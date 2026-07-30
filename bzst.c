// TODO: have a way to set end coord so MT decode doesn't waste cycles.

/* The MIT License

   Copyright (C) 2023, 2026 Genome Research Ltd

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

/*
BZST is a Zstd compatible data file with random access support and
designed for parallel encoding and decoding.

It is a derivation of the earlier BGZF2 draft format.

It is inspired from a union of the zstd seekable format
(https://github.com/facebook/zstd/tree/dev/contrib/seekable_format)
and pzstd (https://github.com/facebook/zstd/tree/dev/contrib/pzstd),
although is not directly compatible with either.

See https://github.com/tfenne/bzst/tree/main/spec for the draft specification.


A Zstd file is a series of frames. Zstd has the notion of data frames
holding compressed data, and skippable frames holding meta-data that
isn't part of the uncompressed output stream.

In BZST
- Each zstd data frame has no back-rerefences to earlier frames and
  can be decoded in isolation.
- Each zstd data frame is immediately preceeded by a block header
  frame (skippable) which holds the size of the adjacent data frame.
  (This is equivalent to pzstd.)
- The file starts with a BZST header frame.
- The file ends with a BZST index frame, which also serves as an EOF
  marker. (This is equivalent to seekable format.)
- Additional format specific skippable frames can occur, either
  prior to the block header frames (eg localised meta-data), or just
  prior to the index frame (meta-data or a second level genomic index
  frame).  This additional layer is built on top of BZST as is its own
  layered format.

TODO: consider also bringing in https://github.com/facebook/zstd/pull/2349
for in-stream dictionaries.  Ideally we'd want multiple dictionaries
and the ability to indicate which dictionary is associated with which
frame.  Maybe the most recent dictionary will always be the one to use?
Combining with parallel decoding this means holding multiple
dictionaries in memory and purging once we know all new blocks relate
to a newer dictionary.  It also means having a random access
capability, so we need to know in the seekable format which entries
are dictionaries.  Eg dsize == 0 and csize != 12, or a separate index
holding the location of dictionary frames so we can use the most
recent one during decode.
*/

/*
Known skippable frame IDs:

0x184D2A50   pzstd, size of next frame
0x184D2A51   aruna footer
0x184D2A52   aruna footer
0x184D2A55   zpkglist - LZ4
0x184D2A56   zpkglist - LZ4
0x184D2A57   zpkglist - LZ4
0x184D2A5B   BZTD (subdivided into BZST frame types)    <----- us
0x184D2A5D   warc-zstd / dict-in-stream
0x184D2A5E   seekable
*/

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <assert.h>
#include <pthread.h>
#include <sys/types.h>
#include <inttypes.h>
#include <zstd.h>
#include <math.h>

#define XXH_STATIC_LINKING_ONLY
#define XXH_IMPLEMENTATION
#include "xxhash.h"
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"
#include "htslib/thread_pool.h"
#include "htslib/bzst.h"
#include "htslib/kstring.h"
#include "cram/pooled_alloc.h"
#include "hts_internal.h" // for hts_bzst_idx_t

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef MAX
#  define MAX(a,b) ((a)>(b)?(a):(b))
#endif

// Mutex debugging
// #define pthread_mutex_lock(x) {fprintf(stderr, "%d: %s:%d lock %p\n", (int)gettid(), __FILE__, __LINE__, x); pthread_mutex_lock(x);}
// #define pthread_mutex_unlock(x) {fprintf(stderr, "%d: %s:%d unlock %p\n", (int)gettid(), __FILE__, __LINE__, x); pthread_mutex_unlock(x);}
// #define pthread_cond_wait(c,m) do {fprintf(stderr, "%d: %s:%d wait %p\n", (int)gettid(), __FILE__, __LINE__, m); pthread_cond_wait(c,m);} while(0)

typedef struct {
    off_t pos;      // cumulative uncompressed position prior to this
    off_t cpos;     // cumulative compression poisition in file
    ssize_t uncomp; // uncompressed size of this block
    ssize_t comp;   // compressed size of this block
} bzst_index_t;

// Genomic index.  We maintain one of these per chromosome, indexed by
// the internal "tid" number.
// TODO: do we use a variable sized integer encoding, or just store in
// 32/64-bit little endian and rely on compression?
typedef struct {
    // TODO: add chr name too so we can do name to tid mappings?
    int tid;             // chromosome+1.  Unused as it's indexed by this?
    hts_pos_t beg, end;  // inclusive range within chromosome

    // or frame number, so we can combine with linear frame index?
    // TODO: rename this.  It's uncomp position and could be within a
    // frame, eg when a frame has multiple references in it.
    off_t frame_start;   // points to pzstd frame prior to data frame
//    size_t frame_offset; // offset within a frame (if multi-ref frames)
//    size_t frame_len;    // length of data from offset onwards

//    // necessary for stats?  Would maybe be better to have nmapped, unmapped,
//    // aligned-bases, etc inline in the pzstd-ish frame, as we can go from
//    // here to there quickly.  We only need one single data point per
//    // chromosome, not per frame, for idxstats meta-data.
//    size_t nmapped;      // number of mapped items in this frame/chr
//    size_t nunmapped;    // number of unmapped items in this frame/chr
} bzst_gindex_t;

// Note this is just a kstring with a linked list.
// We could perhaps embed the link into the ks->s pointer, so empty ones
// are linked together via their data itself.
typedef struct bzst_buffer {
    ssize_t alloc, sz, pos; // allocated, desired size, and curr pos.
    char *buf;
    struct bzst_buffer *next;
    kstring_t meta;
} bzst_buffer;

/*
 * Work-load for single compression or decompression job.
 */
typedef struct bzst_job {
    bzst *fp;
    bzst_buffer *uncomp;
    bzst_buffer *comp;
    int errcode; // +ve is errno, -ve is bzst specific  (TODO)
    int hit_eof; // view from main or view from thread?
    int job_num;
    int known_size;
    struct bzst_job *next;
} bzst_job;

// Message command types for communicating with bzst_mt_reader
enum mtaux_cmd {
    NONE = 0,
    SEEK,
    SEEK_DONE,
    SEEK_FAIL,
    HAS_EOF,
    HAS_EOF_DONE,
    LOAD_SINDEX,
    LOAD_SINDEX_DONE,
    LOAD_GINDEX,
    LOAD_GINDEX_DONE,
    CLOSE,
};

// INTERNAL structure.  Do not use (consider moving to bzst.c and putting
// a blank one here)
struct bzst {
    // Shared with BGZF to provide an easy redirect for code that
    // doesn't want to have lots of conditionals.
    unsigned dummy1:16, is_zstd:1, first_block:1, dummy2:13, is_bzst:1;

    struct hFILE *hfp;    // actual file handle
    int format;           // encoding format (unused, but zlib, zstd, bsc, ...)
    int level;            // compression level
    int is_write;         // open for write
    int block_size;       // ideal block size
    size_t index_frame_sz;// size of seekable index frame
    bzst_index_t *index; // index entries
    int nindex;           // used size of index array
    int aindex;           // allocated size of index array
    int errcode;          // FIXME
    int has_eof;          // Has an EOF block

    off_t frame_pos;      // uncompressed offset of current frame start
    off_t idx_pos;        // as frame_pos, but on explicit flush_try calls only
    off_t tid_pos;        // uncompressed offset of tid within frame.
    off_t last_flush_try; // offset in buffer of last flush attempt
    bzst_buffer *uncomp; // uncompressed data
    bzst_buffer *comp;   // compressed data

    // Genomic indices.
    // We map chr:start-end to uncompressed offsets (use with frame index).
    //
    // OR: CSI style mapping chr/pos to uncompressed offset.
    //     => Same as current, but not a virtual offset. (TODO?)
    int nchr;                // number of chromosomes; size of gindex_sz
    // FIXME: need nchr for 'tid' index, but also number used.
    // Maybe nchr has holes?
    size_t gindex_frame_sz;  // size of zstd frame holding the gindex
    size_t *gindex_sz;       // size of gindex[chr]; consider used+alloc sz
    bzst_gindex_t **gindex; // genomic index per chr / tid
    // For detection of unsorted data.
    int last_used;  // 0 for first rec, 1 for used (sorted), -1 for unsorted
    int last_tid;
    hts_pos_t last_pos;

    // Multi-threading support
    pthread_mutex_t job_pool_m; // when updating this struct
    struct bzst_job *job_free_list;
    int jobs_pending;
    int flush_pending;
    int own_pool;
    hts_tpool *pool;
    hts_tpool_process *out_queue;
    int hit_eof;
    pool_alloc_t *job_pool;
    int job_num; // debugging

    pthread_t io_task;
    pthread_mutex_t command_m; // for messaging with io_task
    pthread_cond_t command_c;
    enum mtaux_cmd command;
    uint64_t seek_to;

    // Shared kstring for optimisation of algorithms
    kstring_t ks;

    // Per frame meta-data callback
    int (*flush_callback)(kstring_t *ks, void *flush_data, int final);
    void *flush_data;
};

static int bzst_add_index(bzst *fp, size_t uncomp, size_t comp);

/*----------------------------------------------------------------------
 * Utility functions shared by both single and multi-threaded implementations
 */

/*
 * Allocates a new bzst buffer capable of holding n bytes.
 *
 * Returns buffer pointer on success,
 *         NULL on failure
 */
bzst_buffer *bzst_buffer_alloc(size_t n) {
    bzst_buffer *b = malloc(sizeof(*b));
    if (!b)
        return NULL;

    if (!(b->buf = malloc(n))) {
        free(b);
        return NULL;
    }

    b->alloc = n;
    b->sz = n;
    b->pos = 0;
    b->next = NULL;

    return b;
}

/*
 * Ensure size is at least n (NB not grow by n).
 * Returns 0 on success,
 *        -1 on failure.
 */
int bzst_buffer_grow(bzst_buffer **bp, size_t n) {
    bzst_buffer *b = *bp;

    if (!*bp) {
        b = *bp = calloc(1, sizeof(*b));
        if (!b)
            return -1;
    }

//    if (n < b->sz)
//      return 0;

    b->sz = n;
    if (n <= b->alloc)
        return 0;

    char *tmp = realloc(b->buf, n);
    if (!tmp)
        return -1;

    b->buf = tmp;
    b->alloc = n;

    return 0;
}

/* Frees memory used by a bzst_buffer */
void bzst_buffer_free(bzst_buffer *b) {
    if (!b)
        return;

    ks_free(&b->meta);
    free(b->buf);
    free(b);
}

static pthread_once_t bzst_comp_once = PTHREAD_ONCE_INIT;
static pthread_key_t bzst_comp_key;

static void bzst_tls_comp_free(void *ptr) {
    ZSTD_CStream *zcs = (ZSTD_CStream *)ptr;
    if (zcs)
        ZSTD_freeCStream(zcs);
}

static void bzst_tls_comp_init(void) {
    pthread_key_create(&bzst_comp_key, bzst_tls_comp_free);
}

/*
 * Compresses a block of data.
 * The output buffer is dynamically grown and can be reused between calls.
 * Initially start with *comp = NULL and *comp_sz = 0.
 *
 * Returns compressed size on success,
 *        -1 on failure.
 */
static ssize_t compress_block(char *uncomp, size_t uncomp_sz,
                              char *comp, size_t comp_alloc,
                              int level) {
    // Check output buffer is large enough
    size_t comp_bound = ZSTD_compressBound(uncomp_sz);
    if (comp_bound > comp_alloc)
        return -1;

    // Currently we use Zstd only.
    // For now we create and configure new streams each time, but we
    // could consider reusing the same zstd stream (NB: needs 1 per thread).
#if 0
    ZSTD_CStream *zcs = ZSTD_createCStream();
    if (!zcs)
        return -1;

    if (ZSTD_initCStream(zcs, level) != 0) {
        ZSTD_freeCStream(zcs);
        return -1;
    }

    ZSTD_CCtx_setParameter(zcs, ZSTD_c_checksumFlag, 1);
    ZSTD_CCtx_setParameter(zcs, ZSTD_c_contentSizeFlag, 1);

    // A bit slower and more system call time.
    ZSTD_inBuffer input = {
        .src = uncomp,
        .size = uncomp_sz,
        .pos = 0
    };
    ZSTD_outBuffer output = {
        .dst = comp,
        .size = comp_alloc,
        .pos = 0
    };

    size_t csize = ZSTD_compressStream2(zcs, &output, &input, ZSTD_e_end);
    return ZSTD_isError(csize) || csize != 0 ? -1 : output.pos;
#elif 1
    // Cache one ZSTD_CStream per running thread.  Little differnce on fast
    // compression, but matters when doing high levels with large blocks as
    // memory usage in the context becomes significant.
    //
    // ./bgzip2 -b8M -11 -@8 enwik9
    // => 18.6s / 2m24s / 1.3s this   -5% real time
    // vs 19.5s / 2m31s / 1.4s below
    //
    // -15
    // => 43.3s / 5m37s / 1.6s this
    // vs 44.3s / 5m41s / 5.0s below
    pthread_once(&bzst_comp_once, bzst_tls_comp_init);
    ZSTD_CStream *zcs = pthread_getspecific(bzst_comp_key);
    if (!zcs) {
        zcs = ZSTD_createCStream();
        pthread_setspecific(bzst_comp_key, zcs);

        if (ZSTD_initCStream(zcs, level) != 0) {
            ZSTD_freeCStream(zcs);
            return -1;
        }
    }
    ZSTD_CCtx_reset(zcs, ZSTD_reset_session_only);
    ZSTD_CCtx_setParameter(zcs, ZSTD_c_checksumFlag, 1);
    ZSTD_CCtx_setParameter(zcs, ZSTD_c_contentSizeFlag, 1);
    // Use a maximum size window log even at low compression levels
//    int wlog = ceil(log(uncomp_sz)/log(2));
//    ZSTD_CCtx_setParameter(zcs, ZSTD_c_windowLog, MIN(24, wlog));

// Helps on bigger buffer sizes (or higher compression levels?)
// Maybe 4-5% smaller for VCF and BCF with similar times, but costs CPU.
// It's essentially just moving the compression level up.
//    ZSTD_CCtx_setParameter(zcs, ZSTD_c_searchLog, 6);
//    ZSTD_CCtx_setParameter(zcs, ZSTD_c_minMatch, 6);

    if (level > 6 && level <= 17) {
        // LDM slows things down and we want to still support the fast
        // modes being fast.  It's often also harmful at the highest
        // compression levels, so we use it for mid-range only.
        //
        // TODO: test novaseq etc.
        ZSTD_CCtx_setParameter(zcs, ZSTD_c_enableLongDistanceMatching, 1);
        ZSTD_CCtx_setParameter(zcs, ZSTD_c_ldmHashLog, 16); // winlog-7 def
        ZSTD_CCtx_setParameter(zcs, ZSTD_c_ldmBucketSizeLog, 3);
        ZSTD_CCtx_setParameter(zcs, ZSTD_c_ldmHashRateLog, 4);
    }

    ZSTD_initCStream(zcs, level);

    size_t csize = ZSTD_compress2(zcs, comp, comp_alloc, uncomp, uncomp_sz);
    return ZSTD_isError(csize) ? -1 : csize;

#else
    ZSTD_CStream *zcs = ZSTD_createCStream();
    if (!zcs)
        return -1;

    if (ZSTD_initCStream(zcs, level) != 0) {
        ZSTD_freeCStream(zcs);
        return -1;
    }

    ZSTD_CCtx_setParameter(zcs, ZSTD_c_checksumFlag, 1);
    ZSTD_CCtx_setParameter(zcs, ZSTD_c_contentSizeFlag, 1);

    size_t csize = ZSTD_compress2(zcs, comp, comp_alloc, uncomp, uncomp_sz);
    ZSTD_freeCStream(zcs);
    return ZSTD_isError(csize) ? -1 : csize;
#endif
}

/*
 * Write a BZST header block as a skippable frame.
 *
 * Format is fixed size:
 * 4:   BGZF2 frame magic number: 0x184D2A5B
 * 4:   remaining frame size
 * 1:   sub-type = 0
 * 1:   format version = 1
 * 4:   "BZST" magic number
 * 4:   format signature (0 == none, but e.g. "BAM\1" etc)
 * 1:   profile (bitmask of abilities used) = 0
 * 1:   flags = 0
 * 8:   XXH64 checksum
 *
 * Returns the number of bytes written on success,
 *         <0 on failure
 */
static int bzst_write_header(bzst *fp) {
    uint8_t buf[28];
    u32_to_le(BZST_SKIPPABLE_ID, buf);
    u32_to_le(20, buf+4);           // size after this
    buf[8] = BZST_HEADER;           // sub-type
    buf[9] = 1;                     // header format
    memcpy(buf+10, "BZST", 4);      // BZST magic
    memcpy(buf+14, "\0\0\0\0", 4);  // parent-format magic
    buf[18] = 0;                    // profile
    buf[19] = 0;                    // flags
    u64_to_le(XXH64(buf, 20, 0), buf+20);

    // Add index entry so offsets work
    if (bzst_add_index(fp, 0, 28) < 0)
        return -1;

    return hwrite(fp->hfp, buf, 28);
}

/*
 * Write a genomic index to be used in conjunction with the seekable index
 * to turn chr:start-end ranges into file offsets.
 *
 * This is more similar to the CRAM index by design, permitting index
 * entries covering all data frames regardless of whether the data is
 * mapped.
 *
 * TODO: maybe this isn't ideal.  We could consider using a CSI index here
 * with virtual offsets be uncompressed positions (to be used with seekable
 * index).
 *
 * Returns 0 on success,
 *        <0 on failure
 */
static int write_genomic_index(bzst *fp) {
    kstring_t ks = {0,0};

    // Header
    ks_resize(&ks, 16); // try 8192
    u32_to_le(BZST_SKIPPABLE_ID, (uint8_t *)ks.s); // BZST skippable frame
    ks.s[8] = GZST_GENOMIC_INDEX;
    ks.s[9] = 0; // index format
    ks.l += 10; // fill out [4..7] later

    // flag
    kputc_(0, &ks); // uncompressed

    // TODO: per file index meta-data.  Basically some bits of "idxstats"

    // Number of chromosomes
    u32_to_le(fp->nchr, (uint8_t *)ks.s + ks.l); ks.l += 4;

    int i;
    for (i = 0; i < fp->nchr; i++) {
        ks_resize(&ks, ks.l + 5 + 20*fp->gindex_sz[i]);

        // flag
        kputc_(0, &ks); // is_aligned, is_sorted... TODO
        // frame count for this chr
        u32_to_le(fp->gindex_sz[i], (uint8_t *)ks.s + ks.l); ks.l += 4;
        // TODO: per-ref meta-data.  Eg other bits of "idxstats"

        bzst_gindex_t *g = fp->gindex[i];
        int j;
        for (j = 0; j < fp->gindex_sz[i]; j++) {
            // Tid isn't needed here.  It belongs out of the loop (so we can
            // have a mismap between tid values to index and array elements).
            u32_to_le(g[j].tid, (uint8_t *)ks.s + ks.l); ks.l += 4;
            // Should we delta these?  Or change to frame no. + offset?
            // FIXME: beg and end are int64.  May want varint.
            u32_to_le(g[j].beg, (uint8_t *)ks.s + ks.l); ks.l += 4;
            u32_to_le(g[j].end, (uint8_t *)ks.s + ks.l); ks.l += 4;
            u64_to_le(g[j].frame_start, (uint8_t *)ks.s + ks.l); ks.l += 8;
        }
    }

    // Footer; used for seeking backwards to start of frame
    ks_resize(&ks, ks.l + 8);
    u32_to_le(ks.l + 8, (uint8_t *)ks.s + ks.l); ks.l += 4;
    u32_to_le(0x8F92EABB, (uint8_t *)ks.s + ks.l); ks.l += 4;

    // Finish up header and write index
    size_t sz = ks.l;
    u32_to_le(sz-8, (uint8_t *)ks.s+4); // size of skippable frame

    //write(3, ks.s, sz);
    int ret = (sz == hwrite(fp->hfp, ks.s, sz) ? 0 : -1);
    free(ks.s);

    return ret;
}

// Submits a command to bzst_mt_reader.  fp->command_m must be unlocked.
static void submit_reader_command(bzst *fp, int command, int done) {
    pthread_mutex_lock(&fp->command_m);

    // Signal the reader thread to execute command
    if (fp->command != CLOSE) {
        fp->command = command;
        fp->errcode = 0;
    }

    pthread_cond_signal(&fp->command_c);
    hts_tpool_wake_dispatch(fp->out_queue);
    do {
        if (fp->command == CLOSE) {
            // possible error in bzst_mt_reader
            pthread_mutex_unlock(&fp->command_m);
            return;
        }
        pthread_cond_signal(&fp->command_c);
        pthread_cond_wait(&fp->command_c, &fp->command_m);
        //fprintf(stderr, "submit_reader_command sees command = %d\n", fp->command);
        if (fp->command == command) {
            // Resend signal intended for bgzf_mt_reader()
            pthread_cond_signal(&fp->command_c);
            break;
        } else if (fp->command == done) {
            // FIXME: can we not use a single DONE flag?
            // Command has completed, check fp->errcode for command status
            break;
        } else if (fp->command == CLOSE) {
            continue;
        } else {
            abort();  // Should not get to any other state
        }
    } while (fp->command != done);
    fp->command = NONE;
    pthread_mutex_unlock(&fp->command_m);
}

static int load_seekable_index_common(bzst *fp);

static int load_seekable_index_mt(bzst *fp) {
    int err = load_seekable_index_common(fp);
    fp->command = LOAD_SINDEX_DONE;
    pthread_cond_signal(&fp->command_c);
    return err;
}

static int load_seekable_index(bzst *fp) {
    int err;
    if (fp->pool) {
        submit_reader_command(fp, LOAD_SINDEX, LOAD_SINDEX_DONE);
        pthread_mutex_lock(&fp->command_m);
        err = fp->errcode;
        pthread_mutex_unlock(&fp->command_m);
    } else {
        err = load_seekable_index_common(fp);
    }
    return err;
}

/*
 * Loads the genomic index is present.
 * Called after reading the seekable index.
 *
 * Returns 0 on success,
 *        -1 on error,
 *        -2 on non-seekable stream,
 *        -3 if no index found.
 */
static int load_genomic_index_common(bzst *fp) {
    uint8_t *buf = NULL;

    if (fp->gindex)
        return 0;

    if (!fp->index) {
        int err;

        if ((err = load_seekable_index_common(fp)) < 0)
            return err;
    }

    // Look for and validate genomic index footer.  This is appear immediately
    // before the seekable index.
    if (hseek(fp->hfp, -(fp->index_frame_sz+8), SEEK_END) < 0)
        return -1 - (errno == ESPIPE);

    uint8_t footer[8];
    if (8 != hread(fp->hfp, footer, 8))
        goto err;

    if (le_to_u32(footer+4) != 0x8F92EABB)
        return -3; // index not found

    // load the index into memory
    uint32_t sz = le_to_u32(footer);
    if (!sz)
        goto err;
    //if (hseek(fp->hfp, -sz, SEEK_CUR) < 0) // why doesn't SEEK_CUR work?
    // FIXME: because it's uint32.  Cast to off_t first.
    if (hseek(fp->hfp, -(fp->index_frame_sz + sz), SEEK_END) < 0)
        goto err;

    if (!(buf = malloc(sz)))
        goto err;
    if (sz != hread(fp->hfp, buf, sz))
        goto err;

    if (le_to_u32(buf) != BZST_SKIPPABLE_ID ||
        buf[8] != GZST_GENOMIC_INDEX ||
        buf[9] != 0) {
        free(buf);
        return -3; // index not found
    }

    // buf[4..7] = skippable frame size, could validate if we wanted to.
    // buf[10] = flag. TODO
    uint8_t *cp = buf+11;
    fp->nchr = le_to_u32(cp);  cp += 4;
//    fprintf(stderr, "Index: nchr %d\n", fp->nchr);

    fp->gindex_frame_sz = sz;
    fp->gindex_sz = calloc(fp->nchr, sizeof(*fp->gindex_sz));
    fp->gindex    = calloc(fp->nchr, sizeof(*fp->gindex));
    if (!fp->gindex_sz || !fp->gindex)
        goto err;

    uint32_t i, j;
    for (i = 0; i < fp->nchr; i++) {
        cp++; //int flag = *cp++; // TODO
        fp->gindex_sz[i] = le_to_u32(cp); cp += 4;
//      fprintf(stderr, "Index: chr %d, nframe %ld\n", i, fp->gindex_sz[i]);

        fp->gindex[i] = calloc(fp->gindex_sz[i], sizeof(*fp->gindex[i]));
        if (!fp->gindex[i])
            goto err;

        bzst_gindex_t *g = fp->gindex[i];
        for (j = 0; j < fp->gindex_sz[i]; j++) {
            g[j].tid = le_to_u32(cp); cp += 4;
            // FIXME: beg and end are int64.  May want varint.
            g[j].beg = le_to_u32(cp); cp += 4;
            g[j].end = le_to_u32(cp); cp += 4;
            g[j].frame_start = le_to_u64(cp); cp += 8;

//          fprintf(stderr, "Index: tid %d, %ld..%ld %ld\n",
//                  g[j].tid, (long)g[j].beg, (long)g[j].end,
//                  g[j].frame_start);
        }
    }

    free(buf);
    return 0;

 err:
    free(buf);
    return -1;
}

static int load_genomic_index_mt(bzst *fp) {
    int err = load_genomic_index_common(fp);
    fp->command = LOAD_GINDEX_DONE;
    pthread_cond_signal(&fp->command_c);
    return err;
}

static int load_genomic_index(bzst *fp) {
    int err;
    if (fp->pool) {
        submit_reader_command(fp, LOAD_GINDEX, LOAD_GINDEX_DONE);
        pthread_mutex_lock(&fp->command_m);
        err = fp->errcode;
        pthread_mutex_unlock(&fp->command_m);
    } else {
        err = load_genomic_index_common(fp);
    }
    return err;
}

/*
 * Finds the uncompressed file offset associated with a specific range.
 * This offset is only suitable for passing into bzst_seek.
 *
 * Note will not likely be the exact location the first record covers this
 * region, but it will be prior to it.  The caller is expected to them
 * discard data out of range.
 *
 * TODO: or should we cache beg/end here and do the discard ourselves?
 * That's how CRAM's API works, and it may be more amenable to a
 * "bzst_query_many" func that can operate on multiple regions in a
 * threaded read-away way.
 *
 * If tid is outside of the bounds of the index, this is an error.
 * If tid is in the index but has no coverage, or we are beyond the end of
 * this reference, then we return the fp offset for the next tid.  This
 * then makes the calling loop immediately hit EOF on range checking.
 *
 * Returns seekable offset on success,
 *        -1 on error,
 *        -2 on non-seekable stream,
 *        -3 if no index found.
 */
int64_t bzst_query(bzst *fp, int tid, hts_pos_t beg, hts_pos_t end) {
    int err;
    if (!fp->gindex) {
        err = load_genomic_index(fp);

        if (err < 0)
            return err;
    }


    // TODO: splitting chromosomes over a frame should provide a new starting
    // point part way into the frame, so we never get the last chromosomes
    // worth of data when doing a tid query. This neatly avoids all questions
    // over whether order of refs in file must match order in header and
    // hence must be incrementing (with -1 wrapping around at end).

    tid++; // unmapped becomes tid 0, even though sorted if diff order.
    if (tid < 0 || tid >= fp->nchr)
        return -1;

    // TODO: Brute force for now, just to get something running.
    // We can replace this with a Nested Containment List (see cram_index.c).
    int i;
    for (i = 0; i < fp->gindex_sz[tid]; i++) {
        if (fp->gindex[tid][i].end < beg)
            continue;

        if (fp->gindex[tid][i].end >= beg)
            return fp->gindex[tid][i].frame_start;
    }

    // Couldn't find it, so use next covered pos
    while (++tid < fp->nchr) {
        if (fp->gindex_sz[tid])
            return fp->gindex[tid][0].frame_start;
    }

    // Could find any, so it's EOF.
    return UINT64_MAX;
}

/*
 * Write an index in the format expected by zstd seekable-format
 * https://github.com/facebook/zstd/blob/dev/contrib/seekable_format/zstd_seekable_compression_format.md
 *
 * We don't add checksums as these are already part of zstd.
 * This index must come last in the file if we wish seekable-format to
 * be able to process this file.
 *
 * Returns 0 on success,
 *        <0 on failure
 */
static int write_seekable_index(bzst *fp) {
    bzst_index_t *idx = fp->index;
    int nidx = fp->nindex;

    uint8_t *buf = malloc(4+4+nidx*8+9);
    int off = 0;
    if (!buf)
        return -1;

    // header
    u32_to_le(0x184D2A5E, buf);
    u32_to_le(nidx*8+9, buf+4);

    // Index entries
    int i;
    for (off = 8, i = 0; i < nidx; i++, off += 8) {
        u32_to_le(idx[i].comp, buf+off);
        u32_to_le(idx[i].uncomp, buf+off+4);
    }

    // Index footer
    u32_to_le(nidx, buf+off); off += 4;
    buf[off++] = 0; // no checksums
    u32_to_le(0x8F92EAB1, buf+off); off += 4;

    int ret = (off == hwrite(fp->hfp, buf, off) ? 0 : -1);
    free(buf);

    return ret;
}

/*
 * Adds an entry to the bzst index
 *
 * Returns 0 on success,
 *        -1 on failure.
 */
static int bzst_add_index(bzst *fp, size_t uncomp, size_t comp) {
    bzst_index_t *idx;

//    static size_t acc = 0; // STATIC but for debugging only
//    fprintf(stderr, "cpos %ld, upos %ld, sz %ld %ld\n",
//          htell(fp->hfp), acc, uncomp, comp);
//    acc += comp;

    // Grow index
    if (fp->nindex >= fp->aindex) {
        size_t n = fp->aindex * 2 + 100;
        if (!(idx = realloc(fp->index, n * sizeof(*idx))))
            return -1;
        fp->aindex = n;
        fp->index = idx;
    }

    // Add
    idx = &fp->index[fp->nindex++];
    idx->comp = comp;
    idx->uncomp = uncomp;
    //idx->pos += uncomp; // not needed while writing

    return 0;
}

/*
 * Writes a per-block meta-data skippable frame holding any format
 * specific additional information per block.  Eg the chromosome range
 * covered by alignments or variants.
 *
 * Format:
 * 4      skippable id: 0x184D2A5B
 * 4      frame size (N+6)
 * 1      sub-type: GZST_BLOCK_META
 * 1      sub-type version (0)
 * N      meta-data
 * 4      checksum (low 32-bit xxhash64)
 *
 * Returns 0 on success,
 *        -1 on failure
 */
static int write_block_metadata(bzst *fp, kstring_t *meta) {
    uint8_t buf[18];

    if (!meta || !meta->l)
	return 0;

    u32_to_le(BZST_SKIPPABLE_ID, buf); // skippable frame id
    u32_to_le(meta->l + 6, buf+4);     // frame size
    buf[8] = GZST_BLOCK_META;          // sub-type
    buf[9] = 0;                        // sub-type version

    int ret = 0;
    ret |= hwrite(fp->hfp, buf, 10) != 10;
    ret |= hwrite(fp->hfp, meta->s, meta->l) != meta->l;

    // or XXH64_state_t *state = XXH64_createState()
    // if we didn't define XXH_STATIC_LINKING_ONLY
    XXH64_state_t state;
    XXH64_reset(&state, 0);
    XXH64_update(&state, buf, 10);
    XXH64_update(&state, meta->s, meta->l);
    u32_to_le(XXH64_digest(&state), buf); // checksum

    ret |= hwrite(fp->hfp, buf, 4) != 4;

    return ret ? -1 : 0;
}

/*
 * Writes a BZST_BLOCK_HEADER skippable frame indicating the size of
 * the next compressed data frame.  We also add it to the index we're
 * building up for zstd seekable.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
static int write_block_header(bzst *fp, ssize_t comp_sz, ssize_t uncomp_sz) {
    uint8_t buf[30];

    u32_to_le(BZST_SKIPPABLE_ID, buf); // skippable ID
    u32_to_le(22, buf + 4);            // frame size
    buf[8] = BZST_BLOCK_HEADER;        // subtype
    u64_to_le(comp_sz, buf+9);         // next block compressed size
    u64_to_le(uncomp_sz, buf+17);      // " uncompressed size
    buf[25] = 0; // flags              // flags (bit 0 = stored; advisory?)
    u32_to_le(XXH64(buf, 0, 26), buf+26);  // low-32 bit checksum 

    int ret = bzst_add_index(fp, 0, 30 + comp_sz);
    ret |= hwrite(fp->hfp, buf, 30) != 30;

    return ret ? -1 : 0;
}

static int write_file_metadata(bzst *fp, kstring_t *meta) {
    uint8_t buf[12];

    u32_to_le(BZST_SKIPPABLE_ID, buf); // pzstd skippable magic no.
    u32_to_le(10 + (meta ? meta->l : 0), buf + 4);
    buf[8] = GZST_FILE_META;
    buf[9] = 0; // meta-data format version

    int ret = bzst_add_index(fp, 0, 18 + (meta ? meta->l : 0));
    ret |= hwrite(fp->hfp, buf, 10) != 10;
    if (meta && meta->l)
        ret |= hwrite(fp->hfp, meta->s, meta->l) != meta->l;

    // pointer back to start for reverse reading
    u32_to_le(10 + meta ? meta->l : 0, buf);
    u32_to_le(0x8F92EA4D, buf+4); // file meta-data magic number
    ret |= hwrite(fp->hfp, buf, 8) != 8;
    // TODO XXhash-64

    return ret ? -1 : 0;
}

// NB: If needed, call this before threading is enabled
static kstring_t *load_file_metadata(bzst *fp, kstring_t *ks) {
    // For ease we load the seekable index and genomic index so we can
    // use the variables there indicating the sizes.
    if (load_genomic_index_common(fp) < 0)
        return NULL;

    off_t off = hseek(fp->hfp, -(fp->index_frame_sz+8 + fp->gindex_frame_sz),
                      SEEK_END);
    if (off < 0)
        return NULL;
    uint8_t footer[8];
    if (8 != hread(fp->hfp, footer, 8))
        return NULL;

    if (le_to_u32(footer+4) != 0x8F92EA4D)
        return NULL; // no meta-data present

    // load the meta-data into memory
    uint32_t sz = le_to_u32(footer);
    if (hseek(fp->hfp, -((off_t)sz+18), SEEK_CUR) < 0)
    if (off < 0)
        return NULL;

    uint8_t header[10];
    if (10 != hread(fp->hfp, header, 10))
        return NULL;
    if (le_to_u32(header) != BZST_SKIPPABLE_ID)
        return NULL;
    if (header[8] != GZST_FILE_META || header[9] != 0)
        return NULL;

    uint32_t mlen = le_to_u32(header+4);
    if (sz+10 != mlen) // not +8? CHECK
        return NULL;

    ks_clear(ks);
    if (ks_resize(ks, mlen+1) < 0)
        return NULL;

    if (mlen != hread(fp->hfp, ks->s, mlen))
        return NULL;
    ks->s[mlen] = 0;
    ks->l = mlen;

    return ks;
}

int bzst_idx_metadata(const hts_idx_t *idx, kstring_t *ks) {
    const hts_bzst_idx_t *bidx = (const hts_bzst_idx_t *)idx;
    return load_file_metadata(bidx->fp, ks) ? 0 : -1;
}

/*----------------------------------------------------------------------
 * Multi-threaded implementation.
 *
 * Encode: the main thread copies data into a block (main thread).  When full
 * this calls bzst_flush to compress and write it out.
 * The bzst_flush function is where things fork, with bzst_write_block
 * performing these operations within the main thread and bzst_write_block_mt
 * creating an encode job to dispatch to the pool of workers.
 *
 * The thread pool workers all call bzst_encode_func, whose only task is
 * the data compression step.
 *
 * bzst_mt_writer is a dedicated I/O thread (not in the pool) whose sole
 * purpoise is to asynchronously take the result of the compression jobs off
 * fp->out_queue and hwrite them to their destination.  This helps remove any
 * potential deadlocks in the main thread and means we can overlap any I/O
 * delays with CPU utilisation.  This is also where we accumulate index
 * entries as it's best performed by a single thread.
 *
 * Finally, we have some tidy up code to drain the queue at the end, as
 * the async nature means the main thread will finish before the workers
 * have completed.
 *
 * Decode: bzst_mt_reader is a dedicated I/O thread (not in the pool) that is
 * reading blocks ahead of time and dispatching decompression jobs to be
 * executed on the main thread pool.  These jobs call bzst_decode_func.
 *
 * We have a "command" channel for communicating with bzst_mt_read, for
 * purposes such as specifying end ranges, closing down, and handling
 * errors.
 */

/*
 * Doesn't explicitly free the memory, but instead puts it on a free list
 * so we can reuse it.
 */
static void bzst_job_free(void *vp) {
    bzst_job *j = (bzst_job *)vp;
    //fprintf(stderr, "bzst_job_free %p %d\n", j, j->job_num);

    bzst *fp = j->fp;
    pthread_mutex_lock(&fp->job_pool_m);
    j->next = fp->job_free_list;
    fp->job_free_list = j;
    pthread_mutex_unlock(&fp->job_pool_m);
}

/*
 * Actually deallocates memory used by a job
 */
static void bzst_job_really_free(bzst_job *j) {
    bzst_buffer_free(j->comp);
    bzst_buffer_free(j->uncomp);
}

/*
 * Fetches a new bzst_job struct.  Reuses previously freed ones if available,
 * or allocates new if not.
 */
static bzst_job *bzst_job_new(bzst *fp) {
    bzst_job *j;

    pthread_mutex_lock(&fp->job_pool_m);
    if (fp->job_free_list) {
        j = fp->job_free_list;
        fp->job_free_list = j->next;
        //fprintf(stderr, "Reuse %d\n", j->job_num);
        pthread_mutex_unlock(&fp->job_pool_m);

    } else {
        pthread_mutex_unlock(&fp->job_pool_m);
        if (!(j = pool_alloc(fp->job_pool)))
            return NULL;
        memset(j, 0, sizeof(*j));
        //fprintf(stderr, "Alloc new block\n");
    }

    bzst_job tmp = *j;
    memset(j, 0, sizeof(*j));
    j->uncomp = tmp.uncomp;
    j->comp   = tmp.comp;
    j->fp     = tmp.fp;

    return j;
}

/*----------
 * MT encoding
 */

/*
 * Single threaded, in a thread of its own (ie async with main).
 *
 * Takes compressed blocks off the results queue and calls hwrite to
 * write them into the output stream.  Also updates any indexing structures.
 *
 * Returns NULL when no more are left, or -1 on error
 */
static void *bzst_mt_writer(void *vp) {
    bzst *fp = (bzst *)vp;
    hts_tpool_result *r;

//    if (fp->idx_build_otf) {
//        fp->idx->moffs = fp->idx->noffs = 1;
//        fp->idx->offs = (bgzidx1_t*) calloc(fp->idx->moffs, sizeof(bgzidx1_t));
//        if (!fp->idx->offs) goto err;
//    }

    // Iterates until result queue is shutdown, where it returns NULL.
    while ((r = hts_tpool_next_result_wait(fp->out_queue))) {
        bzst_job *j = (bzst_job *)hts_tpool_result_data(r);
        assert(j);
//      fprintf(stderr, "job: %d %d->%d %.10s\n",
//              j->job_num, (int)j->uncomp->sz, (int)j->comp->sz,
//              j->uncomp->buf);
//      static int job_last = -1;
//      assert(j->job_num == job_last+1);
//      job_last++;

//        if (fp->idx_build_otf) {
//            fp->idx->noffs++;
//            if ( fp->idx->noffs > fp->idx->moffs )
//            {
//                fp->idx->moffs = fp->idx->noffs;
//                kroundup32(fp->idx->moffs);
//                fp->idx->offs = (bgzidx1_t*) realloc(fp->idx->offs, fp->idx->moffs*sizeof(bgzidx1_t));
//                if ( !fp->idx->offs ) goto err;
//            }
//            fp->idx->offs[ fp->idx->noffs-1 ].uaddr = fp->idx->offs[ fp->idx->noffs-2 ].uaddr + j->uncomp_sz;
//            fp->idx->offs[ fp->idx->noffs-1 ].caddr = fp->idx->offs[ fp->idx->noffs-2 ].caddr + j->comp_sz;
//        }
//
//        // Flush any cached hts_idx_push calls
//        if (bgzf_idx_flush(fp) < 0)
//            goto err;


        if (write_block_metadata(fp, &j->uncomp->meta) < 0)
            goto err;
        if (write_block_header(fp, j->comp->sz, j->uncomp->sz) < 0)
            goto err;

        pthread_mutex_lock(&fp->job_pool_m);
        int ret = bzst_add_index(fp, j->uncomp->pos, j->comp->sz);
        pthread_mutex_unlock(&fp->job_pool_m);
        if (ret < 0)
            goto err;

        if (hwrite(fp->hfp, j->comp->buf, j->comp->sz) != j->comp->sz)
            goto err;

//        // Update our local block_address.  Cannot be fp->block_address due to no
//        // locking in bgzf_tell.
//        pthread_mutex_lock(&mt->idx_m);
//        mt->block_address += j->comp_sz;
//        pthread_mutex_unlock(&mt->idx_m);

        /*
         * Periodically call hflush (which calls fsync when on a file).
         * This avoids the fsync being done at the bgzf_close stage,
         * which can sometimes cause significant delays.  As this is in
         * a separate thread, spreading the sync delays throughout the
         * program execution seems better.
         */
        if (++fp->flush_pending % 32 == 0)
            if (hflush(fp->hfp) != 0)
                goto err;


        hts_tpool_delete_result(r, 0);

        // Also updated by main thread
        pthread_mutex_lock(&fp->job_pool_m);
        fp->jobs_pending--;
        pthread_mutex_unlock(&fp->job_pool_m);
        bzst_job_free(j);
    }

    if (hflush(fp->hfp) != 0)
        goto err;

    hts_tpool_process_destroy(fp->out_queue);
    pthread_mutex_lock(&fp->job_pool_m);
    fp->out_queue = NULL;
    pthread_mutex_unlock(&fp->job_pool_m);

    return NULL;

 err:
    hts_tpool_process_destroy(fp->out_queue);
    return (void *)-1;
}

/*
 * Threaded: runs in one of many thread pool workers.
 *
 * Function called by the bzst encode multi-threaded worker tasks.
 */
static void *bzst_encode_func(void *vp) {
    bzst_job *j = (bzst_job *)vp;

    if ((j->comp->sz = compress_block(j->uncomp->buf, j->uncomp->pos,
                                      j->comp->buf, j->comp->alloc,
                                      j->fp->level)) < 0)
        j->errcode = 1; // TODO design err codes properly

    return vp;
}

/*----------
 * MT decoding
 */

/*
 * Performs the seek (called by reader thread).
 *
 * This simply drains the entire queue, throwing away blocks, seeks,
 * and starts it up again.  Brute force, but maybe sufficient.
 *
 * Returns number of bytes read on success (size of index)
 *        <0 on failure (same codes as load_seekable_index)
 */
static bzst_index_t *index_query(bzst *fp, uint64_t upos);
static void bzst_mt_seek(bzst *fp) {
    if (!fp->index) {
        int err;
        if ((err = load_seekable_index_common(fp)) < 0) {
            pthread_mutex_lock(&fp->job_pool_m);
            fp->errcode = -err;
            fp->command = SEEK_FAIL;
            pthread_mutex_unlock(&fp->job_pool_m);
            return;
        }
    }
    hts_tpool_process_reset(fp->out_queue, 0);
    if (fp->uncomp)
        // triggers a bzst_decode_block_mt call
        fp->uncomp->pos = fp->uncomp->sz;

    pthread_mutex_lock(&fp->job_pool_m);

    bzst_index_t *idx = index_query(fp, fp->seek_to);
    if (!idx) {
        fp->errcode = ENXIO;
        fp->command = SEEK_FAIL;
    } else {
        if (hseek(fp->hfp, idx->cpos, SEEK_SET) < 0) {
            fp->errcode = errno;
            fp->command = SEEK_FAIL;
        } else {
            fp->errcode = 0;
            fp->command = SEEK_DONE;
        }

        // The block is loaded later, as before, but we need to start with
        // pos part way through it.  We do this by modifying seek_to, which
        // is used in bzst_decode_block_mt.
        fp->seek_to -= idx->pos;
        fprintf(stderr, "bgzf seek pos %ld, seek_to %ld\n", idx->cpos, fp->seek_to);
    }
    fp->hit_eof = 0;

    pthread_mutex_unlock(&fp->job_pool_m);

    pthread_cond_signal(&fp->command_c);
}

static int bzst_check_EOF_common(bzst *fp) {
    off_t offset = htell(fp->hfp);
    if (hseek(fp->hfp, -4, SEEK_END) < 0) {
        if (
#ifdef _WIN32
            errno == EINVAL ||
#endif
            errno == ESPIPE) {
            hclearerr(fp->hfp);
            return 2;
        }

        return -1;
    }

    uint8_t buf[4];
    if (hread(fp->hfp, buf, 4) != 4)
        return -1;

    if (hseek(fp->hfp, offset, SEEK_SET) < 0)
        return -1;

    return (memcmp(buf, "\xb1\xea\x92\x8f", 4) == 0) ? 1 : 0;
}

static void bzst_mt_eof(bzst *fp) {
    pthread_mutex_lock(&fp->job_pool_m);
    fp->has_eof = bzst_check_EOF_common(fp);
    pthread_mutex_unlock(&fp->job_pool_m);
    fp->command = HAS_EOF_DONE;
    pthread_cond_signal(&fp->command_c);
}

/*
 * Nul function so we can dispatch a job with the correct serial
 * to mark failure or to indicate an empty read (EOF).
 */
static void *bzst_nul_func(void *arg) {
    return arg;
}

/*
 * Single threaded, called by bzst_mt_reader.
 *
 * Reads a compressed block of data using hread and unwraps the
 * relevant headers to fill out a bzst_job prior to decompression.
 * This doesn't actually dispatch the decompression job.
 *
 * Returns >0 on success (size once uncompressed),
 *          INT_MAX if unknown size (not in frame header),
 *          0 on EOF,
 *         -1 on failure
 */
static int bzst_mt_read_block(bzst *fp, bzst_job *j)
{
    size_t bzst_read_block(bzst *fp, bzst_buffer **comp);

    ssize_t usize = bzst_read_block(fp, &j->comp);
    if (usize == -2)
        return INT_MAX;
    if (usize <= 0)
        return usize;

    // Allocate uncompressed data block
    if (usize != -2) {
        if (usize > BZST_MAX_BLOCK_SIZE)
            return -1; // protect against extreme memory size attacks
        if (bzst_buffer_grow(&j->uncomp, usize) < 0)
            return -1;
        j->known_size = 1;
    } else {
        j->known_size = 0;
    }

    j->fp = fp;
    j->errcode = 0;

    // DEBUGGING ONLY
    static int job_num = 0;
    j->job_num = job_num++;

    return usize;
}

/*
 * Threaded: runs in one of many thread pool workers.
 *
 * Function called by the bzst decode multi-threaded worker tasks.
 */
static void *bzst_decode_func(void *vp) {
    bzst_job *j = (bzst_job *)vp;

    int bzst_decompress_block(bzst_buffer **comp, bzst_buffer **uncomp,
                               size_t usize);

    if (bzst_decompress_block(&j->comp, &j->uncomp,
                               j->known_size ? j->uncomp->sz : -2) < 0)
        j->errcode = 1;

//    if (ZSTD_decompress(j->uncomp->buf, j->uncomp->sz,
//                      j->comp->buf, j->comp->sz) != j->uncomp->sz) {
//      j->errcode = 1; // TODO design err codes properly
//    }

    return vp;
}

/*
 * Single threaded, in a thread of its own (ie async with main).
 *
 * Reads compressed blocks from the file handle and dispatches jobs to the
 * thread pool to decompress them.
 */
static void *bzst_mt_reader(void *vp) {
    bzst *fp = (bzst *)vp;
    bzst_job *j;

    // Note we may have started with load_seekable_index being called by
    // the main application, and we're now in an IO thread. So hread
    // caches represent a data race.  We know however the index isn't
    // loaded after this function starts, so we explicitly lock and unlock
    // to avoid a false positive race.
    //
    // As an alternative, consider making the index load implicit
    // and done on-demand after the initial open, so it's part of this?
    pthread_mutex_lock(&fp->command_m);
    pthread_mutex_unlock(&fp->command_m);

restart:
    if (!(j = bzst_job_new(fp)))
        goto err;
    j->errcode = 0;
    j->hit_eof = 0;
    j->fp = fp;

    while (bzst_mt_read_block(fp, j) > 0) {
        // Dispatch
        if (hts_tpool_dispatch3(fp->pool, fp->out_queue, bzst_decode_func, j,
                                bzst_job_free, bzst_job_free, 0) < 0) {
            bzst_job_free(j);
            goto err;
        }

        // Check for command
        pthread_mutex_lock(&fp->command_m);
        //if (fp->command)
        //    fprintf(stderr, "Reader has command %d\n", fp->command);
        switch (fp->command) {
        case SEEK:
            // Sets fp->command to SEEK_DONE
            bzst_mt_seek(fp);
            pthread_mutex_unlock(&fp->command_m);
            goto restart;

        case HAS_EOF:
            bzst_mt_eof(fp);   // Sets fp->command to HAS_EOF_DONE
            break;

        case LOAD_SINDEX:
            fp->errcode = load_seekable_index_mt(fp);
            break;

        case LOAD_GINDEX:
            fp->errcode = load_genomic_index_mt(fp);
            break;

        case SEEK_DONE:
        case HAS_EOF_DONE:
        case LOAD_SINDEX_DONE:
        case LOAD_GINDEX_DONE:
            pthread_cond_signal(&fp->command_c);
            break;

        case CLOSE:
            pthread_cond_signal(&fp->command_c);
            pthread_mutex_unlock(&fp->command_m);
            hts_tpool_process_destroy(fp->out_queue);
            return NULL;

        default:
            break;
        }
        pthread_mutex_unlock(&fp->command_m);

        // Allocate buffer for next block
        j = bzst_job_new(fp);
        if (!j) {
            hts_tpool_process_destroy(fp->out_queue);
            return NULL;
        }
        j->fp = fp;
    }
    if (j->errcode == 2 /* FIXME */) {
        // Attempt to multi-thread decode a raw gzip stream cannot be done.
        // We tear down the multi-threaded decoder and revert to the old code.
        if (hts_tpool_dispatch3(fp->pool, fp->out_queue, bzst_nul_func, j,
                                bzst_job_free, bzst_job_free, 0) < 0) {
            bzst_job_free(j);
            hts_tpool_process_destroy(fp->out_queue);
            return NULL;
        }
        hts_tpool_process_ref_decr(fp->out_queue);
        return &j->errcode;
    }

    // Dispatch an empty block so EOF is spotted.
    // We also use this mechanism for returning errors, in which case
    // j->errcode is set already.

    j->hit_eof = 1;
    if (hts_tpool_dispatch3(fp->pool, fp->out_queue, bzst_nul_func, j,
                            bzst_job_free, bzst_job_free, 0) < 0) {
        bzst_job_free(j);
        hts_tpool_process_destroy(fp->out_queue);
        return NULL;
    }
    if (j->errcode != 0) {
        hts_tpool_process_destroy(fp->out_queue);
        return &j->errcode;
    }

    // We hit EOF so can stop reading, but we may get a subsequent
    // seek request.  In this case we need to restart the reader.
    //
    // To handle this we wait on a condition variable and then
    // monitor the command. (This could be either seek or close.)
    for (;;) {
        pthread_mutex_lock(&fp->command_m);
        if (fp->command == NONE)
            pthread_cond_wait(&fp->command_c, &fp->command_m);
        switch(fp->command) {
        case SEEK:
            bzst_mt_seek(fp);
            pthread_mutex_unlock(&fp->command_m);
            goto restart;

        case HAS_EOF:
            bzst_mt_eof(fp);   // Sets fp->command to HAS_EOF_DONE
            break;

        case LOAD_SINDEX:
            fp->errcode = load_seekable_index_mt(fp);
            break;

        case LOAD_GINDEX:
            fp->errcode = load_genomic_index_mt(fp);
            break;

        case SEEK_DONE:
        case HAS_EOF_DONE:
        case LOAD_SINDEX_DONE:
        case LOAD_GINDEX_DONE:
            pthread_cond_signal(&fp->command_c);
            break;

        case CLOSE:
            pthread_cond_signal(&fp->command_c);
            pthread_mutex_unlock(&fp->command_m);
            hts_tpool_process_destroy(fp->out_queue);
            return NULL;

        default:
            break;
        }
        pthread_mutex_unlock(&fp->command_m);
    }

 err:
    // We get here if out_queue has been shutdown as we're closing the fp
    pthread_mutex_lock(&fp->command_m);
    fp->command = CLOSE;
    pthread_cond_signal(&fp->command_c);
    pthread_mutex_unlock(&fp->command_m);
    hts_tpool_process_destroy(fp->out_queue);
    return NULL;
}

/*
 * Fetches next decoded block from the thread-pool results queue.
 *
 * Returns >0 on success (number of bytes read),
 *          0 on EOF,
 *         -1 on failure
 *         -2 on switch to stream mode
 */
static int bzst_decode_block_mt(bzst *fp) {
    if (fp->hit_eof)
        return 0;

    hts_tpool_result *r;

    r = hts_tpool_next_result_wait(fp->out_queue);
    if (!r)
        return -1;

    bzst_job *j = (bzst_job *)hts_tpool_result_data(r);
    hts_tpool_delete_result(r, 0);
    if (j->hit_eof) {
        fp->hit_eof = 1;
        bzst_job_free(j);
        return 0;
    }

    // Swap j->uncomp with fp->uncomp to avoid a memcpy
    bzst_buffer *tmp = fp->uncomp;
    fp->uncomp = j->uncomp;
    j->uncomp = tmp;

    bzst_job_free(j);

    // We may need to start part way into this block if we've just done a
    // seek.  Seek_to is reset from absolute uncompressed offset before the
    // seek, to relative offset to seek block after the seek.
    fp->uncomp->pos = fp->seek_to;
    fp->seek_to = 0;

    return fp->uncomp->sz;
}

/*----------
 * Multi-threading code shared by reader and writers
 */

/*
 * Main thread: adds a thread pool to a bzst fp, adding
 * multi-threading capabilities.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int bzst_thread_pool(bzst *fp, hts_tpool *pool, int qsize) {
    fp->own_pool = 0;
    fp->pool = pool;
    if (!qsize)
        qsize = hts_tpool_size(pool)*2;
    if (!(fp->out_queue = hts_tpool_process_init(fp->pool, qsize, 0)))
        return -1;
    hts_tpool_process_ref_incr(fp->out_queue);

    fp->job_pool = pool_create(sizeof(bzst_job));
    if (!fp->job_pool)
        return -1;

    pthread_mutex_init(&fp->job_pool_m, NULL);
    pthread_mutex_init(&fp->command_m, NULL);
//    pthread_mutex_init(&fp->idx_m, NULL);
    pthread_cond_init(&fp->command_c, NULL);
    fp->flush_pending = 0;
    fp->jobs_pending = 0;

    pthread_create(&fp->io_task, NULL,
                   fp->is_write ? bzst_mt_writer : bzst_mt_reader, fp);
    return 0;
}

/*----------------------------------------------------------------------
 * Basic implementation
 */

/*
 * Compress and write a block to disk.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
static int bzst_write_block(bzst *fp, bzst_buffer *buf) {
    if (bzst_buffer_grow(&fp->comp, ZSTD_compressBound(buf->sz)) < 0)
        return -1;

    if ((fp->comp->sz = compress_block(buf->buf, buf->pos,
                                       fp->comp->buf, fp->comp->alloc,
                                       fp->level)) < 0)
        return -1;

    if (write_block_metadata(fp, &buf->meta) < 0)
	return -1;
    if (write_block_header(fp, fp->comp->sz, fp->uncomp->sz) < 0)
        return -1;

    int ret = bzst_add_index(fp, buf->pos, fp->comp->sz);

    if (hwrite(fp->hfp, fp->comp->buf, fp->comp->sz) != fp->comp->sz)
        return -1;

    fp->uncomp->pos = 0; // reset buffered offset

    // We may have temporarily grown it, but keep honouring the user
    // requested block size so one large record doesn't permanently increase
    // block size.
    fp->uncomp->sz = fp->block_size;
    return ret;
}

/*
 * (Runs in the main thread)
 * Asynchronously compress and write a block to disk.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
static int bzst_write_block_mt(bzst *fp, bzst_buffer *buf) {
    bzst_job *j = bzst_job_new(fp);
    if (!j)
        return -1;

    // Fill out job.
    size_t comp_sz = ZSTD_compressBound(buf->sz);
    j->fp = fp;
    j->job_num = fp->job_num++;
    j->next = NULL;

    // We could steal "buf" and optimise this more.
    if (bzst_buffer_grow(&j->uncomp, buf->sz) < 0 ||
        bzst_buffer_grow(&j->comp, comp_sz) < 0)
        return -1;
    memcpy(j->uncomp->buf, buf->buf, buf->pos);
    j->uncomp->pos = buf->pos;
    // It may have grown, but shrink back again
    buf->sz = fp->block_size;

    // Swap buffers.
    kstring_t tmp = j->uncomp->meta;
    j->uncomp->meta = buf->meta;
    buf->meta = tmp;
    ks_clear(&buf->meta);

    // Used by flush function to drain queue.
    // Can we get this via another route?
    pthread_mutex_lock(&fp->job_pool_m);
    fp->jobs_pending++;
    pthread_mutex_unlock(&fp->job_pool_m);

    if (hts_tpool_dispatch3(fp->pool, fp->out_queue, bzst_encode_func,
                            j, bzst_job_free, bzst_job_free, 0) < 0)
        goto err;

    fp->uncomp->pos = 0;
    return 0;

 err:
    bzst_job_free(j);
    pthread_mutex_lock(&fp->job_pool_m);
    fp->jobs_pending--;
    pthread_mutex_unlock(&fp->job_pool_m);
    return -1;
}

/*
 * Flush the bzst stream and ensure we start a new block.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int bzst_flush(bzst *fp) {
    int ret = 0;
    if (!fp->uncomp->pos)
        return 0;

    if (fp->first_block) {
        fp->first_block = 0;
        ret = bzst_write_header(fp) < 0;
    }

    // uncompressed position of next frame
    fp->frame_pos += fp->uncomp->pos;

    // Arbitrary callback to accumulate per-frame meta-data stats
    if (fp->flush_callback)
        ret |= fp->flush_callback(&fp->uncomp->meta, fp->flush_data, 0) < 0;

    if (fp->pool) {
        ret |= bzst_write_block_mt(fp, fp->uncomp);
    } else {
        ret |= bzst_write_block(fp, fp->uncomp);
    }

    fp->last_flush_try = 0;

    return ret;
}

/*
 * Tests whether a write of 'size' would spill over to the next block. If so
 * flush this current one, so we always end blocks on a whole record.
 *
 * Returns 0 on success,
 *        <0 on failure
 */
int bzst_flush_try(bzst *fp, ssize_t size) {
    if (fp->uncomp && fp->uncomp->pos + size > fp->uncomp->sz) {
        int ret = bzst_flush(fp);
        fp->idx_pos = fp->frame_pos;

        // We know we're trying to push size bytes, and we don't want to split
        // write requests, so we'll grow the current buffer if needed.
        if (fp->uncomp && size > fp->uncomp->sz)
            ret |= bzst_buffer_grow(&fp->uncomp, size);

        return ret;
    }

    fp->last_flush_try = fp->uncomp->pos;
    fp->idx_pos = fp->frame_pos;

    return 0;
}

/*
 * Flushes data and drains any asynchronous I/O still waiting to be
 * written.
 *
 * Returns 0 on success,
 *        <0 on failure
 */
int bzst_drain(bzst *fp) {
    if (bzst_flush(fp) < 0)
        return -1;

    if (!fp->pool)
        return 0;

    // Wait for any compression jobs still in the queue or running
    // to complete.  Even after the queue is drained, we may still
    // have jobs_pending > 0 as we're waiting on hwrite I/O too.
    int jp = 0;
    do {
        if (hts_tpool_process_flush(fp->out_queue) < 0)
            return -1;
        pthread_mutex_lock(&fp->job_pool_m);
        jp = fp->jobs_pending;
        pthread_mutex_unlock(&fp->job_pool_m);
    } while (jp);

    // At this point the writer thread should have completed too
    pthread_mutex_lock(&fp->job_pool_m);
    hts_tpool_process_destroy(fp->out_queue);
    pthread_mutex_unlock(&fp->job_pool_m);

    void *retval = NULL;
    pthread_join(fp->io_task, &retval);

    return retval == NULL ? 0 : -1;
}

/*
 * Set the bzst block size.  This can be performed at any point,
 * but it is usually done immediately after opening for write.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int bzst_set_block_size(bzst *fp, size_t sz) {
    if (!fp->is_write)
        return 0;

    fp->block_size = sz;

    if (fp->uncomp)
        if (bzst_flush(fp) < 0)
            return -1;

//    if (fp->first_block)
//      sz = sz < 1024 ? sz : 1024;

    return bzst_buffer_grow(&fp->uncomp, sz);
}

void bzst_set_level(bzst *fp, int level) {
    fprintf(stderr, "bgzf level %d\n", level);
    fp->level = level;
}

static bzst *bzst_open_common(bzst *fp, hFILE *hfp, const char *mode) {
    fp->is_bzst = 1;
    fp->first_block = 1;
    fp->hfp = hfp;
    fp->frame_pos = 0;

    if (*mode == 'w') {
        fp->is_write = 1;
        if (bzst_set_block_size(fp, BZST_DEFAULT_BLOCK_SIZE) < 0) {
            if (hclose(fp->hfp)) // can't ignore :/
                free(fp);
            else
                free(fp);
            return NULL;
        }

        int level = BZST_DEFAULT_LEVEL;
        for (const char *cp = mode; *cp; cp++) {
            if (*cp >= '0' && *cp <= '9') {
                level = atoi(cp);
                break;
            }
        }

        // We could create a ZSTD_CStream here for reuse.
        // This would then also enable us to cache all the params, such
        // as long mode, minMatch, etc.
        fp->level = level;
    } else {
        fp->is_write = 0;
    }

//    // Threading
//    fp->pool = hts_tpool_init(4);
//    fp->out_queue = hts_tpool_process_init(fp->pool, 8, 0);
//    if (bzst_thread_pool(fp, fp->pool, 0) < 0)
//      goto err;

    return fp;
}

/*
 * Opens a bzst file from an existing hfile for read ("r") or write
 * ("w", or "w1" to "w19").
 *
 * Returns bzst file handle on success,
 *         or NULL on failure.
 */
bzst *bzst_hopen(hFILE *hfp, const char *mode) {
    bzst *fp = calloc(1, sizeof(*fp));
    if (!fp)
        return NULL;

    return bzst_open_common(fp, hfp, mode);
}
/*
 * Opens a bzst file 'fn' for read ("r") or write ("w", or "w1" to "w19").
 *
 * Returns bzst file handle on success,
 *         or NULL on failure.
 */
bzst *bzst_open(const char *fn, const char *mode) {
    bzst *fp = calloc(1, sizeof(*fp));
    if (!fp)
        return NULL;

    fp->hfp = hopen(fn, mode);
    if (!fp->hfp) {
        free(fp);
        return NULL;
    }

    return bzst_open_common(fp, fp->hfp, mode);
}


/*
 * Closes a bzst file.
 *
 * Returns 0 on success,
 *        <0 on failure
 */
int bzst_close(bzst *fp) {
    int ret = 0;

    if (fp->is_write) {
        ret |= (bzst_drain(fp) < 0);

        if (fp->flush_callback) {
            kstring_t final_meta = KS_INITIALIZE;
            fp->flush_callback(&final_meta, fp->flush_data, 1);
            ret |= write_file_metadata(fp, &final_meta);
            free(final_meta.s);
        }

        ret |= write_genomic_index(fp);
        ret |= write_seekable_index(fp);
    }

    if (fp->pool) {
        if (!fp->is_write) {
            // Tell the reader to shutdown and wait for it
            pthread_mutex_lock(&fp->command_m);
            fp->command = CLOSE;
            pthread_cond_signal(&fp->command_c);
            hts_tpool_wake_dispatch(fp->out_queue); // unstick the reader
            pthread_mutex_unlock(&fp->command_m);

            if (hts_tpool_process_is_shutdown(fp->out_queue) > 1)
                ret = -1;
        }

        // out_queue has a ref count, so can we this both here and in threads
        hts_tpool_process_destroy(fp->out_queue);

        void *retval = NULL;
        pthread_join(fp->io_task, &retval);
        ret = retval != NULL ? -1 : ret;

        pthread_mutex_destroy(&fp->job_pool_m);
        pthread_mutex_destroy(&fp->command_m);
        pthread_cond_destroy(&fp->command_c);

        bzst_job *j, *n;
        for (j = fp->job_free_list; j; j = n) {
            //fprintf(stderr, "free job %d\n", j->job_num);
            n = j->next;
            bzst_job_really_free(j);
        }

        pool_destroy(fp->job_pool);
        if (fp->own_pool)
            hts_tpool_destroy(fp->pool);
    }

    if (fp->comp)   bzst_buffer_free(fp->comp);
    if (fp->uncomp) bzst_buffer_free(fp->uncomp);

    ret |= hclose(fp->hfp);

    int i;
    for (i = 0; i < fp->nchr; i++)
        free(fp->gindex[i]);
    free(fp->gindex);
    free(fp->gindex_sz);

    if (fp->flush_callback)
        fp->flush_callback(NULL, fp->flush_data, 0);

    free(fp->index);
    free(fp->ks.s);
    free(fp);

    return ret ? -1 : 0;
}

// Registers a callback function and a block of data with the bzst fp.
// We call this callback on each flush and use the results to generate
// arbitrary per-frame meta-data.
void bzst_add_flush_callback(bzst *fp, void *flush_data,
                              int (*flush_callback)(kstring_t *,
                                                    void *,
                                                    int final)) {
    fp->flush_data = flush_data;
    fp->flush_callback = flush_callback;
}

// Returns the registered bzst flush_data
void *bzst_flush_data(bzst *fp) {
    return fp->flush_data;
}

/*
 * Writes a block of data to a bzst file handle.
 *
 * If "can_split" is set then buf may be split into two indexable chunks.
 * Otherwise buf will never span two blocks.
 *
 * Returns number of bytes written on success
 *        -1 on failure
 */
int bzst_write(bzst *fp, const char *buf, size_t buf_sz, int can_split) {
    int r = 0;

    do {
        if (fp->uncomp && fp->uncomp->sz == fp->uncomp->pos) {
            // Full
            if (bzst_flush(fp))
                return -1;
        }

        if (!fp->uncomp)
            if (bzst_buffer_grow(&fp->uncomp, fp->block_size) < 0)
                return -1;

        ssize_t consumes = fp->uncomp->sz - fp->uncomp->pos;
        if (consumes > buf_sz)
            consumes = buf_sz;

        if (consumes == buf_sz || can_split) {
            // Whole or can_split and a portion
            memcpy(fp->uncomp->buf + fp->uncomp->pos, buf, consumes);
            fp->uncomp->pos += consumes;
            buf_sz -= consumes;
            buf    += consumes;
            r      += consumes;
        } else {
            // Can't split and doesn't fit, so flush and write this new item
            // as a new block if it's too large.
            int err = bzst_flush(fp);
            if (err || bzst_buffer_grow(&fp->uncomp,
                                         MAX(buf_sz, fp->block_size)) < 0)
                return -1;
            memcpy(fp->uncomp->buf, buf, buf_sz);

            // else it will fit on the next loop given we've now flushed
        }
    } while (buf_sz);

    return r;
}

/*
 * Reads more compressed data into the bzst_buffer.
 * NB: The return value is the size of the data when uncompressed.
 * The compressed size will be in (*comp)->sz.
 *
 * Returns >0 on success (size when uncompressed),
 *          0 on EOF,
 *         -1 on failure
 *         -2 on unknown sized blocks (eg pzstd format)
 *         -3 on switch to stream mode
 */
/*static*/ size_t bzst_read_block(bzst *fp, bzst_buffer **comp) {
    uint8_t buf[14];
    size_t n;
    uint32_t meta_sz;

 next_block:
    if (8 != (n = hread(fp->hfp, buf, 8)))
        return n == 0 ? 0 : -1;

    uint32_t fsize = le_to_u32(buf+4);

    if (le_to_u32(buf) != BZST_SKIPPABLE_ID) {
        if ((le_to_u32(buf) & 0x184D2A50) == 0x184D2A50) {
            // Another skippable frame, so skip it
            char tmp[8192];
            size_t n, c = fsize;
            while (c > 0 && (n = hread(fp->hfp, tmp, MIN(8192, c))) > 0)
                c -= n;

            if (c)
                return -1;

            goto next_block;
        }
        return -3;
    } else {
        meta_sz = le_to_u32(buf+4);
        if (meta_sz > 100000) {
            fprintf(stderr, "Unexpectedly large meta-data size.  Aborting\n");
            return -3;
        }
        uint8_t *meta = malloc(meta_sz);
        if (!meta)
            return -1;

        // remainder of pzstd skippable frame, now we know it is one
        if (meta_sz != hread(fp->hfp, meta, meta_sz)) {
            free(meta);
            return -1;
        }

        bzst_frame_t type = meta[0];
        if (meta_sz < 6) {
            free(meta);
            return -1;
        }

        memcpy(buf+8, meta, 6); // csize

        // FIXME: put it into buffer kstring
        free(meta);

        if (type != GZST_BLOCK_META)
            goto next_block;
    }

    // Load compressed data
    if (meta_sz < 6)
        return -1;
    size_t csize = le_to_u32(buf+10);
    if (bzst_buffer_grow(comp, csize) < 0)
        return -1;
    if (csize != hread(fp->hfp, (*comp)->buf, csize))
        return -1;

    // Get decompressed size and return it.
    size_t usize = ZSTD_getFrameContentSize((*comp)->buf, csize);
    if (usize == ZSTD_CONTENTSIZE_UNKNOWN)
        return -2;
    if (usize == ZSTD_CONTENTSIZE_ERROR)
        return -1;
    if (usize == 0) {
        return 0; // empty frame => skip to next
        // FIXME: this isn't an EOF, so try again with a goto next_block?
    }

    return usize;
}

static pthread_once_t bzst_once = PTHREAD_ONCE_INIT;
static pthread_key_t bzst_key;

static void bzst_tls_free(void *ptr) {
    ZSTD_DStream *zds = (ZSTD_DStream *)ptr;
    if (zds)
        ZSTD_freeDStream(zds);
}

static void bzst_tls_init(void) {
    pthread_key_create(&bzst_key, bzst_tls_free);
}

/*static*/
int bzst_decompress_block(bzst_buffer **comp, bzst_buffer **uncomp,
                           size_t usize) {
    // Allocate uncompressed data block
    if (usize == -2) {
        // Unknown size, so we have to dynamically grow and do multiple
        // decode calls.  This isn't so efficient, but it's needed if
        // we want to process pzstd output streams as it doesn't add the
        // content size to the frame.
        if (*comp)
            (*comp)->pos = 0;
        if (*uncomp)
            (*uncomp)->pos = 0;
        if (usize != -2) {
            if (bzst_buffer_grow(uncomp, usize) < 0)
                return -1;
        }
        if (!*uncomp || (*uncomp)->alloc == 0)
            // first cautious guess
            if (bzst_buffer_grow(uncomp, (*comp)->sz*4 + 8192) < 0)
                return -1;
        ZSTD_inBuffer input = {
            .src  = (*comp)->buf,
            .size = (*comp)->sz,
            .pos  = 0
        };
        ZSTD_outBuffer output = {
            .dst  = (*uncomp)->buf,
            .size = (*uncomp)->alloc,
            .pos  = 0
        };
        size_t dsize = 1;

        // Cache ZSTD_DStream, one per thread.
        // This is an 8% speed gain in a test, but that may not be
        // representative.
        pthread_once(&bzst_once, bzst_tls_init);
        ZSTD_DStream *zds = pthread_getspecific(bzst_key);
        if (!zds) {
            zds = ZSTD_createDStream();
            pthread_setspecific(bzst_key, zds);
        }
        //ZSTD_DCtx_reset(zds, ZSTD_reset_session_only);

        //ZSTD_DStream *zds = ZSTD_createDStream();
        if (!zds)
            return -1;

        // Keep iterating until all input stream is consumed
        while (input.pos < input.size) {
            //fprintf(stderr, "In %ld out %ld\n", input.size, output.size);
            dsize = ZSTD_decompressStream(zds, &output, &input);
            if (ZSTD_isError(dsize)) {
                fprintf(stderr, "Error: %s\n", ZSTD_getErrorName(dsize));
                return -1;
            }
            // Grow output buffer based on proportion of input, accounting
            // for potential improvements in ratio as we extend.
            if (input.pos < input.size && output.size == output.pos) {
                size_t guess = (input.size+1.0) / (input.pos+1.0) * 1.05 *
                    output.size + 1000;
                guess = guess > (*uncomp)->alloc + 10000
                      ? guess : (*uncomp)->alloc + 10000;
                //fprintf(stderr, "Grow %ld to %ld\n", (*uncomp)->alloc, guess);
                if (bzst_buffer_grow(uncomp, guess) < 0)
                    return -1;

                output.dst  = (*uncomp)->buf;
                output.size = (*uncomp)->alloc;
            }
        }

        // We've used all input, but may still be blocked on output size.
        while (dsize != 0) {
            fprintf(stderr, "Second zstd decode\n");
            if ((*uncomp)->pos == (*uncomp)->sz)
                if (bzst_buffer_grow(uncomp, (*uncomp)->alloc*1.5+100000) < 0)
                    return -1;

            output.dst  = (*uncomp)->buf;
            output.size = (*uncomp)->alloc;

            dsize = ZSTD_decompressStream(zds, &output, &input);
            if (ZSTD_isError(dsize)) {
                fprintf(stderr, "Error: %s\n", ZSTD_getErrorName(dsize));
                return -1;
            }
        }
        (*uncomp)->sz = output.pos;

        //ZSTD_freeDStream(zds);  // uncomment if not using TLS

    } else {
        // Note: we could use the same logic above, but seeding the
        // fp->uncomp buffer size with usize instead of a starting guess.
        //
        // This appears to be similar speed to this code below, but I'm
        // wary as it's much less clear as to how it's working.
        //
        // FIXME: merge at some point when we're more confident in robustness
        // of pzstd support.
        if (usize > BZST_MAX_BLOCK_SIZE)
            return -1; // protect against extreme memory size attacks
        if (bzst_buffer_grow(uncomp, usize) < 0)
            return -1;
        assert((*uncomp)->sz == usize);

        // Uncompress it
        if (ZSTD_decompress((*uncomp)->buf, (*uncomp)->sz,
                            (*comp)->buf, (*comp)->sz) != (*uncomp)->sz) {
            // Embedded size mismatches the stream
            return -1;
        }
    }

    (*uncomp)->pos = 0;

    return (*uncomp)->sz;
}

/*
 * Reads more compressed data and decompresses it.
 * This can work on both non-random access zstd compressed files (TODO)
 * as well as block compressed ones.
 *
 * This fills out fp->uncomp and resets to fp->uncomp->pos to be
 * the start/end of decoded size.
 *
 * TODO: just use the standard zstd API if not pzstd framed
 * We may wish to just switch to streaming mode, but we'll want to push
 * data back to the buffer first?
 *
 * Returns >0 on success (number of bytes read),
 *          0 on EOF,
 *         -1 on failure
 *         -2 on switch to stream mode
 */
static int bzst_decode_block(bzst *fp) {
    ssize_t usize = bzst_read_block(fp, &fp->comp);
    if (usize != -2 && usize <= 0)
        return usize;

    return bzst_decompress_block(&fp->comp, &fp->uncomp, usize);

    // Allocate uncompressed data block
    if (usize == -2) {
        // Unknown size, so we have to dynamically grow and do multiple
        // decode calls.  This isn't so efficient, but it's needed if
        // we want to process pzstd output streams as it doesn't add the
        // content size to the frame.
        if (fp->comp)
            fp->comp->pos = 0;
        if (fp->uncomp)
            fp->uncomp->pos = 0;
        if (!fp->uncomp || fp->uncomp->sz == 0)
            if (bzst_buffer_grow(&fp->uncomp, fp->comp->sz*2) < 0)
                return -1;
        ZSTD_inBuffer input = {
            .src  = fp->comp->buf,
            .size = fp->comp->sz,
            .pos  = 0
        };
        ZSTD_outBuffer output = {
            .dst  = fp->uncomp->buf,
            .size = fp->uncomp->sz,
            .pos  = 0
        };
        size_t dsize;
        ZSTD_DStream *zcs = ZSTD_createDStream();
        if (!zcs)
            return -1;

        // Keep iterating until all input stream is consumed
        while (input.pos < input.size) {
            dsize = ZSTD_decompressStream(zcs, &output, &input);
            if (ZSTD_isError(dsize)) {
                fprintf(stderr, "Error: %s\n", ZSTD_getErrorName(dsize));
                return -1;
            }
            // Grow output buffer based on proportion of input, accounting
            // for potential improvements in ratio as we extend.
            if (input.pos < input.size && output.size == output.pos) {
                size_t guess = (input.size+1.0) / (input.pos+1.0) * 1.05 *
                    output.size + 1000;
                guess = guess > fp->uncomp->sz + 10000
                      ? guess : fp->uncomp->sz + 10000;
                if (bzst_buffer_grow(&fp->uncomp, guess) < 0)
                    return -1;

                output.dst  = fp->uncomp->buf;
                output.size = fp->uncomp->sz;
            }
        }

        // We've used all input, but may still be blocked on output size.
        while (dsize != 0) {
            if (bzst_buffer_grow(&fp->uncomp,
                                  fp->uncomp->sz * 1.5 + 100000) < 0)
                return -1;

            output.dst  = fp->uncomp->buf;
            output.size = fp->uncomp->sz;

            dsize = ZSTD_decompressStream(zcs, &output, &input);
            if (ZSTD_isError(dsize)) {
                fprintf(stderr, "Error: %s\n", ZSTD_getErrorName(dsize));
                return -1;
            }
        }
        fp->uncomp->sz = output.pos;

    } else {
        // FIXME: we could use the same logic above, but seeding the
        // fp->uncomp buffer size with usize instead of a starting guess.
        if (usize > BZST_MAX_BLOCK_SIZE)
            return -1; // protect against extreme memory size attacks
        if (bzst_buffer_grow(&fp->uncomp, usize) < 0)
            return -1;
        assert(fp->uncomp->sz == usize);

        // Uncompress it
        if (ZSTD_decompress(fp->uncomp->buf, fp->uncomp->sz,
                            fp->comp->buf, fp->comp->sz) != fp->uncomp->sz) {
            // Embedded size mismatches the stream
            return -1;
        }
    }

    fp->uncomp->pos = 0;

    return fp->uncomp->sz;
}

/*
 * Refill fp->uncomp when we're out of uncompressed data.
 * Returns >0 on success (number of bytes refilled)
 *        -1 on eof
 *        -2 on failure
 */
static int bzst_refill_uncomp(bzst *fp) {
    ssize_t n = 0;

    if (!fp->uncomp || fp->uncomp->pos == fp->uncomp->sz) {
        // out of buffered content, fetch some more
        n = fp->pool
            ? bzst_decode_block_mt(fp)
            : bzst_decode_block(fp);
        if (n < 0) {
            return -2; // err
        } else if (n == 0) {
            fp->hit_eof = 1;
            return -1; // eof
        }
    }

    return n;
}

/*
 * Reads a data from a bzst file handle.
 *
 * Returns number of bytes read on success
 *        -1 on failure
 */
int bzst_read(bzst *fp, char *buf, size_t buf_sz) {
    if (fp->hit_eof)
        return 0;

    size_t decoded = 0;

    // loop:
    //    fill buffer if empty
    //        read a block to comp, decompress it, store in buffer
    //    copy from buffer to buf
    while (buf_sz) {
#if 1
        switch (bzst_refill_uncomp(fp)) {
        case -1:  return decoded; // EOF
        case -2:  return -1;      // error
        }
#else
        if (!fp->uncomp || fp->uncomp->pos == fp->uncomp->sz) {
            // out of buffered content, fetch some more
            size_t n = fp->pool
                ? bzst_decode_block_mt(fp)
                : bzst_decode_block(fp);
            if (n < 0)
                return -1;
            else if (n == 0) {
                fp->hit_eof = 1;
                return decoded; // EOF
            }
        }
#endif

        size_t n = MIN(buf_sz, fp->uncomp->sz - fp->uncomp->pos);
        memcpy(buf, fp->uncomp->buf + fp->uncomp->pos, n);
        buf += n;
        buf_sz -= n;
        fp->uncomp->pos += n;
        decoded += n;
    }

    return decoded;
}

/*
 * Reads a data from a bzst file handle.  This modifies *buf to
 * point to a block of internal data and returns the size of this data.
 * In will be between 1 and buf_sz bytes long.  This data should not be
 * modified.
 *
 * Returns number of bytes read on success
 *        -1 on failure
 */
int bzst_read_zero_copy(bzst *fp, const char **buf, size_t buf_sz) {
    if (fp->hit_eof)
        return 0;

    size_t decoded = 0;
    *buf = NULL;

    if (!buf_sz)
        return 0;

    if (!fp->uncomp || fp->uncomp->pos == fp->uncomp->sz) {
        // out of buffered content, fetch some more
        ssize_t n = fp->pool
            ? bzst_decode_block_mt(fp)
            : bzst_decode_block(fp);
        if (n < 0)
            return -1;
        else if (n == 0) {
            fp->hit_eof = 1;
            return decoded; // EOF
        }
    }

    size_t n = MIN(buf_sz, fp->uncomp->sz - fp->uncomp->pos);
    *buf = fp->uncomp->buf + fp->uncomp->pos;
    buf += n;
    buf_sz -= n;
    fp->uncomp->pos += n;
    decoded += n;

    return decoded;
}

/*
 * Loads a seekable index from a bzst file.
 *
 * Returns 0 on success,
 *        -1 on error,
 *        -2 on non-seekable stream,
 *        -3 if no index found.
 */
static int load_seekable_index_common(bzst *fp) {
    // Look for and validate index footer
    off_t off = hseek(fp->hfp, -9, SEEK_END);
    if (off < 0)
        return -1 - (errno == ESPIPE);

    uint8_t footer[9];
    if (9 != hread(fp->hfp, footer, 9))
        return -1;

    if (le_to_u32(footer+5) != 0x8F92EAB1 || (footer[4] & 0x7C) != 0)
        return -3;
    int has_chksum = footer[4] & 0x80 ? 1 : 0;

    // Read entire index frame
    uint32_t nframes = le_to_u32(footer);
    size_t sz = 9 + nframes*4*(2+has_chksum) + 8;
    off = hseek(fp->hfp, -sz, SEEK_END);
    if (off < 0)
        return -1;

    bzst_index_t *idx = NULL;
    uint8_t *buf = malloc(sz);
    if (!buf)
        return -1;

    if (sz != hread(fp->hfp, buf, sz))
        goto err;
    if (le_to_u32(buf) != 0x184D2A5E)
        return -3;
    if (le_to_u32(buf+4) != sz-8)
        return -3;

    // Decode index
    fp->index_frame_sz = sz;
    fp->index = idx = malloc(nframes * sizeof(*idx));
    fp->nindex = fp->aindex = nframes;
    if (!idx)
        goto err;

    uint32_t i;
    uint64_t pos = 0, cpos = 0;
    uint8_t *idxp = buf + 8;
    for (i = 0; i < nframes; i++) {
        idx[i].pos    = pos;
        idx[i].cpos   = cpos;
        idx[i].comp   = le_to_u32(idxp);
        idx[i].uncomp = le_to_u32(idxp+4);
//      fprintf(stderr, "frame %ld..%ld  %ld..%ld\n",
//              cpos, cpos+idx[i].comp,
//              pos, pos+idx[i].uncomp);
        idxp += 4*(2+has_chksum);
        pos += idx[i].uncomp;
        cpos += idx[i].comp;
    }

    // rewind
    if (hseek(fp->hfp, 0, SEEK_SET) < 0)
        goto err;

    free(buf);
    return 0;

 err:
    free(buf);
    free(idx);
    return -1;
}

/*
 * Binary search the index for the first block that covers uncompressed
 * position upos.  If the index has pzstd skippable frames, we point to
 * that element.  This ensures we know the size of the next data container.
 *
 * NB: Index entries are non-overlapping.  So A to B, B+1 to C, C+1 to D.
 *
 * Returns an index entry on success,
 *         NULL on failure (errno==ERANGE if beyond end of file)
 */
static bzst_index_t *index_query(bzst *fp, uint64_t upos) {
    int istart = 0;
    int iend = fp->nindex-1;
    int imid;
    bzst_index_t *idx = fp->index;

    // Find approximate location
    for (imid = (iend+1)/2; imid != istart; imid = (istart+iend)/2) {
        if (idx[imid].pos >= upos)
            iend = imid;
        else
            istart = imid;
    }

    // We end with either imid or next (non-skip) imid being correct.
    // Select non-skippable frame
    while (imid+1 < fp->nindex && idx[imid].uncomp == 0)
        imid++;

    // We end with correct being idx[imid] or idx[imid+1]
    if (idx[imid].pos + idx[imid].uncomp <= upos) {
        // Select next non-skippable frame
        if (imid+1 < fp->nindex)
            imid++;
        while (imid+1 < fp->nindex && idx[imid].uncomp == 0)
            imid++;

        // Check for upos being end of index
        if (idx[imid].pos + idx[imid].uncomp <= upos) {
            errno = ERANGE;
            return NULL;
        }
    }

    // Finally walk back back to include skippable frames prior to this
    // offset, as these hold our meta-data.
    while (imid > 0 && idx[imid-1].uncomp == 0)
        imid--;

    return &idx[imid];
}

/*
 * Seeks to uncompressed position upos in a bgzf file opened for read.
 *
 * Returns 0 on success,
 *        -1 on failure
 *
 * TODO: consider "whence" and returning off_t, like lseek?
 */
int bzst_seek(bzst *fp, uint64_t upos) {
    int ret = 0;

    if (fp->pool) {
        // The multi-threaded reader runs asynchronously (fp->io_task).
        // We send it a command to do the seek instead od doing this ourselves.
        pthread_mutex_lock(&fp->command_m);
        fp->command = SEEK;
        fp->seek_to = upos;
        pthread_cond_signal(&fp->command_c);
        hts_tpool_wake_dispatch(fp->out_queue);

        // Now we wait for the io_task to give us the error code back.
        do {
            pthread_cond_wait(&fp->command_c, &fp->command_m);
            switch (fp->command) {
            case SEEK_FAIL:
                ret = -1;
                break;
            case SEEK_DONE:
                break;
            case SEEK:
                // Resend signal intended for bzst_mt_reader(), incase
                // we were woke up by something else.
                pthread_cond_signal(&fp->command_c);
                break;
            default:
                abort();  // Should not get to any other state
            }
        } while (fp->command != SEEK_DONE && fp->command != SEEK_FAIL);
        fp->command = NONE;
        pthread_mutex_unlock(&fp->command_m);

    } else {
        if (!fp->index) {
            int err;
            if ((err = load_seekable_index(fp)) < 0) {
                fp->errcode = -err;
                return -1;
            }
        }

        // Query in index
        bzst_index_t *idx = index_query(fp, upos);
        if (!idx)
            return -1;
        assert(upos >= idx->pos);

        if (idx->cpos != hseek(fp->hfp, idx->cpos, SEEK_SET))
            return -1;
        fprintf(stderr, "Seek to %ld\n", idx->cpos);

        // Load the relevant block
        if (bzst_decode_block(fp) < 0)
            return -1;

        // Skip past any partial data
        fp->uncomp->pos = upos - idx->pos;
    }

    return ret;
}

/*
 * Check for known EOF by detection of seekable index footer.
 *
 * Returns 0 if marker is absent,
 *         1 if present,
 *         2 if unable to check (eg cannot seek),
 *        -1 for I/O error, with errno set.
 */
int bzst_check_EOF(bzst *fp) {
    int has_eof, err = 0;

    if (fp->pool) {
        submit_reader_command(fp, HAS_EOF, HAS_EOF_DONE);
        pthread_mutex_lock(&fp->command_m);
        has_eof = fp->has_eof;
        err = fp->errcode;
        pthread_mutex_unlock(&fp->command_m);
    } else {
        has_eof = bzst_check_EOF_common(fp);
    }

    return err ? -1 : has_eof;
}

/**
 * Read one line from a BGZF file. It is faster than bgzf_getc()
 *
 * @param fp     BGZF file handler
 * @param delim  delimiter
 * @param str    string to write to; must be initialized
 * @return       length of the string (capped at INT_MAX);
 *               -1 on end-of-file; <= -2 on error
 */
int bzst_getline(bzst *fp, int delim, kstring_t *str) {
    int state = 0;
    str->l = 0;
    do {
        ssize_t n = 0;
        if ((n = bzst_refill_uncomp(fp)) < 0)
            return n; // EOF (-1) or error (-2)

        // Find end of line, or point to end of buffer if continues.
        // NB: this latter case shouldn't happen as we should always have
        // whole records in a bzst block, but this is to protect against
        // extreme data (eg >4GB line lengths) and just for extra safety.
        char *eol = memchr(fp->uncomp->buf + fp->uncomp->pos, delim,
                           fp->uncomp->sz - fp->uncomp->pos);
        if (eol)
            state = 1; // any +ve value; found delim
        else
            eol = fp->uncomp->buf + fp->uncomp->sz; // end of buffer
        n = eol - (fp->uncomp->buf + fp->uncomp->pos);

        // Copy the whole or partial line
        if (ks_expand(str, n + 2) < 0) {
            state = -3;
            break;
        }
        memcpy(str->s + str->l, fp->uncomp->buf + fp->uncomp->pos, n);
        fp->uncomp->pos += n+(state==1);
        str->l += n;
    } while (state == 0);

    if (state < -1) return state;
    if (str->l == 0 && state < 0) return state;

    if ( delim=='\n' && str->l>0 && str->s[str->l-1]=='\r' ) str->l--;
    str->s[str->l] = 0;
    return str->l <= INT_MAX ? (int) str->l : INT_MAX;
}

// -1 for EOF, -2 for error, 0-255 for byte.
int bzst_peek(bzst *fp) {
    ssize_t n = 0;
    if ((n = bzst_refill_uncomp(fp)) < 0)
        return n; // EOF (-1) or error (-2)

    return fp->uncomp->buf[fp->uncomp->pos];
}

/*
 * Adds a record to the genomic index.
 * Note this differs from the seekable index which is file offset based.
 *
 * As our index maps genomic coordinates to uncompressed offsets, and the
 * seekable index turns uncompressed offsets to compressed frame positions,
 * we can have multiple index items per frame if we choose, permitting us
 * to skip ahead.
 *
 * Minimally however we require 1 index entry item per chromosome, so this
 * is potentially multiple entries per frame.  Ideally we also require at
 * least index entry per frame to maintain efficiency.  More than one doesn't
 * harm, but is best done through checkpointing in-line with a frame
 * (appending to the pzstd preceeding frame) so that streaming filters are
 * also efficient.
 *
 * Returns 0 on success,
 *        <0 on failure
 */
int bzst_idx_add(bzst *fp, int tid, hts_pos_t beg, hts_pos_t end) {
    // Detect unsorted data
    if (fp->last_used < 0)
        return 0; // unsorted

    // FIXME: if we want to add meta-data like nreads, nmapped, nunmapped,
    // etc then we need a secondary table somewhere.

    if (fp->last_used > 0) {
        if (// unmapped not at end of file
            (fp->last_tid == -1 && tid != -1) ||
            // switching back to a tid we previously used
            (tid != fp->last_tid && tid+1 < fp->nchr && fp->gindex[tid+1]) ||
            // unsorted within a chromosome
            (tid == fp->last_tid && beg < fp->last_pos)) {
            //fprintf(stderr, "File is not sorted\n");
            fp->last_used = -1;
            return 0;
        }
    }
    fp->last_used = 1;
    fp->last_tid = tid;
    fp->last_pos = beg;

    tid++; // so unmapped becomes 0
    if (tid < 0)
        return -1;

    // chr22 10511189 rs1234474560 TTTCTTCCCAAATGTGTATTGATTACAC
    // => bzst_idx_add: 22 10511188 10511216
    //fprintf(stderr, "bzst_idx_add: %d %ld %ld\n", tid, beg, end);

    if (tid == 0)
        beg = end = 0;

    if (tid >= fp->nchr) {
        size_t *t = realloc(fp->gindex_sz, (tid+1)*sizeof(*t));
        if (!t)
            return -1;
        fp->gindex_sz = t;

        bzst_gindex_t **t2 = realloc(fp->gindex, (tid+1)*sizeof(*t2));
        if (!t)
            return -1;
        fp->gindex = t2;

        // Fill in any holes with a while loop, also skips over unmapped
        while (fp->nchr <= tid) {
            fp->gindex_sz[fp->nchr] = 0;
            fp->gindex[fp->nchr]    = NULL;
            fp->nchr++;
        }

        fp->tid_pos = fp->last_flush_try;
    }

    // Check for current entry, or allocate a new one
    bzst_gindex_t *idx = NULL;
    if (fp->gindex_sz[tid] > 0)
        idx = &fp->gindex[tid][fp->gindex_sz[tid]-1];

    if (!idx || fp->last_flush_try <= 0) {
        // A new frame, so add to the index
        bzst_gindex_t *t = realloc(fp->gindex[tid],
                                    (fp->gindex_sz[tid]+1) * sizeof(*t));
        if (!t)
            return -1;
        fp->gindex[tid] = t;

        idx = &fp->gindex[tid][fp->gindex_sz[tid]++];
        idx->frame_start = fp->idx_pos + fp->last_flush_try;
        idx->tid = tid-1;
        idx->beg = beg;
        idx->end = end;

        //fprintf(stderr, "bzst_idx_add: %d %ld %ld: %ld %ld\n", tid, beg, end,
        //      fp->frame_pos, fp->last_flush_try);
    }

    if (idx->beg > beg)
        idx->beg = beg;
    if (idx->end < end)
        idx->end = end;

    return 0;
}

/*
 * Returns the internal temporary kstring associated with this bzst fd.
 */
kstring_t *bzst_ks(bzst *fp) {
   return &fp->ks;
}


/*
 * Returns the next item from a bzst iterator.
 * Returns 0 if found,
 *        -1 at EOF,
 *      <=-2 for error.
 */
int bzst_itr_next(bzst *fp, hts_itr_t *iter, void *r, void *data) {
    if (!iter || iter->finished)
        return -1;

    if (!iter->is_bzst)
        return -2;

    if (iter->n_off) {
        iter->n_off = 0;
        // Check pos and cpos work correctly vs /tmp/_1.bcf
        if (bzst_seek(fp, iter->curr_off) < 0) {
            hts_log_error("Failed to seek to offset %"PRIu64"%s%s",
                          iter->curr_off,
                          errno ? ": " : "", strerror(errno));
            return -2;
        }
    }

    int ret, tid;
    hts_pos_t beg, end;

    // Ensure sort works so -1 is after all others.
    // Or, we could just cast to uint32_t before cmp and force wrap around.
    if (iter->tid == -1)
        iter->tid = INT_MAX;

    do {
        ret = iter->readrec((BGZF *)fp, data, r, &tid, &beg, &end);
        //fprintf(stderr, "bzst_itr_next: %d:%ld-%ld vs %d:%ld-%ld\n",
        //        tid, beg, end, iter->tid, iter->beg, iter->end);
        if (tid == -1)
            tid = INT_MAX;
        if (ret < 0) {
            iter->finished = 1;
            return ret; // EOF or error
        }
    } while (tid <= iter->tid && end <= iter->beg);

    if (tid != iter->tid || beg >= iter->end) {
        //fprintf(stderr, "finished\n");
        iter->finished = 1;
        return -1; // EOF
    }

    // These are part of the public API and things like
    // samtools view can use them in error messages.
    //fprintf(stderr, "found\n");
    iter->curr_tid = tid;
    iter->curr_beg = beg;
    iter->curr_end = end;

    return 0;
}

/*
 * Query a BZST genomic index on region
 *
 * Returns hts_itr_t struct ptr on success,
 *         NULL on failure
 */
hts_itr_t *bzst_itr_query(const hts_idx_t *idx,
			   int tid,
			   hts_pos_t beg,
			   hts_pos_t end,
			   hts_readrec_func *readrec) {
    const hts_bzst_idx_t *bidx = (const hts_bzst_idx_t *)idx;
    hts_itr_t *iter = (hts_itr_t *)calloc(1, sizeof(hts_itr_t));
    if (iter == NULL) return NULL;

    if (tid == HTS_IDX_NOCOOR)
        tid = -1;

    int64_t pos = bzst_query(bidx->fp, tid, beg, end);

    // hts_itr_t is public and extremely BGZF(1) specific.
    // We fill out tid, beg, end from arguments here, and we reuse
    // curr_off to hold the seek position with n_off=1.
    // On first use we seek there and set n_off=0 to indicate we've
    // moved.  We don't seek immediately as we wish to separate querying
    // from seeking.  (We may in theory query multiple places and then seek
    // back to them later on.)

    // An alternative to this would be to put all bzst specific region
    // functionality internal to the bzst fp, as was done with cram_fd.
    // This may be preferable long term for handling multi-region iterators.

    iter->is_bzst = 1;
    iter->tid = tid;
    iter->beg = tid == -1 ? -1 : beg;
    iter->end = tid == -1 ?  0 : end;
    iter->n_off = 1;
    iter->curr_off = pos;
    iter->readrec = readrec;

    return iter;
}
