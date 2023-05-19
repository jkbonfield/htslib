// TODO: have a way to set end coord so MT decode doesn't waste cycles.

// TODO: Make first block fit within 4k so hpeek can work.

// TODO: EOF block - can check via seekable index.

// NB: previous commit was considerably more efficient when overcomitting on
// threads vs CPU count.

/* The MIT License

   Copyright (C) 2023 Genome Research Ltd

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
BGZF2 is a Zstd compatible data file with random access support and
designed for parallel encoding and decoding.

It does this by combining zstd seekable format
(https://github.com/facebook/zstd/tree/dev/contrib/seekable_format)
with pzstd (https://github.com/facebook/zstd/tree/dev/contrib/pzstd).

A Zstd file is a series of frames. Zstd has the notion of data frames
holding compressed data, and skippable frames holding meta-data that
isn't part of the uncompressed output stream.  Both seekable and pzstd
have their own skippable frames.

Zstd format has no indication of the compressed size of a frame, so
it's hard to get this when reading a file.  Pzstd uses a skippable
frame to hold the compressed size of the next data frame, thus
permitting quick read-and-dispatch style decoding.

Seekable-zstd uses a skippable frame at the end of the file holding
the compressed and uncompressed sizes of every frame.  This permits
random access via a trivial binary search.  If we aren't streaming,
this would also be enough to permit parallel decoding, but we cannot
guarantee it.  Our seekable zstd index also includes our pzstd
skippable frames, but these are all listed with an uncompressed size
of zero.  They are required however to get the cumulative compressed
offsets correct.

Note the random access here is purely via uncompressed byte offsets.
The data format is agnostic to any content.  If we wish to store
genomic data sorted by chromosome and position, and to query by
genomic region, then we need an additional index (eg CSI or CRAI).
For now, this is an ancillary file as per existing .bam and vcf.gz.

TODO: consider also bringing in https://github.com/facebook/zstd/pull/2349
for in-stream dictionaries.  Ideally we'd want multiple dictionaries
and the ability to indicate which dictionary is associated with which
frame.  Maybe the most recent dictionary will always be the one to use?
Combining with parallel decoding this means holding multiple
dictionaries in memory and purging once we know all new blocks relate
to a newer dictionary.  It also means having a random access
capability, so we need to know in the seekable format which entries
are dictionaries.  Eg dsize == 0 and csize != 12

TODO: Add XXH64 checksums
*/

/*
Known skippable frame IDs:

0x184D2A50   pzstd, size of next frame
0x184D2A51   aruna footer
0x184D2A52   aruna footer
0x184D2A55   zpkglist - LZ4
0x184D2A56   zpkglist - LZ4
0x184D2A57   zpkglist - LZ4
0x184D2A5D   warc-zstd / dict-in-stream
0x184D2A5E   seekable

Suggests 0x184D2A5B for BGZF2 usage?  Maybe the genomic index, or some
other header meta-data.
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

#include "htslib/hfile.h"
#include "htslib/hts_endian.h"
#include "htslib/thread_pool.h"
#include "htslib/bgzf2.h"
#include "cram/pooled_alloc.h"

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif

// INTERNAL structure.  Do not use (consider moving to bgzf2.c and putting
// a blank one here)
typedef struct {
    off_t pos;     // cumulative uncompressed position prior to this
    size_t uncomp; // uncompressed size of this block
    size_t comp;   // compressed size of this block
    off_t cpos;    // cumulative compression poisition in file
} bgzf2_index_t;

// Note this is just a kstring with a linked list.
// We could perhaps embed the link into the ks->s pointer, so empty ones
// are linked together via their data itself.
typedef struct bgzf2_buffer {
    size_t alloc, sz, pos; // allocated, desired size, and curr pos.
    char *buf;
    struct bgzf2_buffer *next;
} bgzf2_buffer;

/*
 * Work-load for single compression or decompression job.
 */
typedef struct bgzf2_job {
    bgzf2 *fp;
    bgzf2_buffer *uncomp;
    bgzf2_buffer *comp;
    int errcode;
    int hit_eof; // view from main or view from thread?
    int job_num;
    int known_size;
    struct bgzf2_job *next;
} bgzf2_job;

// Message command types for communicating with bgzf2_mt_reader
enum mtaux_cmd {
    NONE = 0,
    SEEK,
    SEEK_DONE,
    SEEK_FAIL,
    HAS_EOF,
    HAS_EOF_DONE,
    CLOSE,
};

// INTERNAL structure.  Do not use (consider moving to bgzf2.c and putting
// a blank one here)
struct bgzf2 {
    // Shared with BGZF to provide an easy redirect for code that
    // doesn't want to have lots of conditionals.
    unsigned dummy1:16, is_zstd:1, first_block:1, dummy2:14;

    struct hFILE *hfp;    // actual file handle
    int format;           // encoding format (unused, but zlib, zstd, bsc, ...)
    int level;            // compression level
    int is_write;         // open for write
    int block_size;       // ideal block size
    bgzf2_index_t *index; // index entries
    int nindex;           // used size of index array
    int aindex;           // allocated size of index array
    int errcode;          // FIXME

    bgzf2_buffer *uncomp; // uncompressed data
    bgzf2_buffer *comp;   // compressed data

    // Multi-threading support
    pthread_mutex_t job_pool_m; // when updating this struct
    struct bgzf2_job *job_free_list;
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
};

/*----------------------------------------------------------------------
 * Utility functions shared by both single and multi-threaded implementations
 */

/*
 * Allocates a new bgzf2 buffer capable of holding n bytes.
 *
 * Returns buffer pointer on success,
 *         NULL on failure
 */
bgzf2_buffer *bgzf2_buffer_alloc(size_t n) {
    bgzf2_buffer *b = malloc(sizeof(*b));
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
int bgzf2_buffer_grow(bgzf2_buffer **bp, size_t n) {
    bgzf2_buffer *b = *bp;

    if (!*bp) {
	b = *bp = calloc(1, sizeof(*b));
	if (!b)
	    return -1;
    }
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

/* Frees memory used by a bgzf2_buffer */
void bgzf2_buffer_free(bgzf2_buffer *b) {
    if (!b)
	return;

    free(b->buf);
    free(b);
}

static pthread_once_t bgzf2_comp_once = PTHREAD_ONCE_INIT;
static pthread_key_t bgzf2_comp_key;

static void bgzf2_tls_comp_free(void *ptr) {
    ZSTD_CStream *zcs = (ZSTD_CStream *)ptr;
    if (zcs)
	ZSTD_freeCStream(zcs);
}

static void bgzf2_tls_comp_init(void) {
    pthread_key_create(&bgzf2_comp_key, bgzf2_tls_comp_free);
}

/*
 * Compresses a block of data.
 * The output buffer is dynamically grown and can be reused between calls.
 * Initially start with *comp = NULL and *comp_sz = 0.
 *
 * Returns compressed size on success,
 *        -1 on failure.
 */
static size_t compress_block(char *uncomp, size_t uncomp_sz,
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
    pthread_once(&bgzf2_comp_once, bgzf2_tls_comp_init);
    ZSTD_CStream *zcs = pthread_getspecific(bgzf2_comp_key);
    if (!zcs) {
	zcs = ZSTD_createCStream();
	pthread_setspecific(bgzf2_comp_key, zcs);

	if (ZSTD_initCStream(zcs, level) != 0) {
	    ZSTD_freeCStream(zcs);
	    return -1;
	}
    }
    ZSTD_CCtx_reset(zcs, ZSTD_reset_session_only);
    ZSTD_CCtx_setParameter(zcs, ZSTD_c_checksumFlag, 1);
    ZSTD_CCtx_setParameter(zcs, ZSTD_c_contentSizeFlag, 1);
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

// Helps on bigger buffer sizes (or higher compression levels?)
//        ZSTD_CCtx_setParameter(zcs, ZSTD_c_searchLog, 6);
//        ZSTD_CCtx_setParameter(zcs, ZSTD_c_minMatch, 6);
//        ZSTD_CCtx_setParameter(zcs, ZSTD_c_enableLongDistanceMatching, 1);
//        ZSTD_CCtx_setParameter(zcs, ZSTD_c_ldmBucketSizeLog, 4);
//        ZSTD_CCtx_setParameter(zcs, ZSTD_c_ldmHashRateLog, 7); // 4 at -3
	
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
static int write_seekable_index(bgzf2 *fp) {
    bgzf2_index_t *idx = fp->index;
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
 * Adds an entry to the bgzf2 index
 *
 * Returns 0 on success,
 *        -1 on failure.
 */
static int bgzf2_add_index(bgzf2 *fp, size_t uncomp, size_t comp) {
    bgzf2_index_t *idx;

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
 * Writes a pzstd compatible skippable frame indicating the size of
 * the next compressed data frame.  We also add it to the index we're
 * building up for zstd seekable.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
static int write_pzstd_skippable(bgzf2 *fp, uint32_t comp_sz) {
    uint8_t buf[12];

    u32_to_le(0x184D2A50, buf); // pzstd skippable magic no.
    u32_to_le(4, buf+4);
    u32_to_le(comp_sz, buf+8);

    int ret = bgzf2_add_index(fp, 0, 12);
    ret |= hwrite(fp->hfp, buf, 12) != 12;

    return ret ? -1 : 0;
}

/*----------------------------------------------------------------------
 * Multi-threaded implementation.
 *
 * Encode: the main thread copies data into a block (main thread).  When full
 * this calls bgzf2_flush to compress and write it out.
 * The bgzf2_flush function is where things fork, with bgzf2_write_block
 * performing these operations within the main thread and bgzf2_write_block_mt
 * creating an encode job to dispatch to the pool of workers.
 *
 * The thread pool workers all call bgzf2_encode_func, whose only task is
 * the data compression step.
 *
 * bgzf2_mt_writer is a dedicated I/O thread (not in the pool) whose sole
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
 * Decode: bgzf2_mt_reader is a dedicated I/O thread (not in the pool) that is
 * reading blocks ahead of time and dispatching decompression jobs to be
 * executed on the main thread pool.  These jobs call bgzf2_decode_func.
 *
 * We have a "command" channel for communicating with bgzf2_mt_read, for
 * purposes such as specifying end ranges, closing down, and handling
 * errors.
 */

/*
 * Doesn't explicitly free the memory, but instead puts it on a free list
 * so we can reuse it.
 */
static void bgzf2_job_free(void *vp) {
    bgzf2_job *j = (bgzf2_job *)vp;
    //fprintf(stderr, "bgzf2_job_free %d\n", j->job_num);

    bgzf2 *fp = j->fp;
    pthread_mutex_lock(&fp->job_pool_m);
    j->next = fp->job_free_list;
    fp->job_free_list = j;
    pthread_mutex_unlock(&fp->job_pool_m);
}

/*
 * Actually deallocates memory used by a job
 */
static void bgzf2_job_really_free(bgzf2_job *j) {
    bgzf2_buffer_free(j->comp);
    bgzf2_buffer_free(j->uncomp);
}

/*
 * Fetches a new bgzf2_job struct.  Reuses previously freed ones if available,
 * or allocates new if not.
 */
static bgzf2_job *bgzf2_job_new(bgzf2 *fp) {
    bgzf2_job *j;

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
static void *bgzf2_mt_writer(void *vp) {
    bgzf2 *fp = (bgzf2 *)vp;
    hts_tpool_result *r;

//    if (fp->idx_build_otf) {
//        fp->idx->moffs = fp->idx->noffs = 1;
//        fp->idx->offs = (bgzidx1_t*) calloc(fp->idx->moffs, sizeof(bgzidx1_t));
//        if (!fp->idx->offs) goto err;
//    }

    // Iterates until result queue is shutdown, where it returns NULL.
    while ((r = hts_tpool_next_result_wait(fp->out_queue))) {
        bgzf2_job *j = (bgzf2_job *)hts_tpool_result_data(r);
        assert(j);
//	fprintf(stderr, "job: %d %d->%d %.10s\n",
//		j->job_num, (int)j->uncomp->sz, (int)j->comp->sz,
//		j->uncomp->buf);
//	static int job_last = -1;
//	assert(j->job_num == job_last+1);
//	job_last++;

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


	if (write_pzstd_skippable(fp, j->comp->sz) < 0)
	    goto err;

	if (bgzf2_add_index(fp, j->uncomp->sz, j->comp->sz) < 0)
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
	bgzf2_job_free(j);
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
 * Function called by the bgzf2 encode multi-threaded worker tasks.
 */
static void *bgzf2_encode_func(void *vp) {
    bgzf2_job *j = (bgzf2_job *)vp;

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
 * Returns 0 on success
 *        <0 on failure (same codes as load_seekabel_index)
 */
int load_seekable_index(bgzf2 *fp);
bgzf2_index_t *index_query(bgzf2 *fp, uint64_t upos);
static void bgzf2_mt_seek(bgzf2 *fp) {
    if (!fp->index) {
	int err;
	if ((err = load_seekable_index(fp)) < 0) {
	    pthread_mutex_lock(&fp->job_pool_m);
	    fp->errcode = -err;
	    fp->command = SEEK_FAIL;
	    pthread_mutex_unlock(&fp->job_pool_m);
	    return;
	}
    }
    hts_tpool_process_reset(fp->out_queue, 0);

    pthread_mutex_lock(&fp->job_pool_m);

    bgzf2_index_t *idx = index_query(fp, fp->seek_to);
    if (!idx) {
	fp->errcode = 1; // FIXME
	fp->command = SEEK_FAIL;
    } else {
	if (hseek(fp->hfp, idx->cpos, SEEK_SET) < 0) {
	    fp->errcode = 99; // BGZF_ERR_IO FIXME
	    fp->command = SEEK_FAIL;
	} else {
	    fp->errcode = 0;
	    fp->command = SEEK_DONE;
	}

	// The block is loaded later, as before, but we need to start with
	// pos part way through it.  We do this by modifying seek_to, which
	// is used in bgzf2_decode_block_mt.
	fp->seek_to -= idx->pos;
    }
    fp->hit_eof = 0;

    pthread_mutex_unlock(&fp->job_pool_m);

    pthread_cond_signal(&fp->command_c);
}

static void bgzf2_mt_eof(bgzf2 *fp) {
    abort();
}

/*
 * Nul function so we can dispatch a job with the correct serial
 * to mark failure or to indicate an empty read (EOF).
 */
static void *bgzf2_nul_func(void *arg) {
    return arg;
}

/*
 * Single threaded, called by bgzf2_mt_reader.
 *
 * Reads a compressed block of data using hread and unwraps the
 * relevant headers to fill out a bgzf2_job prior to decompression.
 * This doesn't actually dispatch the decompression job.
 *
 * Returns >0 on success (size once uncompressed),
 *          INT_MAX if unknown size (not in frame header),
 *          0 on EOF,
 *         -1 on failure
 */
static int bgzf2_mt_read_block(bgzf2 *fp, bgzf2_job *j)
{
    size_t bgzf2_read_block(bgzf2 *fp, bgzf2_buffer **comp);

    ssize_t usize = bgzf2_read_block(fp, &j->comp);
    if (usize == -2)
	return INT_MAX;
    if (usize <= 0)
	return usize;

    // Allocate uncompressed data block
    if (usize != -2) {
	if (usize > BGZF2_MAX_BLOCK_SIZE)
	    return -1; // protect against extreme memory size attacks
	if (bgzf2_buffer_grow(&j->uncomp, usize) < 0)
	    return -1;
	j->known_size = 1;
    } else {
	j->known_size = 0;
    }

    j->fp = fp;
    j->errcode = 0;

    return usize;
}

/*
 * Threaded: runs in one of many thread pool workers.
 *
 * Function called by the bgzf2 decode multi-threaded worker tasks.
 */
static void *bgzf2_decode_func(void *vp) {
    bgzf2_job *j = (bgzf2_job *)vp;

    int bgzf2_decompress_block(bgzf2_buffer **comp, bgzf2_buffer **uncomp,
			       size_t usize);

    if (bgzf2_decompress_block(&j->comp, &j->uncomp,
			       j->known_size ? j->uncomp->sz : -2) < 0)
	j->errcode = 1;

//    if (ZSTD_decompress(j->uncomp->buf, j->uncomp->sz,
//			j->comp->buf, j->comp->sz) != j->uncomp->sz) {
//	j->errcode = 1; // TODO design err codes properly
//    }

    return vp;
}

/*
 * Single threaded, in a thread of its own (ie async with main).
 *
 * Reads compressed blocks from the file handle and dispatches jobs to the
 * thread pool to decompress them.
 */
static void *bgzf2_mt_reader(void *vp) {
    bgzf2 *fp = (bgzf2 *)vp;
    bgzf2_job *j;

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
    if (!(j = bgzf2_job_new(fp)))
	goto err;
    j->errcode = 0;
    j->hit_eof = 0;
    j->fp = fp;

    while (bgzf2_mt_read_block(fp, j) > 0) {
        // Dispatch
        if (hts_tpool_dispatch3(fp->pool, fp->out_queue, bgzf2_decode_func, j,
                                bgzf2_job_free, bgzf2_job_free, 0) < 0) {
            bgzf2_job_free(j);
            goto err;
        }

        // Check for command
        pthread_mutex_lock(&fp->command_m);
        switch (fp->command) {
        case SEEK:
	    // Sets fp->command to SEEK_DONE
            bgzf2_mt_seek(fp);
            pthread_mutex_unlock(&fp->command_m);
            goto restart;

        case HAS_EOF:
            bgzf2_mt_eof(fp);   // Sets fp->command to HAS_EOF_DONE
            break;

        case SEEK_DONE:
        case HAS_EOF_DONE:
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
	j = bgzf2_job_new(fp);
        if (!j) {
            hts_tpool_process_destroy(fp->out_queue);
            return NULL;
        }
	j->fp = fp;
    }

    if (j->errcode == 2 /* FIXME */) {
        // Attempt to multi-thread decode a raw gzip stream cannot be done.
        // We tear down the multi-threaded decoder and revert to the old code.
        if (hts_tpool_dispatch3(fp->pool, fp->out_queue, bgzf2_nul_func, j,
                                bgzf2_job_free, bgzf2_job_free, 0) < 0) {
            bgzf2_job_free(j);
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
    if (hts_tpool_dispatch3(fp->pool, fp->out_queue, bgzf2_nul_func, j,
                            bgzf2_job_free, bgzf2_job_free, 0) < 0) {
        bgzf2_job_free(j);
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
        default:
            pthread_mutex_unlock(&fp->command_m);
            break;

        case SEEK:
            bgzf2_mt_seek(fp);
            pthread_mutex_unlock(&fp->command_m);
            goto restart;

        case HAS_EOF:
            bgzf2_mt_eof(fp);   // Sets fp->command to HAS_EOF_DONE
            pthread_mutex_unlock(&fp->command_m);
            break;

        case SEEK_DONE:
        case HAS_EOF_DONE:
            pthread_cond_signal(&fp->command_c);
            pthread_mutex_unlock(&fp->command_m);
            break;

        case CLOSE:
            pthread_cond_signal(&fp->command_c);
            pthread_mutex_unlock(&fp->command_m);
            hts_tpool_process_destroy(fp->out_queue);
            return NULL;
        }
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
static int bgzf2_decode_block_mt(bgzf2 *fp) {
    if (fp->hit_eof)
	return 0;

    hts_tpool_result *r;

    r = hts_tpool_next_result_wait(fp->out_queue);
    if (!r)
	return -1;

    bgzf2_job *j = (bgzf2_job *)hts_tpool_result_data(r);
    hts_tpool_delete_result(r, 0);
    if (j->hit_eof) {
	fp->hit_eof = 1;
	bgzf2_job_free(j);
	return 0;
    }

    // Swap j->uncomp with fp->uncomp to avoid a memcpy
    bgzf2_buffer *tmp = fp->uncomp;
    fp->uncomp = j->uncomp;
    j->uncomp = tmp;

    bgzf2_job_free(j);

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
 * Main thread: adds a thread pool to a bgzf2 fp, adding
 * multi-threading capabilities.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int bgzf2_thread_pool(bgzf2 *fp, hts_tpool *pool, int qsize) {
    fp->own_pool = 0;
    fp->pool = pool;
    if (!qsize)
        qsize = hts_tpool_size(pool)*2;
    if (!(fp->out_queue = hts_tpool_process_init(fp->pool, qsize, 0)))
	return -1;
    hts_tpool_process_ref_incr(fp->out_queue);

    fp->job_pool = pool_create(sizeof(bgzf2_job));
    if (!fp->job_pool)
	return -1;

    pthread_mutex_init(&fp->job_pool_m, NULL);
    pthread_mutex_init(&fp->command_m, NULL);
//    pthread_mutex_init(&fp->idx_m, NULL);
    pthread_cond_init(&fp->command_c, NULL);
    fp->flush_pending = 0;
    fp->jobs_pending = 0;
    pthread_create(&fp->io_task, NULL,
                   fp->is_write ? bgzf2_mt_writer : bgzf2_mt_reader, fp);

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
static int bgzf2_write_block(bgzf2 *fp, bgzf2_buffer *buf) {
    if (bgzf2_buffer_grow(&fp->comp, ZSTD_compressBound(buf->sz)) < 0)
	return -1;

    if ((fp->comp->sz = compress_block(buf->buf, buf->pos,
				       fp->comp->buf, fp->comp->alloc,
				       fp->level)) < 0)
	return -1;

    if (write_pzstd_skippable(fp, fp->comp->sz) < 0)
	return -1;

    int ret = bgzf2_add_index(fp, buf->pos, fp->comp->sz);

    if (hwrite(fp->hfp, fp->comp->buf, fp->comp->sz) != fp->comp->sz)
	return -1;

    fp->uncomp->pos = 0; // reset buffered offset
    return ret;
}

/*
 * (Runs in the main thread)
 * Asynchronously compress and write a block to disk.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
static int bgzf2_write_block_mt(bgzf2 *fp, bgzf2_buffer *buf) {
    bgzf2_job *j = bgzf2_job_new(fp);
    if (!j)
	return -1;

    // Fill out job.
    size_t comp_sz = ZSTD_compressBound(buf->sz);
    j->fp = fp;
    j->job_num = fp->job_num++;
    j->next = NULL;

    // We could steal "buf" and optimise this more.
    if (bgzf2_buffer_grow(&j->uncomp, buf->sz) < 0 ||
	bgzf2_buffer_grow(&j->comp, comp_sz) < 0)
	return -1;
    memcpy(j->uncomp->buf, buf->buf, buf->pos);
    j->uncomp->pos = buf->pos;

    // Used by flush function to drain queue.
    // Can we get this via another route?
    pthread_mutex_lock(&fp->job_pool_m);
    fp->jobs_pending++;
    pthread_mutex_unlock(&fp->job_pool_m);

    if (hts_tpool_dispatch3(fp->pool, fp->out_queue, bgzf2_encode_func,
			    j, bgzf2_job_free, bgzf2_job_free, 0) < 0)
	goto err;

    fp->uncomp->pos = 0;
    return 0;

 err:
    bgzf2_job_free(j);
    pthread_mutex_lock(&fp->job_pool_m);
    fp->jobs_pending--;
    pthread_mutex_unlock(&fp->job_pool_m);
    return -1;
}

/*
 * Flush the bgzf2 stream and ensure we start a new block.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int bgzf2_flush(bgzf2 *fp) {
    int ret;
    if (!fp->uncomp->pos)
	return 0;

    if (fp->pool) {
	ret = bgzf2_write_block_mt(fp, fp->uncomp);
    } else {
	ret = bgzf2_write_block(fp, fp->uncomp);
    }

    if (fp->first_block) {
	fp->first_block = 0;
	ret = bgzf2_buffer_grow(&fp->uncomp, fp->block_size);
    }

    return ret;
}

/*
 * Flushes data and drains any asynchronous I/O still waiting to be
 * written.
 *
 * Returns 0 on success,
 *        <0 on failure
 */
int bgzf2_drain(bgzf2 *fp) {
    if (bgzf2_flush(fp) < 0)
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
 * Set the bgzf2 block size.  This can be performed at any point,
 * but it is usually done immediately after opening for write.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int bgzf2_set_block_size(bgzf2 *fp, size_t sz) {
    fp->block_size = sz;

    if (fp->uncomp)
	if (bgzf2_flush(fp) < 0)
	    return -1;

    if (fp->first_block)
	sz = sz < 1024 ? sz : 1024;

    return bgzf2_buffer_grow(&fp->uncomp, sz);
}

static bgzf2 *bgzf2_open_common(bgzf2 *fp, hFILE *hfp, const char *mode) { 
    fp->is_zstd = 1;
    fp->first_block = 1;
    fp->hfp = hfp;

    if (*mode == 'w') {
	fp->is_write = 1;
	if (bgzf2_set_block_size(fp, BGZF2_DEFAULT_BLOCK_SIZE) < 0) {
	    if (hclose(fp->hfp)) // can't ignore :/
		free(fp);
	    else
		free(fp);
	    return NULL;
	}

	int level = BGZF2_DEFAULT_LEVEL;
	if (mode[1] >= '0' && mode[1] <= '9')
	    level = atoi(mode+1);

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
//    if (bgzf2_thread_pool(fp, fp->pool, 0) < 0)
//	goto err;

    return fp;
}

/*
 * Opens a bgzf2 file from an existing hfile for read ("r") or write
 * ("w", or "w1" to "w19").
 *
 * Returns bgzf2 file handle on success,
 *         or NULL on failure.
 */
bgzf2 *bgzf2_hopen(hFILE *hfp, const char *mode) {
    bgzf2 *fp = calloc(1, sizeof(*fp));
    if (!fp)
	return NULL;

    return bgzf2_open_common(fp, hfp, mode);
}
/*
 * Opens a bgzf2 file 'fn' for read ("r") or write ("w", or "w1" to "w19").
 *
 * Returns bgzf2 file handle on success,
 *         or NULL on failure.
 */
bgzf2 *bgzf2_open(const char *fn, const char *mode) {
    bgzf2 *fp = calloc(1, sizeof(*fp));
    if (!fp)
	return NULL;

    fp->hfp = hopen(fn, mode);
    if (!fp->hfp) {
	free(fp);
	return NULL;
    }

    return bgzf2_open_common(fp, fp->hfp, mode);
}


/*
 * Closes a bgzf2 file.
 *
 * Returns 0 on success,
 *        <0 on failure
 */
int bgzf2_close(bgzf2 *fp) {
    int ret = 0;

    if (fp->is_write) {
	ret |= (bgzf2_drain(fp) < 0);
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

	bgzf2_job *j, *n;
	for (j = fp->job_free_list; j; j = n) {
	    //fprintf(stderr, "free job %d\n", j->job_num);
	    n = j->next;
	    bgzf2_job_really_free(j);
	}

	pool_destroy(fp->job_pool);
	if (fp->own_pool)
	    hts_tpool_destroy(fp->pool);
    }

    if (fp->comp)   bgzf2_buffer_free(fp->comp);
    if (fp->uncomp) bgzf2_buffer_free(fp->uncomp);

    ret |= hclose(fp->hfp);

    free(fp->index);
    free(fp);

    return ret ? -1 : 0;
}

/*
 * Writes a block of data to a bgzf2 file handle.
 *
 * If "can_split" is set then buf may be split into two indexable chunks.
 * Otherwise buf will never span two blocks.
 *
 * Returns number of bytes written on success
 *        -1 on failure
 */
int bgzf2_write(bgzf2 *fp, char *buf, size_t buf_sz, int can_split) {
    int r = 0;

    do {
	if (fp->uncomp && fp->uncomp->sz == fp->uncomp->pos) {
	    // Full
	    if (bgzf2_flush(fp))
		return -1;
	}

	if (!fp->uncomp)
	    if (bgzf2_buffer_grow(&fp->uncomp, fp->block_size) < 0)
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
	    int err = bgzf2_flush(fp);
	    if (err || bgzf2_buffer_grow(&fp->uncomp, buf_sz) < 0)
		return -1;
	    memcpy(fp->uncomp->buf, buf, buf_sz);

	    // else it will fit on the next loop given we've now flushed
	}
    } while (buf_sz);

    return 0;
}

/*
 * Reads more compressed data into the bgzf2_buffer.
 * NB: The return value is the size of the data when uncompressed.
 * The compressed size will be in (*comp)->sz.
 *
 * Returns >0 on success (size when uncompressed),
 *          0 on EOF,
 *         -1 on failure
 *         -2 on unknown sized blocks (eg pzstd format)
 *         -3 on switch to stream mode
 */
/*static*/ size_t bgzf2_read_block(bgzf2 *fp, bgzf2_buffer **comp) {
    uint8_t buf[12];
    size_t n;

 next_block:
    if (8 != (n = hread(fp->hfp, buf, 8)))
	return n == 0 ? 0 : -1;

    uint32_t fsize = le_to_u32(buf+4);

    // Check if it's a pzstd format frame.
    if (le_to_u32(buf) != 0x184D2A50 || le_to_u32(buf+4) != 4) {
	if (le_to_u32(buf) >= 0x184D2A50 && le_to_u32(buf) <= 0x184D2A5F) {
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
	// remainder of pzstd skippable frame, now we know it is one
	if (4 != hread(fp->hfp, buf+8, 4))
	    return -1;
    }

    // Load compressed data
    size_t csize = le_to_u32(buf+8);
    if (bgzf2_buffer_grow(comp, csize) < 0)
	return -1;
    if (csize != hread(fp->hfp, (*comp)->buf, csize))
	return -1;

    // Get decompressed size and return it.
    size_t usize = ZSTD_getFrameContentSize((*comp)->buf, csize);
    if (usize == ZSTD_CONTENTSIZE_UNKNOWN)
	return -2;
    if (usize == ZSTD_CONTENTSIZE_ERROR)
	return -1;
    if (usize == 0)
	return 0; // empty frame => skip to next

    return usize;
}

static pthread_once_t bgzf2_once = PTHREAD_ONCE_INIT;
static pthread_key_t bgzf2_key;

static void bgzf2_tls_free(void *ptr) {
    ZSTD_DStream *zds = (ZSTD_DStream *)ptr;
    if (zds)
	ZSTD_freeDStream(zds);
}

static void bgzf2_tls_init(void) {
    pthread_key_create(&bgzf2_key, bgzf2_tls_free);
}

/*static*/
int bgzf2_decompress_block(bgzf2_buffer **comp, bgzf2_buffer **uncomp,
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
	    if (bgzf2_buffer_grow(uncomp, usize) < 0)
		return -1;
	}
	if (!*uncomp || (*uncomp)->alloc == 0)
	    // first cautious guess
	    if (bgzf2_buffer_grow(uncomp, (*comp)->sz*4 + 8192) < 0)
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
	pthread_once(&bgzf2_once, bgzf2_tls_init);
	ZSTD_DStream *zds = pthread_getspecific(bgzf2_key);
	if (!zds) {
	    zds = ZSTD_createDStream();
	    pthread_setspecific(bgzf2_key, zds);
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
		if (bgzf2_buffer_grow(uncomp, guess) < 0)
		    return -1;

		output.dst  = (*uncomp)->buf;
		output.size = (*uncomp)->alloc;
	    }
	}

	// We've used all input, but may still be blocked on output size.
	while (dsize != 0) {
	    fprintf(stderr, "Second zstd decode\n");
	    if ((*uncomp)->pos == (*uncomp)->sz)
		if (bgzf2_buffer_grow(uncomp, (*uncomp)->alloc*1.5+100000) < 0)
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
	if (usize > BGZF2_MAX_BLOCK_SIZE)
	    return -1; // protect against extreme memory size attacks
	if (bgzf2_buffer_grow(uncomp, usize) < 0)
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
static int bgzf2_decode_block(bgzf2 *fp) {
    ssize_t usize = bgzf2_read_block(fp, &fp->comp);
    if (usize != -2 && usize <= 0)
	return usize;

    return bgzf2_decompress_block(&fp->comp, &fp->uncomp, usize);

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
	    if (bgzf2_buffer_grow(&fp->uncomp, fp->comp->sz*2) < 0)
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
		if (bgzf2_buffer_grow(&fp->uncomp, guess) < 0)
		    return -1;

		output.dst  = fp->uncomp->buf;
		output.size = fp->uncomp->sz;
	    }
	}

	// We've used all input, but may still be blocked on output size.
	while (dsize != 0) {
	    if (bgzf2_buffer_grow(&fp->uncomp,
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
	if (usize > BGZF2_MAX_BLOCK_SIZE)
	    return -1; // protect against extreme memory size attacks
	if (bgzf2_buffer_grow(&fp->uncomp, usize) < 0)
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
 * Reads a data from a bgzf2 file handle.
 *
 * Returns number of bytes read on success
 *        -1 on failure
 */
int bgzf2_read(bgzf2 *fp, char *buf, size_t buf_sz) {
    if (fp->hit_eof)
	return 0;

    size_t decoded = 0;

    // loop:
    //    fill buffer if empty
    //        read a block to comp, decompress it, store in buffer
    //    copy from buffer to buf
    while (buf_sz) {
	if (!fp->uncomp || fp->uncomp->pos == fp->uncomp->sz) {
	    // out of buffered content, fetch some more
	    size_t n = fp->pool
		? bgzf2_decode_block_mt(fp)
		: bgzf2_decode_block(fp);
	    if (n < 0)
		return -1;
	    else if (n == 0) {
		fp->hit_eof = 1;
		return decoded; // EOF
	    }
	}

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
 * Reads a data from a bgzf2 file handle.  This modifies *buf to
 * point to a block of internal data and returns the size of this data.
 * In will be between 1 and buf_sz bytes long.  This data should not be
 * modified.
 *
 * Returns number of bytes read on success
 *        -1 on failure
 */
int bgzf2_read_zero_copy(bgzf2 *fp, const char **buf, size_t buf_sz) {
    if (fp->hit_eof)
	return 0;

    size_t decoded = 0;
    *buf = NULL;

    if (!buf_sz)
	return 0;

    if (!fp->uncomp || fp->uncomp->pos == fp->uncomp->sz) {
	// out of buffered content, fetch some more
	ssize_t n = fp->pool
	    ? bgzf2_decode_block_mt(fp)
	    : bgzf2_decode_block(fp);
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
 * Loads a seekable index from a bgzf2 file.
 *
 * Returns 0 on success,
 *        -1 on error,
 *        -2 on non-seekable stream.
 *        -3 if no index found.
 */
int load_seekable_index(bgzf2 *fp) {
    // Look for and validate index footer
    off_t off = hseek(fp->hfp, -9, SEEK_END);
    if (off < 0)
	return -1 - (errno == ESPIPE);

    uint8_t footer[9];
    if (9 != hread(fp->hfp, footer, 9))
	return -1;

    if (le_to_u32(footer+5) != 0x8F92EAB1 || (footer[4] & 0x7C) != 0)
	return -3;
    int has_chksum = footer[4] & 0x80;

    // Read entire index frame
    uint32_t nframes = le_to_u32(footer);
    size_t sz = 9 + nframes*4*(2+has_chksum) + 8;
    off = hseek(fp->hfp, -sz, SEEK_END);
    if (off < 0)
	return -1;

    bgzf2_index_t *idx = NULL;
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
bgzf2_index_t *index_query(bgzf2 *fp, uint64_t upos) {
    int istart = 0;
    int iend = fp->nindex-1;
    int imid;
    bgzf2_index_t *idx = fp->index;

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
int bgzf2_seek(bgzf2 *fp, uint64_t upos) {
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
                // Resend signal intended for bgzf2_mt_reader(), incase
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
	bgzf2_index_t *idx = index_query(fp, upos);
	if (!idx)
	    return -1;
	assert(upos >= idx->pos);

	if (idx->cpos != hseek(fp->hfp, idx->cpos, SEEK_SET))
	    return -1;

	// Load the relevant block
	if (bgzf2_decode_block(fp) < 0)
	    return -1;

	// Skip past any partial data
	fp->uncomp->pos = upos - idx->pos;
    }

    return ret;
}

