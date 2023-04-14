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
#include "htslib/bgzf2.h"

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

// INTERNAL structure.  Do not use (consider moving to bgzf2.c and putting
// a blank one here)
struct bgzf2 {
    struct hFILE *hfp;    // actual file handle
    int format;           // encoding format (unused, but zlib, zstd, bsc, ...)
    int level;            // compression level
    int is_write;         // open for write
    int block_size;       // ideal block size
    bgzf2_index_t *index; // index entries
    int nindex;           // used size of index array
    int aindex;           // allocated size of index array

    char *buffer;         // uncompressed data
    size_t buffer_sz;     // used size of buffer
    size_t buffer_alloc;  // allocated size of buffer
    size_t buffer_pos;    // index into current buffer
    char *comp;           // compressed data block
    size_t comp_sz;       // size of compressed data
    size_t comp_alloc;    // allocated size of compressed block
};

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

    return off == hwrite(fp->hfp, buf, off) ? 0 : -1;
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
			     char **comp, size_t *comp_alloc,
			     int level) {
    // Grow output buffer
    size_t comp_bound = ZSTD_compressBound(uncomp_sz);
    if (comp_bound > *comp_alloc) {
	char *c = realloc(*comp, comp_bound);
	if (!c)
	    return -1;
	*comp = c;
	*comp_alloc = comp_bound;
    }

    // Use zstd.
    // For now we create and configure new streams each time, but we
    // could consider reusing the same stream.
    ZSTD_CStream *zcs = ZSTD_createCStream();
    if (!zcs)
	return -1;

    if (ZSTD_initCStream(zcs, level) != 0) {
	ZSTD_freeCStream(zcs);
	return -1;
    }

// Helps on bigger buffer sizes
//        ZSTD_CCtx_setParameter(zcs, ZSTD_c_searchLog, 6);
//        ZSTD_CCtx_setParameter(zcs, ZSTD_c_minMatch, 6);
//        ZSTD_CCtx_setParameter(zcs, ZSTD_c_enableLongDistanceMatching, 1);
//        ZSTD_CCtx_setParameter(zcs, ZSTD_c_ldmBucketSizeLog, 4);
//        ZSTD_CCtx_setParameter(zcs, ZSTD_c_ldmHashRateLog, 7); // 4 at -3
	
    size_t csize = ZSTD_compress2(zcs, *comp, *comp_alloc, uncomp, uncomp_sz);
    ZSTD_freeCStream(zcs);

    return ZSTD_isError(csize) ? -1 : csize;
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

/*
 * Compress and write a block to disk.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
static int bgzf2_write_block(bgzf2 *fp, char *buf, size_t buf_sz) {
    if ((fp->comp_sz = compress_block(buf, buf_sz, &fp->comp, &fp->comp_alloc,
				      fp->level)) < 0)
	return -1;

    if (write_pzstd_skippable(fp, fp->comp_sz) < 0)
	return -1;

    int ret = bgzf2_add_index(fp, buf_sz, fp->comp_sz);

    if (hwrite(fp->hfp, fp->comp, fp->comp_sz) != fp->comp_sz)
	return -1;

    fp->buffer_pos = 0; // reset buffered offset
    return ret;
}

/*
 * Flush the bgzf2 stream and ensure we start a new block.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int bgzf2_flush(bgzf2 *fp) {
    if (!fp->buffer_pos)
	return 0;

    return bgzf2_write_block(fp, fp->buffer, fp->buffer_pos);
}

/*
 * Set the bgzf2 block size.  This can be performed at any point,
 * but it is usually done immediately after opening for write.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int bgzf2_set_block_size(bgzf2 *fp, size_t sz) {
    if (bgzf2_flush(fp) < 0)
	return -1;

    if (fp->buffer_alloc < sz) {
	char *b = realloc(fp->buffer, sz);
	if (!b)
	    return -1;
	fp->buffer = b;
	fp->buffer_alloc = fp->buffer_sz = sz;
    }

    return 0;
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

    return fp;
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
	ret |= (bgzf2_flush(fp) < 0);
	ret |= write_seekable_index(fp);
    }

    ret |= hclose(fp->hfp);

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
	if (fp->buffer_sz == fp->buffer_pos) {
	    // Full
	    if (bgzf2_flush(fp))
		return -1;
	}

	ssize_t consumes = fp->buffer_sz - fp->buffer_pos;
	if (consumes > buf_sz)
	    consumes = buf_sz;

	if (consumes == buf_sz || can_split) {
	    // Whole or can_split and a portion
	    memcpy(fp->buffer + fp->buffer_pos, buf, consumes);
	    fp->buffer_pos += consumes;
	    buf_sz -= consumes;
	    buf    += consumes;
	    r      += consumes;
	} else {
	    // Can't split and doesn't fit, so flush and write this new item
	    // as a new block if it's too large.
	    int err = bgzf2_flush(fp);
	    if (buf_sz >= fp->buffer_alloc) {
		err |= bgzf2_write_block(fp, buf, buf_sz) < 0;

		if (err)
		    return r ? r : -1;
		buf_sz = 0;
	    }
	    // else it will fit on the next loop given we've now flushed
	}
    } while (buf_sz);

    return 0;
}

/*
 * Reads more compressed data and decompresses it.
 * This can work on both non-random access zstd compressed files (TODO)
 * as well as block compressed ones.
 *
 * This fills out fp->buffer and resets to fp->buffer_pos/sz to be
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
static int bgzf2_read_block(bgzf2 *fp) {
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
	return -2;
    } else {
	// remainder of pzstd skippable frame, now we know it is one
	if (4 != hread(fp->hfp, buf+8, 4))
	    return -1;
    }

    // Load compressed data
    fp->comp_sz = le_to_u32(buf+8);
    if (fp->comp_alloc < fp->comp_sz) {
	fp->comp_alloc = fp->comp_sz;
	char *tmp = realloc(fp->comp, fp->comp_alloc);
	if (!tmp)
	    return -1;
	fp->comp = tmp;
    }
    if (fp->comp_sz != hread(fp->hfp, fp->comp, fp->comp_sz))
	return -1;

    // Get decompressed size and allocate
    size_t usize = ZSTD_getFrameContentSize(fp->comp, fp->comp_sz);
    if (usize == ZSTD_CONTENTSIZE_UNKNOWN)
	return -2;
    if (usize == ZSTD_CONTENTSIZE_ERROR)
	return -1;
    if (usize == 0)
	return 0; // empty frame => skip to next

    if (usize > BGZF2_MAX_BLOCK_SIZE)
	return -1; // protect against extreme memory size attacks
    if (fp->buffer_alloc < usize) {
	fp->buffer_alloc = usize;
	char *tmp = realloc(fp->buffer, usize);
	if (!tmp)
	    return -1;
	fp->buffer = tmp;
    }

    // Uncompress
    if (ZSTD_decompress(fp->buffer, usize, fp->comp, fp->comp_sz) != usize) {
	// Embedded size mismatches the stream
	return -1;
    }

    fp->buffer_sz = usize;
    fp->buffer_pos = 0;

    return usize;
}

/*
 * Reads a data from a bgzf2 file handle.
 *
 * Returns number of bytes read on success
 *        -1 on failure
 */
int bgzf2_read(bgzf2 *fp, char *buf, size_t buf_sz) {
    size_t decoded = 0;

    // loop:
    //    fill buffer if empty
    //        read a block to comp, decompress it, store in buffer
    //    copy from buffer to buf
    while (buf_sz) {
	if (fp->buffer_pos == fp->buffer_sz) {
	    // out of buffered content, fetch some more
	    size_t n = bgzf2_read_block(fp);
	    if (n < 0)
		return -1;
	    else if (n == 0)
		return decoded; // EOF
	}

	size_t n = MIN(buf_sz, fp->buffer_sz - fp->buffer_pos);
	memcpy(buf, fp->buffer + fp->buffer_pos, n);
	buf += n;
	buf_sz -= n;
	fp->buffer_pos += n;
	decoded += n;
    }

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
 * Returns an index entry on success,
 *         NULL on failure (errno==ERANGE if beyond end of file)
 */
bgzf2_index_t *index_query(bgzf2 *fp, uint64_t upos) {
    int idx_start = 0;
    int idx_end = fp->nindex;
    bgzf2_index_t *idx = fp->index;

    while (idx_start < idx_end) {
	int mid = (idx_start + idx_end) / 2;
	if (upos < idx[mid].pos)
	    idx_end = mid;
	else if (upos >= idx[mid].pos + idx[mid].uncomp)
	    idx_start = mid;
	else {
	    idx_start = mid;
	    break;
	}
    }

//    fprintf(stderr, "Found %ld at %d: %ld + %ld\n",
//	    (long)upos, idx_start, idx[idx_start].pos, idx[idx_start].uncomp);

    // Skip back including skippable frames
    while (idx_start > 0 && idx[idx_start-1].uncomp == 0)
	idx_start--;

    return &idx[idx_start];
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
    // Query in index
    bgzf2_index_t *idx = index_query(fp, upos);
    if (!idx)
	return -1;
    assert(upos >= idx->pos);
    //assert(upos < idx->pos + idx->uncomp); // can fail due to skippable

    // Load the relevant block
    if (idx->cpos != hseek(fp->hfp, idx->cpos, SEEK_SET))
	return -1;

    if (bgzf2_read_block(fp) < 0)
	return -1;

    // Skip past any partial data
    fp->buffer_pos = upos - idx->pos;

    return 0;
}

