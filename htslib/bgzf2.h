/// @file htslib/bgzf2.h
/// Low-level routines for direct BGZF2 operations.
/*
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

/* BGZF2 is a derivative of BGZF that has a variable block size and uses
 * more modern compression codecs.
 */

/*
 * For now we only support Seekable-Zstd format.  This is compatible with the
 * same format in the zstd source contrib directory, but with some additional
 * headers (TODO) in a skippable frame.
 */

#ifndef HTSLIB_BGZF2_H
#define HTSLIB_BGZF2_H

#include <stdint.h>
#include <sys/types.h>

#include "hts_defs.h"
#include "thread_pool.h"
#include "kstring.h"
#include "hfile.h"

// Ensure ssize_t exists within this header. All #includes must precede this,
// and ssize_t must be undefined again at the end of this header.
#if defined _MSC_VER && defined _INTPTR_T_DEFINED && !defined _SSIZE_T_DEFINED && !defined ssize_t
#define HTSLIB_SSIZE_T
#define ssize_t intptr_t
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct bgzf2 bgzf2;

#include "hts.h"

#define BGZF2_DEFAULT_BLOCK_SIZE 256000
#define BGZF2_DEFAULT_LEVEL 5
#define BGZF2_MAX_BLOCK_SIZE (1<<30)

/*
 * Opens a bgzf2 file from an existing hfile for read ("r") or write
 * ("w", or "w1" to "w19").
 *
 * Returns bgzf2 file handle on success,
 *         or NULL on failure.
 */
bgzf2 *bgzf2_hopen(hFILE *hfp, const char *mode);

/*
 * Opens a bgzf2 file 'fn' for read ("r") or write ("w", or "w1" to "w19").
 *
 * Returns bgzf2 file handle on success,
 *         or NULL on failure.
 */
bgzf2 *bgzf2_open(const char *fn, const char *mode);

/*
 * Closes a bgzf2 file.
 *
 * Returns 0 on success,
 *        <0 on failure
 */
int bgzf2_close(bgzf2 *fp);

/*
 * Set the bgzf2 block size.  This can be performed at any point,
 * but it is usually done immediately after opening for write.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int bgzf2_set_block_size(bgzf2 *fp, size_t sz);

/*
 * Writes a block of data to a bgzf2 file handle.
 *
 * If "can_split" is set then buf may be split into two indexable chunks.
 * Otherwise buf will never span two blocks.
 *
 * Returns number of bytes written on success
 *        -1 on failure
 */
int bgzf2_write(bgzf2 *fp, const char *buf, size_t buf_sz, int can_split);

/*
 * Reads a block of data from a bgzf2 file handle.
 *
 * Returns number of bytes read on success
 *        -1 on failure
 */
int bgzf2_read(bgzf2 *fp, char *buf, size_t buf_sz);

/*
 * Reads a data from a bgzf2 file handle.  This modifies *buf to
 * point to a block of internal data and returns the size of this data.
 * In will be between 1 and buf_sz bytes long.  This data should not be
 * modified.
 *
 * Returns number of bytes read on success
 *        -1 on failure
 */
int bgzf2_read_zero_copy(bgzf2 *fp, const char **buf, size_t buf_sz);

/*
 * Flush the bgzf2 stream and ensure we start a new block.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int bgzf2_flush(bgzf2 *fp);

/*
 * Tests whether a write of 'size' would spill over to the next block. If so
 * flush this current one, so we always end blocks on a whole record.
 *
 * Returns 0 on success,
 *        <0 on failure
 */
int bgzf2_flush_try(bgzf2 *fp, ssize_t size);

/*
 * Finds the uncompressed file offset associated with a specific range.
 * This offset is only suitable for passing into bgzf2_seek.
 *
 * Note will not likely be the exact location the first record covers this
 * region, but it will be prior to it.  The caller is expected to them
 * discard data out of range.
 *
 * TODO: or should we cache beg/end here and do the discard ourselves?
 * That's how CRAM's API works, and it may be more amenable to a
 * "bgzf2_query_many" func that can operate on multiple regions in a
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
int64_t bgzf2_query(bgzf2 *fp, int tid, hts_pos_t beg, hts_pos_t end);

/*
 * Seeks to uncompressed position upos in a bgzf file opened for read.
 *
 * Returns 0 on success,
 *        -1 on failure
 *
 * TODO: consider "whence" and returning off_t, like lseek?
 */
int bgzf2_seek(bgzf2 *fp, uint64_t upos);

/*
 * Check for known EOF by detection of seekable index footer.
 *
 * Returns 0 if marker is absent,
 *         1 if present,
 *         2 if unable to check (eg cannot seek),
 *        -1 for I/O error, with errno set.
 */
int bgzf2_check_EOF(bgzf2 *fp);

/**
 * Enable multi-threading via a shared thread pool.  This means
 * both encoder and decoder can balance usage across a single pool
 * of worker jobs.
 *
 * @param fp          BGZF file handler
 * @param pool        The thread pool (see hts_create_threads)
 * @param qsize       Size of job queue, 0 for auto
 */
int bgzf2_thread_pool(bgzf2 *fp, hts_tpool *pool, int qsize);

/**
 * Read one line from a BGZF file. It is faster than bgzf_getc()
 *
 * @param fp     BGZF file handler
 * @param delim  delimiter
 * @param str    string to write to; must be initialized
 * @return       length of the string (capped at INT_MAX);
 *               -1 on end-of-file; <= -2 on error
 */
int bgzf2_getline(bgzf2 *fp, int delim, kstring_t *str);

/**
 * Returns the next byte in the file without consuming it.
 * @param fp     BGZF file handler
 * @return       -1 on EOF,
 *               -2 on error,
 *               otherwise the unsigned byte value.
 */
int bgzf2_peek(bgzf2 *fp);

/*
 * Adds a record to the index.
 * Returns 0 on success,
 *        <0 on failure
 */
int bgzf2_idx_add(bgzf2 *fp, int tid, hts_pos_t beg, hts_pos_t end);

#ifdef __cplusplus
}
#endif

#ifdef HTSLIB_SSIZE_T
#undef HTSLIB_SSIZE_T
#undef ssize_t
#endif

#endif /* HTSLIB_BGZF2_H */
