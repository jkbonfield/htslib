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

// Ensure ssize_t exists within this header. All #includes must precede this,
// and ssize_t must be undefined again at the end of this header.
#if defined _MSC_VER && defined _INTPTR_T_DEFINED && !defined _SSIZE_T_DEFINED && !defined ssize_t
#define HTSLIB_SSIZE_T
#define ssize_t intptr_t
#endif

#ifdef __cplusplus
extern "C" {
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
typedef struct {
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
} bgzf2;

#define BGZF2_DEFAULT_BLOCK_SIZE 256000
#define BGZF2_DEFAULT_LEVEL 5
#define BGZF2_MAX_BLOCK_SIZE (1<<30)

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
int bgzf2_write(bgzf2 *fp, char *buf, size_t buf_sz, int can_split);

/*
 * Reads a block of data from a bgzf2 file handle.
 *
 * Returns number of bytes read on success
 *        -1 on failure
 */
int bgzf2_read(bgzf2 *fp, char *buf, size_t buf_sz);

/*
 * Flush the bgzf2 stream and ensure we start a new block.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int bgzf2_flush(bgzf2 *fp);

/*
 * Loads a seekable index from a bgzf2 file.
 *
 * Returns 0 on success,
 *        -1 on error,
 *        -2 on non-seekable stream.
 *        -3 if no index found.
 */
int load_seekable_index(bgzf2 *fp);

int bgzf2_seek(bgzf2 *fp, uint64_t upos);

#ifdef __cplusplus
}
#endif

#ifdef HTSLIB_SSIZE_T
#undef HTSLIB_SSIZE_T
#undef ssize_t
#endif

#endif /* HTSLIB_BGZF2_H */
