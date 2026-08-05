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

#include <stdio.h>
#include <stdlib.h>
#include <zstd.h>
#include <inttypes.h>
#include <getopt.h>
#include <unistd.h>
#include <errno.h>

#include "htslib/bzst.h"
#include "htslib/hfile.h"
#include "htslib/thread_pool.h"
#include "htslib/hts_endian.h"

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#define BUFSZ 5000000
static char buff_static[BUFSZ];

enum mode_t {
    MODE_AUTO,
    MODE_TEXT,
    MODE_BIN
};

/* ------------------------------------------------------------------------
 * Compression and decompression.
 */
static int encode(char *in, char *out, int level, size_t block_size,
                  int nthreads, enum mode_t mode, int rec_mod, char hdr_chr) {
    hFILE *fp_in = NULL;
    bzst *fp_out = NULL;
    int ret = -1;
    hts_tpool *pool = NULL;
    char *buffer = buff_static;

    fp_in = hopen(in, "r");
    if (!fp_in) goto err;

    char omode[30];
    snprintf(omode, 30, "w%d", level);
    fp_out = bzst_open(out, omode);
    if (!fp_out) goto err;

    if (nthreads) {
        pool = hts_tpool_init(nthreads);
        if (!pool)
            goto err;
        if (bzst_thread_pool(fp_out, pool, 0) < 0)
            goto err;
    }

    if (bzst_set_block_size(fp_out, block_size))
        goto err;

    ssize_t n;
    if (mode == MODE_TEXT) {
        buffer = malloc(block_size);
        if (!buffer)
            goto err;
        size_t offset = 0;
        size_t orig_block_size = block_size;

        int in_hdr = 1;
    bigger_block:
        while ((n = hread(fp_in, buffer+offset, block_size-offset)) > 0) {
            //fprintf(stderr, "Read %ld bytes to %ld..%ld\n",
            //        n, offset, offset+n);
            int64_t rc = 0; // rec counter in this buffer
            char *cp = buffer, *cp2, *last = buffer;
            size_t remaining = n + offset;
            while ((cp2 = memchr(cp, '\n', remaining))) {
                cp2++;
                remaining -= cp2-cp;
                if (in_hdr && *cp == hdr_chr && *cp2 != hdr_chr) {
                    last = cp2;
                    in_hdr = 0;
                    break;
                } else if (in_hdr) {
                    last = cp2;
                } else if (!in_hdr || *cp != hdr_chr) {
                    in_hdr = 0;
                    if (++rc % rec_mod == 0)
                        last = cp2;
                }
                cp = cp2;
            }

            if (last == buffer) {
                //fprintf(stderr, "Record does not fit in the block size\n");
                block_size *= 1.3;
                //fprintf(stderr, "Resizing block_size to %ld\n", block_size);
                char *buf_new = realloc(buffer, block_size);
                if (!buf_new)
                    goto err;
                buffer = buf_new;
                offset += n;
                goto bigger_block;
            }
            //fprintf(stderr, "bzst_write %ld of %ld\n", last-buffer, n+offset);
            //fprintf(stderr, "<<%.*s>>\n", (int)(last-buffer), buffer);
            if (bzst_write(fp_out, buffer, last-buffer, 0) < 0)
                goto err;
            if (bzst_flush(fp_out) < 0)
                goto err;
            if (orig_block_size > offset)
                block_size = orig_block_size;

            offset = n+offset - (last-buffer);
            //fprintf(stderr, "Move %ld bytes down %ld\n",
            //        offset, (long)(last-buffer));
            memmove(buffer, last, offset);
        }
        //fprintf(stderr, "n=%ld with attempted read of %ld\n", n, block_size-offset);
        if (offset) {
            // remainder if we're not a multiple of rec_mod
            if (bzst_write(fp_out, buffer, offset, 0) < 0)
                goto err;
        }
    } else {
        while ((n = hread(fp_in, buffer, BUFSZ)) > 0) {
            if (bzst_write(fp_out, buffer, n, 1) < 0)
                goto err;
        }
    }

    ret = 0;
 err:
    if (fp_in)
        if (hclose(fp_in))
            fprintf(stderr, "error closing input\n");

    if (fp_out)
        bzst_close(fp_out);

    if (pool)
        hts_tpool_destroy(pool);

    if (buffer != buff_static)
        free(buffer);

    return ret;
}


// Decodes a file, from start to end (specify 0s for all) inclusively.
static int decode(char *in, char *out, uint64_t start, uint64_t end,
                  int nthreads) {
    bzst *fp_in = NULL;
    hFILE *fp_out = NULL;
    int ret = 1;
    size_t remaining = end ? end - start + 1 : INT64_MAX;
    hts_tpool *pool = NULL;

    if (!(fp_in = bzst_open(in, "r"))) {
        perror(in);
        return -1;
    }
    if (!(fp_out = hopen(out, "w"))) {
        perror(out);
        return -1;
    }

    if (nthreads) {
        pool = hts_tpool_init(nthreads);
        if (!pool)
            goto err;
        if (bzst_thread_pool(fp_in, pool, 0) < 0)
            goto err;
    }

    if (end) {
        errno = 0;
        if (bzst_seek(fp_in, start) < 0) {
            if (errno == ERANGE) {
                fprintf(stderr, "Range is beyond end of file\n");
                goto success;
            } else {
                fprintf(stderr, "Failed to seek in bzst file\n");
                goto err;
            }
        }
    }

    ssize_t n = 0;
#if 0
    char buffer[BUFSZ];
    while (remaining > 0 && (n = bzst_read(fp_in, buffer, BUFSZ)) > 0) {
        if (hwrite(fp_out, buffer, MIN(n, remaining)) != n)
            goto err;

        remaining -= n;
    }
#else
    const char *buf0;
    while (remaining > 0 &&
           (n = bzst_read_zero_copy(fp_in, &buf0,
                                     MIN(BUFSZ,remaining))) > 0) {
        if (hwrite(fp_out, buf0, n) != n)
            goto err;

        remaining -= n;
    }
#endif

    if (n == 0 || remaining == 0)
        ret = 0;
 success:
 err:
    if (fp_in)  ret |= (bzst_close(fp_in) < 0);
    if (fp_out) ret |= (hclose(fp_out) < 0);

    if (pool)
        hts_tpool_destroy(pool);

    if (ret)
        fprintf(stderr, "Error decoding file\n");

    return ret ? -1 : 0;
}

/* ------------------------------------------------------------------------
 * ZSTD file structure listing
 */
static int list_bzst_block_metadata(hFILE *fp, uint64_t cpos,
                                     uint32_t *len_p, int level) {
    uint32_t len = *len_p;
    unsigned char buf[4];

    if (level > 1) {
        if (len < 4)
            return -1;
        if (hread(fp, buf, 4) != 4)
            return -1;
        uint32_t csize = le_to_u32(buf);
        printf("BZST block meta skippable, len %d @ %"PRId64
               ", next block csize %u\n", len, cpos, csize);

        len -= 4;
        if (len > 0) {
            char *m = malloc(len);
            if (!m)
                return -1;
            if (hread(fp, m, len) != len)
                return -1;
            printf("    Meta data: %.*s\n", len, m);
            free(m);
        }
        *len_p = 0;
    }

    return 0;
}

static int list_bzst_file_metadata(hFILE *fp, uint64_t cpos,
                                    uint32_t *len_p, int level) {
    uint32_t len = *len_p;

    if (level > 1) {
        printf("File meta skippable, len %d @ %"PRId64"\n",
               len, cpos);
        if (len > 0) {
            char *m = malloc(len);
            if (!m)
                return -1;
            if (hread(fp, m, len) != len)
                return -1;
            printf("    Meta data: %.*s\n", len, m);
            *len_p = 0;
            free(m);
        }
    }

    return 0;
}

static int list_bzst_genomic_index(hFILE *fp, uint64_t cpos,
                                    uint32_t *len_p, int level) {
    uint32_t len = *len_p;

    char *g = malloc(len);
    if (!g)
        return -1;
    if (hread(fp, g, len) != len)
        goto err;

    if (level>1)
        printf("BZST genomic index, len %d, %s @ %"PRId64"\n",
               len, g[0]&1 ? "compressed" : "uncompressed",
               cpos);

    if (level > 2) {
        uint8_t *gp = (uint8_t *)g, *g_end = gp+len;
        gp++; // flag: unused
        if (gp+4 > g_end)
            goto err;
        int nchr = le_to_u32(gp); gp += 4;
        for (int i = 0; i < nchr; i++) {
            if (gp+5 > g_end)
                goto err;
            gp++; // flag
            int index_sz = le_to_u32(gp); gp += 4;
            // Should index_sz be uint64_t?
            printf("    chr-id %d/%d, size %d\n",
                   i+1, nchr, index_sz);
            if (gp+index_sz * 20 > g_end)
                goto err;
            for (int j = 0; j < index_sz; j++) {
                int tid = le_to_u32(gp); gp += 4;
                int beg = le_to_u32(gp); gp += 4;
                int end = le_to_u32(gp); gp += 4;
                uint64_t upos = le_to_u64(gp); gp += 8;
                printf("        %4d: %d, %d..%d at %"PRId64"\n",
                       j, tid, beg, end, upos);
            }
        }
    }

    free(g);
    *len_p = 0;
    return 0;

 err:
    free(g);
    return -1;
}

static int list_bzst_index(hFILE *fp, uint32_t len, uint64_t cpos, int level) {
    uint8_t *g = NULL, buf[22];

    uint8_t flags;
    uint64_t count;
    if (level > 1) {
        if (hread(fp, (char *)buf, 17) != 17)
            return -1;
        flags = buf[0];
        count = le_to_u64(buf+1);
        uint64_t fsize = le_to_u64(buf+9);
        printf("BZST_INDEX, len %d @ %"PRId64", flags=%d, count=%"PRId64
               ", usize=%"PRId64"\n", len+1, cpos, flags, count, fsize);
        len -= 17;
    }

    if (level > 2) {
        // Dump bzst seekable index
        uint8_t *g = malloc(len);
        if (!g)
            goto err;
        if (hread(fp, g, len) != len)
            goto err;

        uint8_t *gp = g, *g_end = gp+len;

        if (count > INT64_MAX/24 || len < 24*count) {
            fprintf(stderr, "BZST index is too small. Aborting\n");
            return -1;
        }

        for (uint64_t i = 0; i < count; i++) {
            uint64_t upos  = le_to_u64(gp);
            uint64_t cpos  = le_to_u64(gp+8);
            uint64_t csize = le_to_u64(gp+16);
            gp += 24;
            printf("    %"PRIu64": upos %"PRIu64", cpos %"PRIu64
                   ", csize %"PRIu64"\n", i, upos, cpos, csize);
        }

        len -= 24*count;
        if (len != 20) {
            fprintf(stderr, "BZST index: expected 20 byte footer, got %d\n",
                    len);
            return -1;
        }
        if (hseek(fp, len, SEEK_CUR) < 0)
            return -1;

        free(g);
    } else {
        if (hseek(fp, len, SEEK_CUR) < 0)
            return -1;
    }

    return 0;

 err:
    free(g);
    return -1;
}


/*
  +--------------------+------------+
  |    Magic_Number    | 4 bytes    |
  +--------------------+------------+
  |    Frame_Header    | 2-14 bytes |
  +--------------------+------------+
  |     Data_Block     | n bytes    |
  +--------------------+------------+
  | [More Data_Blocks] |            |
  +--------------------+------------+
  | [Content_Checksum] | 0-4 bytes  |
  +--------------------+------------+
*/
static int list_zstd_data_frame(hFILE *fp, uint64_t cpos, int level,
                                int64_t *nblock) {
    uint8_t buf[8];

    if (level>1)
        printf("ZStd data frame @ %"PRId64"\n", cpos);

    //-------------------- Frame Header Descriptor
    /*
      +-------------------------+-----------+
      | Frame_Header_Descriptor | 1 byte    |
      +-------------------------+-----------+
      |   [Window_Descriptor]   | 0-1 byte  |
      +-------------------------+-----------+
      |     [Dictionary_ID]     | 0-4 bytes |
      +-------------------------+-----------+
      |  [Frame_Content_Size]   | 0-8 bytes |
      +-------------------------+-----------+
    */
    uint8_t flag;
    if (hread(fp, &flag, 1) < 0)
        goto err;
    int dict_flag = flag&3;
    int checksum_flag = (flag>>2) & 1;
    int single_flag = (flag>>5) & 1;
    int fcs_flag = flag>>6;
    if (level>2)
        printf("    Hdr descriptor: dict=%d, checksum=%d, single=%d, "
               "frame_content_sz_flag=%d\n",
               dict_flag, checksum_flag, single_flag, fcs_flag);
    int fcs_field_size = 1<<fcs_flag;
    if (fcs_field_size == 1 && !single_flag)
        fcs_field_size = 0;

    // Window Descriptor
    if (!single_flag) {
        uint8_t w;
        if (hread(fp, &w, 1) != 1)
            goto err;
        if (level>1) {
            int exponent = w>>3;
            int mantissa = w&7;
            int windowLog = 10 + exponent;
            int windowBase = 1 << windowLog;
            int windowAdd = (windowBase / 8) * mantissa;
            int Window_Size = windowBase + windowAdd;
            printf("    Window size=%d\n", Window_Size);
        }
    }

    // Dictionary ID
    if (dict_flag) {
        int did_field_size = 1<<(dict_flag-1);
        if (hread(fp, buf, did_field_size) != did_field_size)
            goto err;
        if (level>2)
            printf("    DictID field size %d\n", did_field_size);
    }

    // Frame Content Size
    if (fcs_field_size) {
        uint64_t frame_size;
        memset(buf, 0, 8);
        if (hread(fp, buf, fcs_field_size) != fcs_field_size)
            goto err;
        frame_size = le_to_u64(buf);
        if (level>2)
            printf("    Frame size %"PRIu64"\n", frame_size);
    }

    //-------------------- Data blocks
    int last_block = 0;
    do {
        (*nblock)++;
        if (hread(fp, buf, 3) != 3)
            goto err;
        buf[3] = 0;
        uint32_t blk_size = le_to_u32(buf);
        last_block = blk_size & 1;
        int block_type = (blk_size>>1) & 3;
        blk_size >>= 3;
        char *type[] = {"Raw", "RLE", "Compressed", "Reserved"};
        if (level>2)
            printf("    Block last=%d type=%s size=%d\n",
                   last_block, type[block_type], blk_size);
        if (hseek(fp, blk_size, SEEK_CUR) < 0) {
            goto err;
        }
    } while (!last_block);

    //-------------------- Content Checksums
    if (checksum_flag)
        if (hread(fp, buf, 4) != 4)
            goto err;

    return 0;

 err:
    return -1;
}

#define IS_SKIPPABLE(m) ((m & ZSTD_MAGIC_SKIPPABLE_MASK) == ZSTD_MAGIC_SKIPPABLE_START)

static int list_file(char *fn, int level) {
    hFILE *fp = hopen(fn, "r");
    if (!fp)
        return -1;
    int64_t nheader = 0;
    int64_t nblockmeta = 0;
    int64_t nfilemeta = 0;
    int64_t nblockhdr = 0;
    int64_t ndata = 0;
    int64_t nsindex = 0;
    int64_t ngindex = 0;
    int64_t nblock = 0;

    unsigned char buf[20];
    bzst_frame_t type;
    while (hread(fp, (char *)buf, 4) == 4) {
        uint64_t cpos = htell(fp)-4;
        uint32_t magic = le_to_u32(buf), len;
        switch (magic) {
        case BZST_SKIPPABLE_ID: // BZST skippable magic
            if (hread(fp, (char *)buf, 5) != 5)
                goto err;
            len = le_to_u32(buf); // frame size
            type = buf[4];
            len--;

            switch (type) {
            case BZST_HEADER: {
                char dat[100];
                uint8_t version = buf[5];
                char *format = (char *)&buf[6];
                int l = hread(fp, dat, MIN(100, len));
                if (l<0)
                    goto err;

                nheader++;
                if (level>1)
                    printf("BZST_HEADER: len %d, version %d, format \"%.4s\"\n",
                           len+1, version, format);
                len -= l;

                break;
            }

            case BZST_BLOCK_HEADER: {
                char dat[100];
                int l = hread(fp, dat, MIN(100, len));
                if (l<0)
                    goto err;

                uint64_t csize = le_to_u64(dat);
                uint64_t usize = le_to_u64(dat+8);
                uint8_t flags = dat[16];
                // + 4 byte checksum
                
                nblockhdr++;
                if (level > 1)
                    printf("BZST_BLOCK_HEADER: len %d, csize %"PRIu64
                           ", usize %"PRIu64", flags %d\n",
                           len, csize, usize, flags); 

                len -= l;
                break;
            }

            case BZST_INDEX: {
                nsindex++;
                if (list_bzst_index(fp, len, cpos, level) < 0)
                    goto err;
                break;
            }

//            case GZST_BLOCK_META: {
//                nblockmeta++;
//                if (list_bzst_block_metadata(fp, cpos, &len, level) < 0)
//                    goto err;
//                break;
//            }

//            case GZST_FILE_META: {
//                nfilemeta++;
//                if (list_bzst_file_metadata(fp, cpos, &len, level) < 0)
//                    goto err;
//                break;
//            }
//
//            case GZST_GENOMIC_INDEX: {
//                ngindex++;
//                if (list_bzst_genomic_index(fp, cpos, &len, level) < 0)
//                    goto err;
//
//                break;
//            }

            default:
                fprintf(stderr, "Unknown skippable frame with sub-type %d\n",
                        type);
                abort();
            }

            // Consume any remainder of the skippable frame
            if (len > 0)
                if (hseek(fp, len, SEEK_CUR) < 0)
                    goto err;
            break;

        case ZSTD_MAGICNUMBER: // 0xFD2FB528, zstd data frame
            ndata++;
            if (list_zstd_data_frame(fp, cpos, level, &nblock) < 0)
                goto err;
            break;

        default:
            if (IS_SKIPPABLE(magic)) {
                if (hread(fp, (char *)buf, 4) != 4)
                    goto err;
                len = le_to_u32(buf);
                printf("Unknown zstd skippable frame %x, len %d\n",
                       magic, len);
                if (hseek(fp, len, SEEK_CUR) < 0)
                    goto err;
            } else {
                printf("Unknown zstd frame of illegal magic %x\n", magic);
                goto err;
            }
        }
    }

    printf("Frames: %10"PRId64"\tBZST_HEADER\n", nheader);
    printf("Frames: %10"PRId64"\tBZST_BLOCK_HEADER\n", nblockhdr);
    printf("Frames: %10"PRId64"\tZSTD data frames\n", ndata);
    printf("Blocks: %10"PRId64"\tZSTD data blocks\n", nblock);
    printf("Frames: %10"PRId64"\tframe metadata\n", nblockmeta);  
    printf("Frames: %10"PRId64"\tfile metadata\n", nfilemeta); 
    printf("Frames: %10"PRId64"\tgenomic index\n", ngindex);
    printf("Frames: %10"PRId64"\tBZST_INDEX\n", nsindex);

    return hclose(fp);

 err:
    perror("List_file");
    hclose_abruptly(fp);
    return -1;
}

/* ------------------------------------------------------------------------
 */
static void usage(FILE *fp) {
    fprintf(fp, "Usage: bzip2 [opts] [file]\n");
    fprintf(fp, "    -h         show usage\n");
    fprintf(fp, "    -c         output to stdout\n");
    fprintf(fp, "    -o FILE    output to FILE\n");
    fprintf(fp, "    -d         decompress\n");
    fprintf(fp, "    -b SIZE    Specify block size, with optional suffix K, M or G\n");
    fprintf(fp, "    -@ INT     Use INT threads\n");
    fprintf(fp, "    -0 to -19  Specify zstd compression level\n");
    fprintf(fp, "    -c         output to stdout\n");
    fprintf(fp, "    -r RANGE   Uncompressed byte range to decode, eg -r 10M-20M\n");
    fprintf(fp, "    -l         List contents.  Use 2 or 3 times for more verbose output\n");
    fprintf(fp, "    -m MODE    Select input data type; text, bin, auto\n");
    fprintf(fp, "    -L INT     Number of lines per record in text mode\n");
    fprintf(fp, "    -H CHR     Character starting lines to identify file headers\n");
    exit(fp == stderr);
}

int main(int argc, char **argv) {
    int c;
    int level = 0;
    size_t blk_size = BZST_DEFAULT_BLOCK_SIZE;
    int compress = 1;
    char *infn = NULL;
    char *outfn = NULL;
    int64_t start = 0, end = 0;
    int nthreads = 0, list = 0;
    enum mode_t mode = MODE_BIN;
    int rec_mod = 1;
    char hdr_chr = '@';

    while ((c = getopt(argc, argv, "co:dhb:0123456789r:@:lm:L:H:")) >= 0) {
        switch(c) {
        case '@':
            nthreads = atoi(optarg);
            break;

        case 'c': // stdout
            outfn = "-";
            break;

        case 'o':
            outfn = optarg;
            break;

        case 'd':
            compress = 0;
            break;

        case 'b': {
            char *unit;
            blk_size = strtol(optarg, &unit, 0);
            if (*unit == 'k' || *unit == 'K')
                blk_size <<= 10;
            else if (*unit == 'm' || *unit == 'M')
                blk_size <<= 20;
            else if (*unit == 'g' || *unit == 'G')
                blk_size <<= 30;

            if (blk_size > BZST_MAX_BLOCK_SIZE) {
                fprintf(stderr, "Block size is too large, limit is %d bytes\n",
                        BZST_MAX_BLOCK_SIZE);
                return 1;
            }
            if (blk_size < 8)
                blk_size = 8;
            break;
        }
        case '0': case '1': case '2': case '3': case '4':
        case '5': case '6': case '7': case '8': case '9':
            // Also --long option?  Maybe with e.g. -7l ?
            level = level*10 + c-'0';
            break;

        case 'r': {
            // range from X-Y
            // X or Y may be absent which is 0 to Y or X to END
            // Also supports k, m and g suffixes.
            char *endp;
            start = strtol(optarg, &endp, 0);
            if (*endp == 'k' || *endp == 'K')
                start <<= 10, endp++;
            else if (*endp == 'm' || *endp == 'M')
                start <<= 20, endp++;
            else if (*endp == 'g' || *endp == 'G')
                start <<= 30, endp++;

            if (start < 0) {
                // 0 to Y
                end = start;
                start = 0;
            } else if (*endp == '-' && endp[1] != 0) {
                end = strtol(endp+1, &endp, 0);
                if (*endp == 'k' || *endp == 'K')
                    end <<= 10, endp++;
                else if (*endp == 'm' || *endp == 'M')
                    end <<= 20, endp++;
                else if (*endp == 'g' || *endp == 'G')
                    end <<= 30, endp++;
            } else {
                // X to EOF
                end = INT64_MAX;
            }

            if (end < start) {
                fprintf(stderr, "Illegal range '%s'\n", optarg);
                return 1;
            }
            break;
        }

        case 'l':
            list++;
            break;

        case 'm':
            mode = MODE_AUTO;
            if (strcasecmp(optarg, "text") == 0)
                mode = MODE_TEXT;
            else if (strcasecmp(optarg, "bin") == 0)
                mode = MODE_BIN;
            else if (strcasecmp(optarg, "auto") == 0)
                mode = MODE_AUTO;
            else
                fprintf(stderr, "Unknown mode '%s', using 'auto'\n", optarg);
            break;

        case 'L':
            rec_mod = atoi(optarg);
            if (rec_mod < 1)
                rec_mod = 1;
            break;

        case 'H':
            hdr_chr = *optarg;
            break;

        case 'h':
            usage(stdout);
            // fall through
        default:
            usage(stderr);
        }
    }

    if (optind == argc && isatty(fileno((FILE *)stdout)))
        usage(stdout);

    infn = (optind < argc) ? argv[optind++] : "-";
    if (!outfn)
        outfn = (optind < argc) ? argv[optind++] : "-";
    if (!level)
        level = BZST_DEFAULT_LEVEL;

    if (list)
        return list_file(infn, list) ? 1 : 0;

    int ret;
    if (compress) {
        ret = encode(infn, outfn, level, blk_size, nthreads,
                     mode, rec_mod, hdr_chr) ? 1 : 0;
    } else {
        ret = decode(infn, outfn, start, end, nthreads) ? 1 : 0;
    }
    return ret;
}
