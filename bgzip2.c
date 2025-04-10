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

#include "htslib/bgzf2.h"
#include "htslib/hfile.h"
#include "htslib/thread_pool.h"
#include "htslib/hts_endian.h"

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif

//#define BUFSZ 65536
#define BUFSZ 5000000

static int convert(char *in, char *out, int level, long block_size,
                   int nthreads) {
    hFILE *fp_in = NULL;
    bgzf2 *fp_out = NULL;
    char buffer[BUFSZ];
    int ret = -1;
    hts_tpool *pool = NULL;

    fp_in = hopen(in, "r");
    if (!fp_in) goto err;

    char omode[30];
    snprintf(omode, 30, "w%d", level);
    fp_out = bgzf2_open(out, omode);
    if (!fp_out) goto err;

    if (nthreads) {
        pool = hts_tpool_init(nthreads);
        if (!pool)
            goto err;
        if (bgzf2_thread_pool(fp_out, pool, 0) < 0)
            goto err;
    }

    if (bgzf2_set_block_size(fp_out, block_size))
        goto err;

    ssize_t n;
    while ((n = hread(fp_in, buffer, BUFSZ)) > 0) {
        if (bgzf2_write(fp_out, buffer, n, 1) < 0)
            goto err;
    }

    ret = 0;
 err:
    if (fp_in)
        if (hclose(fp_in))
            fprintf(stderr, "error closing input\n");

    if (fp_out)
        bgzf2_close(fp_out);

    if (pool)
        hts_tpool_destroy(pool);

    return ret;
}


// Decodes a file, from start to end (specify 0s for all) inclusively.
static int decode(char *in, char *out, uint64_t start, uint64_t end,
                  int nthreads) {
    bgzf2 *fp_in = NULL;
    hFILE *fp_out = NULL;
    int ret = 1;
    size_t remaining = end ? end - start + 1 : INT64_MAX;
    hts_tpool *pool = NULL;

    if (!(fp_in = bgzf2_open(in, "r"))) {
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
        if (bgzf2_thread_pool(fp_in, pool, 0) < 0)
            goto err;
    }

    if (end) {
        errno = 0;
        if (bgzf2_seek(fp_in, start) < 0) {
            if (errno == ERANGE) {
                fprintf(stderr, "Range is beyond end of file\n");
                goto success;
            } else {
                fprintf(stderr, "Failed to seek in bgzf2 file\n");
                goto err;
            }
        }
    }

    ssize_t n = 0;
#if 0
    char buffer[BUFSZ];
    while (remaining > 0 && (n = bgzf2_read(fp_in, buffer, BUFSZ)) > 0) {
        if (hwrite(fp_out, buffer, MIN(n, remaining)) != n)
            goto err;

        remaining -= n;
    }
#else
    const char *buf0;
    while (remaining > 0 &&
           (n = bgzf2_read_zero_copy(fp_in, &buf0,
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
    if (fp_in)  ret |= (bgzf2_close(fp_in) < 0);
    if (fp_out) ret |= (hclose(fp_out) < 0);

    if (pool)
        hts_tpool_destroy(pool);

    if (ret)
        fprintf(stderr, "Error decoding file\n");

    return ret ? -1 : 0;
}

int list_file(char *fn, int level) {
    hFILE *fp = hopen(fn, "r");
    if (!fp)
        return -1;
    int64_t nmagic = 0;
    int64_t npzstd = 0;
    int64_t ndata = 0;
    int64_t nsindex = 0;
    int64_t ngindex = 0;
    int64_t nblock = 0;

    unsigned char buf[8];
    while (hread(fp, (char *)buf, 4) == 4) {
        uint32_t magic = le_to_u32(buf), len;
        switch (magic) {
        case 0x184D2A5B: // BGZF magic numer
            if (hread(fp, (char *)buf, 4) != 4)
                goto err;
            len = le_to_u32(buf);
            char dat[100];
            int l = hread(fp, dat, len<100?len:100);
            if (l<0)
                goto err;

            if (nmagic == 0) {
                nmagic++;
                if (level>1)
                    printf("BGZF magic, len %d: %.*s\n", len, l, dat);
            } else {
                ngindex++;
                if (level>1)
                    printf("BGZF genomic index, len %d, %s\n", len,
                           dat[0]&1 ? "compressed" : "uncompressed");

                // TODO: decode index too?
            }
            if (len > 100)
                if (hseek(fp, len-100, SEEK_CUR) < 0)
                    goto err;
            break;

        case 0x184D2A50: // pzstd skippable
            npzstd++;
            if (hread(fp, (char *)buf, 4) != 4)
                goto err;
            len = le_to_u32(buf);
            if (level>1)
                printf("PZSTD skippable, len %d\n", len);
            if (hseek(fp, len, SEEK_CUR) < 0)
                goto err;
            break;

        case 0x184D2A5E: // Seekable index
            nsindex++;
            if (hread(fp, (char *)buf, 4) != 4)
                goto err;
            len = le_to_u32(buf);
            if (level>1)
                printf("Seekable Index, len %d\n", len);
            if (hseek(fp, len, SEEK_CUR) < 0)
                goto err;
            break;

        case ZSTD_MAGICNUMBER: // 0xFD2FB528, zstd data frame
            ndata++;
            if (level>1)
                printf("ZStd data frame\n");

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
            if (level>1)
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
                nblock++;
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

            break;

        default:
            if ((magic & ZSTD_MAGIC_SKIPPABLE_MASK)
                == ZSTD_MAGIC_SKIPPABLE_START) {
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

    printf("Frames: %10"PRId64"\tBGZF2 magic number\n", nmagic);
    printf("Frames: %10"PRId64"\tdata\n", ndata);
    printf("Frames: %10"PRId64"\tpzstd header\n", npzstd); 
    printf("Frames: %10"PRId64"\tgenomic index\n", ngindex);
    printf("Frames: %10"PRId64"\tseekable index\n", nsindex);
    printf("Blocks: %10"PRId64"\tdata blocks\n", nblock);

    return hclose(fp);

 err:
    perror("List_file");
    hclose_abruptly(fp);
    return -1;
}

static void usage(FILE *fp) {
    fprintf(fp, "Usage: bzip2 [opts] [file]\n");
    fprintf(fp, "    -h         show usage\n");
    fprintf(fp, "    -c         output to stdout\n");
    fprintf(fp, "    -d         decompress\n");
    fprintf(fp, "    -b SIZE    Specify block size, with optional suffix K, M or G\n");
    fprintf(fp, "    -@ INT     Use INT threads\n");
    fprintf(fp, "    -0 to -19  Specify zstd compression level\n");
    fprintf(fp, "    -c         output to stdout\n");
    fprintf(fp, "    -r RANGE   Uncompressed byte range to decode, eg -r 10M-20M\n");
    fprintf(fp, "    -l         List contents.  Use 2 or 3 times for more verbose output\n");
    exit(fp == stderr);
}

int main(int argc, char **argv) {
    int c;
    int level = 0;
    long blk_size = BGZF2_DEFAULT_BLOCK_SIZE;
    int compress = 1;
    char *infn = NULL;
    char *outfn = NULL;
    int64_t start = 0, end = 0;
    int nthreads = 0, list = 0;

    while ((c = getopt(argc, argv, "cdhb:0123456789r:@:l")) >= 0) {
        switch(c) {
        case '@':
            nthreads = atoi(optarg);
            break;

        case 'c': // stdout
            outfn = "-";
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

            if (blk_size > BGZF2_MAX_BLOCK_SIZE) {
                fprintf(stderr, "Block size is too large, limit is %d bytes\n",
                        BGZF2_MAX_BLOCK_SIZE);
                return 1;
            }
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
    outfn = (optind < argc) ? argv[optind++] : "-";
    if (!level)
        level = BGZF2_DEFAULT_LEVEL;

    if (list)
        return list_file(infn, list) ? 1 : 0;

    int ret;
    if (compress) {
        ret = convert(infn, outfn, level, blk_size, nthreads) ? 1 : 0;
    } else {
        ret = decode(infn, outfn, start, end, nthreads) ? 1 : 0;
    }
    return ret;
}
