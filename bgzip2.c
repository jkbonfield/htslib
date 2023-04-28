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

    fp_in = hopen(in, "r");
    if (!fp_in) goto err;

    char omode[30];
    snprintf(omode, 30, "w%d", level);
    fp_out = bgzf2_open(out, omode);
    if (!fp_out) goto err;

    if (nthreads) {
        hts_tpool *pool = hts_tpool_init(nthreads);
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
    
    if (hclose(fp_in) < 0 || bgzf2_close(fp_out) < 0)
	return -1;

    return 0;

 err:
    if (fp_in)
        if (hclose(fp_in))
            fprintf(stderr, "error closing input\n");

    if (fp_out)
        bgzf2_close(fp_out);

    return -1;
}


// TODO: specify a region
static int decode(char *in, char *out, uint64_t start, uint64_t end,
                  int nthreads) {
    bgzf2 *fp_in = NULL;
    hFILE *fp_out = NULL;
    int ret = 1;
    size_t remaining = end ? end - start : INT64_MAX;

    if (!(fp_in = bgzf2_open(in, "r"))) {
        perror(in);
        return -1;
    }
    if (!(fp_out = hopen(out, "w"))) {
        perror(out);
        return -1;
    }

    if (nthreads) {
        hts_tpool *pool = hts_tpool_init(nthreads);
        if (!pool)
            goto err;
        if (bgzf2_thread_pool(fp_in, pool, 0) < 0)
            goto err;
    }

    if (end) {
        int err = load_seekable_index(fp_in);
        if (err < -2)
            fprintf(stderr, "BGZF2 seekable-index not found\n");

        if (err < 0)
            goto err;

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

    if (ret)
        fprintf(stderr, "Error decoding file\n");

    return ret ? -1 : 0;
}

static void usage(FILE *fp) {
    fprintf(fp, "Usage: bzip2 [opts] [file]\n");
    exit(fp == stderr);
}

int main(int argc, char **argv) {
    int c;
    int level = 0;
    long blk_size = BGZF2_DEFAULT_BLOCK_SIZE;
    int compress = 1;
    char *infn = NULL;
    char *outfn = NULL;
    uint64_t start = 0, end = 0;
    int nthreads = 0;

    while ((c = getopt(argc, argv, "cdhb:0123456789r:@:")) >= 0) {
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

        case 'h':
            usage(stdout);
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

    if (compress) {
        return convert(infn, outfn, level, blk_size, nthreads) ? 1 : 0;
    } else {
        return decode(infn, outfn, start, end, nthreads) ? 1 : 0;
    }
}
