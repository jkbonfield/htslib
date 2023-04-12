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

#include "htslib/bgzf2.h"
#include "htslib/hfile.h"

#define BUFSZ 65536

static int convert(char *in, char *out, int level, long block_size) {
    hFILE *fp_in = NULL;
    bgzf2 *fp_out = NULL;
    char buffer[BUFSZ];

    fp_in = hopen(in, "r");
    if (!fp_in) goto err;

    char omode[30];
    snprintf(omode, 30, "w%d", level);
    fp_out = bgzf2_open(out, omode);
    if (!fp_out) goto err;

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
static int decode(char *in, char *out) {
    char buffer[BUFSZ];
    bgzf2 *fp_in = NULL;
    hFILE *fp_out = NULL;
    int ret = 1;

    fp_in = bgzf2_open(in, "r");
    fp_out = hopen(out, "w");

    size_t n;
    while ((n = bgzf2_read(fp_in, buffer, BUFSZ)) > 0) {
        if (hwrite(fp_out, buffer, n) != n)
            goto err;
    }

    if (n == 0)
        ret = 0;

 err:
    if (fp_in)  ret |= (bgzf2_close(fp_in) < 0);
    if (fp_out) ret |= (hclose(fp_out) < 0);

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

    while ((c = getopt(argc, argv, "cdhb:0123456789")) >= 0) {
        switch(c) {
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
            break;
        }
        case '0': case '1': case '2': case '3': case '4':
        case '5': case '6': case '7': case '8': case '9': 
            // Also --long option?  Maybe with e.g. -7l ?
            level = level*10 + c-'0';
            break;

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
        return convert(infn, outfn, level, blk_size);
    } else {
        return decode(infn, outfn);
    }
}
