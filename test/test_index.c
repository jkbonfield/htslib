/* test_index.c - for testing index generation


  Copyright (C) 2014 Genome Research Ltd.

  Author: Rob Davies <rmd+git@sanger.ac.uk>

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "htslib/sam.h"

#define USAGE "Usage: %s [-bc] [-m INT] <in.bam>\n"

int main(int argc, char **argv) {
    char type = 'b';
    int  min_shift = 14;
    int  opt;
    int  res;
    
    while ((opt = getopt(argc, argv, "bcm:")) >= 0) {
        switch (opt) {
        case 'b':
            type = 'b';
            break;
        case 'c':
            type = 'c';
            break;
        case 'm':
            type = 'c';
            min_shift = atoi(optarg);
            if (min_shift <= 0) {
                fprintf(stderr, "-m option should be > 0\n");
                return EXIT_FAILURE;
            }
            break;
        default:
            fprintf(stderr, USAGE, argv[0]);
            return EXIT_FAILURE;
        }
    }

    if (optind == argc) {
        fprintf(stderr, USAGE, argv[0]);
        return EXIT_FAILURE;
    }

    /* Make the index.  Passing min_shift <= 0 makes a .bai index,
                                           > 0 makes a .csi index */
    res = bam_index_build(argv[optind], type == 'c' ? min_shift : 0);

    return 0 == res ? EXIT_SUCCESS : EXIT_FAILURE;
}
