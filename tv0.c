/*
 * Trivial example of BAM threading in file format decode only.
 * The main thread does a read & ACGT count loop.
 */

#include <config.h>

#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"

int sam_loop(htsFile *in) {
    sam_hdr_t *h = NULL;

    h = sam_hdr_read(in);
    if (h == NULL)
        return -1;

    int r, i;
    int64_t counts[16] = {0}, nr = 0;
    bam1_t *b = bam_init1();

#if 1
    // Simplest.  Slow, but best demonstrates issue of CPU in main thread
    while ((r = sam_read1(in, h, b)) >= 0) {
	nr++;
	uint8_t *seq = bam_get_seq(b);

	for (int k = 0; k < b->core.l_qseq; k++)
	    counts[bam_seqi(seq, k)]++;
    }

#elif 0
    // Faster implementation. Add 2 bases at a time to counts
    int counts2[16] = {0};
    while ((r = sam_read1(in, h, b)) >= 0) {
	nr++;
	uint8_t *seq = bam_get_seq(b);

	int k;
	for (k = 0; k < (b->core.l_qseq&~1); k+=2) {
	    uint8_t s = seq[k/2];
	    counts [s&0xf]++;
	    counts2[s>>4]++;
	}
	if (k < b->core.l_qseq)
	    counts[bam_seqi(seq, k)]++;
    }

    for (i = 0; i < 16; i++)
	counts[i] += counts2[i];

#else
    // Fastest.  Add to a 256 wide array and separate out at end
    int64_t counts256[256] = {0};
    while ((r = sam_read1(in, h, b)) >= 0) {
	nr++;
	uint8_t *seq = bam_get_seq(b);

	int k;
	for (k = 0; k < b->core.l_qseq/2; k++)
	    counts256[seq[k]]++;
	k*=2;
	if (k < b->core.l_qseq)
	    counts[bam_seqi(seq, k)]++;
    }

    for (i = 0; i < 256; i++) {
	counts[i&0xf]+=counts256[i];
	counts[i>>4 ]+=counts256[i];
    }
#endif

    bam_destroy1(b);

    printf("%"PRId64" reads\n", nr);
    for (i = 0; i < 16; i++)
	printf("%c %"PRId64"\n", seq_nt16_str[i], counts[i]);

    sam_hdr_destroy(h);

    return 0;
}

int main(int argc, char *argv[])
{
    htsFile *in;
    int nthreads = 0, c;

    // Parse arguments
    while ((c = getopt(argc, argv, "@:")) >= 0) {
        switch (c) {
        case '@':
	    nthreads = atoi(optarg);
	    break;
        }
    }
    if (argc != optind + 1) {
        fprintf(stderr, "Usage: tv [-@ threads] in.bam\n");
        return 1;
    }

    // Open files
    in = hts_open(argv[optind], "r");
    if (!in) {
        fprintf(stderr, "Error opening \"%s\"\n", argv[optind]);
        return 1;
    }

    if (nthreads > 0)
	hts_set_opt(in,  HTS_OPT_NTHREADS, nthreads);

    // Call main loop
    if (sam_loop(in) < 0)
	exit(1);

    // Close files and free memory
    if (hts_close(in) != 0)
	exit(1);

    return 0;
}
