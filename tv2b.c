/*
 * Trivial example of BAM threading in file format read/write.
 * This is the same task as tv2.c, but we don't have threading for the
 * ACGT count loop.
 *
 * Given we're encoding with bgzf, that is the bottleneck so it makes little
 * difference until we get to a high thread count.
 *
 * Test with 1 million NovaSeq reads.
 * Threads  8     16    32    64
 * tv2:     8.0s  4.2s  2.5s  3.0s  (NBAM=4096)
 * tv2b:    8.1s  6.0s  6.4s  6.4s
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

int sam_loop(htsFile *in, htsFile *out) {
    sam_hdr_t *h = NULL;

    h = sam_hdr_read(in);
    if (h == NULL)
        return -1;
    if (sam_hdr_write(out, h) < 0)
	return -1;

    int r, err = 0;
    int64_t counts[16] = {0}, nr = 0;
    bam1_t *b = bam_init1();

    while ((r = sam_read1(in, h, b)) >= 0) {
	int counts1[16] = {0}, counts2[16] = {0};
	nr++;
	int j;
	uint8_t *seq = bam_get_seq(b);
#if 0
	// Slow implementation is useful for demonstration of main thread
	// becoming bottlenecked.
	// On 1 million novaseq recs with 8 threads,
	// this was 2.9s real 5.8s cpu, vs below 1.8s real 4.6s cpu.
	for (j = 0; j < b->core.l_qseq; j++)
	    counts1[bam_seqi(seq, j)]++;
#else
	// Faster implementation.  We do two bases at a time, and accumulate
	// to two buffers so we don't have successive read incr write clashes.
	for (j = 0; j < (b->core.l_qseq&~1); j+=2) {
	    uint8_t s = seq[j/2];
	    counts1[s&0xf]++;
	    counts2[s>>4]++;
	}
	if (j < b->core.l_qseq)
	    counts1[bam_seqi(seq, j)]++;
#endif

        // Compute and add back GC; A=1, C=2, G=4, T=8
        int AT = counts1[1]+counts2[1]+counts1[8]+counts2[8];
        int CG = counts1[2]+counts2[2]+counts1[4]+counts2[4];
        float gc = (100.0 * CG) / (CG+AT);
        bam_aux_append(b, "cg", 'f', sizeof(float), (uint8_t *)&gc);
	err |= sam_write1(out, h, b) < 0;

	for (j = 0; j < 16; j++)
	    counts[j] += counts1[j]+counts2[j];
    }
    bam_destroy1(b);

    int i;
    printf("%"PRId64" reads\n", nr);
    for (i = 0; i < 16; i++)
	printf("%c %"PRId64"\n", seq_nt16_str[i], counts[i]);

    sam_hdr_destroy(h);

    return -err;
}

int main(int argc, char *argv[])
{
    htsFile *in, *out;
    int nthreads = 0, c;

    // Parse arguments
    while ((c = getopt(argc, argv, "@:")) >= 0) {
        switch (c) {
        case '@':
	    nthreads = atoi(optarg);
	    break;
        }
    }
    if (argc != optind + 2) {
        fprintf(stderr, "Usage: tv [-@ threads] in.bam out.bam\n");
        return 1;
    }

    // Open files
    in = hts_open(argv[optind], "r");
    char out_mode[6] = {'w', 0};
    sam_open_mode(out_mode+1, argv[optind+1], NULL);
    out = hts_open(argv[optind+1], out_mode);

    if (!in) {
        fprintf(stderr, "Error opening \"%s\"\n", argv[optind]);
        return 1;
    }
    if (!out) {
        fprintf(stderr, "Error opening \"%s\"\n", argv[optind+1]);
        return 1;
    }

    // Create and share the thread pool
    htsThreadPool p = {NULL, 0};
    if (nthreads > 0) {
        p.pool = hts_tpool_init(nthreads);
        if (!p.pool) {
            fprintf(stderr, "Error creating thread pool\n");
	    return 1;
	}

	hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
	hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);
    }

    // Call main loop
    if (sam_loop(in, out) < 0)
	exit(1);

    // Close files and free memory
    if (hts_close(in) != 0 || hts_close(out) != 0)
	exit(1);

    if (p.pool)
        hts_tpool_destroy(p.pool);

    return 0;
}
