/*
 * BAM threading example.
 *
 * We load a block of BAM records, run an arbitrary function on that block of
 * records to process them.  Order doesn't matter.
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

// An array of BAM records
#define NBAM 2048
//#define NBAM 512

struct state;

// A worker job
typedef struct bam_job {
    bam1_t ba[NBAM];
    int nbam;
    int status;
    struct state *s;
    struct bam_job *next;
} bam_job;

// State for the threading
typedef struct state {
    // Thread pool details
    hts_tpool *p;
    hts_tpool_process *q;
    htsFile *in, *out;
    sam_hdr_t *h;

    // A linked list of free blocks of bam jobs
    pthread_mutex_t  bamjob_m;
    bam_job *free_jobs;

    // The function that we'll be passing blocks of BAM records to
    void *(*func)(void *arg);

    // Somewhere for us to store the results
    uint64_t nr;
    uint64_t counts[16];
} state;

static bam_job *bam_job_alloc(state *s) {
    bam_job *j = NULL;
    pthread_mutex_lock(&s->bamjob_m);
    if (s->free_jobs) {
	j = s->free_jobs;
	s->free_jobs = s->free_jobs->next;
    }
    pthread_mutex_unlock(&s->bamjob_m);

    if (!j) {
	// Equiv to bam1_init on all elements of j->ba[]
	j = calloc(1, sizeof(*j));
//	int i;
//	for (i = 0; i < NBAM; i++)
//	    // remove BAM_USER_OWNS_DATA
//	    bam_set_mempolicy(&j->ba[i], BAM_USER_OWNS_STRUCT);
    }

    return j;
}

static void bam_job_free(state *s, bam_job *j) {
    pthread_mutex_lock(&s->bamjob_m);
    j->next = s->free_jobs;
    s->free_jobs = j;
    pthread_mutex_unlock(&s->bamjob_m);
}

static void bam_job_destroy(bam_job *j) {
    int i;
    for (i = 0; i < NBAM; i++)
	if (bam_get_mempolicy(&j->ba[i]) & BAM_USER_OWNS_DATA)
	    free(j->ba[i].data);
    free(j);
}

// An arbitrary function to operate on a block of BAM records.
// It's given a bam_job arg
static void *bam_func(void *vp) {
    bam_job *j = (bam_job *)vp;
    j->status = 1; // Mark as OK

    //fprintf(stderr, "Execute bam_func on %d records => %d\n", j->nbam, x);

    // Compute A C G T tallies
    int counts[16] = {0}, i;
#if 1
    // Simplest.  Slow, but best demonstrates issue of CPU in main thread
    for (i = 0; i < j->nbam; i++) {
	bam1_t *b = &j->ba[i];
	uint8_t *seq = bam_get_seq(b);

	for (int k = 0; k < b->core.l_qseq; k++)
	    counts[bam_seqi(seq, k)]++;
    }

#elif 0
    // Faster implementation. Add 2 bases at a time to counts
    int counts2[16] = {0};
    for (i = 0; i < j->nbam; i++) {
	bam1_t *b = &j->ba[i];
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
    for (i = 0; i < j->nbam; i++) {
	bam1_t *b = &j->ba[i];
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

    pthread_mutex_lock(&j->s->bamjob_m);
    j->s->nr += j->nbam;
    for (i = 0; i < 16; i++) {
	j->s->counts[i] += counts[i];
    }
    pthread_mutex_unlock(&j->s->bamjob_m);

    bam_job_free(j->s, j);

    return NULL;
}


int sam_loop(htsFile *in, hts_tpool *p, int nthreads,
	     void *(*func)(void *arg)) {
    sam_hdr_t *h = NULL;
    int ret = 0;

    // Load and copy header
    h = sam_hdr_read(in);
    if (h == NULL)
        return -1;

    // Create the sam state and add thread pool and queues
    state s;
    memset(&s, 0, sizeof(s));
    s.p = p;
    s.q = hts_tpool_process_init(s.p, nthreads*2, 1);
    s.in = in;
    s.h = h;

    pthread_mutex_init(&s.bamjob_m, NULL);
    s.free_jobs = NULL;

    // Main loop 
    for (;;) {
	bam_job *j = bam_job_alloc(&s);
	j->s = &s;

	int i, r, eof = 0;
	for (i = 0; i < NBAM; i++) {

// Attempt to prefetch for writes, but detrimental, esp sam.gz
#if 0
#define N 1
	    if (i+N < NBAM) {
		uint8_t *ptr = (uint8_t *)&j->ba[i+N];
		__builtin_prefetch(ptr,    1);
		__builtin_prefetch(ptr+64, 1);
		if (N>0) {
		    bam1_t *b = &j->ba[i+N-1];
		    if (b->m_data >= 64)
			__builtin_prefetch(b->data, 1);
		    if (b->m_data >= 128)
			__builtin_prefetch(b->data+64, 1);
		    //if (b->m_data >= 192)
		    //	__builtin_prefetch(b->data+128, 1);
		    //if (b->m_data >= 256)
		    //	__builtin_prefetch(b->data+192, 1);
		}
	    }
#elif 0
	    __builtin_prefetch(&j->ba[i], 1);
	    __builtin_prefetch(64+(uint8_t *)&j->ba[i], 1);
#endif

	    if ((r = sam_read1(s.in, s.h, &j->ba[i])) < 0)
		break;
	}
	if (r < -1) {
	    ret = -1;
	    break;
	}

	if (i != NBAM)
	    // eof after this final block
	    eof = 1;

	if (i == 0) {
	    // empty block is eof right now
	    free(j); // NB leaks bam object. TODO: create bam_job_free()
	    break;
	}
	j->nbam = i;
	j->status = 0;

	// We now have nb bam objects in ba, so create and dispatch a job.
	// This can block, but queue shutdown interrupts it and returns 0.
	if (hts_tpool_dispatch(s.p, s.q, func, j) < 0) {
	    ret = -1;
	    break;
	}

	if (eof)
	    break;
    }

    //pthread_join(s.read_t, NULL);

    // Wait for queue to drain as we may still have jobs process
    hts_tpool_process_flush(s.q);
    hts_tpool_process_destroy(s.q);

    int i;
    printf("%"PRId64" reads\n", s.nr);
    for (i = 0; i < 16; i++)
	printf("%c %"PRId64"\n", seq_nt16_str[i], s.counts[i]);

    sam_hdr_destroy(h);

    bam_job *j, *next = NULL;
    for (j = s.free_jobs; j; j = next) {
	next = j->next;
	bam_job_destroy(j);
    }

    return ret;
}

int main(int argc, char *argv[])
{
    htsFile *in;
    int nthreads = 0, c;

    // Parse arguments
    while ((c = getopt(argc, argv, "@:")) >= 0) {
        switch (c) {
        case '@': nthreads = atoi(optarg); break;
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

    // Create and share the thread pool
    htsThreadPool p = {NULL, 0};
    if (nthreads > 0) {
        p.pool = hts_tpool_init(nthreads);
        if (!p.pool) {
            fprintf(stderr, "Error creating thread pool\n");
	    return 1;
	}

	hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
    }

    // Call main loop, which executes bam_func() in parallel
    if (sam_loop(in, p.pool, nthreads, bam_func) < 0)
	exit(1);

    // Close files and free memory
    if (hts_close(in) != 0)
	exit(1);

    if (p.pool)
        hts_tpool_destroy(p.pool);

    return 0;
}
