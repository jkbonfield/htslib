/*
 * BAM threading example.
 *
 * We load a block of BAM records, run an arbitrary function on that block of
 * records to modify them somehow, and then write out the block of BAM records.
 *
 * This therefore requires the correct order and a way to read/write
 * asynchronously.  For that we have separate reader and writer threads
 * to marshall the data, in addition to main, and a queue for result data.
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
#define NBAM 4096

// Commands/options for the reader and writer threads to communicate with
enum sam_cmd {
    SAM_NONE = 0,
    SAM_READ_END,
    SAM_READ_DONE,
};

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

    // A thread to read a block of records along with mutex and condition var
    pthread_t        read_t;
    pthread_mutex_t  read_m;
    pthread_cond_t   read_c;
    enum sam_cmd     read_opt;

    // Similarly to write a block of records.
    pthread_t        write_t;
    pthread_mutex_t  write_m;
    pthread_cond_t   write_c;
    enum sam_cmd     write_opt;

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

    if (!j)
	// Equiv to bam1_init on all elements of j->ba[]
	j = calloc(1, sizeof(*j));

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
    for (i = 0; i < j->nbam; i++) {
	bam1_t *b = &j->ba[i];
	uint8_t *seq = bam_get_seq(b);
        int counts1[16] = {0}, counts2[16] = {0};

#if 0
	// Slow to demo the threading gain in main
	for (int k = 0; k < b->core.l_qseq; k++)
	    counts1[bam_seqi(seq, k)]++;
#else
	// Faster implementation
	int k;
	for (k = 0; k < (b->core.l_qseq&~1); k+=2) {
	    uint8_t s = seq[k/2];
	    counts1[s&0xf]++;
	    counts2[s>>4]++;
	}
	if (k < b->core.l_qseq)
	    counts1[bam_seqi(seq, k)]++;
#endif

        // Compute and add back GC; A=1, C=2, G=4, T=8
        int AT = counts1[1]+counts2[1]+counts1[8]+counts2[8];
        int CG = counts1[2]+counts2[2]+counts1[4]+counts2[4];
        //float gc = (100.0 * CG) / (CG+AT);
        //bam_aux_append(b, "cg", 'f', sizeof(float), (uint8_t *)&gc);

        // For compatibility with qtask_ordered.
        float gc = (double)CG / b->core.l_qseq;
        bam_aux_append(b, "xr", 'f', sizeof(float), (uint8_t *)&gc);

	for (int k = 0; k < 16; k++)
	    counts[k] += counts1[k]+counts2[k];
    }

    pthread_mutex_lock(&j->s->bamjob_m);
    j->s->nr += j->nbam;
    for (i = 0; i < 16; i++)
	j->s->counts[i] += counts[i];
    pthread_mutex_unlock(&j->s->bamjob_m);

    //bam_job_free(j->s, j);
    //return NULL;
    return j;
}

// This thread runs continuously for the duration of the program.
// Its purpose is to read in blocks of data and add to a queue
static void *read_thread(void *vp) {
    state *s = (state *)vp;

    //fprintf(stderr, "read_thread started\n");

    int eof = 0;

    // Loop while data exists
    for (;;) {
	// Check command
	pthread_mutex_lock(&s->read_m);
	int opt = s->read_opt;
	pthread_mutex_unlock(&s->read_m);

	//fprintf(stderr, "opt=%d eof=%d\n", opt, eof);
	if (opt == SAM_READ_END || eof)
	    break;

	// Load a block of data
	bam_job *j = bam_job_alloc(s);
	j->s = s;

	int i, r;
	for (i = 0; i < NBAM; i++) {
	    // bam1_init
	    if ((r = sam_read1(s->in, s->h, &j->ba[i])) < 0)
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
	if (hts_tpool_dispatch(s->p, s->q, s->func, j) < 0)
	    // We could also check read_opt here to see if we've been
	    // explicitly signalled a close, to distinguish errors from
	    // normal shutdown.
	    return NULL;
    }

    // We've finished, so end.  We signal to main thread we're done too.
    // We also submit one EOF job so the writer knows we're complete.
    bam_job *j = bam_job_alloc(s);
    j->s = s;
    j->nbam = 0;
    j->status = 1;
    if (hts_tpool_dispatch(s->p, s->q, s->func, j) < 0)
        return NULL;

    pthread_mutex_lock(&s->read_m);
    s->read_opt = SAM_READ_DONE;
    pthread_cond_signal(&s->read_c);
    pthread_mutex_unlock(&s->read_m);

    //fprintf(stderr, "read_thread exited\n");

    return NULL;
}

static void *write_thread(void *vp) {
    state *s = (state *)vp;

    hts_tpool_result *r;
    int err = 0;
    while ((r = hts_tpool_next_result_wait(s->q))) {
	//fprintf(stderr, "Got result\n");
	bam_job *bj = (bam_job *)hts_tpool_result_data(r);
        if (bj->nbam == 0 && bj->status == 1) {
            // EOF
            bam_job_free(bj->s, bj);
            hts_tpool_delete_result(r, 0);
            break;
        }

        for (int i = 0; i < bj->nbam; i++)
            err |= sam_write1(s->out, s->h, &bj->ba[i]);
        bam_job_free(bj->s, bj);
	hts_tpool_delete_result(r, 0);
	usleep(10);
    }

    return err ? NULL : (void *)1;
}

int sam_loop(htsFile *in, htsFile *out, hts_tpool *p, int nthreads) {
    sam_hdr_t *h = NULL;

    // Load and copy header
    h = sam_hdr_read(in);
    if (h == NULL)
        return -1;

    if (sam_hdr_write(out, h) < 0)
	return -1;

    // Create the sam state and add thread pool and queues
    state s;
    memset(&s, 0, sizeof(s));
    s.p = p;
    s.q = hts_tpool_process_init(s.p, nthreads*2, 0);
    s.in = in;
    s.out = out;
    s.h = h;

    // bam_func is our actual function that operates on blocks of BAM records.
    s.func = bam_func;

    // Create dedicated read and write threads.  These are low CPU and mainly
    // just grouping blocks of data together for CPU operations.  The main
    // decompress and compress is already executed in the shared pool.
    pthread_mutex_init(&s.read_m, NULL);
    pthread_cond_init(&s.read_c,  NULL);
    s.read_opt = 0;

    pthread_mutex_init(&s.write_m, NULL);
    pthread_cond_init(&s.write_c,  NULL);
    s.write_opt = 0;

    pthread_mutex_init(&s.bamjob_m, NULL);
    s.free_jobs = NULL;

    pthread_create(&s.read_t, NULL, read_thread, &s);
    pthread_create(&s.write_t, NULL, write_thread, &s);

    // Now wait on the reader to exit.
    // Main read-write loop
    for (;;) {
	pthread_mutex_lock(&s.read_m);
	pthread_cond_wait(&s.read_c, &s.read_m);
	enum sam_cmd opt = s.read_opt;
	pthread_mutex_unlock(&s.read_m);

	if (opt == SAM_READ_DONE)
	    break;
    }

    // Wait for queue to drain as we may still have jobs process
    hts_tpool_process_flush(s.q);
    pthread_join(s.read_t, NULL);
    pthread_join(s.write_t, NULL);
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

    return 0;
}

int main(int argc, char *argv[])
{
    htsFile *in, *out;
    int nthreads = 0, c;

    // Parse arguments
    while ((c = getopt(argc, argv, "@:")) >= 0) {
        switch (c) {
        case '@': nthreads = atoi(optarg); break;
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
    if (sam_loop(in, out, p.pool, nthreads) < 0)
	exit(1);

    // Close files and free memory
    if (hts_close(in) != 0 || hts_close(out) != 0)
	exit(1);

    if (p.pool)
        hts_tpool_destroy(p.pool);

    return 0;
}
