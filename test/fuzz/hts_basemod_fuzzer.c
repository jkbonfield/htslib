/*  test/fuzz/hts_basemod_fuzzer.c -- Fuzz driver for the base modification API.

    Reads SAM/BAM/CRAM from the fuzz input and exercises the MM/ML parsing
    and iteration functions.  The aim is to find memory errors (out of
    bounds accesses, use after free, leaks, uninitialised data) triggered
    by arbitrary input, not to validate what the functions return.  Return
    values are therefore only checked for initialisation, never for
    plausibility.

    1. Per record, sharing a single hts_base_mod_state across all records:
       - bam_parse_basemod / bam_parse_basemod2
       - bam_mods_recorded, bam_mods_query_type, bam_mods_queryi
       - bam_mods_at_next_pos over every base
       - bam_next_basemod
       - bam_mods_at_qpos with a variable stride

    2. Via mpileup, with per-read state created by a constructor and freed
       by a destructor (as samtools mpileup does):
       - bam_mods_at_qpos at each pileup column
       - bam_plp_insertion and bam_plp_insertion_mod
       - optional overlap detection and a low maxcnt, so that reads get
         modified or discarded while they hold modification state
*/

#include <config.h>

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../htslib/hfile.h"
#include "../../htslib/hts.h"
#include "../../htslib/sam.h"
#include "../../htslib/kstring.h"

// Under MSan, report any returned data that was never initialised,
// regardless of whether we subsequently branch on it.
#if defined(__has_feature)
#  if __has_feature(memory_sanitizer)
#    include <sanitizer/msan_interface.h>
#    define CHECK_INIT_MEM(p, len) __msan_check_mem_is_initialized((p), (len))
#  endif
#endif
#ifndef CHECK_INIT_MEM
#  define CHECK_INIT_MEM(p, len) ((void)0)
#endif
#define CHECK_INIT(x) CHECK_INIT_MEM(&(x), sizeof(x))

#define MAX_MODS      256
#define MAX_REF_SPAN  100000   // skip reads spanning more than this in pileup
#define MAX_COLUMNS   500000   // cap on pileup columns, to avoid timeouts

// Bits of 'param', derived from a hash of the input
#define P_UNCHECKED     (1<<0) // HTS_MOD_REPORT_UNCHECKED in per-record pass
#define P_NMODS_SHIFT   1      // bits 1-2: size of mods[] array passed in
#define P_STRIDE_SHIFT  3      // bits 3-5: qpos stride-1 for bam_mods_at_qpos
#define P_OVERLAPS      (1<<6) // bam_mplp_init_overlaps
#define P_PLP_UNCHECKED (1<<7) // HTS_MOD_REPORT_UNCHECKED in pileup pass
#define P_LOW_MAXCNT    (1<<8) // tiny maxcnt, so pileup drops reads

static const int nmods_sizes[4] = {1, 2, 5, MAX_MODS};

// Only the first n_mods entries are filled out, even if more were found.
static void check_mods(const hts_base_mod *mods, int n, int n_mods) {
    if (n > n_mods)
        n = n_mods;
    if (n > 0)
        CHECK_INIT_MEM(mods, n * sizeof(*mods));
}

static int parse_mods(const bam1_t *b, hts_base_mod_state *m, uint32_t flags) {
    return flags
        ? bam_parse_basemod2(b, m, flags)
        : bam_parse_basemod(b, m);
}

static void query_types(hts_base_mod_state *m) {
    int ntype = 0, strand, implicit;
    char canonical;
    int *types = bam_mods_recorded(m, &ntype);
    CHECK_INIT(ntype);

    for (int i = 0; i < ntype; i++) {
        CHECK_INIT(types[i]);
        if (bam_mods_query_type(m, types[i], &strand, &implicit,
                                &canonical) == 0) {
            CHECK_INIT(strand);
            CHECK_INIT(implicit);
            CHECK_INIT(canonical);
        }
        if (bam_mods_queryi(m, i, &strand, &implicit, &canonical) == 0) {
            CHECK_INIT(strand);
            CHECK_INIT(implicit);
            CHECK_INIT(canonical);
        }
    }

    // Out of range and unknown queries should be handled safely
    (void) bam_mods_queryi(m, ntype, &strand, &implicit, &canonical);
    (void) bam_mods_query_type(m, '#', &strand, &implicit, &canonical);
}

/* ----------------------------------------------------------------------
 * Pass 1: direct per-record API use
 */
static void fuzz_record(const bam1_t *b, hts_base_mod_state *m,
                        uint32_t flags, int n_mods, int stride) {
    hts_base_mod mods[MAX_MODS];
    int len = b->core.l_qseq, i, n;

    if (parse_mods(b, m, flags) < 0)
        return;
    query_types(m);

    // Every base in turn
    for (i = 0; i < len; i++) {
        if ((n = bam_mods_at_next_pos(b, m, mods, n_mods)) < 0)
            break;
        check_mods(mods, n, n_mods);
    }

    // Only the modified bases.  The iteration count is capped, in case
    // a bug leaves it failing to advance.
    if (parse_mods(b, m, flags) < 0)
        return;
    int pos = -1;
    for (i = 0; i <= len; i++) {
        if ((n = bam_next_basemod(b, m, mods, n_mods, &pos)) <= 0)
            break;
        CHECK_INIT(pos);
        check_mods(mods, n, n_mods);
    }

    // Random access with increasing qpos
    if (parse_mods(b, m, flags) < 0)
        return;
    for (i = 0; i < len; i += stride) {
        if ((n = bam_mods_at_qpos(b, i, m, mods, n_mods)) < 0)
            break;
        check_mods(mods, n, n_mods);
    }
}

static void per_record_pass(htsFile *in, sam_hdr_t *hdr, uint32_t param) {
    hts_base_mod_state *m = hts_base_mod_state_alloc();
    bam1_t *b = bam_init1();
    if (!m || !b)
        goto out;

    uint32_t flags = (param & P_UNCHECKED) ? HTS_MOD_REPORT_UNCHECKED : 0;
    int n_mods = nmods_sizes[(param >> P_NMODS_SHIFT) & 3];
    int stride = 1 + ((param >> P_STRIDE_SHIFT) & 7);

    // The state is deliberately reused between records, including after
    // partial iterations and failed parses.
    while (sam_read1(in, hdr, b) >= 0)
        fuzz_record(b, m, flags, n_mods, stride);

 out:
    bam_destroy1(b);
    if (m)
        hts_base_mod_state_free(m);
}

/* ----------------------------------------------------------------------
 * Pass 2: via mpileup
 */
typedef struct {
    htsFile *fp;
    sam_hdr_t *hdr;
    uint32_t flags;
} plp_data;

static int plp_read(void *data, bam1_t *b) {
    plp_data *d = data;
    int ret;
    while ((ret = sam_read1(d->fp, d->hdr, b)) >= 0) {
        if (b->core.tid < 0 || (b->core.flag & BAM_FUNMAP))
            continue;
        // Huge D/N operations generate a column per base; skip them
        if (bam_endpos(b) - b->core.pos > MAX_REF_SPAN)
            continue;
        break;
    }
    return ret;
}

static int plp_construct(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    plp_data *d = data;
    hts_base_mod_state *m = hts_base_mod_state_alloc();
    if (m && bam_parse_basemod2(b, m, d->flags) < 0) {
        hts_base_mod_state_free(m);
        m = NULL;
    }
    cd->p = m;
    return 0;
}

static int plp_destruct(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    if (cd->p)
        hts_base_mod_state_free(cd->p);
    cd->p = NULL;
    return 0;
}

static void pileup_pass(htsFile *in, sam_hdr_t *hdr, uint32_t param) {
    plp_data d = {
        in, hdr, (param & P_PLP_UNCHECKED) ? HTS_MOD_REPORT_UNCHECKED : 0
    };
    void *dp = &d;
    kstring_t ins = KS_INITIALIZE;
    hts_base_mod mods[MAX_MODS];

    bam_mplp_t iter = bam_mplp_init(1, plp_read, &dp);
    if (!iter)
        return;
    bam_mplp_constructor(iter, plp_construct);
    bam_mplp_destructor(iter, plp_destruct);
    bam_mplp_set_maxcnt(iter, (param & P_LOW_MAXCNT) ? 3 : 1000);
    if ((param & P_OVERLAPS) && bam_mplp_init_overlaps(iter) < 0)
        goto out;

    int tid, n_plp, columns = 0;
    hts_pos_t pos;
    const bam_pileup1_t *plp[1];

    while (bam_mplp64_auto(iter, &tid, &pos, &n_plp, plp) > 0) {
        if (++columns > MAX_COLUMNS)
            break; // remaining reads go through the destructor on destroy

        for (int i = 0; i < n_plp; i++) {
            const bam_pileup1_t *p = plp[0] + i;
            hts_base_mod_state *m = p->cd.p;
            int r, del_len;

            if (m && !p->is_del && !p->is_refskip) {
                int n = bam_mods_at_qpos(p->b, p->qpos, m, mods, MAX_MODS);
                check_mods(mods, n, MAX_MODS);
            }

            // Called on every column; non-insertion sites give blank
            // strings.
            if ((r = bam_plp_insertion(p, &ins, &del_len)) >= 0) {
                CHECK_INIT(del_len);
                CHECK_INIT_MEM(ins.s, ins.l);
            }
            if ((r = bam_plp_insertion_mod(p, NULL, &ins, &del_len)) >= 0) {
                CHECK_INIT(del_len);
                CHECK_INIT_MEM(ins.s, ins.l);
            }

            // With state; this also advances it through the inserted bases
            if (m && p->indel > 0) {
                if ((r = bam_plp_insertion_mod(p, m, &ins, &del_len)) >= 0) {
                    CHECK_INIT(del_len);
                    CHECK_INIT_MEM(ins.s, ins.l);
                }
            }
        }
    }

 out:
    ks_free(&ins);
    bam_mplp_destroy(iter);
}

/* ----------------------------------------------------------------------
 * Driver
 */
static htsFile *open_mem(const uint8_t *data, size_t size) {
    uint8_t *copy = malloc(size);
    if (copy == NULL)
        abort();
    memcpy(copy, data, size);

    hFILE *memfile = hopen("mem:", "rb:", copy, size);
    if (memfile == NULL) {
        free(copy);
        return NULL;
    }

    htsFile *in = hts_hopen(memfile, "data", "rb");
    if (in == NULL) {
        if (hclose(memfile) != 0)
            abort();
        return NULL;
    }

    if (hts_get_format(in)->category != sequence_data) {
        hts_close(in);
        return NULL;
    }
    return in;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size == 0)
        return 0;

    hts_set_log_level(HTS_LOG_OFF);

    // Options come from a hash of the input, so each input is deterministic
    uint32_t param = 2166136261u;
    for (size_t i = 0; i < size; i++)
        param = (param ^ data[i]) * 16777619u;

    htsFile *in;
    sam_hdr_t *hdr;

    if (!(in = open_mem(data, size)))
        return 0;
    if ((hdr = sam_hdr_read(in)) != NULL) {
        per_record_pass(in, hdr, param);
        sam_hdr_destroy(hdr);
    }
    hts_close(in);

    // Second pass from a fresh copy, as the first consumed the stream
    if (!(in = open_mem(data, size)))
        return 0;
    if ((hdr = sam_hdr_read(in)) != NULL) {
        pileup_pass(in, hdr, param);
        sam_hdr_destroy(hdr);
    }
    hts_close(in);

    return 0;
}
