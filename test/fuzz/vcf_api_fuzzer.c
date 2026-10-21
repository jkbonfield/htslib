/*  test/fuzz/vcf_api_fuzzer.c -- Fuzz driver for the VCF/BCF record API.

    Reading and writing files is covered by hts_open_fuzzer.  This instead
    exercises the header and record manipulation API on whatever the input
    parses to.  The aim is to find memory errors (out of bounds accesses,
    use after free, double frees, leaks, uninitialised data) triggered by
    arbitrary input, not to validate results: invalid input may legitimately
    produce invalid output.  Return values are therefore only checked for
    initialisation (under MSan), never for plausibility.

    The fuzzer takes care to use the API as documented, so that a failure
    should indicate a problem with htslib's handling of the input rather
    than a misuse by this driver.  Notably:
    - strings passed to update functions are never pointers into the
      record being updated;
    - after any update function fails, that record copy is discarded;
    - a fixed header is used for bcf_translate, and translation and
      subsetting stop if the reading header grows (VCF parsing adds
      undefined tags), as the derived headers no longer match it.

    Header API:
      bcf_hdr_set_samples, bcf_hdr_format, bcf_hdr_get_version,
      bcf_hdr_seqnames, bcf_hrec_format, bcf_hrec_dup, bcf_hrec_find_key,
      bcf_hdr_parse_line, bcf_hdr_get_hrec, bcf_hdr_dup, bcf_hdr_remove,
      bcf_hdr_append, bcf_hdr_sync, bcf_hdr_set_version, bcf_hdr_add_sample,
      bcf_hdr_init, bcf_hdr_merge, bcf_hdr_subset, bcf_hdr_id2* macros

    Record API:
      bcf_unpack, bcf_dup, vcf_format, vcf_parse, bcf_get_info_values and
      bcf_get_format_values (all types), bcf_get_format_string,
      bcf_get_info, bcf_get_fmt, bcf_get_info_id, bcf_get_fmt_id,
      bcf_update_info, bcf_update_format, bcf_update_format_string,
      bcf_update_genotypes, bcf_update_alleles(_str), bcf_update_id,
      bcf_add_id, bcf_update_filter, bcf_add_filter, bcf_remove_filter,
      bcf_has_filter, variant type functions, bcf_format_gt, bcf_translate,
      bcf_subset, and from vcfutils: bcf_calc_ac, bcf_gt_type,
      bcf_trim_alleles, bcf_remove_allele_set
*/

#include <config.h>

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../htslib/hfile.h"
#include "../../htslib/hts.h"
#include "../../htslib/kbitset.h"
#include "../../htslib/kstring.h"
#include "../../htslib/vcf.h"
#include "../../htslib/vcfutils.h"

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
#define CHECK_STR(s) do {                               \
        const char *s_ = (s);                           \
        if (s_) CHECK_INIT_MEM(s_, strlen(s_));         \
    } while (0)

#define MAX_RECORDS  2000  // records processed per input
#define MAX_HDR_ADD  8     // header-defined tags added to each record

// Low bits of 'param'; other uses go via pick()
#define P_SAMPLES_MASK 3      // bcf_hdr_set_samples mode
#define P_PASS         (1<<2) // 'pass' argument to bcf_remove_filter

typedef struct {
    bcf_hdr_t *hdr;      // header records are read with
    bcf_hdr_t *merged;   // extra lines + hdr, for bcf_translate
    bcf_hdr_t *hsub;     // subset of hdr's samples, for bcf_subset
    int *imap, nsub;
    int n_id, n_ctg;     // hdr dictionary sizes when merged/hsub were made
    int sample_subset;   // bcf_hdr_set_samples was applied
    uint32_t param;

    // Buffers reused across calls, as callers normally do.  Each is only
    // ever used with 'hdr', as bcf_get_format_string sizes fstr by the
    // header's sample count.
    int32_t *i32;  int ni32;
    int64_t *i64;  int ni64;
    float   *f;    int nf;
    char    *str;  int nstr;
    char   **fstr; int nfstr;
    int     *ac;   int nac;
    kstring_t ks, ks2;
} fuzz_ctx;

static const int32_t i32_vals[] = {
    0, 1, -1, 127, -120, 128, 32767, -32760, 32768, INT32_MAX, INT32_MIN + 8
};
#define N_I32_VALS (sizeof(i32_vals) / sizeof(*i32_vals))

static const float f_vals[] = { 0.0f, 1.5f, -1e30f, 1e30f };
#define N_F_VALS (sizeof(f_vals) / sizeof(*f_vals))

// Deterministic choice in [0,n) from the input hash and a call-site salt
static unsigned pick(const fuzz_ctx *c, unsigned salt, unsigned n) {
    uint32_t x = c->param ^ (salt * 0x9e3779b9u);
    x ^= x >> 16; x *= 0x45d9f3bu; x ^= x >> 16;
    return n ? x % n : 0;
}

static int grow(void **p, int *m, size_t n, size_t sz) {
    if (n <= (size_t) *m)
        return 0;
    if (n > INT_MAX / 2 || n > SIZE_MAX / sz)
        return -1;
    void *t = realloc(*p, n * sz);
    if (!t)
        return -1;
    *p = t;
    *m = (int) n;
    return 0;
}

static int fill_i32(fuzz_ctx *c, int n, unsigned salt) {
    if (grow((void **) &c->i32, &c->ni32, n, sizeof(*c->i32)) < 0)
        return -1;
    for (int i = 0; i < n; i++)
        c->i32[i] = i32_vals[pick(c, salt + i, N_I32_VALS)];
    if (n > 1) c->i32[n-1] = bcf_int32_vector_end;
    if (n > 2) c->i32[1]   = bcf_int32_missing;
    return 0;
}

static int fill_f(fuzz_ctx *c, int n, unsigned salt) {
    if (grow((void **) &c->f, &c->nf, n, sizeof(*c->f)) < 0)
        return -1;
    for (int i = 0; i < n; i++)
        c->f[i] = f_vals[pick(c, salt + i, N_F_VALS)];
    if (n > 1) bcf_float_set_vector_end(c->f[n-1]);
    if (n > 2) bcf_float_set_missing(c->f[1]);
    return 0;
}

static int header_grew(const fuzz_ctx *c) {
    return c->hdr->n[BCF_DT_ID] != c->n_id
        || c->hdr->n[BCF_DT_CTG] != c->n_ctg;
}

static int format_rec(fuzz_ctx *c, const bcf_hdr_t *h, const bcf1_t *v) {
    c->ks.l = 0;
    int r = vcf_format(h, v, &c->ks);
    if (r == 0)
        CHECK_INIT_MEM(c->ks.s, c->ks.l);
    return r;
}

/* ----------------------------------------------------------------------
 * Header
 */
static int set_samples(fuzz_ctx *c) {
    bcf_hdr_t *h = c->hdr;
    kstring_t *ks = &c->ks;
    int n = bcf_hdr_nsamples(h), i, ret;

    ks->l = 0;
    switch (c->param & P_SAMPLES_MASK) {
    case 0: // all samples
        return 0;

    case 1: // exclude the first
        if (n < 1)
            return 0;
        if (kputc('^', ks) < 0 || kputs(h->samples[0], ks) < 0)
            return -1;
        ret = bcf_hdr_set_samples(h, ks->s, 0);
        break;

    case 2: // every other one
        if (n < 1)
            return 0;
        for (i = 0; i < n; i += 2)
            if ((i && kputc(',', ks) < 0) || kputs(h->samples[i], ks) < 0)
                return -1;
        ret = bcf_hdr_set_samples(h, ks->s, 0);
        break;

    default: // none
        ret = bcf_hdr_set_samples(h, NULL, 0);
        break;
    }

    // ret > 0 means some names weren't found, which is fine
    if (ret < 0)
        return -1;
    c->sample_subset = 1;
    return 0;
}

static void fuzz_header(fuzz_ctx *c) {
    bcf_hdr_t *hdr = c->hdr, *dup;
    kstring_t *ks = &c->ks;
    int i, n;

    for (i = 0; i < 2; i++) {
        ks->l = 0;
        if (bcf_hdr_format(hdr, i, ks) == 0)
            CHECK_INIT_MEM(ks->s, ks->l);
    }

    CHECK_STR(bcf_hdr_get_version(hdr));

    const char **seqs = bcf_hdr_seqnames(hdr, &n);
    if (seqs) {
        for (i = 0; i < n; i++) {
            CHECK_STR(seqs[i]);
            (void) bcf_hdr_name2id(hdr, seqs[i]);
        }
        free(seqs);
    }

    for (i = 0; i < hdr->nhrec; i++) {
        bcf_hrec_t *h = hdr->hrec[i], *h2;

        ks->l = 0;
        if (bcf_hrec_format(h, ks) == 0) {
            CHECK_INIT_MEM(ks->s, ks->l);
            int len = 0;
            if ((h2 = bcf_hdr_parse_line(hdr, ks->s, &len)) != NULL)
                bcf_hrec_destroy(h2);
        }

        if ((h2 = bcf_hrec_dup(h)) != NULL)
            bcf_hrec_destroy(h2);

        int k = bcf_hrec_find_key(h, "ID");
        if (k >= 0 && h->type != BCF_HL_GEN)
            (void) bcf_hdr_get_hrec(hdr, h->type, "ID", h->vals[k],
                                    h->type == BCF_HL_STR ? h->key : NULL);
    }

    // Modify a copy: remove a line, re-add it, change version and samples.
    // This copy is never used with records.
    if (!(dup = bcf_hdr_dup(hdr)))
        return;

    if (dup->nhrec > 0) {
        bcf_hrec_t *h = dup->hrec[pick(c, 1, dup->nhrec)];
        int type = h->type, ok, k;
        char *key = NULL;

        // Take copies first, as removal frees h
        ks->l = 0;
        ok = bcf_hrec_format(h, ks) == 0;
        k = bcf_hrec_find_key(h, "ID");
        if (type == BCF_HL_GEN || type == BCF_HL_STR)
            key = h->key ? strdup(h->key) : NULL;
        else if (k >= 0 && h->vals[k])
            key = strdup(h->vals[k]);

        if (key) {
            bcf_hdr_remove(dup, type, key);
            free(key);
            if (bcf_hdr_sync(dup) == 0 && ok && bcf_hdr_append(dup, ks->s) == 0)
                (void) bcf_hdr_sync(dup);
        }
    }

    (void) bcf_hdr_set_version(dup, "VCFv4.4");
    if (bcf_hdr_add_sample(dup, "fuzz_extra_sample") == 0)
        (void) bcf_hdr_sync(dup);

    ks->l = 0;
    if (bcf_hdr_format(dup, 0, ks) == 0)
        CHECK_INIT_MEM(ks->s, ks->l);

    bcf_hdr_destroy(dup);
}

// A header with extra lines first, so that merged IDs differ from the
// originals and bcf_translate has real work to do.
static bcf_hdr_t *make_merged(fuzz_ctx *c) {
    static const char *extra[] = {
        "##FILTER=<ID=FuzzFlt,Description=\"Fuzz filter\">",
        "##INFO=<ID=FuzzInfo,Number=1,Type=Integer,Description=\"Fuzz\">",
        "##FORMAT=<ID=FuzzFmt,Number=1,Type=Float,Description=\"Fuzz\">",
        "##contig=<ID=fuzz_contig,length=1000>",
    };
    int i;
    bcf_hdr_t *m = bcf_hdr_init("w");
    if (!m)
        return NULL;

    for (i = 0; i < (int) (sizeof(extra) / sizeof(*extra)); i++)
        if (bcf_hdr_append(m, extra[i]) < 0)
            goto fail;
    if (bcf_hdr_merge(m, c->hdr) == NULL)
        goto fail;
    for (i = 0; i < bcf_hdr_nsamples(c->hdr); i++)
        if (bcf_hdr_add_sample(m, c->hdr->samples[i]) < 0)
            goto fail;
    if (bcf_hdr_sync(m) < 0)
        goto fail;
    return m;

 fail:
    bcf_hdr_destroy(m);
    return NULL;
}

// Every other sample, starting from the second (possibly none)
static void make_subset(fuzz_ctx *c) {
    int n = bcf_hdr_nsamples(c->hdr), i, k = 0;
    char **names = malloc((n + 1) * sizeof(*names));
    c->imap = malloc((n + 1) * sizeof(*c->imap));
    if (!names || !c->imap) {
        free(names);
        return;
    }

    for (i = 1; i < n; i += 2)
        names[k++] = c->hdr->samples[i];
    c->hsub = bcf_hdr_subset(c->hdr, k, names, c->imap);
    c->nsub = k;
    free(names);
}

/* ----------------------------------------------------------------------
 * Record queries (read only)
 */
static void variant_types(bcf1_t *v) {
    if (bcf_unpack(v, BCF_UN_STR) < 0)
        return;
    (void) bcf_get_variant_types(v);
    (void) bcf_is_snp(v);
    (void) bcf_has_variant_types(v, VCF_SNP | VCF_INDEL, bcf_match_overlap);
    (void) bcf_has_variant_types(v, VCF_MNP, bcf_match_exact);
    (void) bcf_has_variant_types(v, VCF_SNP | VCF_OTHER, bcf_match_subset);
    // Allele indices must be in range; htslib exits otherwise
    for (int i = 0; i < v->n_allele; i++) {
        (void) bcf_get_variant_type(v, i);
        (void) bcf_has_variant_type(v, i, VCF_SNP | VCF_REF);
    }
}

static void query_info(fuzz_ctx *c, bcf1_t *v, const char *tag) {
    const bcf_hdr_t *h = c->hdr;
    int r;

    if ((r = bcf_get_info_int32(h, v, tag, &c->i32, &c->ni32)) > 0)
        CHECK_INIT_MEM(c->i32, r * sizeof(*c->i32));
    if ((r = bcf_get_info_values(h, v, tag, (void **) &c->i64, &c->ni64,
                                 BCF_HT_LONG)) > 0)
        CHECK_INIT_MEM(c->i64, r * sizeof(*c->i64));
    if ((r = bcf_get_info_float(h, v, tag, &c->f, &c->nf)) > 0)
        CHECK_INIT_MEM(c->f, r * sizeof(*c->f));
    if ((r = bcf_get_info_string(h, v, tag, &c->str, &c->nstr)) > 0)
        CHECK_INIT_MEM(c->str, r);
    (void) bcf_get_info_flag(h, v, tag, &c->i32, &c->ni32);
    (void) bcf_get_info(h, v, tag);
}

static void query_fmt(fuzz_ctx *c, bcf1_t *v, const char *tag) {
    const bcf_hdr_t *h = c->hdr;
    int r, s;

    if ((r = bcf_get_format_int32(h, v, tag, &c->i32, &c->ni32)) > 0)
        CHECK_INIT_MEM(c->i32, r * sizeof(*c->i32));
    if ((r = bcf_get_format_float(h, v, tag, &c->f, &c->nf)) > 0)
        CHECK_INIT_MEM(c->f, r * sizeof(*c->f));
    if ((r = bcf_get_format_char(h, v, tag, &c->str, &c->nstr)) > 0)
        CHECK_INIT_MEM(c->str, r);
    if ((r = bcf_get_format_string(h, v, tag, &c->fstr, &c->nfstr)) > 0)
        for (s = 0; s < bcf_hdr_nsamples(h); s++)
            CHECK_STR(c->fstr[s]);
    (void) bcf_get_fmt(h, v, tag);
}

static void query_gt(fuzz_ctx *c, const bcf_hdr_t *h, bcf1_t *v) {
    int r, s, ial, jal;

    if ((r = bcf_get_genotypes(h, v, &c->i32, &c->ni32)) > 0)
        CHECK_INIT_MEM(c->i32, r * sizeof(*c->i32));

    bcf_fmt_t *fmt = bcf_get_fmt(h, v, "GT");
    if (!fmt)
        return;
    for (s = 0; s < v->n_sample; s++) {
        (void) bcf_gt_type(fmt, s, &ial, &jal);
        c->ks.l = 0;
        if (bcf_format_gt(fmt, s, &c->ks) == 0)
            CHECK_INIT_MEM(c->ks.s, c->ks.l);
    }
}

static void calc_ac(fuzz_ctx *c, const bcf_hdr_t *h, bcf1_t *v) {
    if (v->n_allele == 0)
        return;
    if (grow((void **) &c->ac, &c->nac, v->n_allele, sizeof(*c->ac)) < 0)
        return;
    if (bcf_calc_ac(h, v, c->ac, BCF_UN_INFO | BCF_UN_FMT) > 0)
        CHECK_INIT_MEM(c->ac, v->n_allele * sizeof(*c->ac));
}

/* ----------------------------------------------------------------------
 * Record mutations, each on a private copy
 */

// Rewrite existing INFO fields, sometimes with tweaked values
static void mutate_info(fuzz_ctx *c, bcf1_t *v) {
    const bcf_hdr_t *h = c->hdr;
    bcf1_t *w;
    int i;

    if (v->n_info == 0 || !(w = bcf_dup(v)))
        return;

    for (i = 0; i < v->n_info; i++) {
        int key = v->d.info[i].key, r, u = 0;
        const char *tag = bcf_hdr_int2id(h, BCF_DT_ID, key);
        unsigned p = pick(c, 100 + i, 64);
        if (!tag)
            continue;

        switch (p % 5) {
        case 0:
            if ((r = bcf_get_info_int32(h, w, tag, &c->i32, &c->ni32)) > 0) {
                c->i32[pick(c, 200 + i, r)] = i32_vals[p % N_I32_VALS];
                if (p & 8)
                    c->i32[r-1] = bcf_int32_vector_end;
                u = bcf_update_info_int32(h, w, tag, c->i32, r);
            }
            break;

        case 1:
            // Only a single BCF_HT_LONG value may be stored; htslib aborts
            // on more.
            r = bcf_get_info_values(h, w, tag, (void **) &c->i64, &c->ni64,
                                    BCF_HT_LONG);
            if (r == 1) {
                c->i64[0] = (p & 8) ? ((int64_t) 1 << 40)
                                    : -((int64_t) 1 << 40);
                u = bcf_update_info(h, w, tag, c->i64, 1, BCF_HT_LONG);
            }
            break;

        case 2:
            if ((r = bcf_get_info_float(h, w, tag, &c->f, &c->nf)) > 0) {
                int x = pick(c, 300 + i, r);
                if (p & 8)
                    bcf_float_set_vector_end(c->f[x]);
                else
                    bcf_float_set_missing(c->f[x]);
                u = bcf_update_info_float(h, w, tag, c->f, r);
            }
            break;

        case 3:
            if ((r = bcf_get_info_string(h, w, tag, &c->str, &c->nstr)) > 0)
                u = bcf_update_info_string(h, w, tag, (p & 8) ? "" : c->str);
            break;

        default:
            if (bcf_hdr_idinfo_exists(h, BCF_HL_INFO, key) &&
                bcf_hdr_id2type(h, BCF_HL_INFO, key) == BCF_HT_FLAG)
                u = bcf_update_info_flag(h, w, tag, NULL, (p & 8) ? 1 : 0);
            else
                u = bcf_update_info(h, w, tag, NULL, 0, BCF_HT_INT); // delete
            break;
        }

        if (u < 0)
            goto out;
    }

    format_rec(c, h, w);
    for (i = 0; i < v->n_info; i++) {
        const char *tag = bcf_hdr_int2id(h, BCF_DT_ID, v->d.info[i].key);
        if (tag)
            query_info(c, w, tag); // includes deleted tags
    }
    calc_ac(c, h, w);

 out:
    bcf_destroy(w);
}

// Rewrite existing FORMAT fields, including GT
static void mutate_fmt(fuzz_ctx *c, bcf1_t *v) {
    const bcf_hdr_t *h = c->hdr;
    int nsmpl = bcf_hdr_nsamples(h), i;
    bcf1_t *w;

    if (v->n_fmt == 0 || nsmpl == 0 || !(w = bcf_dup(v)))
        return;

    for (i = 0; i < v->n_fmt; i++) {
        const char *tag = bcf_hdr_int2id(h, BCF_DT_ID, v->d.fmt[i].id);
        unsigned p = pick(c, 400 + i, 64);
        int r, u = 0;
        if (!tag)
            continue;

        if (strcmp(tag, "GT") == 0 && (p & 1)) {
            if ((r = bcf_get_genotypes(h, w, &c->i32, &c->ni32)) > 0) {
                int x = pick(c, 500 + i, r);
                if (c->i32[x] != bcf_int32_missing &&
                    c->i32[x] != bcf_int32_vector_end)
                    c->i32[x] ^= 1; // toggle phasing
                // An out of range allele, as input can also supply
                c->i32[pick(c, 600 + i, r)] = bcf_gt_unphased(w->n_allele + 1);
                u = bcf_update_genotypes(h, w, c->i32, r);
            }
        } else {
            switch (p % 5) {
            case 0:
                if ((r = bcf_get_format_int32(h, w, tag, &c->i32,
                                              &c->ni32)) > 0) {
                    c->i32[pick(c, 700 + i, r)] = i32_vals[p % N_I32_VALS];
                    if (p & 8)
                        c->i32[r-1] = bcf_int32_vector_end;
                    u = bcf_update_format_int32(h, w, tag, c->i32, r);
                }
                break;

            case 1:
                if ((r = bcf_get_format_float(h, w, tag, &c->f, &c->nf)) > 0) {
                    bcf_float_set_vector_end(c->f[pick(c, 800 + i, r)]);
                    u = bcf_update_format_float(h, w, tag, c->f, r);
                }
                break;

            case 2:
                if ((r = bcf_get_format_char(h, w, tag, &c->str,
                                             &c->nstr)) > 0)
                    u = bcf_update_format_char(h, w, tag, c->str, r);
                break;

            case 3:
                if ((r = bcf_get_format_string(h, w, tag, &c->fstr,
                                               &c->nfstr)) > 0)
                    u = bcf_update_format_string(h, w, tag,
                                                 (const char **) c->fstr,
                                                 nsmpl);
                break;

            default:
                u = bcf_update_format(h, w, tag, NULL, 0, BCF_HT_INT);
                break;
            }
        }

        if (u < 0)
            goto out;
    }

    format_rec(c, h, w);
    for (i = 0; i < v->n_fmt; i++) {
        const char *tag = bcf_hdr_int2id(h, BCF_DT_ID, v->d.fmt[i].id);
        if (tag)
            query_fmt(c, w, tag);
    }
    query_gt(c, h, w);
    calc_ac(c, h, w);

 out:
    bcf_destroy(w);
}

// Add tags defined in the header, typed as the header says
static void add_hdr_tags(fuzz_ctx *c, bcf1_t *v) {
    static const char *svals[] = { "a", "", "fuzz,str" };
    const bcf_hdr_t *h = c->hdr;
    int nid = h->n[BCF_DT_ID], nsmpl = bcf_hdr_nsamples(h), j, s;
    const char **sv = NULL;
    bcf1_t *w;

    if (nid <= 0 || nsmpl > INT_MAX / 4 || !(w = bcf_dup(v)))
        return;

    int start = pick(c, 7, nid);
    for (j = 0; j < MAX_HDR_ADD && j < nid; j++) {
        int id = (start + j) % nid, k = 1 + pick(c, 900 + j, 3), u = 0;
        const char *tag = bcf_hdr_int2id(h, BCF_DT_ID, id);
        if (!tag)
            continue; // gap left by removed header lines

        if (bcf_hdr_idinfo_exists(h, BCF_HL_INFO, id)) {
            switch (bcf_hdr_id2type(h, BCF_HL_INFO, id)) {
            case BCF_HT_INT:
                if (fill_i32(c, k, 1000 + j) < 0)
                    goto out;
                u = bcf_update_info_int32(h, w, tag, c->i32, k);
                break;
            case BCF_HT_REAL:
                if (fill_f(c, k, 1100 + j) < 0)
                    goto out;
                u = bcf_update_info_float(h, w, tag, c->f, k);
                break;
            case BCF_HT_STR:
                u = bcf_update_info_string(h, w, tag, svals[k-1]);
                break;
            default:
                u = bcf_update_info_flag(h, w, tag, NULL, 1);
                break;
            }
            if (u < 0)
                goto out;
        }

        if (nsmpl > 0 && bcf_hdr_idinfo_exists(h, BCF_HL_FMT, id)) {
            if (strcmp(tag, "GT") == 0) {
                if (grow((void **) &c->i32, &c->ni32, nsmpl * 2,
                         sizeof(*c->i32)) < 0)
                    goto out;
                for (s = 0; s < nsmpl; s++) {
                    // Allele indices may exceed n_allele
                    c->i32[2*s]   = bcf_gt_unphased(s % 3);
                    c->i32[2*s+1] = (s & 1) ? bcf_int32_vector_end
                                            : bcf_gt_phased(s % 5);
                }
                u = bcf_update_genotypes(h, w, c->i32, nsmpl * 2);
            } else {
                switch (bcf_hdr_id2type(h, BCF_HL_FMT, id)) {
                case BCF_HT_INT:
                    if (fill_i32(c, nsmpl * k, 1200 + j) < 0)
                        goto out;
                    u = bcf_update_format_int32(h, w, tag, c->i32, nsmpl * k);
                    break;
                case BCF_HT_REAL:
                    if (fill_f(c, nsmpl * k, 1300 + j) < 0)
                        goto out;
                    u = bcf_update_format_float(h, w, tag, c->f, nsmpl * k);
                    break;
                case BCF_HT_STR:
                    if (!sv && !(sv = malloc(nsmpl * sizeof(*sv))))
                        goto out;
                    for (s = 0; s < nsmpl; s++)
                        sv[s] = svals[(s + j) % 3];
                    u = bcf_update_format_string(h, w, tag, sv, nsmpl);
                    break;
                default:
                    break; // Flag isn't valid in FORMAT
                }
            }
            if (u < 0)
                goto out;
        }
    }

    format_rec(c, h, w);
    query_gt(c, h, w);
    calc_ac(c, h, w);

 out:
    free(sv);
    bcf_destroy(w);
}

static void mutate_alleles(fuzz_ctx *c, bcf1_t *v) {
    const bcf_hdr_t *h = c->hdr;
    int na = v->n_allele, nals, i;
    char **copies = NULL;
    const char **als = NULL;
    kbitset_t *rm = NULL;
    kstring_t *ks = &c->ks2;
    bcf1_t *w;

    if (na == 0 || !(w = bcf_dup(v)))
        return;

    // String form, built from copies of the alleles
    ks->l = 0;
    for (i = 0; i < na; i++)
        if ((i && kputc(',', ks) < 0) || kputs(v->d.allele[i], ks) < 0)
            goto out;
    if (bcf_update_alleles_str(h, w, ks->s) < 0)
        goto out;
    variant_types(w);

    // Array form: reversed, plus a symbolic allele
    copies = calloc(na, sizeof(*copies));
    als = calloc(na + 1, sizeof(*als));
    if (!copies || !als)
        goto out;
    for (i = 0; i < na; i++) {
        if (!(copies[i] = strdup(v->d.allele[i])))
            goto out;
        als[na - 1 - i] = copies[i];
    }
    nals = na;
    if (na < 0xffff) // n_allele is 16 bits
        als[nals++] = "<*>";
    if (bcf_update_alleles(h, w, als, nals) < 0)
        goto out;
    variant_types(w);
    calc_ac(c, h, w);

    if (bcf_trim_alleles(h, w) < 0)
        goto out;

    if (w->n_allele > 1) {
        if (!(rm = kbs_init(w->n_allele)))
            goto out;
        for (i = 1; i < w->n_allele; i++)
            if (pick(c, 1400 + i, 2))
                kbs_insert(rm, i);
        if (bcf_remove_allele_set(h, w, rm) < 0)
            goto out;
    }
    variant_types(w);
    format_rec(c, h, w);

 out:
    if (copies) {
        for (i = 0; i < na; i++)
            free(copies[i]);
        free(copies);
    }
    free(als);
    if (rm)
        kbs_destroy(rm);
    bcf_destroy(w);
}

static void mutate_filters(fuzz_ctx *c, bcf1_t *v) {
    const bcf_hdr_t *h = c->hdr;
    char pass[] = "PASS", dot[] = ".";
    int n = v->d.n_flt, i;
    int *ids = malloc((n + 1) * sizeof(*ids));
    bcf1_t *w = NULL;

    if (!ids || !(w = bcf_dup(v)))
        goto out;

    for (i = 0; i < n; i++)
        ids[n - 1 - i] = v->d.flt[i];

    (void) bcf_has_filter(h, w, pass);
    (void) bcf_has_filter(h, w, dot);
    for (i = 0; i < n; i++) {
        const char *name = bcf_hdr_int2id(h, BCF_DT_ID, ids[i]);
        if (name)
            (void) bcf_has_filter(h, w, (char *) name);
    }

    if (bcf_update_filter(h, w, ids, n) < 0)
        goto out;
    for (i = 0; i < n; i++)
        if (bcf_remove_filter(h, w, ids[i], !!(c->param & P_PASS)) < 0)
            goto out;
    for (i = 0; i < n; i++)
        if (bcf_add_filter(h, w, ids[i]) < 0)
            goto out;
    format_rec(c, h, w);

    if (bcf_update_filter(h, w, NULL, 0) < 0)
        goto out;
    format_rec(c, h, w);

 out:
    free(ids);
    if (w)
        bcf_destroy(w);
}

static void mutate_id(fuzz_ctx *c, bcf1_t *v) {
    const bcf_hdr_t *h = c->hdr;
    char *id;
    bcf1_t *w;

    if (!v->d.id || !(id = strdup(v->d.id)))
        return;
    if (!(w = bcf_dup(v))) {
        free(id);
        return;
    }

    if (bcf_update_id(h, w, id) < 0
        || bcf_add_id(h, w, "fuzz1;fuzz2") < 0
        || bcf_add_id(h, w, id) < 0)
        goto out;
    format_rec(c, h, w);

    if (bcf_update_id(h, w, ".") < 0)
        goto out;
    format_rec(c, h, w);

 out:
    free(id);
    bcf_destroy(w);
}

// Format to text and parse it back.  Skipped with bcf_hdr_set_samples,
// as vcf_parse would then expect the original sample columns.
static void round_trip(fuzz_ctx *c, bcf1_t *v) {
    bcf1_t *w;

    if (c->sample_subset || format_rec(c, c->hdr, v) < 0)
        return;
    if (c->ks.l && c->ks.s[c->ks.l - 1] == '\n')
        c->ks.s[--c->ks.l] = '\0';

    if (!(w = bcf_init()))
        return;
    if (vcf_parse(&c->ks, c->hdr, w) == 0)
        format_rec(c, c->hdr, w);
    bcf_destroy(w);
}

static void translate(fuzz_ctx *c, bcf1_t *v) {
    bcf1_t *w;

    // bcf_translate exits on records with unchecked errors, and caches a
    // translation table sized by the header at first use.
    if (!c->merged || v->errcode || header_grew(c) || !(w = bcf_dup(v)))
        return;
    if (bcf_translate(c->merged, c->hdr, w) == 0)
        format_rec(c, c->merged, w);
    bcf_destroy(w);
}

static void subset(fuzz_ctx *c, bcf1_t *v) {
    bcf1_t *w;

    if (!c->hsub || header_grew(c) || !(w = bcf_dup(v)))
        return;
    if (bcf_subset(c->hsub, w, c->nsub, c->imap) == 0)
        format_rec(c, c->hsub, w);
    bcf_destroy(w);
}

/* ----------------------------------------------------------------------
 * Per record driver
 */
static void fuzz_record(fuzz_ctx *c, bcf1_t *v) {
    const bcf_hdr_t *h = c->hdr;
    int i;

    if (bcf_unpack(v, BCF_UN_ALL) < 0)
        return;

    CHECK_STR(bcf_seqname_safe(h, v));
    CHECK_STR(v->d.id);
    for (i = 0; i < v->n_allele; i++)
        CHECK_STR(v->d.allele[i]);
    for (i = 0; i < v->d.n_flt; i++)
        CHECK_STR(bcf_hdr_int2id(h, BCF_DT_ID, v->d.flt[i]));
    variant_types(v);

    for (i = 0; i < v->n_info; i++) {
        int key = v->d.info[i].key;
        const char *tag = bcf_hdr_int2id(h, BCF_DT_ID, key);
        (void) bcf_get_info_id(v, key);
        if (!tag)
            continue;
        (void) bcf_hdr_id2int(h, BCF_DT_ID, tag);
        if (bcf_hdr_idinfo_exists(h, BCF_HL_INFO, key)) {
            (void) bcf_hdr_id2type(h, BCF_HL_INFO, key);
            (void) bcf_hdr_id2length(h, BCF_HL_INFO, key);
            (void) bcf_hdr_id2number(h, BCF_HL_INFO, key);
            (void) bcf_hdr_get_hrec(h, BCF_HL_INFO, "ID", tag, NULL);
        }
        query_info(c, v, tag);
    }

    for (i = 0; i < v->n_fmt; i++) {
        int id = v->d.fmt[i].id;
        const char *tag = bcf_hdr_int2id(h, BCF_DT_ID, id);
        (void) bcf_get_fmt_id(v, id);
        if (!tag)
            continue;
        if (bcf_hdr_idinfo_exists(h, BCF_HL_FMT, id)) {
            (void) bcf_hdr_id2type(h, BCF_HL_FMT, id);
            (void) bcf_hdr_id2length(h, BCF_HL_FMT, id);
            (void) bcf_hdr_id2number(h, BCF_HL_FMT, id);
            (void) bcf_hdr_get_hrec(h, BCF_HL_FMT, "ID", tag, NULL);
        }
        query_fmt(c, v, tag);
    }

    query_gt(c, h, v);
    calc_ac(c, h, v);
    format_rec(c, h, v);

    //mutate_info(c, v);
    //mutate_fmt(c, v);
    //add_hdr_tags(c, v);
    //mutate_alleles(c, v);
    //mutate_filters(c, v);
    //mutate_id(c, v);
    round_trip(c, v);
    translate(c, v);
    subset(c, v);
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

    if (hts_get_format(in)->category != variant_data) {
        hts_close(in);
        return NULL;
    }
    return in;
}

// Htslib VCF handling currently has bad error catching, often just giving up
// and calling exit.  This breaks the fuzzer, so for now we catch them and
// keep going.
// Long term we need to fix these, but they're getting in the way of detecting
// real memory corruptions.
//
// Do the final link using -Wl,--wrap=exit -DFUZZ_IGNORE_EXIT
//
// Run with ASAN_OPTIONS=allow_addr2line=1:abort_on_error=1:detect_leaks=0 ./vcf_api_fuzzer -detect_leaks=0 crash-6908485ae0343db75cff385ee2b67fbdb51539b8
// The detect_leaks=0 (both ASAN and fuzzer) is necessary as the exit will be leaking.

#define FUZZ_IGNORE_EXIT

#ifdef FUZZ_IGNORE_EXIT
#include <setjmp.h>

static jmp_buf fuzz_exit_jmp;
static int fuzz_exit_armed;

extern void __sanitizer_print_stack_trace(void) __attribute__((weak));
_Noreturn void __real_exit(int status);

_Noreturn void __wrap_exit(int status) {
    if (fuzz_exit_armed) {
        // Report each distinct call site once, as a to-do list
        static void *seen[64];
        static int nseen;
        void *site = __builtin_return_address(0);
        int i;
        for (i = 0; i < nseen && seen[i] != site; i++)
            ;
        if (i == nseen && nseen < 64) {
            seen[nseen++] = site;
            fprintf(stderr, "htslib called exit(%d); continuing\n", status);
            if (__sanitizer_print_stack_trace)
                __sanitizer_print_stack_trace();
        }
        fuzz_exit_armed = 0;
        longjmp(fuzz_exit_jmp, 1);
    }
    __real_exit(status);
}
#endif

// Mpileup calls exit() directly.  This is problamtic for the fuzzer, but
// we can compile it with -Dexit=fuzz_exit and write our own fuzz_exit
// function that unwinds back to LLVMFuzzerTestOneInput via setjmp
static jmp_buf fuzz_jmp;
void fuzz_exit(int status) {
    fprintf(stderr, "Caught exit of status %d\n", status);
    longjmp(fuzz_jmp, status + 1);
}

int fuzz_one(const uint8_t *data, size_t size) {
    fuzz_ctx c;
    htsFile *in;
    bcf1_t *v = NULL;

    if (size == 0)
        return 0;

    hts_set_log_level(HTS_LOG_OFF);

    memset(&c, 0, sizeof(c));
    c.param = 2166136261u;
    for (size_t i = 0; i < size; i++)
        c.param = (c.param ^ data[i]) * 16777619u;

    if (!(in = open_mem(data, size)))
        return 0;

    if (!(c.hdr = bcf_hdr_read(in)))
        goto out;
    if (set_samples(&c) < 0)
        goto out;

    fuzz_header(&c);
    c.merged = make_merged(&c);
    make_subset(&c);
    c.n_id  = c.hdr->n[BCF_DT_ID];
    c.n_ctg = c.hdr->n[BCF_DT_CTG];

    if (!(v = bcf_init()))
        goto out;
    for (int nrec = 0; nrec < MAX_RECORDS && bcf_read(in, c.hdr, v) == 0;
         nrec++)
        fuzz_record(&c, v);

 out:
    if (v)
        bcf_destroy(v);
    free(c.i32);
    free(c.i64);
    free(c.f);
    free(c.str);
    if (c.fstr) {
        free(c.fstr[0]);
        free(c.fstr);
    }
    free(c.ac);
    free(c.imap);
    ks_free(&c.ks);
    ks_free(&c.ks2);
    if (c.hsub)
        bcf_hdr_destroy(c.hsub);
    if (c.merged)
        bcf_hdr_destroy(c.merged);
    if (c.hdr)
        bcf_hdr_destroy(c.hdr);
    hts_close(in);

    return 0;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
#ifdef FUZZ_IGNORE_EXIT
    if (setjmp(fuzz_exit_jmp))
        return 0; // htslib called exit(); anything this run allocated leaks
    fuzz_exit_armed = 1;
#endif

    int ret = fuzz_one(data, size);

#ifdef FUZZ_IGNORE_EXIT
    fuzz_exit_armed = 0;
#endif
    return ret;
}
