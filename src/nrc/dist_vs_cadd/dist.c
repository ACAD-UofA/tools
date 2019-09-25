/*
 * Copyright (c) 2015, 2016 Graham Gower <graham.gower@gmail.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <errno.h>
#include <signal.h>

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>

// impact
#define I_MODIFIER 0
#define I_LOW 1
#define I_MODERATE 2
#define I_HIGH 3

typedef struct {
	uint32_t d; // pairwise distances
	uint32_t n; // number of loci compared
	uint32_t impact[4];
	uint32_t cadd[100];
} pair_stats_t;

// global data
struct {
	bcf_hdr_t *hdr;
	pair_stats_t *pair_stats;
	unsigned int n_pairs;

	// Genotype array, realloc()'d in bcf_get_genotypes().
	int *gt_arr;
	int *gt_idx;
	tbx_t *cadd_tbx;

	htsFile *cadd_fp;

	int quit_signalled;
} gd;



/*
 * Return the snpEff impact, or -1 if not found.
 * Multiple entries may be found, e.g. for alternate splicing or
 * transcription on both strands.  In this case, take the highest impact.
 */
int
get_snpeff_impact(bcf1_t *rec)
{
	const char *c;
	int n_pipes = 0;
	int impact = -1;
	char *ann = NULL;
	int ann_ndst = 0;

	if (bcf_get_info_string(gd.hdr, rec, "ANN", &ann, &ann_ndst) < 0) {
		/*
		fprintf(stderr, "Error: missing INFO/ANN field for %s:%d.\n",
				bcf_seqname(gd.hdr, rec), rec->pos+1);
		exit(1);
		*/
		return -1;
	}

	for (c=ann; *c!='\0'; c++) {
		switch (*c) {
			case '|':
				n_pipes++;
				if (n_pipes == 2) {
					int tmp = -1;
					if (strncmp(c+1, "LOW", 3))
						tmp = I_LOW;
					else if(strncmp(c+1, "MODIFIER", 8))
						tmp = I_MODIFIER;
					else if(strncmp(c+1, "MODERATE", 8))
						tmp = I_MODERATE;
					else if (strncmp(c+1, "HIGH", 4))
						tmp = I_LOW;
					if (tmp > impact)
						/* take highest impact */
						impact = tmp;
				}
				break;
			case ',':
				/* try the next annotation */
				n_pipes = 0;
				break;
		}

	}

	free(ann);

	return impact;
}

/*
 * Return pointer to field x in a tab delimited string.
 */
const char *
get_tab_field(kstring_t str, int x)
{
	char *p, *end;
	int n;

	if (x < 1)
		return NULL;
	if (x == 1)
		return str.s;

	n = 1;
	end = str.s + str.l;

	for (p=str.s; p<end; p++) {
		if (*p == '\t') {
			if (++n == x) {
				return ++p;
			}
		}
	}
	
	return NULL;
}

int
get_cadd_score(int tid, int pos)
{
	kstring_t str = {0,};
	hts_itr_t *itr;
	const char *s;
	float cadd_f;
	int ret = 0;

	// just use the first entry, scores are the same
	itr = tbx_itr_queryi(gd.cadd_tbx, tid, pos, pos+1);
	if (itr == NULL || (ret = tbx_itr_next(gd.cadd_fp, gd.cadd_tbx, itr, &str)) < 0) {
		fprintf(stderr, "Error: couldn't find %s:%d in CADD index (itr=0x%p, ret=%d).\n",
				bcf_hdr_id2name(gd.hdr, tid), pos+1, itr, ret);
		exit(1);
	}

        tbx_itr_destroy(itr);

	// CADD score is in the 116th column.
	// TODO: parameterise the column number for other fields
	s = get_tab_field(str, 116);
	if (s == NULL) {
		fprintf(stderr, "Error: couldn't get CADD score (column 116) for %s:%d\n",
				bcf_hdr_id2name(gd.hdr, tid), pos);
		exit(1);
	}

	errno = 0;
	cadd_f = strtof(s, NULL);
	if (errno || (cadd_f < 0.0 || cadd_f >= 100.0)) {
		fprintf(stderr, "Error: invalid CADD score (column 116) for %s:%d: ``%s''\n",
				bcf_hdr_id2name(gd.hdr, tid), pos, s);
		exit(1);
	}

	free(str.s);

	//printf("CADD[%s:%d]=%f\n", bcf_hdr_id2name(gd.hdr, tid), pos+1, cadd_f);

	return (int) cadd_f;
}

/*
 * Called for each VCF record.
 * Return -1 for errors, 0 otherwise.
 */
int
process(bcf1_t *rec)
{
	int i, j;
	int idx, jdx; // gt indices for i and j
	unsigned int n = 0;
	int ngt, ngt_arr = 0;
	unsigned int gt_i0, gt_i1, gt_j0, gt_j1;

	int cadd = -1;
	int impact = -1;

	// we only want biallelic SNPs
	if (!bcf_is_snp(rec) || rec->n_allele != 2)
		return 0;

	ngt = bcf_get_genotypes(gd.hdr, rec, &gd.gt_arr, &ngt_arr);
	if (ngt != ngt_arr) {
		fprintf(stderr, "Failed to get genotypes for %s:%d.\n",
				bcf_seqname(gd.hdr, rec), rec->pos+1);
		return -1;
	}

	//impact = get_snpeff_impact(rec);
	cadd = get_cadd_score(rec->rid, rec->pos);

	for (i=0; i<bcf_hdr_nsamples(gd.hdr); i++) {
		idx = gd.gt_idx[i];
		if (bcf_gt_is_missing(gd.gt_arr[idx]) /*|| bcf_gt_is_missing(gd.gt_arr[idx+1])*/) {
			//printf("missing genotype for sample[%d]=%s, idx=%d\n", i, gd.hdr->samples[i], idx);
			n += i;
			continue;
		}

		gt_i0 = bcf_gt_allele(gd.gt_arr[idx]);
		gt_i1 = bcf_gt_allele(gd.gt_arr[idx+1]);

		// skip heterozygote
		if (gt_i0 != gt_i1) {
			n += i;
			continue;
		}

		//printf("sample[%d]=%s, GT=%d/%d, idx=%d\n", i, gd.hdr->samples[i], gt_i0, gt_i1, idx);

		for (j=0; j<i; j++, n++) {
			jdx = gd.gt_idx[j];
			if (bcf_gt_is_missing(gd.gt_arr[jdx]) /*|| bcf_gt_is_missing(gd.gt_arr[jdx+1])*/) {
				continue;
			}

			gt_j0 = bcf_gt_allele(gd.gt_arr[jdx]);
			gt_j1 = bcf_gt_allele(gd.gt_arr[jdx+1]);

			// skip heterozygote
			if (gt_j0 != gt_j1) {
				continue;
			}

			if (gt_i0 != gt_j0) {
				gd.pair_stats[n].d++;
				if (impact != -1)
					gd.pair_stats[n].impact[impact]++;
				if (cadd != -1)
					gd.pair_stats[n].cadd[cadd]++;
			}
			gd.pair_stats[n].n++;
		}
	}

	return 0;
}

/* Comparison function for qsort.
 * Sort by pairwise distance, normalised by number of sites compared.
 */
static int
pair_cmp(const void *p1, const void *p2)
{
	const pair_stats_t *pp1 = p1;
	const pair_stats_t *pp2 = p2;
	return pp1->d * pp2->n - pp2->d * pp1->n;
}

void
handle_sigint(int sig)
{
	gd.quit_signalled = 1;
}


int
main(int argc, char **argv)
{
	int i, j, k;
	unsigned int n;
	bcf1_t *rec;
	htsFile *fp;

	if (argc != 3) {
		fprintf(stderr, "usage %s file.{vcf|bcf}[.gz] 1000G_phase3_inclAnno.tsv.gz\n", argv[0]);
		return -1;
	}

	fp = hts_open(argv[1], "r");
	if (fp == NULL)
		return -1;

	gd.hdr = bcf_hdr_read(fp);
	if (gd.hdr == NULL) {
		hts_close(fp);
		return -1;
	}

	gd.n_pairs = bcf_hdr_nsamples(gd.hdr) * (bcf_hdr_nsamples(gd.hdr)-1) / 2;
	gd.pair_stats = calloc(gd.n_pairs, sizeof(pair_stats_t));
	gd.gt_idx = calloc(bcf_hdr_nsamples(gd.hdr), sizeof(int));
	if (gd.pair_stats == NULL || gd.gt_idx == NULL) {
		perror("calloc");
		hts_close(fp);
		return -1;
	}

	gd.gt_arr = NULL;

	for (i=0; i<bcf_hdr_nsamples(gd.hdr); i++)
		gd.gt_idx[i] = 2*bcf_hdr_id2int(gd.hdr, BCF_DT_SAMPLE, gd.hdr->samples[i]);


	gd.cadd_fp = hts_open(argv[2], "r");
	if (gd.cadd_fp == NULL) {
		hts_close(fp);
		return -1;
	}

	gd.cadd_tbx = tbx_index_load(argv[2]);
	if (gd.cadd_tbx == NULL) {
		hts_close(gd.cadd_fp);
		hts_close(fp);
		return -1;
	}

	gd.quit_signalled = 0;
	signal(SIGINT, handle_sigint);

	rec = bcf_init1();

	while (bcf_read1(fp, gd.hdr, rec) >= 0 && !gd.quit_signalled) {
		if (process(rec) < 0)
			break;
	}

	bcf_destroy1(rec);
	hts_close(fp);
	hts_close(gd.cadd_fp);
	tbx_destroy(gd.cadd_tbx);

	/* sort by pairwise distance */
	qsort(gd.pair_stats, gd.n_pairs, sizeof(pair_stats_t), pair_cmp);

	n = 0;
	printf("sample1\tsample2\tdistance\tn_sites\timpact_MODIFIER\timpact_LOW\timpact_MODERATE\timpact_HIGH");
	for(k=0; k<100; k++)
		printf("\tCADD%d", k);
	printf("\n");
	for (i=0; i<bcf_hdr_nsamples(gd.hdr); i++) {
		for (j=0; j<i; j++) {
			printf("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d",
					gd.hdr->samples[i], gd.hdr->samples[j],
					gd.pair_stats[n].d, gd.pair_stats[n].n,
					gd.pair_stats[n].impact[I_MODIFIER],
					gd.pair_stats[n].impact[I_LOW],
					gd.pair_stats[n].impact[I_MODERATE],
					gd.pair_stats[n].impact[I_HIGH]);
			for(k=0; k<100; k++)
				printf("\t%d", gd.pair_stats[n].cadd[k]);
			printf("\n");
			n++;
		}
	}

	free(gd.pair_stats);
	free(gd.gt_arr);
	free(gd.gt_idx);


	bcf_hdr_destroy(gd.hdr);

	return 0;
}
