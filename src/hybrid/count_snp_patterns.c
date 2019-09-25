/*
 * Copyright (c) 2016 Graham Gower <graham.gower@gmail.com>
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
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <unistd.h>
#include <getopt.h>
#include <signal.h>

#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>


typedef struct {
	char *fn[3]; // filenames
	char *sn[3]; // sample names

	int ld_x; // distance required for loci to be unlinked
	int ignore_transitions;
} opt_t;


typedef struct {
	int ld_x;
	int last_pos;
	uint n_ab, n_bc, n_ac;
} ld_aux;

static int stopped = 0;

void
quit_handler(int sig)
{
	stopped = 1;
}

int
count_snp_patterns(opt_t *opt, ld_aux *lda, int n)
{
	bcf_srs_t *sr;
	bcf_hdr_t *hdr[3];
	bcf1_t *rec[3];
	int ret;
	int i;
	int k;
	int nr;
       
	int ngt_arr = 0;
	int *gt_arr = NULL; // genotype array
	int si[3]; // genotype array indices for each sample
	int gt[3]; // genotypes (number of REF alleles for each sample)
	char ref, alt;

	int last_rid = -1;

	sr = bcf_sr_init();
	if (sr == NULL) {
		ret = -1;
		goto err0;
	}

	sr->require_index = 1;
	sr->collapse = COLLAPSE_ANY;

	for (i=0; i<3; i++)
		bcf_sr_add_reader(sr, opt->fn[i]);

	for (i=0; i<3; i++) {
		hdr[i] = bcf_sr_get_header(sr, i);

		// get sample indicies, if any
		if (opt->sn[i]) {
			si[i] = 2*bcf_hdr_id2int(bcf_sr_get_header(sr, i), BCF_DT_SAMPLE, opt->sn[i]);
			if (si[i] < 0) {
				fprintf(stderr, "Error: %s: no sample ``%s''.\n",
					opt->fn[i], opt->sn[i]);
				ret = -2;
				goto err1;
			}
		} else {
			si[i] = 0;
		}
	}

	while ((nr = bcf_sr_next_line(sr)) && !stopped) {
		if (nr != sr->nreaders)
			continue;

		rec[0] = bcf_sr_get_line(sr, 0);

		//printf("[%s:%d] nr=%d\n", bcf_seqname(hdr[0], rec[0]), rec[0]->pos+1, nr);

		if (last_rid != rec[0]->rid) {
			for (k=0; k<n; k++)
				lda[k].last_pos = -1;
		} else {
			for (k=0; k<n; k++) {
				if (lda[k].last_pos == -1 || rec[0]->pos >= lda[k].last_pos+lda[k].ld_x)
					break;
			}
			if (k == n)
				continue;
		}

		ref = rec[0]->d.allele[0][0];
		alt = '.';

		rec[1] = bcf_sr_get_line(sr, 1);
		rec[2] = bcf_sr_get_line(sr, 2);

		for (i=0; i<3; i++) {

			// filter on QUAL
			if (rec[i]->qual < 25)
				break;

			if (bcf_get_genotypes(hdr[i], rec[i], &gt_arr, &ngt_arr) < 1) {
				fprintf(stderr, "Error: %s: no FORMAT/GT field at %s:%d.\n",
						opt->fn[0], bcf_seqname(hdr[i], rec[i]), rec[i]->pos+1);
				ret = -3;
				goto err2;
			}

			if (bcf_gt_is_missing(gt_arr[si[i]]))
				break;

			int a = bcf_gt_allele(gt_arr[si[i]]);
			int b = bcf_gt_allele(gt_arr[si[i]+1]);

			// skip heterozygotes
			if (a != b)
				break;

			if (a > 0) {
				// only take SNPs
				if (bcf_get_variant_type(rec[i], a) != VCF_SNP)
					break;

				char alt_i = rec[i]->d.allele[a][0];
				if (alt == '.')
					alt = alt_i;
				else if (alt != alt_i) {
					// different ALT genotypes called in different samples
					break;
				}
			}

			// number of REF alleles
			gt[i] = 2 - a - b;
		}

		if (i != 3) {
			// filtered
			continue;
		}

		// skip monomorphic sites
		if (gt[0] == gt[1] && gt[1] == gt[2])
			continue;

		// filter transitions
		if (opt->ignore_transitions &&
		   ((ref == 'C' && alt == 'T') || (alt == 'C' && ref == 'T') ||
		   (ref == 'G' && alt == 'A') || (alt == 'G' && ref == 'A')))
			continue;

		/*
		 * Got an informative site.
		 */
		//printf("[%s:%d] gt={%d,%d,%d}\n", bcf_seqname(hdr[0], rec[0]), rec[0]->pos+1, gt[0], gt[1], gt[2]);

		last_rid = rec[0]->rid;

		for (k=0; k<n; k++) {
			if (lda[k].last_pos != -1 && lda[k].last_pos+lda[k].ld_x > rec[0]->pos)
				continue;

			lda[k].last_pos = rec[0]->pos;

			if (gt[0] == gt[1] && gt[0] != gt[2])
				lda[k].n_ab++;
			else if (gt[1] == gt[2] && gt[0] != gt[1])
				lda[k].n_bc++;
			else if (gt[0] == gt[2] && gt[0] != gt[1])
				lda[k].n_ac++;
		}
	}

	ret = 0;
err2:
	if (gt_arr != NULL)
		free(gt_arr);
err1:
        bcf_sr_destroy(sr);
err0:
	return ret;
}

int
strsplit(char *str, char sep, char ***out)
{
	char *s, *start;
	char **tmp;
	int i, n;

	s = strdup(str);
	if (s == NULL) {
		perror("strsplit: strdup");
		return -1;
	}

	*out = NULL;

	n = 0;
	start = s;
	for (i=0; s[i]; i++) {
		if (s[i] == sep) {
			s[i] = 0;
			n++;
			tmp = realloc(*out, n*sizeof(char *));
			if (tmp == NULL) {
				if (*out)
					free(*out);
				*out = NULL;
				free(s);
				perror("strsplit: realloc");
				return -2;
			}
			tmp[n-1] = strdup(start);
			*out = tmp;
			start = s + i + 1;
		}
	}

	n++;
	tmp = realloc(*out, n*sizeof(char *));
	if (tmp == NULL) {
		if (*out)
			free(*out);
		*out = NULL;
		free(s);
		perror("strsplit: realloc");
		return -2;
	}
	tmp[n-1] = strdup(start);
	*out = tmp;

	free(s);

	return n;
}

int
lda_compar(const void *p1, const void *p2)
{
	const ld_aux *l1 = p1;
	const ld_aux *l2 = p2;
	return l1->ld_x - l2->ld_x;
}

void
usage(char *argv0)
{
	fprintf(stderr, "%s [-t] [-x X] A.vcf[:X] B.vcf[:Y] C.vcf[:Z]\n", argv0);
	exit(1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	ld_aux *lda;
	int n = 0;
	int c;
	int i;

	opt.ignore_transitions = 0;

	while ((c = getopt(argc, argv, "tx:")) != -1) {
		switch (c) {
			case 't':
				opt.ignore_transitions = 1;
				break;
			case 'x':
				{
				char **xx;
				n = strsplit(optarg, ',', &xx);
				lda = calloc(n, sizeof(*lda));
				if (lda == NULL) {
					perror("main: calloc");
					for (i=0; i<n; i++) {
						free(xx[i]);
					}
					free(xx);
					return -1;
				}
				for (i=0; i<n; i++) {
					int mul = 1;
					char *x = xx[i];
					int end = strlen(x)-1;
					if (x[end] == 'k' || x[end] == 'K') {
						mul = 1e3;
					}
					else if (x[end] == 'm' || x[end] == 'M') {
						mul = 1e6;
					}
					lda[i].ld_x = mul * strtol(xx[i], NULL, 0);
					if (lda[i].ld_x < 1 || lda[i].ld_x > 1e7) {
						fprintf(stderr, "-x X=`%s' out of bounds", x);
						for (i=0; i<n; i++) {
							free(xx[i]);
						}
						free(xx);
						usage(argv[0]);
					}
					free(x);
				}
				free(xx);
				qsort(lda, n, sizeof(*lda), lda_compar);
				}
				break;
			default:
				usage(argv[0]);
		}
	}

	if (argc-optind != 3) {
		usage(argv[0]);
	}


	for (i=0; i<3; i++) {
		char *p = strrchr(argv[optind+i], ':');
		opt.fn[i] = argv[optind+i];
		if (p == NULL) {
			opt.sn[i] = NULL;
		} else {
			*p = '\0';
			opt.sn[i] = p+1;
		}
	}

	if (n == 0) {
		n = 1;
		lda = calloc(n, sizeof(*lda));
		if (lda == NULL) {
			perror("main: calloc");
			return -1;
		}
		lda[0].ld_x = 0;
	}

	for (i=0; i<n; i++) {
		//printf("lda[%d].ld_x=%d\n", i, lda[i].ld_x);
		lda[i].n_ab = lda[i].n_bc = lda[i].n_ac = 0;
		lda[i].last_pos = -1;
	}

	signal(SIGINT, quit_handler);
	signal(SIGTERM, quit_handler);

	if (count_snp_patterns(&opt, lda, n) < 0)
		return -2;

	printf("LD_X\tNab\tNbc\tNac\n");
	for (i=0; i<n; i++) {
		printf("%d\t%d\t%d\t%d\n",
			lda[i].ld_x, lda[i].n_ab, lda[i].n_bc, lda[i].n_ac);
	}
	free(lda);

	return 0;
}

