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
#include <errno.h>

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>


/*
 * Identify the most parsimonious topology from allele frequencies.
 */
int
topo_from_af(int n_refs[3], int n_alts[3], int n_tot[3])
{
	double af[3]; // ref allele frequency
	double ab, bc, ac;
	int i;
	int topo;

	for (i=0; i<3; i++) {
		af[i] = (double)n_refs[i] / n_tot[i];
	}

#define sq(x) ((x)*(x))

	ab = sq(af[0]-af[1]);
	bc = sq(af[1]-af[2]);
	ac = sq(af[0]-af[2]);

	if (ab < bc && ab < ac) {
		// ((AB),C)
		topo = 0;
	} else if (bc < ab && bc < ac) {
		// (A,(BC))
		topo = 1;
	} else if (ac < ab && ac < bc) {
		// ((AC),B)
		topo = 2;
	} else {
		// site is uninformative
		topo = -1;
	}

	//printf("topo=%d, A=%.3lf, B=%.3lf, C=%.3lf, ab=%.3lf, bc=%.3lf, ac=%.3lf\n",
	//		topo, af[0], af[1], af[2], ab, bc, ac);

	return topo;
}

int
rec_filter(bcf1_t *rec)
{
	// only take SNPs
	if (!bcf_is_snp(rec))
		return -1;

	// biallelic sites only
	if (rec->n_allele > 2)
		return -2;

	// filter on QUAL
	if (rec->qual < 25)
		return -3;

	// ignore transitions
	if (rec->n_allele == 2 &&
		((rec->d.allele[0][0] == 'C' && rec->d.allele[1][0] == 'T')
		|| (rec->d.allele[1][0] == 'C' && rec->d.allele[0][0] == 'T')))
		return -4;

	return 0;
}

/*
 * Parse vcf files, extracting genotypes for specified samples.
 *
 * fn[3], vcf filenames.
 * n={Nab, Nbc, Nac, Nx}, counts of the observed pattern.
 */
int
count_patterns_from_vcf(char *fn[3], uint (*n)[4])
{
	htsFile *fp[3] = {NULL,};
	tbx_t *tbx[3] = {NULL,};
	bcf_hdr_t *hdr[3] = {NULL,};
	bcf1_t *rec[3] = {NULL, };

	int ret;
	int i;

	for (i=0; i<3; i++) {
		fp[i] = hts_open(fn[i], "r");
		if (fp[i] == NULL) {
			ret = -1;
			goto err0;
		}
		tbx[i] = tbx_index_load(fn[i]);
		if (tbx[i] == NULL) {
			ret = -2;
			goto err0;
		}
		hdr[i] = bcf_hdr_read(fp[i]);
		if (hdr[i] == NULL) {
			ret = -3;
			goto err0;
		}

		rec[i] = bcf_init1();
		if (rec[i] == NULL) {
			perror("bcf_init1: calloc");
			ret = -4;
			goto err0;
		}
	}

	int ngt, ngt_arr = 0;
	int *gt_arr = NULL;

	int skip;
	int ref;
	int n_ref[3], n_alt[3], n_tot[3];

	hts_itr_t *itr = NULL;
	kstring_t line = {0,};

	memset(*n, 0, 4*sizeof(int));

	while (bcf_read1(fp[0], hdr[0], rec[0]) >= 0) {
		if (rec_filter(rec[0]) < 0)
			continue;

		ngt = bcf_get_genotypes(hdr[0], rec[0], &gt_arr, &ngt_arr);
		if (ngt != ngt_arr) {
			fprintf(stderr, "Error: %s: missing genotypes at %s:%d.\n",
					fn[0], bcf_seqname(hdr[0], rec[0]), rec[0]->pos+1);
			ret = -5;
			goto err1;
		}

		if (bcf_gt_is_missing(gt_arr[0]))
			continue;

		ref = bcf_gt_allele(gt_arr[0]) == 0;
		ref += bcf_gt_allele(gt_arr[1]) == 0;

		n_ref[0] = ref;
		n_alt[0] = 2-ref;
		n_tot[0] = 2;

		skip = 0;
		for (i=1; i<3; i++) {
			itr = tbx_itr_queryi(tbx[i], rec[0]->rid, rec[0]->pos, rec[0]->pos+1);
			if (itr == NULL) {
				skip = 1;
				break;
			}

			skip = 1;
			while (tbx_itr_next(fp[i], tbx[i], itr, &line) >= 0) {
				vcf_parse1(&line, hdr[i], rec[i]);
				free(ks_release(&line));
				if (rec_filter(rec[i]) == 0) {
					skip = 0;
					break;
				}
			}

		        tbx_itr_destroy(itr);

			if (skip) {
				free(ks_release(&line));
				break;
			}
			skip = 0;


			ngt = bcf_get_genotypes(hdr[i], rec[i], &gt_arr, &ngt_arr);
			if (ngt != ngt_arr) {
				fprintf(stderr, "Error: %s: missing genotypes at %s:%d.\n",
						fn[i], bcf_seqname(hdr[i], rec[i]), rec[i]->pos+1);
				ret = -5;
				goto err1;
			}

			if (bcf_gt_is_missing(gt_arr[0])) {
				skip = 1;
				break;
			}

			ref = bcf_gt_allele(gt_arr[0]) == 0;
			ref += bcf_gt_allele(gt_arr[1]) == 0;

			n_ref[i] = ref;
			n_alt[i] = 2-ref;
			n_tot[i] = 2;
		}

		if (skip)
			continue;

		//printf("[%s:%d]\n", bcf_seqname(hdr[0], rec[0]), rec[0]->pos+1);
		switch (topo_from_af(n_ref, n_alt, n_tot)) {
			case 0:
				(*n)[0]++;
				break;
			case 1:
				(*n)[1]++;
				break;
			case 2:
				(*n)[2]++;
				break;
			case -1:
				// monomorphic
				(*n)[3]++;
				break;
		}

	}

	ret = 0;
err1:
	if (gt_arr != NULL)
		free(gt_arr);
	//free(line.s);

err0:
	for (i=0; i<3; i++) {
		if (rec[i])
			bcf_destroy1(rec[i]);
		if (hdr[i])
			bcf_hdr_destroy(hdr[i]);
		if (tbx[i])
			tbx_destroy(tbx[i]);
		if (fp[i])
			hts_close(fp[i]);
	}
	return ret;
}


int
main(int argc, char **argv)
{
	uint n[4];

	if (argc != 4) {
		fprintf(stderr, "%s A.vcf B.vcf C.vcf\n", argv[0]);
		return -1;
	}

	if (count_patterns_from_vcf(argv+1, &n) < 0)
		return -2;

	printf("Nab=%d, Nbc=%d, Nac=%d, monomorphic=%d\n",
			n[0], n[1], n[2], n[3]);

	return 0;
}

