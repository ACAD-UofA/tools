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

int
parse_pops(bcf_hdr_t *hdr, char *pops_fn, char *groups[3], int *(*si)[3], int (*ni)[3])
{
	FILE *fp;
	int ret = 0;
	void *tmp;
	int i;
	int idx;

	char *ind;
	char *grp;

	char *p;
	char *line = NULL;
	size_t len = 0;
	int lineno = 0;
	ssize_t read;

	memset(si, 0, sizeof(int *)*3);
	memset(ni, 0, sizeof(int *)*3);

	fp = fopen(pops_fn, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error opening %s: %s\n",
				pops_fn, strerror(errno));
		ret = -1;
		goto err0;
	}

	while ((read = getline(&line, &len, fp)) != -1) {
		lineno++;

		if (len == 0 || line[0] == '#')
			continue;

		if (line[strlen(line)-1] == '\n')
			line[strlen(line)-1] = '\0';

		for (p=line; p-line<len; p++) {
			if (*p == '\t')
				break;
		}
		if (*p != '\t') {
			fprintf(stderr, "Error: %s: line %d: missing tab\n",
					pops_fn, lineno);
			ret = -2;
			goto err1;
		}

		*p = '\0';
		p++;

		ind = line;
		grp = p;

		for (i=0; i<3; i++) {
			if (strcmp(groups[i], grp) == 0)
				break;
		}
		if (i==3)
			// group not in our A,B,C
			continue;

		idx = 2*bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, ind);
		if (idx < 0) {
			fprintf(stderr, "Error: %s: no sample ``%s'' in vcf.\n",
					pops_fn, ind);
			ret = -3;
			goto err1;
		}

		(*ni)[i]++;

		tmp = realloc((*si)[i], sizeof(int) * (*ni)[i]);
		if (tmp == NULL) {
			perror("parse_pops: realloc");
			ret = -4;
			goto err1;
		}
		(*si)[i] = tmp;
		(*si)[i][(*ni)[i] -1] = idx;
	}

	if (!feof(fp)) {
		fprintf(stderr, "Error reading %s: %s\n",
				pops_fn, strerror(errno));
		ret = -4;
		goto err1;
	}

	for (i=0; i<3; i++) {
		if ((*ni)[i] == 0) {
			fprintf(stderr, "Error: %s: no individuals for group ``%s''\n",
					pops_fn, groups[i]);
			ret = -5;
			goto err1;
		}
	}

	ret = 0;
err1:
	if (ret < 0) {
		for (i=0; i<3; i++) {
			if ((*si)[i]) {
				free((*si)[i]);
				(*si)[i] = NULL;
				(*ni)[i] = 0;
			}
		}
	}

	if (line)
		free(line);
	fclose(fp);
err0:
	return ret;
}

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

/*
 * Parse vcf file, extracting genotypes for specified samples.
 *
 * vcf_fn, vcf filename.
 * pops_fn, file containing individual to group mappings.
 * groups[3], list of population group names, for tree tips A, B, C, D.
 * n, counts of the observed gene topologies.
 */
int
count_patterns_from_vcf(char *vcf_fn, char *pops_fn, char *groups[3], uint (*n)[4])
{
	htsFile *fp;
	bcf_hdr_t *hdr;
	bcf1_t *rec;
	int *si[3]; // list of sample ids for each group
	int ni[3]; // number of ids per group
	int ret;
	int ngt, ngt_arr = 0;
	int *gt_arr = NULL;

	fp = hts_open(vcf_fn, "r");
	if (fp == NULL) {
		ret = -1;
		goto err0;
	}

	hdr = bcf_hdr_read(fp);
	if (hdr == NULL) {
		ret = -2;
		goto err1;
	}

	if (parse_pops(hdr, pops_fn, groups, &si, &ni) < 0) {
		ret = -3;
		goto err2;
	}

	rec = bcf_init1();
	if (rec == NULL) {
		perror("bcf_init1: calloc");
		ret = -4;
		goto err2;
	}

	int idx;
	int skip;
	int i, j;
	int ref;
	int n_ref[3], n_alt[3], n_tot[3];

	memset(*n, 0, 3*sizeof(int));

	while (bcf_read1(fp, hdr, rec) >= 0) {
		// only take biallelic SNPs
		if (!bcf_is_snp(rec) || rec->n_allele != 2)
			continue;

		// filter on QUAL
		if (rec->qual < 25)
			continue;

		// ignore transitions
		if ((rec->d.allele[0][0] == 'C' && rec->d.allele[1][0] == 'T')
			|| (rec->d.allele[1][0] == 'C' && rec->d.allele[0][0] == 'T'))
			continue;

		ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
		if (ngt != ngt_arr) {
			fprintf(stderr, "Error: %s: missing genotypes at %s:%d.\n",
					vcf_fn, bcf_seqname(hdr, rec), rec->pos+1);
			ret = -5;
			goto err3;
		}

		skip = 0;
		for (i=0; i<3; i++) {
			n_ref[i] = n_alt[i] = n_tot[i] = 0;
			for (j=0; j<ni[i]; j++) {
				idx = si[i][j];

				if (bcf_gt_is_missing(gt_arr[idx]))
					continue;

				ref = bcf_gt_allele(gt_arr[idx]) == 0;
				ref += bcf_gt_allele(gt_arr[idx+1]) == 0;

				if (ref == 1) {
					skip = 1;
					break;
				}

				n_ref[i] += ref;
				n_alt[i] += 2-ref;
				n_tot[i] += 2;
			}

			if (skip)
				break;

			if (n_tot[i] == 0) {
				// skip missing
				skip = 1;
				break;
			}
		}

		if (skip)
			continue;

		//printf("[%s:%d] ", bcf_seqname(hdr, rec), rec->pos+1);
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
err3:
	if (gt_arr != NULL)
		free(gt_arr);

	bcf_destroy1(rec);
err2:
	bcf_hdr_destroy(hdr);
err1:
	hts_close(fp);
err0:
	return ret;
}


int
main(int argc, char **argv)
{
	uint n[4];

	if (argc != 6) {
		fprintf(stderr, "%s file.vcf pops.txt A B C\n", argv[0]);
		return -1;
	}

	if (count_patterns_from_vcf(argv[1], argv[2], argv+3, &n) < 0)
		return -2;

	printf("Nab=%d, Nbc=%d, Nac=%d, monomorphic=%d\n",
			n[0], n[1], n[2], n[3]);

	return 0;
}

