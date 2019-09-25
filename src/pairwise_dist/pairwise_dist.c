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
#include <getopt.h>
#include <errno.h>
#include <time.h>

#include <htslib/vcf.h>
#include <htslib/khash.h>

KHASH_SET_INIT_INT(int)

typedef struct {
	int min_qual;
	float depth;
	char *contigs_fn;
} opt_t;

khash_t(int) *
parse_contigs(opt_t *opt, bcf_hdr_t *hdr)
{
	khash_t(int) *h;
	FILE *fp;
	char buf[BUFSIZ];
	int i;

	fp = fopen(opt->contigs_fn, "r");
	if (fp == NULL) {
		fprintf(stderr, "%s: fopen: %s\n", opt->contigs_fn, strerror(errno));
		return NULL;
	}

	h = kh_init(int);
	if (h == NULL) {
		perror("kh_init (calloc)");
		fclose(fp);
		return NULL;
	}

	while (fgets(buf, BUFSIZ, fp)) {
		int len = strlen(buf);

		while (len > 0 && (buf[len-1] == '\r' || buf[len-1] == '\n'))
			buf[--len] = '\0';

		if (len == 0)
			continue;

		for (i=0; i<bcf_hdr_nsamples(hdr); i++) {
			khint32_t rid = bcf_hdr_name2id(hdr, buf);
			if (rid != -1) {
				int unused;
				kh_put(int, h, rid, &unused);
				break;
			}
		}
	}

	if (!feof(fp)) {
		fprintf(stderr, "%s: fgets: %s\n", opt->contigs_fn, strerror(errno));
		kh_destroy(int, h);
		h = NULL;
	}

	fclose(fp);

	return h;
}

int
pairwise_dist(opt_t *opt, char *vcf_fn)
{
	htsFile *fp;
	bcf_hdr_t *hdr;
	bcf1_t *rec;
	khash_t(int) *ridset = NULL;
	int ret;
	int i;

	int n_dp = 0;
	int *dp = NULL;
	int ndpr_arr = 0;
	int32_t *dpr_arr = NULL; // DPR array

	float min_depth = opt->depth/2.0;
	float max_depth = opt->depth*2.0;

	uint64_t mm = 0;
	uint64_t sites = 0;

	unsigned short xsubi[3] = {31,41,59}; // random state

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

	if (bcf_hdr_nsamples(hdr) != 1) {
		fprintf(stderr, "Error: %s: Multisample bcf/vcf\n", vcf_fn);
		ret = -3;
		goto err2;
	}

	if (opt->contigs_fn) {
		ridset = parse_contigs(opt, hdr);
		if (ridset == NULL) {
			ret = -4;
			goto err2;
		}
	}

	rec = bcf_init1();
	if (rec == NULL) {
		perror("bcf_init1: calloc");
		ret = -5;
		goto err2;
	}

	while (bcf_read1(fp, hdr, rec) >= 0) {

		if (opt->contigs_fn && kh_get(int, ridset, rec->rid) == kh_end(ridset))
			// filter by contig
			continue;

		if (rec->qual < opt->min_qual)
			continue;

		if (bcf_get_format_int32(hdr, rec, "DP", &dp, &n_dp) < 0) {
			fprintf(stderr, "%s: missing FORMAT/DP field at %s:%d\n",
					vcf_fn, bcf_seqname(hdr, rec), rec->pos+1);
			ret = -6;
			goto err3;
		}

		if (dp[0] < min_depth || dp[0] > max_depth)
			continue;

#if 0
		if (bcf_get_genotypes(hdr, rec, &dpr_arr, &ndpr_arr) < 1) {
                        fprintf(stderr, "%s: missing FORMAT/GT field at %s:%d\n",
                                        vcf_fn, bcf_seqname(hdr, rec), rec->pos+1);
                        ret = -7;
                        goto err3;
                }

                int gt = bcf_gt_allele(dpr_arr[nrand48(xsubi)%1]);
		if (gt)
			mm ++;
		sites++;
#else
		if (bcf_get_format_int32(hdr, rec, "DPR", &dpr_arr, &ndpr_arr) < 1) {
			// DPR is deprecated, try AD
			if (bcf_get_format_int32(hdr, rec, "AD", &dpr_arr, &ndpr_arr) < 1) {
				fprintf(stderr, "Error: %s: no FORMAT/DPR nor FORMAT/AD field at %s:%d.\n",
					vcf_fn, bcf_seqname(hdr, rec), rec->pos+1);
				ret = -7;
				goto err3;
			}
		}

		int32_t dp_ref = dpr_arr[0];
		int32_t dp_alt = 0;
		int gt;

		for (i=1; i<rec->n_allele; i++)
			dp_alt += dpr_arr[i];

		if (dp_alt == 0)
			gt = 0;
		else
			gt = (nrand48(xsubi) % (dp_ref + dp_alt) < dp_ref) ? 0 : 1;

		mm += gt;
		sites++;
#endif
	}

	ret = 0;

	printf("mm=%jd, sites=%jd, mu=%lg\n", mm, sites, (double)mm/sites);

err3:
	if (dp)
		free(dp);
	if (dpr_arr)
		free(dpr_arr);
	bcf_destroy1(rec);
	if (opt->contigs_fn)
		kh_destroy(int, ridset);
err2:
	bcf_hdr_destroy(hdr);
err1:
	hts_close(fp);
err0:
	return ret;
}

void
usage(char *argv0)
{
	fprintf(stderr, "usage: %s [-q MIN_QUAL] [-c contigs.txt] -d MEAN_DEPTH file.{b,v}cf[.gz]\n", argv0);
	exit(1);
}

int
main(int argc, char **argv)
{
	int c;
	char *vcf_fn;
	opt_t opt;

	opt.min_qual = 25;
	opt.depth = -1;
	opt.contigs_fn = NULL;

	while ((c = getopt(argc, argv, "q:d:c:")) != -1) {
		switch (c) {
			case 'q':
				opt.min_qual = strtoul(optarg, NULL, 0);
				if (opt.min_qual < 0 || opt.min_qual > 1000) {
					fprintf(stderr, "min_qual=`%s' out of range\n", optarg);
					return -1;
				}
				break;
			case 'd':
				opt.depth = strtof(optarg, NULL);
				if (opt.depth < 1 || opt.depth > 1000) {
					fprintf(stderr, "depth=`%s' out of range\n", optarg);
					return -1;
				}
				break;
			case 'c':
				opt.contigs_fn = optarg;
				break;
			default:
				usage(argv[0]);
		}
	}

	if (argc-optind != 1 || opt.depth == -1) {
		usage(argv[0]);
	}

	vcf_fn = argv[optind];

	if (pairwise_dist(&opt, vcf_fn) < 0) {
		return -1;
	}

	return 0;
}
