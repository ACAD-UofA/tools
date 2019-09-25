/*
 * Copyright (c) 2016,2017 Graham Gower <graham.gower@gmail.com>
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
#include <ctype.h>

#include <htslib/vcf.h>


// This should be in htslib, but I couldn't find it.
static inline uint32_t bcf_hdr_id2ctglen(const bcf_hdr_t *hdr, int rid) { return hdr->id[BCF_DT_CTG][rid].val->info[0]; }

#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

typedef struct {
	int min_depth;
	int max_depth;
	int window;
	int step;
} opt_t;

typedef struct block {
	uint32_t hom;
	uint32_t n;
} block_t;

int
het_windows(opt_t *opt, char *vcf_fn)
{
	htsFile *fp;
	bcf_hdr_t *hdr;
	bcf1_t *rec;
	int ret;

	int n_dp = 0;
	int *dp = NULL;
	int ngt_arr = 0;
	int32_t *gt_arr = NULL;

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
		fprintf(stderr, "Warning: %s: using first sample of multisample bcf/vcf.\n", vcf_fn);
		// TODO: call bcf_hdr_set_samples() to only parse sample of interest.
	}


	rec = bcf_init1();
	if (rec == NULL) {
		perror("bcf_init1: calloc");
		ret = -3;
		goto err2;
	}

	int32_t last_chrom = -1, win_start = 0, chunk_start = 0;
	int n_chunks = opt->window / opt->step;
	block_t *blocklist;
	block_t *b = blocklist;
	int bi = 0, bcount = 0;
	int window_n = 0, window_hom = 0;

	b = blocklist = calloc(n_chunks, sizeof(*blocklist));
	if (blocklist == NULL) {
		ret = -4;
		goto err3;
	}

	printf("CHROM\tSTART\tEND\tN\tHOM\n");

	while (bcf_read(fp, hdr, rec) >= 0) {

		int var_type = bcf_get_variant_types(rec);
		if (var_type != VCF_REF && var_type != VCF_SNP)
			continue;

		if (bcf_get_format_int32(hdr, rec, "DP", &dp, &n_dp) < 0) {
			fprintf(stderr, "%s: missing FORMAT/DP field at %s:%d\n",
					vcf_fn, bcf_seqname(hdr, rec), rec->pos+1);
			ret = -5;
			goto err4;
		}

		if (dp[0] < opt->min_depth || dp[0] > opt->max_depth)
			continue;

		if (bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr) < 1) {
                        fprintf(stderr, "%s: missing FORMAT/GT field at %s:%d\n",
                                        vcf_fn, bcf_seqname(hdr, rec), rec->pos+1);
                        ret = -6;
                        goto err4;
                }

		int hom = (bcf_gt_allele(gt_arr[0]) == bcf_gt_allele(gt_arr[1]));

		if (rec->rid != last_chrom || rec->pos > chunk_start+opt->step) {
			window_n += b->n;
			window_hom += b->hom;

			if (rec->rid != last_chrom) {
				if (last_chrom != -1) {
					printf("%s\t%d\t%d\t%d\t%d\n",
							bcf_hdr_id2name(hdr, last_chrom),
							win_start,
							min(win_start+opt->window, bcf_hdr_id2ctglen(hdr,last_chrom)),
							window_n,
							window_hom);
				}
				bi = 0;
				b = blocklist;
				bcount = 0;
				window_n = 0;
				window_hom = 0;
				chunk_start = 0;
				win_start = 0;
				memset(blocklist, 0, n_chunks*sizeof(*blocklist));
				last_chrom = rec->rid;
				continue;

			}

			do {
				bi = (bi+1) % n_chunks;
				b = blocklist+bi;
				bcount += (bcount<n_chunks?1:0);
				
				if (rec->pos > win_start+opt->window) {
					printf("%s\t%d\t%d\t%d\t%d\n",
							bcf_hdr_id2name(hdr, last_chrom),
							win_start,
							min(win_start+opt->window, bcf_hdr_id2ctglen(hdr,last_chrom)),
							window_n,
							window_hom);

					window_n -= b->n;
					window_hom -= b->hom;
					assert(window_n >= 0);
					assert(window_hom >= 0);
					b->n = 0;
					b->hom = 0;
					win_start += opt->step;
				}
			} while (rec->rid == last_chrom && rec->pos > win_start+opt->window);

			chunk_start += opt->step;
		}

		b->n++;
		b->hom += hom;
	}

	// mop up
	if (rec->pos > win_start) {
		if (rec->pos > chunk_start) {
			window_n += b->n;
			window_hom += b->hom;
		}
		printf("%s\t%d\t%d\t%d\t%d\n",
			bcf_hdr_id2name(hdr, last_chrom),
			win_start,
			min(win_start+opt->window, bcf_hdr_id2ctglen(hdr,last_chrom)),
			window_n,
			window_hom);
	}

	ret = 0;


err4:
	if (dp)
		free(dp);
	if (gt_arr)
		free(gt_arr);
	free(blocklist);
err3:
	bcf_destroy1(rec);
err2:
	bcf_hdr_destroy(hdr);
err1:
	hts_close(fp);
err0:
	return ret;
}

unsigned long
parse_bp(char *s)
{
	unsigned long x = 0;

again:
	if (strlen(s) <= 1)
		return -1;

	switch (tolower(s[strlen(s)-1])) {
		case 'k':
			x = 1000;
			break;
		case 'm':
			x = 1000000;
			break;
		case 'b':
			if (x == -1)
				return -1;
			x = -1;
			s[strlen(s)-1] = 0;
			goto again;
		default:
			x = 1;
			break;
	}

	return x*strtoul(s, NULL, 0);
}

void
usage(opt_t *opt, char *argv0)
{
	fprintf(stderr, "usage: %s [OPTIONS] file.vcf\n", argv0);
	fprintf(stderr, "   -h MAX_DEPTH   ignore sites with depth higher than MAX_DEPTH [%d]\n", opt->max_depth);
	fprintf(stderr, "   -l MIN_DEPTH   ignore sites with depth lower than MIN_DEPTH [%d]\n", opt->min_depth);
	fprintf(stderr, "   -s STEP        move window along chromosomes in steps of STEP bp [%d]\n", opt->step);
	fprintf(stderr, "   -w WINDOW      output windows of size WINDOW bp [%d]\n", opt->window);
	exit(1);
}

int
main(int argc, char **argv)
{
	int c;
	char *vcf_fn;
	opt_t opt;

	opt.min_depth = 0;
	opt.max_depth = 10000;
	opt.window = 5*1000*1000;
	opt.step = 200*1000;

	while ((c = getopt(argc, argv, "h:l:s:w:")) != -1) {
		switch (c) {
			case 'h':
				opt.max_depth = strtoul(optarg, NULL, 0);
				if (opt.max_depth < 0 || opt.max_depth > 10000) {
					fprintf(stderr, "max_depth=`%s' out of range\n", optarg);
					return -1;
				}
				break;
			case 'l':
				opt.min_depth = strtoul(optarg, NULL, 0);
				if (opt.min_depth < 0 || opt.min_depth > 10000) {
					fprintf(stderr, "min_depth=`%s' out of range\n", optarg);
					return -1;
				}
				break;
			case 's':
				opt.step = parse_bp(optarg);
				if (opt.step < 0 || opt.step > 100*1000*1000) {
					fprintf(stderr, "step=`%s' out of range\n", optarg);
					return -1;
				}
				break;
			case 'w':
				opt.window = parse_bp(optarg);
				if (opt.window < 0 || opt.window > 100*1000*1000) {
					fprintf(stderr, "window=`%s' out of range\n", optarg);
					return -1;
				}
				break;
			default:
				usage(&opt, argv[0]);
		}
	}

	if (argc-optind != 1) {
		usage(&opt, argv[0]);
	}

	if (opt.step > opt.window) {
		fprintf(stderr, "Error: step (%d) > window (%d) doesn't make sense.",
				opt.step, opt.window);
		return -1;
	}

	if (opt.window % opt.step != 0) {
		fprintf(stderr, "Error: step (%d) must divide window (%d) with no remainder.",
				opt.step, opt.window);
		return -1;
	}

	vcf_fn = argv[optind];

	if (het_windows(&opt, vcf_fn) < 0) {
		return -1;
	}

	return 0;
}

