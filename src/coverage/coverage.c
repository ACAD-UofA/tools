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

#include <htslib/vcf.h>

typedef struct {
	int win_sz;
} opt_t;

typedef struct {
	uint64_t *hist;
	size_t slots;

	uint64_t *window;
	uint64_t wcount;

	int min, max;
	uint64_t area;
} hist_t;

int
dump_cov(opt_t *opt, bcf_hdr_t *hdr, hist_t *h, int n)
{
	int i, j;
	FILE *fp;

	for (i=0; i<n; i++) {
		int len = strlen(bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i)) + 64;
		char fn[len+1];
		sprintf(fn, "%s.%d.hist", bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i), opt->win_sz);
		fp = fopen(fn, "w");
		if (fp == NULL) {
			fprintf(stderr, "%s: %s\n", fn, strerror(errno));
			return -1;
		}
		for (j=0; j<h[i].max; j++)
			fprintf(fp, "%zd\n", (intmax_t)h[i].hist[j]);
		fclose(fp);
	}

	return 0;
}

int
coverage(opt_t *opt, char *vcf_fn)
{
	htsFile *fp;
	bcf_hdr_t *hdr;
	bcf1_t *rec;
	int nsamples;
	int ret;

	hist_t *h;
	uint64_t *tmp;
	int i;
	size_t hi; // histogram index

	int n_dp = 0;
	int *dp = NULL;
	int rid = -1;
	size_t wi = 0; // window index

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

	nsamples = bcf_hdr_nsamples(hdr);

	h = calloc(nsamples, sizeof(hist_t));
	if (h == NULL) {
		perror("calloc");
		ret = -4;
		goto err2;
	}

	for (i=0; i<nsamples; i++) {
		h[i].window = calloc(opt->win_sz, sizeof(uint64_t));
		if (h[i].window == NULL) {
			perror("calloc");
			ret = -3;
			goto err3;
		}
	}

	rec = bcf_init1();
	if (rec == NULL) {
		perror("bcf_init1: calloc");
		ret = -3;
		goto err3;
	}

	while (bcf_read1(fp, hdr, rec) >= 0) {
		if (bcf_get_format_int32(hdr, rec, "DP", &dp, &n_dp) < 0) {
			fprintf(stderr, "%s: missing FORMAT/DP field at %s:%d\n",
					vcf_fn, bcf_seqname(hdr, rec), rec->pos+1);
			ret = -5;
			goto err4;
		}

		if (rid != rec->rid) {
			// new chromosome
			for (i=0; i<nsamples; i++) {
				if (rid != -1) {
					// count positions within win_sz of the contig end
					hi = h[i].wcount;
					for (;;) {
						if (++wi == opt->win_sz)
							wi = 0;
						hi -= h[i].window[wi];
						if (hi <= 0)
							break;
						h[i].hist[hi]++;
						h[i].area++;
					}
				}
				// reset window
				h[i].wcount = 0;
				memset(h[i].window, 0, opt->win_sz*sizeof(uint64_t));
			}
			rid = rec->rid;
			wi = 0;
		}

		if (++wi == opt->win_sz)
			wi = 0;

		for (i=0; i<nsamples; i++) {
			int depth = dp[i];
			uint64_t wlast = h[i].window[wi];

			if (depth < 0)
				depth = 0;

			h[i].window[wi] = depth;
			h[i].wcount += depth;
			h[i].wcount -= wlast;
			hi = h[i].wcount;

			if (hi == 0)
				continue;

			if (hi > h[i].max) {
				h[i].max = hi;
				if (hi >= h[i].slots) {
					size_t oldslots = h[i].slots;
					h[i].slots = hi*2;
					tmp = realloc(h[i].hist, sizeof(uint64_t)*h[i].slots);
					if (tmp == NULL) {
						perror("realloc");
						ret = -6;
						goto err4;
					}
					memset(tmp+oldslots, 0, (h[i].slots-oldslots)*sizeof(uint64_t));
					h[i].hist = tmp;
				}
			}

			h[i].hist[hi]++;
			h[i].area++;
		}
	}

	if (rid != -1) {
		for (i=0; i<nsamples; i++) {
			// count positions within win_sz of the contig end
			hi = h[i].wcount;
			for (;;) {
				if (++wi == opt->win_sz)
					wi = 0;
				hi -= h[i].window[wi];
				if (hi <= 0)
					break;
				h[i].hist[hi]++;
				h[i].area++;
			}
		}
	}


	ret = 0;
	if (dump_cov(opt, hdr, h, nsamples) < 0)
		ret = -7;

err4:
	if (dp)
		free(dp);

	for (i=0; i<nsamples; i++) {
		if (h[i].hist)
			free(h[i].hist);
	}
	bcf_destroy1(rec);
err3:
	for (i=0; i<nsamples; i++) {
		if (h[i].window)
			free(h[i].window);
	}
	free(h);
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
	fprintf(stderr, "usage: %s [-w WIN] file.{b,v}cf[.gz]\n", argv0);
	exit(1);
}

int
main(int argc, char **argv)
{
	int c;
	char *vcf_fn;
	opt_t opt;

	opt.win_sz = 1;

	while ((c = getopt(argc, argv, "w:")) != -1) {
		switch (c) {
			case 'w':
				opt.win_sz = strtol(optarg, NULL, 0);
				if (opt.win_sz < 1 || opt.win_sz > 10000) {
					fprintf(stderr, "windowsize=`%s' out of bounds", optarg);
					usage(argv[0]);
				}
				break;
			default:
				usage(argv[0]);
		}
	}

	if (argc-optind != 1) {
		usage(argv[0]);
	}

	vcf_fn = argv[optind];

	if (coverage(&opt, vcf_fn) < 0) {
		return -1;
	}

	return 0;
}
