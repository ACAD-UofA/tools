/*
 * Print length of sequences from input file.
 *
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
#include <stdint.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>

#include <htslib/sam.h>

#include "kseq.h"
#include "kmath.h"
#include "distributions.h"

KSEQ_INIT(int, read);

typedef struct {
	// source of fragments
	int simulate;
	int fasta;
	int bam;

	// input fasta/bam file
	char *fn;

	// bam options
	int mapped;
	int ignore_dupes;
	int ignore_pairs;
	char *contig;
} opt_t;

typedef struct {
	uint64_t *hist;
	int l_min, l_max;
	uint64_t area;
} hist_t;

/*
 * Discretised lognormal density function.
 */
double
lnorm(int x, double mu, double s)
{
	return 0.5*(erf(s*(log(x+0.5)-mu)) -erf(s*(log(x-0.5)-mu)));
}

double
ssqe_exp2(int n, double *x, void *data)
{
	hist_t *h = data;
	double a = x[0];
	double b = x[1];
	double ssqe = 0;
	int i = 0;

	for (i=h->l_min; i<=h->l_max; i++) {
		double diff = exp(-a*i-b/i) - (double)h->hist[i]/h->area;
		ssqe += diff*diff;
		if (isnan(ssqe))
			return INFINITY;
	}

	return ssqe;
}
double
ssqe_exp2mix(int n, double *x, void *data)
{
	hist_t *h = data;
	double a = x[0];
	double b = x[1];
	double c = x[2];
	double y = x[3];
	double ssqe = 0;
	int i = 0;

	if (y < 0 || y > 1)
		return INFINITY;

	if (a < 0 || b < 0 || c < 0)
		return INFINITY;

	for (i=h->l_min; i<=h->l_max; i++) {
		double diff = y*exp(-a*i-b/i)+(1-y)*exp(-c*i-b/i) - (double)h->hist[i]/h->area;
		ssqe += diff*diff;
		if (isnan(ssqe))
			return INFINITY;
	}

	return ssqe;
}

double
ssqe_lnorm(int n, double *x, void *data)
{
	hist_t *h = data;
	double mu = x[0];
	double s = M_SQRT1_2/x[1];
	double ssqe = 0;
	int i = 0;

	for (i=h->l_min; i<=h->l_max; i++) {
		double p_i = lnorm(i, mu, s);
		double diff = p_i - (double)h->hist[i]/h->area;
		ssqe += diff*diff; // / (h->hist[i]*p_i*(1-p_i));
		if (isnan(ssqe))
			return INFINITY;
	}

	return ssqe;
}

double
ssqe_lnormmix(int n, double *x, void *data)
{
	hist_t *h = data;
	double mu1 = x[0];
	double s1 = M_SQRT1_2/x[1];
	double mu2 = x[2];
	double s2 = M_SQRT1_2/x[3];
	double y = x[4];
	double ssqe = 0;
	int i = 0;

	if (y < 0 || y > 1)
		return INFINITY;

	if (mu1 < 0 || mu1 > 10 || mu2 < 0 || mu2 > 10)
		return INFINITY;

	if (x[1] < 0 || x[3] < 0)
		return INFINITY;

	for (i=h->l_min; i<=h->l_max; i++) {
		double p_i1 = lnorm(i, mu1, s1);
		double p_i2 = lnorm(i, mu2, s2);
		double diff = y*p_i1 + (1-y)*p_i2 - (double)h->hist[i]/h->area;
		ssqe += diff*diff; // / (h->hist[i]*(p_i1*(1-p_i1) + p_i2*(1-p_i2)));
		if (isnan(ssqe))
			return INFINITY;
	}

	return ssqe;
}

double
ssqe_lnormgeom(int n, double *x, void *data)
{
	hist_t *h = data;
	double mu = x[0];
	double s = M_SQRT1_2/x[1];
	double p = x[2];
	double ssqe = 0;
	int i = 0;

	for (i=h->l_min; i<=h->l_max; i++) {
		double diff = lnorm(i, mu, s)*pow(1-p, i-1)*p - (double)h->hist[i]/h->area;
		ssqe += diff*diff;
		if (isnan(ssqe))
			return INFINITY;
	}

	return ssqe;
}

double
ssqe_weibull(int n, double *x, void *data)
{
	hist_t *h = data;
	double delta = x[0];
	double eta = x[1];
	double d0 = x[2];
	double ssqe = 0;
	int i = 0;

	for (i=h->l_min; i<=h->l_max; i++) {
		double diff = (delta/eta)*pow((i-d0)/eta,delta-1)*exp(-pow((i-d0)/eta,delta)) - (double)h->hist[i]/h->area;
		ssqe += diff*diff;
		if (isnan(ssqe))
			return INFINITY;
	}

	return ssqe;
}

double
ssqe_rayleigh(int n, double *x, void *data)
{
	hist_t *h = data;
	double loc = x[0];
	double s = 1.0/x[1];
	double ssqe = 0;
	int i = 0;

	for (i=h->l_min; i<=h->l_max; i++) {
		double a = ((double)i-loc)*s;
		double diff =  a*s*exp(-0.5*a*a) - (double)h->hist[i]/h->area;
		ssqe += diff*diff;
		if (isnan(ssqe))
			return INFINITY;
	}

	return ssqe;
}


int
hist_simulate(opt_t *opt, hist_t *h)
{
	int i, j;
	uint64_t k;
	int ret;
	extern uint64_t rpois(krand_t *kr, double lam);

	//const int chrlen = 16500; // mito
	const int chrlen = 1e6; // 1mb chrom
	const int copies = 1;
	const int niter = 100*copies;
	const double a=1e-7, b=2.8, c=0;
	const double mean=0.5, std=0.05;

	rk_state state;
	rk_seed(123, &state);

	h->hist = calloc(chrlen+1, sizeof(uint64_t));
	if (h->hist == NULL) {
		perror("calloc");
		ret = -1;
		goto err0;
	}

	h->hist[chrlen] = copies;
	h->area = copies;

	int brk;
	double p;
	uint64_t nbreaks;

	for (i=0; i<niter; i++) {
		for (j=2; j<=chrlen; j++) {
			if (h->hist[j] == 0)
				continue;
			p = a*pow(j, b)+c;
			if (p > 1)
				p = 1;
			nbreaks = rk_binomial(&state, h->hist[j], p);
			for (k=0; k<nbreaks; k++) {
				do {
					brk = round((rk_normal(&state, mean, std) * j));
				} while (brk<1 || brk>j-1);
				h->hist[brk]++;
				h->hist[j-brk]++;
				h->hist[j]--;
				h->area++;
			}
		}
	}

	/*
	for (j=1; j<=chrlen; j++) {
		if (h->hist[j] == 0)
			continue;
		p = alpha*exp(beta*j);
	}*/

	for (j=chrlen; j>0; j--) {
		if (h->hist[j]) {
			h->l_max = j;
			break;
		}
	}

	ret = 0;
err0:
	return ret;
}

int
hist_from_bam(opt_t *opt, hist_t *h)
{
	samFile *fp;
	bam_hdr_t *hdr;
	bam1_t *b;
	int filter;
	int tid = -1;
	int ret;
	int len;
	size_t h_slots, h_slots2;
	uint64_t *tmp;

	h_slots = 0;
	h->hist = NULL;
	h->l_max = 0;
	h->area = 0;

	fp = sam_open(opt->fn, "r");
	if (fp == NULL) {
		ret = -1;
		goto err0;
	}

	hdr = sam_hdr_read(fp);
	if (hdr == NULL) {
		fprintf(stderr, "%s: no SAM/BAM header found.\n", opt->fn);
		ret = -2;
		goto err1;
	}

	b = bam_init1();
	if (b == NULL) {
		ret = -3;
		goto err2;
	}

	filter = BAM_FSECONDARY|BAM_FQCFAIL|BAM_FSUPPLEMENTARY;
	if (opt->ignore_dupes)
		filter |= BAM_FDUP;
	if (opt->ignore_pairs)
		filter |= BAM_FPAIRED;
	if (opt->mapped == 1)
		filter |= BAM_FUNMAP;

	if (opt->contig) {
		tid = bam_name2id(hdr, opt->contig);
		if (tid == -1) {
			fprintf(stderr, "%s: could not find contig `%s'\n",
					opt->fn, opt->contig);
			return -1;
		}
	}

	while (sam_read1(fp, hdr, b) >= 0) {
		bam1_core_t *c = &b->core;

		if (opt->contig && tid != c->tid)
			continue;

		if (c->flag & filter)
			continue;

		if (opt->mapped == 0 && ((c->flag & BAM_FUNMAP) != BAM_FUNMAP))
			continue;

		if (opt->mapped && (c->flag & BAM_FPAIRED)) {
			if (c->flag & BAM_FREAD2)
				continue;

			if ((c->flag & (BAM_FPROPER_PAIR|BAM_FMUNMAP)) != BAM_FPROPER_PAIR)
				continue;

			len = c->isize; // tlen
			if (len < 0)
				len = -len;
			else if (len == 0)
				continue;
		} else {
			len = c->l_qseq;
		}


		if (len > h->l_max) {
			h->l_max = len;
			if (len >= h_slots) {
				h_slots2 = h_slots;
				h_slots = len*2;
				tmp = realloc(h->hist, sizeof(uint64_t)*h_slots);
				if (tmp == NULL) {
					free(h->hist);
					h->hist = NULL;
					h->l_max = 0;
					ret = -4;
					goto err3;
				}
				h->hist = tmp;
				memset(&tmp[h_slots2], 0, (h_slots-h_slots2)*sizeof(uint64_t));
			}
		}

		h->hist[len]++;
		h->area++;
	}

	ret = 0;

err3:
	bam_destroy1(b);
err2:
	bam_hdr_destroy(hdr);
err1:
	sam_close(fp);
err0:
	return ret;
}

int
hist_from_fasta(opt_t *opt, hist_t *h)
{
	kseq_t *seq;
	int fd;
	int ret;
	int len;
	size_t h_slots, h_slots2;
	uint64_t *tmp;

	h_slots = 0;
	h->hist = NULL;
	h->l_max = 0;
	h->area = 0;

	if (!strcmp(opt->fn, "-")) {
		fd = 0;
	} else {
		fd = open(opt->fn, O_RDONLY);
		if (fd == -1) {
			fprintf(stderr, "%s: %s\n", opt->fn, strerror(errno));
			ret = -1;
			goto err0;
		}
	}

	seq = kseq_init(fd);

	while ((len = kseq_read(seq)) >= 0) {

		if (seq->seq.l > h->l_max) {
			h->l_max = seq->seq.l;
			if (h->l_max >= h_slots) {
				h_slots2 = h_slots;
				h_slots = h->l_max*2;
				tmp = realloc(h->hist, sizeof(uint64_t)*h_slots);
				if (tmp == NULL) {
					free(h->hist);
					h->hist = NULL;
					h->l_max = 0;
					ret = -2;
					goto err1;
				}
				h->hist = tmp;
				memset(&tmp[h_slots2], 0, (h_slots-h_slots2)*sizeof(uint64_t));
			}
		}

		h->hist[seq->seq.l]++;
		h->area++;
	}

	if (len != -1) {
		fprintf(stderr, "%s: unexpected end of file.\n", opt->fn);
	}

	ret = 0;
err1:
	kseq_destroy(seq);
	if (fd != 0)
		close(fd);
err0:
	return ret;
}

/*
 * return 1 if @s ends with @match, 0 otherwise.
 */
int
endswith(char *s, char *match)
{
	int l1, l2;
	int i, j;

	l1 = strlen(s);
	l2 = strlen(match);

	if (l1 < l2)
		return 0;

	for (i=l1-1, j=l2-1; j>=0; i--, j--) {
		if (s[i] != match[j])
			return 0;
	}

	return 1;
}

void
usage(char *argv0)
{
	fprintf(stderr, "usage: %s [-s|-f|-baudp [-c chrN]] in.file\n", argv0);
	exit(1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	int c;
	int i;

	// defaults
	opt.fasta = opt.bam = opt.simulate = 0;
	opt.mapped = 1; // include only mapped reads
	opt.ignore_dupes = 0; // accept reads marked as duplicates
	opt.ignore_pairs = 0; // accept paired reads
	opt.contig = NULL;
	opt.fn = NULL;


	while ((c = getopt(argc, argv, "sfbaudpc:")) != -1)
	{
		switch (c) {
			case 's':
				opt.simulate = 1;
				break;
			case 'f':
				opt.fasta = 1;
				break;
			case 'b':
				opt.bam = 1;
				break;
			case 'a':
				// include mapped and unmapped reads
				opt.mapped = 2;
				break;
			case 'u':
				// include only unmapped reads
				opt.mapped = 0;
				break;
			case 'd':
				// ignore reads marked as duplicates
				opt.ignore_dupes = 1;
				break;
			case 'p':
				// ignore paired reads;
				opt.ignore_pairs = 1;
				break;
			case 'c':
				opt.contig = optarg;
				break;
			default:
				usage(argv[0]);
		}
	}

	if (argc-optind != (1-opt.simulate))
		usage(argv[0]);

	opt.fn = argv[optind];

	if ((opt.simulate && opt.fasta) || (opt.simulate && opt.bam) || (opt.fasta && opt.bam)) {
		fprintf(stderr, "-s, -f and -b flags are mutually exclusive.\n");
		return -1;
	}

	if (!opt.fasta && !opt.bam && !opt.simulate) {
		/* try to detect from filename */
		if (endswith(opt.fn, ".fasta") || endswith(opt.fn, ".fastq") ||
		    endswith(opt.fn, ".fa") || endswith(opt.fn, ".fq")) {
			opt.fasta = 1;
		} else if (endswith(opt.fn, ".bam") || endswith(opt.fn, ".sam")) {
			opt.bam = 1;
		}
	}

	if (!opt.fasta && !opt.bam && !opt.simulate) {
		fprintf(stderr, "%s: unknown input format: "
				"must specify -f for fasta/fastq or -b for bam/sam.\n",
				opt.fn);
		return -2;
	}

	hist_t hist;

	if (opt.fasta) {
		if (hist_from_fasta(&opt, &hist) < 0)
			return -1;
	}
	else if (opt.bam) {
		if (hist_from_bam(&opt, &hist) < 0)
			return -1;
	}
	else if (opt.simulate) {
		if (hist_simulate(&opt, &hist) < 0)
			return -1;
	}

	hist.l_min = hist.l_max;
	for (i=1; i<=hist.l_max; i++) {
		if (hist.hist[i]) {
			hist.l_min = i;
			break;
		}
	}

	uint64_t l_max = hist.l_max;
	//if (hist.l_max > 300)
	//	hist.l_max = 250;

	double ssqe;

	{
	double x[] = {5.0, 0.5}; // initial values for: meanlog, sdlog
	ssqe = kmin_hj(ssqe_lnorm, 2, x, &hist, KMIN_RADIUS, KMIN_EPS, KMIN_MAXCALL);
	printf("#lnorm;SSqE=%lg;meanlog=%lg;sdlog=%lg\n", ssqe, x[0], x[1]);
	}

	{
	double x[] = {0.05, 500}; // a, b
	ssqe = kmin_hj(ssqe_exp2, 2, x, &hist, KMIN_RADIUS, KMIN_EPS, KMIN_MAXCALL);
	printf("#exp2;SSqE=%lg;a=%lg;b=%lg\n", ssqe, x[0], x[1]);
	}

	{
	double x[] = {0.05, 500, 0.05, 0.5}; // a, b, c, y
	ssqe = kmin_hj(ssqe_exp2mix, 4, x, &hist, KMIN_RADIUS, KMIN_EPS, KMIN_MAXCALL);
	printf("#exp2mix;SSqE=%lg;a=%lg;b=%lg;c=%lg;y=%lg\n", ssqe, x[0], x[1], x[2], x[3]);
	}

	{
	double x[] = {5.0, 0.5, 5.0, 1.5, 0.9}; // meanlog1, sdlog1, meanlog2, sdlog2, gamma
	ssqe = kmin_hj(ssqe_lnormmix, 5, x, &hist, KMIN_RADIUS, KMIN_EPS, KMIN_MAXCALL);
	printf("#lnormmix;SSqE=%lg;meanlog1=%lg;sdlog1=%lg;meanlog2=%lg;sdlog2=%lg;gamma=%lg\n",
			ssqe, x[0], x[1], x[2], x[3], x[4]);
	}

	/*
	{
	double ssqe;
	double x[] = {5.0, 0.5, 0.001}; // initial values for: meanlog, sdlog, p
	ssqe = kmin_hj(ssqe_lnormgeom, 3, x, &hist, KMIN_RADIUS, KMIN_EPS, KMIN_MAXCALL);
	printf("#lnormgeom;SSqE=%lg;meanlog=%lg;sdlog=%lg;p=%lg\n", ssqe, x[0], x[1], x[2]);
	}
	*/

	{
	double x[] = {1.5, 30, 20}; // initial values for: delta, eta, d0
	ssqe = kmin_hj(ssqe_weibull, 3, x, &hist, KMIN_RADIUS, KMIN_EPS, KMIN_MAXCALL);
	printf("#weibull;SSqE=%lg;delta=%lg;eta=%lg;d0=%lg\n", ssqe, x[0], x[1], x[2]);
	}

	{
	double x[] = {0.0, 50.0}; // initial values for: loc, scale
	ssqe = kmin_hj(ssqe_rayleigh, 2, x, &hist, KMIN_RADIUS, KMIN_EPS, KMIN_MAXCALL);
	printf("#rayleigh;SSqE=%lg;loc=%lg;scale=%lg\n", ssqe, x[0], x[1]);
	}

	for (i=1; i<=l_max; i++)
		printf("%jd\n", (uintmax_t)hist.hist[i]);

	free(hist.hist);

	return 0;
}
