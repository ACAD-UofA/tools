/*
 * Test if groups belong to the same population, based upon allele
 * frequencies determined by genotyping.
 *
 * Copyright (c) 2015 Graham Gower <graham.gower@gmail.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdint.h>
#include <math.h>
#include <errno.h>
#include <alloca.h>

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <bcftools.h>
#include <version.h> // must be bcftools/version.h

//#include "binom.h"

#define alen(a) (sizeof(a)/sizeof(a[0]))

#define N	1000 // number of values calculated for each distribution
#define NINV	(1.0/((double)N))

/*
 * Determine AF at maximum a posteriori values, 95% credibility intervals,
 * and proportion of overlap between the two distributions.
 */
void
test_af_dists(double dist1[], double dist2[], double dsum1, double dsum2,
		double *p,
		double *af1_map, double *ci1_low, double *ci1_high,
		double *af2_map, double *ci2_low, double *ci2_high)
{
	double diff12; // pairwise distance between a posteriori values
	double p1, p2; // scaled probabilties
	double c1, c2; // cumulative values
	double m1, m2; // maximum values
	int i1, i2; // index at maximum values
	int i;

	*ci1_low = *ci1_high = *ci2_low = *ci2_high = 0;

	diff12 = 0.0;
	c1 = c2 = 0.0;
	m1 = m2 = 0.0;
	i1 = i2 = 0;

	for (i=0; i<N; i++) {

		// scale values such that total area under the curve is 1
		p1 = dist1[i] / dsum1;
		p2 = dist2[i] / dsum2;

		diff12 += fabs(p1 - p2);

		if (p1 > m1) {
			m1 = p1;
			i1 = i;
		}
		if (p2 > m2) {
			m2 = p2;
			i2 = i;
		}

		c1 += p1;
		c2 += p2;
		if (!*ci1_low && c1 >= 0.025 )
			*ci1_low = i*NINV;
		if (!*ci1_high && c1 >= 0.975)
			*ci1_high = i*NINV;
		if (!*ci2_low && c2 >= 0.025 )
			*ci2_low = i*NINV;
		if (!*ci2_high && c2 >= 0.975)
			*ci2_high = i*NINV;
	}

	// proportion of area overlap between the two distributions.
	*p = 1.0 - (diff12/2);

	*af1_map = i1*NINV;
	*af2_map = i2*NINV;
}

/*
 * AF probability distributions of the populations.
 */
void
af_dist(int n_pops, unsigned int *ac, unsigned int *n,
		double **dist, double *dsum)
{
	int i, j, k;
	int sum_ac, sum_n;
	double x, p;
	double p_x_given_ac, prior_p_x;
	double *prior_af;
	
	prior_af = alloca(n_pops*sizeof(double));

	for (j=0; j<n_pops; j++) {
		dsum[j] = 0.0; // areas under the distributions

		/*
		 * Prior expectation of AF from other populations.  Use the
		 * approximate mean, biased towards 0.5 for small sample sizes.
		 * This is our hyperparameter for calculation of the prior.
		 */
		sum_ac = sum_n = 0;
		for (k=0; k<n_pops; k++) {
			if (k != j) {
				sum_ac += ac[k];
				sum_n += n[k];
			}
		}
		prior_af[j] = ((float)(sum_ac + 1)) / (sum_n + 2);
	}

	for (i=0; i<N; i++) {
		x = i*NINV;
		for (j=0; j<n_pops; j++) {

			// weak prior, biased towards prior expectation of AF
//			prior_p_x = binom_lpmf(prior_af[j], 1, x);
			prior_p_x = log(x)*prior_af[j] + log1p(-x)*(1-prior_af[j]);

//			p_x_given_ac = binom_lpmf(ac[j], n[j], x);
			p_x_given_ac = log(x)*ac[j] + log1p(-x)*(n[j]-ac[j]);

			p = exp(p_x_given_ac + prior_p_x);
			if (isnan(p))
				p = 0.0;

			dist[j][i] = p;
			dsum[j] += p;
		}
	}
}




void
version(const char **bcftools_ver, const char **htslib_ver)
{
    *bcftools_ver = BCFTOOLS_VERSION;
    *htslib_ver = hts_version();
}

/*
 * Short description used by 'bcftools plugin -l'.
 */
const char *
about(void)
{
    return "Perform a population test on each site.\n";
}

/*
 * longer description used by 'bcftools +name -h'
 */
const char *
usage(void)
{
    return 
        "\n"
        "About: Perform a two population test on each site.\n"
        "Usage: bcftools +poptest [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -p <file>   Samples in a population. One sample per line."
	"               May be specified multiple times.\n"
        "\n"
        "Example:\n"
        "   bcftools +poptest in.vcf -- -p pop1.list -p pop2.list\n"
        "\n";
}

static bcf_hdr_t *hdr_in;
static unsigned int n_pops = 0; // number of populations
static int *samples_per_pop; // number of samples in each population
static unsigned int **pop_indices; // indices of samples in each population
static double **dist; // AF probability distributions for each population
static double *fst_vec; // fst values from the jacknife

/*
 * Called once at startup, allows to initialize local variables.
 * Return 1 to suppress VCF/BCF header from printing, -1 on critical
 * errors, 0 otherwise.
 */
int
init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
	char **pop_filename = NULL;
	char ***samples;
	int c;
	int i, j;
	int total_samples = 0;

	while ((c = getopt(argc, argv, "p:")) != -1)
	{
		switch (c) {
			case 'p':
				pop_filename = realloc(pop_filename, (n_pops+1)*sizeof(char *));
				if (pop_filename == NULL) {
					perror("realloc");
					return -1;
				}
				pop_filename[n_pops] = optarg;
				n_pops++;
				continue;
			default:
				fprintf(stderr,"%s", usage());
				return -1;
		}
	}

	if (n_pops == 0) {
		fprintf(stderr, "Must specify files containing the samples "
			       "in the populations\nwhich are to be tested.\n");
		fprintf(stderr,"%s", usage());
		return -1;
	}

	samples = malloc(n_pops*sizeof(char **));
	if (samples == NULL) {
		perror("malloc");
		return -1;
	}

	samples_per_pop = malloc(n_pops*sizeof(unsigned int));
	if (samples_per_pop == NULL) {
		perror("malloc");
		return -1;
	}

	size_t tmp_popstr_len = 0;
	char *tmp_popstr;

	for (i=0; i<n_pops; i++) {
		samples[i] = hts_readlist(pop_filename[i], 1, &samples_per_pop[i]);
		if (samples[i] == NULL) {
			fprintf(stderr, "Error reading %s: %s\n",
					pop_filename[i], strerror(errno));
			return -1;
		}

		for (j=0; j<samples_per_pop[i]; j++) 
			tmp_popstr_len += strlen(samples[i][j]);
		tmp_popstr_len += samples_per_pop[i];
		total_samples += samples_per_pop[i];
	}
	tmp_popstr_len++;

	/*
	 * Construct a temporary comma separated string containing samples
	 * from all populations, for passing into bcf_hdr_set_samples().
	 * This is a bit contrived, it would be nice to be able to
	 * pass in each pop_filename to add more samples to the set with
	 * each function call.
	 * Htslib does not appear to work that way however.
	 */
	tmp_popstr = malloc(tmp_popstr_len*sizeof(char));
	if (tmp_popstr == NULL) {
		perror("malloc");
		return -1;
	}
	tmp_popstr[0] = '\0';

	for (i=0; i<n_pops; i++) {
		for (j=0; j<samples_per_pop[i]; j++) {
			strcat(tmp_popstr, samples[i][j]);
			strcat(tmp_popstr, ",");
		}
	}

	// remove trailing comma
	tmp_popstr[tmp_popstr_len-2] = '\0';

	int idx;
	idx = bcf_hdr_set_samples(in, tmp_popstr, 0);
	if (idx != 0) {
		if (idx > 0) {
			idx--;
			/* 
			 * Find the name of the sample that is not in the vcf
			 * and the file in which it was specified.
			 */
			for (i=0; i<n_pops; i++) {
				if (idx < samples_per_pop[i])
					break;
				idx -= samples_per_pop[i];
			}
			fprintf(stderr, "Error: %s:%s not in vcf\n",
					pop_filename[i], samples[i][idx]);
		} else
			perror("malloc");
		return -1;
	}

	free(tmp_popstr);
	free(pop_filename);

	/*
	 * Construct a list of indices for each population.
	 */
	pop_indices = malloc(n_pops*sizeof(unsigned int *));
	if (pop_indices == NULL) {
		perror("malloc");
		return -1;
	}

	for (i=0; i<n_pops; i++) {
		pop_indices[i] = malloc(samples_per_pop[i]*sizeof(unsigned int));
		if (pop_indices[i] == NULL) {
			perror("malloc");
			return -1;
		}
		for (j=0; j<samples_per_pop[i]; j++) {
			pop_indices[i][j] = 2*bcf_hdr_id2int(in, BCF_DT_SAMPLE, samples[i][j]);
			free(samples[i][j]);
		}
		free(samples[i]);
	}

	free(samples);

	dist = malloc(n_pops*sizeof(double*));
	if (dist == NULL) {
		perror("malloc");
		return -1;
	}
	for (i=0; i<n_pops; i++) {
		dist[i] = malloc(N*sizeof(double));
		if (dist[i] == NULL) {
			perror("malloc");
			return -1;
		}
	}

	// preallocate space for the fst jacknife
	fst_vec = malloc(total_samples * sizeof(double));
	if (fst_vec == NULL) {
		perror("malloc");
		return -1;
	}

	hdr_in = in;

	return 1;
}

/*
 * Fst, from Weir & Cockerham (1984).
 */
double
fst(unsigned int n_pops, unsigned int *ac, unsigned int *n, unsigned int *h)
{
	int i;
	unsigned int r = n_pops;
	unsigned int n_sum, n2_sum, p_sum, h_sum;
	double nbar, nc, pbar, s2, hbar;
	double a, b, c;

	assert(n_pops >= 2);

	n_sum = 0;
	n2_sum = 0;
	p_sum = 0;
	h_sum = 0;

	for (i=0; i<r; i++) {
		n_sum += n[i];
		n2_sum += n[i]*n[i];
		p_sum += ac[i]; //n[i] * p[i];
		h_sum += n[i] * h[i];
	}

	nbar = ((float)n_sum) / r;
	nc = (n_sum - ((float)n2_sum)/(r*nbar)) / (r - 1.0);
	pbar = ((float)p_sum) / n_sum;
	hbar = ((float)h_sum) / n_sum;

	s2 = 0.0;
	for (i=0; i<r; i++)
		//s2 += n[i] * (p[i] - pbar);
		s2 += n[i] * ( ((float)ac[i])/n[i] - pbar);
	s2 /= (r - 1) * nbar;

	// Equation (2)
	a = (nbar/nc) * (s2 - (1.0/(nbar-1.0))*(pbar*(1.0-pbar) - s2*(r-1.0)/r - hbar/4.0));
	// Equation (3)
	b = (nbar/(nbar -1.0)) * (pbar*(1.0-pbar) - s2*(r-1.0)/r - hbar*(2.0*nbar -1.0)/(4.0*nbar));
	// Equation (4)
	c = hbar/2.0;

	// return Theta`hat, from Equation (1)
	return a / (a + b + c);
}

void
fst_jackknife(int n_pops, int *gt_arr,
		double *fst_mean, double *ci_low, double *ci_high)
{
	int i, j, k=0;
	int total_obs; // number of samples, across all populations
	int idx;
	unsigned int gt0, gt1; // REF/ALT gentoype values
	unsigned int *ac; // (minor) allele count
	unsigned int *n; // total allele count
	unsigned int *h; // heterozygote count

	ac = alloca(n_pops*sizeof(unsigned int));
	n = alloca(n_pops*sizeof(unsigned int));
	h = alloca(n_pops*sizeof(unsigned int));

	total_obs = 0;
	for (i=0; i<n_pops; i++) {
		ac[i] = 0;
		n[i] = 0;
		h[i] = 0;
		for (j=0; j<samples_per_pop[i]; j++) {
			idx = pop_indices[i][j];
			if (bcf_gt_is_missing(gt_arr[idx]))
				continue;
			gt0 = bcf_gt_allele(gt_arr[idx]);
			gt1 = bcf_gt_allele(gt_arr[idx+1]);
			ac[i] += gt0 + gt1;
			n[i] += 2;
			h[i] += (gt0 == gt1);

			total_obs++;
		}
	}

	*fst_mean = fst(n_pops, ac, n, h);

	for (i=0; i<n_pops; i++) {
		for (j=0; j<samples_per_pop[i]; j++) {
			idx = pop_indices[i][j];
			if (bcf_gt_is_missing(gt_arr[idx]))
				continue;
			gt0 = bcf_gt_allele(gt_arr[idx]);
			gt1 = bcf_gt_allele(gt_arr[idx+1]);

			// subtract this sample from the data set
			ac[i] -= gt0 + gt1;
			n[i] -= 2;
			h[i] -= (gt0 == gt1);

			// find the Fst without the sample
			fst_vec[k++] = fst(n_pops, ac, n, h);

			// add it back in
			ac[i] += gt0 + gt1;
			n[i] += 2;
			h[i] += (gt0 == gt1);
		}
	}

	double sum, mean, var, sd;

	sum = 0.0;
	for (i=0; i< total_obs; i++)
		sum += fst_vec[i];
	mean = sum / total_obs;
	//*fst_mean = mean;

	var = 0.0;
	for (i=0; i< total_obs; i++) {
		double dev = fst_vec[i] - mean;
		var += dev*dev;
	}
	var *= ((float)total_obs -1.0) / total_obs;
	sd = sqrt(var);

	// 95% confidence intervals (or thereabouts, SD*1.97 is good
	// enough, I'm too lazy to figure out a proper t-distribution)
	*ci_low = mean - sd*1.97;
	*ci_high = mean + sd*1.97;
}

/* Genotype array. This is realloc()'d in bcf_get_genotypes(). */
static int *gt_arr = NULL;

/*
 * Called for each VCF record.
 * Return rec to output the line or NULL to suppress output.
 */
bcf1_t *
process(bcf1_t *rec)
{
	static uint64_t n_records = 0;
	int ngt, ngt_arr = 0;
	int i, j, idx;
	unsigned int *ac; // (minor) allele count
	unsigned int *n; // total allele count
	double *dsum;
	double fst_mean, fst_ci_low, fst_ci_high;

	ac = alloca(n_pops*sizeof(unsigned int));
	n = alloca(n_pops*sizeof(unsigned int));
	dsum = alloca(n_pops*sizeof(double));

	n_records++;

	// we only want biallelic SNPs
	if (!bcf_is_snp(rec) || rec->n_allele != 2)
		return NULL;

	ngt = bcf_get_genotypes(hdr_in, rec, &gt_arr, &ngt_arr);
	if (ngt != ngt_arr) {
		fprintf(stderr, "Failed to get genotypes for record %zd.\n",
				(uintmax_t)n_records);
		exit(1);
	}

	for (i=0; i<n_pops; i++) {
		ac[i] = 0;
		n[i] = 0;
		for (j=0; j<samples_per_pop[i]; j++) {
			idx = pop_indices[i][j];
			if (bcf_gt_is_missing(gt_arr[idx]))
				continue;
			ac[i] += bcf_gt_allele(gt_arr[idx]);
			ac[i] += bcf_gt_allele(gt_arr[idx+1]);
			n[i] += 2;
		}
	}

	fst_jackknife(n_pops, gt_arr, &fst_mean, &fst_ci_low, &fst_ci_high);

	af_dist(n_pops, ac, n, dist, dsum);


	double p;
	double af1_map, ci1_low, ci1_high;
	double af2_map, ci2_low, ci2_high;

	test_af_dists(dist[0], dist[1], dsum[0], dsum[1], &p,
			&af1_map, &ci1_low, &ci1_high,
			&af2_map, &ci2_low, &ci2_high);

	printf("%s\t%d\t%s\t%s\t%.3lf (%d/%d)\t%.3lf (%d/%d)\t%lf\t%lf\t%lf\t%lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
			bcf_seqname(hdr_in, rec), rec->pos+1,
			rec->d.allele[0], rec->d.allele[1],
			((float)ac[0])/n[0], ac[0], n[0], ((float)ac[1])/n[1], ac[1], n[1],
			fst_mean, fst_ci_low, fst_ci_high,
			p, // probability that the two allele frequencies overlap
			af1_map, ci1_low, ci1_high,
			af2_map, ci2_low, ci2_high);

	return NULL;
}

/*
 * Called after all lines have been processed to clean up.
 */
void
destroy(void)
{
	int i;

	free(gt_arr);
	free(samples_per_pop);

	for (i=0; i<n_pops; i++) {
		free(dist[i]);
		free(pop_indices[i]);
	}

	free(pop_indices);
	free(dist);
	free(fst_vec);
}
