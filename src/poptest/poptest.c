/*
 * Test if groups belong to the same population, based upon allele
 * frequencies determined by genotyping.
 *
 * Copyright (c) 2015 Graham Gower <graham.gower@gmail.com>
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
#include <stdint.h>
#include <math.h>
#include <errno.h>
#include <alloca.h>

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <bcftools.h>
#include <version.h> // must be bcftools/version.h

#include "binom.h"

#define alen(a) (sizeof(a)/sizeof(a[0]))
#define max(a,b) ((a)>(b)?(a):(b))


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
	int n_jackknives = 1;

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
		n_jackknives *= samples_per_pop[i];
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

	// preallocate space for the fst jacknife
	fst_vec = malloc(n_jackknives * sizeof(double));
	if (fst_vec == NULL) {
		perror("malloc");
		return -1;
	}

	hdr_in = in;

	printf("CHROM\tPOS\tID\tREF\tALT\tLOG LIK");
	for (i=0; i<n_pops; i++) {
		printf("\t[P%d]MAF (MAC/TAC,HET)", i+1);
	}
	printf("\tFst\tFst_CI_low\tFst_CI_high");
	printf("\n");

	return 1;
}

/*
 * Fst, from Weir & Cockerham (1984).
 * For a biallelic locus.
 */
double
fst(unsigned int n_pops, const unsigned int *ac, const unsigned int *n, const unsigned int *h)
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
		h_sum += h[i];
	}

	if (p_sum == 0 || p_sum == n_sum)
		// one allele is fixed in all populations, Fst is undefined
		//return NAN;
		return 0.0;

	pbar = ((double)p_sum) / n_sum;

	// correct for n[i]=number of alleles, not number of individuals
	n_sum /= 2;
	n2_sum /= 4;

	nbar = ((double)n_sum) / r;
	nc = (n_sum - ((double)n2_sum)/n_sum) / (r - 1.0);
	hbar = ((double)h_sum) / n_sum;

	s2 = 0.0;
	for (i=0; i<r; i++) {
		//s2 += n[i] * (p[i] - pbar)*(p[i] - pbar);
		double p_diff = ((double)ac[i])/n[i] - pbar;
		s2 += n[i] * p_diff*p_diff;
	}
	s2 /= (r - 1) * nbar;

	// correct for n[i]=number of alleles, not number of individuals
	s2 /= 2;

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
fst_jackknife(int n_pops, const int *gt_arr,
		double *fst_mean, double *ci_low, double *ci_high)
{
	int i, j, k, l;
	int total_obs; // number of samples, across all populations
	int idx;
	unsigned int gt0, gt1; // REF/ALT gentoype values
	unsigned int *ac; // (minor) allele count
	unsigned int *n; // total allele count
	unsigned int *h; // heterozygote count

	ac = alloca(n_pops*sizeof(unsigned int));
	n = alloca(n_pops*sizeof(unsigned int));
	h = alloca(n_pops*sizeof(unsigned int));

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
			h[i] += (gt0 != gt1);
		}
	}

	//*fst_mean = fst(n_pops, ac, n, h);

	total_obs = 0;
	for (i=0; i<n_pops-1; i++) {
		for (j=0; j<samples_per_pop[i]; j++) {
			idx = pop_indices[i][j];
			if (bcf_gt_is_missing(gt_arr[idx]))
				continue;
			gt0 = bcf_gt_allele(gt_arr[idx]);
			gt1 = bcf_gt_allele(gt_arr[idx+1]);

			// subtract this sample from the data set
			ac[i] -= gt0 + gt1;
			n[i] -= 2;
			h[i] -= (gt0 != gt1);

			for (k=0; k<n_pops; k++) {
				if (i == k)
					continue;
				for (l=0; l<samples_per_pop[k]; l++) {
					idx = pop_indices[k][l];
					if (bcf_gt_is_missing(gt_arr[idx]))
						continue;
					gt0 = bcf_gt_allele(gt_arr[idx]);
					gt1 = bcf_gt_allele(gt_arr[idx+1]);

					// subtract this sample from the data set
					ac[k] -= gt0 + gt1;
					n[k] -= 2;
					h[k] -= (gt0 != gt1);

					// find the Fst without the sample
					fst_vec[total_obs++] = fst(n_pops, ac, n, h);

					// add it back in
					ac[k] += gt0 + gt1;
					n[k] += 2;
					h[k] += (gt0 != gt1);
				}
			}

			idx = pop_indices[i][j];
			gt0 = bcf_gt_allele(gt_arr[idx]);
			gt1 = bcf_gt_allele(gt_arr[idx+1]);

			// add it back in
			ac[i] += gt0 + gt1;
			n[i] += 2;
			h[i] += (gt0 != gt1);
		}
	}

	double sum, mean, var, sd;

	sum = 0.0;
	for (i=0; i<total_obs; i++)
		sum += fst_vec[i];
	mean = sum / total_obs;
	*fst_mean = mean;

	var = 0.0;
	for (i=0; i< total_obs; i++) {
		double dev = fst_vec[i] - mean;
		var += dev*dev;
	}
	var *= ((double)total_obs -1.0) / total_obs;
	sd = sqrt(var);

	// 95% confidence intervals (or thereabouts, SD*1.97 is good
	// enough, I'm too lazy to figure out a proper distribution)
	*ci_low = mean - sd*1.97;
	*ci_high = mean + sd*1.97;
}

void
fst_jackknife2(int n_pops, const int *gt_arr,
		double *fst_mean, double *ci_low, double *ci_high)
{
	int i, j;
	int total_obs; // number of samples, across all populations
	int idx;
	unsigned int gt0, gt1; // REF/ALT gentoype values
	unsigned int *ac; // (minor) allele count
	unsigned int *n; // total allele count
	unsigned int *h; // heterozygote count

	ac = alloca(n_pops*sizeof(unsigned int));
	n = alloca(n_pops*sizeof(unsigned int));
	h = alloca(n_pops*sizeof(unsigned int));

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
			h[i] += (gt0 != gt1);
		}
	}

	//*fst_mean = fst(n_pops, ac, n, h);

	total_obs = 0;
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
			h[i] -= (gt0 != gt1);

			// find the Fst without the sample
			fst_vec[total_obs++] = fst(n_pops, ac, n, h);

			// add it back in
			ac[i] += gt0 + gt1;
			n[i] += 2;
			h[i] += (gt0 != gt1);
		}
	}

	double sum, mean, var, sd;

	sum = 0.0;
	for (i=0; i<total_obs; i++)
		sum += fst_vec[i];
	mean = sum / total_obs;
	*fst_mean = mean;

	var = 0.0;
	for (i=0; i< total_obs; i++) {
		double dev = fst_vec[i] - mean;
		var += dev*dev;
	}
	var *= ((double)total_obs -1.0) / total_obs;
	sd = sqrt(var);

	// 95% confidence intervals (or thereabouts, SD*1.97 is good
	// enough, I'm too lazy to figure out a proper distribution)
	*ci_low = mean - sd*1.97;
	*ci_high = mean + sd*1.97;
}

/*
 * Joint likelihood of observing the allele counts for all the groups,
 * given the allele frequency is specified by the total allele counts.
 * Normalised by the likelihood of the expected allele count.
 */
double
multi_pop_test(unsigned int n_pops, unsigned int *ac, unsigned int *n)
{
	int i;
	int ac_sum, n_sum;
	double af; // total allele frequency
	double ll; // log likelihood

	ac_sum = n_sum = 0;
	for (i=0; i<n_pops; i++) {
		ac_sum += ac[i];
		n_sum += n[i];
	}

	af = ((double)ac_sum) / n_sum;

	ll = 0.0;
	for (i=0; i<n_pops; i++) {
		if (n[i]) {
			/*
			 *  p(E[ac] | af)
			 *  E[ac] is ~ n[i]*af, but may be the integer below or above.
			 *  Take the one with the higher probability.
			 */
			double lp_e_ac = max(binom_lpmf(floor(n[i]*af), n[i], af), binom_lpmf(ceil(n[i]*af), n[i], af));

			ll += binom_lpmf(ac[i], n[i], af) // p(ac | af)
			       	- lp_e_ac; // p(E[ac] | af)
		}
	}

	return ll;
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
	unsigned int gt0, gt1; // REF/ALT gentoype values
	unsigned int *ac; // (minor) allele count
	unsigned int *n; // total allele count
	unsigned int *h; // heterozygote count
	double loglik;
	double fst_mean, fst_ci_low, fst_ci_high;

	ac = alloca(n_pops*sizeof(unsigned int));
	n = alloca(n_pops*sizeof(unsigned int));
	h = alloca(n_pops*sizeof(unsigned int));

	n_records++;

	// we only want biallelic SNPs
	if (!bcf_is_snp(rec) || rec->n_allele != 2)
		return NULL;

	ngt = bcf_get_genotypes(hdr_in, rec, &gt_arr, &ngt_arr);
	if (ngt != ngt_arr) {
		fprintf(stderr, "Failed to get genotypes for record %ju.\n",
				(uintmax_t)n_records);
		exit(1);
	}

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
			h[i] += (gt0 != gt1);
		}
	}

	loglik = multi_pop_test(n_pops, ac, n);

	printf("%s\t%d\t%s\t%s\t%s",
			bcf_seqname(hdr_in, rec), rec->pos+1, rec->d.id,
			rec->d.allele[0], rec->d.allele[1]);
	printf("\t%lf", loglik);
	for (i=0; i<n_pops; i++)
		printf("\t%.3lf (%d/%d,%d)", ((double)ac[i])/n[i], ac[i], n[i], h[i]);

	fst_jackknife(n_pops, gt_arr, &fst_mean, &fst_ci_low, &fst_ci_high);

	printf("\t%lf\t%lf\t%lf", fst_mean, fst_ci_low, fst_ci_high);
	printf("\n");
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
		free(pop_indices[i]);
	}

	free(pop_indices);
	free(fst_vec);
}
