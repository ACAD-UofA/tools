/*
 * Test if population 1 drifted to become population 2, based upon allele
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

#include <pthread.h>
#include <sched.h>

#include "klist.h"


extern double drift_test(unsigned int mac1, unsigned int tac1,
		unsigned int mac2, unsigned int tac2,
		unsigned int t, unsigned int N);

#define _free_no_op(x)
KLIST_INIT(bcf_rec, bcf1_t *, _free_no_op);
typedef struct {
	pthread_mutex_t mutex;
	klist_t(bcf_rec) *q;
	size_t max_size;
	unsigned int done;
} work_queue_t;

// global data
struct {
	bcf_hdr_t *hdr_in;

	unsigned int n_pops; // number of populations
	int *samples_per_pop; // number of samples in each population
	unsigned int **pop_indices; // indices of samples in each population
	unsigned int n_generations;
	unsigned int population_size;

	unsigned int n_threads;
	pthread_t *threads;
	work_queue_t wq;

	pthread_mutex_t print_mutex;
} gd;

void *drift_worker(void *arg);

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
    return "Test allele frequency drift per site.\n";
}

/*
 * longer description used by 'bcftools +name -h'
 */
const char *
usage(void)
{
    return 
        "\n"
        "About: Test allele frequency drift per site.\n"
        "Usage: bcftools +af_drift_test [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -p <file>   Samples in a population. One sample per line."
	"               May be specified multiple times.\n"
	"   -n <num>    Number of threads.\n"
	"   -t <num>    Number of generations.\n"
	"   -N <num>    Population size.\n"
        "\n"
        "Example:\n"
        "   bcftools +af_drift_test in.vcf -- -n 4 -t 50 -N 1000 -p pop1.list -p pop2.list\n"
        "\n";
}

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

	size_t tmp_popstr_len = 0;
	char *tmp_popstr;
	int idx;


	gd.hdr_in = in;
	gd.n_pops = 0;
	gd.n_generations = 0;
	gd.population_size = 0;
	gd.n_threads = 1;

	while ((c = getopt(argc, argv, "p:n:t:N:")) != -1)
	{
		switch (c) {
			case 'p':
				pop_filename = realloc(pop_filename, (gd.n_pops+1)*sizeof(char *));
				if (pop_filename == NULL) {
					perror("realloc");
					return -1;
				}
				pop_filename[gd.n_pops] = optarg;
				gd.n_pops++;
				continue;
			case 'n':
				errno = 0;
				gd.n_threads = strtoul(optarg, NULL, 0);
				if (errno) {
					fprintf(stderr, "Number of threads ``%s'' is invalid.\n", optarg);
					return -1;
				}
				continue;
			case 't':
				errno = 0;
				gd.n_generations = strtoul(optarg, NULL, 0);
				if (errno) {
					fprintf(stderr, "Number of generations ``%s'' is invalid.\n", optarg);
					return -1;
				}
				continue;
			case 'N':
				errno = 0;
				gd.population_size = strtoul(optarg, NULL, 0);
				if (errno) {
					fprintf(stderr, "Population size ``%s'' is invalid.\n", optarg);
					return -1;
				}
				continue;
			default:
				fprintf(stderr,"%s", usage());
				return -1;
		}
	}

	if (gd.n_pops != 2) {
		fprintf(stderr, "Must specify files containing the samples "
			       "in the two populations\nwhich are to be tested.\n");
		fprintf(stderr,"%s", usage());
		return -1;
	}

	if (gd.population_size == 0) {
		fprintf(stderr, "Must specify population size.\n");
		fprintf(stderr,"%s", usage());
		return -1;
	}

	if (gd.n_generations == 0) {
		fprintf(stderr, "Must specify number of generations.\n");
		fprintf(stderr,"%s", usage());
		return -1;
	}

	pthread_mutex_init(&gd.print_mutex, NULL);
	pthread_mutex_init(&gd.wq.mutex, NULL);
	gd.wq.max_size = gd.n_threads * 4;
	gd.wq.done = 0;
	gd.wq.q = kl_init(bcf_rec);


	samples = malloc(gd.n_pops*sizeof(char **));
	if (samples == NULL) {
		perror("malloc");
		return -1;
	}

	gd.samples_per_pop = malloc(gd.n_pops*sizeof(unsigned int));
	if (gd.samples_per_pop == NULL) {
		perror("malloc");
		return -1;
	}

	for (i=0; i<gd.n_pops; i++) {
		samples[i] = hts_readlist(pop_filename[i], 1, &gd.samples_per_pop[i]);
		if (samples[i] == NULL) {
			fprintf(stderr, "Error reading %s: %s\n",
					pop_filename[i], strerror(errno));
			return -1;
		}

		for (j=0; j<gd.samples_per_pop[i]; j++) 
			tmp_popstr_len += strlen(samples[i][j]);
		tmp_popstr_len += gd.samples_per_pop[i];
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

	for (i=0; i<gd.n_pops; i++) {
		for (j=0; j<gd.samples_per_pop[i]; j++) {
			strcat(tmp_popstr, samples[i][j]);
			strcat(tmp_popstr, ",");
		}
	}

	// remove trailing comma
	tmp_popstr[tmp_popstr_len-2] = '\0';

	idx = bcf_hdr_set_samples(gd.hdr_in, tmp_popstr, 0);
	if (idx != 0) {
		if (idx > 0) {
			idx--;
			/* 
			 * Find the name of the sample that is not in the vcf
			 * and the file in which it was specified.
			 */
			for (i=0; i<gd.n_pops; i++) {
				if (idx < gd.samples_per_pop[i])
					break;
				idx -= gd.samples_per_pop[i];
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
	gd.pop_indices = malloc(gd.n_pops*sizeof(unsigned int *));
	if (gd.pop_indices == NULL) {
		perror("malloc");
		return -1;
	}

	for (i=0; i<gd.n_pops; i++) {
		gd.pop_indices[i] = malloc(gd.samples_per_pop[i]*sizeof(unsigned int));
		if (gd.pop_indices[i] == NULL) {
			perror("malloc");
			return -1;
		}
		for (j=0; j<gd.samples_per_pop[i]; j++) {
			gd.pop_indices[i][j] = 2*bcf_hdr_id2int(gd.hdr_in, BCF_DT_SAMPLE, samples[i][j]);
			free(samples[i][j]);
		}
		free(samples[i]);
	}

	free(samples);



	printf("CHROM\tPOS\tID\tREF\tALT\tP");
	for (i=0; i<gd.n_pops; i++) {
		printf("\t[P%d]MAF (MAC/TAC)", i+1);
	}
	printf("\n");


	gd.threads = calloc(gd.n_threads, sizeof(pthread_t));
	if (gd.threads == NULL) {
		perror("calloc");
		return -1;
	}

	for (i=0; i<gd.n_threads; i++) {
		if (pthread_create(&gd.threads[i], NULL, &drift_worker, &gd.wq)) {
			perror("pthread_create");
			return -1;
		}
	}

	return 1;
}

/*
 * Pop a record off the queue.
 * Returns NULL if there are no more records.
 */
bcf1_t *
queue_pop_rec(work_queue_t *wq)
{
	bcf1_t *rec = NULL;

	do {
		if (pthread_mutex_lock(&wq->mutex)) {
			perror("queue_pop_rec: pthread_mutex_lock");
			exit(1);
		}

		kl_shift(bcf_rec, wq->q, &rec);

		if (pthread_mutex_unlock(&wq->mutex)) {
			perror("queue_pop_rec: pthread_mutex_unlock");
			exit(1);
		}

		if (rec == NULL) {
			if (wq->done)
				break;
			sched_yield();
		}

	} while (rec == NULL);

	return rec;
}

/*
 * Push a vcf record on to the queue.
 */
void
queue_push_rec(work_queue_t *wq, bcf1_t *rec)
{
	size_t size;

	do {
		if (pthread_mutex_lock(&wq->mutex)) {
			perror("queue_push_rec: pthread_mutex_trylock");
			exit(1);
		}

		size = wq->q->size;

		if (size < wq->max_size)
			*kl_pushp(bcf_rec, wq->q) = rec;

		if (pthread_mutex_unlock(&wq->mutex)) {
			perror("queue_push_rec: pthread_mutex_unlock");
			exit(1);
		}

		if (size == wq->max_size)
			sched_yield();

	} while (size == wq->max_size);
}


/*
 * Worker thread.
 *
 * Pop a record from the queue, obtain the genotype, run the drift test.
 * Repeat.
 */
void *
drift_worker(void *arg)
{
	work_queue_t *wq = arg;
	bcf1_t *rec;

	/* Genotype array, realloc()'d in bcf_get_genotypes(). */
	int *gt_arr = NULL;
	int ngt, ngt_arr = 0;
	unsigned int gt0, gt1;
	unsigned int ac1, ac2;
	unsigned int n1, n2;
	int idx, j;

	double p;

	while (1) {
		rec = queue_pop_rec(wq);
		if (rec == NULL)
			break;

		bcf_unpack(rec, BCF_UN_STR);

		ngt = bcf_get_genotypes(gd.hdr_in, rec, &gt_arr, &ngt_arr);
		if (ngt != ngt_arr) {
			fprintf(stderr, "Failed to get genotypes for %s:%d.\n",
					bcf_seqname(gd.hdr_in, rec), rec->pos+1);
			exit(1);
		}

		ac1 = 0;
		n1 = 0;
		for (j=0; j<gd.samples_per_pop[0]; j++) {
			idx = gd.pop_indices[0][j];
			if (bcf_gt_is_missing(gt_arr[idx]))
				continue;
			gt0 = bcf_gt_allele(gt_arr[idx]);
			gt1 = bcf_gt_allele(gt_arr[idx+1]);
			ac1 += gt0 + gt1;
			n1 += 2;
		}

		ac2 = 0;
		n2 = 0;
		for (j=0; j<gd.samples_per_pop[1]; j++) {
			idx = gd.pop_indices[1][j];
			if (bcf_gt_is_missing(gt_arr[idx]))
				continue;
			gt0 = bcf_gt_allele(gt_arr[idx]);
			gt1 = bcf_gt_allele(gt_arr[idx+1]);
			ac2 += gt0 + gt1;
			n2 += 2;
		}

		p = drift_test(ac1, n1, ac2, n2, gd.n_generations, gd.population_size);

		if (pthread_mutex_lock(&gd.print_mutex)) {
			perror("drift_worker: pthread_mutex_lock");
			exit(1);
		}

		printf("%s\t%d\t%s\t%s\t%s\t%lf\t%.3lf (%d/%d)\t%.3lf (%d/%d)\n",
			bcf_seqname(gd.hdr_in, rec), rec->pos+1, rec->d.id,
			rec->d.allele[0], rec->d.allele[1], p,
			((double)ac1)/n1, ac1, n1,
			((double)ac2)/n2, ac2, n2);

		if (pthread_mutex_unlock(&gd.print_mutex)) {
			perror("drift_worker: pthread_mutex_unlock");
			exit(1);
		}

		bcf_destroy(rec);
	}

	free(gt_arr);

	return NULL;
}

/*
 * Called for each VCF record.
 * Return rec to output the line or NULL to suppress output.
 */
bcf1_t *
process(bcf1_t *rec)
{
	// we only want biallelic SNPs
	if (!bcf_is_snp(rec) || rec->n_allele != 2)
		return NULL;

	queue_push_rec(&gd.wq, bcf_dup(rec));

	return NULL;
}

/*
 * Called after all lines have been processed to clean up.
 */
void
destroy(void)
{
	int i;

	// signal worker threads to exit
	gd.wq.done = 1;

	// await the completion of worker threads
	for (i=0; i<gd.n_threads; i++) {
		if (pthread_join(gd.threads[i], NULL)) {
			perror("pthread_join");
			exit(1);
		}
	}

	pthread_mutex_destroy(&gd.print_mutex);
	pthread_mutex_destroy(&gd.wq.mutex);

	free(gd.threads);

	free(gd.samples_per_pop);

	for (i=0; i<gd.n_pops; i++) {
		free(gd.pop_indices[i]);
	}

	free(gd.pop_indices);
}
