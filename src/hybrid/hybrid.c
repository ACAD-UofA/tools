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
 *
 *
 * Species tree
 *        /\        ---
 *       /  \        |
 *      /    \       t
 *     /      \      |
 *    /--------\    ---
 *   /  y | 1-y \
 *  /     |      \
 * A      B       C
 *
 * Gene tree 1, occurs with probability y.
 *    /\
 *   /  \
 *  /\   \
 * A  B   C
 *
 * Gene tree 2, occurs with probability 1-y.
 *    /\
 *   /  \
 *  /   /\
 * A   B  C
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <signal.h>
#include <float.h>
#include <math.h>
#include <pthread.h>

#include "kmath.h"

typedef uint gt3_t[3];

static int quit_signalled;


void
rmultinom(krand_t *kr, int N, double *p, uint (*n)[], size_t len)
{
	double p_cum[len];
	double u;
	int i;

	memset(*n, 0, sizeof(int)*len);

	if (N < 1)
		return;

	p_cum[0] = p[0];
	for (i=1; i<len; i++)
		p_cum[i] = p_cum[i-1] + p[i];

	while (--N) {
		u = kr_drand(kr);
		for (i=0; i<len; i++) {
			if (u < p_cum[i]) {
				(*n)[i]++;
				break;
			}
		}
	}
}

void
propose_pmultinom(krand_t *kr, double *last_p, double (*p)[], size_t len)
{
	double x;
	double p_sum;
	int i;
	do {
		p_sum = 0.0;
		for (i=0; i<len-1; i++) {
			do {
				x = kr_normal(kr);
				x = last_p[i] + x*0.01;
			} while (x < 0.0 || x > 1.0);
			(*p)[i] = x;
			p_sum += x;
		}
	} while (p_sum > 1.0);
	(*p)[len-1] = 1 - p_sum;
}

/*
 * pow(p[0],n[0])*...*pow(p[len-1],n[len-1]) * (n[0]+...+n[len-1])! / (n[0]!*...*n[len-1]!)
 */
double
lmultinom(int *n, double *p, size_t len)
{
	double lmn = 0.0;
	double n_sum = 0.0;

	while (--len) {
		lmn += n[len]*log(p[len]) - lgamma(n[len]+1);
		n_sum += n[len];
	}

	lmn += lgamma(n_sum+1);

	return lmn;
}

double
llik(double y, double t, int nab, int nbc, int nac)
{
	double ll;
	double exp_t = exp(-t);
	ll = nab * log(y + (1.0/3.0)*exp_t - y*exp_t)
		+ nbc * log(1-y-(2.0/3.0)*exp_t +y*exp_t)
		+ nac * (log(1.0/3.0) - t);
	return ll;
}

/*
 * Propose new t, time between hybridisation and top of tree,
 * on closed interval [0, TDIFF].
 */
double
t_rand(krand_t *kr, double last_t)
{
#define GEN_TIME 8.0
#define Ne 1e4
// 100kya, lower bound on steppe/wisent divergence
#define TM0_LOWER (1e5/GEN_TIME)/(2.0*Ne)
// 5.9mya, upper bound on steppe/wisent/aurochs divergence
#define TM1_UPPER (5.9e6/GEN_TIME)/(2.0*Ne)
#define TDIFF (TM1_UPPER-TM0_LOWER)
	double t;
	do {
		t = kr_normal(kr);
		t = last_t + t*0.01*TDIFF;
	} while (t < 0.0 || t > TDIFF);
	return t;
}

/*
 * Propose new y, on closed interval [0,1].
 */
double
y_rand(krand_t *kr, double last_y)
{
	double y;
	do {
		y = kr_normal(kr);
		y = last_y + y*0.01;
	} while (y < 0.0 || y > 1.0);
	return y;
}

int
mcmc(FILE *fp, krand_t *kr, int y, gt3_t n, int sampling_period, int samples, double *y_max, double *t_max, double *ll_max)
{
	double y0, y1;
	double t0, t1;
	double ll0, ll1, lldiff;
	int i = 0;

	int n_sum = n[0] + n[1] + n[2];
	double p0[3], p1[3];
	p0[0] = (double)n[0] / n_sum;
	p0[1] = (double)n[1] / n_sum;
	p0[2] = (double)n[2] / n_sum;
	memcpy(p1, p0, sizeof(double)*3);

	y0 = y1 = *y_max = (y == -1) ? 0.5 : y;
	t0 = *t_max = TDIFF/2.0;
	ll0 = *ll_max = -DBL_MAX;

	fprintf(fp, "y_current\ty_proposed\tt_current\tt_proposed\tloglik_current\tloglik_proposed\tp0_ab\tp0_bc\tp0_ac\tp1_ab\tp1_bc\tp1_ac\n");

	do {
		if (y == -1)
			y1 = y_rand(kr, y0);
		t1 = t_rand(kr, t0);

		ll1 = llik(y1, t1, n[0], n[1], n[2]);

		lldiff = ll1 - ll0;
		if (lldiff > 0.0 || log(kr_drand(kr)) < lldiff) {
			// accept new state
			ll0 = ll1;
			y0 = y1;
			t0 = t1;
			memcpy(p1, p0, sizeof(double)*3);

			if (ll1 > *ll_max) {
				*ll_max = ll1;
				*y_max = y1;
				*t_max = t1;
			}
		}

		i++;
		if (i%sampling_period == 0) {
			fprintf(fp, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
					y0, y1, t0, t1, ll0, ll1, p0[0], p0[1], p0[2], p1[0], p1[1], p1[2]);
			i = 0;
			if (samples != -1 && --samples == 0)
				break;
		}
	} while (!quit_signalled);

	return 0;
}

struct thread_stuff {
	pthread_t tid;
	FILE *fp;
	krand_t *kr;
	int y;
	gt3_t n;
	int sampling_period;
	int samples;
	double y_max, t_max, ll_max;
};

void *
mcmc_start(void *arg)
{
	struct thread_stuff *ts = arg;
	intptr_t ret;

	ret = mcmc(ts->fp, ts->kr, ts->y, ts->n, ts->sampling_period, ts->samples,
			&ts->y_max, &ts->t_max, &ts->ll_max);

	return (void *)ret;
}

void
handle_sigint(int sig)
{
	quit_signalled = 1;
}

int
compare_models(FILE *fp, struct thread_stuff ts[3], int *aic_i, int *aicc_i, int *bic_i)
{
	double aic, aicc, bic;
	double aic_min, aicc_min, bic_min;
	int i;
	int k; // number of model parameters
	int n; // number of sites
	char *models[] = {"Hybrid", "(AB)C", "A(BC)"};

	aic_min = aicc_min = bic_min = DBL_MAX;
	*aic_i = *aicc_i = *bic_i -1;

	fprintf(fp, "Model\tk\tN\tNab\tNbc\tNac\ty\tt\tloglik\tAIC\tAICc\tBIC\n");
	for (i=0; i<3; i++) {
		if (i==0)
			k = 2;
		else
			k = 1;
		n = ts[i].n[0] + ts[i].n[1] + ts[i].n[2];
		aic = -2*ts[i].ll_max + 2*k;
		aicc = aic + ((double)(2*k*(k+1))) / ((double)(n-k-1));
		bic = -2*ts[i].ll_max + k*log(n);

		if (aic < aic_min) {
			aic_min = aic;
			*aic_i = i;
		}
		if (aicc < aicc_min) {
			aicc_min = aicc;
			*aicc_i = i;
		}
		if (bic < bic_min) {
			bic_min = bic;
			*bic_i = i;
		}

		fprintf(fp, "%s\t%d\t%d\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
				models[i], k, n,
				ts[i].n[0], ts[i].n[1], ts[i].n[2],
				ts[i].y_max, ts[i].t_max, ts[i].ll_max,
				aic, aicc, bic);
	}

	return 0;
}

int
do_mcmc(FILE *fp, gt3_t n, int sampling_period, int samples, int *aic_i, int *aicc_i, int *bic_i, double *gamma)
{
	char fn[1024];
	struct thread_stuff ts[3] = {[0].y = -1, [1].y = 1, [2].y = 0 };
	pthread_attr_t attr;
	void *res;
	int i;
	int ret = 0;

	for (i=0; i<3; i++) {
		ts[i].kr = NULL;
		memcpy(ts[i].n, n, sizeof(gt3_t));
		ts[i].sampling_period = sampling_period;
		ts[i].samples = samples;

		sprintf(fn, "chain_y%d.txt", ts[i].y);
		ts[i].fp = fopen(fn, "w");
		if (ts[i].fp == NULL) {
			fprintf(stderr, "Error: %s: %s\n", fn, strerror(errno));
			free(ts[i].kr);
			ret = -1*(i+1);
			goto err0;
		}

		//ts[i].kr = kr_srand(getpid()+i);
		ts[i].kr = kr_srand(122+i);
	}

	if (pthread_attr_init(&attr) != 0) {
               perror("pthread_attr_init");
	       ret = -4;
	       goto err0;
	}

	quit_signalled = 0;
	signal(SIGINT, handle_sigint);
	signal(SIGTERM, handle_sigint);

	for (i=0; i<3; i++) {
		if (pthread_create(&ts[i].tid, &attr, mcmc_start, &ts[i]) != 0) {
			perror("pthread_create");
			quit_signalled = 1;
			ret = -2;
			break;
		}
	}

	pthread_attr_destroy(&attr);

	// wait for threads to end
	i--;
	for (; i>=0; i--) {
		if (pthread_join(ts[i].tid, &res) != 0) {
			perror("pthread_join");
			ret = -3;
		}
	}

	if (ret == 0)
		compare_models(fp, ts, aic_i, aicc_i, bic_i);

	*gamma = ts[0].y_max;

err0:
	for (i=0; i<3; i++) {
		if (ts[i].kr){
			fclose(ts[i].fp);
			free(ts[i].kr);
		}
	}

	return ret;
}



// print likelihood surfaces for a y*t grid.
int
llik_surfaces(gt3_t n)
{
	double y;
	double t;

	FILE *fp = fopen("llik.txt", "w");
	if (fp == NULL) {
		fprintf(stderr, "Error: llik.txt: %s\n", strerror(errno));
		return -1;
	}

	fprintf(fp, "y\tt\t(AB,C)\t(A,BC)\t(B,AC)\tData\n");

	for (y=0.0; y<=1.0; y+=0.01) {
		for (t=0; t<=TDIFF; t+=(TDIFF/20.0)) {
			fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
					y,
					t,
					llik(y, t, 1, 0, 0),
					llik(y, t, 0, 1, 0),
					llik(y, t, 0, 0, 1),
					llik(y, t, n[0], n[1], n[2]));
		}
	}

	fclose(fp);

	return 0;
}

#ifdef _HYBRID_MAIN
int
main(int argc, char **argv)
{
	gt3_t n;
	int sampling_period, samples;
	int i;
	int aic, aicc, bic;
	double gamma;

	if (argc != 6) {
		fprintf(stderr, "%s sampling_period samples Nab Nbc Nac\n", argv[0]);
		return -1;
	}

	errno = 0;
	sampling_period = strtol(argv[1], NULL, 0);
	if (errno || sampling_period < 1 || sampling_period > 1e7) {
		fprintf(stderr, "sampling_period=`%s' not on the interval [1,1e7].\n", argv[1]);
		return -3;
	}

	errno = 0;
	samples = strtol(argv[2], NULL, 0);
	if (errno || samples < 0) {
		fprintf(stderr, "samples=`%s' not a positive number.\n", argv[2]);
		return -3;
	}

	/* sample forever */
	if (samples == 0)
		samples = -1;

	for (i=0; i<3; i++) {
		int nn;
		errno = 0;
		nn = strtol(argv[3+i], NULL, 0);
		if (errno || nn < 0) {
			char *sub[] = {"ab", "bc", "ac"};
			fprintf(stderr, "N%s=`%s' not a positive number.\n", sub[i], argv[3+i]);
			return -3;
		}
		n[i] = nn;
	}

	if (llik_surfaces(n) < 0)
		return -5;

	if (do_mcmc(stdout, n, sampling_period, samples, &aic, &aicc, &bic, &gamma) < 0)
		return -6;

	return 0;
}
#endif

#ifdef _POWER_MAIN
int
main()
{
	int i, j, k;
	double gamma, gamma_est;
	krand_t *kr;

	//kr = kr_srand(getpid()+i);
	kr = kr_srand(119);

	gt3_t n1, n2, ntot;
	double p1[3] = {1000.0/1200, 100.0/1200, 100.0/1200};
	double p2[3] = {100.0/1200, 1000.0/1200, 100.0/1200};

	// ensure probabilities add up to 1.0
	p1[2] = 1.0 - (p1[0] + p1[1]);
	p2[2] = 1.0 - (p2[0] + p2[1]);

	const int N = 1200;
	const int Niter = 1000;

	int aic, aicc, bic;
	int aic_sum[3], aicc_sum[3], bic_sum[3];

	double sq_err_sum;

	FILE *null = fopen("/dev/null", "w");

	printf("y\ty_RMSE\tAIC[hybrid]\tAIC[ab]\tAIC[bc]\tAICc[hybrid]\tAICc[ab]\tAICc[bc]\tBIC[hybrid]\tBIC[ab]\tBIC[bc]\n");

	for (gamma=1.0; gamma>=0.80; gamma-=0.01) {
		memset(aic_sum, 0, sizeof(aic_sum));
		memset(aicc_sum, 0, sizeof(aicc_sum));
		memset(bic_sum, 0, sizeof(bic_sum));
		sq_err_sum = 0.0;
		for (i=0; i<Niter; i++) {
			rmultinom(kr, gamma*N, p1, &n1, 3);
			rmultinom(kr, (1-gamma)*N, p2, &n2, 3);
			for (j=0; j<3; j++)
				ntot[j] = n1[j] + n2[j];
			do_mcmc(null, ntot, 100, 1000, &aic, &aicc, &bic, &gamma_est);
			aic_sum[aic]++;
			aicc_sum[aicc]++;
			bic_sum[bic]++;

			sq_err_sum += (gamma-gamma_est)*(gamma-gamma_est);
		}

		printf("%.3lf\t%.3lg", gamma, sqrt(sq_err_sum/Niter));
		for (k=0; k<3; k++)
			printf("\t%.3lf", (double)aic_sum[k]/Niter);
		for (k=0; k<3; k++)
			printf("\t%.3lf", (double)aicc_sum[k]/Niter);
		for (k=0; k<3; k++)
			printf("\t%.3lf", (double)bic_sum[k]/Niter);
		printf("\n");
	}

	fclose(null);

	return 0;
}
#endif
