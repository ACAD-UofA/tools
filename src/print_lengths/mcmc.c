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
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>
#include <float.h>
#include <math.h>
#include <assert.h>

#include "kmath.h"
#include "mcmc.h"

double
propose_normal(krand_t *kr, mcmc_param_t *p)
{
	double x;
	do {
		x = kr_normal(kr);
		x = p->val +x*p->sigma*p->range;
	} while (x < p->range_min || x > p->range_max);
	return x;
}

void param_init(mcmc_param_t *p, const char *name,
		double init, double min, double max, double sigma)
{
	strncpy(p->name, name, 63);
	p->name[63] = '\0';
	p->val = p->x = init;
	p->range_min = min;
	p->range_max = max;
	p->range = max-min;
	assert(p->range > 0);
	p->sigma = sigma;
	p->proposal_func = propose_normal;
}


int
//mcmc(FILE *fp, krand_t *kr, int y, gt3_t n, int sampling_period, int samples, double *y_max, double *t_max, double *ll_max)
mcmc(FILE *fp, krand_t *kr, int sampling_period, int samples, llik_func_t llik, void *data, mcmc_param_t *params, int n_params)
{
	int i = 0;
	double ll0, ll1, ll_max, ll_diff;

	ll0 = ll_max = -DBL_MAX;

	fprintf(fp, "ll0\tll1");
	for (i=0; i<n_params; i++)
		fprintf(fp, "\t%s0\t%s1", params[i].name, params[i].name);
	fprintf(fp, "\n");

	for (;;) {
		for (i=0; i<n_params; i++)
			params[i].x = params[i].proposal_func(kr, &params[i]);

		ll1 = llik(data, params, n_params);

		ll_diff = ll1 - ll0;
		if (ll_diff > 0.0 || log(kr_drand(kr)) < ll_diff) {
			// accept new state
			ll0 = ll1;
			for (i=0; i<n_params; i++) {
				params[i].val = params[i].x;
				if (ll1 > ll_max)
					params[i].val_lmax = params[i].x;
			}
		}

		i++;
		if (i%sampling_period == 0) {
			fprintf(fp, "%lg\t%lg", ll0, ll1);
			for (i=0; i<n_params; i++)
				fprintf(fp, "\t%lg\t%lg", params[i].val, params[i].x);
			fprintf(fp, "\n");
			i = 0;
			if (samples != -1 && --samples == 0)
				break;
		}
	}

	return 0;
}

