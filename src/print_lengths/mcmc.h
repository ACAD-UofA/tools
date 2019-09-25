#ifndef _MCMC_H
#define _MCMC_H

typedef struct mcmc_param {
	char name[64];
	double val, x; // accepted, current
	double val_lmax; // value at maximum likelihood
	double range_min, range_max, range;
	double (*proposal_func)(krand_t *, struct mcmc_param *);
	double sigma; // proposal deviation
	krand_t *kr;
} mcmc_param_t;

typedef double (*llik_func_t)(void *data, mcmc_param_t *params, int n_params);

int mcmc(FILE *fp, krand_t *kr, int sampling_period, int samples,
	llik_func_t llik, void *data, mcmc_param_t *params, int n_params);
void param_init(mcmc_param_t *p, const char *name,
		double init, double min, double max, double sigma);

#endif //_MCMC_H
