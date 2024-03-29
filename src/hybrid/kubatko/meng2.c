/*
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

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>
#include <float.h>
#include <math.h>
#include <pthread.h>

#include <htslib/hts.h>
#include <htslib/vcf.h>

#include "kmath.h"

typedef struct {
	uint nab;
	uint nbc;
	uint nac;
	int nsites;
} gt3_t;

static int quit_signalled;

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
 * Propose new t.
 *
 * Sample a time between hybridisation and top of tree,
 * in units of generations per effective population size.
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
	double t = kr_normal(kr);
	t = last_t + t*0.01*TDIFF;
	if (t < 0.0)
		t = -t;
	if (t > TDIFF)
		t = 2*TDIFF - t;
	return t;
}

/*
 * Propose new y, on closed interval [0,1].
 */
double
y_rand_closed(krand_t *kr, double last_y)
{
	double y = kr_normal(kr);
	y = last_y + y*0.01;
	if (y < 0.0)
		y = -y;
	else if (y > 1.0)
		y = 2.0 -y;
	return y;
}

/*
 * Propose new y, on open interval (0,1).
 */
double
y_rand_open(krand_t *kr, double last_y)
{
	double y = kr_normal(kr);
	y = last_y + y*0.01;
	if (y < 0.0)
		y = 0.0;
	else if (y > 1.0)
		y = 1.0;
	return y;
}

int
mcmc(FILE *fp, krand_t *kr, int y, gt3_t *abc, int sampling_period, int samples, double *y_max, double *t_max, double *ll_max)
{
	double y0, y1;
	double t0, t1;
	double ll0, ll1, lldiff;
	int i = 0;

	int nab = abc->nab;
	int nbc = abc->nbc;
	int nac = abc->nac;

	y0 = y1 = *y_max = (y == -1) ? 0.5 : y;
	t0 = *t_max = TDIFF/2.0;
	ll0 = *ll_max = -DBL_MAX;

	fprintf(fp, "y_current\ty_proposed\tt_current\tt_proposed\tloglik_current\tloglik_proposed\n");

	do {
		if (y == -1)
			y1 = y_rand_closed(kr, y0);
		t1 = t_rand(kr, t0);
		ll1 = llik(y1, t1, nab, nbc, nac);

		lldiff = ll1 - ll0;
		if (lldiff > 0.0 || log(kr_drand(kr)) < lldiff) {
			// accept new state
			ll0 = ll1;
			y0 = y1;
			t0 = t1;

			if (ll1 > *ll_max) {
				*ll_max = ll1;
				*y_max = y1;
				*t_max = t1;
			}
		}

		i++;
		if (i%sampling_period == 0) {
			fprintf(fp, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
					y0, y1, t0, t1, ll0, ll1);
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
	gt3_t *abc;
	int sampling_period;
	int samples;
	double y_max, t_max, ll_max;
};

void *
mcmc_start(void *arg)
{
	struct thread_stuff *ts = arg;
	intptr_t ret;

	ret = mcmc(ts->fp, ts->kr, ts->y, ts->abc, ts->sampling_period, ts->samples,
			&ts->y_max, &ts->t_max, &ts->ll_max);

	return (void *)ret;
}

void
handle_sigint(int sig)
{
	quit_signalled = 1;
}

int
compare_models(struct thread_stuff ts[3])
{
	double aic, aicc, bic;
	int i;
	int k; // number of model parameters
	int n; // number of sites
	char *models[] = {"Hybrid", "(AB)C", "A(BC)"};


	printf("Model\tk\tN\tNab\tNbc\tNac\ty\tt\tloglik\tAIC\tAICc\tBIC\n");
	for (i=0; i<3; i++) {
		if (i==0)
			k = 2;
		else
			k = 1;
		n = ts[i].abc->nab + ts[i].abc->nbc + ts[i].abc->nac;
		aic = -2*ts[i].ll_max + 2*k;
		aicc = aic + ((double)(2*k*(k+1))) / ((double)(n-k-1));
		bic = -2*ts[i].ll_max + k*log(n);
		printf("%s\t%d\t%d\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
				models[i], k, n,
				ts[i].abc->nab, ts[i].abc->nbc, ts[i].abc->nac,
				ts[i].y_max, ts[i].t_max, ts[i].ll_max,
				aic, aicc, bic);
	}

	return 0;
}

int
do_mcmc(gt3_t *abc, int sampling_period, int samples)
{
	char fn[1024];
	struct thread_stuff ts[3] = {[0].y = -1, [1].y = 1, [2].y = 0 };
	pthread_attr_t attr;
	void *res;
	int i;
	int ret = 0;

	for (i=0; i<3; i++) {
		ts[i].kr = NULL;
		ts[i].abc = abc;
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
		compare_models(ts);

err0:
	for (i=0; i<3; i++) {
		if (ts[i].kr){
			fclose(ts[i].fp);
			free(ts[i].kr);
		}
	}

	return ret;
}


int
parse_pops(bcf_hdr_t *hdr, char *pops_fn, char *groups[4], int *(*si)[4], int (*ni)[4])
{
	FILE *fp;
	int ret = 0;
	void *tmp;
	int i;
	int idx;

	char *ind;
	char *grp;

	char *p;
	char *line = NULL;
	size_t len = 0;
	int lineno = 0;
	ssize_t read;

	memset(si, 0, sizeof(int *)*4);
	memset(ni, 0, sizeof(int *)*4);

	fp = fopen(pops_fn, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error opening %s: %s\n",
				pops_fn, strerror(errno));
		ret = -1;
		goto err0;
	}

	while ((read = getline(&line, &len, fp)) != -1) {
		lineno++;

		if (len == 0 || line[0] == '#')
			continue;

		if (line[strlen(line)-1] == '\n')
			line[strlen(line)-1] = '\0';

		for (p=line; p-line<len; p++) {
			if (*p == '\t')
				break;
		}
		if (*p != '\t') {
			fprintf(stderr, "Error: %s: line %d: missing tab\n",
					pops_fn, lineno);
			ret = -2;
			goto err1;
		}

		*p = '\0';
		p++;

		ind = line;
		grp = p;

		for (i=0; i<4; i++) {
			if (strcmp(groups[i], grp) == 0)
				break;
		}
		if (i==4)
			// group not in our A,B,C,D
			continue;

		idx = 2*bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, ind);
		if (idx < 0) {
			fprintf(stderr, "Error: %s: no sample ``%s'' in vcf.\n",
					pops_fn, ind);
			ret = -3;
			goto err1;
		}

		(*ni)[i]++;

		tmp = realloc((*si)[i], sizeof(int) * (*ni)[i]);
		if (tmp == NULL) {
			perror("parse_pops: realloc");
			ret = -4;
			goto err1;
		}
		(*si)[i] = tmp;
		(*si)[i][(*ni)[i] -1] = idx;
	}

	if (!feof(fp)) {
		fprintf(stderr, "Error reading %s: %s\n",
				pops_fn, strerror(errno));
		ret = -4;
		goto err1;
	}

	for (i=0; i<4; i++) {
		if ((*ni)[i] == 0) {
			fprintf(stderr, "Error: %s: no individuals for group ``%s''\n",
					pops_fn, groups[i]);
			ret = -5;
			goto err1;
		}
	}

	ret = 0;
err1:
	if (ret < 0) {
		for (i=0; i<4; i++) {
			if ((*si)[i]) {
				free((*si)[i]);
				(*si)[i] = NULL;
				(*ni)[i] = 0;
			}
		}
	}

	if (line)
		free(line);
	fclose(fp);
err0:
	return ret;
}

/*
 * Use the smallest F2 estimator (Patterson et al. 2012, Appendix A) to
 * inform the gene tree topology.
 */
int
topo_from_f2(int n_refs[4], int n_alts[4], int n_tot[4])
{
	double h[4]; // heterozygosity
	double af[4]; // ref allele frequency
	double f2_ab, f2_bc, f2_ac;
#define USE_OUTGROUP
#ifdef USE_OUTGROUP
	double f2_ad, f2_bd, f2_cd;
#endif
	int i;

	for (i=0; i<4; i++) {
		h[i] = ((double)(n_refs[i]*n_alts[i])) / (n_tot[i] *n_tot[i] * (n_tot[i] -1));
		af[i] = (double)n_refs[i] / n_tot[i];
		printf("h[%d]=%.3lf, af[%d]=%.3lf\n", i, h[i], i, af[i]);
	}

	f2_ab = (af[0]-af[1])*(af[0]-af[1]) - h[0] -h[1];
	f2_bc = (af[1]-af[2])*(af[1]-af[2]) - h[1] -h[2];
	f2_ac = (af[0]-af[2])*(af[0]-af[2]) - h[0] -h[2];

#ifdef USE_OUTGROUP
	f2_ad = (af[0]-af[3])*(af[0]-af[3]) - h[0] -h[3];
	f2_bd = (af[1]-af[3])*(af[1]-af[3]) - h[1] -h[3];
	f2_cd = (af[2]-af[3])*(af[2]-af[3]) - h[2] -h[3];
#endif

	printf("f2_ab=%.3lf, f2_bc=%.3lf, f2_ac=%.3lf, f2_ad=%.3lf, f2_bd=%.3lf, f2_cd=%.3lf\n",
			f2_ab, f2_bc, f2_ac, f2_ad, f2_bd, f2_cd);

	if (f2_ab < f2_bc && f2_ab < f2_ac) {
#ifdef USE_OUTGROUP
		double f2 = f2_ab;
		//if (!(f2 < f2_ad && f2 < f2_bd && f2 < f2_cd))
		if (!(f2 < f2_ad && f2 < f2_bd))
			return -1;
#endif
		// ((AB),C)
		return 0;
	} else if (f2_bc < f2_ab && f2_bc < f2_ac) {
#ifdef USE_OUTGROUP
		double f2 = f2_bc;
		//if (!(f2 < f2_ad && f2 < f2_bd && f2 < f2_cd))
		if (!(f2 < f2_bd && f2 < f2_cd))
			return -1;
#endif
		// (A,(BC))
		return 1;
	} else if (f2_ac < f2_ab && f2_ac < f2_bc) {
#ifdef USE_OUTGROUP
		double f2 = f2_ac;
		//if (!(f2 < f2_ad && f2 < f2_bd && f2 < f2_cd))
		if (!(f2 < f2_ad && f2 < f2_cd))
			return -1;
#endif
		// ((AC),B)
		return 2;
	} else {
		// site is uninformative
		return -1;
	}
}

int
topo_from_af(int n_refs[4], int n_alts[4], int n_tot[4])
{
	double af[4]; // ref allele frequency
	double ab, bc, ac;
//#define USE_OUTGROUP
#ifdef USE_OUTGROUP
	double ad, bd, cd;
#endif
	int i;

	for (i=0; i<4; i++) {
		af[i] = (double)n_refs[i] / n_tot[i];
	}

#define abs(x) ((x) < 0 ? -(x) : (x))
#define sq(x) ((x)*(x))

	ab = sq(af[0]-af[1]);
	bc = sq(af[1]-af[2]);
	ac = sq(af[0]-af[2]);

#ifdef USE_OUTGROUP
	ad = sq(af[0]-af[3]);
	bd = sq(af[1]-af[3]);
	cd = sq(af[2]-af[3]);
#endif

	printf("ab=%.3lf, bc=%.3lf, ac=%.3lf, ad=%.3lf, bd=%.3lf, cd=%.3lf\n",
			ab, bc, ac, ad, bd, cd);

	if (ab < bc && ab < ac) {
#ifdef USE_OUTGROUP
		if (!(ab < ad && ab < bd && ab < cd))
			return -1;
#endif
		// ((AB),C)
		return 0;
	} else if (bc < ab && bc < ac) {
#ifdef USE_OUTGROUP
		if (!(bc < ad && bc < bd && bc < cd))
			return -1;
#endif
		// (A,(BC))
		return 1;
	} else if (ac < ab && ac < bc) {
#ifdef USE_OUTGROUP
		if (!(ac < ad && ac < bd && ac < cd))
			return -1;
#endif
		// ((AC),B)
		return 2;
	} else {
		// site is uninformative
		return -1;
	}
}

/*
 * Parse vcf file, extracting genotypes for specified samples.
 *
 * vcf_fn, filename of vcf or bcf.
 * pops_fn, filename of file containing individual to group mappings.
 * groups[4], list of population group names, for tree tips A, B, C, D.
 * abc, counts of the observed gene topologies.
 */
int
get_gts_from_vcf(char *vcf_fn, char *pops_fn, char *groups[4], gt3_t *abc)
{
	htsFile *fp;
	bcf_hdr_t *hdr;
	bcf1_t *rec;
	int *si[4]; // list of sample ids for each group
	int ni[4]; // number of ids per group
	int ret;
	int ngt, ngt_arr = 0;
	int *gt_arr = NULL;

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

	if (parse_pops(hdr, pops_fn, groups, &si, &ni) < 0) {
		ret = -3;
		goto err2;
	}

	rec = bcf_init1();
	if (rec == NULL) {
		perror("bcf_init1: calloc");
		ret = -4;
		goto err2;
	}

	int idx;
	int skip;
	int i, j;
	int ref;
	int n_ref[4], n_alt[4], n_tot[4];

	int uninf = 0;
	abc->nab = 0;
	abc->nbc = 0;
	abc->nac = 0;

	/*
	for (i=0; i<4; i++) {
		printf("%s[i=%d,ni=%d]: ", groups[i], i, ni[i]);
		for (j=0; j<ni[i]; j++)
			printf("%d,", si[i][j]);
		printf("\n");
	}*/

	while (bcf_read1(fp, hdr, rec) >= 0) {
		// only take biallelic SNPs
		if (!bcf_is_snp(rec) || rec->n_allele != 2)
			continue;

		// filter on QUAL
		if (rec->qual < 25)
			continue;

		// ignore transitions
		if ((rec->d.allele[0][0] == 'C' && rec->d.allele[1][0] == 'T')
			|| (rec->d.allele[1][0] == 'C' && rec->d.allele[0][0] == 'T'))
			continue;

		ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
		if (ngt != ngt_arr) {
			fprintf(stderr, "Error: %s: missing genotypes at %s:%d.\n",
					vcf_fn, bcf_seqname(hdr, rec), rec->pos+1);
			ret = -5;
			goto err3;
		}

		skip = 0;
		for (i=0; i<3; i++) {
			n_ref[i] = n_alt[i] = n_tot[i] = 0;
			for (j=0; j<ni[i]; j++) {
				idx = si[i][j];

				if (bcf_gt_is_missing(gt_arr[idx]))
					continue;

				ref = bcf_gt_allele(gt_arr[idx]) == 0;
				ref += bcf_gt_allele(gt_arr[idx+1]) == 0;

				n_ref[i] += ref;
				n_alt[i] += 2-ref;
				n_tot[i] += 2;
			}

			if (n_tot[i] == 0) {
				// skip missing
				skip = 1;
				break;
			}
		}

		if (skip)
			continue;

		int topo;
		switch ((topo = topo_from_af(n_ref, n_alt, n_tot))) {
			case -1:
				// TODO: add 1/3 to each?
				uninf++;
				break;
			case 0:
				abc->nab++;
				break;
			case 1:
				abc->nbc++;
				break;
			case 2:
				abc->nac++;
				break;
		}

		printf("%s:%d, topo=%d, [A=%.3lf, B=%.3lf, C=%.3lf, D=%.3lf]\n",
				bcf_seqname(hdr, rec), rec->pos+1, topo,
					(double)n_ref[0]/n_tot[0],
					(double)n_ref[1]/n_tot[1],
					(double)n_ref[2]/n_tot[2],
					(double)n_ref[3]/n_tot[3]
				);
	}

	abc->nab += (uninf)/3;
	abc->nbc += (uninf)/3;
	abc->nac += (uninf)/3;


	ret = 0;
err3:
	if (gt_arr != NULL)
		free(gt_arr);

	bcf_destroy1(rec);
err2:
	bcf_hdr_destroy(hdr);
err1:
	hts_close(fp);
err0:
	return ret;
}

// print likelihood surfaces for a y*t grid.
int
llik_surfaces(gt3_t *abc)
{
	double y;
	double t;

	int nab = abc->nab;
	int nbc = abc->nbc;
	int nac = abc->nac;

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
					llik(y, t, nab, nbc, nac));
		}
	}

	fclose(fp);

	return 0;
}

int
main(int argc, char **argv)
{
	gt3_t abc;
	int sampling_period, samples;

	if (argc != 9) {
		fprintf(stderr, "%s sampling_period samples file.vcf pops.txt A B C D\n", argv[0]);
		return -1;
	}

	errno = 0;
	sampling_period = strtoul(argv[1], NULL, 0);
	if (errno || sampling_period < 1 || sampling_period > 1e7) {
		fprintf(stderr, "sampling_period=%s not on the interval [1,1e7].\n", argv[1]);
		return -3;
	}

	errno = 0;
	samples = strtol(argv[2], NULL, 0);
	if (errno) {
		fprintf(stderr, "unable to convert samples (%s) to a number.\n", argv[2]);
		return -3;
	}

	if (get_gts_from_vcf(argv[3], argv[4], argv+5, &abc) < 0)
		return -4;

	//printf("Nab=%d, Nbc=%d, Nac=%d\n", abc.nab, abc.nbc, abc.nac);

	if (llik_surfaces(&abc) < 0)
		return -5;

	if (do_mcmc(&abc, sampling_period, samples) < 0)
		return -6;

	return 0;
}

