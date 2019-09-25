#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <errno.h>

// erand48() seed/state
static unsigned short xsubi[3];

double
af_walk(double p, unsigned int N, unsigned int t)
{
	unsigned int i, j;
	unsigned int mac = 0;

	for (i=0; i<t; i++) {
		mac = 0;
		for (j=0; j<N; j++) {
			if (erand48(xsubi) < p)
				mac++;
		}

		p = ((double)mac) / N;
	}

	return p;
}

int
cmp_d(const void *a, const void *b)
{
	double aa = *(double *)a;
	double bb = *(double *)b;
	return aa < bb ? -1 : 1;
}

#define ITERATIONS 1000
void
af_walk_bootstrap(double p, unsigned int N, unsigned int t,
		double *mean, double *median, double *ci_low, double *ci_high)
{
	double af[ITERATIONS];
	double sum = 0.0;
	int i;

	for (i=0; i<ITERATIONS; i++) {
		af[i] = af_walk(p, N, t);
		sum += af[i];
	}
	
	qsort(af, ITERATIONS, sizeof(double), cmp_d);

	*mean = sum / ITERATIONS;
	*median = (af[ITERATIONS/2] + af[ITERATIONS/2 -1]) / 2;
	*ci_low = af[ITERATIONS/40 -1];
	*ci_high = af[ITERATIONS - ITERATIONS/40 -1];
}

int
main(int argc, char **argv)
{
	double p; // initial allele frequency
	unsigned int N; // population size
	unsigned int t; // time, in generations

	double mean, median, ci_low, ci_high;

	if (argc != 4) {
		fprintf(stderr, "usage: %s p N t\n", argv[0]);
		return 1;
	}

	errno = 0;
	p = strtod(argv[1], NULL);
	if (errno || p < 0 || p > 1) {
		fprintf(stderr, "p=%s is no good.\n", argv[1]);
		return 2;
	}

	errno = 0;
	N = strtoul(argv[2], NULL, 0);
	if (errno) {
		fprintf(stderr, "N=%s is no good.\n", argv[2]);
		return 2;
	}

	
	errno = 0;
	t = strtoul(argv[3], NULL, 0);
	if (errno) {
		fprintf(stderr, "t=%s is no good.\n", argv[3]);
		return 2;
	}

	// seed random number generator
	xsubi[0] = time(NULL) & 0xffff;
	xsubi[1] = time(NULL) >> 16;
	xsubi[2] = getpid();

	af_walk_bootstrap(p, N, t, &mean, &median, &ci_low, &ci_high);

	printf("%lf, %lf, %lf, %lf\n", mean, median, ci_low, ci_high);

	return 0;
}
