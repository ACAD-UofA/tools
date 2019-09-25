/*
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
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <float.h>

#include "cephes.h"

// Precision related epsilon.
#define EPS 1e-30

// Number of terms used before giving up on convergence.
#define TERMS_FOR_INFINITE_SUM 1500

/*
 * Beta probability distribution function, at x.
 */
static double
beta_pdf(unsigned int a, unsigned int b, double x)
{
	return exp((a-1)*log(x) + (b-1)*log1p(-x) - (lgamma(a) + lgamma(b) - lgamma(a+b)));
}

/*
 * Kimura (1955), Eq 5.
 * Kimura (1964), Eq 4.15.
 *
 * Returns the probability of fixation of an allele with frequency p
 * after t generations.
 *
 * For large values of t/N, the probability of fixation approaches p.
 *
 *  p -- Inital allele frequency.
 *  t -- Number of generations.
 *  N -- Population size.
 */
static double
fixation(double p, unsigned int t, unsigned int N)
{
	double u_last, u;
	int i;

	u_last = u = 0.0;

	for (i=1; i<=TERMS_FOR_INFINITE_SUM; i++) {
		u += (i%2?-1:1) * (2*i+1) * cephes_hyp2f1(1-i,i+2,2,p) * exp(-((double)(i*(i+1)*t))/(4*N));
		if (fabs(u - u_last) < EPS)
			// we have converged
			break;
		u_last = u;
	}

	u = p + u * p * (1-p);
	if (u < EPS)
		u = 0.0;

	return u;
}

/*
 * Kimura (1955), Eq 15'.
 * Kimura (1964), Eq 4.10.
 *
 * Returns the probability that a population of size N, with initial
 * allele frequency p, will drift to frequency x after t generations.
 *
 * For large values of t/N the probability approaches 0 for x!=0 && x!=1.
 *
 *  p -- Initial allele frequency.
 *  x -- Final allele frequency.
 *  t -- Number of generations.
 *  N -- Population size.
 */
static double
phi(double p, double x, unsigned int t, unsigned int N)
{
	double phi_last, phi;
	double i_ip1;
	int i;

	/*
	 * Handle the absorbing states.
	 */
	if (x < EPS)
		// allele becomes fixed
		return fixation(p, t, N);
	else if (x > 1.0 - EPS)
		// allele is lost
		return fixation(1-p, t, N);

	phi_last = phi = 0.0;

	for (i=1; i<=TERMS_FOR_INFINITE_SUM; i++) {
		i_ip1 = i * (i+1);
		phi += i_ip1 * (2*i+1) * cephes_hyp2f1(1-i,i+2,2,p) * cephes_hyp2f1(1-i,i+2,2,x) * exp(-(i_ip1*t)/(4*N));
		if (fabs(phi - phi_last) < EPS)
			// we have converged
			break;
		phi_last = phi;
	}

	phi = phi * p * (1-p);

	if (phi < EPS)
		phi = 0.0;

	return phi;
}

/*
 * Given a sampled population of size N, with observed minor allele count ac,
 * and total number of observed alleles n, this function returns the
 * probability that the population's allele frequency will be x after
 * t generations.
 *
 * integral_wrt_p[ phi(p,x,t,N) * Beta_PDF(1+ac,1+(n-ac),p) * dp ]
 *
 * ac -- Observed minor allele count.
 *  n -- Total alleles observed.
 *  x -- Final allele frequency.
 *  t -- Number of generations.
 *  N -- Population size.
 */
static double
drift(unsigned int ac, unsigned int n, double x, unsigned int t, unsigned int N)
{
	double int_phi_p_dp = 0.0;
	double p, pp;

// size of our increment in p for the integral
#define DRIFT_DP 5e-3

	/*
	 * We don't consider p=0, or p=1. The site is assumed biallelic.
	 */
	for (p=DRIFT_DP; p<=1.0-DRIFT_DP; p+=DRIFT_DP) {
		//pp = exp(ac*log(p) + (n-ac)*log1p(-p));
		pp = cephes_btdtr(1+ac, 1+n-ac, p+DRIFT_DP) - cephes_btdtr(1+ac, 1+n-ac, p);
		int_phi_p_dp += phi(p, x, t, N)*pp;
	}

	int_phi_p_dp *= DRIFT_DP;

	// divide by beta function (constant wrt p)
	//int_phi_p_dp /= exp(lgamma(1+ac) + lgamma(1+(n-ac)) - lgamma(2+n));


	return int_phi_p_dp;
}

/*
 * Test if the allele frequency of a bialleleic site in population 1 has
 * drifted to population 2's allele frequency over t generations.
 * True allele frequencies of both populations are unknown, but are
 * specified using Beta distributions from the observed allele counts.
 *
 * Returns the probability of rejecting the null hypothesis.
 *
 * mac1 -- minor allele count for sample of population 1.
 * tac1 -- total allele count for sample of population 1.
 * mac2 -- minor allele count for sample of population 2.
 * tac2 -- total allele count for sample of population 2.
 *    t -- time in generations.
 *    N -- population size (assumed constant w.r.t. time).
 * 
 */
double
drift_test(unsigned int mac1, unsigned int tac1,
		unsigned int mac2, unsigned int tac2,
		unsigned int t, unsigned int N)
{
	double x;
	double p_p1_x, p_p2_x;
	double p_sum = 0.0;

#define DX 5e-3
	for (x=0.0; x<=1.0-DX; x+=DX) {
		//p_p2_x = beta_pdf(1+mac2, 1+tac2-mac2, x);
		p_p2_x = cephes_btdtr(1+mac2, 1+tac2-mac2, x+DX) - cephes_btdtr(1+mac2, 1+tac2-mac2, x);

		if (p_p2_x < EPS)
			continue;
		p_p1_x = drift(mac1, tac1, x, t, N);

		p_sum += (DX * p_p1_x) * (DX * p_p2_x);
	}

	return p_sum;
}

#ifdef TEST1 
int
main(int argc, char **argv)
{
	// minor allele count, total allele count, for populations 1 & 2
	unsigned int mac1, tac1, mac2, tac2;
	unsigned int N; // population size
	unsigned int t; // time, in generations
	double p;

	if (argc != 7) {
		fprintf(stderr, "usage: %s mac1 tac1 mac2 tac2 t N\n", argv[0]);
		return 1;
	}

	errno = 0;
	mac1 = strtoul(argv[1], NULL, 0);
	if (errno) {
		fprintf(stderr, "mac1=%s is no good.\n", argv[1]);
		return 2;
	}

	errno = 0;
	tac1 = strtoul(argv[2], NULL, 0);
	if (errno) {
		fprintf(stderr, "tac1=%s is no good.\n", argv[2]);
		return 2;
	}

	errno = 0;
	mac2 = strtoul(argv[3], NULL, 0);
	if (errno) {
		fprintf(stderr, "mac2=%s is no good.\n", argv[3]);
		return 2;
	}

	errno = 0;
	tac2 = strtoul(argv[4], NULL, 0);
	if (errno) {
		fprintf(stderr, "tac2=%s is no good.\n", argv[4]);
		return 2;
	}

	errno = 0;
	t = strtoul(argv[5], NULL, 0);
	if (errno) {
		fprintf(stderr, "t=%s is no good.\n", argv[5]);
		return 2;
	}
	
	errno = 0;
	N = strtoul(argv[6], NULL, 0);
	if (errno) {
		fprintf(stderr, "N=%s is no good.\n", argv[6]);
		return 2;
	}

	p = drift_test(mac1, tac1, mac2, tac2, t, N);
	printf("%lg\n", p);

	return 0;
}
#endif
