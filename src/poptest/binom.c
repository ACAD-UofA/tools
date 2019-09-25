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
#include <errno.h>
#include <math.h>
#include <assert.h>


/*
 * ln( (n choose k) = n!/(k!*(n-k)!) )
 */
double
lchoose(double k, double n)
{
	return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}

/*
 * ln( Bin(k, n, p) = (n choose k) * p**k * (1-p)**(n-k) )
 */
double
binom_lpmf(double k, double n, double p)
{
	double l_sum = 0.0;

	assert(k >= 0);
	assert(n >= 0);
	assert(p >= 0.0);
	assert(p <= 1.0);

	/*
	 * edge cases
	 */
	if (n == 0)
		return NAN;
	if (k == 0 && p == 0.0)
		return 0.0;
	if (p == 1.0) {
		if (k == n)
			return 0.0;
		else
			return -INFINITY;
	}
	
	l_sum += lchoose(k, n);
	l_sum += log(p)*k;
	l_sum += log1p(-p)*(n-k);

	return l_sum;
}

/*
 * ln( Integral of the binomial pmf with respect to p )
 * = ln( (n choose k) * B(n+1, k+1) ),
 * where B(x,y) = (x-1)!(y-1)!/(x+y-1)!
 */
double
int_binom_wrt_p(int k, int n)
{
	double l_beta = lgamma(k+1) + lgamma(n+1) - lgamma(k+n+1);
	return lchoose(k, n) + l_beta;
}


#ifdef TEST1
/*
 * Evaluate from command line parameters.
 */
int
main(int argc, char **argv)
{
	double n, k;
	double p;

	if (argc != 4) {
		fprintf(stderr, "usage: %s k n p\n", argv[0]);
		return 1;
	}

	errno = 0;
	k = strtod(argv[1], NULL);
	if (errno) {
		fprintf(stderr, "k=%s is no good.\n", argv[1]);
		return 2;
	}

	errno = 0;
	n = strtod(argv[2], NULL);
	if (errno) {
		fprintf(stderr, "n=%s is no good.\n", argv[2]);
		return 2;
	}

	errno = 0;
	p = strtod(argv[3], NULL);
	if (errno) {
		fprintf(stderr, "p=%s is no good.\n", argv[3]);
		return 2;
	}

	printf("ln(Choose(%lf,%lf)) = %lf\n", k, n, lchoose(k, n));
	printf("ln(B(%lf,%lf,%lf)) = %lf\n", k, n, p, binom_lpmf(k, n, p));
	return 0;
}
#endif

#ifdef TEST2
/*
 * For speed comparisons.
 */
int
main(int argc, char **argv)
{
	int n, k;
	double p, sum=0.0;

	for (n=0; n<=1000; n++) {
		for (k=0; k<=n; k++) {
			for (p=0.01; p<=1.0; p+=0.01)
				sum += binom_lpmf(k, n, p);
		}
	}
	printf("%lf\n", sum);
	return 0;
}
#endif
