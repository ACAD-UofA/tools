/*
 * Species tree
 *        /\
 *       /  \
 *      /    \
 *     /      \
 *    /--------\
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

#include <math.h>

/*
 * n choose k = n!/(k!*(n-k)!)
 * n choose 2 = n!/(2*(n-2)!) = n*(n-1)/2
 */
double inline
choose2(int n)
{
	return 0.5 * (n * (n-1));
}

/*
 * t, array of coalesence times measured backwards in time,
 * 	in units of 2N generations, starting with the tip date.
 * T, time at the top of the branch, in units of 2N generations.
 * m, number of lineages at the bottom of the branch.
 * n, number of lineages at the top of the branch.
 * theta = 4*u*N
 */
double
branch_density(double *t, double T, int m, int n, double *theta)
{
	int i;
	int w = m-n; // number of coalesent events
	double f;
	double _2_theta = 2.0/theta;
	double l2_theta = log(2.0) - log(theta);

	if (n == 1 || T-t[w] == 0)
		// root of the tree
		f = 0.0;
	else
		// no coalescence between the top most node and the top of the tree
		f = - _2_theta * choose2(n) * (T - t[w]);

	for (i=0; i<w; i++)
		f += l2_theta - _2_theta * choose2(m-i) * (t[i+1] - t[i]);

	return f;
}

double
tree1_density(double *t)
{
	double f;

	f = branch_density(t, T, 3, 2, theta)
		+ branch_density(t+1, T, 2, 1, theta);

	return f;
}

double
llik(double *y, double **tt, size_t ntrees, size_t treelen, double theta)
{
	int i;
	double ll = 0.0;

	for (i=0; i<ntrees; i++)
		ll += gene_density(tt[i], treelen, theta);

	return ll;
}
