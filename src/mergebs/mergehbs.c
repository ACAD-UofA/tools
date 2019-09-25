/*
 * Merge hairpin ligated or overlapping bisulfite treated reads.
 *
 * Copyright (c) 2017 Graham Gower <graham.gower@gmail.com>
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
#include <string.h>
#include <stdint.h>
#include <getopt.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <zlib.h>
#include <assert.h>
#include "fold.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);

typedef double prob_t;
#define PROB_T_MAX DBL_MAX
#define POW(x,y) pow(x,y)
#define LOG(x) log(x)
#define LOG1P(x) log1p(x)
#define ERF(x) erf(x)
#define EXP(x) exp(x)

#define LL_THRES LOG(20.0)

#define PHRED_MAX 64
//#define PHRED_SHIFT (LOG(PHRED_MAX)/LOG(2.0))
#define PHRED_SHIFT 6

#define READLEN_MAX 1024

#define MOL_UNKNOWN	0x0
#define MOL_HP		0x1
#define MOL_YY		0x2
#define MOL_LONG	0x4 // molecule is longer than read length(s)

typedef struct {
	char *fn1, *fn2;
	char *oprefix; // output prefix
	char *hairpin, *rhairpin; // hairpin and reverse complement
	size_t hlen; // length of hairpin
	char *a1, *a2; // Y adapter sequences
	size_t a1len, a2len; // length of Y adapters

	int phred_scale_in;
	int phred_scale_out;

	prob_t mu, sigma; // mean and std. of log(read length)

	int output_ll;

	FILE *ll_ofp; // log likelihoods
	FILE *hp_ofp; // merged hairpin reads go here
	FILE *yy_ofp; // merged YY reads go here
	FILE *unmerged_r1_ofp,
	     *unmerged_r2_ofp; // unmerged reads go here
} opt_t;


#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

int fit_lognorm(const uint64_t *hist, int len, uint64_t area, double *mu, double *s);

static prob_t _q2p_cache[PHRED_MAX];
/*
 * Covert PHRED scaled Q value to a probability.
 */
static inline prob_t
q2p(char q)
{
	return _q2p_cache[(int)q];
}

static prob_t _match_x_cache[PHRED_MAX*2];
static prob_t _match_cache[PHRED_MAX*PHRED_MAX*2];
static prob_t _lpdf_lognorm_cache[READLEN_MAX];
static prob_t _lsf_lognorm_cache[READLEN_MAX];

prob_t
lpdf_lognorm(size_t x)
{
	if (x >= READLEN_MAX) {
		fprintf(stderr, "reads too long(%zd)! Need to bump READLEN_MAX(%d)\n", x, READLEN_MAX);
		exit(1);
	}
	return _lpdf_lognorm_cache[x];
}

prob_t
lsf_lognorm(size_t x)
{
	if (x >= READLEN_MAX) {
		fprintf(stderr, "reads too long(%zd)! Need to bump READLEN_MAX(%d)\n", x, READLEN_MAX);
		exit(1);
	}
	return _lsf_lognorm_cache[x];
}

static void
init_caches(opt_t *opt)
{
	int x;

	_lpdf_lognorm_cache[0] = -PROB_T_MAX;
	_lsf_lognorm_cache[0] = 0.0;
	for (x=1; x<READLEN_MAX; x++) {

		if (opt->mu > 0) {
			prob_t a, b, c;

			a = ERF((LOG((prob_t)x+0.5) - opt->mu) / (M_SQRT2*opt->sigma));
			b = ERF((LOG((prob_t)x-0.5) - opt->mu) / (M_SQRT2*opt->sigma));
			c = ERF((LOG((prob_t)x) - opt->mu) / (M_SQRT2*opt->sigma));

			_lpdf_lognorm_cache[x] = LOG(a-b) -M_LN2;
			_lsf_lognorm_cache[x] = LOG1P(-c) -M_LN2;
		} else {
			// Uniform prior
			_lpdf_lognorm_cache[x] = 0;
			_lsf_lognorm_cache[x] = 0;
		}
	}

	// probability of error for random nucleotide
	_q2p_cache[0] = _q2p_cache[1] = 0.75; 

	for (x=2; x<PHRED_MAX; x++)
		_q2p_cache[x] = POW(10, -x/10.0);


	for (x=0; x<PHRED_MAX; x++) {
		// match
		_match_x_cache[x<<1|1] = LOG1P(-q2p(x));
		// mismatch
		_match_x_cache[x<<1|0] = LOG(q2p(x)/3.0);
	}


	int x1, x2; // qual scores
	int m; // c1 matches c2?

	for (x1=0; x1<PHRED_MAX; x1++) {
		for (x2=0; x2<PHRED_MAX; x2++) {
			for (m=0; m<2; m++) {
				prob_t l;
				int k = x1<<(PHRED_SHIFT+1) | x2<<1 | m;

				if (m) {
					// match
					l = (1.0-q2p(x1)) * (1.0-q2p(x2));
					l += 3 * (q2p(x1)/3.0) * (q2p(x2)/3.0);
				} else {
					// mismatch
					l = (1.0-q2p(x1)) * (q2p(x2)/3.0);
					l += (q2p(x1)/3.0) * (1.0-q2p(x2));
					l += 2 * (q2p(x1)/3.0) * (q2p(x2)/3.0);
				}

				_match_cache[k] = LOG(l/4.0);
			}
		}
	}
}


static inline prob_t
ll_match_x_1pos(char q, int match)
{
	return _match_x_cache[q<<1 | match];
}

/*
 * loglik that s matches x
 */
prob_t
ll_match_x(char *s, char *q, size_t len, char *x)
{
	prob_t ll = 0.0;
	int i;

	for (i=0; i<len; i++) {
		if (x[i] == 'N')
			ll += LOG(1.0/4.0);
		else
			ll += ll_match_x_1pos(q[i], s[i]==x[i]);
	}

	return ll;
}

static inline prob_t
ll_match_1pos(char q1, char q2, int match)
{
	return _match_cache[q1<<(PHRED_SHIFT+1) | q2<<1 | match];
}

/*
 * loglik that s1 matches s2
 */
prob_t
ll_match(char *s1, char *q1, char *s2, char *q2, size_t len, prob_t *vout)
{
	prob_t x, ll = 0.0;
	int i;

	for (i=0; i<len; i++) {
		x = ll_match_1pos(q1[i], q2[i], s1[i]==s2[i]);
		ll += x;
		if (vout)
			vout[i] = x;
	}

	return ll;
}

/*
 * loglik that s1 matches s2, allowing mismatches due to bisulfite conversion
 * in a hairpin-ligated molecule
 */
prob_t
ll_match_hbs(char *s1, char *q1, char *s2, char *q2, size_t len, prob_t *vout)
{
	prob_t ll = 0.0;
	int i;
	char nt;

	for (i=0; i<len; i++) {
		char c1 = s1[i];
		char c2 = s2[i];
		prob_t e1 = q2p(q1[i]);
		prob_t e2 = q2p(q2[i]);
		prob_t x, l, l_n = 0;

		// P(c1,c2|'A',HP)
		nt = 'A';
		if (c1 == nt)
			l = 1.0-e1;
		else
			l = e1/3.0;
		if (c2 == nt)
			l *= 1.0-e2;
		else
			l *= e2/3.0;
		l_n += l;

		// P(c1,c2|'C',HP)
		nt = 'C';
		if (c1 == nt || c1 == 'T')
			l = (1.0-e1)/2.0;
		else
			l = e1/2.0;
		if (c2 == nt)
			l *= 1.0-e2;
		else
			l *= e2/3.0;
		l_n += l;

		// P(c1,c2|'G',HP)
		nt = 'G';
		if (c1 == nt)
			l = 1.0-e1;
		else
			l = e1/3.0;
		if (c2 == nt || c2 == 'A')
			l *= (1.0-e2)/2.0;
		else
			l *= e2/2.0;
		l_n += l;

		// P(c1,c2|'T',HP)
		nt = 'T';
		if (c1 == nt)
			l = 1.0-e1;
		else
			l = e1/3.0;
		if (c2 == nt)
			l *= 1.0-e2;
		else
			l *= e2/3.0;
		l_n += l;

		// P(c1,c2|HP) = sum_{nt=[A,C,G,T]} P(c1,c2|nt,HP) * P(nt|HP)
		x = LOG(l_n / 4.0);
		ll += x;
		if (vout)
			vout[i] = x;
	}

	return ll;
}

prob_t
sum(prob_t *v, size_t len)
{
	int i;
	prob_t x = 0.0;
	for (i=0; i<len; i++)
		x += v[i];
	return x;
}

/*
 * Y-hairpin ligated molecule.
 * Calculate log[P(data|LEN=i,HP)*P(LEN=i|HP)].
 */
int
ll_given_hp(const opt_t *opt,
		char *s1, char *q1,
		char *s2, char *q2, size_t len,
		prob_t *llvec, prob_t *llsum,
		int *i_max, int *i_2nd)
{
	int ret;
	int i;
	prob_t *_ll_match_hbs;
	prob_t sum_ll_match_hbs;
	prob_t psum;
	prob_t ll_max, ll_2nd;

	ll_max = ll_2nd = -PROB_T_MAX;
	*i_max = *i_2nd = -1;
	llvec[0] = 0;
	psum = 0;

	_ll_match_hbs = calloc(len, sizeof(prob_t));
	if (_ll_match_hbs == NULL) {
		perror("calloc");
		ret = -1;
		goto err0;
	}

	// cache some results
	ll_match_hbs(s1, q1, s2, q2, len, _ll_match_hbs);
	sum_ll_match_hbs = 0;

	llvec[0] = 0;
	psum = 0;

	// hairpin at pos i
	for (i=1; i<=len; i++) {
		prob_t ll;

		sum_ll_match_hbs += _ll_match_hbs[i-1];

		// target sequence upstream of hairpin
		ll = sum_ll_match_hbs;

		// hairpin
		ll += ll_match_x(s1+i, q1+i, min(len-i, opt->hlen), opt->hairpin);
		ll += ll_match_x(s2+i, q2+i, min(len-i, opt->hlen), opt->rhairpin);

		// target sequence downstream of hairpin
		if (i+opt->hlen < len)
			ll += sum(_ll_match_hbs+i+opt->hlen, min(len-(i+opt->hlen), i));

		// Y adapter after target sequence
		if (i+opt->hlen+i < len) {
			ll += ll_match_x(s1+i+opt->hlen+i, q1+i+opt->hlen+i,
					min(len-(i+opt->hlen+i), opt->a1len),
					opt->a1)
				+ ll_match_x(s2+i+opt->hlen+i, q2+i+opt->hlen+i,
					min(len-(i+opt->hlen+i), opt->a2len),
					opt->a2);
		}

		// unknown sequence after the Y adapters
		if (i+opt->hlen+i+opt->a1len < len)
			ll += (len-(i+opt->hlen+i+opt->a1len)) * LOG(1.0/4.0);
		if (i+opt->hlen+i+opt->a2len < len)
			ll += (len-(i+opt->hlen+i+opt->a2len)) * LOG(1.0/4.0);

		if (i<len)
			// log[P(LEN=i|HP)]
			ll += lpdf_lognorm(i);
		else
			// log[P(LEN>=i|HP)]
			ll += lsf_lognorm(i);


		if (ll > ll_max) {
			ll_2nd = ll_max;
			ll_max = ll;
			*i_2nd = *i_max;
			*i_max = i;
		} else if (ll > ll_2nd) {
			ll_2nd = ll;
			*i_2nd = i;
		}

		llvec[i] = ll;
		psum += EXP(ll);
	}

	*llsum = LOG(psum);
	ret = 0;

	free(_ll_match_hbs);
err0:
	return ret;
}

/*
 * Molecule with Y adapter ligated on both ends.
 * Calculate log[P(data|LEN=i,YY)*P(LEN=i|YY)].
 */
int
ll_given_yy(const opt_t *opt,
		char *s1, char *q1,
		char *s2, char *q2, size_t len,
		prob_t *llvec, prob_t *llsum,
		int *i_max, int *i_2nd)
{
	int ret;
	int i;
	char *rs2, *rq2;
	prob_t psum;
	prob_t ll_max, ll_2nd;

	ll_max = ll_2nd = -PROB_T_MAX;
	*i_max = *i_2nd = -1;
	llvec[0] = 0;
	psum = 0;

	rs2 = malloc(len);
	if (rs2 == NULL) {
		perror("malloc");
		ret = -1;
		goto err0;
	}
	memcpy(rs2, s2, len);
	revcomp(rs2, len);

	rq2 = malloc(len);
	if (rq2 == NULL) {
		perror("malloc");
		ret = -2;
		goto err1;
	}
	memcpy(rq2, q2, len);
	reverse(rq2, len);

	/*
	 * Y-Y molecule, reads completely contain the molecule,
	 * Y adapter at pos i, read 1 and read 2 have 'i' bases of sequence identity
	 */
	for (i=1; i<len; i++) {
		prob_t ll;

		// overlapping sequence
		ll = ll_match(s1, q1, rs2+(len-i), rq2+(len-i), i, NULL);

		// adapters
		ll += ll_match_x(s1+i, q1+i, min(len-i, opt->a1len), opt->a1);
		ll += ll_match_x(s2+i, q2+i, min(len-i, opt->a2len), opt->a2);

		// non-matching sequence after the Y adapters
		if (i+opt->a1len < len)
			ll += (len-(i+opt->a1len)) * LOG(1.0/4.0);
		if (i+opt->a2len < len)
			ll += (len-(i+opt->a2len)) * LOG(1.0/4.0);

		// log[P(LEN=i|YY)]
		ll += lpdf_lognorm(i);

		if (ll > ll_max) {
			ll_2nd = ll_max;
			ll_max = ll;
			*i_2nd = *i_max;
			*i_max = i;
		} else if (ll > ll_2nd) {
			ll_2nd = ll;
			*i_2nd = i;
		}

		llvec[i] = ll;
		psum += EXP(ll);
	}

	/*
	 * Y-Y molecule, reads do not completely contain the molecule and no
	 * adapter is observed, read 1 and read 2 have 'len-i' bases of sequence overlap
	 */
	for (i=0; i<=len; i++) {
		prob_t ll;

		// non-overlapping sequence
		ll = 2 * i * LOG(1.0/4.0);

		// overlapping sequence
		ll += ll_match(s1+i, q1+i, rs2, rq2, len-i, NULL);

		if (i<len)
			// log[P(LEN=len+i|YY)]
			ll += lpdf_lognorm(len+i);
		else
			// log[P(LEN>=len+i|YY)]
			ll += lsf_lognorm(len+i);

		if (ll > ll_max) {
			ll_2nd = ll_max;
			ll_max = ll;
			*i_2nd = *i_max;
			*i_max = len+i;
		} else if (ll > ll_2nd) {
			ll_2nd = ll;
			*i_2nd = len+i;
		}

		llvec[len+i] = ll;
		psum += EXP(ll);
	}

	*llsum = LOG(psum);
	ret = 0;

	free(rq2);
err1:
	free(rs2);
err0:
	return ret;
}

/*
 * Maximum posterior classifier of a read as originating from either a
 * Y-Hairpin ligated or a Y-Y ligated molecule.
 */
int
classify(const opt_t *opt, char *s1, char *q1, char *s2, char *q2, size_t len, int *mol_len)
{
	int ret;
	int class;
	int i_max; // length at maximum posterior
	int i_hp_max, i_hp_2nd, i_yy_max, i_yy_2nd;
	prob_t *ll_hp, *ll_yy; // log[P(data|i,HP)], log[P(data|i,YY)]
	prob_t ll_hp_sum, ll_yy_sum; // log[P(data|HP)], log[P(data|YY)]
	prob_t ll_max, ll_2nd; // log[maximum posterior]

	ll_hp = calloc(len+1, sizeof(prob_t));
	if (ll_hp == NULL) {
		perror("calloc");
		ret = -1;
		goto err0;
	}
	if (ll_given_hp(opt, s1, q1, s2, q2, len, ll_hp, &ll_hp_sum, &i_hp_max, &i_hp_2nd) < 0) {
		ret = -2;
		goto err1;
	}

	assert(i_hp_max > 0);
	assert(i_hp_2nd > 0);

	ll_yy = calloc(2*len+1, sizeof(prob_t));
	if (ll_yy == NULL) {
		perror("calloc");
		ret = -3;
		goto err1;
	}

	if (ll_given_yy(opt, s1, q1, s2, q2, len, ll_yy, &ll_yy_sum, &i_yy_max, &i_yy_2nd) < 0) {
		ret = -4;
		goto err2;
	}

	assert(i_yy_max > 0);
	assert(i_yy_2nd > 0);

	ll_max = ll_2nd = -PROB_T_MAX;
	i_max = -1;
	class = MOL_UNKNOWN;

	if (ll_hp_sum > ll_yy_sum && ll_hp_sum -ll_yy_sum > LL_THRES) {
		class = MOL_HP;
		ll_max = ll_hp[i_hp_max];
		i_max = i_hp_max;
		if (i_max == len)
			class |= MOL_LONG;
		ll_2nd = ll_hp[i_hp_2nd];
	} else if (ll_yy_sum > ll_hp_sum && ll_yy_sum -ll_hp_sum > LL_THRES) {
		class = MOL_YY;
		ll_max = ll_yy[i_yy_max];
		i_max = i_yy_max;
		if (i_max == 2*len)
			class |= MOL_LONG;
		ll_2nd = ll_yy[i_yy_2nd];
	}

	if (ll_max - ll_2nd < LL_THRES) {
		// most probable isn't much better than the 2nd most probable
		i_max = -1;
		class = MOL_UNKNOWN;
	}

	if (opt->ll_ofp) {
		int i;
		fprintf(opt->ll_ofp, "HP\tsum\t%lf\n", ll_hp_sum);
		fprintf(opt->ll_ofp, "YY\tsum\t%lf\n", ll_yy_sum);
		for (i=1; i<=len; i++)
			fprintf(opt->ll_ofp, "HP\t%d\t%lf\n", i, ll_hp[i]-ll_hp_sum);
		for (i=1; i<=2*len; i++)
			fprintf(opt->ll_ofp, "YY\t%d\t%lf\n", i, ll_yy[i]-ll_yy_sum);
	}

	ret = class;
	*mol_len = i_max;

err2:
	free(ll_yy);
err1:
	free(ll_hp);
err0:
	return ret;
}

/*
 * Convert ASCII value to phred scaled Q score.
 */
void
clean_quals(const opt_t *opt, const char *s, char *q, size_t len)
{
	int i;

	for (i=0; i<len; i++) {
		char *qi = q+i;

		if (s[i] == 'N') {
			*qi = 0;
			continue;
		}

		*qi -= opt->phred_scale_in;
		if (*q < 2)
			*q = 2;
		if (*q >= PHRED_MAX)
			*q = PHRED_MAX-1;
	}
}

/*
 * Convert phred scaled Q score to ASCII and print to file.
 */
void
fput_qual(FILE *fp, int phred_scale_out, const char *q, size_t len)
{
	int i;
	for (i=0; i<len; i++)
		fputc(q[i]+phred_scale_out, fp);
}

void
write_unmerged(const opt_t *opt, const char *name,
		const char *s1, const char *q1,
		const char *s2, const char *q2,
		size_t len)
{
	fprintf(opt->unmerged_r1_ofp, "@%s\n%s\n+\n", name, s1);
	fput_qual(opt->unmerged_r1_ofp, opt->phred_scale_out, q1, len);
	fputc('\n', opt->unmerged_r1_ofp);

	fprintf(opt->unmerged_r2_ofp, "@%s\n%s\n+\n", name, s2);
	fput_qual(opt->unmerged_r2_ofp, opt->phred_scale_out, q2, len);
	fputc('\n', opt->unmerged_r2_ofp);
}

int
merge_yy(const opt_t *opt, const char *name,
		const char *s1, const char *q1,
		const char *s2, const char *q2,
		size_t len, int mol_len)
{
	int ret;
	int mm;
	int overlap;
	char *rs2, *rq2;
	char *s_out, *q_out;

	rs2 = malloc(len);
	if (rs2 == NULL) {
		perror("malloc");
		ret = -1;
		goto err0;
	}
	memcpy(rs2, s2, len);
	revcomp(rs2, len);

	rq2 = malloc(len);
	if (rq2 == NULL) {
		perror("malloc");
		ret = -2;
		goto err1;
	}
	memcpy(rq2, q2, len);
	reverse(rq2, len);

	s_out = malloc(mol_len+1);
	if (s_out == NULL) {
		perror("malloc");
		ret = -3;
		goto err2;
	}

	q_out = malloc(mol_len+1);
	if (q_out == NULL) {
		perror("malloc");
		ret = -4;
		goto err3;
	}

	if (mol_len <= len) {
		// complete overlap
		overlap = mol_len;
		mm = match2(s1, q1, overlap, rs2+(len-overlap), rq2+(len-overlap), overlap,
				s_out, q_out,
				0, 0, 0);
	} else if (mol_len < 2*len) {
		// partial overlap
		int i;
		overlap = 2*len - mol_len;
		for (i=0; i<len-overlap; i++) {
			s_out[i] = s1[i];
			q_out[i] = q1[i];
		}
		mm = match2(s1+(len-overlap), q1+(len-overlap), overlap, rs2, rq2, overlap,
				s_out+(len-overlap), q_out+(len-overlap),
				0, 0, 0);
		for (i=0; i<len-overlap; i++) {
			s_out[len+i] = rs2[overlap+i];
			q_out[len+i] = rq2[overlap+i];
		}

	} else
		overlap = 0;

	s_out[mol_len] = q_out[mol_len] = '\0';

	if (overlap && mm <= maxdiff(overlap, AVG_ERR, MAXDIFF_THRES)) {
		// merge success
		fprintf(opt->yy_ofp, "@%s\n%s\n+\n", name, s_out);
		fput_qual(opt->yy_ofp, opt->phred_scale_out, q_out, mol_len);
		fputc('\n', opt->yy_ofp);
	} else {
		// failure, too many mismatches
		write_unmerged(opt, name, s1, q1, s2, q2, len);
		ret = 1;
	}

	ret = 0;
//err4:
	free(q_out);
err3:
	free(s_out);
err2:
	free(rq2);
err1:
	free(rs2);
err0:
	return ret;
}

int
merge_hairpin(const opt_t *opt,
		const char *name,
		const char *s1, const char *q1,
		const char *s2, const char *q2,
		size_t len, int mol_len)
{
	int ret;
	int mm;
	char *s_out, *q_out;

	s_out = malloc(mol_len+1);
	if (s_out == NULL) {
		perror("malloc");
		ret = -1;
		goto err0;
	}

	q_out = malloc(mol_len+1);
	if (q_out == NULL) {
		perror("malloc");
		ret = -2;
		goto err1;
	}

	if (mol_len+opt->hlen < len) {
		const char *s3 = s1 + mol_len + opt->hlen;
		const char *q3 = q1 + mol_len + opt->hlen;
		const char *s4 = s2 + mol_len + opt->hlen;
		const char *q4 = q2 + mol_len + opt->hlen;
		int len3 = 2*mol_len+opt->hlen > len ? len - mol_len - opt->hlen : mol_len;
		int len4 = 2*mol_len+opt->hlen > len ? len - mol_len - opt->hlen : mol_len;

		mm = match4(s1, q1, mol_len, s2, q2, mol_len,
				s3, q3, len3, s4, q4, len4,
				s_out, q_out,
				0, 0);
	} else {
		mm = match2(s1, q1, mol_len, s2, q2, mol_len,
				s_out, q_out,
				1, 0, 0);
	}

	s_out[mol_len] = q_out[mol_len] = '\0';

	if (mm <= maxdiff(mol_len*2, AVG_ERR, MAXDIFF_THRES)) {
		// merge success
		fprintf(opt->hp_ofp, "@%s XF:Z:%s|%s|", name, s1, s2);
		fput_qual(opt->hp_ofp, opt->phred_scale_out, q1, len);
		fputc('|', opt->hp_ofp);
		fput_qual(opt->hp_ofp, opt->phred_scale_out, q2, len);
		fprintf(opt->hp_ofp, "\n%s\n+\n", s_out);
		fput_qual(opt->hp_ofp, opt->phred_scale_out, q_out, mol_len);
		fputc('\n', opt->hp_ofp);
	} else {
		// failure, too many mismatches
		write_unmerged(opt, name, s1, q1, s2, q2, len);
		ret = 1;
	}


	ret = 0;
	free(q_out);
err1:
	free(s_out);
err0:
	return ret;
}

/*
 */
int
mergehbs(const opt_t *opt)
{
	gzFile fp1, fp2;
	kseq_t *seq1, *seq2;
	int len1, len2;
	int ret, mret;
	int class;
	int mol_len;

	// count some stuff
	struct {
		uint64_t total_reads;
		uint64_t hp;
		uint64_t yy;
		uint64_t hp_longmol;
		uint64_t yy_longmol;
		uint64_t hp_complete;
		uint64_t yy_complete;
		uint64_t hp_merged;
		uint64_t yy_merged;
		uint64_t unknown;
	} counts = {0,};

	uint64_t *hp_hist = NULL, *yy_hist = NULL;
	int hp_hist_len = 0, yy_hist_len = 0;

	fp1 = gzopen(opt->fn1, "r");
	if (fp1 == NULL) {
		fprintf(stderr, "%s: %s\n", opt->fn1, strerror(errno));
		ret = -1;
		goto err0;
	}

	fp2 = gzopen(opt->fn2, "r");
	if (fp2 == NULL) {
		fprintf(stderr, "%s: %s\n", opt->fn2, strerror(errno));
		ret = -1;
		goto err1;
	}

	seq1 = kseq_init(fp1);
	seq2 = kseq_init(fp2);

	for (;;) {
		len1 = kseq_read(seq1);
		len2 = kseq_read(seq2);
		if (len1 < 0 || len2 < 0)
			break;

		if (seq1->qual.l == 0) {
			fprintf(stderr, "%s: qual scores required.\n", opt->fn1);
			ret = -2;
			goto err3;
		}

		if (seq2->qual.l == 0) {
			fprintf(stderr, "%s: qual scores required.\n", opt->fn2);
			ret = -2;
			goto err3;
		}

		if ((seq1->name.s[seq1->name.l-2] == '/' && seq1->name.s[seq1->name.l-1] == '1' &&
			seq2->name.s[seq2->name.l-2] == '/' && seq2->name.s[seq2->name.l-1] == '2')
		  || (seq1->name.s[seq1->name.l-2] == '.' && seq1->name.s[seq1->name.l-1] == '1' &&
			seq2->name.s[seq2->name.l-2] == '.' && seq2->name.s[seq2->name.l-1] == '2')) {
			seq1->name.l -= 2;
			seq1->name.s[seq1->name.l] = 0;
			seq2->name.l -= 2;
			seq2->name.s[seq2->name.l] = 0;
		}

		if (strcmp(seq1->name.s, seq2->name.s)) {
			fprintf(stderr, "R1 and R2 sequence names don't match:\n%s\n%s\n",
					seq1->name.s, seq2->name.s);
			ret = -3;
			goto err3;
		}

		if (len1 != len2) {
			fprintf(stderr, "R1 and R2 have differing lengths (%d and %d)\n",
					len1, len2);
			ret = -4;
			goto err3;
		}

		counts.total_reads++;

		if (len1 < opt->hlen || len2 < opt->hlen)
			continue;

		if (hp_hist == NULL) {
			hp_hist_len = len1;
			hp_hist = calloc(hp_hist_len, sizeof(*hp_hist));
			if (hp_hist == NULL) {
				perror("malloc");
				ret = -5;
				goto err3;
			}

			yy_hist_len = len1+len2;
			yy_hist = calloc(yy_hist_len, sizeof(*yy_hist));
			if (yy_hist == NULL) {
				perror("malloc");
				ret = -6;
				goto err3;
			}
		} else if (hp_hist_len != len1 || yy_hist_len != len1+len2) {
			fprintf(stderr, "Reads have different length to previous reads: %s\n",
					seq1->name.s);
			ret = -7;
			goto err3;
		}

		clean_quals(opt, seq1->seq.s, seq1->qual.s, seq1->qual.l);
		clean_quals(opt, seq2->seq.s, seq2->qual.s, seq2->qual.l);

		if (opt->ll_ofp)
			fprintf(opt->ll_ofp, "@%s\n", seq1->name.s);

		class = classify(opt,
				seq1->seq.s, seq1->qual.s,
				seq2->seq.s, seq2->qual.s,
				len1,
				&mol_len);

		switch (class & (MOL_HP|MOL_YY)) {
			case MOL_HP:
				counts.hp++;
				if (class & MOL_LONG)
					counts.hp_longmol++;
				if (mol_len+opt->hlen <= len1)
					counts.hp_complete++;

				mret = merge_hairpin(opt, seq1->name.s,
							seq1->seq.s, seq1->qual.s,
							seq2->seq.s, seq2->qual.s,
							len1, mol_len);
				if (mret < 0) {
					ret = -8;
					goto err3;
				}

				if (mret == 0)
					counts.hp_merged++;

				if (mret == 0 && mol_len < hp_hist_len)
					hp_hist[mol_len]++;

				break;
			case MOL_YY:
				counts.yy++;
				if (class & MOL_LONG)
					counts.yy_longmol++;
				if (mol_len+max(opt->a1len,opt->a2len) <= len1)
					counts.yy_complete++;

				mret = merge_yy(opt, seq1->name.s,
						seq1->seq.s, seq1->qual.s,
						seq2->seq.s, seq2->qual.s,
						len1, mol_len);
				if (mret < 0) {
					ret = -9;
					goto err3;
				}

				if (mret == 0)
					counts.yy_merged++;

				if (mret == 0 && mol_len < yy_hist_len)
					yy_hist[mol_len]++;

				break;
			case MOL_UNKNOWN:
				counts.unknown++;
				write_unmerged(opt, seq1->name.s,
						seq1->seq.s, seq1->qual.s,
						seq2->seq.s, seq2->qual.s,
						len1);
				break;
		}
	}

	if (len1 != -1) {
		fprintf(stderr, "%s: unexpected end of file.\n", opt->fn1);
		ret = -10;
		goto err3;
	}

	if (len2 != -1) {
		fprintf(stderr, "%s: unexpected end of file.\n", opt->fn2);
		ret = -11;
		goto err3;
	}

	double hp_mu, hp_s;
	double yy_mu, yy_s;

	fit_lognorm(hp_hist, hp_hist_len, counts.hp_merged, &hp_mu, &hp_s);
	fit_lognorm(yy_hist, yy_hist_len, counts.yy_merged, &yy_mu, &yy_s);

	fprintf(stderr, "counts.total_reads=%jd\n", counts.total_reads);
	fprintf(stderr, "counts.hp=%jd\n", counts.hp);
	fprintf(stderr, "counts.yy=%jd\n", counts.yy);
	fprintf(stderr, "counts.hp_complete=%jd\n", counts.hp_complete);
	fprintf(stderr, "counts.yy_complete=%jd\n", counts.yy_complete);
	fprintf(stderr, "counts.hp_longmol=%jd\n", counts.hp_longmol);
	fprintf(stderr, "counts.yy_longmol=%jd\n", counts.yy_longmol);
	fprintf(stderr, "counts.hp_merged=%jd\n", counts.hp_merged);
	fprintf(stderr, "counts.yy_merged=%jd\n", counts.yy_merged);
	fprintf(stderr, "counts.unknown=%jd\n", counts.unknown);
	fprintf(stderr, "HP lengths ~LogNorm(mu=%lf,sigma=%lf)\n", hp_mu, hp_s);
	fprintf(stderr, "YY lengths ~LogNorm(mu=%lf,sigma=%lf)\n", yy_mu, yy_s);

	/*
	int i;
	for (i=0; i<hp_hist_len; i++)
		fprintf(stderr, "HIST_HP\t%d\t%jd\n", hp_hist[i]);
	for (i=0; i<yy_hist_len; i++)
		fprintf(stderr, "HIST_YY\t%d\t%jd\n", yy_hist[i]);
	*/

	ret = 0;

err3:
	if (yy_hist)
		free(yy_hist);
	if (hp_hist)
		free(hp_hist);
	kseq_destroy(seq2);
	kseq_destroy(seq1);
//err2:
	gzclose(fp2);
err1:
	gzclose(fp1);
err0:
	return ret;
}

// open output files
int
open_ofp(opt_t *opt)
{
	int ret;
	char *tmpfn;

	tmpfn = malloc(strlen(opt->oprefix) + 64);
	if (tmpfn == NULL) {
		perror("malloc");
		ret = -1;
		goto err0;
	}

	sprintf(tmpfn, "%s.hairpin.fq", opt->oprefix);
	if ((opt->hp_ofp = fopen(tmpfn, "w")) == NULL) {
		fprintf(stderr, "%s: %s\n", tmpfn, strerror(errno));
		ret = -2;
		goto err1;
	}

	sprintf(tmpfn, "%s.yy.fq", opt->oprefix);
	if ((opt->yy_ofp = fopen(tmpfn, "w")) == NULL) {
		fprintf(stderr, "%s: %s\n", tmpfn, strerror(errno));
		ret = -3;
		goto err2;
	}

	sprintf(tmpfn, "%s.unmerged_r1.fq", opt->oprefix);
	if ((opt->unmerged_r1_ofp = fopen(tmpfn, "w")) == NULL) {
		fprintf(stderr, "%s: %s\n", tmpfn, strerror(errno));
		ret = -4;
		goto err3;
	}

	sprintf(tmpfn, "%s.unmerged_r2.fq", opt->oprefix);
	if ((opt->unmerged_r2_ofp = fopen(tmpfn, "w")) == NULL) {
		fprintf(stderr, "%s: %s\n", tmpfn, strerror(errno));
		ret = -5;
		goto err4;
	}

	if (opt->output_ll) {
		sprintf(tmpfn, "%s.ll.txt", opt->oprefix);
		if ((opt->ll_ofp = fopen(tmpfn, "w")) == NULL) {
			fprintf(stderr, "%s: %s\n", tmpfn, strerror(errno));
			ret = -6;
			goto err5;
		}
	}

	free(tmpfn);
	return 0;

err5:
	fclose(opt->unmerged_r1_ofp);
err4:
	fclose(opt->unmerged_r2_ofp);
err3:
	fclose(opt->yy_ofp);
err2:
	fclose(opt->hp_ofp);
err1:
	free(tmpfn);
err0:
	return ret;
}

// close output files
void
close_ofp(opt_t *opt)
{
	if (opt->ll_ofp)
		fclose(opt->ll_ofp);
	fclose(opt->hp_ofp);
	fclose(opt->yy_ofp);
	fclose(opt->unmerged_r1_ofp);
	fclose(opt->unmerged_r2_ofp);
}

void
usage(char *argv0, opt_t *opt)
{
	fprintf(stderr, "mergehbs v1\n");
	fprintf(stderr, "usage: %s [...] r1.fq r2.fq\n", argv0);
	fprintf(stderr, " -h STR         Hairpin sequence [%s]\n", opt->hairpin);
	fprintf(stderr, " -l             Output log-likelihoods [%s]\n", ((char *[]){"no", "yes"})[opt->output_ll]);
	fprintf(stderr, " -o STR         Prefix for output files [%s]\n", opt->oprefix);
	fprintf(stderr, " -s FLOAT       Std.dev. of LogNormal molecule length [%.3f]\n", opt->sigma);
	fprintf(stderr, " -u FLOAT       Mean of LogNormal molecule length [%.3f]\n", opt->mu);
	fprintf(stderr, " -1 STR         Adapter 1 sequence [%s]\n", opt->a1);
	fprintf(stderr, " -2 STR         Adapter 2 sequence [%s]\n", opt->a2);
	exit(1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	int c;
	int i;
	int ret;

	memset(&opt, '\0', sizeof(opt_t));
	opt.oprefix = "mergehbs";
	opt.hairpin = "ACGCCGGCGGCAAGTGAAGCCGCCGGCGT";
	opt.a1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
	opt.a2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";
	opt.phred_scale_in = 33;
	opt.phred_scale_out = 33;

	opt.mu = -1;
	opt.sigma = -1;
	// sane aDNA defaults?
	//opt.mu = 4.0;
	//opt.sigma = 0.25;

	while ((c = getopt(argc, argv, "h:lo:s:u:1:2:")) != -1) {
		switch (c) {
			case 'h':
				opt.hairpin = optarg;
				for (i=0; i<strlen(opt.hairpin); i++)
					opt.hairpin[i] = toupper(opt.hairpin[i]);
				break;
			case 'l':
				opt.output_ll = 1;
				break;
			case 'o':
				opt.oprefix = optarg;
				break;
			case 's':
				opt.sigma = atof(optarg);
				break;
			case 'u':
				opt.mu = atof(optarg);
				break;
			case '1':
				opt.a1 = optarg;
				for (i=0; i<strlen(opt.a1); i++)
					opt.a1[i] = toupper(opt.a1[i]);
				break;
			case '2':
				opt.a2 = optarg;
				for (i=0; i<strlen(opt.a2); i++)
					opt.a2[i] = toupper(opt.a2[i]);
				break;
			default:
				usage(argv[0], &opt);
		}
	}

	if (argc-optind != 2) {
		usage(argv[0], &opt);
	}

	opt.fn1 = argv[optind];
	opt.fn2 = argv[optind+1];

	opt.hlen = strlen(opt.hairpin);
	if (opt.hlen < 5 || opt.hlen > 100) {
		fprintf(stderr, "Error: hairpin is too %s (len=%zd).\n",
				opt.hlen<5?"short":"long", opt.hlen);
		usage(argv[0], &opt);
	}

	opt.rhairpin = strdup(opt.hairpin);
	revcomp(opt.rhairpin, opt.hlen);

	opt.a1len = strlen(opt.a1);
	if (opt.a1len < 5 || opt.a1len > 100) {
		fprintf(stderr, "Error: adapter1 is too %s (len=%zd).\n",
				opt.a1len<5?"short":"long", opt.a1len);
		usage(argv[0], &opt);
	}

	opt.a2len = strlen(opt.a2);
	if (opt.a2len < 5 || opt.a2len > 100) {
		fprintf(stderr, "Error: adapter2 is too %s (len=%zd).\n",
				opt.a2len<5?"short":"long", opt.a2len);
		usage(argv[0], &opt);
	}

	if (open_ofp(&opt) < 0) {
		ret = -2;
		goto err1;
	}

	init_caches(&opt);

	ret = mergehbs(&opt);

	close_ofp(&opt);

err1:
	free(opt.rhairpin);

	return ret;
}
