#ifndef _CEPHES_H
#define _CEPHES_H

/* Complex numeral.  */
typedef struct {
    double r;
    double i;
} cephes_cmplx;

extern int cephes_airy(double x, double *ai, double *aip, double *bi,
		double *bip);
extern double cephes_bdtrc(int k, int n, double p);
extern double cephes_bdtr(int k, int n, double p);
extern double cephes_bdtri(int k, int n, double y);
extern double cephes_beta(double a, double b);
extern double cephes_lbeta(double a, double b);
extern double cephes_btdtr(double a, double b, double x);
extern double cephes_cbrt(double x);
extern double cephes_chbevl(double x, double P[], int n);
extern double cephes_chdtrc(double df, double x);
extern double cephes_chdtr(double df, double x);
extern double cephes_chdtri(double df, double y);

extern double cephes_dawsn(double xx);
extern void cephes_eigens(double A[], double RR[], double E[], int N);
extern double cephes_ellie(double phi, double m);
extern double cephes_ellik(double phi, double m);
extern double cephes_ellpe(double x);
extern int cephes_ellpj(double u, double m, double *sn, double *cn, double *dn,
		 double *ph);
extern double cephes_ellpk(double x);
extern double cephes_exp10(double x);
extern double cephes_exp2(double x);
extern double cephes_expn(int n, double x);
extern double cephes_fdtrc(double a, double b, double x);
extern double cephes_fdtr(double a, double b, double x);
extern double cephes_fdtri(double a, double b, double y);

extern int cephes_fresnl(double xxa, double *ssa, double *cca);
extern double cephes_Gamma(double x);
extern double cephes_lgam(double x);
extern double cephes_gdtr(double a, double b, double x);
extern double cephes_gdtrc(double a, double b, double x);
extern int cephes_gels(double A[], double R[], int M, double EPS, double AUX[]);
extern double cephes_hyp2f1(double a, double b, double c, double x);
extern double cephes_hyperg(double a, double b, double x);
extern double cephes_hyp2f0(double a, double b, double x, int type, double *err);
extern double cephes_i0(double x);
extern double cephes_i0e(double x);
extern double cephes_i1(double x);
extern double cephes_i1e(double x);
extern double cephes_igamc(double a, double x);
extern double cephes_igam(double a, double x);
extern double cephes_igami(double a, double y0);
extern double cephes_incbet(double aa, double bb, double xx);
extern double cephes_incbi(double aa, double bb, double yy0);
extern double cephes_iv(double v, double x);
extern double cephes_jv(double n, double x);
extern double cephes_k0(double x);
extern double cephes_k0e(double x);
extern double cephes_k1(double x);
extern double cephes_k1e(double x);
extern double cephes_kn(int nn, double x);

extern void cephes_mmmpy(int r, int c, double *A, double *B, double *Y);
extern int cephes_mtherr(char *name, int code);
extern void cephes_mtransp(int n, double *A, double *T);
extern void cephes_mvmpy(int r, int c, double *A, double *V, double *Y);
extern double cephes_nbdtrc(int k, int n, double p);
extern double cephes_nbdtr(int k, int n, double p);
extern double cephes_nbdtri(int k, int n, double p);
extern double cephes_ndtr(double a);
extern double cephes_erfc(double a);
extern double cephes_erf(double x);
extern double cephes_ndtri(double y0);
extern double cephes_pdtrc(int k, double m);
extern double cephes_pdtr(int k, double m);
extern double cephes_pdtri(int k, double y);
extern double cephes_powi(double x, int nn);
extern double cephes_psi(double x);
extern double cephes_rgamma(double x);
extern double cephes_round(double x);
extern int cephes_sprec(void);
extern int cephes_dprec(void);
extern int cephes_ldprec(void);
extern int cephes_shichi(double x, double *si, double *ci);
extern int cephes_sici(double x, double *si, double *ci);
extern double cephes_simpsn(double f[], double delta);
extern int cephes_simq(double A[], double B[], double X[], int n, int flag,
		int IPS[]);
extern double cephes_radian(double d, double m, double s);

extern void cephes_sincos ( double x, double *s, double *c, int flg);
extern double cephes_sindg(double x);
extern double cephes_cosdg(double x);
extern double cephes_spence(double x);
extern double cephes_stdtr(int k, double t);
extern double cephes_stdtri(int k, double p);
extern double cephes_onef2(double a, double b, double c, double x, double *err);
extern double cephes_threef0(double a, double b, double c, double x, double *err);
extern double cephes_struve(double v, double x);
extern double cephes_tandg(double x);
extern double cephes_cotdg(double x);
extern double cephes_cosm1(double x);
extern double cephes_zeta(double x, double q);
extern double cephes_zetac(double x);

/* polyn.c */
extern void cephes_polini(int maxdeg);
extern void cephes_polprt(double a[], int na, int d);
extern void cephes_polclr(double *a, int n);
extern void cephes_polmov(double *a, int na, double *b);
extern void cephes_polmul(double a[], int na, double b[], int nb, double c[]);
extern void cephes_poladd(double a[], int na, double b[], int nb, double c[]);
extern void cephes_polsub(double a[], int na, double b[], int nb, double c[]);
extern int cephes_poldiv(double a[], int na, double b[], int nb, double c[]);
extern void cephes_polsbt(double a[], int na, double b[], int nb, double c[]);
extern double cephes_poleva(double a[], int na, double x);

/* polmisc.c */
extern void cephes_polatn(double num[], double den[], double ans[], int nn);
extern void cephes_polsqt(double pol[], double ans[], int nn);
extern void cephes_polsin(double x[], double y[], int nn);
extern void cephes_polcos(double x[], double y[], int nn);

/* polrt.c */
int cephes_polrt(double[], double[], int, cephes_cmplx[]);

double cephes_yv(double v, double x);
#endif
