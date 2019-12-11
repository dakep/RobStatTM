/* Copied from robustbase.h */

#include <R.h>
#include <Rinternals.h>
#include <complex.h>

/* For internationalized messages */
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("Matrix", String)
#else
#define _(String) (String)
#define dngettext(pkg, String, StringP, N) (N > 1 ? StringP : String)
#endif

double complex erfz(double complex z);
double complex erfi(double complex z);
SEXP R_erfi(SEXP x);

SEXP R_rho_inf(SEXP cc, SEXP ipsi);

void R_lmrob_S(double *X, double *y, int *n, int *P,
	       int *nRes, double *scale, double *beta_s,
	       double *C, int *iipsi, double *bb,
	       int *best_r, int *Groups, int *N_group,
	       int *K_s, int *max_k, int *max_it_scale,
	       double *rel_tol, double *inv_tol,
               //     ^^^^^^^^^ = refine.tol in R
	       int* converged, int *trace_lev, int *mts, int *ss, int *cutoff);

void R_lmrob_MM(double *X, double *y, int *n, int *P,
		double *beta_initial, double *scale,
		double *beta_m, double *resid,
		int *max_it,
		double *rho_c, int *ipsi, double *loss, double *rel_tol,
		int *converged, int *trace_lev, int *mts, int *ss);

void R_subsample(const double *x, const double *y, int *n, int *m,
		 double *beta, int *ind_space, int *idc, int *idr,
		 double *lu, double *v, int *p,
		 double *_Dr, double *_Dc, int *_rowequ, int *_colequ,
		 int *status, int *sample, int *mts, int *ss, double *tol_inv,
		 int *solve);

SEXP R_psifun(SEXP x_, SEXP c_, SEXP ipsi_, SEXP deriv_);
SEXP R_chifun(SEXP x_, SEXP c_, SEXP ipsi_, SEXP deriv_);
SEXP R_wgtfun(SEXP x_, SEXP c_, SEXP ipsi_);


void R_find_D_scale(double *rr, double *kkappa, double *ttau, int *llength,
		    double *sscale, double *cc, int *iipsi, int *ttype, double *rel_tol,
		    int *max_k, int *converged);

void R_calc_fitted(double *XX, double *bbeta, double *RR, int *nn, int *pp, int *nnrep,
		   int *nnproc, int *nnerr);
