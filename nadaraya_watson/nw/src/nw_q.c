#include <math.h>
#include <R.h>
#include <Rinternals.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("stats", String)
#else
#define _(String) (String)
#endif

static double dokern(double x, int kern)
{
    if(kern == 1) return(1.0);
    if(kern == 2) return(exp(-0.5*x*x));
    return(0.0); /* -Wall */
}

typedef struct {
    double y;
    double w;
} yw_pair;

static int cmp_y(const void *a, const void *b)
{
    double ya = ((yw_pair *)a)->y;
    double yb = ((yw_pair *)b)->y;
    if (ya < yb) return -1;
    if (ya > yb) return 1;
    return 0;
}

static void BDRkquantile(double *x, double *y, R_xlen_t n,
                         double *xp, double *yp, R_xlen_t np,
                         int kern, double bw, double tau)
{
    if (kern == 1) bw *= 0.5;
    if (kern == 2) bw *= 0.3706506;

    yw_pair *arr = (yw_pair *) R_alloc(n, sizeof(yw_pair));

    for (R_xlen_t j = 0; j < np; j++) {
        double x0 = xp[j];
        double wsum = 0.0;

        /* compute weights */
        for (R_xlen_t i = 0; i < n; i++) {
            double w = dokern(fabs(x[i] - x0) / bw, kern);
            arr[i].y = y[i];
            arr[i].w = w;
            wsum += w;
        }

        if (wsum <= 0) {
            yp[j] = NA_REAL;
            continue;
        }

        /* sort by y */
        qsort(arr, n, sizeof(yw_pair), cmp_y);

        /* cumulative weighted quantile */
        double target = tau * wsum;
        double cum = 0.0;

        for (R_xlen_t i = 0; i < n; i++) {
            cum += arr[i].w;
            if (cum >= target) {
                yp[j] = arr[i].y;
                break;
            }
        }
    }
}

SEXP nw_q(SEXP x, SEXP y, SEXP xp,
          SEXP skrn, SEXP sbw, SEXP stau)
{
    int krn = asInteger(skrn);
    double bw = asReal(sbw);
    double tau = asReal(stau);

    if (tau <= 0 || tau >= 1)
        error("tau must be in (0,1)");

    x  = PROTECT(coerceVector(x, REALSXP));
    y  = PROTECT(coerceVector(y, REALSXP));
    xp = PROTECT(coerceVector(xp, REALSXP));

    R_xlen_t nx = XLENGTH(x), np = XLENGTH(xp);
    SEXP yp = PROTECT(allocVector(REALSXP, np));

    BDRkquantile(REAL(x), REAL(y), nx,
                 REAL(xp), REAL(yp), np,
                 krn, bw, tau);

    SEXP ans = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, xp);
    SET_VECTOR_ELT(ans, 1, yp);

    SEXP nm = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(nm, 0, mkChar("x"));
    SET_STRING_ELT(nm, 1, mkChar("y"));
    setAttrib(ans, R_NamesSymbol, nm);

    UNPROTECT(6);
    return ans;
}
