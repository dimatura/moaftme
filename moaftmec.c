
#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>

#include <string.h>

#define DM_SQR(x) ((x)*(x))
#define DM_MULTINOM_EPS (1e-10) 

void foo(double *a, double *b, int *lenptr, double *out) {
    int i;
    int len = *lenptr;
    for (int i=0; i < len; ++i) {
        out[i] = a[i]/b[i];
    }
}

// in fortran ordering, N dimensions d_0...d_(N-1), flat.ix = sum_{0,N-1} n_i prod_{0, i-1} d_j

void foo2(double *a, int *dim0, int *dim1, int *dim2) {
    int i,j,k;
    Rprintf("%d %d %d\n", *dim0, *dim1, *dim2);
    /*
    for(i=0; i < (*dim0)*(*dim1)*(*dim2); ++i) {
        double v = a[i];
        Rprintf("%f\n", v);
    }
    */
    for (k=0; k < *dim2; ++k) {
        for (j=0; j < *dim1; ++j) {
            for (i=0; i < *dim0; ++i) {
                double v = a[i + (*dim0)*j + (*dim0)*(*dim1)*k];
                Rprintf("%f\n", v);
            }
        }
    }
}

// weighted random choice from [0,J)
int dm_rmultinom(double *weights, double sum_weights, int J) {
    double r;
    int j;
    //Rprintf("sum weight: %f, J: %d\n", sum_weights, J);
    if (sum_weights < DM_MULTINOM_EPS) {
        // simply assume equiprobability
        for (j=0; j < J; ++j) {
            weights[j] = 1;
        }
        sum_weights = J;
    }
    r = runif(0, sum_weights);
    for (j=0; j < J; ++j) {
        r -= weights[j];
        if (r < 0) {
            return j;
        }
    }
    return J-1;
}

void dm_sample_zw(double *w,
        double *Xbeta,
        double *Xrho,
        double *sigma,
        int *nptr, 
        int *pptr,
        int *Jptr,
        int *Z
        ) {
    int n = *nptr;
    int p = *pptr;
    int J = *Jptr;
}

void dm_sample_z(double *w,
                 double *Xbeta,
                 double *Xrho,
                 double *sigma,
                 int *nptr, 
                 int *pptr,
                 int *Jptr,
                 int *Z) {
    int n = *nptr;
    int p = *pptr;
    int J = *Jptr;
    int i, j, k;
    // This memory is freed by R at end of .C call
    double *sigma2 = (double *) R_alloc(J, sizeof(double));
    double *logw = (double *) R_alloc(n, sizeof(double));
    Rprintf("n=%d p=%d J=%d\n", n, p, J);
    for (j=0; j < J; ++j) {
        sigma2[j] = DM_SQR(sigma[j]);
    }
    for (i=0; i < n; ++i) {
        logw[i] = log(w[i]);
    }
    // unnormalized row of hij
    double *hi = (double *) R_alloc(J, sizeof(double));
    for (i=0; i < n; ++i) {
        // row sum for hij
        double hi_sum=0;
        // reset hi
        memset(hi, 0, sizeof(double)*J);
        for (j=0; j< J; ++j) {
            hi[j] = (1/(sigma[j]*w[i]))*exp(Xrho[i+n*j]-(0.5/sigma2[j])*DM_SQR(logw[i]-Xbeta[i+n*j]));
            hi_sum += hi[j];
        }
        Z[i] = dm_rmultinom(hi, hi_sum, J);
    }
}

void dm_sample_w(double *tl,
                 double *tu,
                 int *right_censored,
                 int *int_censored,
                 int *Z,
                 double *Xbeta,
                 double *Xrho,
                 double *sigma,
                 int *nptr, 
                 int *pptr,
                 int *Jptr,
                 double *w
                 ) {
    int n = *nptr;
    int p = *pptr;
    int J = *Jptr;
    int i, j, k;
    for (i=0; i < n; ++i) {
        Rprintf("Z[%d] = %d\n", i, Z[i]);
        if (right_censored[i]==0 && int_censored[i]==0) {
            w[i] = tl[i];
            continue;
        }
        double sj = sigma[Z[i]];
        // Xbeta and Xrho come in Fortran order 
        double xibj = Xbeta[i + Z[i]*n];
        double xirj = Xrho[i + Z[i]*n];
        if (int_censored[i]==1) {
            //plnorm args: x, mean, sigma, lowertail=TRUE, log=FALSE
            double Fl = plnorm(tl[i], xibj, sj, 1, 0);
            if (Fl > (1 - 1e-8)) {
                // tl is very large and both f(tl) and f(tu) are very small.
                // qlnorm would return Inf. We sample uniformly from [tl, tu]. 
                w[i] = runif(tl[i], tu[i]);
            }
            double Fu = plnorm(tu[i], xibj, sj, 1, 0);
            // not sure if log helps
            double logFw = log(runif(Fl, Fu));
            w[i] = qlnorm(logFw, xibj, sj, 1, 1);
        } else if (right_censored[i]==1) {
            double Fl = plnorm(tl[i], xibj, sj, 1, 0);
            double logFw = log(runif(Fl, 1));
            w[i] = qlnorm(logFw, xibj, sj, 1, 1);
            // deal with case when lognorm is degenerate (eg very large sigma)
            if (!R_FINITE(w[i])) {
                w[i] = tl[i];
            }
        }
    }
}


