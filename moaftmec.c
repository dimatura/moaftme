
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
    Rprintf("sum weight: %f, J: %d\n", sum_weights, J);
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

void sample_z(double *w,
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
    int i, j, k;
    double *sigma2 = (double *) R_alloc(J, sizeof(double));
    double *logw = (double *) R_alloc(n, sizeof(double));
    Rprintf("n=%d p=%d J=%d\n", n, p, J);
    for (j=0; j < J; ++j) {
        sigma2[j] = DM_SQR(sigma[j]);
    }
    for (i=0; i < n; ++i) {
        logw[i] = log(w[i]);
    }
    Rprintf("before loop\n");
    // unnormalized row of hij
    double *hi = (double *) R_alloc(J, sizeof(double));
    for (i=0; i < n; ++i) {
        Rprintf("i=%d\n", i);
        // row sum for hij
        double hi_sum=0;
        // reset hi
        memset(hi, 0, sizeof(double)*J);
        Rprintf("after memset\n");
        for (j=0; j< J; ++j) {
            Rprintf("j=%d\n", j);
            hi[j] = (1/(sigma[j]*w[i]))*exp(Xrho[i+n*j]-(0.5/sigma2[j])*DM_SQR(logw[i]-Xbeta[i+n*j]));
            hi_sum += hi[j];
            //Rprintf("%f\n", h[i + n*j]);
        }
        Z[i] = dm_rmultinom(hi, hi_sum, J);
    }
}
