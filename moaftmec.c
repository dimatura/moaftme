
#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>

#define DM_SQR(x) ((x)*(x))

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

void sample_z(double *w,
        double *Xbetat,
        double *Xrhot,
        double *sigma,
        int *nptr, 
        int *pptr,
        int *Jptr,
        double *h
        ) {
    int n = *nptr;
    int p = *pptr;
    int J = *Jptr;
    int i, j, k;
    double *sigma2 = (double *) R_alloc(J, sizeof(double));
    double *logw = (double *) R_alloc(n, sizeof(double));
    for (j=0; j < J; ++j) {
        sigma2[j] = sigma[j]*sigma[j];
    }
    for (i=0; i < n; ++i) {
        logw[i] = log(w[i]);
    }
    for (j=0; j < J; ++j) {
        for (i=0; i< n; ++i) {
            h[i + n*j] = (1/(sigma[j]*w[i]))*exp(Xrhot[i+n*j]-(0.5/sigma2[j])*DM_SQR(logw[i]-Xbetat[i+n*j]));
            //Rprintf("%f\n", h[i + n*j]);
        }
    }
}
