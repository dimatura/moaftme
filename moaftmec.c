
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
        //Rprintf("Z[%d] = %d\n", i, Z[i]);
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

double dm_log_post_rho_j(int j,
                     double *rho_j,
                     double *rho,
                     double *X, 
                     int *Z, 
                     double *gamma0,
                     double *tune,
                     int n,
                     int p,
                     int J
                     ) {
    int i, k, l;

    // evaluate logpriori(cand_rho) 
    double log_pri = 0;
    for (k=0; k < p; ++k) {
        // dnorm args are x, mu, sigma, log
        log_pri += dnorm(rho_j[k], 0, *gamma0, 1);
    }
    // evaluate loglik rho j
    // loglik = sum_i [(xi * rho_j) - log(sum_l(exp(xi * rho_k)))] 
    double log_lik = 0;
    for (i=0; i < n; ++i) {
        if (Z[i]!=j) {
            continue; // only care about j candidates
        }
        // calculating log(exp(x_i rho)) = x_i rho
        double xr = 0;
        for (k=0; k < p; ++k) {
            xr += X[i + n*k] * rho_j[k];
        }
        log_lik += xr;
        // calculating sum_l(exp(x_i rho)))
        double sum_exp_xrl = 0;
        for (l = 0; l < J; ++l) {
            if (l==j) {
                // this is the new one
                sum_exp_xrl += exp(xr);
                continue; 
            } 
            double xrl = 0;
            for (k=0; k < p; ++k) {
                xrl += X[i + n*k] * rho[l + J*k];
            }
            sum_exp_xrl += exp(xrl);
        }
        log_lik -= log(sum_exp_xrl);
    }
    // log post = log lik + log pri
    return log_pri + log_lik;
}

void dm_sample_rho(double *X, 
                  int *Z,
                  double *rho, // note: input and output
                  double *gamma0,
                  double *tune,
                  int *nptr, 
                  int *pptr,
                  int *Jptr,
                  int *accept 
                  ) {
    int n = *nptr;
    int p = *pptr;
    int J = *Jptr;
    int i, j, k;
    // for each j
    //    generate candidate rho_j
    //    evaluate logpost(rho_j candidate) 
    //    evaluate logpost(rho_j current) 
    //    evaluate alpha
    //    accept or reject rho_j candidate

    double *cand_rho = (double *) R_alloc(p, sizeof(double));
    double *current_rho = (double *) R_alloc(p, sizeof(double));
    // we start from 1 because first rho is always zero.
    for (j=1; j < J; ++j) {
        // generate candidate rho_j
        for (k=0; k < p; ++k) {
            cand_rho[k] = rnorm(rho[j + J*k], *tune);
        }

        double log_post_cand_rho = dm_log_post_rho_j(j,
                                                     cand_rho,
                                                     rho,
                                                     X,
                                                     Z,
                                                     gamma0,
                                                     tune,
                                                     n,
                                                     p,
                                                     J);
        // extract current rho_j; this is inefficient
        for (k=0; k < p; ++k) {
            current_rho[k] = rho[j + J*k];
        }

        // this should be cached
        double log_post_current_rho = dm_log_post_rho_j(j,
                                                        current_rho,
                                                        rho,
                                                        X,
                                                        Z,
                                                        gamma0,
                                                        tune,
                                                        n,
                                                        p,
                                                        J);
        double ratio = exp(log_post_cand_rho - log_post_current_rho);
        if (runif(0, 1) < ratio) {
            // copy cand_rho into rho
            for (k=0; k < p; ++k) {
                rho[j + J*k] = cand_rho[k];
            }
            *accept = *accept + 1;
        }
        // else rho stays the same.
    }

}

