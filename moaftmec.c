
#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>

#include <string.h>

#define MOAFTME_MULTINOM_EPS (1e-10) 

static R_INLINE double moaftme_sqr(double x) {
    return x*x;
}

static R_INLINE double moaftme_min(double x, double y) {
    return (x < y)? x : y;
}

static R_INLINE double moaftme_max(double x, double y) {
    return (x > y)? x : y;
}

// weighted random choice from [0,J)
int moaftme_rmultinom(double *weights,
        double sum_weights,
        int J) {
    double r;
    int j;
    if (sum_weights < MOAFTME_MULTINOM_EPS) {
        // simply assume equiprobability
        for (j=0; j < J; ++j) {
            weights[j] = 1;
        }
        sum_weights = J;
    }
    r = runif(0., sum_weights);
    for (j=0; j < J; ++j) {
        r -= weights[j];
        if (r < 0) {
            return j;
        }
    }
    return J-1;
}

void moaftme_sample_z(double *w,
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
    // unnormalized row of hij
    double *hi = (double *) R_alloc(J, sizeof(double));
    for (i=0; i < n; ++i) {
        // row sum for hij
        double hi_sum=0;
        // reset hi
        for (j=0; j< J; ++j) {
            hi[j] = (1./(sigma[j]*w[i]))*exp(Xrho[i+n*j]
                    -(0.5/moaftme_sqr(sigma[j]))
                    *moaftme_sqr(log(w[i])-Xbeta[i+n*j]));
            hi_sum += hi[j];
        }
        Z[i] = moaftme_rmultinom(hi, hi_sum, J);
    }
}

double moaftme_sample_int_censored(double tl,
        double tu,
        double mean,
        double sigma) {
    //plnorm args: x, mean, sigma, lowertail=TRUE, log=FALSE
    double Fl = plnorm(tl, mean, sigma, 1, 0);
    if (Fl > (1 - 1e-8)) {
        // tl is very large and both f(tl) and f(tu) are very small.
        // qlnorm would return Inf. We sample uniformly from [tl, tu]. 
        return runif(tl, tu);
    }
    double Fu = plnorm(tu, mean, sigma, 1, 0);
    double Fw = runif(Fl, Fu);
    return qlnorm(Fw, mean, sigma, 1, 0);
    //Rprintf("i %f\n", w[i]);
}

double moaftme_sample_right_censored(double tl,
        double mean,
        double sigma) {
    double Fl = plnorm(tl, mean, sigma, 1, 0);
    double Fw = runif(Fl, 1);
    double w = qlnorm(Fw, mean, sigma, 1, 0);
    // deal with case when lognorm is degenerate (eg very large sigma)
    if (!R_FINITE(w)) {
        return tl;
    }
    return w;
}

void moaftme_sample_w_from_R(double *tl,
        double *tu,
        int *right_censored,
        double *mean,
        double *sigma,
        double *out) {
    GetRNGstate();
    if (*right_censored==1) {
        *out = moaftme_sample_right_censored(*tl, *mean, *sigma);
    } else {
        *out = moaftme_sample_int_censored(*tl, *tu, *mean, *sigma);
    }
    PutRNGstate();
}

void moaftme_sample_w(double *tl,
                 double *tu,
                 int *right_censored,
                 int *int_censored,
                 int *Z,
                 double *Xbeta,
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
        double xibj = Xbeta[i + Z[i]*n];
        if (int_censored[i]==1) {
            w[i] = moaftme_sample_int_censored(tl[i], 
                    tu[i],
                    xibj,
                    sj);
        } else if (right_censored[i]==1) {
            w[i] = moaftme_sample_right_censored(tl[i],
                    xibj,
                    sj);
        }
    }
}

void moaftme_sample_wz(double *tl,
        double *tu,
        int *right_censored,
        int *int_censored,
        double *w,
        double *Xbeta,
        double *Xrho,
        double *sigma,
        int *nptr, 
        int *pptr,
        int *Jptr,
        int *Z
        ) {
    GetRNGstate();
    moaftme_sample_w(tl, tu, right_censored, int_censored,
            Z, Xbeta, sigma, nptr, pptr, Jptr, w);
    moaftme_sample_z(w, Xbeta, Xrho, sigma,
            nptr, pptr, Jptr, Z);
    PutRNGstate();
}

double moaftme_log_post_rho_j(int j,
        double *rho_j,
        double *rho,
        double *X, 
        int *Z, 
        double *gamma0,
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
    // loglik = sum_i [z_ij(xi * rho_j) - log(sum_l(exp(xi * rho_k)))] 
    double log_lik = 0;
    for (i=0; i < n; ++i) {
        // calculating z_ij * x_i rho
        double xr = 0;
        for (k=0; k < p; ++k) {
            xr += X[i + n*k] * rho_j[k];
        }
        if (Z[i]==j) {
            log_lik += xr;
        }
        // calculating sum_l(exp(x_i rho)))
        // TODO
        // See numerical recipes on EM; use 'log-sum-exp' formula
        // log(sum(exp(z_i))) = z_max + log(sum(exp(z_i-z_max)))
        // guarantee at least on exp won't overflow, and those
        // that do could be neglected anyway.
        // Note that when z_i == z_max, exp(z_i-z_max) = 1
        // therefore the calculation will be
        // z_max + log(1 + sum_{z_i != z_max}(exp(z_i-z_max)))
        double sum_exp_xrl = exp(xr);
        for (l = 0; l < J; ++l) {
            if (l==j) { continue; }
            double xrl = 0;
            for (k=0; k < p; ++k) {
                xrl += (X[i + n*k] * rho[l + J*k]);
            }
            sum_exp_xrl += exp(xrl);
        }
        log_lik -= log(sum_exp_xrl);
    }
    return (log_pri + log_lik);
}

void moaftme_sample_rho(double *X, 
                  int *Z,
                  double *rho, // note: input and output
                  double *gamma0,
                  double *tune,
                  int *nptr, 
                  int *pptr,
                  int *Jptr,
                  int *accepts
                  ) {
    int n = *nptr;
    int p = *pptr;
    int J = *Jptr;
    int i, j, k;

    //Rprintf("n=%d p=%d J=%d\n", n, p, J);
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
            cand_rho[k] = rnorm(rho[j + J*k], tune[j]);
        }

        double log_post_cand_rho = moaftme_log_post_rho_j(j,
                                                     cand_rho,
                                                     rho,
                                                     X,
                                                     Z,
                                                     gamma0,
                                                     n,
                                                     p,
                                                     J);
        // a copy of current rho_j
        for (k=0; k < p; ++k) {
            current_rho[k] = rho[j + J*k];
        }

        double log_post_current_rho = moaftme_log_post_rho_j(j,
                                                        current_rho,
                                                        rho,
                                                        X,
                                                        Z,
                                                        gamma0,
                                                        n,
                                                        p,
                                                        J);

        double ratio = exp(log_post_cand_rho - log_post_current_rho);

        if (runif(0, 1) < ratio) {
            // copy cand_rho into rho
            for (k=0; k < p; ++k) {
                rho[j + J*k] = cand_rho[k];
            }
            accepts[j] += 1;
        }
        // else rho stays the same.
    }
}

// in fortran ordering, N dimensions d_0...d_(N-1), flat.ix = sum_{0,N-1} n_i prod_{0, i-1} d_j
void foo(double *a, int *dim0, int *dim1, int *dim2) {
    int i,j,k;
    Rprintf("%d %d %d\n", *dim0, *dim1, *dim2);
    for(i=0; i < (*dim0)*(*dim1)*(*dim2); ++i) {
        double v = a[i];
        Rprintf("%f ", v);
    }
    Rprintf("\n");
    for (k=0; k < *dim2; ++k) {
        for (j=0; j < *dim1; ++j) {
            for (i=0; i < *dim0; ++i) {
                double v = a[i + (*dim0)*j + (*dim0)*(*dim1)*k];
                Rprintf("%f ", v);
            }
        }
    }
}


