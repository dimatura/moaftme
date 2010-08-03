library(mvtnorm)
library(MASS)
options(error=dump.frames)

dyn.load('moaftmec.so')

# t - registered time, no censorship, nx1
# t.l, t.u - lower and upper bounds in interval-censored data, nx1 each 
# if data is uncensored, they should be equal (but only t.l is used)
# right.censored - indicator, nx1
# interv.censored - indicator, nx1
# X covariates - nxp
# w - simulated latent times, nx1
# Z - indicator of latent component, n x J
# beta - coeffcients of AFT models, J x p
# sigma - sigma of AFT models, J x 1 
# rho - coefficients of logistic link, J x p

moaftme.sampler <- function(t.l, t.u, right.censored, int.censored, X, J, M, 
        n0, S0, b0, B0, gamma0, tune) {
    n <- nrow(X)
    p <- ncol(X)

    invB0 <- solve(B0)

    w <- rep(NA, n)
    # 0:(J-1) for use with C indexing
    Z <- sample(0:(J-1), size=n, replace=TRUE)

    sigma.samples <- matrix(nrow=M, ncol=J)
    beta.samples <- array(NA, c(M, J, p))
    rho.samples <- array(NA, c(M, J, p))

    beta.samples[1,,] <- matrix(rnorm(p*J, mean=0, sd=1), ncol=p)
    rho.samples[1,,] <- matrix(rnorm(p*J, mean=0, sd=1), ncol=p)
    # set first col to zero for identifiability
    rho.samples[1,1,] <- matrix(0, nrow=p, ncol=1)
    #sigma.samples[1,] <- matrix(1/rgamma(J, n0/2, n0*S0/2), nrow=1)
    sigma.samples[1,] <- matrix(1, nrow=1, ncol=J)

    accepts <- rep(0, J)

    accepts.rho <- 0
    
    # main loop
    for (m in 2:M) {
        if (m %% 20 == 0) {
            print(m)
        }

        out.wz <- sample.wz(t.l, t.u, 
                beta.samples[m-1,,], 
                rho.samples[m-1,,],
                sigma.samples[m-1,],
                right.censored,
                int.censored,
                Z, X)
        w <- out.wz$w
        Z <- out.wz$Z

        # sample.beta uses log of time
        Y <- log(w)
        out.beta <- sample.beta.sigma(Y, X, Z, beta.samples[m-1,,], sigma.samples[m-1,], n0, S0, b0, B0, invB0, J)
        beta.samples[m,,] <- out.beta$beta
        sigma.samples[m,] <- out.beta$sigma
        
        out.rho <- sample.rho(X, Z, rho.samples[m-1,,], gamma0, tune, J)
        rho.samples[m,,] <- out.rho$rho

        accepts.rho <- accepts.rho + out.rho$accept
        #print(sprintf("m: %d, accepts: %d", m, accepts.rho))
    }

    # take into account that for each iteration we sample rho J-1 times
    accepts.rho <- (accepts.rho*100)/((J-1)*(M-1))

    list(sigma.samples=sigma.samples,
         beta.samples=beta.samples,
         rho.samples=rho.samples,
         accepts.rho=accepts.rho)

}

sample.wz <- function(t.l, t.u, .beta, .rho, .sigma, 
                      right.censored, int.censored, Z, X) {

    Xrho <- tcrossprod(X, .rho)
    Xbeta <- tcrossprod(X, .beta)

    out <- .C("moaftme_sample_wz",
              tl=as.double(t.l),
              tu=as.double(t.u),
              right_censored=as.integer(right.censored),
              int_censored=as.integer(int.censored),
              w=double(nrow(X)),
              Xbeta=as.double(Xbeta),
              Xrho=as.double(Xrho),
              sigma=as.double(.sigma),
              nptr=as.integer(nrow(X)),
              pptr=as.integer(ncol(X)),
              Jptr=as.integer(nrow(.beta)),
              Z=as.integer(Z))
    w <- out$w
    Z <- out$Z
    list(w=w,Z=Z)
}

sample.beta.sigma <- function(Y, X, Z, .beta, .sigma, n0, S0, b0, B0, invB0, J) {
    # TODO fix the t()
    .beta <- t(.beta)
    b0 <- t(b0)
    #J x p
    out.beta <- matrix(NA, nrow=J, ncol=ncol(X))
    #1 x J
    .p <- ncol(X)
    out.sigma <- matrix(NA, nrow=1, ncol=J)
    for (j in 0:(J-1)) {
        X.j <- matrix(X[Z==j,],ncol=.p)
        Y.j <- Y[Z==j]
        XtX <- crossprod(X.j)
        XtY <- crossprod(X.j, Y.j)
        .beta.j <- .beta[,j+1]

        # sample sigma2
        err <- sum((X.j %*% (-.beta.j) + Y.j)^2)
        c.post <- (n0 + nrow(X.j))*0.5
        d.post <- (S0 + err)*0.5
        # inverse gamma
        sigma2.t <- 1./rgamma(1, shape=c.post, rate=d.post)

        # sample beta, given sigma2
        sigma2.inv <- 1/sigma2.t
        sig.beta <- solve(B0 + XtX * sigma2.inv)
        .C <- chol(sig.beta)
        beta.hat <- sig.beta %*% (B0 %*% b0 + (XtY * sigma2.inv))
        beta.t <- (.C %*% rnorm(.p, 0, 1)) + beta.hat

        if (all(is.finite(beta.t)) && is.finite(sigma2.t)) {
            out.sigma[j+1] <- sqrt(sigma2.t)
            # TODO fix t()
            out.beta[j+1,] <- t(beta.t)
        } else {
            # too few points assigned to expert -- sample from prior
            #Zsum <- c()
            #for (k in 0:(J-1)) {
            #    Zsum <- c(Zsum, sum(Z==k))
            #}
            #print(Zsum)
            #out.sigma[j+1] <- sqrt(1./rgamma(1., shape=n0*.5, rate=S0*.5))
            out.sigma[j+1] <- .sigma[j+1]
            # TODO fix t()
            #out.beta[j+1,] <- t(.beta[j+1,])
            out.beta[j+1,] <- rmvnorm(1, b0, invB0)
        }
    }
    list(beta=out.beta,sigma=out.sigma)
}

sample.rho <- function(X, Z, .rho, gamma0, tune, J) {
    #print(.rho)
    n <- nrow(X)
    p <- ncol(X)
    out  <- .C("moaftme_sample_rho", 
              X=as.double(X),
              Z=as.integer(Z),
              rho=as.double(.rho),
              gamma0=as.double(gamma0),
              tune=as.double(tune),
              nptr=as.integer(nrow(X)),
              pptr=as.integer(ncol(X)),
              Jptr=as.integer(J),
              accept=integer(1))
    rho <- matrix(out$rho, nrow=J)

    ##PLOT
    exr <- exp(tcrossprod(X, .rho))
    exr <- exr/repmat(rowSums(exr), 1, ncol(exr))
    matplot(y=exr, type='l', ylim=c(0,1))

    list(rho=rho, accept=out$accept)
}

indicator.to.index <- function(ind) {
    wh <- which(ind==1, arr.ind=TRUE)
    wh[sort(wh[,1],index.return=TRUE)$ix,2]
}

index.to.indicator <- function(Z,J) {
    out <- matrix(0, nrow=length(Z), ncol=J)
    out[cbind(1:nrow(out),Z)] <- 1
    out
}

# matlab-like
repmat <- function(a,n,m) { kronecker(matrix(1, n, m), a) }

# F(t_i | x_i, theta)
sim.data <- function(plot=FALSE) {
    n <- 1000
    X <- matrix(sort(runif(n, -2, 2)),ncol=1)
    .beta <- matrix(c(1,-1), nrow=2)
    .rho <- matrix(c(16, 8), nrow=2)
    .sigma <- matrix(c(.2, .1), nrow=1)
    .p <- 1/(1+exp(-.rho[1]*X))
    z <- rbinom(length(X), 1,  .p) 
    t0 <- exp(.beta[1]*X[z==0])*exp(rnorm(sum(z==0))*.sigma[1]) 
    t1 <- exp(.beta[2]*X[z==1])*exp(rnorm(sum(z==1))*.sigma[2]) 
    if (plot) {
        plot(X, .p, type="l", ylab="t", col="green", xlim=c(-2,2), ylim=c(0, 2))
        points(X[z==0], t0, col="blue", xlim=c(-2,2))
        points(X[z==1], t1, col="red",  xlim=c(-2,2))
        #dev.print(device=pdf,"simdata.pdf")
    }
    t.l <- matrix(c(t0, t1),ncol=1)
    t.u <- t.l
    right.censored <- rep(0, n)
    int.censored <- rep(0,n)
    # censor
    for (i in 1:n) {
        r <- runif(1)
        if (r < 1/3) {
            # interval
            t.l[i] <- ifelse(t.l[i]-0.1 < 0, 0, t.l[i] - 0.1)
            t.u[i] <- t.u[i] + 0.1
            int.censored[i] <- 1
        } else if (r < 2/3) {
            # right
            t.l[i] <- ifelse(t.l[i]-0.1 < 0, 0, t.l[i] - 0.1)
            t.u[i] <- -1
            right.censored[i] <- 1
        }
        # else no censorship
    }
    list(t.l=t.l,t.u=t.u,
            right.censored=right.censored,
            int.censored=int.censored,
            z=z,
            beta=.beta,
            rho=.rho,
            sigma=.sigma,
            X=X)
}


