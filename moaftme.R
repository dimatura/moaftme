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
        n0, S0, b0, B0, tune) {
    n <- nrow(X)
    p <- ncol(X)

    w <- rep(NA, n)
    # we 0:(J-1) for use with C indexing
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

    alpha0 <- 1
    lambda0 <- 1
    gamma0 <- 10

    accepts.rho <- 0
    
    # main loop
    for (m in 2:M) {
        if (m %% 1 == 0) {
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
        #out.beta <- sample.beta(Y, X, Z, beta.samples[m-1,,], sigma.samples[m-1,], n0, S0, b0, B0, J)
        out.beta <- sample.beta.sigma(Y, X, Z, beta.samples[m-1,,], sigma.samples[m-1,], n0, S0, b0, B0, J)
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

    out <- .C("dm_sample_wz",
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
 
# Y: log(w)
# z: column of memberships for jth component
sample.beta <- function(Y, X, Z, .beta, .sigma, n0, S0, b0, B0, J) {
    #J x p
    out.beta <- matrix(NA, nrow=J, ncol=ncol(X))
    #1 x J
    out.sigma <- matrix(NA, nrow=1, ncol=J)
    # note: R uses 1-based indexing, Z is 0-based
    for (j in 0:(J-1)) {
        X.s <- matrix(X[Z==j,],nrow=sum(Z==j), ncol=ncol(X))
        Y.s <- Y[Z==j]
        n <-nrow(X.s)
        p <- ncol(X.s)
        XtX <- crossprod(X.s)
        invXtX <- ginv(XtX)
        XtY <- crossprod(X.s, Y.s)
        B1 <- ginv(ginv(B0) + XtX)
        b1 <- B1 %*% (ginv(B0) %*% b0 + XtY)
        beta.h <- ginv(XtX) %*% XtY 
        n1 <- n0 + n

        S2 <- crossprod(Y.s, (diag(n) - X.s %*% invXtX %*% t(X.s))) %*% Y.s / (n - ncol(X.s))

        n1S1 <- n0*S0 + (n-p)*S2[1,1] + t(beta.h - b0) %*% ginv(B0 + invXtX) %*% (beta.h - b0)
        #print(sprintf("n1S1: %f, S2: %f", n1S1, S2))

        if (is.finite(n1S1) && is.finite(S2)) {
            out.sigma[j+1] <- sqrt(1/rgamma(1, n1/2, n1S1/2))
            out.beta[j+1,] <- rmvnorm(1, b1, out.sigma[j+1] * B1)
        } else {
            out.sigma[j+1] <- .sigma[j+1]
            out.beta[j+1,] <- .beta[j+1,]
        }
    }
    list(beta=out.beta,sigma=out.sigma)
}

sample.beta.sigma <- function(Y, X, Z, .beta, .sigma, n0, S0, b0, B0, J) {
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

        #print(.beta, .sigma)
        #print(beta.t, sigma2.t)

        #if (all(is.finite(beta.t)) && is.finite(sigma2.t)) {
            out.sigma[j+1] <- sqrt(sigma2.t)
        #    # TODO fix t()
            out.beta[j+1,] <- t(beta.t)
        #} else {
        #    out.sigma[j+1] <- .sigma[j+1]
        #    # TODO fix t()
        #    out.beta[j+1,] <- t(.beta[j+1,])
        #}
    }
    list(beta=out.beta,sigma=out.sigma)
}

sample.rho <- function(X, Z, .rho, gamma0, tune, J) {
    #print(.rho)
    n <- nrow(X)
    p <- ncol(X)
    out  <- .C("dm_sample_rho", 
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

test.sampler.1 <- function() {
    n0 <- 0.01
    #S0 <- .5 
    S0 <- 0.01
    b0 <- matrix(0, nrow=1, ncol=2)
    #B0 <- diag(1, 2)
    B0 <- diag(0.01, 2)
    tune <- 0.10
    M <- 1000
    J <- 2

    sim <- sim.data()
    X <- cbind(rep(1,nrow(sim$X)), sim$X)

    rc <- sim$right.censored
    ic <- sim$int.censored

    out <- moaftme.sampler(sim$t.l, sim$t.u, rc, ic, X, J, M, 
        n0=n0, S0=S0, b0=b0, B0=B0, tune)
    out
}

test.sampler.2 <- function() {
    n0 <- 0.001
    #S0 <- diag(100)
    S0 <- 0.001
    b0 <- 1
    #B0 <- 1
    B0 <- diag(0.01, 2)
    tune <- 0.20
    M <- 100000
    J <- 2

    d <- read.table('flourbeetle.txt', header=TRUE)

    X <- model.matrix(~ X1, d)

    rc <- d$rc
    ic <- d$ic
    t.l <- matrix(d$t.l,ncol=1)
    t.u <- matrix(d$t.u,ncol=1)

    out <- moaftme.sampler(t.l, t.u, rc, ic, X, J, M, 
           n0=1, S0=2, b0=matrix(0,2,1), B0=100*diag(2), tune)
    out
}

if (FALSE) {
    ts.plot(out$beta[,1,1])
    plot(density(out$beta[,1,1]))
    ts.plot(out$beta[,1,2])
    plot(density(out$beta[,1,2]))
    ts.plot(out$beta[,2,1])
    plot(density(out$beta[,2,1]))
    ts.plot(out$beta[,2,2])
    plot(density(out$beta[,2,2]))
    ts.plot(out$sigma[,1])
    plot(density(out$sigma[,1]))
    ts.plot(out$sigma[,1])
    ts.plot(out$sigma[,2])
    ts.plot(out$rho[,2,1])
    ts.plot(out$rho[,3,1])
    matplot(cbind(out$beta[,1,1], 10*out$sigma[,1]), type='l')
    matplot(out$beta[,1,1:2], type='l')
    matplot(out$rho[,1,1:2], type='l')
    matplot(out$rho[,2,1:2], type='l')
    matplot(out$rho[,3,1:2], type='l')
    acf(out$beta[,1,1:2])
    acf(out$beta[,2,1:2])
    acf(out$rho[,2,1])
    acf(out$sigma[,1:2])
    
    ts.plot(out$sigma[100:1000,2])
    ts.plot(out$beta[100:1000,1,1:2])
    ts.plot(out$beta[100:1000,2,1:2])
    ts.plot(out$rho[100:1000,2,1])
    mean(out$rho[100:1000,2,2])

}

#sim.data(TRUE)
#package.skeleton(name="moaftme", namespace=TRUE)
out <- test.sampler.1()
#out <- test.sampler.2()

