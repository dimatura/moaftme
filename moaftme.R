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
        print(m)

        out.wz <- sample.wz(t.l, t.u, 
                beta.samples[m-1,,], 
                rho.samples[m-1,,],
                sigma.samples[m-1,],
                right.censored,
                int.censored,
                Z, X)
        w <- out.wz$w
        Z <- out.wz$Z

        #w <- sample.w(t.l, t.u, right.censored, int.censored,
        #        Z.samples[m-1,,],
        #        X,
        #        beta.samples[m-1,,],
        #        sigma.samples[m-1,])

        #Z.samples[m,,] <- sample.z.2(w.samples[m,], X,
        #        beta.samples[m-1,,], 
        #        rho.samples[m-1,,],
        #        sigma.samples[m-1,])
        
        # sample.beta uses log of time
        Y <- log(w)
        out <- sample.beta(Y, X, Z, n0, S0, b0, B0, J)
        beta.samples[m,,] <- out$beta
        sigma.samples[m,] <- out$sigma
        
        out <- sample.rho(X, Z, rho.samples[m-1,,], gamma0, tune$rho)
        rho.samples[m,,] <- out$rho

        accepts.rho <- accepts.rho + out$accept
    }

    accepts.rho <- (accepts.rho*100)/m

    list(sigma.samples=sigma.samples,
         beta.samples=beta.samples,
         rho.samples=rho.samples,
         accepts.rho=accepts.rho)
}

# return .w, 1xp matrix of pseudo-observations
sample.w <- function(t.l, t.u, right.censored, int.censored, Z, X,
        .beta, .sigma) {
    # convert from 0,1 to T,F
    Z <- (Z == 1)
    # by default (no censorship) .w[i] = t.l
    .w <- t.l
    for (i in 1:nrow(t.l)) {
        xtb <- crossprod(X[i,], .beta[Z[i,],])
        s <- .sigma[Z[i,]]
        if (int.censored[i]==1) {
            FlFu <- plnorm(c(t.l[i], t.u[i]), xtb, s)
            # handle case when Fl and Fu ~ 1
            if (FlFu[1] > (1-1e-6)) {
                .w[i] <- runif(1, t.l[i], t.u[i])
            } else {
                .w[i] <- qlnorm(runif(1, FlFu[1], FlFu[2]), xtb, s)
            }
        } else if (right.censored[i]==1) {
            Fl <- plnorm(t.l[i], xtb, s)
            .w[i] <- qlnorm(runif(1, Fl, 1 ), xtb, s)
            .w[i] <- ifelse(is.infinite(.w[i]), t.l[i], .w[i])
        }
        #print(.w[i])
        if (is.infinite(.w[i])) {
            print('Infinite')
        }
        if (is.nan(.w[i])) {
            print('nan')
        }
        if (is.infinite(.w[i]) || is.nan(.w[i])) {
            print(.beta[Z[i,],])
            print(s)
            print(xtb)
            print(sprintf("i=%d, Fl=%f, Fu=%f, rc=%d", i, plnorm(t.l[i], xtb, s), plnorm(t.u[i], xtb, s), right.censored[i]))
        }
 
    }
    .w
}

sample.wz <- function(t.l, t.u, .beta, .rho, .sigma, 
                      right.censored, int.censored, Z, X) {

    Xrho <- tcrossprod(X, .rho)
    Xbeta <- tcrossprod(X, .beta)

    out.w <- .C("dm_sample_w", 
              tl=as.double(t.l),
              tu=as.double(t.u),
              right_censored=as.integer(right.censored),
              int_censored=as.integer(int.censored),
              Z=as.integer(Z),
              Xbeta=as.double(Xbeta),
              Xrho=as.double(Xrho),
              sigma=as.double(.sigma),
              nptr=as.integer(nrow(X)),
              pptr=as.integer(ncol(X)),
              Jptr=as.integer(nrow(.beta)),
              w=double(nrow(X)))
    w <- out.w$w
    #browser()
    
    out.z <- .C("dm_sample_z", 
              w=as.double(w),
              Xbeta=as.double(Xbeta),
              Xrho=as.double(Xrho),
              sigma=as.double(.sigma),
              nptr=as.integer(nrow(X)),
              pptr=as.integer(ncol(X)),
              Jptr=as.integer(nrow(.beta)),
              Z=integer(nrow(X)))

    Z <- out.z$Z
    #browser()

    list(w=w,Z=Z)
}
 
dF.i <- function(t.i, x.i, .beta, .rho, .sigma) {
    .p <- exp(tcrossprod(x.i, .rho))
    .p <- .p/repmat(rowSums(.p), 1, ncol(.p))
    tcrossprod(.p,plnorm(t.i, tcrossprod(x.i,.beta) ,.sigma))
}

# return nxJ matrix Z where Z[i,j]=1 indicates obs. i comes from component j
sample.z <- function(w, X, .beta, .rho, .sigma) {
    n <- nrow(X)
    p <- ncol(X)
    J <- ncol(.beta)

    # .p: nxj matrix where .p(i,j) = p(x.i, rho.j)
    #xtr <- tcrossprod(X, .rho)
    .p <- exp(tcrossprod(X, .rho))
    .f <- dlnorm(repmat(w, 1, J), tcrossprod(X,.beta), repmat(matrix(.sigma,nrow=1), n, 1))

    # logs to (try to) avoid overflow
    #.p <- exp(xtr - repmat(log(rowSums(exp(xtr))), 1, nrow(.rho)))

    # .f: nxj matrix where .f(i,j) = .f_j(x.i, .beta.j) 
    #.f <- dlnorm(repmat(w, 1, ncol(.p)), tcrossprod(X,.beta), repmat(matrix(.sigma,nrow=1), nrow(X), 1))
    # h: nxj matrix where

    .h <- .f * .p

    if (any(is.nan(.h))) {
        #print("nan")
        #print(xtr)
        ##print(rowSums(exp(xtr)))
        #print(.p)
        #print(.f)
    }

    Z <- matrix(nrow=nrow(X), ncol=nrow(.beta))
    for (i in 1:nrow(X)) {
        Z[i,] <- rmultinom(1, 1, .h[i,])
    }

    if (any(colSums(Z)==0)) {
        print("Empty z")
        browser()
    }

    Z
}

# return nxJ matrix Z where Z[i,j]=1 indicates obs. i comes from component j
sample.z.2 <- function(w, X, .beta, .rho, .sigma) {
    Xrhot <- tcrossprod(X, .rho)
    Xbetat <- tcrossprod(X, .beta)
    print ('sameplz')
    out <- .C("sample_z", 
              w=as.double(w),
              Xbetat=as.double(Xbetat),
              Xrhot=as.double(Xrhot),
              sigma=as.double(.sigma),
              nptr=as.integer(nrow(X)),
              pptr=as.integer(ncol(X)),
              Jptr=as.integer(nrow(.beta)),
              Z=integer(nrow(X)))
    browser()
      #print (out)
    h <- matrix(out$h, nrow=nrow(X))
    h <- h / repmat(rowSums(h), 1, ncol(h))
    .p <- exp(tcrossprod(X, .rho))
    .f <- dlnorm(repmat(w, 1, nrow(.beta)), tcrossprod(X,.beta), repmat(matrix(.sigma,nrow=1), nrow(X), 1))
    .h <- .p * .f
    .h <- .h / repmat(rowSums(.h), 1, ncol(.h))
    #browser()

    Z <- matrix(nrow=nrow(X), ncol=nrow(.beta))
    for (i in 1:nrow(X)) {
        Z[i,] <- rmultinom(1, 1, h[i,])
    }
    Z
}

# Y: log(w)
# z: column of memberships for jth component
sample.beta <- function(Y, X, Z, n0, S0, b0, B0, J) {
    #J x p
    out.beta <- matrix(NA, nrow=J, ncol=ncol(X))
    #1 x J
    out.sigma <- matrix(NA, nrow=1, ncol=J)
    # note: R uses 1-based indexing, Z is 0-based
    for (j in 0:(J-1)) {
        #out <- sample.beta.j(Y, X, Z.samples[m,,j], sigma.samples[m-1,j]^2, n0, S0, b0, B0)  
        #sigma.samples[m,j] <- sqrt(out$sig2j)
        #beta.samples[m,j,] <- t(out$beta.j)
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

        #S2 <- t(Y.s) %*% (diag(n) - X.s %*% ginv(t(X.s)%*% X.s)%*%t(X.s) )%*% Y.s / (n - ncol(X.s))
        S2 <- crossprod(Y.s, (diag(n) - X.s %*% invXtX %*% t(X.s))) %*% Y.s / (n - ncol(X.s))

        n1S1 <- n0*S0 + (n-p)*S2[1,1] + t(beta.h - b0) %*% ginv(B0 + invXtX) %*% (beta.h - b0)
        out.sigma[j+1] <- sqrt(1/rgamma(1, n1/2, n1S1/2))
        #out.beta[j,] <- t(rmvt(1,B1,n1)) + b1
        out.beta[j+1,] <- rmvnorm(1, b1, out.sigma[j+1] * B1)
        if (is.nan(out.sigma[j+1])) {
            print("Sigma nan")
            browser()
        }
    }
    list(beta=out.beta,sigma=out.sigma)
}

sample.beta.2 <- function(Y, X, Z) {
    #J x p
    out.beta <- matrix(NA, nrow=ncol(Z), ncol=ncol(X))
    #1 x J
    out.sigma <- matrix(NA, nrow=1, ncol=ncol(Z))

    # 0,1 -> T,F for indexing
    Z <- (Z==1)
    for (j in 1:nrow(out.beta)) {
        X.s <- X[Z[,j],]
        Y.s <- matrix(Y[Z[,j]],ncol=1)
        fit = lm(Y.s ~ 0 + X.s)
        .beta.hat <- matrix(fit$coef, c(1, fit$rank))
        s2 <- sum(fit$residuals^2)/fit$df.residual
        shape <- fit$df.residual/2
        rate <- fit$df.residual/2 * s2
        vbeta <- vcov(fit)/s2
        out.sigma[j] <- sqrt(1/rgamma(1, shape=shape, rate=rate))
        print(.beta.hat)
        print(vbeta)
        .beta <- mvrnorm(1, rep(0, ncol(.beta.hat)), vbeta) 
        out.beta[j,] <- .beta.hat + .beta * out.sigma[j]
    }
    #log.post <- function(.beta) {
    #    P1 <- X %*% solve(t(X.s) %*% X) %*% t(X)
    #    -(p+1)/2 * log(n+1) - n/2 * log(t(y) %*% y -n/(n+1) * t(y) %*% P1 %*% y - 1/(n+1) * t(.beta) %*% t(X) %*% P1 %*% X %*% .beta)
    #}
    list(beta=out.beta, sigma=out.sigma)
}

#tmpitr <- 1

# .rho: rho matrix Jxp
# tune: variance of proposal distrib
#  Z: nxJ membership matrix
# returns 1
sample.rho <- function(X, Z, .rho, gamma0, tune) {
    n <- nrow(X)
    p <- ncol(X)
    J <- ncol(Z)
    # log posteriors of rows of a rho matrix
    # as a vector of length J
    log.post <- function(r) {
        # TODO flatness for log.pri as parameter
        log.pri <- dmvnorm(r, rep(0, p), gamma0*diag(p),log=TRUE)
        exr <- exp(tcrossprod(X, r))
        exr1 <- exr
        xr <- tcrossprod(X, r)
        for(j in 1:J) {
         exr1[,j] <- exp( Z[,j] * xr[,j] )
        }
        exr <- exr1/repmat(rowSums(exr), 1, ncol(exr))
        log.lik <- colSums(log(exr))
        log.lik + log.pri
    }

    cand.rho <- matrix(0, nrow=nrow(.rho), ncol=ncol(.rho))
    # note we leave first row in 0 for identifiablity
    for (j in 2:nrow(.rho)) {
        cand.rho[j,] <- rmvnorm(1, .rho[j,], tune*diag(ncol(.rho)))
    }

    # vector of ratios, 1xJ
    ratio <- exp(log.post(cand.rho) - log.post(.rho)) 

    accept <- 0
    for (j in 2:nrow(.rho)) {
        if (runif(1) < ratio[j]) {
            accept <- accept + 1
            .rho[j,] <- cand.rho[j,] 
        } 
    }

    exr <- exp(tcrossprod(X, .rho))
    exr <- exr/repmat(rowSums(exr), 1, ncol(exr))
    matplot(y=exr, type='l')
    #dev.print(device=pdf,sprintf("itr%d.pdf", tmpitr))
    #tmpitr <<- tmpitr + 1
    #matplot(t(.rho), type='l')
    list(rho=.rho,accept=accept)
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
    list(t.l=t.l,t.u=t.u,
            right.censored=right.censored,
            int.censored=int.censored,
            z=z,
            beta=.beta,
            rho=.rho,
            sigma=.sigma,
            X=X)
}


test.sample.w <- function() {
    n <- 40
    X <- matrix(0.1, nrow=n)
    x.i <- X[1,]
    .beta <- matrix(c(1), nrow=1)
    .sigma <- matrix(c(.2), nrow=1)

    t. <- seq(0, 2, 0.001) 
    Z <- matrix(1,nrow=nrow(X))

    # all int censored
    rc <- rep(0, length(t.))
    ic <- rep(1, length(t.))

    l <- 1.0
    u <- 1.4
    t.l <- matrix(rep(l, nrow(X)), ncol=1)
    t.u <- matrix(rep(u, nrow(X)), ncol=1)
    w <- sample.w(t.l, t.u, rc, ic, Z, X, .beta, .sigma) 

    xtb <- x.i*.beta
    sdlog <- .sigma
    F. <- plnorm(t., xtb, sdlog)
    f. <- dlnorm(t., xtb, sdlog)

    plot(t., F., type='l', col="red", xlim=c(0, max(t.)))
    lines(t., f./2, col="black", lty=4)
    lines(c(l, l), c(0, 10), col="gray")
    lines(c(u, u), c(0, 10), col="gray")
    Fl <- plnorm(l, xtb, sdlog)
    Fu <- plnorm(u, xtb, sdlog)
    Fw <- plnorm(w, xtb, sdlog)
    lines(c(0, max(t.)), c(Fl, Fl), col="green")
    lines(c(0, max(t.)), c(Fu, Fu), col="blue")
    points(w, rep(0, length(w)), pch='|')
    points(rep(0, length(Fw)), Fw, pch='-')
    #dev.print(device=pdf,"samplew.pdf")
}

test.sampler.1 <- function() {
    n0 <- 1
    S0 <- diag(100)
    b0 <- 1
    B0 <- 1
    tune <- list(sigma=0.005, beta=0.005, rho=0.35)
    M <- 1000
    J <- 3

    sim <- sim.data()
    X <- cbind(rep(1,nrow(sim$X)), sim$X)

    rc <- sim$right.censored
    ic <- sim$int.censored

    out <- moaftme.sampler(sim$t.l, sim$t.u, rc, ic, X, J, M, 
        n0=1, S0=2, b0=matrix(0,2,1), B0=100*diag(2), tune)
    #out <- moaftme.MCEM(sim$t.l, sim$t.u, rc, ic, X, J, M)
    out
}

test.sampler.2 <- function() {
    n0 <- 1
    S0 <- diag(100)
    b0 <- 1
    B0 <- 1
    tune <- list(sigma=1, beta=1, rho=0.25)
    M <- 400
    J <- 3

    d <- read.table('flourbeetle.txt', header=TRUE)

    X <- model.matrix(~ X1, d)

    rc <- d$rc
    ic <- d$ic
    t.l <- matrix(d$t.l,ncol=1)
    t.u <- matrix(d$t.u,ncol=1)

    out <- moaftme.sampler(t.l, t.u, rc, ic, X, J, M, 
           n0=1, S0=2, b0=matrix(0,2,1), B0=100*diag(2), tune)
    # out <- moaftme.MCEM(t.l, t.u, rc, ic, X, J, M)
    out
}

if (FALSE) {
    ts.plot(out$beta[,1,1])
    plot(density(out$beta[,1,1]))
    ts.plot(out$beta[,1,2])
    ts.plot(out$beta[,2,1])
    ts.plot(out$beta[,2,2])
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

#test.sample.rho()
#test.sample.w()
#test.sample.z()
#sim.data(TRUE)
#package.skeleton(name="moaftme", namespace=TRUE)
out <- test.sampler.1()
#out <- test.sampler.2()

