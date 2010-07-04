
library(mvtnorm)
library(MASS)
options(error=dump.frames)

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

    # samples
    w.samples <- matrix(NA, nrow=M, ncol=n)
    Z.samples <- array(NA, c(M, n, J))
    sigma.samples <- matrix(nrow=M, ncol=J)
    # usar B0
    beta.samples <- array(NA, c(M, J, p))
    rho.samples <- array(NA, c(M, J, p))

    # initial sampling
    Z.samples[1,,] <- t(rmultinom(n, 1, rep(1/J, J)))
    beta.samples[1,,] <- matrix(rnorm(p*J, mean=0, sd=1), ncol=p)
    rho.samples[1,,] <- matrix(rnorm(p*J, mean=0, sd=1), ncol=p)
    # set first col to zero for identifiability
    rho.samples[1,1,] <- matrix(0, nrow=p, ncol=1)
    #sigma.samples[1,] <- matrix(1/rgamma(J, n0/2, n0*S0/2), nrow=1)
    sigma.samples[1,] <- matrix(1, nrow=1, ncol=J)
    w.samples[1,] <- matrix(0, nrow=n)

    accepts <- rep(0, J)

    # main loop
    for (m in 2:M) {
        print(m)

        w.samples[m,] <- sample.w(t.l, t.u, right.censored, int.censored,
                Z.samples[m-1,,],
                X,
                beta.samples[m-1,,],
                sigma.samples[m-1,])
        #w.samples[m,] <- sample.w.2(t.l, t.u, right.censored, int.censored,
        #        X,
        #        beta.samples[m-1,,],
        #        rho.samples[m-1,,],
        #        sigma.samples[m-1,])

        Z.samples[m,,] <- sample.z(w.samples[m,], X,
                beta.samples[m-1,,], 
                rho.samples[m-1,,],
                sigma.samples[m-1,])
        #print (Z.samples[m,,])
        
        if (any(colSums(Z.samples[m,,])==0)) {
            browser()
        }

        # sample.beta uses log of time
        Y <- log(w.samples[m,])
        #out <- sample.beta.2(Y, X, Z.samples[m,,])
        out <- sample.beta(Y, X, Z.samples[m,,], n0, S0, b0, B0)
        beta.samples[m,,] <- out$beta
        sigma.samples[m,] <- out$sigma
        
        out <- sample.rho(X, Z.samples[m,,], rho.samples[m-1,,], tune)
        rho.samples[m,,] <- out$rho
        accepts <- accepts + out$accept

    }
    acc.rate <- (accepts*100)/m
    list(w.samples=w.samples,Z.samples=Z.samples,sigma.samples=sigma.samples,
            beta.samples=beta.samples,rho.samples=rho.samples,acc.rate=acc.rate)
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

sample.w.2 <- function(t.l, t.u, right.censored, int.censored, X,
        .beta, .rho, .sigma) {
    # convert from 0,1 to T,F
    .w <- matrix(nrow=1,ncol=nrow(X))
    for (i in 1:nrow(t.l)) {
        if (int.censored[i]==1) {
            Fl <- dF.i(t.l[i], X[i,], .beta, .rho, .sigma)
            Fu <- dF.i(t.u[i], X[i,], .beta, .rho, .sigma)
            Fw <- runif(1, Fl, Fu) 
            f.tmp <- function(w) {
                Fw-dF.i(w, X[i,], .beta, .rho, .sigma)
            }
            t. <- seq(0, 1, 0.01) 
            ft. <- apply(as.matrix(t.), 1,f.tmp)
            .w[i,] <- uniroot(f.tmp, lower=t.l[i], upper=t.u[i])
        } else if (right.censored[i]==1) {
            Fl <- dF.i(t.l[i], X[i,], .beta, .rho, .sigma)
            Fw <- function(w) {
                runif(1, Fl, 1)-dF.i(w, X[i,], .beta, .rho, .sigma)
            }
            print(sprintf("%f %f %f %f" , t.l[i], 1, Fw(t.l[i]), Fw(1)))
            .w[i] <- uniroot(Fw, lower=t.l[i], upper=999)
        } else {
            .w[i] <- t.l[i]
        }
    }
    .w
}

# return nxJ matrix Z where Z[i,j]=1 indicates obs. i comes from component j
sample.z <- function(w, X, .beta, .rho, .sigma) {
    # .p: nxj matrix where .p(i,j) = p(x.i, rho.j)
    xtr <- tcrossprod(X, .rho)
    # logs to (try to) avoid overflow
    .p <- exp(xtr - repmat(log(rowSums(exp(xtr))), 1, nrow(.rho)))
    # .f: nxj matrix where .f(i,j) = .f_j(x.i, .beta.j) 

    .f <- dlnorm(repmat(w, 1, ncol(.p)), tcrossprod(X,.beta), repmat(matrix(.sigma,nrow=1), nrow(X), 1))
    # h: nxj matrix where

    #hist(w, 20)
    #print(w)
    #print (.f)

    h <- matrix(1, nrow=nrow(.p), ncol=ncol(.p))
    for (j in 1:ncol(.p)) {
        for (j2 in 1:ncol(.p)) {
            if (j2 != j) {
                h[,j] <- h[,j] + (.p[,j2]*.f[,j2])/(.p[,j]*.f[,j])
                #(.p[,j2]*.f[,j2])
                #(.p[,j]*.f[,j])
                #print((.p[,j2]*.f[,j2])/(.p[,j]*.f[,j]))
                #(1 + (.p[,j2]*.f[,j2])/(.p[,j]*.f[,j]))^(-1)
            }
        }
    }
    h <- 1/h

    if (any(is.nan(h))) {
        browser()
        #print("nan")
        #print(xtr)
        ##print(rowSums(exp(xtr)))
        #print(.p)
        #print(.f)
    }

    Z <- matrix(nrow=nrow(X), ncol=nrow(.beta))
    for (i in 1:nrow(X)) {
        Z[i,] <- rmultinom(1, 1, h[i,])
    }
    Z
}

# Y: log(w)
# z: column of memberships for jth component
sample.beta <- function(Y,X,Z,n0,S0,b0,B0) {

    #J x p
    out.beta <- matrix(NA, nrow=ncol(Z), ncol=ncol(X))
    #1 x J
    out.sigma <- matrix(NA, nrow=1, ncol=ncol(Z))
    for (j in 1:nrow(out.beta)) {
        #out <- sample.beta.j(Y, X, Z.samples[m,,j], sigma.samples[m-1,j]^2, n0, S0, b0, B0)  
        #sigma.samples[m,j] <- sqrt(out$sig2j)
        #beta.samples[m,j,] <- t(out$beta.j)
        X.s<- matrix(X[Z[,j]==1,],nrow=sum(Z[,j]) ,ncol=ncol(X))
        Y.s<- Y[Z[,j]==1]
        .n <-nrow(X.s)
        p <- ncol(X.s)
        b0<-matrix(0,ncol(X.s),1)
        B0<-1000*diag(ncol(X.s))
        B1<- ginv(ginv(B0) + t(X.s) %*% X.s)
        b1<-B1%*%(ginv(B0)%*%b0 + t(X.s)%*%Y.s)
        beta.h<- ginv(t(X.s)%*% X.s ) %*% t(X.s) %*% Y.s
        n1<-n0 + .n
        #beta.j<- t(rmvt(1,B1,n1)) + b1
        beta.j<- t(rmvt(1,B1,n1)) + b1
        S2 <- t(Y.s)%*% ( diag(.n) - X.s%*%ginv(t(X.s)%*% X.s)%*%t(X.s) )%*% Y.s / (.n - ncol(X.s))
        n1S1<-n0*S0 + (.n-p)*S2[1,1] + t(beta.h - b0)%*% ginv( B0 + ginv(t(X.s)%*%X.s) ) %*% (beta.h -b0)
        out.sigma[j] <- sqrt(1/ rgamma(1,n1/2,n1S1/2))
        out.beta[j,] <- beta.j 
        #while (is.nan(sigma.j))  {
        #    sigma.j <- 1/ rgamma(1,n1/2,n1S1/2)
        #}
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

# .rho: rho matrix Jxp
# tune: variance of proposal distrib
#  Z: nxJ membership matrix
# returns 1
sample.rho <- function(X, Z, .rho, tune) {
    # log posteriors of rows of a rho matrix
    # as a vector of length J
    log.post <- function(r) {
        # TODO flatness for log.pri
        log.pri <- dmvnorm(r, rep(0, ncol(X)), 100*diag(ncol(X)),log=TRUE)
        exr <- exp(tcrossprod(X, r))
        exr <- exr/repmat(rowSums(exr), 1, ncol(exr))
        log.lik <- colSums(log(exr))
        log.lik + log.pri
    }

    cand.rho <- matrix(NA, nrow=nrow(.rho), ncol=ncol(.rho))
    # set to zero for identifiability
    cand.rho[1,] <- 0
    for (j in 2:nrow(.rho)) {
        cand.rho[j,] <- rmvnorm(1, .rho[j,], tune*diag(ncol(.rho)))
    }
    # vector of ratios J
    ratio <- exp(log.post(cand.rho) - log.post(.rho)) 

    accept <- rep(0, nrow(.rho))
    for (j in 2:nrow(.rho)) {
        if (runif(1) < ratio[j]) {
            accept[j] <- 1
            .rho[j,] <- cand.rho[j,] 
        } 
    }

    exr <- exp(tcrossprod(X, .rho))
    exr <- exr/repmat(rowSums(exr), 1, ncol(exr))
    matplot(y=exr)

    list(rho=.rho,accept=accept)
}

# matlab-like
repmat <- function(a,n,m) { kronecker(matrix(1, n, m), a) }

# F(t_i | x_i, theta)
dF.i <- function(t.i, x.i, .beta, .rho, .sigma) {
    .p <- exp(tcrossprod(x.i, .rho))
    .p <- .p/repmat(rowSums(.p), 1, ncol(.p))
    sum(.p*plnorm(t.i, meanlog=tcrossprod(x.i,.beta), sdlog=.sigma))
}

df.i <- function(t.i, x.i, .beta, .rho, .sigma) {
    .p <- exp(tcrossprod(x.i, .rho))
    .p <- .p/repmat(rowSums(.p), 1, ncol(.p))
    sum(.p*dlnorm(t.i, meanlog=tcrossprod(x.i,.beta), sdlog=.sigma))
}

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

test.sample.z <- function() {
    d <- sim.data() 
    # use uncensored data
    w <- d$t.l
    z <- sample.z(w, d$X, d$beta, d$rho, d$sigma)
    print(d$z)
    print(z)
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

test.sample.beta.j <- function() {
    z<-sample(c(0,1),size=1000,prob=c(0.4, 0.6) ,replace=T)
    X<-cbind(1,rnorm(1000))
    b<-matrix(c(2,5),2,1)
    Y<-X %*% b + rnorm(1000)

    sample.beta.j(Y,X,z,sig2j=3,n0=1,S0=2,b0=matrix(0,2,1),B0=100*diag(2))
}

test.sample.rho <- function() {
    d <- sim.data() 
    # use uncensored data
    w <- d$t.l
    tune <- 1
    Z <- sample.z(w, d$X, d$beta, d$rho, d$sigma)
    .rho <- matrix(c(2, 3), ncol=1) 
    for (i in 1:100) {
        out <- sample.rho(d$X, Z, .rho, tune)
        .rho <- out$rho
        #print(out$accept)
        #print(.rho)
    }
}

test.sampler.1 <- function() {
    n0 <- 1
    S0 <- diag(100)
    b0 <- 1
    B0 <- 1
    tune <- 1
    M <- 5000
    J <- 2

    sim <- sim.data()
    X <- cbind(rep(1,nrow(sim$X)), sim$X)

    rc <- sim$right.censored
    ic <- sim$int.censored

    out <- moaftme.sampler(sim$t.l, sim$t.u, rc, ic, X, J, M, 
        n0=1, S0=2, b0=matrix(0,2,1), B0=100*diag(2), tune)
}

test.sampler.2 <- function() {
    n0 <- 1
    S0 <- diag(100)
    b0 <- 1
    B0 <- 1
    tune <- 50
    M <- 3000
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

moaftme.MCEM <- function(t.l, t.u, right.censored, int.censored, X, J, M) {

    n <- nrow(X)
    p <- ncol(X)
    .beta <- rnorm(p*J, 0, 1e-6)
    .rho  <- rnorm(p*J, 0, 1e-6)
    .sigma <- 1/rgamma(J, 1/2, 1/2)

    wm <- matrix(NA, nrow=M, ncol=p)
    for (m in 1:M) {

    }

}

#test.sample.beta.j()
#test.sampler.1()
#out <- test.sampler.2()

if (FALSE) {
    ts.plot(out$beta[,1,1])
    plot(density(out$beta[,1,1]))
    ts.plot(out$beta[,1,2])
    ts.plot(out$beta[,2,1])
    ts.plot(out$beta[,2,2])
    ts.plot(out$sigma[,1])
    plot(density(out$sigma[,1]))
    ts.plot(out$sigma[,2])
    ts.plot(out$sigma[,3])
    ts.plot(out$rho[,2,1])
    ts.plot(out$rho[,3,1])
    matplot(cbind(out$beta[,1,1], 10*out$sigma[,1]), type='l')
    matplot(out$beta[,1,1:2], type='l')
    matplot(out$rho[,1,1:2], type='l')
    matplot(out$rho[,2,1:2], type='l')
}

#test.sample.rho()
#test.sample.w()
#test.sample.z()
#sim.data(TRUE)
#package.skeleton(name="moaftme", namespace=TRUE)
test.sampler.1()
#test.sampler.2()
