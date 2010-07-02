
library(mvtnorm)
library(MASS)

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

# matlab-like
repmat <- function(a,n,m) { kronecker(matrix(1, n, m), a) }

# F(t_i | x_i, theta)
dF.i <- function(t.i, x.i, .beta, .rho, .sigma) {
    #browser()
    .p <- exp(x.i %*% t(.rho)) 
    .p <- .p / rowSums(.p)
    #sum(.p*pnorm((log(t.i)-x.i%*%.beta)/.sigma))
    sum(.p*plnorm(t.i, meanlog=x.i%*%t(.beta), sdlog=.sigma))
}

# F(t_i | x_i, theta)
df.i <- function(t.i, x.i, .beta, .rho, .sigma) {
    #browser()
    .p <- exp(x.i %*% t(.rho)) 
    .p <- .p / rowSums(.p)
    #sum(.p*pnorm((log(t.i)-x.i%*%.beta)/.sigma))
    sum(.p*dlnorm(t.i, meanlog=x.i%*%t(.beta), sdlog=.sigma))
}

# return .w, 1xp matrix of pseudo-observations
sample.w <- function(t.l, t.u, right.censored, int.censored, Z, X,
        .beta, .sigma) {
    # convert from 0,1 to T,F
    #browser()
    Z <- Z > 0.5
    .w <- matrix(nrow=1,ncol=nrow(X))
    for (i in 1:nrow(t.l)) {
        if (int.censored[i]==1) {
            #xtb <- crossprod(X[i,], .beta[Z[i,],])
            Fl <- plnorm(t.l[i], X[i,] %*% t(.beta[Z[i,],]), .sigma[Z[i,]])
            Fu <- plnorm(t.u[i], X[i,] %*% t(.beta[Z[i,],]), .sigma[Z[i,]])
            .w[i] <- qlnorm(runif(1, Fl, Fu), X[i,] %*% t(.beta[Z[i,],]), .sigma[Z[i,]])
        } else if (right.censored[i]==1) {
            Fl <- plnorm(t.l[i], X[i,] %*% t(.beta[Z[i,],]), .sigma[Z[i,]])
            .w[i] <- qlnorm(runif(1, Fl, 1 ), X[i,] %*% t(.beta[Z[i,],]), .sigma[Z[i,]])
        } else {
            .w[i] <- t.l[i]
        }
    }
    .w
}

# return nxJ matrix Z where Z[i,j]=1 indicates obs. i comes from component j
sample.z <- function(w, X, .beta, .rho, .sigma) {
    # .p: nxj matrix where .p(i,j) = p(x.i, rho.j)
    xtr <- X %*% t(.rho)
    # logs to (try to) avoid overflow
    .p <- exp(xtr - repmat(log(rowSums(exp(xtr))), 1, length(.rho)))
    # .f: nxj matrix where .f(i,j) = .f_j(x.i, .beta.j) 
    .f <- dlnorm(repmat(w, 1, ncol(.p)),meanlog=X %*% t(.beta),sdlog=repmat(matrix(.sigma,nrow=1), nrow(X), 1))
    # h: nxj matrix where
    # h[i,j] = p_j(x_i,rho_j)*f_j(w_i|x_i,beta_j)/(sum_{j=1}^J(p_j(x_i,rho_j)*f_j(w_i|x_i,beta_j)))
    h <- (.p * .f) / repmat(rowSums(.p * .f), 1, length(.rho))
    Z <- matrix(nrow=nrow(X), ncol=length(.rho))
    for (i in 1:nrow(X)) {
        Z[i,] <- rmultinom(1, 1, h[i,])
    }
    Z
}

# Y: log(w)
# z: column of memberships for jth component
sample.beta.j <- function(Y,X,z,sig2j,n0,S0,b0,B0) {
    X.s<- matrix(X[z==1,],ncol=X)
    Y.s<- Y[z==1]
    .n <-nrow(X.s)
    p <- ncol(X.s)
    b0<-matrix(0,ncol(X.s),1)
    B0<-1000*diag(ncol(X.s))
    B1<- ginv(ginv(B0) + t(X.s) %*% X.s)
    b1<-B1%*%(ginv(B0)%*%b0 + t(X.s)%*%Y.s)
    beta.j<-rmvnorm(1,b1,sig2j*B1)

    beta.h<- ginv(t(X.s)%*% X.s ) %*% t(X.s) %*% Y.s
    n1<-n0 + .n
    S2 <- t(Y.s)%*% ( diag(.n) - X.s%*%ginv(t(X.s)%*% X.s)%*%t(X.s) )%*% Y.s / (.n - ncol(X.s))
    n1S1<-n0*S0 + (.n-p)*S2[1,1] + t(beta.h - b0)%*% ginv( B0 + ginv(t(X.s)%*%X.s) ) %*% (beta.h -b0)
    sigma.j <- 1/ rgamma(1,n1/2,n1S1/2)
    list("beta.j" = beta.j,"sig2j"=sigma.j)
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
        log.pri <- dmvnorm(r, 100*diag(ncol(X)),log=TRUE)
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
    list(rho=.rho,accept=accept)
}

sampler <- function(t.l, t.u, right.censored, int.censored, X, J, M, 
        n0, S0, b0, B0, tune) {
    n <- nrow(X)
    p <- ncol(X)

    # samples
    w.samples <- matrix(nrow=M, ncol=n)
    Z.samples <- array(NA, c(M, n, J))
    sigma.samples <- matrix(nrow=M, ncol=J)
    beta.samples <- array(NA, c(M, J, p))
    rho.samples <- array(NA, c(M, J, p))

    # initial sampling
    beta.samples[1,,] <- matrix(rnorm(p*J, mean=0, sd=10), ncol=p)
    rho.samples[1,,] <- matrix(rnorm(p*J, mean=0, sd=10), ncol=p)
    # set first col to zero for identifiability
    rho.samples[1,1,] <- 0
    sigma.samples[1,] <- matrix(1/rgamma(J, n0/2, n0*S0/2), nrow=1)

    # main loop
    for (m in 2:M) {

        w.samples[m,] <- sample.w(t.l, t.u, right.censored, int.censored,
                Z.samples[m-1,,],
                X,
                beta.samples[m-1,,],
                sigma.samples[m-1,])

        Z.samples[m,,] <- sample.z(w.samples[m,], X,
                beta.samples[m-1,,], 
                rho.samples[m-1,,],
                sigma.samples[m-1,])
        # sample.beta uses log of time
        Y <- log(w.samples[m,])
        for (j in 1:J) {
            out <- sample.beta.j(Y, X, Z.samples[m,,j], sigma.samples[m-1,j]^2, n0, S0, b0, B0)  
            sigma.samples[m,j] <- sqrt(out$sig2j)
            beta.samples[m,j,] <- t(out$beta.j)
        }
        rho.samples[m,,] <- sample.rho(X, Z.samples[m,,], .rho[m-1,,], tune)
    }
    list(w.samples=w,Z.samples=Z.samples,sigma.samples=sigma.samples,
            beta.samples=beta.samples,rho.samples=rho.samples)
}

sim.data <- function(plot=FALSE) {
    n <- 100
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
    #J <- 2
    #M <- 10
    #sampler(t., t.l, t.u, t.p, right.censored, int.censored, X, J, M)
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

test.sampler.1 <- function() {
    n0 <- 1
    S0 <- diag(100)
    b0 <- 1
    B0 <- 1
    tune <- 1
    M <- 10
    J <- 2

    sim <- sim.data()

    rc <- sim$right.censored
    ic <- sim$int.censored

    out <- sampler(sim$t.l, sim$t.u, rc, ic, d$X, J, M, 
        n0=1, S0=2, b0=matrix(0,2,1), B0=100*diag(2), tune)
    print(out)
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

#test.sample.beta.j()
test.sampler.1()
#test.sample.rho()
#test.sample.w()
#test.sample.z()
#sim.data(TRUE)
#package.skeleton(name="moaftme", namespace=TRUE)
#test.sampler.1()
