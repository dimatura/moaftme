
source('moaftme.R')

test.sampler.1 <- function() {
    n0 <- 0.01
    S0 <- 0.01
    b0 <- matrix(0, nrow=1, ncol=2)
    B0 <- diag(0.01, 2)
    tune <- 0.50
    M <- 1000
    J <- 3
    gamma0 <- 10

    sim <- sim.data()
    X <- cbind(rep(1,nrow(sim$X)), sim$X)

    rc <- sim$right.censored
    ic <- sim$int.censored

    out <- moaftme.sampler(sim$t.l, sim$t.u, rc, ic, X, J, M, 
        n0=n0, S0=S0, b0=b0, B0=B0, gamma0=gamma0, tune)
    out
}

test.sampler.2 <- function() {
    n0 <- 0.01
    S0 <- 0.01
    b0 <- matrix(0, nrow=1, ncol=2)
    B0 <- diag(0.1, 2)
    tune <- 8
    M <- 10000
    J <- 2
    gamma0 <- 30

    d <- read.table('flourbeetle.txt', header=TRUE)

    X <- model.matrix(~ X1, d)

    rc <- d$rc
    ic <- d$ic
    t.l <- matrix(d$t.l,ncol=1)
    t.u <- matrix(d$t.u,ncol=1)

    out <- moaftme.sampler(t.l, t.u, rc, ic, X, J, M, 
        n0=n0, S0=S0, b0=b0, B0=B0, gamma0=gamma0, tune)
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
    ts.plot(out$sigma[,2])
    plot(density(out$sigma[,1]))
    ts.plot(out$rho[,2,1])
    ts.plot(out$rho[,2,2])
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
#out <- test.sampler.1()
out <- test.sampler.2()

