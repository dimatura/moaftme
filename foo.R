
dyn.load('moaftmec.so')

a <- c(3,2,1)
b <- c(1,9,8)
c <- rep(NA, 3)

out <- .C("foo", a=as.double(a), b=as.double(b), len=as.integer(length(a)), c=as.double(c), NAOK=TRUE)

A <- array(1:24, c(2,3,4))
out <- .C("foo2", A=as.double(A), as.integer(dim(A)[1]), as.integer(dim(A)[2]), as.integer(dim(A)[3]))
