library(stcos)
library(mvtnorm)

m <- 50

A <- diag(0.5, m)
diag(A[-m,-1]) <- 0.25
diag(A[-1,-m]) <- 0.25
A[1,m] <- 0.25

Sigma <- diag(1, m)
diag(Sigma[-m,-1]) <- 0.5
diag(Sigma[-1,-m]) <- 0.5

# Avoid Kronecker product
system.time(res0 <- covVAR1(A, Sigma, lag_max = 0))

system.time(res <- solve(diag(nrow=m^2, x=1) - A %x% A, as.numeric(Sigma)))
res1 <- matrix(res, m, m)
sum(abs(res0[,,1] - res1))

system.time(res.mc <- sample_Gamma0_VAR1(A, Sigma, R = 1000000, burn = 10000))
sum(abs((res0[,,1] - res.mc)))
