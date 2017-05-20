library(stcos)

m <- 500

A <- diag(0.5, m)
diag(A[-m,-1]) <- 0.25
diag(A[-1,-m]) <- 0.25

Sigma <- diag(1, m)

system.time(out <- solve_Gamma0(A, Sigma))

res <- solve(diag(nrow=m^2, x=1) - A %x% A, as.numeric(Sigma))

cbind(out$x, res)
sum(abs((out$x - res)))
