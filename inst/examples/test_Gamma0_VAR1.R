library(stcos)
library(mvtnorm)

m <- 100

A <- diag(0.5, m)
diag(A[-m,-1]) <- 0.25
diag(A[-1,-m]) <- 0.25

Sigma <- diag(1, m)
diag(Sigma[-m,-1]) <- 0.5
diag(Sigma[-1,-m]) <- 0.5

system.time(out <- solve_Gamma0(A, Sigma))

res <- solve(diag(nrow=m^2, x=1) - A %x% A, as.numeric(Sigma))

cbind(out$x, res)
sum(abs((out$x - res)))

system.time(res.mc <- sample_Gamma0_VAR1(A, Sigma, R = 1000000, burn = 10000))
cbind(out$x, as.numeric(res.mc))
sum(abs((out$x - as.numeric(res.mc))))

sample.Gamma0.VAR1 <- function(A, Sigma, R, burn = 10000)
{
	m <- nrow(Sigma)
	Sigma.half <- chol(Sigma)
	y.t.minus1 <- Sigma.half %*% rnorm(m)
	Gamma0 <- matrix(0, m, m)
	
	for (t in 1:burn) {
		z <- Sigma.half %*% rnorm(m)
		y.t <- A %*% y.t.minus1 + z
		y.t.minus1 <- y.t
	}

	for (t in 1:R) {
		z <- Sigma.half %*% rnorm(m)
		y.t <- A %*% y.t.minus1 + z
		y.t.minus1 <- y.t
		Gamma0 <- Gamma0 + (y.t %*% t(y.t)) / R
	}

	return(Gamma0)
}

system.time(res.mc <- sample.Gamma0.VAR1(A, Sigma, R = 10000000, burn = 10000))
cbind(out$x, as.numeric(res.mc))
sum(abs((out$x - as.numeric(res.mc))))
