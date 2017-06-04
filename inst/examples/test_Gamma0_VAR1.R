library(stcos)
library(mvtnorm)

m <- 3109

A <- diag(0.5, m)
diag(A[-m,-1]) <- 0.25
diag(A[-1,-m]) <- 0.25

Sigma <- diag(1, m)
diag(Sigma[-m,-1]) <- 0.5
diag(Sigma[-1,-m]) <- 0.5

# Avoid Kronecker product
eig <- eigen(A)
V <- eig$vectors
V.inv <- solve(eig$vectors)
C <- V.inv %*% Sigma %*% t(V.inv)
e <- 1 / (1 - eig$values %x% eig$values) * as.numeric(C)
res0 <- V %*% matrix(e, m, m) %*% t(V)

if (FALSE) {
	V <- eig$vectors
	D <- Diagonal(x = eig$values)
	D.tilde <- Diagonal(x = 1 - eig$values)

	eig.comp <- eigen(diag(m) - A)
	eig.comp$values
	eig.comp$vectors
	V %*% D.tilde %*% solve(V)

	eig.kron <- eigen(A %x% A)
	eig.kron$values
	eig.kron$vectors
	idx <- order(eig.kron$values)
	
	eig.kron.comp <- eigen(diag(m^2) - A %x% A)
	eig.kron.comp$values
	eig.kron.comp$vectors
	
	1 - rev(eig.kron$values) - eig.kron.comp$values
	eig.kron$vectors[,idx] - eig.kron.comp$vectors

	V %x% V

	VV <- matrix(NA, m^2, m^2)
	for (j in 1:m) {
		for (k in 1:m) {
			VV[,k + m*(j-1)] <- V[,j] %x% V[,k]
		}
	}
	VV %*% Diagonal(x = 1 - (eig$values %x% eig$values)) %*% solve(VV)

	cbind(A %*% V[,1], eig$values[1] * V[,1])
	cbind((diag(m) - A) %*% V[,1], (1 - eig$values[1]) * V[,1])

	norm(diag(m) - A - V %*% D.tilde %*% solve(V), "F")
	res2 <- solve(diag(nrow=m^2) - A %x% A, as.numeric(Sigma))

	norm(A %x% A - (V %x% V) %*% (D %x% D) %*% (solve(V) %x% solve(V)), type = "F")
	norm((Diagonal(m^2) - A %x% A) - ((V %x% V) %*% Diagonal(x = 1 - (eig$values %x% eig$values)) %*% (solve(V) %x% solve(V))), type = "F")
}

system.time(out <- solve_Gamma0(A, Sigma))

res <- solve(diag(nrow=m^2, x=1) - A %x% A, as.numeric(Sigma))
matrix(res, m, m)

cbind(out$x, res)
sum(abs((out$x - res)))

system.time(res.mc <- sample_Gamma0_VAR1(A, Sigma, R = 1000000, burn = 10000))
cbind(out$x, as.numeric(res.mc))
sum(abs((out$x - as.numeric(res.mc))))

system.time(res.mc <- sample_Gamma0_VAR1_v2(A, Sigma, R = 1000000, burn = 10000))
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
