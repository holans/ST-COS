gibbs.stcos <- function(Z, S, sig2eps, C.inv, H, R,
	report.period = R+1, burn = 0, thin = 1)
{
	Vinv <- 1 / sig2eps
	HpVinv <- t(H) %*% Diagonal(n = length(Vinv), x = Vinv)
	HpinvVH <- HpVinv %*% H

	logger("Begin computing eigenvalues/vectors of HpinvVH\n")
	eig.HpinvVH <- eigen(HpinvVH, symmetric = TRUE)
	logger("Finished computing eigenvalues/vectors of HpinvVH\n")

	r <- ncol(S)
	n <- nrow(Z)
	n_mu <- nrow(HpVinv)

	tt.keep <- 0
	R.keep <- ceiling((R - burn) / thin)
	mu_B.hist <- matrix(NA, R.keep, n_mu)
	xi.hist <- matrix(NA, R.keep, n)
	eta.hist <- matrix(NA, R.keep, r)
	sig2mu.hist <- matrix(NA, R.keep, 1)
	sig2xi.hist <- matrix(NA, R.keep, 1)
	sig2K.hist <- matrix(NA, R.keep, 1)
	Y.hist <- matrix(NA, R.keep, n)

	# Initial values
	mu_B <- rnorm(n_mu)
	eta <- rnorm(r)
	xi <- numeric(n)
	sig2mu <- 1
	sig2xi <- 1
	sig2K <- 1
	Y <- S %*% eta + H %*% mu_B + xi

	logger("Begin computing SpinvV\n")
	SpinvV <- matrix(NA, r, n)
	for (j in 1:r) {
		SpinvV[j,] <- S[,j] / sig2eps
	}
	SpinvVS <- SpinvV %*% S
	logger("Finished computing SpinvV\n")

	logger("Begin Gibbs sampler\n")

	for (tt in 1:R) {
		if (tt %% report.period == 0) {
			logger("Begin iteration %d, using %0.2f GB RAM\n", tt, mem_used() / 2^30)
		}

		# Full Conditional for mu_B, using sparse matrix inverse
		PostLam <- 1 / (eig.HpinvVH$values + 1/sig2mu)
		V.mu_B <- eig.HpinvVH$vectors %*% (PostLam * t(eig.HpinvVH$vectors))
		mean.mu_B <- V.mu_B %*% (HpVinv %*% (Z - S %*% eta - xi))
		mu_B <- mean.mu_B + (eig.HpinvVH$vectors %*% (sqrt(PostLam) * rnorm(n_mu)))

		# Full Conditional for sig2mu
		scale <- as.numeric(0.5 * t(mu_B) %*% mu_B)
		sig2mu <- 1 / rgamma(1, n_mu/2 + 1, 0.000001 + scale)

		# Full Conditional for eta
		zresid <- Z - H %*% mu_B - xi
		V.eta <- solve(as.matrix((1/sig2K * C.inv) + SpinvVS))
		mean.eta <- V.eta %*% (SpinvV %*% zresid)
		eta <- mean.eta + chol(V.eta) %*% rnorm(r)

		# Full Conditional for xi, using sparse matrix inverse
		Sigma.xi.inv <- Vinv + 1/sig2xi
		Sigma.xi <- 1 / Sigma.xi.inv
		mean.xi <- Sigma.xi * Vinv * (Z - S %*% eta - H %*% mu_B)
		xi <- mean.xi + sqrt(Sigma.xi) * rnorm(n)

		# Full Conditional sig2xi
		scale <- as.numeric(0.5 * t(xi) %*% xi)
		sig2xi <- 1 / rgamma(1, n/2 + 1, 0.000001 + scale)

		# Full Conditional for sig2K
		scale <- as.numeric(0.5 * t(eta) %*% C.inv %*% eta)
		sig2K <- 1 / rgamma(1, r/2 + 1, 0.000001 + scale)

		# Update Y
		Y <- S %*% eta + H %*% mu_B + xi

		# Save history
		# TBD: This is very expensive to keep in memory. We should be able to
		# transparently store it on disk, if the user wishes
		if ((tt > burn) & (tt %% thin == 0)) {
			tt.keep <- tt.keep + 1
			mu_B.hist[tt.keep,] <- as.numeric(mu_B)
			eta.hist[tt.keep,] <- as.numeric(eta)
			xi.hist[tt.keep,] <- as.numeric(xi)
			sig2mu.hist[tt.keep] <- sig2mu
			sig2xi.hist[tt.keep] <- sig2xi
			sig2K.hist[tt.keep] <- sig2K
			Y.hist[tt.keep,] <- as.numeric(Y)
		}
	}

	logger("Finished Gibbs sampler\n")

	list(mu_B.hist = mu_B.hist,
		xi.hist = xi.hist,
		eta.hist = eta.hist,
		sig2mu.hist = sig2mu.hist,
		sig2xi.hist = sig2xi.hist,
		sig2K.hist = sig2K.hist,
		Y.hist = Y.hist
	)
}

