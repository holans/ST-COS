gibbs.stcos <- function(Z, S, sigma2.eps, C.inv, H, R = 1000,
	report.period = R+1, burn = 0, thin = 1)
{
	Vinv <- 1 / sig2eps
	HpVinv <- t(H) %*% diag(Vinv) 
	HpinvVH <- HpVinv %*% H
	eig.HpinvVH <- eigen(HpinvVH)

	r <- ncol(S)
	n <- nrow(Z)
	n_mu <- nrow(HpVinv)

	mu_B.hist <- matrix(NA, R, n_mu)
	xi.hist <- matrix(NA, R, n)
	eta.hist <- matrix(NA, R, r)
	sig2mu.hist <- matrix(NA, R, 1)
	sig2xi.hist <- matrix(NA, R, 1)
	sig2K.hist <- matrix(NA, R, 1)
	Y.hist <- matrix(NA, R, n)

	# Initial values
	eta <- rnorm(r)
	xi <- rep(1, n)
	sig2mu <- 1
	sig2xi <- 1
	sig2K <- 1

	for (j in 1:r) {
		SpinvV[j,] <- S[,j] / sigma2.eps
	}

	logger("Gibbs sampler starting\n")

	for (tt in 1:R) {
		# Full Conditional for mu_B
		# sparse matrix inverse
		PostLam <- diag(1 / (eig.HpinvVH$values + 1/sig2mu))
		V.mu_B <- eig.HpinvVH$vectors * PostLam * t(eig.HpinvVH$vectors)
		V.mu_B.half <- eig.HpinvVH$vectors * sqrt(PostLam)
		mean.mu_B <- eig.HpinvVH$vectors * PostLam * V.mu_B * HpVinv * (Z - S*eta - xi)
		mu_B <- mean.mu_B + V.mu_B.half * rnorm(n)

		# Full Conditional for sig2mu
		shapesig2mu <- 0.5 * t(mu_B) %*% mu_B
		sig2mu <- 1/rgamma(1, n_mu/2 + 1, 1 / (0.000001 + shapesig2mu))

		# Full Conditional for eta
		zresid <- Z - H*mu_B - xi
		V.eta <- solve((1/sig2K) * C.inv + SpinvV*S)
		mean.eta <- V.eta * SpinvV * zresid
		eta <- mean.eta + chol(V.eta) * rnorm(r)

		# Full Conditional for xi
		# sparse matrix inverse
		Sigma.xi.inv <- Vinv + 1 / sig2xi
		Sigma.xi <- 1 / Sigma.xi.inv
		mean.xi <- Sigma.xi * Vinv * (Z - S*eta - H*mu_B)
		xi <- mean.xi + sqrt(Sigma.xi) * rnorm(n)

		# Full Conditional for sig2K
		scale <- 0.5 * t(eta) %*% C.inv %*% eta
		sig2K <- 1 / rgamma(1, r/2 + 1, 1 / (0.000001 + scale))

		# Update Y
		Y = S*eta + H*mu_B + xi

		if (tt %% report.period == 0) {
			print(tt)
		}

		# Save history
		if ((tt > burn) & (tt %% thin == 0) {
			mu_B[tt,] <- mu_B
			eta[tt,] <- eta
			xi[tt,] <- xi
			sig2mu[tt] <- sig2mu
			sig2K[tt] <- sig2K
			Y.hist[tt,] <- Y
		}

		logger("Gibbs sampler finished\n")
	}

	list(mu_B.hist = mu_B.hist,
		xi.hist = xi.hist,
		eta.hist = eta.hist,
		sig2mu.hist = sig2mu.hist,
		sig2xi.hist = sig2xi.hist,
		sig2K.hist = sig2K.hist,
		Y.hist = Y.hist
	)
}

