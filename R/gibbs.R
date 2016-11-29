function [Y, S, eta, mu_B, xi, sig2mu, sig2xi, sig2K] = Met_spt_COS_xiismean(Z, T, X, S, sig2eps, Kinv, LamHpVinvH, EigHpVinvV, HpVinv, H)

gibbs.stcos <- function(Z, S, sigma2.eps, Kinv, H, R = 1000, report.period = R+1)
{
	Vinv <- 1 / sig2eps
	HpVinv <- t(H) %*% diag(Vinv) 
	HpinvVH <- HpVinv %*% H
	eig <- eigen(HpinvVH)
	EigHinvVHp <- eig$vectors
	LamHpinvVH <- diag(eig$values)

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

	for (j in 1:r) {
		SpinvV[j,] <- S[,j] / sigma2.eps
	}

	logger("Gibbs sampler starting\n")
	
	for (tt in 1:R) {
		# Full Conditional for mu_B
		# sparse matrix inverse
		PostLaminv <- diag(LamHpVinvH) + 1/sig2mu
		PostLam <- diag(1 / PostLaminv)
		V.mu_B <- EigHpVinvV * PostLam * t(EigHpVinvV)
		V.mu_B.half <- EigHpVinvV * sqrt(PostLam)
		mean.mu_B <- EigHpVinvV * PostLam * V.mu_B * HpVinv * (Z - S*eta - xi)
		mu_B <- mean.mu_B + V.mu_B.half * rnorm(n)

		# Full Conditional for sig2mu
		shapesig2mu <- 0.5 * t(mu_B) %*% mu_B
		sig2mu <- 1/gamrnd(n_mu/2 + 1, 1 / (0.000001 + shapesig2mu))

		# Full Conditional for eta
		zresid <- Z - H*mu_B - xi
		Omega.eta <- solve((1/sig2K) * Kinv + SpinvV*S)
		mean.eta <- Omega.eta * SpinvV * zresid
		eta <- mean.eta + chol(Omega.eta) * rnorm(r)

		# Full Conditional for xi
		# sparse matrix inverse
		Sigma.xi.inv = Vinv + 1 / sig2xi
		Sigma.xi = 1 / Sigma.xi.inv
		mean.xi = PostLam * Vinv * (Z - S*eta - H*mu_B)
		xi <- mean.xi + sqrt(Sigma.xi) * rnorm(n)

		# Full Conditional for sig2K
		scale = 0.5 * t(eta) %*% Kinv %*% eta
		sig2K = 1 / gamrnd(r/2 + 1, 1 / (0.000001 + scale))

		# Update Y
		Y = S*eta + H*mu_B + xi

		if (tt %% report.period == 0) {
			print(tt)
		}

		# Save history
		mu_B[tt,] <- mu_B
		eta[tt,] <- eta
		xi[tt,] <- xi
		sig2mu[tt] <- sig2mu
		sig2K[tt] <- sig2K

		logger("Gibbs sampler finished\n")
	}
}

