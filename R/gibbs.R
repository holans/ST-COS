gibbs.stcos <- function(Z, S, sig2eps, C.inv, H, R,
	report.period = R+1, burn = 0, thin = 1,
	init = NULL, fixed = NULL)
{
	stopifnot(R > burn)

	Vinv <- 1 / sig2eps
	HpVinv <- t(H) %*% Diagonal(n = length(Vinv), x = Vinv)
	HpinvVH <- HpVinv %*% H

	logger("Begin computing eigenvalues/vectors of HpinvVH\n")
	eig.HpinvVH <- eigen(HpinvVH, symmetric = TRUE)
	logger("Finished computing eigenvalues/vectors of HpinvVH\n")

	r <- ncol(S)
	n <- length(Z)
	n_mu <- ncol(H)

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
	if (is.null(init)) { init <- list() }
	if (is.null(init$mu_B)) { init$mu_B <- rnorm(n_mu) }
	if (is.null(init$eta)) { init$eta <- rnorm(r) }
	if (is.null(init$xi)) { init$xi <- numeric(n) }
	if (is.null(init$sig2mu)) { init$sig2mu <- 1 }
	if (is.null(init$sig2xi)) { init$sig2xi <- 1 }
	if (is.null(init$sig2K)) { init$sig2K <- 1 }
	mu_B <- init$mu_B
	eta <- init$eta
	xi <- init$xi
	sig2mu <- init$sig2mu
	sig2xi <- init$sig2xi
	sig2K <- init$sig2K
	Y <- S %*% eta + H %*% mu_B + xi

	# Selection of fixed values
	# mu_b, eta, and xi entries should contain indices of the coordinates to keep fixed
	if (is.null(fixed)) { init <- list() }
	if (is.null(fixed$mu_B)) { fixed$mu_B <- integer(0) }
	if (is.null(fixed$eta)) { fixed$eta <- integer(0) }
	if (is.null(fixed$xi)) { fixed$xi <- integer(0) }
	if (is.null(fixed$sig2mu)) { fixed$sig2mu <- FALSE }
	if (is.null(fixed$sig2xi)) { fixed$sig2xi <- FALSE }
	if (is.null(fixed$sig2K)) { fixed$sig2K <- FALSE }

	logger("Begin computing SpinvV\n")
	SpinvV <- matrix(NA, r, n)
	for (j in 1:r) {
		SpinvV[j,] <- S[,j] / sig2eps
	}
	SpinvVS <- SpinvV %*% S
	logger("Finished computing SpinvV\n")

	timer <- list(mu_B = 0, sig2mu = 0, eta = 0, xi = 0, sig2xi = 0, sig2K = 0, Y = 0)

	logger("Begin Gibbs sampler\n")

	for (tt in 1:R) {
		if (tt %% report.period == 0) {
			logger("Begin iteration %d, using %0.2f GB RAM\n", tt, mem_used() / 2^30)
		}

		# Full Conditional for mu_B, using sparse matrix inverse
		st <- Sys.time()
		PostLam <- 1 / (eig.HpinvVH$values + 1/sig2mu)
		V.mu_B <- eig.HpinvVH$vectors %*% (PostLam * t(eig.HpinvVH$vectors))
		mean.mu_B <- V.mu_B %*% (HpVinv %*% (Z - S %*% eta - xi))
		mu_B.new <- mean.mu_B + (eig.HpinvVH$vectors %*% (sqrt(PostLam) * rnorm(n_mu)))
		idx <- setdiff(1:n_mu, fixed$mu_B)
		mu_B[idx] <- mu_B.new[idx]
		timer$mu_B <- timer$mu_B + as.numeric(Sys.time() - st, units = "secs")

		# Full Conditional for sig2mu
		st <- Sys.time()
		scale <- as.numeric(0.5 * t(mu_B) %*% mu_B)
		sig2mu.new <- 1 / rgamma(1, n_mu/2 + 1, 0.000001 + scale)
		if (!fixed$sig2mu) { sig2mu <- sig2mu.new }
		timer$sig2mu <- timer$sig2mu + as.numeric(Sys.time() - st, units = "secs")

		# Full Conditional for eta
		st <- Sys.time()
		zresid <- Z - H %*% mu_B - xi
		V.eta <- solve(as.matrix((1/sig2K * C.inv) + SpinvVS))
		mean.eta <- V.eta %*% (SpinvV %*% zresid)
		eta.new <- mean.eta + chol(V.eta) %*% rnorm(r)
		idx <- setdiff(1:r, fixed$eta)
		eta[idx] <- eta.new[idx]
		timer$eta <- timer$eta + as.numeric(Sys.time() - st, units = "secs")

		# Full Conditional for xi, using sparse matrix inverse
		st <- Sys.time()
		Sigma.xi.inv <- Vinv + 1/sig2xi
		Sigma.xi <- 1 / Sigma.xi.inv
		mean.xi <- Sigma.xi * Vinv * (Z - S %*% eta - H %*% mu_B)
		xi.new <- mean.xi + sqrt(Sigma.xi) * rnorm(n)
		idx <- setdiff(1:n, fixed$xi)
		xi[idx] <- xi.new[idx]
		timer$xi <- timer$xi + as.numeric(Sys.time() - st, units = "secs")

		# Full Conditional sig2xi
		st <- Sys.time()
		scale <- as.numeric(0.5 * t(xi) %*% xi)
		sig2xi.new <- 1 / rgamma(1, n/2 + 1, 0.000001 + scale)
		if (!fixed$sig2xi) { sig2xi <- sig2xi.new }
		timer$sig2xi <- timer$sig2xi + as.numeric(Sys.time() - st, units = "secs")

		# Full Conditional for sig2K
		st <- Sys.time()
		scale <- as.numeric(0.5 * t(eta) %*% C.inv %*% eta)
		sig2K.new <- 1 / rgamma(1, r/2 + 1, 0.000001 + scale)
		if (!fixed$sig2K) { sig2K <- sig2K.new }
		timer$sig2K <- timer$sig2K + as.numeric(Sys.time() - st, units = "secs")

		# Update Y
		st <- Sys.time()
		Y <- S %*% eta + H %*% mu_B + xi
		timer$Y <- timer$Y + as.numeric(Sys.time() - st, units = "secs")

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
		Y.hist = Y.hist,
		elapsed.sec = timer
	)
}

