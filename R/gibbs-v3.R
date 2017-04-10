# In this version of the MCMC, we have multiplicative errors instead of
# additive ones
gibbs.stcos.v3 <- function(Z, S, sig2eps, C.inv, H, R,
	report.period = R+1, burn = 0, thin = 1,
	init = NULL, fixed = NULL, hyper = NULL)
{
	stopifnot(R > burn)

	# Vinv <- 1 / sig2eps
	# HpVinv <- t(H) %*% Diagonal(n = length(Vinv), x = Vinv)
	# HpinvVH <- HpVinv %*% H

	# logger("Begin computing eigenvalues/vectors of HpinvVH\n")
	# eig.HpinvVH <- eigen(HpinvVH, symmetric = TRUE)
	# logger("Finished computing eigenvalues/vectors of HpinvVH\n")

	r <- ncol(S)
	n <- length(Z)
	n_mu <- ncol(H)

	tt.keep <- 0
	R.keep <- ceiling((R - burn) / thin)
	mu_B.hist <- matrix(NA, R.keep, n_mu)
	gamma.hist <- matrix(NA, R.keep, n)
	eta.hist <- matrix(NA, R.keep, r)
	sig2mu.hist <- matrix(NA, R.keep, 1)
	sig2K.hist <- matrix(NA, R.keep, 1)

	# Initial values
	if (is.null(init)) { init <- list() }
	if (is.null(init$mu_B)) { init$mu_B <- rnorm(n_mu) }
	if (is.null(init$eta)) { init$eta <- rnorm(r) }
	if (is.null(init$gamma)) { init$gamma <- numeric(n) }
	if (is.null(init$sig2mu)) { init$sig2mu <- 1 }
	if (is.null(init$sig2K)) { init$sig2K <- 1 }
	mu_B <- init$mu_B
	eta <- init$eta
	gamma <- init$gamma
	sig2mu <- init$sig2mu
	sig2K <- init$sig2K

	# Selection of fixed values
	# mu_b, eta, and xi entries should contain indices of the coordinates to keep fixed
	if (is.null(fixed)) { init <- list() }
	if (is.null(fixed$mu_B)) { fixed$mu_B <- integer(0) }
	if (is.null(fixed$eta)) { fixed$eta <- integer(0) }
	if (is.null(fixed$gamma)) { fixed$gamma <- integer(0) }
	if (is.null(fixed$sig2mu)) { fixed$sig2mu <- FALSE }
	if (is.null(fixed$sig2K)) { fixed$sig2K <- FALSE }

	# Hyperparameters
	if (is.null(hyper)) { hyper <- list() }
	if (is.null(hyper$a.sig2mu)) { hyper$a.sig2mu <- 2 }
	if (is.null(hyper$a.sig2K)) { hyper$a.sig2K <- 2 }
	if (is.null(hyper$a.gamma)) { hyper$a.gamma <- 2 }
	if (is.null(hyper$b.sig2mu)) { hyper$b.sig2mu <- 2 }
	if (is.null(hyper$b.sig2K)) { hyper$b.sig2K <- 2 }
	if (is.null(hyper$b.gamma)) { hyper$b.gamma <- 2 }

	# logger("Begin computing SpinvV\n")
	# SpinvV <- matrix(NA, r, n)
	# for (j in 1:r) {
	# 	SpinvV[j,] <- S[,j] / sig2eps
	# }
	# SpinvVS <- SpinvV %*% S
	# logger("Finished computing SpinvV\n")

	timer <- list(mu_B = 0, eta = 0, gamma = 0, sig2mu = 0, sig2K = 0)

	logger("Begin Gibbs sampler\n")

	for (tt in 1:R) {
		if (tt %% report.period == 0) {
			logger("Begin iteration %d, using %0.2f GB RAM\n", tt, mem_used() / 2^30)
		}

		# [DONE]
		# Full Conditional for sig2K
		st <- Sys.time()
		scale <- as.numeric(0.5 * t(eta) %*% C.inv %*% eta)
		sig2K.new <- 1 / rgamma(1, r/2 + hyper$a.sig2K, hyper$b.sig2K + scale)
		if (!fixed$sig2K) { sig2K <- sig2K.new }
		timer$sig2K <- timer$sig2K + as.numeric(Sys.time() - st, units = "secs")

		# [DONE]
		# Full Conditional for sig2mu
		st <- Sys.time()
		scale <- as.numeric(0.5 * t(mu_B) %*% mu_B)
		sig2mu.new <- 1 / rgamma(1, n_mu/2 + hyper$a.sig2mu, hyper$b.sig2mu + scale)
		if (!fixed$sig2mu) { sig2mu <- sig2mu.new }
		timer$sig2mu <- timer$sig2mu + as.numeric(Sys.time() - st, units = "secs")

		# [DONE]
		# Full Conditional for gamma
		st <- Sys.time()
		ee <- Z - S %*% eta - H %*% mu_B
		gamma.new <- 1 / rgamma(N, hyper$a.gamma + 1/2, hyper$b.gamma + as.numeric(ee^2 / (2*sig2eps)))
		idx <- setdiff(1:n, fixed$gamma)
		gamma[idx] <- gamma.new[idx]
		timer$gamma <- timer$gamma + as.numeric(Sys.time() - st, units = "secs")

		V.inv <- 1 / (sig2eps*gamma)
		
		# [DONE]
		# Full Conditional for eta
		# st <- Sys.time()
		# Z.resid <- Z - H %*% mu_B
		# V.eta <- solve(as.matrix((1/sig2K * C.inv) + SpinvVS))
		# mean.eta <- V.eta %*% (SpinvV %*% Z.resid)
		# eta.new <- mean.eta + chol(V.eta) %*% rnorm(r)
		# idx <- setdiff(1:r, fixed$eta)
		# eta[idx] <- eta.new[idx]
		# timer$eta <- timer$eta + as.numeric(Sys.time() - st, units = "secs")
browser()
		st <- Sys.time()
		Z.resid <- Z - H %*% mu_B
		V.eta.inv <- t(S) %*% (V.inv * S) + (1/sig2K * C.inv)
		V.eta <- solve(V.eta.inv)
		V.eta <- V.eta/2 + t(V.eta)/2
		mean.eta <- V.eta %*% (t(S) %*% (V.inv * Z.resid))
		eta.new <- mean.eta + chol(V.eta) %*% rnorm(r)
		idx <- setdiff(1:r, fixed$eta)
		eta[idx] <- eta.new[idx]
		timer$eta <- timer$eta + as.numeric(Sys.time() - st, units = "secs")

		# [DONE]
		# Full Conditional for mu_B, using sparse matrix inverse
		# st <- Sys.time()
		# PostLam <- 1 / (eig.HpinvVH$values + 1/sig2mu)
		# V.mu_B <- eig.HpinvVH$vectors %*% (PostLam * t(eig.HpinvVH$vectors))
		# mean.mu_B <- V.mu_B %*% (HpVinv %*% (Z - S %*% eta - xi))
		# mu_B.new <- mean.mu_B + (eig.HpinvVH$vectors %*% (sqrt(PostLam) * rnorm(n_mu)))
		# idx <- setdiff(1:n_mu, fixed$mu_B)
		# mu_B[idx] <- mu_B.new[idx]
		# timer$mu_B <- timer$mu_B + as.numeric(Sys.time() - st, units = "secs")

		# Note that we can do some of this more efficiently because of the matrix forms
		st <- Sys.time()
		Z.resid <- Z - S %*% eta
		V.mu.inv <- t(H) %*% (V.inv * H) + Diagonal(n_mu, sig2mu)
		V.mu <- solve(V.mu.inv)
		mean.mu <- V.mu %*% (t(H) %*% (V.inv * Z.resid))
		system.time(chol(V.mu))
		system.time(mu_B.new <- mean.mu + chol(V.mu) %*% rnorm(n_mu))
		idx <- setdiff(1:n_mu, fixed$mu_B)
		mu_B[idx] <- mu_B.new[idx]
		timer$mu_B <- timer$mu_B + as.numeric(Sys.time() - st, units = "secs")
		
		st <- Sys.time()
		Z.resid <- Z - S %*% eta
		system.time(eig <- eigen(t(H) %*% (V.inv * H)))

		# Save history
		# TBD: This is very expensive to keep in memory. We should be able to
		# transparently store it on disk, if the user wishes
		if ((tt > burn) & (tt %% thin == 0)) {
			tt.keep <- tt.keep + 1
			mu_B.hist[tt.keep,] <- as.numeric(mu_B)
			eta.hist[tt.keep,] <- as.numeric(eta)
			gamma.hist[tt.keep,] <- as.numeric(gamma)
			sig2mu.hist[tt.keep] <- sig2mu
			sig2K.hist[tt.keep] <- sig2K
		}
	}

	logger("Finished Gibbs sampler\n")

	list(mu_B.hist = mu_B.hist,
		gamma.hist = gamma.hist,
		eta.hist = eta.hist,
		sig2mu.hist = sig2mu.hist,
		sig2K.hist = sig2K.hist,
		elapsed.sec = timer
	)
}
