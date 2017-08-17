gibbs.stcos <- function(prep, R, report.period = R+1, burn = 0, thin = 1,
	hyper = NULL, sig2xi.init = NULL)
{
	Z <- prep$get_Z()
	V <- prep$get_V()
	H <- prep$get_H()
	S <- prep$get_reduced_S()
	K.inv <- prep$get_Kinv()
	mle.out <- mle.stcos(Z, S, V, H, init = list(sig2xi.init))

	init <- list(
		sig2xi = mle.out$sig2xi.hat,
		mu_B = mle.out$mu.hat,
		eta = mle.out$eta.hat
	)
	gibbs.out <- gibbs.stcos.raw(Z = Z, S = S, V = V, K.inv = K.inv, H = H,
		R = R, report.period = report.period, burn = burn, thin = thin,
		init = init, hyper = hyper)
}

gibbs.stcos.raw <- function(Z, S, V, K.inv, H, R, report.period = R+1,
	burn = 0, thin = 1, init = NULL, fixed = NULL, hyper = NULL)
{
	timer <- list(mu_B = 0, sig2mu = 0, eta = 0, xi = 0, sig2xi = 0, sig2K = 0,
		Y = 0, pre = 0, post = 0)

	st <- Sys.time()
	stopifnot(R > burn)

	r <- ncol(S)
	n <- length(Z)
	n_mu <- ncol(H)

	Vinv <- 1 / V
	HpVinv <- t(H) %*% Diagonal(n = length(Vinv), x = Vinv)
	HpinvVH <- HpVinv %*% H

	logger("Begin computing eigenvalues/vectors of HpinvVH\n")
	eig.HpinvVH <- eigen(HpinvVH, symmetric = TRUE)
	eig.HpinvVH$vectors <- Matrix(eig.HpinvVH$vectors)
	logger("Finished computing eigenvalues/vectors of HpinvVH\n")

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
	if (is.null(init$sig2mu)) { init$sig2mu <- 1 }
	if (is.null(init$sig2xi)) { init$sig2xi <- 1 }
	if (is.null(init$sig2K)) { init$sig2K <- 1 }
	if (is.null(init$mu_B)) { init$mu_B <- rnorm(n_mu, 0, sqrt(init$sig2mu)) }
	if (is.null(init$eta)) {
		eig <- eigen(K.inv)
		ee <- eig$values
		ee[ee <= 0] <- min(ee[ee > 0])
		ee[ee > 0] <- 1 / ee[ee > 0]
		C <- eig$vectors %*% (ee * t(eig$vectors))
		K.half <- chol(init$sig2K * C)
		init$eta <- K.half %*% rnorm(r)
	}
	if (is.null(init$xi)) { init$xi <- rnorm(n, 0, sqrt(init$sig2xi)) }
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

	# Hyperparameters
	if (is.null(hyper)) { hyper <- list() }
	if (is.null(hyper$a.sig2mu)) { hyper$a.sig2mu <- 2 }
	if (is.null(hyper$a.sig2K)) { hyper$a.sig2K <- 2 }
	if (is.null(hyper$a.sig2xi)) { hyper$a.sig2xi <- 2 }
	if (is.null(hyper$b.sig2mu)) { hyper$b.sig2mu <- 2 }
	if (is.null(hyper$b.sig2K)) { hyper$b.sig2K <- 2 }
	if (is.null(hyper$b.sig2xi)) { hyper$b.sig2xi <- 2 }

	logger("Begin computing SpinvV\n")
	SpinvV <- matrix(NA, r, n)
	for (j in 1:r) {
		SpinvV[j,] <- S[,j] / V
	}
	SpinvVS <- SpinvV %*% S
	logger("Finished computing SpinvV\n")
	timer$pre <- timer$pre + as.numeric(Sys.time() - st, units = "secs")

	logger("Begin Gibbs sampler\n")

	for (tt in 1:R) {
		if (tt %% report.period == 0) {
			logger("Begin iteration %d\n", tt)
		}

		# Draw from [xi | ---] using sparse matrix inverse
		st <- Sys.time()
		Sigma.xi.inv <- Vinv + 1/sig2xi
		Sigma.xi <- 1 / Sigma.xi.inv
		mean.xi <- Sigma.xi * Vinv * (Z - S %*% eta - H %*% mu_B)
		xi.new <- mean.xi + sqrt(Sigma.xi) * rnorm(n)
		idx <- setdiff(1:n, fixed$xi)
		xi[idx] <- xi.new[idx]
		timer$xi <- timer$xi + as.numeric(Sys.time() - st, units = "secs")

		# Draw from [sig2xi | ---]
		st <- Sys.time()
		scale <- as.numeric(0.5 * t(xi) %*% xi)
		sig2xi.new <- 1 / rgamma(1, n/2 + hyper$a.sig2xi, hyper$b.sig2xi + scale)
		if (!fixed$sig2xi) { sig2xi <- sig2xi.new }
		timer$sig2xi <- timer$sig2xi + as.numeric(Sys.time() - st, units = "secs")

		# Draw from [sig2mu | ---]
		st <- Sys.time()
		scale <- as.numeric(0.5 * t(mu_B) %*% mu_B)
		sig2mu.new <- 1 / rgamma(1, n_mu/2 + hyper$a.sig2mu, hyper$b.sig2mu + scale)
		if (!fixed$sig2mu) { sig2mu <- sig2mu.new }
		timer$sig2mu <- timer$sig2mu + as.numeric(Sys.time() - st, units = "secs")

		# Draw from [sig2K | ---]
		st <- Sys.time()
		scale <- as.numeric(0.5 * t(eta) %*% K.inv %*% eta)
		sig2K.new <- 1 / rgamma(1, r/2 + hyper$a.sig2K, hyper$b.sig2K + scale)
		if (!fixed$sig2K) { sig2K <- sig2K.new }
		timer$sig2K <- timer$sig2K + as.numeric(Sys.time() - st, units = "secs")

		# Draw from [eta | ---]
		st <- Sys.time()
		V.eta <- solve(as.matrix((1/sig2K * K.inv) + SpinvVS))
		mean.eta <- V.eta %*% (SpinvV %*% (Z - H %*% mu_B - xi))
		eta.new <- mean.eta + chol(V.eta) %*% rnorm(r)
		idx <- setdiff(1:r, fixed$eta)
		eta[idx] <- eta.new[idx]
		timer$eta <- timer$eta + as.numeric(Sys.time() - st, units = "secs")

		# Draw from [mu_B | ---] using sparse matrix inverse
		# We use the ugly version below because it avoids a large matrix-matrix multiplcation
		st <- Sys.time()
		PostLam <- 1 / (eig.HpinvVH$values + 1/sig2mu)
		# V.mu_B <- eig.HpinvVH$vectors %*% (PostLam * t(eig.HpinvVH$vectors))		# Nice-looking version
		# mean.mu_B <- V.mu_B %*% (HpVinv %*% (Z - S %*% eta - xi))					# Nice-looking version
		mean.mu_B <- eig.HpinvVH$vectors %*% ((PostLam * t(eig.HpinvVH$vectors)) %*% (HpVinv %*% (Z - S %*% eta - xi)))		# Ugly version
		mu_B.new <- mean.mu_B + (eig.HpinvVH$vectors %*% (sqrt(PostLam) * rnorm(n_mu)))
		idx <- setdiff(1:n_mu, fixed$mu_B)
		mu_B[idx] <- mu_B.new[idx]
		timer$mu_B <- timer$mu_B + as.numeric(Sys.time() - st, units = "secs")

		# Update Y (which is also standardized)
		st <- Sys.time()
		Y <- H %*% mu_B + S %*% eta + xi
		timer$Y <- timer$Y + as.numeric(Sys.time() - st, units = "secs")

		# Save history
		# TBD: This can be very expensive to keep in memory. We should be able to
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

	ret <- list(mu_B.hist = mu_B.hist,
		xi.hist = xi.hist,
		eta.hist = eta.hist,
		sig2mu.hist = sig2mu.hist,
		sig2xi.hist = sig2xi.hist,
		sig2K.hist = sig2K.hist,
		Y.hist = Y.hist,
		Z = Z, H = H, S = S, V = V, R.keep = R.keep,
		elapsed.sec = timer
	)
	class(ret) <- "stcos"

	st <- Sys.time()
	ret$loglik <- logLik(ret)
	ret$dic <- DIC(ret)
	timer$post <- timer$post + as.numeric(Sys.time() - st, units = "secs")

	return(ret)
}

logLik.stcos <- function(object, ...)
{
	if (!is.null(object$loglik)) {
		return(object$loglik)
	}

	R.keep <- object$R.keep
	loglik.mcmc <- numeric(R.keep)
	for (r in 1:R.keep) {
		mu_B <- object$mu_B.hist[r,]
		eta <- object$eta.hist[r,]
		sig2xi <- object$sig2xi.hist[r]
		loglik.mcmc[r] <- sum(
			dnorm(object$Z, as.numeric(object$H %*% mu_B + object$S %*% eta),
			sqrt(object$V + sig2xi),
			log = TRUE))
	}
	return(loglik.mcmc)
}

DIC.stcos <- function(object, ...)
{
	if (!is.null(object$dic)) {
		return(object$dic)
	}

	loglik.mcmc <- logLik.stcos(object)
	mu_B.bar <- colMeans(object$mu_B.hist)
	eta.bar <- colMeans(object$eta.hist)
	sig2xi.bar <- mean(object$sig2xi.hist)
	D.thetabar <- -2 * sum(
		dnorm(object$Z, as.numeric(object$H %*% mu_B.bar + object$S %*% eta.bar),
		sqrt(object$V + sig2xi.bar),
		log = TRUE))
	D.bar <- mean(-2*loglik.mcmc)
	D.thetabar + 2*(D.bar - D.thetabar)
}

# TBD: We can compute summaries ourselves and remove dependency on coda
print.stcos <- function (x, ...)
{
	variances.mcmc <- mcmc(cbind(x$sig2mu.hist, x$sig2K.hist, x$sig2xi.hist))
	colnames(variances.mcmc) <- c("sig2mu", "sig2K", "sig2xi")

	total.sec <- sum(unlist(x$elapsed.sec))
	hh <- floor(total.sec / 60^2)
	mm <- floor((total.sec - hh*60^2) / 60)
	ss <- round(total.sec - hh*60^2 - mm*60)

	printf("Fit for STCOS model\n")
	printf("--\n")
	print(summary(variances.mcmc)$statistics)
	printf("--\n")
	printf("Saved %d draws\n", x$R.keep)
	printf("DIC: %f\n", x$dic)
	printf("Elapsed time: %02d:%02d:%02d\n", hh, mm, ss)
}

fitted.stcos <- function (object, H, S, ...)
{
	R.keep <- object$R.keep
	n <- nrow(H)
	E.mcmc <- matrix(NA, R.keep, n)
	for (r in 1:R.keep) {
		mu_B <- object$mu_B.hist[r,]
		eta <- object$eta.hist[r,]
		E.mcmc[r,] <- as.numeric(H %*% mu_B + S %*% eta)
	}
	return(E.mcmc)
}

predict.stcos <- function (object, H, S, ...)
{
	R.keep <- object$R.keep
	n <- nrow(H)
	Y.mcmc <- matrix(NA, R.keep, n)
	for (r in 1:R.keep) {
		mu_B <- object$mu_B.hist[r,]
		eta <- object$eta.hist[r,]
		xi <- object$xi.hist[r,]
		sig2xi <- object$sig2xi.hist[r]
		Y.mcmc[r,] <- rnorm(n, as.numeric(H %*% mu_B + S %*% eta), sqrt(sig2xi))
	}
	return(Y.mcmc)
}
