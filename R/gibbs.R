#' Gibbs Sampler for STCOS Model
#'
#' Run the Gibbs sampling algorithm for the STCOS model. \code{gibbs.stcos}
#' presents a simplified interface, while \code{gibbs.stcos.raw} allows
#' all inputs to be specified separately.
#'
#' @param z Vector which represents the outcome; assumed to be direct
#'        estimates from the survey.
#' @param S Design matrix for basis decomposition.
#' @param v Vector which represents direct variance estimates from the survey.
#' @param K.inv Inverse of the \eqn{K} matrix, which is the covariance of the
#'        random coefficient \eqn{\eta}
#' @param H Matrix of overlaps between source and fine-level supports.
#' @param init A list containing the following initial values for the MCMC:
#' 	      \code{sig2mu}, \code{sig2xi}, \code{sig2K}, \code{mu_B}, \code{eta},
#' 	      \code{xi}. Any values which are not specified are set to arbitrary
#' 	      choices.
#' @param fixed A list specifying which parameters to keep fixed in the MCMC.
#'        This can normally be left blank. If elements \code{sig2mu},
#'        \code{sig2xi}, or \code{sig2K} are specified they should be boolean,
#'        where TRUE means fixed (i.e. not drawn). If elements \code{mu_B},
#'        \code{eta}, or \code{xi} are specified, they should each be a vector
#'        of indicies; the specified indices are to be treated as fixed (i.e.
#'        not drawn).
#' @param prep An \code{STCOSprep} object.
#' @param R Number of MCMC reps.
#' @param report.period Gibbs sampler will report progress each time this many
#'        iterations are completed.
#' @param burn Burn this many of \code{R} the draws, before saving history.
#' @param thin After burn-ikn period, save one out of every \code{thin} draws.
#' @param hyper A list containing the following hyperparameter values:
#' 	      \code{a.sig2mu}, \code{a.sig2K}, \code{a.sig2xi}, \code{b.sig2mu},
#' 	      \code{b.sig2K}, \code{b.sig2xi}. Any hyperparameters which are not
#' 	      specified are set to a default value of 2.
#'
#' @return An \code{stcos} object which contains draws from the sampler.
#'
#' @examples
#' \dontrun{
#' basis <- SpaceTimeBisquareBasis$new(x, y, t, w.s, w.t)
#' sp <- STCOSPrep$new(fine_domain = dom.fine,
#'     fine_domain_geo_name = "GEO_ID",
#'     basis = basis, basis_mc_reps = 500)
#' out1 <- gibbs.stcos(sp, R = 10000, burn = 0, thin = 1)
#' 
#' out2 <- gibbs.stcos.raw(z = sp$get_z(), S = sp$get_reduced_S(),
#'     v = sp$get_v(), K.inv = sp$get_Kinv(), H = sp$get_H(),
#'     R = 10000, burn = 0, thin = 1)
#' }
#' @name gibbs
NULL

#' @name gibbs
#' @export
gibbs.stcos <- function(prep, R, report.period = R+1, burn = 0, thin = 1,
	hyper = NULL, init = NULL)
{
	z <- prep$get_z()
	v <- prep$get_v()
	H <- prep$get_H()
	S <- prep$get_reduced_S()
	K.inv <- prep$get_Kinv()

	gibbs.out <- gibbs.stcos.raw(z = z, v = v, H = H, S = S, K.inv = K.inv,
		R = R, report.period = report.period, burn = burn, thin = thin,
		init = init, hyper = hyper)
}

#' @name gibbs
#' @export
gibbs.stcos.raw <- function(z, v, H, S, K.inv, R, report.period = R+1,
	burn = 0, thin = 1, init = NULL, fixed = NULL, hyper = NULL)
{
	timer <- list(mu_B = 0, sig2mu = 0, eta = 0, xi = 0, sig2xi = 0, sig2K = 0,
		Y = 0, pre = 0, post = 0)

	st <- Sys.time()
	stopifnot(R > burn)

	r <- ncol(S)
	n <- length(z)
	n_mu <- ncol(H)

	vinv <- 1 / v
	HpVinv <- t(H) %*% Diagonal(n = length(vinv), x = vinv)
	HpinvVH <- HpVinv %*% H

	eig.HpinvVH <- eigen(HpinvVH, symmetric = TRUE)
	eig.HpinvVH$vectors <- Matrix(eig.HpinvVH$vectors)

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

	SpinvV <- matrix(NA, r, n)
	for (j in 1:r) {
		SpinvV[j,] <- S[,j] / v
	}
	SpinvVS <- SpinvV %*% S
	timer$pre <- timer$pre + as.numeric(Sys.time() - st, units = "secs")

	logger("Begin Gibbs sampler\n")

	for (tt in 1:R) {
		if (tt %% report.period == 0) {
			logger("Begin iteration %d\n", tt)
		}

		# Draw from [xi | ---] using sparse matrix inverse
		st <- Sys.time()
		Sigma.xi.inv <- vinv + 1/sig2xi
		Sigma.xi <- 1 / Sigma.xi.inv
		mean.xi <- Sigma.xi * vinv * (z - S %*% eta - H %*% mu_B)
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
		mean.eta <- V.eta %*% (SpinvV %*% (z - H %*% mu_B - xi))
		eta.new <- mean.eta + chol(V.eta) %*% rnorm(r)
		idx <- setdiff(1:r, fixed$eta)
		eta[idx] <- eta.new[idx]
		timer$eta <- timer$eta + as.numeric(Sys.time() - st, units = "secs")

		# Draw from [mu_B | ---] using sparse matrix inverse
		# We use the computation below because it avoids a large matrix-matrix multiplcation
		st <- Sys.time()
		PostLam <- 1 / (eig.HpinvVH$values + 1/sig2mu)
		if (FALSE) {
			# Simpler version
			V.mu_B <- eig.HpinvVH$vectors %*% (PostLam * t(eig.HpinvVH$vectors))		
			mean.mu_B <- V.mu_B %*% (HpVinv %*% (z - S %*% eta - xi))
			mu_B.new <- chol(V.mu_B) %*% rnorm(n_mu) + mean.mu_B
		} else {
			# Version that avoids large multiplication
			mean.mu_B <- eig.HpinvVH$vectors %*% ((PostLam * t(eig.HpinvVH$vectors)) %*%
				(HpVinv %*% (z - S %*% eta - xi)))
			mu_B.new <- mean.mu_B + (eig.HpinvVH$vectors %*% (sqrt(PostLam) * rnorm(n_mu)))
		}
		idx <- setdiff(1:n_mu, fixed$mu_B)
		mu_B[idx] <- mu_B.new[idx]
		timer$mu_B <- timer$mu_B + as.numeric(Sys.time() - st, units = "secs")

		# Update Y
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
		z = z, H = H, S = S, v = v, R.keep = R.keep,
		elapsed.sec = timer
	)
	class(ret) <- "stcos"

	st <- Sys.time()
	ret$loglik <- logLik(ret)
	ret$dic <- DIC(ret)
	timer$post <- timer$post + as.numeric(Sys.time() - st, units = "secs")

	return(ret)
}

#' @export
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
			dnorm(object$z, as.numeric(object$H %*% mu_B + object$S %*% eta),
			sqrt(object$v + sig2xi),
			log = TRUE))
	}
	return(loglik.mcmc)
}

#' Deviance Information Criterion
#' 
#' A function to compute the Deviance Information Criterion (DIC) on an
#' \code{stcos} object.
#'
#' @param object A result from the Gibbs sampler.
#'
#' @return DIC computed from saved draws
#' @export
#'
#' @examples
#' \dontrun{
#' out1 <- gibbs.stcos(sp, R = 10000, burn = 0, thin = 1)
#' DIC(out1)
#' }
#' @seealso \code{\link{gibbs.stcos}} \code{\link{gibbs.stcos.raw}}
DIC <- function(object)
{
	if (!is.null(object$dic)) {
		return(object$dic)
	}

	loglik.mcmc <- logLik.stcos(object)
	mu_B.bar <- colMeans(object$mu_B.hist)
	eta.bar <- colMeans(object$eta.hist)
	sig2xi.bar <- mean(object$sig2xi.hist)
	D.thetabar <- -2 * sum(
		dnorm(object$z, as.numeric(object$H %*% mu_B.bar + object$S %*% eta.bar),
		sqrt(object$v + sig2xi.bar),
		log = TRUE))
	D.bar <- mean(-2*loglik.mcmc)
	D.thetabar + 2*(D.bar - D.thetabar)
}

#' @method print stcos
#' @export
print.stcos <- function (x, ...)
{
	# We could compute summaries ourselves and remove dependency on coda ...
	variances.mcmc <- cbind(x$sig2mu.hist, x$sig2K.hist, x$sig2xi.hist)
	colnames(variances.mcmc) <- c("sig2mu", "sig2K", "sig2xi")
	summary.mcmc <- data.frame(
		apply(variances.mcmc, 2, mean),
		apply(variances.mcmc, 2, sd),
		apply(variances.mcmc, 2, quantile, probs = 0.025),
		apply(variances.mcmc, 2, quantile, probs = 0.25),
		apply(variances.mcmc, 2, quantile, probs = 0.75),
		apply(variances.mcmc, 2, quantile, probs = 0.975)
	)
	colnames(summary.mcmc) <- c("Mean", "SD", "2.5%", "25%", "75%", "97.5%")

	total.sec <- sum(unlist(x$elapsed.sec))
	hh <- floor(total.sec / 60^2)
	mm <- floor((total.sec - hh*60^2) / 60)
	ss <- round(total.sec - hh*60^2 - mm*60)

	printf("Fit for STCOS model\n")
	printf("--\n")
	print(summary.mcmc)
	printf("--\n")
	printf("Saved %d draws\n", x$R.keep)
	printf("DIC: %f\n", x$dic)
	printf("Elapsed time: %02d:%02d:%02d\n", hh, mm, ss)

	invisible(x)
}

#' @export
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

#' @export
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
