#' Gibbs Sampler for STCOS Model
#'
#' @param z Vector which represents the outcome; assumed to be direct
#'        estimates from the survey.
#' @param v Vector which represents direct variance estimates from the survey.
#' @param H Matrix of overlaps between source and fine-level supports.
#' @param S Design matrix for basis decomposition.
#' @param Kinv The precision matrix \eqn{\bm{K}^{-1}} of the
#'        random coefficient \eqn{\bm{\eta}}
#' @param init A list containing the following initial values for the MCMC:
#' 	      \code{sig2mu}, \code{sig2xi}, \code{sig2K}, \code{muB}, \code{eta},
#' 	      \code{xi}. Any values which are not specified are set to arbitrary
#' 	      choices.
#' @param fixed A list specifying which parameters to keep fixed in the MCMC.
#'        This can normally be left blank. If elements \code{sig2mu},
#'        \code{sig2xi}, or \code{sig2K} are specified they should be boolean,
#'        where TRUE means fixed (i.e. not drawn). If elements \code{muB},
#'        \code{eta}, or \code{xi} are specified, they should each be a vector
#'        of indicies; the specified indices are to be treated as fixed (i.e.
#'        not drawn).
#' @param R Number of MCMC reps.
#' @param report_period Gibbs sampler will report progress each time this many
#'        iterations are completed.
#' @param burn Number of the \code{R} draws to discard at the beginning of the
#'        chain.
#' @param thin After burn-in period, save one out of every \code{thin} draws.
#' @param hyper A list containing the following hyperparameter values:
#' 	      \code{a_sig2mu}, \code{a_sig2K}, \code{a_sig2xi}, \code{b_sig2mu},
#' 	      \code{b_sig2K}, \code{b_sig2xi}. Any hyperparameters which are not
#' 	      specified are set to a default value of 2.
#' @param object A result from \code{gibbs_stcos}.
#' @param x A result from \code{gibbs_stcos}.
#' @param ... Additional arguments.
#'
#' @return \code{gibbs_stcos} returns an \code{stcos} object which contains
#' draws from the sampler. Helper functions take this object as an input
#' and produce various outputs (see details).
#'
#' @details 
#' Fits the model
#' \deqn{
#'   \bm{Z} = \bm{H} \bm{\mu}_B + \bm{S} \bm{\eta} + \bm{\xi} + \bm{\varepsilon}, \quad
#'   \bm{\varepsilon} \sim \textrm{N}(0, \bm{V}),
#' }
#' \deqn{
#'   \bm{\eta} \sim \textrm{N}(\bm{0}, \sigma_K^2 \bm{K}), \quad
#'   \bm{\xi} \sim \textrm{N}(0, \sigma_{\xi}^2 \bm{I}),
#' }
#' \deqn{
#' \bm{\mu}_B \sim \textrm{N}(\bm{0}, \sigma_\mu^2 \bm{I}), \quad
#' \sigma_\mu^2 \sim \textrm{IG}(a_\mu, b_\mu),
#' }
#' \deqn{
#' \sigma_K^2 \sim \textrm{IG}(a_K, b_K), \quad
#' \sigma_\xi^2 \sim \textrm{IG}(a_\xi, b_\xi),
#' }
#' using a Gibbs sampler with closed-form draws.
#' 
#' Helper functions produce the following outputs:
#' \itemize{
#' \item \code{logLik} computes the log-likelihood for each saved draw.
#' \item \code{DIC} computes the Deviance information criterion for each saved draw.
#' \item \code{print} displays a summary of the draws.
#' \item \code{fitted} computes the mean \eqn{E(Y_i)} for each observation
#'       \eqn{i = 1, \ldots, n}, for each saved draw.
#' \item \code{predict} draws \eqn{Y_i} for each observation
#'       \eqn{i = 1, \ldots, n}, using the parameter values for each saved
#'       Gibbs sampler draw.
#' }
#'
#' @examples
#' \dontrun{
#' demo = prepare_stcos_demo()
#' out = gibbs_stcos(demo$z, demo$v, demo$H, demo$S, solve(demo$K),
#'     R = 100, burn = 0, thin = 1)
#' print(out)
#' logLik(out)
#' DIC(out)
#' fitted(out, demo$H, demo$S)
#' predict(out, demo$H, demo$S)
#' }
#' @name gibbs_stcos
NULL

#' @name gibbs_stcos
#' @export
gibbs_stcos = function(z, v, H, S, Kinv, R, report_period = R+1,
	burn = 0, thin = 1, init = NULL, fixed = NULL, hyper = NULL)
{
	timer = list(muB = 0, sig2mu = 0, eta = 0, xi = 0, sig2xi = 0, sig2K = 0,
		Y = 0, pre = 0, post = 0)

	st = Sys.time()
	stopifnot(R > burn)

	r = ncol(S)
	n = length(z)
	n_mu = ncol(H)

	vinv = 1 / v
	HpVinv = t(H) %*% Diagonal(n = length(vinv), x = vinv)
	HpinvVH = HpVinv %*% H

	eig_HpinvVH = eigen(HpinvVH, symmetric = TRUE)
	eig_HpinvVH$vectors = Matrix(eig_HpinvVH$vectors)

	tt_keep = 0
	R_keep = ceiling((R - burn) / thin)
	muB_hist = matrix(NA, R_keep, n_mu)
	xi_hist = matrix(NA, R_keep, n)
	eta_hist = matrix(NA, R_keep, r)
	sig2mu_hist = matrix(NA, R_keep, 1)
	sig2xi_hist = matrix(NA, R_keep, 1)
	sig2K_hist = matrix(NA, R_keep, 1)
	Y_hist = matrix(NA, R_keep, n)

	# Initial values
	if (is.null(init)) { init = list() }
	if (is.null(init$sig2mu)) { init$sig2mu = 1 }
	if (is.null(init$sig2xi)) { init$sig2xi = 1 }
	if (is.null(init$sig2K)) { init$sig2K = 1 }
	if (is.null(init$muB)) { init$muB = rnorm(n_mu, 0, sqrt(init$sig2mu)) }
	if (is.null(init$eta)) {
		eig = eigen(Kinv)
		ee = eig$values
		ee[ee <= 0] = min(ee[ee > 0])
		ee[ee > 0] = 1 / ee[ee > 0]
		C = eig$vectors %*% (ee * t(eig$vectors))
		K_half = chol(init$sig2K * C)
		init$eta = K_half %*% rnorm(r)
	}
	if (is.null(init$xi)) { init$xi = rnorm(n, 0, sqrt(init$sig2xi)) }
	muB = init$muB
	eta = init$eta
	xi = init$xi
	sig2mu = init$sig2mu
	sig2xi = init$sig2xi
	sig2K = init$sig2K
	Y = S %*% eta + H %*% muB + xi

	# Selection of fixed values
	# muB, eta, and xi entries should contain indices of the coordinates to keep fixed
	if (is.null(fixed)) { init = list() }
	if (is.null(fixed$muB)) { fixed$muB = integer(0) }
	if (is.null(fixed$eta)) { fixed$eta = integer(0) }
	if (is.null(fixed$xi)) { fixed$xi = integer(0) }
	if (is.null(fixed$sig2mu)) { fixed$sig2mu = FALSE }
	if (is.null(fixed$sig2xi)) { fixed$sig2xi = FALSE }
	if (is.null(fixed$sig2K)) { fixed$sig2K = FALSE }

	# Hyperparameters
	if (is.null(hyper)) { hyper = list() }
	if (is.null(hyper$a_sig2mu)) { hyper$a_sig2mu = 2 }
	if (is.null(hyper$a_sig2K)) { hyper$a_sig2K = 2 }
	if (is.null(hyper$a_sig2xi)) { hyper$a_sig2xi = 2 }
	if (is.null(hyper$b_sig2mu)) { hyper$b_sig2mu = 2 }
	if (is.null(hyper$b_sig2K)) { hyper$b_sig2K = 2 }
	if (is.null(hyper$b_sig2xi)) { hyper$b_sig2xi = 2 }

	SpinvV = matrix(NA, r, n)
	for (j in 1:r) {
		SpinvV[j,] = S[,j] / v
	}
	SpinvVS = SpinvV %*% S
	timer$pre = timer$pre + as.numeric(Sys.time() - st, units = "secs")

	logger("Begin Gibbs sampler\n")

	for (tt in 1:R) {
		if (tt %% report_period == 0) {
			logger("Begin iteration %d\n", tt)
		}

		# Draw from [xi | ---] using sparse matrix inverse
		st = Sys.time()
		Sigma_xi_inv = vinv + 1/sig2xi
		Sigma_xi = 1 / Sigma_xi_inv
		mean_xi = Sigma_xi * vinv * (z - S %*% eta - H %*% muB)
		xi.new = mean_xi + sqrt(Sigma_xi) * rnorm(n)
		idx = setdiff(1:n, fixed$xi)
		xi[idx] = xi.new[idx]
		timer$xi = timer$xi + as.numeric(Sys.time() - st, units = "secs")

		# Draw from [sig2xi | ---]
		st = Sys.time()
		scale = as.numeric(0.5 * t(xi) %*% xi)
		sig2xi_new = 1 / rgamma(1, n/2 + hyper$a_sig2xi, hyper$b_sig2xi + scale)
		if (!fixed$sig2xi) { sig2xi = sig2xi_new }
		timer$sig2xi = timer$sig2xi + as.numeric(Sys.time() - st, units = "secs")

		# Draw from [sig2mu | ---]
		st = Sys.time()
		scale = as.numeric(0.5 * t(muB) %*% muB)
		sig2mu_new = 1 / rgamma(1, n_mu/2 + hyper$a_sig2mu, hyper$b_sig2mu + scale)
		if (!fixed$sig2mu) { sig2mu = sig2mu_new }
		timer$sig2mu = timer$sig2mu + as.numeric(Sys.time() - st, units = "secs")

		# Draw from [sig2K | ---]
		st = Sys.time()
		scale = as.numeric(0.5 * t(eta) %*% Kinv %*% eta)
		sig2K_new = 1 / rgamma(1, r/2 + hyper$a_sig2K, hyper$b_sig2K + scale)
		if (!fixed$sig2K) { sig2K = sig2K_new }
		timer$sig2K = timer$sig2K + as.numeric(Sys.time() - st, units = "secs")

		# Draw from [eta | ---]
		st = Sys.time()
		V_eta = solve(as.matrix((1/sig2K * Kinv) + SpinvVS))
		mean_eta = V_eta %*% (SpinvV %*% (z - H %*% muB - xi))
		eta_new = mean_eta + chol(V_eta) %*% rnorm(r)
		idx = setdiff(1:r, fixed$eta)
		eta[idx] = eta_new[idx]
		timer$eta = timer$eta + as.numeric(Sys.time() - st, units = "secs")

		# Draw from [muB | ---] using sparse matrix inverse
		# We use the computation below because it avoids a large matrix-matrix multiplcation
		st = Sys.time()
		PostLam = 1 / (eig_HpinvVH$values + 1/sig2mu)
		if (FALSE) {
			# Simpler version
			V_muB = eig_HpinvVH$vectors %*% (PostLam * t(eig_HpinvVH$vectors))		
			mean_muB = V_muB %*% (HpVinv %*% (z - S %*% eta - xi))
			muB_new = chol(V_muB) %*% rnorm(n_mu) + mean_muB
		} else {
			# Version that avoids large multiplication
			mean_muB = eig_HpinvVH$vectors %*% ((PostLam * t(eig_HpinvVH$vectors)) %*%
				(HpVinv %*% (z - S %*% eta - xi)))
			muB_new = mean_muB + (eig_HpinvVH$vectors %*% (sqrt(PostLam) * rnorm(n_mu)))
		}
		idx = setdiff(1:n_mu, fixed$muB)
		muB[idx] = muB_new[idx]
		timer$muB = timer$muB + as.numeric(Sys.time() - st, units = "secs")

		# Update Y
		st = Sys.time()
		Y = H %*% muB + S %*% eta + xi
		timer$Y = timer$Y + as.numeric(Sys.time() - st, units = "secs")

		# Save history
		if ((tt > burn) & (tt %% thin == 0)) {
			tt_keep = tt_keep + 1
			muB_hist[tt_keep,] = as.numeric(muB)
			eta_hist[tt_keep,] = as.numeric(eta)
			xi_hist[tt_keep,] = as.numeric(xi)
			sig2mu_hist[tt_keep] = sig2mu
			sig2xi_hist[tt_keep] = sig2xi
			sig2K_hist[tt_keep] = sig2K
			Y_hist[tt_keep,] = as.numeric(Y)
		}
	}

	logger("Finished Gibbs sampler\n")

	ret = list(muB_hist = muB_hist,
		xi_hist = xi_hist,
		eta_hist = eta_hist,
		sig2mu_hist = sig2mu_hist,
		sig2xi_hist = sig2xi_hist,
		sig2K_hist = sig2K_hist,
		Y_hist = Y_hist,
		z = z, H = H, S = S, v = v, R_keep = R_keep,
		elapsed_sec = timer
	)
	class(ret) = "stcos_gibbs"

	st = Sys.time()
	ret$loglik = logLik(ret)
	ret$dic = DIC(ret)
	timer$post = timer$post + as.numeric(Sys.time() - st, units = "secs")

	return(ret)
}

#' @name gibbs_stcos
#' @export
logLik.stcos_gibbs = function(object, ...)
{
	if (!is.null(object$loglik)) {
		return(object$loglik)
	}

	R_keep = object$R_keep
	loglik_mcmc = numeric(R_keep)
	for (r in 1:R_keep) {
		muB = object$muB_hist[r,]
		eta = object$eta_hist[r,]
		sig2xi = object$sig2xi_hist[r]
		loglik_mcmc[r] = sum(
			dnorm(object$z, as.numeric(object$H %*% muB + object$S %*% eta),
			sqrt(object$v + sig2xi), log = TRUE)
		)
	}
	return(loglik_mcmc)
}

#' Deviance Information Criterion
#' 
#' Generic function to calculate Deviance Information Criterion (DIC) for a
#' given model object.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments.
#'
#' @return A numeric value of the DIC.
#' @export
DIC = function(object, ...)
{
	UseMethod("DIC")
}

#' @name gibbs_stcos
#' @method DIC stcos_gibbs
#' @export
DIC.stcos_gibbs = function(object, ...)
{
	if (!is.null(object$dic)) {
		return(object$dic)
	}

	loglik_mcmc = logLik(object)
	muB_bar = colMeans(object$muB_hist)
	eta_bar = colMeans(object$eta_hist)
	sig2xi_bar = mean(object$sig2xi_hist)
	D_thetabar = -2 * sum(
		dnorm(object$z, as.numeric(object$H %*% muB_bar + object$S %*% eta_bar),
		sqrt(object$v + sig2xi_bar),
		log = TRUE))
	D_bar = mean(-2*loglik_mcmc)
	D_thetabar + 2*(D_bar - D_thetabar)
}

#' @name gibbs_stcos
#' @method print stcos_gibbs
#' @export
print.stcos_gibbs = function(x, ...)
{
	variances_mcmc = cbind(x$sig2mu_hist, x$sig2K_hist, x$sig2xi_hist)
	colnames(variances_mcmc) = c("sig2mu", "sig2K", "sig2xi")
	summary_mcmc = data.frame(
		apply(variances_mcmc, 2, mean),
		apply(variances_mcmc, 2, sd),
		apply(variances_mcmc, 2, quantile, probs = 0.025),
		apply(variances_mcmc, 2, quantile, probs = 0.25),
		apply(variances_mcmc, 2, quantile, probs = 0.75),
		apply(variances_mcmc, 2, quantile, probs = 0.975)
	)
	colnames(summary_mcmc) = c("Mean", "SD", "2.5%", "25%", "75%", "97.5%")

	total_sec = sum(unlist(x$elapsed_sec))
	hh = floor(total_sec / 60^2)
	mm = floor((total_sec - hh*60^2) / 60)
	ss = round(total_sec - hh*60^2 - mm*60)

	printf("Fit for STCOS model\n")
	printf("--\n")
	print(summary_mcmc)
	printf("--\n")
	printf("Saved %d draws\n", x$R_keep)
	printf("DIC: %f\n", x$dic)
	printf("Elapsed time: %02d:%02d:%02d\n", hh, mm, ss)

	invisible(x)
}

#' @name gibbs_stcos
#' @export
fitted.stcos_gibbs = function(object, H, S, ...)
{
	R_keep = object$R_keep
	n = nrow(H)
	E_mcmc = matrix(NA, R_keep, n)
	for (r in 1:R_keep) {
		muB = object$muB_hist[r,]
		eta = object$eta_hist[r,]
		E_mcmc[r,] = as.numeric(H %*% muB + S %*% eta)
	}
	return(E_mcmc)
}

# @method predict stcos_gibbs
#' @name gibbs_stcos
#' @export
predict.stcos_gibbs = function(object, H, S, ...)
{
	R_keep = object$R_keep
	n = nrow(H)
	Y_mcmc = matrix(NA, R_keep, n)
	for (r in 1:R_keep) {
		muB = object$muB_hist[r,]
		eta = object$eta_hist[r,]
		sig2xi = object$sig2xi_hist[r]
		Y_mcmc[r,] = rnorm(n, as.numeric(H %*% muB + S %*% eta), sqrt(sig2xi))
	}
	return(Y_mcmc)
}
