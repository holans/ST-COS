#' MLE for STCOS Model
#'
#' @param z Vector which represents the outcome; assumed to be direct
#'   estimates from the survey.
#' @param v Vector which represents direct variance estimates from the survey.
#'   The diagonal of the matrix \eqn{\bm{V}} described in the details.
#' @param H Matrix of overlaps between source and fine-level supports.
#' @param S Design matrix for basis decomposition.
#' @param K Variance of the random coefficient \eqn{\bm{\eta}}
#' @param init A list containing the initial values in the MCMC for
#'   \code{sig2xi} and \code{sig2K}. If not specified, we select an
#'   arbitrary initial value.
#' @param optim_control This is passed as the \code{control} argument to
#'   \code{optim}. Note that the value \code{fnscale} is ignored if
#'   specified.
#' @param optim_method Method to be used for likelihood maximization by
#'   \code{optim}. Default is \code{L-BFGS-B}.
#'
#' @return A list containing maximum likelihood estimates.
#'
#' @details Maximize the likelihood of the STCOS model
#' \deqn{
#'   f(\bm{z} \mid \bm{\mu}_B, \sigma_K^2, \sigma_\xi^2)
#'   = \textrm{N}(\bm{z} \mid \bm{H} \bm{\mu}_B, \bm{\Delta}
#'   ),
#'   \quad \bm{\Delta} = \sigma_\xi^2 \bm{I} + \bm{V} + \sigma_K^2 \bm{S} \bm{K} \bm{S}^\top,
#' }
#' by numerical maximization of the profile likelihood
#' \deqn{
#'   \ell(\sigma_K^2, \sigma_\xi^2) =
#'   -\frac{N}{2} \log(2 \pi) -\frac{1}{2} \log |\bm{\Delta}| -\frac{1}{2} (\bm{z} -
#'   \bm{H} \hat{\bm{\mu}}_B)^\top \bm{\Delta}^{-1} (\bm{z} - \bm{H} \hat{\bm{\mu}}_B)
#' }
#' using
#' \eqn{
#'   \hat{\bm{\mu}}_B = (\bm{H}^\top \bm{\Delta}^{-1} \bm{H})^{-1}
#'   \bm{H}^\top \bm{\Delta}^{-1} \bm{z}.
#' }
#'
#' @examples
#' \dontrun{
#' demo = prepare_stcos_demo()
#' mle_out = mle_stcos(demo$z, demo$v, demo$S, demo$H, demo$K)
#' sig2K_hat = mle_out$sig2K_hat
#' sig2xi_hat = mle_out$sig2xi_hat
#' mu_hat = mle_out$mu_hat
#' }
#' @name mle_stcos
#' @export
mle_stcos = function(z, v, H, S, K, init = NULL,
	optim_control = list(), optim_method = "L-BFGS-B")
{
	N = nrow(H)
	n = ncol(H)
	r = ncol(S)

	# Initial values
	if (is.null(init)) { init = list() }
	if (is.null(init$sig2K)) { init$sig2K = 1 }
	if (is.null(init$sig2xi)) { init$sig2xi = 1 }

	SKST = as.matrix(S %*% K %*% t(S))
	loglik = function(par) {
		sig2K = exp(par[1])
		sig2xi = exp(par[2])

		Delta = sig2K * SKST + Diagonal(x = sig2xi + v)
		mu_hat = wls(z, H, Delta)
		z_hat = H %*% mu_hat

		logdet = as.numeric( determinant(Delta)$modulus )
		ll = -N/2 * log(2*pi) - logdet / 2 - sum((z - z_hat) * solve(Delta, z - z_hat)) / 2
		return(ll)
	}

	optim_control$fnscale = -1
	st = Sys.time()
	par_init = c(log(init$sig2K), log(init$sig2xi))
	res = optim(par = par_init, loglik, method = optim_method,
		control = optim_control)
	elapsed_sec = as.numeric(Sys.time() - st, unit = "secs")

	sig2K_hat = exp(res$par[1])
	sig2xi_hat = exp(res$par[2])
	Delta = sig2K_hat * SKST + Diagonal(x = sig2xi_hat + v)
	mu_hat = wls(z, H, Delta)
	z_hat = H %*% mu_hat

	list(sig2K_hat = sig2K_hat, sig2xi_hat = sig2xi_hat,
		mu_hat = as.numeric(mu_hat), res = res, elapsed_sec = elapsed_sec,
		loglik = res$value)
}

# For the model y = X beta + eps, eps ~ N(0, Sigma),
# Compute the Weighted Least Squares estimate for beta, with Sigma given.
# Also compute y_hat = X beta_hat
wls = function(y, X, Sigma)
{
	Sigma_inv_X = solve(Sigma, X)
	pinv(as.matrix(t(X) %*% Sigma_inv_X)) %*% (t(Sigma_inv_X) %*% y)
}
