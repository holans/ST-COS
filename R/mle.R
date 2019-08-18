#' MLE for STCOS Model
#'
#' @param z Vector which represents the outcome; assumed to be direct
#'        estimates from the survey.
#' @param v Vector which represents direct variance estimates from the survey.
#' @param H Matrix of overlaps between source and fine-level supports.
#' @param S Design matrix for basis decomposition.
#' @param K_inv Inverse of the \eqn{K} matrix, which is the covariance of the
#'        random coefficient \eqn{\eta}
#' @param init A list containing the following initial values for the MCMC:
#' 	      \code{sig2xi}. If not specified, we select an arbitrary initial
#' 	      value.
#' @param optim_control This is passed as the \code{control} argument to
#'        \code{optim}. Note that the value \code{fnscale} is ignored if
#'        specified.
#'
#' @return A list containing maximum likelihood estimates.
#'
#' @examples
#' \dontrun{
#' z = sp$get_z()
#' v = sp$get_v()
#' H = sp$get_H()
#' S.reduced = sp$get_reduced_S() 
#' K_inv = sp$get_Kinv(2005:2015)
#' 
#' mle.out = mle_stcos(z, v, S.reduced, H, solve(K_inv))
#' 
#' sig2K_hat = mle.out$sig2K_hat
#' sig2xi_hat = mle.out$sig2xi_hat
#' mu_hat = mle.out$mu_hat
#' }
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

	loglik = function(par) {
		sig2K = exp(par[1])
		sig2xi = exp(par[2])

		Delta = sig2K * S %*% K %*% t(S) + Diagonal(x = sig2xi + v)
		Delta_inv_H = solve(Delta, H)
		mu_hat = pinv(as.matrix(t(H) %*% Delta_inv_H)) %*% (t(Delta_inv_H) %*% z)
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
	Delta_hat = sig2K_hat * S %*% K %*% t(S) + Diagonal(x = sig2xi_hat + v)
	Delta_inv_H_hat = solve(Delta_hat, H)
	mu_hat = pinv(as.matrix(t(H) %*% Delta_inv_H_hat)) %*% (t(Delta_inv_H_hat) %*% z)

	list(sig2K_hat = sig2K_hat, sig2xi_hat = sig2xi_hat,
		mu_hat = as.numeric(mu_hat), res = res, elapsed_sec = elapsed_sec)
}
