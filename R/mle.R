#' MLE for STCOS Model
#'
#' @param z Vector which represents the outcome; assumed to be direct
#'        estimates from the survey.
#' @param v Vector which represents direct variance estimates from the survey.
#' @param S Design matrix for basis decomposition.
#' @param H Matrix of overlaps between source and fine-level supports.
#' @param init A list containing the following initial values for the MCMC:
#' 	      \code{sig2xi}. If not specified, we select an arbitrary initial
#' 	      value.
#' @param optim.control This is passed as the \code{control} argument to
#'        \code{optim}. Note that the value \code{fnscale} is ignored if
#'        specified.
#'
#' @return A list containing maximum likelihood estimates.
#' @export
#'
#' @examples
mle.stcos <- function(z, v, H, S, K, init = NULL,
	optim.control = list())
{
	n <- ncol(H)
	r <- ncol(S)

	# Initial values
	if (is.null(init)) { init <- list() }
	if (is.null(init$sig2K)) { init$sig2K <- 1 }
	if (is.null(init$sig2xi)) { init$sig2xi <- 1 }

	loglik <- function(theta) {
		sig2K <- exp(theta[1])
		sig2xi <- exp(theta[2])

		Sigma <- sig2K * (S %*% K %*% t(S)) + diag(x = sig2xi + v)
		Sigma.inv.H <- solve(Sigma, H)
		mu.hat <- ginv(as.matrix(t(H) %*% Sigma.inv.H)) %*% (t(Sigma.inv.H) %*% z)
		# mu.hat <- as.numeric(solve(t(H) %*% Sigma.inv.H, t(Sigma.inv.H) %*% z))

		logdet <- determinant(Sigma)
		z.star <- z - H %*% mu.hat
		as.numeric(-n/2 * log(2*pi) - (1/2) * logdet$modulus -
			(1/2) * t(z.star) %*% solve(Sigma, z.star))
	}

	optim.control$fnscale <- -1
	st <- Sys.time()
	par.init <- c(log(init$sig2K), log(init$sig2xi))
	res <- optim(par = par.init, loglik, method = "L-BFGS-B",
		control = optim.control)
	elapsed.sec <- as.numeric(Sys.time() - st, unit = "secs")

	sig2K.hat <- exp(res$par[1])
	sig2xi.hat <- exp(res$par[2])
	Sigma.hat <- (sig2K.hat * S %*% K %*% t(S)) + diag(x = sig2xi.hat + v)
	Sigma.inv.H.hat <- solve(Sigma.hat, H)
	# mu.hat <- as.numeric(solve(t(H) %*% Sigma.inv.H.hat, t(Sigma.inv.H.hat) %*% z))
	mu.hat <- ginv(as.matrix(t(H) %*% Sigma.inv.H.hat)) %*% (t(Sigma.inv.H.hat) %*% z)

	list(sig2K.hat = sig2K.hat, sig2xi.hat = sig2xi.hat,
		mu.hat = mu.hat, res = res, elapsed.sec = elapsed.sec)
}
