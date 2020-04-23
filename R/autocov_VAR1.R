#' Compute the autocovariance matrix for a VAR(1) process.
#' 
#' @details
#' Computes the autocovariance matrix \eqn{\bm{\Gamma}(h)} of the
#' \eqn{m}-dimensional VAR(1) process
#' \deqn{
#'   \bm{Y}_t = \bm{A} \bm{Y}_{t-1} + \bm{\epsilon}_t, \quad
#'   \bm{\epsilon}_t \sim \textrm{N}(\bm{0}, \bm{\Sigma})
#' }
#' 
#' For the required computation of \eqn{\bm{\Gamma}(0)}, this function
#' solves the \eqn{m^2 \times m^2} system
#' \deqn{
#' \textrm{vec}(\bm{\Gamma}(0)) = [\bm{I} - \bm{A} \otimes \bm{A}]^{-1} \textrm{vec}(\bm{\Sigma}).
#' }
#' without directly computing \eqn{m^2 \times m^2} matrices.
#' 
#' @param Sigma Covariance matrix  \eqn{\bm{\Sigma}} of the errors.
#' @param A Coefficient matrix \eqn{A} of the autoregression term.
#' @param lag_max maximum number of lags to compute.
#' @return An array \code{Gamma} of dimension \code{c(m, m, lag_max + 1)},
#' where the slice \code{Gamma[,,h]} represents the autocovariance at lag
#' \code{h = 0, 1, ..., lag_max}.
#' 
#' @examples
#' U = matrix(NA, 3, 3)
#' U[,1] = c(1, 1, 1) / sqrt(3)
#' U[,2] = c(1, 0, -1) / sqrt(2)
#' U[,3] = c(0, 1, -1) / sqrt(2)
#' B = U %*% diag(c(0.5, 0.2, 0.1)) %*% t(U)
#' A = (B + t(B)) / 2
#' Sigma = diag(x = 2, nrow = 3)
#' autocov_VAR1(A, Sigma, lag_max = 5)
#' 
#' @export
autocov_VAR1 = function(A, Sigma, lag_max)
{
	A = as.matrix(A)
	m = nrow(Sigma)
	N = m * (lag_max+1)
	Gamma = array(NA, dim = c(m, m, lag_max+1))

	# Avoid the Kronecker product calculation for lag-0, which is infeasible
	# for large m:
	# Gamma[,,1] = solve(Diagonal(m^2,1) - (A %x% A), matrix(Sigma, m^2, 1))
	eig = eigen(A)
	V = eig$vectors
	lambda = eig$values
	if (any(Re(lambda) > 1 - 1e-20)) {
		warning("Unit root detected in autocov_VAR1")
	}
	rm(eig)
	# The columns of V are eigenvectors of A, so V will likely not be symmetric
	if (isSymmetric(V)) {
		V_inv = t(V)
	} else {
		# V_inv = pinv(V)
		V_inv = solve(V)
	}
	# C = V_inv %*% Sigma %*% t(V_inv)
	C = V_inv %*% as.matrix(Sigma) %*% t(V_inv)
	rm(V_inv)
	e = matrix(1 / (1 - lambda %x% lambda), m^2, 1) * matrix(C, m^2, 1)
	rm(C)
	E = matrix(e, m, m)
	rm(e)
	Gamma[,,1] = Re(V %*% E %*% t(V))
	rm(E, V)

	for (h in seq_len(lag_max)) {
		Gamma[,,h+1] = A %*% Gamma[,,h]
	}

	return(Gamma)
}
