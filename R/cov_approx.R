#' Best Approximation to Covariance Structure
#' 
#' Compute the best positive approximant for use in the STCOS
#' model, under several prespecified covariance structures.
#' 
#' @param Delta Covariance (\eqn{n \times n}) for observations within a time
#'     point for the process whose variance we wish to approximate.
#' @param S Design matrix (\eqn{N \times r}) of basis functions evaluated on
#'     the fine-level process over \eqn{T = N / n} time points.
#' 
#' @details
#' Let \eqn{\bm{\Sigma}} be an \eqn{N \times N} symmetric and positive-definite
#' covariance matrix and \eqn{\bm{S}} be an \eqn{N \times r} matrix with
#' rank \eqn{r}. The objective is to compute a matrix \eqn{\bm{K}} which minimizes
#' the Frobenius norm
#' \deqn{
#'   \Vert \bm{\Sigma} - \bm{S} \bm{C} \bm{S}^\top {\Vert}_\textrm{F},
#' }
#' over symmetric positive-definite matrices \eqn{\bm{C}}. The
#' solution is given by
#' \deqn{
#'   \bm{K} = (\bm{S}^\top \bm{S})^{-1} \bm{S}^\top \bm{\Sigma} \bm{S} (\bm{S}^\top \bm{S})^{-1}.
#' }
#' 
#' In the STCOS model, \eqn{\bm{S}} represents the design matrix from a basis
#' function computed from a fine-level support having \eqn{n} areas, using
#' \eqn{T} time steps. Therefore \eqn{N = n T} represents the dimension of
#' covariance for the fine-level support.
#' 
#' We provide functions to handle some possible structures for target
#' covariance matrices of the form
#' \deqn{
#'   \bm{\Sigma} =
#'   \left(
#'   \begin{array}{ccc}
#'   \bm{\Gamma}(1,1) & \cdots & \bm{\Gamma}(1,T) \\
#'   \vdots           & \ddots & \vdots \\
#'   \bm{\Gamma}(T,1) & \cdots & \bm{\Gamma}(T,T)
#'   \end{array}
#'   \right),
#' }
#' where each \eqn{\bm{\Gamma}(s,t)} is an \eqn{n \times n} matrix.
#' 
#' \itemize{
#' \item \code{cov_approx_randwalk} assumes \eqn{\bm{\Sigma}} is based on the
#' autocovariance function of a random walk
#' \deqn{
#'   \bm{Y}_{t+1} = \bm{Y}_{t} + \bm{\epsilon}_t, \quad \bm{\epsilon}_t \sim \textrm{N}(\bm{0}, \bm{\Delta}).
#' }
#' so that
#' \deqn{
#'   \bm{\Gamma}(s,t) = \min(s,t) \bm{\Delta}.
#' }
#' 
#' \item \code{cov_approx_blockdiag} assumes \eqn{\bm{\Sigma}} is based on
#' \deqn{
#'   \bm{Y}_{t+1} = \bm{Y}_{t} + \bm{\epsilon}_t, \quad \bm{\epsilon}_t \sim \textrm{N}(\bm{0}, \bm{\Delta}).
#' }
#' which are independent across \eqn{t}, so that
#' \deqn{
#'   \bm{\Gamma}(s,t) = I(s = t) \bm{\Delta},
#' }
#' }
#' The block structure is used to reduce the computational burden, as \eqn{N}
#' may be large.
#' 
#' @name Covariance Approximation
#' 
#' @export
cov_approx_randwalk = function(Delta, S)
{
	n = nrow(Delta)
	N = nrow(S)
	r = ncol(S)
	lag_max = N / n

	SpS = t(S) %*% S
	if (rankMatrix(SpS) < r) {
		warning("The matrix (S' S) is rank-deficient. Consider reducing the dimension of S")
	}
	SpSinv = pinv(as.matrix(SpS))

	C = matrix(0, r, r)
	for (i in 1:lag_max) {
		for (j in 1:lag_max) {
			idx_i = seq(n*(i-1)+1, n*i)
			idx_j = seq(n*(j-1)+1, n*j)
			C = C + t(S[idx_i,]) %*% (min(i,j)*Delta) %*% S[idx_j,]
		}
	}

	K = SpSinv %*% C %*% SpSinv
	return((K + t(K)) / 2)
}

#' @name Covariance Approximation
#' @export
cov_approx_blockdiag = function(Delta, S)
{
	n = nrow(Delta)
	N = nrow(S)
	r = ncol(S)
	lag_max = N / n

	SpS = t(S) %*% S
	if (rankMatrix(SpS) < r) {
		warning("The matrix (S' S) is rank-deficient. Consider reducing the dimension of S")
	}
	SpSinv = pinv(as.matrix(SpS))
	
	C = matrix(0, r, r)
	for (i in 1:lag_max) {
		idx = seq(n*(i-1)+1, n*i)
		C = C + t(S[idx,]) %*% Delta %*% S[idx,]
	}

	K = SpSinv %*% C %*% SpSinv
	return((K + t(K)) / 2)
}

# This function is currently experimental, and not exported.
# @param M Propagator matrix (\eqn{n \times n}) for VAR(1) covariance
#     structure.
# \item \code{cov_approx_moran} assumes \eqn{\bm{\Sigma}} is based on the
# autocovariance function of a VAR(1) process
# \deqn{
#   \bm{Y}_{t+1} = \bm{M} \bm{Y}_{t} + \bm{\epsilon}_t, \quad \bm{\epsilon}_t \sim \textrm{N}(\bm{0}, \bm{\Delta}).
# }
# using a given \eqn{\bm{M}}, so that
# \deqn{
#   \bm{\Gamma}(t+h,t) = \bm{M}^h \bm{\Gamma}(0), \quad \textrm{if $h \geq 0$},
# }
# \deqn{
#   \bm{\Gamma}(t+h,t) = [\bm{M}^{-h} \bm{\Gamma}(0)]^\top, \quad \textrm{if $h < 0$},
# }
# where
# \deqn{
#   \textrm{vec}(\bm{\Gamma}(0)) = [\bm{I} - \bm{M} \otimes \bm{M}]^{-1} \textrm{vec}(\bm{\Delta}).
# }
cov_approx_varone = function(Delta, M, S)
{
	n = nrow(Delta)
	N = nrow(S)
	r = ncol(S)
	lag_max = N / n

	SpS = t(S) %*% S
	if (rankMatrix(SpS) < r) {
		warning("The matrix (S' S) is rank-deficient. Consider reducing the dimension of S")
	}
	SpSinv = pinv(as.matrix(SpS))

	# Get all the autocovariances we'll need
	G = autocov_VAR1(M, Delta, lag_max)
	Gamma = function(h) {
		if (h >= 0) {
			G[,,h+1]
		} else {
			t(G[,,-h+1])
		}
	}

	C = matrix(0, r, r)
	for (i in 1:lag_max) {
		for (j in 1:lag_max) {
			idx_i = seq(n*(i-1)+1, n*i)
			idx_j = seq(n*(j-1)+1, n*j)
			C = C + t(S[idx_i,]) %*% Gamma(i-j) %*% S[idx_j,]
		}
	}

	K = SpSinv %*% C %*% SpSinv
	return((K + t(K)) / 2)
}
