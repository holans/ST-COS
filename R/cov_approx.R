#' Best Approximation to Covariance Structure
#' 
#' Compute the best positive approximant for use in the STCOS
#' model, under several prespecified covariance structures.
#' 
#' @param Delta Covariance (\eqn{n \times n}) for observations within a time
#'     point for the process whose variance we wish to approximate.
#' @param M Propagator matrix (\eqn{n \times n}) for VAR(1) covariance
#'     structure.
#' @param S Design matrix (\eqn{N \times r}) of basis functions evaluated on
#'     the fine-level process over \eqn{T = N / n} time points.
#' 
#' @details
#' Let \eqn{\bm{\Sigma}} be an \eqn{N \times N} symmetric and positive-definite
#' covariance matrix, which we would like to approximate in the STCOS model.
#' Let \eqn{\bm{S}} be the \eqn{N \times r} design matrix of the basis
#' function computed on \eqn{T} time steps of the fine-level support.
#' The number of observations at each time point is assumed to be \eqn{n},
#' so that the number of total observations is \eqn{N = n T}.
#' 
#' The objective is to compute a symmetric positive-definite matrix
#' \eqn{\bm{K}} which minimizes
#' \deqn{
#'   || \bm{\Sigma} - \bm{S} \bm{C} \bm{S}^\top {||}_\textrm{F},
#' }
#' where \eqn{|| \cdot {||}_\textrm{F}} represents the Frobenius norm. The
#' solution is given by
#' \deqn{
#'   \bm{K} = (\bm{S}^\top \bm{S})^{-1} \bm{S}^\top \bm{\Sigma} \bm{S} (\bm{S}^\top \bm{S})^{-1}.
#' }
#' 
#' We provide functions to handle three possible structures for the target
#' covariance, which are all in the form
#' \deqn{
#'   \bm{\Sigma} =
#'   \left(
#'   \begin{array}{ccc}
#'   \bm{\Gamma}(1,1) & \cdots & \bm{\Gamma}(1,T) \\
#'   \vdots           & \ddots & \vdots \\
#'   \bm{\Gamma}(T,1) & \cdots & \bm{\Gamma}(T,T)
#'   \end{array}
#'   \right).
#' }
#' 
#' \itemize{
#' \item \code{cov_approx_moran} assumes \eqn{\bm{\Sigma}} is based on the
#' autocovariance function of a VAR(1) process
#' \deqn{
#'   \bm{Y}_{t+1} = \bm{M} \bm{Y}_{t} + \bm{\epsilon}_t, \quad \bm{\epsilon}_t \sim \textrm{N}(\bm{0}, \bm{\Delta}).
#' }
#' using a given \eqn{\bm{M}}, so that
#' \deqn{
#'   \bm{\Gamma}(t+h,t) = \bm{M}^h \bm{\Gamma}(0), \quad \textrm{if $h \geq 0$},
#' }
#' \deqn{
#'   \bm{\Gamma}(t+h,t) = [\bm{M}^{-h} \bm{\Gamma}(0)]^\top, \quad \textrm{if $h < 0$},
#' }
#' where
#' \deqn{
#'   \textrm{vec}(\bm{\Gamma}(0)) = [\bm{I} - \bm{M} \otimes \bm{M}]^{-1} \textrm{vec}(\bm{\Delta}).
#' }
#' 
#' \item \code{cov_approx_randwalk} assumes \eqn{\bm{\Sigma}} is based on the
#' autocovariance function of a random walk
#' \deqn{
#'   \bm{Y}_{t+1} = \bm{Y}_{t} + \bm{\epsilon}_t, \quad \bm{\epsilon}_t \sim \textrm{N}(\bm{0}, \bm{\Delta}).
#' }
#' so that
#' \deqn{
#'   \bm{\Gamma}(t+h,t) = t \bm{\Delta}, \quad \textrm{if $h \geq 0$},
#' }
#' \deqn{
#'   \bm{\Gamma}(t+h,t) = (t - |h|) \bm{\Delta}, \quad \textrm{if $-t < h < 0$}.
#' }
#' 
#' \item \code{cov_approx_blockdiag} assumes \eqn{\bm{\Sigma}} is based on
#' \deqn{
#'   \bm{Y}_{t+1} = \bm{Y}_{t} + \bm{\epsilon}_t, \quad \bm{\epsilon}_t \sim \textrm{N}(\bm{0}, \bm{\Delta}).
#' }
#' which are independent across \eqn{t}, so that
#' \deqn{
#'   \bm{\Gamma}(s,s) = \bm{\Delta},
#' }
#' \deqn{
#'   \bm{\Gamma}(s,t) = \bm{0}, \quad \textrm{if $s \neq t$}.
#' }
#' }
#' In any case \eqn{\bm{\Sigma}} may be large and we avoid computing it in its entirety.
#' @name Covariance Approximation
#' @export

#' @name Covariance Approximation
#' @export
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

#' @name Covariance Approximation
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
