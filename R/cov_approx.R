#' Best Approximation to Covariance Structure
#' 
#' Compute the best positive approximant for use in the STCOS
#' model, under several prespecified covariance structures.
#' 
#' @param Qinv
#' @param X
#' @param S
#' @param lag_max The number of lags to compute.
#' 
#' @details
#' Let \eqn{\bm{\Sigma}} be an \eqn{N \times N} symmetric and
#' positive-definite covariance matrix, which we would like to approximate
#' in the STCOS model. Let \eqn{\bm{S}}
#' be the \eqn{N \times r} design matrix of the basis function computed on \code{lag_max}
#' time steps of the fine-level support.
#' 
#' The objective is to compute a symmetric positive-definite matrix
#' \eqn{\bm{K}} which minimizes
#' \deqn{
#'   || \bm{\Sigma} - \bm{S} \bm{C} \bm{S}^\top {||}_\textrm{F}
#' }
#' This is given by 
#' \deqn{
#'   \bm{K} = (\bm{S}^\top \bm{S})^{-1} \bm{S}^\top \bm{\Sigma} \bm{S} (\bm{S}^\top \bm{S})^{-1} 
#' }
#' Sigma is the autocovariance of VAR(1), the first lag_max blocks
#' We don't need to construct Sigma in its entirety to do this operation.
#' @name Covariance Approximation
#' @export

#' @name Covariance Approximation
#' @export
cov_approx_moran = function(Qinv, X, S, lag_max)
{
	n = nrow(Qinv)
	r = ncol(S)

	SpS = t(S) %*% S
	if (rankMatrix(SpS) < r) {
		warning("The matrix (S' S) is rank-deficient. Consider reducing the dimension of S")
	}
	SpSinv = pinv(as.matrix(SpS))

	# Create VAR coefficient matrix from X
	P_perp = Diagonal(n,1) - X %*% solve(t(X) %*% X, t(X))
	eig = eigen(P_perp, symmetric = TRUE)
	M = Re(eig$vectors)
	M = (M + t(M)) / 2

	# Get all the autocovariances we'll need
	G = autocov_VAR1(M, Qinv, lag_max)
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
cov_approx_randwalk = function(Qinv, S, lag_max)
{
	n = nrow(Qinv)
	r = ncol(S)

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
			C = C + t(S[idx_i,]) %*% (min(i,j)*Qinv) %*% S[idx_j,]
		}
	}

	K = SpSinv %*% C %*% SpSinv
	return((K + t(K)) / 2)
}

#' @name Covariance Approximation
#' @export
cov_approx_blockdiag = function(Qinv, S, lag_max)
{
	n = nrow(Qinv)
	r = ncol(S)

	SpS = t(S) %*% S
	if (rankMatrix(SpS) < r) {
		warning("The matrix (S' S) is rank-deficient. Consider reducing the dimension of S")
	}
	SpSinv = pinv(as.matrix(SpS))
	
	C = matrix(0, r, r)
	for (i in 1:lag_max) {
		idx = seq(n*(i-1)+1, n*i)
		C = C + t(S[idx,]) %*% Qinv %*% S[idx,]
	}

	K = SpSinv %*% C %*% SpSinv
	return((K + t(K)) / 2)
}
