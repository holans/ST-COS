# Return K which minimizes || Sigma - S X S^T ||_F over symmetric psd matrices X
# where Sigma is symmetric and pd. This gives
#   K = SpSinvSp %*% Sigma %*% t(SpSinvSp)
# where
#   SpSinvSp = solve(t(S) %*% S) %*% t(S)
#   Sigma is the autocovariance of VAR(1), the first lag_max blocks
# We don't need to construct Sigma in its entirety to do this operation.
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
