# Compute the covariance matrix for a VAR(1) process.
# Sigma is covariance of error term.
# A is VAR(1) coefficient matrix.
#' @export
autocov_VAR1 = function(A, Sigma, lag_max)
{
	A = as.matrix(A)
	m = nrow(Sigma)
	N = m * (lag_max+1)
	Gamma = array(NA, dim = c(m, m, lag_max+1))

	# The simple Kronecker product calculation for lag-0 is infeasible for large m,
	# even if A is fairly sparse.
	# Gamma[,,1] = solve(Diagonal(m^2,1) - (A %x% A), matrix(Sigma, m^2, 1))
	# The following version avoids the Kronecker product.
	eig = eigen(A)
	V = eig$vectors
	lambda = eig$values
	if (any(lambda > 1 - 1e-20)) {
		warning("Unit root detected in autocov_VAR1")
	}
	rm(eig)
	if (isSymmetric(V)) {
		V.inv = t(V)
	} else {
		V.inv = pinv(V)
	}
	C = V.inv %*% Sigma %*% t(V.inv)
	rm(V.inv)
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
