#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::sp_mat compute_basis(const arma::mat& X, const arma::mat& cc, double w_s, double w_t)
{
	size_t N = X.n_rows;
	size_t r = cc.n_rows;
	arma::sp_mat S(N, r);

	double w2_s = pow(w_s, 2.0);
	double w2_t = pow(w_t, 2.0);

	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < r; j++) {
			// Rprintf("compute_basis: i=%d, j=%d\n", i, j);
			double norm2_s = pow(X(i,0) - cc(j,0), 2.0) + pow(X(i,1) - cc(j,1), 2.0);
			double norm2_t = pow(X(i,2) - cc(j,2), 2.0);
			if (norm2_s < w2_s && norm2_t < w2_t) {
				double root_dist = 1 - norm2_s / w2_s - norm2_t / w2_t;
				S(i,j) = pow(root_dist, 2.0);
			}
		}
	}

	return S;
}

