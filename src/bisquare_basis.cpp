#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::sp_mat compute_basis(const arma::mat& X, const arma::mat& cc, const arma::vec& w)
{
	size_t N = X.n_rows;
	size_t r = cc.n_rows;
	arma::sp_mat S(N, r);

	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < r; j++) {
			// Rprintf("compute_basis: i=%d, j=%d\n", i, j);
			arma::vec x = trans(X.row(i) - cc.row(j)) / w;
			double d2 = 1 - dot(x, x);
			if (d2 > 0) {
				S(i,j) = sqrt(d2);
			}
		}
	}

	return S;
}

