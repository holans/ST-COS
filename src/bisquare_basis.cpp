#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::sp_mat compute_basis(const arma::mat& X, const arma::mat& cc, double w_s, double w_t)
{
	size_t N = X.n_rows;
	size_t r = cc.n_rows;
	arma::sp_mat S(N, r);

	double w2_s = w_s * w_s;
	double w2_t = w_t * w_t;

	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < r; j++) {
			// Rprintf("compute_basis: i=%d, j=%d\n", i, j);
			double ds1 = X(i,0) - cc(j,0);
			double ds2 = X(i,1) - cc(j,1);
			double dt = X(i,2) - cc(j,2);
			double norm2_s = ds1*ds1 + ds2*ds2;
			double norm2_t = dt*dt;
			if (norm2_s < w2_s && norm2_t < w2_t) {
				double root_dist = 1 - norm2_s / w2_s - norm2_t / w2_t;
				S(i,j) = root_dist * root_dist;
			}
		}
	}

	return S;
}

