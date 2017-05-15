// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat compute_basis_sp(const arma::mat& X, const arma::mat& cc, double w)
{
	size_t N = X.n_rows;
	size_t r = cc.n_rows;
	arma::mat S(N, r);

	double w2 = w*w;

	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < r; j++) {
			// Rprintf("compute_basis: i=%d, j=%d\n", i, j);
			double ds1 = X(i,0) - cc(j,0);
			double ds2 = X(i,1) - cc(j,1);
			double norm2_s = ds1*ds1 + ds2*ds2;

			double root_dist = 1 - norm2_s / w2;
			S(i,j) = root_dist * root_dist * (norm2_s <= w2);
		}
	}

	return S;
}

// [[Rcpp::export]]
arma::mat compute_basis_spt(const arma::mat& X, const arma::mat& cc, double w_s, double w_t)
{
	// We could use arma::sp_mat here, but it slows down very badly if we set
	// elements the naive way: S(i,j) = el; If memory use becomes a problem, we
	// could try to use batch insertion methods e.g. a row at a time.

	size_t N = X.n_rows;
	size_t r = cc.n_rows;
	arma::mat S(N, r);

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

			double root_dist = 1 - norm2_s / w2_s - norm2_t / w2_t;
			S(i,j) = root_dist * root_dist * (norm2_s <= w2_s && norm2_t <= w2_t);
		}
	}

	return S;
}

