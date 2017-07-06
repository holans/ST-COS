// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat sample_Gamma0_VAR1(const arma::mat& A, const arma::mat& Sigma,
	long unsigned int R, long unsigned int burn = 10000)
{
	size_t m = Sigma.n_rows;
	arma::mat Sigma_half = arma::chol(Sigma);
	arma::vec y_t_minus1 = Sigma_half * arma::randn(m);
	arma::vec y_t;
	arma::vec z;
	arma::mat Gamma0(m, m);
	Gamma0.fill(0);

	for (size_t t = 0; t < burn; t++) {
		z = Sigma_half * arma::randn(m);
		y_t = A * y_t_minus1 + z;
		y_t_minus1 = y_t;
	}
	
	for (size_t t = 0; t < R; t++) {
		z = Sigma_half * arma::randn(m);
		y_t = A * y_t_minus1 + z;
		y_t_minus1 = y_t;
		Gamma0 += (y_t * trans(y_t)) / R;
	}
	
	return Gamma0;
}

