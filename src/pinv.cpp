// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat pinv(const arma::mat& X)
{
	return arma::pinv(X);
}
