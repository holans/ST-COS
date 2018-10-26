// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>

// [[Rcpp::export]]
arma::mat pinv(const arma::mat& X)
{
	return arma::pinv(X);
}
