// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List eigs_sym(const arma::sp_mat& X, unsigned int k)
{
	arma::vec eigval;
	arma::mat eigvec;

	arma::eigs_sym(eigval, eigvec, X, k);
	return Rcpp::List::create(
		Rcpp::Named("values") = eigval,
		Rcpp::Named("vectors") = eigvec
	);
}
