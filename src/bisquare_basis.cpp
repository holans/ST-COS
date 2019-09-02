// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>

// As of the last time I checked, the batch insertion method for arma::sp_mat
// still appears to be a bit faster than initializing S to zeros and then
// setting non-zero elements. We'll keep using batch insertion for now, and
// revisit later...

// [[Rcpp::export]]
arma::sp_mat compute_basis_sp(const arma::mat& X, const arma::mat& cc, double w)
{
	unsigned int N = X.n_rows;
	unsigned int r = cc.n_rows;

	std::vector<unsigned int> ind_row;
	std::vector<unsigned int> ind_col;
	std::vector<double> vals;

	double w2 = w*w;

	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < r; j++) {
			double ds1 = X.at(i,0) - cc.at(j,0);
			double ds2 = X.at(i,1) - cc.at(j,1);
			double norm2 = ds1*ds1 + ds2*ds2;

			if (norm2 <= w2) {
				double root_dist = 1 - norm2 / w2;
				ind_row.push_back(i);
				ind_col.push_back(j);
				vals.push_back(root_dist * root_dist);
			}
		}
	}

	arma::vec values = arma::conv_to<arma::vec>::from(vals);
	arma::umat locations(2, vals.size());
	locations.row(0) = arma::conv_to<arma::uvec>::from(ind_row).t();
	locations.row(1) = arma::conv_to<arma::uvec>::from(ind_col).t();

	arma::sp_mat S(locations, values, N, r);
	return S;
}

// [[Rcpp::export]]
arma::sp_mat compute_basis_spt(const arma::mat& X, const arma::mat& cc, double w_s, double w_t)
{
	unsigned int N = X.n_rows;
	unsigned int r = cc.n_rows;

	std::vector<unsigned int> ind_row;
	std::vector<unsigned int> ind_col;
	std::vector<double> vals;

	double w2_s = w_s * w_s;
	double w2_t = w_t * w_t;

	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < r; j++) {
			double ds1 = X.at(i,0) - cc.at(j,0);
			double ds2 = X.at(i,1) - cc.at(j,1);
			double dt = X.at(i,2) - cc.at(j,2);
			double norm2_s = ds1*ds1 + ds2*ds2;
			double norm2_t = dt*dt;

			if (norm2_s <= w2_s && norm2_t <= w2_t) {
				double root_dist = 2 - norm2_s / w2_s - norm2_t / w2_t;
				ind_row.push_back(i);
				ind_col.push_back(j);
				vals.push_back(root_dist * root_dist / 4);
			}
		}
	}

	arma::vec values = arma::conv_to<arma::vec>::from(vals);
	arma::umat locations(2, vals.size());
	locations.row(0) = arma::conv_to<arma::uvec>::from(ind_row).t();
	locations.row(1) = arma::conv_to<arma::uvec>::from(ind_col).t();

	arma::sp_mat S(locations, values, N, r);
	return S;
}
