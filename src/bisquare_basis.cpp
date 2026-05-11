#include <RcppArmadillo.h>
#include <vector>

// [[Rcpp::export]]
arma::sp_mat compute_basis_sp(const arma::mat& X, const arma::mat& cc, double w)
{
	unsigned int N = X.n_rows;
	unsigned int r = cc.n_rows;

	std::vector<unsigned int> ind_row;
	std::vector<unsigned int> ind_col;
	std::vector<double> vals;

	double w2 = w*w;

	for (unsigned int j = 0; j < r; j++) {
		double cc_x = cc(j,0);
		double cc_y = cc(j,1);

		for (unsigned int i = 0; i < N; i++) {
			double ds1 = X(i,0) - cc_x;
			double ds2 = X(i,1) - cc_y;
			double norm2 = ds1*ds1 + ds2*ds2;

			if (norm2 <= w2) {
				double root_dist = 1 - norm2 / w2;
				ind_row.push_back(i);
				ind_col.push_back(j);
				vals.push_back(root_dist * root_dist);
			}
		}
	}

	// Batch insertion constructor for sparse matrix
	arma::umat loc(2, vals.size());
	loc.row(0) = arma::conv_to<arma::urowvec>::from(ind_row);
	loc.row(1) = arma::conv_to<arma::urowvec>::from(ind_col);
	arma::sp_mat out(loc, arma::vec(vals), N, r);

	return out;
}

// [[Rcpp::export]]
arma::sp_mat compute_basis_spt(const arma::mat& X, const arma::mat& cc,
	double w_s, double w_t)
{
	unsigned int N = X.n_rows;
	unsigned int r = cc.n_rows;

	std::vector<unsigned int> ind_row;
	std::vector<unsigned int> ind_col;
	std::vector<double> vals;

	double w2_s = w_s * w_s;
	double w2_t = w_t * w_t;

	for (unsigned int j = 0; j < r; j++) {
		double cc_x = cc(j,0);
		double cc_y = cc(j,1);
		double cc_t = cc(j,2);

		for (unsigned int i = 0; i < N; i++) {
			double ds1 = X(i,0) - cc_x;
			double ds2 = X(i,1) - cc_y;
			double dt = X(i,2) - cc_t;
			double norm2_s = ds1*ds1 + ds2*ds2;
			double norm2_t = dt*dt;

			if (norm2_s < w2_s && norm2_t < w2_t) {
				double root_dist = 2 - norm2_s / w2_s - norm2_t / w2_t;
				ind_row.push_back(i);
				ind_col.push_back(j);
				vals.push_back(root_dist * root_dist / 4);
			}
		}
	}

	// Batch insertion constructor for sparse matrix
	arma::umat loc(2, vals.size());
	loc.row(0) = arma::conv_to<arma::urowvec>::from(ind_row);
	loc.row(1) = arma::conv_to<arma::urowvec>::from(ind_col);
	arma::sp_mat out(loc, arma::vec(vals), N, r);

	return out;
}
