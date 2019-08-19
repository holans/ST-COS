/*
 * Spatio-Temporal Change of Support model, adapted from <doi:10.1002/sta4.94>.
 */
data {
	int<lower=0> N;   // number of total obs
	int<lower=0> n;   // number of areas in fine-level domain
	int<lower=0> r;   // number of random coefficients
	vector[N] z;      // direct estimates
	vector[N] v;      // direct variance estimates
	matrix[N,n] H;    // design matrix for fixed coeffs
	matrix[N,r] S;    // design matrix for random coeffs
	matrix[r,r] K;
	real alpha_K;
	real beta_K;
	real alpha_xi;
	real beta_xi;
	real alpha_mu;
	real beta_mu;
}
parameters {
	vector[n] mu;
	vector[r] eta;
	vector[N] xi;
	real<lower=0> sig2K;
	real<lower=0> sig2xi;
	real<lower=0> sig2mu;
}
model {
	sig2K ~ inv_gamma(alpha_K, beta_K);
	sig2xi ~ inv_gamma(alpha_xi, beta_xi);
	sig2mu ~ inv_gamma(alpha_mu, beta_mu);
	eta ~ multi_normal(rep_vector(0,r), sig2K * K);
	mu ~ normal(0, sqrt(sig2mu));
	xi ~ normal(0, sqrt(sig2xi));
	z ~ normal(to_vector(H*mu + S*eta + xi), sqrt(v));
}
