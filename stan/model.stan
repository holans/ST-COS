data {
	int<lower=0> N;
	int<lower=0> n;
	int<lower=0> r;
	vector[N] z;
	matrix[N,n] H;
	matrix[N,r] S;
	vector[N] sig2eps;
}
parameters {
	vector[n] mu;
	vector[r] eta;
	// vector[N] xi;
	real<lower=0> sig2mu;
	real<lower=0> sig2xi;
	real<lower=0> sig2K;
}
model {
	vector[N] Hmu;
	vector[N] Seta;

	// prior
	for (j in 1:n) {
		mu[j] ~ normal(0, sqrt(sig2mu));
	}
	for (j in 1:r) {
		eta[j] ~ normal(0, sqrt(sig2K));
	}
	sig2mu ~ inv_gamma(1, 1);
	sig2xi ~ inv_gamma(1, 1);
	sig2K ~ inv_gamma(1, 1);

	// log-likelihood
	Hmu = H*mu;
	Seta = S*eta;
	for (i in 1:N) {
		// print("z[i]=", z[i], " Hmu[i]=", Hmu[i], " Seta[i]=", Seta[i], " sig2eps[i]=", sig2eps[i]);
		z[i] ~ normal(Hmu[i] + Seta[i], sqrt(sig2eps[i] + sig2xi));
	}

	print("sig2mu=", sig2mu, " sig2xi=", sig2xi, " sig2K=", sig2K, " log-likelihood=", target());
}
