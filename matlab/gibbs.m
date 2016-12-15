function [Y, S, eta, mu_B, xi, sig2mu, sig2xi, sig2K] = Met_spt_COS_xiismean(Z, T, X, S, sig2eps, Kinv, LamHpVinvH, EigHpVinvV, HpVinv, H)

p = size(X,2);
r = size(S,2);
n = size(Z,1);
n_mu = size(HpVinv,1);

mu_B = zeros(n_mu,T);
xi = zeros(n,T);

eta = zeros(r,T);
eta(:,1) = randn([r,1]);
sig2mu = ones(1,T);
sig2xi = ones(1,T);
sig2K = ones(T,1);

for j = 1:r
	SpinvV(j,:) = S(:,j)./sig2eps;
end

Vinv = 1./sig2eps;

%% Initialize the Gibbs sampler
%% Start sampling 

for tt = 1:T
	msg = sprintf('Starting iteration %d', tt);
	logger(msg);


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Full Conditional mu_B
	% sparse matrix inverse
	PostLaminv = diag(LamHpVinvH) + (1/sig2mu(tt))*ones(size(LamHpVinvH,1),1);
	PostLam = diag(1./PostLaminv);
	PostCov = EigHpVinvV*PostLam*EigHpVinvV';

	muxi = EigHpVinvV * PostLam * PostCov * (HpVinv*(Z - S*eta(:,tt-1) - xi(:,tt-1)));
	xispatial = EigHpVinvV * sqrt(PostLam) * real(randn([size(HpVinv,1),1]));
	mu_B(:,tt) = muxi + xispatial;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Full Conditional sig2mu
	shapesig2mu = 0.5*(mu_B(:,tt)'*mu_B(:,tt)); 
	sig2mu(tt) = real(1/gamrnd((n_mu/2)+1,(1/(0.000001 + shapesig2mu))));
	   
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Full Conditional eta
	zresid = Z - H*mu_B(:,tt) - xi(:,tt-1);

	tempteta = inv((1/sig2K(tt-1))*Kinv + SpinvV*S);
	choletaforsim = chol(tempteta);
	mueta = tempteta*SpinvV*zresid;
	eta(:,tt) = real(mueta + choletaforsim*mvnrnd(zeros(r,1),eye(r))');

	% Stop if eta has a problem
	if sum(isnan(eta(:,tt)))>0
	   disp('eta problem NaN');
	   break
	end
	if sum(isinf(eta(:,tt)))>0
	   disp('eta problem Inf');
	   break
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Full Conditional xi
	% sparse matrix inverse
	PostLaminv = Vinv + (1/sig2xi(tt-1))*ones(n,1);
	PostLam = (1./PostLaminv);

	muxi = PostLam.*(Vinv.*(Z - S*eta(:,tt-1) - H*mu_B(:,tt)));
	xispatial = sqrt(PostLam).*real(randn([n,1]));
	xi(:,tt) = muxi + xispatial;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Full Conditional xi
	% shapesig2xi = 0.5*(xi(:,tt)'*xi(:,tt)); 
	% sig2xi(tt) = real(1/gamrnd((n/2)+1,(1/(0.000001 + shapesig2xi))));
 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Full Conditional sig2K
	scale = 0.5*eta(:,tt)'*Kinv*eta(:,tt);
	sig2K(tt) = 1/gamrnd((r/2 +1),(1/(0.000001 + scale)));

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Y(:,tt) = S*eta(:,tt) + H*mu_B(:,tt) + xi(:,tt);

	if mod(tt,1)==0
		disp(tt);
	end

end

end
