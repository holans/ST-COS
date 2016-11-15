function [Y, S,eta, xi, xi2, sig2xi,sig2xi2, acceptW, varW,lambda_eta] = Met_spt_COS_xiismean(Z,T,X,S,sig2eps,Kinv,LamHpVinvH,EigHpVinvV,HpVinv,H)

p = size(X,2);
r = size(S,2);
n = size(Z,1);
nxi = size(HpVinv,1);

xi = zeros(nxi,T);
xi2 = zeros(n,T);

eta = zeros(r,T);
eta(:,1) = randn([r,1]);
Lambda = ones(r,T);
sig2xi = ones(1,T);
sig2xi2 = ones(1,T);
lambda_eta =ones(T,1);
beta = zeros(p,T);
W = zeros(2,T);

for j = 1:r
	SpinvV(j,:) = S(:,j)./sig2eps;
end

Vinv = 1./sig2eps;

%% Initialize the Gibbs sampler
%% Start sampling 

acceptW = 0;
varW = 1;
bb = 50;
b = 0;
tt = 1;
while tt < T 
	tt = tt + 1;

	if mod(tt,bb) == 0
		b = 1;
	else
		b = 0;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Full Conditional Xi
	%sparse matrix inverse
	PostLaminv = diag(LamHpVinvH) + (1/sig2xi(tt))*ones(size(LamHpVinvH,1),1);
	PostLam = diag(1./PostLaminv);
	PostCov = EigHpVinvV*PostLam*EigHpVinvV';

	muxi1 = (PostCov*(HpVinv*(Z -X*beta(:,tt)- S*eta(:,tt-1) - xi2(:,tt-1))));
	muxi2 = PostLam*muxi1;
	muxi = EigHpVinvV*muxi2;

	xispatial1 = sqrt(PostLam)*real(randn([size(HpVinv,1),1]));
	xispatial = EigHpVinvV*xispatial1;
	xi(:,tt) = muxi + xispatial;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Full Conditional sig2xi
	bsh = 0.000001;
	shapesig2xi = 0.5*(xi(:,tt)'*xi(:,tt)); 
	sig2xi(tt) = real(1/gamrnd((nxi/2)+1,(1/(bsh+shapesig2xi))));
	   
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Full Conditional Eta
	zresid = Z-X*beta(:,tt)- H*xi(:,tt) - xi2(:,tt-1);

	tempteta = inv((1/lambda_eta(tt-1))*Kinv + SpinvV*S);
	choletaforsim = chol(tempteta);
	mueta = tempteta*SpinvV*zresid;
	eta(:,tt) = real(mueta + choletaforsim*mvnrnd(zeros(r,1),eye(r))');

	%Stop if eta has a problem
	if sum(isnan(eta(:,tt)))>0
	   disp('eta problem NaN');
	   break
	end
	if sum(isinf(eta(:,tt)))>0
	   disp('eta problem Inf');
	   break
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Full Conditional Xi
	%sparse matrix inverse

	PostLaminv = Vinv + (1/sig2xi2(tt-1))*ones(n,1);
	PostLam = (1./PostLaminv);

	muxi1 = (Vinv.*(Z-X*beta(:,tt) - S*eta(:,tt-1) - H*xi(:,tt)));
	muxi = PostLam.*muxi1;

	xispatial = sqrt(PostLam).*real(randn([n,1]));
	xi2(:,tt) = muxi + xispatial;
	 
	scale = 0.5*eta(:,tt)'*Kinv*eta(:,tt);
	lambda_eta(tt) = 1/gamrnd((r/2 +1),(1/(0.000001+scale)));

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Y(:,tt) = X*beta(:,tt)+S*eta(:,tt) + H*xi(:,tt) + xi2(:,tt);

	if mod(tt,1)==0
		disp(tt);
		% disp([acceptW/tt, varW]);
	end

end

end
