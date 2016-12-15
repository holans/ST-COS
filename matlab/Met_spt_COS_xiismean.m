function [Y, S,eta, xi, xi2, sig2xi,sig2xi2, acceptW, varW,lambda_eta] = Met_spt_COS_xiismean(Z,T,X,S,sig2eps,Kinv,LamHpVinvH,EigHpVinvV,HpVinv,H,burn,thin)
% 
% Z = Zagg;
% H = eye(5757);
% PhiQ = real(P);
% Lambda22 = real(diag(D));
% HpinvVH = (1./sigmavar);
% HpVinv = (1./sigmavar);
% Eigit = eye(5757);
% Lamit = (1./sigmavar);
% mub = 0;
% T = 2;
% sig2eps=sigmavar;
% sigma2_beta=10^15;
% LambdaPrior = Dinv;

p=size(X,2);

r = size(S,2);
n = size(Z,1);
nxi = size(HpVinv,1);

T_keep = ceil((T - burn) / thin);
xi_hist = zeros(nxi,T_keep);
xi2_hist = zeros(n,T_keep);
eta_hist = zeros(r,T_keep);
% eta_hist(:,1) = randn([r,1]);
% Lambda = ones(r,T_keep);
sig2xi_hist = ones(1,T_keep);
sig2xi2_hist = ones(1,T_keep);
lambda_eta_hist =ones(T_keep,1);
% beta_hist = zeros(p,T_keep);
% W = zeros(2,T_keep);

% Initial values
xi = zeros(nxi,1);
xi2 = zeros(n,1);
eta = randn([r,1]);
beta = zeros(p,1);
sig2xi = 1;
sig2xi2 = 1;
lambda_eta = 1;

for j = 1:r
SpinvV(j,:) = S(:,j)./sig2eps;
end


Vinv = 1./sig2eps;
%% Initialize the Gibbs sampler
%% Start sampling 

acceptW = 0;
varW = 1;
bb=50;
b=0;
tt = 1;
tt_keep = 0;

msg = sprintf('Using burn = %d and thin = %d', burn, thin);
logger(msg);

while tt < T 
tt = tt + 1;

if mod(tt,bb)==0
    b=1;
else
    b=0;
end

if mod(tt,1) == 0
	msg = sprintf('Starting iteration %d', tt);
	logger(msg);
	% disp(tt);
	% disp([acceptW/tt, varW]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Full Conditional Beta

% tempt = temptbeta;
% mubeta = tempt*XpinvV*(Z - H*xi(:,tt-1) - S*eta(:,tt-1) - xi2(:,tt-1)) + tempt*mub/sigma2_beta;
% beta(:,tt) =real(mubeta + ctemptbeta*mvnrnd(zeros(p,1),eye(p))');
% %Stop if Beta has a problem
% if sum(isnan(beta(:,tt)))>0
%    disp('beta problem NaN');
%    break
% end
% if sum(isinf(beta(:,tt)))>0
%    disp('beta problem Inf');
%    break
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Full Conditional Xi
%sparse matrix inverse

% PostLaminv = diag(LamHpVinvH) + (1/sig2xi(tt-1))*ones(size(LamHpVinvH,1),1);
% PostLam = diag(1./PostLaminv);
% PostCov = EigHpVinvV*PostLam*EigHpVinvV';
% muxi = (PostCov*(HpVinv*(Z -X*beta(:,tt)- S*eta(:,tt-1) - xi2(:,tt-1))));
% xispatial = EigHpVinvV*(sqrt(PostLam)*randn([size(HpVinv,1),1]));
% xi(:,tt) = muxi + xispatial;

PostLaminv = diag(LamHpVinvH) + (1/sig2xi)*ones(size(LamHpVinvH,1),1);
PostLam = 1./PostLaminv;
PostCov = EigHpVinvV * ((PostLam * ones(1, nxi)) .* EigHpVinvV');
muxi = (PostCov*(HpVinv*(Z -X*beta- S*eta - xi2)));
xispatial = EigHpVinvV*(sqrt(PostLam).*randn([size(HpVinv,1),1]));
xi = muxi + xispatial;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Full Conditional sig2xi
% stopxi2nan = 0;
% stopxi2inf = 0;

bsh = 0.000001;
shapesig2xi = 0.5*(xi'*xi);
sig2xi = real(1/gamrnd((nxi/2)+1,(1/(bsh+shapesig2xi))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Full Conditional Eta
zresid = Z-X*beta- H*xi - xi2;

% tempteta = pinv(large_frac(lambda_eta(tt-1))*Kinv + SpinvV*S);
% choletaforsim = cholinv_approxs2(tempteta,1e-4);
tempteta = inv((1/lambda_eta)*Kinv + SpinvV*S);
choletaforsim = chol(tempteta);
mueta = tempteta*(SpinvV*zresid);
eta = real(mueta + choletaforsim*mvnrnd(zeros(1,r),eye(r))');

%Stop if eta has a problem
if sum(isnan(eta))>0
   disp('eta problem NaN');
   break
end
if sum(isinf(eta))>0
   disp('eta problem Inf');
   break
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Full Conditional Xi
%sparse matrix inverse

PostLaminv = Vinv + (1/sig2xi2)*ones(n,1);
PostLam = (1./PostLaminv);
muxi1 = (Vinv.*(Z-X*beta - S*eta - H*xi));
muxi =PostLam.*muxi1;
xispatial = sqrt(PostLam).*real(randn([n,1]));
xi2 = muxi + xispatial;

% % 
%    shapesig2xi2 = 0.5*(xi2(:,tt)'*xi2(:,tt)); 
%    sig2xi2(tt)=real(1/gamrnd((n/2)+1,(1/(0.000001+shapesig2xi2))));

%Metropolis Step for Weight Matrix
% [W(:,tt), varW, tunew,acceptW]= metropolis_block_AOAS_W(2,tt,acceptW,varW,b,W,eta(:,tt),lambda_eta(tt-1),r,LambdaPrior,De);
%     
% WW=inverse_logit_design(W(:,tt),De,r);
% PhiQ = GivensRotation(WW);
% Kinv = PhiQ*diag(LambdaPrior)*PhiQ/lambda_eta(tt-1);

scale = 0.5*eta'*Kinv*eta;
lambda_eta =1/gamrnd((r/2 +1),(1/(0.000001+scale)));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = X*beta+S*eta + H*xi + xi2;


if mod(tt,thin) == 0 && tt > burn
	tt_keep = tt_keep + 1;

	msg = sprintf('\tSaving item %d to history', tt_keep);
	logger(msg);

	xi_hist(:,tt_keep) = xi;
	xi2_hist(:,tt_keep) = xi2;
	eta_hist(:,tt_keep) = eta;
	sig2xi_hist(:,tt_keep) = sig2xi;
	sig2xi2_hist(:,tt_keep) = sig2xi2;
	lambda_eta_hist(tt_keep,:) = lambda_eta;
	Y_hist(:,tt_keep) = Y;
end

end



end
