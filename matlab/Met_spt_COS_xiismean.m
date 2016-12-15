function [Y, S,eta, xi, xi2, sig2xi,sig2xi2, acceptW, varW,lambda_eta] = Met_spt_COS_xiismean(Z,T,X,S,sig2eps,Kinv,LamHpVinvH,EigHpVinvV,HpVinv,H)
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
bb=50;
b=0;
tt = 1;
while tt < T 
tt = tt + 1;

msg = sprintf('Starting iteration %d', tt);
logger(msg);

if mod(tt,bb)==0
    b=1;
else
    b=0;
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

PostLaminv = diag(LamHpVinvH) + (1/sig2xi(tt-1))*ones(size(LamHpVinvH,1),1);
PostLam = 1./PostLaminv;
PostCov = EigHpVinvV * ((PostLam * ones(1, nxi)) .* EigHpVinvV');
muxi = (PostCov*(HpVinv*(Z -X*beta(:,tt)- S*eta(:,tt-1) - xi2(:,tt-1))));
xispatial = EigHpVinvV*(sqrt(PostLam).*randn([size(HpVinv,1),1]));
xi(:,tt) = muxi + xispatial;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Full Conditional sig2xi
% stopxi2nan = 0;
% stopxi2inf = 0;

bsh = 0.000001;
shapesig2xi = 0.5*(xi(:,tt)'*xi(:,tt));
sig2xi(tt)=real(1/gamrnd((nxi/2)+1,(1/(bsh+shapesig2xi))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Full Conditional Eta
zresid = Z-X*beta(:,tt)- H*xi(:,tt) - xi2(:,tt-1);

% tempteta = pinv(large_frac(lambda_eta(tt-1))*Kinv + SpinvV*S);
% choletaforsim = cholinv_approxs2(tempteta,1e-4);
tempteta = inv((1/lambda_eta(tt-1))*Kinv + SpinvV*S);
choletaforsim = chol(tempteta);
mueta = tempteta*(SpinvV*zresid);
eta(:,tt) =real(mueta + choletaforsim*mvnrnd(zeros(1,r),eye(r))');

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Full Conditional Xi
%sparse matrix inverse

PostLaminv = Vinv + (1/sig2xi2(tt-1))*ones(n,1);
PostLam = (1./PostLaminv);
muxi1 = (Vinv.*(Z-X*beta(:,tt) - S*eta(:,tt-1) - H*xi(:,tt)));
muxi =PostLam.*muxi1;
xispatial = sqrt(PostLam).*real(randn([n,1]));
xi2(:,tt) = muxi + xispatial;

% % 
%    shapesig2xi2 = 0.5*(xi2(:,tt)'*xi2(:,tt)); 
%    sig2xi2(tt)=real(1/gamrnd((n/2)+1,(1/(0.000001+shapesig2xi2))));

%Metropolis Step for Weight Matrix
% [W(:,tt), varW, tunew,acceptW]= metropolis_block_AOAS_W(2,tt,acceptW,varW,b,W,eta(:,tt),lambda_eta(tt-1),r,LambdaPrior,De);
%     
% WW=inverse_logit_design(W(:,tt),De,r);
% PhiQ = GivensRotation(WW);
% Kinv = PhiQ*diag(LambdaPrior)*PhiQ/lambda_eta(tt-1);

scale = 0.5*eta(:,tt)'*Kinv*eta(:,tt);
lambda_eta(tt)=1/gamrnd((r/2 +1),(1/(0.000001+scale)));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y(:,tt) = X*beta(:,tt)+S*eta(:,tt) + H*xi(:,tt) + xi2(:,tt);


if mod(tt,1)==0
disp(tt);
% disp([acceptW/tt, varW]);
end

end



end
