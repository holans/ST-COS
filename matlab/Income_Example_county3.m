clear all

%load in everything from data_organization
load prep_sptm_alldataH_7815;
 
betahat = pinv(X0'*X0)*X0'*Zagg;

%MCMC
[Y, S, eta, mu_B, xi, sig2mu, sig2xi, sig2K] = Met_spt_COS_xiismean(Zagg,50,X0,S1,sigmavar,Kinv,LamHpinvVH,EigHinvVHp,HpVinv,H);

%save MCMC results
save new_sptimeCOS_r250_alldata_stpxi_7815 -v7.3 Y S xi eta mu_B sig2xi sig2K;

