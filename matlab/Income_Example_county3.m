clear all

openlogger('output.log');
logger('About to load input data');

%load in everything from data_organization
load prep_sptm_alldataH_7815;

logger('Input data loaded, now running MCMC');
 
betahat = pinv(X0'*X0)*X0'*Zagg;

burn = 0;
thin = 1;

%MCMC
[Y, S,eta, xi, xi2, sig2xi,sig2xi2, acceptW, varW,lambda_eta] = Met_spt_COS_xiismean(Zagg,10000,X0,S1,sigmavar,Kinv,LamHpinvVH,EigHinvVHp,HpVinv,H,burn,thin);

logger('MCMC finished, saving output');

%save MCMC results
save new_sptimeCOS_r250_alldata_stpxi_7815 -v7.3 Y S xi2 eta xi sig2xi acceptW varW lambda_eta;

logger('Output saved, finishing');
closelogger();
