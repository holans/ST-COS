load prep_sptm_alldataH_7815;

% Moran's I Propagator
B = [H, eye(3109)];
[M,LM] = eig(B*pinv(B'*B)*B');
M = real(M);

% Target Covariance
Kinv = make_full_model_sptcovar_9(Q,M,Sconnectorf,3109);
[P,D] = eig(Kinv);
D = real(diag(D));
D(D<0) = 0;
Dinv = D;
Dinv(D>0) = (1./D(D>0));
K = real(P)*diag(Dinv)*real(P');
Kinv = real(P)*diag(D)*real(P');

csvwrite('Kinv_amr.txt', Kinv);
