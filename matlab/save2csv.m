% Sparse matrix format
[I,J] = find(H);
csvwrite('H_sparse.txt', [I, J, nonzeros(H)]);

[I,J] = find(S1);
csvwrite('S1_sparse.txt', [I, J, nonzeros(S1)]);

[I,J] = find(HpinvVH);
csvwrite('HpinvVH_sparse.txt', [I, J, nonzeros(HpinvVH)]);

[I,J] = find(EigHinvVHp);
csvwrite('EigHinvVHp_sparse.txt', [I, J, nonzeros(EigHinvVHp)]);

[I,J] = find(Sconnectorf);
csvwrite('Sconnectorf_sparse.txt', [I, J, nonzeros(Sconnectorf)]);

csvwrite('Zagg.txt', Zagg);
csvwrite('sigmavar.txt', sigmavar);
csvwrite('Kinv.txt', Kinv);
csvwrite('LamHpinvVH.txt', diag(LamHpinvVH));
csvwrite('M.txt', M);
