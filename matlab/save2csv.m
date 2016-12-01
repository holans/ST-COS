% Sparse matrix format
[I,J] = find(H);
csvwrite('H_sparse.txt', [I, J, nonzeros(H)]);

% Sparse matrix format
[I,J] = find(S1);
csvwrite('S1_sparse.txt', [I, J, nonzeros(S1)]);

csvwrite('Zagg.txt', Zagg);
csvwrite('sigmavar.txt', sigmavar);
csvwrite('Kinv.txt', Kinv);

