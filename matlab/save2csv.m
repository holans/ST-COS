% Sparse matrix format
[I,J] = find(H);
csvwrite('H_sparse.txt', [I, J, nonzeros(H)]);

% Sparse matrix format
[I,J] = find(S);
csvwrite('S_sparse.txt', [I, J, nonzeros(S)]);

csvwrite('Zagg.txt', Zagg);
csvwrite('sigmavar.txt', sigmavar);
csvwrite('Kinv.txt', Kinv);

