load prep_sptm_alldataH_7815.mat

% AMR: Not sure exactly what I should put here, but
% trying to match the pattern of what Jon did for his
% basis function computation.
% t2 = [level2, 2012*ones(size(level2,1),1)];
% SGBF2 = ArealBi2_spacetime(county2, 2013, t2, [], [], 100, 0.5, 0.5);

% This returns a number of columns incompatible with existing S matrix
% SGBF2 = ArealBi2_spacetime(county2, 2013, level, [], [], 100, 0.5, 0.5);
% [SGBF2_1, idx] = licols(SGBF2);

% Let's try picking out the same columns as we originally did with S
% SGBF2_1 = SGBF2(:,idx);


SGBF3 = ArealBi2_spacetime(county3, 2013, level, [], [], 100, 0.5, 0.5);
SGBF3_1 = SGBF3(:,idx);

[I,J] = find(SGBF3_1);
csvwrite('basis_2013_1_sparse.dat', [I, J, nonzeros(SGBF3_1)]);

% Also write the overlap matrix H2, which is needed for predictions
% [I,J] = find(H2);
% csvwrite('H_2013_2_sparse.dat', [I, J, nonzeros(H2)]);
