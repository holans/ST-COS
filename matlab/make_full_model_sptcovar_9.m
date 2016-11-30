function K=make_full_model_sptcovar_9(Qinv,M,S,n1)

FullCovar = zeros(size(S,1),size(S,1));
for i = 1:9
    for j = 1:9
        postmult = eye(n1);
        if i==j 
            FullCovar(n1*(i-1)+1:n1*i,n1*(i-1)+1:n1*i) = 0.5*Qinv;
        elseif j>i
                for k = (i+1):j
                    postmult = M*postmult; 
                end
         FullCovar(n1*(i-1)+1:n1*i,n1*(j-1)+1:n1*j) = (postmult*Qinv)*postmult';
        end
    end
	disp(i);
end

FullCovar = FullCovar' + FullCovar;
 
SpS = S'*S;
SpSinv = pinv(SpS);
SpSinvSp = SpSinv*S';

K = (SpSinvSp*FullCovar)*SpSinvSp';

end

