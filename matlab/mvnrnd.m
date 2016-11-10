function x = mvnrnd(mu, Sigma)
    [p,q] = size(Sigma);
    xt = chol(Sigma) * randn([p 1]) + mu';
    x = xt';
end

