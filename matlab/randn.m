% Draw from standard normal without the Statistics Toolbox.
% Uses the Box-Muller transform.
function x = randn(M)
    u1 = rand(M);
    u2 = rand(M);
    r1 = -2 .* log(1 - u1);
    r2 = (1 - u2) .* 2*pi;
    x = sqrt(r1) .* cos(r2);
end

