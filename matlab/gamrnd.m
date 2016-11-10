% Gamma(alpha, beta) generator using Marsaglia and Tsang method
% Algorithm 4.33. This code was copied from
% http://www.hongliangjie.com/2012/12/19/how-to-generate-gamma-random-variables/
function x = gamrnd(alpha, beta)
if alpha > 1
	d = alpha-1/3;
	c = 1/sqrt(9*d);
	flag = 1;
	while flag
		Z = randn(1);
		if Z > -1/c
			V = (1+c*Z)^3;
			U = rand;
			flag = log(U) > (0.5*Z^2+d-d*V+d*log(V));
		end
	end
	x = d*V/beta;
else
	x = gamrnd(alpha+1,beta);
	x = x * rand^(1/alpha);
end

