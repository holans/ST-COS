library(stcos)

cc = 1e3 * matrix(c(
#   -0.1103,    0.0470,    2.0125,
-0.1104,    0.0469,    2.0122,
-0.0783,    0.0361,    2.0125,
-0.0910,    0.0363,    2.0125,
-0.0986,    0.0315,    2.0125,
-0.1147,    0.0424,    2.0130
), 5, 3, byrow = TRUE)

w.s = 1.5
w.t = 1.5

# Jon's code computes basis with rl instead of w.s
# Use type 1 quantile algorithm to match Matlab 
G <- dist(cc)
rl <- w.s * quantile(G[G > 0], prob = 0.05, type = 1)

X = matrix(
	c(-110.3, 47.0, 2012.5,
	  -111.3, 46.0, 2012.5),
	2, 3, byrow = TRUE)

S <- compute_basis(X, cc, rl, w.t)
print(S)
