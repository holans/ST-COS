library(MCMCpack)
library(ks)

V <- cbind(
	c(0,0),
	c(0,1),
	c(1,0),
	c(1,1)
)
d <- nrow(V)
k <- ncol(V)

lambda <- rdirichlet(10000, rep(1,k))
y <- t(V %*% t(lambda))
plot(y)
fhat <- kde(y)
plot(fhat, display="filled.contour2", xlim = c(0,1), ylim = c(0,1))

u <- cbind(
	runif(10000),
	runif(10000))
plot(u)
fhat <- kde(u)
plot(fhat, display="filled.contour2", xlim = c(0,1), ylim = c(0,1))
