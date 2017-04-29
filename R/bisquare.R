BisquareBasis <- R6Class("BisquareBasis",
	public = list(
		initialize = function(cutpoints.x, cutpoints.y, cutpoints.t, w.s, w.t) {
			r <- length(cutpoints.x)
			stopifnot(length(cutpoints.y) == r)
			stopifnot(length(cutpoints.t) == r)
			private$cutpoints <- cbind(cutpoints.x, cutpoints.y, cutpoints.t)
			private$w.s <- w.s
			private$w.t <- w.t
			private$r <- r
		},
		get_dim = function() {
			private$r
		}
	),
	private = list(
		r = NULL,
		cutpoints = NULL,
		w.s = NULL,
		w.t = NULL
	)
)

compute <- function(x, y, time)
{
	X <- cbind(x, y, time)
	cc <- private$cutpoints
	w <- c(private$w.s, private$w.s, private$w.t)
	S <- compute_basis(X, cc, w)
	return(S)
}

compute.old <- function(x, y, time)
{
	X <- cbind(x, y, time)
	N <- nrow(X)

	cc <- private$cutpoints
	w2 <- c(private$w.s, private$w.s, private$w.t)^2
	r <- private$r
	S <- Matrix(0, N, r)

	for (i in 1:N) {
		for (j in 1:r) {
			h2 <- (X[i,] - cc[j,])^2
			d2 <- 1 - sum(h2 / w2)
			if (d2 > 0) {
				S[i,j] <- sqrt(d2)
			}
		}
	}

	return(S)
}

BisquareBasis$set("public", "compute", compute)
BisquareBasis$lock()
