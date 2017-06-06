SpatialBisquareBasis <- R6Class("SpatialBisquareBasis",
	public = list(
		initialize = function(cutpoints.x, cutpoints.y, w) {
			r <- length(cutpoints.x)
			stopifnot(length(cutpoints.y) == r)
			private$cutpoints <- cbind(cutpoints.x, cutpoints.y)
			private$w <- w
			private$r <- r

			# Jon's code computes basis with rl instead of w.s
			# Use type 1 quantile algorithm to match Matlab 
			G <- dist(private$cutpoints)
			private$rl <- w * quantile(G[G > 0], prob = 0.05, type = 1)
			printf("rl = %f\n", private$rl)
		},
		get_dim = function() {
			private$r
		}
	),
	private = list(
		r = NULL,
		cutpoints = NULL,
		w = NULL,
		rl = NULL
	)
)

compute <- function(x, y)
{
	X <- cbind(x, y)
	cc <- private$cutpoints
	S <- compute_basis_sp(X, cc, private$rl)
	return(S)
}

SpatialBisquareBasis$set("public", "compute", compute)
SpatialBisquareBasis$lock()

