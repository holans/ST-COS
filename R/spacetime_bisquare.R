SpaceTimeBisquareBasis <- R6Class("SpaceTimeBisquareBasis",
	private = list(
		r = NULL,
		cutpoints = NULL,
		w.s = NULL,
		w.t = NULL,
		rl = NULL
	),
	public = list(
		initialize = function(cutpoints.x, cutpoints.y, cutpoints.t, w.s, w.t) {
			r <- length(cutpoints.x)
			stopifnot(length(cutpoints.y) == r)
			stopifnot(length(cutpoints.t) == r)
			private$cutpoints <- cbind(cutpoints.x, cutpoints.y, cutpoints.t)
			private$w.s <- w.s
			private$w.t <- w.t
			private$r <- r

			# Jon's code computes basis with rl instead of w.s
			# Use type 1 quantile algorithm to match Matlab 
			G <- dist(private$cutpoints)
			private$rl <- w.s * quantile(G[G > 0], prob = 0.05, type = 1)
			printf("rl = %f\n", private$rl)
		},
		get_dim = function() {
			private$r
		},
		get_cutpoints = function() {
			private$cutpoints
		},
		get_rl = function() {
			private$rl
		},
		get_ws = function() {
			private$w.s
		},
		get_wt = function() {
			private$w.t
		},
		compute = function(x, y, time) {
			X <- cbind(x, y, time)
			S <- compute_basis_spt(X, private$cutpoints, private$rl, private$w.t)
			return(S)
		}
	)
)

# SpaceTimeBisquareBasis$set("public", "compute", compute)
# SpaceTimeBisquareBasis$lock()
