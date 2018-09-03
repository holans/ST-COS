#' Spatial Bisquare Basis
#' 
#' An R6 class for ...
#'
#' @export
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#'
#' @examples
SpatialBisquareBasis <- R6Class("SpatialBisquareBasis",
	private = list(
		r = NULL,
		cutpoints = NULL,
		w = NULL,
		rl = NULL
	),
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
			private$rl <- as.numeric(w * quantile(G[G > 0], prob = 0.05, type = 1))
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
		compute = function(x, y) {
			X <- cbind(x, y)
			S <- compute_basis_sp(X, private$cutpoints, private$rl)
			return(S)
		}
	)
)

# SpatialBisquareBasis$set("public", "compute", compute)
# SpatialBisquareBasis$lock()

