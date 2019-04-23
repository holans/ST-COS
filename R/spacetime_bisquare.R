#' SpaceTime Bisquare Basis
#' 
#' An \code{\link{R6Class}} representing the space-time bisquare basis.
#' 
#' @section Usage:
#' \preformatted{
#' basis <- SpaceTimeBisquareBasis$new(knots.x, knots.y, knots.t,
#'     w.s, w.t)
#' basis$compute(x, y, time)
#' basis$get_dim()
#' basis$get_cutpoints()
#' basis$get_rl()
#' basis$get_ws()
#' basis$get_wt()
#' }
#' 
#' @section Arguments:
#' \itemize{
#' \item \code{knots.x} x-coordinate of knot points.
#' \item \code{knots.y} y-coordinate of knot points.
#' \item \code{knots.t} time coordinate of knot points.
#' \item \code{w.s} (Original, before transformation) spatial radius.
#' \item \code{w.t} Temporal radius.
#' \item \code{x} Vector of x-coordinates for points on which to evaluate the basis.
#' \item \code{y} Vector of y-coordinates for points on which to evaluate the basis.
#' \item \code{time} Vector of time coordinates for points on which to evaluate the basis.
#' }
#' 
#' @section Methods:
#' \itemize{
#' \item \code{new} Create a new \code{SpaceTimeBisquareBasis} object.
#' \item \code{get_dim} Get the number of cutpoints used to construct this basis.
#' \item \code{get_cutpoints} Get the cutpoints used to construct this basis.
#' \item \code{get_rl} Get the transformed spatial radius. The transformation is
#' based on quantiles of distances between knots.
#' \item \code{get_ws} Get the original spatial radius used to construct this
#' basis. This is transformed before it is applied to account for the geography
#' being used.
#' \item \code{get_wt} Get the temporal radius used to construct this basis.
#' \item \code{compute} Evaluate this basis on specific points.
#' }
#'
#' @name SpaceTimeBisquareBasis
#' 
#' @examples
#' set.seed(1234)
#' seq.x <- seq(0, 1, length.out = 3)
#' seq.y <- seq(0, 1, length.out = 3)
#' seq.t <- seq(0, 1, length.out = 3)
#' knots = expand.grid(seq.x, seq.y, seq.t)
#' x <- runif(50)
#' y <- runif(50)
#' t <- sample(1:3, size = 50, replace = TRUE)
#' 
#' basis <- SpaceTimeBisquareBasis$new(knots[,1], knots[,2], knots[,3], w.s = 0.5, w.t = 1)
#' basis$compute(x, y, t)
#' basis$get_dim()
#' basis$get_cutpoints()
#' basis$get_rl()
#' basis$get_ws()
#' basis$get_wt()
#' 
#' # Plot the (spatial) knots and the points at which we evaluated the basis
#' plot(knots[,1], knots[,2], pch = 4, cex = 1.5, col = "red")
#' text(x, y, labels = t, cex = 0.75)
#' 
#' # Draw a circle representing the basis' radius around one of the knot points
#' tseq <- seq(0, 2*pi, length=100) 
#' rad <- basis$get_rl()
#' coords <- cbind(rad * cos(tseq) + seq.x[2], rad * sin(tseq) + seq.y[2])
#' lines(coords, col = "red")
NULL

#' @export
#' @docType class
SpaceTimeBisquareBasis <- R6Class("SpaceTimeBisquareBasis",
	lock_objects = TRUE,
	lock_class = TRUE,
	private = list(
		r = NULL,
		cutpoints = NULL,
		w.s = NULL,
		w.t = NULL,
		rl = NULL
	),
	public = list(
		initialize = function(knots.x, knots.y, knots.t, w.s, w.t) {
			r <- length(knots.x)
			stopifnot(length(knots.y) == r)
			stopifnot(length(knots.t) == r)
			private$cutpoints <- cbind(knots.x, knots.y, knots.t)
			private$w.s <- w.s
			private$w.t <- w.t
			private$r <- r

			# Jon's code computes basis with rl instead of w.s
			# Use type 1 quantile algorithm to match Matlab 
			G <- dist(private$cutpoints)
			private$rl <- as.numeric(w.s * quantile(G[G > 0], prob = 0.05, type = 1))
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
