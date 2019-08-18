#' Areal Spatial Bisquare Basis
#' 
#' An \code{\link{R6Class}} representing the spatial bisquare basis.
#' 
#' @section Usage:
#' \preformatted{
#' basis <- SpatialBisquareBasis$new(knots_x, knots_y, w)
#' basis$compute(x, y, time)
#' basis$get_dim()
#' basis$get_cutpoints()
#' basis$get_rl()
#' basis$get_w()
#' }
#' 
#' @section Arguments:
#' \itemize{
#' \item \code{knots_x} x-coordinate of knot points.
#' \item \code{knots_y} y-coordinate of knot points.
#' \item \code{w} (Original, before transformation) radius.
#' \item \code{x} Vector of x-coordinates for points on which to evaluate the basis.
#' \item \code{y} Vector of y-coordinates for points on which to evaluate the basis.
#' }
#' 
#' @section Methods:
#' \itemize{
#' \item \code{new} Create a new \code{SpatialBisquareBasis} object.
#' \item \code{get_dim} Get the number of cutpoints used to construct this basis.
#' \item \code{get_cutpoints} Get the cutpoints used to construct this basis.
#' \item \code{get_rl} Get the transformed spatial radius. The transformation is
#' based on quantiles of distances between knots.
#' \item \code{get_w} Get the original radius used to construct this
#' basis. This is transformed before it is applied to account for the geography
#' being used.
#' \item \code{compute} Evaluate this basis on specific points.
#' }
#'
#' @name ArealSpatialBisquareBasis
#' 
#' @examples
#' set.seed(1234)
#' seq_x <- seq(0, 1, length.out = 3)
#' seq_y <- seq(0, 1, length.out = 3)
#' knots = merge(seq_x, seq_y)
#' x <- runif(50)
#' y <- runif(50)
#' 
#' basis <- SpatialBisquareBasis$new(knots[,1], knots[,2], w = 0.5)
#' basis$compute(x, y)
#' basis$get_dim()
#' basis$get_cutpoints()
#' basis$get_rl()
#' basis$get_w()
#' 
#' # Plot the knots and the points at which we evaluated the basis
#' plot(knots[,1], knots[,2], pch = 4, cex = 1.5, col = "red")
#' points(x, y, cex = 0.5)
#' 
#' # Draw a circle representing the basis' radius around one of the knot points
#' tseq <- seq(0, 2*pi, length=100) 
#' rad <- basis$get_rl()
#' coords <- cbind(rad * cos(tseq) + seq_x[2], rad * sin(tseq) + seq_y[2])
#' lines(coords, col = "red")
NULL

#' @export
#' @docType class
ArealSpatialBisquareBasis = R6Class("ArealSpatialBisquareBasis",
	lock_objects = TRUE,
	lock_class = TRUE,
	private = list(
		mc_reps = NULL,
		report_period = NULL,
		basis_sp = NULL
	),
	public = list(
		initialize = function(knots_x, knots_y, w, mc_reps) {
			private$mc_reps = mc_reps
			private$basis_sp = SpatialBisquareBasis$new(knots_x, knots_y, w = w)
		},
		get_mc_reps = function() {
			private$mc_reps
		},
		get_basis_sp = function() {
			private$basis_sp
		},
		get_dim = function() {
			private$basis_sp$r
		},
		get_cutpoints = function() {
			private$basis_sp$cutpoints
		},
		get_rl = function() {
			private$basis_sp$rl
		},
		get_w = function() {
			private$basis_sp$w
		},
		compute = function(dom, report_period = nrow(dom) + 1) {
			basis = private$basis_sp
			R = private$mc_reps

			n = nrow(dom)
			r = basis$get_dim()
			S = Matrix(0, n, r)

			for (j in 1:n) {
				if (j %% report_period == 0) {
					logger("Computing basis for area %d of %d\n", j, n)
				}

				# Request a few more samples than we'll need, to prevent the loop in rArea.
				P = rArea(R, dom[j,], blocksize = ceiling(1.2*R))
				S[j,] = S[j,] + colSums(basis$compute(P[,1], P[,2]))
			}

			return(S / R)
		}
	)
)
