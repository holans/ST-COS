#' Areal Spatial Bisquare Basis
#' 
#' An \code{R6Class} representing the spatial bisquare basis.
#' 
#' @section Usage:
#' \preformatted{
#' bs = ArealSpatialBisquareBasis$new(knots_x, knots_y, w, mc_reps)
#' bs$compute(dom)
#' bs$get_dim()
#' bs$get_cutpoints()
#' bs$get_w()
#' }
#' 
#' @section Arguments:
#' \itemize{
#' \item \code{knots_x} numeric vector; x-coordinates of knot points.
#' \item \code{knots_y} numeric vector; y-coordinates of knot points.
#' \item \code{w} numeric; radius for the basis.
#' \item \code{dom} an \code{sf} object; areal units to evaluate.
#' }
#' 
#' @section Methods:
#' \itemize{
#' \item \code{new} Create a new \code{ArealSpatialBisquareBasis} object.
#' \item \code{get_dim} Get the number of cutpoints used to construct this basis.
#' \item \code{get_cutpoints} Get the cutpoints used to construct this basis.
#' \item \code{get_w} Get the radius used to construct this basis.
#' \item \code{compute} Evaluate this basis on specific areal units.
#' }
#'
#' @name ArealSpatialBisquareBasis
#' 
#' @examples
#' set.seed(1234)
#' seq_x = seq(0, 1, length.out = 3)
#' seq_y = seq(0, 1, length.out = 3)
#' knots = merge(seq_x, seq_y)
#' 
#' # Create a simple domain from rectangles
#' shape1 = matrix(c(0.0,0.0, 0.5,0.0, 0.5,0.5, 0.0,0.5, 0.0,0.0), ncol=2, byrow=TRUE)
#' shape2 = shape1 + cbind(rep(0.5,5), rep(0.0,5))
#' shape3 = shape1 + cbind(rep(0.0,5), rep(0.5,5))
#' shape4 = shape1 + cbind(rep(0.5,5), rep(0.5,5))
#' sfc = st_sfc(
#'    st_polygon(list(shape1)),
#'    st_polygon(list(shape2)),
#'    st_polygon(list(shape3)),
#'    st_polygon(list(shape4))
#' )
#' dom = st_sf(data.frame(geoid = 1:length(sfc), geom = sfc))
#' 
#' bs = ArealSpatialBisquareBasis$new(knots[,1], knots[,2], w = 0.5,
#'     mc_reps = 200)
#' bs$compute(dom)
#' bs$get_dim()
#' bs$get_cutpoints()
#' bs$get_w()
#' 
#' # Plot the knots and the points at which we evaluated the basis
#' plot(knots[,1], knots[,2], pch = 4, cex = 1.5, col = "red")
#' plot(dom[,1], col = NA, add = TRUE)
#' 
#' # Draw a circle representing the basis' radius around one of the knot points
#' tseq = seq(0, 2*pi, length=100) 
#' rad = bs$get_w()
#' coords = cbind(rad * cos(tseq) + seq_x[2], rad * sin(tseq) + seq_y[2])
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
			private$basis_sp$get_dim()
		},
		get_cutpoints = function() {
			private$basis_sp$get_cutpoints()
		},
		get_w = function() {
			private$basis_sp$get_w()
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

				# Request a few more samples than we'll need, to prevent the loop in rdomain.
				P = rdomain(R, dom[j,], blocksize = ceiling(1.2*R), itmax = R)
				S[j,] = S[j,] + colSums(basis$compute(P[,1], P[,2]))
			}

			return(S / R)
		}
	)
)
