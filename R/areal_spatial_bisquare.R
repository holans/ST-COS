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
	lock_class = FALSE,
	private = list(
		method = NULL,
		mc_reps = NULL,
		nx = NULL,
		ny = NULL,
		report_period = NULL,
		basis_sp = NULL
	),
	public = list(
		initialize = function(knots_x, knots_y, w) {
			# Default is Monte Carlo with 1000 reps
			self$set_monte_carlo(1000)
			private$report_period = Inf
			private$basis_sp = SpatialBisquareBasis$new(knots_x, knots_y, w = w)
		},
		set_monte_carlo = function(reps) {
			private$method = "MonteCarlo"
			private$mc_reps = reps
			private$nx = NULL
			private$ny = NULL
		},
		set_quad = function(nx, ny) {
			private$method = "Quadrature"
			private$mc_reps = NULL
			private$nx = nx
			private$ny = ny
		},
		set_w = function(w) {
			knots = private$basis_sp$get_cutpoints()
			private$basis_sp = SpatialBisquareBasis$new(knots[,1], knots[,2], w = w)
		},
		set_report_period = function(report_period) {
			private$report_period = report_period
		},
		get_method = function() {
			private$method
		},
		get_mc_reps = function() {
			private$mc_reps
		},
		get_nx = function() {
			private$nx
		},
		get_ny = function() {
			private$nx
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
		get_report_period = function() {
			private$report_period
		},
		compute = function(dom) {
			if (private$method == "MonteCarlo") {
				S = private$compute_mc(dom)
			} else if (private$method == "Quadrature") {
				S = private$compute_quad(dom)
			} else {
				stop("Unknown method")
			}
			return(S)
		}
	)
)

ArealSpatialBisquareBasis$set("private", "compute_mc",
	function(dom) {
		bs = private$basis_sp
		n = nrow(dom)
		r = bs$get_dim()
		S = Matrix(0, n, r)
		R = private$mc_reps
		report_period = private$report_period

		for (j in 1:n) {
			if (j %% report_period == 0) {
				logger("Computing basis for area %d of %d\n", j, n)
			}

			# The factor of 1.2 helps to achieve the desired sample size
			# in rdomain without repeating the loop there.
			P = rdomain(R, dom[j,], blocksize = ceiling(1.2*R), itmax = R)
			B = bs$compute(P$x, P$y)
			S[j,] = S[j,] + colSums(B)
		}

		return(S / R)
	}
)

ArealSpatialBisquareBasis$set("private", "compute_quad",
	function(dom) {
		bs = private$basis_sp
		n = nrow(dom)
		r = bs$get_dim()
		S = Matrix(0, n, r)
		nx = private$nx
		ny = private$ny
		report_period = private$report_period

		for (j in 1:n) {
			if (j %% report_period == 0) {
				logger("Computing basis for area %d of %d\n", j, n)
			}

			out = make_grid(dom[j,], nx, ny)
			X = out$X
			dx = out$dx
			dy = out$dy

			B = bs$compute(X[,1], X[,2])
			area = as.numeric(st_area(dom[j,]))
			S[j,] = colSums(B) * dx * dy / area
		}
		return(S)			
	}
)

ArealSpatialBisquareBasis$lock()
