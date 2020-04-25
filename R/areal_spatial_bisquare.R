#' Areal Spatial Bisquare Basis
#' 
#' An \code{R6Class} representing the spatial bisquare basis.
#' 
#' @section Usage:
#' \preformatted{
#' bs = ArealSpatialBisquareBasis$new(knots_x, knots_y, w)
#' bs$set_monte_carlo(reps)
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
#' \item \code{get_basis_sp} Get the underlying point-level basis function.
#' \item \code{get_dim} Get the number of cutpoints used to construct this basis.
#' \item \code{get_cutpoints} Get the cutpoints used to construct this basis;
#'   the matrix \code{cbind(knots_x,knots_y)}.
#' \item \code{get_w} Get the radius used to construct this basis.
#' \item \code{get_nx} Get the size of the Quadrature grid for the .
#' \item \code{get_ny} Get the size of the Quadrature grid for the y-axis.
#' \item \code{get_mc_reps} Get the number of Monte Carlo reps.
#' \item \code{get_report_period} Get the reporting frequency. A log message
#'   will be printed each time this many areas in a given domain have been
#'   processed.
#' \item \code{get_method} Get the computational method which is currently
#'   set: \code{MonteCarlo} or \code{Quadrature}.
#' \item \code{set_monte_carlo} Set the computational method to Monte Carlo
#'   and use the specified number of reps.
#' \item \code{set_quad} Set the computational method to Quadrature.
#'   A grid of \code{nx} by \code{ny} points will be used to evaluate each
#'   area in the domain.
#' \item \code{set_report_period} Set the reporting frequency.
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
		default_msg = NULL,
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
			private$default_msg = TRUE
			private$report_period = Inf
			private$basis_sp = SpatialBisquareBasis$new(knots_x, knots_y, w = w)
		},
		set_monte_carlo = function(reps) {
			private$method = "MonteCarlo"
			private$mc_reps = reps
			private$nx = NULL
			private$ny = NULL
			private$default_msg = FALSE
		},
		set_quad = function(nx, ny) {
			private$method = "Quadrature"
			private$mc_reps = NULL
			private$nx = nx
			private$ny = ny
			private$default_msg = FALSE
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
			if (private$default_msg) {
				# Print a message the first time compute is called, if the default MC
				# will be used.
				msg = sprintf("Using default method: %s with %d reps", private$method, private$mc_reps)
				message(msg)
				private$default_msg = FALSE
			}

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


#' @export
areal_spatial_bisquare = function(dom, knots, w, control = NULL)
{
	n = nrow(dom)
	r = nrow(knots)
	S = Matrix(0, n, r)

	if (is.null(control)) { control = list() }
	if (is.null(control$mc_reps)) { control$mc_reps = 1000 }
	if (is.null(control$nx)) { control$nx = 50 }
	if (is.null(control$ny)) { control$ny = 50 }
	if (is.null(control$report_period)) { control$report_period = n + 1 }
	if (is.null(control$verbose)) { control$verbose = FALSE }
	if (is.null(control$method)) { control$method = "MC" }
	if (is.null(control$mc_sampling_factor)) { control$mc_sampling_factor = 1.2 }

	R = control$mc_reps
	nx = control$nx
	ny = control$ny
	report_period = control$report_period
	verbose = control$verbose
	method = control$method
	blocksize = ceiling(control$mc_sampling_factor * R)

	if (verbose && method == "MC") {
		printf("Using Monte Carlo method with %d reps\n", R)
		printf("Sampling block size is %d\n", blocksize)
	} else if (verbose && method == "Quad") {
		printf("Using quadrature method with %d x %d grid\n", nx, ny)
	}

	if (verbose) {
		printf("Computing %d areas using %d knots\n", n, r)
		printf("Radius w = %g\n", w)
	}

	for (j in 1:n) {
		if (j %% report_period == 0 && verbose) {
			logger("Computing basis for area %d of %d\n", j, n)
		}

		if (method == "MC") {
			# The blocksize factor of helps to achieve the desired sample size
			# in rdomain without repeating the loop there.
			P = rdomain(R, dom[j,], blocksize = blocksize, itmax = R)
			B = spatial_bisquare(cbind(P$x,P$y), knots, w)
			S[j,] = colSums(B) / R
		} else if (method == "Quad") {
			grid_out = make_grid(dom[j,], nx, ny)
			B = spatial_bisquare(grid_out$X, knots, w)
			area = as.numeric(st_area(dom[j,]))
			S[j,] = colSums(B) * grid_out$dx * grid_out$dy / area
		} else {
			stop("Method not recongnized")
		}
	}

	return(S)
}
