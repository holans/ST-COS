#' Areal Spatial Bisquare Basis
#' 
#' @description
#' Spatial bisquare basis on areal data.
#' 
#' @param dom An \code{sf} or \code{sfc} object with areas
#' \eqn{A_1, \ldots, A_n} to evaluate.
#' @param knots Knots \eqn{\bm{c}_1, \ldots, \bm{c}_r} for the basis.
#' See "Details".
#' @param w Radius for the basis.
#' @param control A \code{list} of control arguments. See "Details".
#'
#' @return A sparse \eqn{n \times r} matrix whose \eqn{i}th row
#' is
#' \eqn{
#' \bm{s}_i^\top =
#' \Big(
#' \bar{\varphi}_1(A_i), \ldots, \bar{\varphi}_r(A_i)
#' \Big).
#' }
#' 
#' @details
#' Notes about arguments:
#' \itemize{
#' \item \code{knots} may be provided as either an \code{sf} or \code{sfc} object, or as a
#'   matrix of points.
#' \item If an \code{sf} or \code{sfc} object is provided for \code{knots}, \eqn{r}
#'   two-dimensional \code{POINT} entries are expected in \code{st_geometry(knots)}.
#'   Otherwise, \code{knots} will be interpreted as an \eqn{r \times 2} numeric matrix.
#' \item If \code{knots} is an \code{sf} or \code{sfc} object, it is checked
#'   to ensure the coordinate system matches \code{dom}.
#' }
#' 
#' For each area \eqn{A} in the given domain, compute an the basis functions
#' \deqn{
#' \bar{\varphi}_j(A) = \frac{1}{|A|} \int_A \varphi_j(\bm{u}) d\bm{u}
#' }
#' for \eqn{j = 1, \ldots, r}. Here, \eqn{\varphi_j(\bm{u})} represent
#' \link{spatial_bisquare} basis functions defined at the point level
#' using \eqn{\bm{c}_j} and \eqn{w}.
#' 
#' The basis requires an integration which may be computed using one
#' of two methods. The \code{mc} method uses
#' \deqn{
#' \bar{\varphi}_j(A) \approx
#' \frac{1}{Q} \sum_{q=1}^Q \varphi_j(\bm{u}_q),
#' }
#' based on a random sample of locations \eqn{\bm{u}_1, \ldots, \bm{u}_Q} from
#' a uniform distribution on area \eqn{A}. The \code{rect} method uses
#' a simple quadrature approximation
#' \deqn{
#' \bar{\varphi}_j(A) \approx
#' \frac{1}{|A|}  \sum_{a=1}^{n_x} \sum_{b=1}^{n_y} \varphi_j(\bm{u}_{ab})
#' I(\bm{u}_{ab} \in A) \Delta_x \Delta_y.
#' }
#' Here, the bounding box \code{st_bbox(A)} is divided evenly into a grid of
#' \eqn{n_x \times n_y} rectangles, each of size \eqn{\Delta_x \times \Delta_y}.
#' Each \eqn{\bm{u}_{ab} = (u_a, u_b)} is a point from the \eqn{(a,b)}th
#' rectangle, for \eqn{a = 1, \ldots, n_x} and \eqn{b = 1, \ldots, n_y}.
#' 
#' Due to the treatment of \eqn{A_i} and \eqn{\bm{c}_j} as objects in a
#' Euclidean space, this basis is more suitable for coordinates from a map
#' projection than coordinates based on a globe representation.
#' 
#' The \code{control} argument is a list which may provide any of the following:
#' \itemize{
#' \item \code{method} specifies computation method: \code{mc} or \code{rect}.
#' Default is \code{mc}.
#' \item \code{mc_reps} is number of repetitions to use for \code{mc}.
#' Default is 1000.
#' \item \code{nx} is number of x-axis points to use for \code{rect}
#' method. Default is 50.
#' \item \code{ny} is number of y-axis oints to use for \code{rect}
#' method. Default is 50.
#' \item \code{report_period} is an integer; print a message with progress each
#' time this many areas are processed. Default is \code{Inf} so that message
#' is suppressed.
#' \item \code{verbose} is a logical; if \code{TRUE} print descriptive
#' messages about the computation. Default is \code{FALSE}.
#' \item \code{mc_sampling_factor} is a positive number; an oversampling factor
#' used to compute \code{blocksize} in the \link{rdomain} function. I.e.,
#' \code{blocksize = ceiling(mc_sampling_factor * mc_reps)}. Default
#' is 1.2.
#' }
#' 
#' @examples
#' set.seed(1234)
#' 
#' # Create knot points
#' seq_x = seq(0, 1, length.out = 3)
#' seq_y = seq(0, 1, length.out = 3)
#' knots = expand.grid(x = seq_x, y = seq_y)
#' knots_sf = st_as_sf(knots, coords = c("x","y"), crs = NA, agr = "constant")
#' 
#' # Create a simple domain (of rectangles) to evaluate
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
#' rad = 0.5
#' areal_spatial_bisquare(dom, knots, rad)
#' areal_spatial_bisquare(dom, knots_sf, rad)
#' 
#' # Plot the knots and the points at which we evaluated the basis
#' plot(knots[,1], knots[,2], pch = 4, cex = 1.5, col = "red")
#' plot(dom[,1], col = NA, add = TRUE)
#' 
#' # Draw a circle representing the basis' radius around one of the knot points
#' tseq = seq(0, 2*pi, length=100) 
#' coords = cbind(rad * cos(tseq) + seq_x[2], rad * sin(tseq) + seq_y[2])
#' lines(coords, col = "red")
#' 
#' @family bisquare
#' @export
areal_spatial_bisquare = function(dom, knots, w, control = NULL)
{
	prep = prepare_bisquare(dom, knots, type = "areal")
	X = prep$X
	knot_mat = prep$knot_mat

	r = nrow(knot_mat)
	n = nrow(dom)
	S = Matrix(0, n, r)

	if (is.null(control)) { control = list() }
	if (is.null(control$mc_reps)) { control$mc_reps = 1000 }
	if (is.null(control$nx)) { control$nx = 50 }
	if (is.null(control$ny)) { control$ny = 50 }
	if (is.null(control$report_period)) { control$report_period = n + 1 }
	if (is.null(control$verbose)) { control$verbose = FALSE }
	if (is.null(control$method)) { control$method = "mc" }
	if (is.null(control$mc_sampling_factor)) { control$mc_sampling_factor = 1.2 }

	reps = control$mc_reps
	nx = control$nx
	ny = control$ny
	report_period = control$report_period
	verbose = control$verbose
	method = control$method
	blocksize = ceiling(control$mc_sampling_factor * reps)

	if (verbose && method == "mc") {
		printf("Using `mc` method with %d reps\n", reps)
		printf("Sampling block size is %d\n", blocksize)
	} else if (verbose && method == "quad") {
		printf("Using 'rect' method with %d x %d grid\n", nx, ny)
	}

	if (verbose) {
		printf("Computing %d areas using %d spatial knots\n", n, r)
		printf("Radius w = %g\n", w)
	}

	dom_area = st_area(dom)

	for (j in 1:n) {
		if (j %% report_period == 0) {
			logger("Computing basis for area %d of %d\n", j, n)
		}

		if (method == "mc") {
			# The blocksize factor of helps to achieve the desired sample size
			# in rdomain without repeating the loop there.
			P = rdomain(reps, dom[j,], blocksize = blocksize, itmax = reps)
			X = st_coordinates(P)
			out = compute_basis_sp(X, knot_mat, w)
			B = sparseMatrix(i = out$ind_row + 1, j = out$ind_col + 1, x = out$values,
				dims = out$dim)
			S[j,] = colSums(B) / reps
		} else if (method == "rect") {
			grid_out = make_grid(dom[j,], nx, ny)
			X = st_coordinates(grid_out$grid)
			out = compute_basis_sp(X, knot_mat, w)
			B = sparseMatrix(i = out$ind_row + 1, j = out$ind_col + 1, x = out$values,
				dims = out$dim)
			S[j,] = colSums(B) * grid_out$dx * grid_out$dy / dom_area[j]
		} else {
			stop("Unrecongnized method. Use `mc` or 'rect'")
		}
	}

	return(S)
}
