#' @export
prepare_bisquare = function(dom, knots_s, knots_t = NULL, type)
{
	if (type == "areal") {
		# here, dom should be an sf whose geometry provides shapes.
		# We won't check for the types of shapes provided.
		# We also leave X as NULL.
		X = NULL
		dom_crs = st_crs(dom)
	} else if (type == "point" && (any(c("sf","sfc") %in% class(dom)))) {
		# Here, dom should be an sf or sfc whose geometry provides POINTs.
		# Points can be 2-d (long,lat) or 3-d (long,lat,time)
		stopifnot(all(st_geometry_type(dom) == "POINT"))
		X = st_coordinates(dom)
		dom_crs = st_crs(dom)
	} else if (type == "point") {
		# Otherwise dom is intepreted as a matrix of points.
		X = as.matrix(dom)
		dom_crs = NULL
	}

	# Make sure X has two columns if knots_t is NULL; otherwise, X should
	# have three columns
	if (is.null(knots_t) && !type == "point") {
		stopifnot(ncol(X) == 2)
	} else if (type == "point") {
		stopifnot(ncol(X) == 3)
	}

	if (any(c("sf","sfc") %in% class(knots_s))) {
		if (!is.null(dom_crs)) {
			# Make sure projections match, if one was given for dom
			stopifnot(dom_crs == st_crs(knots_s))
		}
		stopifnot(all(st_geometry_type(knots_s) == "POINT"))
		knot_mat_s = st_coordinates(knots_s)
	} else {
		knot_mat_s = as.matrix(knots_s)
	}

	R = nrow(knot_mat_s)
	stopifnot(ncol(knot_mat_s) == 2)

	# If temporal knots are given, overall knots should be spatial knots crossed
	# with temporal knots.
	if (is.null(knots_t)) {
		knot_mat = knot_mat_s
	} else {
		T = length(knots_t)
		knot_mat = cbind(knot_mat_s %x% matrix(1,T,1), matrix(1,R,1) %x% knots_t)
	}

	list(X = X, knot_mat = knot_mat, R = R, T = T)
}

