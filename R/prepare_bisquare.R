prepare_bisquare = function(dom, knots, type)
{
	is_dom_sf = any(c("sf","sfc") %in% class(dom))
	if (type == "areal") {
		# here, dom should be an sf whose geometry provides shapes.
		# We won't check for the types of shapes provided.
		# We also leave X as NULL.
		X = NULL
	} else if (type == "point" && is_dom_sf) {
		# Here, dom should be an sf or sfc whose geometry provides POINTs.
		# Points can be 2-d (long,lat) or 3-d (long,lat,time)
		stopifnot(all(st_geometry_type(dom) == "POINT"))
		X = st_coordinates(dom)
	} else if (type == "point") {
		# Otherwise dom is intepreted as a matrix of points.
		X = as.matrix(dom)
	}

	# knots can be either 2d or 3d, and given as a matrix or sf/sfc
	is_knots_sf = any(c("sf","sfc") %in% class(knots))
	if (is_knots_sf) {
		stopifnot(all(st_geometry_type(knots) == "POINT"))
		knot_mat = st_coordinates(knots)
	} else {
		knot_mat = as.matrix(knots)
	}

	# Make sure coordinate systems match, if possible
	if (is_dom_sf && is_knots_sf) {
		stopifnot(st_crs(dom) == st_crs(knots))
	}

	# Make sure dimensions are appropriate, if possible
	if (type == "point") {
		stopifnot(ncol(X) == ncol(knot_mat))
	} else {
		stopifnot(ncol(knot_mat) %in% c(2,3))
	}

	list(X = X, knot_mat = knot_mat)
}
