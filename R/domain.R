#' Draw uniformly distributed points from a set of areas
#' 
#' An alternative to \code{sf::st_sample} which draws uniformly distributed
#' points using a simple accept-reject method.
#' 
#' @param n Number of points desired in the final sample.
#' @param dom An \code{sf} object representing a domain of areal units.
#' @param blocksize Number of candidate points to draw on each pass of
#' accept-reject sampling (see details). Defaults to \code{n}.
#' @param itmax Maximum number of accept-reject samples to attempt. Defaults
#' to \code{Inf}.
#' @return An \code{sf} object with 2-dimensional points.
#' 
#' @details
#' Draws a sample of \code{blocksize} points uniformly from a bounding box on
#' \code{dom}, and accepts only the points which belong to \code{dom}. This
#' yields a uniform sample on \code{dom}. The process is repeated until \code{n}
#' accepted draws are obtained, or until it has been attempted \code{itmax}
#' times. If \code{itmax} iterations are reached without accepting \code{n}
#' draws, an error is thrown.
#' 
#' This seems to be an order of magnitude faster than the current
#' implementation of \code{st_sample}, although the latter can accomplish
#' the same objective and is more general. The improved performance is
#' worthwhile when used in the areal basis functions,
#' which sample repeatedly from the domain.
#' 
#' Performance will degrade when areal units have small area relative to their
#' bounding box, as many candidate points may need to be discarded. For
#' example, this will occur if \code{dom} contains a set of small scattered
#' islands in an ocean. In this case, it would be more efficient to sample
#' from each island at a time.
#' 
#' @examples
#' dom = acs5_2013[c(1,5,8,12),]
#' pts = rdomain(10000, dom)
#' 
#' @export
rdomain = function(n, dom, blocksize = n, itmax = Inf)
{
	bbox = st_bbox(dom)
	itr = 0
	n_accept = 0

	while (n_accept < n && itr < itmax) {
		itr = itr + 1
		u1 = runif(blocksize, bbox[1], bbox[3])
		u2 = runif(blocksize, bbox[2], bbox[4])
		u_df = data.frame(u1 = u1, u2 = u2)
		u_sf = st_as_sf(u_df, coords = c("u1", "u2"), crs = st_crs(dom), agr = "constant")
		idx = unique(unlist(st_contains(dom, u_sf)))

		if (n_accept == 0 && length(idx) > 0) {
			res = u_sf[idx,,drop=FALSE]
		} else {
			res = rbind(res, u_sf[idx,,drop=FALSE])
		}
		n_accept = n_accept + length(idx)
	}

	if (n_accept < n && itr == itmax) {
		stop("Reached maximum number of iterations without accepting n points")
	}

	return(res[1:n,])
}

# Note: the sf::st_make_grid function can similarly be used to make a grid within
# dom[j,], and with more flexibility built in. I.e.
#
# out = st_make_grid(dom[j,], what = "centers", n = c(nx,ny))
# X = matrix(unlist(out), ncol = 2, byrow = TRUE)
#
# As of this writing, the following code is much faster for our purposes.
make_grid = function(dom, nx, ny)
{
	bbox = st_bbox(dom)
	u1 = seq(bbox[1], bbox[3], length.out = nx+1)
	u2 = seq(bbox[2], bbox[4], length.out = ny+1)
	u_df = expand.grid(u1 = u1, u2 = u2)
	u_sf = st_as_sf(u_df, coords = c("u1","u2"), crs = st_crs(dom), agr = "constant")
	idx = unlist(st_contains(dom, u_sf))
	grid = u_sf[idx,]

	dx = (bbox[3] - bbox[1]) / nx
	dy = (bbox[4] - bbox[2]) / ny
	list(grid = grid, dx = dx, dy = dy)
}
