#' Matrix of overlaps between two sets of areas.
#' 
#' A convenience function to convert output from \code{sf::st_intersection}
#' to a sparse matrix as defined in the \code{Matrix} package.
#' 
#' @param dom1 An \code{sf} object representing a domain of areal units.
#' @param dom2 An \code{sf} object representing a domain of areal units.
#' @param proportion Logical; if \code{TRUE}, normalize so that rows sum to 1.
#' Otherwise areas are returned.
#' @return An matrix of overlaps.
#' 
#' @details
#' Returns a matrix \code{H} whose (i,j)th entry represent the area of the overlap
#' between areal units \code{dom1[i,]} and \code{dom2[j,]}.
#' 
#' @examples
#' data("acs_sf")
#' dom1 = acs5_2013[1:10,]
#' dom2 = acs5_2016[1:10,]
#' H1 = overlap_matrix(dom1, dom2)
#' H2 = overlap_matrix(dom1, dom2, proportion = FALSE)
#'
#' @export
overlap_matrix = function(dom1, dom2, proportion = TRUE)
{
	D = dom1 %>%
		mutate(D_row_num = row_number()) %>%
		st_set_agr("constant")
	G = dom2 %>%
		mutate(G_row_num = row_number()) %>%
		st_set_agr("constant")
	INT = st_intersection(D, G)
	H = sparseMatrix(i = INT$D_row_num, j = INT$G_row_num,
		x = as.numeric(st_area(INT)), dims = c(nrow(D), nrow(G)))
	if (proportion) {
		hh = rowSums(H) + (rowSums(H) == 0)
		return(1/hh * H)
	} else {
		return(H)
	}
}