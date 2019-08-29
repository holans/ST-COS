#' Sparse adjacency matrix between two sets of areas.
#' 
#' A convenience function to convert output from \code{sf::st_touches}
#' to a sparse matrix as defined in the \code{Matrix} package.
#' 
#' @param dom An \code{sf} object representing a domain of areal units.
#' @return An adjacency matrix
#' 
#' @details 
#' Returns a matrix \code{A} whose (i,j)th entry contains a 1 if
#' areal units \code{dom[i,]} and \code{dom[j,]} are adjacent;
#' 0 otherwise.
#' 
#' @examples
#' data("acs_sf")
#' dom = acs5_2013[1:4,]
#' A = adjacency_matrix(dom)
#' 
#' @export
adjacency_matrix = function(dom)
{
	adj = st_touches(dom, dom, sparse = TRUE)
	L = length(adj)
	counts = unlist(Map(length, adj))
	n_tuples = sum(counts)
	i = rep(1:L, counts)
	j = unlist(adj)
	sparseMatrix(i = i, j = j, x = 1, dims = c(L, L))
}
