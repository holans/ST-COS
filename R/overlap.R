# Compute matrix of overlaps between two sets of areas (sf objects): D and G.
# If proportion == FALSE, the ith row of the result represents the amount of
# overlap, in area, ith area among each area of dom2. If proportion == TRUE,
# this is normalized to proportions which sum to 1.
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
		hh = rowSums(H)
		hh[hh == 0] = 1
		return( 1/hh * H )
	} else {
		return(H)
	}
}