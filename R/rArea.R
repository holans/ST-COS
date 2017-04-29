rArea <- function(n, area, blocksize = n)
{
	N <- nrow(area)
	res <- matrix(NA, 0, 2)

	# bbox <- st_bbox(area)
	# x.min <- bbox[1]
	# y.min <- bbox[2]
	# x.max <- bbox[3]
	# y.max <- bbox[4]

	done <- FALSE
	while(!done) {
		# u.x <- runif(blocksize, x.min, x.max)
		# u.y <- runif(blocksize, y.min, y.max)
		# u <- st_multipoint(cbind(u.x,u.y))
		# inc <- numeric(blocksize)
		# for (i in 1:blocksize) {
		#	inc[i] <- (length(st_contains(area, st_point(u[i,]))[[1]]) > 0)
		# }

		u <- st_sample(area, size = blocksize)

		# plot(area[,1])
		# plot(u, col = "red", add = TRUE)
		# length(u)

		M <- matrix(unlist(u), length(u), 2, byrow = TRUE)
		res <- rbind(res, M)
		
		done <- (nrow(res) >= n)
	}

	colnames(res) <- c("x", "y")
	res[1:n,]
}
