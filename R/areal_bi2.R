# TBD: Keep drawing until we get B draws
# Also, might want to consider just picking points on a grid
function randomly_generate_Ds <- function(A, B)
{
	stopifnot(length(times) == 1)
	stopifnot(class(A) == "SpatialPolygonsDataFrame")

	lon.low <- A@bbox[1,1]
	lon.hi <- A@bbox[1,2]
	lat.low <- A@bbox[2,1]
	lat.hi <- A@bbox[2,2]

	u.lon <- runif(n, lon.low, lon.hi)
	u.lat <- runif(n, lat.low, lat.hi)
	u <- data.frame(Longitude = u.lon, Latitude = u.lat)
	coordinates(u) <- ~ Longitude + Latitude
	proj4string(u) <- proj4string(A)

	# Now reject the ones outside of the area
	idx <- which(!is.na(over(u, A)$GEO_ID))
	uu <- u@coords[idx,]

	# plot(A)
	# points(u[-idx], col = "black")
	# points(u[idx], col = "red")

	return(uu)
}

function ArealBi2 <- function(D, times, level1, level2, level3, B, srad, trad)
{
	r <- nrow(level1) + nrow(level2) + nrow(level3)
	n <- length(D)
	S <- matrix(0, n, r)
	T <- max(length(times))

	# s1{t} and s2{t} are n x B matrices

	for (t in 1:T) {
		for (j in 1:n) {
			points2 <- randomly_generate_Ds(D[j], B)
			s2{t}(j,:) <- points2(:,1)
			s1{t}(j,:) <- points2(:,2)
		}
	}

	for (t in 1:T) {
		for (i in 1:B) {
			if (i %% 1 == 0) { logger(" Iteration %d", i) }
			Sq <- Create_S([s2{t}(:,i), s1{t}(:,i)], times[t], level1, level2, level3, srad, trad)
			S <- Sq + S
		}
	}

	S <- S / (B*T)
}


function S <- Create_S(crd, times, level1, level2, level3, w, tw)
{
	level{1} = level1
	level{2} = level2
	level{3} = level3

	for i = 1:3
		if isempty(level{i}) == 0
			[B,IX] = sort(level{i},1)
			level{i} = level{i}(IX(:,2),:)
		else
			level{i} = []
		end
	end

	G = dist(level1)
	rl = [w*quantile(G(G>0),0.05), w*quantile(pdist(level2),0.05), w*quantile(pdist(level3),0.05)]
	S = matrix(0, nrow(crd), nrow(level1) + nrow(level2) + nrow(level3))

	count = 0
	for (i in 1:3) {
		if (!isempty(level{i}) == 0) {
			for (j in 1:nrow(crd)) {
				count = count + 1
				h = sqrt((level{i}(:,1) - crd(j,1))^2 + (level{i}(:,2) - crd(j,2))^2)
				ht = sqrt((level{i}(:,3) - times)^2)
				s = (1 - (h/rl(1,i))^2 - (ht/tw)^2)^2
				s[h > rl[1,i]] = 0
				if (i == 1) {
					S(j,1:length(level{i})) = s
				} else if (i == 2) {
					S(j,(length(level{1})+1):(length(level{1})+length(level{2}))) = s
				} else {
					S(j,(length(level{1})+length(level{2})+1):(length(level{1})+length(level{2})+length(level{3}))) = s
				}
			}
		}
	}
}
