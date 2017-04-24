# Current issues
* Some of the shape files don't have their CRS coding set.
* My `Q` matrix seems to be different than Jon's (check the order of areal units).
* My `M` matrix is different than Jon's, even though the matrix `GBI` is the same.
* My sptcovar9 output is different than Jon's.
* `licols` in R and Matlab give very different results; QR decomposition estimates
  rank as 3500 in R and 1250 in Matlab. Maybe we should try a different reduction
  of `S`?

# TODO
- [ ] Can we make the MCMC faster? (Might need to be via linear algebra)
- [ ] Gelman-Rubin test for convergence
- [ ] Work out package design
- [x] Get Scott's feedback on using thinning.
- [ ] Write a `compute.overlaps` function in the package using `sf` code.

# Preprocessing Steps
#### Compute the matrix `H` of overlaps.
1. User provides a set of shape files with their geographies (for observations
   and fine-scale) representing spatio-temporal domains.
2. We also need to make sure we have associated survey estimates `Z` and
   variance estimates `sig2var`. It might help to have them baked into the
   shape file, but I'm not sure yet.

#### Compute the `S` matrix representing the spatio-temporal basis functions.
1. Determine the spatio-temporal cutpoints. In Jon's code, there was a set
   of spatial cutpoints and a set of temporal cutpoints that were determined
   independently.
2. Apply the bisquare basis to obtain `S`.
3. Do some kind of reduction like `licols` or PCA to reduce the columns
   of `S`.

#### Compute the matrix `Kinv`.
1. Compute the adjacency matrix `Q`.
2. Compute the `Sconnector` matrices using another set of cutpoints.
3. Compute the Cov matrix of the target VAR process and combine with
  `Sconnector` matrices.
4. Apply Higham approximation.

# Some Design Notes

#### Inputs
* Shapefiles: for obs'd data and fine scale
* Survey ests and vars
* Time represented by each shape file
* "ACS" period(s) for each shape file

#### Options
* Basis
    - What kind of functions?
    - What should be the knots, cutpoints, weights, etc?
    - Number of reps (if computed with Monte Carlo)
    - How to reduce the dimension?
* Kinv
* MCMC
    - Draws, thin, burn, prior

#### Post-processing
* Access to draws
* loglik, DIC
* Posterior predictive draws

#### Some design questions
* Can we add some top-level "time" and "period" attributes to a shape object?
* Will we blow up the memory usage if we require the user to load all the geographies at once?
* Are all the geographies even needed at once? I think technically not:
    a. Load the fine-level geography
    b. For each obs'd geography, compute its portion of `H`, `S`, `Z` and
       `sig2var`. (It's like we want to build up a persistent state, using
       given choices for basis functions, cut points, etc. Then unload the
       shape file as soon as we can, or at least let the user do it. At the
       end, we should be able to grab all the stuff we need for the MCMC).
    c. 
    
    
```
fg = ... fine-level geography
sp = new STCOSPrep(fg)
for each observed geography
	read shapefile
	add survey estimates to shapefile, if not already there
	sp.AddObs(og, time, period)
	clean up to prevent memory waste
end

# At this point, we may want to reduce the number of basis functions
# For bisquare, we get tons of them because of the number of cutpoints.
# First, user gets the original S to investigate. After finding a good
# transformation, provide it as a function (because we need to do it
# later the exact same way). We should provide a few implementations
# like the licols method and PCA. Maybe the default should be the
# identity transformation?
S = sp$getS()
sp$set_S_reduction(f = function(S) { ... })

# User needs to provide a target covariance matrix. This could either go in
# the constructor, or the call to getCinv could fail if not set. We might be
# able to provide a reasonable default here (VAR(1) of the right size).

# Option 1: Call MCMC with sp
gibbs.out <- gibbs.stcos(sp, R, report.period = 10000, burn = 5000, thin = 10, hyper = NULL)
delete(sp)

# Option 2: Raw interface
# User shouldn't need to do this directly
# But they should be able to if they really want
H = sp.getH()
S = sp.getSreduced()
Z = sp.getZ()
sig2var = sp.getSig2var()
Cinv = sp.getCinv()
gibbs.out <- gibbs.stcos(Z, sig2var, H, S, Cinv, R, report.period = 10000, burn = 5000, thin = 10, hyper = NULL)
```

Okay, now AddObs will look something like this
```
function AddObs(og, time, period, Z_name, sig2var_name)
	m_areas[[n+1]]$H = compute.overlap(og, fg)
	m_areas[[n+1]]$S = compute.basis(og, time, period)
	m_areas[[n+1]]$Z = extract(og, Z_name)
	m_areas[[n+1]]$sig2var = extract(og, sig2var_name)
end
```

Then the accessors could put things together in the right order
```
function getH()
	HH = matrix(NA, N, n)
	for i = 1:n
		HH[idx] = m_areas[[n+1]]$H
	end
	return(HH)
end
```

What if we tried to do this without OOP in plain R? Then the user would have
to pass everything (like fg, all the basis function arguments, etc) each time
anything was called.


Let's see if we can use R's `R6` classes to implement this OOP idea...

```
library(R6)
STCOSPrep <- R6Class("STCOSPrep",
public = list(
	initialize = function(name = NA, hair = NA) {
		private$H_list <- list()
		private$S_list <- list()
		private$Z_list <- list()
		private$sig2var_list <- list()
	},
	add_obs = function(val) {
		...
	},
	get_H = function() {
		...
	}
),
private = list(
	H_list = NULL,
	S_list = NULL,
	Z_list = NULL,
	sig2var_list = NULL,
),
lock_class = TRUE
)
```

TBD: Where to pass all the stuff we need to compute Cinv?