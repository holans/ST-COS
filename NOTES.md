# Notes
* Needed to install udunits library as a prereq for `sf` package. In Ubuntu, I
did this.
``` bash
S sudo apt-get install libudunits2-dev
```
* Needed to install an up-to-date GDAL library as a prereq for `sf` package.
In Ubuntu, I did this.
``` bash
$ sudo add-apt-repository ppa:ubuntugis/ppa
$ sudo apt-get update
$ sudo apt-get upgrade
$ sudo apt-get install libgdal-dev
```

* An explosive `M` might be a problem if we wanted to predict many steps ahead
("long-leap") into the future. But that isn't really a concern in this problem,
so we can leave `M` unconstrained. If we want to discuss this issue, see the
AOAS paper and possibly also the Cressie and Wilke book.

# Current issues
* Some of the shape files don't have their CRS coding set.
* My `Q` matrix seems to be different than Jon's (check the order of areal units).
* My `M` matrix is different than Jon's, even though the matrix `GBI` is the same.
* My sptcovar9 output is different than Jon's.
* `licols` in R and Matlab give very different results; QR decomposition estimates
  rank as 3500 in R and 1250 in Matlab. Maybe we should try a different reduction
  of `S`?

# To Code
- [ ] Can we make the MCMC faster? (Might need to be via linear algebra)
- [ ] Gelman-Rubin test for convergence
- [ ] Work out package design
- [x] Get Scott's feedback on using thinning.
- [x] Write a `compute.overlaps` function in the package using `sf` code.
- [x] Compute basis functions in the code.
- [x] Finish porting the `licols` function to R.
- [x] Add support to reduce the `S` matrix.
- [x] Add code to compute `Kinv`.
- [ ] Write a function to predict on a target geography after the MCMC is
      complete. We'll want predictions, `H %*% mu`, `S %*% eta`, etc. Also
      give direct access to draws, loglik, and DIC

# To Check
- [ ] Correctness of basis calculation vs Matlab
- [ ] Computation of `Kinv`
- [ ] Result of `sptcov` function

# To ask
- [ ] Should we keep the magic factor of 0.9 in `Q = I - 0.9 countAdj`?
- [ ] Was Matlab code switching `Kinv` and `K` incorrectly?
- [ ] 0.05 quantile of pairwise distances used for basis computation
- [ ] Interface to provide reduction of `S`
- [ ] The target VAR(1) process seems like it could become unstable if `M` isn't constrained?
- [ ] The propagator matrix `M` is currently allowed to be dense, but for the ACS example, this
	leads to a huge computational problem to get the autocovariance at lag 0. Could we assume
	some kind of sparsity structure so that not all entries of `Y_t^*` are autocorrelated? For
	example, maybe only adjacent areal units could be nonzero. That might be too strict, so
	we could loosen it to be based on degree in the adjacency graph. But then the norm
	minimization may not be so simple.

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
We might want a simple class hierarchy for the basis functions. the user will
construct one of these and pass it to the `STCOSPrep` object before adding
observed domains. `BisquareBasis` for example might be set up with cut points,
`w_s` and `w_t` parameters, a number of MC draws or gridpoints, etc. The
implementation of `BisquareBasis` could be in Rcpp, if it helps to speed things up.
Maybe the default choice for cutpoints should be called something like
`space_filling`, which takes a domain and a number of spatial points.

The time cutpoints would just be computed from the period. If `period` is
length 1, take cutpoint to be `mean(period-1, period)`. Otherwise take it
to be `mean(period)`. This at least matches Jon's numbers. Or even better,
we can let the user provide a function based on the period.

To predict a new target geography, we could have the user put together a new
`STCOSPrep` object. Except here they should not need to provide `Z` or `V`.
(Perhaps they would use the original object and call a different function
than `add_obs`?). This will allow them to extract the necessary `H` and `S`
matrices. They should then be able to get predictions with a call like
`predict(gibbs.out, H.new, S.new)`. Or they can get the raw MCMC draws
and compute prediction quantities themselves.

Should we somehow add the MCMC draws or summaries to a spatial data frame
for users?? If we wanted, we could copy the given spatial data frame `R.keep`
times, and put the `r`th set of predictions and SDs into the `r`th data frame.
But most users will probably just want to see the summary statistics. Users
might also want to see summaries of `H.new %*% mu_B` or `S.new %*% eta`, which
could also be mapped.

# Scott Meetings
* We should support/discuss reading shapefile and areal data separately, and joining them.
* We could auto-project all the observed domains to the fine-level one.
* We'll want to make some maps. Either by example or by making some helper functions.
* Make sure we give an error if projections aren't compatible
* We should indicate how long some of the longer steps will take