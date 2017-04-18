# Current issues
* Some of the shape files don't have their CRS coding set.
* My Q matrix seems to be different than Jon's (check the order of areal units).
* My M matrix is different than Jon's, even though the matrix GBI is the same.
* My sptcovar9 output is different than Jon's.
* licols in R and Matlab give very different results; QR decomposition estimates
  rank as 3500 in R and 1250 in Matlab. Maybe we should try a different reduction
  of S?

# TODO
- [ ] Can we make the MCMC faster? (Might need to be via linear algebra)
- [ ] Gelman-Rubin test for convergence
- [ ] Work out package design
- [ ] Get Scott's feedback on using thinning.
- [ ] Write a `compute.overlaps` function in the package using `sf` code.

# Preprocessing Steps
Compute the matrix `H` of overlaps.
1. User provides a set of shape files with their geographies (for observations
   and fine-scale) representing spatio-temporal domains.
2. We also need to make sure we have associated survey estimates `Z` and
   variance estimates `sig2var`. It might help to have them baked into the
   shape file, but I'm not sure yet.

Compute the `S` matrix representing the spatio-temporal basis functions.
1. Determine the spatio-temporal cutpoints. In Jon's code, there was a set
   of spatial cutpoints and a set of temporal cutpoints that were determined
   independently.
2. Apply the bisquare basis to obtain `S`.
3. Do some kind of reduction like `licols` or PCA to reduce the columns
   of `S`.

Compute the matrix `Kinv`.
1. Compute the adjacency matrix `Q`.
2. Compute the `Sconnector` matrices using another set of cutpoints.
3. Compute the Cov matrix of the target VAR process and combine with
  `Sconnector` matrices.
4. Apply Higham approximation.
