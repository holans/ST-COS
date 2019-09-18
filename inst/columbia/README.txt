The following files are included in this directory:
- 00-data-prep-aff.R: Assemble analysis datasets which have been downloaded
  from the American FactFinder website. 
- 00-data-prep-api.R: Assemble analysis datasets from the Census Data API.
  (The API is a service provided by the U.S. Census Bureau; therefore, call
  format and data availability are subject to change).
- 01-assemble-data.R: Demonstrates how to assemble inputs for STCOS analysis
  using ACS estimates and shapefiles that can be retrieved publicly on the
  internet.
- 02-prepare-analysis.R: Prepare quantities used in the STCOS model from
  assembled data.
- 03-fit-gibbs.R: Run the Gibbs sampler from the stcos package.
- 04-results.R: Gather results from the Gibbs sampler and prepare outputs.
- 05-fit-stan.R: As an alternative to the stcos Gibbs sampler, use Stan to
  produce MCMC draws from the STCOS model.
- neighborhoods.{dbf,prj,shp,shx}: Shapefiles representing the four City
  of columbia neighborhoods which we take to be target supports.
- stcos.stan: Stan code used by 05-fit-stan.R

Running the R scripts numbered 01-05 above, in the given sequence, should
produce similar results to the article <arXiv:1904.12092>.

The following R packages are needed as prerequisites to run the code in this
directory:
- coda
- dplyr
- fields
- ggforce
- ggplot2
- ggrepel
- gridExtra
- jsonlite
- rstan
- stcos
- sf
- tigris

