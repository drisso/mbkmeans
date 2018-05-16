Goal: Add a function to your package that divides the data into 4 chunks of observations and calculates the mean and standard deviation of each variable for each chunk.

Two functions should be created:
1) A pure R function (i.e. only R code) (`chunkSummaryR` )
2) A R function that calls a C++ function via Rcpp. (`chunkSummaryC` )

Successful completion entails the following:
1) The package must pass R CMD CHECK 
2) The package must include a unit test that checks that the two functions give the same results on a small dataset
3) Appropriate commenting of the code as well creation of meaningful man pages for the function via Roxygen

