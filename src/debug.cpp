#define ARMA_DONT_PRINT_ERRORS
# include <RcppArmadillo.h>
# include <ClusterRHeader.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::depends(ClusterR)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include "functions.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <math.h>

using namespace arma;
using namespace clustR;


// [[Rcpp::export]]
NumericVector debug(NumericVector x) {
  return x * 2;
}
