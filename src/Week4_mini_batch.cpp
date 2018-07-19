#define ARMA_DONT_PRINT_ERRORS
# include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

#ifdef _OPENMP
#include <omp.h>
#endif


#include "utils_rcpp.h"
#include "Week4_functions.h"


#include <algorithm>
#include <iostream>
#include <string>
#include <math.h>

//using Rcpp::Environment;

using namespace arma;
//using Rcpp::Environment;
//using namespace Rcpp;
//using namespace H5;

//' @export
// [[Rcpp::export]]

arma::mat Week4_mini_batch_kmeans(SEXP data){

  SEXP dat1= transfer_data(data);

  arma::mat final_data = Rcpp::as<arma::mat>(dat1);

  return final_data;
}