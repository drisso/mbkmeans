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

#include <algorithm> // std::move_backward
#include <random> // std::default_random_engine
#include <chrono> // std::chrono::system_clock


#include <vector>
#include <stdlib.h>
#include <time.h>
#include<list>

#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"

//using Rcpp::Environment;

using namespace arma;
//using Rcpp::Environment;
//using namespace Rcpp;
//using namespace H5;

//' @export
// [[Rcpp::export]]
SEXP random_choose(SEXP data, double init_fraction){

  auto mattype=beachmat::find_sexp_type(data);
  if (mattype==INTSXP) {
    auto dat1=beachmat::create_integer_matrix(data);

    const size_t& nr=dat1->get_nrow();
    const size_t& nc=dat1->get_ncol();

    int fract = std::ceil(nr * init_fraction);

    arma::uvec init = arma::conv_to< arma::uvec >::from(sample_vec(fract, 0, nr - 1, false));

    arma::uvec samp_init = arma::sort(init);

    Rcpp::IntegerMatrix submat(samp_init.n_rows, nc);
    Rcpp::IntegerVector tmp(nc);

    for(int i=0; i<samp_init.n_rows; i++){
      dat1->get_row(samp_init[i], tmp.begin());
      submat.row(i) = tmp;
    }

    return submat;

  } else if (mattype==REALSXP) {
    auto dat1=beachmat::create_numeric_matrix(data);

    const size_t& nr=dat1->get_nrow();
    const size_t& nc=dat1->get_ncol();

    int fract = std::ceil(nr * init_fraction);

    arma::uvec init = arma::conv_to< arma::uvec >::from(sample_vec(fract, 0, nr - 1, false));

    arma::uvec samp_init = arma::sort(init);

    Rcpp::NumericMatrix submat(samp_init.n_rows, nc);
    Rcpp::NumericVector tmp(nc);

    for(int i=0; i<samp_init.n_rows; i++){
      dat1->get_row(samp_init[i], tmp.begin());
      submat.row(i) = tmp;
    }

    return submat;

  } else {
    throw std::runtime_error("unacceptable matrix type");
  }

}
