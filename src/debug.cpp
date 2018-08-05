#define ARMA_DONT_PRINT_ERRORS
# include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include "utils_rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include "Week4_functions.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <math.h>

using namespace arma;



//shuffle the matrix and select randomly rows(nrows = init_fraction*total)
template<typename T1, typename T2>
SEXP shuffle_matrix(const T1& data, const T2& init_fraction_value){

  const size_t& nc = data->get_ncol();
  const size_t& nr = data->get_nrow();
  int fract = std::ceil(nr*init_fraction_value);

  arma::uvec init = arma::conv_to< arma::uvec >::from(sample_vec(fract, 0, nr - 1, false));

  arma::uvec samp_init = arma::sort(init);

  Rcpp::NumericMatrix submat(samp_init.n_rows, nc);
  Rcpp::NumericVector tmp(nc);

  for(int i=0; i<samp_init.n_rows; i++){
    data->get_row(samp_init[i], tmp.begin());
    submat.row(i) = tmp;
  }

  return submat;

}


//' @export
// [[Rcpp::export]]
arma::mat debug(SEXP data, double init_fraction,int clusters){

  arma::mat update_centroids;

  Rcpp::NumericMatrix tran_data;

  auto matrix_type=beachmat::find_sexp_type(data);

  if(matrix_type == INTSXP){

    auto final_matrix=beachmat::create_integer_matrix(data);

    tran_data =shuffle_matrix(final_matrix,init_fraction);

    //return tran_data;

  }else if(matrix_type ==REALSXP){

    auto final_matrix=beachmat::create_numeric_matrix(data);

    SEXP tran_data = shuffle_matrix(final_matrix,init_fraction);

    //return tran_data;

  }else{

    //Rcpp::stop("The type of matrix is wrong");

    SEXP tran_data = data;

    //return data;
  }

  arma::mat final_data = Rcpp::as<arma::mat>(tran_data);

  update_centroids = kmeans_pp_init(final_data, clusters, false);

  //return tran_data;

  //return 3;



  return update_centroids;
}


