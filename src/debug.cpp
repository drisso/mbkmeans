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


//get the number of rows
template<typename T>
int get_nrow(const T& data){

  auto matrix_type=beachmat::find_sexp_type(data);

  if(matrix_type== INTSXP){
    auto final_matrix=beachmat::create_integer_matrix(data);
    //const size_t& nc = final_matrix->get_ncol();
    const size_t& nr = final_matrix->get_nrow();
    int n_row = nr;
    return n_row;

  }else if(matrix_type== REALSXP){
    auto final_matrix=beachmat::create_numeric_matrix(data);
    //const size_t& nc = final_matrix->get_ncol();
    const size_t& nr = final_matrix->get_nrow();
    int n_row = nr;
    return n_row;
  }else{

    return 0;
  }
}



//get the number of columns
template<typename T>
int get_ncol(const T& data){

  auto matrix_type=beachmat::find_sexp_type(data);

  if(matrix_type== INTSXP){
    auto final_matrix=beachmat::create_integer_matrix(data);
    const size_t& nc = final_matrix->get_ncol();
    //const size_t& nr = final_matrix->get_nrow();
    int n_col = nc;
    return n_col;

  }else if(matrix_type== REALSXP){
    auto final_matrix=beachmat::create_numeric_matrix(data);
    const size_t& nc = final_matrix->get_ncol();
    //const size_t& nr = final_matrix->get_nrow();
    int n_col = nc;
    return n_col;
  }else{

    return 0;
  }
}





// [[Rcpp::export]]
arma::rowvec debug(Rcpp::IntegerMatrix data, Rcpp::NumericMatrix CENTROIDS) {

  ClustHeader clust_header;

  auto final_matrix=beachmat::create_integer_matrix(data);
  int data_n_rows = get_nrow(data);
  int data_n_cols = get_ncol(data);
  Rcpp::NumericMatrix dat_final(1,data_n_cols);
  Rcpp::IntegerVector tmp(data_n_cols);
  arma::rowvec CLUSTERS(data_n_rows);

  int centrod_n_row = get_nrow(CENTROIDS);

  arma::vec tmp_vec(centrod_n_row);

  arma::vec tmp_vec2(centrod_n_row);

  for (unsigned int j = 0; j < data_n_rows; j++) {

    final_matrix->get_row(j, tmp.begin());
    dat_final.row(0) = tmp;
    //arma::mat data_final = Rcpp::as<arma::mat>(dat_final);
    //arma::vec tmp_vec = clust_header.WCSS(arma::conv_to< arma::rowvec >::from(data_final.row(0)), CENTROIDS);                 // returns a rowvec with the SSE for each cluster

    for(unsigned int i =0; i<centrod_n_row;i++){

      tmp_vec(i)= Rcpp::sum(Rcpp::pow(dat_final.row(0) - CENTROIDS.row(i),2));

      tmp_vec2(i) = arma::as_scalar(tmp_vec(i));

    }
    //soft_CLUSTERS.row(j) = arma::conv_to< arma::rowvec >::from(tmp_vec);

    int tmp_idx = clust_header.MinMat(tmp_vec2);                                                                        // returns the index of the tmp_vec with the lowest SSE
    CLUSTERS(j) = tmp_idx+1;
  }
  return CLUSTERS;
}
