#define ARMA_DONT_PRINT_ERRORS
# include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include "functions.h"
#include "utils_rcpp.h"
//#include "ClusterR/utils_rcpp.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <math.h>

using namespace arma;

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

template<typename T>
int get_ncol(const T& data){

  auto matrix_type=beachmat::find_sexp_type(data);

  if(matrix_type== INTSXP){
    auto final_matrix=beachmat::create_integer_matrix(data);
    const size_t& nc = final_matrix->get_ncol();
    //const size_t& nr = final_matrix->get_nrow();
    int n_row = nc;
    return n_row;

  }else if(matrix_type== REALSXP){
    auto final_matrix=beachmat::create_numeric_matrix(data);
    const size_t& nc = final_matrix->get_ncol();
    //const size_t& nr = final_matrix->get_nrow();
    int n_row = nc;
    return n_row;
  }else{

    return 0;
  }
}




template<typename T>
arma::rowvec clusters_WCSS(const T&data,arma::mat CENTROIDS){

  auto matrix_type=beachmat::find_sexp_type(data);

  if(matrix_type == INTSXP){

    auto final_matrix=beachmat::create_integer_matrix(data);
    int data_n_rows = get_nrow(data);
    int data_n_cols = get_ncol(data);
    Rcpp::IntegerMatrix dat_final(data_n_rows,data_n_cols);
    arma::rowvec CLUSTERS(data_n_rows);
    for (unsigned int j = 0; j < data_n_rows; j++) {

      Rcpp::NumericVector tmp(data_n_cols);
      final_matrix->get_row(j, tmp.begin());
      dat_final.row(j) = tmp;
      arma::mat data_final = Rcpp::as<arma::mat>(dat_final);
      arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(data_final.row(j)), CENTROIDS);                 // returns a rowvec with the SSE for each cluster

      //soft_CLUSTERS.row(j) = arma::conv_to< arma::rowvec >::from(tmp_vec);

      int tmp_idx = MinMat(tmp_vec);                                                                        // returns the index of the tmp_vec with the lowest SSE
      CLUSTERS(j) = tmp_idx+1;
    }
    return CLUSTERS;


  }else if(matrix_type ==REALSXP){
    auto final_matrix=beachmat::create_numeric_matrix(data);
    int data_n_rows = get_nrow(data);
    int data_n_cols = get_ncol(data);
    Rcpp::NumericMatrix dat_final(data_n_rows,data_n_cols);
    arma::rowvec CLUSTERS(data_n_rows);


    for (unsigned int j = 0; j < data_n_rows; j++) {

      Rcpp::NumericVector tmp(data_n_cols);
      final_matrix->get_row(j, tmp.begin());
      dat_final.row(j) = tmp;
      arma::mat data_final = Rcpp::as<arma::mat>(dat_final);
      arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(data_final.row(j)), CENTROIDS);                 // returns a rowvec with the SSE for each cluster

      //soft_CLUSTERS.row(j) = arma::conv_to< arma::rowvec >::from(tmp_vec);

      int tmp_idx = MinMat(tmp_vec);                                                                        // returns the index of the tmp_vec with the lowest SSE
      CLUSTERS(j) = tmp_idx+1;
    }
    return CLUSTERS;
  }
}







// [[Rcpp::export]]
arma::rowvec debug(SEXP data, Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue){

  arma::mat CENTROIDS1;

  if (CENTROIDS.isNotNull()) {

    CENTROIDS1 = Rcpp::as<arma::mat>(CENTROIDS);
  }

  int data_n_rows = get_nrow(data);

  int data_n_cols = get_ncol(data);

  Rcpp::NumericMatrix dat_final(data_n_rows,data_n_cols);

  arma::rowvec CLUSTERS(data_n_rows);

  arma::mat soft_CLUSTERS(data_n_rows, CENTROIDS1.n_rows);

  CLUSTERS=clusters_WCSS(data,CENTROIDS1);

  //for (unsigned int j = 0; j < data_n_rows; j++) {

    //Rcpp::NumericVector tmp(data_n_cols);

    //auto final_matrix=beachmat::create_numeric_matrix(data);

    //final_matrix->get_row(j, tmp.begin());

   // dat_final.row(j) = tmp;

    //arma::mat data_final = Rcpp::as<arma::mat>(dat_final);

    //arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(data_final.row(j)), CENTROIDS1);                 // returns a rowvec with the SSE for each cluster

    //soft_CLUSTERS.row(j) = arma::conv_to< arma::rowvec >::from(tmp_vec);

    //int tmp_idx = MinMat(tmp_vec);                                                                        // returns the index of the tmp_vec with the lowest SSE

    //CLUSTERS(j) = tmp_idx+1;
 // }

  return CLUSTERS;


}