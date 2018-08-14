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
#include "functions.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <math.h>

using namespace arma;




//get the number of rows
template<typename T>
int get_nrow(const T& data){
  auto matrix_type=beachmat::find_sexp_type(data);
  if(matrix_type== INTSXP){
    auto final_matrix=beachmat::create_integer_matrix(data);
    //const size_t& nc = final_matrix->get_ncol();
    const size_t& nr = final_matrix->get_nrow();
    int n_col = nr;
    return n_col;

  }else if(matrix_type== REALSXP){
    auto final_matrix=beachmat::create_numeric_matrix(data);
    //const size_t& nc = final_matrix->get_ncol();
    const size_t& nr = final_matrix->get_nrow();
    int n_col = nr;
    return n_col;
  }else{

    return 0;
  }
}


template<typename T1>
arma::vec pred_WCSS1(const T1& data, Rcpp::NumericMatrix CENTROIDS){

    const size_t& nc = data->get_ncol();
    const size_t& nr = data->get_nrow();

    //Rcpp::NumericMatrix submat(samp_init.n_rows, nc);
    Rcpp::NumericMatrix WCSS_matrix;

    Rcpp::NumericVector tmp_vec(nc);

    Rcpp::NumericVector tmp(nc);

    arma::vec final_vec;

    arma::mat CENTROIDS1 = Rcpp::as<arma::mat>(CENTROIDS);

    for(int i=0; i<nr; i++){

      data->get_row(i, tmp.begin());

      //WCSS_matrix[i] = WCSS(tmp[i],CENTROIDS1);

      data->get_row(WCSS_matrix[i],tmp_vec.begin());

      final_vec = tmp_vec;
    }

    return final_vec;
}



//' @export
// [[Rcpp::export]]
arma::mat debug(SEXP data, Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, bool fuzzy = false, double eps = 1.0e-6) {

  arma::mat CENTROIDS1;

  if (CENTROIDS.isNotNull()) {

    CENTROIDS1 = Rcpp::as<arma::mat>(CENTROIDS);
  }

  int data_n_rows = get_nrow(data);

  arma::rowvec CLUSTERS(data_n_rows);

  arma::mat soft_CLUSTERS(data_n_rows, CENTROIDS1.n_rows);

  //arma::vec final = pred_WCSS1(data,CENTROIDS1);

  return 0;

}


