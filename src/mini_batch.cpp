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
Rcpp::List mini_batch(SEXP data, int clusters, int batch_size, int max_iters, int num_init = 1, double init_fraction = 1.0, std::string initializer = "kmeans++",

                      int early_stop_iter = 10, bool verbose = false, Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, double tol = 1e-4, double tol_optimal_init = 0.5, int seed = 1){


  set_seed(seed);             // R's RNG

  int dat_n_rows = get_nrow(data);    // use the template function in this file

  //int dat_n_rows = data.n_rows;

  if (clusters > dat_n_rows - 2 || clusters < 1) { Rcpp::stop("the number of clusters should be at most equal to nrow(data) - 2 and never less than 1"); }

  bool flag = false;

  arma::mat CENTROIDS1;

  if (CENTROIDS.isNotNull()) {

    CENTROIDS1 = Rcpp::as<arma::mat>(CENTROIDS);

    num_init = 1;

    flag = true;
  }

  arma::mat centers_out;

  arma::rowvec iter_before_stop(num_init, arma::fill::zeros);

  int end_init = 0;

  arma::rowvec bst_WCSS(clusters);

  bst_WCSS.fill(arma::datum::inf);           // initialize WCSS to Inf, so that in first iteration it can be compared with the minimum of the 'WSSE'

  if (verbose) { Rcpp::Rcout << " " << std::endl; }

  arma::rowvec flag_exception(num_init, arma::fill::zeros);

  for (int init = 0; init < num_init; init++){

    arma::mat update_centroids;

    if(initializer== "kmeans++"){

      if(init_fraction==1.0){

        //SEXP centroids_matrix = transfer_data(data);

        //update_centroids = Rcpp::as<arma::mat>(centroids_matrix);

        SEXP trans_data = transfer_data(data);
        arma::mat update_centroids = Rcpp::as<arma::mat>(trans_data);


      }else if(init_fraction <1.0 && init_fraction >0.0){

        SEXP tran_data;
        auto matrix_type=beachmat::find_sexp_type(data);

        if(matrix_type == INTSXP){
          auto final_matrix=beachmat::create_integer_matrix(data);
          SEXP tran_data =  shuffle_matrix(final_matrix,init_fraction);
          return tran_data;
        }else if(matrix_type ==REALSXP){
          auto final_matrix=beachmat::create_numeric_matrix(data);
          SEXP tran_data = shuffle_matrix(final_matrix,init_fraction);
          return tran_data;
        }else{
          Rcpp::stop("The type of matrix is wrong");
        }

        arma::mat update_centroids = Rcpp::as<arma::mat>(tran_data);

       // SEXP centroids_matrix = shuffle_matrix(data,init_fraction = init_fraction);


      }else{
        Rcpp::stop("The value of fraction should be larger than 0 and not larger than 1");
      }
    }




  }

  return 0;

}


