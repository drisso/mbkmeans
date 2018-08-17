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


//subset the matrix and select randomly rows(nrows = init_fraction*total)
template<typename T1, typename T2>
SEXP subset_matrix(const T1& data, const T2& init_fraction_value){

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


//subset the matrix and select randomly rows(nrow  = cluster)
template<typename T1>
SEXP subset_matrix_random(const T1& data, int cluster){
  const size_t& nc = data->get_ncol();
  const size_t& nr = data->get_nrow();

  arma::uvec init = arma::conv_to< arma::uvec >::from(sample_vec(cluster, 0, nr - 1, false));

  arma::uvec samp_init = arma::sort(init);

  Rcpp::NumericMatrix submat(samp_init.n_rows, nc);
  Rcpp::NumericVector tmp(nc);

  for(int i=0; i<samp_init.n_rows; i++){
    data->get_row(samp_init[i], tmp.begin());
    submat.row(i) = tmp;
  }

  return submat;

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



//' Predict_mini_batch
//'
//' Prediction function for Mini-batch-k-means for both matrix and HDF5Matrix
//'
//'
//'@param data numeric matrix or integer matrix or HDF5Matrix
//'@param CENTROIDS a matrix of initial cluster centroids. The rows of the CENTROIDS matrix should be equal to the number of clusters and the columns should equal the columns of the data.
//'@return it returns a vector with the clusters.
//'@details
//'
//'This function takes the data and the output centroids and returns the clusters.
//'
//'@references
//'https://github.com/mlampros/ClusterR
//'
//' @export
// [[Rcpp::export]]
Rcpp::List predict_mini_batch(SEXP data, Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, bool fuzzy = false, double eps = 1.0e-6) {

  arma::mat CENTROIDS1;

  if (CENTROIDS.isNotNull()) {

    CENTROIDS1 = Rcpp::as<arma::mat>(CENTROIDS);
  }

  //SEXP trans_data = transfer_data(data);

 // arma::mat dat_final = Rcpp::as<arma::mat>(trans_data);

  int data_n_rows = get_nrow(data);

  int data_n_cols = get_ncol(data);

  arma::rowvec CLUSTERS(data_n_rows);

  arma::mat soft_CLUSTERS(data_n_rows, CENTROIDS1.n_rows);

  CLUSTERS =clusters_WCSS(data,CENTROIDS1);




  //for (unsigned int j = 0; j < data_n_rows; j++) {

  //Rcpp::NumericVector tmp(data_n_cols);

  //auto final_matrix=beachmat::create_numeric_matrix(data);

  //final_matrix->get_row(j, tmp.begin());

  //dat_final.row(j) = tmp;

  //arma::mat data_final = Rcpp::as<arma::mat>(dat_final);

  //arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(data_final.row(j)), CENTROIDS1);

    //arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(dat_final.row(j)), CENTROIDS1);                 // returns a rowvec with the SSE for each cluster

    //soft_CLUSTERS.row(j) = arma::conv_to< arma::rowvec >::from(tmp_vec);

    //int tmp_idx = MinMat(tmp_vec);                                                                        // returns the index of the tmp_vec with the lowest SSE

    //CLUSTERS(j) = tmp_idx+1;
    // }



  if (fuzzy) {
//
//     arma::mat fuzzy_mat(soft_CLUSTERS.n_rows, soft_CLUSTERS.n_cols);
//
//     for (unsigned int i = 0; i < soft_CLUSTERS.n_rows; i++) {
//
//       fuzzy_mat.row(i) = norm_fuzzy(arma::conv_to< arma::rowvec >::from(soft_CLUSTERS.row(i)), eps);
//     }
//
//     return Rcpp::List::create(Rcpp::Named("clusters") = CLUSTERS, Rcpp::Named("fuzzy_clusters") = fuzzy_mat);

  Rcpp::stop("fuzzy clustering is currently not implemented.");
  }


  else {

    return Rcpp::List::create(Rcpp::Named("clusters") = CLUSTERS);
  }
}

