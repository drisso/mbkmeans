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

template<typename T>
arma::rowvec clusters_WCSS(const T&data,arma::mat CENTROIDS){

  ClustHeader clust_header;

  auto matrix_type=beachmat::find_sexp_type(data);

  if(matrix_type == INTSXP){

    auto final_matrix=beachmat::create_integer_matrix(data);
    int data_n_rows = get_nrow(data);
    int data_n_cols = get_ncol(data);
    Rcpp::IntegerMatrix dat_final(data_n_rows,data_n_cols);
    arma::rowvec CLUSTERS(data_n_rows);
    for (unsigned int j = 0; j < data_n_rows; j++) {

      Rcpp::IntegerVector tmp(data_n_cols);
      final_matrix->get_row(j, tmp.begin());
      dat_final.row(j) = tmp;
      arma::mat data_final = Rcpp::as<arma::mat>(dat_final);
      arma::vec tmp_vec = clust_header.WCSS(arma::conv_to< arma::rowvec >::from(data_final.row(j)), CENTROIDS);                 // returns a rowvec with the SSE for each cluster

      //soft_CLUSTERS.row(j) = arma::conv_to< arma::rowvec >::from(tmp_vec);

      int tmp_idx = clust_header.MinMat(tmp_vec);                                                                        // returns the index of the tmp_vec with the lowest SSE
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
      arma::vec tmp_vec = clust_header.WCSS(arma::conv_to< arma::rowvec >::from(data_final.row(j)), CENTROIDS);                 // returns a rowvec with the SSE for each cluster

      //soft_CLUSTERS.row(j) = arma::conv_to< arma::rowvec >::from(tmp_vec);

      int tmp_idx = clust_header.MinMat(tmp_vec);                                                                        // returns the index of the tmp_vec with the lowest SSE
      CLUSTERS(j) = tmp_idx+1;
    }
    return CLUSTERS;
  }

}



//' Predict_mini_batch
//'
//' Prediction function for Mini-batch-k-means for in-memory, delayed, and on-disk matrices
//'
//'
//'@param data matrix, DelayedMatrix, or HDF5Matrix containing numeric or
//'  integer data (obseravtions in rows, variables in columns)
//'@param CENTROIDS a matrix of initial cluster centroids. The rows of the
//'  CENTROIDS matrix should be equal to the number of clusters and the columns
//'  should equal the columns of the data.
//'@return it returns a vector with the clusters.
//'@details
//'
//'This function takes the data and the output centroids and returns the
//'clusters.
//'
//'This implementation relies very heavily on the
//'\code{\link[ClusterR]{MiniBatchKmeans}} implementation. We provide the
//'ability to work with DelayedMatrix and HDF5Matrix through the \code{beachmat}
//'library.
//'
//'@author Lampros Mouselimis and Yuwei Ni
//'
//'@examples
//'data(iris)
//'km = mini_batch(as.matrix(iris[,1:4]), clusters = 3,
//'                batch_size = 10, max_iters = 10)
//'clusters = predict_mini_batch(as.matrix(iris[,1:4]),
//'                              CENTROIDS = km$centroids)
//' @export
// [[Rcpp::export]]
arma::rowvec predict_mini_batch(SEXP data, Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, bool fuzzy = false, double eps = 1.0e-6) {

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
    return CLUSTERS;
  }
}

