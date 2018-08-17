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



template<typename T>
arma::vec distance_wcss(const T&data,arma::mat centroid){

  const size_t& nc = data->get_ncol();
  const size_t& nr = data->get_nrow();

  arma::rowvec CLUSTERS(nr);

  arma::vec tmp_c(centroid.n_rows);

  arma::vec tmp_data(nc);

  Rcpp::NumericMatrix submat(nr,nc);

  arma::vec tmp_vec;

 // arma::rowvec CLUSTERS(nr);

  arma::mat soft_CLUSTERS(nr, centroid.n_rows);


  for(unsigned int i=0;i<nr;i++){


      data->get_row(i,tmp_data.begin());

      //arma::vec tmp_ss = WCSS(tmp_data,centroid);

       tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(tmp_data), centroid);

      //soft_CLUSTERS.row(i) = arma::conv_to< arma::rowvec >::from(tmp_vec);

      int tmp_idx = MinMat(tmp_vec);

      CLUSTERS(i) = tmp_idx;







  }

  return CLUSTERS;

}


// [[Rcpp::export]]
arma::vec debug(Rcpp::IntegerMatrix data, Rcpp::NumericMatrix CENTROIDS){

  arma::mat CENTROIDS1 = Rcpp::as<arma::mat>(CENTROIDS);

  int data_n_rows = get_nrow(data);

  arma::vec tmp_vec_wcss;

  auto final_matrix=beachmat::create_integer_matrix(data);

  tmp_vec_wcss =distance_wcss(final_matrix,CENTROIDS1);


  return tmp_vec_wcss;


}