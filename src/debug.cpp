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

//using Rcpp::Environment;

using namespace arma;
//using Rcpp::Environment;
//using namespace Rcpp;
//using namespace H5;


//template<typename T1, typename T2>
//SEXP get_result(const T1& data, const T2& init_fraction){

//const size_t& nc = data->get_ncol();
//const size_t& nr = data->get_nrow();
//int fract = std::ceil(nr*init_fraction);

//arma::uvec init = arma::conv_to< arma::uvec >::from(sample_vec(fract, 0, nr - 1, false));

//arma::uvec samp_init = arma::sort(init);

//Rcpp::NumericMatrix submat(samp_init.n_rows, nc);
//Rcpp::NumericVector tmp(nc);

//for(int i=0; i<samp_init.n_rows; i++){
//  data->get_row(samp_init[i], tmp.begin());
//  submat.row(i) = tmp;
//}

//return submat;

//}




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


//shuffle the matrix and select randomly rows(nrow  = cluster)
template<typename T1>
SEXP shuffle_matrix_random(const T1& data, int cluster){
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






//' @export
// [[Rcpp::export]]
SEXP debug(SEXP data, double init_fraction, int clusters){

  arma::mat update_centroids;

  SEXP tran_data;

  auto matrix_type=beachmat::find_sexp_type(data);

  if(matrix_type == INTSXP){

    auto final_matrix=beachmat::create_integer_matrix(data);

    SEXP tran_data =  shuffle_matrix(final_matrix,init_fraction);

    //return tran_data;

  }else if(matrix_type ==REALSXP){

    auto final_matrix=beachmat::create_numeric_matrix(data);

    SEXP tran_data = shuffle_matrix(final_matrix,init_fraction);

    //return tran_data;

  }else{

    Rcpp::stop("The type of matrix is wrong");

  }

  //arma::mat final_data = Rcpp::as<arma::mat>(tran_data);

  //update_centroids = kmeans_pp_init(final_data, clusters, false);

  return tran_data;
}




//arma::mat update_centroids;

//SEXP trans_data = transfer_data(data);

//arma::mat final_data = Rcpp::as<arma::mat>(trans_data);

//update_centroids = kmeans_pp_init(final_data, clusters, false);

//return update_centroids;







//debug for shuffle_matrix
//auto matrix_type=beachmat::find_sexp_type(data);
//if(matrix_type==INTSXP){

  //auto final_matrix=beachmat::create_integer_matrix(data);
  //const size_t& nc = final_matrix->get_ncol();
  //const size_t& nr = final_matrix->get_nrow();
  //return get_result(final_matrix,init_fraction);
  //}else if(matrix_type==REALSXP){
  //auto final_matrix=beachmat::create_numeric_matrix(data);
  //const size_t& nc = final_matrix->get_ncol();
  //const size_t& nr = final_matrix->get_nrow();
  //return get_result(final_matrix,init_fraction);

  //}else{
  //return 0;
  //}









//debug for transfer data
  //auto dat = beachmat::create_numeric_matrix(data);

  //const size_t& nc=dat->get_ncol();
  //const size_t& nr=dat->get_nrow();


  //Rcpp::NumericMatrix final(nr,nc);

  //for(int i =0;i<nr;i++){
    //for(int j=0; j<nc;j++){
      //Rcpp::NumericVector tmp1(nr);

      //dat->get(i,j);

      //final(i,j)=dat->get(i,j);
    //}

  //}

  //return final;



  //Rcpp::RObject data1 = dat->yield();

  //SEXP fianl_data = data1.slot("seed");

  //SEXP final_data1 = data1.slot("/HDF5ArrayAUTO00001");

  //std::string final_data = get_class(final_data1);
  //Rcpp::Environment Matrix("package:Matrix"); // Load the Matrix package in R!
  //Rcpp::Function nearPD = Matrix["nearPD"];   // Extract nearPD() R function


  //Rcpp::List PD = nearPD(data);
  //Rcpp::S4 D_s4 = PD["mat"];

  //Rcpp::NumericVector temp = Rcpp::NumericVector(D_s4.slot("x"));
  //Rcpp::NumericVector dims = D_s4.slot("Dim");

  //arma::mat D(temp.begin(), // pointer to NumericVector
  //          dims[0],      // Number of Rows
  //              dims[1],      // Number of Columns
  //                  false,        // Avoid copying by disabling `copy_aux_mem`
  //                  true);        // Bind memory by enabling `strict`


  //SEXP final_matrix = as<NumericMatrix>(dat1);


