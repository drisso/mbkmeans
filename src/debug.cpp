
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

//' @export
// [[Rcpp::export]]
SEXP debug(SEXP data){


  auto dat = beachmat::create_numeric_matrix(data);

  const size_t& nc=dat->get_ncol();
  const size_t& nr=dat->get_nrow();


  Rcpp::NumericMatrix final(nr,nc);

  for(int i =0;i<nr;i++){
    for(int j=0; j<nc;j++){
      Rcpp::NumericVector tmp1(nr);

      dat->get(i,j);

      final(i,j)=dat->get(i,j);
    }

  }




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

  return final;
}

