#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]


#ifdef _OPENMP
#include <omp.h>
#endif

#include "functions.h"


// [[Rcpp::export]]
SEXP transfer_data(SEXP data){

  auto matrix_type=beachmat::find_sexp_type(data);

  if(matrix_type== INTSXP){

    auto dat=beachmat::create_integer_matrix(data);

    const size_t& nc=dat->get_ncol();

    const size_t& nr=dat->get_nrow();

    Rcpp::IntegerMatrix final_matrix(nr,nc);

    for(int i =0;i<nr;i++){

      for(int j=0; j<nc;j++){

        Rcpp::IntegerVector tmp1(nr);

        dat->get(i,j);

        final_matrix(i,j)=dat->get(i,j);
      }

    }

    return final_matrix;

  }else if(matrix_type== REALSXP){

    auto dat = beachmat::create_numeric_matrix(data);

    const size_t& nc=dat->get_ncol();

    const size_t& nr=dat->get_nrow();

    Rcpp::IntegerMatrix final_matrix(nr,nc);

    for(int i =0;i<nr;i++){

      for(int j=0; j<nc;j++){

        Rcpp::IntegerVector tmp1(nr);

        dat->get(i,j);

        final_matrix(i,j)=dat->get(i,j);
      }

    }

    return final_matrix;

  }else{
    Rcpp::stop("The type of matrix is wrong.");
  }
}
