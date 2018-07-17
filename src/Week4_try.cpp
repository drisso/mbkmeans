#include "beachmat/numeric_matrix.h"
#include "Rcpp.h"

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
Rcpp::RObject matrix_mean(SEXP data){

  auto dat = beachmat::create_numeric_matrix(data);


  return dat->yield();
}




