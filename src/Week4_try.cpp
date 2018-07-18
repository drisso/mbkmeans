#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include <typeinfo>
#include "Rinternals.h"
#include "Rcpp.h"

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
String matrix_mean(SEXP data){


  return TYPEOF(data);

}




