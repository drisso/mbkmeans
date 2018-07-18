#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include <typeinfo>
#include "Rinternals.h"
#include "Rcpp.h"

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
int matrix_mean(SEXP data){

  int matrix_type = TYPEOF(data);

  //type of integer matrix is 13
  //type of numeric matrix is 14
  if(matrix_type==13){
    return 1;
  }else if(matrix_type==14){
    return 2;
  }else{
    return 3;
  }

}




