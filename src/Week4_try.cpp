#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include <typeinfo>
#include "Rinternals.h"
#include "Rcpp.h"

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
String matrix_mean(SEXP data){

  int a;

  const std::type_info& type_a = typeid(a);
  const std::type_info& type_data = typeid(data);

  if( type_a.hash_code() == type_data.hash_code()){
    return "same";
  }else{
    return "wrong";
  }

}




