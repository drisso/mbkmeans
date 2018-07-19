#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]


#ifdef _OPENMP
#include <omp.h>
#endif

#include "Week4_functions.h"



// [[Rcpp::export]]
int type_integer_index(SEXP data){
  auto iptr=beachmat::create_integer_matrix(data);
  return 0;
}

// [[Rcpp::export]]
int type_numeric_index(SEXP data){
  auto dptr=beachmat::create_numeric_matrix(data);
  return 0;
}


// [[Rcpp::export]]
std::string make_to_string(const Rcpp::RObject& str) {
  Rcpp::StringVector as_str(str);
  if (as_str.size()!=1) {
    throw std::runtime_error("input RObject should contain a single string");
  }
  return Rcpp::as<std::string>(as_str[0]);
}


// [[Rcpp::export]]
std::string get_class(const Rcpp::RObject& incoming) {
  if (!incoming.isObject()) {
    throw std::runtime_error("object has no class attribute");
  }
  return make_to_string(incoming.attr("class"));
}



// [[Rcpp::export]]
Rcpp::RObject get_safe_slot(const Rcpp::RObject& incoming, const std::string& slotname) {
  if (!incoming.hasSlot(slotname)) {
    std::stringstream err;
    err << "no '" << slotname << "' slot in the " << get_class(incoming) << " object";
    throw std::runtime_error(err.str().c_str());
  }
  return incoming.slot(slotname);
}
