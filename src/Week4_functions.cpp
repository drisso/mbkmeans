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