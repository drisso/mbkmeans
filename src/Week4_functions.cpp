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


// [[Rcpp::export]]
SEXP transfer_data(SEXP data){
  int matrix_type = TYPEOF(data);

  if(matrix_type ==13){

    auto dat=beachmat::create_integer_matrix(data);

    Rcpp::RObject dat1 = dat->yield();

    SEXP final_matrix = as<IntegerMatrix>(dat1);

    return final_matrix;

  }else if(matrix_type==14){

    // returns a std::unique_ptr<beachmat::numeric_matrix> object
    auto dat = beachmat::create_numeric_matrix(data);

    Rcpp::RObject dat1 = dat->yield();

    SEXP final_matrix = as<NumericMatrix>(dat1);

    return final_matrix;

  }else if(matrix_type==25){

    Rcpp::RObject h5seed=get_safe_slot(data, "seed");
    Rcpp::RObject first_val=get_safe_slot(h5seed, "first_val");
    int identify_hdf5matrix= first_val.sexp_type();

    if(identify_hdf5matrix==13){

      auto dat=beachmat::create_integer_matrix(data);

      //Rcpp::RObject dat1 = dat->yield();

      //SEXP final_matrix = as<Rcpp::IntegerMatrix>(dat1);

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

    }else if(identify_hdf5matrix==14){

      // returns a std::unique_ptr<beachmat::numeric_matrix> object
      auto dat = beachmat::create_numeric_matrix(data);

      //Rcpp::RObject dat1 = dat->yield();

      //SEXP final_matrix = as<NumericMatrix>(dat1);

      const size_t& nc=dat->get_ncol();
      const size_t& nr=dat->get_nrow();


      Rcpp::NumericMatrix final_matrix(nr,nc);

      for(int i =0;i<nr;i++){
        for(int j=0; j<nc;j++){
          Rcpp::NumericVector tmp1(nr);

          dat->get(i,j);

          final_matrix(i,j)=dat->get(i,j);
        }

      }


      return final_matrix;
    }else{
      return 0;
    }


return 0;
  }


return 0;
}