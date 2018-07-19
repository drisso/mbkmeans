#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include <typeinfo>
#include "Rinternals.h"

#include "Rcpp.h"


#include "Week4_functions.h"

#include "H5Cpp.h"
using namespace H5;

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
int matrix_mean(SEXP data){

  int matrix_type = TYPEOF(data);

  if(matrix_type==13){
    //type of input is integer matrix

    // creates a std::unique_ptr<beachmat::integer_matrix>
    auto iptr=beachmat::create_integer_matrix(data);

    const size_t& nc=iptr->get_ncol();
    const size_t& nr=iptr->get_nrow();

    Rcpp::NumericVector sums(nc);

    for(int i=0; i<nc; i++) {
      Rcpp::NumericVector tmp(nr);

      iptr->get_col(i, tmp.begin());

      sums[i] = sum(tmp);
    }

    return 1;


  }else if(matrix_type==14){
    //type of input is numeric matrix

    // returns a std::unique_ptr<beachmat::numeric_matrix> object
    auto dptr = beachmat::create_numeric_matrix(data);

    const size_t& nc=dptr->get_ncol();
    const size_t& nr=dptr->get_nrow();

    Rcpp::NumericVector sums(nc);

    for(int i=0; i<nc; i++) {
      Rcpp::NumericVector tmp(nr);

      dptr->get_col(i, tmp.begin());

      sums[i] = sum(tmp);
    }

    return 2;


  }else if(matrix_type==25){

    Rcpp::RObject h5seed=get_safe_slot(data, "seed");
    Rcpp::RObject first_val=get_safe_slot(h5seed, "first_val");
    char m= first_val.sexp_type();

    return 0;



  }else{
    //type os input is other
    int type;

    type=TYPEOF(data);

    return type;
  }

}




