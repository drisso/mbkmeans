#include "beachmat/numeric_matrix.h"
#include "Rcpp.h"

using namespace Rcpp;

//' Compute column sums from matrix using beachmat
//'
//' Input is a numeric matrix of any of the types supported by beachmat.
//'
//' @param dmat a matrix.
//'
//' @return A vector with column means.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector beachmat_colSums(SEXP dmat) {

  // returns a std::unique_ptr<beachmat::numeric_matrix> object
  auto mat = beachmat::create_numeric_matrix(dmat);

  const size_t& nc=mat->get_ncol();
  const size_t& nr=mat->get_nrow();

  Rcpp::NumericVector sums(nc);

  for(int i=0; i<nc; i++) {
    Rcpp::NumericVector tmp(nr);

    mat->get_col(i, tmp.begin());

    sums[i] = sum(tmp);
  }

  return sums;
}
