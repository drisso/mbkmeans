// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// predict_mini_batch
Rcpp::NumericVector predict_mini_batch(SEXP data, Rcpp::NumericMatrix CENTROIDS);
RcppExport SEXP _mbkmeans_predict_mini_batch(SEXP dataSEXP, SEXP CENTROIDSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type CENTROIDS(CENTROIDSSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_mini_batch(data, CENTROIDS));
    return rcpp_result_gen;
END_RCPP
}
// compute_wcss
Rcpp::NumericVector compute_wcss(Rcpp::NumericVector clusters, Rcpp::NumericMatrix cent, SEXP data);
RcppExport SEXP _mbkmeans_compute_wcss(SEXP clustersSEXP, SEXP centSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type clusters(clustersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type cent(centSEXP);
    Rcpp::traits::input_parameter< SEXP >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_wcss(clusters, cent, data));
    return rcpp_result_gen;
END_RCPP
}
// mini_batch
Rcpp::List mini_batch(SEXP data, int clusters, int batch_size, int max_iters, int num_init, double init_fraction, std::string initializer, bool calc_wcss, int early_stop_iter, bool verbose, Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS, double tol);
RcppExport SEXP _mbkmeans_mini_batch(SEXP dataSEXP, SEXP clustersSEXP, SEXP batch_sizeSEXP, SEXP max_itersSEXP, SEXP num_initSEXP, SEXP init_fractionSEXP, SEXP initializerSEXP, SEXP calc_wcssSEXP, SEXP early_stop_iterSEXP, SEXP verboseSEXP, SEXP CENTROIDSSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type clusters(clustersSEXP);
    Rcpp::traits::input_parameter< int >::type batch_size(batch_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type max_iters(max_itersSEXP);
    Rcpp::traits::input_parameter< int >::type num_init(num_initSEXP);
    Rcpp::traits::input_parameter< double >::type init_fraction(init_fractionSEXP);
    Rcpp::traits::input_parameter< std::string >::type initializer(initializerSEXP);
    Rcpp::traits::input_parameter< bool >::type calc_wcss(calc_wcssSEXP);
    Rcpp::traits::input_parameter< int >::type early_stop_iter(early_stop_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type CENTROIDS(CENTROIDSSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(mini_batch(data, clusters, batch_size, max_iters, num_init, init_fraction, initializer, calc_wcss, early_stop_iter, verbose, CENTROIDS, tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mbkmeans_predict_mini_batch", (DL_FUNC) &_mbkmeans_predict_mini_batch, 2},
    {"_mbkmeans_compute_wcss", (DL_FUNC) &_mbkmeans_compute_wcss, 3},
    {"_mbkmeans_mini_batch", (DL_FUNC) &_mbkmeans_mini_batch, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_mbkmeans(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
