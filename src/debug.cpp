#define ARMA_DONT_PRINT_ERRORS
# include <RcppArmadillo.h>
# include <ClusterRHeader.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::depends(ClusterR)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
//#include "functions.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <math.h>

using namespace arma;
using namespace clustR;
template<typename T>
int get_nrow(const T& data){

    auto matrix_type=beachmat::find_sexp_type(data);

    if(matrix_type== INTSXP){
        auto final_matrix=beachmat::create_integer_matrix(data);
        //const size_t& nc = final_matrix->get_ncol();
        const size_t& nr = final_matrix->get_nrow();
        int n_row = nr;
        return n_row;

    }else if(matrix_type== REALSXP){
        auto final_matrix=beachmat::create_numeric_matrix(data);
        //const size_t& nc = final_matrix->get_ncol();
        const size_t& nr = final_matrix->get_nrow();
        int n_row = nr;
        return n_row;
    }else{

        return 0;
    }
}



//get the number of columns
template<typename T>
int get_ncol(const T& data){

    auto matrix_type=beachmat::find_sexp_type(data);

    if(matrix_type== INTSXP){
        auto final_matrix=beachmat::create_integer_matrix(data);
        const size_t& nc = final_matrix->get_ncol();
        //const size_t& nr = final_matrix->get_nrow();
        int n_col = nc;
        return n_col;

    }else if(matrix_type== REALSXP){
        auto final_matrix=beachmat::create_numeric_matrix(data);
        const size_t& nc = final_matrix->get_ncol();
        //const size_t& nr = final_matrix->get_nrow();
        int n_col = nc;
        return n_col;
    }else{

        return 0;
    }
}


template<typename T1>
Rcpp::NumericVector compute_wcss(Rcpp::NumericVector clusters, Rcpp::NumericMatrix cent, const T1& data){

    int nclusters = cent.nrow(); // number of clusters
    int nobs = clusters.size(); // number of obs
    Rcpp::NumericVector labels = Rcpp::sort_unique(clusters); // unique labels

    int data_n_cols = get_ncol(data);
    Rcpp::NumericMatrix wcss_rowdata(1,data_n_cols);
    Rcpp::NumericVector wcss_final(nclusters);

    auto matrix_type=beachmat::find_sexp_type(data);


    if(matrix_type == INTSXP){

        auto wcss_matrix=beachmat::create_integer_matrix(data);
        Rcpp::IntegerVector tmp(data_n_cols);

        for(int i=0; i <nclusters;i++){

            for(int j =0; j<nobs;j++){

                if(clusters[j] == labels[i]){

                    wcss_matrix -> get_row(j,tmp.begin());

                    wcss_rowdata.row(0) = tmp;

                    Rcpp::NumericVector z  = pow((wcss_rowdata.row(0)-cent.row(i)),2);

                    wcss_final[i] += sum(z);
                }
            }
        }
    }else if(matrix_type ==REALSXP){

        auto wcss_matrix=beachmat::create_numeric_matrix(data);
        Rcpp::NumericVector tmp(data_n_cols);

        for(int i=0; i <nclusters;i++){

            for(int j =0; j<nobs;j++){

                if(clusters[j] == labels[i]){

                    wcss_matrix -> get_row(j,tmp.begin());

                    wcss_rowdata.row(0) = tmp;

                    Rcpp::NumericVector z  = pow((wcss_rowdata.row(0)-cent.row(i)),2);

                    wcss_final[i] += sum(z);
                }
            }
        }
    }

    return wcss_final;
}



//' @export
// [[Rcpp::export]]
Rcpp::NumericVector debug(Rcpp::NumericVector clusters, arma::mat cent,SEXP data){

    Rcpp::NumericVector z = compute_wcss(clusters,Rcpp::wrap(cent),data);

    return z;
}
