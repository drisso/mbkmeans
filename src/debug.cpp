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
arma::rowvec wcss_result(Rcpp::NumericVector clusters, arma::mat cent,const T1& data){


    ClustHeader clust_header;

    int centers_row = cent.n_rows;

    //int cluster_length = clusters.size();

    int data_row = get_nrow(data);

    Rcpp::NumericMatrix wcss_tmp(data_row,1);

    arma::rowvec wcss_final(centers_row);

    auto matrix_type=beachmat::find_sexp_type(data);

    if(matrix_type == INTSXP){

        auto wcss_matrix=beachmat::create_integer_matrix(data);

        int data_n_rows = get_nrow(data);

        int data_n_cols = get_ncol(data);

        Rcpp::IntegerVector tmp(data_n_cols);

        Rcpp::IntegerMatrix wcss_rowdata(1,data_n_cols);

        for(int i=0; i <centers_row;i++){

            for(int j =0; j<data_row;j++){

                if(clusters[j] == i+1){

                    wcss_matrix -> get_row(j,tmp.begin());

                    wcss_rowdata.row(0) = tmp;

                    Rcpp::NumericVector z  = pow(wcss_rowdata.row(0)-cent[i],2);

                    wcss_tmp.row(j)= z;

                }else{

                    wcss_tmp.row(j) ==0;
                }

            }
            wcss_final[i] = sum(wcss_tmp);
        }
    }else if(matrix_type ==REALSXP){

        auto wcss_matrix=beachmat::create_numeric_matrix(data);

        int data_n_rows = get_nrow(data);

        int data_n_cols = get_ncol(data);

        Rcpp::NumericVector tmp(data_n_cols);

        Rcpp::NumericMatrix wcss_rowdata(1,data_n_cols);

        for(int i=0; i <centers_row;i++){

            for(int j =0; j<data_row;j++){

                if(clusters[j] == i+1){

                    wcss_matrix -> get_row(j,tmp.begin());

                    wcss_rowdata.row(0) = tmp;

                    Rcpp::NumericVector z  = pow(wcss_rowdata.row(0)-cent[i],2);

                    wcss_tmp.row(j)= z;
                }else{

                    wcss_tmp.row(j) == 0;

                }

            }
            wcss_final[i] = sum(wcss_tmp);
        }
    }

    return wcss_final;


}



//' @export
// [[Rcpp::export]]
arma::rowvec debug(Rcpp::NumericVector clusters, arma::mat cent,SEXP data){

    arma::rowvec z = wcss_result(clusters,cent,data);

    return z;
}
