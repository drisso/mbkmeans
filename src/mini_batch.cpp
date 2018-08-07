#define ARMA_DONT_PRINT_ERRORS
# include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include "utils_rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include "Week4_functions.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <math.h>

using namespace arma;


//get the number of rows
template<typename T>
int get_nrow(const T& data){
  auto matrix_type=beachmat::find_sexp_type(data);
  if(matrix_type== INTSXP){
    auto final_matrix=beachmat::create_integer_matrix(data);
    //const size_t& nc = final_matrix->get_ncol();
    const size_t& nr = final_matrix->get_nrow();
    int n_col = nr;
    return n_col;

  }else if(matrix_type== REALSXP){
    auto final_matrix=beachmat::create_numeric_matrix(data);
    //const size_t& nc = final_matrix->get_ncol();
    const size_t& nr = final_matrix->get_nrow();
    int n_col = nr;
    return n_col;
  }else{

    return 0;
  }
}


//shuffle the matrix and select randomly rows(nrows = init_fraction*total)
template<typename T1, typename T2>
SEXP shuffle_matrix(const T1& data, const T2& init_fraction_value){

  const size_t& nc = data->get_ncol();
  const size_t& nr = data->get_nrow();
  int fract = std::ceil(nr*init_fraction_value);

  arma::uvec init = arma::conv_to< arma::uvec >::from(sample_vec(fract, 0, nr - 1, false));

  arma::uvec samp_init = arma::sort(init);

  Rcpp::NumericMatrix submat(samp_init.n_rows, nc);
  Rcpp::NumericVector tmp(nc);

  for(int i=0; i<samp_init.n_rows; i++){
    data->get_row(samp_init[i], tmp.begin());
    submat.row(i) = tmp;
  }

  return submat;

}


//shuffle the matrix and select randomly rows(nrow  = cluster)
template<typename T1>
SEXP shuffle_matrix_random(const T1& data, int cluster){
  const size_t& nc = data->get_ncol();
  const size_t& nr = data->get_nrow();

  arma::uvec init = arma::conv_to< arma::uvec >::from(sample_vec(cluster, 0, nr - 1, false));

  arma::uvec samp_init = arma::sort(init);

  Rcpp::NumericMatrix submat(samp_init.n_rows, nc);
  Rcpp::NumericVector tmp(nc);

  for(int i=0; i<samp_init.n_rows; i++){
    data->get_row(samp_init[i], tmp.begin());
    submat.row(i) = tmp;
  }

  return submat;

}


//' @export
// [[Rcpp::export]]
Rcpp::List mini_batch(SEXP data, int clusters, int batch_size, int max_iters, int num_init = 1, double init_fraction = 1.0, std::string initializer = "kmeans++",

                      int early_stop_iter = 10, bool verbose = false, Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, double tol = 1e-4, double tol_optimal_init = 0.5, int seed = 1){


  set_seed(seed);             // R's RNG

  int dat_n_rows = get_nrow(data);    // use the template function in this file

  //int dat_n_rows = data.n_rows;

  if (clusters > dat_n_rows - 2 || clusters < 1) { Rcpp::stop("the number of clusters should be at most equal to nrow(data) - 2 and never less than 1"); }

  bool flag = false;

  arma::mat CENTROIDS1;

  if (CENTROIDS.isNotNull()) {

    CENTROIDS1 = Rcpp::as<arma::mat>(CENTROIDS);

    num_init = 1;

    flag = true;
  }

  arma::mat centers_out;

  arma::rowvec iter_before_stop(num_init, arma::fill::zeros);

  int end_init = 0;

  arma::rowvec bst_WCSS(clusters);

  bst_WCSS.fill(arma::datum::inf);           // initialize WCSS to Inf, so that in first iteration it can be compared with the minimum of the 'WSSE'

  if (verbose) { Rcpp::Rcout << " " << std::endl; }

  arma::rowvec flag_exception(num_init, arma::fill::zeros);

  for (int init = 0; init < num_init; init++){

    arma::mat update_centroids;

    if(!flag){

      if(initializer== "kmeans++"){

        if(init_fraction==1.0){

          //SEXP centroids_matrix = transfer_data(data);

          //update_centroids = Rcpp::as<arma::mat>(centroids_matrix);

          SEXP trans_data = transfer_data(data);

          arma::mat final_data = Rcpp::as<arma::mat>(trans_data);

          update_centroids = kmeans_pp_init(final_data, clusters, false);

        }else if(init_fraction <1.0 && init_fraction >0.0){

          Rcpp::NumericMatrix tran_data;

          auto matrix_type=beachmat::find_sexp_type(data);

          if(matrix_type == INTSXP){

            auto final_matrix=beachmat::create_integer_matrix(data);

             tran_data =  shuffle_matrix(final_matrix,init_fraction);

            //return tran_data;

          }else if(matrix_type ==REALSXP){

            auto final_matrix=beachmat::create_numeric_matrix(data);

             tran_data = shuffle_matrix(final_matrix,init_fraction);

            //return tran_data;

          }else{

            Rcpp::stop("The type of matrix is wrong");

          }

          arma::mat final_data = Rcpp::as<arma::mat>(tran_data);

          update_centroids = kmeans_pp_init(final_data, clusters, false);

          // SEXP centroids_matrix = shuffle_matrix(data,init_fraction = init_fraction);

        }else{

          Rcpp::stop("The value of fraction should be larger than 0 and not larger than 1");

        }

      }//kmeans++ ends

      if(initializer == "random"){

        //arma::uvec samp = arma::conv_to< arma::uvec >::from(sample_vec(clusters, 0, data.n_rows - 1, false));

        //update_centroids = data.rows(samp);

        Rcpp::NumericMatrix tran_data_random;

        auto matrix_type=beachmat::find_sexp_type(data);

        if(matrix_type == INTSXP){

          auto final_matrix=beachmat::create_integer_matrix(data);

           tran_data_random =  shuffle_matrix_random(final_matrix,clusters);

          //return tran_data_random;

        }else if(matrix_type ==REALSXP){

          auto final_matrix=beachmat::create_numeric_matrix(data);

           tran_data_random = shuffle_matrix_random(final_matrix,clusters);

          //return tran_data_random;

        }else{

          Rcpp::stop("The type of matrix is wrong");

        }

        update_centroids = Rcpp::as<arma::mat>(tran_data_random);

        //update_centroids = kmeans_pp_init(final_data, clusters, false);

      }//random ends

    }else{

      update_centroids = CENTROIDS1;

    }

    arma::mat previous_centroids = update_centroids;//update the centroids

    arma::mat output_centroids;

    arma::rowvec output_SSE;

    double previous_cost = arma::datum::inf;

    int increment_early_stop = 0;

    int count = 0;

    for (int i = 0; i < max_iters; i++) {

      //select the batch_data

      Rcpp::NumericMatrix batch_data_choose;

      auto matrix_type=beachmat::find_sexp_type(data);

      if(matrix_type == INTSXP){

        auto final_matrix=beachmat::create_integer_matrix(data);

         batch_data_choose =  shuffle_matrix_random(final_matrix,batch_size);

       // return batch_data_choose;

      }else if(matrix_type ==REALSXP){

        auto final_matrix=beachmat::create_numeric_matrix(data);

         batch_data_choose = shuffle_matrix_random(final_matrix,batch_size);

        //return batch_data_choose;

      }else{

        Rcpp::stop("The type of matrix is wrong");

      }


      arma::mat batch_data = Rcpp::as<arma::mat>(batch_data_choose);

      //arma::uvec batch_idx = arma::conv_to< arma::uvec >::from(sample_vec(batch_size, 0, data.n_rows - 1, false));

      //arma::mat batch_data = data.rows(batch_idx);

      arma::rowvec total_SSE(clusters, arma::fill::zeros);

      arma::rowvec CLUSTERS(batch_data.n_rows);

      for (unsigned int j = 0; j < batch_data.n_rows; j++) {

        arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(batch_data.row(j)), update_centroids);         // returns a rowvec with the SSE for each cluster

        int tmp_idx = MinMat(tmp_vec);                                                                              // returns the index of the tmp_vec with the lowest SSE

        total_SSE(tmp_idx) += tmp_vec(tmp_idx);                                                                     // assigns to total_SSE the minimum cost

        CLUSTERS(j) = tmp_idx;
      }

      arma::rowvec cluster_counts(clusters, arma::fill::zeros);

      double eta = 0.0;

      for (unsigned int k = 0; k < batch_data.n_rows; k++) {

        int idx = CLUSTERS(k);

        cluster_counts(idx) += 1;

        eta = 1.0 / cluster_counts(idx);

        arma::rowvec tmp_row = (1.0 - eta) * arma::conv_to< arma::rowvec >::from(update_centroids.row(idx)) + eta * arma::conv_to< arma::rowvec >::from(batch_data.row(k));

        update_centroids.row(idx) = tmp_row;
      }

      double tmp_norm = squared_norm(previous_centroids - update_centroids);

      double calc_cost = arma::accu(total_SSE);

      if (verbose) { Rcpp::Rcout << "iteration: " << i + 1 << "  --> total WCSS: " << calc_cost << "  -->  squared norm: " << tmp_norm << std::endl; }

      count = i;

      if (calc_cost < previous_cost) {

        previous_cost = calc_cost;

        increment_early_stop = 0;

        output_centroids = update_centroids;                // assign end-centroids and SSE when WCSS is minimal

        output_SSE = total_SSE;
      }

      if (calc_cost > previous_cost) {

        increment_early_stop += 1;
      }

      if (tmp_norm < tol || i == max_iters - 1 || increment_early_stop == early_stop_iter - 1) {

        // output_centroids = update_centroids;            // assign end-centroids and SSE when early_stop_iter == increment_early_stop [ repeated calculation of adjusted rand index shows slightly better results for the previous case ]
        //
        // output_SSE = total_SSE;

        break;
      }
    }

    if (arma::accu(output_SSE) < arma::accu(bst_WCSS)) {

      end_init = init + 1;

      bst_WCSS = output_SSE;

      centers_out = output_centroids;
    }

    iter_before_stop(init) = count + 1;

    if (verbose) { Rcpp::Rcout << " " << std::endl; Rcpp::Rcout << "===================== end of initialization " << init + 1 << " =====================" << std::endl; Rcpp::Rcout << " " << std::endl; }

  }

  //int exc = arma::as_scalar(flag_exception(end_init - 1));                        // print warning OR stop function only if duplicates OR NA's are present in relevant output centroids [ end_init - 1 ]

  //if (exc > 0) {

  //if (initializer == "quantile_init") {

  //  std::string message = "the centroid matrix using 'quantile_init' as initializer contains duplicates for number of clusters equal to : " + std::to_string(clusters);

  //  Rcpp::warning(message);}

  //if (initializer == "optimal_init") {

  //  std::string message = "The centroid matrix using 'optimal_init' as initializer contains NA's. Thus, the 'tol_optimal_init' parameter should be (probably) decreased for number of clusters equal to : " + std::to_string(clusters);

  //  Rcpp::stop(message);
  //}
  //}

  return Rcpp::List::create(Rcpp::Named("centroids") = centers_out, Rcpp::Named("WCSS_per_cluster") = bst_WCSS,

                            Rcpp::Named("best_initialization") = end_init, Rcpp::Named("iters_per_initialization") = iter_before_stop);
}


