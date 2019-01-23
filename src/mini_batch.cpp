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

//get the number of rows
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

//calculate the clusters
template<typename T>
arma::rowvec clusters_WCSS(const T&data,Rcpp::NumericMatrix CENTROIDS){

    ClustHeader clust_header;

    auto matrix_type=beachmat::find_sexp_type(data);

    if(matrix_type == INTSXP){

        auto final_matrix=beachmat::create_integer_matrix(data);
        int data_n_rows = get_nrow(data);
        int data_n_cols = get_ncol(data);

        Rcpp::NumericMatrix dat_final(1,data_n_cols);
        Rcpp::IntegerVector tmp(data_n_cols);
        arma::rowvec CLUSTERS(data_n_rows);
        int centrod_n_row = get_nrow(CENTROIDS);

        arma::vec tmp_vec(centrod_n_row);

        for (unsigned int j = 0; j < data_n_rows; j++) {

            final_matrix->get_row(j, tmp.begin());
            dat_final.row(0) = tmp;


            //arma::mat data_final = Rcpp::as<arma::mat>(dat_final);
            //arma::vec tmp_vec = clust_header.WCSS(arma::conv_to< arma::rowvec >::from(data_final.row(0)), CENTROIDS);                 // returns a rowvec with the SSE for each cluster
            for(unsigned int i =0; i<centrod_n_row;i++){

                tmp_vec(i)= Rcpp::sum(Rcpp::pow(dat_final.row(0) - CENTROIDS.row(i),2));

            }
            //soft_CLUSTERS.row(j) = arma::conv_to< arma::rowvec >::from(tmp_vec);

            int tmp_idx = clust_header.MinMat(tmp_vec);                                                                        // returns the index of the tmp_vec with the lowest SSE
            CLUSTERS(j) = tmp_idx+1;
        }
        return CLUSTERS;


    }else if(matrix_type ==REALSXP){
        auto final_matrix=beachmat::create_numeric_matrix(data);
        int data_n_rows = get_nrow(data);
        int data_n_cols = get_ncol(data);
        Rcpp::NumericMatrix dat_final(1,data_n_cols);
        Rcpp::NumericVector tmp(data_n_cols);
        arma::rowvec CLUSTERS(data_n_rows);
        int centrod_n_row = get_nrow(CENTROIDS);

        arma::vec tmp_vec(centrod_n_row);

        //arma::vec tmp_vec2(centrod_n_row);

        for (unsigned int j = 0; j < data_n_rows; j++) {

            final_matrix->get_row(j, tmp.begin());
            dat_final.row(0) = tmp;

            //arma::mat data_final = Rcpp::as<arma::mat>(dat_final);
            //arma::vec tmp_vec = clust_header.WCSS(arma::conv_to< arma::rowvec >::from(data_final.row(0)), CENTROIDS);                 // returns a rowvec with the SSE for each cluster

            for(unsigned int i =0; i<centrod_n_row;i++){

                tmp_vec(i)= Rcpp::sum(Rcpp::pow(dat_final.row(0) - CENTROIDS.row(i),2));

            }

            //soft_CLUSTERS.row(j) = arma::conv_to< arma::rowvec >::from(tmp_vec);

            int tmp_idx = clust_header.MinMat(tmp_vec);                                                                        // returns the index of the tmp_vec with the lowest SSE
            CLUSTERS(j) = tmp_idx+1;
        }
        return CLUSTERS;
    }

}

//' Predict_mini_batch
//'
//' Prediction function for Mini-batch-k-means for in-memory, delayed, and on-disk matrices
//'
//'
//'@param data matrix, DelayedMatrix, or HDF5Matrix containing numeric or
//'  integer data (obseravtions in rows, variables in columns)
//'@param CENTROIDS a matrix of initial cluster centroids. The rows of the
//'  CENTROIDS matrix should be equal to the number of clusters and the columns
//'  should equal the columns of the data.
//'@return it returns a vector with the clusters.
//'@details
//'
//'This function takes the data and the output centroids and returns the
//'clusters.
//'
//'This implementation relies very heavily on the
//'\code{\link[ClusterR]{MiniBatchKmeans}} implementation. We provide the
//'ability to work with DelayedMatrix and HDF5Matrix through the \code{beachmat}
//'library.
//'
//'@author Lampros Mouselimis and Yuwei Ni
//'
//'@examples
//'data(iris)
//'km = mini_batch(as.matrix(iris[,1:4]), clusters = 3,
//'                batch_size = 10, max_iters = 10)
//'clusters = predict_mini_batch(as.matrix(iris[,1:4]),
//'                              CENTROIDS = km$centroids)
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector predict_mini_batch(SEXP data, Rcpp::NumericMatrix CENTROIDS, bool fuzzy = false, double eps = 1.0e-6) {


    //arma::mat CENTROIDS1;


    //if (CENTROIDS.isNotNull()) {

    // CENTROIDS1 = Rcpp::as<arma::mat>(CENTROIDS);

    //}

    //SEXP trans_data = transfer_data(data);

    // arma::mat dat_final = Rcpp::as<arma::mat>(trans_data);

    int data_n_rows = get_nrow(data);

    int data_n_cols = get_ncol(data);

    arma::rowvec CLUSTERS(data_n_rows);

    //arma::mat soft_CLUSTERS(data_n_rows, CENTROIDS1.n_rows);

    CLUSTERS = clusters_WCSS(data,CENTROIDS);

    //for (unsigned int j = 0; j < data_n_rows; j++) {

    //Rcpp::NumericVector tmp(data_n_cols);

    //auto final_matrix=beachmat::create_numeric_matrix(data);

    //final_matrix->get_row(j, tmp.begin());

    //dat_final.row(j) = tmp;

    //arma::mat data_final = Rcpp::as<arma::mat>(dat_final);

    //arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(data_final.row(j)), CENTROIDS1);

    //arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(dat_final.row(j)), CENTROIDS1);                 // returns a rowvec with the SSE for each cluster

    //soft_CLUSTERS.row(j) = arma::conv_to< arma::rowvec >::from(tmp_vec);

    //int tmp_idx = MinMat(tmp_vec);                                                                        // returns the index of the tmp_vec with the lowest SSE

    //CLUSTERS(j) = tmp_idx+1;
    // }

    if (fuzzy) {

        //     arma::mat fuzzy_mat(soft_CLUSTERS.n_rows, soft_CLUSTERS.n_cols);

        //     for (unsigned int i = 0; i < soft_CLUSTERS.n_rows; i++) {

        //       fuzzy_mat.row(i) = norm_fuzzy(arma::conv_to< arma::rowvec >::from(soft_CLUSTERS.row(i)), eps);
        //     }

        //     return Rcpp::List::create(Rcpp::Named("clusters") = CLUSTERS, Rcpp::Named("fuzzy_clusters") = fuzzy_mat);

        Rcpp::stop("fuzzy clustering is currently not implemented.");

    }else {

        //return CLUSTERS;

        Rcpp::NumericVector final_cluster = Rcpp::wrap(CLUSTERS);

        final_cluster.attr("dim") = R_NilValue;

        return final_cluster;
    }
}

//subset the matrix and select randomly rows(nrows = init_fraction*total)
template<typename T1, typename T2>
SEXP subset_matrix(const T1& data, const T2& init_fraction_value){

  ClustHeader clust_header;

  const size_t& nc = data->get_ncol();
  const size_t& nr = data->get_nrow();
  int fract = std::ceil(nr*init_fraction_value);

  arma::uvec init = arma::conv_to< arma::uvec >::from(clust_header.sample_vec(fract, 0, nr - 1, false));

  arma::uvec samp_init = arma::sort(init);

  Rcpp::NumericMatrix submat(samp_init.n_rows, nc);
  Rcpp::NumericVector tmp(nc);

  for(int i=0; i<samp_init.n_rows; i++){
    data->get_row(samp_init[i], tmp.begin());
    submat.row(i) = tmp;
  }

  return submat;

}


//subset the matrix and select randomly rows(nrow  = cluster)
template<typename T1>
SEXP subset_matrix_random(const T1& data, int cluster){

  ClustHeader clust_header;

  const size_t& nc = data->get_ncol();
  const size_t& nr = data->get_nrow();

  arma::uvec init = arma::conv_to< arma::uvec >::from(clust_header.sample_vec(cluster, 0, nr - 1, false));

  arma::uvec samp_init = arma::sort(init);

  Rcpp::NumericMatrix submat(samp_init.n_rows, nc);
  Rcpp::NumericVector tmp(nc);

  for(int i=0; i<samp_init.n_rows; i++){
    data->get_row(samp_init[i], tmp.begin());
    submat.row(i) = tmp;
  }

  return submat;

}


//calculation right WCSS
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


//'
//' Mini_batch
//'
//' Mini-batch-k-means for both matrix and HDF5Matrix
//'
//'@param data numeric matrix or integer matrix or HDF5Matrix
//'@param clusters the number of clusters
//'@param batch_size the size of the mini batches
//'@param num_init number of times the algorithm will be run with different centroid seeds
//'@param max_iters the maximum number of clustering iterations
//'@param init_fraction percentage of data to use for the initialization centroids (applies if initializer is \emph{kmeans++} ). Should be a float number between 0.0 and 1.0.
//'@param initializer the method of initialization. One of \emph{kmeans++} and \emph{random}. See details for more information
//'@param early_stop_iter continue that many iterations after calculation of the best within-cluster-sum-of-squared-error
//'@param verbose either TRUE or FALSE, indicating whether progress is printed during clustering
//'@param CENTROIDS a matrix of initial cluster centroids. The rows of the CENTROIDS matrix should be equal to the number of clusters and the columns should be equal to the columns of the data
//'@param tol a float number. If, in case of an iteration (iteration > 1 and iteration < max_iters) 'tol' is greater than the squared norm of the centroids, then kmeans has converged
//'@param seed integer value for random number generator (RNG)
//'@return a list with the following attributes: centroids, WCSS_per_cluster, best_initialization, iters_per_initialization
//'@details
//'This function performs k-means clustering using mini batches.
//'
//'\strong{kmeans++}: kmeans++ initialization. Reference : http://theory.stanford.edu/~sergei/papers/kMeansPP-soda.pdf AND http://stackoverflow.com/questions/5466323/how-exactly-does-k-means-work
//'
//'\strong{random}: random selection of data rows as initial centroids
//'
//'@references
//'https://github.com/mlampros/ClusterR
//'
//'@examples
//'data = matrix(1:30,nrow = 10)
//'mini_batch(data, 2, 10, 10)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List mini_batch(SEXP data, int clusters, int batch_size, int max_iters, int num_init = 1, double init_fraction = 1.0, std::string initializer = "kmeans++",

                      bool wcss_show = false, int early_stop_iter = 10, bool verbose = false, Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, double tol = 1e-4,  int seed = 1){

  ClustHeader clust_header;

  clust_header.set_seed(seed);             // R's RNG

  int dat_n_rows = get_nrow(data);    // use the template function in this file

  //int dat_n_rows = data.n_rows;

  if (clusters > dat_n_rows - 2 || clusters < 1) { Rcpp::stop("the number of clusters should be at most equal to nrow(data) - 2 and never less than 1"); }

  bool flag = false;

  arma::mat CENTROIDS1;

  if (CENTROIDS.isNotNull()) {

    CENTROIDS1 = Rcpp::as<arma::mat>(CENTROIDS);

    num_init = 1;

    flag = true;                   // get centroids
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

        // if(init_fraction==1.0){
        //
        //   //SEXP centroids_matrix = transfer_data(data);
        //
        //   //update_centroids = Rcpp::as<arma::mat>(centroids_matrix);
        //
        //   SEXP trans_data = transfer_data(data);
        //
        //   arma::mat final_data = Rcpp::as<arma::mat>(trans_data);
        //
        //   update_centroids = clust_header.kmeans_pp_init(final_data, clusters, false);
        //
        // }else
        if(init_fraction <= 1.0 && init_fraction >0.0){

          Rcpp::NumericMatrix tran_data;

          auto matrix_type=beachmat::find_sexp_type(data);

          if(matrix_type == INTSXP){

            auto final_matrix=beachmat::create_integer_matrix(data);

             tran_data =  subset_matrix(final_matrix,init_fraction);

            //return tran_data;

          }else if(matrix_type ==REALSXP){

            auto final_matrix=beachmat::create_numeric_matrix(data);

             tran_data = subset_matrix(final_matrix,init_fraction);

            //return tran_data;

          }else{

            Rcpp::stop("The type of matrix is wrong");

          }

          arma::mat final_data = Rcpp::as<arma::mat>(tran_data);

          update_centroids = clust_header.kmeans_pp_init(final_data, clusters, false);

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

           tran_data_random =  subset_matrix_random(final_matrix,clusters);

          //return tran_data_random;

        }else if(matrix_type ==REALSXP){

          auto final_matrix=beachmat::create_numeric_matrix(data);

           tran_data_random = subset_matrix_random(final_matrix,clusters);

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

         batch_data_choose =  subset_matrix_random(final_matrix,batch_size);

       // return batch_data_choose;

      }else if(matrix_type ==REALSXP){

        auto final_matrix=beachmat::create_numeric_matrix(data);

         batch_data_choose = subset_matrix_random(final_matrix,batch_size);

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

        arma::vec tmp_vec = clust_header.WCSS(arma::conv_to< arma::rowvec >::from(batch_data.row(j)), update_centroids);         // returns a rowvec with the SSE for each cluster

        int tmp_idx = clust_header.MinMat(tmp_vec);                                                                              // returns the index of the tmp_vec with the lowest SSE

        total_SSE(tmp_idx) += tmp_vec(tmp_idx);             //total WCSS                                                          // assigns to total_SSE the minimum cost

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

      double tmp_norm = clust_header.squared_norm(previous_centroids - update_centroids);   //// early-stopping criterium

      double calc_cost = arma::accu(total_SSE);

      if (verbose) { Rcpp::Rcout << "iteration: " << i + 1 << "  --> total WCSS: " << calc_cost << "  -->  squared norm: " << tmp_norm << std::endl; }

      count = i;

      if (calc_cost < previous_cost) {

        previous_cost = calc_cost;     //calc_cost is total WCSS

        increment_early_stop = 0;

        output_centroids = update_centroids;                // assign end-centroids and SSE when WCSS is minimal

        output_SSE = total_SSE;
      }

      if (calc_cost > previous_cost) {

        increment_early_stop += 1;
      }

      if (tmp_norm < tol || i == max_iters - 1 || increment_early_stop == early_stop_iter - 1) {

        output_centroids = update_centroids;

        // assign end-centroids and SSE when early_stop_iter == increment_early_stop [ repeated calculation of adjusted rand index shows slightly better results for the previous case ]

         output_SSE = total_SSE;

        break;
      }
    }

    if (arma::accu(output_SSE) < arma::accu(bst_WCSS)) {   //output_SSE is total wcss, bst_WCSS is WCSS

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

  // Rcpp::Environment package_env("package:mbkmeans");

  // Rcpp::Function rfunction = package_env["predict_mini_batch"];

  // Rcpp::List clusterfinal;

  // clusterfinal= rfunction(data,centers_out);

  Rcpp::NumericVector clusterfinal = predict_mini_batch(data, Rcpp::wrap(centers_out));

  Rcpp::NumericMatrix centrod = Rcpp::wrap(centers_out);

  arma::rowvec wcss_final = compute_wcss(clusterfinal,centrod,data);

  //Rcpp::NumericMatrix wcss_final = wcss_result(clusterfinal,centers_out,data);

  // bool wcss_tmp;
  //
  // int centers_row = centers_out.n_rows;
  //
  // int cluster_length = clusterfinal.size();
  //
  // arma::rowvec wcss_final(centers_row);
  //
  // for(int i =0; i< centers_row; i++){
  //
  //     for(int j =0; j<cluster_length;j++){
  //
  //       wcss_tmp = (clusterfinal[j] = i);
  //
  //       if(wcss_tmp == 1){
  //
  //           wcss_final[i] = sum(pow(clusterfinal[j] - centers_out[i],2));
  //
  //       }
  //
  //     }
  // }
  //
  //
  Rcpp::NumericVector clusterfinal_2 = predict_mini_batch(data, Rcpp::wrap(centers_out));

  //arma::rowvec CLUSTERS;

  //Rcpp::NumericMatrix centers;

  //CLUSTERS = predict_mini_batch(data,centers);

  if(wcss_show == true){

      return Rcpp::List::create(Rcpp::Named("centroids") = centers_out,
                                //Rcpp::Named("WCSS_per_cluster") = bst_WCSS,
                                Rcpp::Named("WCSS_per_cluster") = wcss_final,
                                Rcpp::Named("best_initialization") = end_init,
                                Rcpp::Named("iters_per_initialization") = iter_before_stop,
                                Rcpp::Named("Clusters") = clusterfinal_2);

  }else{

      return Rcpp::List::create(Rcpp::Named("centroids") = centers_out,
                                //Rcpp::Named("WCSS_per_cluster") = bst_WCSS,
                                //Rcpp::Named("WCSS_per_cluster") = wcss_final,
                                Rcpp::Named("best_initialization") = end_init,
                                Rcpp::Named("iters_per_initialization") = iter_before_stop,
                                Rcpp::Named("Clusters") = clusterfinal_2);
  }

}


