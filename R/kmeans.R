#' @importFrom Rcpp sourceCpp
#' @import ClusterR
#' @useDynLib mbkmeans, .registration = TRUE
#'
NULL


#' @title Mini-Batch k-means for large single cell sequencing data
#'
#' @description This is an implementation of the mini-batch k-means algorithm of
#'   Sculley (2010) for large single cell sequencing data with the
#'   dimensionality reduction results as input in in the reducedDim() slot.
#'
#' @details The implementation is largely based on the
#'   \code{\link[ClusterR]{MiniBatchKmeans}} function of the \code{ClusterR}
#'   package. The contribution of this package is to provide support for on-disk
#'   data representations such as HDF5, through the use of \code{DelayedMatrix}
#'   and \code{HDF5Matrix} objects, as well as for sparse data representation
#'   through the classes of the \code{Matrix} package.
#'   We also provide high-level methods for objects of class
#'   \code{SummarizedExperiment}, \code{SingleCellExperiment}, and
#'   \code{LinearEmbeddingMatrix}.
#'
#' @param x The object on which to run mini-batch k-means. It can be a
#'   matrix-like object (e.g., matrix, Matrix, DelayedMatrix, HDF5Matrix) with
#'   genes in the rows and samples in the columns.
#'   Specialized methods are defined for SummarizedExperiment and
#'   SingleCellExperiment.
#' @param ... Arguments to pass to the matrix method.
#' @return A list with the following attributes: centroids, WCSS_per_cluster,
#'   best_initialization, iters_per_initialization.
#' @rdname mbkmeans
#' @importFrom methods is
#' @importFrom SummarizedExperiment assay
#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @references Sculley. Web-Scale K-Means Clustering. WWW 2010, April 26–30, 2010, Raleigh, North Carolina, USA. ACM 978-1-60558-799-8/10/04.
#' @author Lampros Mouselimis and Yuwei Ni
#' @examples
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(matrix(rnorm(100), ncol=10))
#' mbkmeans(se, clusters = 2)
setMethod(
  f = "mbkmeans",
  signature = signature(x = "SummarizedExperiment"),
  definition = function(x, whichAssay = 1, ...){
     mbkmeans(assay(x, whichAssay), ...)
  })

#' @rdname mbkmeans
#' @export
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @param reduceMethod Name of dimensionality reduction results to use as input
#'   to mini-batch k-means. Set to NA to use the full matrix.
#' @param whichAssay The assay to use as input to mini-batch k-means. If x is a SingleCellExperiment, this is ignored unless \code{reduceMethod = NA}.
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(matrix(rnorm(100), ncol=10))
#' mbkmeans(sce, clusters = 2, reduceMethod = NA)
setMethod(
  f = "mbkmeans",
  signature = signature(x = "SingleCellExperiment"),
  definition = function(x, reduceMethod = "PCA", whichAssay = 1, ...)
  {

    if(is.na(reduceMethod)){
      if(NCOL(x)>10000)
        message("Note that you are running kmeans with more than 10,000 cells using all of the dimensions. You might consider running a dimensionality reduction step first.")
      fit <- mbkmeans(assay(x, whichAssay), ...)
    }
    else{
      if(is.null(reducedDimNames(x))){
        stop("There are no dimensionality reduction results
             stored in this object. Use reduceDims() to store
             dimensionality reduction results.")
      }
      if(!(reduceMethod %in% reducedDimNames(x))){
        stop("The argument reduceMethod does not match one
             of the reducedDimNames() in this object. Use
             reducedDimNames() to see what names are in this object.")

      }
      fit <- mbkmeans(t(reducedDim(x, reduceMethod)), ...)

      }

    return(fit)
    })

#' @rdname mbkmeans
#' @export
#' @importClassesFrom SingleCellExperiment LinearEmbeddingMatrix
#' @importFrom SingleCellExperiment sampleFactors
setMethod(
  f = "mbkmeans",
  signature = signature(x = "LinearEmbeddingMatrix"),
  definition = function(x, ...)
  {
    mbkmeans(t(sampleFactors(x)), ...)
  })

#'@rdname mbkmeans
#'@export
#'@importClassesFrom DelayedArray DelayedMatrix
#'@param clusters the number of clusters
#'@param batch_size the size of the mini batches
#'@param num_init number of times the algorithm will be run with different
#'  centroid seeds
#'@param max_iters the maximum number of clustering iterations
#'@param init_fraction percentage of data to use for the initialization
#'  centroids (applies if initializer is \emph{kmeans++} ). Should be a float
#'  number between 0.0 and 1.0.
#'@param initializer the method of initialization. One of \emph{kmeans++} and
#'  \emph{random}. See details for more information
#'@param early_stop_iter continue that many iterations after calculation of the
#'  best within-cluster-sum-of-squared-error
#'@param verbose either TRUE or FALSE, indicating whether progress is printed
#'  during clustering
#'@param CENTROIDS a matrix of initial cluster centroids. The rows of the
#'  CENTROIDS matrix should be equal to the number of clusters and the columns
#'  should be equal to the columns of the data
#'@param tol a float number. If, in case of an iteration (iteration > 1 and
#'  iteration < max_iters) 'tol' is greater than the squared norm of the
#'  centroids, then kmeans has converged
#'@param seed integer value for random number generator (RNG)
#'@return a list with the following attributes: centroids, WCSS_per_cluster,
#'  best_initialization, iters_per_initialization
#'@details This function performs k-means clustering using mini batches.
#'
#'\strong{kmeans++}: kmeans++ initialization. Reference :
#'http://theory.stanford.edu/~sergei/papers/kMeansPP-soda.pdf AND
#'http://stackoverflow.com/questions/5466323/how-exactly-does-k-means-work
#'
#'\strong{random}: random selection of data rows as initial centroids
#'
#'@references https://github.com/mlampros/ClusterR
#'
#'@examples
#'x<-matrix(rnorm(100), ncol=10)
#'mbkmeans(x,clusters = 3)
#'
setMethod(
  f = "mbkmeans",
  signature = signature(x ="ANY"),
  definition = function(x, clusters, batch_size = blocksize(x),
                        max_iters =10, num_init = 1,
                        init_fraction = .25, initializer = "kmeans++",
                        early_stop_iter = 10, verbose = FALSE,
                        CENTROIDS = NULL, tol = 1e-4, seed = 1)
  {

    if(!is(x, "matrix") & !is(x, "Matrix") & !is(x, "HDF5Matrix") &
       !is(x, "DelayedMatrix")) {

      stop("x is not of a supported type")

    } else {

      fit <- mini_batch(t(x), clusters, batch_size, max_iters, num_init,
                        init_fraction, initializer, early_stop_iter,
                        verbose, CENTROIDS, tol, seed)

    }

    return(fit)
})

