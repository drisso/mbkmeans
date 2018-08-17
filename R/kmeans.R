#' @importFrom Rcpp sourceCpp
#' @import ClusterR
#' @useDynLib beachball, .registration = TRUE
#'
NULL


#' @title k-means for large single cell sequencing data
#'
#' @description This is a wrapper for stats::kmeans() for
#' large single cell sequencing data with the dimensionality
#' reduction results as input in in the reducedDim() slot.
#' @param x The object on which to run k-means.
#' @param reduceMethod Name of dimensionality reduction results to use as input
#'   to k-means.
#' @param whichAssay The assay to use as input to k-means. Used only if
#'   \code{reduceMethod = "none"}.
#' @param ... Arguments to pass to \code{\link[stats]{kmeans}}.
#' @return k-means output
#' @name kmeans
#' @rdname kmeans
#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods as
#' @examples
#' se <- SummarizedExperiment(matrix(rnorm(100), ncol=10))
#' kmeans(se, centers = 2, reduceMethod = "none")
setMethod(
  f = "kmeans",
  signature = signature(x = "SummarizedExperiment"),
  definition = function(x, ...){
    kmeans(as(x,"SingleCellExperiment"),...)

  })

#' @rdname kmeans
#' @export
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
setMethod(
  f = "kmeans",
  signature = signature(x = "SingleCellExperiment"),
  definition = function(x, reduceMethod = "PCA", whichAssay = 1, ...)
  {

    if(reduceMethod=="none"){
      if(NCOL(x)>10000)
        message("Note that you are running kmeans with more than 10,000 cells using all of the dimensions. You might consider running a dimensionality reduction step first.")
      fit <- kmeans(assays(x)[[whichAssay]], ...)
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
      fit <- kmeans(reducedDim(x, reduceMethod), ...)

    }

    return(fit)
  })

#' @rdname kmeans
#' @export
#' @importClassesFrom SingleCellExperiment LinearEmbeddingMatrix
#' @importFrom SingleCellExperiment sampleFactors
setMethod(
  f = "kmeans",
  signature = signature(x = "LinearEmbeddingMatrix"),
  definition = function(x, ...)
  {
    kmeans(sampleFactors(x), ...)
  })

#' @rdname kmeans
#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setMethod(
  f = "kmeans",
  signature = signature(x = "DelayedMatrix"),
  definition = function(x, ...)
  {
    stop("kmeans is not yet implemented for a DelayedMatrix/HDF5Matrix object")
  })

#' @rdname kmeans
#' @export
#' @importFrom stats kmeans
# Do we need this function, or will it be called automatically?
# See https://stat.ethz.ch/pipermail/r-devel/2009-March/052680.html and ?Methods_for_S3
setMethod(
  f = "kmeans",
  signature = signature(x = "matrix"),
  definition = function(x, ...)
  {
    stats::kmeans(x,...)
  })
