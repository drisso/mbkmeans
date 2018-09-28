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
#' @name mbkmeans
#' @rdname mbkmeans
#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods as
#' @examples
#' se <- SummarizedExperiment(matrix(rnorm(100), ncol=10))
#' kmeans(se, centers = 2, reduceMethod = "none")
setMethod(
  f = "mbkmeans",
  signature = signature(x = "SummarizedExperiment"),
  definition = function(x, ...){
    mbkmeans(as(x,"SingleCellExperiment"),...)

  })

#' @rdname kmeans
#' @export
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
setMethod(
  f = "mbkmeans",
  signature = signature(x = "SingleCellExperiment"),
  definition = function(x, reduceMethod = "PCA", whichAssay = 1, ...)
  {

    if(reduceMethod=="none"){
      if(NCOL(x)>10000)
        message("Note that you are running kmeans with more than 10,000 cells using all of the dimensions. You might consider running a dimensionality reduction step first.")
      fit <- mbkmeans(assays(x)[[whichAssay]], ...)
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
      fit <- mbkmeans(reducedDim(x, reduceMethod), ...)

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
    mbkmeans(sampleFactors(x), ...)
  })

#' @rdname mbkmeans
#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setMethod(
  f = "mbkmeans",
  signature = signature(x = "ANY"),
  definition = function(x ) #put arguments here
  {
	  mini_batch(x)#put arguments here
  })

