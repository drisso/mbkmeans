#' @title k-means for large single cell sequencing data
#'
#' @description This is a wrapper for stats::kmeans() for
#' large single cell sequencing data with the dimensionality
#' reduction results as input in in the reducedDim() slot.
#' @param reduceMethod Name of dimensionality reduction results
#' to use as input to k-means
#' @return k-means output
#' @name kmeans
#' @rdname kmeans
#' @export
setMethod(
  f = "kmeans",
  signature = signature(x = "SummarizedExperiment"),
  definition = function(x, ...){
    kmeans(as(x,"SingleCellExperiment"),...)

  })

#' @rdname kmeans
#' @export
setMethod(
  f = "kmeans",
  signature = signature(x = "SingleCellExperiment"),
  definition = function(x, reduceMethod = "PCA", which_assay=1,...)
  {
    if(is.null(reducedDimNames(x))){
      kmeans(assays(x)[[1]])
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
      fit = kmeans(reducedDim(x, reduceMethod), ...)

    }

    return(fit)
  })


#' @rdname kmeans
#' @export
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