#' Mini-batch k-means clustering
#'
#' Run the mini-batch k-means \code{\link{mbkmeans}} function with the specified
#' number of centers within \code{clusterRows} from the
#' \code{bluster} Bioconductor package.
#'
#' This function is deprecated. Please use the \code{MbkmeansParam} function in
#' the \code{bluster} Bioconductor package.
#' 
#' @param centers An integer scalar specifying the number of centers.
#' Alternatively, a function that takes the number of observations and returns the number of centers. 
#' Note, the \code{\link{mbkmeans}} function uses the argument \code{clusters} argument to represent this argument. 
#' However, we use \code{centers} to match 
#' @param ... Further arguments to pass to \code{\link{mbkmeans}}.
#' 
#' @name MbkmeansParam
#' @docType class
#' @aliases 
#' show,MbkmeansParam-method
NULL

#' @export
#' @rdname MbkmeansParam
MbkmeansParam <- function(centers, ...) {
  .Deprecated("bluster::MbkmeansParam")
  bluster::MbkmeansParam(centers, ...)
}
