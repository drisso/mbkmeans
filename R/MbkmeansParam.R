#' Mini-batch k-means clustering
#'
#' Run the mini-batch k-means \code{\link{mbkmeans}} function with the specified
#' number of centers within \code{\link{clusterRows}} from the
#' \code{\link{bluster}} Bioconductor package.
#'
#' This function is deprecated. Please use the \code{MbkmeansParam} function in
#' the \code{\link{bluster}} Bioconductor package.
#' @name MbkmeansParam-class
#' @docType class
#' @aliases 
#' show,MbkmeansParam-method
NULL

#' @export
#' @rdname MbkmeansParam-class
MbkmeansParam <- function(centers, ...) {
  .Deprecated("bluster::MbkmeansParam")
  bluster::MbkmeansParam(centers, ...)
}
