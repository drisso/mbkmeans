#' Mini-batch k-means clustering
#'
#' Run the mini-batch k-means \code{\link{mbkmeans}} function with the specified number of centers within \code{\link{clusterRows}}.
#'
#' @param centers An integer scalar specifying the number of centers.
#' Alternatively, a function that takes the number of observations and returns the number of centers.
#' @param ... Further arguments to pass to \code{\link{mbkmeans}}.
#' @inheritParams clusterRows
#' @param BLUSPARAM A \linkS4class{MbkmeansParam} object.
#' @param full Logical scalar indicating whether the full mini-batch k-means statistics should be returned.
#'
#' @author Stephanie Hicks
#'
#' @details
#' This class usually requires the user to specify the number of clusters beforehand.
#' However, we can also allow the number of clusters to vary as a function of the number of observations.
#' The latter is occasionally useful, e.g., to allow the clustering to automatically become more granular for large datasets.
#'
#' To modify an existing MbkmeansParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#'
#' @return 
#' The \code{MbkmeansParam} constructor will return a \linkS4class{MbkmeansParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects};
#' the latter will contain the direct output of \code{\link{mbkmeans}}.
#'
#' @examples
#' clusterRows(iris[,1:4], MbkmeansParam(centers=3))
#' clusterRows(iris[,1:4], MbkmeansParam(centers=3, batch_size=10))
#' clusterRows(iris[,1:4], MbkmeansParam(centers=3, batch_size=10, 
#'             compute_labels=TRUE, calc_wcss=TRUE))
#' @seealso
#' \code{\link{mbkmeans}}, which actually does all the heavy lifting.
#' @name MbkmeansParam-class
#' @docType class
#' @aliases 
#' show,MbkmeansParam-method
NULL

#' @export
#' @rdname MbkmeansParam-class
MbkmeansParam <- function(centers, ...) {
  if (!is.function(centers)) {
    centers <- as.integer(centers)
  }
  new("MbkmeansParam", centers=centers, extra.args=list(...))
}

setMethod(".extras", "MbkmeansParam", function(x) "extra.args")

#' @export
#' @importFrom S4Vectors coolcat
setMethod("show", "MbkmeansParam", function(object) {
  callNextMethod()
  cat(sprintf("centers: %s\n", if (is.function(object@centers)) "variable" else object@centers))
  coolcat("extra.args(%i): %s", names(object@extra.args))
})

#' @export
#' @rdname MbkmeansParam-class
setMethod("clusterRows", c("ANY", "MbkmeansParam"), function(x, BLUSPARAM, full=FALSE) {
  centers <- BLUSPARAM@centers
  if (is.function(centers)) {
    centers <- centers(nrow(x))
  }
  
  args <- c(list(x=as.matrix(t(x)), clusters=centers), BLUSPARAM@extra.args)
  stats <- do.call(mbkmeans, args)
  vec_clusters <- factor(stats$Clusters)
  
  if (full) {
    list(clusters=vec_clusters, objects=stats)
  } else {
    vec_clusters
  }
})
