#' @rdname mbkmeans
setGeneric(
    name = "mbkmeans",
    def = function(x, ...) standardGeneric("mbkmeans")
)

#' Cluster rows of a matrix
#'
#' Cluster rows of a matrix-like object with a variety of algorithms.
#'
#' @param x A numeric matrix-like object where rows represent observations and columns represent variables.
#' @param BLUSPARAM A \linkS4class{BlusterParam} object specifying the algorithm to use.
#' @param full Logical scalar indicating whether the full clustering statistics should be returned for each method.
#'
#' @return
#' By default, a factor of length equal to \code{nrow(x)} containing cluster assignments for each row of \code{x}.
#'
#' If \code{full=TRUE}, a list is returned containing \code{clusters}, a factor as described above;
#' and \code{objects}, an arbitrary object containing algorithm-specific statistics or intermediate objects.
#'
#' @details
#' This generic allows users to write agile code that can use a variety of clustering algorithms.
#' By simply changing \code{BLUSPARAM}, we can tune the clustering procedure in analysis workflows and package functions.
#'
#' @seealso
#' \linkS4class{MbkmeansParam} for some examples of values for \code{BLUSPARAM}.
#' 
#' @author Aaron Lun, Stephanie Hicks
#'
#' @examples
#' m <- matrix(runif(10000), ncol=10)
#'
#' clusterRows(m, MbkmeansParam(2))
#' @export
setGeneric("clusterRows", function(x, BLUSPARAM, full=FALSE) standardGeneric("clusterRows"))

setGeneric(".extras", function(x) standardGeneric(".extras"))