#' Compute labels for mini-batch k-means
#'
#' Given a data matrix and a centroid matrix, it assigns each data point to the
#' closest centroid, using block processing.
#'
#' @param data a matrix-like object with features in row and samples in columns.
#' @param centroids a matrix with the coordinates of the centroids.
#' @param BPPARAM for parallel computations. See the `BiocParallel` package.
#' @param ... passed to `blockApply`.
#'
#' @return a vector of cluster labels for each observation.
#'
#' @importFrom DelayedArray blockApply colSums colAutoGrid
#' @importFrom BiocParallel SerialParam
#' @export
#'
#' @examples
#'
#' data(iris)
#' km <- mini_batch(as.matrix(iris[,1:4]), clusters = 3,
#'                  batch_size = 10, max_iters = 100)
#' predict_mini_batch_r(t(as.matrix(iris[,1:4])), km$centroids)
predict_mini_batch_r <- function(data, centroids,
                                 BPPARAM=BiocParallel::SerialParam(),
                                 ...) {
    unlist(blockApply(data, all_labels, centroids = centroids,
                      grid = colAutoGrid(data),
                      BPPARAM = BPPARAM, ...))
}

one_centroid <- function(x, data) {
    colSums((x - data)^2)
}

all_labels <- function(data, centroids) {
    ss <- apply(centroids, 1, one_centroid, data = data)
    apply(ss, 1, which.min)
}
