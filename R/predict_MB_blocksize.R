#' predict_MB_blocksize
#'
#' Predict Mini-batch-k-means by block for both matrix and HDF5Matrix
#'
#' @param data   numeric matrix or integer matrix or HDF5Matrix
#' @param block_size   the number of rows of data
#' @param clusters  the number of clusters
#' @param batch_size the size of the mini batches
#' @param init_fraction percentage of data to use for the initialization centroids (applies if initializer is \emph{kmeans++} ). Should be a float number between 0.0 and 1.0.
#' @param max_iters the maximum number of clustering iterations
#' @param initializer the method of initialization. One of \emph{kmeans++} and \emph{random}
#'
#' @return  It returns a vector with the clusters
#'
#' @section
#' blocksize: \code{\link[mbkmeans]{blocksize}}
#'
#' @importFrom DelayedArray blockApply RegularArrayGrid
#' @export

predict_MB_blocksize<-function(data,block_size,clusters,batch_size,init_fraction,max_iters,initializer = "kmeans++"){

  if(missing(block_size)) {
    #min(floor(as.numeric(get_ram())/(2*8*ncol(data))), nrow(data))
    block_size = blocksize(data)
  }

  km<-mini_batch(data, clusters = clusters, batch_size = batch_size, init_fraction = init_fraction, max_iters = max_iters,initializer =initializer )

  ##for each row, do predict_mini_batch
    km_block <- blockApply(
    x = data,
    ##FUN = predict_MBatchKMeans,
    FUN = predict_mini_batch,
    CENTROIDS = km$centroids,
    grid = RegularArrayGrid(
      refdim = dim(data),
      spacings = c(block_size, ncol(data))))

    clean_data<-unlist(km_block,use.names = FALSE)

    clean_data
}
