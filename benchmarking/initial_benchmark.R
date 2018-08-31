library(HDF5Array)
library(pryr)
devtools::load_all()

mat <- HDF5Array("benchmarking/pbmc3k_rectangular.h5", name = "counts")

data <- t(mat)[,1:1000]
system.time(km <- mini_batch(data, clusters = 3, batch_size = 100, init_fraction = .1, max_iters = 10))
system.time(clusters <- predict_mini_batch(data, CENTROIDS = km$centroids))

data2 <- as.matrix(data)
system.time(km <- MiniBatchKmeans(data2, clusters = 3, batch_size = 100, init_fraction = .1, max_iters = 10))
system.time(km <- mini_batch(data2, clusters = 3, batch_size = 100, init_fraction = .1, max_iters = 10))
system.time(clusters <- predict_MBatchKMeans(data2, CENTROIDS = km$centroids))
system.time(clusters <- predict_mini_batch(data2, CENTROIDS = km$centroids))

system.time(km <- stats::kmeans(data, centers = 3))

object_size(data)
object_size(data2)

library(DelayedArray)
DelayedArray:::set_verbose_block_processing(FALSE)
block_size <- 1000L
system.time(
km_block <- blockApply(
  x = data,
  FUN = predict_MBatchKMeans,
  CENTROIDS = km$centroids,
  grid = RegularArrayGrid(
    refdim = dim(data),
    spacings = c(block_size, ncol(data))))
)
table(clusters[1,], unlist(km_block))

