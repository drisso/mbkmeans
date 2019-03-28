context("Test mbkmeans function.")
set.seed(13124)

test_that("mini_batch gives the same results on the iris data", {
    data(iris)
    irismat <- as.matrix(iris[,1:4])

    irisres <- list(centroids = matrix(data = c(6.95, 4.983333, 5.8, 3.45,
                                                3.383333, 2.65, 5.9, 1.433333,
                                                4.25, 2.3, 0.2, 1.25),
                                       nrow = 3, ncol=4),
                    WCSS_per_cluster = c(30.235, 15.42333, 46.3725),
                    best_initialization = 1,
                    iters_per_initialization = matrix(10),
                    Clusters =   c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 1, 3, 3, 3, 3,
                                   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                                   3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                                   3, 3, 3, 3, 3, 1, 3, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 3,
                                   3, 1, 1, 1, 1, 3, 1, 3, 1, 3, 1, 1, 3, 3, 1, 1, 1, 1, 1,
                                   3, 3, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 3))

    set.seed(1)
    mb <- mini_batch(irismat, clusters=3, batch_size = 10,
                     max_iters = 10, init_fraction = 0.25)

    expect_true(all(mb$centroids - irisres$centroids <1e4))
    expect_true(all(mb$Clusters - irisres$Clusters <1e4))
    expect_equal(mb$best_initialization, irisres$best_initialization)
    expect_equal(mb$iters_per_initialization, irisres$iters_per_initialization)
})

test_that("WCSS calculation is correct", {
    data(iris)
    irismat <- as.matrix(iris[,1:4])

    mb <- mini_batch(irismat, clusters=3, batch_size = 10,
                     max_iters = 100, init_fraction = 0.25, calc_wcss = TRUE)
    expect_equal(compute_wcss(mb$Clusters, mb$centroids, irismat), mb$WCSS_per_cluster)

    km <- kmeans(irismat, 3)
    expect_equal(compute_wcss(km$cluster, km$centers, irismat), km$withinss)
})

test_that("clustering is accurate", {
    data(iris)
    irismat <- as.matrix(iris[,1:4])

    ## when starting from kmeans solution it should not move

    km <- kmeans(irismat, 3)

    mb <- mini_batch(irismat, clusters=3, batch_size = nrow(irismat),
                     max_iters = 100, init_fraction = 0.25,
                     CENTROIDS = km$centers)
    expect_equal(mb$Clusters, km$cluster)

})

test_that("all mini_batch methods give same result", {
  library(Matrix)
  library(HDF5Array)

  set.seed(123)
  m0 <- matrix(rnorm(100), ncol=10)
  m1 <- as(m0, "sparseMatrix")
  M2 <- as(m0, "HDF5Matrix")          # HDF5Matrix instance
  M3 <- rbind(M2[1:5, ], M2[6:10, ])  # DelayedMatrix instance
  
  set.seed(123)
  res0 <- mini_batch(m0, clusters=3, batch_size=10, max_iters = 10)
  
  set.seed(123)
  res1 <- mini_batch(m1, clusters=3, batch_size=10, max_iters = 10)
  
  set.seed(123)
  res2 <- mini_batch(M2, clusters=3, batch_size=10, max_iters = 10)

  set.seed(123)
  res3 <- mini_batch(M3, clusters=3, batch_size=10, max_iters = 10)

  expect_equal(res0, res1)
  expect_equal(res0, res2)
  expect_equal(res0, res3)
  
})
