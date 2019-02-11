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
