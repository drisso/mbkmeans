context("Test mbkmeans function.")
set.seed(13124)

test_that("all mbkmeans methods give same result", {
    library(SummarizedExperiment)
    library(SingleCellExperiment)
    se <- SummarizedExperiment(matrix(rpois(600, lambda=5), nrow=10, ncol=60),
                               colData = data.frame(bio = gl(2, 30)))

    sce <- as(se, "SingleCellExperiment")

    set.seed(1)
    m_se <- mbkmeans(se, clusters=3)
    set.seed(1)
    m_sce <-  mbkmeans(sce, reduceMethod = NA, clusters=3)
    set.seed(1)
    m_m <- mbkmeans(assay(sce), clusters=3)
    set.seed(1)
    mb <- mini_batch(t(assay(sce)), clusters=3, batch_size = ncol(sce),
                     max_iters = 100, init_fraction = 1)

    expect_equal(m_se, m_sce)
    expect_equal(m_se, m_m)
    expect_equal(m_se, mb)
    
    library(Matrix)
    library(HDF5Array)
    
    set.seed(123)
    m0 <- matrix(rnorm(100), ncol=10)
    m1 <- as(m0, "sparseMatrix")
    M2 <- as(m0, "HDF5Matrix")          # HDF5Matrix instance
    M3 <- M2 + 1 - 1 # DelayedMatrix instance
    
    set.seed(123)
    res0 <- mbkmeans(m0, clusters=3)

    set.seed(123)
    res1 <- mbkmeans(m1, clusters=3)
    
    set.seed(123)
    res2 <- mbkmeans(M2, clusters=3)

    set.seed(123)
    res3 <- mbkmeans(M3, clusters=3)
    
    expect_equal(res0, res1)
    expect_equal(res0, res2)
    expect_equal(res0, res3)
    
})

test_that("Length of cluster labels is the same as number of obs", {
    library(HDF5Array)
    set.seed(123)
    m0 <- matrix(rnorm(2.3e7), ncol=1000)
    M2 <- as(m0, "HDF5Matrix")
    
    k <- 15
    batch <- 500
    
    res_1kcell <- mbkmeans(x=M2,
                           clusters=k, 
                           batch_size = batch, 
                           num_init=1, 
                           max_iters=100,
                           verbose=FALSE,
                           compute_labels=TRUE)
    
    expect_equal(length(res_1kcell$Clusters), NCOL(M2))
    
    fit <- mini_batch(data = t(M2), clusters = k,
                      batch_size = batch, max_iters = 100,
                      num_init = 1,
                      init_fraction = batch/NCOL(M2),
                      initializer = "kmeans++",
                      compute_labels = FALSE,
                      calc_wcss = FALSE,
                      early_stop_iter = 10,
                      verbose = FALSE,
                      CENTROIDS = NULL, tol = 1e-4)
    
    # the following creates 2k labels instead of 1k
    pred_r_1k <- predict_mini_batch_r(data = M2,
                                      centroids = fit$centroids)
    
    expect_equal(length(pred_r_1k), NCOL(M2))
    
    # the following works fine
    pred_c_1k <- predict_mini_batch(data = t(M2),
                                      CENTROIDS = fit$centroids)
    
    expect_equal(length(pred_c_1k), NCOL(M2))
    
    expect_equal(pred_r_1k, pred_c_1k)
    
})
