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

