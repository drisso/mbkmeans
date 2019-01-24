context("Test mbkmeans function.")
set.seed(13124)

test_that("all mbkmeans methods give same result", {
    library(SummarizedExperiment)
    library(SingleCellExperiment)
    se <- SummarizedExperiment(matrix(rpois(600, lambda=5), nrow=10, ncol=60),
                               colData = data.frame(bio = gl(2, 30)))

    sce <- as(se, "SingleCellExperiment")

    ## for now the seed is set internally; remember to set the seed for each run if
    ## removing the internal seed
    m_se <- mbkmeans(se, clusters=3)
    m_sce <-  mbkmeans(sce, reduceMethod = NA, clusters=3)
    m_m <- mbkmeans(assay(sce), clusters=3)
    mb <- mini_batch(t(assay(sce)), clusters=3, batch_size = blocksize(se),
                     max_iters = 10, init_fraction = 0.25)

    expect_equal(m_se, m_sce)
    expect_equal(m_se, m_m)
    expect_equal(m_se, mb)
})
