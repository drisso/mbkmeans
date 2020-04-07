context("Test chunk dims.")

test_that("mini_batch gives the same results on different chunks", {
    library(SummarizedExperiment)
    data(iris)
    irismat <- as.matrix(iris[,1:4])

    sce <- SummarizedExperiment(assays = list(iris = t(irismat)))
    saveHDF5SummarizedExperiment(sce, dir = "tmp1", replace = TRUE)
    saveHDF5SummarizedExperiment(sce, dir = "tmp2", chunkdim = c(nrow(sce), 1),
                                 replace = TRUE)
    saveHDF5SummarizedExperiment(sce, dir = "tmp3", chunkdim = c(1, nrow(sce)),
                                 replace = TRUE)

    tmp1 <- loadHDF5SummarizedExperiment("tmp1")
    tmp2 <- loadHDF5SummarizedExperiment("tmp2")
    tmp3 <- loadHDF5SummarizedExperiment("tmp3")

    set.seed(13124)
    system.time(mbk1 <- mbkmeans(tmp1, clusters = 3))
    set.seed(13124)
    system.time(mbk2 <- mbkmeans(tmp2, clusters = 3))
    set.seed(13124)
    system.time(mbk3 <- mbkmeans(tmp3, clusters = 3))

    expect_equal(mbk1, mbk2)
    expect_equal(mbk1, mbk3)
})
