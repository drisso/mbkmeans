context("Test blocksize function.")
set.seed(13124)

test_that("blocksize works with and without ram argument", {
    data(iris)
    irismat <- as.matrix(iris[,1:4])

    expect_silent(b1 <- blocksize(irismat, ram = 1e9))
    expect_warning(b2 <- blocksize(irismat, ram = NA))

    expect_equal(b1, b2)
})

