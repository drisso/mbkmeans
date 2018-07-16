**Goal**: to modify the `mini_batch_kmeans` function to use the `beachmat` library so that it works with regular matrices as well as HDF5 matrices.

1. Read all the vignettes of the `beachmat` package: https://bioconductor.org/packages/release/bioc/html/beachmat.html
2. Modify the code in `Week3.cpp` to change the current arma matrix into a beachmat matrix -- it should work with both numeric and integer
3. Make sure that the results are the same
4. Write unit test that verifies that we get the same results with a regular matrix and with a HDF5 matrix.

