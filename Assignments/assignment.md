1. Download and install the `rhdf5` package from Bioconductor
2. Read the vignette
3. Simulate a 1000 * 100 random matrix (`rnorm`) 
4. Use `rhdf5` package to write the matrix into hdf5 files
5. Use the `rhdf5` package to read the file
6. Make sure 5 and 3 are same
7. Download the `HDF5Array` Package from Bioconductor
8. Create a `HDF5matrix` with the data from 3
9. Compute column means of 8 and 3, and compare the speed (which is faster)
(`System.time` how much time to compute and get the results) 
10. Create a `analysis` folder in `beachball` repo
11. Add a script for 3-9 in the folder
