
#check the vignette
#rhdf5
#browseVignettes(package = "rhdf5")

#HDF5Array
#https://www.bioconductor.org/packages/release/bioc/html/HDF5Array.html

#DelayedArray--Development of HDF5Array
#https://www.bioconductor.org/packages/release/bioc/vignettes/DelayedArray/inst/doc/01-Working_with_large_arrays.pdf

#SummarizedExperiment
#browseVignettes("SummarizedExperiment")

#SingleCellExperiment
#https://www.bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html


#Generate hdf5 with random data
set.seed(123)
data<-matrix(rnorm(1000*100,mean = 0.5,sd=0.3),1000,100)

#generate hdf5 file with random data
library(rhdf5)
h5createFile("week1.h5")
h5write(data,"week1.h5","matrix1")
h5<-h5read("week1.h5","matrix1")

#calculate the colmean
colMeans(h5)

#SummarizedExperiment with hdf5 object
library(SummarizedExperiment)
matrix1<-SummarizedExperiment(assays = h5)
dir<-sub(getwd(),"hdf5matrix",tempfile())
sematrix<-saveHDF5SummarizedExperiment(matrix1,dir)


#SingleCellExperiment with hdf5 object
library(SingleCellExperiment)
scematrix<-SingleCellExperiment(assays = h5)


#kmeans
sematrix2<-as.matrix(assay(sematrix))
km<-kmeans(sematrix2,5)  ##Kmeans groups data based on the rows
km$centers

#kmean.plus
library(LICORS)
kmpp<-kmeanspp(sematrix2,5)  ##kmean.plus change the way to choose the centers
kmpp$centers

#Question
#1.How to create small shell script that runs R script?
#2.SingleCellExperiment cannot be installed