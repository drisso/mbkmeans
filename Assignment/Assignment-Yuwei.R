#Download the package
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("rhdf5")

#check the vignette
browseVignettes(package = "rhdf5")
setwd("/Users/summer/Desktop/Research Assistant")


#simulate the matrix
set.seed(123)
data<-matrix(rnorm(1000*100,mean = 0.5,sd=0.3),1000,100)


#create the HDF5 file and write the matrix
library(rhdf5)
h5createFile("assignment.h5")
h5write(data,"assignment.h5","matrix1")
h5<-h5read("assignment.h5","matrix1")


#check whether the results are same
identical(data,h5)


#download the HDF5array Package
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
library(HDF5Array)
library(tidyverse)


#save the matrix as an HDF5-based object
library(SummarizedExperiment)
matrix2<-SummarizedExperiment(assays = data)
dir<-sub("/Users/summer/Desktop/Research Assistant","hdf5matrix",tempfile())
hdf5matrix<-saveHDF5SummarizedExperiment(matrix2,dir)



#calculate the mean of each columns

#natural method
colMeans(data)

#hdf5 method
h5_matrix<-loadHDF5SummarizedExperiment(dir)
colMeans(assay(h5_matrix))

#compare the results
identical(colMeans(data),colMeans(assay(h5_matrix)))


#compare the speed
system.time(colMeans(data))
system.time(colMeans(assay(h5_matrix)))
