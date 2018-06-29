#Week 2
#Yuwei Ni
#install all required packaegs
library(SingleCellExperiment)
library(restfulSE)
library(DelayedArray)
library(TENxBrainData)

#import the 1.3 million neurons from 10X Genomics data
my10x = se1.3M()

#transfer the type from SummarizedExperiment to SingleCellExperiment
data<-as(my10x,"SingleCellExperiment")

#choose subset of data
data2 <- data[, seq_len(2000)]
chunksize <- 500
cidx <- snow::splitIndices(ncol(data2), ncol(data2) / chunksize)

lib.sizes <- n.exprs <- numeric(ncol(data2))
tot.exprs <- numeric(nrow(data2))
for (i in head(cidx, 2)) {
  message(".", appendLF=FALSE)
  m <- as.matrix(counts(data2)[,i, drop=FALSE])
  lib.sizes[i] <- colSums(m)
  n.exprs[i] <- colSums(m != 0)
  tot.exprs <- tot.exprs + rowSums(m)
}
ave.exprs <- tot.exprs / ncol(data2)

#calculate and save the library size, number of express genes and the average
colData(data2)$lib.sizes <- lib.sizes
colData(data2)$n.exprs <- n.exprs
rowData(data2)$ave.exprs <- ave.exprs

#filter the rows to remove the 0 counts
feature<- rowSums(counts(data2)>0)>0
dataclean2<-data2[feature,]

#counts for cell to be used a scaling
library(scater)
#dt<-colData(data2)
#dataclean<- dt[which(dt[,"lib.sizes"]>0),]
sizeFactors(data2)
nor<-normalize(data2)
assayNames(nor)


#Cluster
library(ClusterR)



#Appendix

#Example of filtering the rows
y<-matrix(1:4,nrow = 2,ncol = 2)
y1<-SingleCellExperiment(assays = list(counts=y))
feature1<-rowSums(counts(y1)>3)>0
y2<-y1[feature1,]
counts(y2)

#Question:
#1.Q about DropletUtils, see notes
#2.How to filter the rows to remova the genes without expression? Why the code above cannout run?
#3.How to calculate the total number of counts for each cell to be used a scaling factor for normalization ?
#4.