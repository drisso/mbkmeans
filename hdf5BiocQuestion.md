We are working at implementing hdf5 compatible methods in our package, and have some questions about the best practice. 

Here's some start up code that I will use in the following.

```
library(HDF5Array)
library(SummarizedExperiment)
nrows <- 200; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
se0 <- SummarizedExperiment(assays=SimpleList(counts=counts))
dir <- sub("file", "h5_se0_", tempfile())
h5_se0 <- saveHDF5SummarizedExperiment(se0, dir)
```

# Question 1: 

*After writing this, I think this question falls apart into two items 1) Asking Pete why he was opposed to class unions and 2) is there a virtual class for matrices and these objects*

From an interactive point of view, there are many functions that work for users the same for both HDF5Matrix and standard matrices, e.g. rowMeans. From a package developer point of view, however, we do not know the best practices for how to handle this. 

Consider the following simple S3 function. 

```
myNewRowMeans<-function(x,...){
	#a lot of code independent of class of x
	print("This is a lot of code shared regardless of class of x\n")
	out<-rowMeans(x)
	#a lot of code on output of out
	out<-out+1
	return(out)	
}
```

```
testMat<-myNewRowMeans(counts)
testDA<-myNewRowMeans(assay(h5_se0))
```

From an interactive point of view, this function is fine because it doesn't care what class x is. If I want to put it into my package I am not sure if I need to do anything to address the fact that rowMeans for the two classes comes from different packages. I think I should just be able to add

```
#' @importFrom DelayedArray rowMeans
```

Is that sufficient, or is there better practice? (often I use '::' to specificy the package clearly, but that's clearly not possible here).

Furthermore to make this an S4 function, I have options in how I set this up. My preferred choice would be to create a class union -- **Does such a virtual class exist already?**

For example,

```
setGeneric("myNewRowMeans", function(x,...) { standardGeneric("myNewRowMeans")})
setClassUnion("matrixOrDelayed",members=c("matrix", "DelayedArray"))

setMethod("myNewRowMeans", 
          signature = "matrixOrDelayed",
          definition = function(x,...){
		  	print("This is a lot of code shared regardless of class of x\n")
		  	out<-rowMeans(x)
		  	#a lot of code on output of out
		  	out<-out+1
		  	return(out)	
		}
		)
```

However, in a casual conversation with Pete, he seemed negative on this idea (not sure why -- felt not very robust?)

The alternative is to make two different methods, one for `DelayedArray` and another for `matrix`. This requires either copying of shared code (!) or make the shared code an internal function that both call, which is rather annoying if its not very long code. This seems clunky and defeats the point of S4 structure, it seems. 

# Question 2:

*I think this is a real question to be asked*

When the assay slot of a SummarizedExperiment has a HDF5Matrix, it doesn't appear with the class `HDF5Matrix` unless you specify to not have the dimnames. Namely the class of `assay(h5_se0)` varies:

```
> class(assay(h5_se0))
[1] “DelayedMatrix”
attr(,“package”)
[1] “DelayedArray”
> class(assay(h5_se0,withDimnames=FALSE))
[1] “HDF5Matrix”
attr(,“package”)
[1] “HDF5Array”
```
[Note: based on my memory, this is a change, since I thought previously `assay(h5_se0)` was a `HDF5Matrix`].


However, testing the class with `is` always works for a `DelayedMatrix` (i.e. always shows up as a `DelayedMatrix`)

```
> is(assay(h5_se0), "HDF5Matrix")
[1] FALSE
> is(assay(h5_se0,withDimnames=FALSE),"HDF5Matrix")
[1] TRUE
> is(assay(h5_se0), "DelayedMatrix")
[1] TRUE
> is(assay(h5_se0,withDimnames=FALSE),"DelayedMatrix")
[1] TRUE
```


The question is whether in creating methods for these types of objects, is the best practice to use `DelayedArray` (since they both inherit from this class, while `HDF5Matrix` doesn't inherit from `DelayedMatrix`)? The downside is that there might be objects that would be arrays but not matrices that might break some function -- is this a realistic concern?

The other option is to always use `DelayedMatrix` since this class always seems to be recognized, but I don't know if the `setMethod` would recognize this. Also, it makes the function less general if calling directly on `HDF5Matrix`, and not on the `assay` of a `SummarizedExperiment`.

