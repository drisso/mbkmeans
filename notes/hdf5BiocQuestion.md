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

Subject: Virtual class for  `matrix` and `DelayedArray`? (or better strategy for dealing with them both)

I am trying to extend my package to handle `HDF5Matrix` class ( or more generally `DelayedArray`). I currently have S4 functions for `matrix` class. Usually I have a method for `SummarizedExperiment`, which will call call the method on `assay(x)` and I want the method to be able to deal with if `assay(x)` is a `DelayedArray`.

Most of my functions, however, do not require separate code depending on whether `x` is a `matrix` or `DelayedArray`. They are making use of existing functions that will make that choice for me, e.g. rowMeans or subsetting. My goal right now is compatibility, not cleverness, and I'm not creating HDF5 methods to handle other cases. (If something doesn't currently exist, then I just enclose `x` with `data.matrix` or `as.matrix` and call the matrix into memory, though for cleanliness and ease in updating with appropriate methods, I could make separate S4 functions for these specific tasks to dispatch, but that's outside of the scope of my question.) So for simplicity assume I don't really need to dispatch *my code* -- that the methods I'm going to use do that. 

The natural solution for me seem to use `setClassUnion` and I was wondering if such a virtual class already exists? Or is there a better way to handle this?

Here's a simple example, using `rowMeans` as my example:

```
setGeneric("myNewRowMeans", function(x,...) { standardGeneric("myNewRowMeans")})
setClassUnion("matrixOrDelayed",members=c("matrix", "DelayedArray"))

#' @importFrom DelayedArray rowMeans
setMethod("myNewRowMeans", 
          signature = "matrixOrDelayed",
          definition = function(x,...){
		  	# code independent of x
			print("This is a lot of code shared regardless of class of x\n")
			# code that depends on x, but is dispatched by other functions
		  	out<-rowMeans(x)
		  	#a lot of code based on output of out
		  	out<-out+1
		  	return(out)	
		}
)
```



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


The question is whether in creating methods for these types of objects, is the best practice to use `DelayedArray` (since they both inherit from this class, while `HDF5Matrix` doesn't inherit from `DelayedMatrix`)? The conceptual downside is that there might be objects that would be `DelayedArray` but not either of these classes that might break some function -- is this a realistic concern or is the role of `DelayedArray` to be the a virtual class for both?

Another option is to always use `DelayedMatrix` since this class always seems to be recognized, but I don't know if the methods dispatch  would recognize this [I guess I should test it?]. Also, it makes the function less general since it wouldn't work if calling directly on a `HDF5Matrix` object, and not on the `assay` of a `SummarizedExperiment`.

And the last option is to make a class union of `DelayedMatrix` and `HDF5Matrix` (and again, is there already a virtual class for these?).
