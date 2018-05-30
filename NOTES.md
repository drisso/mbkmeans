
# vanilla _k_-means

## Methods 

1. **vanilla _k_-means**: pick _k_, random sample of _k_ "centers", and minimize the average squared distance between points in the same cluster
	- **Good:** simple, fast
	- *Bad:** No accuracy guarantees
2. **_k_-means++** [reference](http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf): _k_-means with randomized seeding for initially picking centers
	- **Good:** faster, more accurate than _k_-means
	- **Bad:** 


# Approximations to _k_-means

## Methods

1.  **Method name** [reference]():
	- **Idea:** 
	- **Assumptions:**
	- **Good:** 	
	- **Bad:** 


# Notes on R and HDF5 files

- [Mike Smith's blog post on parallel processing and HDF5](http://www.msmith.de/2018/05/01/parallel-r-hdf5/)
- [Pete's blog post on working with DelayedArrays (more to come)](https://www.peterhickey.org/2018/05/01/bioc3.7-and-delayedarray/)

# Notes Python and HDF5 files

- [k-NN in python and HDF5](https://labrosa.ee.columbia.edu/millionsong/pages/fast-k-nn-using-hdf5)


# kmeans in parallel
- kmeans in [spark](http://spark.apache.org/docs/latest/mllib-clustering.html#k-means) draws on [kmeans||](http://theory.stanford.edu/~sergei/papers/vldb12-kmpar.pdf) -- a parallel implementation of kmeans++.
- http://www.ece.northwestern.edu/~wkliao/Kmeans/index.html : gives C code for parallel kmeans

# fast kmeans
- MiniBatchKMeans seems popular. Works on random subsamples. Exists implementation in R (ClusterR). Here's some links:
	* https://algorithmicthoughts.wordpress.com/2013/07/26/machine-learning-mini-batch-k-means/
	* http://scikit-learn.org/stable/modules/clustering.html#mini-batch-kmeans
	* https://www.eecs.tufts.edu/~dsculley/papers/fastkmeans.pdf
	* http://theory.stanford.edu/~sergei/papers/vldb12-kmpar.pdf
	* https://cran.r-project.org/web/packages/ClusterR/index.html


