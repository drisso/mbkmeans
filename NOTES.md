
# vanilla _k_-means

## Methods 

1. **vanilla _k_-means**: pick _k_, random sample of _k_ "centers", and minimize the average squared distance between points to the assigned cluster centers. At the minimum, all cluster centers are at the mean of their Voronoi sets (the set of data points which are nearest to the cluster center). 
	- **Good:** simple, fast
	- *Bad:** No accuracy guarantees; falls into local minima so several restarts are useful. 
2. **_k_-means++** ][reference](http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf)]: _k_-means with randomized seeding for initially picking centers
	- **Good:** faster, more accurate than _k_-means
	- **Bad:** 
	
## Algorithms

1. **Hartigan-Wong Algorithm** [[reference](https://www.jstor.org/stable/2346830)]: Default used in the `kmeans()` base R function
2. **Lloyd's Algorithm** [[reference](http://www-evasion.imag.fr/people/Franck.Hetroy/Teaching/ProjetsImage/2007/Bib/lloyd-1982.pdf), [Wikipedia](https://en.wikipedia.org/wiki/Lloyd%27s_algorithm)]:  Implemented in the `kmeans()` base R function
3. **MacQueen Algorithm** [[reference](https://projecteuclid.org/euclid.bsmsp/1200512992)]: a serial _k_-means algorithm; Implemented in the `kmeans()` base R function

 
# What to do if data is too big to be read into memory? 

## Try running _k_-means on all data using HDF5 files

- [`bigkmeans`](https://github.com/argriffing/bigkmeans) -- Python code on GitHub implementing Lloyd's algorithm for _k_-means that allows for HDF5 files
- [scikit-learn](http://scikit-learn.org/stable/datasets/index.html) -- Says "if you manage your own numerical data it is recommended to use an optimized file format such as HDF5 to reduce data load times. Various libraries such as [H5Py](https://www.h5py.org) (access multi-terabyte datasets stored on disk as if they were real NumPy arrays), PyTables and pandas provides a Python interface for reading and writing data in that format."
	- [sklearn.cluster.KMeans()](http://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html): default initialization is _k_-means++; solved using Lloyd's algorithm on entire dataset. 
	- [sklearn.cluster.MiniBatchKMeans()](http://scikit-learn.org/stable/modules/generated/sklearn.cluster.MiniBatchKMeans.html#sklearn.cluster.MiniBatchKMeans): default initialization is _k_-means++; online implementation that does incremental updates of the centers using mini-batches [see notes below on subsampling and [Sculley 2010](http://www.eecs.tufts.edu/~dsculley/papers/fastkmeans.pdf)]

#### Relevant notes on HDF5 files

- [Mike Smith's blog post on parallel processing and HDF5](http://www.msmith.de/2018/05/01/parallel-r-hdf5/)
- [Pete's blog post on working with DelayedArrays (more to come)](https://www.peterhickey.org/2018/05/01/bioc3.7-and-delayedarray/)
- [k-NN in python and HDF5](https://labrosa.ee.columbia.edu/millionsong/pages/fast-k-nn-using-hdf5)


## Consider alternative approaches to speed up _k_-means method

1. Parallelization
	- kmeans in [spark](http://spark.apache.org/docs/latest/mllib-clustering.html#k-means) draws on [kmeans||](http://theory.stanford.edu/~sergei/papers/vldb12-kmpar.pdf) -- a parallel implementation of kmeans++.
	- http://www.ece.northwestern.edu/~wkliao/Kmeans/index.html : gives C code for parallel kmeans
	- Algorithms implementing _k_-means with MapReduce (e.g. Hadoop) -- useful for handling a large volume of data over a distributed computing environment.
		- [Zhao et al. (2009)](https://link.springer.com/chapter/10.1007/978-3-642-10665-1_71): Reference for PKMeans (Parallel K Means)
			- **Bad:** algorithm has an unlimited number of mappers, but at most 1 reducer for each centroid, making it less efficient
		- [Kerdprasop and Kerdprasop 2010](https://pdfs.semanticscholar.org/a76a/f136805e8f535777ad5582b128ae441af75a.pdf): _k_-means clustering on multi-core processors
		- [Anchalia (2014)](https://ieeexplore.ieee.org/document/7046097/): Improves upon Zhao et al. (2009) by noticing the overhead is mostly from the shuffling stage, so they introduce "a combiner in the mapper function to decrease the amount of data to be written by the mapper and the amount of data to be read by the reducer"
		- [Gursoy (2003)](https://link.springer.com/chapter/10.1007/978-3-540-24669-5_31): proposes a _k_-means method using a k-d tree to optimize search for nearest centroids
		- [Jin, Cui, and Yu (2016)](https://arxiv.org/pdf/1608.06347.pdf) (link includes pseudo code implementing _k_-means with MapReduce): Reference for IPKMeans (Improved Parallel K Means). Instead of using k-d tree to optimize for the nearest centroid, uses k-d tree to partition datasets into smaller representative datasets (of the larger dataset). Once the larger dataset has been reduced to smaller representative datasets, each reducer runs one complete _k_-means for each subset, generates 1 group of intermediate centroids, and all intermediate centroids are merged into a final group of centroids. 
			- Does a good job of talking about the high I/O overhead as number of data points increases
			- In contrast to other algorithms above which all need _M_ MapReduce jobs sequentially for _M_ _k_-means iterations, this algorithm only needs 1 MapReduce job (lower I/O overhead compared to PKMeans)
			- Considered different approaches to merge intermediate centroids: hierarchical merging or minimum average sum of square error (ASSE) ? 
			- **Bad:** As one reducer only executes _k_-means on a small portion, can lead to a loss of accuracy
2. Subsampling
	- MiniBatchKMeans seems popular. Works on random subsamples. Exists implementation in R (ClusterR). Here's some links:
		* https://algorithmicthoughts.wordpress.com/2013/07/26/machine-learning-mini-batch-k-means/
		* http://scikit-learn.org/stable/modules/clustering.html#mini-batch-kmeans
		* https://www.eecs.tufts.edu/~dsculley/papers/fastkmeans.pdf
		* http://theory.stanford.edu/~sergei/papers/vldb12-kmpar.pdf
		* https://cran.r-project.org/web/packages/ClusterR/index.html

## Other notes about _k_-means

- [How slow is the _k_-means method?](http://theory.stanford.edu/~sergei/papers/kMeans-socg.pdf)




# Approximations to _k_-means

## Methods

1.  **Method name** [reference]():
	- **Idea:** 
	- **Assumptions:**
	- **Good:** 	
	- **Bad:** 




