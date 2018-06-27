# Kmeans ++

Picks different set of initialization points by iteratively picking next cluster center so that far from the others. 

$X$ set of points in $R^d$. Let $d(x,Y)=\min_{y\in Y} ||x-y||$ (distance of point $x$ to set of points $Y$). 
$C$ is a set of points (usually centroids), and $\phi_Y(C)=\sum_{y\in Y} d^2(y,C)$. 

*Algorithm Initialization of Kmeans++*
1. Pick $c_1$ uniformly from $X$ and set $C=\{c_1\}$
2. while $|C| < k$ do
3.      Sample a single $c_2\in X$, where each $x$ has probability $d^2(x,C)/ \phi_X(C)$ of being selected
4.      Add $c_2$ to $C$.
5. end while.

Then starting with these as the centers of clustesr, do Lloyd iterations

1. Given cluster centers $C$, assign points to cluster whose center is closest
2. Given assignment of points to center, calculate new centers as average of points

## kmeans++||
Total run time of kmeans++ is $O(nkd)$  which is same as Lloyd's algorithm. However, because this initialization is iterative over data, can't parallelize this initialization step. Need $k$ passes (iterations) over the data. 

kmeans++|| instead samples $O(k)$ points in each pass and only $O(logn)$ iterations. So approximately will have fewer iterations if $\logn < k$. For $n=$1million, this is if $k>14$, 10million if $k>16$ etc. 

*Algorithm Initialization of Kmeans++||*

0. Set $\ell=\Omega(k)$ 
1. Pick $c_1$ uniformly from $X$ and set $C=\{c_1\}$
2. Set $r=O(log(\phi_X(C)))$ (the number of iterations to do)
3. for $r$ rounds, do
4.      For each $x\inX$, add $x$ to $C$ with probability $\ell d^2(x,C)/\phi_X(C)$
5. end for
6. For each $c\inC$, assign points to cluster whose center is closest and set weight $w_c$ equal to number of points assigned to cluster $c$
7. Recluster the *weighted* points in $C$ into $k$ clusters (e.g. with kmeans++)
8. Return the centers of these $k$ clusters as initial cluster centers 

Note that step 4 can be done in parallel -- i.e. can calculate probability of point and decide to sample it independently of other points. Also, the size of $C$ after the end of the for-loop is expected to be $>k$: $E(|C|)=r\ell$ which guides the choices of the parameters

They used $\ell=0.1k,0.5k,k,2k,10k$ in their simulations and r=15 for $\ell=0.1k$ and $r=5$ otherwise. 

