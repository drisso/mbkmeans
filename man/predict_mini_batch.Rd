% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{predict_mini_batch}
\alias{predict_mini_batch}
\title{Predict_mini_batch}
\usage{
predict_mini_batch(data, CENTROIDS)
}
\arguments{
\item{data}{matrix-like objectcontaining numeric or
integer data (obseravtions in rows, variables in columns).}

\item{CENTROIDS}{a matrix of initial cluster centroids. The rows of the
CENTROIDS matrix should be equal to the number of clusters and the columns
should equal the columns of the data.}
}
\value{
it returns a vector with the clusters.
}
\description{
Prediction function for mini-batch k-means applied to matrix-like objects.
}
\details{
This function takes the data and the output centroids and returns the
clusters.

This implementation relies very heavily on the
\code{\link[ClusterR]{MiniBatchKmeans}} implementation. We provide the
ability to work with other matrix-like objects other than base matrices (e.g,
DelayedMatrix and HDF5Matrix) through the \code{beachmat}
library.
}
\examples{
data(iris)
km = mini_batch(as.matrix(iris[,1:4]), clusters = 3,
               batch_size = 10, max_iters = 10)
clusters = predict_mini_batch(as.matrix(iris[,1:4]),
                             CENTROIDS = km$centroids)
}
\author{
Yuwei Ni
}
