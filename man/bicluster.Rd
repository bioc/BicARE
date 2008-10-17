\name{bicluster}
\alias{bicluster}
\title{Extract a bicluster}

\description{Extract a bicluster from an object of class biclustering}

\usage{bicluster(biclustering, k, graph=TRUE)}

\arguments{
  \item{biclustering}{an object of class "biclustering" created by
    function \code{\link{FLOC}}}
  \item{k}{the number of the bicluster considered in the "biclustering" object}
  \item{graph}{boolean, indicating whether the graph should be ploted
    or not}
}

\value{Returns the bicluster as a matrix with the genes on rows and the
  samples on columns. Result marix is of class "bicluster". The "graph" option allows to plot the expression profiles of the genes across the conditions in the bicluster.}

\author{Pierre Gestraud}

\examples{
### extract the first bicluster
data(sample.biclustering)
sample.biclustering
bic <- bicluster(sample.biclustering, 1, graph=TRUE)
plot(bic)
}

\keyword{cluster}