\name{testSet}
\alias{testSet}
\title{Find gene sets that are enriched in a bicluster}
\description{Test of the over-representation of gene sets in the biclusters}
\usage{testSet(resBic, geneSetCol)}

\arguments{
  \item{resBic}{a \code{biclustering} object created by \code{\link{FLOC}}}
  \item{geneSetCol}{a \code{\linkS4class{GeneSetCollection}}}
}

\details{
	The over-representation of a gene set in a bicluster is evaluated by an hypergeometric test. 
	    
	The genes identifiers of the gene sets will automaticaly be mapped to
  	the same as those used in the data. 
  
  	Due to the amount of results it is advised to use the \code{\link{makeReport}} function to get a html report. 
}

\value{
  A \code{biclustering} object containing \code{resBic} and updated with
  the results of the tests in \code{resBic$geneSet}.
  
  The results are presented as a list with : 
  \item{GeneSetCollection}{the GeneSetCollection used}
  \item{pvalues}{a matrix containing the pvalues of the tests for each geneSet and each bicluster}
  \item{adjpvalue}{a matrix containing the p-values adjusted by the Benjamini Yekutieli procedure}
}

\author{Pierre Gestraud \email{pierre.gestraud@curie.fr}}

\examples{
data(sample.biclustering)
gss <- GeneSetCollection(sample.biclustering$ExpressionSet[1:50,], setType=GOCollection())
resBic <- testSet(sample.biclustering, gss)
}

\keyword{cluster}