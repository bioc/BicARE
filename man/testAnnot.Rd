\name{testAnnot}

\alias{testAnnot}

\title{Find samples annotations over-represented covariates in biclusters}

\description{Characterisation of the biclusters in term of over-representation of sample covariates.}

\usage{testAnnot(resBic, annot=NULL, covariates="all")}

\arguments{
	\item{resBic}{a biclustering result from \link{FLOC}}
  	\item{annot}{annotation matrix, default value is set to NULL, then phenoData of the ExpressionSet is used}
  	\item{covariates}{the names of the covariates that should be tested,
    default value is set to "all"}
}

\details{
  	For each bicluster and each covariate a chi-squarred test is performed to test the adequation between the distribution of the levels of the covariates in the bicluster and in the original dataset.  
  
	Multiple testing correction is performed by the Benjamini-Yekutieli procedure. The residuals of the tests indicate if the level is over or down represented in the bicluster.  
    
    Due to the amount of results it is advised to use the \code{\link{makeReport}} function to get a html report. 
}

\value{
	A \code{biclustering} object containing \code{resBic} and updated with
  	the results of the tests in \code{resBic$covar}.
 
  	The results are presented as a list with : 
  	\item{covar}{the samples covariates tested}
  	\item{pvalues}{a matrix with the p-values of the tests}
  	\item{adjpvalues}{a matrix with the p-values adjusted by the Benjamini Yekutieli procedure}
  	\item{index}{a list of matrices with the numbers of each level in each bicluster}
  	\item{residuals}{a list of matrices with the residuals of the tests for each modality in each bicluster}
}
\author{Pierre Gestraud}

 	
\examples{
data(sample.biclustering)
resBic <- testAnnot(sample.biclustering, annot=NULL, covariates=c("sex", "type"))
}

\keyword{cluster}