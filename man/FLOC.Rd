\name{FLOC}

\alias{FLOC}

\title{Performs the FLOC algorithm}

\description{Find a given number of biclusters using the a modified version of the FLOC algorithm.}

\usage{FLOC(Data, k = 20, pGene = 0.5, pSample=pGene, r = NULL, N = 8, M = 6, t = 500, blocGene = NULL, blocSample = NULL)}

\arguments{
   \item{Data}{an \code{\linkS4class{ExpressionSet}} or a matrix (with
   genes on rows and conditions on columns)}
   \item{k}{the number of biclusters searched}
   \item{pGene}{genes initial probability of membership to the biclusters}
   \item{pSample}{samples initial probability of membership to the biclusters}
   \item{r}{the residue threshold}
   \item{N}{minimal number of gene per bicluster}
   \item{M}{minimal number of conditions per bicluster}
   \item{t}{number of iterations}
   \item{blocGene}{a matrix indicating the directed initialisation for
   the genes (see details)}
   \item{blocSample}{a matrix indicating the directed initialisation for
   the conditions (see details)}
}

\details{This biclustering algorithm is based on the FLOC algorithm (FLexible
  Overlapped biClustering) defined by Yang et al. (see references). It
  can discover a set of \code{k}, possibly overlapping, biclusters. If \code{r} is set to
  NULL, the residue threshold used in the analysis is the residue of
  \code{Data} divided by 10.
  
  \code{blocGene} and \code{blocSample} are matrix of 0 and 1 with the rows
  representing the features (gene or samples) and the columns the
  biclusters. A 1 on line i and column j indicates that the feature i
  (gene or sample) will be include in the bicluster j during the
  initialisation step and will not be removed from it during the analysis. If the number of columns in these matrices is different from the number of bicluster searched, \code{k} is set to the maximal value
  of these two.
  
  See \code{\link{bicluster}} to extract a bicluster from the biclustering result. 
}

\value{
  Returns an object of class 'biclustering', a list containing at
  least :
  \item{Call}{the matched call.}
  \item{ExpressionSet}{the data used}
  \item{param}{a data.frame with the algorithm parameters}
  \item{bicRow}{a matrix of boolean indicating the belonging of the
    genes to the biclusters}
  \item{bicCol}{the same as for bicRow but for the conditions}
  \item{mat.resvol.bic}{a matrix describing the biclusters}
}

\references{J. Yang, H. Wang, W. Wang, and P.S. Yu. An improved biclustering method for analyzing gene expression. International Journal on Artificial Intelligence Tools, 14(5):771-789, 2005}

\author{
  Pierre Gestraud (\email{pierre.gestraud@curie.fr})
}

\examples{

data(sample.bicData)     ## subset of sample.ExpressionSet from Biobase
residue(sample.bicData)  ##  0.3401921
resBic <- FLOC(sample.bicData, k=10, pGene=0.5,r=0.05,N=8,M=10,t=500)
resBic

## initialising samples of 2 biclusters
iniSample <- matrix(0, ncol=2, nrow=26)
## first bicluster initialised around Female cases
iniSample[pData(sample.bicData)$sex=="Female",1] <- 1
## second bicluster initialised around control cases
iniSample[pData(sample.bicData)$type=="Control",2] <- 1
resBic <- FLOC(sample.bicData, k=10, pGene=0.5, r=0.05, N=8, M=10, t=500, blocSample=iniSample)
resBic
}


\keyword{cluster}