\name{residue}

\alias{residue}

\title{Residue of a matrix}

\description{Returns the residue of a matrix.}

\usage{residue(Data)}

\arguments{
   \item{Data}{an \code{\linkS4class{ExpressionSet}} or a matrix}
}

\details{
  This function computes the residue of a matrix as defined
  by Yang et al (see references). 
}

\references{J. Yang, H. Wang, W. Wang, and P.S. Yu. An improved biclustering method for analyzing gene expression. International Journal on Artificial Intelligence Tools, 14(5):771-789, 2005}

\author{
  Pierre Gestraud
}

\seealso{
  \code{\link{FLOC}}
}

\examples{
data(sample.bicData)
residue(sample.bicData)
}

\keyword{cluster}
