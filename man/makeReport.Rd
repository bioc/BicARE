\name{makeReport}
\alias{makeReport}

\title{Export the results as html files}

\description{Creates a directory with html files containing the
  biclustering results.}

\usage{
makeReport(dirPath, dirName, resBic, browse=TRUE)
}


\arguments{
  \item{dirPath}{path to the directory}
  \item{dirName}{the name of the directory where the report will be created}
  \item{resBic}{a biclustering result}
  \item{browse}{logical. If TRUE the web browser will be opened}
}

\details{
  \code{makeReport} produces a html report of biclustering results in a new directory named \code{dirName}.  
  If the \code{browse} argument is set to TRUE the web browser will be opened on the
  "home.html" file. 
  
  Make sure to have rights to create the result directory.
}

\author{
  Pierre Gestraud \email{pierre.gestraud@curie.fr}
}

\examples{

data(sample.biclustering)
dirPath <- getwd() ## report created in the current working directory
dirName <- "test"
makeReport(dirPath, dirName, sample.biclustering, browse=FALSE)
}

\keyword{cluster}