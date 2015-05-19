\name{data2mat}
\alias{data2mat}
\title{
Convert species list data to species matrix
}
\description{
This function can be used to convert the species list to species matrix. The rows  of the output matrix are plots, or sites. The columns are the species.
}
\usage{
data2mat(data = data)
}
\arguments{
  \item{data}{
The input \code{data}}
}
\details{
The input data will have to include :\code{species},\code{plots} or \code{sites}, \code{abundance}, specifically, a column named \code{"abundance"} must be specified. 
}
\value{
Return a species matrix with each row for each plot, and each column for species.
}
\references{
None
}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}
\examples{
data(testdata)
spmatrix <- data2mat(testdata)
}

\keyword{ species }
\keyword{ matrix }

