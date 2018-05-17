\name{data2mat}
\alias{data2mat}
\title{
Convert field records to community matrix
}
\description{
Convert field records to community matrix
}
\usage{
data2mat(data = data)
}
\arguments{
  \item{data}{
A dataframe with the the following columns: \code{species}, \code{plots} or \code{sites}, \code{abundance}.
}
}

\value{
Return a community matrix ready for computing diversity indices.
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
