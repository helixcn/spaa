\name{freq.calc}
\alias{freq.calc}
\title{
Species relative frequency
}
\description{
This function calculates the species relative frequency which equals to the numbers of occupied plots partitioned by the total number of plots for each species.
}
\usage{
freq.calc(matr)
}

\arguments{
  \item{matr}{ The standard species matrix
}
}
\details{
  The input data is a  standard species matrix with rows for plots and column for species.
}
\value{
Returns a vector that contains relative frequency for each species included in the input matrix.
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
freq.calc(spmatrix)
}


\keyword{ frequency }
\keyword{ species }

