\name{freq.calc}
\alias{freq.calc}
\title{
Compute species' relative frequency
}
\description{
Computing species' relative frequency as defined by the numbers of plots having a species divided by the total number of plots.
}
\usage{
freq.calc(matr)
}

\arguments{
  \item{matr}{
    A community matrix
}
}
\details{
  The input should be a standard community matrix with rows representing sites and columns representing species.
}
\value{
A vector containing relative frequency for each species
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
