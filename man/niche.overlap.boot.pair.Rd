\name{niche.overlap.boot.pair}
\alias{niche.overlap.boot.pair}
\title{
Boostrap the niche overlap indices
}
\description{
bootstrap the niche overlap indices between a pair of species. This is an internal function used by \code{niche.overlap.boot}, use \code{niche.overlap.boot} instead.
}
\usage{
niche.overlap.boot.pair(vectorA, vectorB, method = c("levins",
     "schoener", "petraitis", "pianka", "czech", "morisita"),
      times = 1000, quant = c(0.025, 0.975))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vectorA}{
A numeric vector containing species A's abundance or importance value.
}
  \item{vectorB}{
A numeric vector containing species B's abundance or importance value.
}
  \item{method}{
Name of the index to use.
}
  \item{times}{
Number of bootstraps
}
  \item{quant}{
Confidence intervals of the bootstrap results.
}
}

\value{
This function will return a vector containing:\
\code{"Observed", \
"Boot mean", \
"Boot std", \
"Boot CI1", \
"Boot CI2", \
"times"} \
}
\references{
Zhang Jin-tun,(2004 ) Quantitative Ecology, Science Press, Beijing
}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\note{
This is an internal function, please use \code{niche.overlap.boot}.
}

\seealso{
\code{\link{niche.overlap.boot}}
}
\examples{

data(datasample)
niche.overlap.boot.pair(datasample[,1], datasample[,2], method = "levins")

}

\keyword{ bootstrap }
\keyword{ niche }
\keyword{ overlap }
