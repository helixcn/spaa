\name{niche.overlap.boot.pair}
\alias{niche.overlap.boot.pair}
\title{
Niche overlap boostrap utility function
}
\description{
Compute the bootstrap value between two vectors. This is a internal function called by \code{niche.overlap.boot}, users are encouraged to use the latter function.
}
\usage{
niche.overlap.boot.pair(vectorA, vectorB, method = c("levins", 
     "schoener", "petraitis", "pianka", "czech", "morisita"), 
      times = 999, quant = c(0.025, 0.975))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vectorA}{
A numerical vector including species A's abundance or value of importance.
}
  \item{vectorB}{
A numerical vector including species B's abundance or value of importance.
}
  \item{method}{
Niche overlap indeces to be applied.
}
  \item{times}{
Number of bootstraps
}
  \item{quant}{
Confidence interval of the bootstrap results.
}
}
\details{
To do.
}
\value{
This function will return a vector including the following elements:\
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
Users are ecouraged to call \code{\link{niche.overlap.boot}} rather than this function.
}

\seealso{
\code{\link{niche.overlap.boot}}
}
\examples{

### niche.overlap.boot.pair() example
data(datasample)
niche.overlap.boot.pair(datasample[,1],datasample[,2], method = "levins")

}

\keyword{ bootstrap }
\keyword{ niche }
\keyword{ overlap }
