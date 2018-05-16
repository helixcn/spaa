\name{niche.overlap.pair}
\alias{niche.overlap.pair}
\title{
Compute niche overlap index between two species
}
\description{
Compute niche overlap index between two species. This is an internal function, used \code{\link{niche.overlap}} instead.
}
\usage{
niche.overlap.pair(vectA, vectB, method = c("pianka",
     "schoener","petraitis","czech","morisita", "levins"))
}
\arguments{
  \item{vectA}{
A numeric vector containing species A's abundance or importance value
}
  \item{vectB}{
A numeric vector containing species B's abundance or importance value
}
  \item{method}{
Niche overlap index
}
}

\details{
None
}

\value{
The niche overlap index
}

\references{

Zhang Jin-tun,(2004 ) Quantitative Ecology, Science Press, Beijing

Nicholas J. Gotelli. 2000. Null model analysis of species co-occurrence patterns. Ecology 81:2606-2621. http://esapubs.org/archive/ecol/E081/022/EcoSim%20Help/Niche/NicheOverlap.htm

}

\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\seealso{
\code{\link{niche.overlap}}
}

\examples{

data(datasample)
niche.overlap.pair(datasample[,1],datasample[,2], method = "levins")

}

\keyword{ niche }
\keyword{ overlap }
