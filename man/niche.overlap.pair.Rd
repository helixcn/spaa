\name{niche.overlap.pair}
\alias{niche.overlap.pair}
\title{
Niche overlap between one pair of species.
}
\description{
Compute niche overlap index between one pair of species. Users are encouraged to used \code{\link{niche.overlap}} istead of this function.
}
\usage{
niche.overlap.pair(vectA, vectB, method = c("pianka", 
     "schoener","petraitis","czech","morisita", "levins"))
}
\arguments{
  \item{vectA}{
A numerical vector including species A's abundance or value of importance.
}
  \item{vectB}{
A numerical vector including species B's abundance or value of importance.
}
  \item{method}{
Niche overlap index to be applied.
}
}

\details{
None
}

\value{
Niche overlap index.
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

### niche.overlap.pair() example
data(datasample)
niche.overlap.pair(datasample[,1],datasample[,2], method = "levins")

}

\keyword{ niche }
\keyword{ overlap }
