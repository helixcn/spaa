\name{niche.overlap}
\alias{niche.overlap}
\title{
Niche overlap between each pair of species
}
\description{
Compute niche overlap between each pair of species.
}
\usage{
niche.overlap(mat, method = c("levins", "schoener",
       "petraitis", "pianka", "czech", "morisita"))
}

\arguments{

  \item{mat}{
A community matrix with columns representing species, and rows representing plots.
}
  \item{method}{
A string specifying the name of the index.
}

}
\details{
To add.
}
\value{
A distance matrix contains niche overlap index between each pair of species.
}
\references{

Zhang Jin-tun,(2004 ) Quantitative Ecology, Science Press, Beijing

Nicholas J. Gotelli. 2000. Null model analysis of species co-occurrence patterns. Ecology 81:2606-2621. http://esapubs.org/archive/ecol/E081/022/EcoSim%20Help/Niche/NicheOverlap.htm

}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\seealso{
\code{\link{niche.overlap.pair}}
}
\examples{

data(datasample)
niche.overlap(datasample, method = "levins")

}

\keyword{ niche }
\keyword{ overlap }
