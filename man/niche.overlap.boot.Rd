\name{niche.overlap.boot}
\alias{niche.overlap.boot}
\title{
Boostrap niche overlap indices
}
\description{
Bootstrap niche overlap indices
}
\usage{
niche.overlap.boot(mat, method = c("pianka", "schoener", "petraitis",
    "czech", "morisita", "levins"), times = 1000, quant = c(0.025, 0.975))
}
\arguments{
  \item{mat}{
standard community matrix.
}
  \item{method}{
character string specifying the index.
}
  \item{times}{
Interger, representing the number of bootstrap samples to generate.
}
  \item{quant}{
Quantile of the bootstrap values.
}
}
\details{
This function bootstraps the following niche overlap indices between each pair of species: \
\code{schoener}: Schoener's niche overlap index\
\code{petraitis}: Petraitis' niche overlap index\
\code{czech}: Czechanowski index \
\code{morisita}: Morisita's overlap index\
\code{levins}: Levin's overlap index\
see more information from Gotelli, N(2009).\
}
\value{

a data frame with each row for each pair of species
the columns are \code{"Observed", "Boot mean", "Boot std", "Boot CI1",  "Boot CI2",  "times"} 
}
\references{

Zhang Jin-tun,(2004 ) Quantitative Ecology, Science Press, Beijing\

Nicholas J. Gotelli. 2000. Null model analysis of species co-occurrence patterns. Ecology 81:2606-2621. http://esapubs.org/archive/ecol/E081/022/EcoSim%20Help/Niche/NicheOverlap.htm \

}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\seealso{
\code{\link{niche.overlap.boot.pair}}
}
\examples{

data(datasample)
niche.overlap.boot(datasample[,1:4], method = "pianka")
niche.overlap.boot(datasample[,1:4], method = "schoener")
niche.overlap.boot(datasample[,1:4], method = "czech")
niche.overlap.boot(datasample[,1:4], method = "levins")
}

\keyword{ niche}
\keyword{ bootstrap }
\keyword{ overlap }
