\name{niche.overlap.boot}
\alias{niche.overlap.boot}
\title{
Boostrap of niche overlap
}
\description{
Bootstrap of niche overlap between species.
}
\usage{
niche.overlap.boot(mat, method = c("pianka", "schoener", "petraitis",  
    "czech", "morisita", "levins"), times = 999, quant = c(0.025, 0.975))
}
\arguments{
  \item{mat}{
A community data matrix with columns representing species, and rows representing sites.
}
  \item{method}{
A character string indicating the index of niche overlap to be applied.
}
  \item{times}{
Interger, Number of bootstrap to be implemented.
}
  \item{quant}{
Quantile of bootstrap results, the confidence intervals.
}
}
\details{
This function bootstraps the niche overlap within each pair of species. \
pianka: Pianka's niche overlap index\
schoener: Schoener's niche overlap index\
petraitis: Petraitis' niche overlap index\
czech: Czechanowski index \
morisita: Morisita's overlap index\
levins: Levin's overlap index\
see more information from Gotelli, N(2009).\
}
\value{

a data frame with each row for each pair of species 
the columns are \code{"Observed", \
"Boot mean", \
"Boot std", \
"Boot CI1",  \
"Boot CI2",  \
"times"} \
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

