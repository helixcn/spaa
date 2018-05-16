\name{niche.width}
\alias{niche.width}

\title{
Niche width
}
\description{
Compute niche width for all the species in a community.
}
\usage{
niche.width(mat, method = c("shannon", "levins"))
}
\arguments{
  \item{mat}{
A community matrix with columns representing species, and rows representing plots.
}
  \item{method}{
Character string showing the name of the index.
}
}

\value{
A vetor containing niche width for all the species in the community.
}
\references{
Zhang Jin-tun,(2004 ) Quantitative Ecology, Science Press, Beijing
}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\seealso{
\code{\link{niche.overlap}} for niche overlap
}

\examples{

data(datasample)
niche.width(datasample, method = "levins")
niche.width(datasample, method = "shannon")

}

\keyword{ niche }
