\name{niche.width}
\alias{niche.width}

\title{
Niche width
}
\description{
Compute niche width of the species in a community.
}
\usage{
niche.width(mat, method = c("shannon", "levins"))
}
\arguments{
  \item{mat}{
A community data matrix with each column for each species, and each row for each plot.
}
  \item{method}{
Index of niche width.
}
}
\details{
To be added
}
\value{
A vetor contains niche width index of species in community.
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

# ### niche.width() example 
data(datasample)
niche.width(datasample, method = "levins")
niche.width(datasample, method = "shannon")

}

\keyword{ niche }
