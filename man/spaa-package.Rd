\name{spaa-package}
\alias{spaa-package}
\alias{spaa}
\docType{package}
\title{
\packageTitle{spaa}
}
\description{
\packageDescription{spaa}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{spaa}
\packageIndices{spaa}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
\packageAuthor{spaa}

Maintainer: \packageMaintainer{spaa}
}
\references{
Zhang Jin-tun,(2004 ) Quantitative Ecology, Science Press, Beijing
}

\keyword{ package }
\seealso{
\code{\link[cooccur:cooccur-package]{cooccur}}, \code{\link[ecospat:ecospat-package]{ecospats}}
}
\examples{

  library(vegan)
  data(BCI)
  ## select the top 30 species according to relative frequeny.
  sub <- sub.sp.matrix(BCI, common = 30)
  ## Set the digits to 1
  plotlowertri(cor(sub), size = TRUE, cex = 3, digits = 1)

  #### Niche width and niche overlap
  data(datasample)
  niche.overlap.boot(datasample[,1:3], method = "levins")
  niche.overlap(datasample, method = "levins")
  niche.width(datasample[,1:3], method = "shannon")

}
