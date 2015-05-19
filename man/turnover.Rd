\name{turnover}
\alias{turnover}
\title{
Calculating species turnover
}
\description{
Calculating species turnover based on the mean value between focus quadrat and their neighbours.
}
\usage{
turnover(lab.mat, dist.mat, type = c("quart", "octal"))
}

\arguments{
  \item{lab.mat}{
matrix of quadrat labels.
}
  \item{dist.mat}{
distance matrix between each pair of quadrats
}
  \item{type}{
"quart" indicates four neighbouring quadrats, "octal" indicate eight neighbouring quadrats.
}
}
\details{
species turnover based on the mean value between centred quadrat and its neighbours.
}
\value{
matrix with species turnover.
}
\references{
Lennon J. 2001 The geographical structure of British bird distributions - diversity, spatial turnover and scale Journal of Animal Ecology 70,966-979
}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\seealso{
\code{\link{XYname}} and \code{\link{lab.mat}}
}
\examples{

##example turnover()
plotlab1 <- XYname(4,6)
xxx <- 1:240
dim(xxx) <- c(24, 10)
rownames(xxx) <- plotlab1
### Distance between each pair of plots
ddd <- dist(xxx)
### label matrix
labmat1 <- lab.mat(plotlab1)
yyy <- turnover(labmat1, ddd)

}

\keyword{ turnover }

