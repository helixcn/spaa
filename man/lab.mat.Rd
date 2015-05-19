\name{lab.mat}
\alias{lab.mat}
\title{
Convert vector of XY labels to label matrix
}
\description{
Convert vector of XY labels to label matrix
}
\usage{
lab.mat(plotlab)
}
\arguments{
  \item{plotlab}{
Vector of XY labels
}
}

\value{
XY label matrix
}

\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\seealso{
\code{\link{turnover}}
}
\examples{

### lab.mat() example
plotlab1 <- XYname(4,6)
labmat <- lab.mat(plotlab1)

}

\keyword{ turnover }
