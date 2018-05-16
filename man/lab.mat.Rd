\name{lab.mat}
\alias{lab.mat}
\title{
Convert "mXnY" vector to matrix
}
\description{
In forest dynamic plots, 20m quadrats are often named as "mXnY", for example, "X15Y01", "X10Y08".  This function will convert these labels to a matrix.
}
\usage{
lab.mat(plotlab)
}
\arguments{
  \item{plotlab}{
A character vector contains the names of the 20m quadrats.
}
}

\value{
A matrix with all the labels ordered.
}

\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\seealso{
\code{\link{turnover}}
}
\examples{

plotlab1 <- XYname(4,6)
labmat <- lab.mat(plotlab1)

}

\keyword{ turnover }
