\name{list2dist}
\alias{list2dist}
\title{
Convert pairwise list to distance matrix
}
\description{
Convert pairwise list to distance matrix
}
\usage{
list2dist(dat)
}
\arguments{
  \item{dat}{
dataframe with three columns
}
}
\details{
Dataframe with first column as the column names in the distance matrix, second column as the rownames in the distance matrix, third column the values.
}
\value{
distance matrix
}
\references{
Tuomisto, H. (2003). "Dispersal, Environment, and Floristic Variation of Western Amazonian Forests." Science 299(5604): 241-244.
}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\seealso{
\code{\link{dist2list}}
}
\examples{

##list2dist() example
x <- matrix(rnorm(100), nrow=5)
sampledata <- dist(x)
ddd <- dist2list(sampledata)
list2dist(ddd)

}
\keyword{ distance }

