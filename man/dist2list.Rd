\name{dist2list}
\alias{dist2list}

\title{
Convert distance matrix to pairwised list
}
\description{
Convert distance matrix to pairwised list
}
\usage{
dist2list(dist)
}

\arguments{
  \item{dist}{ distance matrix
}
}
\details{
Pairwise list with first column indicates the rows of the original distance matrix, second column indicates the columns indicates the rows of the original distance matrix, and the third indicats the values.
}
\value{
Dataframe with three columns.
}
\references{
Tuomisto, H. (2003). "Dispersal, Environment, and Floristic Variation of Western Amazonian Forests." Science 299(5604): 241-244.
}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\seealso{
\code{\link{list2dist}}
}
\examples{

##dist2list() example
x <- matrix(rnorm(100), nrow=5)
sampledata <- dist(x)
ddd <- dist2list(sampledata)

}

\keyword{ distance }
\keyword{ list }
