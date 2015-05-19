\name{splist}
\alias{splist}
\docType{data}
\title{
species list used in function add.col()
}
\description{
A species check list example to be used in add.col()
}
\usage{data(splist)}
\format{
  A data frame with 9 observations on the following 3 variables.
  \describe{
    \item{\code{species}}{a factor with levels \code{sp1} to \code{sp8}}
    \item{\code{genera}}{a factor with levels \code{gen1} to \code{gen6}}
    \item{\code{family}}{a factor with levels \code{fam1} to \code{fam5}}
  }
}

\references{
None
}
\examples{

data(splist)
data(testdata)
## add genera from dataframe B to dataframe A.
add.col(inputA = testdata, inputB = splist, add = "genera",
according = "species")
## add family from dataframe B to dataframe A.
add.col(inputA = testdata, inputB = splist, add = "family", 
according = "species")

}
\keyword{datasets}
