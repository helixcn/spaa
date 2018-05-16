\name{sub.sp.matrix}
\alias{sub.sp.matrix}
\title{
Subset species based on relative frequency
}
\description{
Subset species based on relative frequency.}
\usage{
sub.sp.matrix(spmatrix, freq = 0.5, common = NULL)
}

\arguments{
  \item{spmatrix}{
a standard community matrix with rows representing sites and columns representing species.
}
  \item{freq}{
The relative frequency, species with higher relative frequency will be kept in the output.
}
  \item{common}{
The number of most common species to keep.
}
}
\details{
sub.sp.matrix will select the species whose relative frequency above 0.5 (default), or select certain number of species based on relative frequency.
}
\value{
A subset matrix of species with high relative frequency.
}
\references{
None
}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\seealso{
See Also \code{\link{subset}}
}
\examples{

library(vegan)
data(BCI)
## Select the species whose relative frequency
## more than 0.5, from BCI data
sub <- sub.sp.matrix(BCI, freq = 0.5)
## Select the top 30 species according to relative frequency
sub <- sub.sp.matrix(BCI, common = 30)
}

\keyword{ sub }
\keyword{ matrix }
\keyword{ frequency }
