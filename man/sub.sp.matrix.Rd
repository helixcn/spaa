\name{sub.sp.matrix}
\alias{sub.sp.matrix}
\title{
Subset species matrix according to relative frequency
}
\description{
Subset species matrix according to relative frequency.}
\usage{
sub.sp.matrix(spmatrix, freq = 0.5, common = NULL)
}

\arguments{
  \item{spmatrix}{
The spmatrix is a standard species matrix with rows for plots and column for species.
}
  \item{freq}{
The value of relative frequency to be specified, species with higher relative frequency will be reserved in the output matrix. 
}
  \item{common}{
The number of most common (according to relative frequency) species to be specified. 
}
}
\details{
sub.sp.matrix will select the species whose relative frequency above 0.5 (default), or select certain number of species according to relative frequency.  
}
\value{
A subset matrix with certain number of top relative frequency.
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

