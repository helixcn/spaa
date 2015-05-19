\name{geodist}
\alias{geodist}
\title{
Hight precision Great circle distance between two places
}
\description{
Hight precision Great circle distance between two places assuming the earth is elliptic sphere. 
}
\usage{
geodist(L1, phi1, L2, phi2)
}

\arguments{
  \item{L1}{
Longitude of first place in decimal format.
}
  \item{phi1}{
Latitude of first place in decimal format.
}
  \item{L2}{
Longitude of second place in decimal format.
}
  \item{phi2}{
Latitude of second place in decimal format.
}
}
\details{
Hight precision great circle distance between two places assuming the earth is elliptic sphere.
}
\value{
Hight precision great circle distance.
}
\references{
Jean Meeus 1991 Astronomical Algorithms Willmann-Bell 80-83
}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\seealso{
\code{\link{lgeodist}}
}
\examples{

## geodist() example
## Paris
L1 = deg2dec(-2,20,14)
phi1 = deg2dec(48, 50, 11)
## Washington DC
L2 = deg2dec(77,03,56)
phi2 = deg2dec(38,55,17)
##High precision Great Circle distance
geodist(L1, phi1, L2, phi2)

}

\keyword{ distance }
