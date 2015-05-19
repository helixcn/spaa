\name{lgeodist}
\alias{lgeodist}
\title{
Low precision Great circle distance between two places 
}
\description{
Calculating Great circle distance between two places assuming that the earth is sphere.
}
\usage{
lgeodist(L1, phi1, L2, phi2)
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

\value{
Low precision great circle distance between two places.
}
\references{
Jean Meeus 1991 Astronomical Algorithms Willmann-Bell 80-81
}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}
\note{
This function assuming that the earth is sphere.
}

\seealso{
\code{\link{geodist}}
}
\examples{

#lgeodist() example
##Paris
L1 = deg2dec(-2,20,14)
phi1 = deg2dec(48, 50, 11)
##Washington DC
L2 = deg2dec(77,03,56)
phi2 = deg2dec(38,55,17)
#Great circle distance
lgeodist(L1, phi1, L2, phi2)

}

\keyword{ distance }

