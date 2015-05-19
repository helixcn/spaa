\name{deg2dec}
\alias{deg2dec}
\title{
Degree to decimal
}
\description{
Convert latitude or longitude from degree to decimal format
}
\usage{
deg2dec(h, m, s)
}
\arguments{
  \item{h}{
Degree
}
  \item{m}{
Minute
}
  \item{s}{
Second
}
}
\details{
Convert latitude or longitude from degree to decimal format.
}
\value{
Degree of decimal format
}

\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\note{
Places with eastern hemisphere should have longitude and southern hemisphere less than zero. 
}

\examples{

## deg2dec() example
##Paris
L1 = deg2dec(-2,20,14)
phi1 = deg2dec(48, 50, 11)
##Washington DC
L2 = deg2dec(77,03,56)
phi2 = deg2dec(38,55,17)

}

\keyword{ decimal }
