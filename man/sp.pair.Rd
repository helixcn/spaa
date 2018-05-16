\name{sp.pair}
\alias{sp.pair}

\title{
Species association between each pair of species
}
\description{
Compute species association between each pair of species.
}
\usage{
sp.pair(matr)
}

\arguments{
\item{matr}{
Standard community matrix, with rows representing sites and columns representing species.
}
}
\details{
If \code{a}, \code{b}, \code{c}, \code{d} denote the co-occurrence the two species A and B, where:

\code{a} = number of plots occupied both by A and B.

\code{b} = number of plots only have A.

\code{c} = number of plots only have sB.

\code{d} = number of plots without A or B.

\code{N} = \code{a+b+c+d}

Then, it is possible to compute:

Chi square (Yate's correction): \code{chi^{2}=((((a*d-b*c)-0.5*N)^2)*N)/(a+b)*(a+c)*(b+d)*(c+d)}

V ratio: \code{V = ((a+d)-(b+c))/(a + b + c + d)}

Jaccard index: \code{Jaccard =a/(a + b + c)}

Ochiai index: \code{Ochiai = a/sqrt((a+b)*(a+c))}

Dice index: \code{Dice = 2*a/(2*a + b + c)}

The Association Coefficient(\code{AC}):

if \code{a*d>= b*c}:

\code{AC = (a*d - b*c)/((a+b)*(b+d))}

if \code{b*c>= a*d and d>=a}:

\code{AC=(a*d - b*c)/((a+b)*(a+c))}

if \code{b*c>a*d and a<a}:

\code{AC = (a*d - b*c/((b+d)(d+c))}

Point correlation coefficient

(\code{PCC}):

\code{PCC = {a*d-b*c}/{(a+b)*(a+c)*(c+d)*(b+d)}}

}

\value{

\item{chisq}{Chi square matrix}
\item{chisqass}{information for the Chi square matrix }
\item{V}{\code{V} positive or negative association}
\item{Ochiai}{Ochiai's index}
\item{Dice}{Dice's index}
\item{Jaccard}{Jaccard's index}
\item{Pearson}{Pearson's correlation coefficient}
\item{Spearman}{Spearman's rank correlation coefficient}
\item{PCC}{Point correlation coefficient}
\item{AC}{Association coefficient}

}
\references{
Zhang Qiaoying, Peng Shaolin, Zhang Sumei, Zhang Yunchun, Hou Yuping.(2008). Association of dormintant species in Guia hill Municipal Park of Macao.  Ecology and Environment. 17:1541-1547

Zhou XY, Wang BS, Li MG, Zan QJ.(2000). An analysis of interspecific associations in secondary  succession forest communities in Heishiding Nature Reserve, Guangdong Province. Acta Phytoecologica Sinica. 24:332-339.

JIAN Minfei, LIU qijing, ZHU du, YOU hai.(2009). Inter-specific correlations among dorminant populations of tree layer species in evergreen broad-leaved forest in Jiulianshan Mountain of subtropical China. Chinese Journal of Plant Ecology. 33: 672-680
}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\seealso{
See Also as \code{\link{sp.assoc}} for computing association for all the species.
}
\examples{
data(testdata)
spmatrix <- data2mat(testdata)
result <- sp.pair(spmatrix)
}

\keyword{ species }
\keyword{ association }
