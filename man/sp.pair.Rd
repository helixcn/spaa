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

\code{c} = number of plots only have B.

\code{d} = number of plots without A or B.

\code{N} = \code{a+b+c+d}

Then, it is possible to compute:

Chi square (Yate's correction): \code{chi^{2}=((((a*d-b*c)-0.5*N)^2)*N)/(a+b)*(a+c)*(b+d)*(c+d)}

V ratio: \code{V = ((a+d)-(b+c))/(a + b + c + d)}

Jaccard index: \code{Jaccard =a/(a + b + c)}

Ochiai index: \code{Ochiai = a/sqrt((a+b)*(a+c))}

Dice index: \code{Dice = 2*a/(2*a + b + c)}


The Association Coefficient(\code{AC}):

if \code{a*d >= b*c}:

    \code{AC = (a*d - b*c)/((a+b)*(b+d))}

if \code{a*d < b*c  and a <= d}:

    \code{AC = (a*d - b*c)/((a+b)*(a+c))}

if \code{a*d < b*c and a > d}:

    \code{AC = (a*d - b*c)/((b+d)*(c+d))}

Point correlation coefficient

(\code{PCC}):

\code{PCC = {a*d-b*c}/{(a+b)*(a+c)*(c+d)*(b+d)}}

}

\value{

\item{chisq}{Chi Square matrix}
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

HURLBERT, S. H. (1969). A coefficient of interspecific association. Ecology, 50(1), 1-9.

WANG, B. S., & PENG S. L. (1985). Studies on the Measuring Techniques of Interspecific Association of Lower-Subtropical Evergreen-Broadleaved Forests. I. The Exploration and the Revision on the Measuring Formulas of Interspecific Association. Chinese Journal of Plant Ecology, 9(4), 274-285.

JIAN, M. F., LIU, Q. J., ZHU, D., & YOU, H. (2009). Inter-specific correlations among dominant populations of tree layer species in evergreen broad-leaved forest in Jiulianshan Mountain of subtropical China. Chinese Journal of Plant Ecology, 33(4), 672-680.

ZHOU, X. Y., WANG, B. S., LI, M. G., & ZAN, Q. J. (2000). An analysis of interspecific associations in secondary succession forest communities in Heishiding Natural Reserve, Guangdong Province. Chinese Journal of Plant Ecology, 24(3), 332-339

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
