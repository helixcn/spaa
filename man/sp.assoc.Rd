\name{sp.assoc}
\alias{sp.assoc}
\title{
Analyzing species association
}
\description{
Analyzing species association
}
\usage{
sp.assoc(matr)
}
\arguments{
  \item{matr}{
standard community matrix , with rows representing sites and columns representing species.
}
}
\details{
Computations are based on the following formula.

If, \code{N}: Number of plots.

\code{S}: Number of species.

\code{n}: Number of plots occupied by certain species.

\code{Tj}: total number of species for each plot.
s
\code{t}: mean species number for all the plots.


Then:
Variance of species relative frequency is: \code{sigma^{2}{T}= sum{i}=1^{s}P{i}(1-P{i})}.

Variance of species number is: \code{ S^{2}{T}=({1}{N})sum{j=1}^{N}(T{j}-t)^{2} }.

Species relative frequency: \code{P{i}={n{i}}{N}}.

Variance ratio:

If \code{VR > 1} Positively associated,

If \code{VR < 1} Negative associated

\code{VR = {S{T}^{2}}/{sigma{T}^{2}}}

\code{W}: used in comparison with chi square with n degrees of freedom.

\code{ W = VR * N }

}

\value{
Variance ratio, W, Chisq, etc, see \code{details}
\item{pi}{Species frequency}
\item{N}{Number of plots}
\item{S}{Number of species}
\item{Tj}{Total number of species for each plot}
\item{Numspmean}{Mean number of species}
\item{sigmaTsq}{Variance of species relative frequency}
\item{STsq}{Variance of species number}
\item{var.ratio}{Variance ratio}
\item{W}{W statiscit value: used in comparison with chi square.(n)}
}
\references{
Zhang Qiaoying, Peng Shaolin, Zhang Sumei, Zhang Yunchun, Hou Yuping. (2008) Association of dormintant species in Guia hill Municipal Park of Macao. Ecology and Environment. 17:1541-1547

GUO zhongling, MA yuandan, ZHENG Jiping, LIU Wande , JIN Zefeng.(2004)
Biodiversity of tree species,their populations'spatial distribution pattern and interspecific association in mixed deciduous broadleaved forest in Changbai Mountains. Chinese Journal of Applied Ecology. 15:2013-2018

Shi Zuomin, Liu Shirong, Cheng Ruimei, Jiang Youxu.(2001) Interspecific association of plant populations in deciduous broad leaved forest in Baotianman. Scientia Silvae Sinicae. 37:30-35

}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}

\seealso{
See also \code{\link{sp.pair}} for association between each pair of species.
}
\examples{
data(testdata)
spmatrix <- data2mat(testdata)
sp.assoc(spmatrix)
}

\keyword{species}
\keyword{association}
