\name{plotlowertri}
\alias{plotlowertri}
\title{
Generate a graphs showing correlation matrix (lower semi matrix) (Deprecated)
}
\description{
Generate graphs (lower semi matrix) showing lower semi matrix. These graphs are often used to show the structure of a correlation, similarity or dissimilarity matrix.
}
\usage{
plotlowertri(input, valuename = "r",
pchlist = c(19, 17, 15, 1, 5, 2, 7), interval = 6,
cex = 1, ncex = 1, int =1.2, add.number = TRUE,
size = FALSE, add.text = FALSE, show.legend = TRUE,
digits = 2)

}

\arguments{
  \item{input}{
The \code{input}, often a correlation or a distance matrix.
}
  \item{valuename}{
Name of the value that to show in the legend.
}
  \item{pchlist}{
A numberic vector specifying the shapes of points, see \code{pch par()}.
}
  \item{interval}{
Types of point shapes to show
}
  \item{cex}{
A number specifying the text size in the legend
}
  \item{ncex}{
Size of the text shown above each column.
}
  \item{int}{
Space between lines within the legend
}
  \item{add.number}{
If the column number should be shown
}
  \item{size}{
Whether the size of points should change with the value
}
  \item{add.text}{
Logical, If the number should be shown in the grid.
}
  \item{show.legend}{
Logical, If the legend should be appear.
}
  \item{digits}{
Number of digits for the label of each interval.
}
}
\details{
In the legend, space between lines could be adjusted by specifying \code{int}.
}
\value{
lower matrix plot
}


\references{
Zhang Qiaoying, Peng Shaolin, Zhang Sumei, Zhang Yunchun, Hou Yuping.(2008). Association of dormintant species in Guia hill Municipal Park of Macao. Ecology and Environment. 17:1541-1547
}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}
\seealso{
See Also \code{\link{plotnetwork}}, \code{\link[corrplot]{corrplot}}, \code{\link[ggcorrplot]{ggcorrplot}},
}
\examples{
data(testdata)
spmatrix <- data2mat(testdata)
result <- sp.pair(spmatrix)

## Check the legend for 0.00 to 0.33 (Unwanted label)
plotlowertri(result$Pearson, int = 0.5, cex=1.5)
title("Pearson Correlation Lower Matrix Plot")

## Change the size of points and reset the intervals.
## Warning: The lower matrix plot illustrating Pearson
## Correlation between each pair of species. Note the
## triangle didn't appeared in the plots, but have been
## added to the legend. This is due to the distribution
## of data. Be careful in selection of intervals.


plotlowertri(result$Pearson, int = 0.5, cex=1.5,
interval = 4, pchlist = c(19, 17, 15, 1, 5), size = TRUE)

title("Pearson Correlation Lower Matrix Plot")

## "Pure" dots, may have to add legend manually...
plotlowertri(result$Pearson, int = 0.5, cex=2.5,
interval = 4, pchlist = rep(19, 5), size = TRUE,
show.legend = FALSE)

title("Pearson Correlation Lower Matrix Plot")


## Using BCI data
library(vegan)
data(BCI)
## select the top 30 species according to relative frequency.
sub <- sub.sp.matrix(BCI, common = 30)
## Original
plotlowertri(cor(sub))

## Change size
plotlowertri(cor(sub), size = TRUE, cex = 3)

## Set the digits to 1
plotlowertri(cor(sub), size = TRUE, cex = 3, digits = 1,
ncex = 0.7)


}

\keyword{ lower }
\keyword{ matrix }
\keyword{ semimatrix }
