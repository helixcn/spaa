\name{plotlowertri}
\alias{plotlowertri}
\title{
Plotting lower semi matrix
}
\description{
Function for plotting lower semi matrix. These plots are often used to illustrate the relationship in Pearson's correlation, similarity or dissimilarity index between sites or species. 
}
\usage{
plotlowertri(input, valuename = "r", 
pchlist = c(19, 17, 15, 1, 5, 2, 7), interval = 6, 
cex = 1, ncex = 1, int =1.2, add.number = TRUE, 
size = FALSE, add.text = FALSE, show.legend = TRUE, 
digits = 2)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
The \code{input} data, often the results of correlation matrix, can also be class \code{dist} object.
}
  \item{valuename}{
Value name that will be used in the legend.
}
  \item{pchlist}{
The types of point shape to plot see \code{pch par()}. 
}
  \item{interval}{
Number of intervals to illustrate the shape of points
}
  \item{cex}{
Point size
}
  \item{ncex}{
Text size for the add.number, which appeared at the top of each column.
}
  \item{int}{
Row space between lines in legend
}
  \item{add.number}{
Species number for each column.
}
  \item{size}{
Whether the size of points would change according to value.
}
  \item{add.text}{
To add value to the grid.
}
  \item{show.legend}{
Whether the legend should be drawn.
}
  \item{digits}{
Number of digits of for each interval
}
}
\details{
If the matrix contains less than 15 rows/columns, you may have to adjust the row space between the text lines in the legend, using argument \code{int}. Data in class \code{dist} can be include, and will be converted to matrix at first internally.

The lower matrix plot illustrating Pearson's Correlation between each pair of species. Note some value didn't appeared in the plots, may have appeared the legend. This is due to the distribution of data. Be careful in selection of intervals. In this situation you may set \code{show.legend = FALSE}, and add the legend manually. This may be fixed in the future.
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
See Also \code{\link{plotnetwork}}
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
