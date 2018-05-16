\name{plotnetwork}
\alias{plotnetwork}
\title{
Plot correlation network (Deprecated)
}
\description{
Plotting correlation network showing the relationship between each pair of sites connected by segments. The points are arranged in a circle.
}
\usage{
plotnetwork(datainput, interval = 8, xlim = c(-2.5,5),
ylim=c(-3.2,3.2), lty = c(1,2,3,4,4,3,2,1,5), value = "r",
legendx = 3, legendy = 0, right = 1.2, intcept = 0.22,
left = 0.35, linelength = 0.3, cex = 3, lwd = 1.5,
show.legend = TRUE, digits = 2, dit = 1.2,
number.label = TRUE, text.label = TRUE,
linecol = c("red", "black"), ...)
}
\arguments{
  \item{datainput}{
A correlation matrix
}
  \item{interval}{
Number of intervals indicating how to cut the data into different categories
}
  \item{xlim}{
The x range of the plot
}
  \item{ylim}{
The y range of the plot
}
  \item{lty}{
A numeric vector showing the styles of lines used to connect the points, must have n+1 elements.
}
  \item{value}{
Value of the matrix, should be a character, will used in legend.
}
  \item{legendx}{
Coordinate X for the Legend
}
  \item{legendy}{
Coordinate Y for the Legend
}
  \item{right}{
Move legend to the right
}
  \item{intcept}{
Specify the line space in the legend.
}
  \item{left}{
Move legend to the left
}
  \item{linelength}{
Length of the line in the legend.
}
  \item{cex}{
Point size for each circle.
}
  \item{lwd}{
Line width for each circle.
}
  \item{\dots}{
Other arguments to be passed from.
}
  \item{show.legend}{
Whether the legend should be drawn.
}
  \item{digits}{
Number of digits displayed in legend.
}
  \item{dit}{
Distance of text labels from each corner.
}
  \item{number.label}{
Whether the number label should be drawn.
}
  \item{text.label}{
Whether the text label should be drawn.
}
  \item{linecol}{
Colours of the lines specified. The positive correlation will be drawn in the first colour specified.
}

}
\details{
This function could be used to plot the pairwise connections between less than 20 sites
( above 20 is not recommended since there would be too many connections).

The lines will be in red or blue, according the sign of the value of association. Users can adjust the line style and legends based on their requirements.
}
\value{
Correlation network plots.
}
\references{
None
}
\author{
Jinlong Zhang \email{jinlongzhang01@gmail.com}
}
\examples{
data(testdata)
spmatrix <- data2mat(testdata)
result <- sp.pair(spmatrix)
plotnetwork(result$Pearson)

plotnetwork(result$Pearson, linecol = c("orange", "blue"),
number.label = FALSE)

title("Pearson Correlation Network")
}
\keyword{ network }
\keyword{ plot }
