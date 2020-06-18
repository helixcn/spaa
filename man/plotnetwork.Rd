\name{plotnetwork}
\alias{plotnetwork}
\title{
Plot correlation network (Deprecated)
}
\description{
Plotting correlation network showing the relationship between each pair of sites connected by segments. The points are arranged in a circle.
}
\usage{
plotnetwork(datainput, n_levels = 6, xlim = c(-2.5, 5),
    ylim = c(-2.5, 2.5), node_size = 3, lwd.var = TRUE,
    lwd = 4, label_dist = 1.2,
    show.node = TRUE, show.text.label = TRUE,
    linecol = c("orange", "blue"),show.legend = TRUE,
    valuename = "r", legendx = 3,
    legendy = -2, legend_line_space = 1,
    legend_linelength = 0.3, adjust_legend_x = 0,
    adjust_legend_y = 0,digits = 2, ... )
}
\arguments{
  \item{datainput}{
The correlation matrix, usually a lower matrix
}
  \item{n_levels}{
    Number of types of segments to show the correlation, note interval should <12
    }
  \item{xlim}{
    horizontal range of the canvas}
  \item{ylim}{
    vertical range of the canvas}
  \item{node_size}{
    size of the nodes}
  \item{lwd.var}{
    logical, if the segments width should vary with the absolute value
    }
  \item{lwd}{
    width of the segments for connecting the sites, default 1
    }
  \item{label_dist}{
    Distance of text labels from each node.
    }
  \item{show.node}{
    Whether the nodes in the figure should be labeled
    }
  \item{show.text.label}{
    Whether the text label should be drawn.
    }
  \item{linecol}{
    Colours showing positve or negative correlation. The lines representing positive correlations use the first element
    }
  \item{show.legend}{
    If the legend should be shown
    }
  \item{valuename}{
    Name of the variable shown in the legend
    }
  \item{legendx}{
    the starting position of the legend x
    }
  \item{legendy}{
    the starting position of the legend y
    }
  \item{legend_line_space}{
    the space between the lines in the legend
    }
  \item{legend_linelength}{
    length of the segment in the legend
    }
  \item{adjust_legend_x}{
    adjusting the position of legend (x)
    }
  \item{adjust_legend_y}{
    adjusting the position of legend (y)
    }
  \item{digits}{
    Number of disgits shown in the legend, it will be used in generating the breaks as well
    }
  \item{...}{
    other parameters related with plot.default
    }
}
\details{
This function could be used to plot the pairwise connections between less than 20 sites
( above 20 is not recommended since there would be too many connections).

The lines will be in orange or blue, according the sign of the value.
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

plotnetwork(result$Pearson, linecol = c("red", "black"))

plotnetwork(result$Pearson, n_breaks = 7,  node_size = 4,  lwd.var = FALSE, dit = 0.8, show.node = FALSE, show.text.label = FALSE, linecol = c("red", "black"), show.legend = TRUE, valuename = "r", legendx = 3, legend_line_space = 0.5, legend_linelength = 0.5, adjust_legend_x = -1)


title("Pearson Correlation Network")
}
\keyword{ network }
\keyword{ plot }
