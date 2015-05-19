\name{add.col}
\alias{add.col}
\title{
Add one column from A to B, according to one column of common name.
}
\description{
This function can be used to add one column from dataframe B to dataframe A, according to the column names speciefied.
}
\usage{
add.col(inputA, inputB, add, according)
}

\arguments{
  \item{inputA}{
\code{inputA} A dataframe which a column to be added according to one row information. Must contain a column name equals to one column name in dataframe B at the least.
}
  \item{inputB}{
\code{inputB} A dataframe which a column to be abstracted from, must contain a column equals to at least one column name in dataframe A.
}
  \item{add}{
\code{add} the column name in dataframe B to be add to dataframe A. musth be character.
}
  \item{according}{
\code{according} the common column name specified to match the data entries from dataframe A and dataframe B. must be character. 
}
}
\details{
This function can be used to add one column from dataframe B to dataframe A, according to the column names specified. Users have to make sure the to dataframes at least share the common names specified. This function may be an alternative for \code{\link{merge}}
}


\author{
Jinlong Zhang
}

\seealso{
See Also\code{\link{merge}}
}
\examples{

data(splist)
data(testdata)
## add genera from dataframe B to dataframe A.
add.col(inputA = testdata, inputB = splist, add = "genera",
according = "species")
## add family from dataframe B to dataframe A.
add.col(inputA = testdata, inputB = splist, add = "family", 
according = "species")

}

\keyword{ add }

