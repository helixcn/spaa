#' Convert data frame to distance matrix.
#'
#' @param dat A data frame with (at least) three columns.
#'
#' @return A distance matrix.
#' @export
#'
#' @examples
list2dist <-
function (dat) 
{
  dat.name1 <- as.character(dat[, 1])
  dat.name2 <- as.character(dat[, 2])
  dat.value <- dat[, 3]
  #
  total.names = union(dat.name1, dat.name2)
  #
  elements <- matrix(nrow = length(total.names), ncol = length(total.names))
  #
  # Set diagonal elements to zero.
  #
  diag(elements) <- 0
  #
  rownames(elements) <- total.names
  colnames(elements) <- total.names
  #
  rows = factor(dat.name1, levels = total.names)
  cols = factor(dat.name2, levels = total.names)
  #
  for (n in 1:nrow(dat)) {
    irow = as.integer(rows[n])
    jcol = as.integer(cols[n])
    elements[irow, jcol] = dat.value[n]
  }
  print(elements[1:5,1:5])
  #
  as.dist(t(elements), upper = TRUE)
}
