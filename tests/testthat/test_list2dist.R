context("list2dist")

test_that("convert data frame to distance matrix", {
  x <- matrix(rnorm(10), nrow = 5)
  
  # Compute distance matrix directly.
  #
  y2 = dist(x, upper = TRUE)

  # Compute distances in data frame.
  #
  x.dist = lapply(1:5, function(i) {
    lapply(1:5, function(j) {
      orig = x[i,]
      dest = x[j,]
      c(i, j, sqrt(sum((orig - dest)**2)))
    })
  })
  x.dist = data.frame(do.call(rbind, unlist(x.dist, recursive = FALSE)))
  #
  # Convert to distance matrix.
  #
  y1 = list2dist(x.dist)

  expect_equivalent(y1, y2)
})