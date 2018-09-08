#===============================================================================
# region.R
#===============================================================================

#' @title Region of integration for bernoulli sums
#'
#' @param x number of successes
#' @param size vector giving the number of trials for each group.
#' @return a matrix with two columns, the rows of which are the coordinates of
#'   the region of integration.
region <- function(x, size) {
  if (x < 0 || x > sum(size)) stop("provided x is out of bounds")
  coord1 <- min(0, x - size[[2]]):min(x, size[[1]])
  cbind(coord1, x - coord1)
}
