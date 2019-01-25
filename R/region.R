#===============================================================================
# region.R
#===============================================================================

# Imports ======================================================================

#' @import combinat
#' @import partitions




# Functions ====================================================================

#' @title Region of integration for bernoulli sums
#'
#' @param x number of successes
#' @param size vector giving the number of trials for each group.
#' @return a matrix with one column per sample, the rows of which are the
#'   coordinates of the region of integration.
region <- function(x, size) {
  if (x < 0 || x > sum(size)) stop("provided x is out of bounds")
  if (x == 0) return(matrix(0, nrow = 1, ncol = length(size)))
  if (x == sum(size)) return(matrix(size, nrow = 1))
  size_sorted <- sort(size, decreasing = TRUE)

  if (x == 1) {
    part <- matrix(c(1, rep(0, length(size) - 1)), ncol = 1)
  } else {
    if (x < length(size)) {
      print("x < length(size)")
      print(x)
      part <- parts(x)
      part <- rbind(
        part,
        matrix(0, nrow = length(size) - x, ncol = ncol(part))
      )
    } else {
      print("x >= length(size)")
      print(x)
      part <- parts(x)[1:length(size),, drop = FALSE]
    }
    part <- part[,
      apply(part, 2, function(col) all(col <= size_sorted) && (sum(col) >= x))
    ]
  }

  part <- do.call(
    cbind,
    lapply(
      apply(part, 2, function(col) {unique(permn(col))}),
      function(l) {do.call(cbind, l)}
    )
  )
  part <- part[,
    apply(
      part,
      2,
      function(col) {
        all(col <= size)
      }
    )
  ]
  t(part)
}
