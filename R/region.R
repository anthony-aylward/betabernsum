#===============================================================================
# region.R
#===============================================================================

# Imports ======================================================================

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
  size_sorted <- sort(size, decreasing = TRUE)
  parts <- diffparts(x)[1:length(size),, drop = FALSE]
  parts <- parts[,
    apply(
      parts,
      2,
      function(col) {
        all(col <= size_sorted) && (sum(col) >= x)
      }
    )
  ]
  parts <- do.call(
    cbind,
    lapply(
      apply(parts, 2, function(col) {unique(permn(col))}),
      function(l) {do.call(cbind, l)}
    )
  )
  parts <- parts[,
    apply(
      parts,
      2,
      function(col) {
        all(col <= size)
      }
    )
  ]
  t(parts)
}
