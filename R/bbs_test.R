#===============================================================================
# bbs_test.R
#===============================================================================

#' @title Beta Bernoulli Sum Test
#'
#' @description Performs a hypothesis test against a BBS null.
#'
#' @param x number of successes
#' @param size vector giving the number of trials for each group.
#' @param prob vector giving probability of success for each group.
#' @param rho vector giving correlation parameter for each group.
#' @param shape1,shape2 the two (positive) shape parameters of the standard
#'   beta distribution. See the documentation for \code{Betabinom} in the
#'   \code{VGAM} package.
#' @param ... other parameters passed to \code{Betabinom}
#' @return \describe{
#'   \item{statistic}{the number of successes in the input.}
#'   \item{parameter}{the number of trials for each group of the input.}
#'   \item{p_value}{the p-value of the hypothesis test.}
#' }
#' @export
bbs_test <- function(
  x,
  size,
  prob = 0.5,
  rho = 0,
  shape1 = NULL,
  shape2 = NULL,
  alternative = c("two_sided", "less", "greater"),
  ...
) {
  lower_tail_area <- pbbs(
    q = x,
    size = size,
    prob = prob,
    rho = rho,
    shape1 = shape1,
    shape2 = shape2
  )
  upper_tail_area <- pbbs(
    q = x - 1,
    size = size,
    prob = prob,
    rho = rho,
    shape1 = shape1,
    shape2 = shape2,
    lower_tail = FALSE
  )
  if (alternative[[1]] == "two_sided") {
    if (lower_tail_area < upper_tail_area) {
      p_value <- lower_tail_area + pbbs(
        q = sum(size) - x - 1,
        size = size,
        prob = prob,
        rho = rho,
        shape1 = shape1,
        shape2 = shape2,
        lower_tail = FALSE
      )
    } else {
      p_value <- upper_tail_area + pbbs(
        q = sum(size) - x,
        size = size,
        prob = prob,
        rho = rho,
        shape1 = shape1,
        shape2 = shape2
      )
    }
  } else if (alternative[[1]] == "less") {
    p_value <- lower_tail_area
  } else if (alternative[[1]] == "greater") {
    p_value <- upper_tail_area
  }
  list(statistic = x, parameter = size, p_value = p_value)
}