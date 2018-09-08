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
#'   \item{x}{the number of successes in the input.}
#'   \item{size}{the sample size parameter of the input.}
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
  dist_val <- pbbs(
    q = x,
    size = size,
    prob = prob,
    rho = rho,
    shape1 = shape1,
    shape2 = shape2
  )
  if (alternative[[1]] == "two_sided") {
    p_value <- min(1, 2 * min(dist_val, 1 - dist_val))
  } else if (alternative[[1]] == "less") {
    p_value <- dist_val
  } else if (alternative[[1]] == "greater") {
    p_value <- 1 - dist_val
  }
  list(
    statistic = x,
    parameter = size,
    p_value = p_value
  )
}
