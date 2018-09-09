#===============================================================================
# bbs.R
#===============================================================================

#' @title Probability mass function for a sum of beta-bernoulli variables
#'
#' @param x,q vector of quantiles.
#' @param size vector giving the number of trials for each group.
#' @param prob vector giving probability of success for each group.
#' @param rho vector giving correlation parameter for each group.
#' @param shape1,shape2 the two (positive) shape parameters of the standard
#'   beta distribution. See the documentation for \code{Betabinom} in the
#'   \code{VGAM} package.
#' @param lower_tail logical. If TRUE (default), probabilities are
#'   \eqn{P[X <= x]} otherwise, \eqn{P[X > x]}.
#' @param ... other parameters passed to \code{Betabinom}
#' @return numeric, the value of the probability mass or distribution function.
#' @export
dbbs <- function(
  x,
  size,
  prob = 0.5,
  rho = 0,
  shape1 = NULL,
  shape2 = NULL,
  ...
) {
  if (is.null(shape1) && is.null(shape2)) {
    pmf_of_sum_prob_rho(x, size, prob = prob, rho = rho, ...)
  } else if (is.numeric(shape1) && is.numeric(shape2)) {
    pmf_of_sum_shape1_shape2(x, size, shape1 = shape1, shape2 = shape2, ...)
  } else if (is.numeric(shape1) && rho == 0) {
    pmf_of_sum_prob_shape1(x, size, prob = prob, shape1 = shape1, ...)
  } else {
    stop("Invalid arguments")
  }
}

#' @describeIn dbbs Distribution function for a sum of beta-bernoulli variables
#' @export
pbbs <- function(
  q,
  size,
  prob = 0.5,
  rho = 0,
  shape1 = NULL,
  shape2 = NULL,
  lower_tail = TRUE
  ...
) {
  sapply(
    q,
    function(q_i) {
      if (lower_tail) {
        if (q_i < 0) return(0)
        tail <- 0:q_i
      } else if (q_i >= sum(size)) {
        return(0)
      } else {
        tail <- (q_i + 1):sum(size)
      }
      sum(
        sapply(
          tail,
          function(x) {
            dbbs(
              x,
              size = size,
              prob = prob,
              rho = rho,
              shape1 = shape1,
              shape2 = shape2,
              ...
            )
          }
        )
      )
    }
  )
}