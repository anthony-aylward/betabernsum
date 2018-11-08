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
#' @param independent logical. If TRUE (default), assume a sum of two
#'   independent groups of variables. If FALSE, assume all variables are
#'   mutually dependent.
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
  independent = TRUE,
  ...
) {
  if (is.null(shape1) && is.null(shape2)) {
    if (length(prob) == 1) prob <- rep(prob, 2)
    if (length(rho) == 1) rho <- rep(rho, 2)
    if (!all(sapply(list(size, prob, rho), length) == 2)) {
      stop("bad argument lengths")
    }
    if (!independent) {
      stop("provide `shape1` and `shape2` for the dependent case")
    }
  } else if (is.numeric(shape1) && is.numeric(shape2)) {
    if (length(shape1) == 1) shape1 <- rep(shape1, 2)
    if (length(shape2) == 1) shape2 <- rep(shape2, 2)
    if (!all(sapply(list(size, shape1, shape2), length) == 2)) {
      stop("bad argument lengths")
    }
  } else if (is.numeric(shape1) && rho == 0) {
    if (length(prob) == 1) prob <- rep(prob, 2)
    if (length(shape1) == 1) shape1 <- rep(shape1, 2)
    if (!all(sapply(list(size, prob, shape1), length) == 2)) {
      stop("bad argument lengths")
    }
    if (!independent) {
      stop("provide `shape1` and `shape2` for the dependent case")
    }
  }
  sapply(
    x,
    function(x_i) {
      sum(
        apply(
          region(x_i, size),
          1,
          function(row) {
            if (independent) {
              probability_mass_independent(
                row,
                size,
                prob = prob,
                rho = rho,
                shape1 = shape1,
                shape2 = shape2,
                ...
              )
            } else {
              probability_mass_dependent(
                row,
                size,
                shape1 = shape1,
                shape2 = shape2
              )
            }
          }
        )
      )
    }
  )
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
  lower_tail = TRUE,
  independent = TRUE,
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
              independent = independent,
              ...
            )
          }
        )
      )
    }
  )
}
