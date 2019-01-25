#===============================================================================
# bbs.R
#===============================================================================

# Imports ======================================================================

#' @import parallel




# Functions ====================================================================

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
#' @param cores integer. Number of cores to use (default 1).
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
    if (!(length(unique(sapply(list(size, prob, rho), length))) == 1)) {
      stop("bad argument lengths")
    }
    if (!independent) {
      stop("provide `shape1` and `shape2` for the dependent case")
    }
  } else if (is.numeric(shape1) && is.numeric(shape2)) {
    if (!(length(unique(sapply(list(size, shape1, shape2), length))) == 1))  {
      stop("bad argument lengths")
    }
  } else if (is.numeric(shape1) && rho == 0) {
    if (!(length(unique(sapply(list(size, prob, shape1), length))) == 1)) {
      stop("bad argument lengths")
    }
    if (!independent) {
      stop("provide `shape1` and `shape2` for the dependent case")
    }
  }
  sapply(
    x,
    function(x_i) {
      reg <- region(x_i, size)
      summand <- apply(
        reg,
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
      sum(summand)
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
  independent = TRUE,
  lower_tail = TRUE,
  cores = 1,
  ...
) {
  sapply(
    q,
    function(q_i) {
      if (lower_tail) {
        if (q_i < 0) return(0)
        if (q_i <= (sum(size) / 2)) {
          tail <- 0:q_i
          speed_flip <- FALSE
        } else {
          tail <- (q_i + 1):sum(size)
          speed_flip <- TRUE
        }
      } else if (q_i >= sum(size)) {
        return(0)
      } else {
        if (q_i >= (sum(size) / 2)) {
          tail <- (q_i + 1):sum(size)
          speed_flip <- FALSE
        } else {
          tail <- 0:q_i
          speed_flip <- TRUE
        }
      }
      if (any(tail >= 84)) return(NA)
      as.integer(speed_flip) + (1 - 2 * speed_flip) * sum(
        unlist(
          mclapply(
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
            },
            mc.cores = cores
          )
        )
      )
    }
  )
}
