#===============================================================================
# beta_bern_sum.R
#===============================================================================

# Imports ======================================================================

#' @import VGAM




# Functions ====================================================================

#' @title Region of integration for bernoulli sums
#'
#' @param x number of successes
#' @param size vector giving the number of trials for each group.
region <- function(x, size) {
  if (x < 0 || x > sum(size)) stop("provided x is out of bounds")
  coord1 <- min(0, x - size[[2]]):min(x, size[[1]])
  cbind(coord1, x - coord1)
}

#' @title Density of a sum with \code{prob}/\code{rho} parametrization
#'
#' @param x number of successes
#' @param size vector giving the number of trials for each group.
#' @param prob vector giving probability of success for each group.
#' @param rho vector giving correlation parameter for each group.
#' @param ... other parameters passed to dbetabinom
density_of_sum_prob_rho <- function(x, size, prob = 0.5, rho = 0, ...) {
  if (length(prob) == 1) prob <- rep(prob, 2)
  if (length(rho) == 1) rho <- rep(rho, 2)
  if (!all(sapply(list(size, prob, rho), length) == 2)) {
    stop("bad argument lengths")
  }
  sum(
    apply(
      region(x, size),
      1,
      function(row) {
        prod(
          sapply(
            c(1, 2),
            function(index) {
              dbetabinom(
                row[[index]],
                size[[index]],
                prob = prob[[index]],
                rho = rho[[index]],
                ...
              )
            }
          )
        )
      }
    )
  )
}

#' @title Density of a sum with \code{shape1}/\code{shape2} parametrization
#'
#' @param x number of successes
#' @param size vector giving the number of trials for each group.
#' @param shape1,shape2 the two (positive) shape parameters of the standard
#'   beta distribution. See the documentation for \code{Betabinom} in the
#'   \code{VGAM} package.
#' @param ... other parameters passed to dbetabinom.ab
density_of_sum_shape1_shape2 <- function(x, size, shape1, shape2, ...) {
  if (length(shape1) == 1) shape1 <- rep(shape1, 2)
  if (length(shape2) == 1) shape2 <- rep(shape2, 2)
  if (!all(sapply(list(size, shape1, shape2), length) == 2)) {
    stop("bad argument lengths")
  }
  sum(
    apply(
      region(x, size),
      1,
      function(row) {
        prod(
          sapply(
            c(1,2),
            function(index) {
              dbetabinom.ab(
                row[[index]],
                size[[index]],
                shape1 = shape1[[index]],
                shape2 = shape2[[index]],
                ...
              )
            }
          )
        )
      }
    )
  )
}

#' @title Density of a sum with \code{prob}/\code{shape1} parametrization
#'
#' @param x number of successes
#' @param size vector giving the number of trials for each group.
#' @param prob vector giving probability of success for each group.
#' @param shape1 vector giving correlation parameter for each group.
#' @param ... other parameters passed to dbetabinom.ab
density_of_sum_prob_shape1 <- function(x, size, prob, shape1, ...) {
  if (length(prob) == 1) prob <- rep(prob, 2)
  if (length(shape1) == 1) shape1 <- rep(shape1, 2)
  if (!all(sapply(list(size, prob, shape1), length) == 2)) {
    stop("bad argument lengths")
  }
  sum(
    apply(
      region(x, size),
      1,
      function(row) {
        prod(
          sapply(
            c(1,2),
            function(index) {
              dbetabinom.ab(
                row[[index]],
                size[[index]],
                shape1 = shape1[[index]],
                shape2 = shape1[[index]] * (1 - prob[[index]]) / prob[[index]],
                ...
              )
            }
          )
        )
      }
    )
  )
}

#' @title Probability mass function for a sum of beta-bernoulli variables
#'
#' @param x vector of quantiles.
#' @param size vector giving the number of trials for each group.
#' @param prob vector giving probability of success for each group.
#' @param rho vector giving correlation parameter for each group.
#' @param shape1,shape2 the two (positive) shape parameters of the standard
#'   beta distribution. See the documentation for \code{Betabinom} in the
#'   \code{VGAM} package.
#' @param ... other parameters passed to \code{Betabinom}
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
    density_of_sum_prob_rho(x, size, prob = prob, rho = rho, ...)
  } else if (is.numeric(shape1) && is.numeric(shape2)) {
    density_of_sum_shape1_shape2(x, size, shape1 = shape1, shape2 = shape2, ...)
  } else if (is.numeric(shape1) && rho == 0) {
    density_of_sum_prob_shape1(x, size, prob = prob, shape1 = shape1, ...)
  } else {
    stop("Invalid arguments")
  }
}
