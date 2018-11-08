#===============================================================================
# probability_mass_dependent.R
#===============================================================================

# Imports ======================================================================




# Functions ====================================================================

#' @title The integrand for computing probability mass in the dependent case
#'
#' @param t vector of values
#' @param x vector giving the number of successes for each group.
#' @param size vector giving the number of trials for each group.
#' @param shape1,shape2 the two (positive) shape parameters of the standard
#'   beta distribution. 
#' @return numeric, the value of the evaluated integrand.
integrand_dependent <- function(t, x, size, shape1, shape2) {
  sapply(
    t,
    function(t_i) {
      prod(
        sapply(
          1:min(sapply(list(x, size, shape1, shape2), length)),
          function(i) {
            q = qbeta(t_i, shape1[[i]], shape2[[i]])
            if (x[[i]] == 0) {
              (1 - q)^(size[[i]])
            else if (x[[i]] == size[[i]]) {
              q^(x[[i]])
            } else {
              q^x[[i]] * (1 - q)^(size[[i]] - x[[i]])
            }
          }
        )
      )
    }
  )
}

#' @title Logarithm of the integral component of probability mass
#'
#' @param x vector giving the number of successes for each group.
#' @param size vector giving the number of trials for each group.
#' @param shape1,shape2 the two (positive) shape parameters of the standard
#'   beta distribution. 
#' @return numeric, the value of the evaluated integrand.
log_integral_dependent <- function(x, size, shape1, shape2) {
  log(
    integrate(
      integrand,
      0,
      1,
      x = x,
      size = size,
      shape1 = shape1,
      shape2 = shape2
    )[["value"]]
  )
}

#' @title Logarithm of the coefficient component of probability mass
#'
#' @param x vector giving the number of successes for each group.
#' @param size vector giving the number of trials for each group.
#' @param shape1,shape2 the two (positive) shape parameters of the standard
#'   beta distribution. 
#' @return numeric, the value of the log-coefficient.
log_coefficient_dependent <- function(x, size, shape1, shape2) {
  sum(
    sapply(
      1:min(sapply(list(x, size, shape1, shape2), length)),
      function(i) {
        log(choose(size[[i]], x[[i]]))
      }
    )
  )
}

#' @title Probability mass for a certain outcome
#'
#' @param x vector giving the number of successes for each group.
#' @param size vector giving the number of trials for each group.
#' @param shape1,shape2 the two (positive) shape parameters of the standard
#'   beta distribution. 
#' @return numeric, the value of the log-coefficient.
probability_mass_dependent <- function(x, size, shape1, shape2) {
  exp(
    sum(
      log_coefficient(x, size, shape1, shape2),
      log_integral(x, size, shape1, shape2)
    )
  )
}

#' @title Probability mass function of a sum with \code{shape1}/\code{shape2}
#'   parametrization
#'
#' @param x vector of quantiles
#' @param size vector giving the number of trials for each group.
#' @param shape1,shape2 the two (positive) shape parameters of the standard
#'   beta distribution. See the documentation for \code{Betabinom} in the
#'   \code{VGAM} package.
#' @return numeric, the value of the BBS probability mass function.
pmf_of_sum_shape1_shape2_dependent <- function(x, size, shape1, shape2) {
  if (length(shape1) == 1) shape1 <- rep(shape1, 2)
  if (length(shape2) == 1) shape2 <- rep(shape2, 2)
  if (!all(sapply(list(size, shape1, shape2), length) == 2)) {
    stop("bad argument lengths")
  }
  sapply(
    x,
    function(x_i) {
      sum(
        apply(
          region(x_i, size),
          1,
          function(row) {
            probability_mass(row, size, shape1, shape2)
          }
        )
      )
    }
  )
}
