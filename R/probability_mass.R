#===============================================================================
# probability_mass.R
#===============================================================================

# Imports ======================================================================

#' @import VGAM




# Functions ====================================================================

#' @title Probability mass for a bbs in the independent case
#'
#' @param x vector giving the number of successes for each group.
#' @param size vector giving the number of trials for each group.
#' @param shape1,shape2 the two (positive) shape parameters of the standard
#'   beta distribution.
#' @return numeric, the probability mass.
probability_mass_independent <- function(
  x,
  size,
  prob = 0.5,
  rho = 0,
  shape1 = NULL,
  shape2 = NULL,
  ...
) {
  prod(
    sapply(
      1:length(x),
      function(index) {
        if (is.null(shape1) && is.null(shape2)) {
          dbetabinom(
            x[[index]],
            size[[index]],
            prob = prob[[index]],
            rho = rho[[index]],
            ...
          )
        } else if (is.numeric(shape1) && is.numeric(shape2)) {
          dbetabinom.ab(
            x[[index]],
            size[[index]],
            shape1 = shape1[[index]],
            shape2 = shape2[[index]],
            ...
          )
        } else if (is.numeric(shape1) && rho == 0) {
          function(index) {
            dbetabinom.ab(
              x[[index]],
              size[[index]],
              shape1 = shape1[[index]],
              shape2 = shape1[[index]]*(1-prob[[index]])/prob[[index]],
              ...
            )
          }
        } else {
          stop("Invalid arguments")
        }
      }
    )
  )
}

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
          1:length(x),
          function(i) {
            q = qbeta(t_i, shape1[[i]], shape2[[i]])
            if (x[[i]] == 0) {
              (1 - q)^(size[[i]])
            } else if (x[[i]] == size[[i]]) {
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
#' @return numeric, the value of the evaluated integral.
log_integral_dependent <- function(x, size, shape1, shape2) {
  integral <- tryCatch(
    integrate(
      integrand_dependent,
      0,
      1,
      x = x,
      size = size,
      shape1 = shape1,
      shape2 = shape2,
      rel.tol = 1e-15
    )[["value"]],
    error = function(e) NA
  )
  log(integral)
}

#' @title Logarithm of the coefficient component of probability mass
#'
#' @param x vector giving the number of successes for each group.
#' @param size vector giving the number of trials for each group.
#' @param shape1,shape2 the two (positive) shape parameters of the standard
#'   beta distribution. 
#' @return numeric, the value of the log-coefficient.
log_coefficient_dependent <- function(x, size) {
  sum(
    sapply(
      1:length(x),
      function(i) {
        log(choose(size[[i]], x[[i]]))
      }
    )
  )
}

#' @title Probability mass for a bbs in the dependent case
#'
#' @param x vector giving the number of successes for each group.
#' @param size vector giving the number of trials for each group.
#' @param shape1,shape2 the two (positive) shape parameters of the standard
#'   beta distribution. 
#' @return numeric, the probability mass.
probability_mass_dependent <- function(x, size, shape1, shape2) {
  exp(
    sum(
      log_coefficient_dependent(x, size),
      log_integral_dependent(x, size, shape1, shape2)
    )
  )
}
