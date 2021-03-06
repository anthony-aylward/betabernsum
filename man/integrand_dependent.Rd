% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probability_mass.R
\name{integrand_dependent}
\alias{integrand_dependent}
\title{The integrand for computing probability mass in the dependent case}
\usage{
integrand_dependent(t, x, size, shape1, shape2)
}
\arguments{
\item{t}{vector of values}

\item{x}{vector giving the number of successes for each group.}

\item{size}{vector giving the number of trials for each group.}

\item{shape1, shape2}{the two (positive) shape parameters of the standard
beta distribution.}
}
\value{
numeric, the value of the evaluated integrand.
}
\description{
The integrand for computing probability mass in the dependent case
}
