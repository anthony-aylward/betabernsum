#===============================================================================
# probability_mass.py
#===============================================================================

# Imports ======================================================================

from functools import reduce
from math import log
from operator import mul
from scipy.integrate import quad
from scipy.misc import comb
from scipy.special import beta as beta_function
from scipy.stats import beta




# Functions ====================================================================

def prod(iterable):
    """"Product of values in an iterable
    
    Parameters
    ----------
    iterable
        an iterable
    
    Returns
    -------
    the product
    """
    return reduce(mul, iterable, 1)


def betabinom_pmf(k, n, a, b):
    """Probability mass function for a beta binomial distribution

    Parameters
    ----------
    k : int
        number of successes
    n : int
        number of trials
    a : float
        first (positive) shape parameter
    b : float
        second (positive) shape parameter

    Returns
    -------
    float
        the probability mass
    """

    return comb(n, k) * beta_function(k + a, n - k + b) / beta_function(a, b)


def probability_mass_independent(k, n, a, b):
    """Probability mass for a bbs in the independent case

    Parameters
    ----------
    k
        iterable giving the number of successes for each group
    n
        iterable giving the number of trials for each group
    a
        iterable giving the first shape parameter for each group
    b
        iterable giving the second shape parameter for each group
    
    Returns
    -------
    float
        the probability mass
    """
    
    assert len(k) == len(n) == len(a) == len(b)
    return prod(
        betabinom_pmf(k_i, n_i, a_i, b_i)
        for k_i, n_i, a_i, b_i in zip(k, n, a, b)
    )


def evaluate_factor(t, k, n, a, b):
    """Evaluate a factor for the product in `integrand_dependent`

    t : float
        a value on the domain of integration
    k : int
        a number of successes
    n : int
        a number of trials
    a : float
        a first shape parameter
    b : float
        a second shape parameter
    
    Returns
    -------
    float
        value of a factor
    """

    q = beta.ppf(t, a, b)
    if k == 0:
        return (1 - q)**n
    elif k == n:
        return q**k
    else:
        return q**k * (1 - q)**(n - k)
    

def integrand_dependent(t, k, n, a, b):
    """The integrand for computing probability mass in the dependent case

    Parameters
    ----------
    t
        iterable of values on the domain of integration
    k
        iterable giving the number of successes for each group
    n
        iterable giving the number of trials for each group
    a
        iterable giving the first shape parameter for each group
    b
        iterable giving the second shape parameter for each group

    Returns
    -------
    float
    """

    assert len(k) == len(n) == len(a) == len(b)
    return tuple(
        prod(
            evaluate_factor(t_i, k_j, n_j, a_j, b_j)
            for k_j, n_j, a_j, b_j in zip(k, n, a, b)
        )
        for t_i in t
    )


def log_integral_dependent(k, n, a, b):
    """Get logarithm of integral
    
    Parameters
    ----------
    k
    n
    a
    b
    """
    
    return log(quad(integrand_dependent, 0, 1, args = (k, n, a, b))[0])
