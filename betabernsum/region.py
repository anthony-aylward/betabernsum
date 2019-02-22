#===============================================================================
# region.R
#===============================================================================

# Imports ======================================================================

import itertools
from accelasc import accel_asc




# Functions ====================================================================

def region(x, size):
    """Region of integration for bernoulli sums

    Parameters
    ----------
    x : int
        number of successes
    size
        iterable of integers, giving the number of trials for each group
    
    Returns
    -------
    numpy array
        an array with one column per sample, the rows of which are the
        coordinates of the region of integration.
    """

    if x < 0 or x > sum(size):
        raise RuntimeError('provided x is out of bounds')
    if x == 0:
        yield (0,) * len(size)
        return
    if x == sum(size):
        yield tuple(size)
        return
    yield from itertools.chain.from_iterable(
        set(
            permutation for permutation in itertools.permutations(
                partition + [0] * (len(size) - len(partition))
            )
            if all((p <= s for p, s in zip(permutation, size)))
        )
        for partition in accel_asc(x) if len(partition) <= len(size)
    )
