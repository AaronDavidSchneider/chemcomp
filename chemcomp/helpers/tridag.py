import numpy as np
from scipy.linalg import solve_banded

def tridag_components(a, b, c, r):
    """
    The Numerical Recipes tridag() routine for Python, using the SciPy library. 
    This is just a wrapper around the solve_banded function of SciPy, to make it 
    work like the tridag routine.
    """
    n = r.shape
    ab = np.zeros((3, *n))
    ab[0, 1:] = c[:-1]
    ab[1, :] = b[:]
    ab[2, :-1] = a[1:]

    u = np.transpose(
        [solve_banded((1, 1), ab[:, :, i], r[:, i]) for i in range(n[1])])
    return u


def tridag(a, b, c, r):
    """
    The Numerical Recipes tridag() routine for Python, using the SciPy library.
    This is just a wrapper around the solve_banded function of SciPy, to make it
    work like the tridag routine.
    """
    n = r.shape
    ab = np.zeros((3, *n))
    ab[0, 1:] = c[:-1]
    ab[1, :] = b[:]
    ab[2, :-1] = a[1:]

    u = solve_banded((1, 1), ab, r)
    return u


def pentdag(aa, a, b, c, cc, r):
    """
    A 5-diagonal version of tridag.
    """
    n = r.size
    ab = np.zeros((5, n))
    ab[0, 2:] = cc[:-2]
    ab[1, 1:] = c[:-1]
    ab[2, :] = b[:]
    ab[3,:-1] = a[1:]
    ab[4,:-2] = aa[2:]
    u = solve_banded((2,2),ab,r)
    return u
    
