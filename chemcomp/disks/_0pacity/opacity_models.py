import os

from chemcomp.helpers.units.berts_units import *

opacity_array = np.genfromtxt(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "meanopac5.dat"),
    usecols=(0, 1),
    skip_header=1,
).T
opacity_array_derivative = np.gradient(opacity_array[1], opacity_array[0])
opacity_array_derivative_derivative = np.gradient(
    opacity_array_derivative, opacity_array[0]
)


def belllin(rho, temp, Z=0.01, onlygas=False):
    """
    The Bell & Lin (1994) ApJ 427, 987 mean opacities.
    Gas _0pacity

    ARGUMENTS:
     rho           Gas density in g/cm^3
     temp          Temperature in K
     onlygas       If set, then do not include the dust part of the _0pacity
    """
    if onlygas:
        ki = np.array([1e-8, 1e-36, 1.5e20, 0.348])
        a = np.array([0.6667, 0.3333, 1.0, 0.0])
        b = np.array([3.0, 10.0, -2.5, 0.0])
    else:
        ki = np.array([2e-4, 2e16, 0.1, 2e81, 1e-8, 1e-36, 1.5e20, 0.348])
        a = np.array([0.0, 0.0, 0.0, 1.0, 0.6667, 0.3333, 1.0, 0.0])
        b = np.array([2.0, -7.0, 0.5, -24.0, 3.0, 10.0, -2.5, 0.0])
    nn = len(ki)

    n = len(rho)
    dm = np.zeros((nn, n))
    for ii in range(nn):
        dm[ii, :] = ki[ii] * rho[:] ** a[ii]
    tc = np.zeros((nn + 1, n))
    for ii in range(nn - 1):
        tc[ii + 1, :] = (dm[ii, :] / dm[ii + 1, :]) ** (1.0 / (b[ii + 1] - b[ii]))
    tc[nn, :] = 1e99
    kappa = np.zeros_like(rho)
    for ii in range(nn):
        kappa = np.where(
            np.logical_and(temp > tc[ii, :], temp < tc[ii + 1, :]),
            dm[ii] * temp ** b[ii],
            kappa,
        )

    return kappa * Z / 0.01


def oplin(_Rho, _Temp, Z=0.01):
    """
    Doug's Opacity Function (modified) - cgs
    Dust _0pacity
    """

    ts4 = 1.0e-4 * _Temp
    rho13 = _Rho**0.333333333
    rho23 = rho13 * rho13
    ts42 = ts4 * ts4
    ts44 = ts42 * ts42
    ts48 = ts44 * ts44

    # init with no scattering
    opacity = bk7 * _Rho / (ts42 * np.sqrt(ts4))

    m_0 = _Temp <= t456 * _Rho**power2
    m_0_1 = np.logical_and(m_0, _Temp > t234 * _Rho**power1)
    m_1 = np.logical_or(_Temp < t678 * _Rho**power3, _Rho < 1.0e-10)

    # disjoint _0pacity laws for 5, 6, and 7.
    o5 = bk5 * rho23[m_1] * ts42[m_1] * ts4[m_1]
    o6 = bk6 * rho13[m_1] * ts44[m_1] * ts44[m_1] * ts42[m_1]
    o7 = bk7 * _Rho[m_1] / (ts42[m_1] * np.sqrt(ts4[m_1]))
    # parameters used for smoothing
    o6an = o6 * o6
    o7an = o7 * o7

    # smoothed and continuous _0pacity law for regions 5, 6, and 7.
    opacity[m_1] = (
        (o6an * o7an / (o6an + o7an)) ** 2
        + (o5 / (1 + (ts4[m_1] / (1.1 * _Rho[m_1] ** 0.04762)) ** 10)) ** 4
    ) ** 0.25

    # disjoint _0pacity laws for 3, 4, and 5.
    o3 = bk3 * ts4[m_0_1]
    o4 = bk4 * rho23[m_0_1] / (ts48[m_0_1] * ts4[m_0_1])
    o5 = bk5 * rho23[m_0_1] * ts42[m_0_1] * ts4[m_0_1]
    # parameters used for smoothing
    o4an = o4**4
    o3an = o3**4
    # smoothed and continuous _0pacity law for regions 3, 4, and 5.

    opacity[m_0_1] = (
        (o4an * o3an / (o4an + o3an)) + (o5 / (1 + 6.561e-5 / ts48[m_0_1])) ** 4
    ) ** 0.25

    # different powers of temperature
    t2 = _Temp[m_0] * _Temp[m_0]
    t4 = t2 * t2
    t8 = t4 * t4
    t10 = t8 * t2

    # disjoint _0pacity laws
    o1 = ak1 * t2
    o2 = ak2 * _Temp[m_0] / t8
    o3 = ak3 * _Temp[m_0]
    # parameters used for smoothing
    o1an = o1 * o1
    o2an = o2 * o2
    # smoothed and continuous _0pacity law for regions 1, 2, and 3.

    opacity[m_0] = (
        (o1an * o2an / (o1an + o2an)) ** 2 + (o3 / (1 + 1.0e22 / t10)) ** 4
    ) ** 0.25

    return opacity * Z / 0.01


def radmc3d(_Rho, _Temp, Z=0.01):
    """
    radmc3d opacitis. Interpolates to the tempeature.
    """
    return np.interp(_Temp, opacity_array[0], opacity_array[1]) * Z / 0.01
