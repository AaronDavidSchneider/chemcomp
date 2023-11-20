import numpy as np

epsh = 0.4  #   Smoothing in units of H

XMSOL = 1.989e33  # Solar Mass
XMH = 1.67e-24  # Hydrogen Mass
AU = 1.496e13  # Astronomical Unit
BOLTZ = 1.38e-16  # Boltzmann Constant
CLIGHT = 2.998e10  # Light speed
GRAVC = 6.67e-08  # Gravitation Constant
AR = 7.56e-15  # Radiation Constant
RGAS = 8.314e07  # Gas Constant
sr = 5.670e-5  # Stefan Boltzmann
ME = 5.97e27  # Earth mass in gramms
year = 3600.0 * 24.0 * 365.24

pi = np.pi

# Doug's Opacity Function (modified)
power1 = 4.44444444e-2
power2 = 2.381e-2
power3 = 2.267e-1
# data t234,t456,t678/1.6e3, 5.7e3, 2.28e6/
t234 = 1.6e3
t456 = 5.7e3
t678 = 2.28e6
# coefficients for _0pacity laws 1, 2, and 3 in cgs units.
# data ak1,ak2,ak3/2.e-4, 2.e16, 5.e-3/
ak1 = 2.0e-4
ak2 = 2.0e16
ak3 = 5.0e-3
# coefficients for _0pacity laws 3, 4, 5, 6, 7, and 8 in T_4 units.
# data bk3,bk4,bk5,bk6,bk7
# & /50., 2.e-2, 2.e4, 1.d4, 1.5d10/
bk3 = 50.0
bk4 = 2.0e-2
bk5 = 2.0e4
bk6 = 1.0e4
bk7 = 1.5e10
# test T against (T_23 * T_34 * T_34)**0.333333333
