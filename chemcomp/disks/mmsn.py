import astropy.constants as const
import astropy.units as u
import numpy as np

from ._0pacity import oplin

opacity_fct = oplin

from . import Disk

G = const.G.cgs.value
MS = const.M_sun.cgs.value
ME = const.M_earth.cgs.value
m_p = const.m_p.cgs.value
k_B = const.k_B.cgs.value
AU = const.au.cgs.value
year = 3600 * 24 * 364.25


class MMSN(Disk):
    """ docstring for MMSN. Hayashi 1981, Weidenschilling 1977. Likely outdated"""

    def __init__(self,
                 defaults,
                 chemistry,
                 M_STAR=1 * u.M_sun,
                 R0=80.0 * u.au,
                 ALPHA=1e-3,
                 SIGMA_POWER=-15 / 14,
                 ASPECT_POWER=2 / 7,
                 SIGMA_0=1500 * u.g / (u.cm ** 2),
                 ASPECT_0=0.033,
                 **kwargs):
        super().__init__(defaults, M_STAR, ALPHA, chemistry, **kwargs)
        self.sigma_init = SIGMA_0.cgs.value * (self.r / AU) ** SIGMA_POWER * np.exp(- self.r / R0.cgs.value)
        self.aspect_ratio = ASPECT_0 * (self.r / AU) ** ASPECT_POWER  # Hayashi 1981
        self._init_params(**kwargs)
        self.calc_opacity()

    def compute_disktmid(self, viscmodel=None):
        """calculate the disk temperature. """
        self.compute_cs_and_hp()
        self.compute_nu()
        self.calc_opacity()

    def calc_sigma_dust_static(self):
        """
        initializes the dust surface density

        This is equal to the pebble surface density, if DTG: peb has been set. Otherwise it contains the whole dust (including also the small grains)
        """
        vr = 2.0 * self.stokes_number_pebbles / (self.stokes_number_pebbles ** 2. + 1.0) * self.eta * self.v_k
        Pebbleflux = 10e-4 * ME / year
        self.sigma_dust = Pebbleflux / (2. * np.pi * self.r * vr)

        self.chemistry_solid, _ = self._chem.get_composition(self.T)

        self.sigma_dust_components = self.chemistry_solid * self.sigma_dust[:, np.newaxis,
                                                            np.newaxis]

        self.sigma_dust_components[:, 1, 0] = np.zeros_like(self.sigma_g)

        self.sigmin_peb = self.sigmin
        return

    def update_parameters(self):
        """ compute the viscous evolution. This function is called whenever the planet calls the update."""
        self.compute_viscous_evolution()

    def calc_opacity_oplin(self):
        """
        calculate the gas _0pacity - caution value in cgs
        metalicity is dependent of T! However, this is to good approximation constant -> not updated here.
        If T>1500K: No dust -> _0pacity = 0, fix by adding min value to metalicity
        """
        if hasattr(self, "rho_g"):
            self.opacity = 0.5 * (opacity_fct(self.rho_g, self.T))
