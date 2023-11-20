import astropy.constants as const
import astropy.units as u
import numpy as np

from . import Accretion

earth_mass_cgs_value = 1 * const.M_earth.cgs.value
Myr = 1e6 * 3600 * 24 * 365.25
year = 3600 * 24 * 365.25


class GasAccretion(Accretion):
    """ Model for GasAccretion.
    Minimum of Ikoma2000, Machida2010 and disk gas accretion rate (including horseshoe region).
    """

    def __init__(self, config):
        super().__init__()
        self.name = "GasAccretion"
        self.kappa_env = config.get("kappa_env", 0.05 * u.cm ** 2 / u.g)
        self.kappa_env = self.kappa_env.cgs.value
        self.f = config.get("f", 0.2)
        self.f_machida = config.get("f_machida", 1)
        self.f_disk_max = config.get("f_disk_max", 1)

    def _init_params(self):
        self.const_terms = (
                1e3
                * (30*earth_mass_cgs_value)**2.5
                * self.kappa_env/0.05
                * year
        )

    def accretion_domain(self):
        pass

    def update_parameters(self):
        pass

    def remove_mass_from_flux(self):
        return

    def get_min(self):
        """" Get the minimal accretion rate """
        low = (
                0.83
                * self.f_machida
                * self.planet.omega_k
                * self.planet.H_g ** 2
                * self.planet.sigma_g
                * (self.planet.r_h / self.planet.H_g) ** (9 / 2)
        )
        low = np.maximum(low, 1e-60)

        high = (
                0.14
                * self.f_machida
                * self.planet.omega_k
                * self.planet.H_g ** 2
                * self.planet.sigma_g
        )
        high = np.maximum(high, 1e-60)

        disk_max = np.maximum(self.f_disk_max * self.planet.m_dot_disk, 1e-60)
        if self.planet.fSigma < 1:
            dOmegaHS = 1.5 * np.pi * self.planet.omega_k * self.planet._xs
            T_HS = 2 * np.pi / dOmegaHS
            mdot_HS = self.planet._M_HS * self.planet.fSigma / (2 * T_HS)
            mdot_HS = np.maximum(mdot_HS, 1e-60)
            disk_max = disk_max + mdot_HS

        ikoma = self.planet.M/(self.planet.M_c**(-2.5)*self.const_terms)
        ikoma = np.maximum(ikoma, 1e-60)

        valid_accrates = [low, high, disk_max, ikoma]
        minimum = np.argmin(valid_accrates)
        mdot = valid_accrates[minimum]

        assert mdot >= 0, "negative mdot!"
        return str(minimum), mdot

    def calc_m_dot(self):
        """ main routine that calculates the accretion rate."""
        self.update_parameters()
        self.accretion_domain()
        self.m_dot_c = 0
        self.m_dot_c_chem = self.chem_mask_array * self.m_dot_c

        if self.planet.past_pebble_iso:
            self.regime, self.m_dot_a = self.get_min()
            self.m_dot_a_chem = self.planet.chemistry_gas * self.m_dot_a
            self.m_dot = self.m_dot_a + self.m_dot_c

        else:
            self.m_dot_a = 0
            self.m_dot = self.m_dot_a + self.m_dot_c
            self.m_dot_a_chem = self.chem_mask_array * self.m_dot_a


        assert np.all(self.m_dot_a_chem >= 0)
        return