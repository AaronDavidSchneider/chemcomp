import numpy as np


class Accretion(object):
    """Base class of the Accretion model."""

    def __init__(self, config=None):
        self.m_dot_a = 0
        self.m_dot_c = 0
        self.m_dot = self.m_dot_a + self.m_dot_c
        self.M_z = np.zeros(2)
        self.name = "DefaultAccretion"
        self.regime = ""

    def calc_m_dot(self):
        """function that calculates and updates the accretion model"""
        pass

    def accretion_domain(self):
        """function that calculates the regimes of the accretion model"""
        pass

    def init_accretion(self, planet):
        """
        only passes planet and chemistry.
        Is called in __init__ of planet. At the end of init
        """
        self.planet = planet
        self.chem_mask_array = self.planet.disk._chem.mask_array
        self.m_dot_a_chem = self.chem_mask_array * self.m_dot_a
        self.m_dot_c_chem = self.chem_mask_array * self.m_dot_c
        self._init_params()

    def _init_params(self):
        """function that inits parameters needed"""
        pass

    def update_z(self):
        """update the total mass of heavy elements that have been accreted."""
        self.M_z[0] += np.sum(self.m_dot_c_chem[0, :-7] * self.planet.h)
        self.M_z[1] += np.sum(self.m_dot_a_chem[0, :-7] * self.planet.h)

    def remove_mass_from_flux(self):
        pass
