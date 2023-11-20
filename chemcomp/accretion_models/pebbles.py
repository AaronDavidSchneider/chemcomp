import astropy.constants as const
import astropy.units as u
import numpy as np

from . import Accretion

G = 1 * const.G.cgs.value


class PebbleAccretion(Accretion):
    """ Model for PebbleAccretion. """

    def __init__(self, config):
        super().__init__()
        self.name = "PebbleAccretion"
        self.epsilon_p = config.get("epsilon_p", 0.5)
        self.rho_solid = config.get("rho_solid")
        self.u_frag = config.get("u_frag", None)
        self.t_f = config.get("STOKES", None)

        self.factor_peb_accretion = config.get('factor_peb_accretion', 1.0)

        if self.rho_solid is None:
            self._calculate_rho_solid = True
        else:
            self.rho_solid = self.rho_solid.cgs.value
            self._calculate_rho_solid = False

        if self.u_frag is None:
            self._calculate_u_frag = True
            print("Warning: No fixed fragmentation velocity has been set! The code will calculate the fragmentation velocity according to the Drazkowska 2017 approach: 1m/s if silicate, 10m/s if icy.")
        else:
            self.u_frag = self.u_frag.cgs.value
            self._calculate_u_frag = False

        self.regime = "Bondi"
        self.twoD = False

    def _init_params(self):
        """ Pass some pebble parameters to the disk. """
        self.planet.disk.init_pebble_parameters(self.epsilon_p, self.rho_solid, self.t_f, self.u_frag,
                                                self._calculate_rho_solid, self._calculate_u_frag)

    def calc_stokes(self):
        """ Compute the stokes number by interpolating the stokes number of the disk to the planets location. """
        self.t_f = np.interp(self.planet.a_p, self.planet.disk.r, self.planet.disk.stokes_number_pebbles)
        self.tau_f = self.t_f / self.planet.omega_k

    def calc_densities(self):
        """ calculate the pebble scalehight and the volume density of pebbles. """
        self.H_p = np.sqrt(self.planet.disk.alpha_height / self.t_f) * self.planet.H_g
        self.rho_peb_components = self.planet.sigma_peb_components / (np.sqrt(2 * np.pi) * self.H_p)
        return

    def accretion_domain(self):  # checked
        """ decide on the regimes of accretion: bondi vs hill """
        if self.regime == "Bondi":
            transition_mass = (
                    np.sqrt(1 / 3)
                    * self.planet.del_v ** 3
                    / (G * self.planet.omega_k)
            )
            if self.planet.M > transition_mass:
                self.regime = "Hill"

    def twoD_transition(self):  # checked
        """ Transition between 2D and 3D pebble accretion. """
        transition = np.pi * self.R_acc / (2 * np.sqrt(2 * np.pi))
        if transition > self.H_p:
            self.twoD = True
        else:
            self.twoD = False

    def calc_R_acc(self):  # checked
        """Calculate the capture radius of pebble accretion. """
        if self.regime == "Bondi":
            t_b = self.planet.r_b / self.planet.del_v  # checked
            self.R_acc = np.sqrt(4 * self.tau_f / t_b) * self.planet.r_b  # checked
        elif self.regime == "Hill":
            self.R_acc = (self.t_f / 0.1) ** (1 / 3) * self.planet.r_h  # checked
        else:
            print("Error, no correct regime specified: ", self.regime)

        if self.regime in ["Hill", "Bondi"]:
            t_p = (
                    G
                    * self.planet.M
                    / ((np.abs(self.planet.del_v) + self.planet.omega_k * self.planet.r_h) ** 3)
            )  # checked
            self.R_acc = self.R_acc * np.exp(-0.4 * (self.tau_f / t_p) ** 0.65)  # checked

        self.del_v = self.planet.del_v + self.planet.omega_k * self.R_acc  # checked
        return

    def remove_mass_from_flux(self):
        """ Ask the disk to remove the accreted pebbles."""
        self.planet.disk.remove_from_peb(self.m_dot_chem, self.planet.idx)
        return

    def calc_m_dot(self):
        """Calculate the pebble accretion rate."""
        if not self.planet.past_pebble_iso:
            self.calc_stokes()
            self.calc_densities()
            self.accretion_domain()
            self.calc_R_acc()
            self.twoD_transition()

            if self.twoD:
                self.m_dot_chem = 2 * self.R_acc * self.planet.sigma_peb_components * self.del_v  # checked
            else:
                self.m_dot_chem = np.pi * self.R_acc ** 2 * self.rho_peb_components * self.del_v  # checked

            self.m_dot_chem = self.factor_peb_accretion * self.m_dot_chem

            self.m_dot_c_chem = (1 - self.planet.solid_frac) * self.m_dot_chem  # checked
            self.m_dot_a_chem = self.planet.solid_frac * self.m_dot_chem  # checked

            self.m_dot = np.sum(self.m_dot_chem[0])  # sum over all elements

            self.m_dot_c = (1 - self.planet.solid_frac) * self.m_dot  # checked
            self.m_dot_a = self.planet.solid_frac * self.m_dot  # checked


        else:
            self.m_dot = 0
            self.m_dot_c = 1 * self.m_dot
            self.m_dot_a = 1 * self.m_dot
            self.m_dot_chem = 0 * self.planet.chemistry_solid
            self.m_dot_c_chem = 0 * self.planet.chemistry_solid
            self.m_dot_a_chem = 0 * self.planet.chemistry_solid

        return
