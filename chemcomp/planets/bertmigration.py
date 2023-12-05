import numpy as np
from astropy import constants as const
from astropy import units as u

from chemcomp.helpers import eval_kwargs
from . import Planet

sr = const.sigma_sb.cgs.value
G = const.G.cgs.value


class BertPlanet(Planet):
    """Migrating planet. torques adapted from Planettrack.F by Bertram Bitsch. Planet used in Schneider & Bitsch (2021a,b)"""

    def __init__(
        self,
        output,
        disk,
        accretion_models=[],
        rho_c=5.5 * u.g / (u.cm**3),
        a_p=5.2 * u.au,
        t_0=5e5 * u.yr,
        M_c=None,
        M_a=None,
        **kwargs
    ):
        self.gamma_tot = 0.0  # total torque
        self.gamma_norm = 0.0  # torque normalized to gamma_0

        self._use_heat_torque = eval_kwargs(
            kwargs.get("use_heat_torque", True)
        )  # use masset heating torque
        self._use_dynamical_torque = eval_kwargs(
            kwargs.get("use_dynamical_torque", True)
        )  # use paardekooper dynamical torque
        self._apply_gap = eval_kwargs(
            kwargs.get("apply_gap_profile", True)
        )  # if the gap should be added to the gas surface density
        self._migrate = eval_kwargs(
            kwargs.get("migration", True)
        )  # if the planet should migrate

        super().__init__(
            a_p=a_p,
            M_c=M_c,
            M_a=M_a,
            t_0=t_0,
            disk=disk,
            rho_c=rho_c,
            accretion_models=accretion_models,
            output=output,
            **kwargs
        )

        self._old_a_p = 1.0 * self.a_p
        self.init_paardekooper_units()
        print(self.name, (self.M * u.g).to(u.earthMass).value)

    def update_all(self):
        """Update all quantities including disk evolution, migration and planetary quantities.
        Overwrites parent update_all() method"""
        if self._apply_gap:
            self.update(second=False)
            self.calc_gammaeff(migphase=False)
            gap_width = 2.0 * self._xs * self.a_p
            # evaluates true outside the two sigma zone of the gap profile.
            if self.fSigma < 1 and self.past_pebble_iso:
                self.disk.apply_gap_profile(self.a_p, self.fSigma, gap_width)

        interpolation_idx = (
            self.disk.interpolation_idx
            if hasattr(self.disk, "interpolation_idx") and self._apply_gap
            else None
        )
        interpolation_idx_interface = (
            self.disk.interpolation_idx_interface
            if hasattr(self.disk, "interpolation_idx_interface") and self._apply_gap
            else None
        )
        self.disk.update(self.t, self.h)
        self.update(
            second=False,
            interpolation_idx=interpolation_idx,
            interpolation_idx_interface=interpolation_idx_interface,
        )
        self.update_a_p()
        self.update(
            second=True,
            interpolation_idx=interpolation_idx,
            interpolation_idx_interface=interpolation_idx_interface,
        )  # update all quantities effected by a_p

    def init_paardekooper_units(self):
        """Define the paardekooper scaling factors."""
        self._r_0 = 1.0 * self.a_p
        self._omega_0 = 1.0 * self.omega_k
        self._sigma_0 = 1.0 * self.sigma_g
        self._q_d = np.pi * self._r_0**2.0 * self._sigma_0 / self.disk.M_s
        self._nu_0 = 1.0 * self.nu

        # altough mig beta and alpha will change during runtime, we keep them constant since a gap would manipulate this even more!
        self.mig_beta = -self.grad_T  # checked
        self.mig_alpha = -self.grad_sig  # checked
        self.xi = (
            self.mig_beta - (self.disk._chem.gamma - 1.0) * self.mig_alpha
        )  # checked, see equation 5

    def calc_heating_torque(self):
        """
        Calculate the heating torques.
        might be applied afer update_torque to add heating torque
        """

        thermal_crit = self.Vchi * self.c_s / G

        if self.M < thermal_crit:
            L_c_mod_M = 4.0 * np.pi * G * self.Vchi * self.rho_g / self._gammaeff
            L_mod_M = (
                G
                * (self.accretion_models[0].m_dot + self.accretion_models[1].m_dot)
                / self.r_c
            )
            L_norm = L_mod_M / L_c_mod_M

            lambda_c = np.sqrt(self.Vchi / (1.5 * self.omega_k * self._gammaeff))
            eta = self.mig_alpha / 3.0 + (self.mig_beta + 3.0) / 6.0

            gamma_thermal = (
                1.61
                * (self._gammaeff - 1.0)
                / self._gammaeff
                * eta
                * (self.H_g / lambda_c)
                * (L_norm - 1.0)
            )

            self.gamma_tot += self.gamma_zero * gamma_thermal
            self.gamma_norm += gamma_thermal

        return

    def update_a_p(self):
        """
        Migration.
        1. Calculate the torques
        2. Calculate the migration rate
        3. Adapt migration rate to typeII migration if needed
        """
        self.mig_beta = -self.grad_T
        self.mig_alpha = -self.grad_sig
        self.xi = self.mig_beta - (self.disk._chem.gamma - 1.0) * self.mig_alpha

        self.update_torque()
        if self._use_heat_torque:
            self.calc_heating_torque()

        if self._use_dynamical_torque:
            self.tau_m = self.get_dynamical_tau_m()
        else:
            self.tau_m = self.get_static_tau_m()

        if self.fSigma <= 0.53:
            tau_II = (
                -self.a_p**2.0
                / self.nu
                * np.maximum(1.0, self.M / (4 * np.pi * self.sigma_g * self.a_p**2.0))
            )
            if self.fSigma > 0.1:
                self.tau_m = (self.tau_m - tau_II) / (0.53 - 0.1) * (
                    self.fSigma - 0.1
                ) + tau_II
            else:
                self.tau_m = tau_II

        if (self.fSigma >= 0.8 or self.tau_m <= 0.0) and self._migrate:
            # update migration from above if not type 2 and directed inwards
            self.a_p += self.a_p / self.tau_m * self.h

        self.a_p = np.clip(self.a_p, self.r_in, self.disk.r[-1])

    def get_dynamical_tau_m(self):
        """
        Paardekooper 2014 dynamical torques

        Returns
        -------
        tau_m: Migration timescale [s]

        """

        k = (
            8.0
            / 3.0
            / np.pi
            * (3.0 / 2.0 - self.mig_alpha)
            * self.gamma_norm
            * self._q_d**2.0
            * self._xs**3.0
            / (self.aspect_ratio**2.0)
            * self._r_0**2.0
            * self._omega_0
            / self._nu_0
            * (self.a_p / self._r_0) ** (5.0 - 3.0 * self.mig_alpha)
        )
        k = np.minimum(k, 0.5)
        theta = (1.0 - np.sqrt(1.0 - 2.0 * k)) / k
        J_p = self.M * self.a_p**2.0 * self.omega_k
        return 0.5 * J_p / self.gamma_tot / theta

    def get_static_tau_m(self):
        """
        Only static torques

        Returns
        -------
        tau_m: Migration timescale [s]

        """
        J_p = self.M * self.a_p**2.0 * self.omega_k

        return 0.5 * J_p / self.gamma_tot

    @property
    def _kappa(self):
        """Interpolate opacity to planetary position. Adds minimal offset (for stability)"""
        return np.interp(self.a_p, self.disk.r, self.disk.opacity) + 1.0e-300

    def update_torque(self):
        """Calculate Paardekooper2011 torques"""
        self.calc_gammaeff(migphase=True)
        ang_p = self.a_p**2.0 * self.omega_k  # spec. Ang.mom of planet, checked
        pnu = (
            2.0 / 3.0 * np.sqrt((ang_p * self._xs**3.0) / (2.0 * np.pi * self.nu))
        )  # checked
        pxi = np.sqrt((ang_p * self._xs**3.0) / (2.0 * np.pi * self.Vchi))  # checked
        Fpnu = 1.0 / (1.0 + (pnu / 1.3) ** 2.0)  # checked
        Fpxi = 1.0 / (1.0 + (pxi / 1.3) ** 2.0)  # checked

        if pnu < np.sqrt(8.0 / (45.0 * np.pi)):
            Gpnu = (
                16.0 / 25.0 * (45.0 * np.pi / 8.0) ** (0.75) * pnu ** (1.5)
            )  # checked
        else:
            Gpnu = 1.0 - 9.0 / 25.0 * (8.0 / (45.0 * np.pi)) ** (1.33333) * pnu ** (
                -8.0 / 3.0
            )  # checked
        if (pxi) < (np.sqrt(8.0 / (45.0 * np.pi))):
            Gpxi = (
                16.0 / 25.0 * (45.0 * np.pi / 8.0) ** (0.75) * pxi ** (1.5)
            )  # checked
        else:
            Gpxi = 1.0 - 9.0 / 25.0 * (8.0 / (45.0 * np.pi)) ** (1.33333) * pxi ** (
                -8.0 / 3.0
            )  # checked

        if (pnu) < (np.sqrt(28.0 / (45.0 * np.pi))):
            Kpnu = 16.0 / 25.0 * (45.0 * np.pi / 28.0) ** 0.75 * pnu**1.5  # checked
        else:
            Kpnu = 1.0 - 9.0 / 25.0 * (28.0 / (45.0 * np.pi)) ** (1.33333) * pnu ** (
                -8.0 / 3.0
            )  # checked

        if (pxi) < (np.sqrt(28.0 / (45.0 * np.pi))):
            Kpxi = (
                16.0 / 25.0 * (45.0 * np.pi / 28.0) ** (0.75) * pxi ** (1.5)
            )  # checked
        else:
            Kpxi = 1.0 - 9.0 / 25.0 * (28.0 / (45.0 * np.pi)) ** (1.33333) * pxi ** (
                -8.0 / 3.0
            )  # checked

        # Calculate torque contributions
        self.gamma_zero = (
            (self.massratio / self.aspect_ratio) ** 2.0
            * self.sigma_g
            * self.a_p**4.0
            * self.omega_k**2.0
        )  # checked

        if self.fSigma < 1.0 and self.fSigma > 0.1:
            self.gamma_zero = self.gamma_zero * self.fSigma

        self._m_gap = (
            np.sqrt(self.aspect_ratio**5.0 * self.disk.alpha_mig * 1.0 / 0.04)
            * self.disk.M_s
        )
        self._iso_vs_gap = self.peb_iso / self._m_gap

        xfact = self.gamma_zero / self._gammaeff  # checked
        GammaLIN = (
            -(2.5 + 1.7 * self.mig_beta - 0.1 * self.mig_alpha) * xfact
        )  # checked

        GammaHSENT = self.xi / self._gammaeff * 7.9 * xfact  # equation 5, checked
        GammaCLINENT = (
            (2.2 - 1.4 / self._gammaeff) * self.xi * xfact
        )  # equation 7, checked

        GammaCLINBARO = 0.7 * (1.5 - self.mig_alpha) * xfact  # equation 6, checked
        GammaHSBARO = 1.1 * (1.5 - self.mig_alpha) * xfact  # equation 4, checked

        GammaCBARO = GammaHSBARO * Fpnu * Gpnu + (1.0 - Kpnu) * GammaCLINBARO  # checked
        GammaCENT = (
            GammaHSENT * Fpnu * Fpxi * np.sqrt(Gpnu * Gpxi)
            + np.sqrt((1.0 - Kpnu) * (1.0 - Kpxi)) * GammaCLINENT
        )  # checked

        # Normalization factor to Gamma_0 as defined in Paardekooper 2011
        self.gamma_tot = GammaCBARO + GammaCENT + GammaLIN  # checked
        self.gamma_norm = self.gamma_tot / self.gamma_zero  # checked

    def calc_gammaeff(self, migphase):
        """Calculate the effective gamma AND calculate the gapdepth"""
        self.massratio = self.M / self.disk.M_s  # checked
        # 	rename gradients to self.mig_alpha and beta as in Paardekooper 2011
        # 	Calculate Vchi
        self.Vchi = (  # eq. 4 in Paardekooper 2011
            16.0  # faktor of 4 above paardekooper
            * self.disk._chem.gamma
            * (self.disk._chem.gamma - 1.0)
            * sr
            * self.T**4.0
            / (
                3.0
                * self._kappa
                * self.rho_g**2.0
                * self.H_g**2
                * self.omega_k**2
            )
        )  # checked
        # calculate Q
        Q = (
            2.0
            * self.Vchi
            / (3.0 * (self.H_g / self.a_p) ** 3.0)
            / (self.a_p**2.0 * self.omega_k)
        )  # checked
        self._gammaeff = (
            2.0
            * Q
            * self.disk._chem.gamma
            / (
                self.disk._chem.gamma * Q
                + 0.5
                * np.sqrt(
                    2.0
                    * np.sqrt(
                        (self.disk._chem.gamma**2.0 * Q**2.0 + 1.0) ** 2.0
                        - 16.0 * Q**2.0 * (self.disk._chem.gamma - 1.0)
                    )
                    + 2.0 * self.disk._chem.gamma**2.0 * Q**2.0
                    - 2.0
                )
            )
        )  # checked
        #  Calculate functions F,G,K
        self.aspect_ratio = self.H_g / self.a_p  # relative disk thickness
        self._xs = (
            1.11 / self._gammaeff**0.25 * np.sqrt(self.massratio / self.aspect_ratio)
        )  # checked, dimensionless in units of r_p

        #######
        # GAP #
        #######
        if not hasattr(self, "_M_HS"):
            self._M_HS = 4 * np.pi * self._xs * self.a_p**2 * self.sigma_g
            self._sigma_HS = 1 * self.sigma_g

        K = self.massratio**2.0 / (
            self.aspect_ratio**5.0 * self.disk.alpha_mig
        )  # checked
        self.fSigma_kanagawa = 1.0 / (1.0 + 0.04 * K)  # checked
        R = self.a_p**2.0 * self.omega_k / self.nu
        Gap = 3.0 / 4.0 * self.H_g / self.r_h + 50.0 / (self.massratio * R)

        if Gap <= 2.4646:
            f_P = (Gap - 0.541) / 4.0
        else:
            f_P = 1.0 - np.exp(-(Gap**0.75) / 3.0)

        if self.past_pebble_iso:
            # force m_dot_gas = m_dot_a as specified in the gas.py file
            m_dot_gas = self._accretion_models_dict["GasAccretion"].m_dot_a

            if migphase:
                self._sigma_HS = self._sigma_HS * (self._old_a_p / self.a_p) ** 1.5
                self._old_a_p = 1.0 * self.a_p
                M_HS_check = 4.0 * np.pi * self.a_p**2.0 * self._xs * self._sigma_HS

                if M_HS_check < self._M_HS:
                    self._M_HS = self._M_HS + (self.m_dot_disk - m_dot_gas) * self.h
                else:
                    self._M_HS = (
                        self._M_HS
                        + (
                            4.0 * np.pi * self.a_p**2.0 * self._xs
                            - self._M_HS / self._sigma_HS
                        )
                        * self.sigma_g
                    )

                self._sigma_HS = self._M_HS / (4.0 * np.pi * self.a_p**2.0 * self._xs)

                if self.t < self.disk.begin_photevap:
                    assert self._M_HS >= 0.0, "negative mass in horseshoe region"
                else:
                    self._M_HS = np.maximum(self._M_HS, 1e-60)

            f_A = 1.0 - m_dot_gas * self.h / (f_P * self._M_HS)
        else:
            f_A = 1.0

        self.fSigma = np.maximum(f_A * f_P, 1e-7)
        return
