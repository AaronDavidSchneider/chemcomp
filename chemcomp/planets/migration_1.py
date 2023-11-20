import numpy as np
from astropy import units as u

from . import Planet


class LindbladPlanet(Planet):
    """simple migrating planet. torques calculated with viscosity of α = 0.004"""

    def __init__(
            self,
            output,
            disk,
            accretion_models=[],
        rho_c=5.5 * u.g / (u.cm ** 3.0),
        a_p=5.2 * u.au,
        t_0=0.0 * u.yr,
        M_c=None,
        M_a=None,
        **kwargs
    ):

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

    def update_a_p(self):
        self.update_torque()
        J_p = self.M * self.a_p ** 2.0 * self.omega_k
        self.tau_m = 0.5 * J_p / self._gamma_tot
        self.a_p += self.a_p / self.tau_m * self.h

    def update_torque(self):
        adiabatic_index = 1.4
        q = self.M / self.disk.M_s
        h = self.H_g / self.a_p
        self._gamma_0 = (q / h) ** 2 * self.sigma_g * self.a_p ** 4 * self.omega_k ** 2
        tgrad = np.interp(self.a_p, self.disk.r, self.disk.grad_T)  # to be defined
        Siggrad = np.interp(self.a_p, self.disk.r, self.disk.grad_sig)  # checked
        self._gamma_tot = (
            -(2.5 + 1.7 * Siggrad + 0.1 * tgrad) * self._gamma_0 / adiabatic_index
        )
