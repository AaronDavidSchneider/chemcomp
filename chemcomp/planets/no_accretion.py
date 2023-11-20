from astropy import units as u

from . import Planet


class NoAccretion(Planet):
    """NoAccretion planet. Useful if you are just interested in the disk and not in the planet itself"""

    def __init__(
            self,
            output,
            disk,
            accretion_models=[],
            rho_c=5.5 * u.g / (u.cm ** 3),
            a_p=5.2 * u.au,
            t_0=0.0 * u.yr,
            M_c=None,
            M_a=None,
            **kwargs,
    ):
        self.gamma_tot = 0
        self.gamma_norm = 0
        super().__init__(
            a_p=a_p,
            M_c=M_c,
            M_a=M_a,
            t_0=t_0,
            rho_c=rho_c,
            disk=disk,
            accretion_models=accretion_models,
            output=output,
            **kwargs,
        )

    def grow_mass(self):
        self.output.finish()
