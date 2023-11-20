import astropy.constants as const
import astropy.units as u
import numpy as np
from scipy.interpolate import interp1d

from . import Disk
from chemcomp.helpers import eval_kwargs
from ._0pacity import radmc3d

opacity_fct = radmc3d

G = const.G.cgs.value

class TwoPopDisk(Disk):
    """ LBP, Disk of the twopoppy model."""

    def __init__(
            self,
            defaults,
            chemistry,
            M_STAR=0.7 * u.M_sun,
            M0=1e-1 * u.M_sun,
            R0=20 * u.au,
            time=1 * u.yr,
            ALPHA=1e-3,
            output={},
            **kwargs
    ):
        super().__init__(defaults, M_STAR, ALPHA, chemistry, time=time, **kwargs)

        self.tstar = eval_kwargs(kwargs.get("tstar", 4000.0 * u.K))
        self.rstar = eval_kwargs(kwargs.get("rstar", 1.8 * u.R_sun))
        self.rstar = self.rstar.cgs.value
        self.tstar = self.tstar.cgs.value

        self.compute_disktmid()
        self.make_disk_from_lbp_alpha(M0=M0.cgs.value, R0=R0.cgs.value, alpha=ALPHA, time=time.cgs.value)

        self.sigma_g = self.sigma_g / np.trapz(2 * np.pi * self.r * self.sigma_g, x=self.r) * M0.cgs.value

        self.compute_disktmid()
        self.compute_disktmid()
        self._init_params(**kwargs)
        return

    def evolve_init(self):
        pass
        # DEBUGGING
        #mask = np.where(self.r < (5 * u.au).cgs.value)[0]
        #self.sigma_g_components[mask, 1, 7] = 1e-60
        #self.sigma_dust_components[mask, :, :] = 1e-60
        #self.compute_sigma_g_from_mol(self.sigma_g_components[:, 1, :])
        #self.compute_sigma_dust_from_mol(self.sigma_dust_components[:, 1, :])
        #self.sigma_dust_components = self.chemistry_solid * self.sigma_dust[:, np.newaxis, np.newaxis]
        #self.sigma_g_components = self.chemistry_gas * self.sigma_g[:, np.newaxis, np.newaxis]

    def compute_disktmid(self, viscmodel=None):
        self.T = ((0.05 ** 0.25 * self.tstar * (self.r / self.rstar) ** -0.5) ** 4 + (7.) ** 4.0) ** 0.25
        self.tvisc = np.zeros_like(self.T)
        self.tirr = np.zeros_like(self.T)
        self.compute_cs_and_hp()
        if hasattr(self, "rho_g"):
            self.calc_opacity()

    def make_disk_from_lyndenbellpringle(self, M0, R0, nu0, gam, time):
        """
            Make a model of the gas surface density of the disk.
            Here the model is set by the Lynden-Bell & Pringle solution.
            I follow here Lodato, Scardoni, Manara & Testi (2017)
            MNRAS 472, 4700, their Eqs. 2, 3 and 4.

            ARGUMENTS:
              M0      = Initial disk mass [g]
              R0      = Initial disk radius [cm]
              nu0     = Viscosity nu at R0 [cm^2/s]
              gam     = The gamma coefficient of the LBP model
              time    = Time after start [s]
            """
        r = self.r
        nu = nu0 * (r / R0) ** gam  # Eq. 2
        tnu = R0 ** 2.0 / (3.0 * (2. - gam) ** 2.0 * nu0)  # Eq. 4
        T = 1.0 + time / tnu  # Above Eq. 4
        eta = (2.5 - gam) / (2. - gam)  # Below Eq. 3
        self.sigma_g = (M0 / (2.0 * np.pi * R0 ** 2)) * (2. - gam) * \
                       (R0 / r) ** gam * T ** (-eta) * \
                       np.exp(-(r / R0) ** (2. - gam) / T)  # Eq. 3
        self.nu = nu
        self.gam = gam
        self.sigma_g[self.sigma_g < self.sigmin] = self.sigmin


    def make_disk_from_lbp_alpha(self, M0, R0, alpha, time):
        """
        Same as make_disk_from_lyndenbellpringle(), but now assuming a constant
        alpha viscosity model. That is: both nu0 and gamma are internally
        computed.

        ARGUMENTS:
          M0      = Initial disk mass [g]
          R0      = Initial disk radius [cm]
          alpha   = The Shakura and Sunyaev alpha value
          time    = Time after start [s]

        NOTE: self.aspect_ratio must be a powerlaw function!
        """
        assert hasattr(self, "aspect_ratio")

        self.omega_k = np.sqrt(G * self.M_s / self.r ** 3.0)
        self.alpha = alpha
        r = self.r
        aspect_ratio = 1.0 * self.aspect_ratio
        omk = self.omega_k
        cs = aspect_ratio * omk * r
        gam = 1.0
        om1 = self.omega_k[0]
        cs1 = cs[0]
        nu0 = self.alpha * cs1 ** 2.0 / om1

        self.make_disk_from_lyndenbellpringle(M0, R0, nu0, gam, time)

    def update_parameters(self):
        self.compute_viscous_evolution()

