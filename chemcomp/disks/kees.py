import astropy.constants as const
import astropy.units as u
import numpy as np
from scipy.interpolate import interp1d

from . import Disk
from chemcomp.helpers import eval_kwargs

AU = const.au.cgs.value
G = const.G.cgs.value


class DiskLabDisk(Disk):
    """ LBP, viscous disk evolution and viscous heating. Disk of Schneider & Bitsch (2021)a,b."""

    def __init__(
            self,
            defaults,
            chemistry,
            M_STAR=1.0 * u.M_sun,
            M0=1e-2 * u.M_sun,
            R0=200 * u.au,
            time=1.0 * u.yr,
            ALPHA=1e-3,
            **kwargs
    ):
        super().__init__(defaults, M_STAR, ALPHA, chemistry, time=time, **kwargs)

        self.flang = eval_kwargs(kwargs.get("flang", 0.05))
        self.lstar = eval_kwargs(kwargs.get("lstar", 1.0 * const.L_sun).cgs.value)

        self._M0 = M0.cgs.value
        self._R0 = R0.cgs.value

        self.init_T()
        t_iter = 0
        while True:
            self.make_disk_from_lbp_alpha(M0=M0.cgs.value, R0=R0.cgs.value, alpha=float(ALPHA), time=time.cgs.value)
            T0 = self.T.copy()
            self.compute_disktmid()
            self._chem.get_composition(self.T)  # updates mu
            self._init_params(**kwargs)
            if np.all(np.abs((T0-self.T)/T0) < 1e-2):
                break

            assert t_iter < 1000, "We have trouble to converge the initial T-profile"
            t_iter += 1

        self._T_init_finished = True
        return

    def make_disk_from_simplified_lbp(self, Sigc, Rc, gam):
        """
        Make a model of the gas surface density of the disk.
        Here the model is set by the simplified LBP-like powerlaw
        disk with an exponential cut-off with critical radius Rc.
        I follow here Takahashi & Inutsuka (2016)
        AJ 152, 184, their Eq. 18.

        ARGUMENTS:
          Sigc    = Normalization constant [g/cm^2]
          Rc      = Critical disk radius [cm]
          gam     = The gamma coefficient of the LBP model
        """
        r = self.r / AU
        sigma = Sigc * (Rc / r) ** gam * np.exp(-(r / Rc) ** (2. - gam))
        sigma[sigma < self.sigmin] = self.sigmin

        self.sigma_g = sigma

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
        tnu = R0 ** 2 / (3.0 * (2. - gam) ** 2 * nu0)  # Eq. 4
        T = 1.0 + time / tnu  # Above Eq. 4
        eta = (2.5 - gam) / (2. - gam)  # Below Eq. 3
        self.sigma_g = (M0 / (2 * np.pi * R0 ** 2.0)) * (2. - gam) * \
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
        aspect_ratio = 1 * self.aspect_ratio
        omk = self.omega_k
        cs = aspect_ratio * omk * r
        nu = alpha * cs * cs / omk
        gam = (np.log(nu[-1]) - np.log(nu[0])) / (np.log(r[-1]) - np.log(r[0]))
        f = interp1d(r, nu)
        nu0 = f([R0])[0]
        self.make_disk_from_lyndenbellpringle(M0, R0, nu0, gam, time)

    def update_parameters(self):
        """ Things that should be called whenever the update_parameters() function is called in the planet class. """
        self.compute_viscous_evolution()
