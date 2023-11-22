import numpy as np
import scipy.optimize as optimize
from astropy import constants as const
from astropy import units as u
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d

from chemcomp.helpers import solvediffonedee_components, getfluxonedee, eval_kwargs
from . import opacity_fct

sr = const.sigma_sb.cgs.value
year = 3600 * 24 * 365.25
pi = np.pi
k_B = const.k_B.cgs.value
m_p = const.m_p.cgs.value
G = const.G.cgs.value
AU = const.au.cgs.value
ME = const.M_earth.cgs.value

FWHM_TO_SIGMA = 2 * (2 * np.log(2)) ** 0.5


def gaussian(x, max, mu, FWHM):
    """Gaussian helper function."""
    # prevent underflow by maximizing exponent to -1e2 (np.exp(-1e2)=3e-44)
    exp = np.minimum(0.5 * ((x - mu) / FWHM * FWHM_TO_SIGMA) ** 2.0, 1e2)
    return max * np.exp(-exp)


def tvisc_fct(T, sigma_g, tirr, omega_k, alpha, Z, mu):
    """Calculate viscous heating and return the difference between our current guess of the temperature
    to the last guess of the temperature"""
    c_s = (k_B * T / (mu * m_p)) ** 0.5
    h = c_s / omega_k
    rho_g = sigma_g / (np.sqrt(2 * pi) * h)

    nu = alpha * c_s * c_s / omega_k

    opacity = opacity_fct(rho_g, T, Z=Z) * 100.0 * Z

    qvisc = (9.0 / 4.0) * sigma_g * nu * omega_k**2.0
    tau = sigma_g * opacity  # both sides of the disk
    tvisc4 = (
        qvisc * 3.0 * tau / (16.0 * sr)
    )  # factor 0.5 because tau is defined on complete disk

    T_new = np.clip((tvisc4 + tirr**4) ** 0.25, tirr, 1.5e5)

    diff = T - T_new

    return diff


def solve_viscous_heating_globally(
    sigma_g,
    tirr,
    T,
    omega_k,
    alpha,
    Z,
    mu,
    r,
    xtol=1e-1,
    nitermax=100,
    interpolation_idx=None,
):
    """
    Solve the midplane tempeature due to viscous heating (only)
    Note: adapt _0pacity model if needed! only works for low alpha or very smooth _0pacity
    """
    if interpolation_idx is None:
        interpolation_idx = np.ones_like(T, dtype=bool)

    sample = np.linspace(0.1, 300.0, 10000) * AU

    T = np.array(
        [
            optimize.root_scalar(
                tvisc_fct,
                bracket=(tirr[i], 1.5e4),
                x0=T[i],  # estimate T from last iteration
                args=(
                    sigma_g[i],
                    tirr[i],
                    omega_k[i],
                    alpha[i] if hasattr(alpha, "__iter__") else alpha,
                    Z[i] if hasattr(Z, "__iter__") else Z,
                    mu[i] if hasattr(mu, "__iter__") else mu,
                ),
                options={"maxiter": nitermax, "xtol": xtol},
            ).root
            for i in range(len(T))
            if interpolation_idx[i]
        ]
    )
    T_sample = savgol_filter(
        interp1d(r[interpolation_idx], T, fill_value="extrapolate", assume_sorted=True)(
            sample
        ),
        5,
        2,
    )
    T = interp1d(sample, T_sample, fill_value="extrapolate", assume_sorted=True)(r)
    T = np.clip(T, tirr, 1.5e5)
    return T


class Disk(object):
    """Baseclass of the disk."""

    def __init__(
        self,
        defaults,
        M_STAR,
        ALPHA,
        chemistry,
        init_grid=True,
        time=0 * u.yr,
        **kwargs,
    ):
        self.sigmin = 1e-60  # sigma_floor value
        self.defaults = (
            defaults  # dictionary containing the defaults section of the config
        )
        self.M_s = M_STAR.cgs.value  # stellar mass
        self.alpha = float(ALPHA)  # alpha viscosity, fix yaml issue
        self._chem = chemistry  # reference to chemistry object
        self.t_0 = 1 * time.cgs.value  # starting time of the disk
        self.t = 1 * time.cgs.value  # time in the disk (initialised with starting time)

        # init timestep
        self._timestep_input = defaults.get("DEF_TIMESTEP", None)
        if self._timestep_input is not None:
            print(f"WARNING: We are using a custom timestep of {self._timestep_input}")
            self._timestep_input = self._timestep_input.cgs.value

        # init t_end
        self.t_end = defaults.get("DEF_T_END", 100.0 * u.Myr).cgs.value

        # init grid
        if init_grid:
            r_in = defaults.get("DEF_R_IN", 0.2 * u.au)  # inner grid edge
            r_out = defaults.get("DEF_R_OUT", 100.0 * u.au)  # outer grid edge

            self.gridsize = defaults.get("DEF_GRIDSIZE", 1000)  # number of grid cells

            if defaults.get("DEF_LIN_SPACING", True):
                self.r = (
                    np.linspace(r_in, r_out, self.gridsize).to(u.au).value
                )  # linear grid
            else:
                self.r = np.logspace(  # log grid
                    np.log10(r_in.to(u.au).value),
                    np.log10(r_out.to(u.au).value),
                    self.gridsize,
                )

            self.r_i = (
                np.array([self.r[0], *((self.r[1:] + self.r[:-1]) / 2.0), self.r[-1]])
                * u.au
            )
            self.r = self.r * u.au

            self.r_i = self.r_i.cgs.value  # interfacepositions of the radial grid
            self.r = self.r.cgs.value  # cellcenter of the radial grid

        self._init_params(**kwargs)

    def evolve_init(self):
        """function that can be overwritten to be called at the start of the disk evolution"""
        pass

    def _init_params(self, **kwargs):
        """
        Initialize the rest of the parameters.
        Most parameters can only be set if sigma and aspect_ratio have been set before

        Parameters
        ----------
        kwargs: configuration from config_disk section of config.yaml
        """
        # array holding the dust to gas ratio
        self.DTG = {"total": eval_kwargs(kwargs.pop("DTG_total", 0.015))}

        self.DTG_small_grains = self.DTG["total"] * (
            1.0 - 0.75
        )  # DTG small grains is only used for the temperature
        self.conf_static = eval_kwargs(
            kwargs.get("static", False)
        )  # config for static disk
        self.conf_static_stokes = eval_kwargs(
            kwargs.get("static_stokes", False)
        )  # config for static stokes number
        self.conf_evaporation = eval_kwargs(
            kwargs.get("evaporation", True)
        )  # config for static stokes number
        self.conf_temp_evol = eval_kwargs(
            kwargs.get("temp_evol", False)
        )  # use the temperature evolution. Advise: dont do!

        self._evap_width = eval_kwargs(
            kwargs.get("evap_width", 0.1 * u.au)
        ).cgs.value  # config for static stokes number

        self.alpha_height = eval_kwargs(
            kwargs.get("ALPHAHEIGHT", self.alpha)
        )  # vertical mixing alpha
        self.alpha_frag = eval_kwargs(
            kwargs.get("ALPHAFRAG", self.alpha)
        )  # fragmentation alpha, not used in Schneider & Bitsch 2021
        self.alpha_mig = eval_kwargs(
            kwargs.get("ALPHAMIG", self.alpha)
        )  # migration alpha, not used in Schneider & Bitsch 2021

        self._chem.set_z(self.DTG.get("total", 0.015))  # set the solid to gas ratio

        if self.conf_static_stokes:
            self.f_m = np.ones_like(self.r)  # Massdistribution factor for twopop model
            self.stokes_number_small = np.zeros_like(
                self.r
            )  # stokes number of the small population for twopop model
            self.stokes_number_df = np.zeros_like(
                self.r
            )  # stokes number of the drift limited fragmentation for twopop model
            self.stokes_number_drift = np.zeros_like(
                self.r
            )  # stokes number of the drift limited case for twopop model
            self.stokes_number_frag = np.zeros_like(
                self.r
            )  # stokes number of the fragmentation limited case for twopop model

        self.const_f_f = 0.37  # fudge fitting factors from Birnstiel2012
        self.const_f_d = 0.55  # fudge fitting factors from Birnstiel2012
        self.const_N = 0.5  # N of Birnstiel2012
        self.a_0 = eval_kwargs(
            kwargs.get("a_0", 1e-4 * u.cm)
        )  # grainsize of the small population
        self.a_0 = self.a_0.cgs.value  # grainsize of the small population
        self.a_1 = eval_kwargs(kwargs.get("a_0", 1e-4 * u.cm)) * np.ones_like(
            self.r
        )  # grainsize of the large population
        self.a_1 = self.a_1.cgs.value  # grainsize of the large population

        self.omega_k = np.sqrt(G * self.M_s / self.r**3)  # kepler orbitalperiod
        self.v_k = self.omega_k * self.r  # kepler velocity
        self.alpha_factor = np.ones_like(
            self.r
        )  # gap factor that applies the gap to the alphaviscosity

        # Init for output
        self.vr_dust = np.zeros_like(self.r[1:])  # dust velocity
        self.vr_gas = np.zeros_like(self.r[1:])  # gas velocity

        tau_disk = eval_kwargs(kwargs.get("tau_disk", None))  # disk dispersal time
        self.begin_photevap = eval_kwargs(
            kwargs.get("begin_photoevap", 0 * u.Myr).cgs.value
        )  # time until disk disappears

        if tau_disk is not None:
            self.conf_photoevaporation = True
            self.tau_disk = 1 * tau_disk.cgs.value
        else:
            self.conf_photoevaporation = False

        if hasattr(self, "aspect_ratio") and hasattr(self, "sigma_g"):
            # init physical quantities

            # calculate the chemical composition of sigma:
            self.chemistry_solid, self.chemistry_gas = self._chem.get_composition(
                self.T
            )  # chemistry vectors
            dtg = self._chem.get_solid_heavies(self.T)
            self.sigma_g_components = (
                self.chemistry_gas * self.sigma_g[:, np.newaxis, np.newaxis]
            )  # gas surface density as a compositional vector
            self.sigma_g = np.sum(
                self.sigma_g_components[:, 0, :], axis=1
            )  # gas surface density

            self.rho_g = self.sigma_g / (
                np.sqrt(2 * np.pi) * self.aspect_ratio * self.r
            )  # gas volume density
            self.grad_sig = np.gradient(  # gradient of the surfacedensity
                np.log10(self.sigma_g), np.log10(self.r)
            )  # -1.5 for MMSN
            self.c_s = self.aspect_ratio * self.r * self.omega_k  # soundspeed
            self.P = (  # pressure
                self.c_s**2
                * self.sigma_g
                / (np.sqrt(2 * np.pi) * self.aspect_ratio * self.r)
            )
            self.grad_p = np.gradient(  # pressuregradient
                np.log10(self.P), np.log10(self.r)
            )  # -2.75 for MMSN, checked
            self.T = 2.35 * self.c_s**2.0 * m_p / k_B  # Temperature
            if not hasattr(self, "T_irr"):
                self.T_irr = 1.0 * self.T  # effective Temperature from irradiation
            if not hasattr(self, "T_visc"):
                self.T_visc = np.zeros_like(
                    self.T
                )  # viscous Temperature from viscous heating

            self.compute_cs_and_hp()
            self.compute_nu()
            self.calc_opacity()

            self.grad_T = np.gradient(  # Temperature gradient
                np.log10(self.T), np.log10(self.r)
            )  # -0.5 for MMSN
            self.eta = -0.5 * self.aspect_ratio**2 * self.grad_p

            enclosed_mask = self.r < 200.0 * AU
            self.disk_mass = np.trapz(
                2 * np.pi * self.sigma_g[enclosed_mask] * self.r[enclosed_mask],
                self.r[enclosed_mask],
            )
            self.disk_mass_dust = 0

            self.m_dot_components = self.compute_mdot_at_interfaces(
                self.sigma_g_components
            )
            self.m_dot = np.sum(self.m_dot_components[:, 0, :], axis=1)

            self.pebble_flux = np.zeros_like(self.r)
            self.cum_pebble_flux = np.zeros_like(self.r)

            self.Q = self.c_s * self.omega_k / (np.pi * G * self.sigma_g)
            if np.any(self.Q < 1):
                print("caution: Unstable disk: Q_min = {}".format(np.min(self.Q)))

    def calc_pebble_iso(self, diffusion=False):
        """Calculate the pebble isolation mass."""
        f_fit = (self.aspect_ratio / 0.05) ** 3.0 * (
            0.34 * (np.log10(1e-3) / np.log10(self.alpha)) ** 4.0 + 0.66
        )

        if diffusion:
            self.peb_iso_pi_crit = self.alpha / (2.0 * self.stokes_number_pebbles)
            self.peb_iso = (25.0 + self.peb_iso_pi_crit / 0.00476) * ME * f_fit
        else:
            self.peb_iso = 25.0 * ME * f_fit

    def init_pebble_parameters(
        self, epsilon_p, rho_solid, t_f, u_frag, calculate_rho_solid, calculate_u_frag
    ):
        """
        Sets the pebble parameters needed for pebble evaporation model
        all in cgs

        Parameters
        ----------
        epsilon_p: sticking efficiency
        rho_solid: internal dust density
        t_f: Stokes number
        u_frag: fragmentation velocity
        calculate_rho_solid: if rho solid should be calculated according to Drażkowska17

        """
        self.u_frag = u_frag
        self.epsilon_p = epsilon_p
        self.rho_solid = rho_solid
        if self.rho_solid is not None:
            self.rho_solid = np.ones_like(self.sigma_dust) * self.rho_solid

        self._calculate_rho_solid = calculate_rho_solid
        self._calculate_u_frag = calculate_u_frag

        self.calc_sigma_dust()
        self.compute_disktmid()
        self.update_parameters()

        if t_f is not None:
            # If a stokes number is specified in pebble config, this will be used. TwoPop wont be used then!
            self.stokes_number_pebbles = np.ones_like(self.aspect_ratio) * t_f
            self.conf_static_stokes = True
        else:
            self.conf_static_stokes = False
            self.compute_sizes()

        self.calc_sigma_dust()
        self.compute_disktmid()
        self.update_parameters()

        if not self.conf_static_stokes:
            self.compute_sizes()

        self.calc_sigma_dust()
        self.compute_disktmid()
        self.update_parameters()
        return

    def calc_sigma_dust(self):
        """
        initializes the dust surface density

        This is equal to the pebble surface density, if DTG: peb has been set. Otherwise it contains the whole dust (including also the small grains)
        """

        self.solid_fraction = self._chem.get_solid_heavies(self.T)
        self.chemistry_solid, self.chemistry_gas = self._chem.get_composition(self.T)
        sigma_g0 = self.sigma_g / (
            1.0 + (self._chem.dtg0 - self.solid_fraction)
        )  # get the fraction of the background gas

        self.sigma_dust_components = (
            self.chemistry_solid
            * self.solid_fraction[:, np.newaxis, np.newaxis]
            * sigma_g0[:, np.newaxis, np.newaxis]
        )
        self.sigma_g_components = (
            self.chemistry_gas * self.sigma_g[:, np.newaxis, np.newaxis]
        )

        self.sigma_dust = np.sum(self.sigma_dust_components[:, 1, :], axis=1)

        self.sigma_dust_components[:, 1, 0] = np.zeros_like(self.sigma_g)

        assert np.all(
            np.isclose(
                (self.sigma_g + self.sigma_dust), sigma_g0 * (1.0 + self._chem.dtg0)
            )
        ), "error in chemistry: The total masses do not sum up correctly"
        assert np.all(
            np.isclose(
                (
                    (
                        self.sigma_dust_components[:, 1, :]
                        + self.sigma_g_components[:, 1, :]
                    )
                    / (self.sigma_g + self.sigma_dust)[:, np.newaxis]
                ).std(axis=0),
                0,
            )
        ), "error in chemistry: The individual masses do not sum up to a constant value"
        assert np.all(
            np.isclose(
                np.sum(
                    self.sigma_dust_components[:, 1, :]
                    + self.sigma_g_components[:, 1, :],
                    axis=1,
                )
                / sigma_g0,
                (1 + self._chem.dtg0),
            )
        ), "error in chemistry: the dust to gas ratio is not correct"
        assert np.all(
            np.isclose(
                np.sum(
                    self.sigma_dust_components[:, 1, 1:]
                    + self.sigma_g_components[:, 1, 1:],
                    axis=1,
                )
                / sigma_g0,
                self._chem.dtg0,
            )
        ), "error in chemistry: the dust to gas ratio is not correct"

        self.sigmin_peb = 1e-60

        return

    def remove_from_gas(self, mdot, f_red, idx):
        """
        removes m_dot at cell boundaries left and right to planet

        Gas flows from right to left -> accrete from right, remove at cell of planet

        Parameters
        ----------
        mdot: numpy array of length (N+1)x(chem) that holds mdot for every interface
        f_red: reduction rate according to Camille
        idx: mask with indexes at which the mass is removed
        """
        try:
            sum_f_red = np.sum(f_red)
            delta_sigma = f_red / ((self.r_i[1:][idx] - self.r_i[:-1][idx]) * sum_f_red)
            if hasattr(self, "sigdot"):  # remove using visc disk evolution
                self.sigdot_gas_acc = np.zeros_like(self.sigdot)
                self.sigdot_gas_acc[idx] = delta_sigma[:, np.newaxis] * mdot[1]
            else:  # there is no visc disk evolution, remove directly instead
                sum_mdot_el = np.sum(mdot, axis=1)[0]
                self.sigma_g[idx] -= (
                    sum_mdot_el * delta_sigma / (2 * np.pi * self.r[idx]) * self.dt
                )
                self.sigma_g_components[idx] -= (
                    delta_sigma[:, np.newaxis, np.newaxis]
                    / (2.0 * np.pi * self.r[idx, np.newaxis, np.newaxis])
                    * mdot
                    * self.dt
                )
                self.sigma_g[self.sigma_g < self.sigmin] = self.sigmin
                self.sigma_g_components[
                    self.sigma_g_components < self.sigmin
                ] = self.sigmin
            return
        except IndexError:
            return

    def remove_from_peb(self, mdot, idx):
        """
        removes m_dot at cell boundaries left and right to planet
        Pebbles flows from right to left -> accrete from right, remove from right

        Parameters
        ----------
        mdot: numpy array of length (N+1) that holds mdot for every interface
        idx: position of right interface to planet
        """
        try:
            delta_sigma = 1.0 / (self.r_i[1:][idx] - self.r_i[:-1][idx])
            if hasattr(self, "sigdot"):
                self.sigdot_peb_acc = np.zeros_like(self.sigdot)
                self.sigdot_peb_acc[idx] = delta_sigma * mdot[1]  # fill sigdot
            else:
                mdot_sum = np.sum(mdot, axis=1)[0]
                self.sigma_dust[idx] -= (
                    delta_sigma / (2.0 * np.pi * self.r[idx]) * mdot_sum * self.dt
                )
                self.sigma_dust_components[idx] -= (
                    delta_sigma
                    / (2.0 * np.pi * self.r[idx, np.newaxis, np.newaxis])
                    * mdot
                    * self.dt
                )
                self.sigma_dust[self.sigma_dust < self.sigmin_peb] = self.sigmin_peb
                self.sigma_dust_components[
                    self.sigma_dust_components < self.sigmin_peb
                ] = self.sigmin_peb

            return
        except IndexError:
            return

    def apply_gap_profile(self, a_p, fsigma, FWHM):
        # fsigma = np.maximum(fsigma, 1e-7)
        self.alpha_factor = 1 - gaussian(self.r, 1.0 - fsigma, a_p, FWHM)
        self.interpolation_idx = np.abs(self.alpha_factor - 1.0) < 0.01
        self.interpolation_idx[0] = True
        alpha_factor_interface = 0.5 * (self.alpha_factor[1:] + self.alpha_factor[:-1])
        midinterface = np.abs(alpha_factor_interface - 1) < 0.01
        self.interpolation_idx_interface = np.array([True, *midinterface, True])

    def update(self, t, dt):
        """main update routine. This is a wrapper to update_parameters that does the physics of disk evolution."""
        self.t = 1.0 * t
        self.dt = 1.0 * dt
        if not self.conf_static:
            self.update_parameters()

    def update_parameters(self):
        """Update the physics of the disk. See inherited classes."""
        pass

    def compute_sizes(self):
        """
        Compute sizes according to Birnstiel2012

        Returns
        -------

        """
        DENSITY_SILICATES = 3.0
        DENSITY_ICE = 1.0

        if self._calculate_rho_solid:
            msil = np.sum(self.sigma_dust_components[:, 1, 8:], axis=-1)
            mice = np.sum(self.sigma_dust_components[:, 1, 1:8], axis=-1)
            # assert np.all(np.isclose(msil + mice, self.sigma_dust)), 'The internal ice densities do not add up correctly'
            self.rho_solid = (msil + mice) / (
                msil / DENSITY_SILICATES + mice / DENSITY_ICE + 1e-90
            )
            self.rho_solid = np.maximum(self.rho_solid, DENSITY_ICE)

        if self._calculate_u_frag:
            # Calculate the fragmentation velocity according to the Drazkowska 2017 approach: 1m/s if silicate, 10m/s if icy
            mice = np.sum(self.sigma_dust_components[:, 1, 1:8], axis=-1)
            self.u_frag = np.where(mice > 0.01 * self.sigma_dust, 1000, 100)

        a_frag = (
            self.const_f_f
            * 2.0
            / (3.0 * np.pi)
            * self.sigma_g
            / (self.alpha_frag * self.rho_solid)
            * self.u_frag**2.0
            / (self.c_s**2.0)
        )
        tau_grow = self.sigma_g / ((self.omega_k * self.sigma_dust) + 1e-300)

        dlnpdlnr = np.abs(self.grad_p)
        a_drift = (
            self.const_f_d
            * 2
            * self.sigma_dust
            / (np.pi * self.rho_solid)
            * self.v_k**2
            / self.c_s**2
            / dlnpdlnr
        )
        St_df = self.u_frag * self.v_k / (dlnpdlnr * self.c_s**2 * (1 - self.const_N))
        a_df = St_df * 2 * self.sigma_g / (np.pi * self.rho_solid)

        a_1 = np.min([a_drift, a_frag, a_df], axis=0)
        if hasattr(self, "t"):
            val = self.t / tau_grow
            val[val > 700.0] = 700.0  # stabilize numerically
            a_1 = np.min([a_1, self.a_0 * np.exp(val)], axis=0)
        else:
            a_1 = 1 * self.a_0

        stokes_factor = self.rho_solid / self.sigma_g * np.pi / 2.0

        self.a_1 = a_1 * np.ones_like(self.r)  # broadcasting if nescessary
        self.stokes_number_pebbles = a_1 * stokes_factor
        self.stokes_number_small = self.a_0 * stokes_factor

        self.f_m = np.where(
            np.argmin([a_drift, a_frag, a_df], axis=0) == 0,
            0.97 * np.ones_like(self.stokes_number_pebbles),
            0.75 * np.ones_like(self.stokes_number_pebbles),
        )

        self.stokes_number_df = 1.0 * St_df
        self.stokes_number_drift = stokes_factor * a_drift
        self.stokes_number_frag = stokes_factor * a_frag

    def compute_viscous_evolution(self):
        """
        Do the viscous evolution!

        depends on self.sigma_g, self.sigma_g_components -> init will be done, if not been provided

        Returns
        -------
        updates all relevant variables

        """
        if hasattr(self, "dt"):
            self.icelines_mol = self._chem.get_position_of_ice(self.T)

            # calculate Particle sizes according to Birnstiel 2012
            if not self.conf_static_stokes:
                self.compute_sizes()

            # calculate velocities of dust and gas:
            self.compute_interface_velocities()
            self.calc_sigdot(self.icelines_mol)
            self.calc_pebble_flux(noupdate=True)
            if self.conf_photoevaporation:
                self.calc_sigdot_photoevap()

            sigma_g_components_mol = (
                self.get_viscous_evolution_next_timestep_components(
                    self.sigma_g_components, self.dt
                )
            )

            self.compute_dust_next_timestep(self.dt)
            self.compute_sigma_g_from_mol(sigma_g_components_mol)

        else:  # init sigma_g_components
            self.compute_disktmid()
            dtg = self._chem.get_solid_heavies(self.T)
            self.chemistry_solid, self.chemistry_gas = self._chem.get_composition(
                self.T
            )
            self.mu = self._chem.mu
            self.sigma_g_components = (
                self.chemistry_gas * self.sigma_g[:, np.newaxis, np.newaxis]
            )

        self.compute_disktmid()

        if not self.conf_evaporation:
            self.chemistry_solid, self.chemistry_gas = self._chem.get_composition(
                self.T
            )
            self.mu = self._chem.mu

        self.sigma_dust_components = (
            self.chemistry_solid * self.sigma_dust[:, np.newaxis, np.newaxis]
        )
        self.sigma_g_components = (
            self.chemistry_gas * self.sigma_g[:, np.newaxis, np.newaxis]
        )

        if hasattr(self, "vr_peb"):
            self.calc_pebble_flux(noupdate=False)

        enclosed_mask = self.r < 2000.0 * AU

        self._old_disk_mass = 1.0 * self.disk_mass
        self._old_disk_mass_dust = 1.0 * self.disk_mass_dust

        self.disk_mass = np.trapz(
            2.0 * np.pi * self.sigma_g[enclosed_mask] * self.r[enclosed_mask],
            self.r[enclosed_mask],
        )
        self.disk_mass_dust = np.trapz(
            2 * np.pi * self.sigma_dust[enclosed_mask] * self.r[enclosed_mask],
            self.r[enclosed_mask],
        )

        # assert (self.disk_mass + self.disk_mass_dust) <= (
        #            self._old_disk_mass + self._old_disk_mass_dust), "Disk mass is growing"

        self.m_dot_components = self.compute_mdot_at_interfaces(self.sigma_g_components)
        self.m_dot = np.sum(self.m_dot_components[:, 0, :], axis=1)

        self.P = (
            self.c_s**2.0
            * self.sigma_g
            / (np.sqrt(2.0 * np.pi) * self.aspect_ratio * self.r)
        )

        self.grad_sig = np.gradient(np.log10(self.sigma_g), np.log10(self.r))
        self.grad_T = np.gradient(np.log10(self.T), np.log10(self.r))  # -0.5 for MMSN
        self.grad_p = np.gradient(np.log10(self.P), np.log10(self.r))

        self.eta = -0.5 * self.aspect_ratio**2 * self.grad_p

        self.solid_fraction = self.sigma_dust / self.sigma_g

        return

    def get_viscous_evolution_next_timestep_components(self, sigma_g_components, dt):
        """
        Do the gas viscous evolution.
        """
        #
        # Cast into diffusion equation form
        #

        nu = self.nu / self.alpha_factor
        x = self.r
        y = 2.0 * np.pi * x[:, np.newaxis] * sigma_g_components[:, 1, :]
        g = x[:, np.newaxis] ** 0.5 / nu[:, np.newaxis]
        d = 3.0 * nu[:, np.newaxis]
        v = np.zeros(y.shape)

        s = (
            1 * self.sigdot
        )  # if self.t < self.begin_photevap else np.zeros_like(self.sigdot)

        if hasattr(self, "sigdot_gas_acc"):
            s -= self.sigdot_gas_acc
            # checke das vorzeichen!

        if hasattr(self, "sigdot_photoevap"):
            s -= self.sigdot_photoevap
            # checke das vorzeichen!

        clip_val = (
            y
            - 2.0
            * np.pi
            * 0.9
            * x[:, np.newaxis]
            * self.sigmin
            * self.chemistry_gas[:, 1, :]
        ) / dt
        s = np.clip(s, -clip_val, clip_val)
        #
        # Set boundary conditions
        #
        bcl = (1, 0, 0, 0)
        # bcl = (0, 1, 2 * np.pi * x[0] * self.chemistry_gas[0, 1, :] * self.sigmin, 0)
        bcr = (0, 1, 2 * np.pi * x[-1] * self.chemistry_gas[-1, 1, :] * self.sigmin, 0)

        # Get the new value of y after one time step dt
        #
        y = solvediffonedee_components(x, y, v, d, g, s, bcl, bcr, dt=dt, int=False)
        #
        # Obtain new sigma
        #
        # y = np.maximum(y, 2 * np.pi * x[:, np.newaxis] * self.sigmin * self.chemistry_gas[:, 1, :])
        y = np.maximum(y, 2.0 * np.pi * x[:, np.newaxis] * self.sigmin)
        sigma = y / (2.0 * np.pi * x[:, np.newaxis])  # add offset for numerical reasons
        #
        # Return
        #
        return sigma

    def calc_pebble_flux(self, noupdate=False):
        """
        Calculates the pebble_flux at the cell center

        Returns
        -------

        """
        vr_peb = interp1d(
            self.r_i[1:-1], self.vr_peb, assume_sorted=True, fill_value="extrapolate"
        )(self.r)
        self.pebble_flux = np.abs(
            2.0 * np.pi * self.r * vr_peb * self.sigma_dust * self.f_m
        )
        if not noupdate:
            self.cum_pebble_flux += self.pebble_flux * self.dt

        return

    def calc_sigdot(self, icelines_mol):
        """
        Computes the sigdot that is needed for the viscous and dust evolution

        # checke das vorzeichen!

        Parameters
        ----------
        icelines_mol

        Returns
        -------

        """
        x = self.r
        y = 2.0 * np.pi * x[:, np.newaxis] * self.sigma_g_components[:, 1, :]
        self.sigdot = np.zeros(y.shape)
        # interpolate back to cell centers
        vr_dust = interp1d(
            self.r_i[1:-1], self.vr_dust, assume_sorted=True, fill_value="extrapolate"
        )(self.r)
        if not self.conf_evaporation:
            included_icelines = []
        else:
            # curate icelines: (only use icelines that are well enough inside of the grid)
            # included_icelines = np.arange(icelines_mol.shape[0])[
            #    np.logical_and(self.r[icelines_mol] > self.r[5], self.r[icelines_mol] < self.r[-5])]
            included_icelines = np.arange(1, icelines_mol.shape[0])  # exclude rest_mol

            sigdot_cond_const = (
                -3.0
                * self.epsilon_p
                / (2.0 * np.pi * self.rho_solid[:, np.newaxis])
                * y
                * self.sigma_dust[:, np.newaxis]
                * self.omega_k[:, np.newaxis]
                * np.sqrt(self._chem.mu[:, np.newaxis])
            )
            sigdot_cond_const *= (
                self.f_m[:, np.newaxis] / self.a_1[:, np.newaxis]
                + (1.0 - self.f_m[:, np.newaxis]) / self.a_0
            )
        for i in included_icelines:  # exclude rest_mol
            self.sigdot[icelines_mol[i] : -1, i] = sigdot_cond_const[
                icelines_mol[i] : -1, i
            ] / np.sqrt(np.squeeze(self._chem.M_array)[i])
            self.sigdot[: icelines_mol[i], i] = (
                2
                * np.pi
                * x[: icelines_mol[i]]
                * self.sigma_dust_components[: icelines_mol[i], 1, i]
                * np.abs(vr_dust[: icelines_mol[i]])
                / self._evap_width
            )

        min_density = (
            2.0
            * np.pi
            * x[:, np.newaxis]
            * 0.9
            * np.minimum(
                self.sigma_dust_components[:, 1, :], self.sigma_g_components[:, 1, :]
            )
        )

        self.sigdot = np.clip(
            self.sigdot, -min_density / self.dt, min_density / self.dt
        )  # do not transfer more then available in one timestep
        return

    def calc_sigdot_photoevap(self):
        """
        Wraps up calc_sigdot_photoevap_const and determins sigdot

        Returns
        -------

        """
        if self.t > self.begin_photevap:
            if hasattr(self, "tau_disk"):
                self.sigdot_photoevap = (
                    1.0
                    / self.tau_disk
                    * self.sigma_g_components[:, 1, :]
                    * 2.0
                    * np.pi
                    * self.r[:, np.newaxis]
                )
            else:
                self.sigdot_photoevap = (
                    np.zeros_like(self.r) * self.chemistry_gas[:, 1, :]
                )

    def get_dust_radial_drift_next_timestep(self, dt):
        """
        Advance the dust component one time step into the future.
        Radial drift and turbulent mixing included, as well as the
        gas drag as the gas is moving inward.

        ARGUMENTS:
        dt          = Time step in seconds
        alphamodel  = If True, then recompute self.nu from alpha-recipe (default)
        fixgas      = If True, then do *not* include the inward gas motion in dust drift.
        extracond   = (for special purposes only) List of extra internal conditions

        Note: If self.diskmodel.alphamix is present, then this alpha will be used (instead of the
        usual self.alpha) for the turbulent mixing.

        ** BEWARE: **
        Always make sure to have updated the midplane density and temperature,
        and then call the compute_stokes_from_agrain() before calling this subroutine,
        if you have evolved the gas beforehand.
        """
        #
        # Cast into diffusion equation form
        #
        # lmax = np.ones_like(self.r, dtype="bool")
        lmax = self.T < 2000.0
        x = self.r[lmax]
        y = (
            2.0 * np.pi * x[:, np.newaxis] * self.sigma_dust_components[lmax, 1, 1:]
        )  # Dust
        g = 2.0 * np.pi * x * self.sigma_g[lmax]  # Gas
        g = g[:, np.newaxis]
        d = self.nu[lmax]
        di = 0.5 * (d[1:] + d[:-1])[:, np.newaxis]
        vi = self.vr_dust[lmax[:-1], np.newaxis]

        s = -self.sigdot[
            lmax, 1:
        ]  # if self.t < self.begin_photevap else np.zeros_like(self.sigdot[:, 1:])

        if hasattr(self, "sigdot_peb_acc"):
            s -= self.sigdot_peb_acc[lmax, 1:]

        with np.errstate(under="ignore"):
            clip_val = (
                y
                - 2
                * np.pi
                * x[:, np.newaxis]
                * self.sigmin_peb
                * self.chemistry_solid[lmax, 1, 1:]
            ) / dt
            s = np.clip(s, -clip_val, clip_val)

            #
            # Set boundary conditions
            #
            bcl = (1, 0, 0, 1)
            # bcr = (1, 0, 0, 1)
            # bcl = (0, 1, 2 * np.pi * x[0] * self.chemistry_solid[0, 1, 1:] * self.sigmin_peb, 0)
            bcr = (
                0,
                1,
                2 * np.pi * x[-1] * self.chemistry_solid[-1, 1, 1:] * self.sigmin_peb,
                0,
            )
            #
            # Get the new value of y after one time step dt
            #
            y = solvediffonedee_components(
                x, y, vi, di, g, s, bcl, bcr, dt=dt, int=True, upwind=True
            )
            y = np.maximum(
                y,
                2
                * np.pi
                * x[:, np.newaxis]
                * self.chemistry_solid[lmax, 1, 1:]
                * self.sigmin_peb,
            )
            # y = np.maximum(y, 2 * np.pi * x[:, np.newaxis] * 1e-60)
            #
            # Obtain new sigdust
            #
            new_sigma_dust = y / (
                2 * np.pi * x[:, np.newaxis]
            )  # add offset for numerical reasons
        #
        # Return
        #
        sigma_dust = (
            np.ones_like(self.sigma_dust_components[:, 1, 1:]) * self.sigmin_peb
        )
        sigma_dust[lmax] = new_sigma_dust
        return sigma_dust

    def compute_dust_next_timestep(self, dt):
        """
        Wrapper routine around dust_radial_drift_next_timestep(), which stores
        the result. It also recomputes the midplane dust density. If updatestokes
        is False, then the Stokes number is not recomputed, but kept as it is in self.St.
        BEWARE: Make sure that the midplane gas density (self.rhomid) is current.
        NOTE: For computing Stokes number (and stopping time) we keep dv=1e3 for now.
        """
        # treat pebble drift:
        new_sigma_dust_mol = self.sigma_dust_components[:, 1, :]
        new_sigma_dust_mol[:, 1:] = self.get_dust_radial_drift_next_timestep(dt)

        # recalculate chemistry of pebbles:
        self.compute_sigma_dust_from_mol(new_sigma_dust_mol)

    def compute_sigma_g_from_mol(self, sigma_g_components_mol):
        self.sigma_g = np.sum(sigma_g_components_mol, axis=1)
        assert np.all(self.sigma_g > 0.0), "negative surface densities!"
        if self.conf_evaporation:
            self.chemistry_gas = self._chem.get_gas_composition(sigma_g_components_mol)
            self.mu = self._chem.mu

    def compute_sigma_dust_from_mol(self, sigma_dust_components_mol):
        if self.conf_evaporation:
            self.chemistry_solid = self._chem.get_solid_composition(
                sigma_dust_components_mol, self.T
            )

        self.sigma_dust = np.where(
            np.sum(self.chemistry_solid[:, 0], axis=1) > 0.0,
            np.sum(sigma_dust_components_mol, axis=1),
            np.ones_like(self.r) * self.sigmin_peb,
        )
        assert np.all(self.sigma_dust > 0), "negative surface densities!"

    def compute_mdot_at_interfaces(self, sigma_g, oned=False):
        """
        computes selfconsistent with viscous evolution m_dot at the interfaces
        """
        #
        # Cast into diffusion equation form
        #
        x = self.r
        if not oned:
            y = 2.0 * np.pi * self.r[:, np.newaxis, np.newaxis] * sigma_g
        else:
            y = 2.0 * np.pi * self.r * sigma_g
        g = self.r**0.5 / (self.nu / self.alpha_factor)
        d = 3.0 * self.nu / self.alpha_factor
        v = np.zeros(len(x))

        #
        # Compute the mdot
        #
        flux = -getfluxonedee(
            x, y, v, d, g, int=False, oned=oned
        )  # returns flux at interfaces
        return np.array([flux[0], *flux, flux[-1]])

    def compute_vr_at_interfaces(self):
        """
        Computes the radial gas velocity due to the viscous accretion at the interfaces.
        It does so by first calculating the gas Mdot at the interfaces, then dividing by
        2*pi*r*sigma. If upwind==True, then define the vr using the upwind scheme,
        otherwise using the average of 2*pi*r*sigma on both sides of the interface.

        """
        # Compute 2*pi*r*Sigma_g at the interfaces, using the correct averaging
        # or upwinding (note: Mdot has minus sign).
        #
        q = np.zeros(len(self.r) - 1)
        q += (self.m_dot[1:-1] <= 0) * 2.0 * np.pi * self.r[:-1] * self.sigma_g[:-1]
        q += (self.m_dot[1:-1] > 0) * 2.0 * np.pi * self.r[1:] * self.sigma_g[1:]

        #
        # Now compute the radial velocity from the mdot
        #
        self.vr_gas = -self.m_dot[1:-1] / (q + 1e-90)

    def calc_velocity_from_stokes(self, stokes, dlnpdlnr, csi, vki, vrgas):
        """
        Helper function that returns dust velocity from gas velocity and stokes_number
        Parameters
        ----------
        stokes
        dlnpdlnr
        csi
        vki
        vrgas

        Returns
        -------

        """
        stokes_i = (stokes[1:] + stokes[:-1]) / 2.0 + 1e-60
        return vrgas / (1.0 + stokes_i**2) + dlnpdlnr * csi**2 / vki / (
            stokes_i + 1.0 / stokes_i
        )

    def compute_interface_velocities(self):
        """
        Compute the radial drift velocity of the dust at the interfaces. We follow
        Birnstiel 2012
        """

        # Compute v_K and cs at interfaces
        #
        vki = 0.5 * (self.v_k[1:] + self.v_k[:-1])
        csi = 0.5 * (self.c_s[1:] + self.c_s[:-1])

        self.compute_vr_at_interfaces()
        vrgas = self.vr_gas
        #
        # Compute the dln(p)/dln(r) at the interfaces by calling the method
        # that computes omega
        #

        self.P = (
            self.c_s**2.0
            * self.sigma_g
            / (np.sqrt(2 * np.pi) * self.aspect_ratio * self.r)
        )
        self.grad_p = np.gradient(
            np.log10(self.P), np.log10(self.r)
        )  # -2.75 for MMSN, checked

        dlnpdlnr = 0.5 * (self.grad_p[1:] + self.grad_p[:-1])

        #
        # Compute the dust radial drift velocity including the passive advection by the gas
        #
        if not self.conf_static_stokes:
            f_m_i = (self.f_m[1:] + self.f_m[:-1]) / 2.0
            self.vr_peb = self.calc_velocity_from_stokes(
                self.stokes_number_pebbles, dlnpdlnr, csi, vki, vrgas
            )
            self.vr_small = self.calc_velocity_from_stokes(
                self.stokes_number_small, dlnpdlnr, csi, vki, vrgas
            )
            self.vr_dust = f_m_i * self.vr_peb + (1.0 - f_m_i) * self.vr_small
        else:
            self.vr_dust = self.calc_velocity_from_stokes(
                self.stokes_number_pebbles, dlnpdlnr, csi, vki, vrgas
            )
            self.vr_small = 0.0
            self.vr_peb = 1.0 * self.vr_dust

    def init_T(self):
        flux = self.lstar / (4.0 * np.pi * self.r**2.0)
        self.qirr = (
            flux * self.flang
        )  # Factor 2 for two sides, factor 0.5 for half irradiated down
        self.T_irr = (
            0.5 * self.qirr / sr
        ) ** 0.25  # Factor 0.5 because cooling is two-sided
        self.T_irr = np.maximum(self.T_irr, 10.0)
        self.T = 1.0 * self.T_irr

        self._T_init_finished = False
        self.compute_cs_and_hp()
        self.compute_nu()

    def compute_disktmid(self):
        """
        Compute the midplane temperature based on the simple
        irradiation recipe in which half of the irradiated
        energy is radiated away, half is heating the midplane
        of the disk. Will compute self.tmid, the midplane
        temperature, self.cs, the isothermal sound speed,
        and self.hp, the vertical scale height.

        If vischeat==True, then also include the viscous
        heating. But for that we need sigma and mean_opacity_rosseland.

        NOTE: Since this method uses the Rosseland mean _0pacity array
              self.mean_opacity_rosseland[:], if you are not sure if this
              is still up-to-date with the current surface density and
              temperature, then better call self.compute_mean_opacity()

        If keeptvisc==True but vischeat==False, then we do
        not calculate the self.tvisc, but we keep the
        self.tvisc (i.e. we do not put it to zero).

        COMPUTES:
        self.tmid   = Midplane temperature in [K]
        self.cs     = Isothermal sound speed [cm/s]
        self.hp     = Pressure scale height [cm]
        self.qirr   = Irradiative heating rate of the disk [erg/cm^2/s]
        self.qvisc  = (if vischeat==True) Viscous heating rate [erg/cm^2/s]
        """
        if hasattr(self, "alpha_factor"):
            alpha = self.alpha / self.alpha_factor
        else:
            alpha = self.alpha

        if hasattr(self, "T_irr"):
            if self.conf_temp_evol or not self._T_init_finished:
                self.T = solve_viscous_heating_globally(
                    sigma_g=self.sigma_g,
                    tirr=self.T_irr,
                    T=self.T,
                    omega_k=self.omega_k,
                    alpha=alpha,
                    Z=self.DTG_small_grains,
                    r=self.r,
                    mu=self._chem.mu,
                    interpolation_idx=None
                    if not hasattr(self, "interpolation_idx")
                    else self.interpolation_idx,
                )

                self.T_visc = (self.T**4.0 - self.T_irr**4.0) ** 0.25
        else:
            self.init_T()

        self.compute_cs_and_hp()
        self.compute_nu()
        self.calc_opacity()

    def calc_opacity(self):
        """
        calculate the gas _0pacity - caution value in cgs
        metalicity is dependent of T! However, this is to good approximation constant -> not updated here.
        If T>2000: No dust -> _0pacity = 0, fix by adding min value to metalicity
        """
        self.opacity = opacity_fct(self.rho_g, self.T, Z=self.DTG_small_grains)

    def compute_cs_and_hp(self):
        """update all quantities related to T"""
        self.c_s = (k_B * self.T / (self._chem.mu * m_p)) ** 0.5
        self.aspect_ratio = self.c_s / (self.omega_k * self.r)
        if hasattr(self, "sigma_g"):
            self.rho_g = self.sigma_g / (
                np.sqrt(2 * np.pi) * self.aspect_ratio * self.r
            )

    def compute_nu(self):
        """compute the viscosity."""
        self.nu = self.alpha * self.c_s * self.c_s / self.omega_k
