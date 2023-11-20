import os

import matplotlib.pyplot as plt
import numpy as np
import tables
from astropy import constants as const
from astropy import units as u
from scipy.interpolate import interp1d

from chemcomp.helpers import OUTPUT_PATH, eval_kwargs

M_jup = const.M_jup.cgs.value
G = const.G.cgs.value
M_earth = const.M_earth.cgs.value
MS = const.M_sun.cgs.value
year = 3600 * 24 * 365.25
Myr = 1e6 * year
AU = const.au.cgs.value


class DataObject(object):
    def __init__(self, planet, output):
        self.output = output
        name = output.get("name", "planet")
        self.conf_save_disk = output.get("save_disk", False)
        self.conf_save_disk_interval = output.get("save_disk_interval", 100)
        self._save_counter = 0
        filename = output.get("directory", os.path.join(OUTPUT_PATH, "{}.h5".format(name)))
        if os.path.exists(filename):
            os.remove(filename)
        self._h5file = tables.open_file(filename, mode='w', title="Output of {}".format(name))
        self._h5file.set_node_attr("/", "finished", 0)
        self.output_disk = self._h5file.create_group("/", 'disk', 'Output of the Disk')
        self.output_planet = self._h5file.create_group("/", 'planet', 'Output of the Planet')
        self.output_acc = self._h5file.create_group("/", 'acc', 'Output of accretion')

        self.planet = planet
        self.disk = planet.disk

        list_of_disk_quantities = output.get("quantities",
                                             ["a_1", "t", "sigma_g", "sigma_dust", "sigma_dust_components", "T",
                                              "sigma_g_components",
                                              "f_m", "stokes_number_pebbles", "stokes_number_df", "stokes_number_drift",
                                              "stokes_number_frag", "stokes_number_small", "pebble_flux",
                                              "cum_pebble_flux",
                                              "peb_iso", "peb_iso_pi_crit", "T_visc", "T_irr", "mu", "m_dot", "vr_dust", "vr_gas", "m_dot_components"])

        self.dict_of_planet_units = {"t": "u.s", "M": "u.g", "M_a": "u.g", "M_c": "u.g", "a_p": "u.cm",
                                     "T": "u.K", "comp_a": "u.g", "comp_c": "u.g", "tau_m": "u.s",
                                     "sigma_g": "u.g/u.cm**2",
                                     "gamma_tot": "u.g*u.cm**2/u.s**2", "sigma_peb": "u.g/u.cm**2",
                                     "pebble_flux": "u.g/u.s", "peb_iso": "u.g"}
        self.dict_of_acc_units = {"m_dot": "u.g/u.s", "m_dot_a": "u.g/u.s", "m_dot_c": "u.g/u.s",
                                  "m_dot_a_chem": "u.g/u.s", "m_dot_c_chem": "u.g/u.s", "M_z": "u.g"}

        self.list_of_planet_quantities = ["t", "M", "M_a", "M_c", "a_p", "T", "comp_a", "comp_c", "tau_m",
                                          "gamma_tot", "sigma_g", "gamma_norm", "fSigma", "sigma_peb", "pebble_flux",
                                          "peb_iso"]
        self.list_of_acc_quantities = ["m_dot", "m_dot_a", "m_dot_c", "m_dot_a_chem", "m_dot_c_chem", "M_z", "regime"]

        # trim list of quantities to existing:
        _, self.list_of_disk_quantities = [], []
        [self.list_of_disk_quantities.append(quant) if hasattr(self.disk, f"{quant}") else _.append(quant) for quant in
         list_of_disk_quantities]

        self._h5file.create_earray(self.output_disk, "r", obj=self.disk.r)
        self._h5file.create_earray(self.output_disk, "r_i", obj=self.disk.r_i)

        self.list_of_disk_pointer = []
        for quantity in self.list_of_disk_quantities:
            obj = np.array(getattr(self.disk, quantity))[np.newaxis]
            self.list_of_disk_pointer.append(self._h5file.create_earray(self.output_disk, quantity, obj=obj))

        self.list_of_planet_pointer = []
        for quantity in self.list_of_planet_quantities:
            obj = np.array(getattr(self.planet, quantity))[np.newaxis]
            self.list_of_planet_pointer.append(self._h5file.create_earray(self.output_planet, quantity, obj=obj))

        acc_name_dict = {"PebbleAccretion": "peb", "GasAccretion": "gas"}
        # drop acc model if not used:
        self.acc_index_list = []
        self.acc_models = {}
        for i in range(len(self.planet.accretion_models)):
            model = self.planet.accretion_models[i]
            if model.name != "DefaultAccretion":
                self.acc_models[i] = acc_name_dict.get(model.name)
                self.acc_index_list.append(i)

        self.list_of_acc_pointer = {}
        for i in self.acc_index_list:
            self.list_of_acc_pointer[i] = []
            for quantity in self.list_of_acc_quantities:
                obj = np.array(getattr(self.planet.accretion_models[i], quantity))[np.newaxis]
                self.list_of_acc_pointer[i].append(
                    self._h5file.create_earray(self.output_planet, "{}_{}".format(quantity, self.acc_models[i]),
                                               obj=obj))

        self._planet_units = self._h5file.create_group("/planet", 'units', 'units')
        for k, v in self.dict_of_planet_units.items():
            self._h5file.set_node_attr("/planet/units", k, v)

        for _, acc in self.acc_models.items():
            for k, v in self.dict_of_acc_units.items():
                self._h5file.set_node_attr("/planet/units", f"{k}_{acc}", v)

    def save_disk(self):
        for i in range(len(self.list_of_disk_pointer)):
            array = getattr(self.disk, self.list_of_disk_quantities[i])
            array = np.array(array)
            self.list_of_disk_pointer[i].append(array[np.newaxis])

    def save_planet(self):
        for i in range(len(self.list_of_planet_pointer)):
            array = getattr(self.planet, self.list_of_planet_quantities[i])
            array = np.array(array)
            self.list_of_planet_pointer[i].append(array[np.newaxis])

    def save_acc(self):
        for j in self.acc_index_list:
            ptr = self.list_of_acc_pointer[j]
            for i in range(len(ptr)):
                array = getattr(self.planet.accretion_models[j], self.list_of_acc_quantities[i])
                array = np.array(array)
                self.list_of_acc_pointer[j][i].append(array[np.newaxis])

    def save(self, force=False):
        if self.planet.t >= self.planet.save_interval * self._save_counter or force:
            if self.conf_save_disk and self._save_counter % self.conf_save_disk_interval == 0 or force:
                self.save_disk()
            self.save_planet()
            self.save_acc()
            self._save_counter += 1

    def finish(self, error=None):
        if error is None:
            self._h5file.set_node_attr("/", "finished", 1)
        self._h5file.close()


class Planet(object):
    """ docstring for Planet. """

    def __init__(
            self, output, rho_c, a_p, t_0, disk, accretion_models=[], M_c=None, M_a=None, **kwargs
    ):
        self.save_interval = output.get("save_interval", 50000 * u.yr).cgs.value
        self._print_number = 0
        self._total_save_counter = 0
        self.t = 1.0 * disk.t_0
        self.t_0 = t_0.cgs.value + disk.t_0
        self.disk = disk

        self.rho_c = 1.0 * rho_c.cgs.value  # density of the core
        self.a_p = 1.0 * a_p.cgs.value  # semi major axis

        M0_fact = eval_kwargs(kwargs.get("M0_fact", 1.0))
        # Start with pebiso:
        pebiso_start = eval_kwargs(kwargs.get("pebiso_start", False))
        # Fraction of solids accreted into envelope during solid accretion:
        self.solid_frac = eval_kwargs(kwargs.get("solid_frac", 0.1))
        if M_c is None:
            # use Lambrecht2012 transitionmass
            # -> ignored if pebiso_start
            M_c = (1-self.solid_frac)* (
                    np.sqrt(1.0 / 3.0)
                    * (disk.eta * disk.v_k) ** 3.0
                    / (const.G.cgs.value * disk.omega_k)
            )
            M_c = M0_fact * np.interp(a_p.cgs.value, disk.r, M_c) * u.g

        if M_a is None:
            M_a = 0.0 * u.g
            if self.solid_frac != 0.0:
                M_a =  self.solid_frac / (1.0 - self.solid_frac) * M_c

        self.M_c = 1.0 * M_c.cgs.value  # core mass
        self.M_a = 1.0 * M_a.cgs.value  # atmospheric mass
        self.M = self.M_a + self.M_c
        self.accretion_models = accretion_models  # list containing the accretion objects
        self._accretion_models_dict = {acc.name: acc for acc in self.accretion_models}
        self._accrete = eval_kwargs(kwargs.get("accrete", True))
        self.matter_removal = eval_kwargs(kwargs.get("matter_removal", True))
        self.force_insitu = eval_kwargs(kwargs.get("force_insitu", False))
        self.r_in = eval_kwargs(kwargs.get("r_in", self.disk.r_i[0]*u.cm).cgs.value) # where to stop the planet formation

        self.name = output.get("name", "planet")
        self.output_file = output.get("directory", os.path.join(OUTPUT_PATH, "{}.h5").format(self.name))
        self.tau_m = np.inf  # migration rate
        self.fSigma = 1.0    # gap depth

        self.M_end = eval_kwargs(kwargs.get("M_end", 1.0 * u.M_sun).cgs.value)   # final mass at which we interrupt the formation of the planet

        self.past_pebble_iso = False  # boolean flag that stores if we passed the pebiso mass
        self._keep_peb_iso = eval_kwargs(kwargs.get("keep_peb_iso", True))  # keep pebble isolation mass constant when pebble isolation mass is reached (otherwise we might have multiple phases of pebble accretion)
        self._use_pebiso_diffusion = eval_kwargs(kwargs.get("use_pebiso_diffusion", True))  # use the diffusion part of the pebble isolation mass

        self._init_accretion()  # initialise the accretion models

        self.update(second=True)  # first update of all quantities
        self.chem_mask_array = self.disk._chem.mask_array != 0  # see mask_array in chemistry class (needed since there are less elements than mol. species)
        self.comp_a = self.chemistry_solid * self.M_a  # molecular and elemental composition of the planetary envelope
        self.comp_c = self.chemistry_solid * self.M_c  # molecular and elemental composition of the planetary core

        self.output = DataObject(self, output)   # Initialize the IO class
        self.evolve_disk_init()   # start evolving the disk to the point where the planet is placed into the disk

        if pebiso_start:  # if planet is seeded with M=M_iso
            self.calc_pebble_iso()
            print("{}: starting at peb_iso with {:.2e} M_e".format(self.name,(self.peb_iso*u.g).to(u.earthMass)))
            self.M_c = (1.0-self.solid_frac) * self.peb_iso
            self.M_a = self.solid_frac * self.peb_iso
            self.M = self.M_a + self.M_c
            self.past_pebble_iso = True

        self.update(second=True)
        self.comp_a = self.chemistry_solid * self.M_a
        self.comp_c = self.chemistry_solid * self.M_c


    def _init_accretion(self):
        """ Initialise the Accretionmodels. """
        for acc in self.accretion_models:
            acc.init_accretion(self)

    def update_torque(self):
        """ Calculate the Torques"""
        pass

    def update_a_p(self):
        """ Migrate the Planet. """
        pass

    def update(self, second, interpolation_idx=None, interpolation_idx_interface=None):
        """
        Main update routine that ensures that all planetary quantities are interpolated and uptodate.

        Parameters that are dependent on planet Mass need to be kept updated
        second may only be called once after migration

        interpolation_idx can be used if gap should be skipped during the calculation of sigma_gap relevant quantities
        """
        if interpolation_idx is None:  # indexes as which we want to interpolate the disk quantities to the planetary position
            interpolation_idx = np.ones_like(self.disk.r,dtype=bool)
        if interpolation_idx_interface is None:
            interpolation_idx_interface = np.ones_like(self.disk.r_i,dtype=bool)

        self.sigma_g = (  # surface density of the gas phase
            interp1d(self.disk.r[interpolation_idx], self.disk.sigma_g[interpolation_idx], fill_value="extrapolate", assume_sorted=True)(self.a_p)
        )
        self.chemistry_gas = interp1d(self.disk.r, self.disk.chemistry_gas, axis=0, fill_value="extrapolate",
                                      assume_sorted=True)(self.a_p)  # composition of the gas phase
        if self.t < self.disk.begin_photevap:
            assert np.all(self.chemistry_gas >= 0.0), "gas chemistry is negative!"
        elif not np.all(self.chemistry_gas >= 0.0):
            self.chemistry_gas = np.maximum(self.chemistry_gas, np.zeros_like(self.chemistry_gas))   # composition of the gas phase

        if self.t < self.disk.begin_photevap:
            assert self.sigma_g >= 0.0, "negative surface density!"
        else:
            self.sigma_g = np.maximum(self.sigma_g, 1e-60)    # surface density

        self.sigma_g_components = self.sigma_g * self.chemistry_gas    # surface density as a compositional vector

        self.nu = interp1d(self.disk.r, self.disk.nu, fill_value="extrapolate", assume_sorted=True)(self.a_p)  # viscosity

        self.m_dot_disk = (  # viscous accretion rate of the gas phase of the disk
            interp1d(self.disk.r_i[interpolation_idx_interface], np.abs(self.disk.m_dot[interpolation_idx_interface]),fill_value="extrapolate", assume_sorted=True)(self.a_p)
        )

        if self.m_dot_disk <= 0.0:
            # fall back to static mdot if everything else fails
            self.m_dot_disk = 3.0 * np.pi * self.nu * self.sigma_g

        self.m_dot_disk_components = self.chemistry_gas * self.m_dot_disk

        if self.t < self.disk.begin_photevap:
            assert self.m_dot_disk >= 0.0, "the absolute value of m_dot_disk is not positive! Interpolation error"
        else:
            self.m_dot_disk_components = np.maximum(self.m_dot_disk_components, np.ones_like(self.chemistry_gas)*1e-60)
            self.m_dot_disk = np.maximum(1e-60, self.m_dot_disk)

        self.sigma_g_pert = (   # perturbed surface density of the gap
            interp1d(self.disk.r, self.disk.sigma_g, fill_value="extrapolate", assume_sorted=True)(self.a_p)
        )
        self.sigma_g_components_pert = interp1d(self.disk.r,  # perturbed surface density of the gap as a compositional vector
                                           self.disk.sigma_g_components, axis=0,
                                           fill_value="extrapolate", assume_sorted=True)(
            self.a_p)
        self.T = interp1d(self.disk.r, self.disk.T, fill_value="extrapolate", assume_sorted=True)(self.a_p)  # Temperature
        self.H_g = interp1d(self.disk.r, self.disk.aspect_ratio, fill_value="extrapolate", assume_sorted=True)(self.a_p) * self.a_p   # Disk scalehight
        self.c_s = interp1d(self.disk.r, self.disk.c_s, fill_value="extrapolate", assume_sorted=True)(self.a_p)  # soundspeed
        self.rho_g = self.sigma_g / (np.sqrt(2.0 * np.pi) * self.H_g)  # gas volumedensity
        self.grad_sig = interp1d(self.disk.r[interpolation_idx], self.disk.grad_sig[interpolation_idx], fill_value="extrapolate", assume_sorted=True)(self.a_p)   # gradient of the surfacedensity
        self.grad_T = interp1d(self.disk.r, self.disk.grad_T, fill_value="extrapolate", assume_sorted=True)(self.a_p)   # gradient of the Temperature
        self.grad_p = interp1d(self.disk.r, self.disk.grad_p, fill_value="extrapolate", assume_sorted=True)(self.a_p)   # gradient of the pressure
        self.r_h = self.a_p * (self.M / (3.0 * self.disk.M_s)) ** (1.0 / 3.0)  # hill radius
        self.v_k = np.sqrt(G * self.disk.M_s / self.a_p)  # keppler velocity
        self.omega_k = self.v_k / self.a_p  # keppler orbital velocity
        self.del_v = interp1d(self.disk.r, self.disk.eta, fill_value="extrapolate", assume_sorted=True)(
            self.a_p) * self.v_k  # difference between keppler and real gas orbital speed
        self.r_b = (G * self.M) / (self.del_v ** 2.0)  # bondi radius
        self.v_h = self.r_h * self.omega_k  # hill velocity

        self.P = interp1d(self.disk.r, self.disk.P, fill_value="extrapolate", assume_sorted=True)(
            self.a_p)

        if second:
            # This may only be called after migration has been applied
            # migration calculates the gap parameter according to Kanagawa/Ndugu -> self.fsigma

            self.pebble_flux = interp1d(self.disk.r, self.disk.pebble_flux, fill_value="extrapolate", assume_sorted=True)(self.a_p)  # pebble flux
            if hasattr(self.disk, "f_m"):
                self.sigma_peb = interp1d(self.disk.r, self.disk.sigma_dust * self.disk.f_m, fill_value="extrapolate", assume_sorted=True)(self.a_p)  # pebble surface density
                self.sigma_peb_components = interp1d(self.disk.r,  # pebble surface density as a compositional vector
                                                     self.disk.sigma_dust_components * self.disk.f_m[:, np.newaxis,
                                                                                       np.newaxis], axis=0,
                                                     fill_value="extrapolate",
                                                     assume_sorted=True)(self.a_p)

            else:
                self.sigma_peb = interp1d(self.disk.r, self.disk.sigma_dust, fill_value="extrapolate", assume_sorted=True)(self.a_p)  # checked
                self.sigma_peb_components = interp1d(self.disk.r,
                                                     self.disk.sigma_dust_components, axis=0,
                                                     fill_value="extrapolate",
                                                     assume_sorted=True)(self.a_p)

            self.calc_pebble_iso()

            # idx_i_r, idx_i_l and idx:
            # centers:        i-1  | i | i+1
            # interfaces:          i   i+1
            # idx_i_r: interface i+1, idx_i_l: interface i
            # idx: cell center
            # tested -> Check

            self.idx_i_r = np.searchsorted(self.disk.r_i, self.a_p,
                                           side="right")  # interface i+1: idx of right cell interface in self.disk.r_i
            self.idx_i_l = max(self.idx_i_r - 1, 0)  # interface i: idx of left cell interface in self.disk.r_i
            self.idx = 1 * self.idx_i_l  # center i: idx of cell (self.disk.r) planet is located in

            self.chemistry_solid = interp1d(self.disk.r, self.disk.chemistry_solid, axis=0, fill_value="extrapolate",
                                            assume_sorted=True)(self.a_p)   # compositional vector of pebbles

            self.update_parameters()

    def update_parameters(self):
        """ Parameters that are dependent on planet Mass need to be kept updated """
        pass

    @property
    def r_c(self):
        """
        radius of the core
        """
        return (3.0 / 4.0 * self.M_c / (np.pi * self.rho_c)) ** (1.0 / 3.0)  # core radius

    def calc_pebble_iso(self):
        """ calculate the pebble isolation mass."""
        if self._keep_peb_iso and self.past_pebble_iso:
            return

        self.disk.calc_pebble_iso(diffusion=self._use_pebiso_diffusion)

        if hasattr(self.disk,"interpolation_idx"):
            self.peb_iso = interp1d(self.disk.r[self.disk.interpolation_idx], self.disk.peb_iso[self.disk.interpolation_idx], fill_value="extrapolate", assume_sorted=True)(self.a_p)
        else:
            self.peb_iso = interp1d(self.disk.r, self.disk.peb_iso, fill_value="extrapolate", assume_sorted=True)(self.a_p)

        self.past_pebble_iso = True if self.M > self.peb_iso else False
        return

    def update_all(self):
        """ Wrapper that handles the complete update of disk and planet as well as migration."""
        self.disk.update(self.t, self.h)
        self.update(second=False)
        self.update_a_p()
        self.update(second=True)  # update all quantities effected by a_p

    def evolve_disk_init(self):
        """ Evolve the disk before the protoplanetary seed is placed into the disk. """
        self.disk.evolve_init()

        # initialise disk and save
        self.disk.compute_interface_velocities()
        self.h = 0.0
        self.calc_pebble_iso()
        self.print_params(force=True)

        while self.t < self.t_0 and self.keep_running:
            self.disk.compute_interface_velocities()
            self.calc_timestep()
            self.disk.update(self.t, self.h)
            self.calc_pebble_iso()
            self.print_params()

    def calc_timestep(self):
        """ Fix the timestep and update time. Important: This is a fixed timestep currently!"""
        if self.disk._timestep_input:
            self.h = self.disk._timestep_input
        else:
            self.h = 10 * year

        self.t += self.h

    @property
    def keep_running(self):
        """Stopping conditions"""
        if self.disk.disk_mass < 1e-6 * MS:
            return False
        if self.t > self.disk.t_end:
            return False
        if self.a_p < self.r_in and self._migrate:
            return False
        if self.M > self.M_end:
            return False
        return True

    def grow_mass(self) -> None:
        """ Main loop and supervisor of chemcomp.

        dt: the total time that the model should be evolved.
            Needs to be a astropy Quantity
        """

        def increment_mass(acc):
            """ increment the mass of the planet for a certain accretion model """
            acc.calc_m_dot()
            acc.update_z()
            self.total_m_dot += acc.m_dot
            self.comp_a += acc.m_dot_a_chem * self.h
            self.comp_c += acc.m_dot_c_chem * self.h
            self.M_a += acc.m_dot_a * self.h
            self.M_c += acc.m_dot_c * self.h
            if self.matter_removal:
                acc.remove_mass_from_flux()

        def _grow():
            """ wrap around increment_mass and update final mass. """
            self.total_m_dot = 0.0
            list(map(increment_mass, self.accretion_models))
            self.M = self.M_a + self.M_c

        self.print_params()  # print parameters also at beginning
        self.disk.compute_interface_velocities()
        self.calc_timestep()
        self.update_all()  # update first time
        try:
            # main loop
            while self.keep_running:
                self.update_all()  # update for next iteration
                if self.keep_running and self._accrete:
                    _grow()

                self.print_params()

                # for the next timestep:
                self.disk.compute_interface_velocities()
                self.calc_timestep()

            if self.force_insitu and self.a_p <= self.r_in:
                # grow without migration at the inner disk edge
                self._migrate = False
                while self.keep_running:
                    self.update_all()  # update for next iteration
                    if self.keep_running and self._accrete:
                        _grow()

                    # for the next timestep:
                    self.disk.compute_interface_velocities()
                    self.calc_timestep()

            self.print_params(force=True)
            self.output.finish()

        except Exception as e:
            # Handle error and print error
            self.print_params(force=True)
            self.output.finish(error=e)
            print(f"Execution ended with error: {e}")

        return

    def print_params(self, force=False):
        """" Handle the IO using self.output -> DataObject.
        Saves !after! accretion applied
        """
        self.output.save(force)

        # printing some diagnostics to stdio during runtime
        if self.t > self._print_number * 0.5 * Myr or force:
            try:
                CO = self.comp_a[0, 0] / self.comp_a[0, 1] * 16.0 / 12.0
            except FloatingPointError:
                CO = 0.0
            print("{}: t={:.1f}Myr, dt={:.1f}yr, M={:.1f}M_e, r={:.1f}AU, CO={:.2f}".format(self.name, self.t / Myr,
                                                                                            self.h / year,
                                                                                            self.M / M_earth,
                                                                                            self.a_p / AU, CO))
            self._print_number += 1
