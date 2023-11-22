import astropy.constants as const
import astropy.units as u
import h5py
import matplotlib.colors as colors
import numpy as np
import parse
import tables
from adjustText import adjust_text
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

np.seterr(divide="ignore", invalid="ignore")

elements = ["C", "O", "Fe", "S", "Mg", "Si", "Na", "K", "N", "Al", "Ti", "V", "H", "He"]
molecules = [
    "rest",
    "CO",
    "N2",
    "CH4",
    "CO2",
    "NH3",
    "H2S",
    "H2O",
    "Fe3O4",
    "C",
    "FeS",
    "NaAlSi3O8",
    "KAlSi3O8",
    "Mg2SiO4",
    "Fe2O3",
    "VO",
    "MgSiO3",
    "Al2O3",
    "TiO",
]

FeH_dict = {
    "FeH": [-0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4],
    "SiH": [-0.32, -0.25, -0.17, -0.09, 0.0, 0.1, 0.2, 0.32, 0.41],
    "SH": [-0.32, -0.25, -0.17, -0.09, 0.0, 0.1, 0.2, 0.32, 0.41],
    "MgH": [-0.26, -0.19, -0.1, -0.01, 0.08, 0.19, 0.3, 0.42, 0.53],
    "CH": [-0.29, -0.22, -0.14, -0.06, 0.03, 0.11, 0.2, 0.28, 0.36],
    "OH": [-0.19, -0.13, -0.08, -0.02, 0.03, 0.08, 0.13, 0.19, 0.23],
}

OHabu0 = 4.90e-4  # O/H abundance, solar value
CHabu0 = 2.69e-4  # C/H abundance, solar value
SiHabu0 = 3.24e-5  # Si/H adundance, solar value
FeHabu0 = 3.16e-5  # Fe/H adundance, solar value
SHabu0 = 1.32e-5  # S/H adundance, solar value
MgHabu0 = 3.98e-5  # Mg/H adundance, solar value
NHabu0 = 6.76e-5
HeHabu0 = 0.085  # mass ratio of H in protoplanetary disk
AlHabu0 = 2.82e-6
TiHabu0 = 8.91e-8
KHabu0 = 1.07e-7
NaHabu0 = 1.74e-6
VHabu0 = 8.59e-9

MassO = 16.0
MassH = 1.0
MassFe = 56.0
MassMg = 24.3
MassSi = 28.0
MassS = 32.0
MassHe = 4.0
MassTi = 47.867
MassAl = 27
MassK = 39.0983
MassNa = 23
MassN = 14
MassC = 12.0
MassV = 50.9415

el_masses = {
    "C": MassC,
    "O": MassO,
    "H": MassH,
    "Fe": MassFe,
    "Mg": MassMg,
    "Si": MassSi,
    "N": MassN,
    "S": MassS,
    "He": MassHe,
    "Al": MassAl,
    "Ti": MassTi,
    "K": MassK,
    "Na": MassNa,
    "V": MassV,
}

DOT_SPACE = 1e5 * u.yr  # time in years between dots
DOT_SPACE_LARGE = 5e5 * u.yr  # time in years between dots

linestyle_tuple = [
    ("loosely dotted", (0, (1, 10))),
    ("dotted", (0, (1, 1))),
    ("densely dotted", (0, (1, 1))),
    ("loosely dashed", (0, (5, 10))),
    ("dashed", (0, (5, 5))),
    ("densely dashed", (0, (5, 1))),
    ("loosely dashdotted", (0, (3, 10, 1, 10))),
    ("dashdotted", (0, (3, 5, 1, 5))),
    ("densely dashdotted", (0, (3, 1, 1, 1))),
    ("dashdotdotted", (0, (3, 5, 1, 5, 1, 5))),
    ("loosely dashdotdotted", (0, (3, 10, 1, 10, 1, 10))),
    ("densely dashdotdotted", (0, (3, 1, 1, 1, 1, 1))),
]

dict_of_planet_read_names = {
    "t": r"t",
    "M": r"M",
    "M_a": r"M_a",
    "M_c": r"M_c",
    "a_p": r"a_p",
    "r_p": r"R_p",
    "tau_m": r"\tau_m",
    "T": "T",
    "comp_a": "Atmospheric Content",
    "comp_c": "Core Content",
    "sigma_g": r"\Sigma_g",
    "gamma_tot": r"\Gamma",
    "gamma_norm": r"\Gamma / \Gamma_0",
    "fSigma": "fSigma",
}
dict_of_acc_names_template = {
    r"m_dot": r"\dot M_{\mathrm{",
    r"m_dot_a": r"$\dot M_{\mathrm{atm,",
    r"m_dot_c": r"\dot M_{\mathrm{core,",
    r"m_dot_a_chem": r"\dot M_{\mathrm{atm,",
    r"m_dot_c_chem": r"\dot M_{\mathrm{core,",
    r"M_z": r"M_{Z,\mathrm{",
    r"regime": "regime_",
}
dict_of_acc_names = {
    f"{k}_{acc}": r"{}{}".format(v, acc) + r"}}"
    for k, v in dict_of_acc_names_template.items()
    for acc in ["peb", "gas"]
}

dict_of_extra_names = {
    "Z_vs_M": r"\frac{M_{Z}}{M}}",
    "Z_vs_M_a": r"\frac{M_{a, Z}}{M_a}}",
    "Z_vs_M_c": r"\frac{M_{c, Z}}{M_c}}",
}
NAME = {**dict_of_planet_read_names, **dict_of_acc_names, **dict_of_extra_names}


class GrowthPlot:
    def __init__(
        self,
        files,
        out,
        plotlabels=None,
        FeH=0,
        close_plots=False,
        quantities=None,
        dry_run=False,
    ):
        self.files = files
        self.out = out
        self._dry_run = dry_run

        if quantities is not None:
            self._used_quantities = set(quantities)
            self._used_quantities.update(["M", "t"])
        else:
            self._used_quantities = quantities
        self.import_data(files)
        if plotlabels is not None:
            self.plotlabels = plotlabels
        else:
            self.plotlabels = dict(zip(out, out))

        self.init_solar(FeH)
        self.composite_data()
        self.convert_units()
        self.get_len_data()
        self.close_plots = close_plots

    def init_solar(self, FeH):
        def solar_from_FeH(FeH):
            # idx = FeH_dict["FeH"].index(FeH)
            OHabu = OHabu0  # * 10 ** FeH_dict["OH"][idx]
            CHabu = CHabu0  # * 10 ** FeH_dict["CH"][idx]
            SiHabu = SiHabu0  # * 10 ** FeH_dict["SiH"][idx]
            SHabu = SHabu0  # * 10 ** FeH_dict["SH"][idx]
            FeHabu = FeHabu0  # * 10 ** FeH_dict["FeH"][idx]
            MgHabu = MgHabu0  # * 10 ** FeH_dict["MgH"][idx]
            NHabu = 1 * NHabu0
            AlHabu = 1 * AlHabu0
            TiHabu = 1 * TiHabu0
            KHabu = 1 * KHabu0
            NaHabu = 1 * NaHabu0
            VHabu = 1 * VHabu0

            return {
                "O": OHabu * MassO,
                "Si": SiHabu * MassSi,
                "C": CHabu * MassC,
                "S": SHabu * MassS,
                "Fe": FeHabu * MassFe,
                "Mg": MgHabu * MassMg,
                "H": 1,
                "He": 0.085 * MassHe,
                "Na": NaHabu * MassNa,
                "N": NHabu * MassN,
                "Al": AlHabu * MassAl,
                "Ti": TiHabu * MassTi,
                "K": KHabu * MassK,
                "V": MassV * VHabu,
            }

        self.FeH = {}
        for o in self.out:
            if type(FeH) == dict:
                """dict needs to have out as keys"""
                self.data[o]["solar"] = solar_from_FeH(FeH[o])
                self.FeH[o] = FeH[o]
            else:
                self.data[o]["solar"] = solar_from_FeH(FeH)
                self.FeH[o] = FeH

    def import_data(self, files):
        self.data = {}
        self.read_names = {}

        for i in range(len(self.out)):
            try:
                with h5py.File(files[i], "r") as f:
                    units = dict(f["/planet/units"].attrs)
                    # needed for new versions of astropy
                    units = {key: unit.decode('utf-8').replace('u.','') for key, unit in units.items()}
                    keys = list(f["/planet"].keys())
                    keys = (
                        [k for k in keys if k in self._used_quantities]
                        if self._used_quantities is not None
                        else keys
                    )
                    if "units" in keys:
                        keys.remove("units")
                break
            except OSError:
                print("WARNING: non readable file found, trying next one")

        self.faulty_list = []

        if not self._dry_run:
            for i in range(len(self.out)):
                try:
                    with h5py.File(files[i], "r") as f:
                        data = {}
                        for k in keys:
                            d = np.squeeze(np.array(f["/planet"][k]))
                            if d.dtype.type is not np.string_:
                                data[k] = u.Quantity(d, units.get(k))
                            else:
                                data[k] = d
                        data = {
                            q: data[q][np.where(~np.isnan(data["M"]))]
                            for q in data.keys()
                        }
                        self.data[self.out[i]] = data
                except OSError:
                    self.faulty_list.append(i)

            if len(self.faulty_list) > 0:
                print(
                    f"WARNING: There are in total {len(self.faulty_list)} corrupt files!"
                )
        else:
            with h5py.File(files[0], "r") as f:
                fill_data = {}
                for k in keys:
                    fill_data[k] = np.squeeze(
                        u.Quantity(np.array(f["/planet"][k]), units.get(k))
                    )

            for i in range(len(self.out)):
                data = {}
                for k in keys:
                    data[k] = fill_data[k]
                data = {q: data[q] for q in data.keys()}
                self.data[self.out[i]] = data

        self.out = [
            self.out[i] for i in range(len(self.out)) if i not in self.faulty_list
        ]
        self.dt = self.data[self.out[0]]["t"][-1] - self.data[self.out[0]]["t"][-2]

    def composite_data(self):
        def add_to_data(fct, quantity, comp_NAME, **kwargs):
            if self._used_quantities is None or quantity in self._used_quantities:
                for o in self.out:
                    self.data[o][quantity] = fct(o, **kwargs)
                    NAME[quantity] = comp_NAME

        def M_Z_vs_M_a(o):
            return (
                np.sum(self.data[o]["comp_a"][:, 0, :-7], axis=1)
                / self.data[o]["M_a"]
                * u.dimensionless_unscaled
            )

        def M_Z_vs_M(o):
            return (
                np.sum(
                    self.data[o]["comp_a"][:, 0, :-7]
                    + self.data[o]["comp_c"][:, 0, :-7],
                    axis=1,
                )
                / self.data[o]["M"]
                * u.dimensionless_unscaled
            )

        def M_Z(o):
            return np.sum(
                self.data[o]["comp_a"][:, 0, :-7] + self.data[o]["comp_c"][:, 0, :-7],
                axis=1,
            )

        def M_Z_vs_M_c(o):
            return (
                np.sum(self.data[o]["comp_c"][:, 0, :-7], axis=1)
                / self.data[o]["M_c"]
                * u.dimensionless_unscaled
            )

        def elem_a(o):
            return dict(
                zip(
                    elements,
                    (
                        self.data[o]["comp_a"][:, 0, :]
                        / self.data[o]["M_a"][:, np.newaxis]
                    ).T,
                )
            )

        def elem_c(o):
            return dict(
                zip(
                    elements,
                    (
                        self.data[o]["comp_c"][:, 0, :]
                        / self.data[o]["M_c"][:, np.newaxis]
                    ).T,
                )
            )

        def elem_tot(o):
            return dict(
                zip(
                    elements,
                    (
                        (
                            self.data[o]["comp_c"][:, 0, :]
                            + self.data[o]["comp_a"][:, 0, :]
                        )
                        / self.data[o]["M"][:, np.newaxis]
                    ).T,
                )
            )

        def mol_a(o):
            return dict(
                zip(
                    molecules,
                    (
                        self.data[o]["comp_a"][:, 1, :]
                        / self.data[o]["M_a"][:, np.newaxis]
                    ).T,
                )
            )

        def mol_c(o):
            return dict(
                zip(
                    molecules,
                    (
                        self.data[o]["comp_c"][:, 1, :]
                        / self.data[o]["M_c"][:, np.newaxis]
                    ).T,
                )
            )

        def mol_tot(o):
            return dict(
                zip(
                    molecules,
                    (
                        (
                            self.data[o]["comp_c"][:, 1, :]
                            + self.data[o]["comp_a"][:, 1, :]
                        )
                        / self.data[o]["M"][:, np.newaxis]
                    ).T,
                )
            )

        def XH(o, X, orig):
            return (
                self.data[o][f"elem_{orig}"][X]
                / self.data[o][f"elem_{orig}"]["H"]
                / self.data[o]["solar"][X]
            )

        def XY(o, X, Y, orig):
            return (
                self.data[o][f"elem_{orig}"][X]
                / self.data[o][f"elem_{orig}"][Y]
                / self.data[o]["solar"][X]
                * self.data[o]["solar"][Y]
            )

        def XH_not_normed(o, X, orig):
            return (
                self.data[o][f"elem_{orig}"][X]
                / self.data[o][f"elem_{orig}"]["H"]
                / el_masses[X]
            )

        def XY_not_normed(o, X, Y, orig):
            return (
                self.data[o][f"elem_{orig}"][X]
                / self.data[o][f"elem_{orig}"][Y]
                * el_masses[Y]
                / el_masses[X]
            )

        def pebiso(o):
            iso = np.where(
                self.data[o]["m_dot_gas"].value > 0,
                np.ones_like(self.data[o]["t"].value),
                np.zeros_like(self.data[o]["t"].value),
            )
            iso[np.argmax(iso > 0) :] = 1
            iso = iso * u.dimensionless_unscaled
            return iso

        def atmo_frac(o):
            return self.data[o]["M_a"] / self.data[o]["M"]

        def core_frac(o):
            return self.data[o]["M_c"] / self.data[o]["M"]

        add_to_data(elem_c, "elem_c", r"\mathrm{elem}_c")
        add_to_data(elem_a, "elem_a", r"\mathrm{elem}_a")
        add_to_data(elem_tot, "elem_tot", r"\mathrm{elem}")
        add_to_data(mol_c, "mol_c", r"\mathrm{mol}_c")
        add_to_data(mol_a, "mol_a", r"\mathrm{mol}_a")
        add_to_data(mol_tot, "mol_tot", r"\mathrm{mol}")

        for orig in ["a", "c", "tot"]:
            for sp in elements:
                # add orig to names
                add_to_data(XH_not_normed, f"{sp}H_{orig}", f"{sp}/H", X=sp, orig=orig)
                add_to_data(XH, f"{sp}H_normed_{orig}", f"{sp}/H", X=sp, orig=orig)
            add_to_data(XY_not_normed, f"CO_{orig}", "C/O", X="C", Y="O", orig=orig)
            add_to_data(XY_not_normed, f"SiO_{orig}", "Si/O", X="Si", Y="O", orig=orig)
            add_to_data(XY, f"CO_normed_{orig}", "C/O", X="C", Y="O", orig=orig)
            add_to_data(XY, f"SiO_normed_{orig}", "Si/O", X="Si", Y="O", orig=orig)

        add_to_data(M_Z, "M_Z", r"M_{Z}")
        add_to_data(M_Z_vs_M, "Z_vs_M", r"\frac{M_{Z}}{M}")
        add_to_data(M_Z_vs_M_a, "Z_vs_M_a", r"\frac{M_{a,Z}}{M_a}")
        add_to_data(M_Z_vs_M_c, "Z_vs_M_c", r"\frac{M_{c,Z}}{M_c}")

        add_to_data(pebiso, "pebiso", r"M_\mathrm{peb,iso}")
        add_to_data(atmo_frac, "atmo_frac", r"\frac{M_a}{M}")
        add_to_data(core_frac, "core_frac", r"\frac{M_c}{M}")
        return

    def get_len_data(self):
        max_len = np.max([len(data["t"]) for o, data in self.data.items()])
        self.len_data = max_len

    def convert_units(self):
        for o in self.out:
            for k, v in self.data[o].items():
                if hasattr(v, "unit"):
                    if v.unit == u.g:
                        self.data[o][k] = v.to(u.earthMass)
                    if v.unit == u.s:
                        self.data[o][k] = v.to(u.Myr)
                    if v.unit == u.g / u.s:
                        self.data[o][k] = v.to(u.earthMass / u.Myr)
                    if v.unit == u.cm:
                        self.data[o][k] = v.to(u.au)

    def get_label(self, n, quantity):
        if hasattr(quantity, "unit") and quantity.unit != u.dimensionless_unscaled:
            x = r"${}$ / {:latex}".format(n, quantity.unit)
        elif (
            hasattr(self.data[self.out[0]].get(n, None), "unit")
            and self.data[self.out[0]][n].unit != u.dimensionless_unscaled
        ):
            quantity = self.data[self.out[0]][n]
            x = r"${}$ / {:latex}".format(n, quantity.unit)
        elif n == "a_p":
            x = r"${}$ / {:latex}".format(n, (1 * u.au).unit)
        elif n == "t_0":
            x = r"${}$ / {:latex}".format(n, (1 * u.Myr).unit)
        elif n == "M_Z":
            x = r"${}$ / {:latex}".format(n, (1 * u.earthMass).unit)
        elif n == "CO_a":
            x = r"C/O"
        elif n == "CO_tot":
            x = r"C/O"
        else:
            x = r"${}$".format(n)
        return x

    def plot_data(
        self,
        quantities,
        out=None,
        sub_quant=None,
        t_dot=None,
        label=None,
        ulabel=True,
        fig=None,
        units=None,
        ax=None,
        plt=None,
        color_dict=None,
        use_ls=True,
        **kwargs,
    ):
        snapshot = kwargs.get("snapshot", None)
        out = self.out if out is None else out

        if t_dot is not None:
            if hasattr(t_dot, "__iter__"):
                dotspace = t_dot.get("dotspace", DOT_SPACE)
                dotspace_large = t_dot.get("dotspace_large", DOT_SPACE_LARGE)
                use_t_dot = t_dot.get("use_t_dot", True)
                color_t_dot = t_dot.get("color", ["gray"])
                color_t_dot_small = t_dot.get("color_small", ["gray"])
            else:
                dotspace = 1 * DOT_SPACE
                dotspace_large = 1 * DOT_SPACE_LARGE
                use_t_dot = True
                color_t_dot_small = ["gray"]
                color_t_dot = ["gray"]
        else:
            use_t_dot = False

        if snapshot is None:

            def get_dots(x_data, y_data, t):
                t_at_dot = np.arange(
                    dotspace.cgs.value, np.max(t[:-1]).cgs.value, dotspace.cgs.value
                )
                x = np.interp(t_at_dot, t.cgs.value, x_data)
                y = np.interp(t_at_dot, t.cgs.value, y_data)

                t_at_dot = np.arange(
                    dotspace_large.cgs.value,
                    np.max(t[:-1]).cgs.value,
                    dotspace_large.cgs.value,
                )
                x_large = np.interp(t_at_dot, t.cgs.value, x_data)
                y_large = np.interp(t_at_dot, t.cgs.value, y_data)

                x, y = u.Quantity(
                    [xi for xi in x if xi not in x_large], x.unit
                ), u.Quantity([yi for yi in y if yi not in y_large], y.unit)
                return x, y, x_large, y_large

            if not hasattr(self, "_plot_x_data"):
                self._plot_x_data = {}
            if not hasattr(self, "_plot_y_data"):
                self._plot_y_data = {}
            if not hasattr(self, "_plot"):
                self._plot = {}

            self._plot_x_data[ax] = {}
            self._plot_y_data[ax] = {}
            self._plot[ax] = {}
            for o in out:
                if color_dict is not None:
                    kwargs["color"] = color_dict.get(o, None)

                planet = self.data[o].copy()
                if sub_quant is not None:
                    for q in quantities:
                        if type(planet[q]) == dict:
                            planet[q] = planet[q][sub_quant]

                if label is None:
                    legend = self.plotlabels[o]
                else:
                    legend = label

                if units is not None:
                    self._plot_x_data[ax][o], self._plot_y_data[ax][o] = planet[
                        quantities[0]
                    ].to(units[0]), planet[quantities[1]].to(units[1])
                else:
                    self._plot_x_data[ax][o], self._plot_y_data[ax][o] = (
                        planet[quantities[0]],
                        planet[quantities[1]],
                    )

                if ax is not None:
                    if use_ls:
                        ls_mask = self.data[o]["pebiso"] == 0
                        not_ls_mask = np.logical_not(ls_mask)
                        not_ls_mask[np.argmin(ls_mask) - 1] = True
                        (self._plot[ax][o],) = ax.loglog(
                            self._plot_x_data[ax][o][not_ls_mask],
                            self._plot_y_data[ax][o][not_ls_mask],
                            label=legend,
                            ls=":",
                            **kwargs,
                        )
                    else:
                        ls_mask = np.arange(len(self._plot_x_data[ax][o]))

                    (self._plot[ax][o],) = ax.loglog(
                        self._plot_x_data[ax][o][ls_mask],
                        self._plot_y_data[ax][o][ls_mask],
                        label=legend,
                        **kwargs,
                    )

                    ax.set_xlabel(
                        self.get_label(
                            NAME.get(quantities[0], quantities[0]),
                            self._plot_x_data[ax][o],
                        )
                    )
                    ax.set_ylabel(
                        self.get_label(
                            NAME.get(quantities[1], quantities[1]),
                            self._plot_y_data[ax][o],
                        )
                    )

                    if use_t_dot:
                        x, y, x_large, y_large = get_dots(
                            planet[quantities[0]], planet[quantities[1]], planet["t"]
                        )
                        if len(x) > 0:
                            ax.scatter(
                                x,
                                y,
                                marker="o",
                                s=10,
                                color=color_t_dot_small[: len(x)],
                                alpha=1,
                            )
                        if len(x_large) > 0:
                            ax.scatter(
                                x_large,
                                y_large,
                                s=25,
                                marker="o",
                                color=color_t_dot[: len(x_large)],
                                alpha=1,
                            )

                elif plt is not None:
                    (self._plot[o],) = plt.loglog(
                        self._plot_x_data[ax][o],
                        self._plot_y_data[ax][o],
                        label=legend,
                        **kwargs,
                    )

                    plt.xlabel(
                        self.get_label(
                            NAME.get(quantities[0], quantities[0]),
                            self._plot_x_data[ax][o],
                        )
                    )
                    plt.ylabel(
                        self.get_label(
                            NAME.get(quantities[1], quantities[1]),
                            self._plot_y_data[ax][o],
                        )
                    )
                    plt.title(
                        self.get_label(
                            NAME.get(quantities[1], quantities[1]),
                            NAME.get(quantities[0], quantities[0]),
                        )
                    )

                    if use_t_dot:
                        x, y, x_large, y_large = get_dots(
                            planet[quantities[0]], planet[quantities[1]], planet["t"]
                        )
                        plt.plot(x, y, "ko", markersize=4, color="gray", alpha=1)
                        plt.plot(
                            x_large, y_large, "ko", markersize=8, color="gray", alpha=1
                        )

        else:
            for o in out:
                if o in self.out:
                    snapshot = np.minimum(snapshot, len(self._plot_x_data[ax][o]))
                    self._plot[ax][o].set_data(
                        self._plot_x_data[ax][o][:snapshot],
                        self._plot_y_data[ax][o][:snapshot],
                    )

        return self._plot[ax]

    def apply_mask(self, mask):
        """
        filters unused data

        e.g.
        mask = [o for o in Plot.out if Plot.data[o]["a_p"][-1] > 8]
        Plot.apply_mask(mask)

        Parameters
        ----------
        mask: a subset of out

        Returns
        -------


        """
        self.data = {k: v for (k, v) in self.data.items() if k in mask}
        self.out = mask
        self.plotlabels = {k: v for (k, v) in self.plotlabels.items() if k in mask}

    def pqplot_video(
        self,
        plt,
        quantity,
        frames,
        lines=[],
        sub_qaunt=None,
        log_data=False,
        units=None,
        **kwargs,
    ):
        """not working anymore, needs some more fixes"""
        fig, ax = plt.subplots(sharex=True, sharey=True)
        div = make_axes_locatable(ax)
        cax = div.append_axes("right", "5%", "5%")

        def animate(i):
            [ax.clear() for ax in fig.axes]
            self.pqplot(
                plt=plt,
                snapshot=i,
                fig=fig,
                quantity=quantity,
                lines=lines,
                sub_qaunt=sub_qaunt,
                log_data=log_data,
                cax=cax,
                units=units,
                **kwargs,
            )

        return animation.FuncAnimation(fig, animate, frames=frames)

    def pqplot(
        self,
        fig,
        quantity,
        lines=[],
        cax=None,
        sub_quant=None,
        log_data=False,
        symlog_data=False,
        snapshot=-1,
        units=None,
        icelines=None,
        ax=None,
        shape=None,
        cbar_location=None,
        **kwargs,
    ):
        (
            data_all,
            extend_r,
            p_data,
            p_name,
            q_data,
            q_name,
            r_list,
            r_name,
            x_data_all,
            y_data_all,
            z_data,
            z_data_all,
        ) = self.extract_r0t0_data(quantity, snapshot, sub_quant, units)

        if shape is None:
            if len(extend_r) > 2:
                shape = (int(np.ceil(len(extend_r) / 2)), 2)
            else:
                shape = (len(extend_r), 1)

        # begin the plotting
        if ax is None:
            ax = fig.subplots(*shape, sharex=True, sharey=True)
            ax = np.array(ax)

        if log_data:
            if np.array(z_data_all).min() <= 0:
                kwargs["norm"] = colors.SymLogNorm(linthresh=0.03)
            else:
                kwargs["norm"] = colors.LogNorm()
        elif symlog_data:
            kwargs["norm"] = colors.SymLogNorm(linthresh=0.03)
        if kwargs.get("vmin", None) is None:
            kwargs["vmin"] = np.minimum(np.array(z_data_all).min(), 0.03)
        if kwargs.get("vmax", None) is None:
            kwargs["vmax"] = np.maximum(np.array(z_data_all).max(), 0.04)

        for i in range(len(extend_r)):
            if hasattr(ax, "flat"):
                ax_temp = ax.flat[i]
            else:
                ax_temp = ax
            im = ax_temp.pcolor(x_data_all[i], y_data_all[i], z_data_all[i], **kwargs)
            ax_temp.set_xlabel(self.get_label(p_name, x_data_all[i]))
            ax_temp.set_ylabel(self.get_label(q_name, y_data_all[i]))

            for l in lines:
                if isinstance(l[0], str):
                    Z = np.array(
                        [
                            d[l[0]].value[snapshot if snapshot < len(d["t"]) else -1]
                            for d in data_all[i].values()
                        ]
                    )
                else:
                    Z = l[0]

                zi, yi, xi = np.histogram2d(
                    p_data,
                    q_data,
                    bins=(len(np.unique(p_data)), len(np.unique(q_data))),
                    weights=Z,
                    normed=False,
                )
                # zi = gaussian_filter(zi.T, 0.15)
                x_step = l[1].pop("x_step", 1)
                y_step = l[1].pop("y_step", 1)
                zi = zi.T
                zi = np.ma.masked_equal(zi, np.inf)
                CS = ax_temp.contour(
                    x_data_all[i][::x_step, ::y_step],
                    y_data_all[i][::x_step, ::y_step],
                    zi[::x_step, ::y_step],
                    **l[1],
                )
                ax_temp.clabel(CS, l[2])
                l[1]["x_step"] = x_step
                l[1]["y_step"] = y_step

            if icelines is not None:
                self.get_iceline_background(icelines, ax_temp)

            time = 3 * u.Myr if snapshot == -1 else snapshot * self.dt
            if r_name is not None:
                if r_name == "case":
                    ax_temp.set_title(r"{}".format(r_list[i]))
                else:
                    ax_temp.set_title(r"{}={}".format(r_name, r_list[i]))

        shrink = 0.4 if cbar_location == "bottom" else 0.6
        aspect = 6 if cbar_location == "bottom" else 20
        if cax is not None:
            cax.cla()
        if len(extend_r) != shape[0] * shape[1]:
            cbar = fig.colorbar(im, ax=ax[:-1], cax=cax, shrink=shrink, aspect=aspect)
            ax.flat[-1].remove()
        else:
            cbar = fig.colorbar(
                im, ax=ax, cax=cax, location=cbar_location, shrink=shrink, aspect=aspect
            )

        try:
            for temp_ax in ax.flat:
                temp_ax.label_outer()
        except Exception as e:
            print(f"no label outer possible {e}")

        cbar.ax.set_title(self.get_label(quantity, z_data), y=1.03)

        fig.suptitle(
            r"${}$-${}$ Plot for {}, t = {}".format(
                p_name, q_name, quantity, time.to(u.Myr)
            )
        )
        return im

    def extract_r0t0_data(self, quantity, snapshot, sub_quant, units):
        """data in pipline should look like: r0 = [20, 10, 5, 20, 10, 5]), t0 = [1e5, 1e5, 1e5, 1e6, 1e6, 1e6]) - in other words: flattened meshgrids"""
        format_string = "{}:{}, {}:{}"
        format_string_long = "{}:{}, {}:{}, {}:{}"
        if parse.parse(format_string_long, self.plotlabels[self.out[0]]) is not None:
            # -> case of three (or possible more) parameters par_1, par_2, par_3
            pa = [parse.parse(format_string, self.plotlabels[o])[1] for o in self.out]
            try:
                pa = [float(p) for p in pa]
            except:
                pass
            unique = np.unique(pa)
            if np.all(unique == np.array(["evaporation", "plain"])):
                unique = unique[::-1]

            extend_r = [
                [self.out[i] for i in range(len(self.out)) if pa[i] == lsta]
                for lsta in unique
            ]
        elif parse.parse(format_string, self.plotlabels[self.out[0]]) is not None:
            # -> case of only par_1, par_2
            extend_r = [self.out]
        else:
            print("data is not usable")
            quit()
        z_data_all, x_data_all, y_data_all, data_all, r_list = [], [], [], [], []
        for sub in extend_r:
            data = {key: self.data[key] for key in sub}
            p_data, q_data = [], []
            for o in sub:
                if len(extend_r) > 1:
                    r_name, r, p_name, p, q_name, q = parse.parse(
                        format_string_long, self.plotlabels[o]
                    )
                else:
                    p_name, p, q_name, q = parse.parse(
                        format_string, self.plotlabels[o]
                    )
                    r_name = None
                    r = 0
                p_data.append(float(p))
                q_data.append(float(q))

            p_data = np.array(p_data)
            q_data = np.array(q_data)
            if units is not None:
                p_data = (p_data * units[0]).cgs.to(units[0])
                q_data = (q_data * units[1]).cgs.to(units[1])

            if sub_quant is not None:
                z_data = np.array(
                    [
                        data[quantity][sub_quant][
                            snapshot if snapshot < len(data["t"]) else -1
                        ]
                        for data in data.values()
                    ]
                )

            else:
                z_data = np.array(
                    [
                        data[quantity].value[
                            snapshot if snapshot < len(data["t"]) else -1
                        ]
                        for data in data.values()
                    ]
                )

            if units is not None:
                z_data = z_data.to(units[2])

            zi, yi, xi = np.histogram2d(
                p_data,
                q_data,
                bins=(len(np.unique(p_data)), len(np.unique(q_data))),
                weights=z_data,
                normed=False,
            )
            zi = np.ma.masked_equal(zi, 0)
            zi = np.ma.masked_equal(zi, np.inf)
            X, Y = np.meshgrid(np.unique(p_data), np.unique(q_data))
            z_data_all.append(zi.T)
            x_data_all.append(X)
            y_data_all.append(Y)
            r_list.append(r)
            data_all.append(data)
        x_data_all = np.array(x_data_all)
        y_data_all = np.array(y_data_all)
        return (
            data_all,
            extend_r,
            p_data,
            p_name,
            q_data,
            q_name,
            r_list,
            r_name,
            x_data_all,
            y_data_all,
            z_data,
            z_data_all,
        )

    def get_iceline_background(self, file, ax):
        def get_data_iceline():
            data = {}
            with tables.open_file(file, mode="r") as f:
                data["T"] = np.array(f.root.disk.T)
                data["r"] = np.squeeze(np.array(f.root.disk.r) / const.au.cgs.value)
                data["t"] = (
                    np.squeeze(np.array(f.root.disk.t)) / 1e6 / 365.25 / 24 / 3600
                )

            data["icelines"] = np.array(
                [20, 30, 70, 150, 371, 970, 1354, 1423, 1500, 1653, 2000]
            )
            data["iceline_pos"] = np.array(
                [
                    np.searchsorted(-data["T"][i], -data["icelines"])
                    for i in range(len(data["T"]))
                ]
            )
            data["iceline_r"] = np.array(
                [data["r"][data["iceline_pos"][i]] for i in range(len(data["t"]))]
            )
            data["iceline_name"] = [
                "CO & N2",
                "CH4",
                "CO2",
                "H2O & H2S",
                "Fe3O4",
                "NaAlSi3O8 & KAlSi3O8",
                "Mg2SiO4 & Fe2O3",
                "VO",
                "MgSiO3",
                "Al2O3",
                "TiO",
            ]
            return data

        data = get_data_iceline()
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        texts = []
        for i in range(len(data["icelines"])):
            if np.any(
                np.logical_and(
                    data["iceline_r"][:, i] > xlim[0], data["iceline_r"][:, i] < xlim[1]
                )
            ):
                ymask = np.logical_and(data["t"] > ylim[0], data["t"] < ylim[1])
                ax.plot(
                    data["iceline_r"][ymask, i],
                    data["t"][ymask],
                    ls="--",
                    color="gray",
                    alpha=0.8,
                )
                t = ax.text(
                    np.mean(data["iceline_r"][ymask, i]),
                    0.5 * (np.mean(ax.get_ylim())),
                    data["iceline_name"][i],
                    fontsize=10,
                    rotation=90,
                    ha="center",
                    va="center",
                    color="white",
                    backgroundcolor="gray",
                )

                t.set_bbox(dict(alpha=0.0))
                texts.append(t)
        if bool(texts):
            adjust_text(texts, ax=ax)

    def finish_plot(self, plt, file=None, no_legend=False):
        if not no_legend:
            plt.legend()
        if file is not None:
            plt.savefig(file)
        if self.close_plots:
            plt.close()
        else:
            plt.show()
