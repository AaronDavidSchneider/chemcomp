from typing import Tuple
from chemcomp.helpers.units import *

T_array = np.array(
    [
        0,  # background gas (H, He)
        COtemp,
        N2temp,
        CH4temp,
        CO2temp,
        NH3temp,
        H2Stemp,
        H2Otemp,
        Fe3O4temp,
        Ctemp,
        FeStemp,
        NaAlSi3O8temp,
        KAlSi3O8temp,
        Mg2SiO4temp,
        Fe2O3temp,
        VOtemp,
        MgSiO3temp,
        Al2O3temp,
        TiOtemp,
    ]
)
T_array = T_array[:, np.newaxis]

elmasses = np.array(
    [
        MassC,
        MassO,
        MassFe,
        MassS,
        MassMg,
        MassSi,
        MassNa,
        MassK,
        MassN,
        MassAl,
        MassTi,
        MassV,
        MassH,
        MassHe,
        1,  # dummy not used (to match dimensions)
        1,  # dummy not used (to match dimensions)
        1,  # dummy not used (to match dimensions)
        1,  # dummy not used (to match dimensions)
        1,  # dummy not used (to match dimensions)
    ],
)

FeH_dict = {
    "FeH": [-0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4],
    "SiH": [-0.32, -0.25, -0.17, -0.09, 0.0, 0.1, 0.2, 0.32, 0.41],
    "SH": [-0.32, -0.25, -0.17, -0.09, 0.0, 0.1, 0.2, 0.32, 0.41],
    "MgH": [-0.26, -0.19, -0.1, -0.01, 0.08, 0.19, 0.3, 0.42, 0.53],
    "CH": [-0.29, -0.22, -0.14, -0.06, 0.03, 0.11, 0.2, 0.28, 0.36],
    "OH": [-0.19, -0.13, -0.08, -0.02, 0.03, 0.08, 0.13, 0.19, 0.23],
}


def split_molecules(abundances_inp):
    """

    Parameters
    ----------
    abundances_inp: 9xN numpy array of either the gas or solid phase. Can be generated with Chemistry.abundances

    Returns
    -------
    13xN numpy array: For the twelve heavy elements, return the elemental abundances per grid cell, starting from molecular abundances
    The 13th entry gives the abundance of hydrogen contained in the supplied molecular abundances

    """
    (
            rest_mol,
            TotMassCO,
            TotMassN2,
            TotMassCH4,
            TotMassCO2,
            TotMassNH3,
            TotMassH2S,
            TotMassH2O,
            TotMassFe3O4,
            TotMassC,
            TotMassFeS,
            TotMassNaAlSi3O8,
            TotMassKAlSi3O8,
            TotmassMg2SiO4,
            TotMassFe2O3,
            TotMassVO,
            TotMassMgSiO3,
            TotMassAl2O3,
            TotMassTiO,
        ) = abundances_inp

    TotC = (
            TotMassCO * MassC / MassCO
            + TotMassCH4 * MassC / MassCH4
            + TotMassCO2 * MassC / MassCO2
            + TotMassC * MassC / MassC
    )
    TotO = (
            TotMassCO * MassO / MassCO
            + TotMassCO2 * MassO * 2.0 / MassCO2
            + TotMassH2O * MassO / MassH2O
            + TotMassFe3O4 * 4.0 * MassO / MassFe3O4
            + TotmassMg2SiO4 * 4.0 * MassO / MassMg2SiO4
            + TotMassFe2O3 * 3 * MassO / MassFe2O3
            + TotMassMgSiO3 * MassO * 3.0 / MassMgSiO3
            + TotMassTiO * MassO / MassTiO
            + TotMassAl2O3 * 3 * MassO / MassAl2O3
            + TotMassKAlSi3O8 * 8 * MassO / MassKAlSi3O8
            + TotMassNaAlSi3O8 * 8 * MassO / MassNaAlSi3O8
            + TotMassVO * MassO / MassVO
    )
    TotFe = (
            TotMassFe3O4 * 3.0 * MassFe / MassFe3O4
            + TotMassFe2O3 * 2.0 * MassFe / MassFe2O3
            + TotMassFeS * MassFe / MassFeS
    )

    TotS = TotMassFeS * MassS / MassFeS + TotMassH2S * MassS / MassH2S
    TotMg = (
            TotmassMg2SiO4 * MassMg * 2.0 / MassMg2SiO4
            + TotMassMgSiO3 * MassMg / MassMgSiO3
    )
    TotSi = (
            TotmassMg2SiO4 * MassSi / MassMg2SiO4
            + TotMassMgSiO3 * MassSi / MassMgSiO3
            + TotMassKAlSi3O8 * 3 * MassSi / MassKAlSi3O8
            + TotMassNaAlSi3O8 * 3 * MassSi / MassNaAlSi3O8

    )

    TotK = TotMassKAlSi3O8 * MassK / MassKAlSi3O8
    TotNa = TotMassNaAlSi3O8 * MassNa / MassNaAlSi3O8
    TotAl = (
            TotMassKAlSi3O8 * MassAl / MassKAlSi3O8
            + TotMassNaAlSi3O8 * MassAl / MassNaAlSi3O8
            + TotMassAl2O3 * 2 * MassAl / MassAl2O3
    )
    TotTi = TotMassTiO * MassTi / MassTiO
    TotV = TotMassVO * MassV / MassVO
    TotN = TotMassNH3 * MassN / MassNH3 + TotMassN2 * 2 * MassN / MassN2

    TotH_mol = (
            TotMassCH4 * 4.0 * MassH / MassCH4
            + TotMassH2O * 2.0 * MassH / MassH2O
            + TotMassNH3 * 3 * MassH / MassNH3
            + TotMassH2S * 2 * MassH / MassH2S
    )
    # Note how the last entry of the array is TotH_mol
    return np.array([TotC, TotO, TotFe, TotS, TotMg, TotSi, TotNa, TotK, TotN, TotAl, TotTi, TotV, TotH_mol])

def calc_composition(abundances_inp: np.array, solid: bool, He_Mbg: float) -> Tuple[np.array, np.array]:
    """

    Parameters
    ----------
    abundances: 9xN numpy array that can be generated with Chemistry.abundances
    solid: bool that tells if the input array has solid compositions or gas compositions
    He_Mbg: Total mass of He with respect to the total mass of the background H+He gas

    Returns
    -------
    Nx2x10 numpy array: for every gridpoint (T), you get 2 arrays - elemental and molecular abundances
    Nd numpy array: the total mass of Hydrogen in molecules per grid cell

    """
    
    molecule_elements = split_molecules(abundances_inp)
    rest_mol = abundances_inp[0]
    rest = np.sum(molecule_elements[:-1], axis=0)
    TotH_mol = molecule_elements[-1]

    if solid: # These are the solid components, getting the H, He abundance is trivial
        TotH = TotH_mol
        TotHe = np.zeros_like(rest)
    else: # The gas components, where information about the Helium content is required and therefore the supplied He_Mbg argument is used
        TotHe = He_Mbg * rest_mol
        TotH = (1 - rest - TotHe)
        
    # Test about the He+H split
    assert np.isclose(TotH + TotHe - TotH_mol, rest_mol).all(), "Background gas after split does not equal the total background gas from the molecular abundances"

    elements = np.array(
        [
            *molecule_elements[:-1],
            TotH,
            TotHe,
            np.zeros_like(rest),  # dummy not used (to match dimensions)
            np.zeros_like(rest),  # dummy not used (to match dimensions)
            np.zeros_like(rest),  # dummy not used (to match dimensions)
            np.zeros_like(rest),  # dummy not used (to match dimensions)
            np.zeros_like(rest),  # dummy not used (to match dimensions)
        ],
    )

    molecules = abundances_inp.copy() # Avoid unwanted side-effects when returning the supplied molecular abundances

    return np.transpose(np.stack([elements, molecules]), (2, 0, 1)), TotH_mol


class Chemistry:
    """ Chemistry class. Please be careful when changing or adding molecular species,
    some occurancies of the chemistry rely (hardcoded) on the order of the species!
    """
    def __init__(self, config):
        #   Adundances in terms of H (Asplund et al. 2009) + variation from Chiara plot
        if not config.get("use_FeH", True):
            self.OHabu = OHabu0 * 10 ** float(config.get("OH", 0.0))
            self.CHabu = CHabu0 * 10 ** float(config.get("CH", 0.0))
            self.SiHabu = SiHabu0 * 10 ** float(config.get("SiH", 0.0))
            self.SHabu = SHabu0 * 10 ** float(config.get("SH", 0.0))
            self.FeHabu = FeHabu0 * 10 ** float(config.get("FeH", 0.0))
            self.MgHabu = MgHabu0 * 10 ** float(config.get("MgH", 0.0))
        else:
            idx = FeH_dict["FeH"].index(float(config.get("FeH", 0.0)))
            self.OHabu = OHabu0 * 10 ** FeH_dict["OH"][idx]
            self.CHabu = CHabu0 * 10 ** FeH_dict["CH"][idx]
            self.SiHabu = SiHabu0 * 10 ** FeH_dict["SiH"][idx]
            self.SHabu = SHabu0 * 10 ** FeH_dict["SH"][idx]
            self.FeHabu = FeHabu0 * 10 ** FeH_dict["FeH"][idx]
            self.MgHabu = MgHabu0 * 10 ** FeH_dict["MgH"][idx]

        self.NHabu = NHabu0 * 10 ** float(config.get("NH", 0.0))
        self.AlHabu = AlHabu0 * 10 ** float(config.get("AlH", 0.0))
        self.TiHabu = TiHabu0 * 10 ** float(config.get("TiH", 0.0))
        self.KHabu = KHabu0 * 10 ** float(config.get("KH", 0.0))
        self.NaHabu = NaHabu0 * 10 ** float(config.get("NaH", 0.0))
        self.VHabu = VHabu0 * 10 ** float(config.get("VH", 0.0))

        self.HeHabu = HeHabu0 * 10 ** float(config.get("HeH", 0.0))
        self.f_H2 = 1 / (2 * self.HeHabu + 1)

        self.elabu_init = np.array(
            [
                self.CHabu,
                self.OHabu,
                self.FeHabu,
                self.SHabu,
                self.MgHabu,
                self.SiHabu,
                self.NaHabu,
                self.KHabu,
                self.NHabu,
                self.AlHabu,
                self.TiHabu,
                self.VHabu,
                0.0,
                0.0,
                0.0,  # dummy not used (to match dimensions)
                0,  # dummy not used (to match dimensions)
                0,  # dummy not used (to match dimensions)
                0,  # dummy not used (to match dimensions)
                0,  # dummy not used (to match dimensions)
            ],
        )
        # compute temperature independent abu
        self.abuCO2 = float(config.get("CO2_frac", 0.1)) * self.CHabu
        self.abuFe2O3 = 0.25 * (self.FeHabu - 0.9 * self.SHabu)
        self.abuFe3O4 = (1.0 / 6.0) * (self.FeHabu - 0.9 * self.SHabu)
        self.abuFeS = 0.9 * self.SHabu
        self.abuMgSiO3 = self.MgHabu - 2 * (self.MgHabu - (self.SiHabu - 3 * self.KHabu - 3 * self.NaHabu))
        self.abuMg2SiO4 = self.MgHabu - (self.SiHabu - 3 * self.KHabu - 3 * self.NaHabu)
        self.abuCO = float(config.get("CO_frac", 0.45)) * self.CHabu

        if config.get("CH4_frac") is not None:
            self.abuCH4 = config.get("CH4_frac") * self.CHabu
        else:
            self.abuCH4 = (0.45 - float(config.get("C_frac", 0.2))) * self.CHabu

        self.abuC = float(config.get("C_frac", 0.2)) * self.CHabu

        self.abuH2S = 0.1 * self.SHabu
        self.abuNH3 = 0.1 * self.NHabu
        self.abuN2 = 0.45 * self.NHabu
        self.abuNaAlSi3O8 = 1 * self.NaHabu
        self.abuKAlSi3O8 = 1 * self.KHabu
        self.abuTiO = 1 * self.TiHabu
        self.abuVO = 1 * self.VHabu
        self.abuAl2O3 = 0.5 * (self.AlHabu - (self.KHabu + self.NaHabu))

        self.abuH2O = self.OHabu - (
                3 * self.abuMgSiO3
                + 4 * self.abuMg2SiO4
                + self.abuCO
                + 2 * self.abuCO2
                + 3 * self.abuFe2O3
                + 4 * self.abuFe3O4
                + self.abuTiO
                + 3 * self.abuAl2O3
                + 8 * self.abuNaAlSi3O8
                + 8 * self.abuKAlSi3O8
                + self.abuVO
        )

        assert np.isclose(self.abuC + self.abuCO + self.abuCO2 + self.abuCH4, self.CHabu), \
            "Broken chem model, we loose or gain C"

        self.abu_array = np.array(
            [
                0,
                self.abuCO,
                self.abuN2,
                self.abuCH4,
                self.abuCO2,
                self.abuNH3,
                self.abuH2S,
                self.abuH2O,
                self.abuFe3O4,
                self.abuC,
                self.abuFeS,
                self.abuNaAlSi3O8,
                self.abuKAlSi3O8,
                self.abuMg2SiO4,
                self.abuFe2O3,
                self.abuVO,
                self.abuMgSiO3,
                self.abuAl2O3,
                self.abuTiO,
            ]
        )

        M_array = np.array(
            [
                2.3,
                MassCO,
                MassN2,
                MassCH4,
                MassCO2,
                MassNH3,
                MassH2S,
                MassH2O,
                MassFe3O4,
                MassC,
                MassFeS,
                MassNaAlSi3O8,
                MassKAlSi3O8,
                MassMg2SiO4,
                MassFe2O3,
                MassVO,
                MassMgSiO3,
                MassAl2O3,
                MassTiO,
            ]
        )

        self.M_array = M_array[:, np.newaxis]

        self.gamma = 1.4
        # important for outputing:
        self.mask_array = np.ones((2, 19), dtype=int)
        self.mask_array[0, :-5] = 0  # set the unactive dimensions of chemistry model
        # there is one molecules type more then there are elements

        self.mu = 2.3

        return

    def get_position_of_ice(self, T: np.array) -> np.array:
        """
        function that can be used to determine the indicees of the icelines
        This index is the index of the first cell in r/T that has no gas in sigma

        gas | gas | gas | solid | solid | solid
                           idx

        Parameters
        ----------
        T: temperature N-dimensional array in cgs

        Returns
        -------
        idx: position of icelines

        """
        # exclude phantom iceline of rest_mol
        idx = np.squeeze(np.searchsorted(-T, -T_array, side="right"))
        # np.squeeze uses 1e-6 s. Total: 5.3e-6 s
        idx = np.minimum(idx, T.size - 2)

        return idx
    
    def calc_He_Mbg(self, abundances_gas: np.array, abundances_solid: np.array, nHe_nH: np.array, T: np.array) -> float:
        """
        Calculates the mass ratio of Helium with respect to the background gas.
        
        Parameters
        ----------
        abundances_gas: 9xN numpy array containing gas composition. Can be generated with Chemistry.abundances
        abundances_solid: 9xN numpy array containing solid composition. Can be generated with Chemistry.abundances
        nHe_nH: nHe/nH number density ratio of Helium and Hydrogen. This should match with your initial condition
        T: Temperature array
    
        Returns
        -------
        A float describing the mass ratio of Helium w.r.t. the background gas
    
        """
        TotH_solid = split_molecules(abundances_solid)[-1]
        molecule_elements = split_molecules(abundances_gas)
        TotH_mol_gas = molecule_elements[-1]
        rest = np.sum(molecule_elements[:-1], axis=0)
        rest_mol = abundances_gas[0]
        Z = self.get_solid_heavies(T)
        Mdust_Mgas = Z*rest_mol # = Mmol_solid / Mbg * Mbg / Mgas = Mmol_solid/Mgas = Mdust/Mgas
        TotH_solid_Mgas = TotH_solid*Mdust_Mgas # Total mass of Hydrogen in solids over the total gas mass
        
        # This gives the total hydrogen mass in both gaseous and solid phase, normed to the gas mass, (Hgas+Hsolid)/gas
        TotH = (1 - rest + TotH_solid_Mgas) / (1 + nHe_nH * MassHe)
        # This gives Hgas/gas
        TotH_gas = TotH - TotH_solid_Mgas
        # This gives H/gas * He/H = He/gas
        TotHe = TotH * nHe_nH * MassHe
        
        # Some tests to check for some errors for the calculation of TotH and TotHe
        assert np.isclose(1-rest-molecule_elements[-1], rest_mol).all(), "Background gas from the molecular abundance array does not match background gas after splitting the molecules"
        assert np.isclose(TotH_gas + TotHe - TotH_mol_gas, rest_mol).all(), "Splitting the background gas into hydrogen and helium did not work"
        assert np.isclose(TotHe / (TotH*MassHe), nHe_nH).all(), "After separating solid and gaseous H, [He/H] is no longer correct"
        assert np.isclose(TotH_gas+TotHe+rest, 1).all(), "Adding up all gaseous components does not result in the total gas mass"

        
        # It is expected that He_Mbg is constant in space, so it's just saved as a float
        assert np.isclose((TotHe / rest_mol).std(), 0), "He_Mbg is not constant."
        He_Mbg = (TotHe / rest_mol)[0] 
        return He_Mbg

    def get_solid_composition(self, sigma_dust_mol: np.array, T: np.array) -> np.array:
        sigma_dust = np.sum(sigma_dust_mol, axis=1)
        abundances_solid = np.transpose(
            np.divide(sigma_dust_mol, sigma_dust[:, np.newaxis], out=np.zeros_like(sigma_dust_mol),
                      where=np.logical_and(sigma_dust[:, np.newaxis] > 1.0e-60, sigma_dust_mol > 1.0e-60)))
        chem_solid, _ = calc_composition(abundances_solid, solid=True, He_Mbg=self.He_Mbg)
        if np.max(T) > np.max(T_array):
            chem_solid[:np.argmax(T) + 1] = 0

        return chem_solid

    def get_gas_composition(self, sigma_g_components_mol: np.array = None) -> np.array:
        """
        wrapperfunction for use with the chemistry calculations. Also calculates the metalicity.

        Parameters
        ----------
        sigma_g_components: sigma_components in cgs

        Returns
        -------
        chem_solid: joint Array of solid molecular and atomic abundancies
        chem_gas: joint Array of gas molecular and atomic abundancies

        Dimensions of chem_gas and chem_solid:
        N x 2 x 19
        unused elements: see self.mask_array

        """

        sigma = np.sum(sigma_g_components_mol, axis=1)
        abundances_gas = np.transpose(
            np.divide(sigma_g_components_mol, sigma[:, np.newaxis], out=np.zeros_like(sigma_g_components_mol),
                      where=sigma[:, np.newaxis] != 0))
        chem_gas, _ = calc_composition(abundances_gas, solid=False, He_Mbg=self.He_Mbg)
        self.mu = sigma / np.sum(sigma_g_components_mol.T / self.M_array, axis=0)

        return chem_gas

    def set_z(self, Z0):
        """ set the fraction of heavy elements in the disk """
        # The fraction of heavy !molecules! in either gas or solids with respect to the backgroundgas (H+He):
        self.dtg0 = Z0

        # ratio of heavies from the sum over all molecular species:
        self._ZH = np.sum(self.abu_array[1:] * np.squeeze(self.M_array)[1:])
        rezHabu = MassH + self.HeHabu*MassHe + (self.elabu_init*elmasses).sum()
        self._dtg0_of_chem_model = self._ZH/rezHabu

        return

    def get_solid_heavies(self, T):
        """ returns a modified dust to gas ratio according to the temperature of the disk """
        abu_array = np.ones_like(T) * self.abu_array[:, np.newaxis]
        M_masked_solid = np.where(T < T_array, np.ones_like(T) * self.M_array, np.zeros_like(abu_array))
        M_masked_solid[0, :] = 0  # exclude restmol
        Tot_M_solid = np.sum(M_masked_solid * abu_array, axis=0) + 1e-300

        return self.dtg0 * Tot_M_solid / self._ZH

    def get_composition(self, T: np.array) -> np.array:
        """
        wrapperfunction for use with the chemistry calculations. Also calculates the metalicity.
        If the object does not already have the "He_Mbg" attribute, it is calculated

        Parameters
        ----------
        T: Temperature

        Returns
        -------
        chem_solid: joint Array of solid molecular and atomic abundancies
        chem_gas: joint Array of gas molecular and atomic abundancies

        Dimensions of chem_gas and chem_solid:
        N x 2 x 19
        unused elements: see self.mask_array

        """
        abundances_solid, abundances_gas = self.abundances(T)
        self.mu = np.sum(abundances_gas, axis=0) / np.sum(abundances_gas / self.M_array, axis=0)

        # massratio tests for the elemental composition:
        assert np.all(np.isclose(np.sum(abundances_gas, axis=0),
                                 np.ones_like(T))), "error in chemistry, total sum of gas abundancies is not 1"
        assert np.all(np.isclose(np.sum(abundances_solid, axis=0),
                                 np.sum(abundances_solid, axis=0) != 0, np.ones_like(T),
                                 np.zeros_like(T))), "error in chemistry, total sum of solid abundancies is not 1"

        self.He_Mbg = self.calc_He_Mbg(abundances_gas, abundances_solid, self.HeHabu, T)
        chem_solid, _ = calc_composition(abundances_solid, solid=True, He_Mbg=self.He_Mbg)
        chem_gas, _ = calc_composition(abundances_gas, solid=False, He_Mbg=self.He_Mbg)

        if np.max(T) > np.max(T_array):
            chem_solid[:np.argmax(T) + 1] = 0
        else:
            # massratio tests for the molecular composition of the solids:
            # only test, when there are solids!
            assert np.all(np.isclose(np.sum(chem_solid[:, 0, :], axis=1),
                                     np.ones_like(T))), "error in chemistry, total sum of solid abundancies is not 1"

        # massratio tests for the molecular composition of the gas:
        assert np.all(np.isclose(np.sum(chem_gas[:, 0, :], axis=1),
                                 np.ones_like(T))), "error in chemistry, total sum of gas abundancies is not 1"
        ###############
        # Sanity checks
        ###############

        TotH_mol_test = (
                self.abuCH4 * 4.0
                + self.abuH2O * 2.0
                + self.abuNH3 * 3
                + self.abuH2S * 2
        )
        abu_mol_test = (self.abu_array * self.M_array.squeeze()).sum() - TotH_mol_test
        abu_asp_test = (self.elabu_init * elmasses).sum()
        assert np.isclose(abu_asp_test, abu_mol_test), "error in chemistry, total sum of molecular abundancies is not equal to the total sum of the X/H stellar vaules. Partioning is wrong."

        # Test if the heavy elements are identical to the heavy molecules (gas):
        _, abu_gas_test = self.abundances([np.max(T_array)+100])
        chem_gas_test, restH = calc_composition(abu_gas_test, solid=False, He_Mbg=self.He_Mbg)
        elabus_heavy_chem_test = np.sum(chem_gas_test[0,0,:-7])
        assert np.isclose(chem_gas_test[0, 1, 1:].sum(), elabus_heavy_chem_test+restH), "error in chemistry, total sum of elemental mass abundancies is not equal to the total sum of molecular mass abundancies."

        # Test validity of dust to gas ratio:
        assert np.isclose(chem_gas_test[0,1,1:].sum()*(1+self.dtg0),self.dtg0), "Error in dust to gas ratio"

        # Test if the heavy elements are identical to the heavy molecules (solids):
        if np.max(T) < np.max(T_array):
            # only test, when there are solids!
            abu_solid_test, _ = self.abundances([0])
            chem_solid_test, restH = calc_composition(abu_solid_test, solid=True, He_Mbg=self.He_Mbg)
            elabus_heavy_chem_test = np.sum(chem_solid_test[0,0,:-7])
            assert np.isclose(chem_solid_test[0, 1, 1:].sum(), elabus_heavy_chem_test+restH), "error in chemistry, total sum of elemental mass abundancies is not equal to the total sum of molecular mass abundancies."

        return chem_solid, chem_gas

    def abundances(self, T: np.array) -> Tuple[np.array, np.array]:
        """
        Simple chemistry model that calculates abundancies of molecules
        See Bitsch & Battistini 2020

        Parameters
        ----------
        T: temperature array in cgs

        Returns
        -------
        Masses_solid: Array of solid molecular abundancies
        Masses_gas: Array of gas molecular abundancies

        """
        abu_array = np.ones_like(T) * self.abu_array[:, np.newaxis]

        M_masked_solid = np.where(T < T_array, np.ones_like(T) * self.M_array, np.zeros_like(abu_array))
        M_masked_solid[0, :] = 0  # exclude restmol
        Tot_M_solid = np.sum(M_masked_solid * abu_array, axis=0) + 1e-300
        Masses_solid = np.where(T < T_array, abu_array * M_masked_solid / Tot_M_solid,
                                np.zeros_like(abu_array))

        Z = self.get_solid_heavies(T)

        M_masked_gas = np.where(T >= T_array, np.ones_like(T) * self.M_array, np.zeros_like(abu_array))
        M_masked_gas[0, :] = 0  # exclude restmol
        Tot_M_gas = np.sum(M_masked_gas * abu_array, axis=0) + 1e-300
        Masses_gas = np.where(T >= T_array, abu_array * M_masked_gas / Tot_M_gas,
                              np.zeros_like(abu_array))
        Masses_gas[0, :] = 1 / (1 + (self.dtg0 - Z))
        Masses_gas[1:, :] = Masses_gas[1:, :] * (self.dtg0 - Z) / (1 + (self.dtg0 - Z))

        return Masses_solid, Masses_gas
