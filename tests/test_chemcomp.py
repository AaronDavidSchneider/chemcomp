"""
General GCMT tests
"""
import numpy as np
import pytest
from chemcomp.helpers import GrowthPlot
import tables
import astropy.units as u
import astropy.constants as const
import os

from test_chemcomp_common import all_references_testdata


def test_compare_reference_minimal(all_references_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on it."""
    name, expected = all_references_testdata

    files = [f'output/{name}.h5', expected['path']]
    out = ['new', 'old']
    chemcomp = GrowthPlot(files=files, out=out)

    old = chemcomp.data['old']
    new = chemcomp.data['new']

    # planet stuff
    for q in ['M', "M_a", "M_c", "M_z_gas", "M_z_peb", "comp_a", "comp_c", "pebiso"]:
        np.testing.assert_allclose(old[q], new[q], err_msg=f"failed on {q}")

    # disk stuff
    data = {}
    AU = (1 * u.au).cgs.value  # define the astronomical unit conversion from cgs
    Myr = (1 * u.Myr).cgs.value  # define the Megayear unit conversion from cgs

    for o, f in zip(out, files):
        with tables.open_file(f, mode='r') as f:
            data[o] = {}
            data[o]['sigma_g'] = np.array(f.root.disk.sigma_g)
            data[o]['sigma_dust_components'] = np.array(f.root.disk.sigma_dust_components)
            data[o]['t'] = np.array(f.root.disk.t) / Myr  # convert to Myr
            data[o]['r'] = np.array(f.root.disk.r) / AU  # convert to AU

    all_keys = data[out[0]].keys()

    # mask out all the data at which we have infinitesimally small values
    # -> numerical noise may be important there
    threshold = 1e-60
    for key in all_keys:
        old_data = data[out[0]][key]
        new_data = data[out[1]][key]
        mask = np.where(np.logical_and(old_data > threshold, new_data > threshold))
        np.testing.assert_allclose(old_data[mask], new_data[mask], err_msg=f"failed on {key}")










