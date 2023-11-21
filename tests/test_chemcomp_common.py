"""
Common functions for chemcomp testing
"""
import pytest

# mask out all the data at which we have infinitesimally small values
# -> numerical noise may be important there
THRESHOLD_DISK = 1e-30
MAX_RTOL_DISK = 1e-3
MAX_RTOL_PLANET = 1e-4

reference = {
    "Bert": {
        "path": "ci/reference.h5",
        "config": "config/config.yaml"
    },
}

# dictionary of archived experiments and some expected properties
@pytest.fixture
def all_references_testdata():
    yield from reference.items()