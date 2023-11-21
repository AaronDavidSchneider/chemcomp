"""
Common functions for chemcomp testing
"""

import pytest

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