import pytest
from addons import *

@ctest_labeler("quick;cbs")
@pytest.mark.parametrize("distributed", [
    pytest.param(False, id="internal"),
    pytest.param(True,  id="snowflake", marks=using("qcfractal_next")),
])
def test_cbs_xtpl_energy(distributed):
    setenv = ["_PSI4_USE_QCF"] if distributed else None

    ctest_runner(__file__, setenv=setenv)
