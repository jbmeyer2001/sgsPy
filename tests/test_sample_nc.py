import pytest
import sgspy as sgs

class TestNc:
    #TODO add real tests before implementing nc
    @pytest.mark.xfail
    def test_not_implemented(self):
        sgs.sample.nc()
