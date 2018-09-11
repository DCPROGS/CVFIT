import os
import numpy as np

from cvfit import data
from cvfit.equations import Linear
from cvfit.fitting import SingleFitSession

def test_regression_linear_fit():
    filename = "./Example/Example.xlsx"
    set0 = data.read_sets_from_Excel(filename, 2, 0, 2)[0]
    equation = Linear('Linear')
    fsession = SingleFitSession(set0, equation)
    fsession.fit()
    fsession.calculate_errors()
    assert fsession.ndf == 8
    np.testing.assert_almost_equal(fsession.eq.pars[0], 3.666, 3)
    np.testing.assert_almost_equal(fsession.eq.pars[1], 5.204, 3)
    