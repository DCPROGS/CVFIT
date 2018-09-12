import os
import unittest
import numpy as np

from cvfit import data
from cvfit.equations import Linear
from cvfit.fitting import SingleFitSession

class TestRegressionLinearFit:
    """ """
    def setUp(self):
        filename = "./Example/Example.xlsx"
        set0 = data.read_sets_from_Excel(filename, 2, 0, 2)[0]
        equation = Linear('Linear')
        self.fs = SingleFitSession(set0, equation)
        self.fs.fit()
        self.fs.calculate_errors()
    
    def test_degree_of_feedom(self):    
        assert self.fs.ndf == 8

    def test_fitted_parameters(self):
        np.testing.assert_almost_equal(self.fs.eq.pars[0], 3.666, 3)
        np.testing.assert_almost_equal(self.fs.eq.pars[1], 5.204, 3)

    def test_approximate_SD(self):
        np.testing.assert_almost_equal(self.fs.aproxSD[0], 5.216871572087441, 9)
        np.testing.assert_almost_equal(self.fs.aproxSD[1], 1.5729459621912463, 9)

    def test_Smin(self):
        np.testing.assert_almost_equal(self.fs.Smin, 395.8654399999999, 9)

    def test_maxLogLik(self):
        np.testing.assert_almost_equal(self.fs.Lmax, -32.581831644727835, 9)

    def test_lik_limits(self):
        self.fs.Llimits[0][0] is None
        np.testing.assert_almost_equal(self.fs.Llimits[0][1], 16.001557751460275, 9)
        np.testing.assert_almost_equal(self.fs.Llimits[1][0], 1.4776049950866297, 9)
        np.testing.assert_almost_equal(self.fs.Llimits[1][1], 8.923310603763344, 9)
    