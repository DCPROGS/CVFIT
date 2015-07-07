__author__ = "User"
__date__ = "$03-Apr-2015 17:41:27$"

import os
import pylab as plt

from cvfit import data
from cvfit.equations import GHK as eqfit
from cvfit.fitting import SingleFitSession
from cvfit.data import XYDataSet

if __name__ == "__main__":

    filename = "./Example/Example.xlsx"
    set0 = data.read_sets_from_Excel(filename, 2, 0, 3)[0]
    print("Loaded: " + os.path.split(str(filename))[1])
    print (str(set0))
    equation = eqfit('GHK', pars=[1.0, 150.0, 145.0, 5.0])
    fsession = SingleFitSession(set0, equation)
    fsession.fit()
    fsession.calculate_errors()
    print(fsession.string_estimates())
    fsession.calculate_errors()
    print(fsession.string_liklimits())
    
    plX, plY = equation.calculate_plot(set0.X, equation.pars)
    rlim = fsession.Llimits[0]
    plX1, plY1 = equation.calculate_plot(set0.X, [rlim[0], 150.0, 145.0, 5.0])
    plX2, plY2 = equation.calculate_plot(set0.X, [rlim[1], 150.0, 145.0, 5.0])
    
    #avdat = XYDataSet()
    #avdat.pool(set0.X, set0.Y, seto.S)
    set0.average_pooled()
    
    plt.subplot(111)
    plt.plot(set0.X, set0.Y, 'ro') # all data points
    plt.errorbar(set0.avX, set0.avY, yerr=set0.avS, fmt='o')
    plt.plot(plX, plY, 'b-') # fit
    plt.plot(plX1, plY1, 'k--') # lower likelihood limit
    plt.plot(plX2, plY2, 'k--') # higher likelihood limit
    plt.xlabel('Kout, mM')
    plt.ylabel('Erev, mV')
    plt.show()