__author__ = "User"
__date__ = "$03-Apr-2015 17:41:27$"

import os
from cvfit import data
from cvfit.equations import Linear as eqfit
from cvfit.fitting import SingleFitSession

if __name__ == "__main__":

    filename = "./Example/Example.xlsx"
    allsets = data.read_sets_from_Excel(filename, 2, 0, 2)
    print("Loaded: " + os.path.split(str(filename))[1])
    print (str(allsets[0]))
    eqname = 'Linear'
    equation = eqfit(eqname)
    fsession = SingleFitSession(allsets[0], equation)
    fsession.fit()
    fsession.calculate_errors()
    print(fsession.string_estimates())
    fsession.calculate_errors()
    print(fsession.string_liklimits())
