__author__ = "User"
__date__ = "$03-Apr-2015 17:41:27$"

import os
from cvfit import data
from cvfit.equations import Linear as eqfit
from cvfit.fitting import SingleFitSession

if __name__ == "__main__":

    filename = "./Example/Example.xlsx"
    set0 = data.read_sets_from_Excel(filename, 2, 0, 2)[0]
    print("Loaded: " + os.path.split(str(filename))[1])
    print (str(set0))
    equation = eqfit('Linear')
    fsession = SingleFitSession(set0, equation)
    fsession.fit()
    fsession.calculate_errors()
    print(fsession.string_estimates())
    fsession.calculate_errors()
    print(fsession.string_liklimits())
