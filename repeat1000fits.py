#! /usr/bin/python

__author__ = "remis"
__date__ = "$30-May-2015 21:07:31$"

import numpy as np
import pylab as plt

from cvfit.equations import Hill 
from cvfit.data import XYDataSet
from cvfit.fitting import SingleFitSession

if __name__ == "__main__":
    conc = np.array([0.1, 0.3, 1.0, 3.0, 10.0]) # , 30.0, 100.0
    sd = 0.1
    res = []
    for i in range(100):
        eq = Hill('Hill', pars=[0.0, 1.0, 1.0, 1.0])
        resp = eq.calculate_random(conc, sd)
        set = XYDataSet()
        set.from_columns(conc, resp)
        #print (set)

        fs = SingleFitSession(set, eq)
        fs.fit()
        res.append(eq.pars)
    results = np.array(res)
    
    plt.subplot(231)
    plt.hist(results[:, 1])
    plt.xlabel('Ymax')
    plt.ylabel('Number')
    plt.xlim(0,2)
    
    plt.subplot(232)
    plt.hist(results[:, 2])
    plt.xlabel('K')
    plt.ylabel('Number')
    plt.xlim(0,7)
    
    plt.subplot(233)
    plt.hist(results[:, 3])
    plt.xlabel('nH')
    plt.ylabel('Number')
    plt.xlim(0,7)
    
    plt.subplot(234)
    plt.plot(results[:, 1], results[:,2], 'b.')
    plt.xlabel('Ymax')
    plt.xlim(0,2)
    plt.ylabel('K')
    plt.ylim(0,7)
    plt.grid(True)

    plt.subplot(235)
    plt.plot(results[:, 1], results[:,3], 'b.')
    plt.xlabel('Ymax')
    plt.xlim(0,2)
    plt.ylabel('nH')
    plt.ylim(0,7)
    plt.grid(True)
    
    plt.subplot(236)
    plt.plot(results[:, 2], results[:,3], 'b.')
    plt.xlabel('K')
    plt.xlim(0,7)
    plt.ylabel('nH')
    plt.ylim(0,7)
    plt.grid(True)
    
    plt.show()