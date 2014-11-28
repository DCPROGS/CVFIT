#! /usr/bin/python

import numpy as np
from cvfit import errors
from cvfit import cfio
from cvfit import equations
import copy

__author__="User"
__date__ ="$27-Oct-2014 19:04:19$"

if __name__ == "__main__":
    
#    filename = cfio.ask_for_file()
    filename = "Example/Example.csv"
    dataset = cfio.read_sets_from_csv(filename)[4]
    print dataset
    
    # fit: minimise SSD function and get optimal parameters
    #pars = np.array([0.1e-19, 11124.1, 773.322, 1.23857])#0
    #pars = np.array([0.1e-19, 8033.44, 758.689, 1.02105])#1
    #pars = np.array([0.1e-19, 19237, 425.356, 1.35852])#2
    #pars = np.array([0.1e-19, 17073.8, 408.942, 1.34237])#3
    pars = np.array([0.1e-19, 9616.72, 615.511, 1.11186])#4
    names = ['Ymin', 'Ymax', 'EC50', 'nH']
    notfixed = np.array([False, True, True, True])
    fixval = copy.deepcopy(pars)
    minssd = equations.SSD(pars, func, dataset, notfixed, fixval)
    print '\n SSD \n', minssd
    hes = errors.hessian(dataset, equations.hill_equation, pars, notfixed, fixval)
    print '\n Observed information matrix = \n', hes
    covar = errors.covariance_matrix(dataset, equations.hill_equation, pars, notfixed, fixval) 
    print '\n Covariance matrix = \n', covar
    aSD = errors.approximateSD(dataset, equations.hill_equation, pars, notfixed, fixval)
    print '\n aSD \n', aSD
    correl = errors.correlation_matrix(covar)
    print '\n Correlation matrix = \n', correl
    
    
    #lowli = errors.confidence_interval(covar, dataset, equations.hill_equation, pars, notfixed)
    ssec, ssec1, ssec2, count = errors.new_secs3(aSD, covar, dataset, equations.hill_equation, pars, notfixed, hes)
    print '\n ssec = \n', ssec, ssec1, ssec2
    print '\n count = \n', count    
    #print '\n aasdlist = \n', aasdlist  
    fitpars = pars[1:]
    names = names[1:]
    for name, var, sd in zip(names, fitpars, aSD):
        print ('\n{0} = {1:.3g} +/- {2:.3g}'.format(name, var, sd))
    