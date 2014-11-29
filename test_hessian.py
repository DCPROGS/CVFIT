#! /usr/bin/python

from math import pi, log, sqrt, fabs
import numpy as np
from scipy import stats
from scipy import optimize
from cvfit import errors
from cvfit import cfio
from cvfit.equations import hill_equation, SSD, SSDlik, SSDlik_contour
import copy

__author__="User"
__date__ ="$27-Oct-2014 19:04:19$"

if __name__ == "__main__":
    
#    filename = cfio.ask_for_file()
    filename = "Example/Example.csv"
    dataset = cfio.read_sets_from_csv(filename)[2]
    print dataset
    
    # fit: minimise SSD function and get optimal parameters
    #pars = np.array([0.1e-19, 11124.1, 773.322, 1.23857])#0
    #pars = np.array([0.1e-19, 8033.44, 758.689, 1.02105])#1
    pars = np.array([0.1e-19, 19237, 425.356, 1.35852])#2
    #pars = np.array([0.1e-19, 17073.8, 408.942, 1.34237])#3
    #pars = np.array([0.1e-19, 9616.72, 615.511, 1.11186])#4
    names = ['Ymin', 'Ymax', 'EC50', 'nH']
    notfixed = np.array([False, True, True, True])
    fixval = copy.deepcopy(pars)
    Smin = SSD(pars, hill_equation, dataset)
    print '\n SSD \n', Smin
    hes = errors.hessian(dataset, hill_equation, pars, notfixed)
    print '\n Observed information matrix = \n', hes
    covar = errors.covariance_matrix(dataset, hill_equation, pars, notfixed) 
    print '\n Covariance matrix = \n', covar
    aproxSD = errors.approximateSD(dataset, hill_equation, pars, notfixed)
    #print '\n aproximate SD = \n', aSD
    
    print '\n'
    j = 0
    for i in range(len(names)):
        str = 'Parameter {0:d}: {1}\t= {2:.3f}\t'.format(i+1, names[i], pars[i])
        if notfixed[i]:
            str += 'Approx SD = {0:.3f}'.format(aproxSD[j])
            j += 1
        else:
            str += '(fixed)'
        print str
    
    correl = errors.correlation_matrix(covar)
    print '\n Correlation matrix = \n', correl
    
    
    # Number of points fitted = nfit
    nfit = dataset.size()
    # Number of parameters estimated = kfit
    kfit = np.nonzero(notfixed)[0].size
    # Degree of freedom = ndf
    ndf = nfit - kfit
    print('\nNumber of point fitted = {0:d}'.format(nfit))
    print('Number of parameters estimated = {0:d}'.format(kfit))
    print('Degrees of freedom = {0:d}'.format(ndf))
    
    var = Smin / ndf
    Sres, Lmax = sqrt(var), -SSDlik(pars, hill_equation, dataset, notfixed)
    print ('\nResidual error SD = {0:.3f} (variance = {1:.3f})'.format(Sres, var))
    print ('Minimum SSD = {0:.3f}; \tMax log-likelihood = {1:.3f}'.format(Smin, Lmax))
    
    tval = errors.tvalue(ndf)
    m = tval * tval / 2.0
    print "\ntval = {0:.3g}; m = {1:.3g}".format(tval, m)
    clim = sqrt(2. * m)
    print 'Equivalent SD for Gaussian:', clim
    Lcrit = Lmax - m
    print 'Lmax=', Lmax, 'Lcrit=', Lcrit
    
    

    Llimits = errors.lik_intervals(names, pars, aproxSD, 
        notfixed, m, hill_equation, dataset)
    print Llimits