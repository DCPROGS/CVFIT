"""
Equations available to fit data in CVFIT.
"""
from math import sqrt, log, pi
import copy
import numpy as np
from scipy import optimize

def SSDlik(theta, func, data, notfixed, allpars=None):
    """
    Calculate likelihood coresponding to the sum of the squared deviations 
    assuming that errors follow Gausian distribution.
    """

    if allpars != None:
        pars = []
        j = 0
        for i in range(len(allpars)):
            if notfixed[i]:
                pars.append(theta[j])
                j += 1
            else:
                pars.append(allpars[i])
    else:
        pars = theta
    S = 0.0
    for point in data.points:
        S += (point.y - func(pars, point.x)) ** 2
    nfit = data.size()
    kfit = np.nonzero(notfixed)[0].size
    Sres = sqrt(S / (nfit - kfit))
    return nfit * log(sqrt(2 * pi) * Sres) + S / (2 * Sres**2)

def SSDlik_contour(x, num, pars, notfixed, func, dataset):
    fixfix = copy.deepcopy(notfixed)
    fixfix[num] = False
    allpars = pars.copy()
    allpars[num] = x
    theta = allpars[np.nonzero(fixfix)[0]]
    result = optimize.minimize(SSDlik, theta, 
        args=(func, dataset, fixfix, allpars), 
        method='Nelder-Mead', jac=None, hess=None)
    return -result.fun, result.x

def SSD(args, func, data):
    """
    Calculate sum of squared deviations.
    """
    S = 0.0
    for point in data.points:
        S += (point.y - func(args, point.x)) ** 2
    return S

def hill_equation(args, conc):
    '''
    The hill equation.
    '''
    ymin, ymax, ec50, nH = args
    return ymin + (ymax * (conc / ec50) ** nH) / (1 + (conc / ec50) ** nH)
    



