#!/usr/bin/env python

""" CVFIT is a library for weighted least-squares fitting
    of various equations to experimental data, for calculating 
    errors of fitting estimates and for plotting the results. 
    Currently available equations: linear, Hill, ... 
    Calculated are: approximate SD and 95% confidence limits; 
    likelihood intervals will be available too. 
"""

import math
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import binom
from scipy.stats.distributions import t
import matplotlib.pyplot as plt

###################   Equations   ###################################
def linear(x, a, b):
    return a + b * x 

def popen_dCK3(c, E, K):
    """Mechanism of sequential binding of n molecules followed by opening.
    R <-> AR <-> A2R <-> A3R <-> AnR*
    """
    Eterm = E * c**3 / K**3
    Kterm = np.sum([binom(3, r)*(c/K)**r for r in range(1, 4)])
    return  Eterm / (1 + Kterm + Eterm)

def popen_dCK2(c, E, K):
    """Mechanism of sequential binding of n molecules followed by opening.
    R <-> AR <-> A2R <-> AnR*
    """
    Eterm = E * c**2 / K**2
    Kterm = np.sum([binom(2, r)*(c/K)**r for r in range(1, 2)])
    return  Eterm / (1 + Kterm + Eterm)

def GHK(x, r, totOut, In1, In2):
    """Goldman-Hodgkin-Katz equation for bi-ionic condition"""
    return 25 * np.log( (x + r * (totOut - x)) / (In1 + r * In2) )

def Hill(x, Ymin, Ymax, EC50, nH):
    return (Ymin + ((Ymax - Ymin) * (x / EC50) ** nH) / 
        (1 + (x / EC50) ** nH))

def Langmuir(x, Ymin, Ymax, EC50): # Hill slope = 1
    return (Ymin + ((Ymax - Ymin) * (x / EC50)) / (1 + (x / EC50)))

def non_stationary_noise(x, i, N, s0):
    """Parabolic function: ensemble variance vs the mean current- used in 
       non-stationary noise analysis.
       Parameters: i- single channel current amplitude; 
       N- available number of channels; s0- variance of background noise.
    """
    return i * x - x**2 / N + s0

###################   Fit, errors   #################################
def fit(equation, X, Y, theta, S=None): #, limits=None):
    """Use scipy curve_fit."""
    sigma_is_used = True if S is not None else False
    estimates, covariance = curve_fit(equation, X, Y, 
        p0=theta, absolute_sigma=sigma_is_used, sigma = S) #, bounds=limits)
    return estimates, covariance

def errors(estimates, covariance, dof, alpha=0.05):
    """Calculate approximate SD and confidence limits. Default: 95%"""
    tval = t.ppf(1.0-alpha/2., dof) #student-t value
    sigma = np.sqrt(np.diag(covariance))
    confidence_limits = np.array([estimates - sigma * tval, 
                                  estimates + sigma * tval]).T
    return sigma, confidence_limits

def dof(X, theta): # degree of freedom
    return max(0, len(X) - len(theta)) # #datapoints - #parameters

def print_estimates_errors(names, estimates, SD, confidence_limits):
    for name, p, v, c in zip(names, estimates, SD, confidence_limits):
        print('{0}: {1} +/- {2} [{3}  {4}]'.
            format(name, p, v, c[0], c[1]))

def single_curve_fit(equation, X, Y, theta, names):
    estimates, covariance = fit(equation, X, Y, theta) #, limits=(0, [100., 100.]))
    approxSD, conf95 = errors(estimates, covariance, dof(X, theta), alpha=0.05)
    print_estimates_errors(names, estimates, approxSD, conf95)
    return estimates

######################   Examples   #################################
def example_linear():
    X = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5]
    Y = [3.17, 13.25, 19.8, 14.18, 11.43, 25.85, 13.81, 25.49, 26.94, 38.86]
    theta, names = [2, 6], ['a', 'b']
    single_curve_fit(linear, X, Y, theta, names)

def example_popen_dCK3():
    E, K = 10, 10
    X = np.array([1, 2, 5, 10, 50])
    Ytrue = popen_dCK3(X, E, K)
    Y = Ytrue + 0.01 * np.random.normal(size=X.size)
    theta, names = [9, 11], ['E', 'K']
    single_curve_fit(popen_dCK3, X, Y, theta, names)

#####################################################################
if __name__ == '__main__':
    example_linear()
    #example_popen_dCK3()
