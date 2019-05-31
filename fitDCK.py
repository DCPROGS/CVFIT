#!/usr/bin/env python

import math
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import binom
import matplotlib.pyplot as plt

def popen(c, E, K):
    '''
    Mechanism of sequential binding of n molecules followed by opening.
    R <-> AR <-> A2R <-> ... <-> AnR <-> AnR*
    '''
    Eterm = E * c**3 / K**3
    Kterm = np.sum([binom(3, r)*(c/K)**r for r in range(1, 4)])
    return  Eterm / (1 + Kterm + Eterm)

E, K = 10, 10
plotX = np.linspace(0.1, 100, 1000)
plotY = popen(plotX, E, K)
X = np.array([1, 2, 5, 10, 50])
Ytrue = popen(X, E, K)
noise = 0.005 * np.random.normal(size=X.size)
Y = Ytrue + noise
#plt.semilogx(plotX, plotY); # curve
#plt.semilogx(X, Y, 'o'); # points

theta = [10, 10]
popt, pcov = curve_fit(popen, X, Y, p0=theta, bounds=(0, [100., 100.]))
print(popt)
perr = np.sqrt(np.diag(pcov))
print(perr)
print(pcov)