#! /usr/bin/python

__author__="remis"
__date__ ="$05-Nov-2010 21:17:34$"

import numpy as np
from math import*

#import dcSimplex
import fitStats
import fitPlot
#import dialogHill

try:
    from PySide.QtGui import *
    from PySide.QtCore import *
except:
    raise ImportError("pyqt module is missing")

class equationHill(object):
    """
    Hill equation.
    """

    def __init__(self):

        self.opts = self.options()
        #self.curves = []
        #self.guesses = {}
        #self.dicdat = {}
        #self.result = {}


    def setGuesses(self, guesses):
        self.guesses = guesses

    def calcGuesses(self):
        ""

        nsets = self.dicdat['nsets']
        sets = self.dicdat['sets']
        data = self.dicdat['data']

        guesses = {}
        inis = []
        fixed = []
        minmax = []

        ymin = 1e+37
        ymax = 1e-37
        xmin = 1e+37
        xmax = 1e-37
        ymaxsep = []

        for i in range(nsets):
            ymaxsep.append(1e-37)

        for j in range(nsets):
            name = sets[j]
            set = data[name]
            minx = min(set[0])
            if minx < xmin: xmin = minx
            maxx = max(set[0])
            if maxx > xmax: xmax = maxx

            theta, fix = self.guesses(set)

            inis.append(theta[0])
            inis.append(theta[1])
            inis.append(theta[2])
            inis.append(theta[3])
            fixed.append(fix[0])
            fixed.append(fix[1])
            fixed.append(fix[2])
            fixed.append(fix[3])

        minmax.append(xmin)
        minmax.append(xmax)
        minmax.append(ymin/2)
        minmax.append(ymax)

        guesses['minmax'] = minmax
        guesses['ymaxsep'] = ymaxsep
        guesses['inis'] = inis
        guesses['fixed'] = fixed

        return guesses

    def guesses(self, data, verbose=0):
        """
        """

        y0 = 0
        ymax = max(data[1])    # * 1.05

        n = 0
        nH = 1
        K = 1
        X1 = []
        Y1 = []
        for i in range(np.size(data[0])):
            if data[1,i] > 0.1 * ymax and data[1,i] < 0.9 * ymax:
                p = data[1,i] / (ymax - data[1,i])
                n = n + 1
                X1.append(log(data[0,i]))
                Y1.append(log(p))
        if n > 1:
            a, b = fitStats.linreg(X1, Y1, n)
            nH = b
            K = exp(-a/b)
            if verbose: print 'a, b=', a, b
        if verbose: print 'K, nH= ', K, nH

        ymax = ymax * 1.05

        theta = np.array([y0, ymax, K, nH])
        fix = np.array([1, 0, 0, 0])
        return theta, fix

    def setData(self, series_to_fit, data):
        """
        """
        dicdat = {}
        dicdat['data'] = data
        dicdat['nsets'] = len(series_to_fit)
        self.nsets = len(series_to_fit)
        dicdat['sets'] = series_to_fit

        self.dicdat = dicdat

    def testData(self):
        """
        One set of open probabilities at several concentrations.
        Data- wt muscle ACh receptor @ Vm = -100 mV, cell attached (RL)
        Data array structure: row0- x1, row1- y1, row2- sd1, row3- will
        be generated later as weights. Four rows for every set.

        """
        data = np.array([[3.00, 5.00, 6.40, 10.0, 20.0, 30.0, 64.0, 100.0],
                          [0.25, 0.39, 0.55, 0.69, 0.85, 0.92, 0.93, 0.930],
                          [0.045, 0.078, 0.085, 0.052, 0.031, 0.014, 0.012, 0.01],
                          [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.000]])

        dicdat = {}
        dicdat['data'] = data
        dicdat['nsets'] = 1
        dicdat['sets'] = 'set#1'

        return dicdat

    def options(self):
        """
        """
        dict = {}
        dict['nmod'] = 26      # id# of a model in DC's CVFIT
        dict['ncomp'] = 1    # number of sets
        dict['nset'] = 1     # number of sets
        dict['increase'] = True    # function is rising if True; Y increases with X
        dict['fline'] = False    # add linear comp if true
        dict['ifitmode'] = 1     #
        dict['logyfit'] = False

        dict['title'] = "Hill Equation"
        dict['is_y0fixed'] = True
        dict['are_normalised'] = False
        dict['common_list'] = ["Fit separate maximum, separate K values",
            "Fit common maximum, separate K values",
            "Fit common K value, separate maxima",
            "Separate K and Ymax, but common percent Ymax1"]
        dict['common'] = 0

        return dict

    def squeeze(self, guess, fix):

        j = 0
        for i in range(np.size(guess)):
            if fix[i] == 0:
                j += 1

        m = 0
        theta = np.zeros((j))
        for i in range(np.size(guess)):
            if fix[i] == 0:
                theta[m] = guess[i]
                m +=1

        return theta

    def unsqueeze(self, theta, guess, fix):

        j = 0
        for i in range(np.size(guess)):
            if fix[i] == 0:
                guess[i] = theta[j]
                j +=1
        return guess


    def hill(self, theta, xv, opt):
        """
        nmod=26 (hill), nmod=27 (langmuir, nH=1)
        """

        # no common params, one set         ifitmode NOT= 3 (ifitmode.ne.3)
        ycalx = 0.0
        #y0 = theta[0]
        y0 = 0
        ymax = 0.0
        K = 0.0
        nH = 0.0
        j = -1
        n = opt['ncomp']
        for i in range(0, n):
            j = j+1
            ymax = theta[j]
            j = j+1
            K = theta[j]
            c = xv / K
            if opt['nmod'] == 26:
                j = j+1
                nH = theta[j]
                c = pow(c, nH)
            if opt['increase']:    # increasing curve
                x = c / (c + 1.0)
            else:           # decreasing curve
                x = 1.0 / (1.0 + c)
            ycalx = ycalx + ymax * x

        ycalx = y0 + ycalx
        if opt['fline']:
            slope = theta[-1]
            ycalx = ycalx + slope * xv

        return ycalx


    def fitHill(self, theta, data):
        f = fitStats.ssdcv(theta, data, self.opts, self.hill)
        return f

    def initiateFit(self):

        sets = self.dicdat['sets']
        data = self.dicdat['data']
        inis = self.guesses['inis']
        fixed = self.guesses['fixed']
        theta = self.squeeze(inis, fixed)
        result = inis
        values = np.zeros((np.size(theta)))
        sd = []

        for i in range(self.nsets):
            theta1 = np.zeros((3))
            name = sets[i]
            set = data[name]
            theta1 = theta[i*3 : i*3+3]
            values1, smin = fitStats.simplex(theta1, self.fitHill, set)
            values[i*3 : i*3+3] = values1
            sd.append(smin)

        result = self.unsqueeze(values, result, fixed)

        self.result = result
        self.sd = sd

    def calcCurves(self):

        xmin1 = self.guesses['minmax'][0]
        xmax1 = self.guesses['minmax'][1]
        xmin = int(log10(xmin1)-2)
        xmax = int(log10(xmax1)+1)
        nPoint = 1000
        cmin = pow(10, xmin)    #0.001
        cmax = pow(10, xmax)    #1
        dc = (log10(cmax)-log10(cmin)) / float(nPoint-1)

        self.curves = []
        c = np.zeros((nPoint))
        self.curves.append(c)
        for i in range(nPoint):
            c[i] = cmin * pow(10, (i*dc))
        for i in range(self.nsets):
            pars = self.result[i*4+1 : i*4+4]
            f = np.zeros((nPoint))
            for j in range(nPoint):
                f[j] = self.hill(pars, c[j], self.opts)
            self.curves.append(f)

if __name__ == "__main__":

    print "Testing Python implementation of DC's CVFIT program.\n"

    eh = equationHill()
    
    dicdata = eh.testData()
    data = dicdata['data']
    nobs = np.size(data[0])
    data = fitStats.calcSDWeights(data)

    print 'In this test weights from specified SD will be used for fits.\n'
    print '\n Obs#  X value   Y value   sdm(Y)       weight'
    for i in range(0, nobs):
        print ' ',i+1,'   ', data[0,i], '   ', data[1,i], '   ', data[2,i], '   ', data[3,i]

    opts = eh.options()
    print '\nFunction options:', opts

    guess, fix = eh.guesses(data)
    theta = eh.squeeze(guess, fix)

    f = fitStats.ssdcv(theta, data, opts, eh.hill)
    print '\nInitial guesses: SSD =', f
    print theta

    k = np.size(theta)
    nobs = np.size(data[0])
    ndf = nobs - k
    print '\n Number of free parameters =', k
    print ' Number of different x values =', nobs

    values, smin = fitStats.simplex(theta, eh.fitHill, data)

    print ' END OF FITTING'

    print '\n SDmin =', smin
    print ' Final values:', values
    sder = math.sqrt(smin/float(ndf))
    print 'Error S.D.= sqrt(SDmin/df) =', sder, '(', ndf, 'degrees of freedom)'

    thetaDCig = np.array([0.976, 5.49, 1.75])
    f1 = fitStats.ssdcv(thetaDCig, data, opts, eh.hill)
    print "\n DC's initial guesses: SSD =", f1
    print thetaDCig

    thetaDCfv = np.array([0.938078, 5.4319, 1.93519])
    f2 = fitStats.ssdcv(thetaDCfv, data, opts, eh.hill)
    print "\n DC's CVFIT result: SSD =", f2
    print thetaDCfv

    plt1 = fitPlot.plot(values, data, opts, eh.hill)