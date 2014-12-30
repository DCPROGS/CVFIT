#! /usr/bin/python

"""
Python implementations of YCALCV.FOR
CALCULATES I'TH VALUE IN J'TH DATA SET FOR MODEL #NMOD IN CVFIT
function YCALCV(K,THETA,xv1,ival,j,setx,njset)
K- ; THETA- ; xv1- parameter;
ival- if ival<1 then value of XV in common is used;
      if ival>1 then use x=xv1 (parameter)
j- ; setx- ; njset-


Some DC's parameters:
If ifitmode=1 then only one set fitted (#jset, in common)
If ifitmode=2 then all sets fitted separately in loop, the current
    set again being specified by #jset, in common
If ifitmode=3-5 then all sets fitted simultaneously with single eqn

opt- a dictionary; can contain some or all options:
opt['ifitmode'] = 1, 2, 3, 4 or 5, see above;
opt['nmode'] - id# of a model in DC's CVFIT;
opt['ncomp'] - number of components;
opt['increase'] = True or False -function is rising if True;
opt['fline'] = True or False - add linear comp if true;
opt['nset'] - number of sets;
opt['nsfit'] - number of sets to fit;
opt['juse'] - an array containing id# of sets to fit
opt['logyfit'] -



"""

__author__="remis"
__date__ ="$22-Oct-2010 21:08:06$"

import matplotlib.pyplot as plt
import numpy as np
import math
#import Simplex
import dcSimplex
#import fitStats

# Variable scope
global theta
global data
global opt
global func


def testData():
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
    return data

def eqnHill(theta, xv, opt):
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
            c = math.pow(c, nH)
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

def linreg(x, y, n, imode = 0):
    """
    Python implementation of DC's subroutine LINREG(x,y,n,a,b,imode).
    Return a,b for conventional linear regression.
    imode=0 for normal fit.
    imode=1 for fit of slope with a=input value (e.g. a=0 or a=1 for Schild eq).
    """

    a = 1

    en = float(n)
    sx = 0.0
    sxx = 0.0
    sy = 0.0
    syy = 0.0
    sxy = 0.0
    for i in range(n):
        sx = sx + x[i]
        sxx = sxx + x[i] * x[i]
        sy = sy + y[i]
        syy = syy + y[i] * y[i]
        sxy = sxy + x[i] * y[i]
    sxx = sxx - (sx * sx) / en
    sxy = sxy - (sx * sy) / en
    b = sxy / sxx
    if imode == 0: a = (sy / en) - b * sx /en

    return a, b

def ssdcv(theta, data,opt, func):
    """
    Python implementation of SSDCV.FOR
    SUM OF SQUARES FOR FITTING MODEL #NMOD WITH CVFIT
    function SSDCV(K,THETA,Xobs,yobs,w,nj,juse,setx,niobs,njset)
    """

    s = 0.0

    if opt['logyfit']:
        pass
    

    if opt['ifitmode'] < 3:  # Only one set fitted or sets fitted separately
        n = np.size(data[0])     # number of points in a set
        for i in range(0, n):
            x = data[0,i]
            y = func(theta, x, opt)
            w = 1 / (data[2,i] * data[2,i])
            dev = data[1,i] - y
            s = s + w * dev * dev
        return s

    if opt['ifitmode'] >= 3 or opt['ifitmode'] <= 5:
        #fit all sets at once for ifitmode=3
        iset = 0
        for m in range(0, opt['nsfit']):
            j = opt['juse'][m]     # set # used
            iset = iset + 1
            # iset is in COMMON -nec for YCALCV when ifitmode=4,5 if some sets missed so
            # c parameter # in theta() can be matched with proper data set
            for i in range(0, opt['nj'][j]):
                # need some coding if wangt some observations to omit
                # if(nomit.gt.0) then
                #   omit specified obs. Is I equal to any of the elements of JOMIT?
		#   do L=1,nomit
		#	if(i.eq.jomit(L)) goto 4
		#   enddo
		# endif
		# xv1 = Xobs(i,j)
		# Y = ycalcv(kmax,theta,xv1,ival,j,setx,njset)
                # dev = Yobs(i,j) - Y     #normal calc
		#S = S + w(i,j)*dev*dev

                x = data[0,i]
                y = func(theta, x, opt)
                w = 1 / (data[2,i] * data[2,i])
                dev = data[1,i] - y
                s = s + w * dev * dev
		nfit = nfit + 1
        return s

def test(verbose=0):
    """

    """

    global theta
    global data
    global opt
    global func


    dict = {}
    dict['nmod'] = 26      # id# of a model in DC's CVFIT
    dict['ncomp'] = 1    # number of sets
    dict['nset'] = 1     # number of sets
    dict['increase'] = True    # function is rising if True
    dict['fline'] = False    # add linear comp if true
    dict['ifitmode'] = 1     #
    dict['logyfit'] = False

    data = testData()
    k = np.size(data[0])
    for i in range(0, k):
        data[3,i] = 1 / (data[2,i] * data[2,i])
    xmin = int(math.log10(min(data[0]))-2)
    xmax = int(math.log10(max(data[0]))+1)
    y0 = 0
    ymax = max(data[1])
    ymax1 = int(max(data[1])+1)

    

    n = 0
    nH = 1
    K = 1
    X1 = []
    Y1 = []
    for i in range(np.size(data[1])):
        if data[1,i] > 0.1 * ymax and data[1,i] < 0.9 * ymax:
            p = data[1,i] / (ymax - data[1,i])
            n = n + 1
            X1.append(math.log(data[0,i]))
            Y1.append(math.log(p))
    if n > 1:
        a, b = linreg(X1, Y1, n)
        nH = b
        K = math.exp(-a/b)
        if verbose: print 'a, b=', a, b
    if verbose: print 'K, nH= ', K, nH

    ymax = ymax * 1.05
    #ymax = 0.976
    #K = 5.49
    #nH = 1.75

    theta = np.array([ymax, K, nH])

    func = eqnHill

    Y = ssdcv(theta, data, dict, func)
    #print 'Y=', Y



    return data, theta, dict

def plot(theta, data, opt):
    """

    """

    xmin = int(math.log10(min(data[0]))-2)
    xmax = int(math.log10(max(data[0]))+1)
    y0 = 0
    ymax = max(data[1])
    ymax1 = int(max(data[1])+1)


    nPoint = 1000
    cmin = math.pow(10, xmin)    #0.001
    cmax = math.pow(10, xmax)    #1
    dc = (math.log10(cmax)-math.log10(cmin)) / (nPoint-1)


    c = np.zeros((1, nPoint))
    f = np.zeros((2, nPoint))
    for i in range(0, nPoint):
        c[0,i] = cmin * math.pow(10, (i*dc))
        f[0,i] = eqnHill(theta, c[0,i], opt)
        #f[1,i] = eqnHill(theta, c[0,i], 26, 1, 1, False, False)

    line1 = plt.semilogx(c[0], f[0],'r-') #, c[0], f[1], 'b-')
    line2 = plt.semilogx(data[0], data[1], 'bo')
    plt.ylabel('Response')
    plt.xlabel('Dose, M')
    plt.title('Hill/Langmuire equation')
    plt.show()


if __name__ == "__main__":

    print " Python implementation of DC's CVFIT program"

    data, theta, opt = test()
    k = np.size(data[0])
    print '\n Obs#  X value   Y value   sdm(Y)       weight'
    for i in range(0, k):
        print ' ',i+1,'   ', data[0,i], '   ', data[1,i], '   ', data[2,i], '   ', data[3,i]

    nobs = np.size(data[0])
    ndf = nobs - k

    print '\n Number of free parameters =', k
    print ' Number of different x values =', nobs

    f = ssdcv(theta, data, opt, func)
    print '\nInitial guesses: SSD =', f
    print theta

    #theta = np.array([0.938078, 5.4319, 1.93519])   # dc's cvfit result
    values, smin = fitStats.simplex(theta, eqnHill, ssdcv, opt, data)

    print ' END OF FITTING'

    print '\n SDmin =', smin
    print ' Final values:', values
    sder = math.sqrt(smin/float(ndf))
    print 'Error S.D.= sqrt(SDmin/df) =', sder, '(', ndf, 'degrees of freedom)'

    thetaDCig = np.array([0.976, 5.49, 1.75])
    f1 = ssdcv(thetaDCig, data, opt, func)
    print "\n DC's initial guesses: SSD =", f1
    print thetaDCig

    thetaDCfv = np.array([0.938078, 5.4319, 1.93519])
    f2 = ssdcv(thetaDCfv, data, opt, func)
    print "\n DC's CVFIT result: SSD =", f2
    print thetaDCfv

    plt1 = plot(values, data, opt)
