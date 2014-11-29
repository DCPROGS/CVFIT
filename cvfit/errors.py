import numpy as np
from numpy import linalg as nplin
from equations import SSD, SSDlik, SSDlik_contour
from math import fabs, sqrt

def optimise_deltas(dataset, func, args, notfixed):
    """ """
    Lmax = -0.5 * SSD(args, func, dataset)
    Lcrit = 1.005 * Lmax
    deltas = np.zeros((len(args),))
    for index in np.nonzero(notfixed):
        deltas[index] += 0.1 * args[index]
    L = -0.5 * SSD(args + deltas, func, dataset)
    if L > Lcrit:
        count = 0
        while L > Lcrit and count < 100:
            for index in np.nonzero(notfixed):
                deltas[index] *= 2
            L = -0.5 * SSD(args + deltas, func, dataset)
            count += 1
    elif L < Lcrit:
        count = 0
        while L < Lcrit and count < 100:
            for index in np.nonzero(notfixed):
                deltas[index] *= 0.5
            L = -0.5 * SSD(args + deltas, func, dataset)
            count += 1
    #print '\n deltas = \n', deltas
    return deltas

def hessian(dataset, func, pars, notfixed):
    """
    """
    kfit = np.nonzero(notfixed)[0].size
    hessian = np.zeros((kfit, kfit))
    deltas = optimise_deltas(dataset, func, pars, notfixed)
    i = 0
    for i1 in np.nonzero(notfixed)[0]:
        j = 0
        for j1 in np.nonzero(notfixed)[0]:
            coe1, coe2, coe3, coe4 = pars.copy(), pars.copy(), pars.copy(), pars.copy()

            if i == j:
                coe1[j1] += deltas[j1]
                coe3[j1] -= deltas[j1]
                hessian[i, j] = ((
                    SSD(coe1, func, dataset) -
                    2.0 * SSD(pars, func, dataset) +
                    SSD(coe3, func, dataset)) /
                    (deltas[j1]  ** 2))
            else:
                coe1[i1] += deltas[i1]
                coe1[j1] += deltas[j1]
                coe2[i1] += deltas[i1]
                coe2[j1] -= deltas[j1]
                coe3[i1] -= deltas[i1]
                coe3[j1] += deltas[j1]
                coe4[i1] -= deltas[i1]
                coe4[j1] -= deltas[j1]
                hessian[i, j] = ((
                    SSD(coe1, func, dataset) -
                    SSD(coe2, func, dataset) -
                    SSD(coe3, func, dataset) +
                    SSD(coe4, func, dataset)) /
                    (4 * deltas[i1] * deltas[j1]))
            j += 1
        i += 1
    return 0.5 * hessian

def covariance_matrix(dataset, func, pars, notfixed):
    """ """
    cov = nplin.inv(hessian(dataset, func, pars, notfixed))
    minssd = SSD(pars, func, dataset)
    #print '\n SSD \n', minssd
    kfit = np.nonzero(notfixed)[0].size
    #TODO: check if weights were used to calculate SSD. If not then
    # calculate errvar.
    errvar = minssd / (dataset.size() - kfit)
    return cov * errvar

def approximateSD(dataset, func, pars, notfixed):
    """ """
    cov = covariance_matrix(dataset, func, pars, notfixed)
    return np.sqrt(cov.diagonal())
    
def correlation_matrix(covar):
    correl = np.zeros((len(covar),len(covar)))
    
    for i1 in range(len(covar)):
        j1=0
        for j1 in range(len(covar)):
            correl[i1,j1] = covar[i1,j1]/np.sqrt(np.multiply(covar[i1,i1],covar[j1,j1]))
    return correl
            
def tvalue(ndf):
    """
    Return P=0.95 value of Student's t, with ndf degrees of freedom.
    """
    ttable = [12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365,
        2.306, 2.262, 2.228, 2.201, 2.179, 2.160, 2.145, 2.131,
        2.120, 2.210, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069,
        2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042]
    if ndf in range(1, 31):
        tval = ttable[ndf-1]
    elif ndf in range(31, 41):
        frac = (ndf - 30.0) / 10.0 #=0.1 to 1
        tval = 2.042 - frac * (2.042 - 2.021)
    elif ndf in range(41, 61):
        frac = (ndf - 40.0) / 20.0 #=0.05 to 1
        tval = 2.021 - frac * (2.021 - 2.000)
    elif ndf in range(61, 121):
        frac = (ndf - 60.0) / 60.0 #=1/60 to 1
        tval = 2.000 - frac * (2.000 - 1.980)
    elif ndf > 120:
        tval = 1.96
    else:
        print ' ERROR IN TVALUE '
    return tval           

def lik_intervals(names, pars, SD, notfixed, m, func, dataset):
    
    Lmax = -SSDlik(pars, func, dataset, notfixed)
    clim = sqrt(2. * m)
    Lcrit = Lmax - m
    Llimits = []
    
    i = 0
    for j in range(len(pars)):
        if notfixed[j]:
            print('\nCalculating Lik limits for parameter- {0} = {1:.3f}'.format(names[j], pars[j]))
            xhigh1 = pars[j]
            #TODO: if parameter constrained to be positive- ensure that xlow is positive
            xlow1 = pars[j] - 2 * clim * SD[i]
            xlow2 = xhigh1
            xhigh2 = pars[j] + 5 * clim * SD[i]
            print('\tInitial guesses for lower limit: {0:.3f} and {1:.3f}'.format(xlow1, xhigh1))
            print('\tInitial guesses for higher limit: {0:.3f} and {1:.3f}'.format(xlow2, xhigh2))

            found = False
            iter = 0
            xlowlim, xhighlim = None, None
            while not found and iter < 100: 
                L, theta = SSDlik_contour(((xlow1 + xhigh1) / 2), j, pars,
                    notfixed, func, dataset ) 
                if fabs(Lcrit - L) > 0.01:
                    if L < Lcrit:
                        xlow1 = (xlow1 + xhigh1) / 2
                    else:
                        xhigh1 = (xlow1 + xhigh1) / 2
                else:
                    found = True
                    xlowlim = (xlow1 + xhigh1) / 2
                    print 'lower limit found: ', xlowlim
                iter += 1

            found = False
            iter = 0   
            while not found and iter < 100: 
                #L1, L2, L3 = SSDlik_bisect(xlow2, xhigh2, j, pars, notfixed, hill_equation, dataset)
                L, theta = SSDlik_contour(((xlow2 + xhigh2) / 2), j, pars,
                    notfixed, func, dataset ) 
                if fabs(Lcrit - L) > 0.01:
                    if L > Lcrit:
                        xlow2 = (xlow2 + xhigh2) / 2
                    else:
                        xhigh2 = (xlow2 + xhigh2) / 2
                else:
                    found = True
                    xhighlim = (xlow2 + xhigh2) / 2
                    print 'higher limit found: ', xhighlim
                iter += 1


            Llimits.append([xlowlim, xhighlim])
            i += 1

    return Llimits