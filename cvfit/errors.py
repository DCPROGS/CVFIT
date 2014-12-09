import copy
from math import fabs, sqrt, log, pi
from scipy import optimize
import numpy as np
from numpy import linalg as nplin

def residuals(pars, func, X, Y, W):
    '''
    Calculate the weighted residuals.
    '''
    return W * (Y - func.to_fit(pars, X))

def SSD(pars, func, X, Y, W):
    """
    Calculate sum of squared deviations.
    """
    return np.sum(W * (Y - func.to_fit(pars, X))**2)

def SSDlik(theta, func, set):
    """
    Calculate likelihood coresponding to the sum of the squared deviations 
    assuming that errors follow Gausian distribution.
    """
    S = SSD(theta, func, set.X, set.Y, set.W)
    Sres = sqrt(S / (set.size() - len(func.fixed)))
    return set.size() * log(sqrt(2 * pi) * Sres) + S / (2 * Sres**2)

def SSDlik_contour(x, num, theta, func, set):
    functemp = copy.deepcopy(func)
    functemp.fixed[num] = True
    functemp.pars[num] = x
    theta = functemp.get_theta()
    result = optimize.minimize(SSDlik, theta, args=(functemp, set), 
        method='Nelder-Mead', jac=None, hess=None)
    return -result.fun

def optimise_deltas(theta, func, set):
    """ """
    Lmax = -0.5 * SSD(theta, func, set.X, set.Y, set.W)
    Lcrit = 1.005 * Lmax
    deltas = 0.1 * theta
    L = -0.5 * SSD(theta + deltas, func, set.X, set.Y, set.W)
    if L > Lcrit:
        count = 0
        while L > Lcrit and count < 100:
            deltas *= 2
            L = -0.5 * SSD(theta + deltas, func, set.X, set.Y, set.W)
            count += 1
    elif L < Lcrit:
        count = 0
        while L < Lcrit and count < 100:
            deltas *= 0.5
            L = -0.5 * SSD(theta + deltas, func, set.X, set.Y, set.W)
            count += 1
    return deltas

def hessian(theta, func, set):
    """
    """
    kfit = theta.size
    hessian = np.zeros((kfit, kfit))
    deltas = optimise_deltas(theta, func, set)
    i = 0
    for i1 in range(kfit):
        j = 0
        for j1 in range(kfit):
            coe1, coe2, coe3, coe4 = theta.copy(), theta.copy(), theta.copy(), theta.copy()

            if i == j:
                coe1[j1] += deltas[j1]
                coe3[j1] -= deltas[j1]
                hessian[i, j] = ((
                    SSD(coe1, func, set.X, set.Y, set.W) -
                    2.0 * SSD(theta, func, set.X, set.Y, set.W) +
                    SSD(coe3, func, set.X, set.Y, set.W)) /
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
                    SSD(coe1, func, set.X, set.Y, set.W) -
                    SSD(coe2, func, set.X, set.Y, set.W) -
                    SSD(coe3, func, set.X, set.Y, set.W) +
                    SSD(coe4, func, set.X, set.Y, set.W)) /
                    (4 * deltas[i1] * deltas[j1]))
            j += 1
        i += 1
    return 0.5 * hessian

def covariance_matrix(theta, func, set):
    """ """
    cov = nplin.inv(hessian(theta, func, set))
    minssd = SSD(theta, func, set.X, set.Y, set.W)
    kfit = theta.size
    #TODO: check if weights were used to calculate SSD. If not then
    # calculate errvar.
    errvar = minssd / (set.size() - kfit)
    return cov * errvar

def approximateSD(theta, func, set):
    """ """
    cov = covariance_matrix(theta, func, set)
    return np.sqrt(cov.diagonal())
    
def correlation_matrix(covar):
    correl = np.zeros((len(covar),len(covar)))
    for i1 in range(len(covar)):
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

def lik_intervals(theta, SD, m, func, set):
    
    Lmax = -SSDlik(theta, func, set)
    clim = sqrt(2. * m)
    Lcrit = Lmax - m
    Llimits = []
    
    i = 0
    for j in range(len(func.pars)):
        if not func.fixed[j]:
            #print('\nCalculating Lik limits for parameter- {0} = {1:.3f}'.format(func.names[j], theta[i]))
            xhigh1 = theta[i]
            #TODO: if parameter constrained to be positive- ensure that xlow is positive
            xlow1 = theta[i] - 2 * clim * SD[i]
            xlow2 = xhigh1
            xhigh2 = theta[i] + 5 * clim * SD[i]
            #print('\tInitial guesses for lower limit: {0:.3f} and {1:.3f}'.format(xlow1, xhigh1))
            #print('\tInitial guesses for higher limit: {0:.3f} and {1:.3f}'.format(xlow2, xhigh2))

            found = False
            iter = 0
            xlowlim, xhighlim = None, None
            while not found and iter < 100: 
                L = SSDlik_contour(((xlow1 + xhigh1) / 2), j, theta,
                    func, set) 
                if fabs(Lcrit - L) > 0.01:
                    if L < Lcrit:
                        xlow1 = (xlow1 + xhigh1) / 2
                    else:
                        xhigh1 = (xlow1 + xhigh1) / 2
                else:
                    found = True
                    xlowlim = (xlow1 + xhigh1) / 2
                    #print 'lower limit found: ', xlowlim
                iter += 1
            found = False
            iter = 0   
            while not found and iter < 100: 
                #L1, L2, L3 = SSDlik_bisect(xlow2, xhigh2, j, theta, notfixed, hill_equation, dataset)
                L = SSDlik_contour(((xlow2 + xhigh2) / 2), j, theta,
                    func, set) 
                if fabs(Lcrit - L) > 0.01:
                    if L > Lcrit:
                        xlow2 = (xlow2 + xhigh2) / 2
                    else:
                        xhigh2 = (xlow2 + xhigh2) / 2
                else:
                    found = True
                    xhighlim = (xlow2 + xhigh2) / 2
                    #print 'higher limit found: ', xhighlim
                iter += 1
            Llimits.append([xlowlim, xhighlim])
            i += 1
    return Llimits