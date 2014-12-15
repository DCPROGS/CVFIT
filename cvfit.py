#! /usr/bin/python

__author__="remislp"
__date__ ="$29-Nov-2014 15:40:02$"

from math import sqrt
from scipy import optimize
import numpy as np

from cvfit import cfio
from cvfit import errors
from cvfit.errors import residuals, SSD, SSDlik
from cvfit.hill import Hill

def load_data():
    #filename = cfio.ask_for_file()
    filename = "Example/Example.csv"
    try:
        allsets = cfio.read_sets_from_csv(filename, col=2)
        print('File {0} loaded'.format(filename))
        print('{0:d} sets found.'.format(len(allsets)))
        
    except ValueError:
        print('Oops! File did not load properly...')
    
    #TODO: Choose data sets to fit if file contains multiple datasets
    sets = allsets

    # Select weighting method
    print('Please select the weighting method now:')
    print '1: Weights constant; errors from residuals (Default).'
    print '2: Weights from specified s(Y); errors from weights.'
    #weightmode = cfio.check_input('Mode number: ', ['1', '2'], 1)
    weightmode = 1
    for each in sets:
        each.weightmethod = weightmode
    
    for i in range(len(allsets)):
            print '\nSet #{0:d}:'.format(i+1)
            print allsets[i]
    
    return sets
    
        
if __name__ == "__main__":
    
    sets = load_data()
    # TODO: here display data
    
    print '\nPlease, choose between:'
    print '0- fit all sets with the same equation;'
    print '1- fit each set separately [Default].'
    #fit_separate = cfio.check_input('0 or 1: ', ['0', '1'], 1)
    fit_separate = 1
    
    print '\nAvailable equations:'
    print '1. Hill equation'
    #ieq = cfio.check_input('Choose equation to fit [1] :', ['1'], 1)
    #if ieq == 1:
    #    fit_hill_equation(filename, celllist, report)
    
    print 'Do you want to select fit settings separately?'
    print '0- use same settings for all datasets (Default);'
    print '1- set settings for each dataset separately.'
    #fit_setting = cfio.check_input('0 or 1: ', ['0', '1'], 0)
    fit_setting = 0
    
    for set in sets:
        fixed = [True, False, False, False]
        hill = Hill()
        hill.fixed = fixed
        guesses = hill.guesses(set)
        theta = guesses[np.nonzero(np.invert(fixed))[0]]
        
        # Least square fitting
        coeffs, cov, dict, mesg, ier = optimize.leastsq(residuals, theta,
            args=(hill, set.X, set.Y, set.W), full_output=1)
        print '\n\n\t' + set.title + ' fit finished'
        #print 'coeffs=', coeffs
        #print 'cov=', cov
        #print '# fuction evaluation =', dict['nfev']
        #print 'mesg=', mesg
        #print 'ier=', ier
        
        Smin = SSD(coeffs, hill, set.X, set.Y, set.W)
        #print '\n SSD \n', Smin
        hes = errors.hessian(coeffs, hill, set)
        #print '\n Observed information matrix = \n', hes
        covar = errors.covariance_matrix(coeffs, hill, set) 
        #print '\n Covariance matrix = \n', covar
        correl = errors.correlation_matrix(covar)
        
        aproxSD = errors.approximateSD(coeffs, hill, set)
        CVs = 100.0 * aproxSD / coeffs
         
        print 'Number of fuction evaluation =', dict['nfev']
        print('Number of point fitted = {0:d}'.format(set.size()))
        kfit = len(np.nonzero(np.invert(fixed))[0])
        print('Number of parameters estimated = {0:d}'.format(kfit))
        ndf = set.size() - kfit
        print('Degrees of freedom = {0:d}'.format(ndf))
        print cfio.string_estimates(hill, aproxSD, CVs)
        
        var = Smin / ndf
        Sres, Lmax = sqrt(var), -SSDlik(coeffs, hill, set)
        print ('\nResidual error SD = {0:.3f} (variance = {1:.3f})'.format(Sres, var))
        print ('Minimum SSD = {0:.3f}; \tMax log-likelihood = {1:.3f}'.format(Smin, Lmax))
        print '\nCorrelation matrix = \n', correl
        if np.any(np.absolute(correl - np.identity(kfit)) > 0.9):
            print("\nWARNING: SOME PARAMETERS ARE STRONGLY CORRELATED (coeff > 0.9); try different guesses")

        tval = errors.tvalue(ndf)
        m = tval * tval / 2.0
        clim = sqrt(2. * m)
        Lcrit = Lmax - m
        print '\nLIKELIHOOD INTERVALS'
        print ('{0:.3g}-unit Likelihood Intervals'.format(m) + 
            ' (equivalent SD for Gaussian- {0:.3g})'.format(clim))
        print 'Lmax= {0:.6g}; Lcrit= {1:.6g}'.format(Lmax, Lcrit)
        Llimits = errors.lik_intervals(coeffs, aproxSD, m, hill, set)
        print cfio.string_liklimits(hill, Llimits)
