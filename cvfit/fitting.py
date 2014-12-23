import os
from math import sqrt
from scipy import optimize
import numpy as np

import cvfit
from cvfit import cfio
from cvfit import errors
from cvfit.errors import residuals, SSD, SSDlik

class SingleFitSession(object):
    def __init__(self, dataset, equation):
        """
        """
        self.data = dataset
        self.eq = equation
        
        self.eq.propose_guesses(self.data)
        print('\nFitting session for ' + self.data.title + ' initialised!')
        
    def fit(self):
        #theta = equation.get_theta()
        # Least square fitting
        coeffs, cov, dict, mesg, ier = optimize.leastsq(residuals, self.eq.theta,
            args=(self.eq, self.data.X, self.data.Y, self.data.W), full_output=1)
        self.eq.theta = coeffs
        print 'coeffs=', coeffs
        print 'theta=', self.eq.theta
        
        print '\n\n\t' + self.data.title + ' fit finished'
        print 'Number of fuction evaluation =', dict['nfev']
        #print 'mesg=', mesg
        #print 'ier=', ier
        #print 'cov=', cov
        #return coeffs
        
    def calculate_errors(self):

        Smin = SSD(self.eq.theta, self.eq, 
            self.data.X, self.data.Y, self.data.W)
        print '\n SSD \n', Smin
        #hes = errors.hessian(coeffs, eq, set)
        #print '\n Observed information matrix = \n', hes
        covar = errors.covariance_matrix(self.eq.theta, self.eq, self.data)
        #print '\n Covariance matrix = \n', covar
        correl = errors.correlation_matrix(covar)

        aproxSD = errors.approximateSD(self.eq.theta, self.eq, self.data)
        CVs = 100.0 * aproxSD / self.eq.theta

        print('Number of point fitted = {0:d}'.format(self.data.size()))
        kfit = len(np.nonzero(np.invert(self.eq.fixed))[0])
        print('Number of parameters estimated = {0:d}'.format(kfit))
        ndf = self.data.size() - kfit
        print('Degrees of freedom = {0:d}'.format(ndf))
        print cfio.string_estimates(self.eq, aproxSD, CVs)

        var = Smin / ndf
        Sres, Lmax = sqrt(var), -SSDlik(self.eq.theta, self.eq, self.data)
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
        Llimits = errors.lik_intervals(self.eq.theta, aproxSD, m, self.eq, self.data)
        print cfio.string_liklimits(self.eq, Llimits)


def load_data(example=False):

    if example:
        filename = (os.path.dirname(os.path.dirname(cvfit.__file__)) +
            "/Example/Example.csv")
    else:
        filename = cfio.ask_for_file()
    try:
        allsets = cfio.read_sets_from_csv(filename, col=2)
    except ValueError:
        print('fitting.py: WARNING: Oops! File did not load properly...')
    return allsets, filename

def set_weights(sets):
    """ 
    Choose weighting method. 
    """
    
    weightingmodes = ['1'] #, '5']
    mode2 = True
    mode4 = True
    for set in sets:
        if set.S.any() == 0:
            mode2 = False
        for i in np.unique(set.X):
            if len(np.where(set.X == i)[0]) == 1:
                mode4 = False
    
    if mode2:
        weightingmodes.append('2')
        weightingmodes.append('3')
    if mode4:
        weightingmodes.append('4')

    print('\nPlease select the weighting method now:')
    print '1: Weights constant; errors from residuals (Default).'
    if mode2:
        print '2: Weights from specified s(Y); errors from weights.'
        print('3: Weights from specified n, the number of values in the' +
            ' average; errors from weights.')
    else:
        print ('2, 3: s(Y) or n are not specified for some or all pints.' + 
            ' Weights cannot by specified from s(Y) or n.')
    if mode4:
        print('4: Weights from s(Y) calculated from Y repeats at the same X;' +
            ' errors from weights.')
    else:
        print ('4: s(Y) cannot be calculated because some or all X have only' +
            ' one repeat. Weights cannot by specified from s(Y).')
    print '5: Arbitrary weights entered by hand now (NOT IMPLEMENTED YET).'
        
    weightmode = cfio.check_input('Mode number [1]: ', weightingmodes, 1)
    for each in sets:
        each.weightmode = weightmode
    return sets

def general_settings():

    general_settings = {}
    print '\nPlease, choose between:'
    print '0- fit all sets with the same equation;'
    print '1- fit each set separately [Default].'
    #fit_separate = cfio.check_input('0 or 1: ', ['0', '1'], 1)
    general_settings['fit_separate'] = 1
    
    print 'Do you want to select fit settings separately?'
    print '0- use same settings for all datasets (Default);'
    print '1- set settings for each dataset separately.'
    #fit_setting = cfio.check_input('0 or 1: ', ['0', '1'], 0)
    general_settings['same_settings'] = 0

    return general_settings
    

def choose_equation():
    print '\nAvailable equations:'
    print '1. Hill equation'
    print '2. Langmuir equation'
    eq = 'Hill'
    ieq = cfio.check_input('Choose equation to fit [1] : ', ['1', '2'], 1)
    if ieq == 2:
        eq = 'Langmuir'
    return eq

    