import os
import sys
from math import sqrt, fabs
from scipy import optimize
import numpy as np

import cvfit
from cvfit import data
from cvfit import errors
from cvfit.errors import residuals
from cvfit.errors import SSD
from cvfit.errors import SSDlik

class SingleFitSession(object):
    def __init__(self, dataset, equation, output=sys.stdout):
        """
        """
        self.output = output
        self.data = dataset
        self.eq = equation
        self.eq.propose_guesses(self.data)
        self.output.write('Fitting session for ' + self.data.title + ' initialised!')
        
    def fit(self):
        # Least square fitting
        #coeffs, cov, dict, mesg, ier = optimize.leastsq(residuals, self.eq.theta,
        #    args=(self.eq, self.data.X, self.data.Y, self.data.W), full_output=1)
        coeffs, smin = simplex(SSD, self.eq.theta, self.eq, self.data.X, self.data.Y, self.data.W)
        self.eq.theta = coeffs
       
    def calculate_errors(self):

        self.Smin = SSD(self.eq.theta, (self.eq,
            self.data.X, self.data.Y, self.data.W))
        #print '\n SSD \n', Smin
        #hes = errors.hessian(coeffs, eq, set)
        #print '\n Observed information matrix = \n', hes
        self.covar = errors.covariance_matrix(self.eq.theta, self.eq, self.data)
        #print '\n Covariance matrix = \n', covar
        self.correl = errors.correlation_matrix(self.covar)
        self.aproxSD = errors.approximateSD(self.eq.theta, self.eq, self.data)
        self.CVs = 100.0 * self.aproxSD / self.eq.theta
        self.kfit = len(np.nonzero(np.invert(self.eq.fixed))[0])
        self.ndf = self.data.size() - self.kfit
        self.var = self.Smin / self.ndf
        self.Sres, self.Lmax = sqrt(self.var), -SSDlik(self.eq.theta, self.eq, self.data)

        tval = errors.tvalue(self.ndf)
        self.m = tval * tval / 2.0
        self.clim = sqrt(2. * self.m)
        self.Lcrit = self.Lmax - self.m
        self.Llimits = errors.lik_intervals(self.eq.theta, self.aproxSD, self.m, self.eq, self.data)
        
        #self.output.write(self.string_estimates())
        #self.output.write(self.string_liklimits())
        
    def string_estimates(self):
        j = 0
        str = 'Number of point fitted = {0:d}'.format(self.data.size())
        
        str += '\nNumber of parameters estimated = {0:d}'.format(self.kfit)
        str += '\nDegrees of freedom = {0:d}'.format(self.ndf)
        
        str += ('\nResidual error SD = {0:.3f}      (variance = {1:.3f})'.
            format(self.Sres, self.var))
        
        for i in range(len(self.eq.names)):
            str += '\nParameter {0:d}: {1}  \t= {2:.6g}  \t'.format(i+1, self.eq.names[i], self.eq.pars[i])
            if not self.eq.fixed[i]:
                str += '  Approx SD = {0:.6g}\t'.format(self.aproxSD[j])
                str += '  CV = {0:.1f}'.format(self.CVs[j])
                j += 1
            else:
                str += '  (fixed)'

        str += ('\nMinimum SSD = {0:.3f}; \nMax log-likelihood = {1:.3f}'.
            format(self.Smin, self.Lmax))
        str += ('\nCorrelation matrix = ' + 
            '[!!!! PRINTOUT OF CORRELATION MATRIX NOT IMPLEMENTED YET. SORRY.\n')
#        self.output.write(correl)
        if np.any(np.absolute(self.correl - np.identity(self.kfit)) > 0.9):
            str += ("\nWARNING: SOME PARAMETERS ARE STRONGLY CORRELATED (coeff > 0.9); try different guesses")

        if np.any(self.CVs > 33):
            str += "\nWARNING: SOME PARAMETERS POORLY DEFINED (CV > 33%); try different guesses"
        return str

        
    def string_liklimits(self):
        j = 0
        str = '\nLIKELIHOOD INTERVALS\n'
        str += ('{0:.3g}-unit Likelihood Intervals'.format(self.m) +
            '  (equivalent SD for Gaussian- {0:.3g})'.format(self.clim))
        str += '\nLmax= {0:.6g};   Lcrit= {1:.6g}'.format(self.Lmax, self.Lcrit)
        for i in range(len(self.eq.names)):
            str += '\nParameter {0:d}:   {1}\t= {2:.6g}'.format(i+1, self.eq.names[i], self.eq.pars[i])
            if not self.eq.fixed[i]:
                try:
                    str += '\t  LOWER = {0:.6g}'.format(self.Llimits[j][0])
                except:
                    str += '\t  LOWER limit not found'
                try:
                    str += '\t  UPPER = {0:.6g}'.format(self.Llimits[j][1])
                except:
                    str += '\t  UPPER limit not found'
                j += 1
            else:
                str += '\t  (fixed)'
        return str


def load_data(example=False):

    if example:
        filename = (os.path.dirname(os.path.dirname(cvfit.__file__)) +
            "/Example/Example.csv")
    else:
        filename = data.ask_for_file()
    try:
        #allsets = cfio.read_sets_from_csv(filename, col=2)
        allsets = data.read_sets_from_csv(filename, 'csv', col=2, header=0, namesin=False, weight=1)
    except ValueError:
        print('fitting.py: WARNING: Oops! File did not load properly...')
    return allsets, filename

def check_input(text, accept, default):
    '''
    Check if the input is in acceptable range or not
    If not, ask to key in another value
    '''

    inputnumber = raw_input(text)
    if inputnumber:
        while inputnumber not in accept:
            print text
            inputnumber = raw_input(text)
    else:
        inputnumber = default

    return int(inputnumber)

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
        
    weightmode = check_input('Mode number [1]: ', weightingmodes, 1)
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
    ieq = check_input('Choose equation to fit [1] : ', ['1', '2'], 1)
    if ieq == 2:
        eq = 'Langmuir'
    return eq

def simplex(func, theta, *args):

    """
    Python implementation of DC's SIMPLEXV.FOR subroutine.
    """

    verbose = 0

    print '\n USING FAST VERSION OF SIMPLEX'

    # these might come as parameters
    errfac = 1.e-6   #1.e-4
    stpfac= 0.1   #0.1
    reffac = 1    # reflection coeff => 1
    confac = 0.5     # contraction coeff 0 < beta < 1
    extfac = 2     # extension factor > 1
    locfac = 5    # 1

    k = np.size(theta)
    n = k + 1    # # of vertices in simplex
    simp = np.zeros((n, k))
    fval = np.zeros((n))
    step = np.zeros((k))
    crtstp = np.zeros((k))
    pnew = np.zeros((k))
    pnew1 = np.zeros((k))
    thmin = np.zeros((k))

    for j in range(k):
        step[j] = stpfac * theta[j]
	crtstp[j] = errfac * theta[j]

    neval = 0	 # counts function evaluations
    nrestart = 5    # max number of restarts
    irestart = 0    # counts restarts

    while irestart < nrestart:    # 2001	continue	!return here for restart

        irestart += 1
        print 'RESTART#', irestart

        simp[0] = theta
        fval[0] = func(theta, args)
        fsav = fval[0]
        absmin = fval[0]
        thmin = theta
        neval += 1

        # compute offset of the vertices of the starting simplex
	fac = (sqrt(float(n)) - 1.) / (float(k) * sqrt(2.))

        #specify all other vertices of the starting simplex
        for i in range(1, n):
            for j in range(k):
                simp[i,j] = simp[0,j] + step[j] * fac
            simp[i, i-1] = simp[0, i-1] + step[i-1]*(fac+1./sqrt(2.))

            #  and calculate their residuals
            fval[i] = func(simp[i], args)
        neval += k

        if verbose: print '\n simplex at the beginning of restart='
        if verbose: print simp
        if verbose: print ' fval at the begining', fval

        fval, simp = sortShell(fval, simp)
        if fval[0] < absmin:
            absmin = fval[0]
            thmin = simp[0]

        # start iteration loop here at 2000
        L = 0
        niter	= 0
        while L == 0:
            niter += 1
            print 'iter#', niter, 'f=', fval[0], 'theta', simp[0]

            # ----- compute centroid of all vertices except the worst
            centre = np.zeros((k))
            for i in range(k):
                for j in range(n-1):
                    centre[i] = centre[i] + simp[j,i]

            # ----- reflect, with next vertex taken as reflection of worst
            for j in range(k):
                centre[j] = centre[j] / float(k)
                pnew[j] = centre[j] - reffac * (simp[-1,j]-centre[j])
            fnew = func(pnew, args)
            if fnew < absmin:
                absmin = fnew
                thmin = pnew
            neval += 1
            if verbose: print 'reflection: e#', neval, 'f=',fnew, 'pnew=', pnew

            if fnew < fval[0]:
                # ----- new vertex is better than previous best so extend it
                for j in range(k):
                    pnew1[j] = centre[j] + extfac * (pnew[j] - centre[j])
                fnew1 = func(pnew1, args)
                if fnew1 < absmin:
                    absmin = fnew1
                    thmin = pnew1
                neval += 1
                if verbose: print 'extention: e#', neval, 'f1=',fnew1, 'pnew1=', pnew1

                if fnew1 < fnew:     # ----- still better
                    simp[-1] = pnew1
                    fval[-1] = fnew1
                else:
                    simp[-1] = pnew
                    fval[-1] = fnew
                fval, simp = sortShell(fval, simp)
                # goto 1901 for convergence check

            else:     # 112  come here if reflected vertex not
                      # better than best vertex, so no extension wanted

                if fnew < fval[-2]:
                    simp[-1] = pnew
                    fval[-1] = fnew
                else:
                    if fnew < fval[-1]:
                        simp[-1] = pnew
                        fval[-1] = fnew
                    # Contract on the original fval(IHI) side of the centroid
                    for j in range(k):
                        pnew1[j] = centre[j] + confac * (simp[-1, j] - centre[j])
                    fnew1 = func(pnew1, args)
                    if fnew1 < absmin:
                        absmin = fnew1
                        thmin = pnew1
                    neval += 1
                    if verbose: print 'contraction: e#', neval, 'f=',fnew1, 'pnew1=', pnew1

                    # ----- is contracted vertex better than the worst vertex
                    if fnew1 <= fval[-1]:
                        simp[-1] = pnew1
                        fval[-1] = fnew1
                    else:
                        #  ----- no, it is still bad, shrink whole simplex towards best vertex
                        for i in range(n):
                            for j in range(k):
                                if j != i:
                                    simp[i,j] = simp[0,j] + confac * (simp[i,j] - simp[0,j])
                            fval[i] = func(simp[i], args)
                            neval += 1
                            if verbose: print 'reduction: e#', neval, 'f=',fval[i], 'theta=', theta

            fval, simp = sortShell(fval, simp)
            if fval[0] < absmin:
                absmin = fval[0]
                thmin = simp[0]


            # CHECK CONVERGENCE. IF NOT CONVERGED GOTO 2000.
            # This version uses diff between highest and lowest value
            # of parameter of the n values that define a vertex
            # (as in O'Neill version)

            #  ----- order the vertices for all vertices
            # Define L=0 for not converged- do next iteration
            # L=1 for converged via crtstp
            # L=2 for converged via delmin (no restarts)
            # L=3 for abort (no restarts)

            L = 1    #  conv via crtstp
            for j in range(k):     # test each parameter
                if(simp[-1,j] - simp[0,j]) > fabs(crtstp[j]): L = 0    # not conv
        # end of iteration (while L == 0:)

        # ----- convergence attained. Options for ending in this version are:
        # 	(1)look at current best vertex
        # 	(2)look at param values averaged over vertices
        # 	(3)look at absmin,thmin. If better, restart at absmin, as below.
        # 	(4)do local search with each param +/- crtstp, as in O'Neill
        # 	 version, starting at current best vertex. If none are better
        # 	 input current best vertex. If some better restart at better
        # 	 value with crtstp taken as approptiately small initial step.

        if L == 1:
            exvals = []
            exvals.append(fval[0])

            # next average over vertices-put values in pnew()
            for j in range(k):
                pnew[j] = 0.0
                for i in range(n):
                    pnew[j] = pnew[j] + simp[i,j]
                pnew[j] = pnew[j] / float(n)
            fvalav = func(pnew, args)
            exvals.append(fvalav)

            exvals.append(absmin)

            # do local search. Put altered values in pnew1
            for j in range(k):
                pnew1[j] = 0.0
                pnew1[j] = simp[0,j] + locfac * crtstp[j]
            fval1 = func(pnew1, args)
            if fval1 < fval[0]:
                exvals.append(fval1)
            else:
                for j in range(k):
                    pnew1[j] = simp[0,j] - locfac * crtstp[j]     # step in other direction
                fval1 = func(pnew1, args)
                exvals.append(fval1)

            # Test which is best.
            il = 0
            for i in range(1, 4):
                if exvals[i] < exvals[il]: il = i
            if il == 0:
                if irestart == nrestart or fsav == fval[0]:
                    print '\n Returned with best vertex'
                    return simp[0], fval[0]
                else:
                    L = 0
                    theta = simp[0]
                    print '\n Restarted at best vertex'
            elif il == 1:
                if irestart == nrestart or fsav == fvalav:
                    print '\n Returned with averaged vertices'
                    return pnew, fvalav
                else:
                    L = 0
                    theta = pnew
                    print '\n Restarted at averaged vertices'
            elif il == 2:
                if irestart == nrestart or fsav == absmin:
                    print '\n Returned with absolut minimum'
                    return thmin, absmin
                else:
                    L = 0
                    theta = thmin
                    print '\n Restarted at absolut minimum'
            else:
                if irestart == nrestart or fsav == fval1:
                    print '\n Returned with result of local search minimum'
                    return pnew1, fval1
                else:
                    L = 0
                    theta = pnew1
                    print '\n Restarted at result of local search minimum'



def sortShell(vals, simp):
    """
    Shell sort using Shell's (original) gap sequence: n/2, n/4, ..., 1.
    """
    n = np.size(vals)
    gap = n // 2
    while gap > 0:
         # do the insertion sort
         for i in range(gap, n):
             val = vals[i]
             tsimp = np.zeros((n-1))
             for l in range(0, n-1):
                 tsimp[l] = simp[i,l]
             j = i
             while j >= gap and vals[j - gap] > val:
                 vals[j] = vals[j - gap]
                 simp[j] = simp[j - gap]
                 j -= gap
             vals[j] = val
             simp[j] = tsimp
         gap //= 2
    return vals, simp