from __future__ import division
import sys
import os
import codecs


import numpy as np
from scipy import optimize, stats
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import markdown
import prettyplotlib as ppl
import algopy
from collections import defaultdict
import string
import copy


from toolkit import Report, check_input


class Cell(object):

    '''
    A class that store the inportant information for each cell
    '''

    def __init__(self, concentration, response):
        '''
        Make the class by specifying the concentration and response.
        Original data is then sorted according to the concentration.
        '''

        # Keep track of the original data
        self.originalconcentration = concentration
        self.originalresponse = response
        # Sort the data according to concentration
        result = np.vstack((concentration, response))
        result.sort()
        self.concentration = result[0]
        self.response = result[1]
        self.size = len(result[0])


class HillEquation(Cell):

    '''
    A class for fitting the data with Hill equation.
    '''

    def __init__(self, cell):
        Cell.__init__(self, cell.concentration, cell.response)

    def set_condition_hillequation(self, index, rawresult):
        '''
        Setting the conditions for fitting.
        '''

    # Select weighting method
        self.weightmode = None
        for i in np.unique(self.concentration):
            if len(np.where(self.concentration == i)[0]) == 1:
                self.weightmode = 1 # here len is called to check the number of entries for a particular concentration
                print 'Concentration', i, 'has only one value.'
                print 'Using equal weighted method.'
                break

        if self.weightmode is None:
            print 'Please select the weight calculation option.'
            print 'Mode 1: Equally weighted.'
            print 'Mode 2: Weight by the standard error.(Default)'
            self.weightmode = check_input('Mode number:', ['1', '2'], 2)

        # Select the number of Components
        print 'Please key in the number of the Components. Default = 1'
        self.Component = check_input('Number of the Components:',
                                     ['1', '2'], 1)

        # Detect the trend
        # Selecting whether the trend is positive or negative
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            self.concentration, self.response)
        if slope > 0:
            trend = 'positive'
            self.trend = 1
        elif slope < 0:
            trend = 'negative'
            self.trend = -1
        print 'Select the slope.', '(default:', trend, ')'
        print 'Please key in 1 for positive slope and -1 for negative slope'
        self.trend = check_input('-1 or 1:', ['-1', '1'], self.trend)

        if self.trend == 1:
            # Select to fix Ymin/Y(0) or not
            print 'Please select to fix the Y(0) or not.'
            print '0 for False and 1 for True. (Default = 1)'
            self.YminFixed = check_input('0 or 1:', ['0', '1'], 1)
        else:
            self.YminFixed = 0

        # Set the initial value for fixing Ymax
        self.YmaxFixed = False
        
            
        #fix the index and rawresult
        self.rawresult=rawresult
       # self.index=int(index)
          
        return self.weightmode, self.Component, self.trend, \
        self.YminFixed, self.YmaxFixed
        
    def directset_condition_hillequation(self, weightmode, Component, trend,
                                         YminFixed, YmaxFixed):
        self.weightmode = weightmode
        self.Component = Component
        self.trend = trend
        self.YminFixed = YminFixed
        self.YmaxFixed = YmaxFixed

    def hill_equation_guess(self):
        '''
        Calculate the guess for fitting with hill equation.
        '''

        if self.Component == 1:
            self.guess = np.empty(4)
            # Four parameters will be needed to fit one Component Hill
            # equation

            if self.trend == 1:
                # If response increases with concentration

                # Determine Y(0)
                if self.YminFixed:
                    self.guess[0] = 0
                else:
                    self.guess[0] = np.mean(
                        self.response[self.concentration == self.concentration[0]])
                if self.YmaxFixed:
                    self.guess[1] = 1
                else:
                    # Determine Ymax
                    self.guess[1] = np.mean(
                        self.response[self.concentration == self.concentration[-1]]) \
                        - self.guess[0]
            else:
                # If response decreases with concentration

                # Determine Y(0)
                self.guess[0] = np.mean(
                    self.response[self.concentration == self.concentration[-1]])

                # Determine Ymin
                self.guess[1] = np.mean(self.response[self.concentration == self.concentration[0]]) - \
                    self.guess[0]

            # Determine Kr
            Kr = self.guess[2] = 10 ** ((np.log10(self.concentration[0]) +
                                         np.log10(self.concentration[-1])) / 2)

            # Determine nH
            Ymax = np.amax(self.response) 
            LinRegressCon = self.concentration[self.response < Ymax]
            LinRegressResponse = self.response[self.response < Ymax]
            LinRegressX = np.log10(LinRegressCon) - np.log10(self.guess[2])
            LinRegressY = np.log10(
                (LinRegressResponse / Ymax) / (1 - (LinRegressResponse / Ymax)))
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                LinRegressX, LinRegressY) 
            nH = self.guess[3] = slope

        elif self.Component == 2:
            print 'Two Components fitting is not completed.'
            sys.exit(0)
        
    def fitting_hill_equation(self):
        '''
        Use least square fitting method to fit data with hill equation.
        '''

        # Calculate weight
        if self.weightmode == 1:
            self.weight = np.ones(self.size)
        elif self.weightmode == 2:
            self.weight = np.empty(self.size)
            for index, concentration in enumerate(cell.concentration):
                response = cell.response[
                    cell.concentration == concentration] # this is jsut sorting the list values accordingly
                weight[index] = 1 / \
                    (np.std(response) / np.sqrt(len(response)))
                #print index

        # Least square fitting
        coeffs = []
        params = []
        
        coeffs, params = optimize.leastsq(
            hill_equation_residuals, self.guess,
            args=(self.response, self.weight, self.concentration,
                  self.Component, self.YminFixed, self.YmaxFixed))
        print '\n guess, response, weight, concentration, component, yminfixed, ymaxfixed \n', self.guess,self.response, self.weight, self.concentration, self.Component, self.YminFixed, self.YmaxFixed
        self.coeffs = coeffs
        self.params = params
       # print self.params
       # global index
       # for index, indexing in enumerate()
       # hill_parameter_std (self.coeffs, self.response, self.weight, self.concentration,
       #                    self.Component, self.YminFixed, self.YmaxFixed, self.rawresult, self.index)
   # def coeffparamf(coeffparam):
        
      #  coeffs = []        # initialize array
       # params = []        # initialize array

      #  a, b = coeffparam()  # get values

      #  coeffs.append(a)   # store first value in first array
        #params.append(b)   # store second value in first array
      #  return coeffs, params
         
    def normalise_hill_equation(self):
        '''
        Nomalise the coefficients and response.
        '''

        if self.trend == 1:
            # Nomalise the coefficients by fixing the Y(0) and Ymax
            self.normalised_coeffs = self.coeffs.copy()
            self.normalised_coeffs[0] = 0
            self.normalised_coeffs[1] = 1

            # Nomalise the response
            self.normalised_response = \
                (self.response.copy() -
                 self.coeffs[0]) / self.coeffs[1]
        elif self.trend == -1:
            # Nomalise the coefficients by fixing the Y(0) and Ymax
            self.normalised_coeffs = self.coeffs.copy()
            self.normalised_coeffs[0] = 1
            self.normalised_coeffs[1] = -1

            # Nomalise the response
            self.normalised_response = 1 - \
                (self.response.copy() -
                 self.coeffs[0]) / self.coeffs[1]

    def calculate_average_standarderror(self):
        '''
        Calculate the average response in one concentration
        and the stardard deviation of the mean
        '''

        self.uniqueconcentration = np.array([])
        self.averageresponse = np.array([])
        self.standarderror = np.array([])
        for con in np.unique(self.concentration):
            response = self.response[self.concentration == con]
            self.averageresponse = np.append(
                self.averageresponse, np.mean(response))
            self.uniqueconcentration = np.append(
                self.uniqueconcentration, con)
            self.standarderror = np.append(
                self.standarderror, np.std(response) / np.sqrt(len(response)))
        
    
    def plot_originaldata_guess_fit(self, filename, index):
        '''
        Plot the original data with guess and fitted curve
        in log graph.
        '''

        fig, ax = plt.subplots()
        # Plot original data
        ppl.plot(ax, self.concentration, self.response,
                 'o', label='experimental data')

        logplotX = np.log10(self.concentration)
        plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
                                  np.ceil(np.amax(logplotX) + 1), 100)
        self.xmin = plotX[0].copy()
        self.xmax = plotX[-1].copy()

        # Plot guess curve
        plotYguess = hill_equation(
            self.guess, plotX, self.Component, self.YminFixed, False)
        ppl.plot(ax, plotX, plotYguess, label='guess')

        # Plot fitted curve
        plotYfit = hill_equation(
            self.coeffs, plotX, self.Component, self.YminFixed, False)
        ppl.plot(ax, plotX, plotYfit, label='fit')

        ppl.legend(ax, loc='lower right')
        ax.set_xscale('log')
        ax.set_title(
            os.path.split(filename)[1][:-4] + ' ' + str(index + 1))
        fig.savefig(filename[:-4] + '_' + str(index + 1) + '.png')
        plt.close(fig)
        
#    def hill_parameter_std (coeffs, response, weight, concentration,
#                            Component, YminFixed, YmaxFixed):
#                               
#        result = weight * \
#        (response - hill_equation(coeffs, concentration, Component, YminFixed, YmaxFixed))
#        simplehess(coeffs,YminFixed, YmaxFixed)
#        # print weight
#        # print result
#        #calculate error variance for residuals for estimated covariance matrix
#        #SSmin/n(number of points)-p(number of free parameters)
#        #ervar = result ./ (len(rawresult[: , 2*index])-len(coeffs))
#    
#        # j = np.multiply(index,2)
#        # k = rawresult[:, j]
#        # klen = k.shape[1]
#        klen = len(concentration)
#        m = len(coeffs)
#        l = klen-m
#        ervar=np.divide(result,l)
#        # print ervar
#        return ervar
#        


def hill_equation(coeffs, concentration, Component, YminFixed, YmaxFixed):
    '''
    The hill equation
    '''

    if YminFixed:
        ymin = 0
    #    i=0
    else:
        ymin = coeffs[0]
    #    i=1
        
    if YmaxFixed:
        ymax = 1
     #   j=0
    else:
        ymax = coeffs[1]
   #     j=1

    equation = ymin + ymax * (
        ((concentration / coeffs[2]) ** coeffs[3])
        /
        (1 + (concentration / coeffs[2]) ** coeffs[3]))
    return equation
    #print ('equation=' equation)


def hill_equation_residuals(coeffs, response, weight, concentration,
                            Component, YminFixed, YmaxFixed):
    '''
    Calculate the weighted residuals for the hill equation
    '''

    result = weight * \
        (response -
         hill_equation(coeffs, concentration,
                       Component, YminFixed, YmaxFixed))
   # print ('result', result)
    return result
    
def SS(coeffs, response, weight, concentration, Component, YminFixed, YmaxFixed):
    residuals = (response - hill_equation(coeffs, concentration, Component, YminFixed, YmaxFixed))

    if len(response)==100:
        flename = 'C:/Users/Phillip/Desktop/CVFIT/cell5.csv'
        rwresult = np.genfromtxt(flename, delimiter=',')
        print('Rawresult', rwresult)
        SSres = sum(np.power(rwresult,2))
    else:
        weres = weight * residuals
        SSres = sum(np.power(weres,2))
     #   print('weighted residuals', weres)
   # print('response',response)
   # print('Coefffs', coeffs)
   # print('residuals', residuals)
    
    print('ss', SSres)
    
        #print('sum of squares', SS)
    return SSres
        
def simplehess(coeffs, YminFixed, YmaxFixed, concentration, Component, response, weight):
    
   #length of the free parameter array - Y can be fixed by user
    ccoeffs = copy.deepcopy(coeffs)    
    if YminFixed:
        i = 1
        if YmaxFixed:
            j=1
        else:
            j=0
        cdoeffs = ccoeffs[0+i+j:]
    else:
        i= 0
        if YmaxFixed:
            j=1
            cdoeffs = ccoeffs[0]
            cdoeffs=np.append(cdoeffs, ccoeffs[0+j+j:])
        else:
            j=0
            cdoeffs = ccoeffs[0:]
    
    numfree = len(cdoeffs)
    #creating a list to combine with itself to create a look-up "hessian matrix"
    #with variables with respect to which the 1st and 2nd partial derivatives are taken
    # - Hessindma. Hessmat contains the actual values. Hessindma is declared for clarity.
    #freepar = string.lowercase[:numfree]
    freepar=list(range(0,numfree))
    print ('coeffs', coeffs)    
    print ('cdoeffs', cdoeffs)
    print ('freepar', freepar)
    #hessindmat = defaultdict(freepar,freepar)
    hessindma = [[[1,1] for x in xrange(numfree)] for y in xrange(numfree)]
    hessmat = [[1 for g in xrange(numfree)] for h in xrange(numfree)]
    #construct an array of fractions of coeffs such that each individual deltaH does not change likelihood
    #by more than 0.5%
    hdelta = optdel(coeffs, concentration, Component, YminFixed, YmaxFixed, response, weight, freepar, i, j)
    
    for endex in range(len(freepar)):
        undex = 0
       # print endex
        for undex in range(len(freepar)):
            hessindma[endex][undex][0] = freepar[endex]
            hessindma[endex][undex][1] = freepar[undex]
            hessmat[endex][undex] = ddfunk(coeffs, hessindma, undex, endex, concentration, Component, YminFixed, YmaxFixed, i, j, response, weight, hdelta)
      #  hessindmat(index,undex) = 
    #print ('hessindma', hessindma)
    #print ('hessian', hessmat)
    return hessmat
    
def optdel(coeffs, concentration, Component, YminFixed, YmaxFixed, response, weight, freepar, i, j):
    #here the minus sign is just following mathematical convention - no reason for it within the script
    elmax = 0.5 * SS(coeffs, response, weight, concentration, Component, YminFixed, YmaxFixed)
    #practically the delta H is between 0.025 and 0.05 elmax: following the original script by DC
    #hdelta is ony constructed once, and the yndex here doesn't have anything to do with *ndex anywhere else
    hda = copy.deepcopy(coeffs) #just a placeholder
    hdelta = np.zeros(len(coeffs)) #this is the permanent storage
    elcrit = elmax + 0.005 * elmax
    el = SS(coeffs, response, weight, concentration, Component, YminFixed, YmaxFixed)
    for yndex in range(len(freepar)):
        hdel = np.zeros(len(coeffs)) #this is hdelta for within the loop, its flushed and reused
        hdel[yndex + i + j] = 0.1 * hda[yndex + i + j]
        hcoeffs = hda + hdel
        print('hcoeffs', hcoeffs)
        el = 0.5 * SS(hcoeffs, response, weight, concentration, Component, YminFixed, YmaxFixed)
        while el < elcrit:
            hdel[yndex+i+j] = 2*hdel[yndex+i+j]
            hcoeffs[yndex+i+j] = hda[yndex+i+j] + hdel[yndex+i+j]
            el = 0.5 * SS(hcoeffs, response, weight, concentration, Component, YminFixed, YmaxFixed)
        
        while el > elcrit:
      #  else:
            hdel[yndex+i+j] = 0.5*hdel[yndex+i+j]
            hcoeffs[yndex+i+j] = hda[yndex+i+j] + hdel[yndex+i+j]
            el = 0.5 * SS(hcoeffs, response, weight, concentration, Component, YminFixed, YmaxFixed)
        hdelta[yndex + i + j] = hdel[yndex + i + j]
        print('hdel', hdel)
        print('loop hdelta', hdelta)
    print('el', el)
    print('elcrit', elcrit)
    print('hdelta', hdelta)
    return hdelta
        
#def dfunk(coeffs, dcoeffs, ddcoeffs, concentration, Component, YminFixed, YmaxFixed):
def dfunk(dcoeffs, ddcoeffs, concentration, Component, YminFixed, YmaxFixed, response, weight):               
    #dfx = (hill_equation(ddcoeffs, concentration, Component, YminFixed, YmaxFixed) - hill_equation(dcoeffs, concentration, Component, YminFixed, YmaxFixed))
    dfx = (SS(ddcoeffs, response, weight, concentration, Component, YminFixed, YmaxFixed) - SS(dcoeffs, response, weight, concentration, Component, YminFixed, YmaxFixed))
    return dfx    
    
def ddfunk(coeffs, hessindma, undex, endex, concentration, Component, YminFixed, YmaxFixed, i, j, response, weight, hdelta):
    
#    #here, d is a prefix used to mark the number of times +x/step was applied
#    step = 1000
#    dcoeffs = copy.deepcopy(coeffs)
#    ddcoeffs = copy.deepcopy(coeffs)
#    #creating coeffs arrays of the form [x+x/step...]
#    ddcoeffs[int(hessindma[endex][undex][0]+i+j)] = ddcoeffs[int(hessindma[endex][undex][0]+i+j)]+ddcoeffs[int(hessindma[endex][undex][0]+i+j)]/step
#    ddcoeffs[int(hessindma[endex][undex][1]+i+j)] = ddcoeffs[int(hessindma[endex][undex][1]+i+j)]+ddcoeffs[int(hessindma[endex][undex][1]+i+j)]/step
#    dcoeffs[int(hessindma[endex][undex][1]+i+j)] = dcoeffs[int(hessindma[endex][undex][1]+i+j)]+dcoeffs[int(hessindma[endex][undex][1]+i+j)]/step
#    print ('coeffs', coeffs)    
#    print ('dcoeffs', dcoeffs)
#    print ('ddcoeffs', ddcoeffs)
#   # ddfx = (dfunk(coeffs, dcoeffs, ddcoeffs, concentration, Component, YminFixed, YmaxFixed)-dfunk(coeffs, dcoeffs, ddcoeffs, concentration, Component, YminFixed, YmaxFixed))
#    ddfx = (dfunk(dcoeffs, ddcoeffs, concentration, Component, YminFixed, YmaxFixed, response, weight)-dfunk(coeffs, dcoeffs, concentration, Component, YminFixed, YmaxFixed, response, weight))
#   # ddfx2 = (hill_equation(ddcoeffs, concentration, Component, YminFixed, YmaxFixed)-hill_equation(dcoeffs, concentration, Component, YminFixed, YmaxFixed))
#    return ddfx

    #following symmetric numerical approximation formulae: need four different sets of free parameter
    #values: coeff(SSmin)+-delta H. The formulae are different for diagonal and off-diagonal elements
    #den - denominator
    coe1 = copy.deepcopy(coeffs)
    coe2 = copy.deepcopy(coeffs)
    coe3 = copy.deepcopy(coeffs)
    coe4 = copy.deepcopy(coeffs)
    
    if hessindma[endex][undex][0] == hessindma[endex][undex][1]:
        coe1[i+j+hessindma[endex][undex][0]]=coe1[i+j+hessindma[endex][undex][0]] + hdelta[i+j+hessindma[endex][undex][0]]
        coe3[i+j+hessindma[endex][undex][0]]=coe1[i+j+hessindma[endex][undex][0]] - hdelta[i+j+hessindma[endex][undex][0]]
        den = np.power((hdelta[i + j + hessindma[endex][undex][0]]),2)
        ddfx = -0.5 * (SS(coe1, response, weight, concentration, Component, YminFixed, YmaxFixed) - 2 * SS(coeffs, response, weight, concentration, Component, YminFixed, YmaxFixed) + SS(coe3, response, weight, concentration, Component, YminFixed, YmaxFixed)) / den
    else:
        coe1[i+j+hessindma[endex][undex][0]]=coe1[i+j+hessindma[endex][undex][0]] + hdelta[i+j+hessindma[endex][undex][0]]
        coe1[i+j+hessindma[endex][undex][1]]=coe1[i+j+hessindma[endex][undex][1]] + hdelta[i+j+hessindma[endex][undex][1]]
        coe2[i+j+hessindma[endex][undex][0]]=coe2[i+j+hessindma[endex][undex][0]] + hdelta[i+j+hessindma[endex][undex][0]]
        coe2[i+j+hessindma[endex][undex][1]]=coe2[i+j+hessindma[endex][undex][1]] - hdelta[i+j+hessindma[endex][undex][1]]
        coe3[i+j+hessindma[endex][undex][0]]=coe3[i+j+hessindma[endex][undex][0]] - hdelta[i+j+hessindma[endex][undex][0]]
        coe3[i+j+hessindma[endex][undex][1]]=coe3[i+j+hessindma[endex][undex][1]] + hdelta[i+j+hessindma[endex][undex][1]]
        coe4[i+j+hessindma[endex][undex][0]]=coe4[i+j+hessindma[endex][undex][0]] - hdelta[i+j+hessindma[endex][undex][0]]
        coe4[i+j+hessindma[endex][undex][1]]=coe4[i+j+hessindma[endex][undex][1]] - hdelta[i+j+hessindma[endex][undex][1]]
        den = (4 * hdelta[i+j+hessindma[endex][undex][0]] * hdelta[i+j+hessindma[endex][undex][1]])
        ddfx = -0.5 * (SS(coe1, response, weight, concentration, Component, YminFixed, YmaxFixed) - SS(coe2, response, weight, concentration, Component, YminFixed, YmaxFixed) - SS(coe3, response, weight, concentration, Component, YminFixed, YmaxFixed) + SS(coe4, response, weight, concentration, Component, YminFixed, YmaxFixed)) / den
    return ddfx
#def hill_parameter_std (coeffs, response, weight, concentration,
#                    Component, YminFixed, YmaxFixed):
def hill_parameter_std (coeffs, response, weight, concentration, Component, YminFixed, YmaxFixed):
                       
    #result = weight * \
    #(response - hill_equation(coeffs, concentration, Component, YminFixed, YmaxFixed))
    #simplehess(coeffs,YminFixed, YmaxFixed, concentration, Component, response, weight)
    # print weight
    # print result
    #calculate error variance for residuals for estimated covariance matrix
    #SSmin/n(number of points)-p(number of free parameters)
    #ervar = result ./ (len(rawresult[: , 2*index])-len(coeffs))

    # j = np.multiply(index,2)
    # k = rawresult[:, j]
    # klen = k.shape[1]
    ccoeffs = copy.deepcopy(coeffs)    
    if YminFixed:
        i = 1
        if YmaxFixed:
            j=1
        else:
            j=0
        cdoeffs = ccoeffs[0+i+j:]
    else:
        i= 0
        if YmaxFixed:
            j=1
            cdoeffs = ccoeffs[0]
            cdoeffs=np.append(cdoeffs, ccoeffs[0+j+j:])
        else:
            j=0
            cdoeffs = ccoeffs[0:]
    sqsu = SS(coeffs, response, weight, concentration, Component, YminFixed, YmaxFixed)
    df = (len(concentration)-len(cdoeffs))
    ervar = sqsu / df
    erstd = np.sqrt(ervar) 
    print ('SS', sqsu)
    print ('DF', df)
    print ('Error Variance', ervar)
    print ('Error StD', erstd)
    #print('ervar',ervar)
    hessian = simplehess(coeffs, YminFixed, YmaxFixed, concentration, Component, response, weight)
    print ('Hessian', hessian)    
#    reshes = ervar*np.array(hessian)
#    print ('hessian*ervar', reshes)
    infomat = np.linalg.inv(hessian)
    print('information matrix =', infomat)
    resinfomat = ervar*np.array(infomat)
    print('information matrix by res error variance =', resinfomat)
    covmat = np.linalg.inv(infomat)
    print('Covariance Matrix', covmat)
    rescovmat = np.linalg.inv(covmat)
    rescovmat = ervar * rescovmat
    print('covariance with error', rescovmat)
    
#    klen = len(concentration)
#    m = len(coeffs)
#    l = klen-m
#    ervar=np.divide(result,l)
    # print ervar
    return infomat
    
def fitting_curve_hill_equation(filename, celllist, report, rawresult, index):
    print 'Do you want to select fitting conditions separately?'
    print '0 means using one condition for all cells (Default)'
    print '1 means key in condition differently for each cell.'
    separate_condition = check_input('0 or 1:', ['0', '1'], 0)
    processlist = []
    for index, cell in enumerate(celllist):
        print 'Cell:', str(index + 1)
        # Create HillEquation class
        cell = HillEquation(cell)
        processlist.append(cell)
        # Setting fitting conditions
        if (separate_condition == 0) and (index == 0):
            weightmode, Component, trend, YminFixed, YmaxFixed = cell.set_condition_hillequation(rawresult, index)
        elif separate_condition == 0:
            cell.directset_condition_hillequation(
                weightmode, Component, trend, YminFixed, YmaxFixed)
        elif separate_condition == 1:
            cell.set_condition_hillequation(rawresult, index)

        # Calculate the guess
        cell.hill_equation_guess()

        # Fitting the hill equation
        cell.fitting_hill_equation()
        

        hill_parameter_std(cell.coeffs, cell.response, cell.weight, cell.concentration,
                           cell.Component, cell.YminFixed, cell.YmaxFixed)
       
       # Normalise the result
        cell.normalise_hill_equation()

#        hill_parameter_std(cell.normalised_coeffs, cell.normalised_response, cell.weight, cell.concentration,
#                           cell.Component, cell.YminFixed, cell.YmaxFixed)

        # Plot the original data, guess and fit
        cell.plot_originaldata_guess_fit(filename, index)

        print 'Fitting completed'
        print 'Y(0)', cell.coeffs[0], 'Ymax', cell.coeffs[1],
        print 'Kr', cell.coeffs[2], 'nH', cell.coeffs[3]

    # Print all the coefficients into the report
    report.title('Fitting result:', 1)
    report.title('Coefficients:', 2)
    title = ['Y(0)', 'Ymax', 'Kr', 'nH']
    report.tabletitle(title)
    table = np.vstack([cell.coeffs for cell in processlist])
    report.table(table)

    # Plotting three graph
    # Generate the common x axis
    xmin = np.floor(np.log10(min([cell.xmin for cell in processlist])))
    xmax = np.ceil(np.log10(max([cell.xmax for cell in processlist])))
    plotx = 10 ** np.linspace(xmin, xmax, 100)
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    for index, cell in enumerate(processlist):
        # Plot original data and fitted curve
        ppl.plot(ax1, cell.concentration, cell.response,
                 'o', color=ppl.colors.set2[index])
        ppl.plot(ax1, plotx, hill_equation(cell.coeffs, plotx,
                                           cell.Component, YminFixed, False),
                 color=ppl.colors.set2[index], label=str(index + 1))

        # Plot normalised curve only
        ppl.plot(ax2, plotx, hill_equation(cell.normalised_coeffs, plotx,
                                           cell.Component, YminFixed, False),
                 color=ppl.colors.set2[index], label=str(index + 1))

        # Plot normalised curve with normalised data
        ppl.plot(ax3, cell.concentration, cell.normalised_response,
                 'o', color=ppl.colors.set2[index])
        ppl.plot(ax3, plotx, hill_equation(cell.normalised_coeffs, plotx,
                                           cell.Component, YminFixed, False),
                 color=ppl.colors.set2[index], label=str(index + 1))

    ppl.legend(ax1, loc='lower right')
    ax1.set_xscale('log')
    ax1.set_title(os.path.split(filename)[1][:-4])
    fig1.savefig(filename[:-4] + '_originaldata_fittedcurve' + '.png')
    plt.close(fig1)

    ppl.legend(ax2, loc='lower right')
    ax2.set_xscale('log')
    ax2.set_title(os.path.split(filename)[1][:-4])
    fig2.savefig(filename[:-4] + '_normalisedfittedcurve' + '.png')
    plt.close(fig2)

    ppl.legend(ax3, loc='lower right')
    ax3.set_xscale('log')
    ax3.set_title(os.path.split(filename)[1][:-4])
    fig3.savefig(
        filename[:-4] + '_normaliseddata_fittedcurve' + '.png')
    plt.close(fig3)

    # print three plots to report
    report.title('Fitted curve:', 2)
    htmlfilename = os.path.split(filename)[-1]
    report.title('Original data with fitted curve:', 3)
    report.image(
        htmlfilename[:-4] + '_originaldata_fittedcurve' + '.png')
    report.title('Normalised fitted curve:', 3)
    report.image(htmlfilename[:-4] + '_normalisedfittedcurve' + '.png')
    report.title('Normalised data and fitted curve:', 3)
    report.image(
        htmlfilename[:-4] + '_normaliseddata_fittedcurve' + '.png')

    # Pool all the data together and fit
    totalcell = HillEquation(Cell(
        np.hstack([cell.concentration for cell in processlist]),
        np.hstack([cell.normalised_response for cell in processlist])))
    totalcell.directset_condition_hillequation(
        weightmode, Component, trend, YminFixed, True)
    totalcell.hill_equation_guess()
    totalcell.fitting_hill_equation()
    totalcell.calculate_average_standarderror()

    # Plot the fiting of nomalised curve
    fig, ax = plt.subplots()
    ax.errorbar(totalcell.uniqueconcentration,
                totalcell.averageresponse,
                yerr=totalcell.standarderror,
                fmt='o', color=ppl.colors.set2[0],
                ecolor=ppl.colors.set2[1], label='Errorbar')
    ppl.plot(ax, plotx, hill_equation(totalcell.coeffs, plotx,
                                      totalcell.Component, YminFixed, True),

             color=ppl.colors.set2[2], label='Fitting')
    ppl.legend(ax, loc='lower right')
    ax.set_xscale('log')
    ax.set_title(os.path.split(filename)[1][:-4])
    fig.savefig(filename[:-4] + '_allfitting' + '.png')
    plt.close(fig)

    # print the final result
    report.title('Fitting with all normalised data:', 1)
    report.title('Fitting plot:', 2)
    report.image(htmlfilename[:-4] + '_allfitting' + '.png')
    report.title('Fitting cofficients:', 2)
    report.tabletitle(title)
    report.table(totalcell.coeffs)
    report.title('Fitting accuracy:', 2)
    title = ['Concentration', 'Average', 'Standard error']
    report.tabletitle(title)
    table = np.vstack((totalcell.uniqueconcentration,
                       totalcell.averageresponse,
                       totalcell.standarderror))
    table = np.transpose(table)
    report.table(table)

    report.outputhtml()
