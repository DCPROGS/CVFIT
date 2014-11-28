import numpy as np
from numpy import linalg as nplin
from equations import SSD
from scipy import stats
from scipy import optimize
import matplotlib.pyplot as plt
import copy
import math as math

def optimise_deltas(dataset, func, args, notfixed, fixval):
    """ """
    Lmax = -0.5 * SSD(args, func, dataset, notfixed, fixval)
    Lcrit = 1.005 * Lmax
    deltas = np.zeros((len(args),))
    for index in np.nonzero(notfixed):
        deltas[index] += 0.1 * args[index]
    L = -0.5 * SSD(args + deltas, func, dataset, notfixed, fixval)
    if L > Lcrit:
        count = 0
        while L > Lcrit and count < 100:
            for index in np.nonzero(notfixed):
                deltas[index] *= 2
            L = -0.5 * SSD(args + deltas, func, dataset, notfixed, fixval)
            count += 1
    elif L < Lcrit:
        count = 0
        while L < Lcrit and count < 100:
            for index in np.nonzero(notfixed):
                deltas[index] *= 0.5
            L = -0.5 * SSD(args + deltas, func, dataset, notfixed, fixval)
            count += 1
    #print '\n deltas = \n', deltas
    return deltas

def hessian(dataset, func, pars, notfixed, fixval):
    """
    """
    kfit = np.nonzero(notfixed)[0].size
    hessian = np.zeros((kfit, kfit))
    deltas = optimise_deltas(dataset, func, pars, notfixed, fixval)#deltas
   # hayyam = np.zeros(50)#delete
   # deltas = copy.deepcopy(ideltas)# delete
   # for jar in range(-50, 50, 2): #delete
        #print '\n jar \n', jar
         #delete
   #     pos = 2 #delete
   #     deltas[pos] = copy.deepcopy(ideltas[pos])*(10**(jar/10)) #delete
    i = 0
    #print '\n deltas \n', deltas
    for i1 in np.nonzero(notfixed)[0]:
        j = 0
        for j1 in np.nonzero(notfixed)[0]:
            coe1, coe2, coe3, coe4 = pars.copy(), pars.copy(), pars.copy(), pars.copy()

            if i == j:
                coe1[j1] += deltas[j1]
                coe3[j1] -= deltas[j1]
                hessian[i, j] = ((
                    SSD(coe1, func, dataset, notfixed, fixval) -
                    2.0 * SSD(pars, func, dataset, notfixed, fixval) +
                    SSD(coe3, func, dataset, notfixed, fixval)) /
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
                    SSD(coe1, func, dataset, notfixed, fixval) -
                    SSD(coe2, func, dataset, notfixed, fixval) -
                    SSD(coe3, func, dataset, notfixed, fixval) +
                    SSD(coe4, func, dataset, notfixed, fixval)) /
                    (4 * deltas[i1] * deltas[j1]))
            j += 1
        i += 1
      #  hayyam[(jar+50)/2] = copy.deepcopy(hessian[pos-1, pos-1]) #delete
 #   print '\n deltaH SS = \n', hayyam
        #plt.plot(jar, hayyam)#delete
        #plt.show
    return 0.5 * hessian

def covariance_matrix(dataset, func, pars, notfixed, fixval):
    """ """
    cov = nplin.inv(hessian(dataset, func, pars, notfixed, fixval))
    minssd = SSD(pars, func, dataset, notfixed, fixval)
    #print '\n SSD \n', minssd
    kfit = np.nonzero(notfixed)[0].size
    #TODO: check if weights were used to calculate SSD. If not then
    # calculate errvar.
    errvar = minssd / (dataset.size() - kfit)
    return cov * errvar

def approximateSD(dataset, func, pars, notfixed, fixval):
    """ """
    cov = covariance_matrix(dataset, func, pars, notfixed, fixval)
    return np.sqrt(cov.diagonal())
    
def correlation_matrix(covar):
#    for i1 in range(len(covar)):
#        j1=0
#        for j1 in range(len(covar)):
#            if i1==j1:
#                den[i1] = copy.deepcopy(covar[i1,j1])
#                elif i1>j1:
#                    nom [i1,j1] = covar [i1,j1]
#                else:
#                    pass
    correl = np.zeros((len(covar),len(covar)))
    
    for i1 in range(len(covar)):
        j1=0
        for j1 in range(len(covar)):
            correl[i1,j1] = covar[i1,j1]/np.sqrt(np.multiply(covar[i1,i1],covar[j1,j1]))
    return correl
            


    

    
def new_secs3(aSD, covar, dataset, func, pars, notfixed, hes):
    #bear in mind that "true" likelihood distribution should be inverted.Here it isnt
    fixval = copy.deepcopy(pars) #values which the pars are to be fixed at if true
    fixfix = copy.deepcopy(notfixed)
    
    minssd = SSD(pars, func, dataset, notfixed, fixval)
    kfit = np.nonzero(notfixed)[0].size
    #TODO: check if weights were used to calculate SSD. If not then
    # calculate errvar.
    errvar = minssd / (dataset.size() - kfit)
    pi = math.pi
    slik = kfit * np.log(np.sqrt(2*pi)*errvar)
    ELMAX = slik+(minssd/(2*np.sqrt(errvar)**2))
    print '\n errvar, slik, ELMAX \n', errvar, slik, ELMAX
    df = len(covar)
    t = stats.t.ppf(0.95, df)
    aaSD = approximateSD(dataset, func, pars, notfixed, fixval)
    j = len(pars) - len(aaSD) #correcting for discrepancy between parameter array and freepar number
    i = 2 #parameter which the intervals are sought for
    
    Lmaxm = 0.5*SSD(pars, func, dataset, notfixed, fixval)/errvar
    Lmaxn = copy.deepcopy(Lmaxm) #just a deep copy for robustness during testing  
    m = (t**2)*0.5
    print '\n m \n', m
    
    contour = Lmaxn-m
    
    rile = [1,-1] # either go RIght from Smin or LEft, switches by g in the for loop
    #result placeholders
    RES = np.zeros(4)
    RES1 = np.zeros(4)
    RES2 = np.zeros(4)
    fixfix[i+j] = not (i+j) #this tells that parameter 'i' is sought intervals for, and allows to fix it in minimise
    #this is specific to cell, just needed for the "optimize" routine - to check for consistency
    guess = [  0.00000000e+00,   9.48333333e+03,   1.73205081e+03,   9.97865711e-01]
#    response = [  980.,  3193.,  4806.,  4949.,  6162.,  6183.,  6188.,  6190.,  6355.,  8658., 9000.,  9165.,  9276.,  9284.,  9492.,  9674.]
#    weight = [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.]
#    concentration = [   100.,    300.,    600.,   1000.,   1000.,   1000.,   1000.,   1000.,   1000.,   3000.,  10000.,  10000.,  10000.,  30000.,  30000.,  30000.] 
#

    for g in range (0,1):
        SDel = copy.deepcopy(pars)
        print '\n SDel \n', SDel        
        k = 1
        STEP = float(aaSD[i]*2)/float((2**k))
        lohi = copy.deepcopy(rile[g])
        count = 0
        #LGuess = Lmaxn + 3*SD
        SDel[i+j] += lohi * STEP
        print '\n STEP, lohi, SDel, fixfix \n', STEP, lohi, SDel, fixfix
        
        Lcurall = optimize.minimize(SSD, guess, args=(func, dataset, fixfix, SDel), method='Nelder-Mead', jac=None, hess=None)
        Lcurx = copy.deepcopy(Lcurall.x) #this just takes the optimised parameters from Lcurall
        Lcurx = np.asarray(Lcurx)
        Lcur = SSD(Lcurx, func, dataset, fixfix, Lcurx)
        #print '\n Lcur, Lcurx, contour \n', Lcur, Lcurx, contour
        while np.sqrt(((0.5*Lcur / errvar) - contour)**2) >  0.005*Lmaxn/errvar and count<100:
            if (0.5 * Lcur) - contour > 0:
                n = -1
            else:
                n = 1
            k+=1
           # print '\n Step = \n', STEP
            SDel[i+j] += lohi * n * STEP
            STEP = float(2 * aaSD[i])/float((2**k))
            Lcurx[i+j] += lohi * n * STEP
            Lcurall = optimize.minimize(SSD, guess, args=(func, dataset, fixfix, Lcurx), method='Nelder-Mead', jac=None, hess=None)
            Lcurx = copy.deepcopy(Lcurall.x)
            Lcurx = np.asarray(Lcurx)
            Lcur = SSD(Lcurx, func, dataset, fixfix, Lcurx)
                      
            count += 1 
            #this if is just for debug now, basically: if 0=lower, 1=upper
            if g == 0:
                RES1 = copy.deepcopy(Lcurx)
#                for m in range (0,4):
#                    RES1[m] = RES1[m] - pars[m]
                count1 = copy.deepcopy(count)
            else:
                RES2 = Lcurx
            #FIX AASDEL and make a CHANGEABLE ONE
    return RES, RES1, RES2, count1

    

            