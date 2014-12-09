# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 14:14:18 2014

@author: Phillip
"""
import numpy as np

   # for h1 in range(len(covar)):
   #     corlen += h1-1
   # for k1 in range(0,corlen):    
   #     correl[k1] = 
def confidence_interval(covar, dataset, func, pars, notfixed):
    putin = 1 #this is input from user: number of devs from mean
    args = pars
    df = len(covar)
    t = stats.t.ppf(0.975, df)
    aaSD = approximateSD(dataset, func, args, notfixed)
    aSD = aaSD * putin
    Smin = SSD(dataset, func, args)
    contour = Smin+t**2
    lolim = copy.deepcopy(args)
    hilim = copy.deepcopy(args)
    y = 1.071773463
    ym1 = 1
    count1 = 1
#    x = 1
#    deltax = 1
    for i in range(1,4):
        while SSD(dataset, func, lolim) < contour and count1<100:
            lolim[i] += aSD[i]
            lowli, aSD, y, ym1, count1 =  seek_bisect(lolim, i, args, aSD, aaSD, dataset, func, contour)
    print '\n lower interval = \n', lowli
    
    for i in range(1,4):
        while SSD(dataset, func, hilim) < contour and count1<100:
            lolim[i] += aSD[i]
            higli, aSD, y, ym1, count1 =  seek_bisect(hilim, i, args, aSD, aaSD, dataset, func, contour)
    print '\n higher interval = \n', higli
    return

def new_sec(aSD, covar, dataset, func, pars, notfixed):
    df = len(covar)
    t = stats.t.ppf(0.975, df)
    aaSD = approximateSD(dataset, func, pars, notfixed)
    j = len(pars) - len(aaSD)
    i = 2
    Smin = SSD(dataset, func, pars)
    contour = Smin+t**2
    rile = [1,-1] # either go left from Smin or right
    for g in range (0,1):
        STEP = 1.071773463
        lohi = rile[g]
        count = 0
        SDel = np.zeros(len(aaSD)+j)
        SDel[i+j] = SDel[i+j] + (lohi * aaSD[i])
        SDel = SDel + pars # parameters + initial guess for one parameter
        aasdlist = np.zeros(1000) #for diagnostics and finding oscillatory behaviour
        tbl = bisectbl()
        while np.sqrt((SSD(dataset, func, SDel) - contour)**2) >  0.001*Smin and count <1000:
            aaSDel = copy.deepcopy(SDel) #the current estimate
            
            if SSD(dataset, func, aaSDel) - contour > 0:
                n = 1
           # elif SSD(dataset, func, aaSDel) - contour == 0:
           #     break
            else:
                n = -1
            aasdlist[count] = copy.deepcopy(aaSDel[i+j])
            n3 = copy.deepcopy(n2)
            n2 = copy.deepcopy(n)
            if n2 != n and n3 != n:
                STEP /= 1.071773463
            aaSDel[i+j] = aaSDel[i+j] - lohi*n*(aaSDel[i+j]*STEP - aaSDel[i+j])
            STEP *= 1.071773463
            # problem: on reverse of n, STEP still grows, hence oscillation results in huge numbers
            if count in tbl:
                STEP = 1.071773463
            if aasdlist[count] == aasdlist[count-2]:
                SDel[i+j] = copy.deepcopy(aasdlist[count-2])
            count += 1  
            #FIX AASDEL and make a CHANGEABLE ONE
    return SDel, count, aasdlist
    
def bisectbl():
    a = np.zeros(10)
    b = np.zeros((10,10))
    for f in range (0,10):
        a[f] = 100*f
        for k in range (0,10):
            b[f][k] = a[f]+10*k
    return b
    
def new_secs2(aSD, covar, dataset, func, pars, notfixed):
    df = len(covar)
    t = stats.t.ppf(0.975, df)
    print '\n t = \n', t
    aaSD = approximateSD(dataset, func, pars, notfixed)
    j = len(pars) - len(aaSD)
    i = 2
    Smin = 0.5*SSD(dataset, func, pars)
    contour = Smin+t**2
    rile = [1,-1] # either go left from Smin or right
    RES = np.zeros(2)
    for g in range (0,2):
        STEP = aaSD[i]*2
        lohi = copy.deepcopy(rile[g])
        count = 0
        SDel = np.zeros(len(aaSD)+j)
        SDel[i+j] = SDel[i+j] + (lohi * STEP)
        SDel = SDel + pars # parameters + initial guess for one parameter
        print '\n STEP, lohi = \n', STEP, lohi
       
        while np.sqrt(((0.5*SSD(dataset, func, SDel)) - contour)**2) >  0.0005*Smin and count <1000:
            print '\n SS, COntr = \n', 0.5 * SSD(dataset, func, SDel), contour
            #aaSDel = copy.deepcopy(SDel) #the current estimate
            if (0.5 * SSD(dataset, func, SDel)) - contour > 0:
                n = -1
            else:
                n = 1
            STEP /= 2
            print '\n Step = \n', STEP
            SDel[i+j] += n * lohi * STEP
            print '\n SDel = \n', SDel
            count += 1  
            RES[g] = copy.deepcopy(SDel[i+j])
            #FIX AASDEL and make a CHANGEABLE ONE
    return RES, count
    
    def seek_interval_less_sd(lolim, i, aSD, aaSD, y, ym1):
    lolim[i] -= (-aaSD*y + aaSD*ym1)
    #ym1 = copy.deepcopy(y)
    ym1 *= ym1
    return lolim, y, ym1
    
def seek_interval_grth_sd(lolim, i, aSD, aaSD, y, ym1):
    lolim[i] += (aaSD*y - aaSD*ym1)
    #ym1 = copy.deepcopy(y)
    ym1 *= ym1
    return lolim, y, ym1
    
def seek_bisect(lolim, i, args, aSD, aaSD, dataset, func, contour):
    count = 0
    count1 = 0
    y = 1
    ym1 = 1.071773463
    if lolim[i] < args[i]+aSD[i]:
        lolim +=aSD(i)
        SSAD = SSD(dataset, func, lolim)
        while SSAD < contour and count<100:
            lolim, y, ym1 = seek_interval_less_sd(lolim, i, aSD, aaSD, y, ym1)
            count +=1
            SSAD = SSD(dataset, func, lolim)
    else:
        lolim -=aSD[i]
        SSAD = SSD(dataset, func, lolim)
        while SSAD < contour and count<100:
            lolim, y, ym1 = seek_interval_grth_sd(lolim, i, aSD, aaSD, y, ym1)
            count +=1
            SSAD = SSD(dataset, func, lolim)
        count1 += 1
    return lolim, aSD, y, ym1, count1