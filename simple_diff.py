# -*- coding: utf-8 -*-
"""
Created on Wed Oct 01 15:36:35 2014

@author: Phillip
"""

#differentiate some function f(x) with respect to x using arit hmetic approximation

import numpy
import scipy

print 'what is the value of X,Y'
x = float(raw_input('X:'))
y = float(raw_input('Y:'))
step = 1000    
    
    
def funky(x,y):
     
    f = 2*x*x*x*y*y + y*y*y
   # print 'f=', f
    return f
    
funky(x,y)

def dfunk(x,y):
    #f = 
    def dfx (x,y):
        dfx = (funky((x+x/1000), y) - funky(x,y))/(x/1000)
    def dfy (x,y):
        dfy = (funky(x, (y+y/1000)) - funky(x,y))/(y/1000)
   # print "f'x=", dfx
   # print "f'y=", dfy
    return dfunk.dfx
    return dfunk.dfy
    
    def dfunk(cdoeffs, dcoeffs, ddcoeffs, concentration, Component, YminFixed, YmaxFixed):
                
        dfx = (hillequation(ddcoeffs, concentration, Component, YminFixed, YmaxFixed) - hillequation(dcoeffs, concentration, Component, YminFixed, YmaxFixed))
        
    def ddfunk(cdoeffs, hessindma, ondex,endex, concentration, Component, YminFixed, YmaxFixed):
        
        #here, d is a prefix used to mark the number of times +x/step was applied
        dcoeffs = cdoeffs
        ddcoeffs = cdoeffs
        #creating coeffs arrays of the form [x+x/step...]
        ddcoeffs[int(hessindma[endex][undex][0])] = ddcoeffs[int(hessindma[endex][undex][0])]+ddcoeffs[int(hessindma[endex][undex][0])]/step
        ddcoeffs[int(hessindma[endex][undex][1])] = ddcoeffs[int(hessindma[endex][undex][1])]+ddcoeffs[int(hessindma[endex][undex][1])]/step
        dcoeffs[int(hessindma[endex][undex][1])] = dcoeffs[int(hessindma[endex][undex][1])]+dcoeffs[int(hessindma[endex][undex][1])]/step
        
        ddfx = (dfunk(cdoeffs, dcoeffs, ddcoeffs, concentration, Component, YminFixed, YmaxFixed)-dfunk(cdoeffs, dcoeffs, ddcoeffs, concentration, Component, YminFixed, YmaxFixed))
        ddfx2 = (hillequation(ddcoeffs, concentration, Component, YminFixed, YmaxFixed)-hillequation(dcoeffs, concentration, Component, YminFixed, YmaxFixed))
        return ddfx,ddfx2
    
   # print "f'x=", dfx
   # print "f'y=", dfy
    return dfunk.dfx
    return dfunk.dfy
    
dfunk(x,y)

def ddfunk(x,y):
    #f = 
    ddfx = (dfunk.dfx((x+x/1000), y) - dfunk.dfx(x,y))/(x/1000)
    ddfy = (dfunk.dfy(x, (y+y/1000)) - dfunk.dfy(x,y))/(y/1000)
    #ddfxy = (funky((x+x/1000),(y+y/1000)) - funky(x,y))/(x+y)
    ddfxy = (dfunk.dfx(x,(y+y/1000)) - dfunk.dfx(x,y))/(y/1000)
    print "f''x=", ddfx
    print "f''y=", ddfy
    #print "f''xy=", ddfxy
ddfunk(x,y)


