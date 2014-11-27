# -*- coding: utf-8 -*-
"""
Created on Wed Oct 01 15:36:35 2014

@author: Phillip
"""

#differentiate some function f(x) with respect to x using arit hmetic approximation

import numpy
import scipy

coeffs = 
step = 1000
coeffs = np.ndarray(shape=(1,5), dtype = float)
coeffs = (2,3,4,5,6) #feed from Wills
hessian = np.ndarray(shape=(len(coeffs),len(coeffs)), dtype = float) 

def funky(coeffs):
    #defines the function computed with parameters in coeffs
    f = 2*coeffs(0)+3*coeffs(1)+4*coeffs(2)+5*coeffs(3)+6*coeffs(4)
    return f
    
def dfunk(coeffs,i):
    
    #for i in range (1:len(coeffs))
    dcoeffs = coeffs
    dcoeffs(i)= coeffs(i)+coeffs(i)/step
    df = (funky(dcoeffs) - funky(coeffs))/(coeffs(i)/step)
    return df
    
def ddfunk(coeffs,i,j):
    ddcoeffs = coeffs
    ddcoeffs(j) = coeffs(j)+coeffs(j)/step
    ddf = (dfunk(ddcoeffs,i) - dfunk(coeffs,i))/(coeffs(j)/step)
    return ddf
    
for i in range(1:len(coeffs)):
    for j in range(1:len(coeffs)):
        hessian(i,j) = ddfunk(coeffs,i,j)

