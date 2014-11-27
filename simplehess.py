# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 13:06:26 2014

@author: Phillip
"""
def simplehess(coeffs, YminFixed, YmaxFixed)
    
    from collections import defaultdict
    if YMinfixed:
        i = 1
    else:
        i = 0
        if Ymaxfixed:
            j=1
            else:
                j= 0

    freepar = string.lowercase[:len(self.coeffs)-i-j]
    hessindmat = defaultdict(freepar,freepar)
    for ind, index in enumerate(len(self.coeffs)-i-j):
        print index
        for und, undex in enumerate(len(self.coeffs)-i-j):
            hessindmat(index,undex) = [[freepar(index)],[freepar(undex)]]
      #  hessindmat(index,undex) = 
    print hessindmat