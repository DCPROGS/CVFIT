#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="remis"
__date__ ="$05-Nov-2010 21:19:36$"

import numpy as np
from math import*

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

def simplex(theta, func, data, verbose=0):
    """
    Python implementation of DC's SIMPLEXV.FOR subroutine.
    """

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
    nrestart = 100    # max number of restarts
    irestart = 0    # counts restarts

    while irestart < nrestart:    # 2001	continue	!return here for restart

        irestart += 1
        print 'RESTART#', irestart

        simp[0] = theta
        fval[0] = func(theta, data)
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
            fval[i] = func(simp[i], data)
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
            fnew = func(pnew, data)
            if fnew < absmin:
                absmin = fnew
                thmin = pnew
            neval += 1
            if verbose: print 'reflection: e#', neval, 'f=',fnew, 'pnew=', pnew

            if fnew < fval[0]:
                # ----- new vertex is better than previous best so extend it
                for j in range(k):
                    pnew1[j] = centre[j] + extfac * (pnew[j] - centre[j])
                fnew1 = func(pnew1, data)
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
                    fnew1 = func(pnew1, data)
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
                            fval[i] = func(simp[i], data)
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
            fvalav = func(pnew, data)
            exvals.append(fvalav)

            exvals.append(absmin)

            # do local search. Put altered values in pnew1
            for j in range(k):
                pnew1[j] = 0.0
                pnew1[j] = simp[0,j] + locfac * crtstp[j]
            fval1 = func(pnew1, data)
            if fval1 < fval[0]:
                exvals.append(fval1)
            else:
                for j in range(k):
                    pnew1[j] = simp[0,j] - locfac * crtstp[j]     # step in other direction
                fval1 = func(pnew1, data)
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


def linreg(x, y, n, imode = 0):
    """
    Python implementation of DC's subroutine LINREG(x,y,n,a,b,imode).
    Return a,b for conventional linear regression.
    imode=0 for normal fit.
    imode=1 for fit of slope with a=input value (e.g. a=0 or a=1 for Schild eq).
    """

    a = 1

    en = float(n)
    sx = 0.0
    sxx = 0.0
    sy = 0.0
    syy = 0.0
    sxy = 0.0
    for i in range(n):
        sx = sx + x[i]
        sxx = sxx + x[i] * x[i]
        sy = sy + y[i]
        syy = syy + y[i] * y[i]
        sxy = sxy + x[i] * y[i]
    sxx = sxx - (sx * sx) / en
    sxy = sxy - (sx * sy) / en
    b = sxy / sxx
    if imode == 0: a = (sy / en) - b * sx /en

    return a, b

def ssdcv(theta, data,opt, func):
    """
    Python implementation of SSDCV.FOR
    SUM OF SQUARES FOR FITTING MODEL #NMOD WITH CVFIT
    function SSDCV(K,THETA,Xobs,yobs,w,nj,juse,setx,niobs,njset)
    """

    s = 0.0

    if opt['logyfit']:
        pass


    if opt['ifitmode'] < 3:  # Only one set fitted or sets fitted separately
        n = np.size(data[0])     # number of points in a set
        for i in range(0, n):
            x = data[0,i]
            y = func(theta, x, opt)
            w = 1 / (data[2,i] * data[2,i])
            dev = data[1,i] - y
            s = s + w * dev * dev
        return s

    if opt['ifitmode'] >= 3 or opt['ifitmode'] <= 5:
        #fit all sets at once for ifitmode=3
        iset = 0
        for m in range(0, opt['nsfit']):
            j = opt['juse'][m]     # set # used
            iset = iset + 1
            # iset is in COMMON -nec for YCALCV when ifitmode=4,5 if some sets missed so
            # c parameter # in theta() can be matched with proper data set
            for i in range(0, opt['nj'][j]):
                # need some coding if wangt some observations to omit
                # if(nomit.gt.0) then
                #   omit specified obs. Is I equal to any of the elements of JOMIT?
		#   do L=1,nomit
		#	if(i.eq.jomit(L)) goto 4
		#   enddo
		# endif
		# xv1 = Xobs(i,j)
		# Y = ycalcv(kmax,theta,xv1,ival,j,setx,njset)
                # dev = Yobs(i,j) - Y     #normal calc
		#S = S + w(i,j)*dev*dev

                x = data[0,i]
                y = func(theta, x, opt)
                w = 1 / (data[2,i] * data[2,i])
                dev = data[1,i] - y
                s = s + w * dev * dev
		nfit = nfit + 1
        return s

def calcSDWeights(data):
    """
    Calculate weights from specified errors.
    """

    k = np.size(data[0])
    for i in range(0, k):
        data[3,i] = 1 / (data[2,i] * data[2,i])

    return data

if __name__ == "__main__":
    print "Hello World";
