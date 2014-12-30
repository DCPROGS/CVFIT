#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="remis"
__date__ ="$25-Oct-2010 20:29:34$"

import numpy as np
import math


def nwvrtx(i, A, simp, k):
    """
    subroutine NWVRTX(I,A,simp,k,kt1,kt)
    dimension A(kt),simp(kt1,kt)
    copies A() into row I of simp()
    """
    for j in range(0, k):
        simp[i, j] = A[j]
    return simp

def abstst(fval, theta, absmin, k):
    """
    subroutine ABSTST(fval,pnew,absmin,thmin,k,kd1,kd2)
    dimension pnew(kd1),thmin(kd2)
    """
    thmin = np.array((k))
    if fval > absmin:
        return absmin, theta
    else:
        return fval, thmin


def sortsimp(simp, fval, n):
    """
    """

    for i in range(0, n):
        for j in range(0, n-1):
            temp1 = np.array((n-1))
            temp2 = np.array((n-1))
            if fval[j+1] < fval[j]:
                tempv = fval[j]
                fval[j] = fval[j+1]
                fval[j+1] = tempv
                temp1 = simp[j+1]
                temp2 = simp[j]
                print temp1,temp2
                simp[j+1] = temp2
                simp[j] = temp1
        print 'interm simp'
        print simp
        print 'fval=', fval
    return simp, fval


def simplex(theta, func, func1, opt, data, verbose=0):
    """
	call SIMPLEXv(kmax,THETA,stpfac,errfac,nev,nevm,
     & smin,SSDCV,Ndisp,jfix,delmin,confac,irestrt,iconv,
     & Xobs,yobs,w,nj,juse,setx,niobs,njset,ndth)
     subroutine SIMPLEXv(kt,THETA,stpfac,errfac,neval,nevmax,
     & fmin,FUNC,Ndisp,jfix,delmin,confac,irestrt,iconv,
     & Xobs,yobs,w,nj,juse,setx,niobs,njset,kn)

    """

    # these might come as parameters
    irestrt = 3
    errfac = 1.e-4
    delmin = -1.    #		!do not use delmin for convergence
    stpfac=0.1

    nresmax = irestrt    # save input value=max number of restarts
    istart = 0    # counts number of restarts
# Set ndel=no of consec iterations for which reduction in func<delmin for convergence
    ndel = 5
    idel = 0    # counts up to ndel


    print ' USING FAST VERSION OF SIMPLEX'

    reffac = 1.0     # reflection coeff => 1
    confac = 0.5     # contraction coeff 0 < beta < 1
    extfac = 2.0     # extension factor > 1
    neval = 0	 # counts function evaluations

    k = np.size(theta)
    n = k + 1    # # of vertices in simplex
    #simp = np.zeros((n, k))
    simp = np.zeros((n, k))
    fval = np.zeros((n))
    fval1 = np.zeros((n,1))
    step = np.zeros((k))
    crtstp = np.zeros((k))
    thmin = np.zeros((k))
    pnew = np.zeros((k))
    pnew1 = np.zeros((k))
    centre = np.zeros((k))

    for j in range(0, k):
        step[j] = stpfac * theta[j]
	crtstp[j] = errfac * theta[j]

    restart = True
    while restart:    # 2001	continue	!return here for restart

        #simp = nwvrtx(0, theta, simp, k)     #start values=vertex #1
        simp[0] = theta
        fval[0] = func1(theta, data, opt, func)
        fval1[0,0] = fval[0]
        fsav = fval[0]
        neval = neval + 1
        absmin = fval[0]    # starting value

        # compute offset of the vertices of the starting simplex
	fac = (math.sqrt(float(n))-1.) / (float(k)*math.sqrt(2.))

        #specify all other vertices of the starting simplex
        for i in range(1, n):
            for j in range(0, k):
                simp[i,j] = simp[0,j] + step[j]*fac
            simp[i, i-1] = simp[0, i-1] + step[i-1]*(fac+1./math.sqrt(2.))

            #  and calculate their residuals
            i1 = i
            theta = simp[i]
            fval[i] = func1(theta, data, opt, func)
            fval1[i,0] = fval[i]
            if fval[i] < absmin:     #  test for absolut minimum
                absmin = fval[i]
                thmin = theta
        neval = neval + k

        print ' simplex at the beginning of restart='
        print simp
        print ' sds at the begining', fval

	ndsav=0
	niter	= 0
        # start iteration loop here at 2000
        L = 0
        while L == 0:
            niter = niter + 1

            # find best (lowest fval) and worst (highest fval) vertices-
            # indeces of these are ILO,IHI respectively.
            # fnewlo=fval(IHI) corresponds to the vertex to be replaced.

            flo = fval[0]
            fhi = flo
            ilo = 0
            ihi = 0
            for i in range(0, n):
                if fval[i] < flo:
                    flo = fval[i]
                    ilo = i
                if fval[i] > fhi:
                    fhi = fval[i]
                    ihi = i
            theta = simp[ilo]
            # display current best vertex
            print 'iter#', niter, 'best fval=', flo
            print 'theta=', simp[ilo]

            # if(neval.gt.nevmax) goto 1510

            # ----- compute centroid of all vertices except the worst
            for i in range(0, n):
                if i != ihi:
                    for j in range(0, k):
                        centre[j] = centre[j] + simp[i,j]

            # ----- reflect, with next vertex taken as reflection of worst
            # Parameter values that are coord of new vertex in pnew(j)
            for j in range(0, k):
                centre[j] = centre[j] / float(k)
                pnew[j] = centre[j] - reffac * (simp[ihi,j]-centre[j])
            fnew = func1(pnew, data, opt, func)
            neval = neval + 1
            print 'e#', neval, 'val=',fnew, 'reflection'
            print 'pnew=', pnew

            if fnew < absmin:     #  test for absolut minimum
                absmin = fnew
                thmin = pnew

	#if (fnew.ge.flo) goto 112	!fnew worst than prev best

            if fnew < flo:
            # ----- new vertex is better than previous best so extend it
                for j in range(0, k):
                    pnew1[j] = centre[j] + extfac * (pnew[j] - centre[j])
                fnew1 = func1(pnew1, data, opt, func)
                neval = neval + 1
                print 'e#', neval, 'val=',fnew1, 'extention'
                print 'pnew=', pnew1

                if fnew1 < absmin:     #  test for absolut minimum
                    absmin = fnew1
                    thmin = pnew1
                if fnew1 < fnew:     # ----- still better
                    simp = nwvrtx(ihi, pnew1, simp, k)     #start values=vertex #1
                    fval[ihi] = fnew1
                else:
                    simp = nwvrtx(ihi, pnew, simp, k)     #start values=vertex #1
                    fval[ihi] = fnew

                # goto 1901 for convergence check

            else:     # 112  come here if reflected vertex not
                      # better than best vertex, so no extension wanted
                L1 = 0     # number of vertices which are worse than fnew
                for i in range(0, n):
                    if fval[i] > fnew: L1 = L1 + 1
                if L1 == 2:     # reflected vertex better than 2 or more
                                # so input pnew and test conv
                    simp = nwvrtx(ihi, pnew, simp, k)     #start values=vertex #1
                    fval[ihi] = fnew

                if L1 < 2:
                    if L == 1:     # fnew better than one (the worst) vertex so
#c input unextended reflection, pnew, then contract on new (reflected) side
#c (after vertex ihi replaced by reflected vertex loop at 270 does this)
                        simp = nwvrtx(ihi, pnew, simp, k)     #start values=vertex #1
                        fval[ihi] = fnew
                # Contract on the original fval(IHI) side of the centroid
                for j in range(0, k):
                    pnew1[j] = centre[j] +confac * (simp[ihi, j] - centre[j])
                fnew1 = func1(pnew1, data, opt, func)
                neval = neval + 1
                print 'e#', neval, 'val=',fnew1, 'contraction'
                print 'pnew1=', pnew1

                if fnew1 < absmin:     #  test for absolut minimum
                    absmin = fnew1
                    thmin = pnew1
                # ----- is contracted vertex better than the worst vertex
                if fnew1 <= fval[ihi]:
                    simp = nwvrtx(ihi, pnew1, simp, k)     #start values=vertex #1
                    fval[ihi] = fnew1


                else:
                    #  ----- no, it is still bad, shrink whole simplex towards best vertex
                    for i in range(1, n):
                        for j in range(0, k):
                            #simp[i,j] = (simp[i,j] + simp[ilo,j]) * 0.5    # orig version
                            simp[i,j] = simp[ilo,j] + confac * (simp[i,j] - simp[ilo,j])
                            # last line is D.C. version that uses confac rather than 0.5
                            theta[j] = simp[i,j]
                        fval[i] = func1(theta, data, opt, func)
                        neval = neval + 1
                        print 'e#', neval, 'fval=',fval[i], 'pnew1=', theta, 'shrink'
                        if fval[i] < absmin:     #  test for absolut minimum
                            absmin = fval[i]
                            thmin = theta

                print 'fnew1 < fhighest'
                print 'simp='
                print simp
                print 'fvals=', fval



            # CHECK CONVERGENCE. IF NOT CONVERGED GOTO 2000.
            # This version uses diff between highest and lowest value
            # of parameter of the n values that define a vertex
            # (as in O'Neill version)

            #  ----- order the vertices for all vertices
            # Define L=0 for not converged- do next iteration
            # L=1 for converged via crtstp
            # L=2 for converged via delmin (no restarts)
            # L=3 for abort (no restarts)

            if delmin < 0:
                L = 1    #  conv via crtstp
                for j in range(0, k):     # test each parameter
                    il = 0
                    ih = 0
                    for i in range(0, n):    # order values of current param
                        if simp[i,j] < simp[il,j]: il = i
                        if simp[i,j] > simp[ih,j]: ih = i
                    if(simp[ih,j] - simp[il,j]) > abs(crtstp[j]): L = 0    # not conv

            else:    # Use delmin criterion- find f at current best vertex
                il = 1
                for i in range(0,n):
                    if fval[i] < fval[il]: il = i
                del1 = fsav - fval[il]    # reduction in min=pos value
                #IF(DEL1.LT.-1.E-5) GOTO 184	!MIN IS INCREASING! DEL NEGATIVE

                if del1 > -1.0e-5:
                    if del1 >= delmin:
                        idel = 0    # reset when del1 > delmin
                        fsav = fval[il]    # SAVE LAST MINIMUM
                    if del1 < delmin:
                        idel = idel + 1
                        if idel < ndel:    # not yet converged
                            fsav = fval[il]    # SAVE LAST MINIMUM
                        else:
                            print ' Converged via DELMIN '
                            L = 2    # converged


        # end of while L == 0:
	# if (L.eq.0) goto 2000	!next iteration

        # ----- convergence attained. Options for ending in this version are:
        # 	(1)look at current best vertex
        # 	(2)look at param values averaged over vertices
        # 	(3)look at absmin,thmin. If better, restart at absmin, as below.
        # 	(4)do local search with each param +/- crtstp, as in O'Neill
        # 	 version, starting at current best vertex. If none are better
        # 	 input current best vertex. If some better restart at better
        # 	 value with crtstp taken as approptiately small initial step.

        #  get best vertex

        if fval[ihi] > fval[ilo]: ihi = ilo    # addition by I.D.HILL
        temp = np.array((k))
        temp = simp[ihi]
        fhi = fval[ihi]
	fval[0] = fhi    # for test below
        # this is best vertex

        # check absmin
	fval[2] = absmin    # for test below

        if L !=1:
            if fval[0] < fval[2]:     # goto 1254  !aborted- use best vertex
                print ' Return with best vertex'
                iconv = 2    # signals 'return with best vertex'
                theta = temp
                fmin = fval[0]
                return theta, fmin
            else:    #1387
                print ' Return with absmin'
                iconv = 4    # signals 'return with absmin'
                theta = tmin
                fmin = fval[2]
                return theta, fmin

        else:       #if L == 1:      # not aborted       1386
            # next average over vertices-put values in pnew()
            for j in range(0, k):
                pnew[j] = 0.0
                for i in range(0, n):
                    pnew[j] = pnew[j] + simp[i,j]
                pnew[j] = pnew[j] / float(n)
            fval[2] = func1(pnew, data, opt, func)

            # do local search. Temp() already contains best
            # vertex after convergence, with corresp function
            # value in fval(1). Put altered values in pnew1

            for j in range(0,k):
                pnew1[j] = temp[j] + crtstp[j]
                fval[3] = func1(pnew1, data, opt, func)
                if fval[3] > fhi:
                    pnew1[j] = temp[j] - crtstp[j]     # step in other direction
                    fval[3] = func1(pnew1, data, opt, func)
                    if fval[3] > fhi:
                        print ' local search does not improve on theta '
                        fval[3] = fval[0]    # restore orig func value

            # Test which is best. Want to restart at thmin
            # (fval(3))if this is best or restart at pnew1
            # (fval(4))if this is best. Otherwise exit with theta
            # (fval(1)) or with average (fval(2)), whichever is better.
            # For restart new initial values must be in Temp,
            # and step reduced to crtstp.

            # 125	continue
            # find best value

            il = 0
            for i in range(0, 4):
                if fval[0] < fval[il]: il = i

            if istart < nresmax and (il == 2 or il == 3):
                if il == 2:
                    # goto 1252	!restart at thmin
                    print ' Restart at absmin'
                    irestrt = 1     # signals restart at absmin
                    istart = istart + 1     # count number of restarts
                    for j in range(0, k):
                        temp[j] = thmin[j]
                        step[j] = crtstp[j]    # small step size for restart
                    # goto 2001	!restart

                if il == 3:
                    # goto 1253	!restart at pnew1
                    print ' Restart after local search'
                    irestrt = 2    # signals restart after local search
                    istart = istart + 1    # count number of restarts
                    for j in range(0, k):
                        temp[j] = pnew1[j]
                        step[j] = crtstp[j]    # small step size for restart
                    # goto 2001	!restart

                else:
                    if il == 0:
                        # goto 1254	!exit with best vertex=temp
                        print ' Return with best vertex'
                        iconv = 2    # signals 'return with best vertex'
                        theta = temp
                        fmin = fval[0]
                        return

                    if il == 1:
                        #goto 1255	!exit with average
                        print ' Return with averaged vertices'
                        iconv = 3    # signals 'return with average vertex'
                        theta = pnew
                        fmin = fval[1]
                        return

                    if il == 2:
                        #goto 1387	!exit with thmin
                        print ' Return with absmin'
                        iconv = 4    # signals 'return with absmin'
                        theta = thmin
                        fmin = fval[2]
                        return

                    if il == 3:
                        # goto 1388	!exit with pnew1
                        print ' Return with result of local search'
                        iconv = 5    # signals 'return with pnew1'
                        theta = pnew1
                        fmin = fval[3]
                        return


if __name__ == "__main__":
    print "Hello World";


