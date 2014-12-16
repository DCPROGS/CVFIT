	subroutine DO_AoV(xobs,yobs,w,jval,nj,niobs,njset,
     &  kmax,theta,setx,kfit)
c To print analysis of variance fo data set #jval, in cases where there
c are enough d.f. 'within x values' to allow it
c Getgroup already done in CVIN2 but simplest to repeat it here
c This is used only for separate fit to a single set (ifitmode=1 or 2)
c
c NB don't need xcal, ycal here -just need calculated value that
c corresponds to each observed x value
c
	real*4 setx(njset)
	real*4 theta(100)
	integer nj(njset)
c declare MS as real
	real*4 MSerr,MSdev,MSreg,MScol,MStot
c
	real*4 xobs(niobs,njset),yobs(niobs,njset),w(niobs,njset)
	real*4 xgrp(niobs,njset),ybar(niobs,njset),sdm(niobs,njset)
	real*4 sdy(niobs,njset)
	integer ngrp(niobs,njset)
	logical discprt
	common/dp/discprt
c
	jset=jval
	call GETGROUP(xobs,yobs,w,jset,nj,niobs,njset,nxv,
     & xgrp,ybar,sdy,sdm,ngrp,SSwg,ndfwg,SStot,SScol,ymean,n0)
c
	if(ndfwg.lt.1) goto 99
	MSerr=SSwg/float(ndfwg)
c
c The w(i,j) are as in original data, but here we want to use sdy(i,jset)
c as weights, where sdy(k,jset) is the sd for the Y values at the kth X value
c in set #jset.
c Now calculate the SSD due to regression,
c SSDreg=SUM{k=1,nxval(k)}SUM{i=1,nxgrp(k,j)}[w(k)*(Ycalc-ybar)].
c  Rather than summing over all nj(j) values in the set, we sum separately over
c each group of Y values that have the same X.
	if(ndfwg.lt.4) RETURN	!not enough df to be worthwhile
c
	SSDreg=0.0
	SSDdev=0.0
	ival=1	!so xv1 in call to ycalcv is used
	do k=1,nxv	!sum for each different X value (=nxval(j) in main prog)
	   xv1=xgrp(k,jset)
	   Yc=YCALCV(kmax,theta,xv1,ival,j,setx,njset)
	   n=ngrp(k,jset)
	   ybarj=ybar(k,jset)
	   sd=sdy(k,jset)
	   if((sd.gt.0.).and.(sd.lt.3.e18)) then
		wt=1.0/(sd*sd)	!same for all x in the group
		do i=1,n
		   SSDreg=SSDreg + wt*(Yc - ymean)**2
		   SSDdev=SSDdev + wt*(Yc - ybarj)**2
		enddo
	   endif
	enddo
c
	print 1
	if(discprt) write (8,1)
1	format(/,
     & ' ANALYSIS OF VARIANCE',/,
     & '  Source            d.f.   SSD        MS           F       P',/,
     & ' ------------------------------------------------------------')
c
	ndferr=ndfwg
c
	ndfreg=kfit
	MSreg=SSDreg/float(ndfreg)
	F=MSreg/MSerr
	call FDIST(F,ndfreg,ndferr,P)
	print 2, ndfreg,SSDreg,MSreg,F,P
	if(discprt) write (8,2) ndfreg,SSDreg,MSreg,F,P
2	format(
     & '    Due to reg    ',i5,3g12.5,g12.5)
c
	ndf=nxv-1-ndfreg
	if(ndf.lt.1) goto 99
	MSdev=SSDdev/float(ndf)
	F=MSdev/MSerr
	call FDIST(F,ndf,ndferr,P)
	print 3, ndf,SSDdev,MSdev,F,P
	if(discprt) write (8,3) ndf,SSDdev,MSdev,F,P
3	format(
     & '    Dev from reg  ',i5,3g12.5,g12.5,/,
     & ' ------------------------------------------------------------')
c Between cols=X values
	ndf=nxv-1
	if(ndf.lt.1) goto 99
	MScol=SScol/float(ndf)
	F=MScol/MSerr
	call FDIST(F,ndf,ndferr,P)
	print 4, ndf,SScol,MScol,F,P
	if(discprt) write (8,4) ndf,SScol,MScol,F,P
4	format(
     & ' Between X values ',i5,3g12.5,g12.5)
c Error (within x values)
	MSerr=SSwg/float(ndfwg)
	ndf=ndfwg
	if(ndf.lt.1) goto 99
	print 5, ndf,SSwg,MSerr
	if(discprt) write (8,5) ndf,SSwg,MSerr
5	format(
     & ' Within X values  ',i5,2g12.5,/,
     & ' ------------------------------------------------------------')
c Total
	ndf=nj(jset)-1-n0		!n0=number with zero weight
	MStot=SStot/float(ndf)
	if(ndf.lt.1) goto 99
	print 6, ndf,SStot,MStot
	if(discprt) write (8,6) ndf,SStot,MStot
6	format(
     & ' Total            ',i5,2g12.5,/,
     & ' ------------------------------------------------------------')
c
	RETURN
c
99	continue
	print 7
	if(discprt) write (8,7)
7	format(' Cannot complete analysis of varinace (d.f. = 0)')
	RETURN
c
	end
