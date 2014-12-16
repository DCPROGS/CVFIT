	subroutine GETGUESS(nmod,Xobs,Yobs,nj,juse,niobs,njset,nsfit,
     & theta,ndth,thetgues,nset1,kmax,ybarp,ncomp,ifitmode,jset,
     & setx,jfix,nfix,logyfit,norm,sepmax,ymaxsep,guesdone,noguess)
c
c To set initial guesses in CVSIMP
c Initial guesses for set j (separately) kept in thetgues(i,j)
c Initial guesses for simultaneous fit put straight into theta
c IN cases where guesses for the separate fits are already in thetgues()
c enter with guesdone=true, and this routine does not make guesses but
c uses the existing ones to set suitable guesses in theta() for ifitmode
c in use.
	real*4 xobs(niobs,njset),yobs(niobs,njset)
	integer juse(njset),nj(njset)
	real*4 theta(ndth)
	integer JFIX(ndth)
	real*4 thetgues(100,nset1)
	logical noguess(100),logyfit,norm,sepmax,fix0,guesdone
	real setx(njset)
	character defolt*30,cans*30		!to hold text & result of QDIALOG
	character*1 ans
	real*4 ymaxsep(nset1)
	allocatable::x1,y1	!local
	real*4 x1(:),y1(:)
	COMMON/BLOCK2/ castar,xdum1,xdum2,iequiv,ip1,ip2
	common/potrat/jfirst,iset,kmax1			!for pot ratios
c
	njmax=0
	ymin=1.e37
	ymax=-1.e37
	do j=1,nset1
	   noguess(j)=.false.
	   ymaxsep(j)=-1.e37
	enddo
	do m=1,nsfit
	   j=juse(m)		!set # used
	   if(nj(j).gt.njmax) njmax=nj(j)
	   do i=1,nj(j)
		if(Yobs(i,j).lt.ymin) ymin=Yobs(i,j)
		if(Yobs(i,j).gt.ymax) ymax=Yobs(i,j)
		if(Yobs(i,j).gt.ymaxsep(j)) ymaxsep(j)=Yobs(i,j)
	   enddo
	enddo
c  How constant are the maxima
	sy=0.0
	syy=0.0
	if(nsfit.ge.2) then
	   sepmax=.true.
	   do m=1,nsfit
		j=juse(m)		!set # used
		y=ymaxsep(j)
		sy=sy+y
		syy=sy+y*y
	   enddo
	   vary=(syy-sy*sy/float(nsfit))/float(nsfit-1)
	   ybar=sy/float(nsfit)
	   if(vary.gt.0.00001) then
		cv=sqrt(vary)/ybar
		sepmax=cv.gt.1.0	!use sep maxima if CV(ymax)>100%
	   endif
	else
	   sepmax=.true.
	endif
c Following cases are known to have same/diff maxima (for Hill-Langmuir)
	if(ifitmode.eq.4.or.ifitmode.eq.5.or.
     &	(ifitmode.eq.3.and.ip1.eq.0)) then
	   sepmax=.false.
	else if(ifitmode.eq.3.and.(ip1.eq.-1.or.ip1.eq.1)) then
	   sepmax=.true.
	endif
c===next line added 01/27/97 05:43pm
	if(nsfit.eq.1) sepmax=.false.
c
	ALLOCATE(x1(njmax),y1(njmax))
c
c Hill-Langmuir
c For Langmuir use x/y versus x; For Hill have to use Hill plot
c But quite often have y=0, so better have method that avoids 1/y
c eg Eadie-Hofstee  y versus y/x
c (NB for decreasing hyperbola, y*x versus y is similar)
c For ifitmode=3 have
c   ip1=0 -common max, separate K
c   ip1=1 -common K, separate max
c   ip2=0 -separate nH
c   ip2=1 -common nH
c	fix0=.false.
	if((nmod.eq.26.or.nmod.eq.1).and.ncomp.eq.1) then	!Langmuir
	   Y0=ymin/2.		!will this do!!
c	   ans='N'
c	   if(y0.lt.0.05*ymax) ans='Y'
	   ans='Y'
	   call DEFOLTa(ans,defolt)
	   if(iequiv.eq.0) then
		call QDIALOG(1,'Fix Y(0)=0',defolt,11,cans)
	   else
		call QDIALOG(1,'Fix Y(inf)=0',defolt,11,cans)
	   endif
	   call GETINPa(cans,ans)
	   if(ans.eq.'Y') then
c		fix0=.true.
		Y0=0.0
		jfix(1)=1
		nfix=1
	   endif
	   ymm=0.0		!for mean ymax
	   akm=0.0		!for mean K
	   n1=0
	   do m=1,nsfit
		j=juse(m)		!set # used
		if(guesdone) then
		   if(noguess(j)) then
			goto 97				!try next set
		   else
			n1=n1+1
			Y0=thetgues(1,j)	!-guesses for sep fits -see ycalcv
			Ym=thetgues(2,j)
			ak=thetgues(3,j)
		   endif
		else
c This Scatchard type plot can be awful (esp if 1st point too low) -better
c just use observed max!
		 ym=0.0
		 do i=1,nj(j)
		   if(Yobs(i,j).gt.ym) ym=Yobs(i,j)
		 enddo
		 yhalf=0.5*ym
		 if(iequiv.eq.0) then	!increasing hyperbola
		    do i=1,nj(j)-1
			y11=Yobs(i,j)
			y22=Yobs(i+1,j)
			if(y11.le.yhalf.and.y22.ge.yhalf) then
			   aK=0.5*(Xobs(i,j)+Xobs(i+1,j))
			endif
		    enddo
		 else		!decreasing hyperbola
		    do i=1,nj(j)-1
			y11=Yobs(i,j)
			y22=Yobs(i+1,j)
			if(y11.ge.yhalf.and.y22.le.yhalf) then
			   aK=0.5*(Xobs(i,j)+Xobs(i+1,j))
			endif
		    enddo
		 endif
c		 do i=1,nj(j)
c		   x1(i)=Xobs(i,j)
c		   if(iequiv.eq.0) then	!increasing hyperbola
c			y1(i)=Yobs(i,j)
c			x1(i)=Yobs(i,j)/Xobs(i,j)
c		   else			!decreasing hyperbola!
c			y1(i)=Yobs(i,j)*Xobs(i,j)
c			x1(i)=Yobs(i,j)
c		   endif
c		 enddo
		 n=nj(j)
		 if(n.lt.2) then
		   noguess(j)=.true.		!or try another method?
		   goto 97				!try next set
		 endif
		 n1=n1+1
c		 call LINREG(x1,y1,nj(j),a,b,0)
c		 n1=n1+1
c		 if(iequiv.eq.0) then
c		   Ym=a
c		   aK=-b
c		 else
c		   Ym=-a/b
c		   aK=-b
c		 endif
		 thetgues(1,j)=Y0	!-guesses for sep fits -see ycalcv
		 thetgues(2,j)=Ym
		 thetgues(3,j)=aK
		endif
		ymm=ymm+Ym		!for mean ymax
		akm=akm+aK		!for mean K
97		continue
	   enddo
	   ymm=ymm/float(n1)	!mean ymax
	   akm=akm/float(n1)	!mean K
	   if(ifitmode.le.2) then
		theta(1)=Y0
		theta(2)=thetgues(2,jset)	!Ym for 1st set to be fitted, #jset
		theta(3)=thetgues(3,jset)	!K for 1st set to be fitted, #jset
	   else if(ifitmode.eq.3) then
		if(ip1.eq.-1) then		!sep ymax, sep K
		   theta(1)=Y0
		   k=1
		   do m=1,nsfit
			j=juse(m)		!set # used
			k=k+1
			theta(k)=thetgues(3,j)  !K value for each set
			k=k+1
			theta(k)=thetgues(2,j)  !Ym value for each set
		   enddo
		else if(ip1.eq.0) then		!common ymax, sep K
		   theta(1)=Y0
		   theta(2)=ymm 	!mean Ymax for common Ymax guess
		   k=2
		   do m=1,nsfit
			j=juse(m)		!set # used
			k=k+1
			theta(k)=thetgues(3,j)  !K value for each set
		   enddo
		else if(ip1.eq.1) then	!common K, sep ymax
		   theta(1)=Y0
		   k=1
		   do m=1,nsfit
			j=juse(m)		!set # used
			k=k+1
			theta(k)=thetgues(2,j)  !Ym value for each set
		   enddo
		   k=k+1
		   theta(k)=akm	!mean K for common K guess
		   if(logyfit) theta(1)=theta(k)	!only one (see geteqn)
		endif
	   else if(ifitmode.eq.4) then
		theta(1)=Y0
		theta(2)=thetgues(2,jfirst)	!Ym for 1st set to be fitted
		theta(3)=thetgues(3,jfirst)	!K for 1st set to be fitted
		k=3
		do m=2,nsfit
		   j=juse(m)
		   k=k+1
	         theta(k)=thetgues(3,j)/thetgues(3,jfirst)	!approx for pot ratios
		enddo
	   else if(ifitmode.eq.5) then
		theta(1)=Y0
		theta(2)=thetgues(2,jfirst)	!Ym for 1st set to be fitted
		theta(3)=thetgues(3,jfirst)	!K for 1st set to be fitted
		do m=1,nsfit
		   j=juse(m)		!set # used
	         r=thetgues(3,j)/thetgues(3,jfirst)	!approx for pot ratio
		   xB=setx(j)
		   if(r.gt.1.2) then
			n=n+1
			s=s+xB/(r-1.)
		   endif
		enddo
		theta(4)=s/float(n)	!estimate of KB
	   endif
	else if((nmod.eq.27.or.nmod.eq.-1).and.ncomp.eq.1) then	!Hill
c       Try Hill plot between p=0.1 to p=0.9 (approx!)
	   Y0=ymin/2.		!will this do!!
c	   ans='N'
c	   if(y0.lt.0.05*ymax) ans='Y'
	   ans='Y'
	   call DEFOLTa(ans,defolt)
	   if(iequiv.eq.0) then
		call QDIALOG(1,'Fix Y(0)=0',defolt,11,cans)
	   else
		call QDIALOG(1,'Fix Y(inf)=0',defolt,11,cans)
	   endif
	   call GETINPa(cans,ans)
	   if(ans.eq.'Y') then
c		fix0=.true.
		Y0=0.0
		jfix(1)=1
		nfix=1
	   endif
	   ymm=0.0		!for mean ymax
	   akm=0.0		!for mean K
	   enm=0.0		!for mean nH
	   n1=0
	   do m=1,nsfit
		j=juse(m)		!set # used
		if(guesdone) then
		   if(noguess(j)) then
			goto 98				!try next set
		   else
			n1=n1+1
			Y0=thetgues(1,j)	!-guesses for sep fits -see ycalcv
			Ym=thetgues(2,j)
			ak=thetgues(3,j)
			enh=thetgues(4,j)
		   endif
		else
		 n=0
		 if(sepmax) then
		   ymax1=ymaxsep(j)
		 else
		   ymax1=ymax
		 endif
		 do i=1,nj(j)
		   if(Yobs(i,j).gt.0.1*ymax1.and.Yobs(i,j).lt.0.9*Ymax1)then
			p=Yobs(i,j)/(ymax1-Yobs(i,j))
			n=n+1
		      x1(n)=alog(Xobs(i,j))
		      y1(n)=alog(p)
		   endif
		 enddo
		 if(n.lt.2) then
		   noguess(j)=.true.		!or try another method?
		   goto 98				!try next set
c		   goto 99
		 endif
		 call LINREG(x1,y1,n,a,b,0)
		 n1=n1+1
		 enh=b			!nH
		 if(iequiv.eq.1) enh=-b	!decreasing curve
		 aK=exp(-a/b)	!aK
		 Ym=1.05*ymax1
		 thetgues(1,j)=Y0	!-guesses for sep fits -see ycalcv
		 thetgues(2,j)=Ym
		 thetgues(3,j)=aK
		 thetgues(4,j)=enH
		endif
		ymm=ymm+Ym		!for mean ymax
		akm=akm+aK		!for mean K
		enm=enm+enH		!for mean nH
c Could now do modified y vs y/x plot using the value of nH just found, to
c re-estimate K and get proper estimate of max? -or try just using max obs
c value (*1.2?) for now
98		continue	!jump to next set if noguess(j) true
	   enddo
c===change 01/27/97 05:43pm
	   if(n1.ge.1) then
		ymm=ymm/float(n1)	!mean ymax
		akm=akm/float(n1)	!mean K
		enm=enm/float(n1)	!mean nH
	   else
		ymm=ymax1
		enm=1.3	!can't estimate so use 1.3 as guess (K set below)
	   endif
c===end of change 01/27/97 05:43pm
c  Check for sets with only n=<2 so noguess(j)=true
c  For ifitmode=4,5,take one point (nearest to 50% if there are two) and
c using mean ymax, nH, calculate estimate of K from it (only works if curves
c all have similar ymax which may not be true, esp for ifitmode=2) -for now
c use sepmax to judge this (could also do it only for fitmode=4,5?)
	   do m=1,nsfit
	      j=juse(m)		!set # used
	      if(noguess(j)) then
		   if(.not.sepmax) then
			noguess(j)=.false.
c===change 01/27/97 05:43pm
			d1min=1.e37
			do i=1,nj(j)
			   d1=abs(0.5-Yobs(i,j)/ymm)
			   if(d1.lt.d1min) then
				d1min=d1
				imin=i
			   endif
			enddo
			y=Yobs(imin,j)
			x=Xobs(imin,j)
			aK=x*(ymm/y - 1.)**(1./enm)
			thetgues(1,j)=Y0	!-guesses for sep fits -see ycalcv
		      thetgues(2,j)=Ymm	!mean
			thetgues(3,j)=aK	!from above
			thetgues(4,j)=enm	!mean
		   endif
		endif
	   enddo
c===end of change 01/27/97 05:43pm
c
	   if(ifitmode.le.2) then
		theta(1)=Y0
		theta(2)=thetgues(2,jset)	!Ym for 1st set to be fitted, #jset
		theta(3)=thetgues(3,jset)	!K for 1st set to be fitted, #jset
		theta(4)=thetgues(4,jset)	!nH for 1st set to be fitted, #jset
	   else if(ifitmode.eq.3) then
		if(ip1.eq.-1) then		!sep ymax, sep K
		   theta(1)=Y0
		   k=1
		   do m=1,nsfit
			j=juse(m)		!set # used
			k=k+1
			theta(k)=thetgues(3,j)  !K value for each set
			k=k+1
			theta(k)=thetgues(2,j)  !Ym value for each set
		   enddo
		   if(nmod.eq.27) then 	!add nH
			if(ip2.eq.1) then		!common nH
			   k=k+1
			   theta(k)=enm		!mean nH
			else if(ip2.eq.0) then		!nH for each set
			   do m=1,nsfit
			      j=juse(m)		!set # used
			      k=k+1
				theta(k)=thetgues(4,j)  !nH value for each set
			   enddo
			endif
		   endif
		else if(ip1.eq.0) then		!common ymax, sep K
		   theta(1)=Y0
		   theta(2)=ymm 	!mean Ymax for common Ymax guess
		   k=2
		   do m=1,nsfit
			j=juse(m)		!set # used
			k=k+1
			theta(k)=thetgues(3,j)  !K value for each set
		   enddo
		   if(nmod.eq.27) then 	!add nH
			if(ip2.eq.1) then		!common nH
			   k=k+1
			   theta(k)=enm		!mean nH
			else if(ip2.eq.0) then		!nH for each set
			   do m=1,nsfit
			      j=juse(m)		!set # used
			      k=k+1
				theta(k)=thetgues(4,j)  !nH value for each set
			   enddo
			endif
		   endif
		else if(ip1.eq.1) then	!common K, sep ymax
		   theta(1)=Y0
		   k=1
		   do m=1,nsfit
			j=juse(m)		!set # used
			k=k+1
			theta(k)=thetgues(2,j)  !Ym value for each set
		   enddo
		   k=k+1
		   theta(k)=akm		!mean K for common K guess
		   if(nmod.eq.27) then 	!add nH
			if(ip2.eq.1) then		!common nH
			   k=k+1
			   theta(k)=enm		!one nH value
			else
			   do m=1,nsfit
				j=juse(m)		!set # used
				k=k+1
				theta(k)=thetgues(4,j)  !nH value for each set
			   enddo
			endif
		   endif
c
		   if(logyfit) theta(1)=theta(k)	!only one (see geteqn)
		endif
	   else if(ifitmode.eq.4) then
		theta(1)=Y0
		theta(2)=thetgues(2,jfirst)	!Ym for 1st set to be fitted
		theta(3)=thetgues(3,jfirst)	!K for 1st set to be fitted
		theta(4)=thetgues(4,jfirst)	!nH for 1st set to be fitted
		k=4
		do m=2,nsfit
		   j=juse(m)
		   k=k+1
	         theta(k)=thetgues(3,j)/thetgues(3,jfirst)	!approx for pot ratios
		enddo
	   else if(ifitmode.eq.5) then
		theta(1)=Y0
		theta(2)=thetgues(2,jfirst)	!Ym for 1st set to be fitted
		theta(3)=thetgues(3,jfirst)	!K for 1st set to be fitted
		theta(4)=thetgues(4,jfirst)	!nH for 1st set to be fitted
		do m=2,nsfit
		   j=juse(m)		!set # used
	         r=thetgues(3,j)/thetgues(3,jfirst)	!approx for pot ratio
		   xB=setx(j)
		   if(r.gt.1.2) then
			n=n+1
			s=s+xB/(r-1.)
		   endif
		enddo
		theta(5)=s/float(n)	!estimate of KB
	   endif
c
c Power function fit
	else if(nmod.eq.28) then
	   n=0
	   sy=0.0	!overall mean y value for scaling
	   do m=1,nsfit
		j=juse(m)		!set # used
		do i=1,nj(j)
		   n=n+1
		   sy=sy + Yobs(i,j)
		enddo
	   enddo
	   ybarp=sy/float(n)		!for ycalcv via common
c
	   bbar=0.0
	   n=0
	   do m=1,nsfit
	     j=juse(m)		!set # used
	     if(nj(j).gt.1) then
		sy1=0.0
		sx1=0.0
		sxy1=0.0
		sxx1=0.0
		do i=1,nj(j)
		   x=alog10(Xobs(i,j))
		   y=alog10(Yobs(i,j))
		   sy1=sy1 + y
		   sx1=sx1 + x
		   sxy1=sxy1 + x*y
		   sxx1=sxx1 + x*x
		enddo
		en=float(nj(j))
		b=(sxy1-sx1*sy1/en)/(sxx1-sx1*sx1/en)	!slope of log-log plot
		ym1=sy1/en
		xm1=sx1/en
		x0=(alog10(ybarp)-ym1+b*xm1)/b
		x0=10**(x0)
		thetgues(1,j)=x0	!see ycalcv
		thetgues(2,j)=b
		if(norm) thetgues(1,j)=b
		n=n+1
		bbar=bbar+b
	     endif
	   enddo
	   bbar=bbar/float(n)	!mean slope
c  Check for sets with only one observation (i=1)
	   do m=1,nsfit
	      j=juse(m)		!set # used
		i=1
	      if(nj(j).eq.1) then
		   x=alog10(Xobs(i,j))
		   y=alog10(Yobs(i,j))
		   yb=alog10(ybarp)
		   x=x+(yb-y)/bbar
		   thetgues(1,j)=10**x
		   thetgues(2,j)=bbar	!use mean slope
		endif
	   enddo
c  Set theta for the fitmode
	   if(ifitmode.le.2) then
	      theta(1)=thetgues(1,jset)		!x0 for 1st set
		theta(2)=thetgues(2,jset)	!n for 1st set to be fitted, #jset
	   else if(ifitmode.eq.3) then
c/=====to be fixed
	   else if(ifitmode.eq.4) then
c  Now set the extra parameter guesses for the simultaneous fit in theta()
c
		jfirst=juse(1)
		theta(1)=thetgues(1,jfirst)		!x0 for 1st set
		theta(2)=bbar				!mean slope
		k=2
		do m=2,nsfit
		   j=juse(m)
		   k=k+1
	         theta(k)=thetgues(1,j)/thetgues(1,jfirst)	!approx for pot ratio
		enddo
	   else if(ifitmode.eq.5) then
	      theta(1)=thetgues(1,jfirst)		!x0 for 1st set
		theta(2)=thetgues(2,jfirst)	!n for 1st set to be fitted, #jset
		do m=2,nsfit
		   j=juse(m)		!set # used
	         r=thetgues(1,j)/thetgues(1,jfirst)	!approx for pot ratio
		   xB=setx(j)
		   if(r.gt.1.2) then
			n=n+1
			s=s+xB/(r-1.)
		   endif
		enddo
		theta(3)=s/float(n)	!estimate of KB
	   endif
c Straight line fit (nmod=99 and ip1=1)
	else if(nmod.eq.99.and.ip1.eq.1) then
	   do m=1,nsfit
		j=juse(m)		!set # used
		if(guesdone) then
		   if(noguess(j)) then
			goto 96				!try next set
		   else
			n1=n1+1
			Y0=thetgues(1,j)	!-guesses for sep fits -see ycalcv
			slope=thetgues(2,j)
		   endif
		else
		   do i=1,nj(j)
			x1(i)=Xobs(i,j)
			y1(i)=Yobs(i,j)
		   enddo
		   n=nj(j)
		   if(n.lt.2) then
			noguess(j)=.true.		!or try another method?
			goto 96				!try next set
		   endif
		   if(ip2.eq.-1) then	!Y=1+x/KB
			a=1.0
		      call LINREG(x1,y1,nj(j),a,b,1)
			thetgues(1,j)=1.0/b	!-guesses for sep fits -see ycalcv
		   else if(ip2.eq.-2) then	!fit Y=x/KB
			a=0.0
		      call LINREG(x1,y1,nj(j),a,b,1)
			thetgues(1,j)=1.0/b	!-guesses for sep fits -see ycalcv
		   else				!fit Y=a+b*x
		      call LINREG(x1,y1,nj(j),a,b,0)
			Y0=a
			slope=b
			thetgues(1,j)=Y0	!-guesses for sep fits -see ycalcv
			thetgues(2,j)=slope
		   endif
		endif
96		continue
	   enddo
c
	else
	   do j=1,nset1
		noguess(j)=.true.		!guesses not defined here
	   enddo
	endif
c
99	DEALLOCATE(x1,y1)
	RETURN
	end

	subroutine LINREG(x,y,n,a,b,imode)
c To return a,b for conventional linear regression
c imode=0 for normal fit
c imode=1 for fit of slope with a=input value (e.g. a=0 or a=1 for Schild eq)
	real*4 x(n),y(n)
c
	en=float(n)
	sx=0.
	sxx=0.
	sy=0.
	syy=0.
	sxy=0.
	do i=1,n
	   sx=sx+x(i)
	   sxx=sxx+x(i)*x(i)
	   sy=sy+y(i)
	   syy=syy+y(i)*y(i)
	   sxy=sxy+x(i)*y(i)
	enddo
	sxx=sxx - (sx*sx)/en
	sxy=sxy - (sx*sy)/en
c
	b=sxy/sxx
	xbar=sx/en
	if(imode.eq.0) a=(sy/en)-b*xbar
	RETURN
	end


