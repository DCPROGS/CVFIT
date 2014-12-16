	recursive function YCALCV(K,THETA,xv1,ival,j,setx,njset)
C
C CALCULATES I'TH VALUE IN J'TH DATA SET FOR MODEL #NMOD IN CVFIT
c
c With iver=1003 have not got allocatable arrays in common, so
c need to add setx() to parameters
c Xobs cannot be in common now.  Rather than making it a parameter, altered so
c that xv1 is parameter - YCALCJ will use x=xv (from common/block1) if ival<1
c in call, otherwise it will use x=xv1 (parameter) which has been set  to
c xobs(i,j) before call to Ycalcv
c
c Modif 02/01/98 07:11am New option for Hill/Langmuir (ip1=2 for this option)
c Only for ifitmode=3, ncomp=2 -common Ymax1/ymax2 -in this case fit K,
c    Ymax (and nH if req), one value for each set (not separate Ymax1, Ymax2)
c    plus a single value of Ymax1/Ymax for all sets
c Fit percent comp 1, say P1 = Ymax1/(Ymax1+Ymax2) = Ymax1/Ymax.  Then to
c get Ymax1 and Ymax2 (eg in ycalcv) for each set, find
c Ymax1 = P1*Ymax, then Ymax2=Ymax-Ymax1
c
c 12/22/96 02:22pm Power function (nmod=28) altered to have form
c Y=ybar*(x/x0)**n, so Y=ybar when x=x0; adjustable param are x0,n, and
c  ybar is mean of all fitted y values (in common) to give sensible scale factor
c 03/20/96 10:01am Add nmod=28 (imod=23) for power function
c
c Modif 09/30/95 06:51am disjoint polynomial for join-point (x1)
c other than zero by making the polynomial a function of (xv-x1) rather
c than of xv (see Mathcad6 file disjoint.mcd)
c
c Modif 05/05/95 09:20am by adding option for nmod=26,27 to
c (a) normalise all Y values relative to a common conc Y(xnorm)
c (b) in case of common K, sep max, to use log(y) values for fitting as
c suggested by AGH so sep maxima are linear and can be found more easily
c -extra values kept in COMMON/BLOCK3/ logyfit,norm,xnorm
c
c Modif 03/28/95 09:52am Add nmod=26,27 (Langmuir/Hill resp) with ncomp
c components and common param defined via ip1,ip2 thus:
c ip1=0 -common max, separate K
c ip1=1 -common K, separate max
c ip2=0 -separate nH
c ip2=1 -common nH
c Add ival=2 to signify that, for Langmuir/Hill (nmod=26,27) amplitudes
c of components are to be normalised (do this only if ymax(j) defined!)
c
c 09/17/94 11:38am Add nmod=25 =straight lines with common intercept on Y axis
c General exponentials (nmod=19) added 04/25/91 09:57am
c
c Oct 86- two new variables added to common/block2/
c IEQUIV=0 for equiv subunits; =1 for non-equiv (prev signalled via castar)
c IP1,IP2:
c  (1) mods 10,11 IP1=1 if param are K1,K2; =2 if param are K1, K2/K1.
c  (2) mods 9,11  IP2=1 if param is K2 or (K2/K1 if IP1=2)
c		  IP2=2 if param is k(-2) (must have ip2=1 if ip1=2)
c  (3) mods 9,11  IP2 set negative if k(+2) is a param  to be estimated
c  (4) nmod=99 =polynomial. Degree of poly=ip1 (and if ip2>0 then degree
c	is ip2 for x>x1 )
c	If ip2=-1 then Schild eqn is fitted with data Y=dose ratio, X=antag conc
c	If ip2=-2 then Schild eqn is fitted with data Y=dose ratio-1, X=antag conc
c
	real THETA(100)
	LOGICAL CONSTR,FLINE
c=	common/data/Xobs(100,10),Yobs(100,10)
c=	COMMON/SD/w(100,10)
	real*4 setx(njset)
	real*4 ymaxj(100)			!in common/max
c====nj(),jmiss(),juse() no longer in common/block1/
	COMMON/BLOCK1/constr,nset,nfit,nsfit,
     & Xv2,alfa,kmax,ncomp,nmod,fline,nomit,jomit(20),jset,ifitmode
c=	COMMON/BLOCK2/ castar,setx(10),X1,X2,iequiv,ip1,ip2
c====setx no longer in common/block2/
	COMMON/BLOCK2/ castar,X1,X2,iequiv,ip1,ip2
	logical logyfit,norm
	COMMON/BLOCK3/logyfit,norm,xnorm,iostat
	common/potrat/jfirst,iset,kmax1			!for pot ratios etc
c (kmax=actual no of param to be estimated; kmax1=number before extra param
c  added for ifitmode=4,5)
	common/max/ymaxj
	common/pwrfunc/ybarp
	common/binomial/binom(0:10)
c
	call UNDER0('true')	!so underflows set=0 (Lahey manual p12.8)
	ycalcv=0.0
	iostat=0
c	if(i.le.100) XV=Xobs(I,J)
	if(ival.ge.1) then
	   xv=xv1	!parameter
	else
	   xv=xv2	!from common
	endif
c If ival<1 then value of XV in common is used
c
c	if(norm) then
c -calc response at x=xnorm -don't want xv changed after recursive call, so
c if ival=1 for 1st call, make use common for recursive call and vice versa
c Recursive call is too complicated -the call works OK as below, but
c the parameters are numbered differently when norm=false!
c  Better get Ynorm the long way -but need Ymax estimate for this!
c	   ivsav=ival
c	   norm=.false.	!set false for recursive call
c	   if(ival.eq.1) then
c	      ival1=0		!so xvalue in common used
c	      xv2=xnorm	!xV in common
c	      ynorm=YCALCV(K,theta,xvdum,ival1,j,setx,njset)
c	   else
c		ival1=1
c	      ynorm=YCALCV(K,theta,xnorm,ival1,j,setx,njset)
c	   endif
c	   norm=.true.	!reset
c	   ival=ivsav	!reset
c	endif
c
	if(ifitmode.eq.4.and.iset.gt.1) then	!i.e. not for the first set
c	   ratio=theta(kmax-nsfit+iset)
	   ratio=theta(kmax1+iset-1)		!same thing
	   xv=xv/ratio		!divide conc actually used by r1,r2,...
	else if(ifitmode.eq.5) then
	   akb=theta(kmax)	!=KB (see GETEQN)
	   xB=setx(j)		!antag conc
	   xv=xv/(1.0 + xB/aKB)	!mult by dose ratio
	endif
c
c=======temp stuff for Woda problem
	if(nmod.eq.100) then
	   a=theta(1)
	   b=theta(2)
	   c=theta(3)
	   d=0.67
	   ycalcv=a*(1.0-exp(-((xv+b)**d)/c))
	   RETURN
	endif
c=======end of temp
c New versions for Hill/Langmuir nmod=26,27
c theta(1)=Y(0) always.  Then have j2=(ncomp+ncomp*nsfit) parameters for
c the K and Ymax values, so nH values (for nmod=27) start at theta(j2+1)
c NB ymax values for each component are values measured from Y(0)=baseline
c (which are those needed to calc relative amps of components, and
c for normalising)
c When data are normalised relative to Y(xnorm) then Ycal(x) must be
c divided by Ycal(xnorm); when Y(0)=0 and ncomp=1 then Ymax cancels
c and is therefore indeterminate -so Ymax parameter omitted in this case
c (see GETEQN)
	if(nmod.eq.26.or.nmod.eq.27) then
	   if(norm) then
		j0=1+(iset-1)*ncomp	!=1,2,3 ... since ncomp=1 if norm
		j1=nsfit*ncomp+1	!# of params for Y0 and all Ymax values
		ycalx=0.
		ynorm=0.
		Y0=theta(1)
		j0=j0+1
		ak=theta(j0)
	 	c=xv/aK
		cnorm=xnorm/aK
		if(nmod.eq.27) then
		   hill=theta(kmax)
		   c=c**hill
		   cnorm=cnorm**hill
		endif
		if(iequiv.eq.0) then
		   p=c/(c+1.0)	!increasing (normal) curve
		   pnorm=cnorm/(cnorm+1.0)	!increasing (normal) curve
		else
		   p=1.0/(1.0 + c)	!decreasing curve
		   pnorm=1.0/(1.0 + cnorm)	!decreasing curve
		endif
	   	ycalx=ycalx + p
	   	ynorm=ynorm + pnorm
	      ycalcv=ycalx/ynorm
		RETURN
	   endif
c
	   if(logyfit) then	!only for y(0)=0, ncomp=1
		aK=theta(1)
		ymax=theta(j+2)	!ymax for current set, calc in SSDCV
		if(ival.eq.2) ymax=1.0	!for normalised curve
	 	c=xv/aK
		if(nmod.eq.27) then
		  hill=theta(2)
		  c=c**hill
		endif
		if(iequiv.eq.0) then
		   x=c/(c+1.0)	!increasing (normal) curve
		else
		   x=1.0/(1.0 + c)	!decreasing curve
		endif
	   	ycalcv=Ymax*x
		RETURN
	   endif
c
c====	   j2=1+ncomp*(1+nsfit)
	   j2=ncomp*(1+nsfit)
	   Y0=theta(1)
c	   if(ifitmode.le.2) then	!no common params
c  Now ifitmode=4,5 calculated here too in these cases need only first 3/4
c parameters -pot ratio/KB taken into account above
	   if(ifitmode.ne.3) then	!no common params
		ycalx=0.
		j0=1
		do i=1,ncomp
		   j0=j0+1
		   ymax=theta(j0)		!Ymax for each component=theta(2 to ncomp+1)
		   j0=j0+1
		   ak=theta(j0)
	 	   c=xv/aK
		   if(nmod.eq.27) then
		      j0=j0+1
			hill=theta(j0)
			c=c**hill
		   endif
		   if(iequiv.eq.0) then
			x=c/(c+1.0)	!increasing (normal) curve
		   else
			x=1.0/(1.0 + c)	!decreasing curve
		   endif
	   	   ycalx=ycalx + Ymax*x
		enddo
		ycalx=Y0 + ycalx
		if(fline) then
		   slope=theta(kmax)
		   ycalx=ycalx + slope*xv
		endif
	   else if(ifitmode.eq.3.and.ip1.eq.-1) then !sep max, sep K for each set
		j0=1+(iset-1)*2*ncomp	!=1,1+2ncomp,1+4*ncomp for set 1,2,3...
		j2=2+2*nsfit*ncomp	!index for first nH

		ycalx=0.
		jk=j0
		jm=0
		do i=1,ncomp
		   jk=jk+1
		   aK=theta(jk)
		   if(jk.gt.jm) jm=jk
		   jk=jk+1
		   ymax=theta(jk)		!Ymax for each component
		   if(jk.gt.jm) jm=jk
	 	   c=xv/aK
		   if(nmod.eq.27) then
			if(ip2.eq.1) then		!Common nH for all sets
			   m=j2+i-1
			   hill=theta(m)
			   if(m.gt.jm) jm=m
			else
			   m=j2+(iset-1)*2*ncomp+i-1
			   hill=theta(m)
			   if(m.gt.jm) jm=m
			endif
			if(hill*alog10(c).gt.37.or.hill*alog10(c).lt.-37) then
			   iostat=1
			   goto 999
			endif
			c=c**hill
		   endif
		   if(iequiv.eq.0) then
			x=c/(c+1.0)	!increasing (normal) curve
		   else
			x=1.0/(1.0 + c)	!decreasing curve
		   endif
	   	   ycalx=ycalx + Ymax*x
		enddo
		ycalx=Y0 + ycalx
		if(fline) then
		   m1=jm+iset
		   m=kmax-nsfit+iset
		   if(m.ne.m1) then
			call BELL(1)
			print 70,m,m1,iset
70			format(' Index error for slope: m,m1,iset = ',3i5)
		   endif
		   slope=theta(m)
		   ycalx=ycalx + slope*xv
		endif
	   else if(ifitmode.eq.3.and.ip1.eq.0) then !common max, sep K for each set
		j0=1+(iset-1)*ncomp	!=1,1+ncomp,1+2*ncomp for set 1,2,3...
		ycalx=0.
		jm=0
		do i=1,ncomp
		   ymax=theta(i+1)		!Ymax for each component=theta(2 to ncomp+1)
		   m=ncomp+i+j0
		   aK=theta(m)		!K for component #i, set # iset
		   if(m.gt.jm) jm=m	!track max index of theta used so far
	 	   c=xv/aK
		   if(nmod.eq.27) then
			if(ip2.eq.1) then
			   m=j2+1+i
			   hill=theta(m)
			   if(k.gt.jm) jm=m
			else
			   m=j2+j0+i
			   hill=theta(m)
			   if(k.gt.jm) jm=m
			endif
			if(hill*alog10(c).gt.37.or.hill*alog10(c).lt.-37) then
			   iostat=1
			   goto 999
			endif
			c=c**hill
		   endif
		   if(iequiv.eq.0) then
			x=c/(c+1.0)	!increasing (normal) curve
		   else
			x=1.0/(1.0 + c)	!decreasing curve
		   endif
	   	   ycalx=ycalx + Ymax*x
		enddo
		ycalx=Y0 + ycalx
		if(fline) then
		   m1=jm+iset
		   m=kmax-nsfit+iset
		   if(m.ne.m1) then
			call BELL(1)
			print 70,m,m1,iset
c70			format(' Index error for slope: m,m1,iset = ',3i5)
		   endif
		   slope=theta(m)
		   ycalx=ycalx + slope*xv
		endif
c rest is for ifitmode=3
	   else if(ifitmode.eq.3.and.ip1.eq.1) then  !common K, sep Ymax for each set
		j0=1+(iset-1)*ncomp	!=1,1+ncomp,1+2*ncomp for set 1,2,3...
		j1=nsfit*ncomp+1	!# of params for Y0 and all Ymax values
		ycalx=0.
		do i=1,ncomp
		   ymax=theta(j0+i)	!Ymax for component #1, set #iset
c===		   if(ival.eq.2) ymax=ymax/ymaxj(iset)
		   if(ival.eq.2) ymax=ymax/ymaxj(j)
		   aK=theta(i+j1)		!K for component #i
	 	   c=xv/aK
		   if(nmod.eq.27) then
			if(ip2.eq.1) then
			   hill=theta(j2+i+1)
			else
			   hill=theta(j2+j0+i)
			endif
			c=c**hill
		   endif
		   if(iequiv.eq.0) then
			x=c/(c+1.0)	!increasing (normal) curve
		   else
			x=1.0/(c+1.0)	!decreasing curve
		   endif
	   	   ycalx=ycalx + Ymax*x
		enddo
		ycalx=Y0 + ycalx
		if(fline) then
		   m=kmax-nsfit+iset
		   slope=theta(m)
		   ycalx=ycalx + slope*xv
		endif
	   else if(ifitmode.eq.3.and.ncomp.eq.2.and.
     &		(ip1.eq.2.or.ip1.eq.3)) then
c sep K, Ymax common percent comp 1
c  NB this is only for ncomp=2 at present so no need for general ncomp loop
c params are (1) Y(0) (or Y(inf))
c            (2 to nsfit+1)  Ymax(set #1) to Ymax(set #nsfit)   -nsfit values
c            (nsfit+2 to nsfit+3) K1, K2 for set 1
c            (nsfit+4 to nsfit+5) K1, K2 for set 2
c            so K1, for set #iset are params 2+iset*nsfit and K2 is the next
		ycalx=0.
		ymax=theta(1+iset)	!Ymax for set #iset
		if(ival.eq.2) ymax=ymax/ymaxj(j)	!check this!
		if(ip1.eq.2) then
		   rel1=theta(kmax)		!fraction of component 1
		else if(ip1.eq.3) then
		   rel1=1.0-theta(kmax)		!theta=fraction of component 2
		endif
c          NB constrain rel1 between 0 and 1
		if(rel1.lt.1.e-20) then
		   rel1=1.e-20
		   theta(kmax)=rel1
		   if(ip1.eq.3) theta(kmax)=1.0
		else if(rel1.gt.1.0) then
		   rel1=1.0
		   theta(kmax)=rel1
		   if(ip1.eq.3) theta(kmax)=1.e-20
		endif
		ymax1=rel1*ymax
		ymax2=ymax-ymax1
		j1=2 + nsfit + 2*(iset-1)
		aK1=theta(j1)
		aK2=theta(j1+1)
		if(nmod.eq.27) then
		   j2=2 + 3*nsfit 	!index of the first nH
		   if(ip2.eq.1) then
			hill1=theta(j2)	!common nH1
			hill2=theta(j2+1)	!common nH2

		   else
			j2=j2 + 2*(iset-1)
			hill1=theta(j2)
			hill2=theta(j2+1)
		   endif
		endif
		c1=xv/aK1
		c2=xv/aK2
		if(nmod.eq.27) then
		   c1=c1**hill1
		   c2=c2**hill2
		endif
		if(iequiv.eq.0) then
		   x=Ymax1*c1/(c1+1.0) + Ymax2*c2/(c2+1.0) !increasing (normal) curve
		else
		   x=Ymax1/(c1+1.0) + Ymax2/(c2+1.0) !decreasing curve
		endif
c		ycalx=Y0 + ycalx
		ycalx=Y0 + x
		if(fline) then
		   m=kmax-nsfit+iset
		   slope=theta(m)
		   ycalx=ycalx + slope*xv
		endif
	   endif
c All options done -define ycalcv
	   ycalcv=ycalx
c
	   RETURN
	endif		!end of Hill/Langmuir section
c
c n bindings + opening
	if(nmod.eq.30) then
	   nbind=ip1
	   if(ifitmode.le.2) then
		if(iequiv.eq.0) then
		   j1=3*iset
		   aK1=theta(j1-2)		!1,4,7 for set 1,2,3...
		   if(aK1.lt.1.e-10) then
			aK1=1.e-10
			theta(j1-2)=aK1
		   endif
		   E=theta(j1-1)		!2,5,8 for set 1,2,3
		   R=1.0
		   ymax=theta(j1)		!3,6,9
		else if(iequiv.eq.1) then
		   j1=4*iset
		   aK1=theta(j1-3)		!1,5,9 for set 1,2,3...
		   if(aK1.lt.1.e-10) then
			aK1=1.e-10
			theta(j1-3)=aK1
		   endif
		   E=theta(j1-2)		!2,6,10 for set 1,2,3
		   R=theta(j1-1)
		   ymax=theta(j1)		!3,7,11
		endif
	   else if(ifitmode.eq.3) then
		if(iequiv.eq.0) then
		   j1=2*iset
		   aK1=theta(j1-1)		!1,3,5 for set 1,2,3...
c========need lower lim for c (K=1.e-20 causes overflow!)
		   if(aK1.lt.1.e-10) then
			aK1=1.e-10
			theta(j1-1)=aK1
		   endif
		   E=theta(j1)		!2,4,6 for set 1,2,3
		   R=1.0
		   ymax=theta(kmax)
		else if(iequiv.eq.1) then
		   j1=3*iset
		   aK1=theta(j1-2)		!1,4,7 for set 1,2,3...
		   if(aK1.lt.1.e-10) then
			aK1=1.e-10
			theta(j1-2)=aK1
		   else if(aK1.gt.1.e10) then
			aK1=1.e10
			theta(j1-2)=aK1
		   endif
		   E=theta(j1-1)		!2,5,8 for set 1,2,3
		   R=theta(j1)
		   if(R.lt.1.e-10) then
			R=1.e-10
			theta(j1)=R
		   else if(R.gt.1.e10) then
			R=1.e10
			theta(j1)=R
		   endif
		   ymax=theta(kmax)
		endif
	   endif
	   den=1.0
	   aK=aK1	!for first binding
	   do i=1,nbind
		term=binom(i)*(xV**i)/aK
	      den=den+term
		termlast=term
		aK=aK*aK1/(R**i)	!=K1, K1*K2, K1*K2*K3 etc where K2=K1/R etc
	   enddo
	   a=E*termlast
	   ycalcv=ymax*a/(den+a)
	   RETURN
	endif
c MWC(n)
	if(nmod.eq.31) then
	   nbind=ip1
	   an=float(nbind)
c	   if(ifitmode.le.2) then
c		j1=4*iset
c		aK=theta(j1-3)		!1,5,8 for set 1,2,3...
cc========need lower lim for c (K=1.e-20 causes overflow!)
c		if(aK.lt.1.e-10) then
c		   aK=1.e-10
c		   theta(j1-3)=aK
c		endif
c		R=theta(j1-2)		!2,6,10 for set 1,2,3
c		if(fline.and.R.lt.1.) then	!constrain R>1
c		   R=1.001
c		   theta(j1-2)=R
c		endif
c		En=theta(j1-1)		!3,7,11 for set 1,2,3
c		Rn=R**an
c		if(En.gt.1.e-36*Rn) then	!check for underflow of E0
c		   E0=En/Rn
c		else
c		   E0=1.e-36
c		   En=E0*Rn
c		   theta(j1)=En
c		endif
c		ymax=theta(j1)		!4,8,12
c	   else if(ifitmode.eq.3) then
	   if(ifitmode.le.3) then
		if(iequiv.le.0) then	!fit all
		   j1=3*iset
		   aK=theta(j1-2)		!1,4,7 for set 1,2,3...
		   if(aK.lt.1.e-10) then
			aK=1.e-10
			theta(j1-1)=aK
		   endif
		   R=theta(j1-1)		!2,5,8 for set 1,2,3
		   if(fline.and.R.lt.1.) then	!constrain R>1
			R=1.001
			theta(j1-1)=R
		   endif
		   if(iequiv.eq.0) then
			En=theta(j1)			!3,6,9 for set 1,2,3
			Rn=R**an
			if(En.gt.1.e-36*Rn) then	!check for underflow of E0
		         E0=En/(Rn)
			else
			   E0=1.e-36
			   En=E0*Rn
			   theta(j1)=En
			endif
		   else if(iequiv.eq.-1) then
			E0=theta(j1)			!3,6,9 for set 1,2,3
			if(E0.lt.1.e-25) then
		         E0=1.e-25
			   theta(j1)=E0
			endif
		   endif
		   ymax=theta(kmax)
		else if(iequiv.eq.1) then		!common E0
		   j1=2*iset
		   aK=theta(j1-1)		!1,3,5 for set 1,2,3...
		   if(aK.lt.1.e-10) then
			aK=1.e-10
			theta(j1-1)=aK
		   endif
		   R=theta(j1)
		   if(fline.and.R.lt.1.) then	!constrain R>1
			R=1.001
			theta(j1)=R
	         endif
		   E0=theta(kmax-1)
		   ymax=theta(kmax)
		else if(iequiv.eq.2) then		!common R
		   j1=2*iset
		   aK=theta(j1-1)		!1,3,5 for set 1,2,3...
		   if(aK.lt.1.e-10) then
			aK=1.e-10
			theta(j1-1)=aK
		   endif
		   En=theta(j1)
		   R=theta(kmax-1)
		   if(fline.and.R.lt.1.) then	!constrain R>1
			R=1.001
			theta(kmax-1)=R
	         endif
		   Rn=R**an
		   if(En.gt.1.e-36*Rn) then	!check for underflow of E0
		      E0=En/(Rn)
		   else
			E0=1.e-36
			En=E0*Rn
			theta(j1)=En
		   endif
		   ymax=theta(kmax)
		else if(iequiv.eq.3) then		!all common except E0
		   aK=theta(1)
		   if(aK.lt.1.e-10) then
			aK=1.e-10
			theta(1)=aK
		   endif
		   R=theta(2)
		   if(fline.and.R.lt.1.) then	!constrain R>1
			R=1.001
			theta(2)=R
	         endif
		   j1=iset+2
		   E0=theta(j1)
		   ymax=theta(kmax)
		else if(iequiv.eq.4) then		!all common except R=K/K*
		   aK=theta(1)
		   if(aK.lt.1.e-10) then
			aK=1.e-10
			theta(1)=aK
		   endif
		   j1=iset+1
		   R=theta(j1)
		   if(fline.and.R.lt.1.) then	!constrain R>1
			R=1.001
			theta(j1)=R
	         endif
		   E0=theta(kmax-1)
		   ymax=theta(kmax)
		endif
	   endif
	   c=xV/aK
	   cstar=R*c
	   a=E0*(1 + cstar)**nbind
	   a1=(1 + c)**nbind
	   ycalcv=ymax*a/(a1+a)
	   RETURN
	endif		!end of nmod=31
c
c General exponentials (nmod=19)
c amps= theta 2,4,6,8,10; theta(2*i-1) tau=theta 1,3,5,7,9
	if(nmod.eq.19) then
c	  ncomp=kmax/2		!ncomp now in common
	  ycalcv=0.
	  do 53 i=1,ncomp
	   i1=2*i
	   tau=theta(i1-1)
	   amp=theta(i1)
	   ycalcv=ycalcv + amp*exp(-xv/tau)
53	  continue
c	  ycalcv=ycalcv+theta(kmax)		!add Y(inf)
c	  ycalcv=ycalcv+theta(2*ncomp+1)	!add Y(inf)
	  ycalcv=ycalcv+theta(kmax1)		!add Y(inf)
	  RETURN
	endif
c General geometrics (nmod=20)
c amps= theta 2,4,6,8,10; theta(2*i-1) mu=theta 1,3,5,7,9
c mu=1/(1-lambda) ='time constant'
	if(nmod.eq.20) then
	  ycalcv=0.
	  do i=1,ncomp
	   i1=2*i
	   amu=theta(i1-1)
	   if(amu.lt.1.0) then
		amu=1.0001	!lambda must be non-negative
	      theta(i1-1)=amu
	   endif
	   amp=theta(i1)
	   alambda=1.0 -(1.0/amu)
	   ycalcv=ycalcv + amp*(alambda**(xv))
	  enddo
	  ycalcv=ycalcv+theta(kmax1)		!add Y(inf)
	  RETURN
	endif
c
c Lines with common Y-intercept (ith value,jth data set=jth slope)
	if(nmod.eq.25) then
	   slope=theta(iset)
	   xint=theta(kmax)	!intercept on x axis
	   aj=-slope*xint		!intercept on Y axis
	   ycalcv=aj + slope*xv
	   RETURN
	endif
c
	if(nmod.eq.28) then	!power function
	   if(norm) then
		en=theta(1)
		ycalcv=(xv/xnorm)**en
		RETURN
	   endif
	   x0=theta(1)
	   en=theta(2)
	   if(abs(x0).lt.1.e-20) then
		x0=1.e-20
		theta(1)=x0
	   endif
	   ycalcv=(xv/x0)**en
	   ycalcv=ybarp*ycalcv		!ybar is in common
	   RETURN
	endif
c
	if(nmod.eq.16.or.nmod.eq.17.or.nmod.eq.18) goto 329
	if(nmod.eq.21) goto 21
	if(nmod.eq.22) goto 22
c
c polynomial
	if(nmod.ne.99) goto 31
c see NUMERICAL RECIPES for polynomial evaluation
c Also return the slope in X2 (in COMMON)
	if(ip2.gt.0) goto 221		!disjoint polynomial
c Degree of polynomial = ip1 = kmax-1 (but kmax altered for ifitmode=4,5!)
c
c	If ip2=-1 then Schild eqn is fitted with data Y=dose ratio, X=antag conc
c	If ip2=-2 then Schild eqn is fitted with data Y=dose ratio-1, X=antag conc
	if(ip2.eq.-1.or.ip2.eq.-2) then	!fit Y=x/KB or Y=1+x/KB
	   aKB=theta(1)
	   if(ip2.eq.-1) then	!Y=1+x/KB
		ycalcv=1.0 + xv/aKB
	   else if(ip2.eq.-2) then	!fit Y=x/KB
		ycalcv=xv/aKB
	   endif
	   RETURN
c    If ip1=2 (quadratic) and iequiv=1 then fit current-variance parabola
	else if(ip1.eq.2.and.iequiv.eq.1) then
	   cur=theta(1)	!single channel current
	   en=theta(2)
	   ycalcv=cur*xv - xv*xv/en
	   RETURN
	endif
c Now regular polynomial
	x2=0.
	ycalcv=theta(ip1+1)
	do m=ip1,1,-1
	   x2=x2*xv + ycalcv
	   ycalcv=ycalcv*xv + theta(m)
	enddo
	RETURN
221	continue
c For disjoint polynomial degree=ip1 at x<x1; degree=ip2 at x>x1
c Modif 09/30/95 06:51am for x1 other than zero by making the polynomial
c a function of (xv-x1) rather than of xv (see Mathcad6 file disjoint.mcd)
	x2=0.			!initialise for slope
	x3=xv-x1
	if(xv.le.x1) then
	   ycalcv=theta(ip1+1)
 	   do m=ip1,1,-1
		x2=x2*x3+ ycalcv
223		ycalcv=ycalcv*x3 + theta(m)
	   enddo
	else		!xv>x1
	   ycalcv=theta(ip2+1)
 	   do m=ip2,1,-1
		x2=x2*x3 + ycalcv
		ycalcv=ycalcv*x3 + theta(m)
	   enddo
	endif
	RETURN
c
31	if(nmod.le.4) then
c	   NCOMP=KMAX/2		!ncomp now in common
	   IF(NMOD.EQ.-1) ncomp=1
	   AK1=THETA(1)
	   YMAX1=THETA(2)
	   IF(NMOD.EQ.-1) GOTO 11
	   ycalcv=YMAX1*(XV/(XV+AK1))
	   IF(NCOMP.EQ.1) GOTO 10
	   AK2=THETA(3)
	   YMAX2=THETA(4)
	   ycalcv=ycalcv+YMAX2*(XV/(XV+AK2))
c10	   if(fline) ycalcv=ycalcv+theta(kmax)*XV
10	   if(fline) ycalcv=ycalcv+theta(kmax1)*XV
	   RETURN
C HILL EQN
11	   Y0=THETA(3)
	   HILL=THETA(4)
c	   XV=XV**HILL
c	   AK1=AK1**HILL
c	   ycalcv=Y0 + (YMAX1-Y0)*XV/(XV+AK1)
	   c=(xv/aK1)**hill
	   ycalcv=Y0 + (YMAX1-Y0)*c/(c+1.0)
	   RETURN
	endif
	if(nmod.eq.23) then   !single Langmuir per set, with common max
	   y0=theta(1)
	   ymax=theta(2)
	   AK=theta(2+iset)
	   ycalcv=y0 + ymax*xv/(xv+AK)
	   RETURN
	endif
	if(nmod.eq.24) then   !single Hill eqn per set, with common max
	   y0=theta(1)
	   ymax=theta(2)
	   m=2*iset+1
	   AK=theta(m)
	   hill=theta(m+1)
	   c=(xv/AK)**hill
	   ycalcv=y0 + ymax*c/(1.+c)
	   RETURN
	endif
C
C OTHER MODELS (NMOD>4)HERE
300	CONTINUE
	if(nmod.eq.5) then	!binding inhibition with parameter=KB
	   Yinf=theta(1)
	   Y0=theta(2)
	   aKB=theta(3)
	   cB=xv/aKB
	   ycalcv=yinf + y0*(1.0 + castar)/(1.0 + castar + cB)
	   RETURN
	endif
301	IF(NMOD.NE.6) GOTO 302
	AK2=THETA(1)	!DISSOC RATE CONST
	AK1=THETA(2)	!ASSOC RATE CONST
	BINF=THETA(3)	!OVERALL MAX
	X=SETX(J)	!LIGAND CONC FOR SET J
	C=X*AK1/AK2	!NORMALISED CONC
	YMAX=BINF*C/(1.0+C)	!MAX FOR SET J
	T=XV		!TIME
	ycalcv=YMAX*(1.0-EXP(-(AK1*X+AK2)*T))
	RETURN
302	IF(NMOD.NE.7) GOTO 303
	AK2=THETA(1)	!DISSOC RATE CONST
	AK1=THETA(2)	!ASSOC RATE CONST
	BMAX=THETA(2+J)	!MAX FOR SET J
	X=SETX(J)	!LIGAND CONC FOR SET J
	T=XV		!TIME
	ycalcv=BMAX*(1.0-EXP(-(AK1*X+AK2)*T))
	RETURN
303	CONTINUE	!NEXT MODEL
c
C KM2 MODELS
	if(nmod.eq.13) goto 52
	if(nmod.eq.14.or.nmod.eq.15) goto 424
	IF(NMOD.LT.8.OR.NMOD.GT.12) GOTO 304
c Oct 86- two new variables added to common/block2/
c IEQUIV=0 for equiv subunits; =1 for non-equiv (prev signalled via castar)
c  (1) mods 10,11 IP1=1 if param are K1,K2; =2 if param are K1, K2/K1.
c  (2) mods 9,11  IP2=1 if param is K2 or (K2/K1 if IP1=2)
c		  IP2=2 if param is k(-2) (must have ip2=1 if ip1=2)
c  (3) mods 9,11  IP2 set negative if k(+2) is a param  to be estimated
c Value of k(+2) is in X2 if value is fixed, or in theta(kmax) if it
c is to be estimated (ip2 neg)
	if(ip2.lt.0) x2=theta(kmax1)
	ip=iabs(ip2)
	BA=THETA(1)		!B/A FOR N=8,10,12; BETA FOR N=9,11
	if(ip.eq.2) goto 50
c for ip2=+/-1
	AK1=THETA(2)
	AK2=AK1
	IF(ip1.eq.1.and.(NMOD.EQ.10.OR.NMOD.EQ.11)) AK2=THETA(3)
	IF(ip1.eq.2.and.(NMOD.EQ.10.OR.NMOD.EQ.11)) AK2=THETA(3)*ak1
	goto 51
50	continue
c for ip2=+/-2 (cannot have ip1=2 in this case)
	if(nmod.eq.9) AK1=theta(2)/x2	!k-1/k+2
	ak2=ak1
	if(nmod.eq.11) AK1=theta(2)
	if(nmod.eq.11) AK2=theta(3)/x2	!k-2/k+2
51	continue
	AKB=THETA(3)
	IF(NMOD.EQ.10.OR.NMOD.EQ.11) AKB=THETA(4)
	if(nmod.eq.12) cmax=theta(4)	!MAX CURRENT
	if(nmod.eq.12) goto 12
C FOR Nmod=9,11 now calc value for ALPHA from X1=mean open
c time/bst (sec) and value of k(+2) which is in X2 if value
c is fixed, or in theta(kmax) if it is to be estimated (ip2 neg)
c  (1) mods 10,11 IP1=1 if param are K1,K2; =2 if param are K1, K2/K1.
c  (2) mods 9,11  IP2=1 if param is K2 or (K2/K1 if IP1=2)
c		  IP2=2 if param is k(-2) (must have ip2=1 if ip1=2)
c  (3) mods 9,11  IP2 set negative if k(+2) is a param  to be estimated
	X4=X2*AK2			!k(+2)*K2=k(-2) 1/sec
	X5=X2*AK1			!k(+2)*K1=k(-1)
	X6=2.*X4			!DISSOC RATE=2k(-2)
c	IF(NMOD.EQ.11.AND.CASTAR.GT.0.) X6=X4+X5
	IF(NMOD.EQ.11.AND.iequiv.eq.1) X6=X4+X5
	IF(NMOD.EQ.9.OR.NMOD.EQ.11)
     & BA=BA*(X1-(BA/X6)/(BA+X6))/(1.0 + BA/X6)
C NOTE: (X1-(BA/X6)/(BA+X6))=TOTAL OPEN TIME/BURST SO
C NOW BA=BETA/ALPHA FOR ALL MODELS
12	x3=ak1*ak2
c ? fix problem with zero div even when const to be pos (12-nov-85)
	if(x3.gt.1.e-25) goto 120
	ak1=1.e-10
	if(nmod.eq.8.or.nmod.eq.9) ak2=ak1
	goto 12
120	T=(XV*XV)/x3
c	T=(XV*XV)/(AK1*AK2)
	X3=2.*XV/AK1
	IF((NMOD.EQ.10.OR.NMOD.EQ.11).AND.iequiv.eq.1) X3=XV/AK1+XV/AK2
	ycalcv=1.+ X3 +T*(1.+(1.+XV/AKB)*BA)
	ycalcv=BA*T/ycalcv
	if(nmod.eq.12) ycalcv=CMAX*ycalcv
	RETURN
c ind model (13)
52	ba=theta(1)
	ak1=theta(2)
	ba2=ba
	ak2=ak1
	km=2
	if(ip1.eq.1) goto 325
	km=km+1
	ak2=theta(km)	!K2
325	if(ip2.eq.1) goto 326
	km=km+1
	ba2=theta(km)	!beta2/alpha2
326	km=km+1
	akb=theta(km)
	ca1=xv/ak1	!cA1
	ca2=xv/ak2	!cA2
	x3=xv/akb	!cB
	X4=ba*x1	!cA1*b/a
	X5=ba2*x2	!cA2*b2/a2
	x6=x4*x5
	T=1.+ca1+ca2+x4+x5+ ca1*ca2*(1.+ba+ba2) + x6 + x6*x3
	ycalcv=x6/T
	return
C MOD 14 & 15
424	ak11=theta(1)
	ak21=ak11/theta(2)
	ak22=ak21/theta(3)
	ba=theta(4)
	akb=theta(5)
	if(nmod.eq.14) goto 425
c now nmod=15- calc alpha and hence beta/alpha from mean burst length (x1)
c If k(+2) is fixed (kmax=5) it is in x2; if to be estimated kmax=6
	beta=theta(4)
	if(kmax1.eq.5) af2=x2		!k(+2)
	if(kmax1.eq.6) af2=theta(6)	!k(+2)
	ak12=theta(2)*ak22		!K12
	ab1=af2*ak12	!k(-2) for subunit print 1
	ab2=af2*ak22	!k(-2) for subunit print 2
	x3=1./(beta+ab1+ab2)	!mean short gap length
	x4=beta/(ab1+ab2)	!mean short gaps/burst
	alpha=(1.+x4)/(x1-(x3*x4))
	ba=beta/alpha
c
425	continue	!calc Popen from equilib consts
	ca1=xv/ak11
	ca2=xv/ak21
	x3=xv/ak22
	x4=ca1*x3*ba	!c11*c22*b/a
	x5=xv/akb
	T=1.+ca1+ca2+ca1*x3+x4*(1.+x5)
	ycalcv=x4/T
C
304	CONTINUE	!=NEXT MODEL
	RETURN
c
c Model 16
329	continue
	Erev=theta(1)
	g0=theta(2)
	H=theta(3)
	if(nmod.eq.17) goto 330
	if(nmod.eq.18) goto 331
	ycalcv= (xv-Erev)*g0*exp(xv/H)
	RETURN
330	ca=theta(4)		!normalised agonist conc x/K1
c Mod 17: saturation with beta/alpha V-dep. Assume K1 not V-dep
	ba=theta(5)		!beta/alpha at 0 mV
	baV=ba*exp(xv/H)	!beta/alpha at xv mV
	popen=(bav*ca*ca)/(1. + 2.*ca + ca*ca*(1.+bav))
	ycalcv= (xv-Erev)*g0*popen
	RETURN
c Mod 18: saturation with K1 V-dep. Assume beta/alpha not V-dep
331	ca=theta(4)		!normalised agonist conc at E=0
	ba=theta(5)		!beta/alpha
	caV=ca*exp(xv/H)	!normalised at xv mV
	popen=(ba*caV*caV)/(1. + 2.*caV + caV*caV*(1.+ba))
	ycalcv= (xv-Erev)*g0*popen
	RETURN
c Model (21) I/V for GHK with Mg block
21	continue
	gNa=theta(1)		!conductance in pS
	gCs=theta(2)		!conductance in pS
	Eoff=theta(3)	!offset in Erev, relative to GHK value
	aK1=theta(4)		!in mM
	aK2=theta(5)		!in mM
	H1=theta(6)
	H2=theta(7)
	v1=25.249		!mV at 20 deg
	F=96484.56		!Coulomb/mole
	cNa=0.15		!molar
	cCs=0.14		!molar
	PNa=(gNa*v1)/(F*cNa)    !pS*mV=femtoamps
	PCs=(gCs*v1)/(F*cCs)    !pS*mV=femtoamps
	E=xv-Eoff			!xv,Eoff in mV
c problem- if E=0 then e1/(1-exp(-e1)) goes to 0/0; limiting value
c is actually 1.0. If e1<0.0001 then use 1st order approx!
	e1=E/v1
	ex=exp1(-e1)		!use exp1.for to prevent under/overflow
	if(abs(e1).gt.0.0001) then
	   e2=e1/(1.0-ex)
	else
	   e2=1.0/(1.0-0.5*e1)
	endif
c	cur=(0.001*E*F/v1)*(1./(1.-ex))*(PCs*cCs - PNa*cNa*ex)
	cur=(0.001*F*e2)*(PCs*cCs - PNa*cNa*ex)
c For Mg block from inside: k(+1)=---> , k(-1)= <---- so K1= k(-1)/k(+1)
c for dissoc outwards k(-2) = --->.   Define K2=k(-2)/k(+1) (all at 0 mV)
c K1,K2 have V dep given by H1,H2 in mV
	cMg=setx(j)		!Mg conc for set j in mM
	e1=E/H1
	if(e1.gt.87.) e1=87.	!exp1 prevents underflow but not overflow
	e2=E/H2
	if(e2.gt.87.) e2=87.	!exp1 prevents underflow but not overflow
	aKeff=aK1*exp1(e1) + aK2*exp1(e2)
	pMg=cMg/(cMg+aKeff)		!Mg occupancy
	YCALCV=cur*(1.0-pMg)
	RETURN
c Model (22) Exponential voltage dependence of rate constant
22	continue
	tau0=theta(1)
	H=theta(2)
	ycalcv=tau0*exp(xv/H)
c
c Now exit
999	continue
c
	RETURN
	END


