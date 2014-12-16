	function SSDCV(K,THETA,Xobs,yobs,w,nj,juse,setx,niobs,njset)
C SUM OF SQUARES FOR FITTING MODEL #NMOD WITH CVFIT
c With iver=1003 have not got allocatable arrays in common, so
c need to add to parameters Xobs,yobs,w,nj,juse and (for ycalcv) setx.
c Modif 01/17/90 09:59am to allow for new options concerning which
c sets are to be fitted.
c If ifitmode=1 then only one set fitted (#jset, in common)
c If ifitmode=2 then all sets fitted separately in loop, the current
c    set again being specified by #jset, in common
c If ifitmode=3-5 then all sets fitted simultaneously with single eqn
C
	real THETA(100)
	LOGICAL CONSTR,FLINE,fixratio
c=	common/data/Xobs(100,10),Yobs(100,10)
c=	COMMON/SD/w(100,10)
	integer juse(njset),nj(njset)
	real*4 xobs(niobs,njset),yobs(niobs,njset),w(niobs,njset)
	real*4 setx(njset)
	COMMON/BLOCK1/constr,nset,nfit,nsfit,
     & Xv,alfa,kmax,ncomp,nmod,fline,nomit,jomit(20),jset,ifitmode
	COMMON/BLOCK2/castar,X1,X2,iequiv,ip1,ip2
	logical logyfit,norm
	COMMON/BLOCK3/logyfit,norm,xnorm
	common/fix/fixratio,ir1,ir2,rval	!to fix ratio of 2 parameters
	common/potrat/jfirst,iset,kmax1			!for pot ratios
c (kmax=actual no of param to be estimated; kmax1=number before extra param
c  added for ifitmode=4,5)
c
	S=0.0
	IF (.NOT.CONSTR) GOTO 1
	DO 2 I=1,K
2	IF(THETA(I).LT.1.0E-20) THETA(I)=1.0E-20
1	CONTINUE
	NFIT=0
C
	if(fixratio) then
	   theta(ir1)=rval*theta(ir2)
	endif
c
	ival=1	!so value of x in call to Ycalcv is used
c
	if(logyfit) goto 30
	if(ifitmode.ge.3.and.ifitmode.le.5) goto 20  !fit all sets at once
c
c For ifitmode=1,2 fit only set #jset (no omissions at present)
	iset=1	!needed for nmod=26,27
	do i=1,nj(jset)
c	   i1=i
c	   Yc=Ycalcv(kmax,THETA,i1,jset)
	   xv1=Xobs(i,jset)
	   Yc=YCALCV(kmax,theta,xv1,ival,jset,setx,njset)
	   dev=Yobs(i,jset)-Yc
	   S=S + W(i,jset)*dev*dev
	   nfit=nfit+1
	enddo
	goto 99
c
20	continue	!fit all sets at once for ifitmode=3
	iset=0
c	do 3 j=1,nset
c	j1=j
	do m=1,nsfit
	   j=juse(m)		!set # used
c=	   j1=j
	   iset=iset+1
c iset is in COMMON -nec for YCALCV when ifitmode=4,5 if some sets missed so
c parameter # in theta() can be matched with proper data set
	   do 4 i=1,nj(j)
c=		i1=i
		if(nomit.gt.0) then
c omit specified obs. Is I equal to any of the elements of JOMIT?
		   do L=1,nomit
			if(i.eq.jomit(L)) goto 4
		   enddo
		endif
c=		YC=YCALCV(kmax,theta,i1,j1)
		xv1=Xobs(i,j)
		Yc=YCALCV(kmax,theta,xv1,ival,j,setx,njset)
C YCALCV calculates i'th value in j'th data set
		dev=Yobs(i,j)-Yc     !normal calc
		S=S + w(i,j)*dev*dev
		nfit=nfit+1
4	   continue    !end of i=1,nj(j) loop
	enddo    !end of j=1,nset loop
	goto 99
c
c Separate section for logyfit fit
30	continue
	do m=1,nsfit
	   j=juse(m)		!set # used
c In this (logyfit) case the equations are linear in the ymax(j), so they can be
c calculated here for the current K and nH and transmitted to YCALCV in
c theta(3) to theta(nsfit+2)
	   zmax=0.0
	   do i=1,nj(j)
		aK=theta(1)
		c=xobs(i,j)/aK
		if(nmod.eq.27) then
		   hill=theta(2)
		   c=c**hill
		endif
		if(iequiv.eq.0) then
		   x=c/(c+1.0)	!increasing (normal) curve
		else
		   x=1.0/(1.0 + c)	!decreasing curve
		endif
		zmax=zmax + alog(yobs(i,j)) - alog(x)
	   enddo
	   zmax=zmax/nj(j)	!log(Ymax(j))
	   theta(m+2)=exp(zmax)	!Ymax(j)
	enddo
c All log(ymax) now found, so now get SSD
	iset=0
	do m=1,nsfit
	   j=juse(m)		!set # used
	   iset=iset+1
c iset is in COMMON -nec for YCALCV when ifitmode=4,5 if some sets missed so
c parameter # in theta() can be matched with proper data set.
	   do 5 i=1,nj(j)
		if(nomit.gt.0) then
c omit specified obs. Is I equal to any of the elements of JOMIT?
		   do L=1,nomit
			if(i.eq.jomit(L)) goto 5
		   enddo
		endif
		xv1=Xobs(i,j)
		Yc=YCALCV(kmax,theta,xv1,ival,j,setx,njset)
C YCALCV calculates i'th value in j'th data set
		dev=alog(Yobs(i,j)) - alog(Yc)
		S=S + dev*dev			!no weights!
c		dev=Yobs(i,j)-Yc     !normal calc
c		S=S + w(i,j)*dev*dev
		nfit=nfit+1
5	   continue    !end of i=1,nj(j) loop
	enddo    !end of j=1,nset loop
c
99	SSDCV=S
	RETURN
	END

