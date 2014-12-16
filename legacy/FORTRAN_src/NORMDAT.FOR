	subroutine NORMDAT(j1,xval,yval,weight,theta,mode,
     & ncurvd,icurvd,ndat,icurvw,logyfit,juse,nsfit,ymax1,
     & ndv1,ndimd,niobs,njset,ndth,pon)
c To normalise a set of curves relative to the specified maximum for each set
c (in CVDISP) and put individual normalised points in xcal(i,j1), means in
c ycal(j1+1)
c mode=0 does the normalisation first
c mode=1 assumes that individual normalised obs are already in column j1
c	and just calcs mean and SD in col 2
	integer juse(njset)
	real*4 xval(ndv1,ndimd),yval(ndv1,ndimd)
	real*4 weight(niobs,njset)
	real*4 theta(ndth)
	logical discprt,equal,newx,pon
	ALLOCATABLE::nrep,ssy
	real*4 ssy(:)
	integer nrep(:)
c
	dimension ndat(ndimd),icurvd(ndimd)
	integer icurvw(ndimd)	!new for VPLOT5
c for calc curves
	logical logyfit
	real*4 ymaxj(100)			!in common/max
	common/max/ymaxj
	common/dp/discprt
c
c
	ncurvd=2		!one data set only + 1 for means
	icurvd(1)=j1
	icurvd(2)=j1+1
	icurvw(1)=0		!no SD initially
	icurvw(2)=1		!plot SD with means
c
	if(mode.eq.1) goto 1
c
	i1=0
	iset=0
	ymax1=-1.e37	!largest point to be plotted
	if(logyfit) then
	   do m=1,nsfit
		j=juse(m)
		iset=iset+1
		ymaxj(j)=theta(j+2)
		do i=1,ndat(j)
		   i1=i1+1
		   xval(i1,j1)=xval(i,j)
		   yval(i1,j1)=yval(i,j)/ymaxj(j)	!normalised y
		   if(yval(i1,j1).gt.ymax1) ymax1=yval(i1,j1)
		   weight(i1,j1)=weight(i,j) !points skipped if w=0,isdev=0
		enddo
	   enddo
	else
	   Y0=theta(1)
	   do m=1,nsfit
		j=juse(m)
		iset=iset+1
		do i=1,ndat(j)
		   i1=i1+1
		   xval(i1,j1)=xval(i,j)
		   yval(i1,j1)=y0 + (yval(i,j)-Y0)/ymaxj(j)	!normalised y
		   if(yval(i1,j1).gt.ymax1) ymax1=yval(i1,j1)
		   weight(i1,j1)=weight(i,j) !points skipped if w=0,isdev=0
		enddo
	   enddo
	endif
	ndat(j1)=i1
	if(ndat(j1).gt.ndv1) then
	   call BELL(2)
	   print 70,ndat(j1)
70	   format(
     &' To use this option the ''maximum number of obs'''
     &' must be set to at least ',i5)
	   STOP
	endif
c
1	continue
c Sort x values into ascending order
	call SORTr3(Xval,Yval,Weight,j1,ndat(j1),.true.,
     &   ndv1,ndimd,niobs,njset)
c Calc second data set for means (with SD)
	j2=j1+1	!'column' for means
c	First find how many x values are the 'same' (within 1% say)
	ALLOCATE(nrep(ndat(j1)),ssy(ndat(j1)))
	m=0	!counts number of different x values=ndat(j2) -need
c			!the mean and SD of y at each x
c NB omit values when weight=0
	do i=1,ndat(j1)		!for each original x value:
	   x=xval(i,j1)
	   y=yval(i,j1)
	   w=weight(i,j1)
	   newx=.true.
	   if(m.gt.0) then
		do n=1,m		!check if this x already found
		   if(equal(x,xval(n,j2)).and.w.ne.0.0) then	!yes
			yval(n,j2)=yval(n,j2)+y	!accum total
			ssy(n)=ssy(n)+y*y		!accum ssd
			nrep(n)=nrep(n)+1		!1 value at this x
			newx=.false.	!current x value already found
		   endif
		enddo
	  endif
	  if(newx.and.w.ne.0.0) then
	     m=m+1
	     xval(m,j2)=x
	     yval(m,j2)=y		!first y
	     ssy(m)=y*y		!first ssd
	     nrep(m)=1				!1 value at this x
	   endif
	enddo
	ndat(j2)=m
c	calc mean and SD
      print 35
     	if(pon) write(7,35)
      if(discprt) write(8,35)
35	format(/,
     & '    Values for normalised curve:',/,
     & '    x value          n      mean y            SD(mean)')
	do i=1,ndat(j2)
	   en=float(nrep(i))
	   s=ssy(i) - (yval(i,j2)*yval(i,j2)/en)
	   yval(i,j2)=yval(i,j2)/en		!yval=mean now
	   if(nrep(i).ge.2.and.s.gt.1.e-10) then
		weight(i,j2)=(en*(en-1.))/s
	   else
		weight(i,j2)=-1.	!so point plotted with no SD
	   endif
	   if(weight(i,j2).gt.1.e-10) then
		sd=1.0/sqrt(weight(i,j2))
	   else
		sd=-1.
	   endif
c print values
         print 34, xval(i,j2),nrep(i),yval(i,j2),sd
     	   if(pon) write(7,34) xval(i,j2),nrep(i),yval(i,j2),sd
         if(discprt) write(8,34) xval(i,j2),nrep(i),yval(i,j2),sd
34	   format(2x,g13.6,3x,i4,3x,g13.6,3x,g13.6)
	enddo
	DEALLOCATE(nrep,ssy)
c
	RETURN
	end

	logical function EQUAL(x,x1)
c equal (within 1%)
	equal=.false.
	p=0.01
	if(abs(x).lt.1.e-10) then
	   equal=abs(x1).lt.p
	else
	   equal=abs((x-x1)/x).lt.0.01
	endif
	RETURN
	end

