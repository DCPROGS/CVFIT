	subroutine DRAWCONT(aK1,aK2,ncont,svals,jset,xobs,yobs,wt,
     &	niobs,nj,njset,nint,aKint,Ymint)
c To draw SSD contours for fit of hyperbola (called from CVSIMPg)
c Two values of Ymax are calculated for each specified K, by call
c to CONTSSD.  In the call, specify the range of K values to be used
c as aK1 to aK2, and the number of contours to be calculated, as ncont
c contours, for SSD=svals(i), i=1,2,..,ncont
c The path followed towards minimum is plotted as data -nint points
c in Akint,Ymint,
	real*4 svals(ncont)
	allocatable:: yhi,ylo,ycal,xcal,x,y,w
	real*4 yhi(:),ylo(:),ycal(:,:),xcal(:,:),x(:),y(:),w(:)
	allocatable:: XVAL,YVAL
	real XVAL(:,:),YVAL(:,:)		!for VPLOT
	real*4 xobs(niobs,njset),yobs(niobs,njset),wt(niobs,njset)
	integer nj(njset)
	logical complex,allocated,first
	real aKint(nint),Ymint(nint)
c For vplot
	integer ndat(3),icurvd(3),isym(3),ijoin(3)
	integer icurvw(3)		!new for VPLOT5
	real symsiz(3)
c for calc curves
	allocatable:: ncal,icurvc,iline
	integer ncal(:),icurvc(:),iline(:)
	logical fitted,doframe,autplt,interp,landscap
	character*40 titlex,titley
	character*64 title1
	real theta(1)	!not used
	logical mono
	COMMON/cols/icol(100),mono
c
c First calculate the contours
	ncalc=200
	if(aK1.gt.aK2) then
	   temp=aK1
	   aK1=aK2
	   aK2=temp
	endif
	dk=(aK2-aK1)/float(ncalc)
c
	n2=2*ncalc		!get two y values for eah x (=K) value
	Allocate(ylo(ncalc),yhi(ncalc))
	allocate(ycal(n2+1,ncont),xcal(n2+1,ncont))
	allocate(ncal(ncont),icurvc(ncont),iline(ncont))
	allocate(xval(nint,3),yval(nint,3))	!not used
	ndc1=n2+1
	ndimc=ncont
	ndv1=nint
	ndimd=3
c
	nobs=nj(jset)
	Allocate(x(nobs),y(nobs),w(nobs))
c Copy data to x,y
	do i=1,nobs
	   x(i)=xobs(i,jset)
	   y(i)=yobs(i,jset)
	   w(i)=wt(i,jset)
	enddo
	ymin1=1.e37
	ymax1=-1.e37
c Calculate contours
	do n=1,ncont	!contour number n
	   SSD=svals(n)
	   n1=0	!number of pairs of points to be plotted for contour #n
	   first=.false.
	   i=0		!index for x values
	   do k=1,ncalc	!loop through K values
		aK=aK1+float(k-1)*dk	!current K value
		call CONTSSD(y1,y2,complex,x,y,w,nobs,SSD,aK)
c CONTSSD calculates the SSD contours for fit of hyperbola
c Returns the two values of Ymax in y1,y2 (if complex=false), given
c data in x,y,n and specified values of K and SSD.
		if(.not.complex) then
		   first=.true.
		   n1=n1+1
		   yhi(n1)=y1
		   ylo(n1)=y2
		   xcal(n1,n)=aK
		else
		   if(first) goto 1	!end of this contour passed
		endif
	   enddo
1	   continue
c One contour now completed, for current SSD.  Keep it in ycal(i,n).
c To display must have correct sequence (if points are to be joined by lines)
	   ncal(n)=n1
	   do i=1,n1
		ycal(i,n)=ylo(i)
	   enddo
	   n2=2*n1
	   do i=1,n1
		j=n2-i+1
		ycal(j,n)=yhi(i)
		xcal(j,n)=xcal(i,n)
	   enddo
	   ncal(n)=n2		!total number to be plotted
	   do i=1,ncal(n)
		if(ycal(i,n).gt.ymax1) ymax1=ycal(i,n)
		if(ycal(i,n).lt.ymin1) ymin1=ycal(i,n)
	   enddo
c    Add one more point, = first point, so contour closed
	   n2=n2+1
	   ncal(n)=n2
	   xcal(n2,n)=xcal(1,n)
	   ycal(n2,n)=ycal(1,n)
	enddo
c
	deallocate(yhi,ylo,x,y,w)
c
c Now plot them
	ncurvd=3
	ndat(1)=1		!init guess
	ndat(2)=nint-2	!the rest
	ndat(3)=1		!the minimum
	xval(1,1)=aKint(1)
	yval(1,1)=Ymint(1)
	do i=1,nint-2
	   xval(i,2)=aKint(i+1)
	   yval(i,2)=Ymint(i+1)
	enddo
	xval(1,3)=aKint(nint)
	yval(1,3)=Ymint(nint)
	isym(1)=-3	!square
	isym(2)=-1	!triangle
	isym(3)=-7	!circle
	symsiz(1)=2.
	symsiz(2)=1.
	symsiz(3)=2.
	do i=1,3
	   icurvd(i)=i
	   ijoin(i)=0
	   icurvw(i)=-1
	enddo
	isetcol=0		!default colours
	isetcol=1		!default if -1 on entry
	do i=1,100
	   icol(i)=-1
	enddo
c  Calc curves
c  ISCAL=4 if input values of xmin,xmax,ymin,ymax to be used; rest internal
c	iscal=1		!scale internally
	iscal=4
	xmin=0.
	xmax=aK2
	ymin=0.
	ymax=ymax1
c
c	call FIXAX(aK1,aK2,xmin,xmax,xtic,0)	!always non-log
c	call FIXAX(ymin1,ymax1,ymin,ymax,ytic,0)
	ncurvc=ncont
	do i=1,ncont
	   icurvc(i)=i
	   iline(i)=0		!continuous
	enddo
c	ixlo=-1		!whole screen
	ixlo=1300		!square graph
	ixhi=5500
	iylo=1500
	iyhi=5700
	ntx=5
	nty=5
	itx=1		!ticks orientated normally
	ity=1
	symsiz(1)=-1.0 	!default symbol sizes
	doframe=.true.
	landscap=.true.
	autplt=.false.
	ncjump=0
	nvjump=0
	autplt=.false.
	fitted=.false.
	itrace=0
	kwi=niobs   !dimensions of weight()
	kwj=njset
	kth=1		!declared dimension of theta
	cbig=2.5
	ifont=4
	ilog=0
	ilabel=1
	titlex='value of K'
	titley='value of Ymax'
c
	call VPLOT5(XVAL,YVAL,NDAT,icurvd,ncurvd,ijoin,symsiz,ndimd,
     & XCAL,YCAL,NCAL,icurvc,ncurvc,iline,ndimc,ISYM,ILOG,ISCAL,
     & XMIN,XMAX,YMIN,YMAX,XTIC,YTIC,xcross,ycross,ntx,nty,itx,ity,
     & iXLO,iXHI,iYLO,iYHI,y0,yinf,inumx,inumy,ncjump,nvjump,ivplot,
     & titlex,titley,ilabel,doframe,idiskq,autplt,plotonly,itit,title1,
     & cbig,ifont,landscap,fitted,iask,theta,ifitype,ncomp,interp,
     & isetcol,itrace,ndv1,ndc1,weight,kwi,kwj,icurvw,kth)
c
	DEALLOCATE(xcal,ycal,xval,yval)
	DEallocate(ncal,icurvc,iline)
c Leave still in graphics mode
	RETURN
	end



c


