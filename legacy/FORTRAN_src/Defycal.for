	subroutine DEFYcal(xmin,xmax,ncalc,logx,theta,Xcal,Ycal,jval,xgt0,
     & setx,njset,ndimc,ndc1)
c Modif 04/12/92 09:41pm so j=jval taken from arguments, rather than using
c j=jset from common
c To define Ycal,Xcal etc for SET #jset (in COMMON)
c This is for plotting in CVFIT (so when several sets are
c being fitted, possibly with different equations, or different number
c of components, the calc curves can be found while the current model is in
c effect).  The values of ncurvc,icurvc(),ncal(),iline() are now calc in
c CVIN2, after ifitmode is defined.
	real*4 theta(100)
	real*4 xcal(ndc1,ndimc),ycal(ndc1,ndimc)
	real*4 setx(njset)
	logical logx,constr,fline,xgt0
	logical debug,caplock,deb
	COMMON/BLOCK1/constr,nset,nfit,nsfit,
     & Xv,alfa,kmax,ncomp,nmod,fline,nomit,jomit(20),jset,ifitmode
	common/potrat/jfirst,iset,kmax1			!for pot ratios etc
c
	debug()=caplock()
c
	dxcalc=(Xmax-Xmin)/float(ncalc-1)
c
c	j=jset
	j=jval
c
	xmin1=xmin		!lowest point on calc curve
c For general exponentials (nmod=19) prevent calc curve starting below x=0
	if(nmod.eq.19.and.xgt0) then
	   if(xmin1.lt.0.) xmin1=0.
	endif
c
	ival=1	!so value of xv in call to ycalcv is used
	deb=debug()
	do 307 i=1,ncalc
	   xcal(i,j)=Xmin1+float(i-1)*dxcalc
	   if(logx) xcal(i,j)=10.0**xcal(i,j)		!non-log values
	   Xv=xcal(i,j)	!VALUE OF X
c=	   Ycal(i,j)=YCALCV(kmax,theta,L,j)
	   Ycal(i,j)=YCALCV(kmax,theta,xv,ival,j,setx,njset)
c
	   if(deb) then
		print 20,i,j,xcal(i,1),ycal(i,1),xcal(i,1),ycal(i,2)
20	 	format(' i,j,xc(i,j),yc(i,j)= ',2i8,4g13.6)
	   endif
307	continue	!end of i loop
c	if(debug()) print 56,nmod,ncalc,kmax,L,dxcalc,
c     & jsav,jset,ilog,iscal,inumx,inumy,ilabel,
c     & (theta(i),i=1,kmax)
c	if(debug()) pause
c
	RETURN
	end

