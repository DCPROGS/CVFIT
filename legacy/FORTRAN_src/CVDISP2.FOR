	subroutine CVDISP2(Xobs,Yobs,ISYM,ILOG,ISCAL,IJOIN,iline,
     & XMIN,XMAX,YMIN,YMAX,XTIC,YTIC,ixlo,ixhi,iylo,iyhi,y0,yinf,
     & xcross,ycross,titlex,titley,title1,itit,ilabel,inumx,inumy,
     &  xval,yval,ncurvd,ndat,icurvd,ncurvc,ncal,icurvc,Xcal,Ycal,
     &  symsiz,idiskq,ntx,nty,itx,ity,theta,thetsav,fitted,doframe,
     &  dcurve,iw,idest,readini,iflag,ndv1,ndimd,ndc1,ndimc)
c
c Plotting in CVFIT2
c 04/13/92 02:55pm Changed to us VPLOT2 (so weights now in common, rather
c than being argument)
c CVFIT2 07/30/91 04:46pm made for CVFIT2 with xcal,ycal as arguments so
c they are pre-calc during loop for each set
c 12/11/89 10:41am Lahey version for display in CVFIT
c Contrary (I think) to rule book, cannot use Xobs,Yobs directly in call
c to VPLOT1 (OK for set 1 only!), so must copy them to Xval,Yval
c	dimension xcal(2048,20),ycal(2048,20)
c	dimension XVAL(2048,10),YVAL(2048,10)	!can use Xobs directly?
	real*4 xval(ndv1,ndimd),yval(ndv1,ndimd)
	real*4 xcal(ndc1,ndimc),ycal(ndc1,ndimc)
	dimension YOBS(100,10),XOBS(100,10),w(100,10)
	real*4 theta(20),thetsav(20,11)
c	real*4 theta1(20)
c for data
c	dimension ndat(ndimd),icurvd(ndimd),isym(ndimd),ijoin(ndimd)
c for calc curves
c	dimension ncal(ndimc),icurvc(ndimc),iline(ndimc)
	dimension ndat(10),icurvd(10),isym(10),ijoin(10)	!for data
	real symsiz(10)				!for data
	dimension ncal(20),icurvc(20),icurvc1(20),iline(20)	!for calc curve
c	integer*2 videotyp
	LOGICAL pon,debug,fitted,dcurve,doframe,draft,autplt,interp
	character*40 titlex,titley
	character*44 title1
	logical constr,fline		!for COMMON
	logical logx,logy,logity,sqrty,slock,caplock,landscap,readini
	character*1 ans,UC,ans1
	logical discprt
	character defolt*30,cans*30		!to hold text & result of QDIALOG
	common/dp/discprt
	COMMON/BLOCK1/constr,nj(10),nset,nfit,jmiss(10),juse(10),nsfit,
     & Xv,alfa,kmax,ncomp,nmod,fline,nomit,jomit(20),jset,ifitmode
	COMMON/BLOCK2/ CASTAR,SETX(10),X1,X2,iequiv,ip1,ip2
	COMMON/SD/w		  !set isdev=-1 in call if not needed
c
c define function
	pon()=slock()
	debug()=caplock()
c
101	format(a1)
104	FORMAT(I8)
	ncurvc1=ncurvc		!= so input value not altered
	do 1 i=1,20
1	icurvc1(i)=icurvc(i)	!so icurvc in call not altered
c NB weight() always defined in this prog, so set isdev=0 for cases where
c SD bars not defined (e.g. iw=1). If set isdev=-1 then bit in VPLOT2 to
c ignore points with zero weight is ignored
	isdev=0		!no SD available
	if(iw.ge.2.and.iw.le.4) isdev=1	!plot with sd initially
	logx=ilog.eq.1.or.ilog.eq.3.or.ilog.eq.4.or.ilog.eq.6
	logy=ilog.eq.2.or.ilog.eq.3
	logity=ilog.eq.4		!for Hill plot
	sqrty=ilog.eq.5.or.ilog.eq.6
c	ndimd=10		!dimension of Xobs,Yobs	!now a param
c	ndimc=20		!dimension of xcal,ycal
	cbig=2.5
	ifont=4
	ifitype=0		!no display of fitted parameters
	if(nmod.eq.19) then
	   ifitype=4	!exponential fits
	   ncomp=kmax/2
	endif
c	ilabel=-1	!internal labels- now set in CVIN
c	if(debug()) then
c	   ilabel=1
c	   titlex='Conc'
c	   titley='Response'
c	endif
c For Toshiba VIDEOTYP()=18 for graphics mode, =3 for alpha mode
c	if(videotyp().ne.18) call INIPLT(0,.false.,1.0)
c
50	if(fitted.or.dcurve) goto 10
c
	if(idest.eq.0) goto 60		!first time
	print 61
61	format(' Same display? ')
	read 101,ans
	if(UC(ans).eq.'N') goto 60
	iscal=0
	goto 30
c
60	continue
	if(readini) then
	   print 108
108	   format(/)
	   if(ilog.eq.0) then
	      call DCASK('Display Y against X','y',ans)
	   else if(ilog.eq.1) then
	      call DCASK('Display Y against log(X)','y',ans)
	   else if(ilog.eq.2) then
	      call DCASK('Display log(Y) against X','y',ans)
	   else if(ilog.eq.3) then
	      call DCASK('Display log(Y) against log(X)','y',ans)
	   else if(ilog.eq.4) then
	      call DCASK('Display as Hill plot','y',ans)
	   endif
	   if(ans.eq.'N') ilog=0
	else if(nmod.ne.22) then
	   ilog=0
	   ans1='n'
	   if(nmod.ge.8.and.nmod.le.15) ans1='y'
	   if(ifitmode.eq.4.or.ifitmode.eq.5) ans1='y'
	   call DCASK('Display Y against log(X)',ans1,ans)
c	   print 319
c319	   format(' Display Y against log(X) [N] ? ')
c	   read 101,ans
	   if(ans.eq.'Y') ilog=1
	else if(nmod.eq.22) then
	   ilog=0
	   call DCASK('Display log(Y) against X','n',ans)
	   if(ans.eq.'Y') ilog=2
	endif
	logx=ilog.eq.1.or.ilog.eq.3.or.ilog.eq.4.or.ilog.eq.6
	logy=ilog.eq.2.or.ilog.eq.3
	logity=ilog.eq.4		!for Hill plot
	sqrty=ilog.eq.5.or.ilog.eq.6
c	inumx=-1
c	inumy=-1
c	if(ilog.eq.0) goto 40
c	print 41
c41	format('&Use exponent form of numbering for log(x) scale [N]? ')
c	read 101,ans
c	INUMX=-1		X axis in fixed (Fn.d) format
c	INUMX=1		X axis in exponent format (if log)
c	INUMY=-1		Y axis in fixed (Fn.d) format
c	INUMY=1		Y axis in exponent format (if log)
c	if(UC(ans).eq.'Y') inumx=1
c40	continue
c NB axis labels can now be specified in VPLOT
c
c Data sets
c Contrary (I think) to rule book, cannot use Xobs,Yobs directly in call
c to VPLOT1 (OK for set 1 only!), so must copy them to Xval,Yval
	ncurvd=nset
	do 22 j=1,nset
	ndat(j)=nj(j)
	icurvd(j)=j		!display all sets
c	isym(j)=-(j+2)	!filled symbols
	isign=(-1)**j
	isym(j)=isign*(j+2)	!alternate filled/hollow symbols
	symsiz(j)=3.0
	ijoin(j)=-1		!points not joined
	do 22 i=1,nj(j)
	xval(i,j)=xobs(i,j)
22	yval(i,j)=yobs(i,j)
c	ncurvd=1		!1 data set
c	icurvd(1)=1		!- the first
c	ndat(1)=nj(1)
c	isym(1)=3
c	ijoin(1)=-1		!points not joined
	jsav=-1			!display all sets with fitted curve
c
	if(nset.eq.1) goto 27
c
24	continue
	print 23
23	format(' Display:',/,
     & ' (1) All data sets',/,
     & ' (2) Selected data sets',/,
     & ' (3) The fitted data set(s) only',/,
     & '  Option number [3] = ')
	read 104,i
	iopt=3
	if(i.ge.1.and.i.le.3) iopt=i
c
	if(iopt.eq.1) goto 27	!display all already set
	if(iopt.eq.2) then
	   print 14
14	   format('&Number of sets to be displayed= ')
	   read 104,ncurvd
	   do 271 i=1,ncurvd
	   print 272,i
272	   format('&(',i3,'): data set no = ')
	   read 104,icurvd(i)
271	   continue
	else if(iopt.eq.3) then
	   if(ifitmode.eq.1) then
		ncurvd=1
		icurvd(1)=jset		!defined in CVIN
	   else
		ncurvd=0
		do 120 j=1,nset
		if(jmiss(j).eq.1) goto 120	!set omitted
		ncurvd=ncurvd+1
		icurvd(ncurvd)=j
c		if(debug()) print 1201,ncurvd,(jmiss(i),i=1,4),
c     &	(icurvd(i),i=1,4)
c1201		format(' ncurvd=',i8,/,' jmiss=',4i4,/,' icurvd=',4i4)
120		continue
	   endif
	endif
c
27	continue
c
c calc curves:
c NB use ncurvc1 in call to Vplot so ncurvc, now calc in CVIN, not altered
	ncurvc1=0		!none here
	iline(1)=0		!total curve continuous
	iscal=1		!scale internally
	ixlo=-1		!whole screen
	ntx=5
	nty=5
	itx=1		!ticks orientated normally
	ity=1
c	symsiz(1)=-1.0 	!default symbol sizes
	doframe=.true.
	landscap=.true.
	if(nmod.eq.12.or.nmod.eq.16) landscap=.false.
	if(nmod.eq.21) then
	   landscap=.false.
	   iscal=0		!use input values
	   ilog=0
	   logx=.false.
	   logy=.false.
	   call MINMAX2(xval,yval,ndat,xmin1,xmax1,ymin1,
     &   ymax1,logx,logy,ndimd,ncurvd,icurvd)
	   call FIXAX(xmin1,xmax1,xmin,xmax,xtic,0)	!always non-log
	   call FIXAX(ymin1,ymax1,ymin,ymax,ytic,0)
	   xcross=0.
	   ycross=0.
	   itx=0		!ticks central
	   ity=0
	   ntx=5
	   nty=5
	   doframe=.false.
	endif
c	draft=.true.
	autplt=.false.
c
c	if(debug()) print 56,nmod,ncalc,kmax,L,dxcalc,
c     & jsav,jset,ilog,iscal,inumx,inumy,ilabel,
c     & (theta(i),i=1,kmax)
c	if(debug()) pause '0'
c
c
c	if(.not.debug()) goto 30
c	print 31,ncurvd,ndimd,ndimc,(icurvd(i),i=1,8),(ndat(i),i=1,8),
c     & (isym(i),i=1,8)
c31	format(' ncurvd,ndimd,ndimc= ',3i4,/,
c     & ' icurvd= ',8i5,/,
c     & ' ndat  = ',8i5,/,
c     & ' isym  = ',8i5,/)
c	pause '00'
c
30	continue	!return here to plot
c
	iask=2	!do not ask before leaving display; leave graph on screen
c=	iask=-2
	ncjump=0
	nvjump=0
	autplt=.false.
	isetcol=0
	itrace=0
	call VPLOT4(XVAL,YVAL,NDAT,icurvd,ncurvd,ijoin,symsiz,ndimd,
     & XCAL,YCAL,NCAL,icurvc1,ncurvc1,iline,ndimc,ISYM,ILOG,ISCAL,
     & XMIN,XMAX,YMIN,YMAX,XTIC,YTIC,xcross,ycross,ntx,nty,itx,ity,
     & iXLO,iXHI,iYLO,iYHI,y0,yinf,inumx,inumy,ncjump,nvjump,ivplot,
     & titlex,titley,ilabel,doframe,idiskq,autplt,draft,itit,title1,
     & cbig,ifont,landscap,fitted,iask,theta,ifitype,ncomp,isdev,interp,
     & isetcol,itrace,ndv1,ndc1)
c
c Now still in graphics!
	if(pon().and.ilog.eq.4) write(7,25)y0,yinf
      if(discprt.and.ilog.eq.4) write(8,25)y0,yinf
25	format(' For Hill plot y(0), y(inf)= ',2g13.6)
	if(fitted) then
 	   call VIDEOMOD(3)	 !utility lib- back to alpha mode
	   idest=1
	   if(iflag.ne.2) then
	      call DCASK('Repeat same display','N',ans)
	      if(ans.eq.'Y') goto 50
	   else
	      print 232
232	      format(
     &	' (0) Carry on',/,
     &	' (1) Display with all calculated curves',/,
     &	' (2) Display simultaneous fit (parallel curves) only',/,
     &	' (3) Display separate fit to each set only',/,
     &	' Option number [0] = ')
		read 104,iopt
		if(iopt.eq.1) then
		   ncurvc1=ncurvc
		   do 40 i=1,ncurvc1
40			icurvc1(i)=icurvc(i)
		   goto 50
		else if(iopt.eq.2) then
		   ncurvc1=ncurvc/2
		   do 41 i=1,ncurvc1
41		    icurvc1(i)=icurvc(i)
		   goto 50
		else if(iopt.eq.3) then
		   ncurvc1=ncurvc/2
		   do 42 i=1,ncurvc1
42		    icurvc1(i)=icurvc(i+nset)
		   goto 50
		endif
	   endif
	else
c Now still in graphics mode, so write to the dialog box defined in VHIST
	   idest=0
	   ans='Y'
	   call DEFOLTa(ans,defolt)
	   call QDIALOG(1,'Fit this display',defolt,11,cans)
	   call GETINPa(cans,ans)
	   if(ans.eq.'Y') then
		iopt=1
	   else
	      ans='N'
	      call DEFOLTa(ans,defolt)
	      call QDIALOG(1,'Repeat same display',defolt,11,cans)
	      call GETINPa(cans,ans)
	      if(ans.eq.'Y') then
		   iopt=2
		else
		   iopt=3		!end now
		endif
	   endif
c=	   print 231
c=231	   format(
c=     & ' (1) Fit this display',/,
c=     & ' (2) Repeat same display',/,
c=     & ' (3) End now ',/,
c=     & ' Option number [1] = ')
c=	   read 104,iopt
c=	   if(iopt.eq.0) iopt=1
	   if(iopt.eq.2) goto 50
	   if(iopt.eq.3) idest=99
	endif
	RETURN
c
c Now bit done after fit (alpha mode)
10	continue
C
C DISPLAY CALC CURVE AFTER CALC OF ERRORS
c NB use ncurvc1 in call to Vplot so ncurvc, now calc in CVIN, not altered
	iscal=0		!so prev xmin etc used
	itit=1		!so prev title used
	nl=NBLANK1(title1)
	if(nl.eq.1.and.title1(1:1).eq.char(0)) itit=0	!title is blank
c	draft=.true.
c
	if(nmod.eq.-1) then	!problem if diff sets have diff eqns!
	   y0=theta(3)
	   yinf=theta(2)
	endif
c
	if(nmod.eq.21) then
	   print 211
211	   format(' Show:',/,
     &' (0) no extra curves',/,
     &' (1) 1mM + bst length corr',/,
     &' (2) ditto with freq/13 at E>0',/,
     &' (3) both of these',/,' option number= ')
	   read 104,iopt
	   if(iopt.gt.0) then
c calc the curves
		ncalc=201
		do 214 j=5,6
		iline(j)=0
		ncal(j)=ncalc
		do 214 i=1,ncalc
214		xcal(i,j)=xcal(i,4)
c
		do 215 i=1,ncalc
		E=xcal(i,5)
		ycal(i,5)=ycal(i,4)*exp(-E/156.)
		if(E.gt.0.) then
		   ycal(i,6)=ycal(i,5)/13.
		else
		   ycal(i,6)=ycal(i,5)
		endif
c		if(debug()) print 20,i,5,xcal(i,5),ycal(i,5)
c		if(debug()) print 20,i,6,xcal(i,6),ycal(i,6)
215		continue
		print 212
212		format(' Omit data curves ? ')
		read 101,ans
		if(UC(ans).eq.'Y') then
		   ncurvd=0
		   ncurvc1=2
		   if(iopt.eq.1) ncurvc1=1
		   icurvc1(1)=5
		   icurvc1(2)=6
		else
		   ncurvc1=4
		   if(iopt.eq.1) ncurvc1=3
		   icurvc1(1)=3
		   icurvc1(2)=4
		   icurvc1(3)=5
		   icurvc1(4)=6
		endif
         endif
	endif
	if(debug()) then
	   print 217,ncal(1),iline(1),ncurvc1,(icurvc1(i),i=1,ncurvc1)
217	   format(' ncal(1),iline(1),ncurvc1 =',3i4,'  icurvc1= ',10i3)
	   pause
	endif
c
	goto 30		!plot it
c
	end



