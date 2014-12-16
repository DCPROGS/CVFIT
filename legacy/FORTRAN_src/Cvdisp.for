	subroutine CVDISP(Xobs,Yobs,weight,ISYM,ILOG,ISCAL,IJOIN,iline,
     &  XMIN,XMAX,YMIN,YMAX,XTIC,YTIC,xlo,xhi,ylo,yhi,y0,yinf,
     &  xcross,ycross,titlex,titley,title1,itit,ilabel,inumx,inumy,
     &  xval,yval,ncurvd,ndat,icurvd,ncurvc,ncal,icurvc,Xcal,Ycal,
     &  symsiz,idiskq,ntx,nty,itx,ity,theta,fitted,doframe,
     &  dcurve,nodata,iw,idest,readini,iflag,ndv1,ndimd,ndc1,ndimc,
     &  jmiss,nj,niobs,njset,icurvw,juse,ndth,nset1,thetgues,norm1,
     &  plotrue,iver,ibk,
     &  ncurvsav,xcalsav,ycalsav,ncalsav,icrvcsav,ilinesav)
c
c Plotting in CVFIT
c Modif 04/01/95 02:00pm for VPLOT5
c TEMP MODIF 03/30/95 04:25pm -set isdev=-1 if more than 10 sets to
c prevent error with w(100,10) in vplot
c 04/13/92 02:55pm Changed to us VPLOT2 (so weights now in common, rather
c than being argument)
c CVFIT2 07/30/91 04:46pm made for CVFIT2 with xcal,ycal as arguments so
c they are pre-calc during loop for each set
c 12/11/89 10:41am Lahey version for display in CVFIT
c	dimension xcal(2048,20),ycal(2048,20)
c	dimension XVAL(2048,10),YVAL(2048,10)	!can use Xobs directly?
	real*4 xval(ndv1,ndimd),yval(ndv1,ndimd)
	real*4 xcal(ndc1,ndimc),ycal(ndc1,ndimc)
	real*4 xobs(niobs,njset),yobs(niobs,njset),weight(niobs,njset)
	integer juse(njset)
	real*4 thetgues(100,nset1)
c	dimension YOBS(100,10),XOBS(100,10),w(100,10)
	real*4 theta(ndth),theta1(20)
	real*4 ymaxj(100)			!in common/max
	character*20 title 			!for popmenu
	common/popvalm/nxlo,nylo,nyhi,ictm,ibkm,icfm,icupm	!values for popmenu calls
	character*55 strings(8)			!for popmenu
	character*55 helps(8)			!for popmenu
c for data
	dimension ndat(ndimd),icurvd(ndimd),isym(ndimd),ijoin(ndimd)
	integer icurvw(ndimd)	!new for VPLOT5
	real symsiz(ndimd)
c for calc curves
	integer ncal(ndimc),icurvc(ndimc),iline(ndimc)
	allocatable icurvc1
	integer icurvc1(:)
c For last fit
	real*4 xcalsav(ndc1,ndimc),ycalsav(ndc1,ndimc)
	integer ncalsav(ndimc),icrvcsav(ndimc),ilinesav(ndimc)
c
c	dimension ndat(10),icurvd(10),isym(10),ijoin(10)	!for data
c	real symsiz(10)				!for data
c	dimension ncal(20),icurvc(20),icurvc1(20),iline(20)	!for calc curve
c	integer*4 videotyp
	LOGICAL pon,debug,fitted,dcurve,doframe,draft,autplt,interp
	character*40 titlex,titley
	character*64 title1
	logical constr,fline		!for COMMON
	logical logx,logy,logity,slock,caplock,landscap,readini
c	logical sqrty,present
	character*1 ans,UC
	logical discprt,twodisp,nodata,alpha,setwt,super
	character defolt*30,cans*30		!to hold text & result of QDIALOG
	integer jmiss(njset),nj(njset)
	logical cluster,student
	common/dp/discprt
	COMMON/BLOCK1/constr,nset,nfit,nsfit,
     & Xv,alfa,kmax,ncomp,nmod,fline,nomit,jomit(20),jset,ifitmode
	COMMON/BLOCK2/ castar,X1,X2,iequiv,ip1,ip2
	logical logyfit,norm,norm1,plotrue
	COMMON/BLOCK3/logyfit,norm,xnorm,iostat		!CVIN2,YCALCV,CVDISP
c=	COMMON/SD/weight(100,10)	  !set isdev=-1 in call if not needed
	common/max/ymaxj
	COMMON/cols/icol(100),mono
	common/setbig/jrep		!for SETBIGJ in VPLOT5
	common/user/student,cluster,iwindows
c
c define function
	pon()=slock()
	debug()=caplock()
c
101	format(a1)
	ALLOCATE(icurvc1(ndimc))
	do i=1,100
	   icol(i)=-1	!default unless set otherwise below (isetcol=1)
	enddo
	twodisp=.false.
	iopt1=0
	ncurvc1=ncurvc		!= so input value not altered
	do i=1,ndimc
	   icurvc1(i)=icurvc(i)	!so icurvc in call not altered
	enddo
	logx=ilog.eq.1.or.ilog.eq.3.or.ilog.eq.4.or.ilog.eq.6
	logy=ilog.eq.2.or.ilog.eq.3
	logity=ilog.eq.4		!for Hill plot
c	sqrty=ilog.eq.5.or.ilog.eq.6
c	ndimd=10		!dimension of Xobs,Yobs	!now a param
c	ndimc=20		!dimension of xcal,ycal
	cbig=2.5
c=====================azi
	ifont=1
	ifitype=0		!no display of fitted parameters
	if(nmod.eq.19) then
	   ifitype=6	!exponential fits
	   ncomp=kmax/2
	else if(nmod.eq.20) then
	   ifitype=-6	!geometric fits
	   ncomp=kmax/2
	endif
c
c===========================================================
50	if(fitted.or.dcurve.or.plotrue) goto 10
	if(nodata) goto 60
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
	if(.not.readini) ilog=0
	if(.not.nodata) then		!asked before DEFYCAL if nodata
	   print 108
108	   format(/)
	   i=ilog+1
	   print 109,i
109	   format(
     &    ' (1) Display Y against X',/,
     &    ' (2) Display Y against log(X)',/,
     &    ' (3) Display log(Y) against X',/,
     &    ' (4) Display log(Y) against log(X)',/,
     &    ' (5) Display as Hill plot',/,
     &    ' (6) Display as sqrt(y) against X',/,
     &    ' (7) Display as sqrt(y) against log(X)',/,
     &    ' Option number [',i2,'] = ')
	   call INPUTi(i)
	   ilog=i-1
	endif
	logx=ilog.eq.1.or.ilog.eq.3.or.ilog.eq.4.or.ilog.eq.6
	logy=ilog.eq.2.or.ilog.eq.3
	logity=ilog.eq.4		!for Hill plot
c	sqrty=ilog.eq.5.or.ilog.eq.6
c
	if(nodata) goto 10
c
c Data sets
c Copy Xobs,Yobs to Xval,Yval.
c NB weight() always defined in this prog, so set isdev=0 for cases where
c SD bars not defined (e.g. iw=1). If set isdev=-1 then bit in VPLOT to
c ignore points with zero weight is ignored
c Symbol # must be 1-9 (or -1 to -7) -now fixed in SYMBOL()
	isdev=0		!weight(i,j) defined, but not plotted initially
	if(iw.ge.2.and.iw.le.4) isdev=1	!plot with sd initially
	ncurvd=nset
	do j=1,nset
	   ndat(j)=nj(j)
	   icurvd(j)=j		!display all sets
	   icurvw(j)=isdev		!new for VPLOT5!
c	   isym(j)=-(j+2)	!filled symbols
	   isign=(-1)**j
	   isym(j)=isign*(j+2)	!alternate filled/hollow symbols
	   if(ndat(j).le.10) then
	      symsiz(j)=3.0
	   else if(ndat(j).gt.10.and.ndat(j).le.30) then
	      symsiz(j)=2.0
	   else
	      symsiz(j)=1.5
	   endif
	   ijoin(j)=-1		!points not joined
	   do i=1,nj(j)
		xval(i,j)=xobs(i,j)
		yval(i,j)=yobs(i,j)
	   enddo
	enddo
c For normalised curves if pooled (ifitmode=1), nset=3, set 1=all obs,
c set 2=means (former may be unweighted so show SD for latter only initially)
c and set 3=norm point alone
	if(norm.and.ifitmode.eq.1) then			!normalised to common X value
	   ncurvd=2
	   icurvd(1)=2
	   icurvd(2)=3
	   icurvw(1)=0
	   if(iw.eq.1) then
	      icurvw(2)=-1
	   else if(iw.eq.2) then
	      icurvw(2)=1
	   endif
	   icurvw(3)=-1
	   isym(1)=-3		!filled squares for indiv points
	   symsiz(1)=2.0
c	   icol(1)=10		!green for sep points
	   ijoin(1)=-1		!points not joined
	   isym(2)=-7		!filled circles for means
	   symsiz(2)=3.0
c	   icol(2)=12		!red for means
    	   ijoin(2)=-1		!points not joined
	   ndat(3)=1
	   isym(3)=7		!open circles for normalisation point
	   xval(1,3)=xobs(1,3)
	   yval(1,3)=yobs(1,3)
	   symsiz(3)=4.0
	   icol(3)=12		!red for norm point
    	   ijoin(3)=-1		!points not joined
cc	   isetcol=1
	else if(norm.and.ifitmode.ne.1) then	!draw norm point
	   ns1=nsfit+1
	   ncurvd=ns1		!extra set for norm point
	   ndat(ns1)=1 		!one point only
	   icurvd(ns1)=ns1
	   xval(1,ns1)=xobs(1,ns1)
	   yval(1,ns1)=yobs(1,ns1)
	   icurvw(ns1)=-1
	   isym(ns1)=7		!open circles for normalisation point
	   symsiz(ns1)=4.0
	   icol(ns1)=12		!red for norm point
    	   ijoin(ns1)=-1		!points not joined
cc	   isetcol=1
	endif
c Also display mean & SD if normalised to max fitted in sep fits
	if(norm1.and.ifitmode.eq.1) then
	   mode=1		!alread normalised
	   j1=1		!normalised data in yval(i,1) etc
	   call NORMDAT(j1,xval,yval,weight,theta,mode,
     &     ncurvd,icurvd,ndat,icurvw,logyfit,juse,nsfit,ymax1,
     &     ndv1,ndimd,niobs,njset,ndth,pon())
	   isym(1)=-7
	   isym(2)=-3
	   goto 27
	endif
c
	if(nset.eq.1) goto 27

24	continue
	if(norm) goto 121		!ncurvd, icurvd already set
	iopt=3
	print 23,iopt
23	format(' Display:',/,
     & ' (1) All data sets',/,
     & ' (2) Selected data sets',/,
     & ' (3) The fitted data set(s) only',/,
     & '  Option number [',i2,'] = ')
	call INPUTi(iopt)
c
	if(iopt.eq.1) goto 27	!display all already set
	if(iopt.eq.2) then
	   print 14
14	   format('&Number of sets to be displayed = ')
	   call INPUTi(ncurvd)
	   do i=1,ncurvd
	     print 272,i
272	     format('&(',i3,'): data set no = ')
	     call INPUTi(icurvd(i))
	   enddo
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
120		continue
	   endif
	endif
c Reset symbols etc so they refer only to the curves actually plotted
c if nset>9 (NB number of symbols, like number of coloures, is limited!)
121	continue
	if(nset.gt.9) then
	   do j=1,ncurvd
		isign=(-1)**j
		isym(j)=isign*(j+2)	!alternate filled/hollow symbols
		j1=icurvd(j)
		if(ndat(j1).le.10) then
		   symsiz(j)=3.0
		else if(ndat(j1).gt.10.and.ndat(j1).le.30) then
		   symsiz(j)=2.0
		else
		   symsiz(j)=1.5
		endif
	   enddo
	endif
c
27	continue
c
c calc curves:
c NB use ncurvc in call to Vplot, but restore orig values before leaving
c from ncurvc1
	ncurvc=0		!none here
271	continue		!return here if 'nodata', for display
	iline(1)=0		!total curve continuous
	iscal=1		!scale internally
	xlo=-1		!whole screen
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
	autplt=.false.
30	continue	!return here to plot
	autplt=.false.
333   continue
	 if(plotrue) then
	   autplt=.true.
	   do i=1,ndimd
		symsiz(i)=2.5
	   enddo
	   idest=1
	 endif
c
	iask=2	!do not ask before leaving display; leave graph on screen
c=	iask=-2
	if(ilog.ge.5) ilabel=-1		!use internal labels (modif #1)
	ncjump=0
	nvjump=0
	itrace=0
	if(kmax.le.20) then
	   do i=1,kmax
		theta1(i)=theta(i)	!still have theta(20) in VPLOT
	   enddo
	else
	   ifitype=0
	endif
	kwi=niobs   !dimensions of weight()
	kwj=njset
c Want colours same for sim and sep fits -this can be done via new
c setjbig when nsfit>5 so >10 calc curves.  If fewer than this
c setbigj is not invoked so need to set manually
c NB all icol(i) are -1 here, so need to set the actual colours
c as defined in setcol.for, namely
c Example: if sets 2,5 fitted with nset=7, j=icurvc()=2,5,9,12
c and we want 9,12 to be sane colour as 2,5
c Modified 07/08/97 06:55am because in VPLOT setbigj now uses j=1,2,..,ncurv,
c rather than icurv(j)
	if(iflag.eq.2) then
	   icol(1)=9	!plot colours as in \fort90\graph\setcls
	   icol(2)=12
	   icol(3)=10
	   icol(4)=13
	   icol(5)=11
	   icol(6)=0
	   icol(7)=2
	   icol(8)=1
	   icol(9)=6
	   icol(10)=4	!dark red (used for lo/hi bins of histo
	   jrep=nset !so colours repeat for separate fits
	   if(nsfit.gt.10) jrep=10
cc	   isetcol=1
c Modified 07/07/97 10:00pm in line with changes in use of setbigj
c -now just set the first nset line colours equal to the nset data colours
c and with jrep=nset, the next nset line colours will be set correctly in vplot
	   do m=1,ncurvc/2
c		j=icurvc(m) 	!eg if sets 2,5 fitted with nset=7, j=2,5
		icol(m+10)=icol(m)	!sim fits=same col as data
	   enddo
c 	NB must not reset any icol()>20!
	else
	   jrep=10
	endif
c Colours modified again 10/07/99 09:41am so that if all sets plotted, but not
c all are fitted, data and curve colours are the same (done only if nsfit<nset
c and ncurvc=<10, so as not to disturb changes made for case where there
c are more than 10 calc lines)
	if(nsfit.lt.nset.and.ncurvc.le.10) then
	   icol(1)=9	!plot colours as in \fort90\graph\setcls
	   icol(2)=12
	   icol(3)=10
	   icol(4)=13
	   icol(5)=11
	   icol(6)=0
	   icol(7)=2
	   icol(8)=1
	   icol(9)=6
	   icol(10)=4	!dark red (used for lo/hi bins of histo
	   do m=1,ncurvc
		j=icurvc(m) 	!eg if sets 2,5 fitted with nset=7, j=2,5
		icol(m+10)=icol(j)
	   enddo
	   isetcol=1
	endif
c
	super=.false.	!last fit not superimposed
301	continue
	call VPLOT5(XVAL,YVAL,NDAT,icurvd,ncurvd,ijoin,symsiz,ndimd,
     & XCAL,YCAL,NCAL,icurvc,ncurvc,iline,ndimc,ISYM,ILOG,ISCAL,
     & XMIN,XMAX,YMIN,YMAX,XTIC,YTIC,xcross,ycross,ntx,nty,itx,ity,
     & XLO,XHI,YLO,YHI,y0,yinf,inumx,inumy,ncjump,nvjump,ivplot,
     & titlex,titley,ilabel,doframe,idiskq,autplt,draft,itit,title1,
     & cbig,ifont,landscap,fitted,iask,theta1,ifitype,ncomp,interp,
     & isetcol,itrace,ndv1,ndc1,weight,kwi,kwj,icurvw,kmax,iver)
c
	ibk=icol(71)
c
c Superimpose last fit?
c If so, need to decide which columns (j) of xcal,ycal can be overwritten
c Possible problem if we try to normalise after fits superimposed
	if(super) goto 33		!skip save
	if(ncurvsav.ge.1.and.ncurvc.ge.1) then
	   ans='N'
	   call DEFOLTa(ans,defolt)
	   call QDIALOG(1,
     &   'Superimpose the last fit on this graph',defolt,11,cans)
	   call GETINPa(cans,ans)
	   if(ans.eq.'Y') then
		jmax=0
		do j=1,ncurvc
		   if(icurvc(j).gt.jmax) jmax=icurvc(j)
		enddo
		jleft=ndimc-jmax	!#of free cols in xcal, ycal
		if(jleft.lt.ncurvsav) then	!not room for last fits
		   call BELL(1)
		   call WDIALOG(1,' Not room in arrays for last fit',12)
		else
c put last fit into ycal, xcal for display
c NB icurvc not defined properly if have ncurvc=1 and icurvc(1)=3
c (orig code gave ncurvc=2, icurv=[3  0   0  4] instead of icurvc=[3  4]
		   ncurvc=ncurvc+ncurvsav
		   i1=ncurvsav
		   j=jmax
		   do j1=1,ncurvsav
			j=j+1
			i1=i1+1
c			icurvc(j)=j
			icurvc(i1)=j
			j2=icrvcsav(j1)
			ncal(j)=ncalsav(j2)
c			iline(j)=ilinesav(j2)
			iline(j)=2		!make superimposed curves dashed
			do i=1,ncal(j)
			   xcal(i,j)=xcalsav(i,j2)
			   ycal(i,j)=ycalsav(i,j2)
			enddo
		   enddo
		   super=.true.
		   goto 301		!display with last curves superimposed
		endif
	   endif
	endif
c Keep copy of displayed fitted curve
	if(.not.super) then 	!false if last curve has been added
	   if(ncurvc.ge.1.and.(fitted.or.nodata.or.dcurve)) then
		ncurvsav=ncurvc
		do j1=1,ncurvc
		   icrvcsav(j1)=icurvc(j1)
		   j=icurvc(j1)
		   ncalsav(j)=ncal(j)
		   ilinesav(j)=iline(j)
		   do i=1,ncal(j)
			xcalsav(i,j)=xcal(i,j)
			ycalsav(i,j)=ycal(i,j)
		   enddo
		enddo
	   endif
	endif
33	continue	!skip save
c
c Now still in graphics!
	alpha=.false.
c after display of normalised data following ifitmode=2, return to
c cvin2 to fit data again
	if(iflag.eq.4.and.(.not.super)) then
	   ans='Y'
	   call DEFOLTa(ans,defolt)
	   call QDIALOG(1,'Re-fit normalised data',defolt,11,cans)
	   call GETINPa(cans,ans)
	   if(ans.eq.'Y') then
 		call VIDEOMOD(3)	 !utility lib- back to alpha mode
	      idest=2	!return to cvin2 to re-fit
	   else
		idest=99
	   endif
	   RETURN
	endif
	if(pon().and.ilog.eq.4) write(7,25)y0,yinf
      if(discprt.and.ilog.eq.4) write(8,25)y0,yinf
25	format(' For Hill plot y(0), y(inf)= ',2g13.6)
	if(fitted.or.nodata) then
	   if((nmod.eq.26.or.nmod.eq.27).and.ifitmode.eq.2.and.
     &         iflag.eq.3) then	   !separate fits displayed -now normalise?
		if(alpha) then
		   call DCASK(
     &	   'Normalise these curves w.r.t. fitted maxima','Y',ans)
		else
		   ans='Y'
		   call DEFOLTa(ans,defolt)
		   call QDIALOG(1,
     &     'Normalise these curves w.r.t. fitted maxima',defolt,11,cans)
		   call GETINPa(cans,ans)
		endif
		if(ans.eq.'Y') then
c Can overwrite original xval, yval, xobs, yobs with each set, normalised
c to its own maximum.  If absolute weights given (iw=2,3,4) they must be
c scaled too but arbitrary weights (eg all=1) should be left alone.
		   setwt=iw.eq.2.or.iw.eq.3.or.iw.eq.4
c		   j1=nsfit+1
		   j1=1
		   j2=2		!'column' for means
c        if ncomp=1 need to add Ymax1, Ymax2 etc to get max for set j
		   do m=1,nsfit
			j=juse(m)
			ymaxj(j)=0.0	!sum ymax for all components in set j
			k=2
			if(nmod.eq.26) then
			   kdelt=2
			else if(nmod.eq.27) then
			   kdelt=3
			endif
c for ymax k=2,4,6...for Langmuir, k=2,5,8 for Hill, for each component in one set
			do i=1,ncomp
			   ymaxj(j)=ymaxj(j)+thetgues(k,j)
			   k=k+kdelt
			enddo
		   enddo
c            Now normalise data array and vplot arrays
		   ymax1=-1.e37
		   do m=1,nsfit
			j=juse(m)
			ym=ymaxj(j)
			do i=1,nj(j)
			   yval(i,j)=yval(i,j)/ym
			   if(yval(i,j).gt.ymax1) ymax1=yval(i,j)
			   Yobs(i,j)=Yobs(i,j)/ym
			   if(setwt) weight(i,j)=weight(i,j)/(ym*ym)
			enddo
			do i=1,ncal(j)
			   ycal(i,j)=ycal(i,j)/ym
			   if(ycal(i,j).gt.ymax1) ymax1=ycal(i,j)
			enddo
		   enddo
		   ymax=1.1*ymax1		!new maximum for VPLOT for norm data
		   ytic=0.2
c replace ymax params in thetgues with their normalised values!
		   do m=1,nsfit
			j=juse(m)
			ym=ymaxj(j)
			k=2
			if(nmod.eq.26) then
			   kdelt=2
			else if(nmod.eq.27) then
			   kdelt=3
			endif
			do i=1,ncomp
			   thetgues(k,j)=thetgues(k,j)/ym
			   k=k+kdelt
			enddo
		   enddo
c   now display normalised data and then goto cvin2 to ask which mode
c   to fit normalised data in
		   iflag=4	!increment again so no what to do after display
		   goto 301
		endif
	   endif
c
 	   call VIDEOMOD(3)	 !utility lib- back to alpha mode
	   alpha=.true.
	   idest=1
	   if(nodata) idest=99	!skip likelihhod int
	   if(iflag.ne.2) then
	  	if((nmod.eq.26.or.nmod.eq.27).and.ip1.eq.1.and.
     &  	   ifitmode.eq.3) then
		   goto 50
		else
	         call DCASK('Repeat same display','N',ans)
	         if(ans.eq.'Y') goto 50
		endif
	   endif
	   if(iflag.eq.2) then		!after sim and sep fits for ifitmode=4,5
2321	      print 232
232	      format(
     &	' (0) Carry on',/,
     &	' (1) Display with all calculated curves',/,
     &	' (2) Display simultaneous fit (parallel curves) only',/,
     &	' (3) Display separate fit to each set only',/,
     &	' Option number [0] = ')
		iopt=0
	      call INPUTi(iopt)
c NB if some curves omitted (e.g. because there were not enough points for
c a separate fit to be possible) then ncurvc may not be even and ncurvc1/2
c will give wrong answer.
		if(iopt.eq.1) then
		   ncurvc=ncurvc1
		   do 40 i=1,ncurvc
40			icurvc(i)=icurvc1(i)
		   goto 50
		else if(iopt.eq.2) then
c=		   ncurvc=ncurvc1/2
		   ncurvc=nsfit
		   do 41 i=1,ncurvc
41		    icurvc(i)=icurvc1(i)
		   goto 50
		else if(iopt.eq.3) then
c NB if sets omitted, the calc curves are in xcal, ycal(i,j) with j=set
c number eg if nset=4, and fit 1,2,4 then calc curves for simultaneous
c fit are in ycal(i,j=1,2,4), and when sep fits done too they are in
c ycal(i,j=5,6,8) (i.e. nset+ j value for simultaneous fit) in this
c case icurvc1=(1 2 4 5 6 8), and for sep fits only want icurvc=(5 6 8)
c (as long as all sep fits were done -but some may be omitted if there
c are not enough points for a sep fit)
c		   ncurvc=ncurvc1/2
		   ncurvc=ncurvc1-nsfit
		   do i=1,ncurvc
		    if(icurvc1(i+nsfit).le.0) then
			iopt=1
			goto 2321
		    endif
		    icurvc(i)=icurvc1(i+nsfit)
		   enddo
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
	   if(iopt.eq.2) goto 50
	   if(iopt.eq.3) idest=99
	endif
c restore original values
	ncurvc=ncurvc1		!= so input value not altered
	do i=1,ndimc
	   icurvc(i)=icurvc1(i)	!so icurvc in call not altered
	enddo
	DEALLOCATE(icurvc1)
	RETURN
c
c Now bit done after fit (alpha mode)
10	continue
cc	isetcol=1	!so chosen colours can be set
C
C DISPLAY CALC CURVE AFTER CALC OF ERRORS
	iscal=0	   		!so prev xmin etc used
	if(nodata) iscal=1	!no previous display
	itit=1			!so prev title used
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
	   call INPUTi(iopt)
	   if(iopt.gt.0) then
c calc the curves
		ncalc=501
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
		   ncurvc=2
		   if(iopt.eq.1) ncurvc=1
		   icurvc(1)=5
		   icurvc(2)=6
		else
		   ncurvc=4
		   if(iopt.eq.1) ncurvc=3
		   icurvc(1)=3
		   icurvc(2)=4
		   icurvc(3)=5
		   icurvc(4)=6
		endif
         endif
	endif
	if(debug()) then
	   print 217,ncal(1),iline(1),ncurvc,(icurvc(i),i=1,ncurvc)
217	   format(' ncal(1),iline(1),ncurvc1 =',3i4,'  icurvc1= ',10i3)
	   pause
	endif
c Options for fit of Langmuir/Hill with common K and separate Ymax
c Separate curves (option 1) are already in ycal(i,j) for j=1,nsfit
c For normalised curves, nset=2, set 1=all obs, set 2=means -former may
c be unweighted so show SD for latter only initially.
c so make
c	if((nmod.eq.26.or.nmod.eq.27).and.norm.and.
c     &    ifitmode.eq.1) then			!normalised to common X value
	if(norm.and.ifitmode.eq.1) then			!normalised to common X value
	   icurvw(1)=0		!no SD initially
	   if(iw.eq.1) then
	      icurvw(2)=-1
	   else if(iw.eq.2) then
	      icurvw(2)=1		!plot SD with means
	   endif
	   if(alpha) then
	      print 32,iopt1
32	      format(
     &   ' (1) Plot mean and SD of mean only',/,
     &   ' (2) Plot normalised values + their mean and SD',/,
     &   ' (3) No more plots',/,
     &   ' Option number [',i2,'] = ')
	      call INPUTi(iopt1)
	   else
		title='Plot method'
		strings(1)=
     &	  'Plot mean and SD of mean only'
		strings(2)=
     &	  'Plot mean and SD of mean + separate points'
		strings(3)='No more plots'
		nval=3
		jline=iopt1		!default
		nhelp=0
		helps(1)=' '
		call POPMENU(nxlo,nylo,nyhi,strings,nval,ictm,icfm,icupm,
     &	      ibkm,title,helps,nhelp,jline,ans,ival)
		iopt1=jline
	   endif
	   if(iopt1.eq.3) then
 		call VIDEOMOD(3)	 !utility lib- back to alpha mode
		alpha=.true.
		idest=1
c restore original values
		ncurvc=ncurvc1		!= so input value not altered
		do i=1,ndimc
		   icurvc(i)=icurvc1(i)	!so icurvc in call not altered
		enddo
		RETURN
	   else if(iopt1.eq.1) then
		ncurvd=2
		icurvd(1)=2
		icurvd(2)=3
	      ndat(3)=1 		!one point only
	   else if(iopt1.eq.2) then
		ncurvd=3
		icurvd(1)=1
		icurvd(2)=2
		icurvd(3)=3
	      ndat(3)=1 		!one point only
	   endif
	else if(norm.and.ifitmode.ne.1) then	!draw norm point
	   ns1=nsfit+1
	   ncurvd=ns1		!extra set for norm point
	   icurvd(ncurvd)=ns1
	   icurvw(ns1)=-1
	   ndat(ns1)=1 		!one point only
	   xval(1,ns1)=xobs(1,ns1)
	   yval(1,ns1)=yobs(1,ns1)
	   isym(ns1)=7		!open circles for normalisation point
	   symsiz(ns1)=4.0
	   icol(ns1)=12		!red for norm point
    	   ijoin(ns1)=-1		!points not joined
	   isetcol=1
	endif
	if((nmod.eq.26.or.nmod.eq.27).and.ip1.eq.1.and.
     &    ifitmode.eq.3) then
	   if(iopt1.lt.1) iopt1=1
	   if(alpha) then
	   print 31,iopt1
31	   format(
     &   ' (1) Plot separate curves with their estimated maxima',/,
     &   ' (2) Plot single curve with points normalised by maxima',/,
     &   ' (3) No more plots',/,
     &   ' Option number [',i2,'] = ')
	   call INPUTi(iopt1)
	   else
		title='Plot method'
		strings(1)=
     &	  'Plot separate curves with their estimated maxima'
		strings(2)=
     &	  'Plot single curve with points normalised by maxima'
		strings(3)='No more plots'
		nval=3
		jline=iopt1		!default
		nhelp=0
		helps(1)=' '
		call POPMENU(nxlo,nylo,nyhi,strings,nval,ictm,icfm,icupm,
     &	      ibkm,title,helps,nhelp,jline,ans,ival)
		iopt1=jline
	   endif
	   if(iopt1.eq.3) then
 		call VIDEOMOD(3)	 !utility lib- back to alpha mode
		alpha=.true.
		idest=1
c restore original values
		ncurvc=ncurvc1		!= so input value not altered
		do i=1,ndimc
		   icurvc(i)=icurvc1(i)	!so icurvc in call not altered
		enddo
		RETURN
	   else if(iopt1.eq.2) then
c Make new xval,yval, xcal,ycal in j=nsfit+1
c Could also just normalise y but keep as sep curves so can have diff colours
c for each
		j1=nsfit+1
		j2=j1+1	!'column' for means
		if(.not.twodisp) then
		   ndsav=ncurvd	!save to restore option 1
		   idsav1=icurvd(1)	!save to restore option 1
		   idsav2=icurvd(2)	!save to restore option 1
		   idsav3=icurvd(3)	!save to restore option 1
		endif
		mode=0		!do normalisation
		call NORMDAT(j1,xval,yval,weight,theta,mode,
     &	 ncurvd,icurvd,ndat,icurvw,logyfit,juse,nsfit,ymax1,
     &	 ndv1,ndimd,niobs,njset,ndth,pon())
		ymax=1.1*ymax1		!new maximum for VPLOT for norm data
		ytic=0.2
		isym(j1)=-3		!filled squares for indiv points
		symsiz(j1)=2.0
		ijoin(j1)=-1		!points not joined
		isym(j2)=-7		!filled circles for means
		symsiz(j2)=3.0
    		ijoin(j2)=-1		!points not joined
		isym(3)=7		!open circles for normalisation point
		symsiz(3)=3.0
		icol(3)=12		!red for norm point
    		ijoin(3)=-1		!points not joined
c		isdev=1
	      jrep=10	 	!for setbigj
		if(j1.gt.10) then
		   j0=1+mod(j1-1,jrep)	!so j=11-20, 21-30 all -> j0=1-10 etc
		   icol(j0)=10
		else
		   icol(j1)=10		!green for sep points
		endif
		if(j2.gt.10) then
		   j0=1+mod(j2-1,jrep)	!so j=11-20, 21-30 all -> j0=1-10 etc
		   icol(j0)=12
		else
		   icol(j2)=12		!red for means
		endif
		isetcol=1

c    Now calculated curve -need only set ymax=1 and otherwise calculate
c     with fitted values of K and nH (nither j nor iset needed in ycalcv
c	in this case)
c============prev bit->normdat
		ncalc=ncal(1)
		ival=2	!so values normalised by ymaxj(j) in ycalcv
		do i=1,ncalc
		   xcal(i,j1)=xcal(i,1)
		   Xv=xcal(i,j1)
		   Ycal(i,j1)=YCALCV(kmax,theta,xv,ival,j,setx,njset)
               if(iostat.ge.3.and.iostat.le.9) return

            enddo
		ival=1		!reset, in case
		if(.not.twodisp) then
		   ncsav=ncurvc	!save to restore option 1
		   icsav=icurvc(1)	!save to restore option 1
		endif
		twodisp=.true.	!flag that ndsav etc defined
		ncurvc=1		!one calc curve only
		icurvc(1)=j1
		ncal(j1)=ncalc
		iline(j1)=0		!curve continuous
	   else if(iopt1.eq.1) then
		if(twodisp) then
		   ncurvd=ndsav		!restore option 1
		   icurvd(1)=idsav1	!restore option 1
		   icurvd(2)=idsav2	!restore option 1
		   icurvd(3)=idsav3	!restore option 1
		   ncurvc=ncsav		!restore option 1
		   icurvc(1)=icsav	!restore option 1
		endif
	   endif
	   iopt1=iopt1+1
	endif
c
	if(nodata) then
	   goto 271
	else
	   goto 30		!plot it
	endif
c
	end



