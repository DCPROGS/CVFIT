	subroutine CHDISP(Xobs,yval1,nyval,ILOG,ISCAL,IJOIN,iline,
     & XMIN,XMAX,YMIN,YMAX,XTIC,ytic,xlo,xhi,ylo,yhi,
     & xcross,ycross,titlex,titley,title1,itit,ilabel,inumx,inumy,
     & ncurvd,ndat,icurvd,ncurvc,ncal,icurvc,idiskq,
     & ntx,nty,itx,ity,theta,fitted,doframe,idest,
     & xaxis,xsav,nbin,nbw,xwbase,nbdec,logt,ibk)
c
c Version of CVDISP for display of histogram in Xobs
c Can have different bins for diff data sets?- VHIST allows this so
c allow option
c Contrary (I think) to rule book, cannot use Xobs,Yobs directly in call
c to VPLOT1 (OK for set 1 only!), so must copy them to Xval,Yval
	real XVAL(0:511,10),YVAL(0:511,10)
	real XCAL(2048,10),YCAL(2048,10)
	real theta(20)
c for histo data
	integer ndat(10),icurvd(10),ijoin(10)
c for calc curves
	integer ncal(10),icurvc(10),iline(10)
	dimension XOBS(100,10)
	dimension YVAL1(1024)
	dimension FREQ(510),XAXIS(510),xsav(5,4)
c for data
c	dimension ndat(ndimd),icurvd(ndimd),isym(ndimd),ijoin(ndimd)
c for calc curves
c	dimension ncal(ndimc),icurvc(ndimc),iline(ndimc)
	integer*4 videotyp
	LOGICAL pon,debug,fitted,doframe,draft,autplt
	character*40 titlex,titley
	character*64 title1
	logical constr,fline		!for COMMON
	integer jmiss(10),nj(10),juse(10)
	logical logx,logy,sqrty,slock,caplock,landscap,logt,sbin
	character*1 ans,UC
	logical discprt
	common/dp/discprt
	COMMON/cols/icol(100),mono
	COMMON/BLOCK1/constr,nset,nfit,nsfit,
     & Xv,alfa,kmax,ncomp,nmod,fline,nomit,jomit(20),jset,ifitmode
	COMMON/BLOCK2/ CASTAR,X1,X2,iequiv,ip1,ip2
	COMMON/JLOGOS/t1c(10),t2c(10),t1v(10),t2v(10),xoff1,y1v,y2v,
     & y1c,y2c
c
c define function
	pon()=slock()
	debug()=caplock()
c
101	format(a1)
104	FORMAT(I8)
c
	sbin=.false.
	if(idest.eq.0) goto 60		!first time
	print 62
62	format('&Same bins [N] ? ')
	read 101,ans
	if(UC(ans).eq.'Y') then
	   iscal=0
	   sbin=.true.
c	   goto 61		!goto 61 for 'same data'
	endif
c	if(.not.sbin) NBIN=-1		!TO SIGNAL THAT NOT YET SET
c Decide which sets to be displayed if > 1
c Note that if nset=-1 this means that a single set of data has been
c produced in CVIN1, in yval1,nyval rather than in XOBS
60	continue
	ncurvd=1
	icurvd(1)=1
	if(nset.eq.1.or.nset.eq.-1) goto 27
24	print 22
22	format(' Display:',/,
     & ' (1) Selected data sets',/,
     & ' (2) All data sets',/,
     & '  Option number [1] = ')
c     & '',/,
	read 104,i
	if(i.eq.0) iopt=1
c
	print 14
14	format('&Number of sets to be displayed= ')
	read 104,ncurvd
	do 271 j=1,ncurvd
	   print 272,j
272	   format('&(',i3,'): data set no = ')
	   read 104,icurvd(j)
271	continue
c
27	continue
	if(.not.sbin) logt=.false.
	if(sbin) goto 303	!same bins so XAXIS,XWBASE as before
	print 2092
2092	format(
     & ' (1) Distribution of log durations',/,
     & ' (2) Distribution of durations- set bins manually',/,
     & ' (3) Distribution of durations- 20 bins',/,
     & ' (4) Distribution of durations- 40 bins',/,
     & ' (5) Distribution of durations- 60 bins',/,
     & ' (6) Skip histogram',/,
     & ' Option number [2] = ')
c     & '',/,
	read 102,iopt
102	format(i8)
	if(iopt.eq.0) iopt=2
	if(iopt.ne.6) goto 40
90	nbin=0
	RETURN		!skip histo
40	logt=iopt.eq.1
	nbw=1
303	continue
c
c First get bin frequencies for all sets to be displayed
	do 150 j1=1,ncurvd
	j=icurvd(j1)
      print 41,j
      if(pon()) write(7,41) j
      if(discprt) write(8,41) j
41	format(' Data set number ',i4)
c if nset=-1 then yval1,nyval already define in CVIN1
	if(nset.ge.1) then
	   nyval=nj(j)
c Copy data to Yval
	   do 221 i=1,nyval
221	   yval1(i)=Xobs(i,j)	!one-D array for varv1 (orig data in Xobs)
	endif
c -show mean etc
	call VARV1(Yval1,nyval,YBAR,VARY,ymin,ymax)
	print 413,NYVAL,YBAR,sqrt(vary),ymin,ymax
	if(pon()) write(7,413) NYVAL,YBAR,sqrt(vary),ymin,ymax
	if(discprt) write(8,413) NYVAL,YBAR,sqrt(vary),ymin,ymax
413	FORMAT(' Mean and SD of ',I5,' values= ',G13.6,' +/-',g13.6,
     & /,' Range from ',g13.6,' to ',g13.6)
	if(sbin) goto 304		!?
	if(j1.gt.1) then
	   print 151
151	   format(' Same bins as for first data set [Y] ? ')
	   read 101,ans
	   if(UC(ans).ne.'N') then
		do 11 i=1,nbin+1
11	 	xval(i,j)=xval(i,1)
		ndat(j)=ndat(1)
		goto 123	!sort yval to get FREQ, and set Yval
	   endif
	endif
304	continue
c Set bin widths
	call SBINCH(XAXIS,nyval,ymax,ymin,xsav,nbin,nbw,xwbase,
     & nbdec,logt,iopt,sbin)
	ndat(j)=nbin
	do 1 i=1,nbin+1
1	 xval(i,j)=xaxis(i)
	if(j1.eq.1) then
	   xmin=xval(1,j)
	   xmax=xval(nbin+1,j)
	else if(j1.gt.1) then
	   if(xval(1,j).lt.xmin) xmin=xval(1,j)
	   if(xval(nbin+1,j).gt.xmax) xmax=xval(nbin+1,j)
	endif
c
123	continue          !jump here if later sets have same bins
c Sort Yval1 into bins
	call SORTBCH(YVAL1,FREQ,XAXIS,nyval,logt,nbin,xwbase,
     & ymax,flo,fhi)
c     & fmax,flo,fhi)
c
c Xaxis and Freq now defined for the current set. Copy them to Xval,Yval
c (as, eg, in SCVDU)
c For Toshiba VIDEOTYP()=18 for graphics mode, =3 for alpha mode
c histos:
c Now define input data for the current histo (data set #j; see do 150 loop
c above) with current xaxis and freq
	do 21 i=1,nbin
	 yval(i,j)=freq(i)
21	continue
	if(j1.eq.1) then
	   fmax=ymax      !from sortbch
	else if(j1.gt.1) then
	   if(ymax.gt.fmax) fmax=ymax
	endif
	if(flo.gt.0.0001) then
	   yval(0,j)=flo
	   xval(0,j)=xval(1,j)-xwbase
	   if(logt) xval(0,j)=xval(1,j)*(xval(1,j)/xval(2,j))
	   if(j1.eq.1) then
		xmin=xval(0,j)
	   else if(j1.gt.1) then
		if(xval(0,j).lt.xmin) xmin=xval(0,j)
	   endif
	endif
	if(fhi.gt.0.0001) then
	   yval(nbin+2,j)=fhi	!NB Yval(nbin+1,j) not used
	   xval(nbin+2,j)=xval(nbin+1,j)+xwbase
	   if(logt) xval(nbin+2,j)=xval(nbin+1,j)*
     &	(xval(nbin+1,j)/xval(nbin,j))
c	  keep ratio for hi bin as for other bins, so same width on log scale
	   xval(nbin+2,j)=xval(nbin+1,j)*(xval(nbin+1,j)/xval(nbin,j))
	   if(j1.eq.1) then
		xmax=xval(nbin+2,j)
	   else if(j1.gt.1) then
		if(xval(nbin+2,j).gt.xmax) xmax=xval(nbin+2,j)
	   endif
	 endif
c Set other things that need to be specified separately for each data set
	ijoin(j)=0
c
150	continue
c END of loop for each data set
c Xval and Yval defined for all sets to be displayed: now display them
c
	fmin=0.
c xmin,xmax now set during loop for each set (will not include flo,fhi
c if these are zero for all sets)
c	xmin=xaxis(1)
c	xmax=XAXIS(nbin+1)+(xaxis(nbin+1)-xaxis(nbin))	!room for FHI
	fmax=fmax*1.3
	fmax=float(1+IFIX(fmax))	!ROUND UP
	ftic=0.5*10.0**IFIX(ALOG10((fmax-fmin)*0.5))
2102	if((fmax-fmin)/ftic.gt.10.) goto 2101
	ftic=ftic/2.
	goto 2102	!ensure at least 10 tics so get 2 numerical labels
2101	XTIC=2.*XWBASE	!number every 10th bin
c# For Lahey/Hgraph xtic refers to major, not minor tics so make
c bigger ( say 2-fold rather than 5-fold for now?- 5-fold faster!)
	xtic=xtic*5.
	if(.not.logt) ftic=ftic*5.	!smaller tics for sqrt scale for now
61	continue
	if(videotyp().ne.18) call GINO  !INIPLT(0,.false.,1.0)
	autplt=.false.		!for VHIST
	draft=.false.		!for VHIST
	doframe=.true.		!for VHIST
	landscap=.true.		!for VHIST
	isval=0	!no arrow
	xcross=0.0
	ycross=0.0
	inumx=-1	!X axis in fixed format if log scale
	ndimd=10
	ndimc=10
	iask=-2 !do not ask before leaving display; do not leave graph on screen
	lt2=2       !short dash for lo,hi bins
	ilog=0
	if(logt) ilog=6		!plot sqrt(y) vs log(t)
	iscal=0			!use input xmin,xmax etc
	xlo=-1
	if(logt) iscal=1	!scale internally 1st time
	ntx=5		!if not logt set initial input values
	nty=5
	itx=1
	ity=1
	xw=0.0	!to suppress 'per 10 ms' etc
	fitted=.false.
c curves:
	ncurvc=0		!no calc curve yet
	if(debug()) print 77,xmin,xmax,fmin,fmax,xtic,ftic
77	format(' xmin,xmax,fmin,fmax,xtic,ftic=',/,6g13.6)
	call VHIST5(XVAL,YVAL,NDAT,icurvd,ncurvd,ijoin,ndimd,
     & XCAL,YCAL,NCAL,icurvc,ncurvc,iline,ndimc,ILOG,ISCAL,
     & XMIN,XMAX,fMIN,fMAX,XTIC,fTIC,xcross,ycross,ntx,nty,itx,ity,
     & XLO,XHI,YLO,YHI,xw,lt2,inumx,inumy,
     & titlex,titley,ilabel,doframe,idiskq,
     & autplt,draft,itit,title1,csize,ifont,landscap,fitted,
     & theta,ifitype,ncomp,sval,isval,iask,iver)
c
	ibk=icol(71)
	idest=1
	print 231
231	format(' Another display [N] ? ')
	read 101,ans
	if(UC(ans).eq.'Y') goto 61
c	goto 999		!return
c
c Now bit done after fit
c10	continue
c
999	continue
	idest=0
	RETURN
c
	end



