	PROGRAM CVFIT
c
======CVFITG does not yet have hesmat3, or changes in cvfit, cvfout1, cvsimp, 
======and likint1 consequent on change in way error sd incorporated

c GENERAL CURVE FITTING PROGRAM
c Lahey V5.n version of CVFIT2
c NEW ADDITIONS?
c (1) Option to read earlier data from PLOTQ rather than CVDAT?
c (2) ASCII i/o
c (3) Make arrays bigger (and allocatable)
c
c 03/25/95 04:40pm Version =1003 in which xobs, yobs are allocatable and
c there is no restriction on number of sets or of obs per set.  This involves
c a reorganization of the CVDAT.DAT file, so this version should be able to
c read old CVDAT, but can write data only to CVDAT file that has iver=1003
c or greater in record #1 -too complex to mix data with earlier iver?
c   V1003 further modif 03/30/95 09:56pm to take up to 100 parameters (set to
c 100 in main prog at present, but most subroutines now general). NB the array
c titlep() (param names) is at present kept in CVDAT.DAT but not used, so
c remove it from Version 1003.  Also increase titlep() to 20 char.
c
c Modif 09/17/94 07:24pm to add
c  (1) Lines with common Y intercept
c  (2) fix ratio of 2 parameters
c  (3) infil*33 replaces infil*14 to hold whole path
c Modif 06/30/93 02:30pm to cope with old data files that have different iver
c in same data file -check if nset=1-10 when read as if old, and set iver from
c this (overwrites the iver read from record #1)
c Modified 04/10/92 11:39am for relative potency fits (ifitmode=4) or
c dose ratios constrained by Schild eqn (ifitmode=5) -fit and
c  display specified eqn, f(c), (e.g. polynomial or Hill eqn) both
c  (a) separately to each specified set and
c  (b) constrained to be parallel on log(x) scale, so f(c;th) fits 1st set
c f(r1*c;th) fits 2nd set, f(r2*c;th) fits 3rd set  . . . where r1,r2 are the
c nset-1 independent potency ratios relative to set 1 and th = parameters
c for the curve fitted to set 1 (e.g. Hill eqn, polynomial etc)
c   Also option to constrain r(j)=1 + xB(j)/Kb where
c xB(j)=conc of competitive blocker for set j (specified as the 'set variable')
c and Kb is to be estimated (so there is only 1 extra parameter rather
c than nset-1)
c Version 07/30/91 02:52pm that asks which equation to be fitted AFTER the
c initial display. Also allows different eqn to be fitted to each of
c simultaneously display curves; this means that either (a) Ycal, Xcal have
c to be calc during loop for each set (solution adopted here), or (b) a record
c has to be kept of which eqn was fitted to each set to allow xcal, ycal to
c be calc in CVDISP (as it was before).
c
c Uses new data storage defined in CVDISK
c
c MODIFIED 07/09/91 04:26pm so that ifitmode=2 (fit same equation separately
c to several data sets) no longer has to fit all sets
c Lahey version 12/06/89 05:11pm
c (1)SIMULTANEOUS WEIGHTED LS FIT FOR NSET SETS OF OBS WITH
c	NJ(J) OBS IN THE J'TH SET. DIMENSIONS OF WEIGHTS W(*,J)
c	SAME AS THOSE FOR Y(*,J)**2 SO SSD DIMENSIONLESS
c (2)Y(I,J),XOBS(I,J),W(I,J) ARE FOR I'TH OBSERVATION IN J'TH
c	DATA SET
c
c
	ALLOCATABLE:: XVAL,YVAL,XCAL,YCAL
	real XVAL(:,:),YVAL(:,:)		!for VPLOT
	real XCAL(:,:),YCAL(:,:)		!for VPLOT/VHIST
c Problem to allocate xobs,yobs, w for now because they were in common
	ALLOCATABLE:: Xobs,Yobs,w
	real Xobs(:,:),Yobs(:,:),w(:,:)
c	real YOBS(100,10),XOBS(100,10),W(100,10)
c	dimension XVAL(2048,10),YVAL(2048,10)
c	real xcal(2048,20),ycal(2048,20)
	ALLOCATABLE:: titled,jmiss,juse,nj,setx
	character*60 TITLED(:)
	integer jmiss(:),juse(:),nj(:)
	real setx(:)
c	character*60 TITLED(10),titlef	!data set names,file name
c	real setx(10)
c	integer JMISS(10),juse(10),NJ(10)
	real*4 ymaxj(100)
	character*60 titlef		!file name
	character string*50,text*18
	character ascfil*33
	character*60 titw(5)		!names of weighting methods
c	real VAR(20)
	real Sdev(100)
	real*4 theta(100),thetsav(100,11)
	character*20 TITLEP(100)		!names of params for fit
	real yval1(1024)
	real XAXIS(510),xsav(5,4)	!for histo display
	integer JFIX(100),jomit(20)
	character*1 ans,UC
c	character*1 getch		!for debug
c	integer*2 ktype		!for debug
	integer*2 videotyp
	character ndev*2,infil*33
c Allocate other arrays for VPLOT
	ALLOCATABLE:: ndat,icurvd,icurvw,isym,ijoin,symsiz
	ALLOCATABLE:: ncal,icurvc,iline
	integer ndat(:),icurvd(:),icurvw(:),isym(:),ijoin(:)
	integer ncal(:),icurvc(:),iline(:)
	real symsiz(:)				!for data
c	dimension ndat(10),icurvd(10),isym(10),ijoin(10)	!for data
c	real symsiz(10)				!for data
c	dimension ncal(20),icurvc(20),iline(20)		!for calc curve
	integer ncal2(10),icurvc2(10),iline2(10)		!for histo disp
	LOGICAL CONSTR,FLINE,SETVAR,defolt,fitted,dcurve,logt,badset
	logical xgt0,noerr,fixmax,fix0,readini,alpha,allocated,logy
	character*40 titlex,titley
	character*44 title1
	logical discprt,pon,slock,hdisp,fopen,printed,logx,present
	logical debug,caplock
	character*11 cdate,ctime
	character cnum*11
c For mouse
	logical mouse_on,mous_set
	common/mousval/mouse_on,nbutton
c Common block to increase raster vector buffer size for Hgraph
c Maxnum (=2nd dimension of vect) set to 8192 below (rather than default
c which is 1024)- reduces disc work making temp files (eg 05158741.n) which
c have to be made for every maxnum vectors to be plotted (do not alter
c curnum or 1st dimension of VECT (=6))
	integer*2 maxnum,curnum,vect(6,8192)
	COMMON/CHVBLK/maxnum,curnum,vect
	character*40 mtitle*40,filnam*32,prtport*4	!for WINPRINT
	common/dpp/filnam,prtport,ndisc,icol,mtitle !for WINPRINT,ENDPRINT,DISCNUM
	common/dp/discprt
c= xobs, yobs no longer in common
c=	common/data/Xobs,Yobs
	COMMON/LIKBLK/kmax2,kfit,elmax,em,delmin,ndisp,npar,iw,sres
	common/logblk/debug,defolt
	COMMON/BLOCK1/constr,nset,nfit,nsfit,
     & Xv,alfa,kmax,ncomp,nmod,fline,nomit,jomit,jset,ifitmode
c=	COMMON/BLOCK2/ castar,setx(10),X1,X2,iequiv,ip1,ip2
c====setx no longer in common/block2/
	COMMON/BLOCK2/ castar,X1,X2,iequiv,ip1,ip2
c====w(i,j) no longer in common!
c=	COMMON/SD/w		  !for CVDISP2,VPLOT2,SSDCV,YCALCV
	common/potrat/jfirst,iset,kmax1			!for pot ratios
	common/max/ymaxj		!for CVDISP and YCALCV
	logical logyfit,norm,nodata
	COMMON/BLOCK3/logyfit,norm,xnorm
c (kmax=actual no of param to be estimated; kmax1=number before extra param
c  added for ifitmode=4,5)
c
	external SSDCV,ELFUNCV
c
c Oct 86- two new variables added to common/block2/
c IEQUIV=0 for equiv subunits; =1 for non-equiv (prev signalled via castar)
c IP1,IP2:
c  (1) mods 10,11 IP1=1 if param are K1,K2; =2 if param are K1, K2/K1.
c  (2) mods 9,11  IP2=1 if param is K2 or (K2/K1 if IP1=2)
c		  IP2=2 if param is k(-2) (must have ip2=1 if ip1=2)
c  (3) mods 9,11  IP2 set negative if k(+2) is a param  to be estimated
c
C
c define function
	pon()=slock()
	debug()=caplock()
c
101	format(a1)
c	iver=1002		!new version with file title
c	iver=1003		!new version with alloc arrays and new CVDAT org.
c NB iver is initially the version number for the CVDAT being read in
c (may be 1001, 1002, 1003,..).  Writing is done in CVDISK only in the
c format for iver=1003 so set iver=1003 there (and whole CVDAT must be the
c SAME version!)
	ndth=100		!dimension of theta(),jfix(),titlep(),sdev()
	call SETMOUSE()		!define values in common\mousval\ (in IVLIB)
	maxnum=8192		!for Hgraph common/chvblk/
	idest=0
	iflag=0	!for GETEQN
	lflag=0	!for log scale question
	iesc=0	!for DISCNUM abort
	idiskq=-1	!until plot queue disc defined
	isepfit=1	!in case nmod is not 1 or -1
	fitted=.false.
	dcurve=.false.
	fixmax=.false.
	fix0=.false.
c In scan, colours are set in all subroutines via
c common/dialval1/ictd,icfd,icbd  !colours for dialog box #1 (text,frame,back)
c but this not used in CVFIT at present, and ict=11 used only in main prog.
	ict=11	!text colour in dialog box 1
c
      filnam='CVFIT.PRT'
	call WINPRINT	!print-file control
c      OPEN(unit=7,file='PRN',iostat=nerr)		!open printer
      OPEN(unit=7,file=prtport,iostat=nerr)		!open printer
	print 1011
c	if(pon()) write(7,1011)
	if(discprt) write(8,1011)
1011	format(' CVFIT PROGRAM',/)
	call DATE1(cdate)		!DC subroutine
	call TIME(ctime)
	print 1,cdate,ctime(1:8),mtitle
	if(pon()) write(7,1) cdate,ctime(1:8),mtitle
	if(discprt) write(8,1) cdate,ctime(1:8),mtitle
1	format(' Date of analysis: ',a11,/,' Time of analysis: ',a8,/,
     & '   Machine = ',a40)
	print 399
399	format(
     & ' SCROLL LOCK on for printing'/
     & ' CAPS LOCK on for debugging')
c
c
	readini=.false.
	INQUIRE(file='CVFIT.INI',exist=present,flen=len)
	if(present.and.len.gt.0) then
	   if(len.eq.256) then
          OPEN(unit=19,file='CVFIT.INI',status='UNKNOWN',
     &    access='DIRECT',form='UNFORMATTED',recl=256)
	    read(19,rec=1) ndev,ifile1,ifitmode,ilog,idiskq,titlef,infil,
     &	niobs,njset,ascfil
	   else if(len.eq.128) then
          OPEN(unit=19,file='CVFIT.INI',status='UNKNOWN',
     &    access='DIRECT',form='UNFORMATTED',recl=128)
	    read(19,rec=1) ndev,ifile1,ifitmode,ilog,idiskq,titlef,infil,
     &	niobs,njset
	    CLOSE(unit=19)
	   else if(len.eq.64) then		!old version
          OPEN(unit=19,file='CVFIT.INI',status='UNKNOWN',
     &    access='DIRECT',form='UNFORMATTED',recl=64)
	    read(19,rec=1) ndev,ifile1,ifitmode,ilog,idiskq
	    CLOSE(unit=19)
	   endif
	   readini=.true.
	endif
c
c Allocate arrays to size they were in old prog for the moment
c NB need more columns than nset in CVDISP -so use nset+5 to ne on safe side?
12	continue
	print 10,niobs
10	format(/,
     & ' Maximimum number of observations per set [',i5,'] = ')
	call INPUTi(niobs)
	print 11,njset
11	format(
     & '&Maximimum array size for sets (# of sets +5) [',i5,'] = ')
	call INPUTi(njset)	!need extra cols in CVDISP
	if(njset.lt.10) njset=10	!(needed to read old files)
c	ndv1=2048
c	ndimd=10
	ndv1=niobs
	ndimd=njset
c=	ALLOCATE(xval(ndv1,ndimd),yval(ndv1,ndimd))	!later
	ndc1=2048
	ndimc=20
	if(njset.gt.20) ndimc=njset
c	niobs=100
c	njset=10
	ALLOCATE(xobs(niobs,njset),yobs(niobs,njset),w(niobs,njset))
	ALLOCATE(titled(njset),jmiss(njset),juse(njset),
     & nj(njset),setx(njset))
	ALLOCATE(ndat(ndimd),icurvd(ndimd),icurvw(ndimd),isym(ndimd),
     & ijoin(ndimd),symsiz(ndimd))
c=	ALLOCATE(xcal(ndc1,ndimc),ycal(ndc1,ndimc))	!later
	ALLOCATE(ncal(ndimc),icurvc(ndimc),iline(ndimc))
c
	ndev='  '		!before disk # defined
	fopen=.false.	!data file not open yet
	ilabel=-1	!internal axis labels unless otherwise specified
2	continue
      print 39
39	format(' (1) Curve fitting',/,' (2) Histogram display',/,
     & ' Option [1] = ')
	call INPUTi(i)
c	read 104,i
104	format(i8)
	if(i.le.0) i=1
	hdisp=i.eq.2
c	nfirst=1
	nomit=0
c
	call CVIN2(IW,setvar,aval,bval,defolt,ndev,titw,nmiss,
     & titlef,titled,titlex,titley,ilabel,yobs,xobs,w,idest,
     & fopen,ncalc,ncurvc,icurvc,ncal,iline,hdisp,yval1,nyval,printed,
     & readini,ifile,ilog,idiskq,ifitsav,iesc,ndisc,infil,ifile1,
     & niobs,njset,ndimd,ndimc,nj,jmiss,setx,iver,ascfil,nodata)
C/===where next for 'nodata' case; need to get equation/get params/
c/===then calc values and display as for dcurve -skip first display of data!
c
	if(iesc.eq.98) then	!nset too big -reallocate arrays
	   njset=nset+5		!new default size
c=	   DEALLOCATE(xval,yval,xcal,ycal)
	   DEALLOCATE(xobs,yobs,w)
	   DEALLOCATE(titled,jmiss,juse,nj,setx)
	   DEALLOCATE(ndat,icurvd,icurvw,isym,ijoin,symsiz)
	   DEALLOCATE(ncal,icurvc,iline)
	   goto 12
	endif
	if(iesc.eq.99) goto 99	!DISCNUM abort
	if(.not.allocated(xval)) then
	   ALLOCATE(xval(ndv1,ndimd),yval(ndv1,ndimd))
	   ALLOCATE(xcal(ndc1,ndimc),ycal(ndc1,ndimc))
	endif
c Define nsfit=number of sets to be fitted, and juse()=pointer to the set
c numbers (quicker than using jmiss always)
	nsfit=0
	do j=1,nset
	   if(jmiss(j).eq.0) then
	      nsfit=nsfit+1
	      juse(nsfit)=j
	   endif
	enddo
c Now initial display of data before fitting
	if(nodata) then
	   goto 33
	endif
	fitted=.false.
	dcurve=.false.
3062	continue		!return here after fitted for final display
	if(.not.hdisp) then
	  call CVDISP(Xobs,Yobs,w,ISYM,ILOG,ISCAL,IJOIN,iline,
     &  XMIN,XMAX,YMIN,YMAX,XTIC,YTIC,ixlo,ixhi,iylo,iyhi,y0,yinf,
     &  xcross,ycross,titlex,titley,title1,itit,ilabel,inumx,inumy,
     &  xval,yval,ncurvd,ndat,icurvd,ncurvc,ncal,icurvc,Xcal,Ycal,
     &  symsiz,idiskq,ntx,nty,itx,ity,theta,thetsav,fitted,doframe,
     &  dcurve,nodata,iw,idest,readini,iflag,ndv1,ndimd,ndc1,ndimc,
     &  jmiss,nj,niobs,njset,icurvw)
	  logx=ilog.eq.1.or.ilog.eq.3.or.ilog.eq.6	!needed for DEFYcal
	  logy=ilog.eq.2.or.ilog.eq.3
	  if(idest.eq.1) goto 320	!after display with fit, calc likelihood int?
	  if(idest.eq.99) goto 99	!another run?
c       NB if idest=0 (fit next...) then still in graphics mode
	else
	   ncurvc2=ncurvc
	   do 411 i=1,10
	    ncal2(i)=ncal(i)
	    icurvc2(i)=icurvc(i)
411	    iline2(i)=iline(i)
	   call CHDISP(Xobs,yval1,nyval,ILOG,ISCAL,IJOIN,iline2,
     &   XMIN,XMAX,YMIN,YMAX,XTIC,ytic,ixlo,ixhi,iylo,iyhi,
     &   xcross,ycross,titlex,titley,title1,itit,ilabel,inumx,inumy,
     &   ncurvd,ndat,icurvd,ncurvc2,ncal2,icurvc2,idiskq,
     &   ntx,nty,itx,ity,theta,fitted,doframe,idest,
     &   xaxis,xsav,nbin,nbw,xwbase,nbdec,logt)
	   goto 99	!no fitting facility yet
	endif
c
c Do fit
c Fit set #1 (first, if more than one to be fitted); if only one to be
c fitted jset is defined in CVIN (ifitmode=1); if several sets to
c be fitted separately (ifitmode=2,3) then those to be missed have jmiss(j)=1
c Find 1st set to be fitted=#jset
33	continue			!jump here if nodata
	alpha=VIDEOTYP().eq.3		!alpha mode
	if(ifitmode.ge.2) then
c	   j=0
c33	   j=j+1
c	   if(jmiss(j).eq.1) goto 33
c	   jset=j
	   mset=1		!to count sets for ifitmode=2
	   j=juse(1)
	   jset=j
	endif
	jfirst=jset
	if(ifitmode.ge.3.and.ifitmode.le.5) then	   !print titles of all sets
	   if(alpha) then
            print 321
	      if(iver.ge.1002) print 322,titlef
		do j=1,nset
		   if(jmiss(j).ne.1) then
		  	text='FIT TO DATA SET # '
		      print 1092,text,j,titled(j)
		   endif
		enddo
	   else
		n=1
		do j=1,nset
		   if(jmiss(j).ne.1) then
			string(n:n+1)=char(j+48)//','
			n=n+2
		   endif
		enddo
		string(n:n)=' '		!remove final comma
		n=nblank(string)
	      call WDIALOG(1,'SIMULTANEOUS FIT TO SETS '//string(1:n),ict)
	   endif
         if(pon()) write(7,321)
         if(discprt) write(8,321)
321	   format(/,' SIMULTANEOUS FIT TO THE FOLLOWING SETS')
	   if(iver.ge.1002) then
	      if(pon()) write(7,322) titlef
	      if(discprt) write(8,322) titlef
322		format(' File: ',a60)
	   endif
	   do 32 j=1,nset
	    if(jmiss(j).eq.1) goto 32
	    text='FIT TO DATA SET # '
	    if(pon()) write(7,1092) text,j,titled(j)
          if(discprt) write(8,1092) text,j,titled(j)
32	   continue
	endif
	SSDsep=0.0		!initialise for case where summed for each set
	nfittot=0		!ditto
c
31	continue		!return to fit next set
	alpha=VIDEOTYP().eq.3		!alpha mode
	if(ifitmode.eq.1.or.ifitmode.eq.2) then
	   if(iver.ge.1002) then
	      if(pon()) write(7,322) titlef
	      if(discprt) write(8,322) titlef
c322		format(' File: ',a60)
	   endif
	  text='FIT TO DATA SET # '
	  if(dcurve) text='CURVE + DATA SET #'
	  if(nodata) text='CALCULATED CURVE #'
	  if(alpha) then
	     print 322,titlef
	     print 1092,text,jset,titled(jset)
	  else
	     call INTCONV(jset,cnum)
	     string=titled(jset)(1:50)
	     call WDIALOG(1,'File: '//string,ict)
	     call WDIALOG(1,text//charnb(cnum),ict)
	  endif
	  if(pon()) write(7,1092) text,jset,titled(jset)
        if(discprt) write(8,1092) text,jset,titled(jset)
1092	  FORMAT(/,' ',a18,i3,/,3x,a60)
	endif
	if(iflag.ne.2) call GETEQN(IW,titw,titlep,xgt0,iflag,
     &  jmiss,juse,njset,ndth,imod,nodata,iesc)
c If no data, get equation but skip fit
c For ifitmode=4,5: kmax redefined in GETEQN as is kmax1=orig no of param
	if(iesc.eq.99) goto 99	!fitting abort
	call CVSIMPg(kfit,jfix,nfix,defolt,titlep,fitted,dcurve,
     & nodata,theta,thetsav,nev,sres,smin,noerr,iflag,isepfit,idest,
     & fixmax,fix0,Xobs,yobs,w,nj,juse,setx,niobs,njset,ndth,
     & badset)
	if(badset) then
	   if(ifitmode.eq.1) then
		goto 99	!only 1 set and it is bad
	   else if(ifitmode.eq.2) then	!skip to next set
		mset=mset+1
		if(mset.gt.nsfit) goto 35	!display
		j=juse(mset)		!next set # to be fitted
		jset=j
		if(iflag.eq.0) iflag=1		!to tell GETEQN this is not 1st set
		goto 31		!fit next set
	   endif
	endif
c
	if(idest.eq.201) goto 201	!no errors
c
	kfit=kmax-nfix		!redefine (in case kmax redef for ifitmode=4,5)
c	if(debug()) then
	   print 600,kmax,nfix,kfit
600	   format(' kmax,nfix,kfit = ',3i8,/)
c	endif
	call HESMAT1(theta,sdev,nfix,jfix,SMIN,SSDCV,kmax,
     & kfit,100,Xobs,yobs,w,nj,juse,setx,niobs,njset)
c (100=dimension of theta)
C TYPE CALC LS VALUES
201	continue	!jump here for curve with given param
c  After separate fits (following ifitmode=4,5 run) save the SSD and number
c of params fitted, for test of dev from parallism (in CVOUT1)
	if(iflag.eq.2) then
	   if(isepfit.eq.1.or.isepfit.eq.2) then  !one simplex for each set
	      SSDsep=SSDsep+smin		!add values for each set
		nfittot=nfittot+nfit
	      kfitsep=nsfit*(kmax-nfix) !save number of param fitted for ditto
	   else if(isepfit.eq.3) then		!only one call to simplex
		SSDsep=smin
		nfittot=nfit
		kfitsep=kmax-nfix
	   endif
	endif
c
 	if(dcurve.or.nodata) call VIDEOMOD(3) !utility lib- back to alpha mode
	call CVOUT1(theta,noerr,iw,nev,smin,jfix,fitted,titlep,
     & sres,Sdev,titled,yobs,xobs,nj,jmiss,setx,niobs,njset,nodata,
     & dcurve,ndth)
c Calculate here the fitted curve for the model in use for the current set
c (jset is in COMMON) using DEFYcal (NB in NODATA case, ifitmode=2)
c
c	print 108
c108	format(/)
	if(lflag.eq.0) then	!ask once only
	   if(.not.readini) ilog=0
	   i=ilog+1
	   print 109,i
109	   format(' Initial display type:',/,
     &    ' (1) Display Y against X',/,
     &    ' (2) Display Y against log(X)',/,
     &    ' (3) Display log(Y) against X',/,
     &    ' (4) Display log(Y) against log(X)',/,
     &    ' (5) Display as Hill plot',/,
     &    ' Option number [',i2,'] = ')
	   call INPUTi(i)
	   ilog=i-1
	   logx=ilog.eq.1.or.ilog.eq.3.or.ilog.eq.4	!for DEFYCAL
	   logy=ilog.eq.2.or.ilog.eq.3
	endif
c	if(nodata.and.lflag.eq.0) then
c Now ask for xmin, xmax separately for each set in 'nodata' case
	if(nodata) then
	   print 20,jset,xmin,xmax
20	   format(
     &' Calculate curve ',i3,' between Xmin, Xmax [',
     &   g10.3,1x,g10.3,'] = ')
	   call INPUT2r(xmin,xmax)
22	   if(logx.and.((xmin.le.0.).or.(xmax.le.0.))) then
		call BELL(1)
		print 21
21		format(' Xmin, Xmax must be greater than 0 for log scale:',/
     &	,' Xmin, Xmax = ')
		call INPUT2r(xmin,xmax)
		goto 22
	   endif
	   xmin=xmin*1.000001	!so does not round down
	   if(logx) then
		xmin=alog10(xmin)		!for defycal
		xmax=alog10(xmax)
	   endif
	endif
	lflag=1	!scales now defined
c
	if(ifitmode.eq.1) then		!one pass only
	   if(nmod.eq.26.or.nmod.eq.27) then
		iset=1		!needed in commom for ycalcv
	   endif
	   call DEFYcal(xmin,xmax,ncalc,logx,theta,Xcal,Ycal,jset,xgt0,
     &    setx,njset,ndimc,ndc1)
	   goto 35		!display
	else if(ifitmode.eq.2) then	!calc curves for current pass (set #jset)
	   if(nmod.eq.26.or.nmod.eq.27) then
		iset=1		!needed in commom for ycalcv
	   endif
	   jval=jset
	   if(iflag.eq.2) then	!for separate fits
		jval=jset+nset	!sep fits in j=nset+1, nset+2,....,2*nset
		iline(jval)=2
		ncal(jval)=ncalc
		ncurvc=ncurvc+1
		icurvc(ncurvc)=jval
	   endif
	   call DEFYcal(xmin,xmax,ncalc,logx,theta,Xcal,Ycal,jval,xgt0,
     &    setx,njset,ndimc,ndc1)
c	   j=jset	!last set fitted -look for next
c34	   j=j+1
c	   if(j.gt.nset) goto 35	!display
c	   if(jmiss(j).eq.1) goto 34
	   mset=mset+1
	   if(mset.gt.nsfit) goto 35	!display
	   j=juse(mset)		!next set # to be fitted
	   jset=j
	   if(iflag.eq.0) iflag=1		!to tell GETEQN this is not 1st set
c	   text='FIT TO DATA SET # '
c	   print 1092,text,jset,titled(jset)	!say which set next on screen
	   goto 31		!fit next set
	else if(ifitmode.ge.3.and.ifitmode.le.5) then !one pass only -calc all curves
	   iset=0
c	   do 37 j=1,nset
c	    if(jmiss(j).eq.1) goto 37
	   do 37 m=1,nsfit		!go through all fited sets
	    j=juse(m)
	    iset=iset+1	!for common -needed for ifitmode=4,5 if sets omitted
	    jset=j		!for COMMON
	    jval=j
	    if(iflag.eq.2) then	!for sep fits (with common max) after ifit=4,5
		jval=j+nset		!sep fits in j=nset+1, nset+2,....,2*nset
		ncurvc=ncurvc+1
		icurvc(ncurvc)=jval
		iline(jval)=2
		ncal(jval)=ncalc
	    endif
	    call DEFYcal(xmin,xmax,ncalc,logx,theta,Xcal,Ycal,jval,xgt0,
     &    setx,njset,ndimc,ndc1)
37	   continue
c For nmod=26,27, when separate max fitted to each curve (ip1=1) calculate here
c  the total amplitude, measured from Y(0) (see ycalcv), for all components,
c  for set #j and store it in ymaxj(j) (e.g. used for normalising curve in
c  CVDISP) -actual max on graph = Y(0) + SUM(ymaxj(j))
c
	   if((nmod.eq.26.or.nmod.eq.27).and.ip1.eq.1) then
		iset=0
		do j=1,nsfit
		   iset=iset+1
		   j0=1+(iset-1)*ncomp	!=1,1+ncomp,1+2*ncomp for set 1,2,3...
		   ymaxj(j)=0.0	!total ymax for all components in set j
		   do i=1,ncomp
		      ymaxj(j)=ymaxj(j)+theta(j0+i)	!Ymax for component #1, set #iset
		   enddo
		enddo
	   endif
c
	   if(ifitmode.eq.3) goto 35		!display -after separate fits for 4,5
	endif
c For ifitmode=4,5 simultaneous fit now finished, and calc curves (constrained
c to be parallel) calculated.  Before plotting do also the separate fits to
c each curve, and calc the curves for these too (also get reduction in SSD to
c do parallelism test)
46	continue		!return here after likelihood intervals for 4,5 calcs
	if(ifitmode.eq.4.or.ifitmode.eq.5) then
	   SSDpar=smin		!save SSD for parallel (simultaneous) fit
	   kfitpar=kmax-nfix	!save number of param fitted for ditto
	   IFLAG=2		!to indicate last section of ifitmode=4,5 calcs
	   call DCASK(
     &   'Calculate likelihood intervals for simultaneous fit','n',ans)
	   if(ans.eq.'Y') goto 45
         print 38
c         if(pon()) write(7,38)
c         if(discprt) write(8,38)
38	   format(/,' SEPARATE FITS TO EACH SET INDIVIDUALLY')
	   isepfit=1	!in case nmod is not 1 or -1
	   if(iabs(nmod).eq.1) then	!Langmuir (1) or Hill eqn (-1)
48		print 47,isepfit
47		format(' For separate fits:',/,
     &	' (1) Allow all parameters to vary',/,
     &	' (2) Fix all maxima at value for simultaneous fit',/,
     &	' (3) Constrain maximum to be same for all curves',/,
     &	' (4) Omit the separate fits',/,
     &	' Option number [',i2,'] = ')
		call INPUTi(isepfit)
c		read 104,i
		if(isepfit.lt.1.or.isepfit.gt.4) goto 48
		if(isepfit.eq.4) goto 351	!display simultaneous fit only
c            if(pon()) write(7,38)		!sep fit
            if(discprt) write(8,38)
		if(isepfit.eq.1) then
c      	   if(pon()) write(7,471)
	         if(discprt) write(8,471)
471  		   format('  (All parameters allowed to vary)')
		else if(isepfit.eq.2) then
c      	   if(pon()) write(7,473)
	         if(discprt) write(8,473)
473  		   format('  (maximum fixed at value for parallel fit)')
		else if(isepfit.eq.3) then
c      	   if(pon()) write(7,472)
	         if(discprt) write(8,472)
472  		   format('  (maximum constrained to be same for all)')
		endif
	   endif
c	Find first set to be fitted = #jset
	   j=0
40	   j=j+1
	   if(jmiss(j).eq.1) goto 40
	   jset=j
c Set initial guesses (for Hill/Langmuir etc at least) in thetsav(i,j) so need
c not ask for init guesses in CVSIMP
	   if(nmod.le.4.or.nmod.eq.99) then
		fixmax=iabs(nmod).eq.1.and.jfix(2).eq.1   !Ymax fixed for Lang/Hill
		fix0=nmod.eq.-1.and.jfix(3).eq.1		!Y(0) fixed for Hill
		iset=0
		do 42 j=1,nset
	 	if(jmiss(j).eq.1) goto 42
	      iset=iset+1
		if(ifitmode.eq.4.and.iset.gt.1) then
		   ratio=theta(kmax1+iset-1)
		else if(ifitmode.eq.5) then
		   akb=theta(kmax)	!=KB (see GETEQN)
		   xB=setx(j)		!antag conc
		   ratio=(1.0 + xB/aKB)	!mult by dose ratio
		endif
c now multiply all KA values, for iset>1 by 'ratio' to get init guesses
c for each set
		if(iset.eq.1) then		!should have j=jset for iset=1
		   do 41 i=1,kmax1
41		    thetsav(i,jset)=thetsav(i,11)      !values from simultaneous fit
		else
		   do 43 i=1,kmax1
43		    thetsav(i,j)=thetsav(i,11)      !values from simultaneous fit
		   if(nmod.le.4) then
		     thetsav(1,j)=thetsav(1,j)*ratio 			!theta(1)=K1
		     if(ncomp.eq.2) thetsav(3,j)=thetsav(3,j)*ratio	!=K2
		   else if(nmod.eq.99) then		!polynomial
		     do 44 i=2,ip1+1
44			thetsav(i,j)=thetsav(i,j)/(ratio**(i-1))
		   endif
c		   if(debug()) then
c		     print 431,j,(thetsav(i,j),i=1,kmax)
c431		     format(' Set #',i3,' thetsav = ',8g13.6)
c		   endif
		endif
42		continue	!end of j=1,nset loop
	   endif
	   if(isepfit.ne.3) then
	      ifitmode=2		!to do sep fits with original code
	      kmax=kmax1 		!reset to value without extra params (nmod the same)
	      goto 31		!repeat fitting with ifitmode=2, same nmod
	   else if(isepfit.eq.3) then
c		Redefine model and ifitmode, and set initial guesses in thetsav(i,1)
		call NEWEQCV(nmod,kmax,kmax1,titlep,nsfit,juse,ifitmode,
     &	  thetsav,setx,njset,ndth)
		if(debug()) then
		  j=1
		  print 431,j,(thetsav(i,j),i=1,kmax)
431		  format(' Set #',i3,' thetsav = ',20g13.6)
		endif
	      goto 31		!repeat fitting with ifitmode=3, new nmod
	   endif
	endif
c
35	continue
c Now all fitting done -display graph (if relevant, do dev from parallelism
c first)
c Print test of deviation from parallelism following ifitmode=4 or 5
	if(iflag.eq.2) then
	   ndferr=nfittot-kfitsep
	   ndfdev=kfitsep-kfitpar
	   if(ndferr.le.0.or.ndfdev.le.0) then
		call BELL(2)
      	print 53,ndferr,ndfdev,nsfit,nfittot,kfitsep,kfitpar
	      if(pon()) write(7,53) ndferr,ndfdev,nsfit,nfittot,
     &	kfitsep,kfitpar
	      if(discprt) write(8,53) ndferr,ndfdev,nsfit,nfittot,
     &	kfitsep,kfitpar
53		format(' PROBLEM!',/,
     &    ' Degrees of freedom for error = ',i8,/,
     &    ' Degrees of freedom for dev from parallelism = ',i8,/,
     &    ' Number of sets fitted = ',i8,/,
     &    ' Number of observations fitted = ',i8,/,
     &   ' No of parameters fitted for separate, parallel fits = ',2i8)
		goto 351	!skip dev from parallelism test
	   endif
	   aMSerr=SSDsep/float(ndferr)      !error mean square
	   SSDdev=SSDpar-SSDsep			!=improvement for sep fit
	   aMSdev=SSDdev/float(ndfdev)
	   VR=aMSdev/aMSerr			!variance ratio
c	   if(debug()) print 700,nfit,kfitsep,SSDsep,kfitpar,SSDpar
	   print 700,nfit,kfitsep,SSDsep,kfitpar,SSDpar
700	   format(' nfit,kfitsep,SSDsep,kfitpar,SSDpar = ',/,
     &    2i8,g13.6,i8,g13.6)
	   if(vr.ge.0.) then
	      call FDIST(VR,ndfdev,ndferr,P)
	   else
		call BELL(2)
	      print 54,vr,ssdsep,ssdpar
	      if(pon()) write(7,54) vr,ssdsep,ssdpar
	      if(discprt) write(8,54) vr,ssdsep,ssdpar
54		format(' PROBLEM!',/,
     &	' Variance ratio = ',g13.6,/,
     &	' SSD for separate and parallel fits = ',2g13.6)
		goto 351	!skip dev from parallelism test
	   endif
         print 52,aMSerr,ndferr,aMSdev,ndfdev,vr,P
         if(pon()) write(7,52) aMSerr,ndferr,aMSdev,ndfdev,vr,P
         if(discprt) write(8,52) aMSerr,ndferr,aMSdev,ndfdev,vr,P
52	   format(/,' TEST for DEVIATION FROM PARALLELISM',/,
     &  ' Mean square error (from separate fit) = ',g13.6,
     &  ' (',i5,' d.f.)',/,
     &  ' Mean square for diff between parallel and sep fits = ',g13.6,
     &  ' (',i5,' d.f.)',/,
     &  ' Variance ratio = ',g13.6,'   P = ',g13.6)
	endif
351	print 4
4	format(' Type <any key> for display')
	call ANYKEY
	goto 3062	!for display/plot
c
320	continue
	if(iflag.eq.2) goto 99	!already done for ifitmode=4,5
	call DCASK('Calculate likelihood intervals','N',ans)
	if(ans.eq.'N') goto 99
45	continue
	elmax=-smin/2.
	if(iw.eq.1.or.iw.eq.5) elmax=elmax/(sres**2)	!use Sres as error var
	kmax2=kmax	!for likblk
	call LIKINT1(ELFUNCV,noerr,sdev,nmod,constr,
     & theta,jfix,titlep,ndth,Xobs,yobs,w,nj,juse,setx,niobs,njset)
	if(iflag.eq.2) goto 46	!complete ifitmode=4,5 calcs
	goto 320
c
99	continue
	if(VIDEOTYP().ne.3) call VIDEOMOD(3) !utility lib- back to alpha mode
c Write values to .INI
      OPEN(unit=19,file='CVFIT.INI',status='UNKNOWN',
     & access='DIRECT',form='UNFORMATTED',recl=256)
	 write(19,rec=1) ndev,ifile,ifitsav,ilog,idiskq,titlef,infil,
     &	niobs,njset,ascfil
	CLOSE(unit=19)
c
c=	call SAVEDP
	fitted=.false.
	CLOSE(unit=12)
	fopen=.false.
c Rather than call SAVEDP, put its content here, but ask 'another run?'
c while disc file is closed, so when this question is on screen, it
c will be possible to alt-tab to look at output file eg in notepad
	if(discprt) then
	   CLOSE(unit=8)
	   print 3
3	   format(' Another run [Y] ? ')
	   read 101,ans
         OPEN(unit=8,file=FILNAM,status='UNKNOWN',
     &   access='APPEND',form='FORMATTED',carriage control='FORTRAN')
	endif
	idest=1			!not first time in CVIN
	if(UC(ans).eq.'N') goto 999
	fitted=.false.
	dcurve=.false.
	fixmax=.false.
	fix0=.false.
	ncurvc=0
	itit=0
	iflag=0
	isepfit=1
	nfix=0		!no param fixed
	do 62 i=1,20
62	jfix(i)=0
	goto 2		!another run
c
999	continue
 	call VIDEOMOD(3)	 !utility lib- back to alpha mode
	call ENDPRINT
c
	END


