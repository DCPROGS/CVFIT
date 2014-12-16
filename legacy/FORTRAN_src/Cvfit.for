	PROGRAM CVFIT

c
c GENERAL CURVE FITTING PROGRAM
c Lahey 90/Gino version of CVFIT
c
c Modifs 08/24/97 07:14am Major changes to weighting and error options
c (esp to allow use of 'within group' error when there are replicates)
c
c (1) Option to read earlier data from PLOTQ rather than CVDAT?
c
c
c Things still to be done in the new revision
c  (1) Fix initial guesses for straight line/polynomial fits (and other eqns?)
c  (2) Incorp in nmod=27,28
c	(a) fit of binding inhibition curves -param=KB
c	(b) linear component
c
c 01/11/97 06:49am Change method for fixing parameters.  Up to now
c jfix has been zeroed each time CVSIMP entered.  CVSIMP is entered
c more than once only
c  (a) for ifitmode=2 only (entered nsfit times=# of sets fitted)
c  (b) for ifitmode=2 when it FOLLOWS ifitmode=4,5 (entered once for
c  	ifitmode=4 (or 5) and nsfit times for subsequent sep fits
c  (c) for ifitmode=2 when it PRECEDES some other fit mode (for
c	purpose of normalising curves on basis of sep fits)
c  FIXED PARAMS.
c	(1) IFITMODE=2 first:
c	(a) If guesses made and accepted then ask which
c	params to be fixed in getguess (rather than before simplexv), and
c	ask (once) if same eqn to be fitted to all sets -if so then don't
c	reset jfix each time cvsimp is entered (and don't ask again re
c	which eq to be fitted, ncomp etc)
c	(b) If no guesses available, or not accepted, then ask about fixed
c	params in usual place for first set only, if same eqn fitted to all
c	(c) If NOT same eqn for all curves, then must ask everything for
c	every set (guesses no use in this case)
c	(d) If ifitmode=2 is followed by another fit (eg mode=3,4,5,6) then
c	reset jfix when ifitmode=2 completed, and enquire re fixing for
c	the next fitmode -if same eqn is being fitted can give sensible
c	defaults on basis of what was fixed for the sep fits.
c
c	(2) IFITMODE=4,5 first
c	(a) If guesses made and accepted then ask which
c	params to be fixed in getguess (rather than before simplexv).
c	(b) If no guesses available, or not accepted, then ask about fixed
c	params in usual place for first set only, if same eqn fitted to all
c	(c) When moving on to fitmode=2 for sep fits, enquire again re
c	which to be fixed, with sensible defaults on basis of what was
c	fixed during the simultaneous fit (dep on which isepfit chosen)
c	(d) in this case sep fits will always be with same eqn as for
c	the simultaneous fit so needn't ask about this, and don't reset
c	jfix on entering cvsimp (it will have been set during transition
c	from mode=4/5 to mode=2)
c
c
c
c 01/08/97 12:24pm
c (1) Added option after sep fits (ifitmode=2) with Hill/Langmuir to
c	normalise wrt fitted max and refit with another mode (uses iflag=3,4)
c     (prob needs more work if ncomp>1)
c (2) Added ifitmode=6 to pool several sets into set 1 (sets
c  norm1=true so means plotted also in CVDISP -no option yet to use SD
c  calculated here for weighting
c===================
c
c 12/26/96 07:21am Lot of changes in initial guesses (for Hill/Langmuir/
c power function).  Initial guesses for each set separately found in
c thetgues(i,j), and, for ifitmode.ge.3 (simultaneous fits) these used
c to calculate initial guess for the sim fit, in new subroutine, getguess.
c After simultaneous fit, results kept in thetsav(i), and used to get
c initial guesses for subsequent sep fits for ifitmode=4,5
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
c
c arrays to 'superimpose last fit'
	ALLOCATABLE::xcalsav,ycalsav
	real*4 xcalsav(:,:),ycalsav(:,:)
	ALLOCATABLE::ncalsav,icrvcsav,ilinesav
	integer ncalsav(:),icrvcsav(:),ilinesav(:)
c
c	real YOBS(100,10),XOBS(100,10),W(100,10)
c	dimension XVAL(2048,10),YVAL(2048,10)
c	real xcal(2048,20),ycal(2048,20)
	ALLOCATABLE:: titled,jmiss,juse,iuse,nj,setx
	character*60 TITLED(:)
c	real weight(100,10)
c	real weight(kwi,kwj)
	integer jmiss(:),juse(:),iuse(:),nj(:)
	integer icuse(100) 	!see geteqn, cvsimp
	real setx(:)
c	character*60 TITLED(10),titlef	!data set names,file name
c	real setx(10)
	real*4 ymaxj(100)
	character*60 titlef		!file name
	character string*50,text*27
	character ascfil*33
	character*60 titw(5)		!names of weighting methods
c	real VAR(20)
	real Sdev(100),sdevsav(100)
c=	real*4 Sdevr(100)		!for alternative Sdev from residuals
	real*4 theta(100)
	character*20 TITLEP(100)		!names of params for fit
	real yval1(1024)
	real XAXIS(510),xsav(5,4)	!for histo display
	integer JFIX(100),jomit(20),jfixsav(100),iparsav(10)
	character*1 ans,UC
c	character*1 getch		!for debug
c	integer*2 ktype		!for debug
	integer*4 videotyp
	character ndev*2,infil*33
	ALLOCATABLE:: thetsav,thetgues
	real*4 thetsav(:),thetgues(:,:)
	ALLOCATABLE:: cormat
	real*4 cormat(:,:)
	character defolt*30,cans*30		!to hold text & result of QDIALOG
	character*20 title 			!for popmenu
	common/popvalm/nxlo,nylo,nyhi,ictm,ibkm,icfm,icupm	!values for popmenu calls
	character*55 strings(8)			!for popmenu
	character*55 helps(8)			!for popmenu
c Allocate other arrays for VPLOT
	ALLOCATABLE:: ndat,icurvd,icurvw,isym,ijoin,symsiz
	integer ndat(:),icurvd(:),icurvw(:),isym(:),ijoin(:)
	ALLOCATABLE:: ncal,icurvc,iline
	integer ncal(:),icurvc(:),iline(:)
	real symsiz(:)				!for data
c	dimension ndat(10),icurvd(10),isym(10),ijoin(10)	!for data
c	real symsiz(10)				!for data
c	dimension ncal(20),icurvc(20),iline(20)		!for calc curve
	integer ncal2(10),icurvc2(10),iline2(10)		!for histo disp
	LOGICAL CONSTR,FLINE,SETVAR,defolt1,fitted,dcurve,logt,badset
	logical xgt0,noerr,fixmax,fix0,readini,alpha,allocated
	logical nodata,hillang,autosep,neweq,guesdone,norm1,fixset,sameq
	logical noguess(100),doframe,screen,errflag,apfile
	logical logx,logy,logity
	character*40 titlex,titley
	character*64 title1
	logical discprt,pon,slock,hdisp,fopen,printed,present
	logical debug,caplock,plotrue
	character*11 cdate,ctime
	character cnum*11,cnum1*11
c For mouse
	logical mouse_on
c	logical mous_set,eflag
	common/mousval/mouse_on,nbutton
c Also allocate arrays for GETGROUP  -needed for call to CVIN2(can deallocate
c again if not needed)
	ALLOCATABLE::xgrp,ybar,sdy,sdm,ngrp,nxval
	real*4 xgrp(:,:),ybar(:,:),sdy(:,:),sdm(:,:)
	integer ngrp(:,:),nxval(:)
c
	character mtitle*40,filnam*32,prtport*4,path*40	!for WINPRINT
	common/dpp/filnam,prtport,ndisc,icol,mtitle !for WINPRINT,ENDPRINT,DISCNUM
	common/dp/discprt
c= xobs, yobs no longer in common
c=	common/data/Xobs,Yobs
	COMMON/LIKBLK/kmax2,kfit,elmax,em,delmin,ndisp,npar,iw,sres
	common/logblk/debug,defolt1
	COMMON/BLOCK1/constr,nset,nfit,nsfit,
     & Xv,alfa,kmax,ncomp,nmod,fline,nomit,jomit,jset,ifitmode
c=	COMMON/BLOCK2/ castar,setx(10),X1,X2,iequiv,ip1,ip2
c====setx no longer in common/block2/
	COMMON/BLOCK2/ castar,X1,X2,iequiv,ip1,ip2
c====w(i,j) no longer in common!
c=	COMMON/SD/w		  !for CVDISP2,VPLOT2,SSDCV,YCALCV
	common/potrat/jfirst,iset,kmax1			!for pot ratios
	common/max/ymaxj		!for CVDISP and YCALCV
	logical logyfit,norm
	logical cluster,student
	COMMON/BLOCK3/logyfit,norm,xnorm,iostat
c (kmax=actual no of param to be estimated; kmax1=number before extra param
c  added for ifitmode=4,5)
	common/graf/alpha		!for simplexv
	common/popvalt/lxlo,lylo,lyhi,ictt,ibkt,icft,icht	!values for poptable calls
	common/readqplot/plotrue,istrec,iptype,kwi,kwj,ndv1,ndc1
	COMMON/JLOGOS/t1c(10),t2c(10),t1v(10),t2v(10),xoff1,y1v,y2v,
     & y1c,y2c
	common/user/student,cluster,iwindows

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
	cluster=.false.
	student=.false.
	call MYPATH(path)
	if(path(1:2).eq.'O:'.or.path(1:2).eq.'N:') cluster=.true.
	if(cluster) then
	   prtport='lpt1'
	   ndev='C:'
	   ndisc=26
	   icol=1
	   student=.true.
	   iver1=1
	   apfile=.false.
	   discprt=.true.
	   mtitle='Dell'
	   iwindows=2
	   filnam='o:\cvfit\*.prt'
995	   call TITENT0(
     &   'Name for disk print file:Please enter your name in the place
     &of * ',filnam,32,.false.)
	    nr=nblank(filnam)
	   if(filnam(1:9).ne.'o:\cvfit\') then
	      print*,'TRY AGAIN: something like  o:\cvfit\yourname.prt '
            goto 995
	   endif
	   if(filnam(10:11).eq.'*.') then
	      print*,'TRY AGAIN: something like  o:\cvfit\yourname.prt '
            goto 995
	   endif
         if(discprt) OPEN(unit=8,file=FILNAM,status='UNKNOWN',
     &   access='APPEND',form='FORMATTED',carriage control='FORTRAN')
	   if(discprt.and.(.not.apfile)) REWIND(unit=8)
	endif
c
c Use F90 routine to prevent underflow crashes??
	errflag=.true.
c	eflag=.true.
	call UNDFL(errflag)
c	call DVCHK(eflag)
c
c	iver=1002		!new version with file title
c	iver=1003		!new version with alloc arrays and new CVDAT org.
c NB iver is initially the version number for the CVDAT being read in
c (may be 1001, 1002, 1003,..).  Writing is done in CVDISK only in the
c format for iver=1003 so set iver=1003 there (and whole CVDAT must be the
c SAME version!)
	ndth=100		!dimension of theta(),jfix(),titlep(),sdev()
	ncurvsav=-1		!until a fit has been done
	call SETMOUSE()		!define values in common\mousval\ (in IVLIB)
	call GINO
	call vga
	call mode(3)
c In scan, colours are set in all subroutines via
c common/dialval1/ictd,icfd,icbd  !colours for dialog box #1 (text,frame,back)
c but this not used in CVFIT at present, and ict=11 used only in main prog.
	ict=11	!text colour in dialog box 1
c For poptable define
	lxlo=100
	lylo=-1
	lyhi=400
	ictt=1		!text dark blue
	ibkt=7		!light grey background
	ibkt=15		!white background
	icft=3		!frame light blue
	icht=12		!header/footer red
c Values for POPMENU calls
c	common/popvalm/nxlo,nylo,nyhi,ictm,ibkm,icfm,icupm	!values for popmenu calls
	nxlo=100
	nylo=-1
	nyhi=479
c	ictm=15		!text white
c	ibkm=8		!background dark grey
c	icfm=2		!frame/title green
c	icupm=12		!upper case red
c alternative colour scheme (as in VPLOT)
	ictm=15		!text white
	ibkm=1		!background daek blue
	icfm=11		!frame/title light blue
	icupm=12		!upper case red
c
	if(cluster) goto 1122
      filnam='CVFIT.PRT'
	call WINPRINT	!print-file control
1122	continue
      OPEN(unit=7,file=prtport,iostat=nerr)		!open printer
	print 1011
c	if(pon()) write(7,1011)
	if(discprt) write(8,1011)
1011	format(' CVFIT PROGRAM',/)
	call GETSPEC('cvfit.exe')	!print date/time of the program
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
 	if(cluster) then
 	   ndev=' '
 	   ifile1=2
 	   ifitmode=0
 	   ilog=1
 	   idiskq=-1
 	   titlef=' '
 	   infil=' '
 	   niobs=60
 	   njset=35
 	   ascfil=' '
 	   nsfit=1
	   goto 12
 	endif
	INQUIRE(file='CVFIT.INI',exist=present,flen=len)
	if(present.and.len.gt.0) then
	   if(len.eq.256) then
          OPEN(unit=19,file='CVFIT.INI',status='UNKNOWN',
     &    access='DIRECT',form='UNFORMATTED',recl=256)
	    read(19,rec=1) ndev,ifile1,ifitmode,ilog,idiskq,titlef,infil,
     &	niobs,njset,ascfil,nsfit
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
	ALLOCATE(titled(njset),jmiss(njset),juse(njset),iuse(njset),
     & nj(njset),setx(njset))
	ALLOCATE(ndat(ndimd),icurvd(ndimd),icurvw(ndimd),isym(ndimd),
     & ijoin(ndimd),symsiz(ndimd))
c=	ALLOCATE(xcal(ndc1,ndimc),ycal(ndc1,ndimc))	!later
	ALLOCATE(ncal(ndimc),icurvc(ndimc),iline(ndimc))
	ALLOCATE(ncalsav(ndimc),icrvcsav(ndimc),ilinesav(ndimc))
c Also allocate arrays for GETGROUP  -needed for call to CVIN2(can deallocate
c again if not needed)
	ALLOCATE(xgrp(niobs,njset),ybar(niobs,njset),
     &  sdy(niobs,njset),sdm(niobs,njset),ngrp(niobs,njset))
	ALLOCATE(nxval(njset))
c
	ndev='  '		!before disk # defined
	fopen=.false.	!data file not open yet
	ilabel=-1	!internal axis labels unless otherwise specified
      print 39
39	format(' (1) Curve fitting',/,' (2) Histogram display',/,
     & ' Option [1] = ')
	call INPUTi(i)
	idiskq=-1	!until plot queue disc defined
	guesdone=.false.
	norm1=.false.	!for norm after sep fits (see cvin2, cvdisp)
	fixset=.false.	!params to be fixed not set yet
	sameq=.false.	!fit same eq to all sep fits that follow ifitmode=4,5
2	continue
	if(idest.ne.2) idest=0
	iflag=0	!for GETEQN
	lflag=0	!for log scale question
	iesc=0	!for DISCNUM abort
	isepfit=1	!in case nmod is not 1 or -1
	ifitlast=0
	fitted=.false.
	dcurve=.false.
	fixmax=.false.
	fix0=.false.
	autosep=.false.	!separate fits done with defaults -no questions
	do i=1,100
	   noguess(i)=.false.
	enddo
	neweq=.false.
	noerr=.false.
	ncomp=1	!(needs to be defined in cvsimp, even if not used)
	iparsav(1)=-1	!set to ifitmode after ifitmode=4,5
	if(i.le.0) i=1
	hdisp=i.eq.2
c	nfirst=1
	nomit=0
      plotrue=.false.	!set true in CVIN2 if data read from plot queue
c	istrec=-1000
c	 iptype=-1111
	iver=1100
	call CVIN2(IW,setvar,aval,bval,defolt1,ndev,titw,nmiss,
     & titlef,titled,titlex,titley,ilabel,yobs,xobs,w,idest,norm1,
     & fopen,ncalc,ncurvc,icurvc,ncal,iline,hdisp,yval1,nyval,printed,
     & readini,ifile,ilog,idiskq,ifitsav,iesc,ndisc,infil,ifile1,
     & niobs,njset,ndimd,ndimc,nj,jmiss,juse,setx,iver,ascfil,nodata,
     & guesdone,nxval,xgrp,ybar,sdy,sdm,ngrp,SSwg,ndf)
c Now define iuse()= curve number that corresponds to sets specified
c in juse(), which has been defined in CVIN2
	do i=1,njset
	   iuse(i)=0	!zero unused elements
	enddo
	do m=1,nsfit
	   j=juse(m)
	   iuse(j)=m
	enddo
c
	if(plotrue) then
	   if(.not.allocated(xval)) then
	      ALLOCATE(xval(ndv1,ndimd),yval(ndv1,ndimd))
	      ALLOCATE(xcal(ndc1,ndimc),ycal(ndc1,ndimc))
		ALLOCATE(xcalsav(ndc1,ndimc),ycalsav(ndc1,ndimc))
	   endif
	   screen=.true.
	   idev=0
	   call RDVPLQ(istrec,xval,yval,xcal,ycal,ndimd,ndimc,ncurvd,
     & ndat,isym,ijoin,ncurvc,ncal,iline,syms,xmin,xmax,ymin,ymax,
     & xcross,ycross,xtic,ytic,ntx,nty,itx,ity,xlo,xhi,ylo,yhi,
     & itit,title1,ifont,ilog,iscal,doframe,titlex,titley,ilabel,
     & inumx,inumy,sval,theta,ifitype,ncomp,isdev,w,y0,yinf,iptype,
     & ncjump,t1c,t2c,nvjump,t1v,t2v,xoff1,y1v,y2v,y1c,y2c,ivplot,
     & interp,screen,colplotter,
     & itrace,ntrace,ytsep,ndv1,ndc1,
     & kwi,kwj,icurvw,kmax,iver,idev)
	   if(ncurvc.gt.0) fitted=.true.
	   if(ncurvd.gt.0) then
		do i=1,ncurvd
		   icurvd(i)=i
		enddo
	   endif
	   if(ncurvc.gt.0) then
		do i=1,ncurvc
		   icurvc(i)=i
		enddo
	   endif
c	   cbig=2.5
c	   csize=cbig	!from queue
c	   thbig=1.0		!in common for line thickness
c	   isetcol=0		!for VPLOT
c
	   logx=ilog.eq.1.or.ilog.eq.3.or.ilog.eq.4.or.ilog.eq.6
	   logy=ilog.eq.2.or.ilog.eq.3
	   logity=ilog.eq.4
	   if(.not.logx) inumx=1	!see notes in LAXES2
	   if(.not.logy) inumy=1	!see notes in LAXES2
	   goto 3062
	endif
c
C For 'nodata' case; need to get equation/get params/
c then calc values and display as for dcurve -skip first display of data!
c
	if(iesc.eq.98) then	!nset too big -reallocate arrays
	   njset=nset+5		!new default size
c=	   DEALLOCATE(xval,yval,xcal,ycal)
	   DEALLOCATE(xobs,yobs,w)
	   DEALLOCATE(titled,jmiss,juse,iuse,nj,setx)
	   DEALLOCATE(ndat,icurvd,icurvw,isym,ijoin,symsiz)
	   DEALLOCATE(ncal,icurvc,iline)
	   goto 12
	endif
	if(iesc.eq.99) goto 99	!DISCNUM abort
	if(.not.allocated(xval)) then
	   ALLOCATE(xval(ndv1,ndimd),yval(ndv1,ndimd))
	   ALLOCATE(xcal(ndc1,ndimc),ycal(ndc1,ndimc))
	   ALLOCATE(xcalsav(ndc1,ndimc),ycalsav(ndc1,ndimc))
	endif
c
c  Now nset is known can allocate thetsav() and thetgues()
	if(.not.allocated(thetsav)) then
	   ALLOCATE(thetsav(100),thetgues(100,nset))
	endif
c
c Now initial display of data before fitting
	if(nodata) then
	   goto 33
	endif
	fitted=.false.
	dcurve=.false.
3062	continue		!return here after fitted for final display

	if(.not.hdisp) then
	  nset1=nset	!copy for call
	  call CVDISP(Xobs,Yobs,w,ISYM,ILOG,ISCAL,IJOIN,iline,
     &  XMIN,XMAX,YMIN,YMAX,XTIC,YTIC,xlo,xhi,ylo,yhi,y0,yinf,
     &  xcross,ycross,titlex,titley,title1,itit,ilabel,inumx,inumy,
     &  xval,yval,ncurvd,ndat,icurvd,ncurvc,ncal,icurvc,Xcal,Ycal,
     &  symsiz,idiskq,ntx,nty,itx,ity,theta,fitted,doframe,
     &  dcurve,nodata,iw,idest,readini,iflag,ndv1,ndimd,ndc1,ndimc,
     &  jmiss,nj,niobs,njset,icurvw,juse,ndth,nset1,thetgues,norm1,
     &  plotrue,iver,ibk,
     &  ncurvsav,xcalsav,ycalsav,ncalsav,icrvcsav,ilinesav)
	  plotrue=.false.
	  logx=ilog.eq.1.or.ilog.eq.3.or.ilog.eq.4.or.ilog.eq.6	!for DEFYCAL
	  logy=ilog.eq.2.or.ilog.eq.3
	  logity=ilog.eq.4
	  if(idest.eq.1) goto 320	!after display with fit, calc likelihood int?
c Case of refit after normalising (idest=2).  Guesses are already in thetgues(i,j)
c for each of the sets (j) fitted separately before normalisation, but if
c refitted with ifitmode=6,pooled data are put into set 1 (for which a guess
c may not exist, and if it does, may not be appropriate) so in this case
c set guesdone=false in refit section of cvin2, so another guess made in cvsimp
	  if(idest.eq.2) then
		guesdone=.true.	!values for sep fits are in thetgues(i,j)
		autosep=.false.
		mset=1
		goto 2	!refit norm data after ifitmode=2
	  endif
	  if(idest.eq.99) goto 99	!another run?
c       NB if idest=0 (fit next...) then still in graphics mode
	else
	   ncurvc2=ncurvc
	   do 411 i=1,10
	    ncal2(i)=ncal(i)
	    icurvc2(i)=icurvc(i)
411	    iline2(i)=iline(i)
	   call CHDISP(Xobs,yval1,nyval,ILOG,ISCAL,IJOIN,iline2,
     &   XMIN,XMAX,YMIN,YMAX,XTIC,ytic,xlo,xhi,ylo,yhi,
     &   xcross,ycross,titlex,titley,title1,itit,ilabel,inumx,inumy,
     &   ncurvd,ndat,icurvd,ncurvc2,ncal2,icurvc2,idiskq,
     &   ntx,nty,itx,ity,theta,fitted,doframe,idest,
     &   xaxis,xsav,nbin,nbw,xwbase,nbdec,logt,ibk)
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
		string='                                                  '
		do j=1,nset
		   if(jmiss(j).ne.1) then
			string(n:n+1)=char(j+48)//','
			n=n+2
		   endif
		enddo
		n=nblank(string)
		string(n:n)=' '		!remove final comma
	      call WDIALOG(1,'SIMULTANEOUS FIT TO SETS '//string(1:n),ict)
	      call WDIALOG(1,'  (yellow points are set 1)',ict)
	   endif
         if(pon()) write(7,321)
         if(discprt) write(8,321)
321	   format(/,' SIMULTANEOUS FIT TO THE FOLLOWING SETS')
	   if(iver.ge.1002) then
	      if(pon()) write(7,322) titlef
	      if(discprt) write(8,322) titlef
322		format(/,' File: ',a60)
	   endif
	   do 32 j=1,nset
	    if(jmiss(j).eq.1) goto 32
	    text='FIT TO DATA SET # '
	    if(pon()) write(7,1092) text,j,titled(j)
          if(discprt) write(8,1092) text,j,titled(j)
32	   continue
	endif
	do i=1,100
	   icuse(i)=0
	enddo
	SSDsep=0.0		!initialise for case where summed for each set
	nfittot=0		!ditto
c
31	continue		!return to fit next set
	alpha=VIDEOTYP().eq.3		!alpha mode
	if(ifitmode.eq.1.or.ifitmode.eq.2) then
	   if(iver.ge.1002) then
	      if(pon()) write(7,322) titlef
	      if(discprt) write(8,322) titlef
c322		format(/,' File: ',a60)
	   endif
	  text='FIT TO DATA SET # '
	  if(norm1) text='FIT TO POOLED DATA IN SET #'
	  if(dcurve) text='CURVE + DATA SET #'
	  if(nodata) text='CALCULATED CURVE #'
	  if(alpha) then
	     print 322,titlef
	     print 1092,text,jset,titled(jset)
	  else
	     call INTCONV(jset,cnum)
	     string=titled(jset)(1:50)
	     call WDIALOG(1,'File: '//string,ict)
	     call WDIALOG(1,CHARNB(text)//charnb(cnum),ict)
	  endif
	  if(pon()) write(7,1092) text,jset,titled(jset)
        if(discprt) write(8,1092) text,jset,titled(jset)
1092	  FORMAT(/,' ',a27,i3,/,3x,a60)
	endif
13	if(iflag.ne.2) then
c Set default for 'y increases with x' on basis of first set (iequiv in common)
	   j1=juse(1)
c Set default for 'y increases with x' on basis of first set (iequiv in common)
c NB care needed if data are not entered in ascending order of x
         j1=juse(1)
         xmax2=-1.e37
         xmin2=1.e37
         do i=1,nj(j1)
             if(Xobs(i,j1).lt.xmin2) then
                xmin2=Xobs(i,j1)
                imin=i
             else if(Xobs(i,j1).gt.xmax2) then
                xmax2=Xobs(i,j1)
                imax=i
             endif
         enddo
         if(Yobs(imin,j1).lt.Yobs(imax,j1)) then
c=	   if(Yobs(1,j1).lt.Yobs(nj(j1),j1)) then
		iequiv=0       !increasing
	   else
		iequiv=1		!decreasing
	   endif
	   call GETEQN(IW,titw,titlep,xgt0,iflag,sameq,
     & 	 jmiss,juse,iuse,icuse,njset,ndth,imod,nodata,iesc)
c for ifitmode=2 need to reset icuse() for each separate fit
c to value for curve # that corresponds with current set
	   if(ifitmode.eq.2) then
		do m=1,kmax
		   if(icuse(m).ne.0) icuse(m)=mset
		enddo
	   endif
	endif
c If no data, get equation but skip fit
c For ifitmode=4,5: kmax redefined in GETEQN as is kmax1=orig no of param
	if(iesc.eq.99) goto 99	!fitting abort
	nset1=nset	!copy for call
	call CVSIMP(kfit,jfix,nfix,titlep,fitted,dcurve,
     & nodata,theta,thetsav,nev,sres,smin,iflag,fixset,
     & Xobs,yobs,w,nj,juse,icuse,setx,niobs,njset,ndth,
     & badset,nset1,thetgues,xcal,ycal,ndc1,ndimc,ncalc,logx,logy,
     & logity,y0,yinf,xgt0,xmin,xmax,ymin,ymax,autosep,
     & iparsav,jfixsav,sminsav,neweq,guesdone,noguess,mset,ibk,
     & nxval,xgrp,ybar,sdm,ngrp,SSwg,ndf)
	if(neweq) then	!try another equation
	   neweq=.false.
	   goto 13
	endif
	if(badset) then
	   if(ifitmode.eq.1) then
		goto 99	!only 1 set and it is bad
	   else if(ifitmode.eq.2) then	!skip to next set
		mset=mset+1
		if(mset.gt.nsfit) goto 35	!display
		j=juse(mset)		!next set # to be fitted
		jset=j
		if(iflag.eq.0) iflag=1		!to tell GETEQN this is not 1st set
c Set initial guesses (not always defined!)
		do i=1,kmax
		   theta(i)=thetgues(i,jset)
	 	enddo
		goto 31		!fit next set
	   endif
	endif
c
c CALC VARIANCE ETC FOR THE Kfit NON-FIXED PARAMETERS?
c For separate fits following ifitmode=4,5, ask first time only
	ndf=nfit-kfit		!degrees of freedom for error
	if(dcurve.or.nodata.or.ndf.eq.0) then
	   noerr=.true.
	else if(.not.(autosep.and.mset.gt.1)) then
	   ans='Y'
	   if(noerr.or.iflag.eq.2) ans='N'
	   if(alpha) then
		call DCASK('Calculate approximate errors',ans,ans)
	   else
		call DEFOLTa(ans,defolt)
		call QDIALOG(1,'Calculate approximate errors',
     &	     defolt,11,cans)
	      call GETINPa(cans,ans)
	   endif
	   noerr=ans.eq.'N'
	endif
c
c NB must allocate cormat for cvout1, whether used or not
	ALLOCATE(cormat(kmax,kmax))	!for correlation matrix
	if(noerr) goto 201	!no errors
c	if(idest.eq.201) goto 201	!no errors
c
	kfit=kmax-nfix		!redefine (in case kmax redef for ifitmode=4,5)
	if(debug()) then
	   print 600,kmax,nfix,kfit
600	   format(' kmax,nfix,kfit = ',3i8,/)
	endif
c Define serr for scaling finfo etc in hesmat3
c If weights used then calc also Sresiduals using weights=1, so errors
c can be printed ways (use SSCDVr, which is version of SSDCV without weights
	if(iw.eq.2.or.iw.eq.3.or.iw.eq.4) then
	   serr=1.0		!error SD if weights used for errors
	else if(iw.eq.1.or.iw.eq.5) then
	   serr=sres	!error SD if residuals used for errors
	endif
c	call HESMAT2(theta,sdev,nfix,jfix,SMIN,SSDCV,kmax,
c     & kfit,ndth,Xobs,yobs,w,nj,juse,setx,niobs,njset,cormat)
	call HESMAT3(theta,sdev,nfix,jfix,SMIN,SSDCV,kmax,
     & kfit,ndth,Xobs,yobs,w,nj,juse,setx,niobs,njset,cormat,serr)
c (100=dimension of theta)
C TYPE CALC LS VALUES
201	continue	!jump here for curve with given param
c  After separate fits (following ifitmode=4,5 run) save the SSD and number
c of params fitted, for test of dev from parallism (in CVOUT1)
	if(ifitmode.eq.4.or.ifitmode.eq.5) then
	   do i=1,kmax
		sdevsav(i)=sdev(i)	!save for likint
	   enddo
	endif
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
	alpha=VIDEOTYP().eq.3		!alpha mode
	km2=kmax	!copy for call
	call CVOUT1(theta,noerr,iw,nev,smin,jfix,fitted,titlep,
     & sres,Sdev,titled,yobs,xobs,w,nj,jmiss,setx,niobs,njset,nodata,
     & dcurve,ndth,km2,cormat,kfit,
     & nxval,xgrp,ybar,sdm,ngrp,SSwg,ndf)
	DEALLOCATE(cormat)
c Calculate here the fitted curve for the model in use for the current set
c (jset is in COMMON) using DEFYcal (NB in NODATA case, ifitmode=2)
c
c	print 108
c108	format(/)
c No need to ask again unless nodata (in which case this is the first time)
	if(lflag.eq.0.and.nodata) then	!ask once only
	   if(.not.readini) ilog=0
	   i=ilog+1
	   print 109,i
109	   format(' Initial display type:',/,
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
	   logx=ilog.eq.1.or.ilog.eq.3.or.ilog.eq.4.or.ilog.eq.6	!for DEFYCAL
	   logy=ilog.eq.2.or.ilog.eq.3
	   logity=ilog.eq.4
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
c     If there are enough df 'within x values', do analysis of variance for
c      current set and print it
c     Remove for now -still not right
c	   call DO_AoV(xobs,yobs,w,jset,nj,niobs,njset,
c     &  	kmax,theta,setx,kfit)
	   goto 35		!display
c
	else if(ifitmode.eq.2) then	!calc curves for current pass (set #jset)
	   if(nmod.eq.26.or.nmod.eq.27) then
		iset=1		!needed in commom for ycalcv
	   endif
	   jval=jset
c If sep fits done with Hill/Langmuir, with separate Ymax fitted to
c each set, then after the sep fits can ask whether the sets should be
c normalised with respect to the separately fitted maxima, and treated
c as new data to be fitted e.g. with ifitmode=3, 4 or 5
c  No need to go to 35 when all sets done, since can't have iflag=2 here
c so just say 'hit key to display sep fits' and goto 3062 -after display
c of these, ask whether to generate normalised data -if so make it
c print it and specify how to fit it.
	   if(iflag.eq.2) then	!for separate fits AFTER ifitmode=4,5
		jval=jset+nset	!sep fits in j=nset+1, nset+2,....,2*nset
		iline(jval)=2
		ncal(jval)=ncalc
		ncurvc=ncurvc+1
		icurvc(ncurvc)=jval
	   endif
	   if(mset.eq.1) then
		if(alpha) then
		   call DCASK(
     &	   'Fit same equation as before to all curves','Y',ans)
		else
		   ans='Y'
		   call DEFOLTa(ans,defolt)
		   call QDIALOG(1,
     &	    'Fit same equation as before to all curves',
     &	    defolt,11,cans)
		   call GETINPa(cans,ans)
		endif
		sameq=ans.eq.'Y'
	  endif
c       Calculate the curves for the set just fitted
	   call DEFYcal(xmin,xmax,ncalc,logx,theta,Xcal,Ycal,jval,xgt0,
     &    setx,njset,ndimc,ndc1)
c     If there are enough df 'within x values', do analysis of variance for
c      current set and print it
c     Remove for now -still not right
c	   call DO_AoV(xobs,yobs,w,jset,nj,niobs,njset,
c     &  	kmax,theta,setx,kfit)
c       Go on to the next set
	   mset=mset+1
c When all done set iflag=3 here as signal to ask whether to normalise
c after display of separate fits (as for iopt1=2 in cvdisp) (unless
c these are separate fits that FOLLOW ifitmode=4,5)
	   if(mset.gt.nsfit) then
		if(ifitlast.le.3) then
		   iflag=3	!signals that all sep fits done in ifitmode=2
		endif
		if(alpha) then
		   print 4
c4		   format(' Type <any key> for display')
		else
		   call WDIALOG(1,'Hit any key for display',11)
		endif
		call ANYKEY
c==		goto 3062	!for display/plot
		goto 35	!display
	   endif
	   j=juse(mset)		!next set # to be fitted
	   jset=j
	   if(iflag.eq.0) iflag=1		!to tell GETEQN this is not 1st set
c Set guesses for the next set here.  Guesses for each fit separately
c are calculated in thetgues(i,j) (in GETGUESS called in CVSIMP) for
c some models (and in some cases are improved on the basis of the
c simultaneous fit, when ifitmode=4,5 initially)
	   do i=1,kmax
		theta(i)=thetgues(i,jset)
	   enddo
	   if(hillang.and.isepfit.eq.2) then
		theta(2)=thetsav(2)	!value from simultaneous fit must be used
	   endif
	   goto 31		!fit next set
c
	else if(ifitmode.ge.3.and.ifitmode.le.5) then !one pass only -calc all curves
	   iset=0
	   do 37 m=1,nsfit		!go through all fitted sets
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
	   if((nmod.eq.26.or.nmod.eq.27).and.ifitmode.eq.3.and.
     &		ip1.eq.1) then
		iset=0
		do m=1,nsfit
		   j=juse(m)
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
c Need to record ifitmode as ifitlast here, so can distinguish between
c separate fits (fitmode 2) that precede normalisation (iflag=3)
c and separate fits (fitmode 2) that FOLLOW ifitmode=4,5
46	continue		!return here after likelihood intervals for 4,5 calcs
	if(ifitmode.eq.4.or.ifitmode.eq.5) then
	   ifitlast=ifitmode	!record whether 4 or 5
	   SSDpar=smin		!save SSD for parallel (simultaneous) fit
	   kfitpar=kmax-nfix	!save number of param fitted for ditto
	   IFLAG=2		!to indicate last section of ifitmode=4,5 calcs
c  Likelihood intervals for sim fit now done at end
c	   if(alpha) then
c	     call DCASK(
c     &    'Calculate likelihood intervals for simultaneous fit','n',ans)
c	   else
c		ans='N'
c		call DEFOLTa(ans,defolt)
c		call QDIALOG(1,
c     &	'Calculate likelihood intervals for simultaneous fit',
c     &	  defolt,11,cans)
c	      call GETINPa(cans,ans)
c	   endif
c	   if(ans.eq.'Y') goto 45
         if(alpha) print 38
38	   format(/,' SEPARATE FITS TO EACH SET INDIVIDUALLY')
	   isepfit=1	!in case nmod is not 1 or -1
	   if(iabs(nmod).eq.1.or.nmod.eq.26.or.nmod.eq.27) then  !Langmuir/Hill
48		if(alpha) then
		   print 47,isepfit
47		   format(' For separate fits:',/,
     &	' (1) Allow all parameters to vary',/,
     &	' (2) Fix all maxima at value for simultaneous fit',/,
     &	' (3) Constrain maximum to be same for all curves',/,
     &	' (4) Omit the separate fits',/,
     &	' Option number [',i2,'] = ')
		   call INPUTi(isepfit)
		else
		   title='Separate fits:'
		   strings(1)='Allow all parameters to vary'
		   strings(2)='Fix all maxima at value for simultaneous fit'
		   strings(3)='Constrain maximum to be same for all curves'
		   strings(4)='Omit the separate fits'
		   nval=4
		   jline=1		!default
		   nhelp=0
		   helps(1)=' '
		   call POPMENU(nxlo,nylo,nyhi,strings,nval,ictm,icfm,icupm,
     &	      ibkm,title,helps,nhelp,jline,ans,ival)
		   isepfit=jline
		endif
		if(isepfit.lt.1.or.isepfit.gt.4) goto 48
		if(isepfit.eq.4) then
		   call WDIALOG(1,'Hit any key for display',11)
		   goto 351	!display simultaneous fit only
		endif
c           if(pon()) write(7,38)		!sep fit
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
	   else	!if not Hill/Langmuir
		if(alpha) then
		   call DCASK('Omit the separate fits','N',ans)
		else
		   ans='N'
		   call DEFOLTa(ans,defolt)
		   call QDIALOG(1,
     &	'Omit the separate fits',
     &	     defolt,11,cans)
	         call GETINPa(cans,ans)
		endif
		if(ans.eq.'Y') goto 351	!display simultaneous fit only
	   endif
c
	   ans='N'			!default
	   if(autosep) ans='Y'
	   if(alpha) then
		call DCASK('Do separate fits using all defaults',ans,ans)
	   else
		call DEFOLTa(ans,defolt)
		call QDIALOG(1,
     &	'Do the separate fits using all default guesses',
     &	     defolt,11,cans)
	      call GETINPa(cans,ans)
	   endif
	   autosep=ans.eq.'Y'
c
c	Find first set to be fitted = #jset
	   jset=juse(1)
	   mset=1
c Set initial guesses for the first set to be fitted separately -after the
c first separate fit have ifitmode=2 (or 3 if isepfit=3), so do not return here
c Init guesses done for Hill/Langmuir, power function (nmod=28), and
c polynomial (nmod=99) so far -but could be done for any.
	   hillang=nmod.le.4.or.nmod.eq.26.or.nmod.eq.27
	   if(hillang.or.nmod.eq.28.or.nmod.eq.99) then
		fixmax=hillang.and.jfix(2).eq.1   !Ymax fixed for Lang/Hill
		fix0=hillang.and.jfix(1).eq.1		!Y(0) fixed for Hill
		iset=1
		do i=1,kmax1
		   theta(i)=thetgues(i,jset)
		enddo
		titlep(3)='       K1' 	!remove 'set 1' from title
	   else if(nmod.eq.28) then		!power function
		do i=1,kmax1
		   theta(i)=thetgues(i,jset)
		enddo
	   else if(nmod.eq.99.and.ip1.eq.1) then		!polynomial -straight line case
		do i=1,2
		   theta(i)=thetgues(i,jset)
		enddo
	   endif
c
	   hillang=nmod.eq.-1.or.nmod.eq.1.or.nmod.eq.26.or.nmod.eq.27
	   if(isepfit.ne.3) then
	      ifitmode=2		!to do sep fits with original code
	      kmax=kmax1 		!reset to value without extra params (nmod the same)
		call INTCONV(imod,cnum)
		if(alpha) then
		   call DCASK(
     &	   'Fit same equation as before to all curves','Y',ans)
		else
		   ans='Y'
		   call DEFOLTa(ans,defolt)
		   call QDIALOG(1,
     &	    'Fit same equation as before to all curves',
     &	    defolt,11,cans)
		   call GETINPa(cans,ans)
		endif
		sameq=ans.eq.'Y'
		fix0=jfix(1).eq.1	!record before jfix() re-set
		do i=1,kmax		!reset
		   jfix(i)=0
		enddo
		nfix=0
c	   Choose what is to be fixed for separate fits here and set jfix()
		if(hillang) then
		   ans='N'
		   if(fix0) ans='Y'	!Y(0) fixed to 0 for simult fit
		   call DEFOLTa(ans,defolt)
		   if(iequiv.eq.0) then
			call QDIALOG(1,'Fix Y(0)=0',defolt,11,cans)
		   else
			call QDIALOG(1,'Fix Y(inf)=0',defolt,11,cans)
		   endif
		   call GETINPa(cans,ans)
		   if(ans.eq.'Y') then
			theta(1)=0.0
			jfix(1)=1
			nfix=nfix+1
		   endif
		endif
		if(isepfit.eq.2.and.hillang) then
		   if(alpha) then
			ans='Y'
		      print 41,thetsav(2)
41	      	format(
     &        ' Maximum fixed at',g13.6,' (from simult fit): O.K. [Y] ')
			call INPUTa(ans)
	      else
			ans='Y'
			call DEFOLTa(ans,defolt)
		      call REALTOCH(thetsav(2),cnum,11)
			n=nblank(cnum)
			call QDIALOG(1,
     & 	 'Maximum fixed at '//cnum(1:n)//' (from simult fit): O.K.',
     &		  defolt,11,cans)
			call GETINPa(cans,ans)
		   endif
		   if(ans.eq.'Y') then
		      theta(2)=thetsav(2)	!ymax
		      jfix(2)=1
			nfix=nfix+1
		   endif
		endif
		call SETFIX(nfix,jfix,kmax,kfit,theta,titlep,alpha,ndth)
		fixset=.true.	!so not altered in CVSIMP
	      goto 31		!repeat fitting with ifitmode=2, same nmod
	   else if(isepfit.eq.3) then
c		Redefine model and ifitmode, and set initial guesses in thetsav(i,1)
		nset1=nset	!copy for call
		call NEWEQCV(nmod,kmax,kmax1,titlep,nsfit,juse,ifitmode,
     &	  theta,thetsav,setx,njset,ndth)
c Set guesses in neweqcv!
c Need to enquire about fixing?
		ans='N'
		if(fix0) ans='Y'	!Y(0) fixed to 0 for simult fit
		call DEFOLTa(ans,defolt)
		if(iequiv.eq.0) then
		   call QDIALOG(1,'Fix Y(0)=0',defolt,11,cans)
		else
		   call QDIALOG(1,'Fix Y(inf)=0',defolt,11,cans)
		endif
		call GETINPa(cans,ans)
		if(ans.eq.'Y') then
		   theta(1)=0.0
		   jfix(1)=1
		   nfix=nfix+1
		endif
		call SETFIX(nfix,jfix,kmax,kfit,theta,titlep,alpha,ndth)
		fixset=.true.	!so not altered in CVSIMP
		if(debug()) then
		  print 431,j,(thetsav(i),i=1,kmax)
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
		if(alpha) then
      	   print 53,ndferr,ndfdev,nsfit,nfittot,kfitsep,kfitpar
		else
	         call INTCONV(ndferr,cnum)
	         call INTCONV(ndfdev,cnum1)
		   call WDIALOG(1,
     &	 'PROBLEM: d.f. for error = '//charnb(cnum)//
     &	 ' d.f. for non-parallelism = '//charnb(cnum1),12)
	         call INTCONV(nsfit,cnum)
	         call INTCONV(nfittot,cnum1)
		   call WDIALOG(1,
     &	 charnb(cnum)//' sets fitted; '//charnb(cnum1)//
     &	 ' observations fitted',12)
	         call INTCONV(kfitsep,cnum)
	         call INTCONV(kfitpar,cnum1)
		   call WDIALOG(1,
     &	 charnb(cnum)//' parameters for separate fit; '//
     &	 charnb(cnum1)//' params for parallel fit',12)
		endif
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
c	   print 700,nfit,kfitsep,SSDsep,kfitpar,SSDpar
c700	   format(' nfit,kfitsep,SSDsep,kfitpar,SSDpar = ',/,
c     &    2i8,g13.6,i8,g13.6)
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
	   if(alpha) then
            print 52,aMSerr,ndferr,aMSdev,ndfdev,vr,P
	   else
		call WDIALOG(1,
     &	 'TEST for DEVIATION FROM PARALLELISM',11)

	      call REALTOCH(aMSerr,cnum,11)
	      call INTCONV(ndferr,cnum1)
		call WDIALOG(1,
     &	'MS error (from separate fits) = '//charnb(cnum)//
     &	 ' ('//charnb(cnum1)//' d.f.)',11)
	      call REALTOCH(aMSdev,cnum,11)
	      call INTCONV(ndfdev,cnum1)
		call WDIALOG(1,
     &	'MS for diff between sep and parallel = '//charnb(cnum)//
     &	 ' ('//charnb(cnum1)//' d.f.)',11)
	      call REALTOCH(vr,cnum,11)
	      call REALTOCH(P,cnum1,11)
		call WDIALOG(1,
     &	'Var ratio = '//charnb(cnum)//
     &	 '; P = '//charnb(cnum1)//' [hit key for display]',11)
	   endif
         if(pon()) write(7,52) aMSerr,ndferr,aMSdev,ndfdev,vr,P
         if(discprt) write(8,52) aMSerr,ndferr,aMSdev,ndfdev,vr,P
52	   format(/,' TEST for DEVIATION FROM PARALLELISM',/,
     &  ' Mean square error (from separate fit) = ',g13.6,
     &  ' (',i5,' d.f.)',/,
     &  ' Mean square for diff between parallel and sep fits = ',g13.6,
     &  ' (',i5,' d.f.)',/,
     &  ' Variance ratio = ',g13.6,'   P = ',g13.6)
	endif
351	if(alpha) then
	   print 4
4	   format(' Type <any key> for display')
	else if(iflag.ne.2) then
	   call WDIALOG(1,'Hit any key for display',11)
	endif
	call ANYKEY
	goto 3062	!for display/plot
c
320	continue
c==	if(iflag.eq.2) goto 99	!now done here for ifitmode=4,5
	if(iparsav(1).lt.0.) then
	   call DCASK('Calculate likelihood intervals','N',ans)
	   if(ans.eq.'N') goto 99
	else
	   call DCASK(
     &   'Calculate likelihood intervals for simultaneous fit','N',ans)
	   if(ans.eq.'Y') then
		do i=1,kmax
		   theta(i)=thetsav(i)	!restore params for sim fit
		   sdev(i)=sdevsav(i)	!save for likint
		   jfix(i)=jfixsav(i)
		enddo
		smin=sminsav
		elmax=-smin/2.
		ifitmode=iparsav(1)
		nmod=iparsav(2)
		kmax=iparsav(3)
		kmax1=iparsav(4)
		kfit=iparsav(5)
		nfit=iparsav(6)
		ip1=iparsav(7)
		ip2=iparsav(8)
	   else
		goto 99
	   endif
	endif
45	continue
	elmax=-smin/2.
	if(iw.eq.1.or.iw.eq.5) elmax=elmax/(sres**2)	!use Sres as error var
	kmax2=kmax	!for likblk
	call LIKINT1(ELFUNCV,noerr,sdev,nmod,constr,nfit,
     & theta,jfix,titlep,ndth,Xobs,yobs,w,nj,juse,setx,niobs,njset)
c==	if(iflag.eq.2) goto 46	!complete ifitmode=4,5 calcs
	goto 320
c
99	continue
	if(VIDEOTYP().ne.3) call VIDEOMOD(3) !utility lib- back to alpha mode
c Write values to .INI
	if(cluster) goto 997
      OPEN(unit=19,file='CVFIT.INI',status='UNKNOWN',
     & access='DIRECT',form='UNFORMATTED',recl=256)
	 write(19,rec=1) ndev,ifile,ifitsav,ilog,idiskq,titlef,infil,
     &	niobs,njset,ascfil,nsfit
	CLOSE(unit=19)
c
c=	call SAVEDP
997	continue
	fitted=.false.
	CLOSE(unit=12)
	fopen=.false.
c Rather than call SAVEDP, put its content here, but ask 'another run?'
c while disc file is closed, so when this question is on screen, it
c will be possible to alt-tab to look at output file eg in notepad
	if(discprt) then
	   CLOSE(unit=8)
	endif
	ans='Y'
	call DCASK('Another run',ans,ans)
	if(discprt) then
         OPEN(unit=8,file=FILNAM,status='UNKNOWN',
     &   access='APPEND',form='FORMATTED',carriage control='FORTRAN')
	endif
	idest=1			!not first time in CVIN
	if(UC(ans).eq.'N') goto 999
	fitted=.false.
	dcurve=.false.
	fixmax=.false.
	fix0=.false.
	norm1=.false.	!for norm after sep fits (see cvin2, cvdisp)
	fixset=.false.	!params to be fixed not set yet
	sameq=.false.	!fit same eq to all sep fits that follow ifitmode=4,5
	ncurvc=0
	itit=0
	iflag=0
	isepfit=1
	guesdone=.false.
	nfix=0		!no param fixed
	do i=1,100
	   jfix(i)=0
	enddo
 	if(VIDEOTYP().ne.3) call VIDEOMOD(3) !utility lib- back to alpha mode
c	ans='Y'
c	call DCASK('Same data and weighting',ans,ans)
c	if(ans.eq.'Y') then	!this option needs work!
c	   goto 3062
c	else
	   goto 2		!another run
c	endif
c
999	continue
 	call VIDEOMOD(3)	 !utility lib- back to alpha mode
	call ENDPRINT
c
	END


