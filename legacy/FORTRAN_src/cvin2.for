	subroutine CVIN2(IW,SETVAR,AVAL,BVAL,defolt,ndev,titw,nmiss,
     & titlef,titled,titlex,titley,ilabel,yobs,xobs,w,idest,norm1,
     & fopen,ncalc,ncurvc,icurvc,ncal,iline,hdisp,yval1,nyval,printed,
     & readini,ifile,ilog,idiskq,ifitsav,iesc,ndisc,infil,ifile1,
     & niobs,njset,ndimd,ndimc,nj,jmiss,juse,setx,iver,ascfil,nodata,
     & guesdone,nxval,xgrp,ybar,sdy,sdm,ngrp,SSwg,ndf)
c
C INPUT OF DATA FOR CVFIT-window version
c Modif 01/16/96 09:23pm by adding nodata=true when used only to
c plot calculated curves with no data.
c Modif 05/05/95 09:48am by addition of check to see whether same x
c  value occurs in all sets when more than one set is being fitted (if
c  so this value can optionally be used for normalisation (for nmod=26,27
c  at present -see GETEQN, YCALCV)
c  -extra values kept in COMMON/BLOCK3/ logyfit,norm,xnorm
c 03/25/95 06:04pm add nj(),jmiss(),setx() to parameters -no longer
c in common. Also titlep() removed from CVIN, CVDISK (03/30/95 09:59pm) as
c value kept on disk is never used (must still read it from earlier versions
c of cvdat.dat
c CVIN2 07/30/91 02:54pm does not ask which eqn to be fitted -this is now done
c by separate subroutine DEFEQN, after the initial display
c CVIN1 uses new, more compact, data storage in CVDAT (as in RANTEST and VPLQ1)
c 12/07/89 04:24pm Lahey version. Changed to max of 10 data sets, 100 obs/set
c and up to 20 param to be fit
c Also include option to put axis titles in here, and store them with the
c data (may want diff axis labels for diff data sets in a file?- not yet done)
c New data file structure:
c In Lahey fortran file is not of pre-fixed size but has as many records
c as have been written. Also data files quite big, so use one record per
c SET rather than one per file (most files have only one set so will waste
c much space to write 10 sets each time)
c REC #1=nfile,istrec(100)   where istrec(if)=1st rec number for file #if
c Data for file #if is in NSET records as follows:
c istrec(if)=titled(1),nset,nj(10),titlep,setvar,setx(10),iw,
c		    Xobs(i,1),Yobs(i,1),W(i,1),titlex,titley,ilabel
c istrec(if)+1=titled(2),Xobs(i,2),Yobs(i,2),W(i,2)
c and so on up to record number
c istrec(if)+nset-1=titled(nset),Xobs(i,nset),Yobs(i,nset),W(i,nset)
c###PROBLEM is that cannot just overwrite an existing file with a new one
c as new one might have more sets
c
	real*4 xobs(niobs,njset),yobs(niobs,njset),w(niobs,njset)
	real YVAL1(1024),ynorm(niobs)		!data for histogram
	ALLOCATABLE::data			!for INWIND
	real data(:,:)
	ALLOCATABLE::xobs1,yobs1,w1,nj1,xnorms	!for normalised data
	real*4 xobs1(:,:),yobs1(:,:),w1(:,:),xnorms(:)
	integer*4 nj1(:)
	ALLOCATABLE::setx1,titled1
	real*4 setx1(:)
	character*60 titled1(:)
	integer ncal(ndimc),icurvc(ndimc),iline(ndimc)		!for calc curve
	character*79 heading,title
	character*60 titlef	!file name
	character*60 titled(njset)
	character*70 text       !for ascread1
	character*11 cnum
	integer jmiss(njset),juse(njset),nj(njset)
	real setx(njset)
	character*10 TITLEP(20)		!names of params for fit
	character*60 titw(5)		!names of weighting methods
	character ndev*2,infil*33,ascfil*33,filnam1*33,infil1*33
	ALLOCATABLE::colhead
	character*20 colhead(:)		!for ASCREAD1
	logical dcfile,auto,nodata,use1,ascinput,guesdone
	logical shortprt,nodisp,equal1,allocated
	allocatable::Xnum
	real*4 Xnum(:,:)
C THESE ARE ARRAYS TO HOLD TITLES (UP TO 60 CHAR) FOR DATA SETS
C (IN TITLED) AND FOR PARAMETERS (TITLEP)(10 CHAR)
c Also include option to put axis titles in here, and store them with the
c data
	character*40 titlex,titley	!axis labels
	character*1 ans,UC,ans1
	LOGICAL CONSTR,FLINE,SETVAR,setvar1,defolt,specsy,fixwt,addset
	logical discprt,slock,fopen,hdisp,pon,prt,zerow,printed,present
	logical readat,readini,setsd,norm1
	integer*2 jstrec(200),lstrec(200)		!for queue
	integer*4 jstrec1(100),lstrec1(100)		!for old queue
c	character filnam*13
c  Declarations for getgroup (but if these arrays to be available in main
c   program they must be declared there!)
	real*4 xgrp(niobs,njset),ybar(niobs,njset),sdm(niobs,njset)
	real*4 sdy(niobs,njset)
	integer ngrp(niobs,njset),nxval(njset)
	logical cluster,student
c
	common/dp/discprt
	COMMON/BLOCK1/constr,nset,nfit,nsfit,
     & Xv,alfa,kmax,ncomp,nmod,fline,nomit,jomit(20),jset,ifitmode
	logical logyfit,norm,plotrue
	COMMON/BLOCK3/logyfit,norm,xnorm,iostat
	COMMON/BLOCK2/ castar,X1,X2,iequiv,ip1,ip2
	common/readqplot/plotrue,istrec,iptype,kwi,kwj,ndv1,ndc1
	common/user/student,cluster,iwindows
c Oct 86- two new variables added to common/block2/
c IEQUIV=0 for equiv subunits; =1 for non-equiv (prev signalled via castar)
c IP1,IP2:
c  (1) mods 10,11 IP1=1 if param are K1,K2; =2 if param are K1, K2/K1.
c  (2) mods 9,11  IP2=1 if param is K2 or (K2/K1 if IP1=2)
c		  IP2=2 if param is k(-2) (must have ip2=1 if ip1=2)
c  (3) mods 9,11  IP2 set negative if k(+2) is a param  to be estimated
c
c
	pon()=slock()

c
c101	format(a1)
	prt=pon()		!to reduce number of pon() calls (interrupt 11 problems)
c
	plotrue=.false.
	shortprt=.false.	!no print out for very long data sets
	nodisp=.false.
	norm=.false.	!normalised to Y at given X
	iopt2=0
	titw(1)=' (1) Weights constant. Error from residuals.'
	titw(2)=' (2) Weights from specified s(Y): errors from weights'
	titw(3)=' (3) ditto: calculate weights from s(Y) = a + b*X'
	titw(4)=' (4) ditto: calculate weights from s(Y) = a + b*Y '
	titw(5)=
     & ' (5) Use specified weights for fit, but error from residuals'
c
	printed=.false.	!so data not printed twice
c
	if(idest.eq.0) goto 60		!first time
c For refit of norm data after ifitmode=2, get refit mode and start again
c For refit after normalisation, better use the same sets as originally
c fitted so now jump tot 2021, not 202
	if(idest.eq.2) then
	   idest=0
	   idef=6	!default
2031	   print 203,idef
203	   format(/,
     & ' Data for the sets just fitted has been replaced',
     & ' by normalised data'
     & ' Re-fit the normalised curve (s) as follows:',/,
     & ' (3) Fit selected data sets simultaneously with one equation',/,
     & ' (4) Fit selected data sets to estimate relative potencies',/,
     & ' (5) Fit selected data sets to estimate antagonist KB (Schild)',
     & /,'   (antagonist conc specified as ''set variable'')',/,
     & ' (6) Pool selected data sets and fit as one set',/,
     & '  Option number [',i2,'] = ')
	   call INPUTi(idef)
	   if(idef.ge.3.and.idef.le.6) then
	      ifitmode=idef
		idest=0
      	if(pon()) write(7,41)
      	if(discprt) write(8,41)
41		format(/,
     &' Data re-fitted after normalising to maxima from separate fits',
     &	/)
c=		goto 202		!start fitting from scratch
		goto 2021		!start fitting again, using same sets
	   else
		goto 2031
	   endif
	endif
c
c If normalised curve has been done then original data overwritten
c so must re-read data.  Safer to do this anyway (eg weights altered
c in some ifitmode=3 options)
	idest=0              !diff data
c
60	specsy=.true.	!if not reset, specify weight as s(y)
c
	defolt=.true.		!not used now
c
c Check which disc to use (if any) and read data
	nodata=.false.
	ascinput=.false.
200	print 20
20	format(
     & ' (1) Read data from CVDAT data file',/,
     & ' (2) Read data from ASCII file',/,
     & ' (3) Read data from plot queue file [only for new files]',/,
     & ' (4) Input new data from keyboard',/,
     & ' (5) No data (or fit) -plot calculated curve(s) only',/,
     & ' Option number [1] = ')
	iopt=1
	call INPUTi(iopt)
	if(iopt.eq.4) then
	   readat=.false.
	   goto 102
	else if(iopt.eq.3) then		!READ FROM PLOT QUEUE
	   filnam1='\PLOTQ.DAT'
310	   call DISCNUM1(id,ndev,-1,0)
	   if(id.eq.-1) then
		iesc=99
		RETURN			!abort if ESC hit
	   endif
	   INFIL1=ndev//FILNAM1
250	   call TITENT0(
     &   'Name and path of PLOTQ data file:',infil1,33,.false.)
	   INQUIRE(file=infil1,exist=present,flen=len)
	   if(.not.present.or.len.eq.0) then
	      call BELL(3)
	      print 240,infil1
240	      format(1x,a33,' does not exist')
	      present=.false.		!so asks for ndev
	      goto 310
	   else
            OPEN(unit=11,file=INFIL1,status='OLD',
     &      access='DIRECT',form='UNFORMATTED',recl=1024)
	   endif
	   read(11,rec=1) nplot,jstrec,lstrec,iver
	   if(iver.ne.1001.and.iver.ne.1100) then
		call BELL(2)
		print 217
217		format(
     &	'&EITHER ',/,
     &	'  (1) this is not a Plot Queue file',/,
     &	' OR ',/,
     &	'  (2)it is plot queue with old format',/,
     &	' Which option is correct [1] = ')
		iopt1=1
		call INPUTi(iopt1)
		if(iopt1.ne.2) then
 		   goto 250
		else
		   print 2171
2171		   format(' Convert old queue on disk to new format [N] ? ')
		   ans='N'
		   call INPUTa(ans)
	 	   if(UC(ans).ne.'Y') then
		      call CONVQREC(jstrec,jstrec1,lstrec,lstrec1,11,0)
		   else
	  	      call CONVQREC(jstrec,jstrec1,lstrec,lstrec1,11,1)
		   endif
	      endif
	   endif
	   print 5021,nplot
5021     format(' Plot queue has plots=',i3)
	   CLOSE(unit=11)
	   irecl=1024		!for plotq
         OPEN(unit=11,file=INFIL1,status='UNKNOWN',
     &    access='DIRECT',form='UNFORMATTED',recl=irecl)
	   jopt=1
4332	   print 5011
5011     format(' Enter plot number = ')
	   call INPUTi(jopt)
	   iplot=jopt
	   istrec=int4(jstrec(iplot))	!1st record no for plot # iplot
	   nrect=int4(lstrec(iplot))-istrec+1		!total number of records
	   nrect=nrect-2				!number of records for Xval etc
	   read(11,rec=istrec) iptype
	   j=iabs(iptype)
	   if(iptype.eq.14) then
		read(11,rec=istrec) iptype,ndv1,ndimd,ndc1,ndimc
		kwi=100	!dimensions for weight()
		kwj=10
		kmax=20
		if(ndv1.eq.0) ndv1=1
		if(ndimd.eq.0) ndimd=1
		if(ndc1.eq.0) ndc1=1
		if(ndimc.eq.0) ndimc=1
	   else if(iptype.eq.15) then
	     read(11,rec=istrec)iptype,ndv1,ndimd,ndc1,ndimc,kwi,kwj,kmax
	   else
	      ndv1=2048	!dimensions as for earlier versions
		ndimd=10
		ndc1=2048
		ndimc=10
		kwi=100	!dimensions for weight()
		kwj=10
		kmax=20
	   endif
	   if(ndv1.le.0) ndv1=1		!in case no data curves
	   if(ndimd.le.0) ndimd=1		!in case no data curves
	   if(ndc1.le.0) ndc1=1		!in case no calc curves
	   if(ndimc.le.0) ndimc=1		!in case no calc curves
	   if(kwi.le.0) kwi=1
	   if(kwj.le.0) kwj=1
	   if(iptype.ge.15.and.iptype.le.25) then
		CLOSE(unit=11)
            OPEN(unit=11,file=INFIL1,status='UNKNOWN',
     &		access='TRANSPARENT')
		istrec=1 + int4(jstrec(iplot)-1)*1024  !start rec when transparent
	     read(11,rec=istrec) iptype,ndv1,ndimd,ndc1,ndimc,kwi,kwj,kmax
	   endif
	   if(iptype.eq.2) then
	      print 5013
5013         format(' Histogram (not yet done: try again')
		 goto 4332
	   endif
	   plotrue=.true.
	   goto 999
	else if(iopt.eq.5) then
	   nodata=.true.
	   ifitmode=2
	   print 40
40	   format(' Number of curves to be plotted [1] = ')
	   ncurvc=1
	   call INPUTi(ncurvc)
	   nset=ncurvc
	   ncalc=501
	   do i=1,ncurvc
		icurvc(i)=i
	 	ncal(i)=ncalc
	      iline(2)=i-1	!=0 (continuous) for curve 1
	   enddo
	   do i=1,njset
		jmiss(i)=0	!use all
	   enddo
	   goto 999		!return
	else if(iopt.eq.1) then
	   ifile=1		!so CVDISK asks for IFILE and reads it
	   call CVDISK(Xobs,Yobs,W,nset,titlep,setvar,iw,ilabel,
     &    titlef,titled,titlex,titley,ifile,nfile,FOPEN,ndev,hdisp,
     &    prt,readini,iesc,readat,ndisc,infil,ifile1,
     &    niobs,njset,nj,setx,iver)
	   if(.not.readat) goto 102
	   if(iesc.eq.98) RETURN		!reallocate arrays
	   if(iesc.eq.99) RETURN		!discnum abort
	else if(iopt.eq.2) then
	   iver=1003	!otherwise not defined
	   ascinput=.true.
	   ioff=-1
	   ilinhead=-1
c       Get filename, number of rows and cols
	   ncols=3*njset	!safe size for input to ascread1
	   ALLOCATE(colhead(ncols))
	   call ASCREAD1(ioff,ilinhead,ncols,nrows,len,
     &    ascfil,colhead,text,lentext,dcfile)
	   nd1=nrows
	   nd2=ncols
	   ALLOCATE(xnum(nd1,nd2))
	   call ASCREAD2(ioff,ncols,nrows,len,
     & 	ascfil,xnum,nd1,nd2)
	   nc1=ncols
	   if(nc1.gt.6) nc1=6	!print only first 6 cols
	   if(dcfile) then
		if(lentext.gt.0) then
		   print 211,text
211		   format(' Title: ',a70)
		endif
		titlef=text(1:60)
	   else
		do i=1,60
		   titlef(i:i)=' '
		enddo
	   endif
	   print 21,nrows,ncols,nc1
21	   format(1x,i3,' rows of data in ',i3,
     &	' columns: First 3 rows and ',i2,' columns are:')
	   print 212,(colhead(i)(1:13),i=1,nc1)
212	   format(1x,6a13)
	   do i=1,3		!print first 3 lines
		print 22,(xnum(i,j),j=1,nc1)
22		format(6g12.5)
	   enddo
c Print file name
         print 27,ascfil
         if(pon()) write(7,27) ascfil
         if(discprt) write(8,27) ascfil
27	   format(/,' Data from ASCII file: ',a33)
c Check what the ascii data represent!
	   nset=0
	   setsd=.false.
	   if(ncols.eq.1) then
		nset=1
		ncset=1
	   else
30		continue
		if(mod(ncols,2).eq.0) ncset=2		!default
		if(mod(ncols,3).eq.0) ncset=3		!default
		nset=ncols/ncset				!default
		print 28,nset
28		format(
     &     ' Number of sets to be defined from ASCII data: Nset ['
     &	,i3,'] = ')
		call INPUTi(nset)
		ncset=ncols/nset
		print 29,ncset
29		format(
     &	 ' Number of columns per set [',i2,'] = ')
		call INPUTi(ncset)
		if(ncset*nset.ne.ncols) then
		   call BELL(2)
		   goto 30
		endif
	   endif
c
c No need to ask about every column in most cases -do so only if auto=false
	   auto=.false.
	   if(ncset.eq.2) then
		ans='Y'
		call DCASK(
     &	'Columns 1, 2 = X, Y, and so on . . .',ans,ans)
		auto=ans.eq.'Y'
	   else if(ncset.eq.3) then
		ans='Y'
		call DCASK(
     &	'Columns 1, 2, 3 = X, Y, SD and so on . . .',ans,ans)
		auto=ans.eq.'Y'
	   endif

	   do nc=1,ncols
		if(ncols.eq.1) then
		   iopt=2
		   iset=1
		else
c define default as iopt=1,2, 1,2.. when ncols=2 or 4 or 6...;
c			  iopt=1,2,3, 1,2,3,.. for ncols=3 or 6...
		   iset=1 + ((nc-1)/ncset)
		   iopt=1+mod((nc-1),ncset)
		endif
		nj(iset)=nrows	!rows per set (same number for all sets??)
25		continue
		if(.not.auto) then
	         print 23,iset,nc,nc,nc,nc,iopt
23	         format(/,
     &	   ' Set # ',i3,/,
     &	   ' (1) column ',i3,' is x value',/,
     &	   ' (2) column ',i3,' is y value',/,
     &	   ' (3) column ',i3,' is standard deviation',/,
     &	   ' (4) Change set number',/,
     &	   ' (5) Discard column ',i3,/,
     &	   '  Option number [',i2,'] = ')
			call INPUTi(iopt)
		endif
		if(iset.gt.nset) nset=iset	!define nset
		if(iopt.eq.5) goto 26		!next column
		if(iopt.eq.4) then
		   print 24
24		   format(' Set number = ')
		   call INPUTi(iset)
		   goto 25
		else if(iopt.eq.1) then
		   do i=1,nrows
			xobs(i,iset)=xnum(i,nc)
		   enddo
		   if(dcfile) titlex(1:20)=colhead(nc)
		else if(iopt.eq.2) then
		   do i=1,nrows
			yobs(i,iset)=xnum(i,nc)
		   enddo
		   if(dcfile) titley(1:20)=colhead(nc)
		else if(iopt.eq.3) then
		   setsd=.true.
		   do i=1,nrows
			sy=xnum(i,nc)
			if(sy.lt.1.e18) then
			   wt=1.0/(sy*sy)
			else
			   wt=0.0
			endif
			w(i,iset)=wt
		   enddo
		endif
		if(.not.setsd) then
		   do i=1,nrows
			w(i,iset)=1.0
		   enddo
		endif
26		continue
	   enddo		!end of loop for nc=1,ncols
	   DEALLOCATE(xnum)
	   DEALLOCATE(colhead)
	   iw=1		!default weighting method (not specified yet!)
	   if(setsd) iw=2
	   goto 71		!display without asking
	endif
c
	if(nset.gt.15.or.nj(1).gt.100) then
	   ans='Y'
	   call DCASK('Omit printing of observations',ans,ans)
	   shortprt=ans.eq.'Y'
	   if(nset.gt.15) then
		ans='Y'
		call DCASK('Omit also display in blue windows',ans,ans)
		nodisp=ans.eq.'Y'
	   endif
	endif
	print 112
112	FORMAT(' (To omit an observation, set its weight to zero)',/)
	ans1='n'
	call DCASK('View/alter the data',ans1,ans)
	if(ans.eq.'N') then
	   if(iw.ge.1.and.iw.le.5) then
		print 36,titw(iw)
36		format(' Weighting now: ',a60,/,' Change it [N] ? ')
		ans='N'
		call INPUTa(ans)
	   else
		print 37
37		format(' Weighting not defined')
		ans='Y'
	   endif
	   if(UC(ans).eq.'Y') then
		fixwt=.true.
		goto 1027
	   else
		goto 114	!skip view/alter
	   endif
	endif
c
c loop to view each set
c NB with 4 cols window shows s(y) and weights always, and with ICONST=3
c cols 3 are 4 constrained in INWIND.
71	continue
	fixwt=.false.
	if(nodisp) goto 74
	heading=
     & '   X value       Y value    Standard dev   Weight'
	if(hdisp) heading='   X value'
c
c	if(iver.ge.1002.and.(.not.ascinput)) then
	if(iver.ge.1002) then
	   print 322,titlef
	   if(pon()) write(7,322) titlef
	   if(discprt) write(8,322) titlef
322	   format(/,' File: ',a60)
	endif
	n1=1
	n2=nset
711	continue	!jump here with n1=n2=2 to display pooled norm data
	do 113 j=n1,n2
	ncol=4
	if(hdisp) ncol=1
c Allocate 20 extra rows in case new lines added in inwindv
	nd1=nj(j)+20	!extra rows in case lines added in inwindv
	if(allocated(data)) DEALLOCATE(data)
	ALLOCATE(data(nd1,ncol))
	iflag=ncol	!so input data displayed in all cols of window
	call INTCONV(j,cnum)
	title=' SET '//charnb(cnum)//': '//charnb(titled(j))
	if(iopt2.eq.1) then
	   title=' Pooled normalised data (now in set 2)'
	endif
	do i=1,nj(j)
	   data(i,1)=xobs(i,j)
	   if(.not.hdisp) then		!X only for histograms
		data(i,2)=yobs(i,j)
		data(i,4)=w(i,j)
		if(w(i,j).gt.1.e-37) then
		   sd=sqrt(1.0/w(i,j))
		else
		   sd=1.e36
		endif
		data(i,3)=sd
	   endif
	enddo
	nidisp=nj(j)
	if(nidisp.gt.20) nidisp=20	!best value- up to 20?
	iconst=3	!constrain cols 3,4
	if(hdisp) iconst=0	!no constraint
	nr1=-1		!so window initially in default position
c
	call INWINDv(data,nj(j),ncol,nidisp,title,heading,
     & nr1,nc1,iflag,iconst,nd1,ncol)
c
	if(iflag.eq.-1) then
	   DEALLOCATE(data)
	   goto 74		!F5=abort -no more INWIND display
	endif
c
c  Reassign data (if it has been altered)
	if(iflag.eq.1) then
	   do i=1,nj(j)
	    xobs(i,j)=data(i,1)
	    if(.not.hdisp) then		!X only for histograms
	      yobs(i,j)=data(i,2)
	      w(i,j)=data(i,4)
	    endif
	   enddo
	endif
	DEALLOCATE(data)
113	continue		!end of jth set loop
	if(norm.and.iopt2.eq.1) goto 712	!after redisplay of pooled norm data
	if(norm) goto 75	!after redisplay of data norm to Y at given X
c
c ALTERATIONS
c
74	continue
	if(hdisp) goto 1101	!skip alterations,weights for histo
	if(iw.lt.1.or.iw.gt.5) iw=1
	print 1301,titw(iw)
	if(prt) write(7,1301) titw(iw)
      if(discprt) write(8,1301) titw(iw)
1301	format(/,' Weighting method = ',a60,/)
c
	print 72
72	format(' Options to alter data:',/,
     & ' [0] No changes ',/,
     & ' (1) Change weighting method',/,
     & ' (2) Add another data set',/,
     & ' (3) Remove LAST data set',/,
     & ' (4) Alter a set title',/,
     & ' (5) Add another observation',/,
     & ' (6) Remove an observation',/,
     & ' (7) Alter set variables',/,
     & ' (8) Alter file title',/,
     & ' (9) Merge two files in to one',/,
     & ' Option number [0] = ')
	ialt=0
	call INPUTi(ialt)
	if(ialt.eq.0) goto 1101		!print and store on disk?
	if(ialt.lt.0.or.ialt.gt.9) goto 74
	goto (81,82,83,84,85,85,87,88,89) ialt
c
c Change weights
81	continue
	fixwt=.true.
	goto 1027
c
c Add set
82	continue
	print 621
621	format(' (1) Type in data',/,
     & ' (2) Duplicate an existing data set',/,
     & ' Option number = ')
	iopt=1
	call INPUTi(iopt)
	if(iopt.eq.1) then
	   nset=nset+1
	   j1=nset		!so only new set read in below
	   addset=.true.
	   goto 204
	else if(iopt.eq.2) then
	   print 90
90	   format(' Duplicate set number = ')
	   call INPUTi(j)
	   nset=nset+1
	   if(setvar) then
		print 332,nset
		call INPUTr(setx(nset))
	   endif
	   nj(nset)=nj(j)
	   titled(nset)=titled(j)
	   do 91 i=1,nj(j)
		xobs(i,nset)=xobs(i,j)
		yobs(i,nset)=yobs(i,j)
		w(i,nset)=w(i,j)
91	   continue
	endif
	goto 74	!more changes?
c
c Remove last set
83	continue
	nset=nset-1
	goto 74	!more changes?
c option to change set title
84	continue
	print 501
501	format(' Alter title for set number = ')
	call INPUTi(j)
	nlen=60		!actually 60 char, but want not > 40 for iver=1002
	if(iver.eq.1002) nlen=40
	call TITENT0('Enter title for set:',titled(j),nlen,.false.)
	n=NBLANK3(titled(j),60)		!ensure ASCII 0 at end
	goto 74	!more changes?
c
85	continue
	j=1
	if(nset.gt.1) then
	   print 512
512	   format(' Alter set number = ')
	   call INPUTi(j)
	endif
	if(ialt.eq.5) then
	   nj(j)=nj(j)+1
	   if(iw.eq.1) then
		print 515,nj(j),j,nj(j),j
515   	format('&Xobs(',2i3,'), Yobs(',2i3,') = ')
		call INPUT2r(xobs(i,j),yobs(i,j))
		w(i,j)=1.0
	   else
		print 517,nj(j),j,nj(j),j,nj(j),j
517   	format(
     & '&Xobs(',2i3,'), Yobs(',2i3,'),w(',2i3,') (with dec point) = ')
		read 518,xobs(i,j),yobs(i,j),w(i,j)
518		format(3g13.6)
	   endif
	else if(ialt.eq.6) then
520	   print 519,nj(j)
519	   format('  -remove observation # (up to',i3,') = ')
	   call INPUTi(i)
	   if(i.lt.1.or.i.gt.nj(j)) goto 520
         if(prt) write(7,513) i,xobs(i,j)
         if(discprt) write(8,513)  i,xobs(i,j)
513	   format(' Observation #',i4,' omitted: x= ',g13.6)
	   n=0
	   do 514 k=1,nj(j)
		if(k.eq.i) goto 514
		n=n+1
		xobs(n,j)=xobs(k,j)
		yobs(n,j)=yobs(k,j)
		w(n,j)=w(k,j)
514	   continue
	   nj(j)=nj(j)-1
	endif
	goto 74 		!more changes?
c
87	continue
	if(.not.setvar) then
	   call DCASK('Define set variables','y',ans)
	   setvar=ans.eq.'Y'
	   if(setvar) then
		do 872 j=1,nset
		 print 873,j
873		 format(' Set variable for set #',i2,' = ')
		 call INPUTr(setx(j))
872		continue
	   endif
	   goto 74	!more changes?
	endif
	if(setvar) then
	   print 871
871	   format(' Alter set variable for set number [0 to skip] = ')
	   call INPUTi(j)
	   print 873,j
	   call INPUTr(setx(j))
	endif
	goto 74	!more changes?
88	continue
	if(iver.eq.1001) goto 74	!has no file title!
	nlen=60		!default as read from .INI (last titlef used)
	call TITENT0('Enter title for whole file:',titlef,nlen,.false.)
	goto 74	!more changes?
c
c Merge present file (#ifile) with another
89	continue
	print 891,ifile
891	format(' File ',i3,' already read: file # to merge with it = ')
	call INPUTi(ifile2)
	ALLOCATE(xobs1(niobs,njset),yobs1(niobs,njset),w1(niobs,njset))
	ALLOCATE(nj1(njset),setx1(njset),titled1(njset))
	call CVDISK(Xobs1,Yobs1,W1,nset2,titlep,setvar1,iw,ilabel,
     & titlef,titled1,titlex,titley,ifile2,nfile,FOPEN,ndev,hdisp,
     & prt,readini,iesc,readat,ndisc,infil,ifile1,
     & niobs,njset,nj1,setx1,iver)
c
	do j=nset+1,nset+nset2		!new set numbers
	   j1=j-nset			!=1,..,nset2
	   do i=1,nj1(j1)			!obs numbers in new set j
		Xobs(i,j)=Xobs1(i,j1)
		Yobs(i,j)=Yobs1(i,j1)
		w(i,j)=w1(i,j1)
	   enddo
	   setx(j)=setx1(j1)
	   titled(j)=titled1(j1)
	enddo
	do j=nset+1,nset+nset2		!new set numbers
	   j1=j-nset			!=1,..,nset2
	   nj(j)=nj1(j1)
	enddo
	DEALLOCATE(xobs1,yobs1,w1)
	DEALLOCATE(nj1,setx1,titled1)
c Set new values
	nset=nset+nset2
	nlen=60		!default = previous title
	call TITENT0('Enter title for new file:',titlef,nlen,.true.)
	ifile=0		!so CVDISK asks for IFILE and writes it
	call CVDISK(Xobs,Yobs,W,nset,titlep,setvar,iw,ilabel,
     & titlef,titled,titlex,titley,ifile,nfile,FOPEN,ndev,hdisp,
     & prt,readini,iesc,readat,ndisc,infil,ifile1,
     & niobs,njset,nj,setx,iver)
	CLOSE(unit=12)
	fopen=.false.
	goto 74	!more changes?
c
c End of alterations
c
c TYPE YOBS,XOBS,W,NSET,NJ FROM KEYBOARD
c NB all new files typed in here will have iver=1003 and must be written
c to a CVDAT.DAT that has iver=1003 in record #1.
102	continue
	iver=1003
	nlen=60		!default = previous title
	call TITENT0('Enter title for whole file:',titlef,nlen,.true.)
	print 103
103	FORMAT(' Number of data sets in this file = ')
	call INPUTi(nset)
	j1=1		!get all sets from 1 to NSET
c jump to 204 to add a set
204	continue
c
	do 111 j=j1,nset
	print 106,J
106	FORMAT('   Number of obs in set ',i3,' [0 to set in window] = ')
	call INPUTi(nj(j))
	if(nj(j).le.0) nj(j)=-1
c=	nlen=40		!actually 60 char, but want not > 40 for iver=1002
c=	if(iver.ne.1002) nlen=60
	nlen=60
	titled(j)=
     & "                                                           "
	call TITENT0('Enter title for set:',titled(j),nlen,.true.)
	n=NBLANK3(titled(j),60)		!ensure ASCII 0 at end
111	CONTINUE
	SETVAR=.FALSE.
	print 331
331	FORMAT(' Set variable needed [N]? ')
C SET VARIABLE=A NUMBER (EG CONC) FOR EACH SET
	ans='N'
	call INPUTa(ans)
	if(UC(ans).EQ.'Y') SETVAR=.TRUE.
c specify weighting (same for all sets)
1027	continue
	if(hdisp) goto 333	!skip weights
c In orig version iw=1,5,2,3,4 respectively for options below
c	titw(1)=' (1) Weights constant. Error from residuals.'
c	titw(2)=' (2) s(Y) (to be typed in)'
c	titw(3)=' (3) s(Y)=a + b*X    (a,b to be typed in)'
c	titw(4)=' (4) s(Y)=a + b*Y    (a,b to be typed in)'
c	titw(5)=' (5) Arbitrary relative weights. Error from residuals'
	print 1026,(titw(i),i=1,5)
1026	FORMAT(' Weights to be given as:',5(/,1x,a60),/,
     & ' Option number [1] = ')
	iw=1		!default
	call INPUTi(iw)
	if(iw.lt.1.or.iw.gt.5) goto 1027
c
	if(iw.eq.3.or.iw.eq.4) then
	   print 1023
1023	   FORMAT('  values of a,b = ')
	   call INPUT2r(aval,bval)
         if(pon()) write(7,1024)
         if(discprt) write(8,1024)
1024	   format(' a, b = ',2g13.6)
	endif
c
c Section done if existing weighting method is being changed. For iw=1,3,4
c values of s(y) and weights can be calc now and displayed; if iw=2,5 then
c display data with last 2 cols blank so s(y), or w, (whichever preferred)
c can be typed in
c If existing weights being altered then check whether to leave zero weights
c as they are (if there are any)
	if(fixwt) then
	   zerow=.false.
	   do 731 j=1,nset
	   do 731 i=1,nj(j)
		if(w(i,j).eq.0.) zerow=.true.
731	   continue
	   if(zerow) then
		print 732
732		format(' Leave the zero weights unchanged [Y] ? ')
		ans='Y'
		call INPUTa(ans)
		if(UC(ans).eq.'N') zerow=.false.
	   endif
	   if(iw.eq.1.or.iw.eq.3.or.iw.eq.4) then		!set w(i,j) here
	      do 73 j=1,nset
	      do 73 i=1,nj(j)
		if(zerow.and.w(i,j).eq.0.) goto 73		!do not reset w(i,j)
		if(iw.eq.1) sy=1.0
		if(iw.eq.3) sy=aval+bval*xobs(i,j)
		if(iw.eq.4) sy=aval+bval*yobs(i,j)
		if(sy.lt.1.e18) then
		   w(i,j)=1.0/(sy*sy)
		else
		   w(i,j)=0.0
		endif
73	      continue
	      iflag=4	!numbers shown in all 4 cols
	   else if(iw.eq.2.or.iw.eq.5) then
		iflag=2	!last 2 cols blank to enter w or s(y)
	   endif
	   goto 71		!display window
	endif
c
	if(iw.eq.2.or.iw.eq.5) then
33	   print 32
32	   format(' Specify error as',/,
     &    ' (1) standard deviation, s(y)',/,
     &    ' (2) weight, 1/var(y)',/,
     &    ' Option number [1]= ')
	   i=1
	   call INPUTi(i)
	   if(i.eq.0) i=1
	   if(i.ne.1.and.i.ne.2) goto 33
	   specsy=i.eq.1
	endif
c
	if(iw.eq.1.or.iw.eq.3.or.iw.eq.4) then
	   heading='    X value       Y value'
	   ncol=2
	else if(specsy.and.(iw.eq.2.or.iw.eq.5)) then
	   heading='    X value       Y value     Standard dev '
	   ncol=3
	else if(.not.specsy.and.(iw.eq.2.or.iw.eq.5)) then
	   heading='    X value       Y value       Weight     '
	   ncol=3
	endif
c
C  START LOOPS TO READ IN DATA
c
	if(.not.addset) j1=1	!read all sets
333	continue		!jump here for hdisp
c
	print 322,titlef
	if(pon()) write(7,322) titlef
	if(discprt) write(8,322) titlef
c322	format(/,' File: ',a60)
	do 105 j=j1,nset
	print 1093,j
1093	FORMAT(' Data set #',I3,/)
	if(setvar) then
	   print 332,j
332	   format('&Value of set variable for set #',i3,' = ')
	   call INPUTr(setx(j))
	endif
	call INTCONV(j,cnum)
	title=' SET '//charnb(cnum)//': '//charnb(titled(j))
c	heading='    X value       Y value     Standard dev '
c	ncol=3		!defined above
	if(hdisp) then
	   heading='   X value'
	   ncol=1
	endif
	if(nj(j).ge.1) then
	   nd1=nj(j)
	else
	   nd1=niobs	!if nj() not set yet
	endif
	nd1=nd1+20	!extra rows in case lines added in inwindv
	ALLOCATE(data(nd1,ncol))
	iflag=0	!display initially blank
	nidisp=nj(j)
	if(nidisp.gt.20) nidisp=20	!best value- up to 20?
c (NB if nj(j)=-1 so end of data set by # in INWIND, nidisp,ni are reset
c internally in INWIND)
c
	iconst=3	!constrain cols 3,4
	if(hdisp) iconst=0	!no constraint
	nr1=-1		!so window initially in default position
c
	call INWINDv(data,nj(j),ncol,nidisp,title,heading,
     & nr1,nc1,iflag,iconst,nd1,ncol)
c
c allocate data
	do 107 i=1,nj(j)
	xobs(i,j)=data(i,1)
	if(hdisp) goto 107		!X only for histograms
	yobs(i,j)=data(i,2)
c  calc weights
	if(.not.specsy.and.(iw.eq.2.or.iw.eq.5)) then
          w(i,j)=data(i,3)
	else
	   if(iw.eq.1) sy=1.0
	   if(iw.eq.3) sy=aval+bval*xobs(i,j)
	   if(iw.eq.4) sy=aval+bval*yobs(i,j)
	   if(specsy.and.(iw.eq.2.or.iw.eq.5)) sy=data(i,3)
	   if(sy.lt.1.e18) then
		w(i,j)=1.0/(sy*sy)
	   else
		w(i,j)=0.0
	   endif
	endif
107	continue
c Print out values
	if(.not.shortprt) then
	   call PRINTJ(j,titled,setx,hdisp,setvar,prt,nj,w,xobs,yobs,
     &    niobs,njset,norm1)
	endif
	printed=.true.	!so data not printed again below
	DEALLOCATE(data)
105	continue	!end of jth set loop
c
1101	continue
	if(hdisp) titley='  frequency '
	ilabel=1
	print 720,titlex,titley
720	format(' Axis labels:',/,1x,a40,/,1x,a40,/,
     & ' O.K. [Y] ? ')
	ans='Y'
	call INPUTa(ans)
	if(UC(ans).eq.'N') then
	   print 721
721	   format('&Specify X axis label [Y] ? ')
	   ans='Y'
	   call INPUTa(ans)
	   if(UC(ans).ne.'N') call LABENT('Enter label for X axis:',
     &	titlex,40)
	   print 722
722	   format(' Specify Y axis label [Y] ? ')
	   ans='Y'
	   call INPUTa(ans)
	   if(UC(ans).ne.'N') call LABENT('Enter label for Y axis:',
     &	titley,40)
	   ilabel=1
	   print 12
12	   format(/)
	endif
	if(readat) then
	   print 110,ifile
110	   FORMAT(
     &  '&Data is from file ',i3)
	   call DCASK(' : re-write data to disc','n',ans)
	else
	   call DCASK(' Store these data on disc','y',ans)
	endif
	if(ans.eq.'Y') then
	   iesc=0
c	   if(readat) iesc=ifile !default file to write to -ifile1 now a param
	   ifile=0		!so CVDISK asks for IFILE and writes it
	   call CVDISK(Xobs,Yobs,W,nset,titlep,setvar,iw,ilabel,
     &    titlef,titled,titlex,titley,ifile,nfile,FOPEN,ndev,hdisp,
     &    prt,readini,iesc,readat,ndisc,infil,ifile1,
     &    niobs,njset,nj,setx,iver)
	   if(iesc.eq.99) RETURN		!discnum abort
	   close(unit=12)
	   fopen=.false.
	endif
c
C NOW HAVE DATA. START FITTING
114	CONTINUE
c First always sort all sets into increasing order of x values
c  -this both ensures correct plotting (joining of points) and
c   simplifies finding of replicates at 'same' x value
c
	do j=1,nset
	   call SORTr3(Xobs,Yobs,w,j,nj(j),.true.,
     &	   niobs,njset,niobs,njset)
	enddo
c
c	For histogram data have option later to superimpose histos
c	for several sets in the same file. Add option to make single histo
c	but this histo may pool data from seceral sources, with all data
c	put into YVAL1 rather than in XOBS; in this case set NSET=-1 as
c	signal that data is in YVAL1, with totlal number of obs in NYVAL
	if(hdisp) then
116	   print 115
115	   format(
     &	' (1) Use these data set(s) as they are',/,
     &	' (2) Pool with data from 2 or more sets in this file',/,
     &	' (3) Concatenate data from another file',/,
     &	' Option number [1] = ')
	   iopt=1
	   call INPUTi(iopt)
	   if(iopt.le.0) iopt=1
	   if(iopt.eq.1) RETURN		!no fitting for histos yet
	   if(iopt.eq.2) then
		if(nset.eq.1) goto 116		!no other sets to pool
		print 117
117		format('&Number of sets to be pooled (2 or more) = ')
	      call INPUTi(ns1)
		k=0
		do 118 n=1,ns1
		print 146
146		format('& pool set number = ')
	      call INPUTi(j)
      	print 148,nj(j),j
      	if(prt) write(7,148) nj(j),j
      	if(discprt) write(8,148) nj(j),j
148		format(1x,i5,' values concatenated from set ',i3)
		do 118 i=1,nj(j)
		k=k+1
		yval1(k)=xobs(i,j)
118		continue
		nyval=k
		nset=-1
	      RETURN		!no fitting for histos yet
	   else if(iopt.eq.3) then
		j=1
		if(nset.gt.1) then
		   print 149
149		   format('& from present file use set number = ')
	   	   call INPUTi(j)
		endif
		k=0
      	print 148,nj(j),j
      	if(prt) write(7,148) nj(j),j
      	if(discprt) write(8,148) nj(j),j
		do 150 i=1,nj(j)
		k=k+1
		yval1(k)=xobs(i,j)
150		continue
		print 147
147		format(' Number of files to be pooled (1 or more) = ')
	   	call INPUTi(ns1)
		do 144 n=1,ns1
		ifile=1		!so CVDISK asks for IFILE and reads it
	      call CVDISK(Xobs,Yobs,W,nset,titlep,setvar,iw,ilabel,
     &       titlef,titled,titlex,titley,ifile,nfile,FOPEN,ndev,hdisp,
     &       prt,readini,iesc,readat,ndisc,infil,ifile1,
     &       niobs,njset,nj,setx,iver)
	      if(iesc.eq.99) RETURN		!discnum abort
		j=1

		if(nset.gt.1) then
		   print 145
145		   format('& pool set number = ')      !which set in new file?
	 	   call INPUTi(j)
		endif
      	print 148,nj(j),j
      	if(prt) write(7,148) nj(j),j
      	if(discprt) write(8,148) nj(j),j
		do 144 i=1,nj(j)
		k=k+1
		yval1(k)=xobs(i,j)
144		continue
		nyval=k
		nset=-1
	      RETURN		!no fitting for histos yet
	   endif
	endif
c
c Initially all sets used
	do i=1,njset
	   jmiss(i)=0	!use all
	enddo
c
	ncalc=501
	if(nset.eq.1) then
c	   ncal(1)=1
	   ifitmode=1	!if only one set
	   jset=1	!fit set 1 if only one set
	   ncurvc=1
	   ncal(1)=ncalc
	   icurvc(1)=1
	   iline(1)=0	!continuous
	   nsfit=1
	   juse(1)=1
c	   goto 98
	   goto 97	!get nxval
	endif
c
c  For nset>1 ask which to fit
1201	continue
	idef=2	!default
	if(readini) idef=ifitmode
	print 120,idef
120	format(/,' CHOOSE FITTING MODE:',/,
     & ' (1) Fit one data set only',/,
     & ' (2) Fit selected data sets separately',/,
     & ' (3) Fit selected data sets simultaneously with one equation',/,
     & ' (4) Fit selected data sets to estimate relative potencies',/,
     & ' (5) Fit selected data sets to estimate antagonist KB (Schild)',
     & /,'   (antagonist conc specified as ''set variable'')',/,
     & ' (6) Pool selected data sets and fit as one set',/,
     & '  Option number [',i2,'] = ')
c     & '',/,
	call INPUTi(idef)
	if(idef.ge.1.and.idef.le.6) then
	   ifitmode=idef
	else
	   goto 1201
	endif
	ifitsav=ifitmode		!keep initial ifitmode to write to .INI
c For refit after normalisation, better use the same sets as originally
c fitted so now jump tot 2021, not 202
202	continue	!jump here to refit norm data (with new selection of sets)
	do i=1,nset
	   jmiss(i)=1	!miss all
	enddo
	if(ifitmode.eq.1) then
	   jset=1
	   print 124
124	   format('&Fit set number [1] = ')
	   call INPUTi(jset)
	   jmiss(jset)=0	!fit set i
	   nmiss=nset-1
c	else if(ifitmode.ge.2.and.ifitmode.le.5) then
	else
	   if(.not.readini) nsfit=nset	!default=all
130	   print 126,nsfit
126	   format('&Number of sets to be fitted [',i3,'] = ')
	   call INPUTi(nsfit)
	   if(nsfit.lt.0.or.nsfit.gt.nset) goto 130
	   if(nsfit.lt.nset) then
		nmiss=nset-nsfit
		if(nmiss.gt.0) then
		 do 128 j=1,nsfit
		 print 129,j
129		 format('&  (',i2,')  fit set number = ')
	       call INPUTi(j1)
		 if(j1.lt.1.or.j1.gt.nset) goto 130
		 jmiss(j1)=0		!fit set j1
128		 continue
		endif
	   else		!fit all
		nmiss=0
		do 131 i=1,nset
131		jmiss(i)=0	!fit all
	   endif
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
c keep juse in ascending order
	call SORTI(juse,nsfit,.true.,njset)
c
2021	continue	!jump here to refit norm data, using same nsfit, juse etc
c  If several sets to be fitted as one, copy all data into xobs(i,1) etc.
c  and reset ifitmode=1, nsfit=1.  If set #1 is one of those to be fitted
c  data for other sets must be added to it, otherwise copy all data to
c  set 1
c ifitmode=6 to pool several sets into set 1 (sets
c  norm1=true so means plotted also in CVDISP
c -no option yet to use SD calculated here for weighting but new GETGROUP
c call below can be used to help with this
c Case of refit after normalising (idest=2).  Guesses are already in thetgues(i,j)
c for each of the sets (j) fitted separately before normalisation, but if
c refitted with ifitmode=6,pooled data are put into set 1 (for which a guess
c may not exist, and if it does, may not be appropriate) so in this case
c set guesdone=false here, so another guess made in cvsimp
c
	norm1=.false.
	if(ifitmode.eq.6) then
	   guesdone=.false.
	   use1=.false.
	   do m=1,nsfit	!check if set 1 is among those used
		j=juse(m)
		if(j.eq.1) use1=.true.
	   enddo
	   if(use1) then	!add rest of data to set 1
		n=nj(1)
		do m=1,nsfit
		   j=juse(m)
		   if(j.ne.1) then
			do i=1,nj(j)
			   n=n+1
			   Xobs(n,1)=Xobs(i,j)
			   Yobs(n,1)=Yobs(i,j)
			   w(n,1)=w(i,j)
			enddo
		   endif
		enddo
	   else		!copy all data to set 1
		n=0
		do m=1,nsfit
		   j=juse(m)
		   do i=1,nj(j)
			n=n+1
			Xobs(n,1)=Xobs(i,j)
			Yobs(n,1)=Yobs(i,j)
			w(n,1)=w(i,j)
		   enddo
		enddo
	   endif
         print 38
         if(pon()) write(7,38)
         if(discprt) write(8,38)
38	   format(' Data from following sets pooled into set 1')
         print 39,(juse(i),i=1,nsfit)
         if(pon()) write(7,39) (juse(i),i=1,nsfit)
         if(discprt) write(8,39) (juse(i),i=1,nsfit)
39	   format(30(2x,i3))
c      reset values
	   nj(1)=n
	   jset=1
	   nsfit=1
	   ifitmode=1
c  reset jmiss,juse
	   do i=1,nset
		jmiss(i)=1	!miss all
	   enddo
	   jmiss(jset)=0	!fit set i
	   nmiss=nset-1
	   juse(1)=1
	   norm1=.true.	!so means displayed in cvdisp
c       re-sort the arrays into order of ascending x values
	   call SORTr3(Xobs,Yobs,w,1,nj(1),.true.,
     &	   niobs,njset,niobs,njset)
	endif
c
c===
c Check for replicates at each x value and calc SS within groups
c etc if there are replicates
c Now do this befiore asking whether to normalise wrt Y at a specified x
c
c If ndf is big enough, record the error within groups and
c ask if it is to be used for errors -if so, can also
c calculate 'deviations from regression' and test its signif.
c No need to call getgroup in cvsimp (to check that there are enough
c x values to allow fit) if nxval is calculated here for all sets. The
c ybar, sdm also calculated here for each set fitted, and errors
c 'within x values' so calc could be used for error estimation
c -do this now for all sets, whatever fitmode!
97	continue	!jump here if nset=1
	if(ifitmode.ge.1.and.ifitmode.le.5) then
	   do m=1,nsfit	!check if set 1 is among those used
		j=juse(m)
	      call GETGROUP(xobs,yobs,w,j,nj,niobs,njset,nxv,
     &	   xgrp,ybar,sdy,sdm,ngrp,SSwg,ndf,SStot,SScol,ymean,n0)
		nxval(j)=nxv
		if(ndf.ge.1) then
		   Serrwg=sqrt(SSwg/float(ndf))
		   print 421,j
		   if(discprt) write(8,421) j
421		   format(1x,/,' Analysis of replicates for set ',i3)
		   print 42,nj(j),nxval(j),SStot,SScol,
     &		SSwg,ndf,Serrwg
		   if(discprt) write(8,42) nj(j),nxval(j),SStot,SScol,
     &		SSwg,ndf,Serrwg
42		   format(
     & 	   2x,i4,' observations at ',i4,' different X values',/,
     &	   ' Weighted sums of squares are:',/,
     &	   ' Total SSD = ',g13.6,/,
     &	   ' SSD between groups (X values) = ',g13.6,/,
     &	   ' Within group SSD = ',g12.5,' with ',i3,' df.',/,
     &	   ' s(err) = sqrt(SSD/df) = ',g13.6,/,/,
     &      ' Group  n(j)  X value      mean           sd          sdm')
		   do k=1,nxval(j)
			print 43,k,ngrp(k,j),xgrp(k,j),
     &		 ybar(k,j),sdy(k,j),sdm(k,j)
			if(discprt) write(8,43)k,ngrp(k,j),xgrp(k,j),
     &		 ybar(k,j),sdy(k,j),sdm(k,j)
43			format(i3,i8,4g13.6)
		   enddo
		   if(n0.gt.0) then
		   	print 8,n0
		   	if(discprt) write(8,8) n0
8			format(/,' ',i3,
     &		' values omitted because they have zero weight')
		   endif
		endif
	   enddo
	endif
c
c If more than one set to be fitted then check if there is an x value
c that is common to all sets -if so, ask whether to normalise w.r.t.
c y for this x, or whether numbers already normalised (if so, check for 1!)
c This can now be based on results from getgroup
c
	if(norm) goto 11		!already normalised to Y at given X
	if(norm1) goto 11		!already normalised after sep fits
c
	if(nsfit.gt.1) then
	   ALLOCATE(xobs1(niobs,nsfit),yobs1(niobs,nsfit),w1(niobs,nsfit),
     &	nj1(nsfit),xnorms(niobs))
	endif
c Go through each xvalue in the first set -then check if it occurs in all
c other sets -if it does then add it to xnorms() = list of possible x values
c that occur in all sets
	if(nsfit.gt.1) then
	   nnorm=0
	   norm=.true.
	   j1=juse(1)	!first set #

	   do k1=1,nxval(j1)	!go through all x values in set 1
		xv=xgrp(k1,j1)
		do m=2,nsfit	!is xv in all the other sets?
		   j=juse(m)
		   do k=1,nxval(j)
			if(EQUAL1(xv,xgrp(k,j))) then	!value same as xv found -can skip to next set
			   goto 70
			endif
		   enddo
c       if reach here then no xvalue found in set j that was same as xv
		   norm=.false.
70		   continue
		enddo	!end of loop for all sets but first
c If norm is still true here then the current xval is in all sets
		if(norm) then
		   nnorm=nnorm+1
		   xnorms(nnorm)=xv
		endif
	   enddo	!end of loop for all x values in set 1
	   if(nnorm.gt.1) then
		print 108
108		format(/)
		ans='N'
		call DCASK(
     &  'Normalise the observations relative to y at specifed x value',
     &	ans,ans)
		if(ans.eq.'N') then
		   DEALLOCATE(xobs1,yobs1,w1,nj1,xnorms)
		   norm=.false.
		   goto 11
		endif
c
452		continue
c     insert here calc of the # of replicates for each x value
		print 400
400		format(
     &' This option is available only for (a) Langmuir/Hill with one',/,
     &' component per curve, or power function fit at present (and  ',/,
     &' usually better to normalise with respect to fitted maximum  ',/,
     &' when possible) -proceed anyway [N] ? ')
		ans='N'
		call INPUTa(ans)
		if(ans.eq.'N') then
		   DEALLOCATE(xobs1,yobs1,w1,nj1,xnorms)
		   norm=.false.
		   goto 11
		endif
	      print 4
4	      format(
     &      ' The following X values occur in every set:')
		do n=1,nnorm
		   xv=xnorms(n)
		   nrep=0
		   do m=1,nsfit
			j=juse(m)
			do k=1,nxval(j)
			   if(EQUAL1(xv,xgrp(k,j))) nrep=nrep+ngrp(k,j)
			enddo
		   enddo
		   print 44,n,xnorms(n),nrep
44		   format(' (',i3,') X = ',g13.6,' (',i4,' replicates)')
		enddo
451		print 45,nnorm
45		format(' Number of the X value to be used (1 to ',i3,
     &		') = ')
		call INPUTi(in)
		if(in.lt.1.or.in.gt.nnorm) then
		   call BELL(1)
		   goto 451
		endif
		print 46,xnorms(in)
46		format(' Normalise with respect to response at X = ',
     &		g13.6,' [Y] ? ')
		ans='Y'
		call INPUTa(ans)
		if(ans.eq.'N') then
		   ans='Y'
		   call DCASK(
     &  'Cancel the normalisation',ans,ans)
		   if(ans.eq.'Y') then
			DEALLOCATE(xobs1,yobs1,w1,nj1,xnorms)
			norm=.false.
			goto 11
		   else
			goto 452
		   endif
		endif
		norm=.true.
		xnorm=xnorms(in)
	   endif	!end of nnorm.gt.1
	endif		!end of nsfit>1
	do m=1,nsfit
	   j=juse(m)
	   do k=1,nxval(j)
		if(EQUAL1(xgrp(k,j),xnorm)) then
		   ynorm(j)=ybar(k,j)	!normalising mean for set j
		endif
	   enddo
	enddo
	if(norm) then
	   do m=1,nsfit
		i1=0
		j=juse(m)
		do i=1,nj(j)	!go through every indiv observation
		   if(.not.EQUAL1(xobs(i,j),xnorm)) then	!skip points at x=xnorm
			i1=i1+1
			xobs1(i1,j)=xobs(i,j)
			yobs1(i1,j)=yobs(i,j)/ynorm(j)
			if(EQUAL1(1.0,w(i,j))) then
			   w1(i1,j)=1.0
			else
			   w1(i1,j)=w(i,j)/(ynorm(j)**2)	!as though ynorm=constant
			endif
		   endif
		enddo
		nj1(j)=i1
c At this stage still have nsfit separate sets, but data in them are now
c normalised.  In original version these sets were pooled and means put into
c set 2 (indiv values in set 1) -but may be better to keep them separate
c and keep open the options about how to fit them.  Now copy them into
c xobs(), yobs
		nj(j)=nj1(j)
		do i=1,nj(j)
		   xobs(i,j)=xobs1(i,j)
		   yobs(i,j)=yobs1(i,j)
		   w(i,j)=w1(i,j)
		enddo
c now have normalised observations on xobs, yobs, w -do getgroup again
		j=juse(m)
	      call GETGROUP(xobs,yobs,w,j,nj,niobs,njset,nxv,
     &	   xgrp,ybar,sdy,sdm,ngrp,SSwg,ndf,SStot,SScol,ymean,n0)
		nxval(j)=nxv
		if(ndf.ge.1) then
		   Serrwg=sqrt(SSwg/float(ndf))
		   print 422,j
		   if(discprt) write(8,422) j
422		   format(1x,/,
     &	    ' Analysis of replicates for normalised data: set ',i3)
		   print 42,nj(j),nxval(j),SStot,SScol,
     &	    SSwg,ndf,Serrwg
		   if(discprt) write(8,42) nj(j),nxval(j),SStot,SScol,
     &		  SSwg,ndf,Serrwg
		   do k=1,nxval(j)
			print 43,k,ngrp(k,j),xgrp(k,j),
     &		  ybar(k,j),sdy(k,j),sdm(k,j)
			if(discprt) write(8,43)k,ngrp(k,j),xgrp(k,j),
     &		    ybar(k,j),sdy(k,j),sdm(k,j)
c43			format(i3,i8,4g13.6)
		   enddo
		   if(n0.gt.0) then
		   	print 8,n0
		   	if(discprt) write(8,8) n0
c8			format(/,' ',i3,
c     &		' values omitted because they have zero weight')
		   endif
		endif
c===ask for ifitmode again?
	   enddo	!end of do m=1,nsfit
	   ans='Y'
	   call DCASK('View/alter the normalised data',ans,ans)
	   if(ans.eq.'Y') goto 71
	endif		!end of if(norm)
75	continue	!return here after redisplay of norm data
c NB later in prog, option 2 below is signalled by norm=T and ifitmode=1
	if(norm) then
	   print 61,nsfit
61		format(/,
     &' (1) Pool the ',i3,' sets into a single set (only option now)',/,
     &   ' Option number [1] = ')
c	   print 61,nsfit,nsfit
c61		format(/,
c     &   ' (1) Pool the ',i3,' sets into one single set',/,
c     &   ' (2) fit the ',i3,' normalised sets in the usual way',/,
c     &   ' Option number [1] = ')
	   iopt2=1
	   call INPUTi(iopt2)
	   if(iopt2.eq.2) then
1203		print 1202,idef
1202		format(/,' CHOOSE FITTING MODE:',/,
     & ' (2) Fit selected data sets separately',/,
     & ' (3) Fit selected data sets simultaneously with one equation',/,
     & ' (4) Fit selected data sets to estimate relative potencies',/,
     & ' (5) Fit selected data sets to estimate antagonist KB (Schild)',
     & /,'   (antagonist conc specified as ''set variable'')',/,
     & '  Option number [',i2,'] = ')
		call INPUTi(idef)
		if(idef.ge.2.and.idef.le.5) then
		   ifitmode=idef
		else
		   goto 1203
		endif
		ifitsav=ifitmode		!keep initial ifitmode to write to .INI
c Put the normalising point (alone) into an extra set for display
		ns1=nsfit+1
		nj(ns1)=1
		xobs(1,ns1)=xnorm
		yobs(1,ns1)=1.0
		w(1,ns1)=1.0
	   else	!iopt2=1: pool normalised data sets
c		First copy ALL obs into set 1 (values for x=xnorm already omitted)
		if(allocated(xobs1)) DEALLOCATE(xobs1,yobs1,w1)
		n=0
		do j=1,nsfit
		   n=n+nj(j)
		enddo
		ALLOCATE(xobs1(n,njset),yobs1(n,njset),w1(n,njset))
		i1=0
		do m=1,nsfit
		   j=juse(m)
		   do i=1,nj(j)
			i1=i1+1
			xobs1(i1,1)=xobs(i,j)
			yobs1(i1,1)=yobs(i,j)
			w1(i1,1)=w(i,j)
		   enddo
		enddo
		if(i1.gt.niobs) then
		   call BELL(2)
		   print 701,i1
701		   format(
     &' To use this option the ''maximum number of obs'''
     &' must be set to at least ',i5)
		   STOP
		endif
		nj(1)=i1
c        now sort the values, get means and put them in set 2
c Sort x values into ascending order
		call SORTr3(Xobs1,Yobs1,w1,1,nj(1),.true.,
     &	   n,njset,n,njset)
c copy data into Xobs(i,1)
		do i=1,nj(1)
		   xobs(i,1)=xobs1(i,1)
		   yobs(i,1)=yobs1(i,1)
		   w(i,1)=w1(i,1)
		enddo
c Calc second data set for means (with SD)
		j=1
	      call GETGROUP(xobs,yobs,w,j,nj,niobs,njset,nxv,
     &	   xgrp,ybar,sdy,sdm,ngrp,SSwg,ndf,SStot,SScol,ymean,n0)
c      copy means and SD into set 2, and print
		nxval(1)=nxv
		if(ndf.ge.1) then
		   Serrwg=sqrt(SSwg/float(ndf))
		   print 423
		   if(discprt) write(8,423)
423		   format(1x,/,
     & ' Analysis of replicates for pooled normalised data (now set 2)')
		   print 42,nj(1),nxval(1),SStot,SScol,
     &	    SSwg,ndf,Serrwg
		   if(discprt) write(8,42) nj(1),nxval(1),SStot,SScol,
     &		  SSwg,ndf,Serrwg
		   do k=1,nxval(1)
			xobs(k,2)=xgrp(k,1)
			yobs(k,2)=ybar(k,1)
			print 43,k,ngrp(k,1),xgrp(k,1),
     &		  ybar(k,1),sdy(k,1),sdm(k,1)
			if(discprt) write(8,43)k,ngrp(k,1),xgrp(k,1),
     &		    ybar(k,1),sdy(k,1),sdm(k,1)
c43			format(i3,i8,4g13.6)
		   enddo
		   if(n0.gt.0) then
		   	print 8,n0
		   	if(discprt) write(8,8) n0
c8			format(/,' ',i3,
c     &		' values omitted because they have zero weight')
		   endif
		endif
		nj(2)=nxval(1)
		do i=1,nj(2)
		   if(sdm(i,1).gt.1.e-10) then
			w(i,2)=1.0/sdm(i,1)**2
		   else
			w(i,2)=0.0
		   endif
		enddo
		ans='Y'
		call DCASK('View/alter the pooled normalised data',ans,ans)
		if(ans.eq.'Y') then
		   n1=2
		   n2=2
		   goto 711
		endif
	   endif	!temp end of if(iopt2.eq.1) then
	endif		!temp end of if(norm)
712	continue	!jump back here after redisplay of pooled norm data
	if(norm) then
	   if(iopt2.eq.1) then
		print 6
6		format(/,
     &   ' (1) weight the fit with SD calculated from replicates',/,
     &   ' (2) fit with equal weights',/,
     &   ' Option number [1] = ')
		i=1
		call INPUTi(i)
		if(i.eq.2) then
		   iw=1
		   do i=1,nj(2)
			w(i,2)=1.0
		   enddo
		else
c  weights set from sdm() above, but do not reassign here in case values
c  were changed during display in inwind
		   iw=2
c		   do i=1,nj(2)
c			if(sdm(i,1).gt.1.e-10) then
c			   w(i,2)=1.0/sdm(i,1)**2
c			else
c			   w(i,2)=0.0
c			endif
c		   enddo
		endif
		printed=.false.
	      print 7,xnorm
	      if(pon()) write(7,7) xnorm
	      if(discprt) write(8,7) xnorm
7		format(/,
     &      ' Fit to means of data normalised relative to X = ',
     &	f9.3,' (set 2)',/,
     &   '=========================================================',/)
		printed=.true.	!data for normalised fit printed above
c Put the normalising point (alone) into set 3
		nj(3)=1
		xobs(1,3)=xnorm
		yobs(1,3)=1.0
		w(1,3)=1.0
c   set up to fit set 2 only
		ifitmode=1	!if only one set
		nset=3
		jset=2	!fit set 2 =means
		jmiss(1)=1	!miss set 1
		jmiss(2)=0	!fit set 2
		jmiss(3)=1	!miss set 3	!normalisation point for display only
		nsfit=1
		juse(1)=2	!use set 2
		ncurvc=1
		icurvc(1)=2
		ncal(2)=ncalc
		icurvc(2)=1
		iline(2)=0	!continuous
	   endif	!end of if(i.eq.1) then
	endif		!end of if(norm)

	if(allocated(xobs1)) then
	   DEALLOCATE(xobs1,yobs1,w1,nj1,xnorms)
	endif
c Can calculate here values of ncurvc,icurvc etc for calc curves which
c are now calc by DEFYcal() during loop for each fitted set
11	continue
	if(ifitmode.eq.1) then
	   ncurvc=1
	   icurvc(1)=jset
	else
	   ncurvc=0
	   do 121 j=1,nset
	   if(jmiss(j).eq.1) goto 121	!set omitted
	   ncurvc=ncurvc+1
	   icurvc(ncurvc)=j
121	   continue
	endif
	do 122 j1=1,ncurvc
	 j=icurvc(j1)
c	 iline(j)=j1-1		!curve #1 continuous
	 iline(j)=0			!all curves continuous
	 if(ifitmode.eq.4.or.ifitmode.eq.5) iline(j)=0	!all continuous
	 ncal(j)=ncalc
122	continue
c	if(debug()) print 217,ncurvc,(icurvc(i),i=1,ncurvc)
c217	format(' ncurvc=',i4,'  icurvc= ',10i3)
c End of section for normalising wrt Y at specified X
c
c (jump to here if there is only one set)
98	continue
c
c
c If data not printed above then print it here for the fitted sets
	if(.not.printed) then
	   if(iver.ge.1002) then
	      print 322,titlef
	      if(pon()) write(7,322) titlef
	      if(discprt) write(8,322) titlef
c322	      format(/,' File: ',a60)
	   endif
	   do 132 j=1,nset
		if(jmiss(j).eq.1.or.shortprt) goto 132
		call PRINTJ(j,titled,setx,hdisp,setvar,prt,nj,w,xobs,yobs,
     &	 niobs,njset,norm1)
132	   continue
	endif
	if(iw.ge.1.and.iw.le.5) then
	   print 1301,titw(iw)
	   if(prt) write(7,1301) titw(iw)
         if(discprt) write(8,1301) titw(iw)
c1301	   format(/,' Weighting method = ',a60,/)
	endif
c
999	continue
	call flush(7)
	RETURN
	END

	logical function EQUAL1(x1,x2)
c (there is another function EQUAL in normdat)
	if(x1.gt.-1000. .and. x1.lt.1000. .and.
     &	x2.gt.-100. .and. x2.lt.100.) then
	   EQUAL1=abs(x1-x2).lt.1.e-5
	   RETURN
	else
	   if(x1.gt.1.e-30*x2) then
		x=abs(x2/x1) - 1.0
	   else
		x=abs(x1/x2) - 1.0
	   endif
	   EQUAL1=abs(x).lt.1.e-5
	endif
	RETURN
	end

