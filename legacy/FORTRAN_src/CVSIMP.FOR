	subroutine CVSIMP(kfit,jfix,nfix,titlep,fitted,dcurve,
     & nodata,theta,thetsav,nev,sres,smin,iflag,fixset,
     & Xobs,yobs,w,nj,juse,icuse,setx,niobs,njset,ndth,
     & badset,nset1,thetgues,xcal,ycal,ndc1,ndimc,ncalc,logx,logy,
     & logity,y0,yinf,xgt0,xmin,xmax,ymin,ymax,autosep,
     & iparsav,jfixsav,sminsav,neweq,guesdone,noguess,mset,ibk,
     & nxval,xgrp,ybar,sdm,ngrp,SSwg,ndf)
c To fit set #jset in CVFIT (NB jset now in common), or fits all
c sets simultaneously if ifitmode=3,4,5 (see SSCDV)
c In call nset1 is copy of nset (in common) to allow declaration of
c thetsav, thetgues
c
c IFLAG=2 when CVSIMP is doing sep fits for indiv curves (ifitmode=2) as
c the last section of ifitmode=4,5 calcs, after simultaneous fits
c IFLAG=3,4 when separate fits followed by normalistaion and fitmode
c other than 2
c
c Modif 03/27/95 11:34am for allocatable xobs() etc, so all data must
c be passed as params to SIMPLEX and SSDSV
c	function SSDCV(K,THETA,Xobs,yobs,w,nj,juse,setx,niobs,njset)
c
c Modif 09/18/94 06:15pm for VGA mode input
	integer*4 videotyp
	character*20 TITLEP(ndth)		!names of params for fit
	real*4 theta(ndth)
	real*4 thetsav(100),thetgues(100,nset1)
	integer JFIX(ndth),jfixsav(ndth),iparsav(10)
	logical discprt,alpha,skipymax,nodata,xgt0
c	logical fixmax,fix0
	LOGICAL constr,fline,pon,fitted,dcurve,slock,negpar,fixratio
	logical badset,skipq,hillang,autosep,neweq,sepmax,guesdone
	character*1 ans,UC,ans1
	character*11 cnum,cnum1,cnum2,cnum3
	character defolt*30,cans*30		!to hold text & result of QDIALOG
	character*25 text
	logical help
c Params to be passed to SSDCV and Simplex
	real*4 xobs(niobs,njset),yobs(niobs,njset),w(niobs,njset)
	integer juse(njset),nj(njset)
	integer icuse(100) 	!see geteqn, cvsimp
c For guesses etc
	real*4 xcal(ndc1,ndimc),ycal(ndc1,ndimc)
	real setx(njset)
	logical logx,logy,logity,noguess(100),fixset
c
	allocatable::ymaxsep
	real*4 ymaxsep(:)
c
c for getgroup
	real*4 xgrp(niobs,njset),ybar(niobs,njset),sdm(niobs,njset)
	integer ngrp(niobs,njset),nxval(njset)
c
	common/hlp/help				!for QDIALOG
c====nj(),jmiss(),juse() no longer in common/block1/
	COMMON/BLOCK1/constr,nset,nfit,nsfit,
     & Xv,alfa,kmax,ncomp,nmod,fline,nomit,jomit(20),jset,ifitmode
	logical logyfit,norm
	COMMON/BLOCK2/ castar,X1,X2,iequiv,ip1,ip2
	COMMON/BLOCK3/logyfit,norm,xnorm,iostat		!CVIN2,YCALCV,CVDISP,GETEQN,CVSIMP
	common/fix/fixratio,ir1,ir2,rval	!to fix ratio of 2 parameters
	common/potrat/jfirst,iset,kmax1			!for pot ratios
	common/pwrfunc/ybarp
	common/graf/alpha		!for simplexv
	common/dp/discprt
	external SSDCV
c
c
	pon()=slock()
c
101	format(a1)
104	FORMAT(I8)
	alpha=VIDEOTYP().eq.3		!alpha mode
c
	errfac=1.e-4
	delmin=-1.		!do not use delmin for convergence
	stpfac=0.1
	skipq=.false.
c For ifitmode=2, mset>1 keep jfix as for first set (mset=1)
c================
c Also, when (a) doing sep fits (mode 2) after mode 4,5
c and (b) doing another fit after sep fits, may want to set jfix
c elsewhere so not zeroed here
	if(ifitmode.eq.2.and.mset.gt.1) goto 2	!keep jfix()
	if(ifitmode.eq.2.and.iflag.eq.2) goto 2	!keep jfix() after mode 4,5
	if(fixset) goto 2					!keep jfix(); already set
	 nfix=0		!no param fixed
	 do i=1,100
	    jfix(i)=0
	 enddo
	 fixset=.false.
2	continue
c
	print 108
108	format(/)
	negpar=.false.
c Set default initial guesses for some models (for ifitype=4,5 set
c these for the initial simultaneous fit here)
c For nmod=28 calculate init guesses for each set as in pwrfunc.mcd
c and keep them for later in thetgues(i,j).
c If more than one set is being fitted then init guesses for all sets
c (separately) are found here and kept in thetgues(), so do only when
c iflag=0
c Skip unless guess available for nmod!
	if(norm.and.nmod.ne.28) goto 102	!guesses not fixed for Hill yet
	if((iabs(nmod).eq.1
     & .or.(nmod.eq.26.and.(.not.fline))
     & .or.(nmod.eq.27.and.(.not.fline))
     & .or.nmod.eq.28
     & .or.(nmod.eq.99.and.ip1.eq.1))
     & .and.iflag.eq.0.and.ncomp.eq.1.and.
     & (.not.(nodata.or.dcurve))) then
	   ALLOCATE(ymaxsep(nset))	!for guesses for maxima (nmod=27)
	   call GETGUESS(nmod,Xobs,Yobs,nj,juse,niobs,njset,nsfit,
     &    theta,ndth,thetgues,nset1,kmax,ybarp,ncomp,ifitmode,jset,
     &    setx,jfix,nfix,logyfit,norm,sepmax,ymaxsep,guesdone,noguess)
c do fit without asking?  Or use theta values as defaults? -better show
c the guesses on graph and ask if OK to use.  NB for ifitmode=4,5 one
c theta() gives params for all sets, but for ifitmode=2,3 need to load
c appropriate theta() for each set
	   iset=0
	   do m=1,nsfit
		jval=juse(m)
		iset=iset+1	!for common -needed for ifitmode=4,5 if sets omitted
		jset=jval		!for COMMON
		if(ifitmode.le.2) then
		   do i=1,kmax
			theta(i)=thetgues(i,jset)	!params for current set
		   enddo
		endif
c=======	noguess(4)=.true.	!for debug
		if(.not.noguess(jval)) then
		   call DEFYcal(xmin,xmax,ncalc,logx,theta,Xcal,Ycal,
     &		jval,xgt0,setx,njset,ndimc,ndc1)
		endif
	   enddo
	   call DRAWGUES(xmin,xmax,ncalc,logx,logy,logity,y0,yinf,
     &	 Xcal,Ycal,ymin,ymax,nsfit,juse,noguess,njset,ndc1,ndimc,14)
c       For sep fits restore guesses for the first set (after drawing guesses)
	   if(ifitmode.le.2) then
		jset=juse(1)
		do i=1,kmax
		   theta(i)=thetgues(i,jset)	!params for current set
		enddo
	   endif
	   ans='Y'
	   call DEFOLTa(ans,defolt)
	   call QDIALOG(1,
     &	'Use these yellow curves as initial guess for fit',
     &	  defolt,11,cans)
	   call GETINPa(cans,ans)
c (for ifitmode=4,5, asked in cvfit whether to do sep fits with defaults)
	   if(ans.eq.'Y') then
		autosep=.true.
	   else	!guesses not OK
		if(alpha) then
		   call DCASK('Try a different equation','y',ans)
		else
		    ans='Y'
		    call DEFOLTa(ans,defolt)
		    call QDIALOG(1,
     &	      'Try a different equation',
     &	        defolt,11,cans)
	          call GETINPa(cans,ans)
		endif
		if(ans.eq.'Y') then
		   call DRAWGUES(xmin,xmax,ncalc,logx,logy,logity,y0,yinf,
     &	    Xcal,Ycal,ymin,ymax,nsfit,juse,noguess,njset,ndc1,ndimc,ibk)  !delete
		   neweq=.true.
		   RETURN
		else
		   autosep=.false.
		endif
	   endif
c Delete guesses
	   call DRAWGUES(xmin,xmax,ncalc,logx,logy,logity,y0,yinf,
     &    Xcal,Ycal,ymin,ymax,nsfit,juse,noguess,njset,ndc1,ndimc,ibk)	!delete
c  Meed to weight sets differentially if sep max fitted and their
c  maxima are very different
	   if(ifitmode.eq.3.and.nmod.eq.27.and.ip1.eq.1.and.sepmax) then
		if(alpha) then
		   print 109
109		   format(' Maxima are very different',/,
     &	' -suggest weighting sets with (guessed) Ymax: O.K. [Y] ? ')
		   ans='Y'
		   call INPUTa(ans)
		else
		   call WDIALOG(1,'Maxima are very different',11)
		   ans='Y'
		   call DEFOLTa(ans,defolt)
		   call QDIALOG(1,
     &	   '-suggest weighting sets with (guessed) Ymax: O.K.',
     &	   defolt,11,cans)
	         call GETINPa(cans,ans)
		endif
		if(ans.eq.'Y') then
		   do m=1,nsfit
			j=juse(m)
			ym=ymaxsep(j)
			do i=1,nj(j)
			   w(i,j)=w(i,j)/ym**2
			enddo
		   enddo
		   if(pon()) write(7,110)
      	   if(discprt) write(8,110)
110            format(' Each set weighted with guessed Ymax to',/,
     &	     ' compensate for different maxima',/)
		endif
	   endif
	   DEALLOCATE(ymaxsep)
c If autosep then ask which params to be fixed here (and set fixset=true
c so it is not done again) -otherwise ask which to be fixed in normal
c place (which is skipped when autosep true)
	   if(autosep) then
		call SETFIX(nfix,jfix,kmax,kfit,theta,titlep,alpha,ndth)
		fixset=.true.
		fitted=.true.
		dcurve=.false.
		ndisp=20
		skipq=.true.
		negpar=.false.
		do i=1,kmax
		   if(theta(i).lt.0.0) negpar=.true.
		enddo
		if(negpar.or.nmod.eq.21) then
		   constr=.false.
		   ans='N'
		else
		   constr=.true.
		   ans='Y'
		endif
		goto 524				!and skip questions
	   else
		fixset=.false.
	   endif
	endif		!end of section for initial guesses
102	continue
c
c Deal with case where separate fits of Langmuir/Hill are being done
c following simultaneous fit with ifitmode=4,5.  In these cases
c want values that were fixed for simult fit to be still fixed for
c separate fits
c First the case where separate fits have maximum FIXED at the value
c previously found for simultaneous fit
	hillang=nmod.eq.-1.or.nmod.eq.1.or.nmod.eq.26.or.nmod.eq.27
c	if(hillang.and.ifitmode.eq.2) then
c	   if(fixmax) then
c	      jfix(2)=1		!Ymax
c	      nfix=nfix+1
c	   else if(hillang.and.fix0) then	!fix Y(0) also
c		jfix(1)=1
c	      nfix=nfix+1
c	   endif
c	endif
c	if(iflag.eq.2.and.hillang) then
c	 if(isepfit.eq.1) then
c	   if(fixmax) then
c	      jfix(2)=1		!Ymax
c	      nfix=nfix+1
c	   else if(hillang.and.fix0) then	!fix Y(0) also
c		jfix(1)=1
c	      nfix=nfix+1
c	   endif
c	 endif
c	 if(isepfit.eq.2) then
c	   if(alpha) then
c	      print 41,thetsav(2)
c41	      format(
c     &    ' Maximum fixed at',g13.6,' (from simult fit): O.K. [Y] ')
c	      read 101,ans
c	   else
c		ans='Y'
c		call DEFOLTa(ans,defolt)
c	      call REALTOCH(thetsav(2),cnum,11)
c		n=nblank(cnum)
c		call QDIALOG(1,
c     &	'Maximum fixed at '//cnum(1:n)//' (from simult fit): O.K.',
c     &	  defolt,11,cans)
c		call GETINPa(cans,ans)
c	   endif
c	   if(UC(ans).ne.'N') then
c	      theta(2)=thetsav(2)	!ymax
c	      jfix(2)=1
c	      nfix=nfix+1
c	   endif
c	   if(hillang.and.fix0) then	!fix Y(0) also
c		jfix(1)=1
c	      nfix=nfix+1
c	   endif
c	 endif
c	endif
c
c    In following cases init guesses already set in theta (in main prog)
c -fix variables where nec and skip to fit (unless noguess=true for
c current set) (for ifitmode=2 could more simply just leave jfix() as
c it is for first fit, rather than zeroing it at each entry to cvsimp, at least
c as long as same eqn fitted to each curve -if eqn changed need
c to set noguess=true!)
	if(nmod.le.4.or.hillang.or.nmod.eq.99.or.
     &  nmod.eq.23.or.nmod.eq.24.or.nmod.eq.25.or.nmod.eq.28) then
	   j=jset
c
c nmod=23,24 added = Langmuir/Hill (resp) with common
c maximum when more than one set being fitted.  These are defined in ycalcv
c and nmod set to 23,24 in neweqcv (not geteqn); used only if ifitmode=4,5
c and isepfit=3
c nmod=25 =straight lines with common intercept on Y axis
c nmod=28 (imod=23) for power function
c  NB kmax reset to kmax1 (no extra params) in CVFIT after setting iflag=2
c for the separate fits
c	   if(iflag.eq.2.and.(nmod.eq.23.or.nmod.eq.24)) then
c		j=1		!index for thetsav() (see NEWEQCV)
c		if(fixmax) then
c		   jfix(2)=1            !max was fixed for simult fit, so fix now
c	         nfix=nfix+1
c		endif
c		if(fix0.and.nmod.eq.24) then
c		   jfix(1)=1            !Y(0) was fixed for simult fit, so fix now
c	         nfix=nfix+1
c		endif
c	   endif
c don't skip guesses if noguess=true for current set
	   if(autosep.and.(.not.noguess(j))) then
		fitted=.true.
		dcurve=.false.
		ndisp=20
		skipq=.true.
		negpar=.false.
		do i=1,kmax
		   if(theta(i).lt.0.0) negpar=.true.
		enddo
		if(negpar.or.nmod.eq.21) then
		   constr=.false.
		   ans='N'
		else
		   constr=.true.
		   ans='Y'
		endif
		goto 524				!and skip questions
	   endif
	endif
c
c Get the parameter guesses next
	text='Specify initial guesses:'
	dcurve=.false.
	if(nodata) then
	   text='Specify parameter values:'
	   do i=1,ndth
		jfix(i)=1		!for printout in cvout1
	   enddo
	else
	   if(alpha) then
		print 200
200		format(
     &     '&Use parameter values to plot curve without fitting [N] ? ')
		read 101,ans
	   else
		ans='N'
		call DEFOLTa(ans,defolt)
		call QDIALOG(1,
     &	'Use parameters to plot curve without fitting',
     &	defolt,11,cans)
	      call GETINPa(cans,ans)
	   endif
	   if(UC(ans).eq.'Y') then
		fitted=.false.
		dcurve=.true.
		do i=1,ndth
		   jfix(i)=1		!for printout in cvout1
		enddo
		text='Specify parameter values:'
	   endif
	endif
124	continue
c	skipymax=.false.
	if(alpha) then
	   print 1081,text
1081	   format('&',a25,/)
	else
	   call WDIALOG(1,text,11)
	endif
c Set numbering versus curve numbering -would be helpful to specify
c not only the set number but also the curve number (corresponding
c to key on the graph) when asking for parameter values).
c  For ifitmode=2 we fit only one set at a time (#jset) and
c the curve number is already available in mset, but set number
c that is incorporated in parameter names (in titlep()) needs to be
c changed for each set
c  For ifitmode=2,4,5 do only one simultaneous fit, and icuse(k) is the curve
c number that corresponds to set number that is already incorporated
c into titlep(k) in character format (in geteqn).  If icuse()=0 then
c the parameter is a common value that does not refer to
c a particular set
	do 131 i=1,kmax
c For ifitmode=2 must change set number in titlep here -it has already been
c put into titlep() -bit clumsy! -so need to search for it!
	  if(ifitmode.eq.2) then
	   do m=1,18
	    if(titlep(i)(m:m+2).eq.'set') then
		call INTCONV(jset,cnum)
		n=nblank(cnum)
		n2=m+3+n
		if(n2.gt.20) n2=20
		titlep(i)(m:n2)='set '//cnum(1:n)
		goto 421
	    endif
	   enddo
	  endif
421	  continue
	  i1=i
	  if(NMOD.eq.7.and.i.gt.3) i1=3
1321	  continue
c Now ask for guesses!
	  if(icuse(i).ne.0) then
	     call INTCONV(icuse(i),cnum1)
	  endif
71	  if(alpha) then
	     if(icuse(i).eq.0) then
	       print 132,i,TITLEP(i1),theta(i)
132	       format('&Parameter ',i2,': ',a20,' [',g10.3,'] = ')
	       call INPUTr(theta(i))
	     else
	       print 1132,i,TITLEP(i1),icuse(i),theta(i)
1132	       format('&Parameter ',i2,': ',a20,
     &	' (curve ',i2,') [',g10.3,'] = ')
	       call INPUTr(theta(i))
	     endif
	  else
	     call INTCONV(i,cnum)
	     call DEFOLTr(theta(i),defolt)
	     if(icuse(i).eq.0) then
	        call QDIALOG(1,'Param('//CHARNB(cnum)//
     &         '): '//CHARNB(titlep(i1)),defolt,11,cans)
	     else
	        call QDIALOG(1,'Param('//CHARNB(cnum)//
     &         '): '//CHARNB(titlep(i1))//
     &	' (curve '//CHARNB(cnum1)//')',defolt,11,cans)
	     endif
	     call GETINPr(cans,theta(i))
	  endif
	  if(ifitmode.eq.3.and.ip1.eq.2.and.i.eq.kmax) then  !common K, sep Ymax for each set
		if((theta(i).lt.0.).or.theta(i).gt.1.) then
		   call BELL(2)
		   if(alpha) then
			print 70
70		 	format(
     &    ' Fraction of component 1 must be between 0 and 1: try again')
		   else
			call WDIALOG(1,
     &  'Fraction of component 1 must be between 0 and 1: try again',12)
		   endif
		   goto 71
		endif
	  endif
c
	  if(theta(i).lt.0.) negpar=.true.
	  if(nodata.or.dcurve) goto 131
c=	  if(abs(theta(i)).lt.1.e-10) then
c=	    call BELL(2)
	  if(jfix(i).eq.1) then
		if(alpha) then
		   print 1331,i,theta(i)
1331		   format('&Parameter #',i2,' already fixed at ',
     &		g13.6,': O.K. [Y] ? ')
		   ans='Y'
		   call INPUTa(ans)
		else
		   ans='Y'
		   call DEFOLTa(ans,defolt)
	         call INTCONV(i,cnum)
	   	   call REALTOCH(theta(i),cnum1,11)
		   n=nblank(cnum)
	         call QDIALOG(1,'Parameter #'//CHARNB(cnum)//
     &     ' already fixed at '//CHARNB(cnum1)//': O.K.',defolt,12,cans)
	         call GETINPa(cans,ans)
		endif
		if(UC(ans).eq.'N') then
		   goto 1321		!new guess
		endif
	  endif
	  if(abs(theta(i)).lt.1.e-10.and.jfix(i).ne.1) then
	      call BELL(2)
		if(alpha) then
		   print 133,i
133		   format(
     &'&Initial guess of zero not allowed unless you intend to fix',
     &	' this parameter',/,
     &	' Fix parameter number ',i2,' at zero [Y] ? ')
		   ans='Y'
		   call INPUTa(ans)
		else
		   call WDIALOG(1,
     & 	'Initial guess of zero not allowed unless you intend',12)
		   call WDIALOG(1,
     & 	' to fix this parameter:                            ',12)
		   ans='Y'
		   call DEFOLTa(ans,defolt)
	         call INTCONV(i,cnum)
	         call QDIALOG(1,'Fix parameter #'//CHARNB(cnum)//
     &      ' at zero',defolt,12,cans)
	         call GETINPa(cans,ans)
		endif
		if(UC(ans).eq.'N') then
		   goto 1321		!new guess
		else
		   theta(i)=0.0
		   if(jfix(i).ne.1) then	!don't increment nfix if already fixed
			nfix=nfix+1
			jfix(i)=1
		   endif
		endif
	 endif
131	continue
c
	fitted=.true.
	if(dcurve.or.nodata) then
	   fitted=.false.
c Save parameters
	   j=jset
	   if(ifitmode.eq.4.or.ifitmode.eq.5) then	!simultaneous fit
	     do i=1,kmax
		thetsav(i)=theta(i)
	     enddo
	   endif
	   RETURN
	endif
	dcurve=.false.
	if(.not.fixset) then
	   call SETFIX(nfix,jfix,kmax,kfit,theta,titlep,alpha,ndth)
	   fixset=.true.
	endif
C
c  Check initial guesses: if any are neg make default to NOT constrain
c to be positive
	if(negpar.or.nmod.eq.21) then
	   constr=.false.
	   ans='N'
	else
	   constr=.true.
	   ans='Y'
	endif
c
	fixratio=.false.
	if(alpha) then
	   print 134,ans
134	   format('&Constrain parameters to be positive [',a1,'] ? ')
	   read 101,ans1
	   if(UC(ans1).eq.'N'.or.UC(ans1).eq.'Y') ans=ans1
	   constr=ans.eq.'Y'
	   print 135
135	   format('&Fix the ratio of two parameters [N] ? ')
	   read 101,ans
	   if(UC(ans).eq.'Y') then
	      fixratio=.true.
	      print 136
136	      format('&Fix the value of param(i)/param(j):  i,j = ')
		call INPUT2i(ir1,ir2)
	      print 138
138	      format('&Value of param(i)/param(j) required = ')
	      call INPUTr(rval)
	   endif
	else
	   call DEFOLTa(ans,defolt)
	   call QDIALOG(1,'Constrain parameters to be positive',
     &	defolt,11,cans)
	   call GETINPa(cans,ans)
	   constr=ans.eq.'Y'
	   if(nmod.eq.31) then	!MWC
		call DEFOLTa(ans,defolt)
		call QDIALOG(1,'Constrain R to be >1 (i.e. agonist)',
     &	defolt,11,cans)
		call GETINPa(cans,ans)
		fline=ans.eq.'Y'
	   endif
c
	   ans='N'
	   call DEFOLTa(ans,defolt)
	   call QDIALOG(1,'Fix the ratio of two parameters',
     &	defolt,11,cans)
	   call GETINPa(cans,ans)
	   if(UC(ans).eq.'Y') then
	      fixratio=.true.
	     call QDIALOG(1,
     &  'Value of param(i)/param(j) required',' ',11,cans)
	     call GETINPr(cans,rval)
	   endif
	endif
	if(fixratio) then
c Fix parameter #ir1 as its value is set from rval in SSDCV
	   nfix=nfix+1
	   kfit=kfit-1
	   jfix(ir1)=1
	   if(pon()) write(7,139) ir1,ir2,rval
         if(discprt) write(8,139) ir1,ir2,rval
139	   format(
     & ' Ratio of parameter ',i3,' / ',i3,' constrained to ',g13.6)
	endif
	ndisp=20
c
	if(alpha) then
	   call DCASK('Are these guesses O.K.','y',ans)
	else
	   ans='Y'
	   call DEFOLTa(ans,defolt)
	   call QDIALOG(1,'Are these guesses O.K.',
     &	defolt,11,cans)
	   call GETINPa(cans,ans)
	endif
	if(ans.eq.'N') then
	   fixset=.false.
	   goto 124
	endif
c
524	continue
c Check whether fit is possible -if number of different x values
c is less than number of free parameters then fit cannot be done
c (at least when sep fits are being done)
c Now nxval(j) is calculated already in cvin2 for all sets to be fitted
c At the moment getgroup is called only to get nxval, but if there are
c replicate Y at given X, the 'within x' errors calculated here
c could be used for error calcs
	badset=.false.
	if(ifitmode.le.2) then
c
	   kfit=kmax-nfix    !in case jump to 524
	   if(alpha) then
	      print 45,kfit,nxval(jset)
	   else
		call INTCONV(kfit,cnum)
		call INTCONV(nxval(jset),cnum1)
		call WDIALOG(1,
     &'No. of free parameters = '//charnb(cnum)//
     &'  No. of different x values = '//charnb(cnum1),11)
	   endif
	   if(pon()) write(7,45) kfit,nxval(jset)
         if(discprt) write(8,45) kfit,nxval(jset)
45	   format(
     &	' Number of free parameters = ',i3,/,
     &	' Number of different x values = ',i4)
	   if(nxval(jset).lt.kfit) then
		call BELL(1)
		if(alpha) then
	   	   print 46,jset,kfit,nxval(jset)
		else
		   call INTCONV(jset,cnum)
		   call WDIALOG(1,
     &	'Not enough points for fit to set '//charnb(cnum)//
     &	' (hit any key)',12)
		   call ANYKEY()
		endif
	      if(pon()) write(7,46) jset,kfit,nxval(jset)
            if(discprt) write(8,46) jset,kfit,nxval(jset)
46	      format(' Not enough points for fit to set ',i3,/,
     &	' Number of free parameters = ',i3,/,
     &	' Number of different x values = ',i4)
		badset=.true.
		RETURN
	   endif
	endif
c
	irestrt=3
	errfac=1.e-4
	ndisp=20
	if(.not.alpha) ndisp=0
	if(skipq) goto 5241
	if(alpha) then
	   print 312,irestrt
312	   FORMAT(' Limit number of restarts to [',i2,']: n = ')
	   call INPUTi(irestrt)
  	   print 525
525	   format('&Relative error [1.e-4]= ')
	   call INPUTr(errfac)
	   print 523
523	   format(
     & ' Print every Nth estimate [20] :  N= ')
	   call INPUTi(ndisp)
	else
	   call DEFOLTi(irestrt,defolt)
	   call QDIALOG(1,
     &      'Limit number of restarts to',defolt,11,cans)
	   call GETINPi(cans,irestrt)
	   call DEFOLTr(errfac,defolt)
	   call QDIALOG(1,
     &      'Relative error',defolt,11,cans)
	   call GETINPr(cans,errfac)
	   call DEFOLTi(ndisp,defolt)
	   call QDIALOG(1,
     &      'Print every Nth estimate: N ',defolt,11,cans)
	   call GETINPi(cans,ndisp)
	endif
c
5241	continue
	kfit=kmax-nfix    !in case jump to 524
	nevm=-200000		!so iterations not printed to disc
	confac=0.5		!parameter for simplex3
	smin=SSDCV(kmax,theta,Xobs,yobs,w,nj,juse,setx,niobs,njset)
c print initial guesses
	if(pon()) write(7,20) smin,(theta(i),i=1,kmax)
      if(discprt) write(8,20) smin,(theta(i),i=1,kmax)
20	format(/,' Initial guesses: SSD = ',g12.5,/,
     &  10g10.3,/,10g10.3,/,10g10.3,/,10g10.3,/,10g10.3,/)
	if(alpha) then
	  if(kmax.ge.7) print 5231,1,smin,(theta(I),I=1,KMAX)
	  if(kmax.lt.7) print 5232,1,smin,(theta(I),I=1,KMAX)
5231	  format(i8,g10.3,/,10g10.3,/,10g10.3,/,10g10.3,/,10g10.3,/)
5232	  format(i6,7g11.4,/)
	endif
c
	call SIMPLEXv(kmax,THETA,stpfac,errfac,nev,nevm,
     & smin,SSDCV,Ndisp,jfix,delmin,confac,irestrt,iconv,
     & Xobs,yobs,w,nj,juse,setx,niobs,njset,ndth)
C 100=dimension of theta(), jfix().
	sres=-1.
	if(nfit.gt.kfit) then
	   var=smin/float(nfit-kfit)		!error var
	   if(var.gt.1.e-30) sres=sqrt(var)
	else if(kfit.eq.nfit) then
	   var=0.0
	   sres=0.0
	endif
	if(alpha) then
	   print 5233
5233	   FORMAT( ' END OF FITTING')
	   IF(KMAX.GE.7) print 5231,NEV,SMIN,(THETA(I),I=1,KMAX)
	   IF(KMAX.LT.7) print 5232,NEV,SMIN,(THETA(I),I=1,KMAX)
	endif
	ndf=nfit-kfit		!degrees of freedom for error
	if(alpha) then
	   print 5234,nfit,kfit,sres,ndf
	else
	   call INTCONV(nfit,cnum)
	   call INTCONV(kfit,cnum1)
	   call REALTOCH(sres,cnum2,11)
	   call INTCONV(ndf,cnum3)
	   call WDIALOG(1,charnb(cnum)//' points fitted: '//
     &   charnb(cnum1)//' parameters estimated',11)
	   call WDIALOG(1,'Error S.D. ='//charnb(cnum2)//
     &   ' ('//charnb(cnum3)//' degrees of freedom)',11)
	endif
	if(pon()) write(7,5234) nfit,kfit,Smin,sres,ndf
      if(discprt) write(8,5234) nfit,kfit,Smin,sres,ndf
5234	format(/,' Number of points fitted = ',i8,/,
     & ' Number of parameters estimated = ',i8,/,
     & ' Minimum SSD (Smin) = ',g13.6,/,
     & ' Error S.D.= sqrt(Smin/df) = ',g13.6,
     &'  (',i5,' degrees of freedom)',/)
c
c Save parameters in case several fits done
c Need to keep values only when a simultaneous fit is to be followed
c by separate fits, when values from the simultaneous fit may be needed
c (eg to fix maximum for sep fits), and to construct guesses for the sep
c fits).  Keep result of simultaneous fit in thetsav(i,j).  Guesses
c for the sep fits are already in thetgues(i,j), but they could be improved
c on the basis of the simultaneous fit just completed.
	j=jset
	if(ifitmode.eq.4.or.ifitmode.eq.5) then
	   do i=1,kmax
		thetsav(i)=theta(i)
		jfixsav(i)=jfix(i)
	   enddo
	   sminsav=smin
	   iparsav(1)=ifitmode 	!save for lik intervals on sim fit
	   iparsav(2)=nmod
	   iparsav(3)=kmax
	   iparsav(4)=kmax1
	   iparsav(5)=kfit
	   iparsav(6)=nfit
	   iparsav(7)=ip1
	   iparsav(8)=ip2
c
c Revise the guesses for the separate fits using the results of sim fit
c       For power function, indiv slopes should be good, but adjust x0
	   if(nmod.eq.28) then
		thetgues(1,1)=theta(1)		!x0 for 1st set
		do j=2,nsfit
	         thetgues(1,j)=thetgues(1,1)*theta(j+1)
		enddo
	   endif
	endif
c After separate fits (ifitmode=2) keep replace guesses in thetgues()
c by the fitted params for each set, for use if sets are normalised w.r.t.
c their fitted maxima after the sep fits (see cvdisp)
	if((nmod.eq.26.or.nmod.eq.27).and.
     &	ifitmode.eq.2.and.iflag.ne.3) then
	   do i=1,kmax
		thetgues(i,jset)=theta(i)
	   enddo
	endif
c
	RETURN
	end


