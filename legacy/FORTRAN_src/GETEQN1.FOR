	subroutine GETEQN(IW,titw,titlep,xgt0,iflag,sameq,
     & jmiss,juse,iuse,icuse,njset,ndth,imod,nodata,iesc)
c
c To define which equation to be fitted in CVFIT2
c
c NB nmod=1, -1 exist only for ifitmode=4,5; if ifitmode.le.3 then they no longer
c    (except transiently, here) -in this case Langmuir/Hill are now nmod=26,27
c     12/27/96 05:36pm nmod=26,27 when ifitmode=4,5 now
c
c Modif 02/22/98 03:23pm for Hill/Langmuir, ifitmode=3, now has
c ip1=-1 -separate max, separate K
c Modif 02/01/98 03:45pm to define new array iuse(1:kmax) which defines
c for each parameter the curve number (1,2,..., nsfit) to which
c the parameter applies, with iuse(i)=0 if theta(i) is not specific
c for a particular data set.  This array allows curve number (as well as
c set number) to be shown when asking for initial guesses in CVSIMP
c
c Modif 02/01/98 07:11am New option for Hill/Langmuir (ip1=2 for this option)
c Only for ifitmode=3, ncomp=2 -common Ymax1/ymax2 -in this case fit K,
c    Ymax (and nH if req), one value for each set (not separate Ymax1, Ymax2)
c    plus a single value of Ymax1/Ymax for all sets
c Fit percent comp 1, say P1 = Ymax1/(Ymax1+Ymax2) = Ymax1/Ymax.  Then to
c get Ymax1 and Ymax2 (eg in ycalcv) for each set, find
c Ymax1 = P1*Ymax, then Ymax2=Ymax-Ymax1
c
c Modif 04/18/97 06:50am add nmod=31 for MWC(n)
c Modif 04/12/97 add nmod=30 for n bindings + opening
c Modif 03/20/96 10:01am Add nmod=28 (imod=23) for power function
c Modif 09/18/94 06:08pm for graphics mode input (not finished yet!)
c Modif 09/17/94 11:38am Add nmod=25 =straight lines with common
c   intercept on Y axis
c 04/13/92 07:08pm nmod=23,24 added = Langmuir/Hill (resp) with common
c   maximum when more than one set being fitted.  These are defined in ycalcv
c   and nmod set to 23,24 in neweqcv (not here); used only if ifitmode=4,5
c   and isepfit=3.
c 04/11/92 06:12pm ncomp now in common
c
	character*20 TITLEP(ndth)     !names of params for fit
	character*60 titmod(30)		!names of models
	character*60 titw(5)		!names of weighting methods
	character*10 T1(50)		!'standard' param names
	character*11 charstr
c	character*4 tt1(10)
	integer*4 videotyp
C THESE ARE ARRAYS TO HOLD TITLES (UP TO 60 CHAR) FOR DATA SETS
C (IN TITLED) AND FOR PARAMETERS (TITLEP)(10 CHAR)
c Also include option to put axis titles in here, and store them with the
c data
	character*1 ans,UC
	LOGICAL CONSTR,FLINE,alpha
	logical discprt,slock,pon,pprt,xgt0,titprint,sameq
	character cnum*11,cnum1*11
	character defolt*30,cans*30		!to hold text & result of QDIALOG
	logical help,nodata
	integer jmiss(njset),juse(njset)
	integer iuse(njset),icuse(100)
c for popmenu
	character*60 strings(30)
      character*55 helps(5)	!for popmenu
	character*20 title
      character*1 charout 		!for popmenu
	common/popvalm/nxlo,nylo,nyhi,ictm,ibkm,icfm,icupm	!values for popmenu calls
c
	common/hlp/help				!for QDIALOG
	common/dp/discprt
	COMMON/BLOCK1/constr,nset,nfit,nsfit,
     & Xv,alfa,kmax,ncomp,nmod,fline,nomit,jomit(20),jset,ifitmode
	COMMON/BLOCK2/ castar,X1,X2,iequiv,ip1,ip2
	common/potrat/jfirst,iset,kmax1			!for pot ratios
c (kmax=actual no of param to be estimated; kmax1=number before extra param
c  added for ifitmode=4,5)
	logical logyfit,norm
	COMMON/BLOCK3/logyfit,norm,xnorm		!CVIN2,YCALCV,CVDISP,GETEQN
	common/binomial/binom(0:10)
c
c Oct 86- two new variables added to common/block2/
c IEQUIV=0 for equiv subunits; =1 for non-equiv (prev signalled via castar)
c IP1,IP2:
c  (1) mods 10,11 IP1=1 if param are K1,K2; =2 if param are K1, K2/K1.
c  (2) mods 9,11  IP2=1 if param is K2 or (K2/K1 if IP1=2)
c		  IP2=2 if param is k(-2) (must have ip2=1 if ip1=2)
c  (3) mods 9,11  IP2 set negative if k(+2) is a param  to be estimated
c  (4) nmod=99 =polynomial. Degree of poly=ip1, and if ip2>0 then degree
c	is ip2 for x>x1.
c	If ip2=-1 then Schild eqn is fitted with data Y=dose ratio, X=antag conc
c	If ip2=-2 then Schild eqn is fitted with data Y=dose ratio-1, X=antag conc
c  (5) nmod=30  ip1=number of agonist binding sites
c
c	DATA tt1/' b1 ',' b2 ',' b3 ',' b4 ',' b5 ',' b6 ',' b7 ',
c     & ' b8 ',' b9 ','b10 '/
c
	pon()=slock()
c
c101	format(a1)
c
	alpha=VIDEOTYP().eq.3		!alpha mode
c
	pprt=pon()		!to reduce number of pon() calls (interrupt 11 problems)
c	iequiv=0		!initialise in case -no -default value is input now
	titprint=.false.
	T1(1)='       K1 '
	t1(2)='    Ymax1 '
	t1(3)='       K2 '
	t1(4)='    Ymax2 '
	t1(5)='    Slope '
	t1(6)='    Y1(0) '
	t1(7)='   Y(inf) '
	t1(8)='    k(-1) '
	t1(9)='     k(1) '
	t1(10)='     Bmax '
	t1(11)='      B/A '
	t1(12)='     Beta '
	t1(13)='       KB '
	t1(14)='     Y(0) '
	t1(15)='  N(Hill) '
	t1(16)='    K2/K1 '
	t1(17)=' 	k(-2) '
	t1(18)='     k(2) '
	t1(19)='   B2 /A2 '
	t1(20)='   Non-eq '
	t1(21)='     Coop '
	t1(22)=' Erev(mV) '
	t1(23)='   g(E=0) '
	t1(24)='    H (mV)'
	t1(25)=' B/A (E=0)'
	t1(26)='   g(max) '
	t1(27)='   conc/K '
	t1(28)='   tau(0) '
	t1(29)='g(Na) (pS)'
	t1(30)='g(Cs) (pS)'
	t1(31)='E(off)(mV)'
	t1(32)='H(K1) (mV)'
	t1(33)='H(K2) (mV)'
	t1(34)='     nH1  '
	t1(35)='    Y2(0) '
	t1(36)='     nH2  '
	t1(37)='  IC50(1) '
	t1(38)='  IC50(2) '
c
	titmod(1)='(1) Polynomial (inc. straight line)'			!nmod=99
	titmod(2)='(2) Langmuir hyperbola(s) (inc. or dec.) '		!nmod=1/26
	titmod(3)='(3) Hill equation(s) (inc. or dec./common K or max)'!nmod=-1/27
c	titmod(4)='(4) Hyperbola plus straight line'			!nmod=2
c	titmod(5)='(5) (Two hyperbolas: use #2 now)'			!nmod=3
c	titmod(6)='(6) Two hyperbolas plus straight line'		!nmod=4
	titmod(4)='(4) Langmuir hyperbola(s) plus straight line'	!nmod=26			!nmod=2
	titmod(5)='(5) Hill equation(s) plus straight line'  		!nmod=27			!nmod=2
c	titmod(6)='(6)    (not used)                       '
	titmod(6)=
     & '(6) Power function y=ybar*(x/x0)^n (linear on log-log plot)'
	titmod(7)='(7) Binding inhibition curve (parameter=KB)'		!nmod=5
	titmod(8)='(8) Exponential binding onset curves'
c	titmod(8)='(8) Exponential binding onset (common overall max)'
c	titmod(9)='(9) Exponential binding onset (separate max)'
	titmod(9)='(9) General fit of sums of exponentials/geometrics'
	titmod(10)='(10) Popen FOR KM2 MODEL+ BLOCK. K1=K2'
	titmod(11)='(11) Popen FOR KM2 + BLOCK (BETA as param). K1=K2'
	titmod(12)='(12) Popen for KM2 at equilibrium + block: SEP K1,K2'
	titmod(13)='(13) Popen for KM2 + block (BETA as param)'
	titmod(14)=
     & '(14) Total current (I/V) for KM2 at equil + block (K1=K2)'
	titmod(15)='(15) Popen for 2 independent KM1 subunits'
	titmod(16)='(16) Popen for KM2 (Coop and nonequiv) + block'
	titmod(17)=
     & '(17) Popen for KM2 (Coop and nonequiv)+block (beta as param)'
	titmod(18)=
     & '(18) I/V with exponential voltage-dependence of current'
	titmod(19)=
     & '(19) Ditto with saturation (beta/alpha V-dep)'
	titmod(20)=
     & '(20) Ditto with saturation (affinity V-dep)'
	titmod(21)=
     & '(21) I/V for GHK with Mg block'
	titmod(22)=
     & '(22) Exponential voltage dependence of rate constant'
c	titmod(23)=			!not used (now 6)
c     & '(23) Power function y=ybar*(x/x0)^n (linear on log-log plot)'
	titmod(23)='   (not used)                             '
        titmod(24)=
     & '(24) Bind n mols + opening: # of channels same for all sets'
        titmod(25)=
     & '(25) MWC mechanism (any n): # of channels same for all sets'
        ntmod=25                !total number of models defined
c=======temp stuff for Woda problem
        titmod(26)=
     & '(26) Woda equation                                         '
        ntmod=26                !total number of models defined
c=======end of temp
c	titmod(0)=
c     & ''
c
c=	if(sameq) goto 317		!already decided to fit same eq
	if(sameq) goto 316		!already decided to fit same eq
      if(iflag.eq.1) then             !2nd, 3rd,... time
	   call INTCONV(imod,cnum)
	   if(alpha) then
	      call DCASK(
     &	'Fit same equation (#'//charnb(cnum)//')','Y',ans)
	   else
		ans='Y'
		call DEFOLTa(ans,defolt)
		call QDIALOG(1,'Fit same equation (#'//charnb(cnum)//')',
     &	   defolt,11,cans)
		call GETINPa(cans,ans)
	   endif
         if(UC(ans).eq.'Y') goto 316
c=         if(UC(ans).eq.'Y') goto 317
	endif
64    continue
c
	logyfit=.false.
	fline=.false.
c
	if(alpha) then
         print 1151,imod
1151	   FORMAT(
     &' Number of equation to be fitted (-1=abort; 0=HELP) [',i3,'] = ')
	   call INPUTi(imod)
	   if(imod.eq.0) then
		do i=1,ntmod
		   print 66,titmod(i)
66		   format(1x,a60)
		enddo
	      goto 64	!ask for mod # again
	   endif
	else
	   if(imod.lt.1) imod=1
	   call DEFOLTi(imod,defolt)
25	   call QDIALOG(1,
     &     'No. of equation to be fitted (-1=abort) [F1=HELP]',
     &	defolt,11,cans)
	   if(help) then				!F1 hit in qdialog
		call CVHELP(1,titmod,imod)	!outputs imod
		call intconv(imod,charstr)
		call WDIALOG(1,'Fit with equation no. '//charstr,11)
	   else
		call GETINPi(cans,imod)
	   endif
	endif
	if(imod.lt.0) then
	   iesc=99
	   RETURN
	endif
	if(imod.gt.ntmod) goto 64
c
c Define original model number (nmod) from imod to save extensive changes
c later
317	continue			!jump here if 'same equation'
c	if(nmod.eq.0) nmod=99
	if(imod.eq.1) then
	   nmod=99				!polynomial
	else if(imod.eq.2) then
	   nmod=1				!hyperbola
	   fline=.false.
	else if(imod.eq.3) then
	   nmod=-1				!Hill eqn
	   fline=.false.
	else if(imod.eq.4) then
	   nmod=1	!later 26
	   fline=.true.
	else if(imod.eq.5) then
	   nmod=-1	!later 27
	   fline=.true.
	else if(imod.eq.6) then
	   nmod=28	!power function
	else if(imod.eq.8) then       !old imod=8,9; nmod=6,7; exponential onset
24	   continue
	   if(alpha) then
		print 23
23		format(
     &     ' (1) Fit common overall maximum',/,
     &     ' (2) Fit separate maximum for each curve',/,
     &     ' Option number = ')
		call INPUTi(iopt)
	   else
		call WDIALOG(1,'(1) Fit common overall maximum',11)
		call WDIALOG(1,'(2) Fit separate maximum for each curve',11)
		iopt=1
		call DEFOLTi(iopt,defolt)
		call QDIALOG(1,'Option number',defolt,11,cans)
		call GETINPi(cans,iopt)
	   endif
	   if(iopt.lt.1.or.iopt.gt.2) goto 24
	   if(iopt.eq.1) nmod=6
	   if(iopt.eq.2) nmod=7
	else if(imod.eq.9) then		!new general exponentials
	   nmod=19
	   if(alpha) then
	      print 38
38		format(' Fit geometrically decaying curve [N] ')
		ans='N'
		call INPUTa(ans)
	   else
		ans='N'
		call DEFOLTa(ans,defolt)
		call QDIALOG(1,'Fit geometrically decaying curve',
     &		defolt,11,cans)
		call GETINPa(cans,ans)
	   endif
	   if(UC(ans).ne.'Y') then
		nmod=19		!exponentials
	   else
		nmod=20		!geometrics
	   endif
	else if(imod.ge.7.and.imod.le.20.and.imod.ne.8.and.imod.ne.9) then
	   nmod=imod-2		!nmod=5,..,18 for imod=7,...,20 (exc imod=8,9)
c	else if(imod.eq.23) then	!not used now
c	   nmod=28	!power function
	else if(imod.eq.24) then
         nmod=30      !n bindings +opening (nmod=30)
	else if(imod.eq.25) then
         nmod=31      !MWC(n) (nmod=31)
	else
	   nmod=imod			!new ones nmod=21,22 (nmod=20 for geometrics now)
	endif
c=======temp stuff for Woda problem
	if(imod.eq.26) then
	   nmod=100
	   kmax=3
	   titlep(1)='       a    '
	   titlep(2)='       b    '
	   titlep(3)='       c    '
	   goto 316
	endif
c=======end of temp
	if(nodata) then
	   if(pprt) write(7,131) titmod(imod)
	   if(discprt) write(8,131) titmod(imod)
131	   format(/,' Calculated curves:',/,' Equation = ',a60)
	endif
c Common maximum for Hill/Langmuir when more than one set fitted?
c (omit this question if ifitmode=4,5 because these already specified
c to do fit with all param except pot ratio/dose ratio in common)
c For ifitmode=3, choose which param are in common
c For ifitmode=2, no params in common -sep fit
c Use ip1,ip2 to signal options concerning common parameters when ifitmode=3
c ip1=-1 -separate max, separate K
c ip1=0 -common max, separate K
c ip1=1 -common K, separate max
c ip1=2 -separate K, Ymax but common Ymax1/Ymax (ncomp=2 only)
c ip2=0 -separate nH
c ip2=1 -common nH
c If fline set true, then add linear component -this can be done
c for any of the options listed here, except when ifitmode=4 or 5.
c No option yet to have a common slope -if slope fitted there is
c a separate value for each set at present
	if((nmod.eq.1.or.nmod.eq.-1)) then
	   if(ncomp.lt.1) ncomp=1
	   if(alpha) then
		print 541,ncomp
541		format(' Number of components for each set [',i2,'] = ')
		call INPUTi(ncomp)
	   else
		call DEFOLTi(ncomp,defolt)
		call QDIALOG(1,
     &	'Number of components for each set',
     &	defolt,11,cans)
c	      if(help) then				!F1 hit in qdialog
c		   call CVHELP(3,titmod,imod)
c		   goto 28
c	      else
		   call GETINPi(cans,ncomp)
c		endif
	   endif
241	   continue
	   if(ifitmode.eq.4.or.ifitmode.eq.5) then
		if(nmod.eq.1) then
		   nmod=26	!Langmuir
		else
		   nmod=27	!Hill
		endif
		titlep(1)='       Y(0) '  !same for all
		call INTCONV(jset,cnum)
		n=nblank(cnum)
		k=1
c for ifitmode=4,5 Ymax and nH same for all sets
		do i=1,ncomp
		   call INTCONV(i,cnum1)
		   n1=nblank(cnum1)
		   if(.not.norm) then
		      k=k+1
		      titlep(k)='       Ymax'//cnum1(1:n1)
		      icuse(k)=iuse(1)
		   endif
		   k=k+1
		   titlep(k)='   K'//cnum1(1:n1)//' (set 1)'
		   icuse(k)=iuse(1)
		   if(nmod.eq.27) then
			k=k+1
			titlep(k)='       nH'//cnum1(1:n1)
		   endif
		enddo
		kmax=k
		goto 316	!no option to fit line
	   endif
c
	   if(ifitmode.eq.3) then	!one eqn, several sets
	    if(alpha) then
		if(ncomp.eq.1) then
		   print 31
31		   format(
     &     ' (1) Fit separate maximum, separate K values',/,
     &     ' (2) Fit common maximum, separate K values',/,
     &     ' (3) Fit common K value, separate maxima',/,
     &     ' Option number = ')
		else if(ncomp.eq.2) then
		   print 311
311		   format(
     &     ' (1) Fit separate maximum, separate K values',/,
     &     ' (2) Fit common maximum, separate K values',/,
     &     ' (3) Fit common K value, separate maxima',/,
     &     ' (4) Separate K and Ymax, but common percent Ymax1',/,
     &     ' Option number = ')
		endif
		call INPUTi(iopt)
		if(iopt.lt.1.or.iopt.gt.4) goto 241
		ans='Y'	!default
		if(iequiv.eq.1) ans='N'
		print 33,ans
33		format(' Y increasing with X [',a1,'] ? ')
		call INPUTa(ans)
		iequiv=0		     		!increasing
		if(UC(ans).eq.'N') iequiv=1	!decreasing
		if(iequiv.eq.1) titlep(1)='   Y(inf) '  !same for all
		if(nmod.eq.-1) then
		   print 32
32		   format(' Hill slope same for all sets [Y] ? ')
		   ans='Y'
		   call INPUTa(ans)
		   ip2=0		     		!separate nH
		   if(UC(ans).ne.'N') ip2=1	!common nH
		endif
	    else
c		call WDIALOG(1,
c     &	  '(1) Fit separate maximum, separate K values',11)
c		call WDIALOG(1,
c     &	  '(2) Fit common maximum, separate K values',11)
c		call WDIALOG(1,
c     &	  '(3) Fit common K value, separate maxima',11)
c		if(ncomp.eq.1.and.(.not.fline)) then
c		  call WDIALOG(1,
c     &	    '(4) ditto using fast fit to log(y)',11)
c		else if(ncomp.eq.2) then
c		  call WDIALOG(1,
c     &	  '(4) Separate K and Ymax, but common percent Ymax1',11)
c		endif
c		iopt=1
c		call DEFOLTi(iopt,defolt)
c		call QDIALOG(1,'Option number',defolt,11,cans)
c		call GETINPi(cans,iopt)
c		if(iopt.lt.1.or.iopt.gt.4) goto 241
		title='Fit options'
		strings(1)='Fit separate maximum, separate K values'
		strings(2)='Fit common maximum, separate K values'
		strings(3)='Fit common K value, separate maxima'
		if(ncomp.eq.1.and.(.not.fline)) then
		   strings(4)='ditto using fast fit to log(y)'
		   nval=4
		else if(ncomp.eq.2) then
		   strings(4)='Separate K and Ymax, but common % Ymax1'
		   nval=4
		else
		   nval=3
		endif
		helps(1)='  '
		nhelp=1
		iline=1	!default
	      call POPMENU(nxlo,nylo,nyhi,strings,nval,ictm,icfm,icupm,
     &	 ibkm,title,helps,nhelp,iline,charout,ival)
		if(iline.eq.0) then	!iline=0 for ESC=cancel
		   goto 241
		else
		   iopt=iline
		endif
c
		ans='Y'	!default
		if(iequiv.eq.1) ans='N'
		call DEFOLTa(ans,defolt)
		call QDIALOG(1,'Y increasing with X',
     &		defolt,11,cans)
		call GETINPa(cans,ans)
		iequiv=0		     		!increasing
		if(UC(ans).eq.'N') iequiv=1	!decreasing
		if(nmod.eq.-1) then
		   ans='Y'
		   call DEFOLTa(ans,defolt)
		   call QDIALOG(1,'Hill slope same for all sets',
     &		defolt,11,cans)
		   call GETINPa(cans,ans)
		   ip2=0			!separate nH
		   if(UC(ans).eq.'Y') ip2=1
		endif
	    endif
	   endif	!end of ifitmode=3 only bit
c Following bit applies for normalised responses too
c IFITMODE=2 -FITTING ONLY ONE SET AT A TIME so no question of common
c anything -just fit ncomp components separately to each set
	   if(ifitmode.eq.1.or.ifitmode.eq.2) then
		iopt=4	!no common parameters!
		if(alpha) then
		   ans='Y'	!default
		   if(iequiv.eq.1) ans='N'
		   print 33,ans
c33		   format(' Y increasing with X [',a1,'] ? ')
		   call INPUTa(ans)
		   iequiv=0		     		!increasing
		   if(UC(ans).eq.'N') iequiv=1	!decreasing
		else
		   ans='Y'	!default
		   if(iequiv.eq.1) ans='N'
		   call DEFOLTa(ans,defolt)
		   call QDIALOG(1,'Y increasing with X',
     &		defolt,11,cans)
		   call GETINPa(cans,ans)
		   iequiv=0		     		!increasing
		   if(UC(ans).eq.'N') iequiv=1	!decreasing
		endif
	   endif		!end of ifitmode=1,2
c
	   if(iw.ge.1.and.iw.le.5) then
		if(pprt) write(7,1301) titmod(imod),titw(iw)
	      if(discprt) write(8,1301) titmod(imod),titw(iw)
c1301	  format(/,' FITTING:',/,' Equation = ',a60,/,' Weighting = ',a60)
		titprint=.true.
	   endif
	   if(iopt.eq.1) then		!sep max, sep K
		if(pprt) write(7,735)
	      if(discprt) write(8,735)
735	     format(' Separate maximum, separate K values for each set',/)
		ip1=-1
		titlep(1)='     Y(0) ' 			    !same for all
		if(iequiv.eq.1) titlep(1)='   Y(inf) '  !same for all
		k=1
		do j=1,nsfit
		   j1=juse(j)		!set number used
		   call INTCONV(j1,cnum)
		   n=nblank(cnum)
		   do i=1,ncomp
		      call INTCONV(i,cnum1)
		      n1=nblank(cnum1)
			k=k+1
			titlep(k)='K'//cnum1(1:n1)//' (set '//cnum(1:n)//')'
		      icuse(k)=iuse(j1)		!curve # for set j1
			if(.not.norm) then
			  k=k+1
			  titlep(k)=
     &			'Ymax'//charnb(cnum1)//' (set '//cnum(1:n)//')'
		        icuse(k)=iuse(j1)		!curve # for set j1
			endif
		   enddo
		enddo
	   else if(iopt.eq.2) then		!common max, sep K
		if(pprt) write(7,35)
	      if(discprt) write(8,35)
35		format(' Common maximum, separate K values for each set',/)
		ip1=0
		titlep(1)='     Y(0) ' 			    !same for all
		if(iequiv.eq.1) titlep(1)='   Y(inf) '  !same for all
		k=1
		do i=1,ncomp	!ymax(i) for component #i (same for all sets)
		   call INTCONV(i,cnum)
		   n=nblank(cnum)
		   if(.not.norm) then
		      k=k+1
		      titlep(k)='    Ymax'//charnb(cnum)//' '	!='Ymax1','Ymax2',...
		   endif
		enddo
		do j=1,nsfit
		   j1=juse(j)		!set number used
		   call INTCONV(j1,cnum)
		   n=nblank(cnum)
c		   chj=char(j1+48)
		   do i=1,ncomp
		      call INTCONV(i,cnum1)
		      n1=nblank(cnum1)
c			chi=char(i+48)
			k=k+1
			titlep(k)='K'//cnum1(1:n1)//' (set '//cnum(1:n)//')'
		      icuse(k)=iuse(j1)		!curve # for set j1
		   enddo
		enddo
	   else if(iopt.eq.3) then	!common K, sep max
		if(pprt) write(7,36)
	      if(discprt) write(8,36)
36		format(' Common K value, separate maxima for each set',/)
		ip1=1
		titlep(1)='     Y(0) '  !same for all
		if(iequiv.eq.1) titlep(1)='   Y(inf) '  !same for all
		k=1
		do j=1,nsfit
		   j1=juse(j)		!set number used
		   call INTCONV(j1,cnum)
		   n=nblank(cnum)
		   do i=1,ncomp
		      call INTCONV(i,cnum1)
		      n1=nblank(cnum1)
		      if(.not.norm) then
			   k=k+1
		         titlep(k)=
     &			'Ymax'//cnum1(1:n1)//' (set '//cnum(1:n)//')'
		         icuse(k)=iuse(j1)		!curve # for set j1
			endif
		   enddo
		enddo
		do i=1,ncomp	!K(i) for component #i (same for all sets)
		   call INTCONV(i,cnum1)
		   n1=nblank(cnum1)
		   k=k+1
		   titlep(k)='       K'//cnum1(1:n1)//' '	!='K1','K2',...
		enddo
	   else if(iopt.eq.4.and.ifitmode.eq.3) then
		if(ncomp.eq.1) then	!common K, sep max by fitting log(y)
		   if(pprt) write(7,37)
	         if(discprt) write(8,37)
37		   format(' Common K value, separate maxima for each set',/,
     &	' Fit to log(y) values with equal weight for log values',/)
		   ip1=1
		   titlep(1)='     Y(0) '  !same for all
		   if(iequiv.eq.1) titlep(1)='   Y(inf) '  !same for all
		   if(alpha) then
			print 34,titlep(1)(4:9)
34			format(' This method assumes ',a6,
     &		' =0: [Y] ')
			ans='Y'
			call INPUTa(ans)
		   else
			ans='Y'
			call DEFOLTa(ans,defolt)
			call QDIALOG(1,'This method assumes '//
     &		  titlep(1)(4:9)//' = 0: O.K.',
     &		  defolt,11,cans)
		      call GETINPa(cans,ans)
		   endif
		   if(UC(ans).eq.'N') goto 241
		   logyfit=.true.
		   titlep(1)='     K1   '
		   k=1
		else if(ncomp.eq.2) then
c    New option for ncomp=2 only -common Ymax1/ymax2 -in this case fit K,
c    Ymax (and nH if req), one value for each set (not separatt Ymax1, Ymax2)
c    plus a single value of Ymax1/Ymax for all sets
c Fit percent comp 1, say P1 = Ymax1/(Ymax1+Ymax2) = Ymax1/Ymax.  Then to
c get Ymax1 and Ymax2 (eg in ycalcv) for each set, find
c Ymax1 = P1*Ymax, then Ymax2=Ymax-Ymax1
		   if(alpha) then
			print 70
70			format(
     &	     ' (1) Estimate common fraction of component 1',/,
     &	     ' (2) Estimate common fraction of component 2',/,
     &	     ' Option number [1] = ')
			i=1
			call INPUTi(i)
		   else
			call WDIALOG(1,
     &		'(1) Estimate common fraction of component 1',11)
			call WDIALOG(1,
     &		'(2) Estimate common fraction of component 2',11)
			i=1
			call DEFOLTi(i,defolt)
			call QDIALOG(1,'Option number',defolt,11,cans)
			call GETINPi(cans,i)
		   endif
		   if(i.eq.1) then
			ip1=2
		   else
			ip1=3
		   endif
		   if(pprt) write(7,361) i
	         if(discprt) write(8,361) i
361		   format(
     &  ' Separate K and Ymax, but common fraction of component ',i2,/)
		   titlep(1)='     Y(0) '  !same for all
		   if(iequiv.eq.1) titlep(1)='   Y(inf) '  !same for all
		   k=1
		   do j=1,nsfit
			j1=juse(j)		!set number used
			call INTCONV(j1,cnum)
			n=nblank(cnum)
			k=k+1
		      titlep(k)='Ymax (tot) '//'(set '//cnum(1:n)//')' !1 Ymax for each set
		   	icuse(k)=iuse(j1)		!curve # for set j1
		   enddo
		   do j=1,nsfit
			j1=juse(j)		!set number used
			call INTCONV(j1,cnum)
			n=nblank(cnum)
			do i=1,ncomp	!K(i) for component #i (i=1,2 at present)
			   call INTCONV(i,cnum1)
			   n1=nblank(cnum1)
			   k=k+1
			   titlep(k)='   K'//cnum1(1:n1)//' (set '//
     &			cnum(1:n)//')'
		   	   icuse(k)=iuse(j1)		!curve # for set j1
		      enddo
		   enddo
		endif
	   else if(iopt.eq.4) then	!for ifitmode=1, 2: no common params
		if(pprt) write(7,39)
	      if(discprt) write(8,39)
39		format(' Separate fit to each set: no common parameters',/)
		if(nmod.eq.1) then
		   nmod=26	!Langmuir
		else
		   nmod=27	!Hill
		endif
		titlep(1)='     Y(0) '  !same for all
		if(iequiv.eq.1) titlep(1)='   Y(inf) '  !same for all
c	fitting only set #jset at each pass for ifitmode=2
		call INTCONV(jset,cnum)
		n=nblank(cnum)
		k=1
		do i=1,ncomp
		   call INTCONV(i,cnum1)
		   n1=nblank(cnum1)
		   if(.not.norm) then
		      k=k+1
		      titlep(k)=
     &		 'Ymax'//cnum1(1:n1)//' (set '//cnum(1:n)//')'
		      icuse(k)=iuse(jset)		!curve # for set jset
		   endif
		   k=k+1
		   titlep(k)='   K'//cnum1(1:n1)//' (set '//cnum(1:n)//')'
		   icuse(k)=iuse(jset)		!curve # for set jset
		   if(nmod.eq.27) then
			k=k+1
			titlep(k)='nH'//cnum1(1:n1)//' (set '//cnum(1:n)//')'
		      icuse(k)=iuse(jset)		!curve # for set jset
		   endif
		enddo
		if(fline) then
		   k=k+1
		   titlep(k)='slope (set '//cnum(1:n)//')'
		   icuse(k)=iuse(jset)		!curve # for set jset
		endif
		kmax=k
		goto 316	!all set in case of iopt=4 only
	   endif
c Add nH parameters if Hill (not added here if ifitmode=1,2 -they
c were set above, in a preferable order, under iopt=4)
	   if(nmod.eq.1) then
		nmod=26	!Langmuir
		kmax=k
	   else if(nmod.eq.-1) then
		nmod=27	!Hill
		if(ip2.eq.1) then		!common nH (1 for each component)     		!separate nH
		   do i=1,ncomp	!ymax(i) for component #i (same for all sets)
		      call INTCONV(i,cnum)
		      n=nblank(cnum)
			k=k+1
			titlep(k)='     nH'//cnum(1:n)//' '	!='nH1','nH2',...
		   enddo
		else				!sep nH for each set
		   do j=1,nsfit
			j1=juse(j)		!set number used
		      call INTCONV(j1,cnum)
		      n=nblank(cnum)
			do i=1,ncomp
		         call INTCONV(i,cnum1)
		         n1=nblank(cnum1)
			   k=k+1
			 titlep(k)='nH'//cnum1(1:n1)//' (set '//cnum(1:n)//')'
			   icuse(k)=iuse(j1)		!curve # for set j1
			enddo
		   enddo
		endif
		kmax=k
	   endif
c add the percent Ymax1 parameter last if ncomp=2 and ip1=2
	   if(ifitmode.eq.3.and.ncomp.eq.2.and.ip1.eq.2) then
		kmax=kmax+1
		titlep(kmax)='fraction of comp 1 '
	   else if(ifitmode.eq.3.and.ncomp.eq.2.and.ip1.eq.3) then
		kmax=kmax+1
		titlep(kmax)='fraction of comp 2 '
	   endif
c lastly add a straight line for each set
	   if(fline) then
		k=kmax	!add more to kmax!
		do j=1,nsfit
		   j1=juse(j)		!set number used
		   call INTCONV(j1,cnum)
		   n=nblank(cnum)
		   do i=1,ncomp
		      call INTCONV(i,cnum1)
		      n1=nblank(cnum1)
			k=k+1
			titlep(k)='slope'//cnum1(1:n1)//' (set '//
     &		   cnum(1:n)//')'
			icuse(k)=iuse(j1)		!curve # for set j1
		   enddo
		enddo
		kmax=k	!new final value
	   endif
c
	   goto 316
	endif		!end of Langmuir/Hill section
c
	if(nmod.eq.30.or.nmod.eq.31) then	!n bindings + opening/MWCn
	   nbind=3
10	   call DEFOLTi(nbind,defolt)
	   call QDIALOG(1,
     &   'number of agonist binding sites',
     &   defolt,11,cans)
	   call GETINPi(cans,nbind)
	   if(nbind.gt.9) goto 10
	   call INTCONV(nbind,cnum)
	   n=nblank(cnum)
	   ip1=nbind	!keep value for common
	   if(nmod.eq.30) then
	      titmod(24)=
     & 	'(24) Bind '//cnum(1:n)//
     &	' mols + opening: # of channels same for all sets'
		ans='Y'
		call DEFOLTa(ans,defolt)
41		call QDIALOG(1,
     &	   'All binding constants same (R=1: F1 for help)',
     &	   defolt,11,cans)
		if(help) then				!F1 hit in qdialog
		   call CVHELP(4,titmod,imod)		!outputs imod
		   goto 41
		else
		   call GETINPa(cans,ans)
		endif
		iequiv=0
		if(ans.eq.'N') iequiv=1       !keep value in common
	   else if(nmod.eq.31) then
	      titmod(25)=
     & 	'(25) MWC mechanism with '//cnum(1:n)//
     &	' binding sites: # of channels same for all sets'
		iequiv=0
		if(nsfit.gt.1) then
		   title='MWC options'
		   strings(1)='All parameters fitted (use En)'
		   strings(2)='All parameters fitted (use E0)'
		   strings(3)='common E0 for all sets'
		   strings(4)='common R=K/K* for all sets'
		   strings(5)='all common except E0'
		   strings(6)='all common except R=K/K*'
		   nval=6
		   helps(1)=' Full MWC may use either E0 or En  '
		   helps(2)=' as the free parameter (options 1,2)'
		   nhelp=2
		   iline=1	!default
	         call POPMENU(nxlo,nylo,nyhi,strings,nval,ictm,icfm,icupm,
     &	 ibkm,title,helps,nhelp,iline,charout,ival)
		   if(iline.eq.0) then	!iline=0 for ESC=cancel
			goto 64
		   else
			if(iline.eq.1) then
			   iequiv=0       !=0,-1,1,2; keep value in common
			else if(iline.eq.2) then
			   iequiv=-1
			else if(iline.ge.3) then
			   iequiv=iline-2       !1,2,3,4
			endif
		   endif
		endif
	   endif
c      Define binomial coeff here so no need to recalc each time in ycalcv
c      Keep in common/binomial/
	   an1=FACT(nbind)
	   do i=0,nbind
		binom(i)=an1/(FACT(i)*FACT(nbind-i))
	   enddo
c
	   iopt=1	! simple case where all K the same -do others later?
	   if(nmod.eq.30.and.iopt.eq.1) then
		if(ifitmode.le.2) then
		   k=0
		   do j=1,nsfit
			j1=juse(j)		!set number used
			call INTCONV(j1,cnum)
			n=nblank(cnum)
			k=k+1
			titlep(k)='K (set '//cnum(1:n)//')'
			icuse(k)=iuse(j1)		!curve # for set j1
			k=k+1
			titlep(k)='E (set '//cnum(1:n)//')'
			icuse(k)=iuse(j1)		!curve # for set j1
			if(iequiv.eq.1) then
			   k=k+1
			   titlep(k)='r (set '//cnum(1:n)//')'
			   icuse(k)=iuse(j1)		!curve # for set j1
			endif
			k=k+1
			titlep(k)='Ymax (set '//cnum(1:n)//')'
			icuse(k)=iuse(j1)		!curve # for set j1
		   enddo
		else if(ifitmode.eq.3) then
		   k=0
		   do j=1,nsfit
			j1=juse(j)		!set number used
			call INTCONV(j1,cnum)
			n=nblank(cnum)
			k=k+1
			titlep(k)='K (set '//cnum(1:n)//')'
			icuse(k)=iuse(j1)		!curve # for set j1
			k=k+1
			titlep(k)='E (set '//cnum(1:n)//')'
			icuse(k)=iuse(j1)		!curve # for set j1
			if(iequiv.eq.1) then
			   k=k+1
			   titlep(k)='r (set '//cnum(1:n)//')'
			   icuse(k)=iuse(j1)		!curve # for set j1
			endif
		   enddo
		   k=k+1
		   titlep(k)='Ymax [at P(open)=1]'
		endif
		kmax=k
		goto 316
	   else if(nmod.eq.31) then	!MWC
c		if(ifitmode.le.2) then
c		   k=0
c		   do j=1,nsfit
c			j1=juse(j)		!set number used
c			call INTCONV(j1,cnum)
c			n=nblank(cnum)
c			k=k+1
c			titlep(k)='K (shut) (set '//cnum(1:n)//')'
c		      icuse(k)=iuse(j1)		!curve # for set j1
c			k=k+1
c			titlep(k)='R (K/K*) (set '//cnum(1:n)//')'
c		      icuse(k)=iuse(j1)		!curve # for set j1
c			k=k+1
c			titlep(k)='En=E0*R^n (set '//cnum(1:n)//')'
c		      icuse(k)=iuse(j1)		!curve # for set j1
c			k=k+1
c			titlep(k)='Ymax (set '//cnum(1:n)//')'
c		      icuse(k)=iuse(j1)		!curve # for set j1
c		   enddo
c		else if(ifitmode.eq.3) then
		if(ifitmode.le.3) then
		   k=0
		   if(iequiv.le.0) then		!fit all
			do j=1,nsfit
			   j1=juse(j)		!set number used
			   call INTCONV(j1,cnum)
			   n=nblank(cnum)
			   k=k+1
			   titlep(k)='K (shut) (set '//cnum(1:n)//')'
		         icuse(k)=iuse(j1)		!curve # for set j1
			   k=k+1
			   titlep(k)='R (K/K*) (set '//cnum(1:n)//')'
		         icuse(k)=iuse(j1)		!curve # for set j1
			   k=k+1
			   if(iequiv.eq.-1) then
				titlep(k)='E0=En/R^n (set '//cnum(1:n)//')'
		      	icuse(k)=iuse(j1)		!curve # for set j1
			   else if(iequiv.eq.0) then
				titlep(k)='En=E0*R^n (set '//cnum(1:n)//')'
		      	icuse(k)=iuse(j1)		!curve # for set j1
			   endif
			enddo
		   else if(iequiv.eq.1) then		!common E0
			do j=1,nsfit
			   j1=juse(j)		!set number used
			   call INTCONV(j1,cnum)
			   n=nblank(cnum)
			   k=k+1
			   titlep(k)='K (shut) (set '//cnum(1:n)//')'
		         icuse(k)=iuse(j1)		!curve # for set j1
			   k=k+1
			   titlep(k)='R (K/K*) (set '//cnum(1:n)//')'
		         icuse(k)=iuse(j1)		!curve # for set j1
			enddo
			k=k+1
			titlep(k)='E0 (common value)'
		   else if(iequiv.eq.2) then		!common R
			do j=1,nsfit
			   j1=juse(j)		!set number used
			   call INTCONV(j1,cnum)
			   n=nblank(cnum)
			   k=k+1
			   titlep(k)='K (shut) (set '//cnum(1:n)//')'
		         icuse(k)=iuse(j1)		!curve # for set j1
			   k=k+1
			   titlep(k)='En=E0*R^n (set '//cnum(1:n)//')'
		         icuse(k)=iuse(j1)		!curve # for set j1
			enddo
			k=k+1
			titlep(k)='R (K/K*) (common value)'
		   else if(iequiv.eq.3) then		!all common except E0
			k=k+1
			titlep(k)='K (shut) (common value)'
			k=k+1
			titlep(k)='R (K/K*) (common value)'
			do j=1,nsfit
			   j1=juse(j)		!set number used
			   call INTCONV(j1,cnum)
			   n=nblank(cnum)
			   k=k+1
			   titlep(k)='E0=En/R^n (set '//cnum(1:n)//')'
		         icuse(k)=iuse(j1)		!curve # for set j1
			enddo
		   else if(iequiv.eq.4) then		!all common except R=K/K*
			k=k+1
			titlep(k)='K (shut) (common value)'
			do j=1,nsfit
			   j1=juse(j)		!set number used
			   call INTCONV(j1,cnum)
			   n=nblank(cnum)
			   k=k+1
			   titlep(k)='R (K/K*) (set '//cnum(1:n)//')'
		         icuse(k)=iuse(j1)		!curve # for set j1
			enddo
			k=k+1
			titlep(k)='E0=En/R^n (common value)'
		   endif
		   k=k+1
		   titlep(k)='Ymax [at P(open)=1]'
	  	endif
		kmax=k
		goto 316
	   endif
	endif		!end of nmod=30/31
c General exponentials
	if(nmod.eq.19) then
	   xgt0=.true.
	   if(alpha) then
		print 531
531		format('&Allow fitted curves to go below x=0 [N] ? ')
		ans='N'
		call INPUTa(ans)
		if(UC(ans).eq.'Y') xgt0=.false.
		print 53
53		format(
     &	  '&Number of exponential components to be fitted [1] = ')
		call INPUTi(ncomp)
	   else
		ans='N'
		call DEFOLTa(ans,defolt)
		call QDIALOG(1,'Allow fitted curves to go below x=0',
     &		defolt,11,cans)
		call GETINPa(cans,ans)
		if(ans.eq.'Y') xgt0=.false.
		call DEFOLTi(ncomp,defolt)
		call QDIALOG(1,
     &	'Number of exponential components to be fitted',
     &	defolt,11,cans)
		call GETINPi(cans,ncomp)
	   endif
	   if(ncomp.eq.0) ncomp=1
	   kmax=2*ncomp+1		!no of parameters
	   do 5 i=1,ncomp
		ix=2*i-1		!=1,3,5,7..
	      titlep(ix)='  tau('//char(i+48)//') '
	      titlep(ix+1)='  amp('//char(i+48)//') '
5	   continue
	   titlep(kmax)='  Y(inf) '
	   goto 316
	endif
	if(nmod.eq.20) then
	   xgt0=.true.
	   if(alpha) then
		print 55
55		format(
     &	  '&Number of geometric components to be fitted [1] = ')
		call INPUTi(ncomp)
	   else
		call DEFOLTi(ncomp,defolt)
		call QDIALOG(1,
     &	'Number of geometric components to be fitted',
     &	defolt,11,cans)
		call GETINPi(cans,ncomp)
	   endif
	   if(ncomp.eq.0) ncomp=1
	   kmax=2*ncomp+1		!no of parameters
	   do i=1,ncomp
		ix=2*i-1		!=1,3,5,7..
	      titlep(ix)='  ''tau''('//char(i+48)//') '
	      titlep(ix+1)='  amp('//char(i+48)//') '
	   enddo
	   titlep(kmax)='  Y(inf) '
	   goto 316
	endif
c Power function
	if(nmod.eq.28) then
	   kmax=2		!no of parameters
	   titlep(1)=' x0 (=x at y=ybar) '
	   titlep(2)=' n (=log-log slope)'
	   goto 316
	endif
c
	IF(NMOD.EQ.9.OR.NMOD.EQ.11) print 3141
3141	FORMAT(' CALC OF OP TIME/BST CORRECTED 25-JAN-85')
c
	ip1=1
	ip2=1
	if(iw.ge.1.and.iw.le.5) then
	   if(pprt) write(7,1301) titmod(imod),titw(iw)
         if(discprt) write(8,1301) titmod(imod),titw(iw)
c1301	  format(/,' FITTING:',/,' Equation = ',a60,/,' Weighting = ',a60)
	   titprint=.true.
	endif
	if(nmod.eq.23.or.nmod.eq.24) then
	   if(pprt) write(7,1302)
         if(discprt) write(8,1302)
1302	   format(' (common maximum for all sets fitted)')
	endif
c	goto(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
c     & 19,20,21,22),iabs(nmod)
c NMOD now defined
	if(nmod.eq.21) goto 21
	if(nmod.eq.22) goto 22
	if(nmod.eq.16.or.nmod.eq.17.or.nmod.eq.18) goto 329
	if(nmod.ne.99) goto 981		!not polynomial
	ip2=0		!use to signal disjoint polynomial (ip2.ne.0)
	if(alpha) then
	   print 48
48	   format('&Use a disjoint polynomial [N] ? ')
	   ans='N'
	   call INPUTa(ans)
	   if(UC(ans).eq.'Y') then
   		print 50
50   		format(' Make junction at x1 = ')
		call INPUTr(x1)
	      print 471
471	      format(
     &   '&Degree of polynomial below x1 (1=line, 2=quadratic,..) = ')
	      call INPUTi(ip1)
   		print 51,x1
51		format(' Degree of polynomial above x1 = ',g12.4,' = ')
		call INPUTi(ip2)
	   else	!not disjoint
	      print 47
47	      format(
     &	'&Degree of polynomial (1=line, 2=quadratic,...,10) = ')
	      call INPUTi(ip1)
	   endif
c	If ip2=-1 then Schild eqn is fitted with data Y=dose ratio, X=antag conc
c	If ip2=-2 then Schild eqn is fitted with data Y=dose ratio-1, X=antag conc
	   if(ip1.eq.1) then
		ans='N'
		if(ip2.lt.0) ans='Y'
		call DCASK(
     &	'Fit Schild equation, i.e. parameter=KB (slope=1/KB)',
     &	ans,ans)
		if(ans.eq.'Y') then
		   ans='Y'
		   call DCASK(
     &	   'Are the Y values = (dose ratio-1)',ans,ans)
		   if(ans.eq.'Y') then
			ip2=-2
		   else
			ans='Y'
			call DCASK(
     &		'Are the Y values = (dose ratio)',ans,ans)
			if(ans.eq.'Y') then
			   ip2=-1
			else
			   print 58
58			   format(' Invalid Schild data')
			   goto 64
			endif
		   endif
		endif
	   else if(ip1.eq.2) then
		ans='N'
		if(ip2.lt.0) ans='Y'
		call DCASK(
     &	'Fit mean current-variance parabola, param=N (=1/a2)',
     &	ans,ans)
		iequiv=0		!parameters=a1, a2
		if(ans.eq.'Y') then
		   iequiv=1		!parameters=i, N: fit y=ix - x^2/N
		endif
	   endif
	else
	   ans='N'
	   call DEFOLTa(ans,defolt)
27	   call QDIALOG(1,'Use a disjoint polynomial (F1=help)',
     &	defolt,11,cans)
	   if(help) then				!F1 hit in qdialog
		call CVHELP(2,titmod,imod)
		goto 27
	   else
	      call GETINPa(cans,ans)
	   endif
	   if(UC(ans).eq.'N') then
		ip1=1
		call DEFOLTi(ip1,defolt)
26		call QDIALOG(1,
     &	'Degree of polynomial (F1=help)',
     &	 defolt,11,cans)
		if(help) then				!F1 hit in qdialog
		   call CVHELP(2,titmod,imod)
		   goto 26
		else
	         call GETINPi(cans,ip1)
		endif
	   else		!disjoint
		x1=0.0
		call DEFOLTr(x1,defolt)
		call QDIALOG(1,
     &	  'Make junction at x1',
     &	   defolt,11,cans)
		call GETINPr(cans,x1)
		ip1=1
		call DEFOLTi(ip1,defolt)
261		call QDIALOG(1,
     &	'Degree of polynomial below x1 (F1=help)',
     &	 defolt,11,cans)
		if(help) then				!F1 hit in qdialog
		   call CVHELP(2,titmod,imod)
		   goto 261
		else
	         call GETINPi(cans,ip1)
		endif
		ip2=2
		call DEFOLTi(ip2,defolt)
		call QDIALOG(1,
     &     	  'Degree of polynomial above x1',
     &	   defolt,11,cans)
		call GETINPi(cans,ip2)
	   endif
c	If ip2=-1 then Schild eqn is fitted with data Y=dose ratio, X=antag conc
c	If ip2=-2 then Schild eqn is fitted with data Y=dose ratio-1, X=antag conc
	   if(ip1.eq.1) then
		ans='N'
		if(ip2.lt.0) ans='Y'
		call DEFOLTa(ans,defolt)
271		call QDIALOG(1,
     &	   'Fit Schild equation, i.e. parameter=KB (slope=1/KB)',
     &	   defolt,11,cans)
		if(help) then				!F1 hit in qdialog
		   call CVHELP(2,titmod,imod)
		   goto 271
		else
	   	   call GETINPa(cans,ans)
		endif
		if(ans.eq.'Y') then
		   ans='Y'
		   call DEFOLTa(ans,defolt)
		   call QDIALOG(1,'Are the Y values = (dose ratio-1)',
     &	   defolt,11,cans)
	   	   call GETINPa(cans,ans)
		   if(ans.eq.'Y') then
			ip2=-2
		   else
			ans='Y'
			call DEFOLTa(ans,defolt)
			call QDIALOG(1,'Are the Y values = (dose ratio)',
     &		defolt,11,cans)
	   		call GETINPa(cans,ans)
			if(ans.eq.'Y') then
			   ip2=-1
			else
			   call WDIALOG(1,'Invalid Schild data',12)
			   goto 64
			endif
		   endif
		endif
	   else if(ip1.eq.2) then
		ans='N'
		if(ip2.lt.0) ans='Y'
		call DEFOLTa(ans,defolt)
		call QDIALOG(1,
     &	'Fit mean current-variance parabola, param=N (=1/a2)',
     &	  defolt,11,cans)
	   	call GETINPa(cans,ans)
		iequiv=0		!parameters=a1, a2
		if(ans.eq.'Y') then
		   iequiv=1		!parameters=i, N: fit y=ix - x^2/N
		endif
	   endif
	endif
	if(pprt) write(7,52)ip1
      if(discprt) write(8,52)ip1
52	format(' Polynomial of degree= ',i3,' fitted')
	kmax=ip1+1			!polynomial
c Can fit a different degree of polynomial ( eg a line) above x=x1
c as follows. Coeffs are the same for both polynomials up to coeff
c of X**ip1 (or X**ip2 if ip2<ip1)
	if(ip2.gt.0) then
	   if(ip2.gt.ip1) kmax=ip2+1	!to get right no of params
	   if(pprt) write(7,521)x1,ip2
         if(discprt) write(8,521)x1,ip2
521	   format(
     & ' Above x= ',g12.4,' polynomial of degree= ',i3,' fitted')
	else if(ip2.lt.0) then
	   if(pprt) write(7,522)
         if(discprt) write(8,522)
522	   format('  (Fit of Schild equation with KB as parameter)')
	else if(ip1.eq.2.and.iequiv.eq.1) then
	   if(pprt) write(7,523)
         if(discprt) write(8,523)
523	   format('  Fit of mean current-variance parabola')
	endif
	titlep(1)=t1(14)		!Y(0) (=a0)
	do k=2,kmax
	   call INTCONV(k-1,cnum)
c==	   TITLEP(k)='    '//tt1(k)	!b1,b2,...
	   TITLEP(k)='    a'//charnb(cnum)	!a1,a2,...
	enddo
	if(ip2.eq.-1.or.ip2.eq.-2) then	!fit Y=x/KB or Y=1+x/KB
	   kmax=1
	   titlep(1)='                 KB '
	endif
	if(ip1.eq.2.and.iequiv.eq.1) then	!current variance parabola
	   kmax=2
	   titlep(1)='single chan current '
	   titlep(2)='                  N '
	endif
c  Lines with common X intercept (eg reversal potential).  For n lines
c must estimate n slopes + 1 intercept=n+1 parameters
c (only if ifitmode=3)
	if(ip1.eq.1.and.ip2.eq.0.and.nsfit.gt.1.and.ifitmode.eq.3) then
	   if(alpha) then
	      print 40
40	      format(
     &	' Constrain lines to common intercept on x axis [N] ? ')
		ans='N'
		call INPUTa(ans)
	   else
	      ans='N'
	      call DEFOLTa(ans,defolt)
	      call QDIALOG(1,
     &	'Constrain lines to common intercept on x axis',
     &	 defolt,11,cans)
	      call GETINPa(cans,ans)
	   endif
	   if(UC(ans).eq.'Y') then
		nmod=25		!new model #
		kmax=nsfit+1
		do i=1,nsfit
		   call INTCONV(i,cnum)
		   n=nblank(cnum)
		   titlep(i)='slope('//cnum(1:n)//')'	!slope(i)
		enddo
		titlep(kmax)='X (at Y=0)'
	   endif
	endif
	GOTO 316	!MODEL DEFINED-START FITTING
c
981	continue
	if(nmod.eq.13) goto 315
	if(nmod.eq.14.or.nmod.eq.15) goto 424
	if(nmod.ne.10.and.nmod.ne.11) goto 350
351	continue
c  (1) mods 10,11 IP1=1 if param are K1,K2; =2 if param are K1, K2/K1.
c  (2) mods 9,11  IP2=1 if param is K2 or (K2/K1 if IP1=2)
c		  IP2=2 if param is k(-2) (must have ip2=1 if ip1=2)
c  (3) mods 9,11  IP2 set negative if k(+2) is a param  to be estimated
	if(alpha) then
	   print 3184
3184	   format(' Use K1,K2 (Y), [or K1, K2/K1 (N)] as parameters ? ')
	   ans='N'
	   call INPUTa(ans)
	   if(UC(ans).eq.'Y') ip1=1
	   if(UC(ans).ne.'Y') ip1=2
	   print 3183
3183	   format(' Fit for non-equiv subunits (no coop) [N] ? ')
C YES MEANS K1,K2 ARE FOR SUBUNITS print 1,2- ASSUMED SAME
C REGARDLESS OF ORDER OF BINDING (NO COOP). ALSO FORWARD RATE
C FOR NMOD=11 TAKEN AS SAME FOR BOTH SUBUNITS=X2
	   ans='N'
	   call INPUTa(ans)
	else
	   call DEFOLTi(ip1,defolt)
	   call QDIALOG(1,
     &     '(1) Use K1, K2; (2) use K1, K2/K1 as parameters: n',
     &	defolt,11,cans)
	   call GETINPi(cans,ip1)
	   if(ip1.ne.1.and.ip1.ne.2) goto 351
	   ans='N'
	   call DEFOLTa(ans,defolt)
	   call QDIALOG(1,'Fit for non-equiv subunits (no cooperativity)',
     &     defolt,11,cans)
	   call GETINPa(cans,ans)
	endif
	if(UC(ans).ne.'Y'.and.UC(ans).ne.'N') goto 351
	if(UC(ans).eq.'Y') goto 352
c
c	CASTAR=-1.
	iequiv=0
      if(pon()) write(7,353)
      if(discprt) write(8,353)
353	format(' Cooperative binding  (equivalent sites)',/)
	goto 350
c
352	iequiv=1
c352	CASTAR=1.		!USE CASTAR AS FLAG
      if(pon()) write(7,354)
      if(discprt) write(8,354)
354	format(' Nonequivalent sites (non-cooperative binding)',/)
350	continue
c
	IF(NMOD.GT.4) GOTO 300
c32	continue
c Next bit never done now -fit of line done under nmod=26,27 now
	kmax=nmod+1		!for nmod=1-4
	ncomp=kmax/2	!for nmod=1-4
	if(nmod.eq.-1) then
	   ncomp=1
	   kmax=4		!Hill equation
	endif
	FLINE=.false.
	if(nmod.eq.2.or.nmod.eq.4) fline=.true.	!for nmod=2,4 ie imod=4,6
C
	titlep(1)=t1(1)
	titlep(2)=t1(2)
	titlep(3)=t1(5)
	IF(NMOD.NE.-1) GOTO 119
	titlep(3)=t1(14)
	titlep(4)=t1(15)
119	IF(NCOMP.NE.2) GOTO 316
	titlep(3)=t1(3)
	titlep(4)=t1(4)
	titlep(5)=t1(5)
	GOTO 316
C
300	CONTINUE
	if(nmod.eq.5) then	!binding inhibition -1 langmuir only now
	   ncomp=1		!only possibility
	   if(alpha) then
	      print 318
318		format(
     &  	 ' Normalised conc of labelled ligand = ')
	      call INPUTr(castar)
	   else
		castar=1.
		call DEFOLTr(castar,defolt)
29		call QDIALOG(1,
     &	'Normalized conc of labelled ligand (F1=help)',
     &	defolt,11,cans)
	      if(help) then				!F1 hit in qdialog
		   call CVHELP(3,titmod,imod)
		   goto 29
	      else
		   call GETINPr(cans,castar)
		endif
	   endif
	   kmax=3
	   titlep(1)='           Y(inf)   '
	   titlep(2)='           Y(0)     '
	   titlep(3)='           KB       '
	   GOTO 316	!MODEL DEFINED-START FITTING
	endif
C
506	IF(NMOD.NE.6.AND.NMOD.NE.7) GOTO 508
	KMAX=3
	IF(NMOD.EQ.7) KMAX=2+NSET	!K1,K2 + NSET VALUES OF BMAX
	titlep(1)=t1(8)
	titlep(2)=t1(9)
	titlep(3)=t1(7)
	GOTO 316
C
508	CONTINUE	! MODELS 8-12 HERE
	KMAX=3
	IF(NMOD.EQ.10.OR.NMOD.EQ.11) GOTO 510
	titlep(1)=t1(11)
	IF(NMOD.EQ.9) titlep(1)=t1(12)
	titlep(2)=t1(1)
	titlep(3)=t1(13)
	if(nmod.eq.12) titlep(4)=t1(2)
	if(nmod.eq.12) KMAX=4
	GOTO 3161
510	KMAX=4
	IF(NMOD.EQ.10) titlep(1)=t1(11)
	IF(NMOD.EQ.11) titlep(1)=t1(12)
	titlep(2)=t1(1)
	if(ip1.eq.1) titlep(3)=t1(3)
	if(ip1.eq.2) titlep(3)=t1(16)
	titlep(4)=t1(13)
c
c3161	IF(NMOD.EQ.8.OR.NMOD.EQ.10.or.nmod.eq.12) GOTO 316
3161	if(nmod.ne.9.and.nmod.ne.11) goto 316
c Section for NMOD=9 or 11
c  (1) mods 10,11 IP1=1 if param are K1,K2; =2 if param are K1, K2/K1.
c  (2) mods 9,11  IP2=1 if param is K2 or (K2/K1 if IP1=2)
c		  IP2=2 if param is k(-2) (must have ip2=1 if ip1=2)
c  (3) mods 9,11  IP2 set negative if k(+2) is a param  to be estimated
c==============graphics not done from here on
	print 3181
3181	FORMAT(' Mean burst length (ms)= ')
	call INPUTr(x1)
	if(pprt) write(7,3188)x1
      if(discprt) write(8,3188)x1
3188	format(' Mean burst length (ms)= ',g13.6)
	X1=X1/1000.	!SEC
	ip2=1
	if(ip1.eq.2) goto 3191
	print 3190
3190	format(' Use K2 (Y), [or k(-2), (N)] as parameter [Y] ? ')
	ans='Y'
	call INPUTa(ans)
	if(UC(ans).eq.'Y') ip2=1
	if(UC(ans).ne.'Y') ip2=2
	if(ip2.eq.1) goto 3191
	i1=2
	if(nmod.eq.11) i1=3
	titlep(i1)=t1(17)
3191	print 3185
3185	format(' Treat k(+2) as a parameter to be estimated [Y] ? ')
	ans='Y'
	call INPUTa(ans)
	if(UC(ans).eq.'Y') ip2=-ip2
	if(ip2.gt.0) goto 3186
	kmax=kmax+1
	titlep(kmax)=t1(18)
	goto 316
3186	print 3182
3182	FORMAT(' k(+2)= ')
	call INPUTr(x2)
	if(pprt) write(7,3189)x2
      if(discprt) write(8,3189)x2
3189	format(' k(+2) fixed at ',g13.6)
	GOTO 316
C
c model 13
315	print 323
323	format(' Same affinity for both subunits [Y] ? ')
	ans='Y'
	call INPUTa(ans)
	ip1=0
	if(UC(ans).eq.'Y') ip1=1
	print 324
324	format(' Same beta/alpha for both subunits [Y] ? ')
	ans='Y'
	call INPUTa(ans)
	ip2=0
	if(UC(ans).eq.'Y') ip2=1
	titlep(1)=t1(11)
	titlep(2)=t1(1)
	kmax=2
	if(ip1.eq.1) goto 325
	kmax=kmax+1
	titlep(kmax)=t1(3)
325	if(ip2.eq.1) goto 326
	kmax=kmax+1
	titlep(kmax)=t1(19)
326	kmax=kmax+1
	titlep(kmax)=t1(13)
	goto 316
c
c Model 14 & 15
c Kij= ith subunit type, jth (sequential) binding
424	kmax=5
	titlep(1)=t1(1)
	titlep(2)=t1(20)
	titlep(3)=t1(21)
	if(nmod.eq.14) titlep(4)=t1(11)
	if(nmod.eq.15) titlep(4)=t1(12)
	titlep(5)=t1(13)
	if(nmod.eq.14) GOTO 316
c
	print 3181	!' Mean burst length (ms)= ')
	call INPUTr(x1)
	if(pprt) write(7,3188)x1		!Mean burst lengt
      if(discprt) write(8,3188)x1		!Mean burst length (ms)= ,rec= h (ms)= '
	X1=X1/1000.	!SEC
	print 428	!Treat k(+2) as a parameter to be estimated?
428	format(
     & ' Treat k(+2) (assumed to be same for both subunit types)',/,
     & ' as a parameter to be estimated [N] ? ')
	ans='N'
	call INPUTa(ans)
	x2=-1.		!k+2 not to be estimated
	if(UC(ans).eq.'N') goto 426
	kmax=kmax+1
	titlep(kmax)=t1(18)
	goto 316
426	print 3182	!' k(+2)= '
	call INPUTr(x2)
	if(pprt) write(7,429)x2
      if(discprt) write(8,429)x2
429	format(' k(+2) fixed at ',g13.6,/,
     & '   (assumed to be same for both subunit types)')
	goto 316
c
c Model 16
329	continue
	kmax=3
	titlep(1)=t1(22)
	titlep(2)=t1(23)
	titlep(3)=t1(24)
	if(nmod.eq.16) goto 316
	kmax=5		!nmod=17
	titlep(2)=t1(26)
	titlep(4)=t1(27)
	titlep(5)=t1(25)
	if(nmod.eq.18) titlep(5)=t1(11)
	GOTO 316	!MODEL DEFINED-START FITTING
c
c Model 21 ' (21) I/V for GHK with Mg block'
21	continue
	kmax=7
c29	gNa=theta(1)		!conductance in pS
c30	gCs=theta(2)		!conductance in pS
c31	Eoff=theta(3)	!offset in Erev, relative to GHK value
c1	aK1=theta(4)		!in mM
c3	aK2=theta(5)		!in mM
c32	H1=theta(6)
c33	H2=theta(7)
	titlep(1)=t1(29)
	titlep(2)=t1(30)
	titlep(3)=t1(31)
	titlep(4)=t1(1)
	titlep(5)=t1(3)
	titlep(6)=t1(32)
	titlep(7)=t1(33)
	goto 316
c Model 22 ' (22) Exponential voltage dependence of rate constant'
22	continue
	kmax=2
	titlep(1)=t1(28)  !tau
	titlep(2)=t1(24)  !H (mV)
	goto 316
c
c
C Model now defined- start fitting
316	continue
	if(iw.ge.1.and.iw.le.5.and.(.not.titprint)) then
	   if(pprt) write(7,1301) titmod(imod),titw(iw)
         if(discprt) write(8,1301)titmod(imod),titw(iw)
1301	  format(/,' FITTING:',/,' Equation = ',a60,/,' Weighting = ',a60)
	   if(nmod.eq.30.and.iequiv.eq.0) then
	      if(pprt) write(7,42)
            if(discprt) write(8,42)
42	  	format(
     &	' Assume same (microscopic) binding constant for each step')
	   else if(nmod.eq.30.and.iequiv.eq.0) then
	      if(pprt) write(7,43)
            if(discprt) write(8,43)
43	  	format(
     &	'    Cooperativity allowed in binding: index r=K(i-1)/K(i)')
	   endif
	   if(pprt) write(7,108)
         if(discprt) write(8,108)
108	   format(
     &   '===========================================',/)
	endif
c
c Before leaving add extra params for ifitmode=4 (nset-1 potency ratios)
c and for ifitmode=5 (KB)
c For ifitmode=4 jset=set # for the 1st set to be fitted so this is reference
c set for other potency ratios; title of pot ratios = 'rn' where n=data set
c number (n>jset) to which the ratio applies
	kmax1=kmax		!keep number of params before adding extra ones!
	if(ifitmode.eq.4) then
	   i=0
	   do j=jset+1,nset
c	   do 317 j=jset+1,nset
c	    if(jmiss(j).eq.1) goto 317	!set skipped
	    if(jmiss(j).eq.0) then	!set not skipped
	      i=i+1
	      titlep(kmax+i)='       r'//char(j+48)//' '	!='r3','r5',...
	      if(j.eq.10) titlep(kmax+i)='      r10 '
	    endif
	   enddo
c317	   continue
c	   kmax=kmax+i                      !new kmax for simult fit
	   kmax=kmax+nsfit-1                !new kmax for simult fit
	else if(ifitmode.eq.5) then
	   kmax=kmax+1
	   titlep(kmax)=t1(13)		!='       KB '
	endif
c	call flush(7)
	RETURN
	END


