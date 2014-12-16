	subroutine LIKINT1(ELFUNC,noerr,sdev,nmod,constr,nfit,
     & theta,jfix,titlep,ndth,Xobs,yobs,w,nj,juse,setx,niobs,njset)
c LIKINT1 (06/25/92 02:44pm) is version of LIKINT0 (for CVFIT at present)
c that uses SD from HESMAT (if they have been calc) to improve automatic
c initial guesses for likelihood intervals
c
c 04/27/97 05:40pm Many modifs, mainly designed to cope with cases
c where simplex does not converge.  Convergence is helped by use
c of stored initial guesses: problem seems often to be that
c setting init guesses to their ML values (as before) provides poor
c guesses for simplex when one param is fixed at limit value. In bisecv
c now have option to scan a range of limit values separate from bisection
c and each value uses as init guess the params output by the last simplex
c run that converged OK. During this scan, keep the values that were
c closest to the required elfunc value of 0 and use them for initial guess
c when guesses reset. Allocate thetsav(1,j) and xlimsav(j) to store
c output from every convereged fit i=1,...,kmax, j=1,..nfits (eg nfits=200
c to store results of 200 fits) -then when a new putative limit is to be
c tested, look through xlimsav() for closest value that has already been
c tried and use corrsponding thetsav() as initial guess for the non-fixed params
c
c 01/12/96 03:13pm :theta() out of common and fixed for dimension=ndth
c (=100 at present) (also for jfix,sdev,titlep)
c 12/28/89 04:23pm Lahey version of general purpose subroutine
c for calc of likelihood intervals. The subroutine name for calc
c of the likelihood for specified parameters is given in call.
c Version for 20 parameters.
c
	dimension XLL(2,50),XUL(2,50),NPFIL(50)
	real THETA(ndth)
	ALLOCATABLE::thetsav,xlimsav,youtsav
	real*4 thetsav(:,:),xlimsav(:),youtsav(:)
	real*4 sdev(ndth),se(50)
	integer JFIX(ndth)
	character*20 TITLEP(ndth)		!names of params for fit
c	character*1 UC,ans
	logical valid	!true if simplex converged
	LOGICAL DEBUG,LOWER,UPPER,omit,noerr,constr
	logical discprt,slock,caplock,pon,prt,deb,allocated
	common/dp/discprt
	COMMON/LIKBLK/kmax,kfit,elmax,em,delmin,ndisp,npar,iw,sres
c Params to be passed to SSDCV and Simplex
	real*4 xobs(niobs,njset),yobs(niobs,njset),w(niobs,njset)
	integer juse(njset),nj(njset)
	real setx(njset)
	common/simp/iconv,nevmax		!from elfuncv/simplexv
	EXTERNAL ELFUNC
C
C Note that if KFIT=1, i.e. only one param to be estimated then
c when this is fixed at a trial value there are none to be
c estimated so SIMPLEX not needed- just call ELFUNC
c
	pon()=slock()
	debug()=caplock()
c
101	format(a1)
	iopt1=1
	prt=pon()
	deb=debug()
	nfix=0
	do 40 i=1,kmax
	it=i
	if(jfix(i).eq.1) nfix=nfix+1
	if(nmod.eq.7.and.i.gt.3) it=3	!repeat BMAX title for i=3,4,...
	if(jfix(i).eq.0.and.(.not.noerr)) then    !errors calc, not fixed
	   ii=ii+1
	   se(i)=Sdev(i)	!NB Sdev(i) set to -1 for fixed param
c In hesmat3, sres is already incorporated into sdev(i) when residuals used
c for errors (iw=1,5)
c===	   if(iw.eq.1.or.iw.eq.5) se(i)=sres*se(i)
	   print 41,I,TITLEP(IT),THETA(I),se(i)
c	   if(prt) write(7,41)I,TITLEP(it),THETA(I),se(i)
c	   if(discprt) write(8,41)I,TITLEP(it),THETA(I),se(i)
41	   FORMAT(' Param. ',I2,': ',a20,'= ',G13.6,' Approx SD= ',
     &	g13.6)
	else if(jfix(i).eq.1) then	!fixed (whether errors calc or not)
411	   print 412,I,TITLEP(it),THETA(I)
c	   if(prt) write(7,412)I,TITLEP(it),THETA(I)
c         if(discprt) write(8,412)I,TITLEP(it),THETA(I)
412	   FORMAT(' Param. ',I2,': ',a20,'= ',G13.6,' (fixed)')
	else if(noerr.and.jfix(i).eq.0) then		!no errors,not fixed
413	   print 414,I,TITLEP(it),THETA(I)
c	   if(prt) write(7,414)I,TITLEP(it),THETA(I)
c         if(discprt) write(8,414)I,TITLEP(it),THETA(I)
414	   FORMAT(' Param. ',I2,': ',a20,'= ',G13.6)
	endif
40	continue
c
	ndf=nfit-kfit
	print 17,ndf,iopt1
17	FORMAT(
     & ' Compute m-unit intervals:',/,
     & '  (1) Specify m (0.5 for ''1 SD'', 2 for ''2 SD'')',/,
     & '  (2) Use value of m for P=0.95, with ',i3,' d.f.',/,
     & ' Option number [',i2,'] = ')
	call INPUTi(iopt1)
	if(iopt1.eq.1) then
	   em=0.5
	   print 171
171	   format(' m value [0.5] = ')
	   call INPUtr(em)
	else if(iopt1.eq.2) then
	   call TVALUE(ndf,tval)
	   em=tval*tval/2.0
	endif
c	call INPUTr(em)
	clim=sqrt(2.*em)		!equiv SD for Gaussian
	ELCRIT=ELMAX-EM
C INITIALISE
	do 49 i=1,2
	do 49 j=1,20
	XLL(I,J)=-2.E30
49	XUL(I,J)=-2.E30
C
32	print 42
42	FORMAT(
     & '&If calculation does not converge, hit F2 to stop',/,/,
     & ' No of parameters for which limits to be calc [0 to skip] = ')
	np=0
	call INPUTi(np)
	IF(NP.LE.0) GOTO 99
	print 54,em,clim
	if(prt) write(7,54) em,clim
	if(discprt) write(8,54) em,clim
54	format(/,
     & ' LIKELIHOOD INTERVALS',/,
     & ' ',f8.2,'-unit Likelihood Intervals',/,
     & ' (like ',f8.2,'-SD intervals if Gaussian)',/)
c
	DO 43 J=1,NP
131	print 13,j
13	FORMAT('& ',i3,': parameter number = ')
	call INPUTi(i)
	if(i.lt.1.or.i.gt.kmax) goto 131
	NPFIL(J)=i
c Defaults may well not be good for parameters that are neg, or for those
c where limits may include zero
	xhigh=theta(i)
	if(noerr.or.se(i).lt.0.01*theta(i)) then
	   if(theta(i).lt.0.) then
		xlow=5.*theta(i)	!more negative
	   else
		xlow=0.2*theta(i)
	   endif
	else
	   xlow=theta(i) - 2.0*clim*se(i)		!leave 2-fold margin
	endif
	if(constr.and.xlow.lt.1.e-20) xlow=0.05*abs(theta(i))	!when constrained pos
c     NB cannot have guess=0 when stpfac/errfac used in simplex
	if(abs(xlow).lt.0.1*abs(theta(i))) then
	   if(constr) then
		xlow=0.2*abs(theta(i))
	   else
		if(theta(i).lt.0.) then
		   xlow=5.*theta(i)	!more negative
		else
		   xlow=0.2*theta(i)
		endif
	   endif
	endif
	if(xlow.gt.xhigh) then
	   p=xlow
	   xlow=xhigh
	   xhigh=p
	endif
	XLL(1,i)=XLOW
	XLL(2,i)=XHIGH
47	continue
	xlow=theta(i)
	if(noerr.or.se(i).lt.0.01*theta(i)) then
	   xhigh=15.0*abs(theta(i))
	else
	   xhigh=theta(i) + 5.0*clim*se(i)		!leave 5-fold margin
	endif
	if(xlow.gt.xhigh) then
	   p=xlow
	   xlow=xhigh
	   xhigh=p
	endif
	XUL(1,i)=XLOW
	XUL(2,i)=XHIGH
43	continue          !end of param loop
C
C
C PRECISION IN simplex CONTROLLED BY DELMIN=(CHANGE IN MIN)
C SAY ABSOLUTE ERROR OF 0.002 IN L(THETA) IS CLOSE ENOUGH THEN
	nevmax=5000		!max number of iterations for simplex
	ioptsav=1		!initially, for reset of nevmax in bisecv
	omit=.false.	!so BISEC tests and prints init guesses
c	DELMIN=0.005		!defaults used if debug()=false
	delmin=-1.
	yerr=0.01		!max deviation from Lcrit
	NEPSX=1000
	epsy=-1.		!use nstep to control bisection
c	epsy=0.002
	ndisp=-1000		!no print out at all in Simplex
	if(deb) then
	   if(kfit.ne.1) then
		print 5
5 		FORMAT(' Abs accuracy for L(theta) in SIMPLEX: delmin = ')
		call INPUTr(del)
		DELMIN=DEL
		print 3
3		FORMAT(' print every Nth iteration in SIMPLEX. N = ')
		call INPUTi(ndisp)
	   endif
	   print 28
28	   FORMAT(
     & ' Accuracy for param (in Bisection) set to 1 in N. N= ')
	   call INPUTi(nepsx)
	   print 56
56	   FORMAT(
     &   ' Abs accuracy L-Lmax in bisection (ignore if neg) = ')
	   call INPUTr(epsy)	!ignored if neg
	endif
	if(epsy.lt.0.) then
	   print 541,nepsx
	   if(prt) write(7,541) nepsx
	   if(discprt) write(8,541) nepsx
541	   format('& Max error in limit = (ML value)/',i5,/)
	else if(epsy.gt.0.) then
	   print 542,epsy
	   if(prt) write(7,542) epsy
	   if(discprt) write(8,542) epsy
542	   format(' Max error in L-Lmax = ',g13.6)
	endif
C
	nfmax=200	!max # of fits to be stored
	ALLOCATE(thetsav(kmax,nfmax),xlimsav(nfmax),youtsav(nfmax))
	nfits=0	!none stored yet

C NOW START LOOP FOR THE PARAMETERS TO BE DONE
	DO 48 IP1=1,NP
	NPAR=NPFIL(IP1)	!CURRENT PARAMETER NUMBER
	ip=npar
      print 481,npar,titlep(npar),theta(npar)
      if(pon()) write(7,481) npar,titlep(npar),theta(npar)
      if(discprt) write(8,481) npar,titlep(npar),theta(npar)
481	format(' Parameter',i3,': ',a20,' = ',g13.6)
	JFIX(NPAR)=1	!FIX THIS PARAMETER
c Insert check
	call CHECKfix(kmax,jfix,nfix,ndth)
c At first the only available fit is the ML fit itself, so store it
	nfits=1
	do i=1,kmax
	   thetsav(i,1)=theta(i)
	enddo
	xlimsav(1)=theta(npar)
	youtsav(1)=em
c
	EPSX=THETA(NPAR)/FLOAT(NEPSX)	!ABS ACCURACY FOR BISECT
	LOWER=.FALSE.
	UPPER=.FALSE.
	IF(XLL(1,IP).GT.-1.E30) LOWER=.TRUE.	!CALC LOWER LIMIT
	IF(XUL(1,IP).GT.-1.E30) UPPER=.TRUE.	!CALC UPPER LIMIT
C NOW CALC LOWER LIMIT FIRST IF WANTED
553	continue
	IF(.NOT.LOWER) GOTO 50
	XLOW=XLL(1,IP)
	XHIGH=XLL(2,IP)
	if(xlow.gt.xhigh) then
	   x=xhigh
	   xhigh=xlow
	   xlow=x
	endif
	ND1=1	!PRINT EVERY CYCLE IN BISECT
	IF(deb) print 55,XLOW,XHIGH
55	FORMAT( ' Initial guess= ',G13.6,2X,G13.6,/)
	call CHECKfix(kmax,jfix,nfix,ndth)
	call CLRKB()
	call BISECV(ELFUNC,xlow,xhigh,0.0,Xlim,Yout,epsx,epsy,nerr,
     & nd1,omit,ioptsav,kmax,yerr,xlo,xhi,Xobs,yobs,w,nj,juse,
     & setx,niobs,njset,jfix,ndth,thetsav,xlimsav,youtsav,
     & nfits,nfmax,npar,valid)
	if((iabs(nerr).ge.1.and.iabs(nerr).le.4).or.
     & (abs(yout).gt.yerr)) then
	   call BELL(2)
	   if(iabs(nerr).lt.3) then
	      print 551
551	      format(' Guesses for lower limit do not include result:')
	   endif
	   if(iabs(nerr).eq.4) then
	      print 5511
5511	      format(' One or more Simplex runs did not converge:')
	   endif
	   if(nerr.ne.0) then
	     print 71,ip,theta(ip),nerr
71	     format(
     &   ' Param #',i3,'=',g12.5,': ERROR #',i2,' in LOWER limit')
	   else
	     print 711,ip,theta(ip)
711	     format(
     &   ' Param #',i3,'=',g12.5,': ERROR too big in LOWER limit')
	   endif
c     If error in limit, then just check all so far in case good enough
	   ym=1.e37
	   jmin=0
	   do j=1,nfits
		a=abs(youtsav(j))
   		x=xlimsav(j)
		if(a.lt.ym.and.x.lt.theta(ip)) then	!lower lim only
		   ym=a
		   jmin=j
		endif
	   enddo
	   if(jmin.gt.0) then
		xlim=xlimsav(jmin)
		yout=ym
		print 72,xlim,elmax,elcrit,
     &  	 elcrit+yout,yout
72	   	format(
     &    '  Closest estimate of lower limit found so far = ',g12.5,/,
     &    '  [L(max)= ',G13.6,'  L(crit)= ',G13.6,/,
     &    '  L(limit)= ',g13.6,' Dev from L(crit)= ',F10.5,']')
	   endif
	   iopt=1
	   if(abs(yout).lt.2.*yerr.and.jmin.gt.0) iopt=2   !default=accept
	   print 73,iopt
73	   format(
     & ' (1) give new guesses for lower limit and try again',/,
     & ' (2) accept this limit anyway',/,
     & ' (3) go on to next limit',/,
     & ' (4) abandon likelihood interval calculations',/,
     & ' Option number [',i2,'] = ')
	   call INPUTi(iopt)
	   if(iopt.eq.1) then
c		xlow=xll(1,ip)
c		xhigh=xll(2,ip)
		xlow=xlo
		xhigh=xhi
	      print 15,xlow,xhigh
15	      FORMAT('&  New guesses [',g11.4,',',g11.4,'] = ')
	      call INPUT2r(xlow,xhigh)
		if(xlow.gt.xhigh) then
		   p=xlow
		   xlow=xhigh
		   xhigh=p
		endif
		XLL(1,ip)=xlow
		XLL(2,ip)=xhigh
		goto 553
	   endif
	   if(iopt.ne.2) then	!for iopt=2 go straight on
	   	if(iopt.eq.4) jfix(npar)=0	!UNFIX this parameter if abandoned
		if(iopt.eq.3) goto 50
		if(iabs(nerr).eq.3) goto 992	!aborted
		goto 991				!iopt=4
	   endif
	endif
	print 29,EM,XLIM,xlim-theta(npar)
	if(prt) write(7,29) EM,XLIM,xlim-theta(npar)
      if(discprt) write(8,29) EM,XLIM,xlim-theta(npar)
29	FORMAT(
     & '   LOWER ',F6.2,'-unit interval= ',G13.6,' ( dev= ',g13.6,')')
	print 361,ELMAX,ELCRIT,elcrit+yout,YOUT
	if(prt) write(7,361) ELMAX,ELCRIT,elcrit+yout,YOUT
      if(discprt) write(8,361) ELMAX,ELCRIT,elcrit+yout,YOUT
361	FORMAT(
     & '           [L(max)= ',G13.6,'  L(crit)= ',G13.6,/,
     & '            L(limit)= ',g13.6,' Dev from L(crit)= ',F10.5,']',/)
	if(prt) call flush(7)
c	if(upper) then
c	   print 362
c362	   format(' Hit any key to calculate upper limit')
c	   call ANYKEY
c	endif
C
C LOWER LIMIT DONE (IF REQ). NEXT DO UPPER
50	continue			!jump here to do upper limit
	IF(.NOT.UPPER) GOTO 51		!NEXT PARAMETER
555	continue
	print 77,ip,theta(ip)
77	format(' Next do upper limit for parameter ',i3,' = ',g12.5)
	XLOW=XUL(1,IP)
	XHIGH=XUL(2,IP)
	if(xlow.gt.xhigh) then
	   x=xhigh
	   xhigh=xlow
	   xlow=x
	endif
	ND1=1	!PRINT EVERY CYCLE IN BISECT
	IF(deb) print 55,XLOW,XHIGH
	call CHECKfix(kmax,jfix,nfix,ndth)
	call CLRKB()
	call BISECV(ELFUNC,xlow,xhigh,0.0,Xlim,Yout,epsx,epsy,nerr,
     & nd1,omit,ioptsav,kmax,yerr,xlo,xhi,Xobs,yobs,w,nj,juse,
     & setx,niobs,njset,jfix,ndth,thetsav,xlimsav,youtsav,
     & nfits,nfmax,npar,valid)
	if((iabs(nerr).ge.1.and.iabs(nerr).le.4).or.
     & (abs(yout).gt.yerr)) then
	   call BELL(2)
	   if(iabs(nerr).lt.3) then
	      print 554
554	      format(' Guesses for upper limit do not include result:')
	   endif
	   if(iabs(nerr).eq.4) then
	      print 5511
c5511	      format(' One or more Simplex runs did not converge:')
	   endif
	   if(nerr.ne.0) then
	     print 74,ip,theta(ip),nerr
74	     format(
     &   ' Param #',i3,'=',g12.5,': ERROR #',i2,' in UPPER limit')
	   else
	     print 741,ip,theta(ip)
741	     format(
     &   ' Param #',i3,'=',g12.5,': ERROR too big in UPPER limit')
	   endif
c     If error in limit, then just check all so far in case good enough
	   ym=1.e37
	   jmin=0
	   do j=1,nfits
		a=abs(youtsav(j))
   		x=xlimsav(j)
		if(a.lt.ym.and.x.gt.theta(ip)) then	!upper lim only
		   ym=a
		   jmin=j
		endif
	   enddo
	   if(jmin.gt.0) then
		xlim=xlimsav(jmin)
		yout=ym
		print 75,xlim,elmax,elcrit,
     &  	 elcrit+yout,yout
75	   	format(
     &   '  Closest estimate of upper limit found so far = ',g12.5,/,
     &   '  [L(max)= ',G13.6,'  L(crit)= ',G13.6,/,
     &   '  L(limit)= ',g13.6,' Dev from L(crit)= ',F10.5,']')
	   endif
	   iopt=1
	   if(abs(yout).lt.2.*yerr.and.jmin.gt.0) iopt=2   !default=accept
	   print 76,iopt
76	   format(
     & ' (1) give new guesses for upper limit and try again',/,
     & ' (2) accept this limit anyway',/,
     & ' (3) go on to next limit',/,
     & ' (4) abandon likelihood interval calculations',/,
     & ' Option number [',i2,'] = ')
	   iopt=1
	   call INPUTi(iopt)
	   if(iopt.eq.1) then
		xlow=xlo
		xhigh=xhi
	      print 33,xlow,xhigh
33	      FORMAT('&  New guesses [',g11.4,',',g11.4,'] = ')
	      call INPUT2r(xlow,xhigh)
		if(xlow.gt.xhigh) then
		   p=xlow
		   xlow=xhigh
		   xhigh=p
		endif
		XUL(1,ip)=xlow
		XUL(2,ip)=xhigh
		goto 555
	   endif
	   if(iopt.ne.2) then	!for iopt=2 go straight on
	   	if(iopt.eq.4) jfix(npar)=0	!UNFIX this parameter if abandoned
		if(iopt.eq.3) goto 51
		if(iabs(nerr).eq.3) goto 992	!aborted
		goto 991				!iopt=4
	   endif
	endif
	print 37,em,XLIM,xlim-theta(npar)
	if(prt) write(7,37) em,XLIM,xlim-theta(npar)
      if(discprt) write(8,37) em,XLIM,xlim-theta(npar)
37	FORMAT(
     & '   UPPER ',F6.2,'-unit interval= ',G13.6,' ( dev= ',g13.6,')')
	print 361,ELMAX,ELCRIT,elcrit+yout,YOUT
	if(prt) write(7,361) ELMAX,ELCRIT,elcrit+yout,YOUT
      if(discprt) write(8,361) ELMAX,ELCRIT,elcrit+yout,YOUT
c
51	JFIX(npar)=0	!UNFIX THIS PARAM BEFORE GOING TO NEXT
48	continue		!end of loop for each parameter
c
	goto 99		!finished
c
c Error returns
992	print 30
30	format(' Bisection aborted')
991	continue
99	CONTINUE
	if(prt) call flush(7)
	if(allocated(thetsav)) DEALLOCATE(thetsav,xlimsav)
	RETURN
	END

	subroutine CHECKfix(kmax,jfix,nfix,ndth)
	integer JFIX(ndth)
c
	n=0
	do m=1,kmax
	  if(jfix(m).eq.1) n=n+1
	enddo
	if(n.ne.nfix+1) then
	   call BELL(2)
	   print 70,n,nfix+1
70	   format(' Error in LIKINT1:',i2,' fixed param: should be ',i2)
	endif
	RETURN
	end

