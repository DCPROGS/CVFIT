	subroutine BISECV(ELFUNC,x1,x2,y,Xout,Yout,epsx,epsy,nerr,
     & ndisp,omit,ioptsav,kmax,yerr,xlo,xhi,Xobs,yobs,w,nj,juse,
     & setx,niobs,njset,jfix,ndth,thetsav,xlimsav,youtsav,
     & nfits,nfmax,npar,valid)
c
c NB when, after convergence failure, we use 'step through range
c of values' to get idea of where new guesses should be, it may help
c if we always start near to the ML value, using ML estimates as guesses
c for non-fixed params in Simplex, and move away from this value
c progressively so guesses for non-fixed values don't change too suddenly
c If we start at ML value of the fixed param then the selection
c method at start of ELFUNCV will choose theta (stored in thetsav(1)) as
c the initial guess.
	real*4 thetsav(kmax,nfmax),xlimsav(nfmax),youtsav(nfmax)
	integer JFIX(ndth)
	logical valid	!true if simplex converged
c Params to be passed to SSDCV and Simplex
	real*4 xobs(niobs,njset),yobs(niobs,njset),w(niobs,njset)
	integer juse(njset),nj(njset)
	real setx(njset)
	common/simp/iconv,nevmax	!so convergence testable in likint1
c Stuff for BISEC
	LOGICAL OMIT
c	integer*2 ktype
	logical KBHIT		!Utility lib
	character ch,getch

	external ELFUNC
C
c BISECV is special version of BISEC0 which passes necessary parameters
c for ELFUNCv, in CVFIT (3rd line onwards of param)
c
c Solves func(x)=y by bisection, given Y and guesses x1 and
c x2 that include the solution.
c
c Basic bisection subroutine in which all values that are needed
c by FUNC are in common, and FUNC is a function subroutine
c Modif. so nerr=1 if f(x1) and f(x2) both less than Y
c	    nerr=2 if f(x1) and f(x2) both greater than Y
c and nerr positive if f(x) increases with x
c     nerr negative if f(x) decreases with x
c	ABORT with F2: if aborted returns with nerr=3
c
C*** VERSION OF BISECT WHICH,IF OMIT SET TO TRUE, OMITS THE FIRST
C TWO FUNCTION EVALUATIONS AT X1,X2 SO MUST ENSURE BEFORE CALL
C THAT FUNC(X1) < Y AND THAT X1,X2 REALLY INCLUDE ROOT
c EPSx=accuracy for X (ignored if neg- i.e. Nstep not calc so loops
c	until EPSy is satisfied; max Nstep set to 1000)
c EPSY=accuracy for func(X)-Y (ignored if neg)
c
c NERR=1 IF GUESSES BAD,NERR=0 OTHERWISE
C Prints progress every NDISP steps if ndisp>0
c NDISP=-1  for no print of iterations
c NDISP=-2 omits ALL printing
C
	NDSAV=0
	if(ndisp.ge.-1) print 121
121	format(' Bisection may be aborted with F2')
	NERR=0
	yhi=1.e37
	ylo=-1.e37
	xlo=x1	!need defaults for xlo,xhi in case not set
	xhi=x2
	x11=x1	!save input values for debug
	x22=x2
c
	IF(OMIT) GOTO 1
70	continue
c best guess for theta with fixed param now chosen in elfunc
	fx1=ELFUNC(x1,Xobs,yobs,w,nj,juse,setx,niobs,njset,
     & jfix,ndth,thetsav,xlimsav,youtsav,nfits,nfmax)
	if(iconv.ne.6) then
	   valid=.true.
	else
	   nerr=4	!simplex did not converge
	   print 12
	   print 41,1,x1,fx1
41	   format(i5,3x,g13.6,2x,g13.6,'  Simplex not converged')
	   call SETNM(nevmax,ioptsav)		!below
	   if(ioptsav.eq.3) then
		valid=.false.
c		Yout=fx1
c		Xout=x1
		goto 99
	   else if(ioptsav.eq.2) then
		goto 70
	   else if(ioptsav.eq.1) then
		dx=0.1*(abs(x2-x1))
c		  start with x nearest ML value=thetsav(npar,1)
		th=thetsav(npar,1)
		x0=x1
		if(abs(x1-th).gt.abs(x2-th)) then   !x2 nearer ML value
		   x0=x2
		   dx=-dx	!work down from x2
		endif
		print 12
		do i=1,10
		   x=x0+float(i-1)*dx
		   Yout=ELFUNC(x,Xobs,yobs,w,nj,juse,setx,niobs,njset,
     & 		jfix,ndth,thetsav,xlimsav,youtsav,nfits,nfmax)
		   if(iconv.eq.6) then
			valid=.false.
			print 41,i,x,Yout
c41			format(i5,3x,g13.6,2x,g13.6,'  Simplex not converged')
		   else
			valid=.true.
			Xout=x
			print 4,i,x,Yout
			a=abs(yout)
			if(a.lt.yerr) then
			   xout=x
			   print 80
80			   format(' Limit found')
			   goto 99	!result found!
			endif
			if(yout.ge.0.and.yout.lt.yhi) then
			   yhi=yout
			   xhi=x
			else if(yout.lt.0.and.yout.gt.ylo) then
			   ylo=yout
			   xlo=x
			endif
		   endif
		enddo
		goto 99
	   endif
	endif
	if(ndisp.ge.0) print 25,x1,fx1
c=	fx2=func(x2)
71	continue
	fx2=ELFUNC(x2,Xobs,yobs,w,nj,juse,setx,niobs,njset,
     &  jfix,ndth,thetsav,xlimsav,youtsav,nfits,nfmax)
	if(iconv.ne.6) then
	   valid=.true.
	else
	   print 12
	   print 41,2,x2,fx2
c41	   format(i5,3x,g13.6,2x,g13.6,'  Simplex not converged')
	   nerr=4	!simplex did not converge
	   call SETNM(nevmax,ioptsav)		!below
	   if(ioptsav.eq.3) then
		valid=.false.
c		Yout=fx2
c		Xout=x2
		goto 99
	   else if(ioptsav.eq.2) then
		goto 71
	   else if(ioptsav.eq.1) then
		dx=0.1*(abs(x2-x1))
c		  start with x nearest ML value=thetsav(npar,1)
		th=thetsav(npar,1)
		x0=x1
		if(abs(x1-th).gt.abs(x2-th)) then   !x2 nearer ML value
		   x0=x2
		   dx=-dx	!work down from x2
		endif
		print 12
		do i=1,10
c		  start with x nearest ML value!
		   x=x0+float(i-1)*dx
		   Yout=ELFUNC(x,Xobs,yobs,w,nj,juse,setx,niobs,njset,
     & 		 jfix,ndth,thetsav,xlimsav,youtsav,nfits,nfmax)
		   if(iconv.eq.6) then
			valid=.false.
			print 41,i,x,Yout
c41			format(i5,3x,g13.6,2x,g13.6,'  Simplex not converged')
		   else
			print 4,i,x,Yout
			valid=.true.
			Xout=x
			a=abs(yout)
			if(a.lt.yerr) then
			   xout=x
			   print 80
c80			   format(' limit found')
			   goto 99	!result found!
			endif
			if(yout.ge.0.and.yout.lt.yhi) then
			   yhi=yout
			   xhi=x
			else if(yout.lt.0.and.yout.gt.ylo) then
			   ylo=yout
			   xlo=x
			endif
		   endif
		enddo
		goto 99
	   endif
	endif
	if(ndisp.ge.0) print 25,x2,fx2
25	format(' f(',g13.6,')= ',g13.6)
	if((fx1-Y)*(fx2-Y).gt.0.0) then
	   if(ndisp.ge.0) then
	      print 25,x1,fx1
	      print 25,x2,fx2
	   endif
	   goto 98			!SAME SIGN=ERROR; nerr=3
	endif
c	IF(FX1.LT.Y) GOTO 1
	if(fx1.ge.Y) then
	   P=X1
	   X1=X2
	   X2=P
	endif
1	CONTINUE
cd	print 23,x1,x2,epsx
cd23	format(' x1,x2,epsx= ',3g13.6)
	IF(NDISP.GE.0) print 12
c12	FORMAT(' Step #     X(out)      Y(out)',/)
12	FORMAT(' Step #     param val   Output (should be 0)',/)
c
	Nstep=1000	!max if epsy used to control convergence
	a=abs(x1-x2)
	if(epsx.gt.1.e-10*a) then
	   nstep=IFIX(alog(a/epsx)/alog(2.0)+0.5)
	endif
c
	ns=0
22	continue		!return here for next iteration
	if(KBHIT()) then
		ch=GETCH(ktype)
		if(ktype.eq.0.and.ichar(ch).eq.60) goto 999	!abort
	endif
	ns=ns+1
	XOUT=0.5*(X1+X2)
72	Yout=ELFUNC(Xout,Xobs,yobs,w,nj,juse,setx,niobs,njset,
     &     jfix,ndth,thetsav,xlimsav,youtsav,nfits,nfmax)
	if(iconv.ne.6) then
	   valid=.true.
	else
	   nerr=4	!simplex did not converge
	   print 12
	   print 41,ns,xout,yout
c41	   format(i5,3x,g13.6,2x,g13.6,'  Simplex not converged')
	   call SETNM(nevmax,ioptsav)		!below
	   if(ioptsav.eq.3) then
		valid=.false.
		goto 99
	   else if(ioptsav.eq.2) then
		goto 72
	   else if(ioptsav.eq.1) then
		dx=0.1*(abs(x22-x11))
c		  start with x nearest ML value=thetsav(npar,1)
		th=thetsav(npar,1)
		x0=x11
		if(abs(x11-th).gt.abs(x22-th)) then   !x2 nearer ML value
		   x0=x22
		   dx=-dx	!work down from x2
		endif
		print 12
		do i=1,10
		   x=x0+float(i-1)*dx
		   Yout=ELFUNC(x,Xobs,yobs,w,nj,juse,setx,niobs,njset,
     &	     jfix,ndth,thetsav,xlimsav,youtsav,nfits,nfmax)
		   if(iconv.eq.6) then
			valid=.false.
			print 41,i,x,Yout
c41			format(i5,3x,g13.6,2x,g13.6,'  Simplex not converged')
		   else
			print 4,i,x,Yout
			valid=.true.
			Xout=x
			a=abs(yout)
			if(a.lt.yerr) then
			   xout=x
			   print 80
c80			   format(' limit found')
			   goto 99	!result found!
			endif
			if(yout.ge.0.and.yout.lt.yhi) then
			   yhi=yout
			   xhi=x
			else if(yout.lt.0.and.yout.gt.ylo) then
			   ylo=yout
			   xlo=x
			endif
		   endif
		enddo
		goto 99
	   endif
	endif
	IF(NDISP.LT.0.OR.((NS-NDSAV).LT.NDISP)) GOTO 10
	print 4,NS,XOUT,YOUT
4	FORMAT(I5,3X,G13.6,2X,G13.6)
	NDSAV=NS
10	CONTINUE
c check for best values so far
	a=abs(yout)
	if(a.lt.yerr) then
	   xout=xout
	   print 80
c80	   format(' limit found')
	   goto 99	!result found!
	endif
	if(yout.ge.0.and.yout.lt.yhi) then
	   yhi=yout
	   xhi=xout
	else if(yout.lt.0.and.yout.gt.ylo) then
	   ylo=yout
	   xlo=xout
	endif
c
c Now the bisection!
	erry=yout-y
	if(Yout.lt.y) then
	   x1=Xout
	else if(Yout.gt.y) then
	   x2=Xout
	else if(Yout.eq.y) then
	   goto 99
	endif
c	IF(YOUT-Y) 5,6,7
c5	X1=XOUT
c	GOTO 2
c6	GOTO 99
c7	X2=XOUT
c2	CONTINUE
	if(epsy.gt.0.) goto 21		!use epsy to control convergence
	if(ns.lt.nstep) goto 22		! otherwise use nstep
	goto 99
21	if(abs(erry).gt.epsy) goto 22
	goto 99
C
c Set NERR
c Modif. so nerr=1 if f(x1) and f(x2) both less than Y
c	    nerr=2 if f(x1) and f(x2) both greater than Y
c and nerr positive if f(x) increases with x
c     nerr negative if f(x) decreases with x
98	if(fx1.lt.y.and.fx2.lt.y) nerr=1
  	if(fx1.gt.y.and.fx2.gt.y) nerr=2
	if(abs(x2-x1).lt.1.e-30) then
	   if(ndisp.ge.-1) print 26
26	   format(' Initial guesses essentially equal!')
	else
	   if((fx2-fx1)/(x2-x1).lt.0.0) nerr=-nerr
	endif
c	if(ndisp.ge.-1) print 24	!now in likint1
c24	format(' Guesses do not include result')
99	RETURN
c after aborting
999	nerr=3
	RETURN
	END

	subroutine SETNM(nevmax,ioptsav)
c
	if(ioptsav.eq.1.or.ioptsav.eq.2) iopt=ioptsav
	print 1,nevmax,iopt
1	format(' Simplex did not converge after ',i6,' iterations',/,
     & ' (1) Test a range of values for the limit',/,
     & ' (2) Give higher value for number of iterations',/,
     & ' (3) Try new guess for this limit (or abandon)',/,
     & ' Option number [',i2,'] = ')
	call INPUTi(iopt)
	if(iopt.eq.2) then
	   print 2,nevmax
2	   format(' New maximum number of iterations [',i6,'] = ')
	   call INPUTi(nevmax)
	endif
	ioptsav=iopt
	RETURN
	end
