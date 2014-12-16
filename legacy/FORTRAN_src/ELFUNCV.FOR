	function ELFUNCV(X,Xobs,yobs,w,nj,juse,setx,niobs,njset,
     & jfix,ndth,thetsav,xlimsav,youtsav,nfits,nfmax)
c
	real*4 thetsav(kmax,nfmax),xlimsav(nfmax),youtsav(nfmax)
	integer JFIX(ndth)
	real theta1(100)
c	real theta2(100)
c Params to be passed to SSDCV and Simplex
	real*4 xobs(niobs,njset),yobs(niobs,njset),w(niobs,njset)
	integer juse(njset),nj(njset)
	real setx(njset)
	COMMON/LIKBLK/kmax,kfit,elmax,em,delmin,ndisp,npar,iw,sres
	common/simp/iconv,nevmax	!so convergence testable in likint1
	EXTERNAL SSDCV
C
c 12/28/89 04:44pm Lahey version of ELFUNC for CVFIT
c
c 11-73 version Sept 1986.
C Note that if KFIT=1, i.e. only one param to be estimated then
c when this is fixed at a trial value there are none to be
c estimated so SIMPLEX not needed- just call EXPLIK etc
C FUNCTION FOR USE IN CALC OF LIKELIHOOD INTERVALS. INPUT X IS
C VALUE FOR THE FIXED PARAM FOR WHICH LIMITS WANTED-IT IS PARAM NUMBER
C NPAR . JFIX SPECIFIES WHICH PARAMETER HAS BEEN
C FIXED (THAT FOR WHICH LIMITS WANTED) (OTHER PARAM WHICH
C HAVE BEEN FIXED THROUGHOUT STAY FIXED OF COURSE).
C	VALUE OF FUNC IS ELMAX1-ELCRIT WHERE ELMAX1= MAX LIK GIVEN THE
C CONSTRAINED PARAMETERS - ELCRIT= L(THETA HAT)-EM. THUS VALUE OF
C FUNC SHOULD BE ZERO AT SOLUTION FOR EM-UNIT INTERVAL
C	NEED INITIAL GUESSES (IN BP) FOR THE PARAM FOR THIS CONSTRAINED
C FIT. THE FIXED PARAM IS SET TO THE CURRENT GUESS FOR ITS (UPPER OR
C LOWER) LIMIT. FOR THE OTHER PARAM PROB CANNOT DO BETTER THAN
C USE ML VALUES.
C	SPECIFICATION OF ERROR OF MINIMISATION TRICKY SINCE IT IS
C THE MIN THAT WE WANT,NOT THE PARAM VALUES AT THE MIN
C
c
c Choose the best initial guess for the non-fixed parameters from
c those known to converge OK, stored in thetsav()
	dmin=1.e37
	do j=1,nfits
	   d=abs(x - xlimsav(j))
	   if(d.lt.dmin) then
		dmin=d
		jmin=j
	   endif
	enddo
	do i=1,kmax
	   if(i.ne.npar) theta1(i)=thetsav(i,jmin)
	enddo
c
	theta1(NPAR)=X		!SET TO CURRENT GUESS (FIXED)
c
c NB SSDCV gives -2*log likelihood
	if(kfit.eq.1) then
	   elmax1=-0.5*SSDCV(kmax,theta1,Xobs,yobs,w,nj,juse,setx,
     &	niobs,njset)
	   goto 31
	endif
c
c30	continue
	stpfac=0.1
	errfac=1.e-4
	irestrt=3
	confac=0.5		!parameter for simplex3
c	do i=1,kmax
c	   theta2(i)=theta1(i)	!store input guesses
c	enddo
	call SIMPLEXv(kmax,theta1,stpfac,errfac,nev,nevmax,
     & Smin,SSDCV,ndisp,jfix,delmin,confac,irestrt,iconv,
     & Xobs,yobs,w,nj,juse,setx,niobs,njset,ndth)
c NB SSDCV gives -2*log likelihood
	elmax1=-smin/2.0
c
31	continue
c
	if(iw.eq.1.or.iw.eq.5) elmax1=elmax1/(sres**2)	!use Sres as error var
c
	Yout=ELMAX1-(ELMAX-EM)	!elmax1-elcrit=0, ie elmax1=elcrit
c
c Store output to use as initial guesses for next Simplex, IF this one
c converged
c -or better store input? -these are values known to result in convergence!
	if(iconv.ne.6) then
	   nfits=nfits+1
	   if(nfits.le.nfmax) then
		do i=1,kmax
		   thetsav(i,nfits)=theta1(i)		!store output
c		   thetsav(i,nfits)=theta2(i)		!store input
		enddo
		xlimsav(nfits)=theta1(npar)	!the fixed value
		youtsav(nfits)=Yout
	   endif
	endif
c
	ELFUNCV=Yout
	RETURN
	END


