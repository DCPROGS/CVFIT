	subroutine HESMAT2(theta,sdev,nfix,jfix,SMIN,SSD,kmax,
     & kfit,ndth,Xobs,yobs,w,nj,juse,setx,niobs,njset,cormat)
c HESNAT2 is version of HESMAT1 with option for output in graphical
c mode (controlled by alpha, in common).  Assumes dialog box one set up
c already.  Outputs correlation matrix for display in cvout1.
c
c Version of HESMAT with extra parameter, ndth,=declared (or allocated)
c dimension of theta,sdev and jfix in calling prog. And internal arrays
c allocated to kmax.
c
c ### ELIK=NAME OF FUNCTION THAT CALCULATES MINUS LOG-LIKELIHOOD
c ### ARGS are SSD,Smin (not ELIK,ELMAX)
c 06/06/89 08:00pm Lahey version
c Version of HMAT1 that
c (a) takes up to 20 parameters
c (b) scales parameters internally for better calc of derivatives (as
c	in HMAT2.FOR on PDP-11)
c Used with VMAT2 which returns scaled FINFO, COVAR + scaling factors
c in SCAL(i),i-1,...,kfit. (Note that a parameter like g(0)=1.e-30 with
c sd=1.e-31 would have variance=1.e-62 which will underflow, so must keep
c in scaled form until last moment)
c Also returns SDEV= standard dev of params (corrected for scaling) rather
c than variance as in other versions
c
	real*4 sdev(ndth),THETA(ndth)		!param
	integer jfix(ndth)			!param
	real*4 cormat(kmax,kmax)
	ALLOCATABLE::scal,finfo,covar,unit,np
	real*4 scal(:)
	real*8 finfo(:,:),covar(:,:),unit(:,:)
	integer np(:)
	real*8 det,u,den
	LOGICAL DEBUG,pon
	logical discprt,slock,caplock,prt
c Params to be passed to SSDCV and Simplex
	real*4 xobs(niobs,njset),yobs(niobs,njset),w(niobs,njset)
	integer juse(njset),nj(njset)
	real setx(njset)
c For graphical output
c	character defolt*30,cans*30		!to hold text & result of QDIALOG
	logical help
	common/hlp/help				!for QDIALOG
	COMMON/determ/det		!for MATINV
	common/dp/discprt
	logical alpha
	common/graf/alpha		!for simplexv
	common/popvalt/lxlo,lylo,lyhi,ictt,ibkt,icft,icht	!values for poptable calls
	external SSD
c
	pon()=slock()
	debug()=caplock()
c
101   format(a1)
c	deb=caplock()
	prt=slock()
c
	ALLOCATE(scal(kmax),finfo(kmax,kmax),covar(kmax,kmax),
     & unit(kmax,kmax),np(kmax))
c
	fract=-1.0
c	if(alpha) then
c	   print 31
c31	   format('&Use default step size for error calculation [Y] ? ')
c         read 101,ans
c	else
c	   ans='Y'
c	   call DEFOLTa(ans,defolt)
c	   call QDIALOG(1,
c     &	'Use default step size for error calculation',
c     &	  defolt,11,cans)
c	   call GETINPa(cans,ans)
c	endif
c      if(UC(ans).ne.'N') goto 10
cC
c	fract=0.05
c	if(alpha) then
c	   print 11
c11	   FORMAT('&Fraction for calc of Hessian = ')
c	   call INPUTr(fract)
c	else
c	   call DEFOLTr(fract,defolt)
c	   call QDIALOG(1,
c     &	'Fraction for calc of Hessian',
c     &	  defolt,11,cans)
c	   call GETINPr(cans,fract)
c	endif
c10	CONTINUE
c
	call VMAT3(theta,Kmax,fract,finfo,covar,unit,SSD,
     & smin,nfix,jfix,SCAL,
     & ndth,kmax,Xobs,yobs,w,nj,juse,setx,niobs,njset)
c
	if(alpha) print 14,(scal(i),i=1,kfit)
	if(prt) write(7,14) (scal(i),i=1,kfit)
	if(discprt) write(8,14) (scal(i),i=1,kfit)
14	FORMAT(' Scale factors= ',10g9.3)
	if(alpha) print 141
	if(prt) write(7,141)
	if(discprt) write(8,141)
141	format(//,
     & ' Observed information matrix = ')
	do i=1,kfit
	   if(i.le.6) then
		if(alpha) print 45,(FINFO(I,J)/(scal(i)*scal(j)),J=1,I)
		if(prt) write(7,45)(FINFO(I,J)/(scal(i)*scal(j)),J=1,I)
	      if(discprt) write(8,45)(FINFO(I,J)/(scal(i)*scal(j)),J=1,I)
	   else
		if(prt) write(7,45)(FINFO(I,J)/(scal(i)*scal(j)),J=1,I)
      	if(discprt) write(8,45)(FINFO(I,J)/(scal(i)*scal(j)),J=1,I)
		if(alpha) print 45,(FINFO(I,J)/(scal(i)*scal(j)),J=1,6)
		if(prt) write(7,451)(FINFO(I,J)/(scal(i)*scal(j)),J=7,I)
	      if(discprt) write(8,451)(FINFO(I,J)/(scal(i)*scal(j)),J=7,I)
		if(alpha) print 451,(FINFO(I,J)/(scal(i)*scal(j)),J=7,I)
	   endif
	enddo

	if(alpha) print 3711
	if(prt) write(7,3711)
	if(discprt) write(8,3711)
3711	FORMAT(/,' Covariance matrix = ')
	do i=1,kfit
	   if(i.le.6) then
		if(prt) write(7,45)(COVAR(I,J)*scal(i)*scal(j),J=1,I)
	      if(discprt) write(8,45)(COVAR(I,J)*scal(i)*scal(j),J=1,I)
		if(alpha) print 45,(COVAR(I,J)*scal(i)*scal(j),J=1,I)
	   else
		if(prt) write(7,45)(COVAR(I,J)*scal(i)*scal(j),J=1,6)
      	if(discprt) write(8,45)(COVAR(I,J)*scal(i)*scal(j),J=1,6)
		if(alpha) print 45,(COVAR(I,J)*scal(i)*scal(j),J=1,6)
		if(prt) write(7,451)(COVAR(I,J)*scal(i)*scal(j),J=7,I)
	      if(discprt) write(8,451)(COVAR(I,J)*scal(i)*scal(j),J=7,I)
		if(alpha) print 451,(COVAR(I,J)*scal(i)*scal(j),J=7,I)
	   endif
	enddo
c13	CONTINUE
c
C RE-ASSIGN VAR(I) (BEFORE COVAR REDEFINED AS CORRELATION MATRIX)
C SO THAT ALL KMAX VALUES INCLUDED, WITH VAR SET TO -1.0
C FOR FIXED PARAMETERS. ALSO DEFINE LOCAL ARRAY NP(I),I=1,2,..,KFIT
C WITH NP(I)=PARAM NUMBER CORRESPONDING TO I'TH ROW AND COL OF COVAR.
C E.G. IF KMAX=4 AND PARAM 2 IS FIXED NP=(1,3,4)
	i1=0
	do i=1,kmax
	   if(jfix(i).eq.0) then
		i1=i1+1
		c=covar(i1,i1)
		if(c.ge.0.0) then
		   sdev(i)=sqrt(c)*scal(i1)	!unscaled SD
		else
		   sdev(i)=-1.0
		endif
	   else
		sdev(i)=-1.0
	   endif
	enddo
C
C CALC CORRELATION MATRIX IN FINFO AND type it (WITHOUT CHANGING COVAR)
	CONTINUE
	if(alpha) print 421
	if(prt) write(7,421)
	if(discprt) write(8,421)
421	format(/,' Correlation matrix=')
	do i=1,kfit
	   do j=1,kfit
		den=covar(i,i)*covar(j,j)		!still scaled
		if(den.gt.1.e-30) then
		   finfo(i,j)=covar(i,j)/sqrt(den)
		else
		   finfo(i,j)=9999.	!den negative
		endif
	   enddo
	enddo
C
c Print correlation matrix
	DO 43 I=1,KFIT
	IF(I.GT.6) GOTO 56
	if(prt) write(7,45) (FINFO(i,j),j=1,i)	!=correlation matrix
      if(discprt) write(8,45) (FINFO(i,j),j=1,i)
	if(alpha) print 45, (FINFO(i,j),j=1,i)
45	FORMAT(2X,6G13.6)
	GOTO 43
56	if(prt) write(7,45) (FINFO(I,J),J=1,6)
      if(discprt) write(8,45) (FINFO(I,J),J=1,6)
	if(alpha) print 45, (FINFO(I,J),J=1,6)
	if(prt) write(7,451)(FINFO(I,J),J=7,I)
      if(discprt) write(8,451)(FINFO(I,J),J=7,I)
	if(alpha) print 451,(FINFO(I,J),J=7,I)
451	FORMAT(15X,5G13.6,/)	!INDENT REMAINING VALUES IF ANY
43	CONTINUE
	if(prt) write(7,46)DET	!det is in common, from matinv
      if(discprt) write(8,46)DET
	if(det.lt.1.d-10) then
	 if(alpha) then
	   print 46,DET
46	   FORMAT(/' Determinant of INFO matrix=',G13.6//)
	 else
	   call WDIALOG(1,'Determinant of INFO matrix near zero',12)
	 endif
	endif
	DO 471 I=1,KFIT		!IF UNIT O.K. TO 1 IN 1E9 DO NOT PRINT
	DO 471 J=1,KFIT
	U=0.0D0
	IF(I.EQ.J) U=1.0D0
471	IF(DABS(UNIT(I,J)-U).GT.1.0D-9) GOTO 472
	GOTO 50
472	DO 47 I=1,KFIT
	if(prt) write(7,45)(UNIT(I,J),J=1,KFIT)
      if(discprt) write(8,45)(UNIT(I,J),J=1,KFIT)
	if(alpha) print 45,(UNIT(I,J),J=1,KFIT)
47	CONTINUE
50	continue
C
c If in graphics, output correlation matrix for display in cvout1.
	do i=1,kfit
	   do j=1,kfit
		cormat(i,j)=sngl(finfo(i,j))
	   enddo
	enddo
C
	DEALLOCATE(scal,finfo,covar,unit,np)
	RETURN
	END



