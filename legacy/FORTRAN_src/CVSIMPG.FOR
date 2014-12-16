	subroutine CVSIMPG(kfit,jfix,nfix,default,titlep,fitted,dcurve,
     & nodata,theta,thetsav,nev,sres,smin,noerr,iflag,isepfit,idest,
     & fixmax,fix0,Xobs,yobs,w,nj,juse,setx,niobs,njset,ndth,
     & badset)
c Version of CVSIMP that calls SIMPLEXg (draws SSD contours)
c
c To fit set #jset in CVFIT (NB jset now in common), or fits all
c sets simultaneously if ifitmode=3,4,5 (see SSCDV)
c
c IFLAG=2 when CVSIMP is doing sep fits for indiv curves (ifitmode=2) as
c the last section of ifitmode=4,5 calcs, after simultaneous fits
c
c Modif 03/27/95 11:34am for allocatable xobs() etc, so all data must
c be passed as params to SIMPLEX and SSDSV
c	function SSDCV(K,THETA,Xobs,yobs,w,nj,juse,setx,niobs,njset)
c
c Modif 09/18/94 06:15pm for VGA mode input
	integer*2 videotyp
	character*20 TITLEP(ndth)		!names of params for fit
	real*4 theta(ndth),thetsav(100,11)
	integer JFIX(ndth)
	logical discprt,noerr,fixmax,fix0,alpha,skipymax,nodata
	LOGICAL constr,default,pon,fitted,dcurve,slock,negpar,fixratio
	logical badset,skipq
	character*1 ans,UC,ans1
	character cnum*11,cnum1*11
	character defolt*30,cans*30		!to hold text & result of QDIALOG
	character*25 text
	character*55 text1
	logical help
c Params to be passed to SSDCV and Simplex
	real*4 xobs(niobs,njset),yobs(niobs,njset),w(niobs,njset)
	integer juse(njset),nj(njset)
	real setx(njset)
	allocatable:: aKint,Ymint
	real aKint(:),Ymint(:)
c
	allocatable:: xsort	!for sorting
	real*4 xsort(:)
c For contours
	real svals(10)
	logical draw
c
	common/hlp/help				!for QDIALOG
c====nj(),jmiss(),juse() no longer in common/block1/
	COMMON/BLOCK1/constr,nset,nfit,nsfit,
     & Xv,alfa,kmax,ncomp,nmod,fline,nomit,jomit(20),jset,ifitmode
	logical logyfit,norm
	COMMON/BLOCK3/logyfit,norm,xnorm		!CVIN2,YCALCV,CVDISP,GETEQN,CVSIMP
	common/fix/fixratio,ir1,ir2,rval	!to fix ratio of 2 parameters
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
c
	errfac=1.e-4
	delmin=-1.		!do not use delmin for convergence
	stpfac=0.1
	nfix=0		!no param fixed
	do i=1,100
	   jfix(i)=0
	enddo
	skipq=.false.
c
	print 108
108	format(/)
	negpar=.false.
c Deal with case where separate fits of Langmuir/Hill are being done
c following simultaneous fit with ifitmode=4,5.  In these cases
c want values that were fixed for simult fit to be still fixed for
c separate fits
c First the case where separate fits have maximum FIXED at the value
c previously found for simultaneous fit
	if(nmod.eq.1.or.nmod.eq.-1) then
	 if(isepfit.eq.1) then
	   if(fixmax) then
	      jfix(2)=1		!Ymax
	      nfix=nfix+1
	   else if(nmod.eq.-1.and.fix0) then	!fix Y(0) also
		jfix(3)=1
	      nfix=nfix+1
	   endif
	 endif
	 if(isepfit.eq.2) then
	   if(alpha) then
	      print 41,thetsav(2,11)
41	      format(
     &    ' Maximum fixed at',g13.6,' (from simult fit): O.K. [Y] ')
	      read 101,ans
	   else
		ans='Y'
		call DEFOLTa(ans,defolt)
	      call REALTOCH(thetsav(2,11),cnum,11)
		n=nblank(cnum)
		call QDIALOG(1,
     &	'Maximum fixed at '//cnum(1:n)//' (from simult fit): O.K.',
     &	  defolt,11,cans)
		call GETINPa(cans,ans)
	   endif
	   if(UC(ans).ne.'N') then
	      theta(2)=thetsav(2,11)	!ymax
	      jfix(2)=1
	      nfix=nfix+1
	   endif
	   if(nmod.eq.-1.and.fix0) then	!fix Y(0) also
		jfix(3)=1
	      nfix=nfix+1
	   endif
	 endif
	endif
	if(iflag.eq.2.and.(nmod.le.4.or.nmod.eq.99.or.
     &	nmod.eq.24.or.nmod.eq.25)) then
c 	  Init guesses already set in main prog -put them into theta() now
	   j=jset		!index for thetsav()
	   if(nmod.eq.23.or.nmod.eq.24) then
		j=1		!index for thetsav() (see NEWEQCV)
		if(fixmax) then
		   jfix(1)=1            !max was fixed for simult fit, so fix now
	         nfix=nfix+1
		endif
		if(fix0.and.nmod.eq.24) then
		   jfix(2)=1            !Y(0) was fixed for simult fit, so fix now
	         nfix=nfix+1
		endif
	   endif
	   do i=1,kmax
		theta(i)=thetsav(i,j)	!copy init guesses into theta
	   enddo
	   fitted=.true.
	   dcurve=.false.
	   ndisp=20
	   skipq=.true.
	   goto 524				!and skip questions
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
	skipymax=.false.
	if(alpha) then
	   print 1081,text
1081	   format('&',a25,/)
	else
	   call WDIALOG(1,text,11)
	endif
	do 131 i=1,kmax
	  i1=i
	  if(NMOD.eq.7.and.i.gt.3) i1=3
	  if(i.eq.2.and.skipymax) then	!skipymax set below in loop
	     theta(i)=1.0
	     nfix=nfix+1
	     jfix(i)=1
	     if(alpha) then
		print 42
42	 	format(' Ymax will cancel so not needed')
	     else
		call WDIALOG(1,' Ymax will cancel so not needed',11)
	     endif
	     goto 131	!next param	  endif
	  endif
1321	  continue
	  if(alpha) then
	     print 132,i,TITLEP(i1)
132	     format('&Parameter ',i2,': ',a20,'= ')
	     call INPUTr(theta(i))
	  else
	     call INTCONV(i,cnum)
	     call QDIALOG(1,'(Param('//CHARNB(cnum)//
     &      '): '//CHARNB(titlep(i1)),' ',11,cans)
	     call GETINPr(cans,theta(i))
	  endif
	  if(theta(i).lt.0.) negpar=.true.
	  if(nodata.or.dcurve) goto 131
	  if(abs(theta(i)).lt.1.e-20) then
		call BELL(2)
		if(alpha) then
		   print 133,i
133		   format(
     &'&Initial guess of zero not allowed unless you intend to fix',
     &	' this parameter',/,
     &	' Fix parameter number ',i2,' at zero [Y] ? ')
		   read 101,ans
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
		   nfix=nfix+1
		   jfix(i)=1
		   if(nmod.eq.26.or.nmod.eq.27.and.ifitmode.eq.1
     &		.and.i.eq.1.and.ncomp.eq.1) then
			if(norm) skipymax=.true.
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
	   if(ifitmode.eq.4.or.ifitmode.eq.5) j=11	!simultaneous fit
	   do i=1,kmax
		thetsav(i,j)=theta(i)
	   enddo
	   goto 201
	endif
	dcurve=.false.
c	if(alpha) then
c	   print 200
c200	   FORMAT('&Use these values only to plot curve [N] ? ')
c	   read 101,ans
c	else
c	   ans='N'
c	   call DEFOLTa(ans,defolt)
c	   call QDIALOG(1,'Use these values only to plot curve',
c     &	defolt,11,cans)
c	   call GETINPa(cans,ans)
c	endif
c	if(UC(ans).eq.'Y') then
c	   fitted=.false.
c	   dcurve=.true.
cc Save parameters
c	   j=jset
c	   if(ifitmode.eq.4.or.ifitmode.eq.5) j=11	!simultaneous fit
c	   do i=1,kmax
c		thetsav(i,j)=theta(i)
c	   enddo
c	   goto 201
c	endif
c
c	if(default) goto 526
c	print 525
c525	format('&Relative error [1.e-4]= ')
c	call INPUTr(errfac)
c	if(x.gt.1.e-10) errfac=x
cc	if(errfac.lt.0.) delmin=-errfac	!use delmin criterion
c526	continue
c
c Option to omit specified points if only one set being fitted
c ==Omission now  best done by setting weight=0, so point omitted from plot too.
c	if(default.or.nset.ne.1) goto 210
c	print 211
c211	format('&Omit any points [N]? ')
c	read 101,ans
c	if(UC(ans).ne.'Y') goto 210
c	print 217
c217	format('&Number to be omitted= ')
c	read 104,nomit
c	do 218 i=1,nomit
c	print 219,i
c219	format('&Point #',i2,'= ')
c	read 104,jomit(i)
c	if(pon()) write(7,2122)jomit(i)
c      if(discprt) write(8,2122)jomit(i)
c2122	format(' Point #' ,i4,' omitted')
c218	continue
c210	continue
c	if(default) goto 1162
	if(nfix.gt.0) then	!may have fixed some at zero above
	   if(alpha) then
		print 1161,nfix
1161		format('& ',i2,
     & ' parameters already fixed: number of others to be fixed [0] = ')
		i1=0
		call INPUTi(i1)
	   else
	      call INTCONV(nfix,cnum)
		i1=0
	      call DEFOLTi(i1,defolt)
	      call QDIALOG(1,CHARNB(cnum)//
     &      ' parameters already fixed: fix n more; n',defolt,11,cans)
	      call GETINPi(cans,i1)
	   endif
	   if(i1.gt.0) nfix=nfix+i1
	   n=i1		!number to get now
	else
	   i1=0
	   if(alpha) then
		print 116
116		FORMAT('&Number of parameters to be fixed [0]= ')
		call INPUTi(i1)
	   else
	      call DEFOLTi(i1,defolt)
	      call QDIALOG(1,
     &      'Number of parameters to be fixed',defolt,11,cans)
	      call GETINPi(cans,i1)
	   endif
	   if(i1.gt.0) nfix=nfix+i1
	   n=nfix		!number to get now
	endif
	kfit=KMAX-NFIX
	if(i1.eq.0) goto 117    !no more to fix
	do i=1,n
	   if(alpha) then
		print 118,i
118		format('&(',i2,'):   parameter #= ')
		call INPUTi(itemp)
	   else
		call INTCONV(i,cnum1)
		call QDIALOG(1,' '//CHARNB(cnum1)//
     & 	    ': Fix parameter #',' ',11,cans)
		call GETINPi(cans,itemp)
	    endif
	   JFIX(itemp)=1		!ELEMENT OF JFIX=1 IF PARAM FIXED
	enddo
117	continue
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
		call INPUT2r(ir1,ir2)
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
	if(ans.eq.'N') goto 124
c
524	continue
	call VIDEOMOD(3)		!utility lib- this makes graph go for good!
	alpha=.true.
	call LOCATE(0,0)
c Check whether fit is possible -if number of different x values
c is less than number of free parameters then fit cannot be done
c (at least when sep fits are being done)
c First check if there are ANY duplicated x values (fast), and if there
c are then sort x values to count number of different values
	badset=.false.
	if(ifitmode.le.2) then
	   n=nj(jset)
	   nxval=1
	   xv1=Xobs(1,jset)
	   do i=2,n
		xv=Xobs(i,jset)
		if(xv.ne.xv1) then	!at least one duplicate x value
		   nxval=nxval+1
		else
		   goto 43			!sort values
		endif
	   enddo
	   goto 44
43	   continue			!sort values
	   ALLOCATE(xsort(n))
	   call SORTr(xsort,n,n,.true.)
	   nxval=1
	   do i=2,n
		if(xsort(i).ne.xsort(i-1)) nxval=nxval+1
	   enddo
c
44	   continue			!nxval found
	   kfit=kmax-nfix    !in case jump to 524
	   if(nxval.ge.kfit) then
	      print 45,kfit,nxval
	      if(pon()) write(7,45) kfit,nxval
            if(discprt) write(8,45) kfit,nxval
45	      format(
     &	' Number of free parameters = ',i3,/,
     &	' Number of different x values = ',i4)
	   else
		call BELL(1)
	      print 46,jset,kfit,nxval
	      if(pon()) write(7,46) jset,kfit,nxval
            if(discprt) write(8,46) jset,kfit,nxval
46	      format(' Not enough points for fit to set ',i3,/,
     &	' Number of free parameters = ',i3,/,
     &	' Number of different x values = ',i4)
		badset=.true.
		RETURN
	   endif
	endif
c
	if(skipq) goto 5241
	irestrt=3
	print 312,irestrt
312	FORMAT(' Limit number of restarts to [',i2,']: n = ')
	call INPUTi(irestrt)
c
	print 525
525	format('&Relative error [1.e-4]= ')
	call INPUTr(errfac)
c	print 523,errfac
c523	FORMAT(' ERROR=guess*',g13.6,/,
	ndisp=20
	print 523
523	format(
     & ' Print every Nth estimate [20] :  N= ')
	call INPUTi(ndisp)
c
5241	continue
	kfit=kmax-nfix    !in case jump to 524
	nevm=-200000		!so iterations not printed to disc
	confac=0.5		!parameter for simplex3
	smin=SSDCV(kmax,theta,Xobs,yobs,w,nj,juse,setx,niobs,njset)
	if(kmax.ge.7) print 5231,1,smin,(theta(I),I=1,KMAX)
	if(kmax.lt.7) print 5232,1,smin,(theta(I),I=1,KMAX)
5231	format(i8,g10.3,/,10g10.3,/,10g10.3,/)
5232	format(i6,7g11.4,/)
c
c Set last line of parameters for drawing contours
	draw=.false.
	if(nmod.eq.26) then
	   jY=2	!index of Ymax in theta
	   jK=3	!index of K in theta
	   draw=.true.
	   intmax=1000
	   ALLOCATE(akint(intmax),ymint(intmax))
	   nint=0
	   nint=nint+1
	   aKint(nint)=theta(jK)	!initial K value
	   Ymint(nint)=theta(jY)	!initial Ymax value
	   ssdmax=0.
c	   ndisp=-1000
	endif
c
	call SIMPLEXG(kmax,THETA,stpfac,errfac,nev,nevm,
     & smin,SSDCV,Ndisp,jfix,-1.,confac,irestrt,iconv,
     & Xobs,yobs,w,nj,juse,setx,niobs,njset,ndth,
     & draw,akint,ymint,nint,intmax,ssdmax,jK,jY)
c===     & draw,aK1,aK2,ncont,svals,jset,jK,jY)
c
c Draw SSD contours
	if(draw) then
	  print 10
10	  format(' Draw SSD contours [Y] ? ')
	  read 101,ans
	  if(UC(ans).ne.'N') then
	   nint=nint+1
	   aKint(nint)=theta(jK)	!final K value
	   Ymint(nint)=theta(jY)	!final Ymax value
	   if(xobs(1,jset).eq.2.5.and.yobs(1,jset).eq.5.576) then
c Values for LOB example
		aK1=0.
		aK2=40.
		ncont=6
		svals(1)=10.
		svals(2)=20.
		svals(3)=50.
		svals(4)=200.
		svals(5)=500.
		svals(6)=2000.
	   else
		aK=theta(jK)
     		aK1=0.
		aK2=10*aK
		ncont=6
		svals(1)=ifix(smin)+1.
		svals(ncont)=ifix(ssdmax)+10.
		ds=(svals(ncont)-svals(1))/float(ncont)
		do i=2,ncont
		   svals(i)=svals(1)+float(i)*ds
		enddo
	   endif
22	   continue
	   call DRAWCONT(aK1,aK2,ncont,svals,jset,xobs,yobs,w,
     &	niobs,nj,njset,nint,aKint,Ymint)
	   ans='Y'
	   call DEFOLTa(ans,defolt)
	   call QDIALOG(1,'Recalculate contours',defolt,11,cans)
	   call GETINPa(cans,ans)
	   if(ans.eq.'Y') then
		call WDIALOG(1,'Present SSD for contours are:',11)
		j=1
		do i=1,ncont
		   if(pon()) write(7,21) i,svals(i)
		   if(discprt) write(8,21) i,svals(i)
21		   format('SSD for contour #',i3,' = ',g13.6)
	         call REALTOCH(svals(i),cnum,11)
		   n=NBLANK(cnum)
		   if(j+n-1.le.55) then
		      text1(j:j+n-1)=cnum(1:n)
		      j=j+n
		      if(j.lt.55.and.i.lt.ncont) then
			   text1(j:j)=','
		         j=j+1
		      endif
		   else
			call WDIALOG(1,text1,11)
			j=1
			do k=1,55
			   text1(k:k)=char(32)
			enddo
		   endif
		enddo
		call WDIALOG(1,text1,11)
	      call DEFOLTi(ncont,defolt)
	      call QDIALOG(1,'Number of contours',defolt,11,cans)
		call GETINPi(cans,ncont)
		do i=1,ncont
c	         call REALTOCH(svals(i),cnum,11)
		   call DEFOLTr(svals(i),defolt)
	         call INTCONV(i,cnum)
	         call QDIALOG(1,'Contour '//charnb(cnum)//': SDD',
     &	     defolt,11,cans)
	         call GETINPr(cans,svals(i))
		enddo
		goto 22
	   else
		DEALLOCATE(akint,ymint)
		call VIDEOMOD(3)		!utility lib- this makes graph go for good!
	   endif
	  endif
	endif

C 100=dimension of theta(), jfix().
	sres=-1.
	if(nfit.gt.kfit) then
	   var=smin/float(nfit-kfit)		!error var
	   if(var.gt.1.e-30) sres=sqrt(var)
	endif
	print 5233
5233	FORMAT( ' END OF FITTING')
	IF(KMAX.GE.7) print 5231,NEV,SMIN,(THETA(I),I=1,KMAX)
	IF(KMAX.LT.7) print 5232,NEV,SMIN,(THETA(I),I=1,KMAX)
	ndf=nfit-kfit		!degrees of freedom for error
	print 5234,nfit,kfit,sres,ndf
	if(pon()) write(7,5234) nfit,kfit,sres,ndf
      if(discprt) write(8,5234) nfit,kfit,sres,ndf
5234	format(' Number of points fitted = ',i8,/,
     & ' Number of parameters estimated = ',i8,/,
     & ' Error S.D.= ',g13.6,'  (',i5,' degrees of freedom)')
c Save parameters in case several fits done
	j=jset
	if(ifitmode.eq.4.or.ifitmode.eq.5) j=11	!simultaneous fit
	do 20 i=1,kmax
20	thetsav(i,j)=theta(i)
c
C
C CALC VARIANCE ETC FOR THE Kfit NON-FIXED PARAMETERS
	ans1='y'
	if(iflag.eq.2) ans1='n'
	call DCASK('Calculate approximate errors',ans1,ans)
	noerr=ans.eq.'N'
	if(noerr) goto 201
	idest=0
	RETURN
201	idest=201
	RETURN
	end


