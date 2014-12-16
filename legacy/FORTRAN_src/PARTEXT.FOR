	subroutine PARTEXT(text,ifitype,theta,ncomp,ifont,csize,kmax)
c To make the text string to display fitted parameters on a graph
c Call so string constructed only once, not each time graph redrawn.
c  Csize is char size for the text (though will be needed here only if
c size changed eg for subscript while making TEXT; otherwise size is
c needed only at the time TEXT is actually drawn)
c
c Nodif 01/01/97 05:19pm num1 etc made char*9 rather than *8 (and define n1=9
c to facilate changing it again)
c Modif 04/04/95 02:22pm for VPLOT5 by adding kmax=declared dimension of theta
c Ifont (0-11 scale) is the default font for TEXT.
c	To decode THETA need to know what sort of fit was done. This is defined
c in EKDIST by THIST etc, but if queued from other progs to VHIST,VPLOT
c (eg in AUTPLT) better to define this by an integer IFITYPE
c IFITYPE=0 for no display of param
c IFITYPE=1 for time histos (THIST)
c IFITYPE=2 for op/bst histos (OBHIST)
c IFITYPE=3 for amp  histos (AHIST)
c IFITYPE=31 for amp  histos (AHIST) if SD are proportional to mean
c IFITYPE=32 for amp  histos (AHIST) if SD are constrained to be equal
c   (older version used IFITYPE=-3 to indicate equal SD -this still recognized)
c IFITYPE=4 for exponentials fitted to relaxations (CJFIT,VPLOTR)
c IFITYPE=5 for Lorentzians fitted to spectra (PSFIT,VPLOTR)
c IFITYPE=-5 for Lorentzians + line fitted to spectra (PSFIT,VPLOTR)
c IFITYPE=6 for general exponentials (as in CVFIT) (tau, amp)
c IFITYPE=-6 for general geometrics (as in CVFIT) ('mean', amp)
	REAL theta(kmax)
	allocatable:: area,amean,sd
	real*4 area(:),amean(:),sd(:)
c=	REAL AREA(10),AMEAN(10),SD(10)
c	character text*150
	character text*(*)
	character*9 num1,num2,num3
	character*2 fnt,grk,dup
	character*1 nline,sup,sub
c	logical caplock
c
	do 12 i=1,150
12	text(i:i)=char(32)
c
	k2=kmax/2.
	ALLOCATE(area(k2),amean(k2),sd(k2))
c
c Define string control characters (Hgraph p57)
	nline=char(14)		!NB not char(12) as in manual!
	sub=char(4)			!subscript
	sup=char(21)		!subscript
	grk=char(6)//'A'		!to set greek font
	dup=char(6)//'2'		!to set duplex font (eg so x=multiply sign)
	grk=char(6)//'A'		!to set greek font
	fnt=char(6)//char(ifont+48)		!to set default font for string
	call SETFONT1(ifont)	!font for string
	call SETCSIZE(csize,isize)	!this defines ISIZE
c
	isdcon=0		!other than amp histos
	if(ifitype.eq.31) then
	   isdcon=2
	   ifitype=3
	else if(ifitype.eq.32.or.ifitype.eq.-3) then	!include -3 for old queues
	   isdcon=3
	   ifitype=3
	endif
c
	if(ifitype.eq.4) goto 60	!sep section for jumps
	if(iabs(ifitype).eq.6) goto 60	!general exponentials/geometrics
	if(iabs(ifitype).eq.5) goto 70	!sep section for Lorentzians
c
	if(ifitype.ne.1) goto 10
	J=0
	SA=0.
	DO 11 I=1,NCOMP
	J=J+1
	AMEAN(I)=THETA(J)
	IF(I.EQ.NCOMP) GOTO 11
	J=J+1
	AREA(I)=THETA(J)
	SA=SA+AREA(I)
11	CONTINUE
	AREA(NCOMP)=1.0-SA
	goto 40
C
10	IF(ifitype.ne.3) GOTO 30
C NOW AMPLITUDE FIT
	J=0
	SA=0.
	do i=1,ncomp
	   j=j+1
	   amean(i)=theta(j)		! mean(i)
	   if(isdcon.eq.3) then
		if(i.eq.1) j=j+1
		sd(i)=theta(2)		!if SD fixed, otherwise redefined
	   else if(isdcon.eq.2) then
		if(i.eq.1) j=j+1
		sd(i)=theta(2)*amean(i)	!theta(2)=sdfac in this case
	   else				!isdcon=1
		j=j+1
		sd(i)=theta(j)
	   endif
	   if(i.lt.ncomp) then
		j=j+1
		area(i)=theta(j)
		sa=sa+area(i)
	   endif
	enddo
	area(ncomp)=1.0-sa
	goto 40
C
30	CONTINUE
C NOW OP/BURST FIT
	if(ifitype.ne.2) goto 40
	Sa=0.		!total area except for last comp
	j=0
	do 32 i=1,ncomp
	j=j+1
	amean(i)=theta(j)
c	rho(i)=1.0 - (1.0/theta(j))
	j=j+1
	if(i.eq.ncomp) goto 32		!last area not a parameter
	area(i)=theta(j)
	sa=sa+area(i)
32	continue
	area(ncomp)=1.0 - sa
c	goto 40
c
40	continue
c Now define strings to be displayed
c
c Column headings for THIST
c for num1=char*8
	if(ifitype.eq.1) text=grk//'t'//fnt//' (ms) '//' area (%)'
	if(ifitype.eq.2) text='  '//grk//'m'//fnt//'   '//' area (%)'
	if(ifitype.eq.3) text=
     & grk//' m'//fnt//' (pA)'//'  SD   area (%)'
c and now the numbers
	do 1 j=1,ncomp
	call FIXN(amean(j),m)		!m sig figs after dec point
	x=ROUND(amean(j),m)		!round value appropriately
c	call HFORMAT(x,6,m,0,num1)	!0='fixedf' -see \hgraph\consts.for
	n1=9 		!size of num1 etc
	call DCFORMAT(x,n1-2,m,num1)
	y=100.*area(j)
	y=ROUND(y,1)   			!round value appropriately
	call DCFORMAT(y,n1-2,1,num2)
	nl=NBLANK(text)
	if(ifitype.le.2) then
	  text=text(1:nl)//nline//num1//num2
	else if(ifitype.eq.3) then
	  if(isdcon.eq.3.and.j.gt.1) then
		num3='   "    '
	  else
		call FIXN(SD(j),m)
		call DCFORMAT(ROUND(sd(j),m),n1-2,m,num3)
		if(isdcon.eq.2) then
		   n=NBLANK(num3)
		   n=n+1	!room for asterisk
		   if(n.gt.n1) n=n1
		   num3(n:n)='*'	!asterisk to indicate constraint
		endif
	  endif
	  text=text(1:nl)//nline//num1//num3//num2
	endif
1	continue
	goto 99
c
c Section for ifitype=4 (exponentials fitted to jumps)
c and ifitype=6,-6
c Use array amean=tau and area=amplitude
60	continue
	j=0
	do 61 i=1,ncomp
	j=j+1
	amean(i)=theta(j)       !tau
	j=j+1
	area(i)=theta(j)        !amplitude
61	continue
	Yinf=theta(2*ncomp+1)
c
	if(ifitype.eq.4) text=grk//'t'//fnt//' (ms) '//'ampl (pA)'
	if(ifitype.eq.6) text=grk//'t'//fnt//'      '//'amplitude'
	if(ifitype.eq.-6) text='''mean'''//fnt//'   '//'amplitude'
c and now the numbers
	do 62 j=1,ncomp
	call FIXN(amean(j),m)		!m sig figs after dec point
c	if(caplock()) print 701,j,amean(j),m
c701	format(' j,amean(j),m = ',i8,g13.6,i8)
	x=ROUND(amean(j),m)		!round value appropriately
	call DCFORMAT(x,n1-2,m,num1)
	y=area(j)				!=amplitude
	call FIXN(y,m)		!m sig figs after dec point
	y=ROUND(y,m)   			!round value appropriately
	call DCFORMAT(y,n1-1,m,num2)
	nl=NBLANK(text)
	text=text(1:nl)//nline//num1//num2
62	continue
c lastly add Y(inf); NB in Hgraph greek ';' is infinity sign
	nl=NBLANK(text)
	call FIXN(Yinf,m)		!m sig figs after dec point
	x=ROUND(Yinf,m)		!round value appropriately
	call DCFORMAT(x,n1-1,m,num1)
	nl=NBLANK(text)
	text=text(1:nl)//nline// 'Y('//grk//';'//fnt//') = '//num1
	goto 99
c
c Section for Lorentzians ifitype=5,-5
70	continue
c Use array amean=fc, converted to tau, and area=G(0) converted to area
	i=0
	do 71 j=2,2*ncomp,2
	i=i+1
	G0=theta(j-1)*1.e24
	fc=theta(j)
	area(i)=1.570796*G0*fc			!const=pi/2; area (pA*pA)
	amean(i)=1000./(6.283185*fc)		!tau (ms) from fc
71	continue
	if(ifitype.eq.-5) Yline=theta(2*ncomp+1)*1.e24	!line (pA**2/Hz)
c
	text=grk//'t'//fnt//' (ms) '//' Var (pA'//sup//'2'//sub//')'
c and now the numbers
	do 72 j=1,ncomp
	call FIXN(amean(j),m)		!m sig figs after dec point
	x=ROUND(amean(j),m)		!round value appropriately
	call DCFORMAT(x,n1-2,m,num1)
	call FIXN(area(j),m)		!m sig figs after dec point
	x=ROUND(area(j),m)		!round value appropriately
	call DCFORMAT(x,n1-2,m,num2)
	nl=NBLANK(text)
	text=text(1:nl)//nline//num1//num2
72	continue
c lastly add line if ifitype=-5
	if(ifitype.ne.-5) goto 99
	nl=NBLANK(text)
	if(yline.ge.1.e-3) then
	   call FIXN(Yline,m)		!m sig figs after dec point
	   x=ROUND(Yline,m)		!round value appropriately
	  call DCFORMAT(x,n1-1,m,num1)
	   nl=NBLANK(text)
	   text=text(1:nl)//nline//'G('//grk//';'//fnt//') (pA'
     &   //sup//'2'//sub//'s)='//num1
	else
	   do 5 i=1,9
		n=i
		yline=yline*10.
		if(yline.gt.1.0) goto 6
5	   continue
6	   y=ROUND(yline,2)   		!round eg to 2.31
	   call DCFORMAT(y,n1-3,1,num1)
	   nl=NBLANK(text)
	   nl1=NBLANK(num1)
	   text=text(1:nl)//nline//'G('//grk//';'//fnt//') (pA'
     &   //sup//'2'//sub//'s)='//num1(1:nl1)//dup//'x'//fnt//'10'//
     &   sup//'-'//char(n+48)//sub
	endif
	goto 99
c
c	if(debug()) print 2,nl,text(1:80)
c2	format(' nl= ',i8,/,a80)
c	nl=nblank(text)
c	ix=3500
c	iy=3500
c	call JUSTIFYSTRING(ix,iy,text(1:nl),0,isize,center,above)
c	call ENDPLT
c
99	continue
c Restore ifitype (e.g. for plot queue)
	if(isdcon.eq.2) then
	   ifitype=31
	else if(isdcon.eq.3) then
	   ifitype=32
	endif
	i=NBLANK1(text)	!insert ASCII 0 to end string
	DEALLOCATE(area,amean,sd)
	RETURN
	END


	subroutine FIXN(y,m)
c To fix number of sig figs after decimal point in parameter-value
c display in VHIST
	x=abs(y)
	if(x.lt.0.1) m=4	   				!eg for format F7.4
	if((x.ge.0.1).and.(x.lt.1.)) m=3	   		!eg for format F7.3
	if((x.ge.1.).and.(x.lt.10.)) m=2		!eg for format F7.2
c	if((x.ge.10.).and.(x.lt.100.)) m=1		!eg for format F7.1
c	if(x.ge.100.) m=0					!eg for F7.0
	if((x.ge.10.).and.(x.lt.1000.)) m=1		!eg for format F7.1
	if(x.ge.1000.) m=0				!eg for F7.0
	if(x.eq.0.0) m=1
	return
	end


