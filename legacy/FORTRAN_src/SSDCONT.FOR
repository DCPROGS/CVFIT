	program PROG SSDCONT
	logical discprt
	logical pon,slock,debug,caplock,mono
	character*11 cdate,ctime
C For DCMENU
	integer Lb(18)
	character*78 text(18)
c For IV mouse routines
	logical mouse_on,mous_set
	common/mousval/mouse_on,nbutton
c
	common/DMENU/csizb,ifonb,csize,ifont2,nboxlast,nblast           !for DCMENU
	common/dp/discprt
	character*40 mtitle*40,filnam*32,prtport*4      !for WINPRINT
	common/dpp/filnam,prtport,ndisc,jcol,mtitle !for WINPRINT,ENDPRINT,DISCNUM
C
c define function
	pon()=slock()
	debug()=caplock()
c
	filnam='PROG.PRT'
	call WINPRINT   !print file control
      OPEN(unit=7,file=prtport,iostat=nerr)             !open printer
	print 1
	if(pon()) write(7,1)
	if(discprt) write(8,1)
1       FORMAT(' PROG: TO CALCULATE SSD CONTOURS FOR HYPERBOLA',/)
	call DATE1(cdate)               !DC subroutine
	call TIME(ctime)
	print 2,cdate,ctime(1:8),mtitle
	if(pon()) write(7,2) cdate,ctime(1:8),mtitle
	if(discprt) write(8,2) cdate,ctime(1:8),mtitle
2       format(' Date of analysis: ',a11,/,' Time of analysis: ',a8,/,
     & '   Machine = ',a40)
	print 3
3       format(
     & ' SCROLL LOCK on for printing'/
     & ' CAPS LOCK on for debugging')

=============
	call SETMOUSE()         !define values in common\mousval\ (in IVLIB)

============
	mono=jcol.eq.0  !monochrome screen (NB, after WINPRINT)
c defs for DCMENU
c       nbox=10 !number of boxes for DCMENU
	icol1=14        !yellow text/border for boxes in DCMENU
	icol2=8 !grey background for boxes
	icol3=12        !red text/border
	icol4=15        !white fill
	if(mono) then
	   icol3=15             !white text/border
	   icol4=0              !black fill
	endif
	ifonb=0
	csizb=1.7
	ifont=4
	csize=2.
	nbox=10                 !define to prevent error when deleting dcmenu
	nbline=5                !define to prevent error when deleting dcmenu
c
	ivxlo=0         ! for VIEWPORT
	ivxhi=10000
	ivylo=0
	ivyhi=7000
	ixlo=0          ! for graphboundary
	ixhi=10000
	iylo=0
	iyhi=7000
===
	nbox=15 !number of boxes for DCMENU
	call SETLBOX(nbox,Lb,1)
	call NUMSET                     !set num lock on ready for response
	call DCMENU(0,4,Lb,text,0,0)            !delete all
	call DCMENU(nbox,5,Lb,text,icol1,icol2)
===
=========REAL STUFF STARTS HERE (see SSDLOB.mcd)

	subroutine CONTSSD(y1,y2,complex,x,y,n,SSD,aK)
c To calculate the SSD contours for fit of hyperbola
c Returns the two values of Ymax in y1,y2 (if complex=false), given
c data in x,y,n and specified values of K and SSD.
	real*4 x(n),y(n)
	logical complex
c
c  Calculate coefficients of quadratic for give K value (aK) and given SSD value
c
	s=0.
	do i=1,n
	   s=s + (x(i)/(x(i)+aK))**2
	enddo
	den=s
c
	s=0.
	do i=1,n
	   s=s + x(i)*y(i)/(x(i)+aK)
	enddo
	numb=-2*s
c
	s=0.
	do i=1,n
	   s=s + y(i)*y(i)
	enddo
c Coefficients are
	c= (s - SSD)/den
	b=numb/den
c and solutions are
	disc=b*b - 4*c
	if(disc.gt.0) then
	   y1=0.5*(-b+sqrt(disc))
	   y2=0.5*(-b-sqrt(disc))
	   complex=.false.
	else
	   complex=.true.
	endif
	RETURN
	end

=========
	call ENDPRINT
	end

