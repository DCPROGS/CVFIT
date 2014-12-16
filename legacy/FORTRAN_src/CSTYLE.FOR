	subroutine CSTYLE(ncurvd,icurvd,isym,ijoin,
     & ncurvc,icurvc,iline,syms,icol,ndimd,ndimc,imode)
c To read/write curve line types/colours/symbols so that a graph can
c more easily be repeated in the same style in CVFIT etc
c Sometimes there is a very large number of data sets, so keep as in plot queue
c i.e. keep only ncurvd, ncurvc values regardless of which data sets (j values
c in xval(i,j)) are plotted.
c NB call does not alter inoout values of ncurvd,ncurvc,ndimc or ndimc
c	imode=0 to write
c	imode=1 to read
c
	integer icurvd(ndimd),isym(ndimd),ijoin(ndimd)
	real*4 syms(ndimd)
	integer icurvc(ndimc),iline(ndimc)
	integer icol(100)
	ALLOCATABLE::isym1,ijoin1,syms1,iline1
	integer isym1(:),ijoin1(:),iline1(:)
	real*4 syms1(:)
	ALLOCATABLE::icurvd1,icurvc1
	integer icurvd1(:),icurvc1(:)
	logical present
c
	if(imode.eq.1) goto 10
c Write to disc
c Allocate arrays for writing to disc
	if(ncurvd.eq.0) then
	   nd=1			!must allocate anyway
	else
	   nd=ncurvd
	endif
	if(ncurvc.eq.0) then
	   nc=1			!must allocate anyway
	else
	   nc=ncurvc
	endif
	ALLOCATE(isym1(nd),ijoin1(nd),syms1(nd))
	ALLOCATE(iline1(nc))
c Data
	do j1=1,ncurvd
	   j=icurvd(j1)
	   isym1(j1)=isym(j)
	   ijoin1(j1)=ijoin(j)
	   syms1(j1)=syms(j)
	enddo
c Calc curves
	do j1=1,ncurvc
	   j=icurvc(j1)
	   iline1(j1)=iline(j)
	enddo
c Write
	OPEN(unit=13,file='VPLOT.STY',status='UNKNOWN',
     &		access='TRANSPARENT')
	write(13,rec=1) ncurvd,ncurvc,nd,nc,
     & (icurvd(i),i=1,ncurvd),(icurvc(i),i=1,ncurvc),
     &  isym1,ijoin1,syms1,iline1,icol
	CLOSE(unit=13)
	goto 99
c
c Read back from disc
10	continue
	INQUIRE(file='VPLOT.STY',exist=present,flen=nlen)
	if(.not.present) goto 99
	OPEN(unit=13,file='VPLOT.STY',status='UNKNOWN',
     &		access='TRANSPARENT')
	read(13,rec=1) ncurvd1,ncurvc1,nd,nc
	ALLOCATE(isym1(nd),ijoin1(nd),syms1(nd))
	ALLOCATE(iline1(nc))
	ALLOCATE(icurvd1(nd),icurvc1(nc))
	read(13,rec=1) ncurvd1,ncurvc1,nd,nc,
     & (icurvd1(i),i=1,ncurvd1),(icurvc1(i),i=1,ncurvc1),
     &  isym1,ijoin1,syms1,iline1,icol
	CLOSE(unit=13)
c Put data back into original arrays
c Care! If these are from a different data set, there could be too many
c values to fit into the original arrays
	do j1=1,ncurvd1
	   icurvd(j1)=icurvd1(j1)
	   j=icurvd1(j1)
	   if(j.lt.ndimd) then	!do not exceed ndimd for current program
		isym(j)=isym1(j1)
		ijoin(j)=ijoin1(j1)
		syms(j)=syms1(j1)
	   endif
	enddo
c Calc curves
	do j1=1,ncurvc1
	   icurvc(j1)=icurvc1(j1)
	   j=icurvc1(j1)
	   if(j.lt.ndimc) then	!do not exceed ndimd for current program
		iline(j)=iline1(j1)
	   endif
	enddo
	DEALLOCATE(icurvd1,icurvc1)
c
99	continue
	DEALLOCATE(isym1,ijoin1,syms1,iline1)
	RETURN
	end
