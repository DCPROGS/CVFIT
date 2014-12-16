	program  TGRP
c To test GETGROUP
	real xobs(10,2),yobs(10,2),w(10,2)
	ALLOCATABLE::xgrp,ybar,sdy,sdm,ngrp,nxval
	real*4 xgrp(:,:),ybar(:,:),sdy(:,:),sdm(:,:)
	integer ngrp(:,:),nxval(:)
	integer nj(2)
c
	niobs=10
	njset=2

	xobs(1,2)=1
	xobs(2,2)=1
	xobs(3,2)=1
	xobs(4,2)=1
	xobs(5,2)=3
	xobs(6,2)=3
	xobs(7,2)=5
	xobs(8,2)=5
	xobs(9,2)=5
c
	w(1,2)=1.
	w(2,2)=1.
	w(3,2)=1.
	w(4,2)=1.
	w(5,2)=1.
	w(6,2)=1.
	w(7,2)=1.
	w(8,2)=1.
	w(9,2)=1.
	do i=1,9
	   w(i,2)=4.
	enddo
c
	yobs(1,2)=1
	yobs(2,2)=2
	yobs(3,2)=3
	yobs(4,2)=4
	yobs(5,2)=10
	yobs(6,2)=11
	yobs(7,2)=20
	yobs(8,2)=21
	yobs(9,2)=22
c
c Result for above is:
c SStot =   594.222      SScol =   586.722      SSwg =   7.50000
c  group  n(j)  X value      mean       sdm
c   1       4  1.00000      2.50000     0.645497
c   2       2  3.00000      10.5000     0.500000
c   3       3  5.00000      21.0000     0.577350
c
	nj(1)=0
	nj(2)=9
	jset=2
	ALLOCATE(xgrp(niobs,njset),ybar(niobs,njset),sdy(niobs,njset),
     & sdm(niobs,njset),ngrp(niobs,njset),nxval(niobs))
c
	call GETGROUP(xobs,yobs,w,jset,nj,niobs,njset,nxval(jset),
     &	xgrp,ybar,sdy,sdm,ngrp,SSwg,ndf,SStot,SScol,ymean,n0)
	if(n0.gt.0) then
	   print 3,n0
3	   format(' ',i3,' values omitted because they have zero weight')
	endif
	print 2,SStot,SScol,SSwg
2	format(
     & '  SStot = ',g13.6,'  SScol = ',g13.6,'  SSwg = ',g13.6,/,
     & ' group  n(j)  X value      mean(y)     sd            sdm')
	do j=1,nxval(jset)
	   print 1,j,ngrp(j,jset),xgrp(j,jset),ybar(j,jset),
     &    sdy(j,jset),sdm(j,jset)
1	   format(i3,i8,4g13.6)
	enddo
	DEALLOCATE(ybar,sdm,ngrp,nxval)
	end
