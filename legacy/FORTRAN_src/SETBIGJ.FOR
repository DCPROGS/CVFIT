	subroutine SETBIGJ(j,j0)
	common/setbig/jrep
c To recycle colours, line thickness values etc when there more than 10 curves
c  Modified 12/30/96 02:45pm so colours repeat after jrep, rather than 10
c (eg if 6 sets on display then set jrep=6 so if 2 sets of calculated lines
c  are drawn the second set will be same colour as the first)
	if(jrep.gt.0.and.jrep.le.10) then
	   jmax=jrep
	else
	   jmax=10
	endif
c
	if(j.le.jmax) then
	   j0=j
	else
	   j0=1+mod(j-1,jmax)	!so j=11-20, 21-30 all -> j0=1-10 etc
	endif
	RETURN
	end

