	real function FACT(n)
c Returns factorial in floating point
	if(n.eq.0.or.n.eq.1) then
	   FACT=1.0
	else
	   if=1
	   do i=2,n
		if=if*i
	   enddo
	   FACT=float(if)
	endif
c
	RETURN
	end
