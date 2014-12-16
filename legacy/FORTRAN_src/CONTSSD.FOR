	subroutine CONTSSD(y1,y2,complex,x,y,w,n,SSD,aK)
c To calculate the SSD contours for fit of hyperbola
c Returns the two values of Ymax in y1,y2 (if complex=false), given
c data in x,y,n and specified values of K and SSD.
	real*4 x(n),y(n),w(n)
	logical complex
c
c  Calculate coefficients of quadratic for give K value (aK) and given SSD value
c
	s=0.
	do i=1,n
	   s=s + w(i)*(x(i)/(x(i)+aK))**2
	enddo
	den=s
c
	s=0.
	do i=1,n
	   s=s + w(i)*(x(i)*y(i)/(x(i)+aK))
	enddo
	b=-2*s
c
	s=0.
	do i=1,n
	   s=s + w(i)*y(i)*y(i)
	enddo
c Coefficients are
	c= (s - SSD)/den
	b=b/den
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

