	real x(5),y(5),w(5)
	logical complex
c
	n=5
	do i=1,n
	   w(i)=1.
	enddo
	x(1)=2.5
	y(1)=5.576
	x(2)=5.
	y(2)=7.282
	x(3)=10.
	y(3)=12.521
	x(4)=20.
	y(4)=16.138
	x(5)=40.
	y(5)=23.219
c
1	continue
	ssd=20.
c
	call CONTSSD(y1,y2,complex,x,y,w,n,SSD,aK)
c
	goto 1
	end

