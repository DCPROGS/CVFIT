	subroutine DRAWGUES(xmin,xmax,ncalc,logx,logy,logity,y0,yinf,
     & Xcal,Ycal,ymin,ymax,nsfit,juse,noguess,njset,ndc1,ndimc,ic)
c=======================================================================
	real*4 xcal(ndc1,ndimc),ycal(ndc1,ndimc)
	integer juse(njset)
	logical logx,logy,noguess(100),logity
c
c To draw curves for initial guesses on graph after leaving VPLOT
c Code taken from VPLOT5 (after label 303)
	linetype=0		!continuous
	x0=1.e-36		!smallest value for logs
	do m=1,nsfit
	   j=juse(m)
	   if(noguess(j)) goto 221	!skip this set
	   call LINWID(0.20)
	   call lincols(ic,0)		!colour for jth calc curve
	   Jflag=0
	   do k=1,ncalc
		xv=xcal(k,j)
		yv=ycal(k,j)
		if(logx) then
		   if(xv.gt.x0) then
			xv=alog10(xv)
		   else
			goto 220
		   endif
		endif
		if(logy) then
		   if(yv.gt.x0) then
			yv=alog10(yv)
		   else
			goto 220
		   endif
		endif
		x0=1.e-36		!smallest value for taking logs
		if(logity) then
		   if(abs(yinf-yv).lt.1.e-30) goto 220	!for case of yv=yinf
		   if(yv.gt.yinf-x0) goto 220	!skip values above ymax!
		   if(yv.lt.y0+x0) goto 220	!skip values below ymin!
		   yv=(yv-y0)/(yinf-yv)	! Hill scale
		   if(yv.gt.x0) then
			yv=alog10(yv)
		   else
			goto 220
		   endif
		endif
		if(jflag.eq.0) then
		   call gramov(xv,yv)	!1st point in range
		else
		   if(xv.le.xmin.or.xv.ge.xmax.or.
     &		yv.lt.ymin.or.yv.gt.ymax) then
			call graMOV(xv,yv)	!1st point in range
		   else
                  call broken(linetype)
			call gralin(xv,yv)	!1st point in range
		   endif
		endif
		jflag=1
220		continue
	   enddo
221	   continue
	enddo
	RETURN
	end


