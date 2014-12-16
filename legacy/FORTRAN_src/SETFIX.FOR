	subroutine SETFIX(nfix,jfix,kmax,kfit,theta,titlep,alpha,ndth)
c To choose which params to fix in CVFIT -another version
	integer JFIX(ndth)
	real*4 theta(ndth)
	character*20 TITLEP(ndth)		!names of params for fit
	logical alpha,help
	character*11 cnum,cnum1
	character defolt*30,cans*30		!to hold text & result of QDIALOG
	character*1 ans
c For popmenu
	character*27 title 			!for popmenu
	common/popvalm/nxlo,nylo,nyhi,ictm,ibkm,icfm,icupm	!values for popmenu calls
	allocatable strings,helps,thetasc
	character*55 strings(:),helps(:),thetasc(:)
	common/hlp/help				!for QDIALOG
c

	ALLOCATE(strings(kmax+1),helps(5),thetasc(kmax+1))
	if(alpha) goto 1
c
c
c NB can display only 20 parameters at once -otherwise must show popmenu
c nmenu times to show all parameters
	if(kmax.le.20) then
	   nmenu=1
	   kdisp=kmax
	   title='Fix/unfix parameter:'
	else
	   nmenu=1+ (kmax-1)/20
	   kdisp=20
	   kdisp1=20
	   klast=kmax - (nmenu-1)*kdisp
	   title='Fix/unfix parameter: menu 1'
	endif
c
	nmax=0
	npmax=0
	imenu=1	!=1,2,...,nmenu

	do imenu=1,nmenu
2	   call INTCONV(nfix,cnum)
	   call WDIALOG(1,CHARNB(cnum)//
     &  ' parameters already fixed: select more (or ESC when done)',11)
	   if(imenu.gt.1) then
		call INTCONV(imenu,cnum)
		title='Fix/unfix parameter: menu '//charnb(cnum)
		if(imenu.eq.nmenu) kdisp=klast
	   endif
	   do j=1,kdisp
		j1=j + (imenu-1)*kdisp1	!index for parameter (and titlep)
		call REALTOCH(theta(j1),cnum1,11)
		thetasc(j)=cnum1
		n1=NBLANK(thetasc(j))
		np=NBLANK(titlep(j1))
		if(n1.gt.nmax) nmax=n1
		if(np.gt.npmax) npmax=np
	   enddo
	   nc=5+npmax+2+nmax
	   do j=1,kdisp
		j1=j + (imenu-1)*kdisp1	!index for parameter (and titlep)
		call INTCONV(j1,cnum)
		strings(j)=cnum(1:2)//':  '//titlep(j1)(1:npmax)//
     &    '  '//thetasc(j)(1:nmax)
     		if(jfix(j1).eq.1) then
		   strings(j)=strings(j)(1:nc)//'   (fixed)'
		endif
	   enddo
	   if(nmenu.eq.1) then
	      strings(kdisp+1)='All parameters O.K.'
	   else
		j1=1 + (imenu-1)*kdisp1
		j2=j1+kdisp-1
		call INTCONV(j1,cnum)
		call INTCONV(j2,cnum1)
	      strings(kdisp+1)='Parameters '//charnb(cnum)//
     &	' to '//charnb(cnum1)//' O.K.'
	   endif
	   nval=kdisp+1
c Default=accept
	   jline=kdisp+1
	   nhelp=3
	   helps(1)='Select a parameter to fix/unfix.   '
	   helps(2)='When all selected correctly hit ESC'
	   helps(3)='or select the last line, to finish.'
	   call POPMENU(nxlo,nylo,nyhi,strings,nval,ictm,icfm,icupm,
     &      ibkm,title,helps,nhelp,jline,ans,ival)
	   itemp=jline + (imenu-1)*kdisp1 !index for parameter (and titlep)
         if(jline.ge.1.and.jline.le.kdisp) then	!ESC not hit
		if(jfix(itemp).eq.1) then
		   jfix(itemp)=0	!unfix
		   nfix=nfix-1
	      else
	   	   call INTCONV(itemp,cnum)
		   call DEFOLTr(theta(itemp),defolt)
		   call QDIALOG(1,'Fixed value for parameter '//
     &		CHARNB(cnum),defolt,11,cans)
	         call GETINPr(cans,theta(itemp))
		   jfix(itemp)=1
		   nfix=nfix+1
		endif
	      goto 2	!fix another?
	   endif
	enddo		!goto next menu if any when ESC or iline=kmax+1 hit
c
	kfit=kmax-nfix
	DEALLOCATE(strings,helps,thetasc)
	RETURN
c
c Now non-graphics version to be done if nec
1	continue
	kfit=kmax-nfix
	DEALLOCATE(strings,helps,thetasc)
	RETURN
	end
