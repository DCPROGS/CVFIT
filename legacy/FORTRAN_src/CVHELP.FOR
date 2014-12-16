	subroutine CVHELP(ihelp,titmod,imod)
c To define model in CVFIT when in graphics mode (ihelp=1)
c (and other popup help?)
c Lahey/Ioana colours:
c 0=black; 1=dark blue; 2=green; 3=light blue; 4=red; 5=purple
c 6=brown; 7=white; 8=pale white (grey); 9=dark blue -bright; 10=green -bright
c 11=light blue -bright; 12=red -bright; 13=purple -bright; 14=yellow -bright
c 15=white -bright
	character*60 titmod(30)		!names of models
	character*56 helpst(20)		!help strings for poptext
	character*60 strings(30)
      character*55 helps(5)	!for popmenu
	character*20 title
      character*1 charout 		!for popmenu
c
	ixloh=60    !position for poptext box
	iyhih=430
c
2	continue
	SELECT CASE(ihelp)
	CASE(1)
	   nxlo=100
	   nylo=-1
	   nyhi=470
	   ictm=15		!text white
c	   ibkm=8		!background dark grey
	   ibkm=1		!background dark blue
	   icfm=4		!frame/title red
c	   icupm=12		!upper case red
	   icupm=15		!upper case white (too many options to use letter)
	   TITLE='CHOOSE EQUATION: '
c	   do i=1,9
c		strings(i)=titmod(i)
c	   enddo
	   do i=1,3
		strings(i)=titmod(i)
	   enddo
	   strings(4)=' (+option to normalise several sets with common X)'
	   do i=4,9
		strings(i+1)=titmod(i)
	   enddo
	   strings(11)=titmod(22)
	   strings(12)=titmod(23)
	   strings(13)='(*) Equations for P(open) curves . . .'
	   strings(14)='(*) Equations for I/V curves . . .'
	   nval=14
	   nline=nval-2	!number of lines exc (*) options
	   helps(1)='Choose equation to be fitted'
	   nhelp=1
	   call POPMENU(nxlo,nylo,nyhi,strings,nval,ictm,icfm,icupm,
     &	 ibkm,title,helps,nhelp,iline,charout,ival)
	   if(iline.ge.1.and.iline.le.nval) then	!iline=0 for ESC=cancel
		if(iline.ge.1.and.iline.le.10) then
		   if(iline.ge.1.and.iline.le.3) then
			imod=iline
		   else if(iline.ge.4.and.iline.le.10) then
			imod=iline-1
		   endif
		   RETURN
		else if(iline.eq.11) then
		   imod=22
		   RETURN
		else if(iline.eq.12) then
		   imod=23
		   if(imod.eq.23) then	!not in use
			call BELL(1)
			goto 2
		   endif
		   RETURN
		else
		   goto 1	!Popen and i/v options
		endif
	   else
		RETURN		!if ESC hit
	   endif
1	   continue
c Popen and i/v options
	   nxlo=200
	   nylo=-1
	   nyhi=450
	   ictm=15		!white
c	   ibkm=8		!background dark grey
	   ibkm=2		!background dark green
	   icfm=4		!frame/title red
c	   icfm=14		!frame yellow
c	   icupm=12		!upper case red
	   icupm=15		!upper case white (too many options to use letter)
c
	   if(iline.eq.nline+1) then		!Popen options
		i1=1
		do i=10,13
		   strings(i1)=titmod(i)
		   i1=i1+1
		enddo
		do i=15,17
		   strings(i1)=titmod(i)
		   i1=i1+1
		enddo
		strings(i1)=titmod(24)	!n bindings +opening (nmod=30)
		i1=i1+1
		strings(i1)=titmod(25)	!MWCn (nmod=31)
		i1=i1+1
		strings(i1)=' Return to main menu'
		nval=i1
	      call POPMENU(nxlo,nylo,nyhi,strings,nval,ictm,icfm,icupm,
     &	 ibkm,title,helps,nhelp,iline,charout,ival)
		if(iline.ge.1.and.iline.le.nval) then	!iline=0 for ESC=cancel
		   if(iline.ge.1.and.iline.le.4) then
			imod=iline+9	!imod=10-14
		   else if(iline.ge.5.and.iline.le.7) then
			imod=iline+10	!imod=15-17
		   else if(iline.eq.8) then
			imod=24	!n bindings +opening (nmod=30)
		   else if(iline.eq.9) then
			imod=25	!MWCn (nmod=31)
		   else if(iline.eq.nval) then
			goto 2
		   endif
		   RETURN
		endif
	   else if(iline.eq.nline+2) then	!i/v options
		strings(1)=titmod(14)
		i1=2
		do i=18,21
		   strings(i1)=titmod(i)
		   i1=i1+1
		enddo
		strings(i1)=' Return to main menu'
		nval=i1
	      call POPMENU(nxlo,nylo,nyhi,strings,nval,ictm,icfm,icupm,
     &	 ibkm,title,helps,nhelp,iline,charout,ival)
		if(iline.ge.1.and.iline.le.nval) then	!iline=0 for ESC=cancel
		   if(iline.eq.1) then
			imod=14
		   else if(iline.ge.2.and.iline.le.5) then
			imod=iline+16	!imod=18-21
		   else if(iline.eq.nval) then
			goto 2
		   endif
		   RETURN
		endif
	   endif
c
	CASE(2)
	   i=1
	   helpst(i)=
     &   'Degree 1=straight line; 2=quadratic, . . . ,10         '
	   i=i+1
	   helpst(i)=
     &   'Fit Y = Y(0) + a1*x + a2*x^2 +. . .                    '
	   i=i+1
	   helpst(i)=
     &   'Disjoint polynomial has degree above some specified    '
	   i=i+1
	   helpst(i)=
     &   '  value x = x1, different from that below x1: in this '
	   i=i+1
	   helpst(i)=
     &   '  case fits: Y = Y(0) + a1*(x-x1) + a2*(x-x1)^2 +. . .'
	   i=i+1
	   helpst(i)=
     &   ' (To keep 2nd deriv continuous at the join fix a2 = 0)'
	   i=i+1
	   helpst(i)=
     &   'Now specify the value:                                '
	   nhelp=i
c	   subroutine POPTEXT(ixlo,iylo,iyhi,strings,n,ic,ibk,icf)
	   call POPTEXT(ixloh,-1,iyhih,helpst,nhelp,12,15,12)

	CASE(3)
	   i=1
	   helpst(i)=
     &   'Binding inhibition curve with form of a single         '
	   i=i+1
	   helpst(i)=
     &   'declining Lnagmuir component (use Hill equation (eq.3) '
	   i=i+1
	   helpst(i)=
     &   'to check that Hill slope is near 1).                   '
	   i=i+1
	   helpst(i)=
     &   'Estimates the antagonist equilibrium constant directly '
	   i=i+1
	   helpst(i)=
     &   'by fitting Y = Y(inf) + Y(0)(1+cA)/(1+cA+xB/KB),       '
	   i=i+1
	   helpst(i)=
     &   '  where Y=binding, xB=antagonist conc, KB=antagonist   '
	   i=i+1
	   helpst(i)=
     &   '  equilibrium constant, and cA=xA/KA is the normalised '
	   i=i+1
	   helpst(i)=
     &   '  conc of labelled ligand (conc=xA, eq constant=KA).   '
	   nhelp=i
	   call POPTEXT(ixloh,-1,iyhih,helpst,nhelp,12,15,12)

	CASE(4)
	   i=1
	   helpst(i)=
     &   'Define index of cooperativity in binding reaction as   '
	   i=i+1
	   helpst(i)=
     &   'r = K(i-1)/K(i) so R>1 implies positive cooperativity. '
	   nhelp=i
	   call POPTEXT(ixloh,-1,iyhih,helpst,nhelp,12,15,12)

	END SELECT
c
	END

