	subroutine NEWEQCV(nmod,kmax,kmax1,titlep,nsfit,juse,ifitmode,
     & theta,thetsav,setx,njset,ndth)
c Hill/Langmuir -separate fits with common maximum after simultaneous fit
c To alter equation to be fitted to nmod=23 or 24 in CVFIT in case of
c isepfit=3, following simultaneous fit with ifitmode=4 or 5.  Model
c is redefined here, ifitmode changed to 3, and initial guesses, based
c on the simultaneous fit results, put in thetsav(i,1) ready for use
c in CVSIMP (for these models fit done in one pass -only one set of init
c guesses so put in thetsav(i,1))
c Modified 12/31/96 09:00am for new vwesion in which Hill/Langmuir are
c now nmod=26,27 (not 1, -1 as before)
c In this case nmod=1 (Langmuir) or -1 (Hill) and want to have
c separate K1 ( and for Hill, separate NHill) for each curve but a
c common maximum for all, ie nmod=23,24 (see GETEQN)
c
	real*4 theta(100),thetsav(100),setx(njset)
	integer juse(njset)
	character*20 TITLEP(ndth)		!names of params for fit
c
c Put initial guesses into theta(i) before kmax etc altered
c (order of params as specified when titlep() defined below)
	Y0=thetsav(1) 		!Y(0)
	Ymax=thetsav(2)		!value of Ymax from parallel fit
	aka=thetsav(3)		!value of K1 from parallel fit
	hill=thetsav(4) 		!n(Hill)
	do 10 iset=1,nsfit
	 j=juse(iset)
	 ratio=1.0
	 if(ifitmode.eq.4.and.iset.gt.1) then
	   ratio=thetsav(kmax1+iset-1)		!from parallel fit
	 else if(ifitmode.eq.5) then
	   akb=thetsav(kmax)		!=KB (see GETEQN) -value from parallel fit
	   xB=setx(j)		!antag conc
	   ratio=(1.0 + xB/aKB)	!mult by dose ratio
	 endif
	 if(nmod.eq.26) then
	   theta(1)=Y0
	   theta(2)=Ymax		!init guess for nmod=23 fit
	   theta(iset+2)=aka*ratio
	 else if(nmod.eq.27) then
	   theta(1)=Y0		!init guess for nmod=24 fit
	   theta(2)=Ymax		!init guess for nmod=24 fit
	   m1=2*iset+1
	   theta(m1)=aka*ratio	!=theta(3,5,7,...)
	   theta(m1+1)=hill		!=theta(4,6,8,...)
	 endif
10	continue
c
c Redefine model
	if(nmod.eq.26) then
	   nmod=23		!Langmuir with common maximum
	   kmax=nsfit+2	!Y(0), K for each, + common max
	   titlep(1)='     Y(0) '
	   titlep(2)='     Ymax '
	   do 55 m=1,nsfit
		j=juse(m)		!set number used
	      titlep(m+2)='       K'//char(j+48)//' '	!='K3','K5',...
	      if(j.eq.10) titlep(m+1)='      K10 '
55	   continue
	endif
	if(nmod.eq.27) then
	   nmod=24		!Hill with common maximum
	   kmax=2*nsfit+2	!K,nH for each + common max and min
	   titlep(1)='     Y(0) '
	   titlep(2)='     Ymax '
	   do 56 m=1,nsfit
		j=juse(m)		!set number used
		m1=2*m+1
	      titlep(m1)='       K'//char(j+48)//' '	!='K3','K5',...
	      titlep(m1+1)='      nH'//char(j+48)//' '	!='nH3','nH5',...
	      if(j.ge.10.and.j.le.19) then
	         titlep(m1)='      K1'//char(j-10+48)//' '
	         titlep(m1+1)='     nH1'//char(j-10+48)//' '
		endif
56	   continue
	endif
c
c Redefine ifitmode,kmax1
	ifitmode=3				!change ifitmode if not already 3
	kmax1=kmax
c
	RETURN
	end

