      subroutine CVDISK(Xobs,Yobs,W,nset,titlep,setvar,iw,
     & ilabel,titlef,titled,titlex,titley,ifile,nfile,FOPEN,ndev,hdisp,
     & prt,readini,iesc,readat,ndisc,infil,ifile1,
     & niobs,njset,nj,setx,iver)
c  In version 1003 and later titlep() is no longer kept in CVDAT -it is not
c used!
c Modif 09/24/97 04:28pm to make ioffset=2048 (may get more than 1024 bytes
c in the first record if files is big)
c
c=	real YOBS(100,10),XOBS(100,10),W(100,10),DATCOP(3000)
	real*4 xobs(niobs,njset),yobs(niobs,njset),w(niobs,njset)
	ALLOCATABLE datcop,datcop1			!for INWIND
	real datcop(:),datcop1(:)
	ALLOCATABLE titled1,nj1,setx1
	character*60 TITLED1(:)
	integer nj1(:)
	real setx1(:)
c=	real setx(10),setx1(10),datcop1(3000)
	integer jstrec(100),lstrec(100)
	integer*2 jstrec2(200),lstrec2(200)		!for queue
c=	integer nj(10),nj1(10)
	character*60 titled(njset)
	integer nj(njset)
	real setx(njset)
c=	character*60 titled1(10)			!data set names
	character*40 titleds(10),titleds1	!shorter set titles for iver=>1002
	character*60 titlef,titlef1		!file names (for iver=1002)
	character*10 TITLEP(20),titlep1(20)		!names of params for fit
	character*40 titlex,titley,titlex1,titley1
 	character*1 ans,UC,getch
	character ndev*2,ndev1*2,filnam*12,infil*33
	logical slock,pon,debug,caplock,setvar,fopen,filspec,readat
	logical discprt,hdisp,present,hdisp1,setvar1,prt,readini
	logical cluster,student
	common/dp/discprt
	common/user/student,cluster,iwindows
c
c To do disc read/write for CVFIT (based on method in VPLQ1 for minimum
c disk space usage)
c Modif for iver=1003 has quite different (more economical/flexible) file
c format based on transparent file. N.B. CVDAT.DAT can be read for any
c iver number, but is always written in iver=1003 format, so for
c writing to disk the file, if open, is always closed and reopened as
c a transparent file.
c Modified 07/10/92 07:11pm for iver=1002 that has title for whole file
c as well as titles for each set.  To make room in the fixed record length
c must make set titles only 40 characters rather than 60, but at present
c the set titles are copied into titled()*60 to avoid changes in rest
c of program (and to allow iver=1001=orig version to be read in full)
c
c This subroutine does not close unit=12- must set FOPEN and close file
c outside it
c Modif 08/07/91 08:35pm to cope with case where addition of more set(s) to
c  a file (other than the last file) makes that file need an extra record
c  which would overwrite the next file.  In this case all files above that being
c  written to are read from disc, and written back further up.
c  Note that one set may have up to 300 reals=1200 bytes=more than one record.
c Modif 07/26/91 05:32pm so number of file which data comes from (if it was
c  read from disc) is input in IESC for use as default file to write to
c  if adding a set etc). If iesc=0 on entry then NEXT unused file number
c  is default to write to. -NOW ifile1 from .ini is a parameter -iesc not used
c
c Input:
c IFILE=>1 to read values for file #ifile (which is asked for here so set
c	any value eg ifile=1 in call (and also lists all titles if req)
c IFILE=0 to write values to file number decided here (and returned in IFILE)
c IFILE<0 to write values to file number -IFILE (sign restored before return)
c FOPEN=true if data file already open
c
c
	pon()=slock()
	debug()=caplock()
c
101	format(a1)
c	prt=pon()		!to avoid interrupt 11 problems! -now an arg
c Check which disc to use and open file
	iesc=0
	filspec=ifile.lt.0
	if(IFILE.le.0) goto 2	!write the data
	if(FOPEN) goto 5
c
	call DCASK('Read data from disk','Y',ans)
	readat=ans.eq.'Y'
c
c Read initialisation file, cvfit.ini
c ===reading of .ini now moved to main prog, so can keep default array size!
	if(.not.readat) RETURN		!type in the data
	if(.not.readini) then    !.INI not present and so not read in main prog
	   if(cluster) then
		infil='o:\cvfit\*.dat'
	   else if(ndisc.gt.1) then
		filnam='CVDAT.DAT'
	      call DISCNUM1(id,ndev,-1,0)	!ask for winchester or floppy
	      if(id.eq.-1) then
		   iesc=99
		   RETURN			!abort if ESC hit
	      endif
	      INFIL=ndev//FILNAM
	   else
		infil='CVDAT.DAT'		!default for ndisc=1
	   endif
	endif
c
31	continue
	call TITENT0(
     &'Name and path of CVFIT data file:',infil,33,.false.)
	call CHKFILE(infil,nerr,33)
	if(.not.cluster) then
	   if(nerr.ne.0) then
		call BELL(1)
		print 311
311		format(' INVALID file name -try again')
		goto 31
	   endif
	else if(cluster) then
	  if(infil(1:9).ne.'o:\cvfit\') then
	     print*,'TRY AGAIN: something like  o:\cvfit\yourname.dat '
           goto 31
	  endif
	  if(nerr.ne.0.or.infil(10:11).eq.'*.') then
	     print*,'TRY AGAIN: something like  o:\cvfit\yourname.dat '
           goto 31
	  endif
	endif
c
71	continue	!jump here if INFIL already chosen in FINDFIL
	INQUIRE(file=infil,exist=present,flen=len)
	if(.not.present.or.len.eq.0) then
	   call BELL(3)
	   print 24,infil
24	   format(1x,a33,' does not exist')
	   present=.false.		!so asks for ndev
	   goto 31
	else
         OPEN(unit=12,file=INFIL,status='UNKNOWN',
     &   access='DIRECT',form='UNFORMATTED',recl=1024,err=30)
c NB recl is in bytes (4*value for PDP)
	   fopen=.true.
	endif
5	continue	!jump here if already open
c
c
c Section to read data  for file #ifile
1	continue
	read(12,rec=1) nfile,jstrec,lstrec,iver
	if(iver.le.1000.or.iver.gt.1100) iver=1001	!original version
c iver=version # to get format of data files. Original version
c now has iver=1001, If iver=1002 then an extra overall title (titlef(60))
c is stored for each file, as well as the data set titles (titled)
c iver=1003 has new organisation of CVDAT with transparent read/write!
	if(iver.ge.1003) then
	   CLOSE(unit=12)		!then reopen as transparent
         OPEN(unit=12,file=INFIL,status='UNKNOWN',access='TRANSPARENT')
	endif
c
	print 335,nfile
335	format(1x,i5,' files: list file titles [N] ? ')
	read 101,ans
	if(UC(ans).ne.'Y') goto 338
c In Lahey fortran file is not of pre-fixed size but has as many records
c as have been written
	do 518 ifile=1,nfile
	irec=jstrec(ifile)
cc	if(iver.ge.1003) irec=jstrec2(ifile)
c modif to cope with old data files that have different iver in same
c data file -check if nset=1-10 when read as if old, and set iver from this
c (overwrites the iver read from record #1)
	if(iver.lt.1003) then
	   read(12,rec=irec) (titled(k),k=1,10),nset
	   if(nset.ge.1.and.nset.le.10) then
		iver=1001		!old
	   else
		iver=1002		!new
	   endif
	else		!iver=> 1003
	   read(12,rec=irec) nset
	   if(nset.gt.njset-5) then
		call BELL(2)
		print 19,nset,njset
19	     format(' Number of sets = ',i4,' too large for array size',
     &	i4,/,' Allocate a larger size before proceeding')
		iesc=98
		CLOSE(unit=12)
		RETURN
	   endif
	endif
	if(iver.eq.1001) then
	   read(12,rec=irec) (titled(k),k=1,10),nset
         print 520,ifile,nset
         if(prt) write(7,520) ifile,nset
         if(discprt) write(8,520) ifile,nset
520	   format(' File number ',i3,':  contains ',i4,' data sets')
	else if(iver.eq.1002) then
	   read(12,rec=irec) titlef,titleds,nset
	   call SWAPTIT(titled,titleds,nset,njset,1)  !copy titleds into titled
         print 5201,ifile,titlef,nset
         if(prt) write(7,5201) ifile,titlef,nset
         if(discprt) write(8,5201) ifile,titlef,nset
5201	   format(' File no ',i3,': ',a60,/,' contains ',i4,' data sets')
	else if(iver.eq.1003) then
c      N.B. nj() etc already allocated in calling prog to size=njset
	   read(12,rec=irec) nset,titlef,(titled(j),j=1,nset)
         print 5201,ifile,titlef,nset
         if(prt) write(7,5201) ifile,titlef,nset
         if(discprt) write(8,5201) ifile,titlef,nset
	endif
	do 519 j=1,nset
      print 10,j,titled(j)(1:60)
      if(prt) write(7,10) j,titled(j)(1:60)
      if(discprt) write(8,10) j,titled(j)(1:60)
10	format('  data set #',i3,2x,a60)
519	continue
	pause
c      ans=getch(ktype)
518	continue	!end of loop to print titles
c
	print 138
138	format(
     & ' (1) Use this disk file',/,
     & ' (2) Try another disk',/,
     & ' Option number [1] = ')
	i=1
	call INPUTi(i)
c	read 3,i
c	if(i.le.0) i=1
	if(i.eq.2) then
	   close(unit=12)
	   fopen=.false.
	   present=.false.		!so asks for ndev
	   goto 31
	endif
338	continue
c Print file name/path
      if(pon()) write(7,20) infil
      if(discprt) write(8,20) infil
20	format(' Data from: ',a33)
c
	if(nfile.eq.1) then
	   ifile=1
	else if(readini) then
	   if(ifile1.gt.nfile) ifile1=1
	   print 141,nfile,ifile1
141	   FORMAT(' Read file number (1 -',i3,') [',i3,'] = ')
	   ifile=ifile1
	   call INPUTi(ifile)
	else
	   print 140,nfile
140	   FORMAT(' Read file number (1 -',i3,') = ')
	   call INPUTi(ifile)
	endif
c
	ifile1=ifile	!now redefine ifile1 as file actually read
c
	if(ifile.lt.1.or.ifile.gt.nfile) goto 338
	if(prt) write(7,1401)ifile
      if(discprt) write(8,1401)ifile
1401	format(' Data from disk file ',i5)
	nrect=lstrec(ifile)-jstrec(ifile)+1		!total number of records
	nrecd=nrect-1				!number of records for Xval etc
	irec=jstrec(ifile)
ccc	if(iver.ge.1003) irec=jstrec2(ifile)
	if(iver.ge.1003) then	!check array size
	   read(12,rec=irec) nset
	   if(nset.gt.njset-5) then
		call BELL(2)
		print 19,nset,njset
c19	     format(' Number of sets = ',i4,' too large for array size',
c     &	i4,/,' Allocate a larger size before proceeding')
		iesc=98
		CLOSE(unit=12)
		RETURN
	   endif
	endif
	if(iver.eq.1001) then
	   read(12,rec=irec) titled,nset,(nj(k),k=1,10),
     &	(titlep(k),k=1,20),setvar,
     &	(setx(k),k=1,10),iw,titlex,titley,ilabel,hdisp
	else if(iver.eq.1002) then
	   read(12,rec=irec) titlef,titleds,nset,(nj(k),k=1,10),
     &	(titlep(k),k=1,20),setvar,
     &	(setx(k),k=1,10),iw,titlex,titley,ilabel,hdisp
	    call SWAPTIT(titled,titleds,nset,njset,1) !copy titleds into titled
	else if(iver.eq.1003) then
	   read(12,rec=irec) nset,titlef,(titled(j),j=1,nset),
     & 	(nj(j),j=1,nset),(setx(j),j=1,nset),setvar,iw,
     &       titlex,titley,ilabel,hdisp,ioffset,k
c NB add ioffset=byte offset to where data starts, and
c k=number of values in datcop
	   ntot=0
	   do j=1,nset
		ntot=ntot+nj(j)
	   enddo
	endif
	do 521 j=1,nset
      print 10,j,titled(j)
      if(prt) write(7,10) j,titled(j)
      if(discprt) write(8,10) j,titled(j)
521	continue
c Now read 'compressed' data in DATCOP from disk
	if(iver.le.1002) then
	   ALLOCATE(datcop(3000))
	   n=1
	   do j=1,nrecd
		m=n+255
		irec=irec+1
		read(12,rec=irec) (DATCOP(i),i=n,m)
		n=n+256
	   enddo
	else
	   ALLOCATE(datcop(3*ntot))
	   read(12,rec=irec+ioffset) (DATCOP(i),i=1,k)
	endif
c now decode DATCOP into data arrays
	k=0		!index for datcop
	do 47 j=1,nset
	   n=nj(j)
	   do 471 i=1,n
	   k=k+1
471	   Xobs(i,j)=datcop(k)
	   if(hdisp) goto 47	!X only for histograms
	   do 472 i=1,n
	   k=k+1
472	   Yobs(i,j)=datcop(k)
	   do 473 i=1,n
	   k=k+1
473	   w(i,j)=datcop(k)
47	continue
	DEALLOCATE(datcop)
	RETURN
c
c
c Section to WRITE DATA into file number to be specified here (or as
c input if ifile<0
2	continue
c get disc to write to
c NB if file already open then can use as before ONLY if the open file
c is already an iver=1003 file.  Otherwise start a new CVDAT.DAT file, or
c write the data to an existing file that has iver=1003.
c=	if(fopen) goto 501
	if(readat) iver1=iver !version number for file from which data was read
	iver=1003			!version number being written now
	if(iver1.eq.1003.and.fopen) then
	   goto 501
	else if(fopen) then
	   call BELL(3)
	   print 111,iver1-1000,iver-1000,iver-1000
111	   format(
     &' Data was read from a version ',i2,
     &' data file, and this program is version',i2,'.',/,
     &' EITHER write to an existing V',i2,
     &' file, OR start a new CVDAT.DAT file')
	   CLOSE(unit=12)
	endif
c
3	continue
	if(ndisc.gt.1) print 35
35	format(' Specify disk for CVDAT data file (CVDAT.DAT):')
	call DISCNUM1(id,ndev,-1,1)	!ask for hard or floppy disk
	if(id.eq.-1) then
	   iesc=99
	   RETURN			!abort if ESC hit
	endif
	if(cluster) then
		infil='o:\cvfit\*.dat'
	else
	infil=ndev//'\CVDAT\CVDAT.DAT'
	endif
3031	call TITENT0(
     &'Name and path for output CVFIT data file:',infil,33,.false.)
	ifile=iabs(ifile)		!restore sign if nec
	if(cluster) then
	  if(infil(1:9).ne.'o:\cvfit\') then
	     print*,'TRY AGAIN: something like  o:\cvfit\yourname.dat '
           goto 3031
	  endif
	  if(infil(10:11).eq.'*.') then
	     print*,'TRY AGAIN: something like  o:\cvfit\yourname.dat '
           goto 3031
	  endif
	endif
c
	INQUIRE(file=infil,exist=present,flen=len)
	if(.not.present.or.len.eq.0) then
	   call BELL(1)
	   print 241,infil
241	   format(1x,a33,' does not exist: create it [Y] ? ')
	   read 101,ans
	   if(UC(ans).eq.'N') then
		present=.false.		!so asks for ndev
		goto 3    !get new path/filename
	   endif
c	open new data file (code as in CDSKINIT which should not be needed now)
	   do i=1,100
		jstrec(i)=0
		lstrec(i)=0
	   enddo
	   nfile=0
         OPEN(unit=12,file=INFIL,status='UNKNOWN',access='TRANSPARENT')
	   write(12,rec=1) nfile,jstrec,lstrec,iver
	   fopen=.true.
	else
         OPEN(unit=12,file=INFIL,status='UNKNOWN',access='TRANSPARENT')
c NB recl is in bytes (4*value for PDP)
	   read(12,rec=1) nfile,jstrec,lstrec,iver1
	   if(iver1.lt.1003) then
		call BELL(3)
		print 11,iver1-1000,iver-1000,iver-1000
11		format(
     &' This is a version ',i2,
     & ' data file, and this program is version',i2,'.',/,
     &' EITHER write to a V',i2,' file, OR start a new file')
		CLOSE(unit=12)
		goto 3    !get new path/filename
	   endif
	   fopen=.true.
	endif
501	continue
c
c Compress all real*4 data into DATCOP(): Xobs,Yobs,(max=3000 reals)
	ntot=0
	do j=1,nset
	   ntot=ntot+nj(j)
	enddo
	ALLOCATE(datcop(3*ntot))
	k=0		!index for datcop
	do 37 j=1,nset
	   n=nj(j)
	   do 371 i=1,n
	   k=k+1
371	   datcop(k)=Xobs(i,j)
	   if(hdisp) goto 37	!X only for histograms
	   do 372 i=1,n
	   k=k+1
372	   datcop(k)=Yobs(i,j)
	   do 373 i=1,n
	   k=k+1
373	   datcop(k)=w(i,j)
37	continue
	nvalt=k*4			!in bytes
	kdat=k
	ioffset=2048		!in bytes
	nvalt=nvalt+ioffset	!1024 bytes for stuff before datcop
	read(12,rec=1) nfile,jstrec,lstrec
	nfsav=nfile
	if(nfile.eq.0) then
	   jstrec(1)=1025		!now transparent -byte addresses
	   istrec=jstrec(1)
	endif
c
	if(filspec) goto 405	!ifile specified on input
c
44	continue
	if(.not.readat) then
	   ifile=nfile+1		!default file # to write=next unused one
	   print 404,ifile,ifile
404	   format(
     & ' Next unused file number is # ',i3,/,
     & ' File number to write data in (-1 to skip write)  [',i3,'] = ')
	else
	   ifile=ifile1    !default file # to write=one that data was read from
	   print 4041,ifile1,nfile+1,ifile1
4041	   format(
     & ' Data is from file #',i3,'; next unused file is # ',i3,/,
     & ' File number to write data in (0 to skip write)  [',i3,'] = ')
	endif
	call INPUTi(ifile)
	i=ifile
	if(i.le.0) goto 45		!skip write
	if(i.ge.1.and.i.le.nfile+1) ifile=i		!if i=0, ifile=default value
	if(ifile.ge.1.and.ifile.le.nfile) then  !overwrite existing file (nfile unchanged)
	   if(ifile.eq.1) jstrec(1)=1025	!new start record
	   istrec=jstrec(ifile)		!start record for plot to be overwritten
	else if(ifile.eq.nfile+1) then			!add a new file
	   if(nfile.eq.0) then		!when new file created (ifile=1 here)
		jstrec(1)=1025
		istrec=jstrec(1)
	   else
	      istrec=lstrec(nfile)+1	!start record for next file
		if(iver.ge.1003)istrec=lstrec(nfile)+32	!leave 32 byte gap
	   endif
	   nfile=nfile+1			!update nfile
	else
	   goto 45
	endif
405	continue
	if(istrec.lt.1025) then
	   call BELL(2)
	   print 4051,istrec
4051	   format(' ERROR: first record for data = ',i5)
	   goto 44
	endif
c Check that data will fit on floppy: make max file size=711 records which is
c most that will fit on 3.5 inch floppy (1457664 bytes usable, and
c 1457664/2048=711.75)
	lastr=istrec+nvalt-1		!last byte # for current data
	if(lastr.gt.145760) then
	   call BELL(2)
	   print 43,istrec,lastr
43	    format(
     &   ' Not enough room in data file (1st and last record #= ',2i9)
	    nfile=nfsav
	   goto 44
	endif
c If there are files following this one check that they will not be
c overwritten
	if(ifile.lt.nfsav) then
	   if(lastr.ge.jstrec(ifile+1)) then
		call BELL(2)
		iopt=1
		print 50,ifile+1,nfsav+1
50		format(
     &	' This file will overwrite part of the next(#',i3,')',/,
     &	' (1) Move the later files up to make room for this one',/,
     &	' (2) Write this to the next unused file: number',i3,/,
     &	' (3) Do not write to disc at all',/,
     &	' Option number [1] = ')
		i=1
		call INPUTi(i)
		if(i.ge.1.and.i.le.3) iopt=i
		if(iopt.eq.3) then
		   iesc=99
		   RETURN			!abort if ESC hit
		else if(iopt.eq.2) then
		   istrec=lstrec(nfsav)+1	!start record for next file
		   if(iver.ge.1003)istrec=lstrec(nfsav)+32	!leave 32 byte gap
		   nfile=nfsav+1			!update nfile
		else if(iopt.eq.1) then
		   ALLOCATE(titled1(njset),
     &		 nj1(njset),setx1(njset))
		   idrec=lastr-jstrec(ifile+1)+1	!move data up this # of records
		   do 51 if=nfsav,ifile+1,-1		!move all files after #ifile
c read (NB do not overwrite current values!)
			irec=jstrec(if)
			if(iver.eq.1001) then
			   read(12,rec=irec) titled1,nset1,nj1,titlep1,
     &		   setvar1,setx1,iw1,titlex1,titley1,ilabel1,hdisp1
	            else if(iver.eq.1002) then
			 read(12,rec=irec) titlef1,titleds1,nset1,nj1,titlep1,
     &		   setvar1,setx1,iw1,titlex1,titley1,ilabel1,hdisp1
	            else if(iver.ge.1003) then
			   read(12,rec=irec)
     &		    nset1,titlef1,(titled1(j),j=1,nset1),
     &	  	    (nj1(j),j=1,nset1),(setx1(j),j=1,nset1),
     &		    setvar1,iw1,
     &       	    titlex1,titley1,ilabel1,hdisp1,ioffset1,k1
			endif
		      ntot1=0
			do j=1,nset1
			   ntot1=ntot1+nj1(j)
			enddo
			ALLOCATE(datcop1(3*ntot1))
			read(12,rec=irec+ioffset1) (DATCOP1(i),i=1,k1)
c increment 1st,last record #
			jstrec(if)=jstrec(if)+idrec
			lstrec(if)=lstrec(if)+idrec
c and write back in new pos
			irec=jstrec(if)
			if(iver.eq.1001) then
			  write(12,rec=irec) titled1,nset1,nj1,titlep1,
     &		  setvar1,setx1,iw1,titlex1,titley1,ilabel1,hdisp1
			else if(iver.eq.1002) then
			 write(12,rec=irec)titlef1,titleds1,nset1,nj1,titlep1,
     &		  setvar1,setx1,iw1,titlex1,titley1,ilabel1,hdisp1
			else if(iver.eq.1003) then
			   write(12,rec=irec)
     &		    nset1,titlef1,(titled1(j),j=1,nset1),
     &	  	    (nj1(j),j=1,nset1),(setx1(j),j=1,nset1),
     &		    setvar1,iw1,
     &       	    titlex1,titley1,ilabel1,hdisp1,ioffset1,k1
			endif
			write(12,rec=irec+ioffset1) (DATCOP1(i),i=1,k1)
			DEALLOCATE(datcop1)
51		   continue		!end loop for each file to be moved
		   DEALLOCATE(titled1,nj1,setx1)
		endif
	   endif
	endif
c Update jstrec,lstrec
	jstrec(ifile)=istrec	!set value for current plot
	lstrec(ifile)=lastr	!set value for current plot
c
c Now write to disk: update rec #1
c NB Can write only iver=1003!
	iver=1003		!set the value for output
	write(12,rec=1) nfile,jstrec,lstrec,iver
c
c Things to go in 1st record of each plot (nfile copied into rec 1 also
c so it can be retrieved). Write all data set titles etc-(=1004 bytes)
c everything except Xobs,Yobs,SETX
	irec=istrec
	k=kdat
	write(12,rec=irec) nset,titlef,(titled(j),j=1,nset),
     & 	(nj(j),j=1,nset),(setx(j),j=1,nset),setvar,iw,
     &       titlex,titley,ilabel,hdisp,ioffset,k
c
c Now write DATCOP to disk
	write(12,rec=irec+ioffset) (DATCOP(i),i=1,k)
	CLOSE(unit=12)
	print 46,ifile,istrec,lastr,infil
	if(slock()) write(7,46) ifile,istrec,lastr,infil
      if(discprt) write(8,46) ifile,istrec,lastr,infil
46	format
     & (' Data file #',i3,' stored in records ',i9,' to ',i9,' of ',a33)
	call flush(7)
45	continue
	RETURN
c
30	continue	!file open error
	call BELL(2)
	print 32
32	format(' Error in opening file.',/,
     & '  May need initialising by running CDSKINIT')
	present=.false.		!so asks for ndev
	goto 31
	end


	subroutine SWAPTIT(titled,titleds,nset,njset,icop)
c To copy titled into titleds (icop=0) or vice versa (icop=1)
c Titled=60 character set titles used in iver=1001
c Titleds=40 character set titles used in iver=1002
	character*60 titled(njset)
	character*40 titleds(10)
c
	if(icop.eq.0) then
	   do 1 i=1,nset
1	   titleds(i)=titled(i)(1:40)
	else
	   do 2 i=1,nset
	   titled(i)(1:40)=titleds(i)
2	   titled(i)(41:60)='                   '
	endif
	RETURN
	end
