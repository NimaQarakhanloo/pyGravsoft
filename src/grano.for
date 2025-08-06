      program grano
c $Id: grano.for 181 2008-08-14 23:19:51Z tjansson $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                       G R A N O
c
c  program for computing GRS80 free-air and Bouguer anomalies from a
c  list of gravity values (e.g., from gradj) and coordinates.
c  The lists need not be complete.
c  
c  Input:
c
c  <gravityfile>
c  <coordinate file>
c  <elev file>
c  <outputfile>
c  ifmt, lm, density 
c
c  ifmt ...     format of coordinate list: 1 deg, 2 deg min, 3 dms
c
c  lm ...       output format (t = deg min, f = deg)
c
c  density  ... ref.density. For oceans (h < -2m) rho-1.03 is used.
c               if density=0 no Bouguer is output (for ice caps)
c
c  (c) RF/KMS apr 94. Updated oct 94.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8 (a-h,o-z)
      character*72 gfile,cfile,efile,ofile
      parameter (igmax=2000,iemax=2000)
      dimension igstat(igmax),g(igmax)
      dimension iestat(iemax),elev(iemax)
      logical lg(igmax), lm, lelev
      character*12 ch
      radeg = 180/3.14159265d0
c
      write(*,*) 'Input: gravityfile'
      write(*,*) '       coordinatefile'
      write(*,*) '       elevfile (0 if elev in coordinatefile)'
      write(*,*) '       outputfile'
      write(*,*) '       ifmt (1:deg,2:dm,3:dms), lout(t:dm,f:deg),'
      write(*,*) '       density (g/cm**3) '
c
      read(*,'(a)') gfile
      read(*,'(a)') cfile
      read(*,'(a)') efile
      read(*,'(a)') ofile
      read(*,*) ifmt,lm,rho
c
      open(10,file=gfile,status='old')     
      open(20,file=cfile,status='old')
      open(30,file=ofile,status='unknown')
c
      write(*,*)
      write(*,*) '--- G R A N O ---'
      lelev = (efile.ne.'0')
c
c  scan through gravity file
c  skip possible first line with =
c
      ng = 1
10    read(10,'(a12)',end=20) ch 
      if (ch.eq.'            ') goto 10
      if (ch(1:1).eq.'#') goto 10
      if ((ch(1:1).eq.'='.or.ch(2:2).eq.'=').and.ng.eq.1) goto 10
      if (ch(1:1).eq.'='.or.ch(2:2).eq.'=') goto 20
      backspace(10)
      read(10,*,end=20) igstat(ng),g(ng)
      lg(ng) = (.false.)
      ng = ng+1
      if (ng.gt.igmax) stop '*** too many gravity values, > igmax'
      goto 10
c
20    ng = ng-1
      write(*,*) 'Number of gravity values input: ',ng
c
c  read elevation list 
c
24    if (lelev) then
        ne = 1
25      open(15,file=efile,form='formatted',status='old')
        read(15,*,end=26) iestat(ne),elev(ne) 
        ne = ne+1
        if (ne.gt.iemax) stop '*** too many elevations, > iemax'
        goto 25
26      ne = ne-1
        write(*,*) 'Number of elevations input: ',ne
        close(15)
      endif
c
c  coordinate list reading loop
c
      na = 0
      nc = 0
      famin = 9999.d9
      famax = -9999.d9
      bamin = 9999.d9
      bamax = -9999.d9
c
30    goto (31,32,33),ifmt   
c
31    if (lelev) then
        read(20,*,end=40) istat,rlat,rlon
      else
        read(20,*,end=40) istat,rlat,rlon,h     
      endif
      goto 35
c
32    if (lelev) then
        read(20,*,end=40) istat,ilatd,rlat,ilond,rlon
      else
        read(20,*,end=40) istat,ilatd,rlat,ilond,rlon,h
      endif
      if (ilatd.lt.0) rlat = -rlat
      if (ilond.lt.0) rlon = -rlon
      rlat = ilatd + rlat/60
      rlon = ilond + rlon/60
      goto 35
c
33    if (lelev) then
        read(20,*,end=40) istat,ilatd,ilatm,rlat,ilond,ilonm,rlon
      else
        read(20,*,end=40) istat,ilatd,ilatm,rlat,ilond,ilonm,rlon,h
      endif
      if (ilatd.lt.0) ilatm = -ilatm
      if (ilond.lt.0) ilonm = -ilonm
      if (ilatm.lt.0) rlat = -rlat
      if (ilonm.lt.0) rlon = -rlon
      rlat = ilatd + ilatm/60.d0 + rlat/3600
      rlon = ilond + ilonm/60.d0 + rlon/3600
c
c  check if elevation exist
c
35    if (lelev) then
        i = listno(iestat,ne,istat)
        if (i.eq.0) then
          write(*,*) '- no elevation for point: ',istat
          goto 30
        endif
        h = elev(i)
      endif
c
c  find gravity value and compute anomalies
c
      nc = nc+1
      ig = listno(igstat,ng,istat)
      if (ig.eq.0) then
        write(*,*) '- no gravity for point: ',istat
      else
        na = na+1
        lg(ig) = (.true.)
c
        sinfi = sin(rlat/radeg)
	gamma = 978032.677d0*(1+.00193185135d0*sinfi**2)/
     .  sqrt(1 - .00669438002d0*sinfi**2)
        if (h.ge.-2) then
          fa = g(ig)-gamma+0.30877d0*(1-.00142*sinfi**2)*h -.75d-7*h**2
          ba = fa - 0.1119d0*rho/2.67*h            
        else
          fa = g(ig)-gamma
          ba = fa - 0.1119d0*(rho-1.03)/2.67*h
        endif       
c
        if (lm) then
          ilat = rlat
          rlatm = abs(rlat-ilat)*60
          if (ilat.eq.0.and.rlat.lt.0) rlatm = -rlatm
          ilon = rlon
          rlonm = abs(rlon-ilon)*60
          if (ilon.eq.0.and.rlon.lt.0) rlonm = -rlonm         
        endif          
        if (rho.ne.0) then
          if (lm) then
            write(30,37) istat,ilat,rlatm,ilon,rlonm,h,g(ig),fa,ba
          else
            write(30,38) istat,rlat,rlon,h,g(ig),fa,ba
          endif
        else
          if (lm) then
            write(30,37) istat,ilat,rlatm,ilon,rlonm,h,g(ig),fa
          else
            write(30,38) istat,rlat,rlon,h,g(ig),fa
          endif
        endif
37      format(i7,2(i5,f7.3),f9.2,f11.2,2f9.2)
38      format(i7,2f11.5,f9.2,f11.2,2f9.2)
        if (fa.lt.famin) famin=fa
        if (fa.gt.famax) famax=fa
        if (ba.lt.bamin) bamin=ba
        if (ba.gt.bamax) bamax=ba
      endif
      goto 30
c
40    write(*,41) nc,na,famin,famax
41    format(' Number of full coordinates:',i6,', anomalies output:',i6,
     ./' min and max free-air anomaly:',2f11.2)
      if (rho.gt.0) write(*,42) bamin,bamax
42    format('             Bouguer anomaly: ',2f11.2)
c
c  information on gravity points without coordinates
c
      write(*,*)
      write(*,*) 'Stations in gravity file without coordinates: '
      lm = .false.
      do 50 i = 1,ng
        if (.not.lg(i)) then
          write(*,'(i8)') igstat(i)
          lm = .true.
        endif
50    continue
      if (.not.lm) write(*,*) '- none -'
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                    L I S T N O
c
c  subroutine for finding station number in unsorted station list
c  'idx' is zero if 'istat' is not in array 'ia' of 'n' elements
c  search begins at index from previous call
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function listno(ia, n, istat)
      dimension ia(*)
      save idx
      if (n.eq.0) goto 30
      if (idx.le.0.or.idx.gt.n) idx = 1
      j = idx
10    if (ia(j).eq.istat) goto 20
      j = j+1
      if (j.gt.n) j = 1
      if (j.eq.idx) goto 30
      goto 10
20    idx = j
      listno = idx
      return
30    idx = 0
      listno = idx
      return
      end
