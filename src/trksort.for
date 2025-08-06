      program trksort
      implicit double precision(a-h,o-z) 
c $Id: trksort.for 251 2008-12-12 10:26:49Z cct $
c
c  program sorts file with 80-char gravity data where stat no is a line
c  no but order of points within line is random
c  rf jan 1992
c
      character c(500)*80
      character*6 no,oldno
      character*1 dir
      dimension rfi(500),rla(500),ix(500),d(500)
      logical lmore,lfirst
c
      open(10,file='pub29',status='old')
      open(20,file='ofile',status='unknown')
c
      lmore = .true.
      lfirst = .true.
      rfimin = 999
      rfimax =-999
      rlamin = 999
      rlamax =-999
      n = 1
      read(10,'(80a)',end=19) c(n)
c
10    no = c(n)(57:62)
      if (no.ne.oldno.and.(.not.lfirst)) goto 20
      oldno = no
      lfirst = .false. 
      read(c(n)(3:4),'(i2)') ilat
      read(c(n)(5:9),'(f5.2)') rlat
      read(c(n)(12:14),'(i3)') ilon
      read(c(n)(15:19),'(f5.2)') rlon
      rlat = ilat+rlat/60
      rlon = ilon+rlon/60
      if (rlat.lt.rfimin) rfimin = rlat
      if (rlat.gt.rfimax) rfimax = rlat
      if (rlon.lt.rlamin) rlamin = rlon
      if (rlon.gt.rlamax) rlamax = rlon
      rfi(n) = rlat
      rla(n) = rlon
      n = n+1
      if (n.gt.500) stop 'n too large'
      read(10,'(80a)',end=19) c(n)
      goto 10
c
c  sort line
c
19    lmore = .false.
20    n = n-1
      cosfi = cos(rfimin/180*3.1415926)
21    format(' line ',a6,i5,4f10.4,1x,a1)
      if (rfimax-rfimin.gt.(rlamax-rlamin)*cosfi) then
        do 30 k = 1, n
30      d(k) = rfi(k)
        dir = 'N'
      else
        do 31 k = 1, n
31      d(k) = rla(k)
        dir = 'E'
      endif
      write(*,21) oldno,n,rfimin,rfimax,rlamin,rlamax,dir
      call sortin(d,n,ix)
33    format(7i4)
      write(20,'(a80)') (c(ix(j)),j=1,n)
c  
      c(1) = c(n+1)
      n = 1
      lfirst = .true.
      rfimin = 999
      rfimax =-999
      rlamin = 999
      rlamax =-999
      if (lmore) goto 10
      end
c
      subroutine sortin(d,n,ix)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                    S O R T I N
c
c  Subroutine sorts data in array d from 1 to n.
c  Result is given in integer array i, so that d(ix(k)), k = 1..n
c  is the sorted sequence.
c  The routine is fast for near-sorted data, slow otherwise.
c
c  (c) RF/KMS. Jan 1992
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8 (a-h,o-z)
      dimension d(*),ix(*)
c
      ix(1) = 1
      do 20 i = 2, n
        dd = d(i)
        j = i-1
10      if (dd.lt.d(ix(j))) then
          ix(j+1) = ix(j)
        else
          ix(j+1) = i
          goto 20
        endif
        j = j-1
        if (j.gt.0) goto 10
        ix(1) = i
20    continue
      return
      end
