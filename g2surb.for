      program g2sur
c $Id: g2surb.for 181 2008-08-14 23:19:51Z tjansson $
c
c  program to convert gravsoft grids to surfer ascii format (.grd)
c  
c  (c) Rene Forsberg KMS, Denmark
c  Rio sep 11, 1997
c  updated for UTM july 1999
c  updated for surfer binary format apr14, 2002, rf
c
      implicit real*8 (a-h,o-z)
      character*36 ifile,ofile
      real*4 c(8000),hdr(14)
      integer*2 ihdr(4)
      equivalence (hdr,ihdr)
      equivalence (hdr(3),rla1)
      equivalence (hdr(5),rla2)
      equivalence (hdr(7),rfi1)
      equivalence (hdr(9),rfi2)
      equivalence (hdr(11),zmin)
      equivalence (hdr(13),zmax)
c
      write(*,*) 
     .'input name of inputfile (gravsoft) and outputfile (.grd):'
      write(*,*)
      read(*,'(a)') ifile
      read(*,'(a)') ofile
      open(10,file=ifile)
      open(20,file=ofile,status='unknown',form='unformatted',
     .access='direct',recl=4)
c
      read(10,*) rfi1,rfi2,rla1,rla2,dfi,dla
      if (abs(rfi1).gt.100.and.abs(rfi2).gt.100) then
        read(10,*) iell,izone
        write(*,*) '- utm grid assumed -',iell,izone
        rfi1 = rfi1/1000
        rfi2 = rfi2/1000
        rla1 = rla1/1000
        rla2 = rla2/1000
        dfi = dfi/1000
        dla = dla/1000
      endif
      nn = (rfi2-rfi1)/dfi + 1.5
      ne = (rla2-rla1)/dla + 1.5
      if (ne.gt.8000) stop 'e-w rows too long'
      zmin = 9999.d9
      zmax = -9999.d9
c 
      write(*,20) nn,ne
20    format(' grid points n e:',2i5)
c
      ihdr(1) = ichar('S')*256+ichar('D')
      ihdr(2) = ichar('B')*256+ichar('B')
      ihdr(3) = ne
      ihdr(4) = nn
      do ir = 1,10
        write(20,rec=ir) hdr(ir)
      enddo
c
      do 30 i = nn,1,-1
        read(10,*) (c(j),j=1,ne)
        ir = (i-1)*ne + 15
        do 29 j = 1,ne
          if (c(j).lt.zmin) zmin = c(j)
          if (c(j).gt.zmax) zmax = c(j)
          write(20,rec=ir) c(j)
          ir = ir+1
29      continue
30    continue
c
      write(*,40) zmin,zmax 
40    format(' zmin, zmax = ',2f11.3)
      do ir = 11,14
        write(20,rec=ir) hdr(ir)
      enddo
      end
