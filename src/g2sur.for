      program g2sur
c $Id: g2sur.for 181 2008-08-14 23:19:51Z tjansson $
c
c  program to convert gravsoft grids to surfer ascii format (.grd)
c  
c  (c) Rene Forsberg KMS, Denmark
c  Rio sep 11, 1997
c  updated for UTM july 1999
c
      implicit real*8 (a-h,o-z)
      character*36 ifile,ofile
      parameter (nmax=4000000)
      real*4 c(nmax)
c
      write(*,*) 
     .'input name of inputfile (gravsoft) and outputfile (.grd):'
      write(*,*)
      read(*,'(a)') ifile
      read(*,'(a)') ofile
      open(10,file=ifile)
      open(20,file=ofile,status='unknown')
c
      read(10,*) rfi1,rfi2,rla1,rla2,dfi,dla
      if (rfi1.gt.100.and.rfi2.gt.100) then
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
      write(*,*) nn,ne
      k = 0
      zmin = 9999.d9
      zmax = -9999.d9
      do 10 i = 1, nn
        read(10,*) (c(j),j=k+1,k+ne)
        do 11 j=k+1,k+ne
          if (j.gt.nmax) stop '*** too small grid, increase nmax ***'
          if (c(j).ge.9999) then
            c(j) = 1.70141e38
          else
            if (c(j).lt.zmin) zmin = c(j)
            if (c(j).gt.zmax) zmax = c(j)
          endif
11      continue
        k = k+ne
10    continue
c 
      write(*,20) nn,ne,zmin,zmax
20    format(' grid points n e:',2i5,', zmin,zmax =',2f10.2)      
c
      write(20,21) ne,nn,rla1,rla2,rfi1,rfi2,zmin,zmax
21    format('DSAA'/' ',2i5/' ',2f11.5/' ',2f11.5/' ',2f11.2) 
      k = ne*(nn-1)
      do 30 i = 1, nn
        write(20,22) (c(j),j=k+1,k+ne)
22      format(/,50(/' ',6g12.5))
        k = k-ne
30    continue
      end


