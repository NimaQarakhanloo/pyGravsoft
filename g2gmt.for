      program g2gmt
c $Id: g2gmt.for 276 2009-02-19 11:03:09Z cct $
c
c  program to convert gravsoft grid to simple list in gmt format (lon,lat,data)
c  
c  rf, 3-2-2005
c
      implicit real*8 (a-h,o-z)
      dimension c(50000)
      character*256 ifile,ofile
      logical lutm
c
      write(*,*) 
     .'input name of inputfile (gravsoft) and outputfile:'
      write(*,*)
      read(*,'(a)') ifile
      read(*,'(a)') ofile
      open(10,file=ifile)
      open(20,file=ofile,status='unknown')
      rmin = 9.d9
      rmax = -9.d9
c
      read(10,*) rfi1,rfi2,rla1,rla2,dfi,dla
      lutm = .false.
      if (abs(rfi1).gt.100.and.abs(rfi2).gt.100) lutm = .true.
      if (lutm) then 
        read(10,*) iell,izone
        write(*,*) '- utm grid detected zone ',izone
      endif
      nn = (rfi2-rfi1)/dfi+1.5
      ne = (rla2-rla1)/dla+1.5
      write(*,*) nn,ne
      if (ne.gt.5000) stop 'ne too large'
      do 16 i = nn,1,-1
        read(10,*) (c(j),j=1,ne)
        do 15 j= 1,ne
	  rfi = rfi1+(i-1)*dfi
	  rla = rla1+(j-1)*dla
	  if (c(j).lt.rmin) rmin = c(j)
	  if (c(j).gt.rmax) rmax = c(j)
          if (lutm) then
	     stop 'utm not implemented'
          else
	      write(20,14) rla,rfi,c(j)    
14          format(' ',2f10.6,f9.2)
          endif
15      continue
16    continue
c 
      write(*,20) rmin,rmax
 20   format(' min,max =',2f10.2)      
      end
