      program ntozeta    
c $Id: n2zeta.for 273 2017-12-05 14:49:55Z cct $
c the program computes the difference between the orthometric
c height and the normal height as being equal to the Bouger-anomaly
c multiplied by the altitude and divided by normal gravity.
c updated 2017-12-05 to consider complete Bouger-anomaly
c a file with geoid heights may be converted to heightanomalies
c (or the reverse) using the program.
c For the conversion a fixed density of 2.67 g/cm**3 is used,
c corresponding to a Bouger factor of 0.1119 mgal/m.
c
c input files:
c  file with free-air gravity anomalies of the points to be corrected.
c  if corrections are to be applied, file with geoid heights or height
c  anomalies (m)
c output file
c  file to hold corrections
c  
c programmed july 2002 by C.C.Tscherning. Update 2009-02-12.
c Update 2017-12-05 by Sabah
c
      implicit none
      integer i,j,k,n,m
      real*8 rlat,rlatr,rlon,h,rlon1,rlat1,obs1(19),
     *obs(10),d0,dg,ba,bc,corr,corrc,scor,scorc,sscor,sscorc,sd,ssd,so,sso
      real gnormal,pi,degrad
      character*72 ifile,ofile,dfile
      logical ladd,lg2h
c
      d0=0.0d0
      scor=d0
	  scorc=d0
      sscor=d0
	  sscorc=d0
      sd=d0
      ssd=d0
      so=d0
      sso=d0
      pi=atan(1.0d0)*4.0d0
      degrad=pi/180.0d0
c
      write(*,*)' n2zeta, ver. 2017-12-05 '
      write(*,*)
     *' conversion geoid to height anomalies or the reverse ? '
      read(*,*)lg2h
      if (lg2h) then
       write(*,*)' conversion from geoid to height anomalies '
      else
       write(*,*)' conversion from height anomalies to geoid '
      end if
c
      write(*,*)' input gravity and output file names (2 lines) '
      read(*,'(a)')ifile
      read(*,'(a)')ofile
      write(*,*)' files = ',ifile, ofile
      open(12,file=ifile)
      open(14,file=ofile)
      write(*,*)' input data element number (max 10) '
      read(*,*)n
      write(*,*)' correction to be applied on file ? (t/f) '
      read(*,*)ladd
      if (ladd) then
       write(*,*)' input data file name '
       read(*,'(a)')dfile
       open(15,file=dfile)
       write(*,*)dfile
       write(*,*)' input data element number (max 10) '
       read(*,*)m
      end if
c
      if (ladd) then
       write(*,*)' free-air anomalies = fa '
       write(*,*)' Bouger anomalies = ba '
       if (lg2h) then
        write(*,108)
  108   format(' no    lat        long      H         fa       ba   ',
     *  '   corr      N     zeta ') 
       else
        write(*,109)
  109   format(' no    lat        long      H         fa       ba   ',
     *  '   corr    zeta    N   ') 
       end if
      else
      write(*,106)
  106 format(' no    lat        long      H         fa       ba   ',
     *'     corr   ') 
      end if
c 
      j=0
c data on GRAVSOFT format, number, latitude, height, gravity
c anomaly in mgal.
 10   read(12,*,end=99)i,rlat,rlon,h,(obs(k),k=1,n)
      rlatr=rlatr*degrad
      gnormal=9.780327d0*(1.0d0+0.005279*sin(rlatr)**2)-3.088d-6*h
      gnormal=gnormal*1.0d5
      if (ladd) then
       read(15,*,end=99)i,rlat1,rlon1,h,(obs1(k),k=1,m)
       if (abs(rlat-rlat1).gt.0.001.or.abs(rlon-rlon1).gt.0001)
     * write(*,*)' WARNING  out of sequence '
       so=so+obs1(m)
       sso=sso+obs1(m)**2
      end if
      dg=obs(n)
      ba=dg-h*0.1119
	  bc=obs(2)
      if (lg2h) then
       corr=ba*h/gnormal  
	   corrc=bc*h/gnormal
      else
       corr=-ba*h/gnormal  
	   corrc=-bc*h/gnormal
      end if
      scor=scor+corr
      sscor=sscor+corr**2
	  scorc=scorc+corrc
      sscorc=sscorc+corrc**2	  
      if (ladd) then
      write(*,100)i,rlat,rlon,h,dg,ba,corr,obs1(m),obs1(m)-corr
       write(14,100)i,rlat,rlon,h,obs1(m)-corr,corr
       sd=sd+obs1(m)-corr
       ssd=ssd+(obs1(m)-corr)**2
      else
      write(*,100)i,rlat,rlon,h,dg,ba,corr,corrc
      write(14,100)i,rlat,rlon,h,ba,corr,corrc
      end if
 100  format(i6,f9.4,f10.4,f8.2,6f9.3)
      j=j+1
      go to 10
  99  write(*,*)' number of data ',j
      write(*,*)'                         ba: mean      std.  bc: mean      std.'
      sscor=sqrt((sscor-scor**2/j)/(j-1))
	  sscorc=sqrt((sscorc-scorc**2/j)/(j-1))
      write(*,111)scor/j, sscor, scorc/j, sscorc 
  111 format(' corrections             ',2f10.3)
      if (ladd) then
       ssd=sqrt((ssd-sd**2/j)/(j-1))
       sso=sqrt((sso-so**2/j)/(j-1))
       if (lg2h) then
        write(*,112)sd/j, ssd,so/j,sso  
  112   format(' height anomalies (zeta) ',2f10.3,/
     *  ' geoid heights (N)       ',2f10.3)
       else
        write(*,113)sd/j, ssd,so/j,sso  
  113   format(' geoid heights (N)       ',2f10.3,/
     *  ' heigt anomalies (zeta)  ',2f10.3)
       end if
      end if
 
      close(12)
      close(14)
      stop
      end



     
