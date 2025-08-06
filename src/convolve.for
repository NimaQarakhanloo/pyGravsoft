      program convolve
c $Id: convolve.for 181 2008-08-14 23:19:51Z tjansson $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      C O N V O L V E
c
c  this program convolves a file with track data and time (jd)
c  with a kernel (impulse response) given in a separate file.
c  maximal kernel value is assumed to represent t=0, and should
c  be approximately in the center of the file.
c  
c  Input:
c
c  datafile name  (format: trackno*1000+no,lat,lon,h,dat,time)
c  impulseresponse file  (format: equidistant values, 1 pr line)
c  outputfile
c  dtfilt, dtdat  (time spacing of filter coef and data, seconds)
c                 (if dtdat = 0 the data time will be read from the
c                  data file)
c
c  RF/KMS, Greenland aug 92
c  last updated sep 93
c  update to Malaysia 2002 ... change trackno to trackno*10000
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      character*36 dfile,ffile,ofile
      dimension f(5000),d(5000),qf(5000),rf(5000),
     .is(5000),rlat(5000),rlon(5000),rh(5000),t(5000)
      dimension ff(-2500:2500)
      logical lstop
c
      write(*,*) 'Input: datafile name'
      write(*,*) '       kernelfile name'
      write(*,*) '       outputfile'
      write(*,*) '       kernel and data spacing (seconds)' 
c
      read(*,'(a36)') dfile
      read(*,'(a36)') ffile
      read(*,'(a36)') ofile
      read(*,*) dt,dtdat
c
      open(10,file=dfile,form='formatted',status='old')
      open(20,file=ffile,form='formatted',status='old')
      open(30,file=ofile,form='formatted',status='unknown')
c
c  read data kernel
c
      nf = 1
      nfmax = 0
      fsum = 0
      fmax = -99999.9
c
10    read(20,*,end=20) f(nf)
      if (f(nf).gt.fmax) then
        fmax = f(nf)
        nfmax = nf
      endif
      fsum = fsum + f(nf)
      nf = nf+1
      if (nf.gt.5000) stop 'too many kernel points'
      goto 10
20    nf = nf-1
      nc = nf/2
      do 30 i = 1,nf
30    f(i) = f(i)/fsum
c
      call initsp(f,nf,rf,qf)
      rimax = nfmax
      do 32 i = -10,10
        ri = nfmax + i/10.d0
        rr = spline(ri,f,nf,rf)
        if (rr.gt.fmax) then
          rimax = ri
          fmax = rr
        endif
32    continue
      dt2 = (min(rimax-1,nf-rimax)-0.5)*dt
c
      write(*,*) '--- C O N V O L V E ---'
      write(*,31) nf,rimax,fsum,dt2
31    format(' Number of points in filter file: ',i5/
     .' Maximum and filter zero occurs at index: ',f6.1/
     .' Filter sum: ',f9.4,', filter normalized to sum 1'/
     .' Filter window half-length in seconds: ',f7.1/
     .' Every 20th filter values: ')
      do 40 i = 1,nf
40    if (mod(i-nc,20).eq.0) write (*,41) (i-nc)*dt,f(i)
41    format(2f11.4)
c
c  read in data file 
c  loop over lines
c  ---------------
c
      iline = -999999
      ip = 1
      lstop = .false.
c
50    if (dtdat.gt.0) then
        read(10,*,end=51) 
     .  is(ip),rlat(ip),rlon(ip),rh(ip),d(ip)
        t(ip) = dtdat/86400.d0*(ip-1)
      else
        read(10,*,end=51) 
     .  is(ip),rlat(ip),rlon(ip),rh(ip),d(ip),t(ip)
      endif
c  rf malaysia
      ii = is(ip)/10000
      if (iline.eq.-999999) iline = ii
      if (ii.ne.iline) then
        np = ip-1
        goto 60 
      else
        ip = ip+1
        if (ip.gt.5000) stop 'too many points in line'
        goto 50
      endif
51    lstop = .true.
      np = ip-1
c 
c  convolve data track - 1 ... np
c
60    write(*,*) 
      write(*,*) 'Line: ',iline,', No of points: ',np
      if (np.lt.2) stop '*** too few points in line'
      write(*,611) t(1), t(np)
611   format(' times jd: ',2f12.6)
      ddt = (t(np)-t(1))/(np-1)*86400
      radeg = 180.d0/3.14159265
      dist = acos(sin(rlat(1)/radeg)*sin(rlat(np)/radeg) +
     .cos(rlat(1)/radeg)*cos(rlat(np)/radeg)*
     .cos((rlon(1)-rlon(np))/radeg))*radeg*111100
      v = dist/(ddt*np)
      write(*,61) 
     .rlat(1),rlon(1),rh(1),rlat(np),rlon(np),rh(np),ddt,v
61    format(' Start: ',2f11.5,f7.0/' Stop:  ',2f11.5,f7.0,/    
     .' Average data time spacing: ',f9.1,' sec'/
     .' Average velocity (m/s): ',f9.2)
c
c  scale and normalize convolution kernel to current v 
c
      nff = dt2/ddt
      sum = 0
      do 65 i = -nff,nff
        ff(i) = spline(rimax+i*ddt/dt,f,nf,rf)
        sum = sum+ff(i)
65    continue
      do 66 i = -nff,nff
66    ff(i) = ff(i)/sum
c      write(*,*) 'Space domain filter samples:'
c      do 67 i = -10,10
c67    write(*,68) i,ff(i)
c68    format(i5,f11.4)
c
c  do convolution
c
      dsum = 0
      dsum2 = 0
      ssum = 0
      ssum2 = 0
      do 70 i = 1, np
        sum = 0
        ffsum = 0
        do 71 j = -nff,nff
          k = i+j
          if (k.lt.1.or.k.gt.np) goto 71
          sum = sum + d(k)*ff(j)
          ffsum = ffsum+ff(j)
71      continue
        sum = sum/ffsum
c
        write(30,72) is(i),rlat(i),rlon(i),rh(i),sum,t(i)
72      format(i8,2f11.5,2f11.2,f11.5)
        dsum = dsum + d(i)
        dsum2 = dsum2 + d(i)**2
        ssum = ssum + sum
        ssum2 = ssum2 + sum**2
70    continue
c
      write(*,73) dsum/np,sqrt((dsum2-dsum**2/np)/(np-1)),
     .ssum/np,sqrt((ssum2-ssum**2/np)/(np-1))
73    format(' Mean and standard deviations:'/
     .' Before filtering: ',2f10.2/
     .' After      -    : ',2f10.2)      
c
c  end line loop
c  -------------
c  restore data for new line and go back
c
      if (lstop) goto 200
      iline = ii
      is(1) = is(ip)
      rlat(1) = rlat(ip)
      rlon(1) = rlon(ip)
      rh(1) = rh(ip)
      d(1) = d(ip)
      t(1) = t(ip)
      if (dtdat.gt.0) t(1) = 0
      ip = 2
      goto 50
c
200   continue
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  c
c                      i n i t s p                                 c
c                                                                  c
c  initialization procedure for fast 1-dimensional equidistant     c
c  spline interpolation, with free boundary end conditions         c
c  reference: josef stoer: einfuhrung in die numerische mathematik c
c  i, springer 1972.                                               c
c                                                                  c
c  parameters (real):                                              c
c                                                                  c
c  y  given values, y(1), ..., y(n)                                c
c                                                                  c
c  r  spline moments (1 ... n), to be used by function 'spline'    c
c                                                                  c
c  q  work-array, declared at least 1:n                            c
c                                                                  c
c  rene forsberg, july 1983                                        c
c                                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine initsp(y, n, r, q)
c
      implicit double precision(a-h,o-z)
      dimension y(19), r(19), q(19)
c
      q(1) = 0.0
      r(1) = 0.0
      do 11 k = 2, n-1
        p = q(k-1)/2+2
        q(k) = -0.5/p
        r(k) = (3*(y(k+1)-2*y(k)+y(k-1)) - r(k-1)/2)/p
   11 continue
      r(n) = 0.0
      do 12 k = n-1, 2, -1
        r(k) = q(k)*r(k+1)+r(k)
   12 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  c
c                          s p l i n e                             c
c                                                                  c
c  fast one-dimensional equidistant spline interpolation function. c
c                                                                  c
c  parameters:                                                     c
c                                                                  c
c  x   interpolation argument (real), x = 1 first data-point,      c
c      x = n last data-point. outside the range linear extra-      c
c      polation is used.                                           c
c                                                                  c
c  y   real*8 array, 1 .. n : data values                          c
c                                                                  c
c  r   do: spline moments calculated by subroutine 'initsp'        c
c                                                                  c
c  programmer:                                                     c
c  rene forsberg, june 1983                                        c
c                                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function spline(x, y, n, r)
c
      implicit double precision (a-h, o-z)
      dimension y(19), r(19)
c
      if(x.ge.1.0) go to 1
        spline = y(1) + (x-1)*(y(2)-y(1)-r(2)/6)
      return
    1 if(x.le.float(n)) go to 2
        spline = y(n) + (x-n)*(y(n)-y(n-1)+r(n-1)/6)
      return
    2   j = ifrac(x)
        xx = x - j
        spline = y(j) +
     .           xx * ((y(j+1)-y(j)-r(j)/3-r(j+1)/6) +
     .           xx * (r(j)/2 +
     .           xx * (r(j+1)-r(j))/6))
      return
      end
c
      integer function ifrac(r)
      implicit double precision(a-h,o-z)
      if (r.lt.0) goto 10
      ifrac = r
      return
10    i = r
      if (i.eq.r) goto 20
      ifrac = i-1
      return
20    ifrac = i
      return
      end
