      program gpfit
c $Id: gpfit.for 181 2008-08-14 23:19:51Z tjansson $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                            G P F I T
c
c  Program for fitting flat-earth logarithmic covariance function to
c  gravity data in a file.
c
c  input: ifile  (gravity data)
c         ds,smax,d1,d2,t1,t2  (all in km)
c
c  Reference: A new covariance model for inertial gravimetry and
c  gradiometry. Journal of Geophysical Research vol. 92 pp. 1305-10, 1987.
c
c  (c) Rene Forsberg, KMS-Denmark
c  sep 1998
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8 (a-h,o-z)
c
      parameter (maxobs=20000,maxcov=100)
      dimension rfi(maxobs),rla(maxobs),d(maxobs)
      dimension cov(maxcov),ncov(maxcov)
      character*48 ifile
c
      degkm = 111.11
      radeg = 180/3.14159265d0
c
      write(*,*) '--- G P F I T ---'
      write(*,*) 'input: gravityfile '
      read(*,'(a)') ifile
      write(*,*) 'input: ds,smax,d1,d2,t1,t2 (km) '
      read(*,*) ds,smax,d1,d2,t1,t2
      ne = smax/ds+1.5
      if (ne.gt.maxcov) stop 'too many intervals - increase maxcov'
c
      open(10,file=ifile)      
c
c  initialize statistics 
c
      np = 1
      dsum = 0
      dsum2 = 0
      do i = 1, ne
        cov(i) = 0
        ncov(i) = 0
      enddo
c
10    read(10,*,end=20) i,rfi(np),rla(np),r,d(np)
      dsum = dsum + d(np)
      dsum2 = dsum2 + d(np)**2
      np = np+1
      if (np.gt.maxobs) stop 'too many data - increase maxobs'
      goto 10
c
20    np = np-1
      if (np.lt.1) stop 'no data in file'
      c0 = dsum2/np
      write(*,21) np,dsum/np,c0,sqrt(c0)
21    format(' n =',i5,', mean =',f7.2,', C0 =',f9.2,', sqrtC0 =',f7.2)
c
      do 30 i = 1,np
        coskm = cos(rfi(i)/radeg)*degkm
        do 30 j = 1,i
          xx = (rla(i)-rla(j))*coskm
          yy = (rfi(i)-rfi(j))*degkm
          r = sqrt(xx**2+yy**2)
          ir = (r+ds/2)/ds+1
          if (ir.le.ne) then
            cov(ir) = cov(ir) + d(i)*d(j) 
            ncov(ir) = ncov(ir)+1
          endif
30    continue
c
      do i = 1,ne
        if (ncov(i).ne.0) cov(i) = cov(i)/ncov(i)
      enddo
c
c  fit to function - search integer km
c
      idmin = 99999
      itmin = 99999
      summin = 9.9d9
      write(*,35) nint(d1),nint(d2),nint(t1),nint(t2)
35    format(' - searching for best fit in interval:',i3,' -',i4,',', 
     .i3,' -',i4,' km')
c
      do 40 id = nint(d1),nint(d2) 
      do 40 it = nint(t1),nint(t2)
        sum = 0
        dd = id
        t = it
        do 42 k = 1,ne
          s = (k-1)*ds
          c = covp(3,3,s,0.d0,0.d0,dd,t,c0)
          if (ncov(k).eq.0) goto 42
          sum = sum + (cov(k)-c)**2
42      continue
        if (sum.lt.summin) then
          idmin = id
          itmin = it
          summin = sum
        endif
40    continue
      dd = idmin
      t = itmin
      write(*,45) idmin,itmin
45    format(' Best values of fit: D =',i4,' km, T =',i4,' km')
c
      write(*,51) 
51    format(/' Empirical and fitted covariances: '
     ./'    Dist (km)     No of products      Cov     Fitted cov')
      do 55 i = 1,ne
        s = (i-1)*ds
        c = covp(3,3,s,0.d0,0.d0,dd,t,c0)
        write(*,52) s,ncov(i),cov(i),c
52      format('    ',f6.1,'       ',i9,4x,2f12.2)
55    continue
      end
c
      real*8 function covp(ki1,ki2,x,y,zz,d1,d2,var)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                         C O V P 
c
c  function for computation of self-consistent covariance functions
c  for geoid undulations, deflections of the vertical and gravity
c  anomalies, using planar logaritmic covariance functions. 
c  these covariance functions features very good fits to empirical
c  covariance functions, the main decay of the power spectrum 
c  corresponding to kaula's rule.
c
c
c  parameters:
c
c  ki1, ki2    kind of quantities (1: geoid undulations,
c              (3: gravity, 6: N-S deflection, 7: E-W deflection)
c
c  x, y        coordinate differences x2-x1, y2-y1 in kilometer
c              (x positive to the east, y to the north)
c
c  zz          = h1 + h2, sum of heights of points in meter
c
c  d1          depth to top layer (km) - high freq. attenuation
c              (twice the'bjerhammar sphere depth')
c
c  d2          thickness to bottom layer (km) - low freq. attenuation
c
c  var         gravity anomaly variance in mgal
c
c  covariances are output in units  mgal**2, mgal*m or m**2.
c  deflections are given in arcsecs
c
c  programmer: rene forsberg, danish geodetic institute, feb 1986
c  updated nov 1994
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8 (a-h, o-z)
c
      d12 = d1+d2
      d13 = d1+2*d2
      d14 = d1+3*d2
      fc = var/log(d12**3/d13**3*d14/d1)
      zkm = zz/1000
      if (ki1.eq.3.and.ki2.eq.3) goto 20
c
      covp = fc*(plc(ki1,ki2,x,y,d1+zkm)
     *          -3*plc(ki1,ki2,x,y,d12+zkm)
     *          +3*plc(ki1,ki2,x,y,d13+zkm)
     *            -plc(ki1,ki2,x,y,d14+zkm))
      return
c
c  faster computation for gravity autocovariance
c
20    z1 = d1+zkm
      z2 = d12+zkm
      z3 = d13+zkm
      z4 = d14+zkm
      s = x**2 + y**2
      r1 = sqrt(s + z1**2)
      r2 = sqrt(s + z2**2)
      r3 = sqrt(s + z3**2)
      r4 = sqrt(s + z4**2)
      covp = fc*log((z2+r2)**3/(z3+r3)**3*(z4+r4)/(z1+r1))
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                        P L C
c
c  computes basic cross- and auto covariance components 
c  corresponding to  c(gg) = -log(z + r)
c  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function plc(ki1,ki2,x,y,z)
      implicit real*8 (a-h, o-z)
c
      r = sqrt(x**2 + y**2 + z**2)
c
c  table for various kinds 
c
      goto (11, 12, 12, 9, 9, 13, 14,
     *      12, 15, 15, 9, 9, 16, 17,
     *      12, 15, 15, 9, 9, 16, 17,
     *       9,  9,  9, 9, 9,  9,  9,
     *       9,  9,  9, 9, 9,  9,  9,
     *     131,161,161, 9, 9, 18, 19,
     *     141,171,171, 9, 9, 19, 20), (ki1-1)*7+ki2
9     plc = 0
      stop '*** illegal cov codes in pcov'
c
c  geoid * geoid
c
11    plc = .0000010404d0*(.75*z*r + (.25*r**2-.75*z**2)*log(z+r))
      return
c
c  geoid * gravity, ksi, eta (mgal*meter, mgal*arcsec)
c
12    plc = -.00102*(r - z*log(z+r))
      return
13    plc =  .0001052*y*(log(z+r) + .5 + z/(z+r))
      return
131   plc = -.0001052*y*(log(z+r) + .5 + z/(z+r))
      return
14    plc =  .0001052*x*(log(z+r) + .5 + z/(z+r))
      return
141   plc = -.0001052*x*(log(z+r) + .5 + z/(z+r))
      return
c
c  gravity covariance
c
15    plc = -log(z + r)
      return
c
c  gravity/deflection cross covariances
c
16    plc = -y/(z+r) /4.85
      return
161   plc =  y/(z+r) /4.85
      return
17    plc = -x/(z+r) /4.85
      return
171   plc =  x/(z+r) /4.85
      return
c
c  deflection auto- and cross covariances
c
18    plc = -.5*(log(z+r) + .5 + z/(z+r) + y**2/(z+r)**2) /23.5225
      return
19    plc = -.5*x*y/(z+r)**2 /23.5225
      return
20    plc = -.5*(log(z+r) + .5 + z/(z+r) + x**2/(z+r)**2) /23.5225
      return
      end
