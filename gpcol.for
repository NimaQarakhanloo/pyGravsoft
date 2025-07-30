      program gpcol
c $Id: gpcol.for 181 2008-08-14 23:19:51Z tjansson $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                            G P C O L
c
c  Program for flat-earth logarithmic least-squares 
c  collocation. One or more data files are read in and
c  predictions carried out in a grid. 
c  The program may also be used just for interactive listing 
c  of the analytical covariance functions Cgg, CgN, CNN (use nifiles=0).
c  An empirical cov.fkt. will be output for the first data kind.
c
c  input: nifiles            (number of input files, 0=list cov)
c         sqrtC0, D, T       (covariance parameters, mgal and km)
c         ki1, sig1          (input kind and noise)
c         <ifile1>           (first input file)
c         (dlim  -  for ki1 = 5 only)
c         ki2, sig2          (.. second input file ..)
c         <ifile2>
c         ...
c         ko, lhgrid, lblock         
c         (rlat1, rlat2, rlon1, rlon2, dlat, dlon, h  - for lhgrid=false)
c         (<hfile>  - for lhgrid=true)
c         (dlat, dlon, blat, blon, nmin  - for lblock=true)
c         <ofile>
c         <efile>            (only for lblock=f)
c
c
c  ki, ko: input/output type (1:geoid, 2: defl.pair, 3:dg, 6:ksi, 7:eta)
c          deflections given in mgal 
c          specials:
c          ki = 5: input geoid heights, but use slopes (along-track dfv)
c                  additional input: dlim (km) .. in-track max dist
c          ko = 0: just check blocks and emp.cov. of first data
c
c  sig1:   input data std.dev. If sig = -1 then the std.dev. must
c          be given in the datafile.
c
c  lhgrid  use heights in grid file as prediction point heights
c          whole grid predicted, except in blocked mode
c
c  lblock  make blocked computations in (dlat x dlon) blocks with blat,blon
c          borders (degrees). Output is in this case in text line format.
c
c  rlat1, rlat2, rlon1, rlon2, dlat, dlon, h: wanted grid and height
c
c  dlat,dlon,blat,blon   block size and border in degrees
c  nmin                  minimum number of obs in a block for a solution
c
c  Reference: A new covariance model for inertial gravimetry and
c  gradiometry. Journal of Geophysical Research vol. 92 pp. 1305-10, 1987.
c
c  (c) Rene Forsberg, KMS-Denmark
c  f77 version sep 1993
c  last updated dec 1, 1994 (sign error correction)
c               jan 29, 1996 (lhgrid/blocks/slopes - careful, slopes not tested)
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8 (a-h,o-z)
c
c  dimension: icdim = (maxobs*(maxobs+1))/2+maxobs
c
      parameter (maxobs=100,icdim=5015,maxrow=600)
      dimension c(icdim),rfi(maxobs),rla(maxobs),h(maxobs)
     .,d(maxobs),csx(maxobs),ki(maxobs),sig2(maxobs)
     .,pred(maxrow),stdv(maxrow),hinp(maxrow),kiix(6),sigmax(6),dlim(6)
      real*4 cosa(maxobs),sina(maxobs)
      dimension empsum(13),nemp(13)
      equivalence (sig2,csx)
      character*48 ifile(6), ofile, efile, hfile
      logical ldfv,lhgrid,lblock,lslope
c
      degkm = 111.11
      radeg = 180/3.14159265d0
      secfak = 206.265d0/degkm
c
      write(*,*) '--- G P C O L ---'
      write(*,*) 'input no of inputfiles (0 = list covariances) '
      read(*,*) nifiles
      if (nifiles.gt.6) stop '*** too many input files, max 6'
c
1     write(*,*) 'input sqrtC0 (mgal), D, T (km):    (0 0 0 = stop) '
      read(*,*) sqrtc0,d1,d2
      if (sqrtc0.eq.0) stop 
c
      var = sqrtc0**2
      write(*,2) sqrtc0,d1,d2
2     format(' Covariance model: sqrtC0,D,T = ',3f9.2)
c
c  write info on selected planar logarithmic model
c
      write(*,*)
     .'Covariances for gravity, geoid and along-track dfv ',
     .'for selected model:'
      write(*,*)
     .'  dist(km)  Cgg (mgal**2)  Cgn (mgal*m)   Cnn (m**2)',
     .'  Cee (arcsec**2)'
      ddmin = 9.9d9
      ihalf = -1
      do 4 i = 0, 120
        s = i
        cgg = covp(3,3,s,0.d0,0.d0,d1,d2,var)
        dd = abs(cgg-var/2)
        if (dd.lt.ddmin) then
          ihalf = i
          ddmin = dd
        endif
        if (mod(i,10).eq.0)
     .  write(*,3) s,cgg,covp(1,3,s,0.d0,0.d0,d1,d2,var),
     .  covp(1,1,s,0.d0,0.d0,d1,d2,var),
     .  covp(7,7,s,0.d0,0.d0,d1,d2,var)
3       format(f9.0,5f14.3)
4     continue
      write(*,*) 'Model correlation length in km: ',ihalf
      write(*,*)
      if (nifiles.eq.0) goto 1
c
c  read observation file names and remaining input
c  -----------------------------------------------
c
      do 5 i = 1, nifiles
        read(*,*) kiix(i),sigmax(i)
        read(*,'(a)') ifile(i)
        if (kiix(i).eq.5) read(*,*) dlim(i)
        write(*,6) ifile(i)
6       format(' Data input from file: ',a48)
5     continue
      read(*,*) ko,lhgrid,lblock
      if (lhgrid) then
        read(*,'(a)') hfile
        open(19,file=hfile,status='old')
        read(19,*) rfi1,rfi2,rla1,rla2,dfi,dla
      else
        read(*,*) rfi1,rfi2,rla1,rla2,dfi,dla,hh
      endif
      if (lblock) then
        read(*,*) dlat,dlon,blat,blon,nmin
        if (nmin.le.0) nmin = 1
      endif
      read(*,'(a)') ofile
      open(20,file=ofile,status='unknown')
      if (.not.lblock) then
        read(*,'(a)') efile
        open(21,file=efile,status='unknown')
      endif
c
c  initialize statistics and blocks
c  --------------------------------
c
      nptot = 0
      dsum = 0
      dsum2 = 0
      ndsum = 0
      do i = 1, 13
        empsum(i) = 0
        nemp(i) = 0
      enddo
c
      if (lblock) then
        ny = (rfi2-rfi1)/dlat + 0.5
        nx = (rla2-rla1)/dlon + 0.5
        if (ny.lt.1) ny = 1
        if (nx.lt.1) nx = 1
        write(*,7) blat,blon,ny,nx
7       format(' - Blocked computation, overlaps (deg):',2f7.2,
     .  ', blocks:',i4,' x',i4)
      else
        ny = 1
        nx = 1
      endif
c
c  block loop
c  ----------
c
      do 80 iy = ny, 1, -1
      do 80 ix = 1, nx
c
      sum = 0
      n = 1
      nin = 0
      if (lblock) then
        y1 = rfi1 + (iy-1)*dlat
        y2 = rfi1 + iy*dlat
        if (iy.eq.ny) y2 = rfi2+.0001
        ymin = y1-blat
        ymax = y2+blon
        x1 = rla1 + (ix-1)*dlon
        x2 = rla1 + ix*dlon
        if (ix.eq.nx) x2 = rla2+.0001
        xmin = x1-blat
        xmax = x2+blon
        write(*,10) y1,y2,x1,x2
10      format(/' SOLUTION BLOCK ',4f9.3)
      endif
c
c  loop for different input files
c  ------------------------------
c
      do 14 j = 1, nifiles 
        kii = kiix(j)
        if (j.eq.1) kii1 = kii
        sigma = sigmax(j)
        ldfv = (kii.eq.2)
        if (ldfv) kii = 6
        lslope = (kii.eq.5)
c
        open(10,file=ifile(j),status='old')
        dmin = 9.d9
        dmax = -9.d9
c
        if (lslope) then
          if (j.eq.1) write(*,*) 
     .    '- Geoid heights converted to slopes before collocation -'
          dlimd = dlim(j)/degkm
          ntr = 0
          rlatp = 999.99
          rlonp = 999.99
        endif
c
c  point loop entry
c
12      if (sigma.ge.0) then
          if (ldfv) then
            read(10,*,end=13) is,rlat,rlon,hh,dd,eta
          else
            read(10,*,end=13) is,rlat,rlon,hh,dd
          endif
          ss2 = sigma**2
        else
          if (ldfv) then
            read(10,*,end=13) is,rlat,rlon,hh,dd,eta,ss
          else
            read(10,*,end=13) is,rlat,rlon,hh,dd,ss
          endif
          ss2 = ss**2
        endif
c
        if (lblock) then
          if (rlat.lt.ymin.or.rlat.ge.ymax) goto 12
          if (rlon.lt.xmin.or.rlon.ge.xmax) goto 12
          if (y1.le.rlat.and.rlat.le.y2.and.
     .        x1.le.rlon.and.rlon.le.x2) nin = nin+1    
        endif
c
c  geoid heights to slopes
c
        if (lslope) then
          dy = rlat-rlatp
          dx = (rlon-rlonp)*cos(rlat/radeg)
          dr = sqrt(dx**2 + dy**2)
          if (dr.gt.dlimd) then
            ntr = ntr+1
            rlatp = rlat
            rlonp = rlon
            hhp = hh
            ddp = dd
            goto 12
          endif
          if (dr.le.0.0001) then
            write(*,*) '** warning: close points at: ',is,rlat,rlon
            goto 12
          endif
          if (n.le.maxobs) then
            rfi(n) = (rlat+rlatp)/2
            rla(n) = (rlon+rlonp)/2
            h(n) = (hhp+hh)/2
            d(n) = secfak*(dd-ddp)/dr
            cosa(n) = dy/dr
            sina(n) = dx/dr
            sig2(n) = ss2
            ki(n) = kii
            sum = sum+rfi(n)
            rlatp = rlat
            rlonp = rlon
            hhp = hh
            ddp = dd
          endif
c
c  save other obs in arrays
c
        elseif (n.le.maxobs) then
          rfi(n) = rlat
          rla(n) = rlon
          h(n) = hh
          d(n) = dd
          sig2(n) = ss2
          ki(n) = kii
          sum = sum+rfi(n)
        endif
c
        if (d(n).lt.dmin) dmin = d(n)
        if (d(n).gt.dmax) dmax = d(n)
        if (kii.eq.kii1) then
          dsum = dsum + d(n)
          dsum2 = dsum2 + d(n)**2
          ndsum = ndsum+1
        endif
c
        n = n+1
c
        if (ldfv) then 
          if (n.le.maxobs) then
            rfi(n) = rfi(n-1)
            rla(n) = rla(n-1)
            h(n) = h(n-1)
            d(n) = eta
            sig2(n) = sig2(n-1)
            ki(n) = 7 
          endif
          n = n+1
        endif
        goto 12
c
13    close(10)   
      if (lslope) write(*,*) 'Number of breaks/tracks in data: ',ntr
      write(*,131) dmin, dmax
131   format(' Min and max of input data: ',2f10.3)
14    continue
c
c  all files scanned
c
      n = n-1
c     
      if (n.gt.maxobs) then
        write(*,*) '- no computation, too many observations '   
        goto 80
      endif
      if (lblock) then
        write(*,15) n, nin
15      format(' Number of collocation obs.:',i5,
     .  ', in center block:',i5) 
        if (nin.lt.nmin) then
          write(*,*) '- not computed, too few observations in block'
          goto 80
        endif
      else
        write(*,*) 'Total number of observations input: ',n
        if (n.eq.0) stop 'no observations'
      endif 
c
      cosfi = cos(sum/n/57.29578d0)
      coskm = cosfi*degkm
c
c  set up collocation equations - do empcov of 1st data set at same time
c  ---------------------------------------------------------------------
c
      k = 0
      do 20 i = 1,n
      do 20 j = 1,i
        k = k+1
        kii = ki(i)
        kij = ki(j)
        xx = (rla(i)-rla(j))*coskm
        yy = (rfi(i)-rfi(j))*degkm
c
        if (kii.eq.kii1.and.kij.eq.kii1) then
          r = sqrt(xx**2+yy**2)
          ir = (r+5.0)/10+1
          if (ir.le.13) then
            empsum(ir) = empsum(ir) + d(i)*d(j) 
            nemp(ir) = nemp(ir)+1
          endif
          if (ko.eq.0) goto 20
        endif
c
        if (kii.ne.5.and.kij.ne.5) then
          c(k) = covp(kii,kij,xx,yy,h(i)+h(j),d1,d2,var)
        elseif (kii.eq.5) then
          if (kij.eq.5) then
            c(k) = 
     .      cosa(i)*cosa(j)*covp(6,6,xx,yy,h(i)+h(j),d1,d2,var) +
     .      sina(i)*sina(j)*covp(7,7,xx,yy,h(i)+h(j),d1,d2,var) +
     .      2*cosa(i)*sina(j)*covp(6,7,xx,yy,h(i)+h(j),d1,d2,var)
          else
            c(k) = 
     .      cosa(i)*covp(6,kij,xx,yy,h(i)+h(j),d1,d2,var) +
     .      sina(i)*covp(7,kij,xx,yy,h(i)+h(j),d1,d2,var)
          endif
        else 
          c(k) = 
     .    cosa(j)*covp(kii,6,xx,yy,h(i)+h(j),d1,d2,var) +
     .    sina(j)*covp(kii,7,xx,yy,h(i)+h(j),d1,d2,var)
        endif
        if (i.eq.j) c(k) = c(k)+sig2(i)
20    continue
      if (ko.eq.0) goto 80
      irhs = k
      do 21 i = 1,n
21    c(irhs+i) = d(i)
c
c  solve equations
c  ---------------
c
      call chol(c,n,nsing,1)
      if (nsing.ne.0) write(*,*)
     .'*** warning: singular equations - nsing = ',nsing
      write(*,*) '- Collocation equations solved -'
      do 23 i = 1, n
23    d(i) = c(irhs+i)
c
c  find predictions
c  ----------------
c
      nfi = (rfi2-rfi1)/dfi + 1.5
      nla = (rla2-rla1)/dla + 1.5
      if (nla.gt.maxrow) stop '*** too many rows - increase maxrow'
      ldfv = (ko.eq.2)
      if (ldfv) ko = 6
c
c  grid writing - entry for eta
c
28    if (.not.lblock) then
        write(20,29) rfi1,rfi2,rla1,rla2,dfi,dla
        write(21,29) rfi1,rfi2,rla1,rla2,dfi,dla
29      format(6f12.6)
      endif
c
      pmin = 99999.9
      pmax = -99999.9
      hmin = 99999.9
      hmax = -99999.9
      np = 0
c
      do 35 i = nfi,1,-1
      if (lhgrid) read(19,*) (hinp(j),j=1,nla)
      do 34 j = 1, nla
        pfi = rfi1 + (i-1)*dfi
        pla = rla1 + (j-1)*dla
        if (lhgrid) then
          hh = hinp(j)
          if (hh.lt.hmin) hmin = hh
          if (hh.gt.hmax) hmax = hh
        endif
c
c  predictions for point (pfi,pla)
c
        if (lblock) then
          if (pfi.lt.y1.or.pfi.ge.y2.or.pla.lt.x1.or.pla.ge.x2)
     .    goto 34
        endif
c
        sum = 0
        do 31 k = 1, n
          xx = (pla-rla(k))*coskm
          yy = (pfi-rfi(k))*degkm
          kik = ki(k)
          if (kik.ne.5) then
            cov = covp(ko,kik,xx,yy,h(k)+hh,d1,d2,var)
          else
            cov = 
     .      -cosa(k)*covp(ko,6,xx,yy,h(k)+hh,d1,d2,var) -
     .      sina(k)*covp(ko,7,xx,yy,h(k)+hh,d1,d2,var)
          endif
          sum = sum + cov*d(k)
          c(irhs+k) = cov
          csx(k) = cov
31      continue
        ppred = sum
c
c  standard deviations
c
        call chol(c,n,info,3)
        sum = 0
        do 32 k = 1, n
32      sum = sum + csx(k)*c(irhs+k)
        sum = covp(ko,ko,0.d0,0.d0,hh+hh,d1,d2,var) - sum
        if (sum.ge.0) pstdv = sqrt(sum)
        if (sum.lt.0) pstdv = -sqrt(-sum)
c
        np = np+1
        if (ppred.lt.pmin) pmin = ppred
        if (ppred.gt.pmax) pmax = ppred
        if (lblock) then
          is = (i-1)*nla+j
          write(20,33) is, pfi, pla, hh, ppred, pstdv
33        format(i7,2f8.4,f7.1,2f10.3)
        else
          pred(j) = ppred
          stdv(j) = pstdv
        endif
34    continue
      if (.not.lblock) then
        write(20,341) (pred(j),j=1,nla)
        write(21,341) (stdv(j),j=1,nla)
341     format(50(/,8f9.3))
      endif
35    continue
c
      write(*,50) np,pmin,pmax
50    format(' - Predicted:', i6,' points, min max = ',2f12.5)
      if (lhgrid) write(*,51) hmin,hmax
51    format(' - min/max of prediction point heights = ',2f9.1)
      if (ldfv) then
        ko = 7
        ldfv = .false.
        write(20,*)
        write(21,*)
        goto 28
      endif
      nptot = np + nptot
80    if (lhgrid) rewind(19)
      continue
c
      write(*,81) kii1,ndsum,dsum/ndsum,sqrt(dsum2/ndsum),dsum2/ndsum
81    format(/' Empirical covariances at 10 km intervals for ',
     .'data kind ',i1,
     ./' Total no of data:',i6,', mean, rms, C0 = ',2f9.3,f11.3
     ./'    Dist (km)    No of products      Cov')
      do i = 1, 13
        j = (i-1)*10
        r = 0
        if (nemp(i).gt.0) r = empsum(i)/nemp(i) 
        write(*,82) j,nemp(i),r
82      format('    ',i5,'       ',i9,4x,f13.3)
      enddo
      if (ko.eq.0) write(*,*) '- check of data completed w/o coll. -'     
c
      if (lblock) then
        write(*,*)
        write(*,*) '- Total number of predicted points:',nptot
      endif
      if (ko.ne.0) write(*,*) 'Collocation results output on: ',ofile 
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
c
      subroutine chol(c,n,nsing,imode)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      C H O L
c
c  subroutine solves positive definite symmetric linear equations
c  using cholesky decomposition. coefficients stored columnwise
c  in c, followed by righthand side. n is number of unknowns.
c  solution is returned as last column in c.
c  'nsing' is number of singularities. it should be zero for a 
c  succesful solution.
c  
c  special entry points (imode):
c
c  1    normal cholesky equation solver
c  
c  2    performs factorization of c without r.h.s.
c
c  3    factors and solves a particular r.h.s., the
c       r.h.s. is stored and the solution returned
c       in the last column of c, i.e. starting with 
c       element c(n*(n+1)/2 + 1). nsing is not used.
c
c  rf, nov 85. modified jan 86.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8 (a-h,o-z)
      dimension c(*)
c
c  imode: 1: chol, 2: cholfa, 3: cholsl
c
      goto (10,11,12),imode
c
10    nr1 = 1
      nr2 = n+1
      goto 15
c
11    nr1 = 1
      nr2 = n
      goto 15
c
12    nr1 = n+1
      nr2 = n+1
15    if (imode.ne.3) nsing = 0
c
      do 50 nr = nr1, nr2
        i=nr*(nr-1)/2
        ir=i
        do 40 nc = 1,nr
          sum=0
          ic=nc*(nc-1)/2
          i=i+1
          nc1=nc-1
          do 30 np = 1,nc1
30        sum = sum - c(ir+np)*c(ic+np)
          ci = c(i)+sum
          if (nr.eq.nc) goto 31
          c(i) = ci/c(ic+nc)
          goto 40
31        if (nr.gt.n) goto 40
          if (ci.gt.0) goto 32
            nsing = nsing+1
            c(i) = 9.9d39
            goto 40
32        c(i) = sqrt(ci)
40      continue
50    continue
      if (imode.eq.2) return
c
c  back substitution
c
      do 80 nc = n,1,-1
        i=i-1
        ir=i
        ic=nc*(nc+1)/2
        c(i) = c(i)/c(ic)
        do 70 np = nc-1,1,-1
          ir=ir-1
          ic=ic-1
          c(ir) = c(ir) - c(i)*c(ic)
70      continue
80    continue
      return
      end
