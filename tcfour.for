      program tcfour
c $Id: tcfour.for 243 2008-10-29 10:10:19Z cct $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                     t c f o u r
c
c  program for fft analysis of digital terrain models.
c  space domain expressions of the integral kernels are
c  transformed rather than using the fourier domain analytical
c  kernels, in order to allow control over integration
c  radii.
c
c  a subgrid of the grid in 'dtmfile1' is analyzed, specified
c  by its sw-corner (fic,lac) and number of
c  points. both 'in' and 'ie' must be  e v e n .
c  a reference height grid may be subtracted (for rtm effects).
c
c  input to program:
c
c      <dtmfile1>
c      <dtmfile2>
c      <ofile>
c      mode, lref, dist1, dist2 or height,
c      fic1, lac1, in1, ie1
c     (fic2, lac2, in2, ie2 - for mode 2 only)
c
c  'lref' is reference grid specification: f none, t subtract ref.
c      the reference grid name must be specified as 'dtmfile2', this file
c      is only opened if lref is true, and in mode 2.
c
c  'mode' determines the function of the program:
c
c  mode = 0: simple filtering. 'dist1' gives the wavelength
c            for cut-off in units of km.
c            'dist2' > 0: low pass, < 0: high pass.
c
c         1: covariance function and power spectra,
c            (2-d cov.fct. output as dtm with c(0,0)=1000)
c            'dist1, dist2' is not used.
c
c         2: terrain corrections, using grid1 to distance 'dist1',
c            grid2 to distance 'dist2' (both in km).
c            (the detailed grid is used out to a square of half side-
c            length 'dist1', border points in both inner and outer
c            grid being attenuated to compensate edge effects.)
c
c         3: terrain corrections from 'dist1' to 'dist2' (km).
c
c         4: rtm gravity effect, computed as terrain correction
c            to distance 'dist1' and bouguer reduction to reference
c            level.
c
c         5: geoid effect for residual topography, condensation
c            approximation. only one grid, second grid specification
c            is the reference grid (lref = false means no ref used).
c            effects will be computed in elevation 'dist2' (km) above
c            the topography (i.e., above the condensation level).
c            computed to distance 'dist1'.
c
c         6: isostatic geoid effects, second order expansion, in
c            elevation 'dist2' above the geoid. (mixed continental and
c            oceanic area using equivalent sea/rock root conversion).
c            computed to distance 'dist1' (km).
c
c         7: gravity effect of the isostatic compensation,
c            computed to distance 'dist1' (km).
c
c         8: rtm deflections of the vertical (condensation approximation)
c            (written on outfile as two consecutive grids in unit arcsec)
c            effect computed between distances 'dist1' to 'dist2' (km).
c
c  the program recognizes grids in utm by the label on 'dtmfile1'.
c  if utm grids are used (fic, lac) must be northing and easting in km.
c  if (fic, lac) = (0, 0) the sw-corner of the grid is used.
c  the program will automatically zero-pad the grids in not big enough
c  warning: this can give strange effects when reference grids are used!
c
c  overwiev of files:
c
c  file 20   elevation data grid (free format, initiated
c            with label lat1, lat2, lon1, lon2, dlat, dlon
c            specifying grid boundaries)
c  file 21   possible elevation grid for outer zones
c  file 22   possible reference elevation grid
c  file 32   output file (note that results close to the
c            edges are unreliable due to fft-periodicity)
c  file 30,31  scratch files, used for intermediate results
c
c  - note: results for mode 1 are listed with basic distance
c  unit degree.
c
c  (c) rene forsberg, national survey and cadastre of denmark (KMS)
c  first version 1984
c  modified for cdc, december 1984, rf
c  changed and updated for unsw vax, rf mar 1989
c  uppdated with zero padding, kms nov 94
c  updated for error in zero padding with reference dems, rf aug 2008
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      logical lgeog, lrun2, lg1, lref
      real*4 iha
      dimension nnarr(2),hlab(6)
      common /isopar/ d1,d2,diso,alfa
      common /utmpar/ lgeog,iell,izone
      character*128 dfile1,dfile2,ofile
c
c  fixed array bounds   ihadim = 160 x 160 = 25600
c  ipwdim = 160/2 + 1 = 81
c  iwkdim = 160*2 = 320
c  (larger area modification:
c  ihadim = 216 x 192 = 41472)
c
      dimension cha(2, 850000)
      dimension iha(   850000)
      dimension wrk(    2400)
      dimension pw(601), ipw(601), cf(601), icf(601)
      ihadim = 850000
      iwkdim = 2400
      ipwdim = 601
      idim2 = 2*ihadim
c
      read(*,1) dfile1
      read(*,1) dfile2
      read(*,1) ofile
1     format(a128)
c
      read(*,*) mode, lref, dist, dist2
      read(*,*) rfic, rlac, inn, ine
      if (mode.ne.2) goto 14
      read(*,*) rfic2, rlac2, inn2, ine2
14    continue
      iref = 0
      if (lref.and.mode.ne.2.and.mode.ne.3) iref = 1
c
c  open files - unit 30 initially as scratch file
c
      open(20,file=dfile1,status='old')
      if (mode.eq.2) open(21,file=dfile2,status='old')
      if (iref.eq.1) open(22,file=dfile2,status='old')
      open(30,status='scratch',form='unformatted')
      open(31,status='scratch',form='unformatted')
c
      write(*,2) dfile1,mode, dist, dist2, rfic, rlac, inn, ine
2     format(/' ---  T C F O U R  ---'
     ./' dtmfile: ',a36,/' mode = ',i3,', distance(s) = ',2f9.2,
     ./' sw corner: ',2f11.4,', points: ',2i5)
      if (iref.eq.1) write(*,3) dfile2
3     format(' reference grid file: ',a36)
      if (mode.eq.2) write(*,4) dfile2, rfic2, rlac2, inn2, ine2
4     format(' - two grid terrain correction computation -',
     ./' outer zones from file: ',a36,
     ./' sw corner: ',2f11.4,', points: ',2i5)
c
c  constants
c
      rho = 2.67
      rhoiso = 0.4
      diso = 32.0
      g = 0.00667
      grho = g*rho
c
      radeg = 1.0
c
c  check label on dtmfile1 for utm
c
      read(20,*) hlab
      lgeog = (abs(hlab(1)).lt.99.and.abs(hlab(2)).lt.99)
      if (.not.lgeog) radeg = 0.001
      rewind(20)
c
      lrun2 = .false.
c
c  convert degrees to radians
c  if meters, set radeg to convert to degrees
c  basic unit in listings always degrees
c
      rfic = rfic/radeg
      rlac = rlac/radeg
      if (mode.eq.2) rfic2 = rfic2/radeg
      if (mode.eq.2) rlac2 = rlac2/radeg
      if (.not.lgeog) radeg = 1.0/111200.0
c
c  read elevations - entry point for second loop
c  ---------------------------------------------
c
      iunit=20
15    continue
c
      lg1 = lgeog
      call rdgrid(iunit, rfic, rlac, inn, ine, dfi, dla, iha, 
     .ii1z, ii2z, jj1z, jj2z, ihadim)
      if (lg1.neqv.lgeog) stop 'grid1 and grid2 different utm/geog'
c
c  set fft constants, check dimension, find mean
c
      n = inn * ine
      nnarr(1) = ine
      nnarr(2) = inn
      nyqn = inn/2 + 1
      nyqe = ine/2 + 1
      if (n.le.ihadim) goto 17
      write(*,16)
16    format(' *** array size too small: wanted, declared ',2i7)
      stop
c
17    dlacos = dla
      if (lgeog) dlacos = dla * cos((rfic+inn/2*dfi)/
     .57.2957795)
      dydx = dfi/dlacos
      s = 0.0
      do 18 i = 1, n
        s = s + iha(i)
18    continue
      hmean = s/n
c
c  set isostatic root factor to continental or oceanic
c
       alfa = rho/rhoiso
       if (hmean.lt.0) alfa = (rho-1.03)/rhoiso
c
c  transform integral kernel and store on tmpfile
c  -----------------------------------------------
c  only real part needs to be stored except modes ge 6
c
      dfikm = dfi*radeg*111.2
      write(*,24) dfikm, dfikm/dydx
24    format(' gridspacings north, east: ',2f8.3, ' km')
      if (mode.le.1) goto 30
c
      d1 = dist2 - hmean/1000
      d2 = dist2 + diso + alfa*hmean/1000
c
      icase = mode
      if (mode.eq.2.and.lrun2) icase = 1
      j = 0
      if (mode.ge.6) j = 1
c
      call setcha(cha,ihadim,inn,ine,dydx,dfikm,
     .           icase,dist,dist2)
      if (n.le.60) call cprint(cha,ihadim,inn,ine)
      call fourt(cha,nnarr,2,-1,j,wrk,idim2,iwkdim)
      if (n.le.60) call cprint(cha,ihadim,inn,ine)
c
      if (lrun2) rewind(31)
      do 25 i = 1,n
        r = cha(1,i)/n
        s = cha(2,i)/n
        if (mode.lt.6) write(31) r
        if (mode.ge.6) write(31) r,s
25    continue
c
c  transfer elevations to complex array
c  ------------------------------------
c
30    call tofrom(cha,1,n,iha,ihadim)
c
c  read reference grid
c  -------------------
c
      if (iref.eq.0) goto 35
      reffi = 0.0
      refla = 0.0
      irefnn = 0
      irefne = 0
      lg1 = lgeog
      write(*,31)
31    format(/' reference grid:')
      call rdgrid(22, reffi, refla,
     .irefnn, irefne, refdfi, refdla, iha, ii1r,ii2r,jj1r,jj2r,ihadim)
c
      if (lgeog.neqv.lg1) stop 'reference grid and dtm grid utm diff'
      if (ii1r.gt.0.or.ii2r.gt.0.or.jj1r.gt.0.or.jj2r.gt.0) 
     .write(*,*) 'warning: refgrid zero padded: ',ii1r,ii2r,jj1r,jj2r
      call tofrom(cha,2,irefnn*irefne,iha,ihadim)
c
c  subtract reference values or store on unit 30
c
      rs = 0
      rs2 = 0
      do 34 i = inn, 1, -1
      do 34 j = 1, ine
        k = (inn-i)*ine + j
        rfi = rfic + (i-1)*dfi
        rla = rlac + (j-1)*dla
        r = (rfi-reffi)/refdfi + 1.0
        s = (rla-refla)/refdla + 1.0
        if (i.le.inn-ii2z.and.i.ge.1+ii1z.
     .  and.j.le.ine-jj2z.and.j.ge.1+jj1z) then        
          refh = bilinc(r, s, cha, ihadim, irefnn, irefne, 2)
        else
          refh = 0
        endif
c
        if (mode.eq.4) then
          write(30) refh
        else
          rr = cha(1,k) - refh
          cha(1,k) = rr
          rs = rr + rs
          rs2 = rr**2 + rs2
        endif
34    continue
      if (mode.ne.4) write(*,341) rs/n, sqrt((rs2-rs**2/n)/(n-1))
341   format(' mean and std.dev. after ref subtracted: ',2f9.2)
c
c  remove data mean for modes 1 - 4
c  multiply by density for geoid modes
c  set h**2 to obtain both transforms simultanously
c  ------------------------------------------------
c
35    s = 0.0
      do 38 i = 1, n
        r = cha(1,i)
        if (mode.ge.1.and.mode.le.4) r = r - hmean
        cha(1,i) = r
        cha(2,i) = 0.0
        if (mode.ge.2.and.mode.ne.5.and.mode.ne.8) cha(2,i) = r**2
        if (mode.lt.5) goto 37
          r0 = rho
          if (r.lt.0.and.mode.ne.5.and.mode.ne.8) r0 = rho - 1.03
          cha(1,i) = cha(1,i)*r0
          cha(2,i) = cha(2,i)*r0
37      s = s + cha(1,i)**2 + cha(2,i)**2
38    continue
      s = s/n
      write(*, 39) s, hmean
39    format( ' power complex space domain ',e10.4,', mean ',f9.2)
c
c  fourier transformation of elevations
c  ------------------------------------
c  discriminate between complex (j=1) and real (j=0) data
c
      j = 1
      if (mode.le.1.or.mode.eq.5.or.mode.eq.8) j = 0
      call fourt(cha,nnarr,2,-1,j,wrk,idim2,iwkdim)
      if (n.le.60) call cprint(cha,ihadim,inn,ine)
c
      s = 0
      do 41 i = 1,n
        cha(1,i) = cha(1,i)/n
        cha(2,i) = cha(2,i)/n
        s = s + cha(1,i)**2 + cha(2,i)**2
41    continue
      write(*, 42) s
42    format(' power complex freq. domain ',e10.4)
c
c  mode 0
c  ------
c
c  simple filtering of elevations
c
      if (mode.ne.0) goto 50
      if (dist2.ne.0) call setcha(cha,ihadim,inn,ine,
     .dydx,dfikm,0,dist,dist2)
      if (n.le.60) call cprint(cha,
     .ihadim, inn, ine)
      call fourt(cha,nnarr,2,1,1,wrk,idim2,iwkdim)
      if (n.le.90) call cprint(cha, ihadim, inn, ine)
      goto 900
c
c  mode 1
c  ------
c
c  power spectrum and covariance functions
c
50    if (mode.ne.1) goto 80
      do 51 i=1,n
        cha(1,i) = cha(1,i)**2 + cha(2,i)**2
        cha(2,i) = 0.0
51    continue
c
      r = ine/(dydx*inn)
      call azsmoo(cha,ihadim,inn,ine,r,
     .pw,ipw,ipwdim)
c
      call fourt(cha,nnarr,2,1,1,wrk,idim2,iwkdim)
      if (n.le.90) call cprint(cha, ihadim, inn, ine)
      call azsmoo(cha,ihadim,inn,ine,
     .dydx,cf,icf,ipwdim)
c
c  print out covariance function and power spectra
c  power spectra given as logaritmic gravity effect
c
      write(*, 59) dfi*radeg*60, 1/(inn*dfi*radeg), cf(1)
59    format(/' spacing ',f8.3,' arcmin, freq.spacing. ',
     .f8.4,' cycles/degree' / 'c0 = ',f13.4)
      write(*, 60)
60    format(' <   n    psi      power(db)         cov')
      pfakt = inn*dfi*ine*dlacos*(radeg**2)
      pw(1) = pw(1) + hmean**2
      do 62 i = 1, nyqn
        pw(i) = pfakt*pw(i)
        pw(i) = 4.342945 * log(pw(i)) - 19.02
        write(*, 61) i-1, (i-1)*dfi*radeg*60, pw(i),
     .  cf(i)/cf(1), icf(i)
61      format(i6,f7.2,2f14.4,i8)
62    continue
      goto 900
c
c  mode 2, 3 and 4 - terrain corrections
c  -------------------------------------
c
c  terrain correction computation
c
80    if (mode.gt.4) goto 200
      if (mode.eq.2.or.mode.eq.3) write(*,81) dist,dist2
      if (mode.eq.4) write(*, 82) dist
81    format(/' terrain correction zone ',2f7.1, ' km')
82    format(/' terrain correction zone out to ',f6.1,' km')
c
c  input filter from outfile and convolve
c
c  rewind unit 31
c
      rewind(31)
      do 83 i = 1, n
        read(31) r
        if (i.eq.1) r0 = r
        cha(1,i) = cha(1,i)*r
        cha(2,i) = cha(2,i)*r
83    continue
c
c  obtain space domain convolutions for h and h**2
c
      call fourt(cha,nnarr,2,1,1,wrk,idim2,iwkdim)
c
c  compute final corrections in fixed unit microgal
c  add bouguer contribution for mode 4
c  exit for completion of second loop results
c
      if (mode.eq.4) rewind(30)
      cmax = -999999.0
      cmin = 999999.0
      c1max = -999999.0
      c2max = -999999.0
      s = 0.5*grho*dfikm*(dfikm/dydx)*inn*ine
c
      do 84 i = 1, n
        r = iha(i) - hmean
        c1 = s*(cha(2,i) - 2*r*cha(1,i) + r**2*r0)
        if (c1.gt.cmax) cmax = c1
        if (.not.lrun2) cha(1,i) = c1
        if (mode.ne.4) goto 84
        read(30) rr
        r = iha(i)-rr
      write(*,*) 'test i iha ref',i,iha(i),rr
        cha(1,i) = 2000*3.1415926*grho*r - c1
        if (cha(1,i).lt.cmin) cmin = cha(1,i)
        if (cha(1,i).gt.c1max) c1max = cha(1,i)
84    continue
      write(*, 85) cmax/1000
      if (mode.eq.4) write(*, 87) cmin/1000, c1max/1000
87    format(' minimal and maximal rtm effects: ',2f11.6,' mgal')
85    format(' maximal grid terrain correction: ',f11.6,' mgal')
c
c  exit if one-stage terrain correction
c
      if (mode.ne.2) goto 900
      if (lrun2) goto 100
c
c  output inner zone results and heights temporarily on outfile
c
      do 90 i = 1, n
        rr = iha(i)
        write(30) rr, cha(1,i)
90    continue
c
c  set second loop parameters
c
      lrun2 = .true.
      iunit=21
      inn1 = inn
      ine1 = ine
      n1 = n
      dfi1 = dfi
      dla1 = dla
      rfic1 = rfic
      rlac1 = rlac
      inn = inn2
      ine = ine2
      rfic = rfic2
      rlac = rlac2
      ii0 = dist/dfikm
      jj0 = dist/(dfikm/dydx)
      goto 15
c
c  complete terrain corrections loop2
c  rewind temp store on unit 30 and input grid 1 results
c
100    rewind(30)
      do 103 i = 1, inn1
      do 103 j = 1, ine1
        k = (i-1)*ine1 + j
        read(30) rih, c1
        rfi = (inn1-i)*dfi1 + rfic1
        rla = (j-1)*dla1 + rlac1
        ri = (rfi-rfic)/dfi + 1.0
        rj = (rla-rlac)/dla + 1.0
        cha1 = bilinc(ri,rj,cha,ihadim,inn,ine,1)
        cha2 = bilinc(ri,rj,cha,ihadim,inn,ine,2)
        r = rih - hmean
        c2 = s*(cha2 - 2*r*cha1 + r**2*r0)
        iha(k) = c1 + c2
        if (i.lt.ii0.or.i.gt.inn1-ii0) goto 103
        if (j.lt.jj0.or.j.gt.ine1-jj0) goto 103
        if (c1.gt.c1max) c1max = c1
        if (c2.gt.c2max) c2max = c2
        if (iha(k).gt.cmax) cmax = iha(k)
103   continue
      write(*, 104) c1max/1000, c2max/1000, cmax/1000
104   format(' --- maximal terrain corrections reliable area:'/
     .' inner',f10.4,', outer',f10.4,', total',f10.4)
      inn = inn1
      ine = ine1
      dfi = dfi1
      dla = dla1
      rfic = rfic1
      rlac = rlac1
      call tofrom(cha,1,n1,iha,ihadim)
      goto 900
c
c  mode 5, 6, 7 - geoid and gravity
c  --------------------------------
c
200   if (mode.gt.7) goto 300
      if (mode.eq.5) write(*, 201) dist, dist2
201   format(/' rtm geoid effects maxdist, elevation (km): ',2f7.1)
      if (mode.eq.6) write(*, 202) dist, dist2, diso, rhoiso
202   format(/' isostatic geoid effects maxdist, elev, isopar: ',
     .4f6.1)
      if (mode.eq.7) write(*,2021) dist,dist2,diso,rhoiso
2021  format(/' isostatic correction maxdist, elev, isopar: ',
     .4f6.1)
c
c  input filter from tmpfile and convolve
c
      rewind(31)
      if (mode.ge.6) goto 204
      do 203 i = 1, n
        read(31) r
        cha(1,i) = cha(1,i) * r
        cha(2,i) = cha(2,i) * r
203   continue
      goto 210
c
c  resolve complex spectra of h and h**2
c
204   do 205 i = 1, inn
      do 205 j = 1, ine
        k = (i-1)*ine + j
        ii = inn+2-i
        if (ii.gt.inn) ii = ii-inn
        jj = ine+2-j
        if (jj.gt.ine) jj = jj-ine
        kk = (ii-1)*ine + jj
        rr = (cha(1,k)+cha(1,kk))/2
        ri = (cha(2,k)-cha(2,kk))/2
        rr2 = (cha(2,kk)+cha(2,k))/2000
        ri2 = (cha(1,kk)-cha(1,k))/2000
c
        read(31) r, s
        rr = rr*r + rr2/2*s
        ri = ri*r + ri2/2*s
        write(30) rr, ri
205   continue
      rewind(30)
      do 206 i = 1, n
        read(30) rr, ri
        cha(1,i) = rr
        cha(2,i) = ri
206   continue
c
210   call fourt(cha,nnarr,2,1,1,wrk,idim2,iwkdim)
c
      cmax = -999999.9
      c1max = 999999.9
c
c  constant for geoid undulations in mm
c
      s = 1.020*g*dfikm*(dfikm/dydx)*inn*ine
      if (mode.eq.7) s = 1000*g*dfikm*(dfikm/dydx)*inn*ine
      do 220 i = 1, n
        cha(1,i) = cha(1,i)*s
        if (cha(1,i).lt.c1max) c1max = cha(1,i)
        if (cha(1,i).gt.cmax) cmax = cha(1,i)
220   continue
      write(*, 221) c1max/1000, cmax/1000
221   format(' minimum and maximum effects: ',2f12.6)
      goto 900
c
c  mode 8 - deflections of the vertical
c  ------------------------------------
c
300   if (mode.gt.8) goto 900
      write(*, 301) dist, dist2
301   format(/' rtm deflections  zone (km): ',2f7.1)
c
c  input filter from tmpfile and restore components
c
      rewind(31)
      do 310 i = 1, n
        write(30) cha(1,i), cha(2,i)
        read(31) cha(1,i), cha(2,i)
310   continue
      rewind(31)
      do 315 i = 1, inn
      do 315 j = 1, ine
        k = (i-1)*ine+j
        ii = inn+2-i
        if (ii.gt.inn) ii = ii-inn
        jj = ine+2-j
        if (jj.gt.ine) jj = jj-ine
        kk = (ii-1)*ine + jj
        rr = (cha(1,k)+cha(1,kk))/2
        ri = (cha(2,k)-cha(2,kk))/2
        rr2 = (cha(2,kk)+cha(2,k))/2
        ri2 = (cha(1,kk)-cha(1,k))/2
        write(31) rr2, ri2, rr, ri
315   continue
c
c  restore transformed elevations and convolve in two steps
c
      lrun2 = .false.
320   rewind(31)
      rewind(30)
      do 325 i = 1, n
        read(31) rr, ri, rr2, ri2
        read(30) r, s
        if (lrun2) goto 321
        cha(1,i) = rr*r - ri*s
        cha(2,i) = rr*s + ri*r
        goto 325
321     cha(1,i) = rr2*r - ri2*s
        cha(2,i) = rr2*s + ri2*r
325   continue
c
c  fourier transformation
c
      call fourt(cha,nnarr,2,1,1,wrk,idim2,iwkdim)
      cmax = -999999.9
      cmin = 999999.9
      s = 210.34*g*dfikm**2/dydx*inn*ine
      do 330 i = 1, n
        cha(1,i) = cha(1,i)*s
        if (cha(1,i).lt.cmin) cmin = cha(1,i)
        if (cha(1,i).gt.cmax) cmax = cha(1,i)
        if (.not.lrun2) iha(i) = cha(1,i)
330   continue
      if (lrun2) write(*, 332) cmin/1000, cmax/1000
      if (.not.lrun2) write(*, 331) cmin/1000, cmax/1000
331   format(' minimal and maximal effects: eta ',2f12.6)
332   format('                              ksi ',2f12.6)
      if (lrun2) goto 900
      lrun2 = .true.
      goto 320
c
c  write out contents of cha(1,.) on outfile
c  -----------------------------------------
c  center for covariance functions
c
900   continue
c
      close(30)
      close(31)
      open(32,file=ofile,status='unknown')
c
      lrun2 = .false.
c
901   rr = rfic+(inn-1)*dfi
      ri = rlac+(ine-1)*dla
      if (lgeog) write(32, 902) rfic+ii1z*dfi,rr-ii2z*dfi,
     .rlac+jj1z*dla,ri-jj2z*dla,dfi,dla
902   format(' ',4f12.6,2f12.8)
      if (.not.lgeog) write(32, 903) rfic+ii1z*dfi,rr-ii2z*dfi,
     .rlac+jj1z*dla,ri-jj2z*dla,dfi,dla,iell,izone
903   format(' ',6f12.0/' ',2i6)
      ishift = 0
      jshift = 0
      if (mode.eq.1) ishift = inn/2
      if (mode.eq.1) jshift = ine/2+1
      do 910 ii = 1, inn
        do 909 jj = 1, ine
          i = ii + ishift
          j = jj + jshift
          if (i.gt.inn) i = i-inn
          if (j.gt.ine) j = j-ine
          wrk(jj) = cha(1,(i-1)*ine+j)/1000.0
909     continue
        if (ii.gt.ii2z.and.ii.le.inn-ii1z)
     .  write(32, 908) (wrk(jj),jj=1+jj1z,ine-jj2z)
908     format(30(/,' ',9f8.3))
910   continue
      if (mode.ne.8.or.lrun2) goto 912
      lrun2 = .true.
      do 911 i = 1, n
911   cha(1,i) = iha(i)
      goto 901
c
912   end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                    c p r i n t
c
c  prints the contents of complex array cha
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cprint(cha, ihadim, inn, ine)
      implicit double precision(a-h,o-z)
      dimension cha(2,ihadim)
      do 19 k=1,2
        j1=1
        j2=ine
        write(*,11)
11      format(' ')
        do 18 i=1,inn
          write(*,12) (cha(k,j),j=j1,j2)
12        format(8f11.2)
          j1 = j1+ine
          j2 = j2+ine
18      continue
19    continue
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                    t o f r o m
c
c  transfers from integer to complex array
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tofrom(cha,k,n,iha,ihadim)
      implicit double precision(a-h,o-z)
      real*4 iha
      dimension iha(ihadim), cha(2,ihadim)
      do 10 i = 1, n
        cha(k,i) = iha(i)
10    continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                   a z s m o o
c
c  averages power spectrum or covariance function in circles
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine azsmoo(cha,ihadim,inn,ine,
     .                  dydx,cf,icf,icfdim)
      implicit double precision(a-h,o-z)
      dimension cha(2,ihadim)
      dimension icf(icfdim), cf(icfdim)
      nyqn = inn/2+1
      nyqe = ine/2+1
      do 10 i = 1, nyqn
        cf(i) = 0.0
        icf(i) = 0
10    continue
c
      do 18 i = 1, nyqn
      do 18 j = 1, nyqe
        ii = inn+2-i
        jj = ine+2-j
        k11 = (i-1)*ine+j
        k12 = (i-1)*ine+jj
        k21 = (ii-1)*ine+j
        k22 = (ii-1)*ine+jj
        r = sqrt((i-1.0)**2 + ((j-1)/dydx)**2)
        k = int(r+0.5) + 1
        if (k.gt.nyqn) goto 17
c
        if (i.gt.1.or.j.gt.1) goto 11
        cf(k) = cf(k) + cha(1,k11)
        icf(k) = icf(k)  + 1
        goto 17
11      if (i.gt.1) goto 12
        cf(k) = cf(k) + cha(1,k11) + cha(1,k12)
        icf(k) = icf(k) + 2
        goto 17
12      if (j.gt.1) goto 13
        cf(k) = cf(k) + cha(1,k11) + cha(1,k21)
        icf(k) = icf(k) + 2
        goto 17
13      cf(k) = cf(k) + cha(1,k11) + cha(1,k21) +
     .                  cha(1,k12) + cha(1,k22)
        icf(k) = icf(k) + 4
17      continue
18    continue
      do 19 i = 1,nyqn
        cf(i) = cf(i)/icf(i)
19    continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                     b i l i n c
c
c  interpolates values in a dtm array using bilinear
c  (parabolic hyperboloid) interpolation
c  argument ri, rj has (1,1) in lower left corner.
c  interpolation outside grid is only allowed to 0.6 of a
c  gridspacing.
c
c  rf, june 84, updated mar 89.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function bilinc(ri, rj, cha, ihadim, imax, jmax, kk)
      implicit double precision(a-h,o-z)
      dimension cha(2, ihadim)
      in = ri
      ie = rj
      rn = ri - in
      re = rj - ie
      if (in.ge.1) goto 10
      if (ri.lt.0.4) stop 'linear interpolation outside grid (s)'
      in = 1
      rn = 0.0
10    if (in.lt.imax) goto 11
      if (ri.gt.imax+0.6) stop 'linear interpolation outside grid (n)'
      in = imax-1
      rn = 1.0
11    if (ie.ge.1) goto 12
      if (rj.lt.0.4) stop 'linear interpolation outside grid (e)'
      ie = 1
      re = 0.0
12    if (ie.lt.jmax) goto 13
      if (rj.gt.jmax+0.6) stop 'linear interpolation outside grid (w)'
      ie = jmax-1
      re = 1.0
13    k = (imax-in)*jmax+ie
      r = (1-rn)*(1-re)*cha(kk,k) +
     .rn*(1-re)*cha(kk,k-jmax) + (1-rn)*re*cha(kk,k+1) +
     .rn*re*cha(kk,k-jmax+1)
      bilinc = r
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                   s e t c h a
c
c  subroutine for setting array cha with a symmetrically
c  extended argument, determined by 'icase':
c
c  icase:  0:  high pass/low pass to distance 'dist'
c          1:  outer zones tc, square inner border
c          2:  inner zones tc to square
c          3:  1/r**3 out to distance 'dist'
c          4:  1/r**3 from distance 'dist' to 'dist2'
c          5:  1/r in elevation 'dist2'
c          6:  geoid isostasy (complex)
c          7:  gravity isostasy (compensation only)
c          8:  horizontal derivatives (complex)
c
c  'dfi' should be given in units of km.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine setcha(cha,ihadim,inn,ine,dydx,dfikm,
     .                  icase,dist,dist2)
      implicit double precision(a-h,o-z)
      dimension cha(2,ihadim)
      common /isopar/ d1,d2,diso,alfa
c
      z2 = 0.0
      if (icase.le.4.or.icase.eq.8) goto 5
c
c  cylindrical correction radius for distance zero
c
      z2 = dist2**2
      if (icase.eq.6) z2 = d1**2
      if (icase.eq.7) z2 = d2**2
      dlim = sqrt(dist**2 + z2)
      ri = 0.56*dfikm/sqrt(dydx)
      rj = sqrt(ri**2 + z2)
5     nyqn = inn/2 + 1
      nyqe = ine/2 + 1
      if (icase.eq.0) rr = 2.0/dist
c
      do 30 i = 1, nyqn
      do 30 j = 1, nyqe
        ii = inn+2-i
        jj = ine+2-j
        k = (i-1)*ine+j
        k1 = (i-1)*ine+jj
        k2 = (ii-1)*ine+j
        k3 = (ii-1)*ine+jj
        xx = (i-1.0)*dfikm
        yy = (j-1.0)/dydx*dfikm
        if (icase.gt.0)
     .  r = sqrt(xx**2 + yy**2 + z2)
        goto (10,11,12,13,14,15,16,17,18),icase+1
c
c  0: high pass/low pass filtering
c
10      r = sqrt(((i-1)/(inn*dfikm))**2 + ((j-1)/(ine*dfikm/dydx))**2)
        if (dist2.gt.0.and.r.le.rr) goto 30
        if (dist2.lt.0.and.r.gt.rr) goto 30
        goto 28
c
c  1: terrain corrections for coarse grid outer area
c
11      if (r.gt.dist2) goto 28
        ri = dist/dfikm+1.5
        rj = dist/(dfikm/dydx)+1.5
        ii = ri
        jj = rj
        ri = 1.0 - (ri-ii)
        rj = 1.0 - (rj-jj)
        if (i.lt.ii.and.j.lt.jj) goto 28
        rr = 1/r**3
        if (i.eq.ii) rr = ri * rr
        if (j.eq.jj) rr = rj * rr
        goto 27
c
c  2:  terrain corrections for square inner zone
c
12     if (i.eq.1.and.j.eq.1) goto 28
       ri = dist/dfikm + 1.5
       rj = dist/(dfikm/dydx) + 1.5
       ii = ri
       jj = rj
       ri = ri - ii
       rj = rj - jj
       if (i.gt.ii.or.j.gt.jj) goto 28
       rr = 1/r**3
       if (i.eq.ii) rr = ri * rr
       if (j.eq.jj) rr = rj * rr
       goto 27
c
c  3: terrain correction from 'dist1' to 'dist2'
c
13      if (r.le.dist.or.r.gt.dist2) goto 28
        if (i.eq.1.and.j.eq.1) goto 28
        rr = 1/r**3
        goto 27
c
c  4: terrain correction to 'dist'
c
14      if (r.gt.dist) goto 28
        if (i.eq.1.and.j.eq.1) goto 28
        rr = 1/r**3
        goto 27


c
c  5: rtm geoid effects condensation
c
15      if (r.gt.dlim) goto 28
        if (i.eq.1.and.j.eq.1) goto 151
        rr = 1/r
        goto 27
151     rr = 2*3.141593*(rj-dist2)/dfikm/(dfikm/dydx)
        goto 27
c
c  6:  complex geoid isostasy
c
16      if (r.gt.dlim) goto 28
        ri = sqrt(xx**2 + yy**2 + d2**2)
        if (i.eq.1.and.j.eq.1) goto 161
        rr = 1/r - 1/ri - d1*(dist2-d1)/r**3
     .       - d2*(d2-dist2-diso)/ri**3
        ri = d1/r**3 + alfa*d2/ri**3
        goto 26
161     r = 0.0
        if (rj.ne.0) r = 1-d1/rj
        rr = 2*3.141593*(rj - d1 - (dist2-d1)*r)/dfikm
     .       /(dfikm/dydx) - 1/ri - d2*(d2-dist2-diso)/ri**3
        ri = 2*3.141593*r/dfikm/(dfikm/dydx)
     .       + alfa*d2/ri**3
        goto 26
c
c  7:  gravity isostasy
c
17      if (r.gt.dlim) goto 28
        rj = (r**2 - 3*d2**2)/r**5
        rr = -(d2/r**3 + rj*(dist2+diso-d2))
        ri = -alfa*rj
        goto 26
c
c  8:  horizontal gravity, anisotropic filter
c
18      if (r.le.dist) goto 28
        if (r.gt.dist2) goto 28
        if (i.eq.1.and.j.eq.1) goto 28
        rr = -xx/r**3
        ri = yy/r**3
        cha(1,k) = rr
        if (j.gt.1) cha(1,k1) = rr
        if (i.gt.1) cha(1,k2) = -rr
        if (i.gt.1.and.j.gt.1) cha(1,k3) = -rr
        cha(2,k) = ri
        if (j.gt.1) cha(2,k1) = -ri
        if (i.gt.1) cha(2,k2) = ri
        if (i.gt.1.and.j.gt.1) cha(2,k3) = -ri
        goto 30
c
c  extend symmetrically for even function
c
26      cha(2,k) = ri
        if (j.gt.1) cha(2,k1) = ri
        if (i.gt.1) cha(2,k2) = ri
        if (i.gt.1.and.j.gt.1) cha(2,k3) = ri
27      cha(1,k) = rr
        if (j.gt.1) cha(1,k1) = rr
        if (i.gt.1) cha(1,k2) = rr
        if (i.gt.1.and.j.gt.1) cha(1,k3) = rr
        if (icase.lt.6) goto 29
        goto 30
c
c  zero values
c
28      cha(1,k) = 0.0
        if (j.gt.1) cha(1,k1) = 0.0
        if (i.gt.1) cha(1,k2) = 0.0
        if (i.gt.1.and.j.gt.1) cha(1,k3) = 0.0
c
29      cha(2,k) = 0.0
        if (j.gt.1) cha(2,k1) = 0.0
        if (i.gt.1) cha(2,k2) = 0.0
        if (i.gt.1.and.j.gt.1) cha(2,k3) = 0.0
30    continue
      return
      end
c
      subroutine rdgrid(iunit, rfic, rlac, inn, ine, dfi, dla,
     .iha, ii1z,ii2z,jj1z,jj2z,idim)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                       r d g r i d
c
c  subroutine for reading a digital grid file on
c  standard format, i.e. stored rowwise from nw to se, with label.
c
c  as se-corner coordinate 'rfic, rlac' (degrees) will be used
c  (unless they are zero, then the grid corner is used).
c  a grid containing 'inn' x 'ine' points will be put in array
c  'cha' of declared dimension 'idim'.
c  if inn=0 the complete grid will be read.
c  if the wanted grid is too large a zero padding will be done.
c
c  last updated jun 90, rf
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      common /utmpar/ lgeog,iell,izone
      real*4 iha
      dimension iha(idim)
      dimension hlab(6), hrow(2000)
      logical lgeog, lutm
      irdim = 2000
c
c  initialize statistics
c
      nr = 0
      rsum = 0.0
      rsum2 = 0.0
      rmin = 99999
      rmax = -99999
c
      read(iunit,*) (hlab(j),j=1,6)
      lutm = (abs(hlab(1)).ge.400.or.abs(hlab(2)).ge.400)
      lgeog = (.not.lutm)
      if (lgeog) goto 111
      read(iunit,*) iell,izone
      write(*,110) iell,izone
      if (iell.lt.1.or.iell.gt.3.or.izone.lt.1.or.izone.gt.60)
     *stop 'illegal ellipsoid or utm zone'
110   format(' - input grid in utm, ell ',i1,' zone ',i2,' -')
111   dfi = hlab(5)
      dla = hlab(6)
      nn = (hlab(2)-hlab(1))/dfi+1.5
      ne = (hlab(4)-hlab(3))/dla+1.5
      if (nn.gt.irdim.or.ne.gt.irdim) stop 'too long rows/columns'
c
c  find corner indices for wanted subgrid
c
      if (inn.eq.0) then
        inn = nn
        ine = ne
      endif
      if (rfic.eq.0.and.rlac.eq.0) then
        rfic = hlab(1)
        rlac = hlab(3)
      endif
      ifi1 = (rfic-hlab(1))/dfi+1.5
      ila1 = (rlac-hlab(3))/dla+1.5
      ifi2 = ifi1+inn-1
      ila2 = ila1+ine-1
      rfic = (ifi1-1)*dfi + hlab(1)
      rlac = (ila1-1)*dla + hlab(3)
      n = inn*ine
c
c  check boundaries for padding 
c
      ii1z = 0
      ii2z = 0
      jj1z = 0
      jj2z = 0
      if (ifi1.lt.1) ii1z = 1-ifi1
      if (ifi2.gt.nn) ii2z = ifi2-nn
      if (ila1.lt.1) jj1z = 1-ila1
      if (ila2.gt.ne) jj2z = ila2-ne
c
      if (n.gt.idim) then
        write(*, 122) n,idim
122     format(' *** array dim too small - wanted, declared ',2i8)
        stop ' *** sorry ***'
      endif
c
c  read data grid values
c  data in iha array stored with first element at nw corner
c
      ilast = ifi1
      if (ilast.lt.1) ilast = 1
c
      do 130 i = nn,ilast,-1
c
        read(iunit,*,end=131) (hrow(j),j=1,ne)
c
        if (i.lt.ifi1.or.i.gt.ifi2) goto 130
        jj0 = (ifi2-i)*ine - ila1+1
        do 129 j = 1,ne
          r = hrow(j)
          if (j.lt.ila1.or.j.gt.ila2) goto 129
          iha(j+jj0) = r
          nr = nr + 1
          if (r.gt.rmax) rmax = r
          if (r.lt.rmin) rmin = r
          rsum = rsum + r
          rsum2 = rsum2 + r**2
129     continue
130   continue
      goto 133
131     write(*,132) i
132     format(' *** too few data in grid file, lastrow = ',i7)
        stop ' *** check grid label and data ***'
c
c  zero padding
c
133   if (ii1z+ii2z+jj1z+jj2z.gt.0) then
        do 138 i = inn,1,-1
          jj0 = (inn-i)*ine
          if (i.gt.ifi2.or.i.lt.ifi1) then
            do 134 j = 1, ine
134         iha(j+jj0) = 0 
          else
            do 135 j = 1, jj1z
135         iha(j+jj0) = 0
            do 136 j = ine-jj2z+1, ine
136         iha(j+jj0) = 0
          endif
138     continue
      endif
c
c  write information and statistics
c
      if (nr.eq.0) stop '*** no points read from grid, wrong area'
      rfi = hlab(1) + (ifi1-1)*dfi
      rla = hlab(3) + (ila1-1)*dla
      r = rsum/nr
      s = 0.0
      if (n.gt.1)
     .s = sqrt((rsum2 - rsum**2/nr)/(nr-1))
      if (lgeog) write(*,141) (hlab(j),j=1,6),nn,ne
      if (lutm) write(*,142) (hlab(j),j=1,6),nn,ne
141   format(' Gridlab:',4f10.4,2f9.4,i5,i4)
142   format(' Gridlab:',4f10.0,2f8.0,i5,i4)
      if (lgeog) write(*, 143) rfi, rla, inn, ine
143   format(' Selected: sw corner ',2f10.4, ', points ', 2i6,i8)
      if (lutm) write(*, 144) rfi, rla, inn, ine
144   format(' Selected: sw corner ',2f10.0, ', points ', 2i6,i8)
      write(*, 145) nr, r, s, rmin, rmax
145   format(' Statistics of data selected from input grid:'
     ./' pts mean std.dev. min max:',i8,4f9.2)
      if (ii1z+ii2z+jj1z+jj2z.gt.0) write(*,146) ii1z,ii2z,jj1z,jj2z
146   format(' Zero padding done on grid, no of rows/cols S/N/E/W:',4i4)
      return
      end
c
      subroutine fourt(datt,nn,ndim,isign,iform,work,
     .idim1,idim2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                     f o u r t
c
c        version=740301
c        program description norsar n-pd9 dated 1 july 1970
c        author n m brenner
c        further description    three fortran programs etc.
c        issued by lincoln laboratory, mit, july 1967
c        two corrections by hjortenberg 1974
c     the fast fourier transform in usasi basic fortran
c
c     modified to rc fortran rf june 84
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      dimension datt(idim1),nn(ndim),ifact(32),work(idim2)
c
      np0=0
      nprev=0
c
      twopi=6.283185307
      rthlf=.7071067812
      if(ndim-1)920,1,1
1     ntot=2
      do 2 idim=1,ndim
      if(nn(idim))920,920,2
2     ntot=ntot*nn(idim)
c
c     mainloop for each dimension
c
      np1=2
      do 910 idim=1,ndim
      n=nn(idim)
      np2=np1*n
      if(n-1)920,900,5
c
c     is n a power of two and if not, what are its factors
c
5     m=n
      ntwo=np1
      if=1
      idiv=2
10    iquot=m/idiv
      irem=m-idiv*iquot
      if(iquot-idiv)50,11,11
11    if(irem)20,12,20
12    ntwo=ntwo+ntwo
      ifact(if)=idiv
      if=if+1
      m=iquot
      go to 10
20    idiv=3
      inon2=if
30    iquot=m/idiv
      irem=m-idiv*iquot
      if(iquot-idiv)60,31,31
31    if(irem)40,32,40
32    ifact(if)=idiv
      if=if+1
      m=iquot
      go to 30
40    idiv=idiv+2
      go to 30
50    inon2=if
      if(irem)60,51,60
51    ntwo=ntwo+ntwo
      go to 70
60    ifact(if)=m
70    non2p=np2/ntwo
c
c     separate four cases---
c        1. complex transform
c        2. real transform for the 2nd, 3nd, etc. dimension.  method--
c           transform half the datt, supplying the other half by con-
c           jugate symmetry.
c        3. real transform for the 1st dimension,n odd.  method--
c           set the imaginary parts to zero
c        4. real transform for the 1st dimension,n even.method--
c           transform a complex array of lenght n/2 whose real parts
c           are the even numberd real values and whose imaginary parts
c           are the odd numberedreal values.  separate and supply
c           the second half by conjugate summetry.
c
      icase=1
      ifmin=1
      i1rng=np1
      if(idim-4)74,100,100
74    if(iform)71,71,100
71    icase=2
      i1rng=np0*(1+nprev/2)
      if(idim-1)72,72,100
72    icase=3
      i1rng=np1
      if(ntwo-np1)100,100,73
73    icase=4
      ifmin=2
      ntwo=ntwo/2
      n=n/2
      np2=np2/2
      ntot=ntot/2
      i=1
      do 80 j=1,ntot
      datt(j)=datt(i)
80    i=i+2
c
c     shuffle datt by bit reversal, since n=2**k.  as the shuffling
c     can be done by simple interchange, no working array is needed
c
100   if(non2p-1)101,101,200
101   np2hf=np2/2
      j=1
      do 150 i2=1,np2,np1
      if(j-i2)121,130,130
121   i1max=i2+np1-2
      do 125 i1=i2,i1max,2
      do 125 i3=i1,ntot,np2
      j3=j+i3-i2
      tempr=datt(i3)
      tempi=datt(i3+1)
      datt(i3)=datt(j3)
      datt(i3+1)=datt(j3+1)
      datt(j3)=tempr
125   datt(j3+1)=tempi
130   m=np2hf
140   if(j-m)150,150,141
141   j=j-m
      m=m/2
      if(m-np1)150,140,140
150   j=j+m
      go to 300
c
c     shuffle datt by digit reversal for general n
c
200   nwork=2*n
      do 270 i1=1,np1,2
      do 270 i3=i1,ntot,np2
      j=i3
      do 260 i=1,nwork,2
      if(icase-3)210,220,210
210   work(i)=datt(j)
      work(i+1)=datt(j+1)
      go to 240
220   work(i)=datt(j)
      work(i+1)=0.
240   ifp2=np2
      if=ifmin
250   ifp1=ifp2/ifact(if)
      j=j+ifp1
      if(j-i3-ifp2)260,255,255
255   j=j-ifp2
      ifp2=ifp1
      if=if+1
      if(ifp2-np1)260,260,250
260   continue
      i2max=i3+np2-np1
      i=1
      do 270 i2=i3,i2max,np1
      datt(i2)=work(i)
      datt(i2+1)=work(i+1)
270   i=i+2
c
c     main loop for factors of two
c     w=exp(isign*2*pi*sqrt(-1)*m/(4*mmax)).  check for w=isign*sqrt(-1)
c     and repeat for w=w*(1+isign*sqrt(-1))/sqrt(2)
c
300   if(ntwo-np1)600,600,305
305   np1tw=np1+np1
      ipar=ntwo/np1
310   if(ipar-2)350,330,320
320   ipar=ipar/4
      go to 310
330   do 340 i1=1,i1rng,2
      do 340 k1=i1,ntot,np1tw
      k2=k1+np1
      tempr=datt(k2)
      tempi=datt(k2+1)
      datt(k2)=datt(k1)-tempr
      datt(k2+1)=datt(k1+1)-tempi
      datt(k1)=datt(k1)+tempr
340   datt(k1+1)=datt(k1+1)+tempi
350   mmax=np1
360   if(mmax-ntwo/2)370,600,600
370   lmax=max0(np1tw,mmax/2)
      do 570 l=np1,lmax,np1tw
      m=l
      if(mmax-np1)420,420,380
380   theta=-twopi*float(l)/float(4*mmax)
      if(isign)400,390,390
390   theta=-theta
400   wr=cos(theta)
      wi=sin(theta)
410   w2r=wr*wr-wi*wi
      w2i=2.*wr*wi
      w3r=w2r*wr-w2i*wi
      w3i=w2r*wi+w2i*wr
420   do 530 i1=1,i1rng,2
      kmin=i1+ipar*m
      if(mmax-np1)430,430,440
430   kmin=i1
440   kdif=ipar*mmax
450   kstep=4*kdif
      if(kstep-ntwo)460,460,530
460   do 520 k1=kmin,ntot,kstep
      k2=k1+kdif
      k3=k2+kdif
      k4=k3+kdif
      if(mmax-np1)470,470,480
470   u1r=datt(k1)+datt(k2)
      u1i=datt(k1+1)+datt(k2+1)
      u2r=datt(k3)+datt(k4)
      u2i=datt(k3+1)+datt(k4+1)
      u3r=datt(k1)-datt(k2)
      u3i=datt(k1+1)-datt(k2+1)
      if(isign)471,472,472
471   u4r=datt(k3+1)-datt(k4+1)
      u4i=datt(k4)-datt(k3)
      go to 510
472   u4r=datt(k4+1)-datt(k3+1)
      u4i=datt(k3)-datt(k4)
      go to 510
480   t2r=w2r*datt(k2)-w2i*datt(k2+1)
      t2i=w2r*datt(k2+1)+w2i*datt(k2)
      t3r=wr*datt(k3)-wi*datt(k3+1)
      t3i=wr*datt(k3+1)+wi*datt(k3)
      t4r=w3r*datt(k4)-w3i*datt(k4+1)
      t4i=w3r*datt(k4+1)+w3i*datt(k4)
      u1r=datt(k1)+t2r
      u1i=datt(k1+1)+t2i
      u2r=t3r+t4r
      u2i=t3i+t4i
      u3r=datt(k1)-t2r
      u3i=datt(k1+1)-t2i
      if(isign)490,500,500
490   u4r=t3i-t4i
      u4i=t4r-t3r
      go to 510
500   u4r=t4i-t3i
      u4i=t3r-t4r
510   datt(k1)=u1r+u2r
      datt(k1+1)=u1i+u2i
      datt(k2)=u3r+u4r
      datt(k2+1)=u3i+u4i
      datt(k3)=u1r-u2r
      datt(k3+1)=u1i-u2i
      datt(k4)=u3r-u4r
520   datt(k4+1)=u3i-u4i
      kdif=kstep
      kmin=4*(kmin-i1)+i1
      go to 450
530   continue
      m=m+lmax
      if(m-mmax)540,540,570
540   if(isign)550,560,560
550   tempr=wr
      wr=(wr+wi)*rthlf
      wi=(wi-tempr)*rthlf
      go to 410
560   tempr=wr
      wr=(wr-wi)*rthlf
      wi=(tempr+wi)*rthlf
      go to 410
570   continue
      ipar=3-ipar
      mmax=mmax+mmax
      go to 360
c
c     main loop for factoers not equal to two
c     w=exp(isign*2*pi*sqrt(-1)*(j1+j2-i3-1)/ifp2)
c
600   if(non2p-1)700,700,601
601   ifp1=ntwo
      if=inon2
610   ifp2=ifact(if)*ifp1
      theta=-twopi/float(ifact(if))
      if(isign)612,611,611
611   theta=-theta
612   wstpr=cos(theta)
      wstpi=sin(theta)
      do 650 j1=1,ifp1,np1
      thetm=-twopi*float(j1-1)/float(ifp2)
      if(isign)614,613,613
613   thetm=-thetm
614   wminr=cos(thetm)
      wmini=sin(thetm)
      i1max=j1+i1rng-2
      do 650 i1=j1,i1max,2
      do 650 i3=i1,ntot,np2
      i=1
      wr=wminr
      wi=wmini
      j2max=i3+ifp2-ifp1
      do 640 j2=i3,j2max,ifp1
      twowr=wr+wr
      j3max=j2+np2-ifp2
      do 630 j3=j2,j3max,ifp2
      jmin=j3-j2+i3
      j=jmin+ifp2-ifp1
      sr=datt(j)
      si=datt(j+1)
      oldsr=0.
      oldsi=0.
      j=j-ifp1
620   stmpr=sr
      stmpi=si
      sr=twowr*sr-oldsr+datt(j)
      si=twowr*si-oldsi+datt(j+1)
      oldsr=stmpr
      oldsi=stmpi
      j=j-ifp1
      if(j-jmin)621,621,620
621   work(i)=wr*sr-wi*si-oldsr+datt(j)
      work(i+1)=wi*sr+wr*si-oldsi+datt(j+1)
630   i=i+2
      wtemp=wr*wstpi
      wr=wr*wstpr-wi*wstpi
640   wi=wi*wstpr+wtemp
      i=1
      do 650 j2=i3,j2max,ifp1
      j3max=j2+np2-ifp2
      do 650 j3=j2,j3max,ifp2
      datt(j3)=work(i)
      datt(j3+1)=work(i+1)
650   i=i+2
      if=if+1
      ifp1=ifp2
      if(ifp1-np2)610,700,700
c
c     complete areal transform in the 1st dimension, n even, by con-
c     jugate symmetries
c
700   go to (900,800,900,701),icase
701   nhalf=n
      n=n+n
      theta=-twopi/float(n)
      if(isign)703,702,702
702   theta=-theta
703   wstpr=cos(theta)
      wstpi=sin(theta)
      wr=wstpr
      wi=wstpi
      imin=3
      jmin=2*nhalf-1
      go to 725
710   j=jmin
      do 720 i=imin,ntot,np2
      sumr=(datt(i)+datt(j))/2.
      sumi=(datt(i+1)+datt(j+1))/2.
      difr=(datt(i)-datt(j))/2.
      difi=(datt(i+1)-datt(j+1))/2.
      tempr=wr*sumi+wi*difr
      tempi=wi*sumi-wr*difr
      datt(i)=sumr+tempr
      datt(i+1)=difi+tempi
      datt(j)=sumr-tempr
      datt(j+1)=-difi+tempi
720   j=j+np2
      imin=imin+2
      jmin=jmin-2
      wtemp=wr*wstpi
      wr=wr*wstpr-wi*wstpi
      wi=wi*wstpr+wtemp
725   if(imin-jmin)710,730,740
730   if(isign)731,740,740
731   do 735 i=imin,ntot,np2
735   datt(i+1)=-datt(i+1)
740   np2=np2+np2
      ntot=ntot+ntot
      j=ntot+1
      imax=ntot/2+1
745   imin=imax-2*nhalf
      i=imin
      go to 755
750   datt(j)=datt(i)
      datt(j+1)=-datt(i+1)
755   i=i+2
      j=j-2
      if(i-imax)750,760,760
760   datt(j)=datt(imin)-datt(imin+1)
      datt(j+1)=0.
      if(i-j)770,780,780
765   datt(j)=datt(i)
      datt(j+1)=datt(i+1)
770   i=i-2
      j=j-2
      if(i-imin)775,775,765
775   datt(j)=datt(imin)+datt(imin+1)
c  dummy statement for cdc compiler error
      if (i.eq.-9999) write(*,*) -9999
      datt(j+1)=0.
      imax=imin
      go to 745
780   datt(1)=datt(1)+datt(2)
      datt(2)=0.
      go to 900
c
c     complete a real transform for the 2nd, 3rd, etc. dimension by
c     conjugate symmetries.
c
800   if(i1rng-np1)805,900,900
805   do 860 i3=1,ntot,np2
      i2max=i3+np2-np1
      do 860 i2=i3,i2max,np1
      imax=i2+np1-2
      imin=i2+i1rng
      jmax=2*i3+np1-imin
      if(i2-i3)820,820,810
810   jmax=jmax+np2
820   if(idim-2)850,850,830
830   j=jmax+np0
      do 840 i=imin,imax,2
      datt(i)=datt(j)
      datt(i+1)=-datt(j+1)
840   j=j-2
850   j=jmax
      do 860 i=imin,imax,np0
      datt(i)=datt(j)
      datt(i+1)=-datt(j+1)
860   j=j-np0
c
c     end of loop on each dimension
c
900   np0=np1
      np1=np2
910   nprev=n
920   return
      end

