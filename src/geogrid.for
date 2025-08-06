      program geogrid
      implicit double precision(a-h,o-z)
c $Id: geogrid.for 243 2008-10-29 10:10:19Z cct $
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                        g e o g r i d
c
c  program for gridding of irregular distributed data into
c  a regular rectangular grid, for interpolationg data in profiles,
c  or for interpolating individual points using enhanced, fast
c  weighted means interpolation or collocation/kriging.
c  the program may detrend data prior to gridding, or just fit
c  a trend surface to data
c
c  a quadrant search method is used internally in the program
c  to speed up computations. this means that at each prediction
c  point only the 'nqmax' nearest points is used in each of the
c  four quadrants around the point. to speed up the quadrant
c  search an internal data organization is used, where an
c  internal data grid ensuring 'rdat' (p.t. 3) points per compartment
c  in average is used. this grid is not related to the possible
c  wanted prediction grid.
c
c  program input:
c
c  <ifile>
c  <ofile>
c  <efile>
c  nd, idno, nqmax, itrend, ipred,
c  <predpar>,
c  mode, rkm,
c  <pred points>
c
c  where
c
c  nd     ...  data values in line pr. point in <ifile>
c              not used when <ifile> contain a grid (mode 7, see below)
c
c  idno   ...  data value used out of the 'nd' possibilities.
c              if 'idno' is zero the height data are gridded.
c              if 'idno' is negative then the data number abs(idno)
c              is assumed to be followed by a standard deviation.
c              the standard deviation is used in the collocation
c              prediction only if it is greater than the specified
c              sigma, see predpar below.
c              Standard deviations on heights is specified as idno = -99     
c              Special bouguer/free-air option: if idno = 99,
c              free-air is used at sea, bouguer on land, with free-air
c              anomaly position first on line. sea points have negative h.
c
c  nqmax  ...  number of closest points per quadrant used in prediction,
c              if nqmax<=0 all points are used.
c
c  itrend ...  trend surface removal: 0  none (just grid data) 
c                                     1  remove mean
c                                     2  remove linear function in x and y
c                                     3  remove 2nd order polynomium
c                                     4  remove 3rd order polynomial
c                                     5  remove 4-par (geographic only), 
c                                        corresponds to 7-parameter datum shift
c              if 'itrend' is negative only the residuals from the trending
c              are output (i.e. data trend is not restored after prediction)
c
c  ipred  ...  prediction type:  0  none (detrend only)
c                                1  collocation (kriging)
c                                2  weighted means
c
c  <predpar>   depends on prediction type:
c              ipred = 1: xhalf(km), rms(noise)
c                -     2: prediction power (e.g. 2)
c              note: for collocation rms(noise) is the smallest standard
c              deviation used when data is given with sigmas. 
c
c  mode   ...  specifies data to be predicted:
c              mode = 1: grid in geographical coordinates
c                     2, 3: grid in utm projection
c                     4: profile interpolation
c                     5: point interpolation, data points given in
c                        <efile> in standard gi format, i.e. as
c                        gi-id, lat, lon, height, ...
c                     6: do, with input data assumed to be bouguer
c                        anomalies (grs67/2.67). output data will be
c                        gravity values, derived from interpolated
c                        anomalies and prediction point heights.
c                     7: fill-in of unknown (9999) grid values. In this case
c                        <ifile> must contain a grid rather than a point list.
c                        Only the 9999-values are actually predicted.
c              negative: Same as positive value, except only positive
c                        prediction values allowed (all negative values  
c                        set to zero). Use for DTM interpolation only.
c
c  rkm    ...  margin (km) for data selection area surrounding the
c              wanted grid (or rectangular envelope of the profile
c              or wanted points for mode > 3).
c
c  <pred points>  specifies prediction points depending on 'mode':
c
c  mode 1: fi1, fi2, la1, la2, dfi, dla        (geographic grid)
c  -    2: n1, n2, e1, e2, dn, de / iell, zone (utm grid)
c  -    3: fi1, la1, dn, de, nn, ne / iell, zone (utm grid with
c                                               sw-corner lat/lon,
c                                               spacing, and number of
c                                               points)
c  -    4: fi1, la1, fi2, la2, n               (profile with 'n' points)
c  -    5, 6, 7: (nothing)
c
c  the utm projections are defined by their zone no ('uzone') and
c  ellipsoid no ('ie'): ie = 1: wgs84/nad83, ie = 2: hayf/ed50,
c  ie = 3: clarke/nad27.
c
c  the <ifile> contains data points in standard gi-format, i.e. given
c  as lines 'stat-id, lat, lon, height, data(1), .. , data(nd)', where
c  'idno' specifies which data to be used.
c  data values <= -9999 or >= 9999 signals missing data and not used.
c
c  output on grid format is written in 'ofile'. for collocation errors
c  are written on 'ofile'. for weighted mean interpolation 'efile'
c  contains distances in km to the closest data point.
c  for profile and point interpolation all output is on 'ofile'.
c
c  for collocation a second order markov model is used, with c0
c  assigned from the data variance.
c
c  current program limits:
c  maximal number of points: 60000
c  max e-w prediction grid dimension: 1000
c  max selected pr quadrant: 100
c
c  programmer: rene forsberg, danish geodetic institute, nov 85
c  latest updated:  oct 29, 1991, rf/se
c  modified for vax, rf unsw, dec 88
c  inclusion of individual variances and detrending, rf jan 1990
c  grid fill-in included aug 1990, rf
c  negative modes dec 16, 1991, rf, "spar" block changed July 94 
C  by cct. 
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      dimension dl(36), sa(22)
      dimension pred(2500), sdev(2500), isel(60000), row(1000)
      dimension c(5150), cp(100), ov(11), tcoef(11)
      real*4 z,zsig2
      common /dat/ x(250000),y(250000),z(250000),zsig2(250000)
      common /dsort/ ij(250000),ijn(250000),ifc(900),ih(900)
c     common /spar/ xx,yy,nqmax,nyy,nxx,dyy,dxx,yy1,xx1
      common /spar/ xx,yy,dyy,dxx,yy1,xx1,nqmax,nyy,nxx 
      common /tpar/ cosfi,radeg,xx0,yy0
      logical lutm, lcont, lsig, lfaba, lsel, lrst, ligrid, lneg
      character*128 ifile, ofile, efile
c
c  constants
c
      npmax = 250000
      namax = 3600
      irwdim = 2500
      radeg = 180.0/3.14159265
      degkm = 1.852*60
      rdat = 3.0
      eps = (0.001/degkm)**2
c
      write(*,1)
1     format(/,
     .' *********************************************************',
     .'************'/
     .' *    GEOGRID - GRAVSOFT data gridding - vers. SEP 92 - (c',
     .') RF/KMS   *'/
     .' *********************************************************',
     .'************')
      WRITE(6,*)' INPUT NAMES OF IFILE, OFILE, EFILE'
      read(*,2) ifile
      read(*,2) ofile
      read(*,2) efile
2     format(a128)
c
      open(10,file=ifile,status='old')
      open(20,file=ofile,status='unknown')
c
c  input remaining control data
c  ----------------------------
c
      WRITE(6,*)
     *' INPUT NDATA, DATA NO, MAX USED FOR PRED, ITREND, METH 1 OR 2'
      read(*,*) nd,idno,nqmax,itrend,ipred
c
      lrst = (itrend.gt.0)
      itrend = abs(itrend)
      lsig = (idno.lt.0)
      if (idno.eq.-99) idno = 0
      if (lsig) idno = abs(idno)
      lfaba = (idno.eq.99)
      if (lfaba) idno = 1
c
      IF (IPRED.EQ.1) THEN
      WRITE(6,*)' INPUT CORR. LENGTH (KM) AND VAR. NOISE'
      if (ipred.eq.1) read(*,*) xhalf,rmsn
      ELSE
      WRITE(6,*)' INPUT PRED POWER'
      if (ipred.eq.2) read(*,*) ipwr
      END IF
      WRITE(6,*)
     *' INPUT MODE (1-8) AND SEL. RADIUS (KM)'
      read(*,*) mode,rkm
c
      if (mode.eq.5) then
        open(21,file=efile,status='old')
      else
        open(21,file=efile,status='unknown')
      endif
c
      lneg = (mode.lt.0)
      if (lneg) mode = iabs(mode)
      if (mode.ge.1.and.mode.le.7) goto 10
        write(*, 11)
11      format(' *** mode parameter undefined ')
        stop
      if (itrend.gt.5) stop 'itrend parameter too big'
      if (itrend.eq.5.and.(mode.eq.2.or.mode.eq.3.or.mode.eq.7)) 
     .stop 'itrend=3 only allowed for usual geographic grid'
10    lutm = .false.
      if (mode.ge.2.and.mode.le.3) lutm = .true.
      ligrid = (mode.eq.7)
      goto (12,12,13,14,15,15,18),mode
c
c  grid limits. udeg converts to internal program unit: degree
c  mode 1 and 2: grid limit specs
c
12    udeg = 1.0 
      write(*,*)' Input grid specification '
      read(*,*) rfi1,rfi2,rla1,rla2,dfi,dla
      if (lutm) read(*,*) iell,izone
      if (lutm) udeg = 1.0/(degkm*1000)
      rfi1 = rfi1*udeg
      rfi2 = rfi2*udeg
      dfi = dfi*udeg
      rla1 = rla1*udeg
      rla2 = rla2*udeg
      dla = dla*udeg
      nn = (rfi2-rfi1)/dfi+1.499999
      ne = (rla2-rla1)/dla+1.499999
      write(*,*)nn,ne 
      if (lutm) call utmcon(iell, izone, sa)
      goto 20
c
c  mode 3: utm grid spec by sw corner lat/lon
c
13    udeg = 1.0/(degkm*1000)
      read(*,*) rfi,rla,dfi,dla,nn,ne
      read(*,*) iell,izone
      rfi = rfi/radeg
      rla = rla/radeg
      dfi = dfi*udeg
      dla = dla*udeg
      call utmcon(iell, izone, sa)
      call utg(rfi, rla, rn, re, sa, .false., .true.)
      rfi1 = rn*udeg
      rla1 = re*udeg
      rfi2 = rfi1 + (nn-1)*dfi
      rla2 = rla1 + (ne-1)*dla
      goto 20
c
c  mode 4: profile interpolation
c
14    read(*,*) rfi1,rla1,rfi2,rla2,nn
      dfi = .00001
      dla = .00001
      if (nn.gt.1) dfi = (rfi2-rfi1)/(nn-1)
      if (nn.gt.1) dla = (rla2-rla1)/(nn-1)
      goto 20
c
c  mode 5 or 6: point data. find min and max lat and lon.
c
15    rfi1 = 999999.9
      rla1 = 999999.9
      rfi2 = -999999.9
      rla2 = -999999.9
      nn = 0
16    if (mode.eq.5) read(21,*,end=17) id,rlat,rlon,rh
      if (mode.eq.6) read(21,*,end=17) id,rlat,rlon,rh,rba
      nn = nn + 1
      if (rlat.lt.rfi1) rfi1 = rlat
      if (rlat.gt.rfi2) rfi2 = rlat
      if (rlon.lt.rla1) rla1 = rlon
      if (rlon.gt.rla2) rla2 = rlon
      goto 16
17    rewind(21)
      goto 20
c
c  mode 7: grid input
c
18    read(10,*) rfi1,rfi2,rla1,rla2,dfi,dla
      udeg = 1.0 
      lutm = (abs(rfi1).gt.100.or.abs(rfi2).gt.100)
      if (lutm) then 
        read(10,*) iell,izone
        call utmcon(iell, izone, sa)
        udeg = 1.0/(degkm*1000)
      endif
      if (.not.lutm) write(*,1801) rfi1,rfi2,rla1,rla2,dfi,dla
      if (lutm) then
        write(*,1802) rfi1,rfi2,rla1,rla2,dfi,dla,iell,izone
        rfi1 = rfi1*udeg
        rfi2 = rfi2*udeg
        rla1 = rla1*udeg
        rla2 = rla2*udeg
        dfi = dfi*udeg
        dla = dla*udeg
      endif
1801  format(/' input grid: ',6f9.4)
1802  format(/' input grid: ',6f9.0,2i4)
      nn = (rfi2-rfi1)/dfi + 1.5
      ne = (rla2-rla1)/dla + 1.5
      if (ne.gt.irwdim) stop 'too long rows in grid'
c
c  parameter checks
c
20    if (lutm) then
        if (iell.lt.1.or.iell.gt.4.or.
     .  izone.lt.1.or.izone.gt.60)
     .  stop 'utm grid label not ok'
      endif
      lsel = .true.
      if (nqmax.le.0) lsel = .false.
      if (nqmax.gt.100) nqmax = 100
      if (mode.gt.3) goto 22
      if (nn.ge.1.and.ne.ge.1.and.ne.le.irwdim) goto 22
        write(*, 21)
21      format(' *** grid specification illegal')
        stop
22    cosfi = 1.0
      if (.not.lutm) cosfi = cos((rfi1+rfi2)/2/radeg)
      racos = radeg*cosfi
      rdeg = rkm/degkm
      if (rdeg.le.0) rdeg = 0.0001
      if (mode.eq.4) dpr = sqrt(dfi**2+(dla*cosfi)**2)*degkm
c
c   input list of data from ifile
c   -----------------------------
c
      ndp = 0
      np = 0
      nneg = 0
      y1 = min(rfi1,rfi2) - rdeg
      y2 = max(rfi1,rfi2) + rdeg
      x1 = min(rla1,rla2)*cosfi - rdeg
      x2 = max(rla1,rla2)*cosfi + rdeg
      zmin = 9999999.9
      zmax = -9999999.9
      zsmin = zmin
      zsmax = zmax
      zsum = 0
      zsum2 = 0
      xx1 = 9999999.9
      xx2 = -9999999.9
      yy1 = 9999999.9
      yy2 = -9999999.9
      ii = nn+1
      jj = ne+1
c
c  loop entry point for data input
c  -------------------------------
c
25    if (ligrid) then
24      if (jj.eq.ne+1) then 
          if (ii.eq.1) goto 30
          read(10,*) (row(kk),kk=1,ne)
          jj = 1
          ii = ii-1
        endif 
        yy = rfi1 + (ii-1)*dfi
        xx = (rla1 + (jj-1)*dla)*cosfi
        zz = row(jj)
        jj = jj+1
        ndp = ndp + 1
        if (zz.ge.9999) goto 24
      else
        read(10,*,end=30) id,rlat,rlon,rh,(dl(kk+5),kk=1,nd)
        ndp = ndp + 1
        if (.not.lutm) then
          yy = rlat
          xx = rlon*cosfi
        else
          call utg(rlat/radeg,rlon/radeg,rn,re,sa,.false.,.true.)
          yy = rn*udeg
          xx = re*udeg
        endif
        zz = dl(5+idno)
        if (lfaba.and.rh.gt.0.and.rh.lt.9999) zz = dl(5+idno+1)
        if (idno.eq.0) zz = rh
        zzsig = dl(5+idno+1)
c
        if (xx.le.x1.or.xx.ge.x2.or.yy.le.y1.or.yy.ge.y2) goto 25
        if (zz.le.-9999.or.zz.ge.9999) goto 25
      endif
c
      np = np+1
      if (xx.lt.xx1) xx1 = xx
      if (xx.gt.xx2) xx2 = xx
      if (yy.lt.yy1) yy1 = yy
      if (yy.gt.yy2) yy2 = yy
      if (np.le.npmax) goto 29
        write(*, 28) npmax
28      format(' *** too many data points - max: ',i5)
        stop
29    x(np) = xx
      y(np) = yy
      z(np) = zz
      if (lsig) zsig2(np) = zzsig**2
      if (zz.lt.zmin) zmin = zz
      if (zz.gt.zmax) zmax = zz
      if (lsig.and.zzsig.lt.zsmin) zsmin = zzsig
      if (lsig.and.zzsig.gt.zsmax) zsmax = zzsig
      zsum = zz+zsum
      zsum2 = zz**2 +zsum2
c      write(*,*) ' test np = ',np,id
      goto 25
c
c  detrend data if required
c  ------------------------
c
30    if (itrend.eq.0) goto 309
      xx0 = xx1
      yy0 = yy1
      do 301 i = 1, 66
301   c(i) = 0.0
      do 302 k = 1, np
        zz = z(k)
        call setov(x(k),y(k),zz,itrend,nt,ov)
        ii = 0
        do 302 i = 1, nt+1
        do 302 j = 1, i
          ii = ii+1
          c(ii) = c(ii) + ov(i)*ov(j)
302   continue          
      call chol(c,nt,nsing)
      if (nsing.ne.0) write(*,303) nsing
303   format(' *** warning: detrending ',i2,' singularities')
      k = nt*(nt+1)/2
      do 304 i = 1, nt
304   tcoef(i) = c(k+i)
c
      tsum = 0
      tsum2 = 0
      tmin = 999999.9
      tmax = -tmin
      do 305 k = 1, np 
        rr = z(k) - trend(x(k),y(k),itrend,tcoef)
        tsum = rr + tsum
        tsum2 = rr**2 + tsum2
        if (rr.lt.tmin) tmin = rr
        if (rr.gt.tmax) tmax = rr
        z(k) = rr
305   continue
c
c  output info on input, data and selected internal grid
c  -----------------------------------------------------
c
309   if (np.gt.1) goto 31
        write(*, 32)
32      format(' *** too few points in selection area')
        stop
31    rs = sqrt((zsum2 - zsum**2/np)/(np-1))
      rm = zsum/np
      if (itrend.gt.0) then
        ts = sqrt((tsum2 - tsum**2/np)/(np-1))
        tm = tsum/np
      endif
c
c  covariance parameters - second order markov model
c
      if (ipred.ne.1) goto 79
      sqrc0 = rs
      if (itrend.gt.0) sqrc0 = ts
      c0 = sqrc0**2
      cnn = rmsn**2
      alfa = 0.595*xhalf/degkm
c
79    continue
c
      if (.not.lutm) write(*,85) rkm,y1,y2,x1/cosfi,x2/cosfi
85    format(' rkm = ',f6.1,', data selection area: ',4f9.4)
      if (lutm) write(*,86) rkm,y1/udeg,y2/udeg,x1/udeg,x2/udeg
86    format(' rkm = ',f6.1,', data selection area: ',4f10.0)
      if (.not.lutm) goto 89
      call utg(rfi2/udeg,rla1/udeg,yy,xx,sa,.true.,.true.)
      call utg(rfi2/udeg,rla2/udeg,yi,xi,sa,.true.,.true.)
      call utg(rfi1/udeg,rla1/udeg,yj,xj,sa,.true.,.true.)
      call utg(rfi1/udeg,rla2/udeg,yk,xk,sa,.true.,.true.)
      write(*,88) yy*radeg,xx*radeg,yi*radeg,xi*radeg,
     *yj*radeg,xj*radeg,yk*radeg,xk*radeg
88    format(' geographical coordinates of grid corners:',
     .2(/,'   ',2f9.3,'   ',2f9.3))
89    continue
c
      if (ligrid) then
        ndp = nn*ne
        nd = 1
        idno = 1
      endif
      write(*,87) nd,idno,ndp,np,zmin,zmax,rm,rs
87    format(/' data values per point: ',i2,', used no.: ',i2,/,
     *' total points in file:',i6,', selected:',i5,/,
     *' min max mean stddev: ',4f10.3)
      if (lsig) write(*,871) zsmin, zsmax
871   format(' minimal and maximal standard deviations of data: ',2f9.2)
c
      if (itrend.gt.0) write(*,801) itrend,nt,(tcoef(j),j=1,nt)
      if (itrend.gt.0) write(*,802) tmin,tmax,tm,ts
801   format(/' detrending done on data, itrend = ',i1,/
     .' no of trend parameters estimated: ',i2,/
     .' solution: ',5e12.4,/,' ',5e12.4)
802   format(' detrended data (min,max,mean,stddev): ',4f9.3)
c
      if (ipred.eq.1) write(*,81) sqrc0, xhalf, rmsn
81    format(/' collocation prediction - sqrc0,xhalf(km),rmsn = ',3f7.2)
      if (ipred.eq.2) write(*,82) ipwr
82    format(/' weighted means prediction, ipwr =',i2)
      if (.not.lsel) write(*,83)
83    format(' all points used in each prediction')
      if (lsel) write(*,84) nqmax
84    format(' selection:',i3,' closest points per quadrant')
c
c  skip data organization if all points used
c
      if (.not.lsel) goto 100
c
33    r = sqrt((yy2-yy1)*(xx2-xx1)/np*rdat)
      if (r.eq.0) r = 1.0
      nyy = (yy2-yy1)/r+.4999
      nxx = (xx2-xx1)/r+.4999
      if (nyy.lt.1) nyy = 1
      if (nxx.lt.1) nxx = 1
      if (nxx*nyy.le.namax) goto 35
        write(*, 34)
34      format(' - large organization, rdat increased -')
        rdat = 2.0*rdat
        goto 33
35    continue
      dyy = (yy2-yy1)/nyy
      dxx = (xx2-xx1)/nxx
      if (dyy.eq.0) dyy = .00001
      if (dxx.eq.0) dxx = .00001
      nyy1 = nyy+1
      nxx1 = nxx1+1
      n = nxx*nyy
c
      if (.not.lutm) write(*,37) yy1,yy2,xx1/cosfi,xx2/cosfi
37    format(/' data organization limits: ',4f9.4)
      if (lutm) write(*,38) yy1/udeg,yy2/udeg,xx1/udeg,xx2/udeg
38    format(/' data organization limits: ',4f10.0)
      write(*,39) nyy,nxx,n,dyy*degkm,dxx*degkm,np*1.0/n
39    format(' subrectangles (n,e,total): ',3i5,
     */' size (km): ',2f9.1,', average pts per rect (rdat):',f9.3)
c
c  organization of irregular data
c  ------------------------------
c
c  the following piece of code corresponds closely to subroutine
c  'oaf' of 'gspp' by hans sunkel, see osu internal report: gspp
c  manual, ohio state university, 1980
c
c  ij(1...np)   stores element index
c  ijn(1...np)  stores data index in increasing sort element order
c  ifc(1...n1)  stores number of data per sort element
c  ih(1...n1)   pointer vector, data in sort elem k has index ijn(ih(k))
c
      do 41 i = 1,n
41    ifc(i) = 0
      do 45 i = 1,np
        iy = (y(i)-yy1)/dyy
        ix = (x(i)-xx1)/dxx
        if (iy.ge.nyy) iy = nyy-1
        if (ix.ge.nxx) ix = nxx-1
        ixy = iy*nxx + ix+1
        ij(i) = ixy
        ifc(ixy) = ifc(ixy) + 1
45    continue
c
      ih(1) = 1
      do 46 i = 2,n
        i1 = i-1
        ih(i) = ifc(i1)+ih(i1)
46    continue
c
      do 47 i = 1,np
        k = ij(i)
        kk = ih(k)
        ijn(kk) = i
        ih(k) = ih(k) + 1
47    continue
c
      j = 0
      jj = 0
      do 48 i = 1,n
        k = ifc(i)
        ih(i) = ih(i)-k
        if (k.eq.0) j=j+1
        if (k.gt.jj) jj=k
48    continue
c
      write(*, 49) jj, 100.0*j/n
49    format(' max points in subrects:',i5,
     *', percentage with no points: ',f5.1)
c
c  prediction loop
c  ---------------
c
c  the point numbers to be used in prediction is stored in
c  array 'isel', selected by subroutine 'select'.
c
100   zmin = 9999999.9
      zmax = -9999999.9
      zsum = 0
      zsum2 = 0
      if (mode.ge.4.and.mode.le.6) ne = 1
      if (lsel) goto 102
        nsel = np
        do 101 i = 1, nsel
101     isel(i) = i
c
102   if (ligrid) then
        rewind(10)
        read(10,*) r1,r2,r3,r4,r5,r6
        if (lutm) read(10,*) iell,izone
      endif
c
      do 240 ii = nn, 1, -1
      if (ligrid) read(10,*) (row(kk),kk=1,ne)
      do 200 jj = 1, ne
        goto (105,105,105,106,107,107,105),mode
c
c  grid points
c
105     xx = ((jj-1)*dla+rla1)*cosfi
        yy = (ii-1)*dfi+rfi1
        if (ligrid) then
          if (row(jj).lt.9999) then
            zpred = row(jj)
            zstd = 0
            goto 190
          endif  
        endif
        goto 108
c
c  profile points
c
106     xx = ((nn-ii)*dla+rla1)*cosfi
        yy = (nn-ii)*dfi+rfi1
        goto 108
c
c  individual points
c
107     if (mode.eq.5) read(21,*) id,rlat,rlon,rh
        if (mode.eq.6) read(21,*) id,rlat,rlon,rh,rba
        yy = rlat
        xx = rlon*cosfi
c
c  skip predictions for ipred = 0, otherwise select closest data
c
108     if (ipred.eq.0) then
          zpred = 0.0
          zstd = 0.0
          goto 190
        endif
        if (lsel) call select(isel,nsel)
c
        if (ipred.eq.1) goto 150
c
c  prediction by weighted means
c  ----------------------------
c
        wsum = 0
        wfsum = 0
        r2min = 9999999.9
        do 120 i = 1, nsel
          k = isel(i)
          dx = x(k)-xx
          dy = y(k)-yy
          r2 = dx**2 + dy**2
          if (r2.lt.r2min) r2min = r2
          if (r2.gt.eps) goto 110
            zpred = z(k)
            zstd = 0.0
            goto 190
c
c  distance weights
c
110       if (ipwr.ne.2) goto 112
            w = 1/r2
            goto 115
112       w = 1/sqrt(r2)**ipwr
115       wsum = w + wsum
          wfsum = w*z(k) + wfsum
120     continue
        zpred = wfsum/wsum
        zstd = sqrt(r2min)*degkm
        goto 190
c
c  collocation prediction
c  ----------------------
c
150     k = 1
        do 156 i = 1,nsel+1
          if (i.le.nsel) goto 153
          xi = xx
          yi = yy
          i1 = i-1
          goto 154
153       xi = x(isel(i))
          yi = y(isel(i))
          i1 = i
154       do 155 j = 1,i1
            kk = isel(j)
            r = sqrt((x(kk)-xi)**2 + (y(kk)-yi)**2)
c
c  second order markov covariance function
c
            cov = c0*(1 + r/alfa)*exp(-r/alfa)
            if (i.eq.j) cov = cov + max(cnn,dble(zsig2(kk)))
c
               if (k.le.0) write(*,*)k
            if (k.gt.0) c(k) = cov
            k = k+1
          if (i.gt.nsel) cp(j) = cov
155       continue
156     continue
c
        call chol(c,nsel,nsing)
        if (nsing.gt.0) write(*,157) ii,jj,nsing
157     format(' *** cholesky: ',i3,' singularities point: ',2i5)
c
        zpred = 0
        zstd = 0
        i1 = nsel*(nsel+1)/2
        do 158 i = 1,nsel
          r = c(i+i1)
          zpred = zpred + r*z(isel(i))
          zstd = zstd + r*cp(i)
158     continue
        r = c0 - zstd
        if (r.lt.0) zstd = -sqrt(-r)
        if (r.ge.0) zstd = sqrt(r)
c
c  statistics for prediction
c  -------------------------
c
190     if (lrst) zpred = zpred + trend(xx,yy,itrend,tcoef)
        if (lneg) then
          if (zpred.lt.0) zpred = 0
          nneg = nneg+1
        endif
        if (zpred.lt.zmin) zmin = zpred
        if (zpred.gt.zmax) zmax = zpred
        zsum = zsum + zpred
        zsum2 = zsum2 + zpred**2
        pred(jj) = zpred
        sdev(jj) = zstd
200   continue
c
c  output result and possible heading
c  ----------------------------------
c
      mode8 = mode
      if (mode.eq.7.and.lutm) mode8 = 8
      goto (201,211,211,221,231,231,201,211),mode8
c
c  output of geographical grid data
c
201   if (ii.eq.nn) write(20,204) rfi1,rfi1+(nn-1)*dfi,
     *rla1,rla1+(ne-1)*dla,dfi,dla
      write(20,205) (pred(jj),jj=1,ne)
204   format(' ',4f12.6,2f12.8)
205   format(/,50(' ',7f10.3,/))
c
      if (ii.eq.nn.and.ipred.eq.1) write(21,207)
     *rfi1,rfi1+(nn-1)*dfi,rla1,rla1+(ne-1)*dla,dfi,dla
      if (ii.eq.nn.and.ipred.eq.2) write(21,208)
     *rfi1,rfi1+(nn-1)*dfi,rla1,rla1+(ne-1)*dla,dfi,dla
207   format(' ',4f12.6,2f12.7)
208   format(' ',4f12.6,2f12.7)
      if (ipred.eq.1) write(21,205) (sdev(jj),jj=1,ne)
      if (ipred.eq.2) write(21,206) (sdev(jj),jj=1,ne)
206   format(/,40(' ',8f9.1,/))
      goto 240
c
c  utm grid data
c
211   yi = rfi1/udeg
      yj = (rfi1+(nn-1)*dfi)/udeg
      xi = rla1/udeg
      xj = (rla1+(ne-1)*dla)/udeg
      yy = dfi/udeg
      xx = dla/udeg
      if (ii.lt.nn) goto 217
      write(20,212) yi,yj,xi,xj,yy,xx,iell,izone
212   format(' ',6f11.0,/,' ',2i4)
217   write(20,205) (pred(jj),jj=1,ne)
c
      if (ii.lt.nn) goto 220
      write(21,212) yi,yj,xi,xj,yy,xx,iell,izone
220   if (ipred.eq.1) write(21,205) (sdev(jj),jj=1,ne)
      if (ipred.eq.2) write(21,206) (sdev(jj),jj=1,ne)
      goto 240
c
c  profile data
c
221   if (ii.eq.nn) write(*,224) rfi1,rla1,rfi2,rla2,dpr
224   format(/' profile ',2f8.3,'  to ',2f8.3,', spacing ',f9.2,
     *' km')
      write(20,223) (nn-ii)*dpr, zpred, zstd
223   format(' ', f9.2,2f12.3)
      goto 240
c
c  point data
c
231   dl(6) = zpred
      dl(7) = zstd
      if (mode.ne.6) goto 232
c
c  conversion grs67 bouguer anomaly to gravity value
c
      r = sin(rlat/radeg)**2
      gamma67 = 978031.85 * (1 + .005278895*r +
     *                       .000023462*r**2)
      if (rh.ge.0) hfakt = 0.1967
      if (rh.lt.0) hfakt = -0.0687
      dl(6) = zpred + gamma67 - hfakt*rh
c
232   write(20,233) id,rlat,rlon,rh,dl(6),dl(7)
233   format(' ',i9,2f12.5,f9.2,2f9.2)
c
240   continue
      n = nn*ne
      if (n.eq.1) rs = 0.0
      if (n.gt.1) rs = sqrt((zsum2 - zsum**2/n)/(n-1))
      rm = 0.0
      if (n.gt.0) rm = zsum/n
      write(*,250) n,zmin,zmax,rm,rs
250   format(' '/' predicted: ',i7,' points'/
     *' prediction  min max mean std.dev.: ',4f9.2/)
      if (.not.lrst.and.itrend.gt.0) write(*, 251)
251   format(' - trend in data not restored in output -'/)
      if (lneg) write(*,252) nneg
252   format(' - number of negative predictions set to zero:',i7)
c
      end
c
      subroutine select(isel,nsel)
      implicit double precision(a-h,o-z)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                       s e l e c t
c
c  subroutine for selection of 'nqmax' points in each quadrant
c
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension isel(60000)
      real*4 z,zsig2
      common /dat/ x(60000),y(60000),z(60000),zsig2(60000)
      common /dsort/ ij(60000),ijn(60000),ifc(900),ih(900)
c     common /spar/ xx,yy,nqmax,nyy,nxx,dyy,dxx,yy1,xx1
      common /spar/ xx,yy,dyy,dxx,yy1,xx1,nqmax,nyy,nxx 
      common /qpar/ rq2(4,100),iq(4,100),ns(4),r2max(4),iimax(4)
      logical lq1,lq2,lq3,lq4
c
c  initialization of parameters for rect elem sort
c
      lq1 = .false.
      lq2 = .false.
      lq3 = .false.
      lq4 = .false.
      do 11 j = 1,4
        r2max(j) = 0.0
        ns(j) = 0
11    continue
c
      i0 = (yy-yy1)/dyy+1
      j0 = (xx-xx1)/dxx+1
c
c  smallest scan area is 3 x 3
c
      i1 = i0-1
      i2 = i0+1
      j1 = j0-1
      j2 = j0+1
      if (i1.lt.1) i1=1
      if (i2.gt.nyy) i2=nyy
      if (j1.lt.1) j1=1
      if (j2.gt.nxx) j2=nxx
      do 8 i = i1,i2
      do 8 j = j1,j2
8     call rsort(i,j)
c
c  check section
c  check quadrants one by one and scan
c
10    if (ns(1).ge.nqmax) lq1 = .true.
      if (ns(2).ge.nqmax) lq2 = .true.
      if (ns(3).ge.nqmax) lq3 = .true.
      if (ns(4).ge.nqmax) lq4 = .true.
      i1 = i1-1
      i2 = i2+1
      j1 = j1-1
      j2 = j2+1
      if (i2.gt.nyy.and.j1.lt.1) lq1 = .true.
      if (i2.gt.nyy.and.j2.gt.nxx) lq2 = .true.
      if (i1.lt.1.and.j1.lt.1) lq3 = .true.
      if (i1.lt.1.and.j2.gt.nxx) lq4 = .true.
c
      if (lq1.and.lq2.and.lq3.and.lq4) goto 90
c
c  scan edges
c
      if (.not.(lq1.and.lq2)) call rsort(i2,j0)
      if (.not.(lq2.and.lq4)) call rsort(i0,j2)
      if (.not.(lq3.and.lq4)) call rsort(i1,j0)
      if (.not.(lq1.and.lq3)) call rsort(i0,j1)
      if (lq1) goto 30
        do 25 i = i0+1,i2
25      call rsort(i,j1)
        do 26 j = j1+1,j0-1
26      call rsort(i2,j)
30    if (lq2) goto 40
        do 35 i = i0+1,i2
35      call rsort(i,j2)
        do 36 j = j0+1,j2-1
36      call rsort(i2,j)
40    if (lq3) goto 50
        do 45 i = i1,i0-1
45      call rsort(i,j1)
        do 46 j = j1+1,j0-1
46      call rsort(i1,j)
50    if (lq4) goto 60
        do 55 i = i1,i0-1
55      call rsort(i,j2)
        do 56 j = j0+1,j2-1
56      call rsort(i1,j)
c
60    goto 10
c
90    nsel = 0
      k = 0
      do 91 i = 1,4
        k = ns(i)
        do 92 j = 1,k
92      isel(j+nsel) = iq(i,j)
        nsel = nsel + k
91    continue
      return
      end
c
      subroutine rsort(i,j)
      implicit double precision(a-h,o-z)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                       r s o r t
c
c  subroutine sorts 'nqmax' closest points within each quadrant
c  in a subrectangle
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*4 z,zsig2
      common /dat/ x(60000),y(60000),z(60000),zsig2(60000)
      common /dsort/ ij(60000),ijn(60000),ifc(900),ih(900)
c     common /spar/ xx,yy,nqmax,nyy,nxx,dyy,dxx,yy1,xx1
      common /spar/ xx,yy,dyy,dxx,yy1,xx1,nqmax,nyy,nxx 
      common /qpar/ rq2(4,100),iq(4,100),ns(4),r2max(4),iimax(4)
c
      if (i.lt.1.or.i.gt.nyy) return
      if (j.lt.1.or.j.gt.nxx) return
c
      ixy = (i-1)*nxx+j
      k = ifc(ixy)
      if (k.eq.0) return
c
      k1 = ih(ixy)
      k2 = k1+k-1
      do 50 k = k1,k2
        kk = ijn(k)
        dx = x(kk)-xx
        dy = y(kk)-yy
        if (dy.lt.0) goto 15
        iiq = 1
        if (dx.ge.0) iiq = 2
        goto 20
15      iiq = 3
        if (dx.ge.0) iiq = 4
20      r2 = dx**2+dy**2
c
c  add point if not enough in quadrant
c
        ii = ns(iiq)+1
        if (ii.gt.nqmax) goto 30
        ns(iiq) = ii
        iq(iiq,ii) = kk
        rq2(iiq,ii) = r2
        if (r2.lt.r2max(iiq)) goto 50
        r2max(iiq) = r2
        iimax(iiq) = ii
        goto 50
c
c  check if point is among the closest points
c
30      if (r2.gt.r2max(iiq)) goto 50
        ii = iimax(iiq)
        iq(iiq,ii) = kk
        rq2(iiq,ii) = r2
        do 31 jj = 1,nqmax
          rr = rq2(iiq,jj)
          if (rr.le.r2) goto 31
          r2 = rr
          iimax(iiq) = jj
31      continue
        r2max(iiq) = r2
c
50    continue
      end
c
      subroutine setov(xx,yy,obs,itrend,nt,ov)
      implicit double precision(a-h,o-z)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      s e t o t v
c
c  sets observation equation coefficients for detrending
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension ov(11)
      common /tpar/ cosfi,radeg,xx0,yy0
      if (itrend.eq.5) goto 50
c
c  just mean value or polynomials
c
      dx = xx-xx0
      dy = yy-yy0
      nt = 1  
      ov(1) = 1.0
      if (itrend.eq.1) goto 90
      nt = 3
      ov(2) = dx
      ov(3) = dy
      if (itrend.eq.2) goto 90
      nt = 6
      ov(4) = dx**2
      ov(5) = dy**2
      ov(6) = dx*dy
      if (itrend.eq.3) goto 90
      nt = 10
      ov(7) = dx**3
      ov(8) = dy**3
      ov(9) = dx**2*dy
      ov(10) = dx*dy**2
      goto 90
c
c  4-parameter transformation corresponding to 7-par spatial transformation
c
50    nt = 4
      rfi = yy/radeg
      rla = xx/cosfi/radeg
      cosf = cos(rfi)
      ov(1) = cosf*cos(rla)
      ov(2) = cosf*sin(rla)
      ov(3) = sin(rfi)
      ov(4) = 1.0
c
90    ov(nt+1) = obs
      return
      end
c
      function trend(xx,yy,itrend,tsol)
      implicit double precision(a-h,o-z)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                         t r e n d
c
c  finds value of trend function for given point
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension tsol(11),ov(11)
      rr = 0.0
      call setov(xx,yy,rr,itrend,nt,ov)
      do 10 i = 1, nt
10    rr = rr + ov(i)*tsol(i)
      trend = rr
      return
      end
c
      subroutine chol(c,n,nsing)
      implicit double precision(a-h,o-z)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      c h o l
c
c  subroutine solves positive definite symmetric linear equations
c  using cholesky decomposition. coefficients stored columnwise
c  in c, followed by righthand side. n is number of unknowns.
c  solution is returned as last column in c.
c  'nsing' is number of singularities. it should be zero for a
c  succesful solution.
c
c  rf, nov 85
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension c(5150)
      nsing = 0
c
      do 50 nr = 1,n+1
        i=nr*(nr-1)/2
        ir=i
        do 40 nc = 1,nr
          if (nc.gt.n) goto 50
          sum=0
          ic=nc*(nc-1)/2
          i=i+1
          nc1=nc-1
          do 30 np = 1,nc1
30        sum = sum - c(ir+np)*c(ic+np)
          ci = c(i)+sum
          if (nr.eq.nc) goto 31
          cc = c(ic+nc)
          if (cc.eq.0) then
            write(*,*) '*** Cholesky reduction zero-division'
            cc = 1.0e9
          endif
          c(i) = ci/cc
          goto 40
31        if (nr.gt.n) goto 40
          if (ci.gt.0) goto 32
            nsing = nsing+1
            c(i) = 1.0e9
            goto 40
32        c(i) = sqrt(ci)
40      continue
50    continue
c
c  back substitution
c
      do 80 nc = n,1,-1
        ir=i
        ic=nc*(nc+1)/2
        cc = c(ic)
        if (cc.eq.0) then
          write(*,*) '*** Cholesky subsitution zero-division'
          cc = 1.0e9
        endif
        c(i) = c(i)/cc
        do 70 np = nc-1,1,-1
          ir=ir-1
          ic=ic-1
          c(ir) = c(ir) - c(i)*c(ic)
70      continue
        i = i-1
80    continue
      return
      end
c
      subroutine utmcon(isys, izone, sa)
      implicit double precision(a-h,o-z)
      dimension sa(22)
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                         u t m c o n
c
c the procedure produces the constants needed in the transfor-
c mations between transversal mercator and geographical coordina-
c tes and the tolerances needed for the check of the transfor-
c mations. the transformation constants are generated from a
c reg_label defining a transversal mercator system. the formulae
c are taken from k\nig und weise : mathematische grundlagen der
c h\heren geod<sie und kartographie, erster band, berlin 1951.
c
c parameters
c __________
c
c isys, izone         (call)            integers
c specifies ellipsoid and utm zone. the following ellipsoids are
c currently implemented:
c
c     1: wgs84,  2: hayford (ed50),  3: clarke (nad27)
c
c sa                  (return)          array
c the constants needed in the transformation.
c sa(1) =           normalized meridian quadrant  (geotype),
c sa(2) =           easting at the central meridian (geotype),
c sa(3) =           longitude of the central meridian (geotype),
c sa(4)  -  sa(7) = const for ell. geo -> sph. geo (real)
c sa(8)  - sa(11) = const for sph. geo -> ell. geo (real),
c sa(12) - sa(15) = const for sph. n, e -> ell. n, e (real),
c sa(16) - sa(19) = const for ell. n, e -> sph. n, e (real),
c sa(20) =          toler. for utm input, 5 mm.
c sa(21) =          toler. for geo input, do.
c sa(22) =          not used in fortran version, which operates
c                   with  n e g a t i v e  northings s of equator.
c the user may change sa(20) - sa(21) for special checks.
c
c prog: knud poder, danish geodetic institute, 7 nov 1977,
c updated 18 sep 1983;
c rc fortran version alp/rf oct 86
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision n,m
      radeg = 180/3.1415926536
c
c  set ellipsoid parameters
c
      goto (10,20,30),isys
c
c  wgs84 ellipsoid
c
10    a = 6378137.0
      f = 1/298.2572236
      goto 40
c
c  hayford ed 50 ellipsoid
c
20    a = 6378388.0
      f = 1/297.0
      goto 40
c
c  clarke nad 27 ellipsoid
c
30    a = 6378206.4
      f = 1/294.9786982
c
40    eastpr = 500000.0
      dm = 4.0e-4
c
c  normalized meridian quadrant
c  see k\nig und weise p.50 (96), p.19 (38b), p.5 (2)
c
      n=f/(2.0-f)
      m=n**2*(1.0/4.0+n**2/64.0)
      w= a*(-n - dm+m*(1.0-dm))/(1.0+n)
      sa(1)=a + w
c
c  central easting and longitude
c
      sa(2)=eastpr
      sa(3)=((izone - 30)*6 - 3)/radeg
c
c  check-tol for transformation
c  5.0 mm on earth
c
      sa(20) = 0.0050
      sa(21) = sa(20)/a
c
c  coef of trig series
c
c  ell. geo -> sph. geo., kw p186 - 187 (51) - (52)
c
      sa(4) = n*(-2 + n*(2.0/3.0 + n*(4.0/3.0 + n*(-82.0/45.0))))
      sa(5) = n**2*(5.0/3.0 + n*(-16.0/15.0 + n*(-13.0/9.0)))
      sa(6) = n**3*(-26.0/15.0 + n*34.0/21.0)
      sa(7) = n**4*1237.0/630.0
c
c   sph. geo - ell. geo., kw p190 - 191 (61) - (62)
c
      sa(8)=n*(2.0+n*(-2.0/3.0 +n*(-2.0+n*116.0/45.0)))
      sa(9)=n**2*(7.0/3.0+n*(-8.0/5.0+n*(-227.0/45.0)))
      sa(10)=n**3*(56.0/15.0+n*(-136.0)/35.0)
      sa(11)=n**4* (4279.0/630.0)
c
c  sph. n, e -> ell. n, e,  kw p196 (69)
c
      sa(12)=n*(1.0/2.0+n*(-2.0/3.0+
     .                            n*(5.0/16.0+n*41.0/180.0)))
      sa(13)=n**2*(13.0/48.0+
     .                            n*(-3.0/5.0+n*557.0/1440.0))
      sa(14)=n**3*(61.0/240.0+n*(-103.0/140.0))
      sa(15)=n**4*(49561.0/161280.0)
c
c  ell. n, e -> sph. n, e,  kw p194 (65)
c
      sa(16)=n*(-1.0/2.0+n*(2.0/3.0+
     .                           n*(-37.0/96.0+n*1.0/360.0)))
      sa(17)=n**2*(-1.0/48.0+
     .                           n*(-1.0/15.0+n*437.0/1440.0))
      sa(18)=n**3* (-17.0/480.0+n*37.0/840.0)
      sa(19)=n**4*(-4397.0/161280.0)
c
      return
      end
c
      subroutine utg(rn, re, b, l, sa, direct, tcheck)
      implicit double precision(a-h,o-z)
      double precision b, l, sa(22)
      logical direct, tcheck
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                       u t g
c
c dual autochecking transformation procedure for transformation
c between geographical coordinates and transversal mercator co-
c ordinates. the procedure transforms utm->geo when direct is
c true and the reverse when direct is false.
c an alarm is produced when the check by the inverse transforma-
c tion exceeds the tolerance of 5.0 mm or an other value set by
c the user in sa(19) for utm->geo or sa(20) for geo->utm.
c
c n, e              (call)             real
c the utm- or geographical coordinates input for trans-
c formation in meters or radians.
c
c b, l              (return)           real
c the geographical or utm-coordinates output from the procedure
c as transformed and checked coordinates in radians or meters
c
c sa                (call)             array
c transformation constants for direct and inverse transf.
c see fields in sa or set_utm_const for a description
c
c direct            (call)             logical
c direct = true => transformation utm -> geogr.
c direct = false => transformation geogr -> utm
c
c tcheck            (call)             logical
c tcheck = true => check by back transformation
c tcheck = false => no check. this possibility doubles the
c                   speed of the subroutine, with the risk of
c                   obtaining bad results at large (> 60 deg)
c                   distances from the central meridian.
c
c programmer: knud poder, dgi, nov 1977
c rcfortran alp/rf oct 86
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision np, lg, l0, ndif, ncheck
      integer h, s
      double precision n, e
c
      n = rn
      e = re
      qn = sa(1)
      e0 = sa(2)
      l0 = sa(3)
c
c  transformation sequence
c
      if (direct) i=1
      if (.not.direct) i=3
      h = 4-i
      s = 2-i
c
c  check-values
c
      ncheck=n
      echeck=e
c
c  transformation cases
c
      do 100 i=i,h,s
      goto (10,20,30),i
c
c  case 1, utm -> geo
c  ------------------
c
10    np = n/qn
      ep = (e - e0)/qn
c
c  ellip. n, e -> sph. n, e
c
      np = np + clcsin(sa, 15, 4, 2.0*np, 2.0*ep, dn, de)
      ep = ep + de
c
c  sph. n, e = compl. sph. lat -> sph lat, lng
c
      cosbn = cos(np)
      if (cosbn.eq.0) cosbn = 1.0d-33
      snh = (exp(ep) - exp(-ep))/2
      lg   = atan(snh/cosbn)
      bbg  = atan(sin(np)*cos(lg)/cosbn)
c
c  sph. lat, lng -> ell. lat, lng
c
      bbg = bbg + clsin(sa, 7, 4, 2.0*bbg)
      lg = lg + l0
      goto 100
c
c  case 2, transf results
c  ----------------------
c
20    b     = bbg
      n     = bbg
      l     = lg
      e     = lg
      if (tcheck) goto 100
      return
c
c  case 3, geo -> utm
c  ------------------
c
30    bbg   = n + clsin(sa, 3, 4, 2.0*n)
      lg    = e - l0
c
c  sph. lat, lng -> compl. sph. lat = sph n, e
c
      cosbn = cos(bbg)
      if (cosbn.eq.0) cosbn = 1.0d-33
      np    = atan(sin(bbg)/(cos(lg)*cosbn))
      rr = sin(lg)*cosbn
      if (abs(rr).ge.0.95) goto 40
      ep = log((1+rr)/(1-rr))/2
      goto 41
40    ep = 1.0e38
41    continue
c
c  sph. normalized n, e -> ell. n, e
c
      np = np + clcsin(sa, 11, 4, 2.0*np, 2.0*ep, dn, de)
      ep = ep + de
      bbg = qn*np
      lg = qn*ep + e0
c
  100 continue
c
c  in/rev-dif for check
c
      ndif  = bbg - ncheck
      edif  = lg - echeck
      edcos = edif
      if(.not.direct) edcos = edcos*cos(ncheck)
c
c  error actions
c
      if (direct) tol = sa(20)
      if (.not.direct) tol = sa(21)
      if (abs(ndif).lt.tol.and.abs(edcos).lt.tol) return
c
      n = rn
      e = re
      if (direct) n = b
      if (direct) e = l
      ep = 6371000.0
      if (direct) ep = 1.0
      write(*, 90) n*57.29578, e*57.29578, ndif*ep, edcos*ep
90    format(' *** utg coor ',2f7.1,' checkdiff too large: ',
     .2f8.3, ' m')
      return
      end
c
      double precision function clsin(a, i0, g, arg)
      implicit double precision(a-h,o-z)
      dimension a(22)
      integer g
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      c l s i n
c
c computes the sum of a series of a(i+i0)*sin(i*arg) by clenshaw
c summation from g down to 1.
c the sum is the value of the function.
c
c prog.: knud poder 1978
c
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      cosarg = 2*cos(arg)
c
      hr1 = 0.0
      hr  = a(g+i0)
c
      do 10 it = g - 1,1,-1
        hr2 = hr1
        hr1 = hr
        hr = -hr2 + cosarg*hr1 + a(it+i0)
   10 continue
c
      clsin = hr*sin(arg)
c
      return
      end
c
      double precision function clcsin(a, i0, g, argr, argi, r, i)
      implicit double precision(a-h,o-z)
      dimension a(22)
      double precision r, i
      integer g
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                        c l c s i n
c
c computes the sum of a series a(i+i0)*sin(i*arg_r + j*i*arg_i)
c (where j is the imaginary unit) by clenshaw summation
c from g down to 1. the coefficients are here real and
c the argument of sin is complex. the real part of the
c sum is the value of the function.
c
c prog.: knud poder 1978
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision ii
c
      sinar = sin(argr)
      cosar = cos(argr)
      ex = exp(argi)
      emx = exp(-argi)
      sinhar = (ex - emx)/2
      coshar = (ex + emx)/2
      rr = 2*cosar*coshar
      ii = -2*sinar*sinhar
c
      hr1 = 0.0
      hi1 = 0.0
      hi  = 0.0
      hr  = a(g+i0)
c
      do 10 it = g-1, 1, -1
        hr2 = hr1
        hr1 = hr
        hi2 = hi1
        hi1 = hi
        hr  = -hr2 + rr*hr1 - ii*hi1 + a(it+i0)
        hi  = -hi2 + ii*hr1 + rr*hi1
   10 continue
c
      rr = sinar*coshar
      ii = cosar*sinhar
c
      r = rr*hr - ii*hi
      clcsin = r
      i = rr*hi + ii*hr
c
      return
      end

