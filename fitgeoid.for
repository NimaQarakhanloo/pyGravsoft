      program fitgeoid
      implicit double precision(a-h,o-z)
c $Id: fitgeoid.for 181 2008-08-14 23:19:51Z tjansson $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                        f i t g e o i d
c
c  this program calls geogrid and geoip to grid differences between
c  a gravimetric geoid and a set of gps geoid observations
c
c  the geoid must be gravsoft grid format, and the gps levelling
c  data in format
c 
c  id, lat, lon, h, Ngps
c
c  program input:
c  <geoid grid file>
c  <N file>
c  <final geoid file>
c  nqmax, itrend, ipred,
c  <predpar>
c  lat1, lat2, lon1, lon2, dlat, dlon    
c 
c  where
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
c
c  ipred  ...  prediction type:  1  collocation (kriging)
c                                2  weighted means
c
c  <predpar>   depends on prediction type:
c              ipred = 1: xhalf(km), rms(noise)
c                -     2: prediction power (e.g. 2)
c              note: for collocation rms(noise) is the smallest standard
c              deviation used when data is given with sigmas. 
c
c  lat1, lat2, lon1, lon2, dlat, dlon    
c              wanted geographic grid spacing (might be different from original geoid)
c
c  Intermediate results are stored in a number of files with fixed names:
c  
c  "fitgeoid_dif.dat"  Differences Ngps - grid before fitting
c  "fitgeoid_dif.gri"  Grid of correction surface Ngps - grid
c  "fitgeoid_dif.err"  Grid of correction surface errors 
c  "fitgeoid_dif2.dat" Differences Ngps - grid after fitting (as a point list)
c 
c  (c) Rene Forsberg, DTU-Space, Aug 2008
c  Based on gravsoft "geogrid" and "geoip"
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension dl(36), sa(22)
      dimension pred(5000), sdev(5000), isel(100), row(5000)
      dimension c(5150), cp(100), ov(11), tcoef(11)
      real*4 z,zsig2,t
      common /dat/x(79000),y(79000),z(79000),zsig2(79000),
     .t(79000)
      common /dsort/ij(79000),ijn(79000),ifc(8100),ih(8100)
      common /spar/xx,yy,dyy,dxx,yy1,xx1,r2min,dmax,
     .thr,tmin,nqmax,nyy,nxx,ld,lt
      common /tpar/cosfi,radeg,xx0,yy0
      common /g2lpar/
     .re,rfi0,rla0,cosfi0,sinfi0,resin,sinfi2,sinfi12,itype
      logical lutm, lsig, lfaba, lsel, lrst, ligrid, lneg, lt, ld, lthr
      character*72 geoidfile, gpsfile, ofile, file1, file2
c
c  constants
c
      npmax = 79000
      namax = 8100
      irwdim = 5000
      radeg = 180/3.1415926535d0
      degkm = 1.852*60
      rdat = 3.0
c
      write(*,1)
1     format(/,
     .' *********************************************************'/
     .' *    FITGEOID - GRAVSOFT geoid fitter - vers. AUG08     *'/
     .' *********************************************************')
      read(*,2) geoidfile
      read(*,2) gpsfile
      read(*,2) ofile
2     format(a72)
c
c  input remaining control data
c  ----------------------------
c
      read(*,*) nqmax,itrend,ipred
      nd = 1
      idno = 1
c
      lrst = (itrend.gt.0)
      itrend = abs(itrend)
      lsig = (idno.lt.0)
      if (idno.eq.-99) idno = 0
      if (lsig) idno = abs(idno)
      lfaba = (idno.eq.99)
      if (lfaba) idno = 1
c
      if (ipred.eq.1) read(*,*) xhalf,rmsn
      if (ipred.eq.2) read(*,*) ipwr
c
c  fitgeoid
c      read(*,*) mode,rkm
      mode = 1
      rkm = 10.0
c
      lneg = (mode.lt.0)
      if (lneg) mode = iabs(mode)
      if (mode.lt.1.or.mode.gt.9) then
        write(*, 11)
11      format(' *** mode parameter undefined ')
        stop
      endif
      if (itrend.gt.5) stop 'itrend parameter too big'
      if (itrend.eq.5.and.(mode.eq.2.or.mode.eq.3.or.mode.eq.9)) 
     .stop 'itrend=3 only allowed for usual geographic grid'
      lutm = .false.
      if (mode.ge.2.and.mode.le.3) lutm = .true.
      ld = (mode.eq.7.or.mode.eq.8)
      lt = (mode.eq.8)
      ligrid = (mode.eq.9)
      lthr = .false.
c
      goto (12,12,13,14,15,15,15,15,18),mode
c
c  grid limits
c  mode 1 and 2: grid limit specs
c
12    read(*,*) rfi1,rfi2,rla1,rla2,dfi,dla
      if (lutm) read(*,*) iell,izone
      nn = (rfi2-rfi1)/dfi+1.499999
      ne = (rla2-rla1)/dla+1.499999
      if (lutm) call utmcon(iell, izone, sa)
      goto 20
c
c  mode 3: utm grid spec by sw corner lat/lon
c
13    read(*,*) rfi,rla,dfi,dla,nn,ne
      read(*,*) iell,izone
      call utmcon(iell, izone, sa)
      call utg(rfi/radeg, rla/radeg, rn, re, sa, .false., .true.)
      rfi1 = rn
      rla1 = re
      rfi2 = rfi1 + (nn-1)*dfi
      rla2 = rla1 + (ne-1)*dla
      goto 20
c
c  mode 4: profile data
c
14    read(*,*) rfi1,rla1,rfi2,rla2,nn
      dfi = .00001
      dla = .00001
      if (nn.gt.1) dfi = (rfi2-rfi1)/(nn-1)
      if (nn.gt.1) dla = (rla2-rla1)/(nn-1)
      goto 20
c
c  mode 5,6,7,8: point data. find min and max lat and lon.
c
15    rfi1 = 999999.9
      rla1 = 999999.9
      rfi2 = -999999.9
      rla2 = -999999.9
      if (mode.eq.5.or.mode.eq.6) read(*,*) lthr
      if (ld) then
        if (mode.eq.7) read(*,*) dmax,lthr
        if (mode.eq.8) read(*,*) dmax,tmin
        dmax = dmax*1000
        dmax2 = dmax**2
      endif
c
      nn = 0
16    if (mode.le.6) then
        read(21,*,end=17) thr,rlat,rlon,rh
        id = thr
      elseif (mode.eq.7) then
        read(21,*,end=17) thr,rlat,rlon,rh,(dl(j),j=1,nd)
        id = thr
      elseif (lt) then
        read(21,*,end=17) thr,rlat,rlon,rh,(dl(j),j=1,nd)
      endif
c
      nn = nn + 1
      if (rlat.lt.rfi1) rfi1 = rlat
      if (rlat.gt.rfi2) rfi2 = rlat
      if (rlon.lt.rla1) rla1 = rlon
      if (rlon.gt.rla2) rla2 = rlon
      goto 16
17    rewind(21)
      goto 20
c
c  mode 9: grid input
c
18    read(10,*) rfi1,rfi2,rla1,rla2,dfi,dla
      lutm = (abs(rfi1).gt.100.or.abs(rfi2).gt.100)
      if (lutm) then 
        read(10,*) iell,izone
        call utmcon(iell, izone, sa)
        write(*,1802) rfi1,rfi2,rla1,rla2,dfi,dla,iell,izone
1802    format(/' input grid: ',6f9.0,2i4)
      else
        write(*,1801) rfi1,rfi2,rla1,rla2,dfi,dla
1801    format(/' input grid: ',6f9.4)
      endif
      nn = (rfi2-rfi1)/dfi + 1.5
      ne = (rla2-rla1)/dla + 1.5
      if (ne.gt.irwdim) stop 'too long rows in grid'
c
c  parameter checks
c
20    if (lutm) then
        if (iell.lt.1.or.iell.gt.4.or.
     .  izone.lt.1.or.izone.gt.99)
     .  stop 'utm grid label not ok'
      endif
      lsel = .true.
      if (nqmax.le.0) lsel = .false.
      if (nqmax.gt.25) stop '*** nqmax too large'
      if (mode.le.3.and.(nn.lt.1.or.ne.lt.1.or.ne.gt.irwdim)) then
        write(*, 21)
21      format(' *** grid specification illegal')
        stop
      endif
c
c  define internal map projection for geographic 
c  make special case for polar regions
c
      if (lutm) then
        cosfi = 1.0
      else 
        if (rfi2.eq.90.and.abs(rla2-rla1).gt.180) then
          rlat0 = 90
          cosfi = cos(rfi1/radeg)
        elseif (rfi1.eq.-90.and.abs(rla2-rla1).gt.180) then
          rlat0 = -90
          cosfi = cos(rfi2/radeg)
        else
          rlat0 = (rfi1+rfi2)/2
          cosfi = cos(rlat0/radeg)
        endif
        rlon0 = (rla1+rla2)/2
        call setlamb(rlat0,rlon0)
        write(*,22) rlat0,rlon0
22      format(' Internal Lambert map projection center: ',2f9.3)
      endif
      if (mode.eq.4) dpr = sqrt(dfi**2+(dla*cosfi)**2)*degkm
c
c   fitgeoid: open files and subtract (gcomb)
c   -----------------------------------------
      write(*,*)
      write(*,*) '=== Subtraction: Ngps minus geoid grid: ==='
c
      file1 = 'fitgeoid_dif.dat'
      mode = 11
      call geoip(geoidfile,gpsfile,file1,mode)
c
c   now open geogrid files
c   -----------------------
c
      write(*,*) 
      write(*,*) '=== Gridding of GPS corrections ==='
      open(10,file='fitgeoid_dif.dat',status='old')
      open(20,file='fitgeoid_dif.gri',status='unknown')
      open(21,file='fitgeoid_dif.err',status='unknown')
c
c   input list of data from ifile (unit 10)
c   ---------------------------------------
c
      ndp = 0
      np = 0
      nneg = 0
      if (lutm) then
        r = rkm*1000 + 0.1 
        rfimin = min(rfi1,rfi2) - r
        rfimax = max(rfi1,rfi2) + r
        rlamin = min(rla1,rla2) - r
        rlamax = max(rla1,rla2) + r
      else
        r = rkm/degkm + 0.000001
        rfimin = min(rfi1,rfi2) - r
        rfimax = max(rfi1,rfi2) + r
        rlamin = min(rla1,rla2) - r/cosfi
        rlamax = max(rla1,rla2) + r/cosfi
      endif
      zmin = 9999999.9
      zmax = -9999999.9
      zsmin = zmin
      zsmax = zmax
      zpsmin = zmin
      zpsmax = zmax
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
        if (lutm) then
          yy = rfi1 + (ii-1)*dfi
          xx = (rla1 + (jj-1)*dla)
        else
          call geo2lamb(.true.,rfi1+(ii-1)*dfi,rla1+(jj-1)*dla,yy,xx)
        endif
        zz = row(jj)
        jj = jj+1
        ndp = ndp + 1
        if (zz.ge.9999) goto 24
      else
c
c  input point reading
c
        read(10,*,end=30) thr,rlat,rlon,rh,(dl(kk),kk=1,nd)
        id = thr
c
        ndp = ndp + 1
        if (.not.lutm) then
          if (rlat.lt.rfimin.or.rlat.gt.rfimax.or.
     .    rlon.lt.rlamin.or.rlon.gt.rlamax) goto 25
          call geo2lamb(.true.,rlat,rlon,yy,xx)
        else
          call utg(rlat/radeg,rlon/radeg,yy,xx,sa,.false.,.true.)
          if (xx.lt.rlamin.or.xx.gt.rlamax.or.
     .    yy.lt.rfimin.or.yy.gt.rfimax) goto 25
        endif
        zz = dl(idno)
        if (zz.le.-9999.or.zz.ge.9999) goto 25
        if (lfaba.and.rh.gt.0.and.rh.lt.9999) zz = dl(idno+1)
        if (idno.eq.0) zz = rh
        zzsig = dl(idno+1)
      endif
c
      np = np+1
      if (xx.lt.xx1) xx1 = xx
      if (xx.gt.xx2) xx2 = xx
      if (yy.lt.yy1) yy1 = yy
      if (yy.gt.yy2) yy2 = yy
      if (np.gt.npmax) then
        write(*, 28) npmax
28      format(' *** too many data points - max: ',i7)
        stop
      endif
      x(np) = xx
      y(np) = yy
      z(np) = zz
      if (lsig) zsig2(np) = zzsig**2
      if (lt) t(np) = thr
      if (zz.lt.zmin) zmin = zz
      if (zz.gt.zmax) zmax = zz
      if (lsig.and.zzsig.lt.zsmin) zsmin = zzsig
      if (lsig.and.zzsig.gt.zsmax) zsmax = zzsig
      zsum = zz+zsum
      zsum2 = zz**2 +zsum2
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
      trmin = 999999.9
      trmax = -trmin
      do 305 k = 1, np 
        rr = z(k) - trend(x(k),y(k),itrend,tcoef)
        tsum = rr + tsum
        tsum2 = rr**2 + tsum2
        if (rr.lt.trmin) trmin = rr
        if (rr.gt.trmax) trmax = rr
        z(k) = rr
305   continue
c
c  output info on input, data and selected internal grid
c  -----------------------------------------------------
c
309   if (np.lt.1) then
        write(*, 32)
32      format(' *** too few points in selection area')
        write(*,'(4f14.4)') rfimin,rfimax,rlamin,rlamax
        stop
      endif
      if (np.gt.1) then
        rs = sqrt((zsum2 - zsum**2/np)/(np-1))
      else 
        rs = 0
      endif
      rm = zsum/np
      if (itrend.gt.0) then
        if (np.gt.1) then
          ts = sqrt((tsum2 - tsum**2/np)/(np-1))
        else
          ts = 0
        endif
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
      alfa = 0.595*xhalf*1000
c
79    continue
c
      if (.not.lutm) then
c        write(*,85) rkm,rfimin,rfimax,rlamin,rlamax
c85      format(' rkm = ',f6.1,', data selection area: ',4f9.4)
      else
        write(*,86) rkm,
     .  rfimin,rfimax,rlamin,rlamax
86      format(' rkm = ',f6.1,', data selection area: ',4f10.0)
        call utg(rfi2,rla1,yy,xx,sa,.true.,.true.)
        call utg(rfi2,rla2,yi,xi,sa,.true.,.true.)
        call utg(rfi1,rla1,yj,xj,sa,.true.,.true.)
        call utg(rfi1,rla2,yk,xk,sa,.true.,.true.)
        write(*,88) yy*radeg,xx*radeg,yi*radeg,xi*radeg,
     *  yj*radeg,xj*radeg,yk*radeg,xk*radeg
88      format(' geographical coordinates of grid corners:',
     .  2(/,'   ',2f9.3,'   ',2f9.3))
      endif
c
      if (ligrid) then
        ndp = nn*ne
        nd = 1
        idno = 1
      endif
c      write(*,87) nd,idno,ndp,np,zmin,zmax,rm,rs
c87    format(/' data values per point: ',i2,', used no.: ',i2,/,
c     *' total points in file:',i7,', selected:',i7,/,
c     *' min max mean stddev: ',4f10.3)
      if (lsig) write(*,871) zsmin, zsmax
871   format(' minimal and maximal standard deviations of data: ',2f9.2)
c
      if (itrend.gt.0) write(*,801) itrend,nt,(tcoef(j),j=1,nt)
      if (itrend.gt.0) write(*,802) trmin,trmax,tm,ts
801   format(' detrending done on data, itrend = ',i1,/
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
      if (lt) write(*,841) tmin
841   format(' - predictions only using data with time dif ',
     .'greater than: ',f8.3)
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
      n = nxx*nyy
c
      if (.not.lutm) write(*,37) yy1,yy2,xx1,xx2
37    format(/' data organization limits in lambert proj:',4f9.0)
      if (lutm) write(*,38) yy1,yy2,xx1,xx2
38    format(/' data organization limits: ',4f10.0)
      write(*,39) nyy,nxx,n,dyy/1000,dxx/1000,np*1.d0/n
39    format(' subrectangles (n,e,total): ',3i5,
     */' size (km): ',2f9.2,', average pts per rect (rdat):',f9.3)
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
        if (dyy.eq.0) then
          iy = 1
        else
          iy = (y(i)-yy1)/dyy
        endif
        if (dxx.eq.0) then
          ix = 1
        else
          ix = (x(i)-xx1)/dxx
        endif
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
      zvalmin = 9999999.9
      zvalmax = -9999999.9
      zvalsum = 0
      zvalsum2 = 0
      zdifmin = 9999999.9
      zdifmax = -9999999.9
      zdifsum = 0
      zdifsum2 = 0
      no = 0
c
      if (mode.ge.4.and.mode.le.8) ne = 1
      if (lsel) goto 102
        nsel = np
        if (nsel.gt.100) stop '** max 100 pts in collocation'
        do 101 i = 1, nsel
101     isel(i) = i
c
102   if (ligrid) then
        rewind(10)
        read(10,*) r,r,r,r,r,r
        if (lutm) read(10,*) iell,izone
      endif
c
c  main loop entry
c
      do 240 ii = nn, 1, -1
      if (ligrid) read(10,*) (row(kk),kk=1,ne)
      do 200 jj = 1, ne
        goto (105,105,105,106,107,107,107,107,105),mode
c
c  grid points
c
105     if (lutm) then
          xx = (jj-1)*dla+rla1
          yy = (ii-1)*dfi+rfi1
        else
          call geo2lamb(.true.,(ii-1)*dfi+rfi1,(jj-1)*dla+rla1,yy,xx)
        endif
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
106     if (lutm) then
          xx = ((nn-ii)*dla+rla1)
          yy = (nn-ii)*dfi+rfi1
        else
          call geo2lamb(.true.,(nn-ii)*dfi+rfi1,(nn-ii)*dla+rla1,yy,xx)
        endif
        goto 108
c
c  individual points
c
107     if (mode.le.6) then
          read(21,*) thr,rlat,rlon,rh
          id = thr
        elseif (mode.ge.7) then
          read(21,*) thr,rlat,rlon,rh,(dl(j),j=1,nd)
          id = thr
          if (nd.eq.0) then
            zval = rh
          else
            zval = dl(idno)
          endif 
        endif
c
        if (lutm) then
          yy = rlat
          xx = rlon
        else
          call geo2lamb(.true.,rlat,rlon,yy,xx)
        endif
c
c  skip predictions for ipred = 0, otherwise select closest data
c
108     if (ipred.eq.0) then
          zpred = 0.0
          zstd = 0.0
          goto 190
        endif
        if (lsel) call select(isel,nsel)
        if (ld) then
          if (r2min.gt.dmax2) goto 240
        endif
c
        if (ipred.eq.1) goto 150
c
c  prediction by weighted means
c  ----------------------------
c
        wsum = 0
        wfsum = 0
        do 120 i = 1, nsel
          k = isel(i)
          dx = x(k)-xx
          dy = y(k)-yy
          r2 = dx**2 + dy**2
c
c  non-zero value if dist lt 10 cm to avoid singularity
c
          if (r2.lt.0.01) r2 = 0.01     
c
c  distance weights
c
          if (ipwr.eq.2) then
            w = 1/r2
          else          
            w = 1/sqrt(r2)**ipwr
          endif
          wsum = w + wsum
          wfsum = w*z(k) + wfsum
120     continue
        zpred = wfsum/wsum
        zstd = sqrt(r2min)/1000
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
            c(k) = cov
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
        no = no+1
        if (zpred.lt.zmin) zmin = zpred
        if (zpred.gt.zmax) zmax = zpred
        zsum = zsum + zpred
        zsum2 = zsum2 + zpred**2
        if (zstd.lt.zpsmin) zpsmin = zstd
        if (zstd.gt.zpsmax) zpsmax = zstd
        pred(jj) = zpred
        sdev(jj) = zstd
200   continue
c
c  output result and possible heading
c  ----------------------------------
c
      mode10 = mode
      if (mode.eq.9.and.lutm) mode10 = 10
      goto (201,211,211,221,231,231,231,231,201,211),mode10
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
206   format(/,40(' ',8f9.3,/))
      goto 240
c
c  utm grid data
c
211   yi = rfi1
      yj = (rfi1+(nn-1)*dfi)
      xi = rla1
      xj = (rla1+(ne-1)*dla)
      yy = dfi
      xx = dla
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
c  conversion grs80 bouguer anomaly to gravity value
c
231   if (mode.eq.6) then
        r = sin(rlat/radeg)**2
        gamma80 = 978032.677d0*(1+.00193185135d0*r)/
     .  sqrt(1 - .00669438002d0*r)
        if (rh.ge.0) hfakt = 0.1967
        if (rh.lt.0) hfakt = -0.0687
        zpred = zpred + gamma80 - hfakt*rh
      endif
      if (ld) then
c
        zdif = zval - zpred 
c
        if (zval.lt.zvalmin) zvalmin = zval
        if (zval.gt.zvalmax) zvalmax = zval
        zvalsum = zvalsum + zval
        zvalsum2 = zvalsum2 + zval**2
        if (zdif.lt.zdifmin) zdifmin = zdif
        if (zdif.gt.zdifmax) zdifmax = zdif
        zdifsum = zdifsum + zdif
        zdifsum2 = zdifsum2 + zdif**2
c
        if (lt.or.lthr) then
          write(20,234) thr,rlat,rlon,rh,zpred,zstd,zdif
        else
          write(20,233) id,rlat,rlon,rh,zpred,zstd,zdif
        endif
      else
        if (lthr) then        
          write(20,234) thr,rlat,rlon,rh,zpred,zstd
        else
          write(20,233) id,rlat,rlon,rh,zpred,zstd
        endif
      endif
233   format(' ',i9,2f11.5,f9.2,f11.2,2f9.2)
234   format(' ',f12.5,2f11.5,f9.2,f11.3,2f9.3)
c
240   continue
c ------------------------------------------------
c  prediction loop exit
c
      write(*,*)
      if (ld) then
        if (no.le.1) rs = 0.0
        if (no.gt.1) rs = sqrt((zvalsum2 - zvalsum**2/no)/(no-1))
        rm = 0
        if (no.gt.0) rm = zsum/no
        write(*,249) nn,dmax/1000,rm,rs,zvalmin,zvalmax
249     format(' total number of points in datafile (efile): ',i7,/ 
     .  ' prediction pts selected closer than ',f7.3,' km'/
     .  ' selected pts   mean std.dev. min max: ',4f9.2)
      endif

      if (no.le.1) rs = 0.0
      if (no.gt.1) rs = sqrt((zsum2 - zsum**2/no)/(no-1))
      rm = 0
      if (no.gt.0) rm = zsum/no
      write(*,250) no,rm,rs,zmin,zmax,zpsmin,zpsmax
250   format(' predicted: ',i7,' points'/
     *' prediction pts mean std.dev. min max: ',4f9.2/,
     .' prediction error values min max:   ',21x,2f9.2)
c
      if (ld) then
        if (no.le.1) rs = 0.0
        if (no.gt.1) rs = sqrt((zdifsum2 - zdifsum**2/no)/(no-1))
        rm = 0
        if (no.gt.0) rm = zdifsum/no
        write(*,251) rm,rs,zdifmin,zdifmax
251     format(' differences    mean std.dev. min max: ',4f9.2)
      endif
c
      if (.not.lrst.and.itrend.gt.0) write(*, 255)
255   format(' - trend in data not restored in output -'/)
      if (lneg) write(*,256) nneg
256   format(' - number of negative predictions set to zero:',i7)
c
      close(10)
      close(20)
      close(21)
c
c  fitgeoid: now add gridded corrections to original geoid grid
c  ------------------------------------------------------------
c
      write(*,*)
      write(*,*) '=== Addition of gridded corrections to geoid file ==='
      file1 = 'fitgeoid_dif.gri'
      mode = 16
      call geoip(file1,geoidfile,ofile,mode)
c
c  geoidfit: now check fit of geoid 
c  --------------------------------
      write(*,*)
      write(*,*) '=== Fit of GPS levelling data to fitted geoid ==='
      write(*,*) 'Name of fitted geoid: ',ofile
      file2 = 'fitgeoid_dif2.dat'
      mode = 11
      call geoip(ofile,gpsfile,file2,mode)
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
      dimension isel(100)
      real*4 z,zsig2,t
      common /dat/x(79000),y(79000),z(79000),zsig2(79000),
     .t(79000)
      common /dsort/ij(79000),ijn(79000),ifc(8100),ih(8100)
      common /spar/xx,yy,dyy,dxx,yy1,xx1,r2min,dmax,
     .thr,tmin,nqmax,nyy,nxx,ld,lt
      common /qpar/rq2(4,25),iq(4,25),ns(4),r2max(4),iimax(4)
      logical lq1,lq2,lq3,lq4,ld,lt
c
c  initialization of parameters for rect elem sort
c
      lq1 = .false.
      lq2 = .false.
      lq3 = .false.
      lq4 = .false.
      r2min = 9.999d9
      do 11 j = 1,4
        r2max(j) = 0.0
        ns(j) = 0
11    continue
c
      if (dyy.eq.0) then
        i0 = 1
      else
        i0 = (yy-yy1)/dyy+1
      endif
      if (dxx.eq.0) then
        j0 = 1
      else
        j0 = (xx-xx1)/dxx+1
      endif
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
c  check if points in new cells farther away than dmax
c
      if (ld) then
        if ((i0-i1-1)*dyy.gt.dmax.and.(j0-j1-1)*dxx.gt.dmax) goto 90
      endif
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
      real*4 z,zsig2,t
      common /dat/x(79000),y(79000),z(79000),zsig2(79000),
     .t(79000)
      common /dsort/ij(79000),ijn(79000),ifc(8100),ih(8100)
      common /spar/xx,yy,dyy,dxx,yy1,xx1,r2min,dmax,
     .thr,tmin,nqmax,nyy,nxx,ld,lt
      common /qpar/rq2(4,25),iq(4,25),ns(4),r2max(4),iimax(4)
      logical ld,lt
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
        if (lt) then
          if (abs(t(kk)-thr).lt.tmin) goto 50
        endif
        if (dy.lt.0) goto 15
        iiq = 1
        if (dx.ge.0) iiq = 2
        goto 20
15      iiq = 3
        if (dx.ge.0) iiq = 4
20      r2 = dx**2+dy**2
        if (r2.lt.r2min) r2min = r2
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
      common /tpar/cosfi,radeg,xx0,yy0
      common /g2lpar/
     .re,rfi0,rla0,cosfi0,sinfi0,resin,sinfi2,sinfi12,itype
c
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
      call geo2lamb(.false.,yy,xx,rfi,rla)
      cosf = cos(rfi/radeg)
      ov(1) = cosf*cos(rla/radeg)
      ov(2) = cosf*sin(rla/radeg)
      ov(3) = sin(rfi/radeg)
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
c
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
      implicit double precision (a-h,o-z)
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                         u t m c o n
c
c the procedure produces the constants needed in the transfor-
c mations between transversal mercator and geographical coordina-
c tes and the tolerances needed for the check of the transfor-
c mations. the transformation constants are generated from a
c reg_label defining a transversal mercator system. the formulae
c are taken from konig und weise : mathematische grundlagen der
c hoheren geodasie und kartographie, erster band, berlin 1951.
c
c parameters
c __________
c
c isys, izone         (call)            integers
c specifies ellipsoid and utm zone. the following ellipsoids are
c currently implemented:
c
c     1: wgs84,  2: hayford (ed50),  3: clarke (nad27)
c     4: bessel, if utmzone = 99 then bessel and national swedish
c                projection is assumed (nb: this version is only
c                good to approximatively 10 m for sweden, meridian
c                exact longitude unknown)   
c     5, 6: irish and british systems (utmzone 97 and 98)
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
c sa(22) =          north offset 
c                  
c the user may change sa(20) - sa(21) for special checks.
c
c prog: knud poder, danish geodetic institute, 7 nov 1977,
c updated 18 sep 1983;
c rc fortran version alp/rf oct 86, last updated dec 1990
c updated nov 2001/jan 2002 qatar and british systems
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension sa(22)
      double precision n,m,northpr
      dimension aa(6),ff(6)
c
c  ellipsoid parameters
c                  WGS84      Hayford      Clarke       Bessel
c                  Mod.Airy      Airy
      data aa / 6378137.0d0, 6378388.0d0, 6378206.4d0, 6377397.155d0,
     .          6377340.189d0, 6377563.396d0/
      data ff / 298.2572236d0, 297.0d0, 294.9786982d0, 299.153d0,
     .          299.3250d0,  299.3250d0/
      radeg = 180/3.1415926536d0
c
      if (isys.lt.1.or.isys.gt.6) stop 'illegal ellipsoid in utmcon'
      if (izone.lt.1.or.izone.gt.99) stop 'illegal UTM zone in utmcon'
c
      eastpr = 500000.d0
      northpr = 0.d0
      dm = 4.0d-4
      if (izone.eq.95) then
c  qbc
        dm = 1 - 0.999996d0
        eastpr =  400000.d0
        northpr = -2766043.105d0+500000d0
        clong = (50+49/60.d0)/radeg
      elseif (izone.eq.96) then
c  qng
        dm = 1 - 0.99999d0
        eastpr = 200000.d0
        northpr = -2705104.270d0+300000d0
        clong = (51+13/60.d0)/radeg
      elseif (izone.eq.97) then
c  irish national grid
        isys = 5
        dm = -0.000035d0
        eastpr = 200000.d0
        northpr = -5929822.893d0+250000d0
        clong = -8.d0/radeg
      elseif (izone.eq.98) then
c  british osgb36
        isys = 6
        dm = 1 - 0.9996012717d0
        eastpr = 400000.d0 
        northpr = -5427063.818d0-100000d0
        clong = -2.d0/radeg
      elseif (izone.eq.99) then
c  swedish
        isys = 4
        dm = 0.d0
        eastpr = 1500000.d0
        clong = 15.8067/radeg
      else
c  UTM
        clong = ((izone - 30)*6 - 3)/radeg
      endif
c
      a = aa(isys)
      f = 1.d0/ff(isys)
c
c  normalized meridian quadrant
c  see k|nig und weise p.50 (96), p.19 (38b), p.5 (2)
c
      n=f/(2.0-f)
      m=n**2*(1.0/4.0+n**2/64.0)
      w= a*(-n - dm+m*(1.0-dm))/(1.0+n)
      sa(1)=a + w
c
c  central easting and longitude
c
      sa(2)=eastpr
      sa(3)=clong
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
      sa(22)=northpr
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
      double precision n, e, n0
c     radeg = 180/3.1415926536d0
c
      n = rn
      e = re
c
      qn = sa(1)
      e0 = sa(2)
      n0 = sa(22)
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
10    np = (n - n0)/qn
      ep = (e - e0)/qn
c
c  ellip. n, e -> sph. n, e
c
      np = np + clcsin(sa, 15, 4, np*2, ep*2, dn, de)
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
      bbg = bbg + clsin(sa, 7, 4, bbg*2)
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
30    bbg   = n + clsin(sa, 3, 4, n*2)
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
40    ep = 1.0d38
41    continue
c
c  sph. normalized n, e -> ell. n, e
c
      np = np + clcsin(sa, 11, 4, np*2, ep*2, dn, de)
      ep = ep + de
      bbg = qn*np + n0
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
c
      subroutine setlamb(rlat0,rlon0)
      implicit double precision(a-h,o-z)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c  set constants needed for geo2lamb
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common /g2lpar/ 
     .re,rfi0,rla0,cosfi0,sinfi0,resin,sinfi2,sinfi12,itype
c
      radeg = 180/3.1415926536d0
      re = 6371000.d0
c
      if (abs(rlat0).eq.90) then
        itype = 1
      elseif(rlat0.eq.0) then
        itype = 3
      elseif (abs(rlat0).lt.90) then
        itype = 2
      else
        stop 'too large latitude in setlamb'
      endif
      rfi0 = rlat0/radeg
      rla0 = rlon0/radeg
      if (itype.ne.2) return
      cosfi0 = cos(rfi0)
      sinfi0 = sin(rfi0)
      resin = re/sinfi0
      sinfi2 = 2*sinfi0
      sinfi12 = 1 + sinfi0**2
      return
      end
c
      subroutine geo2lamb(ldir,xin,yin,xout,yout)
      implicit double precision(a-h,o-z)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  program implements lambert equivalent projections
c  including special cases on equator and pole.
c  Reference: Ricardus and Adler: Map Projections, 1971.     
c
c  (c) Rene Forsberg, KMS, Dec 2001
c
c  ldir: true: geo->lam, false: lam->geo
c 
c  xin, yin:  lat,lon   or N,E  input
c  xout, yout:  N,E     or lat,lon  output
c
c  (c) Rene Forsberg, KMS, Dec 2001
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical ldir
      common /g2lpar/ 
     .re,rfi0,rla0,cosfi0,sinfi0,resin,sinfi2,sinfi12,itype
      data pi / 3.1415926536d0 /
c
      radeg = 180/pi
      goto (100,200,300),itype
c
c  polar equidistant - R&A p167-13
c
100   if (ldir) then
        rfi = xin/radeg
        rla = yin/radeg-rla0
        if (rfi0.gt.0) then
          r = 2*re*sin(pi/4-rfi/2)
        else
          r = 2*re*sin(pi/4+rfi/2)
        endif
        xout = r*cos(rla)
        yout = r*sin(rla)
      else
        r = sqrt(xin**2+yin**2)/(2*re)
        if (rfi0.gt.0) then
          rfi = pi/2 - 2*asin(r)
        else
          rfi = -(pi/2 - 2*asin(r))
        endif
        rla = atan2(yin,xin) + rla0
        xout = rfi*radeg 
        yout = rla*radeg
      endif
      return
c
c  lambert conical equivalent (albers) - R&A p166-11
c 
200   if (ldir) then
        rfi = xin/radeg
        rla = yin/radeg - rla0
        r = sqrt(sinfi12-sinfi2*sin(rfi))
        xout = resin*(cosfi0 - cos(rla*sinfi0)*r)
        yout = resin*sin(rla*sinfi0)*r
      else
        rx = cosfi0-xin/resin
        ry = yin/resin
        rla = atan2(ry,rx)/sinfi0 + rla0
        rfi = asin((sinfi12-(rx**2+ry**2))/sinfi2) 
        xout = rfi*radeg 
        yout = rla*radeg
      endif
      return
c
c  lambert equator cylindrical equivalent - R&A p166-12
c
300   if (ldir) then
        rfi = xin/radeg
        rla = yin/radeg - rla0
        xout = re*sin(rfi)
        yout = re*rla
      else
        rfi = asin(xin/re)
        rla = yin/re + rla0
        xout = rfi*radeg 
        yout = rla*radeg
      endif
      return
      end
c
c ==============================================================
      subroutine geoip(gfile,sfile,ofile,mode)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                            g e o i p
c
c  program for interpolating values from a grid using bilinear or
c  spline interpolation. the grid or the prediction points may be
c  in either geographical or utm coordinates. the spline prediction
c  is performed in a window of size 'nsp' x 'nsp' points around the
c  wanted points, with typical value of nsp being 8 for a good
c  interpolation.
c
c  the grid file must be in standard format, i.e. scanned in e-w bands
c  from n to s, initiated by a label (lat1,lat2,lon1,lon2,dlat,dlon).
c  for utm grid northing and easting replaces lat and lon, and additio-
c  nally ellipsoid number (1: wgs84, 2:ed50, 3:nad27) and utm zone
c  must be given in label. if utm zone 99 is specified, this signals
c  the swedish national projection rt39 (only an approximative
c  transformation is currently implemented, only for geoid use etc.)
c  the program may interpolate from one utm zone to another.
c
c  the grid file may be in direct access binary format, as produced 
c  by program 'gbin'. the use of direct access format speeds up access
c  time. the program recognizes binary files by a special code (777) 
c  written in the first record.
c
c  the program attempts to interpolate all points with a distance
c  'rmin' or more from the margins (for fft applications, e.g., points
c  near the margin are often useless). the program may interpolate in
c  very big grid files, but assumes the prediction points to be reaso-
c  nably close, reading in the smallest necessary subgrid to perform
c  the wanted interpolations.
c
c  special options:
c  - the program may interpolate in  t w o  grids in the same file. this
c  option is especially designed for deflections of the vertical. the
c  two grids must have identical labels. two grid interpolation is
c  signalled by negative mode or mode > 100, see below. two-grid inter-
c  polation can only be done for grid files in txt format.
c  - two grids representing data in different elevations may be used
c  to interpolate values at some elevation between the two grids. in this
c  case mode = mode+100 must be specified.
c  - the program may also be used to convert a list of terrain corrections
c  into rtm-effects through a bouguer reduction to the interpolated
c  reference level. this approximation is only valid for long-wavelength
c  reference grids. stations at sea (negative heights) will have height
c  set to zero.
c  - in bilinear interpolation mode the grid file may contain 9999-values
c  (signals unknown), interpolated values are assigned to 9999 if any
c  unknown values are encountered in the closest 4 points.
c  For grid interpolation the wanted grid may be to large (9999's will be
c  written at unknown values)
c  - north polar stereographic projection implemented in utm
c  (signalled by negative iell - latitude of true scale = abs iell
c                              - longitude of zero = izone)
c
c  input:
c
c  gridfile,
c  outputfile,
c  mode, nsp, rmin, lsel
c  (lat1,lat2,lon1,lon2 - for lsel true)
c  (pointfile - for mode 1 - 3 and modes > 10)
c  (idno - for modes 11-14,19-22 etc.)
c  (lint,lat1,lat2,lon1,lon2,dlat,dlon - for mode 4 only)
c  (iell,izone - for utm only)
c  (h1,h2 - for modes > 100 only)
c
c  where
c
c  mode ...  1: prediction point list in geographic coordinates (degrees)
c            2: do, with lat and lon given in degrees, minutes, seconds
c            3: prediction point list in utm coordinates
c            4: predictions wanted in grid (geographic or utm)
c            5: individual prediction points in lat, lon (degrees)
c            6: do, with lat, lon in degrees, minutes, seconds
c            7: individual prediction points in utm
c
c           10: like 1, with a data value in file written after predictions
c           11: like 1, predictions  s u b t r a c t e d  from values
c               given in file (data value must follow after the height)
c           12: do, but predictions  a d d e d  to values in file
c           13: like 11, with a second data value for each point not changed
c           14: like 12, do
c           15: 'pointfile' contains a grid, from which the
c               interpolated values from 'gridfile' are  s u b t r a c t e d,
c               i.e. 'outfile' = 'pointfile' - 'gridfile'
c           16: 'pointfile' contains a grid, to which the
c               interpolated values from gridfile are  a d d e d.
c           17: 'pointfile' contains a grid which defines the interpolation 
c               points. the given grid values are only used in two-height
c               interpolation mode (117), see below.
c           18: 'pointfile' contains a grid with unknown values (9999).
c               the unknown values are interpolated from the gridfile,
c               other grid values left untouched.
c           19: list of terrain corrections converted to rtm anomalies
c               through bouguer reduction to reference level 'gfile'
c           20: list of free-air anomaly data converted to rtm-reduced
c               data using a bouguer plate approximation only
c               additional input: density
c           21: conversion of free-air data to Bouguer using grid
c               additional input: density
c           22: migration of ERS-1 data over ice caps. grid is slope
c               grid in degrees
c
c           31, 32, ..: like 11, 12, .. for utm coordinates in 'pointfile'
c
c           - if mode is negative, mode = abs mode and  t w o  grid
c           interpolation is performed (the two grids must be in
c           the same file, e.g. deflections of the vertical from
c           program "geofour")
c           - if mode > 100 then point interpolation is done between
c           two grids in different heights using mode = mode-100.
c           the levels to which the two grids refer must be input.
c           if mode=117 then 'pointfile' is assumed
c           to be a height grid defining the interpolation level.
c
c  nsp .... spline window size. nsp=0 means bilinear interpolation,
c           nsp=1 is equivalent to a 8x8 spline (nsp=8)
c
c  rkm .... minimum required distance to closest edge (km) for
c           interpolation to take place
c
c  lsel  .. select only data wihin a subregion (lat1-lat2, lon1-lon2)
c           this option only works with point lists
c
c  idno  .. data number in line (statno,lat,lon,height,d(1),d(2)..)
c
c  lint ... a logical variable specifying that output should be in
c           integer grid format (only needed for mode 3 and 4)
c
c  files:
c
c  unit10  gridfile   grid file in standard format
c  unit20  pointfile  station list file (no, lat, lon, height, ...)
c  unit30  outfile    output file
c  unit31,32          scratch files for intermediate storage of all
c                     prediction points in binary format
c
c  ------  E X A M P L E -----
c
c  a typical job file, where a file 'geoidfile' containing e.g. a binary
c  geoid grid is used to determine values at some points in 'point-
c  file' through linear interpolation would look like:
c
c  geoidfile
c  outputfile
c  1 0 0
c  pointfile
c
c  the 'pointfile' must contain a listing of form
c     1  56.000  10.000   120.0
c     2  55.894  12.029     0.0
c     .....
c  which gives station number, latitude and longitude in degrees,
c  and height in meter. The outputfile will contain the same data
c  followed by the interpolated value. if utm output is wanted (mode 2),
c  lat/lon must be replaced by northing and easting, and ellipsoid
c  number (1: wgs84, 2: hayford, 3: nad27, 4: bessel) and utm zone
c  must additionally be specified (utm zone = 99 signals swedish
c  national projection approximate transformation)
c
c  -----------------------------
c
c  programmer: rene forsberg, danish geodetic institute.
c  university of calgary, mar 87. modified october 87, honefoss.
c  minor changes, rf thessaloniki nov 88
c  revised for swedish map projection and more points, rf june 89
c  modified aug 8, 1990, rf (output coordinates changed)
c  selection added, greenland aug 12, 1992, rf
c  new modes for ice etc., nov 1992, rf
c  allowing grid interpolation outside (with 9999's) in linear mode, dec95
c  allow reading reals as station id (integer on output), alert may04
c  real as output and polarstereographic grids, Oden Aug 2007 
c  Two-grid error corrected, KL Jun 08
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      dimension sa(22),glab(6),glab2(6),ysp(19),rsp(19),qsp(19),
     .hsp(19),sai(22),dd(20)
      character gfile*72,ofile*72,sfile*72
      logical gutm,ggeo,lsec,ltwo,lskip,lrtm,ldif,ladd,lval,lint,lfa,
     .utu,trutu,lchk,l2h,ligrid,l9999,l2hg,lbin,libin,lstdev,lutm,lgeo,
     .lgrid,lindi,lnindi,ldms,ldos,ltran,ldno,lsel,l2d,l10,l21,l22,
     .l2val,lreal
c
c  dimensions: ihadim  grid value arrays
c              igrowd  max number of points in one prediction grid row
c
      real*4 ha,grow
      dimension ha(1900000),grow(10000)
      ihadim =     1900000
      igrowd =                  10000
      ldos = .true.
      lreal = .false.
c
      call openg(10,gfile,lbin)
      open(30,file=ofile,form='formatted',status='unknown')
c
      radeg = 180/3.1415926536d0 
      degkm = 6371.0/radeg
      lsec = .false.
      twopig = 0.1119
c
c      read(*,*) mode, nsp, rmin, lsel

      nsp = 0
      rmin = 0
      lsel = .false.
c
      if (lsel) then
        write(*,*) 'input selection area lat1,lat2,lon1,lon2: '
        read(*,*) xlat1,xlat2,xlon1,xlon2
        nout = 0
      endif
c
      ltwo = (mode.lt.0.or.mode.gt.100)
      mode = abs(mode)
      l2h = (mode.gt.100)
      if (l2h) mode = mode-100  
      lutm = (mode.eq.3.or.mode.eq.7.or.mode.ge.31)
      lgeo = (.not.lutm)
      if (lutm.and.mode.ge.31) mode = mode-20
      ldms = (mode.eq.2.or.mode.eq.6)
      lval = (mode.ge.10)
      l2val = (lval.and.(ltwo.or.(mode.eq.13.or.mode.eq.14)))
      ldno = (mode.ge.10.and.mode.le.14.or.mode.ge.19)
c
      l10 = (mode.eq.10)
      ldif = (mode.eq.11.or.mode.eq.13.or.mode.eq.15)
      ladd = (mode.eq.12.or.mode.eq.14.or.mode.eq.16)
      lstdev = (mode.eq.13.or.mode.eq.14)
      ligrid = (mode.ge.15.and.mode.le.18)
      l2hg = (mode.eq.17.and.l2h)
      l9999 = (mode.eq.18)
      lrtm = (mode.eq.19)
      lfa =  (mode.eq.20)
      l21 = (mode.eq.21)
      l22 = (mode.eq.22)
      l2d = (l10.or.lstdev.or.l22)
      if (l22) then
        write(*,*) 'Input maximal slope allowed: '
        read(*,*) slopemax
      endif
      if (lrtm.or.lfa.or.l21) then 
        write(*,*) 'Input reduction density: '
        read(*,*) dens
        twopig = 0.1119*dens/2.67
      endif
      if (lval) then
        mode = 1
        if (lutm) mode = 3
        if (ligrid) mode = 4
      endif
      lgrid = (mode.eq.4)
      lindi = (mode.ge.5.and.mode.le.7)
      lnindi = .true.
      if (lindi) lnindi = .false.
c
      if (mode.lt.1.or.mode.gt.7) stop 'illegal mode'
      if ((lgrid.or.lrtm.or.lfa).and.l2h.and.(.not.l2hg)) 
     .stop 'two height grid interpolation not allowed'
      if (ltwo.and.lbin) stop 'two grids not allowed in binary format'
      if (ltwo.and.(.not.l2h).and.lstdev) stop 'two grids not allowed'
      if (nsp.gt.19) stop 'spline window too big'
      if (nsp.eq.1) nsp = 8
      if (nsp.eq.2.or.nsp.eq.3) nsp = 4
c
c  read pointfile name if required 
c
      if (mode.le.3.or.ligrid) then
c        if (ligrid) write(*,104)
c        if (.not.ligrid) write(*,105)
c104     format(' input name of gridfile with interpolation points: ')
c105     format(' input name of file with pointdata: ')
c
c  fitgeoid
c         read(*,101) sfile
        if (.not.ligrid)
     .  open(20,file=sfile,form='formatted',status='old')
        if (ligrid) call openg(20,sfile,libin)
      endif
      if (ldno) then
c        write(*,*) 'input position of data in pointfile line: '
c  fitgeoid
c        read(*,*) idno 
        idno = 1
      endif
c
      if (lnindi) open(31,form='unformatted',status='scratch')
      if (ltwo.and.lnindi)
     .open(32,form='unformatted',status='scratch')
c
c  more input - grid specifications, utm ellipsoid and zone
c
      if (ligrid) goto 131
      if (mode.eq.4) then
        write(*,12)
12      format(' input grid spec: lint (f:real,t:int), ',
     .  'n1,n2,e1,e2,dn,de (deg or meter) ')
        if (ldos) write(*,102)
102     format(' ')

        read(*,*) lint,sfi1,sfi2,sla1,sla2,sdfi,sdla
        lutm = (abs(sfi1).gt.100.or.abs(sfi2).gt.100)
      endif
      if (lutm) then
        write(*,13)
13      format(' input: wanted ellipsoid',
     .  ' (1:wgs84, 2:ed50, 3:nad27) and zone ')
        if (ldos) write(*,102)
        read(*,*) ielli,izonei
        if (ielli.eq.0.or.ielli.gt.4.or.abs(izonei).gt.180)
     .  stop 'wanted utm ellipsoid and zone not valid'
      endif
131   continue
c
c  read grid specification of possible pointfile grid
c
      if (ligrid) then
        lint = .false.
        if (libin) then
          read(20,rec=1)  
     .    icode,sfi1,sfi2,sla1,sla2,sdfi,sdla,lutm,ielli,izonei
          if (icode.ne.777) stop 'pointfile grid lacks 777 code'
        else
          read(20,*) sfi1,sfi2,sla1,sla2,sdfi,sdla
          lutm = (abs(sfi1).gt.100.or.abs(sfi2).gt.100)
          if (lutm) read(20,*) ielli,izonei
        endif
        if (lutm.and.
     .  (ielli.lt.1.or.ielli.gt.4.or.izonei.lt.1.or.izonei.gt.99))
     .  stop 'pointfile grid utm ellipsoid and zone not valid'
      endif
c
c  check main grid file for geog or utm projection
c
      if (lbin) read(10,rec=1) icode,glab,gutm,iell,izone
      if (.not.lbin) then
        read(10,*) glab
        gutm = (abs(glab(1)).ge.100.or.abs(glab(2)).ge.100)
        if (gutm) then
          read(10,*) iell,izone
          if (iell.lt.1.or.iell.gt.4.or.izone.lt.1.or.izone.gt.99)
     .    stop 'utm specification on input grid invalid'
        endif
      endif
      dfi = glab(5)
      dla = glab(6)
      iinn = (glab(2)-glab(1))/glab(5) + 1.5
      iine = (glab(4)-glab(3))/glab(6) + 1.5
      ggeo = .true.
      if (gutm) ggeo = .false.
      ltran = (lutm.or.(lgeo.and.gutm))
c
c  two height grid interpolation
c
      if (l2h) then 
        write(*,132)
  132   format(' input levels of two given grids in m: ')
        read(*,*) rlev1, rlev2
        dlev = rlev2-rlev1
        if (dlev.lt.0.001) stop 'height levels must be different'
      endif
c
c  end of input - header output  
c  ----------------------------
c
c      write(*,140) gfile,ofile,mode,nsp,rmin
c140   format(/' ---  G E O I P  ---',
c     ./' grid file name: ',a36,
c     ./' output file name: ',a36, 
c     ./,' mode = ',i1,', nsp = ',i2,
c     .', minimum edge dist ',f5.1,' km')
c      if (mode.lt.3.or.ligrid) write(*,141) sfile
c141   format(' point file name: ',a36)
c      if (lbin) write(*,142)
c142   format(' - grid file assumed to be in binary format -')
c      if (nsp.le.0) write(*,150)
c150   format(' - bilinear interpolation -')
c      if (nsp.gt.0) write(*,160)
c160   format(' - windowed spline interpolation -')
c      if (lfa) write(*,161)
c161   format(' - rtm bouguer ',
c     .'plate reduction of free-air data to reference level -')
c      if (ligrid) write(*,162)
c162   format(' - pointfile assumed to contain grid data -')
c      if (lrtm) write(*,170)
c170   format(' - conversion terrain corrections to rtm effects -')
c      if (ladd) write(*,171)
c171   format(' - addition of interpolated values to pointfile -')
c      if (ldif) write(*,172)
c172   format(' - subtraction of interpolated values from pointfile -')
c      if (l2h) write(*,173) rlev1,rlev2
c173   format(' - interpolation between height levels ',f8.0,' m and',
c     .f8.0,' m')
c      if (l9999) write(*,174)
c174   format(' - unknown (9999) values interpolated from gridfile -')
c      if (l21) write(*,175)
c175   format(' - conversion of free-air data to Bouguer anomalies -')
c      if (l22) write(*,176)
c176   format(' - ERS-1 migration of altimetry over ice caps -')
c
c  check if utm to utm transformation needed
c  in that case transform from izonei to geographic to izone
c
      utu = lutm.and.gutm
      trutu = (utu.and.(iell.ne.ielli.or.izone.ne.izonei))
      if (trutu) write(*,21)
21    format(' - utm to utm transformation required -')
c
      if (ggeo.and.lutm.or.utu) call utmcon(ielli,izonei,sai)
      if (gutm.and.(.not.lutm).or.trutu) call utmcon(iell,izone,sa)
c
c  set up list of prediction points on scratch file
c  entry point for individual points at 22
c  ------------------------------------------------
c
      istat = 0
22    rfimin = 9.9e9
      rlamin = rfimin
      rfimax = -rfimin
      rlamax = -rfimin
      res = 0.0
      if (lgrid) goto 30
c
c  output point list or individual points wanted
c
      n = 1
      if (lindi) then
        if (lutm) write(*,23)
23      format(' input: N, E (meter - 0,0=exit) ')
        if (mode.eq.5) write(*,24)
24      format(' input: lat, lon (degrees - 0,0=exit) ')
        if (mode.eq.6) write(*,241)
241     format(' input: lat, lon (dms - 0,0..0=exit) ')
        if (ldos) write(*,102)
        if (mode.ne.6) then
          read(*,*) rfi,rla 
        else
          read(*,*) ifid,ifim,rfis,ilad,ilam,rlas
          i = 1
          if (ifid.lt.0) i = -1
          rfi = i*(abs(ifid) + ifim/60.0 + rfis/3600.0)
          i = 1
          if (ilad.lt.0) i = -1
          rla = i*(abs(ilad) + ilam/60.0 + rlas/3600.0)
        endif
        if (rfi.eq.0.and.rla.eq.0) goto 92
        istat = istat+1
        goto 26
      endif
c
c  point list reading loop
c
25    if (lval) then
        if (lstdev.or.(l2val.and.(.not.ltwo))) then
          read(20,*,end=29) rstat,rfi,rla,rh,(dd(j),j=1,idno),stdev
          if (l2val) res2 = stdev
        else
          read(20,*,end=29) rstat,rfi,rla,rh,(dd(j),j=1,idno)
        endif
        istat = rstat
        if (istat.ne.rstat) lreal = .true.
        if (idno.eq.0) then
          res = rh
        else
          res = dd(idno)
        endif
      else
        if (ldms) then
          read(20,*,end=29) rstat,ifid,ifim,rfis,ilad,ilam,rlas,rh
          i = 1
          if (ifid.lt.0) i = -1
          rfi = i*(abs(ifid) + ifim/60.0 + rfis/3600.0)
          i = 1
          if (ilad.lt.0) i = -1
          rla = i*(abs(ilad) + ilam/60.0 + rlas/3600.0)
        else
          read(20,*,end=29) rstat,rfi,rla,rh
        endif
        istat = rstat
        if (istat.ne.rstat) lreal = .true.
      endif
      if (lsel) then 
        if (rfi.lt.xlat1.or.rfi.gt.xlat2.or.rla.lt.xlon1.or.rla.gt.
     .  xlon2) then
          nout = nout+1
          goto 25
        endif
      endif
      if (.not.lutm.and.abs(rfi).gt.100)
     .stop 'point list latitude value too large - use utm mode if utm'
c
c  set ocean height to zero for rtm cases
c
      if ((lrtm.or.lfa).and.rh.lt.0) rh = 0.0
c
c  transform point list to utm
c
26    if (ltran) then
        rfio = rfi
        rlao = rla 
      endif
      if (.not.lutm.and.gutm) then
        call utg(rfi/radeg,rla/radeg,rn,re,sa,.false.,.true.)
        rfi = rn
        rla = re
      elseif (lutm.and.ggeo) then
        call utg(rfi,rla,rn,re,sai,.true.,.true.)
        rfi = rn*radeg
        rla = re*radeg
      endif
      if (trutu) then 
        call utg(rfi,rla,rfi,rla,sai,.true.,.true.)
        call utg(rfi,rla,rfi,rla,sa,.false.,.true.)
      endif
c
      if (rfi.lt.rfimin) rfimin = rfi
      if (rfi.gt.rfimax) rfimax = rfi
      if (rla.lt.rlamin) rlamin = rla
      if (rla.gt.rlamax) rlamax = rla
c
c  output prediction point list and find max and min of coordinates
c
      if (lindi) goto 40
      write(31) rstat,rfi,rla,rh,res
      if (ltran) write(31) rfio,rlao
      if (lstdev) write(31) stdev
      if (l2val) write(31) res2
      n = n+1
      goto 25
29    n = n-1
      goto 36
c
c  output points if grid format wanted
c
30    inn = (sfi2-sfi1)/sdfi+1.5
      ine = (sla2-sla1)/sdla+1.5
      if (ine.gt.igrowd) stop 'too many rows in grid - change igrowd'
      n = inn*ine
      k = 0
c
      do 35 i = inn,1,-1
      if (ligrid) then
        if (.not.libin) read(20,*) (grow(j),j=1,ine)
        if (libin) call inrow(20,grow,i,inn,ine)
      endif
      do 35 j = 1,ine
        lchk = (i.eq.1.or.i.eq.inn).and.(j.eq.1.or.j.eq.ine)
        res = grow(j)
        k = k+1
        rfi = sfi1 + (i-1)*sdfi
        rla = sla1 + (j-1)*sdla
        if (.not.lutm.and.gutm) then
          call utg(rfi/radeg,rla/radeg,rn,re,sa,
     .    .false.,.true.)
          rfi = rn
          rla = re
        elseif (lutm.and.ggeo) then
          call utg(rfi,rla,rn,re,sai,.true.,lchk)
          rfi = rn*radeg
          rla = re*radeg
        endif
        if (trutu) then
          call utg(rfi,rla,rn,re,sai,.true.,lchk)
          call utg(rn,re,rfi,rla,sa,.false.,lchk)
        endif
        if (rfi.lt.rfimin) rfimin = rfi
        if (rfi.gt.rfimax) rfimax = rfi
        if (rla.lt.rlamin) rlamin = rla
        if (rla.gt.rlamax) rlamax = rla
c
c  output prediction grid points on scratch file
c
        write(31) rfi,rla,res
35    continue
c
36    if (ggeo) write(*,37) n,rfimin,rfimax,rlamin,rlamax
37    format(' number of prediction points:',i8,/
     .' within area ',4f10.4)
      if (gutm) write(*,38) n,rfimin,rfimax,rlamin,rlamax,izone
38    format(/' number of prediction points:',i8,/
     .' within area ',4f10.0,' in utm zone ',i2)
      if (lsel) 
     .write(*,*) 'number of points outside selection area',nout
c
c  read in necessary heights for grid interpolation
c  entry point for second loop predictions at 45
c  ------------------------------------------------
c
40    if (nsp.gt.0) then
        rnsp = float(nsp)/2-0.5
        rfimin = rfimin - (rnsp+1)*dfi
        rfimax = rfimax + (rnsp+1)*dfi
        rlamin = rlamin - (rnsp+1)*dla
        rlamax = rlamax + (rnsp+1)*dla
      endif
c
45    lskip = (.not.ltwo.or.lsec)
      if (rfimin.gt.glab(2).or.rfimax.lt.glab(1).
     .or.rlamin.gt.glab(4).or.rlamax.lt.glab(3)) then
        if (lindi) then
          write(*,46)      
46        format(' *** point outside grid')
          if (.not.lbin) rewind(10)
          goto 91 
        endif
        stop 'prediction points completely outside grid'
      endif
      call rdelev(10,rfimin,rfimax,rlamin,rlamax,glab,gutm,iell,izone,
     .lbin,lskip,lnindi,rfi1,rla1,dfi,dla,nfi,nla,ha,ihadim)
      rfi2 = rfi1 + (nfi-1)*dfi
      rla2 = rla1 + (nla-1)*dla
c
c  factors for converting gridunits to km
c
      yfakt = 0.001
      xfakt = 0.001
      if (ggeo) yfakt = degkm
      if (ggeo) xfakt = degkm*cos((rfi2+rfi1)/2/radeg)
c
c  station prediction loop initializations
c  ---------------------------------------
c
      if (lindi.and.istat.gt.1) goto 47
      np = 0
      nrp99 = 0
      rrmin = 9.9e9
      rpmin = 9.9e9
      rpmax = -9.9e9
      rpsum = 0
      rpsum2 = 0
      nd99 = 0
      dmin = 9.9e9
      dmax = -9.9e9
      dsum = 0
      dsum2 = 0
      ndd99 = 0
      ddmin = 9.9e9
      ddmax = -9.9e9
      ddsum = 0
      ddsum2 = 0
      kk = 1
      rh = 0.0
      res = 0.0
c
      if (lindi) goto 47
      rewind(31)
      if (ltwo) rewind(32)
47    continue
c
c  prediction point loop
c  ---------------------
c
      do 90 k = 1, n
        if (lindi) goto 48
        istat = k
        if (lgrid) then
          read(31) rfi,rla,res
        else
          read(31) rstat,rfi,rla,rh,res
          if (ltran) read(31) rfio,rlao
          if (lstdev) read(31) stdev
          if (l2val) read(31) res2
        endif
        if (ltwo.and.lsec) read(32) res1
c       if (l2val.and.lsec) res1 = res2
        if (l2hg) rh = res
        if (l9999) rp = res
        if (l9999.and.int(res).ne.9999) goto 70
c
c  check if prediction within grid
c
48      rn = (rfi-rfi1)/dfi + 1
        re = (rla-rla1)/dla + 1
c
c  for linear int
c
        if (nsp.eq.0) then
	  if ((rn.lt.0.499.or.rn.gt.nfi+0.501 .or.
     .    re.lt.0.499.or.re.gt.nla+0.501).and.(.not.lgrid)) then 
            write(*,481) k,istat,rfi,rla
481         format(' *** station',i6,' id',i9,2f13.4,
     .      ' skipped, outside grid')
	    goto 90
          endif
        else
c
c  for spline int
c
          ii1 = ifrac(rn-rnsp+0.5)
          jj1 = ifrac(re-rnsp+0.5)
          ii2 = ii1 + nsp-1
          jj2 = jj1 + nsp-1
          if (ii1.lt.1.or.ii2.gt.iinn.or.jj1.lt.1.or.jj2.gt.iine) then
            if (lgrid) then
	      write(*,49) rfi,rla,ii1,ii2,jj1,jj2
49            format(' fi, la = ',2f13.3,' select index = ',4i4)
              stop 'prediction grid too close to edge for interpolation'
            else
              write(*,50) k,istat,ii1,ii2,jj1,jj2
50            format(/' *** station',i5,' id',i9,' too close to edge -',
     .        ' window index',4i4)
	    endif
            goto 90
	  endif
        endif
c
c  check radius to edge
c
        rr = min((rfi-glab(1))*yfakt, (glab(2)-rfi)*yfakt,
     .  (rla-glab(3))*xfakt, (glab(4)-rla)*xfakt)
        if (rmin.gt.0.and.rr.lt.rmin-0.001.and.mode.ne.4) then
          write(*,52) k,istat,rr
52        format(' *** station ',2i5,' too close to edge - distance ',
     .    f6.1,' km')
          goto 90
        endif
        if (rr.lt.rrmin) rrmin = rr
c
c  bilinear interpolation
c
        if (nsp.gt.0) goto 60
        rp = bilin(rn,re,ha,nfi,nla,ihadim)
        goto 70
c
c  window spline interpolation
c
60      do 65 ii = 1,nsp
          do 62 jj = 1,nsp
            ysp(jj) = ha((ii1+ii-2)*nla + jj1+jj-1)
            if (int(ysp(jj)).ne.9999) goto 62
              write(*,601) k, istat
601           format(' *** station ',2i6,' skipped, unknown value in',
     .        ' spline zone ***')
              goto 90
62        continue
          call initsp(ysp,nsp,rsp,qsp)
          hsp(ii) = spline(re-jj1+1,ysp,nsp,rsp)
65      continue
        call initsp(hsp,nsp,rsp,qsp)
        rp = spline(rn-ii1+1,hsp,nsp,rsp)
c
c  output prediction results
c
70      np = np+1
        rg = rp
c
        if (.not.lval) goto 704
        if (int(rp).eq.9999) ndd99 = ndd99 + 1
        if (int(rp).eq.9999) goto 701
        ddsum = ddsum + rp
        ddsum2 = ddsum2 + rp**2
        if (rp.lt.ddmin) ddmin = rp
        if (rp.gt.ddmax) ddmax = rp
701     if (int(res).eq.9999) nd99 = nd99 +1
        if (int(res).eq.9999) goto 702
        dsum = dsum + res
        dsum2 = dsum2 + res**2
        if (res.lt.dmin) dmin = res
        if (res.gt.dmax) dmax = res
c
c  difference, sum etc. of predicted quantity rg and original data res
c
702     if (int(rp).eq.9999) goto 705
        if (lrtm) then
	  if (res.ge.9999) res = 0
	  rp = twopig*(rh-rg) - res
	endif    
        if (lfa) rp = res - twopig*(rh-rg)
        if (l21) then
          if (rg.ge.0) rp = res - twopig*rg
          if (rg.lt.0) rp = res - twopig*1.64/2.67*rg
        elseif (l22) then
          if (rg.gt.slopemax) then
            rp = 9999.99
          else
            alfa = rg/57.29578d0
            rp = res - 0.5*alfa**2*802000
          endif
          stdev = rg
        endif
        if (ldif) rp = res - rg
        if (ladd) rp = res + rg
c correction rf nov 95
	if (abs(res).ge.9999) rp = 9999
        if (l10) stdev = res
c
704     continue
        if (lsec.and.l2h) rp = (rp-res1)/dlev*(rh-rlev1)+res1
c
705     if (int(rp).eq.9999) nrp99 = nrp99+1
        if (int(rp).eq.9999) goto 706
        rpsum = rpsum + rp
        rpsum2 = rpsum2 + rp**2
        if (rp.lt.rpmin) rpmin = rp
        if (rp.gt.rpmax) rpmax = rp
c
c  output of point format predictions
c
706     if (lgrid) goto 80
        if (ltran) then
          rfi = rfio
          rla = rlao
        endif
c
        if (ltwo) goto 75
        if (lindi) write(*,708) rp
708     format(' interpolated value =',f10.3)
        if (lreal) then
          if (l22) then
            write(30,73) rstat,rfi,rla,rh,res,rp,stdev
          elseif (l2d) then
            if (lgeo) write(30,73) rstat,rfi,rla,rh,rp,stdev
            if (lutm) write(30,74) rstat,rfi,rla,rh,rp,stdev
          else
            if (lgeo) write(30,73) rstat,rfi,rla,rh,rp
            if (lutm) write(30,74) rstat,rfi,rla,rh,rp
          endif
        else
          istat = rstat
          if (l22) then
            write(30,71) istat,rfi,rla,rh,res,rp,stdev
          elseif (l2d) then
            if (lgeo) write(30,71) istat,rfi,rla,rh,rp,stdev
            if (lutm) write(30,72) istat,rfi,rla,rh,rp,stdev
          else
            if (lgeo) write(30,71) istat,rfi,rla,rh,rp
            if (lutm) write(30,72) istat,rfi,rla,rh,rp
          endif
        endif
71      format(' ',i10,2f11.5,f9.2,f11.3,2f9.4)
72      format(' ',i10,2f10.0,f9.2,f11.3,2f9.4)
73      format(' ',f10.5,2f11.5,f9.2,f11.3,2f9.4)
74      format(' ',f10.5,2f10.0,f9.2,f11.3,2f9.4)
        goto 90
c
c  two grid option 
c
75      istat = rstat
        if (lsec) goto 77
        if (ltwo) write(32) rp
        goto 90
c
77      if (lindi) write(*,708) rp
        if (l2h) goto 78
        if (lgeo) write(30,71) istat,rfi,rla,rh,res,rp
        if (lutm) write(30,72) istat,rfi,rla,rh,res,rp
        goto 90
78      if (l2d) then
          if (lgeo) write(30,71) istat,rfi,rla,rh,rp,stdev
          if (lutm) write(30,72) istat,rfi,rla,rh,rp,stdev
        else
          if (lgeo) write(30,71) istat,rfi,rla,rh,rp
          if (lutm) write(30,72) istat,rfi,rla,rh,rp
        endif
        goto 90
c
c  output of grid format predictions, store one row in grow
c  store first component of two-grid interpolation if needed
c
80      if (.not.l2hg.or.lsec) goto 801
          write(32) rp
          goto 90
801     grow(kk) = rp
        kk = kk+1
        if (k.eq.1.and.mode.eq.4.and.(.not.lutm))
     .  write(30,81) sfi1,sfi1+(inn-1)*sdfi,
     .  sla1,sla1+(ine-1)*sdla,sdfi,sdla
81      format(' ',4f12.6,2f12.7)
        if (k.eq.1.and.mode.eq.4.and.lutm)
     .  write(30,82) sfi1,sfi1+(inn-1)*sdfi,
     .  sla1,sla1+(ine-1)*sdla,sdfi,sdla,ielli,izonei
82      format(' ',6f11.1,/,2i5)
        if (kk.le.ine) goto 90
        kk = 1
        if (.not.lint) write(30,83) (grow(j),j=1,ine)
83      format(/,40(' ',7f10.3,/))
        if (lint) write(30,84) (nint(grow(j)),j=1,ine)
84      format(/,40(' ',12i6,/))
c
90    continue
91    if (lindi) then
        if (.not.lbin) rewind(10)
        if (.not.lbin) read(10,*) glab
        if (.not.lbin.and.gutm) read(10,*) iell,izone
        goto 22
      endif
c
c  statistics for predicted data
c
92    if (np.eq.0) stop 'no prediction points left'
      sdev = 0.0
      k = np-nrp99
      if (k.gt.1) sdev = sqrt((rpsum2 - rpsum**2/k)/(k-1))
      if (k.gt.0) rpsum = rpsum/k
      if (.not.lval) goto 94
      ddev = 0.0
      k = np-nd99
      if (k.gt.1) ddev = sqrt((dsum2 - dsum**2/k)/(k-1))
      if (k.gt.0) dsum = dsum/k
      dddev = 0.0
      k = np-ndd99
      if (k.gt.1.and.lval) dddev = sqrt((ddsum2-ddsum**2/k)/(k-1))
      if (k.gt.0.and.lval) ddsum = ddsum/k
c
94    if (lindi) n = istat
      write(*,95) np,n-np,rrmin
95    format(/,' points predicted:',i7,',  skipped points:',i7,
     ./,' minimum distance to grid edges for predictions: ',f6.1,' km', 
     ./,' statistics:                     mean   std.dev.   min   ',
     .'   max    unknown')
      if (lval) write(*,951) dsum,ddev,dmin,dmax,nd99
951   format(' original data             :',4f9.3,i8)
      if (lval) write(*,953) ddsum,dddev,ddmin,ddmax,ndd99
953   format(' grid interpolation results:',4f9.3,i8)
      write(*,952) rpsum,sdev,rpmin,rpmax,nrp99
952   format(' predicted values output   :',4f9.3,i8)
c
c  read and check label for possible second loop
c
      if (.not.ltwo.or.lsec) goto 99
      write(*,96)
96    format(/' --- second grid prediction ---')
      lsec = .true.
      read(10,*) glab2
      if (glab(1).ne.glab2(1).or.glab(2).ne.glab2(2).
     .or.glab(3).ne.glab2(3).or.glab(4).ne.glab2(4).
     .or.glab(5).ne.glab2(5).or.glab(6).ne.glab2(6))
     .stop 'label on second grid not identical to first'
      if (ggeo) goto 45
      read(10,*) iell2,izone2
      if (iell2.ne.iell.or.izone2.ne.izone)
     .stop 'second grid ellipsoid/zone spec not identical to first'
      goto 45
c
99    close(30)
      close(31)
      close(20)
      close(10)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  
c                       r d e l e v                               
c                                                                  
c  subroutine for reading a digital elevation file, to produce a   
c  subgrid of elevations, covering a given area. the grid may be
c  in binary or txt format, with a label (lat1,lat2,lon1,lon2,
c  dlat,dlon) describing limits of grid and spacing. if lat and 
c  and lon are specified in meters, an utm projection is assumed,
c  and next label line must contain ellipsoid number and utm zone.
c
c  parameters:
c
c  fi1,fi2,la1,la2    (call) limits of wanted area     
c  hlab               (call) ready read grid label
c  lutm,iell,izone    (call) utm information     
c  lbin               (call) binary format
c  lskip              (call) skip remaining parts of grid or not
c  lprint             (call) print information and statistics
c  fi1,la1,dfi,dla,nfi,nla  (return) selected grid, fi1,la1 is sw-point
c  h                  (return) array values read
c  idim               (call) dimension statement of h
c                                                  
c  (c) rene forsberg, kms, new version dec 89, updated dec 95
c                                                 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rdelev(iunit,fi1,fi2,la1,la2,hlab,lutm,iell,izone,
     .lbin,lskip,lprint,fi0,la0,dfi,dla,nfi,nla,h,idim)
c
      implicit double precision(a-h,o-z)
      logical lutm,lskip,lbin,lprint,ltxt
      real*4 h, hh
      dimension h(idim), hh(3000)
      dimension hlab(6)
      double precision la0,la1,la2
c
c  parameters for index coordinate system (0,0) in sw-corner
c
      ltxt = (.not.lbin)
      ii = (hlab(2)-hlab(1))/hlab(5) + 1.5
      jj = (hlab(4)-hlab(3))/hlab(6) + 1.5
      if (lprint) write(*,1) 
1     format(/' grid file information:')           
      if (.not.lutm.and.lprint) write(*,2) (hlab(i),i=1,6),ii,jj
2     format(' gridlab:  ',6f9.4,i5,i4)
      if (lutm.and.lprint) write(*,3) (hlab(i),i=1,6),ii,jj,iell,izone
3     format(' gridlab: ',6f9.0,i5,i4,/,' utm ellipsoid ',
     .i1,' zone ',i2)
c
      hlab1 = hlab(1)
      hlab3 = hlab(3)
      dfi = hlab(5)
      dla = hlab(6)
      i0 = ifrac((fi1-hlab1)/dfi + 0.001)
      j0 = ifrac((la1-hlab3)/dla + 0.001)
      if (i0.lt.0) i0 = 0
      if (j0.lt.0) j0 = 0
      fi0 = dfi*i0 + hlab1
      la0 = dla*j0 + hlab3
      ii0 = ifrac((fi2-hlab1)/dfi + 0.999)
      jj0 = ifrac((la2-hlab3)/dla + 0.999)
      if (ii0.ge.ii) ii0 = ii-1
      if (jj0.ge.jj) jj0 = jj-1
      nfi = ii0-i0+1
      nla = jj0-j0+1
c
c  check array dimension
c
      if (nfi*nla.gt.idim) then
        write(*,10) nfi,nla,nfi*nla,idim
10      format(' *** grid array dimension too small - wanted:',
     .  i4,' x ',i4,' =',i7,', declared:',i7,' ***')
        stop 'decrease area or increase array dimensions'
      endif
c
      nfi0 = nfi + i0
      jj1 = max(1, j0+1)
      jj2 = min(jj, j0+nla)
      if (jj1.gt.jj.or.jj2.lt.1) goto 50
c
c  read standard elevation file
c
      if (lbin) goto 40
      do 13 i = ii, 1, -1
        if (i.gt.nfi0)  then
          read(iunit,*) (hh(j),j=1,jj)
        else
          if (i.le.i0.and.lskip) goto 50
          read(iunit,*) (hh(j),j=1,jj)
          if (i.le.i0) goto 13
          k0 = (i-i0-1)*nla - j0
          do 12 j = jj1, jj2
            h(k0+j) = hh(j)
12        continue
        endif
13    continue
      goto 50
c
c  read binary file
c
40    do 20 i = nfi0, i0+1, -1
        call inrow(iunit,hh,i,ii,jj)
        k0 = (i-i0-1)*nla - j0
        do 20 j = jj1, jj2
          h(k0+j) = hh(j)
20    continue
c
c  statistics of read subgrid
c
50    if (.not.lprint) return
      nz = 0
      n9999 = 0
      k = 0
      sum = 0.0
      sum2 = 0.0
      hmin = 9.9e9
      hmax = -9.9e9
c
      do 60 j = 1, nla
      do 60 i = 1, nfi
        rh = h((i-1)*nla+j)
        if (int(rh).eq.9999) then
          n9999 = n9999+1
          goto 60
        endif
        k = k+1
        if (rh.eq.0) nz = nz+1
        if (rh.lt.hmin) hmin = rh
        if (rh.gt.hmax) hmax = rh
        sum = sum+rh
        sum2 = sum2 + rh**2
60    continue
c
      rm = 0.0
      rs = 0.0
      if (k.gt.0) rm = sum/k
      if (k.gt.1) rs = sqrt((sum2 - sum**2/k)/(k-1))
      if (.not.lutm) write(*,71) fi0,fi0+(nfi-1)*dfi,la0,
     .la0+(nla-1)*dla
71    format(' selected subgrid: ',6f9.4)
      if (lutm) write(*,74) fi0,fi0+(nfi-1)*dfi,la0,
     .la0+(nla-1)*dla
74    format(' selected subgrid: ',6f9.0)
      write(*, 72) nfi, nla, nfi*nla, nz, n9999
72    format(' points: ',i3,' x ',i3,' =',i8,', zero values:',i8,
     .', unknown (9999):',i8)
      write(*, 73) hmin, hmax, rm, rs
73    format(' min  max  mean  std.dev.:',4f10.2)
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
c
      subroutine openg(iunit,name,ldirac)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                          o p e n g
c
c  this subroutine will open a grid file to either a text or
c  binary direct access format file. the code '777' must be first
c  value in binary file, for details see 'gbin' program
c  'ldirac' is true for a direct access file, false for a text file
c
c  (c) rene forsberg, kms, dec 89  
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical ldirac
      character*72 name
      ldirac = .true.
      open(unit=iunit,file=name,status='old',
     .form='unformatted',access='direct',recl=64,err=10)
      read(iunit,rec=1,err=11) icode
      if (icode.eq.777) return 
11    close(iunit)
10    ldirac = .false.
      open(unit=iunit,file=name,status='old',form='formatted',err=12)
      return
12    write(*,*) '*** could not open grid file: ',name
      return
      end
c
      subroutine inrow(iunit,grow,i,nn,ne)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                          i n r o w
c
c  reads one row of the grid in binary direct access format. the row
c  records are starting from the north. 'iunit' must be opened with
c  a record length equal to 16 reals
c
c  (c) rf dec 89
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*4 grow(*)
      jrec = (ne-1)/16
      irec = (nn-i)*(jrec+1)+2
      do 10 j = 0, jrec
10    read(iunit,rec=(irec+j)) (grow(j*16+k),k=1,16)
      return
      end 
c
      subroutine outrow(iunit,grow,i,nn,ne)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c                            o u t r o w
c
c  write records on direct access format, to be input by 'inrow'
c  
c  (c) rf dec 89
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*4 grow(*)
      jrec = (ne-1)/16
      irec = (nn-i)*(jrec+1)+2
      do 10 j = 0, jrec
10    write(iunit,rec=(irec+j)) (grow(j*16+k),k=1,16)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  c
c                           b i l i n                              c
c                                                                  c
c  interpolates values in an array a using bilinear                c
c  (parabolic hyperboloid) interpolation.                          c
c                                                                  c
c  parameters:                                                     c
c                                                                  c
c  bilin       interpolated value                                  c
c                                                                  c
c  ri, rj      interpolation argument, (1,1) in lower left corner, c
c              (imax, jmax) in upper right.                        c
c                                                                  c
c  a           integer*2 array with arguments                      c
c                                                                  c
c  imax, jmax  number of points in grid                            c
c                                                                  c
c  iadim       declared dimension of 'a'                           c
c                                                                  c
c  outside area covered by 'a' the function returns                c
c  an unknown code (9999).
c  if unknown (9999) values are met in file
c  the subroutine returns 9999.                            
c                                                         
c  programmer:                                           
c  rene forsberg, july 1983, modified june 1989, dec 1995
c                                                        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function bilin(ri,rj,a,imax,jmax,iadim)
      implicit double precision(a-h, o-z)
      real*4 a
      dimension a(iadim)
c
      if (ri.lt.0.499.or.ri.gt.imax+0.501.or.
     .rj.lt.0.499.or.rj.gt.jmax+0.501) then
	bilin = 9999
	return
      endif
c     
      in = ifrac(ri)
      ie = ifrac(rj)
      rn = ri - in
      re = rj - ie
c
      if (in.lt.1) then
        in = 1
	rn = 0.0
      elseif (in.ge.imax) then
        in = imax-1
        rn = 1.0
      endif
      if (ie.lt.1) then
        ie = 1
        re = 0.0
      elseif (ie.ge.jmax) then
        ie = jmax-1
        re = 1.0
      endif
c
c  let possible corner value give value if close
c
      k = (in-1)*jmax + ie
      a1 = a(k)
      a2 = a(k+jmax)
      a3 = a(k+1)
      a4 = a(k+jmax+1)
      if (int(a1).eq.9999.or.int(a2).eq.9999.or.
     .int(a3).eq.9999.or.int(a4).eq.9999) then
        bilin = 9999
        return
      endif
c
      bilin = (1-rn)*(1-re)*a1 +
     .rn*(1-re)*a2 + (1-rn)*re*a3 +
     .rn*re*a4
      return
      end
c
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
      subroutine set_pst(rlat_true_scale,rlon_down,false_n,false_e)
      implicit double precision(a-h,o-z)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  The function set_pst receives the ellipsoid
c  parameters and Polar Stereograpic projection parameters as inputs, and
c  sets the corresponding state variables.
c
c  rlat_true_scale  : Latitude of true scale, in deg       (input)
c  rlon_down        : Longitude down from pole, in deg     (input)
c  False_E          : Easting (X) at center of projection, in meters  (input)
c  False_N          : Northing (Y) at center of projection, in meters (input)
c
c  POLAR STEREOGRAPHIC originated from :
c                     U.S. Army Topographic Engineering Center
c                     Geospatial Information Division
c                     7701 Telegraph Road
c                     Alexandria, VA  22310-3864
c                     Provided by Scott Spaunhorst, NGA, Aug 2007
c
c  c-code modified to fortran, rf aug 2007
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common /par/ pi,two_pi,
     .origin_long,origin_lat,es,a_mc,tc,two_a,
     .e4,false_eas,false_nor,radeg 
c
      pi = 3.1415926535897d0
      two_pi = 2*pi
      radeg = 180.d0/pi
c
c  wgs84 parameters
c
      a = 6378137.d0     
      f = 1 / 298.257223563d0 
      es = 0.08181919084262188d0
      rmc = 1.0
      tc = 1.0
      e4 = 1.0033565552493d0
      a_mc = 6378137.d0                 
      two_a = 2*a
c               
c  Polar Stereographic projection Parameters 
c
      False_Eas = false_e              
      False_Nor = false_n  
c
      if (rlon_down.gt.180) rlon_down = rlon_down-360
      if (rlat_true_scale.lt.0) then
        Origin_Lat = -rlat_true_scale/radeg
        Origin_Long = -rlon_down/radeg
      else
        Origin_Lat = rlat_True_Scale/radeg
        Origin_Long = rlon_Down/radeg
      endif
c
      es2 = 2*f - f**2
      es = sqrt(es2)
c
      if (abs(abs(Origin_Lat) - pi/2).gt.1.d-10) then
        slat = sin(Origin_Lat)
        essin = es * slat
        pow_es = ((1.0 - EsSin) / (1.0 + EsSin))**(es/2)
        clat = cos(Origin_Lat)
        rmc = clat / sqrt(1.0 - essin * essin)
        a_mc = a * rmc
        tc = tan(pi/4 - Origin_Lat / 2.0) / pow_es
      else
        one_PLUS_es = 1.0 + es
        one_MINUS_es = 1.0 - es
        e4 = sqrt(one_PLUS_es**one_PLUS_es
     .  * one_MINUS_es**one_MINUS_es)
      endif
      return
      end
c
      subroutine geo_to_pst(rlat,rlon,rn,re)
      implicit double precision(a-h,o-z)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  The function geo_to_pst converts geodetic
c  coordinates (latitude and longitude) to Polar Stereographic coordinates
c  (easting and northing), according to the current ellipsoid
c  and Polar Stereographic projection parameters. 
c
c    Latitude   :  Latitude, in radians                      (input)
c    Longitude  :  Longitude, in radians                     (input)
c    Easting    :  Easting (X), in degrees                   (output)
c    Northing   :  Northing (Y), in degrees                  (output)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common /par/ pi,two_pi,
     .origin_long,origin_lat,es,a_mc,tc,two_a,
     .e4,false_eas,false_nor,radeg
c 
      if (abs(abs(rlat)-pi/2) .lt. 1.0d-10) then
        re = 0.d0
        rn = 0.d0
      else
        dlam = rlon - Origin_Long
        if (dlam.gt.pi) dlam = dlam-two_pi
        if (dlam.lt.-pi) dlam = dlam+two_pi
        slat = sin(rlat)
        essin = es * slat
        pow_es = ((1.0 - EsSin) / (1.0 + EsSin))**(es/2)
        t = tan(pi/4 - rlat/2) / pow_es
c
        if (abs(abs(Origin_Lat) - pi/2).gt.1.d-10) then
          rho = a_mc * t / tc
        else
          rho = two_a * t / e4
        endif
        re = rho * sin(dlam) + False_Eas
        rn = -rho * cos(dlam) + False_Nor
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pst_to_geo(rn,re,rlat,rlon)
      implicit double precision (a-h,o-z)
      common /par/ pi,two_pi,
     .origin_long,origin_lat,es,a_mc,tc,two_a,
     .e4,false_eas,false_nor,radeg 
c
c  lat lon in radians
c 
      tempPHI = 0.0d0
c
      dy = rn - False_Nor
      dx = re - False_Eas
      if (dy.eq.0.and.dx.eq.0) then
        rlat = pi/2
        rlon = Origin_Long
      else
        rho = sqrt(dx**2 + dy**2)
        if (abs(abs(origin_lat)- pi/2).gt.1.d-10) then
          t = rho*tc/a_mc
        else
          t = rho*e4/two_a
        endif
c
        PHI = pi/2 - 2*atan(t)
10      if (abs(PHI - tempPHI).gt.1.0d-10) then
          tempPHI = PHI
          sin_PHI = sin(PHI)
          essin =  es * sin_PHI
          pow_es = ((1.0 - EsSin) / (1.0 + EsSin))**(es/2)
          PHI = pi/2 - 2*atan(t*pow_es)
          goto 10
        endif
      endif
c      
      rlat = PHI
      rlon = (Origin_Long + atan2(dx,-dy))
c
      if (rlon.gt.two_pi) rlon = rlon-two_pi
      if (rlon.lt.-two_pi) rlon = rlon+two_pi
      return
      end

