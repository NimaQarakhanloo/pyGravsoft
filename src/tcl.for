      program tc
      implicit double precision(a-h,o-z)
c $Id: tc.for 267 2009-01-27 14:36:10Z cct $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  
c                               T C                                
c                                                                  
c  Program to compute terrain effects on gravimetric quantities.   
c  various modes of terrain effects may be computed:               
c                                                                  
c  1) Direct topographic effect of all masses above sea-level,     
c     assuming the density to be constant.                         
c                                                                  
c  2) Topographic/isostatic reduction using an airy-model.         
c                                                                  
c  3) Gravimetric terrain correction (effect of topographic        
c     irregularities with respect to a spherical bouguer plate)
c     - in this case with conventional inversion of sign so that       
c     the computed quantities should be a d d e d to observations. 
c                                                                  
c  4) RTM (residual terrain model) effects, i.e. the effects of the 
c     topographic irregularities with respect to a mean surface.   
c     Pre-computed terrain corrections (e.g., from very dense DTM's)
c     can be used, but only for points on the topography.
c
c  The program computes gravity disturbance, deflections of the    
c  the vertical and/or geoid undulations (height anomalies) in     
c  units of mgal, arcsec and m respectively, or second-order       
c  derivatives in eotvos units (1 E = 0.1 mgal/km)                 
c                                                                  
c  The computation is based on two digital elevation models,       
c  a detailed and a coarse, which is used in the inner and outer   
c  zones respectively. the two grids are assumed to have 'common' 
c  boundaries, which is the case if the coarse grid have been con- 
c  structed from the detailed grid by averaging.                   
c                                                                  
c  The integration of the terrain effects is performed using the   
c  formulas for the gravitational effects of a homogeneous rec-    
c  tangular prism. depending on the geometry and accuracy various  
c  formulas are used: exact formulas, spherical harmonic expan-    
c  sion (mc millan) or centered point mass approximation. parame-  
c  ters to determine formula choice and accuracy may be set by the 
c  user in subroutine 'prism1'.                                    
c                                                                  
c  The computation may be done out to a fixed distance from the    
c  computation point or for all masses in a given area. the detai- 
c  led elevation grid is used out to a specified radius, which     
c  should be at least 2 times the gridspacing in the outer grid.   
c  in the local neighbourhood around the point the terrain infor-  
c  mation may be densified using a bicubic spline interpolation.   
c  if the computation point is known to be at the surface of the   
c  topography, the terrain model is modified smoothly in an inner- 
c  zone to give the correct elevation at the computation point.    
c  (a circular area around the computation point may be skipped for
c  for special applications. in this case all prisms with a center 
c  within the specified skip circle is not evaluated. for rtm a  
c  harmonic correction is not done in this case.)               
c  the curvature of the earth is taken into account to the first 
c  order.                                                       
c                                                                
c  If residual terrain effects is wanted, an additional mean     
c  elevation grid (e.g. 30' x 30' mean heights) must be specified.
c  This mean elevation grid should be smooth, e.g. produced by  
c  the program 'tcgrid'                                        
c                                                               
c  The program may use grid files either in utm or geographic    
c  coordinates. If utm grids are used, all grids must be in the   
c  same system (system is specified in gridlabel, see 'rdelev').
c  computation points may be given in geographical coordinates,  
c  but are always output in utm when utm grids are used.        
c
c  The grids must be in standard format, i.e. stored rowwise from N to S
c  with label definining grid, see 'rdelev'
c                                                                
c  Unknown heights in the height data files may be signalled by
c  the value 9999. The program will not output the results for
c  points where the computations have encountered unknown height
c  values. A warning will be written on unit 6, however.
c  NB: a reference grid must not contain unknown height values.
c
c
c  i n p u t                                              
c  *********                                             
c                                                       
c  statfile,           (station list file)
c  dtmfile1,           (detailed elevation grid)
c  dtmfile2,           (coarse elevation grid)
c  dtmfile3,           (reference elevation grid)
c  outfile,            (output file)
c  itype, ikind, izcode, istyp, rho                           
c  fi1, fi2, la1, la2                                         
c  r1, r2,                                                     
c    (rskip - only for r1 < 0),                                  
c    (rtc - only for ikind = 5)
c    (latmin, latmax, lonmin, lonmax, dlat, dlon, elev - for istyp = 0)
c                                                                 
c  dummy file names must be specified if file is not used.        
c  the following codes determine the quantities to be calculated,  
c  the computation kind and the innerzone options:           
c                                                             
c  itype    1  gravity disturbance (dg) - mgal                 
c           2  deflections (ksi, eta) - arcsec                  
c           3  height anomaly (ha) - meter                    
c           4  dg, ksi, eta                                    
c           5  gravity anomaly (dg - ind.eff.)                  
c           6  anomaly, ksi, eta, ha                             
c           7  tzz - vertical gravity gradient (z positive up)
c           8  txx, tyy, tzz
c           9  all gradients
c                                                           
c  ikind    1  topographic effect                            
c           2  isostatic effect                               
c           3  terrain corrections                             
c           4  residual terrain effetcs                         
c           5  do, using precomputed terrain corrections to dist 'rtc' km
c              (if tc-values are missing, ordinary rtm effects are computed)
c                                                            
c  izcode   0  station on terrain, change station elevation   
c           1          do        , change terrain model        
c           2          do        , change terrain model in land points only
c           3  station free                                     
c           4  station free, no spline densification             
c              (dma special - if topography is above computation
c              level, the terrain effects will be computed at the 
c              topography level in both mode 3 and 4)
c                                                               
c  istyp    0  no statfile, compute in grid
c           1  compute effects in statfile
c           2  add effects to value in statfile
c           3  subtract effects     
c           4  statfile with 80-char KMS gravity recs (80-char output
c              for ikind=3, otherwise normal output format)
c         neg  as 1,2 or 3 but UTM data in statfile
c
c  fi1, fi2, la1, la2   fixed maximum area for which ter-     
c              rain effect is computed, irrespectively of the  
c              specified computational radii. unit: deg (utm: m)
c                                                             
c  r1          minimum computation distance of inner grid (km).
c              the inner grid is actually used in the smallest  
c              'subsquare' of the coarse grid, covering a circle 
c              of radius r1. (special option: negative r1 signals
c              computation where topography within a radius rskip 
c              is  n o t  computed. rskip must be input after 'r2')
c                                                              
c  r2          maximal radius of computation (km). if r2 = 0    
c              no outer grid is used. if a fixed-area computation
c              is wanted 'r2' must be sufficiently large.        
c if istyp = 2 or 3 the data column number (nval) must be input.
c                                                                 
c  The station list consists of a                                 
c  listing of no, latitude, longitude (in degrees), elevation, nval
c  data columns.  
c                                                                
c  Marine convention: computations at sea (elevation < 0) are done    
c  at sea level. The depths at a point are shown as a negative height.
c  Caution: Only use izcode=1 for marine data if depths are proper!
c  Use izcode=2 for terrain reductions where sea depths are not used.
c                                                                  
c  Grid option (istyp=0):                                              
c  Computations are performed in the points of a grid, defined by  
c  latitudes 'latmin' to 'latmax' and longitudes 'lonmin' to       
c  'lonmax', with spacing 'dlat' and 'dlon' specified in  
c  degrees.  Elevation of points: 'elev', unless izcode=1 or 2.   
c  Output file is in grid format as well if only one value computed
c                                                                  
c  Note: If the terrain grids does not cover the 'computation      
c  area' as defined by the maximum zone and the radii, elevation   
c  9999 (signalling missing heights) are assumed by the program.
c  Stations encountering missing heights during computations are
c  not computed to end and not output by program.
c  Since the inner grid is always used out to a boundary in the
c  outer grid which is at least 'r1' km away, the 'missing heights'
c  condition may be met even if the point is more than 'r1' away
c  from the inner grid outer boundaries.
c
c  Note: Especially for fixed area reductions, be very careful to  
c  avoid extension of the grids with 9999's. This may be avoided   
c  by specifying a slightly smaller max area, and verified from    
c  the grid specification output of this program.                  
c                                                                  
c  Original version described in ohio state university, dept. of   
c  geodetic science and surveying, report 355, 1984. note that the 
c  grid specification have been  changed  from this report:        
c  limits of grids must allways refer to the  c e n t e r  of the 
c  cells, as natural for point heights.                           
c                                                                  
c  Programmer: Rene Forsberg, Ohio State University / Danish       
c              Geodetic Institute, july 1983.                      
c                                                                  
c  modified for cdc november 1984. grid specification changed      
c  to specify boundaries of grid cell centers. 
c  splineinterpolation have been changed from 5x5 to 7x7 innerzone.
c  an error for rtm-computations of dg has been corrected.         
c  modified for rc8000-fortran, d. arabelos june 86. changes and   
c  inclusion of second order gradients, rf jan 87.                 
c  modified for cdc fortran, rf february 87.                       
c  modified for double precision, 9999-unknown options implemen-
c  ted, rf november 88.
c  modified for swedish national projection, kms unix april 89
c  the swedish national projection is signalled by utm zone = 99 
c  last updated: march 94, rf (80-char input)
c                dec 94, rf (tcnew-dmaac option)
c                jan 96 (add/subtract, inner zone change, exact prism pot)
c                mar 08 (column number add/subtract input).
c                jan 09 error corrected for addition or subtraction.
c                                                                  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      logical dg, dfv, ha, secdev, nogrd2, stgrid, lutm, lskip, sutm
      logical lxx, lyy, lzz, lxy, lxz, lyz, l9999, l80, l80out, lval
      integer*2 h1, h2, href, nval
      real*8 lat,lon,lamin,lamax,la02,la01,la1,la2,
     *lalow,laup,la0ref,t1
      character*1 tccode,tccode2
      character*80 rec
      character*128 statf,dtmf1,dtmf2,dtmf3,outf
      dimension glab(6),sa(22),row(8),t(6),val(6),
     .tmin(6),tmax(6),tsum(6),tsum2(6),
     .vmin(6),vmax(6),vsum(6),vsum2(6),
     .rmin(6),rmax(6),rsum(6),rsum2(6)
c
      data nv,vsum,rsum,tsum /0, 6*0, 6*0, 6*0/
      data vsum2,rsum2,tsum2 /6*0, 6*0, 6*0/
      data vmin,rmin,tmin /6*9.d9, 6*9.d9, 6*9.d9/
      data vmax,rmax,tmax /6*-9.d9, 6*-9.d9, 6*-9.d9/
c
      common /cpar/ itype,ikind,izcode,dg,dfv,ha,r1,r2,r1sq,r2sq,
     .              rski2,rtc2,rho0,rhoiso,dptiso,fktiso,lxx,lyy,lzz,
     .              lxy,lxz,lyz,secdev,lskip,lutm,l9999
      common /result/ sumdg,sumksi,sumeta,sumha,nprism,npriex,
     .              dtmelv,sumxx,sumyy,sumzz,sumxy,sumxz,sumyz
      common /hgridd/ fi01,la01,dfi1,dla1,nfi1,nla1,fi02,la02,dfi2,
     .              dla2,nfi2,nla2
      common /hreff/ fi0ref,la0ref,dfiref,dlaref,nfiref,nlaref
      common /prismp/ r2exac,r2macm,g,gdfv,gha,gsec
c
c  sizes of arrays used to hold elevations adjusted here
c  ----------------------------------------------------- 
c
      dimension h1(3600000), h2(500000), href(360000)
      ih1dm =      3600000
      ih2dm =                  500000
      ihrdm =                                 360000
c
c  read file names
c
      write(*,103)
103   format(' input file names (statfile/dtm1/dtm2/refdtm/outfile):')
      read(*,'(a)') statf
      read(*,102) dtmf1
      read(*,102) dtmf2
      read(*,102) dtmf3
      read(*,102) outf
102   format(a128)
c format statement changer 2008-10-27.
c
      write(*,10)
10    format(' input: itype   (1 dg, 2 defl, 3 ha, 4 g+defl, 5 fa,',
     *' 6 all',/,
     *'                        7 tzz, 8 txx tyy tzz, 9 all grad)',/,
     *'        ikind   (1 topo, 2 iso, 3 tc, 4 rtm, 5 rtm from tc)',/,
     *'        izcode   (0 chg sh, 1 chg ter, 2 chg land, 3 free,',
     *         ' 4 free+nosp)',/,
     *'        istyp   (0 grd, 1 std, 2 sum, 3 dif, 4 80-ch)'/,
     *'        density (g/cm**3) ',/,
     *'        fi1, fi2, la1, la2   (max limits)',/,
     *'        r1, r2   (km, r2=0: no grid2)')
c
      read(*,*) itype, ikind, izcode, istyp, rho0
      read(*,*) fi1, fi2, la1, la2
      read(*,*) r1, r2
c
      lskip = (r1.lt.0)
      r1 = abs(r1)
      if (lskip) then
        write(*,*) 'input: inner skip radius (km)'
        read(*,*) rskip
      endif
      if (ikind.eq.5) then
        write(*,*) 'input: comp. radius of precomputed tc-values (km)'
        read(*,*) rtc
      endif
c
      stgrid = (istyp.eq.0)
      if (stgrid) then
        if (ikind.eq.5) stop '*** grid not allowed for ikind=5'
        write(*, 8001)
8001    format(' input: lat1,lat2,lon1,lon2,dlat,dlon (deg), '
     .  ,'height (m)')
        read(*,*) sfimin, sfimax, slamin, slamax, dfis, dlas, elevs
        nns = (sfimax-sfimin)/dfis+1.5
        nes = (slamax-slamin)/dlas+1.5
        n = nns*nes
      endif
c
      r1 = r1*1000
      r2 = r2*1000
      rski2 = -1.0
      if (lskip) then 
        rskip = rskip*1000
        rski2 = rskip**2
      endif
c
      if (ikind.eq.5) rtc = rtc*1000
c
      dg = (itype.eq.1.or.(itype.ge.4.and.itype.le.6))
      dfv = (itype.eq.2.or.itype.eq.4.or.itype.eq.6)
      ha = (itype.eq.3.or.itype.eq.5.or.itype.eq.6)
      nogrd2 = (r2.le.0)
c
      sutm = (istyp.lt.0)
      istyp = abs(istyp)
      lval = (istyp.eq.2.or.istyp.eq.3.or.ikind.eq.5).and.
     *(itype.eq.1.or.itype.eq.3.or.itype.eq.5.or.itype.eq.7)
C it is only possible to subtractor add a quantity if one value is
c computed.
      if (ikind.eq.5) istyp = 1
      l80 = (istyp.eq.4)
      l80out = (l80.and.ikind.eq.3)
c
      nval = 1
c correction 2009-01-23. Variable itype replaced by istyp.
      if (lval.and.(istyp.eq.3.or.istyp.eq.5)) then
       read(*,*)nval
       write(*,*)' terrain effect added or subtracted from column ',
     * nval
      end if
      if (itype.eq.2) nval = 2
      if (itype.eq.4.or.itype.eq.8) nval = 3
      if (itype.eq.6) nval = 4
      if (itype.eq.9) nval = 6
c
      secdev = (itype.ge.7)
      lxx = (itype.ge.8)
      lyy = (itype.ge.8)
      lzz = (itype.eq.7.or.itype.eq.8)
      lxy = (itype.eq.9)
      lxz = (itype.eq.9)
      lyz = (itype.eq.9)
      lutm = (abs(fi1).ge.100.or.abs(fi2).ge.100)
c
c  density and isostatic paramters
c
      rhoiso = -0.4
      dptiso = 32000.0
      fktiso = abs(1.0/rhoiso)
c
c  prism computation constants
c  r2exac and r2macm determines the shift between different formulas
c
      r2exac = 2.0**2
      r2macm = 9.0**2
      g = 0.00667
      gdfv = 0.001403
      gha = 6.8027d-9
      gsec = 66.7
      radeg = 180.d0/3.141592654d0
c
c  open files
c
      if (.not.stgrid) then
       open(10,file=statf,form='formatted',status='old')
      end if
      open(20,file=dtmf1,form='formatted',status='old')
      if (r2.gt.0) open(21,file=dtmf2,form='formatted',status='old')
      if (ikind.ge.4) open(22,file=dtmf3,form='formatted',status='old')
      open(30,file=outf,form='formatted',status='unknown')
      open(40,status='scratch',form='unformatted')
c
c  read label on detailed file - check for utm
c
      read(20,*) glab
      lutm = (abs(glab(1)).ge.100.or.abs(glab(2)).ge.100)
c
      if (lutm) then
        read(20,*) iell,izone
        if (abs(fi1).lt.100.and.abs(fi2).lt.100) stop
     .  'area specification must be in northing and easting'
        if (iell.lt.1.or.iell.gt.4.or.izone.lt.1.or.izone.gt.99)
     .  stop 'utm specification wrong in dtm file'
        call utmcon(iell,izone,sa)
      endif
c
c  read stations and find max and min of coordinates
c  --------------------------------------------------
c
      fimin = 9.d9
      fimax = -9.d9
      lamin = 9.d9
      lamax = -9.d9
c
      if (.not.stgrid) then
        n = 1
15      if (.not.l80) goto 16
          read(10,151,end=18) rec
          read(rec,152) sfi,sla,elev,istat,tcorr,tccode,ba
          if (tccode.eq.'X') tcorr = 9999.d9
          val(1) = tcorr
          lat = cdeg(sfi)
          lon = cdeg(sla)
151       format(a80)
152       format(1x,f8.2,1x,f9.2,2x,f7.2,27x,i7,f5.1,a1,5x,f7.2)
          goto 17
c
16      if (lval) then
          read(10,*,end=18) istat,lat,lon,elev,(val(j),j=1,nval)
        else
          read(10,*,end=18) istat,lat,lon,elev
        endif
c
17      if (lutm.and.(.not.sutm)) then
          call utg(lat/radeg,lon/radeg,rn,re,sa,.false.,.true.)
          if (i.eq.1) write(*,171) lat,lon,rn,re
171       format(' utm trans first pt: ',2f10.5,1x,2f10.0)
          lat = rn
          lon = re
        endif
        if (lval) then
          write(40) istat,lat,lon,elev,(val(j),j=1,nval)
        elseif (.not.l80out) then
          write(40) istat,lat,lon,elev
        endif
        n = n+1
        if (lat.gt.fimax) fimax = lat
        if (lat.lt.fimin) fimin = lat
        if (lon.gt.lamax) lamax = lon
        if (lon.lt.lamin) lamin = lon
        goto 15
c
18      n = n-1
        if (n.eq.0) stop 'no stations in statfile'
      else
c
c  station grid
c
        istat = 0
        do 20 i = nns, 1, -1
        do 20 j = 1, nes
          istat = istat+1
          lat = sfimin + (i-1)*dfis
          lon = slamin + (j-1)*dlas
          elev = elevs
          if (lutm) then
            call utg(lat/radeg,lon/radeg,rn,re,sa,.false.,.true.)
            if (i.eq.1) write(*,171) lat,lon,rn,re
            lat = rn
            lon = re
            if (lat.gt.fimax) fimax = lat
            if (lat.lt.fimin) fimin = lat
            if (lon.gt.lamax) lamax = lon
            if (lon.lt.lamin) lamin = lon
          endif
          write(40) istat,lat,lon,elev
20      continue
        n = istat
        if (.not.lutm) then
          fimin = sfimin
          fimax = sfimax
          lamin = slamin
          lamax = slamax
        endif
      endif
      continue
c
c  info on parameters and read stations
c  ------------------------------------
c
      write(*, 201)
201   format(/' ===== T C - terrain effect computation =====')
      if (.not.lutm) write(*, 202) itype,ikind,izcode,istyp,
     .fi1,fi2,la1,la2,r1/1000,r2/1000
202   format(' - inputcodes:',4i2,/,'   maxarea:',4f10.2,
     ./,'   computation radii: ', 2f7.1, ' km')
      if (lutm) write(*, 203) itype,ikind,izcode,istyp,fi1,fi2,
     .la1,la2,r1/1000,r2/1000
203   format(' - inputcodes:',4i2,/,'   maxarea:',4f10.0,
     ./,'   computation radii: ', 2f7.1, ' km')
      if (lskip) write(*,204) rskip/1000
204   format(' - no computation within radius ',f6.2,' km -')
      if (.not.lutm) write(*, 205) n, fimin, fimax, lamin, lamax
205   format(' - stations: ',i5,', in area: ',4f9.4)
      if (lutm) write(*, 206) n, fimin, fimax, lamin, lamax
206   format(' - stations: ',i5,', in area: ',4f9.0)
      if (stgrid) write(*, 207) nns, nes, dfis, dlas
207   format('   station grid: ',2i5,', spacing: ',2f9.4)
      if (l80) write(*,*) '- stations input in 80-char format -'
      if (ikind.eq.5) write(*,208) rtc/1000
208   format(' - use of precomputed terrain corrections to',
     .f7.1,' km -')
c
c  read in sufficient elevations to cover station area with the
c  the given calculation radii
c
      f1 = 1.0
      f2 = 1.0
      if (.not.lutm) f1 = 1.0/6371000*radeg
      if (.not.lutm) f2 = f1/cos(fimin/radeg)
c
      do 30 i = 1, 2
        if (i.eq.1) then
          rewind(20)
          dr1 = max(glab(5)/f1, glab(6)/f2)
          if (nogrd2) then
            dr2 = 0
          else
            read(21,*) glab
            if (lutm) read(21,*) k1,k2
            if (lutm.and.(k1.ne.iell.or.k2.ne.izone))
     .      stop 'coarse grid utm spec wrong'
            rewind(21)
            dr2 = max(glab(5)/f1, glab(6)/f2)
          endif
          if (izcode.lt.4) r = max(r1 + dr2, 3.5*dr1)+0.001*dr1
          if (izcode.eq.4) r = max(r1 + dr2, dr1) + 0.001*dr1
        else
          r = r2
          if (nogrd2) then
            fi02 = fi01
            la02 = la01
            dfi2 = dfi1
            dla2 = dla1
            nfi2 = 0
            nla2 = 0
            goto 30
          endif
        endif
        filow = max(fi1, fimin-r*f1)
        fiup  = min(fi2, fimax+r*f1)
        lalow = max(la1, lamin-r*f2)
        laup  = min(la2, lamax+r*f2)
        if (filow.gt.fi2.or.fiup.lt.fi1.or.lalow.gt.la2.
     .  or.laup.lt.la1) 
     .  stop 'wanted computation region outside grid1 coverage' 
c
c  read elevation grids
c  --------------------
c
        if (i.eq.1) then
          write(*, 25)
25        format(/' detailed elevation grid:')
          call rdelev(20, 1, filow, fiup, lalow, laup,
     .    fi01,la01,dfi1,dla1,nfi1,nla1,h1,ih1dm)
          call hinfo(fi01,la01,dfi1,dla1,nfi1,nla1,h1,ih1dm)
        else
          write(*,26)
26        format(/' outer zone elevation grid:')
          call rdelev(21, 1, filow, fiup, lalow, laup,
     .    fi02,la02,dfi2,dla2,nfi2,nla2,h2,ih2dm)
          call hinfo(fi02,la02,dfi2,dla2,nfi2,nla2,h2,ih2dm)
        endif
30    continue
c
c  check relative position of grids
c
      q1 = dfi2/dfi1
      q2 = dla2/dla1
      s1 = (fi01-fi02)/dfi1
      s2 = (la01-la02)/dla1
      frq1=abs(q1-ifrac(q1+0.5))
      frq2=abs(q2-ifrac(q2+0.5))
      frs1=abs(s1-ifrac(s1+0.5))
      frs2=abs(s2-ifrac(s2+0.5))
      if(frq1.gt.0.01.or.frq2.gt.0.01.or.frs1.gt.0.01.or.frs2.gt.0.01)
     . write(*,32) q1,q2,s1,s2
32    format(' *** warning *** outer and inner grid do not fit exactly'
     */,' - grid ratio codes: ',4f12.8)
c
c  read mean elevation grid
c
      if (ikind.ge.4) then
        read(22,*) glab
        if (lutm) read(22,*) k1,k2
        if (lutm.and.(k1.ne.iell.or.k2.ne.izone))
     .  stop 'reference grid utm spec wrong'
        rewind(22)
        write(*, 33)
33      format(/' reference elevation grid:')
        call rdelev(22, -1, filow-glab(5)/2,fiup+glab(5)/2,
     .  lalow-glab(6)/2,laup+glab(6)/2,
     .  fi0ref,la0ref,dfiref,dlaref,nfiref,nlaref,
     .  href,ihrdm)
        call hinfo(fi0ref,la0ref,dfiref,dlaref,nfiref,nlaref,
     .  href,ihrdm)
      endif
c
      write(*,*)
      write(*,*) '========= tc integration results ==========='
      write(*,*)
      if (itype.eq.1.or.itype.eq.5) write(*,*)
     .'     stat    fi        la       h    dtmelv      dg'
      if (itype.eq.2) write(*,*)
     .'     stat    fi        la       h    dtmelv      ksi     eta'
      if (itype.eq.3) write(*,*)
     .'     stat    fi        la       h    dtmelv      ha'
      if (itype.eq.4) write(*,*)
     .'     stat    fi        la       h    dtmelv      ',
     .'dg/anom     ksi       eta'
      if (itype.eq.6) write(*,*)
     .'     stat    fi        la       h    dtmelv      ',
     .'deltag      ksi       eta'
      if (itype.eq.7) write(*,*)
     .'     stat    fi        la       h    dtmelv      tzz'
      if (itype.eq.8) write(*,*)
     .'     stat    fi        la       h    dtmelv      ',
     .'txx      tyy       tzz'
      if (itype.eq.9) write(*,*)
     .'     stat    fi        la       h    dtmelv      '
      if (itype.eq.9) write(*,*)
     .'                txx       tyy       tzz       txy       txz',
     .'       tyz'
c
      n9999 = 0
      ntc9999 = 0
      nprism = 0
      npriex = 0
c
      hdsum = 0
      hdsum2 = 0
      hdmin = 9999999.9
      hdmax = -hdmin
c
      if (l80out) then
        rewind(10) 
        i = r2/1000+.5
        if (nogrd2) i = r1/1000+.5
        tccode = 'D'
        if (i.eq.50) tccode = 'A'
        if (i.eq.166.or.i.eq.167) tccode = 'B'
        if (i.eq.21.or.i.eq.22) tccode = 'C'
      else
        rewind(40)
      endif
c
c  computation station loop
c  -------------------------
c
      do 70 i = 1, n
        if (l80out) then
          read(10,151) rec
          read(rec,152) sfi,sla,elev,istat,tcorr,tccode2,ba
	  if (ikind.ne.3) tccode = tccode2
          lat = cdeg(sfi)
          lon = cdeg(sla)
          if (lutm) then
            call utg(lat/radeg,lon/radeg,rn,re,sa,.false.,.true.)
            lat = rn
            lon = re
          endif
        else
          if (lval) then
            read(40) istat,lat,lon,elev,(val(j),j=1,nval)
          else
            read(40) istat,lat,lon,elev
          endif
        endif
        riih = elev
c
c  set inner radius in mode 5 - if no tc make standard rtm
c
        if (ikind.eq.5) then
          if (val(1).lt.9999) then
            rtc2 = rtc**2
          else
            rtc2 = -1
          endif
        endif
c
c  compute terrain effect
c
        call tcs(lat,lon,elev,h1,h2,href,ih1dm,ih2dm,ihrdm)
c
        if (istyp.ge.2) elev = riih
c
c  write warning for skipped points
c
        if (l9999) then
          n9999 = n9999 + 1
          t(1) = 9999.99d0
          write(*, 60) istat,lat,lon,riih
60        format(' ',i5,2f9.4,' ',f7.1,' - skipped, missing heights')
          goto 65
        endif
c
c  convert gravity disturbance to anomaly for itype 5 and 6
c  small correction oct 90
c
        if (itype.eq.5.or.itype.eq.6) sumdg = sumdg - 0.308*sumha
c
c  complete gradients for itype 9
c
        if (itype.eq.9) sumzz = -(sumxx + sumyy)
c
        goto (601,602,603,604,605,606,607,608,609),itype
601     t(1) = sumdg
        goto 610
602     t(1) = sumksi
        t(2) = sumeta
        goto 610
603     t(1) = sumha
        goto 610
604     t(1) = sumdg
        t(2) = sumksi
        t(3) = sumeta
        goto 610
605     t(1) = sumdg
        goto 610
606     t(1) = sumdg
        t(2) = sumksi
        t(3) = sumeta
        t(4) = sumha
        goto 610
607     t(1) = sumzz
        goto 610
608     t(1) = sumxx
        t(2) = sumyy
        t(3) = sumzz
        goto 610
609     t(1) = sumxx
        t(2) = sumyy
        t(3) = sumzz
        t(4) = sumxy
        t(5) = sumxz
        t(6) = sumyz
610     continue
c
c  write results on * 
c  -------------------
c
        hdif = riih-dtmelv
        hdsum = hdsum + hdif
        hdsum2 = hdsum2 + hdif**2
        if (hdif.lt.hdmin) hdmin = hdif
        if (hdif.gt.hdmax) hdmax = hdif
c
        if (.not.lutm) then
          if (nval.le.4) write(*, 620) istat,lat,lon,
     .    riih,dtmelv,(t(j),j=1,nval)
          if (nval.gt.4) write(*, 621) istat,lat,lon,
     .    riih,dtmelv,(t(j),j=1,nval)
        else
          if (nval.le.4) write(*, 622) istat,lat,lon,
     .    riih,dtmelv,(t(j),j=1,nval)
          if (nval.gt.4) write(*, 623) istat,lat,lon,
     .    riih,dtmelv,(t(j),j=1,nval)
        endif
620     format(' ',i9,2f9.4,' ',2f7.1,' ',6f10.5)
621     format(' ',i9,2f9.4,' ',2f7.1,/,'           ',6f10.5)
622     format(' ',i9,2f9.0,' ',2f7.1,' ',6f10.5)
623     format(' ',i9,2f9.0,' ',2f7.1,/,'           ',6f10.5)
c
c  sum/difference - accumulate statistics
c
        nv = nv + 1
        t1=t(1)
        do j = 1, nval
          if (j.lt.nval.and.lval) then
           t(j)=0.0d0
          end if
          if (j.eq.nval.and.lval) then
           t(j)=t1
          end if
          tsum(j) = tsum(j) + t(j)
          tsum2(j) = tsum2(j) + t(j)**2
          if (t(j).lt.tmin(j)) tmin(j) = t(j)
          if (t(j).gt.tmax(j)) tmax(j) = t(j)
          if (lval) then
c
            if (ikind.eq.5) then
              if (val(j).ge.9999) then
                ntc9999 = ntc9999+1
                val(j) = 0
              endif
              t(j) = t(j)-val(j)
            else
c change 2008-04-08 by cct.
              if (istyp.eq.2) then
               if (j.eq.nval) then
                t(nval) = val(nval)+t1
               else
                 t(j)= val(j)
               end if
              end if
              if (istyp.eq.3) then
                if (j.eq.nval) then
                 t(nval) = val(nval)-t1
                else
                 t(j)= val(j)
                end if
              end if
            endif
c
            vsum(j) = vsum(j) + val(j)
            vsum2(j) = vsum2(j) + val(j)**2
            if (val(j).lt.vmin(j)) vmin(j) = val(j)
            if (val(j).gt.vmax(j)) vmax(j) = val(j)
            rsum(j) = rsum(j) + t(j)
            rsum2(j) = rsum2(j) + t(j)**2
            if (t(j).lt.rmin(j)) rmin(j) = t(j)
            if (t(j).gt.rmax(j)) rmax(j) = t(j)
          endif
        enddo
        if (stgrid.and.nval.eq.1) goto 65
c
c  write results on file
c  ---------------------
c
        if (l80out) then
          read(rec,152) sfi,sla,elev,jstat,tcorr,tccode2,ba
          if (jstat.ne.istat) then
            write (*,*) '*** program error: jstat istat = ',jstat,istat
            stop 'prog error'
          endif
          if (ba.lt.9999) ba = ba-tcorr+sumdg
          tcorr = sumdg
          write(30, 630) rec(1:62),tcorr,tccode,rec(69:73),ba
630       format(a62,f5.1,a1,a5,f7.2)
          goto 70
        else
c
          if (.not.lutm) then
            if (nval.le.4) then
c change 2008-04-02 by cct.
              write(30, 640) istat,lat,lon,elev,(t(j),j=1,nval)
c             write(30, 640) istat,lat,lon,elev,(val(j),j=1,nval)
640           format(' ',i9,f10.5,f11.5,f9.2,4f10.4)
            else
              write(30, 641) istat,lat,lon,elev,(t(j),j=1,nval)
641           format(' ',i9,f10.5,f11.5,f9.2,/,' ',6f10.4)
            endif
          else
            if (nval.le.4) then
              write(30, 642) istat,lat,lon,elev,(t(j),j=1,nval)
642           format(' ',i9,2f10.0,f9.2,4f10.4)
            else
              write(30, 643) istat,lat,lon,elev,(t(j),j=1,nval)
643           format(' ',i9,2f10.0,f9.2,/,' ',6f10.4)
            endif
          endif
        endif
        goto 70
c
c  output on grid format
c
65      if (i.eq.1) then
          jr = 0
          write(30,66) sfimin, sfimax, slamin, slamax, dfis, dlas
66        format(4f11.6,2f11.7)
          write(30,*)
        endif
        jr = jr+1
        row(jr) = t(1)
        if (jr.ge.8.or.mod(i,nes).eq.0) then
          write(30,68) (row(k), k=1,jr)
          if (mod(i,nes).eq.0) write(30,*)
68        format(8f10.3)
          jr = 0
        endif
70    continue
c
c  end of station loop - output statistics
c  ---------------------------------------
c
      if (n9999.gt.0) write(*,89) n9999
89    format(/' --- number of stations skipped due to missing heights: '
     ., i6)
c
      nprism = nprism/nv
      npriex = npriex/nv
      write(*, 90) nprism, npriex
90    format(/' --- average no of prisms/exact formulas pr station: '
     .,2i6)
      if (izcode.le.2.and.nv.gt.1) then
        write(*,91) hdsum/nv,
     .  sqrt((hdsum2-hdsum**2/nv)/(nv-1)),hdmin,hdmax
91      format(' --- difference given - dtm inferred station heights:',
     .  /'     mean stddev min max: ',4f7.0)
      endif
      if (nv.gt.1) then
        write(*,92) nv
92      format (' --- statistics of computed effects, ',
     .  'no of points:', i8/
     .  '       mean    stddev       min       max ')
        do j = 1, nval
          write(*, 93) tsum(j)/nv,
     .    sqrt((tsum2(j)-tsum(j)**2/nv)/(nv-1)),tmin(j),tmax(j)
93        format(' ',4f10.2)
        enddo
        if (lval) then
          write(*,*) '--- statistics of original values in statfile:'
          do j = 1, nval
            write(*, 93) vsum(j)/nv,
     .      sqrt((vsum2(j)-vsum(j)**2/nv)/(nv-1)),vmin(j),vmax(j)
          enddo
          if (istyp.eq.2) write(*,*) '--- sum output on file:'
          if (istyp.eq.3) write(*,*) '--- difference output on file:'
          if (ikind.eq.5) write(*,*) '--- completed RTM effect:'
          do j = 1, nval
            write(*, 93) rsum(j)/nv,
     .      sqrt((rsum2(j)-rsum(j)**2/nv)/(nv-1)),rmin(j),rmax(j)
          enddo
          if (ikind.eq.5) write(*,94) ntc9999
94        format(' --- number of missing, given tc-values (computed ',
     .    'by ordinary rtm):',i7)
        endif
      endif
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  c
c                  s u b r o u t i n e    t c s                    c
c                                                                  c
c  computes terrain effects in individual stations                 c
c                                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine tcs(fi, la, elev, h1, h2, href,ih1dm,ih2dm,ihrdm)
c
      implicit double precision (a-h,o-z)
      logical dg, dfv, ha, splint, lwr1, lwr2, lwr3, lp,
     .lxx, lyy, lzz, lxy, lxz, lyz, l9999, secdev, lutm, lskip
      double precision la,la02,la01,la0ref
      integer*2 h1,h2,href,h,hh
c
      common /cpar/ itype,ikind,izcode,dg,dfv,ha,r1,r2,r1sq,r2sq,
     .              rski2,rtc2,rho0,rhoiso,dptiso,fktiso,lxx,lyy,lzz,
     .              lxy,lxz,lyz,secdev,lskip,lutm,l9999
      common /result/ sumdg,sumksi,sumeta,sumha,nprism,npriex,
     .              dtmelv,sumxx,sumyy,sumzz,sumxy,sumxz,sumyz
      common /hgridd/ fi01,la01,dfi1,dla1,nfi1,nla1,fi02,la02,dfi2,
     .              dla2,nfi2,nla2
      common /hreff/ fi0ref,la0ref,dfiref,dlaref,nfiref,nlaref
      common /dtinfo/ yfakt,xfakt,cfakt,rfi,rla,refdi,refdj,refi0,refj0,
     .              rizsq,stnelv
      dimension h1(ih1dm), h2(ih2dm), href(ihrdm)
c
c  definition parameters for spline densification
c
      dimension  spldef(8), xdiv(16), ydiv(16)
      dimension  a(7), r(7), q(7), splh(15, 7)
c     data  nspl / 8 /
c     data  spldef / 0.0, 0.2, 0.4, 0.55, 0.7, 0.8, 0.9, 0.95 /
      nspl = 8
      spldef(1) = 0.0
      spldef(2) = 0.2
      spldef(3) = 0.4
      spldef(4) = 0.55
      spldef(5) = 0.7
      spldef(6) = 0.8
      spldef(7) = 0.9
      spldef(8) = 0.95
c
      l9999 = .false.
c
      half = 0.5
      one = 1.0
      radeg = 180/3.141592654
c
c  set result variables to zero
c
      sumdg = 0
      sumksi = 0
      sumeta = 0
      sumha = 0
      sumxx = 0
      sumyy = 0
      sumzz = 0
      sumxy = 0
      sumxz = 0
      sumyz = 0
c
c  factors for converting gridunits to meter
c
      yfakt = dfi1
      xfakt = dla1
      cfakt = 1.0
      if (lutm) goto 101
      yfakt = dfi1 * 111195
      cosfi = cos(fi/radeg)
      xfakt = dla1 * 111195 * cosfi
      cfakt = sin(fi/radeg)/cosfi/6371000
c
c  save station elevation, define radii for inner zone and max zone
c  rizsq is the inner zone radius where topography is modified if needed
c  (earlier versions of tc had a 2.0-factor here)
c
  101 stnelv = elev
      r2sq = r2**2
      if (r2.le.0.0) r2sq = r1**2
      rizsq = max(xfakt, yfakt)**2
c
c  grid coordinate system defined with (0, 0) in sw-corner of
c  outer grid, y-axis positive north, x-axis positive east
c  unit is the inner grid unit
c
      iyd = dfi2/dfi1 + 0.5
      jyd = dla2/dla1 + 0.5
      i0 = (fi01-fi02)/dfi1 + 0.5
      j0 = (la01-la02)/dla1 + 0.5
c
      rfi = (fi - fi02)/dfi1
      rla = (la - la02)/dla1
c
      if (ikind.ge.4) then
        refdi = dfiref/dfi1
        refdj = dlaref/dla1
        refi0 = (fi0ref-dfiref/2-fi02)/dfi1
        refj0 = (la0ref-dlaref/2-la02)/dla1
      endif
c
c  grid boundaries: outer zone, inner zone and evt. spline zone
c
      i1y = max(ifrac((rfi - r2/yfakt)/iyd)*iyd, 0)
      i2y = min(ifrac((rfi + r2/yfakt)/iyd+.9999)*iyd, nfi2*iyd)
      j1y = max(ifrac((rla - r2/xfakt)/jyd)*jyd, 0)
      j2y = min(ifrac((rla + r2/xfakt)/jyd+.9999)*jyd, nla2*jyd)
c
      splint = izcode.le.3
      if(.not.splint) go to 2
        ii1 = ifrac(rfi-1.0)
        ii2 = ii1 + 3
        jj1 = ifrac(rla-1.0)
        jj2 = jj1 + 3
  2   continue
c
      ri = r1/yfakt
      rj = r1/xfakt
      if (splint.and.ri.lt.3.5) ri = 3.5
      if (splint.and.rj.lt.3.5) rj = 3.5
      i1 = max(ifrac((rfi-ri)/iyd)*iyd, (i0-1+iyd)/iyd*iyd)
      i2 = min(ifrac((rfi+ri)/iyd+.999999)*iyd, (i0+nfi1)/iyd*iyd)
      j1 = max(ifrac((rla-rj)/jyd)*jyd, (j0-1+jyd)/jyd*jyd)
      j2 = min(ifrac((rla+rj)/jyd+.999999)*jyd, (j0+nla1)/jyd*jyd)
c
c  check station location with respect to grids
c
      nfiyd=nfi2*iyd+0.001
      nljyd=nla2*jyd+0.001
      lwr1=r2.gt.0.and.(rfi.lt.-0.001.or.rfi.gt.nfiyd.or.rla.lt.
     .-0.001.or.rla.gt.nljyd)
      if (lwr1) goto 999
c
      i001=i0-0.001
      i0nfi1=i0+nfi1+0.001
      j001=j0-0.001
      j0nla1=j0+nla1+0.001
      lwr2=rfi.lt.i001.or.
     .rfi.gt.i0nfi1.or.rla.lt.j001.or.rla.gt.j0nla1
      if(.not.lwr2) go to 41
      if (izcode.le.3.or.r2.gt.0) goto 999
      write(*, 9002) fi, la
 9002 format(' ---** tc warning ** station ',2f12.3,
     .' outside innergrid **')
      splint = .false.
      go to 45
   41 continue
c
      ii1m2=ii1-2
      ii2p2=ii2+2
      jj1m2=jj1-2
      jj2p2=jj2+2
      lp=ii1m2.lt.i1.or.ii2p2.gt.i2.or.jj1m2.lt.j1.or.jj2p2.gt.j2
      lwr3=splint.and.lp
      if(.not.lwr3) go to 45
      write(*, 9003) fi,la,ii1,ii2,jj1,jj2
 9003 format(' ---** tc warning ** station ',2f12.3,
     .' too near grid boundary for spline **',4i4)
      splint = .false.
   45 continue
c
c  check boundaries for outside stations
c
      if (r2.gt.0.and.(i1y.ge.i2y.or.j1y.ge.j2y)) return
      if (i1.ge.i2.or.j1.ge.j2) go to 3
      go to 4
   3  continue
        splint = .false.
        i1 = i2y
        i2 = i2y
        j1 = j2y
        j2 = j2y
   4  continue
c
c  computation of inner grid without spline densification
c
      if (splint) go to 32
        ii = ifrac(rfi)+1-i0
        jj = ifrac(rla)+1-j0
        dtmelv=h1((ii-1)*nla1+jj)
        if (dtmelv.eq.9999) goto 999
        if (ii.lt.1.or.ii.gt.nfi1.or.jj.lt.1.or.jj.gt.nla1) dtmelv=0
c
        if (izcode.eq.0.or.izcode.eq.2.and.stnelv.lt.0) 
     .  stnelv = dtmelv
c dmaac special        
        if (izcode.eq.3.or.izcode.eq.4.and.stnelv.lt.dtmelv)
     .  stnelv = dtmelv
c
        do 9 i = i1+1, i2
        do 9 j = j1+1, j2
          hh = h1((i-i0-1)*nla1+j-j0)
          if (hh.eq.9999) goto 999
          call dtc(i-half, j-half, one, one, hh, href,
     .    ihrdm)
   9    continue
        go to 31
   32 continue
c
c  computation of inner grid with spline densification
c
        npoint = 2*nspl - 1
        do 10 i = 1, nspl
          s = spldef(i)
          ydiv(i) = (rfi-ii1)*s
          xdiv(i) = (rla-jj1)*s
          ydiv(npoint+2-i) = 3.0 - (ii2-rfi)*s
          xdiv(npoint+2-i) = 3.0 - (jj2-rla)*s
   10   continue
c
c  compute elevations at horizontal spline lines
c
        do 11 j = 1, 7
          do 111 i = 1, 7
            hh = h1((i+ii1-i0-3)*nla1 + j+jj1-j0-2)
            if (hh.ge.9999) goto 999
            a(i) = hh
  111     continue
          call initsp(a, 7, r, q)
          do 11 i = 1, npoint
            ry = (ydiv(i+1)+ydiv(i))/2
            splh(i,j) = spline(ry+2.5, a, 7, r)
   11   continue
c
c  scan horizontal lines, start with station to get model height
c
        dtmelv = 999999
        do 12 i = nspl, nspl+npoint-1
          ii = i
          if (i.gt.npoint) ii = i - npoint
          ry = (ydiv(ii+1)+ydiv(ii))/2
          do 121 j = 1, 7
            a(j) = splh(ii,j)
            call initsp(a, 7, r, q)
  121     continue
          do 12 j = nspl, nspl+npoint-1
          jj = j
          if (j.gt.npoint) jj = j - npoint
          rx = (xdiv(jj+1)+xdiv(jj))/2
          h = spline(rx+2.5, a, 7, r)
          if (dtmelv.ne.999999) go to 35
            dtmelv = h
            if (izcode.eq.0.or.izcode.eq.2.and.stnelv.lt.0)
     .      stnelv = dtmelv
c dma special   
            if ((izcode.eq.3.or.izcode.eq.4).and.stnelv.lt.dtmelv)
     .      stnelv = dtmelv
c
   35     continue
          call dtc(ry+ii1, rx+jj1, ydiv(ii+1)-ydiv(ii),
     .    xdiv(jj+1)-xdiv(jj), h, href, ihrdm)
   12   continue
c
c  compute rest of innergrid
c
        do 21 i = i1+1, ii1
          k0 = (i-i0-1)*nla1 - j0
          do 21 j = j1+1, j2
            call dtc(i-half, j-half, one, one, h1(k0+j),href,ihrdm)
   21   continue
        if (l9999) goto 999
        do 24 i = ii1+1, ii2
          k0 = (i-i0-1)*nla1 - j0
          do 22 j = j1+1, jj1
            call dtc(i-half, j-half, one, one, h1(k0+j),href,ihrdm)
   22     continue
          do 23 j = jj2+1, j2
            call dtc(i-half, j-half, one, one, h1(k0+j),href,ihrdm)
   23     continue
   24   continue
        if (l9999) goto 999
        do 25 i = ii2+1, i2
          k0 = (i-i0-1)*nla1 - j0
          do 25 j = j1+1, j2
            call dtc(i-half, j-half, one, one, h1(k0+j),href,ihrdm)
   25   continue
   31 continue
      if (l9999) goto 999
c
c  compute outergrid if r2 positive
c
      if (r2.le.0.0) go to 7
        riyd = float(iyd)
        rjyd = float(jyd)
        riyd2 = riyd/2
        rjyd2 = rjyd/2
        do 26 i = i1y+iyd, i1, iyd
          k0 = (i/iyd-1)*nla2
          do 26 j = j1y+jyd, j2y, jyd
            call dtc(i-riyd2,j-rjyd2,riyd,rjyd,h2(k0+j/jyd),
     .      href,ihrdm)
   26   continue
        if (l9999) goto 999
        do 29 i = i1+iyd,i2,iyd
          k0 = (i/iyd-1)*nla2
          do 27 j = j1y+jyd, j1, jyd
            call dtc(i-riyd2,j-rjyd2,riyd,rjyd,h2(k0+j/jyd),
     .      href,ihrdm)
   27     continue
          do 28 j = j2+jyd, j2y, jyd
            call dtc(i-riyd2,j-rjyd2,riyd,rjyd,h2(k0+j/jyd),
     .      href,ihrdm)
   28     continue
   29   continue
        if (l9999) goto 999
        do 30 i = i2+iyd, i2y, iyd
          k0 = (i/iyd-1)*nla2
          do 30 j = j1y+jyd, j2y, jyd
            call dtc(i-riyd2,j-rjyd2,riyd,rjyd,h2(k0+j/jyd),
     .      href,ihrdm)
   30   continue
    7 continue
c
      elev = stnelv
c
c  change sign for conventional gravimetric terrain corrections
c
      if (ikind.ne.3) go to 5
        sumdg = -sumdg
        sumksi = -sumksi
        sumeta = -sumeta
        sumha = -sumha
        if (.not.secdev) goto 5
        sumxx = -sumxx
        sumyy = -sumyy
        sumzz = -sumzz
        sumxy = -sumxy
        sumxz = -sumxz
        sumyz = -sumyz
  5   continue
c
c  harmonic correction for rtm-effects
c  (this correction takes into account that the program so far has
c  computed effects inside the mass when refh > stath)
c
      if (lskip) return
      if (ikind.ge.4.and.dg) go to 6
      return
  6     refh = bilin((rfi-refi0)/refdi, (rla-refj0)/refdj, href,
     .  nfiref, nlaref, ihrdm)
        comph = stnelv
        if (comph.lt.0) comph = 0.0
c       write(*,*)' refh ',refh,comph
        if (refh.gt.comph) sumdg=sumdg-0.083818*rho0*(refh-comph)
c
      return
999   l9999 = .true.
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  c
c                     s u b r o u t i n e   d t c                  c
c                                                                  c
c  computes effect from single integration element                 c
c                                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine dtc(ri, rj, di, dj, hh, href, ihrdm)
c
      implicit double precision (a-h,o-z)
      logical dg, dfv, ha, lxx, lyy, lzz, lxy, lxz, lyz, l9999,
     *secdev, lutm
      double precision la0ref
      integer*2 hh,href
c
      common /cpar/ itype,ikind,izcode,dg,dfv,ha,r1,r2,r1sq,r2sq,
     .              rski2,rtc2,rho0,rhoiso,dptiso,fktiso,lxx,lyy,lzz,
     .              lxy,lxz,lyz,secdev,lskip,lutm,l9999
      common /result/ sumdg,sumksi,sumeta,sumha,nprism,npriex,
     .              dtmelv,sumxx,sumyy,sumzz,sumxy,sumxz,sumyz
      common /hreff/ fi0ref,la0ref,dfiref,dlaref,nfiref,nlaref
      common /dtinfo/ yfakt,xfakt,cfakt,rfi,rla,refdi,refdj,refi0,refj0,
     .              rizsq,stnelv
      dimension href(ihrdm)
c
c  convert from gridunits to meter + meridianconvergence
c
      ym = (ri-rfi)*yfakt
      f = (1.0 - cfakt*ym)
      if (lutm) f = 1.0
      xm = (rj-rla)*xfakt*f
      dist2 = xm**2 + ym**2
c
c  skip far or near sectors
c
      if (dist2.gt.r2sq) return
      if (dist2.le.rski2) return
      if (hh.ge.9999) l9999 = .true.
      if (l9999) return
c
      rho = rho0
      rhow = rho0-1.03
      dx = dj*xfakt*f
      dy = di*yfakt
c
c  change elevations for near sector to match station elevation
c
      h = hh
      if (dist2.lt.rizsq) then
        if (izcode.eq.1 .or. izcode.ge.3.and.stnelv.lt.dtmelv)
     .  h = h + (stnelv-dtmelv)*(rizsq-dist2)/rizsq
        supelv = 0
      else
        supelv = dist2/2/6371000
      endif
c
c  do not allow computations below 0 (neg station height = marine meas)
c
      comph = stnelv+supelv
      if (stnelv.lt.0) comph = supelv
c
c  set mass center and height extent for one or two prisms needed
c
      goto (10,10,20,21,22),ikind
c
c  topography/isostasy
c
10    if (h.lt.0.0) rho = rhow
      dz = h
      zm = dz/2 - comph
      call prism1(rho,xm,ym,zm,dx,dy,dz,dist2+zm**2)
      if (ikind.eq.1) return
      dz = rho*fktiso*h
      zm = -dz/2 - comph - dptiso
      call prism1(rhoiso,xm,ym,zm,dx,dy,dz,dist2+zm**2)
      return
c
c  terrain corrections or rtm effects
c  separate between land and ocean reference surfaces
c
c  tc - marine tc is deviation from Bouguer plate
c
20    refh = stnelv
      if (refh.lt.0) goto 25
      goto 24
c  rtm   
21    refh = bilin((ri-refi0)/refdi,
     .(rj-refj0)/refdj, href, nfiref, nlaref, ihrdm)
c
c   test
c     write(*,777) ri,rj,di,dj,hh,refh
c777  format(/' ri rj di dj hh href =',6f9.4)      
c
      if (refh.lt.0) goto 25
      goto 24
c
c  rtm with tc
c
22    if (dist2.ge.rtc2) goto 21
      refh = stnelv
      h = bilin((ri-refi0)/refdi,
     .(rj-refj0)/refdj, href, nfiref, nlaref, ihrdm)
      rho = -rho
      rhow = -rhow
      if (refh.lt.0) goto 25
c
24    htop = h
      if (htop.lt.0) htop = 0.0
      dz = htop - refh
      zm = dz/2 + refh - comph
      call prism1(rho,xm,ym,zm,dx,dy,dz,dist2+zm**2)
      if (h.ge.0) return
      dz = h
      zm = dz/2 - comph
      call prism1(rhow,xm,ym,zm,dx,dy,dz,dist2+zm**2)
      return
c
25    htop = h
      if (htop.gt.0) htop = 0.0
      dz = htop - refh
      zm = dz/2 + refh - comph
      call prism1(rhow,xm,ym,zm,dx,dy,dz,dist2+zm**2)
      if (h.le.0) return
      dz = h
      zm = dz/2 - comph
      call prism1(rho,xm,ym,zm,dx,dy,dz,dist2+zm**2)
      return
c
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  c
c                s u b r o u t i n e   p r i s m 1                 c
c                                                                  c
c  computes gravity field in orige from rectangular prism with     c
c  coordinates of center (xm, ym, zm) and sidelengths (dx, dy, dz).c
c  exact or approximative formulas are used depending on the prism c
c  geometry.                                                       c
c                                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine prism1(rho, xm, ym, zm, dx, dy, dz, rsq)
c
      implicit double precision (a-h,o-z)
      logical dg, dfv, ha, lxx, lyy, lzz, lxy, lxz, lyz, l9999,
     .secdev, lutm
      common /cpar/ itype,ikind,izcode,dg,dfv,ha,r1,r2,r1sq,r2sq,
     .              rski2,rtc2,rho0,rhoiso,dptiso,fktiso,lxx,lyy,lzz,
     .              lxy,lxz,lyz,secdev,lskip,lutm,l9999
      common /result/ sumdg,sumksi,sumeta,sumha,nprism,npriex,
     .              dtmelv,sumxx,sumyy,sumzz,sumxy,sumxz,sumyz
      common /prismp/ r2exac,r2macm,g,gdfv,gha,gsec
      dimension x(2), y(2), z(2)
      double precision k(2,2,2), h(2,2,2)
c
c  select appropriate formulas
c
      if (abs(dz).le.0.001) return
      nprism = nprism + 1
      dr2 = dx**2 + dy**2 + dz**2
      f2 = rsq/dr2
      if (f2.lt.r2exac) go to 100
c
c  point mass/mac millan
c
      r = sqrt(rsq)
      r3 = r*rsq
      fakt = rho*dx*dy*dz
      if (secdev) goto 30
      if (dg) sumdg = -g*fakt*zm/r3 + sumdg
      if (.not.dfv) go to 1
        sumksi = -fakt*gdfv*ym/r3 + sumksi
        sumeta = -fakt*gdfv*xm/r3 + sumeta
  1   continue
      if (ha) sumha = fakt*gha/r + sumha
      goto 31
30    r5 = r3*rsq
      gg = fakt*gsec
      xm2 = xm**2
      ym2 = ym**2
      zm2 = zm**2
      if (lxx) sumxx = sumxx + gg*(3*xm2 - rsq)/r5
      if (lyy) sumyy = sumyy + gg*(3*ym2 - rsq)/r5
      if (lzz) sumzz = sumzz + gg*(3*zm2 - rsq)/r5
      if (lxy) sumxy = sumxy + gg*3*xm*ym/r5
      if (lxz) sumxz = sumxz + gg*3*xm*zm/r5
      if (lyz) sumyz = sumyz + gg*3*ym*zm/r5
31    if (f2.gt.r2macm) return
c
c mcmillan expansion
c
      alfa = 3*dx**2 - dr2
      beta = 3*dy**2 - dr2
      gamma = 3*dz**2 - dr2
      abg = alfa*xm**2 + beta*ym**2 + gamma*zm**2
      fakt = fakt/24/rsq/r3
      if (secdev) goto 41
      if (dg) sumdg = fakt*g*zm*(gamma - 5*abg/rsq) + sumdg
      if (.not.dfv) go to 2
        sumksi = fakt*gdfv*ym*(beta - 5*abg/rsq) + sumksi
        sumeta = fakt*gdfv*xm*(alfa - 5*abg/rsq) + sumeta
  2   continue
      if (ha) sumha = fakt*gha*abg + sumha
      return
c
41    gg = fakt*gsec
      r4 = rsq**2
      if (lxx) sumxx = sumxx +
     .gg*((35*xm2-5*rsq)/r4*abg - 20*alfa*xm2/rsq + 2*alfa)
      if (lyy) sumyy = sumyy +
     .gg*((35*ym2-5*rsq)/r4*abg - 20*beta*ym2/rsq + 2*beta)
      if (lzz) sumzz = sumzz +
     .gg*((35*zm2-5*rsq)/r4*abg - 20*gamma*zm2/rsq + 2*gamma)
      if (lxy) sumxy = sumxy +
     .gg*(35*xm*ym*abg/r4 - 10*(alfa+beta)*xm*ym/rsq)
      if (lxz) sumxz = sumxz +
     .gg*(35*xm*zm*abg/r4 - 10*(alfa+gamma)*xm*zm/rsq)
      if (lyz) sumyz = sumyz +
     .gg*(35*ym*zm*abg/r4 - 10*(beta+gamma)*ym*zm/rsq)
      return
c
c  exact prism formulas
c
  100 npriex = npriex + 1
      x(1) = xm-dx/2
      x(2) = xm+dx/2
      y(1) = ym-dy/2
      y(2) = ym+dy/2
      z(1) = zm-dz/2
      z(2) = zm+dz/2
      do 11 i1 = 1, 2
      do 11 i2 = 1, 2
      do 11 i3 = 1, 2
        k(i1, i2, i3) = sqrt(x(i1)**2+y(i2)**2+z(i3)**2)
   11 continue
      if (secdev) goto 120
c
      if (dg) sumdg = rho*g*gz(x,y,z,k) + sumdg
c
      if (dfv) then
          do 12 i1 = 1, 2
          do 12 i2 = 1, 2
          do 12 i3 = 1, 2
            h(i1, i2, i3) = k(i2, i3, i1)
   12     continue
          sumksi = rho*gdfv*gz(z, x, y, h) + sumksi
          do 13 i1 = 1, 2
          do 13 i2 = 1, 2
          do 13 i3 = 1, 2
            h(i1, i2, i3) = k(i3, i1, i2)
   13     continue
          sumeta = rho*gdfv*gz(y, z, x, h) + sumeta
      endif
c
      if (.not.ha) return
      ss = 0
      do 14 i1 = 1, 2
      do 14 i2 = 1, 2
      do 14 i3 = 1, 2
        s = 0
        xx = x(i1)
        yy = y(i2)
        zz = z(i3)
        rr = k(i1,i2,i3)
        if (rr.lt.0.05) goto 14
        cc = zz+rr
        if (cc.gt.0) s = s + xx*yy*log(cc)   
        cc = yy+rr 
        if (cc.gt.0) s = s + xx*zz*log(cc)   
        cc = xx+rr
        if (cc.gt.0) s = s + yy*zz*log(cc)   
        if (xx.ne.0) s = s - xx**2/2*atan(yy*zz/xx/rr)
        if (yy.ne.0) s = s - yy**2/2*atan(xx*zz/yy/rr)
        if (zz.ne.0) s = s - zz**2/2*atan(xx*yy/zz/rr)
        ss = ss + (-1)**(i1+i2+i3)*s
   14 continue
      sumha = rho*gha*ss + sumha
      return
c
c  second order derivatives
c
120   gg = rho*gsec
      if (lzz) sumzz = sumzz + gg*tzz(x, y, z, k)
      if (itype.eq.7) return
      if (lxz) sumxz = sumxz + gg*txz(x, y, z, k)
      do 121 i1 = 1,2
      do 121 i2 = 1,2
      do 121 i3 = 1,2
        h(i1,i2,i3) = k(i3,i1,i2)
  121 continue
      if (lxx) sumxx = sumxx + gg*tzz(y, z, x, h)
      if (lxy) sumxy = sumxy + gg*txz(y, z, x, h)
      do 122 i1 = 1,2
      do 122 i2 = 1,2
      do 122 i3 = 1,2
        h(i1,i2,i3) = k(i2,i3,i1)
  122 continue
      if (lyy) sumyy = sumyy + gg*tzz(z, x, y, h)
      if (lyz) sumyz = sumyz + gg*txz(z, x, y, h)
      return
c
      end
c
      double precision function gz(x, y, z, r)
      implicit double precision(a-h,o-z)
      dimension x(2), y(2), z(2), r(2,2,2)
      s = 0.0
      do 17 i1 = 1, 2
        if (abs(x(i1)).le.0.05) go to 17
        s = s+(-1)**i1*x(i1)*log((y(2)+r(i1,2,2))*(y(1)+r(i1,1,1))
     .  /(y(1)+r(i1,1,2))/(y(2)+r(i1,2,1)))
   17 continue
      do 18 i2 = 1, 2
        if (abs(y(i2)).le.0.05) go to 18
        s = s+(-1)**i2*y(i2)*log((x(2)+r(2,i2,2))*(x(1)+r(1,i2,1))
     .  /(x(2)+r(2,i2,1))/(x(1)+r(1,i2,2)))
   18 continue
c
      do 19 i3 = 1,2
        if (abs(z(i3)).le.0.05) go to 19
        do 20 i2 = 1,2
        do 20 i1 = 1,2
          s = s - (-1)**(i1+i2+i3)*z(i3)*atan(x(i1)*y(i2)/z(i3)
     .    /r(i1,i2,i3))
   20   continue
   19 continue
      gz = s
      return
      end
c
      double precision function tzz(x, y, z, r)
      implicit double precision(a-h,o-z)
      dimension x(2), y(2), z(2), r(2,2,2)
      s = 0.0
      do 10 i3 = 1,2
        if (abs(z(i3)).le.0.05) goto 10
        do 11 i2 = 1,2
        do 11 i1 = 1,2
          s = s - (-1)**(i1+i2+i3)*
     .    atan(x(i1)*y(i2)/z(i3)/r(i1,i2,i3))
   11   continue
   10 continue
      tzz = s
      return
      end
c
      double precision function txz(x, y, z, r)
      implicit double precision(a-h,o-z)
      dimension x(2), y(2), z(2), r(2,2,2)
      s1 = (y(2)+r(1,2,1))*(y(1)+r(1,1,2))*
     .     (y(1)+r(2,1,1))*(y(2)+r(2,2,2))
      s2 = (y(1)+r(1,1,1))*(y(2)+r(1,2,2))*
     .     (y(2)+r(2,2,1))*(y(1)+r(2,1,2))
      if (s2.eq.0) goto 10
      s = s1/s2
      if (s.eq.0) goto 10
      txz = log(s1/s2)
      return
   10 txz = 9.99e9
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 
c                       r d e l e v                               
c                                                                 
c  subroutine for reading a digital elevation file, to produce a  
c  subgrid of elevations, covering a given area.                  
c                                                                 
c  elevation be in the format:
c                                                                 
c  1) standard  elevations stored rowvise, from north to south,   
c               with each row scanned from west to east. each     
c               row comprises one record. the first record of     
c               the file is a label record: lat1, lat2, lon1,     
c               lon2, dlat, dlon (degrees), defining              
c               boundary of grid and gridspacing (free format)    
c               for utm northings and eastings are similarly given
c               in meter, followed on a new line by ellipsoid num-
c               ber and utm zone number.                          
c
c  *** subroutine input ***                                       
c                                                                 
c  iunit        fortran unit input number                         
c                                                                 
c  itype        1: grid expanded with 9999's in required
c              -1: do not read grid larger than given
c
c  fi1, fi2, la1, la2   boundaries of wanted area                 
c                                                                 
c  idim1, idim2 maximum (declared) size of h-array                
c                                                                 
c  *** subroutine output ***                                      
c                                                                 
c  fi0, la0     sw-corner of selected grid                        
c                                                                 
c  dfi, dla     grid spacing of selected grid                     
c                                                                 
c  nfi, nla     number of points in latitude and longitude        
c                                                                 
c  h            (integer*2, dimension(idim1,idim2)) elevations,   
c               stored as ('lat index', 'lon index')              
c                                                                 
c  iell, izone  utm ellipsoid no (1: wgs84, 2: int/ed50, 3: nad27)
c               and utm zone no                                   
c                                                                 
c  the subroutine will select the smallest possible grid with the 
c  wanted parameters covering the wanted area. if this area is    
c  larger than the coverage of the elevation file undefined ele-  
c  vations is assumed to be 9999
c                                                                 
c  programmer: rene forsberg, june 1983                           
c  modified grid boundary specification and text format, nov 84   
c  modified for utm, feb 87, u of calgary, rf                     
c  modified for rc8000 binary format and 9999 codes, nov 88, rf   
c  updated Jan 96, rf
c                                                                 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine rdelev(iunit,itype,fi1,fi2,la1,la2,
     .                  fi0,la0,dfi,dla,nfi,nla,h,idim)
c
      implicit double precision (a-h,o-z)
      logical lutm
      integer*2 h,hh
c dimension changed 2009-01-23.
      dimension h(idim), hh(6000),rh(6000)
      double precision hlab(6)
      double precision la1,la2,la0
      data irowdim / 6000 /
c
      read(iunit,*) hlab
      lutm = (abs(hlab(1)).ge.100.or.abs(hlab(2)).ge.100)
      if (lutm) read(iunit,*) iell, izone
c
      if (lutm) then
        if (iell.lt.1.or.iell.gt.4.or.izone.lt.1.or.
     .  izone.gt.99) stop 'utm system specification wrong in file'
      endif
      ii = (hlab(2)-hlab(1))/hlab(5) + 1.5
      jj = (hlab(4)-hlab(3))/hlab(6) + 1.5
      if (.not.lutm) write(*,9000) (hlab(i),i=1,6),ii,jj
9000  format(/' --- gridlab  ',6f9.4,2i5)
      if (lutm) write(*,9010) (hlab(i),i=1,6),ii,jj,iell,izone
9010  format(/' --- gridlab ',6f9.0,2i5,/,'     utm ellipsoid ',
     .i1,' zone ',i2)
c
      dfi = hlab(5)
      dla = hlab(6)
      rfic = hlab(1) - dfi/2
      rlac = hlab(3) - dla/2
      i0 = ifrac((fi1-rfic)/hlab(5)+.001)
      j0 = ifrac((la1-rlac)/hlab(6)+.001)
      if (itype.lt.0) then
        i0 = max(0, i0)
        j0 = max(0, j0)
      endif
      fi0 = dfi*i0 + rfic
      la0 = dla*j0 + rlac
      nfi = ifrac((fi2-fi0)/dfi + .999)
      nla = ifrac((la2-la0)/dla + .999)
      if (itype.lt.0) then
        nfi = max(0, min(nfi, ii-i0))
        nla = max(0, min(nla, jj-j0))
      endif
c
c  check array dimension
c
      if (nfi*nla.gt.idim) write(*,8999) nfi,nla,nfi*nla,idim
 8999 format(' ******** array dimension too small ************',/,
     .' wanted: ',i4,' x ',i4,' =',i7,', declared:',i7)
      if (nfi*nla.gt.idim) stop 'sorry'
      if (nla.gt.irowdim) stop '*** row too long, increase irowdim'
c
c  set elevations to default value 9999
c
      k0 = 0
      do 11 j = 1, nla
      do 11 i = 1, nfi
        k0 = k0+1
        h(k0) = 9999
   11 continue
c
c  read standard elevation file
c
      nfi0 = nfi + i0
      jj1 = max(1, j0+1)
      jj2 = min(jj, j0+nla)
      if (jj1.gt.jj.or.jj2.lt.1) return
c
      do 13 i = ii, 1, -1
        if (i.gt.nfi0) then
          read(iunit,*) (rh(j),j=1,jj)
        else
          if (i.le.i0) return
          read(iunit,*) (rh(j),j=1,jj)
          k0 = (i-i0-1)*nla - j0
          do 12 j = jj1, jj2
            h(k0+j) = nint(rh(j))
 12       continue
        endif
 13   continue
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  c
c                      i f r a c                                   c
c                                                                  c
c  subroutine giving true integer part (entier) of number
c                                                                  c
c  rf, june 1983, modified nov 88                                  c
c                                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer function ifrac(r)
c
      implicit double precision (a-h,o-z)
      if (r.lt.0) go to 1
        ifrac = r
      return
 1      i = r
        if (i.eq.r) go to 2
        ifrac = i-1
        return
 2      ifrac = i
      return
      end
c
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
c  iadim       declared dimension of 'a'
c                                                                  c
c  outside area covered by 'a', the function returns the value of  c
c  the nearest boundary point.                                     c
c                                                                  c
c  programmer:                                                     c
c  rene forsberg, july 1983                                        c
c                                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function bilin(ri,rj,a,imax,jmax,iadim)
c
      implicit double precision (a-h,o-z)
      integer*2 a(iadim)
c
      in = ifrac(ri)
      ie = ifrac(rj)
      rn = ri - in
      re = rj - ie
c
      if(in.lt.1) go to 1
      go to 2
 1      in = 1
        rn = 0.0
      go to 5
 2    if(in.ge.imax) go to 3
      go to 5
 3      in = imax-1
        rn = 1.0
 5    continue
      if(ie.lt.1) go to 4
      go to 7
 4      ie = 1
        re = 0.0
       go to 6
 7    if(ie.ge.jmax) go to 8
      go to 6
 8      ie = jmax-1
        re = 1.0
 6    continue
c
      k = (in-1)*jmax + ie
      bilin = (1-rn)*(1-re)*a(k) +
     .rn*(1-re)*a(k+jmax) + (1-rn)*re*a(k+1) +
     .rn*re*a(k+jmax+1)
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
      dimension y(n), r(n), q(n)
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
      implicit double precision (a-h,o-z)
      dimension y(n), r(n)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  c
c                          h i n f o                               c
c                                                                  c
c  prints information and statistics for an elevation grid, as e.g.c
c  read from a file with 'rdelev'                                  c
c                                                                  c
c  rene forsberg, june 1983                                        c
c                                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine hinfo(fi0,la0,dfi,dla,nfi,nla,h,ihdim)
c
      implicit double precision (a-h,o-z)
      logical lutm
      integer*2 h,hh,hmin,hmax
      double precision la0,la2
      dimension h(ihdim)
c
      fi2 = fi0 + dfi*nfi
      la2 = la0 + dla*nla
      n = nfi*nla
      nz = 0
      n9999 = 0
      sum = 0.0
      sum2 = 0.0
      hmin = 32767
      hmax = -32767
c
      do 10 j = 1, nla
      do 10 i = 1, nfi
        hh = h((i-1)*nla+j)
        if (hh.ge.9999) goto 11
        if (hh.eq.0) nz = nz+1
        if (hh.lt.hmin) hmin = hh
        if (hh.gt.hmax) hmax = hh
        sum = sum+hh
        sum2 = sum2 + (hh*1.0)**2
        goto 10
11      n9999 = n9999+1
10    continue
c
      ns = n - n9999
      if (ns.eq.0) ns = 1
      rm = sum/ns
      if(ns.le.1) go to 1
        rs = sqrt((sum2 - sum**2/ns)/(ns-1))
      go to 2
    1   rs = 0.0
    2 continue
      lutm = (abs(fi0).ge.100.or.abs(fi2).ge.100)
      if (.not.lutm) write(*,9001) fi0+dfi/2,fi2-dfi/2,la0+dla/2,
     .la2-dla/2,dfi,dla
 9001 format(' --- selected ',6f9.4)
      if (lutm) write(*,9010) fi0+dfi/2,fi2-dfi/2,la0+dla/2,
     .la2-dla/2,dfi,dla
9010  format(' --- selected ',6f9.0)
      write(*, 9002) nfi, nla, n, nz, n9999
 9002 format(' points: ',i3,' x ',i3,' = ',i6,', zero values:',i6,
     .', missing/9999:',i6)
      write(*, 9003) hmin, hmax, rm, rs
 9003 format(' hmin  hmax  mean  std.dev.:',2i8,2f11.2)
      write(*, 9004) h((nfi-1)*nla+1),h(nfi*nla),h(1),h(nla)
 9004 format(' corner values: ',4i6)
c
      return
      end
c
      subroutine utmcon(isys, izone, sa)
      implicit double precision (a-h,o-z)
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
c are taken from k|nig und weise : mathematische grundlagen der
c h|heren geod<sie und kartographie, erster band, berlin 1951.
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
c rc fortran version alp/rf oct 86, last updated apr 89
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision n,m
      radeg = 180/3.1415926536
c
c  set ellipsoid parameters
c  check for swedish special projection
c
      if (izone.eq.99) isys = 4
      goto (10,20,30,31),isys
c
c  wgs84 ellipsoid
c
10    a = 6378137.0d0
      f = 1/298.2572236d0
      goto 40
c
c  hayford ed 50 ellipsoid
c
20    a = 6378388.0d0
      f = 1/297.0d0
      goto 40
c
c  clarke nad 27 ellipsoid
c
30    a = 6378206.4d0
      f = 1/294.9786982d0
      goto 40
c
c  bessel ellipsoid
c
31    a = 6377397.155d0
      f = 1/299.153d0
c
40    eastpr = 500000.0
      dm = 4.0e-4
      if (izone.eq.99) dm = 0.0
      if (izone.eq.99) eastpr = 1500000.0
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
      if (izone.ne.99) sa(3)=((izone - 30)*6 - 3)/radeg
      if (izone.eq.99) sa(3)=15.8067/radeg
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
c
      double precision function cdeg(r)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                                 C D E G
c
c  changes number of form -ddmm.mm into true degrees. used for reading
c  numbers from data base.
c  rf jan 1992
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      i = 1
      if (r.lt.0) i = -1
      rr = abs(r)
      ideg = rr/100
      rmin = rr - ideg*100
      cdeg = i*(ideg + rmin/60)
      return
      end
