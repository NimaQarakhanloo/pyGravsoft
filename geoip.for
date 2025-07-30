      program geoip
C $Id: geoip.for 290 2009-09-19 17:38:53Z tjansson $
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
c array size enlarged 2002-03-18 by cct.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      dimension sa(22),glab(6),glab2(6),ysp(19),rsp(19),qsp(19),
     .hsp(19),sai(22),dd(20)
      character gfile*72,ofile*72,sfile*72
      logical gutm,ggeo,lsec,ltwo,lskip,lrtm,ldif,ladd,lval,lint,lfa,
     .utu,trutu,lchk,l2h,ligrid,l9999,l2hg,lbin,libin,lstdev,lutm,lgeo,
     .lgrid,lindi,lnindi,ldms,ldos,ltran,ldno,lsel,l2d,l10,l21,l22,
     .l2val
c
c  dimensions: ihadim  grid value arrays
c              igrowd  max number of points in one prediction grid row
c
      real*4 ha,grow
      dimension ha(400000),grow(9000)
      ihadim =     400000
      igrowd =                  9000
      ldos = .true.
c
      write(*,1)
1     format(
     .' *******************************************************',
     .'**************'/
     .' *   GEOIP - GRAVSOFT GRID INTERPOLATION - vers. DEC95  ',
     .'(c) RF/KMS   *'/
     .' *******************************************************',
     .'**************')
      write(*,2)
2     format(' input name of gridfile to be interpolated: ')
      read(*,101) gfile
      write(*,3)
3     format(' input name of outputfile: ')
      read(*,101) ofile
101   format(a72)
      call openg(10,gfile,lbin)
      open(30,file=ofile,form='formatted',status='unknown')
c
      radeg = 180/3.1415926536d0 
      degkm = 6371.0/radeg
      lsec = .false.
      twopig = 0.1119
c
      write(*,10)
10    format(
     .' input: mode  (1:pointfile deg, 2:ptfile dms,',
     .' 3:ptfile utm, 4:grid,'/
     .'               5:lat/lon deg, 6:lat/lon dms, 7:N/E,'/
     .'               11/13:dif, 12/14:sum, 15:grid dif,',
     .' 16:grid sum ..)'/
     .'        nsp   (interpolation type - 0:linear, 1:spline)'/
     .'        rmin  (minimum dist to grid edge, km) '/
     .'        lsel  (true if subarea to be selected)')
      if (ldos) write(*,102)
102   format(' ')
      read(*,*) mode, nsp, rmin, lsel
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
        if (ligrid) write(*,104)
        if (.not.ligrid) write(*,105)
104     format(' input name of gridfile with interpolation points: ')
105     format(' input name of file with pointdata: ')
        read(*,101) sfile
        if (.not.ligrid)
     .  open(20,file=sfile,form='formatted',status='old')
        if (ligrid) call openg(20,sfile,libin)
      endif
      if (ldno) then
        write(*,*) 'input position of data in pointfile line: '
        read(*,*) idno 
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
        read(*,*) lint,sfi1,sfi2,sla1,sla2,sdfi,sdla
        lutm = (abs(sfi1).gt.100.or.abs(sfi2).gt.100)
      endif
      if (lutm) then
        write(*,13)
13      format(' input: wanted ellipsoid',
     .  ' (1:wgs84, 2:ed50, 3:nad27) and zone ')
        if (ldos) write(*,102)
        read(*,*) ielli,izonei
        if (ielli.lt.1.or.ielli.gt.4.or.izonei.lt.1.or.izonei.gt.99)
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
      write(*,140) gfile,ofile,mode,nsp,rmin
140   format(/' ---  G E O I P  ---',
     ./' grid file name: ',a36,
     ./' output file name: ',a36, 
     ./,' mode = ',i1,', nsp = ',i2,
     .', minimum edge dist ',f5.1,' km')
      if (mode.lt.3.or.ligrid) write(*,141) sfile
141   format(' point file name: ',a36)
      if (lbin) write(*,142)
142   format(' - grid file assumed to be in binary format -')
      if (nsp.le.0) write(*,150)
150   format(' - bilinear interpolation -')
      if (nsp.gt.0) write(*,160)
160   format(' - windowed spline interpolation -')
      if (lfa) write(*,161)
161   format(' - rtm bouguer ',
     .'plate reduction of free-air data to reference level -')
      if (ligrid) write(*,162)
162   format(' - pointfile assumed to contain grid data -')
      if (lrtm) write(*,170)
170   format(' - conversion terrain corrections to rtm effects -')
      if (ladd) write(*,171)
171   format(' - addition of interpolated values to pointfile -')
      if (ldif) write(*,172)
172   format(' - subtraction of interpolated values from pointfile -')
      if (l2h) write(*,173) rlev1,rlev2
173   format(' - interpolation between height levels ',f8.0,' m and',
     .f8.0,' m')
      if (l9999) write(*,174)
174   format(' - unknown (9999) values interpolated from gridfile -')
      if (l21) write(*,175)
175   format(' - conversion of free-air data to Bouguer anomalies -')
      if (l22) write(*,176)
176   format(' - ERS-1 migration of altimetry over ice caps -')
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
          read(20,*,end=29) istat,rfi,rla,rh,(dd(j),j=1,idno),stdev
          if (l2val) res2 = stdev
        else
          read(20,*,end=29) istat,rfi,rla,rh,(dd(j),j=1,idno)
        endif
        if (idno.eq.0) then
          res = rh
        else
          res = dd(idno)
        endif
      else
        if (ldms) then
          read(20,*,end=29) istat,ifid,ifim,rfis,ilad,ilam,rlas,rh
          i = 1
          if (ifid.lt.0) i = -1
          rfi = i*(abs(ifid) + ifim/60.0 + rfis/3600.0)
          i = 1
          if (ilad.lt.0) i = -1
          rla = i*(abs(ilad) + ilam/60.0 + rlas/3600.0)
        else
          read(20,*,end=29) istat,rfi,rla,rh
        endif
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
      write(31) istat,rfi,rla,rh,res
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
37    format(/' number of prediction points:',i8,/
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
          read(31) istat,rfi,rla,rh,res
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
        if (l22) then
          write(30,71) istat,rfi,rla,rh,res,rp,stdev
        elseif (l2d) then
          if (lgeo) write(30,71) istat,rfi,rla,rh,rp,stdev
          if (lutm) write(30,72) istat,rfi,rla,rh,rp,stdev
        else
          if (lgeo) write(30,71) istat,rfi,rla,rh,rp
          if (lutm) write(30,72) istat,rfi,rla,rh,rp
        endif
71      format(' ',i10,2f11.5,f9.2,f11.3,2f9.3)
72      format(' ',i10,2f10.0,f9.2,f11.3,2f9.3)
        goto 90
c
c  two grid option 
c
75      if (lsec) goto 77
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
951   format(' original data (pointfile) :',4f9.3,i8)
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
      write(*,991)
991   format(' ')
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
      open(unit=iunit,file=name,status='old',form='formatted')
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
