      program geoid
c $Id: geoid.for 208 2008-09-08 13:37:44Z cct $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      G E O I D 
c
c  simple linear interpolation and transformation program for
c  geoid gridfiles. This program is an alternative to 'geoip', and
c  differs a.o. in the way that subgrid are are not stored in memory,
c  but accessed everytime from disc. 
c
c  the program may also be used for transformation between UTM, geographic
c  and cartesian coordinate systems, may convert ellipsoidal to orthometric
c  heights and vice versa, and may also be used to fit the geoid onto local
c  control to eliminate residual geoid errors and vertical datum discrepan-
c  cies (included in the ellipsoidal to orthometric height conversion).
c  
c  Input to the program is interactive. 
c
c  The program runs in a number of basic modes:
c
c  0: just coordinate conversion (i.e., geoid height assumed to be zero)
c  1: interpolate geoid heights in same system as in file
c  2: convert ellipsoidal to "normal" heights
c  3: convert "normal" heights to ellipsoidal heights
c  4: interpolate geoid height and convert to other datum (e.g. ED50 geoid)
c  5: deflections of the vertical (simple linear differentiation of geoid)
c  6: deflections of the vertical in local system (e.g. ED50 deflections)
c
c  NOTE: The geoid height given is always assumed to be geocentric  (WGS84).
c  The coordinate transformations in mode 0-3 only relate to the transformation
c  for obtaining latitude and longitude, for interpolating grid values.
c
c  NB: Coordinates are input in the form 56 12 34.0 for geographic
c  coordinates in degrees, minutes and seconds. Longitudes W of
c  Greenwich must be negative (for 0 degree W the sign must be on
c  the minute or second, e.g. 3'53"W = 0 -3 53, 0'53"W = 0 0 -53) 
c  
c  (c) Rene Forsberg, Kort- og Matrikelstyrelsen, dec 1990.
c  Last update feb 1992.
C  Update 2008-09-08.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      character*112 file,gfile,ffile,dfile,ofile
      character*1 char
      logical lutm,linp,lih,lneg,ldos,lok,ltrend,lzero,lfirst,ldfv,lds
      dimension sag(22),sai(22),sao(22),c(66),tcoef(11),ov(11),hh(2)
      common
     *  /gpar/ itask,itask1,rfi1,rla1,dfi,dla,nn,ne,lutm,lutm1,
     *  sag,radeg,id,nrej
      common /tpar/ cosfi0,rfi0,rla0
c
c  default geoid name
c
c     data gfile /'/usr4/geo/rf/GGEOID/geoid92a.bin'/
      data gfile /' '/
      data dfile /' '/
      data ofile /' '/
c
      ldos = .false.
      radeg = 180.0/3.1415926536d0
c
      write(*,2)
2     format(/
     .3x,
     .'***************************************************************'/
     .3x,
     .'*                                                             *'/
     .3x,
     .'*   GEOID - GRAVSOFT geoid interpolation and transformation   *'/
     .3x,
     .'*                                                             *'/
     .3x,
     .'*   vers. MAR92 (c) RF, Kort- og Matrikelstyrelsen, Denmark   *'/
     .3x,
     .'***************************************************************'/
     .' '/
     .' Enter task: 1 = interpolate geoid heights'/
     .'             2 = ellipsoidal to orthometric heights using geoid'/
     .'             3 = orthometric heights to ellipsoidal   -     -'/
     .'             4 = geoid heights in different datum ...')
      if (ldos) write(*,*) '-> '
      read(*,*) itask
c
      ldfv = .false.
      if (itask.eq.5.or.itask.eq.6) then
        write(*,*)
     .  ' --- deflections of the vertical wanted, unit: arcsec ---'
        ldfv = .true.
        if (itask.eq.5) itask = 1
        if (itask.eq.6) itask = 4
      endif
c
      lds = .false.
      if (itask.eq.4) then
        write(*,7)
7       format(' Input wanted datum for geoid prediction:'/,
     .  ' 2 = ED50, 3 = NAD27, 4 = Qornoq: ')
        if (ldos) write(*,*) '-> '
        read(*,*) ids
        if (ids.eq.5) stop 'NWL9D geoid not implemented'
        lds = .true.
        itask = 1  
      endif
c     
      if (itask.gt.0) then
8       write(*,3) gfile
3       format(' Enter binary geoid file name: CR=',a37)
        if (ldos) write(*,*) '-> '
        read(*,1) file
1       format(a72)
        if (file(1:1).ne.' ') gfile = file 
        open(10,file=gfile,status='old',
c     .  iointent='input',
     .  form='unformatted',access='direct',recl=64,err=9)
        goto 10
9         write(*,*) '*** grid file unknown/can not be opened ***'
          goto 8
10      read(10,rec=1,err=11) icode,rfi1,rfi2,rla1,rla2,dfi,dla,
     .  lutm,iellg,izoneg
        goto 12
11        write(*,*) '*** grid file not in binary format, sorry ***'
12      if (icode.ne.777)
     .  stop '*** grid file not in correct format, 777-code missing ***'
        nn = (rfi2-rfi1)/dfi+1.5
        ne = (rla2-rla1)/dla+1.5
        if (nn.lt.2.or.ne.lt.2) stop '*** geoid grid too small ***'
        if (lutm) then
          write(*,13) izoneg,rfi1,rfi2,rla1,rla2,dfi,dla
13        format(' Geoid grid in UTM zone ',i2,', Northing and Easting',
     .    ' limits and spacing:'/' ',6f10.0)
          call utmcon(iellg,izoneg,sag)
        else
          write(*,14) rfi1,rfi2,rla1,rla2,dfi,dla
14        format(' Geoid grid limits and spacing in degrees:'/,
     .    2f12.5,f10.5,2f12.5,f10.5)
        endif  
      endif
c
c  possible trend elimination of the geoid from given heights
c  ----------------------------------------------------------
c
      if (itask.eq.2) then
        write(*,30)
30      format(
     .  ' Do you wish to fit the geoid to local ',
     .  'orthometric heights? (T/F)')
        if (ldos) write(*,*) '-> '
        read(*,*) ltrend
        if (ltrend) then
32        write(*,31)
31        format(' Heights of the local points must be ',
     .    'given in a separate file of form'/
     .    '   1 = statno, lat, lon, h, H (degrees)'/
     .    '   2 = statno, lat, lon, h, H (deg,min,sec)'/
     .    '   3 = statno, lat, lon, H, geoid height (degrees)'/
     .    ' (where h = ellipsoidal height, H = orthometric height)'/
     .    ' File name: ')
          if (ldos) write(*,*) '-> '
          read(*,1) ffile
          open(21,file=ffile,form='formatted',status='old',err=33)
          goto 34
33          write(*,*) '*** file could not be opened'
            goto 32
34        write(*,35)
35        format(' File format: ')
          if (ldos) write(*,*) '-> '
          read(*,*) iform
          if (iform.lt.1.or.iform.gt.3) stop '*** wrong format'
          write(*,36)
36        format(
     .    ' Select type of fit function:'/
     .    '   1 = bias fit (requires at least 1 known point)'/
     .    '   2 = fit linear trend (requires min. 3 points or more)'/
     .    '   3 = fit second order polynomial in x and y',
     .    ' (min. 6 points)'/
     .    '   4 = fit third order polynomial (min. 10 points)'/
     .    '   5 = fit datumshift parameters dX,dY,dZ (min. 3 points)'/
     .    '   6 = fit datumshift dX,dY,dZ,scale (min. 4 points)')
          if (ldos) write(*,*) '-> '
          read(*,*) itrend
          if (itrend.lt.1.or.itrend.gt.6) stop '*** wrong fit type'
c
c  scan through points in the given heights file
c
          np = 0
          do 40 i = 1, 66
40        c(i) = 0.0
c
41        if (iform.eq.1) then
            read(21,*,end=50) id,rfi,rla,rell,rh
            rn = rell-rh
          elseif (iform.eq.2) then
            read(21,*,end=50) id,ifi,mfi,sfi,ila,mla,sla,rell,rh
            ii = 1
            if (ifi.lt.0.or.mfi.lt.0.or.sfi.lt.0) ii = -1
            rfi = ii*(abs(ifi) + abs(mfi)/60.d0 + abs(sfi)/3600)
            ii = 1
            if (ila.lt.0.or.mla.lt.0.or.sla.lt.0) ii = -1
            rla = ii*(abs(ila) + abs(mla)/60.d0 + abs(sla)/3600)
            rn = rell-rh
          else
            read(21,*,end=50) id,rfi,rla,rh,rn
          endif
c
          np = np+1
          if (np.eq.1) then
            rfi0 = rfi
            rla0 = rla
            cosfi0 = cos(rfi0/radeg)
            write(*,42)
42          format(' Apparent geoid errors (GPS-gravimetric) at given',
     .      ' points:')
          endif
c
          rgeoid = geoid1(rfi,rla,lok)
          if (.not.lok)
     .    stop '*** missing interpolated geoid height at given point'
          dn = rn-rgeoid
          call setov(rfi,rla,dn,itrend,nt,ov)
          ii = 0
          do 43 i = 1, nt+1
          do 43 j = 1, i
            ii = ii+1
            c(ii) = c(ii) + ov(i)*ov(j)
43        continue
          write(*,44) id,rfi,rla,rh,dn
44        format(' ',i8,2f12.5,2f11.3)
          goto 41
c
c  solve trend equations
c
50        call chol(c,nt,nsing)
          if (nsing.ne.0) write(*,51) nsing
51        format(' *** warning: detrending gave '
     .    ,i2,' normal equation singularities')
          k = nt*(nt+1)/2
          do 52 i = 1, nt
52        tcoef(i) = c(k+i)
          if (itrend.lt.5) write(*,53) (tcoef(i),i=1,nt)
53        format(' Geoid fit solution vector (all units m): '
     .    /,1x,f9.3,' (bias)',
     .    /,1x,2e14.6,' (slopes in E and N directions)',
     .    /,1x,3e14.6,' (second order terms)',
     .    /,1x,4e14.6,' (third order terms)')
          if (itrend.ge.5) write(*,54) (tcoef(i),i=1,nt)
54        format(' Geoid datum shift fit solution vector:',
     .    /,1x,3f10.3,2x,e12.6)
          close(21)
        endif
      endif
c
c  input data file and format
c  ---------------------------
c
      write(*,15)
15    format(
     .' Do you wish to input data points from a file? (T/F)')
      if (ldos) write(*,*) '-> '
      read(*,*) linp
17    format(a1)
      if (linp) then
18      write(*,16) dfile
16      format(' Enter file name: (CR=',a9,') ')
        if (ldos) write(*,*) '-> '
        read(*,1) file
        if (file(1:1).ne.' ') dfile=file
        open(20,file=dfile,form='formatted',status='old',err=19)
        goto 20 
19        write(*,*) '*** file could not be opened'
          goto 18 
20      continue
      else
        linp = .false.
      endif  
c
      write(*,21) ofile
21    format(' Enter file name for output: (CR=',a9,') ')
      if (ldos) write(*,*) '-> '
      read(*,1) file
      if (file(1:1).ne.' ') ofile=file
      open(30,file=ofile,form='formatted',status='unknown')
c
      if (itask.le.1.and.linp) write(*,23)
23    format(
     .' Type of input: 1 = statno, lat, lon (degrees)'/
     .'                2 = statno, lat, lon (deg,min,sec)'/
     .'                3 = statno, X, Y, Z (meter)'/
     .'                4 = statno, N, E (UTM, meter)')
      if (itask.le.1.and.(.not.linp)) write(*,231)
231   format(
     .' Type of input: 1 = lat, lon (degrees)'/
     .'                2 = lat, lon (deg,min,sec)'/
     .'                3 = X, Y, Z (meter)'/
     .'                4 = N, E (UTM, meter)')
      if (itask.gt.1.and.linp) write(*,24)
24    format(
     .' Type of input: 1 = statno, lat, lon, height (degrees)'/
     .'                2 = statno, lat, lon, height (deg,min,sec)'/
     .'                3 = statno, X, Y, Z (meter)'/
     .'                4 = statno, N, E, height (meter)')
      if (itask.gt.1.and.(.not.linp)) write(*,241)
241   format(
     .' Type of input: 1 = lat, lon, height (degrees)'/
     .'                2 = lat, lon, height (deg,min,sec)'/
     .'                3 = X, Y, Z (meter)'/
     .'                4 = N, E, height (meter)')
      if (ldos) write(*,*) '-> '
      read(*,*) itypi
      if (itypi.eq.3) then
        write(*,25)
25      format(
     .  ' Select datumshift for XYZ input: 1 = None (WGS84 to WGS84)'/
     .  '                                  2 = WGS84 to ED50'/
     .  '                                  3 = WGS84 to NAD27'/
     .  '                                  4 = WGS84 to Qornoq'/
     .  '                                  5 = NWL9D to WGS84')
        if (ldos) write(*,*) '-> '
        read(*,*) idati
        if (idati.lt.1.or.idati.gt.6) stop 'wrong datum'
        ielli = idati
        if (idati.eq.4) ielli = 2
        if (idati.eq.5) ielli = 1
        if (idati.eq.6) then
          ielli = 5
          idati = 1
        endif
      elseif (itypi.eq.4) then
        write(*,26) 
26      format(
     .  ' Enter UTM ellipsoid number 1 = WGS84'/
     .  '                            2 = Hayford (ED50)'/
     .  '                            3 = Clarke (NAD27)'/
     .  '                            4 = Bessel'/
     .  ' and UTM zone number (99 for RT38-approximate): ')
        if (ldos) write(*,*) '-> '
        read(*,*) iell,izone
        if (iell.lt.1.or.iell.gt.4.or.izone.lt.1.or.izone.gt.99)
     .  stop '*** illegal specification of UTM ellipsid and zone'
        call utmcon(iell,izone,sai)
      endif
c
c  output data format
c  -------------------
c
      itypo = 1
      if (itypi.eq.2) itypo = 2
      if (itask.eq.0) then
        write(*,27)
27      format( 
     .  ' Wanted output: 1 = statno, lat, lon (degrees)'/
     .  '                2 = statno, lat, lon (deg,min,sec)'/
     .  '                3 = statno, X, Y, Z (meter)'/
     .  '                4 = statno, N, E (UTM, meter)')
        if (ldos) write(*,*) '-> '
        if (itask.eq.0) read(*,*) itypo
        if (itypo.eq.3.and.itypi.ne.3) write(*,271)
271     format(' - NB: this option requires additional height',
     .  ' given in input -')
      else
        write(*,272) 
272     format(' - output coordinates are geographic degrees -')
      endif
      if (itypo.eq.3) then
        write(*,28)
28      format(
     .  ' Select datumshift for XYZ output: 1 = None (WGS84 to WGS84)'/         
     .  '                                   2 = ED50 to WGS84'/
     .  '                                   3 = NAD27 to WGS84'/
     .  '                                   4 = Qornoq to WGS84'/
     .  '                                   5 = WGS84 to NWL9D')
        if (ldos) write(*,*) '-> '
        read(*,*) idato
        iello = idato
        if (idato.eq.4) iello = 2
        if (idato.eq.5) iello = 1
      elseif (itypo.eq.4) then
        write(*,29) 
29      format(
     .  ' Enter UTM ellipsoid number 1 = WGS84'/
     .  '                            2 = Hayford (ED50)'/
     .  '                            3 = Clarke (NAD27)'/
     .  '                            4 = Bessel'/
     .  ' and UTM zone number (99 for RT38-approximate): ')
        if (ldos) write(*,*) '-> '
        read(*,*) iell,izone
        if (iell.lt.1.or.iell.gt.4.or.izone.lt.1.or.izone.gt.99)
     .  stop '*** illegal specification of UTM ellipsoid and zone'
        call utmcon(iell,izone,sao)
      endif
c
      lih = (itask.gt.1.or.itask.eq.0.and.itypo.eq.3)
c
c  reading loop for individual points
c  ----------------------------------
c
      lfirst = .true.
      n = 0
      nrej = 0  
c
100   n = n+1
      id = n
      h = 0.0
      goto (110,120,130,140),itypi
c
c  degrees
c
110   if (linp) then
        if (lih) read(20,*,end=200) id,rfi,rla,h
        if (.not.lih) read(20,*,end=200) id,rfi,rla
      else
        if (lih) then
          write(*,111)
111       format(' Enter: rfi, rla, h (x=exit)')
          if (ldos) write(*,*) '-> '
          read(*,*,err=200) rfi,rla,h
        else
          write(*,112)
112       format(' Enter: rfi, rla (x=exit)')
          if (ldos) write(*,*) '-> '
          read(*,*,err=200) rfi,rla
        endif
      endif
      goto 150
c
c  deg, min, sec
c
120   if (linp) then
        if (lih) read(20,*,end=200)
     .  id,ifi,mfi,sfi,ila,mla,sla,h
        if (.not.lih) read(20,*,end=200) 
     .  id,ifi,mfi,sfi,ila,mla,sla
      else
        if (lih) then
          write(*,121)
121       format(' Enter: lat, lon, h (deg min sec, x=exit)')
          if (ldos) write(*,*) '-> '
          read(*,*,err=200) ifi,mfi,sfi,ila,mla,sla,h
        else
          write(*,122)
122       format(' Enter: lat, lon (deg min sec, x=exit)')
          if (ldos) write(*,*) '-> '
          read(*,*,err=200) ifi,mfi,sfi,ila,mla,sla
        endif
      endif
      ii = 1
      if (ifi.lt.0.or.mfi.lt.0.or.sfi.lt.0) ii = -1
      rfi = ii*(abs(ifi) + abs(mfi)/60.0 + abs(sfi)/3600)
      ii = 1
      if (ila.lt.0.or.mla.lt.0.or.sla.lt.0) ii = -1
      rla = ii*(abs(ila) + abs(mla)/60.0 + abs(sla)/3600)
      goto 150
c
c  X Y Z, convert to geographic
c
130   if (linp) then
        read(20,*,end=200) id,x,y,z       
      else
        write(*,131)
131     format(' Enter: X Y Z (meter, x=exit)')
        if (ldos) write(*,*) '-> '
        read(*,*,err=200) x,y,z
      endif
      call ctg(x,y,z,rfi,rla,h,ielli,idati,.true.,.true.)
      rfi = rfi*radeg
      rla = rla*radeg
      goto 150
c
c  UTM, convert to geographic
c
140   if (linp) then
        if (lih) read(20,*,end=200) id,rn,re,h
        if (.not.lih) read(20,*,end=200) id,rn,re
      else
        if (lih) then
          write(*,142)
142       format(' Enter: N, E, h (x=exit)')
          if (ldos) write(*,*) '-> '
          read(*,*,err=200) rn,re,h
        else
          write(*,143) 
143       format(' Enter: N, E (x=exit)')
          if (ldos) write(*,*) '-> '
          read(*,*,err=200) rn,re
        endif
      endif
      call utg(rn,re,rfi,rla,sai,.true.,.true.)
      rfi = rfi*radeg
      rla = rla*radeg
      goto 150
c
c  interpolation of geoid value at point rfi,rla
c  ---------------------------------------------
c
150   if (itask.eq.0) goto 169
      rgeoid = geoid1(rfi,rla,lok)
      if (lds) call geoitr(ids,rfi,rla,rgeoid)
c
c  do deflections by taking geoid differences
c
      if (ldfv) then
        rn = geoid1(rfi+0.0001,rla,lok)
        if (lds) call geoitr(ids,rfi+0.0001,rla,rn)
        re = geoid1(rfi,rla+0.0001,lok)
        if (lds) call geoitr(ids,rfi,rla+0.0001,re)
        xi = -(rn-rgeoid)/0.0001*1.85498 
        eta = -(re-rgeoid)/0.0001*1.85498/cos(rfi/radeg)
      endif
      if (.not.lok) goto 100
c
c  convert to/from ellipsoidal height
c
      if (itask.eq.1) then
        h = rgeoid
      elseif (itask.eq.2) then
        h = h-rgeoid
        if (ltrend) then
          dn = trend(rfi,rla,itrend,tcoef)
          write(*,152) dn
152       format('   Geoid trend correction for the following point =',
     .    f9.3)
          h = h-dn
        endif
      else
        h = h+rgeoid
      endif
c
c  output results
c  --------------
c
169   noh = 0
      if (itask.gt.0.or.itypi.eq.3) then
        noh = 1
        hh(1) = h 
      endif
      if (ldfv) then
        noh = 2
        hh(1) = xi
        hh(2) = eta
      endif
c     
      goto (170,175,180,185),itypo
c 
c  geographic degrees
c
170   write(*,171) id,rfi,rla,(hh(j),j=1,noh)
      write(30,171) id,rfi,rla,(hh(j),j=1,noh)
171   format(' ',i8,1x,2f13.8,2f10.3)
      goto 100
c
c  geog deg min sec
c
175   lneg = (rfi.lt.0)
      ifi = abs(rfi)
      mfi = (abs(rfi)-ifi)*60.d0
      sfi = (abs(rfi)-(ifi+mfi/60.d0))*3600.d0
      if (lneg) ifi = -ifi
      if (lneg.and.ifi.eq.0) then
        if (lneg) mfi = -mfi
        if (mfi.eq.0) sfi = -sfi 
      endif
      lneg = (rla.lt.0)
      ila = abs(rla)
      mla = (abs(rla)-ila)*60.d0
      sla = (abs(rla)-(ila+mla/60.d0))*3600.d0
      if (lneg) ila = -ila
      if (lneg.and.ila.eq.0) then
        lzero = .true.
        if (lneg) mla = -mla
        if (mla.eq.0) sla = -sla 
      endif
c
      write(*,177) id,ifi,mfi,sfi,ila,mla,sla,(hh(j),j=1,noh)
      write(30,177) id,ifi,mfi,sfi,ila,mla,sla,(hh(j),j=1,noh)
177   format(' ',i8,1x,2(i6,i3,f9.5),2f10.3)
      goto 100
c
c  X Y Z coordinates
c
180   call ctg(rfi/radeg,rla/radeg,h,x,y,z,iello,idato,.false.,.true.)
      write(*,181) id,x,y,z
      write(30,181) id,x,y,z
181   format(' ',i8,1x,3f13.3)
      goto 100
c
c  utm coordinates
c
185   call utg(rfi/radeg,rla/radeg,rn,re,sao,.false.,.true.)   
      write(*,186) id,rn,re,(hh(j),j=1,noh)
      write(30,186) id,rn,re,(hh(j),j=1,noh)
186   format(' ',i8,1x,2f13.3,2f10.3)
      goto 100
c
c  exit of point reading loop
c  --------------------------
c
200   if (itask.gt.0) write(*,201) n-1-nrej,nrej
      if (itask.eq.0) write(*,202) n-1
201   format(' --- number of points interpolated:',i4,', rejected:',i4)
202   format(' --- number of points transformed:',i4)
      write(*,203) ofile
203   format(' --- outputfile: ',a60)
      end
c
      double precision function geoid1(rfi,rla,lok)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                          g e o i d 1
c
c  this subroutine interpolates geoid value from the geoid grid file,
c  which must be opened in advance, and connected as unit 10.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      logical lutm,lok
      real*4 grec(17)
      dimension sag(22)
      common
     *  /gpar/ itask,itask1,rfi1,rla1,dfi,dla,nn,ne,lutm,lutm1,
     *  sag,radeg,id,nrej
c
      lok = .false.
      xfi = rfi
      xla = rla
      if (lutm)
     .call utg(rfi/radeg,rla/radeg,xfi,xla,sag,.false.,.true.)
      ri = (xfi-rfi1)/dfi+1.0
      rj = (xla-rla1)/dla+1.0
      if (ri.lt.0.999.or.ri.gt.nn+0.001.or.rj.lt.0.999.or.rj.gt.ne
     .+0.001) then
        write(*,*) '*** Point no.',id,' outside grid, ',
     .  'nothing output ***'
        nrej = nrej+1
        return
      endif
      ii = ri
      if (ii.eq.0) ii = 1
      if (ii.ge.nn) ii = nn-1
      jj = rj
      if (jj.eq.0) jj = 1
      if (jj.ge.ne) jj = ne-1
      ri = ri - ii
      rj = rj - jj
c
c  find record number and four surrounding heights
c
      jrec = (ne-1)/16+1
      irec = (nn-ii)*jrec+2
      jjrec = (jj-1)/16
      iirec = irec+jjrec
      kk = jj - jjrec*16
      kk1 = min(kk+1,16)
      read(10,rec=iirec) (grec(k),k=1,kk1)
      if (kk.eq.16) read(10,rec=iirec+1) grec(17)
      g1 = grec(kk)
      g2 = grec(kk+1)
      read(10,rec=iirec-jrec) (grec(k),k=1,kk1)
      if (kk.eq.16) read(10,rec=iirec-jrec+1) grec(17)
      g3 = grec(kk)
      g4 = grec(kk+1)
      if (g1.ge.9999.or.g2.ge.9999.or.g3.ge.9999.or.g4.ge.9999) then
        write(*,*) '*** Point ',id,' skipped, interpolation ',
     .  'involves unknown (9999) value ***'
        nrej = nrej+1
        return
      endif
c
c  linear interpolation
c
      geoid1 = (1-ri)*(1-rj)*g1 + (1-ri)*rj*g2 +
     .ri*(1-rj)*g3 + ri*rj*g4
      lok = .true.
      return
      end
c
      subroutine geoitr(ids,rfi,rla,geoid)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      g e o i t r
c
c  This subroutine will transform a geoid height in WGS84 to a geoid
c  height in the local datum 'ids', see ctg subroutine. 
c
c  rf, oct 91
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      radeg = 180.0/3.1415926536d0
      rn = rfi/radeg
      re = rla/radeg
      h = geoid
      iell = ids
      if (ids.eq.4) iell = 2
c
      call ctg(rn,re,h,x,y,z,1,1,.false.,.true.)
      call ctg(x,y,z,rn,re,h,iell,ids,.true.,.false.)
c
      geoid = h
      return
      end
c
      subroutine ctg(x1,x2,x3,y1,y2,y3,iell,idatum,ldir,lcheck) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                           c t g
c
c  Dual autochecking transformation procedure cartesian to
c  geographic coordinates or reverse. The cartesian system
c  is righthanded with z-axis towards the mean pole and
c  x-axis intersecting the Greenwich meridian. The geographic
c  coordinates are latitude (1), longitude (2) and height
c  above the ellipsoid (3).
c  The subroutine may also do a datum transformation, with
c  the geographical coordinates assumed to be in a specified
c  system in either input or output (XYZ always assumed WGS84)
c
c  parameters:
c
c  ldir = true:  cartesian -> geographic
C  ldir = false: geographic -> cartesian
c
c  x1, x2, x3:   input coordinates
c                (ldir = .true.: X, Y, Z in meter,
c                 ldir = .false.: lat, lon in radians, H in m)
c
c  y1, y2, y3:   output coordinates
c
c  iell:         ellipsoid used in transformation  
c                1  grs80  ellipsoid
c                2  hayford ellipsoid
c                3  nad27 ellipsoid
c                4  bessel ellipsoid
c                5  NWL9D ellipsoid
c
c  idatum:       datum transformation 
c                1  none, XYZ and fi,la,H in same system
c                2  XYZ in WGS84, fi,la,H in ED50
c                3       do     , fi,la,H in NAD27
c                4       do     , fi,la,h in Qornoq Greenland datum
c                5  XYZ in NWL9D, fi,la,h in WGS84
c
c  lcheck:       if true check by reverse transformation
c
c  RF, august 1981.
c  Procedure based on older KP procedures 'trctg' and 'trgtc'
c  Datum shift constants from 'setddatcon' algol procedure
c  Fortran version by RF dec 1990
c  (c) Rene Forsberg, Kort- og Matrikelstyrelsen, Denmark
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h, o-z)
      logical ldir, lcheck
c
c  ellipsoid and datum shift parameters
c   
      dimension aa(5),ff(5),dx(5),dy(5),dz(5),sc(5),rot(5)
c
c       WGS84        Hayford      Clarke      Bessel          NWL9D
      data aa /
     .6378137.0d0, 6378388.0d0, 6378206.4d0, 6377397.155d0, 6378145.0d0/
      data ff / 
     .298.2572236d0, 297.0d0, 294.9786982d0, 299.153d0,     298.25d0/
c
c  datumshifts:
c  2: ED50->WGS84   Ref: Report of an investigation into the use of doppler
c                   satellite positioning to provide coordinates on ED50 in the
c                   area of the North Sea. Prof. Papers new series no 30, 1981.
c                   Followed by DMA NWL9D to WGS84
c  3: NAD27->WGS84  Ref: DMA World Geodetic System 1984 report (CONUS value)
c  4: Qornoq->WGS84 Ref: Anderle 1976, table 7 converted to WGS84, as KMS
c  5: WGS84->NWL9D  Ref: Definition of WGS 84, DMA; Letter from DMA june 1985
c
c               (1)         (2)      (3)     (4)          (5) 
c
      data dx / 0.d0,    -89.5d0,  -8.d0, 163.511d0,      0.d0/  
      data dy / 0.d0,    -93.8d0, 160.d0, 127.553d0,      0.d0/ 
      data dz / 0.d0,   -123.1d0, 176.d0,-159.789d0,    -4.5d0/
      data sc / 0.d0,     1.2d-6,   0.d0,   -0.6d-6,    0.6d-6/
      data rot/ 0.d0,  .75629d-6,   0.d0,-3.9464d-6, 3.9464d-6/
c
      if (iell.lt.1.or.iell.gt.5) stop 'undefined ellipsoid in ctg'
      if (idatum.lt.0.or.idatum.gt.5) stop 'undefined datum in ctg'
      a = aa(iell)
      f = 1/ff(iell)
      f1 = 1-f
      e2 = f*(2-f)  
      em2 = 1.0/f1**2-1 
      c = a/f1
      x1in = x1
      x2in = x2
      x3in = x3  
c
      if (ldir) then
        j = 1
        k = 3
        l = 1
      else
        j = 3  
        k = 1
        l = -1
      endif
      do 50 i = j,k,l
        goto (10,20,30),i
c
c  case 1 - cart->geo 
c
10      if (idatum.gt.1) then
          xx1 = x1in
          x1in = x1in - (dx(idatum)+sc(idatum)*x1in+rot(idatum)*x2in)
          x2in = x2in - (dy(idatum)+sc(idatum)*x2in-rot(idatum)*xx1)       
          x3in = x3in - (dz(idatum)+sc(idatum)*x3in)
        endif  
        p = sqrt(x1in**2+x2in**2)
        rla = atan2(x2in, x1in)
        fi = 0
        w = 0  
        nc = 0
c  repeat loop
12        fi1 = fi  
          fi = atan2(x3in+w, p)
          h = sqrt((x3in+w)**2 + p**2)
          sinfi = sin(fi)
          cos2fi = 1-sinfi**2
          w = c/sqrt(1+em2*cos2fi)*e2*sinfi
          nc = nc+1
          if (abs(fi-fi1).gt.1.0e-9.and.nc.lt.10) goto 12
        x1in = fi
        x2in = rla
        x3in = h - c/sqrt(1+em2*cos2fi)
        goto 50
c 
c  case 2 - dump of result 
c
20      y1 = x1in  
        y2 = x2in  
        y3 = x3in  
        if (.not.lcheck) return
        goto 50
c    
c  case 3 - geo->cart
c
30      cosfi = cos(x1in)        
        sinfi = sin(x1in)
        rn = c/sqrt(1+em2*cosfi**2)
        w = (rn + x3in)*cosfi  
        x1in = w*cos(x2in)  
        x2in = w*sin(x2in)  
        x3in = (rn - rn*e2 + x3in)*sinfi  
        if (idatum.gt.1) then
          xx1 = x1in
          x1in = x1in + (dx(idatum)+sc(idatum)*x1in+rot(idatum)*x2in)
          x2in = x2in + (dy(idatum)+sc(idatum)*x2in-rot(idatum)*xx1)
          x3in = x3in + (dz(idatum)+sc(idatum)*x3in)
        endif  
50    continue   
c
c  tolerance check of results - 2 mm on earth
c
      if (ldir) then
        if (abs(x1-x1in).gt.0.002.or.abs(x2-x2in).gt.0.002.or.
     .  abs(x3-x3in).gt.0.002) write(*,60) x1,x2,x3,
     .  x1-x1in,x2-x2in,x3-x3in
60      format(' *** cart to geo transformation check - error too big:'/
     .  ' X Y Z = ',3f10.1,'  dX dY dZ = ',3f12.3) 
      else
        if (abs(x1-x1in)*6371000.gt.0.002.or.
     .  abs(x2-x2in)*6371000*cosfi.gt.0.002.or.
     .  abs(x3-x3in).gt.0.002) write(*,61) x1*57.29578,x2*57.29578,x3,
     .  (x1-x1in)*6371000,(x2-x2in)*6371000,x3-x3in
61      format(' *** geo to cart transformation check - error too big:'
     .  /' fi la h = ',3f9.4,'  dX dY dZ = ',3f12.3) 
      endif
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
c rc fortran version alp/rf oct 86, last updated dec 1990
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension sa(22)
      double precision n,m
      dimension aa(4),ff(4)
c
c  ellipsoid parameters
c                  WGS84      Hayford      Clarke       Bessel
c
      data aa / 6378137.0d0, 6378388.0d0, 6378206.4d0, 6377397.155d0/
      data ff / 298.2572236d0, 297.0d0, 294.9786982d0, 299.153d0/
      radeg = 180/3.1415926536d0
c
      if (isys.lt.1.or.isys.gt.4) stop 'illegal ellipsoid in utmcon'
      if (izone.lt.1.or.izone.gt.99) stop 'illegal UTM zone in utmcon'
      if (izone.eq.99) isys = 4
      a = aa(isys)
      f = 1.d0/ff(isys)
c
      eastpr = 500000.d0
      dm = 4.0d-4
      if (izone.eq.99) dm = 0.0
      if (izone.eq.99) eastpr = 1500000.d0
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
      radeg = 180/3.1415926536d0
c
      n = rn
      e = re
c
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
      subroutine setov(rfi,rla,obs,itrend,nt,ov)
      implicit double precision(a-h,o-z)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      s e t o t v
c
c  sets observation equation coefficients for detrending
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension ov(11)
      common /tpar/ cosfi0,rfi0,rla0
      radeg = 180/3.1415926536d0
      if (itrend.ge.5) goto 50
c
c  just mean value or polynomials, length unit m
c
      dx = (rla-rla0)*cosfi0*111111
      dy = (rfi-rfi0)*111111
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
50    nt = 3
      if (itrend.eq.6) nt = 4
      cosf = cos(rfi/radeg)
      ov(1) = cosf*cos(rla/radeg)
      ov(2) = cosf*sin(rla/radeg)
      ov(3) = sin(rfi/radeg)
      if (itrend.eq.6) ov(4) = 6371000.0
c
90    ov(nt+1) = obs
      return
      end
c
      function trend(rfi,rla,itrend,tsol)
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
      call setov(rfi,rla,rr,itrend,nt,ov)
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
      dimension c(*)
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
          c(i) = ci/c(ic+nc)
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
        c(i) = c(i)/c(ic)
        do 70 np = nc-1,1,-1
          ir=ir-1
          ic=ic-1
          c(ir) = c(ir) - c(i)*c(ic)
70      continue
        i = i-1
80    continue
      return
      end
