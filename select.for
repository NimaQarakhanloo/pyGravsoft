      program selec
c $Id: select.for 188 2008-08-18 08:33:51Z cct $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                         s e l e c 
c
c  program for selection, averaging or reformatting data.
c  the program may select one data point per pixel, closest to the
c  knots of a given grid, or produce a straight average of all values
c  in the grid cell. the program may also be used only to reformat data
c  into standard format (id, lat, lon, height, data ..). it may also
c  produce a primitive coverage plot.
c
c  input:
c
c  inputfile,
c  outputfile,
c  mode, iang, n,
c     (dataline code spec plus format - only for iang<=0)
c  (rfi1, rfi2, rla1, rla2, dfi, dla (degrees or m) - only for mode>0)
c     (iell, izone - only when limits are northing and easting in meter)
c     (rejection level - for mode 5 and 7 only)
c     (window min., max. latitude and longitude - for mode 6 and 7 only)
c
c  'mode' determines the function of the program:
c   0: reformat all data
c   1: select data closest to knots (if dfi=0 all data in area selected)
c   2: average all data in grid cell and output in list format
c   3: do, output in grid format with 9999.0 in cells without data
c   4: primitive screen plot
c   5: as 1, but data with estimated error larger than 'rejlev' are rejected
c   6: as 1, but a window inside the area is excluded
c   7: as 1, but select minimal value in cell
c   8: as 1, but with random, normal distributed noise added. noise standard
c      deviation must be input after grid specification
c
c  'iang' is a code for input coordinates and format: 
c      1: no, lat/lon in degrees, height, data1, data2  ..
c      2: do, with lat/lon in degrees and minutes 
c      3: do, with lat/lon in deg, min, sec 
c      4: altimetry format (n=1 gives ssh only, 2 ssh+error, 3 ssh+error+t),
c      5: binary format,
c      6: gravity data in kms 80-char format
c         (n=1 fa, n=2 fa+err, n=3 g,fa,ba, n=-3: fa,err,tc)
c      7: grid format (if n=0 a height grid of integers assumed)=
c      if 'iang' = 0 a data line specification must be input (see below),
c      if 'iang' is negative a standard data line format must be input (do)
c
c  'n' is the number of data following no, lat, lon, elev (max 3) in stan-
c      dard format, or number of wanted data in data line format.
c      for averaging only data number 'n' is used and output.
c
c  'rfi1, rfi2, rla1, rla2' (degrees) specifies boundaries of wanted grid
c
c  'dfi, dla' (degrees) specifies wanted grid spacing, if zero all data
c      within area is selected (mode 1), or a default plot is made (mode 3)
c
c  utm option: if rfi1 or rfi2 > 100 then the wanted grid is a utm grid,
c  and ellipsoid number (1: wgs 84, 2: ed50, 3: nad27, ..) and utm zone
c  must be input. If zone is negative polar stereographic assumed.
c  north polar stereographic: -75 0 means true scale 75N and ref lon 0. 
c  if coordinates in file are > 90 or 180 the file is assumed to contain 
c  northings and eastings in same system.
c  
c
c  input data is normally assumed to free format.
c  if 'iang' is negative iang = abs(iang) and a data format in brackets
c  without quotation marks must be specified, e.g. (i4,2f10.6,f8.2,2f7.2)
c
c  dataline option:
c  if 'iang' = 0 then 'n' is the number of data wanted in output.
c  data sequence in a line (free format) is specified by
c  codes 1: station number, 2: latitude (21: lat min, 22: lat sec),
c  3: longitude (31: lon min, 32: lon sec), 4: height, 5, 6, .. data.
c  the data codes must be specified in i3 format, followed by lfree and
c  wanted data, e.g.
c
c    1  2 21  3 31  4  5  6
c    t  6  5
c    (plus format if lfree.eq.false)
c
c  which will read a line with station number, coordinates in degrees
c  and minutes, height followed by two data values, but output data in
c  reverse order. if a station number is not specified (code 1) a sequential
c  number will be given. if heights are not given 0.0 will be used.
c  'lfree' is true if free (*) format may be used. otherwise the
c  format (all reals) must be given on the next line.
c
c  map option: a primitive contour map of
c  data values at pixels will be output for check of station coverage
c  rather than the output list of coordinates. if dfi = 0 a grid will
c  be selected with default options in this case.
c  (the map option should only be attempted for reasonably small grids)
c
c  rejection option: 
c  in mode 5 or 7 the rejection level must be given. if a negative rejection
c  level is given, priority is given in the pixel search to smallest stan-
c  dard deviation (i.e., a low std.dev. is selected rather than the closest
c  point). in standard format the std.dev. must be the second data (n=2).
c
c  programmer: rene forsberg, januar 1985
c  updated jan 89, university of new south wales vax, rf
c  updated june 22, 1989 by cct.
c  last updated jan 16, 1992, rf, mv, hd
c  -            may 25, 1994, rf
c  -            dec 5, 1994, rf
c  -            jan 20, 2001, rf (utm input)
c  -            mar 19, 2002, rf (minimal pixel values)
c  -            dec 30, 2003, rf (mingw update)
c  -            aug 7, 2007, rf (polar stereographic grids)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-h,o-z)
      parameter(iamax=2030000)
      dimension is(iamax)
      real*4 rfi(iamax),rla(iamax)
      real*4 rh(iamax),g(iamax,3),d(iamax)
      dimension sg(9),rdat(40),ssmin(9),ssmax(9),glab(6),
     *g1(9),stdvn(9),s1(9),s2(9),ss1(9),ss2(9),sa(22),sag(22),
     *grow(22000)
      real*4 qfi,qla,qh,qg
      logical all,ldatal,lstd,lfree,lform,lref,lmap,lmean,lrej,lwind,
     *lnoise,lrej1,lutm,lutmo,lgutm,lgrid,lminus,lmin
      integer*4 idat(40),idno(40),ip(40)
      character ifile*72,ofile*72,dfmt*72,hc*1,tccode*1
      data s1/9*0.d0/,s2/9*0.d0/,ss1/9*0.d0/,ss2/9*0.d0/,
     *ssmin/9*9.9d9/,ssmax/9*-9.9d9/
c
      radeg = 180/3.1415926536d0
      lutmo = .false.
c
      write(*,1)
1     format(' input: inputfile and outputfile (two lines): ')
      read(*,2) ifile
      read(*,2) ofile
2     format(a72)
c
c  constants for grid
c
      write(*,20)
20    format(' input: MODE (0:reformat, 1:select, 2:mean, 3:grid, ',
     .'4:plot,',
     ./8x,'      5:sel&rej, 6:sel&wndw, 7:sel_min, 8:sel&noise)',
     ./8x,'IANG (1:deg, 2:dm, 3:dms, 4:alt, 5:bin, 6:80char, 7:grid,'
     ./8x,'      neg:fmt, 0:dline)'
     ./8x,'NDATA'/)
      read(*,*) mode, iang, n
c
      lref = (mode.eq.0)
      lmean = (mode.eq.2.or.mode.eq.3)
      lrej = (mode.eq.5)
      lwind = (mode.eq.6)
      lmin = (mode.eq.7)
      lnoise = (mode.eq.8)
      lminus = (n.lt.0)
c
      n = abs(n)
      if (lnoise) mode = 1
      if (lrej.or.lwind.or.lmin) mode = 1
      if (lrej.and.n.lt.2.and.iang.ne.6)
     *stop 'std.dev. must be data number 2'
      lmap = (mode.eq.4)
c
      lform = (iang.lt.0)
      lfree = (.not.lform)
      if (lform) iang = -iang
      if (lform) write(*,821)
      if (lform) read(*,822) dfmt
      ldatal = (iang.eq.0)
      lstd = (.not.ldatal)
      lgrid = (iang.eq.7)
      if (lgrid.and.n.gt.1) n=1
      if (n.gt.3.and.(.not.lref)) stop 'n limited to max 3'
      if (n.gt.9) stop 'n too large'
c
      if (iang.ne.5) open(10, file=ifile, status='old')
      if (iang.eq.5) open(10, file=ifile, status='old',
     .form='unformatted')
      open(20, file=ofile, status='unknown')
c
      if (lstd) goto 90
c
c  read dataline specification and check if all info is there
c
      write(*,80)
80    format(' input: data line codes in format i3'/
     .' (1:no, 2/21/22:lat, 3/31/32:lon, 4:h, 5,6..:data, -1:end)'/)
      ndl = 1
      read(*,81) idat
81    format(40i3)
      ndl = 1
811   if (idat(ndl).le.0) goto 82
      ndl = ndl+1
      goto 811
82    ndl = ndl-1
      write(*,821) (idat(j),j=1,ndl)
821   format(' ',24i3)
      write(*,824)
824   format(' input: lfree (t:unf, f:form) and wanted data codes: ')
      read(*,*) lfree, (idno(i),i=1,n)
      lform = (.not.lfree)
      if (lform) write(*,825)
825   format(' input format of data line - only real formats allowed')
      if (lform) read(*,822) dfmt
822   format(a72)
c
      ino = 0
      ilat = 0
      ilatm = 0
      ilats = 0
      ilon = 0
      ilonm = 0
      ilons = 0
      ih = 0
      do 85 j= 1, n
85    ip(j) = 0
c
      do 83 i = 1, ndl
        if (idat(i).eq.1) ino = i
        if (idat(i).eq.2) ilat = i
        if (idat(i).eq.21) ilatm = i
        if (idat(i).eq.22) ilats = i
        if (idat(i).eq.3) ilon = i
        if (idat(i).eq.31) ilonm = i
        if (idat(i).eq.32) ilons = i
        if (idat(i).eq.4) ih = i
        do 84 j = 1, n
84      if (idat(i).eq.idno(j)) ip(j) = i
83    continue
c
      if (ilat.eq.0) stop 'no latitude info in file'
      if (ilon.eq.0) stop 'no longitude info in file'
      do 86 j = 1, n
86    if (ip(j).eq.0) stop 'wanted data code not in file'
c
c  read grid and area specification
c  --------------------------------
c
90    if (lref) goto 199 
      write(*,*) 'input: FI1,FI2,LA1,LA2,DFI,DLA (deg or m) '
      read(*,*)  rfi1, rfi2, rla1, rla2, dfi, dla
      if (rfi1.ge.rfi2.or.rla1.ge.rla2) stop 'wrong area specification'
      all = .false.
      if (dfi.eq.0) all = .true.
c
c  check for utm 
c
      lutm = (abs(rfi1).gt.100.or.abs(rfi2).gt.100)
      if (lutm) then
        write(*,*) 'input: utm ellipsoid no and zone: '
        read(*,*) iell,izone
        cosfi = 1.0
        call utmcon(iell,izone,sa)
        call utg(rfi2,rla1,rf1,rl1,sa,.true.,.true.)
        call utg(rfi2,rla2,rf2,rl2,sa,.true.,.true.)
        call utg(rfi1,rla1,rf3,rl3,sa,.true.,.true.)
        call utg(rfi1,rla2,rf4,rl4,sa,.true.,.true.)
        write(*,91) rf1*radeg,rl1*radeg,rf2*radeg,rl2*radeg
     *  ,rf3*radeg,rl3*radeg,rf4*radeg,rl4*radeg
91      format(' corner lat and lon of wanted utm area: '
     *  ,4(/,' ',2f12.5))
      else
        cosfi = cos((rfi2+rfi1)/2/radeg)
      endif
c
      if (lnoise) then
        write(*,*) 'input noise stdv. (1 to n)'
        read(*,*) (stdvn(k),k=1,n)
      endif
      if (lrej) then
        write(*,*) 
     .  'input: rejection level (negative for stddev priority)'
        read(*,*) rejlev
        lrej1 = (rejlev.lt.0)
        if (lrej1) rejlev = abs(rejlev)
      endif
      if (lwind) then
        write(*,*) 
     .  'input: latmin,latmax,lonmin,lonmax for inner window: '
        read(*,*) rfi1i,rfi2i,rla1i,rla2i
      endif
      if (all.and.lmean) stop 'grid spacing must be specified'
c
c  default map parameters
c
      if (.not.(all.and.lmap)) goto 12
      if (rfi1.eq.rfi2.and.rla1.eq.rla2) stop 'area limits wrong'
      all = .false.
      nn = 20
      ne = nn*(rla2-rla1)*cosfi/(rfi2-rfi1)*2.0
      if (ne.gt.78) ne = 78
      if (ne.lt.2) ne = 2
      dfi = (rfi2-rfi1)/(nn-1)
      dla = (rla2-rla1)/(ne-1)
c
12    nn = 0
      ne = 0
      if (.not.all) nn = (rfi2-rfi1)/dfi+1.5
      if (.not.all) ne = (rla2-rla1)/dla+1.5
      nnn = nn*ne
      if (nn*ne.gt.iamax) write(*,11) nn*ne
11    format(' **** specified grid too large ****',i8)
      if (nn*ne.gt.iamax) stop
      if (n.gt.3) n=3
      rfi0 = rfi1-dfi/2
      rla0 = rla1-dla/2
c
c  initialization of arrays storing grid points
c  --------------------------------------------
c
199   nstat = 0
      if (lref) goto 101
      nsel = 0
      do 10 i=1,nnn
        is(i) = 0
        d(i) = 9.9e33
        if (.not.lmean) goto 10
        rfi(i) = 0
        rla(i) = 0
        rh(i) = 0
        d(i) = 0
        g(i,1) = 0
        g(i,2) = 0
        g(i,3) = 9999.0
10    continue
c
      rr2 = dfi**2 + (dla*cosfi)**2
101   nrej = 0
      istat = 0
      nn1 = 0
      nn2 = 0
c
      slatma = -9.9d9
      slatmi = 9.9d9
      slonma = -9.9d9
      slonmi = 9.9d9
c
      write(*,102)
102   format(/' --- S E L E C ---')
      write(*,*) 'output to: ',ofile
c
c  read label for grid data
c 
      if (lgrid) then
        read(10,*) glab 
        lgutm = (abs(glab(1)).gt.100.or.abs(glab(2)).gt.100)
        if (lgutm) read(10,*) iellg,izoneg
        if (lgutm.and.(iellg.lt.0.or.iellg.gt.4.or.
     *  izoneg.lt.1.or.izoneg.gt.99)) stop 'utm ell zone invalid'
        dfig = glab(5)
        dlag = glab(6)
        nng = (glab(2)-glab(1))/dfig + 1.5
        neg = (glab(4)-glab(3))/dlag + 1.5
        if (neg.gt.22000) then 
          write(*,*) neg
          stop 'e-w grid row too long'
        endif
        sh = 0.0
        ing = nng+1
        ieg = neg
        istat = 0
        if (lgutm) call utmcon(iellg,izoneg,sag)
        write(*,14) nng, neg, nng*neg
14      format(' selection from grid data, grid points:',2i7,i9)
      endif
c
c  skip label starting with '<' in 80-char format
c
      if (iang.eq.6) then
        read(10,15) hc
15      format(a1)
        if (hc.ne.'<') rewind(10)
      endif
c
c  scan through gravity data - main loop entry label 200 
c  -----------------------------------------------------
c
200   goto (205,210,220,230,231,233,235,237),iang+1
c
c  iang = 0: read user defined data line
c
205   if (lfree) read(10,*,end=900) (rdat(i),i=1,ndl)
      if (lform) read(10,fmt=dfmt,end=900) (rdat(i),i=1,ndl)
      if (ino.eq.0) istat = nstat+1
      if (ino.gt.0) istat = nint(rdat(ino))
      sfi = rdat(ilat)
      if (ilatm.eq.0) goto 201
        i = 1
        if (sfi.lt.0) i = -1
        sfi = abs(sfi) + rdat(ilatm)/60
        if (ilats.gt.0) sfi = sfi + rdat(ilats)/3600
        sfi = i*sfi
201   sla = rdat(ilon)
      if (ilonm.eq.0) goto 202
        i = 1
        if (sla.lt.0) i = -1
        sla = abs(sla) + rdat(ilonm)/60
        if (ilons.gt.0) sla = sla + rdat(ilons)/3600
        sla = i*sla
202   if (ih.eq.0) sh = 0.0
      if (ih.gt.0) sh = rdat(ih)
      do 203 k = 1, n
203   sg(k) = rdat(ip(k))
      goto 240
c
c  iang = 1,2,3: standard data line format
c
210   if (lfree) read(10,*,end=900) istat,sfi,sla,sh,(sg(k),k=1,n)
      if (lform) read(10,fmt=dfmt,end=900) istat,sfi,sla,sh,
     .(sg(k),k=1,n)
      goto 240
c
220   if (lfree)
     .read(10,*,end=900) istat,idfi,sfi,idla,sla,sh,(sg(k),k=1,n)
      if (lform)
     .read(10,fmt=dfmt,end=900) istat,idfi,sfi,idla,sla,sh,(sg(k),k=1,n)
      i = 1
      if (idfi.lt.0) i = -1
      j = 1
      if (idla.lt.0) j = -1
      sfi = i*(abs(idfi)+sfi/60.d0)
      sla = j*(abs(idla)+sla/60.d0)
      goto 240
c
230   if (lfree) read(10,*,end=900) istat,idfi,imfi,sfi,idla,imla,sla,
     .sh,(sg(k),k=1,n)
      if (lform) read(10,fmt=dfmt,end=900) istat,idfi,imfi,sfi,
     .idla,imla,sla,sh,(sg(k),k=1,n)
      i = 1
      if (idfi.lt.0) i = -1
      j = 1
      if (idla.lt.0) j = -1
      sfi = i*(abs(idfi)+imfi/60.0+sfi/3600)
      sla = j*(abs(idla)+imla/60.0+sla/3600)
      goto 240
c
c  iang = 4: altimeter data format
c
231   continue
c  ers-1 j2-format
c     read(10,*,end=900) irev,time,sfi,sla,sh,ssh,sigma
c     sh = sh+ssh
c  ers-1 t2-format
      read(10,*,end=900) irev,time,sfi,sla,ssh,sigma
      sh = ssh
      istat = irev
      sg(1) = ssh
      sg(2) = sigma
      sg(3) = time
      goto 240
c
c  iang = 5: binary format - australian gravity data
c
233   read(10,end=900) qfi,qla,qg,qh
      istat = istat + 1
      sfi = qfi
      sla = qla
      sh = qh
      gg = qg + 978000.d0
      gg = (gg-979685.74d0)*1.0005118d0 + 979671.86d0
      s = sin(sfi/180*3.14159265)**2
      gamma = 978032.7d0*(1 + .0052790414*s + .0000232718*s**2)
      sg(1) = gg - (gamma - 0.3086*sh)
      sg(2) = sg(1) - 0.1119*sh
      if (sh.lt.-100.or.sh.gt.2500) then
        nn1 = nn1 + 1
        goto 233
      elseif (sg(1).lt.-400.or.sg(1).gt.400) then
        nn2 = nn2 + 1
        goto 233
      endif
      go to 240
c
c  iang = 6: data on 80-character standard gravity data format
c
 235  read(10,236,end=900) sfi,sla,hc,sh,gg,fa,ipub,istat,tc,
     .tccode,err,ba
 236  format(1x,f8.2,1x,f9.2,1x,a1,f7.2,8x,f7.2,1x,f7.2,1x,
     .i3,i7,f5.1,a1,f5.1,f7.2)
      istat = ipub*1000000 + mod(istat,1000000)
      if (hc.eq.'3'.or.hc.eq.'4'.or.hc.eq.'5') sh = 0
      if (tccode.eq.'X') tc = 9999.99
      if (nstat.eq.1)
     .write(*,*) '- sea heights set to zero -'
      sfi = cdeg(sfi) 
      sla = cdeg(sla)
      if (n.eq.1) then
        sg(1) = fa
      elseif (n.eq.2) then
        sg(1) = fa
	if (lminus) then
	  sg(2) = ba
        else
          sg(2) = err
        endif
      elseif (lminus.and.n.eq.3) then
        sg(1) = fa
        sg(2) = err
        sg(3) = tc
      else
        sg(1) = gg+976000
        sg(2) = fa
        sg(3) = ba
	sg(4) = tc
      endif
      goto 240
c
c  iang = 7: data in grid format, values >= 9999 signals unknown
c
237   ieg = ieg+1
      istat = istat+1
      if (ieg.gt.neg) then
        ing = ing-1
        ieg = 1
        if (ing.eq.0) goto 900
        read(10,*) (grow(k),k=1,neg)
        sfi = glab(1) + (ing-1)*dfig
        if (lgutm) ufig = sfi
      endif
      if (grow(ieg).ge.9999) goto 237
      sla = glab(3)+(ieg-1)*dlag
      sg(1) = grow(ieg)
      if (n.eq.0) sh = sg(1)
      if (lgutm) then
        ulag = sla
        call utg(ufig,ulag,sfi,sla,sag,.true.,.false.)
        sfi = sfi*radeg
        sla = sla*radeg
      endif
c
c  select and output data 
c  -----------------------
c
240   nstat = nstat+1
      if (abs(sfi).gt.90.and.abs(sla).gt.360) then
        lutmo = .true.
      else
        if (sla.gt.180.0) sla=sla-360.0
      endif
      if (slatma.lt.sfi) slatma=sfi
      if (slatmi.gt.sfi) slatmi=sfi
      if (slonma.lt.sla) slonma=sla
      if (slonmi.gt.sla) slonmi=sla
      if (lref) goto 241
      if (lwind) then
        if (rfi1i.le.sfi.and.sfi.le.rfi2i.and.
     *  rla1i.le.sla.and.sla.le.rla2i) goto 200
      endif
c
c  convert point to utm if needed
c
      if (lutm) then
        if (lutmo) then
          ufi = sfi
          ula = sla
        else
          call utg(sfi/radeg,sla/radeg,ufi,ula,sa,.false.,.false.)
        endif
      else
        ufi = sfi
        ula = sla
      endif
c
c  output if all points wanted within area
c  ---------------------------------------
c
      if (.not.all) goto 250
c
      if (ufi.lt.rfi1.or.ufi.gt.rfi2.or.
     *ula.lt.rla1.or.ula.gt.rla2) goto 200
      if (lrej.and.sg(2).gt.rejlev) then
        nrej = nrej+1
        goto 200
      endif
241   if (iang.ne.4) then
        if (lutmo) then
          write(20,912) istat,sfi,sla,sh,(sg(k),k=1,n)
        else
          write(20,910) istat,sfi,sla,sh,(sg(k),k=1,n)
        endif
      else
        if (lutmo) then
          write(20,913) istat,sfi,sla,sh,(sg(k),k=1,n)
        else
          write(20,911) istat,sfi,sla,sh,(sg(k),k=1,n)
        endif
      endif
      do 242 j = 1, n
        if (sg(j).lt.ssmin(j)) ssmin(j) = sg(j)
        if (sg(j).gt.ssmax(j)) ssmax(j) = sg(j)
        ss1(j) = ss1(j) + sg(j)
        ss2(j) = ss2(j) + sg(j)**2
242   continue
      nsel = nsel+1
      goto 200
c
c  select and output data - pixel search or averaging
c  --------------------------------------------------
c
250   ii = (ufi-rfi0)/dfi+1.00001
      jj = (ula-rla0)/dla+1.00001
      if (ii.lt.1.or.ii.gt.nn.or.jj.lt.1.or.jj.gt.ne) goto 200
      if (lrej.and.sg(2).gt.rejlev) then
        nrej = nrej+1
        goto 200
      endif
      rf = (ii-1)*dfi+rfi1
      rl = (jj-1)*dla+rla1
      kk = (ii-1)*ne+jj
      if (lmean) goto 300
c
c  pixel search  - use composite stddev+distance with rejlev negative
c
      if (lmin) then
        rr = sg(1)
      else
        rr = (ufi-rf)**2+((ula-rl)*cosfi)**2
        if (lrej1) rr = rr/rr2 + (1 + sg(2))
      endif
      if (rr.ge.d(kk)) goto 200
      d(kk) = rr
      is(kk) = istat
      rfi(kk) = sfi
      rla(kk) = sla
      rh(kk) = sh
      do 275 k=1,n
275   g(kk,k) = sg(k)
      goto 200
c
c  averaging
c
300   is(kk) = is(kk) + 1
      rfi(kk) = rfi(kk) + sfi
      rla(kk) = rla(kk) + sla
      rh(kk) = rh(kk) + sh
      d(kk) = d(kk) + sh**2
      if (n.gt.0) then
        rr = sg(n)
      else
        rr = sg(1)
      endif
      g(kk,1) = g(kk,1) + rr
      g(kk,2) = g(kk,2) + rr**2
      goto 200
c
c  output primitive maps if required
c
900   if (lref) then
        goto 1000
      endif
      if (.not.lmap.or.all) goto 905
      do 901 k = 1, n
901   call pmap(20,g,k,is,nn,ne,rfi1,rfi2,rla1,rla2,dfi,dla,iamax)
      write(*,902)
902   format(' - only pixel coverage map written on outputfile -')
c
c  output selected data in pixel modes
c  -----------------------------------
c
905   if (all) goto 1000
      if (lmean) goto 950
      idum=0
      do 920 i = nn,1,-1
      do 920 j = 1,ne
        kk = (i-1)*ne+j
        if (is(kk).le.0) goto 920
        nsel = nsel+1
        do 906 k = 1,n
          rg = g(kk,k)
          if (rg.lt.ssmin(k)) ssmin(k) = rg
          if (rg.gt.ssmax(k)) ssmax(k) = rg
          ss1(k) = ss1(k) + rg
          ss2(k) = ss2(k) + rg**2
906     continue
c
c  add noise to selected data if wanted
c
        if (lnoise) then
          topi = 6.283185308
          do 916 k=1,n
            r1 = randl(idum)
            r2 = randl(idum)
            sq = sqrt(-2.0*log(r2))
            g1(k)  = cos(topi*r1)*sq*stdvn(k)
            if (abs(g1(k)).gt.(4.0*stdvn(k))) g1(k) = 0.0
            s1(k)=s1(k)+g1(k)
            s2(k)=s2(k)+g1(k)**2
            g(kk,k)=g(kk,k)+g1(k)
916       continue
        endif
c
        if (lmap) goto 920
        if (iang.ne.4) then
          if (lutmo) then
            write(20,912) is(kk),rfi(kk),rla(kk),rh(kk),
     *      (g(kk,k),k=1,n)
          else
            write(20,910) is(kk),rfi(kk),rla(kk),rh(kk),
     *      (g(kk,k),k=1,n)
          endif
        else
          if (lutmo) then
            write(20,911) is(kk),rfi(kk),rla(kk),rh(kk),
     *      (g(kk,k),k=1,n)
          else
            write(20,911) is(kk),rfi(kk),rla(kk),rh(kk),
     *      (g(kk,k),k=1,n)
          endif
        endif
910     format(' ',i9,f10.5,f11.5,f10.3,' ',f10.3,8f9.3)
911     format(' ',i9,f10.5,f11.5,3f9.2,6f13.1)
912     format(' ',i9,f10.1,f11.1,f9.3,' ',f10.3,8f9.3)
913     format(' ',i9,f10.1,f11.1,3f9.2,6f13.1)
920   continue
      goto 1000
c
c  output mean data results
c
950   if (mode.eq.2) write(*,952)
952   format(' - mean height and data generation -'/,
     .' - output data: no, lat, lon, mean h, mean data, np, sig(h), ',
     .'sig(data)')
      do 960 i = nn,1,-1
      do 960 j = 1,ne
        kk = (i-1)*ne+j
        ii = is(kk)
        if (ii.eq.0) goto 960
        nsel = nsel + 1
        rg = g(kk,1)/ii
        if (rg.lt.ssmin(1)) ssmin(1) = rg
        if (rg.gt.ssmax(1)) ssmax(1) = rg
        ss1(1) = ss1(1) + rg
        ss2(1) = ss2(1) + rg**2
        if (mode.eq.3) then
	  g(kk,3) = rg
	  goto 960
        endif
        rr = 0.0
c  sigma may be negative due to rounding errors in real*4
        if (ii.gt.1) rr = (g(kk,2)-g(kk,1)**2/ii)/(ii-1)
        if (rr.lt.0) rr = -sqrt(-rr)
        if (rr.gt.0) rr = sqrt(rr)
        sh = 0.0
        if (ii.gt.1) sh = (d(kk)-rh(kk)**2/ii)/(ii-1)
        if (sh.gt.0) sh = sqrt(sh)
        if (sh.lt.0) sh = -sqrt(-sh)
        if (.not.lutm) 
     *  write(20,956) (nn-i)*ne+j,rfi(kk)/ii,rla(kk)/ii,
     *  rh(kk)/ii,rg,ii,sh,rr
        if (lutm) 
     *  write(20,957) (nn-i)*ne+j,rfi(kk)/ii,rla(kk)/ii,
     *  rh(kk)/ii,rg,ii,sh,rr
956     format(' ',i6,2f12.5,f9.2,f9.2,i7,2f9.2)
957     format(' ',i6,2f12.0,f9.2,f9.2,i7,2f9.2)
960   continue
      if (mode.ne.3) goto 1000
c
c  output mean data grid
c
      write(*,962)
962   format(' - mean data grid written on outputfile -')
      if (.not.lutm) write(20,970) rfi1,rfi2,rla1,rla2,dfi,dla
970   format(' ',4f12.6,2f12.8)
      if (lutm) write(20,971) rfi1,rfi2,rla1,rla2,dfi,dla,iell,izone
971   format(' ',6f12.0,/,' ',2i4)
      do 972 i = nn, 1, -1
        kk = (i-1)*ne
        if (n.ne.0) write(20,973) (g(kk+j,3),j=1,ne)
        if (n.eq.0) write(20,9731) (nint(g(kk+j,3)),j=1,ne)
973     format(20(/,' ',9f8.2))
9731    format(20(/,' ',12i6))
972   continue
c
c  output statistics etc on selected data
c  --------------------------------------
c
1000  write(*,1003) nstat,slatmi,slatma,slonmi,slonma
1003  format(' total points:',i7, 
     */' located within area:',4f13.4)
      all = all.or.lref
      if (.not.all.and.(.not.lutm)) write(*,1004) rfi1,rfi2,rla1,rla2,
     *dfi,dla
1004  format(' wanted pixel grid:  ',4f10.4,2f8.4)
      if (.not.all.and.lutm) write(*,1005) rfi1,rfi2,rla1,rla2,
     *dfi,dla,iell,izone
1005  format(' wanted pixel grid:  ',6f9.0,2i3)
      write(*,1006) nsel, nn*ne
1006  format(' no of output/selected points:',i7,
     *', total poss. pixels:',i7)
      if (lrej) write(*,1010) nrej
1010  format(' number of points in area rejected by too large stddev:',
     *i8)
      if (nsel.gt.1) write(*,1007)
1007  format(' selected data:  mean    std.dev.     min       max')
      if (nsel.gt.1) then
        do 1009 j = 1, n
        write(*,1008) j,ss1(j)/nsel,
     *  sqrt((ss2(j)-ss1(j)**2/nsel)/(nsel-1)),ssmin(j),ssmax(j)
1008    format('     no: ',i1,'  ',4f10.2)
1009    continue
      endif
      if (iang.eq.5) write(*,1020) nn1,nn2
1020  format(' special binary option, nn1, nn2 = ',2i9)
      if (lnoise) then 
        do 1022 k=1,n
        if (nsel.gt.1) s2(k)=sqrt((s2(k)-s1(k)**2/nsel)/(nsel-1))
        write(*,1023) k,s1(k)/nsel,s2(k)
1022     continue
1023    format(' n mean and stdv. of noise',i3,2f10.4)
      endif
      close(20)
      end
c
      subroutine pmap(iunit,g,k,is,nn,ne,rfi1,rfi2,rla1,rla2,dfi,dla,
     *iamax)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                          p m a p
c
c  this subroutine will draw a very primitive map of data given in
c  array g, using different symbols to illustrate the values in a
c  number of levels. values with is.le.0 are treated as unknown.
c
c  rf, dec 88, unsw
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      real*4 g(iamax,3)
      dimension is(iamax)
      character*1 cl(120)
      character*1 clev(10)
      data clev /'c','b','a','=','-','+','*','a','b','c'/
c
      if (nn.gt.120.or.ne.gt.120) stop 'too many pixels for map'
      nlev = 10
c
      rmax = -9999.9
      rmin = 9999.9
      n = nn*ne
      np = 0
      do 5 i = 1, n
        if (is(i).le.0) goto 5
        np = np + 1
        if (g(i,k).lt.rmin) rmin = g(i,k)
        if (g(i,k).gt.rmax) rmax = g(i,k)
5     continue
      if (rmin.eq.9999.9) stop 'no data selected within area'
      dd = (rmax-rmin)/nlev
c
      write(iunit,7) rfi1,rfi2,rla1,rla2,dfi,dla,k
7     format(/' coverage map, area: ',4f8.2,/' grid spacing: ',2f8.4,
     ./' data number',i2,', levels:')
      do 977 j=1,nlev
      rmj=rmin+j*dd
      rmj1=rmin+(j+1)*dd
977   write(iunit,77)clev(j),rmj,rmj1
 77   format(10(/' ',2(a1,' ',f8.2,' to',f8.2,8x)))
      write(iunit,10) (mod(j,10),j=1,ne)
10    format('  ',120i1)
c
      do 20 i = nn, 1, -1
        do 21 j = 1, ne
           kk = (i-1)*ne+j
           if (is(kk).le.0) then
             cl(j) = ' '
             goto 23
           endif
           ip = (g(kk,k)-rmin)/dd + 1.0
           if (ip.lt.1) ip = 1
           if (ip.gt.nlev) ip = nlev
           cl(j) = clev(ip)
23         continue
21       continue
         write(iunit,22) mod(i,10), (cl(j),j=1,ne)
22       format(' ',i1,120a1)
20     continue
       write(iunit,10) (mod(j,10),j=1,ne)
       return
       end
c
      double precision function  randl(j)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                         r a n d l
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c obtained from t.risbo aug. 89.
      integer*4 j
      real   *8 r
c
c     pseudo-random numbers - linear distribution.
c
      real   *8 x0,xs,xm,xp
      data xs /1.d+9/, xm/16807.d+0/,xp /2.147483647 d+9/
c
      if (j .ne. 0) then
        x0 = j
      else
        x0 = xs
      endif
c
      x0 = dmod(xm*x0, xp)
      j = x0
      r = x0/xp
c
      randl = r
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
c special option: polar stereographic. Implemented as negative iell.
c        abs iell = latitude of true scale  
c        izone =    longitude of y-axis
c polar stereographic with scale 1 at pole and oriented to
c Greenwich is thus (iell,izone) = (-90,0)
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision n,m
c
      radeg = 180/3.1415926536d0
c
c  set ellipsoid parameters
c  check for swedish special projection
c  check for polar stereographic, use sa(22) as flag
c
      if (isys.lt.0) then
        write(*,*) '- polarstereographic, true scale at lat: ',
     .  abs(isys),', ref lon: ',izone,' -'
        rlat_true = abs(isys)
        rlon_down = izone
        call set_pst(rlat_true,rlon_down,0.d0,0.d0)
        sa(20) = 1.d0
        sa(21) = 1/111000.d0
        sa(22) = 1.d0 
        return
      endif
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
      sa(22) = 0
      return
      end
c
      subroutine utg(rn, re, b, l, sa, direct, tcheck)
      implicit double precision(a-h,o-z)
      double precision b, l, sa(22)
      logical direct, tcheck
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c stereograpic rf aug 2007
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision np, lg, l0, ndif, ncheck
      integer h, s
      logical lpolar
      double precision n, e
c
      n = rn
      e = re
      qn = sa(1)
      e0 = sa(2)
      l0 = sa(3)
      lpolar = (sa(22).eq.1)
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
10    if (lpolar) then
        call pst_to_geo(rn,re,bbg,lg)
        goto 100
      endif
      np = n/qn
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
30    if (lpolar) then
        call geo_to_pst(n,e,bbg,lg)
        goto 100
      endif
      bbg   = n + clsin(sa, 3, 4, 2.0*n)
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
100   continue
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
      write(*, 90) n, e, ndif*ep, edcos*ep
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
