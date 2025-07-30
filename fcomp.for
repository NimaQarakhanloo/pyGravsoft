      program fcomp
      implicit double precision(a-h,o-z)
c $Id: fcomp.for 267 2009-01-27 14:36:10Z cct $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                             f c o m p
c
c  program reads data from two files and produce sum
c  or difference and writes statistics. the program may
c  also add bias to data, or merge two data files.
c
c  input: inputfile1,
c         inputfile2 (evt. dummy name),
c         outputfile (or '0'),
c         kind, ndata, hist.spacing.
c
c  where kind:  0: just statistics for ifile1
c               1: difference ifile1 minus ifile2
c               2: sum
c               3: factor and bias to ifile1 (one for each data)
c               4: merge file1 and file2 ('ndata' points in file1, 1 in file2)
c               5: subtract file1-file2 only when 2nd field in file1
c                  less than 'trsh' (input), otherwise reject
c               6: as 5, except 2nd field is in file2
c        10,11,..: as 0,1,.. but with long data lines in file 1,
c                  additional input: total number of data and wanted
c                  data positions (file 2 standard format)
c                  wanted data pos = 99 implies BA on land, FA at sea at first
c                  and second field, respectively
c              -1: differences between successive points in ifile1,
c                  profile statistics for both raw differences and ppm,
c                  useful for gps traverses.
c
c  ndata is the number of data following id,lat,lon,elev
c  hist.spacing is histogram spacing, if 0 no histogram is output.
c  if outputfile = '0' no output is written (saves disc space!)
c
c  programmer: rene forsberg, jan 85
c  last updated: jul 93, rf, may 95, rf
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension r(3,9), rold(9)
      logical lcol,ldif,lfile2,lofile,lbias,lmerge,lutm,lgeo,ltrsh
      dimension data(18),bias(9),fact(9)
      dimension idata(9)
      dimension sum(3,9),sum2(3,9),rmin(3,9),rmax(3,9)
      dimension ihist(21,3,9)
      character ifile1*128,ifile2*128,ofile*128
c
      write(*,1)
1     format(' *******************************************************',
     ./,     ' *   FCOMP  -  GRAVSOFT file comparison  - (c) RF/KMS  *',
     ./,     ' *******************************************************',
     ./' input file names (ifile1/ifile2/ofile): ')
      read(*,2) ifile1
      read(*,2) ifile2
      read(*,2) ofile
2     format(a128)
      write(*,12)
12    format(' input kind (0:stat, 1:dif, 2:sum, 3:fak/bias, 4:mrg),',
     .' ndata, hist.sp.: ')
      read(*,*) kind, ndata, rhist
      lcol = .false.
      if (kind.ge.10) then
        if (ldif) stop 'not implemented'
        kind = kind-10
        lcol = .true.
        write(*,901)
901     format(' input no of data in line and wanted position(s): ')
        read(*,*) ndl,(idata(j),j=1,ndata)
      endif
      ltrsh = (kind.eq.5.or.kind.eq.6)
      if (ltrsh) then
        write(*,*) 'input rejection treshold: '
        read(*,*) trsh
      endif
      nn = 3
c
      ldif = (kind.eq.-1)
      lbias = (kind.eq.3)
      lmerge = (kind.eq.4)
      nrej = 0 
      if (lmerge.or.lbias.or.ldif) kind = 0
      lfile2 = (kind.gt.0.or.lmerge)
      if (lbias) lfile2 = .false.
      lofile = (kind.gt.0.or.lcol.or.ldif.or.lbias.or.lmerge)
      if (ofile.eq.'0') lofile = .false.
c
      open(20, file=ifile1, status='old')
      if (lfile2) open(21, file=ifile2, status='unknown')
      if (lofile) open(30, file=ofile, status='unknown')
c
      n = 0
      do 10 i=1,3
      do 10 j=1,9
        sum(i,j) = 0.0
        sum2(i,j) = 0.0
        rmin(i,j) = 999999.99
        rmax(i,j) = -999999.99
        do 10 k=1,21
          ihist(k,i,j) = 0
10    continue
c
      if (ndata.gt.9.or.lmerge.and.ndata.gt.4)
     .stop '*** ndata too big, sorry ***' 
      ndatao = ndata
      if (lmerge.or.ltrsh) ndatao = ndata+1
      if (kind.eq.0.and.(.not.(ldif.or.lbias))) nn = 1
      if (lmerge) nn = 2
c
      if (lbias) then
        write(*,902)
902     format(' input factor and bias(es) to modify data: ')
        read(*,*) (fact(j),bias(j),j=1,ndata)
      endif
c
      ndata1 = ndata
      ndata2 = ndata
      if (kind.eq.5) ndata1 = ndata+1
      if (kind.eq.6) ndata2 = ndata+1
      if (lmerge) ndata2 = 1
c
c  loop for data points
c  --------------------
c
      read(20,*,end=90) istold,rfiold,rlaold,rhold,
     .(rold(i),i=1,ndata1)
      if (.not.ldif) rewind(20)
      lutm = (abs(rfiold).gt.400.or.abs(rlaold).gt.400)
      lgeo = (.not.lutm)
      if (lutm) write(*,19)
19    format(' - file1 and outputfile assumed to be utm -')
c
c  loop entry
c
20    if (lfile2) then
        read(21,*,end=90) ist2,rfi2,rla2,rh2,(r(2,i),i=1,ndata2)
        if (kind.eq.6) tpar = r(2,ndata2)
      endif
c
203   if (lcol) then
        read(20,*,end=90) ist1,rfi1,rla1,rh1,(data(j),j=1,ndl)
        do 201 j=1,ndata
        kk = idata(j)
	if (idata(j).eq.99) then
	  kk = 2
	  if (rh1.lt.0) kk = 1
        endif
201     r(1,j) = data(kk)
      else
        read(20,*,end=90) ist1,rfi1,rla1,rh1,(r(1,i),i=1,ndata1)
        if (kind.eq.5) tpar = r(1,ndata1)
      endif
      if (ltrsh) then
	if (tpar.gt.trsh) then
          nrej = nrej +1 
          goto 20
        endif
      endif
c    
      if (ist1.ne.ist2.and.lfile2) then
        if (abs(rfi1-rfi2).gt.0.1D-8.or.abs(rla1-rla2).gt.0.1D-8) then
         write(*,*) 
     .   '*** station',ist1,' not in file2 - skipped in file 1'
         goto 203
        end if
      endif
c
c  baseline differences
c
      if (ldif) then
        if (lutm) then
          dr = sqrt((rfi1-rfiold)**2 + (rla1-rlaold)**2)*1000
        else 
          dr = sqrt((rfi1-rfiold)**2 +
     .    ((rla1-rlaold)*cos(rfi1/57.29578))**2) * 111.195
        endif
        if (dr.le.0) stop 'identical points on profile'
        do 212 i = 1, ndata
          r(2,i) = r(1,i)-rold(i)
          r(3,i) = abs(r(2,i))/dr*1000
          rold(i) = r(1,i)
          r(1,i) = dr
212     continue
        if (lofile) write(30,213)
     .  istold,ist1,dr,(r(2,i),i=1,ndata),(r(3,i),i=1,ndata)
213     format(' ',2i6,f9.1,3(6f9.2,/))
        istold = ist1
        rfiold = rfi1
        rlaold = rla1
        goto 25
      endif
c
c  difference/sum
c
      do 21 i=1, ndata
        if (kind.eq.1.or.ltrsh) then
          r(3,i) = r(1,i)-r(2,i)
        elseif (kind.eq.2) then
          r(3,i) = r(1,i)+r(2,i)
        else
          r(3,i) = r(1,i)
          if (lbias) then
            r(2,i) = 0.d0
            r(3,i) = fact(i)*r(1,i) + bias(i)
          elseif (lmerge) then
            if (i.eq.1) r(3,ndata+1) = r(2,1)
          endif
        endif
21    continue
      if (ltrsh) r(3,ndatao) = tpar
c
c  write output data
c
      if (lofile) then
        if (lgeo) then
          write(30,22) ist1,rfi1,rla1,rh1,(r(3,i),i=1,ndatao)
22        format(' ',i9,2f11.5,f9.2,9f10.3)
        else
          write(30,221) ist1,rfi1,rla1,rh1,(r(3,i),i=1,ndatao)
221       format(' ',i9,2f11.0,f9.2,9f10.3)
        endif
        if (ist1.ne.ist2.and.lfile2) write(*,24) ist1,ist2
 24     format(' *** warning, different station numbers: ',2i7)
      endif
c
25    n = n+1
      do 30 i= 1,nn
      do 30 j= 1,ndata
        sum(i,j) = sum(i,j) + r(i,j)
        sum2(i,j) = sum2(i,j) + r(i,j)**2
        if (r(i,j).lt.rmin(i,j)) rmin(i,j) = r(i,j)
        if (r(i,j).gt.rmax(i,j)) rmax(i,j) = r(i,j)
        if (rhist.gt.0) ii = r(i,j)/rhist+11.5
        if (ii.lt.1) ii = 1
        if (ii.gt.21) ii = 21
        ihist(ii,i,j) = ihist(ii,i,j)+1
30    continue
c
      goto 20
c
c  exit of reading loop
c  --------------------
c
90    if (kind.eq.0) write(*,91) n
91    format(/' --- fcomp statistics --- number of points: ',i5)
      if (ldif) write(*,911)
911   format(/' successive differences: '/,
     .' file1 = distances in km, file2 = differences, file3 = abs ppm ')
      if (kind.eq.1.or.kind.eq.5) write(*,92) n
92    format(' fcomp statistics, file3 = file1 - file2,',
     *' number of points: ',i5)
      if (kind.eq.2) write(*,93) n
93    format(' fcomp statistics, file3 = file1 + file2,',
     *' number of points: ',i5)
      if (kind.eq.5) write(*,*) 'number of rejected obs: ',nrej
      if (rhist.gt.0) write(*,94) rhist
94    format(' histogram spacing: ',f7.2)
      if (n.eq.0) stop 'no data in file'
c
c
      if (ldif) nn=3
      do 120 i = 1,nn
        write(*,100) i
100     format(' file ',i1,
     *  ': dno   mean  stddev     rms     min     max')
        do 110 j=1,ndata
          rm = sum(i,j)/n
          if (n.eq.1) rstd = 0.0
          if (n.gt.1) rstd = (sum2(i,j)-sum(i,j)**2/n)/(n-1)
        rstd = sqrt(abs(rstd))
        rms = sqrt(sum2(i,j)/n)
        write(*,104) j,rm,rstd,rms,rmin(i,j),rmax(i,j)
104     format('        ',i3,5f8.2)
        if (rhist.gt.0) write(*,105) (ihist(k,i,j),k=1,21)
105   format('               ',21i3/,'                ',
     *  ' * -9 -8 -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7  8',
     *  '  9  *')
110   continue
120   continue
      end
