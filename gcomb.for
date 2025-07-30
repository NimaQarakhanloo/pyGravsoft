      program gcomb
      implicit double precision(a-h,o-z)
c $Id: gcomb.for 181 2008-08-14 23:19:51Z tjansson $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                          g c o m b
c
c  program for combining two grids into one by various manipulations.
c  the grids must have the same grid spacings and relative positions,
c  but not neccessarily the same coverage area.
c  the output grid will cover the same area as the first grid
c  9999-values (unknown) are not added/subtracted.
c  input/output may be either binary or text format
c
c  input:
c
c  gridfile1               or           gridfile
c  gridfile2                            * 
c  outfile                              (lat1,lat2,lon1,lon2 - interactive)
c  mode, iout                           
c
c  where mode determines function of program (except *-mode) 
c
c  mode = 0:  one-grid modification. Additional input: factor, bias
c             (new  = old*factor + bias)
c
c  mode = 1:  subtract 'grid1' minus 'grid2'
c
c  mode = 2:  add 'grid1' plus 'grid2'
c
c  mode = 3:  general grid convert: output = a*grid1 + b*grid2
c             additional input: a, b
c
c  mode = 4:  treshold values. set values in 'grid1' to 9999 when
c             lgt=true:  value in 'grid2' > 'trsh' (input lgt,trsh after iout).
c             lgt=false: value in 'grid2' < 'trsh'  
c
c  mode = 5:  grid overwrite
c             values in 'grid2' is written on top of 'grid1', except when
c             9999 is encountered, then the grid1-values are kept.
c             (special: mode = -5 allows many text grids in 'grid2' file,
c             output must be binary)
c
c  mode = 6:  grid select
c             values in 'grid1' are written only when there is no data in
c             'grid2'.
c
c  mode = 7:  Treshold select. Keep values in grid1 when grid1 > treshold,
c             keep values in grid2 when both grid1 and grid2 < trsh
c
c  mode = 8:  Variable multiplication. Multiply grid 1 by "fak1" when
c             grid 2 .lt. 'treshold', else multiply by "fak2"
c
c  mode = 9:  Grid 1 is a Bouguer anomaly grid, grid 2 a height grid.
c             Output N - zeta (difference geoid minus quasigeoid) in
c             linear approximation. If grid 1 is a free-air grid,
c             the difference zeta* - zeta is obtained.
c
c  mode =10:  Grid 1 is a geoid model with Helmert, grid2 a height grid.
c             Output is a geoid corrected for indirect effect.
c
c  mode =11:  Grid 1 is a quasigeoid at sea level, grid 2 a height grid.
c             Output is a classical geoid N
c
c  iout:      output format:
c             1  binary
c             2  txt, reals
c             3  txt, integers
c
c  special options:
c
c  if 'grid1' = 0 or 9999 then a label must be input, and grid1 will be
c  all zeros or 9999's. 
c
c  if 'outfile' = 0 the output is on the screen without statistics.
c  if both 'grid1' = '0' and 'outfile' = '0' this gives interactive listing
c
c  if 'grid1' and 'outfile' are identical only the relevant rows are
c  updated for binary files (gives fast overwrite)
c
c  Input example - use overwr to extract subgrid from large binary file:
c      0                       
c      binfile
c      subgrid
c      5 3  
c      55.0 55.1  10.0 10.1  .01  .01        (label)
c
c  rene forsberg, thessaloniki/unsw nov 1988, last updated sep 1989
c  updated and changed dec 89, feb 92, rf
c  updated jan 94, chung cheng inst of technology, taiwan, rf
c  update with binary grids, march 1996, rf
c  update with indirect helmert geoid effect, may 21, 2000
c  update with treshold less than, jan 23, 2002
c  geoid update mode11, sep 2003
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension glab1(6),glab2(6)
      real*4 grow1(4000),grow2(4000) 
      character*36 ifile1,ifile2,ofile
      logical lzero, l9999, llab, lint, lutm, libin1, libin2, lobin,
     .lupd, lobinw, ldos, lmulti, lfirst, lgt
c
      ldos = .false.
      write(*,3)
3     format(/
     .' ************************************************************',
     .'********'/
     .' *   GCOMB - GRAVSOFT grid combination - vers. Sep 2003 (c) R',
     .'F/KMS  *'/
     .' ************************************************************',
     .'********')
      write(*,1)
1     format(' input file names - file1/file2/ofile  or  file/*: ')
      if (ldos) write(*,*)
      read(*,2) ifile1
      read(*,2) ifile2
2     format(a36)
c
      if (ifile2.eq.'*') then
        mode = 5
        iout = 3
        ifile2 = ifile1
        ifile1 = '9999'
        ofile = '0'
      else
        read(*,2) ofile
        write(*,4)
4       format(
     .  ' input: mode (0 sgl, 1 dif, 2 sum, 3 faksum, 4 trsh, ',
     .  '5 overwr..)'/,
     .  '        iout (1 bin, 2 txt, 3 txt int)'
     .  )
        if (ldos) write(*,*)
        read(*,*) mode, iout
      endif
c
      if (mode.gt.11) stop 'mode parameter undefined'
      if (mode.eq.4) write(*,5)
5     format(' input lgt,grid2 treshold val for setting grid1 to 9999:')
      if (mode.eq.7) write(*,502)
502   format(' input treshold select value: ')
      if (mode.eq.4) read(*,*) lgt,trsh
      if (mode.eq.7) read(*,*) trsh
      if (mode.eq.8) write(*,*) 'input treshold, fak1, fak2:'
      if (mode.eq.8) read(*,*) trsh,fak1,fak2
      if (mode.eq.0) write(*,501)
501   format(' input factor and bias to be added: ')
      if (mode.eq.0) read(*,*) fact,bias
      if (mode.eq.3) write(*,*) 'input: fak1 and fak2'
      if (mode.eq.3) read(*,*) fak1,fak2
c
      lmulti = (mode.lt.0)
      lfirst = .true.
      if (lmulti) mode = abs(mode)
      lobin = (iout.eq.1)
      lint = (iout.eq.3)
      if (ofile.eq.'0') then
        i20 = 0
        if (lobin) stop 'output must be txt format'
      else
        i20 = 20
      endif
      lzero = (ifile1.eq.'0')
      l9999 = (ifile1.eq.'9999')
      llab = (lzero.or.l9999)
c
c  interactive reentry
c
55    if (llab) then
        write(*,*) 'input wanted area (lat1,lat2,lon1,lon2): '
        read(*,*) (glab1(j),j=1,4)
        if (mode.eq.0) then
          write(*,*) 'input grid spacing: '
          read(*,*) glab1(5),glab1(6)
        endif
      else
        call openg(10,ifile1,libin1)
      endif
c
      lupd = (ofile.eq.ifile1.and.mode.eq.5.and.libin1)
      if (lupd.and.(.not.lobin)) stop 'iout must be 1'
      lobinw = (lobin.and.(.not.lupd))
c
      if (lobin) then
        if (lobinw) 
     .  open(20,file=ofile,form='unformatted',status='unknown',
     .  access='direct',recl=64)
      else
        if (lmulti) stop 'output file must be binary for multi grids'
        if (i20.eq.20) open(20,file=ofile,status='unknown')
      endif
c
      if (.not.llab) then
        if (libin1) then 
          read(10,rec=1) icode,glab1,lutm,iell,izone
          if (icode.ne.777) stop 'binary check code 777 missing'
        else
          read(10,*) glab1
          lutm = (abs(glab1(1)).gt.100.or.abs(glab1(2)).gt.400)
          if (lutm) read(10,*) iell,izone
        endif
      endif
c
c  open grid 2 - entry point for repeated grids
c
8     if (mode.gt.0) then
        if (lfirst) call openg(11,ifile2,libin2)
        if (libin2) then
          read(11,rec=1) icode,glab2,lutm,iell2,izone2
          if (icode.ne.777) stop 'binary code 777 missing on grid2'
          if (lmulti) stop 'multiple grids must be txt format'
        else
          read(11,*,end=90) glab2
          lutm = (abs(glab2(1)).gt.100.or.abs(glab2(2)).gt.400)
          if (lutm) read(11,*) iell2,izone2
        endif
        if (llab) then
          if (lutm) iell = iell2
          if (lutm) izone = izone2
          glab1(5) = glab2(5)
          glab1(6) = glab2(6)
        endif
c
        if (lutm.and.(iell.ne.iell2.or.izone.ne.izone2))
     .  stop 'utm spec difference'
        if (abs(glab1(5)-glab2(5)).gt.0.00001.or.
     .  abs(glab1(6)-glab2(6)).gt.0.00001)
     .  stop 'grid spacings different'
      else
        do 9 j = 1, 6
9       glab2(j) = glab1(j)
      endif
c
      if (lobin) then
        if (lobinw) write(20,rec=1) 
     .  777,(glab1(j),j=1,6),lutm,iell,izone
      else
        if (.not.lutm) write(i20,11) (glab1(j),j=1,6)
11      format(' ',6f12.6)
        if (lutm) write(i20,12) (glab1(j),j=1,6),iell,izone
12      format(' ',6f12.1,/,' ',2i4)
      endif
c
      rfi0 = glab1(1)
      rla0 = glab1(3)
      dfi = glab1(5)
      dla = glab1(6)
      nn = nint((glab1(2)-rfi0)/dfi + 1)
      ne = nint((glab1(4)-rla0)/dla + 1)
      nn2 = nint((glab2(2)-glab2(1))/dfi + 1)
      ne2 = nint((glab2(4)-glab2(3))/dla + 1)
      in1 = nint((glab2(1)-rfi0)/dfi + 1)
      in2 = nint((glab2(2)-rfi0)/dfi + 1)
      ie1 = nint((glab2(3)-rla0)/dla + 1)
      ie2 = nint((glab2(4)-rla0)/dla + 1)
c
      r1min = 9999.99 
      r1max = -r1min
      nr1 = 0
      r1sum = 0.0
      r2sum2 = 0.0
      r2min = 9999.99 
      r2max = -r2min
      nr2 = 0
      r2sum = 0.0
      r2sum2 = 0.0
      rrmin = 9999.99 
      rrmax = -rrmin
      nrr = 0
      rrsum = 0.0
      rrsum2 = 0.0
c
      if (i20.eq.20) write(*,20) 1,nn,1,ne
 20   format(/' ---  G C O M B  ---',/,
     .' number of points in grid1,  north:',2i7,', east:',2i7)
      if (i20.eq.20.and.mode.ne.0) write(*,201) in1,in2,ie1,ie2
201   format(' relative position of grid2, north:',2i7,', east:',2i7) 
      if (libin1.and.(.not.llab)) write(*,*) 'grid1 in binary format'
      if (libin2.and.mode.ne.0) write(*,*) 'grid2 in binary format'
      if (lobinw) write(*,*) '- output grid in binary format'
      if (lupd) write(*,*) '- direct update of binary grid1'
c
c  skip records in grid2 north of grid1
c
      if (in2.gt.nn.and.(.not.libin2).and.mode.gt.0) then 
        do 21 i = 1, in2-nn
21      read(11,*) (grow2(j),j=1,ne2)
      endif
c
      do 50 i = nn, 1, -1
      if (llab) then
        if (lzero) rr = 0.0
        if (l9999) rr = 9999.0
        do 202 j = 1,ne
202     grow1(j) = rr
      elseif (libin1) then
        if (lupd) then
          if (in1.le.i.and.i.le.in2) call inrow(10,grow1,i,nn,ne)
        else
          call inrow(10,grow1,i,nn,ne)
        endif
      else
        read(10,*) (grow1(j),j=1,ne)
      endif
c
      if (in1.le.i.and.i.le.in2.and.mode.gt.0) then
        if (libin2) then
          call inrow(11,grow2,i-in1+1,nn2,ne2)
        else
          read(11,*) (grow2(j),j=1,ne2)
        endif
      endif
c
      do 41 j = 1, ne
        r1 = grow1(j)
        if (i.lt.in1.or.i.gt.in2.or.j.lt.ie1.or.j.gt.ie2) then
          r2 = 9999.0
        elseif (mode.gt.0) then
          r2 = grow2(j-ie1+1)
        endif
c
        goto (30,31,32,33,34,35,36,37,38,381,382,383),mode+1
c
c  one-grid modification
c
 30     if (r1.gt.9998.9) goto 39
        rr = r1*fact + bias
        goto 40
c
c  difference
c
 31     if (r1.gt.9998.9.or.r2.ge.9998.9) goto 39
        rr = r1 - r2
        goto 40
c
c  sum
c
 32     if (r1.gt.9998.9.or.r2.gt.9998.9) goto 39
        rr = r1 + r2
        goto 40
c
c  general sum
c
 33     if (r1.gt.9998.9.or.r2.gt.9998.9) goto 39
        rr = fak1*r1 + fak2*r2
        goto 40
c
c  treshold
c
 34     rr = r1
        if (lgt) then
          if (r2.gt.trsh) rr = 9999.0
        else
          if (r2.lt.trsh) rr = 9999.0
        endif
        goto 40
c
c  overwrite
c
 35     rr = r1
        if (r2.lt.9998.9d0) rr = r2
        goto 40
c
c  select where no data
c
 36     rr = r1
        if (r2.le.9998.9) goto 39
        goto 40
c
c  treshold select
c
 37     rr = 0
        if (r1.gt.trsh) rr = r1
        if (r1.le.trsh.and.r2.lt.trsh) rr = r2
        goto 40
c
c  variable multiplication
c
 38     if (r1.gt.9998.9) goto 39
        if (r2.lt.trsh) rr = r1*fak1
        if (r2.ge.trsh) rr = r1*fak2
        goto 40
c
c  N minus zeta
c
 381    if (r1.gt.9998.9.or.r2.gt.9998.9) goto 39
        rr = r1*r2/981000.d0
        goto 40
c
c  Helmert indirect effect
c  r2 statistics is correction
c
 382    if (r1.gt.9998.9.or.r2.gt.9998.9) goto 39
        r2 = -0.0559*r2**2/981000.d0
        rr = r1+r2
        goto 40
c
c  Geoid-quasigeoid corr
c
 383    if (r1.gt.9998.9.or.r2.gt.9998.9) goto 39
        r2 = -0.1119*r2**2/981000.d0
        rr = r1+r2
        goto 40
c
 39     rr = 9999.0
c
 40     grow1(j) = rr
        if (r1.lt.r1min) r1min = r1
        if (r1.gt.r1max.and.r1.lt.9998.9999) r1max = r1
        if (r1.lt.9998.9999) nr1 = nr1 + 1
        if (r1.lt.9998.9999) r1sum = r1sum + r1
        if (r1.lt.9998.9999) r1sum2 = r1sum2 + r1**2
        if (r2.lt.r2min) r2min = r2
        if (r2.gt.r2max.and.r2.lt.9998.9999) r2max = r2
        if (r2.lt.9998.9999) nr2 = nr2 + 1
        if (r2.lt.9998.9999) r2sum = r2sum + r2
        if (r2.lt.9998.9999) r2sum2 = r2sum2 + r2**2
        if (rr.lt.rrmin) rrmin = rr
        if (rr.gt.rrmax.and.rr.lt.9998.9999) rrmax = rr
        if (rr.lt.9998.9999) nrr = nrr + 1
        if (rr.lt.9998.9999) rrsum = rrsum + rr
        if (rr.lt.9998.9999) rrsum2 = rrsum2 + rr**2
 41   continue
c
c  output results row by row
c
      if (lobin) then
        if (lupd) then
          if (in1.le.i.and.i.le.in2)
     .    call outrow(10,grow1,i,nn,ne)
        else
          call outrow(20,grow1,i,nn,ne) 
        endif
      else
        if (lint) write(i20,51) (nint(grow1(j)),j=1,ne)
 51     format(30(/,12i6))
        if (.not.lint) write(i20,52) (grow1(j),j=1,ne)
 52     format(30(/,8f10.3))
      endif
 50   continue
c
      n = nn*ne
      if (nr1.gt.1) r1sum2 = sqrt((r1sum2-r1sum**2/nr1)/(nr1-1))
      if (nr1.le.1) r1sum2 = 0.0
      if (nr1.gt.0) r1sum = r1sum/nr1
      if (nr2.gt.1) r2sum2 = sqrt((r2sum2-r2sum**2/nr2)/(nr2-1))
      if (nr2.le.1) r2sum2 = 0.0
      if (nr2.gt.0) r2sum = r2sum/nr2
      if (nrr.gt.1) rrsum2 = sqrt((rrsum2-rrsum**2/nrr)/(nrr-1))
      if (nrr.le.1) rrsum2 = 0.0
      if (nrr.gt.0) rrsum = rrsum/nrr
      if (i20.eq.20) write(*,80) 
     .r1sum,r1sum2,r1min,r1max,n-nr1,ifile1,
     .r2sum,r2sum2,r2min,r2max,n-nr2,ifile2,
     .rrsum,rrsum2,rrmin,rrmax,n-nrr,ofile
80    format(' statistics of grids (in grid1 area, '
     .,'grid 2 poss. extended with 9999-values)',/,
     .'            mean   std.dev.    min     max  unknown',
     .'  filename',/,
     .' file1:',4f9.2,i7,3x,a24/,
     .' file2:',4f9.2,i7,3x,a24/' ofile:',4f9.2,i7,3x,a24)
      if (mode.eq.10) write(*,*) 
     .'(2nd row statistics above is ind.corr.)'
c
      if (lmulti) then
        lobinw = .false.
        lfirst = .false.
        close(10)
        open(10,file=ofile,form='unformatted',access='direct',recl=64)
        llab = .false.
        libin1 = .true.
        goto 8
      endif
      if (.not.llab) close(10)
      if (mode.ne.0) close(11)
      close(20)
      if (i20.eq.0.and.ifile1.eq.'9999') goto 55
c
90    end
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
      character*36 name
      ldirac = .true.
      open(unit=iunit,file=name,
     .status='old',form='unformatted',access='direct',recl=64,err=10)
      read(iunit,rec=1,err=11) icode
      if (icode.eq.777) return 
11    close(iunit)
10    ldirac = .false.
      open(unit=iunit,file=name,
     .status='old',form='formatted')
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
