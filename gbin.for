      program gbin
c $Id: gbin.for 243 2008-10-29 10:10:19Z cct $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                          g b i n
c
c  program for converting a text grid into binary or the reverse.
c  the grid file may only contain one grid. 
c
c  a direct access file format is used for binary files.
c  one record consist of 16 floating point (real*4) numbers = 64 bytes.
c  due to the use of real*4, binary files will not work on machines  
c  with a single precision word length different from 4 bytes, in 
c  this case the 'recl' must be changed in the open statements below. 
c
c  input:
c
c      ifile
c      ofile
c      mode (1: txt-> bin, 2: bin->txt, 3: bin->txt integers)
c
c  if ofile = '0' is specified no output is produced (useful for a check).
c
c  (c) rene forsberg oct 89, updated dec 89
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      character*128 ifile,ofile
      character*12 dirall,unfall
      logical libin,lobin,litxt,lotxt,lutm,lint
      real*4 grow(2000)
      data grow/2000*0.0/ 
      igrow =     2000
c
      write(*,10)
10    format(' input: inputfile'/'        outputfile'/
     .'        outputformat (1: bin, 2: txt, 3: txt int)')
      read(*,11) ifile
      read(*,11) ofile
11    format(a128)
      read(*,*) mode
      litxt = (mode.eq.1)         
      libin = (mode.ge.2)
      lotxt = (mode.ge.2)
      lobin = (mode.eq.1)
      lint = (mode.eq.3)
      if (ofile.eq.'0') lotxt = .false.
      if (ofile.eq.'0') lobin = .false.
c
      call openg(10,ifile,libin)
      litxt = (.not.libin)
      if (lobin) open(20,file=ofile,form='unformatted',status='unknown',
     .access='direct',recl=64)
      if (lotxt) open(20,file=ofile,status='unknown')
c
      if (libin) write(*,20)               
20    format(' ---  G B I N  ---'/' conversion from binary format')
      if (litxt) write(*,21)
21    format(' ---  G B I N  ---'/' conversion from text format')
c
      nunk = 0
      rmin = 999999.99
      rmax = -999999.99
      rsum = 0.0
      rsum2 = 0.0
c
      if (libin) then 
        read(10,rec=1) icode,rfi1,rfi2,rla1,rla2,dfi,dla,lutm,iell,izone
        if (icode.ne.777) stop 'binary check code 777 missing'
      else
        read(10,*) rfi1,rfi2,rla1,rla2,dfi,dla
        lutm =  (abs(rfi1).gt.100.or.abs(rfi2).gt.100)
        if (lutm.and.litxt) read(10,*) iell,izone
      endif
      nn = (rfi2-rfi1)/dfi + 1.5
      ne = (rla2-rla1)/dla + 1.5
      write(*,30)
30    format(' grid label:')
      if (.not.lutm) write(*,31) rfi1,rfi2,rla1,rla2,dfi,dla
      if (lutm) write(*,32) rfi1,rfi2,rla1,rla2,dfi,dla,iell,izone
31    format(' ',4f12.6,2f12.7)
32    format(' ',6f12.0,/,' ',2i4)
      write(*,33) nn,ne,nn*ne
33    format(' number of points in grid: ',2i7,i9)
      if (ne.gt.igrow) stop 'rows in grid too long - change igrow'     
c
      if (lobin) write(20,rec=1) 
     .777,rfi1,rfi2,rla1,rla2,dfi,dla,lutm,iell,izone
      if (lotxt.and.(.not.lutm)) write(20,31)
     .rfi1,rfi2,rla1,rla2,dfi,dla
      if (lotxt.and.lutm) write(20,32)
     .rfi1,rfi2,rla1,rla2,dfi,dla,iell,izone
c
      do 40 i = nn, 1, -1
        if (litxt) read(10,*) (grow(j),j=1,ne)
        if (libin) call inrow(10,grow,i,nn,ne)
        if (lotxt.and.lint) write(20,35) (nint(grow(j)),j=1,ne)
        if (lotxt.and.(.not.lint)) write(20,36) (grow(j),j=1,ne)
35      format(100(/,' ',12i6))
36      format(100(/,' ',8f9.3))
        if (lobin) call outrow(20,grow,i,nn,ne)
        do 40 j = 1, ne
          rr = grow(j)
          if (rr.gt.9998.9999) nunk = nunk+1
          if (rr.gt.9998.9999) goto 40 
          if (rr.lt.rmin) rmin = rr
          if (rr.gt.rmax) rmax = rr
          rsum = rsum + rr
          rsum2 = rsum2 + rr**2
40    continue
c
      n = nn*ne - nunk
      if (n.gt.1) rsum2 = sqrt((rsum2 - rsum**2/n)/(n-1))
      if (n.le.1) rsum2 = 0
      if (n.gt.0) rsum = rsum/n
      write(*,50) rsum,rsum2,rmin,rmax,nunk
50    format('      mean   std.dev.   min     max   unknown/9999',/,
     .' ',4f9.2,i9,/)
c
      close(20)
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
 10   write(iunit,rec=(irec+j)) (grow(j*16+k),k=1,16)
      return
      end
