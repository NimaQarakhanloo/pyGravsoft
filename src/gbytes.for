      program gbytes
c $Id: gbytes.for 251 2008-12-12 10:26:49Z cct $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                         g b y t e s
c
c   program for converting a grid in standard format into a file
c   of bytes. the values in the grid are mapped in the range
c   0 to 255 within the interval vmin to vmax. If more than
c   the actual number of points is requested the byte grid is
c   expanded with zeros. the program may also thin out a grid
c   and skip rows and columns to north and east.
c   a list of points may be read in and superimposed on the bytes,
c   list of points must be of form statno, lat, lon ...
c
c   (c) author: rene forsberg, kort- og matrikelstyrelsen, sep 89
c   updated jan 1992, rf
c   thin-out error corrected jun 1992, fr. Minor correction in do
c   loop jj by cct 2008-12-12.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character*1 byte(512,512)
      dimension d(2001)
      character*72 ifile,ofile,pfile
      character*1 ch, cp
      logical lpoint,lskip
c
      write(*,1)
1     format(' input: inputfile'/
     .       '        outputfile'/
     .       '        n_north, n_east (200,320 for PIPS)'/
     .       '        vmin, vmax (min and max range, 0,0 = auto)'/)
      read(*,'(a)') ifile
      read(*,'(a)') ofile
      read(*,*) nwn, nwe, v1, v2
      write(*,*) 'skip or thin data? (y/n) '
      read(*,'(a)') ch
      iskipn = 0
      iskipe = 0
      itn = 1
      ite = 1
      nwnn=nwn
      nwee=nwe
      if (ch.eq.'y'.or.ch.eq.'Y') then
        write(*,2)
2       format(' input: iskipn, iskipe (rows/colums to skip)'/
     .         '        ithinn, ithine (thin factors)'/)
        read(*,*) iskipn,iskipe,itn,ite
        if (itn.eq.0) itn = 1
        if (ite.eq.0) ite = 1
      endif
      write(*,*) 'overwrite a pointfile on grid? (y/n) '
      read(*,'(a)') ch
      lpoint = (ch.eq.'y'.or.ch.eq.'Y')
      if (lpoint) then
        write(*,*) 'point file name: '
        read(*,'(a)') pfile
        write(*,*) 'byte value (0-255)? '
        read(*,*) ib
        cp = char(ib)
        open(30,file=pfile,status='old')
      endif
c
      open(10,file=ifile,status='old')
      open(20,file=ofile,status='unknown')
c
c  find min and max
c
      read(10,*) rfi1,rfi2,rla1,rla2,dfi,dla
      if (abs(rfi1).gt.100.or.abs(rfi2).gt.100) read(10,*) iell,izone
      nn = (rfi2-rfi1)/dfi + 1.5
      ne = (rla2-rla1)/dla + 1.5
      rfiw2 = rfi2-iskipn*dfi
      rfiw1 = rfiw2 - (nwn-1)*dfi
      rlaw1 = rla1+iskipe*dla
      rlaw2 = rlaw1 + (nwe-1)*dla
      write(*,4) nn,ne,rfiw1,rfiw2,rlaw1,rlaw2,nwn,nwe
4     format(' --- GBYTES ---'/' No of points in original grid: ',2i5,
     ./' Byte grid covers area: ',4f12.4,
     ./' No of N-S and E-W pixels: ',2i5)
      if (nn.gt.2001.or.ne.gt.2001) stop 'grid too big - label error?'
c
c  find min and max (optional)
c
      if (v1.eq.0.and.v2.eq.0) then
        v1 = 9.9d9
        v2 = -9.9d9
        do 5 i = 1,nn
          read(10,*) (d(j),j=1,ne)
          do 5 j = 1, ne
            if (d(j).lt.v1) v1 = d(j)
            if (d(j).gt.v2) v2 = d(j)
5       continue
        rewind(10)
        read(10,*) rfi1,rfi2,rla1,rla2,dfi,dla
        if (abs(rfi1).gt.100.or.abs(rfi2).gt.100) read(10,*) iell,izone
        write(*,6) v1,v2
6       format(' Data min and max = ',2f12.2)
      endif
c
      imin = 9999
      imax = -9999
c
c  skip first rows
c
      do 8 i = 1, iskipn
8     read(10,*) (d(j),j=1,ne)
c
c  read wanted rows
c
      irow = (nwn-1)*itn + 1 + iskipn
      ii = 0
      do 11 i = iskipn+1, irow
        if (i.le.nn) then
          read(10,*) (d(j),j=1,ne)
        else
          do 9 j = 1,ne
9         d(j) = 0
        endif
        if(i.gt.(iskipn+nwn)) goto 12
        if (mod(i-iskipn,itn).ne.0) goto 11
        ii = ii+1
c
        lskip=.false.
        do 10 jj = 1, nwe
          j = iskipe + (jj-1)*ite + 1
          if(j.gt.(nwe+iskipe)) then
c changed 2008-12-12 by cct
            lskip=.true.
            goto 11
          endif
          if (j.le.ne) then
            ib = nint((d(j)-v1)/(v2-v1)*256)
	    if (ib.lt.0) ib = 0
	    if (ib.gt.255) ib = 255
	    if (ib.lt.imin) imin = ib
	    if (ib.gt.imax) imax = ib
	    byte(ii,jj) = char(ib)
          else
            byte(ii,jj) = char(0)
          endif
10      continue
11    continue
      if (lskip) jj=jj-1
12    continue
      if(ii.lt.nwn) nwnn=ii
      if(jj.lt.nwe) nwee=jj
c
c  read point list
c
      if (lpoint) then
        j = 0
15      read(30,*,end=16) i,rfi,rla
        ii = (rfiw2-rfi)/(dfi*itn) + 1.5
        jj = (rla-rlaw1)/(dla*ite) + 1.5
        if (1.le.ii.and.ii.le.nwnn.and.1.le.jj.and.jj.le.nwee) then
          byte(ii,jj) = cp
          j = j+1
        endif
        goto 15
16      write(*,*) 'No of points in byte area: ',j
      endif
c
c  output bytes
c
      write(20,19) ((byte(i,j),j=1,nwee),i=1,nwnn)
19    format(512(512a))
      write(*,20) imin,imax
20    format(' Min and max bytes output inside original grid: ',2i5)
      close(20)
      end
