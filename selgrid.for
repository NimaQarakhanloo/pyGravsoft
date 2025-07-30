      program selgrid
c $Id: selgrid.for 181 2008-08-14 23:19:51Z tjansson $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  program for selecting subgrid of binary or txt grid
c  output file is written on txt format (geographic grid only)
c
c  input:
c
c  <ifile> <ofile>
c  rfisw, rlasw, nn, ne
c
c  rf oct 87, cct July 1992 (only Text grids).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character *80 ifile, ofile 
      logical lbin, lgeo
      integer inz,outz,in,out
c
      dimension hrow(1500)
      dimension cha(254000)
      idim = 254000
c
      write(*,*)' SELGRID, Ver. JULY 1992.'
      write(*,*)' Input name of inputfile '
      read(*,'(A)')ifile
      write(*,*)' Input name of outputfile '
      read(*,'(A)')ofile
      inz=20
      outz=30
      in=5
      out=6
      write(*,*)' input latitude and longitude of SW corner (deg.)'
      write(*,*)' input number of steps in north and east ' 
      read(5,*) rfi, rla, nn, ne
c
      open(20,file=ifile,status='old')
      open(30,file=ofile)
c
      call rdgrid(inz,rfi,rla,nn,ne,dfi,dla,lgeo,lbin,cha,idim)
c
      write(outz, 20) rfi, rfi+(nn-1)*dfi, rla, rla+(ne-1)*dla,dfi,dla
c20   format(' < grid selgrid'/' geographic'/
c    .2(2(f14.6,'d'),/),2(f14.7,'d'),/)
 20   format(4f12.6,2f12.7)
c
      do 40 i = 1, nn
        j0 = (i-1)*ne
        do 30 j = 1, ne
30      hrow(j) = cha(j+j0)
        write(outz, 35) (hrow(j),j=1,ne)
35      format(/,65(' ',8f9.2,/))
40    continue
      write(outz, 50) 25
50    format(a6)
      close(outz)
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                       r d g r i d
c
c  subroutine for reading a digital grid file on
c  standard format, i.e. stored rowwise from nw to se, with label.
c  both binary and text files may be read
c
c  as se-corner coordinate 'rfic, rlac' (degrees) will be used
c  (unless they are zero, then the grid corner is used).
c  a grid containing 'inn' x 'ine' points will be put in array
c  'cha' of declared dimension 'idim'.
c  if inn=0 the complete grid will be read.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rdgrid(iunit, rfic, rlac, inn, ine, dfi, dla,
     .                  lgeo, lbin, cha, idim)
      dimension cha(idim)
c
      logical lgeo, lutm, lbin
c     integer*2 hlaB
      dimension hlab(9), hrow(1500)
      integer iunit
      integer *4 hlab1
      equivalence (hlab(1),hlab1)
c
c  initialize statistics
c
      hsum = 0.0
      hsum2 = 0.0
      hmin = 99999
      hmax = -99999
c
c  read label on gravity grid
c  readterlab stores label in hlab(4) ... hlab(9)
c
c  rc8000 special label, convert ellipsoid numbers
c  (1: wgs84, 2: hayford, 3: nad27)
c
c     read(iunit) ii
c     lbin = (ii.eq.7)
c     if (lbin) read(iunit) (hlab(ii), ii=1,9)
c     if (.not.lbin) j = readterlab(iunit, hlab)
      lbin=.false.
c
c     do 106 ii=4,9
c 106 hlab(ii) = hlab(ii)*180/3.14159265
c     j = hlab1.shift.24.shift.(-36)
c     lutm = (j.ne.1)
      lutm=.false. 
      j=100
      k = j/100
      iell = 1
      if (k.eq.1) iell = 2
      if (k.eq.6) iell = 3
      izone = mod(j, 100)
c
      read(20,*) (hlab(j),j=4,9)
c     do 106 ii=4,9
c 106 hlab(ii) = hlab(ii)*180/3.14159265
c     lutm = (abs(hlab(4)).ge.100.or.abs(hlab(5)).ge.100)
      lgeo = (.not.lutm)
      if (lgeo) goto 111
c     read(iunit,*) iell,izone
      write(6,110) iell,izone
      if (iell.lt.1.or.iell.gt.3.or.izone.lt.1.or.izone.gt.60)
     *stop
110   format(' - input grid in utm, ell ',i1,' zone ',i2,' -')
111   nn = (hlab(5)-hlab(4))/hlab(8)+1.5
      ne = (hlab(7)-hlab(6))/hlab(9)+1.5
c
c  find corner indices for wanted subgrid
c
      dfi = hlab(8)
      dla = hlab(9)
      if (inn.eq.0) goto 113
      if (.not.(rfic.eq.0.and.rlac.eq.0)) goto 114
        rfic = hlab(4)
        rlac = hlab(6)
114   ifi1 = (rfic-hlab(4))/dfi+1.5
      ila1 = (rlac-hlab(6))/dla+1.5
      ifi2 = ifi1+inn-1
      ila2 = ila1+ine-1
      rfic = (ifi1-1)*dfi + hlab(4)
      rlac = (ila1-1)*dla + hlab(6)
      goto 115
c
113   ifi1 = 1
      ila1 = 1
      ifi2 = nn
      ila2 = ne
      inn = nn
      ine = ne
      rfic = hlab(4)
      rlac = hlab(6)
c
c  check boundaries and array size
c
115   n = inn*ine
      if (1.le.ifi1.and.ifi2.le.nn.and.1.le.ila1.and.
     .ila2.le.ne) goto 121
      write(6, 119) igrid,nn,ne,ifi1,ifi2,ila1,ila2
119   format(' *** stop: wanted area too large',
     .7i4)
      stop
121   continue
      if (n.le.idim) goto 123
      write(6, 122) n,idim
122   format(' *** stop: array too small - wanted, declared ',2i7)
      stop
123   continue
c
c  read data grid values
c  rc8000 specials due to lack of free format
c
      do 130 i = nn,ifi1,-1
      read(iunit,*,end=199) (hrow(j),j=1,ne)
      do 130 j = 1,ne
        if (lbin) goto 133
c       k = readf(iunit,r)
c       if (k.eq.1) goto 131
        go to 131
199       write(6,132) i,j
132       format(' *** too few data in grid file, i,j = ',2i5)
          stop
131     continue
        r = hrow(j)
        goto 134
133     read(iunit) ii
        r = float(ii)
134     if (i.lt.ifi1.or.i.gt.ifi2) goto 130
        if (j.lt.ila1.or.j.gt.ila2) goto 130
        ii = (ifi2-i)*ine + j-ila1+1
        cha(ii) = r
130   continue
c
      do 135 ii = 1,n
        h = cha(ii)
        if (h.gt.hmax) hmax = h
        if (h.lt.hmin) hmin = h
        hsum = hsum + h
        hsum2 = hsum2 + h**2
135   continue
c
c  write information and statistics
c
      rfi = hlab(4) + (ifi1-1)*hlab(8)
      rla = hlab(6) + (ila1-1)*hlab(9)
      n=inn*ine
      r = hsum/n
      s = 0.0
      if (n.gt.1)
     .s = sqrt((hsum2 - hsum**2/n)/(n-1))
      if (lgeo) write(6,141) (hlab(j),j=4,7),hlab(8)*60,hlab(9)*60
      if (lutm) write(6,142) (hlab(j),j=4,9)
141   format(' input grid: ',4(f10.4,'d'),2(f7.2,'n'))
142   format(' input grid: ', 4f10.0,2f8.0)
      ii=ii-1
      if (lgeo) write(6, 143) rfi, rla, inn, ine, ii
143   format(' selected: sw corner ',2f10.4, ', points ', 3i6 )
      if (lutm) write(6, 144) rfi, rla, inn, ine, ii
144   format(' selected: sw corner ',2f10.0, ', points ', 3i6 )
      write(6, 145) hmin, hmax, r, s
145   format(' statistics: min max mean std.dev. ',4f9.2)
      return
      end
