      program tcgrid
c $Id: tcgrid.for 267 2009-01-27 14:36:10Z cct $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                          t c g r i d
c
c  program to produce a mean  elevation grid of a file containing
c  a digital terrain model (standard format or ngs modified).
c  the mean grid is obtained by simple averaging. areas with no
c  elevations are given value 9999. if an individual grid
c  element is covered partly by elevations, only these will be
c  averaged. if no data is available the unknown flag 9999 will be written.
c
c  if wanted, the mean elevation grid may be low-pass filtered
c  using a moving-average window of 'iffi' x 'ifla' cells.
c  if 'ifla' = -1 the longitude resolution is changed to match
c  the latitude variation
c
c  input:
c
c  inputfile,
c  outputfile,
c  itype, fimin, fimax, lamin, lamax, idfi, idla, iffi, ifla
c
c  where 'itype' specifies dtm-format type (0, 1: std, 2: ngs),
c  'fimin, fimax, lamin, lamax' specifies wanted area (degrees)
c  (if all are zero the entire grid is averaged),
c  'idfi' and 'idla' specifies mean grid average size (in units
c  of given grid - integers only), and
c  'iffi' and 'ifla' (odd numbers) specify averaging window
c  (= 0 or 1 means no low pass filtering).
c
c  note: grids are specified by a label (fi1, fi2, la1, la2,
c  dfi, dla - degrees), where the coordinate limits refer to
c  the  c e n t e r  of the outer blocks. this implies that
c  mean elevation grids will typically have somewhat "ackward"
c  labels.
c
c  values with integer part 9999 signals unknown values.
c  if more than half of the values for an average is unknown,
c  the average value itself is considered unknown and assigned 9999.
c  it is  n o t  allowed with 9999-values in filtering operations.
c
c  rene forsberg, june 1983. updated dec 7, 1983
c  nb: program written in fortran77
c  modified for cdc, with inclusion of filtering, jan 85
c  modified for real input/output and 9999-codes (unknown),
c  november 1988.
c
c  dimension changed from 12500 to 1000. Aug 2008 and to 6000 Jan 09.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision (a-h, l, o-z)
      dimension hhh(251)
      double precision elab(6)
      dimension       h(6000,6000), hc(6000,6000)
      logical lutm,lifla
      character*128 iname,oname
      data ihdim1, ihdim2/6000,6000 /
c
      write(*,1)
1     format(' --- TCGRID, GRAVSOFT grid averager and smoother ---'
     ./' input file name: ')
      read(*,777) iname
      write(*,101)
101   format(' output file name: ')
      read(*,777) oname
777   format(a128)
c
      open(20,file=iname,status='old')
      open(30,file=oname,status='unknown')
c
      write(*,2)
2     format(' input: itype (0 real, 1 integer),',
     .' fimin,fimax,lamin,',
     *'lamax (0 0 0 0 = all area)')
      write(*,2003)
2003  format('        idfi, idla ',
     *'(mean factors), iffi, ifla (filter factors) '/)
c
      read(*, *) itype, fi1, fi2, la1, la2, idfi, idla, iffi, ifla
c
      call relev1(20, itype, fi1, fi2, la1, la2, idfi,idla, lutm,ie,iz,
     .           fi0, la0, dfi, dla, nfi, nla, h, hc, ihdim1, ihdim2)
c
      if (.not.lutm) write(*,3) dfi*60,dla*60,iffi*dfi*60,ifla*dla*60
3     format(/' ofile grid spacing: ',2f8.1,', resolution: ',
     *2f8.1,' (arcmin)')
      if (ifla.eq.-1) write(*,*) 'Longitude resolution changing'
      call hinfo(fi0,la0,dfi,dla,nfi,nla,lutm,h,ihdim1,ihdim2)
      elab(1) = fi0+dfi/2
      elab(3) = la0+dla/2
      elab(2) = fi0 + nfi*dfi - dfi/2
      elab(4) = la0 + nla*dla - dla/2
      elab(5) = dfi
      elab(6) = dla
c
c  filtering of averaged grid
c
      if (iffi.eq.0) iffi = 1
      if (ifla.eq.0) ifla = 1
      if (iffi.eq.1.and.ifla.eq.1) goto 8000
c
      lifla = (ifla.eq.-1)
      ifla0 = ifla
      iffi2 = iffi/2
      ifla2 = ifla/2
      iffi = 2*iffi2+1
      ifla = 2*ifla2+1
c
      do 5 i = 1, nfi
      if (lifla) then
        cosfi = cos(((i-1)*dfi+fi0)/57.29578)
        ifla = iffi*dfi/(dla*cosfi)
        ifla2 = ifla/2
        ifla = 2*ifla2+1
        if (ifla.ne.ifla0) then
          write(*,*) 'Latitude band ',i,' from S set ifla =',ifla
          ifla0 = ifla
        endif
      endif
      do 5 j = 1, nla
        rh = 0
        ic = 0
        do 6 ii = -iffi2,iffi2
        do 6 jj = -ifla2,ifla2
          iii = i+ii
          jjj = j+jj
          if (iii.lt.1.or.iii.gt.nfi) goto 7
          if (jjj.lt.1.or.jjj.gt.nla) goto 7
          rr = h(iii,jjj)
          if (int(rr+0.00001).eq.9999)
     .    stop 'filtering of 9999-codes not allowed'
          rh = rh + rr
          ic = ic + 1
7         continue
6       continue
        hc(i,j) = rh/ic
5     continue
      do 8 i = 1, nfi
      do 8 j = 1, nla
        h(i,j) = hc(i,j)
8     continue
c
8000  if (.not.lutm) write(30,8001) elab
8001  format(' ',6f12.6/)
      if (lutm) write(30,8101) elab,ie,iz
8101  format(' ',6f12.0,/,' ',i5,i5)
      do 10 i = nfi, 1, -1
        do 9 j = 1, nla
          hhh(j) = h(i,j)
    9   continue
        write(*,*)i
        if (itype.eq.0) write(30,8002) (hhh(j),j=1,nla)
        if (itype.eq.1) write(30,8003) (nint(hhh(j)),j=1,nla)
8002    format(8f9.2)
8003    format(10(/,' ',12i6))
   10 continue
c
      end
c
      subroutine relev1(iunit,itype,fi1,fi2,la1,la2,idfi,idla,lutm,ie,
     .                 iz,fi0,la0,dfi,dla,nfi,nla,h,hc,ihdim1,ihdim2)
c
c  subroutine 'relev1' corresponds to subroutine 'rdelev', the
c  difference being that 'relev1' produces mean elevations.
c
      implicit double precision (a-h, l, o-z)
      dimension h(ihdim1, ihdim2), hc(ihdim1, ihdim2), hh(6000)
      double precision hlab(6), lat, lon
      logical lutm
c
c  defining parameters for index coordinate system (0,0) in sw-corner
c
      if (itype.le.1) then
        read(iunit,*) hlab
        lutm = (abs(hlab(1)).gt.100.or.abs(hlab(2)).gt.100)
        if (lutm) read(iunit,*) ie,iz
        if (lutm) write(*,8900) (hlab(i),i=1,6),ie,iz
 8900   format(/' --- utm ',6f10.0,i5,i4)
        if (.not.lutm) write(*,9000) (hlab(i),i=1,6)
 9000   format(/' --- ifile label ',6f10.4)
        hlab(1) = hlab(1)-hlab(5)/2
        hlab(2) = hlab(2)+hlab(5)/2
        hlab(3) = hlab(3)-hlab(6)/2
        hlab(4) = hlab(4)+hlab(6)/2
      else
        hlab(1) = 1.0*ifrac(fi1)-1.0/240
        hlab(3) = 1.0*ifrac(la1)-1.0/240
        hlab(5) = 1.0/120
        hlab(6) = 1.0/120
      endif
c
      if (.not.(
     .fi1.eq.0.and.fi2.eq.0.and.la1.eq.0.and.la2.eq.0)) goto 2
        fi1 = hlab(1)
        fi2 = hlab(2)
        la1 = hlab(3)
        la2 = hlab(4)
 2    dfi = hlab(5)*idfi
      dla = hlab(6)*idla
      i0 = ifrac((fi1-hlab(1))/hlab(5)+.001)
      j0 = ifrac((la1-hlab(3))/hlab(6)+.001)
      fi0 = i0*hlab(5) + hlab(1)
      la0 = j0*hlab(6) + hlab(3)
      nfi = ifrac((fi2-fi0)/dfi + .999)
      nla = ifrac((la2-la0)/dla + .999)
      infi = nfi*idfi
      inla = nla*idla
c
      do 10 j = 1, nla
      do 10 i = 1, nfi
        h(i,j) = 0
        hc(i,j) = 0
   10 continue
c
c  read standard elevation file
c
      if (itype.le.1) then
        ii = (hlab(2)-hlab(1))/hlab(5) + 0.5
        jj = (hlab(4)-hlab(3))/hlab(6) + 0.5
        infi0 = infi + i0
        jj1 = max(1, j0+1)
        jj2 = min(jj, j0+inla)
        if (jj1.gt.jj.or.jj2.lt.1) go to 20
c
        do 13 i = ii, 1, -1
          if (i.gt.infi0) then
            read(iunit,*) (hh(j),j=1,jj)
          else
            if (i.le.i0) go to 20
            read(iunit,*) (hh(j), j=1,jj)
            do 12 j = jj1, jj2
              im = (i-i0-1)/idfi+1
              jm = (j-j0-1)/idla+1
              if (int(hh(j)+0.00001).eq.9999) goto 12
              h(im,jm) = h(im,jm) + hh(j)
              hc(im,jm) = hc(im,jm) + 1.0
   12       continue
          endif
   13   continue
      else
c
c  read ngs format elevations
c
   14   read(iunit, end=20) lat, lon, hh
        ilat = ifrac((lat-fi0)*120 + 1.0d0)
        ilon = ifrac((359-lon-la0)*120 + 2.0d0)
        lon = 360.0 - lon
        write(*, 9001) lat,lon,hh(1),hh(2),hh(3),hh(4)
 9001   format(' --- ngs block ',2f8.3,4f7.0)
c
c  skip if record completely outside area
c
        if ((ilat+14).le.0.or.ilat.gt.infi.or.
     .  (ilon+119).le.0.or.ilon.gt.inla) go to 14
c
        k = 0
        do 15 i = 0, 14
        do 15 j = 1, 120
          k = k + 1
          ifi = ilat + i
          ila = ilon + 120-j
          if (ifi.le.infi.and.ila.le.inla.and.
     .    ifi.ge.1.and.ila.ge.1) then
            im = (ifi-1)/idfi + 1
            jm = (ila-1)/idla + 1
            h(im, jm) = h(im, jm) + hh(k)
            hc(im, jm) = hc(im, jm) + 1
          endif
   15   continue
        go to 14
      endif
c
c  grid scan(s) completed, perform averaging
c
c       nhalf = idfi*idla/2
c       if (nhalf.eq.0) nhalf = 1
   20   nhalf = 1
        write(*,*) 'Min number of data in cell for average:',nhalf
        do 21 i = 1, nfi
        do 21 j = 1, nla
          kk = nint(hc(i,j))
          if (kk.lt.nhalf) h(i,j) = 9999.0
          if (kk.ge.nhalf) h(i,j) = h(i,j)/kk
   21   continue
c
      return
      end
c
      subroutine hinfo(fi0,la0,dfi,dla,nfi,nla,lutm,h,ihdim1,ihdim2)
c
      implicit double precision (a-h, l, o-z)
      dimension h(ihdim1, ihdim2)
      logical lutm
c
      fi2 = fi0 + dfi*nfi
      la2 = la0 + dla*nla
      n = nfi*nla
      nz = 0
      n9999 = 0
      sum = 0
      sum2 = 0
      hmin = 32767
      hmax = -32767
c
      do 10 j = 1, nla
      do 10 i = 1, nfi
        hh = h(i, j)
        if (hh.eq.0) nz = nz+1
        if (int(hh+0.00001).ne.9999) goto 121
          n9999 = n9999 + 1
          goto 10
  121   if (hh.lt.hmin) hmin = hh
        if (hh.gt.hmax) hmax = hh
        sum = sum+hh
        sum2 = sum2 + (hh*1.0)**2
   10 continue
c
      rm = sum/n
      rs = -1.0
      if (n.gt.1) rs = sqrt((sum2 - sum**2/n)/(n-1))
      if (.not.lutm)
     *write(*, 9001) fi0+dfi/2,fi2-dfi/2,la0+dla/2,la2-dla/2,
     *dfi,dla
 9001 format(/' --- ofile gridinfo ',6f9.4)
      if (lutm)
     *write(*, 9101) fi0+dfi/2,fi2-dfi/2,la0+dla/2,la2-dla/2,
     *dfi,dla
 9101 format(/' --- ofile gridinfo ',6f9.0)
      write(*, 9002) nfi, nla, n, nz, n9999
 9002 format(' points: ',i3,' x ',i3,' = ',i6,', zero values:',i5
     .,', missing/9999:',i5)
      write(*, 9003) hmin, hmax, rm, rs
 9003 format(' hmin  hmax  mean  std.dev.:',4f11.2)
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  c
c                      i f r a c                                   c
c                                                                  c
c  subroutine giving true integer part of double precision real    c
c                                                                  c
c  rf, june 1983                                                   c
c                                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer function ifrac(r)
c
      implicit double precision (a-h,l,o-z)
      if (r.ge.0) then
        ifrac = r
      else
        ifrac = r - 0.999999
      endif
      return
      end
