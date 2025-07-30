      program stokes
c $Id: stokes.for 267 2009-01-27 14:36:10Z cct $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                          s t o k e s
c
c  program for integration of gravity data to geoid or deflectios
c  of the vertical. gravity data must be given in gridded format,
c  e.g. produced by 'geogrid'. the integrals are evaluated at given
c  station locations. in a 3 x 3 innerzone around each computation point
c  a local spline densification of the given data is made, analogous
c  to the terrain effect integration programme 'tc'.
c  gravity grid values 9999 or larger signals unknown values. such
c  cells are not taking part in the integration, but 9999's must not be
c  found in the innerzone (7 x 7 subgrid around computation points).
c
c  input:
c
c  gridfile,
c  stationfile,
c  outfile,
c  mode, lmean, psi1, psi2 (degrees)
c
c  where
c
c  mode = 1: stokes integration
c         2: vening-meinesz
c  lmean: true if mean of gravity data to be removed, false otherwise
c
c  the stationfile must contain the prediction point coordinates
c  in standard format (no, lat, lon, height, ..)
c
c  programmer: rene forsberg, danish geodetic institute,
c  vax version, university of new south wales, jan 1989
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      logical first, lmean, lalfa
      dimension g(200000)
      dimension cosfi(800),sinfi(800),cosla(800),sinla(800)
      dimension glab(6)
      dimension spldef(9), xdiv(18), ydiv(18)
      dimension a(7), r(7), q(7), splh(17,7)
c dimension changed 2009-01-25.
      data idim1,idim2 /120000, 800 /
      data  nspl / 9 /
      data  spldef / 0.0, 0.2, 0.4, 0.55, 0.7, 0.8, 0.9, 0.95, 0.98 /
      character*128 gfile,sfile,ofile
c
      pi = 3.141592654
      radeg = 180/pi
      gamma = 981500.0
c
c  input file names and input data
c
      write(*,2)
2     format(' input file names (gfile/sfile/ofile):')
      read(*,1) gfile
      read(*,1) sfile
      read(*,1) ofile
1     format(a128)
c
      open(10,file=gfile,form='formatted',status='old')
      open(20,file=sfile,form='formatted',status='old')
      open(30,file=ofile,form='formatted',status='unknown')
c
      write(*,3)
3     format(' input: mode (1:stoke, 2:vm), lmean, psi1, psi2')
      read(*,*) mode, lmean, psi1, psi2
c
      lalfa = (mode.eq.2)
      if (mode.lt.1.or.mode.gt.2) stop 'mode wrong'
      if (mode.eq.1) sf = 6371000.d0/(4*pi*gamma)
      if (mode.eq.2) sf = 206264.8/(4*pi*gamma)
c
      write(*,71) gfile,sfile,ofile,mode,psi1,psi2
71    format(' --- stokes ---'/,' integral formula evaluation from ',
     .'gridded data'/' grid file: ',a36/' point file: ',a36/
     .' output file: ',a36/
     .' mode: ',i2,', range',f9.4,' to',f9.4,' degrees'/)
c
c  read in gravity grid file
c
c  define internal coordinates with base at (rfic,rlac) at
c  sw corner of grid cell, i.e. one half grid unit sw of corner
c  point
c
      read(10,*) glab
      do 5 i = 1, 6
5     glab(i) = glab(i)/radeg
      rfi1 = glab(1)
      rfi2 = glab(2)
      rla1 = glab(3)
      rla2 = glab(4)
      dfi = glab(5)
      dla = glab(6)
      rfic = rfi1-dfi/2
      rlac = rla1-dla/2
      nn = (rfi2-rfi1)/dfi + 1.5
      ne = (rla2-rla1)/dla + 1.5
      if (nn*ne.gt.idim1) stop 'data array dimension too small'
      if (nn.gt.idim2.or.ne.gt.idim2)
     .stop 'sin and cos arrays too small'
c
      do 7 i = nn,1,-1
        i0 = (i-1)*ne
        read(10,*) (g(j+i0),j=1,ne)
7     continue
c
      call hinfo(rfic*radeg,rlac*radeg,dfi*radeg,dla*radeg,nn,ne,
     .g,idim1,gmean)
c
c  remove mean from gravity data if wanted
c
      n = nn*ne
      if (lmean) then
        do 72 i = 1, n
72      if (g(i).lt.9998.999) g(i) = g(i) - gmean
        write(*,73)
73      format(' mean of gravity data removed prior to integration')
      endif
c
c  convert psi to radians
c
      if (psi2.gt.180) psi2 = 180.0
      psi1 = psi1/radeg
      psi2 = psi2/radeg
      cpsi1 = cos(psi1)
      cpsi2 = cos(psi2)
      if (psi1.lt.1.0e-15) psi1 = -1.0e-15
c
c  set arrays storing cosine and sine of latitude and longitude
c
      do 8 i = 1, nn
        rfi = rfi1 + (i-1)*dfi
        cosfi(i) = cos(rfi)
        sinfi(i) = sin(rfi)
8     continue
      do 9 j = 1, ne
        rla = rla1 + (j-1)*dla
        cosla(j) = cos(rla)
        sinla(j) = sin(rla)
9     continue
c
      if (mode.eq.1) write(*,91)
91    format(/'                 geoid undulations (meter)',/,
     .'   stat     lat       lon       h      inner    total'/)
      if (mode.eq.2) write(*,92)
92    format(/'                 deflections of the vertical (arcsec)',/,
     .'   stat     lat       lon       h           inner'
     .,'             total'/)
c
c  central station loop
c
100   read(20,*,end=200) istat, sfi, sla, sh
      rsum = 0
      rsum1 = 0
      rfi = sfi/radeg
      rla = sla/radeg
      rrfi = (rfi - rfic)/dfi
      rrla = (rla - rlac)/dla
      cosrfi = cos(rfi)
      sinrfi = sin(rfi)
      cosrla = cos(rla)
      sinrla = sin(rla)
      dlacos = dla*cosrfi
c
c  find maximal computation area boundaries
c
      iis = ifrac(rrfi)
      jjs = ifrac(rrla)
      i1 = iis - psi2/dfi
      i2 = iis + psi2/dfi + 1
      if (i1.lt.0) i1 = 0
      if (i2.gt.nn) i2 = nn
      cmin = min(cosfi(i1+1),cosfi(i2))
      j1 = jjs - psi2/(dla*cmin)
      j2 = jjs + psi2/(dla*cmin) + 1
      if (j1.lt.0) j1 = 0
      if (j2.gt.ne) j2 = ne
c
c  find boundaries of spline inner zone
c
      ii1 = iis - 1
      ii2 = ii1 + 3
      jj1 = jjs - 1
      jj2 = jj1 + 3
      if (ii1.lt.2.or.ii2.gt.nn-2.or.jj1.lt.2.or.jj2.gt.ne-2) then
        write(*,101)
101     format(' *** warning: the following station outside area '
     .  ,'or too close to edge ***')
        ii1 = 0
        ii2 = 0
        jj1 = 0
        jj2 = 0
        goto 120
      endif
c
c  skip innerzone logic if psi1 big enough
c
      if (psi1.gt.2.5*sqrt(dfi**2+dlacos**2)) goto 120
c
c  computation of inner grid with spline densification
c
      npoint = 2*nspl - 1
      do 110 i = 1, nspl
        s = spldef(i)
        ydiv(i) = (rrfi-ii1)*s
        xdiv(i) = (rrla-jj1)*s
        ydiv(npoint+2-i) = 3.0 - (ii2-rrfi)*s
        xdiv(npoint+2-i) = 3.0 - (jj2-rrla)*s
  110  continue
c
c  compute data values at horizontal spline lines
c
       do 111 j = 1, 7
         do 112 i = 1, 7
           hh = g((i+ii1-3)*ne + j+jj1-2)
           a(i) = hh
           if (hh.gt.9998.999) then
             write(*,113) istat,sfi,sla
  113        format(' *** warning: stat ',i6,2f12.5,' not computed,'
     .   /'     unknown values in spline zone')
             goto 100
           endif
  112    continue
         call initsp(a, 7, r, q)
         do 111 i = 1, npoint
           ry = (ydiv(i+1)+ydiv(i))/2
           splh(i,j) = spline(ry+2.5, a, 7, r)
  111  continue
c
c  scan horizontal lines, start with center for innermost element
c  innerzone effects computed with psi in planar approximation
c
       first = .true.
       do 115 i = nspl, nspl+npoint-1
         ii = i
         if (i.gt.npoint) ii = i - npoint
         ry = (ydiv(ii+1)+ydiv(ii))/2
         do 116 j = 1, 7
           a(j) = splh(ii,j)
           call initsp(a, 7, r, q)
  116    continue
         do 115 j = nspl, nspl+npoint-1
         jj = j
         if (j.gt.npoint) jj = j - npoint
         rx = (xdiv(jj+1)+xdiv(jj))/2
         dg = spline(rx+2.5, a, 7, r)
c
         y = (ry+ii1 - rrfi)*dfi
         x = (rx+jj1 - rrla)*dlacos
         psi = sqrt(x**2+y**2)
         da = (ydiv(ii+1)-ydiv(ii))*dfi*
     .   (xdiv(jj+1)-xdiv(jj))*dlacos
         if (psi.le.psi1.or.psi.gt.psi2) goto 115
         if (lalfa) then
           if (psi.le.0) psi = 1.0e-37
           cosa = y/psi
           sina = x/psi
         endif
c
c  geoid by stokes, innermost effect by H&M 2-234
c
         if (mode.gt.1) goto 117
         if (first) then
           first = .false.
           rsum = (sqrt(da)*0.56*6371000.0)/gamma*dg/sf
         else
           rsum = rsum + dg*stoke(psi)*da
         endif
         goto 115
c
c  deflections by vening-meinesz, store neighbour values for slopes
c
  117    if (first) goto 118
         rr = dg*vm(psi)*da
         rsum = rsum + rr*cosa
         rsum1 = rsum1 + rr*sina
         goto 115
  118    first = .false.
  115  continue
c
c  compute rest of grid
c
  120 riz = rsum*sf
      if (mode.ge.2) riz1 = rsum1*sf
      do 121 i = i1+1, i2
        k0 = (i-1)*ne
        da = dfi*dla*cosfi(i)
        sinsin = sinfi(i)*sinrfi
        coscos = cosfi(i)*cosrfi
        do 121 j = j1+1, j2
          if (i.gt.ii1.and.i.le.ii2.and.j.gt.jj1.and.j.le.jj2)
     .    goto 121
c
          cosdla = cosla(j)*cosrla + sinla(j)*sinrla
          cpsi = sinsin + coscos*cosdla
          if (cpsi.ge.cpsi1.or.cpsi.lt.cpsi2) goto 121
          psi = acos(cpsi)
          dg = g(k0+j)
          if (dg.gt.9998.999) goto 121
c
          if (mode.eq.2) goto 122
          rsum = rsum + dg*stoke(psi)*da
          goto 121
c
  122     sindla = sinla(j)*cosrla - cosla(j)*sinrla
          sinpsi = sin(psi)
          cosa = (cosrfi*sinfi(i) - sinrfi*cosfi(i)*cosdla)/sinpsi
          sina = cosfi(i)*sindla/sinpsi
c
          rr = dg*vm(psi)*da
          rsum = rsum + rr*cosa
          rsum1 = rsum1 + rr*sina
  121 continue
c
      rsum = rsum*sf
      rsum1 = rsum1*sf
      if (mode.eq.1) write(*,130) istat,sfi,sla,sh,riz,rsum
      if (mode.eq.2) write(*,130) istat,sfi,sla,sh,riz,riz1,rsum,rsum1
      if (mode.eq.1) write(30,130) istat,sfi,sla,sh,rsum
      if (mode.eq.2) write(30,130) istat,sfi,sla,sh,rsum,rsum1
  130 format(' ',i6,2f10.5,f8.2,4f9.3)
      goto 100
c
200   close(30)
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      s t o k e
c
c  stokes function. cf. heiskanen and moritz 1967 (2-164)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function stoke(psi)
      implicit double precision(a-h,o-z)
      spsi2 = sin(psi/2)
      cpsi = cos(psi)
      stoke = 1/spsi2 - 6*spsi2 + 1 - 5*cpsi -
     .3*cpsi*log(spsi2+spsi2**2)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                            v m
c
c  vening-meinesz kernel function, see h&m 2-211
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function vm(psi)
      implicit double precision(a-h,o-z)
      cpsi2 = cos(psi/2)
      spsi2 = sin(psi/2)
      spsi = 2*cpsi2*spsi2
      vm = -0.5*cpsi2/spsi2**2 + 8*spsi - 6*cpsi2 - 3*(1-spsi2)/spsi
     .     + 3*spsi*log(spsi2 + spsi2**2)
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
      subroutine hinfo(fi0,la0,dfi,dla,nfi,nla,h,ihdim,rm)
c
      implicit double precision (a-h,o-z)
      double precision la0,la2
      dimension h(ihdim)
      logical lutm
c
      fi2 = fi0 + dfi*nfi
      la2 = la0 + dla*nla
      n = nfi*nla
      nz = 0
      n9999 = 0
      sum = 0.0
      sum2 = 0.0
      hmin = 9.9e9
      hmax = -9.9e9
c
      do 10 j = 1, nla
      do 10 i = 1, nfi
        hh = h((i-1)*nla+j)
        if (hh.gt.9998.999) goto 11
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
 9001 format(' --- data grid ',6f9.4)
      if (lutm) write(*,9010) fi0+dfi/2,fi2-dfi/2,la0+dla/2,
     .la2-dla/2,dfi,dla
9010  format(' --- selected ',6f9.0)
      write(*, 9002) nfi, nla, n, nz, n9999
 9002 format(' points: ',i3,' x ',i3,' = ',i6,', zero values:',i6,
     .', missing/9999:',i6)
      write(*, 9003) hmin, hmax, rm, rs
 9003 format(' min  max  mean  std.dev.:',4f11.2)
      write(*, 9004) h((nfi-1)*nla+1),h(nfi*nla),h(1),h(nla)
 9004 format(' corner values: ',4f9.2,/)
c
      return
      end
