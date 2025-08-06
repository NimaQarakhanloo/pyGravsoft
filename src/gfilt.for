      program gfilt
      implicit double precision(a-h,o-z)
c $Id: gfilt.for 181 2008-08-14 23:19:51Z tjansson $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                          g f i l t
c
c  program for simple filtering of a grid by a space-domain filter
c
c  input:
c
c  gridfile
c  outfile
c  mode, rdeg, lint
c
c  mode = 1: circular sharp cut-off
c         2: gaussian
c
c  where rdeg is the filtering parameter (full-width resolution)
c  lint: true for integer output
c
c  (c) Rene Forsberg, March 1996
c  min/max and 9999 update june 2004, rf
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character*72 ifile,ofile
      logical lint
      dimension f(1500,1900),g(1500,1900)
      nndim = 1500
      nedim = 1900
c
      write(*,3)
3     format(/
     .' ************************************************************',
     .'**'/
     .' *     GFILT - GRAVSOFT grid filter - vers. JUN04 (c) R',
     .'F/KMS  *'/
     .' ************************************************************',
     .'**')
      write(*,1)
1     format(' input file names - ifile/ofile: ')
      read(*,2) ifile 
      read(*,2) ofile
2     format(a72)
c
      write(*,4)
4     format(' input: mode (1:sharp,2:exp),rdeg, lint')
      read(*,*) mode, rdeg, lint
c
      open(10,file=ifile,status='old')
      open(20,file=ofile,status='unknown')
c
      read(10,*) rfi1,rfi2,rla1,rla2,dfi,dla
      nn = (rfi2-rfi1)/dfi + 1.5
      ne = (rla2-rla1)/dla + 1.5
      write(*,*) nn,ne 
      if (nn.gt.nndim.or.ne.gt.nedim) 
     .stop '*** grid too large, increase nndim or nedim ***'
c
      rdegh = rdeg/2
      rdegh2 = rdegh**2
      ii = 2*rdegh/dfi + 0.5
c    
      write(*,20) nn,ne,rdeg,-ii,ii
 20   format(/' ---  G F I L T  ---',/,
     .' number of points in grid,  north:',i7,', east:',i7/
     .' filter par (degrees): ,',f6.1,', operator length n/s: ',2i5)
      rmin = 9.d9
      rmax = -9.d9
      nr = 0
      rsum = 0
      rsum2 = 0
      n9999 = 0
c
      do 21 i = nn,1,-1
21    read(10,*) (f(i,j),j=1,ne)
c
      do 30 i = nn,1,-1
      rfi = rfi1 + (i-1)*dfi
      if (rfi.gt.89.9) then
        jj = nn
        jj = 75
      else
        jj = 2*rdegh/(dla*cos(rfi/57.29578d0)) + 0.5
      endif
      if (i.eq.nn.or.i.eq.1.or.i.eq.nn/2) then
        write(*,211) i,-jj,jj
211     format(' row ',i4,' operator length e/w ',2i6)
      endif
c
      do 30 j = 1, ne
	  wsum = 0
	  gsum = 0
        if (f(i,j).ge.9999) goto 26
	  cosfi = cos((rfi1+(i-1)*dfi)/57.29578d0)
        
        do 25 ip = -ii,ii
        do 25 jp = -jj,jj
          ik = i+ip 
	    if (ik.lt.1.or.ik.gt.nn) goto 25
	    jk = j+jp 
          if (jk.lt.1.or.jk.gt.ne) goto 25
	    r2 = ((i-ik)*dfi)**2 + ((j-jk)*(dla*cosfi))**2
c
          if (mode.eq.1) then
            if (r2.le.rdegh2) then
              ff = f(ik,jk)
              if (ff.lt.9999) then
	          wsum = wsum + 1
	          gsum = gsum + ff
              endif
            endif
	    else
            w = exp(-r2/rdegh2)
            ff = f(ik,jk)
            if (ff.lt.9999) then
              wsum = wsum + w
              gsum = gsum + ff*w
            endif
          endif
c
25      continue
26      if (wsum.eq.0) then
          g(i,j) = 9999.99
        else
          g(i,j) = gsum/wsum  
        endif
        if (g(i,j).ge.9999) then
          n9999 = n9999+1
        else
          nr = nr+1
          rsum = rsum + g(i,j)
          rsum2 = rsum2 + g(i,j)**2
          if (g(i,j).gt.rmax) rmax = g(i,j)
          if (g(i,j).lt.rmin) rmin = g(i,j)
        endif
30    continue
c
      write(20,32) rfi1,rfi2,rla1,rla2,dfi,dla
32    format(' ',4f12.6,2f12.7)
c       
      do 50 i = nn,1,-1
        if (lint) write(20,51) (nint(g(i,j)),j=1,ne)
51      format(30(/,12i6))
        if (.not.lint) write(20,52) (g(i,j),j=1,ne)
52      format(30(/,8f9.3))
50    continue
c
      close(20)
      sdev = 0.0
      if (nr.gt.1) sdev = sqrt((rsum2 - rsum**2/nr)/(nr-1))
      if (nr.gt.0) rsum = rsum/nr

      write(*,60) nn*ne,n9999,rsum,sdev,rmin,rmax
60    format(' number of points in output grid:',i6,', unknown:',i6,
     ./' mean,stddev,min,max: ',4f9.2) 
      end
