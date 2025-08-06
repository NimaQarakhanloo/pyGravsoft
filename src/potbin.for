      program potbin
      implicit double precision(a-h,o-z)
c $Id: potbin.for 251 2008-12-12 10:26:49Z cct $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                        P O T B I N
c
c  this program converts a potential coefficient file to binary
c  format or vice versa. coefficients are stored as real*8 variables,
c  in a direct access format without any sorting
c  (the direct access format ensures that the binary data may be
c  ported to e.g. dos without modifications)
c
c  input: filenames, lbin (t: txt->bin), lall (t: all coefficents)
c
c  if 'lall' is false only a subset of the coefficients is output,
c  further input (nmin, nmax) must be specified
c
c  (c) rf, national survey and cadastre of denmark, october 1990
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character*72 ifile,ofile
      integer*2 n,m
      logical lbin,ltxt,lall
c
      write(*,*) '--- POTBIN --- '
      write(*,*) 'input: ifile'
      write(*,*) '       ofile'
      write(*,*) '       lbin (t: txt->bin,  f: bin->txt), lall'
      read(*,1) ifile
      read(*,1) ofile
1     format(a72)
      read(*,*) lbin,lall
      ltxt = (.not.lbin)
      nmin = 0
      nmax = 9999
      if (.not.lall) then
        write(*,*) 'input: nmin, nmax'
        read(*,*) nmin,nmax
      endif
c
      if (lbin) open(10,file=ifile,form='formatted',status='old')
      if (ltxt) open(10,file=ifile,form='unformatted',
     .access='direct',recl=20,status='old')
c
      if (lbin) open(20,file=ofile,form='unformatted',
     .access='direct',recl=20,status='unknown')
      if (ltxt) open(20,file=ofile,form='formatted',status='unknown')
c
c  scan through coefficients
c
      nomax = 0
      i = 0
      j = 0
10    i = i+1
      if (lbin) read(10,*) n,m,cnm,snm
C correction 2008-12-12 by cct.
c     if (lbin) read(10,*,end=20) n,m,cnm,snm
      if (ltxt) read(10,rec=i) n,m,cnm,snm
      if (n.le.3) write(*,11) n,m,cnm,snm
      if (n.lt.nmin.or.n.gt.nmax) goto 10
      if (n.gt.nomax) nomax = n
      j = j+1
      if (lbin) write(20,rec=j) n,m,cnm,snm
      if (ltxt) write(20,11) n,m,cnm,snm
11    format(2i4,2g23.15)
      goto 10
c
20    write(*,11) n,m,cnm,snm
      write(*,30) i,j,nomax
30    format(' number of coefficient  pairs read: ',i6,', written: ',
     .i6,', nmax: ',i3)
      end
