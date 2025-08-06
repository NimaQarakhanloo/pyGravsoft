      program glist
      implicit double precision(a-h,o-z)
c $Id: glist.for 243 2008-10-29 10:10:19Z cct $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  converts grid format to line format.
c  values at 9999 or more are not converted (unknown).
c  rf dec 88, unsw, last change 2006-11-30 by cct.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character*128 ifile,ofile
      dimension glab(6)
      dimension grow(50000)
      write(*,5)
5     format(' input file names (ifile/ofile):')
      read(*,10) ifile
      read(*,10) ofile
10    format(a128)
      open(10,file=ifile,status='old')
      open(20,file=ofile,status='unknown')
c
      read(10,*) glab
      nn = (glab(2)-glab(1))/glab(5)+1.5
      ne = (glab(4)-glab(3))/glab(6)+1.5
      k = 0
      do 20 i = nn, 1, -1
        read(10,*) (grow(j),j=1,ne)
        do 20 j = 1, ne
          if (grow(j).gt.9998.99) goto 20
          k = k+1
          rfi = (i-1)*glab(5) + glab(1)
          rla = (j-1)*glab(6) + glab(3)
          write(20,19) k,rfi,rla,grow(j)
19        format(' ',i9,2f12.5,'  0.0',f9.3)
20    continue
      write(*,21) k
21    format(' number of points output: ',i7)
      stop 
      end
