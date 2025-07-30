      program altstack
c $Id: altstack.for 181 2008-08-14 23:19:51Z tjansson $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                        a l t s t a c k
c
c  Program to stack colinear tracks of altimeter data (e.g. Geosat).
c  the observations are interpolated on a time grid defined by
c  the best track and minimum and maximum times relative to a reference
c  time associated with each track. The tracks are merged in a free
c  adjustment using 1 cy/rev cosine and sine minimizing the differences.
c  The residual errors may be analyzed for covariance properties.
c
c  input:
c
c  <ofile>
c  listeach, lstat
c  ntracks, ipr
c  track1 
c
c  if lstat is true a statistical analysis is made of residuals,
c  otherwise only a listing is made in the output file.
c  ipr: print only every 'ipr' ERM in output (gives more clear output)
c
c  (c) programmed by rene forsberg, national survey and cadastre (denmark)
c  march 1989, university of new south wales vax
c
c  modified by Per Knudsen and Maria Brovelli, Nov. 1990.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-h,o-z)
      parameter(maxdat=50000,maxrev=22,maxrv4=26,maxnp=1000)
      integer*2 kerr,kcode
      common/field/ atime(maxdat),krev(maxdat),lat(maxdat),
     .  lon(maxdat),kalt(maxdat),kerr(maxdat),kcode(maxdat)
c
c  simple selection program for altimetry
c
      dimension itrack(maxrev),icount(maxrev),c1(maxrev),c2(maxrev),
     .          sig(maxrev),ts(maxrev),ttide(4)
      dimension hh(maxrev,maxnp),dd(maxrv4,maxnp),rlat(maxrev,maxnp),
     .          rlon(maxrev,maxnp),err(maxrev,maxnp)
      dimension x(maxnp),y(maxnp),ss(maxnp),sv(4,maxnp),ata(990),
     .cov(maxnp),scov(maxnp)
      dimension ic(maxnp),igroup(44)
      character*36 ofile
      logical lskip,lstat,lchang
      DATA DEGRAD/0.017453293D-6/
      igdim = 40000
c
      write(*,2)
2     format(' input: ofile,'/
     .' fi1,fi2,dfi,listeach,lstat',/
     .' ntracks,ipr, and'/' first track number ')
c
      read(*,1) ofile
1     format(a36)
c
c  geosat orbital frequency in seconds
c
      trev = 6037.5
      omega = 2*3.14159265d0 / trev
c
      open(20,file=ofile,status='unknown')
      open(21,file='dat.21',status='unknown')
      open(22,file='dat.22',status='unknown')
      open(23,file='dat.23',status='unknown')
c
      read(*,*) nlist,lstat,ntrack,ipr
      read(*,*) itrack(1)
      do j=2,ntrack
        itrack(j) = itrack(1)+(j-1)*244
      enddo
      irev1 = itrack(1)
      irev2 = itrack(ntrack)
c
      write(*,10)irev1,irev2
      if (lstat) write(20,10) irev1,irev2
10    format(/'1 --- A L T S T A C K ---'//,
     .' input file: erm01-22',/
     .' revolution numbers scanned in file:',i6,' to',i6)
c
c  read geoid grid
c
      do 15 i = 1,ntrack
      icount(i) = 0
      do 15 j = 1,maxnp
      dd(i,j) = 999.99
      err(i,j)= 999.99
15    hh(i,j) = 999.99
      nt = 0
      irevp = 0
      dtime=1.0
      dlim = 1.5
      it = 0
c
c  central loop for reading altimeter data
c  ---------------------------------------
c
      lcode=1
      ireg1=1
      ireg2=22
      ireg3=1
      call rreg(itrack(1),IREG1,IREG2,IREG3,ndat)
      write(*,*) 'ndat=',ndat
c
c find index of min and max latitude
c
      i0min=1
      i0max=1
      do 1001 i=1,ndat
      if(lat(i).lt.lat(i0min)) i0min=i
      if(lat(i).gt.lat(i0max)) i0max=i
 1001 continue
c
c find times for crossing between a nearly orthogonal
c track and each track
c 
      kmax=1
      next=1
      do 1020 k=1,ntrack
      imin=next
      imax=next
      do 1010 i=next,ndat
      if(krev(i).lt.itrack(k)) go to 1010
      if(krev(i).gt.itrack(k)) go to 1011
      icount(k)=icount(k)+1
      if(icount(k).gt.icount(kmax)) kmax=k
      if(lat(i).lt.lat(imin)) imin=i
      if(lat(i).gt.lat(imax)) imax=i
 1010 continue
 1011 continue
      next=i
      if(imin.eq.imax) go to 1020
      ddd=CROSSC(lat(imin),lon(imin),lat(imax),lon(imax),
     .           lat(i0min),lon(i0max),lat(i0max),lon(i0min),
     .           ALFA,BETA,latc,lonc)
      psi=dacos(dsin(lat(imin)*degrad)*dsin(lat(imax)*degrad)+
     .          dcos(lat(imin)*degrad)*dcos(lat(imax)*degrad)*
     .          dcos((lon(imin)-lon(imax))*degrad))
      psic=dacos(dsin(lat(imin)*degrad)*dsin(latc*degrad)+
     .           dcos(lat(imin)*degrad)*dcos(latc*degrad)*
     .           dcos((lon(imin)-lonc)*degrad))
      if(latc.lt.lat(imin)) psic=-psic
      ts(k)=psic/psi*(atime(imax)-atime(imin))+atime(imin)
      if(k.eq.kmax) dfi=psi/(atime(imax)-atime(imin))/degrad*1.0d-6
      write(20,*) k,(latc/1.0e6),(lonc/1.0e6),ts(k)
      write(*,*) k,(latc/1.0e6),(lonc/1.0e6),ts(k)
 1020 continue
      write(20,*) 'best track:',kmax,itrack(kmax),icount(kmax)
c
c  find min and max times relative to ts
c
      time1=1.0d20
      time2=0.0d0
      st=0.0d0
      nt=0
      krevp=0
      do 1030 i=1,ndat
      if(krev(i).ne.krevp) then
        do 1022 k=1,ntrack
        if(krev(i).eq.itrack(k)) go to 1023
 1022   continue
        go to 1030
 1023   continue
        endif
      rtime=atime(i)-ts(k)
      if(rtime.lt.time1) time1=rtime
      if(rtime.gt.time2) time2=rtime
      if(krev(i).eq.krevp.and.krev(i).eq.itrack(kmax)) then
        dt=atime(i)-atime(i-1)
        if(dt.lt.1.5) then
          nt=nt+1
          st=st+dt
          endif
        endif
      krevp=krev(i)
 1030 continue
      dtime=st/nt
      dfi=dabs(dfi*dtime)
      write(20,*) 'dtime=',dtime
c
c find reference time for best track
c
      do 1040 i=1,ndat
      if(krev(i).lt.itrack(kmax)) go to 1040
      if(krev(i).gt.itrack(kmax)) go to 1041
      rtime=atime(i)-ts(kmax)
      ii=rtime/dtime
      if(rtime.lt.0.0d0) ii=ii-1
      dt=rtime-ii*dtime
      if(dt.gt.0.0d0) go to 1041
 1040 continue
 1041 continue
      ii=time1/dtime+1.0
      if(rtime.lt.0.0d0) ii=ii-1
      time1=ii*dtime+dt
      np=(time2-time1)/dtime+1.5
      if (np.gt.maxnp) stop 'number of intervals too big'
      time2 = (np-1)*dtime + time1
      write(20,*) time1,time2,np
c
c start interpolating in time"grid" for each track
c
      write(*,16) nlist
      if (lstat) write(20,16) nlist
16    format(/' (a) listing of first and every',i4,' data points in',
     .' tracks'//,
     .'     lat       lon        dh      err      irev       time',/)
      nt=0
      do 20 ijk=1,(ndat-1)
      irev=krev(ijk)
c
      if (irev.lt.irev1) goto 20
      if (irev.gt.irev2) goto 50
c
c  see if track in list is different from previous
c
      if (irev.eq.irevp) goto 28
c
c  first point in new track
c
        irevp = irev
        do 25 i=1,ntrack
          if (irev.eq.itrack(i)) goto 27
25      continue
        lskip = .true.
        goto 20
c
 27     lskip = .false.
        timep = 999.99
        it = 0
c
 28   continue
      if (lskip) goto 20
      if (irev.ne.itrack(i)) stop 'program error'
c
c  stack geosat data in time intervals
c  -----------------------------------
c  find which 'cell' point belong to
c  and interpolate if the next observation is found
c  within 1.5 sec.
c
      rtime=atime(ijk)-ts(i)
      rfi=lat(ijk)/1.0d6
      rla=lon(ijk)/1.0d6
      if (rla.gt.180) rla = rla-360
      dh=kalt(ijk)/1.0d3
      sigma=kerr(ijk)/1.0d3
      ii = (rtime-time1)/dtime + 1.5
      if (ii.lt.1.or.ii.gt.np) goto 20
      timec = (ii-1)*dtime + time1
      d = timec-rtime
      if(dabs(d).lt.0.05) then
        hh(i,ii) = dh
        err(i,ii)= sigma
        rlat(i,ii) = rfi
        rlon(i,ii) = rla
        endif
      if((atime(ijk+1)-atime(ijk)).gt.dlim) go to 20
      ii = (rtime-time1)/dtime + 2.0d0
      if (ii.lt.1.or.ii.gt.np) goto 20
      rfi2=lat(ijk+1)/1.0d6
      rla2=lon(ijk+1)/1.0d6
      if (rla.gt.180) rla2 = rla2-360
      dh2=kalt(ijk+1)/1.0d3
      sigma2=kerr(ijk+1)/1.0d3
c
      timec = (ii-1)*dtime + time1
      d = timec-rtime
      alfa = d/(atime(ijk+1)-atime(ijk))
      hh(i,ii) = alfa*(dh2-dh) + dh
      err(i,ii)= alfa*(sigma2-sigma) + sigma
      rlat(i,ii) = alfa*(rfi2-rfi) + rfi
      rlon(i,ii) = alfa*(rla2-rla) + rla
c
c  reject outliers
c
      if (hh(i,ii).gt.12) hh(i,ii) = 999.99
c
      if (mod(it,nlist).eq.0)
     .write(*,40) rfi,rla,dh,sigma,irev,atime(ijk)
      if (mod(it,nlist).eq.0)
     .write(20,40) rfi,rla,dh,sigma,irev,atime(ijk)
40    format(' ',f9.4,f10.4,2f8.2,i6,f14.2)
      it = it+1
c
 20   continue
c
c  write out number of points per track
c
50    write(*,51) (itrack(j),icount(j),ts(j),j=1,ntrack)
      if (lstat) write(20,51) (itrack(j),icount(j),ts(j),j=1,ntrack)
51    format(//' (b) number of selected points and ts ',
     .'on wanted tracks: ',
     .25(/,' ',2i8,f11.1))
c
c  adjust cosine and sine terms of selected tracks to geoid
c  --------------------------------------------------------
c
       do 1101 i=1,(ntrack*2)
       y(i)=0.0d0
       ic(i)=1
       ist=(i*(i-1))/2
       do 1101 j=1,i
       ata(ist+j)=0.0d0
 1101 continue
c
      n=0
      sum=0.0d0
      do 1200 i=1,np
      rtime=(i-1)*dtime+time1
      cost=dcos(rtime*omega)
      sint=dsin(rtime*omega)
      cost2=cost*cost
      sint2=sint*sint
      scost=sint*cost
      do 1110 k=1,(ntrack-1)
      if(hh(k,i).lt.999.0) go to 1111
 1110 continue
      go to 1200
 1111 continue
      next=k+1
 1115 do 1120 l=next,ntrack
      if(hh(l,i).lt.999.0) go to 1121
 1120 continue
      go to 1200
 1121 continue
c
c  add observation of difference between track k and l
c  normal equation system
c
      kcos=(k-1)*2+1
      ksin=kcos+1
      lcos=(l-1)*2+1
      lsin=lcos+1
      ist=(kcos*(kcos+1))/2
      ata(ist)=ata(ist)+cost2
      ist=ist+kcos
      ata(ist)=ata(ist)+scost
      ist=ist+1
      ata(ist)=ata(ist)+sint2
      ist=(lcos*(lcos-1))/2+kcos
      ata(ist)=ata(ist)-cost2
      ist=ist+1
      ata(ist)=ata(ist)-scost
      ist=(lcos*(lcos+1))/2
      ata(ist)=ata(ist)+cost2
      ata(ist+kcos)=ata(ist+kcos)-scost
      ata(ist+ksin)=ata(ist+ksin)-sint2
      ata(ist+lcos)=ata(ist+lcos)+scost
      ata(ist+lsin)=ata(ist+lsin)+sint2
      dh=hh(k,i)-hh(l,i)
      y(kcos)=y(kcos)+dh*cost
      y(ksin)=y(ksin)+dh*sint
      y(lcos)=y(lcos)-dh*cost
      y(lsin)=y(lsin)-dh*sint
      n=n+1
      sum=sum+dh**2
      next=l+1
      go to 1115
 1200 continue
      write(20,*)'n and RMS of differences',n,dsqrt(sum/n)
c
c FIND SEPARATE GROUPS OF TRACKS THAT ARE NOT CONNECTED.
C THIS IS DONE BY CHECKING THE NORMAL EQUATIONS.
c
      modpar=2
      do 1407 i=1,(ntrack*MODPAR),MODPAR
 1407 igroup(i)=0
      ngroup=0
c
 1410 ngroup=ngroup+1
      do 1420 i=1,(ntrack*MODPAR),MODPAR
      if(igroup(i).eq.0) go to 1421
 1420 continue
 1421 istart=i
      igroup(istart)=ngroup
c
 1425 do 1441 i=istart,(ntrack*MODPAR),MODPAR
      lchang=.false.
      if(igroup(i).ne.ngroup) go to 1441
      do 1430 j=(i+modpar),(ntrack*MODPAR),MODPAR
      if(igroup(j).gt.0) go to 1430
      ij=(j*(j-1))/2+i
      if(ata(ij).gt.-1.0e-6) go to 1430
      lchang=.true.
      igroup(j)=ngroup
 1430 continue
c
      if (.not.lchang) go to 1441
      lchang=.false.
      do 1440 j=(i+modpar),(ntrack*MODPAR),MODPAR
      if(igroup(j).ne.ngroup) go to 1440
      ij0=(j*(j-1))/2
      do 1435 ii=istart,(j-MODPAR),MODPAR
      if(igroup(ii).gt.0) go to 1435
      ij=ij0+ii
      if(ata(ij).gt.-1.0e-6) go to 1435
      lchang=.true.
      igroup(ii)=ngroup
 1435 continue
 1440 continue
      if(lchang) go to 1425
 1441 continue
c
      do 1450 i=1,(ntrack*MODPAR),MODPAR
      IGROUP(I+1)=IGROUP(I)
      if(igroup(i).eq.0) go to 1410
 1450 continue
c
      write(6,999) ngroup
  999 format('0Number of separate groups',i4)
c
C ADD ZERO-MEAN-VALUE CONSTRAINT IN EACH GROUP (IF IUINC=0)
c
c add constraints: =0 when det=0, sum of c1(.)=0, and sum of c2=0
c
      ist=0
      do 1210 k=1,ntrack
      ist=ist+i
      kcos=(k-1)*2+1
      ksin=kcos+1
      ist=(kcos*(kcos+1))/2
      if(ata(ist).lt.1.0d-10) ata(ist)=ata(ist)+1.0d0
      det=ata(ist)*ata(ist+ksin)-ata(ist+kcos)**2
      if(det.lt.1.0d-10) ata(ist+ksin)=ata(ist+ksin)+1.0d0
      if(det.lt.1.0d-02) ata(ist+ksin)=ata(ist+ksin)+1.0d-3
 1210 continue
C
      DO 1455 igr=1,NGROUP
        do 1220 k=1,ntrack
        if(igroup(2*k).ne.igr) go to 1220
        kcos=(k-1)*2+1
        ksin=kcos+1
        istc=(kcos*(kcos-1))/2
        ists=istc+kcos
          do 1219 l=1,k
          if(igroup(2*l).ne.igr) go to 1219
          lcos=(l-1)*2+1
          lsin=lcos+1
          ata(istc+lcos)=ata(istc+lcos)+1.0d0
          ata(ists+lsin)=ata(ists+lsin)+1.0d0
 1219     continue
 1220   continue
 1455 CONTINUE
c
c solve system
c
      lred=.true.
      lbs=.true.
      call procnl(ata,ic,y,(ntrack*2),dummy,lred,lbs,
     .            990,maxnp,maxnp)
      do 1230 k=1,ntrack
      c1(k)=y((k-1)*2+1)
      c2(k)=y((k-1)*2+2)
 1230 continue
c
ccc      do 60 i = 1, ntrack
ccc        do 61 j  = 1, np
ccc          rtime=(j-1)*dtime+time1
ccc          x(j) = rtime
ccc          y(j) = hh(i,j)
ccc61      continue
ccc        call harm(trev,x,y,np,a,b,sigma)
ccc        c1(i) = a
ccc        c2(i) = b
ccc60    continue
c
c  adjust for cosine and sine terms
c
70    sumr = 0
      sumo = 0
      nr = 0
      do 72 i = 1, ntrack
        sum = 0
        n = 0
        do 8873 j = 1, np
          if (hh(i,j).gt.999) goto 8873
          n = n+1
8873      continue
         write(21,*) n
        n = 0
        do 73 j = 1, np
          tt = omega*((j-1)*dtime+time1)
          rtime = ((j-1)*dtime+time1)
          hmod = c1(i)*cos(tt)+c2(i)*sin(tt)
          if (hh(i,j).gt.999) goto 73
          n = n+1
          sum = sum + (hh(i,j)-hmod)**2
          sumo = sumo + hh(i,j)**2
          hh(i,j) = hh(i,j)-hmod
         write(21,*) rtime,hh(i,j)
73      continue
        if (n.ge.2) sig(i) = sqrt(sum/(n-1))
        sumr = sum + sumr
        nr = n + nr
72    continue
      write(*,75) (itrack(i),c1(i),c2(i),sig(i),
     .i=1,ntrack)
      write(20,75) (itrack(i),c1(i),c2(i),sig(i),
     .i=1,ntrack)
75    format(
     .//' (c) solved cosine and sine terms, sigma (m)'/,
     .40(/' ',i6,3f10.2))
      if (nr.gt.0) write(*,76) sqrt(sumo/nr),sqrt(sumr/nr)
      if (nr.gt.0) write(20,76) sqrt(sumo/nr),sqrt(sumr/nr)
76    format(
     .' overall r.m.s. fit of data to geoid before and',
     .' after slope corr:',/2f11.2)
      n=0
      sum=0.0d0
      do 1300 i=1,np
      rtime=(i-1)*dtime+time1
      cost=dcos(rtime*omega)
      sint=dsin(rtime*omega)
      cost2=cost*cost
      sint2=sint*sint
      scost=sint*cost
      do 1310 k=1,(ntrack-1)
      if(hh(k,i).lt.999.0) go to 1311
 1310 continue
      go to 1300
 1311 continue
      next=k+1
 1315 do 1320 l=next,ntrack
      if(hh(l,i).lt.999.0) go to 1321
 1320 continue
      go to 1300
 1321 continue
c
      dh=hh(k,i)-hh(l,i)
      n=n+1
      sum=sum+dh**2
      next=l+1
      go to 1315
 1300 continue
      write(20,*)'n and RMS of differences',n,dsqrt(sum/n)
c
c  write out results
c
90    write(*,91) time1,time2,dtime
      write(20,91) time1,time2,dtime
91    format(//' (d) interpolated heights (cm) at rel.times from ',
     .       f8.2,' to ',f8.2,/' spacing ',f8.4)
      write(*,92) (itrack(j),j=1,ntrack,ipr)
      write(20,92) (itrack(j),j=1,ntrack,ipr)
92    format('rel.t  ',20i6)
      do 100 i = 1, np
      write(*,101) time1+(i-1)*dtime,
     .             (nint(hh(k,i)*100),k=1,ntrack,ipr)
      write(20,101) time1+(i-1)*dtime,
     .             (nint(hh(k,i)*100),k=1,ntrack,ipr)
101   format(' ',f7.1,20i6)
100   continue
c
c  find residuals at node points relative to mean value
c  ----------------------------------------------------
c
      write(*,119)
      write(20,119)
119   format(/' mean and variance of residuals at same time:'
     ./5x,' time     lat      lon    ',
     .'     n      mean  std.dev. ')
      sums = 0.0
      do 120 j = 1, np
        sum = 0.0
        sum2 = 0.0
        serr = 0.0d0
        sumla = 0.0
        sumlo = 0.0
        n = 0
        do 121 i = 1, ntrack
          if (hh(i,j).gt.999) goto 121
          n = n+1
          sum = hh(i,j)+sum
          sum2 = hh(i,j)**2+sum2
          serr = err(i,j)**2+serr
          sumla = rlat(i,j)+sumla
          sumlo = rlon(i,j)+sumlo
121     continue
        if (n.eq.0) rlam = 999.99
        if (n.gt.0) rlam = sumla/n
        if (n.eq.0) rlom = 999.99
        if (n.gt.0) rlom = sumlo/n
        if (n.eq.0) vmean = 999.99
        if (n.gt.0) vmean = sum/n
        if (n.le.1) vsig = 25.0
        if (n.gt.1) vsig = sqrt((sum2-sum**2/n)/(n-1))
        if (n.le.0) sigma = 0.0
        if (n.gt.0) sigma = sqrt(serr/n)
        if (n.le.0) errm = 0.0
        if (n.gt.0) errm = sqrt((vsig**2+serr)/n)
        rg = vmean
        rtime=(j-1)*dtime+time1+ts(1)
c
       if(vmean.lt.999.0) write(22,*) ((j-1)*dtime+time1),vmean
         write(*,123) rtime,rlam,rlom, n, vmean, vsig,
     .      sigma,errm
         write(20,123) rtime,rlam,rlom, n, vmean, vsig,
     .      sigma,errm
123     format(' ',f10.1,f9.3,f9.3,i8,5f9.2)
        do 125 i = 1, ntrack
125     if (hh(i,j).lt.999) hh(i,j) = hh(i,j) - vmean
        sums = sums + vsig**2
120   continue
      write(*,126) sqrt(sums/np)
      write(20,126) sqrt(sums/np)
126   format(' overall r.m.s. sigma of residuals: ',f9.2)
      if (.not.lstat) goto 900
c
      write(*,128) (itrack(j),j=1,ntrack,ipr)
      write(20,128) (itrack(j),j=1,ntrack,ipr)
128   format(//' (f) residuals (cm) at each time',
     .' relative to mean value',/
     .'   lat  ',20i6)
      do 140 i = 1, np
      write(*,101) 
     .   time1+(i-1)*dtime,(nint(hh(k,i)*100),k=1,ntrack,ipr)
      write(20,101) 
     .   time1+(i-1)*dtime,(nint(hh(k,i)*100),k=1,ntrack,ipr)
140   continue
      do 1400 k=1,ntrack
      n=0
      do 1380 i=1,np
      if(hh(k,i).gt.999.0) go to 1380
      n=n+1
 1380 continue
      write(23,*) n
      do 1390 i=1,np
      if(hh(k,i).gt.999.0) go to 1390
      rtime=(i-1)*dtime+time1
      write(23,*) rtime,hh(k,i)
 1390 continue
 1400 continue
c
c  analyze covariance function of residual errors along track
c
      npc = np/2 + 1
      write(*,200)
      write(20,200)
200   format(//' (g) covariances of errors as a function of distance')
      do 199 j = 1, npc
        scov(j) = -1.0
        ic(j) = 0
199   continue
c
      do 210 i = 1, ntrack
        do 202 j = 1, np
202     x(j) = hh(i,j)*100
        call ecov(x,np,cov)
        do 208 j = 1, npc
          dd(i,j) = cov(j)
          if (cov(j).ge.99999) goto 208
          ic(j) = ic(j) + 1
          scov(j) = scov(j) + cov(j)
208     continue
210   continue
c
      write(*,215) (itrack(i),i=1,ntrac,ipr)
      write(20,215) (itrack(i),i=1,ntrac,ipr)
215   format('        ',20i6)
      do 225 j = 1, npc
        write(*,216) (j-1)*dfi,(nint(dd(i,j)),i=1,ntrack,ipr)
        write(20,216) (j-1)*dfi,(nint(dd(i,j)),i=1,ntrack,ipr)
216     format(1x,f7.1,20i6)
225   continue
c
      write(*,227)
      write(20,227)
227   format(/' average covariance function (cm**2):')
      do 228 j = 1,npc
        if (ic(j).eq.0) cv = 99999
        if (ic(j).gt.0) cv = scov(j)/ic(j)
        write(*,229) (j-1)*dfi,cv
        write(20,229) (j-1)*dfi,cv
229     format(' ',f7.1,f12.0)
228   continue
c
c  analyze time variations
c  -----------------------
c  analyze in multipla of 244 rev, which is the GEOSAT repeat period
c  nt is the number of time steps
c
      nt = ntrack/2-1
c
c  four main tidal periods (m2, s2, k1, o1) in units of days
c
      ttide(1) = 0.51752505d0
      ttide(2) = 0.5d0
      ttide(3) = 0.99726954d0
      ttide(4) = 1.07580587d0
c
      do 349 j = 1, np
        ss(j) = 99999.9
        do 349 i = 1,4
          sv(i,j) = 99999.9
349   continue
c
      do 350 k = 1, nt+4
        tr = ntrack*17.05/k
        if (k.gt.nt) tr = ttide(k-nt)
        do 310 j = 1, np
          n = 0
          sum = 0
          do 320 i = 1, ntrack
            if (hh(i,j).ge.999) goto 320
            n = n+1
            x(n) = ((j-1)*dtime+time1+ts(i))/86400
            y(n) = hh(i,j)*100
            sum = sum + y(n)**2
320       continue
          if (n.eq.0) dd(k,j) = 0.0
          if (n.eq.0) goto 310
          call harm(tr,x,y,n,a,b,vsig)
          dd(k,j) = a**2 + b**2
          if (dd(k,j).gt.99999) dd(k,j) = 99999
          if (k.eq.1) ss(j) = sum/n
          if (k.gt.nt) sv(k-nt,j) = vsig**2
310     continue
350   continue
      write(*,305) (nint(ntrack*17.05/k), k=1,nt)
      write(20,305) (nint(ntrack*17.05/k), k=1,nt)
305   format(//' (h) power spectrum (cm**2)'
     .,' as a function of time period: ',
     ./'   days:',20i6)
      do 353 k = 1,nt
      ic(k) = 0
353   x(k) = 0
      do 352 j = 1,np
        do 354 k = 1,nt
          if (dd(k,j).gt.99998) goto 354
          x(k) = x(k) + dd(k,j)
          ic(k) = ic(k) + 1 
354     continue
        write(*,351) (j-1)*dtime+time1,(nint(dd(k,j)),k=1,nt)
        write(20,351) (j-1)*dtime+time1,(nint(dd(k,j)),k=1,nt)
351     format(' ',f7.1,20i6)
352   continue
c
      do k = 1, nt
        if (ic(k).eq.0) x(k) = 99999.0
        if (ic(k).gt.0) x(k) = x(k)/ic(k)
      enddo 
      write(*,355) (nint(x(k)),k=1,nt)
      write(20,355) (nint(x(k)),k=1,nt)
355   format(/' average',20i6)
c
      write(*,360)
      write(20,360)
360   format(/' variances and power (cm**2) at tidal frequencies:'/
     .'         signal     m2    res     s2    res     k1 ',
     .'   res     o1    res')
      n = 0
      sum = 0
      do 359 i = 1,4
        c1(i) = 0
        c2(i) = 0
359   continue
c
      do 362 j = 1, np
        write(*,361) time1+(j-1)*dtime,
     .  nint(ss(j)),(nint(dd(nt+i,j)),nint(sv(i,j)),i=1,4)
        write(20,361) time1+(j-1)*dtime,
     .  nint(ss(j)),(nint(dd(nt+i,j)),nint(sv(i,j)),i=1,4)
361     format(' ',f7.1,9i7)
        if (ss(j).ge.99999) goto 362
        n = n+1
        sum = sum+ss(j)
        do 364 i=1,4
          c1(i) = c1(i) + dd(nt+i,j)
          c2(i) = c2(i) + sv(i,j)
364     continue
362   continue
      write(*,363) nint(sum/n),(nint(c1(i)/n),nint(c2(i)/n),i=1,4)
      write(20,363) nint(sum/n),(nint(c1(i)/n),nint(c2(i)/n),i=1,4)
363   format(' average',9i7)
c
c  covariance functions in time
c  ----------------------------
c  possibly only some are defined
c
      i1 = nint(itrack(1)/244.0)
      i2 = nint(itrack(ntrack)/244.0)
      ni = i2-i1+1
      nc = ni/2
      do 400 k = 1,nc
        scov(k) = 0.0
        ic(k) = 0
400   continue
c
      write(*,401) (nint((j-1)*17.05),j=1,nc)
      write(20,401) (nint((j-1)*17.05),j=1,nc)
401   format(//' (i) time covariance function (cm**2) '/
     .'   days:',20i6)
      do 410 j = 1, np
        do 420 k = 1,ni
420     x(k) = 999.99
        do 421 i = 1, ntrack
          ik = nint(itrack(i)/244.0)-i1+1
          if (ik.gt.ni) stop 'something wrong in ik'
          x(ik) = hh(i,j)*100
421     continue
        call ecov(x,ni,cov)
        write(*,422) (j-1)*dtime+time1,(nint(cov(k)),k=1,nc)
        write(20,422) (j-1)*dtime+time1,(nint(cov(k)),k=1,nc)
422     format(' ',f7.1,20i6)
        do 430 k = 1,nc
          if (cov(k).gt.99998) goto 430
          scov(k) = scov(k)+cov(k)
          ic(k) = ic(k)+1
430     continue
410   continue
c
      do 451 k = 1,nc
451   if (ic(k).eq.0) ic(k) = 1.0
      write(*,452) (nint(scov(k)/ic(k)),k=1,nc)
      write(20,452) (nint(scov(k)/ic(k)),k=1,nc)
452   format(/' average',20i6)
c
900   end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  c
c                           b i l i n                              c
c                                                                  c
c  interpolates values in an array a using bilinear                c
c  (parabolic hyperboloid) interpolation.                          c
c                                                                  c
c  parameters:                                                     c
c                                                                  c
c  bilin       interpolated value                                  c
c                                                                  c
c  ri, rj      interpolation argument, (1,1) in lower left corner, c
c              (imax, jmax) in upper right.                        c
c                                                                  c
c  a           integer*2 array with arguments                      c
c                                                                  c
c  imax, jmax  number of points in grid                            c
c                                                                  c
c  iadim       declared dimension of 'a'
c                                                                  c
c  outside area covered by 'a' the function returns the value of   c
c  the nearest boundary point, but only to half gridspacing.       c
c                                                                  c
c  programmer:                                                     c
c  rene forsberg, july 1983, updated mar 89                        c
c                                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function bilin(ri,rj,a,imax,jmax,iadim)
      implicit double precision(a-h, o-z)
      dimension a(iadim)
c
      if (ri.lt.0.49) stop 'bilin: point too far south of grid'
      if (ri.gt.imax+0.51) stop 'bilin: point too far north of grid'
      if (rj.lt.0.49) stop 'bilin: point too far west of grid'
      if (rj.gt.jmax+0.51) stop 'bilin: point too far east of grid'
      in = int(ri+1000)-1000
      ie = int(rj+1000)-1000
      rn = ri - in
      re = rj - ie
c
      if(in.lt.1) go to 1
      go to 2
 1      in = 1
        rn = 0.0
      go to 5
 2    if(in.ge.imax) go to 3
      go to 5
 3      in = imax-1
        rn = 1.0
 5    continue
      if(ie.lt.1) go to 4
      go to 7
 4      ie = 1
        re = 0.0
       go to 6
 7    if(ie.ge.jmax) go to 8
      go to 6
 8      ie = jmax-1
        re = 1.0
 6    continue
c
      k = (in-1)*jmax + ie
      bilin = (1-rn)*(1-re)*a(k) +
     .rn*(1-re)*a(k+jmax) + (1-rn)*re*a(k+1) +
     .rn*re*a(k+jmax+1)
      return
      end
c
      subroutine ecov(x,n,cov)
      implicit real*8(a-h,o-z)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                          e c o v
c
c  finds empirical covariances of data in array x(1...n), at equi-
c  distance = data spacing. unknown data are .gt. 999
c  rene forsberg, mar 89, unsw vax
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      parameter(maxnp=1000)
      dimension x(maxnp),cov(maxnp),ic(maxnp)
      do 10 i = 1, n
        ic(i) = 0
        cov(i) = 0.0
10    continue
      do 20 i = 1, n
        if (x(i).ge.999) goto 20
        do 21 j = i, n
          if (x(j).ge.999) goto 21
          id = j-i+1
          ic(id) = ic(id)+1
          cov(id) = cov(id) + x(i)*x(j)
21      continue
20    continue
      do 30 i = 1, n
        if (ic(i).eq.0) cov(i) = 99999.1
        if (ic(i).gt.0) cov(i) = cov(i)/ic(i)
30    continue
      return
      end
c
      subroutine harm(t,x,y,n,a,b,sig)
      implicit real*8(a-h,o-z)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                         h a r m
c
c  fits a harmonic of form a*cos(2*pi/t*x(i)) + b*sin(2*pi/t*x(i)) on
c  data x(i), y(i), i = 1,..n, using least squares.
c  t is the period (in same units as x(i)),
c  sig is the aposteriori variance of fit.
c  y-values at 999 or more are treated as unknown, and not fitted.
c
c  programmer: rene forsberg, march 89, unsw vax
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      parameter(maxnp=1000)
      dimension x(maxnp),y(maxnp)
      dimension c(maxnp),s(maxnp)
c
c  form normals
c
      cc = 0.0
      cs = 0.0
      ss = 0.0
      cy = 0.0
      sy = 0.0
      nn = 0
      do 10 i = 1, n
        if (y(i).ge.999) goto 10
        nn = nn+1
        r = 2*3.14159265d0/t*x(i)
        sinr = sin(r)
        cosr = cos(r)
        cc = cc + cosr**2
        cs = cs + cosr*sinr
        ss = ss + sinr**2
        cy = cy + cosr*y(i)
        sy = sy + sinr*y(i)
        c(i) = cosr
        s(i) = sinr
10    continue
c
      if (nn.gt.0) goto 11
      a = 99999.0
      b = 0.0
      return
c
c  solve
c
11    det = cc*ss - cs**2
      if (det.eq.0) then
        write(*,*) '*** warning *** harmfit: determinant equals zero'
        a = 99999.0
        b = 0.0
        return
      endif
      a = -(cs*sy - ss*cy)/det
      b = (cc*sy - cs*cy)/det
c
c variance
c
      sum = 0
      do 20 i = 1,n
        if (y(i).gt.999) goto 20
        sum = sum + (y(i) - a*c(i)-b*s(i))**2
20    continue
      if (nn.eq.2) sig = 0.0
      if (nn.gt.2) sig = sqrt(sum/(nn-2))
c
      return
      end
c
c*********************************************************
c
      subroutine rreg(itrack1,IREG1,IREG2,IREG3,ndat)
c
c  Subroutine to read altimeter data from alt.reg. belonging
c  to specific revolution numbers, which is itrack1+p*244
c
c  This version is tuned for GEOSAT obtained from OSU.
c  Each ERM is stored in a separate file.
c
c  Per Knudsen, Aug 6, 1990.
c  KMS - Denmark.
c
c  This version 14.11.90, PK.
c
      implicit real*8(a-h)
      parameter(maxdat=50000)
      integer*2 kerr,kcode
      common/field/ atime(maxdat),krev(maxdat),lat(maxdat),
     .  lon(maxdat),kalt(maxdat),kerr(maxdat),kcode(maxdat)
      character*10 rfile(22)
c
c  Altimeter data are stored 'sequentially' on a direct access
c  file.
c  A catalogue and a key in each record describe how data inside
c  subareas are found. The catalogue are stored in the beginning
c  of the file together with:
c
c   nlngth    length of reg.file in records
c   ntotal    total number of data
c   nhead     length of headder in records
c   time1     start time of data
c   time2     end time of data
c   irev1     first revolution number
c   irev2     last revolution number
c   ilat0     lower latitude bound of area
c   ilon0     lower longitude bound of area
c   idlat     size in latitude of subareas
c   idlon     size in longitude of subareas
c   idimi     number of subareas in latitude
c   idimj     number of subareas in longitude
c
c  The catalogue contains for subarea (i,j)
c
c   istore(1,i,j)  first record
c   istore(2,i,j)  last record
c   istore(3,i,j)  number of records in subarea
c
      real*8 time1,time2
      integer*4 irev1,irev2
      integer*4 istore(3,8,15)
c
c  Description of data record:
c
c   time       time in seconds from the beginning of 1985
c   irev       revolution number of track
c   ilat       latitude
c   ilon       longitude
c   issh       sea surface height after orbit and geophys. corr.
c   igeoid     OSU89B geoid heights
c   isstl      Sea surface topography from modified Levitus
c   isstd      Sea surface topography from OSU89 model to deg.10
c   istd       Standard deviation of issh
c   irdr       Orbit correction that has been applied
c   itide      Ocean tide, =32767 when no information
c   ifree1     Not used =32767
c   ifree2     Not used =32767
c   icode      =1 decending
c              =2 ascending
c              +100 std>0.1 m
c   next       number of next record inside sub-area
c              (=0: end)
c
      real*8 time
      integer*4 irev,ilat,ilon,issh,igeoid,next,next2
      integer*2 isstl,isstd,istd,irdr,itide,ifree1,ifree2,icode
      integer*2 i2rec(24)
      equivalence (i2rec(1),time),
     .            (i2rec(5),irev),
     .            (i2rec(7),ilat),
     .            (i2rec(9),ilon),
     .            (i2rec(11),issh),
     .            (i2rec(13),igeoid),
     .            (i2rec(15),next),
     .            (i2rec(17),isstl),
     .            (i2rec(18),isstd),
     .            (i2rec(19),istd),
     .            (i2rec(20),irdr),
     .            (i2rec(21),itide),
     .            (i2rec(22),ifree1),
     .            (i2rec(23),ifree2),
     .            (i2rec(24),icode)
      data rfile/'erm01.reg ',
     .           'erm02.reg ',
     .           'erm03.reg ',
     .           'erm04.reg ',
     .           'erm05.reg ',
     .           'erm06.reg ',
     .           'erm07.reg ',
     .           'erm08.reg ',
     .           'erm09.reg ',
     .           'erm10.reg ',
     .           'erm11.reg ',
     .           'erm12.reg ',
     .           'erm13.reg ',
     .           'erm14.reg ',
     .           'erm15.reg ',
     .           'erm16.reg ',
     .           'erm17.reg ',
     .           'erm18.reg ',
     .           'erm19.reg ',
     .           'erm20.reg ',
     .           'erm21.reg ',
     .           'erm22.reg '/
c
      write(*,555)
  555 FORMAT('1')
c
      ndat=0
c
      do 999 ireg=IREG1,IREG2,IREG3
      irevc=itrack1+(ireg-1)*244
c
      write(*,495) rfile(ireg),irevc
  495 format(/,' Reg.file: ',a10,', revolution no.',i5)
      open(12,file=rfile(ireg),access='direct',recl=48,status='old')
c
c  read headder and catalogue
c
      jrec=1
      read(12,rec=jrec) nlngth,ntotal,nhead,time1,time2,irev1,irev2
c     write(*,*) ' rev.s and ntotal  ',irev1,irev2,ntotal
c     jrec=2
c     read(12,rec=jrec) ilat0,ilon0,idlat,idlon,idimi,idimj
c     do 5 i=1,idimi
c     do 5 j=1,idimj
c     jrec=jrec+1
c     read(12,rec=jrec)(istore(k,i,j),k=1,3)
c   5 continue
c
      idat=0
      icount=0
c
c  Start reading in current register file
c
      do 40 irec=(nhead+1),nlngth
c
      read(12,rec=irec) i2rec
      icount=icount+1
c
      if(irev.lt.irevc) go to 30
      if(irev.gt.irevc) go to 998
      if(icode.gt.2) go to 30
         idat=idat+1
         ndat=ndat+1
         krev(ndat)=irev
         lat(ndat)=ilat
         lon(ndat)=ilon
         kalt(ndat)=issh-igeoid-isstd
         kerr(ndat)=istd
         kcode(ndat)=icode
         atime(ndat)=time
c
   30 continue
   40 continue
c
  998 write(*,510) idat
  510 format(50x,'number of data:',i6)
      close(12)
  999 continue
      return
      end
C
C   *****************************************************************
C
      FUNCTION CROSSC(I1,J1,I2,J2,I3,J3,I4,J4,ALFA,BETA,P1,P2)
C
C    THIS FUNCTION CALCULATES THE CROSS-OVER POSITION (P1,P2)
C    IN X,Y,Z COORDINATES USING A GREAT CIRCLE APPROXIMATION.
C
C    PROGRAMMED BY PER KNUDSEN.
C
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
C
      INTEGER*4 P1,P2
C
      DATA DEGRAD/0.017453293D-6/
C
C  COMPUTE X,Y,Z OF FIRST POINT OF THE NORTH GOING TRACK
C
      COSI1=DCOS(I1*DEGRAD)
      XN0=COSI1*DCOS(J1*DEGRAD)
      YN0=COSI1*DSIN(J1*DEGRAD)
      ZN0=DSIN(I1*DEGRAD)
C
C  AND THE UNIT VECTOR FROM THE FIRST TO THE LAST POINT
C
      COSI2=DCOS(I2*DEGRAD)
      XN1=COSI2*DCOS(J2*DEGRAD)-XN0
      YN1=COSI2*DSIN(J2*DEGRAD)-YN0
      ZN1=DSIN(I2*DEGRAD)-ZN0
      RN1=DSQRT(XN1*XN1+YN1*YN1+ZN1*ZN1)
      XN1=XN1/RN1
      YN1=YN1/RN1
      ZN1=ZN1/RN1
C
C DITTO FOR THE SOUTH GOING TRACK
C
      COSI3=DCOS(I3*DEGRAD)
      XS0=COSI3*DCOS(J3*DEGRAD)
      YS0=COSI3*DSIN(J3*DEGRAD)
      ZS0=DSIN(I3*DEGRAD)
      COSI4=DCOS(I4*DEGRAD)
      XS1=COSI4*DCOS(J4*DEGRAD)-XS0
      YS1=COSI4*DSIN(J4*DEGRAD)-YS0
      ZS1=DSIN(I4*DEGRAD)-ZS0
      RS1=DSQRT(XS1*XS1+YS1*YS1+ZS1*ZS1)
      XS1=XS1/RS1
      YS1=YS1/RS1
      ZS1=ZS1/RS1
C
C  NOW FIND A VECTOR IN THE CROSS SECTION BETWEEN THE TWO PLANES
C  SO (XN0,YN0,ZN0)+A*(XN1,YN1,ZN1)=B*(XS0,YS0,ZS0)+C*(XS1,YS1,ZS1)
C
C  OR  (XS0 XS1 -XN1) (B) (XN0)
C      (YS0 YS1 -YN1)*(C)=(YN0)    OR   FX=Y  AND X=(FTF-1)(FTY)
C      (ZS0 ZS1 -ZN1) (A) (ZX0)
C
C     FTF11=XS0*XS0+YS0*YS0+ZS0*ZS0
      FTF11=1.0D0
      FTF12=XS0*XS1+YS0*YS1+ZS0*ZS1
      FTF13=-XS0*XN1-YS0*YN1-ZS0*ZN1
C     FTF22=XS1*XS1+YS1*YS1+ZS1*ZS1
      FTF22=1.0D0
      FTF23=-XS1*XN1-YS1*YN1-ZS1*ZN1
C     FTF33=XN1*XN1+YN1*YN1+ZN1*ZN1
      FTF33=1.0D0
      FTF21=FTF12
      FTF31=FTF13
      FTF32=FTF23
      FTY1=XS0*XN0+YS0*YN0+ZS0*ZN0
      FTY2=XS1*XN0+YS1*YN0+ZS1*ZN0
      FTY3=-XN1*XN0-YN1*YN0-ZN1*ZN0
C
C REDUCTION OF COLUMNS
C
      CROSSC=0.0
      R=FTF21/FTF11
      FTF21=FTF21-R*FTF11
      FTF22=FTF22-R*FTF12
      FTF23=FTF23-R*FTF13
      FTY2=FTY2-R*FTY1
      R=FTF31/FTF11
      FTF31=FTF31-R*FTF11
      FTF32=FTF32-R*FTF12
      FTF33=FTF33-R*FTF13
      FTY3=FTY3-R*FTY1
      IF(FTF22.EQ.0.0) RETURN
      R=FTF32/FTF22
      FTF32=FTF32-R*FTF22
      FTF33=FTF33-R*FTF23
      FTY3=FTY3-R*FTY2
      IF(FTF33.EQ.0.0) RETURN
      CROSSC=FTF11*FTF22*FTF33
      A=FTY3/FTF33
      C=(FTY2-FTF23*A)/FTF22
C     B=(FTY1-FTF12*C-FTF13*A)/FTF11
C
      XC=XN0+A*XN1
      YC=YN0+A*YN1
      ZC=ZN0+A*ZN1
      R=DSQRT(XC*XC+YC*YC+ZC*ZC)
      COSAP1=DSQRT(1.0D0-(ZC/R)**2)
      AP1=DASIN(ZC/R)/DEGRAD
      AP2=DASIN(YC/(R*COSAP1))/DEGRAD
      ALFA=A/RN1
      BETA=C/RS1
      P1=AP1+0.5
      P2=AP2+0.5
C
      RETURN
      END
C
C   *****************************************************************
C
      SUBROUTINE PROCNL(AN,INUL,H,NT,VAR,LRED,LBS,IANT,INULT,IHT)
C
C     THIS SUBROUTINE USES A CHOLESKY ALGORITHME FOR REDUCING
C     AND SOLVING THE SYSTEM OF LINEAR EQUATIONS
C                   (AT*A)*X=AT*Y
C     WHERE (AT*A) IS SYMETRICAL POSITIV DEFINITE MATRIX OF
C     DIMENSION NT*NT, AND (AT*Y) IS A VECTOR OF DIMENSION NT*1.
C
C     CONTEND OF ARRAYES:
C             AN(.)         THE UPPER PART OF (AT*A), AND RETURNS
C                           WITH LT, WHERE L*LT=(AT*A), IF LRED =
C                           .TRUE.
C             INUL(.)       INDEX OF THE FIRST NON-ZERO ELEMENT
C                           OF EACH ROW.
C             H(.)          THE RIGHT-HANDSIDE (AT*Y), AND RETURNS
C                           WITH X ,IF LBS = .TRUE., ELSE WITH
C                           (L-1)*(AT*Y).
C             VAR           THE PSEUDO DIAGONAL ELEMENT OF (L-1)*
C                           (AT*Y).
C
C
C     PROGRAMMED BY
C                         PER KNUDSEN
C                         GEODETIC INSTITUTE
C                         DK-2920 CHARLOTTENLUND           12.07.85.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AN(IANT),INUL(INULT),H(IHT)
      LOGICAL LRED,LBS
C
C***  THE UPPER PART OF A IS REDUCED INTO LT IF LRED IS TRUE.
C
      IF(.NOT.LRED) GO TO 50
      AN(1)=DSQRT(AN(1))
      DO 25 IS=2,NT
      IST=(IS*(IS-1))/2
      SUM=0.0D0
      IRMIN=INUL(IS)
      IRMAX=IS-1
      DO 10 IR=IRMIN,IRMAX
      IRT=(IR*(IR-1))/2
      SUM=0.0D0
      IIMIN=MAX0(INUL(IS),INUL(IR))
      IIMAX=IR-1
      IF(IIMIN.GT.IIMAX) GO TO 6
      DO 5 II=IIMIN,IIMAX
      SUM=SUM+(AN(IRT+II)*AN(IST+II))
    5 CONTINUE
    6 CONTINUE
      AN(IST+IR)=(AN(IST+IR)-SUM)/AN(IRT+IR)
   10 CONTINUE
      SUM=0.0D0
      IIMIN=INUL(IS)
      IIMAX=IS-1
      DO 15 II=IIMIN,IIMAX
      SUM=SUM+AN(IST+II)**2
   15 CONTINUE
      AN(IST+IS)=DSQRT(AN(IST+IS)-SUM)
   25 CONTINUE
C
C
   50 CONTINUE
C***  SOLVE L-1*H
C
      DO 100 IR=1,NT
      IRT=(IR*(IR-1))/2
      SUM=0.0D0
      IIMIN=INUL(IR)
      IIMAX=IR-1
      IF(IIMIN.GT.IIMAX) GO TO 91
      DO 90 II=IIMIN,IIMAX
      SUM=SUM+(AN(IRT+II)*H(II))
   90 CONTINUE
   91 CONTINUE
      H(IR)=(H(IR)-SUM)/AN(IRT+IR)
  100 CONTINUE
      SUM=0.0D0
      DO 101 II=1,NT
      SUM=SUM+H(II)**2
  101 CONTINUE
      VAR=SUM
C
C***  THE SOLUTION IS FOUND BY BACK SUBSTITUTUION IF LBS IS TRUE.
C
      IF(.NOT.LBS) RETURN
C
      DO 150 IRR=1,NT
      IR=NT+1-IRR
      IRT=(IR*(IR-1))/2
      SUM=0.0D0
      IIMIN=IR+1
      IF(IIMIN.GT.NT) GO TO 141
      DO 140 II=IIMIN,NT
      IIT=(II*(II-1))/2
      SUM=SUM+(AN(IIT+IR)*H(II))
  140 CONTINUE
  141 CONTINUE
      H(IR)=(H(IR)-SUM)/AN(IRT+IR)
  150 CONTINUE
C
      RETURN
      END
