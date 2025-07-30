      program sp1d
c $Id: sp1d.for 243 2008-10-29 10:10:19Z cct $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                    
c                                                                               
c                     S P 1 D
c                                                                               
c  program for spherical fft transformation using 
c  Haagmans 1D FFT method. All data are for speed and simplicity 
c  stored in main memory.
c                                                                               
c  input to program:                                                            
c                                                                               
c      <gridfile>                                                               
c      <ofile>                                                                  
c      lmean, dist1, dist2, nmax, ine
c
c  'lmean'  if true the mean of the data is removed prior to FFT and
c           windowing.
c
c  'dist1, dist2' is the inner and outer range of the kernels in degrees
c
c  'nmax'   is the maximum spherical harmonic wong-gore modification
c           of stokes' integral. Only used for mode=1.
c
c  'ine'    is the number of e-w points actually transformed (may be a subgrid).
c           if the implied subgrid is greater than the
c           actual grid zero padding will be done. 
c           ine must be even, and should be a number with a good prime 
c           factorization to ensure a good fft speed.
c
c  Overwiev of files:
c                                                                               
c  file 20   data grid (free format, initiated                        
c            with label lat1, lat2, lon1, lon2, dlat, dlon                      
c            specifying grid boundaries)                                        
c  file 33   output file (note that results close to the                        
c            edges are unreliable due to fft-periodicity)                       
c                                                                               
c  (c) Rene Forsberg, Kort- og Matrikelstyrelsen,
c      Rentemestervej 8, DK-2400 Copenhagen NV, Denmark.
c      First version dec 1999 - based on 'spfour'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)                                        
      dimension nnarr(2)
      logical lgeog,lmean,liz
      character*128 dfile1,ofile                                          
c                                                                               
c  fixed array bounds   ihadim = 160 x 160 = 25600                              
c  ipwdim = 160/2 + 1 = 81                                                      
c  iwkdim = 160*2 = 320                                                         
c                                                                               
      dimension cha(2, 680000)
      dimension wrk(    9600)
      dimension cosfi(3000),gamma(3000)
      dimension dd(2,5000),sum(2,5000)
      ihadim = 680000
      iwkdim = 9600
      maxn = 3000
      maxe = 5000
      idim2 = 2*ihadim                                                          
c                                                                               
c  input data
c  ----------
c
      read(*,1) dfile1                                                          
      read(*,1) ofile                                                           
1     format(a128)                                                               
c                                                                               
      read(*,*) lmean,psi1,psi2,nmax,ine
c
      liz = (psi1.le.0)
c
c  constants                                                                    
c  ---------
c                                                                               
      pi = 3.141592654d0
      re = 6371000.d0
      radeg = 180.d0/pi
c
      gmin = 9.9d9
      gmax = -9.9d9
      rmsiz = 0
      rmaxiz = -9.9d9
c
c  open files - unit 30 scratch file 
c                                                                               
      open(20,file=dfile1,status='old')                                         
      if (liz) open(30,status='scratch',form='unformatted')
      open(33,file=ofile,status='unknown')
c                                                                               
      write(*,2) dfile1,psi1,psi2,ine
2     format(/,
     .' **************************************************************'/
     .' *   SP1D  -  GRAVSOFT 1D-spherical FFT - (c) RF/KMS dec 1999 *'/
     .' **************************************************************'/
     ./' gridfile: ',a36,
     ./' distance range = ',f9.3,' to ',f9.3,' degrees',
     ./' wanted fft transformation e-w points: ',i6)
      write(*,*) 'GEOID HEIGHT FROM GRAVITY BY 1D-FFT'
      if (nmax.ge.2) WRITE(*,*)
     .'- STOKES KERNEL MODIFIED TO NMAX = ',nmax,' -'
      write(*,*)
      write(*,*) 'Input grid:'
c                                                                               
c  read grid data 
c  ---------------
c  convert angles to radians
c
      call rdgrid(20, rfic, rlac, inn, ine, dfi, dla,
     .lgeog, cha, 1, ii1z, ii2z, jj1z, jj2z, ihadim)
      if (.not.lgeog) stop 'only geographic grids allowed'
      if (inn.gt.maxn) stop 'too many rows - increase maxn'
      if (ine.gt.maxe) stop 'too many e-grid pts - increase maxe' 
      if (ii1z.ne.0.or.ii2z.ne.0.or.jj1z.ne.0) stop 'zero pad error'
      n = inn*ine
      rfic = rfic/radeg
      rlac = rlac/radeg
      dfi = dfi/radeg
      dla = dla/radeg
      psi1 = psi1/radeg
      psi2 = psi2/radeg
      spsi1 = sin(psi1/2)
      spsi2 = sin(psi2/2)
c                                                                               
c  data preparation: find mean, and remove if wanted. 
c                                                                               
      r = 0.0                                                                   
      s = 0.0                                                                   
      do 31 i = 1, n                                                            
        r = r + cha(1,i)                                                        
        s = s + cha(1,i)**2                                                     
        cha(2,i) = 0.0                                                          
31    continue                                                                  
      rmean = r/n                                                               
      s = s/n                                                                   
      write(*, 32) s, rmean                                                     
32    format(/' power space domain ',f10.2,', mean ',f9.2)                      
      if (lmean) then
        do 33 i = 1, n
33      cha(1,i) = cha(1,i)-rmean
        write(*,34)
34      format(' mean value subtracted from input data prior to fft')
      else
        write(*,35)
35      format(' mean value not removed from input data')
      endif
c                                                                               
c  store innerzone correction, multiply data by cosfi and data factor
c  -----------------------------------------------------------------
c                                                                               
      if (liz) write(*,*) '- make innerzone corrections -'
      kk = 0
      ii1 = 1
      ii2 = inn
      jj1 = 1
      jj2 = ine-jj2z
      do 50 ii = inn, 1, -1
        rfi = rfic + (ii-1)*dfi
        cosfi(ii) = cos(rfi)
        gamma(ii) = 978032.d0*(1 + 0.00528*sin(rfi))
        if (liz) then
          s = sqrt(dfi*dla*abs(cosfi(ii))/pi)*re
          s = s/gamma(ii)
        endif
        do 50 jj = 1, ine
          kk = kk + 1
          if (liz) then
            if (ii1.le.ii.and.ii.le.ii2.and
     .      .jj1.le.jj.and.jj.le.jj2)
     .      write(30) s*cha(1,kk)
          endif
          cha(1,kk) = cha(1,kk)*cosfi(ii)
50    continue
      if (liz) rewind(30)
c
c  take fourier transform of data*cos(fi) and store rowwise
c  -----------------------------------------------------------
c
      nnarr(1) = ine
      nyqe = ine/2+1
c
      write(*,*) '- fourier transform of data rows -'
      k = 0
      do 51 i = inn,1,-1
        do j = 1, ine
          dd(1,j) = cha(1,k+j)
          dd(2,j) = 0
        enddo
        call fourt(dd,nnarr,1,-1,0,wrk,idim2,iwkdim)
        do j = 1, ine
          cha(1,k+j) = dd(1,j)
          cha(2,k+j) = dd(2,j)
        enddo
        k = k+ine
51    continue
      write(*,*) '- fourier transformation of data completed -'
c
c  main loop: evaluate rowwise from N to S
c  set stokes operator symmetric
c  ----------------------------------------
c
      rfi2 = rfic + (inn-1)*dfi
      rla2 = rlac + (ine-1-jj2z)*dla
      write(33,55) rfic*radeg, rfi2*radeg,
     .rlac*radeg, rla2*radeg, dfi*radeg, dla*radeg
55    format(' ',4f12.6,2f12.7)
c
      do 80 ip = inn,1,-1
c
        do j = 1,ine
          sum(1,j) = 0
          sum(2,j) = 0
        enddo
c 
        k = 0
        do 60 i = inn,1,-1
          if (abs((ip-i)*dfi).gt.psi2) goto 59
          do j = 1,nyqe
            s = sqrt((sin((ip-i)*dfi/2))**2 + 
     .      (sin((j-1)*dla/2))**2*cosfi(ip)*cosfi(i))		 
            if (s.gt.spsi1.and.s.le.spsi2) then
              dd(1,j) = stokes(s,nmax)
            else
              dd(1,j) = 0
            endif
            dd(2,j) = 0
          enddo
          do j = nyqe+1,ine
            dd(1,j) = dd(1,ine-j+2)
            dd(2,j) = 0
          enddo
          call fourt(dd,nnarr,1,-1,0,wrk,idim2,iwkdim)
c
          do j = 1,ine
            sum(1,j) = sum(1,j) + cha(1,k+j)*dd(1,j)
            sum(2,j) = sum(2,j) + cha(2,k+j)*dd(1,j)
          enddo
59        k = k + ine
60      continue
c
        call fourt(sum,nnarr,1,1,1,wrk,idim2,iwkdim)
c
c  write out result row
c
        fak = re/(4*pi*gamma(ip))/ine*dfi*dla
        do j = 1,ine-jj2z
          sum(1,j) = fak*sum(1,j)
          if (liz) then
            read(30) s
            rmsiz = rmsiz+s**2
            if (abs(s).gt.rmaxiz) rmaxiz = abs(s)
            sum(1,j) = sum(1,j)+s
          endif
          if (sum(1,j).lt.gmin) gmin = sum(1,j)
          if (sum(1,j).gt.gmax) gmax = sum(1,j)
        enddo
        write(33,75) (sum(1,j), j=1,ine-jj2z)
75      format(/,60(/,' ',8f9.3))
        if (mod(ip,20).eq.0) write(*,*) '... completed row ',ip
80    continue
c
c  exit
c  ----
c
      write(*,110) rfic*radeg,rfi2*radeg,rlac*radeg,rla2*radeg,
     .dfi*radeg,dla*radeg,gmin,gmax
      if (liz) write(*,111) sqrt(rmsiz/(inn-ii1z-ii2z)
     ./(ine-jj1z-jj2z)), rmaxiz
110   format(/' Output grid written:'/,
     .' ',4f11.5,2f10.6,/
     .' Minimal and maximal values: ',2f10.3)
111   format(
     .' R.m.s. and max abs value of innerzone correction: ',2f9.3)
c
      if (liz) close(30,status='delete')
      end
c
      subroutine fourt(datt,nn,ndim,isign,iform,work,
     .idim1,idim2)
      implicit double precision(a-h,o-z)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                     f o u r t
c
c        version=740301
c        program description norsar n-pd9 dated 1 july 1970
c        author n m brenner
c        further description    three fortran programs etc.
c        issued by lincoln laboratory, mit, july 1967
c        two corrections by hjortenberg 1974
c     the fast fourier transform in usasi basic fortran
c
c     modified to rc fortran rf june 84
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension datt(idim1),nn(ndim),ifact(32),work(idim2)
c
      np0=0
      nprev=0
c
      twopi=6.283185307d0
      rthlf=.7071067812d0
      if(ndim-1)920,1,1
1     ntot=2
      do 2 idim=1,ndim
      if(nn(idim))920,920,2
2     ntot=ntot*nn(idim)
c
c     mainloop for each dimension
c
      np1=2
      do 910 idim=1,ndim
      n=nn(idim)
      np2=np1*n
      if(n-1)920,900,5
c
c     is n a power of two and if not, what are its factors
c
5     m=n
      ntwo=np1
      if=1
      idiv=2
10    iquot=m/idiv
      irem=m-idiv*iquot
      if(iquot-idiv)50,11,11
11    if(irem)20,12,20
12    ntwo=ntwo+ntwo
      ifact(if)=idiv
      if=if+1
      m=iquot
      go to 10
20    idiv=3
      inon2=if
30    iquot=m/idiv
      irem=m-idiv*iquot
      if(iquot-idiv)60,31,31
31    if(irem)40,32,40
32    ifact(if)=idiv
      if=if+1
      m=iquot
      go to 30
40    idiv=idiv+2
      go to 30
50    inon2=if
      if(irem)60,51,60
51    ntwo=ntwo+ntwo
      go to 70
60    ifact(if)=m
70    non2p=np2/ntwo
c
c     separate four cases---
c        1. complex transform
c        2. real transform for the 2nd, 3nd, etc. dimension.  method--
c           transform half the datt, supplying the other half by con-
c           jugate symmetry.
c        3. real transform for the 1st dimension,n odd.  method--
c           set the imaginary parts to zero
c        4. real transform for the 1st dimension,n even.method--
c           transform a complex array of lenght n/2 whose real parts
c           are the even numberd real values and whose imaginary parts
c           are the odd numberedreal values.  separate and supply
c           the second half by conjugate summetry.
c
      icase=1
      ifmin=1
      i1rng=np1
      if(idim-4)74,100,100
74    if(iform)71,71,100
71    icase=2
      i1rng=np0*(1+nprev/2)
      if(idim-1)72,72,100
72    icase=3
      i1rng=np1
      if(ntwo-np1)100,100,73
73    icase=4
      ifmin=2
      ntwo=ntwo/2
      n=n/2
      np2=np2/2
      ntot=ntot/2
      i=1
      do 80 j=1,ntot
      datt(j)=datt(i)
80    i=i+2
c
c     shuffle datt by bit reversal, since n=2**k.  as the shuffling
c     can be done by simple interchange, no working array is needed
c
100   if(non2p-1)101,101,200
101   np2hf=np2/2
      j=1
      do 150 i2=1,np2,np1
      if(j-i2)121,130,130
121   i1max=i2+np1-2
      do 125 i1=i2,i1max,2
      do 125 i3=i1,ntot,np2
      j3=j+i3-i2
      tempr=datt(i3)
      tempi=datt(i3+1)
      datt(i3)=datt(j3)
      datt(i3+1)=datt(j3+1)
      datt(j3)=tempr
125   datt(j3+1)=tempi
130   m=np2hf
140   if(j-m)150,150,141
141   j=j-m
      m=m/2
      if(m-np1)150,140,140
150   j=j+m
      go to 300
c
c     shuffle datt by digit reversal for general n
c
200   nwork=2*n
      do 270 i1=1,np1,2
      do 270 i3=i1,ntot,np2
      j=i3
      do 260 i=1,nwork,2
      if(icase-3)210,220,210
210   work(i)=datt(j)
      work(i+1)=datt(j+1)
      go to 240
220   work(i)=datt(j)
      work(i+1)=0.
240   ifp2=np2
      if=ifmin
250   ifp1=ifp2/ifact(if)
      j=j+ifp1
      if(j-i3-ifp2)260,255,255
255   j=j-ifp2
      ifp2=ifp1
      if=if+1
      if(ifp2-np1)260,260,250
260   continue
      i2max=i3+np2-np1
      i=1
      do 270 i2=i3,i2max,np1
      datt(i2)=work(i)
      datt(i2+1)=work(i+1)
270   i=i+2
c
c     main loop for factors of two
c     w=exp(isign*2*pi*sqrt(-1)*m/(4*mmax)).  check for w=isign*sqrt(-1)
c     and repeat for w=w*(1+isign*sqrt(-1))/sqrt(2)
c
300   if(ntwo-np1)600,600,305
305   np1tw=np1+np1
      ipar=ntwo/np1
310   if(ipar-2)350,330,320
320   ipar=ipar/4
      go to 310
330   do 340 i1=1,i1rng,2
      do 340 k1=i1,ntot,np1tw
      k2=k1+np1
      tempr=datt(k2)
      tempi=datt(k2+1)
      datt(k2)=datt(k1)-tempr
      datt(k2+1)=datt(k1+1)-tempi
      datt(k1)=datt(k1)+tempr
340   datt(k1+1)=datt(k1+1)+tempi
350   mmax=np1
360   if(mmax-ntwo/2)370,600,600
370   lmax=max0(np1tw,mmax/2)
      do 570 l=np1,lmax,np1tw
      m=l
      if(mmax-np1)420,420,380
380   theta=-twopi*dble(l)/dble(4*mmax)
      if(isign)400,390,390
390   theta=-theta
400   wr=cos(theta)
      wi=sin(theta)
410   w2r=wr*wr-wi*wi
      w2i=2.*wr*wi
      w3r=w2r*wr-w2i*wi
      w3i=w2r*wi+w2i*wr
420   do 530 i1=1,i1rng,2
      kmin=i1+ipar*m
      if(mmax-np1)430,430,440
430   kmin=i1
440   kdif=ipar*mmax
450   kstep=4*kdif
      if(kstep-ntwo)460,460,530
460   do 520 k1=kmin,ntot,kstep
      k2=k1+kdif
      k3=k2+kdif
      k4=k3+kdif
      if(mmax-np1)470,470,480
470   u1r=datt(k1)+datt(k2)
      u1i=datt(k1+1)+datt(k2+1)
      u2r=datt(k3)+datt(k4)
      u2i=datt(k3+1)+datt(k4+1)
      u3r=datt(k1)-datt(k2)
      u3i=datt(k1+1)-datt(k2+1)
      if(isign)471,472,472
471   u4r=datt(k3+1)-datt(k4+1)
      u4i=datt(k4)-datt(k3)
      go to 510
472   u4r=datt(k4+1)-datt(k3+1)
      u4i=datt(k3)-datt(k4)
      go to 510
480   t2r=w2r*datt(k2)-w2i*datt(k2+1)
      t2i=w2r*datt(k2+1)+w2i*datt(k2)
      t3r=wr*datt(k3)-wi*datt(k3+1)
      t3i=wr*datt(k3+1)+wi*datt(k3)
      t4r=w3r*datt(k4)-w3i*datt(k4+1)
      t4i=w3r*datt(k4+1)+w3i*datt(k4)
      u1r=datt(k1)+t2r
      u1i=datt(k1+1)+t2i
      u2r=t3r+t4r
      u2i=t3i+t4i
      u3r=datt(k1)-t2r
      u3i=datt(k1+1)-t2i
      if(isign)490,500,500
490   u4r=t3i-t4i
      u4i=t4r-t3r
      go to 510
500   u4r=t4i-t3i
      u4i=t3r-t4r
510   datt(k1)=u1r+u2r
      datt(k1+1)=u1i+u2i
      datt(k2)=u3r+u4r
      datt(k2+1)=u3i+u4i
      datt(k3)=u1r-u2r
      datt(k3+1)=u1i-u2i
      datt(k4)=u3r-u4r
520   datt(k4+1)=u3i-u4i
      kdif=kstep
      kmin=4*(kmin-i1)+i1
      go to 450
530   continue
      m=m+lmax
      if(m-mmax)540,540,570
540   if(isign)550,560,560
550   tempr=wr
      wr=(wr+wi)*rthlf
      wi=(wi-tempr)*rthlf
      go to 410
560   tempr=wr
      wr=(wr-wi)*rthlf
      wi=(tempr+wi)*rthlf
      go to 410
570   continue
      ipar=3-ipar
      mmax=mmax+mmax
      go to 360
c
c     main loop for factoers not equal to two
c     w=exp(isign*2*pi*sqrt(-1)*(j1+j2-i3-1)/ifp2)
c
600   if(non2p-1)700,700,601
601   ifp1=ntwo
      if=inon2
610   ifp2=ifact(if)*ifp1
      theta=-twopi/dble(ifact(if))
      if(isign)612,611,611
611   theta=-theta
612   wstpr=cos(theta)
      wstpi=sin(theta)
      do 650 j1=1,ifp1,np1
      thetm=-twopi*dble(j1-1)/dble(ifp2)
      if(isign)614,613,613
613   thetm=-thetm
614   wminr=cos(thetm)
      wmini=sin(thetm)
      i1max=j1+i1rng-2
      do 650 i1=j1,i1max,2
      do 650 i3=i1,ntot,np2
      i=1
      wr=wminr
      wi=wmini
      j2max=i3+ifp2-ifp1
      do 640 j2=i3,j2max,ifp1
      twowr=wr+wr
      j3max=j2+np2-ifp2
      do 630 j3=j2,j3max,ifp2
      jmin=j3-j2+i3
      j=jmin+ifp2-ifp1
      sr=datt(j)
      si=datt(j+1)
      oldsr=0.
      oldsi=0.
      j=j-ifp1
620   stmpr=sr
      stmpi=si
      sr=twowr*sr-oldsr+datt(j)
      si=twowr*si-oldsi+datt(j+1)
      oldsr=stmpr
      oldsi=stmpi
      j=j-ifp1
      if(j-jmin)621,621,620
621   work(i)=wr*sr-wi*si-oldsr+datt(j)
      work(i+1)=wi*sr+wr*si-oldsi+datt(j+1)
630   i=i+2
      wtemp=wr*wstpi
      wr=wr*wstpr-wi*wstpi
640   wi=wi*wstpr+wtemp
      i=1
      do 650 j2=i3,j2max,ifp1
      j3max=j2+np2-ifp2
      do 650 j3=j2,j3max,ifp2
      datt(j3)=work(i)
      datt(j3+1)=work(i+1)
650   i=i+2
      if=if+1
      ifp1=ifp2
      if(ifp1-np2)610,700,700
c
c     complete areal transform in the 1st dimension, n even, by con-
c     jugate symmetries
c
700   go to (900,800,900,701),icase
701   nhalf=n
      n=n+n
      theta=-twopi/dble(n)
      if(isign)703,702,702
702   theta=-theta
703   wstpr=cos(theta)
      wstpi=sin(theta)
      wr=wstpr
      wi=wstpi
      imin=3
      jmin=2*nhalf-1
      go to 725
710   j=jmin
      do 720 i=imin,ntot,np2
      sumr=(datt(i)+datt(j))/2.
      sumi=(datt(i+1)+datt(j+1))/2.
      difr=(datt(i)-datt(j))/2.
      difi=(datt(i+1)-datt(j+1))/2.
      tempr=wr*sumi+wi*difr
      tempi=wi*sumi-wr*difr
      datt(i)=sumr+tempr
      datt(i+1)=difi+tempi
      datt(j)=sumr-tempr
      datt(j+1)=-difi+tempi
720   j=j+np2
      imin=imin+2
      jmin=jmin-2
      wtemp=wr*wstpi
      wr=wr*wstpr-wi*wstpi
      wi=wi*wstpr+wtemp
725   if(imin-jmin)710,730,740
730   if(isign)731,740,740
731   do 735 i=imin,ntot,np2
735   datt(i+1)=-datt(i+1)
740   np2=np2+np2
      ntot=ntot+ntot
      j=ntot+1
      imax=ntot/2+1
745   imin=imax-2*nhalf
      i=imin
      go to 755
750   datt(j)=datt(i)
      datt(j+1)=-datt(i+1)
755   i=i+2
      j=j-2
      if(i-imax)750,760,760
760   datt(j)=datt(imin)-datt(imin+1)
      datt(j+1)=0.
      if(i-j)770,780,780
765   datt(j)=datt(i)
      datt(j+1)=datt(i+1)
770   i=i-2
      j=j-2
      if(i-imin)775,775,765
775   datt(j)=datt(imin)+datt(imin+1)
      datt(j+1)=0.
      imax=imin
      go to 745
780   datt(1)=datt(1)+datt(2)
      datt(2)=0.
      go to 900
c
c     complete a real transform for the 2nd, 3rd, etc. dimension by
c     conjugate symmetries.
c
800   if(i1rng-np1)805,900,900
805   do 860 i3=1,ntot,np2
      i2max=i3+np2-np1
      do 860 i2=i3,i2max,np1
      imax=i2+np1-2
      imin=i2+i1rng
      jmax=2*i3+np1-imin
      if(i2-i3)820,820,810
810   jmax=jmax+np2
820   if(idim-2)850,850,830
830   j=jmax+np0
      do 840 i=imin,imax,2
      datt(i)=datt(j)
      datt(i+1)=-datt(j+1)
840   j=j-2
850   j=jmax
      do 860 i=imin,imax,np0
      datt(i)=datt(j)
      datt(i+1)=-datt(j+1)
860   j=j-np0
c
c     end of loop on each dimension
c
900   np0=np1
      np1=np2
910   nprev=n
920   return
      end
c
      subroutine rdgrid(iunit, rfic, rlac, inn, ine, dfi, dla,
     .lgeo, cha, ik, ii1z, ii2z, jj1z, jj2z, idim)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                       r d g r i d
c
c  subroutine for reading a digital grid file on
c  standard format, i.e. stored rowwise from nw to se, with label.
c
c  as se-corner coordinate 'rfic, rlac' (degrees) will be used
c  (unless they are zero, then the grid corner is used).
c  a grid containing 'inn' x 'ine' points will be put in array
c  'cha' of declared dimension 'idim'.
c  if inn=0 the complete lat range will be read.
c  if ine=0 the complete lon range will be read
c  if the wanted grid is too large a zero padding will be done.
c
c  last updated jun 90, rf, updated dec 99
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      common /gridpar/ iell,izone
      dimension cha(2, idim)
      dimension hlab(6), hrow(5000)
      logical lgeo, lutm
      irdim = 5000
c
c  initialize statistics
c
      nr = 0
      rsum = 0.0
      rsum2 = 0.0
      rmin = 99999
      rmax = -99999
c
      read(iunit,*) (hlab(j),j=1,6)
      lutm = (abs(hlab(1)).ge.400.or.abs(hlab(2)).ge.400)
      lgeo = (.not.lutm)
      if (lgeo) goto 111
      read(iunit,*) iell,izone
      write(*,110) iell,izone
      if (iell.lt.1.or.iell.gt.3.or.izone.lt.1.or.izone.gt.60)
     *stop 'illegal ellipsoid or utm zone'
110   format(' - input grid in utm, ell ',i1,' zone ',i2,' -')
111   dfi = hlab(5)
      dla = hlab(6)
      nn = (hlab(2)-hlab(1))/dfi+1.5
      ne = (hlab(4)-hlab(3))/dla+1.5
      if (nn.gt.irdim.or.ne.gt.irdim) stop 'too long rows/columns'
c
c  find corner indices for wanted subgrid
c
      if (inn.eq.0) inn = nn
      if (ine.eq.0) ine = ne
      if (rfic.eq.0.and.rlac.eq.0) then
        rfic = hlab(1)
        rlac = hlab(3)
      endif
      ifi1 = (rfic-hlab(1))/dfi+1.5
      ila1 = (rlac-hlab(3))/dla+1.5
      ifi2 = ifi1+inn-1
      ila2 = ila1+ine-1
      rfic = (ifi1-1)*dfi + hlab(1)
      rlac = (ila1-1)*dla + hlab(3)
      n = inn*ine
c
c  check boundaries for padding 
c
      ii1z = 0
      ii2z = 0
      jj1z = 0
      jj2z = 0
      if (ifi1.lt.1) ii1z = 1-ifi1
      if (ifi2.gt.nn) ii2z = ifi2-nn
      if (ila1.lt.1) jj1z = 1-ila1
      if (ila2.gt.ne) jj2z = ila2-ne
c
      if (n.gt.idim) then
        write(*, 122) n,idim
122     format(' *** array dim too small - wanted, declared ',2i8)
        stop ' *** sorry ***'
      endif
c
c  read data grid values
c  data in cha array stored with first element at nw corner
c
      ilast = ifi1
      if (ilast.lt.1) ilast = 1
c
      do 130 i = nn,ilast,-1
c
        read(iunit,*,end=131) (hrow(j),j=1,ne)
c
        if (i.lt.ifi1.or.i.gt.ifi2) goto 130
        jj0 = (ifi2-i)*ine - ila1+1
        do 129 j = 1,ne
          r = hrow(j)
          if (j.lt.ila1.or.j.gt.ila2) goto 129
          cha(ik,j+jj0) = r
          nr = nr + 1
          if (r.gt.rmax) rmax = r
          if (r.lt.rmin) rmin = r
          rsum = rsum + r
          rsum2 = rsum2 + r**2
129     continue
130   continue
      goto 133
131     write(*,132) i
132     format(' *** too few data in grid file, lastrow = ',i7)
        stop ' *** check grid label and data ***'
c
c  zero padding
c
133   if (ii1z+ii2z+jj1z+jj2z.gt.0) then
        do 138 i = inn,1,-1
          jj0 = (inn-i)*ine
          if (i.gt.ifi2.or.i.lt.ifi1) then
            do 134 j = 1, ine
134         cha(ik,j+jj0) = 0 
          else
            do 135 j = 1, jj1z
135         cha(ik,j+jj0) = 0
            do 136 j = ine-jj2z+1, ine
136         cha(ik,j+jj0) = 0
          endif
138     continue
      endif
c
c  write information and statistics
c
      if (nr.eq.0) stop '*** no points read from grid, wrong area'
      rfi = hlab(1) + (ifi1-1)*dfi
      rla = hlab(3) + (ila1-1)*dla
      r = rsum/nr
      s = 0.0
      if (n.gt.1)
     .s = sqrt((rsum2 - rsum**2/nr)/(nr-1))
      if (lgeo) write(*,141) (hlab(j),j=1,6),nn,ne
      if (lutm) write(*,142) (hlab(j),j=1,6),nn,ne
141   format(' Gridlab:',4f10.4,2f9.4,i5,i4)
142   format(' Gridlab:',4f10.0,2f8.0,i5,i4)
      if (lgeo) write(*, 143) rfi, rla, inn, ine
143   format(' Selected: sw corner ',2f10.4, ', points ', 2i6,i8)
      if (lutm) write(*, 144) rfi, rla, inn, ine
144   format(' Selected: sw corner ',2f10.0, ', points ', 2i6,i8)
      write(*, 145) nr, r, s, rmin, rmax
145   format(' Statistics of data selected from input grid:'
     ./' pts mean std.dev. min max:',i8,4f9.2)
      if (ii1z+ii2z+jj1z+jj2z.gt.0) write(*,146) ii1z,ii2z,jj1z,jj2z
146   format(' Zero padding done on grid, no of rows/cols S/N/E/W:',4i4)
      return
      end
c
      real*8 function stokes(s,nmax)
      implicit real*8(a-h,o-z)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  evaluates stokes function with wong-gore modification of low
c  harmonics using recursion algorithms for pn
c  s is sin(psi/2)
c
c  rf aug 97
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      rr = 1.d0/s-4-6*s+10*s**2-(3-6*s**2)*log(s+s**2)
      if (nmax.le.1) then
        stokes = rr
      else
        t = cos(asin(s)*2)
        pn2 = 1.d0
        pn1 = t
        sum = 0.d0
        do 10 n = 2,nmax
          pn = ((2*n-1)*t*pn1 - (n-1)*pn2)/n
          sum = pn*(2*n+1)/(n-1) + sum
          pn2 = pn1
          pn1 = pn
10      continue
        stokes = rr-sum
      endif
      return
      end
