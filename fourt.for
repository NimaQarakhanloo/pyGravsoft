      stop
      end
c
      subroutine fourt(datt,nn,ndim,isign,iform,work,idim1,idim2)
c $Id: fourt.for 181 2008-08-14 23:19:51Z tjansson $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      implicit double precision(a-h,o-z)
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
