      program spfour 
c $Id: spfour.for 181 2008-08-14 23:19:51Z tjansson $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
c                     S P F O U R
c                   
c  program for spherical fft transformation using a sequence
c  of reference latitude values, with linear interpolation of
c  spherical formulas between latitude parallels.
c
c  the method is based on an expansion of the spherical fft
c  method of strang van hess (manuscripta geodaetica, 1990).
c  the program uses analytical transforms of the fft kernel,
c  this allows additionally control over the effective integration
c  radius.
c                  
c  a subgrid of the grid in 'gridfile' is analyzed. the subgrid
c  is specified by its sw corner (rfic, rlac) and number of points
c  (inn, ine). The numbers should have a good prime factorization
c  with many small factors for the FFT to obtain maximum speed.
c
c  input to program:
c                  
c      <gridfile> 
c      <ofile>   
c      mode, lmean, nref, dist1, dist2
c      (nmod1, nmod2 - only mode 1) 
c      (height - only mode 3 or 11)
c      (<reffile> - only mode 13)
c      fic, lac, inn, ine, iwn
c               
c  'mode' determines the function of the program:
c              
c        1    geoid prediction from gravity data
c        2    (deflections - not implemented yet)
c        3    harmonic upward continuation of gravity data 
c       10    gravity effect of isostasy - 2nd order expansion 
c             (sign: positive on land. add to Bouguer anomalies)
c       11    bouguer effect in airborne gravimetry at "height"
c       12    geoid terrain effect (mass layer approximation)
c       13    geoid RTM terrain effect (3rd order approximation)
c             reffile must be exactly same grid parameters as the basic grid
c
c  'lmean'  if true the mean of the data is removed prior to FFT and
c           windowing.
c
c  'nref'   is the number of reference parallels used. if nref = 1 just
c           one fft is done, nref = 2 yields two ffts with reference on
c           north and south boundary of subgrid, npar = 3 yields three
c           FFT's with reference latitudes on north, central and south
c           parallel, npar = 4 ...
c
c  'dist1, dist2' is the inner and outer range of the kernels in degrees
c
c  'nmod1, nmod2'  is the maximum spherical harmonic modification
c           of stokes' integral. Low harmonics are removed from stokes' 
c           function up to degree 'nmod1' ("wong-gore"), and then linearly
c           tapered to 'nmod2'. Only used for mode=1.
c
c  'height' is the height in m
c
c  'fic, lac' is geographical coordinates of the wanted SW-corner
c
c  'inn, ine' is the number of points actually transformed (may be a subgrid).
c           if the implied subgrid by fic,lac,inn,ine is greater than the
c           actual grid zero padding will be done. The reference latitude
c           will not be affected by the zero padding, and the padded values
c           will not be output (the output grid might thus be smaller than
c           the wanted number of points)
c
c  'iwn'    if iwn.gt.0 the 'iwn' points closed to the edge of the internal
c           selected grid will be windowed by a cosine-tpaer window.
c           note: the windowing is done after taking the data mean!
c
c  Overwiev of files:
c                                                                               
c  file 20   data grid (free format, initiated                        
c            with label lat1, lat2, lon1, lon2, dlat, dlon                      
c            specifying grid boundaries)                                        
c  file 33   output file (note that results close to the                        
c            edges are unreliable due to fft-periodicity)                       
c  file 30,31,32  scratch files, used for intermediate results
c                                                                               
c                                                                               
c  (c) Rene Forsberg, Kort og Matrikelstyrelsen (KMS),
c      Rentemestervej 8, DK-2400 Copenhagen NV, Denmark.
c
c  programmed at u of c, september 1990, based on 'tcfour' and 'geofour'
c  last update feb 1992, july 1992, nov 1992, may 1994
c  wong-gore aug 97
c  mode 13 and updates feb 2002
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)        
      dimension nnarr(2),rfiref(100),irow1(100),irow2(100),wf(50)
      logical lgeog, lgeogr, lmean, liz, lref, lfirst, lzero
      character*72 dfile1,ofile,reffile
c                                                
c  fixed array bounds   ihadim = 160 x 160 = 25600
c  ipwdim = 160/2 + 1 = 81
c  iwkdim = 160*2 = 320
c                    
      dimension cha(2, 2000000)
      dimension wrk(    8000)
      dimension sindl2(4000),sindfi(4000),sindf2(4000)
      ihadim = 2000000
      iwkdim = 8000
      idim2 = 2*ihadim
c                    
c  input data
c  ----------
c
      read(*,1) dfile1
      read(*,1) ofile
1     format(a72)
c
      read(*,*) mode,lmean,nref,psi1,psi2
c
      if (mode.eq.1) read(*,*) nmod1,nmod2
      if (mode.eq.3.or.mode.eq.11) read(*,*) height
      lref = (mode.eq.13)
      if (lref) read(*,1) reffile
c
      read(*,*) rfic,rlac,inn,ine,iwndow
c
      if (nref.lt.1) nref = 1
      if (nref.gt.100) stop 'nref is too big'
      liz = (psi1.le.0.and.mode.ne.10)
c
c  constants                                   
c  ---------
c                                             
      pi = 3.141592654d0
      re = 6371000.d0
      radeg = 180.d0/pi
      gamma = 981500.d0
c
      diso = 32000
      rho = 2.67
      grho = 0.00667d0*rho
      drho = 0.4
c
      gmin = 9.9d9
      gmax = -9.9d9
      rmsiz = 0
      rmaxiz = -9.9d9
c
c  open files - unit 30 and 31 are scratch files for data and results
c                                            
      open(20,file=dfile1,status='old')     
      if (liz) open(30,status='scratch',form='unformatted')
      open(31,status='scratch',form='unformatted')
      open(32,status='scratch',form='unformatted')
      open(33,file=ofile,status='unknown')
      if (lref) then
        open(35,status='scratch',form='unformatted')
        open(36,status='scratch',form='unformatted')
      endif
c                                          
      write(*,2) dfile1,mode,psi1,psi2,rfic,rlac,inn,ine
2     format(/,
     .' **************************************************************'/
     .' *   SPFOUR  -  GRAVSOFT spherical FFT - (c) RF/KMS feb 2002  *'/
     .' **************************************************************'/
     ./' Gridfile: ',a72,/
     .' mode = ',i3,', distance range = ',f9.3,' to ',f9.3,' degrees',
     ./' SW corner: ',2f11.4,', points: ',2i5)
      if (mode.eq.1) then
        write(*,*) 'GEOID HEIGHTS FROM GRAVITY'
        if (nmod2.lt.nmod1) nmod2=nmod1
        if (nmod1.ge.2) write(*,3) nmod1,nmod2 
3       format(' stokes kernel harmonics removed 100% to degree',i5,
     .  /' with linear taper to degree',i5)
      endif
      if (mode.eq.3) write(*,*) 'UPWARD CONTINUATION to height ',
     .int(height)
      if (mode.eq.10) write(*,*) 'ISOSTATIC GRAVITY EFFECT'
      if (mode.eq.11) write(*,*) 'BOUGUER EFFECT IN AIRBORNE GRAVITY',
     .' at height',int(height)
      if (mode.eq.12) write(*,*) 'GEOID TERRAIN EFFECT (condensed)'
      if (mode.eq.13) write(*,*) 'GEOID TERRAIN EFFECT (2nd order RTM)'
      write(*,*)
      write(*,*) 'Input grid:'
c
c  set icomplex: 0 = real data, 1 = complex data
c
      icomplex = 0
      if (mode.ge.10) icomplex = 1
c     
c  read grid data 
c  ---------------
c
      rficr = rfic
      rlacr = rlac
      innr = inn
      iner = ine
c
      call rdgrid(20, rfic, rlac, inn, ine, dfi, dla,
     .lgeog, cha, 1, ii1z, ii2z, jj1z, jj2z, ihadim)
c
      if (.not.lgeog) stop 'only geographic grids allowed'
      n = inn*ine
c
c  read ref grid for mode 13
c
      if (mode.eq.13) then
        open(21,file=reffile,status='unknown')
        write(*,*) 'Reference grid: ',reffile
c
        call rdgrid(21, rficr, rlacr, innr, iner, dfir, dlar,
     .  lgeogr, cha, 2, ii1z, ii2z, jj1z, jj2z, ihadim)
c
        if (.not.lgeogr) stop '*** only geographic grids allowed'
        if (abs(rficr-rfic).gt.0.001.or.abs(rlacr-rlac).gt.0.001.
     .  or.abs(dfir-dfi).gt.0.0001.or.abs(dlar-dla).gt.0.0001.
     .  or.inn.ne.innr.or.ine.ne.iner) stop '*** ref grid not the same'
      endif
c
c  convert angles to radians
c
      rfic = rfic/radeg
      rlac = rlac/radeg
      dfi = dfi/radeg
      dla = dla/radeg
      psi1 = psi1/radeg
      psi2 = psi2/radeg
      spsi1 = sin(psi1/2)
      spsi2 = sin(psi2/2)
c
c  data preparation: find mean, and remove if wanted. windowing.
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
      write(*,32) s, rmean
32    format(/' power space domain ',f10.2,', mean ',f9.2)
      if (lmean) then
        do 33 i = 1, n
33      cha(1,i) = cha(1,i)-rmean
        if (lref) cha(2,i) = cha(2,i)-rmean
        write(*,34)
34      format(' mean value subtracted from input data prior to fft')
      else
        write(*,35)
35      format(' mean value not removed from input data')
      endif
c                                                                               
c  windowing of data grid                                                       
c                                                                               
      if (iwndow.gt.0) then
        if (iwndow.gt.50) stop '*** iwndow too big'
        do 38 i = 1, iwndow  
38      wf(i) = (1 - cos(pi*i/(iwndow+1)))/2
        wsum = 0
        k = 0
        sw = 0.0
        do 40 i = inn, 1, -1
        do 40 j = 1, ine
          k = k+1
          w = 1.0
          if (i.gt.iwndow) goto 41
          w = w*wf(i)
41        if (i.le.inn-iwndow) goto 42
          w = w*wf(inn-i+1)
42        if (j.gt.iwndow) goto 43
          w = w*wf(j)
43        if (j.le.ine-iwndow) goto 44
          w = w*wf(ine-j+1)
44         wsum = wsum + w
          cha(1,k) = cha(1,k)*w
          if (lref) cha(2,k) = cha(2,k)*w
          sw = sw + cha(1,k)**2
40      continue
        wsum = wsum/n
        sw = sw/n
        write(*,45) sw, wsum
45      format(' power after window ',f10.2,', wsum = ', f8.4)
      endif
c
c  data factors for terrain effects - negative heights assumed ocean
c  condense depths to equivalent rock
c  -----------------------------------------------------------------
c
      if (mode.ge.10) then
        kk = 0
        do 46 ii = inn,1,-1
        do 46 jj = 1,ine
          kk = kk+1
c
          h = cha(1,kk)
          if (h.lt.0.and.mode.ne.12) h = h*1.64/rho
          if (mode.eq.10) h = rho/drho*h
          cha(1,kk) = h
46      continue
      endif
c                                                                               
c  store innerzone corrections 
c  ---------------------------
c                                                                               
      if (liz) then
        write(*,*) '- make innerzone corrections -'
        kk = 0
        ii1 = 1+ii1z
        ii2 = inn-ii2z
        jj1 = 1+jj1z
        jj2 = ine-jj2z
        do 47 ii = inn, 1, -1
          rfi = rfic + (ii-1)*dfi
          cosfi = abs(cos(rfi))
c
c  equivalent radius in m for grid cell - make iz in final units
c
          r = sqrt(dfi*dla*cosfi/pi)*re
c
          do 47 jj = 1, ine
            kk = kk + 1
            s = cha(1,kk)
            if (ii1.le.ii.and.ii.le.ii2.
     .      and.jj1.le.jj.and.jj.le.jj2) then
              if (mode.eq.1) then
                write(30) r*s/gamma
              elseif (mode.eq.3) then	
                write(30) s*(1-height/sqrt(r**2+height**2))
              elseif (mode.eq.11) then
                write(30) 2*pi*grho
     .          *(s+sqrt(r**2+(height-s)**2)-sqrt(r**2+height**2))
              elseif (mode.eq.12) then 
                write(30) 2*pi*grho*s*r/gamma
              elseif (mode.eq.13) then
                b = cha(1,kk)-cha(2,kk)
                s = sqrt(r**2+b**2)
                write(30) 
     .          pi*grho*(-b**2 + b*s + r**2*log((b+s)/r))/gamma
              endif 
            endif
47      continue
        rewind(30)
      endif
c
c  multiply data by cosfi - store for later use in mode 13
c  -------------------------------------------------------
c
      lfirst = .true.
c
c  entry point for higher order terms in mode 13
c  lfirst = true in first run, false in second
c
471   kk = 0
      do 48 ii = inn,1,-1
        rfi = rfic + (ii-1)*dfi
        cosfi = abs(cos(rfi))
        do 48 jj = 1,ine
          kk = kk+1  
c
          if (lref) then
            h = cha(1,kk)
            href = cha(2,kk)
            if (lfirst) then
              write(35) h,href
              cha(1,kk) = h-href
              cha(2,kk) = h**3-href**3
            else
              read(35) h,href
              cha(1,kk) = h**2-href**2
              cha(2,kk) = h-href
            endif
          elseif (mode.ge.10) then
            h = cha(1,kk)
            cha(2,kk) = h**2
          endif
c
          cha(1,kk) = cha(1,kk)*cosfi
          if (icomplex.eq.1) cha(2,kk) = cha(2,kk)*cosfi
48    continue
c
      if (icomplex.eq.1)
     .write(*,*) '- complex data used for transform -' 
c
c  take fourier transform of data*cos(fi) and store on unit 31
c  -----------------------------------------------------------
c
      nnarr(1) = ine
      nnarr(2) = inn
      nyqn = inn/2+1
      nyqe = ine/2+1
c
      write(*,*) '- fourier transform of data -'
      call fourt(cha,nnarr,2,-1,icomplex,wrk,idim2,iwkdim)
      write(*,*) '- fourier transformation of data completed -'
c
      do 52 i = 1, n
52    write(31) cha(1,i),cha(2,i)
c                                                                               
c  set tables of the reference latitudes
c  zero padding is not taken into account
c  --------------------------------------
c
      rfi1 = rfic + ii1z*dfi
      rla1 = rlac + jj1z*dla
      rfi2 = rfic + (inn-1-ii2z)*dfi
      rla2 = rlac + (ine-1-jj2z)*dla
      if (nref.eq.1) then
        rfiref(1) = (rfi1+rfi2)/2
      else
        do 53 i = 1, nref
53      rfiref(i) = rfi2 - (rfi2-rfi1)*(i-1)/(nref-1)
      endif
c
c  set band limits of current solution
c
      if (nref.eq.1) then
        irow1(1) = ii1z+1
        irow2(1) = inn-ii2z
      else
        do 54 i = 1,nref-1
        irow1(i) = (rfiref(i+1)-rfic)/dfi+1.5
54      irow2(i) = (rfiref(i)-rfic)/dfi+0.5
        irow2(1) = irow2(1)+1
        irow1(nref) = -999
        irow2(nref) = -999
      endif
c
      if (lfirst) write(*,55) inn
55    format(/' Central latitudes and interpolation bands in ',
     .'the internal grid of',i4,' rows:')
      if (lfirst) 
     .write(*,56) (rfiref(i)*radeg,irow1(i),irow2(i),i=1,nref)
56    format(' ',f10.5,2i8)
c
c  main loop: do fft for each reference latitude
c  merge with previous solution and output interpolated band
c  ---------------------------------------------------------
c
      do 100 iref = 1, nref
        reffi = rfiref(iref)
c
c  set cosine and sine arrays etc
c
        if (iref.eq.1) then
          do 57 i = 1, nyqn
57        sindf2(i) = sin(dfi/2*i)**2
          do 58 i = 1, nyqe
58        sindl2(i) = sin(dla/2*i)**2
          do 59 i = 1, nyqn
59        sindfi(i) = sin(dfi*i)
        endif
c
        cosfir = cos(reffi)
        sinfir = sin(reffi)
        cosfi2 = cosfir**2
        csfir = cosfir*sinfir
c
c  set kernel of transform operator
c  --------------------------------
c
c  utilize that kernel is symmetric around longitude
c  top and left part of kernel array is negative delta fi and delta la
c  arrays are stored from top to bottom
c  zero argument stored at northern border
c
        do 69 i = 1, inn
c
c  set first sin(psi/2)
c
          if (i.eq.inn) then
            sindf = 0
            sinf2 = 0
          elseif (i.lt.nyqn) then
            sindf = sindfi(i)
            sinf2 = sindf2(i)
          else
            sindf = -sindfi(inn-i)
            sinf2 = sindf2(inn-i)
          endif
          do 69 j = 1, nyqe
            jj = ine+2-j
            k = (inn-i)*ine+j
            kk = (inn-i)*ine+jj
            if (j.eq.1) then
              sinl2 = 0
            else
              sinl2 = sindl2(j-1)
            endif
c
            s = sqrt(sinf2 + sinl2*(cosfi2*(1-2*sinf2)+csfir*sindf))
            lzero = (i.eq.inn.and.j.eq.1)
c
c  mode 1: stokes operator
c
            if (mode.eq.1) then
              if (lzero.or.s.le.spsi1.or.s.gt.spsi2) then
                rr = 0
              else
		rr = stokes(s,nmod1,nmod2)
              endif
              ri = 0
            elseif (mode.eq.3) then
c
c  mode 3: upward continuation by poissons integral
c
              if (lzero.or.psi.lt.psi1.or.psi.gt.psi2) then
                rr = 0
              else
                psi = asin(s)*2
                r0 = sqrt((re*psi)**2 + height**2)
                rr = height/r0**3
              endif
              ri = 0
            elseif (mode.ge.10) then
c
c  terrain effects
c
              if (psi.lt.psi1.or.psi.gt.psi2) then
                rr = 0
                ri = 0
              else
                psi = asin(s)*2
                r0 = re*psi
                goto (60,61,62,63),mode-9
c          
c  mode 10 - isostatic root
c
60              r0 = sqrt(r0**2 + diso**2)
                rr = diso/r0**3
                ri = -0.5d0*(r0**2 - 3*diso**2)/r0**5
                goto 68
c
c  mode 11 - airborne bouguer effect
c
61              if (lzero) then
                  rr = 0
                  ri = 0
                else
                  r0 = sqrt(r0**2 + height**2)
                  rr = height/r0**3
                  ri = 0.5d0*(r0**2 - 3*height**2)/r0**5
                endif
                goto 68
c
c  mode 12 - geoid effect of mass layer
c
62              if (lzero) then
                  rr = 0
                else
                  rr = 1/r0
                endif
                ri = 0
                goto 68
c
c  mode 13 - geoid RTM effect to 3rd order 
c
63              if (lzero) then
                  rr = 0
                  ri = 0
                else
                  if (lfirst) then
                    rr = 1/r0 
                    ri = 1.d0/6/r0**3
                  else
                    rr = 1.d0/6/r0**3
                    ri = 0
                  endif
                endif
              endif
            endif
c
c  complete longitude symmetry
c
68          cha(1,k) = rr
            if (j.gt.1) cha(1,kk) = rr
            cha(2,k) = ri
            if (j.gt.1) cha(2,kk) = ri
69      continue
c
c  transform kernel, convolve with data and do inverse transform
c  -------------------------------------------------------------
c
c
        call fourt(cha,nnarr,2,-1,icomplex,wrk,idim2,iwkdim)
c
        rewind(31)
        do 70 i = 1, n
          read(31) rr,ri
          cr = cha(1,i)
          ci = cha(2,i)
          cha(1,i) = (rr*cr - ri*ci)/n
          cha(2,i) = (ri*cr + rr*ci)/n
70      continue
c
        call fourt(cha,nnarr,2,1,1,wrk,idim2,iwkdim)
c
c  test
c        write(*,*) 're' 
c        do 777 ii=1,inn
c        kk=(ii-1)*inn
c777     write(*,778) (cha(1,kk+j),j=1,ine)
c
c  cha now contains the convolution 
c  note: in complex mode real part of (f1+if2)*(g1+ig2) is f1*g1-f2*g2
c  mode 13: store the result and go back to make two more terms
c  ------------------------------------------------------------
c
        if (lref) then
          if (lfirst) then
            rewind(31)
            rewind(35)
            do 72 i = 1,n
              write(36) cha(1,i) 
72          continue
            lfirst = .false.
            goto 471
          else
            rewind(35)
            rewind(36)
            do 73 i = 1,n
              read(35) h,href
              read(36) r
              cha(1,i) = r + 3*h*cha(1,i) - 3*h**2*cha(2,i)
73          enddo 
          endif
        endif
c
c  check solution against possible previous solution
c  output necessary part of solution after poss linear interpolation
c  -----------------------------------------------------------------
c
c  first output label on output grid
c
        if (iref.eq.1)
     .  write(33,80) rfi1*radeg, rfi2*radeg,
     .  rla1*radeg, rla2*radeg, dfi*radeg, dla*radeg
80      format(' ',4f12.6,2f12.7)
c
        jj1 = jj1z+1
        jj2 = ine-jj2z
c
c  merge previous and current solution and output
c
        if (iref.eq.1.and.nref.gt.1) goto 89
        if (iref.gt.1) then
          ii1 = irow1(iref-1)
          ii2 = irow2(iref-1)
        else
          ii1 = irow1(1)
          ii2 = irow2(1)
        endif
        rewind(32)
c
        iihalf = (ii2+ii1)/2
        do 82 i = ii2,ii1,-1
          rfi = rfic + (i-1)*dfi
          k = (inn-i)*ine
          if (iref.gt.1) then
            r2 = (rfiref(iref-1)-rfi)/(rfiref(iref-1)-rfiref(iref))
            r1 = (rfi-rfiref(iref))/(rfiref(iref-1)-rfiref(iref))
            csum = 0.d0
            dsum = 0.d0
            do 83 j = jj1,jj2
              read(32) cc
              cc1 = cha(1,k+j)
              if (i.eq.iihalf) then
                csum = cc**2+csum
                dsum = (cc-cc1)**2+dsum
              endif
              cha(1,k+j) = r1*cc + r2*cc1
83          continue
            if (i.eq.iihalf) then
              j = jj2-jj1+1
              write(*,84) rfi*radeg,sqrt(csum/j),sqrt(dsum/j)
84            format(' Latitude',f9.4,', r.m.s. signal and difference'
     .        ,2e12.4)
            endif
          endif
c
c  write out mixed band, convert units, include innerzone corrections
c  ------------------------------------------------------------------
c
          if (mode.eq.1) then
            sinfir = sin(rfi)**2
            gamma = 978032.d0*(1 + 0.00528*sinfir)
            rr = re*dfi*dla/(4*pi*gamma)
          elseif (mode.eq.3) then
            rr = re**2*dfi*dla/(2*pi)
          elseif (mode.eq.10) then
            rr = re**2*0.00667*drho*dfi*dla
          elseif (mode.eq.11) then
            rr = re**2*grho*dfi*dla
          elseif (mode.eq.12.or.mode.eq.13) then
            rr = re**2*grho*dfi*dla/gamma
          endif
c
          do 86 j=jj1,jj2
            g = cha(1,k+j)*rr
            if (liz) then
              read(30) s
              rmsiz = rmsiz+s**2
              if (abs(s).gt.rmaxiz) rmaxiz = abs(s)
              g = g+s
            endif
            if (g.lt.gmin) gmin = g
            if (g.gt.gmax) gmax = g
            cha(1,k+j) = g
86        continue
          write(33,85) (cha(1,k+j),j=jj1,jj2)
85        format(/,60(/,' ',8f9.3))
82      continue
c
c  store next band for next interpolation
c
89      if (iref.lt.nref) then
          rewind(32)
          do 90 i = irow2(iref),irow1(iref),-1
            k = (inn-i)*ine
            do 90 j = jj1,jj2
              write(32) cha(1,j+k)
90        continue
        endif
100   continue
c
c  exit
c  ----
c
      write(*,110) rfi1*radeg,rfi2*radeg,rla1*radeg,rla2*radeg,
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
      close(31,status='delete')
      close(32,status='delete')
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
      subroutine rdgrid(iunit,rfic,rlac,inn,ine,dfi,dla,
     .lgeo,cha,ik,ii1z,ii2z,jj1z,jj2z,idim)
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
c  if inn=0 the complete grid will be read.
c  if the wanted grid is too large a zero padding will be done.
c
c  last updated jun 90, rf
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      common /gridpar/ iell,izone
      dimension cha(2, idim)
      dimension hlab(6), hrow(2000)
      logical lgeo, lutm
      irdim = 2000
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
      if (inn.eq.0) then
        inn = nn
        ine = ne
      endif
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
      real*8 function stokes(s,nmod1,nmod2)
      implicit real*8(a-h,o-z)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  evaluates stokes function with wong-gore modification of low
c  harmonics, using recursion algorithms for pn
c  100% modification up to nmod1, then linearly tapering to 0% at nmod2
c  s is sin(psi/2)
c
c  rf aug 97 / feb 2002
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      rr = 1.d0/s-4-6*s+10*s**2-(3-6*s**2)*log(s+s**2)
      if (nmod2.le.1) then
        stokes = rr
      else
        t = cos(asin(s)*2)
        pn2 = 1.d0
        pn1 = t
        sum = 0.d0
        alfa = 1
        dr = nmod2-nmod1
        do 10 n = 2,nmod2
          if (n.gt.nmod1) alfa = (nmod2-n)/dr
          pn = ((2*n-1)*t*pn1 - (n-1)*pn2)/n
          sum = alfa*pn*(2*n+1)/(n-1) + sum
          pn2 = pn1
          pn1 = pn
10      continue
        stokes = rr-sum
      endif
      return
      end
