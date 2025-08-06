      program gradj
c $Id: gradj.for 181 2008-08-14 23:19:51Z tjansson $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      G R A D J
c
c  Adjustment program for gravimeter observations.
c  The program uses one observation equations for each
c  reading. Readings are assumed uncorrelated, which means
c  that large adjustments take up only limited space, since first
c  part of normal equations is diagonal.
c  The program solves for bias values, drifts and/or scale factors,
c  with basic obs.eq.
c
c         obs  =  1/scale*(g + bias + drift*t)
c
c  Stations which are only observed once do not take part in the
c  adjustment, but are computed from adjusted parameters.
c
c  Standard deviations of the readings are adjusted to match the
c  estimated standard deviations, while absolute errors are kept
c  fixed. This means that the adjustment requires a number of loops
c  to stabilize.
c
c  Input: 
c
c  <redobs-file>
c  <ofile>
c  <residualfile>
c  <tiefile>
c  dtmax, lsc, ldrift,
c  nfix,
c  stat1, gfix1, sigma,
c  stat2, gfix2, ...
c  ....
c  # G-<inst>
c  code, seqno  
c  ...
c  # G-<inst>
c  code, seqno
c  ...
c
c  where
c
c  <redobs-file>  file of reduced observations with a sequence
c                 number, as provided by 'grredu'.
c                 the file may contain station names.
c
c  <ofile>        outputfile with gravity results and scale factors.
c
c  <residualfile> residual file with tares and drift parameter
c                 only stations observed more than once are listed in this
c                 file.
c
c  <tiefile>      file with ties and residuals on differences
c                 between stations in adjustment (no tie when tare).
c                 this file may be used as input for tie adjustment
c                 programs, or may be sorted to show repeat ties
c                 (dos or unix: sort <tiefile >tiefile.sor)
c
c  dtmax ...  maximum time (hrs) allowed between two measurements
c             a tare is set automatically by program for longer intervals
c             (special: dtmax negative: old format input)
c  lsc  ....  t = solve for scale factors, f = no scale factors
c  ldrift ..  t = solve for linear drift, f = adjustment with bias only
c
c  The headers (# G-...) are mandatory, and the first 7 characters
c  with the inst.no. must correspond exactly to the <redobs-file>
c
c  The optional codes refer to the sequence numbers from 'grredu':
c  d  new drift and bias parameter from observation 'seqno'
c  t  assume tare occurred just before 'seqno' (i.e. new bias)
c  s  skip observation 'seqno'
c  Each 'seqno' must only have one code
c  NB: No blanks are allowed in front of the codes.
c
c  rf june 1993, completed aug 94, Resolute, NWT
c  last updated dec 27, 1994
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
c
c  maxn ... maximal number of unknowns and stations
c  maxc ... max dimension of normal equation matrix
c
      parameter (maxc=8000)
      parameter (maxn=800,maxfix=50,maxlist=50,maxcode=200)
c
      dimension islist(maxn),ia(maxn),lslist(maxn),ix(maxn),
     .g(maxn),sig(maxn)
      character*12 stxt(maxn)
      dimension instlist(maxlist),t0(maxlist)
      dimension isfix(maxfix),gfix(maxfix),sigfix(maxfix),iifix(maxfix)
      dimension iseqc(maxcode),b(maxcode),ib(maxcode)
      dimension c(maxc),sol(maxn)
      logical lslist,lsing(maxn),lsc,ldrift,lold
      character*72 ifile,ofile,rfile,tfile
      character*80 hdr,hdr1
      character*12 obstxt
      character*13 cplot(15)
      character*1 code(maxcode),ch
      common /cpar/ nastat,nbias,ndrift,nsc,iidrift,iisc,iirhs,
     .ldrift,lsc
c
      write(*,*) '--- GRADJ vers. dec 1994 ---' 
      write(*,*) '.. input: red-obs file'
      write(*,*) '          out-tyngde file'
      write(*,*) '          out-resi file'
      write(*,*) '          out-ties file'
      read(*,'(a)') ifile
      open(10,file=ifile,status='old')
      read(*,'(a)') ofile
      read(*,'(a)') rfile
      read(*,'(a)') tfile
      open(30,file=ofile,status='unknown')
      open(31,file=tfile,status='unknown')
      open(32,file=rfile,status='unknown')
c
      write(*,*) '.. dtmax(hr), lsc(f=std.scale), ldrift(t=find drifts)'
      read(*,*) dtmax,lsc,ldrift
      write(*,*) '.. nfix / list of: no,g,sigma (1 .. nfix)'
      read(*,*) nfix
      if (nfix.gt.maxfix) stop 'too many fixed stations'
      do k=1,nfix
        read(*,*) isfix(k),gfix(k),sigfix(k)
      enddo
      write(*,*) '.. list headers in file (e.g. # G-867) + tare list'
      write(*,*)
c
      lold = (dtmax.lt.0)
      dtmax = abs(dtmax)
      iobs = 0
      nobs = 0
      istat = 1
      ibias = 0
      idrift = 0
      nlist = 0
      nsc = 0
      nstat = 0
      nastat = 0
      open(20,form='unformatted',status='scratch')
c
c  reset all gravity by 981000 internally in program
c
      gbias = 981000.d0
      do k = 1,nfix
        gfix(k) = gfix(k) - gbias
      enddo
c
c  read observation file - may be with multiple labels
c  output all observations on tmp unit 20
c  set pointers to scale, bias and drift
c  ----------------------------------------------------
c
10    read(10,'(a)',end=30) hdr
      if (hdr(1:12).eq.'            ') goto 10
c
      if (hdr(1:1).eq.'#') then
c
c  read hdr line and hdr line info in input file
c  ---------------------------------------------
c
        if (nlist.gt.0) write(*,*) 'number of observations: ',iobs
        write(*,*) hdr(1:79)
        if (hdr(3:4).ne.'G-') stop '** header lacks G- in column 3-4'
        inst = iget(hdr,5)
        iobs = 0
        nlist = nlist+1
c
c  find index for instruments and scale factors
c
        if (nsc.eq.0) then
          nsc = 1
          instlist(1) = inst
          isc = 1
        else
          isc = listno(instlist,nsc,inst)
          if (isc.eq.0) then
            nsc = nsc+1
            if (nsc.gt.maxlist) stop '*** too many instruments'
            instlist(nsc) = inst
            isc = nsc
          endif
        endif
c
c  read drift/tare/skip info in input file
c
        if (nlist.eq.1) read(*,'(a)') hdr1
        if (hdr1(1:4).ne.'# G-') 
     .  stop '*** "# G-" not in first 4 char in job input file'
        inst1 = iget(hdr1,5)
        if (inst1.ne.inst)
     .  stop '*** instrument nos different in input and data file'
        ncode = 1
c
15      read(*,'(a)',end=18) hdr1
          ch = hdr1(1:1)
          if (ch.eq.'#') goto 18
          code(ncode) = ch
          iseqc(ncode) = iget(hdr1,2)
          ncode = ncode+1
          if (ncode.gt.maxcode) stop '*** too many codes'
          if (.not.(ch.eq.'t'.or.ch.eq.'d'.or.ch.eq.'s'))
     .    stop '*** code not identical to t, d, or s read from input'
          goto 15
c
18      ncode = ncode-1
        do j = 1, ncode
          write(*,'(1x,a1,i4)') code(j),iseqc(j)
        enddo
      else
c
c  read observation line
c  ---------------------
c
        call gettxt(hdr,obstxt)
        backspace(10)
c20      if (iobs.gt.1) isnoprev = isno
c
21      if (lold) then
          read(10,*,end=30) iseq,isno,i1,i2,i3,t,obs
          iday = i1*10000+i2*100+i3
        else
          read(10,*,end=30) isno,iday,t,iseq,r,r,obs
        endif
c
        iobs = iobs+1
        nobs = nobs+1
c
        if (iobs.gt.1) then
          ttprev = tt
          tt = epoch1900(iday,t) - t0(isc)
          dt = tt-ttprev
          if (dt.lt.0) write(*,22) iday,t
22        format(' *** warning: time reverse at ',i7,f6.2)
        else
          t0(isc) = epoch1900(iday,t)
          tt = 0.0
          dt = 0.0
        endif
c
c  find index for bias, drift etc.
c
        ic = listno(iseqc,ncode,iseq)
        if (ic.gt.0) then
          if (code(ic).eq.'s') then
            iobs = iobs-1
            nobs = nobs-1
            goto 21
          endif
          if (code(ic).eq.'t'.or.code(ic).eq.'d') ibias = ibias+1
          if (code(ic).eq.'d') idrift = idrift+1
        elseif (dt.gt.dtmax/24) then
          ibias = ibias+1
        elseif (iobs.eq.1) then
          ibias = ibias+1
          idrift = idrift+1
        endif
c
c  find station number index
c  islist ... list of all stations (1->nstat)
c  lslist ... t if station in adjustment (1->nstat)
c  check for fixed stations
c  set preliminary g-values from LC&R offset
c
        istat = listno(islist,nstat,isno)
        if (istat.eq.0) then
          nstat = nstat + 1
          if (nstat.gt.maxn) 
     .    stop '*** too many stations, increase maxn'
          islist(nstat) = isno
          stxt(nstat) = obstxt
          g(nstat) = obs - 4520
          lslist(nstat) = .false.
          if (listno(isfix,nfix,isno).gt.0) lslist(nstat) = .true. 
          istat = nstat
        else
          lslist(istat) = .true.
          if (stxt(istat).eq.'            ') stxt(istat) = obstxt
        endif
c
        write(20) iseq,isno,iday,t,obs,tt,isc,ibias,idrift,istat
c  
      endif
      goto 10
30    write(*,*) 'number of observations: ',iobs
c
c  set adjustment station index ia(istat) is pointing to adj.element 
c  -----------------------------------------------------------------
c  
      nastat = 0
      do 32 i = 1,nstat
        if (lslist(i)) then
          nastat = nastat+1
          ia(i) = nastat
        else
          ia(i) = 0
        endif
32    continue
c
c  output adjustment parameter info
c  --------------------------------
c
      ndrift = idrift
      if (.not.ldrift) ndrift = 0
      if (.not.lsc) nsc = 0
      nbias = ibias
c
      n = nastat + nbias + ndrift
      if (lsc) n = n+nsc 
      icdim = (n*(n+1)-nastat*(nastat-1))/2+n
      iirhs = icdim-n
      iidrift = nastat+nbias
      iisc = iidrift+ndrift
c
      write(*,34) nobs,nstat
34    format(' Totally read:',i6,' observations in',i5,' stations')
      write(*,35) nastat,nbias,ndrift,nsc,n
35    format(/' Number of unknowns in adjustment:',
     ./' Station gravity values:  ',i5,
     ./' Bias parameters:         ',i5,
     ./' Drift parameters:        ',i5,
     ./' Scale factors:           ',i5,
     ./' Total number of unknowns:',i5)
      if (icdim.gt.maxc) stop '*** adjustment too large, increase maxc'
      if (n.gt.maxn) stop '*** too many unknowns, increase maxn'
c
c  test output
c
c     rewind(20)
c     do 50 i = 1,nobs
c       if (i.eq.1) write(*,*) 'test: binary obs ofile:'
c       read(20) iseq,isno,iday,t,obs,tt,isc,ibias,idrift,istat
c       write(*,37) iseq,isno,tt,isc,ibias,idrift,istat,ia(istat)
c37     format(2i7,f11.4,5i7)
c50   continue
c
c  set pointers for fixed stations
c  -------------------------------
c
      do 53 i = 1, nfix
        ii = listno(islist,nstat,isfix(i))
        iifix(i) = ii
53    continue
c
c  main loop for normal equations
c  ------------------------------
c
c  first set initial values and loop counter
c  nastat is length of diagonal, n total number of unknowns
c
      maxloop = 10
      if (lsc) maxloop = 20
c
      iloop = 0
      do 55 i = 1,nstat
        if (ia(i).gt.0) sol(ia(i)) = g(i)
55    continue
      do 56 i = nastat+1,iidrift
56    sol(i) = 4520
      do 57 i = iidrift+1,iisc
57    sol(i) = 0.d0
      do 58 i = iisc+1,n
58    sol(i) = 1.d0
c
c  apriori variance of single observation = sigma = 0.02
c  this variance is updated in each loop
c
      sigma = 0.02
      write(*,59) sigma
59    format(/' Adjustment, initial assumed std.dev. of readings:',
     .f7.2)
c
c  adjustment loop entry
c  ---------------------
c
60    rewind(20)
      iloop = iloop+1
c
      naobs = 0
      pvvsum = 0
      gvvsum = 0
      if (iloop.le.maxloop) then
        do 61 i = 1, icdim
61      c(i) = 0
      endif
c
c  scan through reading observations
c  ---------------------------------
c
70    read(20,end=75) iseq,isno,iday,t,obs,tt,isc,ibias,idrift,istat
      iastat = ia(istat)
      if (iastat.eq.0) goto 70
      naobs = naobs+1
c
c  set obs eq and add to normal eq
c
      s0 = 1.d0
      if (lsc) s0 = sol(iisc+isc)
      pobs = predobs(sol,iastat,ibias,idrift,tt,isc)
c
      b(1) = 1/s0
      ib(1) = iastat
      b(2) = 1/s0
      ib(2) = nastat+ibias
      j = 3
      if (ldrift) then
        b(j) = tt/s0
        ib(j) = iidrift+idrift
        j = j+1
      endif
      if (lsc) then
        b(j) = -pobs/s0
        ib(j) = iisc+isc
        j = j+1
      endif
      b(j) = obs - pobs
      ib(j) = n+1
      gvvsum = gvvsum + b(j)**2
c
      pvvsum = pvvsum + b(j)**2
c
      if (iloop.le.maxloop) call addobsd(c, nastat, n, b, ib)
      goto 70
c
75    if (iloop.eq.1) write(*,*)
     .'Number of readings in adjustment:      ',naobs
      if (naobs.gt.0) gvvsum = sqrt(gvvsum/naobs)
c
c  scan through fixed obs
c  ----------------------
c  limit weight when sigma very large
c
      do 77 i = 1, nfix
        if (iifix(i).eq.0) goto 77
        naobs = naobs+1
        weight = sigma/sigfix(i)
        if (weight.gt.100) weight = 100
        j = ia(iifix(i))
        if (j.eq.0.or.j.gt.nastat) stop '*** error in fixed obs prog'
        b(1) = 1.d0 * weight
        ib(1) = j
        b(2) = (gfix(i)-sol(j))*weight
        ib(2) = n+1
        if (iloop.le.maxloop) call addobsd(c, nastat, n, b, ib)
        pvvsum = pvvsum + b(2)**2
77    continue
      if (iloop.eq.1) write(*,*)
     .'Total number of observation equations: ',naobs
      if (naobs.lt.n) stop '*** too few observations'
c
c  statistics and residuals, exit last loop
c  ----------------------------------------
c
      if (nastat.gt.0.and.iloop.gt.1) write(*,78) gvvsum
78    format(' r.m.s. reading residuals (mgal): ',f14.3)
      sigma = -999.999     
      if (naobs-n.gt.0) sigma = sqrt(pvvsum/(naobs-n))
      if (iloop.gt.1) write(*,79) sigma
79    format(' SIGMA single reading = ',f14.3,' mgal')
      if (iloop.gt.maxloop) goto 100
c
      write(*,*)
      write(*,*) 'SOLUTION LOOP: ',iloop
c
c  solve equations
c  ---------------
c
c      write(*,791) (c(kk),kk=1,icdim)
c791   format(' test norm eq matrix'/,8(5e14.6,/))
c
      call chold(c,n,nastat,nsing,cif,lsing)
c
c      write(*,792) iirhs,(c(kk),kk=1,icdim)
c792   format(' test after chold iirhs =',i5,/,8(5e14.6,/))
c
      if (nsing.eq.0) then
        write(*,82) cif
82      format(' - solution OK, max loss of digits:',f5.1)
      else
        write(*,*) '- number of normal equation singularities: ',nsing
        do 83 i= 1,n
83      if (lsing(i)) write(*,*) '  singularity row no.: ',i
      endif
c
c  update solution vector
c  ----------------------
c
      solmax = -9.9d9
      solrms = 0
      do 90 i = 1, n
        r = c(i+iirhs)
        if (r.gt.solmax) solmax = r
        solrms = solrms + r**2
        sol(i) = sol(i) + r
90    continue
      if (iloop.gt.1) write(*,91) sqrt(solrms/n), solmax
91    format(' r.m.s. and max parameter change: ',2f14.7)
      if (solmax.lt.1.d-6) iloop = maxloop+1
c
c  go back to make residuals in start of next loop
c  -----------------------------------------------
c
      goto 60
c
c  exit of loops: make inverse for standard deviations or find single g
c  --------------------------------------------------------------------
c
100   rewind(20)
      do 110 i = 1, nstat
        if (lslist(i)) then
          j = ia(i)
          g(i) = sol(j)
          sig(i) = sqrt(cholinv(c,n,nastat,j))*sigma
        else
105       read(20) iseq,isno,iday,t,obs,tt,isc,ibias,idrift,istat
          if (isno.ne.islist(i)) goto 105
c
          s = 1.0
          if (lsc) s = sol(iisc+isc)
          r = sol(nastat+ibias)
          if (ldrift) r = r + sol(iidrift+idrift)*tt
          g(i) = s*obs - r
          sig(i) = 99.999
        endif
110   continue
c
c  sort station list
c  -----------------
c
      write(*,*)
      write(*,*) '- sorting results, output files: '
      write(*,*) ofile
      write(*,*) rfile
      write(*,*) tfile
      call sortin(islist,nstat,ix)
c
      write(30,*) '#== Fixed stations and adjustment residuals ==='
      write(30,*) '#   stat      fix g     sigma      adj g        v '
      do i = 1, nfix
        if (iifix(i).gt.0) then 
          write(30,140) isfix(i),gfix(i)+gbias,sigfix(i),
     .    sol(ia(iifix(i)))+gbias,gfix(i)-sol(ia(iifix(i))),
     .    stxt(ia(iifix(i)))
140       format(i8,f14.3,f7.3,' ',f12.3,f8.3,2x,a12)
        else
          write(30,141) isfix(i),gfix(i)+gbias,sigfix(i)
141       format(i8,f14.3,f7.3,'   not observed ')
        endif
      enddo
c
      write(30,*)
      write(30,*)
     . '#== Adjusted new gravity values and standard deviations ==='
      do 120 i = 1, nstat
        if (listno(isfix,nfix,islist(ix(i))).gt.0) goto 120
        if (abs(sig(ix(i))-99.999).gt.0.01) then
          write(30,121) islist(ix(i)),g(ix(i))+gbias,sig(ix(i)),
     .    stxt(ix(i))
        else
          write(30,122) islist(ix(i)),g(ix(i))+gbias,stxt(ix(i))
        endif
121     format(i8,f14.3,f9.3,2x,a12)
122     format(i8,f14.3,'      N/A',2x,a12)
120   continue
c
      if (lsc) then
        write(30,*)
        write(30,*) '=== Scale factors of gravimeters ==='
        write(30,*) '  inst.no.   scale      error'
        do i = 1, nsc
          write(30,132) instlist(i),sol(iisc+i),
     .    sqrt(cholinv(c,n,nastat,iisc+i))*sigma
132       format(i8,2f11.6)
        enddo
      endif
c
      write(30,*)
      write(30,*) '=== Statistics of adjustment ==='
      write(30,138) naobs,nastat,n,sigma
138   format(' Adjustment observations:',i6/
     .' Stations:',i6,', total unknowns: ',i6/
     .' SIGMA (single reading at apriori weighting): ',f8.3)
c
c
      write(32,*)
      write(32,*) '=== Reading residuals with tares and drift ==='
      write(32,*) '(excl. single measurements not in adjustment)'
      write(32,*)
      write(32,*)
     .'  statno      t         seqno    obs        v  ',
     .' plot (0.02 mgal/div)'
c 
      cplot(1) = 'M     !      '
      cplot(2) = '*     !      '
      cplot(3) = ' *    !      '
      cplot(4) = '  *   !      '
      cplot(5) = '   *  !      '
      cplot(6) = '    * !      '
      cplot(7) = '     *!      '
      cplot(8) = '      *      '
      cplot(9) = '      !*     '
      cplot(10)= '      ! *    '
      cplot(11)= '      !  *   '
      cplot(12)= '      !   *  '
      cplot(13)= '      !    * '
      cplot(14)= '      !     *'
      cplot(15)= '      !     M'
c
      isc0 = 0
      ibias0 = 0
      idrift0 = 0
      isno00 = 0
      ibias00 = 0
      ibtw = 0
      nv = 0
      gvvsum = 0
      rewind(20)
c
c  output residual loop
c  --------------------
c
      write(31,161)
161   format(4x,
     .' from/to     inst      time   m.betw. dt(hr)    dif        v')
      nsing = 0
c
150   read(20,end=170) iseq,isno,iday,t,obs,tt,isc,ibias,idrift,istat
      iastat = ia(istat)
c
      if (isc.ne.isc0) then
        if (nv.gt.0) then
          write(32,151) sqrt(gvvsum/nv)
151       format(' r.m.s. residuals: ',f8.3)
          nv = 0
          gvvsum = 0
        endif
        write(32,152) instlist(isc)
152     format(/' # G-',i3)
      endif
c
      if (ldrift.and.idrift.ne.idrift0) then
        write(32,155) sol(iidrift+idrift)
155     format(' DRIFT parameter:',f8.3,' mgal/day')
        idrift0 = idrift
      endif
c
      if (ibias.ne.ibias0) then
        if (isc.ne.isc0) write(32,156) sol(nastat+ibias)
        if (isc.eq.isc0) then
          write(32,157) sol(nastat+ibias),
     .    sol(nastat+ibias)-sol(nastat+ibias-1)
        endif
156     format(' BIAS parameter:',f12.3)
157     format(' BIAS parameter:',f12.3, ', TARE:',f9.3)
      endif
c
      if (iastat.ne.0) then
        if (nsing.gt.0) then
          write(32,158) nsing
158       format('         (',i2,' meas.)')
          nsing = 0
        endif
        v = obs-predobs(sol,iastat,ibias,idrift,tt,isc)
        nv = nv+1
        gvvsum = gvvsum + v**2
        i = nint(v/0.02)
        if (i.lt.-7) i = -7
        if (i.gt.7) i = 7
        write(32,159) isno,iday,t,iseq,obs,v,cplot(i+8),stxt(istat)
159     format(i7,i8,',',f6.2,i6,2f10.3,'  ',a13,1x,a12)
      else
        nsing = nsing+1
      endif
c
c  output ties
c
      if (iastat.eq.0) then
        ibtw = ibtw+1
        goto 165
      endif
      if (isno00.ne.isno.and.ibias.eq.ibias00) then
        if (isno00.le.isno) then
          write(31,162)
     .    isno00,isno,instlist(isc),iday,t,ibtw,(tt-tt00)*24,
     .    obs-obs00,v-v00
        else
          write(31,162)
     .    isno,isno00,instlist(isc),iday,t,ibtw,(tt-tt00)*24,
     .    obs00-obs,v00-v
        endif
      endif
162   format(2i7,'  G-',i3,'  ',i6,',',f5.2,i4,f7.1,f10.3,f9.3)
      isno00 = isno
      ibias00 = ibias
      tt00 = tt
      v00 = v
      obs00 = obs
      ibtw = 0
c
165   ibias0 = ibias
      isc0 = isc
      idrift0 = idrift
      goto 150
c
170   if (nv.gt.0) write(32,151) sqrt(gvvsum/nv)
      stop 'GRADJ terminated OK'
      end
c
      real*8 function predobs(sol,iastat,ibias,idrift,tt,isc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                    P R E D O B S
c
c  yields predicted obs based on current contents of sol vector
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      logical ldrift,lsc
      common /cpar/ nastat,nbias,ndrift,nsc,iidrift,iisc,iirhs,
     .ldrift,lsc
      dimension sol(*)
c
      r = sol(iastat) + sol(nastat+ibias)

      if (ldrift) r = r + sol(iidrift+idrift)*tt
      if (lsc) r = r/sol(iisc+isc)
      predobs = r
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                    L I S T N O
c
c  subroutine for finding station number in unsorted station list
c  'idx' is zero if 'istat' is not in array 'ia' of 'n' elements
c  search begins at index from previous call
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function listno(ia, n, istat)
      dimension ia(*)
      save idx
      if (n.eq.0) goto 30
      if (idx.le.0.or.idx.gt.n) idx = 1
      j = idx
10    if (ia(j).eq.istat) goto 20
      j = j+1
      if (j.gt.n) j = 1
      if (j.eq.idx) goto 30
      goto 10
20    idx = j
      listno = idx
      return
30    idx = 0
      listno = idx
      return
      end
c
      subroutine sortin(d,n,ix)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                    S O R T I N
c
c  Subroutine sorts data in array d from 1 to n.
c  Result is given in integer array i, so that d(ix(k)), k = 1..n
c  is the sorted sequence.
c  The routine is fast for near-sorted data, slow otherwise.
c
c  (c) RF/KMS. Jan 1992
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer d,dd
      dimension d(*),ix(*)
c
      ix(1) = 1
      do 20 i = 2, n
        dd = d(i)
        j = i-1
10      if (dd.lt.d(ix(j))) then
          ix(j+1) = ix(j)
        else
          ix(j+1) = i
          goto 20
        endif
        j = j-1
        if (j.gt.0) goto 10
        ix(1) = i
20    continue
      return
      end
c
      real*8 function epoch1900(iday,time)  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                       E P O C H 1 9 0 0
c
c  The procedure compute the time interval in units of days
c  from epoch 1900 (A.D. 1900 January 0.5) to the epoch given
c  by the parameter which must be a GI-standard date.
c  iday is given as an integer in form '130591', meaning May 13, 1991.
c  time is a real giving hr and min (e.g. 17.22), and possible minute
c  decimals.
c
c  Example:
c  The Julian day number of run-time can be written on current
c  output by the call
c  write(out, entier epoch_1900(date_time) + 2415020)
c
c  Willy Lehmann Weng.
c  GI, October 1976, Reg-no 76044.
c  (c) KMS. Fortran version by Rene Forsberg, Oct 91
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 time, m
      integer iday, days, year, month, day, h, i, ii
c
      year  = mod(iday,100)  
      month = mod(iday/100,100)  
      day   = iday/10000  
      if (month.le.0.or.month.gt.12.or.
     .day.le.0.or.day.gt.31) then
        write(*,*) '*****',day,month,year
        stop 'epoch1900 - illegal date spec'
      endif
      i = 0
      if (year.gt.3.and.month.lt.3) i = 1
      goto (1,2,3,4,5,6,7,8,9,10,11,12),month
1     ii = 0
      goto 13
2     ii = 31
      goto 13
3     ii = 59
      goto 13
4     ii = 90
      goto 13
5     ii = 120
      goto 13
6     ii = 151
      goto 13
7     ii = 181
      goto 13
8     ii = 212
      goto 13
9     ii = 243
      goto 13
10    ii = 273
      goto 13
11    ii = 304
      goto 13
12    ii = 334
13    continue
c
      days = year*365 + (year - i)/4 + ii + day
      h = time
      m = (time-h)*100 
c
      epoch1900 = days + h/24.d0 + m/1440.d0 - 0.5d0  
c
      return
      end
c
      subroutine chold(c, n, ndia, nsing, cifmax, lsing)
      implicit double precision(a-h,o-z)
      dimension c(*)
      logical lsing(*)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                           C H O L D
c  
c  subroutine for solving a set of normal equations where the first
c  'ndia' rows are diagonal. the system is solved by choleskys method.
c  the equations are stored columnwise, upper diagonal only, starting
c  with the first 'ndia' diagonal elements only
c
c  
c  parameters:
c  
c  c       array. at call coefficients and right-hand side. On return c con-
c          tains the reduced matrix. normal equations may be build up
c          efficiently with subroutine 'addobsd'. the size of the array
c          must be at least (n*(n+1) - ndia*(ndia-1))/2
c  
c  n       total number of unknowns.
c          if n < 0 is specified an earlier solution
c          is utilized with a new right-hand side.
c
c  ndia    length of initial diagonal.
c
c  nsing   number of singularities, = 0 for ok solution
c
c  cifmax  the actual maximal loss of digits encountered.
c
c  lsing   logical array signalling singular rows.
c
c  (c) rf/kms, june 1990. conversion of algol procedure 'grnll'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer row, col, p, i, ir, ic
        logical findl
c
        cifmax = 0
        nsing = 0
        findl = (n.gt.0)
        n = abs(n)
c
        icol1 = 1
        if (.not.findl) icol1 = n+1
c  
        do 50 row = icol1, n+1
          ir = indx(row,ndia)  
          i = ir
          irow1 = 1
          if (row.le.ndia) irow1 = row
          do 40 col = irow1, row
            sum = 0
            ic = indx(col,ndia)
            i = i+1
            if (col.gt.ndia) then
              do 30 p = 1,col-1
30            sum = sum + c(ir+p)*c(ic+p)
            endif 
c
            ci = c(i)
            if (row.ne.col) then
              c(i)= (ci - sum)/c(indx(col+1,ndia))
            elseif (row.le.n) then
              ci1 = ci - sum
              cif = 0
              if (ci1.eq.0.or.ci.eq.0) then
                cif = 999.d0
              else
                cif = 0.434*log(abs(ci/ci1))
              endif
              if (cif.gt.cifmax) cifmax = cif
              if (ci1.le.0) then
                nsing = nsing+1
                lsing(row) = .true.
                c(i) = 9.9d33
              else
                lsing(row) = .false.
                c(i) = sqrt(ci1)
              endif
            endif
40        continue
50      continue
c
c  back substitution
c
        do 60 col = n, 1, -1
          i = i-1
          ir = i
          ic = indx(col+1,ndia)
          c(i) = c(i)/c(ic)
          if (col.le.ndia) goto 60
          do 58 p = col-1, 1, -1
            ir = ir-1
            ic = ic-1
            c(ir) = c(ir) - c(i)*c(ic)
58        continue
60      continue
        return
        end
c
        integer function indx(m, ndia)  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                          I N D X
c
c  indx is the number of matrix elements prior to column
c  'm', i.e. c(indx(m)+1) is first element of column 'm'
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (m.le.ndia) then 
          indx = m-1 
        else
          indx = (m*(m-1) + ndia*(1-ndia))/2
        endif
        return
        end
c
      double precision function cholinv(c, n, ndia, inv)
      implicit double precision(a-h,o-z) 
      dimension c(*)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                           C H O L I N V
c
c  this subroutine obtains inverse element no. 'inv' from the inverse
c  matrix of the normal equation matrix c. c must be cholesky reduced,
c  i.e. chold must have been called.
c  the right-hand side is not used, but modified on exit
c  
c  (c) rf/kms oct 86, converted to fortran june 90
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer p, col, c1, row
      logical diag
c
      diag = (inv.le.ndia)
      i = indx(n+1, ndia)
c
c  cholesky reduction of right-hand side  0, ..., 1, 0, ..
c
      if (diag) then
        c(i+inv) = 1.0/c(inv)
        row = ndia+1
      else
        c(i+inv) = 1.0/c(indx(inv+1,ndia))
        row = inv
      endif
c
      ir1 = row+1
      if (diag) ir1 = row
      do 20 col = ir1, n
        c1 = col-1
        ic = indx(col,ndia)
        if (diag) then
          sum = -c(i+inv)*c(ic+inv)
        else
          sum = 0
        endif
        do 18 p = row, c1
18      sum = sum - c(i+p)*c(ic+p)
        c(i+col) = sum/c(ic+col)
20    continue
c
c  back substitution
c
      do 30 col = n, row, -1
        ic = indx(col,ndia)
        rr = c(i+col)/c(ic+col)
        c(i+col) = rr
        do 28 p = col-1, row, -1
28      c(i+p) = c(i+p) - rr*c(ic+p)
        if (diag) c(i+inv) = c(i+inv) - rr*c(ic+inv)
30    continue
      if (diag) c(i+inv) = c(i+inv)/c(inv)
      cholinv = c(i+inv)
      return
      end
c
      subroutine addobsd(c, ndia, n, b, ib)
      implicit double precision(a-h,o-z)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                            A D D O B S D
c
c  subroutine adds one observation to normal equations
c
c  parameters:
c
c  c   normal equation matrix
c
c  ndia  length of initial diagonal
c  n   number of unknowns
c
c  b   observation equation coefficients which are non-zero.
c      right hand side (n+1) must be specified as last value.
c
c  ib  column numbers of coefficients in b, terminated by n+1
c      only the first value can be in the range 1 .. ndia
c
c  (c) rf/kms, 1986, converted to fortran june 1990
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension c(*), b(*), ib(*)
      i = 0
c
c  i-loop
c
10    i = i+1
c     write(*,*) 'test addobsd ib,b = ',ib(i),b(i)
      ir = ib(i)
      if (i.gt.1.and.ir.le.ndia) stop '*** obs eqn error ***'
      do 20 j = 1, i
        is = ib(j)
        if (ir.le.ndia) then
          ipos = ir
        else
          ipos = (ir*(ir-1)+ndia*(1-ndia))/2 + is
        endif
        c(ipos) = c(ipos) + b(i)*b(j)
20    continue
      if (ir.le.n) goto 10
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                        g e t t x t
c
c  subroutine read txt portion starting with letter in line,
c  terminates with 2 blanks
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gettxt(ch,txt)
      character*12 txt
      character*80 ch
      character*1 c,c1
c
      txt = '            '
      j = 0
10    j = j+1
      if (j.gt.80) goto 30
      c = ch(j:j)
      if ('a'.le.c.and.c.le.'z'.or.'A'.le.c.and.c.le.'Z') goto 20
      goto 10
c
20    j0 = j
      j12 = j0+13
21    j = j+1
      if (j.le.80) then
        c1 = ch(j:j)
      else
        c1 = ' '
      endif
      if (c.eq.' '.and.c1.eq.' '.or.j.ge.j12) goto 29
      c = c1
      goto 21
c
29    txt = ch(j0:j-2)
30    return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                     i g e t  
c
c  reads integer from array c from position 'ipos',
c  terminated by non-number
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function iget(ch,ipos)
      character ch*80
      character*1 c
c
      j = ipos-1
      isig = 1
c
10    j = j+1
      if (j.gt.80) then
        write(*,*) ch(1:79)
        stop '*** getint: no integer in line'
      endif
      c = ch(j:j)
      if (c.eq.' ') goto 10
      if (c.eq.'-') then
        isig = -1
        goto 10
      endif
      if (c.lt.'0'.or.c.gt.'9') then
        write(*,*) ch(1:79)
        stop '*** getint: wrong character in front of integer'
      endif
c
      i = ichar(c)-48
20    j = j+1
      c = ch(j:j)
      if (j.ge.80.or.c.lt.'0'.or.c.gt.'9') goto 30
      i = i*10 + ichar(c)-48
      goto 20
c
30    iget = isig*i
      return
      end
