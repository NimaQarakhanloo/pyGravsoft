      program crsadj 
c $Id: crsadj.for 248 2008-12-09 15:49:33Z cct $
c
c                            C R S A D J
c
c  Cross-over program for altimeter data, computing
c  cross-over differences, estimating bias/tilt parameters, and
c  correction of the altimeter data.
c
c  Input:
c  lat1,lat2,lon1,lon2,ltilt,mode
c  
c  Input altimeter data format:
c  mode 1:
c  rev.no, lat, lon, dummy, ssh, std.dev.
c  mode 2:
c  rev.no, t, lat, lon, ssh, std.dev.
c
c  Output files:
c  crsdif.dat, crspar.dat
c  
c  The ssh should preferably be minus an EGM like EGM96.
c
c  Programmed by Per Knudsen, KMS - Denmark, Oct 5, 1990.
c  Minor corrections by RF, University of Calgary, March 1993
c  Change 2002-07-12  and 2008-12.10 by cct.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8 (a-h,o-z)
      parameter(maxdat=100000,irest=63250)
      integer*2 kerr,kcode
      character*128 ifile,ofile
C
C COMMON SPACE OF P.T. 563250*4 BYTES (USED BY SUBROUTINE ABIAS)
C
      common/field/ krev(maxdat),lat(maxdat),lon(maxdat),
     .              kalt(maxdat),kerr(maxdat),kcode(maxdat),
     .              i4rest(irest),time(maxdat)
      common /nsave/ndat
      LOGICAL LSAVE,LTILT,LNOERR,LNPASS
c
      write(*,*)' C R S A D J, Ver. DEC 2008.'
      write(*,*)' Input name of input file '
      read(*,'(A)') ifile
      write(*,*)' Data input from file: ',ifile
      write(*,*)' Input name of output file '
      read(*,'(A)') ofile
      write(*,*)' Data output from file: ',ofile
      write(*,*)' input min, max lat & long, ltilt and mode '
      write(*,*)' mode 1: rev., lat, lon, nd, ssh, essh '
      write(*,*)' mode 2: track, time, lat, lon, ssh, essh '
      write(*,*)' mode 3: track, lat, lon, ssh, essh '
      write(*,*)' mode 4: cycle, lat, lon, s0, essh, time '
      write(*,*)' mode 5: cycle, lat, lon, h (=0), ssh '
      read(*,*) alat1,alat2,alon1,alon2,ltilt,mode
      if (mode.lt.1.or.mode.gt.5) stop '*** mode wrong'
c     READ(*,*) iuinc,WEIGHT
      iuinc = 0
      weight = 0.1
c
      write(*,*)' Lat & long boundaries used: '
      write(*,*) alat1,alat2
      write(*,*) alon1,alon2
      open(15,form='formatted',file=ifile)
      open(30,form='formatted',file=ofile)
      open(31,form='formatted',file='crsdif.dat')
      open(33,form='formatted',file='crspar.dat')
c  correction rf
      open(55,status='scratch',form='unformatted')
C
C********************************************************
C  READ DATA AND STORE IN COMMON/FIELD/
c
      ndat=0
      irev0=0
      idat=0
      tt0=0.0
      if (mode.eq.4.or.mode.eq.5) std=0.05
c
c  rf - input of altimeter data
c
  100 if (mode.eq.1) read(15,*,end=101) irev,alat,alon,tt,dh,std
      if (mode.eq.2) read(15,*,end=101) irev,tt,alat,alon,dh,std
      if (mode.eq.3) read(15,*,end=101) irev,alat,alon,dh,std
      if (mode.eq.4) then
        read(15,*,end=101) ir,alat,alon,h0,dh,tt  
        alon=-alon
        if ((tt-tt0).gt.10000) then
         irev=irev+1
        end if
        tt0=tt
      end if
      if (mode.eq.5) then
       read(15,*,end=101)irev,alat,alon,h0,dh
       std=0.05
      end if
      if(alat.lt.alat1.or.alat.gt.alat2) go to 100
      if(alon.lt.alon1.or.alon.gt.alon2) go to 100
c     if (abs(dh).gt.10) go to 100
      if(irev.ne.irev0) then
        irev0=irev
        if(idat.lt.10) ndat=ndat-idat
        idat=0
      endif
      ndat=ndat+1
      if (ndat.gt.maxdat) stop '*** too many data, sorry ***'
      idat=idat+1
      krev(ndat)=irev
      lat(ndat)=alat*1.0d6
      lon(ndat)=alon*1.0d6
      kalt(ndat)=dh*1.0d3
      kerr(ndat)=std*1.0d3
      go to 100
  101 continue
      write(*,*)' passed 101 '
      if(idat.lt.10) ndat=ndat-idat
c     ndat=ndat-1
c
C**********************************************************
C  COMPUTE CROSS-OVER DIFFERENCES AND STORE THEM ON UNIT 31
C
      iuout=31
      DIST0=0.0
      DSTMAX=0.0
      CALL ACRS(IUOUT,ALAT1,ALAT2,ALON1,ALON2,
     *          NDAT,DIST0,DSTMAX)
C
C************************************************************
C  COMPUTE STRAIGHT LINE REGRESSION FOR EACH TRACK.
C  NORMAL MATRIX ELEMENTS ARE STORED ON UNIT 55.
C  COMPUTE STATISTICS.
C  REDUCE DATA.
C
      IUOUT=55
      LSAVE=.TRUE.
      CALL ALINE(IUOUT,NDAT,LSAVE)
C
C**************************************************************
C  CHECK CORSS-OVERS. NO STORING.
C
c     IUOUT=0
c     CALL ACRS(IUOUT,ALAT1,ALAT2,ALON1,ALON2,NDAT,DIST0,DSTMAX)
c
C***************************************************************
C  PERFORME A CROSS-OVER ADJUSTMENT. WRITE BIAS/TILT PARAMENTERS
C  ON UNIT 33.
C
      IUI2=31
      LNOERR=.false.
      hst=0.1
      IUO2=33
      CALL ABIAS(IUI2,IUINC,IUO2,HST,WEIGHT,LTILT,LNOERR)
C
C***************************************************************
c  READ DATA AGAIN SINCE THEY ARE NOT IN COMMON ARRAY ANYMORE.
C
      ndat=0
      rewind(15)
      irev=0
      irev0=-1
      idat=0
c  - rf modification
  200 if (mode.eq.1) read(15,*,end=201) irev,alat,alon,ddd,dh,std
      if (mode.eq.2) read(15,*,end=201) irev,t,alat,alon,dh,std
      if (mode.eq.3) read(15,*,end=201) irev,alat,alon,dh,std
      if (mode.eq.4) then
        read(15,*,end=201) ir,alat,alon,h0,dh,tt  
        alon=-alon
        if ((tt-tt0).gt.10000) then
         irev=irev+1
        end if
        tt0=tt
      end if
      if (mode.eq.5) then
       read(15,*,end=201)irev,alat,alon,h0,dh
       std=0.05
      end if
      if(alat.lt.alat1.or.alat.gt.alat2) go to 200
      if(alon.lt.alon1.or.alon.gt.alon2) go to 200
c     if (abs(dh).gt.10) go to 200
      if(irev.ne.irev0) then
        irev0=irev
        if(idat.lt.10.and.idat.gt.0) then
          ndat=ndat-idat
          write(*,*)' ndat changed ',ndat,idat
        end if
        idat=0
      endif
      ndat=ndat+1
      idat=idat+1
c     time(ndat)= adum
      krev(ndat)=irev
      lat(ndat)=alat*1.0d6
      lon(ndat)=alon*1.0d6
      kalt(ndat)=dh*1.0d3
      kerr(ndat)=std*1.0d3
c     write(*,*)krev(ndat),lat(ndat),lon(ndat),kalt(ndat),
c    *kerr(ndat)
      go to 200
  201 continue
      write(*,*)' label 201 passed '
      if(idat.lt.10) ndat=ndat-idat
C
C***************************************************************
C  CORRECT DATA FOR BIAS/TILT
C
      IUIN2=33
      lnpass=.true.
      call abcor(iuin2,ndat,lnpass,ndat2)
      write(*,*) ' Finished abcor - entering Aplane',ndat2
C
C****************************************************************
C  ESTIMATE A PLANE TROUGH THE CORRECTED OBSERVATIONS
 
c     LSAVE=.TRUE.
c     CALL APLANE(NDAT2,LSAVE)
C
C*****************************************************************
C  IF CORRECTED OBSERVATIONS SHOULD BE WRITTEN, THEN (E.G.)
C
c
      do 10 i=1,ndat2
      iii=(krev(i)*10000)+i
      alat=lat(i)/1.d6
      alon=lon(i)/1.d6
      ssh=kalt(i)/1.d3
      err=kerr(i)/1.d3
c     write(30,500) krev(i),time(i),alat,alon,ssh,err
      if (mode.eq.4)alon=-alon
      write(30,500) krev(i),alat,alon,ssh,err
   10 continue
  500 format(i6,2f10.5,' 0 ',2f10.4)
      rewind(30)
      write(*,*)' CRSADJ TERMINATED '
C
      end
C *******************************************************************
C *****                                                         *****
C *****          SUBROUTINE ALINE                               *****
C *****                                                         *****
C *******************************************************************
C
      SUBROUTINE ALINE(IUOUT,NDAT,LSAVE)
C   *****************************************************************
C
C   PROGRAMS ARE COPYRIGHT BY THE AUTHOR AND KORT- OG MATRIKELSTYRELSEN.
C
C   SUBROUTINE ALINE COMPUTES FOR EACH TRACK BIAS/TILT PARAMETERS FROM
C   RESIDUAL SEA SUEFACE HEIGHTS. THE ELEMENTS IN THE NORMAL MATRIX
C   (2X2) WHICH IS DIAGONAL USING THE MEAN LONGITUDE AS REFERENCE, AND
C   THE RIGHT HAND SIDE ARE STORED (BINARY) ON UNIT=IUOUT FOR LATER USE.
C
C   PROGRAMMED BY         PER KNUDSEN
C                         KORT- OG MATRIKELSTYRELSEN
C                         RENTEMESTERVEJ 8
C                         DK-2400 KOEBENHAVN NV             15.10.90

C   THIS VERSION OCT 15, 1990.      PER KNUDSEN.
C
C   *****************************************************************
C
C   CONNECTED FILES      IUOUT      EQUATION SYSTEM FOR EACH TRACK
C                                    (TO SUBROUTINE ABIAS)
C
C   CONTENT OF ARRAYS    KREV(.)    REVOLUTION NUMBER
C                         LAT(.)     LATITUDE        (*1000000.)
C                         LON(.)     LONGITUDE       (*1000000.)
C                         KALT(.)    GEOID HEIGTH    (* 1000.)
C                         KERR(.)    STD. DEVIATION  (* 1000.)
C                         LNS(.)     = .TRUE.  SOUTH-GOING TRACK
C                                    = .FALSE. NORTH-GOING TRACK.
C
C   MAX NUMBER OF DATA    100000
C
C   *****************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
      parameter(maxdat=100000,MAXREV=100000,irest=63250)
      integer*2 kerr,kcode
      common/field/ krev(maxdat),lat(maxdat),lon(maxdat),
     .              kalt(maxdat),kerr(maxdat),kcode(maxdat),
     .              i4rest(irest),IFIRST(MAXREV),ILAST(MAXREV)
      LOGICAL LSAVE
C
C   ********************************************************************
C
      WRITE(6,490)
  490 FORMAT(   '1SUBROUTINE ALINE, VERSION 15.10.90,',/)
C
      SUM=0.0D0
      SUM2=0.0D0
      DO 10 I=1,NDAT
      DH=KALT(I)/1.D3
      SUM = SUM + DH*1.0D0
      SUM2= SUM2+(DH*1.0D0)**2
   10 CONTINUE
      IF(NDAT.GT.1) SUM2=DSQRT((SUM2-(SUM**2)/NDAT)/(NDAT-1))
      SUM=SUM/NDAT
C
C***  EVALUATION OF NUMBER OF TRACKS AND BEGINNING AND END OF EACH
C     TRACK.
C
      IREV=1
      KREVCC=KREV(1)
      IFIRST(IREV)=1
      DO 40 I=2,NDAT
      IF(KREV(I).EQ.KREVCC) GO TO 40
      ILAST(IREV)=I-1
      IREV=IREV+1
      IFIRST(IREV)=I
      KREVCC=KREV(I)
   40 CONTINUE
      ILAST(IREV)=NDAT
C
      WRITE(6,510) NDAT,IREV,SUM,SUM2
  510 FORMAT(' ENTERED WITH ',I6,' DATA FROM ',I4,' TRACKS.',
     +//,    ' MEAN VALUE OF (OBS-REF)',F8.3,' M',
     +/,     ' STD.DEV.   -      -    ',F8.3,' M.')
C
      WRITE(*,520)
  520 FORMAT(///,
     .'  I   REV   NUMBER     MEAN.LONG    BIAS   TILT     RMS      ',
     .' RMS',/)
C
      SUM=0.0D0
      SUM2=0.0D0
C
C** FOR EACH TRACH DO
C
      DO 200 I=1,IREV
      NDATTR=ILAST(I)-IFIRST(I)+1
      IF(NDATTR.LT.2) GO TO 200
C
C  FIND MEAN LONGITUDE
C
      AMLON=0.0D0
      ERR2=0.0
      DO 100 K=IFIRST(I),ILAST(I)
      ALON=LON(K)/1.0D6
      ERR=KERR(K)/1.0D3+0.0001
      AMLON=AMLON+ALON
      ERR2=ERR2+ERR**2
  100 CONTINUE
      AMLON=AMLON/NDATTR
      ERR2=ERR2/NDATTR
C
C  SET UP NORMAL EQUATION SYSTEM FOR BIAS/TILT ESTIMATION
C  (THE NORMAL MATRIX IS DIAGONAL WHEN THE MEAN LONGITUDE IS
C  USED AS REFERENCE)
C
      DIAG1=NDATTR*1.0D0
      DIAG2=0.0D0
      RIGHT1=0.0D0
      RIGHT2=0.0D0
      DO 110 K=IFIRST(I),ILAST(I)
      DLON=LON(K)/1.0D6-AMLON
      SSH=KALT(K)/1.0D3
      DIAG2=DIAG2+DLON**2
      RIGHT1=RIGHT1+SSH
      RIGHT2=RIGHT2+DLON*SSH
  110 CONTINUE
      BIAS=RIGHT1/DIAG1
      TILT=RIGHT2/DIAG2
C
C REDUCE SEA SURFACE HEIGHT AND SAVE IF LSAVE=TRUE
C
      S=0.0D0
      S2=0.0D0
      DO 120 K=IFIRST(I),ILAST(I)
      DLON=LON(K)/1.0D6-AMLON
      SSH=KALT(K)/1.0D3
      HBT=BIAS+DLON*TILT
      S=S+SSH**2
      SSH=SSH-HBT
      S2=S2+SSH**2
      IF(LSAVE) KALT(K)=SSH*1000.0+0.5
      SUM=SUM+SSH
      SUM2=SUM2+SSH**2
  120 CONTINUE
      S=DSQRT(S/NDATTR)
      S2=DSQRT(S2/NDATTR)
      DIAG1=DIAG1/ERR2
      DIAG2=DIAG2/ERR2
      RIGHT1=RIGHT1/ERR2
      RIGHT2=RIGHT2/ERR2
      WRITE(IUOUT) KREV(IFIRST(I)),AMLON,DIAG1,DIAG2,RIGHT1,RIGHT2
      WRITE(*,530) I,KREV(IFIRST(I)),NDATTR,AMLON,BIAS,TILT,S,S2
  530 FORMAT(I5,2I8,F10.6,4F8.3)
C
  200 CONTINUE
      IF(NDAT.GT.1) SUM2=DSQRT((SUM2-(SUM**2)/NDAT)/(NDAT-1))
      SUM=SUM/NDAT
      WRITE(6,540) SUM,SUM2
  540 FORMAT(//,' AFTER REMOVAL OF TRACK BIAS/TILT',
     +/,        ' MEAN VALUE OF (OBS-REF)',F8.3,' M',
     +/,        ' STD.DEV.   -      -    ',F8.3,' M.')
c
      REWIND IUOUT
C
      RETURN
      END
C *******************************************************************
C *****                                                         *****
C *****          SUBROUTINE ACRS                                *****
C *****                                                         *****
C *******************************************************************
C
      SUBROUTINE ACRS(IUOUT,RLAT1,RLAT2,RLON1,RLON2,
     *                NDAT,DIST0,DSTMAX)
C   *****************************************************************
C
C   SUBROUTINE ACRS AND PROGRAM MODULES CROSS, AND DIST
C   ARE A MODIFICATION OF PROGRAM ALTCROSS (GEOSOFT).
C
C   PROGRAMS ARE COPYRIGHT BY THE AUTHOR AND THE GEODETIC
C   INSTITUTE OF DENMARK, NOW KORT- OG MATRIKELSTYRELSEN.
C
C   SUBROUTINE ACRS COMPUTES THE LOCATIONS OF CROSS-OVER POINTS
C   AND DETERMINES THE DIFFERENCES IN HEIGHT BETWEEN THE ARCS
C   AS THE NORTH-GOING MINUS THE SOUTH-GOING ARC. THIS IS DONE
C   USING A LINEAR INTERPOLATION BETWEEN THE FOUR NEIGHBOURING
C   DATA POINTS.
C
C   PROGRAMMED BY         PER KNUDSEN
C                         GEODETIC INSTITUTE
C                         GAMLEHAVE ALLE 22
C                         DK-2920 CHARLOTTENLUND           27.06.85
C
C   THIS VERSION MAY 17, 1991.      PER KNUDSEN.
C
C   *****************************************************************
C
C   CONNECTED FILES      IUOUT      CROSS-OVER-DATA
C                                    (TO SUBROUTINE ABIAS)
C
C   CONTENT OF ARRAYS    KREV(.)    REVOLUTION NUMBER
C                         LAT(.)     LATITUDE        (*1000000.)
C                         LON(.)     LONGITUDE       (*1000000.)
C                         KALT(.)    GEOID HEIGTH    (* 1000.)
C                         KERR(.)    STD. DEVIATION  (* 1000.)
C                         LNS(.)     = .TRUE.  SOUTH-GOING TRACK
C                                    = .FALSE. NORTH-GOING TRACK.
C
C   MAX NUMBER OF DATA    100000
C
C   *****************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
      parameter(maxdat=100000,MAXREV=100000,irest=63250)
      integer*2 kerr,kcode
      common/field/ krev(maxdat),lat(maxdat),lon(maxdat),
     .              kalt(maxdat),kerr(maxdat),kcode(maxdat),
     .              i4rest(irest),IFIRST(MAXREV),ILAST(MAXREV)
      INTEGER*4 P1,P2,PA,PB
      LOGICAL LNS(maxdat),LCHECK,LNOUT,LNPEAK,LNTWO,LGEO3N,LGEO3S
C
C   ********************************************************************
C
C***  INPUT OF RLAT1,RLAT2    RANGE IN LATITUDE.
C***           RLON1,RLON2      -    - LONGITUDE.
C***           DIST0,DSTMAX   NORMAL AND MAX DISTANCE IN KM.
C***                           BETWEEN SAMPLES. DEFAULTS OF
C***                           DIST0=6.7 AND DSTMAX=2.5*DIST0
C***                           ARE SET IF 0.0 ARE GIVEN. IF
C***                           DIST0=-1 THE AVERAGE IS USED.
C
      WRITE(6,490)
  490 FORMAT(   '1SUBROUTINE ACRS, VERSION 17.5.91,',/)
C
      SUM=0.0D0
      SUM2=0.0D0
      DO 10 I=1,NDAT
      DH=KALT(I)/1.D3
      SUM = SUM + DH*1.0D0
      SUM2= SUM2+(DH*1.0D0)**2
   10 CONTINUE
      IF(NDAT.GT.1) SUM2=DSQRT((SUM2-(SUM**2)/NDAT)/(NDAT-1))
      SUM=SUM/NDAT
C
C***  EVALUATION OF NUMBER OF TRACKS AND DIRECTION; NORTH OR
C***  SOUTH-GOING (LNS=.FALSE. OR .TRUE.).
C
      IREV=0
      I=1
   39 continue
      IREV=IREV+1
      if (irev.gt.maxrev) stop '*** too many revs, increase maxrev'
      IFIRST(IREV)=I
      KREVCC=KREV(I)
      LNS(i)=(LAT(I).GT.LAT(I+1))
      DO 40 I=(ifirst(irev)+1),NDAT
      LNS(i)=(LAT(i-1).GT.LAT(i))
      IF((KREV(I).EQ.KREVCC).and.(lns(i).eqv.lns(i-1))) GO TO 40
      ILAST(IREV)=I-1
      go to 39
   40 CONTINUE
      ILAST(IREV)=NDAT
C
c     DO 42 I=1,IREV
c     IF(IFIRST(I).EQ.ILAST(I)) GO TO 42 
c     LNS(ifirst(i))=(LAT(IFIRST(I)).GT.LAT(ifirst(I)+1))
c     DO 41 J=(IFIRST(I)+1),ILAST(I)
c     LNS(J)=(LAT(j-1).GT.LAT(j))
c  41 CONTINUE
c  42 CONTINUE
C
      WRITE(6,510) NDAT,IREV,SUM,SUM2
  510 FORMAT(' ENTERED WITH ',I6,' DATA FROM ',I4,' TRACKS.',
     +//,    ' MEAN VALUE OF (OBS-REF)',F8.3,' M',
     +/,     ' STD.DEV.   -      -    ',F8.3,' M.')
C
C***  IF DIST0=-1, THEN DIST0 = MEAN SAMPLE DISTANCE IN KM.
C***    ELSE IF DIST0=0, THEN DIST0 = 6.7 KM.
C***    ELSE THE READ VALUE IS USED.
C***  IF DSTMAX=0 THEN DSTMAX = 2.5*DIST0,
C***    ELSE THE READ VALUE IS USED.
C
      IF(DABS(DIST0).LT.1E-6) DIST0=6.7
      IF(DIST0.GT.0.0) GO TO 45
      DIST0=0.0D0
      NONOR=0
      DO 44 I=2,NDAT
      IF(KREV(I-1).NE.KREV(I)) GO TO 44
      DIST0=DIST0+DIST(LAT(I-1),LON(I-1),LAT(I),LON(I))
      NONOR=NONOR+1
   44 CONTINUE
      DIST0=DIST0/NONOR
   45 IF(DSTMAX.LT.1E-6) DSTMAX=2.5*DIST0
      WRITE(6,550) DSTMAX,DIST0
  550 FORMAT(//
     +,      ' CROSS-OVERS ARE DETERMINED IF THE DISTANCE BETWEEN'
     +,/,    ' THE DATA IS LESS THAN',F8.3,' KM. THE ERROR OF'
     +,/,    ' THE VALUE IS MULTIPLIED WITH THE DISTANCES RELA-'
     +,/,    ' TIVE TO AN AVERAGE SAMPLE DISTANCE OF',F7.3,' KM.',//)
C
C***  EACH NORTH-GOING TRACK IS COMPARED WITH EVERY SOUTH-
C***  GOING TRACK AND THE DIFFERENCE IN HEIGTH IS CALCULATED
C***  IF THE CROSS-POINT IS IN THE AREA.
C
      WRITE(6,650)
  650 FORMAT(//,'          REVOLUTION     NUMBER OF',/
     +         ,'            NUMBER      CROSSOVERS'/)
C
      NONOR=0
      NOSOU=0
      IALT =0
      INOR =0
      S1=0.0D0
      S2=0.0D0
      S3=0.0D0
C
      DO 201 IKKK=1,IREV
      IF((ILAST(IKKK)-IFIRST(IKKK)).LT.2) GO TO 201
      DO 200 I=IFIRST(IKKK),ILAST(IKKK)
      IF(LNS(I).OR.KREV(I).EQ.NONOR) GO TO 201
      NONOR=KREV(I)
      INOR=INOR+1
      ISOU=0
C
      DO 101 JJKKK=1,IREV
      IF((ILAST(JJKKK)-IFIRST(JJKKK)).LT.2) GO TO 101
      DO 100 JJ=IFIRST(JJKKK),ILAST(JJKKK)
      J=JJ
      IF(.NOT.LNS(J).OR.KREV(J).EQ.NOSOU) GO TO 101
      NOSOU=KREV(J)
C
      DET=CROSSC(LAT(I),LON(I),LAT(ILAST(IKKK)),LON(ILAST(IKKK)),
     +LAT(J),LON(J),LAT(ILAST(JJKKK)),LON(ILAST(JJKKK)),AL,BE,P1,P2)
C
      IF(DET.EQ.0.0) GO TO 101
      IF(AL.LT.0.0.OR.AL.GT.1.0.OR.BE.LT.0.0.OR.BE.GT.1.0) GO TO 101
      AL=AL*(ILAST(IKKK)-IFIRST(IKKK)+1.0)
      BE=BE*(ILAST(JJKKK)-IFIRST(JJKKK)+1.0)
      K=I
C
C**   THE FOUR NEIGHBOURING DATAPOINTS OF THE CROSS-OVER POINT
C**   ARE FOUND IN MAXIMUM 5 ITERATIONS.
C
      DO 80 NIT=1,5
C
      IA=AL
      IB=BE
      IIA=0
      IIB=0
C
      IF(IA.EQ.0) GO TO 55
      IAD=IA/IABS(IA)
      I1=(1-IAD)/2
      I2=(1+IAD)/2
      IF(LAT(K+I2+IAD).EQ.P1) GO TO 55
      DO 50 IIA=IAD,IA,IAD
      IF((K+IAD+IIA).LT.1.OR.(K+IAD+IIA).GT.NDAT) GO TO 101
      IF(LAT(K+I2+IIA).EQ.P1) GO TO 55
      PA=(LAT(K+I1+IIA)-P1)/(LAT(K+I2+IIA)-P1)
      IF(PA.LE.0) GO TO 55
   50 CONTINUE
   55 CONTINUE
      IF(IB.EQ.0) GO TO 65
      IBD=IB/IABS(IB)
      I1=(1-IBD)/2
      I2=(1+IBD)/2
      IF(LAT(J+I2+IBD).EQ.P1) GO TO 65
      DO 60 IIB=IBD,IB,IBD
      IF((J+IBD+IIB).LT.1.OR.(J+IBD+IIB).GT.NDAT) GO TO 101
      IF(LAT(J+I2+IIB).EQ.P1) GO TO 65
      PB=(LAT(J+I1+IIB)-P1)/(LAT(J+I2+IIB)-P1)
      IF(PB.LE.0) GO TO 65
   60 CONTINUE
   65 CONTINUE
      IA=IIA
      IB=IIB
C
      IF((K+IA+1).GT.NDAT.OR.(J+IB+1).GT.NDAT) GO TO 101
      IF(KREV(K+IA).NE.NONOR.OR.KREV(K+IA+1).NE.NONOR) GO TO 101
      IF(KREV(J+IB).NE.NOSOU.OR.KREV(J+IB+1).NE.NOSOU) GO TO 101
      K=K+IA
      J=J+IB
C
      DET=CROSS(LAT(K),LON(K),LAT(K+1),LON(K+1),
     +LAT(J),LON(J),LAT(J+1),LON(J+1),AL,BE,P1,P2)
C
      IF(DET.EQ.0.0) GO TO 101
      IF(AL.LT.0.)AL=AL-1.
      IF(BE.LT.0.)BE=BE-1.
C
      IF((AL.GE.0..AND.AL.LE.1.).AND.(BE.GE.0..AND.BE.LE.1.)) GO TO 90
C
   80 CONTINUE
      GO TO 101
C
   90 CONTINUE
C
C
C**   THE POSITION AND HEIGHT-DIFFERENCE ARE INTERPOLATED AND
C**   THE RESULTS ARE WRITTEN ON UNIT IUOUT.
C
      ALTK = KALT(K)  /1000.
      ALTK1= KALT(K+1)/1000.
      ALTJ = KALT(J)  /1000.
      ALTJ1= KALT(J+1)/1000.
      DH=(ALTK+AL*(ALTK1-ALTK))-(ALTJ+BE*(ALTJ1-ALTJ))
C
      ERRK = KERR(K)  /1000.
      ERRK1= KERR(K+1)/1000.
      ERRJ = KERR(J)  /1000.
      ERRJ1= KERR(J+1)/1000.
      ERR=(1.-AL)*ERRK**2+AL*ERRK1**2+(1.-BE)*ERRJ**2+BE*ERRJ1**2
      ERR=DSQRT(ERR)
C
      DISTK=DIST(LAT(K),LON(K),LAT(K+1),LON(K+1))
      DISTJ=DIST(LAT(J),LON(J),LAT(J+1),LON(J+1))
      IF(DISTK.GT.DSTMAX.OR.DISTJ.GT.DSTMAX) GO TO 100
      if(Dist0.ne.0.0d0) ERR=ERR*DISTK*DISTJ/(DIST0**2)
C
      ISOU=ISOU+1
      IALT=IALT+1
C     NCRS=NONOR
C     NCRS=100000*NCRS+NOSOU
C
      AP1=P1/1000000.
      AP2=P2/1000000.
      WRITE(IUOUT,700) NONOR,NOSOU,AP1,AP2,DH,ERR
  700 FORMAT(I7,I5,2(F10.4,' '),F8.2,13X,F5.2,10X)
      S1=S1+DH
      S2=S2+DH**2
      S3=S3+(DH/ERR)**2
C
  100 CONTINUE
  101 CONTINUE
C
      WRITE(6,710) INOR,NONOR,ISOU
  710 FORMAT(1X,I5,2I12)
      GO TO 201
C
  200 CONTINUE
  201 CONTINUE
C
      IF(IALT.GT.0) S1=S1/IALT
      IF(IALT.GT.0) S2=DSQRT(S2/IALT)
      IF(IALT.GT.0) S3=DSQRT(S3/IALT)
      WRITE(6,720) IUOUT,IALT,S1,S2,S3
  720 FORMAT(//,' THE NUMBER OF CROSSOVERS',/
     +,         ' WRITTEN ON FILE UNIT',I3,' IS',I6,//
     +,         ' MEAN VALUE         ',F7.3,/
     +,         ' RMS VALUE          ',F7.3,/
     +,         ' RMS OF (DH/ER)      ',F7.3,/)
C
      REWIND IABS(IUOUT)
C
      RETURN
      END
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
      FUNCTION CROSS(I1,J1,I2,J2,I3,J3,I4,J4,ALFA,BETA,P1,P2)
C
C    THIS FUNCTION CALCULATES THE CROSS-OVER POSITION (P1,P2),
C    USING A PLANAR APPROXIMATION, AND RETURNS WITH THE DETERMINANT.
C
C    PROGRAMMED BY PER KNUDSEN.
C
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
C
      INTEGER*4 P1,P2
C
      X1=I1
      X2=I2
      X3=I3
      X4=I4
      Y1=J1
      Y2=J2
      Y3=J3
      Y4=J4
C
      A=X2-X1
      B=X3-X4
      C=Y2-Y1
      D=Y3-Y4
      E=X3-X1
      F=Y3-Y1
C
C***  THE FUNCTION SOLVES     ( A  B ) ( ALFA )   ( E )
C***                           (      )*(      ) = (   )
C***                           ( C  D ) ( BETA )   ( F )
C***
C***  AND RETURNS WITH THE DETERMINANT.
C
      DET=(A*D-B*C)
      CROSS=DET
      IF(DET.lt.1.0d-17) RETURN
      ALFA=(D*E-B*F)/DET
      BETA=(A*F-C*E)/DET
C
      AP1=X1+ALFA*A
      AP2=Y1+ALFA*C
C
      P1=(AP1+.5)
      P2=(AP2+.5)
C
      RETURN
      END
C
C   *****************************************************************
C
      FUNCTION DIST(I1,J1,I2,J2)
C
C    PROGRAMMED BY PER KNUDSEN.
C
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
C
C***  FUNCTION DIST CALCULATES THE DISTANCE BETWEEN TWO POINTS
C***  IN A PLANAR APPROXIMATION GIVEN IN KILOMETERS. IN THE CON-
C***  VERSION FROM RAD. TO KM. A MEAN EARTH RADIUS OF 6371 KM IS
C***  USED.
C***  THE COORDINATES (I1,J1) AND (I2,J2) ARE INTEGERS AND EQUAL
C***  TO THE COORDINATES IN DEGREES*1000000.
C
      DLAT=(I1-I2)/1000000.0
      DLON=(J1-J2)/1000000.0
      RLAT=I1/180.0/1000000.0*3.141592
      DLON=DLON*DCOS(RLAT)
      DIST=DSQRT(DLAT**2+DLON**2)/180.0*3.141592*6371.0
      RETURN
      END
C *******************************************************************
C *****                                                         *****
C *****          SUBROUTINE ABIAS                               *****
C *****                                                         *****
C *******************************************************************
C
      SUBROUTINE ABIAS(IUIN,IUINC,IUOUT,HST,WEIGHT,LTILT,LNOERR)
C   *****************************************************************
C
C   SUBROUTINE ABIAS AND THE PROGRAM MODULES INDI, INDJ,
C   SORTD, AND PROCNL ARE A MODIFICATION OF PROGRAM
C   ALTCOR (GEOSOFT).
C
C   PROGRAMS ARE COPYRIGHT BY THE AUTHOR AND THE GEODETIC
C   INSTITUTE OF DENMARK.
C
C   SUBROUTINE ABIAS COMPUTES A NUMBER OF UNCORRELATED BIAS
C   PARAMETERS FROM THE CROSS-OVER DIFFERENCES PRODUCES BY
C   SUBROUTINE ACRS. ONE PARAMETER IS ASSOCIATED WITH EACH ARC
C   AND THE CONSTRAINT THAT IS NEEDED, IS THAT THE SUM OF THE
C   BIAS PARAMETERS SHALL BE ZERO.
C
C   PROGRAMMED BY         PER KNUDSEN
C                         GEODETIC INSTITUTE
C                         GAMLEHAVE ALLE 22
C                         DK-2920 CHARLOTTENLUND           27.02.85.
C
C   IN THIS VERSION TILT ESTIMATION HAS BEEN INCLUDED (OPTIONAL).
C   THEN 3 MORE CONSTRAINTS ARE NEEDED, BUT 3 CONSTRAINTS SPANNING
C   THIS NULL-SPACE ONLY IS NOT EASYLY FOUND WITHOUT A COORDINATE-
C   TRANSFORMATION (SO TRACKS BECOME PARALLEL). INSTEAD A LITTLE
C   SOMETHING IS ADDED TO DIAGONAL (SUPPRISED?) ASSUMING THAT NO
C   OBSERVABLE HARM IS DONE TO THE SOLUTION. HOWEVER THIS ASSUMPTION
C   IS REASONABLE SINCE ALL EIGENVALUES SPANING ARE LARGE EXCEPT FOR
C   THE LAST FOUR, WHICH ARE VERY SMALL (SEE SCHRAMA'S DISSERTATION).
C   FURTHERMORE ABSOLUTE HEIGHTS (E.G. ABOVE THE GEOID) THROUGH
C   SUBROUTINE ALINE CAN BE USED TOGETHER WITH THE CROSS-OVERS IN
C   ORDER TO FORCE THE ALTIMETRIC SURFACE TO THE GEOID. HENCE THE
C   THE CONSTRAINTS MENTIONED ABOVE ARE NOT USED.
C   IN BOTH CASES SHORT ARCS WITH ONE CROSS POINT ONLY, THEIR TILTS
C   ARE FIXED TO ZERO.
C
C   PROGRAMMED BY         PER KNUDSEN
C                         KORT- OG MATRIKELSTYRELSEN
C                         RENTEMESTERVEJ 8
C                         DK-2400 KBH. NV                    23.10.90
C
C   *****************************************************************
C
C   CONNECTED FILES       IUIN         CROSS-OVER-DATA (FROM ACRS)
C                         IUINC        NORMAL MATRIX ELEMENTS OF
C                                      ABSOLUTE HEIGHTS (FROM ALINE)
C                         IUOUT        RESULT TRACKNO.,BIAS,STD.DEV.,
C                                              MEAN LONGITUDE,TILT,STD.DEV.
C
C   CONTENT OF ARRAYES    NN(.)        NORTH-GOING TRACKNUMBER
C                         NS(.)        SOUTH-GOING      -
C                         ALON(.)      LONGITUDE
C                         HNS(.)       CROSS-OVER DIFFERENCE
C                         ERRY(.)      STD.DEVIATION OF HNS
C                         INOR(.)      LIST OF NORTH-GOING TRACKS
C                         ISOU(.)       -   -  SOUTH-  -     -
C                         AMLON(.)     MEAN LONGITUDE FOR EACH TRACK
C                         IGROUP(.)    SUB GROUP NUMBER
C                         AST(.,3)     USED FOR STATISTICS
C                         NHIST(41,3)   -    -      -
C                         A(.)         UPPER PART OF (AT)*(C-1)*(A)
C                         INUL(.)      INDEX OF FIRST NONZERO ELEMENT
C                         H(.)         (AT)*(C-1)*Y AND RESULT
C                         ERRX(.)      STD.DEVIATION OF RESULT
C                         E(.)         SCRATCH ARRAY
C
C   *****************************************************************
C
C   INPUT ARE     IUIN,IUINC,IUOUT
C                           IF IUINC=0 THEN NO HEIGHTS AND
C                           CONSEQUENTLY USE OF THE ZERO-MEAN-
C                           VALUE CONSTRAINT
C                 HST      INTERVAL IN HISTOGRAMS.
C                 WEIGHT   WEIGHT OF HEIGHT INFORMATIONS
C                           OR QUANTITY TO BE ADDED TO DIAG.
C                 LTILT    TRUE IF TILT PARAMETERS ARE INCLUDED.
C                 LNOERR   IF TRUE, NO ERRORS ARE CALCULATED.
C
C   *****************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(maxdat=100000,maxpar=1500,iat=381625) 
C     COMMON FIELD OF IAT*8 BYTES (PT. .EQ. 563250*4 BYTES)
      COMMON/FIELD/ A(IAT)
      COMMON/FUNC/ IN,JS,INOR(maxdat),ISOU(maxdat)
      REAL*4 HNS,ERRY
      DIMENSION NN(MAXDAT),NS(MAXDAT),HNS(MAXDAT),ERRY(MAXDAT),
     .          ALON(MAXDAT),
     .          AST(MAXPAR,3),AMLON(MAXPAR),H(MAXPAR),E(MAXPAR),
     .          ERRX(MAXPAR),INUL(MAXPAR),IGROUP(MAXPAR),
     .          NHIST(41,3)
      LOGICAL LRED,LBS,LTILT,LNOERR,lchang
      CHARACTER*1 CTYPE(MAXPAR)
      CHARACTER*8 STRING1,STRING2,STRING3
C
      INULT=MAXPAR
      IHT=MAXPAR
C
      MODPAR=1
      IF(LTILT) MODPAR=2
C
C***   INPUT OF CROSS-OVER DATA FROM FILE UNIT=IUIN.
C
      DO 2 I=1,MAXDAT
      READ(IUIN,510,END=4) NN(I),NS(I),ALON(I),HNS(I),ERRY(I)
    2 CONTINUE
      WRITE(*,*) 'NUMBER OF CROSS-OVERS EXCEEDS MAXDAT'
      STOP
    4 CONTINUE
      NCRS=I-1
  510 FORMAT(2X,2I5,11X,F10.4,1X,F8.2,13X,F5.2)
C
      IN=1
      JS=1
      INOR(1)=NN(1)
      ISOU(1)=NS(1)
      DO 5 I=1,41
      NHIST(I,1)=0
      NHIST(I,2)=0
      NHIST(I,3)=0
    5 CONTINUE
      HH=0.
      HHN=0.
      HM=0.
C
C***  EVALUATION OF THE NUMBER OF NORTH-GOING TRACKS (IN) AND SOUTH-
C***  GOING TRACKS (JS) AND THE TOTAL NUMBER (NT). CALCULATION OF MEAN
C***  VALUE AND RMS VALUE OF DIFFERENCES ( N .VS. S ).
C
      DO 50 IJ=1,NCRS
      HHN=HHN+(HNS(IJ)/ERRY(IJ))**2
      HH=HH+(HNS(IJ)**2)
      HM=HM+HNS(IJ)
      AND=HNS(IJ)/HST
      IND=AND
      IF(AND.LT.0.)IND=IND-1
      IND=IND+21
      IF(IND.GT.41.OR.IND.LT.1)IND=41
      NHIST(IND,3)=NHIST(IND,3)+1
      DO 10 I=1,IN
      IF(NN(IJ).EQ.INOR(I)) GO TO 20
   10 CONTINUE
      IN=IN+1
      INOR(IN)=NN(IJ)
   20 CONTINUE
      DO 30 J=1,JS
      IF(NS(IJ).EQ.ISOU(J)) GO TO 40
   30 CONTINUE
      JS=JS+1
      ISOU(JS)=NS(IJ)
   40 CONTINUE
   50 CONTINUE
      NT=IN+JS
      HHN=DSQRT(HHN/NCRS)
      HH=DSQRT(HH/NCRS)
      HM=HM/NCRS
C
      CALL SORTD
C
      WRITE(6,520) NT,NCRS
      I=MODPAR*NT
      IF(LTILT) WRITE(6,530) I
      IF(I.GT.MAXPAR) STOP '*** i .gt. maxpar'
      I=(I*(I+1))/2
      WRITE(6,540) I
      IF(I.GT.IAT) STOP '*** i .gt. iat'
C
C***  INITIALIZING ARRAYS.
C
      DO 90 I=1,(NT*MODPAR)
      IST=(I*(I-1))/2
      DO 80 J=1,I
      A(IST+J)=0.0D0
   80 CONTINUE
      H(I)=0.0D0
      ERRX(I)=0.0D0
      INUL(I)=1
      AMLON(I)=0.0D0
      E(I)=0.0D0
      AST(I,1)=0.0D0
      AST(I,2)=0.0D0
      AST(I,3)=0.0D0
      CTYPE(I)=' '
   90 CONTINUE
C
  520 FORMAT('1',/
     +,' SUBROUTINE ABIAS, VERSION 23.10.90, DETERMINES THE BIAS-',/
     +,' CORRECTION VALUES FOR THE TRACKS FROM THE CROSS-OVER-',/
     +,' DIFFERENCES.',//
     +,' NUMBER OF TRACKS       ',I6,/
     +,' NUMBER OF CROSSOVERS   ',I6)
  530 FORMAT(/
     +,' INCLUDING TILTS THE NUMBER OF PARAMETERS ARE',I6)
  540 FORMAT(/
     +,' NUMBER OF ELEMENTS USED FOR THE EQUATIONS',I8)
  550 FORMAT(/
     +,' MATRIX ELEMENTS ASSOCIATED WITH ABSOLUTE HEIGHT INFORMATION',/
     +,' IS ENTERED FROM UNIT=',I3,' AND USED WITH A WEIGHT OF',D12.4)
  555 FORMAT(/,' NUMBER OF TRACKS INVOLVED',I6)
  560 FORMAT(/,' PARAMETERS ARE CONSTRAINTED BY ADDING TO DIAG',D10.4)
  565 FORMAT(/,' PARAMETERS ARE CONSTRAINTED USING ZERO-MEAN-VALUE.')
C
C  THE MEAN LONGITUDE FOR EACH TRACK IS COMPUTED
C
      DO 100 IJ=1,NCRS
      I=INDI(NN(IJ))
      J=INDJ(NS(IJ))
      IF(I.GT.0.AND.J.GT.0) GO TO 110
      WRITE(6,600)
  600 FORMAT(/,2(/,20(1H*)),' ERROR IN INDEX FUNCTION. STOP')
      STOP
  110 CONTINUE
      J=IN+J
      E(I)=E(I)+1.0/ERRY(IJ)**2
      AMLON(I)=AMLON(I)+ALON(IJ)/ERRY(IJ)**2
      E(J)=E(J)+1.0/ERRY(IJ)**2
      AMLON(J)=AMLON(J)+ALON(IJ)/ERRY(IJ)**2
  100 CONTINUE
      DO 101 I=1,NT
      AMLON(I)=AMLON(I)/E(I)
  101 CONTINUE
C
C  IF ABSOLUTE HEIGHT INFORMATION IS USED TOGETHER WITH THE CROSS-
C  OVERS TO OBTAIN AN ABSOLUTE CONSTRAINT, NORMAL MATRIX ELEMENTS
C  PRODUCED BY SUBROUTINE ALINE ARE ENTERED TOGETHER WITH MEAN LON-
C  GITUDES USED IN THIS COMPUTATION.
C
      NREVC=0
      IF(IUINC.EQ.0) GO TO 1102
      WRITE(*,550) IUINC,WEIGHT
 1100 READ(IUINC,END=1101) KREVC,AMLONC,ATA1C,ATA2C,ATH1C,ATH2C
      I=INDI(KREVC)
      J=INDJ(KREVC)
      IF(I.EQ.0.AND.J.EQ.0) GO TO 1100
      NREVC=NREVC+1
      IF(I.EQ.0) I=IN+J
      AMLON(I)=AMLONC
      IB=(I-1)*MODPAR+1
      H(IB)=ATH1C*WEIGHT
      IST=(IB*(IB+1))/2
      A(IST)=ATA1C*WEIGHT
      IF(LTILT) THEN
        IT=IB+1
        H(IT)=ATH2C*WEIGHT
        IST=(IT*(IT+1))/2
        A(IST)=ATA2C*WEIGHT
        ENDIF
      GO TO 1100
 1101 CONTINUE
      WRITE(*,555) NREVC
 1102 CONTINUE
C
      IF(IUINC.EQ.0) THEN
        IF(LTILT) THEN
          WRITE(*,560) WEIGHT
        ELSE
          WRITE(*,565)
        ENDIF
      ENDIF
C
C***   THE NORMAL EQUATIONS ARE BUILT UP SUCCESSIVELY WITH THE UPPER
C***   PART OF (AT)*(C-1)*(A) IN ARRAY 'A' AND (AT)*(C-1)*Y IN 'H'.
C***   TRACKS WITH THE SAME NUMBER ARE CONSTRAINTED.
C***   SINCE THE RANK OF SYSTEM IS (NT-1) A CONSTRAINT, WHICH ENSURES
C***   THAT THE AVERAGE OF THE CORRECTIONS IS ZERO, IS ADDED (IF IUINC=0).
C***   STATISTICAL PARAMETERS ARE DETERMINED.
C
      DO 120 IJ=1,NCRS
      I=INDI(NN(IJ))
      J=INDJ(NS(IJ))+IN
      IB=(I-1)*MODPAR+1
      JB=(J-1)*MODPAR+1
      H(IB)=H(IB)+HNS(IJ)/(ERRY(IJ)**2)
      H(JB)=H(JB)-HNS(IJ)/(ERRY(IJ)**2)
      IST=(IB*(IB+1))/2
      A(IST)=A(IST)+1.0D0/(ERRY(IJ)**2)
      IST=(JB*(JB+1))/2
      A(IST)=A(IST)+1.0D0/(ERRY(IJ)**2)
      IST=(JB*(JB-1))/2+IB
      A(IST)=-1.0D0/(ERRY(IJ)**2)
      IF(LTILT) THEN
        IT=IB+1
        JT=JB+1
        DLONN=ALON(IJ)-AMLON(I)
        DLONS=ALON(IJ)-AMLON(J)
        H(IT)=H(IT)+DLONN*HNS(IJ)/(ERRY(IJ)**2)
        H(JT)=H(JT)-DLONS*HNS(IJ)/(ERRY(IJ)**2)
        IST=(IT*(IT+1))/2
        A(IST)=A(IST)+(DLONN/ERRY(IJ))**2
        IST=IST-1
        A(IST)=A(IST)+DLONN/(ERRY(IJ)**2)
        IST=(JT*(JT+1))/2
        A(IST)=A(IST)+(DLONS/ERRY(IJ))**2
        IST=IST-1
        A(IST)=A(IST)+DLONS/(ERRY(IJ)**2)
        IST=(JT*(JT-1))/2+IT
        A(IST)=-(DLONN*DLONS)/(ERRY(IJ)**2)
        IST=IST-1
        A(IST)=-DLONS/(ERRY(IJ)**2)
        IST=(JB*(JB-1))/2+IT
        A(IST)=-DLONN/(ERRY(IJ)**2)
      ENDIF
C
      AST(I,1)=AST(I,1)+1.
      AST(J,1)=AST(J,1)+1.
      AST(I,2)=AST(I,2)+HNS(IJ)
      AST(J,2)=AST(J,2)-HNS(IJ)
      AST(I,3)=AST(I,3)+HNS(IJ)**2
      AST(J,3)=AST(J,3)+HNS(IJ)**2
  120 CONTINUE
C
C  CHECK IF A NORTH AND A SOUTH GOING TRACK HAS THE SAME REV.NO.
C  IN THAT CASE MODEL PARAMETERS ARE PUT EQUAL TO EACH OTHER (TAKING
C  THE DIFFERENCE BETWEEN THE AMLON'S INTO ACCOUNT)
C  BIAS(I)+DLONN*TILT(I)=BIAS(J)+DLONS*TILT(J) AND THEN
C  TILT(I)-TILT(J)=0, BIAS(I)-BIAS(J)+(AMLON(J)-AMLON(I))*TILT(I)=0
C
      J=1
      DO 130 I=1,IN
      IF(INOR(I).LT.ISOU(J).OR.J.EQ.JS) GO TO 130
  126 IF(INOR(I).GT.ISOU(J)) J=J+1
      IF(INOR(I).NE.ISOU(J)) GO TO 128
      CTYPE(I)='='
      CTYPE(J)='='
      IB=(I-1)*MODPAR+1
      IST=(IB*(IB+1))/2
      A(IST)=A(IST)+100000.0D0
      JB=(IN+J-1)*MODPAR+1
      IST=(JB*(JB+1))/2
      A(IST)=A(IST)+100000.0D0
      IST=(JB*(JB-1))/2+IB
      A(IST)=-100000.0D0
      IF(LTILT) THEN
        IT=IB+1
        JT=JB+1
        DLONN=AMLON(J)-AMLON(I)
        IST=IST+1
        A(IST)=-DLONN*100000.0D0
        IST=(IT*(IT+1))/2
        A(IST)=A(IST)+(1.0D0+DLONN**2)*100000.0D0
        IST=IST-1
        A(IST)=A(IST)+DLONN*100000.0D0
        IST=(JT*(JT+1))/2
        A(IST)=A(IST)+100000.0D0
        IST=IST-JT+IT
        A(IST)=-100000.0D0
        ENDIF
  128 IF(INOR(I).GT.ISOU(J).AND.J.LT.JS) GO TO 126
  130 CONTINUE
C
c FIND SEPARATE GROUPS OF TRACKS THAT ARE NOT CONNECTED.
C THIS IS DONE BY CHECKING THE NORMAL EQUATIONS.
c
      do 407 i=1,(nt*MODPAR),MODPAR
  407 igroup(i)=0
      ngroup=0
c
  410 ngroup=ngroup+1
      do 420 i=1,(nt*MODPAR),MODPAR
      if(igroup(i).eq.0) go to 421
  420 continue
  421 istart=i
      igroup(istart)=ngroup
c
  425 do 441 i=istart,(in*MODPAR),MODPAR
      lchang=.false.
      if(igroup(i).ne.ngroup) go to 441
      do 430 jj=1,(js*MODPAR),MODPAR
      j=in*MODPAR+jj
      if(igroup(j).gt.0) go to 430
      ij=(j*(j-1))/2+i
      if(a(ij).gt.-1.0e-6) go to 430
      lchang=.true.
      igroup(j)=ngroup
  430 continue
c
      if(.not.lchang) go to 441
      lchang=.false.
      do 440 jj=1,(js*MODPAR),MODPAR
      j=in*MODPAR+jj
      if(igroup(j).ne.ngroup) go to 440
      ij0=(j*(j-1))/2
      do 435 ii=istart,(in*MODPAR),MODPAR
      if(igroup(ii).gt.0) go to 435
      ij=ij0+ii
      if(a(ij).gt.-1.0e-6) go to 435
      lchang=.true.
      igroup(ii)=ngroup
  435 continue
  440 continue
      if(lchang) go to 425
  441 continue
c
      do 450 i=1,(nt*MODPAR),MODPAR
      IGROUP(I+1)=IGROUP(I)
      if(igroup(i).eq.0) go to 410
  450 continue
c
      write(6,999) ngroup
  999 format('0Number of separate groups',i4)
c
C ADD ZERO-MEAN-VALUE CONSTRAINT IN EACH GROUP (IF IUINC=0)
C
      IF(IUINC.EQ.0.AND.(.NOT.LTILT)) THEN
        DO 455 K=1,NGROUP
          do 455 j=1,(nt*MODPAR),MODPAR
            IF(IGROUP(J).NE.K) GO TO 455
            ij0=(j*(j-1))/2
            do 454 i=1,j,MODPAR
              if(igroup(I).ne.k) go to 454
              ij=ij0+i
              a(ij)=a(ij)+1.0
  454       CONTINUE
  455     CONTINUE
        ENDIF
C
      IF(LTILT) THEN
C
C  FIX TILTS TO ZERO FOR SHORT ARCS (FOUND BY EVALUATING SUB-DETERMINANT).
C
        ELIM=1.0
        DO 460 IT=2,(NT*MODPAR),MODPAR
        IST=(IT*(IT+1))/2
        AMS=A(IST)-A(IST-1)*2.0D0/A(IST-IT)
        IF(AMS.LT.ELIM) THEN
          A(IST)=A(IST)+1000.0
          CTYPE((IT-1)/2+1)='-'
        ENDIF
  460   CONTINUE
C
C  ADD A LITTLE SOMETHING (THE PARAMETER WEIGHT) TO THE DIAGONAL
C  (IF IUINC=0)
C
        IF(IUINC.EQ.0) THEN
          DO 470 I=1,(NT*MODPAR)
            IST=(I*(I+1))/2
            A(IST)=A(IST)+WEIGHT
  470     CONTINUE
        ENDIF
      ENDIF
C
C*** THE INDEX OF THE FIRST NON-ZERO ELEMENT IN EACH COLUMN IS FOUND.
C
      DO 146 IS=1,(NT*MODPAR)
      IST=(IS*(IS-1))/2
      DO 144 IR=1,IS
      INUL(IS)=IR
      IF(A(IST+IR).NE.0.0D0) GO TO 146
  144 CONTINUE
  146 CONTINUE
C
C***  WRITE OUT CONDITIONS OF THE SYSTEM
C
      WRITE(6,610) IN,JS
      DO 150 I=1,IN
      IGRP=IGROUP((I-1)*MODPAR+1)
      H1=AST(I,2)/AST(I,1)
      H2=DSQRT(AST(I,3)/AST(I,1))
      WRITE(6,620) I,INOR(I),CTYPE(I),IGRP,AMLON(I),AST(I,1),H1,H2
      AND=H1/HST
      IND=AND
      IF(AND.LT.0.) IND=IND-1
      IND=IND+21
      IF(IND.GT.41.OR.IND.LT.1) IND=41
      NHIST(IND,1)=NHIST(IND,1)+1
      AND=H2/HST
      IND=AND
      IF(AND.LT.0.)IND=IND-1
      IND=IND+21
      IF(IND.GT.41.OR.IND.LT.1) IND=41
      NHIST(IND,2)=NHIST(IND,2)+1
  150 CONTINUE
C
      DO 160 J=1,JS
      H1=AST(IN+J,2)/AST(IN+J,1)
      H2=DSQRT(AST(IN+J,3)/AST(IN+J,1))
      IGRP=IGROUP((IN+J-1)*MODPAR+1)
      WRITE(6,622) J,ISOU(J),CTYPE(IN+J),IGRP,AMLON(IN+J),AST(IN+J,1),
     .             H1,H2
      AND=H1/HST
      IND=AND
      IF(AND.LT.0.) IND=IND-1
      IND=IND+21
      IF(IND.GT.41.OR.IND.LT.1) IND=41
      NHIST(IND,1)=NHIST(IND,1)+1
      AND=H2/HST
      IND=AND
      IF(AND.LT.0.) IND=IND-1
      IND=IND+21
      IF(IND.GT.41.OR.IND.LT.1) IND=41
      NHIST(IND,2)=NHIST(IND,2)+1
  160 CONTINUE
C
  610 FORMAT(/,' NUMBER OF CROSSOVERS WITH MEAN AND RMS VALUE FOR'
     +,/,      ' THE ',I5,' NORTH-GOING TRACKS, AND'
     +,/,      ' THE ',I5,' SOUTH-GOING TRACKS. (IN METERS)'
     +,/,      ' (REV.S MARKED = ARE BOTH NORTH AND SOUTH-GOING)'
     +,/,      ' (REV.S MARKED - WILL HAVE TILT=0)'
     +//,7X
     +,'   REV.NO.  GROUP   MLON    NO.CRS.      MEAN       RMS',/)
  620 FORMAT(1X,I5,'. N',I7,A1,I6,F9.2,F8.0,4X,2F10.3)
  622 FORMAT(1X,I5,'. S',I7,A1,I6,F9.2,F8.0,4X,2F10.3)
  630 FORMAT(/,' THE MEAN VALUE OF ALL CROSS-OVER-DIF.S',F10.3,/
     +,        '  -   RMSVALUE               -         ',F11.3,/
     +,        '  -    -    -   OF CROSS-OVERS/STD.DEV.',F11.3)
C
      WRITE(6,640) HST
      DO 220 J=1,3
      IF(J.NE.2) WRITE(6,645)(NHIST(I,J),I=1,20)
      IF(J.NE.2) WRITE(6,647)
      WRITE(6,645)(NHIST(I,J),I=21,40)
      WRITE(6,645)(I,I=1,20)
      WRITE(6,649)NHIST(41,J)
  220 CONTINUE
C
      WRITE(6,630)HM,HH,HHN
C
  640 FORMAT(/,' DISTRIBUTION WITH INTERVALLENGTH ',F5.2,' M'
     +,/,      ' OF          MEAN DIFFERENCES PER TRACK'
     +,/,      '              RMS        -      -    -  '
     +,/,      '              ALL CROSS-OVER-DIFFERENCES.',/)
  645 FORMAT(1X,25I4)
  647 FORMAT('  -19 -18 -17 -16 -15 -14 -',
     +'13 -12 -11 -10  -9  -8  -7  -6  -5  -4  -3  -2  -1   0')
  649 FORMAT(68X,'OUTSIDE',I4,/)
C
C***  SOLVE THE (NT*MODPAR) LINEAR EQUATIONS WITH CALL OF PROCNL.
C
      LRED=.TRUE.
      LBS=.TRUE.
      CALL PROCNL(A,INUL,H,(NT*MODPAR),E(1),LRED,LBS,IAT,INULT,IHT)
C
      DO 250 I=1,NT
      AST(I,1)=0.
      AST(I,2)=0.
      AST(I,3)=0.
  250 CONTINUE
      HM=0.
      HH=0.
      HHN=0.
      DO 252 I=1,41
      NHIST(I,1)=0
      NHIST(I,2)=0
      NHIST(I,3)=0
  252 CONTINUE
C
C***   CALCULATION OF NEW STATISTICAL PARAMETERS
C
      DO 280 IJ=1,NCRS
      I=INDI(NN(IJ))
      J=INDJ(NS(IJ))+IN
      IB=(I-1)*MODPAR+1
      JB=(J-1)*MODPAR+1
      HNS(IJ)=HNS(IJ)-H(IB)+H(JB)
      IF(LTILT) THEN
        IT=IB+1
        JT=JB+1
        DLONN=ALON(IJ)-AMLON(I)
        DLONS=ALON(IJ)-AMLON(J)
        HNS(IJ)=HNS(IJ)-DLONN*H(IT)+DLONS*H(JT)
        ENDIF
      AST(I,1)=AST(I,1)+1.
      AST(J,1)=AST(J,1)+1.
      AST(I,2)=AST(I,2)+HNS(IJ)
      AST(J,2)=AST(J,2)-HNS(IJ)
      AST(I,3)=AST(I,3)+HNS(IJ)**2
      AST(J,3)=AST(J,3)+HNS(IJ)**2
      HM=HM+HNS(IJ)
      HH=HH+HNS(IJ)**2
      HHN=HHN+(HNS(IJ)/ERRY(IJ))**2
      AND=HNS(IJ)/HST
      IND=AND
      IF(AND.LT.0) IND=IND-1
      IND=IND+21
      IF(IND.GT.41.OR.IND.LT.1) IND=41
      NHIST(IND,3)=NHIST(IND,3)+1
  280 CONTINUE
C
C***   CALCULATION OF STD.DEVIATION OF PARAMETERS (IF LNOERR=FALSE)
C
      LRED=.FALSE.
      LBS=.FALSE.
      DO 289 J=1,(NT*MODPAR)
      ERRX(J)=-1.0
      IF(LNOERR) GO TO 289
      DO 288 K=1,(NT*MODPAR)
      E(K)=0.0D0
  288 CONTINUE
      E(J)=1.0D0
      CALL PROCNL(A,INUL,E,(NT*MODPAR),ERRX(J),LRED,LBS,IAT,INULT,IHT)
      ERRX(J)=DSQRT(ERRX(J))
  289 CONTINUE
C
      DO 290 J=1,NT
      JB=(J-1)*MODPAR+1
      AST(J,2)=AST(J,2)/AST(J,1)
      AST(J,3)=DSQRT(AST(J,3)/AST(J,1))
      AND=H(JB)/HST
      IND=AND
      IF(AND.LT.0.) IND=IND-1
      IND=IND+21
      IF(IND.GT.41.OR.IND.LT.1) IND=41
      NHIST(IND,1)=NHIST(IND,1)+1
      AND=AST(J,3)/HST
      IND=AND
      IF(AND.LT.0.) IND=IND-1
      IND=IND+21
      IF(IND.GT.41.OR.IND.LT.1) IND=41
      NHIST(IND,2)=NHIST(IND,2)+1
  290 CONTINUE
      HM=HM/NCRS
      HH=DSQRT(HH/NCRS)
      HHN=DSQRT(HHN/NCRS)
C
      WRITE(6,650)
      TTT=999.999
      ETT=999.999
      BBBMN=0.0
      TTTMN=0.0
      BBBMS=0.0
      TTTMS=0.0
      STRING1='   -N/A-'
      STRING2='   -N/A-'
      STRING3='   -N/A-'
C
      DO 300 I=1,IN
      IB=(I-1)*MODPAR+1
      IT=IB+1
      BBB=H(IB)
      EBB=ERRX(IB)
      IF(.NOT.LNOERR) WRITE(STRING1,669) EBB
      IF(LTILT) THEN
        TTT=H(IT)
        ETT=ERRX(IT)
        WRITE(STRING2,669) TTT
        IF(.NOT.LNOERR) WRITE(STRING3,669) ETT
        ENDIF
      WRITE(6,660)I,INOR(I),BBB,STRING1,STRING2,STRING3,
     .            AST(I,2),AST(I,3)
      WRITE(IUOUT,720)INOR(I),BBB,EBB,AMLON(I),TTT,ETT
      BBBMN=BBBMN+BBB
      TTTMN=TTTMN+TTT
      STRING1='     -  '
      STRING2='     -  '
      STRING3='     -  '
  300 CONTINUE
      BBBMN=BBBMN/IN
      TTTMN=TTTMN/IN
C
      DO 310 J=1,JS
      JB=(IN+J-1)*MODPAR+1
      JT=JB+1
      BBB=H(JB)
      EBB=ERRX(JB)
      IF(.NOT.LNOERR) WRITE(STRING1,669) EBB
      IF(LTILT) THEN
        TTT=H(JT)
        ETT=ERRX(JT)
        WRITE(STRING2,669) TTT
        IF(.NOT.LNOERR) WRITE(STRING3,669) ETT
        ENDIF
      WRITE(6,662)J,ISOU(J),BBB,STRING1,STRING2,STRING3
     .            ,AST(IN+J,2),AST(IN+J,3)
      WRITE(IUOUT,720)ISOU(J),BBB,EBB,AMLON(IN+J),TTT,ETT
      BBBMS=BBBMS+BBB
      TTTMS=TTTMS+TTT
  310 CONTINUE
      BBBMS=BBBMS/JS
      TTTMS=TTTMS/JS
C
  650 FORMAT(//,1X,75(1H*),//,' RESULT OF ESTIMATION '
     +,//,7X,'    REV.NO.       BIAS  STD.DEV.'
     +,      '    TILT  STD.DEV.    MEAN    RMS',/)
  660 FORMAT(1X,I5,'. N',I8,2X,2X,F8.3,A8,2X,2A8,2X,2F8.3)
  662 FORMAT(1X,I5,'. S',I8,2X,2X,F8.3,A8,2X,2A8,2X,2F8.3)
  669 FORMAT(F8.3)
  720 FORMAT(I10,2F8.3,F10.4,2F8.3)
C
      WRITE(6,670)HST
      DO 320 J=1,3
      IF(J.NE.2) WRITE(6,645)(NHIST(I,J),I=1,20)
      IF(J.NE.2) WRITE(6,647)
      WRITE(6,645)(NHIST(I,J),I=21,40)
      WRITE(6,645)(I,I=1,20)
      WRITE(6,649) NHIST(41,J)
  320 CONTINUE
C
      WRITE(6,630) HM,HH,HHN
      STRING2='   -N/A-'
      IF(LTILT) WRITE(STRING2,669) TTTMN
      WRITE(6,680) BBBMN,STRING2
      STRING2='     -  '
      IF(LTILT) WRITE(STRING2,669) TTTMS
      WRITE(6,682) BBBMS,STRING2
      BBBMN=(IN*BBBMN+JS*BBBMS)/NT
      TTTMN=(IN*TTTMN+JS*TTTMS)/NT
      IF(LTILT) WRITE(STRING2,669) TTTMN
      WRITE(6,684) BBBMN,STRING2
C
  670 FORMAT(/,' DISTRIBUTION WITH INTERVALLENGTH ',F5.2,' M'
     +,/,      ' OF          BIAS ESTIMATIONS'
     +,/,      '              RMS OF DIFF.S PER TRACK'
     +,/,      '              ALL RESIDUALS',/)
  680 FORMAT(//,' MEAN OF NORTH-GOING TRACK BIAS/TILTS',F8.3,A8)
  682 FORMAT(   ' MEAN OF SOUTH-GOING TRACK BIAS/TILTS',F8.3,A8)
  684 FORMAT(   ' MEAN OF ALL TRACK BIASES AND TILTS  ',F8.3,A8)
C
      REWIND IUIN
      IF(IUINC.GT.0) REWIND(IUINC)
      REWIND IUOUT
      write(*,*)' END ABIAS '
C
      RETURN
      END
C   *****************************************************************
C
      FUNCTION INDI(IX)
C
C***  FUNCTION INDI FINDS THE INDEX CORRESPONDING TO A REVOLUTION-
C***  NUMBER, IX=INOR(INDI(IX)), OF NORTH-GOING TRACKS.
C
C    PROGRAMMED BY PER KNUDSEN.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MAXPAR=100000)
      COMMON/FUNC/ IN,JS,INOR(MAXPAR),ISOU(MAXPAR)
      INDI=0
      DO 10 I=1,10
      K=IN/2**I+1
      INDI=MIN0(INDI+K,IN)
      IF(INOR(INDI).GT.IX) INDI=INDI-K
      IF(K.LE.1) GO TO 20
   10 CONTINUE
   20 CONTINUE
      IF(INOR(INDI).NE.IX) INDI=0
      RETURN
      END
C
C   *****************************************************************
C
      FUNCTION INDJ(IY)
C
C***  SIMILAR TO FUNCTION INDI FOR SOUTH-GOING TRACKS.
C
C    PROGRAMMED BY PER KNUDSEN.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MAXPAR=100000)
      COMMON/FUNC/ IN,JS,INOR(MAXPAR),ISOU(MAXPAR)
      INDJ=0
      DO 10 J=1,10
      K=JS/2**J+1
      INDJ=MIN0(INDJ+K,JS)
      IF(ISOU(INDJ).GT.IY) INDJ=INDJ-K
      IF(K.LE.1) GO TO 20
   10 CONTINUE
   20 CONTINUE
      IF(ISOU(INDJ).NE.IY) INDJ=0
      RETURN
      END
C
C   *****************************************************************
C
      SUBROUTINE SORTD
C
C***  SUBROUTINE SORTD SORTS THE REVOLUTIONNUMBERS INCRESING.
C
C    PROGRAMMED BY PER KNUDSEN.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MAXPAR=100000)
      COMMON/FUNC/ IN,JS,INOR(MAXPAR),ISOU(MAXPAR)
      IMAX=IN-1
      DO 100 I=1,IMAX
      II=I
      DO 90 J=I,IN
      IF(INOR(J).LT.INOR(II)) II=J
   90 CONTINUE
      IF(II.EQ.I) GO TO 100
      J=INOR(I)
      INOR(I)=INOR(II)
      INOR(II)=J
  100 CONTINUE
      JMAX=JS-1
      DO 200 J=1,JMAX
      JJ=J
      DO 190 I=J,JS
      IF(ISOU(I).LT.ISOU(JJ)) JJ=I
  190 CONTINUE
      IF(JJ.EQ.J) GO TO 200
      I=ISOU(J)
      ISOU(J)=ISOU(JJ)
      ISOU(JJ)=I
  200 CONTINUE
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
C     CONTEND OF ARRAYES
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
C *******************************************************************
C *****                                                         *****
C *****          SUBROUTINE APLANE                              *****
C *****                                                         *****
C *******************************************************************
C
      SUBROUTINE APLANE(NDAT,LSAVE)
C   *****************************************************************
C
C   PROGRAMS ARE COPYRIGHT BY THE AUTHOR AND KORT- OG MATRIKELSTYRELSEN.
C
C   SUBROUTINE APLANE DETERMINES A SURFACE FROM A NDAT OBSERVATIONS.
C   THE COORDINATES ARE TRANSFORMED INTO COORDINATES PSI1 AND PSI2
C   ALONG THE SATELLITES GROUND TRACKS (PSI1 IN NE DIRECTION AND
C   PSI2 IN NW DIRECTION) AND THE SURFACE IS
C
C              H0 = A+B*PSI1+C*PSI2+D*PSI1*PSI2
C
C   THIS SURFACE REPRESENTS THE ZERO-SPACE IN A CROSS-OVER ADJUSTMENT
C   USING BIAS AND TILT PARAMETERS.
C
C   PROGRAMMED BY      PER KNUDSEN
C                      KORT- OG MATRIKELSTYRELSEN
C                      RENTEMESTERVEJ 8
C                      DK-2400 KBH. NV                       18.10.90
C
C   THIS VERSION  19.10.90 PK.
C
C   *****************************************************************
C
C   CONTENT OF ARRAYS    KREV(.)    REVOLUTION NUMBER
C                         LAT(.)     LATITUDE        (*1000000.)
C                         LON(.)     LONGITUDE       (*1000000.)
C                         KALT(.)    GEOID HEIGTH    (* 1000.)
C                         KERR(.)    STD. DEVIATION  (* 1000.)
C
C   MAX NUMBER OF DATA    100000
C
C   *****************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
      parameter(maxdat=100000,irest=63250)
      integer*2 kerr,kcode
      common/field/ krev(maxdat),lat(maxdat),lon(maxdat),
     .              kalt(maxdat),kerr(maxdat),kcode(maxdat),
     .              i4rest(irest),time(maxdat)
      DIMENSION ATA(10),ATY(4),INUL(4)
      LOGICAL LSAVE,LRED,LBS
C
      PI4TH=DATAN(1.0D0)
      DEGRAD=PI4TH/45.0D6
      DEGRAD2=DEGRAD/2.0D0
      RE=6371.0D0
      RELN10=RE/DLOG(10.0D0)
      DEGKM=DEGRAD*RE
C
C FIND MEAN LATITUDE AND MEAN LONGITUDE
C
      S1=0.0D0
      S2=0.0D0
      DO 10 I=1,NDAT
      A1=LAT(I)/1.0D6
      A2=LON(I)/1.0D6
      S1=S1+A1
      S2=S2+A2
   10 CONTINUE
      AMLAT=S1/NDAT
      AMLON=S2/NDAT
      MLAT=AMLAT*1.0D6+0.5
      MLON=AMLON*1.0D6+0.5
C
      WRITE(*,500) NDAT,AMLAT,AMLON
  500 FORMAT('1SUBROUTINE APLANE, VERSION 19.10.90, DETERMINES A',
     ./,     ' SURFACE FROM A NUMBER OF DATA',I8,',',
     ./,     ' HAVING MEAN LATITUDE AND LONGITUDE OF',2F10.4)
C
C FIND A SOUTH-GOING SATELLITE TRACK CROSSING THE MEAN LATITUDE
C
      DO 20 I=2,NDAT
      IF(LAT(I-1).LT.LAT(I)) GO TO 20
      IF(KREV(I-1).NE.KREV(I)) GO TO 20
      IF(LAT(I-1).GT.MLAT.AND.LAT(I).LT.MLAT) GO TO 21
   20 CONTINUE
      WRITE(*,*) 'NO CROSSING FOUND$$??'
      STOP
   21 CONTINUE
C
C  COMPUTE A VECTOR IN THE PSI1 DIRECTION USING A MERCATOR PROJECTION
C  AND THE RATIO R=DY/DX
C
      DX=(LON(I-1)-LON(I))*DEGKM
      Y1=RELN10*DLOG(DTAN(PI4TH+LAT(I-1)*DEGRAD2))
      Y2=RELN10*DLOG(DTAN(PI4TH+LAT(I)*DEGRAD2))
      DY=Y1-Y2
      R=DY/DX
      SLDEG=DSQRT((DCOS(LAT(I)*DEGRAD)*(LON(I)-LON(I-1)))**2+
     .           (1.0D0*(LAT(I)-LAT(I-1)))**2)
      SY=DSQRT(2.0D0)*DY
C SCALE IS IN 100KM PER UNIT
      SCALE=SLDEG*DEGKM/100.0D0/SY
C      WRITE(*,*) LON(I-1),LON(I),DX
C      WRITE(*,*) LAT(I-1),Y1
C      WRITE(*,*) LAT(I),Y2
C      WRITE(*,*) DY,R,SCALE
C
      X0=MLON*DEGKM
      Y0=RELN10*DLOG(DTAN(PI4TH+MLAT*DEGRAD2))
C
      DO 30 I=1,4
      ATY(I)=0.0D0
      INUL(I)=1
      IST=(I*(I-1))/2
      DO 30 J=1,I
      ATA(IST+J)=0.0D0
   30 CONTINUE
C
C SET MATRIX FOR SCALING X-COORDINATES AND FOR ROTATING 45 DEGREES.
C
      A11=R*DSQRT(2.0D0)/2.0D0
      A21=-A11
      A12=DSQRT(2.0D0)/2.0D0
      A22=A12
C
C  FIND MEAN PSI1 AND MEAN PSI2
C
      PSI10=0.0D0
      PSI20=0.0D0
      S1=0.0D0
      S2=0.0D0
      DO 25 I=1,NDAT
      DX=LON(I)*DEGKM-X0
      DY=RELN10*DLOG(DTAN(PI4TH+LAT(I)*DEGRAD2))-Y0
      PSI1=A11*DX+A12*DY
      PSI2=A21*DX+A22*DY
      PSI10=PSI10+PSI1
      PSI20=PSI20+PSI2
      S1=S1+PSI1**2
      S2=S2+PSI2**2
   25 CONTINUE
      PSI10=PSI10/NDAT
      PSI20=PSI20/NDAT
      S1=DSQRT((S1+S2)*0.5/NDAT)
C      WRITE(*,*)' MEAN PSI',PSI10,PSI20
C      WRITE(*,*)' RMS OF PSI',S1
C
C  RESCALE ROT. MATRIX AND SCALE
C
      S1=S1*2.0D0
      A11=A11/S1
      A12=A12/S1
      A21=A21/S1
      A22=A22/S1
      PSI10=PSI10/S1
      PSI20=PSI20/S1
      SCALE=SCALE*S1
C
C SETUP NORMAL EQUATION SYSTEM FOR H0
C
C   (  ATA        ATY:  )
C   (      1 2 4 7     1 )
C   (        3 5 8     2 )
C   (          6 9     3 )
C   (           10     4 )
C
      ATA(1)=NDAT*1.0D0
C
      DO 40 I=1,NDAT
      DX=LON(I)*DEGKM-X0
      DY=RELN10*DLOG(DTAN(PI4TH+LAT(I)*DEGRAD2))-Y0
      PSI1=A11*DX+A12*DY-PSI10
      PSI2=A21*DX+A22*DY-PSI20
      PPSI=PSI1*PSI2
      H=KALT(I)*1.0D-3
C      IF(I.LT.350)
C     . WRITE(*,555) ((LAT(I)-MLAT)/1.D6),((LON(I)-MLON)/1.D6)
C     .              ,(DX/111111.),(DY/111111.),
C     .               (PSI1/111111.),(PSI2/111111.),H
C  555 FORMAT(4(2F15.5,4X))
      ATY(1)=ATY(1)+H
      ATY(2)=ATY(2)+H*PSI1
      ATY(3)=ATY(3)+H*PSI2
      ATY(4)=ATY(4)+H*PPSI
      ATA(2)=ATA(2)+PSI1
      ATA(3)=ATA(3)+PSI1*PSI1
      ATA(4)=ATA(4)+PSI2
      ATA(5)=ATA(5)+PPSI
      ATA(6)=ATA(6)+PSI2*PSI2
      ATA(7)=ATA(7)+PPSI
      ATA(8)=ATA(8)+PSI1*PPSI
      ATA(9)=ATA(9)+PPSI*PSI2
      ATA(10)=ATA(10)+PPSI*PPSI
   40 CONTINUE
C
C      WRITE(*,*) 'MEAN H',(ATY(1)/ATA(1))
C      WRITE(*,*) ATA(1),ATA(2),ATA(4),ATA(7)
C
C COMPUTE SOLUTION
C
      LRED=.TRUE.
      LBS=.TRUE.
      CALL PROCNL(ATA,INUL,ATY,4,DDD,LRED,LBS,10,4,4)
C
      WRITE(*,520) ATY(1),(ATY(2)/SCALE),(ATY(3)/SCALE),
     .             (ATY(4)/SCALE/SCALE)
  520 FORMAT(/,' THE SOLUTION IS  A=',F7.3,' M',
     ./,       '                   B=',F7.3,' M/100KM ALONG TRACK NE',
     ./,       '                   C=',F7.3,' M/100KM ALONG TRACK NW',
     ./,       '                   D=',F10.6,' M/100KM**2 NORTH')
C
      S1=0.0D0
      S2=0.0D0
      S3=0.0D0
      S4=0.0D0
      DO 50 I=1,NDAT
      DX=LON(I)*DEGKM-X0
      DY=RELN10*DLOG(DTAN(PI4TH+LAT(I)*DEGRAD2))-Y0
      PSI1=A11*DX+A12*DY-PSI10
      PSI2=A21*DX+A22*DY-PSI20
      H0=ATY(1)+ATY(2)*PSI1+ATY(3)*PSI2+ATY(4)*PSI1*PSI2
      H=KALT(I)*1.0D-3
      S1=S1+H
      S2=S2+H**2
      H=H-H0
      S3=S3+H
      S4=S4+H**2
      IF(LSAVE) KALT(I)=H*1.0D3+0.5
   50 CONTINUE
      IF(NDAT.GT.1) S2=DSQRT((S2-(S1**2)/NDAT)/(NDAT-1))
      S1=S1/NDAT
      IF(NDAT.GT.1) S4=DSQRT((S4-(S3**2)/NDAT)/(NDAT-1))
      S3=S3/NDAT
      WRITE(*,560) S1,S2,S3,S4
  560 FORMAT(/,' MEAN AND STD.DEV. BEFORE    ',2F8.3,
     ./,       ' AND AFTER REMOVAL OF SURFACE',2F8.3,/)
C
      IF(LSAVE) THEN
        WRITE(*,*) 'NEW VALUES ARE STORED '
       ELSE
        WRITE(*,*) 'NEW VALUES ARE NOT STORED '
       ENDIF
C
      RETURN
      END
C *******************************************************************
C *****                                                         *****
C *****          SUBROUTINE ABCOR                               *****
C *****                                                         *****
C *******************************************************************
C
      SUBROUTINE ABCOR(IUIN,NDAT,LNPASS,NDAT2)
C   *****************************************************************
C
C    SUBROUTINE ABCOR AND THE PROGRAM MODULES INDX, AND SORTD2
C    ARE A MODIFICATION OF PROGRAM ALTCOR2 (GEOSOFT).
C
C    PROGRAMS ARE COPYRIGHT BY THE AUTHOR AND THE GEODETIC
C    INSTITUTE OF DENMARK.
C
C    SUBROUTINE ABCOR CORRECTS ALTIMETER DATA FOR TRACK
C    RELATED BIASES USING THE PARAMETERS ESTIMATED IN
C    SUBROUTINE ABIAS.
C
C    PROGRAMMED BY        PER KNUDSEN
C                         GEODETIC INSTITUTE
C                         GAMLEHAVE ALLE 22
C                         DK-2920 CHARLOTTENLUND           04.02.86.
C
C    THIS VERSION   OCT 11, 1990.   PER KNUDSEN.
C
C   *****************************************************************
C
C   CONNECTED FILES    IUIN  BIAS/TILT PARAMETERS
C
C   *****************************************************************
C
C    INPUT             LNPASS ALTIMETER DATA WITHOUT A BIAS ARE SKIPPED
C                              IF LNPASS=FALSE.
C
C   *****************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(MAXPAR=100000,maxpa1=1500)
      COMMON/FUNC/ NPAR,IREV(MAXPAR),idum,iirest(maxpar)
      parameter(maxdat=100000,irest=63250)
      integer*2 kerr,kcode
      common/field/ krev(maxdat),lat(maxdat),lon(maxdat),
     .              kalt(maxdat),kerr(maxdat),kcode(maxdat),
     .              i4rest(irest),HB(MAXPA1),AMLON(MAXPA1),
     .              HT(MAXPA1),IREVSK(MAXPA1),NSK(MAXPA1),
     .              DUMMY(94000)                          
      LOGICAL LNPASS
C
C***   INPUT OF TRACKNUMBER, BIAS-PARAMETER, MEAN LONGITUDE, AND
C      TILT PARAMETER FROM UNIT IUIN.
C
      DO 10 I=1,MAXPAR
      READ(IUIN,510,END=20)IREV(I),HB(I),AMLON(I),HT(I)
   10 CONTINUE
  510 FORMAT(I10,F8.3,8X,F10.4,F8.3)
      WRITE(6,520) MAXPAR
      STOP
  520 FORMAT(' NUMBER OF PARAMETERS GREATHER THAN',I5)
   20 CONTINUE
      NPAR=I-1
      CALL SORTD2(HB,AMLON,HT)
C
      WRITE(6,530) NPAR,(IREV(I),HB(I),HT(I),I=1,NPAR)
  530 FORMAT(/'SUBROUTINE ABCOR CORRECTS DATA WITH ',I4,' PARAMETERS',
     +/,     ' OF REVOLUTION NUMBER, BIAS, AND TILT',
     +/,200(/,1X,4(I5,2F6.2,',')))
C
C***   CORRECTION OF ALTIMETER DATA.
C
      S1=0.0D0
      S2=0.0D0
      S3=0.0D0
      S4=0.0D0
      NDAT2=0
      NSKIP=0
      DO 100 I=1,NDAT
      H=KALT(I)/1000.0D0
      S1 = S1 + H
      S2 = S2 + H**2
      IDX=INDX(KREV(I))
      IF(IDX.EQ.0) GO TO 110
      H=H-HB(IDX)
      IF(HT(IDX).LT.999.9) H=H-HT(IDX)*(LON(I)/1.0E6-AMLON(IDX))
      GO TO 120
  110 DO 111 KKK=1,NSKIP
      IF(KREV(I).EQ.IREVSK(KKK)) THEN
        NSK(KKK)=NSK(KKK)+1
        GO TO 112
        ENDIF
  111 CONTINUE
      NSKIP=NSKIP+1
      IREVSK(NSKIP)=KREV(I)
      NSK(NSKIP)=1
  112 CONTINUE
      IF(LNPASS) GO TO 100
  120 NDAT2=NDAT2+1
      S3 = S3 + H
      S4 = S4 + H**2
      KREV(NDAT2)=KREV(I)
      LAT(NDAT2)=LAT(I)
      LON(NDAT2)=LON(I)
      KALT(NDAT2)=H*1000.0
      KERR(NDAT2)=KERR(I)
      KCODE(NDAT2)=KCODE(I)
  100 CONTINUE
      IF(NDAT.GT.1) S2=DSQRT((S2-(S1**2)/NDAT)/(NDAT-1))
      S1=S1/NDAT
      IF(NDAT2.GT.1) S4=DSQRT((S4-(S3**2)/NDAT2)/(NDAT2-1))
      S3=S3/NDAT2
C
      REWIND IUIN
C
      WRITE(*,540) NDAT2,S1,S2,S3,S4
      IF(NSKIP.EQ.0) RETURN
      WRITE(*,542) NSKIP
      WRITE(*,550) (IREVSK(KKK),NSK(KKK),KKK=1,NSKIP)
  540 FORMAT(//,' NUMBER OF CORRECTED VALUES',I10,
     ./,        ' MEAN AND STD.DEV BEFORE',2F8.3,
     ./,        ' AND AFTER CORRECTIONS  ',2F8.3,/)
  542 FORMAT(   ' NUMBER OF UNKNOWN TRACKS',I6,',',
     .//,       ' TRACK NUMBER AND NUMBER OF VALUES',/)
  550 FORMAT(1x,6(i6,'',i5,','))
      IF(LNPASS) WRITE(*,560)
  560 FORMAT(/,' THESE VALUES HAVE BEEN ELIMINATED.')
      RETURN
      END
C
C   *****************************************************************
      FUNCTION INDX(IX)
C
C***  FUNCTION INDX FINDS THE INDEX CORRESPONDING TO A REVOLUTION-
C***  NUMBER, IX=IREV(INDX(IX)).
C
C    PROGRAMMED BY PER KNUDSEN.
C
      PARAMETER(MAXPAR=100000)
      COMMON/FUNC/ NPAR,IREV(MAXPAR),idum,iirest(maxpar)
      INDX=0
      DO 10 I=1,10
      K=NPAR/2**I+1
      INDX=MIN0(INDX+K,NPAR)
      IF(IREV(INDX).GT.IX) INDX=INDX-K
      IF(K.LE.1) GO TO 20
   10 CONTINUE
   20 CONTINUE
      IF(IREV(INDX).NE.IX) INDX=0
      RETURN
      END
C
C   *****************************************************************
C
      SUBROUTINE SORTD2(HB,AMLON,HT)
C
C***  SUBROUTINE SORTD2 SORTS THE REVOLUTIONNUMBERS INCRESING.
C
C    PROGRAMMED BY PER KNUDSEN.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MAXPAR=100000)
      COMMON/FUNC/ NPAR,IREV(MAXPAR),idum,iirest(maxpar)
      DIMENSION HB(NPAR),AMLON(NPAR),HT(NPAR)
      IMAX=NPAR-1
      DO 100 I=1,IMAX
      II=I
      DO 90 J=I,NPAR
      IF(IREV(J).LT.IREV(II)) II=J
   90 CONTINUE
      IF(II.EQ.I) GO TO 100
      J=IREV(I)
      IREV(I)=IREV(II)
      IREV(II)=J
      HH=HB(I)
      HB(I)=HB(II)
      HB(II)=HH
      HH=AMLON(I)
      AMLON(I)=AMLON(II)
      AMLON(II)=HH
      HH=HT(I)
      HT(I)=HT(II)
      HT(II)=HH
  100 CONTINUE
      RETURN
      END
