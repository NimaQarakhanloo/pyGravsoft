      PROGRAM PMITES
C $Id: pmites.for 251 2008-12-12 10:26:49Z cct $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     PROGRAM PMITES, VERSION 870918 FORTRAN 77.                       C
c     program for fitting a spherical harmonic expansion to local data.
c     modified version to accept standard format data files and
c     potential coefficient sets in binary or free format without sdev.
c
c     input:
c
c     gravityfile (unit 1)
c     outputfile (unit 6)
c     coefficient file (unit 10)
c     output coefficients (unit 20)
c     lform, iuni, lmin, lmax, itmax,
c     rdlat, rdlon, rplim
c
c     lform = true: formatted potential coefficients
c     lform = false: binary do
c     modifications by rf, dec 88
c
c     NB: this version must be compiled with the /G_FLOATING option on
c     vax to prevent overflow.
c
C     DESCRIPTION OF INPUT PARAMETERS SEE ROUTINE GEOPMI.              C
C                                                                      C
C     TAPE1...     MEAN FREE AIR ANOMALIES (E.G. F3030UCT86).
C     TAPE5...     input PARAMETERS (pmites.inp)
C     TAPE6...     FORMATTED LINE PRINTER UNIT.                        C
C     TAPE7...     FORMATTED CONSOLE UNIT (DATA TERMINAL SCREEN).      C
C     TAPE10...    FORMATTED UNIT, ON WHICH A START GEOPENTIAL MODEL   C
C                  MAY BE STORED BEFORE THE EXECUTION OF PROGRAM       C
C                  PMITES. TAPE10 MAY BE AN EMPTY TAPE.                C
C     TAPE20...    FORMATTED UNIT, ON WHICH THE NEW MODEL WILL BE      C
C                  STORED DURING THE EXECUTION OF ROUTINE GEOPMI.      C
C     TAPE30...    UNFORMATTED SCRATCH UNIT.                           C
C     TAPE31...    UNFORMATTED SCRATCH UNIT.                           C
C     TAPE32...    UNFORMATTED SCRATCH UNIT.                           C
C                                                                      C
C     USED ROUTINES...                                                 C
C     ================                                                 C
C                                                                      C
C     SUBROUTINE GEOBET                          FROM SCU_LIB LIBRRAY. C
C     SUBROUTINE GEOCTR                          FROM SCU_LIB LIBRRAY. C
C     SUBROUTINE GEOEXT                          FROM SCU_LIB LIBRRAY. C
C     SUBROUTINE GEOHIS                          FROM SCU_LIB LIBRRAY. C
C     SUBROUTINE GEOICS                          FROM SCU_LIB LIBRRAY. C
C     SUBROUTINE GEOINL                          FROM SCU_LIB LIBRRAY. C
C     SUBROUTINE GEOIPL                          FROM SCU_LIB LIBRRAY. C
C     SUBROUTINE GEOPMI                          FROM SCU_LIB LIBRRAY. C
C     SUBROUTINE GEOREF                          FROM SCU_LIB LIBRRAY. C
C                                                                      C
C     PROGRAM CREATION...    870901 BY H.-G. WENZEL,                   C
C                            INSTITUT FUER ERDMESSUNG,                 C
C                            UNIVERSITAET HANNOVER,                    C
C                            NIENBURGER STR. 6,                        C
C                            D-3000 HANNOVER 1,                        C
C                            FEDERAL REPUBLIC OF GERMANY.              C
C     PROGRAM MODIFICATION...870918 BY H.-G. WENZEL.                   C
C**********************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER file1*128,file6*128,file10*128,file20*128
      logical lform
      common /form/lform
      DATA IUN1/1/,IUN5/5/,IUN6/6/,IUN7/7/
      DATA IUN10/10/,IUN20/20/,IUN30/30/,IUN31/31/,IUN32/32/,IUN33/33/
c
      open(5,file='pmites.inp',status='old')
      read(5,10) file1
      read(5,10) file6
      read(5,10) file10
      read(5,10) file20
10    format(a128)
c
c      CDATE=DATE( )
c      CTIME=TIME( )
c      WRITE(IUN6,7000) CDATE,CTIME
c      write(*,7000) CDATE,CTIME
c      CALL GEOEXT(IUN6,IUN7,RJTIME)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     FIRST RECORD ON IUN5...                                          C
C     =======================                                          C
C                                                                      C
C     IUNI...      DATA INPUT UNIT =  1                                C
C     LMIN...      MINIMUM DEGREE  =  1                                C
C     LMAX...      MAXIMUM DEGREE  = 30                                C
C                  AND MODEL       =  1 MGAL.                          C
C     ITMAX...     MAXIMUM NUMBER OF ITERATIONS = 5                    C
C                                                                      C
C     SECOND RECORD ON IUN5...                                         C
C     ========================                                         C
C                                                                      C
C     RDLAT...     BLOCKSIZE IN LATITUDE   = 10. DEGREE.               C
C     RDLON...     BLOCKSIZE IN LONGITUDE  = 10. DEGREE.               C
C     RPLIM...     LIST LIMIT FOR DIFFERENCES BETWEEN OBSERVATIONS     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      READ(5,*)lform,IUNI,LMIN,LMAX,ITMAX
      READ(5,*) RDLAT,RDLON,RPLIM
c
      open(6,file=file6,status='unknown')
      WRITE(IUN6,7003) IUNI,LMIN,LMAX,ITMAX
      write(*,7003) IUNI,LMIN,LMAX,ITMAX
      WRITE(IUN6,7004) RDLAT,RDLON,RPLIM
      write(*,7004) RDLAT,RDLON,RPLIM
c
      OPEN(1,FILE=file1,status='old')
      if (lform) OPEN(10,FILE=file10,status='old')
      if (.not.lform) open(10,file=file10,form='unformatted',
     .status='old')
      OPEN(20,FILE=file20,status='unknown')
      OPEN(30,FORM='UNFORMATTED',status='scratch')
      OPEN(31,FORM='UNFORMATTED',status='scratch')
      OPEN(32,FORM='UNFORMATTED',status='scratch')
      OPEN(33,FORM='UNFORMATTED',status='scratch')
c
      CALL GEOPMI(IUN6,IUN7,IUNI,IUN10,IUN20,IUN30,IUN31,IUN32,IUN33,
     1 LMIN,LMAX,RDLAT,RDLON,RPLIM,ITMAX)
c      CALL GEOEXT(IUN6,IUN7,RJTIME)
      STOP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FORMAT STATEMENTS.                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 7000 FORMAT(/' PROGRAM PMITES, VERSION 870918 FORTRAN 77.'//
     1' STARTED AT ',A10,' ',A10//)
 7001 FORMAT(6I10)
 7002 FORMAT(3F10.5)
 7003 FORMAT(/' PILOT PARAMETERS:'//
     1' INPUT UNIT FOR ANOMALIES        =',I10/
     2' MINIMUM DEGREE                  =',I10/
     3' MAXIMUM DEGREE                  =',I10/
     4' MAXIMUM NUMBER OF ITERATIONS    =',I10/)
 7004 FORMAT(
     1' BLOCKSIZE IN LATITUDE           =',F10.5,' DEGREE'/
     2' BLOCKSIZE IN LONGITUDE          =',F10.5,' DEGREE'/
     3' LIST LIMIT                      =',F10.3,' MGAL  '/)
      END
      SUBROUTINE GEOBET(IUN6,IUN7,RDSIG,L,BETAL)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     ROUTINE GEOBET, VERSION 870518 FORTRAN 77.                       C
C                                                                      C
C     SOURCE CODE LIBRARY...SCU_LIB                                    C
C                                                                      C
C     THE ROUTINE GEOBET COMPUTES THE BETA-DAMPING FUNCTION BY         C
C     A RECURSIVE FORMULA.                                             C
C                                                                      C
C     REFERENCE...SJOEBERG, L. 1980: RECURRENCE RELATION FOR THE       C
C                 BETA-N FUNCTION. BULLETIN GEODESIQUE, VOL. 54,       C
C                 69-72, PARIS 1980.                                   C
C                                                                      C
C     USE OF THE ROUTINE GEOBET...                                     C
C     ============================                                     C
C                                                                      C
C     THERE ARE NO RESTRICTIONS FOR L=0,1 OR 2. BUT FOR L GREATER 2,   C
C     ROUTINE GEOBET HAS TO BE CALLED IN A LOOP  BEGINNING WITH L=0,1  C
C     OR 2 UNTIL THE MAXIMUM DEGREE IS REACHED.                        C
C     A CHANGE OF PARAMETERS RDSIG OR L WITHIN THE LOOP DOES NOT       C
C     INFLUENCE THE COMPUTED DAMPING FUNCTION, EXCEPT FOR L=0,1 OR 2.  C
C     EXAMPLE...   DO 10 LP1=1,LMAX1                                   C
C                  L=LP1-1                                             C
C                  CALL GEOBET(IUN6,RDSIG,L,BETAL)                     C
C                  BL(LP1)=BETAL                                       C
C               10 CONTINUE                                            C
C                                                                      C
C                                                                      C
C     INPUT PARAMETER DESCRIPTION...                                   C
C     ==============================                                   C
C                                                                      C
C     IUN6...      FORMATTED LINE PRINTER UNIT.                        C
C     IUN7...      FORMATTED CONSOLE UNIT (DATA TERMINAL SCREEN).      C
C     RDSIG...     AREA OF THE BLOCK IN RADIAN**2.                     C
C     L...         DEGREE, FOR WHICH THE BETA-DAMPING FUNCTION WILL BE C
C                  COMPUTED.                                           C
C                                                                      C
C     OUTPUT PARAMETER DESCRIPTION...                                  C
C     ===============================                                  C
C                                                                      C
C     BETAL...     BETA DAMPING FUNCTION.                              C
C                                                                      C
C     USED ROUTINES...NONE.                                            C
C     ================                                                 C
C                                                                      C
C                                                                      C
C     NUMERICAL ACCURACY... UNCHECKED                                  C
C     =====================                                            C
C                                                                      C
C     EXECUTION TIME...                                                C
C     =================                                                C
C                                                                      C
C     8.4*10**-6 SEC CPU TIME PER CALL OF ROUTINE GEOBET ON CDC CYBER  C
C     76 OF RRZN HANNOVER, WHICH MEANS 0.084 SEC CPU TIME FOR A COM-   C
C     PLETE LOOP WITH L=0...10 000.                                    C
C     EXECUTION TIME NOT DETERMINED FOR CDC CYBER 990 OF RRZN HANNOVER.C
C                                                                      C
C     PROGRAM TESTS...                                                 C
C     ================                                                 C
C                                                                      C
C     SUCCESSFULLY COMPLETED 830614 BY H.-G.WENZEL.                    C
C                                                                      C
C     ROUTINE CREATION...    830408 BY H.-G.WENZEL,                    C
C                            INSTITUT FUER ERDMESSUNG,                 C
C                            UNIVERSITAET HANNOVER,                    C
C                            NIENBURGER STR.6,                         C
C                            D-3000 HANNOVER 1,                        C
C                            FEDERAL REPUBLIC OF GERMANY.              C
C                            PHONE... 0511-7622796.                    C
C     ROUTINE MODIFICATION...870128 08:00 BY H.-G-WENZEL.              C
C**********************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      save
      BETAL=1.
      IF(RDSIG.EQ.0.) RETURN
      IF(L.GT.2) GOTO 30
      IF(L.EQ.2) GOTO 20
      IF(L.EQ.1) GOTO 10
      IF(L.NE.0) GOTO 5000
    1 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     L EQUAL 0.                                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      BETAL=1.
      RETURN
   10 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     L EQUAL 1.                                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PI=4.*datan(1.d0)
C     PI=4.*ATAN(1.)
      CPSI=1.-RDSIG/(2.*PI)
      BETAL=(1.+CPSI)/2.
      RETURN
   20 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     L EQUAL 2.                                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PI=4.*datan(1.d0)
C     PI=4.*ATAN(1.)
      CPSI=1.-RDSIG/(2.*PI)
      BETAL2=1.
      BETAL1=(1.+CPSI)/2.
      FL=2.
      LN=2
   30 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     RECURSION ALGORITHM FOR L GREATER OR EQUAL 2.                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(L.NE.LN) GOTO 5010
      LN=LN+1
      BETAL=((2.*FL-1.)*CPSI*BETAL1-(FL-2.)*BETAL2)/(FL+1.)
      BETAL2=BETAL1
      BETAL1=BETAL
      FL=FL+1.
      RETURN
 5000 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DEGREE L LESS 0.                                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(IUN6,7001) L
      write(*,7001) L
      L=0
      GOTO 1
 5010 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     WRONG SEQUENCE OF DEGREE L.                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(IUN6,7002) L,LN
      write(*,7002) L,LN
      L=LN
      GOTO 30
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FORMAT STATEMENTS.                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 7001 FORMAT(//' *****ERROR IN ROUTINE GEOBET, VERSION 870128 FTN77.'/
     1' *****THE DEGREE L USED IN THE CALL OF GEOBET HAS TO BE GREATER',
     2' OR QUAL 0 BUT IS',I5/
     3' *****THE DAMPING FUNCTION HAS BEEN PUT TO 1.'/
     4' *****THE DEGREE L HAS BEEN PUT TO 0.'/
     5' *****THE EXECUTION WILL BE CONTINUED.'/)
 7002 FORMAT(//' *****ERROR IN ROUTINE GEOBET, VERSION 870128 FTN77.'/
     1' *****WRONG SEQUENCE OF DEGREE L =',I5/
     2' *****DEGREE L IS MODIFIED TO',I5,' AND THE EXECUTION WILL BE',
     3' CONTINUED.'/)
      END
      SUBROUTINE GEOCTR(IUN6,IUN7,KCI,CLAT,SLAT,CLON,SLON,H,CT,ST,RR,
     1 X,Y,Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     ROUTINE GEOCTR, VERSION 870518 FORTRAN 77.                       C
C                                                                      C
C     SOURCE CODE LIBRARY...SCU_LIB                                    C
C                                                                      C
C     THE ROUTINE GEOCTR COMPUTES TRANSFORMATIONS BETWEEN ELLIPSOIDAL, C
C     SPHERICAL AND CARTESIAN COORDINATES.                             C
C                                                                      C
C     INPUT PARAMETER DESCRIPTION...                                   C
C     ==============================                                   C
C                                                                      C
C     IUN6...      UNIT NUMBER OF LINE PRINTER OUTPUT.                 C
C     IUN7...      FORMATTED CONSOLE UNIT (DATA TERMINAL SCREEN).      C
C     KCI...       DEFINES THE INPUT COORDINATE SYSTEM.                C
C                  FOR KCI=1, ELLIPSOIDAL COORDINATES CLAT, SLAT,      C
C                  CLON, SLON, H ARE ASSUMED TO BE GIVEN. SPHERICAL    C
C                  COORDINATES CT,ST,RR AND CARTESIAN COORDINATES      C
C                  X,Y,Z WILL BE COMPUTED.                             C
C                  FOR KCI=2, SPHERICAL COORDINATES CT,ST,CLON,SLON,RR C
C                  ARE ASSUMED TO BE GIVEN. ELLIPSOIDAL COORDINATES    C
C                  CLAT,SLAT,H AND CARTESIAN COORDINATES X,Y,Z WILL    C
C                  BE COMPUTED.                                        C
C                  FOR KCI=3, CARTESIAN COORDINATES X,Y,Z ARE          C
C                  ASSUMED TO BE GIVEN.ELLIPSOIDAL COORDINATES CLAT,   C
C                  SLAT,CLON,SLON,H AND SPHERICAL COORDINATES CT,ST,   C
C                  CLON,SLON,RR WILL BE COMPUTED.                      C
C                                                                      C
C     INPUT/OUTPUT PARAMETER DESCRIPTION...                            C
C     =====================================                            C
C                                                                      C
C     DEPENDING ON PARAMETER KCI, THE FOLLOWING PARAMETERS ARE         C
C     EITHER INPUT OR OUTPUT PARAMETERS.                               C
C                                                                      C
C     CLAT...      COS OF ELLPSOIDAL LATITUDE.                         C
C     SLAT...      SIN OF ELLPSOIDAL LATITUDE.                         C
C     CLON...      COS OF ELLIPSOIDAL OR SPHERICAL LONGITUDE.          C
C     SLON...      SIN OF ELLIPSOIDAL OR SPHERICAL LONGITUDE.          C
C     H...         ELLPSOIDAL HEIGHT IN METER.                         C
C     CT...        COS OF SPHERICAL POLAR DISTANCE (CO-LATITUDE).      C
C     ST...        SIN OF SPHERICAL POLAR DISTANCE (CO-LATITUDE).      C
C     RR...        RADIUS IN METER.                                    C
C     X,Y,Z...     CARTESIAN COORDINATES IN METER.                     C
C                                                                      C
C     COMMON BLOCK DESCRIPTION...                                      C
C     ===========================                                      C
C                                                                      C
C     COMMON/REFELL/... THE PARAMETERS OF COMMON/REFELL/ HAVE TO BE    C
C     DEFINED BY CALLING ROUTINE GEOREF BEFORE THE EXECUTION OF        C
C     ROUTINE GEOCTR. FOR THE EXECUTION OF ROUTINE GEOCTR, THERE HAVE  C
C     TO BE DEFINED ONLY THE PARAMETERS A AND E2.                      C
C                                                                      C
C     GM...        GEOCENTRIC GRAVITATIONAL CONSTANT IN                C
C                  METER**3/SEC**2.                                    C
C     A...         MAJOR SEMI AXIS IN METER.                           C
C     F...         FLATTENING.                                         C
C     OM...        ROTATION SPEED IN RADIANS/SEC.                      C
C     RMEAN...     MEAN EARTH'S RADIUS IN METER.                       C
C     GAMEAN...    MEAN EARTH'S GRAVITY IN METER/SEC**2.               C
C     GE...        EQUATORIAL NORMAL GRAVITY IN METER/SEC**2.          C
C     RK...        CONSTANT FOR SOMIGLIANA NORMAL GRAVITY FORMULA.     C
C     E2...        SQUARE OF FIRST ECCENTRICITY.                       C
C     ES2...       SQUARE OF SECOND ECCENTRICITY.                      C
C     U0...        NORMAL POTENTIAL OF THE LEVEL ELLIPSOID IN          C
C                  METER**2/SEC**2.                                    C
C     CN...        FULLY NORMALIZED ZONAL HARMONIC COEFFICIENTS OF     C
C                  THE ELLIPSOIDAL NORMAL GRAVITY POTENTIAL UP TO      C
C                  DEGREE 10 . THE COEFFICIENT    C(L,0) IS STORED     C
C                  IN CN(L+1), THE COEFFICIENT C(0,0) IS SET TO 1 AND  C
C                  STORED IN CN(1).                                    C
C                                                                      C
C     COMMON/REFSYS/...NAME OF THE USED REFERENCE ELLIPSOID.           C
C     CREFSY...    NAME OF THE USED REFERENCE ELLIPSOID (CHARACTER*10).C
C                                                                      C
C     USED ROUTINES...  NONE.                                          C
C     ================                                                 C
C                                                                      C
C     NUMERICAL ACCURACY...                                            C
C     =====================                                            C
C                                                                      C
C     FOR KCI= 2 AND 3, ITERATION OF CLAT AND SLAT TO 10**-12 RADIAN   C
C     RESP. 2*10**-7 ARC SEC.                                          C
C                                                                      C
C     EXECUTION TIME...                                                C
C     =================                                                C
C                                                                      C
C     FOR KCI=1, 2.0*10**-5 SEC CPU TIME ON CDC CYBER 76 OF RRZN       C
C                HANNOVER.                                             C
C     FOR KCI=2, 3.2*10**-5 SEC CPU TIME ON CDC CYBER 76 OF RRZN       C
C                HANNOVER FOR HEIGHTS LESS 10 KM.                      C
C     FOR KCI=3, 6.0*10**-5 SEC CPU TIME ON CDC CYBER 76 OF RRZN       C
C                HANNOVER FOR HEIGHTS LESS 10 KM.                      C
C                                                                      C
C     ROUTINE TESTS...                                                 C
C     ================                                                 C
C                                                                      C
C                                                                      C
C     ROUTINE CREATION...    820722 BY H.-G.WENZEL,                    C
C                            INSTITUT FUER ERDMESSUNG,                 C
C                            UNIVERSITAET HANNOVER,                    C
C                            NIENBURGER STR.6,                         C
C                            D-3000 HANNOVER 1,                        C
C                            FEDERAL REPUBLIC OF GERMANY.              C
C     ROUTINE MODIFICATION...870518 BY H.-G.WENZEL.                    C
C**********************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER CREFSY*10
      COMMON/REFELL/GM,A,F,OM,RMEAN,GAMEAN,GE,RK,E2,ES2,U0,CN(11)
      COMMON/REFSYS/ CREFSY
      IF(KCI.EQ.1) GOTO 100
      IF(KCI.EQ.2) GOTO 200
      IF(KCI.EQ.3) GOTO 300
      GOTO 5000
  100 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ELLIPSOIDAL COORDINATES GIVEN, COMPUTE SPHERICAL AND CARTESIAN   C
C     COORDINATES.                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RN=A/DSQRT(1.-E2*SLAT*SLAT)
C     RN=A/SQRT(1.-E2*SLAT*SLAT)
      P=(RN+H)*CLAT
      X=P*CLON
      Y=P*SLON
      Z=(RN*(1.-E2)+H)*SLAT
      RR=DSQRT(X*X+Y*Y+Z*Z)
C     RR=SQRT(X*X+Y*Y+Z*Z)
      CT=Z/RR
      ST=P/RR
      RETURN
  200 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SPHERICAL COORDINATES GIVEN, COMPUTE CARTESIAN AND ELLIPSOIDAL   C
C     COORDINATES.                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      P=RR*ST
      X=P*CLON
      Y=P*SLON
      Z=RR*CT
      GOTO 301
  300 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CARTESIAN COORDINATES GIVEN, COMPUTE ELLIPSOIDAL AND SPHERICAL   C
C     COORDINATES.                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      P=DSQRT(X*X+Y*Y)
      RR=DSQRT(X*X+Y*Y+Z*Z)
C     P=SQRT(X*X+Y*Y)
C     RR=SQRT(X*X+Y*Y+Z*Z)
      CT=Z/RR
      ST=P/RR
      CLON=1.
      SLON=0.
      IF(P.LT.1.E-10) GOTO 301
      CLON=X/P
      SLON=Y/P
  301 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ITERATIVE COMPUTATION OF CLAT,SLAT AND H.                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RX=1./DSQRT(1.-E2*(2.-E2)*ST*ST)
C     RX=1./SQRT(1.-E2*(2.-E2)*ST*ST)
      CLAT0=ST*(1.-E2)*RX
      SLAT0=CT*RX
  302 CONTINUE
      RN=A/DSQRT(1.-E2*SLAT0*SLAT0)
C     RN=A/SQRT(1.-E2*SLAT0*SLAT0)
      H=P*CLAT0+Z*SLAT0-RN*(1.-E2*SLAT0*SLAT0)
      RJ=E2*RN/(RN+H)
      RX=1./DSQRT(1.-RJ*(2.-RJ)*ST*ST)
C     RX=1./SQRT(1.-RJ*(2.-RJ)*ST*ST)
      CLAT=(1.-RJ)*ST*RX
      SLAT=CT*RX
      SDLAT=SLAT*CLAT0-CLAT*SLAT0
      IF(SDLAT*SDLAT.LT.1.E-24) RETURN
      CLAT0=CLAT
      SLAT0=SLAT
      GOTO 302
 5000 WRITE(IUN6,7001) KCI
      write(*,7001) KCI
      STOP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FORMAT STATEMENTS.                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 7001 FORMAT(' ***ERROR IN ROUTINE GEOCTR, VERSION 870518 FTN77.'/
     2' ***WRONG PARAMETER KCI=',I5/
     3' ***ROUTINE GEOCTR STOPS THE EXECUTION.'/)
      END
c
      SUBROUTINE GEOHIS(IUN6,IUN7,IENT,CNAM,XMIN,XMAX,RDX,J,X,NEWP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     ROUTINE GEOHIS, VERSION 870602 FORTRAN 77.                       C
C                                                                      C
C     SOURCE CODE LIBRARY...SCU_LIB                                    C
C     OBJECT CODE LIBRARY...WZLIB 870602                               C
C                                                                      C
C     THE ROUTINE GEOHIS COMPUTES HISTOGRAMS FOR UP TO 20 DIFFERENT    C
C     PARAMETERS IN PARALLEL EXECUTION AND PRINTS THE HISTOGRAMS       C
C     ON LINE PRINTER.                                                 C
C                                                                      C
C     INPUT PARAMETER DESCRIPTION...                                   C
C     ==============================                                   C
C                                                                      C
C     IUN6...      FORMATTED LINE PRINTER UNIT.                        C
C     IUN7...      FORMATTED CONSOLE UNIT (DATA TERMINAL SCREEN).      C
C     IENT...      ENTRY PARAMETER FOR ROUTINE GEOHIS.                 C
C                  FOR IENT=1, ROUTINE GEOHIS WILL BE INITIALIZED      C
C                  FOR PARAMETER NO. J.                                C
C                  FOR IENT=2, THE HISTOGRAM WILL BE UPDATED FOR       C
C                  SAMPLE VALUE X OF PARAMETER NO. J.                  C
C                  FOR IENT=3, THE HISTOGRAM WILL BE PRINTED FOR       C
C                  PARAMETER NO. J AND AGAIN INITIALIZED FOR PARAMETER C
C                  NO. J.                                              C
C     CNAM...      ARRAY OF PARAMETER NAMES (CHARACTER*10).            C
C     XMIN...      ARRAY OF LOWER BORDERS FOR HISTOGRAMS.              C
C     XMAX...      ARRAY OF UPPER BORDERS FOR HISTOGRAMS.              C
C     RDX...       ARRAY OF INCREMENTS FOR HISTOGRAMS.                 C
C     J...         PARAMETER NO.                                       C
C     X...         SAMPLE VALUE OF PARAMETER NO. J.                    C
C     NEWP...      PARAMETER FOR BEGINNING A NEW PAGE.                 C
C                  NEWP=1...HISTOGRAM WILL BE PRINTED ON A NEW PAGE    C
C                  FOR ENTRY PARAMETER IENT=3.                         C
C                  X HAS TO BE DEFINED ONLY FOR ENTRY PARAMETER IENT=2.C
C                                                                      C
C     PROGRAM CREATION...    820420 BY G.WEBER,                        C
C                            INSTITUT FUER ERDMESSUNG,                 C
C                            UNIVERSITAET HANNOVER,                    C
C                            NIENBURGER STR.6,                         C
C                            D-3000 HANNOVER 1,                        C
C                            FEDERAL REPUBLIC OF GERMANY.              C
C     PROGRAM MODIFICATION...830725 BY H.-G.WENZEL                     C
C     PROGRAM MODIFICATION...870521 BY H.-G.WENZEL                     C
C**********************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER CNAM(20)*10,CZEI(40)*1
      DIMENSION XMIN(20),XMAX(20),RDX(20),KMAX(20)
      DIMENSION RMX(20),RMS(20),RXMIN(20),RXMAX(20)
      DIMENSION N(52,20)
c      save
      COMMON/HIS/ N
      IF(IENT.NE.2) GOTO 100
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IENT=2. UPDATING OF ARRAYS N,RMX,RMS FOR PARAMETER NO. J         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      K=(X-XMIN(J))/RDX(J)+2.
      IF(K.LT.1) K=1
      IF(K.GT.KMAX(J)) K=KMAX(J)
      N(K,J)=N(K,J)+1
      RMX(J)=RMX(J)+X
      RMS(J)=RMS(J)+X*X
      IF(X.LT.RXMIN(J)) RXMIN(J)=X
      IF(X.GT.RXMAX(J)) RXMAX(J)=X
      RETURN
  100 IF(IENT.NE.1) GOTO 300
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IENT=1. INITIALIZATION OF ROUTINE GEOHIS FOR PARAMETER NO. J.    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 110 K=1,52
  110 N(K,J)=0
      RMX(J)=0.
      RMS(J)=0.
      RXMIN(J)= 999999.9
      RXMAX(J)=-999999.9
      KMAX(J)=(XMAX(J)-XMIN(J))/RDX(J)+2.
      IF(KMAX(J).GT.52) KMAX(J)=52
      RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IENT=3. PRINT HISTOGRAM FOR PARAMETER NO. J.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  300 NSU=0
      NMA=0
      DO 310 K=1,KMAX(J)
      IF(N(K,J).GT.NMA) NMA=N(K,J)
  310 NSU=NSU+N(K,J)
      IF(NSU.EQ.0) RETURN
      SCALE=40./FLOAT(NMA)
      IF(NEWP.EQ.1) WRITE(IUN6,7006)
      WRITE(IUN6,7000)
      RMX(J)=RMX(J)/FLOAT(NSU)
      RMS(J)=RMS(J)-RMX(J)*RMX(J)*FLOAT(NSU)
      RMS(J)=DSQRT(RMS(J)/FLOAT(NSU))
C     RMS(J)=SQRT(RMS(J)/FLOAT(NSU))
      WRITE(IUN6,7003) CNAM(J),NSU,RMX(J),RMS(J),RXMIN(J),RXMAX(J)
      write(*,311) CNAM(J),NSU,RMX(J),RMS(J),RXMIN(J),RXMAX(J)
  311 format(/' parameter: ',a10,/,
     .' n:',i6,', mean:',f9.2,', rms':,f9.2,', min:',f9.2,', max',f9.2)
      DO 320 K=1,KMAX(J)
      DO 330 L=1,40
  330 CZEI(L)=' '
      CZEI(1)='.'
      IZEIE=FLOAT(N(K,J))*SCALE
      IF(IZEIE.LE.40) GOTO 340
      CZEI(40)='+'
      IZEIE=40
  340 CONTINUE
      IF(IZEIE.LT.1) GOTO 360
      DO 350 L=1,IZEIE
  350 CZEI(L)='*'
  360 CONTINUE
      XMI=FLOAT(K-2)*RDX(J)+XMIN(J)
      XMA=FLOAT(K-1)*RDX(J)+XMIN(J)
      PROZ=FLOAT(N(K,J))/FLOAT(NSU)
      IF(K.NE.1.AND.K.NE.KMAX(J)) WRITE(IUN6,7001) XMI,XMA,PROZ,N(K,J),
     .(CZEI(L),L=1,40)
      IF(K.EQ.1) WRITE(IUN6,7002) XMA,PROZ,N(K,J),(CZEI(L),L=1,40)
      IF(K.EQ.KMAX(J)) WRITE(IUN6,7004) XMI,PROZ,N(K,J),(CZEI(L),L=1,40)
  320 CONTINUE
  370 WRITE(IUN6,7005)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     RE-INITIALIZATION FOR PARAMETER NO.J.                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 380 K=1,52
  380 N(K,J)=0
      RMX(J)=0.
      RMS(J)=0.
      RXMIN(J)=999999.9
      RXMAX(J)=-999999.9
      RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FORMAT STATEMENTS.                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 7000 FORMAT(/' ROUTINE GEOHIS, VERSION 870602 FTN77.'/)
 7001 FORMAT(1X,F10.3,'..',F10.3,F8.3,I8,1X,40A1)
 7002 FORMAT(1X,'-INFINITY ..',F10.3,F8.3,I8,1X,40A1)
 7004 FORMAT(1X,F10.3,'..  INFINITY',F8.3,I8,1X,40A1)
 7003 FORMAT(' HISTOGRAM FOR PARAMETER',2X,A10//
     1' NUMBER OF SAMPLE VALUES',10X,I10/
     2' MEAN OF SAMPLE VALUES  ',10X,F10.3/
     3' RMS  OF SAMPLE VALUES  ',10X,F10.3/
     4' MINIMAL SAMPLE VALUE   ',10X,F10.3/
     5' MAXIMAL SAMPLE VALUE   ',10X,F10.3/)
 7005 FORMAT(' ')
 7006 FORMAT(1H1)
      END
      SUBROUTINE GEOICS(IUN6,IUN7,IENT,RDLONR,CLON,SLON,MAX,CMI,SMI)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     ROUTINE GEOICS, VERSION 870529 FORTRAN 77.                       C
C                                                                      C
C     SOURCE CODE LIBRARY...SCU_LIB                                    C
C                                                                      C
C     THE ROUTINE GEOICS COMPUTES THE INTEGRALS OF COS(M*RLON) AND     C
C     SIN(M*RLON) BY RECURSION FORMULAS.                               C
C                                                                      C
C     FOR THE COMPUTATION OF THE INTEGRALS, THE FOLLOWING RELATIONS    C
C     ARE USED BECAUSE OF THEIR HIGH NUMERICAL STABILITY EVEN FOR      C
C     SMALL LONGITUDE DIFFERENCES...                                   C
C                                                                      C
C     RLON=0.5*(RLONW+RLONE)                                           C
C     RDLONR=(RLONE-RLONW)                                             C
C     CMI=2.*COS(M*RLON)*SIN(0.5*M*RDLONR)/M                           C
C     SMI=2.*SIN(M*RLON)*SIN(0.5*M*RDLONR)/M                           C
C                                                                      C
C     INPUT PARAMETER DESCRIPTION...                                   C
C     ==============================                                   C
C                                                                      C
C     IUN6...      FORMATTED LINE PRINTER UNIT.                        C
C     IUN7...      FORMATTED CONSOLE UNIT (DATA TERMINAL SCREEN).      C
C     IENT...      ENTRY PARAMETER FOR ROUTINE GEOICS.                 C
C                  FOR IENT=1, THE STORAGE ARRAY SMDL, WHICH DEPENDS   C
C                  FROM THE LONGITUDE DIFFERENCE RDLONR ONLY, WILL BE  C
C                  ESTABLISHED.                                        C
C                  FOR IENT=2, THE INTEGRALS FOR A SPECIFIC LONGITUDE  C
C                  WILL BE COMPUTED AND STORED IN ARRAYS CMI ANS SMI.  C
C     RDLONR...    LONGITUDE DIFFERENCE (EAST MINUS WEST BORDER)       C
C                  IN RADIAN.                                          C
C     CLON...      COS OF LONGITUDE OF BLOCK MEAN POINT.               C
C     SLON...      SIN OF LONGITUDE OF BLOCK MEAN POINT.               C
C     MAX...       MAXIMUM ORDER M, FOR WHICH THE INTEGRALS WILL BE    C
C                  COMPUTED.                                           C
C                                                                      C
C     OUTPUT PARAMETER DESCRIPTION...                                  C
C     ===============================                                  C
C                                                                      C
C     CMI...       ARRAY OF INTEGRALS OF COS(M*RLON),                  C
C                  M=0...MAX. THE INTEGRAL OF COS(M*RLON) IS STORED    C
C                  IN CMI(M+1).                                        C
C     SMI...       ARRAY OF INTEGRALS OF SIN(M*RLON),                  C
C                  M=0...MAX. THE INTEGRAL OF SIN(M*RLON) IS STORED    C
C                  IN SMI(M+1).                                        C
C                                                                      C
C     USED ROUTINES...NONE                                             C
C     ================                                                 C
C                                                                      C
C     NUMERICAL ACCURACY...                                            C
C     =====================                                            C
C                                                                      C
C     APPROXIMATELY 5*10**-13 FOR THE COMPUTED INTEGRALS.              C
C                                                                      C
C     EXECUTION TIME...                                                C
C     =================                                                C
C                                                                      C
C     ABOUT 0.001 SEC CPU TIME FOR A CALL WITH IENT=1 AND ABOUT        C
C     0.00041 SEC CPU TIME FOR A CALL WITH IENT=2 AND MAXIMUM ORDER    C
C     500 ON CONTROL DATA CYBER 990 OF RRZN HANNOVER WITH OPTIMIZATION C
C     LEVEL OL=HIGH.                                                   C
C                                                                      C
C     ROUTINE CREATION...    830205 BY H.-G.WENZEL,                    C
C                            INSTITUT FUER ERDMESSUNG                  C
C                            UNIVERSITAET HANNOVER                     C
C                            NIENBURGER STR.6                          C
C                            D-3000 HANNOVER 1                         C
C                            FEDERAL REPUBLIC OF GERMANY               C
C                            PHONE: 0511-7622796.                      C
C     ROUTINE MODIFICATION...870529 BY H.-G.WENZEL.                    C
C**********************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CMI(501),SMI(501),SMDL(501)
c      save
      IF(IENT.EQ.2) GOTO 200
      IF(MAX.GT.500) GOTO 5000
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IENT=1, INITIALIZATION OF ARRAY SMDL.                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CDLONM=1.
      SDLONM=0.
      SMDL(1)=0.
      MAXP1=MAX+1
      CDLON=DCOS(0.5*RDLONR)
      SDLON=DSIN(0.5*RDLONR)
C     CDLON=COS(0.5*RDLONR)
C     SDLON=SIN(0.5*RDLONR)
      RM=1.
      DO 100 MP1=2,MAXP1
      C=CDLON*CDLONM-SDLON*SDLONM
      SDLONM=SDLON*CDLONM+CDLON*SDLONM
      CDLONM=C
      SMDL(MP1)=2./RM*SDLONM
      RM=RM+1.
  100 CONTINUE
      RETURN
  200 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IENT=2.                                                          C
C     M EQUAL ZERO. INITIALIZATION.                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CMI(1)=RDLONR
      SMI(1)=0.
      CLONM=CLON
      SLONM=SLON
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     M GREATER ZERO.                                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 210 MP1=2,MAXP1
      CMI(MP1)=CLONM*SMDL(MP1)
      SMI(MP1)=SLONM*SMDL(MP1)
      C=CLON*CLONM-SLON*SLONM
      SLONM=SLON*CLONM+CLON*SLONM
      CLONM=C
  210 CONTINUE
      RETURN
 5000 WRITE(IUN6,7001) MAX
      write(*,7001) MAX
      STOP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FORMAT STATEMENTS.                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 7001 FORMAT(/,' *** ERROR IN ROUTINE GEOICS, VERSION 870521 FTN77.'/
     1 ' *** MAXIMUM ORDER MAX=',I5,' EXCCEDS THE DIMENSIONS',
     2 ' OF ARRAYS CMI AND SMI.',/,
     3 ' *** ROUTINE GEOICS STOPS THE EXECUTION.')
      END
      SUBROUTINE GEOINL(IUN6,IUN7,IENT,M,CTN,STN,CTS,STS,NDEG,PLM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     ROUTINE GEOINL, VERSION 871026 FORTRAN 77.                       C
C                                                                      C
C     SOURCE CODE LIBRARY...SCU_LIB                                    C
C                                                                      C
C     THE ROUTINE GEOINL COMPUTES THE INTEGRALS OF FULLY NORMALIZED    C
C     LEGENDRE FUNCTIONS FROM RECURSION FORMULAS.                      C
C                                                                      C
C     REFERENCE... PAUL,M. 1978... RECURRENCE RELATIONS FOR THE        C
C                  INTEGRALS OF ASSOCIATED LEGENDRE FUNCTIONS.         C
C                  BULLETIN GEODESIQUE, VOL.52 NO.3, P. 177-190,       C
C                  PARIS 1978.                                         C
C     REFERENCE... WENZEL,H.-G. 1985.. ZUR ANWENDUNG UND BERECHNUNG    C
C                  VON HOCHAUFLOESENDEN KUGELFUNKTIONSMODELLEN FUER    C
C                  DAS GRAVITATIONSPOTENTIAL DER ERDE.                 C
C                  WISSENSCHAFTLICHE ARBEITEN DER FACHRICHTUNG         C
C                  VERMESSUNGSWESEN DER UNIVERSITAET HANNOVER,         C
C                  HANNOVER 1985.                                      C
C                                                                      C
C     USE OF ROUTINE GEOINL...                                         C
C     ========================                                         C
C                                                                      C
C     BEFORE COMPUTING THE INTEGRALS OF FULLY NORMAILZED LEGENDRE      C
C     FUNCTIONS, THE ROUTINE GEOINL  HAS TO BE CALLED WITH IENT=1 IN   C
C     ORDER TO INITIALIZE ARRAYS WITH SQRT(K) AND 1./SQRT(K).          C
C     FOR COMPUTING THE INTEGRALS OF FULLY NORMALIZED LEGENDRE         C
C     FUNCTIONS OF A CERTAIN ORDER, ROUTINE GEOINL  HAS TO BE CALLED   C
C     WITH IENT=2 IN A LOOP BEGINNING WITH M=0 UNTIL THE DESIRED ORDER C
C     M IS REACHED.                                                    C
C                                                                      C
C     EXAMPLE...                                                       C
C                  IENT=1                                              C
C                  CALL GEOINL(IUN6,IUN7,IENT,...,NDEG,....)           C
C                  IENT=2                                              C
C                  NDEG1=NDEG+1                                        C
C                  DO 10 MP1=1,NDEG1                                   C
C                  M=MP1-1                                             C
C                  CALL GEOINL(IUN6,IUN7,IENT,M,............)          C
C                C USE OF COMPUTED INTEGRALS OF ORDER M                C
C               10 CONTINUE                                            C
C                                                                      C
C     INPUT PARAMETER DESCRIPTION...                                   C
C     ==============================                                   C
C                                                                      C
C     IUN6...      FORMATTED LINE PRINTER UNIT.                        C
C     IUN7...      FORMATTED CONSOLE UNIT (DATA TERMINAL SCREEN).      C
C     IENT...      ENTRY PARAMETER FOR ROUTINE GEONLF.                 C
C                  FOR IENT=1, THE ARRAYS RSR AND RISR WILL BE         C
C                  INITIALIZED. IN THAT CASE, ONLY PARAMETER NDEG HAS  C
C                  TO BE DEFINED.                                      C
C                  FOR IENT=2, THE INTEGRALS OF FULLY NORMALIZED       C
C                  LEGENDRE FUNCTIONS OF ORDER M WILL BE COMPUTED FOR  C
C                  DEGREE L=M...NDEG.                                  C
C     M...         ORDER, FOR WHICH THE INTEGRALS OF FULLY NORMALIZED  C
C                  LEGENDRE FUNCTIONS WILL BE COMPUTED FOR DEGREE      C
C                  L=M...NDEG.                                         C
C     CTN...       COS OF POLAR DISTANCE OF NORTHERN BORDER.           C
C     STN...       SIN OF POLAR DISTANCE OF NORTHERN BORDER.           C
C     CTS...       COS OF POLAR DISTANCE OF SOUTHERN BORDER.           C
C     STS...       SIN OF POLAR DISTANCE OF SOUTHERN BORDER.           C
C     NDEG...      MAXIMUM DEGREE AND ORDER, FOR WHICH THE INTEGRALS   C
C                  OF FULLY NORMALIZED LEGENDRE FUNCTIONS WILL BE      C
C                  COMPUTED. NDEG IS  RESTRICTED TO BE LESS OR EQUAL   C
C                  500 DUE TO THE DIMENSION STATEMENTS IN ROUTINE      C
C                  GEOINL.                                             C
C                                                                      C
C     OUTPUT PARAMETER DESCRIPTION...                                  C
C     ===============================                                  C
C                                                                      C
C     PLM...       ARRAY OF INTEGRALS OF FULLY NORMALIZED LEGENDRE     C
C                  FUNCTIONS OF ORDER M. THE INTEGRALS OF FULLY        C
C                  NORMALIZED LEGENDRE FUNCTIONS OF DEGREE L IS STORED C
C                  IN PLM(L+1), BEGINNING WITH L=M.                    C
C                                                                      C
C     COMMON BLOCK DESCRIPTION...                                      C
C     ===========================                                      C
C                                                                      C
C     COMMON/RISR/... SQUARE ROOTS AND INVERSE SQUARE ROOTS.           C
C     THE VARIABLES OF COMMON/RISR/ WILL BE DEFINED BY THE CALL OF     C
C     ROUTINE GEOINL  WITH PARAMETER IENT=1.                           C
C                                                                      C
C     RSR...       ARRAY OF SQUARE ROOTS.                              C
C                  RSR(K) CONTAINS SQRT(K),K=1...2*LMAX+2.             C
C     RISR...      ARRAY OF INVERSE SQUARE ROOTS.                      C
C                  RISR(K) CONTAINS 1./SQRT(K), K=1...2*LMAX+2.        C
C                                                                      C
C     USED ROUTINES...                                                 C
C     ================                                                 C
C                                                                      C
C     SUBROUTINE GEOIPL                          FROM SCU_LIB LIBRARY. C
C                                                                      C
C     NUMERICAL ACCURACY...                                            C
C     =====================                                            C
C                                                                      C
C                                                                      C
C     EXECUTION TIME...                                                C
C     =================                                                C
C                                                                      C
C     ABOUT 0.0025 SEC CPU EXECUTION TIME FOR A CALL OF ROUTINE        C
C     GEOINL WITH IENT=1 AND LMAX=200 ON CDC CYBER 990 OF RRZN         C
C     HANNOVER.                                                        C
C                                                                      C
C     ROUTINE CREATION...    8300510BY H.-G.WENZEL,                    C
C                            INSTITUT FUER ERDMESSUNG,                 C
C                            UNIVERSITAET HANNOVER,                    C
C                            NIENBURGER STR.6,                         C
C                            D-3000 HANNOVER 1,                        C
C                            FEDERAL REPUBLIC OF GERMANY.              C
C                            PHONE...0511-7622796.                     C
C     ROUTINE MODIFICATION...871026 BY H.-G.WENZEL.                    C
C**********************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PLM(501)
      DIMENSION RSR(1002),RISR(1002)
      DIMENSION RIPLL(502)
      COMMON/RISR/ RSR,RISR
c      save
      MAXDEG=500
      IF(IENT.EQ.2) GOTO 100
      IF(IENT.NE.1) GOTO 5000
      IF(NDEG.GT.MAXDEG) GOTO 5010
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IENT=1. STORE SQRT(K) IN DSR(K) AND 1/SQRT(K) IN DISR(K).        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      KM=2*NDEG+2
      IF(KM.LT.15) KM=15
      DO 10 K=1,KM
      FK=FLOAT(K)
      RSR(K)=DSQRT(FK)
C     RSR(K)=SQRT(FK)
      RISR(K)=1./RSR(K)
   10 CONTINUE
      MNEXT=0
      NDEGN=NDEG
      RETURN
  100 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IENT=2.                                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(NDEG.GT.NDEGN) GOTO 5030
      IF(M.NE.0) GOTO 200
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DEFINE START VALUES FOR RECURSION ALGORITM, M=0.                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RCTN=CTN
      RSTN=STN
      RCTS=CTS
      RSTS=STS
      IF(DABS(RCTN).LT.1.D-38) RCTN=0.
      IF(DABS(RSTN).LT.1.D-38) RSTN=0.
      IF(DABS(RCTS).LT.1.D-38) RCTS=0.
      IF(DABS(RSTS).LT.1.D-38) RSTS=0.
C     IF(ABS(RCTN).LT.1.D-38) RCTN=0.
C     IF(ABS(RSTN).LT.1.D-38) RSTN=0.
C     IF(ABS(RCTS).LT.1.D-38) RCTS=0.
C     IF(ABS(RSTS).LT.1.D-38) RSTS=0.
      CALL GEOIPL(IUN6,IUN7,CTN,STN,CTS,STS,NDEG,RIPLL)
      PLM(1)=RIPLL(1)
      RPLM2N=RSR(3)*RCTN
      RPLM2S=RSR(3)*RCTS
      RPMMN=0.
      RPMMS=0.
      RIPLM2=RSR(3)*RISR(4)*(RSTS*RSTS-RSTN*RSTN)
      RPLM1N=RSR(5)*RISR(4)*(RSR(9)*RCTN*RCTN-1.)
      RPLM1S=RSR(5)*RISR(4)*(RSR(9)*RCTS*RCTS-1.)
      RIPLM1=RSR(5)*RISR(4)*(RCTS*RSTS*RSTS-RCTN*RSTN*RSTN)
      PLM(2)=RIPLM2
      PLM(3)=RIPLM1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     RECURSION ALGORITHM FOR M=0.                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RL=3.
      DO 110 L=3,NDEG
      RF1=RSR(2*L+1)/RL
      RF2=(RL-1.)*RISR(2*L-3)
      RPLMN=RF1*(RSR(2*L-1)*RCTN*RPLM1N-RF2*RPLM2N)
      RPLMS=RF1*(RSR(2*L-1)*RCTS*RPLM1S-RF2*RPLM2S)
      RIPLM=RF1*RSR(2*L-1)*(RSTS*RSTS*RPLM1S-RSTN*RSTN*RPLM1N)/(RL+1.)
     1 +(RL-2.)/(RL+1.)*RF1*RF2*RIPLM2
      PLM(L+1)=RIPLM
      RPLM2N=RPLM1N
      RPLM2S=RPLM1S
      RIPLM2=RIPLM1
      RPLM1N=RPLMN
      RPLM1S=RPLMS
      RIPLM1=RIPLM
      RL=RL+1.
  110 CONTINUE
      MNEXT=1
      RETURN
  200 CONTINUE
      IF(M.GT.2) GOTO 230
      IF(M.EQ.2) GOTO 210
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DEFINE START VALUES FOR RECURSION ALGORITHM, M=1.                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RPLM2N=RSR(3)*RSTN
      RPLM2S=RSR(3)*RSTS
      RIPLM2=RIPLL(2)
      RPLM1N=RSR(15)*RCTN*RSTN
      RPLM1S=RSR(15)*RCTS*RSTS
      RIPLM1=RSR(5)*RISR(3)*(RSTS*RSTS*RSTS-RSTN*RSTN*RSTN)
      PLM(2)=RIPLM2
      PLM(3)=RIPLM1
      GOTO 280
  210 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DEFINE START VALUES FOR RECURSION ALGORITHM, M=2.                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RF1=RSR(2*M+1)*RISR(2*M)
      RPLM2N=RF1*RSTN*RPLLN
      RPLM2S=RF1*RSTS*RPLLS
      RIPLM2=RIPLL(3)
      PLM(M+1)=RIPLM2
      RF2=RSR(2*M+3)
      RPLM1N=RF2*RCTN*RPLM2N
      RPLM1S=RF2*RCTS*RPLM2S
      RIPLM1=RSR(2*M+3)*RISR(M+2)*RISR(M+2)*(RSTS*RSTS*RPLM2S-RSTN*RSTN*
     1 RPLM2N)
      PLM(M+2)=RIPLM1
      GOTO 280
  230 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DEFINE START VALUES FOR RECURSION ALGORITHM, M GREATER 2.        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RF1=RSR(2*M+1)*RISR(2*M)
      IF(DABS(RPLLN).LT.1.d-38) RPLLN=0.
      IF(DABS(RPLLS).LT.1.d-38) RPLLS=0.
C     IF(ABS(RPLLN).LT.1.E-100) RPLLN=0.
C     IF(ABS(RPLLS).LT.1.E-100) RPLLS=0.
      RPLM2N=RF1*RSTN*RPLLN
      RPLM2S=RF1*RSTS*RPLLS
      RIPLM2=RIPLL(M+1)
      IF(DABS(RIPLM2).LT.1.d-38) RIPLM2=0.
C     IF(ABS(RIPLM2).LT.1.E-100) RIPLM2=0.
      PLM(M+1)=RIPLM2
      IF(M.EQ.NDEG) GOTO 300
      RF2=RSR(2*M+3)
      RPLM1N=RF2*RCTN*RPLM2N
      RPLM1S=RF2*RCTS*RPLM2S
      IF(DABS(RPLM1N).LT.1.d-38) RPLM1N=0.
      IF(DABS(RPLM1S).LT.1.d-38) RPLM1S=0.
C     IF(ABS(RPLM1N).LT.1.E-100) RPLM1N=0.
C     IF(ABS(RPLM1S).LT.1.E-100) RPLM1S=0.
      RIPLM1=RSR(2*M+3)*RISR(M+2)*RISR(M+2)*(RSTS*RSTS*RPLM2S-RSTN*RSTN*
     1 RPLM2N)
      IF(DABS(RIPLM1).LT.1.d-38) RIPLM1=0.
C     IF(ABS(RIPLM1).LT.1.E-100) RIPLM1=0.
      PLM(M+2)=RIPLM1
      IF(M.EQ.NDEG-1) GOTO 310
  280 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     RECURSION ALGORITHM, M GREATER 0.                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      MP2=M+2
      RPLLN=RPLM2N
      RPLLS=RPLM2S
      RPMM1N=RPMMN
      RPMM1S=RPMMS
      RPMMN=RPLM1N
      RPMMS=RPLM1S
      RL=FLOAT(MP2)
      DO 290 L=MP2,NDEG
      RF1=RSR(2*L+1)*RISR(L+M)*RISR(L-M)
      RF2=RSR(L+M-1)*RSR(L-M-1)*RISR(2*L-3)
      RPLMN=RF1*(RSR(2*L-1)*RCTN*RPLM1N-RF2*RPLM2N)
      RPLMS=RF1*(RSR(2*L-1)*RCTS*RPLM1S-RF2*RPLM2S)
      IF(DABS(RPLMN).LT.1.d-38) RPLMN=0.
      IF(DABS(RPLMS).LT.1.d-38) RPLMS=0.
C     IF(ABS(RPLMN).LT.1.E-100) RPLMN=0.
C     IF(ABS(RPLMS).LT.1.E-100) RPLMS=0.
      RIPLM=(RF1*RSR(2*L-1)*(RSTS*RSTS*RPLM1S-RSTN*RSTN*RPLM1N)+
     1 (RL-2.)*RF1*RF2*RIPLM2)/(RL+1.)
      IF(DABS(RIPLM).LT.1.d-38) RIPLM=0.
C     IF(ABS(RIPLM).LT.1.E-100) RIPLM=0.
      PLM(L+1)=RIPLM
      RPLM2N=RPLM1N
      RPLM2S=RPLM1S
      RIPLM2=RIPLM1
      RPLM1N=RPLMN
      RPLM1S=RPLMS
      RIPLM1=RIPLM
      RL=RL+1.
  290 CONTINUE
      MNEXT=M+1
      RETURN
  300 CONTINUE
      MNEXT=0
      RETURN
  310 CONTINUE
      MNEXT=NDEG
      RPLLN=RPLM2N
      RPLLS=RPLM2S
      RPMM1N=RPMMN
      RPMM1S=RPMMS
      RPMMN=RPLM1N
      RPMMS=RPLM1S
      RETURN
 5000 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     WRONG ENTRY PARAMETER IENT.                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(IUN6,7000) IENT
      write(*,7000) IENT
      STOP
 5010 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     MAXIMUM DEGREE AND ORDER NDEG TOO LARGE.                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(IUN6,7010) NDEG,MAXDEG
      write(*,7010) NDEG,MAXDEG
      STOP
 5020 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     WRONG SEQUENCE OF ORDER M.                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(IUN6,7020) M,MNEXT
      write(*,7020) M,MNEXT
      STOP
 5030 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     MAXIMUM DEGREE AND ORDER CHANGED.                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(IUN6,7030) NDEG,NDEGN,IENT
      write(*,7030) NDEG,NDEGN,IENT
      STOP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FORMAT STATEMENTS.                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 7000 FORMAT(/' ***ERROR IN ROUTINE GEOINL, VERSION 871026 FTN77.'/
     1' ***WRONG ENTRY PARAMETER IENT=',I5/
     2' ***ROUTINE GEOINL STOPS THE EXECUTION.'/)
 7010 FORMAT(/' ***ERROR IN ROUTINE GEOINL, VERSION 871026 FTN77.'/
     1' ***THE WANTED MAXIMUM DEGREE NDEG=',I5,'  EXCEEDS THE ALLOWED MA
     2XIMUM DEGREE MAXDEG=',I5/
     3' ***ROUTINE GEOINL STOPS THE EXECUTION.'/)
 7020 FORMAT(/' ***ERROR IN ROUTINE GEOINL, VERSION 871026 FTN77.'/
     1' ***WRONG SEQUENCE OF ORDER M=',I5,'  WHICH SHOULD BE =',I5/
     2' ***ROUTINE GEOINL STOPS THE EXECUTION.'/)
 7030 FORMAT(/' ***ERROR IN ROUTINE GEOINL, VERSION 871026 FTN77.'/
     1' ***THE MAXIMUM DEGREE AND ORDER NDEG=',I5,' HAS CHANGED FROM NDE
     2GN=',I5,' DURING THE RECURSION LOOP.'/
     3' ***THIS IS NOT ALLOWED FOR ENTRY PARAMETER IENT=',I5/
     4' ***ROUTINE GEOINL STOPS THE EXECUTION.'/)
      END
      SUBROUTINE GEOIPL(IUN6,IUN7,RCTN,RSTN,RCTS,RSTS,LMAX,RIPLL)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     ROUTINE GEOIPL, VERSION 871026 FORTRAN 77.                       C
C                                                                      C
C     SOURCE CODE LIBRARY...SCU_LIB                                    C
C                                                                      C
C     THE ROUTINE GEOIPL COMPUTES THE INTEGRALS OF FULLY               C
C     NORMALIZED LEGENDRE FUNCTIONS OF DEGREE AND ORDER L.             C
C                                                                      C
C     REFERENCE... WENZEL,H.-G. 1985.. ZUR ANWENDUNG UND BERECHNUNG    C
C                  VON HOCHAUFLOESENDEN KUGELFUNKTIONSMODELLEN FUER    C
C                  DAS GRAVITATIONSPOTENTIAL DER ERDE.                 C
C                  WISSENSCHAFTLICHE ARBEITEN DER FACHRICHTUNG         C
C                  VERMESSUNGSWESEN DER UNIVERSITAET HANNOVER,         C
C                  HANNOVER 1985.                                      C
C                                                                      C
C     INPUT PARAMETER DESCRIPTION...                                   C
C     ==============================                                   C
C                                                                      C
C     IUN6...      FORMATTED LINE PRINTER UNIT.                        C
C     IUN7...      FORMATTED CONSOLE UNIT (DATA TERMINAL SCREEN).      C
C     RCTN...      COS OF POLAR DISTANCE OF NORTHERN BORDER.           C
C     RSTN...      SIN OF POLAR DISTANCE OF NORTHERN BORDER.           C
C     RCTS...      COS OF POLAR DISTANCE OF SOUTHERN BORDER.           C
C     RSTS...      SIN OF POLAR DISTANCE OF SOUTHERN BORDER.           C
C     LMAX...      MAXIMUM DEGREE AND ORDER, UP TO WHICH AND           C
C                  INCLUSIVELY THE INTEGRALS OF FULLY NORMALIZED       C
C                  LEGENDRE FUNCTIONS OF DEGREE AND ORDER L WILL BE    C
C                  COMPUTED.                                           C
C                                                                      C
C     OUTPUT PARAMETER DESCRIPTION...                                  C
C     ===============================                                  C
C                                                                      C
C     RIPLL...     ARRAY OF THE INTEGRALS OF FULLY NORMALIZED          C
C                  LEGENDRE FUNCTIONS OF DEGREE AND ORDER L.           C
C                  THE INTEGRAL OF FULLY NORMALIZED LEGENDRE FUNCTION  C
C                  OF DEGREE AND ORDER L IS STORED IN RIPLL(L+1),      C
C                  BEGINNING WITH L=0.                                 C
C                                                                      C
C     COMMON BLOCK DESCRIPTION...                                      C
C     ===========================                                      C
C                                                                      C
C     COMMON/RISR/... SQUARE ROOTS AND INVERSE SQUARE ROOTS.           C
C     THE VARIABLES OF COMMON/RISR/ HAVE TO BE DEFINED BEFORE THE      C
C     CALL OF ROUTINE GEOIPL, BY E.G. CALLING ROUTINE GEOINL.          C
C                                                                      C
C     RSR...       ARRAY OF SQUARE ROOTS.                              C
C                  RSR(K) CONTAINS SQRT(K),K=1...2*LMAX+2.             C
C     RISR...      ARRAY OF INVERSE SQUARE ROOTS.                      C
C                  RISR(K) CONTAINS 1./SQRT(K), K=1...2*LMAX+2.        C
C                                                                      C
C     USED ROUTINES... NONE.                                           C
C     ================                                                 C
C                                                                      C
C     NUMERICAL ACCURACY...                                            C
C     =====================                                            C
C                                                                      C
C     LESS OR EQUAL 1*10**-12 RELATIVE ERROR FOR THE COMPUTED          C
C     INTEGRALS. THE NUMERICAL ACCURACY OF THE COMPLETE SET OF         C
C     COMPUTED INTEGRALS IS CHECKED BY COMPARISON OF THE LAST INTEGRAL C
C     COMPUTED BY DOWNWARD RECURRENCE RELATIONS WITH DIRECT COMPUTED   C
C     INTEGRAL. IN CASE OF RELATIVE ERROR EXCEEDING 1*10**-11, A       C
C     WARNING WILL BE WRITTEN ON LINE PRINTER OUTPUT UNIT IUN6.        C
C                                                                      C
C     EXECUTION TIME...                                                C
C     =================                                                C
C                                                                      C
C     0.040 SEC CPU TIME ON CDC CYBER 76 OF RRZN HANNOVER FOR          C
C     LMAX=180.                                                        C
C                                                                      C
C     ROUTINE CREATION...    830601 BY H.-G.WENZEL,                    C
C                            INSTITUT FUER ERDMESSUNG,                 C
C                            UNIVERSITAET HANNOVER,                    C
C                            NIENBURGER STR.6,                         C
C                            D3000 HANNOVER 1,                         C
C                            FEDERAL REPUBLIC OF GERMANY.              C
C                            PHONE...0511-7622796.                     C
C     ROUTINE MODIFICATION...871026 BY H.-G.WENZEL.                    C
C**********************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RSR(1002),RISR(1002)
      DIMENSION RIPLL(501),RCTNI(502),RSTNI(502),RSTSI(502)
      EQUIVALENCE (RCTNI(1),RSTSI(1))
      COMMON/RISR/ RSR,RISR
c      save
      RSDDT=RSTS*RCTN-RCTS*RSTN
      RSDDT2=RSDDT*RSDDT
      LM=LMAX
      LMP1=LM+1
      LMP2=LM+2
      IF(LMP1.GT.501) GOTO 5000
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     INITIALIZATION OF ARRAYS RCTNI AND RSTNI.                        C
C     RCTNI(I+1) CONTAINS COS(THETA NORTH)**I.                         C
C     RSTNI(I+1) CONTAINS SIN(THETA NORTH)**I.                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RCTNI(1)=1.
      RSTNI(1)=1.
      DO 100 I=1,LMP1
      IP1=I+1
      RCTNI(IP1)=RCTNI(I)*RCTN
      RSTNI(IP1)=RSTNI(I)*RSTN
      IF(DABS(RCTNI(IP1)).LT.1.d-38) RCTNI(IP1)=0.
      IF(DABS(RSTNI(IP1)).LT.1.d-38) RSTNI(IP1)=0.
C     IF(ABS(RCTNI(IP1)).LT.1.E-100) RCTNI(IP1)=0.
C     IF(ABS(RSTNI(IP1)).LT.1.E-100) RSTNI(IP1)=0.
  100 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     COMPUTE THE INTEGRAL OF SIN(THETA)**(L+1) FOR L=LMAX-1 AND       C
C     L=LMAX.                                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 200 I=1,2
      L=LMAX-2+I
      LP1=L+1
      LP2=L+2
      RLP2=FLOAT(LP2)
      RLP1AK=1.
      RZZ=RLP2
      RNN=0.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     LOOP FOR K. START WITH K=L+1.                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RK=RLP2
      RISL=0.
      RISL=0.
      RSDDTK=1.
      DO 210 J=1,LP2
      K=LP2-J
      KP1=K+1
      RK=RK-1.
      RLP2MK=RLP2-RK
      RSERK=RCTNI(J)*RSTNI(KP1)*RLP1AK
      IF(RSERK.EQ.0.) GOTO 230
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     LOOP FOR ITLK. ITLK IS THE INTEGRAL OF SIN(ALPHA)**(L+1) FOR     C
C     ALPHA=0...THETA NORTH - THETA SOUTH.                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RSDI=RSDDTK*RSDDT
      RITLK=RSDI/RLP2MK
      RF=1.
      RZ=RK-1.
      RN=2.
      IT=0
  220 IT=IT+1
      RF=-RF*RZ/RN
      RSDI=RSDI*RSDDT2
      RLP2MK=RLP2MK+2.
      RSER=RF*RSDI/RLP2MK
      RITLK=RITLK+RSER
      IF(ABS(RSER).LT.1.e-30) GOTO 225
C
      RZ=RZ-2.
      RN=RN+2.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CUT OFF OF SERIES DEVELOPMENT FOR ITLK.                          C
C     (1.E-12 RELATIVE ACCURACY).                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RERR=RSER*1.E12
      IF(ABS(RERR).GT.ABS(RITLK)) GOTO 220
      RISL=RISL+RITLK*RSERK
  225 CONTINUE
  230 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     UPDATE L+1 ABOVE K.                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RZZ=RZZ-1.
      RNN=RNN+1.
      RLP1AK=RLP1AK*RZZ/RNN
      RSDDTK=RSDDTK*RSDDT
  210 CONTINUE
      RIPLL(LP1)=RISL
  200 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     COMPUTE RISL(L) BY DOWNWARD RECURRENCE RELATION.                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RSTSI(1)=1.
      DO 300 I=1,LM
      IP1=I+1
      RSTSI(IP1)=RSTSI(I)*RSTS
      IF(DABS(RSTSI(IP1)).LT.1.d-38) RSTSI(IP1)=0.
C     IF(ABS(RSTSI(IP1)).LT.1.E-100) RSTSI(IP1)=0.
  300 CONTINUE
      RZ=FLOAT(LMP2)
      RN=RZ-1.
      DO 400 K=3,LMP1
      L=LMP1-K
      LP1=L+1
      LP3=LP1+2
      RZ=RZ-1.
      RN=RN-1.
      RIPLL(LP1)=(RIPLL(LP3)*RZ+(RSTSI(LP3)*RCTS-RSTNI(LP3)*RCTN))/RN
  400 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NORMALIZATION OF THE INTEGRALS OF  SIN(THETA)**(L+1).            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IZ=1
      IN=0
      RFAC=RSR(2)
      DO 500 L=1,LM
      IZ=IZ+2
      IN=IN+2
      RFAC=RFAC*RSR(IZ)*RISR(IN)
      LP1=L+1
      RIPLL(LP1)=RIPLL(LP1)*RFAC
  500 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CONTROL OF ACCURACY BY DIRECT COMPUTATION OF IP00.               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RIP00=RCTN-RCTS
      RERR=RIPLL(1)-RIP00
      RELER=RERR/RIP00
      RELER2=RELER*RELER
      IF(RELER2.GT.1.E-20) GOTO 5010
      RETURN
 5000 CONTINUE
      WRITE(IUN6,7000) LMAX
      write(*,7000) LMAX
      STOP
 5010 CONTINUE
      WRITE(IUN6,7010) RELER
      write(*,7010) RELER
      RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FORMAT STATEMENTS.                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 7000 FORMAT(/' ***ERROR IN ROUTINE GEOIPL, VERSION 871026 FTN77.'/
     1' ***THE WANTED MAXIMUM DEGREE LMAX=',I5,'  EXCEEDS THE ALLOWED MA
     2XIMUM DEGREE OF 500.'/
     3' ***ROUTINE GEOIPL STOPS THE EXECUTION.'/)
 7010 FORMAT(' ***WARNING FROM ROUTINE GEOIPL, VERSION 871026 FTN77.'/
     1/' ***RELATIVE ERROR OF IP00 IS',E15.5,'  GREATER THAN 1.E-10.')
      END
      SUBROUTINE GEOPMI(IUN6,IUN7,IUNI,IUN10,IUN20,IUN30,IUN31,IUN32,
     1 IUN33,LMIN,LMAX,RDLAT,RDLON,RPLIM,ITMAX)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     ROUTINE GEOPMI, VERSION 871026 FORTRAN 77.                       C
C                                                                      C
C     COMPUTATION OF NEW OR CORRECTION OF EXISTING SPHERICAL HARMONIC  C
C     POTENTIAL COEFFICIENTS BY MEANS OF MEAN FREE AIR ANOMALIES       C
C     USING INTEGRAL FORMULAS WITHIN AN ITERATIV ALGORITHM.            C
C     MEAN FREE AIR ANOMALIES DEFINED IN EQUAL ANGULAR BLOCKS CAN      C
C     BE USED ONLY. A COMPARISON IS CARRIED OUT AFTER EACH ITERATION   C
C     BETWEEN OBSERED MEAN ANOMALIES AND ANOMALIES FROM COMPUTED       C
C     (CORRECTED) SPHERICAL HARMONIC COEFFICIENTS.                     C
C                                                                      C
C     INPUT PARAMETER DESCRIPTION...                                   C
C     ==============================                                   C
C                                                                      C
C     IUN6...    FORMATTED LINE PRINTER UNIT NUMBER.                   C
C     IUN7...    FORMATTED CONSOLE UNIT NUMBER (DATA TERMINAL SCREEN). C
C     IUNI...    FORMATTED UNIT FOR INPUT OF MEAN GRAVITY ANOMALIES    C
C                TO BE USED (E.G. GLWENZX.DATMFA.F3030UCT86).          C
C                IN ORDER TO SAVE COMPUTATION TIME, THE ANOMALIES      C
C                SHOULD BE SORTED TO DECREASING LATITUDE AND           C
C                WITHIN SAME LATITUDE TO INCREASING LONGITUDE.         C
C                THE MEAN FREE AIR GRAVITY ANOMALIES ARE EXPECTED      C
C                TO REFER TO GRS 1980 REFERENCE SYSTEM WITH            C
C                ATMOSPHERIC CORRECTION NOT APPLIED (WILL BE APPLIED   C
C                DURING EXECUTION OF ROUTINE GEOPMI).                  C
C                THE MEAN ELEVATIONS TO BE GIVEN WITH THE ANOMALIES    C
C                ARE EXPECTED TO DESCRIBE THE ELEVATION, TO WHICH      C
C                THE ANOMALIES ARE REFERRING (I.E. ON THE OCEANS       C
C                ELEVATIONS SHOULD BE ZERO). IF ELEVATIONS ARE         C
C                UNDEFINED, THEY SHOULD BE SET TO 9999.99              C
C     IUN10...   FORMATTED INPUT UNIT FOR EXISTING POTENTIAL           C
C                COEFFICIENTS WHICH HAVE TO BE CORRECTED.              C
C                THERE IS NO NEED TO ATTACH IUN10 IN CASE OF NEW       C
C                POTENTIAL COEFFICIENT CALCULATIONS BUT IT HAS TO BE   C
C                OPEND BEFORE CALLING ROUTINE GEOPMI JUST IN ANY CASE. C
C     IUN20...   FORMATTED OUTPUT UNIT FOR NEW/CORRECTED POTENTIAL     C
C                MODEL COEFFICIENTS. IF THERE HAS BEEN USED A START    C
C                MODEL, THE GIVEN STANDARD DEVIATIONS OF THE START     C
C                MODEL WILL ALSO BE USED AS STANDARD DEVIATIONS OF     C
C                THE NEW/CORRECTED MODEL, BECAUSE ROUTINE GEOPMI IS    C
C                NOT ABLE TO ESTIMATE STANDARD DEVIATIONS OF THE       C
C                COEFFICIENTS.                                         C
C     IUN30...   SEQUENTIAL UNFORMATTED SCRATCH TAPE, ON WHICH THE     C
C                COEFFICIENTS OF THE START MODEL WILL BE STORED        C
C                DURING THE EXECUTION OF ROUTINE GEOPMI.               C
C     IUN31...   SEQUENTIAL UNFORMATTED SCRATCH TAPE, ON WHICH THE     C
C                ANOMALIES WILL BE STORED DURING THE EXECUTION OF      C
C                ROUTINE GEOPMI.                                       C
C     IUN32...   SEQUENTIAL UNFORMATTED SCRATCH TAPE, ON WHICH THE     C
C                ANOMALIES WILL BE STORED DURING THE EXECUTION OF      C
C                ROUTINE GEOPMI.                                       C
C     IUN33...   SEQUENTIAL UNFORMATTED SCRATCH TAPE, ON WHICH THE     C
C                COEFFICIENTS AND STANDARD DEVIATIONS OF THE START     C
C                MODEL WILL BE STORED DURING THE EXECUTION OF ROUTINE  C
C                GEOPMI (IF THERE IS A START MODEL AVAILABLE).         C
C     LMIN...    THE COEFFICIENTS OF THE NEW MODEL WILL NOT BE         C
C                COMPUTED BELOW DEGREE LMIN. IF THERE EXIST A START    C
C                MODEL, THE COEFFICIENTS OF THE START MODEL WILL USED  C
C                BELOW DEGREE LMIN.                                    C
C                THE PARAMETER LMIN MAY BE USED TO RESTRICT THE        C
C                DETERMINATION OF CORRECTIONS TO THE COEFFICIENTS      C
C                TO MEDIUM WAVELENGTH PART, DEPENDING ON THE SIZE OF   C
C                THE AREA WHERE MEAN ANOMALIES ARE GIVEN (E.G.         C
C                LMIN = 180/SIZE IN DEGREE).                           C
C     LMAX...    THE COEFFICIENTS OF THE NEW MODEL WILL BE COMPUTED    C
C                UP TO DEGREE AND ORDER LMAX.                          C
C     RDLAT...   BLOCK SIZE OF INPUT ANOMALIES IN LATITUDE  IN DEGREE. C
C     RDLON...   BLOCK SIZE OF INPUT ANOMALIES IN LONGITUDE IN DEGREE. C
C     RPLIM...   LIST LIMIT FOR DIFFERENCES BETWEEN INPUT AND MODEL    C
C                ANOMALIES IN MGAL. ONLY THOSE ANOMALIES WILL BE       C
C                PRINTED, WHERE THE DIFFERENCE BETWEEN INPUT AND       C
C                MODEL ANOMALIES EXCEED RPLIM.                         C
C     ITMAX...   NUMBER OF ITERATIONS.                                 C
C                                                                      C
C     USED ROUTINES...                                                 C
C     ================                                                 C
C                                                                      C
C     SUBROUTINE GEOBET                          FROM SCU_LIB LIBRARY. C
C     SUBROUTINE GEOCTR                          FROM SCU_LIB LIBRARY. C
C     SUBROUTINE GEOHIS                          FROM SCU_LIB LIBRARY. C
C     SUBROUTINE GEOICS                          FROM SCU_LIB LIBRARY. C
C     SUBROUTINE GEOINL                          FROM SCU_LIB LIBRARY. C
C     SUBROUTINE GEOIPL                          FROM SCU_LIB LIBRARY. C
C     SUBROUTINE GEOREF                          FROM SCU_LIB LIBRARY. C
C                                                                      C
C     EXECUTION TIME...                                                C
C     =================                                                C
C                                                                      C
C     USING 1368 30 BY 30 MIN MEAN ANOMALIES, MAXIT = 5, AND           C
C     COMPILER OPTIMIZATION LEVEL OL=HIGH, THE TOTAL CPU EXECUTION     C
C     TIME IS                                                          C
C                                                                      C
C      12.84 SEC FOR LMAX =  30 USING GEM-L2 AS START MODEL,           C
C     278.34 SEC FOR LMAX = 200 USING GPM-2  AS START MODEL            C
C                                                                      C
C     ON CDC CYBER 990 OF RRZN HANNOVER.                               C
C                                                                      C
C     NOTE THAT...  THIS VERSION OF ROUTINE GEOPMI IS LIMITED TO A     C
C     ============  MAXIMUM DEGREE AND ORDER 360 (DEPENDING ON         C
C     PARAMETER IP1 USED IN DIMENSION STATEMENTS.                      C
C                                                                      C
C     ROUTINE CREATION...    860925 BY G. WEBER,                       C
C                            INSTITUT FUER GEODAESIE                   C
C                            UND PHOTOGRAMMETRIE                       C
C                            TECHNISCHE UNIVERSITAET BERLIN            C
C                            STRASSE DES 17. JUNI 135                  C
C                            D-1000 BERLIN 12                          C
C                            FEDERAL REPUBLIC OF GERMANY               C
C     ROUTINE MODIFICATION...871026 BY H.-G. WENZEL,                   C
C                            INSTITUT FUER ERDMESSUNG,                 C
C                            UNIVERSITAET HANNOVER,                    C
C                            NIENBURGER STR. 6,                        C
C                            D-3000 HANNOVER 1,                        C
C                            FEDERAL REPUBLIC OF GERMANY               C
C**********************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical lform
c     real*8 mth$cvt_d_g
      common /form/lform
      CHARACTER CHIS(20)*10
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IP1=20301 FOR LMAX = 200, 65341 FOR LMAX = 360.                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER (IP1=65341)
      DIMENSION XMIN(20),XMAX(20),RDX(20)
      CHARACTER CNAME*10,CREFSN*10
      CHARACTER CTEXT(8)*10,CENDT*10,CFILE*10,CENDF*10
      DIMENSION RCLMO(IP1),RSLMO(IP1),RCLMN(IP1),RSLMN(IP1),PLMSTO(IP1)
      DIMENSION GCMA(501),GSMA(501),GCMAR(501),GSMAR(501)
      DIMENSION RCOLD(501),RSOLD(501),SCOLD(501),SSOLD(501)
      EQUIVALENCE (GCMA(1),RCOLD(1)),(GSMA(1),RSOLD(1))
      EQUIVALENCE (GCMAR(1),SCOLD(1)),(GSMAR(1),SSOLD(1))
      DIMENSION CMI(501),SMI(501),PLM(501),FACTL(501)
      DIMENSION ADVO(501),ADVN(501),ADVD(501)
      COMMON/REFELL/ GM,A,F,OM,RMEAN,GAMEAN,GE,RK,E2,ES2,U0,CN(11)
      DATA CREFSN /'GRS 1980  '/,CENDT/'C*********'/,CENDF/'FILE END  '/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DEFINITION OF CONSTANT PARAMETERS AND DEFINITION OF              C
C     LOGICAL UNIT NUMBERS FOR FIRST ITERATION.                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      MAXDI2=360
      MAXDI2=((MAXDI2+1)*MAXDI2)/2+(MAXDI2+1)
      LDUMMY=99999
      RLIMIT=1.E-14
      PI=4.*datan(1.d0)
C     PI=4.*ATAN(1.)
      RAD=PI/180.
      ITER=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     WRITE FORMAL PARAMETERS, GENERAL INITIALIZATIONS.                C
C     MAXCS IS MAXIMUM SUBSCRIPT OF COEFFICIENTS RCLMO, RSLMO (START   C
C     MODEL RESP. MODEL OF LAST ITERATION) AND RCLMN, RSLMN            C
C     (ACTUAL MODEL DURING ITERATION RESP. NEW MODEL AFTER ITERATION). C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(IUN6,7000) IUNI,LMAX
      write(*,7000) IUNI,LMAX
      LPMAX=11
      MINP1=LMIN+1
      MAXP1=LMAX+1
      MAXP12=MAXP1*2
      MAXCS=MAXP1*(MAXP1-1)/2+MAXP1
      IF(MAXCS.LE.MAXDI2) GOTO 100
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     MAXIMUM DEGREE TOO LARGE.                                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(IUN6,7011) MAXDI2
      write(*,7011) MAXDI2
      GOTO 4000
  100 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     READ MEAN FREE AIR GRAVITY ANOMALIES FROM FORMATTED UNIT IUNI    C
C     (5 OR 1) AND STORE THEM ON UNFORMATTED UNIT IUN31.               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  120 continue
c      READ(IUNI,7015) (CTEXT(J),J=1,8)
c      WRITE(IUN6,7016) (CTEXT(J),J=1,8)
c      write(*,7016) (CTEXT(J),J=1,8)
c      IF(CTEXT(1).NE.CENDT) GOTO 120
c
c  modified input format of gravity data
c
      kk = 0
      fmin = 9999.99
      fmax = -9999.99
  130 READ(IUNI,*,END=140) istat,GLAT,ELON,ELEV,FREAIR
c
      kk = kk+1
      if (freair.lt.fmin) fmin = freair
      if (freair.gt.fmax) fmax = freair
      FMOD=0.
      GLATN=GLAT+RDLAT*0.5
      GLATS=GLAT-RDLAT*0.5
      ELONW=ELON-RDLON*0.5
      ELONE=ELON+RDLON*0.5
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SET UNDEFINED ELEVATION TO ZERO.                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(ELEV.GT.9999.) ELEV=0.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     COMPUTE ATMOSPHERIC CORRECTION.                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FREAIR=FREAIR+0.874-9.9E-5*ELEV+3.5625E-9*ELEV*ELEV
      WRITE(IUN31) GLATN,GLATS,ELONW,ELONE,ELEV,FREAIR,FMOD
      GOTO 130
  140 CONTINUE
      write(*,141) kk,fmin,fmax
  141 format(/' mean anomaly data input: number, min, max =',i6,2f9.2)
      IUNIA=IUN31
      IUNOA=IUN32
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DEFINE OF GEODETIC REFERENCE SYSTEM TO BE USED.                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IPRGRS=1
      CALL GEOREF(IUN6,IUN7,CREFSN,IPRGRS)
      CONST=1.E-5/(GM*4.*PI)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     INITIALIZATION OF ROUTINE ICSM FOR LONGITUDE DIFFERENCE RDLON.   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IENICS=1
      RDLONR=RDLON*RAD
      CALL GEOICS(IUN6,IUN7,IENICS,RDLONR,CLON,SLON,LMAX,CMI,SMI)
      IENICS=2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     INITIALIZATION OF ROUTINE GEOIPL.                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IENINL=1
      CALL GEOINL(IUN6,IUN7,IENINL,M,CTN,STN,CTS,STS,LMAX,PLM)
      IENINL=2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     INITIALIZATION OF COEFFICIENTS ARRAYS RCLMO, SLM, RCLMN AND      C
C     RSLMN. RCLMO AND RSLMO WILL HOLD START MODEL RESP. MODEL         C
C     OF LAST ITERATION.                                               C
C     RCLMN AND RSLMN WILL HOLD ACTUAL MODEL OF THE ITERATION.         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 150 MP1=1,MAXP1
      ADVO(MP1)=0.
      DO 150 LP1=MP1,MAXP1
      LP1MP1=(MP1-1)*(MAXP12-MP1)/2+LP1
      RCLMO(LP1MP1)=0.
      RSLMO(LP1MP1)=0.
      RCLMN(LP1MP1)=0.
  150 RSLMN(LP1MP1)=0.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     READ POTENTIAL COEFFICIENTS OF OLD MODEL.                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ISTAMO=0
c      READ(IUN10,7015,END=250) (CTEXT(J),J=1,8)
c      READ(IUN10,7017) CNAME,MAXO,MAXOO
c      WRITE(IUN6,7018) CNAME,MAXO,MAXOO
c      write(*,7018) CNAME,MAXO,MAXOO
c  160 continue
c      READ(IUN10,7015) (CTEXT(J),J=1,8)
c      WRITE(IUN6,7016) (CTEXT(J),J=1,8)
c      IF(CTEXT(1).NE.CENDT) GOTO 160
c
c  modified reading of potential coefficients - binary fornat
c
      kk = 0
  170 continue
      if (lform) read(iun10,*,end=180) l,m,cii,sii
      IF (l .le. 3) PRINT*,l,m,cii,sii
1701  format(2i4,2d16.9)
c
      if (.not.lform) READ(IUN10,end=180) L,M,CII,SII
c
c  convert to g-floating on vax
c
c     if (.not.lform) cii = mth$cvt_d_g(cii)
c     if (.not.lform) sii = mth$cvt_d_g(sii)
c
      stdcii = 0.0
      stdsii = 0.0
c
      IF(L.GT.1000) GOTO 180
      IF(L.GT.LMAX.OR.M.GT.LMAX) GOTO 170
      if (m.eq.0) kk = kk+1
      if (m.gt.0) kk = kk+2
      LP1=L+1
      MP1=M+1
      LP1MP1=(MP1-1)*(MAXP12-MP1)/2+LP1
      RCLMO(LP1MP1)=CII
      RSLMO(LP1MP1)=SII
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     STANDARD DEVIATIONS WILL BE STORED IN RCLMN, RSLMN, UNITL THE    C
C     START MODEL WILL BE TRANSFERRED TO IUN33.                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RCLMN(LP1MP1)=STDCII
      RSLMN(LP1MP1)=STDSII
      GOTO 170
  180 write(*,181) kk
  181 format(' number of potential coefficients input: ',i6)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ISTAMO=0 IF NO START MODEL IS USED, =1 IF A START MODEL IS USED. C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ISTAMO=1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SUBTRACT NORMAL POTENTIAL.                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 190 LP1=1,LPMAX
  190 RCLMO(LP1)=RCLMO(LP1)-CN(LP1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     STORE START MODEL COEFFICIENTS AND STANDARD DEVIATIONS ON        C
C     SCRATCH UNIT IUN33. (START MODEL COEFFICIENTS ARE REDUCED BY     C
C     NORMAL POTENTIAL).                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REWIND IUN33
      DO 200 MP1=1,MAXP1
      KSA=(MP1-1)*(MAXP12-MP1)/2+MP1
      KSE=(MP1-1)*(MAXP12-MP1)/2+MAXP1
      WRITE(IUN33) (RCLMO(KS),KS=KSA,KSE)
      WRITE(IUN33) (RSLMO(KS),KS=KSA,KSE)
      WRITE(IUN33) (RCLMN(KS),KS=KSA,KSE)
      WRITE(IUN33) (RSLMN(KS),KS=KSA,KSE)
  200 CONTINUE
      REWIND IUN33
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     STORE START MODEL COEFFICIENTS IN RCLMN, RSLMN.                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      KS=0
      DO 210 MP1=1,MAXP1
      DO 210 LP1=MP1,MAXP1
      KS=KS+1
      RCLMN(KS)=RCLMO(KS)
  210 RSLMN(KS)=RSLMO(KS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     INITIALIZATIONS OF ROUTINE GEOHIS.                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  250 CONTINUE
      XMIN(1)=-100.
      XMIN(2)=-100.
      XMIN(3)=-100.
      XMAX(1)= 100.
      XMAX(2)= 100.
      XMAX(3)= 100.
      RDX(1)= 10.
      RDX(2)= 10.
      RDX(3)= 10.
      CHIS(1)='OBS.  FAGA'
      CHIS(2)='MODEL FAGA'
      CHIS(3)='OBS.-MODEL'
      NEWP=0
      IENHIS=1
      CALL GEOHIS(IUN6,IUN7,IENHIS,CHIS,XMIN,XMAX,RDX,1,FREAIR,NEWP)
      CALL GEOHIS(IUN6,IUN7,IENHIS,CHIS,XMIN,XMAX,RDX,2,FMOD,NEWP)
      CALL GEOHIS(IUN6,IUN7,IENHIS,CHIS,XMIN,XMAX,RDX,3,FDIFF,NEWP)
      IENHIS=2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NEW ITERATION. INITIALIZATIONS FOR FIRST OBSERVATION.            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 1000 WRITE(IUN6,7029) ITER
      write(*,7029) ITER
      WRITE(IUN6,7027) RPLIM
c      write(*,7027) RPLIM
      NUMBT=0
      GLATO=99.
      RMSDIF=0.
      REWIND IUNIA
      REWIND IUNOA
      IEND=0
 1010 READ(IUNIA,END=2000) GLATN,GLATS,ELONW,ELONE,ELEV,FREAIR,FMOD
      NUMBT=NUMBT+1
      IF(ELONE.LT.ELONW) ELONE=ELONE+360.
      ELON=(ELONE+ELONW)*0.5
      IF(ELON.GE.360.) ELON=ELON-360.
      CLON=DCOS(ELON*RAD)
      SLON=DSIN(ELON*RAD)
C     CLON=COS(ELON*RAD)
C     SLON=SIN(ELON*RAD)
      IF(DABS(GLATN-GLATO).LT.RLIMIT) GOTO 1100
C     IF(ABS(GLATN-GLATO).LT.RLIMIT) GOTO 1100
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NEW LATITUDE.                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      GLATO=GLATN
      GLAT=0.5*(GLATN+GLATS)
      H=0.
      KCTR=1
      CLATN=DCOS(GLATN*RAD)
      SLATN=DSIN(GLATN*RAD)
C     CLATN=COS(GLATN*RAD)
C     SLATN=SIN(GLATN*RAD)
      CALL GEOCTR(IUN6,IUN7,KCTR,CLATN,SLATN,CLON,SLON,H,CTN,STN,RN,XN,
     1 YN,ZN)
      CLATS=DCOS(GLATS*RAD)
      SLATS=DSIN(GLATS*RAD)
C     CLATS=COS(GLATS*RAD)
C     SLATS=SIN(GLATS*RAD)
      CALL GEOCTR(IUN6,IUN7,KCTR,CLATS,SLATS,CLON,SLON,H,CTS,STS,RS,XS,
     1 YS,ZS)
      R0=0.5*(RN+RS)
      RDSIG=(CTN-CTS)*RDLON*RAD
      CONST2=GM/(R0*R0*RDSIG)*1.E5
      CONST3=GM/(R0*R0*R0*RDSIG)*1.E5
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     STORE FACTL.                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FACTL(1)=0.
      FACTL(2)=0.
      R0IA=R0/A
      R0IAL=R0IA
      DO 1020 LP1=3,MAXP1
      FACTL(LP1)=0.
      L=LP1-1
      RLM1=FLOAT(L-1)
      CALL GEOBET(IUN6,IUN7,RDSIG,L,BETAL)
      R0IAL=R0IAL*R0IA
      IF(LP1.LT.MINP1) GOTO 1020
      FACTL(LP1)=CONST*R0*R0*R0IAL/(RLM1*BETAL)
 1020 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     COMPUTE AND STORE INTEGRALS OVER LEGENDRE FUNCTIONS.             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      KS=0
      DO 1030 MP1=1,MAXP1
      M=MP1-1
      CALL GEOINL(IUN6,IUN7,IENINL,M,CTN,STN,CTS,STS,LMAX,PLM)
      DO 1030 LP1=MP1,MAXP1
      KS=KS+1
 1030 PLMSTO(KS)=PLM(LP1)
      IF(ITER.EQ.0.AND.ISTAMO.EQ.0) GOTO 1100
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     STORE VECTORS GCMA, GSMA, GCMAR,GSMAR.                           C
C     GCMA AND GSMA ARE USED FOR THE COMPUTATION OF MEAN ANOMALIES AT  C
C     ZERO ELEVATION, GCMAR AND GSMAR USE USED FOR THE COMPUTATION OF  C
C     RADIAL DERIVATIVES OF MEAN ANOMALIES.                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      M=0
      MP1=M+1
      GCMA(MP1)=0.
      GSMA(MP1)=0.
      GCMAR(MP1)=0.
      GSMAR(MP1)=0.
      KS=0
      AIRM=1.
      AIRL=AIRM
      AIR=A/R0
      DO 1040 LP1=1,MAXP1
      KS=KS+1
      L=LP1-1
      RLM1=FLOAT(L-1)
      RLP2=FLOAT(L+2)
      GCMA(MP1)=GCMA(MP1)+AIRL*RLM1*RCLMO(KS)*PLMSTO(KS)
      GCMAR(MP1)=GCMAR(MP1)+AIRL*RLM1*RLP2*RCLMO(KS)*PLMSTO(KS)
      AIRL=AIRL*AIR
 1040 CONTINUE
      DO 1050 MP1=2,MAXP1
      AIRM=AIRM*AIR
      AIRL=AIRM
      GCMA(MP1)=0.
      GSMA(MP1)=0.
      GCMAR(MP1)=0.
      GSMAR(MP1)=0.
      DO 1050 LP1=MP1,MAXP1
      KS=KS+1
      RLM1=FLOAT(LP1-2)
      RLP2=RLM1+3.
      GCMA(MP1)=GCMA(MP1)+AIRL*RLM1*RCLMO(KS)*PLMSTO(KS)
      GSMA(MP1)=GSMA(MP1)+AIRL*RLM1*RSLMO(KS)*PLMSTO(KS)
      GCMAR(MP1)=GCMAR(MP1)+AIRL*RLM1*RLP2*RCLMO(KS)*PLMSTO(KS)
      GSMAR(MP1)=GSMAR(MP1)+AIRL*RLM1*RLP2*RSLMO(KS)*PLMSTO(KS)
      AIRL=AIRL*AIR
 1050 CONTINUE
 1100 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NEW LONGITUDE.                                                   C
C     COMPUTE MEAN ANOMALY FROM EXISTING MODEL USING INTEGRALS.        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FMOD=0.
      FMODR=0.
      CALL GEOICS(IUN6,IUN7,IENICS,RDLONR,CLON,SLON,LMAX,CMI,SMI)
      IF(ITER.EQ.0.AND.ISTAMO.EQ.0) GOTO 1115
      DO 1110 MP1=1,MAXP1
      FMOD=FMOD+CMI(MP1)*GCMA(MP1)+SMI(MP1)*GSMA(MP1)
      FMODR=FMODR-CMI(MP1)*GCMAR(MP1)-SMI(MP1)*GSMAR(MP1)
 1110 CONTINUE
      FMOD=FMOD*CONST2
      FMODR=FMODR*CONST3
      FMOD=FMOD+ELEV*FMODR
 1115 FDIFF=FREAIR-FMOD
      RMSDIF=RMSDIF+FDIFF**2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PRINT FREE AIR GRAVITY ANOMALY GRADIENT IN MGAL/KM.              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FMODR=FMODR*1000.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PRINT COMPARISON OF OBSERVED WITH COMPUTED MEAN ANOMALIES.       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(DABS(FDIFF).LT.RPLIM) GOTO 1120
C     IF(ABS(FDIFF).LT.RPLIM) GOTO 1120
      WRITE(IUN6,7003) GLAT,ELON,ELEV,FREAIR,FMOD,FDIFF,FMODR
c      write(*,7003) GLAT,ELON,ELEV,FREAIR,FMOD,FDIFF,FMODR
 1120 WRITE(IUNOA) GLATN,GLATS,ELONW,ELONE,ELEV,FREAIR,FMOD
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     UPDATE HISTOGRAMS.                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(ITER.GT.0) GOTO 1123
      CALL GEOHIS(IUN6,IUN7,IENHIS,CHIS,XMIN,XMAX,RDX,1,FREAIR,NEWP)
 1123 CONTINUE
      CALL GEOHIS(IUN6,IUN7,IENHIS,CHIS,XMIN,XMAX,RDX,2,FMOD,NEWP)
      CALL GEOHIS(IUN6,IUN7,IENHIS,CHIS,XMIN,XMAX,RDX,3,FDIFF,NEWP)
 1125 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     UPDATE NEW SPHERICAL HARMONIC MODEL COEFFICIENTS.                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 1130 MP1=1,MAXP1
      M=MP1-1
      DO 1130 LP1=MP1,MAXP1
      LP1MP1=(MP1-1)*(MAXP12-MP1)/2+LP1
      RCLMN(LP1MP1)=RCLMN(LP1MP1)+FACTL(LP1)*FDIFF*PLMSTO(LP1MP1)*
     1 CMI(MP1)
      RSLMN(LP1MP1)=RSLMN(LP1MP1)+FACTL(LP1)*FDIFF*PLMSTO(LP1MP1)*
     1 SMI(MP1)
 1130 CONTINUE
      GOTO 1010
 2000 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ALL ANOMALIES PROCESSED FOR THIS ITERATION.                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RMSDIF=DSQRT(RMSDIF/FLOAT(NUMBT))
C     RMSDIF=SQRT(RMSDIF/FLOAT(NUMBT))
      WRITE(IUN6,7002) ITER,NUMBT,RMSDIF
      write(*,7002) ITER,NUMBT,RMSDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PRINT HISTOGRAMS FOR THIS ITERATION.                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IENHIS=3
      IF(ITER.GT.0) GOTO 2005
      CALL GEOHIS(IUN6,IUN7,IENHIS,CHIS,XMIN,XMAX,RDX,1,FREAIR,NEWP)
 2005 CONTINUE
      CALL GEOHIS(IUN6,IUN7,IENHIS,CHIS,XMIN,XMAX,RDX,2,FMOD,NEWP)
      CALL GEOHIS(IUN6,IUN7,IENHIS,CHIS,XMIN,XMAX,RDX,3,FDIFF,NEWP)
      IENHIS=2
 2010 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     STORE COEFFICIENTS OF THIS ITERATION IN RCLMO, RSLMO.            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      KS=0
      DO 2020 MP1=1,MAXP1
      DO 2020 LP1=MP1,MAXP1
      KS=KS+1
      RCLMO(KS)=RCLMN(KS)
 2020 RSLMO(KS)=RSLMN(KS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     START NEXT ITERATION OR STOP ITERATIONS.                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ITER=ITER+1
      III=IUNOA
      IUNOA=IUNIA
      IUNIA=III
      IF(ITER.LE.ITMAX) GOTO 1000
 3000 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     END OF ITERATION.                                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 3200 LP1=1,MAXP1
      ADVO(LP1)=0.
      ADVN(LP1)=0.
 3200 ADVD(LP1)=0.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PRINT COEFFICIENTS OF NEW MODEL AND STORE THEM ON TAPE IUN20.    C
C     STANDARD DEVIATIONS OF NEW MODEL ARE SET TO STANDARD DEVIATIONS  C
C     OF START MODEL (IF THERE HAS BE USED ONE).                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REWIND IUN20
      WRITE(*,7013) LMAX,LMAX,numbt
      DO 3260 MP1=1,MAXP1
      M=MP1-1
      IF(ISTAMO.EQ.1) GOTO 3240
      DO 3230 LP1=MP1,MAXP1
      RCOLD(LP1)=0.
      RSOLD(LP1)=0.
      SCOLD(LP1)=0.
      SSOLD(LP1)=0.
 3230 CONTINUE
 3240 CONTINUE
      IF(ISTAMO.EQ.0) GOTO 3250
      READ(IUN33) (RCOLD(LP1),LP1=MP1,MAXP1)
      READ(IUN33) (RSOLD(LP1),LP1=MP1,MAXP1)
      READ(IUN33) (SCOLD(LP1),LP1=MP1,MAXP1)
      READ(IUN33) (SSOLD(LP1),LP1=MP1,MAXP1)
 3250 CONTINUE
      DO 3260 LP1=MP1,MAXP1
      L=LP1-1
      LP1MP1=(MP1-1)*(MAXP12-MP1)/2+LP1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     COMPUTE ANOMALY DEGREE VARIANCES OF NEW MODEL AND OF DIFFERENCES C
C     BETWEEN NEW AND START MODEL.                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ADVO(LP1)=ADVO(LP1)+RCOLD(LP1)**2+RSOLD(LP1)**2
      ADVN(LP1)=ADVN(LP1)+RCLMN(LP1MP1)**2+RSLMN(LP1MP1)**2
      DCLM=RCLMN(LP1MP1)-RCOLD(LP1)
      DSLM=RSLMN(LP1MP1)-RSOLD(LP1)
      ADVD(LP1)=ADVD(LP1)+DCLM*DCLM+DSLM*DSLM
      RCOUT=RCLMN(LP1MP1)
      RSOUT=RSLMN(LP1MP1)
      SCOUT=SCOLD(LP1)
      SSOUT=SSOLD(LP1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ADD NORMAL POTENTIAL.                                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(M.EQ.0.AND.LP1.LE.LPMAX) RCOUT=RCOUT+CN(LP1)
c
c  modified output format of potential coefficients
c
      write(iun20,7001) l,m,rcout,rsout
c
 3260 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PRINT ANOMALY DEGREE VARIANCES.                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(IUN6,7022)
      write(*,7022)
      GMIR2=(GM/RMEAN)**2
      AIR2=(A/RMEAN)**2
      AIRL2=1.
      CF=1.E10/(RMEAN*RMEAN)
      FL=0.
      DO 3300 LP1=1,MAXP1
      L=LP1-1
      ADVO(LP1)=ADVO(LP1)*(FL-1.)*(FL-1.)*GMIR2*AIRL2*CF
      ADVN(LP1)=ADVN(LP1)*(FL-1.)*(FL-1.)*GMIR2*AIRL2*CF
      ADVD(LP1)=ADVD(LP1)*(FL-1.)*(FL-1.)*GMIR2*AIRL2*CF
      DADV=ADVN(LP1)-ADVO(LP1)
      WRITE(IUN6,7020) L,ADVO(LP1),ADVN(LP1),DADV,ADVD(LP1)
      write(*,7020) L,ADVO(LP1),ADVN(LP1),DADV,ADVD(LP1)
      AIRL2=AIRL2*AIR2
 3300 FL=FL+1.
 4000 WRITE(IUN6,7014)
      write(*,7014)
      RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FORMAT STATEMENTS.                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 7000 FORMAT(/ ' ROUTINE GEOPMI ,VERSION 870918 FORTRAN 77.',/,
     1 ' ==========================================',//,
     2 ' OBSERVATIONS WILL BE READ FROM UNIT',I5,/,
     3 ' COEFFICIENTS WILL BE CALCULATED UP TO DEGREE AND ORDER',I5)
 7001 FORMAT(2I5,2d20.12,2d15.8)
 7002 FORMAT(/' ITERATION NUMBER        :',I10/
     1' NUMBER OF OBSERVATIONS  :',I10/
     2' RMS DIFFERENCE          :',F10.3,' MGAL'/)
 7003 FORMAT(' ',2F10.5,5F10.3)
 7011 FORMAT(/,' *** ERROR IN ROUTINE GEOPMI',/,
     1 ' *** MAXIMUM POSSIBLE DEGREE AND ORDER ',I7,' EXCEEDED',/,
     2 ' *** EXECUTION ROUTINE GEOPMI STOPS')
 7013 FORMAT(/' coefficient label - GEOPMI'/
     1 ' max order and degree',2I5,/,
     2 ' SPHERICAL HARMONIC POTENTIAL COEFFICIENTS FROM',/,
     3 ' MEAN FREE AIR ANOMALIES USING INTEGRAL FORMULAS WITH ',/,
     4 ' ROUTINE GEOPMI, VERSION 870914 FORTRAN 77, G. WEBER.',/,
     5 '  ',/,
     6 ' ADDED SPHERICAL HARMONIC COEFFICIENTES OF NORMAL',/,
     7 ' POTENTIAL TAKEN FROM GEODETIC REFERENCE SYSTEM GRS80.',/,
     8 ' NUMBER OF INTRODUCED MEAN FREE AIR ANOMALIES WAS',I7,'.',/,
     9 ' NO ESTIMATION OF COEFFICIENTS STANDARD DEVIATIONS.',/)
 7014 FORMAT(//,' END OF EXECUTION ROUTINE GEOPMI.')
 7015 FORMAT(8A10)
 7016 FORMAT(1X,8A10)
 7017 FORMAT(A10,2I5)
 7018 FORMAT(/' DESCRIPTION OF OLD POTENTIAL MODEL',/,
     1       ' ==================================',//,
     2 ' NAME OF USED SPHERICAL HARMONIC MODEL       ',A10,/,
     3 ' THE MODEL IS COMPLETE TO DEGREE AND ORDER    ',I5,/,
     4 ' COEFFICIENTS ARE GIVEN UP TO DEGREE AND ORDER',I5,//)
 7020 FORMAT(1X,I5,4F12.6)
 7022 FORMAT(/' DIFFERENT ANOMALY DEGREE VARIANCES :',/,
     1       ' ====================================',//,
     2 '     L        ADVO        ADVN   ADVN-ADVO   ADVD(N-O)',/
     3 '         [MGAL**2]   [MGAL**2]   [MGAL**2]   [MGAL**2]',/)
 7026 FORMAT(A10,2F10.5,16X,F10.2,2F8.2)
 7027 FORMAT(/' COMPARISON OF OBSERVED ANOMALIES WITH ANOMALIES FROM',
     2 ' COEFFICIENTS',/,
     3 ' ----------------------------------------------------',
     4 '-------------',//,
     5 ' LIST LIMIT:',F10.3,' MGAL'//
     6 '  LATITUDE LONGITUDE     ELEV.   OBSERV.     MODEL   DIFFER.',
     7 '     GRAD.'/
     8 '  [DEGREE]  [DEGREE]       [M]    [MGAL]    [MGAL]    [MGAL]',
     9 ' [MGAL/KM]'/)
 7029 FORMAT(//,' ITERATION NUMBER IS',I3,/,
     1          ' ----------------------',/)
      END
      SUBROUTINE GEOREF(IUN6,IUN7,CREFSN,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     ROUTINE GEOREF, VERSION 870326 FORTRAN 77.                       C
C                                                                      C
C     SOURCE CODE LIBRARY...SCU_LIB                                    C
C                                                                      C
C     THE ROUTINE GEOREF DEFINES THE PRIMARY CONSTANTS AND COMPUTES    C
C     THE DERIVED CONSTANTS OF THE INTERNATIONAL ELLIPSOID 1930,       C
C                              THE GEODETIC REFERENCE SYSTEM 1967,     C
C                              THE REFERENCE SYSTEM IAG 1975,          C
C                              THE GEODETIC REFERENCE SYSTEM 1980   OR C
C                              THE REFERENCE SYSTEM IAG 1983.          C
C                                                                      C
C     THE PRIMARY PARAMETERS DA, DGM, DJ2 AND DOM FOR THE INTERNATIO-  C
C     NAL ELLIPSOID 1930 HAVE BEEN COMPUTED FROM THE DEFINED PARAME-   C
C     TERS DA, DF, DGE AND DOM WITH PROGRAM INT30 BY H.-G. WENZEL AT   C
C     870724.                                                          C
C                                                                      C
C     REFERENCE... I.A.G. 1967..GEODETIC REFERENCE 1967.               C
C                  PUBLICATION SPECIALE DU BUREAU GEODESIQUE,          C
C                  PARIS 1967.                                         C
C     REFERENCE... MORITZ,H. 1975... REPORT OF SSG NO.5.39 OF I.A.G.   C
C                  PAPER PRESENTED TO 16 TH IUGG GENERAL ASSEMBLY,     C
C                  GRENOBLE 1975.                                      C
C     REFERENCE... MORITZ,H. 1980..GEODETIC REFERENCE SYSTEM 1980.     C
C                  BULLETIN GEODESIQUE, VOL.54 NO.3, P.385-405 ,       C
C                  (THE GEODESIST'S HANDBOOK 1980), PARIS 1980.        C
C     REFERENCE... RAPP,R.H. 1983.. REPORT OF SSG NO. 5.39 OF I.A.G.,  C
C                  FUNDAMENTAL GEODETIC CONSTANTS. PAPER PRESENTED TO  C
C                  18 TH IUGG GENERAL ASSEMBLY, HAMBURG 1983.          C
C                                                                      C
C     INPUT PARAMETER DESCRIPTION...                                   C
C     ==============================                                   C
C                                                                      C
C     IUN6...      FORMATTED LINE PRINTER UNIT.                        C
C     IUN7...      FORMATTED CONSOLE UNIT (DATA TERMINAL SCREEN).      C
C     CREFSN...    DEFINES THE GEODETIC REFERENCE SYSTEM WHICH WILL BE C
C                  USED (CHARACTER*10).                                C
C                  CREFSN='INT 1930  '... THE INTERNATIONAL ELLIPSOID  C
C                                       1930 WILL BE USED.             C
C                  CREFSN='GRS 1967  '... THE GEODETIC REFERENCE       C
C                                       SYSTEM 1967 WILL BE USED.      C
C                  CREFSN='IAG 1975  '... THE REFERENCE SYSTEM IAG     C
C                                       1975 WILL BE USED.             C
C                  CREFSN='GRS 1980  '... THE GEODETIC REFERENCE       C
C                                       SYSTEM 1980 WILL BE USED.      C
C                  CREFSN='IAG 1983  '... THE REFERENCE SYSTEM IAG     C
C                                       1983 WILL BE USED.             C
C                  IF THE INPUT VALUE OF CREFSN DOES NOT AGREE WITH    C
C                  ONE OF THE ABOVE DEFINED STRINGS, THE GEODETIC      C
C                  REFERENCE SYSTEM 1980 WILL BE USED. I.E. IF         C
C                  CREFSN='UNKNOWN ' WILL BE TRANSFERRED TO ROUTINE    C
C                  GEOREF, THE GEODETIC REFERENCE SYSTEM 1980 WILL BE  C
C                  USED.                                               C
C     IPRINT...    LINE PRINTER OUTPUT PARAMETER.                      C
C                  IPRINT=0...NOTHING WILL BE WRITTEN ON UNIT IUN6.    C
C                  IPRINT=1...THE NAME OF THE USED REFERENCE SYSTEM    C
C                             WILL BE WRITTEN ON UNIT IUN6.            C
C                  IPRINT=2...THE PRIMARY AND DERIVED CONSTANTS WILL   C
C                             BE WRITTEN ON UNIT IUN6.                 C
C                                                                      C
C     OUTPUT PARAMETER DESCRIPTION...                                  C
C     ===============================                                  C
C                                                                      C
C     THERE ARE NO OUTPUT PARAMETERS. THE COMPUTED CONSTANTS WILL BE   C
C     TRANSFERRED TO CALLING ROUTINE BY COMMON/REFELL/ AND COMMON      C
C     /REFSYS/.                                                        C
C                                                                      C
C     COMMON BLOCK DESCRIPTION...                                      C
C     ===========================                                      C
C                                                                      C
C     COMMON/REFELL/...  PARAMETERS OF THE REFERENCE ELLIPSOID.        C
C                  THE PARAMETERS OF COMMON/REFELL/ WILL BE DEFINED BY C
C                  CALLING ROUTINE GEOREF.                             C
C                                                                      C
C     GM...        GEOCENTRIC GRAVITATIONAL CONSTANT IN                C
C                  METER**3/SEC**2.                                    C
C     A...         MAJOR SEMI AXIS IN METER.                           C
C     F...         FLATTENING.                                         C
C     OM...        ROTATION SPEED IN RADIANS/SEC.                      C
C     RMEAN...     MEAN EARTH'S RADIUS IN METER.                       C
C     GAMEAN...    MEAN EARTH'S GRAVITY IN METER/SEC**2.               C
C     GE...        EQUATORIAL NORMAL GRAVITY IN METER/SEC**2.          C
C     RK...        CONSTANT FOR SOMIGLIANA NORMAL GRAVITY FORMULA.     C
C     E2...        SQUARE OF FIRST ECCENTRICITY.                       C
C     ES2...       SQUARE OF SECOND ECCENTRICITY.                      C
C     U0...        NORMAL POTENTIAL OF THE LEVEL ELLIPSOID IN          C
C                  METER**2/SEC**2.                                    C
C     CN...        FULLY NORMALIZED ZONAL HARMONIC COEFFICIENTS OF     C
C                  THE ELLIPSOIDAL NORMAL GRAVITY POTENTIAL UP TO      C
C                  DEGREE 10. THE COEFFICIENT C(L, ) IS STORED IN      C
C                  CN(L+1). THE COEFFICIENT C(0,0) IS SET TO 1 AND     C
C                  STORED IN CN(1).                                    C
C                                                                      C
C     COMMON/REFSYS/...NAME OF THE USED REFERENCE ELLIPSOID.           C
C     CREFSY...NAME OF THE USED REFERENCE ELLIPSOID (CHARACTER*10).    C
C                                                                      C
C     USED ROUTINES... NONE.                                           C
C     ================                                                 C
C                                                                      C
C                                                                      C
C     NUMERICAL ACCURACY...                                            C
C     =====================                                            C
C                                                                      C
C     THE RELATIVE ERROR OF ALL VARIABLES IS LESS 10**-14 ON CDC CYBER C
C     76 AND CDC CYBER 990 OF RRZN HANNOVER, RESPECTIVELY.             C
C                                                                      C
C                                                                      C
C     EXECUTION TIME...                                                C
C     =================                                                C
C                                                                      C
C     0.003 SEC CPU TIME WITH IPRINT=0 AND 1, 0.006 SEC CPU TIME WITH  C
C     IPRINT=2 ON CDC CYBER 76 OF RRZN HANNOVER. SAME EXECUTION TIME   C
C     ON CDC CYBER 990 OF RRZN HANNOVER.                               C
C                                                                      C
C     ROUTINE TESTS...                                                 C
C     ================                                                 C
C                                                                      C
C     SUCCESSFULLY COMPLETED 830615 BY F.BOECKMANN                     C
C                                                                      C
C     ROUTINE CREATION...    830408 BY H.-G.WENZEL,                    C
C                            INSTITUT FUER ERDMESSUNG,                 C
C                            UNIVERSITAET HANNOVER,                    C
C                            NIENBURGER STR.6,                         C
C                            D-3000 HANNOVER 1,                        C
C                            FEDERAL REPUBLIC OF GERMANY.              C
C     PROGRAM MODIFICATION...870326 BY H.-G.WENZEL.                    C
C**********************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     IMPLICIT DOUBLE PRECISION (D)
      CHARACTER NRS(5)*10,CREFSN*10,CREFSY*10
      COMMON/REFELL/ GM,A,F,OM,RMEAN,GAMEAN,GE,RK,E2,ES2,U0,CN(11)
      COMMON/REFSYS/ CREFSY
      DATA NRS/'INT 1930  ','GRS 1967  ','IAG 1975  ','GRS 1980  ',
     1'IAG 1983  '/
      KRS=0
      DO 10 I=1,5
      IF(CREFSN.EQ.NRS(I)) KRS=I
   10 CONTINUE
      IF(KRS.EQ.0) GOTO 5000
      GOTO (100,200,300,400,500) KRS
  100 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DEFINE THE PRIMARY PARAMETERS FOR INTERNATIONAL ELLIPSOID 1930.  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DA=6378388.D0
      DGM=398632.904400795D9
      DJ2=1092.03876103097D-6
      DOM=7.292115D-5
      CREFSY=NRS(1)
      GOTO 1000
  200 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DEFINE THE PRIMARY CONSTANTS FOR GRS 1967 REFERENCE ELLIPSOID.   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DA=6378160.D0
      DGM=398603.D9
      DJ2=1082.7D-6
      DOM=7.2921151467D-5
      CREFSY=NRS(2)
      GOTO 1000
  300 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DEFINE THE PRIMARY CONSTANTS FOR IAG 1975 REFERENCE ELLIPSOID.   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DA=6378140.D0
      DGM=398600.5D9
      DJ2=1082.63D-6
      DOM=7.292115D-5
      CREFSY=NRS(3)
      GOTO 1000
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DEFINE THE PRIMARY CONSTANTS FOR GRS 1980 REFERENCE ELLIPSOID.   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  400 CONTINUE
      DA=6378137.D0
      DGM=398600.5D9
      DJ2=1082.63D-6
      DOM=7.292115D-5
      CREFSY=NRS(4)
      GOTO 1000
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DEFINE THE PRIMARY CONSTANTS FOR IAG 1983 REFERENCE ELLIPSOID.   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  500 CONTINUE
      DA=6378136.D0
      DGM=398600.44D9
      DJ2=1082.629D-6
      DOM=7.292115D-5
      CREFSY=NRS(5)
 1000 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     COMPUTE DERIVED CONSTANTS.                                       C
C     ITERATIVE COMPUTATION OF THE SQUARE OF FIRST ECCENTRICITY.       C
C     COMPUTATION OF INITIAL VALUE DE20 FOR THE SQUARE OF FIRST        C
C     ECCENTRICITY ACCORDING TO...                                     C
C                  CHEN,Y.1981...FORMULAE FOR COMPUTING ELLIPSOIDAL    C
C                            PARAMETERS. BULLETIN GEODESIQUE, VOL.55,  C
C                            NO.2, PP.170 - 178, PARIS 1981.           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DQ=DOM*DOM*DA*DA*DA/DGM
      DE20=3.D0*DJ2+DQ-3.D0*DQ*(9.D0*DJ2+3.D0*DQ)/14.D0+(408.D0*DJ2*DQ
     1+149.D0*DQ*DQ-117.D0*DJ2*DJ2)*DQ/392.D0
 1010 CONTINUE
      DES2=DE20/(1.D0-DE20)
      DES=DSQRT(DES2)
      DT=DATAN(DES)
      DM=DQ*DSQRT(1.D0-DE20)
      DQ0=((1.D0+3.D0/DES2)*DT-3.D0/DES)/2.D0
      DE2=3.D0*DJ2+2.D0*DM*DES*DE20/(15.D0*DQ0)
      DERR2=abs(DE2-DE20)
      IF(DERR2.LT.1.0e-15) GOTO 1020
      DE20=DE2
      GOTO 1010
 1020 CONTINUE
      DF=1.D0-DSQRT(1.D0-DE2)
      DES2=DE2/(1.D0-DE2)
      DES=DSQRT(DES2)
      DT=DATAN(DES)
      DE=DSQRT(DE2)
      DU0=DGM*DT/(DA*DE)+DOM*DOM*DA*DA/3.D0
      DM=DQ*(1.D0-DF)
      DQ0=((1.D0+3.D0/DES2)*DT-3.D0/DES)/2.D0
      DQS=3.D0*(1.D0+1.D0/DES2)*(1.D0-DT/DES)-1.D0
      DGE=DGM*1.D0/(DA*DA*(1.D0-DF))*(1.D0-DM-DM*DES/6.D0*DQS/DQ0)
      DGP=DGM*1.D0/(DA*DA)*(1.D0+DM/3.D0*DES*DQS/DQ0)
      DK=(1.D0-DF)*DGP/DGE-1.D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     MEAN RADIUS AND MEAN NORMAL GRAVITY.                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DC=DA/(1.D0-DF)
      DR=DC*(1.D0-2.D0/3.D0*DES2+26.D0/45.D0*DES2*DES2-100.D0/189.D0*
     1 DES2*DES2*DES2+7034.D0/14175.D0*DES2*DES2*DES2*DES2)
      DGAM=DGE*(1.D0+DE2/6.D0+DK/3.D0+59.D0/360.D0*DE2*DE2+5.D0/18.D0*
     1 DE2*DK+2371.D0/15120.D0*DE2*DE2*DE2+259.D0/1080.D0*DE2*DE2*DK+
     2 270229.D0/1814400.D0*DE2*DE2*DE2*DE2+9623.D0/45360.D0*DE2*DE2*
     3 DE2*DK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     COMPUTE FULLY NORMALIZED SPHERICAL HARMONIC COEFFICIENTS.        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CN(1)=1.
      DO 1030 L=1,10
      LP1=L+1
 1030 CN(LP1)=0.
      DEL2=1.D0
      DL2=0.D0
      DO 1040 L2=1,5
      L=2*L2
      DL2=DL2+1.D0
      K=L+1
      DEL2=-DEL2*DE2
      DCN=DEL2*3.D0/(DSQRT(4.D0*DL2+1.D0)*(2.D0*DL2+1.D0)*
     1 (2.D0*DL2+3.D0))
      DCN=DCN*(1.D0-DL2+5.D0*DL2*(1.D0/3.D0-2.D0/45.D0*DM*DES/DQ0))
 1040 CN(K)=DCN
C1040 CN(K)=SNGL(DCN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     TRANSFORMATION FROM DOUBLE PRECISON INTO SINGLE PRECISION.       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      GM=DGM
      A=DA
      F=DF
      OM=DOM
      RMEAN=DR
      RJ2=DJ2
      E2=DE2
      ES2=DES2
      GE=DGE
      GP=DGP
      RK=DK
      GAMEAN=DGAM
      U0=DU0
      IF(IPRINT.EQ.1) WRITE(IUN6,7000) CREFSY
      IF(IPRINT.EQ.1) write(*,7000) CREFSY
      IF(IPRINT.NE.2) RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     OUTPUT ON LINE PRINTER.                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      F1=1./F
      RJ21=RJ2*1.E6
      GM1=GM*1.E-9
      OM1=OM*1.E5
      WRITE(IUN6,7001) CREFSY,A,RJ21,GM1,OM1
      write(*,7001) CREFSY,A,RJ21,GM1,OM1
      WRITE(IUN6,7002) CREFSY,E2,ES2,F1,GE,GP,RK,RMEAN,GAMEAN,U0
      write(*,7002) CREFSY,E2,ES2,F1,GE,GP,RK,RMEAN,GAMEAN,U0
      WRITE(IUN6,7003) CREFSY
      write(*,7003) CREFSY
      DO 1050 L2=1,5
      L=2*L2
      LP1=L+1
      WRITE(IUN6,7004) L,CN(LP1)
 1050 write(*,7004) L,CN(LP1)
      WRITE(IUN6,7006)
      write(*,7006)
      RETURN
 5000 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PARAMETER CREFSN NOT ALLOWED. GEODETIC REFERENCE SYSTEM GRS 1980 C
C     WILL BE USED AND THE EXECUTION WILL BE CONTINUED.                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      KRS=4
      CREFSY=NRS(4)
      WRITE(IUN6,7005) CREFSY
      GOTO 400
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FORMAT STATEMENTS.                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 7000 FORMAT(/' ROUTINE GEOREF, VERSION 870326, FTN 77.'/
     1 ' GEODETIC REFERENCE SYSTEM USED IS',2X,A8/)
 7001 FORMAT('1ROUTINE GEOREF, VERSION 870326, FTN 77.'/
     1 ' PRIMARY CONSTANTS OF GEODETIC REFERENCE SYSTEM',2X,A8//
     2 ' MAJOR SEMI AXIS',15X,F15.1,' METER'/
     3 ' J2',28X,F15.3,' *10**-6'/
     4 ' GEOCENTRIC GRAVITATION', 8X,F15.2,' *10**9 METER**3/SEC**2'/
     5 ' ROTATION SPEED',16X,F15.10,' *10**-5 RADIAN/SEC'/)
 7002 FORMAT(//' DERIVED CONSTANTS OF GEODETIC REFERENCE SYSTEM',2X,
     1 A8//' SQUARE OF FIRST ECCENTRICITY', 1X,F16.14/
     2 ' SQUARE OF SECOND ECCENTRICITY', F16.14/
     3 ' FLATTENING',18X,'1/',F15.10/
     4 ' EQUATORIAL NORMAL GRAVITY', 5X,F15.8,' METER/SEC**2'/
     5 ' POLAR NORMAL GRAVITY',10X,F15.8,' METER/SEC**2'/
     6 ' CONSTANT OF SOMIGLIANA FORMULA',F15.12/
     7 ' MEAN EARTH RADIUS',13X,F15.1,' METER'/
     8 ' MEAN EARTH NORMAL GRAVITY', 5X, F15.8,' METER/SEC**2'/
     9 ' NORMAL POTENTIAL',14X,F15.3,' METER**2/SEC**2'/)
 7003 FORMAT(/' FULLY NORMALIZED ZONAL SPHERICAL HARMONIC',
     1 ' COEFFICIENTS'/
     2 ' OF THE ELLIPSOIDAL NORMAL GRAVITY POTENTIAL FOR THE'/
     3 ' GEODETIC REFERENCE SYSTEM',2X,A8/
     4 '    DEGREE                   C(L,0)'/)
 7004 FORMAT(I10,E25.15)
 7005 FORMAT(//' *****ERROR IN ROUTINE GEOREF, VERSION 870326 FTN 77.'/
     1 ' *****THE PARAMETER CREFSN=',A10,' USED IN THE CALL OF ROUTINE',
     2 ' GEOREF IS NOT ALLOWED.'/
     4 ' *****THE GEODETIC REFERENCE SYSTEM 1980 WILL BE USED AND THE',
     5 ' EXECUTION WILL BE CONTINUED.'//)
 7006 FORMAT(// ' *****EXECUTION OF ROUTINE GEOREF FINISHED.')
       END


