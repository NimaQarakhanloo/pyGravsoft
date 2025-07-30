      PROGRAM ORBIT
C $Id: orbit10.for 181 2008-08-14 23:19:51Z tjansson $
C------------------------------------------------------------------------
C
C THIS PROGRAM COMPUTES A  X,Y,Z,DX/DT,DY/DT,DZ/DT  STATE VECTOR FILE
C GIVEN A KEPLERIAN ORBIT THAT IS PERTURBED ONLY THROUGH A J2 TERM.
C HERE THE J2 SECULAR SOLUTIONS ARE TAKEN INTO ACCOUNT, THEY ARE
C DESCRIBED IN (SCHRAMA,1989), THE ROLE OF RADIAL ORBIT ERRORS IN
C PROCESSING OF SATELLITE ALTIMETER DATA, P.50 EQNS. (4.20).
C
C RECEIVED FROM EOSCHRAMA, JUNE 28, 1989, MODIFIED BY C.C.TSCHERNING. 
C LAST MODIFICATION 2000.02.11
C------------------------------------------------------------------------
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /COMPI/ PI,PI2,DEGRAD,DD,GM 
      LOGICAL LCOMP,LOUTFI
      CHARACTER*72 INAME,ONAME
      REAL*8 S(6),E(6),CH(3),ECOPY(6),AE,C20
      REAL*8 A0,E0,I0,O0,W0,M0
      REAL*8 DODT,DWDT,DMDT
      REAL*8 RTO,RTW,RTM
      REAL*8 T,T0,TB,TE,DT
      REAL*8 M,DUMMY,F,W,O,IIT 
      INTEGER J,K,NOSTPS
      COMMON /ITRANC/SHIFT(10),AX2,E22 
      REAL*8 DXVEC(3)
      EQUIVALENCE (DXVEC(1), DX), (DXVEC(2), DY), (DXVEC(3), DZ)
      REAL*8 DXVECS(3),DDS(6),SI(6),DSI(6),DDSI(6),SMI(6),
     * SSMI(6),SOLD(6),SIOLD(6),DG(6)
      EQUIVALENCE (DXVECS(1), DXS), (DXVECS(2), DYS), (DXVECS(3), DZS)
      REAL*8 ROTMX(3, 3)
C
      WRITE(*,*)' ORBIT * VERS. 9, LAST UPDT. 2000.02.11' 
C
C WHAT YOU ALWAYS NEED
C
      PI  = 4.0D0*DATAN(1D0)
      PI2 = PI*2.0D0 
      DEGRAD=PI/180.0D0
      RADDEG=1.0D0/DEGRAD
      GM  = 3986005D8
      AE  = 6378137D0
      C20 = -0.00108263D0
      AX2=AE
      E22=0.00669438002D0 
      D0=0.0D0
      D1=1.0D0
C
C THE ORBITAL PARAMETERS AT T=TB
C
C A0 --> SEMI MAJOR AXIS IN METERS
C E0 --> ECCENTRICITY
C I0 --> INCLINATION IN DEGREES
C O0 --> RIGHT ASCENSION OF THE ASCENDING NODE IN DEGREES
C W0 --> ARGUMENT OF PERIGEE IN DEGREES
C T0 --> OFFSET TIME AT W0
C
C TIMES (IN SECONDS)
C
C TB --> START TIME OF EPHEMERIS ( M(T)=N*(T-T0), N=DM/DT )
C TE --> END TIME OF EPHEMERIS
C DT --> STEPSIZE
C
C VALUES USED ORIGINALLY BY SHRAMA AS INPUT: 
C     A0 = AE+790D3
C     E0 = 0.001D0
C     I0 = 108D0
      O0 = 45D0
      W0 = 60D0
C     TE = 20000D0
C     DT = 600D0
      T0 = 0D0
      TB = 0D0
      write(*,190)
  190 format(
     *' THE PROGRAM WILL WORK IN FIVE MODES: 0,1,2,3,4.'/
     *' IN MODE 0 THE PROGRAM WILL GIVE THE OUTPUT AS ORIGINALLY'/
     *' DESIGNED C BY SCHRAMA; '/
     *' FOR MODE 1 THE OUTPUT IS TIME, LATITUDE, LONGITUDE (DEC. DEG.)'/
     *' HEIGHT (KM), VELOCITY VECTOR IN INERTIAL AND EARTH-FIXED'/
     *' FRAME, AND INCLINATION/TILT ANGLES WITH RESPECT TO'/
     *' ELLIPSOIDAL NORMAL.'/ 
     *' IN MODE 2 THE OUTPUT IS ONLY THE TIME AND LATITUDE, LONGITUDE'/ 
     *' AND HEIGHT IN M.'/ 
     *' IN MODE 3 THE OUTPUT IS THE TIME AND LATITUDE, LONGITUDE,'/
     *' HEIGHT IN M, AZIMUTH AND TILT OF VELOCITY VECTOR COMPUTED WITH'/
     *' RESPECT TO THE RADIUS VECTOR THROUGH THE SATELLITE.',/
     *' IN MODE 4 IT IS THE KEPLER ELEMENTS '/) 
C 
      WRITE(*,*)
     *' INPUT MODE (0,1,2,3,4) '
      READ(*,*)MODE  
      IF (MODE.LT.0.OR.MODE.GT.4) WRITE(*,*)' WRONG MODE' 
      IF (MODE.NE.4) THEN
      WRITE(*,*)' INPUT HEIGHT (M) or n,EXCENTRICITY,INCL (DEG.)'
      WRITE(*,*)'  CAP. OMEGA, SMALL OMEGA, MEAN ANOMALY  (DEG),' 
      READ(*,*)H,E0,I0,O0,W0,M0  
      IF (H.LT.100.0) THEN
C H IS Mean Motion, n in Rev. per. day. 
      RH=PI2*H/86400.0d0 
      A=(GM/RH**2)**(1.0d0/3.0d0)
      H=A-AE 
      END IF 
      ELSE
      WRITE(*,*)' INPUT STATE-VECTOR '
      READ(*,*)S
      WRITE(*,*)' COMPARISON WITH FILE VALUES ? (T/T) '
      READ(*,*)LCOMP
      IF (LCOMP) THEN
      WRITE(*,*)' INPUT NAME OF FILE '
      READ(*,'(A)')INAME
      OPEN(8,FILE=INAME)
      WRITE(*,*)' FILE FOR COMPARISON: ',INAME
      WRITE(*,*)' OUTPUT IS FILE VALUES AND DIFFERENCES '
      DO 2075,IIN=1,6
      SOLD(IIN)=D0
      SIOLD(IIN)=D0
      SMI(IIN)=D0
 2075 SSMI(IIN)=D0
      END IF
      WRITE(*,*)' OUTPUT TO FILE ? (T/F) '
      READ(*,*)LOUTFI
      IF (LOUTFI) THEN
      WRITE(*,*)' INPUT OUTPUT FILE NAME '
      READ(*,'(A)')ONAME
      OPEN(9,FILE=ONAME)
      END IF
      CALL VECELE(S,E,CH)
      A0=E(1)
      E0=E(2)
      I0=E(3)*RADDEG
      O0=E(4)*RADDEG
      W0=E(5)*RADDEG
      F0=E(6)*RADDEG
      CALL TRANSMF( E0,M0,DUMMY,E(6),-1 )
      M0=M0*RADDEG
      H=A0-AE 
C
      CALL ELEVEC(S,E)
C
      write(*,*)' INITIAL STATE VECTOR (M) '
      write(*,188)S
      END IF
      WRITE(*,*)'  START DAY EPHEMERIS, START DAY ORBIT, END DAY ' 
      WRITE(*,*)'  YY MM DD HH TT SS.SS' 
      READ(*,*)JJ,MO,IDD,IHH,ITT,SS 
      IF (MO.NE.0) THEN 
      CALL GETJUL(JJ,MO,IDD,IHH,ITT,SS,SDEP) 
      ELSE
      EP=DJULIA(JJ,0,0) 
      SDEP=EP+SS
      END IF 
      write(*,2076)JJ,MO,IDD,IHH,ITT,SS 
      READ(*,*)JJ,MO,IDD,IHH,ITT,SS 
      write(*,2076)JJ,MO,IDD,IHH,ITT,SS 
      CALL GETJUL(JJ,MO,IDD,IHH,ITT,SS,SDOR) 
      READ(*,*)JJ,MO,IDD,IHH,ITT,SS 
      write(*,2076)JJ,MO,IDD,IHH,ITT,SS 
 2076 FORMAT(I5,I3,I3,I3,I3,F8.4)
      CALL GETJUL(JJ,MO,IDD,IHH,ITT,SS,EDOR)
      write(*,175)SDEP,SDOR,EDOR
  175 FORMAT(' EPOCH=',D15.7,', START=',D15.7,', END=',D15.7) 
      WRITE(*,*)'  STEPSIZE (S), MIN&MAX LAT, LONG (DEG.)' 
      READ(*,*)DT,RLATMIN, RLATMAX,
     *RLONMIN, RLONMAX 
      A0=AE+H 
      NOBS=0 
      CALL GETJUL(1973,3,21,0,0,0.0d0,D73) 
      CALL GETJUL(1985,1,1,0,0,0.0d0,D85) 
      D85=(D85-D73)*24*3600  
      TB=(SDOR-D73)*24*3600
      T0=(SDEP-D73)*24*3600 
      TE=(EDOR-D73)*24*3600 
      IF (TB.LT.T0.OR.TB.GT.TE) WRITE(*,*)' TIMES WRONG' 
C
C PRINT SOME STATISTICS
C
      WRITE(6,1000) GM,C20,AE
 1000 FORMAT(' GM              ',G20.10,' (M**3/S**2)' /
     .       ' C20             ',F20.10                /
     .       ' AE              ',8X,F12.1,'         (M) '        /)

      WRITE(6,1010) A0,E0,I0,O0,W0,F0,DUMMY*RADDEG,T0
 1010 FORMAT(' SEMI MAJOR AXIS ',F13.3,7X ,' (M) '        /
     .       ' ECCENTRICITY    ',F13.10,7X,                /
     .       ' INCLINATION     ',F15.10,5X,' (DEG) '      /
     .       ' R.A. NODE       ',F15.10,5X,' (DEG) '      /
     .       ' ARG. OF PERIGEE ',F15.10,5x,' (DEG) '      /
     .       ' TRUE INITIAL ANO',F15.10,5x,' (DEG) '      /
     .       ' EXCC INITIAL ANO',F15.10,5x,' (DEG) '      /
     .       ' OFFSET TIME     ',F15.2,5x,' (S)   '      /)

C
C COMPUTE THE SECULAR DRIFT RATES DODT,DWDT,DMDT
C
      I0=I0/180D0*PI
      O0=O0/180D0*PI
      W0=W0/180D0*PI
      M0=M0/180D0*PI 

      CALL SECULA( A0,E0,I0,C20,AE,DODT,DWDT,DMDT )

      RTO=2D0*PI/DODT/86400D0
      RTW=2D0*PI/DWDT/86400D0
      RTM=2D0*PI/(DMDT+DWDT)
      RTDAY=RTM/(24*3600.0) 

      WRITE(6,1020) RTO,RTW,RTM,RTDAY  
 1020 FORMAT(' REVTIM R.A.NODE ',F16.6,' (24H DAYS) ' /
     .       ' REVTIM PERIGEE  ',F16.6,' (24H DAYS) ' /
     .       ' REVOLUTION TIME ',F16.6,' (SEC)      ' /
     .       ' REVOLUTION TIME ',F16.8,' (DAYS)'/)

      TEDAY=(TE-TB)/(24*3600.0) 
      WRITE(6,1030) TB,TE,TEDAY,DT
 1030 FORMAT(' START TIME      ',F16.6,' (SEC)      ' /
     .       ' STOP  TIME      ',F16.6,' (SEC)      ' /
     .       ' STEP  TIME      ',F16.6,' (DAYS)      '/ 
     .       ' STEP SIZE       ',F16.6,' (SEC)      ' /)

C
C AND HERE FOLLOWS THE OUTPUT
C
      E(1)=A0
      E(2)=E0
      E(3)=I0

      IF (MODE.EQ.0) THEN 
      WRITE(6,*) 'OUTPUT AS T,     X,     Y,     Z     (INERTIAL)'
      WRITE(6,*) '             DX/DT, DY/DT, DZ/DT     (INERTIAL)'
      WRITE(6,*) '                 O,     W,     F '
      WRITE(6,*)
      END IF
      IF (MODE.EQ.1) THEN
      WRITE(6,*)
     *'OUTPUT AS T, LAT, LON, H, DX/DT,DY/DT,DZ/DT '
      WRITE(6,*)
     *' VELOCITY IN EARTH-FIX FRAME, & R, LON, PHI, AZ, BETA, PHI' 
      END IF 
      IF (MODE.EQ.3) WRITE(6,*)
     *'OUTPUT AS T, LAT, LON, H, AZ, BETA IN SPHERICAL FRAME ' 
C 
      NOSTPS=NINT((TE-TB)/DT)
      write(*,*)' NUMBER OF STEPS ', nostps 
      IF (NOSTPS.LT.200) THEN
      NSTOUT=0
      ELSE
      WRITE(*,*)' MORE THAN 200 STEPS. INPUT OUTPUT STEPS '
      READ(*,*)NSTOUT
      END IF
C
      DO J=0,NOSTPS
        T=T0+DBLE(J)*DT
        M=DMDT*(T-T0)+M0 
        CALL TRANSMF( E0,M,DUMMY,F,1 )
        W=DWDT*(T-T0)+W0
        O=DODT*(T-T0)+O0
        if (.not.lcomp)
     * write(*,191)DT,DMDT*DT*RADDEG,DWDT*DT*RADDEG,DODT*DT*RADDEG
  191   format(' DT,DM,DW,DO ',4f12.6)
        E(4)=O
        E(5)=W
        E(6)=F
C CONVERSION FACTOR RAD TO DEG:
        DD=180D0/PI 
	DEGRAD=PI/180D0 
        CALL ELEVEC(S,E)
        CALL ANGLE(T,TT)
C TRANSFORM LOCATION TO SOLID-EARTH FRAME:
C NOTE RLAME, PHIE ARE GEOCENTRIC!
        CALL XITOXE(T,S(1),S(2),S(3),XE,YE,ZE,RADE,RLAME,PHIE)
C TRANSFORM VELOCITY TO SOLID-EARTH FRAME:
        CALL XITOXE(T,S(4),S(5),S(6),DX,DY,DZ,RA1,RL1,PH1) 
C DATUM TRANSFORMATION. HERE USED TO COMPUTE GEODETIC COORDS
C (RLA, RLO SATELLITE SOLID EARTH GEODETIC LAT, LON IN RADIANS):
C INTO SLA, CLA, RLA, SLO, CLO, RLO, H.
        CALL TRANS(SLA,CLA,RLA,SLO,CLO,RLO,H,XE,YE,ZE,0) 
	RLAD=RLA*DD
	RLOD=RLO*DD 
        IF (RLAD.GT.RLATMIN.AND.RLAD.LT.RLATMAX.AND.
     *      RLOD.GT.RLONMIN.AND.RLOD.LT.RLONMAX.AND.MODE.GE.1) THEN 
        NOBS=NOBS+1 
        IIT=T-D85 
C DIFF GEOCENTRIC - GEODETIC CCORDS. SOME 11 MIN OF ARC MAX. AT 45 DEG.
        PHIED=PHIE*DD 
	RLAMED=RLAME*DD 
        DPHI=PHIED-RLAD 
	IF (MODE.EQ.3) THEN
	CLO=COS(RLAME)
	SLO=SIN(RLAME)
	CLA=COS(PHIE)
	SLA=SIN(PHIE) 
	END IF 
C       RADE=RADE-6371000.0 
C       WRITE(6,1049) T,PHIED,RLAMED,RADE,(E(K)*180D0/PI,K=4,6) 
C ROTATION MATRIX 
        ROTMX(1,1) = -SLO
	ROTMX(1,2) =  CLO
        ROTMX(1,3) =  0.0D0
	ROTMX(2,1) = -SLA * CLO
	ROTMX(2,2) = -SLA * SLO
        ROTMX(2,3) =  CLA
	ROTMX(3,1) =  CLA * CLO
	ROTMX(3,2) =  CLA * SLO
	ROTMX(3,3) =  SLA
C TRANSFORM NOW DXVEC TO DXVECS, THE SATELLITE COORD. FRAME.
C IN THIS FRAME X IS EAST, Y NORTH, Z UP (ELLIPSOIDAL NORMAL) 
        DO 2010 I = 1,3        
	DXVECS(I) = 0.0
        DO 2010 JJ= 1,3
        DXVECS(I) = DXVECS(I) + ROTMX(I,JJ) * DXVEC(JJ)
 2010   CONTINUE
C USE SATELLITE FRAME DXS, DYS, DZS FOR AZ, BETA COMPUTATION:
C AZIMUTH ZERO AT NORTH, POSITIVE EASTWARD:
        AZ = ATAN2(DXS, DYS) * DD
        IF (AZ .LT. 0.0D0) AZ = AZ + 360.0D0
C BETA POSITIVE FOR CLIMBING SATELLITE:
C NOTE (N.B. N.B. N.B.) THAT BETA IS W.R.T. THE ELLIPSOIDAL
C EARTH SURFACE DUE TO USE OF A ROT. MATRIX FOR GEODETIC LAT, LON!!!!
        BETA = ATAN2(DZS, SQRT(DXS * DXS + DYS * DYS)) * DD 
C
        IF (MODE.EQ.1) 
     *  WRITE(6,1049) IIT,RLAD,RLOD,H/1000,(S(K),K=4,6),
     *  DX,DY,DZ, RA1,RL1,PH1,AZ,BETA,DPHI  
C    *  (E(K)*DD,K=4,6)  
 1049   FORMAT(X,I10,2F11.4,F8.2,3F10.1,/,10X,3F10.1,F8.1,2F10.4,
     *  /,20X,3F10.2) 
        IF (MODE.EQ.2) WRITE(6,1043) IIT,RLAD,RLOD,H
 1043   FORMAT(X,F11.1,2F12.6,F9.1,2F12.6)
        IF (MODE.EQ.3) WRITE(6,1043) IIT,RLAD,RLOD,H,AZ,BETA 
        END IF 
        IF (MODE.EQ.0) 
     *  WRITE(6,1040) T,S,(E(K)*DD,K=4,6),XE,YE,ZE,RADE,RLAMED,PHIED 
 1040   FORMAT(X,F13.1,3F20.4/13X,3F20.4/13X,3F20.4/
     *13X,3F20.4/13X,F20.4,2F20.8/)
        IF (MODE.EQ.4) THEN
        IF (LCOMP) THEN
        READ(8,*)nn,iyy,imm,idd,iti,imi,seci,si
        DO 2041, nn=1,6
        SI(nn)=SI(nn)*1.0D3
        DG(NN)=SI(NN)-SIOLD(NN)
        dsi(nn)=s(nn)-si(nn)
        DDS(NN)=SI(NN)-SIOLD(NN)
        ddsi(nn)=dsi(nn)-SOLD(nn)
        SOLD(NN)=DSI(NN)
        SIOLD(NN)=SI(NN)
        SMI(nn)=SMI(nn)+DSI(nn)
 2041   SSMI(nn)=SSMI(nn)+DSI(nn)**2
c
        SR=SQRT(SI(1)**2+SI(2)**2+SI(3)**2)
        SCN=SQRT(Ch(1)**2+CH(2)**2+CH(3)**2)
        DR=D0
        DH=D0
        GG=D0
        GT=D0
        DO 2042,NN=1,3
        GG=GG+DG(NN+3)*SI(NN)/SR
        GT=GT+DG(NN+3)*CH(NN)/SCN
        DH=DH+DDSI(NN+3)*CH(NN)/SCN
 2042   DR=DR+DDSI(NN+3)*SI(NN)/SR 
        IF (J.EQ.0)DR=0
c
c       IF (MOD(J,NSTOUT).EQ.0.OR.NSTOUT.EQ.0) THEN
        write(*,188)si,DH
        Write(*,188)dsi,DR
c output of gravity in mgal and positive inward.
        WRITE(6,1043) IIT,RLAD,RLOD,H,GG      
        INO=IIT
        IF (LOUTFI.AND.J.GT.1)
     *  WRITE(9,1044) INO,RLAD,RLOD,H,-GG*1.0D5,GT*1.0D5      
 1044   FORMAT(X,I11,2F12.6,F9.1,2F9.1)
c       write(*,199)GG
c 199   format(F8.5) 
c       END IF
        else
        WRITE(*,188)S
  188   FORMAT(3F14.3,3F9.3,F7.4)
        WRITE(*,189)E(1),E(2),E(3)*RADDEG,E(4)*RADDEG
     *  ,E(5)*RADDEG,E(6)*RADDEG
  189 FORMAT(3D16.8)
        END IF
        END IF
      END DO
C
      IF (LOUTFI) CLOSE(9)
      IF (LCOMP) THEN
      CLOSE(8)
      write(*,*)' MEAN AND STDV OF DIFFERENCES (M) '
      DO 2079, J=1,6
      SSMI(J)=SQRT((SSMI(J)-SMI(J)**2/NOSTPS)/(NOSTPS-1))
 2079 SMI(J)=SMI(J)/NOSTPS
C
      WRITE(*,2078)SMI,SSMI
 2078 FORMAT(6F10.2)
      END IF
C
      IF (MODE.EQ.1) WRITE(*,*) ' POINTS',NOBS 
C
      STOP
      END
      SUBROUTINE TRANS(SINLAP,COSLAP,RLATP,SINLOP,COSLOP,RLONGP,HP,
     *X,Y,Z,MODE) 
C ORIGINAL VERSION PROGRAMMED IN 1974 BY C.C.TSCHERNING, GEODAETISK
C INSTITUT. LATEST UPDATE 7 OCT 1987.
C
C THE SUBROUTINE TRANSFORMS THE COORDINATES FROM ONE DATUM TO ANOTHER
C USING THE 7-PARAMETER DATUM-SHIFT GIVEN BY DX,DY,DZ,DL,EPS1,EPS2,
C EPS3. IF MODE=0, NO SHIFT IS DONE. 
C
C
C IF DOUBLE PRECISION IS NEEDED ACTIVATE:
      IMPLICIT INTEGER(I,J,K,N,M),LOGICAL(L),REAL *8(A-H,O-Z)
C AND USE DSIN, DCOS, DATAN2, DABS BELOW.
      COMMON /ITRANC/SINLA0,COSLA0,RLONG0,
     *DX,DY,DZ,EPS3,EPS2,EPS1,S1,AX2,E22
      D0=0.0D0 
      D1=1.0D0 
      X0=X
      Y0=Y
      Z0=Z 
      IF (MODE.NE.0) THEN 
      X1 = X
      Y1 = Y
      X0 = DX+S1*(X+EPS1*Y-EPS2*Z)
      Y0 = DY+S1*(Y-EPS1*X1+EPS3*Z)
      Z0 = DZ+S1*(Z+EPS2*X1-EPS3*Y1)
      END IF 
      XY20= X0*X0+Y0*Y0
      XY0 =  SQRT(XY20)
      DIST20 = XY20+Z0*Z0
      DISTO0 =  SQRT(DIST20)
      RLONG =  ATAN2(Y0,X0)
C
C  COMPUTATION OF THE NEW GEODETIC LATITUDE, CF REF(C) PAGE 183.
      S = AX2
      DH = D0
      RLAT1 = D0 
      COSLA=1.0D0 
   70 RLAT = RLAT1
C
      RLAT1 =  ATAN2(Z0,XY0-E22*S*COSLA)
      COSLA =  COS(RLAT1)
      S = AX2/ SQRT(D1-E22*(D1-COSLA**2))
      DH = XY0/COSLA-S
      IF ( ABS(RLAT1-RLAT).GT.1.0D-10) GO TO 70
C
      RLONGP = RLONG
      RLATP = RLAT1
      SINLOP= SIN(RLONG)
      COSLOP= COS(RLONG)
      COSLAP=COSLA
      SINLAP= SIN(RLATP)
      HP=DH 
      IF (MODE.NE.0) THEN
      X=X0
      Y=Y0
      Z=Z0
      XY=XY0
      XY2=XY20
      DISTO=DISTO0
      DIST2=DIST20
      END IF 
C
      RETURN
      END 
C------------------------------------------------------------------------
C
C SUBROUTINE XITOXE TRANSFORM THE INERTIAL "XI" COORDINATE FRAME TO
C THE EARTH FIXED "XE" COORDINATE FRAME. ADDITIONALLY THE SPHERICAL
C EARTH FIXED COORDINATES RADE "RADIUS" LAME "LONGITUDE" AND PHIE
C "LATITUDE" ARE COMPUTED (THE LATTER TWO IN DEGREES).
C
C------------------------------------------------------------------------
      SUBROUTINE XITOXE( TIME, XI,YI,ZI, XE,YE,ZE, RADE,LAME,PHIE )
      IMPLICIT REAL*8 (A-Z)
      COMMON /COMPI/ PI,PI2,DEGRAD,RADDEG,GM 
      CALL ANGLE(TIME,TT )
C     TTS=SIDTIME(TIME) 
      CT=DCOS(TT)
      ST=DSIN(TT)
C 
      XE =  CT*XI + ST*YI
      YE = -ST*XI + CT*YI
      ZE =  ZI
C
      RADE = DSQRT( XE**2 + YE**2 + ZE**2 )
      LAME = DATAN2( YE,XE )  
      PHIE = DATAN2( ZE,DSQRT(XE**2+YE**2) )  
C 
      END
C------------------------------------------------------------------------
C
C SUBROUTINE ANGLE COMPUTES THE GREENWICH HOUR ANGLE "TT" IN RADIANS
C GIVEN A TIME "TIME" IN SECONDS RELATIVE TO THE DATE 21/3/1973 0H 0M 0S
C
C AS A RESULT:
C
C  EARTH_FIXED_LONGITUDE = INERTIAL_FIXED_LONGITUDE - HOUR_ANGLE
C
C------------------------------------------------------------------------
      SUBROUTINE ANGLE(TIME,TT )
      IMPLICIT REAL*8 (A-Z)

      COMMON /COMPI/ PI,PI2,RAD,DD,GM 
      DATA STD50,OMT50,OMQ50,EPOCH/100.075542D0,360.985647335D0
     ],.29D-12,8480.0D0/

C     RAD=PI/180.0D0
      T=TIME/86400.0D0+EPOCH
      TT=(STD50+(OMT50+OMQ50*T)*T)*RAD
      TT=MOD(TT,PI2) 
      END
C-----------------------------------------------------------------------
C
C SUBROUTINE ELEVEC
C            ------
C
C FROM ORBITAL ELEMENTS IN VECTOR "E" TO THE KEPLERIAN ELEMENTS IN
C THE VECTOR "S"
C
C THE "E" VECTOR CONTAINS:
C
C E(1) SEMI-MAJOR AXIS E(4) ASCENDING NODE      ( =0 IF I=0 OR I=180 )
C E(2) ECCENTRICITY    E(5) ARGUMENT OF PERIGEE ( =0 IF E=0 )
C E(3) INCLINATION     E(6) TRUE ANOMALY
C
C THE "S" VECTOR CONTAINS THE INERTIAL COORDINATES:
C
C S(1) X  E(4) DX/DT
C S(2) Y  E(5) DY/DT
C S(3) Z  E(6) DZ/DT
C
C THE VARIABLE "GM" EQUALS TO THE CANVENDISH CONSTANT TIMES THE MASS
C OF THE PLANET. THE ORIGINAL SOURCE WAS KINDLY PROVIDED BY K. WAKKER.
C
C-----------------------------------------------------------------------
      SUBROUTINE ELEVEC(S,E)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /COMPI/ PI,PI2,DEGRAD,RADDEG,GM 
      DIMENSION E(6),S(6)
      P    = E(1)*(1D0-E(2)**2)
      CV   = DCOS(E(6))
      ECV  = 1D0+E(2)*CV
      R    = P/ECV
      U    = E(5)+E(6)
      CU   = DCOS(U)
      SU   = DSIN(U)
      CO   = DCOS(E(4))
      SO   = DSIN(E(4))
      CI   = DCOS(E(3))
      SI   = DSIN(E(3))
      COCU = CO*CU
      SOSU = SO*SU
      SOCU = SO*CU
      COSU = CO*SU
      FX   = COCU-SOSU*CI
      FY   = SOCU+COSU*CI
      FZ   = SU*SI
      S(1) = R*FX
      S(2) = R*FY
      S(3) = R*FZ
      F    = DSQRT( GM/P )
      VR   = F*E(2)*DSIN(E(6))
      VU   = F*ECV
      S(4) = VR*FX - VU*(COSU+SOCU*CI)
      S(5) = VR*FY - VU*(SOSU-COCU*CI)
      S(6) = VR*FZ + VU*CU*SI
      RETURN
      END
C-----------------------------------------------------------------------
C
C SUBROUTINE VECELE
C            ------
C
C SUBROUTINE FOR CALCULATING THE CLASSICAL ELEMENTS E(I) FROM THE STATE
C VECTOR S(I), SEE ROUTINE ELEVEC FOR A DESCRIPTION OF THE ARGUMENTS OF
C THIS SUBROUTINE.
C
C-----------------------------------------------------------------------
      SUBROUTINE VECELE(S,E,CH)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /COMPI/ PI,PI2,DEGRAD,RADDEG,GM 
      DIMENSION S(6),E(6),CH(3)
      C1  = S(2)*S(6)-S(3)*S(5)
      C2  = S(3)*S(4)-S(1)*S(6)
      C3  = S(1)*S(5)-S(2)*S(4)
      CH(1)=C1
      CH(2)=C2
      CH(3)=C3
      CC  = C1**2 + C2**2 + C3**2
      C   = DSQRT(CC)
      V02 = S(4)**2 + S(5)**2 + S(6)**2
      R0V0= S(1)*S(4) + S(2)*S(5) + S(3)*S(6)
      R02 = S(1)**2 + S(2)**2 + S(3)**2
      R0  = DSQRT(R02)
      X   = R0*V02/GM
      CX  = CC/GM
      STE = R0V0*C/(R0*GM)
      CTE = CX/R0-1D0
      E(1)= R0/(2D0-X)
      E(2)= DSQRT( STE*STE+CTE*CTE )
      E(3)= DACOS(C3/C)
      EPS = 1.0D-20
      IF ( (DABS(E(3)).LT.EPS) .OR. (DABS(PI-E(3)).LT.EPS) ) THEN
        U= DATAN2(S(2),S(1))
        E(4)=0D0
      ELSE
        U=DATAN2(C*S(3),S(2)*C1-S(1)*C2)
        E(4)=DATAN2(C1,-C2)
      ENDIF
 3    IF (E(2)-EPS) 4,4,5
 4    E(6)= U
      E(5)= 0D0
      GOTO 6
 5    E(6)= DATAN2( STE,CTE )
      E(5)= U-E(6)
 6    IF (E(4).LT.0D0) E(4)=E(4)+PI2 
      IF (E(5).LT.0D0) E(5)=E(5)+PI2 
      IF (E(6).LT.0D0) E(6)=E(6)+PI2 
      RETURN
      END
C------------------------------------------------------------------------
C
C SUBROUTINE TRANSMF
C            -------
C
C ALWAYS GIVEN IS E (ECCENTRICITY)
C
C IF IRC=+1 THEN MEAN ANOMALY TRANSFORMED TO ECCENTRIC AND TRUE ANOMALY
C IF IRC=-1 THEN TRUE ANOMALY TRANSFORMED TO ECCENTRIC AND MEAN ANOMALY
C
C ALL ELEMENTS ARE IN R A D I A N S
C
C------------------------------------------------------------------------
      SUBROUTINE TRANSMF( E,INMEAN,ECCEN,TRUE,IRC )
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /COMPI/ PI,PI2,DEGRAD,RADDEG,GM 
      REAL*8 EPS
      INTEGER MAX

      PARAMETER (EPS=1D-15)
      PARAMETER (MAX=1000)

      REAL*8 ECCEN,MEAN,TRUE,RA,E,TEST,DIFF,CF,SF,R,COSE,SINE,INMEAN
      INTEGER*4 ITR,IRC

      IF (IRC.EQ.+1) THEN
        MEAN=MOD(INMEAN,PI*2D0)
        IF (MEAN.LT.0D0) MEAN=MEAN+2D0*PI
        ECCEN=0D0
        ITR=0
 10     ITR=ITR+1
        ECCEN = MEAN+E*SIN(ECCEN)
        TEST = DABS(ECCEN-E*SIN(ECCEN)-MEAN)
        IF ((ITR.LE.MAX).AND.(TEST.GT.EPS)) GOTO 10
        TRUE =DATAN2( DSQRT(1D0-E**2)*DSIN(ECCEN) , DCOS(ECCEN)-E )
        IF (TRUE.LT.0D0) TRUE=TRUE+2D0*PI
      ELSE IF (IRC.EQ.-1) THEN
        CF=DCOS(TRUE)
        SF=DSIN(TRUE)
        RA = (1-E*E)/(1D0+E*CF)
        COSE = (E+RA*CF)
        SINE = RA*SF/DSQRT(1-E*E)
        ECCEN= DATAN2(SINE,COSE)
        IF (ECCEN.LE.0D0) ECCEN=ECCEN+2D0*PI
        MEAN = ECCEN-E*SIN(ECCEN)
        INMEAN=MEAN
      ELSE
        STOP 'ROUTINE TRANSMF: ILLEGAL IRC SELECTED'
      END IF

      END
      DOUBLE PRECISION FUNCTION DJULIA(JJ,MO,DD) 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER DD 
      LOGICAL LEAP
      COMMON /MONTH/MCAL(12)
      DATA MCAL/0,31,59,90,120,151,181,212,243,273,304,334/ 
      J19=JJ-1900
      NLEAP= J19/4
      LEAP= MOD(J19,4).EQ.0
      DJULIA=2415019.5D0+365.0D0*J19+MCAL(MO)+NLEAP+DD
      IF (LEAP.AND.MO.LT.3) DJULIA=DJULIA-1.0D0
      RETURN
      END 
      SUBROUTINE GETJUL(JJ,MO,DD,HOUR,MIN,SEC,EPOCH)
C THE ROUTINE COMPUTES JULIAN DATE. LAST UPDATE 1990.10.30. 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER DD,HOUR
      IF (JJ.LT.1900) JJ=JJ+1900
      IF (JJ.LT.1957) JJ=JJ+100
      EP= DJULIA(JJ,MO,DD)
      EPOCH= EP+HOUR/24.0D0+MIN/1440D0+SEC/86400.0D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION SIDTIME(TIME)
C THE FUNCTION COMPUTES THE HOUR ANGLE IN RADIANS.
C LAST UPDATE 1990.10.30. 
      IMPLICIT REAL *8 (A-H,O-Z)
      COMMON /COMPI/ PI,PI2,DEGRAD,RADDEG,GM 
      JULIA0=(DJULIA(1994,4,1)-DJULIA(1973,3,21))*86400.0D0 
      ST= 3.3008134458717D0+(TIME-JULIA0)*7.29211585531D-5
      SIDTIME=MOD(ST,PI2) 
      IF (SIDTIME.LT.0.0D0) SIDTIME=SIDTIME+PI2 
      RETURN
      END 
      SUBROUTINE SECULA( AM,EM,IM,C20,AE,DODT,DWDT,DMDT )
C------------------------------------------------------------------------
C
C COMPUTE (4.20) IN (SCHRAMA,1989)
C
C INPUT: (ALL IN METERS,SECONDS,KILOGRAM)
C
C AM:   MEAN SEMI MAJOR AXIS
C EM:   MEAN ECCENTRICITY
C IM:   MEAN INCLINATION IN RADIANS
C GM:   THE GRAVITATIONAL CONSTANT
C C20:  THE UNNORMALIZED FLATTENING PARAMETER
C AE:   THE MEAN EQUATORIAL RADIUS
C
C OUTPUT: (ALL IN RADIANS PER SECOND)
C
C DODT: (4.20A)
C DWDT: (4.20B)
C DMDT: (4.20C)
C
C------------------------------------------------------------------------
C     IMPLICIT NONE
      IMPLICIT REAL*8 (A-H,O-Z) 
      COMMON /COMPI/ PI,PI2,DEGRAD,RADDEG,GM 
      REAL*8 AM,EM,IM
      REAL*8 C20,AE
      REAL*8 DWDT,DMDT,DODT
      REAL*8 NM,DUM,CI, NMDAY 

      CI   = DCOS(IM)
      NM   = DSQRT(GM/AM/AM/AM)
      NMDAY=NM*86400/PI2 
      write(*,*)' n=', nmday,' rev/day' 
      DUM  = DSQRT(1D0-EM*EM)
      DUM  = DUM**3
      DMDT = NM-3D0*NM*C20*AE*AE/4D0/AM/AM/DUM*(3D0*CI**2-1D0)
      DUM  = DSQRT(1D0-EM*EM)
      DUM  = DUM**2
      DODT = 3D0*NM*C20*AE*AE/2D0/AM/AM/DUM*CI
      DWDT = 3D0*NM*C20*AE*AE/4D0/AM/AM/DUM*(1D0-5D0*CI**2)

      END
