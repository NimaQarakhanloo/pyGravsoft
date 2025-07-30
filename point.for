      program point
c $Id: point.for 181 2008-08-14 23:19:51Z tjansson $
c Program for determination of a position of a point in a plane
c projection given the coordinates of one or more points as well as
c distances and/or angular measurements.
c it is planned to extend the program with the calculation of the
c point uncertainty.
c Reference: Kahmen, H. & W.Faig: Surveying, 1988.
c Programmed by C.C.Tscherning, Feb. 1998, modified march 20, 1998..
c
      implicit real*8 (a-h,o-z)
      logical lcont,lutm
      D1=1.0D0
      D2=2.0D0
      pi4=atan(1.0d0)
      pi2=pi4*2.0d0
      pi=4.0d0*pi4
      degrad=pi/180.0d0
      radeg=1.0d0/degrad
      write(*,*) ' Program Point, ver. feb. 1998 '
c
      write(*,*)' If coordinates are UTM type t, if Gauss/Kr. type f '
      read(*,*)lutm
c
      write(*,*)' type '
      write(*,*)' 1 for intersection (fremskaering)'
      write(*,*)' 2 for resection (tilbageskaering), '
      write(*,*)' 3 for arc section (bueskaering, trilateration)'
      write(*,*)' 4 polar method (polaer bestemmelse), '
c     write(*,*)' 5 for combined angle and distance method '
  50  read(*,*) method
c
      write(*,*)' Input 1 for angles in ddmmss, 2 for ddmm,'
      write(*,*)' 3 for dd and 4 for cc (gon) '
      read(*,*)iang 
      if (lutm) then
         write(*,*)' scale factor 0.9996 and false E assumed '
         write(*,*)' Distance corrections will be applied '
      else
         write(*,*)' All distances are assumed to be scale-corrected '
      end if
      write(*,*)' All distances must be slope corrected. '
      write(*,*)' '
      write(*,*)' input (N, E) for point A '
      read(*,*)an,ae
      write(*,10)an,ae
      write(*,*)' input (N, E) for point B '
      read(*,*)bn,be
      write(*,10)bn,be
c
      alfa= pi2-atan2(bn-an,be-ae)
      dalfa=alfa*radeg
      c2=((bn-an)**2+(be-ae)**2)
      c=sqrt(c2)
      write(*,60)dalfa,c
  60  format(' azimuth from A to B = ',f12.7,
     *' distance = ',f12.3,' m')
      if (lutm) then
      aeu=ae-500000.0d0
      beu=be-500000.0d0
      scale=(d1-(aeu**2+beu**2+aeu*beu)/(6.0d0*6371000.0d0**2))/0.9996D0
      write(*,80)scale
   80 format(' scale factor used ',f12.8)
      else
      aeu=ae
      beu=be
      scale=D1
      end if
      dh=0.0d0
c
      go to (100,200,300,400,500),method
c
 100  write(*,*)' Intersection (fremskaering)'
      write(*,*)' input <CAB '
      CAB=RA(IDEG,MIN,SEC,IANG)
      write(*,95)CAB
  95  format(f12.7)
      write(*,*)' input <ABC) '
      ABC=RA(IDEG,MIN,SEC,IANG)
      write(*,95)ABC
      RCAB=CAB*degrad
      RABC=ABC*degrad
      write(*,11)dalfa
c test of equations on p. 144 in Borre: Landmaaling, 1990.
c     if (abs(RCAB-RABC).lt.0.0D-7) go to 5
c     de=ae-sin(RABC)*sin(alfa+RCAB)*c/sin(RCAB-RABC)
c     dn=an-sin(RABC)*cos(alfa+RCAB)*c/sin(RCAB-RABC)
c     write(*,10)dn,de,dh
c     dn=an-sin(RABC)*sin(alfa+RCAB)*c/sin(RCAB-RABC)
c     de=ae-sin(RABC)*cos(alfa+RCAB)*c/sin(RCAB-RABC)
c     write(*,10)dn,de,dh
c     dn=an-sin(RCAB)*sin(alfa+RCAB)*c/sin(RCAB-RABC)
c     de=ae-sin(RABC)*cos(alfa+RABC)*c/sin(RCAB-RABC)
c     write(*,10)dn,de,dh
c  5  alfaac=alfa+RCAB
      alfaac=alfa+RCAB
      alfabc=-pi+alfa-RABC
      if (abs(alfaac-pi2).lt.0.0d-7) go to 999
      if (abs(alfabc-pi2).lt.0.0d-7) go to 999
      write(*,75)radeg*alfaac,radeg*alfabc
  75  format(' azimuth AC= ',f12.7,' BC = ',f12.7)
      fn=((be-ae)-(bn-an)*tan(alfabc))/(tan(alfaac)-tan(alfabc))
      dn=an+fn
      de=ae+fn*tan(alfaac)
      go to 999
c
 200  write(*,*)' Resection (tilbageskaering), '
      write(*,*)' input (N, E) for point C '
      read(*,*)cn,ce
      write(*,10)cn,ce
      write(*,*)' unknown point is named D '
      write(*,*)' input <ADB '
      ADB=RA(IDEG,MIN,SEC,IANG)
      write(*,95)ADB
      radb=adb*degrad
      write(*,*)' input < ADC '
      ADC=RA(IDEG,MIN,SEC,IANG)
      write(*,95)ADC
      radc=adc*degrad
      anc=an-cn
      aec=ae-ce
      bnc=bn-cn
      bec=be-ce
      tanacd=(aec/tan(radb)-bec/tan(radc)+(bnc-anc))/
     *(anc/tan(radb)-bnc/tan(radc)-(bec-aec))
      write(*,*)tanacd
      en=((aec+anc/tan(radb))*tanacd+anc-aec/tan(radb))/
     *(D1+tanacd**2)
      dn=en+cn
      ee=en*tanacd
      de=ee+ce
      go to 999
c
 300  write(*,*)' Arc section (bueskaering, trilateration)'
      write(*,*)' input observed distance (BC)=a '
      read(*,*)a
      write(*,*)' input observed distance (AC)=b '
      read(*,*)b
      write(*,96)a,b
      a=a*scale
      b=b*scale
  96  format(' a = ',f12.3,' b = ',f12.3)
      a1=acos((b**2+c2-a**2)/(D2*b*c))
      a2=acos((a**2+c2-b**2)/(D2*a*c))
      write(*,70)a1*radeg,a2*radeg
  70  format(' azimuth AC = ',f12.7,' and BC =',f12.7)
      de=ae+sin(alfa+a1)*b*scale
      dn=an+cos(alfa+a1)*b*scale
      go to 999
c
 400  write(*,*)' Polar Method (polaer bestemmelse), '
      write(*,*)' input observed distance (AC)=b '
      read(*,*)b
      write(*,97)b
  97  format(' b = ',f12.3)
      b=b*scale
      write(*,*)' input <BAC '
      BAC=RA(IDEG,MIN,SEC,IANG)
      write(*,95)BAC
      RBAC=BAC*degrad  
      write(*,11)dalfa
      dn=an+b*scale*cos(alfa+RBAC)
      de=ae+b*scale*sin(alfa+RBAC)
  11  format(' azimuth from A to B = ',f13.8)
      go to 999
c
 500  a=1
c     NOT YET IMPLEMENTED
c
 999  write(*,10)dn,de
  10  format(2f12.3,f8.3)
      write(*,*)' '
c
      write(*,*)' if new calculation type t else f '
      read(*,*)lcont
      if (lcont) go to 50
      stop
      end
c
      DOUBLE PRECISION FUNCTION RA(IDEG,MIN,SEC,IANG)
C THE SUBROUTINE CONVERTS FOR IANG = 1,2,3,4 ANGLES IN (1) DEGREES, MI-
C NUTES, SECONDS, (2) DEGREES, MINUTES, (3) DEGREES AND (4) 400-DEGREES
C TO DECIMAL DEGREES. 
C IF DOUBLE PRECISION IS NEEDED, ACTIVATE:
      IMPLICIT INTEGER(I,J,K,M,N),REAL *8(A-H,O-Z)
      I = 1
      GO TO (1,2,3,4,3),IANG
    1 J = 1
      READ(*,*)IDEG,MIN,SEC
      IF (IDEG .LT. 0 .AND. IANG .LT. 3) I = -1
      IF (MIN.LT.0) J = -1
      SE =I*IDEG*3600+J*MIN*60+SEC
      I = J*I
      GO TO 5
    2 READ(*,*)IDEG,SEC
      IF (IDEG .LT. 0 .AND. IANG .LT. 3) I = -1
      SE=I*IDEG*3600+SEC*60
      GO TO 5
    3 READ(*,*)SEC
      RA = SEC
      RETURN   
    4 READ(*,*)SEC
      SE = SEC*3240
    5 RA= I*SE/3600.0D0  
      RETURN
      END
