        program gra80
c $Id: gra80.for 181 2008-08-14 23:19:51Z tjansson $
c the program computes the gravity anomaly or gravity in GRS80.
c programmed 1999-01-28, last change 1999-05-17 by cct.
        implicit real*8 (a-h,o-z)
      COMMON /CONST/PI,PI4,RADDEG,DEGRAD,D0,D1,D2 
        logical lano,lstop
c
        pi4=atan(1.0d0)
        pi=4.0d0*pi4
        degrad=pi/180.0d0
        raddeg=1.0d0/degrad
        write(*,*)
     *' Calculation of gravity anomaly or gravity, 2003-10-22 '
        write(*,100)
 100    format(' input computation mode: ',/,' 0: normal gravity ',/,
     *  ' 1 gravity from gravity anomaly ',/,
     *  ' 2 gravity anomaly from gravity ')
        read(*,*)mode
        write(*,*)' input 1 for ddmmss, 2 for ddmm, 3 for dd '
        read(*,*)iang
        write(*,*)' Output will be: '
        if (mode.eq.0) then
         write(*,*)' no., lat., long., (deg.), h, '
        write(*,*)
     *' normal gr. at 0 m, H*dg/dh, normal gr. at h (mgal) U m**2/s**2 '
        else
        write(*,*)' no., lat., long., (deg.), h, gravity, '
        write(*,*)
     *' normal gr. at 0 m, delta g, H*dg/dh, (mgal) U m**2/s**2 '
        end if
        if (mode.eq.0) then
         ga=0.0d0
         write(*,*)
     * ' Input no., lat., long., h, lstop '
        else
         if (mode.eq.1) then
          write(*,*)
     *   ' Input no., lat., long., h,  gravity (mgal), lstop '
         else
          write(*,*)
     *   ' Input no., lat., long., h, anomaly (mgal), lstop '
         end if
        end if
c
  10    if (iang.eq.1) then
         if (mode.eq.0) then
          read(*,*)no,nlat,mlat,slat,nlon,mlon,slon,h,lstop
         else
        read(*,*)no,nlat,mlat,slat,nlon,mlon,slon,h,gr,lstop
         end if
        else
        if (iang.eq.2) then
         if (mode.eq.0) then
          read(*,*)no,nlat,slat,nlon,slon,h,lstop
         else
        read(*,*)no,nlat,slat,nlon,slon,h,gr,lstop
         end if
        else
         if (mode.eq.0) then
          read(*,*)no,slat,slon,h,lstop
         else
          read(*,*)no,slat,slon,h,gr,lstop
         end if
        end if
        end if
        call ddeg(nlat,mlat,slat,rlat,IANG)
        call ddeg(nlon,mlon,slon,rlon,IANG)
c
        sinfi = sin(rlat*degrad)
	gamma = 978032.677d0*(1+.00193185135d0*sinfi**2)/
     .  sqrt(1 - .00669438002d0*sinfi**2)
          ha = 0.30877d0*(1-.00142*sinfi**2)*h -.75d-7*h**2
        if (mode.eq.0) then
          ga=gamma-ha
        end if
        if (mode.eq.1) then
          ga = gr-gamma+ha
        end if
        if (mode.eq.2) then
          ga=gr
          gr = gr+gamma-ha
        endif       
        u=6.2636861d7-gr*h*1.0d-5
        if (mode.eq.0) then
         write(*,19)no,rlat,rlon,h,gamma,-ha,ga,u
  19    format(i10,2f12.7,4f10.2,/f12.1)
        else
        write(*,20)no,rlat,rlon,h,gr,gamma,ga,-ha,u
  20    format(i10,2f12.7,4f10.2,/,f7.2,f12.1)
        end if
        if (.not.lstop) go to 10 
        stop
        end
      SUBROUTINE DDEG(IDEG,MIN,SEC,RA,IANG)
C THE SUBROUTINE CONVERTS FOR IANG = 1,2,3,4 ANGLES IN (1) DEGREES, MI-
C NUTES, SECONDS, (2) DEGREES, MINUTES, (3) DEGREES AND (4) 400-DEGREES
C TO DECIMAL DEGREES. 
C IF DOUBLE PRECISION IS NEEDED, ACTIVATE:
      IMPLICIT INTEGER(I,J,K,M,N),REAL *8(A-H,O-Z)
      I = 1
      IF (IDEG .LT. 0 .AND. IANG .LT. 3) I = -1
      GO TO (1,2,3,4,3),IANG
    1 J = 1
      IF (MIN.LT.0) J = -1
      SE =I*IDEG*3600+J*MIN*60+SEC
      I = J*I
      GO TO 5
    2 SE=I*IDEG*3600+SEC*60
      GO TO 5
    3 RA = SEC
      RETURN 
    4 SE = SEC*3240
    5 RA= I*SE/3600.0D0  
      RETURN
      END
