      program degv
c $Id: degv.for 295 2011-07-25 21:36:55Z cct $
c the program computes gravity anomaly degree-variances and error-degree
c variances from a set of spherical harmonic coefficients.
c the degree-variances will refer to the mean-earth radius.
c programmed 2001-02-09 by cct. updated 2009.
c input: file names of input coefficients, output gravity anomaly
c degree-variances and error-degree variances. The degree variances
c are in units of mgal**2 and refer to the mean-earth sphere.
c the coefficients must be on the form: degree (i), order (j),
c Cij, Sij, error of Cij, Sij. The coefficients must be fully
c normalized and unitless.
c 
      implicit  none                
      integer imax, i,j
      real*8 gm,RE,cij,sij,ecij,esij,u20,u40,gvar,evar,hvar,hevar,
     *a,degva(0:2190),edegv(0:2190)
      logical ltab
      character *128 line,egmf,dgvf,edgvf
c
      write(*,*)' Program degv, ver.2011-07-25 '
      write(*,*)' Input name of file with EGM coefficients '
      read(*,'(a)')egmf
      write(*,*)egmf
      open(9,file=egmf)
      write(*,*)' Input name of file to hold degree-variances '
      read(*,'(a)')dgvf
      write(*,*)dgvf
      write(*,*)' Input name of file to hold error-degree variances '
      read(*,'(a)')edgvf
      write(*,*)edgvf
      write(*,*)' output on tabular form ? (T/F) '
      read(*,*)ltab
      write(*,*)ltab
c must be activated if file contains records with text before coefficients.
c     do i=1,12
c     read(9,'(a)')line
c     end do
      gm=3.986005d14
C mean radius.
      RE=6371000.0d0
      a=6378137.0d0
      u20=-0.484165371736E-03 
      u40=0.539873863789E-06 
      gvar=0.0d0
      evar=0.0d0
      hvar=0.0d0
      hevar=0.0d0
      do i=1,2190
       degva(i)=0.0d0
       edegv(i)=0.0d0
      end do
      imax=0
c
  10  read(9,*,end=99)i,j,cij,sij,ecij,esij
c100  format(2i5,2d20.12)
      if (i.gt.imax) imax=i
c we subtract the degree 2,0 and 4,0 terms, since the degree-variances
c are for the anomalous potential (T) and not for the full potential (W).
      if (i.eq.2.and.j.eq.0) then
       cij=cij-u20
       write(*,251)i,j,cij
 251   format(' i,j,cij ',2i3,d15.6)
      end if
      if (i.eq.4.and.j.eq.0) then
       cij=cij-u40
       write(*,251)i,j,cij
      end if
      degva(i)=degva(i)+(cij**2+sij**2)*((a/RE)**i)*(1.0d5*(i-1)
     *  *GM/RE**2)**2
      edegv(i)=edegv(i)+(ecij**2+esij**2)*(1.0d5*(i-1)*GM/RE**2
     **(a/RE)**i)**2
c result in mgal**2.
      go to 10
  99  close(9)
      write(*,*)' coefficients input finished, imax= ',imax
      do i=2,imax
       gvar=gvar+degva(i)
       hvar=hvar+(RE**3/GM)**2/((i-1)**2)*degva(i)*1.0D-10
       evar=evar+edegv(i)
       hevar=hevar+(RE**3/GM)**2/((i-1)**2)*edegv(i)*1.0D-10
      end do
      write(*,250)gvar,evar,hvar,hevar
 250  format(' gravity anomaly variance = ',f12.2,/
     *       ' error-variance           = ',f12.2, ' mgal**2 ',/
     *       ' geoid variance           = ',f12.4, ' m**2 ',/
     *       ' geoid error variance     = ',f12.4, ' m**2 ')
      open(8,file=dgvf)
      open(10,file=edgvf)
      if (ltab) then
       do i=2,imax
        write(8,201)i,degva(i)
        write(10,201)i,edegv(i)
 201    format(i5,e14.6)
       end do
      else
       write(8,200)(degva(i),i=2,imax)
       write(10,200)(edegv(i),i=2,imax)
 200   format(5d14.6)
      end if
      close(8)
      close(10)
      stop
      end

      
