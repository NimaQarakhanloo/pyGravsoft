	program selecgr
c programmed June 2011 by D.Arabelos. Last modification 2011-07-21 by cct.
	implicit real*8 (a-h,o-z)
c selection of satellite data close to the nodes of an arbitrary grid
c the arbitrary grid could be an equal-area 1x1 or 2x2 etc. grid
	dimension phi(20000000),dla(20000000),trr(20000000),h(20000000),
     .            n(20000000),trrd(20000000),err(20000000)
	dimension phip(100000),dlap(100000),trrp(100000),hp(100000),np(100000),
     .            dp(100000),trrdp(100000),errp(100000)
	dimension m(50000),phig(50000),dlag(50000)
        character*72 ifile1,ifile2,ofile
        write(*,*)' selecgr, 2011-06-06 '
        write(*,*)
     *  ' input names of grid-file, data-file and output file '
        read(*,'(a)') ifile1
        read(*,'(a)') ifile2
        read(*,'(a)') ofile
        write(*,*)ifile1,ifile2,ofile
 	open(1,file=ifile1,status='old')
 	open(3,file=ifile2,status='old')
  	open(2,file=ofile,status='unknown')
 	open(4,file='store')
        open(25,file='grerr.dat')
        write(*,*)' suspected errors in grerr.dat '
        mmax=50000
        ngrerr=0
        rlatmin=1.0d10
        rlatmax=-1.0d10
        rlonmin=1.0d10
        rlonmax=-1.0d10
c read the grid
	i=1
1	read(1,*,end=99) m(i),phig(i),dlag(i)
	i=i+1
        if (i.gt.mmax) stop
	goto 1
99	ig=i-1
c ig is the number of the grid points
c read the satellite data
	i=1
2	read(3,*,end=88) n(i),phi(i),dla(i),h(i),trr(i),trrd(i),err(i)
        if (dla(i).lt.0.0) dla(i)=dla(i)+360.0d0
        if (phi(i).gt.rlatmax) rlatmax=phi(i)
        if (phi(i).lt.rlatmin) rlatmin=phi(i)
        if (dla(i).gt.rlonmax) rlonmax=dla(i)
        if (dla(i).lt.rlonmin) rlonmin=dla(i)
        if (abs(trrd(i)).gt.3.0) then
         write(25,110)n(i),trrd(i)
 110     format(i10,f12.3)
         ngrerr=ngrerr+1
        else
	 i=i+1
        end if
	goto 2
88	is=i-1
c is is the number of satellite points
	write(*,*) 'number of grid points = ',ig,' number of 
     . satellite points = ',is,' errors = ',ngrerr
        write(*,*)' min,max lat. lon ',rlatmin,rlatmax,rlonmin,rlonmax
c end of section reading the input
c start the selection section
	it=0
c dlarge defines the max allowed collection radius
c	dlarge=30.d0
 	dlarge=50.d0
	write(*,*) 'no points selected at dist. larger than km ',dlarge
	do i=1,ig
c first define the bounds of each 1-d blocks
 	phi1=phig(i)-0.5 ! degree
	phi2=phig(i)+0.5 ! degree
 	dla1=dlag(i)-360./m(i)/2. ! degree
	dla2=dlag(i)+360./m(i)/2  ! degree
	ic=1
	do 3 j=1,is
	if(phi(j).lt.phi1.or.phi(j).gt.phi2) goto 3
	if(dla(j).lt.dla1.or.dla(j).gt.dla2) goto 3
	phip(ic)=phi(j)
	dlap(ic)=dla(j)
	np(ic)=n(j)
	hp(ic)=h(j)
	trrp(ic)=trr(j)
        trrdp(ic)=trrd(j)
        errp(ic)=err(j)
	ic=ic+1
3	continue
	ic=ic-1
c from all points into the area phi1-dla1-phi2-dla2 select the closest to grid node
	dmin=50.d0
	nmin=0
	do j=1,ic
	 call sphdist(phig(i),phip(j),dlag(i),dlap(j),dd)
	 if(dd.lt.dmin) then
	  dmin=dd
	  nmin=j
	 endif
	enddo
c 
	if(nmin.gt.0.and.dmin.lt.dlarge) then
c it is the total number of the selected points
	it=it+1
c	write(*,100) np(nmin),phip(nmin),dlap(nmin),hp(nmin),trrp(nmin),dmin
	write(4,100) np(nmin),phip(nmin),dlap(nmin),hp(nmin),trrp(nmin),
     *  trrdp(nmin),errp(nmin),dmin
	endif
c
	enddo
c 
	write(*,*) ' the number of selected points is ',it
100	format(i10,2f13.6,f13.1,f12.4,4f10.3)
c arrange the output data along their tracks, according to time (first column)
	rewind (4)
	do i=1,it
 	read(4,*) n(i),phi(i),dla(i),h(i),trr(i),trrd(i),err(i),dmin
	enddo
c
	do i=1,it
	do j=i,it
	if(n(i).gt.n(j)) then
	kk=n(j)
	n(j)=n(i)
	n(i)=kk
c
	pph=phi(j)
	phi(j)=phi(i)
	phi(i)=pph
c
	ddl=dla(j)
	dla(j)=dla(i)
	dla(i)=ddl
c
	hh=h(j)
	h(j)=h(i)
	h(i)=hh
c
	ttrr=trr(j)
	trr(j)=trr(i)
	trr(i)=ttrr
c
        trrdd=trrd(j)
        trrd(j)=trrd(i)
        trrd(i)=trrdd
C
        eerr=err(j)
        err(j)=err(i)
        err(i)=eerr
c
	endif
	enddo
	enddo
c output the shorted data
        write(*,*)' output to file '
	do i=1,it	
	write(2,100) n(i),phi(i),dla(i),h(i),trr(i),trrd(i),err(i)
	enddo
	end	
c
	subroutine sphdist(phi1,phi2,dla1,dla2,d)
c computation of the spherical distance between two points
	implicit real*8(a-h,o-z)
	pi=4.d0*datan(1.d0)
	rho=pi/180.d0
	r=6371.d0
	dphi2=(phi2-phi1)/2.d0
	ddla2=(dla2-dla1)/2.d0
	d=r*dacos(1.d0-2.d0*(dsin(dphi2*rho)**2+
     .    dcos(phi1*rho)*dcos(phi2*rho)*dsin(ddla2*rho)**2))
	return
	end
	
