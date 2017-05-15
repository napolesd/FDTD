PROGRAM FDTDF90
!Author: Jose Manuel Napoles Duarte (napolesd@cifus.uson.mx, napolesd@gmail.com)
IMPLICIT NONE
!--------------------------------------
INTEGER::n,i,j,T
integer,parameter::ie=500,je=500
real*8,parameter::fp=9.37d14,c=3.d8,theta=65.d0,epsi_lon=2.29d0,pi=4.d0*atan(1.d0)
integer::nsteps,npml,ip,jp
real*8::dist,xdist,ydist,radius
real*8::ddx,dt
real*8::freq_in,w_in,lambda_in,periodo,wp,alambdap,thetar,seno,deltaL,deltaX,an,dnu
real*8::spread,T0,pulse,pulset
real*8::curl_e,del_exp
real*8,dimension(0:ie)::gi3,gi2
real*8,dimension(0:je)::gj3,gj2
real*8,dimension(0:ie)::fi3,fi2,fi1
real*8,dimension(0:je)::fj3,fj2,fj1
real*8,dimension(0:ie,0:je)::hz,ex,ey,dx,dy,ihx,ihy,ga,gb,gc,sx,sy,ix,iy
radius=30
!--------------------------------------
!pi=4.d0*atan(1.d0)
ip=ie/2;jp=je/2
!frequency dependent parameters----------------------
freq_in=0.4d0*fp
w_in=2.d0*pi*freq_in
lambda_in=c/freq_in
periodo=lambda_in/c
t0=20.0d0*(periodo)
spread=10.0d0*(periodo)

wp=2.d0*pi*fp
dnu=wp/100.d0
alambdap=2.d0*pi*c/wp
ddx=alambdap/10.d0
dt=ddx/6.d8
del_exp=dexp(-dt*dnu)
thetar=theta*pi/180.d0
seno=dsin(thetar)
lambda_in=c/freq_in
deltaL=2.0d0
deltaX=lambda_in*deltaL
an=dsqrt(epsi_lon)
!frequency dependent parameters----------------------
!###########inicializacion de los campos############
hz=0.d0;dx=0.d0;dy=0.d0
ex=0.d0;ey=0.d0
ihx=0.d0;ihy=0.d0
ix=0.d0;iy=0.d0
sx=0.d0;sy=0.d0
ga=1.d0;gb=0.d0;gc=0.d0
!###################################################
!:::::::::: PML's parameters initialization:::::::::
gi2=1.d0;gi3=1.d0   
fi1=0.d0;fi2=1.d0;fi3=1.d0
!-------------------------
gj2=1.d0;gj3=1.d0   
fj1=0.d0;fj2=1.d0;fj3=1.d0
!::::::Number of PML cells:::::::
npml=7
CALL PMLS(ie,npml,gi2,gi3,fi1,fi2,fi3)
CALL PMLS(je,npml,gj2,gj3,fj1,fj2,fj3)
!###################################################
!---------Plane definition
!..................................................................
do j=10,je-10
do i=ip+10,ie-10
gb(i,j)=wp*wp*dt/dnu
gc(i,j)=-wp*wp*dt/dnu
enddo
enddo
!..................................................................
!-------Cylinder definition
!..................................................................
!do j=1,je
!do i=1,ie
!xdist = (ip+radius-i)
!ydist = (jp-j)
!dist = sqrt(xdist**2+ydist**2)
!if (dist.le.radius) then
!gb(i,j)=wp*wp*dt/dnu
!gc(i,j)=-wp*wp*dt/dnu
!endif
!enddo
!enddo
!..................................................................


nsteps=10000 ! Temporal steps
T=0
do n=1,nsteps
T=T+1
!calculo de Hz
do j=0,je-1
do i=1,ie
hz(i,j)=gi3(i)*gj3(j)*hz(i,j)- &
gi2(i)*gj2(j)*0.5d0*(ey(i,j)-ey(i-1,j) &
-ex(i,j)+ex(i,j-1))
enddo
enddo

! evanescent source definition
pulset=dexp(- (dble(t*dt-t0)/spread)**2 )
!hz(ip-3,jp)=hz(ip-3,jp)+pulset
do j=10,je-10!
pulse=pulset* &
dsin(2.d0*pi*ddx*an*seno*dble(j)/lambda_in+w_in*dt*T)
hz(ip,j)=hz(ip,j)+ &
pulse*dexp(-(dble(j-jp)*ddx/deltaX)**2)
enddo

!-----------------------
!Se hacen cero los extremos de hz: requerido por las PML
!----arriba y abajo
hz(0:ie-2,0)=0.d0
hz(0:ie-2,je-1)=0.d0
!----izquierda y derecha
hz(0,0:je-2)=0.d0
hz(ie-1,0:je-2)=0.d0
!-------el desplazamiento electrico
do j=1,je
do i=0,ie-1
curl_e=hz(i,j+1)-hz(i,j)
ihx(i,j)=ihx(i,j)+fi1(i)*curl_e
dx(i,j)=fj3(j)*dx(i,j)+ &
fj2(j)*0.5d0*curl_e+ihx(i,j)
enddo
enddo
!--------ahora el desplazamiento dy
do j=1,je
do i=0,ie-1
curl_e=hz(i+1,j)-hz(i,j)
ihy(i,j)=ihy(i,j)+fj1(j)*curl_e
dy(i,j)=fi3(i)*dy(i,j) &
 -fi2(i)*0.5d0*curl_e-ihy(i,j)
enddo
enddo
!El campo electrico Ex
do j=1,je
do i=0,ie-1
ex(i,j)=ga(i,j)*(dx(i,j)-ix(i,j)-del_exp*sx(i,j))
ix(i,j)=ix(i,j)+gb(i,j)*ex(i,j)
sx(i,j)=del_exp*sx(i,j)+gc(i,j)*ex(i,j)
enddo
enddo
!El campo electrico Ey
do j=1,je
do i=0,ie-1
ey(i,j)=ga(i,j)*(dy(i,j)-iy(i,j)-del_exp*sy(i,j))
iy(i,j)=iy(i,j)+gb(i,j)*ey(i,j)
sy(i,j)=del_exp*sy(i,j)+gc(i,j)*ey(i,j)
enddo
enddo
if (mod(T,30).eq.0) then
!if (T.eq.nsteps) then
do j=1,je,2
do i=1,ie,2
write(1,*)hz(i,j)
enddo
enddo
endif
!##########FIELD MONITOR########
!----MONITOR(T,Tin,Tfin,Campo)
CALL MONITOR(T,40, 80,  hz(40,40))
!##################################
!........................
!End of temporal loop 
!........................
enddo
END PROGRAM FDTDF90
!monitor: Stores to a file fort.10
!the value of a field at a given point
!during the time interval [Tin,Tfin]
!guarda en el archivo fort.10
!el CAMPO en un punto, en el intevalo Tin-Tfin
SUBROUTINE MONITOR(T,Tin,Tfin,CAMPO)
INTEGER,INTENT(IN)::T,Tin,Tfin
REAL*8,INTENT(IN)::CAMPO
if (T.gt.Tin.and.T.lt.Tfin) then
write(10,*)T,CAMPO
endif
END SUBROUTINE MONITOR

SUBROUTINE PMLS(imax,npmls,g2,g3,f1,f2,f3)
INTEGER,INTENT(IN)::imax,npmls
real*8,DIMENSION(0:imax),intent(inout)::g2,g3,f1,f2,f3
integer::i,xnum,xd
real*8::xxn,xn
do i=0,npmls
xnum=npmls-i
xd=npmls
xxn=dble(xnum)/xd
xn=0.33d0*(xxn)**3
g2(i)=1.d0/(1.d0+xn)
write(2,*)xn,g2(i)
g2(imax-1-i)=1.d0/(1.d0+xn)
g3(i)=(1.d0-xn)/(1.d0+xn)
g3(imax-i-1)=(1.d0-xn)/(1.d0+xn)
xxn=dble(xnum-0.5d0)/xd
xn=0.25d0*xxn**3
f1(i)=xn
f1(imax-2-i)=xn

f2(i)=1.d0/(1.d0+xn)
f2(imax-2-i)=1.d0/(1.d0+xn)
f3(i)=(1.d0-xn)/(1.d0+xn)
f3(imax-2-i)=(1.d0-xn)/(1.d0+xn)
enddo
END SUBROUTINE PMLS
