module interpol_mod

  use par_mod, only: nzmax, maxspec

  implicit none

  real :: uprof(nzmax),vprof(nzmax),wprof(nzmax),wprofeta(nzmax),detaprof(nzmax)
  real :: usigprof(nzmax),vsigprof(nzmax),wsigprof(nzmax),wsigprofeta(nzmax)
  real :: rhoprof(nzmax),rhogradprof(nzmax)

  real :: u,v,w,usig,vsig,wsig,pvi,ueta,veta,weta,wsigeta

  real :: p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2
  integer :: ix,jy,ixp,jyp,ngrid,indz,indzp,indzeta,indzpeta
  integer :: induv,indpuv
  logical :: depoindicator(maxspec)
  logical :: indzindicator(nzmax)

!$OMP THREADPRIVATE(uprof,vprof,wprof,usigprof,vsigprof,wsigprof, &
!$OMP rhoprof,rhogradprof,u,v,w,usig,vsig,wsig,pvi, &
!$OMP p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2,ix,jy,ixp,jyp, &
!$OMP ngrid,indz,indzp,depoindicator,indzindicator, &
!$OMP wprofeta,wsigprofeta,induv,indpuv, &
!$OMP indzeta,indzpeta,ueta,veta,weta,wsigeta,detaprof)

contains

! subroutine interpol_all(itime,xt,yt,zt,zteta)
!   !                          i   i  i  i
!   !*****************************************************************************
!   !                                                                            *
!   !  This subroutine interpolates everything that is needed for calculating the*
!   !  dispersion.                                                               *
!   !                                                                            *
!   !    Author: A. Stohl                                                        *
!   !                                                                            *
!   !    16 December 1997                                                        *
!   !                                                                            *
!   !  Revision March 2005 by AST : all output variables in common block cal-    *
!   !                               culation of standard deviation done in this  *
!   !                               routine rather than subroutine call in order *
!   !                               to save computation time                     *
!   !                                                                            *
!   !*****************************************************************************
!   !                                                                            *
!   ! Variables:                                                                 *
!   ! itime [s]          current temporal position                               *
!   ! memtime(3) [s]     times of the wind fields in memory                      *
!   ! xt,yt,zt           coordinates position for which wind data shall be       *
!   !                    culated                                                 *
!   !                                                                            *
!   ! Constants:                                                                 *
!   !                                                                            *
!   !*****************************************************************************

!   use par_mod
!   use com_mod
!   use interpol_mod
!   use hanna_mod

!   implicit none

!   integer :: itime
!   real :: xt,yt,zt
!   real:: zteta

!   ! Auxiliary variables needed for interpolation
!   real :: ust1(2),wst1(2),oli1(2),oliaux
!   real :: y1(2),y2(2),y3(2),rho1(2),rhograd1(2),psint(2)
!   real :: usl,vsl,wsl,usq,vsq,wsq,xaux,deta1(2),detatemp
!   integer :: i,m,n,indexh
!   real,parameter :: eps=1.0e-30


!   !********************************************
!   ! Multilinear interpolation in time and space
!   !********************************************

!   ! Determine the lower left corner and its distance to the current position
!   !*************************************************************************

!   ddx=xt-real(ix)
!   ddy=yt-real(jy)
!   rddx=1.-ddx
!   rddy=1.-ddy
!   p1=rddx*rddy
!   p2=ddx*rddy
!   p3=rddx*ddy
!   p4=ddx*ddy

!   ! Calculate variables for time interpolation
!   !*******************************************

!   dt1=real(itime-memtime(1))
!   dt2=real(memtime(2)-itime)
!   dtt=1./(dt1+dt2)


!   !*****************************************
!   ! 1. Interpolate u*, w* and Obukhov length
!   !*****************************************

!   ! a) Bilinear horizontal interpolation

!   do m=1,2
!     indexh=memind(m)

!     ust1(m)=p1*ustar(ix ,jy ,1,indexh) &
!          + p2*ustar(ixp,jy ,1,indexh) &
!          + p3*ustar(ix ,jyp,1,indexh) &
!          + p4*ustar(ixp,jyp,1,indexh)
!     wst1(m)=p1*wstar(ix ,jy ,1,indexh) &
!          + p2*wstar(ixp,jy ,1,indexh) &
!          + p3*wstar(ix ,jyp,1,indexh) &
!          + p4*wstar(ixp,jyp,1,indexh)
!     oli1(m)=p1*oli(ix ,jy ,1,indexh) &
!          + p2*oli(ixp,jy ,1,indexh) &
!          + p3*oli(ix ,jyp,1,indexh) &
!          + p4*oli(ixp,jyp,1,indexh)
!   end do

!   ! b) Temporal interpolation

!   ust=(ust1(1)*dt2+ust1(2)*dt1)*dtt
!   wst=(wst1(1)*dt2+wst1(2)*dt1)*dtt
!   oliaux=(oli1(1)*dt2+oli1(2)*dt1)*dtt

!   if (oliaux.ne.0.) then
!     ol=1./oliaux
!   else
!     ol=99999.
!   endif


!   !*****************************************************
!   ! 2. Interpolate vertical profiles of u,v,w,rho,drhodz
!   !*****************************************************


!   ! Determine the level below the current position
!   !***********************************************

!   do i=2,nz
!     if (height(i).gt.zt) then
!       indz=i-1
!       indzp=i
!       exit
!     endif
!   end do

!   !**************************************
!   ! 1.) Bilinear horizontal interpolation
!   ! 2.) Temporal interpolation (linear)
!   !**************************************

!   ! Loop over 2 time steps and indz levels
!   !***************************************

!   do n=indz,indzp
!     wsl=0.
!     wsq=0.
!     do m=1,2
!       indexh=memind(m)

!       y3(m)=p1*ww(ix ,jy ,n,indexh) &
!            +p2*ww(ixp,jy ,n,indexh) &
!            +p3*ww(ix ,jyp,n,indexh) &
!            +p4*ww(ixp,jyp,n,indexh)
!       wsl=wsl+ww(ix ,jy ,n,indexh)+ww(ixp,jy ,n,indexh) &
!            +ww(ix ,jyp,n,indexh)+ww(ixp,jyp,n,indexh)
!       wsq=wsq+ww(ix ,jy ,n,indexh)*ww(ix ,jy ,n,indexh)+ &
!            ww(ixp,jy ,n,indexh)*ww(ixp,jy ,n,indexh)+ &
!            ww(ix ,jyp,n,indexh)*ww(ix ,jyp,n,indexh)+ &
!            ww(ixp,jyp,n,indexh)*ww(ixp,jyp,n,indexh)
!     end do
!     wprof(n)=(y3(1)*dt2+y3(2)*dt1)*dtt
!     indzindicator(n)=.false.

!   ! Compute standard deviations
!   !****************************

!     xaux=wsq-wsl*wsl/8.
!     if (xaux.lt.eps) then
!       wsigprof(n)=0.
!     else
!       wsigprof(n)=sqrt(xaux/7.)
!     endif

!   end do


!   ! Same for zt in eta coordinates
!   !*******************************
!   indzeta=nz-1
!   indzpeta=nz
!   do i=2,nz
!     if (wheight(i).lt.zteta) then
!       indzeta=i-1
!       indzpeta=i
!       exit
!     endif
!   end do

!   induv=nz-1
!   indpuv=nz
!   do i=2,nz
!     if (uvheight(i).lt.zteta) then
!       induv=i-1
!       indpuv=i
!       exit
!     endif
!   end do

!   !**************************************
!   ! 1.) Bilinear horizontal interpolation
!   ! 2.) Temporal interpolation (linear)
!   !**************************************

!   ! Loop over 2 time steps and indz levels
!   !***************************************

!   do n=induv,indpuv
!     usl=0.
!     vsl=0.
!     usq=0.
!     vsq=0.
!     do m=1,2
!       indexh=memind(m)
!       if (ngrid.lt.0) then
!         y1(m)=p1*uupoleta(ix ,jy ,n,indexh) &
!              +p2*uupoleta(ixp,jy ,n,indexh) &
!              +p3*uupoleta(ix ,jyp,n,indexh) &
!              +p4*uupoleta(ixp,jyp,n,indexh)
!         y2(m)=p1*vvpoleta(ix ,jy ,n,indexh) &
!              +p2*vvpoleta(ixp,jy ,n,indexh) &
!              +p3*vvpoleta(ix ,jyp,n,indexh) &
!              +p4*vvpoleta(ixp,jyp,n,indexh)
!         usl=usl+uupoleta(ix ,jy ,n,indexh)+uupoleta(ixp,jy ,n,indexh) &
!              +uupoleta(ix ,jyp,n,indexh)+uupoleta(ixp,jyp,n,indexh)
!         vsl=vsl+vvpoleta(ix ,jy ,n,indexh)+vvpoleta(ixp,jy ,n,indexh) &
!              +vvpoleta(ix ,jyp,n,indexh)+vvpoleta(ixp,jyp,n,indexh)

!         usq=usq+uupoleta(ix ,jy ,n,indexh)*uupoleta(ix ,jy ,n,indexh)+ &
!              uupoleta(ixp,jy ,n,indexh)*uupoleta(ixp,jy ,n,indexh)+ &
!              uupoleta(ix ,jyp,n,indexh)*uupoleta(ix ,jyp,n,indexh)+ &
!              uupoleta(ixp,jyp,n,indexh)*uupoleta(ixp,jyp,n,indexh)
!         vsq=vsq+vvpoleta(ix ,jy ,n,indexh)*vvpoleta(ix ,jy ,n,indexh)+ &
!              vvpoleta(ixp,jy ,n,indexh)*vvpoleta(ixp,jy ,n,indexh)+ &
!              vvpoleta(ix ,jyp,n,indexh)*vvpoleta(ix ,jyp,n,indexh)+ &
!              vvpoleta(ixp,jyp,n,indexh)*vvpoleta(ixp,jyp,n,indexh)
!       else
!         y1(m)=p1*uueta(ix ,jy ,n,indexh) &
!              +p2*uueta(ixp,jy ,n,indexh) &
!              +p3*uueta(ix ,jyp,n,indexh) &
!              +p4*uueta(ixp,jyp,n,indexh)
!         y2(m)=p1*vveta(ix ,jy ,n,indexh) &
!              +p2*vveta(ixp,jy ,n,indexh) &
!              +p3*vveta(ix ,jyp,n,indexh) &
!              +p4*vveta(ixp,jyp,n,indexh)
!         usl=usl+uueta(ix ,jy ,n,indexh)+uueta(ixp,jy ,n,indexh) &
!              +uueta(ix ,jyp,n,indexh)+uueta(ixp,jyp,n,indexh)
!         vsl=vsl+vveta(ix ,jy ,n,indexh)+vveta(ixp,jy ,n,indexh) &
!              +vveta(ix ,jyp,n,indexh)+vveta(ixp,jyp,n,indexh)

!         usq=usq+uueta(ix ,jy ,n,indexh)*uueta(ix ,jy ,n,indexh)+ &
!              uueta(ixp,jy ,n,indexh)*uueta(ixp,jy ,n,indexh)+ &
!              uueta(ix ,jyp,n,indexh)*uueta(ix ,jyp,n,indexh)+ &
!              uueta(ixp,jyp,n,indexh)*uueta(ixp,jyp,n,indexh)
!         vsq=vsq+vveta(ix ,jy ,n,indexh)*vveta(ix ,jy ,n,indexh)+ &
!              vveta(ixp,jy ,n,indexh)*vveta(ixp,jy ,n,indexh)+ &
!              vveta(ix ,jyp,n,indexh)*vveta(ix ,jyp,n,indexh)+ &
!              vveta(ixp,jyp,n,indexh)*vveta(ixp,jyp,n,indexh)
!       endif
!       rhograd1(m)=p1*drhodzeta(ix ,jy ,n,indexh) &
!            +p2*drhodzeta(ixp,jy ,n,indexh) &
!            +p3*drhodzeta(ix ,jyp,n,indexh) &
!            +p4*drhodzeta(ixp,jyp,n,indexh)
!       rho1(m)=p1*rhoeta(ix ,jy ,n,indexh) &
!            +p2*rhoeta(ixp,jy ,n,indexh) &
!            +p3*rhoeta(ix ,jyp,n,indexh) &
!            +p4*rhoeta(ixp,jyp,n,indexh)
!     end do
!     uprof(n)=(y1(1)*dt2+y1(2)*dt1)*dtt
!     vprof(n)=(y2(1)*dt2+y2(2)*dt1)*dtt
!     rhoprof(n)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
!     rhogradprof(n)=(rhograd1(1)*dt2+rhograd1(2)*dt1)*dtt

!     xaux=usq-usl*usl/8.
!     if (xaux.lt.eps) then
!       usigprof(n)=0.
!     else
!       usigprof(n)=sqrt(xaux/7.)
!     endif

!     xaux=vsq-vsl*vsl/8.
!     if (xaux.lt.eps) then
!       vsigprof(n)=0.
!     else
!       vsigprof(n)=sqrt(xaux/7.)
!     endif

!   end do

!   do n=indzeta,indzpeta
!     wsl=0.
!     wsq=0.
!     do m=1,2
!       indexh=memind(m)

!       y3(m)=p1*wweta(ix ,jy ,n,indexh) &
!            +p2*wweta(ixp,jy ,n,indexh) &
!            +p3*wweta(ix ,jyp,n,indexh) &
!            +p4*wweta(ixp,jyp,n,indexh)

!       wsl=wsl+wweta(ix ,jy ,n,indexh)+wweta(ixp,jy ,n,indexh) &
!            +wweta(ix ,jyp,n,indexh)+wweta(ixp,jyp,n,indexh)
!       wsq=wsq+wweta(ix ,jy ,n,indexh)*wweta(ix ,jy ,n,indexh)+ &
!            wweta(ixp,jy ,n,indexh)*wweta(ixp,jy ,n,indexh)+ &
!            wweta(ix ,jyp,n,indexh)*wweta(ix ,jyp,n,indexh)+ &
!            wweta(ixp,jyp,n,indexh)*wweta(ixp,jyp,n,indexh)

!     end do
!     wprofeta(n)=(y3(1)*dt2+y3(2)*dt1)*dtt

!     xaux=wsq-wsl*wsl/8.
!     if (xaux.lt.eps) then
!       wsigprof(n)=0.
!     else
!       wsigprof(n)=sqrt(xaux/7.)
!     endif

!   end do
! end subroutine interpol_all

! subroutine interpol_misslev(n)
!   !                            
!   !*****************************************************************************
!   !                                                                            *
!   !  This subroutine interpolates u,v,w, density and density gradients.        *
!   !                                                                            *
!   !    Author: A. Stohl                                                        *
!   !                                                                            *
!   !    16 December 1997                                                        *
!   !    Update: 2 March 1999                                                    *
!   !                                                                            *
!   !  Revision March 2005 by AST : all output variables in common block cal-    *
!   !                               culation of standard deviation done in this  *
!   !                               routine rather than subroutine call in order *
!   !                               to save computation time                     *
!   !                                                                            *
!   !*****************************************************************************
!   !                                                                            *
!   ! Variables:                                                                 *
!   ! n                  level                                                   *
!   !                                                                            *
!   ! Constants:                                                                 *
!   !                                                                            *
!   !*****************************************************************************

!   use par_mod
!   use com_mod
!   use interpol_mod
!   use hanna_mod

!   implicit none

!   ! Auxiliary variables needed for interpolation
!   real :: y1(2),y2(2),y3(2),rho1(2),rhograd1(2)
!   real :: usl,vsl,wsl,usq,vsq,wsq,xaux,psint(2),psint_t
!   integer :: m,n,indexh
!   real,parameter :: eps=1.0e-30


!   !********************************************
!   ! Multilinear interpolation in time and space
!   !********************************************


!   !**************************************
!   ! 1.) Bilinear horizontal interpolation
!   ! 2.) Temporal interpolation (linear)
!   !**************************************

!   ! Loop over 2 time steps
!   !***********************

!   usl=0.
!   vsl=0.
!   wsl=0.
!   usq=0.
!   vsq=0.
!   wsq=0.
!   do m=1,2
!     indexh=memind(m)

!     y3(m)=p1*ww(ix ,jy ,n,indexh) &
!          +p2*ww(ixp,jy ,n,indexh) &
!          +p3*ww(ix ,jyp,n,indexh) &
!          +p4*ww(ixp,jyp,n,indexh)
!     wsl=wsl+ww(ix ,jy ,n,indexh)+ww(ixp,jy ,n,indexh) &
!          +ww(ix ,jyp,n,indexh)+ww(ixp,jyp,n,indexh)
!     wsq=wsq+ww(ix ,jy ,n,indexh)*ww(ix ,jy ,n,indexh)+ &
!          ww(ixp,jy ,n,indexh)*ww(ixp,jy ,n,indexh)+ &
!          ww(ix ,jyp,n,indexh)*ww(ix ,jyp,n,indexh)+ &
!          ww(ixp,jyp,n,indexh)*ww(ixp,jyp,n,indexh)
!   end do
!   wprof(n)=(y3(1)*dt2+y3(2)*dt1)*dtt
!   indzindicator(n)=.false.


!   ! Compute standard deviations
!   !****************************

!   xaux=wsq-wsl*wsl/8.
!   if (xaux.lt.eps) then
!     wsigprof(n)=0.
!   else
!     wsigprof(n)=sqrt(xaux/7.)
!   endif

!   ! Same for eta coordinates
!   usl=0.
!   vsl=0.
!   wsl=0.
!   usq=0.
!   vsq=0.
!   wsq=0.
!   do m=1,2
!     indexh=memind(m)

!     if (ngrid.lt.0) then
!       y1(m)=p1*uupoleta(ix ,jy ,n,indexh) &
!            +p2*uupoleta(ixp,jy ,n,indexh) &
!            +p3*uupoleta(ix ,jyp,n,indexh) &
!            +p4*uupoleta(ixp,jyp,n,indexh)
!       y2(m)=p1*vvpoleta(ix ,jy ,n,indexh) &
!            +p2*vvpoleta(ixp,jy ,n,indexh) &
!            +p3*vvpoleta(ix ,jyp,n,indexh) &
!            +p4*vvpoleta(ixp,jyp,n,indexh)
!         usl=usl+uupoleta(ix ,jy ,n,indexh)+uupoleta(ixp,jy ,n,indexh) &
!              +uupoleta(ix ,jyp,n,indexh)+uupoleta(ixp,jyp,n,indexh)
!         vsl=vsl+vvpoleta(ix ,jy ,n,indexh)+vvpoleta(ixp,jy ,n,indexh) &
!              +vvpoleta(ix ,jyp,n,indexh)+vvpoleta(ixp,jyp,n,indexh)

!         usq=usq+uupoleta(ix ,jy ,n,indexh)*uupoleta(ix ,jy ,n,indexh)+ &
!              uupoleta(ixp,jy ,n,indexh)*uupoleta(ixp,jy ,n,indexh)+ &
!              uupoleta(ix ,jyp,n,indexh)*uupoleta(ix ,jyp,n,indexh)+ &
!              uupoleta(ixp,jyp,n,indexh)*uupoleta(ixp,jyp,n,indexh)
!         vsq=vsq+vvpoleta(ix ,jy ,n,indexh)*vvpoleta(ix ,jy ,n,indexh)+ &
!              vvpoleta(ixp,jy ,n,indexh)*vvpoleta(ixp,jy ,n,indexh)+ &
!              vvpoleta(ix ,jyp,n,indexh)*vvpoleta(ix ,jyp,n,indexh)+ &
!              vvpoleta(ixp,jyp,n,indexh)*vvpoleta(ixp,jyp,n,indexh)
!     else
!       y1(m)=p1*uueta(ix ,jy ,n,indexh) &
!            +p2*uueta(ixp,jy ,n,indexh) &
!            +p3*uueta(ix ,jyp,n,indexh) &
!            +p4*uueta(ixp,jyp,n,indexh)
!       y2(m)=p1*vveta(ix ,jy ,n,indexh) &
!            +p2*vveta(ixp,jy ,n,indexh) &
!            +p3*vveta(ix ,jyp,n,indexh) &
!            +p4*vveta(ixp,jyp,n,indexh)
!       usl=usl+uueta(ix ,jy ,n,indexh)+uueta(ixp,jy ,n,indexh) &
!            +uueta(ix ,jyp,n,indexh)+uueta(ixp,jyp,n,indexh)
!       vsl=vsl+vveta(ix ,jy ,n,indexh)+vveta(ixp,jy ,n,indexh) &
!            +vveta(ix ,jyp,n,indexh)+vveta(ixp,jyp,n,indexh)

!       usq=usq+uueta(ix ,jy ,n,indexh)*uueta(ix ,jy ,n,indexh)+ &
!            uueta(ixp,jy ,n,indexh)*uueta(ixp,jy ,n,indexh)+ &
!            uueta(ix ,jyp,n,indexh)*uueta(ix ,jyp,n,indexh)+ &
!            uueta(ixp,jyp,n,indexh)*uueta(ixp,jyp,n,indexh)
!       vsq=vsq+vveta(ix ,jy ,n,indexh)*vveta(ix ,jy ,n,indexh)+ &
!            vveta(ixp,jy ,n,indexh)*vveta(ixp,jy ,n,indexh)+ &
!            vveta(ix ,jyp,n,indexh)*vveta(ix ,jyp,n,indexh)+ &
!            vveta(ixp,jyp,n,indexh)*vveta(ixp,jyp,n,indexh)
!     endif
!     y3(m)=p1*wweta(ix ,jy ,n,indexh) &
!          +p2*wweta(ixp,jy ,n,indexh) &
!          +p3*wweta(ix ,jyp,n,indexh) &
!          +p4*wweta(ixp,jyp,n,indexh)
!     rhograd1(m)=p1*drhodzeta(ix ,jy ,n,indexh) &
!          +p2*drhodzeta(ixp,jy ,n,indexh) &
!          +p3*drhodzeta(ix ,jyp,n,indexh) &
!          +p4*drhodzeta(ixp,jyp,n,indexh)
!     rho1(m)=p1*rhoeta(ix ,jy ,n,indexh) &
!          +p2*rhoeta(ixp,jy ,n,indexh) &
!          +p3*rhoeta(ix ,jyp,n,indexh) &
!          +p4*rhoeta(ixp,jyp,n,indexh)

!     wsl=wsl+wweta(ix ,jy ,n,indexh)+wweta(ixp,jy ,n,indexh) &
!          +wweta(ix ,jyp,n,indexh)+wweta(ixp,jyp,n,indexh)
!     wsq=wsq+wweta(ix ,jy ,n,indexh)*wweta(ix ,jy ,n,indexh)+ &
!          wweta(ixp,jy ,n,indexh)*wweta(ixp,jy ,n,indexh)+ &
!          wweta(ix ,jyp,n,indexh)*wweta(ix ,jyp,n,indexh)+ &
!          wweta(ixp,jyp,n,indexh)*wweta(ixp,jyp,n,indexh)
!   end do

!   ! Convert w from Pa/s to eta/s, following FLEXTRA
!   !************************************************

!   uprof(n)=(y1(1)*dt2+y1(2)*dt1)*dtt
!   vprof(n)=(y2(1)*dt2+y2(2)*dt1)*dtt
!   wprofeta(n)=(y3(1)*dt2+y3(2)*dt1)*dtt
!   rhoprof(n)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
!   rhogradprof(n)=(rhograd1(1)*dt2+rhograd1(2)*dt1)*dtt


!   ! Compute standard deviations
!   !****************************

!   xaux=usq-usl*usl/8.
!   if (xaux.lt.eps) then
!     usigprof(n)=0.
!   else
!     usigprof(n)=sqrt(xaux/7.)
!   endif

!   xaux=vsq-vsl*vsl/8.
!   if (xaux.lt.eps) then
!     vsigprof(n)=0.
!   else
!     vsigprof(n)=sqrt(xaux/7.)
!   endif

!   xaux=wsq-wsl*wsl/8.
!   if (xaux.lt.eps) then
!     wsigprofeta(n)=0.
!   else
!     wsigprofeta(n)=sqrt(xaux/7.)
!   endif
! end subroutine interpol_misslev

! subroutine interpol_wind(itime,xt,yt,zt,zteta,pp)
!   !                           i   i  i  i
!   !*****************************************************************************
!   !                                                                            *
!   !  This subroutine interpolates the wind data to current trajectory position.*
!   !                                                                            *
!   !    Author: A. Stohl                                                        *
!   !                                                                            *
!   !    16 December 1997                                                        *
!   !                                                                            *
!   !  Revision March 2005 by AST : all output variables in common block cal-    *
!   !                               culation of standard deviation done in this  *
!   !                               routine rather than subroutine call in order *
!   !                               to save computation time                     *
!   !                                                                            *
!   !*****************************************************************************
!   !                                                                            *
!   ! Variables:                                                                 *
!   ! u,v,w              wind components                                         *
!   ! itime [s]          current temporal position                               *
!   ! memtime(3) [s]     times of the wind fields in memory                      *
!   ! xt,yt,zt           coordinates position for which wind data shall be       *
!   !                    calculated                                              *
!   !                                                                            *
!   ! Constants:                                                                 *
!   !                                                                            *
!   !*****************************************************************************

!   use par_mod
!   use com_mod
!   use interpol_mod

!   implicit none

!   integer :: itime,pp
!   real :: xt,yt,zt
!   real :: zteta

!   ! Auxiliary variables needed for interpolation
!   real :: dz1,dz2,dz,psint(2)
!   real :: u1(2),v1(2),w1(2),dpdeta1(2),uh(2),vh(2),wh(2)
!   real :: usl,vsl,wsl,usq,vsq,wsq,xaux,dpdeta(2),psint_t,dpdetatemp
!   integer :: i,m,n,indexh,indzh
!   real,parameter :: eps=1.0e-30


!   !********************************************
!   ! Multilinear interpolation in time and space
!   !********************************************

!   ! Determine the lower left corner and its distance to the current position
!   !*************************************************************************

!   ddx=xt-real(ix)
!   ddy=yt-real(jy)
!   rddx=1.-ddx
!   rddy=1.-ddy
!   p1=rddx*rddy
!   p2=ddx*rddy
!   p3=rddx*ddy
!   p4=ddx*ddy

!   ! Calculate variables for time interpolation
!   !*******************************************

!   dt1=real(itime-memtime(1))
!   dt2=real(memtime(2)-itime)
!   dtt=1./(dt1+dt2)

!   ! Determine the level below the current position for u,v
!   !*******************************************************
!   indz=nz-1
!   do i=2,nz
!     if (height(i).gt.zt) then
!       indz=i-1
!       exit
!     endif
!   end do

!   ! Vertical distance to the level below and above current position
!   !****************************************************************

!   dz=1./(height(indz+1)-height(indz))
!   dz1=(zt-height(indz))*dz
!   dz2=(height(indz+1)-zt)*dz

!   !**********************************************************************
!   ! 1.) Bilinear horizontal interpolation
!   ! This has to be done separately for 6 fields (Temporal(2)*Vertical(3))
!   !**********************************************************************

!   ! Loop over 2 time steps and 2 levels
!   !************************************

!   wsl=0.
!   wsq=0.
!   do m=1,2
!     indexh=memind(m)
!     do n=1,2
!       indzh=indz+n-1
!       w1(n)=p1*ww(ix ,jy ,indzh,indexh) &
!            +p2*ww(ixp,jy ,indzh,indexh) &
!            +p3*ww(ix ,jyp,indzh,indexh) &
!            +p4*ww(ixp,jyp,indzh,indexh)
!       wsl=wsl+ww(ix ,jy ,indzh,indexh)+ww(ixp,jy ,indzh,indexh) &
!            +ww(ix ,jyp,indzh,indexh)+ww(ixp,jyp,indzh,indexh)
!       wsq=wsq+ww(ix ,jy ,indzh,indexh)*ww(ix ,jy ,indzh,indexh)+ &
!            ww(ixp,jy ,indzh,indexh)*ww(ixp,jy ,indzh,indexh)+ &
!            ww(ix ,jyp,indzh,indexh)*ww(ix ,jyp,indzh,indexh)+ &
!            ww(ixp,jyp,indzh,indexh)*ww(ixp,jyp,indzh,indexh)
!     end do


!   !**********************************
!   ! 2.) Linear vertical interpolation
!   !**********************************

!     wh(m)=dz2*w1(1)+dz1*w1(2)
!   end do


!   !************************************
!   ! 3.) Temporal interpolation (linear)
!   !************************************

!   w=(wh(1)*dt2+wh(2)*dt1)*dtt

!   ! Compute standard deviations
!   !****************************

!   xaux=wsq-wsl*wsl/16.
!   if (xaux.lt.eps) then
!     wsig=0.
!   else
!     wsig=sqrt(xaux/15.)
!   endif

!   ! Same for eta coordinates
!   !*************************

!   induv=nz-1
!   indpuv=nz
!   do i=2,nz
!     if (uvheight(i).lt.zteta) then
!       induv=i-1
!       indpuv=i
!       exit
!     endif
!   end do

!   dz=1./(uvheight(induv+1)-uvheight(induv))
!   dz1=(zteta-uvheight(induv))*dz
!   dz2=(uvheight(induv+1)-zteta)*dz
!   ! if (pp.eq.1) write(*,*) 'uv: ', zteta,induv,uvheight(induv),uvheight(indpuv),dz1,dz2

!   usl=0.
!   vsl=0.
!   usq=0.
!   vsq=0.
!   do m=1,2
!     indexh=memind(m)
!     do n=1,2
!       indzh=induv+n-1

!       if (ngrid.lt.0) then
!         u1(n)=p1*uupoleta(ix ,jy ,indzh,indexh) &
!              +p2*uupoleta(ixp,jy ,indzh,indexh) &
!              +p3*uupoleta(ix ,jyp,indzh,indexh) &
!              +p4*uupoleta(ixp,jyp,indzh,indexh)
!         v1(n)=p1*vvpoleta(ix ,jy ,indzh,indexh) &
!              +p2*vvpoleta(ixp,jy ,indzh,indexh) &
!              +p3*vvpoleta(ix ,jyp,indzh,indexh) &
!              +p4*vvpoleta(ixp,jyp,indzh,indexh)
!         usl=usl+uupoleta(ix ,jy ,indzh,indexh)+ &
!              uupoleta(ixp,jy ,indzh,indexh) &
!              +uupoleta(ix ,jyp,indzh,indexh)+uupoleta(ixp,jyp,indzh,indexh)
!         vsl=vsl+vvpoleta(ix ,jy ,indzh,indexh)+ &
!              vvpoleta(ixp,jy ,indzh,indexh) &
!              +vvpoleta(ix ,jyp,indzh,indexh)+vvpoleta(ixp,jyp,indzh,indexh)

!         usq=usq+uupoleta(ix ,jy ,indzh,indexh)* &
!              uupoleta(ix ,jy ,indzh,indexh)+ &
!              uupoleta(ixp,jy ,indzh,indexh)*uupoleta(ixp,jy ,indzh,indexh)+ &
!              uupoleta(ix ,jyp,indzh,indexh)*uupoleta(ix ,jyp,indzh,indexh)+ &
!              uupoleta(ixp,jyp,indzh,indexh)*uupoleta(ixp,jyp,indzh,indexh)
!         vsq=vsq+vvpoleta(ix ,jy ,indzh,indexh)* &
!              vvpoleta(ix ,jy ,indzh,indexh)+ &
!              vvpoleta(ixp,jy ,indzh,indexh)*vvpoleta(ixp,jy ,indzh,indexh)+ &
!              vvpoleta(ix ,jyp,indzh,indexh)*vvpoleta(ix ,jyp,indzh,indexh)+ &
!              vvpoleta(ixp,jyp,indzh,indexh)*vvpoleta(ixp,jyp,indzh,indexh)
!       else
!         u1(n)=p1*uueta(ix ,jy ,indzh,indexh) &
!              +p2*uueta(ixp,jy ,indzh,indexh) &
!              +p3*uueta(ix ,jyp,indzh,indexh) &
!              +p4*uueta(ixp,jyp,indzh,indexh)
!         v1(n)=p1*vveta(ix ,jy ,indzh,indexh) &
!              +p2*vveta(ixp,jy ,indzh,indexh) &
!              +p3*vveta(ix ,jyp,indzh,indexh) &
!              +p4*vveta(ixp,jyp,indzh,indexh)
!         usl=usl+uueta(ix ,jy ,indzh,indexh)+uueta(ixp,jy ,indzh,indexh) &
!              +uueta(ix ,jyp,indzh,indexh)+uueta(ixp,jyp,indzh,indexh)
!         vsl=vsl+vveta(ix ,jy ,indzh,indexh)+vveta(ixp,jy ,indzh,indexh) &
!              +vveta(ix ,jyp,indzh,indexh)+vveta(ixp,jyp,indzh,indexh)

!         usq=usq+uueta(ix ,jy ,indzh,indexh)*uueta(ix ,jy ,indzh,indexh)+ &
!              uueta(ixp,jy ,indzh,indexh)*uueta(ixp,jy ,indzh,indexh)+ &
!              uueta(ix ,jyp,indzh,indexh)*uueta(ix ,jyp,indzh,indexh)+ &
!              uueta(ixp,jyp,indzh,indexh)*uueta(ixp,jyp,indzh,indexh)
!         vsq=vsq+vveta(ix ,jy ,indzh,indexh)*vveta(ix ,jy ,indzh,indexh)+ &
!              vveta(ixp,jy ,indzh,indexh)*vveta(ixp,jy ,indzh,indexh)+ &
!              vveta(ix ,jyp,indzh,indexh)*vveta(ix ,jyp,indzh,indexh)+ &
!              vveta(ixp,jyp,indzh,indexh)*vveta(ixp,jyp,indzh,indexh)
!       endif
!     end do

!   !**********************************
!   ! 2.) Linear vertical interpolation
!   !**********************************
!     uh(m)=dz2*u1(1)+dz1*u1(2)
!     vh(m)=dz2*v1(1)+dz1*v1(2)
!   end do

!   indzeta=nz-1
!   indzpeta=nz
!   do i=2,nz
!     if (wheight(i).lt.zteta) then
!       indzeta=i-1
!       indzpeta=i
!       exit
!     endif
!   end do

!   dz=1./(wheight(indzeta+1)-wheight(indzeta))
!   dz1=(zteta-wheight(indzeta))*dz
!   dz2=(wheight(indzeta+1)-zteta)*dz
!   ! if (pp.eq.1) write(*,*) 'w: ', zteta,indzeta,wheight(indzeta),wheight(indzpeta),dz1,dz2

!   wsl=0.
!   wsq=0.
!   do m=1,2
!     indexh=memind(m)
!     do n=1,2
!       indzh=indzeta+n-1
!       w1(n)=p1*wweta(ix ,jy ,indzh,indexh) &
!            +p2*wweta(ixp,jy ,indzh,indexh) &
!            +p3*wweta(ix ,jyp,indzh,indexh) &
!            +p4*wweta(ixp,jyp,indzh,indexh)

!       wsl=wsl+wweta(ix ,jy ,indzh,indexh)+wweta(ixp,jy ,indzh,indexh) &
!            +wweta(ix ,jyp,indzh,indexh)+wweta(ixp,jyp,indzh,indexh)
!       wsq=wsq+wweta(ix ,jy ,indzh,indexh)*wweta(ix ,jy ,indzh,indexh)+ &
!            wweta(ixp,jy ,indzh,indexh)*wweta(ixp,jy ,indzh,indexh)+ &
!            wweta(ix ,jyp,indzh,indexh)*wweta(ix ,jyp,indzh,indexh)+ &
!            wweta(ixp,jyp,indzh,indexh)*wweta(ixp,jyp,indzh,indexh)
!     end do

!   !**********************************
!   ! 2.) Linear vertical interpolation
!   !**********************************
!     wh(m)=dz2*w1(1)+dz1*w1(2)
!   end do

!   !************************************
!   ! 3.) Temporal interpolation (linear)
!   !************************************

!   u=(uh(1)*dt2+uh(2)*dt1)*dtt
!   v=(vh(1)*dt2+vh(2)*dt1)*dtt
!   weta=(wh(1)*dt2+wh(2)*dt1)*dtt

!   xaux=wsq-wsl*wsl/16.
!   if (xaux.lt.eps) then
!     wsigeta=0.
!   else
!     wsigeta=sqrt(xaux/15.)
!   endif

!   xaux=usq-usl*usl/16.
!   if (xaux.lt.eps) then
!     usig=0.
!   else
!     usig=sqrt(xaux/15.)
!   endif

!   xaux=vsq-vsl*vsl/16.
!   if (xaux.lt.eps) then
!     vsig=0.
!   else
!     vsig=sqrt(xaux/15.)
!   endif
! end subroutine interpol_wind

! subroutine interpol_wind_short(itime,xt,yt,zt,zteta)
!   !                                 i   i  i  i
!   !*****************************************************************************
!   !                                                                            *
!   !  This subroutine interpolates the wind data to current trajectory position.*
!   !                                                                            *
!   !    Author: A. Stohl                                                        *
!   !                                                                            *
!   !    16 December 1997                                                        *
!   !                                                                            *
!   !  Revision March 2005 by AST : all output variables in common block         *
!   !                                                                            *
!   !*****************************************************************************
!   !                                                                            *
!   ! Variables:                                                                 *
!   ! u,v,w              wind components                                         *
!   ! itime [s]          current temporal position                               *
!   ! memtime(3) [s]     times of the wind fields in memory                      *
!   ! xt,yt,zt           coordinates position for which wind data shall be       *
!   !                    calculated                                              *
!   !                                                                            *
!   ! Constants:                                                                 *
!   !                                                                            *
!   !*****************************************************************************

!   use par_mod
!   use com_mod
!   use interpol_mod

!   implicit none

!   integer, intent(in) :: itime
!   real, intent(in) :: xt,yt,zt
!   real, intent(in) :: zteta

!   ! Auxiliary variables needed for interpolation
!   real :: dz1,dz2,dz
!   real :: u1(2),v1(2),w1(2),uh(2),vh(2),wh(2),dpdeta1(2)
!   integer :: i,m,n,indexh,indzh,psint(2),psint_t,dpdeta


!   !********************************************
!   ! Multilinear interpolation in time and space
!   !********************************************

!   ddx=xt-real(ix)
!   ddy=yt-real(jy)
!   rddx=1.-ddx
!   rddy=1.-ddy
!   p1=rddx*rddy
!   p2=ddx*rddy
!   p3=rddx*ddy
!   p4=ddx*ddy

!   ! Calculate variables for time interpolation
!   !*******************************************

!   dt1=real(itime-memtime(1))
!   dt2=real(memtime(2)-itime)
!   dtt=1./(dt1+dt2)

!   ! Determine the level below the current position for u,v
!   !*******************************************************

!   do i=2,nz
!     if (height(i).gt.zt) then
!       indz=i-1
!       exit
!     endif
!   end do

!   ! Vertical distance to the level below and above current position
!   !****************************************************************

!   dz=1./(height(indz+1)-height(indz))
!   dz1=(zt-height(indz))*dz
!   dz2=(height(indz+1)-zt)*dz


!   !**********************************************************************
!   ! 1.) Bilinear horizontal interpolation
!   ! This has to be done separately for 6 fields (Temporal(2)*Vertical(3))
!   !**********************************************************************

!   ! Loop over 2 time steps and 2 levels
!   !************************************
!   do m=1,2
!     indexh=memind(m)
!     do n=1,2
!       indzh=indz+n-1
!       w1(n)=p1*ww(ix ,jy ,indzh,indexh) &
!            +p2*ww(ixp,jy ,indzh,indexh) &
!            +p3*ww(ix ,jyp,indzh,indexh) &
!            +p4*ww(ixp,jyp,indzh,indexh)
!     end do


!   !**********************************
!   ! 2.) Linear vertical interpolation
!   !**********************************

!     wh(m)=dz2*w1(1)+dz1*w1(2)
!   end do



!   !************************************
!   ! 3.) Temporal interpolation (linear)
!   !************************************
!   w=(wh(1)*dt2+wh(2)*dt1)*dtt


!   ! Same for eta coordinates
!   !*************************

!   induv=nz-1
!   dz1=1.
!   dz2=0.
!   dz=1.
!   do i=2,nz
!     if (uvheight(i).lt.zteta) then
!       induv=i-1
!       dz=1./(uvheight(induv+1)-uvheight(induv))
!       dz1=(zteta-uvheight(induv))*dz
!       dz2=(uvheight(induv+1)-zteta)*dz
!       exit
!     endif
!   end do

!   do m=1,2
!     indexh=memind(m)
!     do n=1,2
!       indzh=induv+n-1

!       if (ngrid.lt.0) then
!         u1(n)=p1*uupoleta(ix ,jy ,indzh,indexh) &
!              +p2*uupoleta(ixp,jy ,indzh,indexh) &
!              +p3*uupoleta(ix ,jyp,indzh,indexh) &
!              +p4*uupoleta(ixp,jyp,indzh,indexh)
!         v1(n)=p1*vvpoleta(ix ,jy ,indzh,indexh) &
!              +p2*vvpoleta(ixp,jy ,indzh,indexh) &
!              +p3*vvpoleta(ix ,jyp,indzh,indexh) &
!              +p4*vvpoleta(ixp,jyp,indzh,indexh)
!       else
!         u1(n)=p1*uueta(ix ,jy ,indzh,indexh) &
!              +p2*uueta(ixp,jy ,indzh,indexh) &
!              +p3*uueta(ix ,jyp,indzh,indexh) &
!              +p4*uueta(ixp,jyp,indzh,indexh)
!         v1(n)=p1*vveta(ix ,jy ,indzh,indexh) &
!              +p2*vveta(ixp,jy ,indzh,indexh) &
!              +p3*vveta(ix ,jyp,indzh,indexh) &
!              +p4*vveta(ixp,jyp,indzh,indexh)
!       endif
!     end do

!   !**********************************
!   ! 2.) Linear vertical interpolation
!   !**********************************

!     uh(m)=dz2*u1(1)+dz1*u1(2)
!     vh(m)=dz2*v1(1)+dz1*v1(2)
!   end do

!   indzeta=nz-1
!   dz1=1.
!   dz2=0.
!   dz=1.
!   do i=2,nz
!     if (wheight(i).lt.zteta) then
!       indzeta=i-1
!       dz=1./(wheight(indzeta+1)-wheight(indzeta))
!       dz1=(zteta-wheight(indzeta))*dz
!       dz2=(wheight(indzeta+1)-zteta)*dz
!       exit
!     endif
!   end do

!   do m=1,2
!     indexh=memind(m)
!     do n=1,2
!       indzh=indzeta+n-1
!       w1(n)=p1*wweta(ix ,jy ,indzh,indexh) &
!            +p2*wweta(ixp,jy ,indzh,indexh) &
!            +p3*wweta(ix ,jyp,indzh,indexh) &
!            +p4*wweta(ixp,jyp,indzh,indexh)
!     end do
!     wh(m)=dz2*w1(1)+dz1*w1(2)
!   end do

!   !************************************
!   ! 3.) Temporal interpolation (linear)
!   !************************************

!   u=(uh(1)*dt2+uh(2)*dt1)*dtt
!   v=(vh(1)*dt2+vh(2)*dt1)*dtt
!   weta=(wh(1)*dt2+wh(2)*dt1)*dtt
! end subroutine interpol_wind_short

! subroutine interpol_rain(yy1,yy2,yy3,nxmax,nymax,nzmax,nx, &
!      ny,iwftouse,xt,yt,level,itime1,itime2,itime,yint1,yint2,yint3)
!   !                          i   i   i    i    i     i   i
!   !i    i    i  i    i     i      i      i     o     o     o
!   !****************************************************************************
!   !                                                                           *
!   !  Interpolation of meteorological fields on 2-d model layers.              *
!   !  In horizontal direction bilinear interpolation interpolation is used.    *
!   !  Temporally a linear interpolation is used.                               *
!   !  Three fields are interpolated at the same time.                          *
!   !                                                                           *
!   !  This is a special version of levlininterpol to save CPU time.            *
!   !                                                                           *
!   !  1 first time                                                             *
!   !  2 second time                                                            *
!   !                                                                           *
!   !                                                                           *
!   !     Author: A. Stohl                                                      *
!   !                                                                           *
!   !     30 August 1996                                                        *
!   !                                                                           *
!   !****************************************************************************
!   !                                                                           *
!   ! Variables:                                                                *
!   !                                                                           *
!   ! dt1,dt2              time differences between fields and current position *
!   ! dz1,dz2              z distance between levels and current position       *
!   ! height(nzmax)        heights of the model levels                          *
!   ! indexh               help variable                                        *
!   ! indz                 the level closest to the current trajectory position *
!   ! indzh                help variable                                        *
!   ! itime                current time                                         *
!   ! itime1               time of the first wind field                         *
!   ! itime2               time of the second wind field                        *
!   ! ix,jy                x,y coordinates of lower left subgrid point          *
!   ! level                level at which interpolation shall be done           *
!   ! iwftouse             points to the place of the wind field                *
!   ! nx,ny                actual field dimensions in x,y and z direction       *
!   ! nxmax,nymax,nzmax    maximum field dimensions in x,y and z direction      *
!   ! xt                   current x coordinate                                 *
!   ! yint                 the final interpolated value                         *
!   ! yt                   current y coordinate                                 *
!   ! yy(0:nxmax,0:nymax,nzmax,3) meteorological field used for interpolation   *
!   ! zt                   current z coordinate                                 *
!   !                                                                           *
!   !****************************************************************************
!   use par_mod, only: numwfmem

!   implicit none

!   integer :: nx,ny,nxmax,nymax,nzmax,memind(numwfmem),m,ix,jy,ixp,jyp
!   integer :: itime,itime1,itime2,level,indexh
!   real :: yy1(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
!   real :: yy2(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
!   real :: yy3(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
!   real :: ddx,ddy,rddx,rddy,dt1,dt2,dt,y1(2),y2(2),y3(2)
!   real :: xt,yt,yint1,yint2,yint3,p1,p2,p3,p4
!   integer :: iwftouse



!   ! If point at border of grid -> small displacement into grid
!   !***********************************************************

!   if (xt.ge.real(nx-1)) xt=real(nx-1)-0.00001
!   if (yt.ge.real(ny-1)) yt=real(ny-1)-0.00001



!   !**********************************************************************
!   ! 1.) Bilinear horizontal interpolation
!   ! This has to be done separately for 2 fields (Temporal)
!   !*******************************************************

!   ! Determine the lower left corner and its distance to the current position
!   !*************************************************************************

!   ix=int(xt)
!   jy=int(yt)
!   ixp=ix+1
!   jyp=jy+1
!   ddx=xt-real(ix)
!   ddy=yt-real(jy)
!   rddx=1.-ddx
!   rddy=1.-ddy
!   p1=rddx*rddy
!   p2=ddx*rddy
!   p3=rddx*ddy
!   p4=ddx*ddy


!   ! Loop over 2 time steps
!   !***********************

!   !  do m=1,2
!   indexh=iwftouse

!   y1(1)=p1*yy1(ix ,jy ,level,indexh) &
!      + p2*yy1(ixp,jy ,level,indexh) &
!      + p3*yy1(ix ,jyp,level,indexh) &
!      + p4*yy1(ixp,jyp,level,indexh)
!   y2(1)=p1*yy2(ix ,jy ,level,indexh) &
!      + p2*yy2(ixp,jy ,level,indexh) &
!      + p3*yy2(ix ,jyp,level,indexh) &
!      + p4*yy2(ixp,jyp,level,indexh)
!   y3(1)=p1*yy3(ix ,jy ,level,indexh) &
!      + p2*yy3(ixp,jy ,level,indexh) &
!      + p3*yy3(ix ,jyp,level,indexh) &
!      + p4*yy3(ixp,jyp,level,indexh)
!   !  end do


!   !************************************
!   ! 2.) Temporal interpolation (linear) - skip to be consistent with clouds
!   !************************************

!   !  dt1=real(itime-itime1)
!   !  dt2=real(itime2-itime)
!   !  dt=dt1+dt2

!   !  yint1=(y1(1)*dt2+y1(2)*dt1)/dt
!   !  yint2=(y2(1)*dt2+y2(2)*dt1)/dt
!   !  yint3=(y3(1)*dt2+y3(2)*dt1)/dt

!    yint1=y1(1)
!    yint2=y2(1)
!    yint3=y3(1)
! end subroutine interpol_rain

! subroutine interpol_vdep(level,vdepo)
!   !                           i     o
!   !****************************************************************************
!   !                                                                           *
!   !  Interpolation of the deposition velocity on 2-d model layer.             *
!   !  In horizontal direction bilinear interpolation interpolation is used.    *
!   !  Temporally a linear interpolation is used.                               *
!   !                                                                           *
!   !  1 first time                                                             *
!   !  2 second time                                                            *
!   !                                                                           *
!   !                                                                           *
!   !     Author: A. Stohl                                                      *
!   !                                                                           *
!   !     30 May 1994                                                           *
!   !                                                                           *
!   !****************************************************************************
!   !                                                                           *
!   ! Variables:                                                                *
!   !                                                                           *
!   ! level                number of species for which interpolation is done    *
!   !                                                                           *
!   !****************************************************************************

!   use par_mod
!   use com_mod
!   use interpol_mod

!   implicit none

!   integer :: level,indexh,m
!   real :: y(2),vdepo

!   ! a) Bilinear horizontal interpolation
!   do m=1,2
!     indexh=memind(m)

!     y(m)=p1*vdep(ix ,jy ,level,indexh) &
!          +p2*vdep(ixp,jy ,level,indexh) &
!          +p3*vdep(ix ,jyp,level,indexh) &
!          +p4*vdep(ixp,jyp,level,indexh)
!   end do



!   ! b) Temporal interpolation

!   vdepo=(y(1)*dt2+y(2)*dt1)*dtt

!   depoindicator(level)=.false.
! end subroutine interpol_vdep

! subroutine interpol_all_nests(itime,xt,yt,zt)
!   !                                i   i  i  i
!   !*****************************************************************************
!   !                                                                            *
!   !  This subroutine interpolates everything that is needed for calculating the*
!   !  dispersion.                                                               *
!   !  Version for interpolating nested grids.                                   *
!   !                                                                            *
!   !    Author: A. Stohl                                                        *
!   !                                                                            *
!   !    9 February 1999                                                         *
!   !    16 December 1997                                                        *
!   !                                                                            *
!   !  Revision March 2005 by AST : all output variables in common block cal-    *
!   !                               culation of standard deviation done in this  *
!   !                               routine rather than subroutine call in order *
!   !                               to save computation time                     *
!   !                                                                            *
!   !*****************************************************************************
!   !                                                                            *
!   ! Variables:                                                                 *
!   ! itime [s]          current temporal position                               *
!   ! memtime(3) [s]     times of the wind fields in memory                      *
!   ! xt,yt,zt           coordinates position for which wind data shall be       *
!   !                    calculated                                              *
!   !                                                                            *
!   ! Constants:                                                                 *
!   !                                                                            *
!   !*****************************************************************************

!   use par_mod
!   use com_mod
!   use interpol_mod
!   use hanna_mod

!   implicit none

!   integer :: itime
!   real :: xt,yt,zt

!   ! Auxiliary variables needed for interpolation
!   real :: ust1(2),wst1(2),oli1(2),oliaux
!   real :: y1(2),y2(2),y3(2),rho1(2),rhograd1(2)
!   real :: usl,vsl,wsl,usq,vsq,wsq,xaux
!   integer :: i,m,n,indexh
!   real,parameter :: eps=1.0e-30


!   !********************************************
!   ! Multilinear interpolation in time and space
!   !********************************************

!   ! Determine the lower left corner and its distance to the current position
!   !*************************************************************************

!   ddx=xt-real(ix)
!   ddy=yt-real(jy)
!   rddx=1.-ddx
!   rddy=1.-ddy
!   p1=rddx*rddy
!   p2=ddx*rddy
!   p3=rddx*ddy
!   p4=ddx*ddy

!   ! Calculate variables for time interpolation
!   !*******************************************

!   dt1=real(itime-memtime(1))
!   dt2=real(memtime(2)-itime)
!   dtt=1./(dt1+dt2)


!   !*****************************************
!   ! 1. Interpolate u*, w* and Obukhov length
!   !*****************************************

!   ! a) Bilinear horizontal interpolation

!   do m=1,2
!     indexh=memind(m)

!     ust1(m)=p1*ustarn(ix ,jy ,1,indexh,ngrid) &
!          + p2*ustarn(ixp,jy ,1,indexh,ngrid) &
!          + p3*ustarn(ix ,jyp,1,indexh,ngrid) &
!          + p4*ustarn(ixp,jyp,1,indexh,ngrid)
!     wst1(m)=p1*wstarn(ix ,jy ,1,indexh,ngrid) &
!          + p2*wstarn(ixp,jy ,1,indexh,ngrid) &
!          + p3*wstarn(ix ,jyp,1,indexh,ngrid) &
!          + p4*wstarn(ixp,jyp,1,indexh,ngrid)
!     oli1(m)=p1*olin(ix ,jy ,1,indexh,ngrid) &
!          + p2*olin(ixp,jy ,1,indexh,ngrid) &
!          + p3*olin(ix ,jyp,1,indexh,ngrid) &
!          + p4*olin(ixp,jyp,1,indexh,ngrid)
!   end do

!   ! b) Temporal interpolation

!   ust=(ust1(1)*dt2+ust1(2)*dt1)*dtt
!   wst=(wst1(1)*dt2+wst1(2)*dt1)*dtt
!   oliaux=(oli1(1)*dt2+oli1(2)*dt1)*dtt

!   if (oliaux.ne.0.) then
!     ol=1./oliaux
!   else
!     ol=99999.
!   endif


!   !*****************************************************
!   ! 2. Interpolate vertical profiles of u,v,w,rho,drhodz
!   !*****************************************************


!   ! Determine the level below the current position
!   !***********************************************

!   do i=2,nz
!     if (height(i).gt.zt) then
!       indz=i-1
!       indzp=i
!       exit
!     endif
!   end do

!   !**************************************
!   ! 1.) Bilinear horizontal interpolation
!   ! 2.) Temporal interpolation (linear)
!   !**************************************

!   ! Loop over 2 time steps and indz levels
!   !***************************************

!   do n=indz,indz+1
!     usl=0.
!     vsl=0.
!     wsl=0.
!     usq=0.
!     vsq=0.
!     wsq=0.
!     do m=1,2
!       indexh=memind(m)
!       y1(m)=p1*uun(ix ,jy ,n,indexh,ngrid) &
!            +p2*uun(ixp,jy ,n,indexh,ngrid) &
!            +p3*uun(ix ,jyp,n,indexh,ngrid) &
!            +p4*uun(ixp,jyp,n,indexh,ngrid)
!       y2(m)=p1*vvn(ix ,jy ,n,indexh,ngrid) &
!            +p2*vvn(ixp,jy ,n,indexh,ngrid) &
!            +p3*vvn(ix ,jyp,n,indexh,ngrid) &
!            +p4*vvn(ixp,jyp,n,indexh,ngrid)
!       y3(m)=p1*wwn(ix ,jy ,n,indexh,ngrid) &
!            +p2*wwn(ixp,jy ,n,indexh,ngrid) &
!            +p3*wwn(ix ,jyp,n,indexh,ngrid) &
!            +p4*wwn(ixp,jyp,n,indexh,ngrid)
!       rhograd1(m)=p1*drhodzn(ix ,jy ,n,indexh,ngrid) &
!            +p2*drhodzn(ixp,jy ,n,indexh,ngrid) &
!            +p3*drhodzn(ix ,jyp,n,indexh,ngrid) &
!            +p4*drhodzn(ixp,jyp,n,indexh,ngrid)
!       rho1(m)=p1*rhon(ix ,jy ,n,indexh,ngrid) &
!            +p2*rhon(ixp,jy ,n,indexh,ngrid) &
!            +p3*rhon(ix ,jyp,n,indexh,ngrid) &
!            +p4*rhon(ixp,jyp,n,indexh,ngrid)

!      usl=usl+uun(ix ,jy ,n,indexh,ngrid)+uun(ixp,jy ,n,indexh,ngrid) &
!           +uun(ix ,jyp,n,indexh,ngrid)+uun(ixp,jyp,n,indexh,ngrid)
!      vsl=vsl+vvn(ix ,jy ,n,indexh,ngrid)+vvn(ixp,jy ,n,indexh,ngrid) &
!           +vvn(ix ,jyp,n,indexh,ngrid)+vvn(ixp,jyp,n,indexh,ngrid)
!      wsl=wsl+wwn(ix ,jy ,n,indexh,ngrid)+wwn(ixp,jy ,n,indexh,ngrid) &
!           +wwn(ix ,jyp,n,indexh,ngrid)+wwn(ixp,jyp,n,indexh,ngrid)

!     usq=usq+uun(ix ,jy ,n,indexh,ngrid)*uun(ix ,jy ,n,indexh,ngrid)+ &
!          uun(ixp,jy ,n,indexh,ngrid)*uun(ixp,jy ,n,indexh,ngrid)+ &
!          uun(ix ,jyp,n,indexh,ngrid)*uun(ix ,jyp,n,indexh,ngrid)+ &
!          uun(ixp,jyp,n,indexh,ngrid)*uun(ixp,jyp,n,indexh,ngrid)
!     vsq=vsq+vvn(ix ,jy ,n,indexh,ngrid)*vvn(ix ,jy ,n,indexh,ngrid)+ &
!          vvn(ixp,jy ,n,indexh,ngrid)*vvn(ixp,jy ,n,indexh,ngrid)+ &
!          vvn(ix ,jyp,n,indexh,ngrid)*vvn(ix ,jyp,n,indexh,ngrid)+ &
!          vvn(ixp,jyp,n,indexh,ngrid)*vvn(ixp,jyp,n,indexh,ngrid)
!     wsq=wsq+wwn(ix ,jy ,n,indexh,ngrid)*wwn(ix ,jy ,n,indexh,ngrid)+ &
!          wwn(ixp,jy ,n,indexh,ngrid)*wwn(ixp,jy ,n,indexh,ngrid)+ &
!          wwn(ix ,jyp,n,indexh,ngrid)*wwn(ix ,jyp,n,indexh,ngrid)+ &
!          wwn(ixp,jyp,n,indexh,ngrid)*wwn(ixp,jyp,n,indexh,ngrid)
!     end do
!     uprof(n)=(y1(1)*dt2+y1(2)*dt1)*dtt
!     vprof(n)=(y2(1)*dt2+y2(2)*dt1)*dtt
!     wprof(n)=(y3(1)*dt2+y3(2)*dt1)*dtt
!     rhoprof(n)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
!     rhogradprof(n)=(rhograd1(1)*dt2+rhograd1(2)*dt1)*dtt
!     indzindicator(n)=.false.

!   ! Compute standard deviations
!   !****************************

!     xaux=usq-usl*usl/8.
!     if (xaux.lt.eps) then
!       usigprof(n)=0.
!     else
!       usigprof(n)=sqrt(xaux/7.)
!     endif

!     xaux=vsq-vsl*vsl/8.
!     if (xaux.lt.eps) then
!       vsigprof(n)=0.
!     else
!       vsigprof(n)=sqrt(xaux/7.)
!     endif


!     xaux=wsq-wsl*wsl/8.
!     if (xaux.lt.eps) then
!       wsigprof(n)=0.
!     else
!       wsigprof(n)=sqrt(xaux/7.)
!     endif

!   end do
! end subroutine interpol_all_nests

! subroutine interpol_misslev_nests(n)
!   !                                  i
!   !*****************************************************************************
!   !                                                                            *
!   !  This subroutine interpolates u,v,w, density and density gradients.        *
!   !                                                                            *
!   !    Author: A. Stohl                                                        *
!   !                                                                            *
!   !    16 December 1997                                                        *
!   !                                                                            *
!   !*****************************************************************************
!   !                                                                            *
!   ! Variables:                                                                 *
!   ! n                  level                                                   *
!   !                                                                            *
!   ! Constants:                                                                 *
!   !                                                                            *
!   !*****************************************************************************

!   use par_mod
!   use com_mod
!   use interpol_mod
!   use hanna_mod

!   implicit none

!   ! Auxiliary variables needed for interpolation
!   real :: y1(2),y2(2),y3(2),rho1(2),rhograd1(2)
!   real :: usl,vsl,wsl,usq,vsq,wsq,xaux
!   integer :: m,n,indexh
!   real,parameter :: eps=1.0e-30


!   !********************************************
!   ! Multilinear interpolation in time and space
!   !********************************************


!   !**************************************
!   ! 1.) Bilinear horizontal interpolation
!   ! 2.) Temporal interpolation (linear)
!   !**************************************

!   ! Loop over 2 time steps
!   !***********************

!   usl=0.
!   vsl=0.
!   wsl=0.
!   usq=0.
!   vsq=0.
!   wsq=0.
!   do m=1,2
!     indexh=memind(m)
!     y1(m)=p1*uun(ix ,jy ,n,indexh,ngrid) &
!          +p2*uun(ixp,jy ,n,indexh,ngrid) &
!          +p3*uun(ix ,jyp,n,indexh,ngrid) &
!          +p4*uun(ixp,jyp,n,indexh,ngrid)
!     y2(m)=p1*vvn(ix ,jy ,n,indexh,ngrid) &
!          +p2*vvn(ixp,jy ,n,indexh,ngrid) &
!          +p3*vvn(ix ,jyp,n,indexh,ngrid) &
!          +p4*vvn(ixp,jyp,n,indexh,ngrid)
!     y3(m)=p1*wwn(ix ,jy ,n,indexh,ngrid) &
!          +p2*wwn(ixp,jy ,n,indexh,ngrid) &
!          +p3*wwn(ix ,jyp,n,indexh,ngrid) &
!          +p4*wwn(ixp,jyp,n,indexh,ngrid)
!     rho1(m)=p1*rhon(ix ,jy ,n,indexh,ngrid) &
!          +p2*rhon(ixp,jy ,n,indexh,ngrid) &
!          +p3*rhon(ix ,jyp,n,indexh,ngrid) &
!          +p4*rhon(ixp,jyp,n,indexh,ngrid)
!     rhograd1(m)=p1*drhodzn(ix ,jy ,n,indexh,ngrid) &
!          +p2*drhodzn(ixp,jy ,n,indexh,ngrid) &
!          +p3*drhodzn(ix ,jyp,n,indexh,ngrid) &
!          +p4*drhodzn(ixp,jyp,n,indexh,ngrid)

!      usl=usl+uun(ix ,jy ,n,indexh,ngrid)+uun(ixp,jy ,n,indexh,ngrid) &
!           +uun(ix ,jyp,n,indexh,ngrid)+uun(ixp,jyp,n,indexh,ngrid)
!      vsl=vsl+vvn(ix ,jy ,n,indexh,ngrid)+vvn(ixp,jy ,n,indexh,ngrid) &
!           +vvn(ix ,jyp,n,indexh,ngrid)+vvn(ixp,jyp,n,indexh,ngrid)
!      wsl=wsl+wwn(ix ,jy ,n,indexh,ngrid)+wwn(ixp,jy ,n,indexh,ngrid) &
!           +wwn(ix ,jyp,n,indexh,ngrid)+wwn(ixp,jyp,n,indexh,ngrid)

!     usq=usq+uun(ix ,jy ,n,indexh,ngrid)*uun(ix ,jy ,n,indexh,ngrid)+ &
!          uun(ixp,jy ,n,indexh,ngrid)*uun(ixp,jy ,n,indexh,ngrid)+ &
!          uun(ix ,jyp,n,indexh,ngrid)*uun(ix ,jyp,n,indexh,ngrid)+ &
!          uun(ixp,jyp,n,indexh,ngrid)*uun(ixp,jyp,n,indexh,ngrid)
!     vsq=vsq+vvn(ix ,jy ,n,indexh,ngrid)*vvn(ix ,jy ,n,indexh,ngrid)+ &
!          vvn(ixp,jy ,n,indexh,ngrid)*vvn(ixp,jy ,n,indexh,ngrid)+ &
!          vvn(ix ,jyp,n,indexh,ngrid)*vvn(ix ,jyp,n,indexh,ngrid)+ &
!          vvn(ixp,jyp,n,indexh,ngrid)*vvn(ixp,jyp,n,indexh,ngrid)
!     wsq=wsq+wwn(ix ,jy ,n,indexh,ngrid)*wwn(ix ,jy ,n,indexh,ngrid)+ &
!          wwn(ixp,jy ,n,indexh,ngrid)*wwn(ixp,jy ,n,indexh,ngrid)+ &
!          wwn(ix ,jyp,n,indexh,ngrid)*wwn(ix ,jyp,n,indexh,ngrid)+ &
!          wwn(ixp,jyp,n,indexh,ngrid)*wwn(ixp,jyp,n,indexh,ngrid)
!   end do
!   uprof(n)=(y1(1)*dt2+y1(2)*dt1)*dtt
!   vprof(n)=(y2(1)*dt2+y2(2)*dt1)*dtt
!   wprof(n)=(y3(1)*dt2+y3(2)*dt1)*dtt
!   rhoprof(n)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
!   rhogradprof(n)=(rhograd1(1)*dt2+rhograd1(2)*dt1)*dtt
!   indzindicator(n)=.false.

!   ! Compute standard deviations
!   !****************************

!   xaux=usq-usl*usl/8.
!   if (xaux.lt.eps) then
!     usigprof(n)=0.
!   else
!     usigprof(n)=sqrt(xaux/7.)
!   endif

!   xaux=vsq-vsl*vsl/8.
!   if (xaux.lt.eps) then
!     vsigprof(n)=0.
!   else
!     vsigprof(n)=sqrt(xaux/7.)
!   endif


!   xaux=wsq-wsl*wsl/8.
!   if (xaux.lt.eps) then
!     wsigprof(n)=0.
!   else
!     wsigprof(n)=sqrt(xaux/7.)
!   endif
! end subroutine interpol_misslev_nests

! subroutine interpol_wind_short_nests(itime,xt,yt,zt)
!   !                                       i   i  i  i
!   !*****************************************************************************
!   !                                                                            *
!   !  This subroutine interpolates the wind data to current trajectory position.*
!   !                                                                            *
!   !    Author: A. Stohl                                                        *
!   !                                                                            *
!   !    16 December 1997                                                        *
!   !                                                                            *
!   !  Revision March 2005 by AST : all output variables in common block         *
!   !                                                                            *
!   !*****************************************************************************
!   !                                                                            *
!   ! Variables:                                                                 *
!   ! u,v,w              wind components                                         *
!   ! itime [s]          current temporal position                               *
!   ! memtime(3) [s]     times of the wind fields in memory                      *
!   ! xt,yt,zt           coordinates position for which wind data shall be       *
!   !                    calculated                                              *
!   !                                                                            *
!   ! Constants:                                                                 *
!   !                                                                            *
!   !*****************************************************************************

!   use par_mod
!   use com_mod
!   use interpol_mod

!   implicit none

!   integer :: itime
!   real :: xt,yt,zt

!   ! Auxiliary variables needed for interpolation
!   real :: dz1,dz2,dz
!   real :: u1(2),v1(2),w1(2),uh(2),vh(2),wh(2)
!   integer :: i,m,n,indexh,indzh


!   !********************************************
!   ! Multilinear interpolation in time and space
!   !********************************************

!   ddx=xt-real(ix)
!   ddy=yt-real(jy)
!   rddx=1.-ddx
!   rddy=1.-ddy
!   p1=rddx*rddy
!   p2=ddx*rddy
!   p3=rddx*ddy
!   p4=ddx*ddy

!   ! Calculate variables for time interpolation
!   !*******************************************

!   dt1=real(itime-memtime(1))
!   dt2=real(memtime(2)-itime)
!   dtt=1./(dt1+dt2)

!   ! Determine the level below the current position for u,v
!   !*******************************************************

!   do i=2,nz
!     if (height(i).gt.zt) then
!       indz=i-1
!       exit
!     endif
!   end do

!   ! Vertical distance to the level below and above current position
!   !****************************************************************

!   dz=1./(height(indz+1)-height(indz))
!   dz1=(zt-height(indz))*dz
!   dz2=(height(indz+1)-zt)*dz


!   !**********************************************************************
!   ! 1.) Bilinear horizontal interpolation
!   ! This has to be done separately for 6 fields (Temporal(2)*Vertical(3))
!   !**********************************************************************

!   ! Loop over 2 time steps and 2 levels
!   !************************************

!   do m=1,2
!     indexh=memind(m)
!     do n=1,2
!       indzh=indz+n-1

!       u1(n)=p1*uun(ix ,jy ,indzh,indexh,ngrid) &
!            +p2*uun(ixp,jy ,indzh,indexh,ngrid) &
!            +p3*uun(ix ,jyp,indzh,indexh,ngrid) &
!            +p4*uun(ixp,jyp,indzh,indexh,ngrid)
!       v1(n)=p1*vvn(ix ,jy ,indzh,indexh,ngrid) &
!            +p2*vvn(ixp,jy ,indzh,indexh,ngrid) &
!            +p3*vvn(ix ,jyp,indzh,indexh,ngrid) &
!            +p4*vvn(ixp,jyp,indzh,indexh,ngrid)
!       w1(n)=p1*wwn(ix ,jy ,indzh,indexh,ngrid) &
!            +p2*wwn(ixp,jy ,indzh,indexh,ngrid) &
!            +p3*wwn(ix ,jyp,indzh,indexh,ngrid) &
!            +p4*wwn(ixp,jyp,indzh,indexh,ngrid)

!     end do


!   !**********************************
!   ! 2.) Linear vertical interpolation
!   !**********************************

!     uh(m)=dz2*u1(1)+dz1*u1(2)
!     vh(m)=dz2*v1(1)+dz1*v1(2)
!     wh(m)=dz2*w1(1)+dz1*w1(2)
!   end do


!   !************************************
!   ! 3.) Temporal interpolation (linear)
!   !************************************

!   u=(uh(1)*dt2+uh(2)*dt1)*dtt
!   v=(vh(1)*dt2+vh(2)*dt1)*dtt
!   w=(wh(1)*dt2+wh(2)*dt1)*dtt
! end subroutine interpol_wind_short_nests

! subroutine interpol_wind_nests(itime,xt,yt,zt)
!   !                                 i   i  i  i
!   !*****************************************************************************
!   !                                                                            *
!   !  This subroutine interpolates the wind data to current trajectory position.*
!   !                                                                            *
!   !    Author: A. Stohl                                                        *
!   !                                                                            *
!   !    16 December 1997                                                        *
!   !    16 December 1997                                                        *
!   !                                                                            *
!   !  Revision March 2005 by AST : all output variables in common block cal-    *
!   !                               culation of standard deviation done in this  *
!   !                               routine rather than subroutine call in order *
!   !                               to save computation time                     *
!   !                                                                            *
!   !*****************************************************************************
!   !                                                                            *
!   ! Variables:                                                                 *
!   ! u,v,w              wind components                                         *
!   ! itime [s]          current temporal position                               *
!   ! memtime(3) [s]     times of the wind fields in memory                      *
!   ! xt,yt,zt           coordinates position for which wind data shall be       *
!   !                    calculated                                              *
!   !                                                                            *
!   ! Constants:                                                                 *
!   !                                                                            *
!   !*****************************************************************************

!   use par_mod
!   use com_mod
!   use interpol_mod

!   implicit none

!   integer :: itime
!   real :: xt,yt,zt

!   ! Auxiliary variables needed for interpolation
!   real :: dz1,dz2,dz
!   real :: u1(2),v1(2),w1(2),uh(2),vh(2),wh(2)
!   real :: usl,vsl,wsl,usq,vsq,wsq,xaux
!   integer :: i,m,n,indexh,indzh
!   real,parameter :: eps=1.0e-30


!   !********************************************
!   ! Multilinear interpolation in time and space
!   !********************************************

!   ! Determine the lower left corner and its distance to the current position
!   !*************************************************************************

!   ddx=xt-real(ix)
!   ddy=yt-real(jy)
!   rddx=1.-ddx
!   rddy=1.-ddy
!   p1=rddx*rddy
!   p2=ddx*rddy
!   p3=rddx*ddy
!   p4=ddx*ddy

!   ! Calculate variables for time interpolation
!   !*******************************************

!   dt1=real(itime-memtime(1))
!   dt2=real(memtime(2)-itime)
!   dtt=1./(dt1+dt2)

!   ! Determine the level below the current position for u,v
!   !*******************************************************

!   do i=2,nz
!     if (height(i).gt.zt) then
!       indz=i-1
!       exit
!     endif
!   end do

!   ! Vertical distance to the level below and above current position
!   !****************************************************************

!   dz=1./(height(indz+1)-height(indz))
!   dz1=(zt-height(indz))*dz
!   dz2=(height(indz+1)-zt)*dz


!   !**********************************************************************
!   ! 1.) Bilinear horizontal interpolation
!   ! This has to be done separately for 6 fields (Temporal(2)*Vertical(3))
!   !**********************************************************************

!   ! Loop over 2 time steps and 2 levels
!   !************************************

!   usl=0.
!   vsl=0.
!   wsl=0.
!   usq=0.
!   vsq=0.
!   wsq=0.
!   do m=1,2
!     indexh=memind(m)
!     do n=1,2
!       indzh=indz+n-1

!       u1(n)=p1*uun(ix ,jy ,indzh,indexh,ngrid) &
!            +p2*uun(ixp,jy ,indzh,indexh,ngrid) &
!            +p3*uun(ix ,jyp,indzh,indexh,ngrid) &
!            +p4*uun(ixp,jyp,indzh,indexh,ngrid)
!       v1(n)=p1*vvn(ix ,jy ,indzh,indexh,ngrid) &
!            +p2*vvn(ixp,jy ,indzh,indexh,ngrid) &
!            +p3*vvn(ix ,jyp,indzh,indexh,ngrid) &
!            +p4*vvn(ixp,jyp,indzh,indexh,ngrid)
!       w1(n)=p1*wwn(ix ,jy ,indzh,indexh,ngrid) &
!            +p2*wwn(ixp,jy ,indzh,indexh,ngrid) &
!            +p3*wwn(ix ,jyp,indzh,indexh,ngrid) &
!            +p4*wwn(ixp,jyp,indzh,indexh,ngrid)

!       usl=usl+uun(ix ,jy ,indzh,indexh,ngrid)+ &
!            uun(ixp,jy ,indzh,indexh,ngrid) &
!            +uun(ix ,jyp,indzh,indexh,ngrid)+ &
!            uun(ixp,jyp,indzh,indexh,ngrid)
!       vsl=vsl+vvn(ix ,jy ,indzh,indexh,ngrid)+ &
!            vvn(ixp,jy ,indzh,indexh,ngrid) &
!            +vvn(ix ,jyp,indzh,indexh,ngrid)+ &
!            vvn(ixp,jyp,indzh,indexh,ngrid)
!       wsl=wsl+wwn(ix ,jy ,indzh,indexh,ngrid)+ &
!            wwn(ixp,jy ,indzh,indexh,ngrid) &
!            +wwn(ix ,jyp,indzh,indexh,ngrid)+ &
!            wwn(ixp,jyp,indzh,indexh,ngrid)

!       usq=usq+uun(ix ,jy ,indzh,indexh,ngrid)* &
!            uun(ix ,jy ,indzh,indexh,ngrid)+ &
!            uun(ixp,jy ,indzh,indexh,ngrid)*uun(ixp,jy ,indzh,indexh,ngrid)+ &
!            uun(ix ,jyp,indzh,indexh,ngrid)*uun(ix ,jyp,indzh,indexh,ngrid)+ &
!            uun(ixp,jyp,indzh,indexh,ngrid)*uun(ixp,jyp,indzh,indexh,ngrid)
!       vsq=vsq+vvn(ix ,jy ,indzh,indexh,ngrid)* &
!            vvn(ix ,jy ,indzh,indexh,ngrid)+ &
!            vvn(ixp,jy ,indzh,indexh,ngrid)*vvn(ixp,jy ,indzh,indexh,ngrid)+ &
!            vvn(ix ,jyp,indzh,indexh,ngrid)*vvn(ix ,jyp,indzh,indexh,ngrid)+ &
!            vvn(ixp,jyp,indzh,indexh,ngrid)*vvn(ixp,jyp,indzh,indexh,ngrid)
!       wsq=wsq+wwn(ix ,jy ,indzh,indexh,ngrid)* &
!            wwn(ix ,jy ,indzh,indexh,ngrid)+ &
!            wwn(ixp,jy ,indzh,indexh,ngrid)*wwn(ixp,jy ,indzh,indexh,ngrid)+ &
!            wwn(ix ,jyp,indzh,indexh,ngrid)*wwn(ix ,jyp,indzh,indexh,ngrid)+ &
!            wwn(ixp,jyp,indzh,indexh,ngrid)*wwn(ixp,jyp,indzh,indexh,ngrid)
!     end do


!   !**********************************
!   ! 2.) Linear vertical interpolation
!   !**********************************

!     uh(m)=dz2*u1(1)+dz1*u1(2)
!     vh(m)=dz2*v1(1)+dz1*v1(2)
!     wh(m)=dz2*w1(1)+dz1*w1(2)
!   end do


!   !************************************
!   ! 3.) Temporal interpolation (linear)
!   !************************************

!   u=(uh(1)*dt2+uh(2)*dt1)*dtt
!   v=(vh(1)*dt2+vh(2)*dt1)*dtt
!   w=(wh(1)*dt2+wh(2)*dt1)*dtt


!   ! Compute standard deviations
!   !****************************

!   xaux=usq-usl*usl/16.
!   if (xaux.lt.eps) then
!     usig=0.
!   else
!     usig=sqrt(xaux/15.)
!   endif

!   xaux=vsq-vsl*vsl/16.
!   if (xaux.lt.eps) then
!     vsig=0.
!   else
!     vsig=sqrt(xaux/15.)
!   endif


!   xaux=wsq-wsl*wsl/16.
!   if (xaux.lt.eps) then
!     wsig=0.
!   else
!     wsig=sqrt(xaux/15.)
!   endif
! end subroutine interpol_wind_nests

! subroutine interpol_rain_nests(yy1,yy2,yy3,nxmaxn,nymaxn,nzmax, &
!        maxnests,ngrid,nxn,nyn,iwftouse,xt,yt,level,itime1,itime2,itime, &
!        yint1,yint2,yint3)
!   !                                i   i   i    i      i      i
!   !   i       i    i   i    i    i  i    i     i      i      i
!   !  o     o     o
!   !****************************************************************************
!   !                                                                           *
!   !  Interpolation of meteorological fields on 2-d model layers for nested    *
!   !  grids. This routine is related to levlin3interpol.f for the mother domain*
!   !                                                                           *
!   !  In horizontal direction bilinear interpolation interpolation is used.    *
!   !  Temporally a linear interpolation is used.                               *
!   !  Three fields are interpolated at the same time.                          *
!   !                                                                           *
!   !  This is a special version of levlininterpol to save CPU time.            *
!   !                                                                           *
!   !  1 first time                                                             *
!   !  2 second time                                                            *
!   !                                                                           *
!   !                                                                           *
!   !     Author: A. Stohl                                                      *
!   !                                                                           *
!   !     15 March 2000                                                         *
!   !                                                                           *
!   !****************************************************************************
!   !                                                                           *
!   ! Variables:                                                                *
!   !                                                                           *
!   ! dt1,dt2              time differences between fields and current position *
!   ! dz1,dz2              z distance between levels and current position       *
!   ! height(nzmax)        heights of the model levels                          *
!   ! indexh               help variable                                        *
!   ! indz                 the level closest to the current trajectory position *
!   ! indzh                help variable                                        *
!   ! itime                current time                                         *
!   ! itime1               time of the first wind field                         *
!   ! itime2               time of the second wind field                        *
!   ! ix,jy                x,y coordinates of lower left subgrid point          *
!   ! level                level at which interpolation shall be done           *
!   ! iwftouse             points to the place of the wind field                *
!   ! nx,ny                actual field dimensions in x,y and z direction       *
!   ! nxmax,nymax,nzmax    maximum field dimensions in x,y and z direction      *
!   ! xt                   current x coordinate                                 *
!   ! yint                 the final interpolated value                         *
!   ! yt                   current y coordinate                                 *
!   ! yy(0:nxmax,0:nymax,nzmax,3) meteorological field used for interpolation   *
!   ! zt                   current z coordinate                                 *
!   !                                                                           *
!   !****************************************************************************
!   use par_mod, only: numwfmem

!   implicit none

!   integer :: maxnests,ngrid
!   integer :: nxn(maxnests),nyn(maxnests),nxmaxn,nymaxn,nzmax,iwftouse
!   integer :: m,ix,jy,ixp,jyp,itime,itime1,itime2,level,indexh
!   real :: yy1(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,maxnests)
!   real :: yy2(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,maxnests)
!   real :: yy3(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,maxnests)
!   real :: ddx,ddy,rddx,rddy,dt1,dt2,dt,y1(2),y2(2),y3(2)
!   real :: xt,yt,yint1,yint2,yint3,p1,p2,p3,p4



!   ! If point at border of grid -> small displacement into grid
!   !***********************************************************

!   ! if (xt.ge.(real(nxn(ngrid)-1)-0.0001)) &
!   !      xt=real(nxn(ngrid)-1)-0.0001
!   ! if (yt.ge.(real(nyn(ngrid)-1)-0.0001)) &
!   !      yt=real(nyn(ngrid)-1)-0.0001

!   ! ESO make it consistent with interpol_rain
!   if (xt.ge.(real(nxn(ngrid)-1))) xt=real(nxn(ngrid)-1)-0.00001
!   if (yt.ge.(real(nyn(ngrid)-1))) yt=real(nyn(ngrid)-1)-0.00001



!   !**********************************************************************
!   ! 1.) Bilinear horizontal interpolation
!   ! This has to be done separately for 2 fields (Temporal)
!   !*******************************************************

!   ! Determine the lower left corner and its distance to the current position
!   !*************************************************************************

!   ix=int(xt)
!   jy=int(yt)

!   ixp=ix+1
!   jyp=jy+1
!   ddx=xt-real(ix)
!   ddy=yt-real(jy)
!   rddx=1.-ddx
!   rddy=1.-ddy
!   p1=rddx*rddy
!   p2=ddx*rddy
!   p3=rddx*ddy
!   p4=ddx*ddy


!   ! Loop over 2 time steps
!   !***********************

!   !  do m=1,2
!   !    indexh=memind(m)
!     indexh=iwftouse

!     y1(1)=p1*yy1(ix ,jy ,level,indexh,ngrid) &
!          + p2*yy1(ixp,jy ,level,indexh,ngrid) &
!          + p3*yy1(ix ,jyp,level,indexh,ngrid) &
!          + p4*yy1(ixp,jyp,level,indexh,ngrid)
!     y2(1)=p1*yy2(ix ,jy ,level,indexh,ngrid) &
!          + p2*yy2(ixp,jy ,level,indexh,ngrid) &
!          + p3*yy2(ix ,jyp,level,indexh,ngrid) &
!          + p4*yy2(ixp,jyp,level,indexh,ngrid)
!     y3(1)=p1*yy3(ix ,jy ,level,indexh,ngrid) &
!          + p2*yy3(ixp,jy ,level,indexh,ngrid) &
!          + p3*yy3(ix ,jyp,level,indexh,ngrid) &
!          + p4*yy3(ixp,jyp,level,indexh,ngrid)
!   !  end do


!   !************************************
!   ! 2.) Temporal interpolation (linear)
!   !************************************

!   ! dt1=real(itime-itime1)
!   ! dt2=real(itime2-itime)
!   ! dt=dt1+dt2

!   ! yint1=(y1(1)*dt2+y1(2)*dt1)/dt
!   ! yint2=(y2(1)*dt2+y2(2)*dt1)/dt
!   ! yint3=(y3(1)*dt2+y3(2)*dt1)/dt

!    yint1=y1(1)
!    yint2=y2(1)
!    yint3=y3(1)
! end subroutine interpol_rain_nests

! subroutine interpol_vdep_nests(level,vdepo)
!   !                                 i     o
!   !****************************************************************************
!   !                                                                           *
!   !  Interpolation of the deposition velocity on 2-d model layer.             *
!   !  In horizontal direction bilinear interpolation interpolation is used.    *
!   !  Temporally a linear interpolation is used.                               *
!   !                                                                           *
!   !  1 first time                                                             *
!   !  2 second time                                                            *
!   !                                                                           *
!   !                                                                           *
!   !     Author: A. Stohl                                                      *
!   !                                                                           *
!   !     30 May 1994                                                           *
!   !                                                                           *
!   !****************************************************************************
!   !                                                                           *
!   ! Variables:                                                                *
!   !                                                                           *
!   ! level                number of species for which interpolation is done    *
!   !                                                                           *
!   !****************************************************************************

!   use par_mod
!   use com_mod
!   use interpol_mod

!   implicit none

!   integer :: level,indexh,m
!   real :: y(2),vdepo

!   ! a) Bilinear horizontal interpolation

!   do m=1,2
!     indexh=memind(m)

!     y(m)=p1*vdepn(ix ,jy ,level,indexh,ngrid) &
!          +p2*vdepn(ixp,jy ,level,indexh,ngrid) &
!          +p3*vdepn(ix ,jyp,level,indexh,ngrid) &
!          +p4*vdepn(ixp,jyp,level,indexh,ngrid)
!   end do


!   ! b) Temporal interpolation

!   vdepo=(y(1)*dt2+y(2)*dt1)*dtt

!   depoindicator(level)=.false.
! end subroutine interpol_vdep_nests

end module interpol_mod



