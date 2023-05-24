! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine interpol_all(itime,xt,yt,zt,zteta)
  !                          i   i  i  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates everything that is needed for calculating the*
  !  dispersion.                                                               *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
  !                                                                            *
  !  Revision March 2005 by AST : all output variables in common block cal-    *
  !                               culation of standard deviation done in this  *
  !                               routine rather than subroutine call in order *
  !                               to save computation time                     *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          current temporal position                               *
  ! memtime(3) [s]     times of the wind fields in memory                      *
  ! xt,yt,zt           coordinates position for which wind data shall be       *
  !                    culated                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use interpol_mod
  use hanna_mod

  implicit none

  integer :: itime
  real :: xt,yt,zt
  real:: zteta

  ! Auxiliary variables needed for interpolation
  real :: ust1(2),wst1(2),oli1(2),oliaux
  real :: y1(2),y2(2),y3(2),rho1(2),rhograd1(2),psint(2)
  real :: usl,vsl,wsl,usq,vsq,wsq,xaux,deta1(2),detatemp
  integer :: i,m,n,indexh
  real,parameter :: eps=1.0e-30


  !********************************************
  ! Multilinear interpolation in time and space
  !********************************************

  ! Determine the lower left corner and its distance to the current position
  !*************************************************************************

  ddx=xt-real(ix)
  ddy=yt-real(jy)
  rddx=1.-ddx
  rddy=1.-ddy
  p1=rddx*rddy
  p2=ddx*rddy
  p3=rddx*ddy
  p4=ddx*ddy

  ! Calculate variables for time interpolation
  !*******************************************

  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)


  !*****************************************
  ! 1. Interpolate u*, w* and Obukhov length
  !*****************************************

  ! a) Bilinear horizontal interpolation

  do m=1,2
    indexh=memind(m)

    ust1(m)=p1*ustar(ix ,jy ,1,indexh) &
         + p2*ustar(ixp,jy ,1,indexh) &
         + p3*ustar(ix ,jyp,1,indexh) &
         + p4*ustar(ixp,jyp,1,indexh)
    wst1(m)=p1*wstar(ix ,jy ,1,indexh) &
         + p2*wstar(ixp,jy ,1,indexh) &
         + p3*wstar(ix ,jyp,1,indexh) &
         + p4*wstar(ixp,jyp,1,indexh)
    oli1(m)=p1*oli(ix ,jy ,1,indexh) &
         + p2*oli(ixp,jy ,1,indexh) &
         + p3*oli(ix ,jyp,1,indexh) &
         + p4*oli(ixp,jyp,1,indexh)
  end do

  ! b) Temporal interpolation

  ust=(ust1(1)*dt2+ust1(2)*dt1)*dtt
  wst=(wst1(1)*dt2+wst1(2)*dt1)*dtt
  oliaux=(oli1(1)*dt2+oli1(2)*dt1)*dtt

  if (oliaux.ne.0.) then
    ol=1./oliaux
  else
    ol=99999.
  endif


  !*****************************************************
  ! 2. Interpolate vertical profiles of u,v,w,rho,drhodz
  !*****************************************************


  ! Determine the level below the current position
  !***********************************************

  do i=2,nz
    if (height(i).gt.zt) then
      indz=i-1
      indzp=i
      exit
    endif
  end do

  !**************************************
  ! 1.) Bilinear horizontal interpolation
  ! 2.) Temporal interpolation (linear)
  !**************************************

  ! Loop over 2 time steps and indz levels
  !***************************************

  do n=indz,indzp
    wsl=0.
    wsq=0.
    do m=1,2
      indexh=memind(m)

      y3(m)=p1*ww(ix ,jy ,n,indexh) &
           +p2*ww(ixp,jy ,n,indexh) &
           +p3*ww(ix ,jyp,n,indexh) &
           +p4*ww(ixp,jyp,n,indexh)
      wsl=wsl+ww(ix ,jy ,n,indexh)+ww(ixp,jy ,n,indexh) &
           +ww(ix ,jyp,n,indexh)+ww(ixp,jyp,n,indexh)
      wsq=wsq+ww(ix ,jy ,n,indexh)*ww(ix ,jy ,n,indexh)+ &
           ww(ixp,jy ,n,indexh)*ww(ixp,jy ,n,indexh)+ &
           ww(ix ,jyp,n,indexh)*ww(ix ,jyp,n,indexh)+ &
           ww(ixp,jyp,n,indexh)*ww(ixp,jyp,n,indexh)
    end do
    wprof(n)=(y3(1)*dt2+y3(2)*dt1)*dtt
    indzindicator(n)=.false.

  ! Compute standard deviations
  !****************************

    xaux=wsq-wsl*wsl/8.
    if (xaux.lt.eps) then
      wsigprof(n)=0.
    else
      wsigprof(n)=sqrt(xaux/7.)
    endif

  end do


  ! Same for zt in eta coordinates
  !*******************************
  indzeta=nz-1
  indzpeta=nz
  do i=2,nz
    if (wheight(i).lt.zteta) then
      indzeta=i-1
      indzpeta=i
      exit
    endif
  end do

  induv=nz-1
  indpuv=nz
  do i=2,nz
    if (uvheight(i).lt.zteta) then
      induv=i-1
      indpuv=i
      exit
    endif
  end do

  !**************************************
  ! 1.) Bilinear horizontal interpolation
  ! 2.) Temporal interpolation (linear)
  !**************************************

  ! Loop over 2 time steps and indz levels
  !***************************************

  do n=induv,indpuv
    usl=0.
    vsl=0.
    usq=0.
    vsq=0.
    do m=1,2
      indexh=memind(m)
      if (ngrid.lt.0) then
        y1(m)=p1*uupoleta(ix ,jy ,n,indexh) &
             +p2*uupoleta(ixp,jy ,n,indexh) &
             +p3*uupoleta(ix ,jyp,n,indexh) &
             +p4*uupoleta(ixp,jyp,n,indexh)
        y2(m)=p1*vvpoleta(ix ,jy ,n,indexh) &
             +p2*vvpoleta(ixp,jy ,n,indexh) &
             +p3*vvpoleta(ix ,jyp,n,indexh) &
             +p4*vvpoleta(ixp,jyp,n,indexh)
        usl=usl+uupoleta(ix ,jy ,n,indexh)+uupoleta(ixp,jy ,n,indexh) &
             +uupoleta(ix ,jyp,n,indexh)+uupoleta(ixp,jyp,n,indexh)
        vsl=vsl+vvpoleta(ix ,jy ,n,indexh)+vvpoleta(ixp,jy ,n,indexh) &
             +vvpoleta(ix ,jyp,n,indexh)+vvpoleta(ixp,jyp,n,indexh)

        usq=usq+uupoleta(ix ,jy ,n,indexh)*uupoleta(ix ,jy ,n,indexh)+ &
             uupoleta(ixp,jy ,n,indexh)*uupoleta(ixp,jy ,n,indexh)+ &
             uupoleta(ix ,jyp,n,indexh)*uupoleta(ix ,jyp,n,indexh)+ &
             uupoleta(ixp,jyp,n,indexh)*uupoleta(ixp,jyp,n,indexh)
        vsq=vsq+vvpoleta(ix ,jy ,n,indexh)*vvpoleta(ix ,jy ,n,indexh)+ &
             vvpoleta(ixp,jy ,n,indexh)*vvpoleta(ixp,jy ,n,indexh)+ &
             vvpoleta(ix ,jyp,n,indexh)*vvpoleta(ix ,jyp,n,indexh)+ &
             vvpoleta(ixp,jyp,n,indexh)*vvpoleta(ixp,jyp,n,indexh)
      else
        y1(m)=p1*uueta(ix ,jy ,n,indexh) &
             +p2*uueta(ixp,jy ,n,indexh) &
             +p3*uueta(ix ,jyp,n,indexh) &
             +p4*uueta(ixp,jyp,n,indexh)
        y2(m)=p1*vveta(ix ,jy ,n,indexh) &
             +p2*vveta(ixp,jy ,n,indexh) &
             +p3*vveta(ix ,jyp,n,indexh) &
             +p4*vveta(ixp,jyp,n,indexh)
        usl=usl+uueta(ix ,jy ,n,indexh)+uueta(ixp,jy ,n,indexh) &
             +uueta(ix ,jyp,n,indexh)+uueta(ixp,jyp,n,indexh)
        vsl=vsl+vveta(ix ,jy ,n,indexh)+vveta(ixp,jy ,n,indexh) &
             +vveta(ix ,jyp,n,indexh)+vveta(ixp,jyp,n,indexh)

        usq=usq+uueta(ix ,jy ,n,indexh)*uueta(ix ,jy ,n,indexh)+ &
             uueta(ixp,jy ,n,indexh)*uueta(ixp,jy ,n,indexh)+ &
             uueta(ix ,jyp,n,indexh)*uueta(ix ,jyp,n,indexh)+ &
             uueta(ixp,jyp,n,indexh)*uueta(ixp,jyp,n,indexh)
        vsq=vsq+vveta(ix ,jy ,n,indexh)*vveta(ix ,jy ,n,indexh)+ &
             vveta(ixp,jy ,n,indexh)*vveta(ixp,jy ,n,indexh)+ &
             vveta(ix ,jyp,n,indexh)*vveta(ix ,jyp,n,indexh)+ &
             vveta(ixp,jyp,n,indexh)*vveta(ixp,jyp,n,indexh)
      endif
      rhograd1(m)=p1*drhodzeta(ix ,jy ,n,indexh) &
           +p2*drhodzeta(ixp,jy ,n,indexh) &
           +p3*drhodzeta(ix ,jyp,n,indexh) &
           +p4*drhodzeta(ixp,jyp,n,indexh)
      rho1(m)=p1*rhoeta(ix ,jy ,n,indexh) &
           +p2*rhoeta(ixp,jy ,n,indexh) &
           +p3*rhoeta(ix ,jyp,n,indexh) &
           +p4*rhoeta(ixp,jyp,n,indexh)
    end do
    uprof(n)=(y1(1)*dt2+y1(2)*dt1)*dtt
    vprof(n)=(y2(1)*dt2+y2(2)*dt1)*dtt
    rhoprof(n)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
    rhogradprof(n)=(rhograd1(1)*dt2+rhograd1(2)*dt1)*dtt

    xaux=usq-usl*usl/8.
    if (xaux.lt.eps) then
      usigprof(n)=0.
    else
      usigprof(n)=sqrt(xaux/7.)
    endif

    xaux=vsq-vsl*vsl/8.
    if (xaux.lt.eps) then
      vsigprof(n)=0.
    else
      vsigprof(n)=sqrt(xaux/7.)
    endif

  end do

  do n=indzeta,indzpeta
    wsl=0.
    wsq=0.
    do m=1,2
      indexh=memind(m)

      y3(m)=p1*wweta(ix ,jy ,n,indexh) &
           +p2*wweta(ixp,jy ,n,indexh) &
           +p3*wweta(ix ,jyp,n,indexh) &
           +p4*wweta(ixp,jyp,n,indexh)

      wsl=wsl+wweta(ix ,jy ,n,indexh)+wweta(ixp,jy ,n,indexh) &
           +wweta(ix ,jyp,n,indexh)+wweta(ixp,jyp,n,indexh)
      wsq=wsq+wweta(ix ,jy ,n,indexh)*wweta(ix ,jy ,n,indexh)+ &
           wweta(ixp,jy ,n,indexh)*wweta(ixp,jy ,n,indexh)+ &
           wweta(ix ,jyp,n,indexh)*wweta(ix ,jyp,n,indexh)+ &
           wweta(ixp,jyp,n,indexh)*wweta(ixp,jyp,n,indexh)

    end do
    wprofeta(n)=(y3(1)*dt2+y3(2)*dt1)*dtt

    xaux=wsq-wsl*wsl/8.
    if (xaux.lt.eps) then
      wsigprof(n)=0.
    else
      wsigprof(n)=sqrt(xaux/7.)
    endif

  end do
end subroutine interpol_all
