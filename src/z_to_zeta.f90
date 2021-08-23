! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine z_to_zeta(itime,xt,yt,zold,zteta)
  !                        i    i   o  o  o
  !        o       o       o    i  i  i   o
  !*****************************************************************************
  ! Converting z from eta coordinates to meters                                *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          current temporal position                               *
  ! xteta,yteta,zteta                   spatial position of trajectory         *
  !                                                                            *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use interpol_mod

  implicit none
  integer, intent(in) ::          &
    itime                           ! time index
  integer ::                      &
    i,m,indexh                  ! loop indices
  real(kind=dp), intent(in) ::    &
    xt,yt                           ! particle position
  real, intent(in) ::             &
    zold                            ! particle verticle position in eta coordinates
  real, intent(inout) ::          &
    zteta                           ! converted output z in meters
  real ::                         &
    ttemp_old,ttemp1(2),ttemp_new,& ! storing virtual temperature
    ew,                           & ! why does this function need to be declared here?
    ztemp1,ztemp2,                & ! z positions of the two encompassing levels
    frac,                         & ! fraction between z levels
    psint1(2),psint                 ! pressure of encompassing levels

  call determine_grid_coordinates(real(xt),real(yt))
  call find_grid_distances(real(xt),real(yt))
  call find_time_variables(itime)

  call bilinear_horizontal_interpolation(ps,psint1,1,1)
  call temporal_interpolation(psint1(1),psint1(2),psint)

  call bilinear_horizontal_interpolation(tvirtual,ttemp1,1,nzmax)
  call temporal_interpolation(ttemp1(1),ttemp1(2),ttemp_old)

  ! Integration method as used in the original verttransform_ecmwf.f90
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ztemp1 = 0.
  do i=2,nz-1

    call bilinear_horizontal_interpolation(tvirtual,ttemp1,i,nzmax)
    call temporal_interpolation(ttemp1(1),ttemp1(2),ttemp_new)

    if (abs(ttemp_new-ttemp_old).gt.0.2) then
      ztemp2=ztemp1+r_air/ga*log((akz(i-1)+bkz(i-1)*psint)/(akz(i)+bkz(i)*psint))* &
        (ttemp_new-ttemp_old)/log(ttemp_new/ttemp_old)
    else
      ztemp2=ztemp1+r_air/ga*log((akz(i-1)+bkz(i-1)*psint)/(akz(i)+bkz(i)*psint))*ttemp_new
    endif

    if (ztemp2.gt.zold) then
      frac = (zold-ztemp1)/(ztemp2-ztemp1)
      exit
    else if (i.eq.nz-1) then
      frac = 1.
      exit
    endif
    ttemp_old=ttemp_new
    ztemp1=ztemp2
  end do

  zteta=uvheight(i-1)*(1.-frac)+uvheight(i)*frac

end subroutine z_to_zeta