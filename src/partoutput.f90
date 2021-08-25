! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine partoutput(itime)
  !                        i
  !*****************************************************************************
  !                                                                            *
  !     Dump all particle positions                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     12 March 1999                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use interpol_mod
  use coordinates_ecmwf

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,j,jjjjmmdd,ihmmss
  !integer :: ix,jy,ixp,jyp,indexh,m,il,ind,indz,indzp
  real :: xlon,ylat,ztemp
  real :: topo,hm(2),hmixi,qvi
  real :: tti,rhoi,pvi
  real :: tr(2),tri
  character :: adate*8,atime*6


  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss


  ! Some variables needed for temporal interpolation
  !*************************************************

  ! dt1=real(itime-memtime(1))
  ! dt2=real(memtime(2)-itime)
  ! dtt=1./(dt1+dt2)
  call find_time_variables(itime)

  ! Open output file and write the output
  !**************************************

  if (ipout.eq.1.or.ipout.eq.3) then
    open(unitpartout,file=path(2)(1:length(2))//'partposit_'//adate// &
         atime,form='unformatted')
  else
    open(unitpartout,file=path(2)(1:length(2))//'partposit_end', &
         form='unformatted')
  endif

  ! Write current time to file
  !***************************

  write(unitpartout) itime
  do i=1,numpart

  ! Take only valid particles
  !**************************

    if (itra1(i).eq.itime) then
      xlon=xlon0+xtra1(i)*dx
      ylat=ylat0+ytra1(i)*dy

  !*****************************************************************************
  ! Interpolate several variables (PV, specific humidity, etc.) to particle position
  !*****************************************************************************
      call determine_grid_coordinates(real(xtra1(i)),real(ytra1(i)))
      call find_grid_distances(real(xtra1(i)),real(ytra1(i)))
  ! Topography
  !***********
      call bilinear_horizontal_interpolation_2dim(oro,topo)

      ! First set dz1out from interpol_mod to -1 so it only is calculated once per particle
      !************************************************************************************
      dz1out=-1
      ! Potential vorticity
      call interpol_partoutput_value('PV',pvi,i)
      ! Specific humidity
      call interpol_partoutput_value('QV',qvi,i)
      ! Temperature
      call interpol_partoutput_value('TT',tti,i)
      ! Density
      call interpol_partoutput_value('RH',rhoi,i)
      ! Reset dz1out
      !*************
      dz1out=-1

  ! Tropopause and PBL height
  !**************************
  ! Tropopause
      call bilinear_horizontal_interpolation(tropopause,tr,1,1)
      call temporal_interpolation(tr(1),tr(2),tri)
  ! PBL height
      call bilinear_horizontal_interpolation(hmix,hm,1,1)
      call temporal_interpolation(hm(1),hm(2),hmixi)


  ! Convert eta z coordinate to meters if necessary
  !************************************************
      if (wind_coord_type.eq.'ETA') then
        call zeta_to_z(itime,xtra1(i),ytra1(i),ztra1eta(i),ztemp)
      else
        ztemp=ztra1(i)
      endif
  ! Write the output
  !*****************      
      write(unitpartout) npoint(i),xlon,ylat,ztemp, &
           itramem(i),topo,pvi,qvi,rhoi,hmixi,tri,tti, &
           (xmass1(i,j),j=1,nspec)
    endif
  end do
  write(unitpartout) -99999,-9999.9,-9999.9,-9999.9,-99999, &
       -9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9, &
       (-9999.9,j=1,nspec)


  close(unitpartout)

end subroutine partoutput
