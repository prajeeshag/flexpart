! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

  !*****************************************************************************
  !                                                                            *
  ! This module contains all subroutines computing the vertical coordinate     *
  ! transformation of the meteorological input data (L. Bakels 2021)           *
  !                                                                            *
  !*****************************************************************************

module verttransform_mod
  use par_mod
  use com_mod
  use qvsat_mod
  use cmapf_mod, only: cc2gll
  use windfields_mod

  implicit none

contains


subroutine verttransform_init(n)
  implicit none

  integer, intent(in) :: n
  real ::  tvold,pold,pint,tv
  integer :: ix,jy,kz,ixm,jym
  real,parameter :: const=r_air/ga

  if ((ipin.eq.1).or.(ipin.eq.4)) then
    call read_heightlevels(height,nmixz)
    return
  endif

  jym=nymin1
  ixm=0
  loop1: do jy=0,nymin1
    do ix=0,nxmin1
      if (ps(ix,jy,1,n).gt.100000.) then
        ixm=ix
        jym=jy
        exit loop1
      endif
    end do
  end do loop1

  tvold=tt2(ixm,jym,1,n)*(1.+0.378*ew(td2(ixm,jym,1,n),ps(ixm,jym,1,n))/ &
       ps(ixm,jym,1,n))
  pold=ps(ixm,jym,1,n)
  height(1)=0.

  do kz=2,nuvz
    pint=akz(kz)+bkz(kz)*ps(ixm,jym,1,n)
    tv=tth(ixm,jym,kz,n)*(1.+0.608*qvh(ixm,jym,kz,n))

    if (abs(tv-tvold).gt.0.2) then
      height(kz)= height(kz-1)+const*log(pold/pint)*(tv-tvold)/log(tv/tvold)
    else
      height(kz)=height(kz-1)+const*log(pold/pint)*tv
    endif

    tvold=tv
    pold=pint
  end do

  ! Determine highest levels that can be within PBL
  !************************************************

  do kz=1,nz
    if (height(kz).gt.hmixmax) then
      nmixz=kz
      exit
    endif
  end do

  if (loutrestart.ne.-1) then
    call output_heightlevels(height,nmixz)
  endif
end subroutine verttransform_init

subroutine output_heightlevels(height_tmp,nmixz_tmp)
  implicit none

  real,intent(in) :: height_tmp(nzmax)
  integer,intent(in) :: nmixz_tmp
  integer :: kz
  character(len=256) :: heightlevels_filename

  heightlevels_filename = path(2)(1:length(2))//'heightlevels.bin'

  write(*,*) 'Writing Initialised heightlevels to file:', trim(heightlevels_filename)

  open(unitheightlevels,file=trim(heightlevels_filename),form='unformatted')

  write(unitheightlevels) nmixz_tmp

  do kz=1,nz
    write(unitheightlevels) height_tmp(kz)
  end do
  close(unitheightlevels)
end subroutine output_heightlevels

subroutine read_heightlevels(height_tmp,nmixz_tmp)
  implicit none

  real,intent(out) :: height_tmp(nzmax)
  integer,intent(out) :: nmixz_tmp
  integer :: kz,ios
  character(len=256) :: heightlevels_filename

  heightlevels_filename = path(2)(1:length(2))//'heightlevels.bin'

  write(*,*) 'Reading heightlevels from file:', trim(heightlevels_filename)

  open(unitheightlevels,file=trim(heightlevels_filename),form='unformatted',err=9988)

  read(unitheightlevels,iostat=ios) nmixz_tmp

  do kz=1,nz
    read(unitheightlevels) height_tmp(kz)
  end do
  close(unitheightlevels)

  return

9988   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE             #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'heightlevels.bin'//'    #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS        #### '
  write(*,*) ' #### NAME DOES NOT EXISTS, REMOVE call read_heightlevels #### '
  write(*,*) ' #### FROM VERTTRANSFORM_MOD.                 #### '
end subroutine read_heightlevels



subroutine convert_cloud_params(ix,jy,nxlim,nylim,max_cloudthck_eta,conv_clrange_eta, &
  highconvp_clrange_eta,lowconvp_clrange_eta,uvzlev)

  implicit none

  integer, intent(in) :: ix,jy,nxlim,nylim
  real,intent(in),dimension(0:nxlim,0:nylim,nzmax) :: uvzlev

  ! converted parameters for eta coordinates:
  integer,intent(out) :: max_cloudthck_eta
  integer, dimension(2),intent(out) :: conv_clrange_eta,highconvp_clrange_eta, &
    lowconvp_clrange_eta
  integer :: i, kz


  ! Convert cloud parameters to eta coords.
  ! Reverse sign when using eta (eta:1-0, meter:0-max)
  max_cloudthck_eta=int(uvheight(nz)*eta_convert)
  conv_clrange_eta=int(uvheight(nz)*eta_convert)
  highconvp_clrange_eta=int(uvheight(nz)*eta_convert)
  lowconvp_clrange_eta=int(uvheight(nz)*eta_convert)
  do kz=1,nz
    if (uvzlev(ix,jy,kz).gt.max_cloudthck) then
      max_cloudthck_eta=int(uvheight(kz)*eta_convert)
      exit
    endif
  end do
  do i=1,2
    do kz=1,nz
      if (uvzlev(ix,jy,kz).gt.conv_clrange(i)) then
        conv_clrange_eta(i)=int(uvheight(kz)*eta_convert)
        exit
      endif
    end do
    do kz=1,nz
      if (uvzlev(ix,jy,kz).gt.highconvp_clrange(i)) then
        highconvp_clrange_eta(i)=int(uvheight(kz)*eta_convert)
        exit
      endif
    end do
    do kz=1,nz
      if (uvzlev(ix,jy,kz).gt.lowconvp_clrange(i)) then
        lowconvp_clrange_eta(i)=int(uvheight(kz)*eta_convert)
        exit
      endif
    end do
  end do
end subroutine convert_cloud_params

subroutine identify_cloud(ix,jy,lcw_tmp,lcwsum_tmp,nxlim,nylim, &
  ctwc_tmp,clwc_tmp,ciwc_tmp,icloudbot_tmp,icloudtop_tmp,rho_tmp, &
  tt_tmp,qv_tmp,uvzlev,wzlev)

  implicit none

  logical,intent(in) :: lcw_tmp,lcwsum_tmp
  integer, intent(in) :: ix,jy,nxlim,nylim
  real,intent(out) :: ctwc_tmp(0:nxlim,0:nylim)

  real,intent(inout) :: clwc_tmp(0:nxlim,0:nylim,nzmax)
  real,intent(in) :: ciwc_tmp(0:nxlim,0:nylim,nzmax)
  real,intent(in),dimension(0:nxlim,0:nylim,nzmax) :: rho_tmp,tt_tmp,qv_tmp
  real,intent(in),dimension(0:nxlim,0:nylim,nzmax) :: uvzlev,wzlev

  integer,intent(out) :: icloudbot_tmp(0:nxlim,0:nylim), icloudtop_tmp(0:nxlim,0:nylim)

  integer :: kz
  real :: pressure,rh
  real :: clw

  !*******************************************************************************
  if (lcw_tmp) then ! identify clouds based on cloud water content
  !*******************************************************************************

    ctwc_tmp(ix,jy) = 0. ! initialise cloud total water content
    if (.not.lcwsum_tmp) clwc_tmp(ix,jy,:) = clwc_tmp(ix,jy,:) + ciwc_tmp(ix,jy,:)
    do kz = 1,nz-1 ! Changed order of loop to prevent ETA computation to be done unnecessarily

      ! vertically integrate cloud water and determine cloud bottom, top
      ! cloud water per cell in kg / m2
      ! calculate cloud water mass per area: kgCW/kgAIR * kgAIR/m3 * m = kgCW/m2

      ! assuming rho is in kg/m3 and hz in m gives: kg/kg * kg/m3 *m3/kg /m = m2/m3



      clw = clwc_tmp(ix,jy,kz)*rho_tmp(ix,jy,kz)*(height(kz+1)-height(kz))

      ! Add this layer to column cloud water [m3/m3]
      ctwc_tmp(ix,jy) = ctwc_tmp(ix,jy)+clw ! kg / m2 (or eta) in column

      if (clw .gt. 0.) then ! cloud layer - maybe use threshold?





        if (icloudbot_tmp(ix,jy) .eq. icmv) &
            icloudbot_tmp(ix,jy) = (height(kz))
          icloudtop_tmp(ix,jy) = (height(kz))

      endif
    end do

  !**************************************************************************
  else       ! identify clouds using relative humidity
  !**************************************************************************
    do kz = 1,nz-1 ! Changed order of loop to prevent ETA computation to be done unnecessarily
      pressure=rho_tmp(ix,jy,kz)*r_air*tt_tmp(ix,jy,kz)
      rh=qv_tmp(ix,jy,kz)/f_qvsat(pressure,tt_tmp(ix,jy,kz))
      ! PS if (prec.gt.0.01) print*,'relhum',prec,kz,rh,height(kz)
      if (rh .ge. rhmin) then






        if (icloudbot_tmp(ix,jy) .eq. icmv) then
          icloudbot_tmp(ix,jy)=(height(kz))
        endif
        icloudtop_tmp(ix,jy)=(height(kz))


      endif
    end do
!**************************************************************************
  endif ! lcw true/false
!**************************************************************************
end subroutine identify_cloud

subroutine apply_cloud_bounds(ix,jy,nxlim,nylim,lsprec_tmp,convprec_tmp,uvzlev, &
  icloudbot_tmp,icloudtop_tmp,max_cloudthck_eta,conv_clrange_eta, &
  highconvp_clrange_eta,lowconvp_clrange_eta)

  implicit none

  integer, intent(in) :: ix,jy,nxlim,nylim

  real,intent(in) :: lsprec_tmp(0:nxlim,0:nylim,numpf),convprec_tmp(0:nxlim,0:nylim,numpf)
  real,intent(in),dimension(0:nxlim,0:nylim,nzmax) :: uvzlev

  ! converted parameters for eta coordinates:
  integer,intent(in) :: max_cloudthck_eta
  integer,intent(in),dimension(2) :: conv_clrange_eta,highconvp_clrange_eta,lowconvp_clrange_eta

  integer,intent(inout) :: icloudbot_tmp(0:nxlim,0:nylim), icloudtop_tmp(0:nxlim,0:nylim)

  integer :: icloudtop_old, min_cloudtop
  integer :: kz
  real :: lsp,convp,prec
  logical lconvectprec


  real,parameter :: precmin = 0.002 ! minimum prec in mm/h for cloud diagn.


  ! memorise icloudtop
  icloudtop_old = icloudtop_tmp(ix,jy)
  ! top level, kz=nz-1
  ! limit cloud top to 19 km:
  if (icloudbot_tmp(ix,jy) .eq. icmv) then ! if no bottom found, no top either
    icloudtop_tmp(ix,jy)=icmv
  else if (icloudtop_tmp(ix,jy) .gt. max_cloudthck) then ! max cloud height
    icloudtop_tmp(ix,jy) = max_cloudthck
  endif
 ! PS  get rid of too thin clouds
  if (icloudtop_tmp(ix,jy) .lt. icloudbot_tmp(ix,jy) + min_cloudthck) then
    icloudbot_tmp(ix,jy)=icmv
    icloudtop_tmp(ix,jy)=icmv
  endif
  ! PS implement a rough fix for badly represented convection
  ! PS is based on looking at a limited set of comparison data
  lsp=  sum(   lsprec_tmp(ix,jy,:) )
  convp=sum( convprec_tmp(ix,jy,:) )
  prec=lsp+convp
  if (lsp.gt.convp) then !  prectype='lsp'
    lconvectprec = .false.
  else ! prectype='cp '
    lconvectprec = .true.
  endif
  if (lconvectprec .and. prec .gt. precmin .and.  &
    (icloudtop_old .lt. conv_clrange(2) .or. &
      icloudbot_tmp(ix,jy) .gt. conv_clrange(1)) ) then
    if (convp .lt. 0.1) then
      icloudbot_tmp(ix,jy) = lowconvp_clrange(1)
      icloudtop_tmp(ix,jy) = lowconvp_clrange(2)
    else
      icloudbot_tmp(ix,jy) = highconvp_clrange(1)
      icloudtop_tmp(ix,jy) = highconvp_clrange(2)
    endif
  endif

end subroutine apply_cloud_bounds

subroutine verttransform_gfs(n,uuh,vvh,wwh,pvh)
  !                          i  i   i   i   i
  !*****************************************************************************
  !                                                                            *
  !     This subroutine transforms temperature, dew point temperature and      *
  !     wind components from eta to meter coordinates.                         *
  !     The vertical wind component is transformed from Pa/s to m/s using      *
  !     the conversion factor pinmconv.                                        *
  !     In addition, this routine calculates vertical density gradients        *
  !     needed for the parameterization of the turbulent velocities.           *
  !                                                                            *
  !     Author: A. Stohl, G. Wotawa                                            *
  !                                                                            *
  !     12 August 1996                                                         *
  !     Update: 16 January 1998                                                *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  !  CHANGES                                                                   *
  !     Major update: 17 February 1999                                         *
  !     by G. Wotawa                                                           *
  !                                                                            *
  !     - Vertical levels for u, v and w are put together                      *
  !     - Slope correction for vertical velocity: Modification of calculation  *
  !       procedure                                                            *
  !                                                                            *
  !  Bernd C. Krueger, Feb. 2001:                                              *
  !   Variables tth and qvh (on eta coordinates) from common block             *
  !                                                                            *
  ! Sabine Eckhardt, March 2007:                                               *
  ! added the variable cloud for use with scavenging - descr. in com_mod       *
  ! PS/AT 2018/-21: variable "cloud" is replaced by quickfix, see below        *
  !                                                                            *
  ! Unified ECMWF and GFS builds                                               *
  ! Marian Harustak, 12.5.2017                                                 *
  !     - Renamed from verttransform to verttransform_ecmwf                    *
  !                                                                            *
  !  undocumented modifications by NILU for v10                                *
  !                                                                            *
  !  Petra Seibert, 2018-06-13:                                                *
  !   - put back SAVE attribute for INIT, just to be safe                      *
  !   - minor changes, most of them just cosmetics                             *
  !   for details see changelog.txt in branch unive                            *
  !                                                                            *
  !  Petra Seibert, Anne Philipp, 2019-05-02: implement wetdepo quickfix       *
  !  Petra Seibert, Anne Tipka, 2020-11-19: reimplement in latest version      *
  !                                                                            *
  ! ****************************************************************************

  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! Note PS, AT 2021-01-29: all these fields are 0:nxmax-1,0:nymax-1 !!        *
  ! nx,ny,nz                        field dimensions in x,y and z direction    *
  ! icloudbot(0:nxmax,0:nymax,numwfmem) cloud bottom field for wet deposition  *
  ! icloudtop(0:nxmax,0:nymax,numwfmem) cloud thickness for wet deposition    *
  ! uu(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in x-direction [m/s]*
  ! vv(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in y-direction [m/s]*
  ! ww(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in z-direction      *
  !                                          [deltaeta/s]                      *
  ! tt(0:nxmax,0:nymax,nzmax,numwfmem)     temperature [K]                     *
  ! pv(0:nxmax,0:nymax,nzmax,numwfmem)     potential voriticity (pvu)          *
  ! ps(0:nxmax,0:nymax,numwfmem)           surface pressure [Pa]               *
  !                                                                            *
  !*****************************************************************************

  implicit none

  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: uuh,vvh,pvh
  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nwzmax)  :: wwh

  real,dimension(0:nxmax-1,0:nymax-1,nzmax) :: uvwzlev
  real,dimension(nuvzmax) :: rhoh
  real,dimension(nwzmax) :: wzlev
  real,dimension(nzmax) :: pinmconv

! local array introduced in v10 by ?? to achieve better loop order (PS)
  integer,dimension(0:nxmax-1,0:nymax-1) :: idx

  integer :: icloudtop_old
  integer :: ix,jy,kz,iz,n,kmin,kl,klp,ix1,jy1,ixp,jyp,ixm,jym,kz_inv
  real :: clw ! cloud water in kg / m2 in a grid cell
  real :: pressure,rh,lsp,convp,prec

  real :: pint,tv,tvold,pold,dz1,dz2,dz,ui,vi
  real :: xlon,ylat,xlonr,dzdx,dzdy
  real :: dzdx1,dzdx2,dzdy1,dzdy2,cosf
  real :: uuaux,vvaux,uupolaux,vvpolaux,ddpol,ffpol,wdummy

  real,parameter :: const=r_air/ga
  real,parameter :: precmin = 0.002 ! minimum prec in mm/h for cloud diagnostics

  logical lconvectprec
  logical, save :: init = .true. !PS 2018-06-13: add back save attribute

  ! NCEP version
  integer :: llev, i

  integer :: rank_thread , nr_threads


  !*************************************************************************
  ! If verttransform is called the first time, initialize heights of the   *
  ! z levels in meter. The heights are the heights of model levels, where  *
  ! u,v,T and qv are given.                                                *
  !*************************************************************************

  if (init) then

  ! Search for a point with high surface pressure (i.e. not above significant topography)
  ! Then, use this point to construct a reference z profile, to be used at all times
  !*****************************************************************************
    call verttransform_init(n)

  ! Do not repeat initialization of the Cartesian z grid
  !*****************************************************

    init=.false.

  endif


  ! Loop over the whole grid
  !*************************
!$OMP DO
  do jy=0,nymin1
    do ix=0,nxmin1


!   if ((jy.eq.0).and.(ix.eq.0)) print*, 'in loop 1'

  ! NCEP version: find first level above ground
      llev = 0
      do i=1,nuvz
        if (ps(ix,jy,1,n).lt.akz(i)) llev=i
      enddo
       llev = llev+1
       if (llev.gt.nuvz-2) llev = nuvz-2
  !     if (llev.eq.nuvz-2) write(*,*) 'verttransform
  !    +WARNING: LLEV eq NUZV-2'
  ! NCEP version


  ! compute height of pressure levels above ground
  !***********************************************

      tvold=tth(ix,jy,llev,n)*(1.+0.608*qvh(ix,jy,llev,n))
      pold=akz(llev)
      wzlev(llev)=0.
      uvwzlev(ix,jy,llev)=0.
      rhoh(llev)=pold/(r_air*tvold)

!    if ((jy.eq.0).and.(ix.eq.0)) print*, 'in loop 2'

      do kz=llev+1,nuvz
        pint=akz(kz)+bkz(kz)*ps(ix,jy,1,n)
        tv=tth(ix,jy,kz,n)*(1.+0.608*qvh(ix,jy,kz,n))
        rhoh(kz)=pint/(r_air*tv)

        if (abs(tv-tvold).gt.0.2) then
          uvwzlev(ix,jy,kz)=uvwzlev(ix,jy,kz-1)+const*log(pold/pint)* &
          (tv-tvold)/log(tv/tvold)
        else
          uvwzlev(ix,jy,kz)=uvwzlev(ix,jy,kz-1)+const*log(pold/pint)*tv
        endif
        wzlev(kz)=uvwzlev(ix,jy,kz)

        tvold=tv
        pold=pint
      enddo

  ! pinmconv=(h2-h1)/(p2-p1)
! if ((jy.eq.0).and.(ix.eq.0)) print*, 'in loop 3'

      pinmconv(llev)=(uvwzlev(ix,jy,llev+1)-uvwzlev(ix,jy,llev))/ &
           ((aknew(llev+1)+bknew(llev+1)*ps(ix,jy,1,n))- &
           (aknew(llev)+bknew(llev)*ps(ix,jy,1,n)))
      do kz=llev+1,nz-1
        pinmconv(kz)=(uvwzlev(ix,jy,kz+1)-uvwzlev(ix,jy,kz-1))/ &
             ((aknew(kz+1)+bknew(kz+1)*ps(ix,jy,1,n))- &
             (aknew(kz-1)+bknew(kz-1)*ps(ix,jy,1,n)))
      enddo
      pinmconv(nz)=(uvwzlev(ix,jy,nz)-uvwzlev(ix,jy,nz-1))/ &
           ((aknew(nz)+bknew(nz)*ps(ix,jy,1,n))- &
           (aknew(nz-1)+bknew(nz-1)*ps(ix,jy,1,n)))


  ! Levels, where u,v,t and q are given
  !************************************
! if ((jy.eq.0).and.(ix.eq.0)) print*, 'in loop 4'


      uu(ix,jy,1,n)=uuh(ix,jy,llev)
      vv(ix,jy,1,n)=vvh(ix,jy,llev)
      tt(ix,jy,1,n)=tth(ix,jy,llev,n)
      qv(ix,jy,1,n)=qvh(ix,jy,llev,n)


  ! IP & SEC, 201812 add clouds
      if (lcw) then
         clwc(ix,jy,1,n)=clwch(ix,jy,llev,n)
      endif
      pv(ix,jy,1,n)=pvh(ix,jy,llev)
      rho(ix,jy,1,n)=rhoh(llev)
      pplev(ix,jy,1,n)=akz(llev)
      uu(ix,jy,nz,n)=uuh(ix,jy,nuvz)
      vv(ix,jy,nz,n)=vvh(ix,jy,nuvz)
      tt(ix,jy,nz,n)=tth(ix,jy,nuvz,n)
      qv(ix,jy,nz,n)=qvh(ix,jy,nuvz,n)
  ! IP & SEC, 201812 add clouds


      if (lcw) then
         clwc(ix,jy,nz,n)=clwch(ix,jy,nuvz,n)
      endif
      pv(ix,jy,nz,n)=pvh(ix,jy,nuvz)
      rho(ix,jy,nz,n)=rhoh(nuvz)
      pplev(ix,jy,nz,n)=akz(nuvz)
      kmin=llev+1


      do iz=2,nz-1
        do kz=kmin,nuvz
          ! print*, 'in loop 4.3.1', jy, ix, iz, kz
          if (height(iz).gt.uvwzlev(ix,jy,nuvz)) then
            uu(ix,jy,iz,n)=uu(ix,jy,nz,n)
            vv(ix,jy,iz,n)=vv(ix,jy,nz,n)
            tt(ix,jy,iz,n)=tt(ix,jy,nz,n)
            qv(ix,jy,iz,n)=qv(ix,jy,nz,n)
  ! IP & SEC, 201812 add clouds
            if (lcw) then
               clwc(ix,jy,iz,n)=clwc(ix,jy,nz,n)
            endif
            pv(ix,jy,iz,n)=pv(ix,jy,nz,n)
            rho(ix,jy,iz,n)=rho(ix,jy,nz,n)
            pplev(ix,jy,iz,n)=pplev(ix,jy,nz,n)
            exit
          endif
          if ((height(iz).gt.uvwzlev(ix,jy,kz-1)).and. &
           (height(iz).le.uvwzlev(ix,jy,kz))) then
            !real,dimension(0:nxmax-1,0:nymax-1,nzmax) :: uvwzlev
            dz1=height(iz)-uvwzlev(ix,jy,kz-1)
            dz2=uvwzlev(ix,jy,kz)-height(iz)
            dz=dz1+dz2
            uu(ix,jy,iz,n)=(uuh(ix,jy,kz-1)*dz2+uuh(ix,jy,kz)*dz1)/dz
            vv(ix,jy,iz,n)=(vvh(ix,jy,kz-1)*dz2+vvh(ix,jy,kz)*dz1)/dz
            tt(ix,jy,iz,n)=(tth(ix,jy,kz-1,n)*dz2+tth(ix,jy,kz,n)*dz1)/dz
            qv(ix,jy,iz,n)=(qvh(ix,jy,kz-1,n)*dz2+qvh(ix,jy,kz,n)*dz1)/dz

  ! IP & SEC, 201812 add clouds
            if (lcw) then
              clwc(ix,jy,iz,n)= &
                (clwch(ix,jy,kz-1,n)*dz2+clwch(ix,jy,kz,n)*dz1)/dz
            endif
            pv(ix,jy,iz,n)=(pvh(ix,jy,kz-1)*dz2+pvh(ix,jy,kz)*dz1)/dz

            rho(ix,jy,iz,n)=(rhoh(kz-1)*dz2+rhoh(kz)*dz1)/dz

            pplev(ix,jy,iz,n)=(akz(kz-1)*dz2+akz(kz)*dz1)/dz

          endif
        enddo
      enddo


  ! Interpolation of vertical motion (levels where w is given)
  !***********************************************************

      ww(ix,jy,1,n)=wwh(ix,jy,llev)*pinmconv(llev)
      ww(ix,jy,nz,n)=wwh(ix,jy,nwz)*pinmconv(nz)
      kmin=llev+1
      do iz=2,nz
        do kz=kmin,nwz
          if ((height(iz).gt.wzlev(kz-1)).and. &
           (height(iz).le.wzlev(kz))) then
            dz1=height(iz)-wzlev(kz-1)
            dz2=wzlev(kz)-height(iz)
            dz=dz1+dz2
            ww(ix,jy,iz,n)=(wwh(ix,jy,kz-1)*pinmconv(kz-1)*dz2 &
              +wwh(ix,jy,kz)*pinmconv(kz)*dz1)/dz
          endif
        enddo
      enddo

!if ((jy.eq.0).and.(ix.eq.0)) print*, 'in loop 6'
  ! Compute density gradients at intermediate levels
  !*************************************************

      drhodz(ix,jy,1,n)=(rho(ix,jy,2,n)-rho(ix,jy,1,n))/(height(2)-height(1))
      do kz=2,nz-1
        drhodz(ix,jy,kz,n)=(rho(ix,jy,kz+1,n)-rho(ix,jy,kz-1,n))/ &
          (height(kz+1)-height(kz-1))
      enddo
      drhodz(ix,jy,nz,n)=drhodz(ix,jy,nz-1,n)

    !if ((jy.eq.0).and.(ix.eq.0))


    enddo
  enddo
!$OMP END DO

  !****************************************************************
  ! Compute slope of eta levels in windward direction and resulting
  ! vertical wind correction
  !****************************************************************


!$OMP DO
  do jy=1,ny-2
    cosf=cos((real(jy)*dy+ylat0)*pi180)

    do ix=1,nx-2
   ! print*, '1] slope of eta levels jy, ix=',jy, ix

   !if ((jy.le.2).and.(ix.eq.1)) print*, 'in eta loop 1'

  ! NCEP version: find first level above ground
      llev = 0
      do i=1,nuvz
       if (ps(ix,jy,1,n).lt.akz(i)) llev=i
      end do
      llev = llev+1

!if ((jy.le.2).and.(ix.eq.1)) print*, 'in eta loop 2'


      if (llev.gt.nuvz-2) llev = nuvz-2
  !     if (llev.eq.nuvz-2) write(*,*) 'verttransform
  !    +WARNING: LLEV eq NUZV-2'
  ! NCEP version

      kmin=llev+1

      !if ((jy.le.3).and.(ix.gt.250)) print*, 'in eta loop 3'

      do iz=2,nz-1

        ui=uu(ix,jy,iz,n)*dxconst/cosf
        vi=vv(ix,jy,iz,n)*dyconst

        klp=nz+1
        do kz=kmin,nz

          if ((height(iz).gt.uvwzlev(ix,jy,kz-1)).and. &
          (height(iz).le.uvwzlev(ix,jy,kz))) then
            dz1=height(iz)-uvwzlev(ix,jy,kz-1)
            dz2=uvwzlev(ix,jy,kz)-height(iz)
            dz=dz1+dz2
            kl=kz-1
            klp=kz
            exit
          endif

        enddo

        if (klp.eq.nz+1) then
          klp=nz
          kl=nz-1

         ! real,dimension(0:nxmax-1,0:nymax-1,nzmax) :: uvwzlev

          dz1=uvwzlev(ix,jy,klp)-uvwzlev(ix,jy,kl)


          dz2=0.
        endif

        ix1=ix-1
        jy1=jy-1
        ixp=ix+1
        jyp=jy+1

        dzdx1=(uvwzlev(ixp,jy,kl)-uvwzlev(ix1,jy,kl))*0.5
        dzdx2=(uvwzlev(ixp,jy,klp)-uvwzlev(ix1,jy,klp))*0.5
        dzdx=(dzdx1*dz2+dzdx2*dz1)/dz

        dzdy1=(uvwzlev(ix,jyp,kl)-uvwzlev(ix,jy1,kl))*0.5
        dzdy2=(uvwzlev(ix,jyp,klp)-uvwzlev(ix,jy1,klp))*0.5
        dzdy=(dzdy1*dz2+dzdy2*dz1)/dz

        ww(ix,jy,iz,n)=ww(ix,jy,iz,n)+(dzdx*ui+dzdy*vi)

      enddo ! z

      !if ((jy.le.3).and.(ix.eq.1)) print*, 'in eta loop end z jy=',jy


    enddo


  enddo
!$OMP END DO


  ! If north pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************

  if (nglobal) then
    do iz=1,nz
      do jy=int(switchnorthg)-2,nymin1
        ylat=ylat0+real(jy)*dy
        do ix=0,nxmin1
          xlon=xlon0+real(ix)*dx
          call cc2gll(northpolemap,ylat,xlon,uu(ix,jy,iz,n), &
               vv(ix,jy,iz,n),uupol(ix,jy,iz,n),vvpol(ix,jy,iz,n))
        enddo
      enddo
    enddo


    do iz=1,nz

      ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,nymin1,iz,n)**2+vv(nx/2-1,nymin1,iz,n)**2)
      if (vv(nx/2-1,nymin1,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,nymin1,iz,n)/vv(nx/2-1,nymin1,iz,n))-xlonr
      elseif (vv(nx/2-1,nymin1,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,nymin1,iz,n)/vv(nx/2-1,nymin1,iz,n))-xlonr
      else
        ddpol=pi*0.5-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.*pi+ddpol
      if(ddpol.gt.2.*pi) ddpol=ddpol-2.*pi

      ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.
      xlonr=xlon*pi/180.
      ylat=90.
      uuaux=-ffpol*sin(xlonr+ddpol)
      vvaux=-ffpol*cos(xlonr+ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux,vvpolaux)
      jy=nymin1
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      enddo

    enddo


    ! Fix: Set W (vertical wind) at pole to the zonally averaged W of the next
    ! equator-ward parallel


    do iz=1,nz
      wdummy=0.
      jy=ny-2
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      enddo
      wdummy=wdummy/real(nx)
      jy=nymin1
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      enddo
    enddo

  endif


  ! If south pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************

  if (sglobal) then
    do iz=1,nz
      do jy=0,int(switchsouthg)+3
        ylat=ylat0+real(jy)*dy
        do ix=0,nxmin1
          xlon=xlon0+real(ix)*dx
          call cc2gll(southpolemap,ylat,xlon,uu(ix,jy,iz,n), &
            vv(ix,jy,iz,n),uupol(ix,jy,iz,n),vvpol(ix,jy,iz,n))
        enddo
      enddo
    enddo

    do iz=1,nz

      ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,0,iz,n)**2+vv(nx/2-1,0,iz,n)**2)
      if (vv(nx/2-1,0,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,0,iz,n)/vv(nx/2-1,0,iz,n))+xlonr
      elseif (vv(nx/2-1,0,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,0,iz,n)/vv(nx/2-1,0,iz,n))-xlonr
      else
        ddpol=pi*0.5-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.*pi+ddpol
      if(ddpol.gt.2.*pi) ddpol=ddpol-2.*pi

      ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.
      xlonr=xlon*pi/180.
      ylat=-90.
      uuaux=+ffpol*sin(xlonr-ddpol)
      vvaux=-ffpol*cos(xlonr-ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux,vvpolaux)

      jy=0
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      enddo

    enddo


    ! Fix: Set W at pole to the zonally averaged W of the next equator-
    ! ward parallel of latitude

    do iz=1,nz
      wdummy=0.
      jy=1
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=0
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      enddo
    enddo

  endif


  ! PS, AT: for v10.5, we add back the quick fix to interpolate clouds in
  !   interpol_rain.f90 developed by PS for v8 and extend it to using
  !   cloud water fields

  !*******************************************************************************
  if (lcw) then ! identify clouds based on cloud water content
  !*******************************************************************************

    write(*,*) 'Global NCEP fields: using cloud water'

    ctwc(:,:,n)=0. ! initialise cloud total water content

    ! If water/ice are read separately into clwc and ciwc, store sum in clwc
    if (.not. lcwsum) clwc(:,:,:,n) = clwc(:,:,:,n) + ciwc(:,:,:,n)

    do kz = 1,nz-1
      do jy=0,nymin1
        do ix=0,nxmin1
          if (kz .eq. 1) then
            icloudbot(ix,jy,n) = icmv
!!          icloudtop=icmv ! this is just a local variable
!           we will use icloudtop as workspace for cloud top
          endif

          ! vertically integrate cloud water and determine cloud bottom, top
          ! cloud water per cell in kg / m2
          ! calculate cloud water mass per area: kgCW/kgAIR * kgAIR/m3 * m = kgCW/m2

          clw = clwc(ix,jy,kz,n)*rho(ix,jy,kz,n)*(height(kz+1)-height(kz))
          ! Add this layer to column cloud water [m3/m3]
          ctwc(ix,jy,n) = ctwc(ix,jy,n)+clw ! kg / m2 in column

          if (clw .gt. 0.) then ! cloud layer - maybe use threshold?
            if (icloudbot(ix,jy,n) .eq. icmv) &
              icloudbot(ix,jy,n) = nint(height(kz))
            icloudtop(ix,jy,n) = nint(height(kz))
          endif

          if (kz .eq. nz-1) then ! top level
            ! memorise icloudtop
            icloudtop_old = icloudtop(ix,jy,n)
            ! limit cloud top to 19 km:
            if (icloudtop(ix,jy,n) .gt. 19000) icloudtop(ix,jy,n) = 19000
            if (icloudbot(ix,jy,n) .eq. icmv) then
              icloudtop(ix,jy,n) = icmv
            endif

           ! PS  get rid of too thin clouds
            if (icloudtop(ix,jy,n) .lt. 50) then
              icloudbot(ix,jy,n)=icmv
              icloudtop(ix,jy,n)=icmv
            endif

            ! PS implement a rough fix for badly represented convection
            ! PS is based on looking at a limited set of comparison data
            lsp=  sum(   lsprec(ix,jy,1,:,n) )
            convp=sum( convprec(ix,jy,1,:,n) )
            prec=lsp+convp
            if (lsp.gt.convp) then !  prectype='lsp'
              lconvectprec = .false.
            else ! prectype='cp '
              lconvectprec = .true.
            endif
            if (lconvectprec .and. prec .gt. precmin .and.  &
              (icloudtop_old .lt. 6000 .or. icloudbot(ix,jy,n) .gt. 3000) ) then
              if (convp .lt. 0.1) then
                icloudbot(ix,jy,n) = 500
                icloudtop(ix,jy,n) =         8000
              else
                icloudbot(ix,jy,n) = 0
                icloudtop(ix,jy,n) =      10000
              endif
            endif
          endif ! end top level

        enddo ! ix loop
      enddo ! jy loop
    enddo ! kz loop


!**************************************************************************
  else       ! identify clouds using relative humidity
!**************************************************************************
!   clouds occur where rh>90% (using rh_ice for T<-20 deg C)

    write(*,*) 'NCEP fields: using relative humidity for cloud &
        &identification'
    do kz = 1,nz-1
      do jy=0,nymin1
        do ix=0,nxmin1

!PS       note that original by Sabine Eckhart was 80%
!PS       however, for T<-20 C we consider saturation over ice
!PS       so I think 90% should be enough
          if (kz .eq. 1) then
            icloudbot(ix,jy,n) = icmv
!!          icloudtop=icmv ! this is just a local variable
!           we will use icloudtop as workspace for cloud top
          endif
!98        do kz=1,nz
          pressure=rho(ix,jy,kz,n)*r_air*tt(ix,jy,kz,n)
          rh=qv(ix,jy,kz,n)/f_qvsat(pressure,tt(ix,jy,kz,n))
!PS            if (prec.gt.0.01) print*,'relhum',prec,kz,rh,height(kz)
          if (rh .ge. rhmin) then
            if (icloudbot(ix,jy,n) .eq. icmv) then
              icloudbot(ix,jy,n)=nint(height(kz))! use int to save memory
            endif
            icloudtop(ix,jy,n)=nint(height(kz)) ! use int to save memory
          endif
!          enddo

!PS/AT 2021: in this version, we skip the iteration with smaller rhmin
!PS try to get a cloud thicker than 50 m
!PS if there is at least .01 mm/h  - changed to 0.002 and put into
!PS parameter precpmin
!          if ((icloudbot(ix,jy,n) .eq. icmv .or. &
!            icloudtop-icloudbot(ix,jy,n) .lt. 50) .and. prec .gt. precmin) then
!              rhmin = rhmin - 0.05
!              if (rhmin .ge. 0.30) goto 98 ! give up for <= 25% rel.hum.
!          endif
!          if (icloudtop .ne. icmv) then
!            icloudtop(ix,jy,n) = icloudtop-icloudbot(ix,jy,n)
!          else
!            icloudtop(ix,jy,n) = icmv
!          endif

          if (kz .eq. nz-1) then ! top level

            ! memorise icloudtop
            icloudtop_old = icloudtop(ix,jy,n)
            ! limit cloud top to 19 km:
            if (icloudtop(ix,jy,n) .gt. 19000) icloudtop(ix,jy,n) = 19000
            if (icloudbot(ix,jy,n) .ne. icmv) then
              icloudtop(ix,jy,n) = icloudtop(ix,jy,n)-icloudbot(ix,jy,n)
            else
              icloudtop(ix,jy,n) = icmv
            endif

           ! PS  get rid of too thin clouds
            if (icloudtop(ix,jy,n) .lt. 50) then
              icloudbot(ix,jy,n)=icmv
              icloudtop(ix,jy,n)=icmv
            endif

            ! PS implement a rough fix for badly represented convection
            ! PS is based on looking at a limited set of comparison data
            lsp=  sum(   lsprec(ix,jy,1,:,n) )
            convp=sum( convprec(ix,jy,1,:,n) )
            prec=lsp+convp
            if (lsp.gt.convp) then !  prectype='lsp'
              lconvectprec = .false.
            else ! prectype='cp '
              lconvectprec = .true.
            endif
            if (lconvectprec .and. prec .gt. precmin .and.  &
              (icloudtop_old .lt. 6000 .or. icloudbot(ix,jy,n) .gt. 3000) ) then
              if (convp .lt. 0.1) then
                icloudbot(ix,jy,n) = 500
                icloudtop(ix,jy,n) = 8000
              endif
            else
              icloudbot(ix,jy,n) = 0
              icloudtop(ix,jy,n) = 10000
            endif

          endif ! end top level

        enddo ! ix loop
      enddo ! jy loop
    enddo ! kz loop

!**************************************************************************
  endif ! lcw true/false
!**************************************************************************

end subroutine verttransform_gfs

subroutine verttransform_ecmwf_heights(nxlim,nylim, &
  tt2_tmp,td2_tmp,ps_tmp,qvh_tmp,tth_tmp,prsh_tmp, &
  rhoh_tmp,pinmconv,uvzlev,wzlev)

  implicit none

  integer, intent(in) :: nxlim,nylim
  real,intent(in),dimension(0:nxlim,0:nylim) :: tt2_tmp,td2_tmp,ps_tmp
  real,intent(in),dimension(0:nxlim,0:nylim,nuvzmax) :: qvh_tmp,tth_tmp
  real,intent(out),dimension(0:nxlim,0:nylim,nuvzmax) :: rhoh_tmp,prsh_tmp
  real,intent(out),dimension(0:nxlim,0:nylim,nzmax) :: pinmconv
  real,intent(out),dimension(0:nxlim,0:nylim,nuvzmax) :: uvzlev,wzlev
  !real,dimension(0:nxlim,0:nylim) :: tvold,pold,pint,tv
  real :: tvold,pold,pint,tv
  real,parameter :: const=r_air/ga
  integer :: ix,jy,kz

  ! Loop over the whole grid
  !*************************
!$OMP PARALLEL PRIVATE(jy,ix,pint,tv,tvold,pold,kz)
!$OMP DO
  do jy=0,nylim
    do ix=0,nxlim
      tvold=tt2_tmp(ix,jy)*(1.+0.378*ew(td2_tmp(ix,jy),ps_tmp(ix,jy))/ &
           ps_tmp(ix,jy))
      pold=ps_tmp(ix,jy)
      uvzlev(ix,jy,1)=0.
      wzlev(ix,jy,1)=0.
      rhoh_tmp(ix,jy,1)=pold/(r_air*tvold)
      prsh_tmp(ix,jy,1)=ps_tmp(ix,jy)

      do kz=2,nuvz
        pint=akz(kz)+bkz(kz)*ps_tmp(ix,jy)
        prsh_tmp(ix,jy,kz)=pint
        tv=tth_tmp(ix,jy,kz)*(1.+0.608*qvh_tmp(ix,jy,kz))
        rhoh_tmp(ix,jy,kz)=pint/(r_air*tv)

        if (abs(tv-tvold).gt.0.2) then
          uvzlev(ix,jy,kz)=uvzlev(ix,jy,kz-1)+const* &
               log(pold/pint)*(tv-tvold)/log(tv/tvold)
        else
          uvzlev(ix,jy,kz)=uvzlev(ix,jy,kz-1)+const* &
               log(pold/pint)*tv
        endif

        tvold=tv
        pold=pint

      end do

      do kz=2,nwz-1
        wzlev(ix,jy,kz)=(uvzlev(ix,jy,kz+1)+uvzlev(ix,jy,kz))*0.5
      end do
      wzlev(ix,jy,nwz)=wzlev(ix,jy,nwz-1)+ &
           uvzlev(ix,jy,nuvz)-uvzlev(ix,jy,nuvz-1)


      pinmconv(ix,jy,1)=(uvzlev(ix,jy,2))/ &
           ((aknew(2)+bknew(2)*ps_tmp(ix,jy))- &
           (aknew(1)+bknew(1)*ps_tmp(ix,jy)))
      do kz=2,nz-1
        pinmconv(ix,jy,kz)=(uvzlev(ix,jy,kz+1)-uvzlev(ix,jy,kz-1))/ &
             ((aknew(kz+1)+bknew(kz+1)*ps_tmp(ix,jy))- &
             (aknew(kz-1)+bknew(kz-1)*ps_tmp(ix,jy)))
      end do
      pinmconv(ix,jy,nz)=(uvzlev(ix,jy,nz)-uvzlev(ix,jy,nz-1))/ &
           ((aknew(nz)+bknew(nz)*ps_tmp(ix,jy))- &
           (aknew(nz-1)+bknew(nz-1)*ps_tmp(ix,jy)))
    end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  ! pold(:,:)=ps_tmp(:,:)
  ! uvzlev(:,:,1)=0.
  ! wzlev(:,:,1)=0.
  ! rhoh_tmp(:,:,1)=pold(:,:)/(r_air*tvold(:,:))
  ! prsh_tmp(:,:,1)=ps_tmp(:,:)

  ! Compute heights of eta levels
  !******************************

  ! do kz=2,nuvz
  !   pint(:,:)=akz(kz)+bkz(kz)*ps_tmp(:,:)
  !   prsh_tmp(:,:,kz)=pint(:,:)
  !   tv(:,:)=tth_tmp(:,:,kz)*(1.+0.608*qvh_tmp(:,:,kz))
  !   rhoh_tmp(:,:,kz)=pint(:,:)/(r_air*tv(:,:))

  !   where (abs(tv(:,:)-tvold(:,:)).gt.0.2)
  !     uvzlev(:,:,kz)=uvzlev(:,:,kz-1)+const*&
  !          &log(pold(:,:)/pint(:,:))* &
  !          (tv(:,:)-tvold(:,:))/&
  !          &log(tv(:,:)/tvold(:,:))
  !   elsewhere
  !     uvzlev(:,:,kz)=uvzlev(:,:,kz-1)+const*&
  !          &log(pold(:,:)/pint(:,:))*tv(:,:)
  !   endwhere

  !   tvold(:,:)=tv(:,:)
  !   pold(:,:)=pint(:,:)

  ! end do

  ! do kz=2,nwz-1
  !   wzlev(:,:,kz)=(uvzlev(:,:,kz+1)+uvzlev(:,:,kz))/2.
  ! end do
  ! wzlev(:,:,nwz)=wzlev(:,:,nwz-1)+ &
  !      uvzlev(:,:,nuvz)-uvzlev(:,:,nuvz-1)


  ! pinmconv(:,:,1)=(uvzlev(:,:,2))/ &
  !      ((aknew(2)+bknew(2)*ps_tmp(:,:))- &
  !      (aknew(1)+bknew(1)*ps_tmp(:,:)))
  ! do kz=2,nz-1
  !   pinmconv(:,:,kz)=(uvzlev(:,:,kz+1)-uvzlev(:,:,kz-1))/ &
  !        ((aknew(kz+1)+bknew(kz+1)*ps_tmp(:,:))- &
  !        (aknew(kz-1)+bknew(kz-1)*ps_tmp(:,:)))
  ! end do
  ! pinmconv(:,:,nz)=(uvzlev(:,:,nz)-uvzlev(:,:,nz-1))/ &
  !      ((aknew(nz)+bknew(nz)*ps_tmp(:,:))- &
  !      (aknew(nz-1)+bknew(nz-1)*ps_tmp(:,:)))
end subroutine verttransform_ecmwf_heights


end module verttransform_mod
