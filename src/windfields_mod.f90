! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

!*****************************************************************************
!                                                                            *
! This module stores all meteorological input data and contains routines     *
! reading and allocating this data                                           *
!                                                                            *
! L. Bakels 2023: Dynamical allocation of all windfields                     *
!                                                                            *
!*****************************************************************************

module windfields_mod
   use par_mod
   use com_mod
   use point_mod
   use pbl_profile_mod

   implicit none

   !******************************************************************************
   ! Variables associated with the ECMWF meteorological input data ("wind fields")
   !******************************************************************************
   !
   logical :: debug = .true.

   integer :: &
      numbwf                ! actual number of wind fields

   integer, allocatable, dimension(:) :: &
      wftime         ! times relative to beginning time of wind fields [s]

   character(len=255), allocatable, dimension(:) :: &
      wfname        ! file names of wind fields

   !Windfield parameters
   !********************
   integer :: nxmax, nymax, nuvzmax, nwzmax, nzmax !Size of windfield

   ! Fixed fields, unchangeable with time
   !*************************************
   real, allocatable, dimension(:, :) :: &
      oro, & ! orography of the ECMWF model
      excessoro, & ! excess orography mother domain
      lsm                                 ! land sea mask of the ECMWF model

   ! 3d fields
   !**********
   real, allocatable, dimension(:, :, :, :) :: &
      uu, vv, ww, & ! wind components in x,y and z direction [m/s]
      uupol, vvpol, & ! wind components in polar stereographic projection [m/s]
      tt, tth, & ! temperature data on internal and half model levels [K]
      qv, qvh, & ! specific humidity data on internal and half model levels (eta if 'ETA')
      pv, & ! potential vorticity
      rho, & ! air density [kg/m3]
      drhodz, & ! vertical air density gradient [kg/m2]
      pplev, & ! Pressure on half model levels
      prs, & ! air pressure RLT
      rho_dry                                 ! dry air density RLT Only printed out in binary mode???

   ! Cloud properties
   !*****************************************
   real, allocatable, dimension(:, :, :, :) :: &
      clwc, & ! liquid   [kg/kg] ZHG
      ciwc, & ! ice      [kg/kg] ZHG
      clwch, & ! original eta level liquid [kg/kg] ZHG
      ciwch                                   ! original eta level ice [kg/kg] ZHG
   real, allocatable, dimension(:, :, :) :: &
      ctwc                                    ! ESO: =icloud_stats(:,:,4,:) total cloud water content
   integer, allocatable, dimension(:, :, :) :: & ! new scavenging AT 2021
      icloudbot, & ! cloud bottom height [m/eta]
      icloudtop                               ! cloud top [m/eta]

   ! 2d fields
   !**********
   real, allocatable, dimension(:, :, :, :) :: &
      ps, & ! surface pressure
      sd, & ! snow depth
      msl, & ! mean sea level pressure
      tcc, & ! total cloud cover
      u10, & ! 10 meter u
      v10, & ! 10 meter v
      tt2, & ! 2 meter temperature
      td2, & ! 2 meter dew point
      sshf, & ! surface sensible heat flux
      ssr, & ! surface solar radiation
      sfcstress, & ! surface stress
      ustar, & ! friction velocity [m/s]
      wstar, & ! convective velocity scale [m/s]
      hmix, & ! mixing height [m]
      tropopause, & ! altitude of thermal tropopause [m]
      oli                                   ! inverse Obukhov length (1/L) [m]

   ! 2d fields
   !**********
   real, allocatable, dimension(:, :, :, :, :) :: & ! newWetDepoScheme, extra precip dimension AT 2021
      lsprec, & ! large scale total precipitation [mm/h]
      convprec                                  ! convective precipitation [mm/h]

   integer :: metdata_format  ! storing the input data type (ECMWF/NCEP)

   !****************************************************************************
   ! Variables defining actual size and geographical location of the wind fields
   !****************************************************************************

   integer :: &
      nx, ny, nz, & ! actual dimensions of wind fields in x,y and z direction, respectively
      nxmin1, & ! nx-1
      nymin1, & ! ny-1
      nxfield, & ! same as nx for limited area fields, but for global fields nx=nxfield+1
      nuvz, nwz, & ! vertical dimension of original ECMWF data (u,v components/ w components(staggered grid))
      nmixz, & ! number of levels up to maximum PBL height (3500 m)
      nlev_ec      ! number of levels ECMWF model
   real :: &
      dxconst, & ! auxiliary variables for utransform
      dyconst      ! auxiliary variables for vtransform
   integer :: nconvlevmax !maximum number of levels for convection
   integer :: na ! parameter used in Emanuel's convect subroutine

   !*************************************************
   ! Variables used for vertical model discretization
   !*************************************************
   real, allocatable, dimension(:) :: &
      height, & ! heights of all levels [m]
      wheight, & ! model level heights [m]
      uvheight, & ! half-model level heights [m]
      akm, bkm, & ! coefficients which regulate vertical discretization of ecmwf model levels
      akz, bkz, & ! model discretization coefficients at the centre of the layers
      aknew, bknew                      ! model discretization coefficients at the interpolated levels

contains

   subroutine detectformat

      !*****************************************************************************
      !                                                                            *
      !   This routine reads the 1st file with windfields to determine             *
      !   the format.                                                              *
      !                                                                            *
      !     Authors: M. Harustak                                                   *
      !                                                                            *
      !     6 May 2015                                                             *
      !                                                                            *
      !   Unified ECMWF and GFS builds                                             *
      !   Marian Harustak, 12.5.2017                                               *
      !     - Added routine to FP10 Flexpart distribution                          *
      !*****************************************************************************
      !                                                                            *
      ! Variables:                                                                 *
      ! fname                file name of file to check                            *
      !                                                                            *
      !*****************************************************************************

      use par_mod
      use com_mod
      use class_gribfile_mod

      implicit none

      character(len=255) :: filename

      ! If no file is available
      if (numbwf .le. 0) then
         print *, 'No wind file available'
         metdata_format = GRIBFILE_CENTRE_UNKNOWN
         return
      end if

      ! construct filename
      filename = path(3) (1:length(3))//trim(wfname(1))

      ! get format
      metdata_format = gribfile_centre(TRIM(filename))
   end subroutine detectformat

   subroutine gridcheck_gfs

      !**********************************************************************
      !                                                                     *
      !             FLEXPART MODEL SUBROUTINE GRIDCHECK                     *
      !                                                                     *
      !**********************************************************************
      !                                                                     *
      !             AUTHOR:      G. WOTAWA                                  *
      !             DATE:        1997-08-06                                 *
      !             LAST UPDATE: 1997-10-10                                 *
      !                                                                     *
      !             Update:      1999-02-08, global fields allowed, A. Stohl*
      !             CHANGE: 17/11/2005, Caroline Forster, GFS data          *
      !             CHANGE: 11/01/2008, Harald Sodemann, GRIB1/2 input with *
      !                                 ECMWF grib_api                      *
      !             CHANGE: 03/12/2008, Harald Sodemann, update to f90 with *
      !                                 ECMWF grib_api                      *
      !                                                                     *
      !   Unified ECMWF and GFS builds                                      *
      !   Marian Harustak, 12.5.2017                                        *
      !     - Renamed routine from gridcheck to gridcheck_gfs               *
      !                                                                     *
      !                                                                     *
      !  Anne Tipka, Petra Seibert 2021-02: implement new interpolation     *
      !    for precipitation according to #295 using 2 additional fields    *
      !                                                                     *
      !**********************************************************************
      !                                                                     *
      ! DESCRIPTION:                                                        *
      !                                                                     *
      ! THIS SUBROUTINE DETERMINES THE GRID SPECIFICATIONS (LOWER LEFT      *
      ! LONGITUDE, LOWER LEFT LATITUDE, NUMBER OF GRID POINTS, GRID DIST-   *
      ! ANCE AND VERTICAL DISCRETIZATION OF THE ECMWF MODEL) FROM THE       *
      ! GRIB HEADER OF THE FIRST INPUT FILE. THE CONSISTANCY (NO CHANGES    *
      ! WITHIN ONE FLEXPART RUN) IS CHECKED IN THE ROUTINE "READWIND" AT    *
      ! ANY CALL.                                                           *
      !                                                                     *
      ! XLON0                geographical longitude of lower left gridpoint *
      ! YLAT0                geographical latitude of lower left gridpoint  *
      ! NX                   number of grid points x-direction              *
      ! NY                   number of grid points y-direction              *
      ! DX                   grid distance x-direction                      *
      ! DY                   grid distance y-direction                      *
      ! NUVZ                 number of grid points for horizontal wind      *
      !                      components in z direction                      *
      ! NWZ                  number of grid points for vertical wind        *
      !                      component in z direction                       *
      ! sizesouth, sizenorth give the map scale (i.e. number of virtual grid*
      !                      points of the polar stereographic grid):       *
      !                      used to check the CFL criterion                *
      ! UVHEIGHT(1)-         heights of gridpoints where u and v are        *
      ! UVHEIGHT(NUVZ)       given                                          *
      ! WHEIGHT(1)-          heights of gridpoints where w is given         *
      ! WHEIGHT(NWZ)                                                        *
      !                                                                     *
      !**********************************************************************

      use grib_api
      use cmapf_mod, only: stlmbr, stcm2p

      implicit none

      !HSO  parameters for grib_api
      integer :: ifile
      integer :: iret
      integer :: igrib, stat, size1
      real(kind=4) :: xaux1, xaux2, yaux1, yaux2
      real(kind=8) :: xaux1in, xaux2in, yaux1in, yaux2in
      integer :: gribVer, parCat, parNum, typSfc, valSurf, discipl
      !HSO  end
      integer :: ix, jy, i, ifn, ifield, j, k, iumax, iwmax, numskip
      real :: sizesouth, sizenorth, xauxa
      real :: xsec18 !ip
      real, allocatable, dimension(:) :: akm_usort, pres, tmppres
      real, parameter :: eps = 0.0001

      ! NCEP GFS
      real :: help

      integer :: i179, i180, i181

      ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

      integer :: isec1(8), isec2(3)
      real(kind=4), allocatable, dimension(:) :: zsec4
      character(len=1) :: opt

      !HSO  grib api error messages
      character(len=24) :: gribErrorMsg = 'Error reading grib file'
      character(len=20) :: gribFunction = 'gridcheckwind_gfs'

      iumax = 0
      iwmax = 0

      if (ideltas .gt. 0) then
         ifn = 1
      else
         ifn = numbwf
      end if
      !
      ! OPENING OF DATA FILE (GRIB CODE)
      !
5     call grib_open_file(ifile, path(3) (1:length(3)) &
                          //trim(wfname(ifn)), 'r', iret)
      if (iret .ne. GRIB_SUCCESS) then
         goto 999   ! ERROR DETECTED
      end if
      !turn on support for multi fields messages
      call grib_multi_support_on

      ifield = 0
      do
         ifield = ifield + 1
         !
         ! GET NEXT FIELDS
         !
         call grib_new_from_file(ifile, igrib, iret)
         if (iret .eq. GRIB_END_OF_FILE) then
            exit    ! EOF DETECTED
         elseif (iret .ne. GRIB_SUCCESS) then
            goto 999   ! ERROR DETECTED
         end if

         !first see if we read GRIB1 or GRIB2
         call grib_get_int(igrib, 'editionNumber', gribVer, iret)
         call grib_check(iret, gribFunction, gribErrorMsg)

         if (gribVer .eq. 1) then ! GRIB Edition 1

            !read the grib1 identifiers
            call grib_get_int(igrib, 'indicatorOfParameter', isec1(6), iret)
            call grib_check(iret, gribFunction, gribErrorMsg)
            call grib_get_int(igrib, 'indicatorOfTypeOfLevel', isec1(7), iret)
            call grib_check(iret, gribFunction, gribErrorMsg)
            call grib_get_int(igrib, 'level', isec1(8), iret)
            call grib_check(iret, gribFunction, gribErrorMsg)

            xsec18 = real(isec1(8))

         else ! GRIB Edition 2

            !read the grib2 identifiers
            call grib_get_int(igrib, 'discipline', discipl, iret)
            call grib_check(iret, gribFunction, gribErrorMsg)
            call grib_get_int(igrib, 'parameterCategory', parCat, iret)
            call grib_check(iret, gribFunction, gribErrorMsg)
            call grib_get_int(igrib, 'parameterNumber', parNum, iret)
            call grib_check(iret, gribFunction, gribErrorMsg)
            call grib_get_int(igrib, 'typeOfFirstFixedSurface', typSfc, iret)
            call grib_check(iret, gribFunction, gribErrorMsg)
            call grib_get_int(igrib, 'scaledValueOfFirstFixedSurface', &
                              valSurf, iret)
            call grib_check(iret, gribFunction, gribErrorMsg)

            !convert to grib1 identifiers
            isec1(6) = -1
            isec1(7) = -1
            isec1(8) = -1
            xsec18 = -1.0

            if ((parCat .eq. 2) .and. (parNum .eq. 2) .and. (typSfc .eq. 100)) then ! U
               isec1(6) = 33          ! indicatorOfParameter
               isec1(7) = 100         ! indicatorOfTypeOfLevel
               isec1(8) = valSurf/100 ! level, convert to hPa
               xsec18 = valSurf/100.0 ! level, convert to hPa

               ! fixgfs11
               call grib_get_size(igrib, 'values', size1, iret)
               allocate (zsec4(size1), stat=stat)
               if (stat .ne. 0) error stop "Could not allocate zsec4"
               call grib_get_real4_array(igrib, 'values', zsec4, iret)
               call grib_check(iret, gribFunction, gribErrorMsg)
               !PRINT*,zsec4(1:15)
               !stop    'MIP2 in gridcheck'
               deallocate (zsec4)

            elseif ((parCat .eq. 3) .and. (parNum .eq. 5) .and. (typSfc .eq. 1)) then ! TOPO
               isec1(6) = 7           ! indicatorOfParameter
               isec1(7) = 1           ! indicatorOfTypeOfLevel
               isec1(8) = 0
               xsec18 = real(0)

            elseif ((parCat .eq. 0) .and. (parNum .eq. 0) .and. (typSfc .eq. 1) &
                    .and. (discipl .eq. 2)) then ! LSM
               isec1(6) = 81          ! indicatorOfParameter
               isec1(7) = 1           ! indicatorOfTypeOfLevel
               isec1(8) = 0
               xsec18 = real(0)
            end if

         end if ! gribVer

         if (ifield .eq. 1) then

            !get the required fields from section 2
            !store compatible to gribex input
            call grib_get_int(igrib, 'numberOfPointsAlongAParallel', &
                              isec2(2), iret)
            call grib_check(iret, gribFunction, gribErrorMsg)
            call grib_get_int(igrib, 'numberOfPointsAlongAMeridian', &
                              isec2(3), iret)
            call grib_check(iret, gribFunction, gribErrorMsg)
            call grib_get_real8(igrib, 'longitudeOfFirstGridPointInDegrees', &
                                xaux1in, iret)
            call grib_check(iret, gribFunction, gribErrorMsg)
            call grib_get_real8(igrib, 'longitudeOfLastGridPointInDegrees', &
                                xaux2in, iret)
            call grib_check(iret, gribFunction, gribErrorMsg)
            call grib_get_real8(igrib, 'latitudeOfLastGridPointInDegrees', &
                                yaux1in, iret)
            call grib_check(iret, gribFunction, gribErrorMsg)
            call grib_get_real8(igrib, 'latitudeOfFirstGridPointInDegrees', &
                                yaux2in, iret)
            call grib_check(iret, gribFunction, gribErrorMsg)

            ! Fix for flexpart.eu ticket #48
            if (xaux2in .lt. 0) xaux2in = 359.0

            xaux1 = real(xaux1in)
            xaux2 = real(xaux2in)
            yaux1 = real(yaux1in)
            yaux2 = real(yaux2in)

            nxfield = isec2(2)
            ny = isec2(3)
            if ((abs(xaux1) .lt. eps) .and. (xaux2 .ge. 359)) then ! NCEP DATA FROM 0 TO

               ! fixgfs11
               ! xaux1=-179.0                             ! 359 DEG EAST ->
               ! xaux2=-179.0+360.-360./real(nxfield)    ! TRANSFORMED TO -179
               ! reset to working v10 settings
               xaux1 = -180.0                             ! 359 DEG EAST ->
               xaux2 = -180.0 + 360.-360./real(nxfield)    ! TRANSFORMED TO -179
            end if                                      ! TO 180 DEG EAST

            if (xaux1 .gt. 180) xaux1 = xaux1 - 360.0
            if (xaux2 .gt. 180) xaux2 = xaux2 - 360.0
            if (xaux1 .lt. -180) xaux1 = xaux1 + 360.0
            if (xaux2 .lt. -180) xaux2 = xaux2 + 360.0
            if (xaux2 .lt. xaux1) xaux2 = xaux2 + 360.
            xlon0 = xaux1
            ylat0 = yaux1
            dx = (xaux2 - xaux1)/real(nxfield - 1)
            dy = (yaux2 - yaux1)/real(ny - 1)
            dxconst = 180./(dx*r_earth*pi)
            dyconst = 180./(dy*r_earth*pi)
            !HSO end edits

            ! Check whether fields are global
            ! If they contain the poles, specify polar stereographic map
            ! projections using the stlmbr- and stcm2p-calls
            !***********************************************************

            xauxa = abs(xaux2 + dx - 360.-xaux1)
            if (xauxa .lt. 0.001) then
               nx = nxfield + 1                 ! field is cyclic
               xglobal = .true.
               if (abs(nxshift) .ge. nx) &
                  error stop 'nxshift in file par_mod is too large'
               xlon0 = xlon0 + real(nxshift)*dx
            else
               nx = nxfield
               xglobal = .false.
               if (nxshift .ne. 0) &
                  error stop 'nxshift (par_mod) must be zero for non-global domain'
            end if
            nxmin1 = nx - 1
            nymin1 = ny - 1
            if (xlon0 .gt. 180.) xlon0 = xlon0 - 360.
            xauxa = abs(yaux1 + 90.)
            if (xglobal .and. xauxa .lt. 0.001) then
               sglobal = .true.               ! field contains south pole
               ! Enhance the map scale by factor 3 (*2=6) compared to north-south
               ! map scale
               sizesouth = 6.*(switchsouth + 90.)/dy
               call stlmbr(southpolemap, -90., 0.)
               call stcm2p(southpolemap, 0., 0., switchsouth, 0., sizesouth, &
                           sizesouth, switchsouth, 180.)
               switchsouthg = (switchsouth - ylat0)/dy
            else
               sglobal = .false.
               switchsouthg = 999999.
            end if
            xauxa = abs(yaux2 - 90.)
            if (xglobal .and. xauxa .lt. 0.001) then
               nglobal = .true.               ! field contains north pole
               ! Enhance the map scale by factor 3 (*2=6) compared to north-south
               ! map scale
               sizenorth = 6.*(90.-switchnorth)/dy
               call stlmbr(northpolemap, 90., 0.)
               call stcm2p(northpolemap, 0., 0., switchnorth, 0., sizenorth, &
                           sizenorth, switchnorth, 180.)
               switchnorthg = (switchnorth - ylat0)/dy
            else
               nglobal = .false.
               switchnorthg = 999999.
            end if
            ! Set nxmax and nymax and allocate the fields for oro lsm and excessoro
            nxmax = nx
            nymax = ny
            call alloc_fixedfields

         end if ! ifield.eq.1

         if (nxshift .lt. 0) error stop 'nxshift (par_mod) must not be negative'
         if (nxshift .ge. nxfield) error stop 'nxshift (par_mod) too large'

         if (isec1(6) .ne. -1) then
            !  get the size and data of the values array
            !get the size and data of the values array
            call grib_get_size(igrib, 'values', size1, iret)
            call grib_check(iret, gribFunction, gribErrorMsg)
            allocate (zsec4(size1), stat=stat)
            if (stat .ne. 0) error stop "Could not allocate zsec4"
            call grib_get_real4_array(igrib, 'values', zsec4, iret)
            call grib_check(iret, gribFunction, gribErrorMsg)
         end if

         ! NCEP ISOBARIC LEVELS
         !*********************

         if ((isec1(6) .eq. 33) .and. (isec1(7) .eq. 100)) then ! check for U wind
            iumax = iumax + 1
            ! fixgfs11
            allocate (tmppres(iumax), stat=stat)
            if (stat .ne. 0) error stop "Could not allocate tmppres"
            if (iumax .gt. 1) tmppres(1:iumax - 1) = pres
            !pres(iumax)=real(isec1(8))*100.0
            call move_alloc(tmppres, pres)
            pres(iumax) = xsec18*100.0
            ! ip 30.1.24 fix vertical coordinate reading bug
         end if

         ! fixgfs11 TODO: finish cleanup
         i179 = nint(179./dx)
         if (dx .lt. 0.7) then
            i180 = nint(180./dx) + 1    ! 0.5 deg data
         else
            i180 = nint(179./dx) + 1    ! 1 deg data
         end if
         i181 = i180 + 1

         ! fixgfs11 -- revert to working v10.4 setting
         i180 = nint(180./dx)    ! 0.5 deg data
         i181 = i180
         i179 = i180

         ! NCEP TERRAIN
         !*************

         if ((isec1(6) .eq. 007) .and. (isec1(7) .eq. 001)) then

            ! IP 8/5/23
            do jy = 0, ny - 1
               do ix = 0, nxfield - 1
                  help = zsec4(nxfield*(ny - jy - 1) + ix + 1)
                  if (ix .le. i180) then
                     oro(i179 + ix, jy) = help
                     excessoro(i179 + ix, jy) = 0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
                  else
                     oro(ix - i181, jy) = help
                     excessoro(ix - i181, jy) = 0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
                  end if
               end do
            end do
         end if

         ! NCEP LAND SEA MASK
         !*******************

         if ((isec1(6) .eq. 081) .and. (isec1(7) .eq. 001)) then
            do jy = 0, ny - 1
               do ix = 0, nxfield - 1
                  help = zsec4(nxfield*(ny - jy - 1) + ix + 1)
                  if (ix .le. i180) then
                     lsm(i179 + ix, jy) = help
                  else
                     lsm(ix - i181, jy) = help
                  end if
               end do
            end do
         end if

         call grib_release(igrib)
         if (isec1(6) .ne. -1) deallocate (zsec4) !IP 28/11/23 fix to run GFS tests
      end do                      !! READ NEXT LEVEL OR PARAMETER
      !
      ! CLOSING OF INPUT DATA FILE
      !

      ! HSO
      call grib_close_file(ifile)
      ! HSO end edits

      nuvz = iumax
      nwz = iumax
      nlev_ec = iumax

      ! Allocate memory for windfields
      !*******************************
      nwzmax = nwz
      nuvzmax = nuvz
      nzmax = nuvz
      nconvlevmax = iumax
      na = nuvzmax
      call alloc_windfields

      if (nx .gt. nxmax) then
         write (*, *) 'FLEXPART error: Too many grid points in x direction.'
         write (*, *) 'Reduce resolution of wind fields.'
         write (*, *) 'Or change parameter settings in file par_mod.'
         write (*, *) nx, nxmax
         error stop
      end if

      if (ny .gt. nymax) then
         write (*, *) 'FLEXPART error: Too many grid points in y direction.'
         write (*, *) 'Reduce resolution of wind fields.'
         write (*, *) 'Or change parameter settings in file par_mod.'
         write (*, *) ny, nymax
         error stop
      end if

      if (nuvz .gt. nuvzmax) then
         write (*, *) 'FLEXPART error: Too many u,v grid points in z '// &
            'direction.'
         write (*, *) 'Reduce resolution of wind fields.'
         write (*, *) 'Or change parameter settings in file par_mod.'
         write (*, *) nuvz, nuvzmax
         error stop
      end if

      if (nwz .gt. nwzmax) then
         write (*, *) 'FLEXPART error: Too many w grid points in z '// &
            'direction.'
         write (*, *) 'Reduce resolution of wind fields.'
         write (*, *) 'Or change parameter settings in file par_mod.'
         write (*, *) nwz, nwzmax
         error stop
      end if

      ! If desired, shift all grids by nxshift grid cells
      !**************************************************

      if (xglobal) then
         call shift_field_0(oro, nxfield, ny)
         call shift_field_0(lsm, nxfield, ny)
         call shift_field_0(excessoro, nxfield, ny)
      end if

      ! Output of grid info
      !********************

      if (lroot) then
         write (*, *)
         write (*, *)
         write (*, '(a,2i7)') 'Vertical levels in NCEP data: ', &
            nuvz, nwz
         write (*, *)
         write (*, '(a)') 'Mother domain:'
         write (*, '(a,f10.2,a1,f10.2,a,f10.2)') '  Longitude range: ', &
            xlon0, ' to ', xlon0 + (nx - 1)*dx, '   Grid distance: ', dx
         write (*, '(a,f10.2,a1,f10.2,a,f10.2)') '  Latitude range : ', &
            ylat0, ' to ', ylat0 + (ny - 1)*dy, '   Grid distance: ', dy
         write (*, *)
      end if

      ! CALCULATE VERTICAL DISCRETIZATION OF ECMWF MODEL
      ! PARAMETER akm,bkm DESCRIBE THE HYBRID "ETA" COORDINATE SYSTEM

      allocate (akm_usort(nwzmax), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate akm_usort"
      numskip = nlev_ec - nuvz  ! number of ecmwf model layers not used
      ! by trajectory model
      do i = 1, nwz
         j = numskip + i
         k = nlev_ec + 1 + numskip + i
         akm_usort(nwz - i + 1) = pres(nwz - i + 1)
         bkm(nwz - i + 1) = 0.0
      end do

      !******************************
      ! change Sabine Eckhardt: akm should always be in descending order ... readwind adapted!
      !******************************
      do i = 1, nwz
         if (akm_usort(1) .gt. akm_usort(2)) then
            akm(i) = akm_usort(i)
         else
            akm(i) = akm_usort(nwz - i + 1)
         end if
      end do

      !
      ! CALCULATION OF AKZ, BKZ
      ! AKZ,BKZ: model discretization parameters at the center of each model
      !     layer
      !
      ! Assign the 10 m winds to an artificial model level with akz=0 and bkz=1.0,
      ! i.e. ground level
      !*****************************************************************************

      do i = 1, nuvz
         akz(i) = akm(i)
         bkz(i) = bkm(i)
      end do

      ! NOTE: In FLEXPART versions up to 4.0, the number of model levels was doubled
      ! upon the transformation to z levels. In order to save computer memory, this is
      ! not done anymore in the standard version. However, this option can still be
      ! switched on by replacing the following lines with those below, that are
      ! currently commented out. For this, similar changes are necessary in
      ! verttransform.f and verttranform_nest.f
      !*****************************************************************************

      nz = nuvz
      if (nz .gt. nzmax) error stop 'nzmax too small'
      do i = 1, nuvz
         aknew(i) = akz(i)
         bknew(i) = bkz(i)
      end do

      ! Switch on following lines to use doubled vertical resolution
      !*************************************************************
      !nz=nuvz+nwz-1
      !if (nz.gt.nzmax) stop 'nzmax too small'
      !do 100 i=1,nwz
      !  aknew(2*(i-1)+1)=akm(i)
      !00     bknew(2*(i-1)+1)=bkm(i)
      !do 110 i=2,nuvz
      !  aknew(2*(i-1))=akz(i)
      !10     bknew(2*(i-1))=bkz(i)
      ! End doubled vertical resolution
      deallocate (akm_usort, pres)
      return

999   write (*, *)
      write (*, *) ' ###########################################'// &
         '###### '
      write (*, *) '       TRAJECTORY MODEL SUBROUTINE GRIDCHECK:'
      write (*, *) ' CAN NOT OPEN INPUT DATA FILE '//wfname(ifn)
      write (*, *) ' ###########################################'// &
         '###### '
      write (*, *)
      write (*, '(a)') '!!! PLEASE INSERT A NEW CD-ROM AND   !!!'
      write (*, '(a)') '!!! PRESS ANY KEY TO CONTINUE...     !!!'
      write (*, '(a)') '!!! ...OR TERMINATE FLEXPART PRESSING!!!'
      write (*, '(a)') '!!! THE "X" KEY...                   !!!'
      write (*, *)
      read (*, '(a)') opt
      if (opt .eq. 'X') then
         error stop
      else
         goto 5
      end if

   end subroutine gridcheck_gfs

   subroutine readwind_gfs(indj, n, uuh, vvh, wwh)

      !***********************************************************************
      !                                                                      *
      !              TRAJECTORY MODEL SUBROUTINE READWIND                    *
      !                                                                      *
      !***********************************************************************
      !                                                                      *
      !              AUTHOR:      G. WOTAWA                                  *
      !              DATE:        1997-08-05                                 *
      !              LAST UPDATE: 2000-10-17, Andreas Stohl                  *
      !              CHANGE: 01/02/2001, Bernd C. Krueger, Variables tth and *
      !                      qvh (on eta coordinates) in common block        *
      !              CHANGE: 16/11/2005, Caroline Forster, GFS data          *
      !              CHANGE: 11/01/2008, Harald Sodemann, Input of GRIB1/2   *
      !                      data with the ECMWF grib_api library            *
      !              CHANGE: 03/12/2008, Harald Sodemann, update to f90 with *
      !                                  ECMWF grib_api                      *
      !                                                                      *
      !   Unified ECMWF and GFS builds                                       *
      !   Marian Harustak, 12.5.2017                                         *
      !     - Renamed routine from readwind to readwind_gfs                  *
      !                                                                      *
      !   Petra Seibert, Anne Tipka, 2021-02: implement new interpolation    *
      !   just catch numpf>1 and produce error msg, adjust rank of precip    *
      !   and correct some loops in bad order                                *
      !                                                                      *
      !                                                                      *
      !  Anne Tipka, Petra Seibert 2021-02: implement new interpolation      *
      !    for precipitation according to #295 using 2 additional fields     *
      !    change some double loops in wrong order to forall constructs      *
      !                                                                      *
      !***********************************************************************
      !*                                                                     *
      !* DESCRIPTION:                                                        *
      !*                                                                     *
      !* READING OF ECMWF METEOROLOGICAL FIELDS FROM INPUT DATA FILES. THE   *
      !* INPUT DATA FILES ARE EXPECTED TO BE AVAILABLE IN GRIB CODE          *
      !*                                                                     *
      !* INPUT:                                                              *
      !* indj               indicates number of the wind field to be read in *
      !* n                  temporal index for meteorological fields (1 to 3)*
      !*                                                                     *
      !* IMPORTANT VARIABLES FROM COMMON BLOCK:                              *
      !*                                                                     *
      !* wfname             File name of data to be read in                  *
      !* nx,ny,nuvz,nwz     expected field dimensions                        *
      !* nlev_ec            number of vertical levels ecmwf model            *
      !* uu,vv,ww           wind fields                                      *
      !* tt,qv              temperature and specific humidity                *
      !* ps                 surface pressure                                 *
      !*                                                                     *
      !***********************************************************************

      use grib_api
      use qvsat_mod

      implicit none

      !HSO  new parameters for grib_api
      integer :: ifile
      integer :: iret, size1, size2, stat
      integer :: igrib
      integer :: ipf
      integer :: gribVer, parCat, parNum, typSfc, valSurf, discipl
      !HSO end edits
      real :: uuh(0:nxmax - 1, 0:nymax - 1, nuvzmax)
      real :: vvh(0:nxmax - 1, 0:nymax - 1, nuvzmax)
      real :: wwh(0:nxmax - 1, 0:nymax - 1, nwzmax)
      integer :: ii, indj, i, j, k, n, levdiff2, ifield, iumax, iwmax

      ! NCEP
      integer :: numpt, numpu, numpv, numpw, numprh, numpclwch
      real :: help, temp
      real :: elev
      real :: ulev1(0:nxmax - 1, 0:nymax - 1), vlev1(0:nxmax - 1, 0:nymax - 1)
      real :: tlev1(0:nxmax - 1, 0:nymax - 1)
      real :: qvh2(0:nxmax - 1, 0:nymax - 1)

      integer :: i179, i180, i181

      ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING
      !HSO kept isec1, isec2 and zsec4 for consistency with gribex GRIB input

      integer :: isec1(8), isec2(3)
      real           :: xsec18  ! IP 29.1.24
      real(kind=4), allocatable, dimension(:) :: zsec4
      real(kind=4) :: xaux, yaux, xaux0, yaux0
      real(kind=8) :: xauxin, yauxin
      real, parameter :: eps = 1.e-4
      real(kind=4) :: ewss(0:nxmax - 1, 0:nymax - 1), nsss(0:nxmax - 1, 0:nymax - 1)
      real :: plev1, hlev1, ff10m, fflev1

      logical :: hflswitch, strswitch

      !HSO  for grib api error messages
      character(len=24) :: gribErrorMsg = 'Error reading grib file'
      character(len=20) :: gribFunction = 'readwind_gfs'
      character(len=20) :: shortname
      character(len=64) :: typeOfLevel, name, stepRange

      if (numpf .gt. 1) goto 777 ! additional precip fields not implemented in GFS

      hflswitch = .false.
      strswitch = .false.
      levdiff2 = nlev_ec - nwz + 1
      iumax = 0
      iwmax = 0

      ! OPENING OF DATA FILE (GRIB CODE)

      !HSO
      call grib_open_file(ifile, path(3) (1:length(3)) &
                          //trim(wfname(indj)), 'r', iret)
      if (iret .ne. GRIB_SUCCESS) then
         goto 888   ! ERROR DETECTED
      end if
      !turn on support for multi fields messages
      call grib_multi_support_on

      numpt = 0
      numpu = 0
      numpv = 0
      numpw = 0
      numprh = 0
      numpclwch = 0
      ifield = 0
      do
         ifield = ifield + 1
         !
         ! GET NEXT FIELDS
         !
         call grib_new_from_file(ifile, igrib, iret)
         if (iret .eq. GRIB_END_OF_FILE) then
            exit   ! EOF DETECTED
         elseif (iret .ne. GRIB_SUCCESS) then
            goto 888   ! ERROR DETECTED
         end if

         !first see if we read GRIB1 or GRIB2
         call grib_get_int(igrib, 'editionNumber', gribVer, iret)
         !  call grib_check(iret,gribFunction,gribErrorMsg)

         if (gribVer .eq. 1) then ! GRIB Edition 1

            !read the grib1 identifiers
            call grib_get_int(igrib, 'indicatorOfParameter', isec1(6), iret)
            !  call grib_check(iret,gribFunction,gribErrorMsg)
            call grib_get_int(igrib, 'indicatorOfTypeOfLevel', isec1(7), iret)
            !  call grib_check(iret,gribFunction,gribErrorMsg)
            call grib_get_int(igrib, 'level', isec1(8), iret)
            !  call grib_check(iret,gribFunction,gribErrorMsg)

! IP 01.24: port form nilu dev branch
!JMA / SH: isec1(8) not evaluated any more below
!b/c with GRIB 2 this may be a real variable
            xsec18 = real(isec1(8))

         else ! GRIB Edition 2

            !read the grib2 identifiers
            call grib_get_string(igrib, 'shortName', shortname, iret)
            call grib_get_string(igrib, 'typeOfLevel', typeOfLevel, iret)
            call grib_get_string(igrib, 'name', name, iret)
            call grib_get_string(igrib, 'stepRange', stepRange, iret)
            if (debug) print *, "Shortname -> ", trim(shortname)//":"//trim(typeOfLevel)//":"//trim(name)//":"//trim(stepRange)

            call grib_get_int(igrib, 'discipline', discipl, iret)
            !  call grib_check(iret,gribFunction,gribErrorMsg)
            call grib_get_int(igrib, 'parameterCategory', parCat, iret)
            !  call grib_check(iret,gribFunction,gribErrorMsg)
            call grib_get_int(igrib, 'parameterNumber', parNum, iret)
            !  call grib_check(iret,gribFunction,gribErrorMsg)
            call grib_get_int(igrib, 'typeOfFirstFixedSurface', typSfc, iret)
            !  call grib_check(iret,gribFunction,gribErrorMsg)
            call grib_get_int(igrib, 'scaledValueOfFirstFixedSurface', &
                              valSurf, iret)
            !  call grib_check(iret,gribFunction,gribErrorMsg)

            !  write(*,*) 'Field: ',ifield,parCat,parNum,typSfc,shortname
            !convert to grib1 identifiers
            isec1(6) = -1
            isec1(7) = -1
            isec1(8) = -1

            xsec18 = -1.0

            if ((parCat .eq. 0) .and. (parNum .eq. 0) .and. (typSfc .eq. 100)) then ! T
               isec1(6) = 11          ! indicatorOfParameter
               isec1(7) = 100         ! indicatorOfTypeOfLevel
               isec1(8) = valSurf/100 ! level, convert to hPa
               xsec18 = valSurf/100.0 ! level, convert to hPa

               ! IPfixgfs11
               call grib_get_size(igrib, 'values', size1, iret)
               allocate (zsec4(size1), stat=stat)
               call grib_get_real4_array(igrib, 'values', zsec4, iret)
               call grib_check(iret, gribFunction, gribErrorMsg)
               deallocate (zsec4)

            elseif ((parCat .eq. 2) .and. (parNum .eq. 2) .and. (typSfc .eq. 100)) then ! U
               isec1(6) = 33          ! indicatorOfParameter
               isec1(7) = 100         ! indicatorOfTypeOfLevel
               isec1(8) = valSurf/100 ! level, convert to hPa
               xsec18 = valSurf/100.0 ! level, convert to hPa

               ! IPfixgfs11
               call grib_get_size(igrib, 'values', size1, iret)
               allocate (zsec4(size1), stat=stat)
               call grib_get_real4_array(igrib, 'values', zsec4, iret)
               call grib_check(iret, gribFunction, gribErrorMsg)
               deallocate (zsec4)

            elseif ((parCat .eq. 2) .and. (parNum .eq. 3) .and. (typSfc .eq. 100)) then ! V
               isec1(6) = 34          ! indicatorOfParameter
               isec1(7) = 100         ! indicatorOfTypeOfLevel
               isec1(8) = valSurf/100 ! level, convert to hPa
               xsec18 = valSurf/100.0 ! level, convert to hPa
            elseif ((parCat .eq. 2) .and. (parNum .eq. 8) .and. (typSfc .eq. 100)) then ! W
               isec1(6) = 39          ! indicatorOfParameter
               isec1(7) = 100         ! indicatorOfTypeOfLevel
               isec1(8) = valSurf/100 ! level, convert to hPa
               xsec18 = valSurf/100.0 ! level, convert to hPa
            elseif ((parCat .eq. 1) .and. (parNum .eq. 1) .and. (typSfc .eq. 100)) then ! RH
               isec1(6) = 52          ! indicatorOfParameter
               isec1(7) = 100         ! indicatorOfTypeOfLevel
               isec1(8) = valSurf/100 ! level, convert to hPa
               xsec18 = valSurf/100.0 ! level, convert to hPa
            elseif ((parCat .eq. 1) .and. (parNum .eq. 1) .and. (typSfc .eq. 103)) then ! RH2
               isec1(6) = 52          ! indicatorOfParameter
               isec1(7) = 105         ! indicatorOfTypeOfLevel
               isec1(8) = 2
               xsec18 = real(2)
            elseif ((parCat .eq. 0) .and. (parNum .eq. 0) .and. (typSfc .eq. 103)) then ! T2
               isec1(6) = 11          ! indicatorOfParameter
               isec1(7) = 105         ! indicatorOfTypeOfLevel
               isec1(8) = 2
               xsec18 = real(2)
            elseif ((parCat .eq. 2) .and. (parNum .eq. 2) .and. (typSfc .eq. 103)) then ! U10
               isec1(6) = 33          ! indicatorOfParameter
               isec1(7) = 105         ! indicatorOfTypeOfLevel
               isec1(8) = 10
               xsec18 = real(10)
            elseif ((parCat .eq. 2) .and. (parNum .eq. 3) .and. (typSfc .eq. 103)) then ! V10
               isec1(6) = 34          ! indicatorOfParameter
               isec1(7) = 105         ! indicatorOfTypeOfLevel
               isec1(8) = 10
               xsec18 = real(10)
            elseif ((parCat .eq. 1) .and. (parNum .eq. 22) .and. (typSfc .eq. 100)) then ! CLWMR Cloud Mixing Ratio [kg/kg]:
               isec1(6) = 153         ! indicatorOfParameter
               isec1(7) = 100         ! indicatorOfTypeOfLevel
               isec1(8) = valSurf/100 ! level, convert to hPa
               xsec18 = valSurf/100.0 ! level, convert to hPa
            elseif ((parCat .eq. 3) .and. (parNum .eq. 1) .and. (typSfc .eq. 101)) then ! SLP
               isec1(6) = 2           ! indicatorOfParameter
               isec1(7) = 102         ! indicatorOfTypeOfLevel
               isec1(8) = 0
               xsec18 = real(0)
            elseif ((parCat .eq. 3) .and. (parNum .eq. 0) .and. (typSfc .eq. 1) .and. (discipl .eq. 0)) then ! SP
               isec1(6) = 1           ! indicatorOfParameter
               isec1(7) = 1           ! indicatorOfTypeOfLevel
               isec1(8) = 0
               xsec18 = real(0)
            elseif ((parCat .eq. 1) .and. (parNum .eq. 13) .and. (typSfc .eq. 1)) then ! SNOW
               isec1(6) = 66          ! indicatorOfParameter
               isec1(7) = 1           ! indicatorOfTypeOfLevel
               isec1(8) = 0
               xsec18 = real(0)
            elseif ((parCat .eq. 0) .and. (parNum .eq. 0) .and. (typSfc .eq. 104)) then ! T sigma 0
               isec1(6) = 11          ! indicatorOfParameter
               isec1(7) = 107         ! indicatorOfTypeOfLevel
               isec1(8) = 0.995       ! lowest sigma level !LB: isec1 is an integer array!!!
               xsec18 = 0.995         ! lowest sigma level
            elseif ((parCat .eq. 2) .and. (parNum .eq. 2) .and. (typSfc .eq. 104)) then ! U sigma 0
               isec1(6) = 33          ! indicatorOfParameter
               isec1(7) = 107         ! indicatorOfTypeOfLevel
               isec1(8) = 0.995       ! lowest sigma level !LB: isec1 is an integer array!!!
               xsec18 = 0.995         ! lowest sigma level
            elseif ((parCat .eq. 2) .and. (parNum .eq. 3) .and. (typSfc .eq. 104)) then ! V sigma 0
               isec1(6) = 34          ! indicatorOfParameter
               isec1(7) = 107         ! indicatorOfTypeOfLevel
               isec1(8) = 0.995       ! lowest sigma level !LB: isec1 is an integer array!!!
               xsec18 = 0.995         ! lowest sigma level
            elseif ((parCat .eq. 3) .and. (parNum .eq. 5) .and. (typSfc .eq. 1)) then ! TOPO
               isec1(6) = 7           ! indicatorOfParameter
               isec1(7) = 1           ! indicatorOfTypeOfLevel
               isec1(8) = 0
               xsec18 = real(0)
            elseif ((parCat .eq. 0) .and. (parNum .eq. 0) .and. (typSfc .eq. 1) &
                    .and. (discipl .eq. 2)) then ! LSM
               isec1(6) = 81          ! indicatorOfParameter
               isec1(7) = 1           ! indicatorOfTypeOfLevel
               isec1(8) = 0
               xsec18 = real(0)
            elseif ((parCat .eq. 3) .and. (parNum .eq. 196) .and. (typSfc .eq. 1)) then ! BLH
               isec1(6) = 221         ! indicatorOfParameter
               isec1(7) = 1           ! indicatorOfTypeOfLevel
               isec1(8) = 0
               xsec18 = real(0)
            elseif ((parCat .eq. 1) .and. (parNum .eq. 7) .and. (typSfc .eq. 1)) then ! LSP/TP
               isec1(6) = 62          ! indicatorOfParameter
               isec1(7) = 1           ! indicatorOfTypeOfLevel
               isec1(8) = 0
               xsec18 = real(0)
            elseif ((parCat .eq. 1) .and. (parNum .eq. 196) .and. (typSfc .eq. 1)) then ! CP
               isec1(6) = 63          ! indicatorOfParameter
               isec1(7) = 1           ! indicatorOfTypeOfLevel
               isec1(8) = 0
               xsec18 = real(0)
            end if

         end if ! gribVer

         if (ifield .eq. 1) then

            !get the required fields from section 2
            !store compatible to gribex input
            call grib_get_int(igrib, 'numberOfPointsAlongAParallel', &
                              isec2(2), iret)
            !  call grib_check(iret,gribFunction,gribErrorMsg)
            call grib_get_int(igrib, 'numberOfPointsAlongAMeridian', &
                              isec2(3), iret)
            !  call grib_check(iret,gribFunction,gribErrorMsg)
            call grib_get_real8(igrib, 'longitudeOfFirstGridPointInDegrees', &
                                xauxin, iret)
            !  call grib_check(iret,gribFunction,gribErrorMsg)
            call grib_get_real8(igrib, 'latitudeOfLastGridPointInDegrees', &
                                yauxin, iret)
            !  call grib_check(iret,gribFunction,gribErrorMsg)
            xaux = real(xauxin) + real(nxshift)*dx
            yaux = real(yauxin)

            ! CHECK GRID SPECIFICATIONS

            if (isec2(2) .ne. nxfield) error stop 'READWIND: NX NOT CONSISTENT'
            if (isec2(3) .ne. ny) error stop 'READWIND: NY NOT CONSISTENT'
            ! if(xaux.eq.0.) xaux=-179.0     ! NCEP DATA
            ! IPfixgfs11: revert to working v10.4 settings

            if (xaux .eq. 0.) xaux = -180.0     ! NCEP DATA

            xaux0 = xlon0
            yaux0 = ylat0
            if (xaux .lt. 0.) xaux = xaux + 360.
            if (yaux .lt. 0.) yaux = yaux + 360.
            if (xaux0 .lt. 0.) xaux0 = xaux0 + 360.
            if (yaux0 .lt. 0.) yaux0 = yaux0 + 360.
            if (abs(xaux - xaux0) .gt. eps) &
               error stop 'READWIND GFS: LOWER LEFT LONGITUDE NOT CONSISTENT'
            if (abs(yaux - yaux0) .gt. eps) &
               error stop 'READWIND GFS: LOWER LEFT LATITUDE NOT CONSISTENT'
         end if
         !HSO end of edits

         if (isec1(6) .ne. -1) then
            !  get the size and data of the values array
            call grib_get_size(igrib, 'values', size1, iret)
            allocate (zsec4(size1), stat=stat)
            if (stat .ne. 0) error stop "Could not allocate zsec4"
            call grib_get_real4_array(igrib, 'values', zsec4, iret)
            call grib_check(iret, gribFunction, gribErrorMsg)
         end if

         i179 = nint(179./dx)
         if (dx .lt. 0.7) then
            i180 = nint(180./dx) + 1    ! 0.5 deg data
         else
            i180 = nint(179./dx) + 1    ! 1 deg data
         end if
         i181 = i180 + 1

         ! IPfixgfs11: revert to v10.4 working settings
         i180 = nint(180./dx)
         i181 = i180
         i179 = i180

         if (isec1(6) .ne. -1) then

            do j = 0, nymin1
               do i = 0, nxfield - 1
                  if ((isec1(6) .eq. 011) .and. (isec1(7) .eq. 100)) then
                     if (debug .and. (i .eq. 0 .and. j .eq. 0)) print *, "Reading -> temperature (tth)"
                     ! TEMPERATURE
                     if ((i .eq. 0) .and. (j .eq. 0)) then
                        !do ii=1,nuvz
                        !  if ((isec1(8)*100.0).eq.akz(ii)) numpt=ii
                        !end do
                        numpt = minloc(abs(xsec18*100.0 - akz), dim=1) ! IP 29.1.24
                        ! IPfixgfs11
                        ! numpt was const 1, and akzs were from not initialized allocation
                     end if
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (help .le. 0) then
                        write (*, *) 'i, j: ', i, j
                        stop 'help <= 0.0 from zsec4'
                     end if
!          if(i.le.i180) then ! 1==180 fills missing 0 lines in tth
                     if (i .lt. i180) then
                        tth(i179 + i, j, numpt, n) = help
                     else
                        tth(i - i181, j, numpt, n) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 033) .and. (isec1(7) .eq. 100)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> u-velocity (uuh)"
                     ! U VELOCITY
                     if ((i .eq. 0) .and. (j .eq. 0)) then
                        ! do ii=1,nuvz
                        !   if ((isec1(8)*100.0).eq.akz(ii)) numpu=ii
                        ! end do
                        numpu = minloc(abs(xsec18*100.0 - akz), dim=1)
                     end if
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        uuh(i179 + i, j, numpu) = help
                     else
                        uuh(i - i181, j, numpu) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 034) .and. (isec1(7) .eq. 100)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> v-velocity (vvh)"
                     ! V VELOCITY
                     if ((i .eq. 0) .and. (j .eq. 0)) then
                        !do ii=1,nuvz
                        !  if ((isec1(8)*100.0).eq.akz(ii)) numpv=ii
                        !end do
                        numpv = minloc(abs(xsec18*100.0 - akz), dim=1)
                     end if
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        vvh(i179 + i, j, numpv) = help
                     else
                        vvh(i - i181, j, numpv) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 052) .and. (isec1(7) .eq. 100)) then
                     ! RELATIVE HUMIDITY -> CONVERT TO SPECIFIC HUMIDITY LATER
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> relative humidity (qvh)"
                     if ((i .eq. 0) .and. (j .eq. 0)) then
!              do ii=1,nuvz
!                if ((isec1(8)*100.0).eq.akz(ii)) numprh=ii
!              end do
                        numprh = minloc(abs(xsec18*100.0 - akz), dim=1)

                     end if
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        qvh(i179 + i, j, numprh, n) = help
                     else
                        qvh(i - i181, j, numprh, n) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 001) .and. (isec1(7) .eq. 001)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> surface pressure (ps)"
                     ! SURFACE PRESSURE
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        ps(i179 + i, j, 1, n) = help
                     else
                        ps(i - i181, j, 1, n) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 039) .and. (isec1(7) .eq. 100)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> vertical velocity (wwh)"
                     ! W VELOCITY
!           if((i.eq.0).and.(j.eq.0)) then
!              do ii=1,nuvz
!                if ((isec1(8)*100.0).eq.akz(ii)) numpw=ii
!              end do
!          endif
                     if ((i .eq. 0) .and. (j .eq. 0)) numpw = minloc(abs(xsec18*100.0 - akz), dim=1)

                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        wwh(i179 + i, j, numpw) = help
                     else
                        wwh(i - i181, j, numpw) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 066) .and. (isec1(7) .eq. 001)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> snow depth (sd)"
                     ! SNOW DEPTH
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        sd(i179 + i, j, 1, n) = help
                     else
                        sd(i - i181, j, 1, n) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 002) .and. (isec1(7) .eq. 102)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> mean sea level pressure (mslp)"
                     ! MEAN SEA LEVEL PRESSURE
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        msl(i179 + i, j, 1, n) = help
                     else
                        msl(i - i181, j, 1, n) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 071) .and. (isec1(7) .eq. 244)) then
                     ! TOTAL CLOUD COVER
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        tcc(i179 + i, j, 1, n) = help
                     else
                        tcc(i - i181, j, 1, n) = help
                     end if
                  end if
                  ! if((isec1(6).eq.033).and.(isec1(7).eq.105).and. &
                  !      (nint(xsec18).eq.10)) then
!             (isec1(8).eq.10)) then
                  if (trim(shortName) .eq. '10u') then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> 10m u velocity (u10)"
                     ! 10 M U VELOCITY
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        u10(i179 + i, j, 1, n) = help
                     else
                        u10(i - i181, j, 1, n) = help
                     end if
                  end if
                  ! if ((isec1(6) .eq. 034) .and. (isec1(7) .eq. 105) .and. &
                  !     (nint(xsec18) .eq. 10)) then
                  if (trim(shortName) .eq. '10v') then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> 10m v velocity (v10)"
                     ! 10 M V VELOCITY
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        v10(i179 + i, j, 1, n) = help
                     else
                        v10(i - i181, j, 1, n) = help
                     end if
                  end if
                  ! if ((isec1(6) .eq. 011) .and. (isec1(7) .eq. 105) .and. &
                  !     (nint(xsec18) .eq. 2)) then
!             (isec1(8).eq.02)) then
                  if (trim(shortName) .eq. '2t') then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> 2m temperature (tt2)"
                     ! 2 M TEMPERATURE
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        tt2(i179 + i, j, 1, n) = help
                     else
                        tt2(i - i181, j, 1, n) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 017) .and. (isec1(7) .eq. 105) .and. &
                      (nint(xsec18) .eq. 2)) then
                     !      (isec1(8).eq.02)) then

                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> 2m dew point temperature (td2)"
                     ! 2 M DEW POINT TEMPERATURE
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        td2(i179 + i, j, 1, n) = help
                     else
                        td2(i - i181, j, 1, n) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 062) .and. (isec1(7) .eq. 001)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> large scale precipitation (lsprec)"
                     ! LARGE SCALE PREC.
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        lsprec(i179 + i, j, 1, 1, n) = help
                     else
                        lsprec(i - i181, j, 1, 1, n) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 063) .and. (isec1(7) .eq. 001)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> convective precipitation (convprec)"
                     ! CONVECTIVE PREC.
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        convprec(i179 + i, j, 1, 1, n) = help
                     else
                        convprec(i - i181, j, 1, 1, n) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 007) .and. (isec1(7) .eq. 001)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> topography (oro)"
                     ! TOPOGRAPHY
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        oro(i179 + i, j) = help
                        excessoro(i179 + i, j) = 0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
                     else
                        oro(i - i181, j) = help
                        excessoro(i - i181, j) = 0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
                     end if
                  end if
                  if ((isec1(6) .eq. 081) .and. (isec1(7) .eq. 001)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> land/sea mask (lsm)"
                     ! LAND SEA MASK
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        lsm(i179 + i, j) = help
                     else
                        lsm(i - i181, j) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 221) .and. (isec1(7) .eq. 001)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> mixing height (hmix)"
                     ! MIXING HEIGHT
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        hmix(i179 + i, j, 1, n) = help
                     else
                        hmix(i - i181, j, 1, n) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 052) .and. (isec1(7) .eq. 105) .and. &
                      (isec1(8) .eq. 02)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> 2 m relative humidity (qvh2)"
                     ! 2 M RELATIVE HUMIDITY
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        qvh2(i179 + i, j) = help
                     else
                        qvh2(i - i181, j) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 011) .and. (isec1(7) .eq. 107)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> temperature at lowest sigma level (tlev1)"
                     ! TEMPERATURE LOWEST SIGMA LEVEL
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        tlev1(i179 + i, j) = help
                     else
                        tlev1(i - i181, j) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 033) .and. (isec1(7) .eq. 107)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> u-velocity at lowest sigma level (ulev1)"
                     ! U VELOCITY LOWEST SIGMA LEVEL
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        ulev1(i179 + i, j) = help
                     else
                        ulev1(i - i181, j) = help
                     end if
                  end if
                  if ((isec1(6) .eq. 034) .and. (isec1(7) .eq. 107)) then
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> v-velocity at lowest sigma level (vlev1)"
                     ! V VELOCITY LOWEST SIGMA LEVEL
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        vlev1(i179 + i, j) = help
                     else
                        vlev1(i - i181, j) = help
                     end if
                  end if
                  ! SEC & IP 12/2018 read GFS clouds
                  if ((isec1(6) .eq. 153) .and. (isec1(7) .eq. 100)) then  !! CLWCR  Cloud liquid water content [kg/kg]
                     if (debug .and. i .eq. 0 .and. j .eq. 0) print *, "Reading -> cloud liquid water content (clwch)"
                     if ((i .eq. 0) .and. (j .eq. 0)) then
                        !  do ii=1,nuvz
                        !1    if ((isec1(8)*100.0).eq.akz(ii)) numpclwch=ii
                        ! end do
                        numpclwch = minloc(abs(xsec18*100.0 - akz), dim=1)
                     end if
                     help = zsec4(nxfield*(ny - j - 1) + i + 1)
                     if (i .lt. i180) then
                        clwch(i179 + i, j, numpclwch, n) = help
                     else
                        clwch(i - i181, j, numpclwch, n) = help
                     end if
                     lcw = .true.
                     lcwsum = .true.
                  end if

               end do
            end do

         end if

         if ((isec1(6) .eq. 33) .and. (isec1(7) .eq. 100)) then
            ! NCEP ISOBARIC LEVELS
            iumax = iumax + 1
         end if

         call grib_release(igrib)

         if (isec1(6) .ne. -1) deallocate (zsec4) !IP 28/11/23 fix deallocation error
      end do                      !! READ NEXT LEVEL OR PARAMETER
      !
      ! CLOSING OF INPUT DATA FILE
      !

      !HSO close grib file
      call grib_close_file(ifile)

      ! SENS. HEAT FLUX
      sshf(:, :, 1, n) = 0.0     ! not available from gfs.tccz.pgrbfxx files
      hflswitch = .false.    ! Heat flux not available
      ! SOLAR RADIATIVE FLUXES
      ssr(:, :, 1, n) = 0.0      ! not available from gfs.tccz.pgrbfxx files
      ! EW SURFACE STRESS
      ewss = 0.0         ! not available from gfs.tccz.pgrbfxx files
      ! NS SURFACE STRESS
      nsss = 0.0         ! not available from gfs.tccz.pgrbfxx files
      strswitch = .false.    ! stress not available

      ! CONVERT TP TO LSP (GRIB2 only)
      if (gribVer .eq. 2) then
         do j = 0, nymin1
         do i = 0, nxfield - 1
            if (i .le. i180) then
            if (convprec(i179 + i, j, 1, 1, n) .lt. lsprec(i179 + i, j, 1, 1, n)) then ! neg precip would occur
               lsprec(i179 + i, j, 1, 1, n) = &
                  lsprec(i179 + i, j, 1, 1, n) - convprec(i179 + i, j, 1, 1, n)
            else
               lsprec(i179 + i, j, 1, 1, n) = 0
            end if
            else
            if (convprec(i - i181, j, 1, 1, n) .lt. lsprec(i - i181, j, 1, 1, n)) then
               lsprec(i - i181, j, 1, 1, n) = &
                  lsprec(i - i181, j, 1, 1, n) - convprec(i - i181, j, 1, 1, n)
            else
               lsprec(i - i181, j, 1, 1, n) = 0
            end if
            end if
         end do
         end do
      end if
      !HSO end edits

      ! TRANSFORM RH TO SPECIFIC HUMIDITY

      do j = 0, ny - 1
         do i = 0, nxfield - 1
            do k = 1, nuvz
               help = qvh(i, j, k, n)
               temp = tth(i, j, k, n)
               if (temp .le. 0.0) then
                  write (*, *) 'STOP: TRANSFORM RH TO SPECIFIC HUMIDITY: temp, i, j, k, n'
                  write (*, *) temp, i, j, k, n
!          temp = 273.0
                  stop
               end if

               plev1 = akm(k) + bkm(k)*ps(i, j, 1, n)
               !print*, temp,plev1
               elev = ew(temp, plev1)*help/100.0
               qvh(i, j, k, n) = xmwml*(elev/(plev1 - ((1.0 - xmwml)*elev)))
            end do
         end do
      end do

      ! CALCULATE 2 M DEW POINT FROM 2 M RELATIVE HUMIDITY
      ! USING BOLTON'S (1980) FORMULA
      ! BECAUSE td2 IS NOT AVAILABLE FROM NCEP GFS DATA
      k = 2 ! CHECK THIS!!!
      do j = 0, ny - 1
         do i = 0, nxfield - 1
            help = qvh2(i, j)
            temp = tt2(i, j, 1, n)
            if (temp .le. 0.0) then
               write (*, *) 'STOP: CALCULATE 2 M DEW POINT FROM 2 M RELATIVE HUMIDITY: temp, i, j, k, n'
               write (*, *) temp, i, j, k, n
!          temp = 273.0
               stop
            end if

            plev1 = akm(k) + bkm(k)*ps(i, j, 1, n)
            elev = ew(temp, plev1)/100.*help/100.   !vapour pressure in hPa
            td2(i, j, 1, n) = 243.5/(17.67/log(elev/6.112) - 1) + 273.
            if (help .le. 0.) td2(i, j, 1, n) = tt2(i, j, 1, n)
         end do
      end do

      if (levdiff2 .eq. 0) then
         iwmax = nlev_ec + 1
         do i = 0, nxmin1
            do j = 0, nymin1
               wwh(i, j, nlev_ec + 1) = 0.
            end do
         end do
      end if

      ! For global fields, assign the leftmost data column also to the rightmost
      ! data column; if required, shift whole grid by nxshift grid points
      !*************************************************************************

      if (xglobal) then
         call shift_field_0(ewss, nxfield, ny)
         call shift_field_0(nsss, nxfield, ny)
         call shift_field_0(oro, nxfield, ny)
         call shift_field_0(excessoro, nxfield, ny)
         call shift_field_0(lsm, nxfield, ny)
         call shift_field_0(ulev1, nxfield, ny)
         call shift_field_0(vlev1, nxfield, ny)
         call shift_field_0(tlev1, nxfield, ny)
         call shift_field_0(qvh2, nxfield, ny)
         call shift_field(ps, nxfield, ny, 1, 1, 2, n)
         call shift_field(sd, nxfield, ny, 1, 1, 2, n)
         call shift_field(msl, nxfield, ny, 1, 1, 2, n)
         call shift_field(tcc, nxfield, ny, 1, 1, 2, n)
         call shift_field(u10, nxfield, ny, 1, 1, 2, n)
         call shift_field(v10, nxfield, ny, 1, 1, 2, n)
         call shift_field(tt2, nxfield, ny, 1, 1, 2, n)
         call shift_field(td2, nxfield, ny, 1, 1, 2, n)
         do ipf = 1, numpf
            call shift_field(lsprec(:, :, :, ipf, n), nxfield, ny, 1, 1, 1, 1)
            call shift_field(convprec(:, :, :, ipf, n), nxfield, ny, 1, 1, 1, 1)
         end do
         call shift_field(sshf, nxfield, ny, 1, 1, 2, n)
         call shift_field(ssr, nxfield, ny, 1, 1, 2, n)
         call shift_field(hmix, nxfield, ny, 1, 1, 2, n)
         call shift_field(tth, nxfield, ny, nuvzmax, nuvz, 2, n)
         call shift_field(qvh, nxfield, ny, nuvzmax, nuvz, 2, n)
         call shift_field(uuh, nxfield, ny, nuvzmax, nuvz, 1, 1)
         call shift_field(vvh, nxfield, ny, nuvzmax, nuvz, 1, 1)
         call shift_field(wwh, nxfield, ny, nwzmax, nwz, 1, 1)
         ! IP & SEC adding GFS Clouds 20181205
         call shift_field(clwch, nxfield, ny, nuvzmax, nuvz, 2, n)
      end if

      do i = 0, nxmin1
         do j = 0, nymin1
            ! Convert precip. from mm/s -> mm/hour
            convprec(i, j, 1, 1, n) = convprec(i, j, 1, 1, n)*3600.
            lsprec(i, j, 1, 1, n) = lsprec(i, j, 1, 1, n)*3600.
            sfcstress(i, j, 1, n) = sqrt(ewss(i, j)**2 + nsss(i, j)**2)
         end do
      end do

      if ((.not. hflswitch) .or. (.not. strswitch)) then
         !  write(*,*) 'WARNING: No flux data contained in GRIB file ',
         !    +  wfname(indj)

         ! CALCULATE USTAR AND SSHF USING THE PROFILE METHOD
         !***************************************************************************

         do j = 0, nymin1
            do i = 0, nxmin1
               hlev1 = 30.0                     ! HEIGHT OF FIRST MODEL SIGMA LAYER
               ff10m = sqrt(u10(i, j, 1, n)**2 + v10(i, j, 1, n)**2)
               fflev1 = sqrt(ulev1(i, j)**2 + vlev1(i, j)**2)
               call pbl_profile(ps(i, j, 1, n), td2(i, j, 1, n), hlev1, &
                                tt2(i, j, 1, n), tlev1(i, j), ff10m, fflev1, &
                                sfcstress(i, j, 1, n), sshf(i, j, 1, n))
               if (sshf(i, j, 1, n) .gt. 200.) sshf(i, j, 1, n) = 200.
               if (sshf(i, j, 1, n) .lt. -400.) sshf(i, j, 1, n) = -400.
            end do
         end do
      end if

      if (iumax .ne. nuvz) error stop 'READWIND: NUVZ NOT CONSISTENT'
      if (iumax .ne. nwz) error stop 'READWIND: NWZ NOT CONSISTENT'

      return

777   continue
      write (*, *) ' #### FLEXPART MODEL ERROR!                       ####'
      write (*, *) ' #### Additional precip fields not implemented in ####'
      write (*, *) ' #### GFS version                                 ####'
      write (*, *) ' #### Set numpf=1 in par_mod.f90 and recompile!   ####'
      stop 'Execution terminated'

888   write (*, *) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
      write (*, *) ' #### ', wfname(indj), '                    #### '
      write (*, *) ' #### IS NOT GRIB FORMAT !!!                  #### '
      error stop 'Execution terminated'

   end subroutine readwind_gfs

   subroutine shift_field_0(field, nxf, nyf)
      !                          i/o   i   i
      !*****************************************************************************
      !                                                                            *
      !  This subroutine shifts global fields by nxshift grid cells, in order to   *
      !  facilitate all sorts of nested wind fields, or output grids, which,       *
      !  without shifting, would overlap with the domain "boundary".               *
      !                                                                            *
      !    Author: A. Stohl                                                        *
      !                                                                            *
      !    3 July 2002                                                             *
      !                                                                            *
      !*****************************************************************************
      !                                                                            *
      ! Variables:                                                                 *
      !                                                                            *
      ! Constants:                                                                 *
      !                                                                            *
      !*****************************************************************************

      implicit none

      integer :: nxf, nyf, ix, jy, ixs
      real :: field(0:nxmax - 1, 0:nymax - 1), xshiftaux(0:nxmax - 1)

      ! Loop over y and z
      !******************

      do jy = 0, nyf - 1

         ! Shift the data
         !***************

         if (nxshift .ne. 0) then
            do ix = 0, nxf - 1
               if (ix .ge. nxshift) then
                  ixs = ix - nxshift
               else
                  ixs = nxf - nxshift + ix
               end if
               xshiftaux(ixs) = field(ix, jy)
            end do
            do ix = 0, nxf - 1
               field(ix, jy) = xshiftaux(ix)
            end do
         end if

         ! Repeat the westernmost grid cells at the easternmost domain "boundary"
         !***********************************************************************

         field(nxf, jy) = field(0, jy)
      end do

      return
   end subroutine shift_field_0

   subroutine shift_field(field, nxf, nyf, nzfmax, nzf, nmax, n)
      !                        i/o   i   i    i     i   i   i
      !*****************************************************************************
      !                                                                            *
      !  This subroutine shifts global fields by nxshift grid cells, in order to   *
      !  facilitate all sorts of nested wind fields, or output grids, which,       *
      !  without shifting, would overlap with the domain "boundary".               *
      !                                                                            *
      !    Author: A. Stohl                                                        *
      !                                                                            *
      !    3 July 2002                                                             *
      !                                                                            *
      !*****************************************************************************
      !                                                                            *
      ! Variables:                                                                 *
      !                                                                            *
      ! Constants:                                                                 *
      !                                                                            *
      !*****************************************************************************

      implicit none

      integer :: nxf, nyf, nzf, n, ix, jy, kz, ixs, nzfmax, nmax
      real :: field(0:nxmax - 1, 0:nymax - 1, nzfmax, nmax), xshiftaux(0:nxmax - 1)

      ! Loop over y and z
      !******************

      do kz = 1, nzf
         do jy = 0, nyf - 1

            ! Shift the data
            !***************

            if (nxshift .ne. 0) then
               do ix = 0, nxf - 1
                  if (ix .ge. nxshift) then
                     ixs = ix - nxshift
                  else
                     ixs = nxf - nxshift + ix
                  end if
                  xshiftaux(ixs) = field(ix, jy, kz, n)
               end do
               do ix = 0, nxf - 1
                  field(ix, jy, kz, n) = xshiftaux(ix)
               end do
            end if

            ! Repeat the westernmost grid cells at the easternmost domain "boundary"
            !***********************************************************************

            field(nxf, jy, kz, n) = field(0, jy, kz, n)
         end do
      end do
   end subroutine shift_field

   subroutine alloc_fixedfields
      implicit none
      integer :: stat

      allocate (oro(0:nxmax - 1, 0:nymax - 1), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate oro"
      allocate (excessoro(0:nxmax - 1, 0:nymax - 1), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate excessoro"
      allocate (lsm(0:nxmax - 1, 0:nymax - 1), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate lsm"
   end subroutine alloc_fixedfields

   subroutine alloc_windfields
      implicit none
      integer :: stat
      ! Eta coordinates
      !****************
      ! Intrinsic coordinates
      !**********************
      allocate (uu(0:nxmax - 1, 0:nymax - 1, nzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate uu"
      allocate (vv(0:nxmax - 1, 0:nymax - 1, nzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate vv"
      allocate (ww(0:nxmax - 1, 0:nymax - 1, nzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate ww"
      allocate (uupol(0:nxmax - 1, 0:nymax - 1, nzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate uupol"
      allocate (vvpol(0:nxmax - 1, 0:nymax - 1, nzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate vvpol"
      allocate (tt(0:nxmax - 1, 0:nymax - 1, nzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate tt"
      allocate (tth(0:nxmax - 1, 0:nymax - 1, nuvzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate tth"
      allocate (pv(0:nxmax - 1, 0:nymax - 1, nzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate pv"
      allocate (qv(0:nxmax - 1, 0:nymax - 1, nzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate qv"
      allocate (qvh(0:nxmax - 1, 0:nymax - 1, nuvzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate qvh"
      allocate (rho(0:nxmax - 1, 0:nymax - 1, nzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate rho"
      allocate (drhodz(0:nxmax - 1, 0:nymax - 1, nzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate drhodz"
      allocate (pplev(0:nxmax - 1, 0:nymax - 1, nuvzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate pplev"
      allocate (prs(0:nxmax - 1, 0:nymax - 1, nzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate prs"
      allocate (rho_dry(0:nxmax - 1, 0:nymax - 1, nzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate rho_dry"

      ! Cloud data
      !***********
      allocate (clwc(0:nxmax - 1, 0:nymax - 1, nzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate clwc"
      allocate (ciwc(0:nxmax - 1, 0:nymax - 1, nzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate ciwc"
      allocate (clwch(0:nxmax - 1, 0:nymax - 1, nuvzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate clwch"
      allocate (ciwch(0:nxmax - 1, 0:nymax - 1, nuvzmax, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate ciwch"

      clwc = 0.0
      ciwc = 0.0
      clwch = 0.0
      ciwch = 0.0

      allocate (ctwc(0:nxmax - 1, 0:nymax - 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate ctwc"
      allocate (icloudbot(0:nxmax - 1, 0:nymax - 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate icloudbot"
      allocate (icloudtop(0:nxmax - 1, 0:nymax - 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate icloudtop"

      ! 2d fields
      !**********
      allocate (ps(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate ps"
      allocate (sd(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate sd"
      allocate (msl(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate msl"
      allocate (tcc(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate tcc"
      allocate (u10(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate u10"
      allocate (v10(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate v10"
      allocate (tt2(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate tt2"
      allocate (td2(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate td2"
      allocate (lsprec(0:nxmax - 1, 0:nymax - 1, 1, numpf, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate lsprec"
      allocate (convprec(0:nxmax - 1, 0:nymax - 1, 1, numpf, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate convprec"
      allocate (sshf(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate sshf"
      allocate (ssr(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate ssr"
      allocate (sfcstress(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate sfcstress"
      allocate (ustar(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate ustar"
      allocate (wstar(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate wstar"
      allocate (hmix(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate hmix"
      allocate (tropopause(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate tropopause"
      allocate (oli(0:nxmax - 1, 0:nymax - 1, 1, numwfmem), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate oli"

      ! Vertical descritisation arrays
      !*******************************
      allocate (height(nzmax), wheight(nzmax), uvheight(nzmax), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate height arrays"
      allocate (akm(nwzmax), bkm(nwzmax), akz(nuvzmax), bkz(nuvzmax), &
                aknew(nzmax), bknew(nzmax), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate model level parameters"
   end subroutine alloc_windfields

   subroutine dealloc_windfields
      implicit none

      deallocate (wftime, wfname)
      deallocate (oro, excessoro, lsm)
      deallocate (uu, vv, ww, uupol, vvpol, tt, tth, qv, qvh, pv, rho, drhodz, pplev, prs, rho_dry)
      deallocate (clwc, ciwc, clwch, ciwch, ctwc)
      deallocate (ps, sd, msl, tcc, u10, v10, tt2, td2, lsprec, convprec, sshf, ssr, sfcstress, &
                  ustar, wstar, hmix, tropopause, oli)
      deallocate (height, wheight, uvheight, akm, bkm, akz, bkz, aknew, bknew)
   end subroutine dealloc_windfields

end module windfields_mod
