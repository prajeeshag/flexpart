module interpol_mod

  !use par_mod, only: nzmax, maxspec

  use par_mod
  use com_mod

  implicit none

  real :: uprof(nzmax),vprof(nzmax),wprof(nzmax),wprofeta(nzmax),detaprof(nzmax)
  real :: usigprof(nzmax),vsigprof(nzmax),wsigprof(nzmax),wsigprofeta(nzmax)
  real :: rhoprof(nzmax),rhogradprof(nzmax)

  real :: u,v,w,usig,vsig,wsig,pvi,ueta,veta,weta,wsigeta

  real :: p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2
  real :: xtn,ytn 
  integer :: nix,njy
  integer :: ix,jy,ixp,jyp,ngrid,indz,indzp,indzeta,indzpeta
  integer :: induv,indpuv
  logical :: depoindicator(maxspec)
  logical :: indzindicator(nzmax)
  logical :: lbounds(2),lbounds_w(2),lbounds_uv(2) ! marking particles below or above bounds

!$OMP THREADPRIVATE(uprof,vprof,wprof,usigprof,vsigprof,wsigprof, &
!$OMP rhoprof,rhogradprof,u,v,w,usig,vsig,wsig,pvi, &
!$OMP p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2,ix,jy,ixp,jyp, &
!$OMP ngrid,indz,indzp,depoindicator,indzindicator, &
!$OMP wprofeta,wsigprofeta,induv,indpuv,lbounds,lbounds_w,lbounds_uv, &
!$OMP indzeta,indzpeta,ueta,veta,weta,wsigeta,detaprof, &
!$OMP xtn,ytn,nix,njy)

contains

subroutine initialise_interpol_mod(itime,xt,yt,zt,zteta)
  ! This routine initialises all important values used in the interpol module
  ! This includes:
  ! - The current grid number in which the particle is positioned
  ! - The interpolation fractions of the grid (x,y,z) and of time

  implicit none

  integer, intent(in) :: itime             ! time step
  real, intent(in)    :: xt,yt             ! particle positions
  real, intent(in)    :: zt                ! height in meters
  real, intent(in)    :: zteta             ! height in eta coordinates

  call determine_grid_coordinates(xt,yt)
  call find_grid_distances(xt,yt)
  call find_time_variables(itime)
  call find_z_level_meters(zt)
  call find_z_level_eta(zteta)
end subroutine initialise_interpol_mod

subroutine determine_grid_coordinates(xt,yt)
  implicit none 

  real, intent(in) :: xt,yt                 ! particle positions

  if (ngrid.gt.0) then
    xtn=(xt-xln(ngrid))*xresoln(ngrid)
    ytn=(yt-yln(ngrid))*yresoln(ngrid)
    ix=int(xtn)
    jy=int(ytn)
    nix=nint(xtn)
    njy=nint(ytn)
  else
    ix=int(xt)
    jy=int(yt)
    nix=nint(xt)
    njy=nint(yt)
  endif
  ixp=ix+1
  jyp=jy+1

  ! eso: Temporary fix for particle exactly at north pole
  if (jyp >= nymax) then
    ! write(*,*) 'WARNING: advance.f90 jyp >= nymax. xt,yt:',xt,yt
    jyp=jyp-1
  end if

  if (jyp >= nymax) then
    jyp=jyp-1
  end if
end subroutine determine_grid_coordinates

subroutine find_grid_distances(xt,yt)

  implicit none 

  real, intent(in) :: xt,yt                 ! particle positions

  ddx=xt-real(ix)
  ddy=yt-real(jy)
  rddx=1.-ddx
  rddy=1.-ddy
  p1=rddx*rddy
  p2=ddx*rddy
  p3=rddx*ddy
  p4=ddx*ddy
end subroutine find_grid_distances

subroutine find_time_variables(itime)
  
  implicit none  

  integer, intent(in) :: itime             ! time step
  
  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)  
end subroutine find_time_variables

subroutine find_z_level_meters(zt)
  implicit none
  real, intent(in)     :: zt       ! height in meters
  integer              :: i

  indz=nz-1
  indzp=nz
  if (zt.le.height(1)) then
    lbounds(1)=.true.
    lbounds(2)=.false.
    indz=1
    indzp=2
  else if (zt.ge.height(nz)) then
    lbounds(1)=.false.
    lbounds(2)=.true.
  else
    lbounds(1)=.false.
    lbounds(2)=.false.
    do i=2,nz
      if (height(i).gt.zt) then
        indz=i-1
        indzp=i
        exit
      endif
    end do
  endif
end subroutine find_z_level_meters

subroutine find_z_level_eta(zteta)
  implicit none
  real, intent(in)       :: zteta    ! height in eta coordinates
  integer                :: i        ! loop variable

  indzeta=nz-1
  indzpeta=nz
  ! Flag particles that are above or below bounds
  if (zteta.ge.wheight(1)) then
    lbounds_w(1)=.true.
    lbounds_w(2)=.false.
    indzeta=1
    indzpeta=2
  else if (zteta.le.wheight(nz)) then
    lbounds_w(1)=.false.
    lbounds_w(2)=.true.
  else
    lbounds_w(1)=.false.
    lbounds_w(2)=.false.
    do i=2,nz
      if (wheight(i).lt.zteta) then
        indzeta=i-1
        indzpeta=i
        exit
      endif
    end do
  endif

  induv=nz-1
  indpuv=nz
  if (zteta.gt.uvheight(1)) then
    lbounds_uv(1)=.true.
    lbounds_uv(2)=.false.
    induv=1 
    indpuv=2
  else if (zteta.lt.uvheight(nz)) then
    lbounds_uv(1)=.false.
    lbounds_uv(2)=.true.
  else
    lbounds_uv(1)=.false.
    lbounds_uv(2)=.false.
    do i=2,nz
      if (uvheight(i).lt.zteta) then
        induv=i-1
        indpuv=i
        exit
      endif
    end do
  endif
end subroutine find_z_level_eta

subroutine find_vertical_variables(vertlevels,zpos,zlevel,dz1,dz2,bounds)
  implicit none
  real, intent(in)    :: vertlevels(:)     ! vertical levels in coordinate system
  real, intent(in)    :: zpos              ! verticle particle position
  integer, intent(in) :: zlevel            ! vertical level of interest
  logical, intent(in) :: bounds(2)         ! flag marking if particles are outside bounds  
  real, intent(inout) :: dz1,dz2
  real                :: dz

  ! If the particle is below bounds (bounds(1)==.true.):
  if (bounds(1)) then
    dz1=0.
    dz2=1.
  ! If above bounds (bounds(2)==.true.):
  else if (bounds(2)) then
    dz1=1.
    dz2=0.
  else
    dz=1./(vertlevels(zlevel+1)-vertlevels(zlevel))
    dz1=(zpos-vertlevels(zlevel))*dz
    dz2=(vertlevels(zlevel+1)-zpos)*dz
  endif
end subroutine find_vertical_variables

subroutine temporal_interpolation(time1,time2,output)

  implicit none

  real, intent(in)    :: time1,time2     ! input data at two timesteps 
  real, intent(inout) :: output          ! interpolated data

  output=(time1*dt2+time2*dt1)*dtt
end subroutine temporal_interpolation

subroutine bilinear_horizontal_interpolation(field,output,zlevel,ztot)

  implicit none
  integer, intent(in) :: zlevel,ztot                              ! interpolation z level, z
  real, intent(in)    :: field(0:nxmax-1,0:nymax-1,ztot,numwfmem)    ! input field to interpolate over
  real, intent(inout) :: output(2)                                ! interpolated values
  integer             :: m, indexh

  do m=1,2
    indexh=memind(m)

    output(m)=p1*field(ix ,jy ,zlevel,indexh) &
         + p2*field(ixp,jy ,zlevel,indexh) &
         + p3*field(ix ,jyp,zlevel,indexh) &
         + p4*field(ixp,jyp,zlevel,indexh)
  end do
end subroutine bilinear_horizontal_interpolation

subroutine bilinear_spatial_interpolation(field,output,zlevel,dz1,dz2,ztot)
  implicit none
  integer, intent(in) :: zlevel,ztot                               ! interpolation z level
  real, intent(in)    :: field(0:nxmax-1,0:nymax-1,ztot,numwfmem)  ! input field to interpolate over
  real, intent(in)    :: dz1,dz2
  real, intent(inout) :: output(2)                                 ! interpolated values
  integer             :: m,n,indexh,indzh
  real                :: output1(2)

  do m=1,2
    indexh=memind(m)
    do n=1,2
      indzh=zlevel+n-1
      output1(n)=p1*field(ix ,jy ,indzh,indexh) &
           + p2*field(ixp,jy ,indzh,indexh) &
           + p3*field(ix ,jyp,indzh,indexh) &
           + p4*field(ixp,jyp,indzh,indexh)
    end do
  !**********************************
  ! 2.) Linear vertical interpolation
  !**********************************
    output(m)=dz2*output1(1)+dz1*output1(2)    
  end do
end subroutine bilinear_spatial_interpolation
  
subroutine compute_standard_deviation(field,output,zlevel1,zlevel2,ztot)
  implicit none
  real,parameter      :: eps=1.0e-30
  integer, intent(in) :: zlevel1,zlevel2,ztot                     ! interpolation z levels
  real, intent(in)    :: field(0:nxmax-1,0:nymax-1,ztot,numwfmem) ! input field to interpolate over
  real, intent(inout) :: output                                   ! standard deviation
  real                :: xaux,ndivide
  real                :: sl, sq
  integer             :: m, indexh,zlevel

  sl=0.
  sq=0.
  do m=1,2
    indexh=memind(m)
    do zlevel=zlevel1,zlevel2
      sl=sl+field(ix ,jy ,zlevel,indexh)+field(ixp,jy ,zlevel,indexh) &
           +field(ix ,jyp,zlevel,indexh)+field(ixp,jyp,zlevel,indexh)
      sq=sq+field(ix ,jy ,zlevel,indexh)*field(ix ,jy ,zlevel,indexh)+ &
           field(ixp,jy ,zlevel,indexh)*field(ixp,jy ,zlevel,indexh)+ &
           field(ix ,jyp,zlevel,indexh)*field(ix ,jyp,zlevel,indexh)+ &
           field(ixp,jyp,zlevel,indexh)*field(ixp,jyp,zlevel,indexh)
    end do
  end do

  if (zlevel1.eq.zlevel2) then
    ndivide=8.
  else
    ndivide=16.
  endif

  xaux=sq-sl*sl/ndivide
  if (xaux.lt.eps) then
    output=0.
  else
    output=sqrt(xaux/(ndivide-1.))
  endif  
end subroutine compute_standard_deviation

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

  use hanna_mod

  implicit none

  integer, intent(in) :: itime
  real, intent(in)    :: xt,yt,zt,zteta

  ! Auxiliary variables needed for interpolation
  real :: ust1(2),wst1(2),oli1(2),oliaux
  real :: y1(2),y2(2),y3(2),rho1(2),rhograd1(2)
  integer :: i,m,n,indexh
  real,parameter :: eps=1.0e-30

  !********************************************
  ! Multilinear interpolation in time and space
  !********************************************

  ! Determine the lower left corner and its distance to the current position
  !*************************************************************************
  call find_grid_distances(xt,yt)

  ! Calculate variables for time interpolation
  !*******************************************
  call find_time_variables(itime)

  !*****************************************
  ! 1. Interpolate u*, w* and Obukhov length
  !*****************************************

  ! a) Bilinear horizontal interpolation
  call bilinear_horizontal_interpolation(ustar,ust1,1,1)
  call bilinear_horizontal_interpolation(wstar,wst1,1,1)
  call bilinear_horizontal_interpolation(oli,oli1,1,1)

  ! b) Temporal interpolation
  call temporal_interpolation(ust1(1),ust1(2),ust)
  call temporal_interpolation(wst1(1),wst1(2),wst)
  call temporal_interpolation(oli1(1),oli1(2),oliaux)

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
  call find_z_level_meters(zt)

  !**************************************
  ! 1.) Bilinear horizontal interpolation
  ! 2.) Temporal interpolation (linear)
  !**************************************

  ! Loop over 2 time steps and indz levels
  !***************************************

  do n=indz,indzp
    call bilinear_horizontal_interpolation(ww,y3,n,nwzmax)
    call temporal_interpolation(y3(1),y3(2),wprof(n))
    indzindicator(n)=.false.

  ! Compute standard deviations
  !****************************
    call compute_standard_deviation(ww,wsigprof(n),n,n,nwzmax)
  end do

  ! Same for zt in eta coordinates
  !*******************************
  call find_z_level_eta(zteta)

  !**************************************
  ! 1.) Bilinear horizontal interpolation
  ! 2.) Temporal interpolation (linear)
  !**************************************

  ! Loop over 2 time steps and indz levels
  !***************************************

  do n=induv,indpuv
    if (ngrid.lt.0) then
      call bilinear_horizontal_interpolation(uupoleta,y1,n,nzmax)
      call bilinear_horizontal_interpolation(vvpoleta,y2,n,nzmax)
      call compute_standard_deviation(uupoleta,usigprof(n),n,n,nzmax)
      call compute_standard_deviation(vvpoleta,vsigprof(n),n,n,nzmax)
    else
      call bilinear_horizontal_interpolation(uueta,y1,n,nzmax)
      call bilinear_horizontal_interpolation(vveta,y2,n,nzmax)
      call compute_standard_deviation(uueta,usigprof(n),n,n,nzmax)
      call compute_standard_deviation(vveta,vsigprof(n),n,n,nzmax)
    endif
    call bilinear_horizontal_interpolation(drhodzeta,rhograd1,n,nzmax)
    call bilinear_horizontal_interpolation(rhoeta,rho1,n,nzmax+1)

    call temporal_interpolation(y1(1),y1(2),uprof(n))
    call temporal_interpolation(y2(1),y2(2),vprof(n))
    call temporal_interpolation(rho1(1),rho1(2),rhoprof(n))
    call temporal_interpolation(rhograd1(1),rhograd1(2),rhogradprof(n))
  end do

  do n=indzeta,indzpeta
    call bilinear_horizontal_interpolation(wweta,y3,n,nzmax)
    call compute_standard_deviation(wweta,wsigprofeta(n),n,n,nzmax)
    call temporal_interpolation(y3(1),y3(2),wprofeta(n))
  end do
end subroutine interpol_all

subroutine interpol_misslev(n)
  !                            
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates u,v,w, density and density gradients.        *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
  !    Update: 2 March 1999                                                    *
  !                                                                            *
  !  Revision March 2005 by AST : all output variables in common block cal-    *
  !                               culation of standard deviation done in this  *
  !                               routine rather than subroutine call in order *
  !                               to save computation time                     *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! n                  level                                                   *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use hanna_mod

  implicit none

  ! Auxiliary variables needed for interpolation
  real :: y1(2),y2(2),y3(2),rho1(2),rhograd1(2)
  integer :: n


  !********************************************
  ! Multilinear interpolation in time and space
  !********************************************

  !**************************************
  ! 1.) Bilinear horizontal interpolation
  ! 2.) Temporal interpolation (linear)
  !**************************************

  call bilinear_horizontal_interpolation(ww,y3,n,nwzmax)
  call compute_standard_deviation(ww,wsigprof(n),n,n,nwzmax)

  indzindicator(n)=.false.

  if (ngrid.lt.0) then
    call bilinear_horizontal_interpolation(uupoleta,y1,n,nzmax)
    call bilinear_horizontal_interpolation(vvpoleta,y2,n,nzmax)
    call compute_standard_deviation(uupoleta,usigprof(n),n,n,nzmax)
    call compute_standard_deviation(vvpoleta,vsigprof(n),n,n,nzmax)
  else
    call bilinear_horizontal_interpolation(uueta,y1,n,nzmax+1)
    call bilinear_horizontal_interpolation(vveta,y2,n,nzmax+1)
    call compute_standard_deviation(uueta,usigprof(n),n,n,nzmax+1)
    call compute_standard_deviation(vveta,vsigprof(n),n,n,nzmax+1)
  endif

  call bilinear_horizontal_interpolation(wweta,y3,n,nzmax)
  call compute_standard_deviation(wweta,wsigprofeta(n),n,n,nzmax)
 
  call bilinear_horizontal_interpolation(drhodzeta,rhograd1,n,nzmax)
  call bilinear_horizontal_interpolation(rhoeta,rho1,n,nzmax+1)

  call temporal_interpolation(y1(1),y1(2),uprof(n))
  call temporal_interpolation(y2(1),y2(2),vprof(n))
  call temporal_interpolation(y3(1),y3(2),wprofeta(n))
  call temporal_interpolation(rho1(1),rho1(2),rhoprof(n))
  call temporal_interpolation(rhograd1(1),rhograd1(2),rhogradprof(n)) 
end subroutine interpol_misslev

subroutine interpol_wind(itime,xt,yt,zt,zteta,pp)
  !                           i   i  i  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates the wind data to current trajectory position.*
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
  ! u,v,w              wind components                                         *
  ! itime [s]          current temporal position                               *
  ! memtime(3) [s]     times of the wind fields in memory                      *
  ! xt,yt,zt           coordinates position for which wind data shall be       *
  !                    calculated                                              *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************


  implicit none

  integer :: itime,pp
  real :: xt,yt,zt
  real :: zteta

  ! Auxiliary variables needed for interpolation
  real :: dz1,dz2,dz
  real :: uh(2),vh(2),wh(2)
  ! ! Multilinear interpolation in time and space
  ! !********************************************

  ! Determine the lower left corner and its distance to the current position
  !*************************************************************************
  call find_grid_distances(xt,yt)

  ! Calculate variables for time interpolation
  !*******************************************
  call find_time_variables(itime)

  ! Determine the level below the current position for u,v
  !*******************************************************
  call find_z_level_meters(zt)

  ! Vertical distance to the level below and above current position
  !****************************************************************
  call find_vertical_variables(height,zt,indz,dz1,dz2,lbounds)

  !**********************************************************************
  ! 1.) Bilinear horizontal interpolation
  ! This has to be done separately for 6 fields (Temporal(2)*Vertical(3))
  !**********************************************************************

  ! Loop over 2 time steps and 2 levels
  !************************************
  call compute_standard_deviation(ww,wsig,indz,indz+1,nwzmax)
  call bilinear_spatial_interpolation(ww,wh,indz,dz1,dz2,nwzmax)
  call temporal_interpolation(wh(1),wh(2),w)

  ! Same for eta coordinates
  !*************************
  ! First the half levels
  !**********************
  call find_z_level_eta(zteta)
  call find_vertical_variables(uvheight,zteta,induv,dz1,dz2,lbounds_uv)

  if (ngrid.lt.0) then
    call compute_standard_deviation(uupoleta,usig,induv,induv+1,nzmax)
    call bilinear_spatial_interpolation(uupoleta,uh,induv,dz1,dz2,nzmax)
    call compute_standard_deviation(vvpoleta,vsig,induv,induv+1,nzmax)
    call bilinear_spatial_interpolation(vvpoleta,vh,induv,dz1,dz2,nzmax)
  else
    call compute_standard_deviation(uueta,usig,induv,induv+1,nzmax)
    call bilinear_spatial_interpolation(uueta,uh,induv,dz1,dz2,nzmax)
    call compute_standard_deviation(vveta,vsig,induv,induv+1,nzmax)
    call bilinear_spatial_interpolation(vveta,vh,induv,dz1,dz2,nzmax)
  endif
  call temporal_interpolation(uh(1),uh(2),u)
  call temporal_interpolation(vh(1),vh(2),v)

  ! Then for the model levels
  !**************************
  call find_vertical_variables(wheight,zteta,indzeta,dz1,dz2,lbounds_w)

  call compute_standard_deviation(wweta,wsigeta,indzeta,indzeta+1,nzmax)
  call bilinear_spatial_interpolation(wweta,wh,indzeta,dz1,dz2,nzmax)
  call temporal_interpolation(wh(1),wh(2),weta)
end subroutine interpol_wind

subroutine interpol_wind_short(itime,xt,yt,zt,zteta)
  !                                 i   i  i  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates the wind data to current trajectory position.*
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
  !                                                                            *
  !  Revision March 2005 by AST : all output variables in common block         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! u,v,w              wind components                                         *
  ! itime [s]          current temporal position                               *
  ! memtime(3) [s]     times of the wind fields in memory                      *
  ! xt,yt,zt           coordinates position for which wind data shall be       *
  !                    calculated                                              *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************


  implicit none

  integer, intent(in) :: itime
  real, intent(in) :: xt,yt,zt
  real, intent(in) :: zteta

  ! Auxiliary variables needed for interpolation
  real :: dz1,dz2,dz
  real :: uh(2),vh(2),wh(2)

  if (ngrid.gt.0) then
    call interpol_wind_short_nests(itime,xtn,ytn,zt)
    return
  endif
  !********************************************
  ! Multilinear interpolation in time and space
  !********************************************
  call find_grid_distances(xt,yt)

  ! Calculate variables for time interpolation
  !*******************************************
  call find_time_variables(itime)

  ! Determine the level below the current position for u,v
  !*******************************************************
  call find_z_level_meters(zt)
  call find_z_level_eta(zteta)

  ! Vertical distance to the level below and above current position
  !****************************************************************
  call find_vertical_variables(height,zt,indz,dz1,dz2,lbounds)

  !**********************************************************************
  ! 1.) Bilinear horizontal interpolation
  ! This has to be done separately for 6 fields (Temporal(2)*Vertical(3))
  !**********************************************************************

  ! Loop over 2 time steps and 2 levels
  !************************************
  call bilinear_spatial_interpolation(ww,wh,indz,dz1,dz2,nwzmax)
  call temporal_interpolation(wh(1),wh(2),w)

  ! Same for eta coordinates
  !*************************
  ! U,V level
  !**********
  call find_vertical_variables(uvheight,zteta,induv,dz1,dz2,lbounds_uv)

  if (ngrid.lt.0) then
    call bilinear_spatial_interpolation(uupoleta,uh,induv,dz1,dz2,nzmax)
    call bilinear_spatial_interpolation(vvpoleta,vh,induv,dz1,dz2,nzmax)
  else
    call bilinear_spatial_interpolation(uueta,uh,induv,dz1,dz2,nzmax)
    call bilinear_spatial_interpolation(vveta,vh,induv,dz1,dz2,nzmax)
  endif
  call temporal_interpolation(uh(1),uh(2),u)
  call temporal_interpolation(vh(1),vh(2),v)

  ! W level
  !**********
  call find_vertical_variables(wheight,zteta,indzeta,dz1,dz2,lbounds_w)
  call bilinear_spatial_interpolation(wweta,wh,indzeta,dz1,dz2,nzmax)
  call temporal_interpolation(wh(1),wh(2),weta)
end subroutine interpol_wind_short

subroutine interpol_rain(yy1,yy2,yy3,nxmax,nymax,nzmax,nx, &
     ny,iwftouse,xt,yt,level,itime1,itime2,itime,yint1,yint2,yint3)
  !                          i   i   i    i    i     i   i
  !i    i    i  i    i     i      i      i     o     o     o
  !****************************************************************************
  !                                                                           *
  !  Interpolation of meteorological fields on 2-d model layers.              *
  !  In horizontal direction bilinear interpolation interpolation is used.    *
  !  Temporally a linear interpolation is used.                               *
  !  Three fields are interpolated at the same time.                          *
  !                                                                           *
  !  This is a special version of levlininterpol to save CPU time.            *
  !                                                                           *
  !  1 first time                                                             *
  !  2 second time                                                            *
  !                                                                           *
  !                                                                           *
  !     Author: A. Stohl                                                      *
  !                                                                           *
  !     30 August 1996                                                        *
  !                                                                           *
  !****************************************************************************
  !                                                                           *
  ! Variables:                                                                *
  !                                                                           *
  ! dt1,dt2              time differences between fields and current position *
  ! dz1,dz2              z distance between levels and current position       *
  ! height(nzmax)        heights of the model levels                          *
  ! indexh               help variable                                        *
  ! indz                 the level closest to the current trajectory position *
  ! indzh                help variable                                        *
  ! itime                current time                                         *
  ! itime1               time of the first wind field                         *
  ! itime2               time of the second wind field                        *
  ! ix,jy                x,y coordinates of lower left subgrid point          *
  ! level                level at which interpolation shall be done           *
  ! iwftouse             points to the place of the wind field                *
  ! nx,ny                actual field dimensions in x,y and z direction       *
  ! nxmax,nymax,nzmax    maximum field dimensions in x,y and z direction      *
  ! xt                   current x coordinate                                 *
  ! yint                 the final interpolated value                         *
  ! yt                   current y coordinate                                 *
  ! yy(0:nxmax,0:nymax,nzmax,3) meteorological field used for interpolation   *
  ! zt                   current z coordinate                                 *
  !                                                                           *
  !****************************************************************************
  use par_mod, only: numwfmem

  implicit none

  integer :: nx,ny,nxmax,nymax,nzmax,memind(numwfmem),m,ix,jy,ixp,jyp
  integer :: itime,itime1,itime2,level,indexh
  real :: yy1(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: yy2(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: yy3(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: ddx,ddy,rddx,rddy,dt1,dt2,dt,y1(2),y2(2),y3(2)
  real :: xt,yt,yint1,yint2,yint3,p1,p2,p3,p4
  integer :: iwftouse



  ! If point at border of grid -> small displacement into grid
  !***********************************************************

  if (xt.ge.real(nx-1)) xt=real(nx-1)-0.00001
  if (yt.ge.real(ny-1)) yt=real(ny-1)-0.00001



  !**********************************************************************
  ! 1.) Bilinear horizontal interpolation
  ! This has to be done separately for 2 fields (Temporal)
  !*******************************************************

  ! Determine the lower left corner and its distance to the current position
  !*************************************************************************

  ix=int(xt)
  jy=int(yt)
  ixp=ix+1
  jyp=jy+1
  ddx=xt-real(ix)
  ddy=yt-real(jy)
  rddx=1.-ddx
  rddy=1.-ddy
  p1=rddx*rddy
  p2=ddx*rddy
  p3=rddx*ddy
  p4=ddx*ddy


  ! Loop over 2 time steps
  !***********************

  !  do m=1,2
  indexh=iwftouse

  y1(1)=p1*yy1(ix ,jy ,level,indexh) &
     + p2*yy1(ixp,jy ,level,indexh) &
     + p3*yy1(ix ,jyp,level,indexh) &
     + p4*yy1(ixp,jyp,level,indexh)
  y2(1)=p1*yy2(ix ,jy ,level,indexh) &
     + p2*yy2(ixp,jy ,level,indexh) &
     + p3*yy2(ix ,jyp,level,indexh) &
     + p4*yy2(ixp,jyp,level,indexh)
  y3(1)=p1*yy3(ix ,jy ,level,indexh) &
     + p2*yy3(ixp,jy ,level,indexh) &
     + p3*yy3(ix ,jyp,level,indexh) &
     + p4*yy3(ixp,jyp,level,indexh)
  !  end do


  !************************************
  ! 2.) Temporal interpolation (linear) - skip to be consistent with clouds
  !************************************

  !  dt1=real(itime-itime1)
  !  dt2=real(itime2-itime)
  !  dt=dt1+dt2

  !  yint1=(y1(1)*dt2+y1(2)*dt1)/dt
  !  yint2=(y2(1)*dt2+y2(2)*dt1)/dt
  !  yint3=(y3(1)*dt2+y3(2)*dt1)/dt

   yint1=y1(1)
   yint2=y2(1)
   yint3=y3(1)
end subroutine interpol_rain

subroutine interpol_vdep(level,vdepo)
  !                           i     o
  !****************************************************************************
  !                                                                           *
  !  Interpolation of the deposition velocity on 2-d model layer.             *
  !  In horizontal direction bilinear interpolation interpolation is used.    *
  !  Temporally a linear interpolation is used.                               *
  !                                                                           *
  !  1 first time                                                             *
  !  2 second time                                                            *
  !                                                                           *
  !                                                                           *
  !     Author: A. Stohl                                                      *
  !                                                                           *
  !     30 May 1994                                                           *
  !                                                                           *
  !****************************************************************************
  !                                                                           *
  ! Variables:                                                                *
  !                                                                           *
  ! level                number of species for which interpolation is done    *
  !                                                                           *
  !****************************************************************************

  implicit none

  integer :: level,indexh,m
  real :: y(2),vdepo

  ! a) Bilinear horizontal interpolation
  do m=1,2
    indexh=memind(m)

    y(m)=p1*vdep(ix ,jy ,level,indexh) &
         +p2*vdep(ixp,jy ,level,indexh) &
         +p3*vdep(ix ,jyp,level,indexh) &
         +p4*vdep(ixp,jyp,level,indexh)
  end do



  ! b) Temporal interpolation

  vdepo=(y(1)*dt2+y(2)*dt1)*dtt

  depoindicator(level)=.false.
end subroutine interpol_vdep

subroutine interpol_all_nests(itime,xt,yt,zt)
  !                                i   i  i  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates everything that is needed for calculating the*
  !  dispersion.                                                               *
  !  Version for interpolating nested grids.                                   *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    9 February 1999                                                         *
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
  !                    calculated                                              *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use hanna_mod

  implicit none

  integer :: itime
  real :: xt,yt,zt

  ! Auxiliary variables needed for interpolation
  real :: ust1(2),wst1(2),oli1(2),oliaux
  real :: y1(2),y2(2),y3(2),rho1(2),rhograd1(2)
  real :: usl,vsl,wsl,usq,vsq,wsq,xaux
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

    ust1(m)=p1*ustarn(ix ,jy ,1,indexh,ngrid) &
         + p2*ustarn(ixp,jy ,1,indexh,ngrid) &
         + p3*ustarn(ix ,jyp,1,indexh,ngrid) &
         + p4*ustarn(ixp,jyp,1,indexh,ngrid)
    wst1(m)=p1*wstarn(ix ,jy ,1,indexh,ngrid) &
         + p2*wstarn(ixp,jy ,1,indexh,ngrid) &
         + p3*wstarn(ix ,jyp,1,indexh,ngrid) &
         + p4*wstarn(ixp,jyp,1,indexh,ngrid)
    oli1(m)=p1*olin(ix ,jy ,1,indexh,ngrid) &
         + p2*olin(ixp,jy ,1,indexh,ngrid) &
         + p3*olin(ix ,jyp,1,indexh,ngrid) &
         + p4*olin(ixp,jyp,1,indexh,ngrid)
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

  do n=indz,indz+1
    usl=0.
    vsl=0.
    wsl=0.
    usq=0.
    vsq=0.
    wsq=0.
    do m=1,2
      indexh=memind(m)
      y1(m)=p1*uun(ix ,jy ,n,indexh,ngrid) &
           +p2*uun(ixp,jy ,n,indexh,ngrid) &
           +p3*uun(ix ,jyp,n,indexh,ngrid) &
           +p4*uun(ixp,jyp,n,indexh,ngrid)
      y2(m)=p1*vvn(ix ,jy ,n,indexh,ngrid) &
           +p2*vvn(ixp,jy ,n,indexh,ngrid) &
           +p3*vvn(ix ,jyp,n,indexh,ngrid) &
           +p4*vvn(ixp,jyp,n,indexh,ngrid)
      y3(m)=p1*wwn(ix ,jy ,n,indexh,ngrid) &
           +p2*wwn(ixp,jy ,n,indexh,ngrid) &
           +p3*wwn(ix ,jyp,n,indexh,ngrid) &
           +p4*wwn(ixp,jyp,n,indexh,ngrid)
      rhograd1(m)=p1*drhodzn(ix ,jy ,n,indexh,ngrid) &
           +p2*drhodzn(ixp,jy ,n,indexh,ngrid) &
           +p3*drhodzn(ix ,jyp,n,indexh,ngrid) &
           +p4*drhodzn(ixp,jyp,n,indexh,ngrid)
      rho1(m)=p1*rhon(ix ,jy ,n,indexh,ngrid) &
           +p2*rhon(ixp,jy ,n,indexh,ngrid) &
           +p3*rhon(ix ,jyp,n,indexh,ngrid) &
           +p4*rhon(ixp,jyp,n,indexh,ngrid)

     usl=usl+uun(ix ,jy ,n,indexh,ngrid)+uun(ixp,jy ,n,indexh,ngrid) &
          +uun(ix ,jyp,n,indexh,ngrid)+uun(ixp,jyp,n,indexh,ngrid)
     vsl=vsl+vvn(ix ,jy ,n,indexh,ngrid)+vvn(ixp,jy ,n,indexh,ngrid) &
          +vvn(ix ,jyp,n,indexh,ngrid)+vvn(ixp,jyp,n,indexh,ngrid)
     wsl=wsl+wwn(ix ,jy ,n,indexh,ngrid)+wwn(ixp,jy ,n,indexh,ngrid) &
          +wwn(ix ,jyp,n,indexh,ngrid)+wwn(ixp,jyp,n,indexh,ngrid)

    usq=usq+uun(ix ,jy ,n,indexh,ngrid)*uun(ix ,jy ,n,indexh,ngrid)+ &
         uun(ixp,jy ,n,indexh,ngrid)*uun(ixp,jy ,n,indexh,ngrid)+ &
         uun(ix ,jyp,n,indexh,ngrid)*uun(ix ,jyp,n,indexh,ngrid)+ &
         uun(ixp,jyp,n,indexh,ngrid)*uun(ixp,jyp,n,indexh,ngrid)
    vsq=vsq+vvn(ix ,jy ,n,indexh,ngrid)*vvn(ix ,jy ,n,indexh,ngrid)+ &
         vvn(ixp,jy ,n,indexh,ngrid)*vvn(ixp,jy ,n,indexh,ngrid)+ &
         vvn(ix ,jyp,n,indexh,ngrid)*vvn(ix ,jyp,n,indexh,ngrid)+ &
         vvn(ixp,jyp,n,indexh,ngrid)*vvn(ixp,jyp,n,indexh,ngrid)
    wsq=wsq+wwn(ix ,jy ,n,indexh,ngrid)*wwn(ix ,jy ,n,indexh,ngrid)+ &
         wwn(ixp,jy ,n,indexh,ngrid)*wwn(ixp,jy ,n,indexh,ngrid)+ &
         wwn(ix ,jyp,n,indexh,ngrid)*wwn(ix ,jyp,n,indexh,ngrid)+ &
         wwn(ixp,jyp,n,indexh,ngrid)*wwn(ixp,jyp,n,indexh,ngrid)
    end do
    uprof(n)=(y1(1)*dt2+y1(2)*dt1)*dtt
    vprof(n)=(y2(1)*dt2+y2(2)*dt1)*dtt
    wprof(n)=(y3(1)*dt2+y3(2)*dt1)*dtt
    rhoprof(n)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
    rhogradprof(n)=(rhograd1(1)*dt2+rhograd1(2)*dt1)*dtt
    indzindicator(n)=.false.

  ! Compute standard deviations
  !****************************

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


    xaux=wsq-wsl*wsl/8.
    if (xaux.lt.eps) then
      wsigprof(n)=0.
    else
      wsigprof(n)=sqrt(xaux/7.)
    endif

  end do
end subroutine interpol_all_nests

subroutine interpol_misslev_nests(n)
  !                                  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates u,v,w, density and density gradients.        *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! n                  level                                                   *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use hanna_mod

  implicit none

  ! Auxiliary variables needed for interpolation
  real :: y1(2),y2(2),y3(2),rho1(2),rhograd1(2)
  real :: usl,vsl,wsl,usq,vsq,wsq,xaux
  integer :: m,n,indexh
  real,parameter :: eps=1.0e-30


  !********************************************
  ! Multilinear interpolation in time and space
  !********************************************


  !**************************************
  ! 1.) Bilinear horizontal interpolation
  ! 2.) Temporal interpolation (linear)
  !**************************************

  ! Loop over 2 time steps
  !***********************

  usl=0.
  vsl=0.
  wsl=0.
  usq=0.
  vsq=0.
  wsq=0.
  do m=1,2
    indexh=memind(m)
    y1(m)=p1*uun(ix ,jy ,n,indexh,ngrid) &
         +p2*uun(ixp,jy ,n,indexh,ngrid) &
         +p3*uun(ix ,jyp,n,indexh,ngrid) &
         +p4*uun(ixp,jyp,n,indexh,ngrid)
    y2(m)=p1*vvn(ix ,jy ,n,indexh,ngrid) &
         +p2*vvn(ixp,jy ,n,indexh,ngrid) &
         +p3*vvn(ix ,jyp,n,indexh,ngrid) &
         +p4*vvn(ixp,jyp,n,indexh,ngrid)
    y3(m)=p1*wwn(ix ,jy ,n,indexh,ngrid) &
         +p2*wwn(ixp,jy ,n,indexh,ngrid) &
         +p3*wwn(ix ,jyp,n,indexh,ngrid) &
         +p4*wwn(ixp,jyp,n,indexh,ngrid)
    rho1(m)=p1*rhon(ix ,jy ,n,indexh,ngrid) &
         +p2*rhon(ixp,jy ,n,indexh,ngrid) &
         +p3*rhon(ix ,jyp,n,indexh,ngrid) &
         +p4*rhon(ixp,jyp,n,indexh,ngrid)
    rhograd1(m)=p1*drhodzn(ix ,jy ,n,indexh,ngrid) &
         +p2*drhodzn(ixp,jy ,n,indexh,ngrid) &
         +p3*drhodzn(ix ,jyp,n,indexh,ngrid) &
         +p4*drhodzn(ixp,jyp,n,indexh,ngrid)

     usl=usl+uun(ix ,jy ,n,indexh,ngrid)+uun(ixp,jy ,n,indexh,ngrid) &
          +uun(ix ,jyp,n,indexh,ngrid)+uun(ixp,jyp,n,indexh,ngrid)
     vsl=vsl+vvn(ix ,jy ,n,indexh,ngrid)+vvn(ixp,jy ,n,indexh,ngrid) &
          +vvn(ix ,jyp,n,indexh,ngrid)+vvn(ixp,jyp,n,indexh,ngrid)
     wsl=wsl+wwn(ix ,jy ,n,indexh,ngrid)+wwn(ixp,jy ,n,indexh,ngrid) &
          +wwn(ix ,jyp,n,indexh,ngrid)+wwn(ixp,jyp,n,indexh,ngrid)

    usq=usq+uun(ix ,jy ,n,indexh,ngrid)*uun(ix ,jy ,n,indexh,ngrid)+ &
         uun(ixp,jy ,n,indexh,ngrid)*uun(ixp,jy ,n,indexh,ngrid)+ &
         uun(ix ,jyp,n,indexh,ngrid)*uun(ix ,jyp,n,indexh,ngrid)+ &
         uun(ixp,jyp,n,indexh,ngrid)*uun(ixp,jyp,n,indexh,ngrid)
    vsq=vsq+vvn(ix ,jy ,n,indexh,ngrid)*vvn(ix ,jy ,n,indexh,ngrid)+ &
         vvn(ixp,jy ,n,indexh,ngrid)*vvn(ixp,jy ,n,indexh,ngrid)+ &
         vvn(ix ,jyp,n,indexh,ngrid)*vvn(ix ,jyp,n,indexh,ngrid)+ &
         vvn(ixp,jyp,n,indexh,ngrid)*vvn(ixp,jyp,n,indexh,ngrid)
    wsq=wsq+wwn(ix ,jy ,n,indexh,ngrid)*wwn(ix ,jy ,n,indexh,ngrid)+ &
         wwn(ixp,jy ,n,indexh,ngrid)*wwn(ixp,jy ,n,indexh,ngrid)+ &
         wwn(ix ,jyp,n,indexh,ngrid)*wwn(ix ,jyp,n,indexh,ngrid)+ &
         wwn(ixp,jyp,n,indexh,ngrid)*wwn(ixp,jyp,n,indexh,ngrid)
  end do
  uprof(n)=(y1(1)*dt2+y1(2)*dt1)*dtt
  vprof(n)=(y2(1)*dt2+y2(2)*dt1)*dtt
  wprof(n)=(y3(1)*dt2+y3(2)*dt1)*dtt
  rhoprof(n)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
  rhogradprof(n)=(rhograd1(1)*dt2+rhograd1(2)*dt1)*dtt
  indzindicator(n)=.false.

  ! Compute standard deviations
  !****************************

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


  xaux=wsq-wsl*wsl/8.
  if (xaux.lt.eps) then
    wsigprof(n)=0.
  else
    wsigprof(n)=sqrt(xaux/7.)
  endif
end subroutine interpol_misslev_nests

subroutine interpol_wind_short_nests(itime,xt,yt,zt)
  !                                       i   i  i  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates the wind data to current trajectory position.*
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
  !                                                                            *
  !  Revision March 2005 by AST : all output variables in common block         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! u,v,w              wind components                                         *
  ! itime [s]          current temporal position                               *
  ! memtime(3) [s]     times of the wind fields in memory                      *
  ! xt,yt,zt           coordinates position for which wind data shall be       *
  !                    calculated                                              *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: itime
  real :: xt,yt,zt

  ! Auxiliary variables needed for interpolation
  real :: dz1,dz2,dz
  real :: u1(2),v1(2),w1(2),uh(2),vh(2),wh(2)
  integer :: i,m,n,indexh,indzh


  !********************************************
  ! Multilinear interpolation in time and space
  !********************************************

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

  ! Determine the level below the current position for u,v
  !*******************************************************

  do i=2,nz
    if (height(i).gt.zt) then
      indz=i-1
      exit
    endif
  end do

  ! Vertical distance to the level below and above current position
  !****************************************************************

  dz=1./(height(indz+1)-height(indz))
  dz1=(zt-height(indz))*dz
  dz2=(height(indz+1)-zt)*dz


  !**********************************************************************
  ! 1.) Bilinear horizontal interpolation
  ! This has to be done separately for 6 fields (Temporal(2)*Vertical(3))
  !**********************************************************************

  ! Loop over 2 time steps and 2 levels
  !************************************

  do m=1,2
    indexh=memind(m)
    do n=1,2
      indzh=indz+n-1

      u1(n)=p1*uun(ix ,jy ,indzh,indexh,ngrid) &
           +p2*uun(ixp,jy ,indzh,indexh,ngrid) &
           +p3*uun(ix ,jyp,indzh,indexh,ngrid) &
           +p4*uun(ixp,jyp,indzh,indexh,ngrid)
      v1(n)=p1*vvn(ix ,jy ,indzh,indexh,ngrid) &
           +p2*vvn(ixp,jy ,indzh,indexh,ngrid) &
           +p3*vvn(ix ,jyp,indzh,indexh,ngrid) &
           +p4*vvn(ixp,jyp,indzh,indexh,ngrid)
      w1(n)=p1*wwn(ix ,jy ,indzh,indexh,ngrid) &
           +p2*wwn(ixp,jy ,indzh,indexh,ngrid) &
           +p3*wwn(ix ,jyp,indzh,indexh,ngrid) &
           +p4*wwn(ixp,jyp,indzh,indexh,ngrid)

    end do


  !**********************************
  ! 2.) Linear vertical interpolation
  !**********************************

    uh(m)=dz2*u1(1)+dz1*u1(2)
    vh(m)=dz2*v1(1)+dz1*v1(2)
    wh(m)=dz2*w1(1)+dz1*w1(2)
  end do


  !************************************
  ! 3.) Temporal interpolation (linear)
  !************************************

  u=(uh(1)*dt2+uh(2)*dt1)*dtt
  v=(vh(1)*dt2+vh(2)*dt1)*dtt
  w=(wh(1)*dt2+wh(2)*dt1)*dtt
end subroutine interpol_wind_short_nests

subroutine interpol_wind_nests(itime,xt,yt,zt)
  !                                 i   i  i  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates the wind data to current trajectory position.*
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
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
  ! u,v,w              wind components                                         *
  ! itime [s]          current temporal position                               *
  ! memtime(3) [s]     times of the wind fields in memory                      *
  ! xt,yt,zt           coordinates position for which wind data shall be       *
  !                    calculated                                              *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: itime
  real :: xt,yt,zt

  ! Auxiliary variables needed for interpolation
  real :: dz1,dz2,dz
  real :: u1(2),v1(2),w1(2),uh(2),vh(2),wh(2)
  real :: usl,vsl,wsl,usq,vsq,wsq,xaux
  integer :: i,m,n,indexh,indzh
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

  ! Determine the level below the current position for u,v
  !*******************************************************

  do i=2,nz
    if (height(i).gt.zt) then
      indz=i-1
      exit
    endif
  end do

  ! Vertical distance to the level below and above current position
  !****************************************************************

  dz=1./(height(indz+1)-height(indz))
  dz1=(zt-height(indz))*dz
  dz2=(height(indz+1)-zt)*dz


  !**********************************************************************
  ! 1.) Bilinear horizontal interpolation
  ! This has to be done separately for 6 fields (Temporal(2)*Vertical(3))
  !**********************************************************************

  ! Loop over 2 time steps and 2 levels
  !************************************

  usl=0.
  vsl=0.
  wsl=0.
  usq=0.
  vsq=0.
  wsq=0.
  do m=1,2
    indexh=memind(m)
    do n=1,2
      indzh=indz+n-1

      u1(n)=p1*uun(ix ,jy ,indzh,indexh,ngrid) &
           +p2*uun(ixp,jy ,indzh,indexh,ngrid) &
           +p3*uun(ix ,jyp,indzh,indexh,ngrid) &
           +p4*uun(ixp,jyp,indzh,indexh,ngrid)
      v1(n)=p1*vvn(ix ,jy ,indzh,indexh,ngrid) &
           +p2*vvn(ixp,jy ,indzh,indexh,ngrid) &
           +p3*vvn(ix ,jyp,indzh,indexh,ngrid) &
           +p4*vvn(ixp,jyp,indzh,indexh,ngrid)
      w1(n)=p1*wwn(ix ,jy ,indzh,indexh,ngrid) &
           +p2*wwn(ixp,jy ,indzh,indexh,ngrid) &
           +p3*wwn(ix ,jyp,indzh,indexh,ngrid) &
           +p4*wwn(ixp,jyp,indzh,indexh,ngrid)

      usl=usl+uun(ix ,jy ,indzh,indexh,ngrid)+ &
           uun(ixp,jy ,indzh,indexh,ngrid) &
           +uun(ix ,jyp,indzh,indexh,ngrid)+ &
           uun(ixp,jyp,indzh,indexh,ngrid)
      vsl=vsl+vvn(ix ,jy ,indzh,indexh,ngrid)+ &
           vvn(ixp,jy ,indzh,indexh,ngrid) &
           +vvn(ix ,jyp,indzh,indexh,ngrid)+ &
           vvn(ixp,jyp,indzh,indexh,ngrid)
      wsl=wsl+wwn(ix ,jy ,indzh,indexh,ngrid)+ &
           wwn(ixp,jy ,indzh,indexh,ngrid) &
           +wwn(ix ,jyp,indzh,indexh,ngrid)+ &
           wwn(ixp,jyp,indzh,indexh,ngrid)

      usq=usq+uun(ix ,jy ,indzh,indexh,ngrid)* &
           uun(ix ,jy ,indzh,indexh,ngrid)+ &
           uun(ixp,jy ,indzh,indexh,ngrid)*uun(ixp,jy ,indzh,indexh,ngrid)+ &
           uun(ix ,jyp,indzh,indexh,ngrid)*uun(ix ,jyp,indzh,indexh,ngrid)+ &
           uun(ixp,jyp,indzh,indexh,ngrid)*uun(ixp,jyp,indzh,indexh,ngrid)
      vsq=vsq+vvn(ix ,jy ,indzh,indexh,ngrid)* &
           vvn(ix ,jy ,indzh,indexh,ngrid)+ &
           vvn(ixp,jy ,indzh,indexh,ngrid)*vvn(ixp,jy ,indzh,indexh,ngrid)+ &
           vvn(ix ,jyp,indzh,indexh,ngrid)*vvn(ix ,jyp,indzh,indexh,ngrid)+ &
           vvn(ixp,jyp,indzh,indexh,ngrid)*vvn(ixp,jyp,indzh,indexh,ngrid)
      wsq=wsq+wwn(ix ,jy ,indzh,indexh,ngrid)* &
           wwn(ix ,jy ,indzh,indexh,ngrid)+ &
           wwn(ixp,jy ,indzh,indexh,ngrid)*wwn(ixp,jy ,indzh,indexh,ngrid)+ &
           wwn(ix ,jyp,indzh,indexh,ngrid)*wwn(ix ,jyp,indzh,indexh,ngrid)+ &
           wwn(ixp,jyp,indzh,indexh,ngrid)*wwn(ixp,jyp,indzh,indexh,ngrid)
    end do


  !**********************************
  ! 2.) Linear vertical interpolation
  !**********************************

    uh(m)=dz2*u1(1)+dz1*u1(2)
    vh(m)=dz2*v1(1)+dz1*v1(2)
    wh(m)=dz2*w1(1)+dz1*w1(2)
  end do


  !************************************
  ! 3.) Temporal interpolation (linear)
  !************************************

  u=(uh(1)*dt2+uh(2)*dt1)*dtt
  v=(vh(1)*dt2+vh(2)*dt1)*dtt
  w=(wh(1)*dt2+wh(2)*dt1)*dtt


  ! Compute standard deviations
  !****************************

  xaux=usq-usl*usl/16.
  if (xaux.lt.eps) then
    usig=0.
  else
    usig=sqrt(xaux/15.)
  endif

  xaux=vsq-vsl*vsl/16.
  if (xaux.lt.eps) then
    vsig=0.
  else
    vsig=sqrt(xaux/15.)
  endif


  xaux=wsq-wsl*wsl/16.
  if (xaux.lt.eps) then
    wsig=0.
  else
    wsig=sqrt(xaux/15.)
  endif
end subroutine interpol_wind_nests

subroutine interpol_rain_nests(yy1,yy2,yy3,nxmaxn,nymaxn,nzmax, &
       maxnests,ngrid,nxn,nyn,iwftouse,xt,yt,level,itime1,itime2,itime, &
       yint1,yint2,yint3)
  !                                i   i   i    i      i      i
  !   i       i    i   i    i    i  i    i     i      i      i
  !  o     o     o
  !****************************************************************************
  !                                                                           *
  !  Interpolation of meteorological fields on 2-d model layers for nested    *
  !  grids. This routine is related to levlin3interpol.f for the mother domain*
  !                                                                           *
  !  In horizontal direction bilinear interpolation interpolation is used.    *
  !  Temporally a linear interpolation is used.                               *
  !  Three fields are interpolated at the same time.                          *
  !                                                                           *
  !  This is a special version of levlininterpol to save CPU time.            *
  !                                                                           *
  !  1 first time                                                             *
  !  2 second time                                                            *
  !                                                                           *
  !                                                                           *
  !     Author: A. Stohl                                                      *
  !                                                                           *
  !     15 March 2000                                                         *
  !                                                                           *
  !****************************************************************************
  !                                                                           *
  ! Variables:                                                                *
  !                                                                           *
  ! dt1,dt2              time differences between fields and current position *
  ! dz1,dz2              z distance between levels and current position       *
  ! height(nzmax)        heights of the model levels                          *
  ! indexh               help variable                                        *
  ! indz                 the level closest to the current trajectory position *
  ! indzh                help variable                                        *
  ! itime                current time                                         *
  ! itime1               time of the first wind field                         *
  ! itime2               time of the second wind field                        *
  ! ix,jy                x,y coordinates of lower left subgrid point          *
  ! level                level at which interpolation shall be done           *
  ! iwftouse             points to the place of the wind field                *
  ! nx,ny                actual field dimensions in x,y and z direction       *
  ! nxmax,nymax,nzmax    maximum field dimensions in x,y and z direction      *
  ! xt                   current x coordinate                                 *
  ! yint                 the final interpolated value                         *
  ! yt                   current y coordinate                                 *
  ! yy(0:nxmax,0:nymax,nzmax,3) meteorological field used for interpolation   *
  ! zt                   current z coordinate                                 *
  !                                                                           *
  !****************************************************************************
  use par_mod, only: numwfmem

  implicit none

  integer :: maxnests,ngrid
  integer :: nxn(maxnests),nyn(maxnests),nxmaxn,nymaxn,nzmax,iwftouse
  integer :: m,ix,jy,ixp,jyp,itime,itime1,itime2,level,indexh
  real :: yy1(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,maxnests)
  real :: yy2(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,maxnests)
  real :: yy3(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,maxnests)
  real :: ddx,ddy,rddx,rddy,dt1,dt2,dt,y1(2),y2(2),y3(2)
  real :: xt,yt,yint1,yint2,yint3,p1,p2,p3,p4



  ! If point at border of grid -> small displacement into grid
  !***********************************************************

  ! if (xt.ge.(real(nxn(ngrid)-1)-0.0001)) &
  !      xt=real(nxn(ngrid)-1)-0.0001
  ! if (yt.ge.(real(nyn(ngrid)-1)-0.0001)) &
  !      yt=real(nyn(ngrid)-1)-0.0001

  ! ESO make it consistent with interpol_rain
  if (xt.ge.(real(nxn(ngrid)-1))) xt=real(nxn(ngrid)-1)-0.00001
  if (yt.ge.(real(nyn(ngrid)-1))) yt=real(nyn(ngrid)-1)-0.00001



  !**********************************************************************
  ! 1.) Bilinear horizontal interpolation
  ! This has to be done separately for 2 fields (Temporal)
  !*******************************************************

  ! Determine the lower left corner and its distance to the current position
  !*************************************************************************

  ix=int(xt)
  jy=int(yt)

  ixp=ix+1
  jyp=jy+1
  ddx=xt-real(ix)
  ddy=yt-real(jy)
  rddx=1.-ddx
  rddy=1.-ddy
  p1=rddx*rddy
  p2=ddx*rddy
  p3=rddx*ddy
  p4=ddx*ddy


  ! Loop over 2 time steps
  !***********************

  !  do m=1,2
  !    indexh=memind(m)
    indexh=iwftouse

    y1(1)=p1*yy1(ix ,jy ,level,indexh,ngrid) &
         + p2*yy1(ixp,jy ,level,indexh,ngrid) &
         + p3*yy1(ix ,jyp,level,indexh,ngrid) &
         + p4*yy1(ixp,jyp,level,indexh,ngrid)
    y2(1)=p1*yy2(ix ,jy ,level,indexh,ngrid) &
         + p2*yy2(ixp,jy ,level,indexh,ngrid) &
         + p3*yy2(ix ,jyp,level,indexh,ngrid) &
         + p4*yy2(ixp,jyp,level,indexh,ngrid)
    y3(1)=p1*yy3(ix ,jy ,level,indexh,ngrid) &
         + p2*yy3(ixp,jy ,level,indexh,ngrid) &
         + p3*yy3(ix ,jyp,level,indexh,ngrid) &
         + p4*yy3(ixp,jyp,level,indexh,ngrid)
  !  end do


  !************************************
  ! 2.) Temporal interpolation (linear)
  !************************************

  ! dt1=real(itime-itime1)
  ! dt2=real(itime2-itime)
  ! dt=dt1+dt2

  ! yint1=(y1(1)*dt2+y1(2)*dt1)/dt
  ! yint2=(y2(1)*dt2+y2(2)*dt1)/dt
  ! yint3=(y3(1)*dt2+y3(2)*dt1)/dt

   yint1=y1(1)
   yint2=y2(1)
   yint3=y3(1)
end subroutine interpol_rain_nests

subroutine interpol_vdep_nests(level,vdepo)
  !                                 i     o
  !****************************************************************************
  !                                                                           *
  !  Interpolation of the deposition velocity on 2-d model layer.             *
  !  In horizontal direction bilinear interpolation interpolation is used.    *
  !  Temporally a linear interpolation is used.                               *
  !                                                                           *
  !  1 first time                                                             *
  !  2 second time                                                            *
  !                                                                           *
  !                                                                           *
  !     Author: A. Stohl                                                      *
  !                                                                           *
  !     30 May 1994                                                           *
  !                                                                           *
  !****************************************************************************
  !                                                                           *
  ! Variables:                                                                *
  !                                                                           *
  ! level                number of species for which interpolation is done    *
  !                                                                           *
  !****************************************************************************


  implicit none

  integer :: level,indexh,m
  real :: y(2),vdepo

  ! a) Bilinear horizontal interpolation

  do m=1,2
    indexh=memind(m)

    y(m)=p1*vdepn(ix ,jy ,level,indexh,ngrid) &
         +p2*vdepn(ixp,jy ,level,indexh,ngrid) &
         +p3*vdepn(ix ,jyp,level,indexh,ngrid) &
         +p4*vdepn(ixp,jyp,level,indexh,ngrid)
  end do


  ! b) Temporal interpolation

  vdepo=(y(1)*dt2+y(2)*dt1)*dtt

  depoindicator(level)=.false.
end subroutine interpol_vdep_nests

end module interpol_mod



