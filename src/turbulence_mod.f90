module turbulence_mod
  use par_mod
  use com_mod
  use interpol_mod
  use hanna_mod
  use random_mod, only: ran3
  use coordinates_ecmwf
  use particle_mod

  implicit none

contains

subroutine turbulence_boundarylayer(ipart,nrand,dt,zts,rhoa,rhograd)

  implicit none 

  integer, intent(in) ::          &
    ipart                           ! particle index
  integer, intent(inout) ::       &
    nrand                           ! random number used for turbulence
  real,intent(in) ::              &
    dt,                           & ! real(ldt)
    rhoa,                         & ! air density, used in CBL
    rhograd                         ! vertical gradient of the air density, used in CBL
  real,intent(inout) ::           &
    zts                             ! local 'real' copy of the particle position
  real ::                         &
    delz,                         & ! change in vertical position due to turbulence
    ru,rv,rw,                     & ! used for computing turbulence
    dtf,rhoaux,dtftlw,ath,bth,    & ! CBL related
    ptot_lhh,Q_lhh,phi_lhh,       & ! CBL related
    old_wp_buf,dcas,dcas1,        & ! CBL related
    del_test                        ! CBL related
  integer ::                      &
    flagrein,                     & ! flag used in CBL scheme
    i                               ! ĺoop variable


    if (turbswitch) then
      call hanna(zts)
    else
      call hanna1(zts)
    endif

  !*****************************************
  ! Determine the new diffusivity velocities
  !*****************************************

  ! Horizontal components
  !**********************
  ! sigu,sigv, tlv and tlu are defined in hanna_mod, would be better to move below there
    if (nrand+1.gt.maxrand) nrand=1
    if (dt/tlu.lt..5) then
      part(ipart)%turbvel%u=(1.-dt/tlu)*part(ipart)%turbvel%u+rannumb(nrand)*sigu*sqrt(2.*dt/tlu)
    else
      ru=exp(-dt/tlu)
      part(ipart)%turbvel%u=ru*part(ipart)%turbvel%u+rannumb(nrand)*sigu*sqrt(1.-ru**2)
    endif
    if (dt/tlv.lt..5) then
      part(ipart)%turbvel%v=(1.-dt/tlv)*part(ipart)%turbvel%v+rannumb(nrand+1)*sigv*sqrt(2.*dt/tlv)
    else
      rv=exp(-dt/tlv)
      part(ipart)%turbvel%v=rv*part(ipart)%turbvel%v+rannumb(nrand+1)*sigv*sqrt(1.-rv**2)
    endif
    nrand=nrand+2


    if (nrand+ifine.gt.maxrand) nrand=1
    rhoaux=rhograd/rhoa
    dtf=dt*fine

    dtftlw=dtf/tlw

  ! Loop over ifine short time steps for vertical component
  !********************************************************
  ! tlw,dsigwdz and dsigw2dz is defined in hanna_mod, maybe move some below there
    do i=1,ifine

  ! Determine the drift velocity and density correction velocity
  !*************************************************************
      
      if (turbswitch) then
        if (dtftlw.lt..5) then
  !*************************************************************
  !************** CBL options added by mc see routine cblf90 ***
          ! LB needs to be checked if this works with openmp
          if (cblflag.eq.1) then  !modified by mc
            if (-h/ol.gt.5) then  !modified by mc
                flagrein=0
                nrand=nrand+1
                old_wp_buf=part(ipart)%turbvel%w
                call cbl(part(ipart)%turbvel%w,zts,ust,wst,h,rhoa,rhograd,&
                  sigw,dsigwdz,tlw,ptot_lhh,Q_lhh,phi_lhh,ath,bth,ol,flagrein) !inside the routine for inverse time
                part(ipart)%turbvel%w=(part(ipart)%turbvel%w+ath*dtf+&
                  bth*rannumb(nrand)*sqrt(dtf))*real(part(ipart)%icbt) 
                delz=part(ipart)%turbvel%w*dtf
                if (flagrein.eq.1) then
                    call re_initialize_particle(zts,ust,wst,h,sigw,old_wp_buf,nrand,ol)
                    part(ipart)%turbvel%w=old_wp_buf
                    delz=part(ipart)%turbvel%w*dtf
                    nan_count=nan_count+1
                end if             
            else 
                nrand=nrand+1
                old_wp_buf=part(ipart)%turbvel%w
                ath=-part(ipart)%turbvel%w/tlw+sigw*dsigwdz+&
                  part(ipart)%turbvel%w*part(ipart)%turbvel%w/sigw*dsigwdz+sigw*sigw/rhoa*rhograd  !1-note for inverse time should be -wp/tlw*ldirect+... calculated for wp=-wp
                                                                                    !2-but since ldirect =-1 for inverse time and this must be calculated for (-wp) and
                                                                                    !3-the gaussian pdf is symmetric (i.e. pdf(w)=pdf(-w) ldirect can be discarded
                bth=sigw*rannumb(nrand)*sqrt(2.*dtftlw)
                part(ipart)%turbvel%w=(part(ipart)%turbvel%w+ath*dtf+bth)*real(part(ipart)%icbt)  
                delz=part(ipart)%turbvel%w*dtf
                del_test=(1.-part(ipart)%turbvel%w)/part(ipart)%turbvel%w !catch infinity value
                if (isnan(part(ipart)%turbvel%w).or.isnan(del_test)) then 
                    nrand=nrand+1                      
                    part(ipart)%turbvel%w=sigw*rannumb(nrand)
                    delz=part(ipart)%turbvel%w*dtf
                    nan_count2=nan_count2+1
                end if  
            end if
  !******************** END CBL option *******************************            
  !*******************************************************************            
          else
               part(ipart)%turbvel%w=((1.-dtftlw)*part(ipart)%turbvel%w+rannumb(nrand+i)*sqrt(2.*dtftlw) &
               +dtf*(dsigwdz+rhoaux*sigw))*real(part(ipart)%icbt) 
               delz=part(ipart)%turbvel%w*sigw*dtf
          end if
        else
          rw=exp(-dtftlw)
          part(ipart)%turbvel%w=(rw*part(ipart)%turbvel%w+rannumb(nrand+i)*sqrt(1.-rw**2) &
               +tlw*(1.-rw)*(dsigwdz+rhoaux*sigw))*real(part(ipart)%icbt)
          delz=part(ipart)%turbvel%w*sigw*dtf
        endif
        
      else
        rw=exp(-dtftlw)
        part(ipart)%turbvel%w=(rw*part(ipart)%turbvel%w+rannumb(nrand+i)*sqrt(1.-rw**2)*sigw &
             +tlw*(1.-rw)*(dsigw2dz+rhoaux*sigw**2))*real(part(ipart)%icbt)
        delz=part(ipart)%turbvel%w*dtf
      endif

  !****************************************************
  ! Compute turbulent vertical displacement of particle
  !****************************************************

      if (abs(delz).gt.h) delz=mod(delz,h)

  ! Determine if particle transfers to a "forbidden state" below the ground
  ! or above the mixing height
  !************************************************************************

      if (delz.lt.-zts) then         ! reflection at ground
        part(ipart)%icbt=-1
        call set_z(ipart,-zts-delz)
      else if (delz.gt.(h-zts)) then ! reflection at h
        part(ipart)%icbt=-1
        call set_z(ipart,-zts-delz+2.*h)
      else                         ! no reflection
        part(ipart)%icbt=1
        call set_z(ipart,zts+delz)
      endif

      if (i.ne.ifine) then
        zeta=zts/h
        call hanna_short(zts)
      endif
      zts=real(part(ipart)%z)
    end do
    if (cblflag.ne.1) nrand=nrand+i
end subroutine turbulence_boundarylayer

subroutine turbulence_stratosphere(dt,nrand,ux,vy,wp,tropop,zts)
  use par_mod
  use com_mod

  implicit none
  
  integer, intent(inout) ::       &
    nrand                           ! random number used for turbulence
  real, intent(inout) ::          &
    ux,vy,wp                        ! random turbulent velocities above PBL
  real, intent(in) ::             &
    tropop,                       & ! height of troposphere
    zts,                          & ! height of particle
    dt                              ! real(ldt)
  real ::                         &
    uxscale,wpscale,              & ! factor used in calculating turbulent perturbations above PBL
    weight                          ! transition above the tropopause

  if (zts.lt.tropop) then  ! in the troposphere
    uxscale=sqrt(2.*d_trop/dt)
    if (nrand+1.gt.maxrand) nrand=1
    ux=rannumb(nrand)*uxscale
    vy=rannumb(nrand+1)*uxscale
    nrand=nrand+2
    wp=0.
  else if (zts.lt.tropop+1000.) then     ! just above the tropopause: make transition
    weight=(zts-tropop)/1000.
    uxscale=sqrt(2.*d_trop/dt*(1.-weight))
    if (nrand+2.gt.maxrand) nrand=1
    ux=rannumb(nrand)*uxscale
    vy=rannumb(nrand+1)*uxscale
    wpscale=sqrt(2.*d_strat/dt*weight)
    wp=rannumb(nrand+2)*wpscale+d_strat/1000.
    nrand=nrand+3
  else                 ! in the stratosphere
    if (nrand.gt.maxrand) nrand=1
    ux=0.
    vy=0.
    wpscale=sqrt(2.*d_strat/dt)
    wp=rannumb(nrand)*wpscale
    nrand=nrand+1
  endif
end subroutine turbulence_stratosphere

subroutine turbulence_mesoscale(nrand,dxsave,dysave,ipart,eps_eta)
  use point_mod
  use par_mod
  use com_mod
  use interpol_mod
  use hanna_mod
  use cmapf_mod
  use coordinates_ecmwf
  use particle_mod

  implicit none
  integer, intent(inout) ::       &
    nrand                           ! random number used for turbulence
  integer, intent(in) ::          &
    ipart                              ! particle index
  real, intent(in) ::             &
    eps_eta
  real, intent(inout) ::          &
    dxsave,dysave                   ! accumulated displacement in long and lat
  real ::                         &
    r,rs,                         & ! mesoscale related
    ux,vy                           ! random turbulent velocities above PBL

  r=exp(-2.*real(abs(lsynctime))/real(lwindinterv))
  rs=sqrt(1.-r**2)
  if (nrand+2.gt.maxrand) nrand=1
  part(ipart)%mesovel%u=r*part(ipart)%mesovel%u+rs*rannumb(nrand)*usig*turbmesoscale
  part(ipart)%mesovel%v=r*part(ipart)%mesovel%v+rs*rannumb(nrand+1)*vsig*turbmesoscale
  dxsave=dxsave+part(ipart)%mesovel%u*real(lsynctime)
  dysave=dysave+part(ipart)%mesovel%v*real(lsynctime)

  select case (wind_coord_type)
    case ('ETA')
      part(ipart)%mesovel%w=r*part(ipart)%mesovel%w+rs*rannumb(nrand+2)*wsigeta*turbmesoscale
      call update_zeta(ipart,part(ipart)%mesovel%w*real(lsynctime))
      if (part(ipart)%zeta.ge.1.) call set_zeta(ipart,1.-(part(ipart)%zeta-1.))
      if (part(ipart)%zeta.eq.1.) call update_zeta(ipart,-eps_eta)

    case ('METER')
      part(ipart)%mesovel%w=r*part(ipart)%mesovel%w+rs*rannumb(nrand+2)*wsig*turbmesoscale
      call update_z(ipart,part(ipart)%mesovel%w*real(lsynctime))
      if (part(ipart)%z.lt.0.) call set_z(ipart,-1.*part(ipart)%z)    ! if particle below ground -> refletion

    case default
      part(ipart)%mesovel%w=r*part(ipart)%mesovel%w+rs*rannumb(nrand+2)*wsig*turbmesoscale
      call update_z(ipart,part(ipart)%mesovel%w*real(lsynctime))
      if (part(ipart)%z.lt.0.) call set_z(ipart,-1.*part(ipart)%z)    ! if particle below ground -> refletion
  end select
end subroutine turbulence_mesoscale


end module turbulence_mod