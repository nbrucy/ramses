!!! FlorentR - PATCH GROUP
!========================================================================================
!== Initial conditions to setup 1 or more disk galaxy - Florent Renaud - March 17, 2011
!== Adapted from Frederic Bournaud's and Damien Chapon's patch "merger"
!========================================================================================
!==  Namelist settings: each variable is a list (one item per galaxy)
!==
!== 	gal_center_x :		x-coordinates of the centers of mass of the galaxies [kpc]
!== 	gal_center_y :		y-coordinates of the centers of mass of the galaxies [kpc]
!== 	gal_center_z :		z-coordinates of the centers of mass of the galaxies [kpc]
!== 	Vgal_x :		vx-coordinates of the centers of mass of the galaxies [km/s]
!== 	Vgal_y :		vy-coordinates of the centers of mass of the galaxies [km/s]
!== 	Vgal_z :		vz-coordinates of the centers of mass of the galaxies [km/s]
!== 	gal_axis_x : 		x-components of the angular momentum of the disks
!== 	gal_axis_y : 		x-components of the angular momentum of the disks
!== 	gal_axis_z : 		x-components of the angular momentum of the disks
!==     Mgas_disk :		gaseous masses of the disks [10^9 Msun]
!== 	typ_radius :		caracteristic radii of the disks [kpc]
!== 	cut_radius :		truncation radii of the disks [kpc]
!== 	typ_height :		caracteristic heights of the disks [kpc]
!== 	cut_height :		truncation heights of the disks [kpc]
!== 	rad_profile :		radial density profile of the disks ['exponential' or 'Toomre']
!== 	z_profile :		vertical density profile of the disks ['exponential' or 'gaussian']
!==	Vcirc_dat_file :	path to the files containing the velocity curves (col1 = radii in pc, col2 = Vcirc in km/s) [ORDERED IN RADIUS, NO DUPLICATED LINES]
!==	ic_part_file :		file in the IC directory that contains the particle data (see also init_part)
!==
!==  Global scalar:
!== 	IG_density_factor :	density contrast for the intergalactic medium
!==
!========================================================================================


module group_commons
  use amr_commons

  integer,parameter::MAXNGAL=10 ! maximum number of galaxies
  integer::ngal ! actual number of galaxies

  ! Galactic data
  real(dp), dimension(1:MAXNGAL)::Mgas_disk
  real(dp), dimension(1:MAXNGAL)::typ_radius
  real(dp), dimension(1:MAXNGAL)::cut_radius
  real(dp), dimension(1:MAXNGAL)::typ_height
  real(dp), dimension(1:MAXNGAL)::cut_height
  character(len=16), dimension(1:MAXNGAL)::rad_profile
  character(len=16), dimension(1:MAXNGAL)::z_profile
  real(dp), dimension(1:3,1:MAXNGAL)::gal_pos ! x, y, z
  real(dp), dimension(1:3,1:MAXNGAL)::gal_vel ! vx, vy, vz
  real(dp), dimension(1:3,1:MAXNGAL)::gal_axis ! axis_x, axis_y, axis_z
  real(dp), dimension(1:MAXNGAL)::gal_rho ! central_density
  real(dp)::dmin ! intergalactic density

  ! velocity data
  integer, dimension(1:MAXNGAL)::vcirc_nsample
  real(dp), dimension(:,:,:), allocatable::vcirc_dat

  ! particle data
  character(len=512), dimension(1:MAXNGAL)::ic_part_file

!!! TEMP
  ! if compatibility_vfactor == .false. (i.e NOT defined in the GROUP_PARAMS in the namelist), velocities in 'ic_part' are in km/s (This is the new version). See also init_part.f90
  ! if compatibility_vfactor == .true., velocities in 'ic_part' are in code units and no scaling will be done (This is the old version). See also init_part.f90
  ! THIS IS ONLY FOR THE PARTICLES VELOCITIES, NOT THE GAS, NOR THE GALAXIES
  logical::compatibility_vfactor=.false.
!!!

end module group_commons

!==================================================================================
!==================================================================================
!==================================================================================

subroutine condinit_group(x,u,dx,nn)

  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================

  use group_commons
  use amr_parameters
  use hydro_parameters
  implicit none

  ! amr data
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Position of cell center

  ! Namelist definitions
  namelist/group_params/ gal_center_x, gal_center_y, gal_center_z, &
       & Vgal_x, Vgal_y, Vgal_z, gal_axis_x, gal_axis_y, gal_axis_z, &
       & Mgas_disk, typ_radius, cut_radius, typ_height, cut_height, &
       & rad_profile, z_profile, Vcirc_dat_file, IG_density_factor, &
       & ic_part_file, &
!!! TEMP
       & compatibility_vfactor
!!!

  ! galaxy data
  real(dp), dimension(1:MAXNGAL)::gal_center_x
  real(dp), dimension(1:MAXNGAL)::gal_center_y
  real(dp), dimension(1:MAXNGAL)::gal_center_z
  real(dp), dimension(1:MAXNGAL)::Vgal_x
  real(dp), dimension(1:MAXNGAL)::Vgal_y
  real(dp), dimension(1:MAXNGAL)::Vgal_z
  real(dp), dimension(1:MAXNGAL)::gal_axis_x
  real(dp), dimension(1:MAXNGAL)::gal_axis_y
  real(dp), dimension(1:MAXNGAL)::gal_axis_z

  ! misc
  character(len=512), dimension(1:MAXNGAL)::Vcirc_dat_file
  real(dp)::IG_density_factor = 1.0D-5
  integer::membership ! index of the galaxy containing the current cell
  integer:: i, k, ind_gal, ierr, d ! counters, error handlers
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2, a2
  real(dp)::mfactor, kpcfactor, pcfactor, vfactor, pi
  real(dp)::pot, pot_tmp, norm, r, rr, abs_z, vcirc, weight, vrot, dmin_tmp
  real(dp),dimension(3,MAXNGAL)::xx
  real(dp),dimension(3)::xx_rad

  ! utils
  logical,save:: init_nml=.false.
  logical::nml_ok=.true.
  character(LEN=80)::infile
  logical::file_exists


  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  a2 = T2_star / scale_T2 ! sound speed squared

  ! Init: read user-defined group parameters from the namelist, scale parameters, compute central and intergalactic densities
  if (.not. init_nml) then

    pi = dacos(-1.0D0)

    call getarg(1,infile) ! get the name of the namelist
    open(1,file=infile)
    rewind(1) ! really needed?
    read(1,NML=group_params,END=106) ! open namelist and read the group_params data
    goto 107
    106 if(myid==1)write(*,*)' You have to set up namelist &GROUP_PARAMS in parameter file'
    call clean_stop
    107 continue
    close(1)

!!! TEMP
    if(compatibility_vfactor .and. myid==1) write(*,*)' Velocities of particles will not be scaled.'
!!!

    ngal = size(gal_center_x) ! get the number of galaxies from the namelist
    if(ngal>MAXNGAL) then
      if(myid==1)write(*,*) 'Error: number of galaxies > MAXNGAL (=',MAXNGAL, ')'
      call clean_stop
    endif

    ! store galaxy data
    do ind_gal=ngal, 1, -1 ! Loop over galaxies (REVERSE ORDER)
      gal_pos(:,ind_gal) = (/ gal_center_x(ind_gal), gal_center_y(ind_gal), gal_center_z(ind_gal) /)
      gal_vel(:,ind_gal) = (/ Vgal_x(ind_gal), Vgal_y(ind_gal), Vgal_z(ind_gal) /)
      gal_axis(:,ind_gal) = (/ gal_axis_x(ind_gal), gal_axis_y(ind_gal), gal_axis_z(ind_gal) /)

      ! remove zero-mass galaxies
      if(Mgas_disk(ind_gal)==0.0D0) then
        ngal = ngal - 1 ! decrease counter
        if(ind_gal==ngal) goto 108 ! last galaxy in the array: do nothing
        ! left-shift the rest of the array (i.e. non-zero-mass galaxies)
        gal_pos(:,ind_gal:MAXNGAL-1) = gal_pos(:,ind_gal+1:MAXNGAL)
        gal_vel(:,ind_gal:MAXNGAL-1) = gal_vel(:,ind_gal+1:MAXNGAL)
        gal_axis(:,ind_gal:MAXNGAL-1) = gal_axis(:,ind_gal+1:MAXNGAL)
        Mgas_disk(ind_gal:MAXNGAL-1) = Mgas_disk(ind_gal+1:MAXNGAL)
        typ_radius(ind_gal:MAXNGAL-1) = typ_radius(ind_gal+1:MAXNGAL)
        cut_radius(ind_gal:MAXNGAL-1) = cut_radius(ind_gal+1:MAXNGAL)
        typ_height(ind_gal:MAXNGAL-1) = typ_height(ind_gal+1:MAXNGAL)
        cut_height(ind_gal:MAXNGAL-1) = cut_height(ind_gal+1:MAXNGAL)
        rad_profile(ind_gal:MAXNGAL-1) = rad_profile(ind_gal+1:MAXNGAL)
        z_profile(ind_gal:MAXNGAL-1) = z_profile(ind_gal+1:MAXNGAL)
        Vcirc_dat_file(ind_gal:MAXNGAL-1) = Vcirc_dat_file(ind_gal+1:MAXNGAL)
        ic_part_file(ind_gal:MAXNGAL-1) = ic_part_file(ind_gal+1:MAXNGAL)
      endif
      108 continue
    enddo ! End loop over galaxies

    do ind_gal=1, ngal ! Loop over galaxies
      ! Normalization of axis
      norm = sqrt(dot_product(gal_axis(:,ind_gal),gal_axis(:,ind_gal)))
      if(norm==0.0D0) then
        if(myid==1)write(*,*)'Error: Galactic axis(#', ind_gal,') is zero '
        nml_ok=.false.
      endif
      if(norm.NE.1.0D0) then
        gal_axis(:,ind_gal) = gal_axis(:,ind_gal) / norm
      end if

      ! Check for circular velocity files
!      inquire(file=trim(Vcirc_dat_file(ind_gal)), exist=file_exists)
      inquire(file=trim(initfile(levelmin))//'/'//trim(Vcirc_dat_file(ind_gal)), exist=file_exists)
      if(.NOT. file_exists) then
        if(myid==1)write(*,*)"Error: Vcirc_dat_file ", trim(Vcirc_dat_file(ind_gal)), " doesn't exist "
        nml_ok=.false.
      end if

      ! Check for particle files
      inquire(file=trim(initfile(levelmin))//'/'//trim(ic_part_file(ind_gal)), exist=file_exists)
      if(.NOT. file_exists) then
        if(myid==1)write(*,*)"Error: ic_part_file ", trim(ic_part_file(ind_gal)), " doesn't exist "
        nml_ok=.false.
      end if

      ! Count the number of data points in the vcirc file
      vcirc_nsample(ind_gal) = 0 ! init
      open(unit=123, file=TRIM(initfile(levelmin))//'/'//TRIM(Vcirc_dat_file(ind_gal)), iostat=ierr)
!      open(unit=123, file=trim(Vcirc_dat_file(ind_gal)), iostat=ierr)
      do while(ierr==0)
        read(123,*,iostat=ierr)
        if(ierr==0) then
          vcirc_nsample(ind_gal) = vcirc_nsample(ind_gal) + 1  ! Number of data-points
        end if
      end do
      close(123)
    enddo ! End loop over galaxies

    if(.not. nml_ok)then
      if(myid==1)write(*,*)'Too many errors in the namelist'
      if(myid==1)write(*,*)'Aborting...'
      call clean_stop
    end if

    ! read the velocity curves data
    allocate(vcirc_dat(1:maxval(vcirc_nsample),2,ngal)) ! warning: never deallocated
    Vcirc_dat = 0.0D0 ! init

    do ind_gal=1, ngal ! Loop over galaxies
      open(unit=123, file=TRIM(initfile(levelmin))//'/'//TRIM(Vcirc_dat_file(ind_gal)), iostat=ierr)
!      open(unit=123, file=trim(Vcirc_dat_file(ind_gal)), iostat=ierr)
      do i=1, vcirc_nsample(ind_gal)
        read(123,*) vcirc_dat(i,:,ind_gal) ! read the velocity data
      end do
      close(123)
    enddo

    ! data has been read and understood
    init_nml = .true.

    ! Scaling (done only once!)
    mfactor= 1.0D9 * 1.9891D33  / (scale_d * scale_l**3)
    kpcfactor = 3.085677581282D21 / scale_l
    pcfactor = kpcfactor/1.0D3
    vfactor = 1.0D5 / scale_v

    vcirc_dat(:,1,:) = vcirc_dat(:,1,:) * pcfactor
    vcirc_dat(:,2,:) = vcirc_dat(:,2,:) * vfactor

    do ind_gal=1, ngal ! Loop over galaxies
      Mgas_disk(ind_gal) = Mgas_disk(ind_gal) * mfactor
      typ_radius(ind_gal) = typ_radius(ind_gal) * kpcfactor
      cut_radius(ind_gal) = cut_radius(ind_gal) * kpcfactor
      typ_height(ind_gal) = typ_height(ind_gal) * kpcfactor
      cut_height(ind_gal) = cut_height(ind_gal) * kpcfactor
      gal_pos(:,ind_gal) = gal_pos(:,ind_gal) * kpcfactor
      gal_vel(:,ind_gal) = gal_vel(:,ind_gal) * vfactor

      ! Central gas density and intergalactic gas density
      select case (rad_profile(ind_gal))
        case ('Toomre')
          !write(*,*) 'Galaxy #', ind_gal, 'rad_profile TOOMRE'
          gal_rho(ind_gal) = dsqrt(1.0D0 + cut_radius(ind_gal)**2/typ_radius(ind_gal)**2) - 1.0D0
          dmin_tmp = 1.0D0 / dsqrt(1.0D0 + cut_radius(ind_gal)**2/typ_radius(ind_gal)**2)
        case default ! exponential disk
          !write(*,*) 'Galaxy #', ind_gal, 'rad_profile EXP'
          gal_rho(ind_gal) = 1.0D0 - exp(-cut_radius(ind_gal) / typ_radius(ind_gal)) * (1.0D0 + cut_radius(ind_gal) / typ_radius(ind_gal))
          dmin_tmp = exp(-cut_radius(ind_gal) / typ_radius(ind_gal))
      end select

      !! JFENSCH
      !select case (z_profile(ind_gal))
      !  case ('exponential')
      !write(*,*) 'Galaxy #', ind_gal, 'z_profile EXP'
      gal_rho(ind_gal) = gal_rho(ind_gal) * (1.0D0 - exp(-cut_height(ind_gal) / typ_height(ind_gal)))
      dmin_tmp = dmin_tmp * exp(-cut_height(ind_gal) / typ_height(ind_gal))
      !  case default ! gaussian
      !    write(*,*) 'Galaxy #', ind_gal, 'z_profile Gaussian'
      !    gal_rho(ind_gal) = gal_rho(ind_gal) * (dsqrt(pi/2.0D0) * erf(cut_height(ind_gal) / (dsqrt(2.0D0)*typ_height(ind_gal))))
      !    dmin_tmp = dmin_tmp * exp(-0.5D0 * cut_height(ind_gal)**2 / typ_height(ind_gal)**2 )
      !end select
      !! JFENSCH

      gal_rho(ind_gal) = Mgas_disk(ind_gal) / (4.0D0 * pi * typ_radius(ind_gal)**2 * typ_height(ind_gal) * gal_rho(ind_gal))

      !write(*,*)'Galaxy #', ind_gal, 'gal_rho =', gal_rho
      !write(*,*)'Galaxy #', ind_gal, 'dmin temp =', dmin_tmp

      if(ind_gal==1) then
        dmin = gal_rho(ind_gal) * dmin_tmp
        !write(*,*)'dmin 1st galaxy = ', dmin
      else
        dmin = min(dmin, gal_rho(ind_gal) * dmin_tmp)
        !write(*,*)'dmin 2nd galaxy = ', dmin
      endif
    enddo ! End loop over galaxies

    !write(*,*)'dmin after loop =', dmin
    !write(*,*)'IG_density_factor =', IG_density_factor

    dmin = IG_density_factor * dmin

    !write(*,*)'Final dmin =', dmin

  end if
  ! end of initialization


  ! Loop over cells
  do i=1,nn
    do ind_gal=1, ngal ! Loop over galaxies
      do d=1,3
        xx(d,ind_gal) = x(i,d) - (gal_pos(d,ind_gal) + boxlen / 2.0D0) ! distance to galaxy center
      enddo
      ! determine the galactic membership of the current cell
      pot_tmp = dot_product(xx(:,ind_gal),xx(:,ind_gal)) / Mgas_disk(ind_gal)
      if(ind_gal==1) then
        pot = pot_tmp
        membership = ind_gal
      else
        if(pot_tmp<pot) then ! deeper in the potential well
          pot = pot_tmp ! set new minimum
          membership = ind_gal ! set new galaxy membership
        endif
      endif
    enddo ! End loop over galaxies


    xx_rad = xx(:,membership) - dot_product(xx(:,membership),gal_axis(:,membership)) * gal_axis(:,membership)
    r = dsqrt(dot_product(xx_rad,xx_rad)) ! cylindric radius: distance between the cell and the galactic rotation axis
    abs_z = dsqrt(dot_product(xx(:,membership),xx(:,membership)) - r**2) ! distance to the galactic plane


    if( ((r-dx/2.0D0)<cut_radius(membership)) .AND. ((abs_z-dx/2.0D0)<cut_height(membership)) )then ! cell in the disk : analytical density profile + rotation velocity
      weight = ( min(r+dx/2.0D0, cut_radius(membership)) - (r-dx/2.0D0) ) / dx
      if (weight.NE.1.0D0) r = r + (weight - 1.0D0) * dx/2.0D0 ! cell partially outside the disk
      weight = weight * ( min(abs_z + dx / 2.0D0, cut_height(membership)) - (abs_z - dx / 2.0D0) ) / dx

      ! Circular velocity
      k=2
      do while (r>vcirc_dat(k,1,membership) .AND. k.NE.vcirc_nsample(membership))
        k = k + 1
      end do
      ! if radius > largest radius defined: keep the last slope available (should be ~flat).
      ! 1st order interpolation
      vcirc = vcirc_dat(k-1,2,membership) + (r-vcirc_dat(k-1,1,membership)) / (vcirc_dat(k,1,membership)-vcirc_dat(k-1,1,membership)) * (vcirc_dat(k,2,membership)-vcirc_dat(k-1,2,membership))

      ! Density
      select case (rad_profile(membership))
        case ('Toomre')
          q(i,1) = gal_rho(membership) / dsqrt(1.0D0 + r**2/typ_radius(membership)**2)
        case default ! exponential disk
          q(i,1) = gal_rho(membership) * exp(-r / typ_radius(membership))
      end select

      select case (z_profile(membership))
        case ('exponential')
          q(i,1) = q(i,1) * exp(-abs_z / typ_height(membership))
        case default ! gaussian
          q(i,1) = q(i,1) * exp(-0.5D0 * abs_z**2 / typ_height(membership)**2)
      end select

      q(i,1) = max(weight * q(i,1), dmin)

      ! V = vrot * (u_rot^xx_rad)/r + Vx_gal    -> vrot = sqrt(Vcirc**2 - 3*Cs**2 + r/rho * grad(rho) * Cs**2)
      select case (rad_profile(membership))
        case ('Toomre')
          vrot = dsqrt(max(vcirc**2 - 3.0D0*a2 - r**2/(r**2+typ_radius(membership)**2)*a2, 0.0D0))
        case default ! exponential disk
          vrot = dsqrt(max(vcirc**2 - 3.0D0*a2 - r/typ_radius(membership) * a2, 0.0D0))
      end select

      q(i,2) = weight * vrot * ( gal_axis(2,membership) * xx_rad(3) - gal_axis(3,membership) * xx_rad(2) ) /r + gal_vel(1,membership)
      q(i,3) = weight * vrot * ( gal_axis(3,membership) * xx_rad(1) - gal_axis(1,membership) * xx_rad(3) ) /r + gal_vel(2,membership)
      q(i,4) = weight * vrot * ( gal_axis(1,membership) * xx_rad(2) - gal_axis(2,membership) * xx_rad(1) ) /r + gal_vel(3,membership)

    else ! Cell out of the gaseous disk : density = IG density, velocity = v_gal
      q(i,1) = dmin
      q(i,2:4) = gal_vel(:,membership)
    endif

    q(i,5)=a2*q(i,1) ! P = rho * a**2 = rho * Cs**2
  enddo ! End loop over cells


  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum : Omega = rho * V
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0D0
  u(1:nn,5)=u(1:nn,5)+0.5D0*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5D0*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5D0*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy (E = Ec + P / (gamma - 1))
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! passive scalars
  do i=6,nvar
     u(1:nn,i)=q(1:nn,1)*q(1:nn,i)
  end do

end subroutine condinit_group

!==================================================================================
!==================================================================================
!==================================================================================
!!! FRenaud


!==================================================================================
!==================================================================================
!==================================================================================
