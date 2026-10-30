
module deltaE_module

  use amr_commons

  implicit none

  ! DeltaE params
  logical::deltaE_enable=.true. ! Whether to enable deltaE output
  logical::deltaE_correct_pressure_fix=.true. ! Whether to be "pressure_fix" aware when computing energy
  logical::deltaE_debug=.false. ! print intermediate energies
  logical::deltaE_force_all_levels=.false. ! Force computation on all levels for all processes

  ! Zone settings (similar to geometrical refine zone definition)
  ! The regions are mutually exclusive: a cell/particles can belong on to one region, with the lowest index
  integer,parameter :: DELTAE_MAX_REGIONS=99
  integer::deltaE_nregions=0 ! Number of deltaE zones
  real(dp),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_size =-1       ! Size of refinement region
  real(dp),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_x = 0          ! Center of deltaE region
  real(dp),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_y = 0          ! Center of deltaE region
  real(dp),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_z = 0          ! Center of deltaE region
  real(dp),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_exp = 2        ! Exponent for distance
  real(dp),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_a = 1          ! Ellipticity (Y/X)
  real(dp),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_b = 1          ! Ellipticity (Z/X)
  character(len=20),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_name = "" ! Name of the region

  ! DeltaE status flags
  logical, dimension(:), allocatable :: deltaE_use_unew
  logical, dimension(:), allocatable :: deltaE_use_phi_old

  ! Arrays
  ! Energies. To add an energy, add its index and name at the same position
  enum, bind(c)
  enumerator :: iekin = 1, iepot, ieint, iemag, iekin_gas, &
                        & iekin_part, iepot_gas, iepot_part, &
                        & iekin_part_star, iekin_part_dm, &
                        & iepot_part_star, iepot_part_dm, &
                        nb_energy_end
  end enum
  integer, parameter :: nb_energy_fixed = nb_energy_end - 1
  character(len=14), dimension(1:nb_energy_fixed) :: energy_names = [character(len=14) :: &
                        & "ekin", "epot", "eint", &
                        & "emag", "ekin_gas", "ekin_part",&
                        & "epot_gas", "epot_part", &
                        & "ekin_part_star", "ekin_part_dm", &
                        & "epot_part_star", "epot_part_dm"]

  ! Turbulent kinetic energy, one reservoir per level: index iekin_turb0 + ilevel
  integer, parameter :: iekin_turb0 = nb_energy_fixed
  integer, parameter :: nb_energy_kind = nb_energy_fixed + MAXLEVEL

  ! Per grid totals for the turbulent energy
  real(dp), allocatable, dimension(:,:) :: grid_u  ! (1:ngridmax, 1:ndim+1)

  type process
  real(dp), dimension(1:nb_energy_kind, 0:DELTAE_MAX_REGIONS) :: v
  character(len=20) :: name
  logical :: compute_eturb
  logical :: compute_gas
  logical :: compute_part
end type

! Process indices. To add a process: add its index here and one line in initialize_processes.
enum, bind(c)
  enumerator :: iflux_gas = 1, iflux_part, icooling, igravity_gas, igravity_part, &
              & istar_formation, ifeedback, iturb_driving, ipoisson, &
              & imagnetic_diffusion, irefine, icorrections, nb_process_end
end enum
integer, parameter :: nb_process = nb_process_end - 1

type(process), dimension(nb_process) :: deltaE

contains

  subroutine read_deltaE_params(namelist_unit,nml_ok)

    integer,intent(in)::namelist_unit
    logical,intent(inout)::nml_ok
    integer::nml_err

    namelist/deltaE_params/deltaE_correct_pressure_fix,deltaE_debug,deltaE_enable,deltaE_force_all_levels, &
    & deltaE_nregions,deltaE_region_size,deltaE_region_x,deltaE_region_y,deltaE_region_z, &
    & deltaE_region_exp,deltaE_region_a,deltaE_region_b,deltaE_region_name


   ! Go to the beginning of the file
    rewind(namelist_unit)

    ! Read namelist
    read(namelist_unit,NML=deltaE_params,IOSTAT=nml_err)

    if(nml_err>0)then
      if(myid==1)write(*,*)'Error reading namelist &DELTAE_PARAMS. Check formatting.'
      nml_ok=.false.
    end if

    if(deltaE_nregions < 0 .or. deltaE_nregions > DELTAE_MAX_REGIONS) then
      if(myid==1)write(*,*)'Error in &DELTAE_PARAMS. deltaE_nregions should be positive and below DELTAE_MAX_REGIONS'
      nml_ok=.false.
    end if

    allocate(deltaE_use_unew(levelmin:nlevelmax))
    allocate(deltaE_use_phi_old(levelmin:nlevelmax))

    deltaE_use_phi_old = .false.
    deltaE_use_unew = .false.

    call initialize_processes

  end subroutine read_deltaE_params

  subroutine initialize_processes
    !      index                          value   name                  eturb    gas      part
    deltaE(iflux_gas)           = process(0.0d0, "flux_gas",           .true.,  .true.,  .false.)
    deltaE(iflux_part)          = process(0.0d0, "flux_part",          .false., .false., .true.)
    deltaE(icooling)            = process(0.0d0, "cooling",            .false., .true.,  .false.)
    deltaE(igravity_gas)        = process(0.0d0, "gravity_gas",        .true.,  .true.,  .false.)
    deltaE(igravity_part)       = process(0.0d0, "gravity_part",       .false., .false., .true.)
    deltaE(istar_formation)     = process(0.0d0, "star_formation",     .true.,  .true.,  .true.)
    deltaE(ifeedback)           = process(0.0d0, "feedback",           .true.,  .true.,  .true.)
    deltaE(iturb_driving)       = process(0.0d0, "turb_driving",       .true.,  .true.,  .false.)
    deltaE(ipoisson)            = process(0.0d0, "poisson",            .false., .true.,  .true.)
    deltaE(imagnetic_diffusion) = process(0.0d0, "magnetic_diffusion", .false., .true.,  .false.)
    deltaE(irefine)             = process(0.0d0, "refine",             .true.,  .true.,  .true.)
    deltaE(icorrections)        = process(0.0d0, "corrections",        .true.,  .true.,  .true.)
  end subroutine

  subroutine reset_processes
    integer :: j
    do j = 1, nb_process
      deltaE(j)%v = 0.0d0
    end do
  end subroutine

  function energy_name(i) result(name)
    integer, intent(in) :: i
    character(len=20) :: name
    if (i <= nb_energy_fixed) then
      name = energy_names(i)
    else
      write(name, '(A,I2.2)') "ekin_gas_turb_", i - iekin_turb0
    end if
  end function

  function region_name(k) result(name)
    integer, intent(in) :: k
    character(len=20) :: name
    if (k == 0) then
      name = ""
    else
      if (deltaE_region_name(k) == "") then
        write(name, '(A,I2.2)') "region_", k
      else
        name = deltaE_region_name(k)
      end if
    end if
  end function

  logical function is_energy_printed(i)
    integer, intent(in) :: i
    is_energy_printed = (i <= nb_energy_fixed) .or. &
         & (i - iekin_turb0 >= levelmin .and. i - iekin_turb0 <= nlevelmax)
  end function

  subroutine print_energies(energies)
    implicit none
    real(dp), dimension(1:nb_energy_kind, 0:deltaE_nregions), intent(in) :: energies
    integer :: i, k

    write (*, *) "Energy output"
    do k = 0, deltaE_nregions
      do i = 1, nb_energy_kind
        if (.not. is_energy_printed(i)) cycle
        write (*, 996) region_name(k), energy_name(i), energies(i, k)
      end do
    end do
    write (*, *) "End energy output"

996 format(A20, A20, 1pe20.7)

  end subroutine

  subroutine print_processes
    implicit none
    integer :: i, j, k

    write (*, *) "DeltaE output "
    write (*, 997) "", "", (deltaE(j)%name, j = 1, nb_process)
    do k = 0, deltaE_nregions
      do i = 1, nb_energy_kind
        if (.not. is_energy_printed(i)) cycle
        write (*, 998) region_name(k), energy_name(i), (deltaE(j)%v(i,k), j = 1, nb_process)
      end do
    end do

    write (*, *) "End deltaE output"

997 format(A20, A20, *(A20))
998 format(A20, A20, *(1pe20.5))
  end subroutine

  subroutine compute_energies(ilevel, energies, compute_eturb, compute_gas, compute_part)
    implicit none

    integer, intent(in):: ilevel
    real(dp), dimension(1:nb_energy_kind, 0:deltaE_nregions), intent(inout) :: energies
    logical, intent(in), optional :: compute_eturb, compute_gas, compute_part
    logical compute_eturb_loc, compute_gas_loc, compute_part_loc
    real(dp), dimension(1:nlevelmax, 0:deltaE_nregions) :: eturb

    compute_eturb_loc = .true. ; compute_gas_loc = .true. ; compute_part_loc = .true.
    if (present(compute_eturb)) compute_eturb_loc = compute_eturb
    if (present(compute_gas)) compute_gas_loc = compute_gas
    if (present(compute_part)) compute_part_loc = compute_part

    energies = 0.0d0
    if (hydro .and. compute_gas_loc) then
      call compute_energy_gas(ilevel, energies)
    end if
    if (pic .and. compute_part_loc) then
      call compute_energy_part(ilevel, energies)
    end if
    if (hydro .and. compute_eturb_loc) then
      eturb = 0.0d0
      call compute_ekin_turb_all_levels(eturb)
      energies(iekin_turb0 + levelmin:iekin_turb0 + nlevelmax, :) = eturb(levelmin:nlevelmax,:)
    end if
    energies(iekin,:) = energies(iekin_gas,:) + energies(iekin_part,:)
    energies(iepot,:) = energies(iepot_gas,:) + energies(iepot_part,:)
  end subroutine

  subroutine compute_transfer(levelstart, levelend, deltaE_process, step)

    implicit none
    integer, intent(in):: levelstart, levelend
    type(process), intent(inout) :: deltaE_process
    integer, intent(in) :: step
    integer :: act_levelstart, act_levelend

    integer :: ilevel

    real(dp), dimension(1:nb_energy_kind,0:DELTAE_MAX_REGIONS), save :: energy_before
    real(dp), dimension(1:nb_energy_kind,0:deltaE_nregions) ::energy_level, energy_after

    if (deltaE_force_all_levels) then
      act_levelstart = levelmin
      act_levelend = nlevelmax
    else
      act_levelstart = levelstart
      act_levelend = levelend
    end if

    if (step == 1) then
      energy_before = 0.0d0
    else
      energy_after = 0.0d0
    end if

    do ilevel = act_levelstart, act_levelend
      call compute_energies(ilevel, energy_level, &
      & compute_eturb = (ilevel == act_levelstart .and. deltaE_process%compute_eturb), &
      & compute_gas = deltaE_process%compute_gas, &
      & compute_part = deltaE_process%compute_part)
      if (step == 1) then
        energy_before(:, 0:deltaE_nregions) = energy_before(:,0:deltaE_nregions) + energy_level
      else
        energy_after = energy_after + energy_level
      end if
    end do


    if (step == 1) then
      call reduce_energies(energy_before(:, 0:deltaE_nregions))
    else
      call reduce_energies(energy_after)
    end if

    if (step == 2) then

      deltaE_process%v(:,0:deltaE_nregions) = deltaE_process%v(:,0:deltaE_nregions) &
         &  + energy_after(:, 0:deltaE_nregions) - energy_before(:, 0:deltaE_nregions)

      if (deltaE_debug .and. myid == 1) then
        write(*,*) "DeltaE debug before ", deltaE_process%name
        call print_energies(energy_before)
        write(*,*) "DeltaE debug after ", deltaE_process%name
        call print_energies(energy_after)
      end if
    end if

  end subroutine

  subroutine reduce_energies(energies)
    ! Sum over ranks of the local energies
    use mpi_mod
    implicit none
    real(dp), dimension(1:nb_energy_kind, 0:deltaE_nregions), intent(inout) :: energies
#ifndef WITHOUTMPI
    real(dp), dimension(1:nb_energy_kind, 0:deltaE_nregions) :: energies_loc
    integer :: info
    energies_loc = energies
    call MPI_ALLREDUCE(energies_loc, energies, nb_energy_kind*(deltaE_nregions + 1), MPI_DOUBLE_PRECISION, &
         & MPI_SUM, MPI_COMM_WORLD, info)
#endif
  end subroutine reduce_energies

  subroutine compute_energy_gas(ilevel, energies)
    use amr_commons
    use hydro_commons
    use poisson_commons
    use mpi_mod
    implicit none
#ifndef WITHOUTMPI
    integer::info
    real(kind=8), dimension(4,0:deltaE_nregions)::comm_buffin, comm_buffout
#endif
    integer, intent(in)::ilevel
    real(dp), dimension(1:nb_energy_kind, 0:deltaE_nregions), intent(inout) :: energies

    integer::i, ivar, ind, ncache, igrid, ix, iy, iz, iskip, idim, k
    integer::nleaf, ngrid, nx_loc
    integer, dimension(1:nvector), save::ind_grid, ind_cell, ind_leaf, ind_grid_of_leaf

    integer, dimension(1:nvector) :: region_index = 0
    real(dp)::er,xr,yr,zr,rr,xn,yn,zn,r,aa,bb
    real(dp),dimension(1:3)::skip_loc
    real(dp),dimension(1:twotondim,1:3)::xc
    real(dp),dimension(1:nvector,1:ndim),save::xx
    real(dp):: dx_code, dx, vol, scale


    real(kind=8)::mass_loc, ekin_leaf, emag_leaf, epot_leaf
    real(kind=8), dimension(0:deltaE_nregions+1) :: ekin_loc, eint_loc, emag_loc, epot_loc
    real(dp), dimension(1:nvector, 1:nvar_all), save::uu
    real(dp), dimension(1:nvector), save::phi_level
    real(dp), dimension(1:nvector), save::rhoinv
    integer, dimension(1:nvar_all)::ivar_used
    integer::kv, nused



    real(dp):: e_cons, e_prim, e_trunc, div

    mass_loc = 0.0d0
    ekin_loc = 0.0d0; ekin_leaf = 0.0d0
    emag_loc = 0.0d0; emag_leaf = 0.0d0
    eint_loc = 0.0d0
    epot_loc = 0.0d0; epot_leaf = 0.0d0

    energies(iekin_gas,:) = 0.0d0
    energies(ieint,:) = 0.0d0
    energies(iemag,:) = 0.0d0
    energies(iepot_gas,:) = 0.0d0

    if (numbtot(1, ilevel) == 0) return

    ! Variables used below: rho, rho*v, E, non-thermal energies, B left/right
    nused = 0
    do ivar = 1, ndim + 1
      nused = nused + 1; ivar_used(nused) = ivar
    end do
    nused = nused + 1; ivar_used(nused) = neul
#if NENER>0
    do ivar = 1, nener
      nused = nused + 1; ivar_used(nused) = nhydro + ivar
    end do
#endif
#ifdef SOLVERmhd
    do ivar = 1, ndim
      nused = nused + 1; ivar_used(nused) = neul + ivar
      nused = nused + 1; ivar_used(nused) = nvar + ivar
    end do
#endif

    ! Mesh spacing at that level
    nx_loc = icoarse_max - icoarse_min + 1
    scale = boxlen/dble(nx_loc)
    if(ndim>0)skip_loc(1)=dble(icoarse_min)
    if(ndim>1)skip_loc(2)=dble(jcoarse_min)
    if(ndim>2)skip_loc(3)=dble(kcoarse_min)
    dx_code = 0.5D0**ilevel
    dx = dx_code*scale
    vol = dx**ndim

    ! Set position of cell centers relative to grid center
    do ind=1,twotondim
      iz=(ind-1)/4
      iy=(ind-1-4*iz)/2
      ix=(ind-1-2*iy-4*iz)
      if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx_code
      if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx_code
      if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx_code
    end do


    ! Loop over active grids by vector sweeps
    ncache = active(ilevel)%ngrid
    do igrid = 1, ncache, nvector
      ngrid = MIN(nvector, ncache - igrid + 1)
      do i = 1, ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
      end do

      ! Loop over cells
      do ind = 1, twotondim
        iskip = ncoarse + (ind - 1)*ngridmax
        do i = 1, ngrid
          ind_cell(i) = ind_grid(i) + iskip
        end do

        ! Gather leaf cells
        nleaf = 0
        do i = 1, ngrid
          if (son(ind_cell(i)) == 0) then
            nleaf = nleaf + 1
            ind_leaf(nleaf) = ind_cell(i)
            ind_grid_of_leaf(nleaf) = ind_grid(i)
          end if
        end do

        region_index = deltaE_nregions + 1
        if (deltaE_nregions > 0) then
          ! Compute cell center in code units
          do idim=1,ndim
            do i=1,nleaf
              xx(i,idim)=xg(ind_grid_of_leaf(i),idim)+xc(ind,idim)
            end do
          end do
          ! Rescale position from code units to user units
          do idim=1,ndim
            do i=1,nleaf
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
            end do
          end do

          do k = deltaE_nregions, 1, -1
            er = deltaE_region_exp       (k) ! Exponent defining norm
            xr = deltaE_region_x         (k) ! Region centre
            yr = deltaE_region_y         (k)
            zr = deltaE_region_z         (k)
            rr = deltaE_region_size  (k) ! Region diameter
            aa = deltaE_region_a         (k) ! Ellipticity (Y/X)
            bb = deltaE_region_b         (k) ! Ellipticity (Z/X)
            do i=1,nleaf
              xn = 0.0d0; yn=0.0d0; zn=0.0d0
              xn = abs(xx(i,1)-xr)
              xn = 2.0d0*xn/rr
#if NDIM > 1
              yn = abs(xx(i,2)-yr)
              yn = 2.0d0*yn/(aa*rr)
#endif
#if NDIM >2
              zn = abs(xx(i,3)-zr)
              zn = 2.0d0*zn/(bb*rr)
#endif
              if(er<10)then
                r = (xn**er+yn**er+zn**er)**(1/er)
              else
                r = max(xn,yn,zn)
              end if
              if (r < 1.0) region_index(i) = k
            end do
          end do
        end if


        ! Gather the hydro variables used below only
        do kv = 1, nused
          ivar = ivar_used(kv)
          if (deltaE_use_unew(ilevel)) then
            do i = 1, nleaf
              uu(i, ivar) = unew(ind_leaf(i), ivar)
            end do
          else
            do i = 1, nleaf
              uu(i, ivar) = uold(ind_leaf(i), ivar)
            end do
          end if
        end do
        do i = 1, nleaf
          rhoinv(i) = 1.0d0/uu(i, 1)
        end do


        if (poisson) then
          if (deltaE_use_phi_old(ilevel)) then
            do i = 1, nleaf
              phi_level(i) = phi_old(ind_leaf(i))
            end do
          else
            do i = 1, nleaf
              phi_level(i) = phi(ind_leaf(i))
            end do
          end if

          do i = 1, nleaf
            epot_leaf = 0.5*vol*uu(i, 1)*scale*phi_level(i)
            epot_loc(0) = epot_loc(0) + epot_leaf
            epot_loc(region_index(i)) = epot_loc(region_index(i)) + epot_leaf
          end do
        end if

        ! Compute total internal energy step 1
        do i = 1, nleaf
          eint_loc(0) = eint_loc(0) + uu(i, neul)*vol
          eint_loc(region_index(i)) = eint_loc(region_index(i)) + uu(i, neul)*vol
        end do

        ! Compute total energies
        do ivar = 1, ndim
          do i = 1, nleaf
            ! Compute total kinetic energy
            ekin_leaf = 0.5d0*vol*uu(i, 1 + ivar)**2*rhoinv(i)
            ekin_loc(0) = ekin_loc(0)  + ekin_leaf
            ekin_loc(region_index(i)) = ekin_loc(region_index(i)) + ekin_leaf
#ifdef SOLVERmhd
            ! Compute total magnetic energy
            emag_leaf = 0.125d0*(uu(i, neul + ivar) + uu(i, nvar + ivar))**2*vol
            emag_loc(0) = emag_loc(0) + emag_leaf
            emag_loc(region_index(i)) = emag_loc(region_index(i)) + emag_leaf
#endif

            ! Compute total internal energy step 2
            eint_loc(0) = eint_loc(0) - ekin_leaf - emag_leaf
            eint_loc(region_index(i)) = eint_loc(region_index(i)) - ekin_leaf - emag_leaf
          end do
        end do

        ! Compute total internal energy step 3
#if NENER>0
        do ivar = 1, nener
          do i = 1, nleaf
            eint_loc(0) = eint_loc(0) - uu(i, nhydro + ivar)*vol
            eint_loc(region_index(i)) = eint_loc(region_index(i)) - uu(i, nhydro + ivar)*vol
          end do
        end do
#endif

        if(pressure_fix .and. deltaE_correct_pressure_fix .and. deltaE_use_unew(ilevel)) then
          ! Correct internal energy if too small
          do i=1, nleaf
            ekin_leaf = 0.
            do ivar = 1, ndim
              ekin_leaf = ekin_leaf +  0.5d0*uu(i, 1 + ivar)**2*rhoinv(i)
            end do
#if NENER>0
            do ivar=1,nener
              ekin_leaf = ekin_leaf + uu(i, nhydro + ivar)
            end do
#endif
#ifdef SOLVERmhd
            do ivar = 1, ndim
              ekin_leaf = ekin_leaf + 0.125d0*(uu(i, neul + ivar) + uu(i, nvar + ivar))**2
            end do
#endif

            e_cons = uu(i, neul) - ekin_leaf
            e_prim = enew(ind_leaf(i))
            ! Note: here divu=-div.u*dt
            div = abs(divu(ind_leaf(i)))*dx/dtnew(ilevel)
            e_trunc = beta_fix*uu(i, 1)*max(div,3.0d0*hexp*dx)**2
            if(e_cons<e_trunc)then
                eint_loc(0) = eint_loc(0) - e_cons*vol + e_prim*vol
                eint_loc(region_index(i)) = eint_loc(region_index(i)) - e_cons*vol + e_prim*vol
            end if
          end do
        end if
      end do
      ! End loop over cells

    end do
    ! End loop over grids

    energies(iekin_gas,:) = ekin_loc(0:deltaE_nregions)
    energies(ieint,:) = eint_loc(0:deltaE_nregions)
    energies(iemag,:) = emag_loc(0:deltaE_nregions)
    energies(iepot_gas,:) = epot_loc(0:deltaE_nregions)
  end subroutine compute_energy_gas



  subroutine compute_ekin_turb_all_levels(eturb)
    ! Turbulent kinetic energy per level. For every grid at level L:
    !   eturb(L) += sum_i 1/2 m_i |v_i - V|^2,  V = mass-weighted mean velocity of the grid.
    ! Refined cells use the total mass and momentum of their child grid.
    use amr_commons
    use hydro_commons
    use mpi_mod
    implicit none
    real(dp), dimension(1:nlevelmax, 0:deltaE_nregions), intent(out) :: eturb
    real(dp), dimension(1:nlevelmax, 0:deltaE_nregions+1) :: eturb_loc
    integer, dimension(1:nvector) :: ind_grid
    integer :: ilevel, igrid, ngrid, ncache, i
#ifndef WITHOUTMPI
    integer :: info
#endif

    if (.not. allocated(grid_u)) allocate(grid_u(1:ngridmax, 1:ndim+1))

    eturb_loc = 0.0d0
    do ilevel = nlevelmax, levelmin, -1
      ncache = active(ilevel)%ngrid
      do igrid = 1, ncache, nvector
        ngrid = min(nvector, ncache - igrid + 1)
        do i = 1, ngrid
          ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
        end do
        call compute_eturb_aux(ilevel, ngrid, ind_grid, eturb_loc)
      end do
      ! Totals of level ilevel grids to their virtual copies, for parents owned by another rank
      if (ilevel > levelmin) call make_virtual_grid_totals(ilevel)
    end do

    eturb = eturb_loc(:, 0:deltaE_nregions)
  end subroutine compute_ekin_turb_all_levels

  subroutine compute_eturb_aux(ilevel, ngrid, ind_grid, eturb)
    ! For a batch of at most nvector grids owned by this rank at level ilevel: stores their total
    ! mass and momentum in grid_u and adds the fluctuation energy of their cells to eturb.
    use amr_commons
    use hydro_commons
    implicit none
    integer, intent(in) :: ilevel, ngrid
    integer, dimension(1:nvector), intent(in) :: ind_grid
    real(dp), dimension(1:nlevelmax,  0:deltaE_nregions+1), intent(inout) :: eturb

    real(dp), dimension(1:nvector, 1:twotondim) :: mass
    real(dp), dimension(1:nvector, 1:ndim, 1:twotondim) :: momentum
    real(dp), dimension(1:nvector) :: mass_total
    real(dp), dimension(1:nvector, 1:ndim) :: momentum_total

    integer, dimension(1:nvector) :: region_index = 0
    real(dp)::er,xr,yr,zr,rr,xn,yn,zn,r,aa,bb
    real(dp),dimension(1:3)::skip_loc
    real(dp),dimension(1:twotondim,1:3)::xc
    real(dp),dimension(1:nvector,1:ndim),save::xx
    real(dp):: dx_code, dx, vol, scale, eturb_cell
    integer :: nx_loc, k, ix, iy, iz


    integer :: i, ind, iskip, icell, ison, idim
    logical :: use_unew

    ! Mesh spacing at that level
    nx_loc = icoarse_max - icoarse_min + 1
    scale = boxlen/dble(nx_loc)
    if(ndim>0)skip_loc(1)=dble(icoarse_min)
    if(ndim>1)skip_loc(2)=dble(jcoarse_min)
    if(ndim>2)skip_loc(3)=dble(kcoarse_min)
    dx_code = 0.5D0**ilevel
    dx = dx_code*scale
    vol = dx**ndim

    ! Set position of cell centers relative to grid center
    do ind=1,twotondim
      iz=(ind-1)/4
      iy=(ind-1-4*iz)/2
      ix=(ind-1-2*iy-4*iz)
      if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx_code
      if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx_code
      if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx_code
    end do

    use_unew = deltaE_use_unew(ilevel)

    ! Leaves: own cells, from the current u array. Refined cells: totals of the child grid
    do ind = 1, twotondim
      iskip = ncoarse + (ind - 1)*ngridmax
      do i = 1, ngrid
        icell = ind_grid(i) + iskip
        ison = son(icell)
        if (ison > 0) then
          mass(i, ind) = grid_u(ison, 1)
          momentum(i, 1:ndim, ind) = grid_u(ison, 2:ndim+1)
        else if (use_unew) then
          mass(i, ind) = unew(icell, 1)*vol
          momentum(i, 1:ndim, ind) = unew(icell, 2:ndim+1)*vol
        else
          mass(i, ind) = uold(icell, 1)*vol
          momentum(i, 1:ndim, ind) = uold(icell, 2:ndim+1)*vol
        end if
      end do
    end do

    ! Grid totals, stored for the parent level
    mass_total(1:ngrid) = 0.0d0
    momentum_total(1:ngrid, 1:ndim) = 0.0d0
    do ind = 1, twotondim
      do i = 1, ngrid
        mass_total(i) = mass_total(i) + mass(i, ind)
      end do
      do idim = 1, ndim
        do i = 1, ngrid
          momentum_total(i, idim) = momentum_total(i, idim) + momentum(i, idim, ind)
        end do
      end do
    end do
    do i = 1, ngrid
      grid_u(ind_grid(i), 1) = mass_total(i)
    end do
    do idim = 1, ndim
      do i = 1, ngrid
        grid_u(ind_grid(i), idim + 1) = momentum_total(i, idim)
      end do
    end do


    do ind = 1, twotondim
      region_index = deltaE_nregions + 1
      if (deltaE_nregions > 0) then
        ! Compute cell center in code units
        do idim=1,ndim
          do i=1,ngrid
            xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
          end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
          do i=1,ngrid
            xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
          end do
        end do

        do k = deltaE_nregions, 1, -1
          er = deltaE_region_exp       (k) ! Exponent defining norm
          xr = deltaE_region_x         (k) ! Region centre
          yr = deltaE_region_y         (k)
          zr = deltaE_region_z         (k)
          rr = deltaE_region_size  (k) ! Region diameter
          aa = deltaE_region_a         (k) ! Ellipticity (Y/X)
          bb = deltaE_region_b         (k) ! Ellipticity (Z/X)
          do i=1,ngrid
            xn = 0.0d0; yn=0.0d0; zn=0.0d0
            xn = abs(xx(i,1)-xr)
            xn = 2.0d0*xn/rr
#if NDIM > 1
            yn = abs(xx(i,2)-yr)
            yn = 2.0d0*yn/(aa*rr)
#endif
#if NDIM >2
            zn = abs(xx(i,3)-zr)
            zn = 2.0d0*zn/(bb*rr)
#endif
            if(er<10)then
              r = (xn**er+yn**er+zn**er)**(1/er)
            else
              r = max(xn,yn,zn)
            end if
            if (r < 1.0) region_index(i) = k
          end do
        end do
      end if

    ! Fluctuations of the 8 cells around the grid mean
      do idim = 1, ndim
        do i = 1, ngrid
          eturb_cell =  0.5d0*mass(i, ind) &
          & *(momentum(i, idim, ind)/max(mass(i, ind), smallr*vol) &
          &   - momentum_total(i, idim)/max(mass_total(i), smallr*vol))**2
          eturb(ilevel, 0) = eturb(ilevel, 0) + eturb_cell
          eturb(ilevel, region_index(i)) = eturb(ilevel, region_index(i)) + eturb_cell
        end do
      end do
    end do
  end subroutine compute_eturb_aux

  subroutine make_virtual_grid_totals(ilevel)
    ! Sends grid_u of the grids owned at level ilevel to their virtual copies.
    ! Same emission/reception lists as make_virtual_fine_dp, with ndim+1 values per grid
    ! instead of one value per cell. Each cpu block is laid out like make_virtual_fine_dp:
    ! component ivar of grid i at position i + (ivar-1)*n.
    use amr_commons
    use mpi_mod
    implicit none
    integer, intent(in) :: ilevel
#ifndef WITHOUTMPI
    integer, parameter :: nv = ndim + 1, tag = 1101
    real(dp), allocatable, dimension(:) :: sendbuf, recvbuf
    integer, dimension(1:ncpu) :: nsend, nrecv, offsend, offrecv, reqsend, reqrecv
    integer, dimension(MPI_STATUS_SIZE, 1:ncpu) :: statuses
    integer :: icpu, i, ivar, step, countsend, countrecv, info
#ifdef LIGHT_MPI_COMM
    integer :: idx, offset
#endif

    if (numbtot(1, ilevel) == 0) return

    ! Number of grids sent to / received from each cpu, and offsets of the cpu blocks
    nsend = 0
#ifdef LIGHT_MPI_COMM
    do idx = 1, emission(ilevel)%nactive
      nsend(emission(ilevel)%cpuid(idx)) = emission(ilevel)%ngrids(idx)
    end do
#else
    do icpu = 1, ncpu
      nsend(icpu) = emission(icpu, ilevel)%ngrid
    end do
#endif
    do icpu = 1, ncpu
      nrecv(icpu) = reception(icpu, ilevel)%ngrid
    end do
    offsend(1) = 0
    offrecv(1) = 0
    do icpu = 2, ncpu
      offsend(icpu) = offsend(icpu - 1) + nv*nsend(icpu - 1)
      offrecv(icpu) = offrecv(icpu - 1) + nv*nrecv(icpu - 1)
    end do
    allocate(sendbuf(1:max(nv*sum(nsend), 1)), recvbuf(1:max(nv*sum(nrecv), 1)))

    ! Receive all messages
    countrecv = 0
    do icpu = 1, ncpu
      if (nrecv(icpu) > 0) then
        countrecv = countrecv + 1
        call MPI_IRECV(recvbuf(offrecv(icpu) + 1), nv*nrecv(icpu), MPI_DOUBLE_PRECISION, &
             & icpu - 1, tag, MPI_COMM_WORLD, reqrecv(countrecv), info)
      end if
    end do

    ! Gather the emission buffer
#ifdef LIGHT_MPI_COMM
    offset = 0
    do idx = 1, emission(ilevel)%nactive
      icpu = emission(ilevel)%cpuid(idx)
      do ivar = 1, nv
        step = offsend(icpu) + (ivar - 1)*nsend(icpu)
        do i = 1, nsend(icpu)
          sendbuf(step + i) = grid_u(emission(ilevel)%igrid(offset + i), ivar)
        end do
      end do
      offset = offset + nsend(icpu)
    end do
#else
    do icpu = 1, ncpu
      do ivar = 1, nv
        step = offsend(icpu) + (ivar - 1)*nsend(icpu)
        do i = 1, nsend(icpu)
          sendbuf(step + i) = grid_u(emission(icpu, ilevel)%igrid(i), ivar)
        end do
      end do
    end do
#endif

    ! Send all messages
    countsend = 0
    do icpu = 1, ncpu
      if (nsend(icpu) > 0) then
        countsend = countsend + 1
        call MPI_ISEND(sendbuf(offsend(icpu) + 1), nv*nsend(icpu), MPI_DOUBLE_PRECISION, &
             & icpu - 1, tag, MPI_COMM_WORLD, reqsend(countsend), info)
      end if
    end do

    ! Scatter the reception buffer into the virtual grids
    call MPI_WAITALL(countrecv, reqrecv, statuses, info)
    do icpu = 1, ncpu
      do ivar = 1, nv
        step = offrecv(icpu) + (ivar - 1)*nrecv(icpu)
        do i = 1, nrecv(icpu)
#ifdef LIGHT_MPI_COMM
          grid_u(reception(icpu, ilevel)%pcomm%igrid(i), ivar) = recvbuf(step + i)
#else
          grid_u(reception(icpu, ilevel)%igrid(i), ivar) = recvbuf(step + i)
#endif
        end do
      end do
    end do

    call MPI_WAITALL(countsend, reqsend, statuses, info)
    deallocate(sendbuf, recvbuf)
#endif
  end subroutine make_virtual_grid_totals


  subroutine compute_energy_part(ilevel, energies)
    use amr_commons
    use hydro_commons
    use mpi_mod
    use pm_commons

    implicit none

    ! Parameters & outputs
    integer, intent(in):: ilevel
    real(dp), dimension(1:nb_energy_kind, 0:deltaE_nregions), intent(inout) :: energies

    real(dp), dimension(-NFAMILIES:NFAMILIES,0:deltaE_nregions+1) :: ekin_families_loc
    real(dp), dimension(-NFAMILIES:NFAMILIES,0:deltaE_nregions+1) :: epot_families_loc

    real(kind=8), dimension(0:deltaE_nregions+1)::ekin_part_loc, epot_part_loc

    ! AMR variables
    integer :: igrid, jgrid

    ! Particle variable
    integer::ipart, jpart, next_part, ip, ig, npart1, local_counter

    integer, dimension(1:nvector), save:: ind_grid, ind_part, ind_grid_part

    ! MPI variables
#ifndef WITHOUTMPI
    integer::info
    real(kind=8), dimension(6,0:deltaE_nregions)::comm_buffin, comm_buffout
#endif


    energies(iepot_part,:) = 0.d0
    energies(iepot_part_dm,:) = 0.d0
    energies(iepot_part_star,:) = 0.d0
    energies(iekin_part,:) = 0.d0
    energies(iekin_part_dm,:) = 0.d0
    energies(iekin_part_star,:) = 0.d0

    ekin_part_loc = 0.d0
    epot_part_loc = 0.d0
    ekin_families_loc = 0.d0
    epot_families_loc = 0.d0

    if (pic) then
      ig = 0
      ip = 0
      ! Loop over particles that are not tracers (same buffering as move_fine)
      do jgrid = 1, active(ilevel)%ngrid
        igrid = active(ilevel)%igrid(jgrid)
        npart1 = numbp(igrid)  ! Number of particles in the grid
        if (npart1 > 0) then
          ig = ig + 1
          ind_grid(ig) = igrid
          ipart = headp(igrid)
          local_counter = 0
          ! Loop over particles
          do jpart = 1, npart1
            ! Save next particle  <---- Very important !!!
            next_part = nextp(ipart)
            if (ig == 0) then
              ig = 1
              ind_grid(ig) = igrid
            end if
            ! Skip tracers
            if (.not. is_tracer(typep(ipart))) then
              local_counter = local_counter + 1
              ip = ip + 1
              ind_part(ip) = ipart
              ind_grid_part(ip) = ig
              if (ip == nvector) then
                call energy_part_helper(ind_grid, ind_part, ind_grid_part, ig, ip, &
                     & ekin_part_loc, ekin_families_loc, epot_part_loc, epot_families_loc, ilevel)
                local_counter = 0
                ip = 0
                ig = 0
              end if
            end if

            ipart = next_part  ! Go to next particle
          end do
          ! End loop over particles
          ! If there was no particle in the grid, remove the grid from the buffer
          if (local_counter == 0 .and. ig > 0) ig = ig - 1
        end if
      end do
      ! End loop over grids
      if (ip > 0) then
        call energy_part_helper(ind_grid, ind_part, ind_grid_part, ig, ip, &
             & ekin_part_loc, ekin_families_loc, epot_part_loc, epot_families_loc, ilevel)
      end if

      energies(iekin_part,:) = ekin_part_loc(0:deltaE_nregions)
      energies(iekin_part_dm,:) = ekin_families_loc(FAM_DM,0:deltaE_nregions)
      energies(iekin_part_star,:) = ekin_families_loc(FAM_STAR,0:deltaE_nregions)
      energies(iepot_part,:) = epot_part_loc(0:deltaE_nregions)
      energies(iepot_part_dm,:) = epot_families_loc(FAM_DM,0:deltaE_nregions)
      energies(iepot_part_star,:) = epot_families_loc(FAM_STAR,0:deltaE_nregions)

    end if

  end subroutine compute_energy_part

  subroutine energy_part_helper(ind_grid, ind_part, ind_grid_part, ng, nn, &
    & ekin_loc, ekin_families_loc, epot_loc, epot_families_loc, ilevel)
    use amr_commons
    use pm_commons
    use hydro_commons
    implicit none

    integer::ng, nn, ilevel
    integer, dimension(1:nvector)::ind_grid, ind_part, ind_grid_part

    integer::i, idim, nx_loc, k
    real(dp)::scale, epot_part, ekin_part
    real(dp), dimension(1:nvector), save:: mass_part, phi_part
    integer, dimension(1:nvector) :: type_part
    real(dp), dimension(1:nvector, 1:ndim)::vel_part
    real(dp), dimension(-NFAMILIES:NFAMILIES, 0:1 + deltaE_nregions), intent(inout) :: ekin_families_loc
    real(dp), dimension(-NFAMILIES:NFAMILIES, 0:1 + deltaE_nregions), intent(inout) :: epot_families_loc
    real(dp), intent(inout), dimension(0:1 + deltaE_nregions):: ekin_loc, epot_loc
    integer, dimension(1:nvector) :: region_index
    real(dp)::er, xr, yr, zr, rr, xn, yn, zn, r, aa, bb

    real(dp)::dx

    ! Compute time step
    dx = 0.5D0**ilevel
    nx_loc = (icoarse_max - icoarse_min + 1)
    scale = boxlen/dble(nx_loc)

    do idim = 1, ndim
      do i = 1, nn
        vel_part(i, idim) = vp(ind_part(i), idim)
      end do
    end do

    ! Get region index
    region_index = deltaE_nregions + 1
    if (deltaE_nregions > 0) then
      do k = deltaE_nregions, 1, -1
        er = deltaE_region_exp(k) ! Exponent defining norm
        xr = deltaE_region_x(k) ! Region centre
        yr = deltaE_region_y(k)
        zr = deltaE_region_z(k)
        rr = deltaE_region_size(k) ! Region size (diameter for a circle)
        aa = deltaE_region_a(k) ! Ellipticity (Y/X)
        bb = deltaE_region_b(k) ! Ellipticity (Z/X)
        do i = 1, nn
          xn = 0.0d0; yn = 0.0d0; zn = 0.0d0
          xn = abs(xp(ind_part(i), 1) - xr)
          xn = 2.0d0*xn/rr
#if NDIM > 1
          yn = abs(xp(ind_part(i), 2) - yr)
          yn = 2.0d0*yn/(aa*rr)
#endif
#if NDIM >2
          zn = abs(xp(ind_part(i), 3) - zr)
          zn = 2.0d0*zn/(bb*rr)
#endif
          if (er < 10) then
            r = (xn**er + yn**er + zn**er)**(1/er)
          else
            r = max(xn, yn, zn)
          end if
          if (r < 1.0) region_index(i) = k
        end do
      end do
    end if

    ! Fetch mass and type
    do i = 1, nn
      mass_part(i) = mp(ind_part(i))
      type_part(i) = typep(ind_part(i))%family
    end do

    ! Potential at the current particle position, from the current grid potential
    if (poisson) then
      call cic_phi_part(ind_grid, ind_part, ind_grid_part, ng, nn, ilevel, phi_part)
      do i = 1, nn
        epot_part = 0.5*mass_part(i)*phi_part(i)*scale
      epot_loc(0) = epot_loc(0) + epot_part
      epot_loc(region_index(i)) = epot_loc(region_index(i)) + epot_part
      epot_families_loc(type_part(i), 0) = epot_families_loc(type_part(i), 0) + epot_part
      epot_families_loc(type_part(i), region_index(i)) = epot_families_loc(type_part(i), region_index(i)) + epot_part
    end do
    end if

    ! Compute kinetic energy
    do idim = 1, ndim
      do i = 1, nn
        ekin_part = 0.5D0*mass_part(i)*vel_part(i, idim)**2
        ekin_loc(0) = ekin_loc(0) + ekin_part
        ekin_loc(region_index(i)) = ekin_loc(region_index(i)) + ekin_part
        ekin_families_loc(type_part(i), 0) = ekin_families_loc(type_part(i), 0) + ekin_part
        ekin_families_loc(type_part(i), region_index(i)) = ekin_families_loc(type_part(i), region_index(i)) + ekin_part
      end do
    end do

  end subroutine energy_part_helper

  subroutine cic_phi_part(ind_grid, ind_part, ind_grid_part, ng, np, ilevel, phi_part)
    !------------------------------------------------------------
    ! CIC interpolation of the potential at the particle positions,
    ! Cloud at level ilevel if the cloud lies
    ! entirely in it, level ilevel-1 otherwise.
    ! Uses phi_old instead of phi for the cells of a level whose
    ! deltaE_use_phi_old flag is set.
    !------------------------------------------------------------
    use amr_commons
    use pm_commons
    use poisson_commons
    implicit none
    integer, intent(in)::ng, np, ilevel
    integer, dimension(1:nvector), intent(in)::ind_grid, ind_grid_part, ind_part
    real(dp), dimension(1:nvector), intent(out)::phi_part

    integer::j, ind, idim, nx_loc
    real(dp)::dx, scale
    logical::use_old_fine, use_old_coarse
    integer, dimension(1:nvector), save::father_cell
    real(dp), dimension(1:nvector, 1:ndim), save::x0
    integer, dimension(1:nvector, 1:threetondim), save::nbors_father_cells
    logical, dimension(1:nvector), save::ok
    real(dp), dimension(1:nvector, 1:ndim), save::x, dd, dg
    integer, dimension(1:nvector, 1:ndim), save::ig, id, igg, igd, icg, icd
    real(dp), dimension(1:nvector, 1:twotondim), save::vol
    integer, dimension(1:nvector, 1:twotondim), save::igrid, icell, indp, kg
    real(dp), dimension(1:3)::skip_loc

    dx = 0.5D0**ilevel
    nx_loc = (icoarse_max - icoarse_min + 1)
    skip_loc = (/0.0d0, 0.0d0, 0.0d0/)
    if (ndim > 0) skip_loc(1) = dble(icoarse_min)
    if (ndim > 1) skip_loc(2) = dble(jcoarse_min)
    if (ndim > 2) skip_loc(3) = dble(kcoarse_min)
    scale = boxlen/dble(nx_loc)

    ! Lower left corner of 3x3x3 grid-cube
    do idim = 1, ndim
      do j = 1, ng
        x0(j, idim) = xg(ind_grid(j), idim) - 3.0D0*dx
      end do
    end do
    do j = 1, ng
      father_cell(j) = father(ind_grid(j))
    end do
    call get3cubefather(father_cell, nbors_father_cells, ng, ilevel)

    ! Particle position in units of dx relative to the 3x3x3 cube
    do idim = 1, ndim
      do j = 1, np
        x(j, idim) = (xp(ind_part(j), idim)/scale + skip_loc(idim) - x0(ind_grid_part(j), idim))/dx
      end do
    end do

    ! CIC at level ilevel
    do idim = 1, ndim
      do j = 1, np
        dd(j, idim) = x(j, idim) + 0.5D0
        id(j, idim) = int(dd(j, idim))
        dd(j, idim) = dd(j, idim) - id(j, idim)
        dg(j, idim) = 1.0D0 - dd(j, idim)
        ig(j, idim) = id(j, idim) - 1
        igg(j, idim) = ig(j, idim)/2
        igd(j, idim) = id(j, idim)/2
      end do
    end do
#if NDIM==1
    do j = 1, np
      kg(j, 1) = 1 + igg(j, 1)
      kg(j, 2) = 1 + igd(j, 1)
    end do
#endif
#if NDIM==2
    do j = 1, np
      kg(j, 1) = 1 + igg(j, 1) + 3*igg(j, 2)
      kg(j, 2) = 1 + igd(j, 1) + 3*igg(j, 2)
      kg(j, 3) = 1 + igg(j, 1) + 3*igd(j, 2)
      kg(j, 4) = 1 + igd(j, 1) + 3*igd(j, 2)
    end do
#endif
#if NDIM==3
    do j = 1, np
      kg(j, 1) = 1 + igg(j, 1) + 3*igg(j, 2) + 9*igg(j, 3)
      kg(j, 2) = 1 + igd(j, 1) + 3*igg(j, 2) + 9*igg(j, 3)
      kg(j, 3) = 1 + igg(j, 1) + 3*igd(j, 2) + 9*igg(j, 3)
      kg(j, 4) = 1 + igd(j, 1) + 3*igd(j, 2) + 9*igg(j, 3)
      kg(j, 5) = 1 + igg(j, 1) + 3*igg(j, 2) + 9*igd(j, 3)
      kg(j, 6) = 1 + igd(j, 1) + 3*igg(j, 2) + 9*igd(j, 3)
      kg(j, 7) = 1 + igg(j, 1) + 3*igd(j, 2) + 9*igd(j, 3)
      kg(j, 8) = 1 + igd(j, 1) + 3*igd(j, 2) + 9*igd(j, 3)
    end do
#endif
    do ind = 1, twotondim
      do j = 1, np
        igrid(j, ind) = son(nbors_father_cells(ind_grid_part(j), kg(j, ind)))
      end do
    end do

    ! Cloud entirely in level ilevel ?
    ok(1:np) = .true.
    do ind = 1, twotondim
      do j = 1, np
        ok(j) = ok(j) .and. igrid(j, ind) > 0
      end do
    end do

    ! If not, redo CIC at level ilevel-1
    do idim = 1, ndim
      do j = 1, np
        if (.not. ok(j)) then
          x(j, idim) = x(j, idim)/2.0D0
          dd(j, idim) = x(j, idim) + 0.5D0
          id(j, idim) = int(dd(j, idim))
          dd(j, idim) = dd(j, idim) - id(j, idim)
          dg(j, idim) = 1.0D0 - dd(j, idim)
          ig(j, idim) = id(j, idim) - 1
        end if
      end do
    end do

    ! Parent cell position
    do idim = 1, ndim
      do j = 1, np
        if (ok(j)) then
          icg(j, idim) = ig(j, idim) - 2*igg(j, idim)
          icd(j, idim) = id(j, idim) - 2*igd(j, idim)
        else
          icg(j, idim) = ig(j, idim)
          icd(j, idim) = id(j, idim)
        end if
      end do
    end do
#if NDIM==1
    do j = 1, np
      icell(j, 1) = 1 + icg(j, 1)
      icell(j, 2) = 1 + icd(j, 1)
    end do
#endif
#if NDIM==2
    do j = 1, np
      if (ok(j)) then
        icell(j, 1) = 1 + icg(j, 1) + 2*icg(j, 2)
        icell(j, 2) = 1 + icd(j, 1) + 2*icg(j, 2)
        icell(j, 3) = 1 + icg(j, 1) + 2*icd(j, 2)
        icell(j, 4) = 1 + icd(j, 1) + 2*icd(j, 2)
      else
        icell(j, 1) = 1 + icg(j, 1) + 3*icg(j, 2)
        icell(j, 2) = 1 + icd(j, 1) + 3*icg(j, 2)
        icell(j, 3) = 1 + icg(j, 1) + 3*icd(j, 2)
        icell(j, 4) = 1 + icd(j, 1) + 3*icd(j, 2)
      end if
    end do
#endif
#if NDIM==3
    do j = 1, np
      if (ok(j)) then
        icell(j, 1) = 1 + icg(j, 1) + 2*icg(j, 2) + 4*icg(j, 3)
        icell(j, 2) = 1 + icd(j, 1) + 2*icg(j, 2) + 4*icg(j, 3)
        icell(j, 3) = 1 + icg(j, 1) + 2*icd(j, 2) + 4*icg(j, 3)
        icell(j, 4) = 1 + icd(j, 1) + 2*icd(j, 2) + 4*icg(j, 3)
        icell(j, 5) = 1 + icg(j, 1) + 2*icg(j, 2) + 4*icd(j, 3)
        icell(j, 6) = 1 + icd(j, 1) + 2*icg(j, 2) + 4*icd(j, 3)
        icell(j, 7) = 1 + icg(j, 1) + 2*icd(j, 2) + 4*icd(j, 3)
        icell(j, 8) = 1 + icd(j, 1) + 2*icd(j, 2) + 4*icd(j, 3)
      else
        icell(j, 1) = 1 + icg(j, 1) + 3*icg(j, 2) + 9*icg(j, 3)
        icell(j, 2) = 1 + icd(j, 1) + 3*icg(j, 2) + 9*icg(j, 3)
        icell(j, 3) = 1 + icg(j, 1) + 3*icd(j, 2) + 9*icg(j, 3)
        icell(j, 4) = 1 + icd(j, 1) + 3*icd(j, 2) + 9*icg(j, 3)
        icell(j, 5) = 1 + icg(j, 1) + 3*icg(j, 2) + 9*icd(j, 3)
        icell(j, 6) = 1 + icd(j, 1) + 3*icg(j, 2) + 9*icd(j, 3)
        icell(j, 7) = 1 + icg(j, 1) + 3*icd(j, 2) + 9*icd(j, 3)
        icell(j, 8) = 1 + icd(j, 1) + 3*icd(j, 2) + 9*icd(j, 3)
      end if
    end do
#endif

    ! Parent cell addresses
    do ind = 1, twotondim
      do j = 1, np
        if (ok(j)) then
          indp(j, ind) = ncoarse + (icell(j, ind) - 1)*ngridmax + igrid(j, ind)
        else
          indp(j, ind) = nbors_father_cells(ind_grid_part(j), icell(j, ind))
        end if
      end do
    end do

    ! Cloud volumes
#if NDIM==1
    do j = 1, np
      vol(j, 1) = dg(j, 1)
      vol(j, 2) = dd(j, 1)
    end do
#endif
#if NDIM==2
    do j = 1, np
      vol(j, 1) = dg(j, 1)*dg(j, 2)
      vol(j, 2) = dd(j, 1)*dg(j, 2)
      vol(j, 3) = dg(j, 1)*dd(j, 2)
      vol(j, 4) = dd(j, 1)*dd(j, 2)
    end do
#endif
#if NDIM==3
    do j = 1, np
      vol(j, 1) = dg(j, 1)*dg(j, 2)*dg(j, 3)
      vol(j, 2) = dd(j, 1)*dg(j, 2)*dg(j, 3)
      vol(j, 3) = dg(j, 1)*dd(j, 2)*dg(j, 3)
      vol(j, 4) = dd(j, 1)*dd(j, 2)*dg(j, 3)
      vol(j, 5) = dg(j, 1)*dg(j, 2)*dd(j, 3)
      vol(j, 6) = dd(j, 1)*dg(j, 2)*dd(j, 3)
      vol(j, 7) = dg(j, 1)*dd(j, 2)*dd(j, 3)
      vol(j, 8) = dd(j, 1)*dd(j, 2)*dd(j, 3)
    end do
#endif

    ! Interpolate the potential: the cloud cells are at level ilevel if ok,
    ! at level ilevel-1 otherwise
    use_old_fine = deltaE_use_phi_old(ilevel)
    use_old_coarse = .false.
    if (ilevel > levelmin) use_old_coarse = deltaE_use_phi_old(ilevel - 1)
    phi_part(1:np) = 0.0d0
    do ind = 1, twotondim
      do j = 1, np
        if ((ok(j) .and. use_old_fine) .or. (.not. ok(j) .and. use_old_coarse)) then
          phi_part(j) = phi_part(j) + phi_old(indp(j, ind))*vol(j, ind)
        else
          phi_part(j) = phi_part(j) + phi(indp(j, ind))*vol(j, ind)
        end if
      end do
    end do

  end subroutine cic_phi_part

  !#########################################################################

end module deltaE_module
