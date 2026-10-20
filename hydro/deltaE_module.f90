
module deltaE_module

  use amr_commons

  implicit none

  ! DeltaE params
  logical::deltaE_enable=.true. ! Whether to enable deltaE output
  logical::deltaE_correct_pressure_fix=.true. ! Whether to be "pressure_fix" aware when computing energy
  logical::deltaE_debug=.false. ! print intermediate energies
  logical::deltaE_force_all_levels=.false. ! Force computation on all levels for all processes


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

  type process
  real(dp), dimension(1:nb_energy_kind) :: v
  character(len=20) :: name
  logical :: compute_eturb
end type

! Process indices. To add a process: add its index here and one line in initialize_processes.
enum, bind(c)
  enumerator :: iflux_gas = 1, iflux_part, icooling, igravity_gas, igravity_part, &
              & istar_formation, ifeedback, iturb_driving, ipoisson, &
              & imagnetic_diffusion, irefine, icorrections, nb_process_end
end enum
integer, parameter :: nb_process = nb_process_end - 1

type(process), dimension(1:nb_process) :: deltaE

contains

  subroutine read_deltaE_params(namelist_unit,nml_ok)
   
    integer,intent(in)::namelist_unit
    logical,intent(inout)::nml_ok
    integer::nml_err
   
    namelist/deltaE_params/deltaE_correct_pressure_fix,deltaE_debug,deltaE_enable,deltaE_force_all_levels
   ! Go to the beginning of the file
    rewind(namelist_unit)

    ! Read namelist
    read(namelist_unit,NML=deltaE_params,IOSTAT=nml_err)

    if(nml_err>0)then
      if(myid==1)write(*,*)'Error reading namelist &DELTAE_PARAMS. Check formatting.'
      nml_ok=.false.
    end if

    allocate(deltaE_use_unew(levelmin:nlevelmax))
    allocate(deltaE_use_phi_old(levelmin:nlevelmax))

    deltaE_use_phi_old = .false.
    deltaE_use_unew = .false.

    call initialize_processes

  end subroutine read_deltaE_params

  subroutine initialize_processes
    deltaE(iflux_gas)           = process(0.0d0, "flux_gas",           .true.)
    deltaE(iflux_part)          = process(0.0d0, "flux_part",          .false.)
    deltaE(icooling)            = process(0.0d0, "cooling",            .false.)
    deltaE(igravity_gas)        = process(0.0d0, "gravity_gas",        .true.)
    deltaE(igravity_part)       = process(0.0d0, "gravity_part",       .false.)
    deltaE(istar_formation)     = process(0.0d0, "star_formation",     .true.)
    deltaE(ifeedback)           = process(0.0d0, "feedback",           .true.)
    deltaE(iturb_driving)       = process(0.0d0, "turb_driving",       .true.)
    deltaE(ipoisson)            = process(0.0d0, "poisson",            .false.)
    deltaE(imagnetic_diffusion) = process(0.0d0, "magnetic_diffusion", .false.)
    deltaE(irefine)             = process(0.0d0, "refine",             .true.)
    deltaE(icorrections)        = process(0.0d0, "corrections",        .true.)
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

  logical function is_energy_printed(i)
    integer, intent(in) :: i
    is_energy_printed = (i <= nb_energy_fixed) .or. &
         & (i - iekin_turb0 >= levelmin .and. i - iekin_turb0 <= nlevelmax)
  end function

  subroutine print_energies(energies)
    implicit none
    real(dp), dimension(1:nb_energy_kind), intent(in) :: energies
    integer :: i

    write (*, *) "Energy output"
    do i = 1, nb_energy_kind
      if (.not. is_energy_printed(i)) cycle
      write (*, 996) " ", energy_name(i), energies(i)
    end do
    write (*, *) "End energy output"

996 format(A1, A20, 1pe20.7)

  end subroutine

  subroutine print_processes
    implicit none
    integer :: i, j

    write (*, *) "DeltaE output"
    write (*, 997) " ", (deltaE(j)%name, j = 1, nb_process)
    do i = 1, nb_energy_kind
      if (.not. is_energy_printed(i)) cycle
      write (*, 998) energy_name(i), (deltaE(j)%v(i), j = 1, nb_process)
    end do
    write (*, *) "End deltaE output"

997 format(A20, *(A20))
998 format(A20, *(1pe20.5))
  end subroutine

  subroutine compute_energies(ilevel, energies, compute_eturb)
    implicit none

    integer, intent(in):: ilevel
    real(dp), dimension(1:nb_energy_kind), intent(inout) :: energies
    logical, intent(in) :: compute_eturb
    real(dp), dimension(1:nlevelmax) :: eturb

    energies = 0.0d0
    if (hydro) then
      call compute_energy_gas(ilevel, energies)
    end if
    if (pic) then
      call compute_energy_part(ilevel, energies)
    end if
    if (hydro .and. compute_eturb) then
      eturb = 0.0d0
      call compute_ekin_turb_all_levels(eturb)
      energies(iekin_turb0 + levelmin:iekin_turb0 + nlevelmax) = eturb(levelmin:nlevelmax)
    end if
    energies(iekin) = energies(iekin_gas) + energies(iekin_part)
    energies(iepot) = energies(iepot_gas) + energies(iepot_part)
  end subroutine

  subroutine compute_transfer(levelstart, levelend, deltaE_process, step)

    implicit none
    integer, intent(in):: levelstart, levelend
    integer :: act_levelstart, act_levelend
    type(process), intent(inout) :: deltaE_process
    integer, intent(in) :: step

    integer :: ilevel

    real(dp), dimension(1:nb_energy_kind), save :: energy_before
    real(dp), dimension(1:nb_energy_kind) ::energy_level, energy_after

    if (deltaE_force_all_levels) then
      act_levelstart = levelmin
      act_levelend = nlevelmax
    else 
      act_levelstart = levelstart
      act_levelend = min(levelend + 1, nlevelmax)
    end if

    if (step == 1) then
      energy_before = 0.0d0
    else
      energy_after = 0.0d0
    end if

    do ilevel = act_levelstart, act_levelend
      call compute_energies(ilevel, energy_level, ilevel == act_levelstart .and. deltaE_process%compute_eturb)
      if (step == 1) then
        energy_before = energy_before + energy_level
      else
        energy_after = energy_after + energy_level
      end if
    end do


    if (step == 2) then

      deltaE_process%v = deltaE_process%v + energy_after - energy_before

      if (deltaE_debug .and. myid == 1) then
        write(*,*) "DeltaE debug before ", deltaE_process%name
        call print_energies(energy_before)
        write(*,*) "DeltaE debug after ", deltaE_process%name
        call print_energies(energy_after)
      end if
    end if

  end subroutine

  subroutine compute_energy_gas(ilevel, energies)
    use amr_commons
    use hydro_commons
    use poisson_commons
    use mpi_mod
    implicit none
#ifndef WITHOUTMPI
    integer::info
    real(kind=8), dimension(4)::comm_buffin, comm_buffout
#endif
    integer, intent(in)::ilevel
    real(dp), dimension(1:nb_energy_kind), intent(inout) :: energies

    integer::i, ivar,  ind, ncache, igrid, iskip
    integer::nleaf, ngrid, nx_loc
    integer, dimension(1:nvector), save::ind_grid, ind_cell, ind_leaf

    real(dp):: dx, vol, scale
    real(kind=8)::mass_loc, ekin_loc, eint_loc, emag_loc, epot_loc, ekin_leaf, emag_leaf, epot_leaf
    real(kind=8)::mass_all
    real(dp), dimension(1:nvector, 1:nvar_all), save::uu
    real(dp), dimension(1:nvector), save::phi_level

    real(dp), dimension(1:nvector, 1:ndim), save::gg

    real(dp):: e_cons, e_prim, e_trunc, div

    mass_loc = 0.0d0
    ekin_loc = 0.0d0; ekin_leaf = 0.0d0
    emag_loc = 0.0d0; emag_leaf = 0.0d0
    eint_loc = 0.0d0
    epot_loc = 0.0d0; epot_leaf = 0.0d0

    energies(iekin_gas) = 0.0d0
    energies(ieint) = 0.0d0
    energies(iemag) = 0.0d0
    energies(iepot_gas) = 0.0d0

    if (numbtot(1, ilevel) == 0) return

    ! Mesh spacing at that level
    nx_loc = icoarse_max - icoarse_min + 1
    scale = boxlen/dble(nx_loc)
    dx = 0.5D0**ilevel*scale
    vol = dx**ndim

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
          end if
        end do

        ! Gather hydro variables
        do ivar = 1, nvar_all
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
            epot_loc = epot_loc + epot_leaf
          end do
        end if

        ! Compute total internal energy step 1
        do i = 1, nleaf
          eint_loc = eint_loc + uu(i, neul)*vol
        end do

        ! Compute total energies
        do ivar = 1, ndim
          do i = 1, nleaf
            ! Compute total kinetic energy
            ekin_leaf = 0.5d0*vol*uu(i, 1 + ivar)**2/uu(i, 1)
            ekin_loc = ekin_loc + ekin_leaf
#ifdef SOLVERmhd
            ! Compute total magnetic energy
            emag_leaf = 0.125d0*(uu(i, neul + ivar) + uu(i, nvar + ivar))**2*vol
            emag_loc = emag_loc + emag_leaf
#endif

            ! Compute total internal energy step 2
            eint_loc = eint_loc - ekin_leaf - emag_leaf
          end do
        end do

        ! Compute total internal energy step 3 
#if NENER>0
        do ivar = 1, nener
          do i = 1, nleaf
            eint_loc = eint_loc - uu(i, nhydro + ivar)*vol
          end do
        end do
#endif

        if(pressure_fix .and. deltaE_correct_pressure_fix .and. deltaE_use_unew(ilevel)) then
          ! Correct internal energy if too small
          do i=1, nleaf
            ekin_leaf = 0.
            do ivar = 1, ndim
              ekin_leaf = ekin_leaf +  0.5d0*uu(i, 1 + ivar)**2/uu(i, 1)
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
                eint_loc = eint_loc - e_cons*vol + e_prim*vol
            end if
          end do
        end if
      end do
      ! End loop over cells

    end do
    ! End loop over grids

    ! Compute global quantities
#ifndef WITHOUTMPI
    comm_buffin(1) = ekin_loc
    comm_buffin(2) = eint_loc
    comm_buffin(3) = emag_loc
    comm_buffin(4) = epot_loc
    call MPI_ALLREDUCE(comm_buffin, comm_buffout, 4, MPI_DOUBLE_PRECISION, MPI_SUM,&
            &MPI_COMM_WORLD, info)
    energies(iekin_gas) = comm_buffout(1)
    energies(ieint) = comm_buffout(2)
    energies(iemag) = comm_buffout(3)
    energies(iepot_gas) = comm_buffout(4)
#endif
#ifdef WITHOUTMPI
    energies(iekin_gas) = ekin_loc
    energies(ieint) = eint_loc
    energies(iemag) = emag_loc
    energies(iepot_gas) = epot_loc
#endif
  end subroutine compute_energy_gas


  subroutine compute_ekin_turb_all_levels(eturb)
    ! Turbulent kinetic energy per level. For every grid i at level L:
    !   eturb(L) += sum_i 1/2 m_i |v_i - V|^2,  V = mass-weighted mean velocity of the grid.
    ! Refined cells use the total mass and momentum of their child grid,
    ! since some updates are done only on leaf cells (eg. feedback, star formation, etc.)
    use amr_commons
    use hydro_commons
    use mpi_mod
    implicit none
    real(dp), dimension(1:nlevelmax), intent(out) :: eturb
    real(dp), dimension(1:nlevelmax) :: eturb_loc
    integer, dimension(1:nvector) :: ind_grid
    real(dp), dimension(1:nvector) :: mass_total
    real(dp), dimension(1:ndim, 1:nvector) :: momentum_total
    integer :: igrid, ngrid, ncache, i
#ifndef WITHOUTMPI
    integer :: info
#endif

    eturb_loc = 0.0d0
    ncache = active(levelmin)%ngrid
    do igrid = 1, ncache, nvector
      ngrid = min(nvector, ncache - igrid + 1)
      do i = 1, ngrid
        ind_grid(i) = active(levelmin)%igrid(igrid + i - 1)
      end do
      call compute_eturb_aux(levelmin, ngrid, ind_grid, mass_total, momentum_total, eturb_loc)
    end do

#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(eturb_loc, eturb, nlevelmax, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#else
    eturb = eturb_loc
#endif
  end subroutine compute_ekin_turb_all_levels

  recursive subroutine compute_eturb_aux(ilevel, ngrid, ind_grid, mass_total, momentum_total, eturb)
    ! For a batch of at most nvector grids at level ilevel: returns their total mass and
    ! momentum, and adds the fluctuation energy of their cells around the grid mean to
    ! eturb. Refined cells are handled by recursive calls on vectors of at most nvector grids.
    use amr_commons
    use hydro_commons
    implicit none
    integer, intent(in) :: ilevel, ngrid
    integer, dimension(1:ngrid), intent(in) :: ind_grid
    real(dp), dimension(1:ngrid), intent(out) :: mass_total
    real(dp), dimension(1:ndim, 1:ngrid), intent(out) :: momentum_total
    real(dp), dimension(1:nlevelmax), intent(inout) :: eturb

    real(dp), dimension(1:nvector, 1:twotondim) :: mass
    real(dp), dimension(1:ndim, 1:ngrid, 1:twotondim) :: momentum
    integer, dimension(1:twotondim*nvector) :: child, child_igrid, child_ind
    real(dp), dimension(1:twotondim*nvector) :: child_mass
    real(dp), dimension(1:ndim, 1:twotondim*nvector) :: child_momentum
    real(dp), dimension(1:ndim) :: vmean
    real(dp) :: vol
    integer :: igrid, ind, iskip, icell, nchild, k, nk, idim
    logical :: use_unew

    vol = (0.5d0**ilevel*boxlen/dble(icoarse_max - icoarse_min + 1))**ndim
    use_unew = deltaE_use_unew(ilevel)

    ! Leaves: take from the current u array. 
    ! Refined cells: takes from their child
    nchild = 0
    do ind = 1, twotondim
      iskip = ncoarse + (ind - 1)*ngridmax
      do igrid = 1, ngrid
        icell = ind_grid(igrid) + iskip
        if (son(icell) > 0) then
          nchild = nchild + 1
          child(nchild) = son(icell)
          child_igrid(nchild) = igrid
          child_ind(nchild) = ind
        else if (use_unew) then
          mass(igrid, ind) = unew(icell, 1)*vol
          momentum(:, igrid, ind) = unew(icell, 2:ndim+1)*vol
        else
          mass(igrid, ind) = uold(icell, 1)*vol
          momentum(:, igrid, ind) = uold(icell, 2:ndim+1)*vol
        end if
      end do
    end do

    ! Recursive call to level ilevel + 1 only for children of this vector
    do k = 1, nchild, nvector
      nk = min(nvector, nchild - k + 1)
      call compute_eturb_aux(ilevel + 1, nk, child(k:k+nk-1), child_mass(k:k+nk-1), child_momentum(:, k:k+nk-1), eturb)
    end do
    do k = 1, nchild
      mass(child_igrid(k), child_ind(k)) = child_mass(k)
      momentum(:, child_igrid(k), child_ind(k)) = child_momentum(:, k)
    end do

    ! Grid totals, then fluctuations of the 8 cells around the grid mean
    do igrid = 1, ngrid
      mass_total(igrid) = 0.0d0
      momentum_total(:, igrid) = 0.0d0
    end do
    do ind = 1, twotondim
      do igrid = 1, ngrid
        mass_total(igrid) = mass_total(igrid) + mass(igrid, ind)
        momentum_total(:, igrid) = momentum_total(:, igrid) + momentum(:, igrid, ind)
      end do
    end do
    do igrid = 1, ngrid
      vmean = momentum_total(:, igrid)/max(mass_total(igrid), smallr*vol)
      do ind = 1, twotondim
        do idim = 1, ndim
          eturb(ilevel) = eturb(ilevel) + 0.5d0*mass(igrid, ind) &
               & *(momentum(idim, igrid, ind)/max(mass(igrid, ind), smallr*vol) - vmean(idim))**2
        end do
      end do
    end do
  end subroutine compute_eturb_aux


  subroutine compute_energy_part(ilevel, energies)
    use amr_commons
    use hydro_commons
    use mpi_mod
    use pm_commons

    implicit none

    ! Parameters & outputs
    integer, intent(in):: ilevel
    real(dp), dimension(1:nb_energy_kind), intent(inout) :: energies

    real(dp), dimension(-NFAMILIES:NFAMILIES) :: ekin_families_loc
    real(dp), dimension(-NFAMILIES:NFAMILIES) :: epot_families_loc

    real(kind=8)::ekin_part_loc, epot_part_loc

    ! AMR variables
    integer :: igrid, jgrid

    ! Particle variable
    integer::ipart, jpart, next_part, ip, npart1, local_counter

    integer, dimension(1:nvector), save:: ind_part

    ! MPI variables
#ifndef WITHOUTMPI
    integer::info
#endif

    energies(iepot_part) = 0.d0
    energies(iepot_part_dm) = 0.d0
    energies(iepot_part_star) = 0.d0
    energies(iekin_part) = 0.d0
    energies(iekin_part_dm) = 0.d0
    energies(iekin_part_star) = 0.d0

    ekin_part_loc = 0.d0
    epot_part_loc = 0.d0
    ekin_families_loc = 0.d0
    epot_families_loc = 0.d0

    if (pic) then
      ip = 0
      ! Loop over particles that are not tracers
      do jgrid = 1, active(ilevel)%ngrid
        igrid = active(ilevel)%igrid(jgrid)
        npart1 = numbp(igrid)  ! Number of particles in the grid
        if (npart1 > 0) then
          ipart = headp(igrid)
          ! Loop over particles
          do jpart = 1, npart1
            ! Save next particle  <---- Very important !!!
            next_part = nextp(ipart)
            ! Skip tracers
            if (.not. is_tracer(typep(ipart))) then
              ip = ip + 1
              ind_part(ip) = ipart
              if (ip == nvector) then
                call ekin_part_helper(ind_part, ekin_part_loc, ekin_families_loc, epot_part_loc, epot_families_loc, ip, ilevel)
                ip = 0
              end if
            end if

            ipart = next_part  ! Go to next particle
          end do
          ! End loop over particles
        end if
      end do
      ! End loop over grids
      if (ip > 0) then
        call ekin_part_helper(ind_part, ekin_part_loc, ekin_families_loc, epot_part_loc, epot_families_loc,  ip, ilevel)
      end if

#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(ekin_part_loc, energies(iekin_part), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
      call MPI_ALLREDUCE(ekin_families_loc(FAM_DM), energies(iekin_part_dm), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(ekin_families_loc(FAM_STAR), energies(iekin_part_star), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
      call MPI_ALLREDUCE(epot_part_loc, energies(iepot_part), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
      call MPI_ALLREDUCE(epot_families_loc(FAM_DM), energies(iepot_part_dm), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(epot_families_loc(FAM_STAR), energies(iepot_part_star), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#endif
#ifdef WITHOUTMPI
      energies(iekin_part) = ekin_part_loc
      energies(iekin_part_dm) = ekin_families_loc(FAM_DM)
      energies(iekin_part_star) = ekin_families_loc(FAM_STAR)
      energies(iepot_part) = epot_part_loc
      energies(iepot_part_dm) = epot_families_loc(FAM_DM)
      energies(iepot_part_star) = epot_families_loc(FAM_STAR)
#endif

    end if

  end subroutine compute_energy_part

  subroutine ekin_part_helper(ind_part, ekin_loc, ekin_families_loc, epot_loc, epot_families_loc, nn, ilevel)
    use amr_commons
    use pm_commons
    use hydro_commons
    implicit none

    integer::nn, ilevel
    integer, dimension(1:nvector)::ind_part

    integer::i, idim, nx_loc
    real(dp)::scale
    real(dp), dimension(1:nvector), save:: mass_part, phi_part
    integer, dimension(1:nvector) :: type_part
    real(dp), dimension(1:nvector, 1:ndim)::vel_part
    real(dp), dimension(-NFAMILIES:NFAMILIES), intent(inout) :: ekin_families_loc
    real(dp), dimension(-NFAMILIES:NFAMILIES), intent(inout) :: epot_families_loc
    real(dp), intent(inout):: ekin_loc, epot_loc

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

    ! Fetch mass and type
    do i = 1, nn
      mass_part(i) = mp(ind_part(i))
      type_part(i) = typep(ind_part(i))%family
#ifdef OUTPUT_PARTICLE_POTENTIAL
      phi_part(i) = ptcl_phi(ind_part(i))
      epot_loc = epot_loc + 0.5 * mass_part(i) * phi_part(i) * scale
      epot_families_loc(type_part(i)) = epot_families_loc(type_part(i)) + 0.5 * mass_part(i) *  phi_part(i) * scale
#endif
    end do

    ! Compute kinetic energy
    do idim = 1, ndim
      do i = 1, nn
        ekin_loc = ekin_loc + 0.5D0*mass_part(i)*vel_part(i, idim)**2
        ekin_families_loc(type_part(i)) = ekin_families_loc(type_part(i)) + 0.5D0*mass_part(i)*vel_part(i, idim)**2
      end do
    end do

  end subroutine ekin_part_helper

  !#########################################################################

end module deltaE_module
