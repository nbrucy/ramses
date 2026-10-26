
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
  real(dp),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_diameter =-1   ! Diameter of refinement region
  real(dp),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_x = 0          ! Center of deltaE region
  real(dp),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_y = 0          ! Center of deltaE region
  real(dp),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_z = 0          ! Center of deltaE region
  real(dp),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_exp = 2        ! Exponent for distance
  real(dp),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_a = 1          ! Ellipticity (Y/X)
  real(dp),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_b = 1          ! Ellipticity (Z/X)
  character(len=20),dimension(1:DELTAE_MAX_REGIONS)::deltaE_region_names = "" ! Names of the region

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
  real(dp), dimension(1:nb_energy_kind, 0:DELTAE_MAX_REGIONS) :: v
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

type(process), dimension(nb_process) :: deltaE

contains

  subroutine read_deltaE_params(namelist_unit,nml_ok)
   
    integer,intent(in)::namelist_unit
    logical,intent(inout)::nml_ok
    integer::nml_err
   
    namelist/deltaE_params/deltaE_correct_pressure_fix,deltaE_debug,deltaE_enable,deltaE_force_all_levels, &
    & deltaE_nregions,deltaE_region_diameter,deltaE_region_x,deltaE_region_y,deltaE_region_z, &
    & deltaE_region_exp,deltaE_region_a,deltaE_region_b

  
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

  function region_name(k) result(name)
    integer, intent(in) :: k
    character(len=20) :: name
    if (k == 0) then
      name = ""
    else
      if (deltaE_region_names(k) == "") then
        write(name, '(A,I2.2)') "region_", k
      else
        name = deltaE_region_names(k)
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
    write (*, 997) " ", (deltaE(j)%name, j = 1, nb_process)
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

  subroutine compute_energies(ilevel, energies, compute_eturb)
    implicit none

    integer, intent(in):: ilevel
    real(dp), dimension(1:nb_energy_kind, 0:deltaE_nregions), intent(inout) :: energies
    logical, intent(in) :: compute_eturb
    real(dp), dimension(1:nlevelmax) :: eturb

    energies = 0.0d0
    if (hydro) then
      call compute_energy_gas(ilevel, energies)
    end if
    if (pic) then
      call compute_energy_part(ilevel, energies(:,0))
    end if
    if (hydro .and. compute_eturb) then
      eturb = 0.0d0
      call compute_ekin_turb_all_levels(eturb)
      energies(iekin_turb0 + levelmin:iekin_turb0 + nlevelmax, 0) = eturb(levelmin:nlevelmax)
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
        energy_before(:, 0:deltaE_nregions) = energy_before(:,0:deltaE_nregions) + energy_level
      else
        energy_after = energy_after + energy_level
      end if
    end do


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

    real(dp):: dx, vol, scale
    real(dp),dimension(1:3)::skip_loc
    real(dp),dimension(1:twotondim,1:3)::xc
    real(dp),dimension(1:nvector,1:ndim),save::xx
  
    real(kind=8)::mass_loc, ekin_leaf, emag_leaf, epot_leaf
    real(kind=8), dimension(0:deltaE_nregions+1) :: ekin_loc, eint_loc, emag_loc, epot_loc
    real(dp), dimension(1:nvector, 1:nvar_all), save::uu
    real(dp), dimension(1:nvector), save::phi_level
    real(dp)::er,xr,yr,zr,rr,xn,yn,zn,r,aa,bb


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

    ! Mesh spacing at that level
    nx_loc = icoarse_max - icoarse_min + 1
    scale = boxlen/dble(nx_loc)
    if(ndim>0)skip_loc(1)=dble(icoarse_min)
    if(ndim>1)skip_loc(2)=dble(jcoarse_min)
    if(ndim>2)skip_loc(3)=dble(kcoarse_min)
    dx = 0.5D0**ilevel*scale
    vol = dx**ndim

    ! Set position of cell centers relative to grid center
    do ind=1,twotondim
      iz=(ind-1)/4
      iy=(ind-1-4*iz)/2
      ix=(ind-1-2*iy-4*iz)
      if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
      if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
      if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
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
            rr = deltaE_region_diameter  (k) ! Region diameter
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
            ekin_leaf = 0.5d0*vol*uu(i, 1 + ivar)**2/uu(i, 1)
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
                eint_loc(0) = eint_loc(0) - e_cons*vol + e_prim*vol
                eint_loc(region_index(i)) = eint_loc(region_index(i)) - e_cons*vol + e_prim*vol
            end if
          end do
        end if
      end do
      ! End loop over cells

    end do
    ! End loop over grids

    ! Compute global quantities
#ifndef WITHOUTMPI
    comm_buffin(1,:) = ekin_loc(0:deltaE_nregions)
    comm_buffin(2,:) = eint_loc(0:deltaE_nregions)
    comm_buffin(3,:) = emag_loc(0:deltaE_nregions)
    comm_buffin(4,:) = epot_loc(0:deltaE_nregions)
    call MPI_ALLREDUCE(comm_buffin, comm_buffout, 4*deltaE_nregions, MPI_DOUBLE_PRECISION, MPI_SUM,&
            &MPI_COMM_WORLD, info)
    energies(iekin_gas,:) = comm_buffout(1,:)
    energies(ieint,:) = comm_buffout(2,:)
    energies(iemag,:) = comm_buffout(3,:)
    energies(iepot_gas,:) = comm_buffout(4,:)
#endif
#ifdef WITHOUTMPI
    energies(iekin_gas,:) = ekin_loc(0:deltaE_nregions)
    energies(ieint,:) = eint_loc(0:deltaE_nregions)
    energies(iemag,:) = emag_loc(0:deltaE_nregions)
    energies(iepot_gas,:) = epot_loc(0:deltaE_nregions)
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
    real(dp), dimension(1:nb_energy_kind, 0:deltaE_nregions), intent(inout) :: energies

    real(dp), dimension(-NFAMILIES:NFAMILIES,0:deltaE_nregions+1) :: ekin_families_loc
    real(dp), dimension(-NFAMILIES:NFAMILIES,0:deltaE_nregions+1) :: epot_families_loc

    real(kind=8), dimension(0:deltaE_nregions+1)::ekin_part_loc, epot_part_loc

    ! AMR variables
    integer :: igrid, jgrid

    ! Particle variable
    integer::ipart, jpart, next_part, ip, npart1

    integer, dimension(1:nvector), save:: ind_part

    ! MPI variables
#ifndef WITHOUTMPI
    integer::info
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
      call MPI_ALLREDUCE(ekin_part_loc, energies(iekin_part,0:deltaE_nregions), 1+deltaE_nregions, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
      call MPI_ALLREDUCE(ekin_families_loc(FAM_DM,0:deltaE_nregions), energies(iekin_part_dm,0:deltaE_nregions), 1+deltaE_nregions, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(ekin_families_loc(FAM_STAR,0:deltaE_nregions), energies(iekin_part_star,0:deltaE_nregions), 1+deltaE_nregions, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
      call MPI_ALLREDUCE(epot_part_loc, energies(iepot_part,0:deltaE_nregions), 1+deltaE_nregions, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
      call MPI_ALLREDUCE(epot_families_loc(FAM_DM,0:deltaE_nregions), energies(iepot_part_dm,0:deltaE_nregions), 1+deltaE_nregions, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(epot_families_loc(FAM_STAR,0:deltaE_nregions), energies(iepot_part_star,0:deltaE_nregions), 1+deltaE_nregions, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#endif
#ifdef WITHOUTMPI
      energies(iekin_part,0:deltaE_nregions) = ekin_part_loc(0:deltaE_nregions)
      energies(iekin_part_dm,0:deltaE_nregions) = ekin_families_loc(FAM_DM,0:deltaE_nregions)
      energies(iekin_part_star,0:deltaE_nregions) = ekin_families_loc(FAM_STAR,0:deltaE_nregions)
      energies(iepot_part,0:deltaE_nregions) = epot_part_loc(0:deltaE_nregions)
      energies(iepot_part_dm,0:deltaE_nregions) = epot_families_loc(FAM_DM,0:deltaE_nregions)
      energies(iepot_part_star,0:deltaE_nregions) = epot_families_loc(FAM_STAR,0:deltaE_nregions)
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

    integer::i, idim, nx_loc, k
    real(dp)::scale, epot_part, ekin_part
    real(dp), dimension(1:nvector), save:: mass_part, phi_part
    integer, dimension(1:nvector) :: type_part
    real(dp), dimension(1:nvector, 1:ndim)::vel_part
    real(dp), dimension(-NFAMILIES:NFAMILIES,0:1+deltaE_nregions), intent(inout) :: ekin_families_loc
    real(dp), dimension(-NFAMILIES:NFAMILIES,0:1+deltaE_nregions), intent(inout) :: epot_families_loc
    real(dp), intent(inout), dimension(0:1+deltaE_nregions):: ekin_loc, epot_loc
    integer, dimension(1:nvector) :: region_index
    real(dp)::er,xr,yr,zr,rr,xn,yn,zn,r,aa,bb


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
        er = deltaE_region_exp       (k) ! Exponent defining norm
        xr = deltaE_region_x         (k) ! Region centre
        yr = deltaE_region_y         (k)
        zr = deltaE_region_z         (k)
        rr = deltaE_region_diameter  (k) ! Region diameter
        aa = deltaE_region_a         (k) ! Ellipticity (Y/X)
        bb = deltaE_region_b         (k) ! Ellipticity (Z/X)
        do i=1,nn
          xn = 0.0d0; yn=0.0d0; zn=0.0d0
          xn = abs(xp(ind_part(i), 1)-xr)
          xn = 2.0d0*xn/rr
#if NDIM > 1
          yn = abs(xp(ind_part(i),2)-yr)
          yn = 2.0d0*yn/(aa*rr)
#endif
#if NDIM >2
          zn = abs(xp(ind_part(i),3)-zr)
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



    ! Fetch mass and type
    do i = 1, nn
      mass_part(i) = mp(ind_part(i))
      type_part(i) = typep(ind_part(i))%family
#ifdef OUTPUT_PARTICLE_POTENTIAL
      phi_part(i) = ptcl_phi(ind_part(i))
      epot_part = 0.5 * mass_part(i) * phi_part(i) * scale
      epot_loc(0) = epot_loc(0) + epot_part
      epot_loc(region_index(i)) = epot_loc(region_index(i)) + epot_part
      epot_families_loc(type_part(i), 0) = epot_families_loc(type_part(i), 0) + epot_part
      epot_families_loc(type_part(i), region_index(i)) = epot_families_loc(type_part(i), region_index(i)) + epot_part
#endif
    end do

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

  end subroutine ekin_part_helper

  !#########################################################################

end module deltaE_module
