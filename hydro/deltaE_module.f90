
module deltaE_module

    use amr_commons

    implicit none

    integer, parameter :: nb_energy_kind = 13
    integer, parameter :: iekin = 1, iepot = 2, ieint = 3, iemag = 4, iekin_gas = 5, iekin_part = 6, iepot_gas = 7, iepot_part = 8, &
                         & iekin_gas_turb = 9,  iekin_part_star = 10, iekin_part_dm = 11, iepot_part_star = 12, iepot_part_dm = 13
    character(len=14),  dimension(1:nb_energy_kind) :: energy_names = [character(len=14) ::"ekin", "epot", "eint", "emag", "ekin_gas", "ekin_part", "epot_gas", "epot_part", &
                                & "ekin_gas_turb", "ekin_part_star", "ekin_part_dm", "epot_part_star", "epot_part_dm"]
   
    type processes
        real(dp), dimension(1:nb_energy_kind) :: flux_gas, flux_part, cooling, gravity_gas, gravity_part, &
                                                &  star_formation, feedback, magnetic_diffusion, corrections
    contains
        procedure :: initialize_processes
        procedure :: print_processes
    end type


    type(processes) :: deltaE 


    contains


subroutine initialize_processes(this)

    implicit none
    class(processes), intent(inout) :: this

    this%flux_gas = 0.0d0
    this%flux_part = 0.0d0
    this%cooling = 0.0d0
    this%gravity_gas = 0.0d0
    this%gravity_part = 0.0d0
    this%star_formation = 0.0d0
    this%feedback = 0.0d0
    this%magnetic_diffusion = 0.0d0
    this%corrections = 0d0

end subroutine


subroutine print_energies(energies)
    implicit none
    real(dp), dimension(1:nb_energy_kind), intent(in) :: energies
    integer :: i

    write(*,*) "Energy output"
    do i=1, nb_energy_kind 
        write(*,998) " ", energy_names(i), energies(i)
    end do
    write(*,*) "End energy output"
    
    998 format(A1,A14, 1pe20.7)


end subroutine

subroutine print_processes(this)

    implicit none
    class(processes), intent(in) :: this
    integer :: i


    write(*,*) "DeltaE output"
    write(*,997) " ", "flux_gas", "flux_part", "cooling", "gravity_gas", "gravity_part", "star_formation", "feedback", "magnetic_diffusion", "corrections"
    do i=1, nb_energy_kind 
        write(*,998) energy_names(i), this%flux_gas(i), this%flux_part(i),  this%cooling(i), this%gravity_gas(i), &
                     this%gravity_part(i), this%star_formation(i), this%feedback(i), this%magnetic_diffusion(i), &
                     & this%corrections(i)
    end do
    write(*,*) "End deltaE output"
    
    997 format(A20, A20, A20, A20, A20, A20, A20, A20, A20, A20)
    998 format(A20,1pe20.5, 1pe20.5,1pe20.5,1pe20.5,1pe20.5,1pe20.5,1pe20.5,1pe20.5,1pe20.5)
    
end subroutine

subroutine compute_energies(ilevel, use_unew, energies)

    implicit none

    integer, intent(in):: ilevel
    logical, intent(in):: use_unew  
    real(dp), dimension(1:nb_energy_kind), intent(inout) :: energies

    energies = 0.0d0
    call compute_energy_gas(ilevel, energies(iekin_gas), energies(ieint), energies(iemag), use_unew)
    call compute_ekin_part(ilevel, energies(iekin_part), energies(iekin_part_dm), energies(iekin_part_star))
    energies(iekin) = energies(iekin_gas) + energies(iekin_part)

    energies(iekin_gas_turb) = compute_ekin_turb(ilevel)
end subroutine

subroutine compute_transfer(levelstart, levelend, use_unew, deltaE_process, step)

    implicit none 
    integer, intent(in):: levelstart, levelend
    logical, intent(in):: use_unew
    real(dp), dimension(1:nb_energy_kind), intent(inout) :: deltaE_process
    integer, intent(in) :: step 

    integer :: ilevel

    real(dp), dimension(1:nb_energy_kind), save :: energy_before
    real(dp), dimension(1:nb_energy_kind) ::energy_level, energy_after

    if (step == 1) then 
        energy_before = 0.0d0
    else
        energy_after = 0.0d0
    end if

    do ilevel=levelstart, levelend
        call compute_energies(ilevel, use_unew, energy_level)
        if (step == 1) then 
            energy_before = energy_before + energy_level
        else
            energy_after = energy_after + energy_level
        end if
    end do

    if (step == 2) then 
        deltaE_process = deltaE_process + energy_after - energy_before
    end if

end subroutine



subroutine compute_ekin_part(ilevel, ekin_part, ekin_part_dm, ekin_part_star)
    use amr_commons
    use hydro_commons
    use mpi_mod
    use pm_commons

    implicit none

    ! Parameters & outputs
    integer, intent(in):: ilevel
    real(dp), intent(inout):: ekin_part, ekin_part_dm, ekin_part_star
    real(dp), dimension(-NFAMILIES:NFAMILIES) :: ekin_families_loc


    real(kind=8)::ekin_loc

    ! AMR variables
    integer :: igrid, jgrid

    ! Particle variable
    integer :: ip, npart1, ipart, jpart
    integer,dimension(1:nvector),save::ind_part
    logical :: ok


    ! MPI variables
#ifndef WITHOUTMPI
    integer::info
#endif

    ekin_part=0; ekin_part_dm=0; ekin_part_star=0
    ekin_loc=0
    ekin_families_loc = 0

    if (pic) then 

        ! Compute maximum time step on active region
        if(numbl(myid,ilevel)>0)then
        ! Loop over grids
        ip=0
        do jgrid=1,active(ilevel)%ngrid
            igrid=active(ilevel)%igrid(jgrid)
            npart1=numbp(igrid)   ! Number of particles in the grid
            if(npart1>0)then
                ! Loop over particles
                ipart=headp(igrid)
                do jpart=1,npart1
                    ! MC_tracer ================
                    ! skip tracer particles here
                    if (tracer) then
                    ok = is_not_tracer(typep(ipart))
                    else
                    ok = .true.
                    end if
                    ! End MC Tracer ============

                    if (ok) then
                    ip=ip+1
                    ind_part(ip)=ipart
                    if(ip==nvector)then
                        call ekin_part_helper(ind_part,ekin_loc,ekin_families_loc,ip,ilevel)
                        ip=0
                    end if
                    end if
                    ipart=nextp(ipart)    ! Go to next particle
                end do
                ! End loop over particles
            end if
        end do
        ! End loop over grids
        if(ip>0)call ekin_part_helper(ind_part,ekin_loc,ekin_families_loc,ip,ilevel)
        end if

#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(ekin_loc,ekin_part,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(ekin_families_loc(FAM_DM),ekin_part_dm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(ekin_families_loc(FAM_STAR),ekin_part_star,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
        ekin_part = ekin_loc
        ekin_part_dm = ekin_families_loc(FAM_DM)
        ekin_part_star = ekin_families_loc(FAM_STAR)
#endif

    end if

end subroutine compute_ekin_part

subroutine ekin_part_helper(ind_part,ekin_loc,ekin_families_loc,nn,ilevel)
    use amr_commons
    use pm_commons
    use hydro_commons
    implicit none
    real(kind=8)::ekin_loc

    integer::nn,ilevel
    integer,dimension(1:nvector)::ind_part
  
    integer::i,idim,nx_loc
    real(dp)::scale
    real(dp),dimension(1:nvector),save:: mass_part, type_part
    real(dp),dimension(1:nvector,1:ndim)::vel_part
    real(dp), dimension(-NFAMILIES:NFAMILIES) :: ekin_families_loc
    real(dp)::dx

    ! Compute time step
    dx=0.5D0**ilevel
    nx_loc=(icoarse_max-icoarse_min+1)
    scale=boxlen/dble(nx_loc)
  
    do idim=1,ndim
       do i=1,nn
          vel_part(i, idim) = vp(ind_part(i), idim)
       end do
    end do
  
    ! Fetch mass and type
    do i = 1, nn
       mass_part(i) = mp(ind_part(i))
       type_part(i) = typep(ind_part(i))%family
    end do
  
    ! Compute kinetic energy
    do idim=1,ndim
       do i=1,nn
          ekin_loc=ekin_loc+0.5D0*mass_part(i)*vel_part(i, idim)**2
          ekin_families_loc(type_part(i)) = ekin_families_loc(type_part(i)) + 0.5D0*mass_part(i)*vel_part(i, idim)**2
       end do
    end do
  
  end subroutine ekin_part_helper
  



function compute_ekin_turb(ilevel) result(ekin_turb)
    use amr_commons
    use hydro_commons
    use mpi_mod
    use amr_constants, only:i1min,i1max,j1min,j1max,k1min,k1max

    implicit none

    ! Parameters and outputs
    integer, intent(in) :: ilevel
    real(dp) :: ekin_turb

    ! AMR Variables
    integer :: i, ind, ncache, igrid, iskip, idim
    integer :: nleaf, ngrid,  nx_loc
    integer, dimension(nvector) :: ind_grid, ind_cell, ind_leaf, index_current_grid
    integer,dimension(1:nvector,1:threetondim)::nbor_cells
    integer::i1,j1,k1,ind_father


    ! Physical variables
    real(dp) :: scale, dx, vol, ekin_turb_loc
    real(dp), dimension(nvector, ndim):: velocity_mean ! Mean velocity in the grid
    real(dp), dimension(nvector):: mass
    real(dp), dimension(ndim) :: mass_velocity ! Mass* velocity of the current leaf cell
    real(dp) :: vel

    ! MPI variables
#ifndef WITHOUTMPI
        integer::info
#endif
    
    ! Initialize arrays
    ekin_turb_loc = 0.0d0
    
    ! Mesh spacing at that level
    nx_loc = icoarse_max-icoarse_min+1
    scale = boxlen/dble(nx_loc)
    dx = 0.5D0**ilevel*scale
    vol = dx**ndim
    
    
    ! Loop over active grids by vector sweeps
    ncache=active(ilevel)%ngrid
    do igrid=1,ncache,nvector

        velocity_mean = 0.0d0
        mass = 0.0d0

        ngrid=MIN(nvector,ncache-igrid+1)

        do i=1,ngrid
            ind_grid(i)=active(ilevel)%igrid(igrid + i - 1)
            ind_cell(i)=father(ind_grid(i)) ! also gather father cells
        end do

        ! Collect neighbor cells 
        call get3cubefather(ind_cell,nbor_cells,ngrid,ilevel-1)

        ! Loop over 3x3x3 neighboring father cells 
        do k1=k1min,k1max
            do j1=j1min,j1max
                do i1=i1min,i1max

                    ! Get neighbor cell index
                    ind_father=1+i1+3*j1+9*k1

                    do i=1,ngrid
                        ! Get mass and velocity and aggregate
                        mass_velocity = 2**(ndim)*vol*uold(nbor_cells(i,ind_father), 2:1+ndim)
                        velocity_mean(i, :) = velocity_mean(i, :) + mass_velocity
                        mass(i) = mass(i) + 2**(ndim)*vol * uold(nbor_cells(i,ind_father), 1)
                    end do
                end do
            end do
        end do
        
        ! Normalize
        do i=1,ngrid
            velocity_mean(i, :) = velocity_mean(i, :) /  mass(i) 
        end do
    
           ! Loop over cells
        do ind=1,twotondim
            iskip=ncoarse+(ind-1)*ngridmax
            do i=1,ngrid
                ind_cell(i)=ind_grid(i)+iskip
            end do

            ! Gather leaf cells
            nleaf=0
            do i=1,ngrid
                if(son(ind_cell(i))==0)then
                    nleaf=nleaf+1
                    ind_leaf(nleaf)=ind_cell(i)
                    index_current_grid(nleaf) = i
                end if
            end do
    
            ! Compute variance
            do i=1,nleaf
                do idim=1,ndim
                    vel = uold(ind_leaf(i), 1 + idim) / max(uold(ind_leaf(i), 1), smallr)
                    ekin_turb_loc = ekin_turb_loc + 0.5* (vol * uold(ind_leaf(i), 1)*(vel - velocity_mean(index_current_grid(i),idim))**2)
                end do
            end do
        end do
    end do
    
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(ekin_turb_loc,ekin_turb,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
        ekin_turb = ekin_turb_loc
#endif

end function


subroutine compute_energy_gas(ilevel, ekin_all, eint_all, emag_all, use_unew)
    use amr_commons
    use hydro_commons
    use poisson_commons
    use mpi_mod
    implicit none
#ifndef WITHOUTMPI
    integer::info
    real(kind=8),dimension(4)::comm_buffin,comm_buffout
#endif
    integer, intent(in)::ilevel
    real(dp), intent(inout)::ekin_all,eint_all,emag_all
    logical, intent(in)::use_unew


    integer::i,ivar,idim,ind,ncache,igrid,iskip
    integer::nleaf,ngrid,nx_loc
    integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

    real(dp)::dt_lev,dx,vol,scale
    real(kind=8)::mass_loc,ekin_loc,eint_loc,emag_loc, ekin_leaf, emag_leaf
    real(kind=8)::mass_all
    real(dp),dimension(1:nvector,1:nvar_all),save::uu
    real(dp),dimension(1:nvector,1:ndim),save::gg


    mass_all=0.0d0; mass_loc=0.0d0
    ekin_all=0.0d0; ekin_loc=0.0d0
    emag_all=0.0d0; emag_loc=0.0d0
    eint_all=0.0d0; eint_loc=0.0d0
    ekin_leaf=0.0d0; emag_leaf=0.0d0


    if(numbtot(1,ilevel)==0)return


    ! Mesh spacing at that level
    nx_loc=icoarse_max-icoarse_min+1
    scale=boxlen/dble(nx_loc)
    dx=0.5D0**ilevel*scale
    vol=dx**ndim


    ! Loop over active grids by vector sweeps
    ncache=active(ilevel)%ngrid
    do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
            ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim
            iskip=ncoarse+(ind-1)*ngridmax
            do i=1,ngrid
                ind_cell(i)=ind_grid(i)+iskip
            end do

            ! Gather leaf cells
            nleaf=0
            do i=1,ngrid
                if(son(ind_cell(i))==0)then
                    nleaf=nleaf+1
                    ind_leaf(nleaf)=ind_cell(i)
                end if
            end do

            ! Gather hydro variables
            do ivar=1,nvar_all
                if (use_unew) then
                    do i=1,nleaf
                        uu(i,ivar)=unew(ind_leaf(i),ivar)
                    end do
                else 
                    do i=1,nleaf
                        uu(i,ivar)=uold(ind_leaf(i),ivar)
                    end do
                end if
            end do

            ! Compute total internal energy step 1
            do i=1,nleaf
                eint_loc=eint_loc + uu(i,neul)*vol
            end do

            ! Compute total energies
            do ivar=1,ndim
                do i=1,nleaf
                ! Compute total kinetic energy
                    ekin_leaf = 0.5d0*vol*uu(i, 1 + ivar)**2 / uu(i,1)        
                    ekin_loc = ekin_loc + ekin_leaf
#ifdef SOLVERmhd
                ! Compute total magnetic energy
                    emag_leaf = 0.125d0*(uu(i,neul+ivar)+uu(i,nvar+ivar))**2*vol
                    emag_loc = emag_loc + emag_leaf
#endif

                ! Compute total internal energy step 2
                    eint_loc = eint_loc - ekin_leaf - emag_leaf
                end do
            end do

            ! Compute total internal energy step 3
#if NENER>0
            do ivar=1,nener
                do i=1,nleaf
                    eint_loc=eint_loc-uu(i,nhydro+ivar)*vol
                end do
            end do
#endif

        end do
        ! End loop over cells

    end do
    ! End loop over grids

    ! Compute global quantities
#ifndef WITHOUTMPI
    comm_buffin(1)=ekin_loc
    comm_buffin(2)=eint_loc
    comm_buffin(3)=emag_loc
    call MPI_ALLREDUCE(comm_buffin,comm_buffout,4,MPI_DOUBLE_PRECISION,MPI_SUM,&
            &MPI_COMM_WORLD,info)
    ekin_all=comm_buffout(1)
    eint_all=comm_buffout(2)
    emag_all=comm_buffout(3)
#endif
#ifdef WITHOUTMPI
    ekin_all=ekin_loc
    eint_all=eint_loc
    emag_all=emag_loc
#endif

end subroutine compute_energy_gas


end module deltaE_module