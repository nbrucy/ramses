
module deltaE_module
    implicit none

    contains

function compute_total_energy(ilevel, use_unew) result(etot_all)
    use amr_commons
    use hydro_commons
    use mpi_mod

    implicit none

    ! Parameters & outputs
    integer, intent(in):: ilevel
    logical, intent(in):: use_unew
    real(dp):: etot_all

    ! Specific variables
    real(dp) ::  etot_loc=0

    ! AMR variables
    integer :: i,ind,ncache,igrid,iskip
    integer :: nleaf,ngrid,nx_loc
    integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf
    real(dp)::dx,vol,scale


    ! MPI variables
#ifndef WITHOUTMPI
    integer::info
#endif
    
    etot_loc = 0.0d0
    etot_all = 0.0d0

    ! Mesh spacing at that level
    nx_loc = icoarse_max-icoarse_min+1
    scale = boxlen/dble(nx_loc)
    dx = 0.5D0**ilevel*scale
    vol = dx**ndim


    ! Loop over active grids by vector sweeps
    ncache=active(ilevel)%ngrid
    do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
            ind_grid(i)=active(ilevel)%igrid(igrid + i - 1)
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
            

            if(use_unew) then 
                do i=1,nleaf
                    etot_loc = etot_loc + unew(ind_leaf(i), neul) * vol
                end do
            else 
                do i=1,nleaf
                    etot_loc = etot_loc + uold(ind_leaf(i), neul) * vol
                end do
            end if
        end do
    end do


#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(etot_loc,etot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        & MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
    etot_all = etot_loc
#endif

end function compute_total_energy


function compute_kinetic_energy_part(ilevel) result(ekin_part)
    use amr_commons
    use hydro_commons
    use mpi_mod
    use pm_commons

    implicit none

    ! Parameters & outputs
    integer, intent(in):: ilevel
    real(dp):: ekin_part


    real(kind=8)::dt_loc,ekin_loc

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
    
    ekin_part=0; ekin_loc=0

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
                        call newdt2(ind_part,dt_loc,ekin_loc,ip,ilevel)
                        ip=0
                    end if
                    end if
                    ipart=nextp(ipart)    ! Go to next particle
                end do
                ! End loop over particles
            end if
        end do
        ! End loop over grids
        if(ip>0)call newdt2(ind_part,dt_loc,ekin_loc,ip,ilevel)
        end if

#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(ekin_loc,ekin_part,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
            & MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
        ekin_part = ekin_loc
#endif

    end if

end function compute_kinetic_energy_part


end module deltaE_module