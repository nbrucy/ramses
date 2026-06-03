
module deltaE_module
    implicit none

    contains

function compute_total_energy_gas(ilevel, use_unew) result(etot_all)
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

end function compute_total_energy_gas


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
    call MPI_ALLREDUCE(ekin_loc,ekin_part,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
    ekin_part = ekin_loc
#endif

end if

end function compute_kinetic_energy_part

function compute_turbulent_energy(ilevel) result(ekin_turb)
use amr_commons
use hydro_commons
use mpi_mod

implicit none


! Parameters and outputs
integer :: ilevel
real(dp) :: ekin_turb

! AMR Variables
integer :: i, ind, ncache, igrid, iskip, idim
integer :: nleaf, ngrid, nx_loc
integer, dimension(nvector) :: ind_grid, ind_cell, ind_leaf

! Physical variables
real(dp) :: scale, dx, vol, ekin_turb_loc
real(dp), dimension(nvector, ndim):: mass_velocity_mean ! Mean velocity in the grid
real(dp), dimension(nvector):: mass, mass_velocity_variance ! Variance of the 3D velocity in the grid
real(dp) :: vel, mass_velocity ! Mass* velocity of the current leaf cell

! MPI variables
#ifndef WITHOUTMPI
    integer::info
#endif

! Initialize arrays
mass_velocity_mean = 0.0d0
mass_velocity_variance = 0.0d0
ekin_turb_loc = 0.0d0
mass = 0.0d0

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
        
        ! Compute mean velocities in grids
        do i=1,nleaf
            ! Get velocities and aggregateekin_turb
            do idim=1,ndim
                mass_velocity =  vol * uold(ind_leaf(i), 1 + idim)
                mass_velocity_mean(i, idim) = mass_velocity_mean(i, idim) + mass_velocity
                mass(i) = mass(i)  + vol * uold(ind_leaf(i), 1)
            end do 
        end do   

        ! Compute variance
        do i=1,nleaf
            do idim=1,ndim
                vel = uold(ind_leaf(i), 1 + idim) / max(uold(ind_leaf(i), 1), smallr)
                mass_velocity_variance(i) = mass_velocity_variance(i) + vol * uold(ind_leaf(i), 1)*(vel - mass_velocity_mean(i, idim)/mass(i))**2 
            end do
        end do
    end do

    ! Loo


end function

#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(ekin_turb_loc,ekin_turb,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
    ekin_turb = ekin_turb_loc
#endif


end function


! subroutine compute_energy_splitted(ilevel, ekin_all, eint_all, emag_all)
! use amr_commons
! use hydro_commons
! use poisson_commons
! use mpi_mod
! implicit none
! #ifndef WITHOUTMPI
! integer::info
! real(kind=8),dimension(4)::comm_buffin,comm_buffout
! #endif
! integer::ilevel

! integer::i,ivar,idim,ind,ncache,igrid,iskip
! integer::nleaf,ngrid,nx_loc
! integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

! real(dp)::dt_lev,dx,vol,scale
! real(kind=8)::mass_loc,ekin_loc,eint_loc,emag_loc
! real(kind=8)::mass_all,ekin_all,eint_all,emag_all
! real(dp),dimension(1:nvector,1:nvar_all),save::uu
! real(dp),dimension(1:nvector,1:ndim),save::gg


! mass_all=0.0d0; mass_loc=0.0d0
! ekin_all=0.0d0; ekin_loc=0.0d0
! emag_all=0.0d0; emag_loc=0.0d0
! eint_all=0.0d0; eint_loc=0.0d0

! if(numbtot(1,ilevel)==0)return


! ! Mesh spacing at that level
! nx_loc=icoarse_max-icoarse_min+1
! scale=boxlen/dble(nx_loc)
! dx=0.5D0**ilevel*scale
! vol=dx**ndim


! ! Loop over active grids by vector sweeps
! ncache=active(ilevel)%ngrid
! do igrid=1,ncache,nvector
!     ngrid=MIN(nvector,ncache-igrid+1)
!     do i=1,ngrid
!         ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
!     end do

!     ! Loop over cells
!     do ind=1,twotondim
!         iskip=ncoarse+(ind-1)*ngridmax
!         do i=1,ngrid
!             ind_cell(i)=ind_grid(i)+iskip
!         end do

!         ! Gather leaf cells
!         nleaf=0
!         do i=1,ngrid
!             if(son(ind_cell(i))==0)then
!             nleaf=nleaf+1
!             ind_leaf(nleaf)=ind_cell(i)
!             end if
!         end do

!         ! Gather hydro variables
!         do ivar=1,nvar_all
!             do i=1,nleaf
!             uu(i,ivar)=uold(ind_leaf(i),ivar)
!             end do
!         end do

!         ! Compute total energy
!         do i=1,nleaf
!             ekin_loc=ekin_loc+uu(i,neul)*vol
!         end do

!         ! Compute total magnetic energy
!         do ivar=1,3
!             do i=1,nleaf
!             emag_loc=emag_loc+0.125d0*(uu(i,5+ivar)+uu(i,nvar+ivar))**2*vol
!             end do
!         end do

!         ! Compute total internal energy
!         do i=1,nleaf
!             eint_loc=eint_loc+uu(i,neul)*vol
!         end do
!         do ivar=1,3
!             do i=1,nleaf
!             eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/uu(i,1)*vol &
!                     & -0.125d0*(uu(i,neul+ivar)+uu(i,nvar+ivar))**2*vol
!             end do
!         end do
! #if NENER>0
!         do ivar=1,nener
!             do i=1,nleaf
!             eint_loc=eint_loc-uu(i,nhydro+ivar)*vol
!             end do
!         end do
! #endif

!     end do
!     ! End loop over cells

! end do
! ! End loop over grids

! ! Compute global quantities
! #ifndef WITHOUTMPI
! comm_buffin(1)=ekin_loc
! comm_buffin(2)=eint_loc
! comm_buffin(3)=emag_loc
! call MPI_ALLREDUCE(comm_buffin,comm_buffout,4,MPI_DOUBLE_PRECISION,MPI_SUM,&
!         &MPI_COMM_WORLD,info)
! ekin_all=comm_buffout(1)
! eint_all=comm_buffout(2)
! emag_all=comm_buffout(3)
! #endif
! #ifdef WITHOUTMPI
! ekin_all=ekin_loc
! eint_all=eint_loc
! emag_all=emag_loc
! #endif

! end subroutine compute_energy_splitted


end module deltaE_module