module alpha_disk_module
   use amr_parameters
   !================================================================
   ! This module contains the variable needed for disks
   !================================================================

   ! Protostellar disks
   real(dp) :: disk_radius = 1.    ! Radius of the disk
   real(dp) :: disk_density = 1.0    ! Disk density at the limit of the disk
   real(dp) :: inner_boundary = 0.25    ! Inner boundary, in unit of disk_radius
   real(dp) :: outer_boundary = 3    ! Outer boundary, in unit of disk_radius
   real(dp) :: h_over_r = 0.1 ! h/r param in lisa setup
   logical  :: damping= .true. ! Whether to use damping boudary conditions
   real(dp) :: damping_time = 1 ! Damping at inner  boundary, in units of rotation time (2 pi * omega^-1)
   real(dp) :: damping_inner_boundary = 0.1 ! Inner boundary of the damping region, in unit of disk_radius (should be less than inner_boundary)
   logical  :: disk_local_isothermal = .false. ! Whether to use local isothermal disk

contains

subroutine read_alpha_disk_params()
   implicit none

   character(LEN=80)::infile

   !--------------------------------------------------
   ! Namelist definitions
   !--------------------------------------------------
   namelist/disk_params/disk_radius, disk_density, inner_boundary, outer_boundary, h_over_r&
   & ,damping, damping_time, damping_inner_boundary &
   & ,disk_local_isothermal

   ! Read namelist file
   call getarg(1, infile) ! get the name of the namelist
   open (1, file=infile)
   rewind (1)
   read (1, NML=disk_params)
   close (1)

   if (damping .and. damping_inner_boundary > inner_boundary) then
      write(*,*) "error, damping_inner_boundary should be less than inner_boundary (in disk_radius units)"
      call clean_end()
   end if

end subroutine read_alpha_disk_params

end module alpha_disk_module


!#########################################################
!#########################################################
!#########################################################
subroutine boundary_alpha_disk(ilevel)
   Use amr_commons      !, ONLY: dp,ndim,nvector,boxlen,t
   !  use hydro_parameters !, ONLY: nvar,boundary_var,gamma,bx_bound,by_bound,bz_bound,turb,dens0,V0
   use hydro_commons
   use alpha_disk_module
   use poisson_parameters
   use gravana_utils
   implicit none
   integer::ilevel
   !----------------------------------------------------------
   ! This routine reset a part of the box to its initial value
   !----------------------------------------------------------
   integer::igrid, ngrid, ncache, i, ind, iskip, ix, iy, iz
   integer::nx_loc, idim
   real(dp)::dx, dx_loc, scale
   real(dp), dimension(1:twotondim, 1:3)::xc
   real(dp), dimension(1:3)::skip_loc

   integer, dimension(1:nvector), save::ind_grid, ind_cell
   real(dp), dimension(1:nvector, 1:ndim), save::x
   real(dp):: x0, y0, xx, yy, cs, omega, xx_softened, yy_softened, ur
   real(dp):: rc, rc_softened, damping_time_inner_boundary, damping_factor
   real(dp):: gmass, emass, r0, d0, density, ekin, eint, softened_inner_radius
   real(dp), dimension(1:nvar) :: u0 = 0.! initial condtions
   real(dp),parameter::pi = acos(-1.0d0)

   if (numbtot(1, ilevel) == 0 .or. t == 0) return

   ! Mesh size at level ilevel in coarse cell units
   dx = 0.5D0**ilevel

   ! Rescaling factors
   nx_loc = (icoarse_max - icoarse_min + 1)
   skip_loc = (/0.0d0, 0.0d0, 0.0d0/)
   if (ndim > 0) skip_loc(1) = dble(icoarse_min)
   if (ndim > 1) skip_loc(2) = dble(jcoarse_min)
   if (ndim > 2) skip_loc(3) = dble(kcoarse_min)
   scale = boxlen/dble(nx_loc)
   dx_loc = dx*scale

   ! Position of the point mass
   x0 = gravity_params(3)
   y0 = gravity_params(4)

   ! Central mass
   gmass = gravity_params(1)
   ! Softening coefficient
   emass = gravity_params(2)

   ! Density reference
   d0 = disk_density
   ! Outer limit of the disk
   r0 = disk_radius

   ! Set position of cell centers relative to grid center
   do ind = 1, twotondim
      iz = (ind - 1)/4
      iy = (ind - 1 - 4*iz)/2
      ix = (ind - 1 - 2*iy - 4*iz)
      if (ndim > 0) xc(ind, 1) = (dble(ix) - 0.5D0)*dx
      if (ndim > 1) xc(ind, 2) = (dble(iy) - 0.5D0)*dx
      if (ndim > 2) xc(ind, 3) = (dble(iz) - 0.5D0)*dx
   end do

   !---------------------------------------------------------
   ! Compute analytical velocity field for the external cells
   !---------------------------------------------------------
   ncache = active(ilevel)%ngrid

   ! Loop over grids by vector sweeps
   do igrid = 1, ncache, nvector
      ngrid = MIN(nvector, ncache - igrid + 1)
      do i = 1, ngrid
         ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
      end do

      ! Loop over cells
      do ind = 1, twotondim

         ! Gather cell indices
         iskip = ncoarse + (ind - 1)*ngridmax
         do i = 1, ngrid
            ind_cell(i) = iskip + ind_grid(i)
         end do

         ! Gather cell centre positions
         do idim = 1, ndim
            do i = 1, ngrid
               x(i, idim) = xg(ind_grid(i), idim) + xc(ind, idim)
            end do
         end do
         ! Rescale position from code units to user units
         do idim = 1, ndim
            do i = 1, ngrid
               x(i, idim) = (x(i, idim) - skip_loc(idim))*scale
            end do
         end do

         do i = 1, ngrid

            ! shift coordinate system
            xx = x(i, 1) - x0
#if NDIM > 1
            yy = x(i, 2) - y0
#endif

            ! cylindrical radius
            rc = sqrt(xx**2 + yy**2)
            if (cubic_spline_kernel) then
               rc_softened = sqrt(1d0/cubic_spline(rc, cubic_kernel_rsoft))
            else
               rc_softened = sqrt(xx**2 + yy**2 + emass**2)
            end if

            ! softened coordinates
            xx_softened = xx*(rc_softened/rc)
            yy_softened = yy*(rc_softened/rc)

            ! Reinitialize the density for the internal and external border (cylindrical)

            if (rc < r0*inner_boundary .or. rc > r0*outer_boundary) then

               density = d0*(rc_softened/r0)**(-1/2.)
               ! density
               u0(1) = density

               omega = sqrt((gmass / rc_softened**3 ) * (1 - (3/2.)*h_over_r**2))
               cs = h_over_r * sqrt(gmass / rc_softened)

               ! momentum
               u0(2) = -u0(1)*omega*yy_softened
               u0(3) = u0(1)*omega*xx_softened

               ! Also add radial velocity
               if (add_viscosity .and. alpha_viscosity > 0) then
                  ur = - (3/2.) * alpha_viscosity * cs * h_over_r
                  u0(2)  = u0(2) + u0(1)*ur*xx_softened/rc_softened
                  u0(3) = u0(3) + u0(1)*ur*yy_softened/rc_softened
               end if

               ! internal energy
               eint = u0(1)*cs**2/(gamma - 1)

               ! kinetic energy
               ekin = 0
               do idim = 2, ndim + 1
                  ekin = ekin + 0.5*u0(idim)**2/u0(1)
               end do

               ! energy
               u0(neul) = eint + ekin

               ! Apply damping
               if (damping .and. rc < r0*inner_boundary .and. rc >= r0*damping_inner_boundary  ) then
                  if (cubic_spline_kernel) then
                     softened_inner_radius = sqrt(1d0/cubic_spline(r0*inner_boundary, cubic_kernel_rsoft))
                  else
                     softened_inner_radius = sqrt((r0*inner_boundary)**2 + emass**2)
                  end if

                  damping_time_inner_boundary = damping_time * 2 * pi / sqrt(gmass / softened_inner_radius**3 ) 
                  damping_factor = max(0., 1. - ((rc/r0 - damping_inner_boundary)/(inner_boundary - damping_inner_boundary))**2)
                  uold(ind_cell(i), 1:neul) = uold(ind_cell(i), 1:neul) - (uold(ind_cell(i), 1:neul) - u0(1:neul)) &
                                                & * (damping_factor * dtold(ilevel) / damping_time_inner_boundary)
               else   ! No damping for the outer boundary
                  uold(ind_cell(i), 1:neul) = u0(1:neul)
               end if
            end if

         end do

      end do
      ! End loop over cells

   end do
   ! End loop over grids
end subroutine boundary_alpha_disk
  
!================================================================
!================================================================
!================================================================
!================================================================
  subroutine condinit_alpha_disk(x,q,dx,nn)
   use amr_parameters
   use hydro_parameters
   use poisson_parameters
   use alpha_disk_module
   use gravana_utils
 
   implicit none
   integer ::nn                            ! Number of cells
   real(dp)::dx                            ! Cell size
   real(dp),dimension(1:nvector,1:nvar_all)::q ! Primitive variables

   real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
   !================================================================
   ! This routine generates an analytical disk potential initial conditions for RAMSES.
   !================================================================
   integer ::i
   real(dp):: x0, y0, xx, yy, cs, omega, ur, xx_softened, yy_softened
   real(dp):: rc, rc_softened, gm_r3
   real(dp):: gmass, emass, r0, d0, d

   logical,save:: first_call = .true.       ! True if this is the first call to condinit

   if (first_call) then
      call read_alpha_disk_params()
      first_call = .false.
   end if
 
   ! Position of the point mass
   x0 = gravity_params(3)
   y0 = gravity_params(4)
 
   
   ! Central mass
   gmass = gravity_params(1)
   ! Softening coefficient
   emass = gravity_params(2)
 
   ! Density reference
   d0 = disk_density
 
   ! Outer limit of the disk
   r0 = disk_radius
 
   do i=1,nn
      ! shift coordinate system
      xx = x(i,1) - x0
#if NDIM > 1
      yy = x(i,2) - y0
#endif
 
      ! cylindrical radius
      rc = sqrt(xx**2 + yy**2)

      if (cubic_spline_kernel) then
         rc_softened = sqrt(1d0/cubic_spline(rc, cubic_kernel_rsoft))
      else
         rc_softened = sqrt(xx**2 + yy**2 + emass**2)
      end if
 
      ! softened coordinates
      xx_softened =  xx * (rc_softened / rc);
      yy_softened =  yy * (rc_softened / rc);
 
      ! Here d is the column density - lisa SETUP
      d = d0 * (rc_softened / r0)**(-1/2.) 
 
      d = max(d, smallr)
      q(i, 1) = d
 
      ! angular velocity
      omega = sqrt((gmass/rc_softened**3) * (1 - (3/2.)*h_over_r**2))
      cs = h_over_r * sqrt(gmass / rc_softened)
  
      ! velocity
      q(i, 2) = - omega * yy_softened
      q(i, 3) =   omega * xx_softened
 
      ! Also add radial velocity
      if (add_viscosity .and. alpha_viscosity > 0) then
         ur = - (3/2.) * alpha_viscosity * cs * h_over_r
         q(i, 2) = q(i, 2) +  ur * xx_softened / rc_softened
         q(i, 3) = q(i, 3) +  ur * yy_softened / rc_softened
      end if
 
      ! pressure
      q(i,neul) =  q(i, 1)*cs**2
      

   end do
 end subroutine condinit_alpha_disk