module disk_module
   use amr_parameters
   !================================================================
   ! This module contains the variable needed for disks
   !================================================================

   ! Protostellar disks
   real(dp) :: disk_radius = 0.25    ! Radius of the disk
   real(dp) :: disk_density = 1.0    ! Disk density at the limit of the disk
   real(dp) :: radius_min_factor = 0.5    ! Radius of the inner isothermal zone
   real(dp) :: radius_max_factor = 3    ! Outer limit of the simulation
   real(dp) :: h_over_r=0.1 ! h/r param in lisa setup
   logical  :: damping= .true. ! Whether to use damping boudary conditions
   real(dp) :: damping_time=1 ! Damping at inner and outer boundary, in units of rotation time (2 pi * omega^-1)


end module disk_module

subroutine read_disk_params()
   use disk_module
   implicit none

   character(LEN=80)::infile

   !--------------------------------------------------
   ! Namelist definitions
   !--------------------------------------------------
   namelist/disk_params/disk_radius, disk_density, radius_min_factor, radius_max_factor, h_over_r, damping, damping_time

   ! Read namelist file
   call getarg(1, infile) ! get the name of the namelist
   open (1, file=infile)
   rewind (1)
   read (1, NML=disk_params)
   close (1)

end subroutine read_disk_params

!#########################################################
!#########################################################
!#########################################################
subroutine boundary_disk(ilevel)
   Use amr_commons      !, ONLY: dp,ndim,nvector,boxlen,t
   !  use hydro_parameters !, ONLY: nvar,boundary_var,gamma,bx_bound,by_bound,bz_bound,turb,dens0,V0
   use hydro_commons
   use disk_module
   use poisson_parameters
   implicit none
   integer::ilevel
   !----------------------------------------------------------
   ! This routine reset a part of the box to its initial value
   !----------------------------------------------------------
   integer::igrid, ngrid, ncache, i, ind, iskip, ix, iy, iz
   integer::nx_loc, idim, neul = 5
   real(dp)::dx, dx_loc, scale, u, v, w, A, B, C
   real(dp), dimension(1:twotondim, 1:3)::xc
   real(dp), dimension(1:3)::skip_loc

   integer, dimension(1:nvector), save::ind_grid, ind_cell
   real(dp), dimension(1:nvector, 1:ndim), save::x
   real(dp):: x0, y0, z0, xx, yy, cs, omega, xx_soft, yy_soft, ur
   real(dp):: rc, rc_soft, damping_time_r
   real(dp):: mass, emass, r0, d0, density, rin, ekin, eint
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
   scale = dble(nx_loc)/boxlen
   dx_loc = dx/scale

   ! Position of the point mass
   x0 = gravity_params(3)
   y0 = gravity_params(4)

   ! Central mass
   mass = gravity_params(1)
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
               x(i, idim) = (x(i, idim) - skip_loc(idim))/scale
            end do
         end do

         do i = 1, ngrid

            ! shift coordinate system
            xx = x(i, 1) - x0
            yy = x(i, 2) - y0

            ! cylindrical radius
            rc = sqrt(xx**2 + yy**2)
            rc_soft = sqrt(xx**2 + yy**2 + emass**2)

            ! softened coordinates
            xx_soft = xx*(rc_soft/rc)
            yy_soft = yy*(rc_soft/rc)

           


            ! Reinitialize the density for the internal and external border (cylindrical)

            if (rc < r0*radius_min_factor .or. rc > r0*radius_max_factor) then


               density = d0*(rc_soft/r0)**(-1/2.)
               ! density
               u0(1) = density

               omega = sqrt((mass / rc_soft**3 ) * (1 - (3/2.)*h_over_r**2))
               cs = h_over_r * sqrt(mass / rc_soft)

               ! momentum
               u0(2) = -u0(1)*omega*yy_soft
               u0(3) = u0(1)*omega*xx_soft

               ! Also add radial velocity
               if (alpha_viscosity > 0) then
                  ur = - (3/2.) * alpha_viscosity * cs * h_over_r
                  u0(2)  = u0(2) + u0(1)*ur*xx_soft/rc_soft
                  u0(3) = u0(3) + u0(1)*ur*yy_soft/rc_soft
               end if

               ! internal energy
               eint = u0(1)*cs**2/(gamma - 1)

               ! kinetic energy
               ekin = 0
               do idim = 2, ndim + 1
                  ekin = ekin + 0.5*u0(idim)**2/u0(1)
               end do

               ! energy
               u0(ndim + 2) = eint + ekin

               ! Apply damping
               if (damping .and. rc < r0*radius_min_factor) then
                  damping_time_r = damping_time * 2 * pi / omega
                  uold(ind_cell(i), 1:ndim + 2) = uold(ind_cell(i), 1:ndim + 2) - (uold(ind_cell(i), 1:ndim + 2) - u0(1:ndim + 2)) * (dtold(ilevel) / damping_time_r)
               else   ! No damping for the outer boundary
                  uold(ind_cell(i), 1:ndim + 2) = u0(1:ndim + 2)
               end if
            end if

         end do

      end do
      ! End loop over cells

   end do
   ! End loop over grids

end subroutine boundary_disk

