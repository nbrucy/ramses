!#########################################################
!#########################################################
!#########################################################
!#########################################################
module gravana_utils
   implicit none

   contains

   pure function cubic_spline(rr, rsoft)
      use amr_commons
      implicit none
      real(dp), intent(in)::rr, rsoft
      real(dp)::cubic_spline
      !----------------------------------------------------------------
      ! This function returns the value of the cubic spline kernel
      ! at a given distance rr and softening length rsoft
      ! The result is homogeneous to 1/rr^2
      !----------------------------------------------------------------
      if (rr < 0.5d0*rsoft) then
         cubic_spline = (32./3.*rr/rsoft**3.-192./5.*rr**3./rsoft**5.+32.*rr**4./rsoft**6.)
      else if (rr < rsoft) then
         cubic_spline = (-1./(15.*rr**2.) + 64./3.*rr/rsoft**3.-48.*rr**2./rsoft**4.+ &
                        192./5.*rr**3./rsoft**5.-160./15.*rr**4./rsoft**6.)
      else
         cubic_spline = 1/rr**2
      end if
   end function cubic_spline


   pure function grav_factor(rx, ry, rz, gmass, emass)
      use amr_commons
      use poisson_parameters
      implicit none
      real(dp), intent(in)::rx, ry, rz, gmass, emass
      real(dp)::rr, rsoft
      real(dp)::grav_factor
      !----------------------------------------------------------------
      ! This function returns the gravitational factor GM/r^3,
      ! depending on the softening strategy
      ! If cubic_spline_kernel is true, it uses the cubic spline kernel
      ! Otherwise, it uses the a classical softening length
      !----------------------------------------------------------------
      if (cubic_spline_kernel) then
         rr = sqrt(rx**2 + ry**2 + rz**2)
         rsoft = cubic_kernel_rsoft
         grav_factor = -gmass*cubic_spline(rr, rsoft)
      else
         rr = sqrt(rx**2 + ry**2 + rz**2 + emass**2)
         grav_factor = -gmass/rr**3
      end if
   end function grav_factor

end module gravana_utils

subroutine gravana(x, f, dx, ncell)
   use amr_commons
   use poisson_parameters
   use constants
   use gravana_utils
   implicit none

   integer ::ncell                         ! Size of input arrays
   real(dp)::dx                            ! Cell size
   real(dp), dimension(1:nvector, 1:ndim)::f ! Gravitational acceleration
   real(dp), dimension(1:nvector, 1:ndim)::x ! Cell center position.
   !================================================================
   ! This routine computes the acceleration using analytical models.
   ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
   ! f(i,1:ndim) is the gravitational acceleration in user units.
   ! Only if there is no self-gravity
   !================================================================
   integer::idim, i
   real(dp)::gmass, emass, xmass, ymass, zmass, rr, rx, ry, rz
   real(dp)::xmass1, ymass1, zmass1, fact, fact1, fact2
   real(dp)::gmass2, xmass2, ymass2, zmass2, emass2, omega, separation
   real(dp):: a1, a2, z0
   real(dp)::scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2

   
   select case (gravity_type)

      ! Constant vector
   case (1)
      do idim = 1, ndim
         do i = 1, ncell
            f(i, idim) = gravity_params(idim)
         end do
      end do

      ! Point mass
   case (2)
      gmass = gravity_params(1) ! GM
      emass = dx
      emass = gravity_params(2) ! Softening length
      xmass = gravity_params(3) ! Point mass coordinates
      ymass = gravity_params(4)
      zmass = gravity_params(5)

      do i = 1, ncell
         rx = 0.0d0; ry = 0.0d0; rz = 0.0d0
         rx = x(i, 1) - xmass
#if NDIM>1
         ry = x(i, 2) - ymass
#endif
#if NDIM>2
         rz = x(i, 3) - zmass
#endif
         fact = grav_factor(rx, ry, rz, gmass, emass)         

         f(i, 1) = fact*rx
#if NDIM>1
         f(i, 2) = fact*ry
#endif
#if NDIM>2
         f(i, 3) = fact*rz
#endif
      end do

   case (3)
      ! vertical galactic gravitational field
      ! Kuijken & Gilmore 1989 taken from Joung & MacLow (2006)
      ! g = -a1 z / sqrt(z^2+z0^2) - a2 z
      a1 = gravity_params(1) ! Star potential coefficient in kpc Myr-2
      a2 = gravity_params(2) ! DM potential coefficient in Myr-2
      z0 = gravity_params(3) ! Scale height in pc
      ! standard values are: a1 = 1.42d-3, a2 = 5.49d-4, z0 = 0.18d3 pc

      ! The gravitational field is given by
      ! g = -a1 z / sqrt(z^2+z0^2) - a2 z
      ! rho = [(a1 / z0) ( (z/z0)^2 + 1)^(-3/2) + a2] / (4piG)

      ! convert to code units
      call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
      a1 = a1*kpc2cm/Myr2sec**2/scale_l*scale_t**2
      a2 = a2/Myr2sec**2*scale_t**2
      z0 = z0*pc2cm/scale_l

      do i = 1, ncell
         ! the last dimension is vertical (1D -> x, 2D -> y, 3D -> z)
         rz = x(i, ndim) - 0.5d0*boxlen
         f(i, ndim) = -a1*rz/(rz**2 + z0**2)**0.5 - a2*rz
      end do

   case (4)
      ! Proper 2 body problem
      gmass = gravity_params(1) ! GM
      emass = dx
      emass = gravity_params(2) ! Softening length
      xmass = gravity_params(3) ! center of mass coordinates
      ymass = gravity_params(4)
      zmass = gravity_params(5)

      gmass2 = gmass_secondary  ! GM of the second point mass
      emass2 = soft_secondary

      omega = sqrt((gmass + gmass2)/separation**3) ! Keplerian rotation speed

      fact1 = gmass2/(gmass + gmass2)
      fact2 = gmass/(gmass + gmass2)

      xmass1 = xmass - fact1*separation*cos(omega*t)
      ymass1 = ymass - fact1*separation*sin(omega*t)
      zmass1 = zmass

      xmass2 = xmass + fact2*separation*cos(omega*t)
      ymass2 = ymass + fact2*separation*sin(omega*t)
      zmass2 = zmass

      do i = 1, ncell
         rx = 0.0d0; ry = 0.0d0; rz = 0.0d0
         rx = x(i, 1) - xmass1
#if NDIM>1
         ry = x(i, 2) - ymass1
#endif
#if NDIM>2
         rz = x(i, 3) - zmass1
#endif
         fact = grav_factor(rx, ry, rz, gmass, emass)

         f(i, 1) = -fact*rx
#if NDIM>1
         f(i, 2) = -fact*ry
#endif
#if NDIM>2
         f(i, 3) = -fact*rz
#endif

         ! redo for the second mass
         rx = 0.0d0; ry = 0.0d0; rz = 0.0d0
         rx = x(i, 1) - xmass2
#if NDIM>1
         ry = x(i, 2) - ymass2
#endif
#if NDIM>2
         rz = x(i, 3) - zmass2
#endif
         rr = sqrt(rx**2 + ry**2 + rz**2 + emass**2)
         f(i, 1) = f(i, 1) - gmass2*rx/rr**3
#if NDIM>1
         f(i, 2) = f(i, 2) - gmass2*ry/rr**3
#endif
#if NDIM>2
         f(i, 3) = f(i, 3) - gmass2*rz/rr**3
#endif
      end do

   end select
end subroutine gravana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine phi_ana(rr, pp, ngrid)
   use amr_commons
   use poisson_commons
   use constants, only: twopi
   implicit none
   integer::ngrid
   real(dp), dimension(1:nvector)::rr, pp
   ! -------------------------------------------------------------------
   ! This routine set up boundary conditions for fine levels.
   ! -------------------------------------------------------------------

   integer :: i
   real(dp):: fourpi

   fourpi = 2*twopi

  do i=1,ngrid
#if NDIM==1
     pp(i)=multipole(1)*fourpi/2*rr(i)
#elif NDIM==2
     pp(i)=multipole(1)*2*log(rr(i))
#elif NDIM==3
     pp(i)=-multipole(1)/rr(i)
#endif
  end do
end subroutine phi_ana
