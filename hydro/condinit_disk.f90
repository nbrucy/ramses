!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit_disk(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use disk_module
  implicit none
  integer ::nn                              ! Number of cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x   ! Cell center position.
  real(dp),parameter::pi = acos(-1.0d0)     ! pi ~ 3.14

  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft,
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft,
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================

  !================================================================
  ! Initial conditions for a disk orbiting around a central mass
  ! The disk is locally isothermal
  ! WARNING : All modifications there should be repeated in velocity_fine
  ! in the file courant_init.f90
  !================================================================
  integer ::i, idim
  real(dp):: x0, y0, z0, xx, yy, cs, omega, ur, xx_soft, yy_soft
  real(dp):: rc, rc_soft
  real(dp):: mass, emass, r0, cs0, d0, d, rin, ekin, eint


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

  do i=1,nn
     ! shift coordinate system
     xx = x(i,1) - x0
     yy = x(i,2) - y0

     ! cylindrical radius
     rc = sqrt(xx**2 + yy**2)
     rc_soft = sqrt(xx**2 + yy**2 + emass**2)

     ! softened coordinates
     xx_soft =  xx * (rc_soft / rc);
     yy_soft =  yy * (rc_soft / rc);

     ! Here d is the column density - lisa SETUP
     d = d0 * (rc_soft / r0)**(-1/2.) 

     d = max(d, smallr)
     u(i, 1) = d

     ! angular velocity
     omega = sqrt((mass / rc_soft**3 ) * (1 - (3/2.)*h_over_r**2))
     cs = h_over_r * omega * rc_soft
 
     ! momentum
     u(i, 2) = - u(i, 1) * omega * yy_soft
     u(i, 3) =  u(i, 1) * omega * xx_soft

     ! Also add radial velocity
     if (alpha_viscosity > 0) then
        ur = - (3/2.) * alpha_viscosity * cs * h_over_r
        u(i, 2) = u(i, 2) + u(i, 1) * ur * xx_soft / rc_soft
        u(i, 3) = u(i, 3) + u(i, 1) * ur * yy_soft / rc_soft
     end if

     ! internal energy
     eint = + u(i, 1)*cs**2 /(gamma -1)
     ! kinetic energy
     ekin = 0
     do idim=2,ndim+1
        ekin =  ekin + 0.5*u(i, idim)**2 / u(i, 1)
     end do
     ! energy
     u(i, ndim+2) = eint + ekin

  end do
end subroutine condinit_disk