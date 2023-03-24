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
  real(dp):: x0, y0, z0, xx, yy, zz=0, cs, omega, ur, xx_soft, yy_soft
  real(dp):: rc, rc_soft, rs,  rs_soft
  real(dp):: mass, emass, r0, cs0, d0, d, rin, ekin, eint


  ! Position of the point mass
  x0 = gravity_params(3)
  y0 = gravity_params(4)
  z0 = gravity_params(5)

  ! Central mass
  mass = gravity_params(1)
  ! Softening coefficient
  emass = gravity_params(2)

  ! Density reference
  d0 = disk_density
  ! Sound of speed reference
  cs0 = sqrt(temper_iso)
  ! Outer limit of the disk
  r0 = disk_radius


  do i=1,nn
     ! shift coordinate system
     xx = x(i,1) - x0
     yy = x(i,2) - y0

#if NDIM>2
     zz = x(i,3) - z0
#endif
     ! cylindrical radius
     rc = sqrt(xx**2 + yy**2)
     rc_soft = sqrt(xx**2 + yy**2 + emass**2)

     ! spherical radius
     rs = sqrt(xx**2 + yy**2 + zz**2)
     rs_soft = sqrt(xx**2 + yy**2 + zz**2 + emass**2)

     ! softened coordinates
     xx_soft =  xx * (rs_soft / rs);
     yy_soft =  yy * (rs_soft / rs);

     ! Inner limit of the isothermal zone
#if NDIM>2
     rin = sqrt(r0*radius_min_factor*(r0*radius_min_factor + inner_iso_z_flaring*abs(zz)))
#else
    rin = r0*radius_min_factor
#endif

     ! sound velocity
     if (rc_soft > rin) then
        cs = cs0*(rc_soft/r0)**(-temper_expo/2.)
     else
        cs = cs0*(rin/r0)**(-temper_expo/2.)
     end if

     ! density. The exponent of the central radial profile (3 - temper_expo/2.)
     ! is chosen in order to kill dependency of temper_expo in the mass profile.
     ! In Meru & Bate condition, the exponent is chosen so that the column density
     ! is inversely proportional to rc
#if NDIM>2
     if (merubate) then
        d = d0 * (r0 / rc_soft)**((temper_expo - 5.)/ 2.) * exp((mass/cs**2)*(1./rs_soft - 1./rc_soft))
     else
        d = d0 * (r0 / rc_soft)**(3 - temper_expo/2.) * exp((mass/cs**2)*(1./rs_soft - 1./rc_soft))
     end if
#else
     ! Here d is the column density
     d = d0 * (rc_soft / r0)**(-1/2.) 
#endif


     if(rc_soft > r0 .or. abs(zz) > 0.5 * r0) then
          d = d / contrast_factor
     end if

     d = max(d, smallr)
     u(i, 1) = d

     ! angular velocity
#if NDIM>2
     if (temper_expo == 1. .or. rc_soft > rin) then
        omega = sqrt(max(mass/((rc_soft**2)*rs_soft) - (4. - temper_expo/2.)*(cs**2/rc_soft**2), 0.0))
     else
        omega = sqrt(max(mass/(rs_soft**3) - (3. - temper_expo/2.)*(cs**2/rc_soft**2), 0.0))
     end if
#else
    omega = sqrt(mass / rc_soft**3 - (3/2.)*(cs**2/rc_soft**2) )
#endif

     ! momentum
     u(i, 2) = - u(i, 1) * omega * yy_soft
     u(i, 3) =  u(i, 1) * omega * xx_soft
#if NDIM>2
     u(i, 4) = 0.
#endif 

     ! Also add radial velocity
     if (alpha > 0) then
        ur = - (3/2.) * alpha * cs**2 * sqrt(rc_soft)
        u(i, 2) = u(i, 2) + u(i, 1) * ur * xx_soft / rc_soft
        u(i, 3) =  u(i, 3) + u(i, 1) * ur * yy_soft / rc_soft
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

#ifdef SOLVERMHD
     ! magnetic field
     u(i, 6:8) = (/0., 0., 0./)
     u(i, nvar+1:nvar+3) =  (/0., 0., 0./)
#endif

  end do
end subroutine condinit_disk