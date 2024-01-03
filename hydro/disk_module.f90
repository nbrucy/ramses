module disk_module
    use amr_parameters
    !================================================================
    ! This module contains the variable needed for disks
    !================================================================
  
    ! Protostellar disks
    real(dp) :: disk_radius         = 0.25   ! Radius of the disk
    real(dp) :: disk_density        = 1.0    ! Disk density at the limit of the disk
    real(dp) :: temper_iso          = 0.1    ! Sound velocity at the limit of the disk
    real(dp) :: temper_expo         = 1.0    ! Exponent used in local isothermal profile
    real(dp) :: radius_min_factor   = 0.1    ! Radius of the inner isothermal zone
    real(dp) :: contrast_factor     = 1000.  ! magnitude between the disk and the rest of the simulation
    real(dp) :: inner_iso_z_flaring = 1.0    ! Flaring of the inner isothermal zone in z to prevent resolution effects
    real(dp) :: radius_max_factor   = 3    ! Outer limit of the simulation
    !real(dp) :: gravity_threshold_factor = 1. ! Cell with a density under gravity_threshold_factor * smallr will not undergo gravity kick
    !logical  :: local_cooling       = .false. ! Disk cooling with cooling timae proportional to angular speed
    !real(dp) :: beta_cool           = 10.     ! tcool = beta_cool / omega
    logical  :: merubate           = .false.  ! whether to use Meru & Bate initial setup
  
end module disk_module
  
subroutine read_disk_params()
    use disk_module
    implicit none
  
    character(LEN=80)::infile
  
    !--------------------------------------------------
    ! Namelist definitions
    !--------------------------------------------------
    namelist/disk_params/disk_radius, disk_density, temper_iso, temper_expo, radius_min_factor, contrast_factor, inner_iso_z_flaring, radius_max_factor
  
    ! Read namelist file
    call getarg(1, infile) ! get the name of the namelist
    open (1, file=infile)
    rewind(1)
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
    integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
    integer::nx_loc,idim,neul=5
    real(dp)::dx,dx_loc,scale,u,v,w,A,B,C
    real(dp),dimension(1:twotondim,1:3)::xc
    real(dp),dimension(1:3)::skip_loc
  
    integer ,dimension(1:nvector),save::ind_grid,ind_cell
    real(dp),dimension(1:nvector,1:ndim),save::x
    real(dp):: x0, y0, z0, xx, yy, zz = 0., cs, omega, xx_soft, yy_soft, ur
    real(dp):: rc, rc_soft, rs,  rs_soft
    real(dp):: mass, emass, r0, cs0, d0, density, rin, ekin, eint
    
    if (numbtot(1, ilevel) == 0) return
  
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
  
  
    ! Set position of cell centers relative to grid center
    do ind=1,twotondim
       iz=(ind-1)/4
       iy=(ind-1-4*iz)/2
       ix=(ind-1-2*iy-4*iz)
       if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
       if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
       if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
    end do
  
    !---------------------------------------------------------
    ! Compute analytical velocity field for the external cells
    !---------------------------------------------------------
    ncache=active(ilevel)%ngrid
  
    ! Loop over grids by vector sweeps
    do igrid=1,ncache,nvector
       ngrid=MIN(nvector,ncache-igrid+1)
       do i=1,ngrid
          ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
       end do
  
       ! Loop over cells
       do ind=1,twotondim
  
          ! Gather cell indices
          iskip=ncoarse+(ind-1)*ngridmax
          do i=1,ngrid
             ind_cell(i)=iskip+ind_grid(i)
          end do
  
          ! Gather cell centre positions
          do idim=1,ndim
             do i=1,ngrid
                x(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
             end do
          end do
          ! Rescale position from code units to user units
          do idim=1,ndim
             do i=1,ngrid
                x(i,idim)=(x(i,idim)-skip_loc(idim))/scale
             end do
          end do
  
  
  
          do i=1,ngrid
  
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
             xx_soft =  xx * (rs_soft / rs)
             yy_soft =  yy * (rs_soft / rs)

             ! Inner limit of the isothermal zone
#if NDIM>2
             rin = sqrt(r0*radius_min_factor*(r0*radius_min_factor + inner_iso_z_flaring*abs(zz)))
#else
             rin = r0*radius_min_factor
#endif
     
             

             ! Reinitialize the density for the internal and external border (cylindrical)
             if (rc < r0*radius_min_factor .or. rc > r0*radius_max_factor .or. abs(zz) > 0.45 * boxlen) then
  
                ! sound velocity
                if (rc_soft > rin) then
                    cs = cs0*(rc_soft/r0)**(-temper_expo/2.)
                else
                    cs = cs0*(rin/r0)**(-temper_expo/2.)
                end if
#if NDIM>2
                if (merubate) then
                    densityd = d0 * (r0 / rc_soft)**((temper_expo - 5.)/ 2.) * exp((mass/cs**2)*(1./rs_soft - 1./rc_soft))
                else
                    density = d0 * (r0 / rc_soft)**(3 - temper_expo/2.) * exp((mass/cs**2)*(1./rs_soft - 1./rc_soft))
                end if
#else
                density  = d0 * (rc_soft / r0)**(-1/2.) 
#endif
                if(rc_soft > r0 .or. abs(zz) > 0.5 * r0) then
                    density = density / contrast_factor
                end if
                density = max(density, smallr)
  
              ! angular velocity
#if NDIM>2
                omega = sqrt(max(mass/((rc_soft**2)*rs_soft) - (4. - temper_expo/2.)*(cs**2/rc_soft**2), 0.0))
#else
                omega = sqrt(mass / rc_soft**3 - (3/2.)*(cs**2/rc_soft**2) )
#endif
                ! density
                uold(ind_cell(i), 1) = density

                     ! Also add radial velocity
               if (alpha_viscosity > 0) then
                  ur = - (3/2.) * alpha_viscosity * cs**2 * sqrt(rc_soft)
                  uold(i, 2) = uold(i, 2) + uold(i, 1) * ur * xx_soft / rc_soft
                  uold(i, 3) =  uold(i, 3) + uold(i, 1) * ur * yy_soft / rc_soft
               end if
  
                ! momentum
                uold(ind_cell(i), 2) = - uold(ind_cell(i), 1) * omega * yy_soft
                uold(ind_cell(i), 3) =  uold(ind_cell(i), 1) * omega * xx_soft
#if NDIM>2
                uold(ind_cell(i), 4) = 0.
#endif
  
                ! internal energy
                eint = + uold(ind_cell(i), 1)*cs**2 /(gamma -1)
                ! kinetic energy
                ekin = 0
                do idim=2,ndim+1
                   ekin =  ekin + 0.5*uold(ind_cell(i), idim)**2 / uold(ind_cell(i), 1)
                end do
                ! energy
                uold(ind_cell(i), ndim+2) = eint + ekin

#ifdef SOLVERMHD  
                ! magnetic field
                uold(ind_cell(i), 6:8) = (/0., 0., 0./)
                uold(ind_cell(i), nvar+1:nvar+3) =  (/0., 0., 0./)
#endif

             end if
          end do
  
       end do
       ! End loop over cells
  
    end do
    ! End loop over grids  
  
  end subroutine boundary_disk
  
  