SUBROUTINE vaccafit(M,S)
  ! This routine is called in sink_RT_feedback
  ! perform a fit of the Vacca et al. 96 ionising flux
  ! M - stellar mass / solar masses
  ! S - photon emission rate in / s

  use amr_parameters,only:dp
  use sink_feedback_parameters
  implicit none

  real(dp),intent(in)::M
  real(dp),intent(out)::S
  
  S = stf_K * (M / stf_m0)**stf_a / (1. + (M / stf_m0)**stf_b)**stf_c

END SUBROUTINE
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_sn_stellar
  use pm_commons
  use amr_commons
  use hydro_commons
  use sink_feedback_parameters
  use constants, only:pi,pc2cm
  use mpi_mod
  implicit none

  integer:: ivar
  integer:: ilevel, ind, ix, iy, iz, ngrid, iskip, idim
  integer:: i, nx_loc, igrid, ncache
  integer, dimension(1:nvector), save:: ind_grid, ind_cell
  real(dp):: dx, scale, dx_min, dx_loc, vol_loc
  real(dp), dimension(1:3):: skip_loc
  real(dp), dimension(1:twotondim, 1:3):: xc
  logical, dimension(1:nvector), save:: ok
  real(dp), dimension(1:nvector, 1:ndim), save:: xx
  real(dp):: sn_r, sn_m, sn_p_local, sn_e_local, sn_vol, sn_d, sn_ed
  real(dp):: rr,pgas,dgas,ekin,mass_sn_tot,mass_sn_tot_all
  integer:: info
  real(dp),dimension(1:nvector,1:ndim)::x
  real(dp),dimension(1:3):: xshift, x_sn
  logical, save:: first = .true.
  real(dp), save:: xseed
  real(dp) ::dens_max_loc,dens_max_loc_all
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  logical, dimension(1:nstellarmax):: mark_del
  integer:: istellar,isink
  real(dp)::T_sn,sn_ed_lim,pnorm_sn,pnorm_sn_all,vol_sn,vol_sn_all,vol_rap, mass_sn, mass_sn_all,dens_moy
  real(dp)::pgas_check,pgas_check_all
  integer, parameter:: navg = 3
  integer, parameter:: nsph = 1
  real(dp), dimension(1:nsph, 1:3):: avg_center
  real(dp), dimension(1:nsph):: avg_radius
  real(dp), dimension(1:navg):: avg_rpow
  real(dp), dimension(1:navg, 1:nvar+3):: avg_upow
  real(dp), dimension(1:navg, 1:nsph):: avg
  real(dp):: norm, rad_sn

  if(.not. hydro)return
  if(ndim .ne. 3)return

  if(verbose)write(*,*)'Entering make_sn_stellar'

  if (first) then
     xseed = 0.5
     call random_number(xseed)
     first = .false.
  endif

  ! Mesh spacing in that level
  nx_loc = icoarse_max - icoarse_min + 1
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale = boxlen / dble(nx_loc)
  dx_min = scale * 0.5d0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  sn_r = 3.0d0*(0.5d0**levelmin)*scale
  if(sn_r_sat .ne. 0) sn_r = max(sn_r, sn_r_sat * pc2cm / scale_l) !impose a minimum size of 12 pc for the radius
!  sn_r = 2.*(0.5**levelmin)*scale
  sn_m = sn_mass_ref !note this is replaced later
  sn_p_local = sn_p_ref
  sn_e_local = sn_e_ref

!  if(sn_r /= 0.0) then
    sn_vol = 4. / 3. * pi * sn_r**3
!    sn_d = sn_m / sn_vol
!    sn_ed = sn_e_local / sn_vol
!  end if

  !we loop over stellar objects to determine whether one is turning supernovae
  !after it happens, the object is removed from the list
  mark_del = .false.
  do istellar = 1, nstellar
    if(t - tstellar(istellar) < ltstellar(istellar)) cycle
    mark_del(istellar) = .true.

    ! find correct index in sink array which will be equal or lower than id_sink due to sink merging
    isink = id_stellar(istellar)
    do while (id_stellar(istellar) .ne. idsink(isink))
      isink = isink - 1
    end do

    ! the mass of the massive stars 
    sn_m = mstellar(istellar) 

    !remove the mass that is dumped in the grid
    msink(isink) = msink(isink) - sn_m

    !the velocity dispersion times the life time of the object
    rad_sn = ltstellar(istellar)*Vdisp

    !find a random point within a sphere of radius < 1
    norm=2.
    do while (norm .gt. 1)
       norm=0.
       do idim = 1, ndim
          call random_number(xseed)
          xshift(idim) = (xseed-0.5)*2.
          norm = norm + xshift(idim)**2
       end do
    end do
    do idim = 1, ndim - 1 
        xshift(idim) = xshift(idim) * rad_sn
    end do
    !special treatment for the z-coordinates to maintain it at low altitude
    xshift(3) = xshift(3) * min(rad_sn,100.)

    ! place the supernovae around sink particles
    x_sn(:) = xsink(isink, :) + xshift(:)

    !apply periodic boundary conditions (only along x and y)
    !PH note that this should also be modified for the shearing box 24/01/2017
    ! TC: TODO: check which BC!
    ! if BC -> aply, else: move SN back in the box
    if( x_sn(1) .lt. 0) x_sn(1) = boxlen - x_sn(1) 
    if( x_sn(2) .lt. 0) x_sn(2) = boxlen - x_sn(2) 
    if( x_sn(1) .gt. boxlen) x_sn(1) = - boxlen + x_sn(1) 
    if( x_sn(2) .gt. boxlen) x_sn(2) = - boxlen + x_sn(2) 

    if(.true.) then
      avg_center(1, :) = x_sn(:)
      avg_radius = sn_r
      avg_rpow = 0.0d0
      avg_upow = 0.0d0
      ! avg_rpow(1) = 0 ; avg_upow(1, :) = 0 -> integrand(1) = 1
      ! avg_rpow(2) = 0 ; avg_upow(2, 1) = 1 -> integrand(2) = density
      avg_upow(2, 1) = 1.0d0
      ! avg_rpow(3) = 1 ; avg_upow(3, :) = 0 -> integrand(3) = radius
      avg_rpow(3) = 1.0d0
      call sphere_average(navg, nsph, avg_center, avg_radius, avg_rpow, avg_upow, avg)
      vol_sn = avg(1, 1)
      mass_sn = avg(2, 1) + vol_sn * sn_d ! region average + ejecta
      pnorm_sn = avg(3, 1)
    else
      !do a first path to compute the volume of the cells that are enclosed in the supernovae radius
      !this is to correct for the grid effects
      vol_sn = 0. ; vol_sn_all=0.
      mass_sn = 0. ; mass_sn_all=0.
      pnorm_sn = 0. ; pnorm_sn_all=0.

      do ilevel = levelmin, nlevelmax
        ! Computing local volume (important for averaging hydro quantities)
        dx = 0.5d0**ilevel
        dx_loc = dx * scale
        vol_loc = dx_loc**ndim

        ! Cell center position relative to grid center position
        do ind=1,twotondim
          iz = (ind - 1) / 4
          iy = (ind - 1 - 4 * iz) / 2
          ix = (ind - 1 - 2 * iy - 4 * iz)
          if(ndim>0) xc(ind,1) = (dble(ix) - 0.5d0) * dx
          if(ndim>1) xc(ind,2) = (dble(iy) - 0.5d0) * dx
          if(ndim>2) xc(ind,3) = (dble(iz) - 0.5d0) * dx
        end do

        ! Loop over grids
        ncache=active(ilevel)%ngrid
        do igrid = 1, ncache, nvector
          ngrid = min(nvector, ncache - igrid + 1)
          do i = 1, ngrid
            ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
          end do

          ! Loop over cells
          do ind = 1, twotondim
            ! Gather cell indices
            iskip = ncoarse + (ind - 1) * ngridmax
            do i = 1, ngrid
              ind_cell(i) = iskip + ind_grid(i)
            end do

            ! Gather cell center positions
            do i = 1, ngrid
              xx(i, :) = xg(ind_grid(i), :) + xc(ind, :)
            end do
            ! Rescale position from coarse grid units to code units
            do idim=1,ndim
               do i=1,ngrid
                  xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
               end do
            end do

            ! Flag leaf cells
            do i = 1, ngrid
              ok(i) = (son(ind_cell(i)) == 0)
            end do

            do i = 1, ngrid
              if(ok(i)) then
                rr = 0.
                do idim=1,ndim
                   rr = rr + ( (xx(i,idim) - x_sn(idim)) / sn_r)**2
                enddo

                if(rr < 1.) then
                  vol_sn = vol_sn + vol_loc
                  mass_sn = mass_sn + vol_loc * (uold(ind_cell(i), 1) + sn_d)
                  pnorm_sn = pnorm_sn + vol_loc * sqrt(rr)
                endif
              endif
            end do
            !  End loop over sublist of cells
          end do
          ! End loop over cells
        end do
        ! End loop over grids
      end do
      ! End loop over levels

#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(vol_sn,vol_sn_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(mass_sn,mass_sn_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(pnorm_sn,pnorm_sn_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)

      vol_sn  = vol_sn_all
      mass_sn = mass_sn_all
      pnorm_sn = pnorm_sn_all
#endif
    end if

    !compute energy and mass density
    sn_d = sn_m / vol_sn
    sn_ed = sn_e_local / vol_sn

    dens_moy = mass_sn / vol_sn

    dens_max_loc = 0.
    mass_sn_tot = 0.
    dens_max_loc_all = 0.
    mass_sn_tot_all = 0.
    pgas_check=0.

    !now loop over cells again and damp energies, mass and momentum
    !loop over levels 
    do ilevel = levelmin, nlevelmax
      ! Computing local volume (important for averaging hydro quantities)
      dx = 0.5d0**ilevel
      dx_loc = dx * scale
      vol_loc = dx_loc**ndim

      ! Cell center position relative to grid center position
      do ind=1,twotondim
        iz = (ind - 1) / 4
        iy = (ind - 1 - 4 * iz) / 2
        ix = (ind - 1 - 2 * iy - 4 * iz)
        xc(ind,1) = (dble(ix) - 0.5d0) * dx
        xc(ind,2) = (dble(iy) - 0.5d0) * dx
        xc(ind,3) = (dble(iz) - 0.5d0) * dx
      end do

      ! Loop over grids
      ncache=active(ilevel)%ngrid
      do igrid = 1, ncache, nvector
        ngrid = min(nvector, ncache - igrid + 1)
        do i = 1, ngrid
          ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
        end do

        ! Loop over cells
        do ind = 1, twotondim
          ! Gather cell indices
          iskip = ncoarse + (ind - 1) * ngridmax
          do i = 1, ngrid
            ind_cell(i) = iskip + ind_grid(i)
          end do

          ! Gather cell center positions
          do i = 1, ngrid
            xx(i, :) = xg(ind_grid(i), :) + xc(ind, :)
          end do
          ! Rescale position from coarse grid units to code units
          do idim=1,ndim
             do i=1,ngrid
                xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
             end do
          end do

          ! Flag leaf cells
          do i = 1, ngrid
            ok(i) = (son(ind_cell(i)) == 0)
          end do

          do i = 1, ngrid
            if(ok(i)) then
!                rr = sqrt(sum(((xx(i,:) - x_sn(:)) / sn_r)**2))
                rr = 0.
                do idim=1,ndim
                   rr = rr + ((xx(i,idim) - x_sn(idim)) / sn_r)**2
                enddo
                rr = sqrt(rr)

                if(rr < 1.) then
                  uold(ind_cell(i), 1) = uold(ind_cell(i), 1) + sn_d
                  dgas = uold(ind_cell(i), 1)

                  if(dgas .gt. dens_max_loc) dens_max_loc = dgas

                  mass_sn_tot = mass_sn_tot + dgas*vol_loc


                  !compute velocity of the gas within this cell assuming 
                  !energy equipartition
                  !pgas = sqrt(eff_sn*sn_ed / max(dens_moy,dgas) ) * dgas 
                  !for cells where dgas < dens_moy, take the same velocity as if
                  !the density were equal to dens_moy
                  pgas = min(sn_p_local / pnorm_sn * rr /  dgas , Vsat) * dgas
                  pgas_check = pgas_check + pgas * vol_loc

                  ekin=0.
                  do idim=1,ndim
                     ekin = ekin + ( (uold(ind_cell(i),idim+1))**2 ) / dgas / 2.
                  enddo

                  uold(ind_cell(i), 2+ndim) = uold(ind_cell(i), 2+ndim) - ekin

                  do idim=1,ndim
                     uold(ind_cell(i),idim+1) = uold(ind_cell(i),idim+1) + pgas * (xx(i,idim) - x_sn(idim)) / (rr * sn_r)
                  enddo

                  ekin=0.
                  do idim=1,ndim
                     ekin = ekin + ( (uold(ind_cell(i),idim+1))**2 ) / dgas / 2.
                  enddo

                  !before adding thermal energy make sure the temperature is not too high (too small timesteps otherwise)
                  T_sn = (sn_ed / dgas * (gamma-1.) ) * scale_T2
                  T_sn = min( T_sn , Tsat) / scale_T2
                  sn_ed_lim = T_sn * dgas / (gamma-1.)

                  uold(ind_cell(i), 2+ndim) = uold(ind_cell(i), 2+ndim) + ekin + sn_ed_lim

                end if

            end if
          end do
          ! End loop over sublist of cells
        end do
        ! End loop over cells
      end do
      ! End loop over grids
    end do
    ! End loop over levels

#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(mass_sn_tot,mass_sn_tot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(dens_max_loc,dens_max_loc_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(pgas_check,pgas_check_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
    mass_sn_tot_all = mass_sn_tot
    dens_max_loc_all = dens_max_loc
    pgas_check_all = pgas_check
#endif

    if(myid == 1) write(*, *) "SN momentum (injected, expected):", pgas_check_all, sn_p_local
    if(myid == 1) write(*, *) "Physical units:", pgas_check_all * scale_d * scale_l**3 * scale_v, sn_p_local * scale_d * scale_l**3 * scale_v

    !calculate grid effect
    vol_rap = vol_sn / sn_vol

    !TC: should be outputted to the log
    if(myid .eq. 1) then 
       open(103,file='supernovae2.txt',form='formatted',status='unknown',access='append')
         write(103,112) t,x_sn(1),x_sn(2),x_sn(3),dens_max_loc_all,mass_sn_tot_all,vol_rap,pgas_check_all,sn_p_local
       close(103)
    endif

112 format(9e12.4)

  end do ! end of the loop over stellar objects

  call delete_stellar(mark_del)

  ! Update hydro quantities for split cells
  do ilevel = nlevelmax, levelmin, -1
    call upload_fine(ilevel)
    do ivar = 1, nvar
      call make_virtual_fine_dp(uold(1, ivar), ilevel)
    enddo
  enddo

end subroutine make_sn_stellar
!################################################################
!################################################################
!################################################################
!################################################################
subroutine sphere_average(navg, nsph, center, radius, rpow, upow, avg)
    use amr_parameters, only: boxlen, dp, hydro, icoarse_max, icoarse_min &
        & , jcoarse_min, kcoarse_min, levelmin, ndim, ngridmax, nlevelmax &
        & , nvector, twotondim, verbose
    use amr_commons, only: active, ncoarse, son, xg, myid
    use hydro_parameters, only: nvar
    use hydro_commons, only: uold
    use mpi_mod
    implicit none

    ! Integrate quantities over spheres
    ! The integrand is (r / radius(isph))**rpow(iavg) * product(u(ivar)**upow(iavg, ivar), ivar=1:nvar)

    integer, intent(in):: navg                               ! Number of quantities
    integer, intent(in):: nsph                               ! Number of spheres
    real(dp), dimension(1:nsph, 1:ndim), intent(in):: center ! Sphere centers
    real(dp), dimension(1:nsph), intent(in):: radius         ! Sphere radii
    real(dp), dimension(1:navg), intent(in):: rpow           ! Power of radius in the integral
#ifdef SOLVERmhd
    real(dp), dimension(1:navg, 1:nvar+3), intent(in):: upow ! Power of hydro variables in the integral
#else
    real(dp), dimension(1:navg, 1:nvar), intent(in):: upow   ! Power of hydro variables in the integral
#endif
    real(dp), dimension(1:navg, 1:nsph), intent(out):: avg   ! Averages

    integer:: i, ivar, ilevel, igrid, ind, ix, iy, iz, iskip, isph, idim
    integer:: nx_loc, ncache, ngrid
    integer, dimension(1:nvector):: ind_grid, ind_cell
    logical, dimension(1:nvector):: ok

    real(dp):: scale, dx, dx_loc, vol_loc, rr
    real(dp), dimension(1:3):: skip_loc
    real(dp), dimension(1:twotondim, 1:3):: xc
    real(dp), dimension(1:nvector, 1:ndim):: xx

    integer:: info
    real(dp), dimension(1:navg, 1:nsph):: avg_loc
    real(dp), dimension(1:navg):: integrand
    real(dp), dimension(1:navg):: utemp

    if(.not. hydro)return
    if(ndim .ne. 3)return

    if(verbose .and. myid == 1) write(*, *) 'Entering sphere_average'

    ! Mesh spacing in that level
    nx_loc = icoarse_max - icoarse_min + 1
    skip_loc=(/0.0d0,0.0d0,0.0d0/)
    if(ndim>0)skip_loc(1)=dble(icoarse_min)
    if(ndim>1)skip_loc(2)=dble(jcoarse_min)
    if(ndim>2)skip_loc(3)=dble(kcoarse_min)
    scale = boxlen / dble(nx_loc)

    avg_loc = 0.0d0

    do ilevel = levelmin, nlevelmax
        ! Computing local volume (important for averaging hydro quantities)
        dx = 0.5d0**ilevel
        dx_loc = dx * scale
        vol_loc = dx_loc**ndim

        ! Cell center position relative to grid center position
        do ind = 1, twotondim
            iz = (ind - 1) / 4
            iy = (ind - 1 - 4 * iz) / 2
            ix = (ind - 1 - 2 * iy - 4 * iz)
            if(ndim>0) xc(ind, 1) = (dble(ix) - 0.5d0) * dx
            if(ndim>1) xc(ind, 2) = (dble(iy) - 0.5d0) * dx
            if(ndim>2) xc(ind, 3) = (dble(iz) - 0.5d0) * dx
        end do

        ! Loop over grids
        ncache = active(ilevel)%ngrid
        do igrid = 1, ncache, nvector
            ngrid = min(nvector, ncache - igrid + 1)
            do i = 1, ngrid
                ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
            end do

            ! Loop over cells
            do ind = 1, twotondim
                ! Gather cell indices
                iskip = ncoarse + (ind - 1) * ngridmax
                do i = 1, ngrid
                    ind_cell(i) = iskip + ind_grid(i)
                end do

                ! Gather cell center positions
                do i = 1, ngrid
                    xx(i, :) = xg(ind_grid(i), :) + xc(ind, :)
                end do

                ! Rescale position from coarse grid units to code units
                do idim=1,ndim
                   do i=1,ngrid
                      xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
                   end do
                end do

                ! Flag leaf cells
                do i = 1, ngrid
                    ok(i) = (son(ind_cell(i)) == 0)
                end do

                do i = 1, ngrid
                    if(ok(i)) then
                        do isph = 1, nsph
                            rr = sqrt(sum(((xx(i, :) - center(isph, :)) / radius(isph))**2))

                            if(rr < 1.) then
                                integrand = rr**rpow
                                where(abs(rpow) < 1.0d-10) ! Avoid NaNs of the form 0**0
                                    integrand = 1.0d0
                                end where
#ifdef SOLVERmhd
                                do ivar = 1, nvar + 3
#else
                                do ivar = 1, nvar
#endif
                                    utemp(:) = uold(ind_cell(i), ivar)
                                    where(abs(upow(:, ivar)) < 1.0d-10) ! Avoid NaNs of the form 0**0
                                        utemp = 1.0d0
                                    end where
                                    integrand = integrand * utemp**upow(:, ivar)
                                end do
                                avg_loc(:, isph) = avg_loc(:, isph) + vol_loc * integrand
                            endif
                            ! End test on radius
                        end do
                        ! End loop over spheres
                    endif
                    ! End test on leaf cells
                end do
                ! End loop over sublist of cells
            end do
            ! End loop over cells
        end do
        ! End loop over grids
    end do
    ! End loop over levels

#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(avg_loc, avg, navg * nsph, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#else
    avg = avg_loc
#endif
end subroutine sphere_average
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! THIS SECTION DEALS WITH INDIVIDUAL FIXED SOURCES IN THE NAMELIST
! TC: currently not used but we leave it for future reference

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine feedback_fixed(ilevel)
  use amr_commons
  use hydro_parameters
  use sink_feedback_parameters
  implicit none
  integer::ilevel,isn
  !---------------------
  ! Dump a supernova/winds into cells
  ! NOTE: variables internally still labelled "sn", but valid for winds too
  !---------------------
  ! This check should be done already in amr_step, but best to be sure I suppose
  if(FB_on .and. ilevel == levelmin) then
     ! Reset per-step flags/variables used to signal to other routines
     FB_sourceactive = .false.
     FB_dtnew = 0d0
     do isn=1,FB_nsource
        ! Deal with supernovae
        if(FB_sourcetype(isn) .eq. 'supernova') then
          if(t >= FB_start(isn) .and. .not. FB_done(isn)) then
             ! NOTE: FB_done NOT SAVED BY OUTPUT, SO RELAUNCHING AFTER FB_start WILL LAUNCH SUPERNOVA
             if(myid == 1) write (*,*) 'Supernova ',isn,' @ t = ', t, &
                  & 'FB_start =', FB_start(isn)
             call make_fb_fixed(ilevel,isn)
             FB_done(isn) = .true.
          endif
        endif
        ! Deal with winds
        if(FB_sourcetype(isn) .eq. 'wind') then
          if(t >= FB_start(isn) .and. t <= FB_end(isn)) then
            ! Inject winds
            call make_fb_fixed(ilevel,isn)
            ! Log file book-keeping, wind starts
            if(.not. FB_done(isn)) then
              if(myid == 1) write (*,*) 'Wind started ',isn,' @ t = ', t, &
                    & 'FB_start =', FB_start(isn)
              FB_done(isn) = .true.
            endif
          endif
          ! Log file book-keeping, wind ends
          if(t > FB_end(isn) .and. FB_done(isn)) then
            if(myid == 1) write (*,*) 'Wind stopped ',isn,' @ t = ', t, &
                  & 'FB_end =', FB_end(isn)
            FB_done(isn) = .false.
          endif
        endif
     end do
  endif
end subroutine feedback_fixed

subroutine courant_fb_fixed(dtout)
  ! Find fixed source timestep
  ! NOTE: This basically just ensures that the source starts on time
  ! It's more efficient to let the other courant limiters work after that
  ! dtout - output parameter
  use amr_commons
  use hydro_parameters
  use amr_parameters,only:dp
  use sink_feedback_parameters
  implicit none
  
  integer::isn
  real(dp),intent(out)::dtout
  
  ! Loop through each source
  do isn=1,FB_nsource
     ! Check to see whether we're going to overshoot the supernova
!     if (t+dtout > FB_start(isn)) then
!        dtout = FB_start(isn) - t
!     endif
     ! OLD CODE FOR LIMITING TIMING, USE OTHER COURANT LIMITERS INSTEAD
     ! TODO: REMOVE
     ! Is the source active this timestep?
     !if (FB_sourceactive(isn)) then
     !   ! Set the timestep to the current minimum
     !   if (dtout .gt. 0d0) then
     !      dtout = min(dtout,FB_dtnew(isn))
     !   else
     !      dtout = FB_dtnew(isn)
     !   end if
     !end if
  end do

end subroutine courant_fb_fixed

subroutine make_fb_fixed(currlevel,isn)
  ! Adapted from O. Iffrig's make_sn_blast
  use amr_commons
  use hydro_commons
  use amr_parameters,only:dp
  use sink_feedback_parameters
  implicit none

  integer, intent(in) :: currlevel,isn

  integer:: ivar
  integer:: ilevel, ind, ix, iy, iz, ngrid, iskip, idim
  integer:: i, nx_loc, igrid, ncache
  integer, dimension(1:nvector), save:: ind_grid, ind_cell
  real(dp):: dx, dt
  real(dp):: scale, dx_min, dx_loc, vol_loc

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_msun, scale_ecgs

  real(dp), dimension(1:3):: skip_loc
  real(dp), dimension(1:twotondim, 1:3):: xc
  logical, dimension(1:nvector), save:: ok

  real(dp),dimension(1:ndim):: sn_cent
  real(dp), dimension(1:nvector, 1:ndim), save:: xx
  real(dp):: sn_r, sn_m, sn_e_local, sn_vol, sn_d, sn_ed, dx_sel, sn_p_local, sn_v
  real(dp):: rr, pi
  real(dp), dimension(1:ndim)::rvec
  logical:: sel = .false.
  real(dp),parameter::m_sun=1.9891d33  ! Solar mass [g]
  real(dp),parameter::year2=3.154d7  ! 1 year [s]

  if(.not. hydro)return
  if(ndim .ne. 3)return

  if(verbose)write(*,*)'Entering make_fb_fixed'

  pi = acos(-1.0)

  ! Make source active this timestep
  FB_sourceactive(isn) = .true.

  ! Mesh spacing in that level
  nx_loc = icoarse_max - icoarse_min + 1
  skip_loc = (/ 0.0d0, 0.0d0, 0.0d0 /)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale = boxlen / dble(nx_loc)
  dx_min = scale * 0.5d0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_d * scale_l**ndim / m_sun
  scale_ecgs = scale_d * scale_v**2 * scale_l**ndim

  ! Hard-code sn properties to centre of box
  sn_r = FB_radius(isn)*(0.5**nlevelmax)*scale

  !Set up injection energy, mass, velocity and position
  sn_m = FB_mejecta(isn) / scale_msun ! Put in 10 solar masses
  sn_e_local = FB_energy(isn) / scale_ecgs
  sn_v = sqrt(2.0*(sn_e_local*scale_ecgs)/(sn_m*scale_msun*m_sun))
  sn_v = sn_v / scale_v
  sn_cent(1)= FB_pos_x(isn)*boxlen
  sn_cent(2)= FB_pos_y(isn)*boxlen
  sn_cent(3)= FB_pos_z(isn)*boxlen

  ! HACK - force Courant condition on winds before first hydro step
  ! TODO: think about how this is synchronised more carefully
  FB_dtnew(isn) = courant_factor*dx_min/sn_v
  dt = min(dtnew(currlevel),dt)
  dt = dtnew(currlevel)

  ! If this is a wind, scale the luminosity by the timestep
  if (FB_sourcetype(isn) .eq. 'wind') then
    sn_m = sn_m * dt*scale_t/year2 ! dt in years
    sn_e_local = sn_e_local * dt*scale_t/year2
  endif

  ! HACK !!! - KINETIC BLAST ONLY WORKS FOR sn_r > 0.0 !!!
  if(sn_r /= 0.0) then
     sn_vol = 4. / 3. * pi * sn_r**3
     sn_d = sn_m / sn_vol
     sn_ed = sn_e_local / sn_vol
     sn_p_local = sn_d*sn_v ! uniform momentum of blast ejecta
  end if
     
  if(myid .eq. 1 .and. FB_sourcetype(isn) .eq. 'supernova') then
     write(*,*) 'Supernova blast! Wow!'
     write(*,*) 'x_sn, y_sn, z_sn, ',sn_cent(1),sn_cent(2),sn_cent(3)
  endif

  ! Loop over levels
  do ilevel = levelmin, nlevelmax
    ! Computing local volume (important for averaging hydro quantities)
    dx = 0.5d0**ilevel
    dx_loc = dx * scale
    vol_loc = dx_loc**ndim

    ! Cell center position relative to grid center position
    do ind=1,twotondim
      iz = (ind - 1) / 4
      iy = (ind - 1 - 4 * iz) / 2
      ix = (ind - 1 - 2 * iy - 4 * iz)
      if(ndim>0) xc(ind,1) = (dble(ix) - 0.5d0) * dx
      if(ndim>1) xc(ind,2) = (dble(iy) - 0.5d0) * dx
      if(ndim>2) xc(ind,3) = (dble(iz) - 0.5d0) * dx
    end do

    ! Loop over grids
    ncache=active(ilevel)%ngrid
    do igrid = 1, ncache, nvector
      ngrid = min(nvector, ncache - igrid + 1)
      do i = 1, ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
      end do

      ! Loop over cells
      do ind = 1, twotondim
        ! Gather cell indices
        iskip = ncoarse + (ind - 1) * ngridmax
        do i = 1, ngrid
          ind_cell(i) = iskip + ind_grid(i)
        end do

        ! Gather cell center positions
        do i = 1, ngrid
          xx(i, :) = xg(ind_grid(i), :) + xc(ind, :)
        end do
        ! Rescale position from coarse grid units to code units
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
           end do
        end do
        ! Flag leaf cells
        do i = 1, ngrid
          ok(i) = (son(ind_cell(i)) == 0)
        end do

        do i = 1, ngrid
          if(ok(i)) then
            if(sn_r == 0.0) then
              sn_d = sn_m / vol_loc ! XXX
              sn_ed = sn_e_local / vol_loc ! XXX
              rr = 1.0
              do idim = 1, ndim
                !rr = rr * max(1.0 - abs(xx(i, idim) - sn_center(sn_i, idim)) / dx_sel, 0.0)
                rr = rr * max(1.0 - abs(xx(i, idim) - sn_cent(idim)) / dx_loc, 0.0)
              end do
              !if(rr > 0.0) then
                !if(.not. sel) then
                  !! We found a leaf cell near the supernova center
                  !sel = .true.
                  !sn_d = sn_m / sn_vol
                  !sn_ed = sn_e_local / sn_vol
                !end if
                uold(ind_cell(i), 1) = uold(ind_cell(i), 1) + sn_d * rr
                uold(ind_cell(i), 2+ndim) = uold(ind_cell(i), 2+ndim) + sn_ed * rr
              !end if
            else
               ! Get direction to point the explosion in
               do idim = 1, ndim
                  rvec(idim) = (xx(i, idim) - sn_cent(idim)) / sn_r
               enddo
               rr = sqrt(sum(rvec**2))
               rvec = rvec/rr
   
               if(rr < 1.) then
                  uold(ind_cell(i), 1) = uold(ind_cell(i), 1) + sn_d
                  ! If not entirely thermal injection, add some velocity
                  if (.not.FB_thermal(isn)) then
                     do idim=1,ndim 
                        uold(ind_cell(i),1+idim) = uold(ind_cell(i),1+idim) + &
                          & sn_p_local * rvec(idim)
                     enddo
                  end if
                  uold(ind_cell(i), 2+ndim) = uold(ind_cell(i), 2+ndim) + sn_ed
               endif
            endif
          endif
        end do
      end do
      ! End loop over cells
    end do
    ! End loop over grids
  end do
  ! End loop over levels

  ! Update hydro quantities for split cells
  do ilevel = nlevelmax, levelmin, -1
    call upload_fine(ilevel)
    do ivar = 1, nvar
      call make_virtual_fine_dp(uold(1, ivar), ilevel)
    enddo
  enddo
end subroutine make_fb_fixed

!################################################################
!################################################################
!################################################################
!################################################################

! TC: currently not used but we leave it for future reference
SUBROUTINE feedback_refine(xx,ok,ncell,ilevel)

! This routine flags cells immediately around SN sources to the finest
! level of refinement. The criteria for refinement at a point are:
! a) The point is less than one ilevel cell width from an SN source.
! b) The point is within FB_r_wind finest level cell widths from
!    the SN source.
!-------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use amr_parameters, only:dp
  use sink_feedback_parameters
  implicit none
  integer::ncell,ilevel,i,k,nx_loc,isn
  real(dp),dimension(1:nvector,1:ndim)::xx
  logical ,dimension(1:nvector)::ok
  real(dp)::dx_loc,rvec(ndim),w,rmag,rFB
!-------------------------------------------------------------------------
  nx_loc=(icoarse_max-icoarse_min+1)
  dx_loc = boxlen*0.5D0**ilevel/dble(nx_loc)
  ! Loop over regions
#if NDIM==3
  do isn=1,FB_nsource
     do i=1,ncell
        rFB = FB_r_refine(isn)*boxlen*0.5D0**nlevelmax
        rvec(1)=xx(i,1)-FB_pos_x(isn)*boxlen
        rvec(2)=xx(i,2)-FB_pos_y(isn)*boxlen
        rvec(3)=xx(i,3)-FB_pos_z(isn)*boxlen
        rmag=sqrt(sum(rvec**2))
        if(rmag .le. 2*rFB+dx_loc) then
           ok(i)=.true.
        endif
     end do
  end do
#endif
  
END SUBROUTINE feedback_refine
