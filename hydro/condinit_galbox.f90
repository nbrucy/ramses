module galbox_module
  use amr_parameters
  !================================================================
  ! This module contains the variable needed for galbox IC
  !================================================================

  real(dp),save::turb=0.
  real(dp),save::dens0=0.
  real(dp),save::Height0=0.
  real(dp),save::Bx=0.,By=0.,Bz=0.

end module galbox_module


subroutine read_galbox_params()
  use galbox_module
  implicit none

  character(LEN=80)::infile

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/galbox_params/turb, Height0, dens0,Bx,By,Bz

  ! Read namelist file
  call getarg(1,infile) ! get the name of the namelist
  open(1, file=infile)
  read(1, NML=galbox_params)
  close(1)

end subroutine read_galbox_params


subroutine condinit_galbox(x,u,dx,nn)

  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================

  use amr_parameters
  use amr_commons
  use hydro_parameters
  use galbox_module
  implicit none

  ! amr data
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size

#ifdef SOLVERmhd
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
#else
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
#endif

  real(dp),dimension(1:nvector,1:ndim)::x ! Position of cell center

  integer::ivar,i,j,k
  real(dp)::pi,xx,yy
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,mag_norm

  logical,save:: first_call = .true.
  real(dp),dimension(1:3,1:100,1:100,1:100),save::q_idl
  real(dp),save::vx_tot,vy_tot,vz_tot,vx2_tot,vy2_tot,vz2_tot,vx,vy,vz,v_rms
  integer,save:: n_size
  integer:: ind_i, ind_j, ind_k
  real(dp),save:: ind,seed1,seed2,seed3,xi,yi,zi
  real(dp):: n_total
  real(dp):: temper


  !start initialising q and u to zero
  u=0.
  q=0.

  !for mass_sph refinnement - 10. as a reference for coarse level (convenient when in cc for instance)
  mass_sph = 10. * (boxlen*(0.5**levelmin))**3

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  mag_norm = sqrt(1.*8000./scale_T2*2.*1.5)



!!! Step 1 : Read params and initialize turbulent velocity field
if(first_call) then
     
     if(myid .eq. 1) write(*,*) 'mag_norm', mag_norm

     call read_galbox_params()

     if(myid .eq. 1) then 
        write(*,*) 'turb, Height0, dens0,Bx,By,Bz',turb, Height0, dens0,Bx,By,Bz
     endif

     ! Read the turbulent velocity field used as initial condition
     if(myid ==1) write(*,*) '[condinit] Read the file which contains the initial turbulent velocity field'
     open(20, file='ramses.data', form='formatted')
     read(20, *) n_size, ind, seed1, seed2, seed3

     if(n_size .ne. 100) then
        write(*,*) '[Error] [condinit] Unexpected field size in initial turbulence'
        call clean_stop
     endif

     v_rms=0.

     vx_tot  = 0.
     vy_tot  = 0.
     vz_tot  = 0.
     vx2_tot = 0.
     vy2_tot = 0.
     vz2_tot = 0.

     do k=1,n_size
        do j=1,n_size
           do i=1,n_size
              read(20, *) xi, yi, zi, vx, vy, vz
              q_idl(1,i,j,k) = vx
              q_idl(2,i,j,k) = vy
              q_idl(3,i,j,k) = vz

              xi = boxlen*((i-0.5) / n_size-0.5)
              yi = boxlen*((j-0.5) / n_size-0.5)
              zi = boxlen*((k-0.5) / n_size-0.5)

              vx_tot = vx_tot + vx
              vy_tot = vy_tot + vy
              vz_tot = vz_tot + vz

              vx2_tot = vx2_tot + vx**2
              vy2_tot = vy2_tot + vy**2
              vz2_tot = vz2_tot + vz**2

           enddo
        enddo
     enddo
     close(20)

     n_total = n_size**3

     vx_tot = vx_tot / n_total
     vy_tot = vy_tot / n_total
     vz_tot = vz_tot / n_total

     vx2_tot = vx2_tot / n_total
     vy2_tot = vy2_tot / n_total
     vz2_tot = vz2_tot / n_total

     v_rms = sqrt(vx2_tot-vx_tot**2 + vy2_tot-vy_tot**2 + vz2_tot-vz_tot**2)

     ! Calculate the coefficient by which the turbulence velocity needs
     ! to be multiplied
     ! turb is in km/s ,  1.d5 converts it in cm/s
     v_rms =  turb*1.d5 / scale_v / v_rms

     if (myid ==1 ) write(*,*) 'turb ', turb, ', v_rms ', v_rms , 'first_call ',first_call

100  format(i5,4e10.5)
101  format(6e10.5)
102  format(i5)

     if (myid ==1)  write(*,*) '[condinit] Reading achieved'
     first_call = .false.
  end if

!!! Step 2: Initialize primitive field
  do i=1,nn

     x(i,1) = x(i,1) - 0.5 * boxlen
     x(i,2) = x(i,2) - 0.5 * boxlen
     x(i,3) = x(i,3) - 0.5 * boxlen


#ifdef SOLVERmhd
     !Bx component 
     q(i,6     ) = Bx * mag_norm * exp(-x(i,3)**2/(2.*Height0**2)) !exponential profile along z
     q(i,nvar+1) = q(i,6)

     !By component
     q(i,7     ) =  By * mag_norm
     q(i,nvar+2) =  q(i,7)

     !Bz component
     q(i,8     ) =  Bz * mag_norm
     q(i,nvar+3) =  q(i,8)
#endif


     ! in cgs

     ! density
     q(i,1) = dens0 * max(exp(-x(i,3)**2 / (2.*Height0**2)), 1.d-2) ! exponential profile along z
     ! pressure
     temper = (8000. / scale_T2 ) / dens0
     q(i,5) =  q(i,1)*temper

     ! initialise the turbulent velocity field
     ! make a zero order interpolation (should be improved)
     ind_i = int((x(i,1) / boxlen + 0.5) * n_size) + 1
     ind_j = int((x(i,2) / boxlen + 0.5) * n_size) + 1
     ind_k = int((x(i,3) / boxlen + 0.5) * n_size) + 1


     if( ind_i .lt. 1 .or. ind_i .gt. n_size) write(*,*) 'ind_i ',ind_i,boxlen,x(i,1),n_size
     if( ind_j .lt. 1 .or. ind_j .gt. n_size) write(*,*) 'ind_j ',ind_j
     if( ind_k .lt. 1 .or. ind_k .gt. n_size) write(*,*) 'ind_k ',ind_k

     q(i,2) = v_rms*(q_idl(1,ind_i,ind_j,ind_k)-vx_tot)
     q(i,3) = v_rms*(q_idl(2,ind_i,ind_j,ind_k)-vy_tot)
     q(i,4) = v_rms*(q_idl(3,ind_i,ind_j,ind_k)-vz_tot)
  end do





!!! Step 3: Convert primitive to conservative variables

  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)


#ifdef SOLVERmhd
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  u(1:nn,6:8)=q(1:nn,6:8)
#endif


#ifdef SOLVERmhd
  ! passive scalars
  do ivar=9,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#else 
  ! passive scalars
  do ivar=6,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit_galbox

