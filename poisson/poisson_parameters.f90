module poisson_parameters
  use amr_parameters

  ! Convergence criterion for Poisson solvers
  real(dp)::epsilon=1.0D-4

  ! Type of force computation
  integer ::gravity_type=0

  ! Gravity parameters
  real(dp),dimension(1:10)::gravity_params=0

  ! Cubic spline kernel softening by RMR
  logical :: cubic_spline_kernel = .false.  ! whether or not to use a cubic kernel that reduces to a Newtonian potential outside a given radius
  real(dp):: cubic_kernel_rsoft = 0.0d0  ! softening length for the cubic kernel
  
  ! 2 body problem
  real(dp):: gmass_secondary = 0.
  real(dp):: soft_secondary = 0.

  ! Maximum level for CIC dark matter interpolation
  integer :: cic_levelmax=0

  ! Min level for CG solver
  ! level < cg_levelmin uses fine multigrid
  ! level >=cg_levelmin uses conjugate gradient
  integer :: cg_levelmin=999

  ! Gauss-Seidel smoothing sweeps for fine multigrid
  integer, parameter :: ngs_fine   = 2
  integer, parameter :: ngs_coarse = 2

  ! Number of multigrid cycles for coarse levels *in safe mode*
  !   1 is the fastest,
  !   2 is slower but can give much better convergence in some cases
  integer, parameter :: ncycles_coarse_safe = 1

end module poisson_parameters
