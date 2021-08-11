module set_inputs

  use set_precision,   only : prec
  use set_constants,   only : zero, one, half, pi
  use fluid_constants, only : R_gas, gamma

  implicit none

  private

  public :: imax, neq, xmin, xmax, n_ghost_cells
  public :: i_high, i_low, ig_high, ig_low
  public :: areaStar, area, darea
  public :: CFL, eps, tol, eps_roe, beta_lim, epsM, kappaM
  public :: max_iter, max_newton_iter, newton_tol, counter
  public :: soln_save, res_save, res_out
  public :: p0, T0, a0, rho0
  public :: set_derived_inputs, flux_scheme, limiter_scheme, cons
  public :: leftV, rightV, leftU, rightU, limiter_freeze, psi_plus, psi_minus

  integer :: imax    = 128
  integer :: i_low   = 10
  integer :: i_high  = 10
  integer :: ig_low  = 10
  integer :: ig_high = 10
  integer :: neq  = 3
  integer :: max_iter = 150000
  integer :: n_ghost_cells   = 2
  integer :: max_newton_iter = 1000

  real(prec) :: newton_tol = 1.0e-15_prec
  real(prec) :: tol        = 1.0e-11_prec
  real(prec) :: eps        = 1.0e-3_prec
  real(prec) :: p0         = one
  real(prec) :: T0         = one
  real(prec) :: a0         = one
  real(prec) :: rho0       = one
  real(prec) :: areaStar   = 0.2_prec
  real(prec) :: xmin       = -one
  real(prec) :: xmax       = one
  real(prec) :: CFL        = 0.1_prec
  integer :: flux_scheme   = 1
  integer :: limiter_scheme = 2
  real(prec) :: beta_lim = 2
  integer :: soln_save     = 150000
  integer :: res_save      = 10
  integer :: res_out       = 1000
  real(prec) :: eps_roe    = 0.1_prec
  real(prec) :: epsM       = zero
  real(prec) :: kappaM     = -one
  logical :: limiter_freeze = .false.
  logical :: cons           = .true.
  real(prec), dimension(:,:), allocatable :: leftV, rightV, leftU, rightU
  real(prec), dimension(:,:), allocatable :: psi_plus, psi_minus
  integer :: counter = 1
  contains

  !=================================== area ==================================80
  !>
  !! Description: Calculates area distribution for nozzle.
  !!
  !! Inputs:      x:   Position coordinate along x-axis.
  !!
  !! Outputs:     area: Area at specified coordinate.
  !<
  !===========================================================================80
  function area(x)

    real(prec) :: area
    real(prec), intent(in) :: x

    area = 0.1_prec !0.2_prec + 0.4_prec*( one + sin(pi*(x-0.5_prec)))

  end function area

  !=================================== darea =================================80
  !>
  !! Description: Calculates area derivative distribution for nozzle.
  !!
  !! Inputs:      x:   Position coordinate along x-axis.
  !!
  !! Outputs:     darea: analytic dA/dx evaluated at specified coordinate.
  !<
  !===========================================================================80
  function darea(x)

    real(prec) :: darea
    real(prec), intent(in) :: x

    darea = zero !0.4_prec*pi*cos(pi*(x-0.5_prec))

  end function darea

  !=========================== set_derived_inputs ============================80
  !>
  !! Description: Sets derived quantities and prints to STDOUT.
  !<
  !===========================================================================80
  subroutine set_derived_inputs

    a0   = sqrt(gamma*R_gas*T0)
    rho0 = 1000.0_prec*p0/(R_gas*T0)
    i_low = 1
    i_high = imax
    ig_low  = 1 - n_ghost_cells
    ig_high = imax + n_ghost_cells
    write(*,'(A8,F20.14,A13)') 'R     = ', R_gas, ' [J/(kmol*K)]'
    write(*,'(A8,F20.14)')     'gamma = ', gamma
    write(*,'(A8,F20.14,A6)')  'a_0   = ', a0, ' [m/s]'
    write(*,'(A8,F20.14,A9)')  'rho_0 = ', rho0, ' [kg/m^3]'
    write(*,'(A8,F20.14,A6)')  'P_0   = ', p0, ' [kPa]'
    write(*,'(A8,F20.14,A4)')  'T_0   = ', T0, ' [K]'
    write(*,'(A8,F20.14,A6)')  'A*    = ', areaStar, ' [m^2]'

    allocate(leftV(i_low-1:i_high,1:neq))
    allocate(rightV(i_low-1:i_high,1:neq))
    allocate(leftU(i_low-1:i_high,1:neq))
    allocate(rightU(i_low-1:i_high,1:neq))
    allocate(psi_plus(ig_low:ig_high,1:neq))
    allocate(psi_minus(ig_low:ig_high,1:neq))

  end subroutine set_derived_inputs

end module set_inputs
