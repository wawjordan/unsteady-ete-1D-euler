module fluid_constants

  use set_precision, only : prec
  use set_constants, only : zero, one

  implicit none

  private

  public :: M_air, R_u, R_gas, gamma
  public :: set_fluid_constants

  real(prec) :: M_air = 28.96_prec
  real(prec) :: R_u   = 8314.0_prec
  real(prec) :: R_gas = zero
  real(prec) :: gamma = 1.4_prec

  contains

  !=========================== set_fluid_constants ===========================80
  !>
  !! Description: Sets derived quantities.
  !!
  !<
  !===========================================================================80
  subroutine set_fluid_constants

    R_gas = R_u/M_air

  end subroutine set_fluid_constants

end module fluid_constants
