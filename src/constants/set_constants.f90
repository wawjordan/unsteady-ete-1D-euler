module set_constants

  use set_precision, only : prec

  implicit none

  private

  public :: zero, one, two, three, four, six
  public :: half, third, fourth, fifth, sixth, tenth
  public :: pi, set_derived_constants

  real(prec), parameter :: zero   = 0.0_prec
  real(prec), parameter :: tenth  = 0.1_prec
  real(prec), parameter :: sixth  = 1.0_prec/6.0_prec
  real(prec), parameter :: fifth  = 0.2_prec
  real(prec), parameter :: fourth = 0.25_prec
  real(prec), parameter :: third  = 1.0_prec/3.0_prec
  real(prec), parameter :: half   = 0.5_prec
  real(prec), parameter :: one    = 1.0_prec
  real(prec), parameter :: two    = 2.0_prec
  real(prec), parameter :: three  = 3.0_prec
  real(prec), parameter :: four   = 4.0_prec
  real(prec), parameter :: six    = 6.0_prec
  real(prec)            :: pi     = 3.0_prec

  contains

    !======================== set_derived_constants ==========================80
    !>
    !! Description: Sets derived quantities to specified precision.
    !<
    !=========================================================================80
    subroutine set_derived_constants

      implicit none

      pi = acos(-1.0_prec)

    end subroutine set_derived_constants

end module set_constants
