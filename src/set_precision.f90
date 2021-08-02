!============================== set_precision ================================80
!>
!! Provides IEEE 754 compliant kinds via the iso_c_binding module and,
!! by extension, C interoperability.
!! The r* naming convention means real and rd is real_default.
!! The c* naming convention means complex and cd is complex_default.
!! The i* naming convention means integer and id is integer_default.
!! * == number of bytes, the defaults are all double precision
!<
!=============================================================================80
module set_precision

  use iso_c_binding, only : c_float, c_double, c_long_double, c_float_complex, &
                            c_double_complex, c_long_double_complex,           &
                            c_short, c_int, c_long_long

  implicit none

  private

  public :: r4, r8, r16, rd, dp, prec
  public :: c8, c16, c32, cd
  public :: i2, i4, i8, id

! Real
  integer, parameter :: r4  = c_float
  integer, parameter :: r8  = c_double
  integer, parameter :: r16 = c_long_double
  integer, parameter :: dp  = r8
  integer, parameter :: rd  = r8
  integer, parameter :: prec  = r8

! Complex
  integer, parameter :: c8  = c_float_complex
  integer, parameter :: c16 = c_double_complex
  integer, parameter :: c32 = c_long_double_complex
  integer, parameter :: cd  = c16

! Integer
  integer, parameter :: i2 = c_short
  integer, parameter :: i4 = c_int
  integer, parameter :: i8 = c_long_long
  integer, parameter :: id = i4

end module set_precision
!=============================================================================80
