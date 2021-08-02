module init_problem

  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, half
  use fluid_constants, only : R_gas, gamma
  use set_inputs,      only : p0, T0, eps, ig_low, ig_high
  use soln_type
  use grid_type
  use variable_conversion

  implicit none
  private
  public :: initialize

  contains

  !================================ initialize  ==============================80
  !>
  !! Description:
  !!
  !! Inputs:      grid :
  !!              soln :
  !!
  !! Outputs:     grid :
  !!              soln :
  !<
  !===========================================================================80
  subroutine initialize( grid, soln )

    implicit none

    integer :: i
    type( soln_t ), intent(inout) :: soln
    type( grid_t ), intent(inout) :: grid

    where(grid%xc<=0.0_prec)
      soln%V(:,1) = one
      soln%V(:,2) = zero
      soln%V(:,3) = 100000.0_prec
    elsewhere
      soln%V(:,1) = 0.125_prec
      soln%V(:,2) = zero
      soln%V(:,3) = 10000.0_prec
    end where
    !do i = ig_low,ig_high
    !  if ( soln%mach(i) < eps ) then
    !    soln%mach(i) = eps
    !  end if
    !end do

    !call isentropic_relations(soln%mach,soln%V)
    call prim2cons(soln%U,soln%V)
    call update_states( soln )

  end subroutine initialize

end module init_problem
