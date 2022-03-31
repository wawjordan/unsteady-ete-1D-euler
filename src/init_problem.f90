module init_problem

  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, half
  use fluid_constants, only : R_gas, gamma
  use set_inputs,      only : p0, T0, eps, ig_low, ig_high
  use soln_type
  use grid_type
  use variable_conversion
  use riemann_problem, only : init_bracket, solve_riemann, cell_avg_riemann

  implicit none
  private
  public :: initialize

  contains

  !================================ initialize  ==============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine initialize( grid, soln, U_L, U_R )

    implicit none

    integer :: i
    type( soln_t ), intent(inout) :: soln
    type( grid_t ), intent(inout) :: grid
    real(prec), intent(in) :: U_L(neq), U_R(neq)
    real(prec) :: p0, bnd_a, bnd_b, pStar, uStar, xstar, t

    where(grid%xc<=0.0_prec)
      soln%V(1,:) = U_L(1)
      soln%V(2,:) = U_L(2)
      soln%V(3,:) = U_L(3)
    elsewhere
      soln%V(1,:) = U_R(1)
      soln%V(2,:) = U_R(2)
      soln%V(3,:) = U_R(3)
    end where
  !  soln%V(1,:) = one-0.1_prec*grid%xc**2
  !  soln%V(2,:) = grid%xc**2-1
  !  soln%V(3,:) = one-0.2_prec*grid%xc**2
    !do i = ig_low, ig_high
    !  soln%V(1,i) = real(i,prec)
    !  soln%V(2,i) = real(i,prec)
    !  soln%V(3,i) = real(i,prec)
    !end do
    !do i = ig_low, ig_high
    !  soln%U(1,i) = real(i,prec)
    !  soln%U(2,i) = real(i,prec)
    !  soln%U(3,i) = real(i,prec)
    !end do

    call prim2cons(soln%U,soln%V)
    call update_states( soln )


  end subroutine initialize

end module init_problem
