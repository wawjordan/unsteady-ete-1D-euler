module time_integration

  use set_precision, only : prec
  use set_constants, only : zero, one, half, third, fourth
  use set_inputs,    only : neq, i_low, i_high, ig_low, ig_high
  use grid_type, only : grid_t
  use soln_type, only : soln_t
  use variable_conversion

  implicit none

  private

  public :: calc_time_step, explicit_RK, calc_residual, residual_norms
  public :: explicit_Euler
  contains

  !============================= calc_time_step ==============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine calc_time_step( dx, V, asnd, lambda, dt )

    use set_inputs,          only : CFL

    real(prec), dimension(neq,ig_low:ig_high), intent(in)  :: V
    real(prec), dimension(ig_low:ig_high), intent(in)  :: dx
    real(prec), dimension(ig_low:ig_high),   intent(out) :: lambda
    real(prec), dimension(i_low:i_high),   intent(out) :: dt

    real(prec), dimension(ig_low:ig_high), intent(in)    :: asnd

    lambda(:) = abs(V(2,:)) + asnd
    dt(:) = CFL*dx(i_low:i_high)/lambda(i_low:i_high)
    dt(:) = minval(dt)

  end subroutine calc_time_step

  !============================== explicit_RK ================================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine explicit_RK( grid, soln, VL, VR )
    use flux_calc, only : calc_flux_1D
    use variable_conversion, only : update_states, prim2cons
    use other_subroutines, only : MUSCL_extrap
    use basic_boundaries, only : explicit_characteristic_bndry
    type(grid_t), intent(inout) :: grid
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq), intent(in) :: VL, VR
    real(prec), dimension(neq,i_low-1:i_high) :: leftV, rightV, leftU, rightU
    real(prec), dimension(neq,ig_low:ig_high) :: old_U
    real(prec), dimension(4) :: k
    integer :: i, j, n

    old_U = soln%U

    k = (/ fourth, third, half, one /)
    !k = one
    do j = 1,4
      call explicit_characteristic_bndry(soln, VL, VR)
      call prim2cons(soln%U,soln%V)
      call calc_flux_1D(grid,soln)
      call calc_residual(grid,soln)

      do i = 1,neq
        soln%U(i,i_low:i_high) = &
          old_U(i,i_low:i_high) - k(j)* &
          soln%R(i,i_low:i_high)*soln%dt(i_low:i_high)/ &
          ( grid%Ac(i_low:i_high)*grid%dx(i_low:i_high) )
      end do
      call update_states(soln)
    end do
  end subroutine explicit_RK

  !============================== calc_residual ==============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine calc_residual(grid,soln)

    type(grid_t), intent(inout) :: grid
    type(soln_t), intent(inout) :: soln
    integer :: i

    do i = i_low,i_high
      soln%R(:,i) = grid%Ai(i)*soln%F(:,i) &
                  - grid%Ai(i-1)*soln%F(:,i-1) &
                  - grid%Ac(i)*grid%dx(i)*soln%S(:,i)
    end do
  end subroutine calc_residual


  !============================= explicit_euler ==============================80
  !>
  !! Description:
  !<
  !===========================================================================80

  subroutine explicit_euler( grid, src, dt, F, U, R )

    type(grid_t), intent(in) :: grid

    real(prec), dimension(neq,ig_low:ig_high), intent(inout) :: U
    real(prec), dimension(neq,i_low:i_high),   intent(out)   :: R
    real(prec), dimension(neq,i_low-1:i_high), intent(in)    :: F
    real(prec), dimension(i_low:i_high)  ,   intent(in)    :: src, dt

    !integer :: i

    R(1,i_low:i_high) = F(1,i_low:i_high)*grid%Ai(i_low:i_high) &
           - F(1,i_low-1:i_high-1)*grid%Ai(i_low-1:i_high-1)
    R(2,i_low:i_high) = F(2,i_low:i_high)*grid%Ai(i_low:i_high) &
           - F(2,i_low-1:i_high-1)*grid%Ai(i_low-1:i_high-1) &
           - 0.0_prec*src(i_low:i_high)*grid%dx(i_low:i_high)
    R(3,i_low:i_high) = F(3,i_low:i_high)*grid%Ai(i_low:i_high) &
           - F(3,i_low-1:i_high-1)*grid%Ai(i_low-1:i_high-1)


    U(1,i_low:i_high) = U(1,i_low:i_high) &
                      - dt(i_low:i_high)/ &
                        (grid%Ac(i_low:i_high)*grid%dx(i_low:i_high))*R(1,:)
    U(2,i_low:i_high) = U(2,i_low:i_high) &
                      - dt(i_low:i_high)/ &
                        (grid%Ac(i_low:i_high)*grid%dx(i_low:i_high))*R(2,:)
    U(3,i_low:i_high) = U(3,i_low:i_high) &
                      - dt(i_low:i_high)/ &
                        (grid%Ac(i_low:i_high)*grid%dx(i_low:i_high))*R(3,:)

  end subroutine explicit_euler

  !============================= residual_norms ==============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine residual_norms( R, Rnorm, rinit )

    real(prec), dimension(neq,i_low:i_high), intent(in) :: R
    real(prec), dimension(neq), intent(in)  :: rinit
    real(prec), dimension(neq), intent(out) :: Rnorm
    real(prec) :: Linv
    integer :: i

    Linv = one/real(size(R,2),prec)

    do i = 1,neq
      Rnorm(i) = sqrt(Linv*sum(R(i,:)**2))
    end do
    Rnorm = Rnorm/rinit

  end subroutine residual_norms

end module time_integration
