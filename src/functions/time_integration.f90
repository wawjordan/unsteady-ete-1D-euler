module time_integration

  use set_precision, only : prec
  use set_constants, only : one, half
  use set_inputs,    only : neq, i_low, i_high, ig_low, ig_high
  use grid_type
  use soln_type
  use variable_conversion

  implicit none

  contains

  !============================= calc_time_step ==============================80
  !>
  !! Description:
  !!
  !! Inputs:      dx :
  !!              V  :
  !!
  !! Outputs:     lambda :
  !!              dt     :
  !<
  !===========================================================================80
  subroutine calc_time_step( dx, V, asnd, lambda, dt )

    use set_inputs,          only : CFL

    real(prec), dimension(ig_low:ig_high,neq), intent(in)  :: V
    real(prec), intent(in) :: dx
    real(prec), dimension(ig_low:ig_high),   intent(out) :: lambda
    real(prec), dimension(i_low:i_high),   intent(out) :: dt

    real(prec), dimension(ig_low:ig_high), intent(in)    :: asnd

    lambda(:) = abs(V(:,2)) + asnd
    dt(:) = CFL*dx/lambda(i_low:i_high)
    dt(:) = minval(CFL*dx/lambda(i_low:i_high))

  end subroutine calc_time_step

  !============================= explicit_euler ==============================80
  !>
  !! Description:
  !!
  !! Inputs:      grid :
  !!              S    :
  !!              F    :
  !!              dt   :
  !!              U    :
  !!
  !! Outputs:     lambda :
  !!              dt     :
  !!              R      :
  !!              U      :
  !<
  !===========================================================================80

  subroutine explicit_euler( grid, src, dt, F, U, R )

    type(grid_t), intent(in) :: grid

    real(prec), dimension(ig_low:ig_high,neq), intent(inout) :: U
    real(prec), dimension(i_low:i_high,neq),   intent(out)   :: R
    real(prec), dimension(i_low-1:i_high,neq), intent(in)    :: F
    real(prec), dimension(i_low:i_high)  ,   intent(in)    :: src, dt

    !integer :: i

    R(i_low:i_high,1) = F(i_low:i_high,1)*grid%Ai(i_low:i_high) &
           - F(i_low-1:i_high-1,1)*grid%Ai(i_low-1:i_high-1)
    R(i_low:i_high,2) = F(i_low:i_high,2)*grid%Ai(i_low:i_high) &
           - F(i_low-1:i_high-1,2)*grid%Ai(i_low-1:i_high-1) &
           - src(i_low:i_high)*grid%dx
    R(i_low:i_high,3) = F(i_low:i_high,3)*grid%Ai(i_low:i_high) &
           - F(i_low-1:i_high-1,3)*grid%Ai(i_low-1:i_high-1)


    U(i_low:i_high,1) = U(i_low:i_high,1) &
                      - dt(i_low:i_high)/ &
                        (grid%Ac(i_low:i_high)*grid%dx)*R(:,1)
    U(i_low:i_high,2) = U(i_low:i_high,2) &
                      - dt(i_low:i_high)/ &
                        (grid%Ac(i_low:i_high)*grid%dx)*R(:,2)
    U(i_low:i_high,3) = U(i_low:i_high,3) &
                      - dt(i_low:i_high)/ &
                        (grid%Ac(i_low:i_high)*grid%dx)*R(:,3)

  end subroutine explicit_euler

  !============================= residual_norms ==============================80
  !>
  !! Description:
  !!
  !! Inputs:      R :
  !!
  !! Outputs:     Rnorm :
  !<
  !===========================================================================80
  subroutine residual_norms( R, Rnorm, pnorm, rinit )

    real(prec), dimension(i_low:i_high,1:neq), intent(in) :: R
    real(prec), dimension(1,1:neq), intent(in)  :: rinit
    real(prec),  intent(inout) :: Rnorm(1,1:neq)
    integer :: pnorm

    if (pnorm == 0) then
      Rnorm(1,1:neq) = maxval(abs(R),1)
    elseif (pnorm == 1) then
      Rnorm(1,1) = (one/real(size(R,1)))*sum(abs(R(:,1)),1)
      Rnorm(1,2) = (one/real(size(R,1)))*sum(abs(R(:,2)),1)
      Rnorm(1,3) = (one/real(size(R,1)))*sum(abs(R(:,3)),1)
    elseif (pnorm == 2) then
      Rnorm(1,1) = sqrt((one/real(size(R,1)))*sum(R(:,1)**2,1))
      Rnorm(1,2) = sqrt((one/real(size(R,1)))*sum(R(:,2)**2,1))
      Rnorm(1,3) = sqrt((one/real(size(R,1)))*sum(R(:,3)**2,1))
    else
      Rnorm(1,1) = maxval(abs(R(:,1)))
      Rnorm(1,2) = maxval(abs(R(:,2)))
      Rnorm(1,3) = maxval(abs(R(:,3)))
    end if
     Rnorm(1,1) = Rnorm(1,1)/rinit(1,1)
     Rnorm(1,2) = Rnorm(1,2)/rinit(1,2)
     Rnorm(1,3) = Rnorm(1,3)/rinit(1,3)

  end subroutine residual_norms

!subroutine advance_solution(grid,soln)
!
!  type(grid_t) :: grid
!  type(soln_t) :: soln
!
!  call calculate_time_step(soln%lambda,soln%dt,soln%V)
!
!  call enforce_bndry(soln)
!
!  call extrapolate_variables(soln%U,Uextrap,...)
!
!  call central_flux(Uextrap,Fextrap)
!
!  call jst_damping(soln%lambda,soln%D,soln%E,Fextrap,...)
!
!  call calculate_sources(soln%S,...)
!
!  call explicit_euler(grid,soln%S,soln%dt,Fextrap,soln%U,soln%R)
!
!  call cons2prim(soln%U,soln%V)
!
!  call isentropic_relations(soln%M,soln%V,soln%T) ! calc_mach_number?
!
!  call output_soln()
!
!  call enforce_bndry(soln)
!
!end subroutine advance_solution
end module time_integration
