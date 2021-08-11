module exact_soln_type

  use set_precision      , only : prec
  use set_constants      , only : zero, one, two, half
  use fluid_constants    , only : R_gas, gamma
  use set_inputs         , only : imax, neq, n_ghost_cells, eps, ig_low, ig_high
  use set_inputs         , only : max_newton_iter, newton_tol, i_low, i_high
  use variable_conversion, only : isentropic_relations, prim2cons
  use riemann_problem    , only : init_bracket, solve_riemann, &
              cell_avg_riemann, sample_riemann
  use grid_type

  implicit none

  private

  public :: exact_soln_t, setup_exact_soln, destroy_exact_soln, &
            sample_exact_soln

!============================== exact_soln_t =================================80
!>
!! Description: Derived type for analytic solutions.
!<
!=============================================================================80
  type exact_soln_t

    real(prec), allocatable, dimension(:)   :: aq ! speed of sound
    real(prec), allocatable, dimension(:)   :: Mq ! Mach #
    real(prec), allocatable, dimension(:)   :: Tq ! Temperature
    real(prec), allocatable, dimension(:,:) :: Vq ! (rho,u,p)^T
    real(prec), allocatable, dimension(:,:) :: Uq ! (rho,u,p)^T
    real(prec), allocatable, dimension(:) :: left_state, right_state
    real(prec) :: pStar, uStar, aStar, areaStar

  end type exact_soln_t
  contains

    !========================= setup_exact_soln ==============================80
    !>
    !! Description: Sets up exact_soln_t type.
    !<
    !=========================================================================80
    subroutine setup_exact_soln( this, grid, U_L, U_R )

      use set_constants, only : zero, one
      use riemann_problem, only : solve_riemann

      type(exact_soln_t), intent(inout) :: this
      type(grid_t), intent(in) :: grid
      real(prec), intent(in) :: U_L(:), U_R(:)

      allocate( this%aq(ig_low:ig_high), &
                this%Mq(ig_low:ig_high), &
                this%Tq(ig_low:ig_high), &
                this%Vq(ig_low:ig_high,neq), &
                this%Uq(ig_low:ig_high,neq) )
      this%aq = zero
      this%Mq = zero
      this%Tq = zero
      this%Vq = zero
      this%Uq = zero

      this%pStar = one
      this%uStar = one
      this%aStar = one
      this%areaStar = one
      this%left_state  = U_L
      this%right_state = U_R

      call solve_riemann_problem(this, grid)

    end subroutine setup_exact_soln

    subroutine sample_exact_soln(this, grid, t)
      type(exact_soln_t), intent(inout) :: this
      type(grid_t), intent(in) :: grid
      real(prec), intent(in) :: t
      real(prec) :: xlocs(5)
      integer :: xmsk(5)
      !call sample_riemann(this%uStar, this%pStar, &
      !         this%left_state, this%right_state, &
      !         grid%xc, &
      !         t,newton_tol,this%Vq,xlocs,xmsk)
      call cell_avg_riemann(this%uStar, this%pStar, &
               this%left_state, this%right_state, &
               grid%xi, &
               t,newton_tol,this%Vq,xlocs,xmsk)
    end subroutine sample_exact_soln

    !======================== destroy_exact_soln =============================80
    !>
    !! Description: Deallocates exact_soln_t type.
    !<
    !=========================================================================80
    subroutine destroy_exact_soln( this )

      type(exact_soln_t) :: this

      deallocate( this%aq, this%Mq, this%Tq, this%Vq, this%Uq )

    end subroutine destroy_exact_soln

    subroutine solve_riemann_problem( this, grid )

      type(exact_soln_t), intent(inout) :: this
      type(grid_t), intent(in) :: grid
      real(prec) :: p0, bnd_a, bnd_b
      integer :: i

     call init_bracket( this%left_state, this%right_state, &
                        1.0e-6_prec, p0, bnd_a, bnd_b )
     call solve_riemann( p0, bnd_a, bnd_b, newton_tol, max_newton_iter, &
                        this%left_state, this%right_state, this%pStar,  &
                        this%uStar )
     !write(*,*) this%pStar, this%uStar
     !stop

    end subroutine solve_riemann_problem

end module exact_soln_type
