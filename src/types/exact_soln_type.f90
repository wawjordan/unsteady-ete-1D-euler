module exact_soln_type

  use set_precision      , only : prec
  use set_constants      , only : zero, one, two, half
  use fluid_constants    , only : R_gas, gamma
  use set_inputs         , only : imax, neq, n_ghost_cells, eps
  use set_inputs         , only : max_newton_iter, newton_tol, i_low, i_high
  use variable_conversion, only : isentropic_relations, prim2cons
  use grid_type

  implicit none

  private

  public :: exact_soln_t
  public :: allocate_exact_soln, deallocate_exact_soln
  public :: solve_exact_soln

!============================== exact_soln_t =================================80
!>
!! Description: Derived type for analytic solutions.
!<
!=============================================================================80
  type :: exact_soln_t

    real(prec), allocatable, dimension(:)   :: aq ! speed of sound
    real(prec), allocatable, dimension(:)   :: Mq ! Mach #
    real(prec), allocatable, dimension(:)   :: Tq ! Temperature
    real(prec), allocatable, dimension(:,:) :: Vq ! (rho,u,p)^T
    real(prec), dimension(neq) :: left_state, right_state
    real(prec) :: pStar, uStar, aStar, areaStar

  contains

    procedure, public, pass :: setup_exact_soln
    procedure, public, pass :: destroy_exact_soln
    procedure(solve_exact_soln_i),  pointer, public, nopass :: solve_exact_soln
    procedure(sample_exact_soln_i), pointer, public, nopass :: sample_exact_soln

  end type exact_soln_t

  abstract interface

    !========================== solve_exact_soln_i ===========================80
    !>
    !! Description: Interface for subroutine to solve for exact solution.
    !<
    !=========================================================================80
    subroutine solve_exact_soln_i(this,grid,t)

      use set_precision, only : prec
      use grid_type,     only : grid_t

      import exact_soln_t

      implicit none

      type(exact_soln_t), intent(inout) :: this
      type(grid_t), intente(in) :: grid
      real(prec), intent(in) :: t

    end subroutine solve_exact_soln_i

    !========================== sample_exact_soln_i ==========================80
    !>
    !! Description: Interface for subroutine to sample exact solution values.
    !<
    !=========================================================================80
    subroutine sample_exact_soln_i(this,grid,t)

      use set_precision, only : prec
      use grid_type,     only : grid_t

      import exact_soln_t

      implicit none

      type(exact_soln_t), intent(inout) :: this
      type(grid_t), intente(in) :: grid
      real(prec), intent(in) :: t

    end subroutine sample_exact_soln_i

  end interface

  contains

    !========================= setup_exact_soln ==============================80
    !>
    !! Description: Sets up exact_soln_t type.
    !<
    !=========================================================================80
    subroutine setup_exact_soln( this )

      use set_constants, only : zero, one, SOLUTION_TYPE
      use exact_solvers, only : solve_q1d_isentropic_nozzle, &
                                solve_riemann_problem

      type(exact_soln_t) :: this

      allocate( this%aq(i_low-1:i_high), &
                this%Mq(i_low-1:i_high), &
                this%Tq(i_low-1:i_high), &
                this%Vq(i_low-1:i_high,neq) )
      this%aq = zero
      this%Mq = zero
      this%Tq = zero
      this%Vq = zero

      this%pStar = one
      this%uStar = one
      this%aStar = one
      this%areaStar = one
      this%left_state  = (/ one, one, one /)
      this%right_state = (/ one, one, one /)

      select case( SOLUTION_TYPE )

        case( Q1D )

          this%solve_exact_soln => solve_q1d_isentropic_nozzle

        case( SHOCK_TUBE )

          this%solve_exact_soln => solve_riemann_problem

        case default

      end select

    end subroutine setup_exact_soln

    !======================== destroy_exact_soln =============================80
    !>
    !! Description: Deallocates exact_soln_t type.
    !<
    !=========================================================================80
    subroutine destroy_exact_soln( this )

      type(exact_soln_t) :: this

      deallocate( this%aq, this%Mq, this%Tq, this%Vq )

    end subroutine destroy_exact_soln

    !========================== solve_exact_soln =============================80
    !>
    !! Description: Calculates analytic/quasi-analytic solution.
    !!
    !! Inputs:      soln:    exact_q1d_t type.
    !!
    !! Outputs:     soln:    exact_q1d_t type.
    !<
    !=========================================================================80
  subroutine solve_q1d_isentropic_nozzle( this, grid, t )

    type(exact_q1d_t), intent(inout) :: this
    type(grid_t),      intent(in) :: grid
    real(prec),        intent(in) :: t
    real(prec) :: M0
    real(prec) :: M1
    real(prec), dimension(2) :: Mk
    real(prec) :: err
    integer :: i
!====================== Solution at cell centers =========================80

    Mk = -9999.99_prec
    err  = -9999.99_prec

    M0 = eps
    M1 = one
    call solve_Q1D(grid%Ac(i_low),areaStar,M0,M0,M1,tol,max_iter,this%Mq)
    call newton_safe( grid%Ac(i_low), fun, dfun, M0, M1, soln%Mc(i_low), Mk, err )

    do i = i_low+1,i_high
      if ( (iSS==1).and.(grid%Ai(i) > grid%Ai(i-1)) ) then
        M0 = one - eps
        M1 = 10.0_prec
      else
        M0 = eps
        M1 = one+eps
      endif

      call Solve( grid%Ac(i), fun, dfun, M0, M1, soln%Mc(i), Mk, err )

    end do

    call isentropic_relations( soln%Mc, soln%Vc )
    call prim2cons(soln%Uc,soln%Vc)

  end subroutine solve_exact_q1d

  !================================== fun ====================================80
  !>
  !! Description: Fixed-point form of Area-Mach number relation as a function
  !!              of Mach number and area.
  !!
  !! Inputs:      M:    Mach number.
  !!              A:   Local area.
  !!
  !! Outputs:     fun:    f(M,A).
  !<
  !===========================================================================80
  function fun (M,A)

    real(prec) :: fun
    real(prec), intent (in) :: M
    real(prec), intent (in) :: A
    fun = ((two/(gamma+1))*(one+half*(gamma-1)*M**2))**((gamma+1)/(gamma-1)) - ((A/Astar)**2)*M**2
    return

  end function fun

  !================================= dfun ====================================80
  !>
  !! Description: Derivative of fixed-point form of Area-Mach number relation
  !!              as a function of Mach number and area.
  !!
  !! Inputs:      M:    Mach number.
  !!              A:   Local area.
  !!
  !! Outputs:     dfun:    df(M,A)/dM.
  !<
  !===========================================================================80
  function dfun (M,A)

    real(prec) :: dfun
    real(prec), intent (in) :: M
    real(prec), intent (in) :: A
    dfun = two*M*( ((two/(gamma+1))*(one+half*(gamma-1)*M**2))**(two/(gamma-1)) - (A/Astar)**2 )
    return

  end function dfun


  !============================== newton_safe3 ===============================80
  !>
  !! Description: Safe-guarded 1D Newton's method, but with hard-coded input
  !!              parameter 'A1' and reduced memory overhead because I'm dumb
  !!              and hard-coded array sizes in a previous version of the code.
  !!
  !! Inputs:      A      :  Area.
  !!              fun    :  two variable function f(x,A1).
  !!              dfun   :  Derivative of f(x) w.r.t. x.
  !!              bnd_a1 :  Beginning of bracketing interval.
  !!              bnd_b1 :  End of bracketing interval.
  !!
  !! Outputs:     M      :  zero of f(x).
  !!              Mk     :  Array of last two iterates.
  !!              err    :  approximate relative iterative error.
  !<
  !===========================================================================80
  subroutine newton_safe( A, fun, dfun, bnd_a1, bnd_b1, M, Mk, err )

    real(prec), external :: fun, dfun

    real(prec), intent(in)  :: bnd_a1, bnd_b1, A
    real(prec), intent(out) :: M
    real(prec), intent(out) :: err

    real(prec), dimension(2), intent(out) :: Mk

    real(prec) :: M_new
    real(prec) :: bnd_a, bnd_b
    integer    :: k = 1

    Mk = zero

    bnd_a = bnd_a1
    bnd_b = bnd_b1
    Mk(1)  = bnd_a
    Mk(2)  = Mk(1) - fun(Mk(1),A)/dfun(Mk(1),A)

    if ( (Mk(2) < bnd_a).or.(Mk(2) > bnd_b) ) then
      Mk(2) = bnd_a + half*(bnd_b-bnd_a)
      if ( int(sign(one,fun(bnd_a,A))) == int(sign(one,fun(Mk(2),A))) ) then
        bnd_a = Mk(2)
      else
        bnd_b = Mk(2)
      endif
    endif

    err = abs(Mk(2)-Mk(1))/abs(Mk(2))
    do k = 2, max_newton_iter
      if (err < newton_tol) exit
      M_new = Mk(2) - fun(Mk(2),A)/dfun(Mk(2),A)
      if ( (M_new < bnd_a).or.(M_new > bnd_b) ) then
        M_new = bnd_a + half*(bnd_b-bnd_a)
        if ( int(sign(one,fun(bnd_a,A))) == int(sign(one,fun(M_new,A))) ) then
          bnd_a = M_new
        else
          bnd_b = M_new
        endif
      endif
      Mk(1) = Mk(2)
      Mk(2) = M_new
      err = abs(Mk(2)-Mk(1))/abs(Mk(2))
    end do

    M = Mk(2)

  end subroutine newton_safe

end module exact_q1d_type
