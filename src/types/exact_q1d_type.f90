module exact_q1d_type

  use set_precision      , only : prec
  use set_constants      , only : zero, one, two, half
  use fluid_constants    , only : R_gas, gamma
  use set_inputs         , only : imax, neq, n_ghost_cells, areaStar, iSS, eps
  use set_inputs         , only : max_newton_iter, newton_tol, i_low, i_high
  use variable_conversion, only : isentropic_relations, prim2cons
  use grid_type

  implicit none

  private

  public :: exact_q1d_t
  public :: allocate_exact_q1d, deallocate_exact_q1d
  public :: solve_exact_q1d

!============================== exact_q1d_t ==================================80
!>
!! Description: Derived type for analytic solution for quasi-1D nozzle flow.
!<
!=============================================================================80
  type exact_q1d_t

    real(prec), allocatable, dimension(:)   :: Mc, Mi
    real(prec), allocatable, dimension(:)   :: Tc, Ti
    real(prec), allocatable, dimension(:,:) :: Vc, Vi
    real(prec), allocatable, dimension(:,:) :: Uc

  end type exact_q1d_t


contains

  !======================== allocate_exact_q1d ===============================80
  !>
  !! Description: Allocates exact_q1d_t type.
  !!
  !! Inputs:      soln:    exact_q1d_t type (unallocated).
  !!
  !! Outputs:     soln:    exact_q1d_t type (allocated).
  !<
  !===========================================================================80
  subroutine allocate_exact_q1d( soln )

    type(exact_q1d_t), intent(inout) :: soln

    allocate( soln%Mc(i_low:i_high), &
              soln%Tc(i_low:i_high), &
              soln%Vc(i_low:i_high,neq), &
              soln%Uc(i_low:i_high,neq) )
    allocate( soln%Mi(i_low-1:i_high), &
              soln%Ti(i_low-1:i_high), &
              soln%Vi(i_low-1:i_high,neq) )

    soln%Mc = zero
    soln%Tc = zero
    soln%Uc = zero
    soln%Vc = zero

    soln%Mi = zero
    soln%Ti = zero
    soln%Vi = zero

  end subroutine allocate_exact_q1d

  !====================== deallocate_exact_q1d ===============================80
  !>
  !! Description: Deallocates exact_q1d_t type.
  !!
  !! Inputs:      soln:    exact_q1d_t type (allocated).
  !!
  !! Outputs:     soln:    exact_q1d_t type (deallocated).
  !<
  !===========================================================================80
  subroutine deallocate_exact_q1d( soln )

    type(exact_q1d_t), intent(inout) :: soln

    deallocate( soln%Mc, soln%Tc, soln%Vc, soln%Uc )
    deallocate( soln%Mi, soln%Ti, soln%Vi )

  end subroutine deallocate_exact_q1d

  !=========================== solve_exact_q1d ===============================80
  !>
  !! Description: Calculates Mach number in nozzle from area relation.
  !!
  !! Inputs:      soln:    exact_q1d_t type.
  !!
  !! Outputs:     soln:    exact_q1d_t type.
  !<
  !===========================================================================80
  subroutine solve_exact_q1d( soln, grid )

    type(exact_q1d_t), intent(inout) :: soln
    type(grid_t), intent(in) :: grid
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
    call newton_safe( grid%Ac(i_low), fun, dfun, M0, M1, soln%Mc(i_low), Mk, err )

    do i = i_low+1,i_high
      if ( (iSS==1).and.(grid%Ai(i) > grid%Ai(i-1)) ) then
        M0 = one - eps
        M1 = 10.0_prec
      else
        M0 = eps
        M1 = one+eps
      endif

      call newton_safe( grid%Ac(i), fun, dfun, M0, M1, soln%Mc(i), Mk, err )

    end do

    call isentropic_relations( soln%Mc, soln%Vc )
    call prim2cons(soln%Uc,soln%Vc)

!====================== Solution at cell interfaces =========================80

    Mk   = -9999.99_prec
    err  = -9999.99_prec

    M0 = eps
    M1 = one
    call newton_safe( grid%Ai(i_low-1), fun, dfun, M0, M1, soln%Mi(i_low-1), Mk, err )

    do i = i_low,i_high
      if ( (iSS==1).and.(grid%Ai(i) > grid%Ai(i-1)) ) then
        M0 = one - eps
        M1 = 10.0_prec
      else
        M0 = eps
        M1 = one+eps
      endif

      call newton_safe( grid%Ai(i), fun, dfun, M0, M1, soln%Mi(i), Mk, err )

    end do

    call isentropic_relations( soln%Mi, soln%Vi )

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
    fun = ((two/(gamma+1))*(one+half*(gamma-1)*M**2))**((gamma+1)/(gamma-1)) &
            - ((A/areaStar)**2)*M**2
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
    dfun = two*M*( ((two/(gamma+1))*(one+half*(gamma-1)*M**2))**(two/(gamma-1))&
            - (A/areaStar)**2 )
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
