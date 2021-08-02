module q1d_problem

  use set_precision, only : prec
  use set_constants, only : zero, one, two, half
  use set_inputs, only : neq
  use fluid_constants, only : gamma
  use quadrature, only : gauss_integ
  use variable_conversion, only : speed_of_sound
  use utilities, only : newton_safe, hunt

  implicit none

  private

  public :: solve_Q1D, solve_q1d_isentropic_nozzle

contains

  subroutine solve_Q1D(A,areaStar,M0,bnd_a1,bnd_b1,tol,max_iter,M)
    real(prec), intent(in) :: A, areaStar, M0, bnd_a1, bnd_b1, tol
    real(prec), intent(out) :: M
    real(prec) :: err
    integer, intent(in) :: max_iter

    call newton_safe( f_Q1D, df_Q1D, M0, bnd_a1, bnd_b1, &
                                      tol, max_iter, M, err )
  contains

    function f_Q1D(M)
      use set_inputs
      real(prec) :: f_Q1D
      real(prec), intent(in) :: M

      f_Q1D = ((two/(gamma+one))*(one+half*(gamma-one)*M**2)) &
            **((gamma+one)/(gamma-one)) - ((A/areaStar)**2)*M**2
    end function f_Q1D

    function df_Q1D(M)
      use set_inputs
      real(prec) :: df_Q1D
      real(prec), intent(in) :: M

      df_Q1D = two*M*( ((two/(gamma+one))*(one+half*(gamma-one)*M**2)) &
            **(two/(gamma-one)) - (A/areaStar)**2 )
    end function df_Q1D

  end subroutine solve_Q1D

  subroutine solve_q1d_isentropic_nozzle(x1,x2,xi,w,Mach)
    use set_inputs, only : areaStar, area
    !real(prec), external :: area
    real(prec), intent(in) :: x1, x2
    real(prec), intent(in) :: xi(:), w(:)
    real(prec), intent(out) :: Mach
    integer :: i
    real(prec) :: M0,Ma,Mb

    M0 = 0.1_prec
    Ma = 0.001_prec
    Mb = 1.0_prec
    call gauss_integ(Mach,f,x1,x2,xi,w)
    Mach = Mach/abs(x2-x1)

  contains

    function f(x)
      real(prec) :: f
      real(prec), intent(in) :: x
      real(prec) :: atmp

      atmp = area(x)
      call solve_Q1D(atmp,areaStar,M0,Ma,Mb,1.0e-16_prec,100,f)
      write(*,*) x, atmp, f
    end function f

  end subroutine solve_q1d_isentropic_nozzle

end module q1d_problem
