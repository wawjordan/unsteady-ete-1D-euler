module limiter_type

  use set_precision, only : prec
  use set_constants, only : small, zero, one, two, three, four, half, fourth
  use set_inputs, only : neq, beta_lim, n_ghost
  implicit none

  private

  public :: limiter_t, select_limiter
  type limiter_t
    procedure( calc_limiter ), nopass, pointer :: ptr => null()
  contains
    procedure, pass :: limiter_fun
    procedure, pass :: select_limiter
  end type limiter_t

  abstract interface
  !============================== calc_limiter ===============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  pure subroutine calc_limiter(r,psi,dpsi)

    import :: prec
    real(prec), intent(in) :: r
    real(prec), intent(out) :: psi
    real(prec), intent(out) :: dpsi

  end subroutine calc_limiter

  end interface

contains
  ! type-bound elemental
  elemental subroutine limiter_fun(this,r,psi,dpsi)
    class(limiter_t), intent(in) :: this
    real(prec), intent(in) :: r
    real(prec), intent(out) :: psi
    real(prec), intent(out) :: dpsi
    !if(associated(this%ptr)) call this%ptr(r,psi,dpsi)
    call this%ptr(r,psi,dpsi)
  end subroutine limiter_fun

  !============================== select_limiter =============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine select_limiter(this,limiter_scheme)

    integer, intent(in) :: limiter_scheme
    class(limiter_t), intent(inout) :: this

    select case(limiter_scheme)
      case(0)
        this%ptr => null_limiter
        write(*,*) 'soln%lim%limiter_fun = null_limiter'
      case(1)
        this%ptr => van_leer_limiter
        write(*,*) 'soln%lim%limiter_fun = van_leer_limiter'
      case(2)
        this%ptr => van_albada_limiter
        write(*,*) 'soln%lim%limiter_fun = van_albada_limiter'
      case(3)
        this%ptr => minmod_limiter
        write(*,*) 'soln%lim%limiter_fun = minmod_limiter'
      case(4)
        this%ptr => beta_limiter
        write(*,*) 'soln%lim%limiter_fun = beta_limiter'
      case default
    end select

  end subroutine select_limiter

  !============================= null_limiter ================================80
  !>
  !! Description:
  !!
  !! Inputs:      r   :
  !!
  !! Outputs:     psi :
  !<
  !===========================================================================80
  pure subroutine null_limiter(r,psi,dpsi)
    real(prec), intent(in) :: r
    real(prec), intent(out) :: psi
    real(prec), intent(out) :: dpsi
    real(prec) :: dum

    dum = r ! dummy variable to avoid compiler warning

    psi  = one
    dpsi = zero

  end subroutine null_limiter

  !=========================== van_leer_limiter ==============================80
  !>
  !! Description:
  !!
  !! Inputs:      r   :
  !!
  !! Outputs:     psi :
  !<
  !===========================================================================80
  pure subroutine van_leer_limiter(r,psi,dpsi)
    real(prec), intent(in) :: r
    real(prec), intent(out) :: psi
    real(prec), intent(out) :: dpsi

    psi = one + r
    psi = sign(one, psi)*max(abs(psi),small)
    psi = (r + abs(r))/psi
    psi = half*(one+sign(one,r))*psi

    dpsi = abs(r)*(one + r)**2
    dpsi = sign(one, dpsi)*max(abs(dpsi),small)
    dpsi = (r + abs(r))/dpsi
    dpsi = half*(one+sign(one,r))*dpsi

  end subroutine van_leer_limiter

  !======================== van_albada_limiter ===============================80
  !>
  !! Description:
  !!
  !! Inputs:      r   :
  !!
  !! Outputs:     psi :
  !<
  !===========================================================================80
  pure subroutine van_albada_limiter(r,psi,dpsi)
    real(prec), intent(in) :: r
    real(prec), intent(out) :: psi
    real(prec), intent(out) :: dpsi

    psi = one + r**2
    psi = sign(one, psi)*max(abs(psi),small)
    psi = ( (r + one)*r )/psi
    psi = half*(one+sign(one,r))*psi

    dpsi = (one + r**2)**2
    dpsi = sign(one, dpsi)*max(abs(dpsi),small)
    dpsi = ( one - (r - two)*r )/dpsi
    dpsi = half*(one+sign(one,r))*dpsi

  end subroutine van_albada_limiter

  !============================== minmod_limiter =============================80
  !>
  !! Description:
  !!
  !! Inputs:      r   :
  !!
  !! Outputs:     psi :
  !<
  !===========================================================================80
  pure subroutine minmod_limiter(r,psi,dpsi)
    real(prec), intent(in) :: r
    real(prec), intent(out) :: psi
    real(prec), intent(out) :: dpsi
    real(prec) :: tmp1, tmp2

    tmp1 = half*(one+sign(one,r))
    tmp2 = half*(one+sign(one,one-r))

    psi = min(r,one)
    psi = tmp1*psi

    dpsi = tmp1*tmp2

  end subroutine minmod_limiter

  !============================== beta_limiter ===============================80
  !>
  !! Description:
  !!
  !! Inputs:      r   :
  !!
  !! Outputs:     psi :
  !<
  !===========================================================================80
  pure subroutine beta_limiter(r,psi,dpsi)
    real(prec), intent(in) :: r
    real(prec), intent(out) :: psi
    real(prec), intent(out) :: dpsi
    logical :: tmp1, tmp2

    psi = max(zero, min(beta_lim*r,one), min(r,beta_lim))
    psi = half*(one+sign(one,r))*psi

    tmp1 = (r<one).and.(r>=zero)
    tmp2 = (r>=one).and.(r<beta_lim)

    dpsi = zero + beta_lim*merge(one,zero,tmp1) + merge(one,zero,tmp2)

  end subroutine beta_limiter

end module limiter_type
