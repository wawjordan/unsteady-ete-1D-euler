module limiter_calc

  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, four, half, fourth
  use set_inputs, only : i_low, i_high, ig_low, ig_high, imax
  use set_inputs, only : neq, beta_lim, n_ghost
  implicit none

  private

  public :: calculate_limiters
  type limiter_t
    procedure( calc_limiter ), nopass, pointer :: ptr => null()
  contains
    procedure, pass :: limiter_fun
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
    if(associated(this%ptr)) call this%ptr(r,psi,dpsi)
  end subroutine limiter_fun


  !============================== select_limiter =============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine select_limiter()

    use set_inputs, only : limiter_scheme

    type(limiter_t) :: limiter

    select case(limiter_scheme)
      case(0)
        limiter%ptr => null_limiter
      case(1)
        limiter%ptr => van_leer_limiter
      case(2)
        limiter%ptr => van_albada_limiter
      case(3)
        limiter%ptr => minmod_limiter
      case(4)
        limiter%ptr => beta_limiter
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
    psi = sign(one, psi)*max(abs(psi),1.0e-12_prec)
    psi = (r + abs(r))/psi
    psi = half*(one+sign(one,r))*psi

    dpsi = abs(r)*(one + r)**2
    dpsi = sign(one, dpsi)*max(abs(dpsi),1.0e-12_prec)
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
    psi = sign(one, psi)*max(abs(psi),1.0e-12_prec)
    psi = ( (r + one)*r )/psi
    psi = half*(one+sign(one,r))*psi

    dpsi = (one + r**2)**2
    dpsi = sign(one, dpsi)*max(abs(dpsi),1.0e-12_prec)
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

  !========================== calculate_limiter ==============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine calculate_limiters(soln)

    use soln_type, only : soln_t
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq,i_low-1:i_high) :: &
            psi_p, psi_m, dpsi_p, dpsi_m, r_plus, r_minus
    real(prec) :: den(neq)
    integer :: i

    !den = soln%V(:,i_low  :i_high+1)       &
    !       - soln%V(:,i_low-1:I_high)
    !den = sign(one,den)*max(abs(den),1e-6_prec)
    do i = i_low, i_high
      den = soln%V(:,i+1) - soln%V(:,i)
      den = sign(one,den)*max(abs(den),1e-6_prec)
      r_plus(:,i)  = ( soln%V(:,i+2) &
                - soln%V(:,i+1) )/den
      r_minus(:,i) = ( soln%V(:,i) &
                - soln%V(:,i-1) )/den
    end do
    !r_plus  = ( soln%V(:,i_low+1:i_high+2) &
    !             - soln%V(:,i_low  :i_high+1) )/den

    !r_minus = ( soln%V(:,i_low-1:i_high ) &
    !             - soln%V(:,i_low-2:i_high-1) )/den
    call soln%lim%limiter_fun(r_plus,psi_p,dpsi_p)
    call soln%lim%limiter_fun(r_minus,psi_m,dpsi_m)


    soln%psi_p( :,i_low-1:i_high) = psi_p
    soln%psi_m( :,i_low-1:i_high) = psi_m

    do i = 1,n_ghost
    soln%psi_p( :,i_low-1-i) = one!&
         !soln%psi_p(:,i_low-i)
    soln%psi_m( :,i_low-1-i) = one!&
         !soln%psi_m(:,i_low-i)
    soln%psi_p( :,i_high+i) = one!&
         !soln%psi_p(:,i_high-1+i)
    soln%psi_m( :,i_high+i) = one!&
         !soln%psi_m(:,i_high-1+i)
    end do

  end subroutine calculate_limiters

end module limiter_calc
