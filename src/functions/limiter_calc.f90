module limiter_calc

  use set_precision, only : prec
  use set_constants, only : small, zero, one
  use set_inputs, only : i_low, i_high, neq, n_ghost
  use soln_type, only : soln_t

  implicit none

  private

  !public :: calculate_limiters

contains


  elemental subroutine drplus(um1,ui,up1,dm1,di,dp1)

    real(prec), intent(in)  :: um1, ui, up1
    real(prec), intent(out) :: dm1, di, dp1
    real(prec) :: den

    den = ui - um1
    den = sign(one,den)*max(abs(den),small)

    dm1 =  ( up1 - ui  )/( den**2 )
    di  = -( up1 - um1 )/( den**2 )
    dp1 =              one/den

  end subroutine drplus

  elemental subroutine drminus(um1,ui,up1,dm1,di,dp1)

    real(prec), intent(in)  :: um1, ui, up1
    real(prec), intent(out) :: dm1, di, dp1
    real(prec) :: den

    den = up1 - ui
    den = sign(one,den)*max(abs(den),small)

    dm1 =             -one/den
    di =   ( up1 - um1 )/( den**2 )
    dp1 = -( ui - um1  )/( den**2 )

  end subroutine drminus

  !========================== calculate_limiter ==============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine calculate_limiters(soln)

    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq,i_low-1:i_high) :: &
            psi_p, psi_m, dpsi_p, dpsi_m, r_plus, r_minus
    real(prec) :: den(neq)
    integer :: i

    !den = soln%V(:,i_low  :i_high+1)       &
    !       - soln%V(:,i_low-1:I_high)
    !den = sign(one,den)*max(abs(den),1e-6_prec)
    !r_plus  = ( soln%V(:,i_low+1:i_high+2) &
    !             - soln%V(:,i_low  :i_high+1) )/den

    !r_minus = ( soln%V(:,i_low-1:i_high ) &
    !             - soln%V(:,i_low-2:i_high-1) )/den
    do i = i_low, i_high
      den = soln%V(:,i+1) - soln%V(:,i)
      den = sign(one,den)*max(abs(den),small)
      r_plus(:,i)  = ( soln%V(:,i+2) &
                - soln%V(:,i+1) )/den
      r_minus(:,i) = ( soln%V(:,i) &
                - soln%V(:,i-1) )/den
    end do

    soln%dpsi_p = one
    soln%dpsi_m = one
    soln%psi_p = zero
    soln%psi_m = zero


    call soln%lim%limiter_fun(r_plus,psi_p,dpsi_p)
    call soln%lim%limiter_fun(r_minus,psi_m,dpsi_m)

    soln%psi_p( :,i_low-1:i_high)  = psi_p
    soln%psi_m( :,i_low-1:i_high)  = psi_m
    soln%dpsi_p( :,i_low-1:i_high) = dpsi_p
    soln%dpsi_m( :,i_low-1:i_high) = dpsi_m

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
