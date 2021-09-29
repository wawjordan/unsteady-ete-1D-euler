module solution_reconstruction

  use set_precision,       only : prec
  use set_constants,       only : small, zero, one, fourth
  use set_inputs,          only : neq, i_low, i_high, ig_low, ig_high, &
                                  epsM, kappaM
  use fluid_constants,     only : gamma
  use soln_type,           only : soln_t
  use grid_type,           only : grid_t

  implicit none

  private

  public :: MUSCL_extrap, MUSCL_extrap_lite

  contains

  !================================ MUSCL_extrap =============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine MUSCL_extrap(soln,UL,UR)
    type(soln_t), intent(inout) :: soln
    real(prec), intent(out), dimension(neq,i_low-1:i_high) :: UL, UR
    real(prec), dimension(neq,i_low-1:i_high) ::    r_plus,    r_minus, &
                                                  psi_plus,  psi_minus, &
                                                 dpsi_plus, dpsi_minus
    real(prec) :: den(neq), e4, km1, kp1
    integer :: i, j, k

    e4  = fourth*epsM
    km1 = one - kappaM
    kp1 = one + kappaM

    ! Calculate consecutive variations
    i = i_low-1
    den = soln%U(:,i+1) - soln%U(:,i)
    den = sign(one,den)*max(abs(den),small)
    r_plus(:,i) = ( soln%U(:,i+2) - soln%U(:,i+1) )/den

    do  i = i_low,i_high
      den = soln%U(:,i+1) - soln%U(:,i)
      den = sign(one,den)*max(abs(den),small)
      r_plus(:,i)  = ( soln%U(:,i+2) - soln%U(:,i+1) )/den
      r_minus(:,i) = ( soln%U(:,i) - soln%U(:,i-1) )/den
    end do

    i = i_high+1
    den = soln%U(:,i+1) - soln%U(:,i)
    den = sign(one,den)*max(abs(den),small)
    r_minus(:,i) = ( soln%U(:,i) - soln%U(:,i-1) )/den

    ! Calculate limiters
    call soln%lim%limiter_fun(r_plus,psi_plus,dpsi_plus)
    call soln%lim%limiter_fun(r_minus,psi_minus,dpsi_minus)

    ! Left states & derivatives
    do  i = i_low,i_high+1
      UL(:,i-1) = soln%U(:,i) + e4*(                                         &
           km1*psi_plus(:,i-1)*( soln%U(:,i) - soln%U(:,i-1) )               &
         + kp1*psi_minus(:,i)*( soln%U(:,i+1) - soln%U(:,i) ) )
      soln%duduL(1,:,i-1) = e4*(                                             &
           km1*( dpsi_plus(:,i-1)*r_plus(:,i-1) - psi_plus(:,i-1) )          &
         - kp1*( dpsi_minus(:,i) ) )
      soln%duduL(2,:,i-1) = one - e4*(                                       &
           km1*( dpsi_plus(:,i-1)*(r_plus(:,i-1)+one) - psi_plus(:,i-1) )    &
         - kp1*( dpsi_minus(:,i)*(r_minus(:,i)+one) - psi_minus(:,i) ) )
      soln%duduL(3,:,i-1) = e4*(                                             &
           km1*( dpsi_plus(:,i-1) )                                          &
         - kp1*( dpsi_minus(:,i)*r_minus(:,i) - psi_minus(:,i) ) )
    end do

    ! Right states & derivatives
    do i = i_low-1,i_high
      UR(:,i) = soln%U(:,i+1) - e4*(                                         &
           km1*psi_minus(:,i+1)*( soln%U(:,i+2) - soln%U(:,i+1) )            &
         + kp1*psi_plus(:,i)*( soln%U(:,i+1) - soln%U(:,i) ) )
      soln%duduR(1,:,i) = e4*(                                               &
           km1*( dpsi_minus(:,i+1) )                                         &
         - kp1*( dpsi_plus(:,i)*r_plus(:,i) - psi_plus(:,i) ) )
      soln%duduR(2,:,i) = one - e4*(                                         &
           km1*( dpsi_minus(:,i+1)*(r_minus(:,i+1)+one) - psi_minus(:,i+1) ) &
         - kp1*( dpsi_plus(:,i)*(r_plus(:,i)+one) - psi_plus(:,i) ) )
      soln%duduR(3,:,i) = e4*(                                               &
           km1*( dpsi_minus(:,i+1)*r_minus(:,i+1) - psi_minus(:,i+1) )       &
         - kp1*( dpsi_plus(:,i) ) )
    end do

  end subroutine MUSCL_extrap


  !============================ MUSCL_extrap_lite ============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine MUSCL_extrap_lite(soln,UL,UR)
    type(soln_t), intent(inout) :: soln
    real(prec), intent(out), dimension(neq,i_low-1:i_high) :: UL, UR
    real(prec), dimension(neq,i_low-1:i_high) ::    r_plus,    r_minus, &
                                                  psi_plus,  psi_minus, &
                                                 dpsi_plus, dpsi_minus
    real(prec) :: den(neq), e4, km1, kp1
    integer :: i, j, k

    e4  = fourth*epsM
    km1 = one - kappaM
    kp1 = one + kappaM

    ! Calculate consecutive variations
    i = i_low-1
    den = soln%U(:,i+1) - soln%U(:,i)
    den = sign(one,den)*max(abs(den),small)
    r_plus(:,i) = ( soln%U(:,i+2) - soln%U(:,i+1) )/den

    do  i = i_low,i_high
      den = soln%U(:,i+1) - soln%U(:,i)
      den = sign(one,den)*max(abs(den),small)
      r_plus(:,i)  = ( soln%U(:,i+2) - soln%U(:,i+1) )/den
      r_minus(:,i) = ( soln%U(:,i) - soln%U(:,i-1) )/den
    end do

    i = i_high+1
    den = soln%U(:,i+1) - soln%U(:,i)
    den = sign(one,den)*max(abs(den),small)
    r_minus(:,i) = ( soln%U(:,i) - soln%U(:,i-1) )/den

    ! Calculate limiters
    call soln%lim%limiter_fun(r_plus,psi_plus,dpsi_plus)
    call soln%lim%limiter_fun(r_minus,psi_minus,dpsi_minus)

    ! Left states & derivatives
    do  i = i_low,i_high+1
      UL(:,i-1) = soln%U(:,i) + e4*(                                         &
           km1*psi_plus(:,i-1)*( soln%U(:,i) - soln%U(:,i-1) )               &
         + kp1*psi_minus(:,i)*( soln%U(:,i+1) - soln%U(:,i) ) )
    end do

    ! Right states & derivatives
    do i = i_low-1,i_high
      UR(:,i) = soln%U(:,i+1) - e4*(                                         &
           km1*psi_minus(:,i+1)*( soln%U(:,i+2) - soln%U(:,i+1) )            &
         + kp1*psi_plus(:,i)*( soln%U(:,i+1) - soln%U(:,i) ) )
    end do

  end subroutine MUSCL_extrap_lite

end module solution_reconstruction
