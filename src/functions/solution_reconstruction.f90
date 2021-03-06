module solution_reconstruction

  use set_precision,       only : prec
  use set_constants,       only : small, zero, one, two, fourth
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
    real(prec), dimension(neq,i_low-1:i_high+1) ::    r_plus,    r_minus, &
                                                  psi_plus,  psi_minus, &
                                                 dpsi_plus, dpsi_minus
    !real(prec), dimension(neq,i_low-1:i_high) ::    r_plus,    r_minus, &
    !                                              psi_plus,  psi_minus, &
    !                                             dpsi_plus, dpsi_minus
    real(prec) :: den(neq), e4, km1, kp1
    integer :: i, j, k

    e4  = fourth*epsM
    km1 = one - kappaM
    kp1 = one + kappaM

    r_plus     = zero
    r_minus    = zero
    psi_plus   = zero
    psi_minus  = zero
    dpsi_plus  = zero
    dpsi_minus = zero

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
    ! Extrapolate consecutive variations at endpoints
    r_plus(:,i) = two*r_plus(:,i-1) - r_plus(:,i-2)
    i = i_low-1
    r_minus(:,i) = two*r_minus(:,i+1) - r_minus(:,i+2)

    ! Calculate limiters
    call soln%lim%limiter_fun(r_plus,psi_plus,dpsi_plus)
    call soln%lim%limiter_fun(r_minus,psi_minus,dpsi_minus)


    !do i = i_low-1,i_high+1
    !  write(*,*) i, psi_minus(:,i)
    !end do
    !write(*,*)
    !do i = i_low-1,i_high+1
    !  write(*,*) i, psi_plus(:,i)
    !end do
    !stop


    ! Left states & derivatives
    do  i = i_low,i_high+1
      UL(:,i-1) = soln%U(:,i) + e4*(                                         &
           km1*psi_plus(:,i-1)*( soln%U(:,i) - soln%U(:,i-1) )               &
         + kp1*psi_minus(:,i)*( soln%U(:,i+1) - soln%U(:,i) ) )
      soln%duduL(:,1,i-1) = e4*(                                             &
           km1*( dpsi_plus(:,i-1)*r_plus(:,i-1) - psi_plus(:,i-1) )          &
         - kp1*( dpsi_minus(:,i) ) )
      soln%duduL(:,2,i-1) = one - e4*(                                       &
           km1*( dpsi_plus(:,i-1)*(r_plus(:,i-1)+one) - psi_plus(:,i-1) )    &
         - kp1*( dpsi_minus(:,i)*(r_minus(:,i)+one) - psi_minus(:,i) ) )
      soln%duduL(:,3,i-1) = e4*(                                             &
           km1*( dpsi_plus(:,i-1) )                                          &
         - kp1*( dpsi_minus(:,i)*r_minus(:,i) - psi_minus(:,i) ) )
    end do

    ! Right states & derivatives
    do i = i_low-1,i_high
      UR(:,i) = soln%U(:,i+1) - e4*(                                         &
           km1*psi_minus(:,i+1)*( soln%U(:,i+2) - soln%U(:,i+1) )            &
         + kp1*psi_plus(:,i)*( soln%U(:,i+1) - soln%U(:,i) ) )
      soln%duduR(:,1,i) = e4*(                                               &
           km1*( dpsi_minus(:,i+1) )                                         &
         - kp1*( dpsi_plus(:,i)*r_plus(:,i) - psi_plus(:,i) ) )
      soln%duduR(:,2,i) = one - e4*(                                         &
           km1*( dpsi_minus(:,i+1)*(r_minus(:,i+1)+one) - psi_minus(:,i+1) ) &
         - kp1*( dpsi_plus(:,i)*(r_plus(:,i)+one) - psi_plus(:,i) ) )
      soln%duduR(:,3,i) = e4*(                                               &
           km1*( dpsi_minus(:,i+1)*r_minus(:,i+1) - psi_minus(:,i+1) )       &
         - kp1*( dpsi_plus(:,i) ) )
    end do

  end subroutine MUSCL_extrap


  !============================ MUSCL_extrap_lite ============================80
  !>
  !! Description: MUSCL extrapolation for primitive variables
  !<
  !===========================================================================80
  subroutine MUSCL_extrap_lite(soln,VL,VR)
    type(soln_t), intent(inout) :: soln
    real(prec), intent(out), dimension(neq,i_low-1:i_high) :: VL, VR
    real(prec), dimension(neq,i_low-1:i_high) ::    r_plus,    r_minus, &
                                                  psi_plus,  psi_minus, &
                                                 dpsi_plus, dpsi_minus
    !real(prec), dimension(neq,i_low-1:i_high) ::    r_plus, &
    !                                              psi_plus, &
    !                                             dpsi_plus
    !real(prec), dimension(neq,i_low:i_high+1) ::   r_minus, &
    !                                             psi_minus, &
    !                                            dpsi_minus
    real(prec) :: den(neq), e4, km1, kp1
    integer :: i, j, k

    e4  = fourth*epsM
    km1 = one - kappaM
    kp1 = one + kappaM

    ! Calculate consecutive variations
    i = i_low-1
    den = soln%V(:,i+1) - soln%V(:,i)                      ! u_{i+1} - u_{i}
    den = sign(one,den)*max(abs(den),small)
    r_plus(:,i) = ( soln%V(:,i+2) - soln%V(:,i+1) )/den    ! r+_{i+1/2}

    do  i = i_low,i_high
      den = soln%V(:,i+1) - soln%V(:,i)                    ! u_{i+1} - u_{i}
      den = sign(one,den)*max(abs(den),small)
      r_plus(:,i)  = ( soln%V(:,i+2) - soln%V(:,i+1) )/den ! r+_{i+1/2}
      r_minus(:,i-1) = ( soln%V(:,i) - soln%V(:,i-1) )/den   ! r-_{i+1/2}
    end do

    i = i_high+1
    den = soln%V(:,i+1) - soln%V(:,i)                      ! u_{i+1} - u_{i}
    den = sign(one,den)*max(abs(den),small)
    r_minus(:,i-1) = ( soln%V(:,i) - soln%V(:,i-1) )/den     ! r-_{i+1/2}

    ! Calculate limiters
    call soln%lim%limiter_fun(r_plus,psi_plus,dpsi_plus)
    call soln%lim%limiter_fun(r_minus,psi_minus,dpsi_minus)

    ! Left states & derivatives
    do  i = i_low,i_high+1
      VL(:,i-1) = soln%V(:,i) + e4*(                                         &
           km1*psi_plus(:,i-1)*( soln%V(:,i) - soln%V(:,i-1) )               &
         + kp1*psi_minus(:,i)*( soln%V(:,i+1) - soln%V(:,i) ) )
    end do

    ! Right states & derivatives
    do i = i_low-1,i_high
      VR(:,i) = soln%V(:,i+1) - e4*(                                         &
           km1*psi_minus(:,i+1)*( soln%V(:,i+2) - soln%V(:,i+1) )            &
         + kp1*psi_plus(:,i)*( soln%V(:,i+1) - soln%V(:,i) ) )
    end do

    !do i = i_low-1,i_high
    !  write(*,*) (UL(j,i),j=1,neq), '  ', (UR(j,i),j=1,neq)
    !end do
    !stop

  end subroutine MUSCL_extrap_lite

end module solution_reconstruction
