module residuals

  use set_precision, only : prec
  use set_constants, only : zero, one, half, third, fourth
  use set_inputs,    only : neq, i_low, i_high, ig_low, ig_high
  use grid_type, only : grid_t
  use soln_type, only : soln_t
  use variable_conversion

  implicit none

  private

  public :: calc_residual, residual_norms, build_LHS_matrix
  contains

  !============================== calc_residual ==============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine calc_residual(grid,soln)

    type(grid_t), intent(inout) :: grid
    type(soln_t), intent(inout) :: soln
    integer :: i

    do i = i_low,i_high
      soln%R(:,i) = grid%Ai(i)*soln%F(:,i) &
                  - grid%Ai(i-1)*soln%F(:,i-1) &
                  - grid%Ac(i)*grid%dx(i)*soln%S(:,i)
    end do
  end subroutine calc_residual

  !============================= residual_norms ==============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine residual_norms( R, Rnorm, rinit )

    real(prec), dimension(neq,i_low:i_high), intent(in) :: R
    real(prec), dimension(neq), intent(in)  :: rinit
    real(prec), dimension(neq), intent(out) :: Rnorm
    real(prec) :: Linv
    integer :: i

    Linv = one/real(size(R,2),prec)

    do i = 1,neq
      Rnorm(i) = sqrt(Linv*sum(R(i,:)**2))
    end do
    Rnorm = Rnorm/rinit

  end subroutine residual_norms

  subroutine build_LHS_matrix(soln,grid)
    type(soln_t), intent(inout) :: soln
    type(grid_t), intent(in)    :: grid
    real(prec), dimension(neq,i_low-1:i_high) :: Left, Right
    real(prec), dimension(neq,neq) :: dfduim1, dfdui, dfduip1, dfduip2
    real(prec), dimension(neq,neq) :: dgduim1, dgdui, dgduip1, dgduip2
    real(prec) :: vt
    integer :: i, j

    !call MUSCL_extrap(soln,Left,Right)
    ! left boundary -> calculate dfdu_{i-1/2} wrt i, i+1, & i+2
!    i = i_low-1
!    call VL_flux_jac( Left(:,i), Right(:,i), &
!            soln%U(:,i-1), soln%U(:,i), soln%U(:,i+1), soln%U(:,i+2), &
!            soln%psi_p(:,i-1), soln%psi_p(:,i), &
!            soln%psi_m(:,i), soln%psi_m(:,i+1), &
!            soln%dpsi_p(:,i-1), soln%dpsi_p(:,i), &
!            soln%dpsi_m(:,i), soln%dpsi_m(:,i+1), &
!            dfduim1, dfdui, dfduip1, dfduip2 )
!    soln%LHS(:,:,i+1,3) = -grid%Ai(i)*dfduim1    !A*df{i-1/2}/du{1}
!    soln%LHS(:,:,i+1,4) = -grid%Ai(i)*dfdui  !A*df{i-1/2}/du{2}
!    soln%LHS(:,:,i+1,5) = -grid%Ai(i)*dfduip1  !A*df{i-1/2}/du{3}
!    do i = i_low,i_high-1
!      !write(*,*) soln%psi_p(1,i), soln%psi_m(1,i)
!      call VL_flux_jac( Left(:,i), Right(:,i), &
!              soln%U(:,i-1), soln%U(:,i), soln%U(:,i+1), soln%U(:,i+2), &
!              soln%psi_p(:,i-1), soln%psi_p(:,i), &
!              soln%psi_m(:,i), soln%psi_m(:,i+1), &
!              soln%dpsi_p(:,i-1), soln%dpsi_p(:,i), &
!              soln%dpsi_m(:,i), soln%dpsi_m(:,i+1), &
!              dfduim1, dfdui, dfduip1, dfduip2 )
!      do j = 1,neq
!        write(*,*) dfdui(j,:)
!      end do
!      write(*,*)
!      vt = grid%Ac(i)*grid%dx(i)/soln%dt(i)
!      soln%LHS(:,:,i,1) = soln%LHS(:,:,i,1) + grid%Ai(i)*dfduim1
!      soln%LHS(:,:,i,2) = soln%LHS(:,:,i,2) + grid%Ai(i)*dfdui
!      soln%LHS(:,:,i,3) = vt + soln%LHS(:,:,i,3) + grid%Ai(i)*dfduip1
!      soln%LHS(:,:,i,4) = soln%LHS(:,:,i,4) + grid%Ai(i)*dfduip2
!
!      soln%LHS(:,:,i+1,2) = -grid%Ai(i)*dfduim1
!      soln%LHS(:,:,i+1,3) = -grid%Ai(i)*dfdui
!      soln%LHS(:,:,i+1,4) = -grid%Ai(i)*dfduip1
!      soln%LHS(:,:,i+1,5) = -grid%Ai(i)*dfduip2
!    end do
!    i = i_high+1
!    call VL_flux_jac( Left(:,i), Right(:,i), &
!            soln%U(:,i-1), soln%U(:,i), soln%U(:,i+1), soln%U(:,i+2), &
!            soln%psi_p(:,i-1), soln%psi_p(:,i), &
!            soln%psi_m(:,i), soln%psi_m(:,i+1), &
!            soln%dpsi_p(:,i-1), soln%dpsi_p(:,i), &
!            soln%dpsi_m(:,i), soln%dpsi_m(:,i+1), &
!            dfduim1, dfdui, dfduip1, dfduip2 )
!      soln%LHS(:,:,i,1) = soln%LHS(:,:,i,1) + grid%Ai(i)*dfdui
!      soln%LHS(:,:,i,2) = soln%LHS(:,:,i,2) + grid%Ai(i)*dfduip1
!      soln%LHS(:,:,i,3) = soln%LHS(:,:,i,3) + grid%Ai(i)*dfduip2

  end subroutine build_LHS_matrix

end module residuals
