module residuals

  use set_precision, only : prec
  use set_constants, only : zero, one, half, third, fourth
  use set_inputs,    only : neq, i_low, i_high, ig_low, ig_high
  use grid_type, only : grid_t
  use soln_type, only : soln_t
  use solution_reconstruction, only : MUSCL_extrap
  use flux_jacobians, only : calc_vl_dfdu

  implicit none

  private

  public :: calc_residual, residual_norms, build_LHS_matrix, build_LHS_test
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

  subroutine mark_diag(M,val)
    real(prec), dimension(:,:), intent(inout) :: M
    real(prec), intent(in) :: val
    integer :: n, i, j

    n = size(M,1)

    do j = 1,n
      do i = 1,n
        if (i == j) then
          M(i,j) = val
        else
          M(i,j) = -val
        end if
      end do
    end do
  end subroutine mark_diag

  subroutine build_LHS_test(soln)
    type(soln_t), intent(inout) :: soln
    integer :: i, j, ii, jj
    !------------------------------ Row 1 -------------------------------------
    i = i_low
      !soln%LHS(:,:,i,3) = soln                       ! dRdu_(i)
      call mark_diag( soln%LHS(:,:,i,3), 3000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,4) = tmp                        ! dRdu_(i+1)
      call mark_diag( soln%LHS(:,:,i,4), 4000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,5) = tmp                        ! dRdu_(i+2)
      call mark_diag( soln%LHS(:,:,i,5), 5000.0_prec + real(i,prec) )
    !------------------------------ Row 2 -------------------------------------
    i = i_low+1
      !soln%LHS(:,:,i,2) = tmp                        ! dRdu_(i-1)
      call mark_diag( soln%LHS(:,:,i,2), 2000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,3) = tmp                        ! dRdu_(i)
      call mark_diag( soln%LHS(:,:,i,3), 3000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,4) = tmp                        ! dRdu_(i+1)
      call mark_diag( soln%LHS(:,:,i,4), 4000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,5) = tmp                        ! dRdu_(i+2)
      call mark_diag( soln%LHS(:,:,i,5), 5000.0_prec + real(i,prec) )
    !------------------------------ Rows 3:N-2 --------------------------------
    do i = i_low+2,i_high-2
      !soln%LHS(:,:,i,1) = tmp                        ! dRdu_(i-2)
      call mark_diag( soln%LHS(:,:,i,1), 1000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,2) = tmp                       ! dRdu_(i-1)
      call mark_diag( soln%LHS(:,:,i,2), 2000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,3) = tmp                        ! dRdu_(i)
      call mark_diag( soln%LHS(:,:,i,3), 3000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,4) = tmp                        ! dRdu_(i+1)
      call mark_diag( soln%LHS(:,:,i,4), 4000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,5) = tmp                        ! dRdu_(i+2)
      call mark_diag( soln%LHS(:,:,i,5), 5000.0_prec + real(i,prec) )
    end do
    !------------------------------ Row N-1 -----------------------------------
    i = i_high-1
      !soln%LHS(:,:,i,1) = tmp                        ! dRdu_(i-2)
      call mark_diag( soln%LHS(:,:,i,1), 1000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,2) = tmp                        ! dRdu_(i-1)
      call mark_diag( soln%LHS(:,:,i,2), 2000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,3) = tmp                        ! dRdu_(i)
      call mark_diag( soln%LHS(:,:,i,3), 3000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,4) = tmp                        ! dRdu_(i+1)
      call mark_diag( soln%LHS(:,:,i,4), 4000.0_prec + real(i,prec) )
    !------------------------------ Row N -------------------------------------
    i = i_high
      !soln%LHS(:,:,i,1) = tmp                        ! dRdu_(i-2)
      call mark_diag( soln%LHS(:,:,i,1), 1000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,2) = tmp                        ! dRdu_(i-1)
      call mark_diag( soln%LHS(:,:,i,2), 2000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,3) = tmp                        ! dRdu_(i)
      call mark_diag( soln%LHS(:,:,i,3), 3000.0_prec + real(i,prec) )
!stop
      do j = 1,5
        do i = i_low,i_high
          do jj = 1,neq
            do ii = 1,neq
              write(*,*) soln%LHS(ii,jj,i,j)
            end do
          end do
        end do
      end do

      stop

  end subroutine build_LHS_test


  subroutine build_LHS_matrix(soln,grid)
    type(soln_t), intent(inout) :: soln
    type(grid_t), intent(in)    :: grid
    real(prec), dimension(neq,i_low-1:i_high) :: UL,UR
    real(prec), dimension(neq) :: FL0, FR0, FL1, FR1, &
                                  DtmpL0, DtmpR0, DtmpL1, DtmpR1
    real(prec), dimension(neq,neq) :: tmp, dfduL0, dfduR0, dfduL1, dfduR1
    integer :: ii, jj, i, j

    call MUSCL_extrap(soln,UL,UR)

    do i = i_low-1,i_high
    write(*,*) i, soln%U(:,i)
    end do
    write(*,*)
    do i = i_low,i_high+1
    write(*,*) i, soln%U(:,i)
    end do
    write(*,*) 'MUSCL:'

    do i = i_low-1,i_high
    write(*,*) i, UL(:,i)
    end do
    write(*,*)
    do i = i_low-1,i_high
    write(*,*) i, UL(:,i)
    end do
    stop

    !==========================================================================
    !------------------------------ Row 1 -------------------------------------
    !==========================================================================
    i = i_low
      call calc_vl_dfdu(UL(:,i-1),FL0,dfduL0,.true.) ! dfdu_L_(i-1/2)
      call calc_vl_dfdu(UR(:,i-1),FR0,dfduR0,.false.)! dfdu_R_(i-1/2)
      call calc_vl_dfdu(UL(:,i),FL1,dfduL1,.true.)   ! dfdu_L_(i+1/2)
      call calc_vl_dfdu(UR(:,i),FR1,dfduR1,.false.)  ! dfdu_R_(i+1/2)

      soln%F(:,i-1) = FL0 + FR0                      ! f_(i-1/2)
      soln%F(:,i)   = FL1 + FR1                      ! f_(i+1/2)

      DtmpL0 = grid%Ai(i-1)*soln%duduL(:,3,i-1)      ! A_(i-1/2)*duLdu_(i)
      DtmpR0 = grid%Ai(i-1)*soln%duduR(:,2,i-1)      ! A_(i-1/2)*duRdu_(i)

      DtmpL1 = grid%Ai(i)*soln%duduL(:,2,i)          ! A_(i+1/2)*duLdu_(i)
      DtmpR1 = grid%Ai(i)*soln%duduR(:,1,i)          ! A_(i+1/2)*duRdu_(i)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduL1(ii,jj)*DtmpL1(ii) + dfduR1(ii,jj)*DtmpR1(ii) &
                   - (dfduL0(ii,jj)*DtmpL0(ii) + dfduR0(ii,jj)*DtmpR0(ii))
        end do
      end do
      soln%LHS(:,:,i,3) = tmp                        ! dRdu_(i)
      do ii = 1,neq
        write(*,*) (tmp(ii,jj),jj=1,3)
      end do
      stop

      DtmpR0 = grid%Ai(i-1)*soln%duduR(:,3,i-1)      ! A_(i-1/2)*duRdu_(i+1)

      DtmpL1 = grid%Ai(i)*soln%duduL(:,3,i)          ! A_(i+1/2)*duLdu_(i+1)
      DtmpR1 = grid%Ai(i)*soln%duduR(:,2,i)          ! A_(i+1/2)*duRdu_(i+1)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduL1(ii,jj)*DtmpL1(ii) + dfduR1(ii,jj)*DtmpR1(ii) &
                     - dfduR0(ii,jj)*DtmpR0(ii)
        end do
      end do
      soln%LHS(:,:,i,4) = tmp                        ! dRdu_(i+1)

      DtmpR1 = grid%Ai(i)*soln%duduR(:,3,i)          ! A_(i+1/2)*duRdu_(i+2)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduR1(ii,jj)*DtmpR1(ii)
        end do
      end do
      soln%LHS(:,:,i,5) = tmp                        ! dRdu_(i+2)

    !==========================================================================
    !------------------------------ Row 2 -------------------------------------
    !==========================================================================
    i = i_low+1
      FL0 = FL1
      FR0 = FR1
      dfduL0 = dfduL1                                ! dfdu_L_(i-1/2)
      dfduR0 = dfduR1                                ! dfdu_R_(i-1/2)
      call calc_vl_dfdu(UL(:,i),FL1,dfduL1,.true.)   ! dfdu_L_(i+1/2)
      call calc_vl_dfdu(UR(:,i),FR1,dfduR1,.false.)  ! dfdu_R_(i+1/2)

      soln%F(:,i-1) = FL0 + FR0                      ! f_(i-1/2)
      soln%F(:,i)   = FL1 + FR1                      ! f_(i+1/2)

      DtmpL0 = grid%Ai(i-1)*soln%duduL(:,2,i-1)      ! A_(i-1/2)*duLdu_(i-1)
      DtmpR0 = grid%Ai(i-1)*soln%duduR(:,1,i-1)      ! A_(i-1/2)*duRdu_(i-1)

      DtmpL1 = grid%Ai(i)*soln%duduL(:,1,i)          ! A_(i+1/2)*duLdu_(i-1)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduL1(ii,jj)*DtmpL1(ii) &
                   - (dfduL0(ii,jj)*DtmpL0(ii) + dfduR0(ii,jj)*DtmpR0(ii))
        end do
      end do
      soln%LHS(:,:,i,2) = tmp                        ! dRdu_(i-1)

      DtmpL0 = grid%Ai(i-1)*soln%duduL(:,3,i-1)      ! A_(i-1/2)*duLdu_(i)
      DtmpR0 = grid%Ai(i-1)*soln%duduR(:,2,i-1)      ! A_(i-1/2)*duRdu_(i)

      DtmpL1 = grid%Ai(i)*soln%duduL(:,2,i)          ! A_(i+1/2)*duLdu_(i)
      DtmpR1 = grid%Ai(i)*soln%duduR(:,1,i)          ! A_(i+1/2)*duRdu_(i)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduL1(ii,jj)*DtmpL1(ii) + dfduR1(ii,jj)*DtmpR1(ii) &
                   - (dfduL0(ii,jj)*DtmpL0(ii) + dfduR0(ii,jj)*DtmpR0(ii))
        end do
      end do
      soln%LHS(:,:,i,3) = tmp                        ! dRdu_(i)

      DtmpR0 = grid%Ai(i-1)*soln%duduR(:,3,i-1)      ! A_(i-1/2)*duRdu_(i+1)

      DtmpL1 = grid%Ai(i)*soln%duduL(:,3,i)          ! A_(i+1/2)*duLdu_(i+1)
      DtmpR1 = grid%Ai(i)*soln%duduR(:,2,i)          ! A_(i+1/2)*duRdu_(i+1)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduL1(ii,jj)*DtmpL1(ii) + dfduR1(ii,jj)*DtmpR1(ii) &
                     - dfduR0(ii,jj)*DtmpR0(ii)
        end do
      end do
      soln%LHS(:,:,i,4) = tmp                        ! dRdu_(i+1)

      DtmpR1 = grid%Ai(i)*soln%duduR(:,3,i)          ! A_(i+1/2)*duRdu_(i+2)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduR1(ii,jj)*DtmpR1(ii)
        end do
      end do
      soln%LHS(:,:,i,5) = tmp                        ! dRdu_(i+2)

    !==========================================================================
    !------------------------------ Rows 3:N-2 --------------------------------
    !==========================================================================
    do i = i_low+2,i_high-2
      FL0 = FL1
      FR0 = FR1
      dfduL0 = dfduL1                                ! dfdu_L_(i-1/2)
      dfduR0 = dfduR1                                ! dfdu_R_(i-1/2)
      call calc_vl_dfdu(UL(:,i),FL1,dfduL1,.true.)   ! dfdu_L_(i+1/2)
      call calc_vl_dfdu(UR(:,i),FR1,dfduR1,.false.)  ! dfdu_R_(i+1/2)

      DtmpL0 = grid%Ai(i-1)*soln%duduL(:,1,i-1)      ! A_(i-1/2)*duLdu_(i-2)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = - dfduL0(ii,jj)*DtmpL0(ii)
        end do
      end do
      soln%LHS(:,:,i,1) = tmp                        ! dRdu_(i-2)

      DtmpL0 = grid%Ai(i-1)*soln%duduL(:,2,i-1)      ! A_(i-1/2)*duLdu_(i-1)
      DtmpR0 = grid%Ai(i-1)*soln%duduR(:,1,i-1)      ! A_(i-1/2)*duRdu_(i-1)

      DtmpL1 = grid%Ai(i)*soln%duduL(:,1,i)          ! A_(i+1/2)*duLdu_(i-1)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduL1(ii,jj)*DtmpL1(ii) &
                   - (dfduL0(ii,jj)*DtmpL0(ii) + dfduR0(ii,jj)*DtmpR0(ii))
        end do
      end do
      soln%LHS(:,:,i,2) = tmp                        ! dRdu_(i-1)


      DtmpL0 = grid%Ai(i-1)*soln%duduL(:,3,i-1)      ! A_(i-1/2)*duLdu_(i)
      DtmpR0 = grid%Ai(i-1)*soln%duduR(:,2,i-1)      ! A_(i-1/2)*duRdu_(i)

      DtmpL1 = grid%Ai(i)*soln%duduL(:,2,i)          ! A_(i+1/2)*duLdu_(i)
      DtmpR1 = grid%Ai(i)*soln%duduR(:,1,i)          ! A_(i+1/2)*duRdu_(i)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduL1(ii,jj)*DtmpL1(ii) + dfduR1(ii,jj)*DtmpR1(ii) &
                   - (dfduL0(ii,jj)*DtmpL0(ii) + dfduR0(ii,jj)*DtmpR0(ii))
        end do
      end do
      soln%LHS(:,:,i,3) = tmp                        ! dRdu_(i)


      DtmpR0 = grid%Ai(i-1)*soln%duduR(:,3,i-1)      ! A_(i-1/2)*duRdu_(i+1)

      DtmpL1 = grid%Ai(i)*soln%duduL(:,3,i)          ! A_(i+1/2)*duLdu_(i+1)
      DtmpR1 = grid%Ai(i)*soln%duduR(:,2,i)          ! A_(i+1/2)*duRdu_(i+1)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduL1(ii,jj)*DtmpL1(ii) + dfduR1(ii,jj)*DtmpR1(ii) &
                     - dfduR0(ii,jj)*DtmpR0(ii)
        end do
      end do
      soln%LHS(:,:,i,4) = tmp                        ! dRdu_(i+1)

      DtmpR1 = grid%Ai(i)*soln%duduR(:,3,i)          ! A_(i+1/2)*duRdu_(i+2)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduR1(ii,jj)*DtmpR1(ii)
        end do
      end do
      soln%LHS(:,:,i,5) = tmp                        ! dRdu_(i+2)
    end do

    !==========================================================================
    !------------------------------ Row N-1 -----------------------------------
    !==========================================================================
    i = i_high-1
      FL0 = FL1
      FR0 = FR1
      dfduL0 = dfduL1                                ! dfdu_L_(i-1/2)
      dfduR0 = dfduR1                                ! dfdu_R_(i-1/2)
      call calc_vl_dfdu(UL(:,i),FL1,dfduL1,.true.)   ! dfdu_L_(i+1/2)
      call calc_vl_dfdu(UR(:,i),FR1,dfduR1,.false.)  ! dfdu_R_(i+1/2)

      DtmpL0 = grid%Ai(i-1)*soln%duduL(:,1,i-1)      ! A_(i-1/2)*duLdu_(i-2)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = - dfduL0(ii,jj)*DtmpL0(ii)
        end do
      end do
      soln%LHS(:,:,i,1) = tmp                        ! dRdu_(i-2)

      DtmpL0 = grid%Ai(i-1)*soln%duduL(:,2,i-1)      ! A_(i-1/2)*duLdu_(i-1)
      DtmpR0 = grid%Ai(i-1)*soln%duduR(:,1,i-1)      ! A_(i-1/2)*duRdu_(i-1)

      DtmpL1 = grid%Ai(i)*soln%duduL(:,1,i)          ! A_(i+1/2)*duLdu_(i-1)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduL1(ii,jj)*DtmpL1(ii) &
                   - (dfduL0(ii,jj)*DtmpL0(ii) + dfduR0(ii,jj)*DtmpR0(ii))
        end do
      end do
      soln%LHS(:,:,i,2) = tmp                        ! dRdu_(i-1)


      DtmpL0 = grid%Ai(i-1)*soln%duduL(:,3,i-1)      ! A_(i-1/2)*duLdu_(i)
      DtmpR0 = grid%Ai(i-1)*soln%duduR(:,2,i-1)      ! A_(i-1/2)*duRdu_(i)

      DtmpL1 = grid%Ai(i)*soln%duduL(:,2,i)          ! A_(i+1/2)*duLdu_(i)
      DtmpR1 = grid%Ai(i)*soln%duduR(:,1,i)          ! A_(i+1/2)*duRdu_(i)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduL1(ii,jj)*DtmpL1(ii) + dfduR1(ii,jj)*DtmpR1(ii) &
                   - (dfduL0(ii,jj)*DtmpL0(ii) + dfduR0(ii,jj)*DtmpR0(ii))
        end do
      end do
      soln%LHS(:,:,i,3) = tmp                        ! dRdu_(i)


      DtmpR0 = grid%Ai(i-1)*soln%duduR(:,3,i-1)      ! A_(i-1/2)*duRdu_(i+1)

      DtmpL1 = grid%Ai(i)*soln%duduL(:,3,i)          ! A_(i+1/2)*duLdu_(i+1)
      DtmpR1 = grid%Ai(i)*soln%duduR(:,2,i)          ! A_(i+1/2)*duRdu_(i+1)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduL1(ii,jj)*DtmpL1(ii) + dfduR1(ii,jj)*DtmpR1(ii) &
                     - dfduR0(ii,jj)*DtmpR0(ii)
        end do
      end do
      soln%LHS(:,:,i,4) = tmp                        ! dRdu_(i+1)

    !==========================================================================
    !------------------------------ Row N -------------------------------------
    !==========================================================================
    i = i_high
      FL0 = FL1
      FR0 = FR1
      dfduL0 = dfduL1                                ! dfdu_L_(i-1/2)
      dfduR0 = dfduR1                                ! dfdu_R_(i-1/2)
      call calc_vl_dfdu(UL(:,i),FL1,dfduL1,.true.)   ! dfdu_L_(i+1/2)
      call calc_vl_dfdu(UR(:,i),FR1,dfduR1,.false.)  ! dfdu_R_(i+1/2)

      DtmpL0 = grid%Ai(i-1)*soln%duduL(:,1,i-1)      ! A_(i-1/2)*duLdu_(i-2)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = - dfduL0(ii,jj)*DtmpL0(ii)
        end do
      end do
      soln%LHS(:,:,i,1) = tmp                        ! dRdu_(i-2)

      DtmpL0 = grid%Ai(i-1)*soln%duduL(:,2,i-1)      ! A_(i-1/2)*duLdu_(i-1)
      DtmpR0 = grid%Ai(i-1)*soln%duduR(:,1,i-1)      ! A_(i-1/2)*duRdu_(i-1)

      DtmpL1 = grid%Ai(i)*soln%duduL(:,1,i)          ! A_(i+1/2)*duLdu_(i-1)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduL1(ii,jj)*DtmpL1(ii) &
                   - (dfduL0(ii,jj)*DtmpL0(ii) + dfduR0(ii,jj)*DtmpR0(ii))
        end do
      end do
      soln%LHS(:,:,i,2) = tmp                        ! dRdu_(i-1)


      DtmpL0 = grid%Ai(i-1)*soln%duduL(:,3,i-1)      ! A_(i-1/2)*duLdu_(i)
      DtmpR0 = grid%Ai(i-1)*soln%duduR(:,2,i-1)      ! A_(i-1/2)*duRdu_(i)

      DtmpL1 = grid%Ai(i)*soln%duduL(:,2,i)          ! A_(i+1/2)*duLdu_(i)
      DtmpR1 = grid%Ai(i)*soln%duduR(:,1,i)          ! A_(i+1/2)*duRdu_(i)

      do jj = 1,3
        do ii = 1,neq
          tmp(ii,jj) = dfduL1(ii,jj)*DtmpL1(ii) + dfduR1(ii,jj)*DtmpR1(ii) &
                   - (dfduL0(ii,jj)*DtmpL0(ii) + dfduR0(ii,jj)*DtmpR0(ii))
        end do
      end do
      soln%LHS(:,:,i,3) = tmp                        ! dRdu_(i)
!stop
      do j = 1,5
        do i = i_low,i_high
          do jj = 1,neq
            do ii = 1,neq
              write(*,*) soln%LHS(ii,jj,i,j)
            end do
          end do
        end do
      end do

      stop

  end subroutine build_LHS_matrix

end module residuals
