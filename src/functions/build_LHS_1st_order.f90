module build_LHS_1st_order
  use dispmodule
  use set_precision, only : prec
  use set_constants, only : zero, one, half, third, fourth
  use set_inputs,    only : neq, i_low, i_high, ig_low, ig_high
  use grid_type, only : grid_t
  use soln_type, only : soln_t
  use solution_reconstruction, only : MUSCL_extrap
  use flux_jacobians, only : calc_vl_dfdu

  implicit none

  private

  public :: build_LHS_matrix_1st_order, build_tridiag_test
  contains

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

  subroutine build_tridiag_test(soln)
    type(soln_t), intent(inout) :: soln
    integer :: i, j, ii, jj
    !------------------------------ Row 1 -------------------------------------
    i = i_low
      !soln%LHS(:,:,i,2) = soln                       ! dRdu_(i)
      call mark_diag( soln%LHS(:,:,i,2), 2000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,3) = tmp                        ! dRdu_(i+1)
      call mark_diag( soln%LHS(:,:,i,3), 3000.0_prec + real(i,prec) )
    !------------------------------ Rows 2:N-1 --------------------------------
    do i = i_low+1,i_high-1
      !soln%LHS(:,:,i,1) = tmp                        ! dRdu_(i-1)
      call mark_diag( soln%LHS(:,:,i,1), 1000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,3) = tmp                        ! dRdu_(i)
      call mark_diag( soln%LHS(:,:,i,2), 2000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,4) = tmp                        ! dRdu_(i+1)
      call mark_diag( soln%LHS(:,:,i,3), 3000.0_prec + real(i,prec) )
    end do
    !------------------------------ Row N -------------------------------------
    i = i_high
      !soln%LHS(:,:,i,1) = tmp                        ! dRdu_(i-2)
      call mark_diag( soln%LHS(:,:,i,1), 1000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,2) = tmp                        ! dRdu_(i-1)
      call mark_diag( soln%LHS(:,:,i,2), 2000.0_prec + real(i,prec) )
      !soln%LHS(:,:,i,3) = tmp                        ! dRdu_(i)
      call mark_diag( soln%LHS(:,:,i,3), 3000.0_prec + real(i,prec) )
!stop
      do j = 1,3
        do i = i_low,i_high
          do jj = 1,neq
            do ii = 1,neq
              write(*,*) soln%LHS(ii,jj,i,j)
            end do
          end do
        end do
      end do

      stop

  end subroutine build_tridiag_test


  subroutine build_LHS_matrix_1st_order(soln,grid)
    type(soln_t), intent(inout) :: soln
    type(grid_t), intent(in)    :: grid
    real(prec), dimension(neq) :: FLm,      &
                                  FLi, FRi, &
                                       FRp, &
                                  DtmpLm,   &
                                  DtmpLi, DtmpRi, &
                                          DtmpRp
    real(prec), dimension(neq,neq) :: tmp, &
                                  dfduLm,         &
                                  dfduLi, dfduRi, &
                                          dfduRp
    integer :: ii, jj, i, j
    ! m --> i-1 , i --> i  , p --> i+1
    ! L --> +   , R --> -

      !do i = i_low,i_high
      !  do j = 0,1
      !    call disp('    ',soln%LHS(:,:,i,j), advance='NO')
      !  end do
      ! call disp('    ',soln%LHS(:,:,i,j+1))
      !end do
      !stop

    !==========================================================================
    !------------------------------ Row 1 -------------------------------------
    !==========================================================================
    i = i_low
      call calc_vl_dfdu(soln%U(:,i-1),FLm,dfduLm,.true. ) ! dfdu_+_(i-1)
      call calc_vl_dfdu(soln%U(:,i  ),FLi,dfduLi,.true. ) ! dfdu_+_(i)
      call calc_vl_dfdu(soln%U(:,i  ),FRi,dfduRi,.false.) ! dfdu_-_(i)
      call calc_vl_dfdu(soln%U(:,i+1),FRp,dfduRp,.false.) ! dfdu_-_(i+1)
      soln%F(:,i-1) = FLm + FRi                      ! f_(i-1/2)
      soln%F(:,i)   = FLi + FRp                      ! f_(i+1/2)

      DtmpLm = grid%Ai(i-1)        ! A_(i-1/2)*duLdu_(i-1)
      DtmpRi = grid%Ai(i-1)        ! A_(i-1/2)*duRdu_(i-1)
      DtmpLi = grid%Ai(i)          ! A_(i+1/2)*duLdu_(i)
      DtmpRp = grid%Ai(i)          ! A_(i+1/2)*duRdu_(i)

      do jj = 1,neq
        do ii = 1,neq
          tmp(ii,jj) = dfduLi(ii,jj)*DtmpLi(ii) - dfduRi(ii,jj)*DtmpRi(ii)
        end do
      end do
      soln%LHS(:,:,i,2) = tmp                        ! dRdu_(i)

      do jj = 1,neq
        do ii = 1,neq
          tmp(ii,jj) = dfduRp(ii,jj)*DtmpRp(ii)
        end do
      end do
      soln%LHS(:,:,i,3) = tmp                        ! dRdu_(i+1)


    !==========================================================================
    !------------------------------ Rows 2:N-1 --------------------------------
    !==========================================================================
    do i = i_low+1,i_high-1

      FLm = FLi
      FRi = FRp
      dfduLm = dfduLi                                     ! dfdu_+_(i-1)
      dfduRi = dfduRp                                     ! dfdu_-_(i)
      call calc_vl_dfdu(soln%U(:,i  ),FLi,dfduLi,.true. ) ! dfdu_+_(i)
      call calc_vl_dfdu(soln%U(:,i+1),FRp,dfduRp,.false.) ! dfdu_-_(i+1)

      !soln%F(:,i-1) = FLm + FRi                      ! f_(i-1/2)
      soln%F(:,i)   = FLi + FRp                      ! f_(i+1/2)

      DtmpLm = grid%Ai(i-1)        ! A_(i-1/2)*duLdu_(i-1)
      DtmpRi = grid%Ai(i-1)        ! A_(i-1/2)*duRdu_(i-1)
      DtmpLi = grid%Ai(i)          ! A_(i+1/2)*duLdu_(i)
      DtmpRp = grid%Ai(i)          ! A_(i+1/2)*duRdu_(i)

      do jj = 1,neq
        do ii = 1,neq
          tmp(ii,jj) = dfduLm(ii,jj)*DtmpLm(ii)
        end do
      end do
      soln%LHS(:,:,i,1) = tmp                        ! dRdu_(i-1)

      do jj = 1,neq
        do ii = 1,neq
          tmp(ii,jj) = dfduLi(ii,jj)*DtmpLi(ii) + dfduRi(ii,jj)*DtmpRi(ii)
        end do
      end do
      soln%LHS(:,:,i,2) = tmp                        ! dRdu_(i)

      do jj = 1,neq
        do ii = 1,neq
          tmp(ii,jj) = dfduRp(ii,jj)*DtmpRp(ii)
        end do
      end do
      soln%LHS(:,:,i,3) = tmp                        ! dRdu_(i+1)

    end do

    !==========================================================================
    !------------------------------ Row N -------------------------------------
    !==========================================================================
    i = i_high

      FLm = FLi
      FRi = FRp
      dfduLm = dfduLi                                     ! dfdu_+_(i-1)
      dfduRi = dfduRp                                     ! dfdu_-_(i)
      call calc_vl_dfdu(soln%U(:,i  ),FLi,dfduLi,.true. ) ! dfdu_+_(i)
      call calc_vl_dfdu(soln%U(:,i+1),FRp,dfduRp,.false.) ! dfdu_-_(i+1)

      !soln%F(:,i-1) = FLm + FRi                      ! f_(i-1/2)
      soln%F(:,i)   = FLi + FRp                      ! f_(i+1/2)

      DtmpLm = grid%Ai(i-1)        ! A_(i-1/2)*duLdu_(i-1)
      DtmpRi = grid%Ai(i-1)        ! A_(i-1/2)*duRdu_(i-1)
      DtmpLi = grid%Ai(i)          ! A_(i+1/2)*duLdu_(i)
      DtmpRp = grid%Ai(i)          ! A_(i+1/2)*duRdu_(i)

      do jj = 1,neq
        do ii = 1,neq
          tmp(ii,jj) = dfduLm(ii,jj)*DtmpLm(ii)
        end do
      end do
      soln%LHS(:,:,i,1) = tmp                        ! dRdu_(i-1)

      do jj = 1,neq
        do ii = 1,neq
          tmp(ii,jj) = dfduLi(ii,jj)*DtmpLi(ii) + dfduRi(ii,jj)*DtmpRi(ii)
        end do
      end do
      soln%LHS(:,:,i,2) = tmp                        ! dRdu_(i)

      !do j = 1,3
      !  do i = i_low,i_high
      !    do jj = 1,neq
      !      do ii = 1,neq
      !        write(*,*) soln%LHS(ii,jj,i,j)
      !      end do
      !    end do
      !  end do
      !end do
      !do i = i_low,i_high
      !  do j = 1,2
      !    call disp('    ',soln%LHS(:,:,i,j), advance='NO')
      !  end do
      ! call disp('    ',soln%LHS(:,:,i,3))
      !end do
      !stop

  end subroutine build_LHS_matrix_1st_order

end module build_LHS_1st_order
