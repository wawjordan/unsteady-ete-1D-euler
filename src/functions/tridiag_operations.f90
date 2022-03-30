module tridiag_operations

  use set_precision, only : prec
  use set_constants, only : zero, one

  implicit none

  public :: trivec, tprec

contains

!  subroutine trivec(jj,neq,da,dd,dc)
!    implicit doubleprecision(a-h,o-z)
!    dimension da(201,3,3),dd(201,3,3),dc(201,3,3)
!    do 1 j=1,jj
!    do 2 k=1,neq-1
!    do 3 l=k+1,neq
!      dd(j,l,k) = dd(j,l,k)/dd(j,k,k)
!3   continue
!    do 4 l=k+1,neq
!    do 4 m=k+1,neq
!      dd(j,l,m) = dd(j,l,m) - dd(j,l,k)*dd(j,k,m)
!4   continue
!2   continue
!    if(j .ne. jj) then
!    do 5 l=1,neq-1
!    do 5 m=l+1,neq
!    do 5 mm=1,neq
!      dc(j,m,mm) = dc(j,m,mm) - dd(j,m,l)*dc(j,l,mm)
!5   continue
!    do 6 l=neq,1,-1
!    do 7 mm=1,neq
!7     dc(j,l,mm) = dc(j,l,mm)/dd(j,l,l)
!    do 8 m=1,l-1
!    do 8 mm=1,neq
!      dc(j,m,mm) = dc(j,m,mm) - dd(j,m,l)*dc(j,l,mm)
!8   continue
!6   continue
!    do 9 k=1,neq
!    do 9 l=1,neq
!    do 9 mm=1,neq
!      dd(j+1,l,k) = dd(j+1,l,k) - da(j+1,l,mm)*dc(j,mm,k)
!9   continue
!    endif
!1   continue
!    return
!  end subroutine trivec

  ! **********************************************************************
  !
  !     Subroutine TRIVEC does the LU decomposition of the block
  !     tridiagonal matrices DA,DD,DC
  !
  ! **********************************************************************
  subroutine trivec(jj,neq,da,dd,dc)
    implicit none
    integer, intent(in) :: jj, neq
    real(prec), intent(inout) :: da(jj,neq,neq), dd(jj,neq,neq), dc(jj,neq,neq)
    integer :: j, k, l, m, mm

    do j = 1,jj           ! 1 *********************************************** 1
      do k = 1,neq-1      ! 2 ******************************************** 2  |
        do l = k+1,neq    ! 3 ***************************************** 3  |  |
          dd(j,l,k) = dd(j,l,k)/dd(j,k,k) !                             |  |  |
        end do            ! 3 continue ******************************** 3  |  |
        do l = k+1,neq    ! 4 ***************************************** 4  |  |
          do m = k+1,neq  ! 4 *************************************** 4 |  |  |
            dd(j,l,m) = dd(j,l,m) - dd(j,l,k)*dd(j,k,m) !             | |  |  |
          end do          ! 4 continue ****************************** 4 |  |  |
        end do            ! 4 continue ******************************** 4  |  |
      end do              ! 2 continue *********************************** 2  |
      if (j /= jj) then   ! ================================================= |
        do l = 1,neq-1    ! 5 ******************************************** 5  |
          do m = l+1,neq  ! 5 ****************************************** 5 |  |
            do mm = 1,neq ! 5 **************************************** 5 | |  |
              dc(j,m,mm) = dc(j,m,mm) - dd(j,m,l)*dc(j,l,mm) !         | | |  |
            end do        ! 5 continue ******************************* 5 | |  |
          end do          ! 5 continue ********************************* 5 |  |
        end do            ! 5 continue *********************************** 5  |
        do l = neq,1,-1   ! 6 ******************************************** 6  |
          do mm = 1,neq   ! 7 ***************************************** 7  |  |
            dc(j,l,mm) = dc(j,l,mm)/dd(j,l,l) !                         |  |  |
          end do          ! 7 continue (?) **************************** 7  |  |
          do m = 1,l-1    ! 8 ***************************************** 8  |  |
            do mm = 1,neq ! 8 *************************************** 8 |  |  |
              dc(j,m,mm) = dc(j,m,mm) - dd(j,m,l)*dc(j,l,mm) !        | |  |  |
            end do        ! 8 continue ****************************** 8 |  |  |
          end do          ! 8 continue ******************************** 8  |  |
        end do            ! 6 continue *********************************** 6  |
        do k = 1,neq      ! 9 ******************************************** 9  |
          do l = 1,neq    ! 9 ****************************************** 9 |  |
            do mm = 1,neq ! 9 **************************************** 9 | |  |
              dd(j+1,l,k) = dd(j+1,l,k) - da(j+1,l,mm)*dc(j,mm,k) !    | | |  |
            end do        ! 9 continue ******************************* 9 | |  |
          end do          ! 9 continue ********************************* 9 |  |
        end do            ! 9 continue *********************************** 9  |
      endif               ! ================================================= |
    end do                ! 1 continue ************************************** 1
  end subroutine trivec


!  subroutine tprec(jj,neq,da,dd,dc,c)
!
!      Use Select_Precision
!
!      implicit none
!
!      Integer :: jj
!      Integer :: neq
!      Integer :: i = -99
!      Integer :: k = -99
!      Integer :: l = -99
!      Integer :: m = -99
!!     Integer :: mm = -99
!      Real(kind=Prec) :: sum  = -99.9_Prec
!      Real(kind=Prec),Dimension(jj,neq,neq) :: da
!      Real(kind=Prec),Dimension(jj,neq,neq) :: dd
!      Real(kind=Prec),Dimension(jj,neq,neq) :: dc
!      Real(kind=Prec),Dimension(jj,neq)     :: c
!      do 1 k=1,jj
!        if(k .ne. 1) then
!          do 2 i=1,neq
!           sum = da(k,i,1)*c(k-1,1)
!          do 22 m=2,neq
!    22     sum = sum + da(k,i,m)*c(k-1,m)
!           c(k,i) = c(k,i) - sum
!    2     continue
!        endif
! ! ----- forward solve
!        do 4 l=1,neq-1
!        do 4 m=l+1,neq
!          c(k,m) = c(k,m) - dd(k,m,l)*c(k,l)
!  4     continue
! ! ----- back solve
!        do 5 l=neq,1,-1
!          c(k,l) = c(k,l)/dd(k,l,l)
!          do 6 m=l-1,1,-1
!    6       c(k,m) = c(k,m) - dd(k,m,l)*c(k,l)
!  5     continue
! 1    continue
!      do 44 k=jj-1,1,-1
!      do 44 i=1,neq
!        sum = dc(k,i,1)*c(k+1,1)
!        do 54 m=2,neq
!   54    sum = sum + dc(k,i,m)*c(k+1,m)
!        c(k,i) = c(k,i) - sum
! 44   continue
!      return
!      end


  ! *******************************************************************
  !
  !     Subroutine TPREC is a non-vectorizable backsolve routine for
  !     solving the block-tridiagonal system LUx=b given a factorization
  !
  ! *******************************************************************
  subroutine tprec(jj,neq,da,dd,dc,c)
    implicit none
    integer :: jj
    integer :: neq
    integer :: i = -99
    integer :: k = -99
    integer :: l = -99
    integer :: m = -99
!   integer :: mm = -99
    real(prec) :: sum  = -99.9_prec
    real(prec) :: da(jj,neq,neq)
    real(prec) :: dd(jj,neq,neq)
    real(prec) :: dc(jj,neq,neq)
    real(prec) :: c(jj,neq)

    do k = 1,jj           ! 1 *********************************************** 1
      if (k /= 1) then    ! ================================================= |
        do i = 1,neq      ! 2 ******************************************** 2  |
          sum = da(k,i,1)*c(k-1,1) !                                       |  |
          do m = 2,neq    ! 22 **************************************** 22 |  |
            sum = sum + da(k,i,m)*c(k-1,m) !                            |  |  |
          end do          ! 22 continue ******************************* 22 |  |
          c(k,i) = c(k,i) - sum !                                          |  |
        end do            ! 2 continue *********************************** 2  |
      endif               ! ================================================= |
 ! ----- forward solve                                                        |
      do l = 1,neq-1      ! 4 ******************************************** 4  |
        do m = l+1,neq    ! 4 ****************************************** 4 |  |
          c(k,m) = c(k,m) - dd(k,m,l)*c(k,l) !                           | |  |
        end do            ! 4 continue ********************************* 4 |  |
      end do              ! 4 continue *********************************** 4  |
 ! ----- back solve                                                           |
      do l = neq,1,-1     ! 5 ******************************************** 5  |
        c(k,l) = c(k,l)/dd(k,l,l) !                                        |  |
        do m = l-1,1,-1   ! 6 ***************************************** 6  |  |
          c(k,m) = c(k,m) - dd(k,m,l)*c(k,l) !                          |  |  |
        end do            ! 6 continue ******************************** 6  |  |
      end do              ! 5 continue *********************************** 5  |
    end do                ! 1 continue ************************************** 1
    do k = jj-1,1,-1      ! 44 ********************************************* 44
      do i = 1,neq        ! 44 ******************************************* 44 |
        sum = dc(k,i,1)*c(k+1,1) !                                          | |
        do m = 2,neq      ! 54 **************************************** 54  | |
          sum = sum + dc(k,i,m)*c(k+1,m) !                              |   | |
        end do            ! 54 continue ******************************* 54  | |
        c(k,i) = c(k,i) - sum !                                             | |
      end do              ! 44 continue ********************************** 44 |
    end do                ! 44 continue ************************************ 44
  end subroutine tprec

end module tridiag_operations
