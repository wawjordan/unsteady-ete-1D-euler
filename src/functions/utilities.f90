module utilities

  use set_precision, only : prec
  use set_constants, only : zero, one, two, half

  implicit none

  public :: newton_safe, hunt, block_band2full, &
          block_band2sparse1, block_band2sparse2

contains

  subroutine block_band2sparse2(A1,M1,M2,q,N)
    real(prec), intent(in)  :: A1(q,q,N,M1+M2+1)
    integer, intent(in) :: M1, M2, q, N
    real(prec) :: val
    integer :: MM, i, j, k, ii, jj, row, col
    MM = M1 + M2 + 1    ! block bandwidth
    do k = 1,N
    do j = 1,MM
    do ii = 1,q
    do jj = 1,q
    row = (k-1)*q + ii
    col = (j-1)*q + jj
    val = A1(ii,jj,k,j)
    write(*,*) row, col, val
    end do
    end do
    end do
    end do
  end subroutine block_band2sparse2

  subroutine block_band2sparse1(A1,M1,M2,q,N)
    real(prec), intent(in)  :: A1(q,q,N,M1+M2+1)
    integer, intent(in) :: M1, M2, q, N
    real(prec) :: val, A2(q,q,N,M1+M2+1)
    integer :: L, MM, i, j, k, ii, jj, row, col

    MM = M1 + M2 + 1    ! block bandwidth
    A2 = A1
    L = M1
    do i = 1,M1
      do j = M1+2-i,MM
        A2(:,:,i,j-L)  = A1(:,:,i,j)
      end do
      L = L - 1
      do j = MM-L,MM
        A2(:,:,i,j) = zero
      end do
    end do
    L = 0
    do k = 1,N
      if (k>M1+1) then
      if (L<N) then
        L = L + 1
      end if
      end if
      do j = 1,min(MM,N-(k-M2-1))
        do ii = 1,q
          do jj = 1,q
            ! row, column, value
            row = (k-1)*q + ii
            col = (L+j-1)*q + jj
            val = A2(ii,jj,k,j)
            write(*,*) row, col, val
          end do
        end do
      end do
    end do

  end subroutine block_band2sparse1



  subroutine block_band2full(A1,A2,M1,M2,q,N)
    real(prec), intent(in)  :: A1(q,q,N,M1+M2+1)
    real(prec), intent(out) :: A2(q*N,q*N)
    integer, intent(in) :: M1, M2, q, N
    integer :: MM, i, j, k, ii, jj, h

    MM = M1 + M2 + 1    ! block bandwidth

    A2 = zero

    do ii = 1,N
      h = (ii-1)*q
      do j = 1,q
        do i = 1,q
          A2(h+i,h+j)  = A1(i,j,ii,M1+1)
        end do
      end do
    end do

    do jj = 1,M1
      k = M1 + 1 - jj
      do ii = M1+2-jj,N
        h = (ii-1)*q
        do j = 1,q
          do i = 1,q
            A2(h+i,(h-k*q)+j) = A1(i,j,ii,jj)
          end do
        end do
      end do
    end do

    do jj= M1+2,MM
      k = MM + 1 - jj
      do ii = 1,N-k
        h = (ii-1)*q
        do j = 1,q
          do i = 1,q
            A2(h+i,(h+k*q)+j) = A1(i,j,ii,jj)
          end do
        end do
      end do
    end do

  end subroutine block_band2full

  subroutine newton_safe( fun, dfun, guess, bnd_a1, bnd_b1, &
                                      tol, max_iter, x, err )

    real(prec), external :: fun, dfun
    real(prec), intent(in)  :: guess, bnd_a1, bnd_b1, tol
    integer,    intent(in)  :: max_iter
    real(prec), intent(out) :: x, err
    real(prec), dimension(2) :: xk
    real(prec) :: a, b, fk, dfk, tmp
    integer :: k

    xk = zero

    a = bnd_a1
    b = bnd_b1
    xk(1) = guess
    fk = fun(xk(1))
    dfk = dfun(xk(1))
    xk(2) = xk(1) - fk/dfk

    if ( (xk(2)<a).or.(xk(2)>b) ) then
      !xk(2) = a + (b-a)/2
      xk(2) = guess
      tmp = fun(a)
      if ( int(sign(one,tmp)) == int(sign(one,fk)) ) then
        a = xk(2)
      else
        b = xk(2)
      end if
    end if

    err = one

    do k = 2, max_iter
      if ( (fk < tol).and.(err < tol) ) exit
      fk = fun(xk(2))
      dfk = dfun(xk(2))
      tmp = xk(2) - fk/dfk
      if ( (tmp<a).or.(tmp>b) ) then
        xk(1) = xk(2)
        xk(2) = a + (b-a)/2
        tmp = fun(a)
        if ( int(sign(one,tmp)) == int(sign(one,fk)) ) then
          a = xk(2)
        else
          b = xk(2)
        end if
      else
        xk(1) = xk(2)
        xk(2) = tmp
      end if
      err = abs(xk(2)-xk(1))/abs(xk(2))
    end do

    x = xk(2)

  end subroutine newton_safe

  subroutine hunt(xx,x,jlo)
    real(prec), intent(in) :: xx(:)
    real(prec), intent(in) :: x
    integer, intent(inout) :: jlo
    integer :: n, inc, jhi, jm
    logical :: ascnd

    n = size(xx,1)

    ascnd=(xx(n) >= xx(1))
    if (jlo <= 0 .or. jlo > n) then      ! check if guess is useful
      jlo = 0
      jhi = n + 1
    else
      inc = 1                            ! set hunting increment
      if (x >= xx(jlo) .eqv. ascnd) then ! hunt up
        do
          jhi = jlo + inc
          if ( jhi > n ) then            ! done hunting (off end of table)
            jhi = n + 1
            exit
          else
            if (x < xx(jhi) .eqv. ascnd) exit
            jlo = jhi                    ! not done hunting
            inc = inc + inc              ! double increment
          end if
        end do                           ! try again
      else                               ! hunt down
        jhi = jlo
        do
          jlo = jhi - inc
          if (jlo < 1) then              ! done hunting (off end of table)
            jlo = 0
            exit
          else
            if (x >= xx(jlo) .eqv. ascnd ) exit
            jhi = jlo                    ! not done hunting
            inc = inc + inc              ! double increment
          end if
        end do                           ! try again
      end if
    end if                               ! done hunting (value bracketed)
    do                                   ! bisection
      if (jhi - jlo <= 1) then
        if (x == xx(n)) jlo = n - 1
        if (x == xx(n)) jlo = 1
        exit
      else
        jm = (jhi + jlo)/2
        if (x >= xx(jm) .eqv. ascnd) then
          jlo = jm
        else
          jhi = jm
        end if
      end if
    end do

  end subroutine hunt

end module utilities
