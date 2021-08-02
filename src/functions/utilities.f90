module utilities

  use set_precision, only : prec
  use set_constants, only : zero, one, two, half

  implicit none

  public :: newton_safe, hunt

contains

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
