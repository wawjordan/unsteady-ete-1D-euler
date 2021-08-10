module extrapolation

  use set_precision, only : prec
  use set_constants, only : one, two, three, four, five, six, seven, eight, &
                            nine, ten, half, third, fourth, fifth, sixth,   &
                            eighth, ninth, tenth, small
  use set_inputs,    only : neq, ig_low, ig_high
  use soln_type,     only : soln_t

  implicit none

  private

  public :: polint

  contains

  subroutine polint(xa,ya,n,x,y,dy)

    integer,    intent(in) :: n
    real(prec), intent(in) :: xa(n), ya(n)
    real(prec), intent(in) :: x
    real(prec), intent(inout) :: y, dy
    integer :: i, m, ns
    real(prec) :: den, dif, dift, ho, hp, w
    real(prec), dimension(ubound(xa,1)-lbound(xa,1)) :: c, d

    ns  = 1
    dif = abs(x - xa(1))
    do i = 1,n
      dift = abs(x-xa(i))
      if (dift < dif) then
        ns  = i
        dif = dift
      end if
      c(i) = ya(i)
      d(i) = ya(i)
    end do
    y = ya(ns)
    ns = ns - 1
    do m = 1,n-1
      do i = 1,n-m
        ho  = xa(i) - x
        hp  = xa(i+m) - x
        w   = c(i+1) - d(i)
        den = ho - hp
        if (den <= small) then
          write(*,*) 'failure in polint'
        end if
        den = w/den
        d(i) = hp*den
        c(i) = ho*den
      end do
      if (2*ns < n-m) then
        dy = c(ns+1)
      else
        dy = d(ns)
        ns = ns - 1
      end if
      y = y + dy
    end do

  end subroutine polint

end module extrapolation


