module quadrature

  use set_precision, only : prec
  use set_constants, only : pi, zero, one, two, three, four, half, fourth
  use set_inputs,    only : imax, i_low, i_high
 ! use grid_type,     only : grid_t

  implicit none

  private

  public :: gauleg, gauss_integ

  contains

    subroutine gauleg(x,w,tol,n,max_iter)
      real(prec), intent(in) :: tol
      integer, intent(in) :: n, max_iter
      real(prec), intent(out) :: x(n), w(n)
      integer :: i, j, k, m
      real(prec) :: x1, x2, p1, p2, p3, pp, xl, xm, z, z1
      x1 = -one
      x2 = one
      m = (n+1)/2
      xm = half*(x2+x1)
      xl = half*(x2-x1)
      do i = 1,m
        z = cos(pi*(real(i,prec)-fourth)/(real(n,prec)+half))
        do k = 1,max_iter
          p1 = one
          p2 = zero
          do j = 1,n
            p3 = p2
            p2 = p1
            p1 = ( ( two*real(j,prec) - one )*z*p2 &
                 - ( real(j,prec) - one )*p3 )/real(j,prec)
          end do
          pp = real(n,prec)*(z*p1-p2)/(z*z-one)
          z1 = z
          z = z1 - p1/pp
          if ((abs(p1) < tol).and.(abs(z-z1)/abs(z) < tol)) exit
        end do
        x(i) = xm - xl*z
        x(n+1-i) = xm + xl*z
        w(i) = two*xl/((one-z*z)*pp*pp)
        w(n+1-i) = w(i)
      end do
    end subroutine gauleg

    subroutine gauss_integ(I,f,a,b,xi,w)
      real(prec), intent(out) :: I
      real(prec), external :: f
      real(prec), intent(in) :: a,b
      real(prec), intent(in) :: xi(:),w(:)
      real(prec), allocatable :: x(:)
      real(prec) :: c1, c2, dxdxi
      integer :: n, j
      n = size(xi,1)
      allocate(x(n))
      c1 = (b-a)/two
      c2 = (a+b)/two
      dxdxi = (b-a)/two

      x = c1*xi + c2
      I = zero
      do j = 1,n
        I = I + w(j)*f(x(j))
      end do
      I = dxdxi*I

      deallocate(x)
    end subroutine gauss_integ

end module quadrature
