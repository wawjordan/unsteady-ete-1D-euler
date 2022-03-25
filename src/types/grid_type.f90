module grid_type

  use set_precision, only : prec

  implicit none

  type grid_t

    integer :: imax

    real(prec) :: xmin
    real(prec) :: xmax

    real(prec), allocatable, dimension(:) :: xi
    real(prec), allocatable, dimension(:) :: xc
    real(prec), allocatable, dimension(:) :: dx
    real(prec), allocatable, dimension(:) :: Ai
    real(prec), allocatable, dimension(:) :: Ac
    real(prec), allocatable, dimension(:) :: dAc

  end type grid_t

  contains

  !============================= allocate_grid ===============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine allocate_grid( grid )

    use set_constants, only : zero, one, half
    use set_inputs   , only : imax, ig_low, ig_high, i_low, i_high, &
                              n_ghost, xmin, xmax

    type( grid_t ), intent( inout ) :: grid
    integer :: i

    allocate( grid%xi(ig_low-1:ig_high), &
              grid%xc(ig_low:ig_high)  , &
              grid%dx(ig_low:ig_high)  , &
              grid%Ai(ig_low-1:ig_high), &
              grid%Ac(ig_low:ig_high)  , &
              grid%dAc(i_low:i_high)   )

    grid%dx = (xmax - xmin)/real(imax,prec)


    grid%xi = [ (xmin + real(i,prec)/real(imax,prec)*(xmax-xmin), &
               i=ig_low-1,ig_high) ]
    grid%xc = [ (half*(grid%xi(i) + grid%xi(i+1)),i=ig_low-1,ig_high-1) ]
    !do i = ig_low,ig_high
    !  write(*,'(I3,3G12.4)') i, grid%xi(i), grid%xc(i), grid%dx(i)
    !end do

    grid%Ai  = one
    grid%Ac  = one
    grid%dAc = zero

  end subroutine allocate_grid

! ==> Grid using specified node locations
!  subroutine allocate_grid( grid, xi )
!
!    use set_constants, only : zero, one, half
!    use set_inputs   , only : imax, ig_low, ig_high, i_low, i_high, &
!                              n_ghost, xmin, xmax
!
!    type( grid_t ), intent( inout ) :: grid
!    real(prec), dimension(:), intent(in) :: xi
!    real(prec) :: dx_left, dx_right
!    integer :: i, nmax
!
!    nmax = size(xi)
!    if (i /= imax) then
!      imax    = nmax-1
!      i_low   = 1
!      i_high  = imax
!      ig_low  = 1 - n_ghost
!      ig_high = imax + n_ghost
!    end if
!    xmin = minval(xi)
!    xmax = maxval(xi)
!
!    allocate( grid%xi(ig_low-1:ig_high), &
!              grid%xc(ig_low:ig_high)  , &
!              grid%dx(ig_low:ig_high)  , &
!              grid%Ai(ig_low-1:ig_high), &
!              grid%Ac(ig_low:ig_high)  , &
!              grid%dAc(i_low:i_high)   )
!
!    grid%xi = zero
!
!    grid%xi(i_low-1:i_high) = xi
!
!    grid%dx(i_low:i_high) = xi(2:nmax) - xi(1:nmax-1)
!
!    grid%dx(i_low-1:ig_low:-1) = grid%dx(i_low)
!    grid%dx(i_high+1:ig_high) = grid%dx(i_high)
!
!
!    do i = i_low-1,ig_low,-1
!      grid%xi(i-1) = grid%xi(i) - grid%dx(i_low)
!    end do
!    do i = i_high+1,ig_high
!      grid%xi(i) = grid%xi(i-1) + grid%dx(i_high)
!    end do
!
!    !grid%dx = (xmax - xmin)/real(imax,prec)
!
!
!    !grid%xi = [ (xmin + real(i,prec)/real(imax,prec)*(xmax-xmin), &
!    !           i=ig_low-1,ig_high) ]
!    grid%xc = [ (half*(grid%xi(i) + grid%xi(i+1)),i=ig_low-1,ig_high-1) ]
!    !do i = ig_low,ig_high
!    !  write(*,'(I3,3G12.4)') i, grid%xi(i), grid%xc(i), grid%dx(i)
!    !end do
!
!    grid%Ai = one
!    grid%Ac = one
!    grid%dAc = zero
!
!  end subroutine allocate_grid

  !============================= deallocate_grid =============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine deallocate_grid( grid )

    type( grid_t ), intent( inout ) :: grid

    deallocate( grid%xi, grid%xc, grid%Ai, grid%Ac, grid%dAc )

  end subroutine deallocate_grid

end module grid_type

