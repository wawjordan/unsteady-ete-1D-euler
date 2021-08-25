module geometry

  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, half, pi
  use fluid_constants, only : R_gas, gamma
  use set_inputs,      only : area, darea, eps
  use soln_type
  use grid_type

  implicit none
  private
  public :: setup_geometry, teardown_geometry
  contains

  !============================= setup geometry  =============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine setup_geometry( grid, soln )

    use set_inputs, only : i_low, i_high

    type( soln_t ), intent(inout) :: soln
    type( grid_t ), intent(inout) :: grid

    integer :: i

    call allocate_grid( grid )
    call allocate_soln( soln )

    do i = i_low-1,i_high
      grid%Ai(i) = area( grid%xi(i) )
    end do

    do i = i_low,i_high
      grid%Ac(i) = area( grid%xc(i) )
    end do

    !do i = i_low,i_high
    !  grid%dAc(i) = darea( grid%xc(i) )
    !end do
    !grid%dAc(ig_low:i_low-1) = zero
    do i = i_low,i_high
      grid%dAc(i) = ( grid%Ai(i) - grid%Ai(i-1) )/grid%dx(i)
    end do
    !grid%dAc(i_high+1:ig_high) = zero

  end subroutine setup_geometry
!  subroutine setup_geometry( grid, soln, xi )
!
!    use set_inputs, only : i_low, i_high
!
!    type( soln_t ), intent(inout) :: soln
!    type( grid_t ), intent(inout) :: grid
!    real(prec), dimension(:), intent(in) :: xi
!
!    integer :: i
!
!    call allocate_grid( grid, xi )
!    call allocate_soln( soln )
!
!    do i = i_low-1,i_high
!      grid%Ai(i) = area( grid%xi(i) )
!    end do
!
!    do i = i_low,i_high
!      grid%Ac(i) = area( grid%xc(i) )
!    end do
!
!    !do i = i_low,i_high
!    !  grid%dAc(i) = darea( grid%xc(i) )
!    !end do
!    !grid%dAc(ig_low:i_low-1) = zero
!    do i = i_low,i_high
!      grid%dAc(i) = ( grid%Ai(i) - grid%Ai(i-1) )/grid%dx(i)
!    end do
!    !grid%dAc(i_high+1:ig_high) = zero
!
!  end subroutine setup_geometry

  !========================== teardown geometry  =============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine teardown_geometry( grid, soln )

    type( soln_t ) :: soln
    type( grid_t ) :: grid

    call deallocate_grid( grid )
    call deallocate_soln( soln )

  end subroutine teardown_geometry

end module geometry
