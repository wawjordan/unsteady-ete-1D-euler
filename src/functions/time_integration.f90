module time_integration

  use set_precision, only : prec
  use set_constants, only : zero, one, half, third, fourth
  use set_inputs,    only : neq, i_low, i_high, ig_low, ig_high
  use grid_type, only : grid_t
  use soln_type, only : soln_t
  use variable_conversion
  use residuals, only : calc_residual
  use build_LHS_1st_order, only : build_LHS_matrix_1st_order

  implicit none

  private

  public :: calc_time_step, explicit_RK
  public :: explicit_Euler, implicit_Euler
  contains

  !============================= calc_time_step ==============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine calc_time_step( dx, V, asnd, lambda, dt )

    use set_inputs,          only : CFL

    real(prec), dimension(neq,ig_low:ig_high), intent(in)  :: V
    real(prec), dimension(ig_low:ig_high),     intent(in)  :: dx
    real(prec), dimension(ig_low:ig_high),     intent(out) :: lambda
    real(prec), dimension(i_low:i_high),       intent(out) :: dt
    real(prec), dimension(ig_low:ig_high),     intent(in)  :: asnd

    lambda(:) = abs(V(2,:)) + asnd
    dt(:) = CFL*dx(i_low:i_high)/lambda(i_low:i_high)
    dt(:) = minval(dt)

  end subroutine calc_time_step

  !============================== explicit_RK ================================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine explicit_RK( grid, soln, VL, VR )
    use flux_calc, only : calc_flux_1D
    use variable_conversion, only : update_states, prim2cons
    use basic_boundaries, only : explicit_characteristic_bndry
    type(grid_t), intent(inout) :: grid
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq), intent(in) :: VL, VR
    real(prec), dimension(neq,i_low-1:i_high) :: leftV, rightV, leftU, rightU
    real(prec), dimension(neq,ig_low:ig_high) :: old_U
    real(prec), dimension(4) :: k
    integer :: i, j, n

    old_U = soln%U

    k = (/ fourth, third, half, one /)
    !k = one
    do j = 1,4
      call explicit_characteristic_bndry(soln, VL, VR)
      call prim2cons(soln%U,soln%V)
      call calc_flux_1D(grid,soln)
      call calc_residual(grid,soln)

      do i = 1,neq
        soln%U(i,i_low:i_high) = &
          old_U(i,i_low:i_high) - k(j)* &
          soln%R(i,i_low:i_high)*soln%dt(i_low:i_high)/ &
          ( grid%Ac(i_low:i_high)*grid%dx(i_low:i_high) )
      end do
      call update_states(soln)
    end do
  end subroutine explicit_RK



  subroutine implicit_euler(grid,soln)
    use dispmodule
    use tridiag_operations, only : trivec_alt, tprec_alt
    use flux_calc, only : calc_flux_1D
    type(grid_t), intent(inout) :: grid
    type(soln_t), intent(inout) :: soln
    integer :: N
    !real(prec), allocatable :: du(:,:)
    real(prec), allocatable :: AL(:,:,:,:), du(:,:)
    real(prec) :: d
    integer :: i, j, ii, jj

    N = i_high - i_low + 1

    !allocate( du(neq,N) )
    allocate( AL(neq,neq,N,3), du(neq,N) )

    call prim2cons(soln%U,soln%V)
    !call calc_flux_1D(grid,soln)
    call calc_residual(grid,soln)
    call build_LHS_matrix_1st_order(soln,grid)
    AL = soln%LHS
    do i = i_low, i_high
      do j = 1,neq
      AL(j,j,i,2) = AL(j,j,i,2) + &
          ( grid%Ac(i)*grid%dx(i)/soln%dt(i) )
      end do
      !call disp('AL = ',AL(:,:,i,2))
    end do

    !do j = 1,3
    !  do i = i_low, i_high
    !    do jj = 1,neq
    !      do ii = 1,neq
    !        write(*,*) AL(ii,jj,i,j)
    !      end do
    !    end do
    !  end do
    !end do
    !stop

    du = -soln%R(:,:)

    !do i = i_low,i_high
    !  do j = 1,neq
    !    write(*,*) du(j,i)
    !  end do
    !end do
    !stop

    call trivec_alt(N,neq,AL)
    call tprec_alt(N,neq,AL,du)

    do i = i_low, i_high
      soln%U(:,i) = soln%U(:,i) + du(:,i)
      !call disp(soln%U(:,i))
    end do

    deallocate( AL, du )

  end subroutine implicit_euler

!  subroutine implicit_euler_old(grid,soln)
!    use block_matrix_operations, only : blk_bandec, blk_banbks
!    use residuals, only : build_LHS_matrix
!    use flux_calc, only : calc_flux_1D
!    type(grid_t), intent(inout) :: grid
!    type(soln_t), intent(inout) :: soln
!    integer :: M1, M2, MP, MPL, N, NP, q
!    integer :: indx(1:i_high-i_low+1)
!    real(prec), allocatable :: AL(:,:,:,:), du(:,:)
!    real(prec) :: d
!    integer :: i
!
!    M1 = 2
!    M2 = 2
!    MP = M1 + M2 + 1
!    MPL = M1
!    N = i_high - i_low + 1
!    NP = N
!    q = neq
!
!    allocate( AL(q,q,N,M1), du(q,N) )
!
!    !call prim2cons(soln%U,soln%V)
!    call calc_flux_1D(grid,soln)
!    call calc_residual(grid,soln)
!    call build_LHS_matrix(soln,grid)
!
!    AL = zero
!    du = -soln%R(:,:)
!
!    call blk_bandec(soln%LHS,N,q,M1,M2,NP,MP,AL,MPL,indx,d)
!
!    call blk_banbks(soln%LHS,N,q,M1,M2,NP,MP,AL,MPL,indx,du)
!
!    do i = i_low, i_high
!      soln%U(:,i) = soln%U(:,i) + du(:,i)
!    end do
!
!    deallocate( AL, du )
!
!  end subroutine implicit_euler_old

  !============================= explicit_euler ==============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine explicit_euler( grid, src, dt, F, U, R )

    type(grid_t), intent(in) :: grid

    real(prec), dimension(neq,ig_low:ig_high), intent(inout) :: U
    real(prec), dimension(neq,i_low:i_high),   intent(out)   :: R
    real(prec), dimension(neq,i_low-1:i_high), intent(in)    :: F
    real(prec), dimension(i_low:i_high)  ,   intent(in)    :: src, dt

    !integer :: i

    R(1,i_low:i_high) = F(1,i_low:i_high)*grid%Ai(i_low:i_high) &
           - F(1,i_low-1:i_high-1)*grid%Ai(i_low-1:i_high-1)
    R(2,i_low:i_high) = F(2,i_low:i_high)*grid%Ai(i_low:i_high) &
           - F(2,i_low-1:i_high-1)*grid%Ai(i_low-1:i_high-1) &
           - 0.0_prec*src(i_low:i_high)*grid%dx(i_low:i_high)
    R(3,i_low:i_high) = F(3,i_low:i_high)*grid%Ai(i_low:i_high) &
           - F(3,i_low-1:i_high-1)*grid%Ai(i_low-1:i_high-1)


    U(1,i_low:i_high) = U(1,i_low:i_high) &
                      - dt(i_low:i_high)/ &
                        (grid%Ac(i_low:i_high)*grid%dx(i_low:i_high))*R(1,:)
    U(2,i_low:i_high) = U(2,i_low:i_high) &
                      - dt(i_low:i_high)/ &
                        (grid%Ac(i_low:i_high)*grid%dx(i_low:i_high))*R(2,:)
    U(3,i_low:i_high) = U(3,i_low:i_high) &
                      - dt(i_low:i_high)/ &
                        (grid%Ac(i_low:i_high)*grid%dx(i_low:i_high))*R(3,:)

  end subroutine explicit_euler

end module time_integration
