module soln_type

  use set_precision, only : prec
  use set_constants, only : zero, one
  use set_inputs,    only : i_low, i_high, ig_low, ig_high
  use set_inputs,    only : neq, max_iter, limiter_scheme
  use limiter_type,  only : limiter_t, select_limiter

  implicit none

  type soln_t

    real(prec), allocatable, dimension(:,:,:,:) :: LHS
    real(prec), allocatable, dimension(:,:) :: V
    real(prec), allocatable, dimension(:,:) :: U
    real(prec), allocatable, dimension(:,:) :: R
    real(prec), allocatable, dimension(:,:) :: S
    real(prec), allocatable, dimension(:,:) :: F
    real(prec), allocatable, dimension(:,:) :: psi_p
    real(prec), allocatable, dimension(:,:) :: psi_m
    real(prec), allocatable, dimension(:,:) :: dpsi_p
    real(prec), allocatable, dimension(:,:) :: dpsi_m
    real(prec), allocatable, dimension(:,:) :: DE
    real(prec), allocatable, dimension(:)   :: asnd
    real(prec), allocatable, dimension(:)   :: mach
    real(prec), allocatable, dimension(:)   :: temp
    real(prec), allocatable, dimension(:)   :: lambda
    real(prec), allocatable, dimension(:)   :: dt
    real(prec), allocatable, dimension(:)   :: rinit
    real(prec), allocatable, dimension(:)   :: rold
    real(prec), allocatable, dimension(:)   :: rnorm
    real(prec), allocatable, dimension(:)   :: DEnorm
    real(prec) :: time
    type(limiter_t) :: lim

  end type soln_t

  contains


  !============================= allocate_soln ===============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine allocate_soln( soln )

    type(soln_t), intent(inout) :: soln

    allocate( soln%LHS( neq, neq, i_low:i_high, 3 ) )

    allocate( soln%V(  neq,  ig_low:ig_high ), &
              soln%U(  neq,  ig_low:ig_high ), &
              soln%R(  neq,   i_low:i_high ),  &
              soln%S(  neq,   i_low:i_high ),  &
              soln%F(  neq, i_low-1:i_high ),  &
              soln%DE( neq,   i_low:i_high ),  &
              soln%asnd( ig_low:ig_high ),     &
              soln%mach( ig_low:ig_high ),     &
              soln%temp( ig_low:ig_high ),     &
              soln%lambda( ig_low:ig_high ),   &
              soln%dt( i_low:i_high ),         &
              soln%rinit( neq ),               &
              soln%rold( neq ),                &
              soln%rnorm( neq ),               &
              soln%DEnorm( neq ) )
    allocate( soln%psi_p(  neq, ig_low-1:ig_high ), &
              soln%psi_m(  neq, ig_low-1:ig_high ), &
              soln%dpsi_p( neq, ig_low-1:ig_high ), &
              soln%dpsi_m( neq, ig_low-1:ig_high ) )

    soln%LHS    = zero
    soln%V      = zero
    soln%U      = zero
    soln%R      = zero
    soln%S      = zero
    soln%F      = zero
    soln%DE     = zero
    soln%asnd   = zero
    soln%mach   = zero
    soln%temp   = zero
    soln%lambda = zero
    soln%dt     = zero
    soln%rinit  = zero
    soln%rold   = zero
    soln%rnorm  = zero
    soln%DEnorm = zero
    soln%psi_p  = zero
    soln%psi_m  = zero
    soln%dpsi_p = zero
    soln%dpsi_m = zero
    soln%time   = zero

    call soln%lim%select_limiter(limiter_scheme)

  end subroutine allocate_soln


  !=========================== deallocate_soln ===============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine deallocate_soln( soln )

    implicit none

    type(soln_t), intent(inout) :: soln

    deallocate(soln%LHS,    &
               soln%V,      &
               soln%U,      &
               soln%R,      &
               soln%S,      &
               soln%F,      &
               soln%DE,     &
               soln%asnd,   &
               soln%mach,   &
               soln%temp,   &
               soln%lambda, &
               soln%dt,     &
               soln%rinit,  &
               soln%rold,   &
               soln%rnorm,  &
               soln%DEnorm, &
               soln%psi_p,  &
               soln%psi_m,  &
               soln%dpsi_p, &
               soln%dpsi_m  )

  end subroutine deallocate_soln

end module soln_type
