module soln_type

  use set_precision, only : prec
  use set_constants, only : zero, one
  use set_inputs,    only : i_low, i_high, ig_low, ig_high
  use set_inputs,    only : neq, max_iter, shock
  
  implicit none

  type soln_t
    
    real(prec), allocatable, dimension(:,:) :: V
    real(prec), allocatable, dimension(:,:) :: U
    real(prec), allocatable, dimension(:,:) :: R
    real(prec), allocatable, dimension(:,:) :: F
    real(prec), allocatable, dimension(:,:) :: D
    real(prec), allocatable, dimension(:)   :: asnd
    real(prec), allocatable, dimension(:)   :: mach
    real(prec), allocatable, dimension(:)   :: temp
    real(prec), allocatable, dimension(:)   :: dt
    real(prec), allocatable, dimension(:)   :: src
    real(prec), allocatable, dimension(:)   :: lambda
    real(prec), allocatable, dimension(:,:) :: DE
    real(prec), allocatable, dimension(:)   :: DEnorm
    real(prec), allocatable, dimension(:)   :: rnorm
    real(prec), allocatable, dimension(:)   :: rold
    real(prec), allocatable, dimension(:)   :: rinit
    
  end type soln_t

  contains
  
  
  !============================= allocate_soln ===============================80
  !>
  !! Description: 
  !!
  !! Inputs:      soln : 
  !!
  !! Outputs:     soln : 
  !<
  !===========================================================================80
  subroutine allocate_soln( soln )

    type(soln_t), intent(inout) :: soln
    
    allocate( soln%V( ig_low:ig_high, neq ), &
              soln%U( ig_low:ig_high, neq ), &
              soln%R( i_low:i_high,   neq ), &
              soln%F( i_low-1:i_high, neq ), &
              soln%D( i_low-1:i_high, neq ), &
              soln%asnd( ig_low:ig_high ),      &
              soln%mach( ig_low:ig_high ),      &
              soln%temp( ig_low:ig_high ),      &
              soln%src( i_low:i_high ),      &
              soln%dt( i_low:i_high ),     &
              soln%lambda( ig_low:ig_high ), &
              soln%rnorm( 1:neq ),        &
              soln%rold( 1:neq ),        &
              soln%rinit( 1:neq ) )
    
    if (shock.eq.0) then
      allocate( soln%DE( i_low:i_high, neq ), &
                soln%DEnorm( 1:neq )  )
      soln%DE = zero
      soln%DE = zero
    end if
    
    soln%V   = zero
    soln%U   = zero
    soln%R   = zero
    soln%F   = zero
    soln%D   = zero
    soln%asnd   = zero
    soln%mach   = zero
    soln%temp   = zero
    soln%src   = zero
    soln%dt  = zero
    soln%lambda = zero
    soln%rnorm = zero
    soln%rold = zero
    soln%rinit = zero

  end subroutine allocate_soln
  
  
  !=========================== deallocate_soln ===============================80
  !>
  !! Description: 
  !!
  !! Inputs:      soln : 
  !!
  !! Outputs:     soln : 
  !<
  !===========================================================================80
  subroutine deallocate_soln( soln )
  
    implicit none
    
    type(soln_t), intent(inout) :: soln
    
    deallocate(soln%V,      &
               soln%U,      &
               soln%R,      &
               soln%F,      &
               soln%D,      &
               soln%asnd,      &
               soln%mach,      &
               soln%temp,      &
               soln%src,      &
               soln%dt,     &
               soln%lambda, &
               soln%rnorm,  &
               soln%rold,  &
               soln%rinit  )
    
    if (shock.eq.0) then
      deallocate( soln%DE, soln%DEnorm )
    end if
    
  end subroutine deallocate_soln

end module soln_type
