module flux_calc
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, four, half, fourth
  use fluid_constants, only : gamma
  use set_inputs, only : neq, imax, i_low, i_high, ig_low, ig_high, eps_roe
  use variable_conversion, only : cons2prim, speed_of_sound
  
  implicit none
  
  private
  
  public :: flux_fun, select_flux
  
  procedure( calc_flux ), pointer :: flux_fun
   
  abstract interface
    
  !================================ calc_flux ================================80
  !>
  !! Description:
  !!
  !! Inputs:      left_state  :
  !!              right_state : 
  !!
  !! Outputs:     F     :
  !<
  !===========================================================================80
  subroutine calc_flux(left_state, right_state, F)
    
    import :: prec, i_low, i_high, neq
    real(prec), dimension(:,:), intent(in) :: left_state, right_state
    real(prec), dimension(i_low-1:i_high,neq), intent(out) :: F
    
  end subroutine calc_flux
    
  end interface
  
contains
  
  !================================ select_flux ==============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine select_flux()
    
    use set_inputs, only : flux_scheme
    
    flux_fun => null()
    
    select case(flux_scheme)
    
    case(1)
      flux_fun => central_flux
    case(2)
      flux_fun => van_leer_flux
    case(3)
      flux_fun => roe_flux
    case default
    
    end select
  
  end subroutine select_flux
  
  !================================ central_flux =============================80
  !>
  !! Description:
  !!
  !! Inputs:      left  :
  !!              right : 
  !!
  !! Outputs:     F     :
  !<
  !===========================================================================80
  subroutine central_flux(left, right, F)
    
    real(prec), dimension(:,:), intent(in) :: left, right
    real(prec), dimension(i_low-1:i_high,neq), intent(out)   :: F
    real(prec), dimension(i_low-1:i_high,neq) :: Ui
    
    
    Ui(i_low-1:i_high,:) = half*(left + right)
    
    F(:,1) = Ui(:,2) 
    F(:,2) = half*(three-gamma)*( Ui(:,2)**2 )/Ui(:,1) &
           + (gamma-one)*Ui(:,3)
    F(:,3) = Ui(:,3)*Ui(:,2)/Ui(:,1) &
           + Ui(:,2)/Ui(:,1)*( (gamma-one)*Ui(:,3) &
           - half*(gamma-one)*Ui(:,2)**2/Ui(:,1) )
  
  end subroutine central_flux
  
  !============================== van_leer_flux ==============================80
  !>
  !! Description:
  !!
  !! Inputs:      left  :
  !!              right : 
  !!
  !! Outputs:     F     :
  !<
  !===========================================================================80
  subroutine van_leer_flux(left, right, F)
    
    real(prec), dimension(:,:)               , intent(in)  :: left, right
    real(prec), dimension(i_low-1:i_high,neq), intent(out) :: F
    real(prec), dimension(i_low-1:i_high,neq) :: VL, VR
    real(prec), dimension(i_low-1:i_high)     ::  aL, aR
    real(prec) :: rhoL, rhoR
    real(prec) ::   uL,   uR
    real(prec) ::   pL,   pR
    real(prec) ::   ML,   MR
    real(prec) ::  htL,  htR
    real(prec) :: M_plus, M_minus
    real(prec) :: p_plus, p_minus
    real(prec) :: c_plus, c_minus
    real(prec) :: d_plus, d_minus
    real(prec) :: alpha_plus, alpha_minus
    real(prec) :: beta_L, beta_R
    real(prec) :: Fc1, Fc2, Fc3, Fp
    integer :: i
    
    call cons2prim(left,VL)
    call cons2prim(right,VR)
    call speed_of_sound(VL(:,3),VL(:,1),aL)
    call speed_of_sound(VR(:,3),VR(:,1),aR)
    
    do i = i_low-1,i_high
      rhoL = VL(i,1)
      rhoR = VR(i,1)
      uL   = VL(i,2)
      uR   = VR(i,2)
      pL   = VL(i,3)
      pR   = VR(i,3)
      
      ML = uL/aL(i)
      MR = uR/aR(i)
      M_plus  =  fourth*(ML+one)**2
      M_minus = -fourth*(MR-one)**2
      beta_L = -max(zero,one-int(abs(ML)))
      beta_R = -max(zero,one-int(abs(MR)))
      alpha_plus  = half*(one+sign(one,ML))
      alpha_minus = half*(one-sign(one,MR))
      c_plus  = alpha_plus*(one+beta_L)*ML - beta_L*M_plus
      c_minus = alpha_minus*(one+beta_R)*MR - beta_R*M_minus
      htL = aL(i)**2/(gamma-one) + half*uL**2
      htR = aR(i)**2/(gamma-one) + half*uR**2
      
      Fc1 = rhoL*aL(i)*c_plus + rhoR*aR(i)*c_minus
      Fc2 =  rhoL*aL(i)*c_plus*uL + rhoR*aR(i)*c_minus*uR
      Fc3 =  rhoL*aL(i)*c_plus*htL + rhoR*aR(i)*c_minus*htR
      
      p_plus  = M_plus*(-ML + two)
      p_minus = M_minus*(-MR - two)
      d_plus  = alpha_plus*(one+beta_L) - beta_L*p_plus
      d_minus = alpha_minus*(one+beta_R) - beta_R*p_minus
      
      Fp = d_plus*pL + d_minus*pR
      
      F(i,1) = Fc1
      F(i,2) = Fc2 + Fp
      F(i,3) = Fc3
    end do
    
  end subroutine van_leer_flux
  
  !================================ roe_flux =================================80
  !>
  !! Description:
  !!
  !! Inputs:      left  :
  !!              right : 
  !!
  !! Outputs:     F     :
  !<
  !===========================================================================80
  subroutine roe_flux( left, right, F )
    
    real(prec), dimension(:,:)               , intent(in)  :: left, right
    real(prec), dimension(i_low-1:i_high,neq), intent(out) :: F
    real(prec), dimension(i_low-1:i_high,neq) :: VL, VR
    real(prec), dimension(i_low-1:i_high)     ::  aL, aR
    real(prec), dimension(i_low-1:i_high)     ::  R
    real(prec), dimension(3) :: FL, FR, rvec1, rvec2, rvec3, lambda
    real(prec) :: rhoL, rhoR, rho2
    real(prec) ::   uL,   uR,   u2
    real(prec) ::   pL,   pR,   a2
    real(prec) ::  htL,  htR,  ht2
    real(prec) :: dw1, dw2, dw3
    integer :: i
    
    call cons2prim(left,VL)
    call cons2prim(right,VR)
    call speed_of_sound(VL(:,3),VL(:,1),aL)
    call speed_of_sound(VR(:,3),VR(:,1),aR)
    
    do i = i_low-1,i_high
      rhoL = VL(i,1)
      rhoR = VR(i,1)
      uL   = VL(i,2)
      uR   = VR(i,2)
      pL   = VL(i,3)
      pR   = VR(i,3)
      htL  = aL(i)**2/(gamma-one) + half*uL**2
      htR  = aR(i)**2/(gamma-one) + half*uR**2
      
      R(i) = sqrt(rhoR/rhoL)
      rho2 = R(i)*rhoL
      u2   = (R(i)*uR + uL)/(R(i) + one)
      ht2  = (R(i)*htR + htL)/(R(i) + one)
      a2   = sqrt((gamma-one)*(ht2-half*u2**2))
      
      lambda = (/ u2, u2 + a2, u2 - a2 /)
      
      rvec1 = (/ one, lambda(1), half*u2**2 /)
      rvec2 = half*(rho2/a2)*(/ one, lambda(2), ht2 + u2*a2 /)
      rvec3 = -half*(rho2/a2)*(/ one, lambda(3), ht2 - u2*a2 /)
      
      lambda = abs(lambda)
      lambda = half*(one+sign(one,lambda-two*eps_roe*a2))*lambda &
           & + half*(one-sign(one,lambda-two*eps_roe*a2))*&
           & (lambda**2/(four*eps_roe*a2) + eps_roe*a2)
      
      dw1 = (rhoR - rhoL) - (pR - pL)/a2**2
      dw2 = (uR - uL) + (pR - pL)/(rho2*a2)
      dw3 = (uR - uL) - (pR - pL)/(rho2*a2)
      
      FL = (/ rhoL*uL, rhoL*uL**2 + pL, rhoL*uL*htL /)
      FR = (/ rhoR*uR, rhoR*uR**2 + pR, rhoR*uR*htR /)
      
      F(i,:) = half*(FL+FR) - half*(lambda(1)*dw1*rvec1 + &
             & lambda(2)*dw2*rvec2 + lambda(3)*dw3*rvec3)
    end do
  end subroutine roe_flux
  
end module flux_calc
