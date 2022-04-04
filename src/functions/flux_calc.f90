module flux_calc

  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, four, half, fourth
  use fluid_constants, only : gamma
  use set_inputs, only : neq, i_low, i_high, ig_low, ig_high, eps_roe
  use set_inputs, only : n_ghost
  use variable_conversion, only : cons2prim, speed_of_sound, limit_primitives
  use solution_reconstruction, only : MUSCL_extrap, MUSCL_extrap_lite
  use soln_type, only : soln_t
  use grid_type, only : grid_t

  implicit none

  private

  public :: flux_fun, select_flux, calc_flux_1D, exact_flux_cons

  procedure( calc_flux ), pointer :: flux_fun

  abstract interface

  !================================ calc_flux ================================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine calc_flux(left_state, right_state, F)

    import :: prec, neq
    real(prec), dimension(neq), intent(in) :: left_state, right_state
    real(prec), dimension(neq), intent(out) :: F

  end subroutine calc_flux

  end interface

contains

  !========================= calc_flux_1D_explicit ===========================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine calc_flux_1D(grid,soln)

    type(grid_t), intent(in)    :: grid
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq,i_low-1:i_high) :: Vleft, Vright, Uleft, Uright
    integer :: i

    call MUSCL_extrap_lite(soln,Vleft,Vright)
    call cons2prim(Uleft,Vleft)
    call cons2prim(Uright,Vright)

    do i = i_low-1,i_high
      call flux_fun(Uleft(:,i),Uright(:,i),soln%F(:,i))
    end do

end subroutine calc_flux_1D


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
      flux_fun => van_leer_flux
    case(2)
      flux_fun => roe_flux
    case default

    end select

  end subroutine select_flux

  subroutine exact_flux_cons(U,F)

    real(prec), dimension(neq), intent(in)  :: U
    real(prec), dimension(neq), intent(out) :: F

    F(1) = U(2)
    F(2) = U(2)**2/U(1) + (gamma-one)*( U(3) - half*U(2)**2/U(1) )
    F(3) = gamma*U(2)*U(3)/U(1) - half*(gamma-one)*U(2)**3/U(1)**2

  end subroutine exact_flux_cons

  !============================== van_leer_flux ==============================80
  !>
  !! Description: van Leer flux in conserved variables
  !<
  !===========================================================================80
  subroutine van_leer_flux(left, right, F)

    real(prec), dimension(neq), intent(in)  :: left, right
    real(prec), dimension(neq), intent(out) :: F
    real(prec), dimension(neq) :: VL, VR

    call cons2prim(left,VL)
    call cons2prim(right,VR)

    call van_leer_flux_prim(VL,VR,F)

  end subroutine van_leer_flux



  !============================ van_leer_flux_prim ===========================80
  !>
  !! Description: van Leer flux in primitive variables
  !<
  !===========================================================================80
  subroutine van_leer_flux_prim(left, right, F)

    real(prec), dimension(neq), intent(in)  :: left, right
    real(prec), dimension(neq), intent(out) :: F
    real(prec) ::   aL,   aR
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

    rhoL = left(1)
    uL   = left(2)
    pL   = left(3)
    aL   = speed_of_sound(pL,rhoL)

    rhoR = right(1)
    uR   = right(2)
    pR   = right(3)
    aR   = speed_of_sound(pR,rhoR)

    ML = uL/aL
    MR = uR/aR
    M_plus  =  fourth*(ML+one)**2
    M_minus = -fourth*(MR-one)**2
    beta_L = -max(zero,one-int(abs(ML)))
    beta_R = -max(zero,one-int(abs(MR)))
    alpha_plus  = half*(one+sign(one,ML))
    alpha_minus = half*(one-sign(one,MR))
    c_plus  = alpha_plus*(one+beta_L)*ML - beta_L*M_plus
    c_minus = alpha_minus*(one+beta_R)*MR - beta_R*M_minus
    htL = aL**2/(gamma-one) + half*uL**2
    htR = aR**2/(gamma-one) + half*uR**2

    Fc1 =  rhoL*aL*c_plus + rhoR*aR*c_minus
    Fc2 =  rhoL*aL*c_plus*uL + rhoR*aR*c_minus*uR
    Fc3 =  rhoL*aL*c_plus*htL + rhoR*aR*c_minus*htR

    p_plus  = M_plus*(-ML + two)
    p_minus = M_minus*(-MR - two)
    d_plus  = alpha_plus*(one+beta_L) - beta_L*p_plus
    d_minus = alpha_minus*(one+beta_R) - beta_R*p_minus

    Fp = d_plus*pL + d_minus*pR

    F(1) = Fc1
    F(2) = Fc2 + Fp
    F(3) = Fc3

  end subroutine van_leer_flux_prim

  !================================ roe_flux =================================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine roe_flux( left, right, F )

    real(prec), dimension(neq), intent(in)  :: left, right
    real(prec), dimension(neq), intent(out) :: F
    real(prec), dimension(neq) :: FL, FR, rvec1, rvec2, rvec3, lambda
    real(prec) ::  R
    real(prec) ::  aL, aR, a2
    real(prec) :: rhoL, rhoR, rho2
    real(prec) ::   uL,   uR,   u2
    real(prec) ::   pL,   pR
    real(prec) ::  htL,  htR,  ht2
    real(prec) :: dw1, dw2, dw3

    rhoL = left(1)
    uL   = left(2)
    pL   = left(3)
    aL   = speed_of_sound(pL,rhoL)

    rhoR = right(1)
    uR   = right(2)
    pR   = right(3)
    aR   = speed_of_sound(pR,rhoR)

    htL  = aL**2/(gamma-one) + half*uL**2
    htR  = aR**2/(gamma-one) + half*uR**2

    R = sqrt(rhoR/rhoL)
    rho2 = R*rhoL
    u2   = (R*uR + uL)/(R + one)
    ht2  = (R*htR + htL)/(R + one)
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

    F(:) = half*(FL+FR) - half*(lambda(1)*dw1*rvec1 + &
           & lambda(2)*dw2*rvec2 + lambda(3)*dw3*rvec3)
  end subroutine roe_flux

end module flux_calc
