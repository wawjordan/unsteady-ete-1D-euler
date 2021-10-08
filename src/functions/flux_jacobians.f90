module flux_jacobians

  use set_precision,   only : prec
  use set_constants,   only : small, zero, one, two, three, eight, &
                              half, fourth, eighth
  use fluid_constants, only : gamma
  use set_inputs,      only : neq, i_low, i_high, ig_low, ig_high, epsM, kappaM
  use soln_type,       only : soln_t
  use variable_conversion, only : cons2prim_1d, speed_of_sound
  use flux_calc,       only : exact_flux_cons

  implicit none

  private

  public :: calc_vl_dfdu, flux_jac_cons1D, flux_jac_prim1D

  contains



!  subroutine VL_flux_jac(UL,UR,Uim1,Ui,Uip1,Uip2,psipm1,psip,psim,psimp1, &
!                    dpsipm1,dpsim,dpsip,dpsimp1, dfduim1,dfdui,dfduip1,dfduip2)
!    real(prec), intent(in), dimension(neq) :: UL, UR, Uim1, Ui, Uip1, Uip2, &
!                                              psipm1, psip, psim, psimp1, &
!                                              dpsipm1, dpsip, dpsim, dpsimp1
!    real(prec), intent(out), dimension(neq,neq) :: dfduim1, dfdui, dfduip1, &
!                                                   dfduip2
!    real(prec) :: fL(neq), fR(neq), dfduL(neq,neq), dfduR(neq,neq)
!    real(prec), dimension(neq) :: duLim1, duLi, duLip1, duRi, duRip1, duRip2
!    real(prec) :: rhoL, rhoR, uvelL, uvelR, presL, presR, avelL, avelR, ML, MR
!    integer :: i
!
!
!    call cons2prim_1D(UL(1),UL(2),UL(3),rhoL,uvelL,presL)
!    call cons2prim_1D(UR(1),UR(2),UR(3),rhoR,uvelR,presR)
!    call speed_of_sound(presL,rhoL,avelL)
!    call speed_of_sound(presR,rhoR,avelR)
!
!    ML = uvelL/avelL
!    MR = uvelR/avelR
!
!    dfduL = zero
!    dfduR = zero
!    dfduim1 = zero
!    dfdui   = zero
!    dfduip1 = zero
!    dfduip2 = zero
!
!    duLim1 = one
!    duLi   = one
!    duLip1 = one
!    duRi   = one
!    duRip1 = one
!    duRip2 = one
!
!
!
!    !call duL_MUSCL_im1(Uim1,Ui,Uip1,psipm1,dpsipm1,dpsim,duLim1)
!    !call duL_MUSCL_i(Uim1,Ui,Uip1,psipm1,psim,dpsipm1,dpsim,duLi)
!    !call duL_MUSCL_ip1(Uim1,Ui,Uip1,psim,dpsipm1,dpsim,duLip1)
!    !call duR_MUSCL_i(Ui,Uip1,Uip2,psip,dpsip,dpsimp1,duRi)
!    !call duR_MUSCL_ip1(Ui,Uip1,Uip2,psip,psimp1,dpsip,dpsimp1,duRip1)
!    !call duR_MUSCL_ip2(Ui,Uip1,Uip2,psimp1,dpsip,dpsimp1,duRip2)
!
!    call calc_vl_dfdu(UL,fL,dfduL,.true.)
!    call calc_vl_dfdu(UR,fR,dfduR,.false.)
!
!    do i = 1,neq
!      dfduim1(:,i) = dfduL(:,i)*duLim1(i)
!    end do
!    do i = 1,neq
!      dfdui(:,i) = dfduL(:,i)*duLi(i) + dfduR(:,i)*duRi(i)
!    end do
!    do i = 1,neq
!      dfduip1(:,i) = dfduL(:,i)*duLip1(i) + dfduR(:,i)*duRip1(i)
!    end do
!    do i = 1,neq
!      dfduip2(:,i) = dfduR(:,i)*duRip2(i)
!    end do
!
!  end subroutine VL_flux_jac


  subroutine calc_vl_dfdu(U,f,dfdu,mask)

    real(prec), intent(in) :: U(neq)
    real(prec), intent(out) :: f(neq), dfdu(neq,neq)
    real(prec) :: f1(neq), dfdu1(neq,neq)
    logical, intent(in) :: mask ! .true./.false. = L/R
    real(prec) :: pm, g1, gg1, gg2, isgg1, phi, phi32, ggphi
    real(prec) :: P, uvel, rho, a, M, M1
    real(prec) :: U1, U2, U3
    real(prec) :: fa, fb, sub, sup
    real(prec) :: da1, da2, da3, db1, db2, db3
    logical :: supL, supR

    pm = merge(one,-one,mask)

    U1 = U(1)
    U2 = U(2)
    U3 = U(3)

    call cons2prim_1D(U1,U2,U3,rho,uvel,P)
    call speed_of_sound(P,rho,a)

    M = uvel/a
    M1 = M + pm

    g1  = gamma - one
    gg1 = gamma*(gamma-one)
    gg2 = gamma**2 - one
    isgg1 = one/sqrt(gg1)

    phi = U1*U3 - half*U2**2
    phi32 = one/(sqrt(phi)*phi)
    ggphi = sqrt(gg1/phi)

    fa = pm*fourth*U1*a*M1**2
    fb = pm*( g1*M*a + pm*two*a)

    f(1) = fa
    f(2) = fa*fb/gamma
    f(3) = half/gg2*fa*fb**2

    da1 = pm*( eighth*(M1**2)*U3*ggphi - fourth*U1*a*M1*(U2*U3)*isgg1*phi32)
    da2 = pm*(-eighth*(M1**2)*U2*ggphi +   half*U1*a*M1*(U1*U3)*isgg1*phi32)
    da3 = pm*( eighth*(M1**2)*U1*ggphi - fourth*U1*a*M1*(U1*U2)*isgg1*phi32)

    db1 = pm/U1*(pm*(U3*ggphi - two*a) - g1*U2/U1 )
    db2 = pm/U1*(g1 - pm*U2*ggphi)
    db3 = ggphi

    dfdu(1,1) = da1
    dfdu(2,1) = pm*(one/gamma)*(fb*da1 + fa*db1)
    dfdu(3,1) = (one/gg2)*(half*fb**2*da1 + fa*fb*db1)

    dfdu(1,2) = da2
    dfdu(2,2) = pm*(one/gamma)*(fb*da2 + fa*db2)
    dfdu(3,2) = (one/gg2)*(half*fb**2*da2 + fa*fb*db2)

    dfdu(1,3) = da3
    dfdu(2,3) = pm*(one/gamma)*(fb*da3 + fa*db3)
    dfdu(3,3) = (one/gg2)*(half*fb**2*da3 + fa*fb*db3)

    call exact_flux_cons(U,f1)
    call flux_jac_cons1D(U,dfdu1)

    sub  = merge( one, zero, ( abs(M)<one ) )
    supL = (M>=one).and.mask
    supR = (M<=-one).and.(.not.mask)
    sup  = merge( one, zero, ( supL.or.supR ) )

    f = f*sub + f1*sup
    dfdu = dfdu*sub + dfdu1*sup

  end subroutine calc_vl_dfdu

  subroutine flux_jac_cons1D(U,A)

    real(prec), intent(in)  :: U(3)
    real(prec), intent(out) :: A(3,3)

    A(1,1) = zero
    A(2,1) = half*(gamma-three)*(U(2)/U(1))**2
    A(3,1) = (gamma - one)*(U(2)/U(1))**3 - gamma*(U(2)*U(3))/U(1)
    A(1,2) = one
    A(2,2) = (three - gamma)*(U(2)/U(1))
    A(3,2) = gamma*(U(3)/U(1)) - half*three*(gamma-one)*(U(2)/U(1))**2
    A(1,3) = zero
    A(2,3) = gamma - one
    A(3,3) = gamma*(U(2)/U(1))

  end subroutine flux_jac_cons1D

  subroutine flux_jac_prim1D(V,A)

    real(prec), intent(in)  :: V(3)
    real(prec), intent(out) :: A(3,3)

    A(1,1) = V(2)
    A(2,1) = zero
    A(3,1) = zero
    A(1,2) = V(1)
    A(2,2) = V(2)
    A(3,2) = gamma*V(3)
    A(1,3) = zero
    A(2,3) = one/V(1)
    A(3,3) = V(2)

  end subroutine flux_jac_prim1D

end module flux_jacobians
