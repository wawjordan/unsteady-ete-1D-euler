module riemann_problem

  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, half
  use set_inputs,    only : neq
  use fluid_constants, only : gamma
  use quadrature, only : gauss_integ
  use variable_conversion, only : speed_of_sound
  use utilities, only : newton_safe, hunt

  implicit none

  private

  public :: solve_riemann, sample_riemann, cell_avg_riemann, init_bracket

contains

subroutine pressure_function( pres, U_K, f_out )
    real(prec), intent(in)  :: pres
    real(prec), intent(in)  :: U_K(:)
    real(prec), intent(out) :: f_out
    real(prec) :: alpha, rho_K, pres_K, avel_K, A_K, B_K, &
                            f1, f2, df1, df2, msk1, msk2

    alpha = (gamma - one)/(gamma + one)

    rho_K  = U_K(1)
    pres_K = U_K(3)
    call speed_of_sound( pres_K, rho_K, avel_K)

    A_K = two/( (gamma + one)*rho_K )
    B_K = alpha*pres_K

    f1 = (pres - pres_K)*sqrt( A_K/(pres + B_K) ) !shock
    f2 = two*avel_K/(gamma - one)*( (pres/pres_K) &
           **( (gamma - one)/(two*gamma) ) - one) !expansion

    msk1 = merge(one,zero,(pres>pres_K))
    msk2 = merge(one,zero,(pres<=pres_K))

    f_out = f1*msk1 + f2*msk2
  end subroutine pressure_function

  subroutine pressure_derivative_function(pres,U_K,df_out)
    real(prec), intent(in)  :: pres
    real(prec), intent(in)  :: U_K(:)
    real(prec), intent(out) :: df_out
    real(prec) :: alpha, rho_K, pres_K, avel_K, A_K, B_K, &
                            f1, f2, df1, df2, msk1, msk2
    alpha = (gamma - one)/(gamma + one)

    rho_K  = U_K(1)
    pres_K = U_K(3)
    call speed_of_sound( pres_K, rho_K, avel_K )

    A_K = two/( (gamma + one)*rho_K )
    B_K = alpha*pres_K

    df1 = sqrt( A_K/(B_K + pres) )*( one - &
              ( (pres - pres_K)/(two*(B_K + pres) ) ) )
    df2 = (one/(rho_K*avel_K) )*(pres/pres_K)**(-(gamma + one)/(two*gamma) )

    msk1 = merge(one,zero,(pres>pres_K))
    msk2 = merge(one,zero,(pres<=pres_K))

    df_out = df1*msk1 + df2*msk2
  end subroutine pressure_derivative_function

  subroutine init_bracket( U_L, U_R, tol, p0, bnd_a, bnd_b )
    real(prec), intent(in)  :: U_L(:), U_R(:)
    real(prec), intent(in)  :: tol
    real(prec), intent(out) :: p0, bnd_a, bnd_b
    real(prec) :: rho_L, rho_R, pres_L, pres_R, delta_u, avel_L, avel_R,     &
                  pmin, pmax, fmin, fmax, A_L, B_L, A_R, B_R, pTR, pPV, pTS, &
                  gL, gR, p01, p02, p03, p04

    rho_L   = U_L(1)
    rho_R   = U_R(1)
    delta_u = U_R(2) - U_L(2)
    pres_L  = U_L(3)
    pres_R  = U_R(3)

    call speed_of_sound(pres_L, rho_L, avel_L)
    call speed_of_sound(pres_R, rho_R, avel_R)

    pmin = min(pres_L,pres_R)
    pmax = max(pres_L,pres_R)

    call pressure_function(pmin,U_L,fmin)
    call pressure_function(pmax,U_R,fmax)

!    if ( (fmin>zero).and.(fmax>zero) ) then
!      bnd_a = zero
!      bnd_b = pmin
!    elseif ( (fmin<=zero).and.(fmax>=zero) ) then
!      bnd_a = pmin
!      bnd_b = pmax
!    elseif ( (fmin<zero).and.(fmax<zero) ) then
!      bnd_a = pmax
!      !bnd_b = large
!      bnd_b = 1.0e99_prec
!    elseif &
!      ( (two*avel_L/(gamma - one) + two*avel_R/(gamma - one) ) > delta_u ) &
!        then
!      write(*,*) 'Initial data violates pressure positivity condition'
!    end if
    bnd_a = zero
    bnd_b = 1.0e99_prec
    A_L = two/( (gamma + one)*rho_L )
    B_L = (gamma - one)/(gamma + one)*pres_L

    A_R = two/( (gamma + one)*rho_R )
    B_R = (gamma - one)/(gamma + one)*pres_R

    pTR = ( ( avel_L + avel_R - half*(gamma - one)*delta_u ) / &
          ( (avel_L/pres_L)**( (gamma - one)/(two*gamma) ) +   &
            (avel_R/pres_R)**( (gamma - one)/(two*gamma) ) ) ) &
            **( two*gamma/(gamma - one) )
    p01 = pTR

    pPV = half*(pres_L + pres_R) &
            - 0.125_prec*delta_u*(rho_L + rho_R)*(avel_L + avel_R)
    p02 = max(tol,pPV)

    gL = sqrt( A_L/(p02 + B_L) )
    gR = sqrt( A_R/(p02 + B_R) )
    pTS = (gL*pres_L + gR*pres_R - delta_u)/(gL + gR)
    p03 = max(tol,pTS)

    p04 = half*(pres_L + pres_R)

    !p0 = fourth*(p01 + p02 + p03 + p04)
    p0 = p03
  end subroutine init_bracket

  subroutine solve_riemann(p0,a,b,tol,max_iter,U_L,U_R,pStar,uStar)

    real(prec), intent(in)  :: p0, a, b, tol
    real(prec), intent(in)  :: U_L(:), U_R(:)
    real(prec), intent(out) :: pStar, uStar
    integer, intent(in) :: max_iter
    real(prec) :: delta_u, fL, fR, err

    delta_u = U_R(2) - U_L(2)

    call newton_safe( f, df, p0, a, b, tol, max_iter, pStar, err)

    call pressure_function(pStar,U_L,fL)
    call pressure_function(pStar,U_R,fR)
    uStar = half*(U_L(2) + U_R(2)) + half*(fR - fL)

  contains
    function f(pres)
      real(prec) :: f, f1, f2
      real(prec), intent(in) :: pres
      call pressure_function(pres,U_L,f1)
      call pressure_function(pres,U_R,f2)
      f = f1 + f2 + delta_u
!      write(*,*) pres
!      write(*,*) U_L
!      write(*,*) U_R
!      write(*,*) f1
!      write(*,*) f2
!      write(*,*) f
    end function f

    function df(pres)
      real(prec) :: df, df1, df2
      real(prec), intent(in) :: pres
      call pressure_derivative_function(pres,U_L,df1)
      call pressure_derivative_function(pres,U_R,df2)
      df = df1 + df2
    end function df

  end subroutine solve_riemann

  subroutine sample_riemann(uStar,pStar,U_L,U_R,x,t,tol,U,xlocs,xmsk)
    real(prec), intent(in)  :: uStar, pStar, t, tol
    real(prec), intent(in)  :: U_L(:), U_R(:)
    real(prec), intent(in)  :: x(:)
    real(prec), intent(inout) :: U(:,:)
    real(prec), intent(out)   :: xlocs(5)
    real(prec), dimension(lbound(U,1):ubound(U,1)) :: S
    integer, dimension(5) :: xmsk
    logical, dimension(lbound(U,1):ubound(U,1)) :: &
    msk_L1, msk_L2, msk_L3, msk_L4, msk_R1, msk_R2, msk_R3, msk_R4

    !logical, allocatable    :: mask1(:)
    !real(prec), allocatable :: S(:)
    real(prec) :: U_star_shk_L(neq), U_star_shk_R(neq)
    real(prec) :: U_star_fan_L(neq), U_star_fan_R(neq)
    real(prec) :: rho_L, rho_R, uvel_L, uvel_R, pres_L, pres_R, alpha,     &
                 avel_L, avel_R, shkvel_L, shkvel_R, rho_shk_L, rho_shk_R, &
                 rho_fan_L, rho_fan_R, aStarL, aStarR, S_HL, S_TL, S_HR, S_TR, &
                 t2, xoffset
    integer :: i, j, k, n, mL, mR


    n = size(x,1)
    t2 = t
    t2 = max(tol,t2)
    xoffset = zero !half*(maxval(x)-minval(x))
    S = x/t2

    xmsk  = 0
    xlocs = zero

    rho_L  = U_L(1)
    rho_R  = U_R(1)
    uvel_L = U_L(2)
    uvel_R = U_R(2)
    pres_L = U_L(3)
    pres_R = U_R(3)

    alpha = (gamma - one)/(gamma + one)

    call speed_of_sound(pres_L,rho_L,avel_L)
    call speed_of_sound(pres_R,rho_R,avel_R)

    shkvel_L = uvel_L - avel_L*sqrt( ( (gamma + one)/(two*gamma) )* &
                   (pStar/pres_L) + ( (gamma - one)/(two*gamma) ) )
    shkvel_R = uvel_R + avel_R*sqrt( ( (gamma + one)/(two*gamma) )* &
                   (pStar/pres_R) + ( (gamma - one)/(two*gamma) ) )
    rho_shk_L = rho_L*( (pStar/pres_L) + alpha )/(alpha*(pStar/pres_L) + one)
    rho_shk_R = rho_R*( (pStar/pres_R) + alpha )/(alpha*(pStar/pres_R) + one)

    rho_fan_L = rho_L*(pStar/pres_L)**(one/gamma)
    rho_fan_R = rho_R*(pStar/pres_R)**(one/gamma)

    aStarL = avel_L*(pStar/pres_L)**( (gamma-one)/(two*gamma) )
    aStarR = avel_R*(pStar/pres_R)**( (gamma-one)/(two*gamma) )

    S_HL = uvel_L - avel_L
    S_TL = uStar - aStarL

    S_HR = uvel_R + avel_R
    S_TR = uStar + aStarR

    U_Star_shk_L = (/ rho_shk_L, uStar, pStar /)
    U_Star_shk_R = (/ rho_shk_R, uStar, pStar /)

    U_Star_fan_L = (/ rho_fan_L, uStar, pStar /)
    U_Star_fan_R = (/ rho_fan_R, uStar, pStar /)

    msk_L1 = (S < uStar)
    msk_L2 = (S < shkvel_L)
    msk_L3 = (S < S_HL)
    msk_L4 = (S > S_TL)

    msk_R1 = .not.(msk_L1)
    msk_R2 = (S > shkvel_R)
    msk_R3 = (S > S_HR)
    msk_R4 = (S < S_TR)

    xmsk(3)  = 2
    xlocs(3) = ustar*t + xoffset

    if (pStar > U_L(3)) then
      xlocs(2) = shkvel_L*t + xoffset
      xmsk(2)  = 2
      forall (i=1:n,(msk_L1(i).and.msk_L2(i)) ) U(i,:) = U_L
      forall (i=1:n,(msk_L1(i).and.(.not.msk_L2(i)) ) ) U(i,:) = U_Star_shk_L
    else
      xlocs(1) = S_HL*t + xoffset
      xlocs(2) = S_TL*t + xoffset
      xmsk(1)  = 1
      xmsk(2)  = 1
      forall (i=1:n,(msk_L1(i).and.msk_L3(i)) ) U(i,:) = U_L
      forall (i=1:n,(msk_L1(i).and.msk_L4(i)) ) U(i,:) = U_Star_fan_L
      forall (i=1:n,(msk_L1(i).and.( &
             (.not.msk_L3(i)).and.(.not.msk_L4(i) ) ) ) ) &
           U(i,:) = (/ rho_in_fan_L(rho_L,uvel_L,avel_L,x(i),t), &
                     uvel_in_fan_L(uvel_L,avel_L,x(i),t), &
                      pres_in_fan_L(pres_L,uvel_L,avel_L,x(i),t) /)
    end if
    if (pStar > U_R(3)) then
      xlocs(4) = shkvel_R*t + xoffset
      xmsk(4)  = 2
      forall (i=1:n,(msk_R1(i).and.msk_R2(i)) ) U(i,:) = U_R
      forall (i=1:n,(msk_R1(i).and.(.not.msk_R2(i)) ) ) U(i,:) = U_Star_shk_R
    else
      xlocs(4) = S_TR*t + xoffset
      xlocs(5) = S_HR*t + xoffset
      xmsk(4)  = 1
      xmsk(5)  = 1
      forall (i=1:n,(msk_R1(i).and.msk_R3(i)) ) U(i,:) = U_R
      forall (i=1:n,(msk_R1(i).and.msk_R4(i)) ) U(i,:) = U_Star_fan_R
      forall (i=1:n,(msk_R1(i).and.( &
             (.not.msk_R3(i)).and.(.not.msk_R4(i) ) ) ) ) &
           U(i,:) = (/ rho_in_fan_R(rho_R,uvel_R,avel_R,x(i),t), &
                      uvel_in_fan_R(uvel_R,avel_R,x(i),t), &
                      pres_in_fan_R(pres_R,uvel_R,avel_R,x(i),t) /)
    end if

  end subroutine sample_riemann

  subroutine cell_avg_riemann(uStar,pStar,U_L,U_R,x,t,tol,U,xlocs,xmsk)
    real(prec), intent(in)  :: uStar, pStar, t, tol
    real(prec), intent(in)  :: U_L(:), U_R(:)
    real(prec), intent(in)  :: x(:)
    real(prec), intent(inout) :: U(:,:)
    real(prec), intent(out)   :: xlocs(5)
    real(prec), dimension(lbound(x,1):ubound(x,1),lbound(U,2):ubound(U,2)) :: IU
    real(prec), dimension(lbound(x,1):ubound(x,1)) :: S
    real(prec) :: U_L2(5,neq), U_R2(5,neq)
    integer, dimension(5) :: xmsk, ilocs
    logical, dimension(lbound(x,1):ubound(x,1)) :: &
    msk_L1, msk_L2, msk_L3, msk_L4, msk_R1, msk_R2, msk_R3, msk_R4

    !logical, allocatable    :: mask1(:)
    !real(prec), allocatable :: S(:)
    real(prec) :: U_star_shk_L(neq), U_star_shk_R(neq), U_star_L(neq)
    real(prec) :: U_star_fan_L(neq), U_star_fan_R(neq), U_star_R(neq)
    real(prec) :: rho_L, rho_R, uvel_L, uvel_R, pres_L, pres_R, alpha,     &
                 avel_L, avel_R, shkvel_L, shkvel_R, rho_shk_L, rho_shk_R, &
                 rho_fan_L, rho_fan_R, aStarL, aStarR, S_HL, S_TL, S_HR, S_TR, &
                 t2, xoffset
    integer :: i, j, k, n, mL, mR, lbnd, hbnd


    n = size(x,1)
    t2 = t
    t2 = max(tol,t2)
    xoffset = zero !half*(maxval(x)-minval(x))
    S = x/t2
    IU = zero

    xmsk  = 0
    xlocs = zero
    ilocs = n/2

    rho_L  = U_L(1)
    rho_R  = U_R(1)
    uvel_L = U_L(2)
    uvel_R = U_R(2)
    pres_L = U_L(3)
    pres_R = U_R(3)

    alpha = (gamma - one)/(gamma + one)

    call speed_of_sound(pres_L,rho_L,avel_L)
    call speed_of_sound(pres_R,rho_R,avel_R)

    shkvel_L = uvel_L - avel_L*sqrt( ( (gamma + one)/(two*gamma) )* &
                   (pStar/pres_L) + ( (gamma - one)/(two*gamma) ) )
    shkvel_R = uvel_R + avel_R*sqrt( ( (gamma + one)/(two*gamma) )* &
                   (pStar/pres_R) + ( (gamma - one)/(two*gamma) ) )
    rho_shk_L = rho_L*( (pStar/pres_L) + alpha )/(alpha*(pStar/pres_L) + one)
    rho_shk_R = rho_R*( (pStar/pres_R) + alpha )/(alpha*(pStar/pres_R) + one)

    rho_fan_L = rho_L*(pStar/pres_L)**(one/gamma)
    rho_fan_R = rho_R*(pStar/pres_R)**(one/gamma)

    aStarL = avel_L*(pStar/pres_L)**( (gamma-one)/(two*gamma) )
    aStarR = avel_R*(pStar/pres_R)**( (gamma-one)/(two*gamma) )

    S_HL = uvel_L - avel_L
    S_TL = uStar - aStarL

    S_HR = uvel_R + avel_R
    S_TR = uStar + aStarR

    U_Star_shk_L = (/ rho_shk_L, uStar, pStar /)
    U_Star_shk_R = (/ rho_shk_R, uStar, pStar /)

    U_Star_fan_L = (/ rho_fan_L, uStar, pStar /)
    U_Star_fan_R = (/ rho_fan_R, uStar, pStar /)

    msk_L1 = (S < uStar)
    msk_L2 = (S < shkvel_L)
    msk_L3 = (S < S_HL)
    msk_L4 = (S > S_TL)

    msk_R1 = .not.(msk_L1)
    msk_R2 = (S > shkvel_R)
    msk_R3 = (S > S_HR)
    msk_R4 = (S < S_TR)

    xmsk(3)  = 2
    xlocs(3) = ustar*t

    if (pStar > U_L(3)) then
      xlocs(2) = shkvel_L*t
      xmsk(2)  = 2
      U_star_L = U_Star_shk_L
      forall (i=1:n,(msk_L1(i).and.msk_L2(i)) ) IU(i,:) = U_L*x(i)
      forall (i=1:n,(msk_L1(i).and.(.not.msk_L2(i)) ) ) IU(i,:) = &
                                                         U_Star_shk_L*x(i)
    else
      xlocs(1) = S_HL*t
      xlocs(2) = S_TL*t
      xmsk(1)  = 1
      xmsk(2)  = 1
      U_Star_L = U_star_fan_L
      forall (i=1:n,(msk_L1(i).and.msk_L3(i)) ) IU(i,:) = U_L*x(i)
      forall (i=1:n,(msk_L1(i).and.msk_L4(i)) ) IU(i,:) = U_Star_fan_L*x(i)
      forall (i=1:n,(msk_L1(i).and.( &
             (.not.msk_L3(i)).and.(.not.msk_L4(i) ) ) ) ) &
           IU(i,:) = (/ int_rho_L(rho_L,uvel_L,avel_L,x(i),t), &
                        int_uvel_L(uvel_L,avel_L,x(i),t), &
                        int_pres_L(pres_L,uvel_L,avel_L,x(i),t) /)
    end if
    if (pStar > U_R(3)) then
      xlocs(4) = shkvel_R*t
      xmsk(4)  = 2
      U_Star_R = U_star_shk_R
      forall (i=1:n,(msk_R1(i).and.msk_R2(i)) ) IU(i,:) = U_R*x(i)
      forall (i=1:n,(msk_R1(i).and.(.not.msk_R2(i)) ) ) IU(i,:) = &
                                                         U_Star_shk_R*x(i)
    else
      xlocs(4) = S_TR*t
      xlocs(5) = S_HR*t
      xmsk(4)  = 1
      xmsk(5)  = 1
      U_Star_R = U_Star_fan_R
      forall (i=1:n,(msk_R1(i).and.msk_R3(i)) ) IU(i,:) = U_R*x(i)
      forall (i=1:n,(msk_R1(i).and.msk_R4(i)) ) IU(i,:) = U_Star_fan_R*x(i)
      forall (i=1:n,(msk_R1(i).and.( &
             (.not.msk_R3(i)).and.(.not.msk_R4(i) ) ) ) ) &
             IU(i,:) = (/ int_rho_R(rho_R,uvel_R,avel_R,x(i),t), &
                          int_uvel_R(uvel_R,avel_R,x(i),t), &
                          int_pres_R(pres_R,uvel_R,avel_R,x(i),t) /)
    end if

    do i = 1,5
      call hunt(x,xlocs(i),ilocs(i))
    end do

    write(*,*) ilocs

    do i = 1,5
      if (xmsk(i) == 0) then  ! no feature detected (left)
        U_L2(i,:) = U_L*xlocs(i)
        U_R2(i,:) = U_R*xlocs(i)
      elseif (i == 1 .and. xmsk(i) == 1) then ! left expansion (head)
        U_L2(i,:) = U_L*xlocs(i)
        U_R2(i,:) = (/ int_rho_L(rho_L,uvel_L,avel_L,xlocs(i),t), &
                      int_uvel_L(uvel_L,avel_L,xlocs(i),t), &
                      int_pres_L(pres_L,uvel_L,avel_L,xlocs(i),t) /)
      elseif (i == 2 .and. xmsk(i) == 1) then ! left expansion (tail)
        U_L2(i,:) = (/ int_rho_L(rho_L,uvel_L,avel_L,xlocs(i),t), &
                      int_uvel_L(uvel_L,avel_L,xlocs(i),t), &
                      int_pres_L(pres_L,uvel_L,avel_L,xlocs(i),t) /)
        U_R2(i,:) = U_Star_fan_L*xlocs(i)
      elseif (i == 2 .and. xmsk(i) == 2) then ! left shock
        U_L2(i,:) = U_L*xlocs(i)
        U_R2(i,:) = U_Star_shk_L*xlocs(i)
      elseif (i == 3) then ! contact discontinuity
        U_L2(i,:) = U_Star_L*xlocs(i)
        U_R2(i,:) = U_Star_R*xlocs(i)
      elseif (i == 4 .and. xmsk(i) == 2) then ! right shock
        U_L2(i,:) = U_star_shk_R*xlocs(i)
        U_R2(i,:) = U_R*xlocs(i)
      elseif (i == 4 .and. xmsk(i) == 1) then ! right expansion (tail)
        U_L2(i,:) = U_Star_fan_R*xlocs(i)
        U_R2(i,:) = (/ int_rho_R(rho_R,uvel_R,avel_R,xlocs(i),t), &
                      int_uvel_R(uvel_R,avel_R,xlocs(i),t), &
                      int_pres_R(pres_R,uvel_R,avel_R,xlocs(i),t) /)
      elseif (i == 5 .and. xmsk(i) == 1) then ! right expansion (head)
        U_L2(i,:) = (/ int_rho_R(rho_R,uvel_R,avel_R,xlocs(i),t), &
                      int_uvel_R(uvel_R,avel_R,xlocs(i),t), &
                      int_pres_R(pres_R,uvel_R,avel_R,xlocs(i),t) /)
        U_R2(i,:) = U_R*xlocs(i)
      end if
    end do


    !write(*,*) lbound(U,1), ubound(U,1)
    !write(*,*) lbound(IU,1), ubound(IU,1)
    !write(*,*) lbound(x,1), ubound(x,1)

    !stop

    do i = 1,n-1
      U(i,:) = ( IU(i+1,:) - IU(i,:) )/( x(i+1) - x(i) )
    end do
    !do i = 1,n-1
    !  U(i,:) = ( IU(i+1,:) - IU(i,:) )/( x(i+1) - x(i) )
    !end do

    do i = 1,5
      U(ilocs(i),:) = ( U_L2(i,:) - IU(ilocs(i),:)  + &
                      IU(ilocs(i)+1,:) - U_R2(i,:) )/ &
                        ( x(ilocs(i)+1) - x(ilocs(i)) )
    end do

    xlocs = xlocs + xoffset

  end subroutine cell_avg_riemann

  elemental function rho_in_fan_L(rho_L,uvel_L,avel_L,x,t)
    real(prec) :: rho_in_fan_L
    real(prec), intent(in) :: rho_L, uvel_L, avel_L, x, t
    real(prec) :: alpha
    alpha  = (gamma - one)/(gamma + one)
    rho_in_fan_L = rho_L*( two/(gamma + one) &
                 + alpha/avel_L*(uvel_L - x/t) )**( two/(gamma - one) )
  end function rho_in_fan_L

  elemental function uvel_in_fan_L(uvel_L,avel_L,x,t)
    real(prec) :: uvel_in_fan_L
    real(prec), intent(in) :: uvel_L, avel_L, x, t
    real(prec) :: alpha
    alpha  = (gamma - one)/(gamma + one)
    uvel_in_fan_L = two/(gamma + one)* &
                    (avel_L + half*(gamma - one)*uvel_L + x/t)
  end function uvel_in_fan_L

  elemental function pres_in_fan_L(pres_L,uvel_L,avel_L,x,t)
    real(prec) :: pres_in_fan_L
    real(prec), intent(in) :: pres_L, uvel_L, avel_L, x, t
    real(prec) :: alpha
    alpha  = (gamma - one)/(gamma + one)
    pres_in_fan_L = pres_L*( two/(gamma + one) &
                  + alpha/avel_L*(uvel_L - x/t) )**( two*gamma/(gamma - one) )
  end function pres_in_fan_L

  elemental function rho_in_fan_R(rho_R,uvel_R,avel_R,x,t)
    real(prec) :: rho_in_fan_R
    real(prec), intent(in) :: rho_R, uvel_R, avel_R, x, t
    real(prec) :: alpha
    alpha  = (gamma - one)/(gamma + one)
    rho_in_fan_R = rho_R*( two/(gamma + one) &
                 - alpha/avel_R*(uvel_R - x/t) )**( two/(gamma - one) )
  end function rho_in_fan_R

  elemental function uvel_in_fan_R(uvel_R,avel_R,x,t)
    real(prec) :: uvel_in_fan_R
    real(prec), intent(in) :: uvel_R, avel_R, x, t
    real(prec) :: alpha
    alpha  = (gamma - one)/(gamma + one)
    uvel_in_fan_R = two/(gamma + one)* &
                    (-avel_R + half*(gamma - one)*uvel_R + x/t)
  end function uvel_in_fan_R

  elemental function pres_in_fan_R(pres_R,uvel_R,avel_R,x,t)
    real(prec) :: pres_in_fan_R
    real(prec), intent(in) :: pres_R, uvel_R, avel_R, x, t
    real(prec) :: alpha
    alpha  = (gamma - one)/(gamma + one)
    pres_in_fan_R = pres_R*( two/(gamma + one) &
                  - alpha/avel_R*(uvel_R - x/t) )**( two*gamma/(gamma - one) )
  end function pres_in_fan_R

!========================================================================
  elemental function int_rho_L(rho_L,uvel_L,avel_L,x,t)
    real(prec) :: int_rho_L
    real(prec), intent(in) :: rho_L, uvel_L, avel_L, x, t

    int_rho_L = -( rho_L*t*( two/(gamma + one) + &
                 ( (uvel_L - x/t)*(gamma - one) )/ &
                 ( avel_L*(gamma + one) ) )** &
                 ( two/(gamma - one) )* &
                 ( two*avel_L - uvel_L + gamma*(uvel_L - x/t) + x/t ) )/ &
                 (gamma + one)
  end function int_rho_L

  elemental function int_uvel_L(uvel_L,avel_L,x,t)
    real(prec) :: int_uvel_L
    real(prec), intent(in) :: uvel_L, avel_L, x, t

    int_uvel_L = ( t*( avel_L + x/t + half*uvel_L*(gamma - one) )**2 )/ &
                 (gamma + one)
  end function int_uvel_L

  elemental function int_pres_L(pres_L,uvel_L,avel_L,x,t)
    real(prec) :: int_pres_L
    real(prec), intent(in) :: pres_L, uvel_L, avel_L, x, t

    int_pres_L = t*( two/(gamma + one) + ( (uvel_L - x/t)*(gamma - one) )/ &
                   ( avel_L*(gamma + one) ) )** &
                   ( (two*gamma)/(gamma - one) )* &
                   ( ( (pres_L - gamma*pres_L)*(uvel_L - x/t) )/&
                   (three*gamma - one) - (two*avel_L*pres_L)/ &
                   (three*gamma - one) )
  end function int_pres_L

  elemental function int_rho_R(rho_R,uvel_R,avel_R,x,t)
    real(prec) :: int_rho_R
    real(prec), intent(in) :: rho_R, uvel_R, avel_R, x, t

    int_rho_R = ( rho_R*t*( two/(gamma + one) - &
                ( (uvel_R - x/t)*(gamma - one) )/ &
                ( avel_R*(gamma + one) ) )** &
                ( two/(gamma - one) )* &
                ( two*avel_R + uvel_R - gamma*(uvel_R - x/t) - x/t ) )/ &
                (gamma + one)
  end function int_rho_R

  elemental function int_uvel_R(uvel_R,avel_R,x,t)
    real(prec) :: int_uvel_R
    real(prec), intent(in) :: uvel_R, avel_R, x, t

    int_uvel_R = ( t*( x/t - avel_R + half*uvel_R*(gamma - one) )**2 )/ &
                 (gamma + one)
  end function int_uvel_R

  elemental function int_pres_R(pres_R,uvel_R,avel_R,x,t)
    real(prec) :: int_pres_R
    real(prec), intent(in) :: pres_R, uvel_R, avel_R, x, t

    int_pres_R = t*( two/(gamma + one) - ( (uvel_R - x/t)*(gamma - one) )/ &
                   ( avel_R*(gamma + one) ) )** &
                   ( (two*gamma)/(gamma - one) )* &
                   ( ( (pres_R - gamma*pres_R)*(uvel_R - x/t) )/ &
                   (three*gamma - one) + (two*avel_R*pres_R)/ &
                   (three*gamma - one) )
  end function int_pres_R

end module riemann_problem
