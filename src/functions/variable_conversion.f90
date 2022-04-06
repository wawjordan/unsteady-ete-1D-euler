module variable_conversion

  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, three, eight, half, fourth, eighth
  use fluid_constants, only : gamma, R_gas
  use set_inputs,      only : p0,T0, neq, ig_low, ig_high
  use soln_type,       only : soln_t

  implicit none

  private

  public :: riemann_1, riemann_2, riemann_3, riem2prim, &
            update_states, update_mach, speed_of_sound, &
            prim2cons, cons2prim, limit_primitives,     &
            isentropic_relations, cons2prim_1D

  interface cons2prim
    module procedure cons2prim_scalar
    module procedure cons2prim_array
  end interface

  interface prim2cons
    module procedure prim2cons_scalar
    module procedure prim2cons_array
  end interface

  interface update_mach
    module procedure update_mach_scalar
    module procedure update_mach_array
  end interface

  contains

  !============================== speed_of_sound  ============================80
  !>
  !! Description:
  !<
  !===========================================================================80
!  elemental subroutine speed_of_sound( pressure, rho, sound_speed )
!
!    real(prec), intent(in)  :: pressure, rho
!    real(prec), intent(out) :: sound_speed
!
!    sound_speed = sqrt(gamma*pressure/rho)
!
!  end subroutine speed_of_sound
  elemental function speed_of_sound( pressure, rho )

    real(prec) :: speed_of_sound
    real(prec), intent(in) :: pressure, rho

    speed_of_sound = sqrt(gamma*pressure/rho)

  end function speed_of_sound

  !=============================== cons2prim =================================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine cons2prim_scalar( U, V )

    real(prec), dimension(:), intent(in) :: U
    real(prec), dimension(:), intent(out) :: V

    V(1) = U(1)
    V(2) = U(2)/U(1)
    V(3) = (gamma - one)*U(3) - half*(gamma - one)*U(2)**2/U(1)

  end subroutine cons2prim_scalar
  subroutine cons2prim_array( U, V )

    real(prec), dimension(:,:), intent(in) :: U
    real(prec), dimension(:,:), intent(out) :: V
    !integer :: i,j
    !do i = lbound(U,2),ubound(U,2)
    !  write(*,*) (U(j,i),j=1,neq)
    !end do

    V(1,:) = U(1,:)
    V(2,:) = U(2,:)/U(1,:)
    V(3,:) = (gamma - one)*U(3,:) - half*(gamma - one)*U(2,:)**2/U(1,:)

  end subroutine cons2prim_array
  elemental subroutine cons2prim_1D( U1, U2, U3, V1, V2, V3 )
    real(prec), intent(in)  :: U1, U2, U3
    real(prec), intent(out) :: V1, V2, V3

    V1 = U1
    V2 = U2/U1
    V3 = (gamma - one)*U3 - half*(gamma - one)*U2**2/U1
  end subroutine cons2prim_1D

  !=============================== prim2cons =================================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine prim2cons_scalar( U, V )

    real(prec), dimension(:), intent(out)  :: U
    real(prec), dimension(:), intent(in)   :: V

    U(1) = V(1)
    U(2) = V(1)*V(2)
    U(3) = ( V(3)/(gamma - one) ) + half*V(1)*V(2)**2

  end subroutine prim2cons_scalar
  subroutine prim2cons_array( U, V )

    real(prec), dimension(:,:), intent(out)  :: U
    real(prec), dimension(:,:), intent(in)   :: V

    U(1,:) = V(1,:)
    U(2,:) = V(1,:)*V(2,:)
    U(3,:) = ( V(3,:)/(gamma - one) ) + half*V(1,:)*V(2,:)**2

  end subroutine prim2cons_array
  elemental subroutine prim2cons_1D( U1, U2, U3, V1, V2, V3 )
    real(prec), intent(out) :: U1, U2, U3
    real(prec), intent(in)  :: V1, V2, V3

    U1 = V1
    U2 = V1*V2
    U3 = ( V3/(gamma - one) ) + half*V1*V2**2
  end subroutine prim2cons_1D

  !================================ update_mach ==============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine update_mach_scalar( V, M )

    real(prec), dimension(:),  intent(in)    :: V
    real(prec),                intent(inout) :: M

    M = abs(V(2))/speed_of_sound(V(3),V(1))

  end subroutine update_mach_scalar
  subroutine update_mach_array( V, M )

    real(prec), dimension(:,:),  intent(in) :: V
    real(prec), dimension(:), intent(inout) :: M

    M = abs(V(2,:))/speed_of_sound(V(3,:),V(1,:))

  end subroutine update_mach_array


  elemental function riemann_1(rho, pres)
    real(prec) :: riemann_1
    real(prec), intent(in) :: rho, pres
    riemann_1 = pres/rho**gamma
  end function riemann_1

  elemental function riemann_2(rho, uvel, avel)
    real(prec) :: riemann_2
    real(prec), intent(in) :: rho, uvel, avel
    riemann_2 = uvel + two*avel/(gamma - one)
  end function riemann_2

  elemental function riemann_3(rho, uvel, avel)
    real(prec) :: riemann_3
    real(prec), intent(in) :: rho, uvel, avel
    riemann_3 = uvel - two*avel/(gamma - one)
  end function riemann_3

  elemental subroutine riem2prim(R_plus, R_minus, R_zero, rho, uvel, avel, pres)
    real(prec), intent(in) :: R_plus, R_minus, R_zero
    real(prec), intent(out) :: rho, uvel, avel, pres

    uvel = half*(R_plus + R_minus)
    avel = fourth*(gamma-one)*(R_plus - R_minus)
    rho = ( avel**2/(gamma*R_zero) )**( one/(gamma - one) )
    pres = rho*avel**2/gamma
  end subroutine riem2prim

  !============================== update_states  =============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine update_states( soln )

    type(soln_t) :: soln
    call cons2prim(soln%U,soln%V)
    call limit_primitives(soln%V)
    soln%asnd = speed_of_sound(soln%V(3,:),soln%V(1,:))

    soln%mach = abs(soln%V(2,:))/soln%asnd
    soln%temp = T0/(one + half*(gamma - one)*soln%mach**2)

  end subroutine update_states




  !========================= limit_primitives ================================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine limit_primitives(V)

    !real(prec), dimension(:,:), intent(inout) :: U
    real(prec), dimension(:,:), intent(inout) :: V
    integer :: i
    !logical, dimension(lbound(U,1):ubound(U,1)) :: mask

    !mask = .false.

    !where ( (V(:,1)<0.001_prec).or.(V(:,2)<10.0_prec).or.(V(:,3)<500.0_prec) )
    !  mask = .true.
    !end where

    do i = lbound(V,2),ubound(V,2)
      if (V(1,i)<1.0e-12_prec) then
        V(1,i) = 1.0e-12_prec
      end if
      if (V(3,i)<1.0e-12_prec) then
        V(3,i) = 1.0e-12_prec
      end if
      !if (mask(i)) then
      !  call prim2cons(U,V)
      !end if
    end do
    !call prim2cons(U,V)

  end subroutine limit_primitives

  !========================= isentropic_relations ============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine isentropic_relations(M,V)

    real(prec), dimension(:,:), intent(inout) :: V
    real(prec), dimension(:),   intent(in) :: M

    real(prec), dimension(lbound(V,2):ubound(V,2))  :: T
    real(prec), dimension(lbound(V,2):ubound(V,2)) :: a

    T = T0/(one + half*(gamma - one)*M(:)**2)

    V(3,:) = 1000.0_prec*p0/(one + half*(gamma - one)*M**2)**(gamma/(gamma-1))
    V(1,:) = V(3,:)/(R_gas*T)
    a = speed_of_sound(V(3,:),V(1,:))
    V(2,:) = M*a

  end subroutine isentropic_relations

end module variable_conversion
