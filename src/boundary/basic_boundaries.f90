module basic_boundaries

  use set_precision, only : prec
  use set_constants, only : zero, one, two
  use set_inputs,    only : imax, neq, eps, i_low, i_high, ig_low, ig_high
  use variable_conversion, only : prim2cons, cons2prim, riem2prim, &
                                  speed_of_sound, isentropic_relations, &
                                  update_mach, riemann_1, riemann_2, riemann_3
  use soln_type, only : soln_t

  implicit none

  private
  public :: explicit_characteristic_bndry

contains

  subroutine explicit_characteristic_bndry(soln, VL, VR)
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq), intent(in) :: VL, VR
    real(prec) :: rho, uvel, avel, pres
    real(prec) :: rho_ghost, uvel_ghost, avel_ghost, pres_ghost
    real(prec) :: Rp_ghost, Rm_ghost, R0_ghost
    real(prec), dimension(2) :: tmp
    logical :: flag
    integer :: i

    character(9) :: str1, str2

    flag = .false.

    ! Left boundary
    rho  = soln%V(1,i_low)
    uvel = soln%V(2,i_low)
    pres = soln%V(3,i_low)
    avel = soln%asnd(i_low)
    if ( uvel < zero ) then ! inflow
      if ( abs(uvel) >= avel ) then ! supersonic
        call Left_sup_in_bndry(soln,VL)
        str1 = 'SUP_L_IN'
      else ! subsonic
        call Left_sub_in_bndry(soln,VL)
        flag = .true.
        str1 = 'SUB_L_IN'
      end if
    else ! outflow
      if ( abs(uvel) >= avel ) then ! supersonic
        call Left_sup_out_bndry(soln,VL)
        str1 = 'SUP_L_OUT'
      else ! subsonic
        call Left_sub_out_bndry(soln,VL)
        flag = .true.
        str1 = 'SUB_L_OUT'
      end if
    end if

    ! Right boundary
    rho  = soln%V(1,i_high)
    uvel = soln%V(2,i_high)
    pres = soln%V(3,i_high)
    avel = soln%asnd(i_high)
    if ( uvel < zero ) then ! inflow
      if ( abs(uvel) >= avel ) then ! supersonic
        call Right_sup_in_bndry(soln,VR)
        str2 = 'SUP_R_IN'
      else ! subsonic
        call Right_sub_in_bndry(soln,VR,flag)
        str2 = 'SUB_R_IN'
      end if
    else ! outflow
      if ( abs(uvel) >= avel ) then ! supersonic
        call Right_sup_out_bndry(soln,VR)
        str2 = 'SUP_R_OUT'
      else ! subsonic
        call Right_sub_out_bndry(soln,VR,flag)
        str2 = 'SUB_R_OUT'
      end if
    end if

    !write(*,'(F10.6,A3,A9,A2,A9)') soln%time,' : ',str1,', ', str2

  end subroutine explicit_characteristic_bndry

  subroutine Left_sub_in_bndry(soln,Vspec)
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq), intent(in) :: Vspec
    real(prec) :: rho, uvel, avel, pres
    real(prec) :: rho_ghost, uvel_ghost, avel_ghost, pres_ghost
    real(prec) :: Rp_ghost, Rm_ghost, R0_ghost
    real(prec), dimension(2) :: tmp
    integer :: i

    rho  = soln%V(1,i_low)
    uvel = soln%V(2,i_low)
    pres = soln%V(3,i_low)
    avel = soln%asnd(i_low)

    ! left subsonic inflow
    R0_ghost = riemann_1(Vspec(1),Vspec(3))
    Rp_ghost = riemann_2(Vspec(1),uvel,avel)
    tmp(1)   = riemann_3(rho,uvel,avel)
    tmp(2)   = riemann_3(soln%V(1,i_low+1),soln%V(2,i_low+1),soln%asnd(i_low+1))
    do i = i_low-1,ig_low,-1
      Rm_ghost = two*tmp(1) - tmp(2)
      tmp(2) = tmp(1)
      tmp(1) = Rm_ghost

      call riem2prim(Rp_ghost,Rm_ghost,R0_ghost, &
            rho_ghost,uvel_ghost,avel_ghost,pres_ghost)

      soln%V(:,i) = (/ rho_ghost, uvel_ghost, pres_ghost /)
      soln%asnd(i) = avel_ghost
    end do
  end subroutine Left_sub_in_bndry

  subroutine Right_sub_in_bndry(soln,Vspec,flag)
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq), intent(in) :: Vspec
    real(prec) :: rho, uvel, avel, pres
    real(prec) :: rho_ghost, uvel_ghost, avel_ghost, pres_ghost
    real(prec) :: Rp_ghost, Rm_ghost, R0_ghost
    real(prec), dimension(2) :: tmp
    logical :: flag
    integer :: i

    rho  = soln%V(1,i_high)
    uvel = soln%V(2,i_high)
    pres = soln%V(3,i_high)
    avel = soln%asnd(i_high)

    ! right subsonic inflow
    if (flag.eqv..false.) then
      R0_ghost = riemann_1(Vspec(1),Vspec(3))
      Rp_ghost = riemann_2(Vspec(1),uvel,avel)
    else
      R0_ghost = riemann_1(Vspec(1),pres)
      Rp_ghost = riemann_2(Vspec(1),Vspec(2),avel)
    end if
    tmp(1)   = riemann_3(rho,uvel,avel)
    tmp(2)   = riemann_3(soln%V(1,i_high-1),&
                         soln%V(2,i_high-1),soln%asnd(i_high-1))
    do i = i_high+1,ig_high
      Rm_ghost = two*tmp(1) - tmp(2)
      tmp(2) = tmp(1)
      tmp(1) = Rm_ghost

      call riem2prim(Rp_ghost,Rm_ghost,R0_ghost, &
            rho_ghost,uvel_ghost,avel_ghost,pres_ghost)

      soln%V(:,i) = (/ rho_ghost, uvel_ghost, pres_ghost /)
      soln%asnd(i) = avel_ghost
    end do
  end subroutine Right_sub_in_bndry

  subroutine Left_sub_out_bndry(soln,Vspec)
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq), intent(in) :: Vspec
    real(prec) :: rho, uvel, avel, pres, aspec
    real(prec) :: rho_ghost, uvel_ghost, avel_ghost, pres_ghost
    real(prec) :: Rp_ghost, Rm_ghost, R0_ghost
    real(prec), dimension(2) :: tmp1, tmp2
    logical :: flag
    integer :: i

    rho  = soln%V(1,i_low)
    uvel = soln%V(2,i_low)
    pres = soln%V(3,i_low)
    avel = soln%asnd(i_low)
    aspec = speed_of_sound(Vspec(3),rho)

    ! left subsonic outflow
    Rm_ghost = riemann_3(rho,uvel,aspec)

    tmp1(1)  = riemann_1(rho,pres)
    tmp1(2)  = riemann_1(soln%V(1,i_low+1),soln%V(3,i_low+1))
    tmp2(1)  = riemann_2(rho,uvel,avel)
    tmp2(2)  = riemann_2(soln%V(1,i_low+1),soln%V(2,i_low+1),soln%asnd(i_low+1))

    do i = i_low-1,ig_low,-1
      R0_ghost = two*tmp1(1) - tmp1(2)
      Rp_ghost = two*tmp2(1) - tmp2(2)
      tmp1(2) = tmp1(1)
      tmp1(1) = R0_ghost
      tmp2(2) = tmp2(1)
      tmp2(1) = Rp_ghost

      call riem2prim(Rp_ghost,Rm_ghost,R0_ghost, &
            rho_ghost,uvel_ghost,avel_ghost,pres_ghost)

      soln%V(:,i) = (/ rho_ghost, uvel_ghost, pres_ghost /)
      soln%asnd(i) = avel_ghost
    end do
  end subroutine Left_sub_out_bndry

  subroutine Right_sub_out_bndry(soln,Vspec,flag)
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq), intent(in) :: Vspec
    real(prec) :: rho, uvel, avel, pres, aspec
    real(prec) :: rho_ghost, uvel_ghost, avel_ghost, pres_ghost
    real(prec) :: Rp_ghost, Rm_ghost, R0_ghost
    real(prec), dimension(2) :: tmp1, tmp2
    logical :: flag
    integer :: i

    rho  = soln%V(1,i_high)
    uvel = soln%V(2,i_high)
    pres = soln%V(3,i_high)
    avel = soln%asnd(i_high)
    aspec = speed_of_sound(Vspec(3),rho)

    ! right subsonic outflow
    if (flag.eqv..false.) then
      Rm_ghost = riemann_3(rho,uvel,aspec)
    else
      Rm_ghost = riemann_3(rho,Vspec(2),avel)
    end if

    tmp1(1)  = riemann_1(rho,pres)
    tmp1(2)  = riemann_1(soln%V(1,i_high-1),&
                         soln%V(3,i_high-1))
    tmp2(1)  = riemann_2(rho,uvel,avel)
    tmp2(2)  = riemann_2(soln%V(1,i_high-1),&
                         soln%V(2,i_high-1),soln%asnd(i_high-1))

    do i = i_high+1,ig_high
      R0_ghost = two*tmp1(1) - tmp1(2)
      Rp_ghost = two*tmp2(1) - tmp2(2)
      tmp1(2) = tmp1(1)
      tmp1(1) = R0_ghost
      tmp2(2) = tmp2(1)
      tmp2(1) = Rp_ghost

      call riem2prim(Rp_ghost,Rm_ghost,R0_ghost, &
            rho_ghost,uvel_ghost,avel_ghost,pres_ghost)

      soln%V(:,i) = (/ rho_ghost, uvel_ghost, pres_ghost /)
      soln%asnd(i) = avel_ghost
    end do
  end subroutine Right_sub_out_bndry

  subroutine Left_sup_in_bndry(soln,Vspec)
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq), intent(in) :: Vspec
    real(prec) :: rho, uvel, avel, pres
    real(prec) :: rho_ghost, uvel_ghost, avel_ghost, pres_ghost
    real(prec) :: Rp_ghost, Rm_ghost, R0_ghost
    integer :: i

    rho  = Vspec(1)
    uvel = Vspec(2)
    pres = Vspec(3)
    avel = speed_of_sound(pres,rho)

    ! left supersonic inflow
    R0_ghost = riemann_1(rho,pres)
    Rp_ghost = riemann_2(rho,uvel,avel)
    Rm_ghost = riemann_3(rho,uvel,avel)
    call riem2prim(Rp_ghost,Rm_ghost,R0_ghost, &
          rho_ghost,uvel_ghost,avel_ghost,pres_ghost)

    do i = i_low-1,ig_low,-1
      soln%V(:,i) = (/ rho_ghost, uvel_ghost, pres_ghost /)
      soln%asnd(i) = avel_ghost
    end do
  end subroutine Left_sup_in_bndry

  subroutine Right_sup_in_bndry(soln,Vspec)
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq), intent(in) :: Vspec
    real(prec) :: rho, uvel, avel, pres
    real(prec) :: rho_ghost, uvel_ghost, avel_ghost, pres_ghost
    real(prec) :: Rp_ghost, Rm_ghost, R0_ghost
    integer :: i

    rho  = Vspec(1)
    uvel = Vspec(2)
    pres = Vspec(3)
    avel = speed_of_sound(pres,rho)

    ! left supersonic inflow
    R0_ghost = riemann_1(rho,pres)
    Rp_ghost = riemann_2(rho,uvel,avel)
    Rm_ghost = riemann_3(rho,uvel,avel)
    call riem2prim(Rp_ghost,Rm_ghost,R0_ghost, &
          rho_ghost,uvel_ghost,avel_ghost,pres_ghost)

    do i = i_high+1,ig_high
      soln%V(:,i) = (/ rho_ghost, uvel_ghost, pres_ghost /)
      soln%asnd(i) = avel_ghost
    end do
  end subroutine Right_sup_in_bndry

  subroutine Left_sup_out_bndry(soln,Vspec)
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq), intent(in) :: Vspec
    real(prec) :: rho, uvel, avel, pres
    real(prec) :: rho_ghost, uvel_ghost, avel_ghost, pres_ghost
    real(prec) :: Rp_ghost, Rm_ghost, R0_ghost
    real(prec), dimension(2) :: tmp1, tmp2, tmp3
    integer :: i

    rho  = Vspec(1)
    uvel = Vspec(2)
    pres = Vspec(3)
    avel = speed_of_sound(pres,rho)

    ! left supersonic outflow
    tmp1(1)  = riemann_1(soln%V(1,i_low),&
                         soln%V(3,i_low))
    tmp1(2)  = riemann_1(soln%V(1,i_low+1),&
                         soln%V(3,i_low+1))
    tmp2(1)  = riemann_2(soln%V(1,i_low),&
                         soln%V(2,i_low),soln%asnd(i_low))
    tmp2(2)  = riemann_2(soln%V(1,i_low+1),&
                         soln%V(2,i_low+1),soln%asnd(i_low+1))
    tmp3(1)  = riemann_3(soln%V(1,i_low),&
                         soln%V(2,i_low),soln%asnd(i_low))
    tmp3(2)  = riemann_3(soln%V(1,i_low+1),&
                         soln%V(2,i_low+1),soln%asnd(i_low+1))

    do i = i_low-1,ig_low,-1
      R0_ghost = two*tmp1(1) - tmp1(2)
      Rp_ghost = two*tmp2(1) - tmp2(2)
      Rm_ghost = two*tmp3(1) - tmp3(2)
      tmp1(2) = tmp1(1)
      tmp1(1) = R0_ghost
      tmp2(2) = tmp2(1)
      tmp2(1) = Rp_ghost
      tmp3(2) = tmp3(1)
      tmp3(1) = Rm_ghost

      call riem2prim(Rp_ghost,Rm_ghost,R0_ghost, &
            rho_ghost,uvel_ghost,avel_ghost,pres_ghost)

      soln%V(:,i) = (/ rho_ghost, uvel_ghost, pres_ghost /)
      soln%asnd(i) = avel_ghost
    end do
  end subroutine Left_sup_out_bndry

  subroutine Right_sup_out_bndry(soln,Vspec)
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq), intent(in) :: Vspec
    real(prec) :: rho, uvel, avel, pres
    real(prec) :: rho_ghost, uvel_ghost, avel_ghost, pres_ghost
    real(prec) :: Rp_ghost, Rm_ghost, R0_ghost
    real(prec), dimension(2) :: tmp1, tmp2, tmp3
    integer :: i

    rho  = Vspec(1)
    uvel = Vspec(2)
    pres = Vspec(3)
    avel = speed_of_sound(pres,rho)

    ! left supersonic inflow
    tmp1(1)  = riemann_1(soln%V(1,i_high),&
                         soln%V(3,i_high))
    tmp1(2)  = riemann_1(soln%V(1,i_high-1),&
                         soln%V(3,i_high-1))
    tmp2(1)  = riemann_2(soln%V(1,i_high),&
                         soln%V(2,i_high),soln%asnd(i_high))
    tmp2(2)  = riemann_2(soln%V(1,i_high-1),&
                         soln%V(2,i_high-1),soln%asnd(i_high-1))
    tmp3(1)  = riemann_3(soln%V(1,i_high),&
                         soln%V(2,i_high),soln%asnd(i_high))
    tmp3(2)  = riemann_3(soln%V(1,i_high-1),&
                         soln%V(2,i_high-1),soln%asnd(i_high-1))

    do i = i_high+1,ig_high
      R0_ghost = two*tmp1(1) - tmp1(2)
      Rp_ghost = two*tmp2(1) - tmp2(2)
      Rm_ghost = two*tmp3(1) - tmp3(2)
      tmp1(2) = tmp1(1)
      tmp1(1) = R0_ghost
      tmp2(2) = tmp2(1)
      tmp2(1) = Rp_ghost
      tmp3(2) = tmp3(1)
      tmp3(1) = Rm_ghost

      call riem2prim(Rp_ghost,Rm_ghost,R0_ghost, &
            rho_ghost,uvel_ghost,avel_ghost,pres_ghost)

      soln%V(:,i) = (/ rho_ghost, uvel_ghost, pres_ghost /)
      soln%asnd(i) = avel_ghost
    end do
  end subroutine Right_sup_out_bndry

end module basic_boundaries
