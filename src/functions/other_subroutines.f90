module other_subroutines

  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, half, fourth
  use set_inputs, only : imax, neq, i_low, i_high, ig_low, ig_high
  use set_inputs, only : epsM, kappaM
  use fluid_constants, only : gamma
  use variable_conversion
  use flux_jacobians, only : VL_flux_jac
  !use limiter_calc, only : limiter_fun!, calc_consecutive_variations
  use soln_type, only : soln_t
  use exact_soln_type, only : exact_soln_t
  use grid_type, only : grid_t

  implicit none

  contains

  !============================= calculate_sources ===========================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine calculate_sources(soln,grid)

    type(soln_t), intent(inout) :: soln
    type(grid_t), intent(in)    :: grid
    real(prec), dimension(i_low:i_high) :: pres, darea

    pres  = soln%V(3,i_low:i_high)
    darea = grid%dAc(i_low:i_high)

    soln%S(2,:) = source_term_1D_nozzle(pres,darea)

  end subroutine calculate_sources

  elemental function source_term_1D_nozzle(pres,darea)
    real(prec) :: source_term_1D_nozzle
    real(prec), intent(in)  :: pres, darea

    source_term_1D_nozzle = pres*darea
  end function source_term_1D_nozzle

  subroutine build_LHS_matrix(soln,grid)
    type(soln_t), intent(inout) :: soln
    type(grid_t), intent(in)    :: grid
    real(prec), dimension(neq,i_low-1:i_high) :: Left, Right
    real(prec), dimension(neq,neq) :: dfduim1, dfdui, dfduip1, dfduip2
    real(prec), dimension(neq,neq) :: dgduim1, dgdui, dgduip1, dgduip2
    real(prec) :: vt
    integer :: i, j

    call MUSCL_extrap(soln,Left,Right)
    ! left boundary -> calculate dfdu_{i-1/2} wrt i, i+1, & i+2
!    i = i_low-1
!    call VL_flux_jac( Left(:,i), Right(:,i), &
!            soln%U(:,i-1), soln%U(:,i), soln%U(:,i+1), soln%U(:,i+2), &
!            soln%psi_p(:,i-1), soln%psi_p(:,i), &
!            soln%psi_m(:,i), soln%psi_m(:,i+1), &
!            soln%dpsi_p(:,i-1), soln%dpsi_p(:,i), &
!            soln%dpsi_m(:,i), soln%dpsi_m(:,i+1), &
!            dfduim1, dfdui, dfduip1, dfduip2 )
!    soln%LHS(:,:,i+1,3) = -grid%Ai(i)*dfduim1    !A*df{i-1/2}/du{1}
!    soln%LHS(:,:,i+1,4) = -grid%Ai(i)*dfdui  !A*df{i-1/2}/du{2}
!    soln%LHS(:,:,i+1,5) = -grid%Ai(i)*dfduip1  !A*df{i-1/2}/du{3}
    do i = i_low,i_high-1
      !write(*,*) soln%psi_p(1,i), soln%psi_m(1,i)
      call VL_flux_jac( Left(:,i), Right(:,i), &
              soln%U(:,i-1), soln%U(:,i), soln%U(:,i+1), soln%U(:,i+2), &
              soln%psi_p(:,i-1), soln%psi_p(:,i), &
              soln%psi_m(:,i), soln%psi_m(:,i+1), &
              soln%dpsi_p(:,i-1), soln%dpsi_p(:,i), &
              soln%dpsi_m(:,i), soln%dpsi_m(:,i+1), &
              dfduim1, dfdui, dfduip1, dfduip2 )
      do j = 1,neq
        write(*,*) dfdui(j,:)
      end do
      write(*,*)
      vt = grid%Ac(i)*grid%dx(i)/soln%dt(i)
      soln%LHS(:,:,i,1) = soln%LHS(:,:,i,1) + grid%Ai(i)*dfduim1
      soln%LHS(:,:,i,2) = soln%LHS(:,:,i,2) + grid%Ai(i)*dfdui
      soln%LHS(:,:,i,3) = vt + soln%LHS(:,:,i,3) + grid%Ai(i)*dfduip1
      soln%LHS(:,:,i,4) = soln%LHS(:,:,i,4) + grid%Ai(i)*dfduip2

      soln%LHS(:,:,i+1,2) = -grid%Ai(i)*dfduim1
      soln%LHS(:,:,i+1,3) = -grid%Ai(i)*dfdui
      soln%LHS(:,:,i+1,4) = -grid%Ai(i)*dfduip1
      soln%LHS(:,:,i+1,5) = -grid%Ai(i)*dfduip2
    end do
!    i = i_high+1
!    call VL_flux_jac( Left(:,i), Right(:,i), &
!            soln%U(:,i-1), soln%U(:,i), soln%U(:,i+1), soln%U(:,i+2), &
!            soln%psi_p(:,i-1), soln%psi_p(:,i), &
!            soln%psi_m(:,i), soln%psi_m(:,i+1), &
!            soln%dpsi_p(:,i-1), soln%dpsi_p(:,i), &
!            soln%dpsi_m(:,i), soln%dpsi_m(:,i+1), &
!            dfduim1, dfdui, dfduip1, dfduip2 )
!      soln%LHS(:,:,i,1) = soln%LHS(:,:,i,1) + grid%Ai(i)*dfdui
!      soln%LHS(:,:,i,2) = soln%LHS(:,:,i,2) + grid%Ai(i)*dfduip1
!      soln%LHS(:,:,i,3) = soln%LHS(:,:,i,3) + grid%Ai(i)*dfduip2

  end subroutine build_LHS_matrix


  !================================ MUSCL_extrap =============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine MUSCL_extrap(soln,Left,Right)
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq,i_low-1:i_high), intent(out) :: Left, Right
    integer :: i
    do i = i_low-1,i_high
      Left(:,i) = soln%V(:,i) + fourth*epsM*( &
            (one-kappaM)*soln%psi_p(:,i-1)*( soln%V(:,i) - soln%V(:,i-1) ) + &
            (one+kappaM)*soln%psi_m(:,i)  *( soln%V(:,i+1) - soln%V(:,i) ) )

      Right(:,i) = soln%V(:,i+1) - fourth*epsM*( &
            (one-kappaM)*soln%psi_m(:,i+1)*( soln%V(:,i+2) - soln%V(:,i+1) ) + &
            (one+kappaM)*soln%psi_p(:,i)  *( soln%V(:,i+1) - soln%V(:,i) ) )
    end do

   ! Left  = soln%V(:,i_low-1:i_high) + fourth*epsM*(       &
   !         (one-kappaM)*soln%psi_p(:,i_low-2:i_high-1)* &
   !                           ( soln%V(:,i_low-1:i_high)  &
   !                           - soln%V(:,i_low-2:i_high-1) )&
   !       + (one+kappaM)*soln%psi_m(:,i_low-1:i_high)* &
   !                           ( soln%V(:,i_low  :i_high+1)  &
   !                           - soln%V(:,i_low-1:i_high) ))
   ! Right  = soln%V(:,i_low  :i_high+1) - fourth*epsM*(       &
   !         (one-kappaM)*soln%psi_m(:,i_low  :i_high+1)* &
   !                           ( soln%V(:,i_low+1:i_high+2)  &
   !                           - soln%V(:,i_low  :i_high+1) )&
   !       + (one+kappaM)*soln%psi_p(:,i_low-1:i_high)* &
   !                           ( soln%V(:,i_low  :i_high+1)  &
   !                           - soln%V(:,i_low-1:i_high) ))

  end subroutine MUSCL_extrap

  !================================== calc_de ==== ===========================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine calc_de( soln, exact_soln, DE, DEnorm, pnorm, cons )

    type(soln_t), intent(inout) :: soln
    type(exact_soln_t), intent(in) :: exact_soln
    logical, intent(in) :: cons
    real(prec), dimension(neq,i_low:i_high), intent(out) :: DE
    real(prec), dimension(neq), intent(out) :: DEnorm
    integer, intent(in) :: pnorm
    real(prec) :: Linv
    Linv = one/real(i_high-i_low)
    if (cons) then
      DE = soln%U(:,i_low:i_high) - exact_soln%Uq(:,i_low:i_high)
    else
      DE = soln%V(:,i_low:i_high) - exact_soln%Vq(:,i_low:i_high)
    end if

    if (pnorm == 0) then
      DEnorm(1) = maxval(abs(DE(1,:)))
      DEnorm(2) = maxval(abs(DE(2,:)))
      DEnorm(3) = maxval(abs(DE(3,:)))
    elseif (pnorm == 1) then
      DEnorm(1) = Linv*sum(abs(DE(1,:)))
      DEnorm(2) = Linv*sum(abs(DE(2,:)))
      DEnorm(3) = Linv*sum(abs(DE(3,:)))
    elseif (pnorm == 2) then
      DEnorm(1) = sqrt(Linv*sum(DE(1,:)**2))
      DEnorm(2) = sqrt(Linv*sum(DE(2,:)**2))
      DEnorm(3) = sqrt(Linv*sum(DE(3,:)**2))
    else
      DEnorm(1) = maxval(abs(DE(1,:)))
      DEnorm(2) = maxval(abs(DE(2,:)))
      DEnorm(3) = maxval(abs(DE(3,:)))
    end if

  end subroutine calc_de

  !========================== output_file_headers ============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine output_file_headers

    use set_inputs, only : imax
    use set_inputs, only : flux_scheme, limiter_scheme, beta_lim
    use set_inputs, only : CFL, epsM

    character(len=1024) :: dirname_hist
    character(len=1024) :: dirname_field
    character(len=1024) :: filename
    character(len=64) ::   ncells_str
    character(len=64) ::   CFL_str
    character(len=64) ::   flux_str
    character(len=64) ::   order_str
    character(len=64) ::   limiter_str
    character(len=64) ::   kappa_str

    ! Set up output directories
    write (ncells_str  , "(A1,I0.3,A1)") "N"  , imax   , "_"
    write (CFL_str  , "(A4,I0.3)") "CFL-", int(1000*cfl)

    select case(flux_scheme)
    case(1)
      write (flux_str,"(A3)") "V-L"
    case(2)
      write (flux_str,"(A3)") "ROE"
    case default
    end select
    if (nint(epsM).eq.0) then
      write (order_str,"(A1)") "1"
    else
      write (order_str,"(A1)") "2"
    end if

    select case(limiter_scheme)
    case(0)
      write (limiter_str,"(A3)") "_NA"
    case(1)
      write (limiter_str,"(A3)") "_VL"
    case(2)
      write (limiter_str,"(A3)") "_VA"
    case(3)
      write (limiter_str,"(A3)") "_MM"
    case(4)
      write (limiter_str,"(A2,I0.2,A1)") "_B",int(10*beta_lim),"_"
    case default
    end select

    write (kappa_str, "(A2,I3.2,A1)") "_K"  , int(10*kappaM),"_"

    !write(*,*) ncells_str
    !write(*,*) CFL_str
    !write(*,*) flux_str
    !write(*,*) order_str
    !write(*,*) limiter_str
    !write(*,*) kappa_str
    !stop



    write (dirname_hist, *) "/hist/"
    write (dirname_field, *) "/field/"
    write (filename,*)  trim(ncells_str)//  &
    &                   trim(flux_str)//    &
    &                   trim(order_str)//   &
    &                   trim(limiter_str)// &
    &                   trim(kappa_str)//   &
    &                   trim(CFL_str)
    !write(*,*) filename
    !stop

    call execute_command_line ('mkdir -p ../results/hist/')

    call execute_command_line ('mkdir -p ../results/field/')
    ! Set up output files (history and solution)
    open(30,file= '../results/hist/'//trim(adjustl(filename))// &
    &                        '_history.dat',status='unknown')

    write(30,*) 'TITLE = "Q1D res hist ('// &
    & trim(adjustl(filename))//')"'

    write(30,*) 'variables="Iteration""R<sub>1</sub>""R<sub>2</sub>"&
    & "R<sub>3</sub>""DE<sub>1</sub>""DE<sub>2</sub>""DE<sub>3</sub>"'

    open(40,file= '../results/field/'//trim(adjustl(filename))// &
    &                        '_field.dat',status='unknown')

    write(40,*) 'TITLE = "Q1D soln ('// &
    & trim(adjustl(filename))//')"'

    write(40,*) 'variables="x (m)""A (m<sup>2</sup>)"&
    & "<greek>r</greek> (kg/m<sup>3</sup>)""u (m/s)"&
    & "p (N/m<sup>2</sup>)"  &
    & "M""U<sub>1</sub>""U<sub>2</sub>""U<sub>3</sub>"&
    & "M<sub>exact</sub>""<greek>r</greek><sub>exact</sub>"&
    & "u<sub>exact</sub>""p<sub>exact</sub>"&
    & "DE<sub>1</sub>""DE<sub>2</sub>""DE<sub>3</sub>"'

  end subroutine output_file_headers

  !============================= output_soln =================================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine output_soln(grid,soln,ex_soln,num_iter)

    use set_inputs, only : counter

    type( grid_t ), intent(in) :: grid
    type( soln_t ), intent(in) :: soln
    type( exact_soln_t ), intent(in) :: ex_soln
    integer,        intent(in) :: num_iter

    integer :: i, j, low, high

    low = i_low
    high = i_high

    open(40,status='unknown')
    write(40,*) 'zone T="',num_iter,'" '
    write(40,*) 'I=',high-low+2
    write(40,*) 'J=',1
    write(40,*) 'DT=(DOUBLE DOUBLE DOUBLE &
                   & DOUBLE DOUBLE DOUBLE &
                   & DOUBLE DOUBLE DOUBLE &
                   & DOUBLE DOUBLE DOUBLE &
                   & DOUBLE DOUBLE DOUBLE &
                   & DOUBLE)'
    write(40,*) 'DATAPACKING=BLOCK'
    write(40,*) 'VARLOCATION=([1-16]=CELLCENTERED)'
    write(40,*) 'STRANDID=',1
    write(40,*) 'SOLUTIONTIME=',real(num_iter,prec)
    !write(40,*) ( ( grid%xi(i) ,  i = low-1,high), NEW_LINE('a'), j = 1,2 )
    !write(40,*) ( ( real(j-1,prec)*grid%Ai(i), i = low-1,high), &
    !            NEW_LINE('a'), j = 1,2 )
    write(40,*) ( grid%xc(i),  i = low,high )
    write(40,*) ( grid%Ac(i),  i = low,high )
    write(40,*) ( soln%V(1,i),  i = low,high )
    write(40,*) ( soln%V(2,i),  i = low,high )
    write(40,*) ( soln%V(3,i),  i = low,high )
    write(40,*) ( soln%mach(i), i = low,high )
    write(40,*) ( soln%U(1,i),  i = low,high )
    write(40,*) ( soln%U(2,i),  i = low,high )
    write(40,*) ( soln%U(3,i),  i = low,high )
    write(40,*) ( ex_soln%Mq(i)   , i = low,high )
    write(40,*) ( ex_soln%Vq(1,i) , i = low,high )
    write(40,*) ( ex_soln%Vq(2,i) , i = low,high )
    write(40,*) ( ex_soln%Vq(3,i) , i = low,high )
    write(40,*) ( soln%DE(1,i), i = low,high )
    write(40,*) ( soln%DE(2,i), i = low,high )
    write(40,*) ( soln%DE(3,i), i = low,high )
    !do i = i_low,i_high
    !  write(40,*) grid%xc(i),grid%Ac(i),soln%V(1,i),soln%V(2,i),soln%V(3,i),&
    !         & soln%mach(i),soln%U(1,i),soln%U(2,i),soln%U(3,i), &
    !         & ex_soln%Mq(i), ex_soln%Vq(1,i), ex_soln%Vq(2,i),  &
    !         & ex_soln%Vq(3,i), soln%DE(1,i), soln%DE(2,i), soln%DE(3,i)
    !end do

  end subroutine output_soln

  !============================= output_res ==================================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine output_res(soln,num_iter)

    type( soln_t), intent(in) :: soln
    integer, intent(in) :: num_iter
    integer :: i
    write(30,*) num_iter,(soln%rnorm(i),i=1,neq),(soln%DEnorm(i),i=1,neq)

  end subroutine output_res

end module other_subroutines
