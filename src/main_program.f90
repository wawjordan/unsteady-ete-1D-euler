program main_program

  use set_constants, only : set_derived_constants, zero, one
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs, xmin, xmax, neq
  use set_inputs, only : max_iter, tol, soln_save, res_save
  use set_inputs, only : leftV, rightV, leftU, rightU, flux_scheme
  use set_inputs, only : limiter_scheme, limiter_freeze, res_out, cons
  use variable_conversion
  use time_integration, only : explicit_Euler, explicit_RK, calc_time_step, residual_norms
  use basic_boundaries, only : explicit_characteristic_bndry
  use limiter_type, only : select_limiter
  use flux_calc, only : select_flux, flux_fun, calc_flux_1D
  use other_subroutines
  use geometry, only : setup_geometry, teardown_geometry
  use init_problem, only : initialize
  use namelist, only : read_namelist
  use grid_type
  use soln_type
  use exact_soln_type, only : exact_soln_t, setup_exact_soln, &
           destroy_exact_soln, sample_exact_soln

 ! use dns_matrix_derived_type, only : dns_matrix_t, allocate_dns_matrix

 ! use iterative_solvers, only : gmres
  implicit none

  external dgesv
  character(len=100) :: header_str1
  character(len=100) :: niter
  integer :: i, j, pnorm
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  type( exact_soln_t ) :: ex_soln
  real(prec), dimension(3) :: VL, VR
  real(prec), dimension(:), allocatable :: xi
  real(prec), dimension(3,3) :: A
  real(prec), dimension(3) :: B
  integer :: IPIV(3),INFO

  !A = transpose( reshape( (/ 8.0_prec, 1.0_prec, 6.0_prec, &
  !                           3.0_prec, 5.0_prec, 7.0_prec, &
  !                           4.0_prec, 9.0_prec, 2.0_prec /), shape(A) ) )
  !B = (/1.0_prec, 3.0_prec, 4.0_prec /)

  !call dgesv(neq,neq,A,neq,IPIV,B,neq,INFO)

  !write(*,*) 'A = '
  !do j = 1,neq
  !  do i = 1,neq
  !    write(*,'(F3.1,A1)') A(i,j)
  !  end do
  !  write(*,*)
  !end do
  !write(*,*) 'B = '
  !do j = 1,neq
    !do i = 1,neq
    !  write(*,'(F3.1,A1)') B(i,j)
    !end do
  !  write(*,'(F10.6,A1)') B(j)
  !  write(*,*)
  !end do

  !stop

  pnorm = 1
  j = 0

  VL = (/one,zero,one/)
  VR = (/0.125_prec,zero,0.1_prec/)
  !VL = (/one,-two,0.4_prec/)
  !VR = (/one,two,0.4_prec/)

  !call allocate_dns_matrix(A_mat,3)

  !call A_mat%add_array(1,1,3,3,array1)

  !call gmres(A_mat,b,x)

  !call A_mat%extract_array((/1,1/),(/3,3/),array2)

  !do i = 1,3
  !  write(*,*) array2(i,1)
  !end do
  !write(*,*)
  !do i = 1,3
  !  write(*,*) b(i)
  !end do
  !write(*,*)
  !do i = 1,3
  !  write(*,*) x(i)
  !end do


  !call A_mat%destroy

  !stop



  open(50,file='temp.txt',status='unknown')

  100 format(1(I0.8),3(G20.7))

  write(header_str1,*) " Iter  |   ||density||    |"// &
  & "   ||velocity||    |   ||pressure||    |"
  call set_derived_constants
  call set_fluid_constants
  call read_namelist('input.nml')
  call output_file_headers
  call set_derived_inputs
  allocate(xi(imax))

  xi = [ (xmin + real(i,prec)/real(imax,prec)*(xmax-xmin), &
          i=i_low-1,i_high) ]
  call setup_geometry(grid,soln)
  deallocate(xi)
  call select_flux()
  !call select_limiter(limiter_scheme)

  call setup_exact_soln( ex_soln, grid, VL, VR)
  call initialize(grid,soln,VL,VR)
  call output_soln(grid,soln,ex_soln,0)


  !call calculate_sources(soln%V(3,:),grid%dAc,soln%src)
!
!  call MUSCL_extrap( soln%V, leftV, rightV )
!  call prim2cons(leftU,leftV)
!  call prim2cons(rightU,rightV)
!  call flux_fun(leftU,rightU,soln%F)

  call calc_time_step(grid%dx,soln%V,soln%asnd,soln%lambda,soln%dt)
  !call calc_flux_1D(grid,soln)
  !call explicit_euler(grid,soln%S,soln%dt,soln%F,soln%U,soln%R)
  !call explicit_characteristic_bndry(soln, VL, VR)
  call explicit_RK(grid,soln,VL,VR)
  soln%time = soln%time + minval(soln%dt)
  call sample_exact_soln(ex_soln,grid,soln%time)
  call output_soln(grid,soln,ex_soln,1)
  call residual_norms(soln%R,soln%rinit,(/one,one,one/))
  soln%rold = soln%rinit
  write(*,*) 'Residual Norms: Iteration 0'
  write(*,*) header_str1
  write(*,100) j, soln%rinit(1), soln%rinit(2), soln%rinit(3)
  write(*,*)
  write(*,*) 'Relative Residual Norms:'
!
!
  do j = 2,max_iter

    !call calculate_sources(soln,grid)
    call calc_time_step(grid%dx,soln%V,soln%asnd,soln%lambda,soln%dt)
    call explicit_RK(grid,soln,VL,VR)

    !call explicit_characteristic_bndry(soln, VL, VR)
    !call update_states( soln )

    !call MUSCL_extrap( soln%V, leftV, rightV )
    !call prim2cons(leftU,leftV)
    !call prim2cons(rightU,rightV)
    !call flux_fun(leftU,rightU,soln%F)
    !call calc_flux_1D(grid,soln)

    !call explicit_euler(grid,soln%S,soln%dt,soln%F,soln%U,soln%R)
    !call update_states( soln )

    soln%time = soln%time + minval(soln%dt)
    call sample_exact_soln(ex_soln,grid,soln%time)
    !VL = ex_soln%Vq(:,i_low)
    !VR = ex_soln%Vq(:,i_high)

   if (mod(j,soln_save)==0) then
      call calc_de( soln, ex_soln, soln%DE, soln%DEnorm, pnorm, cons )
      call output_soln(grid,soln,ex_soln,j)
    end if
    call residual_norms(soln%R,soln%rnorm,soln%rinit)

    !if (all(soln%rnorm<tol) ) then
    !  exit
    !elseif (any(soln%rnorm>soln%rold)) then
    !    limiter_freeze = .true.
    !end if
    soln%rold = soln%rnorm
    if (mod(j,10*res_out)==0) then
      write(*,*) header_str1
    end if
    if (mod(j,res_out)==0) then
      write(*,100) j, soln%rnorm(1), soln%rnorm(2), soln%rnorm(3)
    end if
    if (mod(j,res_save)==0) then
      call output_res(soln,j)
    end if
  end do

  if (mod(j,res_out)/=0) then
    write(*,100) j, soln%rnorm(1), soln%rnorm(2), soln%rnorm(3)
  end if
  if (all(soln%rnorm<tol) ) then
    write(*,*) 'Solution converged!'
  else
    write(niter,'(I7)') max_iter
    write(*,*) 'Solution failed to converge in ',&
       & trim(adjustl(niter)),' iterations'
  end if
  write(50,*) imax, j, soln%rnorm(1), soln%rnorm(2), soln%rnorm(3)

  if (mod(j,res_out)/=0) then
    call output_soln(grid,soln,ex_soln,j+1)
    call output_res(soln,j)
  end if
  call destroy_exact_soln( ex_soln )
  call teardown_geometry(grid,soln)
  close(50)
end program main_program
