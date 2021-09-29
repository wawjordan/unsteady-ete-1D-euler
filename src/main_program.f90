program main_program

  use set_constants, only : set_derived_constants, zero, one
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs, xmin, xmax, neq
  use set_inputs, only : max_iter, tol, soln_save, res_save
  use set_inputs, only : leftV, rightV, leftU, rightU, flux_scheme
  use set_inputs, only : limiter_scheme, limiter_freeze, res_out, cons
  use variable_conversion
  !use time_integration, only : explicit_Euler, explicit_RK, calc_time_step, &
  !                             residual_norms, implicit_Euler
  !use basic_boundaries, only : explicit_characteristic_bndry
  use limiter_type, only : select_limiter
  use flux_calc, only : select_flux, flux_fun, calc_flux_1D
  use flux_jacobians, only : calc_vl_dfdu, flux_jac_cons1D
  use other_subroutines
  use geometry, only : setup_geometry, teardown_geometry
  use init_problem, only : initialize
  use namelist, only : read_namelist
  use grid_type
  use soln_type
  use exact_soln_type, only : exact_soln_t, setup_exact_soln, &
           destroy_exact_soln, sample_exact_soln

  implicit none

  character(len=100) :: header_str1
  character(len=100) :: niter
  integer :: i, j, pnorm
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  type( exact_soln_t ) :: ex_soln
  real(prec), dimension(3) :: VL, VR, f
  real(prec), dimension(:), allocatable :: xi
  real(prec), dimension(3,3) :: A, B, C
  !real(prec), dimension(3) :: B
  integer :: IPIV(3),INFO

  pnorm = 1
  j = 0

  VL = (/one,zero,one/)
  VR = (/0.125_prec,zero,0.1_prec/)

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
  call setup_exact_soln( ex_soln, grid, VL, VR)
  call initialize(grid,soln,VL,VR)

  open(50,file='temp.txt',status='unknown')

  100 format(1(I0.8),3(G20.7))

  write(header_str1,*) " Iter  |   ||density||    |"// &
  & "   ||velocity||    |   ||pressure||    |"
  call output_soln(grid,soln,ex_soln,0)
  call calc_vl_dfdu(soln%U(:,1),f,A,.true.)
  write(*,*) 'A+ = '
  do i = 1,neq
      write(*,'(*(F10.4))') (A(i,j), j = 1,3)
  end do
  write(*,*)

  call calc_vl_dfdu(soln%U(:,i_high),f,B,.false.)
  write(*,*) 'A- = '
  do i = 1,neq
      write(*,'(*(F10.4))') (B(i,j), j = 1,3)
  end do
  write(*,*)

  call flux_jac_cons1D(soln%U(:,1),C)
  write(*,*) 'A = '
  do i = 1,neq
      write(*,'(*(F10.4))') (C(i,j), j = 1,3)
  end do
  write(*,*)
  write(*,*) soln%U(:,1)

  stop

!!!  call calc_time_step(grid%dx,soln%V,soln%asnd,soln%lambda,soln%dt)
!!!  call build_LHS_matrix(soln,grid)
!!!  write(*,*) 'LHS(:,:,1,1)'
!!!  do i = 1,neq
!!!    do j = 1,neq
!!!      write(*,*) soln%LHS(i,j,1,1)
!!!    end do
!!!  end do
!!!  write(*,*)
!!!  write(*,*) 'LHS(:,:,1,2)'
!!!  do i = 1,neq
!!!    do j = 1,neq
!!!      write(*,*) soln%LHS(i,j,1,2)
!!!    end do
!!!  end do
!!!  write(*,*)
!!!  write(*,*) 'LHS(:,:,1,3)'
!!!  do i = 1,neq
!!!    do j = 1,neq
!!!      write(*,*) soln%LHS(i,j,1,3)
!!!    end do
!!!  end do
!!!  write(*,*)
!!!  write(*,*) 'LHS(:,:,1,4)'
!!!  do i = 1,neq
!!!    do j = 1,neq
!!!      write(*,*) soln%LHS(i,j,1,4)
!!!    end do
!!!  end do
!!!  write(*,*)
!!!  write(*,*) 'LHS(:,:,1,5)'
!!!  do i = 1,neq
!!!    do j = 1,neq
!!!      write(*,*) soln%LHS(i,j,1,5)
!!!    end do
!!!  end do
!!!  stop
!!!  !call explicit_euler(grid,soln%S,soln%dt,soln%F,soln%U,soln%R)
!!!  !call explicit_characteristic_bndry(soln, VL, VR)
!!!  !call explicit_RK(grid,soln,VL,VR)
!!!  call implicit_euler(grid,soln)
!!!  call update_states( soln )
!!!  write(*,*) minval(soln%dt)
!!!  soln%time = soln%time + minval(soln%dt)
!!!  call sample_exact_soln(ex_soln,grid,soln%time)
!!!  call output_soln(grid,soln,ex_soln,1)
!!!  call residual_norms(soln%R,soln%rinit,(/one,one,one/))
!!!  soln%rold = soln%rinit
!!!  write(*,*) 'Residual Norms: Iteration 0'
!!!  write(*,*) header_str1
!!!  write(*,100) j, soln%rinit(1), soln%rinit(2), soln%rinit(3)
!!!  write(*,*)
!!!  write(*,*) 'Relative Residual Norms:'
!!!  do j = 2,max_iter
!!!
!!!    !call calculate_sources(soln,grid)
!!!    call calc_time_step(grid%dx,soln%V,soln%asnd,soln%lambda,soln%dt)
!!!    call explicit_RK(grid,soln,VL,VR)
!!!
!!!    call implicit_euler(grid,soln)
!!!
!!!    soln%time = soln%time + minval(soln%dt)
!!!    call sample_exact_soln(ex_soln,grid,soln%time)
!!!
!!!   if (mod(j,soln_save)==0) then
!!!      call calc_de( soln, ex_soln, soln%DE, soln%DEnorm, pnorm, cons )
!!!      call output_soln(grid,soln,ex_soln,j)
!!!    end if
!!!    call residual_norms(soln%R,soln%rnorm,soln%rinit)
!!!
!!!    soln%rold = soln%rnorm
!!!    if (mod(j,10*res_out)==0) then
!!!      write(*,*) header_str1
!!!    end if
!!!    if (mod(j,res_out)==0) then
!!!      write(*,100) j, soln%rnorm(1), soln%rnorm(2), soln%rnorm(3)
!!!    end if
!!!    if (mod(j,res_save)==0) then
!!!      call output_res(soln,j)
!!!    end if
!!!  end do
!!!
!!!  if (mod(j,res_out)/=0) then
!!!    write(*,100) j, soln%rnorm(1), soln%rnorm(2), soln%rnorm(3)
!!!  end if
!!!  if (all(soln%rnorm<tol) ) then
!!!    write(*,*) 'Solution converged!'
!!!  else
!!!    write(niter,'(I7)') max_iter
!!!    write(*,*) 'Solution failed to converge in ',&
!!!       & trim(adjustl(niter)),' iterations'
!!!  end if
!!!  write(50,*) imax, j, soln%rnorm(1), soln%rnorm(2), soln%rnorm(3)
!!!
!!!  if (mod(j,res_out)/=0) then
!!!    call output_soln(grid,soln,ex_soln,j+1)
!!!    call output_res(soln,j)
!!!  end if
  call destroy_exact_soln( ex_soln )
  call teardown_geometry(grid,soln)
  close(50)
end program main_program
