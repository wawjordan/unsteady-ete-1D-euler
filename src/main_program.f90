program main_program

  use set_constants, only : set_derived_constants, zero, one
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs, xmin, xmax, neq
  use set_inputs, only : max_iter, tol, soln_save, res_save
  use set_inputs, only : leftV, rightV, leftU, rightU, flux_scheme
  use set_inputs, only : limiter_scheme, limiter_freeze, res_out, cons
  use variable_conversion
  use time_integration, only : explicit_Euler, calc_time_step, implicit_Euler
  use basic_boundaries, only : explicit_characteristic_bndry
  use limiter_type, only : select_limiter
  use flux_calc, only : select_flux, flux_fun, calc_flux_1D
  use flux_jacobians, only : calc_vl_dfdu, flux_jac_cons1D
  use other_subroutines
  use residuals, only : calc_residual, residual_norms
  use build_LHS_1st_order, only : build_LHS_matrix_1st_order
  use tridiag_operations, only : trivec_alt, tprec_alt
  use geometry, only : setup_geometry, teardown_geometry
  use init_problem, only : initialize
  use namelist, only : read_namelist
  !use utilities, only : block_band2full, block_band2sparse1
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

  i = 0
  j = 0
  pnorm = 1

  VL = (/one,0.1_prec,one/)
  VR = (/0.125_prec,0.1_prec,0.1_prec/)

  call set_derived_constants
  call set_fluid_constants
  call read_namelist('input.nml')
  call output_file_headers
  call set_derived_inputs
  call setup_geometry(grid,soln)
  call select_flux()
  call setup_exact_soln( ex_soln, grid, VL, VR)
  call initialize(grid,soln,VL,VR)

  open(50,file='temp.txt',status='unknown')

  100 format(1(I0.8),3(G20.7))

  write(header_str1,*) " Iter  |   ||density||    |"// &
  & "   ||velocity||    |   ||pressure||    |"
  !call output_soln(grid,soln,ex_soln,0)
  !call build_tridiag_test(soln)
  !call build_LHS_test(soln)
  call build_LHS_matrix_1st_order(soln,grid)
  !stop
 ! call block_band2full(soln%LHS,A_out,2,2,3,11)
 ! do i = 1,imax*neq
 !     write(*,*) (A_out(i,j), j = 1,imax*neq)
 ! end do

  call calc_residual(grid,soln)

  call calc_time_step(grid%dx,soln%V,soln%asnd,soln%lambda,soln%dt)
  !stop
!!!  !call explicit_euler(grid,soln%S,soln%dt,soln%F,soln%U,soln%R)
  call explicit_characteristic_bndry(soln, VL, VR)
  !stop
!!!  !call explicit_RK(grid,soln,VL,VR)
  call implicit_euler(grid,soln)
  !stop
  call update_states( soln )
  !stop
!!!  write(*,*) minval(soln%dt)
!!!  soln%time = soln%time + minval(soln%dt)
  !call sample_exact_soln(ex_soln,grid,soln%time)
  !stop
  call output_soln(grid,soln,ex_soln,1)
  call residual_norms(soln%R,soln%rinit,(/one,one,one/))
  soln%rold = soln%rinit
  write(*,*) 'Residual Norms: Iteration 0'
  write(*,*) header_str1
  write(*,100) j, soln%rinit(1), soln%rinit(2), soln%rinit(3)
  write(*,*)
  !stop
  write(*,*) 'Relative Residual Norms:'
  do j = 2,5 !max_iter
!!!    !call calculate_sources(soln,grid)
    call build_LHS_matrix_1st_order(soln,grid)
    call calc_time_step(grid%dx,soln%V,soln%asnd,soln%lambda,soln%dt)
    call implicit_euler(grid,soln)
    soln%time = soln%time + minval(soln%dt)
    !call sample_exact_soln(ex_soln,grid,soln%time)
   if (mod(j,soln_save)==0) then
      call calc_de( soln, ex_soln, soln%DE, soln%DEnorm, pnorm, cons )
      call output_soln(grid,soln,ex_soln,j)
    end if
    call residual_norms(soln%R,soln%rnorm,soln%rinit)
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
!!!
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
