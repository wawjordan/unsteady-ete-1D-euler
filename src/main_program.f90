program main_program

  use set_constants, only : set_derived_constants, zero, one
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs
  use set_inputs, only : max_iter, tol, soln_save, res_save
  use set_inputs, only : leftV, rightV, leftU, rightU, flux_scheme
  use set_inputs, only : limiter_freeze, res_out, cons
  use variable_conversion
  use time_integration
  use basic_boundaries, only : explicit_characteristic_bndry
  use limiter_calc, only : select_limiter
  use flux_calc, only : select_flux, flux_fun
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


  character(len=100) :: header_str1
  character(len=100) :: niter
  integer :: i, j, pnorm
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  type( exact_soln_t ) :: ex_soln
  real(prec), dimension(3) :: VL, VR
 ! type(dns_matrix_t) :: A_mat
 ! real(prec), dimension(3,3) :: array1, array2
 ! real(prec), dimension(3) :: b, x


  !array1 = transpose(reshape((/8.,1.,6.,3.,5.,7.,4.,9.,2./),shape(array1)))
  !array2 = one
  !x = zero

  !do i = 1,3
  !  b(i) = real(i,prec)
  !end do

  pnorm = 1
  j = 0

  VL = (/one,zero,one/)
  VR = (/0.125_prec,zero,0.1_prec/)

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
  call setup_geometry(grid,soln)

  call select_flux()
  call select_limiter()

  call setup_exact_soln( ex_soln, grid, VL, VR)
  call initialize(grid,soln,VL,VR)
  call output_soln(grid,soln,ex_soln,0)

  call explicit_characteristic_bndry(soln, VL, VR)
  call update_states( soln )
  call calculate_sources(soln%V(:,3),grid%dAc,soln%src)
  call calc_time_step(grid%dx,soln%V,soln%asnd,soln%lambda,soln%dt)

  call MUSCL_extrap( soln%V, leftV, rightV )
  call prim2cons(leftU,leftV)
  call prim2cons(rightU,rightV)
  call flux_fun(leftU,rightU,soln%F)

  call explicit_euler(grid,soln%src,soln%dt,soln%F,soln%U,soln%R)
  call update_states( soln )

  soln%time = soln%time + minval(soln%dt)

  call sample_exact_soln(ex_soln,grid,soln%time)

  call output_soln(grid,soln,ex_soln,1)

  call residual_norms(soln%R,soln%rinit,pnorm,(/one,one,one/))
  soln%rold = soln%rinit

  write(*,*) 'Residual Norms: Iteration 0'
  write(*,*) header_str1
  write(*,100) j, soln%rinit(1), soln%rinit(2), soln%rinit(3)
  write(*,*)
  write(*,*) 'Relative Residual Norms:'


  do j = 2,max_iter

    call explicit_characteristic_bndry(soln, VL, VR)
    call update_states( soln )
    call calculate_sources(soln%V(:,3),grid%dAc,soln%src)
    call calc_time_step(grid%dx,soln%V,soln%asnd,soln%lambda,soln%dt)

    call MUSCL_extrap( soln%V, leftV, rightV )
    call prim2cons(leftU,leftV)
    call prim2cons(rightU,rightV)
    call flux_fun(leftU,rightU,soln%F)

    call explicit_euler(grid,soln%src,soln%dt,soln%F,soln%U,soln%R)
    call update_states( soln )

    soln%time = soln%time + minval(soln%dt)
    call sample_exact_soln(ex_soln,grid,soln%time)

    if (mod(j,soln_save)==0) then
      call calc_de( soln, ex_soln, soln%DE, soln%DEnorm, pnorm, cons )
      call output_soln(grid,soln,ex_soln,j)
    end if

    call residual_norms(soln%R,soln%rnorm,pnorm,soln%rinit)

    if (all(soln%rnorm<tol) ) then
      exit
    elseif (any(soln%rnorm>soln%rold)) then
        limiter_freeze = .true.
    end if
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
