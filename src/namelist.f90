module namelist

  use fluid_constants, only : R_gas, gamma
  use set_inputs, only : imax, i_high, i_low, ig_high, ig_low
  use set_inputs, only : neq, xmin, xmax, n_ghost
  use set_inputs, only : areaStar, area, darea
  use set_inputs, only : CFL, eps, tol, eps_roe, beta_lim, epsM, kappaM
  use set_inputs, only : max_iter, max_newton_iter, newton_tol
  use set_inputs, only : soln_save, res_save, res_out
  use set_inputs, only : p0, T0, a0, rho0
  use set_inputs, only : set_derived_inputs, flux_scheme, limiter_scheme, cons
  use set_inputs, only : leftV, rightV, leftU, rightU

  implicit none
  private

  public :: read_namelist

contains

  subroutine read_namelist(file_path)

    !use set_precision, only : prec
    !use set_constants, only : zero, one

    character(len=*), intent(in) :: file_path
    integer :: fstat
    integer :: funit
    logical :: fexist
    logical :: fopen = .false.
    namelist /grid/ imax, xmin, xmax, n_ghost
    namelist /geometry/ areaStar
    namelist /constants/ R_gas, gamma
    namelist /initial/ p0, T0
    namelist /numerical/ CFL, eps, tol, max_iter, &
            & max_newton_iter, newton_tol
    namelist /flux/ flux_scheme, limiter_scheme, &
            & eps_roe, beta_lim
    namelist /output/ soln_save, res_save, res_out, cons
    namelist /reconstruction/ epsM, kappaM

    inquire( file=file_path,exist=fexist )
    if ( .not. fexist ) then
      write(*,*) 'ERROR: ',trim(file_path),' does not exist.'
      return
    end if

    open(file=file_path,status='old',action='read',iostat=fstat,newunit=funit)
    fopen = .true.
    ! grid
    rewind(funit)
    read(funit,nml=grid,iostat=fstat)
    if (fstat > 0) then
      write(*,*) 'ERROR: error in namelist "grid".'
      stop
    end if

    ! geometry
    rewind(funit)
    read(funit,nml=geometry,iostat=fstat)
    if (fstat>0) then
      write(*,*) 'ERROR: error in namelist "geometry".'
      stop
    end if

    ! constants
    rewind(funit)
    read(funit,nml=constants,iostat=fstat)
    if (fstat>0) then
      write(*,*) 'ERROR: error in namelist "constants".'
      stop
    end if

    ! initial
    rewind(funit)
    read(funit,nml=initial,iostat=fstat)
    if (fstat>0) then
      write(*,*) 'ERROR: error in namelist "initial".'
      stop
    end if


    ! numerical
    rewind(funit)
    read(funit,nml=numerical,iostat=fstat)
    if (fstat>0) then
      write(*,*) 'ERROR: error in namelist "numerical".'
      stop
    end if

    ! flux
    rewind(funit)
    read(funit,nml=flux,iostat=fstat)
    if (fstat>0) then
      write(*,*) 'ERROR: error in namelist "flux".'
      stop
    end if

    ! output
    rewind(funit)
    read(funit,nml=output,iostat=fstat)
    if (fstat>0) then
      write(*,*) 'ERROR: error in namelist "output".'
      stop
    end if

    ! reconstruction
    rewind(funit)
    read(funit,nml=reconstruction,iostat=fstat)
    if (fstat>0) then
      write(*,*) 'ERROR: error in namelist "reconstruction".'
      stop
    end if

    close(funit)

  end subroutine read_namelist



end module namelist
