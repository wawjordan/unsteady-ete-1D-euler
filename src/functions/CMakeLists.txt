# Propagate lib_srcs list to parent
add_lib_sources(
     variable_conversion.f90
     time_integration.f90
     flux_calc.f90
     flux_jacobians.f90
     limiter_calc.f90
     other_subroutines.f90
     quadrature.f90
     riemann_problem.f90
     q1d_problem.f90
     utilities.f90
     extrapolation.f90
     block_matrix_operations.f90
)
