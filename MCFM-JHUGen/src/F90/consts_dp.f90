  module consts_dp
    implicit none

    double complex, public, parameter :: czero = (0d0,0d0)
    double complex, public, parameter :: ci = (0d0,1d0)
    double complex, public, parameter :: cone = (1d0,0d0)
    double complex, public, parameter :: two = 2d0
    double complex, public, parameter :: chalf = (0.5d0,0d0)
    double precision, public, parameter :: sqrt2 = &
       &1.4142135623730950488016887242096980785696718753769d0

    logical, public :: case_b2,qbq_and_gluons,qbq_WW_and_gluons,WWqqqq, &
       & case_a3,extra_ferm_pair1,case_b1

    double precision, public, parameter :: propcut = 1d-10
    double precision, public, parameter :: tol = 1d-14

  end module consts_dp