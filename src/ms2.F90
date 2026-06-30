!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 5.0                !
!  (c) 2026 by OVGU Magdeburg / TU Berlin                      !
!==============================================================!
!  Program ms2                                                 !
!  This file contains the main routine                         !
!==============================================================!

!****************************************************************
!* Source code and latest version:                              *
!* https://github.com/moltherm-computation/ms2                  *
!****************************************************************

#ifndef ARCH
#define ARCH    0
#define FORTRAN 90
#define MPI_VER 0
#endif

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2.F90...'
#endif


program ms2

#if MPI_VER > 0 && defined(MPI_USE_MODULE)
  use mpi_f08
#endif

  use ms2_simulation
  use ms2_global

  implicit none

    ! Include MPI header
#if MPI_VER > 0 && !defined(MPI_USE_MODULE)
    include 'mpif.h'
#endif

  ! Declare local variables
  type(TSimulation) :: Simulation

  ! Initialize program
  call InitializeProgram

  ! Create simulation object
  call Construct( Simulation )

  ! Run simulation
  call Run( Simulation )

  ! Destroy simulation object
  call Destruct( Simulation )

  ! Finalize program
  call FinalizeProgram

end program ms2
