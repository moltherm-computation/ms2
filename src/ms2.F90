!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 4.0                !
!  (c) 2020 by TU Kaiserslautern / TU Berlin                   !
!      P.O. Box 67653                                          !
!      67653 Kaiserslautern                                    !
!==============================================================!
!  Program ms2                                                 !
!  This file contains the main routine                         !
!==============================================================!

!****************************************************************
!* Updates and auxiliary routines available from            *   
!* http://www.ms-2.de                                           *   
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

#if MPI_VER > 0
  use mpi_f08
#endif

  use ms2_simulation
  use ms2_global

  implicit none

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
