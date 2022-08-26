!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 2.0 + IDF          !
!  (c) 2014 by TU Kaiserslautern                               !
!      P.O. Box 67653                                          !
!      67653 Kaiserslautern                                    !
!==============================================================!
!  Program ms2                                                 !
!  This file contains the main routine                         !
!==============================================================!

!****************************************************************
!* Updates and auxiliary routines are available from            *   
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

  use ms2_simulation
  use ms2_global

  implicit none

    ! Include MPI header
#if MPI_VER > 0
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
