!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 1.0                !
!  (c) 2011 by TU Kaiserslautern                               !
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

#ifndef TRANS
#define TRANS 0
#endif

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2.F90...'
#endif



program ms2

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
