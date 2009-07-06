!==============================================================!
!  MOLECULAR SIMULATION PROGRAM MS2 Version 1.1 v12            !
!  (c) 2001 by Sergey Lishchuk, ITT                            !
!  (c) 2007 by Bernhard Eckl, ITT                              !
!  (c) 2007 by Ekaterina Elts, TUM                             !
!==============================================================!
!  Program ms2                                                 !
!  This file contains the main routine                         !
!==============================================================!

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
