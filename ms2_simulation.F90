!==============================================================!
!  MOLECULAR SIMULATION PROGRAM MS2 Version 1.1 v12            !
!  (c) 2001 by Sergey Lishchuk, ITT                            !
!  (c) 2007 by Bernhard Eckl, ITT                              !
!==============================================================!
!  Module ms2_simulation                                       !
!  Contains TSimulation object                                 !
!==============================================================!

#ifndef ARCH
#define ARCH    0
#define FORTRAN 90
#define MPI_VER 0
#endif

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2_simulation.F90...'
#endif

module ms2_simulation

  use ms2_ensemble
  use ms2_global



!==============================================================!
!  Type TSimulation                                            !
!==============================================================!

  type TSimulation

    ! Number of ensembles
    integer :: NEnsembles

    ! Ensembles
    type(TEnsemble), pointer :: Ensemble(:)

    ! I/O unit for result file
    integer :: iounit_result

    ! I/O unit for running average result file
    integer :: iounit_runave

    ! I/O unit for final result file
    integer :: iounit_errors

  end type TSimulation

  interface Construct
    module procedure TSimulation_Construct
  end interface

  interface Destruct
    module procedure TSimulation_Destruct
  end interface

  interface CreateAccumulators
    module procedure TSimulation_CreateAccumulators
  end interface

  interface DestroyAccumulators
    module procedure TSimulation_DestroyAccumulators
  end interface

  interface Run
    module procedure TSimulation_Run
  end interface

  interface RunSteps
    module procedure TSimulation_RunSteps
  end interface

  interface RunMDStep
    module procedure TSimulation_RunMDStep
  end interface

  interface RunMCStep
    module procedure TSimulation_RunMCStep
  end interface

  interface RunSVCStep
    module procedure TSimulation_RunSVCStep
  end interface

  interface CheckNPart
    module procedure TSimulation_CheckNPart
  end interface

  interface ResetEnsembles
    module procedure TSimulation_ResetEnsembles
  end interface

  interface ResultOpen
    module procedure TSimulation_ResultOpen
  end interface

  interface ResultUpdate
    module procedure TSimulation_ResultUpdate
  end interface

  interface ResultClose
    module procedure TSimulation_ResultClose
  end interface

  interface ErrorsUpdate
    module procedure TSimulation_ErrorsUpdate
  end interface

  interface SVCOutput
    module procedure TSimulation_SVCOutput
  end interface

  interface VisualOpen
    module procedure TSimulation_VisualOpen
  end interface

  interface VisualUpdate
    module procedure TSimulation_VisualUpdate
  end interface

  interface VisualClose
    module procedure TSimulation_VisualClose
  end interface

  interface RestartSave
    module procedure TSimulation_RestartSave
  end interface

  interface RestartRead
    module procedure TSimulation_RestartRead
  end interface



contains



!==============================================================!
!  Subroutine TSimulation_Construct                            !
!==============================================================!

  subroutine TSimulation_Construct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    character( IOBufferLength ) :: str
    integer                     :: i
    integer                     :: stat
    real(RK)                    :: KappaL_h
    real(RK)                    :: debyelen_h
    integer                     :: grid_h,spline_h
    integer                     :: nvecmax_h,nsqmax_h,nmax_h

    ! Read configuration file
#if ARCH == 1 || ARCH == 2 || ARCH == 3
    if( Restart ) then
      write( IOBuffer, '("Restarting from file: ", A)' ) &
&       trim( RestartFileName )
      call LogWrite
    else
      ParameterFileName = trim( OutputNameTag )//ParameterFileExtension
    end if
#else
    call FileReset( iounit_config, ProgramFileName//ConfigFileExtension )
    call FileReadParameter( iounit_config, IdRestart )
    read( IOBuffer, * ) str
    select case( str )
    case( 'YES', 'Yes', 'yes' )
      call FileReadParameter( iounit_config, IdRestartFileName )
      read( IOBuffer, * ) RestartFileName
      Restart = .true.
      write( IOBuffer, '("Restarting from file: ", A)' ) RestartFileName
      call LogWrite
      call FileReset( iounit_restart, RestartFileName )
      read( iounit_restart, '(A128)' ) ParameterFileName
    case( 'NO', 'No', 'no' )
      call FileReadParameter( iounit_config, IdParamsFileName )
      read( IOBuffer, * ) ParameterFileName
    case default
      call Error( 'Select yes/no for restart in file '// &
&       ProgramFileName//ConfigFileExtension )
    end select
#endif
    write( IOBuffer, '("Parameter file name: ", A)' ) trim( ParameterFileName )
    call LogWrite
#if ARCH != 1 && ARCH != 2 && ARCH != 3
    call FileClose( iounit_config )
#endif
    call LogWriteBlank

    ! Open parameter file for reading
    call FileReset( iounit_params, ParameterFileName )
    call LogWriteBlank
    write( IOBuffer, '("Reading parameters of simulation")' )
    call LogWrite

    ! Read name tag for output files
#if ARCH != 1 && ARCH != 2 && ARCH != 3
    call FileReadParameter( iounit_params , IdOutputNameTag )
    read( IOBuffer, * ) OutputNameTag
#endif
    write( IOBuffer, '("Name tag for output files: ", A)' ) &
&     trim( OutputNameTag )
    call LogWrite

    ! Read type of units
    call FileReadParameter( iounit_params , IdUseReducedUnits )
    read( IOBuffer, * ) str
    select case( str )
    case( 'REDUCED', 'Reduced', 'reduced' )
      UseReducedUnits = .true.
      str = 'reduced'
    case( 'SI', 'si' )
      UseReducedUnits = .false.
      str = 'SI'
    case default
      call Error( trim( str )//' system of units is not implemented' )
    end select
    write( IOBuffer, '("System of units: ", A)' ) trim( str )
    call LogWrite

    ! Read unit of length
    call FileReadParameter( iounit_params, IdUnitLength )
    read( IOBuffer, * ) UnitLength
    UnitLength = UnitLength * Angstroem
    write( IOBuffer, '("Unit of length: ", F6.3, " A")' ) &
&     UnitLength / Angstroem
    call LogWrite

    ! Read unit of energy
    call FileReadParameter( iounit_params, IdUnitEnergy )
    read( IOBuffer, * ) UnitEnergy
    UnitEnergy = UnitEnergy * kBoltzmann
    write( IOBuffer, '("Unit of energy: ", F8.3, " K")' ) &
&     UnitEnergy / kBoltzmann
    call LogWrite

    ! Read unit of mass
    call FileReadParameter( iounit_params, IdUnitMass )
    read( IOBuffer, * ) UnitMass
    UnitMass = UnitMass * .001_RK / NAvogadro
    write( IOBuffer, '("Unit of mass: ", F8.3, " a.u.")' ) &
&     UnitMass * NAvogadro * 1000._RK
    call LogWrite

    ! Calculate derived reduced units
    UnitVolume = UnitLength**3
    UnitTemperature = UnitEnergy / kBoltzmann
    UnitDensity = .001_RK / NAvogadro / UnitVolume
    UnitTime = sqrt(UnitMass / UnitEnergy) * UnitLength
    UnitForce = UnitEnergy / UnitLength
    UnitTorque = UnitEnergy
    UnitPressure = UnitForce / UnitLength**2
    UnitInertia = UnitMass * UnitLength**2
    UnitCharge = sqrt( 4._RK * Pi * VacuumPermittivity &
&     * UnitLength * UnitEnergy )
    UnitDipole = UnitCharge * UnitLength
    UnitQuadrupole = UnitCharge * UnitLength**2

    ! Read type of simulation
    call FileReadParameter( iounit_params , IdSimulationType )
    read( IOBuffer, * ) str
    select case( str )
    case( 'MD', 'md' )
      SimulationType = MolecularDynamics
      SimulationTypeString = 'Molecular Dynamics'
    case( 'MC', 'mc' )
      SimulationType = MonteCarlo
      SimulationTypeString = 'Monte-Carlo'
    case( 'SVC', 'svc', '2VC', '2vc' )
      SimulationType = SecondVirialCoeff
      SimulationTypeString = 'Second Virial Coefficient'
    case default
      call Error( trim( str )//' simulation is not implemented' )
    end select
    write( IOBuffer, '("Simulation type: ", A)' ) trim( SimulationTypeString )
    call LogWrite

    ! Read parameters specific to given simulation type
    if( SimulationType .eq. SecondVirialCoeff ) then

      ! Read number of orientations
      call FileReadParameter( iounit_params , IdNOrient )
      read( IOBuffer, * ) NOrient
      write( IOBuffer, '("Number of orientations: ", I7)' ) NOrient
      call LogWrite

      ! Read number of steps
      call FileReadParameter( iounit_params , IdRSteps )
      read( IOBuffer, * ) NSteps
      write( IOBuffer, '("Number of radial steps: ", I7)' ) NSteps
      call LogWrite

      ! Read minimum radius
      call FileReadParameter( iounit_params , IdMinRadius )
      read( IOBuffer, * ) MinRadius
      if( .not. UseReducedUnits ) then
        MinRadius = MinRadius / UnitLength * Angstroem
      end if
      write( IOBuffer, '("Minimum radius: ", F8.3, " A")' ) &
&       MinRadius * UnitLength / Angstroem
      call LogWrite

      ! Read maximum radius
      call FileReadParameter( iounit_params , IdMaxRadius )
      read( IOBuffer, * ) MaxRadius
      if( .not. UseReducedUnits ) then
        MaxRadius = MaxRadius / UnitLength * Angstroem
      end if
      write( IOBuffer, '("Maximum radius: ", F8.3, " A")' ) &
&       MaxRadius * UnitLength / Angstroem
      call LogWrite

      ! Set output frequencies
      BlockSize = 0
      ErrorsUpdateFrequency = NSteps
      VisualUpdateFrequency = 0

      ! Set cutoff mode
      CutoffMode = CenterofMass

    else

      ! Read MD simulation parameters
      if( SimulationType .eq. MolecularDynamics ) then

        ! Type of integrator
        call FileReadParameter( iounit_params , IdIntegratorType )
        read( IOBuffer, * ) str
        select case( str )
        case( 'GEAR', 'Gear', 'gear' )
          IntegratorType = IntegratorTypeGear
          IntegratorTypeString = 'Gear predictor-corrector'
        case( 'LEAPFROG', 'LeapFrog', 'Leapfrog', 'leapfrog' )
          IntegratorType = IntegratorTypeLeapFrog
          IntegratorTypeString = 'LeapFrog'
        case( 'VERLET', 'Verlet', 'verlet' )
          IntegratorType = IntegratorTypeVerlet
          IntegratorTypeString = 'Verlet'
        case( 'VV', 'Vv', 'vV', 'vv', &
&         'VELOCITY VERLET', 'Velocity Verlet', 'velocity Verlet', &
&         'velocity verlet', 'VELOCITY-VERLET', 'Velocity-Verlet', &
&         'velocity-Verlet', 'velocity-verlet', 'VELOCITYVERLET', &
&         'VelocityVerlet', 'velocityVerlet', 'velocityverlet' )
          IntegratorType = IntegratorTypeVV
          IntegratorTypeString = 'Velocity-Verlet'
        case default
          call Error( trim( str )//' integrator is not implemented' )
        end select
        write( IOBuffer, '("Integrator type: ", A)' ) &
&         trim( IntegratorTypeString )
        call LogWrite

        ! Time step
        call FileReadParameter( iounit_params , IdTimeStep )
        read( IOBuffer, * ) TimeStep
        write( IOBuffer, '("Time step: ", F9.6, " fs")' ) &
&         TimeStep * UnitTime * 1E15_RK
        call LogWrite
        write( IOBuffer, '("Reduced time step: ", F8.6)' ) TimeStep
        call LogWrite
        TimeStep2 = .5_RK * TimeStep
        TimeStepSquared = TimeStep**2
        TimeStepSquared2 = .5_RK * TimeStepSquared
        TimeStepSquaredInv2 = .5_RK / TimeStepSquared

      ! Read MC simulation parameters
      else

        ! Acceptance rate
        call FileReadParameter( iounit_params , IdAcceptance )
        read( IOBuffer, * ) Acceptance
        if( Acceptance < 0.05_RK ) then
          Acceptance = 0.05_RK
        else if( Acceptance > 0.95_RK ) then
          Acceptance = 0.95_RK
        end if
        write( IOBuffer, '("Acceptance rate: ", F6.2, "%")' ) &
&         Acceptance * 100._RK
        call LogWrite
        AccUpperLimit = Acceptance * 1.1_RK
        AccLowerLimit = Acceptance * 0.9_RK

      end if

      ! Read type of ensembles
      call FileReadParameter( iounit_params , IdEnsembleType )
      read( IOBuffer, * ) str
      select case( str )
      case( 'NVE', 'nve' )
        EnsembleType = EnsembleTypeNVE
        ConstantTemperature = .false.
        ConstantPressure = .false.
        EnsembleTypeString = 'NVE'
      case( 'NVT', 'nvt' )
        EnsembleType = EnsembleTypeNVT
        ConstantTemperature = .true.
        ConstantPressure = .false.
        EnsembleTypeString = 'NVT'
      case( 'NPH', 'nph' )
        EnsembleType = EnsembleTypeNPH
        ConstantTemperature = .false.
        ConstantPressure = .true.
        EnsembleTypeString = 'NPH'
      case( 'NPT', 'npt' )
        EnsembleType = EnsembleTypeNPT
        ConstantTemperature = .true.
        ConstantPressure = .true.
        EnsembleTypeString = 'NPT'
      case( 'GE', 'ge' )
        EnsembleType = EnsembleTypeGE
        ConstantTemperature = .true.
        ConstantPressure = .false.
        EnsembleTypeString = 'GrandEquilibrium'
      case( 'HA', 'ha' )
        EnsembleType = EnsembleTypeHA
        ConstantTemperature = .true.
        ConstantPressure = .true.
        EnsembleTypeString = 'Humid Air'
      case default
        call Error( trim( str )//' ensemble is not implemented' )
      end select
      write( IOBuffer, '("Ensemble type: ", A)' ) trim( EnsembleTypeString )
      call LogWrite

      ! Check whether simulation type is applicable to ensemble type
      if( SimulationType .eq. MonteCarlo .and. .not. ConstantTemperature ) &
&       call Error( trim( SimulationTypeString )//" simulation of " &
&         //trim( EnsembleTypeString )//" ensemble is not implemented" )

      if( (EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA) &
&         .and. .not. SimulationType .eq. MonteCarlo ) &
&       call Error( trim( SimulationTypeString )//" simulation of " &
&         //trim( EnsembleTypeString )//" ensemble is not implemented" )

      ! Read number of MC overlap reduction steps
      if( SimulationType .eq. MolecularDynamics ) then
        call FileReadParameter( iounit_params , IdNStepsMC )
        read( IOBuffer, * ) NStepsMC
        if( NStepsMC > 0 ) then
          write( IOBuffer, '("Number of MC overlap reduction steps: ", I7)' ) &
&           NStepsMC
          call LogWrite
          MCOverlapReduction = .true.
          Acceptance = .5_RK
          AccUpperLimit = Acceptance * 1.1_RK
          AccLowerLimit = Acceptance * 0.9_RK
        else
          write( IOBuffer, '("No MC overlap reduction")' )
          call LogWrite
          MCOverlapReduction = .false.
        end if
      else
        MCOverlapReduction = .false.
      end if

      ! Read number of NVT equilibration steps
      call FileReadParameter( iounit_params , IdNStepsV )
      read( IOBuffer, * ) NStepsV
      write( IOBuffer, '("Number of NVT equilibration steps: ", I7)' ) &
&       NStepsV
      call LogWrite

      ! Read number of NPT equilibration steps
      if( ConstantPressure ) then
        if( EnsembleType .eq. EnsembleTypeHA ) then
          call FileReadParameter( iounit_params , IdNStepsMueP )
          read( IOBuffer, * ) NStepsP
          write( IOBuffer, '("Number of HA equilibration steps: ", I7)' ) &
&           NStepsP
          call LogWrite
        else
          call FileReadParameter( iounit_params , IdNStepsP )
          read( IOBuffer, * ) NStepsP
          write( IOBuffer, '("Number of NPT equilibration steps: ", I7)' ) &
&           NStepsP
          call LogWrite
        end if
      else if( EnsembleType .eq. EnsembleTypeGE ) then
        call FileReadParameter( iounit_params , IdNStepsMue )
        read( IOBuffer, * ) NStepsP
        write( IOBuffer, '("Number of GE equilibration steps: ", I7)' ) &
&         NStepsP
        call LogWrite
      else
        NStepsP = 0
      end if

      ! Read number of production steps
      call FileReadParameter( iounit_params , IdNSteps )
      read( IOBuffer, * ) NSteps
      write( IOBuffer, '("Number of production steps: ", I7)' ) NSteps
      call LogWrite

      ! Read frequency of updating result file
      call FileReadParameter( iounit_params , IdBlockSize )
      read( IOBuffer, * ) BlockSize
      if( BlockSize > 0 ) then
        write( IOBuffer, &
&         '("Result files will be updated each", I7, " time steps")' ) &
&         BlockSize
      else
        write( IOBuffer, '("Result files will not be created")' )
      end if
      call LogWrite

      ! Calculate number of blocks and block sizes
      if( BlockSize > 0 ) then
        NBlocksMax = max( NStepsV, NStepsP, NSteps ) / BlockSize
        NBlockSizesMax = int( sqrt( real( NSteps / BlockSize, RK ) ) )
      else
        NBlocksMax = 0
        NBlockSizesMax = 0
      end if

      ! Read frequency of updating final result file
      call FileReadParameter( iounit_params , IdErrorsUpdateFrequency )
      read( IOBuffer, * ) ErrorsUpdateFrequency
      if( ErrorsUpdateFrequency < 1 ) then
        ErrorsUpdateFrequency = NSteps
      else if( ErrorsUpdateFrequency < BlockSize * 4 ) then
        ErrorsUpdateFrequency = BlockSize * 4
      end if
      if( ErrorsUpdateFrequency < NSteps ) then
        write( IOBuffer, &
&        '("Final result files will be updated each", I7, " time steps")' ) &
&         ErrorsUpdateFrequency
      else
        write( IOBuffer, '("Final result files will be created at the end")' )
      end if
      call LogWrite

      ! Read frequency of updating visualisation file
      call FileReadParameter( iounit_params , IdVisualUpdateFrequency )
      read( IOBuffer, * ) VisualUpdateFrequency
      if( VisualUpdateFrequency > 0 ) then
        write( IOBuffer, &
&        '("Visualization files will be updated each", I7, " time steps")' ) &
&         VisualUpdateFrequency
      else
        write( IOBuffer, '("Visualization files will not be created")' )
      end if
      call LogWrite

      ! Read cutoff mode
      call FileReadParameter( iounit_params , IdCutoffMode )
      read( IOBuffer, * ) str
      select case( str )
      case( 'COM', 'com', 'CenterOfMass', 'CenterofMass', 'centerofmass' )
        CutoffMode = CenterofMass
        CutoffModeString = 'Center of Mass'
        write( IOBuffer, '("Cutoff mode: ", A)' ) trim( CutoffModeString )
      case( 'Site', 'site', 'Site-Site', 'site-site' )
        CutoffMode = SiteSite
        CutoffModeString = 'Site-Site'
      case default
        call Error( trim( str )//' is not a valid cutoff mode' )
      end select
      call LogWrite

      ! Read LongRange mode
      call FileReadParameter( iounit_params , IdLongRange )
      read( IOBuffer, * ) str
      select case( str )
        case( 'Ewald', 'ew', 'ewald', 'EWALD')
            LongRange = Ewald
            LongRangeString = 'EwaldSum'
            write( IOBuffer, '("Long Range Correction: ", A)' ) trim( LongRangeString )
            call LogWrite
            ! Read Ewald Parameters
            call FileReadParameter( iounit_params , IdKappa )
            read( IOBuffer, * ) KappaL_h
            write( IOBuffer, '("Ewald Parameter KappaL:", F8.3)' )KappaL_h
            call LogWrite

            call FileReadParameter( iounit_params , Idnsqmax )
            read( IOBuffer, * ) nsqmax_h
            write( IOBuffer, '("Ewald Parameter NsqMax:", I7)' ) nsqmax_h
            call LogWrite

            call FileReadParameter( iounit_params , IdNVecMax )
            read( IOBuffer, * ) nvecmax_h

            write( IOBuffer, '("Ewald Parameter NVecMax:", I7)' ) nvecmax_h
            call LogWrite
            call FileReadParameter( iounit_params , IdNMax )
            read( IOBuffer, * ) nmax_h

            write( IOBuffer, '("Ewald Parameter NMax:", I7)' ) nmax_h
            call LogWrite

        case( 'PME', 'pme', 'SPME', 'spme')
            LongRange = PME
            LongRangeString = 'Smooth Particle Mesh Ewald Summation'
            write( IOBuffer, '("Long Range Correction: ", A)' ) trim( LongRangeString )
            call LogWrite
            ! Read SPM Ewald Parameters
            call FileReadParameter( iounit_params , IdKappa )
            read( IOBuffer, * ) KappaL_h
            write( IOBuffer, '("Ewald Parameter KappaL:", F8.3)' )KappaL_h
            call LogWrite

            call FileReadParameter( iounit_params , IdGrid )
            read( IOBuffer, * ) grid_h
            write( IOBuffer, '("Grid Space SPME:", I7)' ) grid_h
            call LogWrite

            call FileReadParameter( iounit_params , IdSpline )
            read( IOBuffer, * ) spline_h
            write( IOBuffer, '("order of SPME Spline:", I7)' ) spline_h

      case( 'ReactionField', 'RF', 'reactionfield', 'rf' )
         LongRange = RField
        LongRangeString = 'Reaction Field'
        write( IOBuffer, '("Long Range Correction: ", A)' ) trim( LongRangeString )

      case( 'ExtReactionField', 'ExtRF', 'extreactionfield', 'extrf' )
         LongRange = ExtRField
        LongRangeString = 'Extended Reaction Field by Tironi et al.'
        write( IOBuffer, '("Long Range Correction: ", A)' ) trim( LongRangeString )
            ! Read extended Reaction Field Parameters
            call FileReadParameter( iounit_params , IdDebyeLen )
            read( IOBuffer, * ) debyelen_h
            write( IOBuffer, '("Debye Length [A]:", F8.3)' )debyelen_h
            call LogWrite

      case default
        call Error( trim( str )//' is not a valid cutoff mode' )
      end select
      call LogWrite

    end if

    ! Read number of ensembles
    call FileReadParameter( iounit_params , IdNEnsembles )
    read( IOBuffer, * ) this%NEnsembles
    write( IOBuffer, '("Number of ensembles:", I3)' ) this%NEnsembles
    call LogWrite

    ! Create ensembles
    allocate( this%Ensemble(this%NEnsembles), STAT = stat )
    call AllocationError( stat, 'ensembles', this%NEnsembles )
    do i = 1, this%NEnsembles
      if (LongRange .eq. Ewald) then
            this%ensemble(i)%KappaL = KappaL_h
            this%ensemble(i)%nsqmax = nsqmax_h
            this%ensemble(i)%nvecmax = nvecmax_h
            this%ensemble(i)%nmax = nmax_h
      else if (LongRange .eq. PME) then
            this%ensemble(i)%KappaL = KappaL_h
            this%ensemble(i)%gridx  = grid_h
            this%ensemble(i)%gridy  = grid_h
            this%ensemble(i)%gridz  = grid_h
            this%ensemble(i)%splineorder = spline_h
            allocate(this%ensemble(i)%qgrida(2,(grid_h)**3+1),STAT=stat)
            if(stat >0) write(*,*) 'Allocation Error grida'
            allocate(this%ensemble(i)%qgrida_old(2,(grid_h)**3+1),STAT=stat)
            if(stat >0) write(*,*) 'Allocation Error grida_old'
            allocate(this%ensemble(i)%qgridb(2,(grid_h)**3+1),STAT=stat)
            if(stat >0) write(*,*) 'Allocation Error gridb'
      else if (LongRange .eq. ExtRField) then
            this%ensemble(i)%DebyeLen = debyelen_h * Angstroem / UnitLength

      end if
      if( SimulationType .eq. SecondVirialCoeff ) then
        call ConstructSVC( this%Ensemble(i), i )
      else
        call Construct( this%Ensemble(i), i )
      end if
!       if (LongRange .eq. PME) then
!         call PMESelfTerm(this%Ensemble(i))
!       endif
    end do

!DEBUG
!   if( any( this%Ensemble(:)%Component(:)%ChemPotMethod .eq. ChemPotMethodGradIns ) then
!     call FileReadParameter( iounit_params, IDFluctFreq )
!     read( IOBuffer, * ) GradInsFrequency
!     call FileReadParameter( iounit_params, IDNFullFluct )
!     read( IOBuffer, * ) NFullFluct
!     call FileReadParameter( iounit_params, IDMaxCounter )
!     read( IOBuffer, * ) maxcounter
!   end if
  GradInsFrequency = BlockSize
  NFullFluct = 20
  maxcounter = 0
!DEBUG

    ! Close parameter file
    call LogWriteBlank
    call FileClose( iounit_params )

    ! Create accumulators
    call CreateAccumulators( this )

    ! Set I/O unit numbers
    this%iounit_result = iounit_result
    this%iounit_runave = iounit_runave
    this%iounit_errors = iounit_errors

    ! Open result and visualisation files
    call LogWriteBlank
    call ResultOpen( this )
    call VisualOpen( this )

  end subroutine TSimulation_Construct



!==============================================================!
!  Subroutine TSimulation_Destruct                             !
!==============================================================!

  subroutine TSimulation_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Close result and visualisation files
    call LogWriteBlank
    call ResultClose( this )
    call VisualClose( this )

    ! Destroy accumulators
    call DestroyAccumulators( this )

    ! Destroy ensembles
    if( associated( this%Ensemble ) ) then
      do i = 1, this%NEnsembles
        call Destruct( this%Ensemble(i) )
      end do
      deallocate( this%Ensemble )
    end if

  end subroutine TSimulation_Destruct



!==============================================================!
!  Subroutine TSimulation_CreateAccumulators                   !
!==============================================================!

  subroutine TSimulation_CreateAccumulators( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Construct accumulators
    do i = 1, this%NEnsembles
      call CreateAccumulators( this%Ensemble(i) )
    end do

  end subroutine TSimulation_CreateAccumulators



!==============================================================!
!  Subroutine TSimulation_DestroyAccumulators                  !
!==============================================================!

  subroutine TSimulation_DestroyAccumulators( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Destruct accumulators
    do i = 1, this%NEnsembles
      call DestroyAccumulators( this%Ensemble(i) )
    end do


  end subroutine TSimulation_DestroyAccumulators



!==============================================================!
!  Subroutine TSimulation_Run                                  !
!==============================================================!

  subroutine TSimulation_Run( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: StepStart, StepEnd
    integer :: i
    logical :: NPartsOk

    tooManyParticles = .false.

    if( Restart ) then
      call RestartRead( this )
      StepStart = Step + 1
      MCOverlapReduction = .false.                                             ! no MC overlap reduction in case of restart
    else
      StepTotal = 0
      StepStart = 1
      if( SimulationType .eq. SecondVirialCoeff ) then
        Equilibration = .false.
        NVTEquilibration = .false.
      else
        Equilibration = .true.
        NVTEquilibration = .true.
      end if
    end if

    ! Run MC overlap reduction
    if( MCOverlapReduction .and. .not. TerminateProgram ) then
      StepEnd = NStepsMC
      call LogWriteBlank
      if( Restart ) then
        write( IOBuffer, '("Resuming MC overlap reduction")' )
        Restart = .false.
      else
        write( IOBuffer, '("Starting MC overlap reduction")' )
      end if
      call LogWriteTime
      SimulationType = MonteCarlo

      call RunSteps( this, StepStart, StepEnd )

      if( .not. TerminateProgram ) then
        write( IOBuffer, '("MC overlap reduction completed")' )
        MCOverlapReduction = .false.
        SimulationType = MolecularDynamics
        do i = 1, this%NEnsembles
          call InitMolecularDynamics( this%Ensemble(i), .true. )
        end do
      else
        write( IOBuffer, '("MC overlap reduction terminated")' )
      end if
      call LogWriteTime
      StepStart = 1
    end if

eqloop: do
      ! Run NVT equilibration
      if( NVTEquilibration .and. .not. TerminateProgram ) then
        StepEnd = NStepsV
        call LogWriteBlank
        if( Restart ) then
          write( IOBuffer, '("Resuming NVT equilibration")' )
          Restart = .false.
        else
          write( IOBuffer, '("Starting NVT equilibration")' )
        end if
        call LogWriteTime

        call RunSteps( this, StepStart, StepEnd )

        if( .not. TerminateProgram ) then
          write( IOBuffer, '("NVT equilibration completed")' )
          NVTEquilibration = .false.
        else
          write( IOBuffer, '("NVT equilibration terminated")' )
        end if
        call LogWriteTime
        StepStart = 1
      end if

      ! Run GE or NPT equilibration
      if( Equilibration .and. .not. TerminateProgram ) then
        StepEnd = NStepsP
        if( EnsembleType .eq. EnsembleTypeGE ) then
          call LogWriteBlank
          if( Restart ) then
            write( IOBuffer, '("Resuming GE equilibration")' )
            Restart = .false.
          else
            write( IOBuffer, '("Starting GE equilibration")' )
          end if
          call LogWriteTime

          call RunSteps( this, StepStart, StepEnd )

          if( .not. TerminateProgram ) then
            call CheckNPart( this, NPartsOk )
            if( NPartsOk ) then
              write( IOBuffer, '("GE equilibration completed")' )
              Equilibration = .false.
            else
              write( IOBuffer, &
&               '("GE equilibration ended with too many/too less particles")' )
              call LogWriteTime
              write( IOBuffer, '("Restarting equilibration")' )
              call LogWrite
              call ResetEnsembles( this )
              tooManyParticles = .false.
              NVTEquilibration = .true.
              StepStart = 1
              cycle eqloop
            end if
          else
            write( IOBuffer, '("GE equilibration terminated")' )
          end if
          call LogWriteTime

        else if( EnsembleType .eq. EnsembleTypeHA ) then
          call LogWriteBlank
          if( Restart ) then
            write( IOBuffer, '("Resuming HA equilibration")' )
            Restart = .false.
          else
            write( IOBuffer, '("Starting HA equilibration")' )
          end if
          call LogWriteTime

          call RunSteps( this, StepStart, StepEnd )

          if( .not. TerminateProgram ) then
            call CheckNPart( this, NPartsOk )
            if( NPartsOk ) then
              write( IOBuffer, '("HA equilibration completed")' )
              Equilibration = .false.
            else
              write( IOBuffer, &
&               '("HA equilibration ended with too many/too less particles")' )
              call LogWriteTime
              write( IOBuffer, '("Restarting equilibration")' )
              call LogWrite
              call ResetEnsembles( this )
              tooManyParticles = .false.
              NVTEquilibration = .true.
              StepStart = 1
              cycle eqloop
            end if
          else
            write( IOBuffer, '("HA equilibration terminated")' )
          end if
          call LogWriteTime

        else if( ConstantPressure ) then
          call LogWriteBlank
          if( Restart ) then
            write( IOBuffer, '("Resuming NPT equilibration")' )
            Restart = .false.
          else
            write( IOBuffer, '("Starting NPT equilibration")' )
          end if
          call LogWriteTime

          call RunSteps( this, StepStart, StepEnd )

          if( .not. TerminateProgram ) then
            write( IOBuffer, '("NPT equilibration completed")' )
            Equilibration = .false.
          else
            write( IOBuffer, '("NPT equilibration terminated")' )
          end if
          call LogWriteTime

        else
          Equilibration = .false.
        end if
        StepStart = 1
      end if

      exit eqloop
    end do eqloop

    ! GradIns initialization
!     if( GradInsInitialization .and. .not. TerminateProgram ) then
!       StepEnd = NStepsG
!       call LogWriteBlank
!       if( Restart ) then
!         write( IOBuffer, '("Resuming GradIns initialization")' )
!         Restart = .false.
!       else
!         write( IOBuffer, '("Starting GradIns initialization")' )
!         call LogWrite
!         write( IOBuffer, '("  (adjustment of weigthing factors)")' )
!       end if
!       call LogWriteTime
! 
!       do i = 1, this%NEnsembles
!         call GradInsInit( this%Ensemble(i) )
!       end do
! 
!       if( .not. TerminateProgram ) then
!         write( IOBuffer, '("GradIns initialization completed")' )
!         GradInsInitialization = .false.
!       else
!         write( IOBuffer, '("GradIns initialization terminated")' )
!       end if
!       call LogWriteTime
!       StepStart = 1
!     end if

    ! Run production
    if( .not. TerminateProgram ) then
      StepEnd = NSteps
      call LogWriteBlank
      if( Restart ) then
        write( IOBuffer, '("Resuming simulation")' )
        Restart = .false.
      else
        write( IOBuffer, '("Starting simulation")' )
      end if
      call LogWriteTime

      call RunSteps( this, StepStart, StepEnd )

      if( .not. TerminateProgram ) then
        write( IOBuffer, '("Simulation completed")' )
!         GradInsInitialization = .false.
      else
        write( IOBuffer, '("Simulation terminated")' )
      end if
      call LogWriteTime
    end if

    ! Output for second virial coefficient run
    if( SimulationType .eq. SecondVirialCoeff ) call SVCOutput( this )

    ! Save restart file
    call LogWriteBlank
    write( IOBuffer, '("Saving simulation restart file")' )
    call LogWrite
    if( Step > StepEnd ) then
      Step = StepEnd
      if( BlockSize > 0 ) NBlocks = 1 + (Step - 1) / BlockSize
    end if
    call RestartSave( this )

  end subroutine TSimulation_Run



!==============================================================!
!  Subroutine TSimulation_RunSteps                             !
!==============================================================!

  subroutine TSimulation_RunSteps( this, StepStart, StepEnd )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation)   :: this
    integer, intent(in) :: StepStart, StepEnd

    ! Declare local variables
#if MPI_VER > 0
    logical :: AnyTerminateProgram
#endif

    ! Run simulation steps
    do Step = StepStart, StepEnd

      ! Update total number of steps
      StepTotal = StepTotal + 1

      ! Set current block number
      if( BlockSize > 0 ) then
        NBlocks = 1 + (Step - 1) / BlockSize
        NBlockSizes = int( sqrt( real( Step / BlockSize, RK ) ) )
      end if

      ! Run simulation step
      select case( SimulationType )
      case( MolecularDynamics )
        call RunMDStep( this )
      case( MonteCarlo )
        call RunMCStep( this )
      case( SecondVirialCoeff )
        call RunSVCStep( this )
      end select

      ! Update result and visualisation files
      call ResultUpdate( this )
      call VisualUpdate( this )

      ! Update log and result files
      if( mod( Step, LogUpdateFrequency ) == 0 .or. Step == StepEnd ) &
&       call LogWriteStep
      if( .not. Equilibration .and. &
&       ( mod( Step, ErrorsUpdateFrequency ) == 0 .or. Step == StepEnd )) &
&       call ErrorsUpdate( this )

      ! Check for termination request (caused by signal handler)
#if MPI_VER > 0 && ( ARCH == 1 || ARCH == 2 )
      call MPI_Allreduce( TerminateProgram, AnyTerminateProgram, 1, &
&       MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror )
      if( AnyTerminateProgram ) then
        TerminateProgram = .true.
        exit
      end if
#else
      if( TerminateProgram ) exit
#endif

      ! Check for too many particles (GE only)
      if( tooManyParticles ) exit

    end do

  end subroutine TSimulation_RunSteps



!==============================================================!
!  Subroutine TSimulation_RunMDStep                            !
!==============================================================!

  subroutine TSimulation_RunMDStep( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Run MD simulation step
    do i = 1, this%NEnsembles
      call RunMDStep( this%Ensemble(i) )
    end do

  end subroutine TSimulation_RunMDStep



!==============================================================!
!  Subroutine TSimulation_RunMCStep                            !
!==============================================================!

  subroutine TSimulation_RunMCStep( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Run MC simulation step
    do i = 1, this%NEnsembles
      call RunMCStep( this%Ensemble(i) )
    end do

  end subroutine TSimulation_RunMCStep



!==============================================================!
!  Subroutine TSimulation_RunSVCStep                           !
!==============================================================!

  subroutine TSimulation_RunSVCStep( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Run SVC simulation step
    do i = 1, this%NEnsembles
      call RunSVCStep( this%Ensemble(i) )
    end do

  end subroutine TSimulation_RunSVCStep



!==============================================================!
!  Subroutine TSimulation_CheckNPart                           !
!==============================================================!

  subroutine TSimulation_CheckNPart( this, NPartsOk )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation)    :: this
    logical, intent(out) :: NPartsOk

    ! Declare local variables
    integer :: i

    ! Check number of particles in ensembles
    NPartsOk = .true.
    do i = 1, this%NEnsembles
      call CheckNPart( this%Ensemble(i), NPartsOk )
    end do

  end subroutine TSimulation_CheckNPart



!==============================================================!
!  Subroutine TSimulation_ResetEnsembles                       !
!==============================================================!

  subroutine TSimulation_ResetEnsembles( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Run MC simulation step
    do i = 1, this%NEnsembles
      call ResetEnsemble( this%Ensemble(i) )
    end do

  end subroutine TSimulation_ResetEnsembles



!==============================================================!
!  Subroutine TSimulation_ResultOpen                           !
!==============================================================!

  subroutine TSimulation_ResultOpen( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Return if no output
    if( BlockSize < 1 .and. .not. SimulationType .eq. SecondVirialCoeff ) return

!     if( this%NEnsembles > 1 ) then
! 
!       if( Restart ) then
!         ! Open result file
!         call FileAppend( this%iounit_result, &
! &         trim( OutputNameTag )//ResultFileExtension )
! 
!         ! Open running average result file
!         call FileAppend( this%iounit_runave, &
! &         trim( OutputNameTag )//RunAveFileExtension )
! 
!       else
!         ! Open result file
!         call FileRewrite( this%iounit_result, &
! &         trim( OutputNameTag )//ResultFileExtension )
! 
!         ! Open running average result file
!         call FileRewrite( this%iounit_runave, &
! &         trim( OutputNameTag )//RunAveFileExtension )
! 
!       end if
! 
!     end if

    ! Open ensemble result files
    do i = 1, this%NEnsembles
      call ResultOpen( this%Ensemble(i) )
    end do

  end subroutine TSimulation_ResultOpen



!==============================================================!
!  Subroutine TSimulation_ResultUpdate                         !
!==============================================================!

  subroutine TSimulation_ResultUpdate( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Return if no output
    if( BlockSize < 1 ) return

    ! No output for MCOverlapReduction
    if( MCOverlapReduction ) return

!     ! Generate simulation result files
!     if( this%NEnsembles > 1 ) then
! 
!       ! Update result header
!       if( Step == 1 ) then
!         call FileWriteBlank( this%iounit_result )
!         call FileWriteBlank( this%iounit_runave )
! 
!         ! Number of steps
!         write( IOBuffer, '("     NR")' )
!         call FileWriteNoAdvance( this%iounit_result )
!         call FileWriteNoAdvance( this%iounit_runave )
! 
!         call FileWriteBlank( this%iounit_result )
!         call FileWriteBlank( this%iounit_runave )
!       end if
! 
!       ! Update result files
!       if( mod( Step, BlockSize ) == 0 ) then
! 
!         ! Number of steps
!         write( IOBuffer, '(I7)' ) Step
!         call FileWriteNoAdvance( this%iounit_result )
!         call FileWriteNoAdvance( this%iounit_runave )
! 
!         call FileWriteBlank( this%iounit_result )
!         call FileWriteBlank( this%iounit_runave )
!       end if
! 
!     end if

    ! Update ensemble result files
    do i = 1, this%NEnsembles
      call ResultUpdate( this%Ensemble(i) )
    end do

  end subroutine TSimulation_ResultUpdate



!==============================================================!
!  Subroutine TSimulation_ResultClose                          !
!==============================================================!

  subroutine TSimulation_ResultClose( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Return if no output
    if( BlockSize < 1 .and. .not. SimulationType .eq. SecondVirialCoeff ) return

    ! Close ensemble result files
    do i = 1, this%NEnsembles
      call ResultClose( this%Ensemble(i) )
    end do

!     if( this%NEnsembles > 1 ) then
! 
!       ! Close running average result file
!       call FileClose( this%iounit_runave )
! 
!       ! Close result file
!       call FileClose( this%iounit_result )
! 
!     end if

  end subroutine TSimulation_ResultClose



!==============================================================!
!  Subroutine TSimulation_ErrorsUpdate                         !
!==============================================================!

  subroutine TSimulation_ErrorsUpdate( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Return if no output
    if( BlockSize < 1 ) return

    ! Update log file
    write( IOBuffer, '("Saving simulation results")' )
    call LogWrite

!     if( this%NEnsembles > 1 ) then
! 
!       ! Open final result file
!       call FileRewrite( this%iounit_errors, &
! &       trim( OutputNameTag )//ErrorsFileExtension )
! 
!       ! Separator
!       write( IOBuffer, '(76("="))' )
!       call FileWrite( this%iounit_errors )
!       call FileWriteBlank( this%iounit_errors )
!       write( IOBuffer, '(T24, "SIMULATION RESULT FILE")' )
!       call FileWrite( this%iounit_errors )
!       write( IOBuffer, '(T24, "----------------------")' )
!       call FileWrite( this%iounit_errors )
!       call FileWriteBlank( this%iounit_errors )
! 
!       ! Simulation type
!       write( IOBuffer, '("Simulation type", T36, ":", 9X, A)' ) &
! &       trim( SimulationTypeString )
!       call FileWrite( this%iounit_errors )
!       write( IOBuffer, '("Ensemble type", T36, ":", 9X, A)' ) &
! &       trim( EnsembleTypeString )
!       call FileWrite( this%iounit_errors )
!       if( SimulationType .eq. MolecularDynamics ) then
!         write( IOBuffer, '("Integrator type", T36, ":", 9X, A)' ) &
! &         trim( IntegratorTypeString )
!         call FileWrite( this%iounit_errors )
!       end if
!       call FileWriteBlank( this%iounit_errors )
! 
!       ! Number of steps
!       write( IOBuffer, '("Number of NVT equilibration steps", T36, ":", I10)' ) &
! &       NStepsV
!       call FileWrite( this%iounit_errors )
!       write( IOBuffer, '("Number of NPT equilibration steps", T36, ":", I10)' ) &
! &       NStepsP
!       call FileWrite( this%iounit_errors )
!       write( IOBuffer, '("Number of production steps", T36, ":", I10)' ) &
! &       Step
!       call FileWrite( this%iounit_errors )
!       call FileWriteBlank( this%iounit_errors )
! 
!       ! Time step
!       if( SimulationType .eq. MolecularDynamics ) then
!         write( IOBuffer, '("Time step", T29, "reduced:", F20.9)' ) &
! &         TimeStep
!         call FileWrite( this%iounit_errors )
!         write( IOBuffer, '(T31, "in fs:", F20.9)' ) &
! &         TimeStep * UnitTime * 1E15_RK
!         call FileWrite( this%iounit_errors )
!         call FileWriteBlank( this%iounit_errors )
!       end if
! 
!       ! Acceptance rate
!       if( SimulationType .eq. MonteCarlo ) then
!         write( IOBuffer, '("Acceptance rate", T36, ":", F20.9)' ) &
! &         Acceptance
!         call FileWrite( this%iounit_errors )
!         call FileWriteBlank( this%iounit_errors )
!       end if
! 
!       ! Number of ensembles
!       write( IOBuffer, '("Number of ensembles", T36, ":", I10)' ) &
! &       this%NEnsembles
!       call FileWrite( this%iounit_errors )
!       call FileWriteBlank( this%iounit_errors )
! 
!       ! System of units
!       write( IOBuffer, '("Unit of length", T36, ":", F20.9, " A")' ) &
! &       UnitLength / Angstroem
!       call FileWrite( this%iounit_errors )
!       write( IOBuffer, '("Unit of energy", T36, ":", F20.9, " K")' ) &
! &       UnitEnergy / kBoltzmann
!       call FileWrite( this%iounit_errors )
!       write( IOBuffer, '("Unit of mass", T36, ":", F20.9, " a.u.")' ) &
! &       UnitMass * NAvogadro * 1000._RK
!       call FileWrite( this%iounit_errors )
!       call FileWriteBlank( this%iounit_errors )
! 
!       ! Separator
!       write( IOBuffer, '(76("="))' )
!       call FileWrite( this%iounit_errors )
!       call FileWriteBlank( this%iounit_errors )
!       write( IOBuffer, '("VALUE", T31, "UNITS", T46, "AVERAGE", T66, "ERROR")' )
!       call FileWrite( this%iounit_errors )
!       write( IOBuffer, '(76("-"))' )
!       call FileWrite( this%iounit_errors )
!       call FileWriteBlank( this%iounit_errors )
! 
!       ! Separator
!       write( IOBuffer, '(76("="))' )
!       call FileWrite( this%iounit_errors )
!       call FileWriteBlank( this%iounit_errors )
! 
!       ! Close final result file
!       call FileClose( this%iounit_errors )
! 
!     end if

    ! Save ensemble results
    do i = 1, this%NEnsembles
      call ErrorsUpdate( this%Ensemble(i) )
    end do

  end subroutine TSimulation_ErrorsUpdate



!==============================================================!
!  Subroutine TSimulation_SVCOutput                            !
!==============================================================!

  subroutine TSimulation_SVCOutput( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Update log file
    write( IOBuffer, '("Saving simulation results")' )
    call LogWrite

    ! Save ensemble results
    do i = 1, this%NEnsembles
      call SVCOutput( this%Ensemble(i) )
    end do

  end subroutine TSimulation_SVCOutput



!==============================================================!
!  Subroutine TSimulation_VisualOpen                           !
!==============================================================!

  subroutine TSimulation_VisualOpen( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Return if no output
    if( VisualUpdateFrequency < 1 ) return

    ! Open ensemble visualisation files
    do i = 1, this%NEnsembles
      call VisualOpen( this%Ensemble(i) )
    end do

  end subroutine TSimulation_VisualOpen



!==============================================================!
!  Subroutine TSimulation_VisualUpdate                         !
!==============================================================!

  subroutine TSimulation_VisualUpdate( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Return if no output
    if( VisualUpdateFrequency < 1 ) return

    ! Return if equilibration
    if( Equilibration ) return

    ! Update ensemble visualisation files
    if( mod( StepTotal - 1, VisualUpdateFrequency ) == 0 ) then
      do i = 1, this%NEnsembles
        call VisualUpdate( this%Ensemble(i) )
      end do
    end if

  end subroutine TSimulation_VisualUpdate



!==============================================================!
!  Subroutine TSimulation_VisualClose                          !
!==============================================================!

  subroutine TSimulation_VisualClose( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Return if no output
    if( VisualUpdateFrequency < 1 ) return

    ! Close ensemble visualisation files
    do i = 1, this%NEnsembles
      call VisualClose( this%Ensemble(i) )
    end do

  end subroutine TSimulation_VisualClose



!==============================================================!
!  Subroutine TSimulation_RestartSave                          !
!==============================================================!

  subroutine TSimulation_RestartSave( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

!DEBUG
  if( SimulationType .eq. SecondVirialCoeff ) return
!DEBUG

    ! Open restart file for writing
    call FileRewrite( iounit_restart, &
&     trim( OutputNameTag )//RestartFileExtension )

    ! Save contents to restart file
    write( iounit_restart, '(A)' ) trim( ParameterFileName )
    write( iounit_restart, '(2I10)' ) Step, StepTotal
    write( iounit_restart, '(2L5)' ) Equilibration, NVTEquilibration

    ! Save ensembles
    do i = 1, this%NEnsembles
      call RestartSave( this%Ensemble(i) )
    end do

    ! Close restart file
    call FileClose( iounit_restart )

  end subroutine TSimulation_RestartSave



!==============================================================!
!  Subroutine TSimulation_RestartRead                          !
!==============================================================!

  subroutine TSimulation_RestartRead( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation)   :: this

    ! Declare local variables
    integer :: i

    if( RootProc ) then

      ! Read contents from restart file
      read( iounit_restart, '(2I10)' ) Step, StepTotal
      read( iounit_restart, '(2L5)' ) Equilibration, NVTEquilibration

    end if

#if MPI_VER > 0
    call MPI_Bcast( Step, 1, MPI_INTEGER, &
&     NRootProc, MPI_COMM_WORLD, ierror )
    call MPI_Bcast( StepTotal, 1, MPI_INTEGER, &
&     NRootProc, MPI_COMM_WORLD, ierror )
    call MPI_Bcast( Equilibration, 1, MPI_LOGICAL, &
&     NRootProc, MPI_COMM_WORLD, ierror )
    call MPI_Bcast( NVTEquilibration, 1, MPI_LOGICAL, &
&     NRootProc, MPI_COMM_WORLD, ierror )
#endif

    ! Set current block number
    if( BlockSize > 0 ) then
      NBlocks = 1 + (Step - 1) / BlockSize
      NBlockSizes = int( sqrt( real( Step / BlockSize, RK ) ) )
    end if

    ! Read ensembles
    do i = 1, this%NEnsembles
      call RestartRead( this%Ensemble(i) )
    end do

    ! Close restart file
    call FileClose( iounit_restart )

 end subroutine TSimulation_RestartRead






end module ms2_simulation
