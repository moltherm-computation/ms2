!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 1.0                !
!  (c) 2011 by TU Kaiserslautern                               !
!      P.O. Box 67653                                          !
!      67653 Kaiserslautern                                    !
!==============================================================!
!  Module ms2_simulation                                       !
!  Contains TSimulation object                                 !
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
!DEC$ MESSAGE:'Compiling ms2_simulation.F90...'
#endif

module ms2_simulation

  use ms2_global
  use ms2_ensemble
  use ms2_stopwatch



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

#if  TRANS == 1
!TRANSPORT_start
    ! I/O unit for correlation function
    integer :: iounit_rescf
!TRANSPORT_END
#endif


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

  interface RunMCStep_Gibbs
    module procedure TSimulation_RunMCStep_Gibbs
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
  
  interface RDFOpen
    module procedure TSimulation_RDFOpen
  end interface

  interface RDFUpdate
    module procedure TSimulation_RDFUpdate
  end interface

  interface RDFClose
    module procedure TSimulation_RDFClose
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
      write( IOBuffer, '("Restarting from file: ", A)' ) trim( RestartFileName )
      call LogWrite

    else
      write( IOBuffer, '("Using parameters from file: ", A)' ) trim( ParameterFileName )
      call LogWrite
    end if

#else
    call FileReset( iounit_config, ProgramFileName//ConfigFileExtension )
    call FileReadParameter( str, iounit_config, IdRestart, .true., "NO" )
    select case( str )

    case( 'YES', 'Yes', 'yes' )
      Restart = .true.
      call FileReadParameter( RestartFileName, iounit_config, IdRestartFileName, .true. )
      write( IOBuffer, '("Restarting from file: ", A)' ) RestartFileName
      call LogWrite
      call FileReset( iounit_restart, RestartFileName )
      read( iounit_restart, '(A128)' ) ParameterFileName

    case( 'NO', 'No', 'no' )
      call FileReadParameter( ParameterFileName, iounit_config, IdParamsFileName, .true. )

    case default
      call Error( 'Select yes/no for restart in file '// ProgramFileName//ConfigFileExtension )
    end select
#endif

    write( IOBuffer, '(72(1H*))')
    call LogWrite
    write( IOBuffer, '(T24, "Reading Simulation Input")')
    call LogWrite
    write( IOBuffer, '(72(1H*))')
    call LogWrite
    write( IOBuffer, '("Parameter file name: ", A)' ) trim( ParameterFileName )
    call LogWrite

#if ARCH != 1 && ARCH != 2 && ARCH != 3
    call FileClose( iounit_config )
#endif
!    call LogWriteBlank

    ! Open parameter file for reading
    call FileReset( iounit_params, ParameterFileName )
    call LogWriteBlank
    write( IOBuffer, '(72(1H-))')
    call LogWrite
    write( IOBuffer, '(T20, "Reading parameters of simulation")' )
    call LogWrite

    ! Read name tag for output files
    call FileReadParameter( str, iounit_params , IdOutputNameTag, .true., status=stat )
    if ( OutputNameTagfromCommandline ) then

      if ( RootProc .and. stat .eq. 0 ) then
        print *,"INFO: output prefix from command line (", trim(OutputNameTag) &
&              ,") overwrites the one from the parameter file (", trim(str) ,")"
      end if

      str = "(from command line)"

    else if ( stat .eq. 0 ) then
      OutputNameTag = trim(str)                  ! possible truncation
      str = "(from parameter file)"

    else
      str = "(default)"
    end if
    write( IOBuffer, '("Name tag for output ",A,": ",T44, A)' ) trim( str ), trim( OutputNameTag )
    call LogWrite

    ! Read type of units
    call FileReadParameter( str, iounit_params , IdUseReducedUnits, .true. )
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
    write( IOBuffer, '("System of units: ",T26, A)' ) trim( str )
    call LogWrite

    ! Read unit of length
    call FileReadParameter( UnitLength, iounit_params, IdUnitLength, .true., 3.5_RK )
    UnitLength = UnitLength * Angstroem
    write( IOBuffer, '("Unit of length: ",T23, F8.3, " A")' ) UnitLength / Angstroem
    call LogWrite

    ! Read unit of energy
    call FileReadParameter( UnitEnergy, iounit_params, IdUnitEnergy, .true., 100.0_RK )
    UnitEnergy = UnitEnergy * kBoltzmann
    write( IOBuffer, '("Unit of energy: ",T23, F8.3, " K")' ) UnitEnergy / kBoltzmann
    call LogWrite

    ! Read unit of mass
    call FileReadParameter( UnitMass, iounit_params, IdUnitMass, .true., 40.0_RK )
    UnitMass = UnitMass * .001_RK / NAvogadro
    write( IOBuffer, '("Unit of mass:   ",T23, F8.3, " a.u.")' ) UnitMass * NAvogadro * 1000._RK
    call LogWrite
    call LogWriteBlank

    ! Calculate derived reduced units
    UnitVolume = UnitLength**3
    UnitTemperature = UnitEnergy / kBoltzmann
    UnitDensity = .001_RK / NAvogadro / UnitVolume
    UnitTime = sqrt(UnitMass / UnitEnergy) * UnitLength
    UnitForce = UnitEnergy / UnitLength
    UnitTorque = UnitEnergy
    UnitPressure = UnitForce / UnitLength**2
    UnitInertia = UnitMass * UnitLength**2
    UnitCharge = sqrt( 4._RK * Pi * VacuumPermittivity * UnitLength * UnitEnergy )
    UnitDipole = UnitCharge * UnitLength
    UnitQuadrupole = UnitCharge * UnitLength**2

    ! Read type of simulation
    call FileReadParameter( str, iounit_params , IdSimulationType, .true. )
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

    case( 'GibbsMC', 'gibbsmc', 'gibbs', 'Gibbs' )
      SimulationType = Gibbs
      SimulationTypeString = 'Gibbs-Monte-Carlo'

    case default
      call Error( trim( str )//' simulation is not implemented' )
    end select
    write( IOBuffer, '("Simulation type: ",T26, A)' ) trim( SimulationTypeString )
    call LogWrite

    ! Read parameters specific to given simulation type
    if( SimulationType .eq. SecondVirialCoeff ) then

      ! Read number of orientations
      call FileReadParameter( NOrient, iounit_params , IdNOrient, .true. )
      write( IOBuffer, '("Number of orientations: ",T24, I7)' ) NOrient
      call LogWrite

      ! Read number of steps
      call FileReadParameter( NSteps, iounit_params , IdRSteps, .true. )
      write( IOBuffer, '("Number of radial steps: ",T24, I7)' ) NSteps
      call LogWrite

      ! Read minimum radius
      call FileReadParameter( MinRadius, iounit_params , IdMinRadius, .true. )
      if( .not. UseReducedUnits ) then
        MinRadius = MinRadius / UnitLength * Angstroem
      end if
      write( IOBuffer, '("Minimum radius: ",T27, F8.3, " A")' ) MinRadius * UnitLength / Angstroem
      call LogWrite

      ! Read maximum radius
      call FileReadParameter( MaxRadius, iounit_params , IdMaxRadius, .true. )
      if( .not. UseReducedUnits ) then
        MaxRadius = MaxRadius / UnitLength * Angstroem
      end if
      write( IOBuffer, '("Maximum radius: ",T27, F8.3, " A")' ) MaxRadius * UnitLength / Angstroem
      call LogWrite

      ! Set output frequencies
      BlockSize = 0
      ErrorsUpdateFrequency = NSteps
      VisualUpdateFrequency = 0
      RDFUpdateFrequency = 0
      ! Set cutoff mode
      CutoffMode = CenterofMass

    else

      ! Read MD simulation parameters
      if( SimulationType .eq. MolecularDynamics ) then

        ! Type of integrator
        call FileReadParameter( str, iounit_params , IdIntegratorType, .true., "GEAR" )
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
        write( IOBuffer, '("Integrator type: ",T26, A)' ) trim( IntegratorTypeString )
        call LogWrite

        ! Time step
        call FileReadParameter( TimeStep, iounit_params , IdTimeStep, .true., 5.0E-4_RK )
        write( IOBuffer, '("Time step: ",T26, F9.6, " fs")' ) TimeStep * UnitTime * 1E15_RK
        call LogWrite
        write( IOBuffer, '("Reduced time step: ",T26, F9.6)' ) TimeStep
        call LogWrite
        TimeStep2 = .5_RK * TimeStep
        TimeStepSquared = TimeStep**2
        TimeStepSquared2 = .5_RK * TimeStepSquared
        TimeStepSquaredInv2 = .5_RK / TimeStepSquared

      ! Read MC simulation parameters
      else

        ! Acceptance rate
        call FileReadParameter( Acceptance, iounit_params , IdAcceptance, .true., 0.5_RK )
        if( Acceptance < 0.05_RK ) then
          Acceptance = 0.05_RK
        else if( Acceptance > 0.95_RK ) then
          Acceptance = 0.95_RK
        end if
        write( IOBuffer, '("Acceptance rate: ",T24, F6.2, "%")' ) Acceptance * 100._RK
        call LogWrite
        AccUpperLimit = Acceptance * 1.1_RK
        AccLowerLimit = Acceptance * 0.9_RK

      end if

      ! Read type of ensembles
      call FileReadParameter( str, iounit_params , IdEnsembleType, .true. )
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
      call LogWriteBlank
      write( IOBuffer, '("Ensemble type: ",T26, A)' ) trim( EnsembleTypeString )
      call LogWrite

      ! Check whether simulation type is applicable to ensemble type
      if( SimulationType .eq. MonteCarlo .and. .not. ConstantTemperature ) &
&         call Error( trim( SimulationTypeString )//" simulation of " &
&         //trim( EnsembleTypeString )//" ensemble is not implemented" )

      if( (EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA) &
&         .and. .not. SimulationType .eq. MonteCarlo ) &
&         call Error( trim( SimulationTypeString )//" simulation of " &
&         //trim( EnsembleTypeString )//" ensemble is not implemented" )

      ! Read number of MC overlap reduction steps
      call LogWriteBlank
      if( SimulationType .eq. MolecularDynamics ) then
        call FileReadParameter( NStepsMC, iounit_params , IdNStepsMC, .true., 0 )
        if( NStepsMC > 0 ) then
          write( IOBuffer, '("Number of MC overlap reduction steps: ",T40, I7)' ) NStepsMC
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
        call FileReadParameter( NStepsMC, iounit_params , IdNStepsMC, .true., 0 )
        GradInsInit = NStepsMC
        write( IOBuffer, '("Grad. Ins. initialization (if needed): ", T40, I7)' ) GradInsInit
        call LogWrite
      end if

      ! Read number of NVT equilibration steps
      call FileReadParameter( NStepsV, iounit_params , IdNStepsV, .true., 0 )
      write( IOBuffer, '("Number of NVT equilibration steps: ",T40, I7)' ) NStepsV
      call LogWrite

      ! Read number of NPT equilibration steps
      if( ConstantPressure ) then
        if( EnsembleType .eq. EnsembleTypeHA ) then
          call FileReadParameter( NStepsP, iounit_params , IdNStepsMueP, .true., 0 )
          write( IOBuffer, '("Number of HA equilibration steps: ",T40, I7)' ) NStepsP
          call LogWrite

        else
          call FileReadParameter( NStepsP, iounit_params , IdNStepsP, .true., 0 )
          write( IOBuffer, '("Number of NPT equilibration steps: ",T40, I7)' ) NStepsP
          call LogWrite
        end if

      else if( EnsembleType .eq. EnsembleTypeGE ) then
        call FileReadParameter( NStepsP, iounit_params , IdNStepsMue, .true., 0 )
        write( IOBuffer, '("Number of GE equilibration steps: ",T40, I7)' ) NStepsP
        call LogWrite

      else
        NStepsP = 0
      end if

      ! Read number of production steps
      call FileReadParameter( NSteps, iounit_params , IdNSteps, .true., 0 )
      write( IOBuffer, '("Number of production steps: ",T40, I7)' ) NSteps
      call LogWrite
      call LogWriteBlank

#if MPI_VER > 0
      if ( SimulationType .eq. MonteCarlo ) then
        NSteps = ceiling(real(NSteps)/NProcs)
      endif
#endif

      ! Read frequency of updating result file
      call FileReadParameter( BlockSize, iounit_params , IdBlockSize, .true., NSteps )
      if( BlockSize > 0 ) then
        write( IOBuffer, '("Result files will be updated each", I7, " time steps")' ) BlockSize
      else
        write( IOBuffer, '("All result files will not be created")' )
      end if
      call LogWrite

      ! Calculate number of blocks and block sizes
      if( BlockSize > 0 ) then
        NBlocksMax = ceiling(max( NStepsV, NStepsP, NSteps ) / real(BlockSize))

        ! Warning, if simulation is extended
        if ( mod(NSteps,BlockSize) .ne. 0._RK) then
          if (NSteps > BlockSize) then
            NSteps = ceiling(real(NSteps)/BlockSize)*BlockSize

            if ( SimulationType .eq. MonteCarlo ) then
              write( IOBuffer, '("Production steps are extended to",T40, I7, " cycles")' ) NSteps*NProcs
            else
              write( IOBuffer, '("Production steps are extended to",T40, I7, " time steps")' ) NSteps
            end if
            call LogWrite

          else
            BlockSize = NSteps
            write( IOBuffer, '("BlockSize is reduced to ",T40, I7, " due to small number of steps")' ) BlockSize
            call LogWrite
            NBlocksMax = ceiling(max( NStepsV, NStepsP, NSteps ) / real(BlockSize))
          endif
        end if

        NBlockSizesMax = int( sqrt( real( NSteps / BlockSize, RK ) ) )
        if (NBlocksMax .eq. 0) then
            call Error( 'ResultFreq < RunSteps, please change input variables' )
        end if

      else
        NBlocksMax = 0
        NBlockSizesMax = 0
      end if
      
      if ( NBlocksMax .lt. 10) then
        call LogWriteBlank
        write(IOBuffer, '("!!! WARNING !!!")')
        call LogWrite
        write(IOBuffer, '("Underestimated variances expected due to the small number of blocks.")' )
        call LogWrite
        call LogWriteBlank
      end if

      ! Read frequency of updating final result file
      if ( BlockSize > 0 ) then
        call FileReadParameter( ErrorsUpdateFrequency, iounit_params , IdErrorsUpdateFrequency, .true., 0 )
        if( ErrorsUpdateFrequency < 1 ) then
          ErrorsUpdateFrequency = NSteps
        else if( ErrorsUpdateFrequency < BlockSize * 4 ) then
          ErrorsUpdateFrequency = BlockSize * 4
        end if

        if( ErrorsUpdateFrequency < NSteps ) then
          write( IOBuffer, '("Final result files will be updated each", I7, " time steps")' ) ErrorsUpdateFrequency
        else
          write( IOBuffer, '("Final result files will be created at the end")' )
        end if

      call LogWrite
      else
        ErrorsUpdateFrequency = NSteps
      end if

      ! Read frequency of updating visualisation file
      call FileReadParameter( VisualUpdateFrequency, iounit_params , IdVisualUpdateFrequency, .true., 0 )
      if( VisualUpdateFrequency > 0 ) then
        write( IOBuffer, '("Visualization files will be updated each", I7, " time steps")' ) VisualUpdateFrequency
      else
        write( IOBuffer, '("Visualization files will not be created")' )
      end if
      call LogWrite
      call LogWriteBlank
      
      ! Read frequency of updating visualisation file
      call FileReadParameter( RDFUpdateFrequency, iounit_params , IdRDFUpdateFrequency, .true., 0 )
      if( RDFUpdateFrequency > 0 ) then
        write( IOBuffer, '("RDF files will be updated each", I7, " time steps")' ) RDFUpdateFrequency
      else
        write( IOBuffer, '("RDF files will not be created")' )
      end if
      call LogWrite
      call LogWriteBlank
      
      if( RDFUpdateFrequency > 0 ) then
      call FileReadParameter( RDFNumberShells, iounit_params , IdRDFNumberShells, .true., 200 )
        write( IOBuffer, '("RDF will operate with", I7, " shells")' ) RDFNumberShells
      call LogWrite
      call LogWriteBlank
      end if

      ! Read cutoff mode
      call FileReadParameter( str, iounit_params , IdCutoffMode, .true., "COM" )
      select case( str )

      case( 'COM', 'com', 'CenterOfMass', 'CenterofMass', 'centerofmass' )
        CutoffMode = CenterofMass
        CutoffModeString = 'Center of Mass'
        write( IOBuffer, '("Cutoff mode: ",T26, A)' ) trim( CutoffModeString )

      case( 'Site', 'site', 'Site-Site', 'site-site' )
        CutoffMode = SiteSite
        CutoffModeString = 'Site-Site'
        write( IOBuffer, '("Cutoff mode: ",T26, A)' ) trim( CutoffModeString )

      case default
        call Error( trim( str )//' is not a valid cutoff mode' )
      end select
      call LogWrite
      call LogWriteBlank

      ! Read LongRange mode
      call FileReadParameter( str, iounit_params , IdLongRange, .true., "rf" )
      select case( str )
        case( 'Ewald', 'ew', 'ewald', 'EWALD')
            LongRange = Ewald
            LongRangeString = 'EwaldSum'
            write( IOBuffer, '("Long Range Correction: ", A)' ) trim( LongRangeString )
            call LogWrite
            call FileReadParameter( KappaL_h, iounit_params , IdKappa, .true., 5.6_RK )
            write( IOBuffer, '("Ewald: KappaL:", T23, F8.3)' ) KappaL_h
            call LogWrite

            call FileReadParameter( nsqmax_h, iounit_params , Idnsqmax, .true., 27 )
            write( IOBuffer, '("Ewald: NsqMax:",T20, I7)' ) nsqmax_h
            call LogWrite

            call FileReadParameter( nvecmax_h, iounit_params , IdNVecMax, .true., 1000 )
            write( IOBuffer, '("Ewald: NVecMax:",T20, I7)' ) nvecmax_h
            call LogWrite

            call FileReadParameter( nmax_h, iounit_params , IdNMax, .true., 5 )
            write( IOBuffer, '("Ewald: NMax:",T20, I7)' ) nmax_h
            call LogWrite

        case( 'PME', 'pme', 'SPME', 'spme')
            LongRange = PME
            LongRangeString = 'Smooth Particle Mesh Ewald Summation'
            write( IOBuffer, '("Long Range Correction: ", A)' ) trim( LongRangeString )
            call LogWrite
            ! Read SPM Ewald Parameters
            call FileReadParameter( KappaL_h, iounit_params , IdKappa, .true., 5.6_RK )
            write( IOBuffer, '("Ewald: KappaL:", F8.3)' )KappaL_h
            call LogWrite

            call FileReadParameter( grid_h, iounit_params , IdGrid, .true. )
            write( IOBuffer, '("Grid Space SPME:", I7)' ) grid_h
            call LogWrite

            call FileReadParameter( spline_h, iounit_params , IdSpline, .true. )
            write( IOBuffer, '("order of SPME Spline:", I7)' ) spline_h

        case( 'ReactionField', 'RF', 'reactionfield', 'rf' )
            LongRange = RField
            LongRangeString = 'Reaction Field'
            write( IOBuffer, '("Long Range Correction: ", A)' ) trim( LongRangeString )

        case( 'ExtReactionField', 'ExtRF', 'extreactionfield', 'extrf' )
            LongRange = ExtRField
            LongRangeString = 'Extended Reaction Field by Tironi et al.'
            write( IOBuffer, '("Long Range Correction: ", A)' ) trim( LongRangeString )
            call LogWrite
            ! Read extended Reaction Field Parameters
            call FileReadParameter( debyelen_h, iounit_params , IdDebyeLen, .true.)
            write( IOBuffer, '("Debye Length [A]:", F8.3)' )debyelen_h

        case( 'Rodgers', 'rodgers' )
            LongRange = rodgers
            LongRangeString = 'Rodgers'
            write( IOBuffer, '("Long Range Correction: ", A)' ) trim( LongRangeString )
            call LogWrite
            ! Read Rodgers Parameters
            call FileReadParameter( KappaL_h, iounit_params , IdKappa, .true., 0.15_RK )
            write( IOBuffer, '("Rodgers Parameter KappaL:", F8.3)' )KappaL_h
            call LogWrite

      case default
        call Error( trim( str )//' is not a valid longrange correction' )
      end select
      call LogWrite

    end if

    ! Read type of simulation with/without internal degree of freedom
    ! Read type of MD simulation with/without internal degree of freedom  
    call FileReadParameter( str, iounit_params , IdUseIntDegFreed, .true., "off" )
    select case( str )
    case( 'ON', 'On', 'on', 'YES', 'yes' )
       UseIntDegFreed = .true.
       str = 'Flexible molecules'
    case( 'OFF', 'off', 'no', 'No' )
       UseIntDegFreed = .false.
       str = 'Rigid molecules'
    case default
       call Error( trim( str )//'To switch on internal degree of freedom use on or yes' )
    end select
    write( IOBuffer, '("Using internal degree of freedom: ", A)' ) trim( str )
    call LogWrite

    ! Read printIDF parameter - to print all contributions to inramolecular energy if need
    if (UseIntDegFreed) then
      call FileReadParameter( str, iounit_params , IdPrintIDF, .true., "off" )
      select case( str )
      case( 'ON', 'On', 'on', 'YES', 'yes' )
         printIDF = .true.
      case( 'OFF', 'off', 'no', 'No' )
         printIDF = .false.
      case default
         call Error( trim( str )//'To print contributions to intramolecular energy use on or yes' )
      end select
    ! Read tolerance for Shake/QShake algorithm, if < 0, then no constraint dynamics is used and all bond lengths can vibrate
      call FileReadParameter( Shake, iounit_params , IdShake, .true., 0.0_RK )
      if ( Shake > 0 ) then 
        str = 'yes'
      else 
        str = 'no, all bonds can vibrate' 
      end if
      write( IOBuffer, '("Using Shake algorithm for bonds: ", A)' ) trim( str )
      call LogWrite
      if (str == 'yes') then 
        write( IOBuffer, '("Shake tolerance: ", F9.6)' ) Shake
        call LogWrite
      end if   
    ! Read parameters for intramolecular nonbonded interactions
      call FileReadParameter( str, iounit_params , IdIntraLJEl, .true., "no" )
      select case( str )
      case( 'ON', 'On', 'on', 'YES', 'yes' ) ! include all intramolecular 1-5 electrostatic & LJ interaction 
         IntraLJEl = .true.
         str = 'Include all intramolecular 1-5 nonbonded interactions'
      case( 'OFF', 'off', 'no', 'No' )
         IntraLJEl = .false.
         str = 'No intramolecular nonbonded interactions'
      case default
         call Error( trim( str )//'To switch on intramolecular 1-5 nonbonded interactions use on or yes' )
      end select
      write( IOBuffer, '("Intramolecular nonbonded interactions: ", A)' ) trim( str )
      call LogWrite

      if (IntraLJEl) then 
        call FileReadParameter( str, iounit_params , IdLJEl14, .true., "no" )
        select case( str )
        case( 'ON', 'On', 'on', 'YES', 'yes' ) ! include all intramolecular 1-4 electrostatic & LJ interaction 
           LJEl14 = .true.
           str = 'Include all intramolecular 1-4 nonbonded interactions'
        case( 'OFF', 'off', 'no', 'No' )
           LJEl14 = .false.
           str = 'No intramolecular 1-4 nonbonded interactions'
        case default
           call Error( trim( str )//'To switch on intramolecular 1-4 nonbonded interactions use on or yes' )
        end select
        write( IOBuffer, '("Intramolecular nonbonded interactions: ", A)' ) trim( str )
        call LogWrite
      end if
    end if
    
    ! Read number of ensembles
    call FileReadParameter( this%NEnsembles, iounit_params , IdNEnsembles, .true., 1 )
    write( IOBuffer, '("Number of ensembles:",T24, I3)' ) this%NEnsembles
    call LogWrite

    
#if  TRANS == 1
!TRANSPORT_start
    ! Read correlation function mode
    call FileReadParameter( str , iounit_params , IdCorrFun, .true. , 'no' )
    select case( str )

    case( 'yes' , 'ok', 'ja' )
      CorrfunMode = active
      CorrfunModeString = 'Include transport properties'

    case( 'no', 'nein' )
      CorrfunMode = inactive
      CorrfunModeString = 'No transport properties'
      call Error( 'Use a binary compiled without -DTRANS if you do not &
&                  wish to calculate transport properties. If you do, set CorrFunMode = yes ' )

    case default
      call Error( 'Unknown transport properties ('//trim(IdCorrFun)//'='//trim(str)//')' )
    end select
    write( IOBuffer, '("Transport properties:",T26, A)' ) trim(CorrfunModeString)
    call LogWrite
!TRANSPORT_END
#endif

    if( (SimulationType .eq. Gibbs .and. this%NEnsembles .ne. 2) )  &
&       call Error( trim( SimulationTypeString )//" simulation of " &
&       //trim( SimulationTypeString )//" needs 2 Ensembles" )

    ! Create ensembles
    allocate( this%Ensemble(this%NEnsembles), STAT = stat )
    call AllocationError( stat, 'ensembles', this%NEnsembles )
    do i = 1, this%NEnsembles
      if (LongRange .eq. Ewald) then
            this%Ensemble(i)%KappaL = KappaL_h
            this%Ensemble(i)%nsqmax = nsqmax_h
            this%Ensemble(i)%nvecmax = nvecmax_h
            this%Ensemble(i)%nmax = nmax_h

#ifdef SPME
      else if (LongRange .eq. PME) then
            this%Ensemble(i)%KappaL = KappaL_h
            this%Ensemble(i)%gridx  = grid_h
            this%Ensemble(i)%gridy  = grid_h
            this%Ensemble(i)%gridz  = grid_h
            this%Ensemble(i)%splineorder = spline_h
            allocate(this%Ensemble(i)%qgrida(2,(grid_h)**3+1),STAT=stat)
            if(stat >0) write(*,*) 'Allocation Error grida'
            allocate(this%Ensemble(i)%qgrida_old(2,(grid_h)**3+1),STAT=stat)
            if(stat >0) write(*,*) 'Allocation Error grida_old'
            allocate(this%Ensemble(i)%qgridb(2,(grid_h)**3+1),STAT=stat)
            if(stat >0) write(*,*) 'Allocation Error gridb'
#endif

      else if (LongRange .eq. ExtRField) then
            this%Ensemble(i)%DebyeLen = debyelen_h / Angstroem * UnitLength
      else if (LongRange .eq. Rodgers) then
            this%Ensemble(i)%KappaL = KappaL_h
      end if
      if( SimulationType .eq. SecondVirialCoeff ) then
        call ConstructSVC( this%Ensemble(i), i )
      else
        call Construct( this%Ensemble(i), i )
      end if
    end do

  GradInsFrequency = BlockSize
  NFullFluct = 20
  maxcounter = 0

    ! Close parameter file
    call FileClose( iounit_params )
    write( IOBuffer, '(T18, "Reading Simulation Input successful")')
    call LogWrite
    write( IOBuffer, '(72(1H*))')
    call LogWrite
    call LogWriteBlank

    ! Create accumulators
    call CreateAccumulators( this )

    ! Set I/O unit numbers
    this%iounit_result = iounit_result
    this%iounit_runave = iounit_runave
    this%iounit_errors = iounit_errors
#if  TRANS == 1
    this%iounit_rescf  = iounit_rescf  !TRANSPORT_thisline
#endif

    ! Open result and visualisation files
    call LogWriteBlank
    call LogWriteBlank
    write( IOBuffer, '(72(1H*))')
    call LogWrite
    write( IOBuffer, '(T28, "Start Simulation")')
    call LogWrite
    write( IOBuffer, '(72(1H*))')
    call LogWrite
    call ResultOpen( this )
    call VisualOpen( this )
    call RDFOpen( this )

  end subroutine TSimulation_Construct



!==============================================================!
!  Subroutine TSimulation_Destruct                             !
!==============================================================!

  subroutine TSimulation_Destruct( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Close result and visualisation files
    call LogWriteBlank
    call ResultClose( this )
    call VisualClose( this )
    call RDFClose( this )
    
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
    integer :: i, j, s, t
    logical :: NPartsOk
    type(TStopwatch) :: RunTimer,RunStepsTimer
    integer :: k

#if MPI_VER > 0
    type(TComponent), pointer :: pc
    type(TInteraction), pointer :: pi
    integer :: n1, n2
    
    integer :: color, NGroups, Proc_Max_Eff
    integer :: statusHost, lengthHost, tmpVal
    character(255) :: hostnameStr
    logical :: multNodes

#endif 

    tooManyParticles = .false.
    call Construct(RunTimer,"TSimulation_Run",CStopwatch_doMPIStartBarrier)
    call Construct(RunStepsTimer)

    call start_Timer(RunTimer)
    call logwritestart_Timer(RunTimer)
    
#if MPI_VER > 0
    ! This is for the restart - in case there is a restart, the root reads and communicates
    if (SimulationType .eq. MonteCarlo) then 
        RootProc = NProc==NRootProc
    endif
#endif

    if( Restart ) then
      call RestartRead( this )
      StepStart = Step + 1
      MCOverlapReduction = .false.    ! no MC overlap reduction in case of restart
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

#if MPI_VER > 0 
    ! For MC parallelization: if we have common equilibration 
    ! active, we revert to one rootproc
    if (SimulationType .eq. MonteCarlo) then 
      if (CommonEqui) then
        NProcs_W = NProcs
        NProc_W = NProc

        multNodes = .false.

        call MPI_GET_PROCESSOR_NAME(hostnameStr,lengthHost, ierror)
        if (len(trim(hostnameStr))==0) then
          statusHost = 1    
        else
          statusHost = 0         
        endif    
        
        if (statusHost >0) then
          
          write( IOBuffer, '("WARNING: This platform/compiler does not support MPI_GET_PROCESSOR_NAME")' ) 
          call LogWrite
          write( IOBuffer, '("WARNING: properly, therefore equilibration is split arbitrarily.")' ) 
          call LogWrite
          write( IOBuffer, '("WARNING: This may result in poor performance during equilibration")' ) 
          call LogWrite
          
          !The maximum number of processes
          Proc_Max_Eff = 8
            
          if (NProcs .gt. Proc_Max_Eff) then
            multNodes = .true.
            NGroups = NProcs/Proc_Max_Eff          
            color=mod(NProc,NGroups)  
        
            if (NProc .ge. NGroups*Proc_Max_Eff) then
              color = 1000000              
            endif
          
            call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,NProc,Communicator,ierror) 
            call SetCommunicator( Communicator )             
          endif    
        else
          if (NProcs .gt. 4) then
            color = 0
            do i=1,lengthHost
              tmpVal = ichar (trim(hostnameStr(i:i)))
              color = color + (tmpVal**2)*i
            enddo

            call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,NProc,Communicator,ierror) 
            call SetCommunicator( Communicator )
       
              
            call MPI_ALLREDUCE(NProcs, Proc_Max_Eff, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror )
          
            if (Proc_Max_Eff .lt. 3) then
               write( IOBuffer, '("WARNING: MPI_GET_PROCESSOR_NAME may have given a processor specific name")' ) 
               call LogWrite
               write( IOBuffer, '("WARNING: if you have more than 2 PE per node, something is wrong. Try a different")' ) 
               call LogWrite
               write( IOBuffer, '("WARNING: compiler. P.s: Due to that, the equilibration is slow, sorry!")' ) 
               call LogWrite
            endif
              
              
            if (Proc_Max_Eff .lt. NProcs_W) then
              multNodes = .true.
              if (Proc_Max_Eff .gt. NProcs) then
                color = 1000000
                tmpVal = 0
              else
                tmpVal = 1
              endif
              call MPI_ALLREDUCE(tmpVal, NGroups, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror )
              NGroups = NGroups/Proc_Max_Eff
            endif
          endif
        endif
        
        if (multNodes) then
          
           
           if (color == 1000000) then
              if (MCOverlapReduction) then
                NStepsMC = 1
              endif

              if (NVTEquilibration) then
                NStepsV = 1
              endif

              if (Equilibration) then
                NStepsP = 1
              endif
           endif
           if (RootProc) then
             if (NProc_W .ne. NRootProc) then
               write( IOBuffer, '(I16)' ) NProc_W  
               call FileRewrite( iounit_log, trim( OutputNameTag )//'_Equi_'//trim( adjustl( IOBuffer ) )//LogFileExtension )
               
               do j = 1, this%NEnsembles

                 ! Open running average result file
                 write( IOBuffer, '(I16)' ) NProc_W
                 call FileRewrite( this%Ensemble(j)%iounit_runave, &
&                     trim( OutputNameTag )//'_Equi_'//trim( adjustl( IOBuffer ) )//RunAveFileExtension )

                 ! Open result file
                 write( IOBuffer, '(I16)' ) NProc_W
                 call FileRewrite( this%Ensemble(j)%iounit_result, &
&                     trim( OutputNameTag )//'_Equi_'//trim( adjustl( IOBuffer ) )//ResultFileExtension )
               enddo
             endif 
           endif   
           
           call Randomize( seed = (5333+(color+1)) )
 
        endif

      else
        call Randomize( seed = (5333*(NProc+1)) )
      endif

      ! adapt procrange for to the given equilibration scheme
      do j = 1, this%NEnsembles
        do i = 1, this%Ensemble(j)%NComponents
           pc => this%Ensemble(j)%Component(i)
           pc%NPart1 = ProcRange( pc%NPart, pc%NPart0, pc%NPart2 )
        end do
        
        ! Recalculate Energies to avoid energy artefacts 
        call Mol2Atom( this%Ensemble(j) )
        ! Recalculate LongRange Correction
        call CalculateCorr( this%Ensemble(j) )
        if ( (LongRange .eq. Ewald) .or. (LongRange .eq. PME) ) then
          this%Ensemble(j)%NBox1 = ProcRange( this%Ensemble(j)%BoxenAnzahlMax, this%Ensemble(j)%NBox0, this%Ensemble(j)%NBox2 )
        end if

         ! Set all potential energy matrices
         call Energy( this%Ensemble(j), this%Ensemble(j)%EPot )
         call UpdateEnergy( this%Ensemble(j) )

      end do
    endif 
#endif

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
      SimulationType = MonteCarlo
      
      if ( UseIntDegFreed ) then
        do i=1,this%NEnsembles
          call Flex2Rigid( this%Ensemble(i) )
        end do
      end if

      call Timer_setTag(RunStepsTimer,"MC overlap reduction")
      call start_Timer(RunStepsTimer)
      call logwritestart_Timer(RunStepsTimer)

      call RunSteps( this, StepStart, StepEnd )
      
      if ( UseIntDegFreed ) then
        do i=1,this%NEnsembles
          call Rigid2Flex( this%Ensemble(i) )
        end do
      end if
      
      call stop_Timer(RunStepsTimer)
      call logwritestop_Timer(RunStepsTimer)

      if( .not. TerminateProgram ) then
        write( IOBuffer, '("MC overlap reduction completed")' )
        MCOverlapReduction = .false.
        SimulationType = MolecularDynamics

        do i = 1, this%NEnsembles
          call InitMolecularDynamics( this%Ensemble(i), .true. )
        end do
      else
        write( IOBuffer, '("MC overlap reduction TERMINATED")' )
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

        call Timer_setTag(RunStepsTimer,"NVT equilibration")
        call start_Timer(RunStepsTimer)
        call logwritestart_Timer(RunStepsTimer)
        call RunSteps( this, StepStart, StepEnd )
        call stop_Timer(RunStepsTimer)
        call logwritestop_Timer(RunStepsTimer)

        if( .not. TerminateProgram ) then
          write( IOBuffer, '("NVT equilibration completed")' )
          NVTEquilibration = .false.
        else
          write( IOBuffer, '("NVT equilibration TERMINATED")' )
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

          call Timer_setTag(RunStepsTimer,"GE equilibration")
          call start_Timer(RunStepsTimer)
          call logwritestart_Timer(RunStepsTimer)
          call RunSteps( this, StepStart, StepEnd )
          call stop_Timer(RunStepsTimer)
          call logwritestop_Timer(RunStepsTimer)

          if( .not. TerminateProgram ) then
            call CheckNPart( this, NPartsOk )

            if( NPartsOk ) then
              write( IOBuffer, '("GE equilibration completed")' )
              Equilibration = .false.

            else
              write( IOBuffer, '("GE equilibration ended with too many/too less particles")' )
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
            write( IOBuffer, '("GE equilibration TERMINATED")' )
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

          call Timer_setTag(RunStepsTimer,"HA equilibration")
          call start_Timer(RunStepsTimer)
          call logwritestart_Timer(RunStepsTimer)
          call RunSteps( this, StepStart, StepEnd )
          call stop_Timer(RunStepsTimer)
          call logwritestop_Timer(RunStepsTimer)

          if( .not. TerminateProgram ) then
            call CheckNPart( this, NPartsOk )
            if( NPartsOk ) then
              write( IOBuffer, '("HA equilibration completed")' )
              Equilibration = .false.
            else
              write( IOBuffer, '("HA equilibration ended with too many/too less particles")' )
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
            write( IOBuffer, '("HA equilibration TERMINATED")' )
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

          call Timer_setTag(RunStepsTimer,"NPT equilibration")
          call start_Timer(RunStepsTimer)
          call logwritestart_Timer(RunStepsTimer)
          call RunSteps( this, StepStart, StepEnd )
          call stop_Timer(RunStepsTimer)
          call logwritestop_Timer(RunStepsTimer)

          if( .not. TerminateProgram ) then
            write( IOBuffer, '("NPT equilibration completed")' )
            Equilibration = .false.
          else
            write( IOBuffer, '("NPT equilibration TERMINATED")' )
          end if
          call LogWriteTime

        else if( SimulationType .eq. Gibbs ) then
          StepEnd = NStepsV
          call LogWriteBlank

          if( Restart ) then
            write( IOBuffer, '("Resuming Gibbs equilibration")' )
            Restart = .false.
          else
            write( IOBuffer, '("Starting Gibbs equilibration")' )
          end if

          call LogWriteTime
          call Timer_setTag(RunStepsTimer,"NPT equilibration")
          call start_Timer(RunStepsTimer)
          call logwritestart_Timer(RunStepsTimer)
          call RunSteps( this, StepStart, StepEnd )
          call stop_Timer(RunStepsTimer)
          call logwritestop_Timer(RunStepsTimer)

          if( .not. TerminateProgram ) then
            call CheckNPart( this, NPartsOk )

            if( NPartsOk ) then
              write( IOBuffer, '("Gibbs equilibration completed")' )
              Equilibration = .false.

            else
              write( IOBuffer, '("Gibbs equilibration ended with too many/too less particles")' )
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
            write( IOBuffer, '("Gibbs equilibration TERMINATED")' )
          end if
          call LogWriteTime


        else
          Equilibration = .false.
        end if
        StepStart = 1
      end if

      exit eqloop
    end do eqloop

    ! In the MC parallelization, every process is regarded as its own root from here 
    ! (the equilibration is finished. From now on, every process runs its own simulation etc.)
#if MPI_VER > 0 
    if (SimulationType .eq. MonteCarlo .and. CommonEqui) then 
      
      do k = 1, this%NEnsembles
          do i = 1, this%Ensemble(k)%NRealComponents
            do j = 1, this%Ensemble(k)%NRealComponents
              pi => this%Ensemble(k)%Interaction(j, i)
              n1 = pi%NPart1
              n2 = pi%NPart2
        
              call MPI_Allreduce( pi%EPot(1:n1, 1:n2), pi%EPotNew(1:n1, 1:n2), n1*n2 , &
&                  MPI_RK, MPI_SUM, Communicator, ierror )
              pi%EPot(1:n1, 1:n2) =  pi%EPotNew(1:n1, 1:n2)
       
              if ( this%Ensemble(k)%OptPressure ) then
                call MPI_Allreduce( pi%Virial(1:n1, 1:n2) ,pi%VirialNew(1:n1, 1:n2), n1*n2 , &
&                    MPI_RK, MPI_SUM, Communicator, ierror )
                pi%Virial(1:n1, 1:n2)  =  pi%VirialNew(1:n1, 1:n2)
              endif
            end do
          end do
      end do
      
      
      if (NProcs_W .gt. Proc_Max_Eff) then
      
        if (RootProc) then
          if (NProc_W .ne. NRootProc) then
            ! Close all files keeping track of the equilibration
            do j = 1, this%NEnsembles
              call FileClose( this%Ensemble(j)%iounit_runave )
              call FileClose( this%Ensemble(j)%iounit_result )
            enddo
            call LogClose
          endif
        endif
         
          if (NProcs_W .gt. NGroups*Proc_Max_Eff) then
           ! build new communicator, including the Root and all processes not having
           ! equilibrated
           if (NProc_W == NRootProc) then
            ! the Root receives the corresponding color (see above)
            color = 1000000
           endif 
          
           call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,NProc_W,Communicator,ierror) 
           call SetCommunicator( Communicator )
           
           ! only these processes are involved in the communication
           if  ((NProc_W .ge. NGroups*Proc_Max_Eff) .or. (NProc_W == NRootProc)) then
             do j = 1, this%NEnsembles
               call MPI_Bcast( this%Ensemble(j)%EPot, 1, MPI_RK, NRootProc, Communicator, ierror )
               call MPI_Bcast( this%Ensemble(j)%DispVol, 1, MPI_RK, NRootProc, Communicator, ierror )            
               do i = 1, this%Ensemble(j)%NComponents
                 call MPI_Bcast( this%Ensemble(j)%Component(i)%Pm0(:, :), size( this%Ensemble(j)%Component(i)%P0 ), &
&                     MPI_RK, NRootProc, Communicator, ierror )

                 if( this%Ensemble(j)%Component(i)%Molecule%isElongated ) then
                    call MPI_Bcast( this%Ensemble(j)%Component(i)%Qm0(:, :), size( this%Ensemble(j)%Component(i)%Q0 ), &
&                        MPI_RK, NRootProc, Communicator, ierror )
                 endif 
               enddo

               do i = 1,  this%Ensemble(j)%NRealComponents
                 call MPI_Bcast( this%Ensemble(j)%Component(i)%DispTran, 1, MPI_RK, NRootProc, Communicator, ierror )
                 call MPI_Bcast( this%Ensemble(j)%Component(i)%DispRot, 1, MPI_RK, NRootProc, Communicator, ierror )
               enddo
             enddo
           endif
          endif
          
          ! Set Communicator to COMM_WORLD
          call SetCommunicator (MPI_COMM_WORLD)
              
      endif      
      
      ! New random number seed for different simulations (distinct simulation in every process)
      call Randomize( seed = (5333*(NProc+1)) )

      ! adapt procrange such that each simulation calculates all its interactions from now on
      do j = 1, this%NEnsembles
        do i = 1, this%Ensemble(j)%NComponents
          pc => this%Ensemble(j)%Component(i)
          pc%NPart1 = ProcRange( pc%NPart, pc%NPart0, pc%NPart2 )
        end do
        ! Convert molecular coordinates to atom positions
        call Mol2Atom( this%Ensemble(j) )
        
        ! Recalculate LongRange Correction
        call CalculateCorr( this%Ensemble(j) )
        if ( (LongRange .eq. Ewald) .or. (LongRange .eq. PME) ) then
          this%Ensemble(j)%NBox1 = ProcRange( this%Ensemble(j)%BoxenAnzahlMax, this%Ensemble(j)%NBox0, this%Ensemble(j)%NBox2 )
        end if
        
        ! Set all potential energy matrices
        call Energy( this%Ensemble(j), this%Ensemble(j)%EPot )
        call UpdateEnergy( this%Ensemble(j) )
      end do


    endif 
#endif

     GradInsInitialization = .false.
     do j = 1, this%NEnsembles
       do i = 1, this%Ensemble(j)%NComponents
           if( (this%Ensemble(j)%Component(i)%WFMethod .eq. WFMethodGuess) .and. &
&              (this%Ensemble(j)%Component(i)%ChemPotMethod .eq. ChemPotMethodGradIns) ) then
             GradInsInitialization = .true.
           endif
       enddo
     enddo
     
      if( GradInsInitialization) then
       call LogWriteBlank

       if( Restart ) then
         write( IOBuffer, '("Resuming GradIns initialization")' )
         Restart = .false.
         StepStart = Step + 1

       else
         StepStart = 1
         write( IOBuffer, '("Starting GradIns initialization")' )
         call LogWrite
         write( IOBuffer, '("  (adjustment of weighting factors)")' )
       end if
       call LogWriteTime
       
       do Step = StepStart, GradInsInit
         do i = 1, this%NEnsembles
           call ChemicalPotential( this%Ensemble(i) )
         end do
       end do
       
       
       write( IOBuffer, '("Number of GradIns initialization iterations: ",T40, I7)' ) max(NStepsMC,1)*this%NEnsembles
       call LogWrite
       
       Step = 1
       if( .not. TerminateProgram ) then
         write( IOBuffer, '("GradIns initialization completed")' )
         GradInsInitialization = .false.
       else
         write( IOBuffer, '("GradIns initialization TERMINATED")' )
       end if

       call LogWriteTime
       StepStart = 1
     end if

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

      call Timer_setTag(RunStepsTimer,"simulation")
      call start_Timer(RunStepsTimer)
      call logwritestart_Timer(RunStepsTimer)
      call RunSteps( this, StepStart, StepEnd )
      call stop_Timer(RunStepsTimer)
      call logwritestop_Timer(RunStepsTimer)

      if( .not. TerminateProgram ) then
        write( IOBuffer, '("Simulation completed")' )
      else
        write( IOBuffer, '("Simulation TERMINATED")' )
      end if
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

    call stop_Timer(RunTimer)
    call logwritestop_Timer(RunTimer)

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
#if MPI_VER > 0 && ( ARCH == 1 || ARCH == 2 )
    logical :: AnyTerminateProgram
#endif

#if TRANS==1
    integer:: StepCF
#endif

    integer:: o, i, j, t, s

    ! Run simulation steps
    do Step = StepStart, StepEnd

      ! Update total number of steps
      StepTotal = StepTotal + 1

      ! Set current block number
      if( BlockSize > 0 ) then
        NBlocks = 1 + (Step - 1) / BlockSize
        NBlockSizes = int( sqrt( real( Step / BlockSize, RK ) ) )

#if TRANS==1
      ! Run simulation step
        if(( CorrfunMode == active ).and.(.not. Equilibration )) then
          do i = 1, this%Nensembles
            if(mod((Step+NStepCorr-1),NStepCorr) .eq. 0) then
              StepCF = (Step + NStepCorr -1) / NStepCorr
              if ( StepCF >= this%Ensemble(i)%Ncorr )then
                NBlocksCF = 1 + ( StepCF - 1 - this%Ensemble(i)%Ncorr ) / ( BlockSizeCF * this%Ensemble(i)%NSpancf )
                NBlockSizesCF = int( sqrt( real(( StepCF - this%Ensemble(i)%Ncorr) / (BlockSizeCF * this%Ensemble(i)%NSpancf ), RK)))
              else
                NBlocksCF     = 0
                NBlockSizesCF = 0
              end if
          
            end if
          end do
        end if
#endif
      end if

      ! Run simulation step
      select case( SimulationType )
      case( MolecularDynamics )
        call RunMDStep( this )
      case( MonteCarlo )
        call RunMCStep( this )
      case( SecondVirialCoeff )
        call RunSVCStep( this )
      case( Gibbs )
        call RunMCStep_Gibbs( this )
      end select

      ! Update result and visualisation files
      call ResultUpdate( this )
      call VisualUpdate( this )
      call RDFUpdate ( this )

      ! Update log and result files
      if( mod( Step, LogUpdateFrequency ) == 0 .or. Step == StepEnd ) call LogWriteStep
      if( .not. Equilibration .and. ( mod( Step, ErrorsUpdateFrequency ) == 0 .or. Step == StepEnd )) call ErrorsUpdate( this )

      ! Check for termination request (caused by signal handler)
#if MPI_VER > 0 && ( ARCH == 1 || ARCH == 2 )
      if (SimulationType .eq. MonteCarlo) then
        if (Step == StepEnd .and. .not. Equilibration) then
          call MPI_Allreduce( TerminateProgram, AnyTerminateProgram, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror )
          if( AnyTerminateProgram ) then
            TerminateProgram = .true.
            exit
          end if
        else
          call MPI_Allreduce( TerminateProgram, AnyTerminateProgram, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror )
          if( AnyTerminateProgram ) then
            TerminateProgram = .true.
            exit
          end if
        endif 

      else
        call MPI_Allreduce( TerminateProgram, AnyTerminateProgram, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror )
        if( AnyTerminateProgram ) then
          TerminateProgram = .true.
          exit
        end if
      endif
      
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
!  Subroutine TSimulation_RunMCStep_GibbsEnsemble              !
!==============================================================!

  subroutine TSimulation_RunMCStep_Gibbs( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i
    real(RK):: EPotOldliq
    real(RK):: VolOldliq
    real(RK):: EPotDelta
    real(RK):: dv, NEns_h
    logical :: accept

    integer :: NEns
    integer :: nc,np
    integer :: NTransfer


    NTransfer = 100
    ! Simulations Setup Check for Gibbs
    if ( (SimulationType .eq. Gibbs) .and. ConstantPressure ) then
      call Error( 'Gibbs Ensemble only implemented for NVT')
    end if
    if ( (SimulationType .eq. Gibbs) .and. (this%NEnsembles .ne. 2) ) then
      call Error( 'Gibbs Ensemble needs two SimBoxes: one liquid and one vapor SimBox')
    end if

    ! Run MC simulation step
    do i = 1, this%NEnsembles
      call RunMCStep( this%Ensemble(i) )
    end do

! Volume Change in both boxes
    if (.not. NVTEquilibration) then
      accept = .false.
      EPotOldliq = this%Ensemble(1)%EPot
      VolOldliq = this%Ensemble(1)%Volume0

      call Resize_Gibbs( this%Ensemble(1),dv,EPotDelta )
      call Resize_Gibbs( this%Ensemble(2),dv,EPotDelta,accept )

    ! Accept volume change
      call Update_Gibbs ( this%Ensemble(1),accept,EPotOldliq,VolOldliq )

! Particle change in both boxes
      DO i=1,NTransfer
        accept = .false.
        NEns_h = rnd(0._RK,1._RK)
        NEns   = int(anint(NEns_h) + 1._RK)

        call Remove_Gibbs( this%Ensemble(NEns),nc,np,EPotDelta )
        call Insert_Gibbs( this%Ensemble(3-NEns),nc,EPotDelta, accept )
        call Update_Gibbs ( this%Ensemble(NEns),nc,np,NTransfer,accept )
      END DO
    end if


  end subroutine TSimulation_RunMCStep_Gibbs



!==============================================================!
!  Subroutine TSimulation_CheckNPart                           !
!==============================================================!

  subroutine TSimulation_CheckNPart( this, NPartsOk )

    implicit none

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

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Return if no output
    if( BlockSize < 1 .and. .not. SimulationType .eq. SecondVirialCoeff ) return

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

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( SimulationType .eq. MolecularDynamics ) then
       if( .not. RootProc ) return
    endif

    ! Return if no output
    if( BlockSize < 1 ) return

    ! No output for MCOverlapReduction
    if( MCOverlapReduction ) return

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

  end subroutine TSimulation_ResultClose



!==============================================================!
!  Subroutine TSimulation_ErrorsUpdate                         !
!==============================================================!

  subroutine TSimulation_ErrorsUpdate( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( SimulationType .eq. MolecularDynamics ) then
       if( .not. RootProc ) return
    endif

    ! Return if no output
    if( BlockSize < 1 ) return

    ! Update log file
    write( IOBuffer, '("Saving simulation results")' )
    call LogWrite

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
!  Subroutine TSimulation_RDFOpen                           !
!==============================================================!

  subroutine TSimulation_RDFOpen( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Return if no output
    if( RDFUpdateFrequency < 1 ) return

    ! Open ensemble visualisation files
    do i = 1, this%NEnsembles
      call RDFOpen( this%Ensemble(i) )
    end do

  end subroutine TSimulation_RDFOpen

!==============================================================!
!  Subroutine TSimulation_RDFUpdate                         !
!==============================================================!

  subroutine TSimulation_RDFUpdate( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Return if no output
    if( RDFUpdateFrequency < 1 ) return

    ! Return if equilibration
    if( Equilibration ) return

    ! Update ensemble visualisation files
    if( mod( StepTotal - 1, RDFUpdateFrequency ) == 0 ) then
      do i = 1, this%NEnsembles
        call RDFUpdate( this%Ensemble(i) )
      end do
    end if

  end subroutine TSimulation_RDFUpdate


!==============================================================!
!  Subroutine TSimulation_RDFClose                          !
!==============================================================!

  subroutine TSimulation_RDFClose( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Return if no output
    if( RDFUpdateFrequency < 1 ) return

  end subroutine TSimulation_RDFClose
  
!==============================================================!
!  Subroutine TSimulation_RestartSave                          !
!==============================================================!

  subroutine TSimulation_RestartSave( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

  if( SimulationType .eq. SecondVirialCoeff ) return

    ! Open restart file for writing
    call FileRewrite( iounit_restart, trim( OutputNameTag )//RestartFileExtension )

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
    call MPI_Bcast( Step, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
    call MPI_Bcast( StepTotal, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
    call MPI_Bcast( Equilibration, 1, MPI_LOGICAL, NRootProc, Communicator, ierror )
    call MPI_Bcast( NVTEquilibration, 1, MPI_LOGICAL, NRootProc, Communicator, ierror )
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
