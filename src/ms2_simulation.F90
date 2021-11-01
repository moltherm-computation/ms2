!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 3.0                !
!  (c) 2017 by TU Kaiserslautern / U Paderborn                 !
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

#ifndef OSMOP
#define OSMOP 0
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
    ! first and last ensemble to be processed
    integer :: firstEnsembleIdx, lastEnsembleIdx
    ! Number of MPI ensemble groups (only relevant for MPI version, set to 0 otherwise)
    integer :: mpiEnsembleGroups

    ! Ensembles
    type(TEnsemble), pointer, contiguous :: Ensemble(:)

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

#if OSMOP > 0
  interface ProfileOpen
    module procedure TSimulation_ProfileOpen
  end interface

  interface ProfileUpdate
    module procedure TSimulation_ProfileUpdate
  end interface

  interface ProfileClose
    module procedure TSimulation_ProfileClose
  end interface
#endif

  interface RDFOpen
    module procedure TSimulation_RDFOpen
  end interface

  interface RDFUpdate
    module procedure TSimulation_RDFUpdate
  end interface

  interface RDFClose
    module procedure TSimulation_RDFClose
  end interface
  
  interface KBIOpen
    module procedure TSimulation_KBIOpen
  end interface

  interface KBIUpdate
    module procedure TSimulation_KBIUpdate
  end interface

  interface KBIClose
    module procedure TSimulation_KBIClose
  end interface
  
  interface ALPHA2Update
    module procedure TSimulation_ALPHA2Update
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

#if MPI_VER > 0
    integer  :: icommunicator
    real(RK) :: dummyR
    !integer  :: dummyI
#endif

    ! Read configuration file
#if ARCH == 1 || ARCH == 2 || ARCH == 3
    if( Restart ) then
      write( IOBuffer, '("Restarting ",A," using ",A,"_*.",A," files")' ) &
&            trim( ParameterFileName ),trim(OutputNameTag),RestartFileExtension
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
    ! Read parVersionNr
    call FileReadParameter( parVersionNr, iounit_params , IdparVersionNr, .true., 1.0_RK )
    write( IOBuffer, '("File created with/for ms2-version: ",T38, F6.3)' ) parVersionNr
    call LogWrite
    if ( parVersionNr .lt. ms2VersionNr ) then
      write( IOBuffer, '("Hint: Your ms2-version is newer than your parameter file, consider updating it.")' )
      call LogWrite
    endif

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
    call LogWriteBlank

    call FileReadParameter( max_time , iounit_params , IdWallTime , .true., 20160  )
    write( IOBuffer, '("Specified walltime: ",T23, I5, " m")' ) max_time
    call LogWrite

    call FileReadParameter( time_limit , iounit_params , IdTimeLimit , .true., 60  )
    write( IOBuffer, '("Specified time limit: ",T23, I5, " m")' ) time_limit
    call LogWrite
    call LogWriteBlank

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
      write( IOBuffer, '("Number of radial steps: ",T23, I8)' ) NSteps
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
      KBIUpdateFrequency = 0
      BlockSizeKBI = 0
      ALPHA2UpdateFrequency = 0
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
        if (.not. UseReducedUnits ) then
          if ( parVersionNr .ge. 2.0_RK ) then
            TimeStep = TimeStep / UnitTime
          else
            write( IOBuffer, '("WARNING: Time step in SI-Units was not implemented for your version of the par-file.")' )
            call LogWrite
          endif
        endif
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
      if( SimulationType .eq. MonteCarlo .and. EnsembleType .eq. EnsembleTypeHA) &
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
      end if

      ! Read number of NVT equilibration steps
      call FileReadParameter( NStepsV, iounit_params , IdNStepsV, .true., 0 )
      write( IOBuffer, '("Number of NVT equilibration steps: ",T40, I7)' ) NStepsV
      call LogWrite

      ! Read number of NVE equilibration steps
      if( EnsembleType .eq. EnsembleTypeNVE ) then
        call FileReadParameter( NStepsE, iounit_params , IdNStepsE, .true., 0 )
        write( IOBuffer, '("Number of NVE equilibration steps: ",T40, I7)' ) NStepsE
        call LogWrite
      else
        NStepsE = 0
      end if

      ! Read number of constant pressure equilibration steps
      if( ConstantPressure ) then
        if( EnsembleType .eq. EnsembleTypeHA ) then
          call FileReadParameter( NStepsP, iounit_params , IdNStepsMueP, .true., 0 )
          write( IOBuffer, '("Number of HA equilibration steps: ",T40, I7)' ) NStepsP
          call LogWrite
          
        else if( EnsembleType .eq. EnsembleTypeNPH ) then
          call FileReadParameter( NStepsH, iounit_params , IdNStepsH, .true., 0 )
          write( IOBuffer, '("Number of NPH equilibration steps: ",T40, I7)' ) NStepsH
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
      write( IOBuffer, '("Number of production steps: ",T39, I8)' ) NSteps
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
        write( IOBuffer, '("Result files will be updated each", T40, I7, " time steps")' ) BlockSize
      else
        write( IOBuffer, '("All result files will not be created")' )
      end if
      call LogWrite

      ! Calculate number of blocks and block sizes
      if( BlockSize > 0 ) then
        NBlocksMax = ceiling(max( NStepsV, NStepsE, NStepsP, NStepsH, NSteps ) / real(BlockSize))

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
            NBlocksMax = ceiling(max( NStepsV, NStepsE, NStepsP, NStepsH, NSteps ) / real(BlockSize))
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
          write( IOBuffer, '("Final result files will be updated each ", T40, I7, " time steps")' ) ErrorsUpdateFrequency
        else
          ErrorsUpdateFrequency = NSteps
          write( IOBuffer, '("Final result files will be created at the end")' )
        end if

      call LogWrite
      else
        ErrorsUpdateFrequency = NSteps
      end if

      ! Read frequency of updating visualisation file
      call FileReadParameter( VisualUpdateFrequency, iounit_params , IdVisualUpdateFrequency, .true., 0 )
      if( VisualUpdateFrequency > 0 ) then
        write( IOBuffer, '("Visualization files will be updated each", T40, I7, " time steps")' ) VisualUpdateFrequency
      else
        write( IOBuffer, '("Visualization files will not be created")' )
      end if
      call LogWrite
      call LogWriteBlank

      ! Read frequency of updating visualisation file
      call FileReadParameter( RDFUpdateFrequency, iounit_params , IdRDFUpdateFrequency, .true., 0 )
      if( RDFUpdateFrequency > 0 ) then
        write( IOBuffer, '("RDF files will be updated each", T40, I7, " time steps")' ) RDFUpdateFrequency
      else
        write( IOBuffer, '("RDF files will not be created")' )
      end if
      call LogWrite
      
      if( RDFUpdateFrequency > 0 ) then
      call FileReadParameter( RDFNumberShells, iounit_params , IdRDFNumberShells, .true., 200 )
        write( IOBuffer, '("RDF will operate with", I7, " shells")' ) RDFNumberShells
      call LogWrite
      end if
      call LogWriteBlank
      
      ! Read frequency of updating KBI file
      call FileReadParameter( KBIUpdateFrequency, iounit_params , IdKBIUpdateFrequency, .true., 0 )
      if( KBIUpdateFrequency > 0 ) then
        if( .not. EnsembleType .eq. EnsembleTypeNVT) then 
            call Error( trim( str )//' -> Kirkwood-Buff integration is in the NVT ensemble only defined' )
        else
            if (SimulationType .eq. MolecularDynamics ) KBIUpdateFrequency=1 !with MD and KBI -> KBISum is calculated while traversing the interaction matrix with RunMDStep            
            write( IOBuffer, '("RDF for KBI will be updated each", T40, I7, " time steps")' ) KBIUpdateFrequency
        end if
      else
        write( IOBuffer, '("KBI files will not be created")' )
      end if
      call LogWrite
      
      if( KBIUpdateFrequency > 0 ) then
        call FileReadParameter( BlockSizeKBI, iounit_params , IdKBIResetFrequency, .true., 10000 )
        !rounding up if KBIResetFreq is not a multiple of KBIUpdateFreq
        BlockSizeKBI = KBIUpdateFrequency*ceiling(real(BlockSizeKBI,RK)/real(KBIUpdateFrequency,RK))
        write( IOBuffer, '("RDF for KBI will be reset each", T40, I7, " time steps")' ) BlockSizeKBI
        call LogWrite
        call FileReadParameter( KBINumberShells, iounit_params , IdKBINumberShells, .true., 200 )
        write( IOBuffer, '("RDF for KBI will operate with", I7, " shells")' ) KBINumberShells
        call LogWrite
        KBINumberShellsMax=ceiling(sqrt(3*real(KBINumberShells,RK)**2))
        KBINShellsCubeEdge=floor(sqrt(2*real(KBINumberShells,RK)**2))
#if MPI_VER > 0     
        if (SimulationType .eq. MonteCarlo) then 
            BlockSizeKBI=int(BlockSizeKBI/NProcs) !KBIBlockSize per process
            !rounding up if KBIResetFreq is not a multiple of KBIUpdateFreq
            BlockSizeKBI=KBIUpdateFrequency*int(BlockSizeKBI/KBIUpdateFrequency)
        end if
#endif
        ! Calculate number of blocks and block sizes for KBI
        NBlocksMaxKBI = ceiling(max( NStepsV, NStepsE, NStepsP, NStepsH, NSteps ) / real(BlockSizeKBI))
        NBlockSizesMaxKBI = int( sqrt( real( NSteps / BlockSizeKBI, RK ) ) )
      end if
      call LogWriteBlank
      
      ! Read frequency of updating Alpha2 correlation function
      call FileReadParameter( ALPHA2UpdateFrequency, iounit_params, IdALPHA2UpdateFrequency, .true., 0 )
      if ( ALPHA2UpdateFrequency > 0 ) then
        if ( SimulationType .eq. MolecularDynamics ) then 
            call FileReadParameter( ALPHA2Length, iounit_params, IdALPHA2Length, .true., 10000 )
            call FileReadParameter( ALPHA2Shift,  iounit_params, IdALPHA2Shift,  .true., 1000  )
            write( IOBuffer, '("Alpha2 will be updated each", T40, I7, " time steps")' ) ALPHA2UpdateFrequency
            call LogWrite
            write( IOBuffer, '("Alpha2 correlation length: ", T40, I7, " time steps")' ) ALPHA2Length
            call LogWrite
            write( IOBuffer, '("Alpha2 correlation shift each", T40, I7, " time steps")' ) ALPHA2Shift
            call LogWrite           
        else
            call Error( trim( str )//' -> Alpha2 correlation function is defined for MD only' )
        end if      
      end if
      call LogWriteBlank
      
#if OSMOP > 0
      if ( SimulationType .eq. MonteCarlo ) then
        write( IOBuffer, '("Osmotic Pressure calculation with in Monte-Carlo not possible. Continuing without")' )
        call LogWrite
        call LogWriteBlank
      else
        !Number of Bins for the Density, Chem. Potential and Pressure 
        call FileReadParameter( NBinsDen, iounit_params , IdNBinsDen, .true., 500 )
        write( IOBuffer, '("Osmotic Pressure calculation with ", I7, " Bins")' ) NBinsDen
        call LogWrite
        call FileReadParameter( kForceOsmoticPressure, iounit_params , IdWallForce, .true., 41868._RK )
        if( .not. UseReducedUnits ) then
          kForceOsmoticPressure = kForceOsmoticPressure/(NAvogadro*Angstroem**2)*(UnitLength**2)/UnitEnergy
        end if
        write( IOBuffer, '("Forceconstant of the wall: ",T26, F10.5, " ?")' ) kForceOsmoticPressure
        call LogWrite
        call LogWriteBlank
      end if
#endif

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

            call FileReadParameter( nsqmax_h, iounit_params , Idnsqmax, .true. )
            write( IOBuffer, '("Ewald: NsqMax:",T20, I7)' ) nsqmax_h
            call LogWrite

            call FileReadParameter( nvecmax_h, iounit_params , IdNVecMax, .true. )
            write( IOBuffer, '("Ewald: NVecMax:",T20, I7)' ) nvecmax_h
            call LogWrite

            call FileReadParameter( nmax_h, iounit_params , IdNMax, .true. )
            write( IOBuffer, '("Ewald: NMax:",T20, I7)' ) nmax_h
            call LogWrite

!         case( 'PME', 'pme', 'SPME', 'spme')
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

    ! Read number of ensembles
    call FileReadParameter( this%NEnsembles, iounit_params , IdNEnsembles, .true., 1 )
    write( IOBuffer, '("Number of ensembles:",T24, I3)' ) this%NEnsembles
    call LogWrite
    
    this%firstEnsembleIdx=1
    this%lastEnsembleIdx=this%NEnsembles

    if( (SimulationType .eq. Gibbs .and. this%NEnsembles .ne. 2) )  &
&       call Error( trim( SimulationTypeString )//" simulation of " &
&       //trim( SimulationTypeString )//" needs 2 Ensembles" )

    ! Read number of MPI ensemble groups
    call FileReadParameter( this%mpiEnsembleGroups, iounit_params , IdmpiEnsembleGroups, .true., 0 )
#if MPI_VER > 0
    write( IOBuffer, '("mpiEnsembleGroups:",T24, I3)' ) this%mpiEnsembleGroups
    call LogWrite
    if ( this%mpiEnsembleGroups .eq. 1 ) this%mpiEnsembleGroups=this%NEnsembles
    if ( this%mpiEnsembleGroups .gt. this%NEnsembles .or. this%mpiEnsembleGroups .gt. NProcs_W ) &
&      this%mpiEnsembleGroups=min(this%NEnsembles,NProcs_W)
    if ( this%mpiEnsembleGroups .le. 1 ) this%mpiEnsembleGroups=0
    
    if (this%mpiEnsembleGroups .gt. 1) then
      ! Close the ParameterFile to reopen it within the subcommunicators
      call FileClose( iounit_params )
      call MPI_Bcast( ParameterFileName, FileNameLength, MPI_CHARACTER, NRootProc, Communicator, ierror )
      ! create subcommunicators to process subranges of the ensembles ++++++++++++++++++++++++++++++
      call SplitCommunicator(this%mpiEnsembleGroups)    ! setting NCommunicator, NCommunicators and Communicator etc
      ! 1-index based
      this%firstEnsembleIdx=this%NEnsembles*NCommunicator/NCommunicators+1
      this%lastEnsembleIdx=this%NEnsembles*(NCommunicator+1)/NCommunicators
      write( IOBuffer, '("MPI communicator",I3," (out of",I3,") with ",I3," PEs computes ensemble",I3," -",I3)' ) &
&            NCommunicator+1,NCommunicators,NProcs,this%firstEnsembleIdx,this%lastEnsembleIdx
      ! be aware that e.g. the random number generator calls might be different
      call LogWrite
      ! Reopen the ParameterFile (dirty hack) for each communicator
      call FileReset( iounit_params, ParameterFileName )
      !call FileReadParameter( dummyI, iounit_params , IdNEnsembles, .true., 1 )
    endif
#else
    if ( this%mpiEnsembleGroups /= 0 ) then
      write( IOBuffer, '("Warning: mpiEnsembleGroups only supported in MPI version")')
      call LogWrite
      write( IOBuffer, '("         neglecting mpiEnsembleGroups =",I6)') this%mpiEnsembleGroups
      call LogWrite
      this%mpiEnsembleGroups=0
    end if
#endif

    ! Create ensembles
    allocate( this%Ensemble(this%firstEnsembleIdx:this%lastEnsembleIdx), STAT = stat )
    call AllocationError( stat, 'ensembles', this%NEnsembles )

#if  TRANS == 1
!TRANSPORT_start
    ! Read correlation function mode
    if ( parVersionNr .lt. 2.0_RK ) then
      call FileReadParameter( str , iounit_params , IdCorrFun, .true. , 'no' )
      select case( str )

      case( 'yes' , 'ok', 'ja' )
        this%Ensemble(:)%CorrFunMode = .true.
        str = 'Include transport properties for all ensembles'

      case( 'no', 'nein' )
        this%Ensemble(:)%CorrFunMode = .false.
        str = 'No transport properties for any ensemble'
        call Error( 'Use a binary compiled without -DTRANS if you do not &
&                    wish to calculate transport properties. If you do, set CorrFunMode = yes ' )

      case default
        call Error( 'Unknown transport properties ('//trim(IdCorrFun)//'='//trim(str)//')' )
      end select
      write( IOBuffer, '("Transport properties:",T26, A)' ) trim(str)
      call LogWrite
    endif
!TRANSPORT_END
#endif

#if MPI_VER > 0
    ! force sequential reading of parameter file (within Ensemble Construct)    better use MPI-IO!
    do icommunicator = 0,NCommunicators-1
      if (icommunicator==NCommunicator) then
#endif
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
      if (LongRange .eq. Ewald) then
            this%Ensemble(i)%KappaL = KappaL_h
            this%Ensemble(i)%nsqmax = nsqmax_h
            this%Ensemble(i)%nvecmax = nvecmax_h
            this%Ensemble(i)%nmax = nmax_h

#if SPME > 0
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
#if MPI_VER > 0
      else
        ! dirty hack to move the file pointer to the next ensemble
        call FileReadParameter( dummyR, iounit_params , IdRefTemperature, .false. )
        !call FileReadParameter( dummyR, iounit_params , IdRefDensity, .false. )
      end if
      call MPI_Barrier( MPI_COMM_WORLD, ierror )
    end do
#endif

#if  TRANS == 1
!TRANSPORT_start
    ! Read correlation function mode
    if ( parVersionNr .ge. 2.0_RK ) then
      if ( .not. ANY(this%Ensemble(:)%CorrFunMode) ) then
        str = 'No transport properties for any ensemble'
        call Error( 'Use a binary compiled without -DTRANS if you do not &
&                    wish to calculate transport properties. If you do, set CorrFunMode = yes for one ensemble ' )

        write( IOBuffer, '("Transport properties:",T26, A)' ) trim(str)
        call LogWrite
      endif
    endif
!TRANSPORT_END
#endif

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
    call KBIOpen( this )
#if OSMOP > 0
    if ( SimulationType .ne. MonteCarlo ) call ProfileOpen(this )
#endif

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
    call KBIClose( this )
    !call RDFClose( this ) ! file is closed after updating
#if OSMOP > 0
    if ( SimulationType .ne. MonteCarlo ) call ProfileClose(this )
#endif
    ! Destroy accumulators
    call DestroyAccumulators( this )

    ! Destroy ensembles
    if( associated( this%Ensemble ) ) then
      do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
    integer :: i, j, s, t, NGradInsInit
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
    logical :: AnyNPartOk = .false.
#endif 

    tooManyParticles = .false.
    call Construct(RunTimer,"TSimulation_Run",CStopwatch_doMPIStartBarrier)
    call Construct(RunStepsTimer)
#if MPI_VER > 0
    call Timer_SetMPIcommunicator(RunTimer,MPI_COMM_WORLD)
    call Timer_SetMPIcommunicator(RunStepsTimer,Communicator)
#endif

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
          
            call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,NProc_W,Communicator,ierror) 
            call SetCommunicator( Communicator )             
          endif    
        else
          if (NProcs .gt. 4) then
            color = 0
            do i=1,lengthHost
              tmpVal = ichar (trim(hostnameStr(i:i)))
              color = color + (tmpVal**2)*i
            enddo

            call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,NProc_W,Communicator,ierror) 
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
                if( EnsembleType .eq. EnsembleTypeGE ) NStepsP = 1
                if( ConstantPressure ) then 
                  if(EnsembleType .eq. EnsembleTypeNPH ) then
                    NStepsH = 1
                  else 
                    NStepsP =  1
                  end if
                end if
                if( EnsembleType .eq. EnsembleTypeNVE ) NStepsE = 1
              endif
           endif
           if (RootProc) then
             if (NProc_W .ne. NRootProc) then
               write( IOBuffer, '(I16)' ) NProc_W  
               call FileRewrite( iounit_log, trim( OutputNameTag )//'_Equi_'//trim( adjustl( IOBuffer ) )//LogFileExtension )
               
               do j = this%firstEnsembleIdx, this%lastEnsembleIdx

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
      do j = this%firstEnsembleIdx, this%lastEnsembleIdx
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
      call Timer_setTag(RunStepsTimer,"MC overlap reduction")
      call start_Timer(RunStepsTimer)
      call logwritestart_Timer(RunStepsTimer)

      call RunSteps( this, StepStart, StepEnd )
      
      call stop_Timer(RunStepsTimer)
      call logwritestop_Timer(RunStepsTimer)

      if( .not. TerminateProgram ) then
        write( IOBuffer, '("MC overlap reduction completed")' )
        MCOverlapReduction = .false.
        SimulationType = MolecularDynamics

        do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
      ! Run GE, NpT or NVE equilibration
      if( Equilibration .and. .not. TerminateProgram ) then
        if( EnsembleType .eq. EnsembleTypeGE ) then
          StepEnd = NStepsP

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
#if MPI_VER > 0 && ( ARCH == 1 || ARCH == 2 )
            call MPI_Allreduce( NPartsOk, AnyNPartOk, 1, MPI_LOGICAL, MPI_LAND, Communicator, ierror )
            if ( .not. AnyNPartOk) then
                NPartsOk = .false.
            endif
#endif

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
          StepEnd = NStepsP
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

        else if( EnsembleType .eq. EnsembleTypeNPT ) then
          StepEnd = NStepsP
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
  
        else if( EnsembleType .eq. EnsembleTypeNPH ) then
          StepEnd = NStepsH
          call LogWriteBlank
          if( Restart ) then
            write( IOBuffer, '("Resuming NPH equilibration")' )
            Restart = .false.
          else
            write( IOBuffer, '("Starting NPH equilibration")' )
          end if

          call Timer_setTag(RunStepsTimer,"NPH equilibration")
          call start_Timer(RunStepsTimer)
          call logwritestart_Timer(RunStepsTimer)
          call RunSteps( this, StepStart, StepEnd )
          call stop_Timer(RunStepsTimer)
          call logwritestop_Timer(RunStepsTimer)

          if( .not. TerminateProgram ) then
            write( IOBuffer, '("NPH equilibration completed")' )
            Equilibration = .false.
          else
            write( IOBuffer, '("NPH equilibration TERMINATED")' )
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

        else if( EnsembleType .eq. EnsembleTypeNVE ) then
          StepEnd = NStepsE
          call LogWriteBlank
          if( Restart ) then
            write( IOBuffer, '("Resuming NVE equilibration")' )
            Restart = .false.
          else
            write( IOBuffer, '("Starting NVE equilibration")' )
          end if

          call Timer_setTag(RunStepsTimer,"NVE equilibration")
          call start_Timer(RunStepsTimer)
          call logwritestart_Timer(RunStepsTimer)
          call RunSteps( this, StepStart, StepEnd )
          call stop_Timer(RunStepsTimer)
          call logwritestop_Timer(RunStepsTimer)

          if( .not. TerminateProgram ) then
            write( IOBuffer, '("NVE equilibration completed")' )
            Equilibration = .false.
          else
            write( IOBuffer, '("NVE equilibration TERMINATED")' )
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
      
      do k = this%firstEnsembleIdx, this%lastEnsembleIdx
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
      
      
      if (multNodes) then
      
        if (RootProc) then
          if (NProc_W .ne. NRootProc) then
            ! Close all files keeping track of the equilibration
            do j = this%firstEnsembleIdx, this%lastEnsembleIdx
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
             do j = this%firstEnsembleIdx, this%lastEnsembleIdx
               call MPI_Bcast( this%Ensemble(j)%EPot, 1, MPI_RK, NRootProc, Communicator, ierror )
               call MPI_Bcast( this%Ensemble(j)%DispVol, 1, MPI_RK, NRootProc, Communicator, ierror )            
               do i = 1, this%Ensemble(j)%NComponents
                 call MPI_Bcast( this%Ensemble(j)%Component(i)%P0(:, :), size( this%Ensemble(j)%Component(i)%P0 ), &
&                     MPI_RK, NRootProc, Communicator, ierror )

                 if( this%Ensemble(j)%Component(i)%Molecule%isElongated ) then
                    call MPI_Bcast( this%Ensemble(j)%Component(i)%Q0(:, :), size( this%Ensemble(j)%Component(i)%Q0 ), &
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
      do j = this%firstEnsembleIdx, this%lastEnsembleIdx
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
     do j = this%firstEnsembleIdx, this%lastEnsembleIdx
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
       
       NGradInsInit = 1      
       do j= this%firstEnsembleIdx, this%lastEnsembleIdx  
        do i = 1, this%Ensemble(j)%NComponents
         NGradInsInit = NGradInsInit + this%Ensemble(j)%Component(i)%GradInsInit
        end do 
       end do
      
       do j= this%firstEnsembleIdx, this%lastEnsembleIdx
        do Step = StepStart, NGradInsInit
        call ChemicalPotential( this%Ensemble(j) )
        end do 
       end do
             
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
      if( BlockSizeKBI > 0 ) NBlocksKBI = 1 + (Step - 1) / BlockSizeKBI
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
#if MPI_VER > 0
    integer :: mpistatus(MPI_STATUS_SIZE)
    integer :: TerminateStatus
    integer :: mpireqbcastTerm, mpireqmsgTerm
    logical :: doneBcastTerm, doneMsgTerm
    integer :: numMsgTerm_send, numMsgTerm_recv
#endif

#if TRANS==1
    integer:: StepCF
#endif

    integer:: o, i, j, t, s

#if MPI_VER > 0
   if (NCommunicators > 1 ) then
     TerminateStatus=0
     doneBcastTerm=.false.
     doneMsgTerm=.false.
     numMsgTerm_send=0
     numMsgTerm_recv=0
     if ( RootProc) then
       if ( RootProc_R ) then
         ! RootProc_W subcommunicator root starts receiving a TerminateStatus message
         call MPI_Irecv(TerminateStatus, 1, MPI_INTEGER, MPI_ANY_SOURCE, mpimsgtag_simTerm, Communicator_R, mpireqmsgTerm, ierror)
       else ! (RootProc.and.).not.RootProc_R
         ! non_RootProc_R subcommunicator roots start receiving TerminateStatus broadcast of TerminateStatus before the loop
         call MPI_Ibcast(TerminateStatus, 1, MPI_INTEGER, NRootProc_R, Communicator_R, mpireqbcastTerm, ierror)
       end if
     end if
   end if
#endif

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
        if ( .not. Equilibration ) then
          do i = this%firstEnsembleIdx, this%lastEnsembleIdx
            if ( this%Ensemble(i)%CorrFunMode ) then
              if (mod((Step+this%Ensemble(i)%NStepCorr-1),this%Ensemble(i)%NStepCorr) .eq. 0) then
                StepCF = (Step + this%Ensemble(i)%NStepCorr -1) / this%Ensemble(i)%NStepCorr
                if ( StepCF >= this%Ensemble(i)%Ncorr )then
                  NBlocksCF = 1 + int(( StepCF - 1 - this%Ensemble(i)%Ncorr ) / &
&                                            ( BlockSizeCF * this%Ensemble(i)%NSpancf ))
                  NBlockSizesCF = int( sqrt( real(( StepCF - this%Ensemble(i)%Ncorr) / &
&                                                ( BlockSizeCF * this%Ensemble(i)%NSpancf ), RK)))
                else
                  NBlocksCF     = 0
                  NBlockSizesCF = 0
                end if
              end if

            end if
          end do
        end if
#endif
      end if
      
      ! Set current block number KBI
      if( BlockSizeKBI > 0 ) then
        NBlocksKBI = 1 + (Step - 1) / BlockSizeKBI
        NBlockSizesKBI = int( sqrt( real( Step / BlockSizeKBI, RK ) ) )
      end if
      

      ! Run simulation step
      select case( SimulationType )
      case( MolecularDynamics )
        call RunMDStep( this )
        call ALPHA2Update( this )
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
      call KBIUpdate ( this )

      ! Update log and result files
      if( mod( Step, LogUpdateFrequency ) == 0 .or. Step == StepEnd ) call LogWriteStep
      if( .not. Equilibration .and. ( mod( Step, ErrorsUpdateFrequency ) == 0 .or. Step == StepEnd )) then
        call ErrorsUpdate( this )
#if OSMOP > 0
        if ( SimulationType .ne. MonteCarlo ) call ProfileUpdate(this )
#endif
      endif

      ! Check for termination request (caused by signal handler)
#if MPI_VER > 0
      if (NCommunicators > 1 ) then
        ! transfer termination information to TerminateStatus, delete the flags and wait for the broadcast...
        if (TerminateProgram) TerminateStatus=IOR(TerminateStatus,1)
        if (tooManyParticles) TerminateStatus=IOR(TerminateStatus,2)
        ! terminate solely after the terminate broadcast was received
        TerminateProgram= .false.
        tooManyParticles= .false.
        !call MPI_Allreduce( MPI_IN_PLACE, TerminateStatus, 1, MPI_INTEGER, MPI_BOR, Communicator, ierror )
        if ( RootProc ) then
          call MPI_Reduce( MPI_IN_PLACE, TerminateStatus, 1, MPI_INTEGER, MPI_BOR, NRootProc, Communicator, ierror )
          if ( .not. doneMsgTerm ) then
            if ( RootProc_W ) then
              call MPI_Test(mpireqmsgTerm, doneMsgTerm, mpistatus, ierror)
              if ( doneMsgTerm ) then
                write( IOBuffer, '("received message with termination status (",B0,") within step ",I0,"/",I0)' ) &
&                      TerminateStatus, Step, StepTotal
                call LogWriteTime
                doneMsgTerm=.true.
                numMsgTerm_recv = numMsgTerm_recv + 1
              end if
            else ! (RootProc.and.).not.RootProc_W
              if (TerminateStatus /= 0) then
                write( IOBuffer, '("sending message with termination status (",B0,") from PE",I0," within step ",I0,"/",I0)' ) &
&                      NProc_W, TerminateStatus, Step, StepTotal
                call LogWriteTime
                call MPI_ISend(TerminateStatus, 1, MPI_INTEGER, NRootProc_R, mpimsgtag_simTerm, Communicator_R, mpireqmsgTerm, ierror)
                doneMsgTerm=.true.
                numMsgTerm_send = numMsgTerm_send + 1
              end if
            end if
          end if
          if ( .not. doneBcastTerm ) then
            if ( RootProc_R ) then
              if (TerminateStatus /= 0) then
                write( IOBuffer, '("broadcasting termination status (",B0,") within step ",I0,"/",I0)' ) &
&                      TerminateStatus, Step, StepTotal
                call LogWriteTime
                call MPI_Ibcast(TerminateStatus, 1, MPI_INTEGER, NRootProc_R, Communicator_R, mpireqbcastTerm, ierror)
                doneBcastTerm = .true.
              end if
            else
              call MPI_Test(mpireqbcastTerm, doneBcastTerm, mpistatus, ierror)
              if (doneBcastTerm .and. TerminateStatus>0) then
                !TerminateProgram=.true.
                !if (IAND(TerminateStatus,2).eq.2) tooManyParticles=.true.
                write( IOBuffer, '("received broadcast with termination status (",B0,") within step ",I0,"/",I0)' ) &
&                      TerminateStatus, Step, StepTotal
                call LogWriteTime
              end if
            end if
          end if
        else !.not.RootProc
          call MPI_Reduce( TerminateStatus, TerminateStatus, 1, MPI_INTEGER, MPI_BOR, NRootProc, Communicator, ierror )
        end if
        
        ! broadcast TerminateStatus within the subcommunicator
        call MPI_Bcast(TerminateStatus, 1, MPI_INTEGER, NRootProc, Communicator, ierror)
        if (TerminateStatus > 0) TerminateProgram=.true.
        if (IAND(TerminateStatus,2).eq.2) tooManyParticles=.true.
      else
        !                                                                            Communicator
        call MPI_Allreduce( MPI_IN_PLACE, TerminateProgram, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror )
        call MPI_Allreduce( MPI_IN_PLACE, tooManyParticles, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror )
      end if
#endif

      if ( TerminateProgram ) then
        write( IOBuffer, '("terminating program after step ",I0,"/",I0)' ) Step,StepTotal
        call LogWriteTime
        exit
      end if
      ! Check for too many particles (GE only)
      if ( tooManyParticles ) exit

    end do

#if MPI_VER > 0
    if (NCommunicators > 1 ) then
      ! clean up (but don't use MPI_Cancel)
      if ( RootProc ) then
        if ( RootProc_R ) then
          call MPI_Reduce( MPI_IN_PLACE, numMsgTerm_send, 1, MPI_INTEGER, MPI_SUM, NRootProc_R, Communicator_R, ierror )
!          if ( .not. doneMsgTerm ) then
!            ! check again, if terminate message was received
!            call MPI_Test(mpireqmsgTerm, doneMsgTerm, mpistatus, ierror)
!            if ( doneMsgTerm ) then
!              write( IOBuffer, '("received message with termination status (",B0,") after step ",I0,"/",I0)' ) &
!&                    TerminateStatus, Step, StepTotal
!              call LogWriteTime
!              doneMsgTerm=.true.
!              numMsgTerm_recv = numMsgTerm_recv + 1
!            end if
!          end if
          !                             1 irecv is received or pending
          do i = 1, numMsgTerm_send-max(numMsgTerm_recv,1)
            call MPI_Recv(TerminateStatus, 1, MPI_INTEGER, MPI_ANY_SOURCE, mpimsgtag_simTerm, Communicator_R, ierror)
            if (IAND(TerminateStatus,1).eq.1) TerminateProgram=.true.
            if (IAND(TerminateStatus,2).eq.2) tooManyParticles=.true.
          end do
        else ! .not.RootProc_R
          !if ( .not. doneMsgTerm .and. NProc_R.eq.1 ) then ! only works if NRootProc_R.ne.1 (NRootProc_R==0)
          if ( .not. doneMsgTerm .and. NProc_R.eq.mod(NRootProc_R+1,NProcs_R) ) then    ! should work for NProcs_R.gt.1
            ! at least one terminate message should be sent to serve the RootProc_R irecv - e.g. NProc==1
              write( IOBuffer, '("sending message with termination status (",B0,") from PE",I0," after step ",I0,"/",I0)' ) &
&                    NProc_W, TerminateStatus, Step, StepTotal
              call LogWriteTime
              call MPI_ISend(TerminateStatus, 1, MPI_INTEGER, NRootProc_R, mpimsgtag_simTerm, Communicator_R, mpireqmsgTerm, ierror)
              doneMsgTerm=.true.
              numMsgTerm_send = numMsgTerm_send + 1
          end if
          call MPI_Reduce( numMsgTerm_send, numMsgTerm_send, 1, MPI_INTEGER, MPI_SUM, NRootProc_R, Communicator_R, ierror )
        end if
        if ( doneMsgTerm ) then
          call MPI_Wait(mpireqmsgTerm, mpistatus, ierror)
        end if
        
        if ( .not. doneBcastTerm ) then
          if ( RootProc_R ) then
!            write( IOBuffer, '("broadcasting termination status (",B0,") after step ",I0,"/",I0)' ) &
!&                  TerminateStatus, Step, StepTotal
!            call LogWriteTime
            call MPI_Ibcast(TerminateStatus, 1, MPI_INTEGER, NRootProc_R, Communicator_R, mpireqbcastTerm, ierror)
            doneBcastTerm = .true.
!          else
!            call MPI_Test(mpireqbcastTerm, doneBcastTerm, mpistatus, ierror)
!            if (doneBcastTerm .and. TerminateStatus>0) then
!              write( IOBuffer, '("received broadcast with termination status (",B0,") after step ",I0,"/",I0)' ) &
!&                    TerminateStatus, Step, StepTotal
!              call LogWriteTime
!            end if
          end if
          call MPI_Wait(mpireqbcastTerm, mpistatus, ierror)
        end if
      end if
      
      if (TerminateProgram) TerminateStatus=IOR(TerminateStatus,1)
      if (tooManyParticles) TerminateStatus=IOR(TerminateStatus,2)
      call MPI_Allreduce( MPI_IN_PLACE, TerminateStatus, 1, MPI_INTEGER, MPI_BOR, MPI_COMM_WORLD, ierror )
      if (IAND(TerminateStatus,1).eq.1) TerminateProgram=.true.
      if (IAND(TerminateStatus,2).eq.2) tooManyParticles=.true.
    end if
#endif

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
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
    if ( (SimulationType .eq. Gibbs) .and. ConstantPressure ) &
&       call Error( 'Gibbs Ensemble only implemented for NVT')
    if ( (SimulationType .eq. Gibbs) .and. (this%NEnsembles .ne. 2) ) &
&       call Error( 'Gibbs Ensemble needs two SimBoxes: one liquid and one vapor SimBox')

    ! Run MC simulation step
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
#if MPI_VER > 0
    if(SimulationType .ne. MonteCarlo) then
      if( .not. RootProc ) return
    endif
#else
    ! Check for root process
    if( .not. RootProc ) return
#endif
    ! Return if no output
    if( BlockSize < 1 .and. .not. SimulationType .eq. SecondVirialCoeff ) return

    ! Open ensemble result files
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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

#if MPI_VER > 0
    if(SimulationType .ne. MonteCarlo) then
      if( .not. RootProc ) return
    endif
#else
    ! Check for root process
    if( .not. RootProc ) return
#endif
    ! Return if no output
    if( BlockSize < 1 .and. .not. SimulationType .eq. SecondVirialCoeff ) return

    ! Close ensemble result files
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
      do i = this%firstEnsembleIdx, this%lastEnsembleIdx
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
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
      call VisualClose( this%Ensemble(i) )
    end do

  end subroutine TSimulation_VisualClose


#if OSMOP > 0
!==============================================================!
!  Subroutine TSimulation_ProfileOpen                          !
!==============================================================!

  subroutine TSimulation_ProfileOpen( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Open ensemble visualisation files
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
      call ProfileOpen( this%Ensemble(i) )
    end do

  end subroutine TSimulation_ProfileOpen


!==============================================================!
!  Subroutine TSimulation_ProfileUpdate                        !
!==============================================================!

  subroutine TSimulation_ProfileUpdate( this )

    implicit none

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

    ! Update ensemble visualisation files
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
       call ProfileUpdate( this%Ensemble(i) )
    end do

  end subroutine TSimulation_ProfileUpdate


!==============================================================!
!  Subroutine TSimulation_ProfileClose                         !
!==============================================================!

  subroutine TSimulation_ProfileClose( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Close ensemble visualisation files
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
      call ProfileClose( this%Ensemble(i) )
    end do

  end subroutine TSimulation_ProfileClose
#endif


!==============================================================!
!  Subroutine TSimulation_RDFOpen                              !
!==============================================================!

  subroutine TSimulation_RDFOpen( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Return if no output
    if( RDFUpdateFrequency < 1 ) return

    ! Open ensemble visualisation files
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
      call RDFOpen( this%Ensemble(i) )
    end do

  end subroutine TSimulation_RDFOpen

!==============================================================!
!  Subroutine TSimulation_RDFUpdate                            !
!==============================================================!

  subroutine TSimulation_RDFUpdate( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Return if no output
    if( RDFUpdateFrequency < 1 ) return

    ! Return if equilibration
    if( Equilibration ) return

    ! Update ensemble visualisation files
    if( mod( Step - 1, RDFUpdateFrequency ) == 0 ) then
      do i = this%firstEnsembleIdx, this%lastEnsembleIdx
        call RDFUpdate( this%Ensemble(i) )
      end do
    end if

  end subroutine TSimulation_RDFUpdate


!==============================================================!
!  Subroutine TSimulation_RDFClose                             !
!==============================================================!

  subroutine TSimulation_RDFClose( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Return if no output
    if( RDFUpdateFrequency < 1 ) return

    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
      call RDFClose( this%Ensemble(i) )
    end do

  end subroutine TSimulation_RDFClose


!==============================================================!
!  Subroutine TSimulation_KBIOpen                              !
!==============================================================!

  subroutine TSimulation_KBIOpen( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Return if no output
    if( KBIUpdateFrequency < 1 ) return

    ! Open ensemble visualisation files
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
      call KBIOpen( this%Ensemble(i) )
    end do

  end subroutine TSimulation_KBIOpen

!==============================================================!
!  Subroutine TSimulation_KBIUpdate                            !
!==============================================================!

  subroutine TSimulation_KBIUpdate( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Return if no output
    if( KBIUpdateFrequency < 1 ) return

    ! Return if equilibration
    if( Equilibration ) return

    ! Update ensemble visualisation files
    if( mod( Step - 1, KBIUpdateFrequency ) == 0 ) then
      do i = this%firstEnsembleIdx, this%lastEnsembleIdx
        call KBIUpdate( this%Ensemble(i) )
      end do
    end if

  end subroutine TSimulation_KBIUpdate


!==============================================================!
!  Subroutine TSimulation_KBIClose                             !
!==============================================================!

  subroutine TSimulation_KBIClose( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Return if no output
    if( KBIUpdateFrequency < 1 ) return

    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
      call KBIClose( this%Ensemble(i) )
    end do

  end subroutine TSimulation_KBIClose
  
!==============================================================!
!  Subroutine TSimulation_ALPHA2Update                         !
!==============================================================!

  subroutine TSimulation_ALPHA2Update( this )

    implicit none

    ! Declare arguments
    type(TSimulation) :: this

    ! Declare local variables
    integer :: i

    ! Return if no output
    if( ALPHA2UpdateFrequency < 1 ) return

    ! Return if equilibration
    if( Equilibration ) return

    ! Update ensemble visualisation files
    if( mod( Step - 1, ALPHA2UpdateFrequency ) == 0 ) then
      do i = this%firstEnsembleIdx, this%lastEnsembleIdx
        call ALPHA2Update( this%Ensemble(i) )
      end do
    end if

  end subroutine TSimulation_ALPHA2Update
  
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
    integer :: i,j

    ! Check for root process
    if( RootProc ) then

        if( SimulationType .eq. SecondVirialCoeff ) return

        write( RestartFileName, '(A,A)' ) trim(OutputNameTag),RestartFileExtension
#if MPI_VER > 0
        if ( NCommunicators .gt. 1 ) then
          write( RestartFileName, '(A,"_",I0,A)' ) trim(OutputNameTag),NCommunicator+1,RestartFileExtension
        endif
#endif

        write( IOBuffer, '("Saving restart file ", A)' ) trim( RestartFileName )
        call LogWriteTime

        ! Open restart file for writing
        call FileRewrite( iounit_restart, trim(RestartFileName) )

        ! Save contents to restart file
        write( iounit_restart, '(A)' ) trim( ParameterFileName )
        write( iounit_restart, '(2I10)' ) Step, StepTotal
        write( IOBuffer, '("saving restart data at step",I10," (of",I10,")")' ) Step, StepTotal
        call LogWrite
        write( iounit_restart, '(2L5)' ) Equilibration, NVTEquilibration

    end if
    
    ! Save ensembles
    do i = this%firstEnsembleIdx, this%lastEnsembleIdx
        if( RootProc ) then
            write( IOBuffer, '("writing ensemble",I7)' ) i
            call LogWriteTime
            write( iounit_restart, '(A,":",I0)' ) RstEnsembleMarker,i
        end if
        ! saving ensemble data
        call RestartSave( this%Ensemble(i) )
    end do  
    
    ! Close restart file
    call FileClose( iounit_restart )

    write( IOBuffer, '("Finished saving restart file ", A)' ) trim( RestartFileName )
    call LogWriteTime

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
    character(FileNameLength) :: parfilename
    character(IOBufferLength) :: ensemblemarker
    integer :: pos
#if MPI_VER > 0
    integer :: filepos
    integer :: stat
#endif

    write( IOBuffer, '("Reading restart file ")' )
    call LogWriteTime

    if( RootProc ) then

      write( RestartFileName, '(A,A)' ) trim(OutputNameTag),RestartFileExtension
#if MPI_VER > 0
      if ( NCommunicators .gt. 1 ) then
        write( RestartFileName, '(A,"_",I0,A)' ) trim(OutputNameTag),NCommunicator+1,RestartFileExtension
      endif
#endif
      call FileReset( iounit_restart, trim(RestartFileName) )

      ! Read non-ensemble specifif contents from restart file first
      read( iounit_restart, '(A128)' ) parfilename
      if (trim(parfilename) /= trim(ParameterFileName)) then
        call LogWriteBlank
        write( IOBuffer, '("WARNING: ",A," was created with par-file ",A," and NOT ",A)' ) &
&                        trim(RestartFileName), trim(parfilename), trim(ParameterFileName)
        call LogWrite
        call LogWriteBlank
      endif
      read( iounit_restart, '(2I10)' ) Step, StepTotal
      write( IOBuffer, '("restarting at step",I10," (of",I10,")")' ) Step, StepTotal
      call LogWrite
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
    
    ! Set current block number KBI
    if( BlockSizeKBI > 0 ) then
      NBlocksKBI = 1 + (Step - 1) / BlockSizeKBI
      NBlockSizesKBI = int( sqrt( real( Step / BlockSizeKBI, RK ) ) )
    end if

      ! Read ensembles
      do i = this%firstEnsembleIdx, this%lastEnsembleIdx
        if( RootProc ) then
          read( iounit_restart, '(A)' ) ensemblemarker
          pos = index( ensemblemarker,':')
          if( pos<=1 .or. trim(ensemblemarker(1:pos-1))/=trim(RstEnsembleMarker) ) then
            call LogWriteBlank
            write( IOBuffer, '("WARNING: expected marker ",A," but read ",A)' ) trim(RstEnsembleMarker), trim(ensemblemarker)
            call LogWrite
            call LogWriteBlank
          end if
          write( IOBuffer, '("reading ensemble",I6," (marker ",A,")")' ) i, trim(ensemblemarker)
          call LogWriteTime
        end if
        ! reading ensemble data
        call RestartRead( this%Ensemble(i) )
      end do
      
      ! Close restart file
      call FileClose( iounit_restart )

    write( IOBuffer, '("Finished reading restart file", A)' ) trim( RestartFileName )
    call LogWriteTime
    
 end subroutine TSimulation_RestartRead



end module ms2_simulation
