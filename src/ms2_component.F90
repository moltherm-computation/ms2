!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 3.0                !
!  (c) 2017 by TU Kaiserslautern / U Paderborn                 !
!      P.O. Box 67653                                          !
!      67653 Kaiserslautern                                    !
!==============================================================!
!  Module ms2_component                                        !
!  Contains TComponent object                                  !
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
!DEC$ MESSAGE:'Compiling ms2_component.F90...'
#endif

module ms2_component

  use ms2_accumulator
  use ms2_global
  use ms2_molecule
  use ms2_site



!==============================================================!
!  Type TComponent                                             !
!==============================================================!

  type TComponent

    ! Charged component
    logical           :: charged

#if OSMOP > 0
    ! Permeability
    logical           :: permeable
#endif

    ! Positions and orientations of test particles
    real(RK), pointer, contiguous :: P0Test(:, :), Q0Test(:, :)

    ! Centers of mass positions and their derivatives
    real(RK), pointer, contiguous :: P0(:, :)
    real(RK), pointer, contiguous :: P0Save(:, :)
    real(RK), pointer, contiguous :: P0old(:, :)
    real(RK), pointer, contiguous :: P1(:, :)
    real(RK), pointer, contiguous :: P2(:, :)
    real(RK), pointer, contiguous :: P3(:, :)
    real(RK), pointer, contiguous :: P4(:, :)
    real(RK), pointer, contiguous :: P5(:, :)

    ! Quaternion parameters and their derivatives
    real(RK), pointer, contiguous :: Q0(:, :)
    real(RK), pointer, contiguous :: Q0Save(:, :)
    real(RK), pointer, contiguous :: Q0tmp(:, :)
    real(RK), pointer, contiguous :: Q1(:, :)
    real(RK), pointer, contiguous :: Q2(:, :)
    real(RK), pointer, contiguous :: Q3(:, :)
    real(RK), pointer, contiguous :: Q4(:, :)

    ! Angular velocities and their derivatives
    real(RK), pointer, contiguous :: W0(:, :)
    real(RK), pointer, contiguous :: W1(:, :)
    real(RK), pointer, contiguous :: W2(:, :)
    real(RK), pointer, contiguous :: W3(:, :)
    real(RK), pointer, contiguous :: W4(:, :)

    ! Displacement
    real(RK), pointer, contiguous :: Disp(:, :)

    ! Alpha2 matrix
    real(RK), pointer, contiguous :: ri0_x(:, :)
    real(RK), pointer, contiguous :: ri0_y(:, :)
    real(RK), pointer, contiguous :: ri0_z(:, :)

#if TRANS==1
    !EinsteinCoef  ri0_E Component
    real(RK), pointer, contiguous :: ri0_E_x(:, :)
    real(RK), pointer, contiguous :: ri0_E_y(:, :)
    real(RK), pointer, contiguous :: ri0_E_z(:, :)
    integer :: NEinstein !  NCorr / NSpanCF
#endif

    ! Total forces
    real(RK), pointer, contiguous :: F(:, :)
#if MPI_VER > 0
    real(RK), pointer, contiguous :: FAll(:, :)
#endif

    ! Total torques
    real(RK), pointer, contiguous :: T(:, :)
#if MPI_VER > 0
    real(RK), pointer, contiguous :: TAll(:, :)
#endif

#if OSMOP > 0
    ! Force for osmotic pressure
    real(RK), pointer, contiguous :: FOsmoticPressure(:)

    ! Density profile
    integer, pointer, contiguous  :: DensityProfileN(:)
#if OSMOP == 2
    real(RK), pointer, contiguous :: ChemPotProfile(:)
#endif
#endif

    ! NPart additional
#if MPI_VER > 0
    logical, pointer, contiguous :: NAdd(:)
#endif

#if  TRANS == 1
!TRANSPORT_start
    real(RK), pointer, contiguous :: KinETran(:,:)
    real(RK), pointer, contiguous :: KinEPart(:)
    real(RK) :: KinETranTotal(3)
    real(RK) :: PartialMolarEnthalpy

    real(RK), pointer, contiguous :: FS(:,:)
    real(RK), pointer, contiguous :: FB(:,:)
    real(RK), pointer, contiguous :: FTC(:,:)
    real(RK), pointer, contiguous :: FRC(:,:)

    real(RK), pointer, contiguous :: FTC1(:,:)
    real(RK), pointer, contiguous :: FTC2(:,:)
    real(RK), pointer, contiguous :: FTC3(:,:)

    real(RK), pointer, contiguous :: FRC1(:,:)
    real(RK), pointer, contiguous :: FRC2(:,:)
    real(RK), pointer, contiguous :: FRC3(:,:)
#if MPI_VER > 0
    real(RK), pointer, contiguous :: FSAll(:,:)
    real(RK), pointer, contiguous :: FBAll(:,:)
    real(RK), pointer, contiguous :: FRCAll(:,:)

! Components of the FTC Tensor(3)
    real(RK), pointer, contiguous :: FTC1All(:,:)
    real(RK), pointer, contiguous :: FTC2All(:,:)
    real(RK), pointer, contiguous :: FTC3All(:,:)

! Components of the FRC Tensor(3)
    real(RK), pointer, contiguous :: FRC1All(:,:)
    real(RK), pointer, contiguous :: FRC2All(:,:)
    real(RK), pointer, contiguous :: FRC3All(:,:)
#endif

#endif

    ! Total dipole moment of molecule for reaction field
    real(RK), pointer, contiguous :: MueX(:), MueY(:), MueZ(:)

    ! Torques from reaction field, space fixed
    real(RK), pointer, contiguous :: tRFX(:), tRFY(:), tRFZ(:)

    ! Total dipole moment of test particles for reaction field
    real(RK), pointer, contiguous :: MueXTest(:), MueYTest(:), MueZTest(:)

    ! Gear corrector local arrays
    real(RK), pointer, contiguous :: Corr0(:, :)
    real(RK), pointer, contiguous :: Corr1(:, :)

    ! Length of simulation box
    real(RK), pointer :: BoxLength

    ! Mole fraction of this component
    real(RK) :: Fraction

    ! Maximum number of particles in component
    integer, pointer :: NPartMax

    ! Number of particles in component
    integer, pointer :: NPart

    ! Number of particles in process
    ! Starting position, Number of Particles, Endposition
    integer, pointer :: NPart0, NPart1, NPart2

    ! Number of test particles
    integer, pointer :: NTest
    integer          :: NTestAll

    ! Number of degrees of freedom
    integer :: NDFTran, NDFRot, NDF

    ! Maximum allowed MC displacements
    real(RK), pointer :: DispTran, DispRot

    ! Number of MC attempts and successes
    integer :: NMoveAttempts, NMoveSuccesses
    integer :: NRotateAttempts, NRotateSuccesses
    integer :: NMoveBiasedAttempts, NMoveBiasedSuccesses
    integer :: NRotateBiasedAttempts, NRotateBiasedSuccesses

    ! Kinetic energy
    real(RK) :: EKinTran, EKinRot

    ! Chemical potential
    logical  :: CalcChemPot
    integer  :: ChemPotMethod, WFMethod, NGradThis
    integer  :: FluctState
    integer  :: GradInsInit
    real(RK) :: ChemPot, WidomContribution
    real (RK) :: HW_counter, HW_denom
!DEBUG
    real(RK) :: ChemPot1, ChemPot2
!DEBUG
    real(RK) :: ChemPot0, PartialMolarVolume
    real(RK) :: VarChemPot, VarPartialMolarVolume

    integer  :: BiasedPartners
    integer  :: BiasedPartnersNum

! Ewald
    real(RK) :: EPotTestSelf

    ! Fluctuating components and weighting factors
    integer           :: NFluctState, NFluctMax
    integer, pointer, contiguous  :: NState(:), NStateWF(:), NFluctComp(:)
    real(RK), pointer, contiguous :: WF(:)
    real(RK)          :: ProbW0, ProbW1, ProbW0V, ProbW1Rho
!DEBUG
    integer, pointer, contiguous  :: NFluctUpAttempts(:), NFluctUpSuccesses(:)
    integer, pointer, contiguous  :: NFluctDownAttempts(:), NFluctDownSuccesses(:)
!     integer, pointer  :: NStateBF(:)
!     real(RK), pointer :: BFSumState(:)
!DEBUG

    ! Variables for Thermodynamic Integration
    integer           :: NBins
    integer, pointer, contiguous  :: BinsVisit(:)
    real(RK)          :: Lambda, LambdaExponent, LaMin, LaMax, deltaLa, LaStepMax, ExpMinusBetaEnLaMin
    real(RK), pointer, contiguous :: BinsEn(:), BinsdEndLa(:), BinsIntdEndLa(:), BinsdEndLaV(:), BinsdEndLaH(:), BinsIntVW(:), BinsIntHW(:)

    ! Mole fraction in corresponding liquid simulation (for GE ensemble only)
    real(RK) :: LiqFraction

    ! Long-range corrections
    real(RK) :: EPotTestCorrMIE
    real(RK) :: EPotTestCorrTT68
    real(RK) :: EPotTestCorrRF

    ! Accumulated sums, averages and errors
    type(TAccumulator) :: SumInvChemPotRho
    type(TAccumulator) :: SumInvChemPot
    type(TAccumulator) :: SumChemPotV
    type(TAccumulator) :: SumChemPotVV
    type(TAccumulator) :: SumHW_counter
    type(TAccumulator) :: SumHW_denom
    type(TAccumulator) :: SumVW
    type(TAccumulator) :: SumHM
    type(TAccumulator) :: SumFraction
    type(TAccumulator) :: SumChemPotThermoIntWidom
    type(TAccumulator) :: SumChemPotThermoIntWidomV
#if OSMOP > 0
    type(TAccumulator),pointer, contiguous :: SumDenProfile(:)
#if OSMOP == 2
    type(TAccumulator),pointer, contiguous :: SumChemPotProfile(:)
#endif
#endif

    ! Potential model for this component
    type(TMolecule) :: Molecule

    ! File name for potential model
    character(FileNameLength) :: PotModFileName

  end type TComponent

  interface Construct
    module procedure TComponent_Construct
    module procedure TComponent_ConstructSVC
    module procedure TComponent_ConstructFluct
    module procedure TComponent_ConstructThermoInt
  end interface

  interface Destruct
    module procedure TComponent_Destruct
  end interface

  interface DestructFluct
    module procedure TComponent_DestructFluct
  end interface

  interface Allocate
    module procedure TComponent_Allocate
  end interface

  interface Deallocate
    module procedure TComponent_Deallocate
  end interface

  interface CreateAccumulators
    module procedure TComponent_CreateAccumulators
  end interface

  interface DestroyAccumulators
    module procedure TComponent_DestroyAccumulators
  end interface

  interface LongRangeCheck
    module procedure TComponent_LongRangeCheck
  end interface

  interface InitVelocities
    module procedure TComponent_InitVelocities
  end interface

  interface InitIntegratorGear
    module procedure TComponent_InitIntegratorGear
  end interface

  interface InitIntegratorLeapFrog
    module procedure TComponent_InitIntegratorLeap
  end interface

  interface InitIntegratorVerlet
    module procedure TComponent_InitIntegratorVerlet
  end interface

  interface InitIntegratorVV
    module procedure TComponent_InitIntegratorVV
  end interface

  interface RemoveNetMomentum
    module procedure TComponent_RemoveNetMomentum
  end interface

  interface CalculateEKin
    module procedure TComponent_CalculateEKin
  end interface

  interface Mol2Atom
    module procedure TComponent_Mol2Atom
  end interface

  interface Mol2Atom1
    module procedure TComponent_Mol2Atom1
  end interface

  interface Mol2AtomTest
    module procedure TComponent_Mol2AtomTest
  end interface

  interface Atom2Mol
    module procedure TComponent_Atom2Mol
  end interface

  interface Atom2Mol1
    module procedure TComponent_Atom2Mol1
  end interface

  interface Atom2Mol_Trans
    module procedure TComponent_Atom2Mol_Trans
  end interface

#if OSMOP > 0
  interface DensityProfile
    module procedure TComponent_DensityProfile
  end interface
#endif

  interface PredictGear
    module procedure TComponent_PredictGear
  end interface

  interface CorrectGear
    module procedure TComponent_CorrectGear
  end interface

  interface PredictLeapFrog
    module procedure TComponent_PredictLeapFrog
  end interface

  interface CorrectLeapFrog
    module procedure TComponent_CorrectLeapFrog
  end interface

  interface PredictVerlet
    module procedure TComponent_PredictVerlet
  end interface

  interface CorrectVerlet
    module procedure TComponent_CorrectVerlet
  end interface

  interface PredictVV
    module procedure TComponent_PredictVV
  end interface

  interface CorrectVV
    module procedure TComponent_CorrectVV
  end interface

  interface ZeroNAttempts
    module procedure TComponent_ZeroNAttempts
  end interface

  interface UpdateDisplacements
    module procedure TComponent_UpdateDisplacements
  end interface

  interface AddParticle
    module procedure TComponent_AddParticle
  end interface

  interface RemoveParticle
    module procedure TComponent_RemoveParticle
  end interface

  interface UpdateChemPot
    module procedure TComponent_UpdateChemPot
  end interface

  interface SaveState
    module procedure TComponent_SaveState
  end interface

  interface RestoreState
    module procedure TComponent_RestoreState
  end interface

  interface RestartSave
    module procedure TComponent_RestartSave
  end interface

  interface RestartRead
    module procedure TComponent_RestartRead
  end interface

#if CONSTR > 0
  interface CorrectGear_Constraint
    module procedure TComponent_CorrectGear_Constraint
  end interface
#endif

#if  TRANS == 1
!TRANSPORT_start
  interface ForceTransport
    module procedure TComponent_ForceTransport
  end interface
!TRANSPORT_END
#endif

contains



!==============================================================!
!  Subroutine TComponent_Construct                             !
!==============================================================!

  subroutine TComponent_Construct( this, comp )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent) :: this
    integer          :: comp

    ! Declare local variables
    character( IOBufferLength ) :: str
    integer                     :: stat

    ! Allocate number of particles in component
    allocate( this%NPart, STAT = stat )
    call AllocationError( stat, 'number of particles' )

    ! Allocate number of particles in process
    allocate( this%NPart0, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NPart1, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NPart2, STAT = stat )
    call AllocationError( stat, 'number of particles' )

    ! Allocate number of test particles
    allocate( this%NTest, STAT = stat )
    call AllocationError( stat, 'number of test particles' )

    ! Read file name for potential model
    call FileReadParameter( this%PotModFileName, iounit_params , IdPotModFileName, .false. )

    ! Read mole fraction of this component
    write( IOBuffer, '(72(1H-))')
    call LogWrite
    write( IOBuffer, '(T13, "Reading component", I3," for ensemble")') comp
    call LogWrite
    call FileReadParameter( this%Fraction, iounit_params , IdFraction, .false. )
    write( IOBuffer, '("Mole fraction of component ", A, ": ", F9.6)' ) trim( this%PotModFileName ), this%Fraction
    call LogWrite

#if TRANS==1
 ! Read partial molar enthalpy from the paremeters file
    call FileReadParameter( this%PartialMolarEnthalpy, iounit_params , IdPartialMolarEnthalpy, .false., 0._RK )

    if (this%PartialMolarEnthalpy .ne. 0._RK) then
      write( IOBuffer,'("Reduced PartMolEnt of component ", A, ": ", F12.8 )' ) &
&       trim( this%PotModFilename ), this%PartialMolarEnthalpy
      call LogWrite
    end if

#endif


    ! Initialize flag for calculation of chemical potential
    this%CalcChemPot = .false.
    this%NTest = 0

    ! Initialize fluctuating state (for GradIns)
    this%FluctState = -1



    if( EnsembleType .eq. EnsembleTypeGE ) then
      ! Read mole fraction of liquid simulation
      call FileReadParameter( this%LiqFraction, iounit_params , IdLiqFraction, .false. )

      ! Read chemical potential and partial molar volume and their
      ! uncertainties for Grand Equilibrium
      call FileReadParameter( this%ChemPot0, iounit_params , IdChemPot, .false. )
      call FileReadParameter( this%VarChemPot, iounit_params , IdVarChemPot, .false. )
      call FileReadParameter( this%PartialMolarVolume, iounit_params , IdPartialMolarVolume, .false. )
      call FileReadParameter( this%VarPartialMolarVolume, iounit_params , IdVarPartialMolarVolume, .false. )
      write( IOBuffer,'("Reduced ChemPot0 of component ", A, ": ", F9.6, " (", F9.6, ")")' ) &
&       trim( this%PotModFileName ), this%ChemPot0, this%VarChemPot
      call LogWrite
      write( IOBuffer,'("Reduced PartMolVol of component ", A, ": ", F9.6, " (", F9.6, ")")' ) &
&       trim( this%PotModFilename ), this%PartialMolarVolume, this%VarPartialMolarVolume
      call LogWrite

    else if( EnsembleType .eq. EnsembleTypeHA ) then
      if( comp == 1 ) then
        ! Read chemical potential of phase changing component (first one)
        call FileReadParameter( this%ChemPot, iounit_params , IdChemPot, .false. )
        call FileReadParameter( this%VarChemPot, iounit_params , IdVarChemPot, .false. )
        write( IOBuffer, '("Reduced ChemPot of component ", A, ": ", F9.6, " (", F9.6, ")")' ) &
&         trim( this%PotModFileName ), this%ChemPot0, this%VarChemPot
        call LogWrite
      end if

    else
      ! Read method for calculation of chemical potential
      call FileReadParameter( str, iounit_params, IdChemPotMethod, .false., "NONE" )
      select case( str )
      case( 'NONE', 'None', 'none' )
        this%ChemPotMethod = ChemPotMethodNone
        str = 'no calculation'
      case( 'WIDOM', 'Widom', 'widom' )
        this%ChemPotMethod = ChemPotMethodWidom
        str = 'Widom''s test particle method'
      case( 'GRADINS', 'GradIns', 'Gradins', 'gradins' )
        this%ChemPotMethod = ChemPotMethodGradIns
        this%FluctState = 0
        str = 'gradual insertion'
      case( 'THERMOINT', 'ThermoInt', 'Thermoint', 'thermoint' )
        this%ChemPotMethod = ChemPotMethodThermoInt
        str = 'thermodynamic integration'
      case default
        call Error( trim( str )//  ' method for calculation of chemical potential is not implemented' )
      end select
      if( this%ChemPotMethod .eq. ChemPotMethodGradIns .and. .not. SimulationType .eq. MonteCarlo ) &
&       call Error( 'Gradual insertion is only allowed for MonteCarlo simulation' )
      write( IOBuffer, '("Chemical potential of ", A, " will be calculated by: ", A)' ) &
&       trim( this%PotModFilename )
      call LogWrite
      write( IOBuffer, '(T10, "-> ", A)' ) trim( str )
      call LogWrite

      ! Read Gradual Insertion Initialization Steps
      if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
        call FileReadParameter( this%GradInsInit, iounit_params , IdGradInsInit, .false., 0 )
        write( IOBuffer, '("Grad. Ins. initialization Steps: ", T40, I7)' ) this%GradInsInit
        call LogWrite
      end if

      ! Read number of test particles
      if( this%ChemPotMethod .eq. ChemPotMethodWidom ) then
        call FileReadParameter( this%NTest, iounit_params, IdNTest, .false. )
        if( this%NTest <= 0 ) call Error( 'Number of test particles need to be > 0' )
        write( IOBuffer, '(T10, "-> Number of test particles:", I11 )' ) this%NTest

#if MPI_VER>0
        if (SimulationType .eq. MolecularDynamics) then
           this%NTest = ((this%NTest -1)/NProcs +1)
        endif
#endif
      end if

      ! Read weighting factors method
      this%WFMethod = WFMethodNone
      if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
        call FileReadParameter( str, iounit_params, IdWeightFactors, .false. )
        select case(str)
        case( 'auto', 'Auto' )
          call Error( 'Method "auto" for weighting factors is not implemented' )
          str = 'automatic method'
        case( 'guess', 'Guess' )
          this%WFMethod = WFMethodGuess
          str = 'first guess'
        case( 'optset', 'Optset', 'OptSet' )
          this%WFMethod = WFMethodOptSet
          str = 'optimized set'
        case default
          call Error( trim( str )// ' method for weighting factors is not implemented' )
        end select
        write( IOBuffer, '("Estimation of weighting factors: using ", A )' ) trim( str )
        call LogWrite
      end if
      if (this%ChemPotMethod .eq. ChemPotMethodThermoInt ) then
        call FileReadParameter( this%LaMin, iounit_params , IdLambdaMin, .false., 0.2_RK )
        write( IOBuffer, '("Thermo. Int. LambdaMin: ", T40, F7.5)' ) this%LaMin
        call LogWrite
        call FileReadParameter( this%LaMax, iounit_params , IdLambdaMax, .false., 1.0_RK )
        write( IOBuffer, '("Thermo. Int. LambdaMax: ", T40, F7.5)' ) this%LaMax
        call LogWrite
        call FileReadParameter( this%NBins, iounit_params , IdNBins, .false., 100 )
        write( IOBuffer, '("Thermo. Int. NBins: ", T40, I7)' ) this%NBins
        call LogWrite
        if (SimulationType .eq. MolecularDynamics) then
          call FileReadParameter( this%LaStepMax, iounit_params , IdLambdaStepMax, .false., 0.01_RK)
        else
          call FileReadParameter( this%LaStepMax, iounit_params , IdLambdaStepMax, .false., 0.1_RK)
        end if
        write( IOBuffer, '("Thermo. Int. LambdaStepMax: ", T40, F7.5)' ) this%LaStepMax
        call LogWrite
        call FileReadParameter( this%LambdaExponent, iounit_params , IdLambdaExponent, .false., 4.0_RK)
        write( IOBuffer, '("Thermo. Int. LambdaExponent: ", T40, F7.5)' ) this%LambdaExponent
        call LogWrite
        call FileReadParameter( this%NTest, iounit_params, IdNTest, .false., 100 )
        write( IOBuffer, '(T10, "-> Number of test particles:", I11 )' ) this%NTest
        call LogWrite
#if MPI_VER>0
        if (SimulationType .eq. MolecularDynamics) then
          this%NTest = ((this%NTest-1)/NProcs +1)
        endif
#endif
        if (this%LaMin**this%LambdaExponent .lt. 1E-30_RK) then
          this%LaMin = 1E-30_RK**(1._RK/this%LambdaExponent)
          write( IOBuffer, '("LambdaMin too low for simulation and was changed!")')
          call LogWrite
        endif
        this%deltaLa=(this%LaMax-this%LaMin)/this%NBins
      end if

    end if

#if OSMOP > 0
    call FileReadParameter( str, iounit_params, IdPermeability, .false.)
    select case( str )
      case( 'OFF', 'Off', 'off', 'NO', 'No', 'no', 'false', 'False', 'FALSE' )
        this%permeable = .false.
        write( IOBuffer, '("component ", A, " is not permeable.")' ) trim( this%PotModFilename )
        call LogWrite
      case('YES', 'Yes', 'yes', 'ON', 'On', 'on', 'right', 'Right', 'RIGHT' )
        this%permeable = .true.
        write( IOBuffer, '("component ", A, " is permeable.")' ) trim( this%PotModFilename )
        call LogWrite
      case default
        call Error( trim( str )//  '  unknown. Set value to Yes or On.' )
      end select
#endif

    ! Create potential model
    call Construct( this%Molecule, this%PotModFileName, &
&     merge(0, -1, this%ChemPotMethod .eq. ChemPotMethodGradIns) )
    this%NFluctMax = this%Molecule%NFluct
    this%NFluctState = 0

    ! Allocate maximum allowed MC displacements
    if( (SimulationType .eq. MonteCarlo) .or. (SimulationType .eq. Gibbs) .or. MCOverlapReduction ) then
      allocate( this%DispTran, STAT = stat )
      call AllocationError( stat, 'maximum MC displacement' )
      allocate( this%DispRot, STAT = stat )
      call AllocationError( stat, 'maximum MC displacement' )
    end if

    ! Allocate and read weighting factors
    if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
      allocate( this%WF( 0:this%NFluctMax ), STAT = stat )
      call AllocationError( stat, 'fluctuating particle states', &
&       this%NFluctMax + 1 )
      if( this%WFMethod .eq. WFMethodGuess .or. this%WFMethod .eq. WFMethodOptSet ) then
        if( RootProc ) read( iounit_params, * ) this%WF
#if MPI_VER > 0
        call MPI_Bcast( this%WF, size( this%WF ), MPI_RK, &
&         NRootProc, Communicator, ierror )
#endif
      end if
    end if

    ! Allocate fluctuating particle components vector
    nullify( this%NFluctComp )
    if( this%NFluctMax > 0 ) then
      allocate( this%NFluctComp( 0:this%NFluctMax ), STAT = stat )
      call AllocationError( stat, 'fluctuating particle components', this%NFluctMax + 1 )
    end if

    write( IOBuffer, '(T8, "Reading component", I3," for ensemble successful")') comp
    call LogWrite
    write( IOBuffer, '(72(1H-))')
    call LogWrite

  end subroutine TComponent_Construct


!==============================================================!
!  Subroutine TComponent_ConstructSVC                          !
!==============================================================!

  subroutine TComponent_ConstructSVC( this, PotModFileName )

    implicit none

    ! Declare arguments
    type(TComponent)                      :: this
    character(FileNameLength), intent(in) :: PotModFileName

    ! Declare local variables
    integer :: stat

    ! Allocate number of particles in component
    allocate( this%NPart, STAT = stat )
    call AllocationError( stat, 'number of particles' )

    ! Allocate number of particles in process
    allocate( this%NPart0, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NPart1, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NPart2, STAT = stat )
    call AllocationError( stat, 'number of particles' )

    ! Allocate number of test particles
    allocate( this%NTest, STAT = stat )
    call AllocationError( stat, 'number of test particles' )
    this%NTest = 0

    ! Set file name for potential model
    this%PotModFileName = PotModFileName

    ! Set mole fraction of this component
    this%Fraction = 1._RK

    ! Create potential model
    call Construct( this%Molecule, this%PotModFileName, -1 )

  end subroutine TComponent_ConstructSVC

!==============================================================!
!  Subroutine TComponent_ConstructThermoInt                    !
!==============================================================!

  subroutine TComponent_ConstructThermoInt( this, comp0 )

    implicit none

    ! Declare arguments
    type(TComponent)             :: this
    type(TComponent), intent(in) :: comp0

    ! Declare local variables
    integer :: stat

    ! Allocate number of particles in component
    allocate( this%NPart, STAT = stat )
    call AllocationError( stat, 'number of particles' )

    ! Allocate number of particles in process
    allocate( this%NPart0, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NPart1, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NPart2, STAT = stat )
    call AllocationError( stat, 'number of particles' )

    ! Allocate number of particles in component
    allocate( this%NTest, STAT = stat )
    call AllocationError( stat, 'number of particles' )

    ! Copy file name for potential model
    this%PotModFileName = comp0%PotModFileName

    ! Set number of particles and mole fraction of this component
    this%NPart = 0
    this%NTest = 0
    this%Fraction = 0._RK
    this%NBins = 0

    this%Lambda = 1.0_RK - Zero

    ! Set fluctuating state (for GradIns)
    this%FluctState = 0
    this%ChemPotMethod = ChemPotMethodNone

    ! Set chemical potential flag
    this%CalcChemPot = .false.

    ! Create potential model
    call Construct( this%Molecule, this%PotModFileName, -1 )

    ! Set maximum allowed MC displacements
    this%DispTran => comp0%DispTran
    this%DispRot => comp0%DispRot

  end subroutine TComponent_ConstructThermoInt

!==============================================================!
!  Subroutine TComponent_ConstructFluct                        !
!==============================================================!

  subroutine TComponent_ConstructFluct( this, comp0, state )

    implicit none

    ! Declare arguments
    type(TComponent)             :: this
    type(TComponent), intent(in) :: comp0
    integer, intent(in)          :: state

    ! Declare local variables
    integer :: stat

    ! Allocate number of particles in component
    allocate( this%NPart, STAT = stat )
    call AllocationError( stat, 'number of particles' )

    ! Allocate number of particles in process
    allocate( this%NPart0, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NPart1, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NPart2, STAT = stat )
    call AllocationError( stat, 'number of particles' )

    ! Allocate number of particles in component
    allocate( this%NTest, STAT = stat )
    call AllocationError( stat, 'number of particles' )

    ! Copy file name for potential model
    this%PotModFileName = comp0%PotModFileName

    ! Set number of particles and mole fraction of this component
    this%NPart = 0
    this%NTest = 0
    this%Fraction = 0._RK

    ! Set fluctuating state (for GradIns)
    this%FluctState = state
    this%ChemPotMethod = ChemPotMethodNone

    ! Set chemical potential flag
    this%CalcChemPot = .false.

    ! Create potential model
    call Construct( this%Molecule, this%PotModFileName, state )

    ! Set maximum allowed MC displacements
    this%DispTran => comp0%DispTran
    this%DispRot => comp0%DispRot

  end subroutine TComponent_ConstructFluct



!==============================================================!
!  Subroutine TComponent_Destruct                              !
!==============================================================!

  subroutine TComponent_Destruct( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Destroy potential model
    call Destruct( this%Molecule )

    ! Deallocate number of test particles
    if( associated( this%NTest ) ) then
      deallocate( this%NTest )
    end if

    ! Deallocate number of particles in process
    if( associated( this%NPart0 ) ) then
      deallocate( this%NPart0 )
    end if
    if( associated( this%NPart1 ) ) then
      deallocate( this%NPart1 )
    end if
    if( associated( this%NPart2 ) ) then
      deallocate( this%NPart2 )
    end if

    ! Deallocate number of particles in component
    if( associated( this%NPart ) ) then
      deallocate( this%NPart )
    end if

    if( SimulationType .eq. MonteCarlo .or. MCOverlapReduction ) then
      ! Deallocate maximum allowed MC displacements
      if( associated( this%DispTran ) ) then
        deallocate( this%DispTran )
      end if
      if( associated( this%DispRot ) ) then
        deallocate( this%DispRot )
      end if
    end if

  end subroutine TComponent_Destruct


!==============================================================!
!  Subroutine TComponent_DestructFluct                         !
!==============================================================!

  subroutine TComponent_DestructFluct( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Destroy potential model
    call Destruct( this%Molecule )

    ! Deallocate number of test particles
    if( associated( this%NTest ) ) then
      deallocate( this%NTest )
    end if

    ! Deallocate number of particles in process
    if( associated( this%NPart0 ) ) then
      deallocate( this%NPart0 )
    end if
    if( associated( this%NPart1 ) ) then
      deallocate( this%NPart1 )
    end if
    if( associated( this%NPart2 ) ) then
      deallocate( this%NPart2 )
    end if

    ! Deallocate number of particles in component
    if( associated( this%NPart ) ) then
      deallocate( this%NPart )
    end if

  end subroutine TComponent_DestructFluct


!==============================================================!
!  Subroutine TComponent_CreateAccumulators                    !
!==============================================================!

  subroutine TComponent_CreateAccumulators( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

#if OSMOP > 0
    ! Declare local variables
    integer          :: i, stat

    allocate( this%SumDenProfile(NBinsDen ), STAT = stat )
    call AllocationError( stat, 'NBinsDen', NBinsDen )
#if OSMOP == 2
    allocate( this%SumChemPotProfile(NBinsDen ), STAT = stat )
    call AllocationError( stat, 'NBinsDen', NBinsDen )
#endif
#endif

    ! Construct accumulators
    select case( this%ChemPotMethod )
    case( ChemPotMethodGradIns )
      call Construct( this%SumInvChemPotRho, .true. )
      call Construct( this%SumInvChemPot, .true. )
      call Construct( this%SumVW, .true. )
      call Construct( this%SumHM, .true. )
    case( ChemPotMethodWidom )
      call Construct( this%SumChemPotV, .false. )
      call Construct( this%SumChemPotVV, .false. )
      call Construct( this%SumHW_counter, .false. )
      call Construct( this%SumHW_denom, .false. )
      call Construct( this%SumVW, .true. )
      call Construct( this%SumHM, .true. )
    case( ChemPotMethodThermoInt )
      call Construct( this%SumChemPotV, .true. )
      call Construct( this%SumChemPotVV, .false. )
      call Construct( this%SumChemPotThermoIntWidom, .false. )
      call Construct( this%SumChemPotThermoIntWidomV, .false. )
      call Construct( this%SumHW_counter, .false. )
      call Construct( this%SumHW_denom, .false. )
      call Construct( this%SumVW, .true. )
      call Construct( this%SumHM, .true. )
    end select

    if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
      call Construct( this%SumFraction, .false. )
    end if

#if OSMOP > 0
    do i = 1, NBinsDen
      call Construct( this%SumDenProfile(i), .false. )
#if OSMOP == 2
      if ( this%ChemPotMethod .eq. ChemPotMethodWidom ) call Construct( this%SumChemPotProfile(i), .false. )
#endif
    end do
#endif

  end subroutine TComponent_CreateAccumulators



!==============================================================!
!  Subroutine TComponent_DestroyAccumulators                   !
!==============================================================!

  subroutine TComponent_DestroyAccumulators( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

#if OSMOP > 0
    ! Declare local variables
    integer          :: i
#endif

    ! Destruct accumulators
    select case( this%ChemPotMethod )
    case( ChemPotMethodGradIns )
      call Destruct( this%SumInvChemPotRho )
      call Destruct( this%SumInvChemPot )
      call Destruct( this%SumVW )
      call Destruct( this%SumHM )
    case( ChemPotMethodWidom )
      call Destruct( this%SumChemPotV )
      call Destruct( this%SumChemPotVV )
      call Destruct( this%SumHW_counter )
      call Destruct( this%SumHW_denom )
      call Destruct( this%SumVW )
      call Destruct( this%SumHM )
    case( ChemPotMethodThermoInt )
      call Destruct( this%SumChemPotV )
      call Destruct( this%SumChemPotVV )
      call Destruct( this%SumChemPotThermoIntWidom )
      call Destruct( this%SumChemPotThermoIntWidomV )
      call Destruct( this%SumHW_counter )
      call Destruct( this%SumHW_denom )
      call Destruct( this%SumVW )
      call Destruct( this%SumHM )
    end select

    if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
      call Destruct( this%SumFraction )
    end if

#if OSMOP > 0
    do i = 1, NBinsDen
      call Destruct( this%SumDenProfile(i) )
#if OSMOP == 2
      if ( this%ChemPotMethod .eq. ChemPotMethodWidom ) call Destruct( this%SumChemPotProfile(i) )
#endif
    end do
    if( associated( this%SumDenProfile ) ) then
      deallocate( this%SumDenProfile )
    end if
#if OSMOP == 2
    if( associated( this%SumChemPotProfile ) ) then
      deallocate( this%SumChemPotProfile )
    end if
#endif
#endif

  end subroutine TComponent_DestroyAccumulators

!==============================================================!
!  Subroutine TComponent_Allocate                              !
!==============================================================!

  subroutine TComponent_Allocate( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Declare local variables
    integer :: np, ntest, nf
    integer :: i
    integer :: stat

    ! Set maximum number of particles and number of test particles
    np = this%NPartMax
    ntest = this%NTest

    ! Nullify pointers
    nullify( this%P0 )
    nullify( this%P0Save )
    nullify( this%P0old )
    nullify( this%P1 )
    nullify( this%P2 )
    nullify( this%P3 )
    nullify( this%P4 )
    nullify( this%P5 )
    nullify( this%Disp )
    nullify( this%ri0_x )
    nullify( this%ri0_y )
    nullify( this%ri0_z )
    nullify( this%F )
#if MPI_VER > 0
    nullify( this%FAll )
    nullify( this%NAdd )
#endif
    nullify( this%Q0 )
    nullify( this%Q0Save )
    nullify( this%Q0tmp )
    nullify( this%Q1 )
    nullify( this%Q2 )
    nullify( this%Q3 )
    nullify( this%Q4 )
    nullify( this%W0 )
    nullify( this%W1 )
    nullify( this%W2 )
    nullify( this%W3 )
    nullify( this%W4 )
    nullify( this%T )
#if MPI_VER > 0
    nullify( this%TAll )
#endif
    nullify( this%MueX )
    nullify( this%MueY )
    nullify( this%MueZ )
    nullify( this%tRFX )
    nullify( this%tRFY )
    nullify( this%tRFZ )
    nullify( this%MueXTest )
    nullify( this%MueYTest )
    nullify( this%MueZTest )
    nullify( this%Corr0 )
    nullify( this%Corr1 )
    nullify( this%NState )
    nullify( this%NStateWF )
#if  TRANS == 1
!  Transport  !TRANSPORT_start
    nullify(this%KinETran)
    nullify(this%KinEPart)
    nullify( this%FS )
    nullify( this%FB )
    nullify( this%FRC )
    nullify( this%FTC )
    nullify( this%FTC1)
    nullify( this%FTC2 )
    nullify( this%FTC3 )
    nullify( this%FRC1)
    nullify( this%FRC2 )
    nullify( this%FRC3 )

    nullify( this%ri0_E_x ) !EinsteinCoef ri0_E nullify
    nullify( this%ri0_E_y )
    nullify( this%ri0_E_z )

#if OSMOP > 0
    nullify( this%FOsmoticPressure )
    nullify( this%DensityProfileN )
#if OSMOP == 2
    nullify( this%ChemPotProfile )
#endif
#endif

#if MPI_VER > 0
    nullify( this%FSAll )
    nullify( this%FBAll )
    nullify( this%FRCAll )

    nullify( this%FTC1All )
    nullify( this%FTC2All )
    nullify( this%FTC3All )

    nullify( this%FRC1All )
    nullify( this%FRC2All )
    nullify( this%FRC3All )
#endif

    nullify( this%BinsVisit )
    nullify( this%BinsEn )
    nullify( this%BinsdEndLa )
    nullify( this%BinsIntdEndLa )
    nullify( this%BinsdEndLaV )
    nullify( this%BinsdEndLaH )
    nullify( this%BinsIntVW )
    nullify( this%BinsIntHW )

    ! Transport
    allocate( this%KinETran( np, 3 ), STAT = stat )
    call AllocationError( stat, '3*particles', np )
    allocate( this%KinEPart( np), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%FS( np, 3 ), STAT = stat )
    call AllocationError( stat, '3*particles', np )
    allocate( this%FB( np, 3 ), STAT = stat )
    call AllocationError( stat, '3*particles', np )
    allocate( this%FTC( np, 3 ), STAT = stat )
    call AllocationError( stat, '3*particles', np )
    allocate( this%FRC( np, 3 ), STAT = stat )
    call AllocationError( stat, '3*particles', np )
    allocate( this%FTC1( np, 3 ), STAT = stat )
    call AllocationError( stat, '3*particles', np )
    allocate( this%FTC2( np, 3 ), STAT = stat )
    call AllocationError( stat, '3*particles', np )
    allocate( this%FTC3( np, 3 ), STAT = stat )
    call AllocationError( stat, '3*particles', np )
    allocate( this%FRC1( np, 3 ), STAT = stat )
    call AllocationError( stat, '3*particles', np )
    allocate( this%FRC2( np, 3 ), STAT = stat )
    call AllocationError( stat, '3*particles', np )
    allocate( this%FRC3( np, 3 ), STAT = stat )
    call AllocationError( stat, '3*particles', np )
    allocate( this%Q0( np, 4 ), STAT = stat )
    call AllocationError( stat, '4*particles', np )

#if MPI_VER > 0
    allocate( this%FSAll( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%FBAll( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%FRCAll( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%FTC1All( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%FTC2All( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%FTC3All( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%FRC1All( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%FRC2All( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%FRC3All( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

#endif

    this%FS(: , :)   = 0._RK
    this%FB(: , :)   = 0._RK
    this%FTC(: , :)  = 0._RK
    this%FRC(: , :)  = 0._RK

    this%FTC1(:,:)  = 0._RK
    this%FTC2(:,:)  = 0._RK
    this%FTC3(:,:)  = 0._RK

    this%FRC1(:,:)  = 0._RK
    this%FRC2(:,:)  = 0._RK
    this%FRC3(:,:)  = 0._RK
!TRANSPORT_END
#endif

    ! Centers of mass positions
    allocate( this%P0( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%P0Save( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%P0old( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    if( SimulationType .eq. MolecularDynamics ) then

      ! Centers of mass positions' derivatives
      allocate( this%P1( np, 3 ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%P2( np, 3 ), STAT = stat )
      call AllocationError( stat, 'particles', np )

      if( IntegratorType .eq. IntegratorTypeGear ) then
        allocate( this%P3( np, 3 ), STAT = stat )
        call AllocationError( stat, 'particles', np )
        allocate( this%P4( np, 3 ), STAT = stat )
        call AllocationError( stat, 'particles', np )
        allocate( this%P5( np, 3 ), STAT = stat )
        call AllocationError( stat, 'particles', np )
      end if

      ! Displacement
      allocate( this%Disp( np, 3 ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      this%Disp(:, :) = 0._RK

      if( ALPHA2UpdateFrequency > 0 ) then
          ! Alpha2
          allocate( this%ri0_x( np, 0:ALPHA2Length/ALPHA2Shift-1 ), STAT = stat )
          call AllocationError( stat, 'particles', np )
          allocate( this%ri0_y( np, 0:ALPHA2Length/ALPHA2Shift-1 ), STAT = stat )
          call AllocationError( stat, 'particles', np )
          allocate( this%ri0_z( np, 0:ALPHA2Length/ALPHA2Shift-1 ), STAT = stat )
          call AllocationError( stat, 'particles', np )
          this%ri0_x(:, :) = 0._RK
          this%ri0_y(:, :) = 0._RK
          this%ri0_z(:, :) = 0._RK
      end if

#if TRANS==1
      !EinsteinCoef allocate ri0_E
      if( EinsteinCoefCalc) then
          allocate( this%ri0_E_x( np, 0:this%NEinstein-1 ), STAT = stat )
          call AllocationError( stat, 'particles', np )
          allocate( this%ri0_E_y( np, 0:this%NEinstein-1 ), STAT = stat )
          call AllocationError( stat, 'particles', np )
          allocate( this%ri0_E_z( np, 0:this%NEinstein-1 ), STAT = stat )
          call AllocationError( stat, 'particles', np )
          this%ri0_E_x(:, :) = 0._RK
          this%ri0_E_y(:, :) = 0._RK
          this%ri0_E_z(:, :) = 0._RK
      end if
#endif

      ! Total forces
      allocate( this%F( np, 3 ), STAT = stat )
      call AllocationError( stat, 'particles', np )
#if MPI_VER > 0
      allocate( this%FAll( np, 3 ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%NAdd( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
#endif

#if OSMOP > 0
      ! Force for osmotic pressure
      allocate( this%FOsmoticPressure( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      ! Density Profile
      allocate( this%DensityProfileN( NBinsDen ), STAT = stat )
      call AllocationError( stat, 'BinsDensity', NBinsDen )
#if OSMOP == 2
      allocate( this%ChemPotProfile( NBinsDen ), STAT = stat )
      call AllocationError( stat, 'BinsChemPot', NBinsDen )
#endif
#endif
    end if

    if( this%Molecule%isElongated ) then

#if  TRANS != 1
! For the calculation of transport properties, the necessary quaternion matrix has
! already been allocated in this subroutine!
      ! Quaternion parameters
      allocate( this%Q0( np, 4 ), STAT = stat )
      call AllocationError( stat, 'particles', np )
#endif
      allocate( this%Q0Save( np, 4 ), STAT = stat )
      call AllocationError( stat, 'particles', np )

      if( SimulationType .eq. MolecularDynamics ) then

        if( IntegratorType .eq. IntegratorTypeLeapFrog ) then
          allocate( this%Q0tmp( np, 4 ), STAT = stat )
          call AllocationError( stat, 'particles', np )
        end if

        ! Quaternion parameters' derivatives
        allocate( this%Q1( np, 4 ), STAT = stat )
        call AllocationError( stat, 'particles', np )

        if( IntegratorType .eq. IntegratorTypeGear ) then
          allocate( this%Q2( np, 4 ), STAT = stat )
          call AllocationError( stat, 'particles', np )
          allocate( this%Q3( np, 4 ), STAT = stat )
          call AllocationError( stat, 'particles', np )
          allocate( this%Q4( np, 4 ), STAT = stat )
          call AllocationError( stat, 'particles', np )
        end if

        ! Angular velocities and their derivatives
        allocate( this%W0( np, 3 ), STAT = stat )
        call AllocationError( stat, 'particles', np )
        allocate( this%W1( np, 3 ), STAT = stat )
        call AllocationError( stat, 'particles', np )

        if( IntegratorType .eq. IntegratorTypeGear ) then
          allocate( this%W2( np, 3 ), STAT = stat )
          call AllocationError( stat, 'particles', np )
          allocate( this%W3( np, 3 ), STAT = stat )
          call AllocationError( stat, 'particles', np )
          allocate( this%W4( np, 3 ), STAT = stat )
          call AllocationError( stat, 'particles', np )
        end if

        ! Total torques
        allocate( this%T( np, 3 ), STAT = stat )
        call AllocationError( stat, 'particles', np )
#if MPI_VER > 0
        allocate( this%TAll( np, 3 ), STAT = stat )
        call AllocationError( stat, 'particles', np )
#endif
        ! Torques from reaction field
        allocate( this%tRFX( np ), STAT = stat )
        call AllocationError( stat, 'particles', np )
        allocate( this%tRFY( np ), STAT = stat )
        call AllocationError( stat, 'particles', np )
        allocate( this%tRFZ( np ), STAT = stat )
        call AllocationError( stat, 'particles', np )
      end if

      ! Total dipole moment of molecules for reaction field
      if( CutoffMode .eq. CenterofMass ) then
        allocate( this%MueX( np ), STAT = stat )
        call AllocationError( stat, 'particles', np )
        allocate( this%MueY( np ), STAT = stat )
        call AllocationError( stat, 'particles', np )
        allocate( this%MueZ( np ), STAT = stat )
        call AllocationError( stat, 'particles', np )

        if( ntest > 0 ) then
          allocate( this%MueXTest( ntest ), STAT = stat )
          call AllocationError( stat, 'particles', ntest )
          allocate( this%MueYTest( ntest ), STAT = stat )
          call AllocationError( stat, 'particles', ntest )
          allocate( this%MueZTest( ntest ), STAT = stat )
          call AllocationError( stat, 'particles', ntest )
        end if
      end if

    end if

    ! Gear corrector local arrays
    if( SimulationType .eq. MolecularDynamics .and. IntegratorType .eq. IntegratorTypeGear ) then
      allocate( this%Corr0( np, merge( 4, 3, this%Molecule%isElongated ) ),STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%Corr1( np, merge( 4, 3, this%Molecule%isElongated ) ),STAT = stat )
      call AllocationError( stat, 'particles', np )
    end if

    ! Site positions, orientations, forces and torques
    do i = 1, this%Molecule%NMIEnm
      this%Molecule%SiteMIEnm(i)%NPartMax => this%NPartMax
      this%Molecule%SiteMIEnm(i)%NPart => this%NPart
      this%Molecule%SiteMIEnm(i)%NTest => this%NTest
      this%Molecule%SiteMIEnm(i)%NPart0 => this%NPart0
      this%Molecule%SiteMIEnm(i)%NPart1 => this%NPart1
      this%Molecule%SiteMIEnm(i)%NPart2 => this%NPart2

      call Allocate( this%Molecule%SiteMIEnm(i) )
      this%Molecule%SiteMIEnm(i)%PX => this%P0(:, 1)
      this%Molecule%SiteMIEnm(i)%PY => this%P0(:, 2)
      this%Molecule%SiteMIEnm(i)%PZ => this%P0(:, 3)

      if( ntest > 0 ) then
        this%Molecule%SiteMIEnm(i)%PXTest => this%P0Test(:, 1)
        this%Molecule%SiteMIEnm(i)%PYTest => this%P0Test(:, 2)
        this%Molecule%SiteMIEnm(i)%PZTest => this%P0Test(:, 3)
      end if

#if TRANS==1
      this%Molecule%SiteMIEnm(i)%Q0r => this%Q0
#endif
    end do

    do i = 1, this%Molecule%NTT68
      this%Molecule%SiteTT68(i)%NPartMax => this%NPartMax
      this%Molecule%SiteTT68(i)%NPart => this%NPart
      this%Molecule%SiteTT68(i)%NTest => this%NTest
      this%Molecule%SiteTT68(i)%NPart0 => this%NPart0
      this%Molecule%SiteTT68(i)%NPart1 => this%NPart1
      this%Molecule%SiteTT68(i)%NPart2 => this%NPart2

      call Allocate( this%Molecule%SiteTT68(i) )
      this%Molecule%SiteTT68(i)%PX => this%P0(:, 1)
      this%Molecule%SiteTT68(i)%PY => this%P0(:, 2)
      this%Molecule%SiteTT68(i)%PZ => this%P0(:, 3)

      if( ntest > 0 ) then
        this%Molecule%SiteTT68(i)%PXTest => this%P0Test(:, 1)
        this%Molecule%SiteTT68(i)%PYTest => this%P0Test(:, 2)
        this%Molecule%SiteTT68(i)%PZTest => this%P0Test(:, 3)
      end if

#if TRANS==1
      this%Molecule%SiteTT68(i)%Q0r => this%Q0
#endif
    end do
    do i = 1, this%Molecule%NCharge
      this%Molecule%SiteCharge(i)%NPartMax => this%NPartMax
      this%Molecule%SiteCharge(i)%NPart => this%NPart
      this%Molecule%SiteCharge(i)%NTest => this%NTest
      this%Molecule%SiteCharge(i)%NPart0 => this%NPart0
      this%Molecule%SiteCharge(i)%NPart1 => this%NPart1
      this%Molecule%SiteCharge(i)%NPart2 => this%NPart2

      call Allocate( this%Molecule%SiteCharge(i) )
      this%Molecule%SiteCharge(i)%PX => this%P0(:, 1)
      this%Molecule%SiteCharge(i)%PY => this%P0(:, 2)
      this%Molecule%SiteCharge(i)%PZ => this%P0(:, 3)

      if( ntest > 0 ) then
        this%Molecule%SiteCharge(i)%PXTest => this%P0Test(:, 1)
        this%Molecule%SiteCharge(i)%PYTest => this%P0Test(:, 2)
        this%Molecule%SiteCharge(i)%PZTest => this%P0Test(:, 3)
      end if

#if TRANS==1
      this%Molecule%SiteCharge(i)%Q0r => this%Q0
#endif
    end do

    do i = 1, this%Molecule%NDipole
      this%Molecule%SiteDipole(i)%NPartMax => this%NPartMax
      this%Molecule%SiteDipole(i)%NPart => this%NPart
      this%Molecule%SiteDipole(i)%NTest => this%NTest
      this%Molecule%SiteDipole(i)%NPart0 => this%NPart0
      this%Molecule%SiteDipole(i)%NPart1 => this%NPart1
      this%Molecule%SiteDipole(i)%NPart2 => this%NPart2

      call Allocate( this%Molecule%SiteDipole(i) )
      this%Molecule%SiteDipole(i)%PX => this%P0(:, 1)
      this%Molecule%SiteDipole(i)%PY => this%P0(:, 2)
      this%Molecule%SiteDipole(i)%PZ => this%P0(:, 3)

      if( ntest > 0 ) then
        this%Molecule%SiteDipole(i)%PXTest => this%P0Test(:, 1)
        this%Molecule%SiteDipole(i)%PYTest => this%P0Test(:, 2)
        this%Molecule%SiteDipole(i)%PZTest => this%P0Test(:, 3)
      end if

#if TRANS==1
      this%Molecule%SiteDipole(i)%Q0r => this%Q0
#endif
    end do

    do i = 1, this%Molecule%NQuadrupole
      this%Molecule%SiteQuadrupole(i)%NPartMax => this%NPartMax
      this%Molecule%SiteQuadrupole(i)%NPart => this%NPart
      this%Molecule%SiteQuadrupole(i)%NTest => this%NTest
      this%Molecule%SiteQuadrupole(i)%NPart0 => this%NPart0
      this%Molecule%SiteQuadrupole(i)%NPart1 => this%NPart1
      this%Molecule%SiteQuadrupole(i)%NPart2 => this%NPart2

      call Allocate( this%Molecule%SiteQuadrupole(i) )
      this%Molecule%SiteQuadrupole(i)%PX => this%P0(:, 1)
      this%Molecule%SiteQuadrupole(i)%PY => this%P0(:, 2)
      this%Molecule%SiteQuadrupole(i)%PZ => this%P0(:, 3)

      if( ntest > 0 ) then
        this%Molecule%SiteQuadrupole(i)%PXTest => this%P0Test(:, 1)
        this%Molecule%SiteQuadrupole(i)%PYTest => this%P0Test(:, 2)
        this%Molecule%SiteQuadrupole(i)%PZTest => this%P0Test(:, 3)
      end if

#if TRANS==1
      this%Molecule%SiteQuadrupole(i)%Q0r => this%Q0
#endif
    end do

    ! Fluctuating particle states
    if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
      nf = this%NFluctMax
      allocate( this%NState( 0: nf ), STAT = stat )
      call AllocationError( stat, 'fluctuating particle states', nf + 1 )
      allocate( this%NStateWF( 0: nf ), STAT = stat )
      call AllocationError( stat, 'fluctuating particle states', nf + 1 )

!DEBUG
      allocate( this%NFluctUpAttempts( nf ), STAT = stat )
      call AllocationError( stat, 'fluctuating particle states', nf )
      allocate( this%NFluctUpSuccesses( nf ), STAT = stat )
      call AllocationError( stat, 'fluctuating particle states', nf )
      allocate( this%NFluctDownAttempts( nf ), STAT = stat )
      call AllocationError( stat, 'fluctuating particle states', nf )
      allocate( this%NFluctDownSuccesses( nf ), STAT = stat )
      call AllocationError( stat, 'fluctuating particle states', nf )
!DEBUG
    end if

    if( this%ChemPotMethod .eq. ChemPotMethodThermoInt ) then
      allocate( this%BinsVisit( 0: this%NBins-1 ), STAT = stat )
      call AllocationError( stat, 'Number of Bins', this%NBins )
      allocate( this%BinsEn( 0: this%NBins-1 ), STAT = stat )
      call AllocationError( stat, 'En', this%NBins )
      allocate( this%BinsdEndLa( 0: this%NBins-1 ), STAT = stat )
      call AllocationError( stat, 'dEndLa', this%NBins )
      allocate( this%BinsIntdEndLa( 0: this%NBins-1 ), STAT = stat )
      call AllocationError( stat, 'IntdEndLa', this%NBins )
      allocate( this%BinsdEndLaV( 0: this%NBins-1 ), STAT = stat )
      call AllocationError( stat, 'dEndLaV', this%NBins )
      allocate( this%BinsdEndLaH( 0: this%NBins-1 ), STAT = stat )
      call AllocationError( stat, 'dEndLaH', this%NBins )
      allocate( this%BinsIntVW( 0: this%NBins-1 ), STAT = stat )
      call AllocationError( stat, 'IntVW', this%NBins )
      allocate( this%BinsIntHW( 0: this%NBins-1 ), STAT = stat )
      call AllocationError( stat, 'IntHW', this%NBins )
    end if

    ! Update log file
    write( IOBuffer, '("Memory for ", A, " allocated successfully")' ) trim( this%PotModFileName )
    call LogWrite

  end subroutine TComponent_Allocate



!==============================================================!
!  Subroutine TComponent_Deallocate                            !
!==============================================================!

  subroutine TComponent_Deallocate( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Declare local variables
    integer :: i

    ! Centers of mass positions and their derivatives
    if( associated( this%P0 ) ) then
      deallocate( this%P0 )
    end if
    if( associated( this%P0Save ) ) then
      deallocate( this%P0Save )
    end if
    if( associated( this%P0old ) ) then
      deallocate( this%P0old )
    end if
    if( associated( this%P1 ) ) then
      deallocate( this%P1 )
    end if
    if( associated( this%P2 ) ) then
      deallocate( this%P2 )
    end if
    if( associated( this%P3 ) ) then
      deallocate( this%P3 )
    end if
    if( associated( this%P4 ) ) then
      deallocate( this%P4 )
    end if
    if( associated( this%P5 ) ) then
      deallocate( this%P5 )
    end if

    ! Displacement
    if( associated( this%Disp ) ) then
      deallocate( this%Disp )
    end if
    ! Alpha2
    if( associated( this%ri0_x ) ) then
      deallocate( this%ri0_x )
    end if
    if( associated( this%ri0_y ) ) then
      deallocate( this%ri0_y )
    end if
    if( associated( this%ri0_z ) ) then
      deallocate( this%ri0_z )
    end if

    ! Total forces
    if( associated( this%F ) ) then
      deallocate( this%F )
    end if

    ! Quaternion parameters and their derivatives
    if( associated( this%Q0 ) ) then
      deallocate( this%Q0 )
    end if
    if( associated( this%Q0Save ) ) then
      deallocate( this%Q0Save )
    end if
    if( associated( this%Q0tmp ) ) then
      deallocate( this%Q0tmp )
    end if
    if( associated( this%Q1 ) ) then
      deallocate( this%Q1 )
    end if
    if( associated( this%Q2 ) ) then
      deallocate( this%Q2 )
    end if
    if( associated( this%Q3 ) ) then
      deallocate( this%Q3 )
    end if
    if( associated( this%Q4 ) ) then
      deallocate( this%Q4 )
    end if

    ! Angular velocities and their derivatives
    if( associated( this%W0 ) ) then
      deallocate( this%W0 )
    end if
    if( associated( this%W1 ) ) then
      deallocate( this%W1 )
    end if
    if( associated( this%W2 ) ) then
      deallocate( this%W2 )
    end if
    if( associated( this%W3 ) ) then
      deallocate( this%W3 )
    end if
    if( associated( this%W4 ) ) then
      deallocate( this%W4 )
    end if

    ! Total torques
    if( associated( this%T ) ) then
      deallocate( this%T )
    end if

    ! Total dipole moment of molecules for reaction field
    if( associated( this%MueX ) ) then
      deallocate( this%MueX )
    end if
    if( associated( this%MueY ) ) then
      deallocate( this%MueY )
    end if
    if( associated( this%MueZ ) ) then
      deallocate( this%MueZ )
    end if

    ! Torques from reaction field
    if( associated( this%tRFX ) ) then
      deallocate( this%tRFX )
    end if
    if( associated( this%tRFY ) ) then
      deallocate( this%tRFY )
    end if
    if( associated( this%tRFZ ) ) then
      deallocate( this%tRFZ )
    end if

    ! Total dipole moment of test particles for reaction field
    if( associated( this%MueXTest ) ) then
      deallocate( this%MueXTest )
    end if
    if( associated( this%MueYTest ) ) then
      deallocate( this%MueYTest )
    end if
    if( associated( this%MueZTest ) ) then
      deallocate( this%MueZTest )
    end if

    ! Gear corrector local arrays
    if( associated( this%Corr0 ) ) then
      deallocate( this%Corr0 )
    end if
    if( associated( this%Corr1 ) ) then
      deallocate( this%Corr1 )
    end if

#if OSMOP > 0
    ! Total forces for osmotic pressure
    if( associated( this%FOsmoticPressure ) ) then
      deallocate( this%FOsmoticPressure )
    end if
    ! Density Profile
    if( associated( this%DensityProfileN ) ) then
      deallocate( this%DensityProfileN )
    end if
#if OSMOP == 2
    if( associated( this%ChemPotProfile ) ) then
      deallocate( this%ChemPotProfile )
    end if
#endif
#endif

#if  TRANS == 1
! Transport !TRANSPORT_start
    !EinsteinCoef ri0_E deallocate
    if( associated( this%ri0_E_x ) ) then
      deallocate( this%ri0_E_x )
    end if
    if( associated( this%ri0_E_y ) ) then
      deallocate( this%ri0_E_y )
    end if
    if( associated( this%ri0_E_z ) ) then
      deallocate( this%ri0_E_z )
    end if


    if( associated( this%KinETran) ) then
      deallocate( this%KinETran )
    end if
    if( associated( this%KinEPart) ) then
      deallocate( this%KinEPart )
    end if
    if( associated( this%FS ) ) then
      deallocate( this%FS )
    end if
    if( associated( this%FB ) ) then
      deallocate( this%FB )
    end if
    if( associated( this%FTC ) ) then
      deallocate( this%FTC )
    end if
    if( associated( this%FRC ) ) then
      deallocate( this%FRC )
    end if
    if( associated( this%FTC1 ) ) then
      deallocate( this%FTC1 )
    end if
    if( associated( this%FTC2 ) ) then
      deallocate( this%FTC2 )
    end if
    if( associated( this%FTC3 ) ) then
      deallocate( this%FTC3 )
    end if
    if( associated( this%FRC1 ) ) then
      deallocate( this%FRC1 )
    end if
    if( associated( this%FRC2 ) ) then
      deallocate( this%FRC2 )
    end if
    if( associated( this%FRC3 ) ) then
      deallocate( this%FRC3 )
    end if
#if MPI_VER > 0
    if( associated( this%FBAll ) ) then
      deallocate( this%FBAll )
    end if
    if( associated( this%FSAll ) ) then
      deallocate( this%FSAll )
    end if
    if( associated( this%FRCAll ) ) then
      deallocate( this%FRCAll )
    end if
    if( associated( this%FTC1All ) ) then
      deallocate( this%FTC1All )
    end if
    if( associated( this%FTC2All ) ) then
      deallocate( this%FTC2All )
    end if
    if( associated( this%FTC3All ) ) then
      deallocate( this%FTC3All )
    end if
    if( associated( this%FRC1All ) ) then
      deallocate( this%FRC1All )
    end if
    if( associated( this%FRC2All ) ) then
      deallocate( this%FRC2All )
    end if
    if( associated( this%FRC3All ) ) then
      deallocate( this%FRC3All )
    end if
#endif
!TRANSPORT_END
#endif

    ! Site positions, orientations, forces and torques
    do i = 1, this%Molecule%NMIEnm
      call Deallocate( this%Molecule%SiteMIEnm(i) )
    end do
    do i = 1, this%Molecule%NTT68
      call Deallocate( this%Molecule%SiteTT68(i) )
    end do
    do i = 1, this%Molecule%NCharge
      call Deallocate( this%Molecule%SiteCharge(i) )
    end do
    do i = 1, this%Molecule%NDipole
      call Deallocate( this%Molecule%SiteDipole(i) )
    end do
    do i = 1, this%Molecule%NQuadrupole
      call Deallocate( this%Molecule%SiteQuadrupole(i) )
    end do

    if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
     ! Fluctuating particle states
      if( associated( this%NState ) ) then
        deallocate( this%NState )
      end if
      if( associated( this%NStateWF ) ) then
        deallocate( this%NStateWF )
      end if
      if( associated( this%NFluctComp ) ) then
        deallocate( this%NFluctComp )
      end if
      if( associated( this%WF ) ) then
        deallocate( this%WF )
      end if
!DEBUG
      if( associated( this%NFluctUpAttempts ) ) then
        deallocate( this%NFluctUpAttempts )
      end if
      if( associated( this%NFluctUpSuccesses ) ) then
        deallocate( this%NFluctUpSuccesses )
      end if
      if( associated( this%NFluctDownAttempts ) ) then
        deallocate( this%NFluctDownAttempts )
      end if
      if( associated( this%NFluctDownSuccesses ) ) then
        deallocate( this%NFluctDownSuccesses )
      end if
!DEBUG
    end if

    if( this%ChemPotMethod .eq. ChemPotMethodThermoInt ) then
      if( associated( this%BinsVisit ) ) then
        deallocate( this%BinsVisit )
      end if
      if( associated( this%BinsEn ) ) then
        deallocate( this%BinsEn )
      end if
      if( associated( this%BinsdEndLa ) ) then
        deallocate( this%BinsdEndLa )
      end if
      if( associated( this%BinsIntdEndLa ) ) then
        deallocate( this%BinsIntdEndLa )
      end if
      if( associated( this%BinsdEndLaV ) ) then
        deallocate( this%BinsdEndLaV )
      end if
      if( associated( this%BinsdEndLaH ) ) then
        deallocate( this%BinsdEndLaH )
      end if
      if( associated( this%BinsIntVW ) ) then
        deallocate( this%BinsIntVW )
      end if
      if( associated( this%BinsIntHW ) ) then
        deallocate( this%BinsIntHW )
      end if
    end if

#if MPI_VER > 0
    if( associated( this%FAll ) ) then
      deallocate( this%FAll )
    end if
    if( associated( this%TAll ) ) then
      deallocate( this%TAll )
    end if
#endif

  end subroutine TComponent_Deallocate


!==============================================================!
!  Subroutine TComponent_LongRangeCheck                        !
!==============================================================!

   subroutine TComponent_LongRangeCheck ( this, q )

   implicit none
   type(TComponent)           :: this
   real(RK),intent(in out)    :: q
   integer                    :: i

   do i=1,this%Molecule%NCharge
     q = q + this%NPart * this%Molecule%SiteCharge(i)%e
   end do

   if (abs(q) .ge. 1e-1) this%charged = .true.

! Reaction Field Check
     if ((LongRange .eq. RField) .and. (abs(this%Molecule%Charge) .ge. 1e-1)) then
       write (ErrorBuffer,'("You have a non-neutral component.\n NetCharge norm&
&      red = ", F15.10, "\n Conflicts with ReactionField")') this%Molecule%Charge
       call Error
     end if

     if ( ((EnsembleType .eq. EnsembleTypeGE) .or. (EnsembleType .eq. EnsembleTypeHA)) .and. (abs(q) .ge. 1e-1) ) then
       write (ErrorBuffer,'("GrandEquilibrium not possible in a charged system")') q
       call Error
     end if

   end subroutine TComponent_LongRangeCheck


!==============================================================!
!  Subroutine TComponent_InitVelocities                        !
!==============================================================!

  subroutine TComponent_InitVelocities( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Declare local variables
    integer :: i, j

    ! Set random linear velocities
    do i = 1, 3
      do j = 1, this%NPart
        this%P1(j, i) = rnd( -1._RK, 1._RK )
      end do
    end do

    ! Normalize translational velocity vectors (only done once - needs not to be efficient)
    do j = 1, this%NPart
      this%P1(j, :) = this%P1(j, :) / sqrt( dot_product( this%P1(j, :), this%P1(j, :) ))
    end do

    ! Nullify angular velocities
    if( this%Molecule%isElongated ) this%W0(:, :) = 0._RK

  end subroutine TComponent_InitVelocities


!==============================================================!
!  Subroutine TComponent_InitIntegratorGear                    !
!==============================================================!

  subroutine TComponent_InitIntegratorGear( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Zero accelerations
    this%P2(:, :) = 0._RK
    this%P3(:, :) = 0._RK
    this%P4(:, :) = 0._RK
    this%P5(:, :) = 0._RK
    if( this%Molecule%isElongated ) then
      this%Q1(:, :) = 0._RK
      this%Q2(:, :) = 0._RK
      this%Q3(:, :) = 0._RK
      this%Q4(:, :) = 0._RK
      this%W1(:, :) = 0._RK
      this%W2(:, :) = 0._RK
      this%W3(:, :) = 0._RK
      this%W4(:, :) = 0._RK
    end if

  end subroutine TComponent_InitIntegratorGear


!==============================================================!
!  Subroutine TComponent_InitIntegratorLeap                    !
!==============================================================!

  subroutine TComponent_InitIntegratorLeap( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Zero accelerations
    this%P2(:, :) = 0._RK

    if( this%Molecule%isElongated ) then
      this%Q1(:, :) = 0._RK
      this%W1(:, :) = 0._RK
    end if

  end subroutine TComponent_InitIntegratorLeap


!==============================================================!
!  Subroutine TComponent_InitIntegratorVerlet                  !
!==============================================================!

  subroutine TComponent_InitIntegratorVerlet( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Issue error
    call Error( 'Subroutine TComponent_InitIntegratorVerlet is not implemented' )

  end subroutine TComponent_InitIntegratorVerlet


!==============================================================!
!  Subroutine TComponent_InitIntegratorVV                      !
!==============================================================!

  subroutine TComponent_InitIntegratorVV( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Issue error
    call Error( 'Subroutine TComponent_InitIntegratorVV is not implemented' )

  end subroutine TComponent_InitIntegratorVV


!==============================================================!
!  Subroutine TComponent_RemoveNetMomentum                     !
!==============================================================!

  subroutine TComponent_RemoveNetMomentum( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Declare local variables
    real(RK) :: P(3), L(3)
    integer :: i, j
    real(RK) :: Pim

    ! Return if zero particles in component
    if( this%NPart == 0 ) return

    ! Calculate net momentum
    P(:) = 0._RK
    L(:) = 0._RK
    do i = 1, 3
      P(i) = P(i) + this%Molecule%Mass * sum( this%P1(1:this%NPart, i) )

      if( i <= this%Molecule%NDFRot ) L(i) = L(i) + this%Molecule%MOI(i) * sum( this%W0(1:this%NPart, i) )

    end do
    P(:) = P(:) / this%NPart
    L(:) = L(:) / this%NPart

    ! Remove net momentum
    do i = 1, 3
      Pim = P(i) / this%Molecule%Mass
      do j = 1, this%NPart
        this%P1(j, i) = this%P1(j, i) - Pim
      end do

      if( i <= this%Molecule%NDFRot ) then
        Pim = L(i) / this%Molecule%MOI(i)
        do j = 1, this%NPart
          this%W0(j, i) = this%W0(j, i) - Pim
        end do
      end if
    end do

  end subroutine TComponent_RemoveNetMomentum


!==============================================================!
!  Subroutine TComponent_CalculateEKin                         !
!==============================================================!

  subroutine TComponent_CalculateEKin( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Declare local variables
    integer :: i

    ! Calculate translational kinetic energy
    this%EKinTran = this%Molecule%Mass * TimeStepSquaredInv2 &
&     * sum( this%P1(1:this%NPart, :)**2 ) * this%BoxLength**2

    ! Calculate rotational kinetic energy
    this%EKinRot = 0._RK
    do i = 1, this%Molecule%NDFRot
      this%EKinRot = this%EKinRot + this%Molecule%MOI(i) * .5_RK &
&       * sum( this%W0(1:this%NPart, i)**2 )
    end do

  end subroutine TComponent_CalculateEKin


!==============================================================!
!  Subroutine TComponent_Mol2Atom                              !
!==============================================================!

!   subroutine TComponent_Mol2Atom( this, i0, l )
    subroutine TComponent_Mol2Atom( this, l )
    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent)    :: this
!     integer, intent(in) :: i0
    integer, intent(in) :: l

    ! Declare local variables
    real(RK)                       :: BoxLengthInv
    real(RK)                       :: PX(l), PY(l), PZ(l)
    real(RK)                       :: q1, q2, q3, q4, qinv
    real(RK)                       :: A11(l), A12(l), A13(l)
    real(RK)                       :: A21(l), A22(l), A23(l)
    real(RK)                       :: A31(l), A32(l), A33(l)
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    real(RK)                       :: mue1, mue2, mue3
    type(TSiteMIEnm), pointer      :: pMIEnm
    type(TSiteTT68), pointer       :: pTT68
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i0, i1, i, j

!     ! i0 startpos of proc, i1 endpos of proc; l = difference
!     i1 = l + i0 - 1
    i0 = 1
    i1 = l

    ! Broadcast positions and orientations to all processes
#if MPI_VER > 0
    ! in MC simulations, we only communicate during common equilibration
    if ( SimulationType .ne. MonteCarlo .or. ((Equilibration .and. CommonEqui) )) then
      call MPI_Bcast( this%P0(:, :), size( this%P0 ), MPI_RK, NRootProc, Communicator, ierror )
      if( this%Molecule%isElongated ) call MPI_Bcast( this%Q0(:, :), size( this%Q0 ),&
      & MPI_RK, NRootProc, Communicator, ierror )
    endif
#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    ! Check number of rotation axes
    if( this%Molecule%isElongated ) then

      ! Loop over molecules
      do i = 1, l
        ! Positions and quaternions of particle i
        PX(i) = this%P0(i-1+i0, 1)
        PY(i) = this%P0(i-1+i0, 2)
        PZ(i) = this%P0(i-1+i0, 3)
        q1 = this%Q0(i-1+i0, 1)
        q2 = this%Q0(i-1+i0, 2)
        q3 = this%Q0(i-1+i0, 3)
        q4 = this%Q0(i-1+i0, 4)

        ! Normalise quaternions
        qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
        q1 = q1 * qinv
        q2 = q2 * qinv
        q3 = q3 * qinv
        q4 = q4 * qinv
        this%Q0(i-1+i0, 1) = q1
        this%Q0(i-1+i0, 2) = q2
        this%Q0(i-1+i0, 3) = q3
        this%Q0(i-1+i0, 4) = q4

        ! Calculate rotation matrix elements
        A11(i) = q1**2 + q2**2 - q3**2 - q4**2
        A12(i) = 2._RK * (q2 * q3 + q1 * q4)
        A13(i) = 2._RK * (q2 * q4 - q1 * q3)
        A21(i) = 2._RK * (q2 * q3 - q1 * q4)
        A22(i) = q1**2 - q2**2 + q3**2 - q4**2
        A23(i) = 2._RK * (q3 * q4 + q1 * q2)
        A31(i) = 2._RK * (q2 * q4 + q1 * q3)
        A32(i) = 2._RK * (q3 * q4 - q1 * q2)
        A33(i) = q1**2 - q2**2 - q3**2 + q4**2
      end do

      ! Loop over MIEnm sites in molecule
      do j = 1, this%Molecule%NMIEnm
        pMIEnm => this%Molecule%SiteMIEnm(j)
        r1 = pMIEnm%r(1) * BoxLengthInv
        r2 = pMIEnm%r(2) * BoxLengthInv
        r3 = pMIEnm%r(3) * BoxLengthInv
        do i = 1, l
          pMIEnm%RX(i-1+i0) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pMIEnm%RY(i-1+i0) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pMIEnm%RZ(i-1+i0) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
        end do
      end do

      ! Loop over TT68 sites in molecule
      do j = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(j)
        r1 = pTT68%r(1) * BoxLengthInv
        r2 = pTT68%r(2) * BoxLengthInv
        r3 = pTT68%r(3) * BoxLengthInv
        do i = 1, l
          pTT68%RX(i-1+i0) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pTT68%RY(i-1+i0) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pTT68%RZ(i-1+i0) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
        end do
      end do

      ! Loop over charge sites in molecule
      do j = 1, this%Molecule%NCharge
        pCharge => this%Molecule%SiteCharge(j)
        r1 = pCharge%r(1) * BoxLengthInv
        r2 = pCharge%r(2) * BoxLengthInv
        r3 = pCharge%r(3) * BoxLengthInv
        do i = 1, l
          pCharge%RX(i-1+i0) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pCharge%RY(i-1+i0) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pCharge%RZ(i-1+i0) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
        end do
      end do

      ! Loop over dipole sites in molecule
      do j = 1, this%Molecule%NDipole
        pDipole => this%Molecule%SiteDipole(j)
        r1 = pDipole%r(1) * BoxLengthInv
        r2 = pDipole%r(2) * BoxLengthInv
        r3 = pDipole%r(3) * BoxLengthInv
        or1 = pDipole%or(1)
        or2 = pDipole%or(2)
        or3 = pDipole%or(3)
        do i = 1, l
          pDipole%RX(i-1+i0) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pDipole%RY(i-1+i0) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pDipole%RZ(i-1+i0) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
          pDipole%OX(i-1+i0) = or1 * A11(i) + or2 * A21(i) + or3 * A31(i)
          pDipole%OY(i-1+i0) = or1 * A12(i) + or2 * A22(i) + or3 * A32(i)
          pDipole%OZ(i-1+i0) = or1 * A13(i) + or2 * A23(i) + or3 * A33(i)
        end do
      end do

      ! Loop over quadrupole sites in molecule
      do j = 1, this%Molecule%NQuadrupole
        pQuadrupole => this%Molecule%SiteQuadrupole(j)
        r1 = pQuadrupole%r(1) * BoxLengthInv
        r2 = pQuadrupole%r(2) * BoxLengthInv
        r3 = pQuadrupole%r(3) * BoxLengthInv
        or1 = pQuadrupole%or(1)
        or2 = pQuadrupole%or(2)
        or3 = pQuadrupole%or(3)
        do i = 1, l
          pQuadrupole%RX(i-1+i0) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pQuadrupole%RY(i-1+i0) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pQuadrupole%RZ(i-1+i0) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
          pQuadrupole%OX(i-1+i0) = or1 * A11(i) + or2 * A21(i) + or3 * A31(i)
          pQuadrupole%OY(i-1+i0) = or1 * A12(i) + or2 * A22(i) + or3 * A32(i)
          pQuadrupole%OZ(i-1+i0) = or1 * A13(i) + or2 * A23(i) + or3 * A33(i)
        end do
      end do

      ! Rotate total dipole moment
      if( CutoffMode .eq. CenterofMass ) then
        mue1 = this%Molecule%Mue(1)
        mue2 = this%Molecule%Mue(2)
        mue3 = this%Molecule%Mue(3)
        do i = 1, l
          this%MueX(i-1+i0) = mue1 * A11(i) + mue2 * A21(i) + mue3 * A31(i)
          this%MueY(i-1+i0) = mue1 * A12(i) + mue2 * A22(i) + mue3 * A32(i)
          this%MueZ(i-1+i0) = mue1 * A13(i) + mue2 * A23(i) + mue3 * A33(i)
        end do
      end if

    else

      ! Loop over MIEnm sites in molecule
      do j = 1, this%Molecule%NMIEnm
        pMIEnm => this%Molecule%SiteMIEnm(j)
        do i = 1, l
          pMIEnm%RX(i-1+i0) = this%P0(i, 1)
          pMIEnm%RY(i-1+i0) = this%P0(i, 2)
          pMIEnm%RZ(i-1+i0) = this%P0(i, 3)
        end do
      end do

      ! Loop over TT68 sites in molecule
      do j = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(j)
        do i = 1, l
          pTT68%RX(i-1+i0) = this%P0(i, 1)
          pTT68%RY(i-1+i0) = this%P0(i, 2)
          pTT68%RZ(i-1+i0) = this%P0(i, 3)
        end do
      end do

      ! Loop over charge sites in molecule
      if (LongRange .ne. RField) then
        do j = 1, this%Molecule%NCharge
          pCharge => this%Molecule%SiteCharge(j)
          do i = 1, l
            pCharge%RX(i-1+i0) = this%P0(i,1)
            pCharge%RY(i-1+i0) = this%P0(i,2)
            pCharge%RZ(i-1+i0) = this%P0(i,3)
          end do
        end do
      end if
    end if

  end subroutine TComponent_Mol2Atom


!==============================================================!
!  Subroutine TComponent_Mol2Atom1                             !
!==============================================================!

  subroutine TComponent_Mol2Atom1( this, np )

    implicit none

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: np

    ! Declare local variables
    real(RK)                       :: BoxLengthInv
    real(RK)                       :: PXi, PYi, PZi
    real(RK)                       :: q1, q2, q3, q4, qinv
    real(RK)                       :: A11, A12, A13, A21, A22, A23, A31, A32, A33
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    real(RK)                       :: mue1, mue2, mue3
    type(TSiteMIEnm), pointer      :: pMIEnm
    type(TSiteTT68), pointer       :: pTT68
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    ! Positions of particle n
    PXi = this%P0(np, 1)
    PYi = this%P0(np, 2)
    PZi = this%P0(np, 3)

    ! Check number of rotation axes
    if( this%Molecule%isElongated ) then

      ! Normalise quaternions
      q1 = this%Q0(np, 1)
      q2 = this%Q0(np, 2)
      q3 = this%Q0(np, 3)
      q4 = this%Q0(np, 4)
      qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
      q1 = q1 * qinv
      q2 = q2 * qinv
      q3 = q3 * qinv
      q4 = q4 * qinv
      this%Q0(np, 1) = q1
      this%Q0(np, 2) = q2
      this%Q0(np, 3) = q3
      this%Q0(np, 4) = q4

      ! Calculate rotation matrix elements
      A11 = q1**2 + q2**2 - q3**2 - q4**2
      A12 = 2._RK * (q2 * q3 + q1 * q4)
      A13 = 2._RK * (q2 * q4 - q1 * q3)
      A21 = 2._RK * (q2 * q3 - q1 * q4)
      A22 = q1**2 - q2**2 + q3**2 - q4**2
      A23 = 2._RK * (q3 * q4 + q1 * q2)
      A31 = 2._RK * (q2 * q4 + q1 * q3)
      A32 = 2._RK * (q3 * q4 - q1 * q2)
      A33 = q1**2 - q2**2 - q3**2 + q4**2

      ! Loop over MIEnm sites in molecule
      do i = 1, this%Molecule%NMIEnm
        pMIEnm => this%Molecule%SiteMIEnm(i)
        r1 = pMIEnm%r(1) * BoxLengthInv
        r2 = pMIEnm%r(2) * BoxLengthInv
        r3 = pMIEnm%r(3) * BoxLengthInv
        pMIEnm%RX(np) = PXi + r1 * A11 + r2 * A21 + r3 * A31
        pMIEnm%RY(np) = PYi + r1 * A12 + r2 * A22 + r3 * A32
        pMIEnm%RZ(np) = PZi + r1 * A13 + r2 * A23 + r3 * A33
      end do

      ! Loop over TT68 sites in molecule
      do i = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(i)
        r1 = pTT68%r(1) * BoxLengthInv
        r2 = pTT68%r(2) * BoxLengthInv
        r3 = pTT68%r(3) * BoxLengthInv
        pTT68%RX(np) = PXi + r1 * A11 + r2 * A21 + r3 * A31
        pTT68%RY(np) = PYi + r1 * A12 + r2 * A22 + r3 * A32
        pTT68%RZ(np) = PZi + r1 * A13 + r2 * A23 + r3 * A33
      end do

      ! Loop over charge sites in molecule
      do i = 1, this%Molecule%NCharge
        pCharge => this%Molecule%SiteCharge(i)
        r1 = pCharge%r(1) * BoxLengthInv
        r2 = pCharge%r(2) * BoxLengthInv
        r3 = pCharge%r(3) * BoxLengthInv
        pCharge%RX(np) = PXi + r1 * A11 + r2 * A21 + r3 * A31
        pCharge%RY(np) = PYi + r1 * A12 + r2 * A22 + r3 * A32
        pCharge%RZ(np) = PZi + r1 * A13 + r2 * A23 + r3 * A33
      end do

      ! Loop over dipole sites in molecule
      do i = 1, this%Molecule%NDipole
        pDipole => this%Molecule%SiteDipole(i)
        r1 = pDipole%r(1) * BoxLengthInv
        r2 = pDipole%r(2) * BoxLengthInv
        r3 = pDipole%r(3) * BoxLengthInv
        or1 = pDipole%or(1)
        or2 = pDipole%or(2)
        or3 = pDipole%or(3)
        pDipole%RX(np) = PXi + r1 * A11 + r2 * A21 + r3 * A31
        pDipole%RY(np) = PYi + r1 * A12 + r2 * A22 + r3 * A32
        pDipole%RZ(np) = PZi + r1 * A13 + r2 * A23 + r3 * A33
        pDipole%OX(np) = or1 * A11 + or2 * A21 + or3 * A31
        pDipole%OY(np) = or1 * A12 + or2 * A22 + or3 * A32
        pDipole%OZ(np) = or1 * A13 + or2 * A23 + or3 * A33
      end do

      ! Loop over quadrupole sites in molecule
      do i = 1, this%Molecule%NQuadrupole
        pQuadrupole => this%Molecule%SiteQuadrupole(i)
        r1 = pQuadrupole%r(1) * BoxLengthInv
        r2 = pQuadrupole%r(2) * BoxLengthInv
        r3 = pQuadrupole%r(3) * BoxLengthInv
        or1 = pQuadrupole%or(1)
        or2 = pQuadrupole%or(2)
        or3 = pQuadrupole%or(3)
        pQuadrupole%RX(np) = PXi + r1 * A11 + r2 * A21 + r3 * A31
        pQuadrupole%RY(np) = PYi + r1 * A12 + r2 * A22 + r3 * A32
        pQuadrupole%RZ(np) = PZi + r1 * A13 + r2 * A23 + r3 * A33
        pQuadrupole%OX(np) = or1 * A11 + or2 * A21 + or3 * A31
        pQuadrupole%OY(np) = or1 * A12 + or2 * A22 + or3 * A32
        pQuadrupole%OZ(np) = or1 * A13 + or2 * A23 + or3 * A33
      end do

      ! Rotate total dipole moment
      if( CutoffMode .eq. CenterofMass ) then
        mue1 = this%Molecule%Mue(1)
        mue2 = this%Molecule%Mue(2)
        mue3 = this%Molecule%Mue(3)
        this%MueX(np) = mue1 * A11 + mue2 * A21 + mue3 * A31
        this%MueY(np) = mue1 * A12 + mue2 * A22 + mue3 * A32
        this%MueZ(np) = mue1 * A13 + mue2 * A23 + mue3 * A33
      end if

    else

      ! Loop over MIEnm sites in molecule
      do i = 1, this%Molecule%NMIEnm
        pMIEnm => this%Molecule%SiteMIEnm(i)
        pMIEnm%RX(np) = PXi
        pMIEnm%RY(np) = PYi
        pMIEnm%RZ(np) = PZi
      end do

      ! Loop over TT68 sites in molecule
      do i = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(i)
        pTT68%RX(np) = PXi
        pTT68%RY(np) = PYi
        pTT68%RZ(np) = PZi
      end do

      ! Loop over charge sites in molecule
      if (LongRange .ne. RField) then
        do i = 1, this%Molecule%NCharge
          pCharge => this%Molecule%SiteCharge(i)
            pCharge%RX(np) = PXi
            pCharge%RY(np) = PYi
            pCharge%RZ(np) = PZi
        end do
      end if

    end if

  end subroutine TComponent_Mol2Atom1


!==============================================================!
!  Subroutine TComponent_Mol2AtomTest                          !
!==============================================================!

  subroutine TComponent_Mol2AtomTest( this, np )

    implicit none

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: np

    ! Declare local variables
    real(RK)                       :: BoxLengthInv
    real(RK)                       :: PX(np), PY(np), PZ(np)
    real(RK)                       :: q1, q2, q3, q4
    real(RK)                       :: A11(np), A12(np), A13(np)
    real(RK)                       :: A21(np), A22(np), A23(np)
    real(RK)                       :: A31(np), A32(np), A33(np)
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    real(RK)                       :: mue1, mue2, mue3
    type(TSiteMIEnm), pointer      :: pMIEnm
    type(TSiteTT68), pointer       :: pTT68
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    ! Check number of rotation axes
    if( this%Molecule%isElongated ) then

      ! Loop over molecules
      do i = 1, np
        ! Positions and quaternions of test particle i
        PX(i) = this%P0Test(i, 1)
        PY(i) = this%P0Test(i, 2)
        PZ(i) = this%P0Test(i, 3)
        q1 = this%Q0Test(i, 1)
        q2 = this%Q0Test(i, 2)
        q3 = this%Q0Test(i, 3)
        q4 = this%Q0Test(i, 4)

        ! Calculate rotation matrix elements
        A11(i) = q1**2 + q2**2 - q3**2 - q4**2
        A12(i) = 2._RK * (q2 * q3 + q1 * q4)
        A13(i) = 2._RK * (q2 * q4 - q1 * q3)
        A21(i) = 2._RK * (q2 * q3 - q1 * q4)
        A22(i) = q1**2 - q2**2 + q3**2 - q4**2
        A23(i) = 2._RK * (q3 * q4 + q1 * q2)
        A31(i) = 2._RK * (q2 * q4 + q1 * q3)
        A32(i) = 2._RK * (q3 * q4 - q1 * q2)
        A33(i) = q1**2 - q2**2 - q3**2 + q4**2
      end do

      ! Loop over MIEnm sites in molecule
      do j = 1, this%Molecule%NMIEnm
        pMIEnm => this%Molecule%SiteMIEnm(j)
        r1 = pMIEnm%r(1) * BoxLengthInv
        r2 = pMIEnm%r(2) * BoxLengthInv
        r3 = pMIEnm%r(3) * BoxLengthInv
        do i = 1, np
          pMIEnm%RXTest(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pMIEnm%RYTest(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pMIEnm%RZTest(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
        end do
      end do

      ! Loop over TT68 sites in molecule
      do j = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(j)
        r1 = pTT68%r(1) * BoxLengthInv
        r2 = pTT68%r(2) * BoxLengthInv
        r3 = pTT68%r(3) * BoxLengthInv
        do i = 1, np
          pTT68%RXTest(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pTT68%RYTest(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pTT68%RZTest(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
        end do
      end do

      ! Loop over charge sites in molecule
      do j = 1, this%Molecule%NCharge
        pCharge => this%Molecule%SiteCharge(j)
        r1 = pCharge%r(1) * BoxLengthInv
        r2 = pCharge%r(2) * BoxLengthInv
        r3 = pCharge%r(3) * BoxLengthInv
        do i = 1, np
          pCharge%RXTest(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pCharge%RYTest(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pCharge%RZTest(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
        end do
      end do

      ! Loop over dipole sites in molecule
      do j = 1, this%Molecule%NDipole
        pDipole => this%Molecule%SiteDipole(j)
        r1 = pDipole%r(1) * BoxLengthInv
        r2 = pDipole%r(2) * BoxLengthInv
        r3 = pDipole%r(3) * BoxLengthInv
        or1 = pDipole%or(1)
        or2 = pDipole%or(2)
        or3 = pDipole%or(3)
        do i = 1, np
          pDipole%RXTest(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pDipole%RYTest(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pDipole%RZTest(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
          pDipole%OXTest(i) = or1 * A11(i) + or2 * A21(i) + or3 * A31(i)
          pDipole%OYTest(i) = or1 * A12(i) + or2 * A22(i) + or3 * A32(i)
          pDipole%OZTest(i) = or1 * A13(i) + or2 * A23(i) + or3 * A33(i)
        end do
      end do

      ! Loop over quadrupole sites in molecule
      do j = 1, this%Molecule%NQuadrupole
        pQuadrupole => this%Molecule%SiteQuadrupole(j)
        r1 = pQuadrupole%r(1) * BoxLengthInv
        r2 = pQuadrupole%r(2) * BoxLengthInv
        r3 = pQuadrupole%r(3) * BoxLengthInv
        or1 = pQuadrupole%or(1)
        or2 = pQuadrupole%or(2)
        or3 = pQuadrupole%or(3)
        do i = 1, np
          pQuadrupole%RXTest(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pQuadrupole%RYTest(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pQuadrupole%RZTest(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)

          pQuadrupole%OXTest(i) = or1 * A11(i) + or2 * A21(i) + or3 * A31(i)
          pQuadrupole%OYTest(i) = or1 * A12(i) + or2 * A22(i) + or3 * A32(i)
          pQuadrupole%OZTest(i) = or1 * A13(i) + or2 * A23(i) + or3 * A33(i)
        end do
      end do

      if( CutoffMode .eq. CenterofMass ) then
        mue1 = this%Molecule%Mue(1)
        mue2 = this%Molecule%Mue(2)
        mue3 = this%Molecule%Mue(3)
        do i = 1, np
          this%MueXTest(i) = mue1 * A11(i) + mue2 * A21(i) + mue3 * A31(i)
          this%MueYTest(i) = mue1 * A12(i) + mue2 * A22(i) + mue3 * A32(i)
          this%MueZTest(i) = mue1 * A13(i) + mue2 * A23(i) + mue3 * A33(i)
        end do
      end if

    else

      ! Loop over MIEnm sites in molecule
      do i = 1, this%Molecule%NMIEnm
        pMIEnm => this%Molecule%SiteMIEnm(i)
        do j = 1, np
          pMIEnm%RXTest(j) = this%P0Test(j, 1)
          pMIEnm%RYTest(j) = this%P0Test(j, 2)
          pMIEnm%RZTest(j) = this%P0Test(j, 3)
        end do
      end do

      ! Loop over TT68 sites in molecule
      do i = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(i)
        do j = 1, np
          pTT68%RXTest(j) = this%P0Test(j, 1)
          pTT68%RYTest(j) = this%P0Test(j, 2)
          pTT68%RZTest(j) = this%P0Test(j, 3)
        end do
      end do

      ! Loop over charge sites in molecule
      if (LongRange .ne. RField) then
        do i = 1, this%Molecule%NCharge
          pCharge => this%Molecule%SiteCharge(i)
          do j = 1, np
            pCharge%RX(j) = this%P0Test(j,1)
            pCharge%RY(j) = this%P0Test(j,2)
            pCharge%RZ(j) = this%P0Test(j,3)
          end do
        end do
      end if

    end if

  end subroutine TComponent_Mol2AtomTest


!==============================================================!
!  Subroutine TComponent_Atom2Mol                              !
!==============================================================!

!   subroutine TComponent_Atom2Mol( this, i0, l )
    subroutine TComponent_Atom2Mol( this, l )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent)    :: this
!     integer, intent(in) :: i0
    integer, intent(in) :: l

    ! Declare local variables
    real(RK)                       :: BoxLength
    real(RK)                       :: rx(l), ry(l), rz(l), r1x, r1y, r1z
    real(RK)                       :: q1(l), q2(l), q3(l), q4(l)
    real(RK)                       :: fx, fy, fz, tx, ty, tz
    real(RK)                       :: A11, A12, A13, A21, A22, A23, A31, A32, A33
    type(TSiteMIEnm), pointer      :: pMIEnm
    type(TSiteTT68), pointer       :: pTT68
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i0, i1, i, j
#if MPI_VER > 0
    real(RK),allocatable           :: OsmoPAll(:)

    allocate( OsmoPAll(this%NPart) )
    OsmoPAll(:) = 0._RK
#endif

!     ! i0 startpos of proc, i1 endpos of proc; i = difference
!     i1 = l + i0 - 1
    i0 = 1
    i1 = l

    ! Assign local variables
    BoxLength = this%BoxLength

    ! Initialize forces
    this%F(:, :) = 0._RK
#if OSMOP > 0
    !if (RootProc) then
      this%FOsmoticPressure(:) = 0._RK
      if ( .not. this%permeable ) then
        do i = i0, i1 ! 1, this%NPart
          if( this%P0(i,1) .ge. 0.25_RK ) then
            this%FOsmoticPressure(i) = kForceOsmoticPressure*( this%P0(i,1)-.25_RK)*BoxLength
            this%F(i,1) = this%F(i,1) - this%FOsmoticPressure(i)
          elseif( this%P0(i,1) .le. -0.25_RK ) then
            this%FOsmoticPressure(i) = kForceOsmoticPressure*(-this%P0(i,1)-.25_RK)*BoxLength
            this%F(i,1) = this%F(i,1) + this%FOsmoticPressure(i)
          end if
        end do
      end if
    !end if
#if MPI_VER > 0
    call MPI_Reduce( this%FOsmoticPressure(:), OsmoPAll(:), size( this%FOsmoticPressure), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    if (RootProc) this%FOsmoticPressure(:) = OsmoPAll(:)
#endif
#endif

    ! Check number of rotation axes
    if( this%Molecule%isElongated ) then

      ! Initialize torques
      this%T(:, :) = 0._RK

      ! Initialize local arrays
      rx(:) = this%P0(i0:i1, 1)
      ry(:) = this%P0(i0:i1, 2)
      rz(:) = this%P0(i0:i1, 3)
      q1(:) = this%Q0(i0:i1, 1)
      q2(:) = this%Q0(i0:i1, 2)
      q3(:) = this%Q0(i0:i1, 3)
      q4(:) = this%Q0(i0:i1, 4)

      ! Loop over MIEnm sites in molecule
      do j = 1, this%Molecule%NMIEnm
        pMIEnm => this%Molecule%SiteMIEnm(j)
        do i = 1, l
          fx = pMIEnm%FX(i-1+i0)
          fy = pMIEnm%FY(i-1+i0)
          fz = pMIEnm%FZ(i-1+i0)
          r1x = ( pMIEnm%RX(i-1+i0) - rx(i) ) * BoxLength
          r1y = ( pMIEnm%RY(i-1+i0) - ry(i) ) * BoxLength
          r1z = ( pMIEnm%RZ(i-1+i0) - rz(i) ) * BoxLength
          this%F(i-1+i0, 1) = this%F(i-1+i0, 1) + fx
          this%F(i-1+i0, 2) = this%F(i-1+i0, 2) + fy
          this%F(i-1+i0, 3) = this%F(i-1+i0, 3) + fz
          this%T(i-1+i0, 1) = this%T(i-1+i0, 1) + r1y * fz - r1z * fy
          this%T(i-1+i0, 2) = this%T(i-1+i0, 2) + r1z * fx - r1x * fz
          this%T(i-1+i0, 3) = this%T(i-1+i0, 3) + r1x * fy - r1y * fx
        end do
      end do

      ! Loop over TT68 sites in molecule
      do j = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(j)
        do i = 1, l
          fx = pTT68%FX(i-1+i0)
          fy = pTT68%FY(i-1+i0)
          fz = pTT68%FZ(i-1+i0)
          r1x = ( pTT68%RX(i-1+i0) - rx(i) ) * BoxLength
          r1y = ( pTT68%RY(i-1+i0) - ry(i) ) * BoxLength
          r1z = ( pTT68%RZ(i-1+i0) - rz(i) ) * BoxLength
          this%F(i-1+i0, 1) = this%F(i-1+i0, 1) + fx
          this%F(i-1+i0, 2) = this%F(i-1+i0, 2) + fy
          this%F(i-1+i0, 3) = this%F(i-1+i0, 3) + fz
          this%T(i-1+i0, 1) = this%T(i-1+i0, 1) + r1y * fz - r1z * fy
          this%T(i-1+i0, 2) = this%T(i-1+i0, 2) + r1z * fx - r1x * fz
          this%T(i-1+i0, 3) = this%T(i-1+i0, 3) + r1x * fy - r1y * fx
        end do
      end do

      ! Loop over charge sites in molecule
      do j = 1, this%Molecule%NCharge
        pCharge => this%Molecule%SiteCharge(j)
        do i = 1, l
          fx = pCharge%FX(i-1+i0)
          fy = pCharge%FY(i-1+i0)
          fz = pCharge%FZ(i-1+i0)
          r1x = ( pCharge%RX(i-1+i0) - rx(i) ) * BoxLength
          r1y = ( pCharge%RY(i-1+i0) - ry(i) ) * BoxLength
          r1z = ( pCharge%RZ(i-1+i0) - rz(i) ) * BoxLength
          this%F(i-1+i0, 1) = this%F(i-1+i0, 1) + fx
          this%F(i-1+i0, 2) = this%F(i-1+i0, 2) + fy
          this%F(i-1+i0, 3) = this%F(i-1+i0, 3) + fz
          this%T(i-1+i0, 1) = this%T(i-1+i0, 1) + r1y * fz - r1z * fy
          this%T(i-1+i0, 2) = this%T(i-1+i0, 2) + r1z * fx - r1x * fz
          this%T(i-1+i0, 3) = this%T(i-1+i0, 3) + r1x * fy - r1y * fx
        end do
      end do

      ! Loop over dipole sites in molecule
      do j = 1, this%Molecule%NDipole
        pDipole => this%Molecule%SiteDipole(j)
        do i = 1, l
          fx = pDipole%FX(i-1+i0)
          fy = pDipole%FY(i-1+i0)
          fz = pDipole%FZ(i-1+i0)
          r1x = ( pDipole%RX(i-1+i0) - rx(i) ) * BoxLength
          r1y = ( pDipole%RY(i-1+i0) - ry(i) ) * BoxLength
          r1z = ( pDipole%RZ(i-1+i0) - rz(i) ) * BoxLength
          this%F(i-1+i0, 1) = this%F(i-1+i0, 1) + fx
          this%F(i-1+i0, 2) = this%F(i-1+i0, 2) + fy
          this%F(i-1+i0, 3) = this%F(i-1+i0, 3) + fz
          this%T(i-1+i0, 1) = this%T(i-1+i0, 1) + pDipole%OY(i) * pDipole%TZ(i) &
&                                    - pDipole%OZ(i) * pDipole%TY(i) + r1y * fz - r1z * fy
          this%T(i-1+i0, 2) = this%T(i-1+i0, 2) + pDipole%OZ(i) * pDipole%TX(i) &
&                                    - pDipole%OX(i) * pDipole%TZ(i) + r1z * fx - r1x * fz
          this%T(i-1+i0, 3) = this%T(i-1+i0, 3) + pDipole%OX(i) * pDipole%TY(i) &
&                                    - pDipole%OY(i) * pDipole%TX(i) + r1x * fy - r1y * fx
        end do
      end do

      ! Loop over quadrupole sites in molecule
      do j = 1, this%Molecule%NQuadrupole
        pQuadrupole => this%Molecule%SiteQuadrupole(j)
        do i = 1, l
          fx = pQuadrupole%FX(i-1+i0)
          fy = pQuadrupole%FY(i-1+i0)
          fz = pQuadrupole%FZ(i-1+i0)
          r1x = ( pQuadrupole%RX(i-1+i0) - rx(i) ) * BoxLength
          r1y = ( pQuadrupole%RY(i-1+i0) - ry(i) ) * BoxLength
          r1z = ( pQuadrupole%RZ(i-1+i0) - rz(i) ) * BoxLength
          this%F(i-1+i0, 1) = this%F(i-1+i0, 1) + fx
          this%F(i-1+i0, 2) = this%F(i-1+i0, 2) + fy
          this%F(i-1+i0, 3) = this%F(i-1+i0, 3) + fz
          this%T(i-1+i0, 1) = this%T(i-1+i0, 1) + pQuadrupole%OY(i) * pQuadrupole%TZ(i) &
&                                    - pQuadrupole%OZ(i) * pQuadrupole%TY(i) + r1y * fz - r1z * fy
          this%T(i-1+i0, 2) = this%T(i-1+i0, 2) + pQuadrupole%OZ(i) * pQuadrupole%TX(i) &
&                                    - pQuadrupole%OX(i) * pQuadrupole%TZ(i) + r1z * fx - r1x * fz
          this%T(i-1+i0, 3) = this%T(i-1+i0, 3) + pQuadrupole%OX(i) * pQuadrupole%TY(i) &
&                                    - pQuadrupole%OY(i) * pQuadrupole%TX(i) + r1x * fy - r1y * fx
        end do
      end do

      do i = 1, l
        ! Add torques from reaction field
        tx = this%T(i-1+i0, 1) + this%tRFX(i-1+i0)
        ty = this%T(i-1+i0, 2) + this%tRFY(i-1+i0)
        tz = this%T(i-1+i0, 3) + this%tRFZ(i-1+i0)

        ! Convert torque to body-fixed coordinates
        A11 = q1(i)**2 + q2(i)**2 - q3(i)**2 - q4(i)**2
        A12 = 2._RK * (q2(i) * q3(i) + q1(i) * q4(i))
        A13 = 2._RK * (q2(i) * q4(i) - q1(i) * q3(i))
        A21 = 2._RK * (q2(i) * q3(i) - q1(i) * q4(i))
        A22 = q1(i)**2 - q2(i)**2 + q3(i)**2 - q4(i)**2
        A23 = 2._RK * (q3(i) * q4(i) + q1(i) * q2(i))
        A31 = 2._RK * (q2(i) * q4(i) + q1(i) * q3(i))
        A32 = 2._RK * (q3(i) * q4(i) - q1(i) * q2(i))
        A33 = q1(i)**2 - q2(i)**2 - q3(i)**2 + q4(i)**2
        this%T(i-1+i0, 1) = A11 * tx + A12 * ty + A13 * tz
        this%T(i-1+i0, 2) = A21 * tx + A22 * ty + A23 * tz
        this%T(i-1+i0, 3) = A31 * tx + A32 * ty + A33 * tz
      end do

    else

      ! Loop over MIEnm sites in molecule
      do j = 1, this%Molecule%NMIEnm
        pMIEnm => this%Molecule%SiteMIEnm(j)
        do i = 1, l
          this%F(i-1+i0, 1) = this%F(i-1+i0, 1) + pMIEnm%FX(i-1+i0)
          this%F(i-1+i0, 2) = this%F(i-1+i0, 2) + pMIEnm%FY(i-1+i0)
          this%F(i-1+i0, 3) = this%F(i-1+i0, 3) + pMIEnm%FZ(i-1+i0)
        end do
      end do

      ! Loop over TT68 sites in molecule
      do j = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(j)
        do i = 1, l
          this%F(i-1+i0, 1) = this%F(i-1+i0, 1) + pTT68%FX(i-1+i0)
          this%F(i-1+i0, 2) = this%F(i-1+i0, 2) + pTT68%FY(i-1+i0)
          this%F(i-1+i0, 3) = this%F(i-1+i0, 3) + pTT68%FZ(i-1+i0)
        end do
      end do

      ! Loop over charge sites in molecule
      if (LongRange .ne. RField) then
        do j = 1, this%Molecule%NCharge
          pCharge => this%Molecule%SiteCharge(j)
          do i = 1, l
            this%F(i-1+i0, 1) = this%F(i-1+i0, 1) + pCharge%FX(i)
            this%F(i-1+i0, 2) = this%F(i-1+i0, 2) + pCharge%FY(i)
            this%F(i-1+i0, 3) = this%F(i-1+i0, 3) + pCharge%FZ(i)
          end do
        end do
      end if

    end if

    ! Add forces and torques by demand
#if MPI_VER > 0
!     do i = 1, i0-1
!       if ( this%NAdd(i) ) call Atom2Mol1(this,i)
!     end do
!     do i = i1, this%NPart
!       if ( this%NAdd(i) ) call Atom2Mol1(this,i)
!     end do

    ! Reduce forces and torques from all processes
    call MPI_Reduce( this%F(:, :), this%FAll(:, :), size( this%F), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    if( this%Molecule%isElongated ) &
&     call MPI_Reduce( this%T(:, :), this%TAll(:, :), size( this%T ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
#if  TRANS == 1

! Transport  !TRANSPORT_start
    call MPI_Reduce( this%FB(:, :), this%FBAll(:, :), size( this%FB ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    call MPI_Reduce( this%FS(:, :), this%FSAll(:, :), size( this%FS ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )

  !  if (this%Conductivity) then
      call MPI_Reduce( this%FTC1(:, :), this%FTC1All(:, :), size( this%FTC1 ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( this%FTC2(:, :), this%FTC2All(:, :), size( this%FTC2 ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( this%FTC3(:, :), this%FTC3All(:, :), size( this%FTC3 ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )

      call MPI_Reduce( this%FRC1(:, :), this%FRC1All(:, :), size( this%FRC1 ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( this%FRC2(:, :), this%FRC2All(:, :), size( this%FRC2 ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( this%FRC3(:, :), this%FRC3All(:, :), size( this%FRC3 ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
   ! end if
!TRANSPORT_END
#endif
#endif

  end subroutine TComponent_Atom2Mol


!==============================================================!
!  Subroutine TComponent_Atom2Mol1                             !
!==============================================================!

  subroutine TComponent_Atom2Mol1( this, np )

    implicit none

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: np

    ! Declare local variables
    real(RK)                       :: BoxLength
    real(RK)                       :: rx, ry, rz, r1x, r1y, r1z
    real(RK)                       :: q1, q2, q3, q4
    real(RK)                       :: fx, fy, fz, tx, ty, tz
    real(RK)                       :: A11, A12, A13, A21, A22, A23, A31, A32, A33
    type(TSiteMIEnm), pointer      :: pMIEnm
    type(TSiteTT68), pointer       :: pTT68
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j


    ! Assign local variables
    BoxLength = this%BoxLength

    ! Check number of rotation axes
    if( this%Molecule%isElongated ) then


      ! Initialize local arrays
      rx = this%P0(np, 1)
      ry = this%P0(np, 2)
      rz = this%P0(np, 3)
      q1 = this%Q0(np, 1)
      q2 = this%Q0(np, 2)
      q3 = this%Q0(np, 3)
      q4 = this%Q0(np, 4)

      ! Loop over MIEnm sites in molecule
      do j = 1, this%Molecule%NMIEnm
        pMIEnm => this%Molecule%SiteMIEnm(j)
        fx = pMIEnm%FX(np)
        fy = pMIEnm%FY(np)
        fz = pMIEnm%FZ(np)
        r1x = ( pMIEnm%RX(np) - rx ) * BoxLength
        r1y = ( pMIEnm%RY(np) - ry ) * BoxLength
        r1z = ( pMIEnm%RZ(np) - rz ) * BoxLength
        this%F(np, 1) = this%F(np, 1) + fx
        this%F(np, 2) = this%F(np, 2) + fy
        this%F(np, 3) = this%F(np, 3) + fz
        this%T(np, 1) = this%T(np, 1) + r1y * fz - r1z * fy
        this%T(np, 2) = this%T(np, 2) + r1z * fx - r1x * fz
        this%T(np, 3) = this%T(np, 3) + r1x * fy - r1y * fx
      end do

      ! Loop over TT68 sites in molecule
      do j = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(j)
        fx = pTT68%FX(np)
        fy = pTT68%FY(np)
        fz = pTT68%FZ(np)
        r1x = ( pTT68%RX(np) - rx ) * BoxLength
        r1y = ( pTT68%RY(np) - ry ) * BoxLength
        r1z = ( pTT68%RZ(np) - rz ) * BoxLength
        this%F(np, 1) = this%F(np, 1) + fx
        this%F(np, 2) = this%F(np, 2) + fy
        this%F(np, 3) = this%F(np, 3) + fz
        this%T(np, 1) = this%T(np, 1) + r1y * fz - r1z * fy
        this%T(np, 2) = this%T(np, 2) + r1z * fx - r1x * fz
        this%T(np, 3) = this%T(np, 3) + r1x * fy - r1y * fx
      end do

      ! Loop over charge sites in molecule
      do j = 1, this%Molecule%NCharge
        pCharge => this%Molecule%SiteCharge(j)
        fx = pCharge%FX(np)
        fy = pCharge%FY(np)
        fz = pCharge%FZ(np)
        r1x = ( pCharge%RX(np) - rx ) * BoxLength
        r1y = ( pCharge%RY(np) - ry ) * BoxLength
        r1z = ( pCharge%RZ(np) - rz ) * BoxLength
        this%F(np, 1) = this%F(np, 1) + fx
        this%F(np, 2) = this%F(np, 2) + fy
        this%F(np, 3) = this%F(np, 3) + fz
        this%T(np, 1) = this%T(np, 1) + r1y * fz - r1z * fy
        this%T(np, 2) = this%T(np, 2) + r1z * fx - r1x * fz
        this%T(np, 3) = this%T(np, 3) + r1x * fy - r1y * fx
      end do

      ! Loop over dipole sites in molecule
      do j = 1, this%Molecule%NDipole
        pDipole => this%Molecule%SiteDipole(j)
        fx = pDipole%FX(np)
        fy = pDipole%FY(np)
        fz = pDipole%FZ(np)
        r1x = ( pDipole%RX(np) - rx ) * BoxLength
        r1y = ( pDipole%RY(np) - ry ) * BoxLength
        r1z = ( pDipole%RZ(np) - rz ) * BoxLength
        this%F(np, 1) = this%F(np, 1) + fx
        this%F(np, 2) = this%F(np, 2) + fy
        this%F(np, 3) = this%F(np, 3) + fz
        this%T(np, 1) = this%T(np, 1) + pDipole%OY(np) * pDipole%TZ(np) &
&                                - pDipole%OZ(np) * pDipole%TY(np) + r1y * fz - r1z * fy
        this%T(np, 2) = this%T(np, 2) + pDipole%OZ(np) * pDipole%TX(np) &
&                                - pDipole%OX(np) * pDipole%TZ(np) + r1z * fx - r1x * fz
        this%T(np, 3) = this%T(np, 3) + pDipole%OX(np) * pDipole%TY(np) &
&                                - pDipole%OY(np) * pDipole%TX(np) + r1x * fy - r1y * fx
      end do

      ! Loop over quadrupole sites in molecule
      do j = 1, this%Molecule%NQuadrupole
        pQuadrupole => this%Molecule%SiteQuadrupole(j)
        fx = pQuadrupole%FX(np)
        fy = pQuadrupole%FY(np)
        fz = pQuadrupole%FZ(np)
        r1x = ( pQuadrupole%RX(np) - rx ) * BoxLength
        r1y = ( pQuadrupole%RY(np) - ry ) * BoxLength
        r1z = ( pQuadrupole%RZ(np) - rz ) * BoxLength
        this%F(np, 1) = this%F(np, 1) + fx
        this%F(np, 2) = this%F(np, 2) + fy
        this%F(np, 3) = this%F(np, 3) + fz
        this%T(np, 1) = this%T(np, 1) + pQuadrupole%OY(np) * pQuadrupole%TZ(np) &
&                                - pQuadrupole%OZ(np) * pQuadrupole%TY(np) + r1y * fz - r1z * fy
        this%T(np, 2) = this%T(np, 2) + pQuadrupole%OZ(np) * pQuadrupole%TX(np) &
&                                - pQuadrupole%OX(np) * pQuadrupole%TZ(np) + r1z * fx - r1x * fz
        this%T(np, 3) = this%T(np, 3) + pQuadrupole%OX(np) * pQuadrupole%TY(np) &
&                                - pQuadrupole%OY(np) * pQuadrupole%TX(np) + r1x * fy - r1y * fx
      end do

      ! Add torques from reaction field
      tx = this%T(np, 1) + this%tRFX(np)
      ty = this%T(np, 2) + this%tRFY(np)
      tz = this%T(np, 3) + this%tRFZ(np)

      ! Convert torque to body-fixed coordinates
      A11 = q1**2 + q2**2 - q3**2 - q4**2
      A12 = 2._RK * (q2 * q3 + q1 * q4)
      A13 = 2._RK * (q2 * q4 - q1 * q3)
      A21 = 2._RK * (q2 * q3 - q1 * q4)
      A22 = q1**2 - q2**2 + q3**2 - q4**2
      A23 = 2._RK * (q3 * q4 + q1 * q2)
      A31 = 2._RK * (q2 * q4 + q1 * q3)
      A32 = 2._RK * (q3 * q4 - q1 * q2)
      A33 = q1**2 - q2**2 - q3**2 + q4**2
      this%T(np, 1) = A11 * tx + A12 * ty + A13 * tz
      this%T(np, 2) = A21 * tx + A22 * ty + A23 * tz
      this%T(np, 3) = A31 * tx + A32 * ty + A33 * tz

    else

      ! Loop over MIEnm sites in molecule
      do j = 1, this%Molecule%NMIEnm
        pMIEnm => this%Molecule%SiteMIEnm(j)
        this%F(np, 1) = this%F(np, 1) + pMIEnm%FX(np)
        this%F(np, 2) = this%F(np, 2) + pMIEnm%FY(np)
        this%F(np, 3) = this%F(np, 3) + pMIEnm%FZ(np)
      end do

      ! Loop over TT68 sites in molecule
      do j = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(j)
        this%F(np, 1) = this%F(np, 1) + pTT68%FX(np)
        this%F(np, 2) = this%F(np, 2) + pTT68%FY(np)
        this%F(np, 3) = this%F(np, 3) + pTT68%FZ(np)
      end do

      ! Loop over charge sites in molecule
      if (LongRange .ne. RField) then
        do j = 1, this%Molecule%NCharge
          pCharge => this%Molecule%SiteCharge(j)
          this%F(np, 1) = this%F(np, 1) + pCharge%FX(np)
          this%F(np, 2) = this%F(np, 2) + pCharge%FY(np)
          this%F(np, 3) = this%F(np, 3) + pCharge%FZ(np)
        end do
      end if

    end if

  end subroutine TComponent_Atom2Mol1


!==============================================================!
!  Subroutine TComponent_Atom2Mol_Trans                        !
!==============================================================!

  subroutine TComponent_Atom2Mol_Trans( this, np )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: np

    ! Declare local variables
    real(RK)                       :: BoxLength
    real(RK)                       :: rx(np), ry(np), rz(np), r1x, r1y, r1z
    real(RK)                       :: q1(np), q2(np), q3(np), q4(np)
    real(RK)                       :: fx, fy, fz, tx, ty, tz
    real(RK)                       :: A11, A12, A13, A21, A22, A23, A31, A32, A33
    type(TSiteMIEnm), pointer      :: pMIEnm
    type(TSiteTT68), pointer       :: pTT68
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j

#if  TRANS == 1
    !TRANSPORT_start
    real(RK)                       :: vsx,vsy,vsz
    real(RK)                       :: vsux,vsuy,vsuz
    real(RK)                       :: vbx,vby,vbz
    real(RK)                       :: cx, cy, cz
    real(RK)                       :: tux, tuy, tuz, tlx, tly, tlz, tdx, tdy, tdz
    real(RK)                       :: BoxLength_dt
#endif
!TRANSPORT_END

    ! Assign local variables
    BoxLength = this%BoxLength
#if  TRANS == 1
    !TRANSPORT_start
    BoxLength_dt = this%BoxLength/TimeStep !TRANSPORT_thisline
    this%FS(:,:) = 0._RK
    this%FB(:,:) = 0._RK
  !  if (this%Conductivity) then
      this%FTC(:,:) = 0._RK
      this%FRC(:,:) = 0._RK
      this%FTC1(:,:) = 0._RK
      this%FTC2(:,:) = 0._RK
      this%FTC3(:,:) = 0._RK
      this%FRC1(:,:) = 0._RK
      this%FRC2(:,:) = 0._RK
      this%FRC3(:,:) = 0._RK
   ! end if
    !TRANSPORT_END
#endif

    ! Initialize forces
    this%F(1:np, :) = 0._RK

    ! Check number of rotation axes
    if( this%Molecule%isElongated ) then

      ! Initialize torques
      this%T(1:np, :) = 0._RK

      ! Initialize local arrays
      rx(:) = this%P0(:, 1)
      ry(:) = this%P0(:, 2)
      rz(:) = this%P0(:, 3)
      q1(:) = this%Q0(:, 1)
      q2(:) = this%Q0(:, 2)
      q3(:) = this%Q0(:, 3)
      q4(:) = this%Q0(:, 4)

      ! Loop over MIEnm sites in molecule
      do j = 1, this%Molecule%NMIEnm
        pMIEnm => this%Molecule%SiteMIEnm(j)
        do i = 1, np
          fx = pMIEnm%FX(i)
          fy = pMIEnm%FY(i)
          fz = pMIEnm%FZ(i)
#if  TRANS == 1
          !TRANSPORT_start
          vsx = pMIEnm%vsMIEx(i)
          vsy = pMIEnm%vsMIEy(i)
          vsz = pMIEnm%vsMIEz(i)
          vbx = pMIEnm%vbMIEx(i)
          vby = pMIEnm%vbMIEy(i)
          vbz = pMIEnm%vbMIEz(i)
       !   if (this%Conductivity) then
            vsux= pMIEnm%vsuMIEx(i)
            vsuy= pMIEnm%vsuMIEy(i)
            vsuz= pMIEnm%vsuMIEz(i)
            cx  = pMIEnm%cMIEx(i)
            cy  = pMIEnm%cMIEy(i)
            cz  = pMIEnm%cMIEz(i)
            tux = pMIEnm%tuMIEx(i)
            tuy = pMIEnm%tuMIEy(i)
            tuz = pMIEnm%tuMIEz(i)
            tlx = pMIEnm%tlMIEx(i)
            tly = pMIEnm%tlMIEy(i)
            tlz = pMIEnm%tlMIEz(i)
            tdx = pMIEnm%tdMIEx(i)
            tdy = pMIEnm%tdMIEy(i)
            tdz = pMIEnm%tdMIEz(i)
       !   end if
          !TRANSPORT_END
#endif
          r1x = ( pMIEnm%RX(i) - rx(i) ) * BoxLength
          r1y = ( pMIEnm%RY(i) - ry(i) ) * BoxLength
          r1z = ( pMIEnm%RZ(i) - rz(i) ) * BoxLength
          this%F(i, 1) = this%F(i, 1) + fx
          this%F(i, 2) = this%F(i, 2) + fy
          this%F(i, 3) = this%F(i, 3) + fz
          this%T(i, 1) = this%T(i, 1) + r1y * fz - r1z * fy
          this%T(i, 2) = this%T(i, 2) + r1z * fx - r1x * fz
          this%T(i, 3) = this%T(i, 3) + r1x * fy - r1y * fx
#if  TRANS == 1
          !TRANSPORT_start
          this%FS(i, 1)= this%FS(i, 1)+ vsx
          this%FS(i, 2)= this%FS(i, 2)+ vsy
          this%FS(i, 3)= this%FS(i, 3)+ vsz
          this%FB(i, 1)= this%FB(i, 1)+ vbx
          this%FB(i, 2)= this%FB(i, 2)+ vby
          this%FB(i, 3)= this%FB(i, 3)+ vbz

         ! if (this%Conductivity) then
            this%FTC1(i, 1)= this%FTC1(i, 1) +(cx+vbx)
            this%FTC1(i, 2)= this%FTC1(i, 2) + vsux
            this%FTC1(i, 3)= this%FTC1(i, 3) + vsuy
            this%FTC2(i, 1)= this%FTC2(i, 1) + vsx
            this%FTC2(i, 2)= this%FTC2(i, 2) +(cy+vby)
            this%FTC2(i, 3)= this%FTC2(i, 3) + vsuz
            this%FTC3(i, 1)= this%FTC3(i, 1) + vsy
            this%FTC3(i, 2)= this%FTC3(i, 2) + vsz
            this%FTC3(i, 3)= this%FTC3(i, 3) +(cz+vbz)

            this%FRC1(i,1) = this%FRC1(i,1) + tdx
            this%FRC1(i,2) = this%FRC1(i,2) + tux
            this%FRC1(i,3) = this%FRC1(i,3) + tuy
            this%FRC2(i,1) = this%FRC2(i,1) + tlx
            this%FRC2(i,2) = this%FRC2(i,2) + tdy
            this%FRC2(i,3) = this%FRC2(i,3) + tuz
            this%FRC3(i,1) = this%FRC3(i,1) + tly
            this%FRC3(i,2) = this%FRC3(i,2) + tlz
            this%FRC3(i,3) = this%FRC3(i,3) + tdz
         ! end if
           !TRANSPORT_END
#endif
        end do
      end do

      ! Loop over TT68 sites in molecule
      do j = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(j)
        do i = 1, np
          fx = pTT68%FX(i)
          fy = pTT68%FY(i)
          fz = pTT68%FZ(i)
#if  TRANS == 1
          !TRANSPORT_start
          vsx = pTT68%vsTTx(i)
          vsy = pTT68%vsTTy(i)
          vsz = pTT68%vsTTz(i)
          vbx = pTT68%vbTTx(i)
          vby = pTT68%vbTTy(i)
          vbz = pTT68%vbTTz(i)
       !   if (this%Conductivity) then
            vsux= pTT68%vsuTTx(i)
            vsuy= pTT68%vsuTTy(i)
            vsuz= pTT68%vsuTTz(i)
            cx  = pTT68%cTTx(i)
            cy  = pTT68%cTTy(i)
            cz  = pTT68%cTTz(i)
            tux = pTT68%tuTTx(i)
            tuy = pTT68%tuTTy(i)
            tuz = pTT68%tuTTz(i)
            tlx = pTT68%tlTTx(i)
            tly = pTT68%tlTTy(i)
            tlz = pTT68%tlTTz(i)
            tdx = pTT68%tdTTx(i)
            tdy = pTT68%tdTTy(i)
            tdz = pTT68%tdTTz(i)
       !   end if
          !TRANSPORT_END
#endif
          r1x = ( pTT68%RX(i) - rx(i) ) * BoxLength
          r1y = ( pTT68%RY(i) - ry(i) ) * BoxLength
          r1z = ( pTT68%RZ(i) - rz(i) ) * BoxLength
          this%F(i, 1) = this%F(i, 1) + fx
          this%F(i, 2) = this%F(i, 2) + fy
          this%F(i, 3) = this%F(i, 3) + fz
          this%T(i, 1) = this%T(i, 1) + r1y * fz - r1z * fy
          this%T(i, 2) = this%T(i, 2) + r1z * fx - r1x * fz
          this%T(i, 3) = this%T(i, 3) + r1x * fy - r1y * fx
#if  TRANS == 1
          !TRANSPORT_start
          this%FS(i, 1)= this%FS(i, 1)+ vsx
          this%FS(i, 2)= this%FS(i, 2)+ vsy
          this%FS(i, 3)= this%FS(i, 3)+ vsz
          this%FB(i, 1)= this%FB(i, 1)+ vbx
          this%FB(i, 2)= this%FB(i, 2)+ vby
          this%FB(i, 3)= this%FB(i, 3)+ vbz

         ! if (this%Conductivity) then
            this%FTC1(i, 1)= this%FTC1(i, 1) +(cx+vbx)
            this%FTC1(i, 2)= this%FTC1(i, 2) + vsux
            this%FTC1(i, 3)= this%FTC1(i, 3) + vsuy
            this%FTC2(i, 1)= this%FTC2(i, 1) + vsx
            this%FTC2(i, 2)= this%FTC2(i, 2) +(cy+vby)
            this%FTC2(i, 3)= this%FTC2(i, 3) + vsuz
            this%FTC3(i, 1)= this%FTC3(i, 1) + vsy
            this%FTC3(i, 2)= this%FTC3(i, 2) + vsz
            this%FTC3(i, 3)= this%FTC3(i, 3) +(cz+vbz)

            this%FRC1(i,1) = this%FRC1(i,1) + tdx
            this%FRC1(i,2) = this%FRC1(i,2) + tux
            this%FRC1(i,3) = this%FRC1(i,3) + tuy
            this%FRC2(i,1) = this%FRC2(i,1) + tlx
            this%FRC2(i,2) = this%FRC2(i,2) + tdy
            this%FRC2(i,3) = this%FRC2(i,3) + tuz
            this%FRC3(i,1) = this%FRC3(i,1) + tly
            this%FRC3(i,2) = this%FRC3(i,2) + tlz
            this%FRC3(i,3) = this%FRC3(i,3) + tdz
         ! end if
           !TRANSPORT_END
#endif
        end do
      end do

      ! Loop over charge sites in molecule
      do j = 1, this%Molecule%NCharge
        pCharge => this%Molecule%SiteCharge(j)
        do i = 1, np
          fx = pCharge%FX(i)
          fy = pCharge%FY(i)
          fz = pCharge%FZ(i)
#if  TRANS == 1
          !TRANSPORT_start
          vsx = pCharge%vsCx(i)
          vsy = pCharge%vsCy(i)
          vsz = pCharge%vsCz(i)
          vbx = pCharge%vbCx(i)
          vby = pCharge%vbCy(i)
          vbz = pCharge%vbCz(i)
      !    if (this%Conductivity) then
            vsux= pCharge%vsuCx(i)
            vsuy= pCharge%vsuCy(i)
            vsuz= pCharge%vsuCz(i)
            cx  = pCharge%cCx(i)
            cy  = pCharge%cCy(i)
            cz  = pCharge%cCz(i)
            tux = pCharge%tuCx(i)
            tuy = pCharge%tuCy(i)
            tuz = pCharge%tuCz(i)
            tlx = pCharge%tlCx(i)
            tly = pCharge%tlCy(i)
            tlz = pCharge%tlCz(i)
            tdx = pCharge%tdCx(i)
            tdy = pCharge%tdCy(i)
            tdz = pCharge%tdCz(i)
       !   end if
          !TRANSPORT_END
#endif
          r1x = ( pCharge%RX(i) - rx(i) ) * BoxLength
          r1y = ( pCharge%RY(i) - ry(i) ) * BoxLength
          r1z = ( pCharge%RZ(i) - rz(i) ) * BoxLength
          this%F(i, 1) = this%F(i, 1) + fx
          this%F(i, 2) = this%F(i, 2) + fy
          this%F(i, 3) = this%F(i, 3) + fz
          this%T(i, 1) = this%T(i, 1) + r1y * fz - r1z * fy
          this%T(i, 2) = this%T(i, 2) + r1z * fx - r1x * fz
          this%T(i, 3) = this%T(i, 3) + r1x * fy - r1y * fx
#if  TRANS == 1
          !TRANSPORT_start
          this%FS(i, 1)= this%FS(i, 1)+ vsx
          this%FS(i, 2)= this%FS(i, 2)+ vsy
          this%FS(i, 3)= this%FS(i, 3)+ vsz
          this%FB(i, 1)= this%FB(i, 1)+ vbx
          this%FB(i, 2)= this%FB(i, 2)+ vby
          this%FB(i, 3)= this%FB(i, 3)+ vbz

        !  if (this%Conductivity) then
            this%FTC1(i, 1)= this%FTC1(i, 1) +(cx+vbx)
            this%FTC1(i, 2)= this%FTC1(i, 2) + vsux
            this%FTC1(i, 3)= this%FTC1(i, 3) + vsuy
            this%FTC2(i, 1)= this%FTC2(i, 1) + vsx
            this%FTC2(i, 2)= this%FTC2(i, 2) +(cy+vby)
            this%FTC2(i, 3)= this%FTC2(i, 3) + vsuz
            this%FTC3(i, 1)= this%FTC3(i, 1) + vsy
            this%FTC3(i, 2)= this%FTC3(i, 2) + vsz
            this%FTC3(i, 3)= this%FTC3(i, 3) +(cz+vbz)

            this%FRC1(i,1) = this%FRC1(i,1) + tdx
            this%FRC1(i,2) = this%FRC1(i,2) + tux
            this%FRC1(i,3) = this%FRC1(i,3) + tuy
            this%FRC2(i,1) = this%FRC2(i,1) + tlx
            this%FRC2(i,2) = this%FRC2(i,2) + tdy
            this%FRC2(i,3) = this%FRC2(i,3) + tuz
            this%FRC3(i,1) = this%FRC3(i,1) + tly
            this%FRC3(i,2) = this%FRC3(i,2) + tlz
            this%FRC3(i,3) = this%FRC3(i,3) + tdz
         ! end if
          !TRANSPORT_END
#endif
        end do
      end do

      ! Loop over dipole sites in molecule
      do j = 1, this%Molecule%NDipole
        pDipole => this%Molecule%SiteDipole(j)
        do i = 1, np
          fx = pDipole%FX(i)
          fy = pDipole%FY(i)
          fz = pDipole%FZ(i)
#if  TRANS == 1
          !TRANSPORT_start
          vsx = pDipole%vsDx(i)
          vsy = pDipole%vsDy(i)
          vsz = pDipole%vsDz(i)
          vbx = pDipole%vbDx(i)
          vby = pDipole%vbDy(i)
          vbz = pDipole%vbDz(i)
         ! if (this%Conductivity) then
            vsux= pDipole%vsuDx(i)
            vsuy= pDipole%vsuDy(i)
            vsuz= pDipole%vsuDz(i)
            cx  = pDipole%cDx(i)
            cy  = pDipole%cDy(i)
            cz  = pDipole%cDz(i)
            tux = pDipole%tuDx(i)
            tuy = pDipole%tuDy(i)
            tuz = pDipole%tuDz(i)
            tlx = pDipole%tlDx(i)
            tly = pDipole%tlDy(i)
            tlz = pDipole%tlDz(i)
            tdx = pDipole%tdDx(i)
            tdy = pDipole%tdDy(i)
            tdz = pDipole%tdDz(i)
        !  end if
          !TRANSPORT_END
#endif
          r1x = ( pDipole%RX(i) - rx(i) ) * BoxLength
          r1y = ( pDipole%RY(i) - ry(i) ) * BoxLength
          r1z = ( pDipole%RZ(i) - rz(i) ) * BoxLength
          this%F(i, 1) = this%F(i, 1) + fx
          this%F(i, 2) = this%F(i, 2) + fy
          this%F(i, 3) = this%F(i, 3) + fz
          this%T(i, 1) = this%T(i, 1) + pDipole%OY(i) * pDipole%TZ(i) &
&                                     - pDipole%OZ(i) * pDipole%TY(i) &
&                                     + r1y * fz - r1z * fy
          this%T(i, 2) = this%T(i, 2) + pDipole%OZ(i) * pDipole%TX(i) &
&                                     - pDipole%OX(i) * pDipole%TZ(i) &
&                                     + r1z * fx - r1x * fz
          this%T(i, 3) = this%T(i, 3) + pDipole%OX(i) * pDipole%TY(i) &
&                                     - pDipole%OY(i) * pDipole%TX(i) &
&                                     + r1x * fy - r1y * fx
#if  TRANS == 1
         !TRANSPORT_start
          this%FS(i, 1)= this%FS(i, 1)+ vsx
          this%FS(i, 2)= this%FS(i, 2)+ vsy
          this%FS(i, 3)= this%FS(i, 3)+ vsz
          this%FB(i, 1)= this%FB(i, 1)+ vbx
          this%FB(i, 2)= this%FB(i, 2)+ vby
          this%FB(i, 3)= this%FB(i, 3)+ vbz

       !   if (this%Conductivity) then
            this%FTC1(i, 1)= this%FTC1(i, 1) +(cx+vbx)
            this%FTC1(i, 2)= this%FTC1(i, 2) + vsux
            this%FTC1(i, 3)= this%FTC1(i, 3) + vsuy
            this%FTC2(i, 1)= this%FTC2(i, 1) + vsx
            this%FTC2(i, 2)= this%FTC2(i, 2) +(cy+vby)
            this%FTC2(i, 3)= this%FTC2(i, 3) + vsuz
            this%FTC3(i, 1)= this%FTC3(i, 1) + vsy
            this%FTC3(i, 2)= this%FTC3(i, 2) + vsz
            this%FTC3(i, 3)= this%FTC3(i, 3) +(cz+vbz)

            this%FRC1(i,1) = this%FRC1(i,1) + tdx
            this%FRC1(i,2) = this%FRC1(i,2) + tux
            this%FRC1(i,3) = this%FRC1(i,3) + tuy
            this%FRC2(i,1) = this%FRC2(i,1) + tlx
            this%FRC2(i,2) = this%FRC2(i,2) + tdy
            this%FRC2(i,3) = this%FRC2(i,3) + tuz
            this%FRC3(i,1) = this%FRC3(i,1) + tly
            this%FRC3(i,2) = this%FRC3(i,2) + tlz
            this%FRC3(i,3) = this%FRC3(i,3) + tdz
        !  end if
!TRANSPORT_END
#endif
        end do
      end do

      ! Loop over quadrupole sites in molecule
      do j = 1, this%Molecule%NQuadrupole
        pQuadrupole => this%Molecule%SiteQuadrupole(j)
        do i = 1, np
          fx = pQuadrupole%FX(i)
          fy = pQuadrupole%FY(i)
          fz = pQuadrupole%FZ(i)
#if  TRANS == 1
          !TRANSPORT_start
          vsx = pQuadrupole%vsQx(i)
          vsy = pQuadrupole%vsQy(i)
          vsz = pQuadrupole%vsQz(i)
          vbx = pQuadrupole%vbQx(i)
          vby = pQuadrupole%vbQy(i)
          vbz = pQuadrupole%vbQz(i)
     !     if (this%Conductivity) then
            vsux= pQuadrupole%vsuQx(i)
            vsuy= pQuadrupole%vsuQy(i)
            vsuz= pQuadrupole%vsuQz(i)
            cx  = pQuadrupole%cQx(i)
            cy  = pQuadrupole%cQy(i)
            cz  = pQuadrupole%cqz(i)
            tux = pQuadrupole%tuQx(i)
            tuy = pQuadrupole%tuQy(i)
            tuz = pQuadrupole%tuQz(i)
            tlx = pQuadrupole%tlQx(i)
            tly = pQuadrupole%tlQy(i)
            tlz = pQuadrupole%tlQz(i)
            tdx = pQuadrupole%tdQx(i)
            tdy = pQuadrupole%tdQy(i)
            tdz = pQuadrupole%tdQz(i)
      !    end if
          !TRANSPORT_END
#endif
          r1x = ( pQuadrupole%RX(i) - rx(i) ) * BoxLength
          r1y = ( pQuadrupole%RY(i) - ry(i) ) * BoxLength
          r1z = ( pQuadrupole%RZ(i) - rz(i) ) * BoxLength
          this%F(i, 1) = this%F(i, 1) + fx
          this%F(i, 2) = this%F(i, 2) + fy
          this%F(i, 3) = this%F(i, 3) + fz
          this%T(i, 1) = this%T(i, 1) + pQuadrupole%OY(i) * pQuadrupole%TZ(i) &
&                                     - pQuadrupole%OZ(i) * pQuadrupole%TY(i) &
&                                     + r1y * fz - r1z * fy
          this%T(i, 2) = this%T(i, 2) + pQuadrupole%OZ(i) * pQuadrupole%TX(i) &
&                                     - pQuadrupole%OX(i) * pQuadrupole%TZ(i) &
&                                     + r1z * fx - r1x * fz
          this%T(i, 3) = this%T(i, 3) + pQuadrupole%OX(i) * pQuadrupole%TY(i) &
&                                     - pQuadrupole%OY(i) * pQuadrupole%TX(i) &
&                                     + r1x * fy - r1y * fx
#if  TRANS == 1
  !TRANSPORT_start
          this%FS(i, 1)= this%FS(i, 1)+ vsx
          this%FS(i, 2)= this%FS(i, 2)+ vsy
          this%FS(i, 3)= this%FS(i, 3)+ vsz
          this%FB(i, 1)= this%FB(i, 1)+ vbx
          this%FB(i, 2)= this%FB(i, 2)+ vby
          this%FB(i, 3)= this%FB(i, 3)+ vbz

      !    if (this%Conductivity) then
            this%FTC1(i, 1)= this%FTC1(i, 1) +(cx+vbx)
            this%FTC1(i, 2)= this%FTC1(i, 2) + vsux
            this%FTC1(i, 3)= this%FTC1(i, 3) + vsuy
            this%FTC2(i, 1)= this%FTC2(i, 1) + vsx
            this%FTC2(i, 2)= this%FTC2(i, 2) +(cy+vby)
            this%FTC2(i, 3)= this%FTC2(i, 3) + vsuz
            this%FTC3(i, 1)= this%FTC3(i, 1) + vsy
            this%FTC3(i, 2)= this%FTC3(i, 2) + vsz
            this%FTC3(i, 3)= this%FTC3(i, 3) +(cz+vbz)

            this%FRC1(i,1) = this%FRC1(i,1) + tdx
            this%FRC1(i,2) = this%FRC1(i,2) + tux
            this%FRC1(i,3) = this%FRC1(i,3) + tuy
            this%FRC2(i,1) = this%FRC2(i,1) + tlx
            this%FRC2(i,2) = this%FRC2(i,2) + tdy
            this%FRC2(i,3) = this%FRC2(i,3) + tuz
            this%FRC3(i,1) = this%FRC3(i,1) + tly
            this%FRC3(i,2) = this%FRC3(i,2) + tlz
            this%FRC3(i,3) = this%FRC3(i,3) + tdz
       !   end if
!TRANSPORT_END
#endif
        end do
      end do

      do i = 1, np
        ! Add torques from reaction field
        tx = this%T(i, 1) + this%tRFX(i)
        ty = this%T(i, 2) + this%tRFY(i)
        tz = this%T(i, 3) + this%tRFZ(i)

        ! Convert torque to body-fixed coordinates
        A11 = q1(i)**2 + q2(i)**2 - q3(i)**2 - q4(i)**2
        A12 = 2._RK * (q2(i) * q3(i) + q1(i) * q4(i))
        A13 = 2._RK * (q2(i) * q4(i) - q1(i) * q3(i))
        A21 = 2._RK * (q2(i) * q3(i) - q1(i) * q4(i))
        A22 = q1(i)**2 - q2(i)**2 + q3(i)**2 - q4(i)**2
        A23 = 2._RK * (q3(i) * q4(i) + q1(i) * q2(i))
        A31 = 2._RK * (q2(i) * q4(i) + q1(i) * q3(i))
        A32 = 2._RK * (q3(i) * q4(i) - q1(i) * q2(i))
        A33 = q1(i)**2 - q2(i)**2 - q3(i)**2 + q4(i)**2
        this%T(i, 1) = A11 * tx + A12 * ty + A13 * tz
        this%T(i, 2) = A21 * tx + A22 * ty + A23 * tz
        this%T(i, 3) = A31 * tx + A32 * ty + A33 * tz
      end do

    else

      ! Loop over MIEnm sites in molecule
      do j = 1, this%Molecule%NMIEnm
        pMIEnm => this%Molecule%SiteMIEnm(j)
        do i = 1, np
#if  TRANS == 1
        !TRANSPORT_start
          vsx = pMIEnm%vsMIEx(i)
          vsy = pMIEnm%vsMIEy(i)
          vsz = pMIEnm%vsMIEz(i)
          vbx = pMIEnm%vbMIEx(i)
          vby = pMIEnm%vbMIEy(i)
          vbz = pMIEnm%vbMIEz(i)
     !     if (this%Conductivity) then
            vsux= pMIEnm%vsuMIEx(i)
            vsuy= pMIEnm%vsuMIEy(i)
            vsuz= pMIEnm%vsuMIEz(i)
            cx  = pMIEnm%cMIEx(i)
            cy  = pMIEnm%cMIEy(i)
            cz  = pMIEnm%cMIEz(i)
      !    end if
          !TRANSPORT_END
#endif
          this%F(i, 1) = this%F(i, 1) + pMIEnm%FX(i)
          this%F(i, 2) = this%F(i, 2) + pMIEnm%FY(i)
          this%F(i, 3) = this%F(i, 3) + pMIEnm%FZ(i)
#if  TRANS == 1
          !TRANSPORT_start
          this%FS(i, 1) = this%FS(i, 1) + vsx
          this%FS(i, 2) = this%FS(i, 2) + vsy
          this%FS(i, 3) = this%FS(i, 3) + vsz
          this%FB(i, 1) = this%FB(i, 1) + vbx
          this%FB(i, 2) = this%FB(i, 2) + vby
          this%FB(i, 3) = this%FB(i, 3) + vbz

       !   if (this%Conductivity) then
            this%FTC1(i, 1)= this%FTC1(i, 1) +(cx+vbx)
            this%FTC1(i, 2)= this%FTC1(i, 2) + vsux
            this%FTC1(i, 3)= this%FTC1(i, 3) + vsuy
            this%FTC2(i, 1)= this%FTC2(i, 1) + vsx
            this%FTC2(i, 2)= this%FTC2(i, 2) +(cy+vby)
            this%FTC2(i, 3)= this%FTC2(i, 3) + vsuz
            this%FTC3(i, 1)= this%FTC3(i, 1) + vsy
            this%FTC3(i, 2)= this%FTC3(i, 2) + vsz
            this%FTC3(i, 3)= this%FTC3(i, 3) +(cz+vbz)
        !  end if
!TRANSPORT_END
#endif
        end do
      end do

      ! Loop over TT68 sites in molecule
      do j = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(j)
        do i = 1, np
#if  TRANS == 1
        !TRANSPORT_start
          vsx = pTT68%vsTTx(i)
          vsy = pTT68%vsTTy(i)
          vsz = pTT68%vsTTz(i)
          vbx = pTT68%vbTTx(i)
          vby = pTT68%vbTTy(i)
          vbz = pTT68%vbTTz(i)
     !     if (this%Conductivity) then
            vsux= pTT68%vsuTTx(i)
            vsuy= pTT68%vsuTTy(i)
            vsuz= pTT68%vsuTTz(i)
            cx  = pTT68%cTTx(i)
            cy  = pTT68%cTTy(i)
            cz  = pTT68%cTTz(i)
      !    end if
          !TRANSPORT_END
#endif
          this%F(i, 1) = this%F(i, 1) + pTT68%FX(i)
          this%F(i, 2) = this%F(i, 2) + pTT68%FY(i)
          this%F(i, 3) = this%F(i, 3) + pTT68%FZ(i)
#if  TRANS == 1
          !TRANSPORT_start
          this%FS(i, 1) = this%FS(i, 1) + vsx
          this%FS(i, 2) = this%FS(i, 2) + vsy
          this%FS(i, 3) = this%FS(i, 3) + vsz
          this%FB(i, 1) = this%FB(i, 1) + vbx
          this%FB(i, 2) = this%FB(i, 2) + vby
          this%FB(i, 3) = this%FB(i, 3) + vbz

       !   if (this%Conductivity) then
            this%FTC1(i, 1)= this%FTC1(i, 1) +(cx+vbx)
            this%FTC1(i, 2)= this%FTC1(i, 2) + vsux
            this%FTC1(i, 3)= this%FTC1(i, 3) + vsuy
            this%FTC2(i, 1)= this%FTC2(i, 1) + vsx
            this%FTC2(i, 2)= this%FTC2(i, 2) +(cy+vby)
            this%FTC2(i, 3)= this%FTC2(i, 3) + vsuz
            this%FTC3(i, 1)= this%FTC3(i, 1) + vsy
            this%FTC3(i, 2)= this%FTC3(i, 2) + vsz
            this%FTC3(i, 3)= this%FTC3(i, 3) +(cz+vbz)
        !  end if
!TRANSPORT_END
#endif
        end do
      end do

      ! Loop over charge sites in molecule
      if (LongRange .ne. RField) then
        do j = 1, this%Molecule%NCharge
          pCharge => this%Molecule%SiteCharge(j)
          do i = 1, np
            this%F(i, 1) = this%F(i, 1) + pCharge%FX(i)
            this%F(i, 2) = this%F(i, 2) + pCharge%FY(i)
            this%F(i, 3) = this%F(i, 3) + pCharge%FZ(i)
#if  TRANS == 1
            !TRANSPORT_start
            vsx = pCharge%vsCx(i)
            vsy = pCharge%vsCy(i)
            vsz = pCharge%vsCz(i)
            vbx = pCharge%vbCx(i)
            vby = pCharge%vbCy(i)
            vbz = pCharge%vbCz(i)
            this%FS(i, 1)= this%FS(i, 1)+ vsx
            this%FS(i, 2)= this%FS(i, 2)+ vsy
            this%FS(i, 3)= this%FS(i, 3)+ vsz
            this%FB(i, 1)= this%FB(i, 1)+ vbx
            this%FB(i, 2)= this%FB(i, 2)+ vby
            this%FB(i, 3)= this%FB(i, 3)+ vbz
            !TRANSPORT_END
#endif
          end do
        end do
      end if

    end if

    ! Reduce forces and torques from all processes
#if MPI_VER > 0
    call MPI_Reduce( this%F(:, :), this%FAll(:, :), size( this%F ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    if( this%Molecule%isElongated ) call MPI_Reduce( this%T(:, :), this%TAll(:, :), size( this%T ), &
&     MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )

#if  TRANS == 1
! Transport  !TRANSPORT_start
    call MPI_Reduce( this%FB(:, :), this%FBAll(:, :), size( this%FB ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    call MPI_Reduce( this%FS(:, :), this%FSAll(:, :), size( this%FS ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )

 !   if (this%Conductivity) then
      call MPI_Reduce( this%FTC1(:, :), this%FTC1All(:, :), size( this%FTC1 ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( this%FTC2(:, :), this%FTC2All(:, :), size( this%FTC2 ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( this%FTC3(:, :), this%FTC3All(:, :), size( this%FTC3 ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )

      call MPI_Reduce( this%FRC1(:, :), this%FRC1All(:, :), size( this%FRC1 ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( this%FRC2(:, :), this%FRC2All(:, :), size( this%FRC2 ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( this%FRC3(:, :), this%FRC3All(:, :), size( this%FRC3 ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
  !  end if
!TRANSPORT_END
#endif
#endif

  end subroutine TComponent_Atom2Mol_Trans


#if OSMOP > 0
!==============================================================!
!  Subroutine TComponent_DensityProfile                        !
!==============================================================!

  subroutine TComponent_DensityProfile( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent)    :: this

    ! Declare local variables
    integer             :: i, j
#if MPI_VER > 0
    integer,allocatable :: DensityN(:)

    allocate( DensityN(NBinsDen) )
    DensityN(:) = 0
#endif

    ! Initialize local arrays
    this%DensityProfileN(:) = 0

    ! Loop over MIEnm or TT68 sites in molecule
#if MPI_VER > 0
loop1:do i = this%NPart0, this%NPart2
#else
loop1:do i = 1, this%NPart
#endif
      do j = 1, NBinsDen
        if (this%P0(i,1) .ge. real(j-1)/NBinsDen-.5_RK) then
          if (this%P0(i,1) < real(j)/NBinsDen-.5_RK) then
            this%DensityProfileN(j) = this%DensityProfileN(j)+1
            cycle loop1
          end if
        end if
      end do
    end do loop1

#if MPI_VER > 0
    call MPI_Reduce( this%DensityProfileN(:), DensityN(:), size(this%DensityProfileN), MPI_INTEGER, MPI_SUM, NRootProc, Communicator, ierror )
    if (RootProc) this%DensityProfileN(:) = DensityN(:)
#endif

  end subroutine TComponent_DensityProfile
#endif


!==============================================================!
!  Subroutine TComponent_PredictGear                           !
!==============================================================!

  subroutine TComponent_PredictGear( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Declare local variables
    integer :: np, nra
    integer :: i, j

    ! Assign local variables
    np = this%NPart
    nra = this%Molecule%NDFRot

    ! Predict COM positions and their derivatives
    do j = 1, 3
      do i = 1, np
        this%P0(i, j) = this%P0(i, j) + this%P1(i, j) + this%P2(i, j) + this%P3(i, j) + this%P4(i, j) + this%P5(i, j)
        this%P1(i, j) = this%P1(i, j) + 2._RK * this%P2(i, j) + 3._RK * this%P3(i, j) + 4._RK * this%P4(i, j) + 5._RK * this%P5(i, j)
        this%P2(i, j) = this%P2(i, j) + 3._RK * this%P3(i, j) + 6._RK * this%P4(i, j) + 10._RK * this%P5(i, j)
        this%P3(i, j) = this%P3(i, j) + 4._RK * this%P4(i, j) + 10._RK * this%P5(i, j)
        this%P4(i, j) = this%P4(i, j) + 5._RK * this%P5(i, j)

        ! Calculate displacement
        this%Disp(i, j) = this%Disp(i, j) + this%P0(i, j) - this%P0old(i, j)

        ! Check for conservation of particles in primary cell
#if ARCH == 1
        if( this%P0(i, j) < -.5_RK ) then
          this%P0(i, j) = this%P0(i, j) + 1._RK
        elseif( this%P0(i, j) > .5_RK ) then
          this%P0(i, j) = this%P0(i, j) - 1._RK
        end if
#else
        this%P0(i, j) = this%P0(i, j) - anint( this%P0(i, j) )
#endif
        this%P0old(i, j) = this%P0(i, j)
      end do
    end do

    if( this%Molecule%IsElongated ) then

      ! Predict quaternion parameters and their derivatives
      do j = 1, 4
        do i = 1, np
          this%Q0(i, j) = this%Q0(i, j) + this%Q1(i, j) + this%Q2(i, j) + this%Q3(i, j) + this%Q4(i, j)
          this%Q1(i, j) = this%Q1(i, j) + 2._RK * this%Q2(i, j) + 3._RK * this%Q3(i, j) + 4._RK * this%Q4(i, j)
          this%Q2(i, j) = this%Q2(i, j) + 3._RK * this%Q3(i, j) + 6._RK * this%Q4(i, j)
          this%Q3(i, j) = this%Q3(i, j) + 4._RK * this%Q4(i, j)
        end do
      end do

      ! Predict angular velocities and their derivatives
      do j = 1, nra
        do i = 1, np
          this%W0(i, j) = this%W0(i, j) + this%W1(i, j) + this%W2(i, j) + this%W3(i, j) + this%W4(i, j)
          this%W1(i, j) = this%W1(i, j) + 2._RK * this%W2(i, j) + 3._RK * this%W3(i, j) + 4._RK * this%W4(i, j)
          this%W2(i, j) = this%W2(i, j) + 3._RK * this%W3(i, j) + 6._RK * this%W4(i, j)
          this%W3(i, j) = this%W3(i, j) + 4._RK * this%W4(i, j)
        end do
      end do

    end if

  end subroutine TComponent_PredictGear


!==============================================================!
!  Subroutine TComponent_CorrectGear                           !
!==============================================================!

  subroutine TComponent_CorrectGear( this, dLogVolumeThird )

    implicit none

    ! Declare arguments
    type(TComponent)  :: this
    real(RK)          :: dLogVolumeThird

    ! Declare local variables
    real(RK)          :: BoxLengthInv
    real(RK)          :: MassInv
    real(RK)          :: Moi23, Moi31, Moi12
    real(RK)          :: TMoi1, TMoi2, TMoi3
    real(RK), pointer, contiguous :: pF(:, :), pT(:, :)
    integer           :: np, nra
    integer           :: i, j

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength
    MassInv = 1._RK / this%Molecule%Mass
    np = this%NPart
    nra = this%Molecule%NDFRot

    ! Correct COM positions and their derivatives
#if MPI_VER > 0
    pF => this%FAll(:, :)
#else
    pF => this%F(:, :)
#endif

    do j = 1, 3
      do i = 1, np
        this%Corr0(i, j) = pF(i, j) * TimeStepSquared2 * BoxLengthInv * MassInv

        if( ConstantPressure .and. .not. NVTEquilibration ) this%Corr0(i, j) = this%Corr0(i, j) &
&           - this%P1(i, j) * dLogVolumeThird

        this%Corr1(i, j) = this%Corr0(i, j) - this%P2(i, j)
        this%P0(i, j) = this%P0(i, j) + this%Corr1(i, j) * Gear20
        this%P1(i, j) = this%P1(i, j) + this%Corr1(i, j) * Gear21
        this%P2(i, j) =                    this%Corr0(i, j)
        this%P3(i, j) = this%P3(i, j) + this%Corr1(i, j) * Gear23
        this%P4(i, j) = this%P4(i, j) + this%Corr1(i, j) * Gear24
        this%P5(i, j) = this%P5(i, j) + this%Corr1(i, j) * Gear25

        ! Calculate displacement
        this%Disp(i, j) = this%Disp(i, j) + this%P0(i, j) - this%P0old(i, j)

        ! Check for conservation of particles in primary cell
#if ARCH == 1
        if( this%P0(i, j) < -.5_RK ) then
          this%P0(i, j) = this%P0(i, j) + 1._RK
        elseif( this%P0(i, j) > .5_RK ) then
          this%P0(i, j) = this%P0(i, j) - 1._RK
        end if
#else
        this%P0(i, j) = this%P0(i, j) - anint( this%P0(i, j) )
#endif
        this%P0old(i, j) = this%P0(i, j)
      end do
    end do

    if( this%Molecule%isElongated ) then

      ! Correct quaternion parameters and their derivatives
      do i = 1, np
        this%Corr0(i, 1) = TimeStep2 * ( - this%Q0(i, 2) * this%W0(i, 1) - this%Q0(i, 3) * this%W0(i, 2) &
&                                        - this%Q0(i, 4) * this%W0(i, 3))

        this%Corr0(i, 2) = TimeStep2 * ( + this%Q0(i, 1) * this%W0(i, 1) - this%Q0(i, 4) * this%W0(i, 2) &
&                                        + this%Q0(i, 3) * this%W0(i, 3))

        this%Corr0(i, 3) = TimeStep2 * ( + this%Q0(i, 4) * this%W0(i, 1) + this%Q0(i, 1) * this%W0(i, 2) &
&                                        - this%Q0(i, 2) * this%W0(i, 3))

        this%Corr0(i, 4) = TimeStep2 * ( - this%Q0(i, 3) * this%W0(i, 1) + this%Q0(i, 2) * this%W0(i, 2) &
&                                        + this%Q0(i, 1) * this%W0(i, 3))

      end do
      do j = 1, 4
        do i = 1, np
          this%Corr1(i, j) = this%Corr0(i, j) - this%Q1(i, j)
          this%Q0(i, j) = this%Q0(i, j) + this%Corr1(i, j) * Gear10
          this%Q1(i, j) =                 this%Corr0(i, j)
          this%Q2(i, j) = this%Q2(i, j) + this%Corr1(i, j) * Gear12
          this%Q3(i, j) = this%Q3(i, j) + this%Corr1(i, j) * Gear13
          this%Q4(i, j) = this%Q4(i, j) + this%Corr1(i, j) * Gear14
        end do
      end do

      ! Correct angular velocities and their derivatives
#if MPI_VER > 0
      pT => this%TAll(:, :)
#else
      pT => this%T(:, :)
#endif
      TMoi1 = TimeStep / this%Molecule%MOI(1)
      TMoi2 = TimeStep / this%Molecule%MOI(2)

      if( this%Molecule%is3D ) then
        Moi23 = this%Molecule%MOI(2) - this%Molecule%MOI(3)
        Moi31 = this%Molecule%MOI(3) - this%Molecule%MOI(1)
        Moi12 = this%Molecule%MOI(1) - this%Molecule%MOI(2)
        TMoi3 = TimeStep / this%Molecule%MOI(3)
        do i = 1, np
          this%Corr0(i, 1) = (pT(i, 1) + this%W0(i, 2) * this%W0(i, 3) * Moi23) * TMoi1
          this%Corr0(i, 2) = (pT(i, 2) + this%W0(i, 3) * this%W0(i, 1) * Moi31) * TMoi2
          this%Corr0(i, 3) = (pT(i, 3) + this%W0(i, 1) * this%W0(i, 2) * Moi12) * TMoi3
        end do

      else
        do i = 1, np
          this%Corr0(i, 1) = pT(i, 1) * TMoi1
          this%Corr0(i, 2) = pT(i, 2) * TMoi2
        end do
      end if

      do j = 1, nra
        do i = 1, np
          this%Corr1(i, j) = this%Corr0(i, j) - this%W1(i, j)
          this%W0(i, j) = this%W0(i, j) + this%Corr1(i, j) * Gear10
          this%W1(i, j) = this%Corr0(i, j)
          this%W2(i, j) = this%W2(i, j) + this%Corr1(i, j) * Gear12
          this%W3(i, j) = this%W3(i, j) + this%Corr1(i, j) * Gear13
          this%W4(i, j) = this%W4(i, j) + this%Corr1(i, j) * Gear14
        end do
      end do

    end if

  end subroutine TComponent_CorrectGear



!==============================================================!
!  Subroutine TComponent_PredictLeapFrog                       !
!==============================================================!

  subroutine TComponent_PredictLeapFrog( this, scale )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Declare arguments
    real(RK), intent(in) :: scale

    ! Declare local variables
    real(RK) :: Korr
    integer  :: np, nra
    integer  :: i, j

    Korr = 2._RK - 1._RK / scale
    np = this%NPart

    do j = 1, 3
      do i = 1, np
        this%P1(i, j) = Korr * this%P1(i, j) + this%P2(i, j)
        this%P0(i, j) = this%P0(i, j) + this%P1(i, j)

        ! Calculate displacement
        this%Disp(i, j) = this%Disp(i, j) + this%P0(i, j) - this%P0old(i, j)

        ! Check for conservation of particles in primary cell
#if ARCH == 1
        if( this%P0(i, j) < -.5_RK ) then
          this%P0(i, j) = this%P0(i, j) + 1._RK
        elseif( this%P0(i, j) > .5_RK ) then
          this%P0(i, j) = this%P0(i, j) - 1._RK
        end if
#else
        this%P0(i, j) = this%P0(i, j) - anint( this%P0(i, j) )
#endif
        this%P0old(i, j) = this%P0(i, j)
      end do
    end do

    if( this%Molecule%IsElongated ) then
      nra = this%Molecule%NDFRot
      do i = 1, np
        do j = 1, 4
          this%Q0tmp(i, j) = this%Q0(i, j) + .5_RK * this%Q1(i, j)
        end do

        do j = 1, nra
          this%W0(i, j) = Korr * this%W0(i, j) + .5_RK * this%W1(i, j)
        end do

        this%Q1(i, 1) = TimeStep2 * ( - this%Q0tmp(i, 2) * this%W0(i, 1) - this%Q0tmp(i, 3) * this%W0(i, 2) &
&                                     - this%Q0tmp(i, 4) * this%W0(i, 3))

        this%Q1(i, 2) = TimeStep2 * ( + this%Q0tmp(i, 1) * this%W0(i, 1) - this%Q0tmp(i, 4) * this%W0(i, 2) &
&                                     + this%Q0tmp(i, 3) * this%W0(i, 3))

        this%Q1(i, 3) = TimeStep2 * ( + this%Q0tmp(i, 4) * this%W0(i, 1) + this%Q0tmp(i, 1) * this%W0(i, 2) &
&                                     - this%Q0tmp(i, 2) * this%W0(i, 3))

        this%Q1(i, 4) = TimeStep2 * ( - this%Q0tmp(i, 3) * this%W0(i, 1) + this%Q0tmp(i, 2) * this%W0(i, 2) &
&                                     + this%Q0tmp(i, 1) * this%W0(i, 3))

        do j = 1, 4
          this%Q0(i, j) = this%Q0(i, j) + this%Q1(i, j)
        end do
      end do
    end if

  end subroutine TComponent_PredictLeapFrog



!==============================================================!
!  Subroutine TComponent_CorrectLeapFrog                       !
!==============================================================!

  subroutine TComponent_CorrectLeapFrog( this, dLogVolumeThird )

    implicit none

    ! Declare arguments
    type(TComponent)  :: this
    real(RK)          :: dLogVolumeThird

    ! Declare local variables
    real(RK), pointer, contiguous :: pF(:, :), pT(:, :)
    real(RK)          :: BoxLengthInv, MassInv
    real(RK)          :: Moi23, Moi31, Moi12
    real(RK)          :: TMoi1, TMoi2, TMoi3
    integer           :: np, nra
    integer           :: i, j

    BoxLengthInv = 1._RK / this%BoxLength
    MassInv = 1._RK / this%Molecule%Mass
    np = this%NPart
#if MPI_VER > 0
    pF => this%FAll(:, :)
#else
    pF => this%F(:, :)
#endif

    do j = 1, 3
      do i = 1, np
        this%P2(i, j) = pF(i, j) * TimeStepSquared2 * BoxLengthInv * MassInv &
&                      - (this%P1(i, j) + this%P2(i, j)) * dLogVolumeThird

        this%P1(i, j) = this%P1(i, j) + this%P2(i, j)
      end do
    end do

    if( this%Molecule%IsElongated ) then
      nra = this%Molecule%NDFRot
      TMoi1 = TimeStep / this%Molecule%MOI(1)
      TMoi2 = TimeStep / this%Molecule%MOI(2)
#if MPI_VER > 0
      pT => this%TAll(:, :)
#else
      pT => this%T(:, :)
#endif
      if( this%Molecule%is3D ) then
        TMoi3 = TimeStep / this%Molecule%MOI(3)
        Moi23 = this%Molecule%MOI(2) - this%Molecule%MOI(3)
        Moi31 = this%Molecule%MOI(3) - this%Molecule%MOI(1)
        Moi12 = this%Molecule%MOI(1) - this%Molecule%MOI(2)
        do i = 1, np
          this%W1(i, 1) = (pT(i, 1) + this%W0(i, 2) * this%W0(i, 3) * Moi23) * TMoi1
          this%W1(i, 2) = (pT(i, 2) + this%W0(i, 3) * this%W0(i, 1) * Moi31) * TMoi2
          this%W1(i, 3) = (pT(i, 3) + this%W0(i, 1) * this%W0(i, 2) * Moi12) * TMoi3
        end do

      else
        do i = 1, np
          this%W1(i, 1) = pT(i, 1) * TMoi1
          this%W1(i, 2) = pT(i, 2) * TMoi2
        end do
      end if

      do j = 1, nra
        do i = 1, np
          this%W0(i, j) = this%W0(i, j) + .5_RK * this%W1(i, j)
        end do
      end do
      do i = 1, np
        this%Q1(i, 1) = TimeStep2 * ( - this%Q0(i, 2) * this%W0(i, 1) - this%Q0(i, 3) * this%W0(i, 2) &
&                                     - this%Q0(i, 4) * this%W0(i, 3))

        this%Q1(i, 2) = TimeStep2 * ( + this%Q0(i, 1) * this%W0(i, 1) - this%Q0(i, 4) * this%W0(i, 2) &
&                                     + this%Q0(i, 3) * this%W0(i, 3))

        this%Q1(i, 3) = TimeStep2 * ( + this%Q0(i, 4) * this%W0(i, 1) + this%Q0(i, 1) * this%W0(i, 2) &
&                                     - this%Q0(i, 2) * this%W0(i, 3))

        this%Q1(i, 4) = TimeStep2 * ( - this%Q0(i, 3) * this%W0(i, 1) + this%Q0(i, 2) * this%W0(i, 2) &
&                                     + this%Q0(i, 1) * this%W0(i, 3))

      end do
    end if

  end subroutine TComponent_CorrectLeapFrog


!==============================================================!
!  Subroutine TComponent_PredictVerlet                         !
!==============================================================!

  subroutine TComponent_PredictVerlet( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Issue error
    call Error( 'Subroutine TComponent_PredictVerlet is not implemented' )

  end subroutine TComponent_PredictVerlet



!==============================================================!
!  Subroutine TComponent_CorrectVerlet                         !
!==============================================================!

  subroutine TComponent_CorrectVerlet( this, dLogVolumeThird )

    implicit none

    ! Declare arguments
    type(TComponent) :: this
    real(RK)         :: dLogVolumeThird

    ! Issue error
    call Error( 'Subroutine TComponent_CorrectVerlet is not implemented' )

  end subroutine TComponent_CorrectVerlet



!==============================================================!
!  Subroutine TComponent_PredictVV                             !
!==============================================================!

  subroutine TComponent_PredictVV( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Issue error
    call Error( 'Subroutine TComponent_PredictVV is not implemented' )

  end subroutine TComponent_PredictVV



!==============================================================!
!  Subroutine TComponent_CorrectVV                             !
!==============================================================!

  subroutine TComponent_CorrectVV( this, dLogVolumeThird )

    implicit none

    ! Declare arguments
    type(TComponent) :: this
    real(RK)         :: dLogVolumeThird

    ! Issue error
    call Error( 'Subroutine TComponent_CorrectVV is not implemented' )

  end subroutine TComponent_CorrectVV



!==============================================================!
!  Subroutine TComponent_ZeroNAttempts                         !
!==============================================================!

  subroutine TComponent_ZeroNAttempts( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Zero number of MC attempts and successes
    this%NMoveAttempts    = 0
    this%NMoveSuccesses   = 0
    this%NRotateAttempts = 0
    this%NRotateSuccesses = 0

    if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
      this%NMoveBiasedAttempts = 0
      this%NMoveBiasedSuccesses = 0
      this%NRotateBiasedAttempts = 0
      this%NRotateBiasedSuccesses = 0
!DEBUG
      this%NFLuctUpAttempts(:) = 0
      this%NFluctUpSuccesses(:) = 0
      this%NFluctDownAttempts(:) = 0
      this%NFluctDownSuccesses(:) = 0
!DEBUG
    end if

  end subroutine TComponent_ZeroNAttempts



!==============================================================!
!  Subroutine TComponent_UpdateDisplacements                   !
!==============================================================!

  subroutine TComponent_UpdateDisplacements( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Update translational displacement
    if( this%NMoveSuccesses < this%NMoveAttempts * Acceptance ) then
      this%DispTran = this%DispTran * .95_RK

    else if( this%DispTran < DispTranLimit ) then
      this%DispTran = this%DispTran * 1.05_RK

    end if

    ! Update rotational displacement
    if( this%NRotateSuccesses < this%NRotateAttempts * Acceptance ) then
      this%DispRot = this%DispRot * .95_RK

    else if( this%DispRot < DispRotLimit ) then
      this%DispRot = this%DispRot * 1.05_RK

    end if

  end subroutine TComponent_UpdateDisplacements



!==============================================================!
!  Subroutine TComponent_AddParticle                           !
!==============================================================!

  subroutine TComponent_AddParticle( this, r, q )

    implicit none

    ! Declare arguments
    type(TComponent)               :: this
    real(RK), intent(in)           :: r(3)
    real(RK), intent(in), optional :: q(4)

    ! Test boundaries of particle arrays
    if( this%NPart >= this%NPartMax .and. EnsembleType .eq. EnsembleTypeGE ) then
      tooManyParticles = .true.
      return
    end if

    if( this%NPart > this%NPartMax .and. EnsembleType .ne. EnsembleTypeGE) then
      tooManyParticles = .true.
      return
    end if

    ! Increase NPart
    this%NPart = this%NPart + 1
#if MPI_VER > 0
    this%NPart1 = ProcRange( this%NPart, this%NPart0, this%NPart2 )
#endif

    ! Set coordinates and orientation of new particle
    this%P0(this%NPart, :) = r(:)
    if( this%Molecule%isElongated ) then
      this%Q0(this%NPart, :) = q(:)
    end if

  end subroutine TComponent_AddParticle


!==============================================================!
!  Subroutine TComponent_RemoveParticle                        !
!==============================================================!

  subroutine TComponent_RemoveParticle( this, np )

    implicit none

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: np

    if( np .ne. this%NPart ) then

      ! Copy coordinates and orientation of last particle
      this%P0(np, :) = this%P0(this%NPart, :)
      if( this%Molecule%isElongated ) then
        this%Q0(np, :) = this%Q0(this%NPart, :)
      end if
      call Mol2Atom1( this, np )

    end if

    ! Remove last particle
    this%NPart = this%NPart - 1
#if MPI_VER > 0
    this%NPart1 = ProcRange( this%NPart, this%NPart0, this%NPart2 )
#endif

  end subroutine TComponent_RemoveParticle



!==============================================================!
!  Subroutine TComponent_UpdateChemPot                         !
!==============================================================!

  subroutine TComponent_UpdateChemPot( this, diffpressure )

    implicit none

    ! Declare arguments
    type(TComponent)     :: this
    real(RK), intent(in) :: diffpressure

    ! Update chemical potential
    this%ChemPot = this%ChemPot0 + this%PartialMolarVolume * diffpressure

  end subroutine TComponent_UpdateChemPot



!==============================================================!
!  Subroutine TComponent_SaveState                             !
!==============================================================!

  subroutine TComponent_SaveState( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Save current state
    this%P0Save = this%P0
    if( this%Molecule%IsElongated ) this%Q0Save = this%Q0

  end subroutine TComponent_SaveState



!==============================================================!
!  Subroutine TComponent_RestoreState                          !
!==============================================================!

  subroutine TComponent_RestoreState( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Restore saved state
    this%P0 = this%P0Save
    if( this%Molecule%IsElongated ) this%Q0 = this%Q0Save

    ! Calculate site positions
! #if MPI_VER > 0
!     call Mol2Atom( this, this%NPart0, this%NPart2-this%NPart+1 )
! #else
!     call Mol2Atom( this, 1, this%NPart )
! #endif
    call Mol2Atom( this, this%NPart )

  end subroutine TComponent_RestoreState



!==============================================================!
!  Subroutine TComponent_RestartSave                           !
!==============================================================!

  subroutine TComponent_RestartSave( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Declare local variables
    integer  :: np, i, j
    real(RK) :: pos(3), quat(4)

    ! Assign local variables
    np = this%NPart

    ! Check for root process
    if( .not. RootProc ) return

    ! Save contents to restart file
    write( iounit_restart, '(I10)' ) np

    ! Centers of mass positions
    do i = 1, np
      pos(:) = this%P0(i,:)
      write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) pos(:)
    end do

    if( SimulationType .eq. MolecularDynamics ) then
      ! Centers of mass positions' derivatives
      do i = 1, np
        pos(:) = this%P1(i,:)
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) pos(:)
      end do
      do i = 1, np
        pos(:) = this%P2(i,:)
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) pos(:)
      end do

      if( IntegratorType .eq. IntegratorTypeGear ) then
        do i = 1, np
          pos(:) = this%P3(i,:)
          write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) pos(:)
        end do
        do i = 1, np
          pos(:) = this%P4(i,:)
          write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) pos(:)
        end do
        do i = 1, np
          pos(:) = this%P5(i,:)
          write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) pos(:)
        end do
      end if

      do i = 1, np
        pos(:) = this%Disp(i,:)
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) pos(:)
      end do

      if( ALPHA2UpdateFrequency > 0 ) then
        do i = 1, np
          do j = 0, ALPHA2Length/ALPHA2Shift-1
            write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%ri0_x(i,j),this%ri0_y(i,j),this%ri0_z(i,j)
          end do
        end do
      end if

#if TRANS == 1
      if( EinsteinCoefCalc ) then  !EinsteinCoef ri0_E rest write
            do i = 1, np
              do j = 0, this%NEinstein-1
                write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%ri0_E_x(i,j),this%ri0_E_y(i,j),this%ri0_E_z(i,j)
              end do
            end do
       end if
#endif

    else
      write( iounit_restart, '(ES20.12E3)' ) this%DispTran
      write( iounit_restart, '(2I10)' ) this%NMoveAttempts, this%NMoveSuccesses
      write( iounit_restart, '(2I10)' ) this%NMoveBiasedAttempts, &
&       this%NMoveBiasedSuccesses
    end if

    if( this%Molecule%isElongated ) then
      ! Quaternion parameters
      do i = 1, np
        quat(:) = this%Q0(i,:)
        write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) quat(:)
      end do

      if( SimulationType .eq. MolecularDynamics ) then
        ! Quaternion parameters' derivatives
        do i = 1, np
          quat(:) = this%Q1(i,:)
          write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) quat(:)
        end do

        if( IntegratorType .eq. IntegratorTypeGear ) then
          do i = 1, np
            quat(:) = this%Q2(i,:)
            write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) quat(:)
          end do
          do i = 1, np
            quat(:) = this%Q3(i,:)
            write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) quat(:)
          end do
          do i = 1, np
            quat(:) = this%Q4(i,:)
            write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) quat(:)
          end do
        end if

        ! Angular velocities and their derivatives
        do i = 1, np
          pos(:) = this%W0( i, : )
          write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) pos(:)
        end do
        do i = 1, np
          pos(:) = this%W1( i, : )
          write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) pos(:)
        end do

        if( IntegratorType .eq. IntegratorTypeGear ) then
          do i = 1, np
            pos(:) = this%W2( i, : )
            write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) pos(:)
          end do
          do i = 1, np
            pos(:) = this%W3( i, : )
            write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) pos(:)
          end do
          do i = 1, np
            pos(:) = this%W4( i, : )
            write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) pos(:)
          end do
        end if
      else
        write( iounit_restart, '(ES20.12E3)' ) this%DispRot
        write( iounit_restart, '(2I10)' ) this%NRotateAttempts, this%NRotateSuccesses
        write( iounit_restart, '(2I10)' ) this%NRotateBiasedAttempts, this%NRotateBiasedSuccesses
      end if

    end if

    if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
      write( iounit_restart, '(ES20.12E3)' ) this%WF(:)
      write( iounit_restart, '(I10)' ) this%NState(:)
      write( iounit_restart, '(I10)' ) this%NStateWF(:)
    end if

  end subroutine TComponent_RestartSave



!==============================================================!
!  Subroutine TComponent_RestartRead                           !
!==============================================================!

  subroutine TComponent_RestartRead( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent) :: this

    ! Declare local variables
    integer :: i, j, np

    if( RootProc ) then

      ! Read contents from restart file
      read( iounit_restart, '(I10)' ) np
      if( np > this%NPartMax ) call Error( 'Not enough memory to read particles from restart file' )
      this%NPart = np

      ! Centers of mass positions
      do i = 1, np
        read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%P0( i, : )
      end do

      if( SimulationType .eq. MolecularDynamics ) then
        ! Centers of mass positions' derivatives
        do i = 1, np
          read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%P1( i, : )
        end do

        do i = 1, np
          read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%P2( i, : )
        end do

        if( IntegratorType .eq. IntegratorTypeGear ) then
          do i = 1, np
            read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%P3( i, : )
          end do

          do i = 1, np
            read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%P4( i, : )
          end do

          do i = 1, np
            read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%P5( i, : )
          end do
        end if

        do i = 1, np
          read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%Disp( i, : )
        end do

        if( ALPHA2UpdateFrequency > 0 ) then
          do i = 1, np
            do j = 0, ALPHA2Length/ALPHA2Shift-1
              read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%ri0_x(i,j),this%ri0_y(i,j),this%ri0_z(i,j)
            end do
          end do
        end if
#if TRANS==1
        !EinsteinCoef rest read ri0_E_x
         if( EinsteinCoefCalc ) then
              do i = 1, np
                do j = 0, this%NEinstein-1
                  read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%ri0_E_x(i,j),this%ri0_E_y(i,j),this%ri0_E_z(i,j)
                end do
              end do
         end if
#endif
      else
        read( iounit_restart, '(ES20.12E3)' ) this%DispTran
        read( iounit_restart, '(2I10)' ) this%NMoveAttempts, this%NMoveSuccesses
        read( iounit_restart, '(2I10)' ) this%NMoveBiasedAttempts, this%NMoveBiasedSuccesses
      end if

      if( this%Molecule%isElongated ) then
        ! Quaternion parameters
        do i = 1, np
          read( iounit_restart, '(4(ES20.12E3, :, X))' ) this%Q0( i, : )
        end do

        if( SimulationType .eq. MolecularDynamics ) then
          ! Quaternion parameters' derivatives
          do i = 1, np
            read( iounit_restart, '(4(ES20.12E3, :, X))' ) this%Q1( i, : )
          end do

          if( IntegratorType .eq. IntegratorTypeGear ) then
            do i = 1, np
              read( iounit_restart, '(4(ES20.12E3, :, X))' ) this%Q2( i, : )
            end do

            do i = 1, np
              read( iounit_restart, '(4(ES20.12E3, :, X))' ) this%Q3( i, : )
            end do

            do i = 1, np
              read( iounit_restart, '(4(ES20.12E3, :, X))' ) this%Q4( i, : )
            end do
          end if

          ! Angular velocities and their derivatives
          do i = 1, np
            read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%W0( i, : )
          end do

          do i = 1, np
            read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%W1( i, : )
          end do

          if( IntegratorType .eq. IntegratorTypeGear ) then
            do i = 1, np
              read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%W2( i, : )
            end do

            do i = 1, np
              read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%W3( i, : )
            end do

            do i = 1, np
              read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%W4( i, : )
            end do

          end if

        else
          read( iounit_restart, '(ES20.12E3)' ) this%DispRot
          read( iounit_restart, '(2I10)' ) this%NRotateAttempts, this%NRotateSuccesses
          read( iounit_restart, '(2I10)' ) this%NRotateBiasedAttempts, this%NRotateBiasedSuccesses
        end if
      end if

      if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
        read( iounit_restart, '(ES20.12E3)' ) this%WF(:)
        read( iounit_restart, '(I10)' ) this%NState(:)
        read( iounit_restart, '(I10)' ) this%NStateWF(:)
      end if

    end if

#if MPI_VER > 0
    call MPI_Bcast( this%NPart, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
    call MPI_Bcast( this%P0(:, :), size( this%P0 ), MPI_RK, NRootProc, Communicator, ierror )

    if( this%Molecule%isElongated ) call MPI_Bcast( this%Q0(:, :), size( this%Q0 ), MPI_RK, &
&       NRootProc, Communicator, ierror )

    if( (SimulationType .eq. MonteCarlo) .or. (SimulationType .eq. Gibbs) ) then

      call MPI_Bcast( this%DispTran, 1, MPI_RK, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%NMoveAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%NMoveSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%NMoveBiasedAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%NMoveBiasedSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )

      if( this%Molecule%isElongated ) then
        call MPI_Bcast( this%DispRot, 1, MPI_RK, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NRotateAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NRotateSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NRotateBiasedAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NRotateBiasedSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
      end if
    end if

    if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
      call MPI_Bcast( this%WF, size( this%WF ), MPI_RK, NRootProc, Communicator, ierror )
      call MPI_BCast( this%NState, size( this%NState ), MPI_INTEGER, NRootProc, Communicator, ierror )
      call MPI_BCast( this%NStateWF, size( this%NStateWF ), MPI_INTEGER, NRootProc, Communicator, ierror )
    end if
#endif

    ! Update old positions
    this%P0old = this%P0

  end subroutine TComponent_RestartRead


!TRANSPORT_start
!==============================================================!
!  Subroutine TComponent_ForceTransport                        !
!==============================================================!

subroutine TComponent_ForceTransport( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent)  :: this
#if TRANS==1
    integer           :: i, j, k, nra
    real(RK), pointer, contiguous :: pFTC1(:,:), pFTC2(:,:), pFTC3(:,:)
    real(RK), pointer, contiguous :: pFRC1(:,:), pFRC2(:,:), pFRC3(:,:)
    real(RK)          :: BoxLength_dt
    real(RK)          :: BoxLength_dt2

    !declare local variables
    BoxLength_dt = this%BoxLength/TimeStep
    BoxLength_dt2 = BoxLength_dt*BoxLength_dt


    this%FTC(:,:) = 0._RK
    this%FRC(:,:) = 0._RK
    this%KinETran(:,:) = 0._RK
    this%KinEPart(:) = 0._RK

    if (RootProc) then

   !   if (this%Conductivity) then
#if MPI_VER > 0
        pFTC1 => this%FTC1All(:,:)
        pFTC2 => this%FTC2All(:,:)
        pFTC3 => this%FTC3All(:,:)

        pFRC1 => this%FRC1All(:,:)
        pFRC2 => this%FRC2All(:,:)
        pFRC3 => this%FRC3All(:,:)
#else
        pFTC1 => this%FTC1(:,:)
        pFTC2 => this%FTC2(:,:)
        pFTC3 => this%FTC3(:,:)

        pFRC1 => this%FRC1(:,:)
        pFRC2 => this%FRC2(:,:)
        pFRC3 => this%FRC3(:,:)
#endif
        do k = 1, 3
          do i= 1, this%Npart
            this%FTC(i,1)= this%FTC(i,1)+ pFTC1(i,k)*this%P1(i,k)
            this%FTC(i,2)= this%FTC(i,2)+ pFTC2(i,k)*this%P1(i,k)
            this%FTC(i,3)= this%FTC(i,3)+ pFTC3(i,k)*this%P1(i,k)
          end do
        end do

        nra = this%Molecule%NDFRot
        do k= 1, nra
          do i= 1, this%Npart
            this%FRC(i,1)= this%FRC(i,1)+ pFRC1(i,k)*this%W0(i,k)
            this%FRC(i,2)= this%FRC(i,2)+ pFRC2(i,k)*this%W0(i,k)
            this%FRC(i,3)= this%FRC(i,3)+ pFRC3(i,k)*this%W0(i,k)
          end do
        end do

        this%FTC(:,:) = this%FTC(:,:)*BoxLength_dt


        ! Calculate kinetic energy / molecule
        do j = 1, 3
          this%KinETran(:,j) = this%P1(:,j)*this%P1(:,j)
        end do

        this%KinETran(:,:) = this%KinETran(:,:)* this%Molecule%Mass*BoxLength_dt2

        do j = 1,3
          this%KinEPart(:) = this%KinEPart(:) + this%KinETran(:,j)
        end do

        this%KinETranTotal(1) = sum(this%KinETran(:,1))
        this%KinETranTotal(2) = sum(this%KinETran(:,2))
        this%KinETranTotal(3) = sum(this%KinETran(:,3))

    end if ! RootProc
#endif

  end subroutine TComponent_ForceTransport
!TRANSPORT_END


#if CONSTR > 0
!==============================================================!
!  Subroutine TComponent_CorrectGear                           !
!==============================================================!

  subroutine TComponent_CorrectGear_Constraint(this,aa,dLogVolumeThird,Forc,drx,dry,drz )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent)         :: this
    real(RK),intent(in)      :: dLogVolumeThird
    integer,intent (in)      :: aa
    real(RK), intent(in out) :: Forc
    real(RK),intent(in)      :: drx,dry,drz

    ! Declare local variables
    real(RK)          :: BoxLength
    real(RK)          :: Mass
    real(RK)          :: np
    real(RK)          :: ff
    real(RK)          :: Corr0,Corr0ff,Corr1
    real(RK)          :: dr(3)
    integer           :: i, j

    ! Assign local variables
    BoxLength = this%BoxLength
    Mass = this%Molecule%Mass
    np = 2
    dr(1) = drx
    dr(2) = dry
    dr(3) = drz

    ! Correct COM positions and their derivatives
    do j = 1, 3,1

      Corr1 = + dr(j) / Gear20
      Corr0 = Corr1 + this%P2(aa,j)

      Corr0ff = Corr0
      if (ConstantPressure .and. .not. NVTEquilibration) Corr0ff = Corr0ff + this%P1(aa,j)*dLogVolumeThird

        ff = Corr0ff * BoxLength* Mass / TimeStepSquared2
        Forc = Forc + ff
        this%P0(aa, j) = this%P0(aa, j) + Corr1 * Gear20

        ! Check for conservation of particles in primary cell

#if ARCH == 1
        if( this%P0(aa, j) < -.5_RK ) then
          this%P0(aa, j) = this%P0(aa, j) + 1._RK
        elseif( this%P0(i, j) > .5_RK ) then
          this%P0(aa, j) = this%P0(aa, j) - 1._RK
        end if
#else
        this%P0(aa, j) = this%P0(aa, j) - anint( this%P0(aa, j) )
#endif
        this%P0old(aa, j) = this%P0(aa, j)
    end do

  end subroutine TComponent_CorrectGear_Constraint

#endif



end module ms2_component
