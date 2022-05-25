!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 4.0                !
!  (c) 2020 by TU Kaiserslautern / TU Berlin                   !
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

#include "mathMacros.F90"

module ms2_component

#if MPI_VER > 0 && defined(MPI_USE_MODULE)
  use mpi_f08
#endif

  use ms2_accumulator
  use ms2_global
  use ms2_molecule
  use ms2_site
  use math_types



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

    ! Positions and orientations for units of test particles
    real(RK), pointer :: P0Test(:, :, :), Q0Test(:, :, :)

    ! Centers of mass positions for molecules
    real(RK), pointer, contiguous :: Pm0(:, :)
    real(RK), pointer, contiguous :: P0Save( :, :, :)
    ! Centers of mass positions and their derivatives for Units
    real(RK), pointer, contiguous :: P0(:, :, :)
    real(RK), pointer, contiguous :: P1(:, :, :)
    real(RK), pointer, contiguous :: P2(:, :, :)
    real(RK), pointer, contiguous :: P3(:, :, :)
    real(RK), pointer, contiguous :: P4(:, :, :)
    real(RK), pointer, contiguous :: P5(:, :, :)

    ! Quaternion parameters and their derivatives for Units
    real(RK), pointer, contiguous :: Q0(:, :, :)
    real(RK), pointer, contiguous :: Q0Save(:, :, :)
    real(RK), pointer, contiguous :: Q0tmp(:, :, :)
    real(RK), pointer, contiguous :: Q1(:, :, :)
    real(RK), pointer, contiguous :: Q2(:, :, :)
    real(RK), pointer, contiguous :: Q3(:, :, :)
    real(RK), pointer, contiguous :: Q4(:, :, :)

    ! Angular velocities and their derivatives for units
    real(RK), pointer, contiguous :: W0(:, :, :)
    real(RK), pointer, contiguous :: W1(:, :, :)
    real(RK), pointer, contiguous :: W2(:, :, :)
    real(RK), pointer, contiguous :: W3(:, :, :)
    real(RK), pointer, contiguous :: W4(:, :, :)

    ! Intramolecular energy of test particles
    real(RK), pointer, contiguous :: EPotTestIntra(:)

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
    ! Total forces acting on units
    real(RK), pointer, contiguous :: F(:, :, :)
#if MPI_VER > 0
    real(RK), pointer, contiguous :: FAll(:, :, :)
#endif

    ! Total torques acting on units
    real(RK), pointer, contiguous :: T(:, :, :)
#if MPI_VER > 0
    real(RK), pointer, contiguous :: TAll(:, :, :)
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

    ! Total dipole moment of units of a molecule for reaction field
    real(RK), pointer, contiguous :: MueX(:, :), MueY(:, :), MueZ(:, :)

    ! Torques from reaction field, space fixed
    real(RK), pointer, contiguous :: tRFX(:, :), tRFY(:, :), tRFZ(:, :)

    ! Total dipole moment of test particles for reaction field
    real(RK), pointer, contiguous :: MueXTest(:, :), MueYTest(:, :), MueZTest(:, :)

    ! Gear corrector local arrays
    real(RK), pointer, contiguous :: Corr0(:, :, :)
    real(RK), pointer, contiguous :: Corr1(:, :, :)

    ! Length of simulation box
    real(RK), pointer :: BoxLength

    ! Mole fraction of this component
    real(RK) :: Fraction

    ! Maximum number of particles in component
    integer, pointer :: NPartMax

    ! Maximum number of units in component
    integer, pointer :: NUnitMax

    ! Number of particles in component
    integer, pointer :: NPart

    ! Number of particles in process
    ! Starting position, Number of Particles, Endposition
    integer, pointer :: NPart0, NPart1, NPart2

    ! Number of test particles
    integer, pointer :: NTest
    integer          :: NTestAll
    integer, pointer :: NTest0, NTest1, NTest2

    ! Number of degrees of freedom
    integer :: NDFTran, NDFRot, NDF

    ! Maximum allowed MC displacements
    real(RK), pointer :: DispTran, DispRot
    real(RK), pointer :: DispMolTran, DispMolRot

    ! Number of MC attempts and successes
    integer :: NMoveAttempts, NMoveSuccesses
    integer :: NRotateAttempts, NRotateSuccesses
    integer :: NMoveBiasedAttempts, NMoveBiasedSuccesses
    integer :: NRotateBiasedAttempts, NRotateBiasedSuccesses

    ! Number of MC attempts and successes for IDF
    integer :: NMoveMolAttempts, NMoveMolSuccesses
    integer :: NRotateMolAttempts, NRotateMolSuccesses
    integer :: NMoveBiasedMolAttempts, NMoveBiasedMolSuccesses
    integer :: NRotateBiasedMolAttempts, NRotateBiasedMolSuccesses

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

    ! IDF
    integer, pointer, contiguous :: UnitLJ(:),UnitC(:),UnitDP(:),UnitQP(:)

! Ewald
    real(RK) :: EPotTestSelf

    ! Fluctuating components and weighting factors
    integer           :: NFluctState, NFluctMax
    integer, pointer, contiguous  :: NState(:), NStateWF(:), NFluctComp(:)
    real(RK), pointer, contiguous :: WF(:)
    real(RK)          :: ProbW0, ProbW1, ProbW0V, ProbW1Rho
!DEBUG
    integer, pointer  :: NFluctUpAttempts(:), NFluctUpSuccesses(:)
    integer, pointer  :: NFluctDownAttempts(:), NFluctDownSuccesses(:)
!     integer, pointer  :: NStateBF(:)
!     real(RK), pointer :: BFSumState(:)
!DEBUG

    ! Variables for Thermodynamic Integration
    integer           :: NBins
    integer           :: changeLaFreq, changeLaPart, forfeitLaSampl
    integer, pointer, contiguous  :: BinsVisit(:)
    real(RK)          :: Lambda, LambdaExponent, LaMin, LaMax, deltaLa, LaStepMax, ExpMinusBetaEnLaMin
    real(RK)          :: currentBinsEn
    real(RK), pointer, contiguous :: BinsEn(:), BinsdEndLa(:), BinsIntdEndLa(:), BinsdEndLaV(:), BinsdEndLaH(:), BinsIntVW(:), BinsIntHW(:)

    ! Mole fraction in corresponding liquid simulation (for GE ensemble only)
    real(RK) :: LiqFraction

    ! Long-range corrections
    real(RK) :: EPotTestCorrMIE
    real(RK) :: EPotTestCorrTT68
    real(RK) :: EPotTestCorrRF
    
    ! Internal degrees of freedom
    integer,pointer, contiguous :: BondCount(:)
    integer,pointer, contiguous :: BoPartner(:,:)
    integer,pointer, contiguous :: AngleCount(:)
    integer,pointer, contiguous :: AnglePartner(:,:)
    integer,pointer, contiguous :: DihedralCount(:)
    integer,pointer, contiguous :: DihedralPartner(:,:)

    ! Accumulated sums, averages and errors
    type(TAccumulator) :: SumInvChemPotRho
    type(TAccumulator) :: SumInvChemPot
    type(TAccumulator) :: SumChemPotV
    type(TAccumulator) :: SumChemPotVV
    type(TAccumulator) :: SumChemPotGE
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

  interface RemoveNetMomentum
    module procedure TComponent_RemoveNetMomentum
  end interface

  interface CalculateEKin
    module procedure TComponent_CalculateEKin
  end interface

  interface ResizeMol
    module procedure TComponent_ResizeMol
  end interface

  interface RotateMol
    module procedure TComponent_RotateMol
  end interface
  
  interface RotateTest
    module procedure TComponent_RotateTest
  end interface

  interface InitUnit
    module procedure TComponent_InitUnit
  end interface

  interface Unit2Atom
    module procedure TComponent_Unit2Atom
  end interface
  
  interface Unit2Atom1
    module procedure TComponent_Unit2Atom1Mol
    module procedure TComponent_Unit2Atom1
    module procedure TComponent_Unit2Atom1Root
  end interface

  interface Unit2AtomShake
    module procedure TComponent_Unit2AtomShake
  end interface

  interface Unit2AtomTest
    module procedure TComponent_Unit2AtomTest
  end interface

  interface Atom2Unit
    module procedure TComponent_Atom2Unit
  end interface

  interface Atom2Unit_Trans
    module procedure TComponent_Atom2Unit_Trans
  end interface

  interface Unit2Mol
    module procedure TComponent_Unit2Mol
    module procedure TComponent_Unit2Mol1
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
  
  interface ReverseLeapFrog
    module procedure TComponent_ReverseLeapFrog
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

  interface Constraints
    module procedure TComponent_Constraints
  end interface

  interface AddParticle
    module procedure TComponent_AddParticle
  end interface
  
  interface DuplicateParticle
    module procedure TComponent_DuplicateParticle
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
#if MPI_VER > 0 && !defined(MPI_USE_MODULE)
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
    allocate( this%NTest0, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NTest1, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NTest2, STAT = stat )
    call AllocationError( stat, 'number of particles' )

    ! Read file name for potential model
    call FileReadParameter( this%PotModFileName, paramsFile%iounit , IdPotModFileName, .false. )

    ! Read mole fraction of this component
    write( IOBuffer, '(72("-"))')
    call LogWrite
    write( IOBuffer, '(T13, "Reading component", I3," for ensemble")') comp
    call LogWrite
    call FileReadParameter( this%Fraction, paramsFile%iounit , IdFraction, .false. )
    write( IOBuffer, '("Mole fraction of component ", A, ": ", F9.6)' ) trim( this%PotModFileName ), this%Fraction
    call LogWrite

#if TRANS==1
 ! Read partial molar enthalpy from the paremeters file
    call FileReadParameter( this%PartialMolarEnthalpy, paramsFile%iounit , IdPartialMolarEnthalpy, .false., 0._RK )

    if (this%PartialMolarEnthalpy .ne. 0._RK) then
      write( IOBuffer,'("Reduced PartMolEnt of component ", A, ": ", F9.6 )' ) &
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
      call FileReadParameter( this%LiqFraction, paramsFile%iounit , IdLiqFraction, .false. )

      ! Read chemical potential and partial molar volume and their
      ! uncertainties for Grand Equilibrium
      call FileReadParameter( this%ChemPot0, paramsFile%iounit , IdChemPot, .false. )
      call FileReadParameter( this%VarChemPot, paramsFile%iounit , IdVarChemPot, .false. )
      call FileReadParameter( this%PartialMolarVolume, paramsFile%iounit , IdPartialMolarVolume, .false. )
      call FileReadParameter( this%VarPartialMolarVolume, paramsFile%iounit , IdVarPartialMolarVolume, .false. )
      write( IOBuffer,'("Reduced ChemPot0 of component ", A, ": ", F9.6, " (", F9.6, ")")' ) &
&       trim( this%PotModFileName ), this%ChemPot0, this%VarChemPot
      call LogWrite
      write( IOBuffer,'("Reduced PartMolVol of component ", A, ": ", F9.6, " (", F9.6, ")")' ) &
&       trim( this%PotModFilename ), this%PartialMolarVolume, this%VarPartialMolarVolume
      call LogWrite

    else if( EnsembleType .eq. EnsembleTypeHA ) then
      if( comp == 1 ) then
        ! Read chemical potential of phase changing component (first one)
        call FileReadParameter( this%ChemPot, paramsFile%iounit , IdChemPot, .false. )
        call FileReadParameter( this%VarChemPot, paramsFile%iounit , IdVarChemPot, .false. )
        write( IOBuffer, '("Reduced ChemPot of component ", A, ": ", F9.6, " (", F9.6, ")")' ) &
&         trim( this%PotModFileName ), this%ChemPot0, this%VarChemPot
        call LogWrite
      end if

    else if( EnsembleType .eq. EnsembleTypeMUVT ) then
      ! Read chemical potential
      call FileReadParameter( this%ChemPot0, paramsFile%iounit , IdChemPot, .false. )
      write( IOBuffer,'("Reduced ChemPot0 of component ", A, ": ", F9.6, " (", F9.6, ")")' ) &
      &       trim( this%PotModFileName ), this%ChemPot0, this%VarChemPot
      call LogWrite

    else
      ! Read method for calculation of chemical potential
      call FileReadParameter( str, paramsFile%iounit, IdChemPotMethod, .false., "NONE" )
      select case( str )
      case( 'NONE', 'None', 'none' )
        this%ChemPotMethod = ChemPotMethodNone
        str = 'no calculation'
            if (EnsembleType .eq. EnsembleTypeNPTSVC) then
        this%ChemPotMethod = ChemPotMethodWidom
        end if
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
        call FileReadParameter( this%GradInsInit, paramsFile%iounit , IdGradInsInit, .false., 0 )
        write( IOBuffer, '("Grad. Ins. initialization Steps: ", T40, I7)' ) this%GradInsInit
        call LogWrite
      end if

      ! Read number of test particles
      if( this%ChemPotMethod .eq. ChemPotMethodWidom ) then
        call FileReadParameter( this%NTest, paramsFile%iounit, IdNTest, .false. )
        if( this%NTest <= 0 ) call Error( 'Number of test particles need to be > 0' )
        write( IOBuffer, '(T10, "-> Number of test particles:", I11 )' ) this%NTest


#if MPI_VER > 0
        if (( (SimulationType .eq. MolecularDynamics) .or. (mpiMCCommonGroups > 0) ) .and. .not. UseIntDegFreed) then
           this%NTest = ((this%NTest -1)/NProcs +1) !for mpiMCCommonGroups > 0: inside a group chemical potential calculation is parallelized by this%NTest
        endif
#endif
      end if

      ! Read weighting factors method
      this%WFMethod = WFMethodNone
      if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
        call FileReadParameter( str, paramsFile%iounit, IdWeightFactors, .false. )
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
        if (UseIntDegFreed) then
            call FileReadParameter( this%LaMin, paramsFile%iounit , IdLambdaMin, .false., 0.1_RK )
        else
            call FileReadParameter( this%LaMin, paramsFile%iounit , IdLambdaMin, .false., 0.2_RK )
        end if
        write( IOBuffer, '("Thermo. Int. LambdaMin: ", T40, F8.5)' ) this%LaMin
        call LogWrite
        call FileReadParameter( this%LaMax, paramsFile%iounit , IdLambdaMax, .false., 1.0_RK )
        write( IOBuffer, '("Thermo. Int. LambdaMax: ", T40, F8.5)' ) this%LaMax
        call LogWrite
        call FileReadParameter( this%NBins, paramsFile%iounit , IdNBins, .false., 100 )
        write( IOBuffer, '("Thermo. Int. NBins: ", T40, I8)' ) this%NBins
        call LogWrite
        if (SimulationType .eq. MolecularDynamics) then
            if (.not. UseIntDegFreed) then
                call FileReadParameter( this%LaStepMax, paramsFile%iounit , IdLambdaStepMax, .false., 0.01_RK)
            else
                write( IOBuffer, '("In MD simulations LambdaStepMax is determined by NBins, LambdaMax and LambdaMin.")' )
                call LogWrite
                this%LaStepMax = -(this%LaMax-this%LaMin) / (this%NBins-1)
                this%deltaLa = -this%LaStepMax
            end if
        else
          call FileReadParameter( this%LaStepMax, paramsFile%iounit , IdLambdaStepMax, .false., 0.1_RK)
          this%deltaLa = (this%LaMax-this%LaMin) / this%NBins
        end if
        write( IOBuffer, '("Thermo. Int. LambdaStepMax: ", T40, F8.5)' ) this%LaStepMax
        call LogWrite
        call FileReadParameter( this%LambdaExponent, paramsFile%iounit , IdLambdaExponent, .false., 4.0_RK)
        write( IOBuffer, '("Thermo. Int. LambdaExponent: ", T40, F8.5)' ) this%LambdaExponent
        call LogWrite
        if (UseIntDegFreed) then
            call FileReadParameter( this%NTest, paramsFile%iounit, IdNTest, .false., 25 ) ! Michael Sch.: changed from 250
        else
            call FileReadParameter( this%NTest, iounit_params, IdNTest, .false., 100 )
        end if
        write( IOBuffer, '(T10, "-> Number of test particles:", I11 )' ) this%NTest
        call LogWrite

        if (UseIntDegFreed) then
            ! Michael Sch.: new for enhanced possibilites for E(la) sampling
            if (SimulationType .eq. MolecularDynamics) then
              call FileReadParameter( this%changeLaFreq, iounit_params , IdchangeLaFreq, .false., 5000)
            else
              call FileReadParameter( this%changeLaFreq, iounit_params , IdchangeLaFreq, .false., 100)
            end if
            write( IOBuffer, '("Frequency for lambda changes:", I11 )' ) this%changeLaFreq
            call LogWrite
            this%changeLaPart = 1.8*this%changeLaFreq*(this%LaMax-this%LaMin)/this%LaStepMax
            if (SimulationType .eq. MonteCarlo) this%changeLaPart=10*this%changeLaPart
            ! read parameter in case the user wants no or more different particles to be sampled...not advised too much^^
            call FileReadParameter( this%changeLaPart, iounit_params , IdchangeLaPart, .false., this%changeLaPart)
            write( IOBuffer, '("Frequency for changing the fluctuating particle:", I11 )' ) this%changeLaPart
            call LogWrite
            if (SimulationType .eq. MolecularDynamics) then
              call FileReadParameter( this%forfeitLaSampl, iounit_params , IdforfeitLaSampl, .false., 2000)
            else
              call FileReadParameter( this%forfeitLaSampl, iounit_params , IdforfeitLaSampl, .false., 40)
            end if
            write( IOBuffer, '("Forfeit steps of each lambda value regarding the energy sampling:", I11 )' ) this%forfeitLaSampl
            call LogWrite
        end if

#if MPI_VER>0
        if ( (SimulationType .eq. MolecularDynamics) .or. (mpiMCCommonGroups > 0) ) then
          this%NTest = ((this%NTest-1)/NProcs +1) !for mpiMCCommonGroups > 0: inside a group chemical potential calculation is parallelized by this%NTest
        endif
#endif
        if (this%LaMin**this%LambdaExponent .lt. 1E-30_RK) then
          this%LaMin = 1E-30_RK**(1._RK/this%LambdaExponent)
          if (.not. UseIntDegFreed) then
              write( IOBuffer, '("LambdaMin too low for simulation and was changed!")')
          else
              write( IOBuffer, '("LambdaMin too low for simulation! Value was changed to: ", F8.5)' ) this%LaMin
          end if
          call LogWrite
        endif
        if (.not. UseIntDegFreed) this%deltaLa=(this%LaMax-this%LaMin)/this%NBins
      end if

    end if

#if OSMOP > 0
    call FileReadParameter( str, paramsFile%iounit, IdPermeability, .false.)
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

    ! Set Unit Borders
    this%UnitLJ => this%Molecule%UnitLJ
    this%UnitC  => this%Molecule%UnitC
    this%UnitDP => this%Molecule%UnitDP
    this%UnitQP => this%Molecule%UnitQP

    ! Allocate maximum allowed MC displacements
    if( (SimulationType .eq. MonteCarlo) .or. (SimulationType .eq. Gibbs) .or. MCOverlapReduction ) then
      allocate( this%DispTran, STAT = stat )
      call AllocationError( stat, 'maximum MC displacement' )
      allocate( this%DispRot, STAT = stat )
      call AllocationError( stat, 'maximum MC displacement' )
      allocate( this%DispMolTran, STAT = stat )
      call AllocationError( stat, 'maximum MC molecule displacement' )
      allocate( this%DispMolRot, STAT = stat )
      call AllocationError( stat, 'maximum MC molecule displacement' )
    end if

    ! Allocate and read weighting factors
    if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
      allocate( this%WF( 0:this%NFluctMax ), STAT = stat )
      call AllocationError( stat, 'fluctuating particle states', &
&       this%NFluctMax + 1 )
      if( this%WFMethod .eq. WFMethodGuess .or. this%WFMethod .eq. WFMethodOptSet ) then
        if( RootProc ) read( paramsFile%iounit, * ) this%WF
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
    write( IOBuffer, '(72("-"))')
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
    allocate( this%NTest0, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NTest1, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NTest2, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    this%NTest = 0

    ! Set file name for potential model
    this%PotModFileName = PotModFileName

    ! Set mole fraction of this component
    this%Fraction = 1._RK

    ! Create potential model
    call Construct( this%Molecule, this%PotModFileName, -1 )

    ! Set Unit Borders
    this%UnitLJ => this%Molecule%UnitLJ
    this%UnitC  => this%Molecule%UnitC
    this%UnitDP => this%Molecule%UnitDP
    this%UnitQP => this%Molecule%UnitQP

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
    allocate( this%NTest0, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NTest1, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NTest2, STAT = stat )
    call AllocationError( stat, 'number of particles' )

    ! Copy file name for potential model
    this%PotModFileName = comp0%PotModFileName

    ! Set number of particles and mole fraction of this component
    this%NPart = 0
    this%NTest = 0
    this%Fraction = 0._RK
    this%NBins = 0

    if (UseIntDegFreed) then
        this%Lambda = 1.0_RK !- Zero ! test Minh
    else
        this%Lambda = 1.0_RK - Zero ! test Minh
    end if

    ! Set fluctuating state (for GradIns)
    this%FluctState = 0
    this%ChemPotMethod = ChemPotMethodNone

    ! Set chemical potential flag
    this%CalcChemPot = .false.

    ! Create potential model
    call Construct( this%Molecule, this%PotModFileName, -1 )

    ! Set Unit Borders
    this%UnitLJ => this%Molecule%UnitLJ
    this%UnitC  => this%Molecule%UnitC
    this%UnitDP => this%Molecule%UnitDP
    this%UnitQP => this%Molecule%UnitQP

    ! Set maximum allowed MC displacements
    this%DispTran => comp0%DispTran
    this%DispRot => comp0%DispRot
    this%DispMolTran => comp0%DispMolTran
    this%DispMolRot => comp0%DispMolRot

    ! Set Degrees of Freedom
    this%Molecule%Unit(1:this%Molecule%NUnit)%NDF = comp0%Molecule%Unit(1:comp0%Molecule%NUnit)%NDF

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
    allocate( this%NTest0, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NTest1, STAT = stat )
    call AllocationError( stat, 'number of particles' )
    allocate( this%NTest2, STAT = stat )
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

    ! Set Unit Borders
    this%UnitLJ => this%Molecule%UnitLJ
    this%UnitC  => this%Molecule%UnitC
    this%UnitDP => this%Molecule%UnitDP
    this%UnitQP => this%Molecule%UnitQP

    ! Set maximum allowed MC displacements
    this%DispTran => comp0%DispTran
    this%DispRot => comp0%DispRot
    this%DispMolTran => comp0%DispMolTran
    this%DispMolRot => comp0%DispMolRot

    ! Set Degrees of Freedom
    this%Molecule%Unit(1:this%Molecule%NUnit)%NDF = comp0%Molecule%Unit(1:comp0%Molecule%NUnit)%NDF

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
    if( associated( this%NTest0 ) ) then
      deallocate( this%NTest0 )
    end if
    if( associated( this%NTest1 ) ) then
      deallocate( this%NTest1 )
    end if
    if( associated( this%NTest2 ) ) then
      deallocate( this%NTest2 )
    end if
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
      if( associated( this%DispMolTran ) ) then
        deallocate( this%DispMolTran )
      end if
      if( associated( this%DispMolRot ) ) then
        deallocate( this%DispMolRot )
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
    if( associated( this%NTest0 ) ) then
      deallocate( this%NTest0 )
    end if
    if( associated( this%NTest1 ) ) then
      deallocate( this%NTest1 )
    end if
    if( associated( this%NTest2 ) ) then
      deallocate( this%NTest2 )
    end if
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
    if (EnsembleType .eq. EnsembleTypeGE) then
      call Construct( this%SumChemPotGE, .false. )
    endif

    if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeMUVT .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
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
      call Destruct( this%SumChemPotThermoIntWidom )
      call Destruct( this%SumChemPotThermoIntWidomV )
      call Destruct( this%SumHW_counter )
      call Destruct( this%SumHW_denom )
      call Destruct( this%SumVW )
      call Destruct( this%SumHM )
    end select
    if (EnsembleType .eq. EnsembleTypeGE) then
      call Destruct( this%SumChemPotGE )
    endif

    if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeMUVT .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
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
    integer :: nu, nup, nlj, nch, ndi, nqu, j
    integer :: i
    integer :: stat
    logical :: Site1, Site2, Site3, Site4
    integer :: SiteId1, SiteId2, SiteId3, SiteId4

    ! Set maximum number of particles and number of test particles
    np = this%NPartMax
    nu = this%Molecule%NUnit
    nup = nu*np

    ntest = this%NTest

    ! Nullify pointers
    nullify( this%P0 )
    nullify( this%EPotTestIntra )
    nullify( this%Pm0 )
    nullify( this%P0Save )
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

    !EinsteinCoef ri0_E nullify
    if (( TransMethod .eq. Einstein) .or. (TransMethod .eq. GKEinstein)) then
      nullify( this%ri0_E_x ) 
      nullify( this%ri0_E_y )
      nullify( this%ri0_E_z )
    end if  

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
    allocate( this%Q0( np, 4, nu ), STAT = stat )
    call AllocationError( stat, 'units*4*particles', nup )

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
    if (ntest > 0) then
      allocate( this%EPotTestIntra( ntest ), STAT = stat )
      call AllocationError( stat, 'testparticles', ntest )
    end if
    allocate( this%Pm0( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%P0Save( np, 3, nu ), STAT = stat )
    call AllocationError( stat, 'units*particles', nup )
    ! Centers of mass positions for Units
    allocate( this%P0( np, 3, nu ), STAT = stat )
    call AllocationError( stat, 'units*particles', nup )

    if( SimulationType .eq. MolecularDynamics ) then

      ! Centers of mass positions' derivatives
      allocate( this%P1( np, 3, nu ), STAT = stat )
      call AllocationError( stat, 'units*particles', nup )
      allocate( this%P2( np, 3, nu ), STAT = stat )
      call AllocationError( stat, 'units*particles', nup )

      if( IntegratorType .eq. IntegratorTypeGear ) then
        allocate( this%P3( np, 3, nu ), STAT = stat )
        call AllocationError( stat, 'units*particles', nup )
        allocate( this%P4( np, 3, nu ), STAT = stat )
        call AllocationError( stat, 'units*particles', nup )
        allocate( this%P5( np, 3, nu ), STAT = stat )
        call AllocationError( stat, 'units*particles', nup )
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
      if( (TransMethod .eq. Einstein) .or. (TransMethod .eq. GKEinstein))then
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
      allocate( this%F( np, 3, nu ), STAT = stat )
      call AllocationError( stat, 'units*particles', nup )
#if MPI_VER > 0
      allocate( this%FAll( np, 3, nu ), STAT = stat )
      call AllocationError( stat, 'units*particles', nup )
      allocate( this%NAdd( np ), STAT = stat )
      call AllocationError( stat, 'units*particles', nup )
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
      allocate( this%Q0( np, 4, nu ), STAT = stat )
      call AllocationError( stat, 'units*particles', nup )
#endif
      allocate( this%Q0Save( np, 4, nu ), STAT = stat )
      call AllocationError( stat, 'units*particles', nup )

      if( SimulationType .eq. MolecularDynamics ) then

        if( IntegratorType .eq. IntegratorTypeLeapFrog ) then
          allocate( this%Q0tmp( np, 4, nu ), STAT = stat )
          call AllocationError( stat, 'units*particles', nup )
        end if

        ! Quaternion parameters' derivatives
        allocate( this%Q1( np, 4, nu ), STAT = stat )
        call AllocationError( stat, 'units*particles', nup )

        if( IntegratorType .eq. IntegratorTypeGear ) then
          allocate( this%Q2( np, 4, nu ), STAT = stat )
          call AllocationError( stat, 'units*particles', nup )
          allocate( this%Q3( np, 4, nu ), STAT = stat )
          call AllocationError( stat, 'units*particles', nup )
          allocate( this%Q4( np, 4, nu ), STAT = stat )
          call AllocationError( stat, 'units*particles', nup )
        end if

        ! Angular velocities and their derivatives
        allocate( this%W0( np, 3, nu ), STAT = stat )
        call AllocationError( stat, 'units*particles', nup )
        allocate( this%W1( np, 3, nu ), STAT = stat )
        call AllocationError( stat, 'units*particles', nup )

        if( IntegratorType .eq. IntegratorTypeGear ) then
          allocate( this%W2( np, 3, nu ), STAT = stat )
          call AllocationError( stat, 'units*particles', nup )
          allocate( this%W3( np, 3, nu ), STAT = stat )
          call AllocationError( stat, 'units*particles', nup )
          allocate( this%W4( np, 3, nu ), STAT = stat )
          call AllocationError( stat, 'units*particles', nup )
        end if

        ! Total torques
        allocate( this%T( np, 3, nu ), STAT = stat )
        call AllocationError( stat, 'units*particles', nup )
#if MPI_VER > 0
        allocate( this%TAll( np, 3, nu ), STAT = stat )
        call AllocationError( stat, 'units*particles', nup )
#endif
        ! Torques from reaction field
        allocate( this%tRFX( np,nu ), STAT = stat )
        call AllocationError( stat, 'units*particles', nup )
        allocate( this%tRFY( np,nu ), STAT = stat )
        call AllocationError( stat, 'units*particles', nup )
        allocate( this%tRFZ( np,nu ), STAT = stat )
        call AllocationError( stat, 'units*particles', nup )
      end if

      ! Total dipole moment of units for reaction field
      if( CutoffMode .eq. CenterofMass ) then
        allocate( this%MueX( np, nu ), STAT = stat )
        call AllocationError( stat, 'units*particles', nup )
        allocate( this%MueY( np, nu ), STAT = stat )
        call AllocationError( stat, 'units*particles', nup )
        allocate( this%MueZ( np, nu ), STAT = stat )
        call AllocationError( stat, 'units*particles', nup )

        if( ntest > 0 ) then
          allocate( this%MueXTest( ntest,nu ), STAT = stat )
          call AllocationError( stat, 'units*particles', nu*ntest )
          allocate( this%MueYTest( ntest,nu ), STAT = stat )
          call AllocationError( stat, 'units*particles', nu*ntest )
          allocate( this%MueZTest( ntest,nu ), STAT = stat )
          call AllocationError( stat, 'units*particles', nu*ntest )
        end if
      end if

    end if

    ! Gear corrector local arrays
    if( SimulationType .eq. MolecularDynamics .and. IntegratorType .eq. IntegratorTypeGear ) then
      allocate( this%Corr0( np, merge( 4, 3, this%Molecule%isElongated ), nu),STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%Corr1( np, merge( 4, 3, this%Molecule%isElongated ), nu),STAT = stat )
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
      this%Molecule%SiteMIEnm(i)%NTest0 => this%NTest0
      this%Molecule%SiteMIEnm(i)%NTest1 => this%NTest1
      this%Molecule%SiteMIEnm(i)%NTest2 => this%NTest2

      call Allocate( this%Molecule%SiteMIEnm(i) )
      this%Molecule%SiteMIEnm(i)%PX => this%Pm0(:, 1)
      this%Molecule%SiteMIEnm(i)%PY => this%Pm0(:, 2)
      this%Molecule%SiteMIEnm(i)%PZ => this%Pm0(:, 3)

#if TRANS==1
      this%Molecule%SiteMIEnm(i)%Qm0r => this%Q0
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
      this%Molecule%SiteTT68(i)%PX => this%P0(:, 1, 1)
      this%Molecule%SiteTT68(i)%PY => this%P0(:, 2, 1)
      this%Molecule%SiteTT68(i)%PZ => this%P0(:, 3, 1)

      if( ntest > 0 ) then
        this%Molecule%SiteTT68(i)%PXTest => this%P0Test(:, 1, 1)
        this%Molecule%SiteTT68(i)%PYTest => this%P0Test(:, 2, 1)
        this%Molecule%SiteTT68(i)%PZTest => this%P0Test(:, 3, 1)
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
      this%Molecule%SiteCharge(i)%NTest0 => this%NTest0
      this%Molecule%SiteCharge(i)%NTest1 => this%NTest1
      this%Molecule%SiteCharge(i)%NTest2 => this%NTest2

      call Allocate( this%Molecule%SiteCharge(i) )
      this%Molecule%SiteCharge(i)%PX => this%Pm0(:, 1)
      this%Molecule%SiteCharge(i)%PY => this%Pm0(:, 2)
      this%Molecule%SiteCharge(i)%PZ => this%Pm0(:, 3)

#if TRANS==1
      this%Molecule%SiteCharge(i)%Qm0r => this%Q0
#endif
    end do

    do i = 1, this%Molecule%NDipole
      this%Molecule%SiteDipole(i)%NPartMax => this%NPartMax
      this%Molecule%SiteDipole(i)%NPart => this%NPart
      this%Molecule%SiteDipole(i)%NTest => this%NTest
      this%Molecule%SiteDipole(i)%NPart0 => this%NPart0
      this%Molecule%SiteDipole(i)%NPart1 => this%NPart1
      this%Molecule%SiteDipole(i)%NPart2 => this%NPart2
      this%Molecule%SiteDipole(i)%NTest0 => this%NTest0
      this%Molecule%SiteDipole(i)%NTest1 => this%NTest1
      this%Molecule%SiteDipole(i)%NTest2 => this%NTest2

      call Allocate( this%Molecule%SiteDipole(i) )
      this%Molecule%SiteDipole(i)%PX => this%Pm0(:, 1)
      this%Molecule%SiteDipole(i)%PY => this%Pm0(:, 2)
      this%Molecule%SiteDipole(i)%PZ => this%Pm0(:, 3)

#if TRANS==1
      this%Molecule%SiteDipole(i)%Qm0r => this%Q0
#endif
    end do

    do i = 1, this%Molecule%NQuadrupole
      this%Molecule%SiteQuadrupole(i)%NPartMax => this%NPartMax
      this%Molecule%SiteQuadrupole(i)%NPart => this%NPart
      this%Molecule%SiteQuadrupole(i)%NTest => this%NTest
      this%Molecule%SiteQuadrupole(i)%NPart0 => this%NPart0
      this%Molecule%SiteQuadrupole(i)%NPart1 => this%NPart1
      this%Molecule%SiteQuadrupole(i)%NPart2 => this%NPart2
      this%Molecule%SiteQuadrupole(i)%NTest0 => this%NTest0
      this%Molecule%SiteQuadrupole(i)%NTest1 => this%NTest1
      this%Molecule%SiteQuadrupole(i)%NTest2 => this%NTest2

      call Allocate( this%Molecule%SiteQuadrupole(i) )
      this%Molecule%SiteQuadrupole(i)%PX => this%Pm0(:, 1)
      this%Molecule%SiteQuadrupole(i)%PY => this%Pm0(:, 2)
      this%Molecule%SiteQuadrupole(i)%PZ => this%Pm0(:, 3)

#if TRANS==1
      this%Molecule%SiteQuadrupole(i)%Qm0r => this%Q0
#endif
    end do

    ! Internal degrees of freedom

    ! Units
    nlj=0
    nch=0
    ndi=0
    nqu=0
    do i = 1, this%Molecule%NUnit
      this%Molecule%Unit(i)%NPartMax => this%NPartMax
      this%Molecule%Unit(i)%NPart => this%NPart
      this%Molecule%Unit(i)%NPart0 => this%NPart0
      this%Molecule%Unit(i)%NPart1 => this%NPart1
      this%Molecule%Unit(i)%NPart2 => this%NPart2
      this%Molecule%Unit(i)%PX => this%P0(:, 1, i)
      this%Molecule%Unit(i)%PY => this%P0(:, 2, i)
      this%Molecule%Unit(i)%PZ => this%P0(:, 3, i)

      if (this%Molecule%Unit(i)%NMIEnm > 0) then
        do j = 1, this%Molecule%Unit(i)%NMIEnm
          nlj = nlj+1
          this%Molecule%Unit(i)%SiteMIEnm(j)%r=this%Molecule%SiteMIEnm(nlj)%r
          this%Molecule%Unit(i)%SiteMIEnm(j)%RX=>this%Molecule%SiteMIEnm(nlj)%RX
          this%Molecule%Unit(i)%SiteMIEnm(j)%RY=>this%Molecule%SiteMIEnm(nlj)%RY
          this%Molecule%Unit(i)%SiteMIEnm(j)%RZ=>this%Molecule%SiteMIEnm(nlj)%RZ
          if (ntest>0) then
            this%Molecule%Unit(i)%SiteMIEnm(j)%RXTest=>this%Molecule%SiteMIEnm(nlj)%RXTest
            this%Molecule%Unit(i)%SiteMIEnm(j)%RYTest=>this%Molecule%SiteMIEnm(nlj)%RYTest
            this%Molecule%Unit(i)%SiteMIEnm(j)%RZTest=>this%Molecule%SiteMIEnm(nlj)%RZTest
            this%Molecule%Unit(i)%SiteMIEnm(j)%PXTest => this%P0Test(:, 1, i)
            this%Molecule%Unit(i)%SiteMIEnm(j)%PYTest => this%P0Test(:, 2, i)
            this%Molecule%Unit(i)%SiteMIEnm(j)%PZTest => this%P0Test(:, 3, i)
            this%Molecule%SiteMIEnm(nlj)%PXTest => this%Molecule%Unit(i)%SiteMIEnm(j)%PXTest
            this%Molecule%SiteMIEnm(nlj)%PYTest => this%Molecule%Unit(i)%SiteMIEnm(j)%PYTest
            this%Molecule%SiteMIEnm(nlj)%PZTest => this%Molecule%Unit(i)%SiteMIEnm(j)%PZTest
          endif
          this%Molecule%Unit(i)%SiteMIEnm(j)%FX=>this%Molecule%SiteMIEnm(nlj)%FX
          this%Molecule%Unit(i)%SiteMIEnm(j)%FY=>this%Molecule%SiteMIEnm(nlj)%FY
          this%Molecule%Unit(i)%SiteMIEnm(j)%FZ=>this%Molecule%SiteMIEnm(nlj)%FZ
          this%Molecule%SiteMIEnm(nlj)%PX =>this%Molecule%Unit(i)%PX
          this%Molecule%SiteMIEnm(nlj)%PY =>this%Molecule%Unit(i)%PY
          this%Molecule%SiteMIEnm(nlj)%PZ =>this%Molecule%Unit(i)%PZ
        end do
      end if
      if (this%Molecule%Unit(i)%NCharge > 0) then
        do j = 1, this%Molecule%Unit(i)%NCharge
          nch = nch+1
          this%Molecule%Unit(i)%SiteCharge(j)%r=this%Molecule%SiteCharge(nch)%r
          this%Molecule%Unit(i)%SiteCharge(j)%RX=>this%Molecule%SiteCharge(nch)%RX
          this%Molecule%Unit(i)%SiteCharge(j)%RY=>this%Molecule%SiteCharge(nch)%RY
          this%Molecule%Unit(i)%SiteCharge(j)%RZ=>this%Molecule%SiteCharge(nch)%RZ
          if (ntest>0) then
            this%Molecule%Unit(i)%SiteCharge(j)%RXTest=>this%Molecule%SiteCharge(nch)%RXTest
            this%Molecule%Unit(i)%SiteCharge(j)%RYTest=>this%Molecule%SiteCharge(nch)%RYTest
            this%Molecule%Unit(i)%SiteCharge(j)%RZTest=>this%Molecule%SiteCharge(nch)%RZTest
            this%Molecule%Unit(i)%SiteCharge(j)%PXTest => this%P0Test(:, 1, i)
            this%Molecule%Unit(i)%SiteCharge(j)%PYTest => this%P0Test(:, 2, i)
            this%Molecule%Unit(i)%SiteCharge(j)%PZTest => this%P0Test(:, 3, i)
            this%Molecule%SiteCharge(nch)%PXTest => this%Molecule%Unit(i)%SiteCharge(j)%PXTest
            this%Molecule%SiteCharge(nch)%PYTest => this%Molecule%Unit(i)%SiteCharge(j)%PYTest
            this%Molecule%SiteCharge(nch)%PZTest => this%Molecule%Unit(i)%SiteCharge(j)%PZTest
          endif
          this%Molecule%Unit(i)%SiteCharge(j)%FX=>this%Molecule%SiteCharge(nch)%FX
          this%Molecule%Unit(i)%SiteCharge(j)%FY=>this%Molecule%SiteCharge(nch)%FY
          this%Molecule%Unit(i)%SiteCharge(j)%FZ=>this%Molecule%SiteCharge(nch)%FZ
          this%Molecule%SiteCharge(nch)%PX => this%Molecule%Unit(i)%PX
          this%Molecule%SiteCharge(nch)%PY => this%Molecule%Unit(i)%PY
          this%Molecule%SiteCharge(nch)%PZ => this%Molecule%Unit(i)%PZ
        end do
      end if
      if (this%Molecule%Unit(i)%NDipole > 0) then
        do j = 1, this%Molecule%Unit(i)%NDipole
          ndi = ndi+1
          this%Molecule%Unit(i)%SiteDipole(j)%r=this%Molecule%SiteDipole(ndi)%r
          this%Molecule%Unit(i)%SiteDipole(j)%or=this%Molecule%SiteDipole(ndi)%or
          this%Molecule%Unit(i)%SiteDipole(j)%RX=>this%Molecule%SiteDipole(ndi)%RX
          this%Molecule%Unit(i)%SiteDipole(j)%RY=>this%Molecule%SiteDipole(ndi)%RY
          this%Molecule%Unit(i)%SiteDipole(j)%RZ=>this%Molecule%SiteDipole(ndi)%RZ
          this%Molecule%Unit(i)%SiteDipole(j)%OX=>this%Molecule%SiteDipole(ndi)%OX
          this%Molecule%Unit(i)%SiteDipole(j)%OY=>this%Molecule%SiteDipole(ndi)%OY
          this%Molecule%Unit(i)%SiteDipole(j)%OZ=>this%Molecule%SiteDipole(ndi)%OZ
          if (ntest>0) then
            this%Molecule%Unit(i)%SiteDipole(j)%RXTest=>this%Molecule%SiteDipole(ndi)%RXTest
            this%Molecule%Unit(i)%SiteDipole(j)%RYTest=>this%Molecule%SiteDipole(ndi)%RYTest
            this%Molecule%Unit(i)%SiteDipole(j)%RZTest=>this%Molecule%SiteDipole(ndi)%RZTest
            this%Molecule%Unit(i)%SiteDipole(j)%OXTest=>this%Molecule%SiteDipole(ndi)%OXTest
            this%Molecule%Unit(i)%SiteDipole(j)%OYTest=>this%Molecule%SiteDipole(ndi)%OYTest
            this%Molecule%Unit(i)%SiteDipole(j)%OZTest=>this%Molecule%SiteDipole(ndi)%OZTest
            this%Molecule%Unit(i)%SiteDipole(j)%PXTest => this%P0Test(:, 1, i)
            this%Molecule%Unit(i)%SiteDipole(j)%PYTest => this%P0Test(:, 2, i)
            this%Molecule%Unit(i)%SiteDipole(j)%PZTest => this%P0Test(:, 3, i)
            this%Molecule%SiteDipole(ndi)%PXTest => this%Molecule%Unit(i)%SiteDipole(j)%PXTest
            this%Molecule%SiteDipole(ndi)%PYTest => this%Molecule%Unit(i)%SiteDipole(j)%PYTest
            this%Molecule%SiteDipole(ndi)%PZTest => this%Molecule%Unit(i)%SiteDipole(j)%PZTest
          endif
          this%Molecule%Unit(i)%SiteDipole(j)%FX=>this%Molecule%SiteDipole(ndi)%FX
          this%Molecule%Unit(i)%SiteDipole(j)%FY=>this%Molecule%SiteDipole(ndi)%FY
          this%Molecule%Unit(i)%SiteDipole(j)%FZ=>this%Molecule%SiteDipole(ndi)%FZ
          this%Molecule%Unit(i)%SiteDipole(j)%TX=>this%Molecule%SiteDipole(ndi)%TX
          this%Molecule%Unit(i)%SiteDipole(j)%TY=>this%Molecule%SiteDipole(ndi)%TY
          this%Molecule%Unit(i)%SiteDipole(j)%TZ=>this%Molecule%SiteDipole(ndi)%TZ
          this%Molecule%SiteDipole(ndi)%PX=> this%Molecule%Unit(i)%PX
          this%Molecule%SiteDipole(ndi)%PY=> this%Molecule%Unit(i)%PY
          this%Molecule%SiteDipole(ndi)%PZ=> this%Molecule%Unit(i)%PZ
        end do
      end if
      if (this%Molecule%Unit(i)%NQuadrupole > 0) then
        do j = 1, this%Molecule%Unit(i)%NQuadrupole
          nqu = nqu+1
          this%Molecule%Unit(i)%SiteQuadrupole(j)%r=this%Molecule%SiteQuadrupole(nqu)%r
          this%Molecule%Unit(i)%SiteQuadrupole(j)%or=this%Molecule%SiteQuadrupole(nqu)%or
          this%Molecule%Unit(i)%SiteQuadrupole(j)%RX=>this%Molecule%SiteQuadrupole(nqu)%RX
          this%Molecule%Unit(i)%SiteQuadrupole(j)%RY=>this%Molecule%SiteQuadrupole(nqu)%RY
          this%Molecule%Unit(i)%SiteQuadrupole(j)%RZ=>this%Molecule%SiteQuadrupole(nqu)%RZ
          this%Molecule%Unit(i)%SiteQuadrupole(j)%OX=>this%Molecule%SiteQuadrupole(nqu)%OX
          this%Molecule%Unit(i)%SiteQuadrupole(j)%OY=>this%Molecule%SiteQuadrupole(nqu)%OY
          this%Molecule%Unit(i)%SiteQuadrupole(j)%OZ=>this%Molecule%SiteQuadrupole(nqu)%OZ
          if (ntest>0) then
            this%Molecule%Unit(i)%SiteQuadrupole(j)%RXTest=>this%Molecule%SiteQuadrupole(nqu)%RXTest
            this%Molecule%Unit(i)%SiteQuadrupole(j)%RYTest=>this%Molecule%SiteQuadrupole(nqu)%RYTest
            this%Molecule%Unit(i)%SiteQuadrupole(j)%RZTest=>this%Molecule%SiteQuadrupole(nqu)%RZTest
            this%Molecule%Unit(i)%SiteQuadrupole(j)%OXTest=>this%Molecule%SiteQuadrupole(nqu)%OXTest
            this%Molecule%Unit(i)%SiteQuadrupole(j)%OYTest=>this%Molecule%SiteQuadrupole(nqu)%OYTest
            this%Molecule%Unit(i)%SiteQuadrupole(j)%OZTest=>this%Molecule%SiteQuadrupole(nqu)%OZTest
            this%Molecule%Unit(i)%SiteQuadrupole(j)%PXTest => this%P0Test(:, 1, i)
            this%Molecule%Unit(i)%SiteQuadrupole(j)%PYTest => this%P0Test(:, 2, i)
            this%Molecule%Unit(i)%SiteQuadrupole(j)%PZTest => this%P0Test(:, 3, i)
            this%Molecule%SiteQuadrupole(nqu)%PXTest => this%Molecule%Unit(i)%SiteQuadrupole(j)%PXTest
            this%Molecule%SiteQuadrupole(nqu)%PYTest => this%Molecule%Unit(i)%SiteQuadrupole(j)%PYTest
            this%Molecule%SiteQuadrupole(nqu)%PZTest => this%Molecule%Unit(i)%SiteQuadrupole(j)%PZTest
          endif
          this%Molecule%Unit(i)%SiteQuadrupole(j)%FX=>this%Molecule%SiteQuadrupole(nqu)%FX
          this%Molecule%Unit(i)%SiteQuadrupole(j)%FY=>this%Molecule%SiteQuadrupole(nqu)%FY
          this%Molecule%Unit(i)%SiteQuadrupole(j)%FZ=>this%Molecule%SiteQuadrupole(nqu)%FZ
          this%Molecule%Unit(i)%SiteQuadrupole(j)%TX=>this%Molecule%SiteQuadrupole(nqu)%TX
          this%Molecule%Unit(i)%SiteQuadrupole(j)%TY=>this%Molecule%SiteQuadrupole(nqu)%TY
          this%Molecule%Unit(i)%SiteQuadrupole(j)%TZ=>this%Molecule%SiteQuadrupole(nqu)%TZ
          this%Molecule%SiteQuadrupole(nqu)%PX=> this%Molecule%Unit(i)%PX
          this%Molecule%SiteQuadrupole(nqu)%PY=> this%Molecule%Unit(i)%PY
          this%Molecule%SiteQuadrupole(nqu)%PZ=> this%Molecule%Unit(i)%PZ
        end do
      end if
    end do

    ! Idf Site positions and  forces
    do i = 1, this%Molecule%NBond
      this%Molecule%IdfBond(i)%NPartMax => this%NPartMax
      this%Molecule%IdfBond(i)%NPart => this%NPart
      this%Molecule%IdfBond(i)%NPart0 => this%NPart0
      this%Molecule%IdfBond(i)%NPart1 => this%NPart1
      this%Molecule%IdfBond(i)%NPart2 => this%NPart2
      SiteId1 = this%Molecule%IdfBond(i)%SiteId1
      SiteId2 = this%Molecule%IdfBond(i)%SiteId2
      Site1 = .false.
      Site2 = .false.
      if ( this%Molecule%NMIEnm>0 ) then
        do j = 1, this%Molecule%NMIEnm
          if (this%Molecule%SiteMIEnm(j)%SiteId==SiteId1) then
            this%Molecule%IdfBond(i)%RX1=>this%Molecule%SiteMIEnm(j)%RX(:)
            this%Molecule%IdfBond(i)%RY1=>this%Molecule%SiteMIEnm(j)%RY(:)
            this%Molecule%IdfBond(i)%RZ1=>this%Molecule%SiteMIEnm(j)%RZ(:)
            this%Molecule%IdfBond(i)%FX1=>this%Molecule%SiteMIEnm(j)%FX(:)
            this%Molecule%IdfBond(i)%FY1=>this%Molecule%SiteMIEnm(j)%FY(:)
            this%Molecule%IdfBond(i)%FZ1=>this%Molecule%SiteMIEnm(j)%FZ(:)
            this%Molecule%IdfBond(i)%PX1=>this%Molecule%SiteMIEnm(j)%PX(:)
            this%Molecule%IdfBond(i)%PY1=>this%Molecule%SiteMIEnm(j)%PY(:)
            this%Molecule%IdfBond(i)%PZ1=>this%Molecule%SiteMIEnm(j)%PZ(:)
            Site1 = .true.
          else if (this%Molecule%SiteMIEnm(j)%SiteId==SiteId2) then
            this%Molecule%IdfBond(i)%RX2=>this%Molecule%SiteMIEnm(j)%RX(:)
            this%Molecule%IdfBond(i)%RY2=>this%Molecule%SiteMIEnm(j)%RY(:)
            this%Molecule%IdfBond(i)%RZ2=>this%Molecule%SiteMIEnm(j)%RZ(:)
            this%Molecule%IdfBond(i)%FX2=>this%Molecule%SiteMIEnm(j)%FX(:)
            this%Molecule%IdfBond(i)%FY2=>this%Molecule%SiteMIEnm(j)%FY(:)
            this%Molecule%IdfBond(i)%FZ2=>this%Molecule%SiteMIEnm(j)%FZ(:)
            this%Molecule%IdfBond(i)%PX2=>this%Molecule%SiteMIEnm(j)%PX(:)
            this%Molecule%IdfBond(i)%PY2=>this%Molecule%SiteMIEnm(j)%PY(:)
            this%Molecule%IdfBond(i)%PZ2=>this%Molecule%SiteMIEnm(j)%PZ(:)
            Site2 = .true.
          end if
          if (Site1 .and. Site2) exit
        end do
      end if
      if((.not.Site1 .or. .not. Site2) .and. (this%Molecule%NCharge > 0) ) then
        do j = 1, this%Molecule%NCharge
          if (this%Molecule%SiteCharge(j)%SiteId==SiteId1) then
            this%Molecule%IdfBond(i)%RX1=>this%Molecule%SiteCharge(j)%RX(:)
            this%Molecule%IdfBond(i)%RY1=>this%Molecule%SiteCharge(j)%RY(:)
            this%Molecule%IdfBond(i)%RZ1=>this%Molecule%SiteCharge(j)%RZ(:)
            this%Molecule%IdfBond(i)%FX1=>this%Molecule%SiteCharge(j)%FX(:)
            this%Molecule%IdfBond(i)%FY1=>this%Molecule%SiteCharge(j)%FY(:)
            this%Molecule%IdfBond(i)%FZ1=>this%Molecule%SiteCharge(j)%FZ(:)
            this%Molecule%IdfBond(i)%PX1=>this%Molecule%SiteCharge(j)%PX(:)
            this%Molecule%IdfBond(i)%PY1=>this%Molecule%SiteCharge(j)%PY(:)
            this%Molecule%IdfBond(i)%PZ1=>this%Molecule%SiteCharge(j)%PZ(:)
            Site1 = .true.
          else if (this%Molecule%SiteCharge(j)%SiteId==SiteId2) then
            this%Molecule%IdfBond(i)%RX2=>this%Molecule%SiteCharge(j)%RX(:)
            this%Molecule%IdfBond(i)%RY2=>this%Molecule%SiteCharge(j)%RY(:)
            this%Molecule%IdfBond(i)%RZ2=>this%Molecule%SiteCharge(j)%RZ(:)
            this%Molecule%IdfBond(i)%FX2=>this%Molecule%SiteCharge(j)%FX(:)
            this%Molecule%IdfBond(i)%FY2=>this%Molecule%SiteCharge(j)%FY(:)
            this%Molecule%IdfBond(i)%FZ2=>this%Molecule%SiteCharge(j)%FZ(:)
            this%Molecule%IdfBond(i)%PX2=>this%Molecule%SiteCharge(j)%PX(:)
            this%Molecule%IdfBond(i)%PY2=>this%Molecule%SiteCharge(j)%PY(:)
            this%Molecule%IdfBond(i)%PZ2=>this%Molecule%SiteCharge(j)%PZ(:)
            Site2 = .true.
          end if
          if (Site1 .and. Site2) exit
        end do
      end if
      if((.not.Site1 .or. .not. Site2) .and. (this%Molecule%NDipole > 0) ) then
        do j = 1, this%Molecule%NDipole
          if (this%Molecule%SiteDipole(j)%SiteId==SiteId1) then
            this%Molecule%IdfBond(i)%RX1=>this%Molecule%SiteDipole(j)%RX(:)
            this%Molecule%IdfBond(i)%RY1=>this%Molecule%SiteDipole(j)%RY(:)
            this%Molecule%IdfBond(i)%RZ1=>this%Molecule%SiteDipole(j)%RZ(:)
            this%Molecule%IdfBond(i)%FX1=>this%Molecule%SiteDipole(j)%FX(:)
            this%Molecule%IdfBond(i)%FY1=>this%Molecule%SiteDipole(j)%FY(:)
            this%Molecule%IdfBond(i)%FZ1=>this%Molecule%SiteDipole(j)%FZ(:)
            this%Molecule%IdfBond(i)%PX1=>this%Molecule%SiteDipole(j)%PX(:)
            this%Molecule%IdfBond(i)%PY1=>this%Molecule%SiteDipole(j)%PY(:)
            this%Molecule%IdfBond(i)%PZ1=>this%Molecule%SiteDipole(j)%PZ(:)
            Site1 = .true.
          else if (this%Molecule%SiteDipole(j)%SiteId==SiteId2) then
            this%Molecule%IdfBond(i)%RX2=>this%Molecule%SiteDipole(j)%RX(:)
            this%Molecule%IdfBond(i)%RY2=>this%Molecule%SiteDipole(j)%RY(:)
            this%Molecule%IdfBond(i)%RZ2=>this%Molecule%SiteDipole(j)%RZ(:)
            this%Molecule%IdfBond(i)%FX2=>this%Molecule%SiteDipole(j)%FX(:)
            this%Molecule%IdfBond(i)%FY2=>this%Molecule%SiteDipole(j)%FY(:)
            this%Molecule%IdfBond(i)%FZ2=>this%Molecule%SiteDipole(j)%FZ(:)
            this%Molecule%IdfBond(i)%PX2=>this%Molecule%SiteDipole(j)%PX(:)
            this%Molecule%IdfBond(i)%PY2=>this%Molecule%SiteDipole(j)%PY(:)
            this%Molecule%IdfBond(i)%PZ2=>this%Molecule%SiteDipole(j)%PZ(:)
            Site2 = .true.
          end if
          if (Site1 .and. Site2) exit
        end do
      end if
      if((.not.Site1 .or. .not. Site2) .and. (this%Molecule%NQuadrupole > 0) ) then
        do j = 1, this%Molecule%NQuadrupole
          if (this%Molecule%SiteQuadrupole(j)%SiteId==SiteId1) then
            this%Molecule%IdfBond(i)%RX1=>this%Molecule%SiteQuadrupole(j)%RX(:)
            this%Molecule%IdfBond(i)%RY1=>this%Molecule%SiteQuadrupole(j)%RY(:)
            this%Molecule%IdfBond(i)%RZ1=>this%Molecule%SiteQuadrupole(j)%RZ(:)
            this%Molecule%IdfBond(i)%FX1=>this%Molecule%SiteQuadrupole(j)%FX(:)
            this%Molecule%IdfBond(i)%FY1=>this%Molecule%SiteQuadrupole(j)%FY(:)
            this%Molecule%IdfBond(i)%FZ1=>this%Molecule%SiteQuadrupole(j)%FZ(:)
            this%Molecule%IdfBond(i)%PX1=>this%Molecule%SiteQuadrupole(j)%PX(:)
            this%Molecule%IdfBond(i)%PY1=>this%Molecule%SiteQuadrupole(j)%PY(:)
            this%Molecule%IdfBond(i)%PZ1=>this%Molecule%SiteQuadrupole(j)%PZ(:)
            Site1 = .true.
          else if (this%Molecule%SiteQuadrupole(j)%SiteId==SiteId2) then
            this%Molecule%IdfBond(i)%RX2=>this%Molecule%SiteQuadrupole(j)%RX(:)
            this%Molecule%IdfBond(i)%RY2=>this%Molecule%SiteQuadrupole(j)%RY(:)
            this%Molecule%IdfBond(i)%RZ2=>this%Molecule%SiteQuadrupole(j)%RZ(:)
            this%Molecule%IdfBond(i)%FX2=>this%Molecule%SiteQuadrupole(j)%FX(:)
            this%Molecule%IdfBond(i)%FY2=>this%Molecule%SiteQuadrupole(j)%FY(:)
            this%Molecule%IdfBond(i)%FZ2=>this%Molecule%SiteQuadrupole(j)%FZ(:)
            this%Molecule%IdfBond(i)%PX2=>this%Molecule%SiteQuadrupole(j)%PX(:)
            this%Molecule%IdfBond(i)%PY2=>this%Molecule%SiteQuadrupole(j)%PY(:)
            this%Molecule%IdfBond(i)%PZ2=>this%Molecule%SiteQuadrupole(j)%PZ(:)
            Site2 = .true.
          end if
          if (Site1 .and. Site2) exit
        end do
      end if
    end do
    do i = 1, this%Molecule%NAngle
      this%Molecule%IdfAngle(i)%NPartMax => this%NPartMax
      this%Molecule%IdfAngle(i)%NPart => this%NPart
      this%Molecule%IdfAngle(i)%NPart0 => this%NPart0
      this%Molecule%IdfAngle(i)%NPart1 => this%NPart1
      this%Molecule%IdfAngle(i)%NPart2 => this%NPart2
      SiteId1 = this%Molecule%IdfAngle(i)%SiteId1
      SiteId2 = this%Molecule%IdfAngle(i)%SiteId2
      SiteId3 = this%Molecule%IdfAngle(i)%SiteId3
      Site1 = .false.
      Site2 = .false.
      Site3 = .false.
      if ( this%Molecule%NMIEnm>0 ) then
        do j = 1, this%Molecule%NMIEnm
          if (this%Molecule%SiteMIEnm(j)%SiteId==SiteId1) then
            this%Molecule%IdfAngle(i)%RX1=>this%Molecule%SiteMIEnm(j)%RX(:)
            this%Molecule%IdfAngle(i)%RY1=>this%Molecule%SiteMIEnm(j)%RY(:)
            this%Molecule%IdfAngle(i)%RZ1=>this%Molecule%SiteMIEnm(j)%RZ(:)
            this%Molecule%IdfAngle(i)%FX1=>this%Molecule%SiteMIEnm(j)%FX(:)
            this%Molecule%IdfAngle(i)%FY1=>this%Molecule%SiteMIEnm(j)%FY(:)
            this%Molecule%IdfAngle(i)%FZ1=>this%Molecule%SiteMIEnm(j)%FZ(:)
            Site1 = .true.
          else if (this%Molecule%SiteMIEnm(j)%SiteId==SiteId2) then
            this%Molecule%IdfAngle(i)%RX2=>this%Molecule%SiteMIEnm(j)%RX(:)
            this%Molecule%IdfAngle(i)%RY2=>this%Molecule%SiteMIEnm(j)%RY(:)
            this%Molecule%IdfAngle(i)%RZ2=>this%Molecule%SiteMIEnm(j)%RZ(:)
            this%Molecule%IdfAngle(i)%FX2=>this%Molecule%SiteMIEnm(j)%FX(:)
            this%Molecule%IdfAngle(i)%FY2=>this%Molecule%SiteMIEnm(j)%FY(:)
            this%Molecule%IdfAngle(i)%FZ2=>this%Molecule%SiteMIEnm(j)%FZ(:)
            Site2 = .true.
          else if (this%Molecule%SiteMIEnm(j)%SiteId==SiteId3) then
            this%Molecule%IdfAngle(i)%RX3=>this%Molecule%SiteMIEnm(j)%RX(:)
            this%Molecule%IdfAngle(i)%RY3=>this%Molecule%SiteMIEnm(j)%RY(:)
            this%Molecule%IdfAngle(i)%RZ3=>this%Molecule%SiteMIEnm(j)%RZ(:)
            this%Molecule%IdfAngle(i)%FX3=>this%Molecule%SiteMIEnm(j)%FX(:)
            this%Molecule%IdfAngle(i)%FY3=>this%Molecule%SiteMIEnm(j)%FY(:)
            this%Molecule%IdfAngle(i)%FZ3=>this%Molecule%SiteMIEnm(j)%FZ(:)
            Site3 = .true.
          end if
          if (Site1 .and. Site2 .and. Site3) exit
        end do
      end if
      if((.not.Site1.or. .not.Site2 .or. .not.Site3) .and. (this%Molecule%NCharge>0) ) then
        do j = 1, this%Molecule%NCharge
          if (this%Molecule%SiteCharge(j)%SiteId==SiteId1) then
            this%Molecule%IdfAngle(i)%RX1=>this%Molecule%SiteCharge(j)%RX(:)
            this%Molecule%IdfAngle(i)%RY1=>this%Molecule%SiteCharge(j)%RY(:)
            this%Molecule%IdfAngle(i)%RZ1=>this%Molecule%SiteCharge(j)%RZ(:)
            this%Molecule%IdfAngle(i)%FX1=>this%Molecule%SiteCharge(j)%FX(:)
            this%Molecule%IdfAngle(i)%FY1=>this%Molecule%SiteCharge(j)%FY(:)
            this%Molecule%IdfAngle(i)%FZ1=>this%Molecule%SiteCharge(j)%FZ(:)
            Site1 = .true.
          else if (this%Molecule%SiteCharge(j)%SiteId==SiteId2) then
            this%Molecule%IdfAngle(i)%RX2=>this%Molecule%SiteCharge(j)%RX(:)
            this%Molecule%IdfAngle(i)%RY2=>this%Molecule%SiteCharge(j)%RY(:)
            this%Molecule%IdfAngle(i)%RZ2=>this%Molecule%SiteCharge(j)%RZ(:)
            this%Molecule%IdfAngle(i)%FX2=>this%Molecule%SiteCharge(j)%FX(:)
            this%Molecule%IdfAngle(i)%FY2=>this%Molecule%SiteCharge(j)%FY(:)
            this%Molecule%IdfAngle(i)%FZ2=>this%Molecule%SiteCharge(j)%FZ(:)
            Site2 = .true.
          else if (this%Molecule%SiteCharge(j)%SiteId==SiteId3) then
            this%Molecule%IdfAngle(i)%RX3=>this%Molecule%SiteCharge(j)%RX(:)
            this%Molecule%IdfAngle(i)%RY3=>this%Molecule%SiteCharge(j)%RY(:)
            this%Molecule%IdfAngle(i)%RZ3=>this%Molecule%SiteCharge(j)%RZ(:)
            this%Molecule%IdfAngle(i)%FX3=>this%Molecule%SiteCharge(j)%FX(:)
            this%Molecule%IdfAngle(i)%FY3=>this%Molecule%SiteCharge(j)%FY(:)
            this%Molecule%IdfAngle(i)%FZ3=>this%Molecule%SiteCharge(j)%FZ(:)
            Site3 = .true.
          end if
          if (Site1 .and. Site2 .and. Site3) exit
        end do
      end if
      if((.not.Site1.or. .not.Site2 .or. .not.Site3) .and. (this%Molecule%NDipole>0) ) then
        do j = 1, this%Molecule%NDipole
          if (this%Molecule%SiteDipole(j)%SiteId==SiteId1) then
            if (SiteId1 == SiteId2) then
              this%Molecule%IdfAngle(i)%RX1=>this%Molecule%SiteDipole(j)%OX(:)
              this%Molecule%IdfAngle(i)%RY1=>this%Molecule%SiteDipole(j)%OY(:)
              this%Molecule%IdfAngle(i)%RZ1=>this%Molecule%SiteDipole(j)%OZ(:)
            else
              this%Molecule%IdfAngle(i)%RX1=>this%Molecule%SiteDipole(j)%RX(:)
              this%Molecule%IdfAngle(i)%RY1=>this%Molecule%SiteDipole(j)%RY(:)
              this%Molecule%IdfAngle(i)%RZ1=>this%Molecule%SiteDipole(j)%RZ(:)
            end if
            this%Molecule%IdfAngle(i)%FX1=>this%Molecule%SiteDipole(j)%FX(:)
            this%Molecule%IdfAngle(i)%FY1=>this%Molecule%SiteDipole(j)%FY(:)
            this%Molecule%IdfAngle(i)%FZ1=>this%Molecule%SiteDipole(j)%FZ(:)
            Site1 = .true.
          else if (this%Molecule%SiteDipole(j)%SiteId==SiteId2) then
            this%Molecule%IdfAngle(i)%RX2=>this%Molecule%SiteDipole(j)%RX(:)
            this%Molecule%IdfAngle(i)%RY2=>this%Molecule%SiteDipole(j)%RY(:)
            this%Molecule%IdfAngle(i)%RZ2=>this%Molecule%SiteDipole(j)%RZ(:)
            this%Molecule%IdfAngle(i)%FX2=>this%Molecule%SiteDipole(j)%FX(:)
            this%Molecule%IdfAngle(i)%FY2=>this%Molecule%SiteDipole(j)%FY(:)
            this%Molecule%IdfAngle(i)%FZ2=>this%Molecule%SiteDipole(j)%FZ(:)
            Site2 = .true.
          else if (this%Molecule%SiteDipole(j)%SiteId==SiteId3) then
            if (SiteId3 == SiteId2) then
              this%Molecule%IdfAngle(i)%RX3=>this%Molecule%SiteDipole(j)%OX(:)
              this%Molecule%IdfAngle(i)%RY3=>this%Molecule%SiteDipole(j)%OY(:)
              this%Molecule%IdfAngle(i)%RZ3=>this%Molecule%SiteDipole(j)%OZ(:)
            else
              this%Molecule%IdfAngle(i)%RX3=>this%Molecule%SiteDipole(j)%RX(:)
              this%Molecule%IdfAngle(i)%RY3=>this%Molecule%SiteDipole(j)%RY(:)
              this%Molecule%IdfAngle(i)%RZ3=>this%Molecule%SiteDipole(j)%RZ(:)
            end if
            this%Molecule%IdfAngle(i)%FX3=>this%Molecule%SiteDipole(j)%FX(:)
            this%Molecule%IdfAngle(i)%FY3=>this%Molecule%SiteDipole(j)%FY(:)
            this%Molecule%IdfAngle(i)%FZ3=>this%Molecule%SiteDipole(j)%FZ(:)
            Site3 = .true.
          end if
          if (Site1 .and. Site2 .and. Site3) exit
        end do
      end if
      if((.not.Site1.or. .not.Site2 .or. .not.Site3) .and. (this%Molecule%NQuadrupole>0) ) then
        do j = 1, this%Molecule%NQuadrupole
          if (this%Molecule%SiteQuadrupole(j)%SiteId==SiteId1) then
            if (SiteId1 == SiteId2) then
              this%Molecule%IdfAngle(i)%RX1=>this%Molecule%SiteQuadrupole(j)%OX(:)
              this%Molecule%IdfAngle(i)%RY1=>this%Molecule%SiteQuadrupole(j)%OY(:)
              this%Molecule%IdfAngle(i)%RZ1=>this%Molecule%SiteQuadrupole(j)%OZ(:)
            else
              this%Molecule%IdfAngle(i)%RX1=>this%Molecule%SiteQuadrupole(j)%RX(:)
              this%Molecule%IdfAngle(i)%RY1=>this%Molecule%SiteQuadrupole(j)%RY(:)
              this%Molecule%IdfAngle(i)%RZ1=>this%Molecule%SiteQuadrupole(j)%RZ(:)
            end if
            this%Molecule%IdfAngle(i)%FX1=>this%Molecule%SiteQuadrupole(j)%FX(:)
            this%Molecule%IdfAngle(i)%FY1=>this%Molecule%SiteQuadrupole(j)%FY(:)
            this%Molecule%IdfAngle(i)%FZ1=>this%Molecule%SiteQuadrupole(j)%FZ(:)
            Site1 = .true.
          else if (this%Molecule%SiteQuadrupole(j)%SiteId==SiteId2) then
            this%Molecule%IdfAngle(i)%RX2=>this%Molecule%SiteQuadrupole(j)%RX(:)
            this%Molecule%IdfAngle(i)%RY2=>this%Molecule%SiteQuadrupole(j)%RY(:)
            this%Molecule%IdfAngle(i)%RZ2=>this%Molecule%SiteQuadrupole(j)%RZ(:)
            this%Molecule%IdfAngle(i)%FX2=>this%Molecule%SiteQuadrupole(j)%FX(:)
            this%Molecule%IdfAngle(i)%FY2=>this%Molecule%SiteQuadrupole(j)%FY(:)
            this%Molecule%IdfAngle(i)%FZ2=>this%Molecule%SiteQuadrupole(j)%FZ(:)
            Site2 = .true.
          else if (this%Molecule%SiteQuadrupole(j)%SiteId==SiteId3) then
            if (SiteId3 == SiteId2) then
              this%Molecule%IdfAngle(i)%RX3=>this%Molecule%SiteQuadrupole(j)%OX(:)
              this%Molecule%IdfAngle(i)%RY3=>this%Molecule%SiteQuadrupole(j)%OY(:)
              this%Molecule%IdfAngle(i)%RZ3=>this%Molecule%SiteQuadrupole(j)%OZ(:)
            else
              this%Molecule%IdfAngle(i)%RX3=>this%Molecule%SiteQuadrupole(j)%RX(:)
              this%Molecule%IdfAngle(i)%RY3=>this%Molecule%SiteQuadrupole(j)%RY(:)
              this%Molecule%IdfAngle(i)%RZ3=>this%Molecule%SiteQuadrupole(j)%RZ(:)
            end if
            this%Molecule%IdfAngle(i)%FX3=>this%Molecule%SiteQuadrupole(j)%FX(:)
            this%Molecule%IdfAngle(i)%FY3=>this%Molecule%SiteQuadrupole(j)%FY(:)
            this%Molecule%IdfAngle(i)%FZ3=>this%Molecule%SiteQuadrupole(j)%FZ(:)
            Site3 = .true.
          end if
          if (Site1 .and. Site2 .and. Site3) exit
        end do
      end if
    end do
    do i = 1, this%Molecule%NDihedral
      this%Molecule%IdfDihedral(i)%NPartMax => this%NPartMax
      this%Molecule%IdfDihedral(i)%NPart => this%NPart
      this%Molecule%IdfDihedral(i)%NPart0 => this%NPart0
      this%Molecule%IdfDihedral(i)%NPart1 => this%NPart1
      this%Molecule%IdfDihedral(i)%NPart2 => this%NPart2
      SiteId1 = this%Molecule%IdfDihedral(i)%SiteId1
      SiteId2 = this%Molecule%IdfDihedral(i)%SiteId2
      SiteId3 = this%Molecule%IdfDihedral(i)%SiteId3
      SiteId4 = this%Molecule%IdfDihedral(i)%SiteId4
      Site1 = .false.
      Site2 = .false.
      Site3 = .false.
      Site4 = .false.
      if ( this%Molecule%NMIEnm>0 ) then
        do j = 1, this%Molecule%NMIEnm
          if (this%Molecule%SiteMIEnm(j)%SiteId==SiteId1) then
            this%Molecule%IdfDihedral(i)%RX1=>this%Molecule%SiteMIEnm(j)%RX(:)
            this%Molecule%IdfDihedral(i)%RY1=>this%Molecule%SiteMIEnm(j)%RY(:)
            this%Molecule%IdfDihedral(i)%RZ1=>this%Molecule%SiteMIEnm(j)%RZ(:)
            this%Molecule%IdfDihedral(i)%FX1=>this%Molecule%SiteMIEnm(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY1=>this%Molecule%SiteMIEnm(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ1=>this%Molecule%SiteMIEnm(j)%FZ(:)
            Site1 = .true.
          else if (this%Molecule%SiteMIEnm(j)%SiteId==SiteId2) then
            this%Molecule%IdfDihedral(i)%RX2=>this%Molecule%SiteMIEnm(j)%RX(:)
            this%Molecule%IdfDihedral(i)%RY2=>this%Molecule%SiteMIEnm(j)%RY(:)
            this%Molecule%IdfDihedral(i)%RZ2=>this%Molecule%SiteMIEnm(j)%RZ(:)
            this%Molecule%IdfDihedral(i)%FX2=>this%Molecule%SiteMIEnm(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY2=>this%Molecule%SiteMIEnm(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ2=>this%Molecule%SiteMIEnm(j)%FZ(:)
            Site2 = .true.
          else if (this%Molecule%SiteMIEnm(j)%SiteId==SiteId3) then
            this%Molecule%IdfDihedral(i)%RX3=>this%Molecule%SiteMIEnm(j)%RX(:)
            this%Molecule%IdfDihedral(i)%RY3=>this%Molecule%SiteMIEnm(j)%RY(:)
            this%Molecule%IdfDihedral(i)%RZ3=>this%Molecule%SiteMIEnm(j)%RZ(:)
            this%Molecule%IdfDihedral(i)%FX3=>this%Molecule%SiteMIEnm(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY3=>this%Molecule%SiteMIEnm(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ3=>this%Molecule%SiteMIEnm(j)%FZ(:)
            Site3 = .true.
          else if (this%Molecule%SiteMIEnm(j)%SiteId==SiteId4) then
            this%Molecule%IdfDihedral(i)%RX4=>this%Molecule%SiteMIEnm(j)%RX(:)
            this%Molecule%IdfDihedral(i)%RY4=>this%Molecule%SiteMIEnm(j)%RY(:)
            this%Molecule%IdfDihedral(i)%RZ4=>this%Molecule%SiteMIEnm(j)%RZ(:)
            this%Molecule%IdfDihedral(i)%FX4=>this%Molecule%SiteMIEnm(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY4=>this%Molecule%SiteMIEnm(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ4=>this%Molecule%SiteMIEnm(j)%FZ(:)
            Site4 = .true.
          end if
          if (Site1 .and. Site2 .and. Site3 .and. Site4) exit
        end do
      end if
      if((.not.Site1 .or. .not. Site2 .or. .not. Site3 .or. .not. Site4) &
&              .and. (this%Molecule%NCharge > 0) ) then
        do j = 1, this%Molecule%NCharge
          if (this%Molecule%SiteCharge(j)%SiteId==SiteId1) then
            this%Molecule%IdfDihedral(i)%RX1=>this%Molecule%SiteCharge(j)%RX(:)
            this%Molecule%IdfDihedral(i)%RY1=>this%Molecule%SiteCharge(j)%RY(:)
            this%Molecule%IdfDihedral(i)%RZ1=>this%Molecule%SiteCharge(j)%RZ(:)
            this%Molecule%IdfDihedral(i)%FX1=>this%Molecule%SiteCharge(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY1=>this%Molecule%SiteCharge(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ1=>this%Molecule%SiteCharge(j)%FZ(:)
            Site1 = .true.
          else if (this%Molecule%SiteCharge(j)%SiteId==SiteId2) then
            this%Molecule%IdfDihedral(i)%RX2=>this%Molecule%SiteCharge(j)%RX(:)
            this%Molecule%IdfDihedral(i)%RY2=>this%Molecule%SiteCharge(j)%RY(:)
            this%Molecule%IdfDihedral(i)%RZ2=>this%Molecule%SiteCharge(j)%RZ(:)
            this%Molecule%IdfDihedral(i)%FX2=>this%Molecule%SiteCharge(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY2=>this%Molecule%SiteCharge(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ2=>this%Molecule%SiteCharge(j)%FZ(:)
            Site2 = .true.
          else if (this%Molecule%SiteCharge(j)%SiteId==SiteId3) then
            this%Molecule%IdfDihedral(i)%RX3=>this%Molecule%SiteCharge(j)%RX(:)
            this%Molecule%IdfDihedral(i)%RY3=>this%Molecule%SiteCharge(j)%RY(:)
            this%Molecule%IdfDihedral(i)%RZ3=>this%Molecule%SiteCharge(j)%RZ(:)
            this%Molecule%IdfDihedral(i)%FX3=>this%Molecule%SiteCharge(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY3=>this%Molecule%SiteCharge(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ3=>this%Molecule%SiteCharge(j)%FZ(:)
            Site3 = .true.
          else if (this%Molecule%SiteCharge(j)%SiteId==SiteId4) then
            this%Molecule%IdfDihedral(i)%RX4=>this%Molecule%SiteCharge(j)%RX(:)
            this%Molecule%IdfDihedral(i)%RY4=>this%Molecule%SiteCharge(j)%RY(:)
            this%Molecule%IdfDihedral(i)%RZ4=>this%Molecule%SiteCharge(j)%RZ(:)
            this%Molecule%IdfDihedral(i)%FX4=>this%Molecule%SiteCharge(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY4=>this%Molecule%SiteCharge(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ4=>this%Molecule%SiteCharge(j)%FZ(:)
            Site4 = .true.
          end if
          if (Site1 .and. Site2 .and. Site3 .and. Site4) exit
        end do
      end if
      if((.not.Site1 .or. .not. Site2 .or. .not. Site3 .or. .not. Site4) &
&              .and. (this%Molecule%NDipole > 0) ) then
        do j = 1, this%Molecule%NDipole
          if (this%Molecule%SiteDipole(j)%SiteId==SiteId1) then
            if ( SiteId1 == SiteId2 ) then
              this%Molecule%IdfDihedral(i)%RX1=>this%Molecule%SiteDipole(j)%OX(:)
              this%Molecule%IdfDihedral(i)%RY1=>this%Molecule%SiteDipole(j)%OY(:)
              this%Molecule%IdfDihedral(i)%RZ1=>this%Molecule%SiteDipole(j)%OZ(:)
            else
              this%Molecule%IdfDihedral(i)%RX1=>this%Molecule%SiteDipole(j)%RX(:)
              this%Molecule%IdfDihedral(i)%RY1=>this%Molecule%SiteDipole(j)%RY(:)
              this%Molecule%IdfDihedral(i)%RZ1=>this%Molecule%SiteDipole(j)%RZ(:)
            end if
            this%Molecule%IdfDihedral(i)%FX1=>this%Molecule%SiteDipole(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY1=>this%Molecule%SiteDipole(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ1=>this%Molecule%SiteDipole(j)%FZ(:)
            Site1 = .true.
          else if (this%Molecule%SiteDipole(j)%SiteId==SiteId2) then
            this%Molecule%IdfDihedral(i)%RX2=>this%Molecule%SiteDipole(j)%RX(:)
            this%Molecule%IdfDihedral(i)%RY2=>this%Molecule%SiteDipole(j)%RY(:)
            this%Molecule%IdfDihedral(i)%RZ2=>this%Molecule%SiteDipole(j)%RZ(:)
            this%Molecule%IdfDihedral(i)%FX2=>this%Molecule%SiteDipole(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY2=>this%Molecule%SiteDipole(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ2=>this%Molecule%SiteDipole(j)%FZ(:)
            Site2 = .true.
          else if (this%Molecule%SiteDipole(j)%SiteId==SiteId3) then
            this%Molecule%IdfDihedral(i)%RX3=>this%Molecule%SiteDipole(j)%RX(:)
            this%Molecule%IdfDihedral(i)%RY3=>this%Molecule%SiteDipole(j)%RY(:)
            this%Molecule%IdfDihedral(i)%RZ3=>this%Molecule%SiteDipole(j)%RZ(:)
            this%Molecule%IdfDihedral(i)%FX3=>this%Molecule%SiteDipole(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY3=>this%Molecule%SiteDipole(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ3=>this%Molecule%SiteDipole(j)%FZ(:)
            Site3 = .true.
          else if (this%Molecule%SiteDipole(j)%SiteId==SiteId4) then
            if ( SiteId4 == SiteId3 ) then
              this%Molecule%IdfDihedral(i)%RX4=>this%Molecule%SiteDipole(j)%OX(:)
              this%Molecule%IdfDihedral(i)%RY4=>this%Molecule%SiteDipole(j)%OY(:)
              this%Molecule%IdfDihedral(i)%RZ4=>this%Molecule%SiteDipole(j)%OZ(:)
            else
              this%Molecule%IdfDihedral(i)%RX4=>this%Molecule%SiteDipole(j)%RX(:)
              this%Molecule%IdfDihedral(i)%RY4=>this%Molecule%SiteDipole(j)%RY(:)
              this%Molecule%IdfDihedral(i)%RZ4=>this%Molecule%SiteDipole(j)%RZ(:)
            end if
            this%Molecule%IdfDihedral(i)%FX4=>this%Molecule%SiteDipole(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY4=>this%Molecule%SiteDipole(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ4=>this%Molecule%SiteDipole(j)%FZ(:)
            Site4 = .true.
          end if
          if (Site1 .and. Site2 .and. Site3 .and. Site4) exit
        end do
      end if
      if((.not.Site1 .or. .not. Site2 .or. .not. Site3 .or. .not. Site4) &
&              .and. (this%Molecule%NQuadrupole > 0) ) then
        do j = 1, this%Molecule%NQuadrupole
          if (this%Molecule%SiteQuadrupole(j)%SiteId==SiteId1) then
            if ( SiteId1 == SiteId2 ) then
              this%Molecule%IdfDihedral(i)%RX1=>this%Molecule%SiteQuadrupole(j)%OX(:)
              this%Molecule%IdfDihedral(i)%RY1=>this%Molecule%SiteQuadrupole(j)%OY(:)
              this%Molecule%IdfDihedral(i)%RZ1=>this%Molecule%SiteQuadrupole(j)%OZ(:)
            else
              this%Molecule%IdfDihedral(i)%RX1=>this%Molecule%SiteQuadrupole(j)%RX(:)
              this%Molecule%IdfDihedral(i)%RY1=>this%Molecule%SiteQuadrupole(j)%RY(:)
              this%Molecule%IdfDihedral(i)%RZ1=>this%Molecule%SiteQuadrupole(j)%RZ(:)
            end if
            this%Molecule%IdfDihedral(i)%FX1=>this%Molecule%SiteQuadrupole(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY1=>this%Molecule%SiteQuadrupole(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ1=>this%Molecule%SiteQuadrupole(j)%FZ(:)
            Site1 = .true.
          else if (this%Molecule%SiteQuadrupole(j)%SiteId==SiteId2) then
            this%Molecule%IdfDihedral(i)%RX2=>this%Molecule%SiteQuadrupole(j)%RX(:)
            this%Molecule%IdfDihedral(i)%RY2=>this%Molecule%SiteQuadrupole(j)%RY(:)
            this%Molecule%IdfDihedral(i)%RZ2=>this%Molecule%SiteQuadrupole(j)%RZ(:)
            this%Molecule%IdfDihedral(i)%FX2=>this%Molecule%SiteQuadrupole(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY2=>this%Molecule%SiteQuadrupole(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ2=>this%Molecule%SiteQuadrupole(j)%FZ(:)
            Site2 = .true.
          else if (this%Molecule%SiteQuadrupole(j)%SiteId==SiteId3) then
            this%Molecule%IdfDihedral(i)%RX3=>this%Molecule%SiteQuadrupole(j)%RX(:)
            this%Molecule%IdfDihedral(i)%RY3=>this%Molecule%SiteQuadrupole(j)%RY(:)
            this%Molecule%IdfDihedral(i)%RZ3=>this%Molecule%SiteQuadrupole(j)%RZ(:)
            this%Molecule%IdfDihedral(i)%FX3=>this%Molecule%SiteQuadrupole(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY3=>this%Molecule%SiteQuadrupole(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ3=>this%Molecule%SiteQuadrupole(j)%FZ(:)
            Site3 = .true.
          else if (this%Molecule%SiteQuadrupole(j)%SiteId==SiteId4) then
            if ( SiteId4 == SiteId3 ) then
              this%Molecule%IdfDihedral(i)%RX4=>this%Molecule%SiteQuadrupole(j)%OX(:)
              this%Molecule%IdfDihedral(i)%RY4=>this%Molecule%SiteQuadrupole(j)%OY(:)
              this%Molecule%IdfDihedral(i)%RZ4=>this%Molecule%SiteQuadrupole(j)%OZ(:)
            else
              this%Molecule%IdfDihedral(i)%RX4=>this%Molecule%SiteQuadrupole(j)%RX(:)
              this%Molecule%IdfDihedral(i)%RY4=>this%Molecule%SiteQuadrupole(j)%RY(:)
              this%Molecule%IdfDihedral(i)%RZ4=>this%Molecule%SiteQuadrupole(j)%RZ(:)
            end if
            this%Molecule%IdfDihedral(i)%FX4=>this%Molecule%SiteQuadrupole(j)%FX(:)
            this%Molecule%IdfDihedral(i)%FY4=>this%Molecule%SiteQuadrupole(j)%FY(:)
            this%Molecule%IdfDihedral(i)%FZ4=>this%Molecule%SiteQuadrupole(j)%FZ(:)
            Site4 = .true.
          end if
          if (Site1 .and. Site2 .and. Site3 .and. Site4) exit
        end do
      end if
    end do


    if (UseIntDegFreed) then
      this%BondCount => this%Molecule%BondCount
      this%BoPartner => this%Molecule%BoPartner
      this%AngleCount => this%Molecule%AngleCount
      this%AnglePartner => this%Molecule%AnglePartner
      this%DihedralCount => this%Molecule%DihedralCount
      this%DihedralPartner => this%Molecule%DihedralPartner
    end if

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

    if( associated( this%EPotTestIntra ) ) then
      deallocate( this%EPotTestIntra )
    end if

    ! Centers of mass positions and their derivatives
    if( associated( this%P0 ) ) then
      deallocate( this%P0 )
    end if
    if( associated( this%Pm0 ) ) then
      deallocate( this%Pm0 )
    end if
    if( associated( this%P0Save ) ) then
      deallocate( this%P0Save )
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

     if ( ((EnsembleType .eq. EnsembleTypeGE) .or. (EnsembleType .eq. EnsembleTypeMUVT) .or. (EnsembleType .eq. EnsembleTypeHA)) .and. (abs(q) .ge. 1e-1) ) then
       write (ErrorBuffer,'("GrandEquilibrium not possible in a charged system, q=",G16.9)') q
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
    integer :: iUnit

    ! Set random linear velocities
    if (EMinimizationIDF) then ! Michaael Sch.: changed here for preEquilibration
      do i = 1, 3
        do j = 1, this%NPart
          this%P1(j, i, :) = rnd( -1._RK, 1._RK )
        end do
      end do
    else
      do iUnit = 1, this%Molecule%NUnit
        do i = 1, 3
          do j = 1, this%NPart
            this%P1(j, i, iUnit) = rnd( -1._RK, 1._RK )
          end do
        end do
      end do
    end if

    ! Normalize translational velocity vectors (only done once - needs not to be efficient)
    do iUnit = 1, this%Molecule%NUnit
      do j = 1, this%NPart
        this%P1(j, :, iUnit) = this%P1(j, :, iUnit) / sqrt( dot_product( this%P1(j, :, iUnit), this%P1(j, :, iUnit) ))
      end do
    end do

    ! Nullify angular velocities
    if( this%Molecule%isElongated ) this%W0(:, :, :) = 0._RK

  end subroutine TComponent_InitVelocities


!==============================================================!
!  Subroutine TComponent_InitIntegratorGear                    !
!==============================================================!

  subroutine TComponent_InitIntegratorGear( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Zero accelerations
    this%P2(:, :, :) = 0._RK
    this%P3(:, :, :) = 0._RK
    this%P4(:, :, :) = 0._RK
    this%P5(:, :, :) = 0._RK
    if( this%Molecule%isElongated ) then
      this%Q1(:, :, :) = 0._RK
      this%Q2(:, :, :) = 0._RK
      this%Q3(:, :, :) = 0._RK
      this%Q4(:, :, :) = 0._RK
      this%W1(:, :, :) = 0._RK
      this%W2(:, :, :) = 0._RK
      this%W3(:, :, :) = 0._RK
      this%W4(:, :, :) = 0._RK
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
    this%P2(:, :, :) = 0._RK

    if( this%Molecule%isElongated ) then
      this%Q1(:, :, :) = 0._RK
      this%W1(:, :, :) = 0._RK
    end if

  end subroutine TComponent_InitIntegratorLeap


!==============================================================!
!  Subroutine TComponent_RemoveNetMomentum                     !
!==============================================================!

  subroutine TComponent_RemoveNetMomentum( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Declare local variables
    real(RK) :: P(3), L(3)
    integer :: i, j, k
    real(RK) :: Pim

    ! Return if zero particles in component
    if( this%NPart == 0 ) return

    ! Calculate net momentum
    do k = 1, this%Molecule%NUnit
      P(:) = 0._RK
      L(:) = 0._RK
      do i = 1, 3
        if (.not. UseIntDegFreed) then
            P(i) = P(i) + this%Molecule%Mass * sum( this%P1(1:this%NPart, i, k) )
        else
            P(i) = P(i) + this%Molecule%Unit(k)%Mass * sum( this%P1(1:this%NPart, i, k) )
        end if

        if( i <= this%Molecule%Unit(k)%NDFRot ) L(i) = L(i) + this%Molecule%Unit(k)%MOI(i) * sum( this%W0(1:this%NPart, i, k) )

      end do
      P(:) = P(:) / this%NPart
      L(:) = L(:) / this%NPart

      ! Remove net momentum
      do i = 1, 3
        if (.not. UseIntDegFreed) then
            Pim = P(i) / this%Molecule%Mass
        else
            Pim = P(i) / this%Molecule%Unit(k)%Mass
        end if

        do j = 1, this%NPart
          this%P1(j, i, k) = this%P1(j, i, k) - Pim
        end do

        if( i <= this%Molecule%Unit(k)%NDFRot ) then
          Pim = L(i) / this%Molecule%Unit(k)%MOI(i)
          do j = 1, this%NPart
            this%W0(j, i, k) = this%W0(j, i, k) - Pim
          end do
        end if
      end do
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
    integer :: i, iUnit

    this%EKinTran = 0._RK

    ! Calculate translational kinetic energy
    do  iUnit = 1,  this%Molecule%NUnit
      if (.not. UseIntDegFreed .and. .not. EMinimizationIDF) then
          this%EKinTran = this%EkinTran+this%Molecule%Mass * TimeStepSquaredInv2 &
&           * sum( this%P1(1:this%NPart, :, iUnit)**2 ) * this%BoxLength**2
      else
          this%EKinTran = this%EkinTran+this%Molecule%Unit(iUnit)%Mass * TimeStepSquaredInv2 &
&           * sum( this%P1(1:this%NPart, :, iUnit)**2 ) * this%BoxLength**2
      end if
    end do

    ! Calculate rotational kinetic energy
    this%EKinRot = 0._RK
    do iUnit = 1, this%Molecule%NUnit
      do i = 1, this%Molecule%Unit(iUnit)%NDFRot
        this%EKinRot = this%EKinRot + this%Molecule%Unit(iUnit)%MOI(i) * .5_RK &
&         * sum( this%W0(1:this%NPart, i, iUnit)**2 )
      end do
    end do

  end subroutine TComponent_CalculateEKin


!==============================================================!
!  Subroutine TComponent_Unit2Atom                             !
!==============================================================!

!   subroutine TComponent_Mol2Atom( this, i0, l )
  subroutine TComponent_Unit2Atom( this, l )

    implicit none

    ! Include MPI header
#if MPI_VER > 0 && !defined(MPI_USE_MODULE)
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
    integer                        :: i, j, iUnit

!     ! i0 startpos of proc, i1 endpos of proc; l = difference
!     i1 = l + i0 - 1

    ! Broadcast positions and orientations to all processes
#if MPI_VER > 0
    ! in MC simulations, we only communicate during common equilibration
    if ( SimulationType .ne. MonteCarlo .or. (Equilibration .and. CommonEqui) .or. (mpiMCCommonGroups > 0)) then
      call MPI_Bcast( this%P0(:, :, :), size( this%P0 ), MPI_RK, NRootProc, Communicator, ierror )
      if( this%Molecule%isElongated ) then
        call MPI_Bcast( this%Q0(:, :, :), size( this%Q0 ), MPI_RK, NRootProc, Communicator, ierror )
      end if
    end if
#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    ! Loop over all units in Molecule
    do iUnit = 1, this%Molecule%NUnit
      ! Check number of rotation axes
      if( this%Molecule%Unit(iUnit)%isElongated ) then
        ! Loop over molecules
        do i = 1, l
          PX(i) = this%P0(i, 1, iUnit)
          PY(i) = this%P0(i, 2, iUnit)
          PZ(i) = this%P0(i, 3, iUnit)
          q1 = this%Q0(i, 1, iUnit)
          q2 = this%Q0(i, 2, iUnit)
          q3 = this%Q0(i, 3, iUnit)
          q4 = this%Q0(i, 4, iUnit)

          ! Normalise quaternions
          qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
          q1 = q1 * qinv
          q2 = q2 * qinv
          q3 = q3 * qinv
          q4 = q4 * qinv
          this%Q0(i, 1, iUnit) = q1
          this%Q0(i, 2, iUnit) = q2
          this%Q0(i, 3, iUnit) = q3
          this%Q0(i, 4, iUnit) = q4

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

        ! Loop over MIEnm sites in unit
        do j = 1, this%Molecule%Unit(iUnit)%NMIEnm
          pMIEnm => this%Molecule%Unit(iUnit)%SiteMIEnm(j)
          r1 = pMIEnm%r(1) * BoxLengthInv
          r2 = pMIEnm%r(2) * BoxLengthInv
          r3 = pMIEnm%r(3) * BoxLengthInv
          do i = 1, l
            pMIEnm%RX(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
            pMIEnm%RY(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
            pMIEnm%RZ(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
          end do
        end do

      ! Loop over TT68 sites in molecule
      do j = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(j)
        r1 = pTT68%r(1) * BoxLengthInv
        r2 = pTT68%r(2) * BoxLengthInv
        r3 = pTT68%r(3) * BoxLengthInv
        do i = 1, l
          pTT68%RX(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pTT68%RY(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pTT68%RZ(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
        end do
      end do

        ! Loop over charge sites in molecule
        do j = 1, this%Molecule%Unit(iUnit)%NCharge
          pCharge => this%Molecule%Unit(iUnit)%SiteCharge(j)
          r1 = pCharge%r(1) * BoxLengthInv
          r2 = pCharge%r(2) * BoxLengthInv
          r3 = pCharge%r(3) * BoxLengthInv
          do i = 1, l
            pCharge%RX(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
            pCharge%RY(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
            pCharge%RZ(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
          end do
        end do

        ! Loop over dipole sites in molecule
        do j = 1, this%Molecule%Unit(iUnit)%NDipole
          pDipole => this%Molecule%Unit(iUnit)%SiteDipole(j)
          r1 = pDipole%r(1) * BoxLengthInv
          r2 = pDipole%r(2) * BoxLengthInv
          r3 = pDipole%r(3) * BoxLengthInv
          or1 = pDipole%or(1)
          or2 = pDipole%or(2)
          or3 = pDipole%or(3)
          do i = 1, l
            pDipole%RX(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
            pDipole%RY(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
            pDipole%RZ(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
            pDipole%OX(i) = or1 * A11(i) + or2 * A21(i) + or3 * A31(i)
            pDipole%OY(i) = or1 * A12(i) + or2 * A22(i) + or3 * A32(i)
            pDipole%OZ(i) = or1 * A13(i) + or2 * A23(i) + or3 * A33(i)
          end do
        end do

        ! Loop over quadrupole sites in molecule
        do j = 1, this%Molecule%Unit(iUnit)%NQuadrupole
          pQuadrupole => this%Molecule%Unit(iUnit)%SiteQuadrupole(j)
          r1 = pQuadrupole%r(1) * BoxLengthInv
          r2 = pQuadrupole%r(2) * BoxLengthInv
          r3 = pQuadrupole%r(3) * BoxLengthInv
          or1 = pQuadrupole%or(1)
          or2 = pQuadrupole%or(2)
          or3 = pQuadrupole%or(3)
          do i = 1, l
            pQuadrupole%RX(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
            pQuadrupole%RY(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
            pQuadrupole%RZ(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
            pQuadrupole%OX(i) = or1 * A11(i) + or2 * A21(i) + or3 * A31(i)
            pQuadrupole%OY(i) = or1 * A12(i) + or2 * A22(i) + or3 * A32(i)
            pQuadrupole%OZ(i) = or1 * A13(i) + or2 * A23(i) + or3 * A33(i)
          end do
        end do

      ! Rotate total dipole moment
        if( CutoffMode .eq. CenterofMass ) then
          if (.not. UseIntDegFreed .and. ( SimulationType .ne. MonteCarlo )) then
              mue1 = this%Molecule%Mue(1)
              mue2 = this%Molecule%Mue(2)
              mue3 = this%Molecule%Mue(3)
          else
              mue1 = this%Molecule%Unit(iUnit)%Mue(1)
              mue2 = this%Molecule%Unit(iUnit)%Mue(2)
              mue3 = this%Molecule%Unit(iUnit)%Mue(3)
          end if
          do i = 1, l
            this%MueX(i, iUnit) = mue1 * A11(i) + mue2 * A21(i) + mue3 * A31(i)
            this%MueY(i, iUnit) = mue1 * A12(i) + mue2 * A22(i) + mue3 * A32(i)
            this%MueZ(i, iUnit) = mue1 * A13(i) + mue2 * A23(i) + mue3 * A33(i)
          end do
        end if

      else

        ! Loop over MIEnm sites in molecule
        do i = 1, this%Molecule%Unit(iUnit)%NMIEnm
          pMIEnm => this%Molecule%Unit(iUnit)%SiteMIEnm(i)
          do j = 1, l
            pMIEnm%RX(j) = this%P0(j, 1, iUnit)
            pMIEnm%RY(j) = this%P0(j, 2, iUnit)
            pMIEnm%RZ(j) = this%P0(j, 3, iUnit)
          end do
        end do

      ! Loop over TT68 sites in molecule
      do j = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(j)
        do i = 1, l
          pTT68%RX(i) = this%P0(i, 1, iUnit)
          pTT68%RY(i) = this%P0(i, 2, iUnit)
          pTT68%RZ(i) = this%P0(i, 3, iUnit)
        end do
      end do

        ! Loop over charge sites in molecule
        if ((LongRange .ne. RField) .or. UseIntDegFreed) then
            do i = 1, this%Molecule%Unit(iUnit)%NCharge
              pCharge => this%Molecule%Unit(iUnit)%SiteCharge(i)
              do j = 1, l
                pCharge%RX(j) = this%P0(j, 1, iUnit)
                pCharge%RY(j) = this%P0(j, 2, iUnit)
                pCharge%RZ(j) = this%P0(j, 3, iUnit)
              end do
            end do
        end if
      end if
    end do

  end subroutine TComponent_Unit2Atom


!==============================================================!
!  Subroutine TComponent_Unit2Atom1 (per unit)                 !
!==============================================================!

  subroutine TComponent_Unit2Atom1( this, np, nu )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: np
    integer, intent(in) :: nu

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

    ! Broadcast positions and orientations to all processes
#if MPI_VER > 0
    ! in MC simulations, we only communicate during common equilibration
    if ( SimulationType .ne. MonteCarlo .or. ((Equilibration .and. CommonEqui) )) then
      call MPI_Bcast( this%P0(np, :, nu), 3, MPI_RK, NRootProc, Communicator, ierror )
      if( this%Molecule%isElongated ) then
        call MPI_Bcast( this%Q0(np, :, nu), 4, MPI_RK, NRootProc, Communicator, ierror )
      end if
    end if
#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    ! Positions and quaternions of unit k in particle i
    PXi = this%P0(np, 1, nu)
    PYi = this%P0(np, 2, nu)
    PZi = this%P0(np, 3, nu)

    ! Check number of rotation axes
    if( this%Molecule%Unit(nu)%isElongated ) then

      ! Normalise quaternions
      q1 = this%Q0(np, 1, nu)
      q2 = this%Q0(np, 2, nu)
      q3 = this%Q0(np, 3, nu)
      q4 = this%Q0(np, 4, nu)
      qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
      q1 = q1 * qinv
      q2 = q2 * qinv
      q3 = q3 * qinv
      q4 = q4 * qinv
      this%Q0(np, 1, nu) = q1
      this%Q0(np, 2, nu) = q2
      this%Q0(np, 3, nu) = q3
      this%Q0(np, 4, nu) = q4

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

      ! Loop over MIEnm sites in unit
      do i = 1, this%Molecule%Unit(nu)%NMIEnm
        pMIEnm => this%Molecule%Unit(nu)%SiteMIEnm(i)
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
      do i = 1, this%Molecule%Unit(nu)%NCharge
        pCharge => this%Molecule%Unit(nu)%SiteCharge(i)
        r1 = pCharge%r(1) * BoxLengthInv
        r2 = pCharge%r(2) * BoxLengthInv
        r3 = pCharge%r(3) * BoxLengthInv
        pCharge%RX(np) = PXi + r1 * A11 + r2 * A21 + r3 * A31
        pCharge%RY(np) = PYi + r1 * A12 + r2 * A22 + r3 * A32
        pCharge%RZ(np) = PZi + r1 * A13 + r2 * A23 + r3 * A33
      end do

      ! Loop over dipole sites in molecule
      do i = 1, this%Molecule%Unit(nu)%NDipole
        pDipole => this%Molecule%Unit(nu)%SiteDipole(i)
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
      do i = 1, this%Molecule%Unit(nu)%NQuadrupole
        pQuadrupole => this%Molecule%Unit(nu)%SiteQuadrupole(i)
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
        mue1 = this%Molecule%Unit(nu)%Mue(1)
        mue2 = this%Molecule%Unit(nu)%Mue(2)
        mue3 = this%Molecule%Unit(nu)%Mue(3)
        this%MueX(np, nu) = mue1 * A11 + mue2 * A21 + mue3 * A31
        this%MueY(np, nu) = mue1 * A12 + mue2 * A22 + mue3 * A32
        this%MueZ(np, nu) = mue1 * A13 + mue2 * A23 + mue3 * A33
      end if

    else

      ! Loop over MIEnm sites in molecule
      do i = 1, this%Molecule%Unit(nu)%NMIEnm
        pMIEnm => this%Molecule%Unit(nu)%SiteMIEnm(i)
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
      if ((LongRange .ne. RField) .or. UseIntDegFreed) then
        do i = 1, this%Molecule%Unit(nu)%NCharge
          pCharge => this%Molecule%Unit(nu)%SiteCharge(i)
          pCharge%RX(np) = PXi
          pCharge%RY(np) = PYi
          pCharge%RZ(np) = PZi
        end do
      end if

    end if

  end subroutine TComponent_Unit2Atom1


!==============================================================!
!  Subroutine TComponent_Unit2AtomTest                         !
!==============================================================!

  subroutine TComponent_Unit2AtomTest( this, np, nu )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: np
    integer, intent(in) :: nu

    ! Declare local variables
    real(RK)                       :: BoxLengthInv
    real(RK)                       :: PX(np), PY(np), PZ(np)
    real(RK)                       :: q1, q2, q3, q4, qinv
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
    integer                        :: i, i0, i1, j, k

#if MPI_VER > 0
    if (UseIntDegFreed) then
      i0 = this%NTest0
      i1 = this%NTest2
    else
      i0 = 1
      i1 = this%NTest
    end if
#else
    i0 = 1
    i1 = this%NTest
#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    ! Loop over all units in Molecule
    do k = 1, nu
      ! Check number of rotation axes
      if( this%Molecule%Unit(k)%isElongated ) then

        ! Loop over molecules
        do i = i0, i1
          ! Positions and quaternions of test particle i
          PX(i) = this%P0Test(i, 1, k)
          PY(i) = this%P0Test(i, 2, k)
          PZ(i) = this%P0Test(i, 3, k)
          q1 = this%Q0Test(i, 1, k)
          q2 = this%Q0Test(i, 2, k)
          q3 = this%Q0Test(i, 3, k)
          q4 = this%Q0Test(i, 4, k)

          ! Normalise quaternions
#if ARCH == 3
!          qinv = rsqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#else
!          qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#endif
!          q1 = q1 * qinv
!          q2 = q2 * qinv
!          q3 = q3 * qinv
!          q4 = q4 * qinv
!          this%Q0Test(i, 1, k) = q1
!          this%Q0Test(i, 2, k) = q2
!          this%Q0Test(i, 3, k) = q3
!          this%Q0Test(i, 4, k) = q4

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

        ! Loop over MIEnm sites in unit
        do j = 1, this%Molecule%Unit(k)%NMIEnm
          pMIEnm => this%Molecule%Unit(k)%SiteMIEnm(j)
          r1 = pMIEnm%r(1) * BoxLengthInv
          r2 = pMIEnm%r(2) * BoxLengthInv
          r3 = pMIEnm%r(3) * BoxLengthInv
          do i = i0, i1
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
        do j = 1, this%Molecule%Unit(k)%NCharge
          pCharge => this%Molecule%Unit(k)%SiteCharge(j)
          r1 = pCharge%r(1) * BoxLengthInv
          r2 = pCharge%r(2) * BoxLengthInv
          r3 = pCharge%r(3) * BoxLengthInv
          do i = i0, i1
            pCharge%RXTest(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
            pCharge%RYTest(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
            pCharge%RZTest(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
          end do
        end do

        ! Loop over dipole sites in molecule
        do j = 1, this%Molecule%Unit(k)%NDipole
          pDipole => this%Molecule%Unit(k)%SiteDipole(j)
          r1 = pDipole%r(1) * BoxLengthInv
          r2 = pDipole%r(2) * BoxLengthInv
          r3 = pDipole%r(3) * BoxLengthInv
          or1 = pDipole%or(1)
          or2 = pDipole%or(2)
          or3 = pDipole%or(3)
          do i = i0, i1
            pDipole%RXTest(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
            pDipole%RYTest(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
            pDipole%RZTest(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
            pDipole%OXTest(i) = or1 * A11(i) + or2 * A21(i) + or3 * A31(i)
            pDipole%OYTest(i) = or1 * A12(i) + or2 * A22(i) + or3 * A32(i)
            pDipole%OZTest(i) = or1 * A13(i) + or2 * A23(i) + or3 * A33(i)
          end do
        end do

        ! Loop over quadrupole sites in molecule
        do j = 1, this%Molecule%Unit(k)%NQuadrupole
          pQuadrupole => this%Molecule%Unit(k)%SiteQuadrupole(j)
          r1 = pQuadrupole%r(1) * BoxLengthInv
          r2 = pQuadrupole%r(2) * BoxLengthInv
          r3 = pQuadrupole%r(3) * BoxLengthInv
          or1 = pQuadrupole%or(1)
          or2 = pQuadrupole%or(2)
          or3 = pQuadrupole%or(3)
          do i = i0, i1
            pQuadrupole%RXTest(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
            pQuadrupole%RYTest(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
            pQuadrupole%RZTest(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)

            pQuadrupole%OXTest(i) = or1 * A11(i) + or2 * A21(i) + or3 * A31(i)
            pQuadrupole%OYTest(i) = or1 * A12(i) + or2 * A22(i) + or3 * A32(i)
            pQuadrupole%OZTest(i) = or1 * A13(i) + or2 * A23(i) + or3 * A33(i)
          end do
        end do

        if( CutoffMode .eq. CenterofMass ) then
          mue1 = this%Molecule%Unit(k)%Mue(1)
          mue2 = this%Molecule%Unit(k)%Mue(2)
          mue3 = this%Molecule%Unit(k)%Mue(3)
          do i = i0, i1
            this%MueXTest(i, k) = mue1 * A11(i) + mue2 * A21(i) + mue3 * A31(i)
            this%MueYTest(i, k) = mue1 * A12(i) + mue2 * A22(i) + mue3 * A32(i)
            this%MueZTest(i, k) = mue1 * A13(i) + mue2 * A23(i) + mue3 * A33(i)
          end do
        end if

      else

        ! Loop over MIEnm sites in molecule
        do i = 1, this%Molecule%Unit(k)%NMIEnm
          pMIEnm => this%Molecule%Unit(k)%SiteMIEnm(i)
          do j = i0, i1
            pMIEnm%RXTest(j) = this%P0Test(j, 1, k)
            pMIEnm%RYTest(j) = this%P0Test(j, 2, k)
            pMIEnm%RZTest(j) = this%P0Test(j, 3, k)
          end do
        end do

      ! Loop over TT68 sites in molecule
      do i = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(i)
        do j = 1, np
          pTT68%RXTest(j) = this%P0Test(j, 1, k)
          pTT68%RYTest(j) = this%P0Test(j, 2, k)
          pTT68%RZTest(j) = this%P0Test(j, 3, k)
        end do
      end do

        ! Loop over charge sites in molecule
        if ((LongRange .ne. RField) .or. UseIntDegFreed) then
          do i = 1, this%Molecule%Unit(k)%NCharge
            pCharge => this%Molecule%Unit(k)%SiteCharge(i)
            do j = i0, i1
              pCharge%RXTest(j) = this%P0Test(j, 1, k)
              pCharge%RYTest(j) = this%P0Test(j, 2, k)
              pCharge%RZTest(j) = this%P0Test(j, 3, k)
            end do
          end do
        end if

      end if
    end do

  end subroutine TComponent_Unit2AtomTest


!==============================================================!
!  Subroutine TComponent_Atom2Unit                             !
!==============================================================!

!   subroutine TComponent_Atom2Mol( this, i0, l )
  subroutine TComponent_Atom2Unit( this, l )

    implicit none

    ! Include MPI header
#if MPI_VER > 0 && !defined(MPI_USE_MODULE)
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
    integer                        :: i0, i1, i, j, iUnit
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
    this%F(:, :, :) = 0._RK
#if OSMOP > 0
    !if (RootProc) then
      this%FOsmoticPressure(:) = 0._RK
      if ( .not. this%permeable ) then
        do i = i0, i1 ! 1, this%NPart
          do iUnit = 1, this%Molecule%NUnit
            if( this%P0(i, iUnit,1) .ge. 0.25_RK ) then
              this%FOsmoticPressure(i) = kForceOsmoticPressure*( this%P0(i, iUnit,1)-.25_RK)*BoxLength
              this%F(i, iUnit,1) = this%F(i, iUnit,1) - this%FOsmoticPressure(i)
            elseif( this%P0(i,iUnit,1) .le. -0.25_RK ) then
              this%FOsmoticPressure(i) = kForceOsmoticPressure*(-this%P0(i, iUnit,1)-.25_RK)*BoxLength
              this%F(i, iUnit,1) = this%F(i, iUnit,1) + this%FOsmoticPressure(i)
            end if
          end do
        end do
      end if
    !end if
#if MPI_VER > 0
    call MPI_Reduce( this%FOsmoticPressure(:), OsmoPAll(:), size( this%FOsmoticPressure), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    if (RootProc) this%FOsmoticPressure(:) = OsmoPAll(:)
#endif
#endif

    ! Loop over all Units in Molecule
    do iUnit = 1, this%Molecule%NUnit
      ! Check number of rotation axes
      if( this%Molecule%Unit(iUnit)%isElongated ) then

        ! Initialize torques
        this%T(:, :, iUnit) = 0._RK

        ! Initialize local arrays
        rx(:) = this%P0(:, 1, iUnit)
        ry(:) = this%P0(:, 2, iUnit)
        rz(:) = this%P0(:, 3, iUnit)
        q1(:) = this%Q0(:, 1, iUnit)
        q2(:) = this%Q0(:, 2, iUnit)
        q3(:) = this%Q0(:, 3, iUnit)
        q4(:) = this%Q0(:, 4, iUnit)

        ! Loop over MIEnm sites in unit
        do j = 1, this%Molecule%Unit(iUnit)%NMIEnm
          pMIEnm => this%Molecule%Unit(iUnit)%SiteMIEnm(j)
          do i = 1, l
            fx = pMIEnm%FX(i-1+i0)
            fy = pMIEnm%FY(i-1+i0)
            fz = pMIEnm%FZ(i-1+i0)
            r1x = ( pMIEnm%RX(i-1+i0) - rx(i) ) * BoxLength
            r1y = ( pMIEnm%RY(i-1+i0) - ry(i) ) * BoxLength
            r1z = ( pMIEnm%RZ(i-1+i0) - rz(i) ) * BoxLength
            this%F(i-1+i0, 1, iUnit) = this%F(i-1+i0, 1, iUnit) + fx
            this%F(i-1+i0, 2, iUnit) = this%F(i-1+i0, 2, iUnit) + fy
            this%F(i-1+i0, 3, iUnit) = this%F(i-1+i0, 3, iUnit) + fz
            this%T(i-1+i0, 1, iUnit) = this%T(i-1+i0, 1, iUnit) + r1y * fz - r1z * fy
            this%T(i-1+i0, 2, iUnit) = this%T(i-1+i0, 2, iUnit) + r1z * fx - r1x * fz
            this%T(i-1+i0, 3, iUnit) = this%T(i-1+i0, 3, iUnit) + r1x * fy - r1y * fx
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
          this%F(i-1+i0, 1, iUnit) = this%F(i-1+i0, 1, iUnit) + fx
          this%F(i-1+i0, 2, iUnit) = this%F(i-1+i0, 2, iUnit) + fy
          this%F(i-1+i0, 3, iUnit) = this%F(i-1+i0, 3, iUnit) + fz
          this%T(i-1+i0, 1, iUnit) = this%T(i-1+i0, 1, iUnit) + r1y * fz - r1z * fy
          this%T(i-1+i0, 2, iUnit) = this%T(i-1+i0, 2, iUnit) + r1z * fx - r1x * fz
          this%T(i-1+i0, 3, iUnit) = this%T(i-1+i0, 3, iUnit) + r1x * fy - r1y * fx
        end do
      end do

        ! Loop over charge sites in unit
        do j = 1, this%Molecule%Unit(iUnit)%NCharge
          pCharge => this%Molecule%Unit(iUnit)%SiteCharge(j)
          do i = 1, l
            fx = pCharge%FX(i-1+i0)
            fy = pCharge%FY(i-1+i0)
            fz = pCharge%FZ(i-1+i0)
            r1x = ( pCharge%RX(i-1+i0) - rx(i) ) * BoxLength
            r1y = ( pCharge%RY(i-1+i0) - ry(i) ) * BoxLength
            r1z = ( pCharge%RZ(i-1+i0) - rz(i) ) * BoxLength
            this%F(i-1+i0, 1, iUnit) = this%F(i-1+i0, 1, iUnit) + fx
            this%F(i-1+i0, 2, iUnit) = this%F(i-1+i0, 2, iUnit) + fy
            this%F(i-1+i0, 3, iUnit) = this%F(i-1+i0, 3, iUnit) + fz
            this%T(i-1+i0, 1, iUnit) = this%T(i-1+i0, 1, iUnit) + r1y * fz - r1z * fy
            this%T(i-1+i0, 2, iUnit) = this%T(i-1+i0, 2, iUnit) + r1z * fx - r1x * fz
            this%T(i-1+i0, 3, iUnit) = this%T(i-1+i0, 3, iUnit) + r1x * fy - r1y * fx
          end do
        end do

        ! Loop over dipole sites in unit
        do j = 1, this%Molecule%Unit(iUnit)%NDipole
          pDipole => this%Molecule%Unit(iUnit)%SiteDipole(j)
          do i = 1, l
            fx = pDipole%FX(i-1+i0)
            fy = pDipole%FY(i-1+i0)
            fz = pDipole%FZ(i-1+i0)
            r1x = ( pDipole%RX(i-1+i0) - rx(i) ) * BoxLength
            r1y = ( pDipole%RY(i-1+i0) - ry(i) ) * BoxLength
            r1z = ( pDipole%RZ(i-1+i0) - rz(i) ) * BoxLength
            this%F(i-1+i0, 1, iUnit) = this%F(i-1+i0, 1, iUnit) + fx
            this%F(i-1+i0, 2, iUnit) = this%F(i-1+i0, 2, iUnit) + fy
            this%F(i-1+i0, 3, iUnit) = this%F(i-1+i0, 3, iUnit) + fz
            this%T(i-1+i0, 1, iUnit) = this%T(i-1+i0, 1, iUnit) + pDipole%OY(i) * pDipole%TZ(i) &
&                                       - pDipole%OZ(i) * pDipole%TY(i) + r1y * fz - r1z * fy
            this%T(i-1+i0, 2, iUnit) = this%T(i-1+i0, 2, iUnit) + pDipole%OZ(i) * pDipole%TX(i) &
&                                       - pDipole%OX(i) * pDipole%TZ(i) + r1z * fx - r1x * fz
            this%T(i-1+i0, 3, iUnit) = this%T(i-1+i0, 3, iUnit) + pDipole%OX(i) * pDipole%TY(i) &
&                                       - pDipole%OY(i) * pDipole%TX(i) + r1x * fy - r1y * fx
          end do
        end do

        ! Loop over quadrupole sites in unit
        do j = 1, this%Molecule%Unit(iUnit)%NQuadrupole
          pQuadrupole => this%Molecule%Unit(iUnit)%SiteQuadrupole(j)
          do i = 1, l
            fx = pQuadrupole%FX(i-1+i0)
            fy = pQuadrupole%FY(i-1+i0)
            fz = pQuadrupole%FZ(i-1+i0)
            r1x = ( pQuadrupole%RX(i-1+i0) - rx(i) ) * BoxLength
            r1y = ( pQuadrupole%RY(i-1+i0) - ry(i) ) * BoxLength
            r1z = ( pQuadrupole%RZ(i-1+i0) - rz(i) ) * BoxLength
            this%F(i-1+i0, 1, iUnit) = this%F(i-1+i0, 1, iUnit) + fx
            this%F(i-1+i0, 2, iUnit) = this%F(i-1+i0, 2, iUnit) + fy
            this%F(i-1+i0, 3, iUnit) = this%F(i-1+i0, 3, iUnit) + fz
            this%T(i-1+i0, 1, iUnit) = this%T(i-1+i0, 1, iUnit) + pQuadrupole%OY(i) * pQuadrupole%TZ(i) &
&                                       - pQuadrupole%OZ(i) * pQuadrupole%TY(i) + r1y * fz - r1z * fy
            this%T(i-1+i0, 2, iUnit) = this%T(i-1+i0, 2, iUnit) + pQuadrupole%OZ(i) * pQuadrupole%TX(i) &
&                                       - pQuadrupole%OX(i) * pQuadrupole%TZ(i)+ r1z * fx - r1x * fz
            this%T(i-1+i0, 3, iUnit) = this%T(i-1+i0, 3, iUnit) + pQuadrupole%OX(i) * pQuadrupole%TY(i) &
&                                       - pQuadrupole%OY(i) * pQuadrupole%TX(i) + r1x * fy - r1y * fx
          end do
        end do

        do i = 1, l
          ! Add torques from reaction field
          tx = this%T(i-1+i0, 1, iUnit) + this%tRFX(i-1+i0, iUnit)
          ty = this%T(i-1+i0, 2, iUnit) + this%tRFY(i-1+i0, iUnit)
          tz = this%T(i-1+i0, 3, iUnit) + this%tRFZ(i-1+i0, iUnit)

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
          this%T(i-1+i0, 1, iUnit) = A11 * tx + A12 * ty + A13 * tz
          this%T(i-1+i0, 2, iUnit) = A21 * tx + A22 * ty + A23 * tz
          this%T(i-1+i0, 3, iUnit) = A31 * tx + A32 * ty + A33 * tz
        end do

      else

        ! Loop over MIEnm sites in unit
        do j = 1, this%Molecule%Unit(iUnit)%NMIEnm
          pMIEnm => this%Molecule%Unit(iUnit)%SiteMIEnm(j)
          do i = 1, l
            this%F(i-1+i0, 1, iUnit) = this%F(i-1+i0, 1, iUnit) + pMIEnm%FX(i-1+i0)
            this%F(i-1+i0, 2, iUnit) = this%F(i-1+i0, 2, iUnit) + pMIEnm%FY(i-1+i0)
            this%F(i-1+i0, 3, iUnit) = this%F(i-1+i0, 3, iUnit) + pMIEnm%FZ(i-1+i0)
          end do
        end do

      ! Loop over TT68 sites in molecule
      do j = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(j)
        do i = 1, l
          this%F(i-1+i0, 1, iUnit) = this%F(i-1+i0, 1, iUnit) + pTT68%FX(i-1+i0)
          this%F(i-1+i0, 2, iUnit) = this%F(i-1+i0, 2, iUnit) + pTT68%FY(i-1+i0)
          this%F(i-1+i0, 3, iUnit) = this%F(i-1+i0, 3, iUnit) + pTT68%FZ(i-1+i0)
        end do
      end do

        ! Loop over charge sites in molecule
        if ((LongRange .ne. RField) .or. UseIntDegFreed) then
          do j = 1, this%Molecule%Unit(iUnit)%NCharge
            pCharge => this%Molecule%Unit(iUnit)%SiteCharge(j)
            do i = 1, l
              this%F(i-1+i0, 1, iUnit) = this%F(i-1+i0, 1, iUnit) + pCharge%FX(i)
              this%F(i-1+i0, 2, iUnit) = this%F(i-1+i0, 2, iUnit) + pCharge%FY(i)
              this%F(i-1+i0, 3, iUnit) = this%F(i-1+i0, 3, iUnit) + pCharge%FZ(i)
            end do
          end do
        end if

      end if

    end do

#if MPI_VER > 0
    ! Reduce forces and torques from all processes
    call MPI_Reduce( this%F(:, :, :), this%FAll(:, :, :), size( this%F ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    if( this%Molecule%isElongated ) &
&     call MPI_Reduce( this%T(:, :, :), this%TAll(:, :, :), size( this%T ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
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

  end subroutine TComponent_Atom2Unit

#if 0
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
#endif


!==============================================================!
!  Subroutine TComponent_Atom2Unit_Trans                       !
!==============================================================!

  subroutine TComponent_Atom2Unit_Trans( this, np, nu )

    implicit none

    ! Include MPI header
#if MPI_VER > 0 && !defined(MPI_USE_MODULE)
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: np
    integer, intent(in) :: nu

    ! Declare local variables
    real(RK)                       :: BoxLength
    real(RK)                       :: rx(np,nu), ry(np,nu), rz(np,nu), r1x, r1y, r1z
    real(RK)                       :: q1(np,nu), q2(np,nu), q3(np,nu), q4(np,nu)
    real(RK)                       :: fx, fy, fz, tx, ty, tz
    real(RK)                       :: A11, A12, A13, A21, A22, A23, A31, A32, A33
    type(TSiteMIEnm), pointer      :: pMIEnm
    type(TSiteTT68), pointer       :: pTT68
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j, k

#if  TRANS == 1
    real(RK)                       :: vsx,vsy,vsz
    real(RK)                       :: vsux,vsuy,vsuz
    real(RK)                       :: vbx,vby,vbz
    real(RK)                       :: cx, cy, cz
    real(RK)                       :: tux, tuy, tuz, tlx, tly, tlz, tdx, tdy, tdz
    real(RK)                       :: BoxLength_dt

    BoxLength_dt = this%BoxLength/TimeStep

    this%FS(:,:) = 0._RK
    this%FB(:,:) = 0._RK
    this%FTC(:,:) = 0._RK
    this%FRC(:,:) = 0._RK
    this%FTC1(:,:) = 0._RK
    this%FTC2(:,:) = 0._RK
    this%FTC3(:,:) = 0._RK
    this%FRC1(:,:) = 0._RK
    this%FRC2(:,:) = 0._RK
    this%FRC3(:,:) = 0._RK
#endif
    
    ! Assign local variables
    BoxLength = this%BoxLength

    k = nu   ! Michael Sch.: add do loop over molecules when implementing ms2 with transportproperties for flexible molecules

    ! Initialize forces
    this%F(1:np, :, k) = 0._RK

    ! Check number of rotation axes
    if( this%Molecule%Unit(k)%isElongated ) then

      ! Initialize torques
      this%T(1:np, :, k) = 0._RK

      ! Initialize local arrays
      rx(:,k) = this%P0(:, 1, k)
      ry(:,k) = this%P0(:, 2, k)
      rz(:,k) = this%P0(:, 3, k)
      q1(:,k) = this%Q0(:, 1, k)
      q2(:,k) = this%Q0(:, 2, k)
      q3(:,k) = this%Q0(:, 3, k)
      q4(:,k) = this%Q0(:, 4, k)

      ! Loop over MIEnm sites in unit
      do j = 1, this%Molecule%NMIEnm
        pMIEnm => this%Molecule%SiteMIEnm(j)
        do i = 1, np
          fx = pMIEnm%FX(i)
          fy = pMIEnm%FY(i)
          fz = pMIEnm%FZ(i)
#if  TRANS == 1
      
          vsx = pMIEnm%vsMIEx(i)
          vsy = pMIEnm%vsMIEy(i)
          vsz = pMIEnm%vsMIEz(i)
          vbx = pMIEnm%vbMIEx(i)
          vby = pMIEnm%vbMIEy(i)
          vbz = pMIEnm%vbMIEz(i)
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
#endif
          r1x = ( pMIEnm%RX(i) - rx(i,k) ) * BoxLength
          r1y = ( pMIEnm%RY(i) - ry(i,k) ) * BoxLength
          r1z = ( pMIEnm%RZ(i) - rz(i,k) ) * BoxLength
          this%F(i, 1, k) = this%F(i, 1, k) + fx
          this%F(i, 2, k) = this%F(i, 2, k) + fy
          this%F(i, 3, k) = this%F(i, 3, k) + fz
          this%T(i, 1, k) = this%T(i, 1, k) + r1y * fz - r1z * fy
          this%T(i, 2, k) = this%T(i, 2, k) + r1z * fx - r1x * fz
          this%T(i, 3, k) = this%T(i, 3, k) + r1x * fy - r1y * fx
#if  TRANS == 1

          this%FS(i, 1)= this%FS(i, 1)+ vsx
          this%FS(i, 2)= this%FS(i, 2)+ vsy
          this%FS(i, 3)= this%FS(i, 3)+ vsz
          this%FB(i, 1)= this%FB(i, 1)+ vbx
          this%FB(i, 2)= this%FB(i, 2)+ vby
          this%FB(i, 3)= this%FB(i, 3)+ vbz

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
  
          vsx = pTT68%vsTTx(i)
          vsy = pTT68%vsTTy(i)
          vsz = pTT68%vsTTz(i)
          vbx = pTT68%vbTTx(i)
          vby = pTT68%vbTTy(i)
          vbz = pTT68%vbTTz(i)
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
#endif
          r1x = ( pTT68%RX(i) - rx(i, 1) ) * BoxLength
          r1y = ( pTT68%RY(i) - ry(i, 1) ) * BoxLength
          r1z = ( pTT68%RZ(i) - rz(i, 1) ) * BoxLength
          this%F(i, 1, k) = this%F(i, 1, k) +  fx
          this%F(i, 2, k) = this%F(i, 2, k) + fy
          this%F(i, 3, k) = this%F(i, 3, k) + fz
          this%T(i, 1, k) = this%T(i, 1, k) + r1y * fz - r1z * fy
          this%T(i, 2, k) = this%T(i, 2, k) + r1z * fx - r1x * fz
          this%T(i, 3, k) = this%T(i, 3, k) + r1x * fy - r1y * fx
#if  TRANS == 1
    
          this%FS(i, 1)= this%FS(i, 1)+ vsx
          this%FS(i, 2)= this%FS(i, 2)+ vsy
          this%FS(i, 3)= this%FS(i, 3)+ vsz
          this%FB(i, 1)= this%FB(i, 1)+ vbx
          this%FB(i, 2)= this%FB(i, 2)+ vby
          this%FB(i, 3)= this%FB(i, 3)+ vbz

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
#endif
        end do
      end do

      ! Loop over charge sites in unit
      do j = 1, this%Molecule%NCharge
        pCharge => this%Molecule%SiteCharge(j)
        do i = 1, np
          fx = pCharge%FX(i)
          fy = pCharge%FY(i)
          fz = pCharge%FZ(i)
#if  TRANS == 1

          vsx = pCharge%vsCx(i)
          vsy = pCharge%vsCy(i)
          vsz = pCharge%vsCz(i)
          vbx = pCharge%vbCx(i)
          vby = pCharge%vbCy(i)
          vbz = pCharge%vbCz(i)
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
         
#endif
          r1x = ( pCharge%RX(i) - rx(i,k) ) * BoxLength
          r1y = ( pCharge%RY(i) - ry(i,k) ) * BoxLength
          r1z = ( pCharge%RZ(i) - rz(i,k) ) * BoxLength
          this%F(i, 1, k) = this%F(i, 1, k) + fx
          this%F(i, 2, k) = this%F(i, 2, k) + fy
          this%F(i, 3, k) = this%F(i, 3, k) + fz
          this%T(i, 1, k) = this%T(i, 1, k) + r1y * fz - r1z * fy
          this%T(i, 2, k) = this%T(i, 2, k) + r1z * fx - r1x * fz
          this%T(i, 3, k) = this%T(i, 3, k) + r1x * fy - r1y * fx
#if  TRANS == 1
     
          this%FS(i, 1)= this%FS(i, 1)+ vsx
          this%FS(i, 2)= this%FS(i, 2)+ vsy
          this%FS(i, 3)= this%FS(i, 3)+ vsz
          this%FB(i, 1)= this%FB(i, 1)+ vbx
          this%FB(i, 2)= this%FB(i, 2)+ vby
          this%FB(i, 3)= this%FB(i, 3)+ vbz

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
#endif
        end do
      end do

      ! Loop over dipole sites in unit
      do j = 1, this%Molecule%NDipole
        pDipole => this%Molecule%SiteDipole(j)
        do i = 1, np
          fx = pDipole%FX(i)
          fy = pDipole%FY(i)
          fz = pDipole%FZ(i)
#if  TRANS == 1
        
          vsx = pDipole%vsDx(i)
          vsy = pDipole%vsDy(i)
          vsz = pDipole%vsDz(i)
          vbx = pDipole%vbDx(i)
          vby = pDipole%vbDy(i)
          vbz = pDipole%vbDz(i)
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
         
#endif
          r1x = ( pDipole%RX(i) - rx(i,k) ) * BoxLength
          r1y = ( pDipole%RY(i) - ry(i,k) ) * BoxLength
          r1z = ( pDipole%RZ(i) - rz(i,k) ) * BoxLength
          this%F(i, 1, k) = this%F(i, 1, k) + fx
          this%F(i, 2, k) = this%F(i, 2, k) + fy
          this%F(i, 3, k) = this%F(i, 3, k) + fz
          this%T(i, 1, k) = this%T(i, 1, k) + pDipole%OY(i) * pDipole%TZ(i) &
&                                     - pDipole%OZ(i) * pDipole%TY(i) &
&                                     + r1y * fz - r1z * fy
          this%T(i, 2, k) = this%T(i, 2, k) + pDipole%OZ(i) * pDipole%TX(i) &
&                                     - pDipole%OX(i) * pDipole%TZ(i) &
&                                     + r1z * fx - r1x * fz
          this%T(i, 3, k) = this%T(i, 3, k) + pDipole%OX(i) * pDipole%TY(i) &
&                                     - pDipole%OY(i) * pDipole%TX(i) &
&                                     + r1x * fy - r1y * fx
#if  TRANS == 1
         
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
#endif
        end do
      end do

      ! Loop over quadrupole sites in unit
      do j = 1, this%Molecule%NQuadrupole
        pQuadrupole => this%Molecule%SiteQuadrupole(j)
        do i = 1, np
          fx = pQuadrupole%FX(i)
          fy = pQuadrupole%FY(i)
          fz = pQuadrupole%FZ(i)
#if  TRANS == 1
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
#endif
          r1x = ( pQuadrupole%RX(i) - rx(i,k) ) * BoxLength
          r1y = ( pQuadrupole%RY(i) - ry(i,k) ) * BoxLength
          r1z = ( pQuadrupole%RZ(i) - rz(i,k) ) * BoxLength
          this%F(i, 1, k) = this%F(i, 1, k) + fx
          this%F(i, 2, k) = this%F(i, 2, k) + fy
          this%F(i, 3, k) = this%F(i, 3, k) + fz
          this%T(i, 1, k) = this%T(i, 1, k) + pQuadrupole%OY(i) * pQuadrupole%TZ(i) &
&                                     - pQuadrupole%OZ(i) * pQuadrupole%TY(i) &
&                                     + r1y * fz - r1z * fy
          this%T(i, 2, k) = this%T(i, 2, k) + pQuadrupole%OZ(i) * pQuadrupole%TX(i) &
&                                     - pQuadrupole%OX(i) * pQuadrupole%TZ(i) &
&                                     + r1z * fx - r1x * fz
          this%T(i, 3, k) = this%T(i, 3, k) + pQuadrupole%OX(i) * pQuadrupole%TY(i) &
&                                     - pQuadrupole%OY(i) * pQuadrupole%TX(i) &
&                                     + r1x * fy - r1y * fx
#if  TRANS == 1
  
          this%FS(i, 1)= this%FS(i, 1)+ vsx
          this%FS(i, 2)= this%FS(i, 2)+ vsy
          this%FS(i, 3)= this%FS(i, 3)+ vsz
          this%FB(i, 1)= this%FB(i, 1)+ vbx
          this%FB(i, 2)= this%FB(i, 2)+ vby
          this%FB(i, 3)= this%FB(i, 3)+ vbz

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
#endif
        end do
      end do

      do i = 1, np
        ! Add torques from reaction field
        tx = this%T(i, 1, k) + this%tRFX(i, k)
        ty = this%T(i, 2, k) + this%tRFY(i, k)
        tz = this%T(i, 3, k) + this%tRFZ(i, k)

        ! Convert torque to body-fixed coordinates
        A11 = q1(i,k)**2 + q2(i,k)**2 - q3(i,k)**2 - q4(i,k)**2
        A12 = 2._RK * (q2(i,k) * q3(i,k) + q1(i,k) * q4(i,k))
        A13 = 2._RK * (q2(i,k) * q4(i,k) - q1(i,k) * q3(i,k))
        A21 = 2._RK * (q2(i,k) * q3(i,k) - q1(i,k) * q4(i,k))
        A22 = q1(i,k)**2 - q2(i,k)**2 + q3(i,k)**2 - q4(i,k)**2
        A23 = 2._RK * (q3(i,k) * q4(i,k) + q1(i,k) * q2(i,k))
        A31 = 2._RK * (q2(i,k) * q4(i,k) + q1(i,k) * q3(i,k))
        A32 = 2._RK * (q3(i,k) * q4(i,k) - q1(i,k) * q2(i,k))
        A33 = q1(i,k)**2 - q2(i,k)**2 - q3(i,k)**2 + q4(i,k)**2
        this%T(i, 1, k) = A11 * tx + A12 * ty + A13 * tz
        this%T(i, 2, k) = A21 * tx + A22 * ty + A23 * tz
        this%T(i, 3, k) = A31 * tx + A32 * ty + A33 * tz
      end do

    else

      ! Loop over MIEnm sites in unit
      do j = 1, this%Molecule%NMIEnm
        pMIEnm => this%Molecule%SiteMIEnm(j)
        do i = 1, np
#if  TRANS == 1
  
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
#endif
          this%F(i, 1, k) = this%F(i, 1, k) + pMIEnm%FX(i)
          this%F(i, 2, k) = this%F(i, 2, k) + pMIEnm%FY(i)
          this%F(i, 3, k) = this%F(i, 3, k) + pMIEnm%FZ(i)
#if  TRANS == 1

          this%FS(i, 1) = this%FS(i, 1) + vsx
          this%FS(i, 2) = this%FS(i, 2) + vsy
          this%FS(i, 3) = this%FS(i, 3) + vsz
          this%FB(i, 1) = this%FB(i, 1) + vbx
          this%FB(i, 2) = this%FB(i, 2) + vby
          this%FB(i, 3) = this%FB(i, 3) + vbz

          this%FTC1(i, 1)= this%FTC1(i, 1) +(cx+vbx)
          this%FTC1(i, 2)= this%FTC1(i, 2) + vsux
          this%FTC1(i, 3)= this%FTC1(i, 3) + vsuy
          this%FTC2(i, 1)= this%FTC2(i, 1) + vsx
          this%FTC2(i, 2)= this%FTC2(i, 2) +(cy+vby)
          this%FTC2(i, 3)= this%FTC2(i, 3) + vsuz
          this%FTC3(i, 1)= this%FTC3(i, 1) + vsy
          this%FTC3(i, 2)= this%FTC3(i, 2) + vsz
          this%FTC3(i, 3)= this%FTC3(i, 3) +(cz+vbz)
#endif
        end do
      end do

      ! Loop over TT68 sites in molecule
      do j = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(j)
        do i = 1, np
#if  TRANS == 1
          vsx = pTT68%vsTTx(i)
          vsy = pTT68%vsTTy(i)
          vsz = pTT68%vsTTz(i)
          vbx = pTT68%vbTTx(i)
          vby = pTT68%vbTTy(i)
          vbz = pTT68%vbTTz(i)
          vsux= pTT68%vsuTTx(i)
          vsuy= pTT68%vsuTTy(i)
          vsuz= pTT68%vsuTTz(i)
          cx  = pTT68%cTTx(i)
          cy  = pTT68%cTTy(i)
          cz  = pTT68%cTTz(i)
#endif
          this%F(i, 1, k) = this%F(i, 1, k) + pTT68%FX(i)
          this%F(i, 2, k) = this%F(i, 2, k) + pTT68%FY(i)
          this%F(i, 3, k) = this%F(i, 3, k) + pTT68%FZ(i)

#if  TRANS == 1

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
        
#endif
        end do
      end do

      ! Loop over TT68 sites in molecule
      do j = 1, this%Molecule%NTT68
        pTT68 => this%Molecule%SiteTT68(j)
        do i = 1, np
#if  TRANS == 1
          vsx = pTT68%vsTTx(i)
          vsy = pTT68%vsTTy(i)
          vsz = pTT68%vsTTz(i)
          vbx = pTT68%vbTTx(i)
          vby = pTT68%vbTTy(i)
          vbz = pTT68%vbTTz(i)
          vsux= pTT68%vsuTTx(i)
          vsuy= pTT68%vsuTTy(i)
          vsuz= pTT68%vsuTTz(i)
          cx  = pTT68%cTTx(i)
          cy  = pTT68%cTTy(i)
          cz  = pTT68%cTTz(i)
#endif
          this%F(i, 1, k) = this%F(i, 1, k) + pTT68%FX(i)
          this%F(i, 2, k) = this%F(i, 2, k) + pTT68%FY(i)
          this%F(i, 3, k) = this%F(i, 3, k) + pTT68%FZ(i)

#if  TRANS == 1

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
        
#endif
        end do
      end do

      ! Loop over charge sites in molecule
      do j = 1, this%Molecule%NCharge
        pCharge => this%Molecule%SiteCharge(j)
          do i = 1, np
            this%F(i, 1, k) = this%F(i, 1, k) + pCharge%FX(i)
            this%F(i, 2, k) = this%F(i, 2, k) + pCharge%FY(i)
            this%F(i, 3, k) = this%F(i, 3, k) + pCharge%FZ(i)
#if  TRANS == 1
       
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
           
#endif
          end do
      end do
    end if

    ! Reduce forces and torques from all processes
#if MPI_VER > 0
    call MPI_Reduce( this%F(:, :, :), this%FAll(:, :, :), size( this%F ), MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    if( this%Molecule%isElongated ) call MPI_Reduce( this%T(:, :, :), this%TAll(:, :, :), size( this%T ), &
&     MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )

#if  TRANS == 1
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
#endif
#endif

  end subroutine TComponent_Atom2Unit_Trans


#if OSMOP > 0
!==============================================================!
!  Subroutine TComponent_DensityProfile                        !
!==============================================================!

  subroutine TComponent_DensityProfile( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0 && !defined(MPI_USE_MODULE)
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
        if (this%Pm0(i,1) .ge. real(j-1)/NBinsDen-.5_RK) then
          if (this%Pm0(i,1) < real(j)/NBinsDen-.5_RK) then
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
    integer :: np, nra, iUnit
    integer :: i, j
    real(RK) :: r(this%NPart, 3), P0old

    ! Assign local variables
    np = this%NPart

    ! Predict COM positions and their derivatives
    do j = 1, 3
      do i = 1, np

        if (.not. UseIntDegFreed) then ! there should exist only one unit
            P0old = this%P0(i, j, 1)
        end if

        do iUnit = 1, this%Molecule%NUnit
          this%P0(i, j, iUnit) = this%P0(i, j, iUnit) + this%P1(i, j, iUnit) + this%P2(i, j, iUnit) + this%P3(i, j, iUnit) + this%P4(i, j, iUnit) + this%P5(i, j, iUnit)
          this%P1(i, j, iUnit) = this%P1(i, j, iUnit) + 2._RK * this%P2(i, j, iUnit) + 3._RK * this%P3(i, j, iUnit) + 4._RK * this%P4(i, j, iUnit) + 5._RK * this%P5(i, j, iUnit)
          this%P2(i, j, iUnit) = this%P2(i, j, iUnit) + 3._RK * this%P3(i, j, iUnit) + 6._RK * this%P4(i, j, iUnit) + 10._RK * this%P5(i, j, iUnit)
          this%P3(i, j, iUnit) = this%P3(i, j, iUnit) + 4._RK * this%P4(i, j, iUnit) +10._RK * this%P5(i, j, iUnit)
          this%P4(i, j, iUnit) = this%P4(i, j, iUnit) + 5._RK * this%P5(i, j, iUnit)
        end do

        ! Calculate displacement
        if (.not. UseIntDegFreed) then ! there should exist only one unit
            this%Disp(i, j) = this%Disp(i, j) + this%P0(i, j, 1) - P0old
        end if

        r(i, j) = 0._RK

        do iUnit = 1, this%Molecule%NUnit
          ! Check for conservation of particles in primary cell
#if ARCH == 1
          if( this%P0(i, j, iUnit) < -.5_RK ) then
            this%P0(i, j, iUnit) = this%P0(i, j, iUnit) + 1._RK
          elseif( this%P0(i, j, iUnit) > .5_RK ) then
            this%P0(i, j, iUnit) = this%P0(i, j, iUnit) - 1._RK
          end if
#else
          this%P0(i, j, iUnit) = this%P0(i, j, iUnit) - anint( this%P0(i, j, iUnit) )
#endif
          ! Calculate new positions of COM for molecules from new COM of units
          r(i, j) = r(i, j) + this%Molecule%Unit(iUnit)%Mass*(this%P0(i,j, iUnit)-anint(this%P0(i,j, iUnit)-this%Pm0(i,j)))
        end do

        if (UseIntDegFreed) then
          P0old = this%Pm0(i, j)
          this%Pm0(i,j) = r(i, j)/this%Molecule%Mass ! (unit-)mass averaged coordinates
          ! Calculate displacement of molecules
          this%Disp(i, j) = this%Disp(i, j) + this%Pm0(i, j) - P0old
          this%Pm0(i,j) = this%Pm0(i,j) - anint(this%Pm0(i,j))
        else ! there should exist only one unit
          this%Pm0(i,j) = this%P0(i, j, 1)
        end if
      end do
    end do

    do iUnit = 1, this%Molecule%NUnit
      if( this%Molecule%Unit(iUnit)%IsElongated ) then

        ! Predict quaternion parameters and their derivatives
        do j = 1, 4
          do i = 1, np
            this%Q0(i, j, iUnit) = this%Q0(i, j, iUnit) + this%Q1(i, j, iUnit) + this%Q2(i, j, iUnit) + this%Q3(i, j, iUnit) + this%Q4(i, j, iUnit)
            this%Q1(i, j, iUnit) = this%Q1(i, j, iUnit) + 2._RK * this%Q2(i, j, iUnit) + 3._RK * this%Q3(i, j, iUnit) + 4._RK * this%Q4(i, j, iUnit)
            this%Q2(i, j, iUnit) = this%Q2(i, j, iUnit) + 3._RK * this%Q3(i, j, iUnit) + 6._RK * this%Q4(i, j, iUnit)
            this%Q3(i, j, iUnit) = this%Q3(i, j, iUnit) + 4._RK * this%Q4(i, j, iUnit)
          end do
        end do

        nra = this%Molecule%Unit(iUnit)%NDFRot
        ! Predict angular velocities and their derivatives
        do j = 1, nra
          do i = 1, np
            this%W0(i, j, iUnit) = this%W0(i, j, iUnit) + this%W1(i, j, iUnit) + this%W2(i, j, iUnit) + this%W3(i, j, iUnit) + this%W4(i, j, iUnit)
            this%W1(i, j, iUnit) = this%W1(i, j, iUnit) + 2._RK * this%W2(i, j, iUnit) + 3._RK * this%W3(i, j, iUnit) + 4._RK * this%W4(i, j, iUnit)
            this%W2(i, j, iUnit) = this%W2(i, j, iUnit) + 3._RK * this%W3(i, j, iUnit) + 6._RK * this%W4(i, j, iUnit)
            this%W3(i, j, iUnit) = this%W3(i, j, iUnit) + 4._RK * this%W4(i, j, iUnit)
          end do
        end do

      end if

    end do

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
    real(RK)          :: MassInv, P0old
    real(RK)          :: Moi23, Moi31, Moi12
    real(RK)          :: TMoi1, TMoi2, TMoi3
    real(RK), pointer, contiguous :: pF(:, :, :), pT(:, :, :)
    integer           :: np, nra,iUnit
    integer           :: i, j
    real(RK)          :: r(this%NPart, 3)

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength
    np = this%NPart

    ! Correct COM positions and their derivatives
#if MPI_VER > 0
    pF => this%FAll(:, :, :)
#else
    pF => this%F(:, :, :)
#endif

    do j = 1, 3
      do i = 1, np
        r(i, j) = 0._RK

        if (.not. UseIntDegFreed) then
            P0old = this%P0(i, j, 1)
        end if

        do iUnit= 1, this%Molecule%NUnit
          if (.not. UseIntDegFreed) then
              MassInv = 1._RK / this%Molecule%Mass
          else
              MassInv = 1._RK / this%Molecule%Unit(iUnit)%Mass
          end if
          this%Corr0(i, j, iUnit) = pF(i, j, iUnit) * TimeStepSquared2 * BoxLengthInv * MassInv

          if( ConstantPressure .and. .not. NVTEquilibration ) this%Corr0(i, j, iUnit) = this%Corr0(i, j, iUnit) &
&             - this%P1(i, j, iUnit) * dLogVolumeThird

          this%Corr1(i, j, iUnit) = this%Corr0(i, j, iUnit) - this%P2(i, j, iUnit)
          this%P0(i, j, iUnit) = this%P0(i, j, iUnit) + this%Corr1(i, j, iUnit) * Gear20
          this%P1(i, j, iUnit) = this%P1(i, j, iUnit) + this%Corr1(i, j, iUnit) * Gear21
          this%P2(i, j, iUnit) =                    this%Corr0(i, j, iUnit)
          this%P3(i, j, iUnit) = this%P3(i, j, iUnit) + this%Corr1(i, j, iUnit) * Gear23
          this%P4(i, j, iUnit) = this%P4(i, j, iUnit) + this%Corr1(i, j, iUnit) * Gear24
          this%P5(i, j, iUnit) = this%P5(i, j, iUnit) + this%Corr1(i, j, iUnit) * Gear25

          if (.not. UseIntDegFreed) then ! there should exist only one unit
              this%Disp(i, j) = this%Disp(i, j) + this%P0(i, j, 1) - P0old
          end if

          ! Check for conservation of particles in primary cell
#if ARCH == 1
          if( this%P0(i, j, iUnit) < -.5_RK ) then
            this%P0(i, j, iUnit) = this%P0(i, j, iUnit) + 1._RK
          elseif( this%P0(i, j, iUnit) > .5_RK ) then
            this%P0(i, j, iUnit) = this%P0(i, j, iUnit) - 1._RK
          end if
#else
          this%P0(i, j, iUnit) = this%P0(i, j, iUnit) - anint( this%P0(i, j, iUnit) )
#endif
          ! Calculate new positions of COM for molecules from new COM of units
          if (.not. UseIntDegFreed) then
              r(i, j) = r(i, j) + this%Molecule%Mass*(this%P0(i,j,iUnit)-anint(this%P0(i,j,iUnit)-this%Pm0(i,j)))
          else
              r(i, j) = r(i, j) + this%Molecule%Unit(iUnit)%Mass*(this%P0(i,j,iUnit)-anint(this%P0(i,j,iUnit)-this%Pm0(i,j)))
          end if
        end do

        if (UseIntDegFreed) then
            P0old = this%Pm0(i, j)
            this%Pm0(i, j) = r(i, j)/this%Molecule%Mass
            ! Calculate displacement of molecules
            this%Disp(i, j) = this%Disp(i, j) + this%Pm0(i, j) - P0old
            this%Pm0(i, j) = this%Pm0(i, j) - anint(this%Pm0(i, j))
        else ! there should exist only one unit
            this%Pm0(i, j) = this%P0(i, j, 1)
        end if
      end do
    end do

    do iUnit = 1, this%Molecule%NUnit
      if( this%Molecule%Unit(iUnit)%isElongated ) then
      nra = this%Molecule%Unit(iUnit)%NDFRot

    ! Correct quaternion parameters and their derivatives
        do i = 1, np
          this%Corr0(i, 1, iUnit) = TimeStep2 * ( - this%Q0(i, 2, iUnit) * this%W0(i, 1, iUnit) - this%Q0(i, 3, iUnit) * this%W0(i, 2, iUnit) &
&                                             - this%Q0(i, 4, iUnit) * this%W0(i, 3, iUnit))

          this%Corr0(i, 2, iUnit) = TimeStep2 * ( + this%Q0(i, 1, iUnit) * this%W0(i, 1, iUnit) - this%Q0(i, 4, iUnit) * this%W0(i, 2, iUnit) &
&                                             + this%Q0(i, 3, iUnit) * this%W0(i, 3, iUnit))

          this%Corr0(i, 3, iUnit) = TimeStep2 * ( + this%Q0(i, 4, iUnit) * this%W0(i, 1, iUnit) + this%Q0(i, 1, iUnit) * this%W0(i, 2, iUnit) &
&                                             - this%Q0(i, 2, iUnit) * this%W0(i, 3, iUnit))

          this%Corr0(i, 4, iUnit) = TimeStep2 * ( - this%Q0(i, 3, iUnit) * this%W0(i, 1, iUnit) + this%Q0(i, 2, iUnit) * this%W0(i, 2, iUnit) &
&                                             + this%Q0(i, 1, iUnit) * this%W0(i, 3, iUnit))

        end do
        do j = 1, 4
          do i = 1, np
            this%Corr1(i, j, iUnit) = this%Corr0(i, j, iUnit) - this%Q1(i, j, iUnit)
            this%Q0(i, j, iUnit) = this%Q0(i, j, iUnit) + this%Corr1(i, j, iUnit) * Gear10
            this%Q1(i, j, iUnit) =                    this%Corr0(i, j, iUnit)
            this%Q2(i, j, iUnit) = this%Q2(i, j, iUnit) + this%Corr1(i, j, iUnit) * Gear12
            this%Q3(i, j, iUnit) = this%Q3(i, j, iUnit) + this%Corr1(i, j, iUnit) * Gear13
            this%Q4(i, j, iUnit) = this%Q4(i, j, iUnit) + this%Corr1(i, j, iUnit) * Gear14
          end do
        end do

        ! Correct angular velocities and their derivatives
#if MPI_VER > 0
        pT => this%TAll(:, :, :)
#else
        pT => this%T(:, :, :)
#endif
        TMoi1 = TimeStep / this%Molecule%Unit(iUnit)%MOI(1)
        TMoi2 = TimeStep / this%Molecule%Unit(iUnit)%MOI(2)

        if( this%Molecule%Unit(iUnit)%is3D ) then
          Moi23 = this%Molecule%Unit(iUnit)%MOI(2) - this%Molecule%Unit(iUnit)%MOI(3)
          Moi31 = this%Molecule%Unit(iUnit)%MOI(3) - this%Molecule%Unit(iUnit)%MOI(1)
          Moi12 = this%Molecule%Unit(iUnit)%MOI(1) - this%Molecule%Unit(iUnit)%MOI(2)
          TMoi3 = TimeStep / this%Molecule%Unit(iUnit)%MOI(3)
          do i = 1, np
            this%Corr0(i, 1, iUnit) = (pT(i, 1, iUnit) + this%W0(i, 2, iUnit) * this%W0(i, 3, iUnit) * Moi23) * TMoi1
            this%Corr0(i, 2, iUnit) = (pT(i, 2, iUnit) + this%W0(i, 3, iUnit) * this%W0(i, 1, iUnit) * Moi31) * TMoi2
            this%Corr0(i, 3, iUnit) = (pT(i, 3, iUnit) + this%W0(i, 1, iUnit) * this%W0(i, 2, iUnit) * Moi12) * TMoi3
          end do
        else
          do i = 1, np
            this%Corr0(i, 1, iUnit) = pT(i, 1, iUnit) * TMoi1
            this%Corr0(i, 2, iUnit) = pT(i, 2, iUnit) * TMoi2
          end do
        end if

        do j = 1, nra
          do i = 1, np
            this%Corr1(i, j, iUnit) = this%Corr0(i, j, iUnit) - this%W1(i, j, iUnit)
            this%W0(i, j, iUnit) = this%W0(i, j, iUnit) + this%Corr1(i, j, iUnit) * Gear10
            this%W1(i, j, iUnit) = this%Corr0(i, j, iUnit)
            this%W2(i, j, iUnit) = this%W2(i, j, iUnit) + this%Corr1(i, j, iUnit) * Gear12
            this%W3(i, j, iUnit) = this%W3(i, j, iUnit) + this%Corr1(i, j, iUnit) * Gear13
            this%W4(i, j, iUnit) = this%W4(i, j, iUnit) + this%Corr1(i, j, iUnit) * Gear14
          end do
        end do

      end if
    end do

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
    real(RK) :: Korr, P0old
    integer  :: np, nra, iUnit
    integer  :: i, j
    real(RK) :: r(this%NPart, 3)

    Korr = 2._RK - 1._RK / scale
    np = this%NPart

    do j = 1, 3
      do i = 1, np

        if (.not. UseIntDegFreed) then
            P0old = this%P0(i, j, 1)
        end if

        do iUnit = 1, this%Molecule%NUnit
          this%P1(i, j, iUnit) = Korr * this%P1(i, j, iUnit) + this%P2(i, j, iUnit)
          this%P0(i, j, iUnit) = this%P0(i, j, iUnit) + this%P1(i, j, iUnit)
        end do

        if (.not. UseIntDegFreed) then
            ! Calculate displacement
            this%Disp(i, j) = this%Disp(i, j) + this%P0(i, j, 1) - P0old
        end if

        r(i, j) = 0._RK

        do iUnit= 1, this%Molecule%NUnit
          ! Check for conservation of particles in primary cell
#if ARCH == 1
          if( this%P0(i, j, iUnit) < -.5_RK ) then
            this%P0(i, j, iUnit) = this%P0(i, j, iUnit) + 1._RK
          elseif( this%P0(i, j, iUnit) > .5_RK ) then
            this%P0(i, j, iUnit) = this%P0(i, j, iUnit) - 1._RK
          end if
#else
          this%P0(i, j, iUnit) = this%P0(i, j, iUnit) - anint( this%P0(i, j, iUnit) )
#endif
          ! Calculate new positions of COM for molecules from new COM of units
          r(i, j) = r(i, j) + this%Molecule%Unit(iUnit)%Mass*(this%P0(i,j,iUnit)-anint(this%P0(i,j,iUnit)-this%Pm0(i,j)))


          if (UseIntDegFreed) then
              P0old = this%Pm0(i, j)
              this%Pm0(i, j) = r(i, j)/this%Molecule%Mass
              ! Calculate displacement of molecules
              this%Disp(i, j) = this%Disp(i, j) + this%Pm0(i, j) - P0old
              this%Pm0(i, j) = this%Pm0(i,j) - anint(this%Pm0(i,j))
          else
              this%Pm0(i, j) = this%P0(i,j,1)
          end if
        end do
      end do
    end do

    do iUnit = 1, this%Molecule%NUnit
      if( this%Molecule%Unit(iUnit)%IsElongated ) then
        nra = this%Molecule%Unit(iUnit)%NDFRot
        do i = 1, np
          do j = 1, 4
            this%Q0tmp(i, j, iUnit) = this%Q0(i, j, iUnit) + .5_RK * this%Q1(i, j, iUnit)
          end do

          do j = 1, nra
            this%W0(i, j, iUnit) = Korr * this%W0(i, j, iUnit) + .5_RK * this%W1(i, j, iUnit)
          end do

          this%Q1(i, 1, iUnit) = TimeStep2 * ( - this%Q0tmp(i, 2, iUnit) * this%W0(i, 1, iUnit) - this%Q0tmp(i, 3, iUnit) * this%W0(i, 2, iUnit) &
&                                          - this%Q0tmp(i, 4, iUnit) * this%W0(i, 3, iUnit))

          this%Q1(i, 2, iUnit) = TimeStep2 * ( + this%Q0tmp(i, 1, iUnit) * this%W0(i, 1, iUnit) - this%Q0tmp(i, 4, iUnit) * this%W0(i, 2, iUnit) &
&                                          + this%Q0tmp(i, 3, iUnit) * this%W0(i, 3, iUnit))

          this%Q1(i, 3, iUnit) = TimeStep2 * ( + this%Q0tmp(i, 4, iUnit) * this%W0(i, 1, iUnit) + this%Q0tmp(i, 1, iUnit) * this%W0(i, 2, iUnit) &
&                                          - this%Q0tmp(i, 2, iUnit) * this%W0(i, 3, iUnit))

          this%Q1(i, 4, iUnit) = TimeStep2 * ( - this%Q0tmp(i, 3, iUnit) * this%W0(i, 1, iUnit) + this%Q0tmp(i, 2, iUnit) * this%W0(i, 2, iUnit) &
&                                          + this%Q0tmp(i, 1, iUnit) * this%W0(i, 3, iUnit))

          do j = 1, 4
            this%Q0(i, j, iUnit) = this%Q0(i, j, iUnit) + this%Q1(i, j, iUnit)
          end do
        end do
      end if
    end do

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
    real(RK), pointer, contiguous :: pF(:, :, :), pT(:, :, :)
    real(RK)          :: BoxLengthInv, MassInv
    real(RK)          :: Moi23, Moi31, Moi12
    real(RK)          :: TMoi1, TMoi2, TMoi3
    integer           :: np, nra, iUnit
    integer           :: i, j

    BoxLengthInv = 1._RK / this%BoxLength
    np = this%NPart
#if MPI_VER > 0
    pF => this%FAll(:, :, :)
#else
    pF => this%F(:, :, :)
#endif

    do iUnit = 1, this%Molecule%NUnit
      MassInv = 1._RK / this%Molecule%Unit(iUnit)%Mass
      do j = 1, 3
        do i = 1, np
          this%P2(i, j, iUnit) = pF(i, j, iUnit) * TimeStepSquared2 * BoxLengthInv * MassInv &
&                            - (this%P1(i, j, iUnit) + this%P2(i, j, iUnit)) * dLogVolumeThird

          this%P1(i, j, iUnit) = this%P1(i, j, iUnit) + this%P2(i, j, iUnit)
        end do
      end do
    end do

    if( this%Molecule%isElongated ) then
#if MPI_VER > 0
      pT => this%TAll(:, :, :)
#else
      pT => this%T(:, :, :)
#endif
    end if

    do iUnit = 1, this%Molecule%NUnit
      if( this%Molecule%Unit(iUnit)%IsElongated ) then
        nra = this%Molecule%Unit(iUnit)%NDFRot
        TMoi1 = TimeStep / this%Molecule%Unit(iUnit)%MOI(1)
        TMoi2 = TimeStep / this%Molecule%Unit(iUnit)%MOI(2)

        if( this%Molecule%Unit(iUnit)%is3D ) then
          TMoi3 = TimeStep / this%Molecule%Unit(iUnit)%MOI(3)
          Moi23 = this%Molecule%Unit(iUnit)%MOI(2) - this%Molecule%Unit(iUnit)%MOI(3)
          Moi31 = this%Molecule%Unit(iUnit)%MOI(3) - this%Molecule%Unit(iUnit)%MOI(1)
          Moi12 = this%Molecule%Unit(iUnit)%MOI(1) - this%Molecule%Unit(iUnit)%MOI(2)
          do i = 1, np
            this%W1(i, 1, iUnit) = (pT(i, 1, iUnit) + this%W0(i, 2, iUnit) * this%W0(i, 3, iUnit) * Moi23) * TMoi1
            this%W1(i, 2, iUnit) = (pT(i, 2, iUnit) + this%W0(i, 3, iUnit) * this%W0(i, 1, iUnit) * Moi31) * TMoi2
            this%W1(i, 3, iUnit) = (pT(i, 3, iUnit) + this%W0(i, 1, iUnit) * this%W0(i, 2, iUnit) * Moi12) * TMoi3
          end do

        else
          do i = 1, np
            this%W1(i, 1, iUnit) = pT(i, 1, iUnit) * TMoi1
            this%W1(i, 2, iUnit) = pT(i, 2, iUnit) * TMoi2
          end do
        end if

        do j = 1, nra
          do i = 1, np
            this%W0(i, j, iUnit) = this%W0(i, j, iUnit) + .5_RK * this%W1(i, j, iUnit)
          end do
        end do
        do i = 1, np
          this%Q1(i, 1, iUnit) = TimeStep2 * ( - this%Q0(i, 2, iUnit) * this%W0(i, 1, iUnit) - this%Q0(i, 3, iUnit) * this%W0(i, 2, iUnit) &
&                                              - this%Q0(i, 4, iUnit) * this%W0(i, 3, iUnit))

          this%Q1(i, 2, iUnit) = TimeStep2 * ( + this%Q0(i, 1, iUnit) * this%W0(i, 1, iUnit) - this%Q0(i, 4, iUnit) * this%W0(i, 2, iUnit) &
&                                              + this%Q0(i, 3, iUnit) * this%W0(i, 3, iUnit))

          this%Q1(i, 3, iUnit) = TimeStep2 * ( + this%Q0(i, 4, iUnit) * this%W0(i, 1, iUnit) + this%Q0(i, 1, iUnit) * this%W0(i, 2, iUnit) &
&                                              - this%Q0(i, 2, iUnit) * this%W0(i, 3, iUnit))

          this%Q1(i, 4, iUnit) = TimeStep2 * ( - this%Q0(i, 3, iUnit) * this%W0(i, 1, iUnit)  + this%Q0(i, 2, iUnit) * this%W0(i, 2, iUnit) &
&                                              + this%Q0(i, 1, iUnit) * this%W0(i, 3, iUnit))

        end do
      end if
    end do

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
!  Subroutine TComponent_AddParticle                           !
!==============================================================!

  subroutine TComponent_AddParticle( this, r, q )

    implicit none

    ! Declare arguments
    type(TComponent)               :: this
    real(RK), intent(in)           :: r(3)
    real(RK), intent(in), optional :: q(4)

    ! Declare local variables
    integer            :: selected, i

    ! Test boundaries of particle arrays
    if( this%NPart >= this%NPartMax .and. EnsembleType .eq. EnsembleTypeGE ) then
      tooManyParticles = .true.
      return
    end if

    if( this%NPart >= this%NPartMax .and. EnsembleType .eq. EnsembleTypeMUVT ) then
      tooManyParticles = .true.
      return
    end if

    if( this%NPart > this%NPartMax .and. EnsembleType .ne. EnsembleTypeGE) then
      tooManyParticles = .true.
      return
    end if

    if( this%NPart > this%NPartMax .and. EnsembleType .ne. EnsembleTypeMUVT) then
      tooManyParticles = .true.
      return
    end if

    ! Increase NPart
    if (UseIntDegFreed) selected = rnd( this%NPart )
    !IOBuffer = '' ! Michael DEBUG
    !write( IOBuffer, '("Duplicating particle number " I4," for insertion at end. ")' ) selected
    this%NPart = this%NPart + 1
#if MPI_VER > 0
    this%NPart1 = ProcRange( this%NPart, this%NPart0, this%NPart2 )
#endif

    ! Set coordinates and orientation of new particle by an representative of the configuration
    do i = 1, this%Molecule%NUnit
        if (.not. UseIntDegFreed) then
            this%P0(this%NPart,1:3,i) = r(1:3)
            this%Pm0(this%NPart,1:3) = this%P0(this%NPart,1:3,i)
        else
            this%P0(this%NPart,1:3,i) = this%P0(selected,1:3,i) + r(1:3)
            this%P0(this%NPart,1:3,i) = this%P0(this%NPart,1:3,i) - anint(this%P0(this%NPart,1:3,i))
        end if
    end do

    if (SimulationType .eq. MolecularDynamics ) then
      this%P1(this%NPart,:,:) = this%P1(selected,:,:)
      this%P2(this%NPart,:,:) = 0._RK
      if( this%Molecule%isElongated ) then
        this%W0(this%NPart,:,:) = this%W0(selected,:,:)
        this%Q1(this%NPart,:,:) = 0._RK
        this%W1(this%NPart,:,:) = 0._RK
      end if
      if ( IntegratorType .eq. IntegratorTypeGear) then
        this%P3(this%NPart,:,:) = 0._RK
        this%P4(this%NPart,:,:) = 0._RK
        this%P5(this%NPart,:,:) = 0._RK
        if( this%Molecule%isElongated ) then
          this%Q2(this%NPart,:,:) = 0._RK
          this%Q3(this%NPart,:,:) = 0._RK
          this%Q4(this%NPart,:,:) = 0._RK
          this%W2(this%NPart,:,:) = 0._RK
          this%W3(this%NPart,:,:) = 0._RK
          this%W4(this%NPart,:,:) = 0._RK
        end if
      end if
    end if

    if (UseIntDegFreed) then
        this%Pm0(this%NPart,:) = 0._RK
        call Unit2Mol(this, this%NPart) 
    end if

    if (this%Molecule%isElongated) then
        if (.not. UseIntDegFreed) then
            this%Q0(this%NPart,:,1) = q(:)
        else
            this%Q0(this%NPart,:,:) = this%Q0(selected,:,:)
            call RotateMol( this, this%NPart, q)
        end if
    end if
    if (UseIntDegFreed) call Unit2Atom1( this, this%NPart, RootProc )

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
      this%Pm0(np, :) = this%Pm0(this%NPart, :)
      this%P0(np, :, :) = this%P0(this%NPart, :, :)
      if( this%Molecule%isElongated ) then
        this%Q0(np, :, :) = this%Q0(this%NPart, :, :)
      end if

      if (SimulationType .eq. MolecularDynamics ) then
        this%P1(np,:,:) = this%P1(this%NPart,:,:)
        this%P2(np,:,:) = this%P2(this%NPart,:,:)
        if( this%Molecule%isElongated ) then
          this%Q1(np,:,:) = this%Q1(this%NPart,:,:)
          this%W1(np,:,:) = this%W1(this%NPart,:,:)
        end if
        if ( IntegratorType .eq. IntegratorTypeGear) then
          this%P3(np,:,:) = this%P3(this%NPart,:,:)
          this%P4(np,:,:) = this%P4(this%NPart,:,:)
          this%P5(np,:,:) = this%P5(this%NPart,:,:)
          if( this%Molecule%isElongated ) then
            this%Q2(np,:,:) = this%Q2(this%NPart,:,:)
            this%Q3(np,:,:) = this%Q3(this%NPart,:,:)
            this%Q4(np,:,:) = this%Q4(this%NPart,:,:)
            this%W2(np,:,:) = this%W2(this%NPart,:,:)
            this%W3(np,:,:) = this%W3(this%NPart,:,:)
            this%W4(np,:,:) = this%W4(this%NPart,:,:)
          end if
        end if
      end if

      ! Calculate Unit / Atom positions
      call Unit2Atom1( this, np, RootProc )

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
    integer          :: i

    ! Save current state
    this%P0Save = this%P0
    do i = 1, this%Molecule%NUnit
      if( this%Molecule%Unit(i)%IsElongated ) this%Q0Save = this%Q0
    end do

  end subroutine TComponent_SaveState



!==============================================================!
!  Subroutine TComponent_RestoreState                          !
!==============================================================!

  subroutine TComponent_RestoreState( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this
    integer          :: i

    ! Restore saved state
    this%P0 = this%P0Save
    do i = 1, this%Molecule%NUnit
      if( this%Molecule%Unit(i)%IsElongated ) this%Q0 = this%Q0Save
    end do

    ! Calculate site positions
! #if MPI_VER > 0
!     call Mol2Atom( this, this%NPart0, this%NPart2-this%NPart+1 )
! #else
!     call Mol2Atom( this, 1, this%NPart )
! #endif
    call Unit2Atom( this, this%NPart)

  end subroutine TComponent_RestoreState



!==============================================================!
!  Subroutine TComponent_RestartSave                           !
!==============================================================!

  subroutine TComponent_RestartSave( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Declare local variables
    integer  :: np, i, k, nu, j
    real(RK) :: pos(3), quat(4)

    ! Assign local variables
    np = this%NPart
    nu = this%Molecule%NUnit

    ! Check for root process
    if( .not. RootProc ) return

    ! Save contents to restart file
    write( restartFile%iounit, '(I10)' ) np
    if ( UseIntDegFreed ) write( restartFile%iounit, '(I10)' ) nu

    ! Centers of mass positions
    do i = 1, np
      do k = 1, nu
        pos(:) = this%P0(i,:, k)
        write( restartFile%iounit, '(3(ES20.12E3, :, ";"))' ) pos(:)
      end do
    end do

    if( SimulationType .eq. MolecularDynamics ) then
      ! Centers of mass positions' derivatives
      do i = 1, np
        do k = 1, nu
          pos(:) = this%P1(i,:, k)
          write( restartFile%iounit, '(3(ES20.12E3, :, ";"))' ) pos(:)
        end do
      end do
      do i = 1, np
        do k = 1, nu
          pos(:) = this%P2(i,:, k)
          write( restartFile%iounit, '(3(ES20.12E3, :, ";"))' ) pos(:)
        end do
      end do

      if( IntegratorType .eq. IntegratorTypeGear ) then
        do i = 1, np
          do k = 1, nu
            pos(:) = this%P3(i,:, k)
            write( restartFile%iounit, '(3(ES20.12E3, :, ";"))' ) pos(:)
          end do
        end do
        do i = 1, np
          do k = 1, nu
            pos(:) = this%P4(i,:, k)
            write( restartFile%iounit, '(3(ES20.12E3, :, ";"))' ) pos(:)
          end do
        end do
        do i = 1, np
          do k = 1, nu
            pos(:) = this%P5(i,:, k)
            write( restartFile%iounit, '(3(ES20.12E3, :, ";"))' ) pos(:)
          end do
        end do
      end if

      if (.not. printIDF) then
          do i = 1, np
            pos(:) = this%Disp(i,:)
            write( restartFile%iounit, '(3(ES20.12E3, :, ";"))' ) pos(:)
          end do

          if( ALPHA2UpdateFrequency > 0 ) then
            do i = 1, np
              do j = 0, ALPHA2Length/ALPHA2Shift-1
                write( restartFile%iounit, '(3(ES20.12E3, :, ";"))' ) this%ri0_x(i,j),this%ri0_y(i,j),this%ri0_z(i,j)
              end do
            end do
          end if
      end if

#if TRANS == 1
      if( (TransMethod .eq. Einstein) .or. (TransMethod .eq. GKEinstein) ) then  !EinsteinCoef ri0_E rest write
            do i = 1, np
              do j = 0, this%NEinstein-1
                write( restartFile%iounit, '(3(ES20.12E3, :, ";"))' ) this%ri0_E_x(i,j),this%ri0_E_y(i,j),this%ri0_E_z(i,j)
              end do
            end do
       end if
#endif

    else
      write( restartFile%iounit, '(ES20.12E3)' ) this%DispTran
      write( restartFile%iounit, '(2I10)' ) this%NMoveAttempts, this%NMoveSuccesses
      write( restartFile%iounit, '(2I10)' ) this%NMoveBiasedAttempts, &
&       this%NMoveBiasedSuccesses
      if ( UseIntDegFreed ) then
        write( iounit_restart, '(ES20.12E3)' ) this%DispMolTran
        write( iounit_restart, '(2I10)' ) this%NMoveMolAttempts, this%NMoveMolSuccesses
        write( iounit_restart, '(2I10)' ) this%NMoveBiasedMolAttempts, this%NMoveBiasedMolSuccesses
      end if
    end if

    if( this%Molecule%isElongated ) then
      ! Quaternion parameters
      do i = 1, np
        do k = 1, nu
          quat(:) = this%Q0(i,:, k)
          write( restartFile%iounit, '(4(ES20.12E3, :, ";"))' ) quat(:)
        end do
      end do

      if( SimulationType .eq. MolecularDynamics ) then
        ! Quaternion parameters' derivatives
        do i = 1, np
          do k = 1, nu
            quat(:) = this%Q1(i,:, k)
            write( restartFile%iounit, '(4(ES20.12E3, :, ";"))' ) quat(:)
          end do
        end do

        if( IntegratorType .eq. IntegratorTypeGear ) then
          do i = 1, np
            do k = 1, nu
              quat(:) = this%Q2(i,:, k)
              write( restartFile%iounit, '(4(ES20.12E3, :, ";"))' ) quat(:)
            end do
          end do
          do i = 1, np
            do k = 1, nu
              quat(:) = this%Q3(i,:, k)
              write( restartFile%iounit, '(4(ES20.12E3, :, ";"))' ) quat(:)
            end do
          end do
          do i = 1, np
            do k = 1, nu
              quat(:) = this%Q4(i,:, k)
              write( restartFile%iounit, '(4(ES20.12E3, :, ";"))' ) quat(:)
            end do
          end do
        end if

        ! Angular velocities and their derivatives
        do i = 1, np
          do k = 1, nu
            pos(:) = this%W0(i,:, k)
            write( restartFile%iounit, '(3(ES20.12E3, :, ";"))' ) pos(:)
          end do
        end do
        do i = 1, np
          do k = 1, nu
            pos(:) = this%W1(i,:, k)
            write( restartFile%iounit, '(3(ES20.12E3, :, ";"))' ) pos(:)
          end do
        end do

        if( IntegratorType .eq. IntegratorTypeGear ) then
          do i = 1, np
            do k = 1, nu
              pos(:) = this%W2(i,:, k)
              write( restartFile%iounit, '(3(ES20.12E3, :, ";"))' ) pos(:)
            end do
          end do
          do i = 1, np
            do k = 1, nu
              pos(:) = this%W3(i,:, k)
              write( restartFile%iounit, '(3(ES20.12E3, :, ";"))' ) pos(:)
            end do
          end do
          do i = 1, np
            do k = 1, nu
              pos(:) = this%W4(i,:, k)
              write( restartFile%iounit, '(3(ES20.12E3, :, ";"))' ) pos(:)
            end do
          end do
        end if
      else
        write( restartFile%iounit, '(ES20.12E3)' ) this%DispRot
        write( restartFile%iounit, '(2I10)' ) this%NRotateAttempts, this%NRotateSuccesses
        write( restartFile%iounit, '(2I10)' ) this%NRotateBiasedAttempts, this%NRotateBiasedSuccesses
        if ( UseIntDegFreed ) then
          write( restartFile%iounit, '(ES20.12E3)' ) this%DispMolRot
          write( restartFile%iounit, '(2I10)' ) this%NRotateMolAttempts, this%NRotateMolSuccesses
          write( restartFile%iounit, '(2I10)' ) this%NRotateBiasedMolAttempts, this%NRotateBiasedMolSuccesses
        end if
      end if

    end if

    if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
      write( restartFile%iounit, '(ES20.12E3)' ) this%WF(:)
      write( restartFile%iounit, '(I10)' ) this%NState(:)
      write( restartFile%iounit, '(I10)' ) this%NStateWF(:)
    end if

  end subroutine TComponent_RestartSave



!==============================================================!
!  Subroutine TComponent_RestartRead                           !
!==============================================================!

  subroutine TComponent_RestartRead( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0 && !defined(MPI_USE_MODULE)
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent) :: this

    ! Declare local variables
    integer :: i, np, nu, k, j
    real(RK):: r(3)

    if( RootProc ) then

      ! Read contents from restart file
      read( restartFile%iounit, '(I10)' ) np
      if( np > this%NPartMax ) call Error( 'Not enough memory to read particles from restart file' )
      this%NPart = np
      if (UseIntDegFreed) then
        read( iounit_restart, '(I10)' ) nu
        this%Molecule%NUnit = nu
      else
        nu = 1
        this%Molecule%NUnit = nu
      end if

      ! Centers of mass positions
      do i = 1, np
        do k = 1, nu
          read( restartFile%iounit, '(3(ES20.12E3, :, 1X))' ) (this%P0( i, j, k ),j=1,3)
        end do
      end do

      ! Calculate positions of COM for molecules from  COM of units...needed for intialization
      this%Pm0(:,:) = 0._RK
      do j = 1, 3 ! 3 iterations to find close to real COM...3 is randomly chosen to account for high impulses (without 1 should sufficie)
        do i = 1, np
          r(:) = 0._RK
          do k= 1, nu
            ! Calculate new positions of COM for molecules from new COM of units
            r(1:3) = r(1:3) + this%Molecule%Unit(k)%Mass*this%P0(i, 1:3, k)
          end do
          this%Pm0(i,:) = r(:)/this%Molecule%Mass
          ! Calculate displacement of molecules
          this%Pm0(i,:) = this%Pm0(i,:) - anint(this%Pm0(i,:))
        end do
      end do

      if( SimulationType .eq. MolecularDynamics ) then
        ! Centers of mass positions' derivatives
        do i = 1, np
          do k = 1, nu
            read( restartFile%iounit, '(3(ES20.12E3, :, 1X))' ) (this%P1( i, j, k),j=1,3)
          end do
        end do

        do i = 1, np
          do k = 1, nu
            read( restartFile%iounit, '(3(ES20.12E3, :, 1X))' ) (this%P2( i, j, k),j=1,3)
          end do
        end do

        if( IntegratorType .eq. IntegratorTypeGear ) then
          do i = 1, np
            do k = 1, nu
              read( restartFile%iounit, '(3(ES20.12E3, :, 1X))' ) (this%P3( i, j, k),j=1,3)
            end do
          end do

          do i = 1, np
            do k = 1, nu
              read( restartFile%iounit, '(3(ES20.12E3, :, 1X))' ) (this%P4( i, j, k),j=1,3)
            end do
          end do

          do i = 1, np
            do k = 1, nu
              read( restartFile%iounit, '(3(ES20.12E3, :, 1X))' ) (this%P5( i, j, k),j=1,3)
            end do
          end do
        end if

        if (.not. printIDF) then
            do i = 1, np
              read( restartFile%iounit, '(3(ES20.12E3, :, 1X))' ) (this%Disp( i, j ),j=1,3)
            end do

            if( ALPHA2UpdateFrequency > 0 ) then
              do i = 1, np
                do j = 0, ALPHA2Length/ALPHA2Shift-1
                  read( restartFile%iounit, '(3(ES20.12E3, :, 1X))' ) this%ri0_x(i,j),this%ri0_y(i,j),this%ri0_z(i,j)
                end do
              end do
            end if
        end if
#if TRANS==1
        !EinsteinCoef rest read ri0_E_x
         if( (TransMethod .eq. GKEinstein) .or. (TransMethod .eq. Einstein) ) then
             do i = 1, np
               do j = 0, this%NEinstein-1
                 read( restartFile%iounit, '(3(ES20.12E3, :, 1X))' ) this%ri0_E_x(i,j),this%ri0_E_y(i,j),this%ri0_E_z(i,j)
               end do
             end do
         end if
#endif
      else
        read( restartFile%iounit, '(ES20.12E3)' ) this%DispTran
        read( restartFile%iounit, '(2I10)' ) this%NMoveAttempts, this%NMoveSuccesses
        read( restartFile%iounit, '(2I10)' ) this%NMoveBiasedAttempts, this%NMoveBiasedSuccesses
        if ( UseIntDegFreed ) then
          read( restartFile%iounit, '(ES20.12E3)' ) this%DispMolTran
          read( restartFile%iounit, '(2I10)' ) this%NMoveMolAttempts, this%NMoveMolSuccesses
          read( restartFile%iounit, '(2I10)' ) this%NMoveBiasedMolAttempts, this%NMoveBiasedMolSuccesses
        end if
      end if

      if( this%Molecule%isElongated ) then
        ! Quaternion parameters
        do i = 1, np
          do k = 1, nu
            read( restartFile%iounit, '(4(ES20.12E3, :, 1X))' ) (this%Q0( i, j, k),j=1,4)
          end do
        end do

        if( SimulationType .eq. MolecularDynamics ) then
          ! Quaternion parameters' derivatives
          do i = 1, np
            do k = 1, nu
              read( restartFile%iounit, '(4(ES20.12E3, :, 1X))' ) (this%Q1( i, j, k),j=1,4)
            end do
          end do

          if( IntegratorType .eq. IntegratorTypeGear ) then
            do i = 1, np
              do k = 1, nu
                read( restartFile%iounit, '(4(ES20.12E3, :, 1X))' ) (this%Q2( i, j, k),j=1,4)
              end do
            end do

            do i = 1, np
              do k = 1, nu
                read( restartFile%iounit, '(4(ES20.12E3, :, 1X))' ) (this%Q3( i, j, k),j=1,4)
              end do
            end do

            do i = 1, np
              do k = 1, nu
                read( restartFile%iounit, '(4(ES20.12E3, :, 1X))' ) (this%Q4( i, j, k),j=1,4)
              end do
            end do
          end if

          ! Angular velocities and their derivatives
          do i = 1, np
            do k = 1, nu
              read( restartFile%iounit, '(3(ES20.12E3, :, 1X))' ) (this%W0( i, j, k),j=1,3)
            end do
          end do

          do i = 1, np
            do k = 1, nu
              read( restartFile%iounit, '(3(ES20.12E3, :, 1X))' ) (this%W1( i, j, k),j=1,3)
            end do
          end do

          if( IntegratorType .eq. IntegratorTypeGear ) then
            do i = 1, np
              do k = 1, nu
                read( restartFile%iounit, '(3(ES20.12E3, :, 1X))' ) (this%W2( i, j, k),j=1,3)
              end do
            end do

            do i = 1, np
              do k = 1, nu
                read( restartFile%iounit, '(3(ES20.12E3, :, 1X))' ) (this%W3( i, j, k),j=1,3)
              end do
            end do

            do i = 1, np
              do k = 1, nu
                read( restartFile%iounit, '(3(ES20.12E3, :, 1X))' ) (this%W4( i, j, k),j=1,3)
              end do
            end do

          end if

        else
          read( restartFile%iounit, '(ES20.12E3)' ) this%DispRot
          read( restartFile%iounit, '(2I10)' ) this%NRotateAttempts, this%NRotateSuccesses
          read( restartFile%iounit, '(2I10)' ) this%NRotateBiasedAttempts, this%NRotateBiasedSuccesses
          if ( UseIntDegFreed ) then
            read( restartFile%iounit, '(ES20.12E3)' ) this%DispMolRot
            read( restartFile%iounit, '(2I10)' ) this%NRotateMolAttempts, this%NRotateMolSuccesses
            read( restartFile%iounit, '(2I10)' ) this%NRotateBiasedMolAttempts, this%NRotateBiasedMolSuccesses
          end if
        end if
      end if

      if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
        read( restartFile%iounit, '(ES20.12E3)' ) this%WF(:)
        read( restartFile%iounit, '(I10)' ) this%NState(:)
        read( restartFile%iounit, '(I10)' ) this%NStateWF(:)
      end if

    end if

#if MPI_VER > 0
    call MPI_Bcast( this%NPart, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
    call MPI_Bcast( this%Pm0(:, :), size( this%Pm0 ), MPI_RK, NRootProc, Communicator, ierror )
    call MPI_Bcast( this%P0(:, :, :), size( this%P0 ), MPI_RK, NRootProc, Communicator, ierror )
    if( this%Molecule%isElongated ) call MPI_Bcast( this%Q0(:, :, :), size( this%Q0 ), MPI_RK, &
&       NRootProc, Communicator, ierror )

    if( (SimulationType .eq. MonteCarlo) .or. (SimulationType .eq. Gibbs) ) then

      call MPI_Bcast( this%DispTran, 1, MPI_RK, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%NMoveAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%NMoveSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%NMoveBiasedAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%NMoveBiasedSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
      if ( UseIntDegFreed ) then
        call MPI_Bcast( this%DispMolTran, 1, MPI_RK, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NMoveMolAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NMoveMolSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NMoveBiasedMolAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NMoveBiasedMolSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
      end if

      if( this%Molecule%isElongated ) then
        call MPI_Bcast( this%DispRot, 1, MPI_RK, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NRotateAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NRotateSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NRotateBiasedAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NRotateBiasedSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        if ( UseIntDegFreed ) then
          call MPI_Bcast( this%DispRot, 1, MPI_RK, NRootProc, Communicator, ierror )
          call MPI_Bcast( this%NRotateAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
          call MPI_Bcast( this%NRotateSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
          call MPI_Bcast( this%NRotateBiasedAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
          call MPI_Bcast( this%NRotateBiasedSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        end if
      end if
    end if

    if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
      call MPI_Bcast( this%WF, size( this%WF ), MPI_RK, NRootProc, Communicator, ierror )
      call MPI_BCast( this%NState, size( this%NState ), MPI_INTEGER, NRootProc, Communicator, ierror )
      call MPI_BCast( this%NStateWF, size( this%NStateWF ), MPI_INTEGER, NRootProc, Communicator, ierror )
    end if
#endif

  end subroutine TComponent_RestartRead


!TRANSPORT_start
!==============================================================!
!  Subroutine TComponent_ForceTransport                        !
!==============================================================!

subroutine TComponent_ForceTransport( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0 && !defined(MPI_USE_MODULE)
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
       !Michael Sch.: this only works for 1unit per molecule...no P1 and W0 defined for COM only for the units !!!
        do k = 1, 3
          do i= 1, this%Npart  
            this%FTC(i,1)= this%FTC(i,1)+ pFTC1(i,k)*this%P1(i,k, 1)
            this%FTC(i,2)= this%FTC(i,2)+ pFTC2(i,k)*this%P1(i,k, 1)
            this%FTC(i,3)= this%FTC(i,3)+ pFTC3(i,k)*this%P1(i,k, 1)
          end do
        end do

        nra = this%Molecule%Unit(1)%NDFRot
        do k= 1, nra
          do i= 1, this%Npart
            this%FRC(i,1)= this%FRC(i,1)+ pFRC1(i,k)*this%W0(i,k, 1)
            this%FRC(i,2)= this%FRC(i,2)+ pFRC2(i,k)*this%W0(i,k, 1)
            this%FRC(i,3)= this%FRC(i,3)+ pFRC3(i,k)*this%W0(i,k, 1)
          end do
        end do

        this%FTC(:,:) = this%FTC(:,:)*BoxLength_dt


        ! Calculate kinetic energy / molecule
        do j = 1, 3
          this%KinETran(:,j) = this%P1(:,j, 1)*this%P1(:,j, 1)
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


!==============================================================!
!  Subroutine TComponent_ReverseLeapFrog                       !
!==============================================================!

  subroutine TComponent_ReverseLeapFrog( this, oldF, dLogVolumeThird )

    implicit none

    ! Declare arguments
    type(TComponent)       :: this
    real(RK), intent( in ) :: oldF(this%NPart,3,this%Molecule%NUnit)
    real(RK), intent( in ) :: dLogVolumeThird

    ! Declare local variables
    integer           :: nu, np
    integer           :: i, j, k
    real(RK)          :: r(3), BoxLengthInv, P0old(3)

    BoxLengthInv = 1._RK / this%BoxLength
    np = this%NPart
    nu = this%Molecule%NUnit

    do k = 1, nu
      do j = 1, 3
        do i = 1, np
          this%P0(i, j, k) = this%P0(i, j, k) - this%P1(i, j, k)
          this%P1(i, j, k) = this%P1(i, j, k) - 2._RK*this%P2(i, j, k)
          if (abs(dLogVolumeThird) > 0._RK) then
            this%P2(i, j, k) = (oldF(i, j, k) * TimeStepSquared2 * BoxLengthInv / this%Molecule%Unit(k)%Mass &
&                               - this%P2(i, j, k) ) / dLogVolumeThird - this%P1(i, j, k)
          ! else oldP2 is not needed and new value can correctly be calculated without previous one
          endif
        end do
      end do
    end do

    do i = 1, np
      r(:) = 0._RK
      do k= 1, nu
        do j = 1, 3
          ! Check for conservation of particles in primary cell
#if ARCH == 1
          if( this%P0(i, j, k) < -.5_RK ) then
            this%P0(i, j, k) = this%P0(i, j, k) + 1._RK
          elseif( this%P0(i, j, k) > .5_RK ) then
            this%P0(i, j, k) = this%P0(i, j, k) - 1._RK
          end if
#else
          this%P0(i, j, k) = this%P0(i, j, k) - anint( this%P0(i, j, k) )
#endif
          ! Calculate new positions of COM for molecules from new COM of units
          r(j) = r(j) + this%Molecule%Unit(k)%Mass*(this%P0(i,j,k)-anint(this%P0(i,j,k)-this%Pm0(i,j)))
        end do
      end do

      P0old = this%Pm0(i, :)

      this%Pm0(i,:) = r(:)/this%Molecule%Mass
      ! Calculate displacement of molecules
      this%Disp(i, :) = this%Disp(i, :) + this%Pm0(i, :) - P0old
      this%Pm0(i,:) = this%Pm0(i,:) - anint(this%Pm0(i,:))
    end do

    do k = 1, nu
      if( this%Molecule%Unit(k)%IsElongated ) then
        do j = 1, 4
          this%Q0(1:np, j, k) = this%Q0(1:np, j, k) - this%Q1(1:np, j, k)
        end do
        do j = 1, this%Molecule%Unit(k)%NDFRot
          this%W0(1:np, j, k) = this%W0(1:np, j, k) - this%W1(1:np, j, k)
        end do
      end if
    end do

  end subroutine TComponent_ReverseLeapFrog


!==============================================================!
!  Subroutine TComponent_Constraints                           !
!==============================================================!

 subroutine TComponent_Constraints( this, VirialShake)

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent)        :: this
    real(RK),intent(in out) :: VirialShake

    ! Declare local variables
    logical                 :: need
    type(TIdfBond), pointer :: pBond
    integer                 :: it, itmax, unstableMol
    integer                 :: i, i0, i1, j, k
    integer                 :: np, nu, Unit1, Unit2
    real(RK)                :: BoxLength, Shake2, Shake002
    real(RK)                :: RX1, RY1, RZ1, RX2, RY2, RZ2
    real(RK)                :: PX1, PY1, PZ1, PX2, PY2, PZ2
    real(RK)                :: RXij(this%Molecule%NBond), RYij(this%Molecule%NBond), RZij(this%Molecule%NBond)
    real(RK)                :: R0Xij(this%Molecule%NBond), R0Yij(this%Molecule%NBond), R0Zij(this%Molecule%NBond)
    real(RK)                :: P0Xij(this%Molecule%NBond), P0Yij(this%Molecule%NBond), P0Zij(this%Molecule%NBond)
    real(RK)                :: dRmax, dRmaxold, RSquared, R0Sq, dRSquared(this%Molecule%NBond)
    real(RK)                :: PR1(this%Molecule%NBond,3), PR2(this%Molecule%NBond,3)
    real(RK)                :: e(this%Molecule%NBond,3), EffM(this%Molecule%NBond)
    real(RK)                :: fx, fy, fz, Fijconstr, EMass1, EMass2, Coeff, Coeffalt, DotProd
    real(RK)                :: Term1(3), Term2(3), Term3(3), MOI(3), intermedQ0(4), addW1(3), addQ1(4)
    real(RK)                :: tempP0(3,this%Molecule%NUnit), tempQ0(4,this%Molecule%NUnit), tempW0(3)
    real(RK)                :: tempF(3,this%Molecule%NUnit), tempT(3,this%Molecule%NUnit)
#if MPI_VER > 0
    integer                 :: itRoot, unstableMolRoot

    call MPI_Bcast( this%P0(:, :, :), size( this%P0 ), MPI_RK, NRootProc, Communicator, ierror )
    call MPI_Bcast( this%P1(:, :, :), size( this%P1 ), MPI_RK, NRootProc, Communicator, ierror )
    if( this%Molecule%isElongated ) then
      call MPI_Bcast( this%Q0(:, :, :), size( this%Q0 ), MPI_RK, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%Q1(:, :, :), size( this%Q1 ), MPI_RK, NRootProc, Communicator, ierror )
    end if
#endif

    ! Assign local variables
    BoxLength = this%BoxLength
    Shake2 = Shake**2
    Shake002 = Shake2/100
    np = this%NPart
    nu = this%Molecule%NUnit
    unstableMol = 0
    itmax = 9999

#if MPI_VER > 0
    i0 = this%NPart0
    i1 = this%NPart2
molloop: do i = i0, i1
#else
    i1 = this%NPart
molloop: do i = 1, i1
#endif
!     do i = 1, np
      ! get/save old site and COM (of unit) positions
      tempP0(:,:) = this%P0(i,:,:) - this%P1(i,:,:)
      tempQ0(:,:) = this%Q0(i,:,:) - this%Q1(i,:,:)
      need = .true.
      it = 0
      dRmax = 0._RK

      do j = 1, this%Molecule%NBond
        pBond => this%Molecule%IdfBond(j)
        Unit1 = pBond%UnitId1
        Unit2 = pBond%UnitId2

        RX1 = pBond%RX1(i)
        RY1 = pBond%RY1(i)
        RZ1 = pBond%RZ1(i)
        RX2 = pBond%RX2(i)
        RY2 = pBond%RY2(i)
        RZ2 = pBond%RZ2(i)

        PX1 = tempP0(1,Unit1)
        PY1 = tempP0(2,Unit1)
        PZ1 = tempP0(3,Unit1)
        PX2 = tempP0(1,Unit2)
        PY2 = tempP0(2,Unit2)
        PZ2 = tempP0(3,Unit2)

        ! Calculate unit-unit and site-site distance vector at begin of this timestep
        R0Xij(j) = RX1 - RX2
        R0Yij(j) = RY1 - RY2
        R0Zij(j) = RZ1 - RZ2
        P0Xij(j) = PX1 - PX2
        P0Yij(j) = PY1 - PY2
        P0Zij(j) = PZ1 - PZ2
        R0Xij(j) = (R0Xij(j) - anint( R0Xij(j) )) * BoxLength
        R0Yij(j) = (R0Yij(j) - anint( R0Yij(j) )) * BoxLength
        R0Zij(j) = (R0Zij(j) - anint( R0Zij(j) )) * BoxLength
        P0Xij(j) = (P0Xij(j) - anint( P0Xij(j) )) * BoxLength
        P0Yij(j) = (P0Yij(j) - anint( P0Yij(j) )) * BoxLength
        P0Zij(j) = (P0Zij(j) - anint( P0Zij(j) )) * BoxLength

        e(j,1) = R0Xij(j)/pBond%R0
        e(j,2) = R0Yij(j)/pBond%R0
        e(j,3) = R0Zij(j)/pBond%R0

        ! Calculate unit-site distance vector's at start of this timestep
        if (this%Molecule%Unit(Unit1)%IsElongated) then
         PR1(j,1) = (RX1 - PX1) * BoxLength
         PR1(j,2) = (RY1 - PY1) * BoxLength
         PR1(j,3) = (RZ1 - PZ1) * BoxLength
        endif
        if (this%Molecule%Unit(Unit2)%IsElongated) then
         PR2(j,1) = (RX2 - PX2) * BoxLength
         PR2(j,2) = (RY2 - PY2) * BoxLength
         PR2(j,3) = (RZ2 - PZ2) * BoxLength
        endif

        EMass1 = 1._RK/this%Molecule%Unit(Unit1)%Mass
        EMass2 = 1._RK/this%Molecule%Unit(Unit2)%Mass

        ! Contribution of the units to the effective mass; changing Shake to QShake (Forester & Smith, 1998)
        ! e_ = d_/d = (RXij/R0, RYij/R0, RZij/R0)
        ! sA_ = PR1(i,j,:)
        ! IA_ = (MOI(1) 0 0; 0 MOI(2) 0; 0 0 MOI(3) )
        if (this%Molecule%Unit(Unit1)%IsElongated) then
          MOI(:) = this%Molecule%Unit(Unit1)%MOI(:)

          Term3(1) = PR1(j,2)*e(j,3) - PR1(j,3)*e(j,2)
          Term3(2) = PR1(j,3)*e(j,1) - PR1(j,1)*e(j,3)
          Term3(3) = PR1(j,1)*e(j,2) - PR1(j,2)*e(j,1)

          call qTerm( tempQ0(:, Unit1), MOI(:), Term3(:), Term2(:) )

          Term1(1) = Term2(2)*PR1(j,3) - Term2(3)*PR1(j,2)
          Term1(2) = Term2(3)*PR1(j,1) - Term2(1)*PR1(j,3)
          Term1(3) = Term2(1)*PR1(j,2) - Term2(2)*PR1(j,1)

          EMass1 = EMass1 + Term1(1)*e(j,1) + Term1(2)*e(j,2) + Term1(3)*e(j,3)
        endif
        if (this%Molecule%Unit(Unit2)%IsElongated) then
          MOI(:) = this%Molecule%Unit(Unit2)%MOI(:)

          Term3(1) = PR2(j,2)*e(j,3) - PR2(j,3)*e(j,2)
          Term3(2) = PR2(j,3)*e(j,1) - PR2(j,1)*e(j,3)
          Term3(3) = PR2(j,1)*e(j,2) - PR2(j,2)*e(j,1)

          call qTerm( tempQ0(:, Unit2), MOI(:), Term3(:), Term2(:) )

          Term1(1) = Term2(2)*PR2(j,3) - Term2(3)*PR2(j,2)
          Term1(2) = Term2(3)*PR2(j,1) - Term2(1)*PR2(j,3)
          Term1(3) = Term2(1)*PR2(j,2) - Term2(2)*PR2(j,1)

          EMass2 = EMass2 + Term1(1)*e(j,1) + Term1(2)*e(j,2) + Term1(3)*e(j,3)
        endif

        EffM(j)=1._RK/(EMass1+EMass2)/TimeStepSquared

      end do ! bond loop for t0 calculation

      !calculates site positions of unconstrained timestep
      call Unit2AtomShake(this, i, this%P0(i,1:3,1:nu), this%Q0(i,1:4,1:nu))

      ! Loop over all bonds in molecule
      do j = 1, this%Molecule%NBond
        pBond => this%Molecule%IdfBond(j)

        R0Sq = pBond%R0**2 ! squared equlibrium bond length
        RX1 = pBond%RX1(i)
        RY1 = pBond%RY1(i)
        RZ1 = pBond%RZ1(i)
        RX2 = pBond%RX2(i)
        RY2 = pBond%RY2(i)
        RZ2 = pBond%RZ2(i)

        ! Calculate temporary bond vector
        RXij(j) = RX1 - RX2
        RYij(j) = RY1 - RY2
        RZij(j) = RZ1 - RZ2
        RXij(j) = (RXij(j) - anint( RXij(j) )) * BoxLength
        RYij(j) = (RYij(j) - anint( RYij(j) )) * BoxLength
        RZij(j) = (RZij(j) - anint( RZij(j) )) * BoxLength

        RSquared=RXij(j)**2+RYij(j)**2+RZij(j)**2

        ! Deviation from equilibrium
        dRSquared(j) = (R0Sq - RSquared)
        dRmax=max(dRmax,abs(dRSquared(j)/R0Sq))
      end do
      if (dRmax < Shake2) then ! molecule already inside of constraints even without QShake
        need = .false.
      else
        dRmaxold = 1.1 * dRmax
      end if

      do while (need .and. ( it < itmax )) ! calculate shake-force iteratively

        it = it+1
        tempF(:,:) = 0._RK
        tempT(:,:) = 0._RK

        do j= 1, this%Molecule%NBond
          pBond => this%Molecule%IdfBond(j)

          Unit1 = pBond%UnitId1
          Unit2 = pBond%UnitId2
          R0Sq = pBond%R0**2
          DotProd = R0Xij(j)*RXij(j)+R0Yij(j)*RYij(j)+R0Zij(j)*RZij(j)

          ! firs order equation of QShake
          ! (not monoton decreasing for all conformations, therefore not used)
          !Fijconstr = 0.5_RK*EffM(j)*dRSquared(j)/DotProd

          ! second order equation
          Coeff = dRSquared(j) + DotProd**2 / R0Sq
          Coeff = Coeff / R0Sq
          Fijconstr = - DotProd / R0Sq
          if ( Coeff > 0._RK ) then
            Coeff = sqrt(Coeff)
            Coeffalt  = Fijconstr - Coeff
            Fijconstr = Fijconstr + Coeff
            if ( abs(Coeffalt) < abs(Fijconstr) ) Fijconstr = Coeffalt
          end if
          Fijconstr = Fijconstr * EffM(j)

          fx = Fijconstr * R0Xij(j)
          fy = Fijconstr * R0Yij(j)
          fz = Fijconstr * R0Zij(j)

          VirialShake = VirialShake + ( fx*P0Xij(j) + fy*P0Yij(j) + fz*P0Zij(j) )

          tempF(1,Unit1) = tempF(1,Unit1) + fx
          tempF(2,Unit1) = tempF(2,Unit1) + fy
          tempF(3,Unit1) = tempF(3,Unit1) + fz
          tempF(1,Unit2) = tempF(1,Unit2) - fx
          tempF(2,Unit2) = tempF(2,Unit2) - fy
          tempF(3,Unit2) = tempF(3,Unit2) - fz
          if (this%Molecule%Unit(Unit1)%IsElongated) then
            ! Torque
            tempT(1,Unit1) = tempT(1,Unit1) + PR1(j,2)*fz - PR1(j,3)*fy
            tempT(2,Unit1) = tempT(2,Unit1) + PR1(j,3)*fx - PR1(j,1)*fz
            tempT(3,Unit1) = tempT(3,Unit1) + PR1(j,1)*fy - PR1(j,2)*fx
          end if
          if (this%Molecule%Unit(Unit2)%IsElongated) then
            ! Torque
            tempT(1,Unit2) = tempT(1,Unit2) - PR2(j,2)*fz + PR2(j,3)*fy
            tempT(2,Unit2) = tempT(2,Unit2) - PR2(j,3)*fx + PR2(j,1)*fz
            tempT(3,Unit2) = tempT(3,Unit2) - PR2(j,1)*fy + PR2(j,2)*fx
          end if

        end do ! constraint force calculation for all bonds

        do j=1,nu ! displacement of all units due to forces
          Coeff = TimeStepSquared / BoxLength / this%Molecule%Unit(j)%Mass
          ! Translational Correction
          this%P0(i, 1, j) = this%P0(i, 1, j) + Coeff*tempF(1,j)
          this%P0(i, 2, j) = this%P0(i, 2, j) + Coeff*tempF(2,j)
          this%P0(i, 3, j) = this%P0(i, 3, j) + Coeff*tempF(3,j)
          this%F(i,1:3,j) = this%F(i,1:3,j) + tempF(1:3,j)

          ! Rotational Correction
          if (this%Molecule%Unit(j)%IsElongated) then ! new term
            addW1(:) = 0._RK
            ! Changes to Rotational Matrix due to QShake
            this%T(i,:,j) = this%T(i,:,j) + tempT(:,j)
            addW1(1) = tempT(1,j) / this%Molecule%Unit(j)%MOI(1)
            addW1(2) = tempT(2,j) / this%Molecule%Unit(j)%MOI(2)
            if( this%Molecule%Unit(j)%is3D ) then
              addW1(3) = tempT(3,j) / this%Molecule%Unit(j)%MOI(3)
            end if
            addW1(:) = 0.5_RK*TimeStep * addW1(:)
            tempW0(:) = addW1(:)
            addQ1(1) = TimeStep2 * ( - this%Q0(i, 2, j) * tempW0(1) &
&                                    - this%Q0(i, 3, j) * tempW0(2) - this%Q0(i, 4, j) * tempW0(3))
            addQ1(2) = TimeStep2 * ( + this%Q0(i, 1, j) * tempW0(1) &
&                                    - this%Q0(i, 4, j) * tempW0(2) + this%Q0(i, 3, j) * tempW0(3))
            addQ1(3) = TimeStep2 * ( + this%Q0(i, 4, j) * tempW0(1) &
&                                    + this%Q0(i, 1, j) * tempW0(2) - this%Q0(i, 2, j) * tempW0(3))
            addQ1(4) = TimeStep2 * ( - this%Q0(i, 3, j) * tempW0(1) &
&                                    + this%Q0(i, 2, j) * tempW0(2) + this%Q0(i, 1, j) * tempW0(3))
            intermedQ0(:) = this%Q0(i, :, j) + 0.5_RK*addQ1(:)
            tempW0(:) = tempW0(:) + addW1(:)
            addQ1(1) = TimeStep2 * ( - intermedQ0(2) * tempW0(1) &
&                                    - intermedQ0(3) * tempW0(2) - intermedQ0(4) * tempW0(3))
            addQ1(2) = TimeStep2 * ( + intermedQ0(1) * tempW0(1) &
&                                    - intermedQ0(4) * tempW0(2) + intermedQ0(3) * tempW0(3))
            addQ1(3) = TimeStep2 * ( + intermedQ0(4) * tempW0(1) &
&                                    + intermedQ0(1) * tempW0(2) - intermedQ0(2) * tempW0(3))
            addQ1(4) = TimeStep2 * ( - intermedQ0(3) * tempW0(1) &
&                                    + intermedQ0(2) * tempW0(2) + intermedQ0(1) * tempW0(3))
            this%Q0(i, :, j) = this%Q0(i, :, j) + addQ1(:)
          end if

        end do ! unit loop

        call Unit2AtomShake(this, i, this%P0(i,1:3,1:nu), this%Q0(i,1:4,1:nu))
        dRmaxold = dRmax
        dRmax = 0._RK

        !Loop over all bonds in molecule
        do j = 1, this%Molecule%NBond
          pBond => this%Molecule%IdfBond(j)

          R0Sq = pBond%R0**2
          RX1 = pBond%RX1(i)
          RY1 = pBond%RY1(i)
          RZ1 = pBond%RZ1(i)
          RX2 = pBond%RX2(i)
          RY2 = pBond%RY2(i)
          RZ2 = pBond%RZ2(i)

          ! Calculate temporary bond vector
          RXij(j) = RX1 - RX2
          RYij(j) = RY1 - RY2
          RZij(j) = RZ1 - RZ2
          RXij(j) = (RXij(j) - anint( RXij(j) )) * BoxLength
          RYij(j) = (RYij(j) - anint( RYij(j) )) * BoxLength
          RZij(j) = (RZij(j) - anint( RZij(j) )) * BoxLength

          RSquared=RXij(j)**2+RYij(j)**2+RZij(j)**2

          ! Deviation from equilibrium
          dRSquared(j) = (R0Sq - RSquared)
          dRmax=max(dRmax,abs(dRSquared(j)/R0Sq))
        end do

        if (dRmax < Shake2) then ! constraints inside precision
          need = .false.
        elseif ( dRmax >= dRmaxold .or. Shake002 > (dRmaxold - dRmax) ) then
          dRmax = 0.1_RK*Shake2
          need = .false.
          unstableMol = i
          ! happens for Force*timestep too high; QShake can't resolve constraints if displacements each timestep are very large
          ! checks to-do: timestep below 10fs and able to display dynamics? Average displacement < sigma?
          !dRmax = 0._RK
          exit molloop
        end if

      end do ! shake-force calculation for molecule i finished

      do j=1,nu
        do k=1,3
          this%P0(i,k,j) = tempP0(k,j) + this%P1(i,k,j)
        end do
        if (this%Molecule%Unit(j)%IsElongated) then
          do k=1,4
            this%Q0(i,k,j) = tempQ0(k,j) + this%Q1(i,k,j)
          end do
        endif
      end do

    end do molloop ! molecule loop

#if MPI_VER > 0
    call MPI_Reduce( this%F(:, :, :), this%FAll(:, :, :), size( this%F ), &
&     MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    if( this%Molecule%isElongated ) call MPI_Reduce( this%T(:, :, :), this%TAll(:, :, :), size( this%T ), &
&     MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    call MPI_Reduce( it, itRoot, 1, MPI_INTEGER, MPI_MAX, NRootProc, Communicator, ierror )
    call MPI_Reduce( unstableMol, unstableMolRoot, 1, MPI_INTEGER, MPI_MAX, NRootProc, Communicator, ierror )
    if (RootProc) then
      it = itRoot
      unstableMol = unstableMolRoot
#endif

      if (unstableMol > 0 .or. it >= itmax) then ! write debug information (pos, vel, force) in restart file
        call FileRewrite( iounit_restart, trim(RestartFileName) )
        write( iounit_restart, '(A)' ) trim( ParameterFileName )
        write( iounit_restart, '(2I10)' ) Step, StepTotal
        write( iounit_restart, '(2L5)' ) Equilibration, NVTEquilibration
        call RestartSave( this )
        do i = 1, np
          do k = 1, nu
            write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%F(np,:,nu)
          end do
        end do
        call FileClose( iounit_restart )
      end if

      if (unstableMol > 0) then
        write( IOBuffer, '("QShake was not converging to zero for molecule", I6, &
&                          " in step", I10, ". Stop at iteration: ", I4)' ) unstableMol, Step, it
        call LogWrite
        call Error( 'QShake was not converging. Occuring forces to high for integration.' )
      end if
      if ( it >= itmax ) then
        call Error( 'Initial density to high for QShake. Probably to many overlaps.' )
      end if
#if MPI_VER > 0
    end if
#endif


contains

    subroutine qTerm (q, MOI, term3, term2)

    ! Declare arguments
    real(RK), intent( in )  :: q(4)
    real(RK), intent( in )  :: MOI(3)
    real(RK), intent( in )  :: term3(3)
    real(RK), intent( out ) :: term2(3)

    ! Declare local variables
    real(RK)                :: q1, q2, q3, q4, qinv
    real(RK)                :: A11, A12, A13, A21, A22, A23, A31, A32, A33
    real(RK)                :: M11, M12, M13, M21, M22, M23, M31, M32, M33

    term2(:) = 0._RK

    ! Normalise quaternions
    q1 = q(1)
    q2 = q(2)
    q3 = q(3)
    q4 = q(4)
#if ARCH == 3
    qinv = rsqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#else
    qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#endif
    q1 = q1 * qinv
    q2 = q2 * qinv
    q3 = q3 * qinv
    q4 = q4 * qinv

    A11 = q1**2 + q2**2 - q3**2 - q4**2
    A12 = 2._RK * ( q2 * q3 + q2 * q4)
    A13 = 2._RK * ( q2 * q4 - q1 * q3)
    A21 = 2._RK * ( q2 * q3 - q1 * q4)
    A22 = q1**2 -  q2**2 + q3**2 - q4**2
    A23 = 2._RK * (q3 * q4 + q1 * q2)

    M11 = A11**2/MOI(1) + A21**2/MOI(2)
    M12 = A11*A12/MOI(1) + A21*A22/MOI(2)
    M13 = A11*A13/MOI(1) + A21*A23/MOI(2)
    M21 = A12*A11/MOI(1) + A22*A21/MOI(2)
    M22 = A12**2/MOI(1) + A22**2/MOI(2)
    M23 = A12*A13/MOI(1) + A22*A23/MOI(2)
    M31 = A13*A11/MOI(1) + A23*A21/MOI(2)
    M32 = A13*A12/MOI(1) + A23*A22/MOI(2)
    M33 = A13**2/MOI(1) + A23**2/MOI(2)

    if (MOI(3) > Zero) then
      A31 = 2._RK * (q2 * q4 + q1 * q3)
      A32 = 2._RK * (q3 * q4 - q1 * q2)
      A33 = q1**2 -  q2**2 - q3**2 + q4**2
      M11 = M11 + A31**2/MOI(3)
      M12 = M12 + A31*A32/MOI(3)
      M13 = M13 + A31*A33/MOI(3)
      M21 = M21 + A32*A31/MOI(3)
      M22 = M22 + A32**2/MOI(3)
      M23 = M23 + A32*A33/MOI(3)
      M31 = M31 + A33*A31/MOI(3)
      M32 = M32 + A33*A32/MOI(3)
      M33 = M33 + A33**2/MOI(3)
    end if

    term2(1) = M11*term3(1) + M12*term3(2) + M13*term3(3)
    term2(2) = M21*term3(1) + M22*term3(2) + M23*term3(3)
    term2(3) = M31*term3(1) + M32*term3(2) + M33*term3(3)

    end subroutine qTerm


 end subroutine TComponent_Constraints



!==============================================================!
!  Subroutine TComponent_DuplicateParticle                     !
!==============================================================!

  subroutine TComponent_DuplicateParticle( comp1, comp2, np)

    implicit none

    ! Declare arguments
    type(TComponent)       :: comp1, comp2
    integer, intent(in)    :: np

    ! Test boundaries of particle arrays
    if( comp1%NPart > comp1%NPartMax ) then
      tooManyParticles = .true.
      return
    elseif( comp1%NPart == comp1%NPartMax .and. EnsembleType .eq. EnsembleTypeGE ) then
      tooManyParticles = .true.
      return
    end if

    ! Increase NPart
    comp1%NPart = comp1%NPart + 1
#if MPI_VER > 0
    comp1%NPart1 = ProcRange( comp1%NPart, comp1%NPart0, comp1%NPart2 )
#endif

    ! Set coordinates and orientation of new particle
    comp1%Pm0(comp1%NPart, :) = comp2%Pm0(np, :)
    comp1%P0(comp1%NPart, :, :) = comp2%P0(np, :, :)
    if (comp1%Molecule%isElongated) then
      comp1%Q0(comp1%NPart, :, :) = comp2%Q0(np, :, :)
    end if

  end subroutine TComponent_DuplicateParticle


!===============================================================!
!  Subroutine TComponent_Unit2AtomShake (per molecule, w/o MPI) !
!===============================================================!

  subroutine TComponent_Unit2AtomShake( this, np, P0, Q0)

    implicit none

    ! Declare arguments
    type(TComponent)     :: this
    integer, intent(in)  :: np
    real(RK), intent(in) :: P0(3,this%Molecule%NUnit)
    real(RK), intent(in) :: Q0(4,this%Molecule%NUnit)

    ! Declare local variables
    real(RK)                       :: BoxLengthInv
    real(RK)                       :: q1, q2, q3, q4, qinv
    real(RK)                       :: A11, A12, A13
    real(RK)                       :: A21, A22, A23
    real(RK)                       :: A31, A32, A33
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    real(RK)                       :: mue1, mue2, mue3
    type(TSiteMIEnm), pointer      :: pLJ126
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j, k
    integer                        :: nu

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength
    nu = this%Molecule%NUnit

    do k = 1, nu
      ! Check number of rotation axes
      if( this%Molecule%Unit(k)%isElongated ) then
        ! Positions and quaternions of unit k in particle i
        q1 = Q0(1, k)
        q2 = Q0(2, k)
        q3 = Q0(3, k)
        q4 = Q0(4, k)

        ! Normalise quaternions
#if ARCH == 3
        qinv = rsqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#else
        qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#endif
        q1 = q1 * qinv
        q2 = q2 * qinv
        q3 = q3 * qinv
        q4 = q4 * qinv
        this%Q0(np,1,k) = q1
        this%Q0(np,2,k) = q2
        this%Q0(np,3,k) = q3
        this%Q0(np,4,k) = q4

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

        ! Loop over LJ126 sites in unit
        do j = 1, this%Molecule%Unit(k)%NMIEnm
          pLJ126 => this%Molecule%Unit(k)%SiteMIEnm(j)
          r1 = pLJ126%r(1) * BoxLengthInv
          r2 = pLJ126%r(2) * BoxLengthInv
          r3 = pLJ126%r(3) * BoxLengthInv
          pLJ126%RX(np) = P0(1, k) + r1 * A11 + r2 * A21 + r3 * A31
          pLJ126%RY(np) = P0(2, k) + r1 * A12 + r2 * A22 + r3 * A32
          pLJ126%RZ(np) = P0(3, k) + r1 * A13 + r2 * A23 + r3 * A33
        end do

        ! Loop over charge sites in molecule
        do j = 1, this%Molecule%Unit(k)%NCharge
          pCharge => this%Molecule%Unit(k)%SiteCharge(j)
          r1 = pCharge%r(1) * BoxLengthInv
          r2 = pCharge%r(2) * BoxLengthInv
          r3 = pCharge%r(3) * BoxLengthInv
          pCharge%RX(np) = P0(1, k) + r1 * A11 + r2 * A21 + r3 * A31
          pCharge%RY(np) = P0(2, k) + r1 * A12 + r2 * A22 + r3 * A32
          pCharge%RZ(np) = P0(3, k) + r1 * A13 + r2 * A23 + r3 * A33
        end do

        ! Loop over dipole sites in molecule
        do j = 1, this%Molecule%Unit(k)%NDipole
          pDipole => this%Molecule%Unit(k)%SiteDipole(j)
          r1 = pDipole%r(1) * BoxLengthInv
          r2 = pDipole%r(2) * BoxLengthInv
          r3 = pDipole%r(3) * BoxLengthInv
          or1 = pDipole%or(1)
          or2 = pDipole%or(2)
          or3 = pDipole%or(3)
          pDipole%RX(np) = P0(1, k) + r1 * A11 + r2 * A21 + r3 * A31
          pDipole%RY(np) = P0(2, k) + r1 * A12 + r2 * A22 + r3 * A32
          pDipole%RZ(np) = P0(3, k) + r1 * A13 + r2 * A23 + r3 * A33
          pDipole%OX(np) = or1 * A11 + or2 * A21 + or3 * A31
          pDipole%OY(np) = or1 * A12 + or2 * A22 + or3 * A32
          pDipole%OZ(np) = or1 * A13 + or2 * A23 + or3 * A33
        end do

        ! Loop over quadrupole sites in molecule
        do j = 1, this%Molecule%Unit(k)%NQuadrupole
          pQuadrupole => this%Molecule%Unit(k)%SiteQuadrupole(j)
          r1 = pQuadrupole%r(1) * BoxLengthInv
          r2 = pQuadrupole%r(2) * BoxLengthInv
          r3 = pQuadrupole%r(3) * BoxLengthInv
          or1 = pQuadrupole%or(1)
          or2 = pQuadrupole%or(2)
          or3 = pQuadrupole%or(3)
          pQuadrupole%RX(np) = P0(1, k) + r1 * A11 + r2 * A21 + r3 * A31
          pQuadrupole%RY(np) = P0(2, k) + r1 * A12 + r2 * A22 + r3 * A32
          pQuadrupole%RZ(np) = P0(3, k) + r1 * A13 + r2 * A23 + r3 * A33
          pQuadrupole%OX(np) = or1 * A11 + or2 * A21 + or3 * A31
          pQuadrupole%OY(np) = or1 * A12 + or2 * A22 + or3 * A32
          pQuadrupole%OZ(np) = or1 * A13 + or2 * A23 + or3 * A33
        end do

        if( CutoffMode .eq. CenterofMass ) then
          mue1 = this%Molecule%Unit(k)%Mue(1)
          mue2 = this%Molecule%Unit(k)%Mue(2)
          mue3 = this%Molecule%Unit(k)%Mue(3)
          this%MueX(np, k) = mue1 * A11 + mue2 * A21 + mue3 * A31
          this%MueY(np, k) = mue1 * A12 + mue2 * A22 + mue3 * A32
          this%MueZ(np, k) = mue1 * A13 + mue2 * A23 + mue3 * A33
        end if

      else !If unit is not elongated

        ! Loop over LJ126 sites in molecule
        do i = 1, this%Molecule%Unit(k)%NMIEnm
          pLJ126 => this%Molecule%Unit(k)%SiteMIEnm(i)
          pLJ126%RX(np) = P0(1, k)
          pLJ126%RY(np) = P0(2, k)
          pLJ126%RZ(np) = P0(3, k)
        end do

        ! Loop over charge sites in molecule
        do i = 1, this%Molecule%Unit(k)%NCharge
          pCharge => this%Molecule%Unit(k)%SiteCharge(i)
          pCharge%RX(np) = P0(1, k)
          pCharge%RY(np) = P0(2, k)
          pCharge%RZ(np) = P0(3, k)
        end do

      end if
    end do

  end subroutine TComponent_Unit2AtomShake


!==============================================================!
!  Subroutine TComponent_InitUnit                              !
!==============================================================!

subroutine TComponent_InitUnit( this, np, dq )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: np
    real(RK),intent(in) :: dq(3)

    ! Declare local variables
    type(TUnit), pointer :: pUnit
    real(RK)             :: BoxLengthInv
    real(RK)             :: A11, A12, A13
    real(RK)             :: A21, A22, A23
    real(RK)             :: A31, A32, A33
    real(RK)             :: q1, q2, q3, q4, qinv
    integer              :: i

    ! Broadcast center of mass (molecule) positions to all processes
#if MPI_VER > 0
    ! in MC simulations, we only communicate during common equilibration
    if ( SimulationType .ne. MonteCarlo .or. ((Equilibration .and. CommonEqui) )) then
      call MPI_Bcast( this%Pm0(np, :), 3, MPI_RK, NRootProc, Communicator, ierror )
    end if
#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    ! Calculate rotation matrix elements
    q1 = 1._RK
    q2 = dq(1)
    q3 = dq(2)
    q4 = dq(3)
    ! Normalise quaternions
#if ARCH == 3
    qinv = rsqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#else
    qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#endif
    q1 = q1 * qinv
    q2 = q2 * qinv
    q3 = q3 * qinv
    q4 = q4 * qinv

    A11 = q2**2 - q3**2 - q4**2 + q1**2
    A12 = 2._RK * (q2 * q3 + q4*q1)
    A13 = 2._RK * (q2 * q4 - q3*q1)
    A21 = 2._RK * (q2 * q3 - q4*q1)
    A22 = - q2**2 + q3**2 - q4**2 + q1**2
    A23 = 2._RK * (q3 * q4 + q2*q1)
    A31 = 2._RK * (q2 * q4 + q3*q1)
    A32 = 2._RK * (q3 * q4 - q2*q1)
    A33 = - q2**2 - q3**2 + q4**2 + q1**2

    do i=1,this%Molecule%NUnit
      pUnit => this%Molecule%Unit(i)

      ! Calculating new Positions and quaternions of unit i after rotation
      this%P0(np,1,i) = this%Pm0(np,1) + (pUnit%P0(1) * A11 + pUnit%P0(2) * A21 + pUnit%P0(3) * A31) * BoxLengthInv
      this%P0(np,2,i) = this%Pm0(np,2) + (pUnit%P0(1) * A12 + pUnit%P0(2) * A22 + pUnit%P0(3) * A32) * BoxLengthInv
      this%P0(np,3,i) = this%Pm0(np,3) + (pUnit%P0(1) * A13 + pUnit%P0(2) * A23 + pUnit%P0(3) * A33) * BoxLengthInv

      ! Unit%Q0*dq w/o norm
      this%Q0(np,1,i) = pUnit%Q0(1) - dq(1)*pUnit%Q0(2) - dq(2)*pUnit%Q0(3) - dq(3)*pUnit%Q0(4)
      this%Q0(np,2,i) = pUnit%Q0(2) + dq(1)*pUnit%Q0(1) - dq(3)*pUnit%Q0(3) + dq(2)*pUnit%Q0(4)
      this%Q0(np,3,i) = pUnit%Q0(3) + dq(2)*pUnit%Q0(1) + dq(3)*pUnit%Q0(2) - dq(1)*pUnit%Q0(4)
      this%Q0(np,4,i) = pUnit%Q0(4) + dq(3)*pUnit%Q0(1) - dq(2)*pUnit%Q0(2) + dq(1)*pUnit%Q0(3)
    end do

  end subroutine TComponent_InitUnit


!==============================================================!
!  Subroutine TComponent_RotateTest                            !
!==============================================================!

subroutine TComponent_RotateTest( this, np, dq )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: np
    real(RK),intent(in) :: dq(3)

    ! Declare local variables
    real(RK)         :: BoxLengthInv
    real(RK)         :: PX, PY, PZ
    real(RK)         :: A11, A12, A13
    real(RK)         :: A21, A22, A23
    real(RK)         :: A31, A32, A33
    real(RK)         :: r1, r2, r3, rm(3)
    real(RK)         :: q1, q2, q3, q4, qinv
    integer          :: i

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    ! Calculate rotation matrix elements
    q1 = 1._RK
    q2 = dq(1)
    q3 = dq(2)
    q4 = dq(3)
    ! Normalise quaternions
#if ARCH == 3
    qinv = rsqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#else
    qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#endif
    q1 = q1 * qinv
    q2 = q2 * qinv
    q3 = q3 * qinv
    q4 = q4 * qinv

    A11 = q2**2 - q3**2 - q4**2 + q1**2
    A12 = 2._RK * (q2 * q3 + q4*q1)
    A13 = 2._RK * (q2 * q4 - q3*q1)
    A21 = 2._RK * (q2 * q3 - q4*q1)
    A22 = - q2**2 + q3**2 - q4**2 + q1**2
    A23 = 2._RK * (q3 * q4 + q2*q1)
    A31 = 2._RK * (q2 * q4 + q3*q1)
    A32 = 2._RK * (q3 * q4 - q2*q1)
    A33 = - q2**2 - q3**2 + q4**2 + q1**2

    rm(:) = 0._RK
    do i = 1, this%Molecule%NUnit
      rm(1:3) = rm(1:3) + this%Molecule%Unit(i)%Mass*this%P0Test(np,1:3,i)
    end do
    rm(:) = rm(:)/this%Molecule%Mass

    do i=1,this%Molecule%NUnit
      ! Positions and quaternions of unit i in particle np
      PX = this%P0Test(np, 1, i)
      PY = this%P0Test(np, 2, i)
      PZ = this%P0Test(np, 3, i)
      q1 = this%Q0Test(np, 1, i)
      q2 = this%Q0Test(np, 2, i)
      q3 = this%Q0Test(np, 3, i)
      q4 = this%Q0Test(np, 4, i)

      ! Distance unit-COM
      r1 = (PX-rm(1))
      r2 = (PY-rm(2))
      r3 = (PZ-rm(3))

      ! Calculating new Positions and quaternions of unit i after rotation
      this%P0Test(np,1,i) = rm(1) + r1 * A11 + r2 * A21 + r3 * A31
      this%P0Test(np,2,i) = rm(2) + r1 * A12 + r2 * A22 + r3 * A32
      this%P0Test(np,3,i) = rm(3) + r1 * A13 + r2 * A23 + r3 * A33

      ! this%Q0*dq
      this%Q0Test(np,1,i) = q1 - dq(1)*q2 - dq(2)*q3 - dq(3)*q4
      this%Q0Test(np,2,i) = q2 + dq(1)*q1 - dq(3)*q3 + dq(2)*q4
      this%Q0Test(np,3,i) = q3 + dq(2)*q1 + dq(3)*q2 - dq(1)*q4
      this%Q0Test(np,4,i) = q4 + dq(3)*q1 - dq(2)*q2 + dq(1)*q3

    end do

  end subroutine TComponent_RotateTest


!==============================================================!
!  Subroutine TComponent_RotateMol                             !
!==============================================================!

subroutine TComponent_RotateMol( this, np, dq )

    implicit none

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: np
    real(RK),intent(in) :: dq(3)

    ! Declare local variables
    real(RK)         :: BoxLengthInv
    real(RK)         :: PX, PY, PZ
    real(RK)         :: A11, A12, A13
    real(RK)         :: A21, A22, A23
    real(RK)         :: A31, A32, A33
    real(RK)         :: r1, r2, r3
    real(RK)         :: q1, q2, q3, q4, qinv
    integer          :: i

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    ! Calculate rotation matrix elements
    q1 = 1._RK
    q2 = dq(1)
    q3 = dq(2)
    q4 = dq(3)
    ! Normalise quaternions
#if ARCH == 3
    qinv = rsqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#else
    qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#endif
    q1 = q1 * qinv
    q2 = q2 * qinv
    q3 = q3 * qinv
    q4 = q4 * qinv

    A11 = q2**2 - q3**2 - q4**2 + q1**2
    A12 = 2._RK * (q2 * q3 + q4*q1)
    A13 = 2._RK * (q2 * q4 - q3*q1)
    A21 = 2._RK * (q2 * q3 - q4*q1)
    A22 = - q2**2 + q3**2 - q4**2 + q1**2
    A23 = 2._RK * (q3 * q4 + q2*q1)
    A31 = 2._RK * (q2 * q4 + q3*q1)
    A32 = 2._RK * (q3 * q4 - q2*q1)
    A33 = - q2**2 - q3**2 + q4**2 + q1**2

    do i=1,this%Molecule%NUnit
      ! Positions and quaternions of unit i in particle np
      PX = this%P0(np, 1, i)
      PY = this%P0(np, 2, i)
      PZ = this%P0(np, 3, i)
      q1 = this%Q0(np, 1, i)
      q2 = this%Q0(np, 2, i)
      q3 = this%Q0(np, 3, i)
      q4 = this%Q0(np, 4, i)

      ! Distance unit-COM
      r1 = (PX-this%Pm0(np,1))
      r2 = (PY-this%Pm0(np,2))
      r3 = (PZ-this%Pm0(np,3))
      r1 = r1 - anint(r1)
      r2 = r2 - anint(r2)
      r3 = r3 - anint(r3)

      ! Calculating new Positions and quaternions of unit i after rotation
      this%P0(np,1,i) = this%Pm0(np,1) + r1 * A11 + r2 * A21 + r3 * A31
      this%P0(np,2,i) = this%Pm0(np,2) + r1 * A12 + r2 * A22 + r3 * A32
      this%P0(np,3,i) = this%Pm0(np,3) + r1 * A13 + r2 * A23 + r3 * A33
      this%P0(np,1,i) = this%P0(np,1,i) - anint(this%P0(np,1,i))
      this%P0(np,2,i) = this%P0(np,2,i) - anint(this%P0(np,2,i))
      this%P0(np,3,i) = this%P0(np,3,i) - anint(this%P0(np,3,i))

      ! this%Q0*dq
      this%Q0(np,1,i) = q1 - dq(1)*q2 - dq(2)*q3 - dq(3)*q4
      this%Q0(np,2,i) = q2 + dq(1)*q1 - dq(3)*q3 + dq(2)*q4
      this%Q0(np,3,i) = q3 + dq(2)*q1 + dq(3)*q2 - dq(1)*q4
      this%Q0(np,4,i) = q4 + dq(3)*q1 - dq(2)*q2 + dq(1)*q3

    end do

  end subroutine TComponent_RotateMol


!==============================================================!
!  Subroutine TComponent_ResizeMol                             !
!==============================================================!

  subroutine TComponent_ResizeMol( this, DelBoxFrac )

    implicit none

    ! Declare arguments
    real(RK),intent(in) :: DelBoxFrac
    type(TComponent)    :: this

    ! Declare local variables
    real(RK)            :: PXij, PYij,PZij
    integer             :: nu, np
    integer             :: i, j


    ! Calculate positions of units after global resize
    nu = this%Molecule%NUnit
    np = this%NPart

    do i=1, np
      do j=1,nu
        PXij = this%P0(i,1,j) - this%Pm0(i,1)
        PYij = this%P0(i,2,j) - this%Pm0(i,2)
        PZij = this%P0(i,3,j) - this%Pm0(i,3)

        this%P0(i,1,j) = ( PXij - anint(PXij) ) / DelBoxFrac + this%Pm0(i,1)
        this%P0(i,2,j) = ( PYij - anint(PYij) ) / DelBoxFrac + this%Pm0(i,2)
        this%P0(i,3,j) = ( PZij - anint(PZij) ) / DelBoxFrac + this%Pm0(i,3)

        this%P0(i,1,j) = this%P0(i,1,j) - anint(this%P0(i,1,j))
        this%P0(i,2,j) = this%P0(i,2,j) - anint(this%P0(i,2,j))
        this%P0(i,3,j) = this%P0(i,3,j) - anint(this%P0(i,3,j))
      end do
    end do

    end subroutine TComponent_ResizeMol


!==============================================================!
!  Subroutine TComponent_Unit2Atom1Root (per mol + if root)    !
!==============================================================!

  subroutine TComponent_Unit2Atom1Root( this, np, Proc )

    implicit none

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: np
    logical, intent(in) :: Proc

    ! Declare local variables
    real(RK)                       :: BoxLengthInv
    real(RK)                       :: PX, PY, PZ
    real(RK)                       :: q1, q2, q3, q4, qinv
    real(RK)                       :: A11, A12, A13
    real(RK)                       :: A21, A22, A23
    real(RK)                       :: A31, A32, A33
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    real(RK)                       :: mue1, mue2, mue3
    type(TSiteMIEnm), pointer      :: pLJ126
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j, k
    integer                        :: nu

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength
    nu = this%Molecule%NUnit

    ! Loop over all units in Molecule
    do k = 1, nu
      ! Check number of rotation axes
      if( this%Molecule%Unit(k)%isElongated ) then
        ! Positions and quaternions of unit k in particle i
        PX = this%P0(np, 1, k)
        PY = this%P0(np, 2, k)
        PZ = this%P0(np, 3, k)
        q1 = this%Q0(np, 1, k)
        q2 = this%Q0(np, 2, k)
        q3 = this%Q0(np, 3, k)
        q4 = this%Q0(np, 4, k)

        ! Normalise quaternions
#if ARCH == 3
        qinv = rsqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#else
        qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#endif
        q1 = q1 * qinv
        q2 = q2 * qinv
        q3 = q3 * qinv
        q4 = q4 * qinv
        this%Q0(np, 1, k) = q1
        this%Q0(np, 2, k) = q2
        this%Q0(np, 3, k) = q3
        this%Q0(np, 4, k) = q4

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

        ! Loop over LJ126 sites in unit
        do j = 1, this%Molecule%Unit(k)%NMIEnm
          pLJ126 => this%Molecule%Unit(k)%SiteMIEnm(j)
          r1 = pLJ126%r(1) * BoxLengthInv
          r2 = pLJ126%r(2) * BoxLengthInv
          r3 = pLJ126%r(3) * BoxLengthInv
          pLJ126%RX(np) = PX + r1 * A11 + r2 * A21 + r3 * A31
          pLJ126%RY(np) = PY + r1 * A12 + r2 * A22 + r3 * A32
          pLJ126%RZ(np) = PZ + r1 * A13 + r2 * A23 + r3 * A33
        end do

        ! Loop over charge sites in molecule
        do j = 1, this%Molecule%Unit(k)%NCharge
          pCharge => this%Molecule%Unit(k)%SiteCharge(j)
          r1 = pCharge%r(1) * BoxLengthInv
          r2 = pCharge%r(2) * BoxLengthInv
          r3 = pCharge%r(3) * BoxLengthInv
          pCharge%RX(np) = PX + r1 * A11 + r2 * A21 + r3 * A31
          pCharge%RY(np) = PY + r1 * A12 + r2 * A22 + r3 * A32
          pCharge%RZ(np) = PZ + r1 * A13 + r2 * A23 + r3 * A33
        end do

        ! Loop over dipole sites in molecule
        do j = 1, this%Molecule%Unit(k)%NDipole
          pDipole => this%Molecule%Unit(k)%SiteDipole(j)
          r1 = pDipole%r(1) * BoxLengthInv
          r2 = pDipole%r(2) * BoxLengthInv
          r3 = pDipole%r(3) * BoxLengthInv
          or1 = pDipole%or(1)
          or2 = pDipole%or(2)
          or3 = pDipole%or(3)
          pDipole%RX(np) = PX + r1 * A11 + r2 * A21 + r3 * A31
          pDipole%RY(np) = PY + r1 * A12 + r2 * A22 + r3 * A32
          pDipole%RZ(np) = PZ + r1 * A13 + r2 * A23 + r3 * A33
          pDipole%OX(np) = or1 * A11 + or2 * A21 + or3 * A31
          pDipole%OY(np) = or1 * A12 + or2 * A22 + or3 * A32
          pDipole%OZ(np) = or1 * A13 + or2 * A23 + or3 * A33
        end do

        ! Loop over quadrupole sites in molecule
        do j = 1, this%Molecule%Unit(k)%NQuadrupole
          pQuadrupole => this%Molecule%Unit(k)%SiteQuadrupole(j)
          r1 = pQuadrupole%r(1) * BoxLengthInv
          r2 = pQuadrupole%r(2) * BoxLengthInv
          r3 = pQuadrupole%r(3) * BoxLengthInv
          or1 = pQuadrupole%or(1)
          or2 = pQuadrupole%or(2)
          or3 = pQuadrupole%or(3)
          pQuadrupole%RX(np) = PX + r1 * A11 + r2 * A21 + r3 * A31
          pQuadrupole%RY(np) = PY + r1 * A12 + r2 * A22 + r3 * A32
          pQuadrupole%RZ(np) = PZ + r1 * A13 + r2 * A23 + r3 * A33
          pQuadrupole%OX(np) = or1 * A11 + or2 * A21 + or3 * A31
          pQuadrupole%OY(np) = or1 * A12 + or2 * A22 + or3 * A32
          pQuadrupole%OZ(np) = or1 * A13 + or2 * A23 + or3 * A33
        end do

        if( CutoffMode .eq. CenterofMass ) then
          mue1 = this%Molecule%Unit(k)%Mue(1)
          mue2 = this%Molecule%Unit(k)%Mue(2)
          mue3 = this%Molecule%Unit(k)%Mue(3)
          this%MueX(np, k) = mue1 * A11 + mue2 * A21 + mue3 * A31
          this%MueY(np, k) = mue1 * A12 + mue2 * A22 + mue3 * A32
          this%MueZ(np, k) = mue1 * A13 + mue2 * A23 + mue3 * A33
        end if

      else !If unit is not elongated

        ! Loop over LJ126 sites in molecule
        do i = 1, this%Molecule%Unit(k)%NMIEnm
          pLJ126 => this%Molecule%Unit(k)%SiteMIEnm(i)
          pLJ126%RX(np) = this%P0(np, 1, k)
          pLJ126%RY(np) = this%P0(np, 2, k)
          pLJ126%RZ(np) = this%P0(np, 3, k)
        end do

        ! Loop over charge sites in molecule
        do i = 1, this%Molecule%Unit(k)%NCharge
          pCharge => this%Molecule%Unit(k)%SiteCharge(i)
          pCharge%RX(np) = this%P0(np, 1, k)
          pCharge%RY(np) = this%P0(np, 2, k)
          pCharge%RZ(np) = this%P0(np, 3, k)
        end do

      end if
    end do

  end subroutine TComponent_Unit2Atom1Root


!==============================================================!
!  Subroutine TComponent_Unit2Atom1 (per molecule)             !
!==============================================================!

  subroutine TComponent_Unit2Atom1Mol( this, np )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: np

    ! Declare local variables
    real(RK)                       :: BoxLengthInv
    type(vector)                   :: P, siteCoord
    real(RK)                       :: q1, q2, q3, q4, qinv
    real(RK)                       :: A11, A12, A13
    real(RK)                       :: A21, A22, A23
    real(RK)                       :: A31, A32, A33
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    real(RK)                       :: mue1, mue2, mue3
    type(TSiteMIEnm), pointer      :: pLJ126
    type(TSiteTT68), pointer      :: pTT68
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j, k
    integer                        :: nu

    ! Broadcast positions and orientations to all processes
#if MPI_VER > 0
    ! in MC simulations, we only communicate during common equilibration
    if ( SimulationType .ne. MonteCarlo .or. ((Equilibration .and. CommonEqui) )) then
      call MPI_Bcast( this%P0(np, :, :), this%Molecule%NUnit*3, MPI_RK, NRootProc, Communicator, ierror )
      if( this%Molecule%isElongated ) then
        call MPI_Bcast( this%Q0(np, :, :), this%Molecule%NUnit*4, MPI_RK, NRootProc, Communicator, ierror )
      end if
    end if
#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength
    nu = this%Molecule%NUnit

    ! Loop over all units in Molecule
    do k = 1, nu
      ! Check number of rotation axes
      if( this%Molecule%Unit(k)%isElongated ) then
        ! Positions and quaternions of unit k in particle i
        P = vector(this%P0(np, 1, k), this%P0(np, 2, k), this%P0(np, 3, k))
        q1 = this%Q0(np, 1, k)
        q2 = this%Q0(np, 2, k)
        q3 = this%Q0(np, 3, k)
        q4 = this%Q0(np, 4, k)

        ! Normalise quaternions
#if ARCH == 3
        qinv = rsqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#else
        qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#endif
        q1 = q1 * qinv
        q2 = q2 * qinv
        q3 = q3 * qinv
        q4 = q4 * qinv
        this%Q0(np, 1, k) = q1
        this%Q0(np, 2, k) = q2
        this%Q0(np, 3, k) = q3
        this%Q0(np, 4, k) = q4

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

        ! Loop over LJ126 sites in unit
        do j = 1, this%Molecule%Unit(k)%NMIEnm
          pLJ126 => this%Molecule%Unit(k)%SiteMIEnm(j)
          r1 = pLJ126%r(1) * BoxLengthInv
          r2 = pLJ126%r(2) * BoxLengthInv
          r3 = pLJ126%r(3) * BoxLengthInv
          pLJ126%RX(np) = P%x + r1 * A11 + r2 * A21 + r3 * A31
          pLJ126%RY(np) = P%y + r1 * A12 + r2 * A22 + r3 * A32
          pLJ126%RZ(np) = P%z + r1 * A13 + r2 * A23 + r3 * A33
        end do

        ! Loop over TT sites in unit
        do j = 1, this%Molecule%NTT68
          pTT68 => this%Molecule%SiteTT68(j)
          r1 = pTT68%r(1) * BoxLengthInv
          r2 = pTT68%r(2) * BoxLengthInv
          r3 = pTT68%r(3) * BoxLengthInv
          pTT68%RX(np) = P%x + r1 * A11 + r2 * A21 + r3 * A31
          pTT68%RY(np) = P%y + r1 * A12 + r2 * A22 + r3 * A32
          pTT68%RZ(np) = P%z + r1 * A13 + r2 * A23 + r3 * A33
        end do

        ! Loop over charge sites in molecule
        do j = 1, this%Molecule%Unit(k)%NCharge
          pCharge => this%Molecule%Unit(k)%SiteCharge(j)
          siteCoord = vector(pCharge%r(1), pCharge%r(2), pCharge%r(3))
          siteCoord = SCALE_VECTOR(siteCoord, BoxLengthInv)
          pCharge%RX(np) = P%x + siteCoord%x * A11 + siteCoord%y * A21 + siteCoord%z * A31
          pCharge%RY(np) = P%y + siteCoord%x * A12 + siteCoord%y * A22 + siteCoord%z * A32
          pCharge%RZ(np) = P%z + siteCoord%x * A13 + siteCoord%y * A23 + siteCoord%z * A33
        end do

        ! Loop over dipole sites in molecule
        do j = 1, this%Molecule%Unit(k)%NDipole
          pDipole => this%Molecule%Unit(k)%SiteDipole(j)
          r1 = pDipole%r(1) * BoxLengthInv
          r2 = pDipole%r(2) * BoxLengthInv
          r3 = pDipole%r(3) * BoxLengthInv
          or1 = pDipole%or(1)
          or2 = pDipole%or(2)
          or3 = pDipole%or(3)
          pDipole%RX(np) = P%x + r1 * A11 + r2 * A21 + r3 * A31
          pDipole%RY(np) = P%y + r1 * A12 + r2 * A22 + r3 * A32
          pDipole%RZ(np) = P%z + r1 * A13 + r2 * A23 + r3 * A33
          pDipole%OX(np) = or1 * A11 + or2 * A21 + or3 * A31
          pDipole%OY(np) = or1 * A12 + or2 * A22 + or3 * A32
          pDipole%OZ(np) = or1 * A13 + or2 * A23 + or3 * A33
        end do

        ! Loop over quadrupole sites in molecule
        do j = 1, this%Molecule%Unit(k)%NQuadrupole
          pQuadrupole => this%Molecule%Unit(k)%SiteQuadrupole(j)
          r1 = pQuadrupole%r(1) * BoxLengthInv
          r2 = pQuadrupole%r(2) * BoxLengthInv
          r3 = pQuadrupole%r(3) * BoxLengthInv
          or1 = pQuadrupole%or(1)
          or2 = pQuadrupole%or(2)
          or3 = pQuadrupole%or(3)
          pQuadrupole%RX(np) = P%x + r1 * A11 + r2 * A21 + r3 * A31
          pQuadrupole%RY(np) = P%y + r1 * A12 + r2 * A22 + r3 * A32
          pQuadrupole%RZ(np) = P%z + r1 * A13 + r2 * A23 + r3 * A33
          pQuadrupole%OX(np) = or1 * A11 + or2 * A21 + or3 * A31
          pQuadrupole%OY(np) = or1 * A12 + or2 * A22 + or3 * A32
          pQuadrupole%OZ(np) = or1 * A13 + or2 * A23 + or3 * A33
        end do

        if( CutoffMode .eq. CenterofMass ) then
          mue1 = this%Molecule%Unit(k)%Mue(1)
          mue2 = this%Molecule%Unit(k)%Mue(2)
          mue3 = this%Molecule%Unit(k)%Mue(3)
          this%MueX(np, k) = mue1 * A11 + mue2 * A21 + mue3 * A31
          this%MueY(np, k) = mue1 * A12 + mue2 * A22 + mue3 * A32
          this%MueZ(np, k) = mue1 * A13 + mue2 * A23 + mue3 * A33
        end if

      else !If unit is not elongated

        ! Loop over LJ126 sites in molecule
        do i = 1, this%Molecule%Unit(k)%NMIEnm
          pLJ126 => this%Molecule%Unit(k)%SiteMIEnm(i)
          pLJ126%RX(np) = this%P0(np, 1, k)
          pLJ126%RY(np) = this%P0(np, 2, k)
          pLJ126%RZ(np) = this%P0(np, 3, k)
        end do

        ! Loop over TT68 sites in molecule
        do i = 1, this%Molecule%NTT68
          pTT68 => this%Molecule%SiteTT68(i)
          pTT68%RX(np) = this%P0(np, 1, k)
          pTT68%RY(np) = this%P0(np, 2, k)
          pTT68%RZ(np) = this%P0(np, 3, k)
        end do

        ! Loop over charge sites in molecule
        do i = 1, this%Molecule%Unit(k)%NCharge
          pCharge => this%Molecule%Unit(k)%SiteCharge(i)
          pCharge%RX(np) = this%P0(np, 1, k)
          pCharge%RY(np) = this%P0(np, 2, k)
          pCharge%RZ(np) = this%P0(np, 3, k)
        end do

      end if
    end do


  end subroutine TComponent_Unit2Atom1Mol


!==============================================================!
!  Subroutine TComponent_Unit2Mol                              !
!==============================================================!

  subroutine TComponent_Unit2Mol( this )

    implicit none

    ! Declare arguments
    type(TComponent)   :: this

    ! Declare local variables
    real(RK)                 :: mass, avg, dist(this%Molecule%NUnit)
    real(RK)                 :: PX(this%NPart),PY(this%NPart),PZ(this%NPart)
    integer                  :: i, j
    integer                  :: np

    np = this%NPart
    mass = 0._RK
    PX(1:np)   = 0._RK
    PY(1:np)   = 0._RK
    PZ(1:np)   = 0._RK

    do i=1,this%Molecule%NUnit
      mass = this%Molecule%Unit(i)%Mass
      PX(1:np)   = PX(1:np)   + (this%P0(1:np,1,i)-&
&          anint(this%P0(1:np,1,i)-this%Pm0(1:np,1)) )*mass
      PY(1:np)   = PY(1:np)   + (this%P0(1:np,2,i)-&
&          anint(this%P0(1:np,2,i)-this%Pm0(1:np,2)) )*mass
      PZ(1:np)   = PZ(1:np)   + (this%P0(1:np,3,i)-&
&          anint(this%P0(1:np,3,i)-this%Pm0(1:np,3)) )*mass
    end do

    mass = this%Molecule%Mass
    this%Pm0(1:np,1) = PX / mass
    this%Pm0(1:np,2) = PY / mass
    this%Pm0(1:np,3) = PZ / mass

  end subroutine TComponent_Unit2Mol


!==============================================================!
!  Subroutine TComponent_Unit2Mol1                             !
!==============================================================!

  subroutine TComponent_Unit2Mol1( this, np )

    implicit none

    ! Declare arguments
    type(TComponent)   :: this
    integer,intent (in) :: np

    ! Declare local variables
    real(RK)                 :: mass, avg
    real(RK)                 :: PX,PY,PZ
    integer                  :: i

    mass = 0._RK
    PX   = 0._RK
    PY   = 0._RK
    PZ   = 0._RK

    do i=1,this%Molecule%NUnit
      mass = this%Molecule%Unit(i)%Mass
      PX   = PX   + (this%P0(np,1,i)-&
&          anint(this%P0(np,1,i)-this%Pm0(np,1)) )*mass
      PY   = PY   + (this%P0(np,2,i) - &
&          anint(this%P0(np,2,i)-this%Pm0(np,2)) )*mass
      PZ   = PZ   + ( this%P0(np,3,i) - &
&          anint(this%P0(np,3,i)-this%Pm0(np,3) ) )*mass
    end do

    mass = this%Molecule%Mass
    this%Pm0(np,1) = PX / mass
    this%Pm0(np,2) = PY / mass
    this%Pm0(np,3) = PZ / mass

  end subroutine TComponent_Unit2Mol1


end module ms2_component
