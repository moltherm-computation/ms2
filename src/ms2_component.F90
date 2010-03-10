!==============================================================!
!  MOLECULAR SIMULATION PROGRAM MS2 Version 1.1 v12            !
!  (c) 2001 by Sergey Lishchuk, ITT                            !
!  (c) 2007 by Bernhard Eckl, ITT                              !
!==============================================================!
!  Module ms2_component                                        !
!  Contains TComponent object                                  !
!==============================================================!

#ifndef ARCH
#define ARCH    0
#define FORTRAN 90
#define MPI_VER 0
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

    ! Positions and orientations of test particles
    real(RK), pointer :: Pm0Test(:, :), Qm0Test(:, :)

    ! Positions and orientations for units of test particles
    real(RK), pointer :: P0Test(:, :, :), Q0Test(:, :, :)

    ! Centers of mass positions and their derivatives for molecules
    real(RK), pointer :: Pm0(:, :)
    real(RK), pointer :: P0Save( :, :, :)
    real(RK), pointer :: Pm0old(:, :)


    ! Centers of mass positions and their derivatives for Units
    real(RK), pointer :: P0(:, :, :)
    real(RK), pointer :: P1(:, :, :)
    real(RK), pointer :: P2(:, :, :)
    real(RK), pointer :: P3(:, :, :)
    real(RK), pointer :: P4(:, :, :)
    real(RK), pointer :: P5(:, :, :)


    ! Quaternion parameters for molecules - only to calculate the initial orientation
    real(RK), pointer :: Qm0(:, :)



    ! Quaternion parameters and their derivatives for Units
    real(RK), pointer :: Q0(:, :, :)
     real(RK), pointer :: Q0Save(:, :, :)
    real(RK), pointer :: Q0tmp(:, :, :)
    real(RK), pointer :: Q1(:, :, :)
    real(RK), pointer :: Q2(:, :, :)
    real(RK), pointer :: Q3(:, :, :)
    real(RK), pointer :: Q4(:, :, :)


    ! Angular velocities and their derivatives for units
    real(RK), pointer :: W0(:, :, :)
    real(RK), pointer :: W1(:, :, :)
    real(RK), pointer :: W2(:, :, :)
    real(RK), pointer :: W3(:, :, :)
    real(RK), pointer :: W4(:, :, :)

    ! Displacement for molecules
    real(RK), pointer :: Disp(:, :)

    ! Total forces acting on units
    real(RK), pointer :: F(:, :, :)
#if MPI_VER > 0
    real(RK), pointer :: FAll(:, :, :)
#endif

    ! Total torques acting on units
    real(RK), pointer :: T(:, :, :)
#if MPI_VER > 0
    real(RK), pointer :: TAll(:, :, :)
#endif

    ! Total dipole moment of units of molecule for reaction field
    real(RK), pointer :: MueX(:, :), MueY(:, :), MueZ(:, :)

    ! Torques from reaction field, space fixed
    real(RK), pointer :: tRFX(:, :), tRFY(:, :), tRFZ(:, :)

    ! Total dipole moment of test particles for reaction field
    real(RK), pointer :: MueXTest(:, :), MueYTest(:, :), MueZTest(:, :)

    ! Gear corrector local arrays
    real(RK), pointer :: Corr0(:, :, :)
    real(RK), pointer :: Corr1(:, :, :)

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
    integer  :: ChemPotMethod, WFMethod
    integer  :: FluctState
    real(RK) :: ChemPot, WidomContribution
!DEBUG
    real(RK) :: ChemPot1, ChemPot2
!DEBUG
    real(RK) :: ChemPot0, PartialMolarVolume
    real(RK) :: VarChemPot, VarPartialMolarVolume

    ! IDF
    integer, pointer :: UnitLJ(:),UnitC(:),UnitDP(:),UnitQP(:)

    ! Ewald
    real(RK) :: EPotTestSelf

    ! Fluctuating components and weighting factors
    integer           :: NFluctState, NFluctMax
    integer, pointer  :: NState(:), NStateWF(:), NFluctComp(:)
    real(RK), pointer :: WF(:)
    real(RK)          :: ProbW0, ProbW1, ProbW0V, ProbW1Rho
!DEBUG
    integer, pointer  :: NFluctUpAttempts(:), NFluctUpSuccesses(:)
    integer, pointer  :: NFluctDownAttempts(:), NFluctDownSuccesses(:)
!     integer, pointer  :: NStateBF(:)
!     real(RK), pointer :: BFSumState(:)
!DEBUG

    ! Mole fraction in corresponding liquid simulation (for GE ensemble only)
    real(RK) :: LiqFraction

    ! Long-range corrections
    real(RK) :: EPotTestCorrLJ
    real(RK) :: EPotTestCorrRF

    ! Inner Degrees of Freedom
    integer,pointer :: BondCount(:)
    integer,pointer :: BoPartner(:,:)
    integer,pointer :: AngleCount(:)
    integer,pointer :: AnglePartner(:,:)
    integer,pointer :: DihedralCount(:)
    integer,pointer :: DihedralPartner(:,:)


    ! Accumulated sums, averages and errors
    type(TAccumulator) :: SumInvChemPotRho
    type(TAccumulator) :: SumInvChemPot
!DEBUG
    type(TAccumulator) :: SumInvChemPotRho1
    type(TAccumulator) :: SumInvChemPot1
    type(TAccumulator) :: SumInvChemPotRho2
    type(TAccumulator) :: SumInvChemPot2
!DEBUG
    type(TAccumulator) :: SumChemPotV
    type(TAccumulator) :: SumChemPotVV
    type(TAccumulator) :: SumVW
!DEBUG
    type(TAccumulator) :: SumVW1
    type(TAccumulator) :: SumVW2
!DEBUG
    type(TAccumulator) :: SumFraction

    ! Potential model for this component
    type(TMolecule) :: Molecule

    ! File name for potential model
    character(FileNameLength) :: PotModFileName

  end type TComponent

  interface Construct
    module procedure TComponent_Construct
    module procedure TComponent_ConstructSVC
    module procedure TComponent_ConstructFluct
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

  !interface Mol2Atom
  !  module procedure TComponent_Mol2Atom
  !end interface

  !interface Mol2Atom1
  !  module procedure TComponent_Mol2Atom1
  !end interface

  interface Mol2AtomTest
    module procedure TComponent_Mol2AtomTest
  end interface

!  interface Atom2Mol
!    module procedure TComponent_Atom2Mol
!  end interface

  interface Mol2Resize
    module procedure TComponent_Mol2Resize
  end interface

  interface Mol2Unit
    module procedure TComponent_Mol2Unit
    module procedure TComponent_Mol2UnitRotate
end interface

  interface Mol2Unit1
    module procedure TComponent_Mol2Unit1
  end interface

  interface Mol2Unit1Test
    module procedure TComponent_Mol2Unit1Test
  end interface

  interface Unit2Atom
    module procedure TComponent_Unit2Atom
  end interface

  interface Unit2Atom1
    module procedure TComponent_Unit2Atom1Mol
    module procedure TComponent_Unit2Atom1
  end interface

  interface Unit2AtomTest
    module procedure TComponent_Unit2Atom1Test
  end interface

  interface Atom2Unit
    module procedure TComponent_Atom2Unit
  end interface

  interface Unit2Mol
    module procedure TComponent_Unit2Mol
    module procedure TComponent_Unit2Mol1
  end interface
  
  interface Flex2Rigid
    module procedure TComponent_Flex2Rigid
  end interface
  
  interface Rigid2Flex
    module procedure TComponent_Rigid2Flex
  end interface

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
    write( IOBuffer, '("-----------------------------------------------------------")')
    call LogWrite
    write( IOBuffer, '(T13, "Reading components for ensemble", I3)') comp
    call LogWrite
    call FileReadParameter( this%Fraction, iounit_params , IdFraction, .false. )
    write( IOBuffer, '("Mole fraction of component ", A, ": ", F9.6)' ) &
&     trim( this%PotModFileName ), this%Fraction
    call LogWrite

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
      write( IOBuffer, &
&       '("Reduced ChemPot0 of component ", A, ": ", F9.6, " (", F9.6, ")")' ) &
&       trim( this%PotModFileName ), this%ChemPot0, this%VarChemPot
      call LogWrite
      write( IOBuffer, &
&       '("Reduced PartMolVol of component ", A, ": ", F9.6, " (", F9.6, ")")' ) &
&       trim( this%PotModFilename ), this%PartialMolarVolume, &
&       this%VarPartialMolarVolume
      call LogWrite

    else if( EnsembleType .eq. EnsembleTypeHA ) then
      if( comp == 1 ) then
        ! Read chemical potential of phase changing component (first one)
        call FileReadParameter( this%ChemPot, iounit_params , IdChemPot, .false. )
        call FileReadParameter( this%VarChemPot, iounit_params , IdVarChemPot, .false. )
        write( IOBuffer, &
&         '("Reduced ChemPot of component ", A, ": ", F9.6, " (", F9.6, ")")' ) &
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
      case default
        call Error( trim( str )// &
&         ' method for calculation of chemical potential is not implemented' )
      end select
      if( this%ChemPotMethod .eq. ChemPotMethodGradIns .and. &
&         .not. SimulationType .eq. MonteCarlo ) &
&       call Error( 'Gradual insertion is only allowed for MonteCarlo simulation' )
      write( IOBuffer, &
&       '("Chemical potential of ", A, " will be calculated by: ", A)' ) &
&       trim( this%PotModFilename )
      call LogWrite
      write( IOBuffer, '(T10, "-> ", A)' ) trim( str )
      call LogWrite

      ! Read number of test particles
      if( this%ChemPotMethod .eq. ChemPotMethodWidom ) then
        call FileReadParameter( this%NTest, iounit_params, IdNTest, .false. )
        if( this%NTest <= 0 ) &
&         call Error( 'Number of test particles need to be > 0' )
        write( IOBuffer, '(T10, "-> Number of test particles:", I11 )' ) this%NTest
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
          call Error( trim( str )// &
&           ' method for weighting factors is not implemented' )
        end select
        write( IOBuffer, '("Estimation of weighting factors: using ", A )' ) &
&         trim( str )
        call LogWrite
      end if

    end if

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
      nullify ( this%DispTran )
      nullify ( this%DispRot )
      allocate( this%DispTran, STAT = stat )
      call AllocationError( stat, 'maximum MC displacement' )
      allocate( this%DispRot, STAT = stat )
      call AllocationError( stat, 'maximum MC displacement' )
      nullify ( this%DispMolTran )
      nullify ( this%DispMolRot )
      allocate( this%DispMolTran, STAT = stat )
      call AllocationError( stat, 'maximum MC molecule displacement' )
      allocate( this%DispMolRot, STAT = stat )
      call AllocationError( stat, 'maximum MC molecule displacement' )
    end if

    ! Allocate and read weighting factors
    if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
      nullify( this%WF )
      allocate( this%WF( 0:this%NFluctMax ), STAT = stat )
      call AllocationError( stat, 'fluctuating particle states', &
&       this%NFluctMax + 1 )
      if( this%WFMethod .eq. WFMethodGuess .or. &
&         this%WFMethod .eq. WFMethodOptSet ) then
        if( RootProc ) read( iounit_params, * ) this%WF
#if MPI_VER > 0
        call MPI_Bcast( this%WF, size( this%WF ), MPI_DOUBLE_PRECISION, &
&         NRootProc, MPI_COMM_WORLD, ierror )
#endif
      end if
    end if

    ! Allocate fluctuating particle components vector
    nullify( this%NFluctComp )
    if( this%NFluctMax > 0 ) then
      allocate( this%NFluctComp( 0:this%NFluctMax ), STAT = stat )
      call AllocationError( stat, 'fluctuating particle components', &
&       this%NFluctMax + 1 )
    end if

    write( IOBuffer, '(T8, "Reading components for ensemble", I3, " successful")') comp
    call LogWrite
    write( IOBuffer, '("-----------------------------------------------------------")')
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
    nullify( this%NPart )
    nullify( this%NPart1 )
    nullify( this%NPart2 )
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

    ! Deallocation of MC vectors - even in MD
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
!  Subroutine TComponent_Destruct                              !
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

    ! Construct accumulators
    select case( this%ChemPotMethod )
    case( ChemPotMethodGradIns )
      call Construct( this%SumInvChemPotRho, .true. )
      call Construct( this%SumInvChemPot, .true. )
!DEBUG
      call Construct( this%SumInvChemPotRho1, .true. )
      call Construct( this%SumInvChemPot1, .true. )
      call Construct( this%SumInvChemPotRho2, .true. )
      call Construct( this%SumInvChemPot2, .true. )
!DEBUG
      call Construct( this%SumVW, .true. )
!DEBUG
      call Construct( this%SumVW1, .true. )
      call Construct( this%SumVW2, .true. )
!DEBUG
    case( ChemPotMethodWidom )
      call Construct( this%SumChemPotV, .false. )
      call Construct( this%SumChemPotVV, .false. )
      call Construct( this%SumVW, .true. )
    end select

    if( EnsembleType .eq. EnsembleTypeGE .or. &
&       EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
      call Construct( this%SumFraction, .false. )
    end if

  end subroutine TComponent_CreateAccumulators



!==============================================================!
!  Subroutine TComponent_DestroyAccumulators                   !
!==============================================================!

  subroutine TComponent_DestroyAccumulators( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Destruct accumulators
    select case( this%ChemPotMethod )
    case( ChemPotMethodGradIns )
      call Destruct( this%SumInvChemPotRho )
      call Destruct( this%SumInvChemPot )
!DEBUG
      call Destruct( this%SumInvChemPotRho1 )
      call Destruct( this%SumInvChemPot1 )
      call Destruct( this%SumInvChemPotRho2 )
      call Destruct( this%SumInvChemPot2 )
!DEBUG
      call Destruct( this%SumVW )
!DEBUG
      call Destruct( this%SumVW1 )
      call Destruct( this%SumVW2 )
!DEBUG
    case( ChemPotMethodWidom )
      call Destruct( this%SumChemPotV )
      call Destruct( this%SumChemPotVV )
      call Destruct( this%SumVW )
    end select

    if( EnsembleType .eq. EnsembleTypeGE .or. &
&       EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
      call Destruct( this%SumFraction )
    end if

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
    integer :: nu, nup, neu, neup
    integer :: i
    integer :: j, k, l, index
    integer :: stat
    logical :: Site1, Site2, Site3, Site4
    integer :: SiteId1, SiteId2, SiteId3, SiteId4
    logical :: ok


    ! Set maximum number of particles and number of test particles
    np = this%NPartMax
    nu = this%Molecule%NUnit
    neu = this%Molecule%NEUnit ! number of elongated Units in molecule
    nup = nu*np
    neup = neu*np
    ntest = this%NTest

    ! Nullify pointers
    nullify( this%Pm0 )
    nullify( this%P0Save )
    nullify( this%Pm0old )
    nullify( this%P0 )
    nullify( this%P1 )
    nullify( this%P2 )
    nullify( this%P3 )
    nullify( this%P4 )
    nullify( this%P5 )
    nullify( this%Disp )
    nullify( this%F )
#if MPI_VER > 0
    nullify( this%FAll )
#endif
    nullify( this%Qm0 )
    nullify( this%Q0Save )
    nullify( this%Q0tmp )
    nullify( this%Q0 )
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

    ! Centers of mass positions
    allocate( this%Pm0( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%P0Save( np, 3, nu ), STAT = stat )
    call AllocationError( stat, 'units*particles', nup )
    allocate( this%Pm0old( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

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

      ! Total forces
      allocate( this%F( np, 3, nu ), STAT = stat )
      call AllocationError( stat, 'units*particles', nup )
#if MPI_VER > 0
      allocate( this%FAll( np, 3, nu ), STAT = stat )
      call AllocationError( stat, 'units*particles', nup )
#endif

    end if

    if( this%Molecule%isElongated ) then

      ! Quaternion parameters
      allocate( this%Qm0( np, 4 ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%Q0Save( np, 4, nu ), STAT = stat )
      call AllocationError( stat, 'units*particles', nup )

      ! Quaternion parameters for Units
      allocate( this%Q0( np, 4, nu ), STAT = stat )
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
    if( SimulationType .eq. MolecularDynamics &
&     .and. IntegratorType .eq. IntegratorTypeGear ) then
      allocate( this%Corr0( np, merge( 4, 3, this%Molecule%isElongated ),nu ), &
&               STAT = stat )
      call AllocationError( stat, 'units*particles', nup )
      allocate( this%Corr1( np, merge( 4, 3, this%Molecule%isElongated ),nu ), &
&               STAT = stat )
      call AllocationError( stat, 'units*particles', nup )
    end if

    ! Site positions, orientations, forces and torques
    do i = 1, this%Molecule%NLJ126
      this%Molecule%SiteLJ126(i)%NPartMax => this%NPartMax
      this%Molecule%SiteLJ126(i)%NPart => this%NPart
      this%Molecule%SiteLJ126(i)%NTest => this%NTest
      this%Molecule%SiteLJ126(i)%NPart0 => this%NPart0
      this%Molecule%SiteLJ126(i)%NPart1 => this%NPart1
      this%Molecule%SiteLJ126(i)%NPart2 => this%NPart2
      call Allocate( this%Molecule%SiteLJ126(i) )
!      this%Molecule%SiteLJ126(i)%PX => this%P0(:, 1)
!      this%Molecule%SiteLJ126(i)%PY => this%P0(:, 2)
!      this%Molecule%SiteLJ126(i)%PZ => this%P0(:, 3)
      if( ntest > 0 ) then
       this%Molecule%SiteLJ126(i)%PXTest => this%Pm0Test(:, 1)
       this%Molecule%SiteLJ126(i)%PYTest => this%Pm0Test(:, 2)
       this%Molecule%SiteLJ126(i)%PZTest => this%Pm0Test(:, 3)
      end if
    end do
    do i = 1, this%Molecule%NCharge
      this%Molecule%SiteCharge(i)%NPartMax => this%NPartMax
      this%Molecule%SiteCharge(i)%NPart => this%NPart
      this%Molecule%SiteCharge(i)%NTest => this%NTest
      this%Molecule%SiteCharge(i)%NPart0 => this%NPart0
      this%Molecule%SiteCharge(i)%NPart1 => this%NPart1
      this%Molecule%SiteCharge(i)%NPart2 => this%NPart2
      call Allocate( this%Molecule%SiteCharge(i) )
!      this%Molecule%SiteCharge(i)%PX => this%P0(:, 1)
!      this%Molecule%SiteCharge(i)%PY => this%P0(:, 2)
!      this%Molecule%SiteCharge(i)%PZ => this%P0(:, 3)
      if( ntest > 0 ) then
       this%Molecule%SiteCharge(i)%PXTest => this%Pm0Test(:, 1)
       this%Molecule%SiteCharge(i)%PYTest => this%Pm0Test(:, 2)
       this%Molecule%SiteCharge(i)%PZTest => this%Pm0Test(:, 3)
      end if
    end do
    do i = 1, this%Molecule%NDipole
      this%Molecule%SiteDipole(i)%NPartMax => this%NPartMax
      this%Molecule%SiteDipole(i)%NPart => this%NPart
      this%Molecule%SiteDipole(i)%NTest => this%NTest
      this%Molecule%SiteDipole(i)%NPart0 => this%NPart0
      this%Molecule%SiteDipole(i)%NPart1 => this%NPart1
      this%Molecule%SiteDipole(i)%NPart2 => this%NPart2
      call Allocate( this%Molecule%SiteDipole(i) )
!      this%Molecule%SiteDipole(i)%PX => this%P0(:, 1)
!      this%Molecule%SiteDipole(i)%PY => this%P0(:, 2)
!      this%Molecule%SiteDipole(i)%PZ => this%P0(:, 3)
      if( ntest > 0 ) then
       this%Molecule%SiteDipole(i)%PXTest => this%Pm0Test(:, 1)
       this%Molecule%SiteDipole(i)%PYTest => this%Pm0Test(:, 2)
       this%Molecule%SiteDipole(i)%PZTest => this%Pm0Test(:, 3)
      end if
    end do
    do i = 1, this%Molecule%NQuadrupole
      this%Molecule%SiteQuadrupole(i)%NPartMax => this%NPartMax
      this%Molecule%SiteQuadrupole(i)%NPart => this%NPart
      this%Molecule%SiteQuadrupole(i)%NTest => this%NTest
      this%Molecule%SiteQuadrupole(i)%NPart0 => this%NPart0
      this%Molecule%SiteQuadrupole(i)%NPart1 => this%NPart1
      this%Molecule%SiteQuadrupole(i)%NPart2 => this%NPart2
      call Allocate( this%Molecule%SiteQuadrupole(i) )
!      this%Molecule%SiteQuadrupole(i)%PX => this%P0(:, 1)
!      this%Molecule%SiteQuadrupole(i)%PY => this%P0(:, 2)
!      this%Molecule%SiteQuadrupole(i)%PZ => this%P0(:, 3)
      if( ntest > 0 ) then
       this%Molecule%SiteQuadrupole(i)%PXTest => this%Pm0Test(:, 1)
       this%Molecule%SiteQuadrupole(i)%PYTest => this%Pm0Test(:, 2)
       this%Molecule%SiteQuadrupole(i)%PZTest => this%Pm0Test(:, 3)
      end if
    end do

    ! Inner degrees of freedom

    ! Units
    do i = 1, this%Molecule%NUnit
      this%Molecule%Unit(i)%NPartMax => this%NPartMax
      this%Molecule%Unit(i)%NPart => this%NPart
!      this%Molecule%Unit(i)%NTest => this%NTest
      this%Molecule%Unit(i)%NPart0 => this%NPart0
      this%Molecule%Unit(i)%NPart1 => this%NPart1
      this%Molecule%Unit(i)%NPart2 => this%NPart2
      this%Molecule%Unit(i)%PX => this%P0(:, 1, i)
      this%Molecule%Unit(i)%PY => this%P0(:, 2, i)
      this%Molecule%Unit(i)%PZ => this%P0(:, 3, i)
      if (this%Molecule%Unit(i)%NLJ126 > 0) then
        do j = 1, this%Molecule%Unit(i)%NLJ126
          if (UseIntDegFreed) then
              call binar_search(this%Molecule%SiteLJ126%SiteId, this%Molecule%Unit(i)%SiteLJ126(j)%SiteId, ok, index )
              if (ok) then
                  this%Molecule%Unit(i)%SiteLJ126(j)%r=>this%Molecule%SiteLJ126(index)%r
                  this%Molecule%Unit(i)%SiteLJ126(j)%RX=>this%Molecule%SiteLJ126(index)%RX
                  this%Molecule%Unit(i)%SiteLJ126(j)%RY=>this%Molecule%SiteLJ126(index)%RY
                  this%Molecule%Unit(i)%SiteLJ126(j)%RZ=>this%Molecule%SiteLJ126(index)%RZ
                  this%Molecule%Unit(i)%SiteLJ126(j)%FX=>this%Molecule%SiteLJ126(index)%FX
                  this%Molecule%Unit(i)%SiteLJ126(j)%FY=>this%Molecule%SiteLJ126(index)%FY
                  this%Molecule%Unit(i)%SiteLJ126(j)%FZ=>this%Molecule%SiteLJ126(index)%FZ
                  this%Molecule%SiteLJ126(index)%PX => this%Molecule%Unit(i)%PX
                  this%Molecule%SiteLJ126(index)%PY => this%Molecule%Unit(i)%PY
                  this%Molecule%SiteLJ126(index)%PZ => this%Molecule%Unit(i)%PZ
              end if
           else
                 this%Molecule%Unit(i)%SiteLJ126(j)%r=>this%Molecule%SiteLJ126(j)%r
                 this%Molecule%Unit(i)%SiteLJ126(j)%RX=>this%Molecule%SiteLJ126(j)%RX
                 this%Molecule%Unit(i)%SiteLJ126(j)%RY=>this%Molecule%SiteLJ126(j)%RY
                 this%Molecule%Unit(i)%SiteLJ126(j)%RZ=>this%Molecule%SiteLJ126(j)%RZ
                 this%Molecule%Unit(i)%SiteLJ126(j)%FX=>this%Molecule%SiteLJ126(j)%FX
                 this%Molecule%Unit(i)%SiteLJ126(j)%FY=>this%Molecule%SiteLJ126(j)%FY
                 this%Molecule%Unit(i)%SiteLJ126(j)%FZ=>this%Molecule%SiteLJ126(j)%FZ
                 this%Molecule%SiteLJ126(j)%PX =>this%Molecule%Unit(i)%PX
                 this%Molecule%SiteLJ126(j)%PY =>this%Molecule%Unit(i)%PY
                 this%Molecule%SiteLJ126(j)%PZ =>this%Molecule%Unit(i)%PZ
            end if
         end do
       end if
       if (this%Molecule%Unit(i)%NCharge > 0) then
         do j = 1, this%Molecule%Unit(i)%NCharge
            if (UseIntDegFreed) then
                call binar_search(this%Molecule%SiteCharge%SiteId, this%Molecule%Unit(i)%SiteCharge(j)%SiteId, ok, index )
                if (ok) then
                    this%Molecule%Unit(i)%SiteCharge(j)%r=>this%Molecule%SiteCharge(index)%r
                    this%Molecule%Unit(i)%SiteCharge(j)%RX=>this%Molecule%SiteCharge(index)%RX
                    this%Molecule%Unit(i)%SiteCharge(j)%RY=>this%Molecule%SiteCharge(index)%RY
                    this%Molecule%Unit(i)%SiteCharge(j)%RZ=>this%Molecule%SiteCharge(index)%RZ
                    this%Molecule%Unit(i)%SiteCharge(j)%FX=>this%Molecule%SiteCharge(index)%FX
                    this%Molecule%Unit(i)%SiteCharge(j)%FY=>this%Molecule%SiteCharge(index)%FY
                    this%Molecule%Unit(i)%SiteCharge(j)%FZ=>this%Molecule%SiteCharge(index)%FZ
                    this%Molecule%SiteCharge(index)%PX =>this%Molecule%Unit(i)%PX
                    this%Molecule%SiteCharge(index)%PY =>this%Molecule%Unit(i)%PY
                    this%Molecule%SiteCharge(index)%PZ =>this%Molecule%Unit(i)%PZ
                end if
            else
                this%Molecule%Unit(i)%SiteCharge(j)%r=>this%Molecule%SiteCharge(j)%r
                this%Molecule%Unit(i)%SiteCharge(j)%RX=>this%Molecule%SiteCharge(j)%RX
                this%Molecule%Unit(i)%SiteCharge(j)%RY=>this%Molecule%SiteCharge(j)%RY
                this%Molecule%Unit(i)%SiteCharge(j)%RZ=>this%Molecule%SiteCharge(j)%RZ
                this%Molecule%Unit(i)%SiteCharge(j)%FX=>this%Molecule%SiteCharge(j)%FX
                this%Molecule%Unit(i)%SiteCharge(j)%FY=>this%Molecule%SiteCharge(j)%FY
                this%Molecule%Unit(i)%SiteCharge(j)%FZ=>this%Molecule%SiteCharge(j)%FZ
                this%Molecule%SiteCharge(j)%PX => this%Molecule%Unit(i)%PX
                this%Molecule%SiteCharge(j)%PY => this%Molecule%Unit(i)%PY
                this%Molecule%SiteCharge(j)%PZ => this%Molecule%Unit(i)%PZ
          end if
         end do
       end if
       if (this%Molecule%Unit(i)%NDipole > 0) then
         do j = 1, this%Molecule%Unit(i)%NDipole
            if (UseIntDegFreed) then
                call binar_search(this%Molecule%SiteDipole%SiteId, this%Molecule%Unit(i)%SiteDipole(j)%SiteId, ok, index )
                if (ok) then
                    this%Molecule%Unit(i)%SiteDipole(j)%r=>this%Molecule%SiteDipole(index)%r
                    this%Molecule%Unit(i)%SiteDipole(j)%or=>this%Molecule%SiteDipole(index)%or
                    this%Molecule%Unit(i)%SiteDipole(j)%RX=>this%Molecule%SiteDipole(index)%RX
                    this%Molecule%Unit(i)%SiteDipole(j)%RY=>this%Molecule%SiteDipole(index)%RY
                    this%Molecule%Unit(i)%SiteDipole(j)%RZ=>this%Molecule%SiteDipole(index)%RZ
                    this%Molecule%Unit(i)%SiteDipole(j)%OX=>this%Molecule%SiteDipole(index)%OX
                    this%Molecule%Unit(i)%SiteDipole(j)%OY=>this%Molecule%SiteDipole(index)%OY
                    this%Molecule%Unit(i)%SiteDipole(j)%OZ=>this%Molecule%SiteDipole(index)%OZ
                    this%Molecule%Unit(i)%SiteDipole(j)%FX=>this%Molecule%SiteDipole(index)%FX
                    this%Molecule%Unit(i)%SiteDipole(j)%FY=>this%Molecule%SiteDipole(index)%FY
                    this%Molecule%Unit(i)%SiteDipole(j)%FZ=>this%Molecule%SiteDipole(index)%FZ
                    this%Molecule%Unit(i)%SiteDipole(j)%TX=>this%Molecule%SiteDipole(index)%TX
                    this%Molecule%Unit(i)%SiteDipole(j)%TY=>this%Molecule%SiteDipole(index)%TY
                    this%Molecule%Unit(i)%SiteDipole(j)%TZ=>this%Molecule%SiteDipole(index)%TZ
                    this%Molecule%SiteDipole(index)%PX => this%Molecule%Unit(i)%PX
                    this%Molecule%SiteDipole(index)%PY => this%Molecule%Unit(i)%PY
                    this%Molecule%SiteDipole(index)%PZ => this%Molecule%Unit(i)%PZ
                end if
            else
                this%Molecule%Unit(i)%SiteDipole(j)%r=>this%Molecule%SiteDipole(j)%r
                this%Molecule%Unit(i)%SiteDipole(j)%or=>this%Molecule%SiteDipole(j)%or
                this%Molecule%Unit(i)%SiteDipole(j)%RX=>this%Molecule%SiteDipole(j)%RX
                this%Molecule%Unit(i)%SiteDipole(j)%RY=>this%Molecule%SiteDipole(j)%RY
                this%Molecule%Unit(i)%SiteDipole(j)%RZ=>this%Molecule%SiteDipole(j)%RZ
                this%Molecule%Unit(i)%SiteDipole(j)%OX=>this%Molecule%SiteDipole(j)%OX
                this%Molecule%Unit(i)%SiteDipole(j)%OY=>this%Molecule%SiteDipole(j)%OY
                this%Molecule%Unit(i)%SiteDipole(j)%OZ=>this%Molecule%SiteDipole(j)%OZ
                this%Molecule%Unit(i)%SiteDipole(j)%FX=>this%Molecule%SiteDipole(j)%FX
                this%Molecule%Unit(i)%SiteDipole(j)%FY=>this%Molecule%SiteDipole(j)%FY
                this%Molecule%Unit(i)%SiteDipole(j)%FZ=>this%Molecule%SiteDipole(j)%FZ
                this%Molecule%Unit(i)%SiteDipole(j)%TX=>this%Molecule%SiteDipole(j)%TX
                this%Molecule%Unit(i)%SiteDipole(j)%TY=>this%Molecule%SiteDipole(j)%TY
                this%Molecule%Unit(i)%SiteDipole(j)%TZ=>this%Molecule%SiteDipole(j)%TZ
                this%Molecule%SiteDipole(j)%PX=> this%Molecule%Unit(i)%PX
                this%Molecule%SiteDipole(j)%PY=> this%Molecule%Unit(i)%PY
                this%Molecule%SiteDipole(j)%PZ=> this%Molecule%Unit(i)%PZ
            end if
         end do
       end if
       if (this%Molecule%Unit(i)%NQuadrupole > 0) then
         do j = 1, this%Molecule%Unit(i)%NQuadrupole
            if (UseIntDegFreed) then
                call binar_search(this%Molecule%SiteQuadrupole%SiteId, this%Molecule%Unit(i)%SiteQuadrupole(j)%SiteId, ok, index )
                if (ok) then
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%r=>this%Molecule%SiteQuadrupole(index)%r
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%or=>this%Molecule%SiteQuadrupole(index)%or
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%RX=>this%Molecule%SiteQuadrupole(index)%RX
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%RY=>this%Molecule%SiteQuadrupole(index)%RY
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%RZ=>this%Molecule%SiteQuadrupole(index)%RZ
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%OX=>this%Molecule%SiteQuadrupole(index)%OX
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%OY=>this%Molecule%SiteQuadrupole(index)%OY
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%OZ=>this%Molecule%SiteQuadrupole(index)%OZ
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%FX=>this%Molecule%SiteQuadrupole(index)%FX
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%FY=>this%Molecule%SiteQuadrupole(index)%FY
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%FZ=>this%Molecule%SiteQuadrupole(index)%FZ
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%TX=>this%Molecule%SiteQuadrupole(index)%TX
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%TY=>this%Molecule%SiteQuadrupole(index)%TY
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%TZ=>this%Molecule%SiteQuadrupole(index)%TZ
                    this%Molecule%SiteQuadrupole(index)%PX=>this%Molecule%Unit(i)%PX
                    this%Molecule%SiteQuadrupole(index)%PY=>this%Molecule%Unit(i)%PY
                    this%Molecule%SiteQuadrupole(index)%PZ=>this%Molecule%Unit(i)%PZ
                end if
            else
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%r=>this%Molecule%SiteQuadrupole(j)%r
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%or=>this%Molecule%SiteQuadrupole(j)%or
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%RX=>this%Molecule%SiteQuadrupole(j)%RX
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%RY=>this%Molecule%SiteQuadrupole(j)%RY
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%RZ=>this%Molecule%SiteQuadrupole(j)%RZ
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%OX=>this%Molecule%SiteQuadrupole(j)%OX
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%OY=>this%Molecule%SiteQuadrupole(j)%OY
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%OZ=>this%Molecule%SiteQuadrupole(j)%OZ
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%FX=>this%Molecule%SiteQuadrupole(j)%FX
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%FY=>this%Molecule%SiteQuadrupole(j)%FY
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%FZ=>this%Molecule%SiteQuadrupole(j)%FZ
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%TX=>this%Molecule%SiteQuadrupole(j)%TX
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%TY=>this%Molecule%SiteQuadrupole(j)%TY
                    this%Molecule%Unit(i)%SiteQuadrupole(j)%TZ=>this%Molecule%SiteQuadrupole(j)%TZ
                    this%Molecule%SiteQuadrupole(j)%PX=> this%Molecule%Unit(i)%PX
                    this%Molecule%SiteQuadrupole(j)%PY=> this%Molecule%Unit(i)%PY
                    this%Molecule%SiteQuadrupole(j)%PZ=> this%Molecule%Unit(i)%PZ
            end if
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
       if ( this%Molecule%NLJ126>0 ) then
          do j = 1, this%Molecule%NLJ126
             if (this%Molecule%SiteLJ126(j)%SiteId==SiteId1) then
                 this%Molecule%IdfBond(i)%RX1=>this%Molecule%SiteLJ126(j)%RX(:)
                 this%Molecule%IdfBond(i)%RY1=>this%Molecule%SiteLJ126(j)%RY(:)
                 this%Molecule%IdfBond(i)%RZ1=>this%Molecule%SiteLJ126(j)%RZ(:)
                 this%Molecule%IdfBond(i)%FX1=>this%Molecule%SiteLJ126(j)%FX(:)
                 this%Molecule%IdfBond(i)%FY1=>this%Molecule%SiteLJ126(j)%FY(:)
                 this%Molecule%IdfBond(i)%FZ1=>this%Molecule%SiteLJ126(j)%FZ(:)
                 this%Molecule%IdfBond(i)%PX1=>this%Molecule%SiteLJ126(j)%PX(:)   ! For calculation of virial
                 this%Molecule%IdfBond(i)%PY1=>this%Molecule%SiteLJ126(j)%PY(:)
                 this%Molecule%IdfBond(i)%PZ1=>this%Molecule%SiteLJ126(j)%PZ(:)
                 Site1 = .true.
             else if (this%Molecule%SiteLJ126(j)%SiteId==SiteId2) then
                 this%Molecule%IdfBond(i)%RX2=>this%Molecule%SiteLJ126(j)%RX(:)
                 this%Molecule%IdfBond(i)%RY2=>this%Molecule%SiteLJ126(j)%RY(:)
                 this%Molecule%IdfBond(i)%RZ2=>this%Molecule%SiteLJ126(j)%RZ(:)
                 this%Molecule%IdfBond(i)%FX2=>this%Molecule%SiteLJ126(j)%FX(:)
                 this%Molecule%IdfBond(i)%FY2=>this%Molecule%SiteLJ126(j)%FY(:)
                 this%Molecule%IdfBond(i)%FZ2=>this%Molecule%SiteLJ126(j)%FZ(:)
                 this%Molecule%IdfBond(i)%PX2=>this%Molecule%SiteLJ126(j)%PX(:)
                 this%Molecule%IdfBond(i)%PY2=>this%Molecule%SiteLJ126(j)%PY(:)
                 this%Molecule%IdfBond(i)%PZ2=>this%Molecule%SiteLJ126(j)%PZ(:)
                 Site2 = .true.
              end if
              if (Site1 .and. Site2) then
                  exit
              end if
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
             if (Site1 .and. Site2) then
                  exit
             end if
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
       if ( this%Molecule%NLJ126>0 ) then
          do j = 1, this%Molecule%NLJ126
             if (this%Molecule%SiteLJ126(j)%SiteId==SiteId1) then
                 this%Molecule%IdfAngle(i)%RX1=>this%Molecule%SiteLJ126(j)%RX(:)
                 this%Molecule%IdfAngle(i)%RY1=>this%Molecule%SiteLJ126(j)%RY(:)
                 this%Molecule%IdfAngle(i)%RZ1=>this%Molecule%SiteLJ126(j)%RZ(:)
                 this%Molecule%IdfAngle(i)%FX1=>this%Molecule%SiteLJ126(j)%FX(:)
                 this%Molecule%IdfAngle(i)%FY1=>this%Molecule%SiteLJ126(j)%FY(:)
                 this%Molecule%IdfAngle(i)%FZ1=>this%Molecule%SiteLJ126(j)%FZ(:)
                 Site1 = .true.
             else if (this%Molecule%SiteLJ126(j)%SiteId==SiteId2) then
                 this%Molecule%IdfAngle(i)%RX2=>this%Molecule%SiteLJ126(j)%RX(:)
                 this%Molecule%IdfAngle(i)%RY2=>this%Molecule%SiteLJ126(j)%RY(:)
                 this%Molecule%IdfAngle(i)%RZ2=>this%Molecule%SiteLJ126(j)%RZ(:)
                 this%Molecule%IdfAngle(i)%FX2=>this%Molecule%SiteLJ126(j)%FX(:)
                 this%Molecule%IdfAngle(i)%FY2=>this%Molecule%SiteLJ126(j)%FY(:)
                 this%Molecule%IdfAngle(i)%FZ2=>this%Molecule%SiteLJ126(j)%FZ(:)
                 Site2 = .true.
             else if (this%Molecule%SiteLJ126(j)%SiteId==SiteId3) then
                 this%Molecule%IdfAngle(i)%RX3=>this%Molecule%SiteLJ126(j)%RX(:)
                 this%Molecule%IdfAngle(i)%RY3=>this%Molecule%SiteLJ126(j)%RY(:)
                 this%Molecule%IdfAngle(i)%RZ3=>this%Molecule%SiteLJ126(j)%RZ(:)
                 this%Molecule%IdfAngle(i)%FX3=>this%Molecule%SiteLJ126(j)%FX(:)
                 this%Molecule%IdfAngle(i)%FY3=>this%Molecule%SiteLJ126(j)%FY(:)
                 this%Molecule%IdfAngle(i)%FZ3=>this%Molecule%SiteLJ126(j)%FZ(:)
                 Site3 = .true.
              end if
              if (Site1 .and. Site2 .and. Site3) then
                  exit
              end if
            end do
       end if
       if((.not.Site1 .or. .not. Site2 .or. .not. Site3) .and. (this%Molecule%NCharge > 0) ) then
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
             if (Site1 .and. Site2 .and. Site3) then
                  exit
             end if
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
       if ( this%Molecule%NLJ126>0 ) then
          do j = 1, this%Molecule%NLJ126
             if (this%Molecule%SiteLJ126(j)%SiteId==SiteId1) then
                 this%Molecule%IdfDihedral(i)%RX1=>this%Molecule%SiteLJ126(j)%RX(:)
                 this%Molecule%IdfDihedral(i)%RY1=>this%Molecule%SiteLJ126(j)%RY(:)
                 this%Molecule%IdfDihedral(i)%RZ1=>this%Molecule%SiteLJ126(j)%RZ(:)
                 this%Molecule%IdfDihedral(i)%FX1=>this%Molecule%SiteLJ126(j)%FX(:)
                 this%Molecule%IdfDihedral(i)%FY1=>this%Molecule%SiteLJ126(j)%FY(:)
                 this%Molecule%IdfDihedral(i)%FZ1=>this%Molecule%SiteLJ126(j)%FZ(:)
!                 this%Molecule%IdfDihedral(i)%PX1=>this%Molecule%SiteLJ126(j)%PX(:) !for 1,4 intramolecular interaction virial
!                 this%Molecule%IdfDihedral(i)%PY1=>this%Molecule%SiteLJ126(j)%PY(:)
!                 this%Molecule%IdfDihedral(i)%PZ1=>this%Molecule%SiteLJ126(j)%PZ(:)
!                 this%Molecule%IdfDihedral(i)%Sigma1=this%Molecule%SiteLJ126(j)%sig
!                 this%Molecule%IdfDihedral(i)%Epsilon1=this%Molecule%SiteLJ126(j)%eps
                 Site1 = .true.
             else if (this%Molecule%SiteLJ126(j)%SiteId==SiteId2) then
                 this%Molecule%IdfDihedral(i)%RX2=>this%Molecule%SiteLJ126(j)%RX(:)
                 this%Molecule%IdfDihedral(i)%RY2=>this%Molecule%SiteLJ126(j)%RY(:)
                 this%Molecule%IdfDihedral(i)%RZ2=>this%Molecule%SiteLJ126(j)%RZ(:)
                 this%Molecule%IdfDihedral(i)%FX2=>this%Molecule%SiteLJ126(j)%FX(:)
                 this%Molecule%IdfDihedral(i)%FY2=>this%Molecule%SiteLJ126(j)%FY(:)
                 this%Molecule%IdfDihedral(i)%FZ2=>this%Molecule%SiteLJ126(j)%FZ(:)
                 Site2 = .true.
             else if (this%Molecule%SiteLJ126(j)%SiteId==SiteId3) then
                 this%Molecule%IdfDihedral(i)%RX3=>this%Molecule%SiteLJ126(j)%RX(:)
                 this%Molecule%IdfDihedral(i)%RY3=>this%Molecule%SiteLJ126(j)%RY(:)
                 this%Molecule%IdfDihedral(i)%RZ3=>this%Molecule%SiteLJ126(j)%RZ(:)
                 this%Molecule%IdfDihedral(i)%FX3=>this%Molecule%SiteLJ126(j)%FX(:)
                 this%Molecule%IdfDihedral(i)%FY3=>this%Molecule%SiteLJ126(j)%FY(:)
                 this%Molecule%IdfDihedral(i)%FZ3=>this%Molecule%SiteLJ126(j)%FZ(:)
                 Site3 = .true.
             else if (this%Molecule%SiteLJ126(j)%SiteId==SiteId4) then
                 this%Molecule%IdfDihedral(i)%RX4=>this%Molecule%SiteLJ126(j)%RX(:)
                 this%Molecule%IdfDihedral(i)%RY4=>this%Molecule%SiteLJ126(j)%RY(:)
                 this%Molecule%IdfDihedral(i)%RZ4=>this%Molecule%SiteLJ126(j)%RZ(:)
                 this%Molecule%IdfDihedral(i)%FX4=>this%Molecule%SiteLJ126(j)%FX(:)
                 this%Molecule%IdfDihedral(i)%FY4=>this%Molecule%SiteLJ126(j)%FY(:)
                 this%Molecule%IdfDihedral(i)%FZ4=>this%Molecule%SiteLJ126(j)%FZ(:)
!                 this%Molecule%IdfDihedral(i)%PX4=>this%Molecule%SiteLJ126(j)%PX(:)!for 1,4 intramolecular interaction virial
!                 this%Molecule%IdfDihedral(i)%PY4=>this%Molecule%SiteLJ126(j)%PY(:)
!                 this%Molecule%IdfDihedral(i)%PZ4=>this%Molecule%SiteLJ126(j)%PZ(:)
!                 this%Molecule%IdfDihedral(i)%Sigma4=this%Molecule%SiteLJ126(j)%sig
!                 this%Molecule%IdfDihedral(i)%Epsilon4=this%Molecule%SiteLJ126(j)%eps
                 Site4 = .true.
              end if
              if (Site1 .and. Site2 .and. Site3 .and. Site4) then
                  exit
              end if
            end do
       end if
       if((.not.Site1 .or. .not. Site2 .or. .not. Site3 .or. .not. Site4) .and. (this%Molecule%NCharge > 0) ) then
          do j = 1, this%Molecule%NCharge
             if (this%Molecule%SiteCharge(j)%SiteId==SiteId1) then
                 this%Molecule%IdfDihedral(i)%RX1=>this%Molecule%SiteCharge(j)%RX(:)
                 this%Molecule%IdfDihedral(i)%RY1=>this%Molecule%SiteCharge(j)%RY(:)
                 this%Molecule%IdfDihedral(i)%RZ1=>this%Molecule%SiteCharge(j)%RZ(:)
                 this%Molecule%IdfDihedral(i)%FX1=>this%Molecule%SiteCharge(j)%FX(:)
                 this%Molecule%IdfDihedral(i)%FY1=>this%Molecule%SiteCharge(j)%FY(:)
                 this%Molecule%IdfDihedral(i)%FZ1=>this%Molecule%SiteCharge(j)%FZ(:)
!                 this%Molecule%IdfDihedral(i)%PX1=>this%Molecule%SiteCharge(j)%PX(:)!for 1,4 intramolecular interaction virial
!                 this%Molecule%IdfDihedral(i)%PY1=>this%Molecule%SiteCharge(j)%PY(:)
!                 this%Molecule%IdfDihedral(i)%PZ1=>this%Molecule%SiteCharge(j)%PZ(:)
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
!                 this%Molecule%IdfDihedral(i)%PX4=>this%Molecule%SiteCharge(j)%PX(:)
!                 this%Molecule%IdfDihedral(i)%PY4=>this%Molecule%SiteCharge(j)%PY(:)
!                 this%Molecule%IdfDihedral(i)%PZ4=>this%Molecule%SiteCharge(j)%PZ(:)
                 Site4 = .true.
             end if
             if (Site1 .and. Site2 .and. Site3 .and. Site4) then
                  exit
             end if
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
!       allocate( this%NStateBF( nf ), STAT = stat )
!       call AllocationError( stat, 'fluctuating particle states', nf )
!       allocate( this%BFSumState( nf ), STAT = stat )
!       call AllocationError( stat, 'fluctuating particle states', nf )
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

    ! Update log file
    write( IOBuffer, '("Memory for ", A, " allocated successfully")' ) &
&     trim( this%PotModFileName )
    call LogWrite

   contains

    subroutine binar_search (array, Id, treffer, index)

      ! Declare arguments
      integer, dimension(:), intent( in ) :: array
      integer, intent( in )               :: Id
      logical, intent( out )              :: treffer
      integer, intent( out )              :: index

      ! Declare local variables
      integer                             :: anfang, ende, mitte

      anfang = 1
      ende = size (array)
      do
         if ( anfang == ende ) exit
         mitte = (anfang + ende)*0.5
         if ( id <= array(mitte) ) then
           ende = mitte
         else
           anfang = mitte + 1
         end if
      end do
      index = anfang
      treffer = (id == array(index))

    end subroutine binar_search


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
    if( associated( this%Pm0 ) ) deallocate( this%Pm0 )
    if( associated( this%Pm0old ) ) deallocate( this%Pm0old )
    if( associated( this%P0Save ) ) deallocate( this%P0Save )
    if( associated( this%P0 ) ) deallocate( this%P0 )
    if( associated( this%P1 ) ) deallocate( this%P1 )
    if( associated( this%P2 ) ) deallocate( this%P2 )
    if( associated( this%P3 ) ) deallocate( this%P3 )
    if( associated( this%P4 ) ) deallocate( this%P4 )
    if( associated( this%P5 ) ) deallocate( this%P5 )

    ! Displacement
    if( associated( this%Disp ) ) deallocate( this%Disp )

    ! Total forces
    if( associated( this%F ) ) deallocate( this%F )

    ! Quaternion parameters and their derivatives
    if( associated( this%Qm0 ) ) deallocate( this%Qm0 )
    if( associated( this%Q0Save ) ) deallocate( this%Q0Save )
    if( associated( this%Q0 ) ) deallocate( this%Q0 )
    if( associated( this%Q0tmp ) ) deallocate( this%Q0tmp )
    if( associated( this%Q1 ) ) deallocate( this%Q1 )
    if( associated( this%Q2 ) ) deallocate( this%Q2 )
    if( associated( this%Q3 ) ) deallocate( this%Q3 )
    if( associated( this%Q4 ) ) deallocate( this%Q4 )

    ! Angular velocities and their derivatives
    if( associated( this%W0 ) ) deallocate( this%W0 )
    if( associated( this%W1 ) ) deallocate( this%W1 )
    if( associated( this%W2 ) ) deallocate( this%W2 )
    if( associated( this%W3 ) ) deallocate( this%W3 )
    if( associated( this%W4 ) ) deallocate( this%W4 )

    ! Total torques
    if( associated( this%T ) ) deallocate( this%T )

    ! Total dipole moment of molecules for reaction field
    if( associated( this%MueX ) ) deallocate( this%MueX )
    if( associated( this%MueY ) ) deallocate( this%MueY )
    if( associated( this%MueZ ) ) deallocate( this%MueZ )

    ! Torques from reaction field
    if( associated( this%tRFX ) ) deallocate( this%tRFX )
    if( associated( this%tRFY ) ) deallocate( this%tRFY )
    if( associated( this%tRFZ ) ) deallocate( this%tRFZ )

    ! Total dipole moment of test particles for reaction field
    if( associated( this%MueXTest ) ) deallocate( this%MueXTest )
    if( associated( this%MueYTest ) ) deallocate( this%MueYTest )
    if( associated( this%MueZTest ) ) deallocate( this%MueZTest )

    ! Gear corrector local arrays
    if( associated( this%Corr0 ) ) deallocate( this%Corr0 )
    if( associated( this%Corr1 ) ) deallocate( this%Corr1 )

    ! Site positions, orientations, forces and torques
    do i = 1, this%Molecule%NLJ126
      call Deallocate( this%Molecule%SiteLJ126(i) )
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
!     if( associated( this%NStateBF ) ) then
!       deallocate( this%NStateBF )
!     end if
!     if( associated( this%BFSumState ) ) then
!       deallocate( this%BFSumState )
!     end if
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
    end if
!DEBUG

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

! Reaction Field Check
   if ((LongRange .eq. RField) .and. (CutOffMode .eq. SiteSite) ) then
     if (this%Molecule%NCharge > 0 ) then
       write (ErrorBuffer,'("Reaction Field in combination with SiteSite-Cutoff doesnot support partial charges")')
     end if
   end if


   if ((LongRange .eq. RField) .and. (abs(q) .ge. 1e-1)) then
     write (ErrorBuffer,'("You have a non-neutral component.\n NetCharge norm&red = ", F15.10, "\n Conflicts with ReactionField")') q
     call Error
   end if

   if ( ((EnsembleType .eq. EnsembleTypeGE) .or. (EnsembleType .eq. EnsembleTypeHA)) &
&         .and. (abs(q) .ge. 1e-4) ) then
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
    integer :: i, j, k
    integer :: nu

    nu = this%Molecule%NUnit
    ! Set random linear velocities
    do k=1,nu
      do i = 1, 3
        do j = 1, this%NPart
          this%P1(j, i,k) = rnd( -1._RK, 1._RK )
        end do
      end do
    end do

    ! Nullify angular velocities
    if( this%Molecule%isElongated ) this%W0(:, :,:) = 0._RK

  end subroutine TComponent_InitVelocities



!==============================================================!
!  Subroutine TComponent_InitIntegratorGear                    !
!==============================================================!

  subroutine TComponent_InitIntegratorGear( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Local variables
    integer          :: i
    integer          :: nu

    nu = this%Molecule%NUnit

    ! Zero accelerations
    this%P2(:, :, :) = 0._RK
    this%P3(:, :, :) = 0._RK
    this%P4(:, :, :) = 0._RK
    this%P5(:, :, :) = 0._RK
    do i=1,nu
!      if( this%Molecule%Unit(i)%isElongated ) then
        this%Q1(:, :, i) = 0._RK
        this%Q2(:, :, i) = 0._RK
        this%Q3(:, :, i) = 0._RK
        this%Q4(:, :, i) = 0._RK
        this%W1(:, :, i) = 0._RK
        this%W2(:, :, i) = 0._RK
        this%W3(:, :, i) = 0._RK
        this%W4(:, :, i) = 0._RK
!    end if
   end do

  end subroutine TComponent_InitIntegratorGear



!==============================================================!
!  Subroutine TComponent_InitIntegratorLeap                    !
!==============================================================!

  subroutine TComponent_InitIntegratorLeap( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Local variables
    integer          :: i, nu
    nu = this%Molecule%NUnit

    ! Zero accelerations
    this%P2(:, :, :) = 0._RK
    do i = 1, nu
!     if( this%Molecule%Unit(i)%isElongated ) then
        this%Q1(:, :, i) = 0._RK
        this%W1(:, :, i) = 0._RK
!      end if
    end do

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
!  Subroutine TComponent_RemoveNetMomentum                      !
!==============================================================!

  subroutine TComponent_RemoveNetMomentum( this,nu )

    implicit none

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: nu

    ! Declare local variables
    real(RK) :: P(3,nu), L(3,nu)
    integer :: i, j,k
    real(RK) :: Pim(nu)

    ! Return if zero particles in component
    if( this%NPart == 0 ) return

!    ! Calculate net momentum
!    P(:) = 0._RK
!    L(:) = 0._RK
!    do i = 1, 3
!      P(i) = P(i) &
!&       + this%Molecule%Mass * sum( this%P1(1:this%NPart, i) )
!      if( i <= this%Molecule%NDFRot ) &
!&       L(i) = L(i) &
!&         + this%Molecule%MOI(i) * sum( this%W0(1:this%NPart, i) )
!    end do
!    P(:) = P(:) / this%NPart
!    L(:) = L(:) / this%NPart

    ! Remove net momentum
    do k = 1, nu
      P(:, k) = 0._RK
      L(:, k) = 0._RK
      ! Calculate net momentum
      do i = 1, 3
        P(i, k) = P(i, k) &
&         + this%Molecule%Unit(k)%Mass * sum( this%P1(1:this%NPart, i, k) )
        if( i <= this%Molecule%Unit(k)%NDFRot ) &
&         L(i, k) = L(i, k) &
&           + this%Molecule%Unit(k)%MOI(i) * sum( this%W0(1:this%NPart, i, k) )
      end do
      P(:, k) = P(:, k) / this%NPart
      L(:, k) = L(:, k) / this%NPart

      ! Remove net momentum
      do i = 1, 3
        Pim(k) = P(i, k) / this%Molecule%Unit(k)%Mass




        do j = 1, this%NPart
          this%P1(j, i, k) = this%P1(j, i, k) - Pim(k)
        end do
        if( i <= this%Molecule%Unit(k)%NDFRot ) then
          Pim(k) = L(i, k) / this%Molecule%Unit(k)%MOI(i)
          do j = 1, this%NPart
            this%W0(j, i, k) = this%W0(j, i, k) - Pim(k)
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
    integer :: i, k, nu

    nu = this%Molecule%NUnit

    ! Calculate translational kinetic energy
    this%EKinTran = 0._RK
    do  k = 1,  nu
      this%EKinTran = this%EkinTran+this%Molecule%Unit(k)%Mass * TimeStepSquaredInv2 &
&       * sum( this%P1(1:this%NPart, :, k)**2 ) * this%BoxLength**2
    end do

    ! Calculate rotational kinetic energy
    this%EKinRot = 0._RK
    do k = 1, nu
      do i = 1, this%Molecule%Unit(k)%NDFRot
        this%EKinRot = this%EKinRot + this%Molecule%Unit(k)%MOI(i) * .5_RK &
&         * sum( this%W0(1:this%NPart, i, k)**2 )
      end do
    end do

  end subroutine TComponent_CalculateEKin



! ! ! !==============================================================!
! ! ! !  Subroutine TComponent_Mol2Atom                              !
! ! ! !==============================================================!
! ! !
! ! !   subroutine TComponent_Mol2Atom( this, np )
! ! !
! ! !     implicit none
! ! !
! ! !     ! Include MPI header
! ! ! #if MPI_VER > 0
! ! !     include 'mpif.h'
! ! ! #endif
! ! !
! ! !     ! Declare arguments
! ! !     type(TComponent)    :: this
! ! !     integer, intent(in) :: np
! ! !
! ! !     ! Declare local variables
! ! !     real(RK)                       :: BoxLengthInv
! ! !     real(RK)                       :: PX(np), PY(np), PZ(np)
! ! !     real(RK)                       :: q1, q2, q3, q4, qinv
! ! !     real(RK)                       :: A11(np), A12(np), A13(np)
! ! !     real(RK)                       :: A21(np), A22(np), A23(np)
! ! !     real(RK)                       :: A31(np), A32(np), A33(np)
! ! !     real(RK)                       :: r1, r2, r3, or1, or2, or3
! ! !     real(RK)                       :: mue1, mue2, mue3
! ! !     type(TSiteLJ126), pointer      :: pLJ126
! ! !     type(TSiteCharge), pointer     :: pCharge
! ! !     type(TSiteDipole), pointer     :: pDipole
! ! !     type(TSiteQuadrupole), pointer :: pQuadrupole
! ! !     integer                        :: i, j
! ! !
! ! !     ! Broadcast positions and orientations to all processes
! ! ! #if MPI_VER > 0
! ! !     call MPI_Bcast( this%P0(:, :), size( this%P0 ), &
! ! ! &     MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
! ! !     if( this%Molecule%isElongated ) &
! ! ! &     call MPI_Bcast( this%Q0(:, :), size( this%Q0 ), &
! ! ! &       MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
! ! ! #endif
! ! !
! ! !     ! Assign local variables
! ! !     BoxLengthInv = 1._RK / this%BoxLength
! ! !
! ! !     ! Check number of rotation axes
! ! !     if( this%Molecule%isElongated ) then
! ! !
! ! !       ! Loop over molecules
! ! !       do i = 1, np
! ! !         ! Positions and quaternions of particle i
! ! !         PX(i) = this%P0(i, 1)
! ! !         PY(i) = this%P0(i, 2)
! ! !         PZ(i) = this%P0(i, 3)
! ! !         q1 = this%Q0(i, 1)
! ! !         q2 = this%Q0(i, 2)
! ! !         q3 = this%Q0(i, 3)
! ! !         q4 = this%Q0(i, 4)
! ! !
! ! !         ! Normalise quaternions
! ! ! #if ARCH == 3
! ! !         qinv = rsqrt( q1**2 + q2**2 + q3**2 + q4**2 )
! ! ! #else
! ! !         qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
! ! ! #endif
! ! !         q1 = q1 * qinv
! ! !         q2 = q2 * qinv
! ! !         q3 = q3 * qinv
! ! !         q4 = q4 * qinv
! ! !         this%Q0(i, 1) = q1
! ! !         this%Q0(i, 2) = q2
! ! !         this%Q0(i, 3) = q3
! ! !         this%Q0(i, 4) = q4
! ! !
! ! !         ! Calculate rotation matrix elements
! ! !         A11(i) = q1**2 + q2**2 - q3**2 - q4**2
! ! !         A12(i) = 2._RK * (q2 * q3 + q1 * q4)
! ! !         A13(i) = 2._RK * (q2 * q4 - q1 * q3)
! ! !         A21(i) = 2._RK * (q2 * q3 - q1 * q4)
! ! !         A22(i) = q1**2 - q2**2 + q3**2 - q4**2
! ! !         A23(i) = 2._RK * (q3 * q4 + q1 * q2)
! ! !         A31(i) = 2._RK * (q2 * q4 + q1 * q3)
! ! !         A32(i) = 2._RK * (q3 * q4 - q1 * q2)
! ! !         A33(i) = q1**2 - q2**2 - q3**2 + q4**2
! ! !       end do
! ! !
! ! !       ! Loop over LJ126 sites in molecule
! ! !       do j = 1, this%Molecule%NLJ126
! ! !         pLJ126 => this%Molecule%SiteLJ126(j)
! ! !         r1 = pLJ126%r(1) * BoxLengthInv
! ! !         r2 = pLJ126%r(2) * BoxLengthInv
! ! !         r3 = pLJ126%r(3) * BoxLengthInv
! ! !         do i = 1, np
! ! !           pLJ126%RX(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
! ! !           pLJ126%RY(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
! ! !           pLJ126%RZ(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
! ! !         end do
! ! !       end do
! ! !
! ! !       ! Loop over charge sites in molecule
! ! !       do j = 1, this%Molecule%NCharge
! ! !         pCharge => this%Molecule%SiteCharge(j)
! ! !         r1 = pCharge%r(1) * BoxLengthInv
! ! !         r2 = pCharge%r(2) * BoxLengthInv
! ! !         r3 = pCharge%r(3) * BoxLengthInv
! ! !         do i = 1, np
! ! !           pCharge%RX(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
! ! !           pCharge%RY(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
! ! !           pCharge%RZ(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
! ! !         end do
! ! !       end do
! ! !
! ! !       ! Loop over dipole sites in molecule
! ! !       do j = 1, this%Molecule%NDipole
! ! !         pDipole => this%Molecule%SiteDipole(j)
! ! !         r1 = pDipole%r(1) * BoxLengthInv
! ! !         r2 = pDipole%r(2) * BoxLengthInv
! ! !         r3 = pDipole%r(3) * BoxLengthInv
! ! !         or1 = pDipole%or(1)
! ! !         or2 = pDipole%or(2)
! ! !         or3 = pDipole%or(3)
! ! !         do i = 1, np
! ! !           pDipole%RX(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
! ! !           pDipole%RY(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
! ! !           pDipole%RZ(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
! ! !           pDipole%OX(i) = or1 * A11(i) + or2 * A21(i) + or3 * A31(i)
! ! !           pDipole%OY(i) = or1 * A12(i) + or2 * A22(i) + or3 * A32(i)
! ! !           pDipole%OZ(i) = or1 * A13(i) + or2 * A23(i) + or3 * A33(i)
! ! !         end do
! ! !       end do
! ! !
! ! !       ! Loop over quadrupole sites in molecule
! ! !       do j = 1, this%Molecule%NQuadrupole
! ! !         pQuadrupole => this%Molecule%SiteQuadrupole(j)
! ! !         r1 = pQuadrupole%r(1) * BoxLengthInv
! ! !         r2 = pQuadrupole%r(2) * BoxLengthInv
! ! !         r3 = pQuadrupole%r(3) * BoxLengthInv
! ! !         or1 = pQuadrupole%or(1)
! ! !         or2 = pQuadrupole%or(2)
! ! !         or3 = pQuadrupole%or(3)
! ! !         do i = 1, np
! ! !           pQuadrupole%RX(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
! ! !           pQuadrupole%RY(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
! ! !           pQuadrupole%RZ(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
! ! !           pQuadrupole%OX(i) = or1 * A11(i) + or2 * A21(i) + or3 * A31(i)
! ! !           pQuadrupole%OY(i) = or1 * A12(i) + or2 * A22(i) + or3 * A32(i)
! ! !           pQuadrupole%OZ(i) = or1 * A13(i) + or2 * A23(i) + or3 * A33(i)
! ! !         end do
! ! !       end do
! ! !
! ! !       ! Rotate total dipole moment
! ! !       if( CutoffMode .eq. CenterofMass ) then
! ! !         mue1 = this%Molecule%Mue(1)
! ! !         mue2 = this%Molecule%Mue(2)
! ! !         mue3 = this%Molecule%Mue(3)
! ! !         do i = 1, np
! ! !           this%MueX(i) = mue1 * A11(i) + mue2 * A21(i) + mue3 * A31(i)
! ! !           this%MueY(i) = mue1 * A12(i) + mue2 * A22(i) + mue3 * A32(i)
! ! !           this%MueZ(i) = mue1 * A13(i) + mue2 * A23(i) + mue3 * A33(i)
! ! !         end do
! ! !       end if
! ! !
! ! !     else
! ! !
! ! !       ! Loop over LJ126 sites in molecule
! ! !       do i = 1, this%Molecule%NLJ126
! ! !         pLJ126 => this%Molecule%SiteLJ126(i)
! ! !         do j = 1, np
! ! !           pLJ126%RX(j) = this%P0(j, 1)
! ! !           pLJ126%RY(j) = this%P0(j, 2)
! ! !           pLJ126%RZ(j) = this%P0(j, 3)
! ! !         end do
! ! !       end do
! ! !
! ! !       ! Loop over charge sites in molecule
! ! !       if (LongRange .ne. RField) then
! ! !         do i = 1, this%Molecule%NCharge
! ! !           pCharge => this%Molecule%SiteCharge(i)
! ! !           do j = 1, np
! ! !             pCharge%RX(j) = this%P0(j,1)
! ! !             pCharge%RY(j) = this%P0(j,2)
! ! !             pCharge%RZ(j) = this%P0(j,3)
! ! !           end do
! ! !         end do
! ! !       end if
! ! !
! ! !     end if
! ! !
! ! !   end subroutine TComponent_Mol2Atom
! ! !
! ! !
! ! !
! ! ! !==============================================================!
! ! ! !  Subroutine TComponent_Mol2Atom1                             !
! ! ! !==============================================================!
! ! !
! ! !   subroutine TComponent_Mol2Atom1( this, n )
! ! !
! ! !     implicit none
! ! !
! ! !     ! Declare arguments
! ! !     type(TComponent)    :: this
! ! !     integer, intent(in) :: n
! ! !
! ! !     ! Declare local variables
! ! !     real(RK)                       :: BoxLengthInv
! ! !     real(RK)                       :: PXi, PYi, PZi
! ! !     real(RK)                       :: q1, q2, q3, q4, qinv
! ! !     real(RK)                       :: A11, A12, A13, A21, A22, A23, &
! ! ! &                                     A31, A32, A33
! ! !     real(RK)                       :: r1, r2, r3, or1, or2, or3
! ! !     real(RK)                       :: mue1, mue2, mue3
! ! !     type(TSiteLJ126), pointer      :: pLJ126
! ! !     type(TSiteCharge), pointer     :: pCharge
! ! !     type(TSiteDipole), pointer     :: pDipole
! ! !     type(TSiteQuadrupole), pointer :: pQuadrupole
! ! !     integer                        :: i
! ! !
! ! !     ! Assign local variables
! ! !     BoxLengthInv = 1._RK / this%BoxLength
! ! !
! ! !     ! Positions of particle n
! ! !     PXi = this%P0(n, 1)
! ! !     PYi = this%P0(n, 2)
! ! !     PZi = this%P0(n, 3)
! ! !
! ! !     ! Check number of rotation axes
! ! !     if( this%Molecule%isElongated ) then
! ! !
! ! !       ! Normalise quaternions
! ! !       q1 = this%Q0(n, 1)
! ! !       q2 = this%Q0(n, 2)
! ! !       q3 = this%Q0(n, 3)
! ! !       q4 = this%Q0(n, 4)
! ! ! #if ARCH == 3
! ! !       qinv = rsqrt( q1**2 + q2**2 + q3**2 + q4**2 )
! ! ! #else
! ! !       qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
! ! ! #endif
! ! !       q1 = q1 * qinv
! ! !       q2 = q2 * qinv
! ! !       q3 = q3 * qinv
! ! !       q4 = q4 * qinv
! ! !       this%Q0(n, 1) = q1
! ! !       this%Q0(n, 2) = q2
! ! !       this%Q0(n, 3) = q3
! ! !       this%Q0(n, 4) = q4
! ! !
! ! !       ! Calculate rotation matrix elements
! ! !       A11 = q1**2 + q2**2 - q3**2 - q4**2
! ! !       A12 = 2._RK * (q2 * q3 + q1 * q4)
! ! !       A13 = 2._RK * (q2 * q4 - q1 * q3)
! ! !       A21 = 2._RK * (q2 * q3 - q1 * q4)
! ! !       A22 = q1**2 - q2**2 + q3**2 - q4**2
! ! !       A23 = 2._RK * (q3 * q4 + q1 * q2)
! ! !       A31 = 2._RK * (q2 * q4 + q1 * q3)
! ! !       A32 = 2._RK * (q3 * q4 - q1 * q2)
! ! !       A33 = q1**2 - q2**2 - q3**2 + q4**2
! ! !
! ! !       ! Loop over LJ126 sites in molecule
! ! !       do i = 1, this%Molecule%NLJ126
! ! !         pLJ126 => this%Molecule%SiteLJ126(i)
! ! !         r1 = pLJ126%r(1) * BoxLengthInv
! ! !         r2 = pLJ126%r(2) * BoxLengthInv
! ! !         r3 = pLJ126%r(3) * BoxLengthInv
! ! !         pLJ126%RX(n) = PXi + r1 * A11 + r2 * A21 + r3 * A31
! ! !         pLJ126%RY(n) = PYi + r1 * A12 + r2 * A22 + r3 * A32
! ! !         pLJ126%RZ(n) = PZi + r1 * A13 + r2 * A23 + r3 * A33
! ! !       end do
! ! !
! ! !       ! Loop over charge sites in molecule
! ! !       do i = 1, this%Molecule%NCharge
! ! !         pCharge => this%Molecule%SiteCharge(i)
! ! !         r1 = pCharge%r(1) * BoxLengthInv
! ! !         r2 = pCharge%r(2) * BoxLengthInv
! ! !         r3 = pCharge%r(3) * BoxLengthInv
! ! !         pCharge%RX(n) = PXi + r1 * A11 + r2 * A21 + r3 * A31
! ! !         pCharge%RY(n) = PYi + r1 * A12 + r2 * A22 + r3 * A32
! ! !         pCharge%RZ(n) = PZi + r1 * A13 + r2 * A23 + r3 * A33
! ! !       end do
! ! !
! ! !       ! Loop over dipole sites in molecule
! ! !       do i = 1, this%Molecule%NDipole
! ! !         pDipole => this%Molecule%SiteDipole(i)
! ! !         r1 = pDipole%r(1) * BoxLengthInv
! ! !         r2 = pDipole%r(2) * BoxLengthInv
! ! !         r3 = pDipole%r(3) * BoxLengthInv
! ! !         or1 = pDipole%or(1)
! ! !         or2 = pDipole%or(2)
! ! !         or3 = pDipole%or(3)
! ! !         pDipole%RX(n) = PXi + r1 * A11 + r2 * A21 + r3 * A31
! ! !         pDipole%RY(n) = PYi + r1 * A12 + r2 * A22 + r3 * A32
! ! !         pDipole%RZ(n) = PZi + r1 * A13 + r2 * A23 + r3 * A33
! ! !         pDipole%OX(n) = or1 * A11 + or2 * A21 + or3 * A31
! ! !         pDipole%OY(n) = or1 * A12 + or2 * A22 + or3 * A32
! ! !         pDipole%OZ(n) = or1 * A13 + or2 * A23 + or3 * A33
! ! !       end do
! ! !
! ! !       ! Loop over quadrupole sites in molecule
! ! !       do i = 1, this%Molecule%NQuadrupole
! ! !         pQuadrupole => this%Molecule%SiteQuadrupole(i)
! ! !         r1 = pQuadrupole%r(1) * BoxLengthInv
! ! !         r2 = pQuadrupole%r(2) * BoxLengthInv
! ! !         r3 = pQuadrupole%r(3) * BoxLengthInv
! ! !         or1 = pQuadrupole%or(1)
! ! !         or2 = pQuadrupole%or(2)
! ! !         or3 = pQuadrupole%or(3)
! ! !         pQuadrupole%RX(n) = PXi + r1 * A11 + r2 * A21 + r3 * A31
! ! !         pQuadrupole%RY(n) = PYi + r1 * A12 + r2 * A22 + r3 * A32
! ! !         pQuadrupole%RZ(n) = PZi + r1 * A13 + r2 * A23 + r3 * A33
! ! !         pQuadrupole%OX(n) = or1 * A11 + or2 * A21 + or3 * A31
! ! !         pQuadrupole%OY(n) = or1 * A12 + or2 * A22 + or3 * A32
! ! !         pQuadrupole%OZ(n) = or1 * A13 + or2 * A23 + or3 * A33
! ! !       end do
! ! !
! ! !       ! Rotate total dipole moment
! ! !       if( CutoffMode .eq. CenterofMass ) then
! ! !         mue1 = this%Molecule%Mue(1)
! ! !         mue2 = this%Molecule%Mue(2)
! ! !         mue3 = this%Molecule%Mue(3)
! ! !         this%MueX(n) = mue1 * A11 + mue2 * A21 + mue3 * A31
! ! !         this%MueY(n) = mue1 * A12 + mue2 * A22 + mue3 * A32
! ! !         this%MueZ(n) = mue1 * A13 + mue2 * A23 + mue3 * A33
! ! !       end if
! ! !
! ! !     else
! ! !
! ! !       ! Loop over LJ126 sites in molecule
! ! !       do i = 1, this%Molecule%NLJ126
! ! !         pLJ126 => this%Molecule%SiteLJ126(i)
! ! !         pLJ126%RX(n) = PXi
! ! !         pLJ126%RY(n) = PYi
! ! !         pLJ126%RZ(n) = PZi
! ! !       end do
! ! !
! ! !       ! Loop over charge sites in molecule
! ! !       if (LongRange .ne. RField) then
! ! !         do i = 1, this%Molecule%NCharge
! ! !           pCharge => this%Molecule%SiteCharge(i)
! ! !             pCharge%RX(n) = PXi
! ! !             pCharge%RY(n) = PYi
! ! !             pCharge%RZ(n) = PZi
! ! !         end do
! ! !       end if
! ! !
! ! !     end if
! ! !
! ! !   end subroutine TComponent_Mol2Atom1
! ! !
! ! !
! ! !
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
    type(TSiteLJ126), pointer      :: pLJ126
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j, nu

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    nu=this%Molecule%NUnit

    ! Check number of rotation axes
    if( this%Molecule%isElongated ) then

      ! Loop over molecules
      do i = 1, np
        ! Positions and quaternions of test particle i
        PX(i) = this%Pm0Test(i, 1)
        PY(i) = this%Pm0Test(i, 2)
        PZ(i) = this%Pm0Test(i, 3)
        q1 = this%Qm0Test(i, 1)
        q2 = this%Qm0Test(i, 2)
        q3 = this%Qm0Test(i, 3)
        q4 = this%Qm0Test(i, 4)

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

      ! Loop over LJ126 sites in molecule
      do j = 1, this%Molecule%NLJ126
        pLJ126 => this%Molecule%SiteLJ126(j)
        r1 = pLJ126%r(1) * BoxLengthInv
        r2 = pLJ126%r(2) * BoxLengthInv
        r3 = pLJ126%r(3) * BoxLengthInv
        do i = 1, np
          pLJ126%RXTest(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pLJ126%RYTest(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pLJ126%RZTest(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
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
          pQuadrupole%RXTest(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + &
&                                 r3 * A31(i)
          pQuadrupole%RYTest(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + &
&                                 r3 * A32(i)
          pQuadrupole%RZTest(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + &
&                                 r3 * A33(i)
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
          this%MueXTest(i,1) = mue1 * A11(i) + mue2 * A21(i) + mue3 * A31(i)
          this%MueYTest(i,1) = mue1 * A12(i) + mue2 * A22(i) + mue3 * A32(i)
          this%MueZTest(i,1) = mue1 * A13(i) + mue2 * A23(i) + mue3 * A33(i)
        end do
      end if

    else

      ! Loop over LJ126 sites in molecule
      do i = 1, this%Molecule%NLJ126
        pLJ126 => this%Molecule%SiteLJ126(i)
        do j = 1, np
          pLJ126%RXTest(j) = this%Pm0Test(j, 1)
          pLJ126%RYTest(j) = this%Pm0Test(j, 2)
          pLJ126%RZTest(j) = this%Pm0Test(j, 3)
        end do
      end do

      ! Loop over charge sites in molecule
      if (LongRange .ne. RField) then
        do i = 1, this%Molecule%NCharge
          pCharge => this%Molecule%SiteCharge(i)
          do j = 1, np
            pCharge%RX(j) = this%Pm0Test(j,1)
            pCharge%RY(j) = this%Pm0Test(j,2)
            pCharge%RZ(j) = this%Pm0Test(j,3)
          end do
        end do
      end if

    end if

  end subroutine TComponent_Mol2AtomTest
! ! !
! ! !
! ! !
! ! ! !==============================================================!
! ! ! !  Subroutine TComponent_Atom2Mol                              !
! ! ! !==============================================================!
! ! !
! ! !   subroutine TComponent_Atom2Mol( this, np )
! ! !
! ! !     implicit none
! ! !
! ! !     ! Include MPI header
! ! ! #if MPI_VER > 0
! ! !     include 'mpif.h'
! ! ! #endif
! ! !
! ! !     ! Declare arguments
! ! !     type(TComponent)    :: this
! ! !     integer, intent(in) :: np
! ! !
! ! !     ! Declare local variables
! ! !     real(RK)                       :: BoxLength
! ! !     real(RK)                       :: rx(np), ry(np), rz(np), r1x, r1y, r1z
! ! !     real(RK)                       :: q1(np), q2(np), q3(np), q4(np)
! ! !     real(RK)                       :: fx, fy, fz, tx, ty, tz
! ! !     real(RK)                       :: A11, A12, A13, A21, A22, A23, &
! ! ! &                                     A31, A32, A33
! ! !     type(TSiteLJ126), pointer      :: pLJ126
! ! !     type(TSiteCharge), pointer     :: pCharge
! ! !     type(TSiteDipole), pointer     :: pDipole
! ! !     type(TSiteQuadrupole), pointer :: pQuadrupole
! ! !     integer                        :: i, j
! ! !
! ! !     ! Assign local variables
! ! !     BoxLength = this%BoxLength
! ! !
! ! !     ! Initialize forces
! ! !     this%F(1:np, :) = 0._RK
! ! !
! ! !     ! Check number of rotation axes
! ! !     if( this%Molecule%isElongated ) then
! ! !
! ! !       ! Initialize torques
! ! !       this%T(1:np, :) = 0._RK
! ! !
! ! !       ! Initialize local arrays
! ! !       rx(:) = this%P0(:, 1)
! ! !       ry(:) = this%P0(:, 2)
! ! !       rz(:) = this%P0(:, 3)
! ! !       q1(:) = this%Q0(:, 1)
! ! !       q2(:) = this%Q0(:, 2)
! ! !       q3(:) = this%Q0(:, 3)
! ! !       q4(:) = this%Q0(:, 4)
! ! !
! ! !       ! Loop over LJ126 sites in molecule
! ! !       do j = 1, this%Molecule%NLJ126
! ! !         pLJ126 => this%Molecule%SiteLJ126(j)
! ! !         do i = 1, np
! ! !           fx = pLJ126%FX(i)
! ! !           fy = pLJ126%FY(i)
! ! !           fz = pLJ126%FZ(i)
! ! !           r1x = ( pLJ126%RX(i) - rx(i) ) * BoxLength
! ! !           r1y = ( pLJ126%RY(i) - ry(i) ) * BoxLength
! ! !           r1z = ( pLJ126%RZ(i) - rz(i) ) * BoxLength
! ! !           this%F(i, 1) = this%F(i, 1) + fx
! ! !           this%F(i, 2) = this%F(i, 2) + fy
! ! !           this%F(i, 3) = this%F(i, 3) + fz
! ! !           this%T(i, 1) = this%T(i, 1) + r1y * fz - r1z * fy
! ! !           this%T(i, 2) = this%T(i, 2) + r1z * fx - r1x * fz
! ! !           this%T(i, 3) = this%T(i, 3) + r1x * fy - r1y * fx
! ! !         end do
! ! !       end do
! ! !
! ! !       ! Loop over charge sites in molecule
! ! !       do j = 1, this%Molecule%NCharge
! ! !         pCharge => this%Molecule%SiteCharge(j)
! ! !         do i = 1, np
! ! !           fx = pCharge%FX(i)
! ! !           fy = pCharge%FY(i)
! ! !           fz = pCharge%FZ(i)
! ! !           r1x = ( pCharge%RX(i) - rx(i) ) * BoxLength
! ! !           r1y = ( pCharge%RY(i) - ry(i) ) * BoxLength
! ! !           r1z = ( pCharge%RZ(i) - rz(i) ) * BoxLength
! ! !           this%F(i, 1) = this%F(i, 1) + fx
! ! !           this%F(i, 2) = this%F(i, 2) + fy
! ! !           this%F(i, 3) = this%F(i, 3) + fz
! ! !           this%T(i, 1) = this%T(i, 1) + r1y * fz - r1z * fy
! ! !           this%T(i, 2) = this%T(i, 2) + r1z * fx - r1x * fz
! ! !           this%T(i, 3) = this%T(i, 3) + r1x * fy - r1y * fx
! ! !         end do
! ! !       end do
! ! !
! ! !       ! Loop over dipole sites in molecule
! ! !       do j = 1, this%Molecule%NDipole
! ! !         pDipole => this%Molecule%SiteDipole(j)
! ! !         do i = 1, np
! ! !           fx = pDipole%FX(i)
! ! !           fy = pDipole%FY(i)
! ! !           fz = pDipole%FZ(i)
! ! !           r1x = ( pDipole%RX(i) - rx(i) ) * BoxLength
! ! !           r1y = ( pDipole%RY(i) - ry(i) ) * BoxLength
! ! !           r1z = ( pDipole%RZ(i) - rz(i) ) * BoxLength
! ! !           this%F(i, 1) = this%F(i, 1) + fx
! ! !           this%F(i, 2) = this%F(i, 2) + fy
! ! !           this%F(i, 3) = this%F(i, 3) + fz
! ! !           this%T(i, 1) = this%T(i, 1) + pDipole%OY(i) * pDipole%TZ(i) &
! ! ! &                                     - pDipole%OZ(i) * pDipole%TY(i) &
! ! ! &                                     + r1y * fz - r1z * fy
! ! !           this%T(i, 2) = this%T(i, 2) + pDipole%OZ(i) * pDipole%TX(i) &
! ! ! &                                     - pDipole%OX(i) * pDipole%TZ(i) &
! ! ! &                                     + r1z * fx - r1x * fz
! ! !           this%T(i, 3) = this%T(i, 3) + pDipole%OX(i) * pDipole%TY(i) &
! ! ! &                                     - pDipole%OY(i) * pDipole%TX(i) &
! ! ! &                                     + r1x * fy - r1y * fx
! ! !         end do
! ! !       end do
! ! !
! ! !       ! Loop over quadrupole sites in molecule
! ! !       do j = 1, this%Molecule%NQuadrupole
! ! !         pQuadrupole => this%Molecule%SiteQuadrupole(j)
! ! !         do i = 1, np
! ! !           fx = pQuadrupole%FX(i)
! ! !           fy = pQuadrupole%FY(i)
! ! !           fz = pQuadrupole%FZ(i)
! ! !           r1x = ( pQuadrupole%RX(i) - rx(i) ) * BoxLength
! ! !           r1y = ( pQuadrupole%RY(i) - ry(i) ) * BoxLength
! ! !           r1z = ( pQuadrupole%RZ(i) - rz(i) ) * BoxLength
! ! !           this%F(i, 1) = this%F(i, 1) + fx
! ! !           this%F(i, 2) = this%F(i, 2) + fy
! ! !           this%F(i, 3) = this%F(i, 3) + fz
! ! !           this%T(i, 1) = this%T(i, 1) + pQuadrupole%OY(i) * pQuadrupole%TZ(i) &
! ! ! &                                     - pQuadrupole%OZ(i) * pQuadrupole%TY(i) &
! ! ! &                                     + r1y * fz - r1z * fy
! ! !           this%T(i, 2) = this%T(i, 2) + pQuadrupole%OZ(i) * pQuadrupole%TX(i) &
! ! ! &                                     - pQuadrupole%OX(i) * pQuadrupole%TZ(i) &
! ! ! &                                     + r1z * fx - r1x * fz
! ! !           this%T(i, 3) = this%T(i, 3) + pQuadrupole%OX(i) * pQuadrupole%TY(i) &
! ! ! &                                     - pQuadrupole%OY(i) * pQuadrupole%TX(i) &
! ! ! &                                     + r1x * fy - r1y * fx
! ! !         end do
! ! !       end do
! ! !
! ! !       do i = 1, np
! ! !         ! Add torques from reaction field
! ! !         tx = this%T(i, 1) + this%tRFX(i)
! ! !         ty = this%T(i, 2) + this%tRFY(i)
! ! !         tz = this%T(i, 3) + this%tRFZ(i)
! ! !
! ! !         ! Convert torque to body-fixed coordinates
! ! !         A11 = q1(i)**2 + q2(i)**2 - q3(i)**2 - q4(i)**2
! ! !         A12 = 2._RK * (q2(i) * q3(i) + q1(i) * q4(i))
! ! !         A13 = 2._RK * (q2(i) * q4(i) - q1(i) * q3(i))
! ! !         A21 = 2._RK * (q2(i) * q3(i) - q1(i) * q4(i))
! ! !         A22 = q1(i)**2 - q2(i)**2 + q3(i)**2 - q4(i)**2
! ! !         A23 = 2._RK * (q3(i) * q4(i) + q1(i) * q2(i))
! ! !         A31 = 2._RK * (q2(i) * q4(i) + q1(i) * q3(i))
! ! !         A32 = 2._RK * (q3(i) * q4(i) - q1(i) * q2(i))
! ! !         A33 = q1(i)**2 - q2(i)**2 - q3(i)**2 + q4(i)**2
! ! !         this%T(i, 1) = A11 * tx + A12 * ty + A13 * tz
! ! !         this%T(i, 2) = A21 * tx + A22 * ty + A23 * tz
! ! !         this%T(i, 3) = A31 * tx + A32 * ty + A33 * tz
! ! !       end do
! ! !
! ! !     else
! ! !
! ! !       ! Loop over LJ126 sites in molecule
! ! !       do j = 1, this%Molecule%NLJ126
! ! !         pLJ126 => this%Molecule%SiteLJ126(j)
! ! !         do i = 1, np
! ! !           this%F(i, 1) = this%F(i, 1) + pLJ126%FX(i)
! ! !           this%F(i, 2) = this%F(i, 2) + pLJ126%FY(i)
! ! !           this%F(i, 3) = this%F(i, 3) + pLJ126%FZ(i)
! ! !         end do
! ! !       end do
! ! !
! ! !       ! Loop over charge sites in molecule
! ! !       if (LongRange .ne. RField) then
! ! !         do j = 1, this%Molecule%NCharge
! ! !           pCharge => this%Molecule%SiteCharge(j)
! ! !           do i = 1, np
! ! !           this%F(i, 1) = this%F(i, 1) + pCharge%FX(i)
! ! !           this%F(i, 2) = this%F(i, 2) + pCharge%FY(i)
! ! !           this%F(i, 3) = this%F(i, 3) + pCharge%FZ(i)
! ! !           end do
! ! !         end do
! ! !       end if
! ! !
! ! !     end if
! ! !
! ! !     ! Reduce forces and torques from all processes
! ! ! #if MPI_VER > 0
! ! !     call MPI_Reduce( this%F(:, :), this%FAll(:, :), size( this%F ), &
! ! ! &     MPI_DOUBLE_PRECISION, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
! ! !     if( this%Molecule%isElongated ) &
! ! ! &     call MPI_Reduce( this%T(:, :), this%TAll(:, :), size( this%T ), &
! ! ! &       MPI_DOUBLE_PRECISION, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
! ! ! #endif
! ! !
! ! !   end subroutine TComponent_Atom2Mol
! ! !

!==============================================================!
!  Subroutine TComponent_Mol2Unit                              !
!==============================================================!

  subroutine TComponent_Mol2Unit( this, np, nu )

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
    real(RK)                       :: PmX(np), PmY(np), PmZ(np)
    real(RK)                       :: q1, q2, q3, q4, qinv
    real(RK)                       :: A11(np), A12(np), A13(np)
    real(RK)                       :: A21(np), A22(np), A23(np)
    real(RK)                       :: A31(np), A32(np), A33(np)
    integer                        :: nup
    real(RK)                       :: U11(nu), U12(nu), U13(nu)
    real(RK)                       :: U21(nu), U22(nu), U23(nu)
    real(RK)                       :: U31(nu), U32(nu), U33(nu)
    real(RK)                       :: UN11(nu), UN12(nu), UN13(nu)
    real(RK)                       :: UN21(nu), UN22(nu), UN23(nu)
    real(RK)                       :: UN31(nu), UN32(nu), UN33(nu)
    real(RK)                       :: UA11(nu*np), UA12(nu*np), UA13(nu*np)
    real(RK)                       :: UA21(nu*np), UA22(nu*np), UA23(nu*np)
    real(RK)                       :: UA31(nu*np), UA32(nu*np), UA33(nu*np)
    real(RK)                       :: ort(3,3), AUP(3,3)
    real(RK)                       :: E(3,3), C(3,3)
    real(RK)                       :: determinant
    real(RK)                       :: PX(nu*np), PY(nu*np), PZ(nu*np)
    real(RK)                       :: T, S, SInv
    real(RK)                       :: qu01, qu02, qu03, qu04
    real(RK)                       :: qu1, qu2, qu3, qu4, quinv
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    real(RK)                       :: mue1, mue2, mue3
    type(TSiteLJ126), pointer      :: pLJ126
    type(TSiteCharge), pointer     :: pCharge
!     type(TSiteDipole), pointer     :: pDipole
!     type(TSiteQuadrupole), pointer :: pQuadrupole
    type(TUnit), pointer           :: pUnit
    integer                        :: i, j, ij, k

    ! Broadcast positions and orientations to all processes
#if MPI_VER > 0
    call MPI_Bcast( this%Pm0(:, :), size( this%Pm0 ), &
&     MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
    if( this%Molecule%isElongated ) &
&     call MPI_Bcast( this%Qm0(:, :), size( this%Qm0 ), &
&       MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength
    nup = nu*np

    ! Check number of rotation axes
    if( this%Molecule%isElongated ) then
       ! Loop over molecules
       do i = 1, np
         ! Positions and quaternions of particle i
          PmX(i) = this%Pm0(i, 1)
          PmY(i) = this%Pm0(i, 2)
          PmZ(i) = this%Pm0(i, 3)
          q1 = this%Qm0(i, 1)
          q2 = this%Qm0(i, 2)
          q3 = this%Qm0(i, 3)
          q4 = this%Qm0(i, 4)

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
          this%Qm0(i, 1) = q1
          this%Qm0(i, 2) = q2
          this%Qm0(i, 3) = q3
          this%Qm0(i, 4) = q4

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

        ! Calculate initial COM position and quartenions for Units
        ! Loop over Units in molecule
        do j = 1, nu
           pUnit => this%Molecule%Unit(j)
           do i = 1, np
             this%P0(i, 1, j) = PmX(i) + (pUnit%P0(1)*A11(i)+pUnit%P0(2)*A21(i)+pUnit%P0(3)*A31(i)) * BoxLengthInv ! COM of Unit in space-fixed system
             this%P0(i, 2, j) = PmY(i) + (pUnit%P0(1)*A12(i)+pUnit%P0(2)*A22(i)+pUnit%P0(3)*A32(i)) * BoxLengthInv
             this%P0(i, 3, j) = PmZ(i) + (pUnit%P0(1)*A13(i)+pUnit%P0(2)*A23(i)+pUnit%P0(3)*A33(i)) * BoxLengthInv

             this%Q0(i,1,j) = this%Qm0(i,1)*pUnit%Q0(1) - this%Qm0(i,2)*pUnit%Q0(2) - &
&                              this%Qm0(i,3)*pUnit%Q0(3) - this%Qm0(i,4)*pUnit%Q0(4)
             this%Q0(i,2,j) = this%Qm0(i,1)*pUnit%Q0(2) + this%Qm0(i,2)*pUnit%Q0(1) + &
&                              this%Qm0(i,3)*pUnit%Q0(4) - this%Qm0(i,4)*pUnit%Q0(3)
             this%Q0(i,3,j) = this%Qm0(i,1)*pUnit%Q0(3) + this%Qm0(i,3)*pUnit%Q0(1) - &
&                              this%Qm0(i,2)*pUnit%Q0(4) + this%Qm0(i,4)*pUnit%Q0(2)
             this%Q0(i,4,j) = this%Qm0(i,1)*pUnit%Q0(4) + this%Qm0(i,4)*pUnit%Q0(1) - &
&                              this%Qm0(i,2)*pUnit%Q0(3) - this%Qm0(i,3)*pUnit%Q0(2)
           end do
         end do

      else    ! if Molecule is not Elongated
       do i = 1, np
         ! Positions and quaternions of particle i
          PmX(i) = this%Pm0(i, 1)
          PmY(i) = this%Pm0(i, 2)
          PmZ(i) = this%Pm0(i, 3)
          do j = 1, nu
            this%P0(i, 1, j) = PmX(i) ! COM of Unit in space-fixed system
            this%P0(i, 2, j) = PmY(i)
            this%P0(i, 3, j) = PmZ(i)
          end do
       end do
      end if

  end subroutine TComponent_Mol2Unit



!==============================================================!
!  Subroutine TComponent_Mol2Unit1                             !
!==============================================================!

  subroutine TComponent_Mol2Unit1( this, np, nu )

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
    real(RK)                       :: PmX, PmY, PmZ
    real(RK)                       :: q1, q2, q3, q4, qinv
    real(RK)                       :: A11, A12, A13
    real(RK)                       :: A21, A22, A23
    real(RK)                       :: A31, A32, A33
!    integer                        :: nu
    real(RK)                       :: U11(nu), U12(nu), U13(nu)
    real(RK)                       :: U21(nu), U22(nu), U23(nu)
    real(RK)                       :: U31(nu), U32(nu), U33(nu)
    real(RK)                       :: UA11(nu*np), UA12(nu*np), UA13(nu*np)
    real(RK)                       :: UA21(nu*np), UA22(nu*np), UA23(nu*np)
    real(RK)                       :: UA31(nu*np), UA32(nu*np), UA33(nu*np)
    real(RK)                       :: T, S, SInv
    real(RK)                       :: qu01, qu02, qu03, qu04
    real(RK)                       :: qu1, qu2, qu3, qu4, quinv
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    type(TUnit), pointer           :: pUnit
    integer                        :: i, j, ij, k

    ! Broadcast positions and orientations to all processes
#if MPI_VER > 0
    call MPI_Bcast( this%Pm0(:, :), size( this%Pm0 ), &
&     MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
    if( this%Molecule%isElongated ) &
&     call MPI_Bcast( this%Qm0(:, :), size( this%Qm0 ), &
&       MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength
!    nu  = this%Molecule%NUnit

    ! Check number of rotation axes
    if( this%Molecule%isElongated ) then
       ! Loop over molecules
         ! Positions and quaternions of particle i
          PmX = this%Pm0(np, 1)
          PmY = this%Pm0(np, 2)
          PmZ = this%Pm0(np, 3)
          q1 = this%Qm0(np, 1)
          q2 = this%Qm0(np, 2)
          q3 = this%Qm0(np, 3)
          q4 = this%Qm0(np, 4)

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
          this%Qm0(np, 1) = q1
          this%Qm0(np, 2) = q2
          this%Qm0(np, 3) = q3
          this%Qm0(np, 4) = q4

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


        ! Calculate initial COM position and quartenions for Units
        ! Loop over Units in molecule
        do j = 1, nu
           pUnit => this%Molecule%Unit(j)
           do i = 1, np
             this%P0(np, 1, j) = PmX + (pUnit%P0(1)*A11+pUnit%P0(2)*A21+pUnit%P0(3)*A31) * BoxLengthInv ! COM of Unit in space-fixed system
             this%P0(np, 2, j) = PmY + (pUnit%P0(1)*A12+pUnit%P0(2)*A22+pUnit%P0(3)*A32) * BoxLengthInv
             this%P0(np, 3, j) = PmZ + (pUnit%P0(1)*A13+pUnit%P0(2)*A23+pUnit%P0(3)*A33) * BoxLengthInv


             this%Q0(i,1,j) = this%Qm0(i,1)*pUnit%Q0(1) - this%Qm0(i,2)*pUnit%Q0(2) - &
&                              this%Qm0(i,3)*pUnit%Q0(3) - this%Qm0(i,4)*pUnit%Q0(4)
             this%Q0(i,2,j) = this%Qm0(i,1)*pUnit%Q0(2) + this%Qm0(i,2)*pUnit%Q0(1) + &
&                              this%Qm0(i,3)*pUnit%Q0(4) - this%Qm0(i,4)*pUnit%Q0(3)
             this%Q0(i,3,j) = this%Qm0(i,1)*pUnit%Q0(3) + this%Qm0(i,3)*pUnit%Q0(1) - &
&                              this%Qm0(i,2)*pUnit%Q0(4) + this%Qm0(i,4)*pUnit%Q0(2)
             this%Q0(i,4,j) = this%Qm0(i,1)*pUnit%Q0(4) + this%Qm0(i,4)*pUnit%Q0(1) - &
&                              this%Qm0(i,2)*pUnit%Q0(3) - this%Qm0(i,3)*pUnit%Q0(2)
           end do

        end do

    else    ! if Molecule is not Elongated
        ! Positions and quaternions of particle np
        PmX = this%Pm0(np, 1)
        PmY = this%Pm0(np, 2)
        PmZ = this%Pm0(np, 3)
        do j = 1, nu
          this%P0(np, 1, j) = PmX ! COM of Unit in space-fixed system
          this%P0(np, 2, j) = PmY
          this%P0(np, 3, j) = PmZ
        end do
    end if

  end subroutine TComponent_Mol2Unit1




!==============================================================!
!  Subroutine TComponent_Mol2Unit1Test                             !
!==============================================================!

  subroutine TComponent_Mol2Unit1Test( this, np, nu )

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
    real(RK)                       :: PmX, PmY, PmZ
    real(RK)                       :: q1, q2, q3, q4, qinv
    real(RK)                       :: A11, A12, A13
    real(RK)                       :: A21, A22, A23
    real(RK)                       :: A31, A32, A33
!    integer                        :: nu
    real(RK)                       :: U11(nu), U12(nu), U13(nu)
    real(RK)                       :: U21(nu), U22(nu), U23(nu)
    real(RK)                       :: U31(nu), U32(nu), U33(nu)
    real(RK)                       :: UA11(nu*np), UA12(nu*np), UA13(nu*np)
    real(RK)                       :: UA21(nu*np), UA22(nu*np), UA23(nu*np)
    real(RK)                       :: UA31(nu*np), UA32(nu*np), UA33(nu*np)
    real(RK)                       :: T, S, SInv
    real(RK)                       :: qu01, qu02, qu03, qu04
    real(RK)                       :: qu1, qu2, qu3, qu4, quinv
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    type(TUnit), pointer           :: pUnit
    integer                        :: i, j, ij, k

    ! Broadcast positions and orientations to all processes
#if MPI_VER > 0
    call MPI_Bcast( this%Pm0Test(:, :), size( this%Pm0Test ), &
&     MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
    if( this%Molecule%isElongated ) &
&     call MPI_Bcast( this%Qm0Test(:, :), size( this%Qm0Test ), &
&       MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength
!    nu  = this%Molecule%NUnit

    ! Check number of rotation axes
    if( this%Molecule%isElongated ) then
       ! Loop over molecules
         ! Positions and quaternions of particle i
          PmX = this%Pm0Test(np, 1)
          PmY = this%Pm0Test(np, 2)
          PmZ = this%Pm0Test(np, 3)
          q1 = this%Qm0Test(np, 1)
          q2 = this%Qm0Test(np, 2)
          q3 = this%Qm0Test(np, 3)
          q4 = this%Qm0Test(np, 4)

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
          this%Qm0Test(np, 1) = q1
          this%Qm0Test(np, 2) = q2
          this%Qm0Test(np, 3) = q3
          this%Qm0Test(np, 4) = q4

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


        ! Calculate initial COM position and quartenions for Units
        ! Loop over Units in molecule
        do j = 1, nu
           pUnit => this%Molecule%Unit(j)
           do i = 1, np
             this%P0Test(np, 1, j) = PmX + (pUnit%P0(1)*A11+pUnit%P0(2)*A21+pUnit%P0(3)*A31) * BoxLengthInv ! COM of Unit in space-fixed system
             this%P0Test(np, 2, j) = PmY + (pUnit%P0(1)*A12+pUnit%P0(2)*A22+pUnit%P0(3)*A32) * BoxLengthInv
             this%P0Test(np, 3, j) = PmZ + (pUnit%P0(1)*A13+pUnit%P0(2)*A23+pUnit%P0(3)*A33) * BoxLengthInv

             this%Q0Test(i,1,j) = this%Qm0(i,1)*pUnit%Q0(1) - this%Qm0(i,2)*pUnit%Q0(2) - &
&                                   this%Qm0(i,3)*pUnit%Q0(3) - this%Qm0(i,4)*pUnit%Q0(4)
             this%Q0Test(i,2,j) = this%Qm0(i,1)*pUnit%Q0(2) + this%Qm0(i,2)*pUnit%Q0(1) + &
&                                   this%Qm0(i,3)*pUnit%Q0(4) - this%Qm0(i,4)*pUnit%Q0(3)
             this%Q0Test(i,3,j) = this%Qm0(i,1)*pUnit%Q0(3) + this%Qm0(i,3)*pUnit%Q0(1) - &
&                                   this%Qm0(i,2)*pUnit%Q0(4) + this%Qm0(i,4)*pUnit%Q0(2)
             this%Q0Test(i,4,j) = this%Qm0(i,1)*pUnit%Q0(4) + this%Qm0(i,4)*pUnit%Q0(1) - &
&                                   this%Qm0(i,2)*pUnit%Q0(3) - this%Qm0(i,3)*pUnit%Q0(2)
           end do

        end do

    else    ! if Molecule is not Elongated
        ! Positions and quaternions of particle np
        PmX = this%Pm0Test(np, 1)
        PmY = this%Pm0Test(np, 2)
        PmZ = this%Pm0Test(np, 3)
        do j = 1, nu
          this%P0Test(np, 1, j) = PmX ! COM of Unit in space-fixed system
          this%P0Test(np, 2, j) = PmY
          this%P0Test(np, 3, j) = PmZ
        end do
    end if

  end subroutine TComponent_Mol2Unit1Test



!==============================================================!
!  Subroutine TComponent_Mol2Resize                            !
!==============================================================!

  subroutine TComponent_Mol2Resize( this, DelBoxFrac )

    implicit none

    ! Declare arguments
    real(RK),intent(in) :: DelBoxFrac
    type(TComponent)    :: this
!     integer, intent(in) :: np
!     integer, intent(in) :: nu

    ! Declare local variables
    real(RK)            :: BoxLengthInv
    integer             :: nu, np
    integer             :: i, j

    ! Calculate positions of units after global resize
    nu = this%Molecule%NUnit
    np = this%NPart
    if (nu .eq. 1) then
      do i=1, np
        this%P0(i,1,1) = this%Pm0(i,1)
        this%P0(i,2,1) = this%Pm0(i,2)
        this%P0(i,3,1) = this%Pm0(i,3)
      end do
    else
      do i=1, np
!         scalex = (DelBoxFac-1._RK ) * this%Pm0(i,1)
!         scaley = (DelBoxFac-1._RK )* this%Pm0(i,2)
!         scalez = (DelBoxFac-1._RK )* this%Pm0(i,3)

        ! New molecular center of masses
!         this%Pm0(i,1) = this%Pm0(i,1) * DelBoxFrac
!         this%Pm0(i,2) = this%Pm0(i,2) * DelBoxFrac
!         this%Pm0(i,3) = this%Pm0(i,3) * DelBoxFrac
        do j=1,nu
!           this%P0(i,1,j) = (this%P0(i,1,j) - this%Pm0(i,1) ) / DelBoxFrac + this%Pm0(i,1)
!           this%P0(i,2,j) = (this%P0(i,2,j) - this%Pm0(i,2) ) / DelBoxFrac + this%Pm0(i,2)
!           this%P0(i,3,j) = (this%P0(i,3,j) - this%Pm0(i,3) ) / DelBoxFrac + this%Pm0(i,3)
          this%P0(i,1,j) = (this%P0(i,1,j) + this%Pm0(i,1) ) - this%Pm0(i,1)/DelBoxFrac
          this%P0(i,2,j) = (this%P0(i,2,j) + this%Pm0(i,2) ) - this%Pm0(i,2)/DelBoxFrac 
          this%P0(i,3,j) = (this%P0(i,3,j) + this%Pm0(i,3) ) - this%Pm0(i,3)/DelBoxFrac
        end do
      end do
    end if

    end subroutine TComponent_Mol2Resize




!==============================================================!
!  Subroutine TComponent_Unit2Atom                             !
!==============================================================!

  subroutine TComponent_Unit2Atom( this, np, nu )

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
    integer                        :: nup
    real(RK)                       :: BoxLengthInv
    real(RK)                       :: PX(nu*np), PY(nu*np), PZ(nu*np)
    real(RK)                       :: q1, q2, q3, q4, qinv
    real(RK)                       :: A11(nu*np), A12(nu*np), A13(nu*np)
    real(RK)                       :: A21(nu*np), A22(nu*np), A23(nu*np)
    real(RK)                       :: A31(nu*np), A32(nu*np), A33(nu*np)
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    real(RK)                       :: mue1, mue2, mue3
    type(TSiteLJ126), pointer      :: pLJ126
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j, k, ik, jk

    ! Broadcast positions and orientations to all processes
#if MPI_VER > 0
    call MPI_Bcast( this%P0(:, :, :), size( this%P0 ), &
&     MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
    if( this%Molecule%isElongated ) &
&     call MPI_Bcast( this%Q0(:, :, :), size( this%Q0 ), &
&       MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength
    nup = nu*np

    if ( this%Molecule%isElongated ) then

      ! Loop over all units in Molecule
      do k = 1, nu
        ! Check number of rotation axes
!        if( this%Molecule%Unit(k)%isElongated ) then
!           print *, 'Unit is elongated'
          ! Loop over molecules
          do i = 1, np
            ik = (i-1)*nu+k
            ! Positions and quaternions of unit k in particle i
            PX(ik) = this%P0(i, 1, k)
            PY(ik) = this%P0(i, 2, k)
            PZ(ik) = this%P0(i, 3, k)
            q1 = this%Q0(i, 1, k)
            q2 = this%Q0(i, 2, k)
            q3 = this%Q0(i, 3, k)
            q4 = this%Q0(i, 4, k)

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
            this%Q0(i, 1, k) = q1
            this%Q0(i, 2, k) = q2
            this%Q0(i, 3, k) = q3
            this%Q0(i, 4, k) = q4

          ! Calculate rotation matrix elements
            A11(ik) = q1**2 + q2**2 - q3**2 - q4**2
            A12(ik) = 2._RK * (q2 * q3 + q1 * q4)
            A13(ik) = 2._RK * (q2 * q4 - q1 * q3)
            A21(ik) = 2._RK * (q2 * q3 - q1 * q4)
            A22(ik) = q1**2 - q2**2 + q3**2 - q4**2
            A23(ik) = 2._RK * (q3 * q4 + q1 * q2)
            A31(ik) = 2._RK * (q2 * q4 + q1 * q3)
            A32(ik) = 2._RK * (q3 * q4 - q1 * q2)
            A33(ik) = q1**2 - q2**2 - q3**2 + q4**2
          end do

          ! Loop over LJ126 sites in unit
          do j = 1, this%Molecule%Unit(k)%NLJ126
            pLJ126 => this%Molecule%Unit(k)%SiteLJ126(j)
            r1 = pLJ126%r(1) * BoxLengthInv
            r2 = pLJ126%r(2) * BoxLengthInv
            r3 = pLJ126%r(3) * BoxLengthInv
            do i = 1, np
              ik = (i-1)*nu+k
              pLJ126%RX(i) = PX(ik) + r1 * A11(ik) + r2 * A21(ik) + r3 * A31(ik)
              pLJ126%RY(i) = PY(ik) + r1 * A12(ik) + r2 * A22(ik) + r3 * A32(ik)
              pLJ126%RZ(i) = PZ(ik) + r1 * A13(ik) + r2 * A23(ik) + r3 * A33(ik)
            end do
          end do

          ! Loop over charge sites in molecule
          do j = 1, this%Molecule%Unit(k)%NCharge
            pCharge => this%Molecule%Unit(k)%SiteCharge(j)
            r1 = pCharge%r(1) * BoxLengthInv
            r2 = pCharge%r(2) * BoxLengthInv
            r3 = pCharge%r(3) * BoxLengthInv
            do i = 1, np
              ik = (i-1)*nu+k
              pCharge%RX(i) = PX(ik) + r1 * A11(ik) + r2 * A21(ik) + r3 * A31(ik)
              pCharge%RY(i) = PY(ik) + r1 * A12(ik) + r2 * A22(ik) + r3 * A32(ik)
              pCharge%RZ(i) = PZ(ik) + r1 * A13(ik) + r2 * A23(ik) + r3 * A33(ik)
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
            do i = 1, np
              ik = (i-1)*nu+k
              pDipole%RX(i) = PX(ik) + r1 * A11(ik) + r2 * A21(ik) + r3 * A31(ik)
              pDipole%RY(i) = PY(ik) + r1 * A12(ik) + r2 * A22(ik) + r3 * A32(ik)
              pDipole%RZ(i) = PZ(ik) + r1 * A13(ik) + r2 * A23(ik) + r3 * A33(ik)
              pDipole%OX(i) = or1 * A11(ik) + or2 * A21(ik) + or3 * A31(ik)
              pDipole%OY(i) = or1 * A12(ik) + or2 * A22(ik) + or3 * A32(ik)
              pDipole%OZ(i) = or1 * A13(ik) + or2 * A23(ik) + or3 * A33(ik)
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
            do i = 1, np
              ik = (i-1)*nu+k
              pQuadrupole%RX(i) = PX(ik) + r1 * A11(ik) + r2 * A21(ik) + r3 * A31(ik)
              pQuadrupole%RY(i) = PY(ik) + r1 * A12(ik) + r2 * A22(ik) + r3 * A32(ik)
              pQuadrupole%RZ(i) = PZ(ik) + r1 * A13(ik) + r2 * A23(ik) + r3 * A33(ik)
              pQuadrupole%OX(i) = or1 * A11(ik) + or2 * A21(ik) + or3 * A31(ik)
              pQuadrupole%OY(i) = or1 * A12(ik) + or2 * A22(ik) + or3 * A32(ik)
              pQuadrupole%OZ(i) = or1 * A13(ik) + or2 * A23(ik) + or3 * A33(ik)
            end do
          end do

        if( CutoffMode .eq. CenterofMass ) then
          mue1 = this%Molecule%Unit(k)%Mue(1)
          mue2 = this%Molecule%Unit(k)%Mue(2)
          mue3 = this%Molecule%Unit(k)%Mue(3)
         do i = 1, np
            ik = (i-1)*nu+k
            this%MueX(i, k) = mue1 * A11(ik) + mue2 * A21(ik) + mue3 * A31(ik)
            this%MueY(i, k) = mue1 * A12(ik) + mue2 * A22(ik) + mue3 * A32(ik)
            this%MueZ(i, k) = mue1 * A13(ik) + mue2 * A23(ik) + mue3 * A33(ik)
         end do
        end if

        end do
      else ! If molecule is not elongated
        do k = 1, nu
          ! Loop over LJ126 sites in molecule
         do i = 1, this%Molecule%Unit(k)%NLJ126
           pLJ126 => this%Molecule%Unit(k)%SiteLJ126(i)
           do j = 1, np
             pLJ126%RX(j) = this%P0(j, 1, k)
             pLJ126%RY(j) = this%P0(j, 2, k)
             pLJ126%RZ(j) = this%P0(j, 3, k)
           end do
         end do

         ! Loop over LJ126 sites in molecule
        do i = 1, this%Molecule%Unit(k)%NCharge
          pCharge => this%Molecule%Unit(k)%SiteCharge(i)
          do j = 1, np
            pCharge%RX(j) = this%P0(j, 1, k)
            pCharge%RY(j) = this%P0(j, 2, k)
            pCharge%RZ(j) = this%P0(j, 3, k)
          end do
        end do
      end do
    end if


  end subroutine TComponent_Unit2Atom


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
    real(RK)                       :: PX, PY, PZ
    real(RK)                       :: q1, q2, q3, q4, qinv
    real(RK)                       :: A11, A12, A13
    real(RK)                       :: A21, A22, A23
    real(RK)                       :: A31, A32, A33
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    real(RK)                       :: mue1, mue2, mue3
    type(TSiteLJ126), pointer      :: pLJ126
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j, k, ik, jk
    integer                        :: nu

    ! Broadcast positions and orientations to all processes
#if MPI_VER > 0
    call MPI_Bcast( this%P0(:, :, :), size( this%P0 ), &
&     MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
    if( this%Molecule%isElongated ) &
&     call MPI_Bcast( this%Q0(:, :, :), size( this%Q0 ), &
&       MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    if ( this%Molecule%isElongated ) then
      nu = this%Molecule%NUnit
      ! Loop over all units in Molecule
      do k = 1, nu
        ! Check number of rotation axes
!        if( this%Molecule%Unit(k)%isElongated ) then
!           print *, 'Unit is elongated'
            ik = (np-1)*nu+k
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
          do j = 1, this%Molecule%Unit(k)%NLJ126
            pLJ126 => this%Molecule%Unit(k)%SiteLJ126(j)
            r1 = pLJ126%r(1) * BoxLengthInv
            r2 = pLJ126%r(2) * BoxLengthInv
            r3 = pLJ126%r(3) * BoxLengthInv
            ik = (np-1)*nu+k
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
            ik = (np-1)*nu+k
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
            ik = (np-1)*nu+k
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
            ik = (np-1)*nu+k
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
            ik = (np-1)*nu+k
            this%MueX(np, k) = mue1 * A11 + mue2 * A21 + mue3 * A31
            this%MueY(np, k) = mue1 * A12 + mue2 * A22 + mue3 * A32
            this%MueZ(np, k) = mue1 * A13 + mue2 * A23 + mue3 * A33
        end if

        end do
      else ! If molecule is not elongated
        do k = 1, nu
          ! Loop over LJ126 sites in molecule
         do i = 1, this%Molecule%Unit(k)%NLJ126
           pLJ126 => this%Molecule%Unit(k)%SiteLJ126(i)
             pLJ126%RX(np) = this%P0(np, 1, k)
             pLJ126%RY(np) = this%P0(np, 2, k)
             pLJ126%RZ(np) = this%P0(np, 3, k)
         end do

         ! Loop over LJ126 sites in molecule
        do i = 1, this%Molecule%Unit(k)%NCharge
          pCharge => this%Molecule%Unit(k)%SiteCharge(i)
          pCharge%RX(np) = this%P0(np, 1, k)
          pCharge%RY(np) = this%P0(np, 2, k)
          pCharge%RZ(np) = this%P0(np, 3, k)
        end do
      end do
    end if


  end subroutine TComponent_Unit2Atom1Mol




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
    real(RK)                       :: PX, PY, PZ
    real(RK)                       :: q1, q2, q3, q4, qinv
    real(RK)                       :: A11, A12, A13
    real(RK)                       :: A21, A22, A23
    real(RK)                       :: A31, A32, A33
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    real(RK)                       :: mue1, mue2, mue3
    type(TSiteLJ126), pointer      :: pLJ126
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j, k, ik, jk
    integer                        :: nup

    ! Broadcast positions and orientations to all processes
#if MPI_VER > 0
    call MPI_Bcast( this%P0(:, :, :), size( this%P0 ), &
&     MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
    if( this%Molecule%isElongated ) &
&     call MPI_Bcast( this%Q0(:, :, :), size( this%Q0 ), &
&       MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    if ( this%Molecule%isElongated ) then
      nup = this%Molecule%NUnit
      ik = (np-1)*nup+nu
      ! Positions and quaternions of unit k in particle i
      PX = this%P0(np, 1, nu)
      PY = this%P0(np, 2, nu)
      PZ = this%P0(np, 3, nu)
      q1 = this%Q0(np, 1, nu)
      q2 = this%Q0(np, 2, nu)
      q3 = this%Q0(np, 3, nu)
      q4 = this%Q0(np, 4, nu)

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

          ! Loop over LJ126 sites in unit
      do j = 1, this%Molecule%Unit(nu)%NLJ126
        pLJ126 => this%Molecule%Unit(nu)%SiteLJ126(j)
        r1 = pLJ126%r(1) * BoxLengthInv
        r2 = pLJ126%r(2) * BoxLengthInv
        r3 = pLJ126%r(3) * BoxLengthInv
        ik = (np-1)*nup+nu
        pLJ126%RX(np) = PX + r1 * A11 + r2 * A21 + r3 * A31
        pLJ126%RY(np) = PY + r1 * A12 + r2 * A22 + r3 * A32
        pLJ126%RZ(np) = PZ + r1 * A13 + r2 * A23 + r3 * A33
      end do

          ! Loop over charge sites in molecule
      do j = 1, this%Molecule%Unit(nu)%NCharge
        pCharge => this%Molecule%Unit(nu)%SiteCharge(j)
        r1 = pCharge%r(1) * BoxLengthInv
        r2 = pCharge%r(2) * BoxLengthInv
        r3 = pCharge%r(3) * BoxLengthInv
        ik = (np-1)*nup+nu
        pCharge%RX(np) = PX + r1 * A11 + r2 * A21 + r3 * A31
        pCharge%RY(np) = PY + r1 * A12 + r2 * A22 + r3 * A32
        pCharge%RZ(np) = PZ + r1 * A13 + r2 * A23 + r3 * A33
      end do

          ! Loop over dipole sites in molecule
      do j = 1, this%Molecule%Unit(nu)%NDipole
        pDipole => this%Molecule%Unit(nu)%SiteDipole(j)
        r1 = pDipole%r(1) * BoxLengthInv
        r2 = pDipole%r(2) * BoxLengthInv
        r3 = pDipole%r(3) * BoxLengthInv
        or1 = pDipole%or(1)
        or2 = pDipole%or(2)
        or3 = pDipole%or(3)
        ik = (np-1)*nup+nu
        pDipole%RX(np) = PX + r1 * A11 + r2 * A21 + r3 * A31
        pDipole%RY(np) = PY + r1 * A12 + r2 * A22 + r3 * A32
        pDipole%RZ(np) = PZ + r1 * A13 + r2 * A23 + r3 * A33
        pDipole%OX(np) = or1 * A11 + or2 * A21 + or3 * A31
        pDipole%OY(np) = or1 * A12 + or2 * A22 + or3 * A32
        pDipole%OZ(np) = or1 * A13 + or2 * A23 + or3 * A33
      end do

      ! Loop over quadrupole sites in molecule
      do j = 1, this%Molecule%Unit(nu)%NQuadrupole
        pQuadrupole => this%Molecule%Unit(nu)%SiteQuadrupole(j)
        r1 = pQuadrupole%r(1) * BoxLengthInv
        r2 = pQuadrupole%r(2) * BoxLengthInv
        r3 = pQuadrupole%r(3) * BoxLengthInv
        or1 = pQuadrupole%or(1)
        or2 = pQuadrupole%or(2)
        or3 = pQuadrupole%or(3)
        ik = (np-1)*nup+nu
        pQuadrupole%RX(np) = PX + r1 * A11 + r2 * A21 + r3 * A31
        pQuadrupole%RY(np) = PY + r1 * A12 + r2 * A22 + r3 * A32
        pQuadrupole%RZ(np) = PZ + r1 * A13 + r2 * A23 + r3 * A33
        pQuadrupole%OX(np) = or1 * A11 + or2 * A21 + or3 * A31
        pQuadrupole%OY(np) = or1 * A12 + or2 * A22 + or3 * A32
        pQuadrupole%OZ(np) = or1 * A13 + or2 * A23 + or3 * A33
      end do

      if( CutoffMode .eq. CenterofMass ) then
        mue1 = this%Molecule%Unit(nu)%Mue(1)
        mue2 = this%Molecule%Unit(nu)%Mue(2)
        mue3 = this%Molecule%Unit(nu)%Mue(3)
        this%MueX(np, nu) = mue1 * A11 + mue2 * A21 + mue3 * A31
        this%MueY(np, nu) = mue1 * A12 + mue2 * A22 + mue3 * A32
        this%MueZ(np, nu) = mue1 * A13 + mue2 * A23 + mue3 * A33
      end if

    else ! If molecule is not elongated
       do k = 1, nu
         ! Loop over LJ126 sites in molecule
         do i = 1, this%Molecule%Unit(nu)%NLJ126
           pLJ126 => this%Molecule%Unit(nu)%SiteLJ126(i)
           pLJ126%RX(np) = this%P0(np, 1, nu)
           pLJ126%RY(np) = this%P0(np, 2, nu)
           pLJ126%RZ(np) = this%P0(np, 3, nu)
         end do

         ! Loop over LJ126 sites in molecule
         do i = 1, this%Molecule%Unit(nu)%NCharge
           pCharge => this%Molecule%Unit(nu)%SiteCharge(i)
           pCharge%RX(np) = this%P0(np, 1, nu)
           pCharge%RY(np) = this%P0(np, 2, nu)
           pCharge%RZ(np) = this%P0(np, 3, nu)
         end do
       end do
    end if


  end subroutine TComponent_Unit2Atom1



!==============================================================!
!  Subroutine TComponent_Unit2Atom1Test                        !
!==============================================================!

  subroutine TComponent_Unit2Atom1Test( this, np )

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
    real(RK)                       :: PX, PY, PZ
    real(RK)                       :: q1, q2, q3, q4, qinv
    real(RK)                       :: A11, A12, A13
    real(RK)                       :: A21, A22, A23
    real(RK)                       :: A31, A32, A33
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    real(RK)                       :: mue1, mue2, mue3
    type(TSiteLJ126), pointer      :: pLJ126
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j, k, ik, jk
    integer                        :: nu

    ! Broadcast positions and orientations to all processes
#if MPI_VER > 0
    call MPI_Bcast( this%P0(:, :, :), size( this%P0 ), &
&     MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
    if( this%Molecule%isElongated ) &
&     call MPI_Bcast( this%Q0(:, :, :), size( this%Q0 ), &
&       MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    if ( this%Molecule%isElongated ) then
      nu = this%Molecule%NUnit
      ! Loop over all units in Molecule
      do k = 1, nu
        ! Check number of rotation axes
!        if( this%Molecule%Unit(k)%isElongated ) then
!           print *, 'Unit is elongated'
            ik = (np-1)*nu+k
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
          do j = 1, this%Molecule%Unit(k)%NLJ126
            pLJ126 => this%Molecule%Unit(k)%SiteLJ126(j)
            r1 = pLJ126%r(1) * BoxLengthInv
            r2 = pLJ126%r(2) * BoxLengthInv
            r3 = pLJ126%r(3) * BoxLengthInv
            ik = (np-1)*nu+k
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
            ik = (np-1)*nu+k
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
            ik = (np-1)*nu+k
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
            ik = (np-1)*nu+k
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
            ik = (np-1)*nu+k
            this%MueX(np, k) = mue1 * A11 + mue2 * A21 + mue3 * A31
            this%MueY(np, k) = mue1 * A12 + mue2 * A22 + mue3 * A32
            this%MueZ(np, k) = mue1 * A13 + mue2 * A23 + mue3 * A33
        end if

        end do
      else ! If molecule is not elongated
        do k = 1, nu
          ! Loop over LJ126 sites in molecule
         do i = 1, this%Molecule%Unit(k)%NLJ126
           pLJ126 => this%Molecule%Unit(k)%SiteLJ126(i)
             pLJ126%RX(np) = this%P0(np, 1, k)
             pLJ126%RY(np) = this%P0(np, 2, k)
             pLJ126%RZ(np) = this%P0(np, 3, k)
         end do

         ! Loop over LJ126 sites in molecule
        do i = 1, this%Molecule%Unit(k)%NCharge
          pCharge => this%Molecule%Unit(k)%SiteCharge(i)
          pCharge%RX(np) = this%P0(np, 1, k)
          pCharge%RY(np) = this%P0(np, 2, k)
          pCharge%RZ(np) = this%P0(np, 3, k)
        end do
      end do
    end if

  end subroutine TComponent_Unit2Atom1Test




!==============================================================!
!  Subroutine TComponent_Atom2Unit                             !
!==============================================================!

  subroutine TComponent_Atom2Unit( this, np, nu )

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
    integer                        :: nup, neu
    real(RK)                       :: BoxLength
    real(RK)                       :: rx(np, nu), ry(np, nu), rz(np, nu), r1x, r1y, r1z
    real(RK)                       :: q1(np, nu), q2(np, nu), q3(np, nu), q4(np, nu)
    real(RK)                       :: fx, fy, fz, tx, ty, tz
    real(RK)                       :: A11, A12, A13, A21, A22, A23, &
&                                     A31, A32, A33
    type(TSiteLJ126), pointer      :: pLJ126
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j, k, ik

    ! Assign local variables
    BoxLength = this%BoxLength
    nup = nu*np
    neu = this%Molecule%NEUnit

    ! Initialize forces
    this%F(1:np, :, :) = 0._RK

    ! Loop over all Units in Molecule
    do k = 1, nu
      ! Check number of rotation axes
      if( this%Molecule%Unit(k)%isElongated ) then

        ! Initialize torques
        this%T(1:np, :, k) = 0._RK

         ! Initialize local arrays
         rx(:, k) = this%P0(:, 1, k)
         ry(:, k) = this%P0(:, 2, k)
         rz(:, k) = this%P0(:, 3, k)
         q1(:, k) = this%Q0(:, 1, k)
         q2(:, k) = this%Q0(:, 2, k)
         q3(:, k) = this%Q0(:, 3, k)
         q4(:, k) = this%Q0(:, 4, k)

         ! Loop over LJ126 sites in unit
         do j = 1, this%Molecule%Unit(k)%NLJ126
           pLJ126 => this%Molecule%Unit(k)%SiteLJ126(j)
           do i = 1, np
!             ik = (i-1)*nu+k
             fx = pLJ126%FX(i)
             fy = pLJ126%FY(i)
             fz = pLJ126%FZ(i)
             r1x = ( pLJ126%RX(i) - rx(i, k) ) * BoxLength
             r1y = ( pLJ126%RY(i) - ry(i, k) ) * BoxLength
             r1z = ( pLJ126%RZ(i) - rz(i, k) ) * BoxLength
             this%F(i, 1, k) = this%F(i, 1, k) + fx
             this%F(i, 2, k) = this%F(i, 2, k) + fy
             this%F(i, 3, k) = this%F(i, 3, k) + fz
             this%T(i, 1, k) = this%T(i, 1, k) + r1y * fz - r1z * fy
             this%T(i, 2, k) = this%T(i, 2, k) + r1z * fx - r1x * fz
             this%T(i, 3, k) = this%T(i, 3, k) + r1x * fy - r1y * fx
           end do
         end do

         ! Loop over charge sites in unit
         do j = 1, this%Molecule%Unit(k)%NCharge
           pCharge => this%Molecule%Unit(k)%SiteCharge(j)
           do i = 1, np
!             ik = (i-1)*nu+k
             fx = pCharge%FX(i)
             fy = pCharge%FY(i)
             fz = pCharge%FZ(i)
             r1x = ( pCharge%RX(i) - rx(i, k) ) * BoxLength
             r1y = ( pCharge%RY(i) - ry(i, k) ) * BoxLength
             r1z = ( pCharge%RZ(i) - rz(i, k) ) * BoxLength
             this%F(i, 1, k) = this%F(i, 1, k) + fx
             this%F(i, 2, k) = this%F(i, 2, k) + fy
             this%F(i, 3, k) = this%F(i, 3, k) + fz
             this%T(i, 1, k) = this%T(i, 1, k) + r1y * fz - r1z * fy
             this%T(i, 2, k) = this%T(i, 2, k) + r1z * fx - r1x * fz
             this%T(i, 3, k) = this%T(i, 3, k) + r1x * fy - r1y * fx
           end do
         end do

         ! Loop over dipole sites in unit
         do j = 1, this%Molecule%Unit(k)%NDipole
           pDipole => this%Molecule%Unit(k)%SiteDipole(j)
           do i = 1, np
!             ik = (i-1)*nu+k
             fx = pDipole%FX(i)
             fy = pDipole%FY(i)
             fz = pDipole%FZ(i)
             r1x = ( pDipole%RX(i) - rx(i, k) ) * BoxLength
             r1y = ( pDipole%RY(i) - ry(i, k) ) * BoxLength
             r1z = ( pDipole%RZ(i) - rz(i, k) ) * BoxLength
             this%F(i, 1, k) = this%F(i, 1, k) + fx
             this%F(i, 2, k) = this%F(i, 2, k) + fy
             this%F(i, 3, k) = this%F(i, 3, k) + fz
             this%T(i, 1, k) = this%T(i, 1, k) + pDipole%OY(i) * pDipole%TZ(i) &
&                                        - pDipole%OZ(i) * pDipole%TY(i) &
&                                        + r1y * fz - r1z * fy
             this%T(i, 2, k) = this%T(i, 2, k) + pDipole%OZ(i) * pDipole%TX(i) &
&                                        - pDipole%OX(i) * pDipole%TZ(i) &
&                                        + r1z * fx - r1x * fz
             this%T(i, 3, k) = this%T(i, 3, k) + pDipole%OX(i) * pDipole%TY(i) &
&                                        - pDipole%OY(i) * pDipole%TX(i) &
&                                        + r1x * fy - r1y * fx
           end do
         end do

         ! Loop over quadrupole sites in unit
         do j = 1, this%Molecule%Unit(k)%NQuadrupole
           pQuadrupole => this%Molecule%Unit(k)%SiteQuadrupole(j)
           do i = 1, np
             fx = pQuadrupole%FX(i)
             fy = pQuadrupole%FY(i)
             fz = pQuadrupole%FZ(i)
             r1x = ( pQuadrupole%RX(i) - rx(i, k) ) * BoxLength
             r1y = ( pQuadrupole%RY(i) - ry(i, k) ) * BoxLength
             r1z = ( pQuadrupole%RZ(i) - rz(i, k) ) * BoxLength
             this%F(i, 1, k) = this%F(i, 1, k) + fx
             this%F(i, 2, k) = this%F(i, 2, k) + fy
             this%F(i, 3, k) = this%F(i, 3, k) + fz
             this%T(i, 1, k) = this%T(i, 1, k) + pQuadrupole%OY(i) * pQuadrupole%TZ(i) &
&                                        - pQuadrupole%OZ(i) * pQuadrupole%TY(i) &
&                                        + r1y * fz - r1z * fy
             this%T(i, 2, k) = this%T(i, 2, k) + pQuadrupole%OZ(i) * pQuadrupole%TX(i) &
&                                        - pQuadrupole%OX(i) * pQuadrupole%TZ(i) &
&                                        + r1z * fx - r1x * fz
             this%T(i, 3, k) = this%T(i, 3, k) + pQuadrupole%OX(i) * pQuadrupole%TY(i) &
&                                        - pQuadrupole%OY(i) * pQuadrupole%TX(i) &
&                                        + r1x * fy - r1y * fx
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

        ! Loop over LJ126 sites in unit
         do j = 1, this%Molecule%Unit(k)%NLJ126
           pLJ126 => this%Molecule%Unit(k)%SiteLJ126(j)
           do i = 1, np
             ik = (i-1)*nu+k
             this%F(i, 1, k) = this%F(i, 1, k) + pLJ126%FX(i)
             this%F(i, 2, k) = this%F(i, 2, k) + pLJ126%FY(i)
             this%F(i, 3, k) = this%F(i, 3, k) + pLJ126%FZ(i)
           end do
         end do
       end if

    end do



       ! Reduce forces and torques from all processes
#if MPI_VER > 0
    call MPI_Reduce( this%F(:, :, :), this%FAll(:, :, :), size( this%F ), &
&     MPI_DOUBLE_PRECISION, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
    if( this%Molecule%isElongated ) &
&     call MPI_Reduce( this%T(:, :, :), this%TAll(:, :, :), size( this%T ), &
&       MPI_DOUBLE_PRECISION, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
#endif

  end subroutine TComponent_Atom2Unit


!==============================================================!
!  Subroutine TComponent_Unit2Mol                              !
!==============================================================!

  subroutine TComponent_Unit2Mol( this )

    implicit none

    ! Declare arguments
    type(TComponent)   :: this

    ! Declare local variables
    real(RK)                 :: mass
    real(RK)                 :: PX(this%NPart),PY(this%NPart),PZ(this%NPart)
    integer                  :: i
    integer                  :: np
    type(TMolecule), pointer :: pm

    np = this%NPart
    mass = 0._RK
    PX(1:np)   = 0._RK
    PY(1:np)   = 0._RK
    PZ(1:np)   = 0._RK

    do i=1,this%Molecule%NUnit
      mass = mass + this%Molecule%Unit(i)%Mass
      PX(1:np)   = PX(1:np)   + this%P0(1:np,1,i)*this%Molecule%Unit(i)%Mass
      PY(1:np)   = PY(1:np)   + this%P0(1:np,2,i)*this%Molecule%Unit(i)%Mass
      PZ(1:np)   = PZ(1:np)   + this%P0(1:np,3,i)*this%Molecule%Unit(i)%Mass
    end do

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
    real(RK)                 :: mass
    real(RK)                 :: PX,PY,PZ
    integer                  :: i
    type(TMolecule), pointer :: pm

    mass = 0._RK
    PX   = 0._RK
    PY   = 0._RK
    PZ   = 0._RK

!    pm => this%Molecule

    do i=1,this%Molecule%NUnit
      mass = mass + this%Molecule%Unit(i)%Mass
      PX   = PX   + this%P0(np,1,i)*this%Molecule%Unit(i)%Mass
      PY   = PY   + this%P0(np,2,i)*this%Molecule%Unit(i)%Mass
      PZ   = PZ   + this%P0(np,3,i)*this%Molecule%Unit(i)%Mass
    end do

    this%Pm0(np,1) = PX / mass
    this%Pm0(np,2) = PY / mass
    this%Pm0(np,3) = PZ / mass

  end subroutine TComponent_Unit2Mol1



!==============================================================!
!  Subroutine TComponent_Unit2Atom1 (per molecule)             !
!==============================================================!

  subroutine TComponent_Mol2UnitRotate( this, np, dq )

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
    real(RK)                       :: BoxLengthInv
    real(RK)                       :: PX, PY, PZ
    real(RK)                       :: A11, A12, A13
    real(RK)                       :: A21, A22, A23
    real(RK)                       :: A31, A32, A33
    real(RK)                       :: r1, r2, r3
    real(RK)                       :: q1, q2, q3,q4
    integer                        :: i, ik
    integer                        :: nup

    ! Broadcast positions and orientations to all processes
#if MPI_VER > 0
    call MPI_Bcast( this%P0(:, :, :), size( this%P0 ), &
&     MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
    if( this%Molecule%isElongated ) &
&     call MPI_Bcast( this%Q0(:, :, :), size( this%Q0 ), &
&       MPI_DOUBLE_PRECISION, NRootProc, MPI_COMM_WORLD, ierror )
#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    ! Calculate rotation matrix elements
    q1 = 1._RK
    q2 = dq(1)
    q3 = dq(2)
    q4 = dq(3)

    A11 = q2**2 - q3**2 - q4**2 + q1**2
    A12 = 2._RK * (q2 * q3 + q4*q1)
    A13 = 2._RK * (q2 * q4 - q3*q1)
    A21 = 2._RK * (q2 * q3 - q4*q1)
    A22 = - q2**2 + q3**2 - q4**2 + q1**2
    A23 = 2._RK * (q3 * q4 + q2*q1)
    A31 = 2._RK * (q2 * q4 + q3*q1)
    A32 = 2._RK * (q3 * q4 - q2*q1)
    A33 = - q2**2 - q3**2 + q4**2 + q1**2

    nup = this%Molecule%NUnit

    do i=1,nup
      ! Check number of rotation axes
      ik = (np-1)*nup+i
      ! Positions and quaternions of unit k in particle i
      PX = this%P0(np, 1, i)
      PY = this%P0(np, 2, i)
      PZ = this%P0(np, 3, i)
      
      q1 = this%Q0(np, 1, i)
      q2 = this%Q0(np, 2, i)
      q3 = this%Q0(np, 3, i)
      q4 = this%Q0(np, 4, i)

      ! Loop over LJ126 sites in unit
      r1 = (PX-this%Pm0(np,1)) * BoxLengthInv
      r2 = (PY-this%Pm0(np,2)) * BoxLengthInv
      r3 = (PZ-this%Pm0(np,3)) * BoxLengthInv

      this%P0(np,1,i) = this%Pm0(np,1) + r1 * A11 + r2 * A21 + r3 * A31
      this%P0(np,2,i) = this%Pm0(np,2) + r1 * A12 + r2 * A22 + r3 * A32
      this%P0(np,3,i) = this%Pm0(np,3) + r1 * A13 + r2 * A23 + r3 * A33
      
      this%Q0(np, 1, i) = q1 - dq(1) * q2 - dq(2) * q3 - dq(3) * q4
      this%Q0(np, 2, i) = q2 + dq(1) * q1 - dq(2) * q4 + dq(3) * q3
      this%Q0(np, 3, i) = q3 + dq(1) * q4 + dq(2) * q1 - dq(3) * q2
      this%Q0(np, 4, i) = q4 - dq(1) * q3 + dq(2) * q2 + dq(3) * q1

    end do

  end subroutine TComponent_Mol2UnitRotate






!==============================================================!
!  Subroutine TComponent_Flex2Rigid                            !
!==============================================================!

  subroutine TComponent_Flex2Rigid( this )

    implicit none

    ! Declare arguments
    type(TComponent)     :: this

    ! Declare local variables
!     type(TMolecule), pointer :: pm
    integer :: i
    
!     pm => this%Molecule
    
    do i=1, this%Molecule%NBond
      this%Molecule%IDFBond(i)%ForConst = this%Molecule%IDFBond(i)%ForConst * 1e10_RK
    end do
    do i=1, this%Molecule%NAngle
      this%Molecule%IDFAngle(i)%ForConst = this%Molecule%IDFAngle(i)%ForConst * 1e10_RK
    end do
    do i=1, this%Molecule%NDihedral
      this%Molecule%IDFDihedral(i)%ForConst = this%Molecule%IDFDihedral(i)%ForConst * 1e10_RK
    end do
    
  end subroutine TComponent_Flex2Rigid


!==============================================================!
!  Subroutine TComponent_Rigid2Flex                            !
!==============================================================!

  subroutine TComponent_Rigid2Flex( this )

    implicit none

    ! Declare arguments
    type(TComponent)     :: this

    ! Declare local variables
!     type(TMolecule), pointer :: pm
    integer :: i
    
!     pm => this%Molecule
    
    do i=1, this%Molecule%NBond
      this%Molecule%IDFBond(i)%ForConst = this%Molecule%IDFBond(i)%ForConst / 1e10_RK
    end do
    do i=1, this%Molecule%NAngle
      this%Molecule%IDFAngle(i)%ForConst = this%Molecule%IDFAngle(i)%ForConst / 1e10_RK
    end do
    do i=1, this%Molecule%NDihedral
      this%Molecule%IDFDihedral(i)%ForConst = this%Molecule%IDFDihedral(i)%ForConst / 1e10_RK
    end do
    
  end subroutine TComponent_Rigid2Flex



!==============================================================!
!  Subroutine TComponent_PredictGear                           !
!==============================================================!

  subroutine TComponent_PredictGear( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Declare local variables
    integer :: np, nu
    integer :: i, j, k

    ! Assign local variables
    np = this%NPart
!    nra = this%Molecule%NDFRot
    nu = this%Molecule%NUnit

    ! Predict COM positions and their derivatives
    do k=1,nu
     do j = 1, 3
      do i = 1, np
        this%P0(i, j, k) = this%P0(i, j, k) &
&                        + this%P1(i, j, k) &
&                        + this%P2(i, j, k) &
&                        + this%P3(i, j, k) &
&                        + this%P4(i, j, k) &
&                        + this%P5(i, j, k)
        this%P1(i, j, k) = this%P1(i, j, k) &
&                + 2._RK * this%P2(i, j, k) &
&                + 3._RK * this%P3(i, j, k) &
&                + 4._RK * this%P4(i, j, k) &
&                + 5._RK * this%P5(i, j, k)
        this%P2(i, j, k) = this%P2(i, j, k) &
&                + 3._RK * this%P3(i, j, k) &
&                + 6._RK * this%P4(i, j, k) &
&                +10._RK * this%P5(i, j, k)
        this%P3(i, j, k) = this%P3(i, j, k) &
&                + 4._RK * this%P4(i, j, k) &
&                +10._RK * this%P5(i, j, k)
        this%P4(i, j, k) = this%P4(i, j, k) &
&                + 5._RK * this%P5(i, j, k)
      end do
     end do

     if( this%Molecule%Unit(k)%IsElongated ) then

      ! Predict quaternion parameters and their derivatives
       do j = 1, 4
         do i = 1, np
           this%Q0(i, j, k) = this%Q0(i, j, k) &
&                           + this%Q1(i, j, k) &
&                           + this%Q2(i, j, k) &
&                           + this%Q3(i, j, k) &
&                           + this%Q4(i, j, k)
           this%Q1(i, j, k) = this%Q1(i, j, k) &
&                   + 2._RK * this%Q2(i, j, k) &
&                   + 3._RK * this%Q3(i, j, k) &
&                   + 4._RK * this%Q4(i, j, k)
           this%Q2(i, j, k) = this%Q2(i, j, k) &
&                   + 3._RK * this%Q3(i, j, k) &
&                   + 6._RK * this%Q4(i, j, k)
           this%Q3(i, j, k) = this%Q3(i, j, k) &
&                   + 4._RK * this%Q4(i, j, k)
         end do
       end do

      ! Predict angular velocities and their derivatives
       do j = 1, this%Molecule%Unit(k)%NDFRot
         do i = 1, np
           this%W0(i, j, k) = this%W0(i, j, k) &
&                           + this%W1(i, j, k) &
&                           + this%W2(i, j, k) &
&                           + this%W3(i, j, k) &
&                           + this%W4(i, j, k)
           this%W1(i, j, k) = this%W1(i, j, k) &
&                   + 2._RK * this%W2(i, j, k) &
&                   + 3._RK * this%W3(i, j, k) &
&                   + 4._RK * this%W4(i, j, k)
           this%W2(i, j, k) = this%W2(i, j, k) &
&                   + 3._RK * this%W3(i, j, k) &
&                   + 6._RK * this%W4(i, j, k)
           this%W3(i, j, k) = this%W3(i, j, k) &
&                   + 4._RK * this%W4(i, j, k)
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
    real(RK)          :: MassInv
    real(RK)          :: Moi23, Moi31, Moi12
    real(RK)          :: TMoi1, TMoi2, TMoi3
    real(RK), pointer :: pF(:, :, :), pT(:, :, :)
    integer           :: np, nra, nu
    integer           :: i, j, k
    real(RK)          :: r(3)

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength
    np = this%NPart
    nu = this%Molecule%NUnit

    ! Correct COM positions and their derivatives
#if MPI_VER > 0
    pF => this%FAll(:, :, :)
#else
    pF => this%F(:, :, :)
#endif
    do k= 1, nu
      MassInv = 1._RK / this%Molecule%Unit(k)%Mass
      do j = 1, 3
        do i = 1, np
          this%Corr0(i, j, k) = pF(i, j, k) &
&           * TimeStepSquared2 * BoxLengthInv * MassInv
          if( ConstantPressure .and. .not. NVTEquilibration ) &
&           this%Corr0(i, j, k) = this%Corr0(i, j, k) &
&             - this%P1(i, j, k) * dLogVolumeThird
          this%Corr1(i, j, k) = this%Corr0(i, j, k) - this%P2(i, j, k)
          this%P0(i, j, k) = this%P0(i, j, k) + this%Corr1(i, j, k) * Gear20
          this%P1(i, j, k) = this%P1(i, j, k) + this%Corr1(i, j, k) * Gear21
          this%P2(i, j, k) =                    this%Corr0(i, j, k)
          this%P3(i, j, k) = this%P3(i, j, k) + this%Corr1(i, j, k) * Gear23
          this%P4(i, j, k) = this%P4(i, j, k) + this%Corr1(i, j, k) * Gear24
          this%P5(i, j, k) = this%P5(i, j, k) + this%Corr1(i, j, k) * Gear25
        end do
      end do
    end do

    ! Calculate new positions of COM for molecules from new COM of units
    do i = 1, np
      r(:) = 0._RK
      do k= 1, nu
         r(:) = r(:) + this%Molecule%Unit(k)%Mass*this%P0(i,:,k)
      end do
      this%Pm0(i,:) = r(:)/this%Molecule%Mass
    end do

    ! Calculate displacement of molecules
    do i = 1, np
      do j = 1, 3
        this%Disp(i, j) = this%Disp(i, j) + this%Pm0(i, j) - this%Pm0old(i, j)

    ! Check for conservation of particles in primary cell
#if ARCH == 1
        if( this%Pm0(i, j) < -.5_RK ) then
          do k = 1, nu
             this%P0(i, j, k) = this%P0(i, j, k) + 1._RK
          end do
        elseif( this%Pm0(i, j) > .5_RK ) then
          do k = 1, nu
             this%P0(i, j, k) = this%P0(i, j, k) - 1._RK
          end do
        end if
#else
        do k = 1, nu
           this%P0(i, j, k) = this%P0(i, j, k) - anint( this%Pm0(i, j) )
        end do
      end do
    end do
#endif
    ! Calculate new positions of COM for molecules from new COM of units
    do i = 1, np
      r(:) = 0._RK
      do k= 1, nu
         r(:) = r(:) + this%Molecule%Unit(k)%Mass*this%P0(i,:,k)
      end do
      this%Pm0(i,:) = r(:)/this%Molecule%Mass
      this%Pm0old(i,:) = this%Pm0(i, :)
    end do


    if( this%Molecule%isElongated ) then

      ! Correct quaternion parameters and their derivatives
    do k = 1, nu
      do i = 1, np
        this%Corr0(i, 1, k) = TimeStep2 * ( - this%Q0(i, 2, k) * this%W0(i, 1, k) &
&                                        - this%Q0(i, 3, k) * this%W0(i, 2, k) &
&                                        - this%Q0(i, 4, k) * this%W0(i, 3, k))
        this%Corr0(i, 2, k) = TimeStep2 * ( + this%Q0(i, 1, k) * this%W0(i, 1, k) &
&                                        - this%Q0(i, 4, k) * this%W0(i, 2, k) &
&                                        + this%Q0(i, 3, k) * this%W0(i, 3, k))
        this%Corr0(i, 3, k) = TimeStep2 * ( + this%Q0(i, 4, k) * this%W0(i, 1, k) &
&                                        + this%Q0(i, 1, k) * this%W0(i, 2, k) &
&                                        - this%Q0(i, 2, k) * this%W0(i, 3, k))
        this%Corr0(i, 4, k) = TimeStep2 * ( - this%Q0(i, 3, k) * this%W0(i, 1, k) &
&                                        + this%Q0(i, 2, k) * this%W0(i, 2, k) &
&                                        + this%Q0(i, 1, k) * this%W0(i, 3, k))
      end do
      do j = 1, 4
        do i = 1, np
          this%Corr1(i, j, k) = this%Corr0(i, j, k) - this%Q1(i, j, k)
          this%Q0(i, j, k) = this%Q0(i, j, k) + this%Corr1(i, j, k) * Gear10
          this%Q1(i, j, k) =                 this%Corr0(i, j, k)
          this%Q2(i, j, k) = this%Q2(i, j, k) + this%Corr1(i, j, k) * Gear12
          this%Q3(i, j, k) = this%Q3(i, j, k) + this%Corr1(i, j, k) * Gear13
          this%Q4(i, j, k) = this%Q4(i, j, k) + this%Corr1(i, j, k) * Gear14
        end do
      end do
    end do
      ! Correct angular velocities and their derivatives
#if MPI_VER > 0
      pT => this%TAll(:, :, :)
#else
      pT => this%T(:, :, :)
#endif
    do k = 1, nu
      TMoi1 = TimeStep / this%Molecule%Unit(k)%MOI(1)
      TMoi2 = TimeStep / this%Molecule%Unit(k)%MOI(2)
      if( this%Molecule%Unit(k)%is3D ) then
        Moi23 = this%Molecule%Unit(k)%MOI(2) - this%Molecule%Unit(k)%MOI(3)
        Moi31 = this%Molecule%Unit(k)%MOI(3) - this%Molecule%Unit(k)%MOI(1)
        Moi12 = this%Molecule%Unit(k)%MOI(1) - this%Molecule%Unit(k)%MOI(2)
        TMoi3 = TimeStep / this%Molecule%Unit(k)%MOI(3)
        do i = 1, np
          this%Corr0(i, 1, k) = (pT(i, 1, k) + this%W0(i, 2, k) * this%W0(i, 3, k) * &
&                             Moi23) * TMoi1
          this%Corr0(i, 2, k) = (pT(i, 2, k) + this%W0(i, 3, k) * this%W0(i, 1, k) * &
&                             Moi31) * TMoi2
          this%Corr0(i, 3, k) = (pT(i, 3, k) + this%W0(i, 1, k) * this%W0(i, 2, k) * &
&                             Moi12) * TMoi3
        end do
      else
        do i = 1, np
          this%Corr0(i, 1, k) = pT(i, 1, k) * TMoi1
          this%Corr0(i, 2, k) = pT(i, 2, k) * TMoi2
        end do
      end if

      do j = 1, this%Molecule%Unit(k)%NDFRot
        do i = 1, np
          this%Corr1(i, j, k) = this%Corr0(i, j, k) - this%W1(i, j, k)
          this%W0(i, j, k) = this%W0(i, j, k) + this%Corr1(i, j, k) * Gear10
          this%W1(i, j, k) = this%Corr0(i, j, k)
          this%W2(i, j, k) = this%W2(i, j, k) + this%Corr1(i, j, k) * Gear12
          this%W3(i, j, k) = this%W3(i, j, k) + this%Corr1(i, j, k) * Gear13
          this%W4(i, j, k) = this%W4(i, j, k) + this%Corr1(i, j, k) * Gear14
        end do
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
    integer  :: np, nra, nu
    integer  :: i, j, k
    real(RK)          :: r(3)


    Korr = 2._RK - 1._RK / scale
    np = this%NPart
    nu = this%Molecule%NUnit

    do k = 1, nu
      do j = 1, 3
        do i = 1, np
          this%P1(i, j, k) = Korr * this%P1(i, j, k) + this%P2(i, j, k)
          this%P0(i, j, k) = this%P0(i, j, k) + this%P1(i, j, k)
        end do
      end do
    end do

    ! Calculate new positions of COM for molecules from new COM of units
    do i = 1, np
      r(:) = 0._RK
      do k= 1, nu
         r(:) = r(:) + this%Molecule%Unit(k)%Mass*this%P0(i,:,k)
      end do
      this%Pm0(i,:) = r(:)/this%Molecule%Mass
    end do

    ! Calculate displacement of molecules
    do i = 1, np
      do j = 1, 3
        this%Disp(i, j) = this%Disp(i, j) + this%Pm0(i, j) - this%Pm0old(i, j)

        ! Check for conservation of particles in primary cell
#if ARCH == 1
        if( this%Pm0(i, j) < -.5_RK ) then
          do k = 1, nu
             this%P0(i, j, k) = this%P0(i, j, k) + 1._RK
          end do
        elseif( this%Pm0(i, j) > .5_RK ) then
          do k = 1, nu
             this%P0(i, j, k) = this%P0(i, j, k) - 1._RK
          end do
        end if
#else
        do k = 1, nu
           this%P0(i, j, k) = this%P0(i, j, k) - anint( this%Pm0(i, j) )
        end do
      end do
    end do
#endif
    ! Calculate new positions of COM for molecules from new COM of units
    do i = 1, np
      r(:) = 0._RK
      do k= 1, nu
         r(:) = r(:) + this%Molecule%Unit(k)%Mass*this%P0(i,:,k)
      end do
      this%Pm0(i,:) = r(:)/this%Molecule%Mass
      this%Pm0old(i,:) = this%Pm0(i, :)
    end do

    if( this%Molecule%IsElongated ) then
      do k = 1, nu
        nra = this%Molecule%Unit(k)%NDFRot
        do i = 1, np
          do j = 1, 4
            this%Q0tmp(i, j, k) = this%Q0(i, j, k) + .5_RK * this%Q1(i, j, k)
          end do
          do j = 1, nra
            this%W0(i, j, k) = Korr * this%W0(i, j, k) + .5_RK * this%W1(i, j, k)
          end do
          this%Q1(i, 1, k) = TimeStep2 * ( - this%Q0tmp(i, 2, k) * this%W0(i, 1, k) &
&                                       - this%Q0tmp(i, 3, k) * this%W0(i, 2, k) &
&                                       - this%Q0tmp(i, 4, k) * this%W0(i, 3, k))
          this%Q1(i, 2, k) = TimeStep2 * ( + this%Q0tmp(i, 1, k) * this%W0(i, 1, k) &
&                                       - this%Q0tmp(i, 4, k) * this%W0(i, 2, k) &
&                                       + this%Q0tmp(i, 3, k) * this%W0(i, 3, k))
          this%Q1(i, 3, k) = TimeStep2 * ( + this%Q0tmp(i, 4, k) * this%W0(i, 1, k) &
&                                       + this%Q0tmp(i, 1, k) * this%W0(i, 2, k) &
&                                       - this%Q0tmp(i, 2, k) * this%W0(i, 3, k))
          this%Q1(i, 4, k) = TimeStep2 * ( - this%Q0tmp(i, 3, k) * this%W0(i, 1, k) &
&                                       + this%Q0tmp(i, 2, k) * this%W0(i, 2, k) &
&                                       + this%Q0tmp(i, 1, k) * this%W0(i, 3, k))
          do j = 1, 4
            this%Q0(i, j, k) = this%Q0(i, j, k) + this%Q1(i, j, k)
          end do
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
    real(RK), pointer :: pF(:, :, :), pT(:, :, :)
    real(RK)          :: BoxLengthInv, MassInv
    real(RK)          :: Moi23, Moi31, Moi12
    real(RK)          :: TMoi1, TMoi2, TMoi3
    integer           :: np, nra, nu
    integer           :: i, j, k

    BoxLengthInv = 1._RK / this%BoxLength
    np = this%NPart
    nu = this%Molecule%NUnit
#if MPI_VER > 0
    pF => this%FAll(:, :, :)
#else
    pF => this%F(:, :, :)
#endif
   do k = 1, nu
     MassInv = 1._RK / this%Molecule%Unit(k)%Mass
     do j = 1, 3
       do i = 1, np
         this%P2(i, j, k) = pF(i, j, k) * TimeStepSquared2 * BoxLengthInv * MassInv &
&                       - (this%P1(i, j, k) + this%P2(i, j, k)) * dLogVolumeThird
         this%P1(i, j, k) = this%P1(i, j, k) + this%P2(i, j, k)
       end do
     end do
   end do

    if( this%Molecule%IsElongated ) then
       do k = 1, nu
         nra = this%Molecule%Unit(k)%NDFRot
         TMoi1 = TimeStep / this%Molecule%Unit(k)%MOI(1)
         TMoi2 = TimeStep / this%Molecule%Unit(k)%MOI(2)
#if MPI_VER > 0
      pT => this%TAll(:, :, :)
#else
      pT => this%T(:, :, :)
#endif
      if( this%Molecule%Unit(k)%is3D ) then
        TMoi3 = TimeStep / this%Molecule%Unit(k)%MOI(3)
        Moi23 = this%Molecule%Unit(k)%MOI(2) - this%Molecule%Unit(k)%MOI(3)
        Moi31 = this%Molecule%Unit(k)%MOI(3) - this%Molecule%Unit(k)%MOI(1)
        Moi12 = this%Molecule%Unit(k)%MOI(1) - this%Molecule%Unit(k)%MOI(2)
        do i = 1, np
          this%W1(i, 1, k) = (pT(i, 1, k) + this%W0(i, 2, k) * this%W0(i, 3, k) * Moi23) * &
&                         TMoi1
          this%W1(i, 2, k) = (pT(i, 2, k) + this%W0(i, 3, k) * this%W0(i, 1, k) * Moi31) * &
&                         TMoi2
          this%W1(i, 3, k) = (pT(i, 3, k) + this%W0(i, 1, k) * this%W0(i, 2, k) * Moi12) * &
&                         TMoi3
        end do
      else
        do i = 1, np
          this%W1(i, 1, k) = pT(i, 1, k) * TMoi1
          this%W1(i, 2, k) = pT(i, 2, k) * TMoi2
        end do
      end if

      do j = 1, nra
        do i = 1, np
          this%W0(i, j, k) = this%W0(i, j, k) + .5_RK * this%W1(i, j, k)
        end do
      end do
      do i = 1, np
        this%Q1(i, 1, k) = TimeStep2 * ( - this%Q0(i, 2, k) * this%W0(i, 1, k) &
&                                     - this%Q0(i, 3, k) * this%W0(i, 2, k) &
&                                     - this%Q0(i, 4, k) * this%W0(i, 3, k))
        this%Q1(i, 2, k) = TimeStep2 * ( + this%Q0(i, 1, k) * this%W0(i, 1, k) &
&                                     - this%Q0(i, 4, k) * this%W0(i, 2, k) &
&                                     + this%Q0(i, 3, k) * this%W0(i, 3, k))
        this%Q1(i, 3, k) = TimeStep2 * ( + this%Q0(i, 4, k) * this%W0(i, 1, k) &
&                                     + this%Q0(i, 1, k) * this%W0(i, 2, k) &
&                                     - this%Q0(i, 2, k) * this%W0(i, 3, k))
        this%Q1(i, 4, k) = TimeStep2 * ( - this%Q0(i, 3, k) * this%W0(i, 1, k) &
&                                     + this%Q0(i, 2, k) * this%W0(i, 2, k) &
&                                     + this%Q0(i, 1, k) * this%W0(i, 3, k))
      end do
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

    if (UseIntDegFreed) then 
      if( this%NMoveMolSuccesses < this%NMoveMolAttempts * Acceptance ) then
        this%DispMolTran = this%DispMolTran * .95_RK
      else if( this%DispMolTran < DispMolTranLimit ) then
        this%DispTran = this%DispTran * 1.05_RK
      end if

      ! Update rotational displacement
      if( this%NRotateSuccesses < this%NRotateAttempts * Acceptance ) then
        this%DispMolRot = this%DispMolRot * .95_RK
      else if( this%DispMolRot < DispMolRotLimit ) then
        this%DispMolRot = this%DispMolRot * 1.05_RK
      end if
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
    real(RK), intent(in), optional :: q(4, 1:this%Molecule%NUnit)

    integer :: k

    ! Test boundaries of particle arrays
    if( this%NPart >= this%NPartMax ) then
      tooManyParticles = .true.
      return
    end if

    ! Increase NPart
    this%NPart = this%NPart + 1
#if MPI_VER > 0
    this%NPart1 = ProcRange( this%NPart, this%NPart0, this%NPart2 )
!     this%NPart1 = 1 + (this%NPart - 1) / NProcs
!     this%NPart0 = 1 + this%NPart1 * NProc
!     this%NPart2 = min( this%NPart0 + this%NPart1 - 1, this%NPart )
!     this%NPart1 = this%NPart2 - this%NPart0 + 1
#endif

    ! Set coordinates and orientation of new particle
    if ( this%Molecule%NUnit .ne. 1) then
        write(*,*) 'Adding particles not supported for flexible molecules!'
        STOP
    end if

    ! Add Particle
    this%Pm0(this%NPart, :) = r(:)
    this%P0(this%NPart, :,1) = r(:)
    if( this%Molecule%isElongated ) then
!       this%Qm0(this%NPart, :) = q(:,1)
      this%Q0(this%NPart, :,:) = q(:,:)
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

    ! Declare local variables
    integer :: k, nu

    if( np .ne. this%NPart ) then

      ! Copy coordinates and orientation of last particle
      this%Pm0(np, :) = this%Pm0(this%NPart, :)
      if( this%Molecule%isElongated ) then
        this%Qm0(np, :) = this%Qm0(this%NPart, :)
      end if


      nu=this%Molecule%NUnit

      ! Calculate Unit / Atom positions
      call Mol2Unit1 ( this, np, nu )
      do k=1,nu
        call Unit2Atom1( this, np, k )
      end do

    end if

    ! Remove last particle
    this%NPart = this%NPart - 1
#if MPI_VER > 0
    this%NPart1 = 1 + (this%NPart - 1) / NProcs
    this%NPart0 = 1 + this%NPart1 * NProc
    this%NPart2 = min( this%NPart0 + this%NPart1 - 1, this%NPart )
    this%NPart1 = this%NPart2 - this%NPart0 + 1
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
    integer          :: k

    ! Save current state
    this%P0Save = this%P0
    do k = 1, this%Molecule%NUnit
      if( this%Molecule%Unit(k)%IsElongated ) then
        this%Q0Save = this%Q0
        exit
      end if
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
      if( this%Molecule%Unit(i)%IsElongated ) then
        this%Q0 = this%Q0Save
        exit
      end if
    end do

    ! Calculate site positions
    do i=1,this%Molecule%NUnit
      call Unit2Atom( this, this%NPart,1 )
    end do

  end subroutine TComponent_RestoreState



!==============================================================!
!  Subroutine TComponent_RestartSave                           !
!==============================================================!

  subroutine TComponent_RestartSave( this )

    implicit none

    ! Declare arguments
    type(TComponent) :: this

    ! Declare local variables
    integer :: np, i
    integer :: k,nu

    ! Assign local variables
    np = this%NPart
    nu = this%Molecule%NUnit

    ! Check for root process
    if( .not. RootProc ) return

    ! Save contents to restart file
    write( iounit_restart, '(I10)' ) np
    write( iounit_restart, '(I10)' ) nu


    ! Centers of mass positions for units
    do i = 1, np
      do k = 1, nu
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P0( i, :, k )
      end do
    end do

    if( SimulationType .eq. MolecularDynamics ) then
      ! Centers of mass positions' derivatives
      do i = 1, np
        do k = 1, nu
          write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P1( i, :, k )
        end do
      end do
      do i = 1, np
        do k = 1, nu
          write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P2( i, :, k )
        end do
      end do

      if( IntegratorType .eq. IntegratorTypeGear ) then
        do i = 1, np
          do k = 1, nu
             write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P3( i, :, k )
          end do
        end do
        do i = 1, np
          do k = 1, nu
             write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P4( i, :, k )
          end do
        end do
        do i = 1, np
          do k = 1, nu
             write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P5( i, :, k )
          end do
        end do
      end if
    else
      write( iounit_restart, '(ES20.12E3)' ) this%DispTran
      write( iounit_restart, '(2I10)' ) this%NMoveAttempts, this%NMoveSuccesses
      write( iounit_restart, '(2I10)' ) this%NMoveBiasedAttempts, &
&       this%NMoveBiasedSuccesses
      if ( UseIntDegFreed ) then
        write( iounit_restart, '(ES20.12E3)' ) this%DispMolTran
        write( iounit_restart, '(2I10)' ) this%NMoveMolAttempts, this%NMoveMolSuccesses
        write( iounit_restart, '(2I10)' ) this%NMoveBiasedMolAttempts, &
&         this%NMoveBiasedMolSuccesses
      end if
    end if

    if( this%Molecule%isElongated ) then
      ! Quaternion parameters
      do i = 1, np
         do  k = 1, nu
            write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q0( i, :, k )
         end do
      end do

      if( SimulationType .eq. MolecularDynamics ) then
        ! Quaternion parameters' derivatives
        do i = 1, np
         do  k = 1, nu
            write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q1( i, :, k )
         end do
        end do

        if( IntegratorType .eq. IntegratorTypeGear ) then
          do i = 1, np
            do  k = 1, nu
              write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q2( i, :, k )
            end do
          end do
          do i = 1, np
            do  k = 1, nu
              write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q3( i, :, k )
            end do
          end do
          do i = 1, np
            do  k = 1, nu
              write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q4( i, :, k )
            end do
          end do
        end if

        ! Angular velocities and their derivatives
        do i = 1, np
          do k = 1, nu
            write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W0( i, :, k )
          end do
        end do
        do i = 1, np
          do k = 1, nu
            write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W1( i, :, k )
          end do
        end do

        if( IntegratorType .eq. IntegratorTypeGear ) then
          do i = 1, np
            do k= 1, nu
              write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W2( i, :, k )
            end do
          end do
          do i = 1, np
            do k= 1, nu
              write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W3( i, :, k )
            end do
          end do
          do i = 1, np
            do k= 1, nu
              write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W4( i, :, k )
            end do
          end do
        end if
      else
        write( iounit_restart, '(ES20.12E3)' ) this%DispRot
        write( iounit_restart, '(2I10)' ) this%NRotateAttempts, &
&         this%NRotateSuccesses
        write( iounit_restart, '(2I10)' ) this%NRotateBiasedAttempts, &
&         this%NRotateBiasedSuccesses
        if ( UseIntDegFreed ) then
          write( iounit_restart, '(ES20.12E3)' ) this%DispMolRot
          write( iounit_restart, '(2I10)' ) this%NRotateMolAttempts, &
&           this%NRotateMolSuccesses
          write( iounit_restart, '(2I10)' ) this%NRotateBiasedMolAttempts, &
&           this%NRotateBiasedMolSuccesses
        end if
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
    integer :: i, np
    integer :: nu, k
    real(RK):: r(3)

    if( RootProc ) then

      ! Read contents from restart file
      read( iounit_restart, '(I10)' ) np
      if( np > this%NPartMax ) &
&       call Error( 'Not enough memory to read particles from restart file' )
      this%NPart = np
      read( iounit_restart, '(I10)' ) nu
      this%Molecule%NUnit = nu

      ! Centers of mass positions
      do i = 1, np
        do k = 1, nu
          read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P0( i, :, k )
        end do
      end do

      ! Calculate positions of COM for molecules from  COM of units
      do i = 1, np
        r(:) = 0._RK
        do k= 1, nu
           r(:) = r(:) + this%Molecule%Unit(k)%Mass*this%P0(i,:,k)
        end do
        this%Pm0(i,:) = r(:)/this%Molecule%Mass
        this%Pm0old(i,:) = this%Pm0(i, :)
      end do


      if( SimulationType .eq. MolecularDynamics ) then
        ! Centers of mass positions' derivatives
        do i = 1, np
          do k = 1, nu
            read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P1( i, : , k )
          end do
        end do
        do i = 1, np
          do k = 1, nu
            read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P2( i, : , k )
          end do
        end do

        if( IntegratorType .eq. IntegratorTypeGear ) then
          do i = 1, np
            do k = 1, nu
              read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P3( i, :, k )
            end do
          end do
          do i = 1, np
            do k = 1, nu
              read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P4( i, :, k )
            end do
          end do
          do i = 1, np
            do k = 1, nu
              read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P5( i, :, k )
            end do
          end do
        end if
      else
        read( iounit_restart, '(ES20.12E3)' ) this%DispTran
        read( iounit_restart, '(2I10)' ) this%NMoveAttempts, this%NMoveSuccesses
        read( iounit_restart, '(2I10)' ) this%NMoveBiasedAttempts, &
&         this%NMoveBiasedSuccesses
        if ( UseIntDegFreed ) then
          read( iounit_restart, '(ES20.12E3)' ) this%DispMolTran
          read( iounit_restart, '(2I10)' ) this%NMoveMolAttempts, this%NMoveMolSuccesses
          read( iounit_restart, '(2I10)' ) this%NMoveBiasedMolAttempts, &
&           this%NMoveBiasedMolSuccesses
        end if
      end if

      if( this%Molecule%isElongated ) then
        ! Quaternion parameters
        do i = 1, np
          do k = 1, nu
            read( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q0( i, :, k )
          end do
        end do

        if( SimulationType .eq. MolecularDynamics ) then
          ! Quaternion parameters' derivatives
          do i = 1, np
            do k = 1, nu
              read( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q1( i, :, k )
            end do
          end do

          if( IntegratorType .eq. IntegratorTypeGear ) then
            do i = 1, np
              do k = 1, nu
                read( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q2( i, :, k )
              end do
            end do
            do i = 1, np
              do k = 1, nu
                read( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q3( i, :, k )
              end do
            end do
            do i = 1, np
              do k = 1, nu
                read( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q4( i, :, k )
              end do
            end do
          end if

          ! Angular velocities and their derivatives
          do i = 1, np
            do k = 1, nu
              read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W0( i, :, k )
            end do
          end do
          do i = 1, np
            do k = 1, nu
              read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W1( i, :, k )
            end do
          end do

          if( IntegratorType .eq. IntegratorTypeGear ) then
            do i = 1, np
              do k = 1, nu
                read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W2( i, : , k)
              end do
            end do
            do i = 1, np
              do k = 1, nu
                read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W3( i, : , k)
              end do
            end do
            do i = 1, np
              do k = 1, nu
                read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W4( i, : , k)
              end do
            end do
          end if
        else
          read( iounit_restart, '(ES20.12E3)' ) this%DispRot
          read( iounit_restart, '(2I10)' ) this%NRotateAttempts, &
&           this%NRotateSuccesses
          read( iounit_restart, '(2I10)' ) this%NRotateBiasedAttempts, &
&           this%NRotateBiasedSuccesses
          if ( UseIntDegFreed ) then
            read( iounit_restart, '(ES20.12E3)' ) this%DispMolRot
            read( iounit_restart, '(2I10)' ) this%NRotateMolAttempts, &
&             this%NRotateMolSuccesses
            read( iounit_restart, '(2I10)' ) this%NRotateBiasedMolAttempts, &
&             this%NRotateBiasedMolSuccesses
          end if

        end if
      end if

      if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
        read( iounit_restart, '(ES20.12E3)' ) this%WF(:)
        read( iounit_restart, '(I10)' ) this%NState(:)
        read( iounit_restart, '(I10)' ) this%NStateWF(:)
      end if

    end if

#if MPI_VER > 0
    call MPI_Bcast( this%NPart, 1, MPI_INTEGER, NRootProc, &
&     MPI_COMM_WORLD, ierror )
    call MPI_Bcast( this%P0(:, :, :), size( this%P0 ), MPI_DOUBLE_PRECISION, &
&     NRootProc, MPI_COMM_WORLD, ierror )
    call MPI_Bcast( this%Pm0(:, :), size( this%Pm0 ), MPI_DOUBLE_PRECISION, &
&     NRootProc, MPI_COMM_WORLD, ierror )
    if( this%Molecule%isElongated ) &
&     call MPI_Bcast( this%Q0(:, :, :), size( this%Q0 ), MPI_DOUBLE_PRECISION, &
&       NRootProc, MPI_COMM_WORLD, ierror )
    if( this%Molecule%isElongated ) &
&     call MPI_Bcast( this%Qm0(:, :), size( this%Qm0 ), MPI_DOUBLE_PRECISION, &
&       NRootProc, MPI_COMM_WORLD, ierror )
    if( (SimulationType .eq. MonteCarlo) .or. (SimulationType .eq. Gibbs) ) then
      call MPI_Bcast( this%DispTran, 1, MPI_DOUBLE_PRECISION, NRootProc, &
&       MPI_COMM_WORLD, ierror )
      call MPI_Bcast( this%NMoveAttempts, 1, MPI_INTEGER, NRootProc, &
&       MPI_COMM_WORLD, ierror )
      call MPI_Bcast( this%NMoveSuccesses, 1, MPI_INTEGER, NRootProc, &
&       MPI_COMM_WORLD, ierror )
      call MPI_Bcast( this%NMoveBiasedAttempts, 1, MPI_INTEGER, NRootProc, &
&       MPI_COMM_WORLD, ierror )
      call MPI_Bcast( this%NMoveBiasedSuccesses, 1, MPI_INTEGER, NRootProc, &
&       MPI_COMM_WORLD, ierror )
      if ( UseIntDegFreed ) then
        call MPI_Bcast( this%DispMolTran, 1, MPI_DOUBLE_PRECISION, NRootProc, &
&         MPI_COMM_WORLD, ierror )
        call MPI_Bcast( this%NMoveMolAttempts, 1, MPI_INTEGER, NRootProc, &
&         MPI_COMM_WORLD, ierror )
        call MPI_Bcast( this%NMoveMolSuccesses, 1, MPI_INTEGER, NRootProc, &
&         MPI_COMM_WORLD, ierror )
        call MPI_Bcast( this%NMoveBiasedMolAttempts, 1, MPI_INTEGER, NRootProc, &
&         MPI_COMM_WORLD, ierror )
        call MPI_Bcast( this%NMoveBiasedMolSuccesses, 1, MPI_INTEGER, NRootProc, &
&         MPI_COMM_WORLD, ierror )
      end if
      if( this%Molecule%isElongated ) then
        call MPI_Bcast( this%DispRot, 1, MPI_DOUBLE_PRECISION, NRootProc, &
&         MPI_COMM_WORLD, ierror )
        call MPI_Bcast( this%NRotateAttempts, 1, MPI_INTEGER, NRootProc, &
&         MPI_COMM_WORLD, ierror )
        call MPI_Bcast( this%NRotateSuccesses, 1, MPI_INTEGER, NRootProc, &
&         MPI_COMM_WORLD, ierror )
        call MPI_Bcast( this%NRotateBiasedAttempts, 1, MPI_INTEGER, NRootProc, &
&         MPI_COMM_WORLD, ierror )
        call MPI_Bcast( this%NRotateBiasedSuccesses, 1, MPI_INTEGER, &
&         NRootProc, MPI_COMM_WORLD, ierror )
        if ( UseIntDegFreed ) then
          call MPI_Bcast( this%DispRot, 1, MPI_DOUBLE_PRECISION, NRootProc, &
&           MPI_COMM_WORLD, ierror )
          call MPI_Bcast( this%NRotateAttempts, 1, MPI_INTEGER, NRootProc, &
&           MPI_COMM_WORLD, ierror )
          call MPI_Bcast( this%NRotateSuccesses, 1, MPI_INTEGER, NRootProc, &
&           MPI_COMM_WORLD, ierror )
          call MPI_Bcast( this%NRotateBiasedAttempts, 1, MPI_INTEGER, NRootProc, &
&           MPI_COMM_WORLD, ierror )
          call MPI_Bcast( this%NRotateBiasedSuccesses, 1, MPI_INTEGER, &
&           NRootProc, MPI_COMM_WORLD, ierror )
        end if
      end if
    end if
    if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
      call MPI_Bcast( this%WF, size( this%WF ), MPI_DOUBLE_PRECISION, &
&       NRootProc, MPI_COMM_WORLD, ierror )
      call MPI_BCast( this%NState, size( this%NState ), MPI_INTEGER, &
&       NRootProc, MPI_COMM_WORLD, ierror )
      call MPI_BCast( this%NStateWF, size( this%NStateWF ), MPI_INTEGER, &
&       NRootProc, MPI_COMM_WORLD, ierror )
    end if
#endif

    ! Update old positions
    this%Pm0old = this%Pm0

  end subroutine TComponent_RestartRead





end module ms2_component

