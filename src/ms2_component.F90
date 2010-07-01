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

#ifndef TRANS
#define TRANS 0
#endif

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2_component.F90...'
#endif

module ms2_component

!#ifdev MPI_VER > 0
!  use mpi
!#endif

  use ms2_accumulator
  use ms2_global
  use ms2_molecule
  use ms2_site



!==============================================================!
!  Type TComponent                                             !
!==============================================================!

  type TComponent

    ! Positions and orientations of test particles
    real(RK), pointer :: P0Test(:, :), Q0Test(:, :)

    ! Centers of mass positions and their derivatives
    real(RK), pointer :: P0(:, :)
    real(RK), pointer :: P0Save(:, :)
    real(RK), pointer :: P0old(:, :)
    real(RK), pointer :: P1(:, :)
    real(RK), pointer :: P2(:, :)
    real(RK), pointer :: P3(:, :)
    real(RK), pointer :: P4(:, :)
    real(RK), pointer :: P5(:, :)

    ! Quaternion parameters and their derivatives
    real(RK), pointer :: Q0(:, :)
    real(RK), pointer :: Q0Save(:, :)
    real(RK), pointer :: Q0tmp(:, :)
    real(RK), pointer :: Q1(:, :)
    real(RK), pointer :: Q2(:, :)
    real(RK), pointer :: Q3(:, :)
    real(RK), pointer :: Q4(:, :)

    ! Angular velocities and their derivatives
    real(RK), pointer :: W0(:, :)
    real(RK), pointer :: W1(:, :)
    real(RK), pointer :: W2(:, :)
    real(RK), pointer :: W3(:, :)
    real(RK), pointer :: W4(:, :)

    ! Displacement
    real(RK), pointer :: Disp(:, :)

    ! Total forces
    real(RK), pointer :: F(:, :)
#if MPI_VER > 0
    real(RK), pointer :: FAll(:, :)
#endif

    ! Total torques
    real(RK), pointer :: T(:, :)
#if MPI_VER > 0
    real(RK), pointer :: TAll(:, :)
#endif
#if  TRANS == 1
!TRANSPORT_start
 ! Transport
    real(RK), pointer :: KinETran(:,:)

    real(RK), pointer :: fs(:,:)
    real(RK), pointer :: fb(:,:)
    real(RK), pointer :: ftc(:,:)
    real(RK), pointer :: frc(:,:)

    real(RK), pointer :: ftc1(:,:)
    real(RK), pointer :: ftc2(:,:)
    real(RK), pointer :: ftc3(:,:)

    real(RK), pointer :: frc1(:,:)
    real(RK), pointer :: frc2(:,:)
    real(RK), pointer :: frc3(:,:)
#if MPI_VER > 0
    real(RK), pointer :: fsAll(:,:)
    real(RK), pointer :: fbAll(:,:)
    real(RK), pointer :: frcAll(:,:)

! Components of the FTC Tensor(3)
    real(RK), pointer :: ftc1All(:,:)
    real(RK), pointer :: ftc2All(:,:)
    real(RK), pointer :: ftc3All(:,:)

! Components of the FRC Tensor(3)
    real(RK), pointer :: frc1All(:,:)
    real(RK), pointer :: frc2All(:,:)
    real(RK), pointer :: frc3All(:,:)

#endif
!TRANSPORT_END
#endif

    ! Total dipole moment of molecule for reaction field
    real(RK), pointer :: MueX(:), MueY(:), MueZ(:)

    ! Torques from reaction field, space fixed
    real(RK), pointer :: tRFX(:), tRFY(:), tRFZ(:)

    ! Total dipole moment of test particles for reaction field
    real(RK), pointer :: MueXTest(:), MueYTest(:), MueZTest(:)

    ! Gear corrector local arrays
    real(RK), pointer :: Corr0(:, :)
    real(RK), pointer :: Corr1(:, :)

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
    integer  :: ChemPotMethod, WFMethod
    integer  :: FluctState
    real(RK) :: ChemPot, WidomContribution
!DEBUG
    real(RK) :: ChemPot1, ChemPot2
!DEBUG
    real(RK) :: ChemPot0, PartialMolarVolume
    real(RK) :: VarChemPot, VarPartialMolarVolume

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

    ! Allocate maximum allowed MC displacements
    if( SimulationType .eq. MonteCarlo .or. MCOverlapReduction ) then
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
      if( this%WFMethod .eq. WFMethodGuess .or. &
&         this%WFMethod .eq. WFMethodOptSet ) then
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
      call AllocationError( stat, 'fluctuating particle components', &
&       this%NFluctMax + 1 )
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
&       EnsembleType .eq. EnsembleTypeHA) then
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
&       EnsembleType .eq. EnsembleTypeHA ) then
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
    integer :: i
    integer :: stat
#if  TRANS == 1
    real(RK), pointer:: Q00(: , :)
#endif

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
    nullify( this%F )
#if MPI_VER > 0
    nullify( this%FAll )
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
    nullify( this%fs )
    nullify( this%fb )
    nullify( this%frc )
    nullify( this%ftc )

    nullify( this%ftc1)
    nullify( this%ftc2 )
    nullify( this%ftc3 )

    nullify( this%frc1)
    nullify( this%frc2 )
    nullify( this%frc3 )

#if MPI_VER > 0
    nullify( this%fsAll )
    nullify( this%fbAll )
    nullify( this%frcAll )

    nullify( this%ftc1All )
    nullify( this%ftc2All )
    nullify( this%ftc3All )

    nullify( this%frc1All )
    nullify( this%frc2All )
    nullify( this%frc3All )
#endif

    ! Transport

    allocate( this%KinETran( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%fs( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%fb( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%ftc( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%frc( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%ftc1( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%ftc2( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%ftc3( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%frc1( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%frc2( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%frc3( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )


    allocate( Q00( np, 4 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

#if MPI_VER > 0
    allocate( this%fsAll( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%fbAll( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%frcAll( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%ftc1All( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%ftc2All( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%ftc3All( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%frc1All( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%frc2All( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    allocate( this%frc3All( np, 3 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

#endif

    allocate( Q00( np, 4 ), STAT = stat )
    call AllocationError( stat, 'particles', np )

    this%fs(: , :)   = 0._RK
    this%fb(: , :)   = 0._RK
    this%ftc(: , :)  = 0._RK
    this%frc(: , :)  = 0._RK
    Q00(: , :)       = 0._RK

    this%ftc1(:,:)  = 0._RK
    this%ftc2(:,:)  = 0._RK
    this%ftc3(:,:)  = 0._RK

    this%frc1(:,:)  = 0._RK
    this%frc2(:,:)  = 0._RK
    this%frc3(:,:)  = 0._RK
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

      ! Total forces
      allocate( this%F( np, 3 ), STAT = stat )
      call AllocationError( stat, 'particles', np )
#if MPI_VER > 0
      allocate( this%FAll( np, 3 ), STAT = stat )
      call AllocationError( stat, 'particles', np )
#endif

    end if

    if( this%Molecule%isElongated ) then

      ! Quaternion parameters
      allocate( this%Q0( np, 4 ), STAT = stat )
      call AllocationError( stat, 'particles', np )
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
    if( SimulationType .eq. MolecularDynamics &
&     .and. IntegratorType .eq. IntegratorTypeGear ) then
      allocate( this%Corr0( np, merge( 4, 3, this%Molecule%isElongated ) ), &
&               STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%Corr1( np, merge( 4, 3, this%Molecule%isElongated ) ), &
&               STAT = stat )
      call AllocationError( stat, 'particles', np )
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
      this%Molecule%SiteLJ126(i)%PX => this%P0(:, 1)
      this%Molecule%SiteLJ126(i)%PY => this%P0(:, 2)
      this%Molecule%SiteLJ126(i)%PZ => this%P0(:, 3)
      if( ntest > 0 ) then
        this%Molecule%SiteLJ126(i)%PXTest => this%P0Test(:, 1)
        this%Molecule%SiteLJ126(i)%PYTest => this%P0Test(:, 2)
        this%Molecule%SiteLJ126(i)%PZTest => this%P0Test(:, 3)
        end if
#if TRANS==1
        if (this%Molecule%isElongated) then
          this%Molecule%SiteLJ126(i)%Q0r => this%Q0
        else
          this%Molecule%SiteLJ126(i)%Q0r => Q00
        end if
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
      if (this%Molecule%isElongated) then
        this%Molecule%SiteCharge(i)%Q0r => this%Q0
      else
        this%Molecule%SiteCharge(i)%Q0r => Q00
      end if
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
      if (this%Molecule%isElongated) then
        this%Molecule%SiteDipole(i)%Q0r => this%Q0
      else
        this%Molecule%SiteDipole(i)%Q0r => Q00
      end if
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
      if (this%Molecule%isElongated) then
        this%Molecule%SiteQuadrupole(i)%Q0r => this%Q0
      else
        this%Molecule%SiteQuadrupole(i)%Q0r => Q00
      end if
#endif
    end do

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
#if  TRANS == 1
! Transport !TRANSPORT_start

    if( associated( this%KinETran) ) then
      deallocate( this%KinETran )
    end if
    if( associated( this%fs ) ) then
      deallocate( this%fs )
    end if
    if( associated( this%fb ) ) then
      deallocate( this%fb )
    end if
    if( associated( this%ftc ) ) then
      deallocate( this%ftc )
    end if
    if( associated( this%frc ) ) then
      deallocate( this%frc )
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

     ! Fluctuating particle states
    if( associated( this%NFluctComp ) ) then
      deallocate( this%NFluctComp )
    end if
    if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
      if( associated( this%WF ) ) then
        deallocate( this%WF )
      end if
      if( associated( this%NState ) ) then
        deallocate( this%NState )
      end if
      if( associated( this%NStateWF ) ) then
        deallocate( this%NStateWF )
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
!DEBUG
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
!  Subroutine TComponent_RemoveNetMomentum                      !
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
      P(i) = P(i) &
&       + this%Molecule%Mass * sum( this%P1(1:this%NPart, i) )
      if( i <= this%Molecule%NDFRot ) &
&       L(i) = L(i) &
&         + this%Molecule%MOI(i) * sum( this%W0(1:this%NPart, i) )
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

  subroutine TComponent_Mol2Atom( this, np )

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
    real(RK)                       :: PX(np), PY(np), PZ(np)
    real(RK)                       :: q1, q2, q3, q4, qinv
    real(RK)                       :: A11(np), A12(np), A13(np)
    real(RK)                       :: A21(np), A22(np), A23(np)
    real(RK)                       :: A31(np), A32(np), A33(np)
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    real(RK)                       :: mue1, mue2, mue3
    type(TSiteLJ126), pointer      :: pLJ126
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j

    ! Broadcast positions and orientations to all processes
#if MPI_VER > 0
    ! in MC simulations, we only communicate during common equilibration
    if ( SimulationType .ne. MonteCarlo .or. ((Equilibration .and. CommonEqui) )) then
      call MPI_Bcast( this%P0(:, :), size( this%P0 ), &
&      MPI_RK, NRootProc, Communicator, ierror )
      if( this%Molecule%isElongated ) &
&      call MPI_Bcast( this%Q0(:, :), size( this%Q0 ), &
&       MPI_RK, NRootProc, Communicator, ierror )
    endif

#endif

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    ! Check number of rotation axes
    if( this%Molecule%isElongated ) then

      ! Loop over molecules
      do i = 1, np
        ! Positions and quaternions of particle i
        PX(i) = this%P0(i, 1)
        PY(i) = this%P0(i, 2)
        PZ(i) = this%P0(i, 3)
        q1 = this%Q0(i, 1)
        q2 = this%Q0(i, 2)
        q3 = this%Q0(i, 3)
        q4 = this%Q0(i, 4)

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
        this%Q0(i, 1) = q1
        this%Q0(i, 2) = q2
        this%Q0(i, 3) = q3
        this%Q0(i, 4) = q4

        ! Calculate rotation matrix elements   ! FIXME Hier ist was doppelt, diese Berechnung wird nochmals gemacht
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
          pLJ126%RX(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pLJ126%RY(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pLJ126%RZ(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
        end do
      end do

      ! Loop over charge sites in molecule
      do j = 1, this%Molecule%NCharge
        pCharge => this%Molecule%SiteCharge(j)
        r1 = pCharge%r(1) * BoxLengthInv
        r2 = pCharge%r(2) * BoxLengthInv
        r3 = pCharge%r(3) * BoxLengthInv
        do i = 1, np
          pCharge%RX(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pCharge%RY(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pCharge%RZ(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
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
          pDipole%RX(i) = PX(i) + r1 * A11(i) + r2 * A21(i) + r3 * A31(i)
          pDipole%RY(i) = PY(i) + r1 * A12(i) + r2 * A22(i) + r3 * A32(i)
          pDipole%RZ(i) = PZ(i) + r1 * A13(i) + r2 * A23(i) + r3 * A33(i)
          pDipole%OX(i) = or1 * A11(i) + or2 * A21(i) + or3 * A31(i)
          pDipole%OY(i) = or1 * A12(i) + or2 * A22(i) + or3 * A32(i)
          pDipole%OZ(i) = or1 * A13(i) + or2 * A23(i) + or3 * A33(i)
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
        mue1 = this%Molecule%Mue(1)
        mue2 = this%Molecule%Mue(2)
        mue3 = this%Molecule%Mue(3)
        do i = 1, np
          this%MueX(i) = mue1 * A11(i) + mue2 * A21(i) + mue3 * A31(i)
          this%MueY(i) = mue1 * A12(i) + mue2 * A22(i) + mue3 * A32(i)
          this%MueZ(i) = mue1 * A13(i) + mue2 * A23(i) + mue3 * A33(i)
        end do
      end if

    else

      ! Loop over LJ126 sites in molecule
      do i = 1, this%Molecule%NLJ126
        pLJ126 => this%Molecule%SiteLJ126(i)
        do j = 1, np
          pLJ126%RX(j) = this%P0(j, 1)
          pLJ126%RY(j) = this%P0(j, 2)
          pLJ126%RZ(j) = this%P0(j, 3)
        end do
      end do

    end if

  end subroutine TComponent_Mol2Atom



!==============================================================!
!  Subroutine TComponent_Mol2Atom1                             !
!==============================================================!

  subroutine TComponent_Mol2Atom1( this, n )

    implicit none

    ! Declare arguments
    type(TComponent)    :: this
    integer, intent(in) :: n

    ! Declare local variables
    real(RK)                       :: BoxLengthInv
    real(RK)                       :: PXi, PYi, PZi
    real(RK)                       :: q1, q2, q3, q4, qinv
    real(RK)                       :: A11, A12, A13, A21, A22, A23, &
&                                     A31, A32, A33
    real(RK)                       :: r1, r2, r3, or1, or2, or3
    real(RK)                       :: mue1, mue2, mue3
    type(TSiteLJ126), pointer      :: pLJ126
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i

    ! Assign local variables
    BoxLengthInv = 1._RK / this%BoxLength

    ! Positions of particle n
    PXi = this%P0(n, 1)
    PYi = this%P0(n, 2)
    PZi = this%P0(n, 3)

    ! Check number of rotation axes
    if( this%Molecule%isElongated ) then

      ! Normalise quaternions
      q1 = this%Q0(n, 1)
      q2 = this%Q0(n, 2)
      q3 = this%Q0(n, 3)
      q4 = this%Q0(n, 4)
#if ARCH == 3
      qinv = rsqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#else
      qinv = 1._RK / sqrt( q1**2 + q2**2 + q3**2 + q4**2 )
#endif
      q1 = q1 * qinv
      q2 = q2 * qinv
      q3 = q3 * qinv
      q4 = q4 * qinv
      this%Q0(n, 1) = q1
      this%Q0(n, 2) = q2
      this%Q0(n, 3) = q3
      this%Q0(n, 4) = q4

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

      ! Loop over LJ126 sites in molecule
      do i = 1, this%Molecule%NLJ126
        pLJ126 => this%Molecule%SiteLJ126(i)
        r1 = pLJ126%r(1) * BoxLengthInv
        r2 = pLJ126%r(2) * BoxLengthInv
        r3 = pLJ126%r(3) * BoxLengthInv
        pLJ126%RX(n) = PXi + r1 * A11 + r2 * A21 + r3 * A31
        pLJ126%RY(n) = PYi + r1 * A12 + r2 * A22 + r3 * A32
        pLJ126%RZ(n) = PZi + r1 * A13 + r2 * A23 + r3 * A33
      end do

      ! Loop over charge sites in molecule
      do i = 1, this%Molecule%NCharge
        pCharge => this%Molecule%SiteCharge(i)
        r1 = pCharge%r(1) * BoxLengthInv
        r2 = pCharge%r(2) * BoxLengthInv
        r3 = pCharge%r(3) * BoxLengthInv
        pCharge%RX(n) = PXi + r1 * A11 + r2 * A21 + r3 * A31
        pCharge%RY(n) = PYi + r1 * A12 + r2 * A22 + r3 * A32
        pCharge%RZ(n) = PZi + r1 * A13 + r2 * A23 + r3 * A33
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
        pDipole%RX(n) = PXi + r1 * A11 + r2 * A21 + r3 * A31
        pDipole%RY(n) = PYi + r1 * A12 + r2 * A22 + r3 * A32
        pDipole%RZ(n) = PZi + r1 * A13 + r2 * A23 + r3 * A33
        pDipole%OX(n) = or1 * A11 + or2 * A21 + or3 * A31
        pDipole%OY(n) = or1 * A12 + or2 * A22 + or3 * A32
        pDipole%OZ(n) = or1 * A13 + or2 * A23 + or3 * A33
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
        pQuadrupole%RX(n) = PXi + r1 * A11 + r2 * A21 + r3 * A31
        pQuadrupole%RY(n) = PYi + r1 * A12 + r2 * A22 + r3 * A32
        pQuadrupole%RZ(n) = PZi + r1 * A13 + r2 * A23 + r3 * A33
        pQuadrupole%OX(n) = or1 * A11 + or2 * A21 + or3 * A31
        pQuadrupole%OY(n) = or1 * A12 + or2 * A22 + or3 * A32
        pQuadrupole%OZ(n) = or1 * A13 + or2 * A23 + or3 * A33
      end do

      ! Rotate total dipole moment
      if( CutoffMode .eq. CenterofMass ) then
        mue1 = this%Molecule%Mue(1)
        mue2 = this%Molecule%Mue(2)
        mue3 = this%Molecule%Mue(3)
        this%MueX(n) = mue1 * A11 + mue2 * A21 + mue3 * A31
        this%MueY(n) = mue1 * A12 + mue2 * A22 + mue3 * A32
        this%MueZ(n) = mue1 * A13 + mue2 * A23 + mue3 * A33
      end if

    else

      ! Loop over LJ126 sites in molecule
      do i = 1, this%Molecule%NLJ126
        pLJ126 => this%Molecule%SiteLJ126(i)
        pLJ126%RX(n) = PXi
        pLJ126%RY(n) = PYi
        pLJ126%RZ(n) = PZi
      end do

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
    type(TSiteLJ126), pointer      :: pLJ126
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
          this%MueXTest(i) = mue1 * A11(i) + mue2 * A21(i) + mue3 * A31(i)
          this%MueYTest(i) = mue1 * A12(i) + mue2 * A22(i) + mue3 * A32(i)
          this%MueZTest(i) = mue1 * A13(i) + mue2 * A23(i) + mue3 * A33(i)
        end do
      end if

    else

      ! Loop over LJ126 sites in molecule
      do i = 1, this%Molecule%NLJ126
        pLJ126 => this%Molecule%SiteLJ126(i)
        do j = 1, np
          pLJ126%RXTest(j) = this%P0Test(j, 1)
          pLJ126%RYTest(j) = this%P0Test(j, 2)
          pLJ126%RZTest(j) = this%P0Test(j, 3)
        end do
      end do

    end if

  end subroutine TComponent_Mol2AtomTest



!==============================================================!
!  Subroutine TComponent_Atom2Mol                              !
!==============================================================!

  subroutine TComponent_Atom2Mol( this, np )

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
  
    real(RK)                       :: A11, A12, A13, A21, A22, A23, &
&                                     A31, A32, A33
    type(TSiteLJ126), pointer      :: pLJ126
    type(TSiteCharge), pointer     :: pCharge
    type(TSiteDipole), pointer     :: pDipole
    type(TSiteQuadrupole), pointer :: pQuadrupole
    integer                        :: i, j

    ! Assign local variables
    BoxLength = this%BoxLength
#if  TRANS == 1
    BoxLength_dt = this%BoxLength/TimeStep !TRANSPORT_thisline
#endif
    ! Initialize forces
    this%F(1:np, :) = 0._RK
#if  TRANS == 1
    !TRANSPORT_start
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
!TRANSPORT_END
#endif

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
      ! Loop over LJ126 sites in molecule
      do j = 1, this%Molecule%NLJ126
        pLJ126 => this%Molecule%SiteLJ126(j)
        do i = 1, np
          fx = pLJ126%FX(i)
          fy = pLJ126%FY(i)
          fz = pLJ126%FZ(i)
#if  TRANS == 1
          !TRANSPORT_start
          vsx = pLJ126%vsLJx(i)
          vsy = pLJ126%vsLJy(i)
          vsz = pLJ126%vsLJz(i)
          vsux= pLJ126%vsuLJx(i)
          vsuy= pLJ126%vsuLJy(i)
          vsuz= pLJ126%vsuLJz(i)
          vbx = pLJ126%vbLJx(i)
          vby = pLJ126%vbLJy(i)
          vbz = pLJ126%vbLJz(i)
          cx  = pLJ126%cLJx(i)
          cy  = pLJ126%cLJy(i)
          cz  = pLJ126%cLJz(i)
          tux = pLJ126%tuLJx(i)
          tuy = pLJ126%tuLJy(i)
          tuz = pLJ126%tuLJz(i)
          tlx = pLJ126%tlLJx(i)
          tly = pLJ126%tlLJy(i)
          tlz = pLJ126%tlLJz(i)
          tdx = pLJ126%tdLJx(i)
          tdy = pLJ126%tdLJy(i)
          tdz = pLJ126%tdLJz(i)
          !TRANSPORT_END
#endif
          r1x = ( pLJ126%RX(i) - RX(i) ) * BoxLength
          r1y = ( pLJ126%RY(i) - RY(i) ) * BoxLength
          r1z = ( pLJ126%RZ(i) - RZ(i) ) * BoxLength
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
          vsux= pCharge%vsuCx(i)
          vsuy= pCharge%vsuCy(i)
          vsuz= pCharge%vsuCz(i)
          vbx = pCharge%vbCx(i)
          vby = pCharge%vbCy(i)
          vbz = pCharge%vbCz(i)
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
          vsux= pDipole%vsuDx(i)
          vsuy= pDipole%vsuDy(i)
          vsuz= pDipole%vsuDz(i)
          vbx = pDipole%vbDx(i)
          vby = pDipole%vbDy(i)
          vbz = pDipole%vbDz(i)
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
          vsux= pQuadrupole%vsuQx(i)
          vsuy= pQuadrupole%vsuQy(i)
          vsuz= pQuadrupole%vsuQz(i)
          vbx = pQuadrupole%vbQx(i)
          vby = pQuadrupole%vbQy(i)
          vbz = pQuadrupole%vbQz(i)
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

      ! Loop over LJ126 sites in molecule
      do j = 1, this%Molecule%NLJ126
        pLJ126 => this%Molecule%SiteLJ126(j)
        do i = 1, np
#if  TRANS == 1
        !TRANSPORT_start
          vsx = pLJ126%vsLJx(i)
          vsy = pLJ126%vsLJy(i)
          vsz = pLJ126%vsLJz(i)
          vsux= pLJ126%vsuLJx(i)
          vsuy= pLJ126%vsuLJy(i)
          vsuz= pLJ126%vsuLJz(i)
          vbx = pLJ126%vbLJx(i)
          vby = pLJ126%vbLJy(i)
          vbz = pLJ126%vbLJz(i)
          cx  = pLJ126%cLJx(i)
          cy  = pLJ126%cLJy(i)
          cz  = pLJ126%cLJz(i)
          !TRANSPORT_END
#endif
          this%F(i, 1) = this%F(i, 1) + pLJ126%FX(i)
          this%F(i, 2) = this%F(i, 2) + pLJ126%FY(i)
          this%F(i, 3) = this%F(i, 3) + pLJ126%FZ(i)
#if  TRANS == 1
          !TRANSPORT_start
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
!TRANSPORT_END
#endif
        end do
      end do

    end if

    ! Reduce forces and torques from all processes
#if MPI_VER > 0
    call MPI_Reduce( this%F(:, :), this%FAll(:, :), size( this%F ), &
&     MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    if( this%Molecule%isElongated ) then
     call MPI_Reduce( this%T(:, :), this%TAll(:, :), size( this%T ), &
&      MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    end if
#if  TRANS == 1
! Transport  !TRANSPORT_start
    call MPI_Reduce( this%FB(:, :), this%FBAll(:, :), size( this%FB ), &
&     MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    call MPI_Reduce( this%FS(:, :), this%FSAll(:, :), size( this%FS ), &
&       MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )

    call MPI_Reduce( this%FTC1(:, :), this%FTC1All(:, :), size( this%FTC1 ), &
&     MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    call MPI_Reduce( this%FTC2(:, :), this%FTC2All(:, :), size( this%FTC2 ), &
&     MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    call MPI_Reduce( this%FTC3(:, :), this%FTC3All(:, :), size( this%FTC3 ), &
&     MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )

    call MPI_Reduce( this%FRC1(:, :), this%FRC1All(:, :), size( this%FRC1 ), &
&       MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    call MPI_Reduce( this%FRC2(:, :), this%FRC2All(:, :), size( this%FRC2 ), &
&       MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    call MPI_Reduce( this%FRC3(:, :), this%FRC3All(:, :), size( this%FRC3 ), &
&       MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
!TRANSPORT_END
#endif

#endif

  end subroutine TComponent_Atom2Mol



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
        this%P0(i, j) = this%P0(i, j) &
&                     + this%P1(i, j) &
&                     + this%P2(i, j) &
&                     + this%P3(i, j) &
&                     + this%P4(i, j) &
&                     + this%P5(i, j)
        this%P1(i, j) = this%P1(i, j) &
&             + 2._RK * this%P2(i, j) &
&             + 3._RK * this%P3(i, j) &
&             + 4._RK * this%P4(i, j) &
&             + 5._RK * this%P5(i, j)
        this%P2(i, j) = this%P2(i, j) &
&             + 3._RK * this%P3(i, j) &
&             + 6._RK * this%P4(i, j) &
&             +10._RK * this%P5(i, j)
        this%P3(i, j) = this%P3(i, j) &
&             + 4._RK * this%P4(i, j) &
&             +10._RK * this%P5(i, j)
        this%P4(i, j) = this%P4(i, j) &
&             + 5._RK * this%P5(i, j)
      end do
    end do

    if( this%Molecule%IsElongated ) then

      ! Predict quaternion parameters and their derivatives
      do j = 1, 4
        do i = 1, np
          this%Q0(i, j) = this%Q0(i, j) &
&                       + this%Q1(i, j) &
&                       + this%Q2(i, j) &
&                       + this%Q3(i, j) &
&                       + this%Q4(i, j)
          this%Q1(i, j) = this%Q1(i, j) &
&               + 2._RK * this%Q2(i, j) &
&               + 3._RK * this%Q3(i, j) &
&               + 4._RK * this%Q4(i, j)
          this%Q2(i, j) = this%Q2(i, j) &
&               + 3._RK * this%Q3(i, j) &
&               + 6._RK * this%Q4(i, j)
          this%Q3(i, j) = this%Q3(i, j) &
&               + 4._RK * this%Q4(i, j)
        end do
      end do

      ! Predict angular velocities and their derivatives
      do j = 1, nra
        do i = 1, np
          this%W0(i, j) = this%W0(i, j) &
&                       + this%W1(i, j) &
&                       + this%W2(i, j) &
&                       + this%W3(i, j) &
&                       + this%W4(i, j)
          this%W1(i, j) = this%W1(i, j) &
&               + 2._RK * this%W2(i, j) &
&               + 3._RK * this%W3(i, j) &
&               + 4._RK * this%W4(i, j)
          this%W2(i, j) = this%W2(i, j) &
&               + 3._RK * this%W3(i, j) &
&               + 6._RK * this%W4(i, j)
          this%W3(i, j) = this%W3(i, j) &
&               + 4._RK * this%W4(i, j)
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
    real(RK), pointer :: pF(:, :), pT(:, :)
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
        this%Corr0(i, j) = pF(i, j) &
&         * TimeStepSquared2 * BoxLengthInv * MassInv
        if( ConstantPressure .and. .not. NVTEquilibration ) &
&         this%Corr0(i, j) = this%Corr0(i, j) &
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
        this%Corr0(i, 1) = TimeStep2 * ( - this%Q0(i, 2) * this%W0(i, 1) &
&                                        - this%Q0(i, 3) * this%W0(i, 2) &
&                                        - this%Q0(i, 4) * this%W0(i, 3))
        this%Corr0(i, 2) = TimeStep2 * ( + this%Q0(i, 1) * this%W0(i, 1) &
&                                        - this%Q0(i, 4) * this%W0(i, 2) &
&                                        + this%Q0(i, 3) * this%W0(i, 3))
        this%Corr0(i, 3) = TimeStep2 * ( + this%Q0(i, 4) * this%W0(i, 1) &
&                                        + this%Q0(i, 1) * this%W0(i, 2) &
&                                        - this%Q0(i, 2) * this%W0(i, 3))
        this%Corr0(i, 4) = TimeStep2 * ( - this%Q0(i, 3) * this%W0(i, 1) &
&                                        + this%Q0(i, 2) * this%W0(i, 2) &
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
          this%Corr0(i, 1) = (pT(i, 1) + this%W0(i, 2) * this%W0(i, 3) * &
&                             Moi23) * TMoi1
          this%Corr0(i, 2) = (pT(i, 2) + this%W0(i, 3) * this%W0(i, 1) * &
&                             Moi31) * TMoi2
          this%Corr0(i, 3) = (pT(i, 3) + this%W0(i, 1) * this%W0(i, 2) * &
&                             Moi12) * TMoi3
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
        this%Q1(i, 1) = TimeStep2 * ( - this%Q0tmp(i, 2) * this%W0(i, 1) &
&                                     - this%Q0tmp(i, 3) * this%W0(i, 2) &
&                                     - this%Q0tmp(i, 4) * this%W0(i, 3))
        this%Q1(i, 2) = TimeStep2 * ( + this%Q0tmp(i, 1) * this%W0(i, 1) &
&                                     - this%Q0tmp(i, 4) * this%W0(i, 2) &
&                                     + this%Q0tmp(i, 3) * this%W0(i, 3))
        this%Q1(i, 3) = TimeStep2 * ( + this%Q0tmp(i, 4) * this%W0(i, 1) &
&                                     + this%Q0tmp(i, 1) * this%W0(i, 2) &
&                                     - this%Q0tmp(i, 2) * this%W0(i, 3))
        this%Q1(i, 4) = TimeStep2 * ( - this%Q0tmp(i, 3) * this%W0(i, 1) &
&                                     + this%Q0tmp(i, 2) * this%W0(i, 2) &
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
    real(RK), pointer :: pF(:, :), pT(:, :)
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
          this%W1(i, 1) = (pT(i, 1) + this%W0(i, 2) * this%W0(i, 3) * Moi23) * &
&                         TMoi1
          this%W1(i, 2) = (pT(i, 2) + this%W0(i, 3) * this%W0(i, 1) * Moi31) * &
&                         TMoi2
          this%W1(i, 3) = (pT(i, 3) + this%W0(i, 1) * this%W0(i, 2) * Moi12) * &
&                         TMoi3
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
        this%Q1(i, 1) = TimeStep2 * ( - this%Q0(i, 2) * this%W0(i, 1) &
&                                     - this%Q0(i, 3) * this%W0(i, 2) &
&                                     - this%Q0(i, 4) * this%W0(i, 3))
        this%Q1(i, 2) = TimeStep2 * ( + this%Q0(i, 1) * this%W0(i, 1) &
&                                     - this%Q0(i, 4) * this%W0(i, 2) &
&                                     + this%Q0(i, 3) * this%W0(i, 3))
        this%Q1(i, 3) = TimeStep2 * ( + this%Q0(i, 4) * this%W0(i, 1) &
&                                     + this%Q0(i, 1) * this%W0(i, 2) &
&                                     - this%Q0(i, 2) * this%W0(i, 3))
        this%Q1(i, 4) = TimeStep2 * ( - this%Q0(i, 3) * this%W0(i, 1) &
&                                     + this%Q0(i, 2) * this%W0(i, 2) &
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
    if( this%NPart >= this%NPartMax ) then
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
    integer :: np, i

    ! Assign local variables
    np = this%NPart

    ! Check for root process
    if( .not. RootProc ) return

    ! Save contents to restart file
    write( iounit_restart, '(I10)' ) np

    ! Centers of mass positions
    do i = 1, np
      write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P0( i, : )
    end do

    if( SimulationType .eq. MolecularDynamics ) then
      ! Centers of mass positions' derivatives
      do i = 1, np
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P1( i, : )
      end do
      do i = 1, np
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P2( i, : )
      end do

      if( IntegratorType .eq. IntegratorTypeGear ) then
        do i = 1, np
          write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P3( i, : )
        end do
        do i = 1, np
          write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P4( i, : )
        end do
        do i = 1, np
          write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P5( i, : )
        end do
      end if
    else
      write( iounit_restart, '(ES20.12E3)' ) this%DispTran
      write( iounit_restart, '(2I10)' ) this%NMoveAttempts, this%NMoveSuccesses
      write( iounit_restart, '(2I10)' ) this%NMoveBiasedAttempts, &
&       this%NMoveBiasedSuccesses
    end if

    if( this%Molecule%isElongated ) then
      ! Quaternion parameters
      do i = 1, np
        write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q0( i, : )
      end do

      if( SimulationType .eq. MolecularDynamics ) then
        ! Quaternion parameters' derivatives
        do i = 1, np
          write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q1( i, : )
        end do

        if( IntegratorType .eq. IntegratorTypeGear ) then
          do i = 1, np
            write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q2( i, : )
          end do
          do i = 1, np
            write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q3( i, : )
          end do
          do i = 1, np
            write( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q4( i, : )
          end do
        end if

        ! Angular velocities and their derivatives
        do i = 1, np
          write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W0( i, : )
        end do
        do i = 1, np
          write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W1( i, : )
        end do

        if( IntegratorType .eq. IntegratorTypeGear ) then
          do i = 1, np
            write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W2( i, : )
          end do
          do i = 1, np
            write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W3( i, : )
          end do
          do i = 1, np
            write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W4( i, : )
          end do
        end if
      else
        write( iounit_restart, '(ES20.12E3)' ) this%DispRot
        write( iounit_restart, '(2I10)' ) this%NRotateAttempts, &
&         this%NRotateSuccesses
        write( iounit_restart, '(2I10)' ) this%NRotateBiasedAttempts, &
&         this%NRotateBiasedSuccesses
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

    if( RootProc ) then

      ! Read contents from restart file
      read( iounit_restart, '(I10)' ) np
      if( np > this%NPartMax ) &
&       call Error( 'Not enough memory to read particles from restart file' )
      this%NPart = np

      ! Centers of mass positions
      do i = 1, np
        read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P0( i, : )
      end do

      if( SimulationType .eq. MolecularDynamics ) then
        ! Centers of mass positions' derivatives
        do i = 1, np
          read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P1( i, : )
        end do
        do i = 1, np
          read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P2( i, : )
        end do

        if( IntegratorType .eq. IntegratorTypeGear ) then
          do i = 1, np
            read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P3( i, : )
          end do
          do i = 1, np
            read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P4( i, : )
          end do
          do i = 1, np
            read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%P5( i, : )
          end do
        end if
      else
        read( iounit_restart, '(ES20.12E3)' ) this%DispTran
        read( iounit_restart, '(2I10)' ) this%NMoveAttempts, this%NMoveSuccesses
        read( iounit_restart, '(2I10)' ) this%NMoveBiasedAttempts, &
&         this%NMoveBiasedSuccesses
      end if

      if( this%Molecule%isElongated ) then
        ! Quaternion parameters
        do i = 1, np
          read( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q0( i, : )
        end do

        if( SimulationType .eq. MolecularDynamics ) then
          ! Quaternion parameters' derivatives
          do i = 1, np
            read( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q1( i, : )
          end do

          if( IntegratorType .eq. IntegratorTypeGear ) then
            do i = 1, np
              read( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q2( i, : )
            end do
            do i = 1, np
              read( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q3( i, : )
            end do
            do i = 1, np
              read( iounit_restart, '(4(ES20.12E3, :, ";"))' ) this%Q4( i, : )
            end do
          end if

          ! Angular velocities and their derivatives
          do i = 1, np
            read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W0( i, : )
          end do
          do i = 1, np
            read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W1( i, : )
          end do

          if( IntegratorType .eq. IntegratorTypeGear ) then
            do i = 1, np
              read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W2( i, : )
            end do
            do i = 1, np
              read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W3( i, : )
            end do
            do i = 1, np
              read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%W4( i, : )
            end do
          end if
        else
          read( iounit_restart, '(ES20.12E3)' ) this%DispRot
          read( iounit_restart, '(2I10)' ) this%NRotateAttempts, &
&           this%NRotateSuccesses
          read( iounit_restart, '(2I10)' ) this%NRotateBiasedAttempts, &
&           this%NRotateBiasedSuccesses
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
&     Communicator, ierror )
    call MPI_Bcast( this%P0(:, :), size( this%P0 ), MPI_RK, &
&     NRootProc, Communicator, ierror )
    if( this%Molecule%isElongated ) &
&     call MPI_Bcast( this%Q0(:, :), size( this%Q0 ), MPI_RK, &
&       NRootProc, Communicator, ierror )
    if( SimulationType .eq. MonteCarlo ) then
      call MPI_Bcast( this%DispTran, 1, MPI_RK, NRootProc, &
&       Communicator, ierror )
      call MPI_Bcast( this%NMoveAttempts, 1, MPI_INTEGER, NRootProc, &
&       Communicator, ierror )
      call MPI_Bcast( this%NMoveSuccesses, 1, MPI_INTEGER, NRootProc, &
&       Communicator, ierror )
      call MPI_Bcast( this%NMoveBiasedAttempts, 1, MPI_INTEGER, NRootProc, &
&       Communicator, ierror )
      call MPI_Bcast( this%NMoveBiasedSuccesses, 1, MPI_INTEGER, NRootProc, &
&       Communicator, ierror )
      if( this%Molecule%isElongated ) then
        call MPI_Bcast( this%DispRot, 1, MPI_RK, NRootProc, &
&         Communicator, ierror )
        call MPI_Bcast( this%NRotateAttempts, 1, MPI_INTEGER, NRootProc, &
&         Communicator, ierror )
        call MPI_Bcast( this%NRotateSuccesses, 1, MPI_INTEGER, NRootProc, &
&         Communicator, ierror )
        call MPI_Bcast( this%NRotateBiasedAttempts, 1, MPI_INTEGER, NRootProc, &
&         Communicator, ierror )
        call MPI_Bcast( this%NRotateBiasedSuccesses, 1, MPI_INTEGER, &
&         NRootProc, Communicator, ierror )
      end if
    end if
    if( this%ChemPotMethod .eq. ChemPotMethodGradIns ) then
      call MPI_Bcast( this%WF, size( this%WF ), MPI_RK, &
&       NRootProc, Communicator, ierror )
      call MPI_BCast( this%NState, size( this%NState ), MPI_INTEGER, &
&       NRootProc, Communicator, ierror )
      call MPI_BCast( this%NStateWF, size( this%NStateWF ), MPI_INTEGER, &
&       NRootProc, Communicator, ierror )
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
    real(RK), pointer :: pftc1(:,:), pftc2(:,:), pftc3(:,:)
    real(RK), pointer :: pfrc1(:,:), pfrc2(:,:), pfrc3(:,:)
    real(RK)          :: BoxLength_dt
    real(RK)          :: BoxLength_dt2

    !declare local variables
    BoxLength_dt = this%BoxLength/TimeStep
    BoxLength_dt2 = BoxLength_dt*BoxLength_dt


    this%FTC(:,:) = 0._RK
    this%FRC(:,:) = 0._RK
    this%KinETran(:,:) = 0._RK

 if (RootProc) then

#if MPI_VER > 0
    pftc1 => this%FTC1All(:,:)
    pftc2 => this%FTC2All(:,:)
    pftc3 => this%FTC3All(:,:)

    pfrc1 => this%FRC1All(:,:)
    pfrc2 => this%FRC2All(:,:)
    pfrc3 => this%FRC3All(:,:)
#else
    pftc1 => this%FTC1(:,:)
    pftc2 => this%FTC2(:,:)
    pftc3 => this%FTC3(:,:)

    pfrc1 => this%FRC1(:,:)
    pfrc2 => this%FRC2(:,:)
    pfrc3 => this%FRC3(:,:)
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
     do i = 1, this%NPart

       this%KinETran(i,j) = this%P1(i,j)*this%P1(i,j)

     end do
   end do

      this%KinETran(:,:) = this%KinETran(:,:)* this%Molecule%Mass*BoxLength_dt2

   end if

#endif
  end subroutine TComponent_ForceTransport

!TRANSPORT_END

end module ms2_component

