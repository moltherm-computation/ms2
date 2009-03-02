!==============================================================!
!  MOLECULAR SIMULATION PROGRAM MS2 Version 1.1 v12            !
!  (c) 2001 by Sergey Lishchuk, ITT                            !
!  (c) 2007 by Bernhard Eckl, ITT                              !
!==============================================================!
!  Module ms2_ensemble                                         !
!  Contains TEnsemble object                                   !
!==============================================================!

#ifndef ARCH
#define ARCH    0
#define FORTRAN 90
#define MPI_VER 0
#endif

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2_ensemble.F90...'
#endif

module ms2_ensemble

  use ms2_accumulator
  use ms2_component
  use ms2_global
  use ms2_potential
  use ms2_interaction
  use ms2_site



!==============================================================!
!  Type TEnsemble                                              !
!==============================================================!

  type TEnsemble

    ! Number of ensemble
    integer :: EnsembleNumber

    ! I/O unit for result file
    integer :: iounit_result

    ! I/O unit for running average result file
    integer :: iounit_runave

    ! I/O unit for final result file
    integer :: iounit_errors

    ! I/O unit for visualization file
    integer :: iounit_visual

    ! Maximum number of particles
    integer, pointer :: NPartMax

    ! Number of particles in ensemble
    integer :: NPart, NPartInitial
    integer :: NPartLBound, NPartUBound

    ! Maximum number of test particles
    integer :: NTestMax

    ! Maximum number of fluctuating states
    integer :: NFluctMax

    ! Number of degrees of freedom
    integer :: NDFTran, NDFRot, NDF

    ! Mass of piston
    real(RK) :: PistonMass

    ! Positions and orientations of test particles
    real(RK), pointer :: P0Test(:, :), Q0Test(:, :)

    ! Number of unit cells in one dimension of lattice
    integer :: NCells

    ! Number of components in ensemble
    integer :: NComponents, NRealComponents

    ! Maximum numbers of sites in components
    integer :: NLJ126Max, NChargeMax, NDipoleMax, NQuadrupoleMax

    ! Components
    type(TComponent), pointer :: Component(:)

    ! Interactions
    type(TInteraction), pointer :: Interaction(:, :)

    ! Initial values of temperature, pressure, density
    real(RK) :: RefTemperature, RefPressure, RefDensity

    ! Values of density, enthalpy, betaT, dHdp and their uncertainties
    ! in corresponding liquid simulation (for GE ensemble only)
    real(RK) :: LiqDensity, VarLiqDensity, LiqEnthalpy, VarLiqEnthalpy
    real(RK) :: LiqBetaT, VarLiqBetaT, LiqdHdP, VarLiqdHdP

    ! Current values of temperature, pressure, density
    real(RK) :: Temperature, Pressure, Density

    ! Velocity scaling factor for temperature control
    real(RK) :: scale

    ! Virial
    real(RK) :: Virial

    ! Scale coefficients for LJ126 epsilon and sigma
    real(RK), pointer :: ScaleEpsilon(:, :), ScaleSigma(:, :)

    ! Cutoff radii
    real(RK) :: RCutoffLJ126LJ126
    real(RK) :: RCutoffDipoleDipole
    real(RK) :: RCutoffDipoleQuadrupole
    real(RK) :: RCutoffQuadrupoleQuadrupole

    ! Characteristic dielectric constant for reaction field method
    real(RK) :: RFEpsilon

    ! Maximum cutoff radius
    real(RK) :: RCutoffMax2
    integer  :: NRCutoffMax

    ! Volume of simulation box and its derivatives
    real(RK) :: Volume0, Volume1, Volume2, Volume3, Volume4, Volume5

    ! Length of simulation box
    real(RK), pointer :: BoxLength

    ! Maximum allowed MC displacements
    real(RK) :: DispVol

    ! Number of MC attempts and successes
    integer :: NResizeAttempts, NResizeSuccesses
    integer :: NInsertAttempts, NInsertSuccesses
    integer :: NDeleteAttempts, NDeleteSuccesses

    ! Kinetic energy
    real(RK) :: EKin, EKinTran, EKinRot

    ! Potential energy
    real(RK) :: EPot

    ! Potential energy of test particles
    real(RK), pointer :: EPotTest(:)

    ! Long-range corrections
    real(RK) :: EPotCorrLJ, VirialCorrLJ
    real(RK) :: EPotCorrRF

    ! Accumulated sums, averages and errors
    ! 1.) Basic sums
    type(TAccumulator) :: SumPressure
    type(TAccumulator) :: SumDensity
    type(TAccumulator) :: SumTemperature
    type(TAccumulator) :: SumEPot
    type(TAccumulator) :: SumEnthalpy
    type(TAccumulator) :: SumVolume
    type(TAccumulator) :: SumVirial
    type(TAccumulator) :: SumNPart

    ! 2.) Combined sums
    type(TAccumulator) :: SumEPotSquared
    type(TAccumulator) :: SumEPotV
    type(TAccumulator) :: SumEPotVirial
    type(TAccumulator) :: SumEnthalpySquared
    type(TAccumulator) :: SumEnthalpyV
    type(TAccumulator) :: SumVolumeSquared

    ! 3.) Derived sums
    type(TAccumulator) :: SumBetaT
    type(TAccumulator) :: SumdHdP
    type(TAccumulator) :: SumdUdV
    type(TAccumulator) :: SumCV
    type(TAccumulator) :: SumCP
    type(TAccumulator) :: SumAlphaP


  ! Ewald Parameters
  real(RK),pointer :: Ewald_Prefac(:)
  real(RK),pointer :: Ewald_Vec(:,:)
  integer          :: NVecMax, NSQMAX, NMAX,Boxenanzahl
  real(RK)         :: Kappa, KappaL
  real(RK)         :: USelbstterm
  real(RK)         :: UIntra
  real(RK)         :: UFourier, EVirial
!   real(RK),pointer :: SSin_Fac(:), SCos_Fac(:)
  integer,pointer ::  NBox0(:),NBox1(:),NBox2(:)
   real(RK),pointer :: U_fourierLocal(:)
   real(RK),pointer:: SSin(:),SCos(:)
   real(RK),pointer:: sinfac_s(:,:,:), cosfac_s(:,:,:)
   real(RK),pointer:: sinfac(:),cosfac(:)
!   real(RK),pointer:: sinfac_s_old(:),cosfac_s_old(:)
   real(RK),pointer:: rold(:,:)
   real(RK),pointer:: SSin_Vec(:),SCos_Vec(:)
   real(RK),pointer:: Vec2(:)
   real(RK),pointer:: Faktor(:)
   real(RK),pointer:: HFac(:)
   real(RK),pointer:: distx(:,:,:),disty(:,:,:),distz(:,:,:)
   real(RK),pointer:: VirIntra(:)

  ! SPME parameters
   integer         :: splineorder
   integer         :: gridx
   integer         :: gridy
   integer         :: gridz
   real(RK),pointer::   qgrida(:,:), qgrida_old(:,:), qgridb(:,:)
   integer*8       :: qgrid_forward, qgrid_backward
   real(RK),pointer:: bbtot(:)
   real(RK),pointer:: bsp_arr(:), bsp_modx(:), bsp_mody(:), bsp_modz(:)
   real(RK)        :: EVirialIntra
   real(RK),pointer:: EPotPME(:), VirialPME(:), mm2(:)
!    real(RK),pointer:: spline(:), dspline(:)

   ! Extended ReactionField Method
   real(RK)        :: DebyeLen


#if CONSTR > 0
   integer         :: NCons
   integer,pointer :: Cons1Comp(:)
   integer,pointer :: Cons2Comp(:)
   integer,pointer :: Cons1(:)
   integer,pointer :: Cons2(:)
   real(RK),pointer:: ConsR(:)
   real(RK),pointer:: FCons(:)
   real(RK),pointer:: UCons(:)
   logical         :: consup
#endif

#ifdef ABL
   real(RK),pointer:: AblPS(:,:)
   real(RK),pointer:: AblPE(:,:)
   real(RK),pointer:: AblRhoS(:,:)
   real(RK),pointer:: AblRhoE(:,:)
#endif

  end type TEnsemble

  interface Construct
    module procedure TEnsemble_Construct
  end interface

  interface ConstructSVC
    module procedure TEnsemble_ConstructSVC
  end interface

  interface Destruct
    module procedure TEnsemble_Destruct
  end interface

  interface CreateComponents
    module procedure TEnsemble_CreateComponents
  end interface

  interface DestroyComponents
    module procedure TEnsemble_DestroyComponents
  end interface

  interface CreatePotentials
    module procedure TEnsemble_CreatePotentials
  end interface

  interface DestroyPotentials
    module procedure TEnsemble_DestroyPotentials
  end interface

  interface CreateAccumulators
    module procedure TEnsemble_CreateAccumulators
  end interface

  interface DestroyAccumulators
    module procedure TEnsemble_DestroyAccumulators
  end interface

  interface CalculateNPart
    module procedure TEnsemble_CalculateNPart
  end interface

  interface UpdateBoxLength
    module procedure TEnsemble_UpdateBoxLength
  end interface

  interface UpdateFractions
    module procedure TEnsemble_UpdateFractions
  end interface

  interface FindNSiteMax
    module procedure TEnsemble_FindNSiteMax
  end interface

  interface Allocate
    module procedure TEnsemble_Allocate
  end interface

  interface DeallocateEPot
    module procedure TEnsemble_DeallocateEPot
  end interface

  interface Deallocate
    module procedure TEnsemble_Deallocate
  end interface

  interface CalculateCorr
    module procedure TEnsemble_CalculateCorr
  end interface

  interface InitPositions
    module procedure TEnsemble_InitPositions
  end interface

  interface InitOrientations
    module procedure TEnsemble_InitOrientations
  end interface

  interface InitMolecularDynamics
    module procedure TEnsemble_InitMolecularDynamics
  end interface

  interface InitVelocities
    module procedure TEnsemble_InitVelocities
  end interface

  interface InitIntegrator
    module procedure TEnsemble_InitIntegrator
  end interface

  interface InitIntegratorGear
    module procedure TEnsemble_InitIntegratorGear
  end interface

  interface InitIntegratorLeapFrog
    module procedure TEnsemble_InitIntegratorLeap
  end interface

  interface InitIntegratorVerlet
    module procedure TEnsemble_InitIntegratorVerlet
  end interface

  interface InitIntegratorVV
    module procedure TEnsemble_InitIntegratorVV
  end interface

  interface RemoveNetMomentum
    module procedure TEnsemble_RemoveNetMomentum
  end interface

  interface CalculateEKin
    module procedure TEnsemble_CalculateEKin
  end interface

  interface CheckNPart
    module procedure TEnsemble_CheckNPart
  end interface

  interface ResetEnsemble
    module procedure TEnsemble_ResetEnsemble
  end interface

  interface RunMDStep
    module procedure TEnsemble_RunMDStep
  end interface

  interface RunMCStep
    module procedure TEnsemble_RunMCStep
  end interface

  interface RunSVCStep
    module procedure TEnsemble_RunSVCStep
  end interface

  interface Mol2Atom
    module procedure TEnsemble_Mol2Atom
  end interface

  interface Atom2Mol
    module procedure TEnsemble_Atom2Mol
  end interface

  interface Predict
    module procedure TEnsemble_Predict
  end interface

  interface Correct
    module procedure TEnsemble_Correct
  end interface

  interface PredictGear
    module procedure TEnsemble_PredictGear
  end interface

  interface CorrectGear
    module procedure TEnsemble_CorrectGear
  end interface

  interface PredictLeapFrog
    module procedure TEnsemble_PredictLeapFrog
  end interface

  interface CorrectLeapFrog
    module procedure TEnsemble_CorrectLeapFrog
  end interface

  interface PredictVerlet
    module procedure TEnsemble_PredictVerlet
  end interface

  interface CorrectVerlet
    module procedure TEnsemble_CorrectVerlet
  end interface

  interface PredictVV
    module procedure TEnsemble_PredictVV
  end interface

  interface CorrectVV
    module procedure TEnsemble_CorrectVV
  end interface

  interface Force
    module procedure TEnsemble_Force
  end interface

  interface ChemicalPotential
    module procedure TEnsemble_ChemicalPotential
  end interface

  interface UpdateEnergy
    module procedure TEnsemble_UpdateEnergy
    module procedure TEnsemble_UpdateEnergy1
  end interface

  interface Energy
    module procedure TEnsemble_Energy
    module procedure TEnsemble_Energy1
    module procedure TEnsemble_Energy1_CF
  end interface

  interface GetEnergy
    module procedure TEnsemble_GetEnergy
    module procedure TEnsemble_GetEnergy1
  end interface

  interface GetVirial
    module procedure TEnsemble_GetVirial
  end interface

  interface Move
    module procedure TEnsemble_Move
  end interface

  interface Rotate
    module procedure TEnsemble_Rotate
  end interface

  interface MoveBiased
    module procedure TEnsemble_MoveBiased
  end interface

  interface RotateBiased
    module procedure TEnsemble_RotateBiased
  end interface

  interface ChangeFluct
    module procedure TEnsemble_ChangeFluct
  end interface

  interface Insert
    module procedure TEnsemble_Insert
  end interface

  interface Delete
    module procedure TEnsemble_Delete
  end interface

  interface Move2End
    module procedure TEnsemble_Move2End
  end interface

  interface Resize
    module procedure TEnsemble_Resize
  end interface

  interface ZeroNAttempts
    module procedure TEnsemble_ZeroNAttempts
  end interface

  interface UpdateDisplacements
    module procedure TEnsemble_UpdateDisplacements
  end interface

  interface SaveState
    module procedure TEnsemble_SaveState
  end interface

  interface RestoreState
    module procedure TEnsemble_RestoreState
  end interface

  interface ResultOpen
    module procedure TEnsemble_ResultOpen
  end interface

  interface ResultUpdate
    module procedure TEnsemble_ResultUpdate
  end interface

  interface ResultClose
    module procedure TEnsemble_ResultClose
  end interface

  interface ErrorsUpdate
    module procedure TEnsemble_ErrorsUpdate
  end interface

  interface SVCOutput
    module procedure TEnsemble_SVCOutput
  end interface

  interface VisualOpen
    module procedure TEnsemble_VisualOpen
  end interface

  interface VisualUpdate
    module procedure TEnsemble_VisualUpdate
  end interface

  interface VisualClose
    module procedure TEnsemble_VisualClose
  end interface

  interface RestartSave
    module procedure TEnsemble_RestartSave
  end interface

  interface RestartRead
    module procedure TEnsemble_RestartRead
  end interface

  interface EwaldSelfTerm
    module procedure TEnsemble_EwaldSelfTerm
  end interface

  interface EwaldFourierTerm
    module procedure TEnsemble_EwaldFourierTerm
  end interface

  interface EwaldFourierEnergy
    module procedure TEnsemble_EwaldFourierEnergy
    module procedure TEnsemble_EwaldFourierEnergy1
    module procedure TEnsemble_EwaldFourierEnergy_CF
  end interface

  interface EwaldSelfTerm_Energy
    module procedure TEnsemble_EwaldSelf_Energy
  end interface

  interface Ewald_ChemPotSelf
    module procedure TEnsemble_Ewald_ChemPotSelf
  end interface
! ! ! 
! ! !   interface Ewald_ChemPotFour
! ! !     module procedure TEnsemble_Ewald_ChemPotFour
! ! !   end interface

!   interface PMESelfTerm
!     module procedure TEnsemble_PMESelfTerm
!   end interface

  interface PMEFourierTerm
    module procedure TEnsemble_PMEFourierTerm
  end interface

  interface PMESetup
    module procedure TEnsemble_PME_Setup
  end interface

  interface PMESelfTermMC
    module procedure TEnsemble_PMESelfTerm_MC
  end interface

  interface PMEFourierTermMC
    module procedure TEnsemble_PMEFourierTerm_MC
  end interface

  interface chargegrid_min
   module procedure TEnsemble_PMChargeGrid_min
   end interface

  interface chargegrid_plus
   module procedure TEnsemble_PMChargeGrid_plus
   end interface

#if CONSTR > 0
  interface Constraints
    module procedure TEnsemble_Constraints
  end interface
#endif


contains



!==============================================================!
!  Subroutine TEnsemble_Construct                              !
!==============================================================!

  subroutine TEnsemble_Construct( this, ne )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: ne

    ! Declare local variables
    integer :: i, j
    integer :: stat

    ! Allocate simulation box length
    allocate( this%BoxLength, STAT = stat )
    call AllocationError( stat, 'simulation box length' )

    ! Allocate maximum number of particles
    allocate( this%NPartMax, STAT = stat )
    call AllocationError( stat, 'maximum number of particles' )

    ! Set number of ensemble
    this%EnsembleNumber = ne
    call LogWriteBlank
    write( IOBuffer, '("Reading parameters of ensemble", I3)' ) &
&     this%EnsembleNumber
    call LogWrite

    ! Read temperature
    call FileReadParameter( iounit_params , IdRefTemperature )
    read( IOBuffer, * ) this%RefTemperature
    if( .not. UseReducedUnits ) then
      this%RefTemperature = this%RefTemperature / UnitTemperature
    end if

    ! Read pressure
    if( EnsembleType .eq. EnsembleTypeGE ) then
      call FileReadParameter( iounit_params , IdPressure0 )
      read( IOBuffer, * ) this%RefPressure
      if( .not. UseReducedUnits ) then
        this%RefPressure = this%RefPressure * 1E6_RK / UnitPressure
      end if
    end if
    if( ConstantPressure ) then
      call FileReadParameter( iounit_params , IdRefPressure )
      read( IOBuffer, * ) this%RefPressure
      if( .not. UseReducedUnits ) then
        this%RefPressure = this%RefPressure * 1E6_RK / UnitPressure
      end if
    end if

    ! Read liquid simulation data
    if( EnsembleType .eq. EnsembleTypeGE ) then
      call FileReadParameter( iounit_params , IdLiqDensity )
      read( IOBuffer, * ) this%LiqDensity
      if( .not. UseReducedUnits ) then
        this%LiqDensity = this%LiqDensity / UnitDensity
      end if
      call FileReadParameter( iounit_params , IdVarLiqDensity )
      read( IOBuffer, * ) this%VarLiqDensity
      if( .not. UseReducedUnits ) then
        this%VarLiqDensity = this%VarLiqDensity / UnitDensity
      end if
      call FileReadParameter( iounit_params , IdLiqEnthalpy )
      read( IOBuffer, * ) this%LiqEnthalpy
      if( .not. UseReducedUnits ) then
        this%LiqEnthalpy = this%LiqEnthalpy / ( UnitEnergy * NAvogadro )
      end if
      call FileReadParameter( iounit_params , IdVarLiqEnthalpy )
      read( IOBuffer, * ) this%VarLiqEnthalpy
      if( .not. UseReducedUnits ) then
        this%VarLiqEnthalpy = this%VarLiqEnthalpy / ( UnitEnergy * NAvogadro )
      end if
      call FileReadParameter( iounit_params , IdLiqBetaT )
      read( IOBuffer, * ) this%LiqBetaT
      if( .not. UseReducedUnits ) then
        this%LiqBetaT = this%LiqBetaT * UnitPressure * 1E-6_RK
      end if
      call FileReadParameter( iounit_params , IdVarLiqBetaT )
      read( IOBuffer, * ) this%VarLiqBetaT
      if( .not. UseReducedUnits ) then
        this%VarLiqBetaT = this%VarLiqBetaT * UnitPressure * 1E-6_RK
      end if
      call FileReadParameter( iounit_params , IdLiqdHdP )
      read( IOBuffer, * ) this%LiqdHdP
      if( .not. UseReducedUnits ) then
        this%LiqdHdP = this%LiqdHdP * UnitDensity
      end if
      call FileReadParameter( iounit_params , IdVarLiqdHdP )
      read( IOBuffer, * ) this%VarLiqdHdP
      if( .not. UseReducedUnits ) then
        this%VarLiqdHdP = this%VarLiqdHdP * UnitDensity
      end if
    end if

    ! Read density
    call FileReadParameter( iounit_params , IdRefDensity )
    read( IOBuffer, * ) this%RefDensity
    if( .not. UseReducedUnits ) then
      this%RefDensity = this%RefDensity / UnitDensity
    end if

    ! Update log file
    write( IOBuffer, '("Temperature: ", F9.3, " K")' ) &
&     this%RefTemperature * UnitTemperature
    call LogWrite
    if( ConstantPressure ) then
      write( IOBuffer, '("Pressure: ", F9.3, " MPa")' ) &
&       this%RefPressure * UnitPressure * 1E-6_RK
      call LogWrite
    end if
    if( EnsembleType .eq. EnsembleTypeGE ) then
      write( IOBuffer, '("Pressure0: ", F9.6, " MPa")' ) &
&       this%RefPressure * UnitPressure * 1E-6_RK
      call LogWrite
      write( IOBuffer, '("Liquid density: ", F9.6, " (", F9.6, ") mol/l")' ) &
&       this%LiqDensity * UnitDensity, this%VarLiqDensity * UnitDensity
      call LogWrite
      write( IOBuffer, '("Liquid enthalpy: ", F9.2, " (", F9.2, ") J/mol")' ) &
&       this%LiqEnthalpy * UnitEnergy * NAvogadro, &
&       this%VarLiqEnthalpy * UnitEnergy * NAvogadro
      call LogWrite
      write( IOBuffer, '("Liquid betaT: ", F8.6, "( ", F8.6, ") 1/MPa")' ) &
&       this%LiqBetaT / UnitPressure * 1E6_RK, &
&       this%VarLiqBetaT / UnitPressure * 1E6_RK
      call LogWrite
      write( IOBuffer, '("Liquid dHdP: ", F8.6, "( ", F8.6, ") l/mol")' ) &
&       this%LiqdHdP / UnitDensity, this%VarLiqdHdP / UnitDensity
      call LogWrite
    end if
    write( IOBuffer, '("Density: ", F9.3, " mol/l")' ) &
&     this%RefDensity * UnitDensity
    call LogWrite
    write( IOBuffer, '("Reduced temperature: ", F9.6)' ) this%RefTemperature
    call LogWrite
    if( ConstantPressure ) then
      write( IOBuffer, '("Reduced pressure: ", F9.6)' ) this%RefPressure
      call LogWrite
    end if
    if( EnsembleType .eq. EnsembleTypeGE ) then
      write( IOBuffer, '("Reduced pressure0: ", F9.6)' ) this%RefPressure
      call LogWrite
      write( IOBuffer, &
&       '("Reduced liquid density: ", F9.6, " (", F9.6, ")")' ) &
&       this%LiqDensity, this%VarLiqDensity
      call LogWrite
      write( IOBuffer, &
&       '("Reduced liquid enthalpy: ", F9.4, " (", F9.4, ")")' ) &
&       this%LiqEnthalpy, this%VarLiqEnthalpy
      call LogWrite
      write( IOBuffer, '("Reduced liquid betaT: ", F8.6, "( ", F8.6, ")")' ) &
&       this%LiqBetaT, this%VarLiqBetaT
      call LogWrite
      write( IOBuffer, '("Reduced liquid dHdP: ", F8.4, "( ", F8.4, ")")' ) &
&       this%LiqdHdP, this%VarLiqdHdP
      call LogWrite
    end if
    write( IOBuffer, '("Reduced density: ", F9.6)' ) this%RefDensity
    call LogWrite

    ! Read mass of piston
    if( SimulationType .eq. MolecularDynamics .and. ConstantPressure ) then
      call FileReadParameter( iounit_params , IdPistonMass )
      read( IOBuffer, * ) this%PistonMass
!       if( .not. UseReducedUnits ) then
!         this%PistonMass = this%PistonMass / UnitMass * UnitLength**4
!       end if

      write( IOBuffer, '("Mass of piston: ", F12.9)' ) this%PistonMass
      call LogWrite
    end if

    ! Read initial number of particles in ensemble
    call FileReadParameter( iounit_params , IdNPart )
    read( IOBuffer, * ) this%NPart
    if( EnsembleType .eq. EnsembleTypeGE .or. &
&       EnsembleType .eq. EnsembleTypeHA ) then
      this%NPartInitial = this%NPart
      this%NPartLBound = int( real( this%NPart, RK ) / 1.2_RK )
      this%NPartUBound = int( real( this%NPart, RK ) * 1.2_RK )
    end if

    ! Read number of components in ensemble
    call FileReadParameter( iounit_params , IdNComponents )
    read( IOBuffer, * ) this%NComponents
    write( IOBuffer, '("Number of components:", I3)' ) this%NComponents
    call LogWrite
    if( this%NComponents <= 0 ) then
      write( ErrorBuffer, &
&       '("There must be at least 1 component in ensemble", I2)' ) &
&       this%EnsembleNumber
      call Error
    end if
    if( this%NComponents > 999 ) &
&     call Error( 'Cannot work with more than 999 components on '//Hardware )

#if MPI_VER > 0
    if ((LongRange .eq. Ewald) .or. (LongRange .eq. PME))then
      allocate(this%NBox0(NProcs),STAT=stat)
      allocate(this%NBox1(NProcs),STAT=stat)
      allocate(this%NBox2(NProcs),STAT=stat)
    end if
#endif
    ! Create components
    call CreateComponents( this )

    ! Calculate initial number of particles of each component
    call CalculateNPart( this )

    ! Calculate maximum numbers of sites in components
    call FindNSiteMax( this )

    ! Charge sites need center of mass cutoff
    if(( CutoffMode .ne. CenterofMass ) .and. ( this%NChargeMax .gt. 0 )) &
&     call Error( 'Center of mass correction is necessary when using charges')

    ! Allocate arrays
    call Allocate( this )

    ! Read scale coefficients for LJ126 epsilon and sigma
    this%ScaleSigma(:, :) = 1._RK
    this%ScaleEpsilon(:, :) = 1._RK
    do i = 1, this%NRealComponents - 1
      do j = i + 1, this%NRealComponents
        call FileReadParameter( iounit_params , IdScaleSigma )
        read( IOBuffer, * ) this%ScaleSigma(i, j)
        if( i /= j ) this%ScaleSigma(j, i) = this%ScaleSigma(i, j)
        call FileReadParameter( iounit_params , IdScaleEpsilon )
        read( IOBuffer, * ) this%ScaleEpsilon(i, j)
        if( i /= j ) this%ScaleEpsilon(j, i) = this%ScaleEpsilon(i, j)
        write( IOBuffer, &
&         '(A, "-", A, " Lennard-Jones interaction:  eta =", F6.3, ", xi =", F6.3)' ) &
&         trim( this%Component(i)%PotModFileName ), &
&         trim( this%Component(j)%PotModFileName ), &
&         this%ScaleSigma(i, j), this%ScaleEpsilon(i, j)
        call LogWrite
      end do
    end do

    ! Read cutoff radii
    this%RCutoffLJ126LJ126 = 0._RK
    this%RCutoffDipoleDipole = 0._RK
    this%RCutoffDipoleQuadrupole = 0._RK
    this%RCutoffQuadrupoleQuadrupole = 0._RK
    if( CutoffMode .eq. CenterofMass ) then
      call FileReadParameter( iounit_params , IdRCutoffCOM )
      read( IOBuffer, * ) this%RCutoffLJ126LJ126
      if (this%RCutoffLJ126LJ126 < 0._RK) then
        this%RCutoffLJ126LJ126 = 0.9*0.5*(this%NPart / (NAvogadro*this%RefDensity*UnitDensity*1000))**(1._RK/3._RK)/UnitLength
      end if
      write( IOBuffer, '("Reduced center of mass cutoff radius: ", F6.3)' ) &
&       this%RCutoffLJ126LJ126
      call LogWrite
      this%RCutoffDipoleDipole = this%RCutoffLJ126LJ126
      this%RCutoffDipoleQuadrupole = this%RCutoffLJ126LJ126
      this%RCutoffQuadrupoleQuadrupole = this%RCutoffLJ126LJ126
    else
      if( this%NLJ126Max > 0 ) then
        call FileReadParameter( iounit_params , IdRCutoffLJ126LJ126 )
        read( IOBuffer, * ) this%RCutoffLJ126LJ126
        write( IOBuffer, &
&         '("Lennard-Jones cutoff radius: ", F6.3, " sigma")' ) &
&         this%RCutoffLJ126LJ126
        call LogWrite
      end if
      if( this%NDipoleMax > 0 ) then
        call FileReadParameter( iounit_params , IdRCutoffDipoleDipole )
        read( IOBuffer, * ) this%RCutoffDipoleDipole
        write( IOBuffer, '("Reduced dipole-dipole cutoff radius: ", F8.3)' ) &
&         this%RCutoffDipoleDipole
        call LogWrite
        if( this%NQuadrupoleMax > 0 ) then
          call FileReadParameter( iounit_params , IdRCutoffDipoleQuadrupole )
          read( IOBuffer, * ) this%RCutoffDipoleQuadrupole
          write( IOBuffer, &
&           '("Reduced dipole-quadrupole cutoff radius: ", F8.3)' ) &
&           this%RCutoffDipoleQuadrupole
          call LogWrite
        end if
      end if
      if( this%NQuadrupoleMax > 0 ) then
        call FileReadParameter( iounit_params , IdRCutoffQuadrupoleQuadrupole )
        read( IOBuffer, * ) this%RCutoffQuadrupoleQuadrupole
        write( IOBuffer, &
&         '("Reduced quadrupole-quadrupole cutoff radius: ", F8.3)' ) &
&         this%RCutoffQuadrupoleQuadrupole
        call LogWrite
      end if
    end if

    ! Read characteristic dielectric constant
    this%RFEpsilon = 0._RK
    if(( this%NDipoleMax > 0 ) .or. ( this%NChargeMax > 0 )) then
      call FileReadParameter( iounit_params , IdRFEpsilon )
      read( IOBuffer, * ) this%RFEpsilon
      write( IOBuffer, '("Characteristic dielectric constant: ", E16.9)' ) &
&       this%RFEpsilon
      call LogWrite
    end if


#if CONSTR > 0
    write( IOBuffer, '("CONSTRAINED DYNAMICS")' )
    call LogWrite

    call FileReadParameter( iounit_params , IdNCons )
    read( IOBuffer, * ) this%NCons
    write( IOBuffer, '("Number of Constrained Molecules:", I3)' ) this%NCons
    call LogWrite

    allocate( this%Cons1Comp( this%NCons), STAT = stat )
       if(stat >0) write(*,*) 'Allocation Error Cons1Comp'
!     call AllocationError( stat, 'Constraint Molecules', this%NCons )
    allocate( this%Cons2Comp( this%NCons), STAT = stat )
       if(stat >0) write(*,*) 'Allocation Error Cons2Comp'
!     call AllocationError( stat, 'Constraint Molecules', this%NCons )
    allocate( this%Cons1( this%NCons), STAT = stat )
           if(stat >0) write(*,*) 'Allocation Error Cons1'
! call AllocationError( stat, 'Constraint Molecules', this%NCons )
    allocate( this%Cons2( this%NCons), STAT = stat )
         if(stat >0) write(*,*) 'Allocation Error Cons2'
!   call AllocationError( stat, 'Constraint Molecules', this%NCons )
    allocate( this%ConsR( this%NCons), STAT = stat )
         if(stat >0) write(*,*) 'Allocation Error ConsR'
!   call AllocationError( stat, 'Constraint Molecules', this%NCons )
    allocate( this%FCons( this%NCons), STAT = stat )
         if(stat >0) write(*,*) 'Allocation Error FCons'
!   call AllocationError( stat, 'Constraint Molecules', this%NCons )
    allocate( this%UCons( this%NCons), STAT = stat )
         if(stat >0) write(*,*) 'Allocation Error UCons'
!   call AllocationError( stat, 'Constraint Molecules', this%NCons )

   DO i=1,this%NCons,1
    call FileReadParameter( iounit_params , IdCons1Comp )
    read( IOBuffer, * ) this%Cons1Comp(i)
    write( IOBuffer, '("Constrained Mol Typ 1:", I3)' ) this%Cons1Comp(i)
    call LogWrite

    call FileReadParameter( iounit_params , IdCons1 )
    read( IOBuffer, * ) this%Cons1(i)
    write( IOBuffer, '("Constrained Mol 1:", I3)' ) this%Cons1(i)
    call LogWrite

    call FileReadParameter( iounit_params , IdCons2Comp )
    read( IOBuffer, * ) this%Cons2Comp(i)
    write( IOBuffer, '("Constrained Mol Typ 2:", I3)' ) this%Cons2Comp(i)
    call LogWrite

    call FileReadParameter( iounit_params , IdCons2 )
    read( IOBuffer, * ) this%Cons2(i)
    write( IOBuffer, '("Constrained Mol 2:", I3)' ) this%Cons2(i)
    call LogWrite

    call FileReadParameter( iounit_params , IdConsR )
    read( IOBuffer, * ) this%ConsR(i)
    write( IOBuffer, '("Constrained Mol Distance:", F6.3)' ) this%ConsR(i)
    call LogWrite

    this%ConsR(i) = this%ConsR(i) * Angstroem
    this%ConsR(i) = this%ConsR(i) / UnitLength
! Use only squared!
    this%ConsR(i) = this%ConsR(i) * this%ConsR(i)

   END DO

! REduce the number of degrees of freedom in the system
    this%NDF = this%NDF - this%NCons

#endif


    ! Create potentials
    call CreatePotentials( this )

    ! Calculate long-range corrections
    call CalculateCorr( this )
    write( IOBuffer, &
&     '("Cutoff correction to potential energy from LJ:", F12.8)' ) &
&     this%EPotCorrLJ * NProcs / this%NPart
    call LogWrite
    write( IOBuffer, &
&     '("Cutoff correction to pressure from LJ        :", F12.8)' ) &
&     this%VirialCorrLJ * NProcs / this%NPart
    call LogWrite

    do i = 1, this%NRealComponents
      write( IOBuffer, &
&       '("Cutoff correction to chem. pot. of ", A, " from LJ:", F12.8)' ) &
&       trim( this%Component(i)%PotModFileName ), &
&       this%Component(i)%EPotTestCorrLJ
      call LogWrite
    end do

    write( IOBuffer, &
&     '("Cutoff correction to potential energy from reaction field:", F12.8)' ) &
&     this%EPotCorrRF * NProcs / this%NPart
    call LogWrite

    do i = 1, this%NRealComponents
      write( IOBuffer, &
&       '("Cutoff correction to chem. pot. of ", A, " from reaction field:", F12.8)' ) &
&       trim( this%Component(i)%PotModFileName ), &
&       this%Component(i)%EPotTestCorrRF
      call LogWrite
    end do

    ! Calculate maximum cutoff radius
!     this%RCutoffMax2 = 0._RK       ! Bug
    if( this%NDipoleMax > 0 ) then
      this%RCutoffMax2 = max(this%RCutoffMax2, &
&       2._RK * this%RCutoffDipoleDipole )
      if( this%NQuadrupoleMax > 0 ) then
        this%RCutoffMax2 = max(this%RCutoffMax2, &
&         2._RK * this%RCutoffDipoleQuadrupole )
      end if
    end if
    if( this%NQuadrupoleMax > 0 ) then
      this%RCutoffMax2 = max(this%RCutoffMax2, &
&       2._RK * this%RCutoffQuadrupoleQuadrupole )
    end if

    if( .not. Restart ) then

      ! Update all BoxLength-dependent constants
      call UpdateBoxLength( this )

      ! Set initial positions of particles in simulation box
      call InitPositions( this )

      ! Set initial orientations of particles in simulation box
      call InitOrientations( this )

      ! Calculate initial energies for the Ewald Summation
      if (LongRange .eq. Ewald) then
         if (this%KappaL .eq. 0.) then
            this%Kappa = sqrt(PI) * (4.0*this%NPart / this%Volume0**2)**(1./6.)
            this%KappaL = this%Kappa*this%BoxLength
         else
            this%Kappa = this%KappaL/this%BoxLength   !Boxlength bereits normiert
         end if
         do i=1,this%NComponents
           do j=1,this%NCOmponents
             this%Interaction(i,j)%Kappa = this%Kappa
           end do
         end do
         call EwaldSelfTerm( this )
        ! Memory Allocation for Ewald Summation
         allocate(this%Faktor(this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error Faktor'
         allocate(this%U_fourierLocal(this%BoxenAnzahl),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error U_fourier'
         allocate(this%SSin(this%BoxenAnzahl),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error SSin'
         allocate(this%SCos(this%BoxenAnzahl),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error SCos'
         allocate(this%sinfac(this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error sinfac'
         allocate(this%cosfac(this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error cosfac'
         allocate(this%SSin_Vec(this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error SSin_Vec'
         allocate(this%SCos_Vec(this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error SCos_Vec'
         allocate(this%sinfac_s(this%NComponents,5,this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error sin_facs'
         allocate(this%cosfac_s(this%NComponents,5,this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error cosfac_s'
!          allocate(this%sinfac_s_old(5),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error sinfac'
!          allocate(this%cosfac_s_old(5),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error cosfac'
         allocate(this%rold(5,3),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error rold'
         allocate(this%Vec2(this%BoxenAnzahl),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error Vec2'
         allocate(this%HFac(this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error HFac'
!          allocate(this%distx(this%NComponents,5,this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error distx'
!          allocate(this%disty(this%NComponents,5,this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error disty'
!          allocate(this%distz(this%NComponents,5,this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error distz'
!          allocate(this%VirIntra(this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error VirIntra'
         allocate(this%distx(this%NComponents,5,this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error distx'
         allocate(this%disty(this%NComponents,5,this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error disty'
         allocate(this%distz(this%NComponents,5,this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error distz'
         allocate(this%VirIntra(this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error VirIntra'

! Chemical Potential
         DO i=1,this%NComponents
           call Ewald_ChemPotSelf(this,i)
         END DO

      else if (LongRange .eq. PME) then

         if (this%KappaL .eq. 0.) then
            this%Kappa = sqrt(PI) * (4.0*this%NPart / this%Volume0**2)**(1./6.)
            this%KappaL = this%Kappa*this%BoxLength
         else
            this%Kappa = this%KappaL/this%BoxLength   !Boxlength bereits normiert
         end if
         do i=1,this%NComponents
           do j=1,this%NCOmponents
             this%Interaction(i,j)%Kappa = this%Kappa
           end do
         end do
         allocate(this%bsp_arr(this%gridx+5),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error bsp_arr'
         allocate(this%bsp_modx(this%gridx+5),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error bsp_modx'
         allocate(this%bsp_mody(this%gridy+5),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error bsp_mody'
         allocate(this%bsp_modz(this%gridz+5),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error bsp_modz'

! Setup SPME calculation
         call PMESetup (this)
       else if (LongRange .eq. ExtRField) then
         do i=1,this%NComponents
           do j=1,this%NCOmponents
             this%Interaction(i,j)%DebyeLen = this%DebyeLen
           end do
         end do
      end if

      if( SimulationType .eq. MolecularDynamics &
&       .and. .not. MCOverlapReduction ) then

        ! Initialize molecular dynamics simulation
        call InitMolecularDynamics( this, .false. )

      else

        ! Set temperature
        this%Temperature = this%RefTemperature

        ! Convert molecular coordinates to atom positions
        call Mol2Atom( this )

        ! Set all potential energy matrices
        call Energy( this, this%EPot )
        call UpdateEnergy( this )

        ! Set initial values of maximum allowed MC displacements
        this%DispVol = DispVolStart
        do i = 1, this%NRealComponents
          this%Component(i)%DispTran = DispTranStart
          this%Component(i)%DispRot = DispRotStart
        end do

      end if
    else        ! Restart
      ! Calculate initial energies for the Ewald Summation
      if (LongRange .eq. Ewald) then
!          write (*,*) 'this%BoxLength wurde festgelegt in Zeile ms2_ensemble 958'
! 	 this%BoxLength=3.6445654373
!          this%BoxLength = 1.0_RK
!          this%Kappa = 0.0_RK
!          this%Kappa = this%KappaL/this%BoxLength   !Boxlength bereits normiert
!          do i=1,this%NCOmponents
!            do j=1,this%NCOmponents
!              this%Interaction(i,j)%Kappa = this%Kappa
!            end do
!          end do
!          call EwaldSelfTerm( this )
        ! Memory Allocation for Ewald Summation
!          allocate(this%Faktor(this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error Faktor'
!          allocate(this%U_fourierLocal(this%BoxenAnzahl),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error U_fourier'
!          allocate(this%SSin(this%BoxenAnzahl),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error SSin'
!          allocate(this%SCos(this%BoxenAnzahl),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error SCos'
!          allocate(this%sinfac(this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error sinfac'
!          allocate(this%cosfac(this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error cosfac'
!          allocate(this%SSin_Vec(this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error SSin_Vec'
!          allocate(this%SCos_Vec(this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error SCos_Vec'
!          allocate(this%sinfac_s(this%NComponents,5,this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error sin_facs'
!          allocate(this%cosfac_s(this%NComponents,5,this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error cosfac_s'
!          allocate(this%Vec2(this%BoxenAnzahl),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error Vec2'
!          allocate(this%HFac(this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error HFac'
!          allocate(this%distx(this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error distx'
!          allocate(this%disty(this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error disty'
!          allocate(this%distz(this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error distz'
!          allocate(this%VirIntra(this%NPart),STAT=stat)
!          if(stat >0) write(*,*) 'Allocation Error VirIntra'

      end if
    end if         ! Restart

    ! Set I/O unit numbers
    i = FilesPerEnsemble * this%EnsembleNumber
    this%iounit_result = iounit_result + i
    this%iounit_runave = iounit_runave + i
    this%iounit_errors = iounit_errors + i
    this%iounit_visual = iounit_visual + i


  end subroutine TEnsemble_Construct



!==============================================================!
!  Subroutine TEnsemble_ConstructSVC                           !
!==============================================================!

  subroutine TEnsemble_ConstructSVC( this, ne )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: ne

    ! Declare local variables
    integer                   :: i, j
    integer                   :: stat
    type(TComponent), pointer :: pc
    character(FileNameLength) :: PotModFileName
    real(RK)                  :: scaleSigma, scaleEpsilon

    ! Allocate simulation box length
    allocate( this%BoxLength, STAT = stat )
    call AllocationError( stat, 'simulation box length' )

    ! Allocate maximum number of particles
    allocate( this%NPartMax, STAT = stat )
    call AllocationError( stat, 'maximum number of particles' )

    ! Set number of ensemble
    this%EnsembleNumber = ne
    call LogWriteBlank
    write( IOBuffer, '("Reading parameters of ensemble", I3)' ) &
&     this%EnsembleNumber
    call LogWrite

    ! Read temperature
    call FileReadParameter( iounit_params , IdRefTemperature )
    read( IOBuffer, * ) this%Temperature
    if( .not. UseReducedUnits ) then
      this%Temperature = this%Temperature / UnitTemperature
    end if

    ! Update log file
    write( IOBuffer, '("Temperature: ", F9.3, " K")' ) &
&     this%Temperature * UnitTemperature
    call LogWrite
    write( IOBuffer, '("Reduced temperature: ", F9.6)' ) this%Temperature
    call LogWrite

    ! Read number of components in ensemble
    call FileReadParameter( iounit_params , IdNComponents )
    read( IOBuffer, * ) this%NComponents
    write( IOBuffer, '("Number of components:", I3)' ) this%NComponents
    call LogWrite
    if( this%NComponents <= 0 ) then
      write( ErrorBuffer, &
&       '("There must be at least 1 component in ensemble", I2)' ) &
&       this%EnsembleNumber
      call Error
    end if
    if( this%NComponents > 999 ) &
&     call Error( 'Cannot work with more than 999 components on '//Hardware )

    ! Create components
    this%NComponents = 2 * this%NComponents
    allocate( this%Component(this%NComponents), STAT = stat )
    call AllocationError( stat, 'components', this%NComponents )
    do i = 1, this%NComponents, 2

      ! Read file name for potential model
      call FileReadParameter( iounit_params , IdPotModFileName )
      read( IOBuffer, * ) PotModFileName

      call Construct( this%Component(i), PotModFileName )
      call Construct( this%Component(i+1), PotModFileName )

    end do

    ! Calculate number of particles in each process
    do i = 1, this%NComponents
      pc => this%Component(i)
      pc%NPart = NOrient
      pc%NPart1 = 1 + (pc%NPart - 1) / NProcs
      pc%NPart0 = 1 + pc%NPart1 * NProc
      pc%NPart2 = min( pc%NPart0 + pc%NPart1 - 1, pc%NPart )
      pc%NPart1 = pc%NPart2 - pc%NPart0 + 1
    end do
    this%NPartMax = NOrient
    this%NTestMax = 0

    ! Calculate maximum numbers of sites in components
    call FindNSiteMax( this )

    ! Allocate arrays
    call Allocate( this )

    ! Read scale coefficients for LJ126 epsilon and sigma
    this%ScaleSigma(:, :) = 1._RK
    this%ScaleEpsilon(:, :) = 1._RK
    do i = 1, this%NComponents - 2, 2
      do j = i + 2, this%NComponents, 2
        call FileReadParameter( iounit_params , IdScaleSigma )
        read( IOBuffer, * ) scaleSigma
        this%ScaleSigma(i:i+1, j:j+1) = scaleSigma
        if( i /= j ) this%ScaleSigma(j:j+1, i:i+1) = scaleSigma
        call FileReadParameter( iounit_params , IdScaleEpsilon )
        read( IOBuffer, * ) scaleEpsilon
        this%ScaleEpsilon(i:i+1, j:j+1) = scaleEpsilon
        if( i /= j ) this%ScaleEpsilon(j:j+1, i:i+1) = scaleEpsilon
        write( IOBuffer, &
&         '(A, "-", A, " Lennard-Jones interaction:  eta =", F6.3, ", xi =", F6.3)' ) &
&         trim( this%Component(i)%PotModFileName ), &
&         trim( this%Component(j)%PotModFileName ), &
&         this%ScaleSigma(i, j), this%ScaleEpsilon(i, j)
        call LogWrite
      end do
    end do

    ! Set cutoff radii
    this%RCutoffLJ126LJ126 = MaxRadius

    ! Disable reaction field
    this%RFEpsilon = 0._RK

    ! Create potentials
    call CreatePotentials( this )

    ! Calculate long-range corrections
    do i = 1, this%NComponents, 2
      do j = i + 1, this%NComponents, 2
        this%Interaction(i, j)%EPotCorrLJ = &
&         sum( this%Interaction(i, j)%PotLJ126LJ126(:, :)%EPotCorr )
        write( IOBuffer, &
&         '("Cutoff correction to SVC of ", A, "-", A, " from LJ:", F12.8)' ) &
&         trim( this%Component(i)%Molecule%PotModFileName ), &
&         trim( this%Component(j)%Molecule%PotModFileName ), &
&         .5_RK * this%Interaction(i, j)%EPotCorrLJ / this%Temperature
        call LogWrite
      end do
    end do

    ! Update all BoxLength-dependent constants
    call UpdateBoxLength( this )

    ! Set initial positions of particles
    do i = 1, this%NComponents
      this%Component(i)%P0 = 0._RK
    end do

    ! Set initial orientations of particles
    call InitOrientations( this )

    ! Convert molecular coordinates to atom positions
    call Mol2Atom( this )

    ! Set I/O unit numbers
    i = FilesPerEnsemble * this%EnsembleNumber
    this%iounit_result = iounit_result + i
    this%iounit_runave = iounit_runave + i
    this%iounit_errors = iounit_errors + i
    this%iounit_visual = iounit_visual + i

  end subroutine TEnsemble_ConstructSVC



!==============================================================!
!  Subroutine TEnsemble_Destruct                               !
!==============================================================!

  subroutine TEnsemble_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif
    include 'fftw3.f'

    ! Declare arguments
    type(TEnsemble) :: this

    ! Destroy potentials
    call DestroyPotentials( this )

    ! Deallocate arrays
    call Deallocate( this )

    ! Destroy components
    call DestroyComponents( this )

    ! Deallocate maximum number of particles
    if( associated( this%NPartMax ) ) then
      deallocate( this%NPartMax )
    end if

    ! Deallocate simulation box length
    if( associated( this%BoxLength ) ) then
      deallocate( this%BoxLength )
    end if

    if (LongRange .eq. Ewald) then
      if ( associated ( this%Ewald_Prefac ) ) then 
         deallocate( this%Ewald_Prefac )
      end if
      if ( associated ( this%Ewald_Vec ) ) then 
         deallocate( this%Ewald_Vec )
      end if
#if MPI_VER > 0
      if ( associated ( this%NBox0 ) ) then 
         deallocate( this%NBox0 )
      end if
      if ( associated ( this%NBox1 ) ) then 
         deallocate( this%NBox1 )
      end if
      if ( associated ( this%NBox2 ) ) then 
         deallocate( this%NBox2 )
      end if
#endif
      if ( associated ( this%Faktor ) ) then 
         deallocate( this%Faktor )
      end if
      if ( associated ( this%U_fourierLocal ) ) then 
         deallocate( this%U_fourierLocal )
      end if
      if ( associated ( this%SSin ) ) then 
         deallocate( this%SSin )
      end if
      if ( associated ( this%SCos ) ) then 
         deallocate( this%SCos )
      end if
      if ( associated ( this%sinfac ) ) then 
         deallocate( this%sinfac )
      end if
      if ( associated ( this%cosfac ) ) then 
         deallocate( this%cosfac )
      end if
      if ( associated ( this%SSin_Vec ) ) then 
         deallocate( this%SSin_Vec )
      end if
      if ( associated ( this%SCos_Vec ) ) then 
         deallocate( this%SCos_Vec )
      end if
      if ( associated ( this%sinfac_s ) ) then 
         deallocate( this%sinfac_s )
      end if
      if ( associated ( this%cosfac_s ) ) then 
         deallocate( this%cosfac_s )
      end if
      if ( associated ( this%rold ) ) then 
         deallocate( this%rold )
      end if
      if ( associated ( this%Vec2 ) ) then 
         deallocate( this%Vec2 )
      end if
      if ( associated ( this%HFac ) ) then 
         deallocate( this%HFac )
      end if
      if ( associated ( this%distx ) ) then 
         deallocate( this%distx )
      end if
      if ( associated ( this%disty ) ) then 
         deallocate( this%disty )
      end if
      if ( associated ( this%distz ) ) then 
         deallocate( this%distz )
      end if
      if ( associated ( this%VirIntra ) ) then 
         deallocate( this%VirIntra )
      end if
    end if

    if (LongRange .eq. PME) then
      call dfftw_destroy_plan(this%qgrid_forward)
      call dfftw_destroy_plan(this%qgrid_backward)
      if ( associated ( this%qgrida ) ) then 
         deallocate( this%qgrida )
      end if
      if ( associated ( this%qgrida_old ) ) then 
         deallocate( this%qgrida_old )
      end if
      if ( associated ( this%bsp_arr ) ) then 
         deallocate( this%bsp_arr )
      end if
      if ( associated ( this%bsp_modx ) ) then 
         deallocate( this%bsp_modx )
      end if
      if ( associated ( this%bsp_mody ) ) then 
         deallocate( this%bsp_mody )
      end if
      if ( associated ( this%bsp_modz ) ) then 
         deallocate( this%bsp_modz )
      end if
    endif

  end subroutine TEnsemble_Destruct



!==============================================================!
!  Subroutine TEnsemble_CreateComponents                       !
!==============================================================!

  subroutine TEnsemble_CreateComponents( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i, j, nfluct, ncomp
    integer :: stat
    type(TComponent), pointer :: reallocate(:)

    ! Create components
    ncomp = this%NComponents
    allocate( this%Component(ncomp), STAT = stat )
    call AllocationError( stat, 'components', ncomp )
    do i = 1, ncomp
      call Construct( this%Component(i), i )
    end do

    ! Create components for fluctuating particle states
    this%NRealComponents = ncomp
    do i = 1, this%NRealComponents
      if( this%Component(i)%ChemPotMethod .eq. ChemPotMethodGradIns ) then
        nfluct = this%Component(i)%Molecule%NFluct
        ncomp = ncomp + nfluct

        ! Reallocate component array
        allocate( reallocate(ncomp), STAT = stat )
        call AllocationError( stat, 'components', ncomp )
        reallocate( 1:size(this%Component) ) = this%Component(:)
        deallocate( this%Component )
        this%Component => reallocate

        do j = 1, nfluct
          call Construct( this%Component(ncomp - nfluct + j), &
&           this%Component(i), j )
          this%Component(i)%NFluctComp(j) = ncomp - nfluct + j
        end do
        this%Component(i)%NFluctComp(0) = i
      end if
    end do

    ! Set new number of components (including fluctuating particles)
    this%NComponents = ncomp

  end subroutine TEnsemble_CreateComponents



!==============================================================!
!  Subroutine TEnsemble_DestroyComponents                      !
!==============================================================!

  subroutine TEnsemble_DestroyComponents( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Destroy components
    do i = 1, this%NComponents
      call Destruct( this%Component(i) )
    end do
    if( associated( this%Component ) ) then
      deallocate( this%Component )
    end if

  end subroutine TEnsemble_DestroyComponents



!==============================================================!
!  Subroutine TEnsemble_CreatePotentials                       !
!==============================================================!

  subroutine TEnsemble_CreatePotentials( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i, j
    integer :: stat

    ! Create interactions
    allocate( this%Interaction( &
&     this%NComponents, this%NComponents ), STAT = stat )
    call AllocationError( stat, 'components', this%NComponents )
    if( SimulationType .eq. SecondVirialCoeff ) then
      do i = 1, this%NComponents, 2
        do j = i + 1, this%NComponents, 2
          call Construct( &
&           this%Interaction(i, j), i, j, &
&           this%Component(i), this%Component(j), &
&           this%RCutoffLJ126LJ126, &
&           this%RCutoffDipoleDipole, &
&           this%RCutoffDipoleQuadrupole, &
&           this%RCutoffQuadrupoleQuadrupole, &
&           this%ScaleSigma(i, j), &
&           this%ScaleEpsilon(i, j), &
&           this%RFEpsilon )
          call UpdateBoxLength( this%Interaction( j, i ), this%BoxLength )
        end do
      end do
    else
      do i = 1, this%NComponents
        do j = 1, this%NComponents
          call Construct( &
&           this%Interaction(i, j), i, j, &
&           this%Component(i), this%Component(j), &
&           this%RCutoffLJ126LJ126, &
&           this%RCutoffDipoleDipole, &
&           this%RCutoffDipoleQuadrupole, &
&           this%RCutoffQuadrupoleQuadrupole, &
&           this%ScaleSigma(i, j), &
&           this%ScaleEpsilon(i, j), &
&           this%RFEpsilon )
        end do
      end do
    end if

  end subroutine TEnsemble_CreatePotentials



!==============================================================!
!  Subroutine TEnsemble_DestroyPotentials                      !
!==============================================================!

  subroutine TEnsemble_DestroyPotentials( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i, j

    ! Destroy interactions
    if( SimulationType .eq. SecondVirialCoeff ) then
      do i = 1, this%NComponents, 2
        do j = i + 1, this%NComponents, 2
          call Destruct( this%Interaction(i, j) )
        end do
      end do
    else
      do i = 1, this%NComponents
        do j = 1, this%NComponents
          call Destruct( this%Interaction(i, j) )
        end do
      end do
    end if
    if( associated( this%Interaction ) ) then
      deallocate( this%Interaction )
    end if

  end subroutine TEnsemble_DestroyPotentials



!==============================================================!
!  Subroutine TEnsemble_CreateAccumulators                     !
!==============================================================!

  subroutine TEnsemble_CreateAccumulators( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Construct accumulators
    if( .not. SimulationType .eq. SecondVirialCoeff ) then

      ! 1.) Basic sums
      call Construct( this%SumPressure, .false. )
      call Construct( this%SumDensity, .false. )
      call Construct( this%SumTemperature, .false. )
      call Construct( this%SumEPot, .false. )
      call Construct( this%SumEnthalpy, .false. )
      call Construct( this%SumVolume, .false. )
      call Construct( this%SumVirial, .false. )
      if( EnsembleType .eq. EnsembleTypeGE .or. &
&         EnsembleType .eq. EnsembleTypeHA ) then
        call Construct( this%SumNPart, .false. )
      end if

      ! 2.) Combined sums
      call Construct( this%SumEPotSquared, .false. )
      call Construct( this%SumEPotV, .false. )
      call Construct( this%SumEPotVirial, .false. )
      call Construct( this%SumEnthalpySquared, .false. )
      call Construct( this%SumEnthalpyV, .false. )
      call Construct( this%SumVolumeSquared, .false. )

      ! 3.) Derived sums
      call Construct( this%SumBetaT, .true. )
      call Construct( this%SumdHdP, .true. )
      call Construct( this%SumdUdV, .true. )
      call Construct( this%SumCV, .true. )
      call Construct( this%SumCP, .true. )
      call Construct( this%SumAlphaP, .true. )

      do i = 1, this%NRealComponents
        call CreateAccumulators( this%Component(i) )
      end do

    end if

  end subroutine TEnsemble_CreateAccumulators



!==============================================================!
!  Subroutine TEnsemble_DestroyAccumulators                    !
!==============================================================!

  subroutine TEnsemble_DestroyAccumulators( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Destruct accumulators
    ! 1.) Basic sums
    call Destruct( this%SumPressure )
    call Destruct( this%SumDensity )
    call Destruct( this%SumTemperature )
    call Destruct( this%SumEPot )
    call Destruct( this%SumEnthalpy )
    call Destruct( this%SumVolume )
    call Destruct( this%SumVirial )
    if( EnsembleType .eq. EnsembleTypeGE .or. &
&       EnsembleType .eq. EnsembleTypeHA ) then
      call Destruct( this%SumNPart )
    end if

    ! 2.) Combined sums
    call Destruct( this%SumEPotSquared )
    call Destruct( this%SumEPotV )
    call Destruct( this%SumEPotVirial )
    call Destruct( this%SumEnthalpySquared )
    call Destruct( this%SumEnthalpyV )
    call Destruct( this%SumVolumeSquared )

    ! 3.) Derived sums
    call Destruct( this%SumBetaT )
    call Destruct( this%SumdHdP )
    call Destruct( this%SumdUdV )
    call Destruct( this%SumCV )
    call Destruct( this%SumCP )
    call Destruct( this%SumAlphaP )

    do i = 1, this%NRealComponents
      call DestroyAccumulators( this%Component(i) )
    end do


  end subroutine TEnsemble_DestroyAccumulators



!==============================================================!
!  Subroutine TEnsemble_CalculateNPart                         !
!==============================================================!

  subroutine TEnsemble_CalculateNPart( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK)                  :: s
    type(TComponent), pointer :: pc
    integer                   :: i

    ! Adjust number of cells to cube of integer
    if( this%NPart < NPartInCell ) this%NPart = NPartInCell
    this%NCells = ceiling( (real( this%NPart, RK ) &
&                             / real( NPartInCell, RK ))**Third )

    ! Set maximum number of particles
    if( EnsembleType .eq. EnsembleTypeGE .or. &
&       EnsembleType .eq. EnsembleTypeHA ) then
      this%NPartMax = 2 * this%NPart
    else
      this%NPartMax = this%NPart
    end if

    ! Normalize molar fractions
    s = 0
    do i = 1, this%NRealComponents
      pc => this%Component(i)
      s = s + pc%Fraction
    end do
    do i = 1, this%NRealComponents
      pc => this%Component(i)
      pc%Fraction = pc%Fraction / s
    end do

    ! Calculate number of particles in each component
    ! according to molar fraction
    this%Component(this%NRealComponents)%NPart = this%NPart
    do i = 1, this%NRealComponents - 1
      pc => this%Component(i)
      pc%NPart = nint( this%NPart * pc%Fraction )
      this%Component(this%NRealComponents)%NPart = &
&       this%Component(this%NRealComponents)%NPart - pc%NPart
    end do

    ! Set molar fractions according to real number of particles
    ! and calculate number of degrees of freedom
    this%NDFTran = 0
    this%NDFRot = 0
    do i = 1, this%NComponents
      pc => this%Component(i)
      pc%Fraction = real( pc%NPart, RK ) / real( this%NPart, RK )
      pc%NDFTran = pc%NPart * 3
      pc%NDFRot = pc%NPart * pc%Molecule%NDFRot
      pc%NDF = pc%NDFTran + pc%NDFRot
      this%NDFTran = this%NDFTran + pc%NDFTran
      this%NDFRot = this%NDFRot + pc%NDFRot
    end do
    this%NDF = this%NDFTran + this%NDFRot

    ! Calculate number of particles and test particles in each process
    ! and maximum number of test particles and fluctuating states
    this%NTestMax = 0
    this%NFluctMax = 0
    do i = 1, this%NComponents
      pc => this%Component(i)
      pc%NPart1 = 1 + (pc%NPart - 1) / NProcs
      pc%NPart0 = 1 + pc%NPart1 * NProc
      pc%NPart2 = min( pc%NPart0 + pc%NPart1 - 1, pc%NPart )
      pc%NPart1 = pc%NPart2 - pc%NPart0 + 1
      if( pc%NTest > 0 ) pc%NTest = 1 + (pc%NTest - 1) / NProcs
      pc%NTestAll = NProcs * pc%NTest
      this%NTestMax = max( pc%NTest, this%NTestMax )
      this%NFluctMax = max( pc%NFluctMax, this%NFluctMax )
    end do

    ! Calculate volume of simulation box
    this%Volume0 = this%NPart / this%RefDensity

    ! Update log file
    write( IOBuffer, '("Number of particles:", I11)' ) this%NPart
    call LogWrite
    do i = 1, this%NRealComponents
      pc => this%Component(i)
      write( IOBuffer, '("Molar fraction of ", A, ":", F6.3, &
&         ";  number of particles:", I11)' ) &
&       trim( pc%PotModFileName ), pc%Fraction, pc%NPart
      call LogWrite
      if( pc%NTestAll > 0 ) then
        write( IOBuffer, '("Number of test particles:", I11)' ) pc%NTestAll
        call LogWrite
      end if
    end do

  end subroutine TEnsemble_CalculateNPart



!==============================================================!
!  Subroutine TEnsemble_UpdateBoxLength                        !
!==============================================================!

  subroutine TEnsemble_UpdateBoxLength( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i, j

    ! Update size of the simulation box
    if( SimulationType .eq. SecondVirialCoeff ) then
      this%BoxLength = 100._RK
    else
      this%Density = this%NPart / this%Volume0
#if ARCH == 3
      this%BoxLength = cbrt( this%Volume0 )
#else
      this%BoxLength = this%Volume0**Third
#endif
      if (LongRange .eq. Ewald) then
        this%Kappa = this%KappaL / this%BoxLength
        do i=1,this%NComponents
          do j=1,this%NComponents
            this%Interaction(i,j)%Kappa = this%Kappa
          end do
        end do
        call EwaldSelfTerm(this)
      else if (LongRange .eq. PME ) then
        this%Kappa = this%KappaL / this%BoxLength
        do i=1,this%NComponents
          do j=1,this%NComponents
            this%Interaction(i,j)%Kappa = this%Kappa
          end do
        end do
        call PMESelfTermMC ( this )
      end if
    end if

    do i = 1, this%NComponents
      do j = 1, this%NComponents
        call UpdateBoxLength( this%Interaction( j, i ), this%BoxLength )
      end do
    end do

! Debug
!    write(*,*) 'this%BoxLength ist manuell festgelegt in ms2_ensemble.F90:1685 (UpdateBoxLength!)'
!    this%BoxLength = 3.644565432
! Debug
  end subroutine TEnsemble_UpdateBoxLength



!==============================================================!
!  Subroutine TEnsemble_UpdateFractions                        !
!==============================================================!

  subroutine TEnsemble_UpdateFractions( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer                   :: i
    type(TComponent), pointer :: pc

    ! Set molar fractions according to real number of particles
    ! and calculate number of degrees of freedom
    this%NDFTran = 0
    this%NDFRot = 0
    do i = 1, this%NComponents
      pc => this%Component(i)
      this%Component(i)%Fraction = real( pc%NPart, RK ) / &
&       real( this%NPart, RK )
      pc%NDFTran = pc%NPart * 3
      pc%NDFRot = pc%NPart * pc%Molecule%NDFRot
      pc%NDF = pc%NDFTran + pc%NDFRot
      this%NDFTran = this%NDFTran + pc%NDFTran
      this%NDFRot = this%NDFRot + pc%NDFRot
    end do
    this%NDF = this%NDFTran + this%NDFRot

  end subroutine TEnsemble_UpdateFractions



!==============================================================!
!  Subroutine TEnsemble_FindNSiteMax                           !
!==============================================================!

  subroutine TEnsemble_FindNSiteMax( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer                   :: i
    type(TComponent), pointer :: pc

    ! Calculate maximum numbers of sites in components
    this%NLJ126Max      = 0
    this%NChargeMax     = 0
    this%NDipoleMax     = 0
    this%NQuadrupoleMax = 0
    do i = 1, this%NComponents
      pc => this%Component(i)
      if( pc%Molecule%NLJ126 > this%NLJ126Max ) &
&       this%NLJ126Max = pc%Molecule%NLJ126
      if( pc%Molecule%NCharge > this%NChargeMax ) &
&       this%NChargeMax = pc%Molecule%NCharge
      if( pc%Molecule%NDipole > this%NDipoleMax ) &
&       this%NDipoleMax = pc%Molecule%NDipole
      if( pc%Molecule%NQuadrupole > this%NQuadrupoleMax ) &
&       this%NQuadrupoleMax = pc%Molecule%NQuadrupole
    end do

  end subroutine TEnsemble_FindNSiteMax



!==============================================================!
!  Subroutine TEnsemble_Allocate                               !
!==============================================================!

  subroutine TEnsemble_Allocate( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i
    integer :: stat
    integer :: number


    ! Nullify pointers
    nullify( this%P0Test )
    nullify( this%Q0Test )
    nullify( this%EPotTest )

    ! Allocate scale coefficients for sigma and epsilon
    allocate( this%ScaleSigma(this%NComponents, this%NComponents), &
&     STAT = stat )
    call AllocationError( stat, 'components', this%NComponents )
    allocate( this%ScaleEpsilon(this%NComponents, this%NComponents), &
&     STAT = stat )
    call AllocationError( stat, 'components', this%NComponents )

    ! Allocate test particles
    if( this%NTestMax > 0 ) then
      allocate( this%P0Test( this%NTestMax, 3 ), STAT = stat )
      call AllocationError( stat, 'test particles', this%NTestMax )
      allocate( this%Q0Test( this%NTestMax, 4 ), STAT = stat )
      call AllocationError( stat, 'test particles', this%NTestMax )
      allocate( this%EPotTest( this%NTestMax ), STAT = stat )
      call AllocationError( stat, 'test particles', this%NTestMax )
      write( IOBuffer, '("Memory for test particles allocated successfully")' )
      call LogWrite
    end if

    ! Allocate components
    if( EnsembleType .eq. EnsembleTypeGE .or. &
&       EnsembleType .eq. EnsembleTypeHA ) then
       do i = 1, this%NComponents
         this%Component(i)%NPartMax => this%NPartMax
         if( this%Component(i)%NTest > 0 ) then
           this%Component(i)%P0Test => this%P0Test
           this%Component(i)%Q0Test => this%Q0Test
         end if
         this%Component(i)%BoxLength => this%BoxLength
         call Allocate( this%Component(i) )
       end do
    else
       do i = 1, this%NComponents
         this%Component(i)%NPartMax => this%NPartMax
         if( this%Component(i)%NTest > 0 ) then
           this%Component(i)%P0Test => this%P0Test
           this%Component(i)%Q0Test => this%Q0Test
         end if
         this%Component(i)%BoxLength => this%BoxLength
         call Allocate( this%Component(i) )
       end do
    end if


    ! Ewald Summation - we allocate  - look above in the variable declaration
    if (LongRange .eq. Ewald) then
     allocate(this%Ewald_Prefac(this%NVecMax),STAT=stat)
     allocate(this%Ewald_Vec(3,this%NVecMax),STAT=stat)
     number = 0
     do i=1,this%NComponents
       number = number + this%Component(i)%NPart * this%Component(i)%Molecule%NCharge
     end do
!      allocate(this%SSin_fac(number),STAT=stat)
!      allocate(this%SCos_fac(number),STAT=stat)
    end if 

#ifdef ABL
     allocate(this%AblPS(this%NComponents,10),STAT=stat)
     allocate(this%AblPE(this%NComponents,10),STAT=stat)
     nullify(this%AblPS);
     nullify(this%AblPE);
     allocate(this%AblRhoS(this%NComponents,10),STAT=stat)
     allocate(this%AblRhoE(this%NComponents,10),STAT=stat)
     nullify(this%AblRhoS);
     nullify(this%AblRhoE);
#endif
  end subroutine TEnsemble_Allocate



!==============================================================!
!  Subroutine TEnsemble_DeallocateEPot                         !
!==============================================================!

  subroutine TEnsemble_DeallocateEPot( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i, j

    ! Deallocate potential energy matrix
    do i = 1, this%NComponents
      do j = 1, this%NComponents
        call DeallocateEPot( this%Interaction(i, j) )
      end do
    end do

  end subroutine TEnsemble_DeallocateEPot



!==============================================================!
!  Subroutine TEnsemble_Deallocate                             !
!==============================================================!

  subroutine TEnsemble_Deallocate( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Deallocate arrays
    do i = 1, this%NComponents
      call Deallocate( this%Component(i) )
    end do
    if( associated( this%P0Test ) ) then
      deallocate( this%P0Test )
    end if
    if( associated( this%Q0Test ) ) then
      deallocate( this%Q0Test )
    end if
    if( associated( this%EPotTest ) ) then
      deallocate( this%EPotTest )
    end if
    if( associated( this%ScaleEpsilon ) ) then
      deallocate( this%ScaleEpsilon )
    end if
    if( associated( this%ScaleSigma ) ) then
      deallocate( this%ScaleSigma )
    end if

  end subroutine TEnsemble_Deallocate



!==============================================================!
!  Subroutine TEnsemble_CalculateCorr                          !
!==============================================================!

  subroutine TEnsemble_CalculateCorr( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK)                        :: NPartInv, Scale, RFConst
    type(TComponent), pointer       :: pc
    type(TPotLJ126LJ126), pointer   :: plj
    integer                         :: i1, i2, j1, j2

    ! Assign local variables
    NPartInv = 1._RK / this%NPart
    if (LongRange .eq. ExtRField) then
      RFConst = -1._RK / this%RCutoffDipoleDipole**3 &
&       * ((this%RFEpsilon - 1._RK)*(1._RK+this%DebyeLen*this%RCutoffDipoleDipole)+&
&          this%RFEpsilon*(this%DebyeLen*this%RCutoffDipoleDipole)**2)&
&       / ((2._RK * this%RFEpsilon+1._RK)*(1._RK+this%DebyeLen*this%RCutoffDipoleDipole)+&
&          this%RFEpsilon*(this%DebyeLen*this%RCutoffDipoleDipole)**2)
    else 
      RFConst = -1._RK / this%RCutoffDipoleDipole**3 &
&       * (this%RFEpsilon - 1._RK) / (2._RK * this%RFEpsilon + 1._RK)
    endif

    ! Set maximum cutoff radius
    this%NRCutoffMax = 0
    this%RCutoffMax2 = 0._RK

    ! Zero long-range corrections
    this%EPotCorrLJ   = 0._RK
    this%EPotCorrRF   = 0._RK
    this%VirialCorrLJ = 0._RK
    do i1 = 1, this%NComponents
      this%Component(i1)%EPotTestCorrLJ = 0._8
    end do

    ! Calculate Lennard-Jones long-range corrections
    if( this%NLJ126Max > 0 ) then
      do i1 = 1, this%NComponents
        do i2 = 1, this%NComponents
          Scale = this%Component(i1)%NPart * this%Component(i2)%NPart &
&                 * NPartInv
          do j1 = 1, this%Component(i1)%Molecule%NLJ126
            do j2 = 1, this%Component(i2)%Molecule%NLJ126
              plj => this%Interaction(i1, i2)%PotLJ126LJ126(j1, j2)
              this%EPotCorrLJ = this%EPotCorrLJ + Scale * plj%EPotCorr
              this%VirialCorrLJ = this%VirialCorrLJ + Scale * plj%VirialCorr
              this%Component(i1)%EPotTestCorrLJ &
&               = this%Component(i1)%EPotTestCorrLJ &
&                 + this%Component(i2)%Fraction * plj%EPotTestCorr
              this%RCutoffMax2 = &
&               max( this%RCutoffMax2, 2._RK * sqrt( plj%RCutoffSquared ) )
            end do
          end do
        end do
      end do
      this%EPotCorrLJ = this%EPotCorrLJ / NProcs
      this%VirialCorrLJ = this%VirialCorrLJ / NProcs
    end if

    ! Calculate electrostatic long-range corrections
    ! This is the self term of the reaction field
    if(((this%NChargeMax>0).AND.((LongRange.eq.RField).or.(LongRange.eq.ExtRField))) &
&              .or.(this%NDipoleMax > 0) ) then
      do i1 = 1, this%NComponents
        pc => this%Component(i1)
        pc%EPotTestCorrRF = pc%Molecule%MueSquared * 2._RK * RFConst
        this%EPotCorrRF = this%EPotCorrRF + pc%Molecule%MueSquared * pc%NPart
      end do
      this%EPotCorrRF = this%EPotCorrRF * RFConst / NProcs
    end if

    if ((this%NChargeMax > 0) .and. ((LongRange .eq. Ewald) .or. (LongRange .eq. PME))) then
      this%EPotCorrRF = 0.0
      do i1 = 1, this%NComponents
        pc => this%Component(i1)
        pc%EPotTestCorrRF = 0._RK
      end do

    end if
  end subroutine TEnsemble_CalculateCorr



!==============================================================!
!  Subroutine TEnsemble_InitPositions                          !
!==============================================================!

  subroutine TEnsemble_InitPositions( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer                   :: i, j, k, l, n, nc, nm
    real(RK)                  :: xl
    integer                   :: comp(this%NComponents)
    type(TComponent), pointer :: pc

    ! Initialize comp array
    do i = 1, this%NComponents
      comp(i) = this%Component(i)%NPart
    end do

    ! Create FCC lattice
    n = 0
    xl = 1._RK / real( this%NCells, RK )
loop:do l = 1, NPartInCell
      do i = 1, this%NCells
        do j = 1, this%NCells
          do k = 1, this%NCells
            nc = select_component( comp )
            nm = comp(nc) + 1
            pc => this%Component(nc)
            pc%P0(nm, 1) = xl * (CellX(l) + i - 1) - 0.5_RK
            pc%P0(nm, 2) = xl * (CellY(l) + j - 1) - 0.5_RK
            pc%P0(nm, 3) = xl * (CellZ(l) + k - 1) - 0.5_RK
            n = n + 1
            if( n == this%NPart ) exit loop
          end do
        end do
      end do
    end do loop

    ! Save old positions
    do i = 1, this%NComponents
      this%Component(i)%P0old = this%Component(i)%P0
    end do


  contains



    ! Select component randomly
    function select_component(comp) result(nc)

      ! Declare arguments
      integer :: comp(:)

      ! Declare result
      integer :: nc

      ! Declare local variables
      real(RK) :: harvest
      real(RK) :: level(this%NComponents)
      integer  :: i

      ! Initialize level array
      level(1) = comp(1)
      do i = 2, this%NComponents
        level(i) = level(i - 1) + comp(i)
      end do
      level = level / sum(comp)

      ! Generate random number
      harvest = rnd( 0._RK, 1._RK )

      ! Select number of component
      nc = 1
      do while( harvest > level(nc) )
        nc = nc + 1
      end do
      comp(nc) = comp(nc) - 1

    end function select_component



  end subroutine TEnsemble_InitPositions



!==============================================================!
!  Subroutine TEnsemble_InitOrientations                       !
!==============================================================!

  subroutine TEnsemble_InitOrientations( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer                   :: i, j, k
    real(RK)                  :: r
    type(TComponent), pointer :: pc

    ! Set random orientations of particles
    do i = 1, this%NComponents
      pc => this%Component(i)
      if(pc%Molecule%isElongated ) then
        do j = 1, pc%NPart
          do
            do k = 1, 4
              pc%Q0(j, k) = rnd( -1._RK, 1._RK )
            end do
            r = sum( pc%Q0(j, :)**2 )
            if( r <= 1._RK ) exit
          end do
          pc%Q0(j, :) = pc%Q0(j, :) / sqrt( r )
        end do
      end if
    end do

  end subroutine TEnsemble_InitOrientations



!==============================================================!
!  Subroutine TEnsemble_InitMolecularDynamics                  !
!==============================================================!

  subroutine TEnsemble_InitMolecularDynamics( this, dealloc )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    logical, intent(in) :: dealloc

    ! Declare local variables
    integer :: i

    ! Reallocate ensemble
    if( dealloc ) call DeallocateEPot( this )

    ! Set initial velocities of particles
    call InitVelocities( this )

    ! Initialize integrator
    call InitIntegrator( this )

    ! Set old positions for displacement
    do i = 1, this%NComponents
      this%Component(i)%P0old = this%Component(i)%P0
    end do

    ! Remove net momentum
    call RemoveNetMomentum( this )

    ! Rescale velocities
    call CalculateEKin( this, .true. )
    call CalculateEKin( this, .false. )

    ! Set velocity scaling factor
    this%scale = 1._RK

  end subroutine TEnsemble_InitMolecularDynamics



!==============================================================!
!  Subroutine TEnsemble_InitVelocities                         !
!==============================================================!

  subroutine TEnsemble_InitVelocities( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Init velocities of each component
    do i = 1, this%NComponents
      call InitVelocities( this%Component(i) )
    end do

  end subroutine TEnsemble_InitVelocities



!==============================================================!
!  Subroutine TEnsemble_InitIntegrator                         !
!==============================================================!

  subroutine TEnsemble_InitIntegrator( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Check for root process
    if( .not. RootProc ) return

    ! Call InitIntegrator
    select case( IntegratorType )
    case( IntegratorTypeGear )
      call InitIntegratorGear( this )
    case( IntegratorTypeLeapFrog )
      call InitIntegratorLeapFrog( this )
    case( IntegratorTypeVerlet )
      call InitIntegratorVerlet( this )
    case( IntegratorTypeVV )
      call InitIntegratorVV( this )
    end select

  end subroutine TEnsemble_InitIntegrator



!==============================================================!
!  Subroutine TEnsemble_InitIntegratorGear                     !
!==============================================================!

  subroutine TEnsemble_InitIntegratorGear( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Call InitIntegrator for each component
    do i = 1, this%NComponents
      call InitIntegratorGear( this%Component(i) )
    end do

    ! Zero volume time derivatives
    this%Volume1 = 0._RK
    this%Volume2 = 0._RK
    this%Volume3 = 0._RK
    this%Volume4 = 0._RK
    this%Volume5 = 0._RK

  end subroutine TEnsemble_InitIntegratorGear



!==============================================================!
!  Subroutine TEnsemble_InitIntegratorLeap                     !
!==============================================================!

  subroutine TEnsemble_InitIntegratorLeap( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Call InitIntegrator for each component
    do i = 1, this%NComponents
      call InitIntegratorLeapFrog( this%Component(i) )
    end do

    ! Zero volume time derivatives
    this%Volume1 = 0._RK
    this%Volume2 = 0._RK

  end subroutine TEnsemble_InitIntegratorLeap



!==============================================================!
!  Subroutine TEnsemble_InitIntegratorVerlet                   !
!==============================================================!

  subroutine TEnsemble_InitIntegratorVerlet( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Call InitIntegrator for each component
    do i = 1, this%NComponents
      call InitIntegratorVerlet( this%Component(i) )
    end do

  end subroutine TEnsemble_InitIntegratorVerlet



!==============================================================!
!  Subroutine TEnsemble_InitIntegratorVV                       !
!==============================================================!

  subroutine TEnsemble_InitIntegratorVV( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Call InitIntegrator for each component
    do i = 1, this%NComponents
      call InitIntegratorVV( this%Component(i) )
    end do

  end subroutine TEnsemble_InitIntegratorVV



!==============================================================!
!  Subroutine TEnsemble_RemoveNetMomentum                      !
!==============================================================!

  subroutine TEnsemble_RemoveNetMomentum( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Remove net momentum of each component
    do i = 1, this%NComponents
      call RemoveNetMomentum( this%Component(i) )
    end do

  end subroutine TEnsemble_RemoveNetMomentum



!==============================================================!
!  Subroutine TEnsemble_CalculateEKin                          !
!==============================================================!

  subroutine TEnsemble_CalculateEKin( this, rescale )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this
    logical         :: rescale

    ! Declare local variables
    integer                   :: i
    integer                   :: np
    real(RK)                  :: scale
    type(TComponent), pointer :: pc

    ! Check for root process
    if( RootProc ) then

      ! Nullify kinetic energies
      this%EKinTran = 0._RK
      this%EKinRot = 0._RK

      ! Loop over components
      do i = 1, this%NComponents
        pc => this%Component(i)
        call CalculateEKin( pc )
        this%EKinTran = this%EKinTran + pc%EKinTran
        this%EKinRot = this%EKinRot + pc%EKinRot
      end do
      this%EKin = this%EKinTran + this%EKinRot

      ! Calculate temperature
      this%Temperature = 2._RK * this%EKin / this%NDF

      ! Rescale velocities
      if( rescale ) then
      scale = sqrt( this%RefTemperature / this%Temperature )
        do i = 1, this%NComponents
          pc => this%Component(i)
          np = pc%NPart
          pc%P1(:, :) = pc%P1(:, :) * scale
          if( pc%Molecule%isElongated ) then
            pc%W0(:, :) = pc%W0(:, :) * scale
          end if
        end do
      else
        scale = 1._RK
      end if

      ! Save scaling factor
      this%scale = scale

    end if

    ! Broadcast temperature
#if MPI_VER > 0
    call MPI_Bcast( this%Temperature, 1, MPI_DOUBLE_PRECISION, &
&     NRootProc, MPI_COMM_WORLD, ierror )
#endif

  end subroutine TEnsemble_CalculateEKin



!==============================================================!
!  Subroutine TEnsemble_CheckNPart                             !
!==============================================================!

  subroutine TEnsemble_CheckNPart( this, NPartOk )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)        :: this
    logical, intent(inout) :: NPartOk

    ! Test if number of particles is ok
    NPartOk = NPartOk .and. this%NPart > this%NPartLBound &
&                     .and. this%NPart < this%NPartUBound

  end subroutine TEnsemble_CheckNPart



!==============================================================!
!  Subroutine TEnsemble_ResetEnsemble                          !
!==============================================================!

  subroutine TEnsemble_ResetEnsemble( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
!     integer :: i

    ! Calculate new initial density
    this%RefDensity = this%RefDensity * real( this%NPart, RK ) &
&     / real( this%NPartInitial, RK )
    this%NPart = this%NPartInitial

    ! Reset fractions
!     do i = 1, this%NComponents
!       this%Component(i)%Fraction = this%Component(i)%StartFraction
!     end do

    ! Calculate number of particles of each component
    call CalculateNPart( this )

    ! Calculate long-range corrections
    call CalculateCorr( this )

    ! Update all BoxLength-dependent constants
    call UpdateBoxLength( this )

    ! Set initial positions of particles in simulation box
    call InitPositions( this )

    ! Set initial orientations of particles in simulation box
    call InitOrientations( this )

    ! Convert molecular coordinates to atom positions
    call Mol2Atom( this )

    ! Set all potential energy matrices
    call Energy( this, this%EPot )
    call UpdateEnergy( this )

  end subroutine TEnsemble_ResetEnsemble



!==============================================================!
!  Subroutine TEnsemble_RunMDStep                              !
!==============================================================!

  subroutine TEnsemble_RunMDStep( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Zero displacement
    if( Step == 1 ) then
      do i = 1, this%NComponents
        this%Component(i)%Disp(:, :) = 0._RK
      end do
    end if

    ! Run MD simulation step
    call Predict( this )
    call Mol2Atom( this )

#if DEBUG_4Part > 0
! DEBUG
! DEBUG STEPHAN
    this%Component(1)%Molecule%SiteCharge(1)%RX(1) = 0.5041221762148172/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RY(1) = 0.8107396528595158/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(1) = 0.1894949962472363/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RX(1) = 0.5022685541517008/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RY(1) = 0.8117808162647254/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(1) = 0.1869276678455920/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RX(2) = 2.099380236705095 /this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RY(2) = 1.860649777419017/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(2) = 1.313235759816588/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RX(2) = 2.097708401492369/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RY(2) = 1.858592437418999/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(2) = 1.311215009892196/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RX(3) = 1.317611713633559 /this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RY(3) = 2.087579656401032/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(3) = 0.4484239832120850/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RX(3) = 1.315500144816407/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RY(3) = 2.085384351686846/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(3) = 0.4470700905082808/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RX(4) = 1.482005925799646/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RY(4) = 0.7714424297398743/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(4) = 1.556938840875683/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RX(4) = 1.484582190214423/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RY(4) = 0.7716945714654903/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(4) = 1.554838746276274/this%Boxlength

    this%Component(1)%Molecule%SiteCharge(1)%RX(1) = 0.5041221762148172 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(1) = 0.8107396528595158 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(1) = 0.1894949962472363 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(1) = 0.5022685541517008 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(1) = 0.8117808162647254 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(1) = 0.1869276678455920 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(2) = 2.099380236705095  / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(2) = 1.860649777419017 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(2) = 1.313235759816588 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(2) = 2.097708401492369 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(2) = 1.858592437418999 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(2) = 1.311215009892196 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(3) = 1.317611713633559 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(3) = 2.087579656401032 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(3) = 0.448423983212085 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(3) = 1.315500144816407 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(3) = 2.085384351686846 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(3) = 0.447070090508280 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(4) = 1.482005925799646 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(4) = 0.771442429739874 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(4) = 1.556938840875683 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(4) = 1.484582190214423 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(4) = 0.771694571465490 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(4) = 1.554838746276274 / this%BoxLength
#endif

    call Force( this )
    call ChemicalPotential( this )
    call Atom2Mol( this )
    call Correct( this )
#if CONSTR > 0
    call Constraints(this)
#endif

    call CalculateEKin( this, ConstantTemperature .or. Equilibration )
    if( .not. Equilibration .and. this%RCutoffMax2 > this%BoxLength ) &
&     this%NRCutoffMax = this%NRCutoffMax + 1

  end subroutine TEnsemble_RunMDStep



!==============================================================!
!  Subroutine TEnsemble_RunMCStep                              !
!==============================================================!

  subroutine TEnsemble_RunMCStep( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer  :: r, s
    integer  :: nc, np, ndf
    integer  :: i
    real(RK) :: rx, sx
    real(RK) :: diffpressure
    real(RK) :: NDFsystem

    ! Zero number of MC attempts and successes
    if( Step == 1 ) call ZeroNAttempts( this )

#if DEBUG_4Part > 0
! DEBUG STEPHAN
    this%Component(1)%Molecule%SiteCharge(1)%RX(1) = 0.5041221762148172/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RY(1) = 0.8107396528595158/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(1) = 0.1894949962472363/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RX(1) = 0.5022685541517008/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RY(1) = 0.8117808162647254/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(1) = 0.1869276678455920/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RX(2) = 2.099380236705095 /this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RY(2) = 1.860649777419017/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(2) = 1.313235759816588/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RX(2) = 2.097708401492369/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RY(2) = 1.858592437418999/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(2) = 1.311215009892196/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RX(3) = 1.317611713633559 /this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RY(3) = 2.087579656401032/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(3) = 0.4484239832120850/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RX(3) = 1.315500144816407/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RY(3) = 2.085384351686846/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(3) = 0.4470700905082808/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RX(4) = 1.482005925799646/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RY(4) = 0.7714424297398743/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(4) = 1.556938840875683/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RX(4) = 1.484582190214423/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RY(4) = 0.7716945714654903/this%Boxlength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(4) = 1.554838746276274/this%Boxlength

    this%Component(1)%Molecule%SiteCharge(1)%RX(1) = 0.5041221762148172 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(1) = 0.8107396528595158 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(1) = 0.1894949962472363 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(1) = 0.5022685541517008 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(1) = 0.8117808162647254 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(1) = 0.1869276678455920 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(2) = 2.099380236705095  / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(2) = 1.860649777419017 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(2) = 1.313235759816588 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(2) = 2.097708401492369 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(2) = 1.858592437418999 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(2) = 1.311215009892196 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(3) = 1.317611713633559 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(3) = 2.087579656401032 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(3) = 0.448423983212085 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(3) = 1.315500144816407 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(3) = 2.085384351686846 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(3) = 0.447070090508280 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(4) = 1.482005925799646 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(4) = 0.771442429739874 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(4) = 1.556938840875683 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(4) = 1.484582190214423 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(4) = 0.771694571465490 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(4) = 1.554838746276274 / this%BoxLength
#endif

    ! Outer loop
    NDFsystem = this%NDF
    do i = 1, NDFsystem / 3

      ! Choose particle randomly
      s = 0
      r = rnd( this%NDF )
loop1:do nc = 1, this%NComponents
        s = s + this%Component(nc)%NDF
        if( r <= s ) exit loop1
      end do loop1
      ndf = this%Component(nc)%Molecule%NDF
      np = 1 + (s - r) / ndf

      ! Move or rotate
      if( mod( s - r, ndf ) < 3 ) then
        call Move( this, nc, np )
      else
        call Rotate( this, nc, np )
      end if

    end do

    ! Calculate potential energy and virial
#if MPI_VER > 0
    call MPI_Allreduce( GetEnergy( this ), this%EPot, 1 , &
&     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
    call MPI_Allreduce( GetVirial( this ), this%Virial, 1 , &
&     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
#else
    this%EPot = GetEnergy( this )
    this%Virial = GetVirial( this )
#endif

    ! Resize simulation box
    if( ConstantPressure .and. .not. NVTEquilibration ) then
      call Resize( this )

      ! Check whether cutoff radius is too large
      if( this%RCutoffMax2 > this%BoxLength ) &
&       this%NRCutoffMax = this%NRCutoffMax + 1
    end if

    ! Calculate pressure
    this%Pressure = this%Density * this%Temperature &
&                   + this%Virial / this%Volume0

    if( EnsembleType .eq. EnsembleTypeGE ) then

      if( .not. NVTEquilibration ) then
        ! Update chemical potentials
        diffpressure = ( this%Pressure - this%RefPressure ) / this%Temperature
        do nc = 1, this%NComponents
          call UpdateChemPot( this%Component(nc), diffpressure )
        end do

        ! Attempt inserts and deletes
        do i = 1, 2
          sx = 0._RK
          rx = rnd( 0._RK, 1._RK )
loop2:    do nc = 1, this%NComponents
            sx = sx + this%Component(nc)%Fraction
            if( rx <= sx ) exit loop2
          end do loop2
          call Insert( this, nc )

          s = 0._RK
          r = rnd( this%NPart )
loop3:    do nc = 1, this%NComponents
            s = s + this%Component(nc)%NPart
            if( r <= s ) exit loop3
          end do loop3
          np = 1 + s - r
          call Delete( this, nc, np )

        end do
      end if

    elseif( EnsembleType .eq. EnsembleTypeHA ) then

      if( .not. NVTEquilibration ) then
        ! Attempt inserts and deletes on phase changing component (first one)
        do i = 1, 2
          if( rnd( 0._RK, 1._RK ) < this%Component(1)%Fraction ) &
&           call Insert( this, 1 )

          np = rnd( this%NPart )
          if( np <= this%Component(1)%NPart ) &
&           call Delete( this, 1, np )
        end do
      end if

    else

      ! Calculate chemical potential
      call ChemicalPotential( this )

    end if

    ! Update MC displacements
    if( Equilibration .and. mod( Step, DispUpdateFrequency ) == 0 ) &
&     call UpdateDisplacements( this )

  end subroutine TEnsemble_RunMCStep



!==============================================================!
!  Subroutine TEnsemble_RunSVCStep                             !
!==============================================================!

  subroutine TEnsemble_RunSVCStep( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK)                    :: r, rdist, rsquared(NSteps)
    real(RK)                    :: betaneg, betaneg1, betaneg2, Bij0
    integer                     :: i, j, n, np
    type(TInteraction), pointer :: pi

    ! Inverse temperature
    betaneg = -1._RK / this%Temperature
    betaneg1 = -1._RK / ( .9999_RK * this%Temperature )
    betaneg2 = -1._RK / ( 1.0001_RK * this%Temperature )

    ! Constant for integration
    Bij0 = Pi23 * MinRadius**3

    ! Calculate distance
    rdist = (MaxRadius - MinRadius) / real(NSteps - 1, RK)
    r = MinRadius + (Step - 1) * rdist
    do i = 1, Step
      rsquared(i) = ( MinRadius + (i - 1) * rdist ) ** 2
    end do

    ! Set distance
    do i = 2, this%NComponents, 2
      this%Component(i)%P0(:, 1) = r / this%BoxLength
      call Mol2Atom( this%Component(i), this%Component(i)%NPart )
    end do

    ! Loop over components
    do i = 1, this%NComponents, 2
      do j = i + 1, this%NComponents, 2
        pi => this%Interaction(i, j)
        n = pi%NPart2
        pi%MayerFFunction(Step) = 0._RK

        ! Loop over particles
        do np = 1, this%Component(i)%NPart
          call Energy( pi, np, this%BoxLength )

          ! Sum Mayer f-function
          pi%MayerFFunction(Step) = pi%MayerFFunction(Step) + &
&           sum( exp( betaneg * pi%EPot1(1:n) ) - 1._RK )
          pi%MayerFFunction1(Step) = pi%MayerFFunction1(Step) + &
&           sum( exp( betaneg1 * pi%EPot1(1:n) ) - 1._RK )
          pi%MayerFFunction2(Step) = pi%MayerFFunction2(Step) + &
&           sum( exp( betaneg2 * pi%EPot1(1:n) ) - 1._RK )

        end do

        ! Average Mayer f-function
        pi%MayerFFunction(Step) = pi%MayerFFunction(Step) / &
&         real( n * this%Component(i)%NPart, RK )
        pi%MayerFFunction1(Step) = pi%MayerFFunction1(Step) / &
&         real( n * this%Component(i)%NPart, RK )
        pi%MayerFFunction2(Step) = pi%MayerFFunction2(Step) / &
&         real( n * this%Component(i)%NPart, RK )

        ! Integrate Mayer f-function
        pi%IntFFunction = Bij0 + Piminus2 * &
&         simpson( pi%MayerFFunction * rsquared, rdist, Step)
        pi%IntFFunction1 = Bij0 + Piminus2 * &
&         simpson( pi%MayerFFunction1 * rsquared, rdist, Step)
        pi%IntFFunction2 = Bij0 + Piminus2 * &
&         simpson( pi%MayerFFunction2 * rsquared, rdist, Step)

      end do
    end do

    ! Update result header
    if( Step == 1 ) then
      call FileWriteBlank( this%iounit_result )

      ! Number of steps
      write( IOBuffer, '("     NR")' )
      call FileWriteNoAdvance( this%iounit_result )

      ! Radius
      write( IOBuffer, '("         R")' )
      call FileWriteNoAdvance( this%iounit_result )

      ! Mayer f-function and integral
      do i = 1, this%NComponents, 2
        do j = i + 1, this%NComponents, 2
          write( IOBuffer, '("      F", I2, "-", I2)' ) i, j - 1
          call FileWriteNoAdvance( this%iounit_result )
          write( IOBuffer, '("    IntF", I2, "-", I2)' ) i, j - 1
          call FileWriteNoAdvance( this%iounit_result )
        end do
      end do
      call FileWriteBlank( this%iounit_result )
    end if

    ! Number of steps
    write( IOBuffer, '(I7)' ) Step
    call FileWriteNoAdvance( this%iounit_result )

    ! Radius
    write( IOBuffer, '(F10.5)' ) r
    call FileWriteNoAdvance( this%iounit_result )

    ! Mayer f-function and integral
    do i = 1, this%NComponents, 2
      do j = i + 1, this%NComponents, 2
        write( IOBuffer, '(F12.4)' ) this%Interaction(i, j)%MayerFFunction(Step)
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(F13.4)' ) this%Interaction(i, j)%IntFFunction(Step)
        call FileWriteNoAdvance( this%iounit_result )
      end do
    end do

    call FileWriteBlank( this%iounit_result )
#if ARCH == 2
    call flush( this%iounit_result )
#endif



  contains



    function simpson(values, stepsize, n) result(integral)

      ! Declare arguments
      real(RK), intent(in) :: values(:), stepsize
      integer, intent(in)  :: n

      ! Declare result
      real(RK) :: integral(n)

      ! Declare local variables
      integer :: i

      ! Return if no values to integrate
      if( n < 1 ) return

      ! Initialize result
      integral = 0._RK

      ! Calculate integral via Simpson's rule
      do i = 3, n, 2
        integral(i) = integral(i-2) + values(i) + 4._RK * values(i-1) &
&         + values(i-2)
        integral(i-1) = .5 * (integral(i) + integral(i-2))
      end do
      if( mod(n, 2) == 0 .and. n > 2 ) integral(n) = integral(n-1) &
&       + .5_RK * values(n) + 2._RK * values(n-1) + .5_RK * values(n-2)
      integral = integral * stepsize / 3._RK

    end function



  end subroutine TEnsemble_RunSVCStep



!==============================================================!
!  Subroutine TEnsemble_Mol2Atom                               !
!==============================================================!

  subroutine TEnsemble_Mol2Atom( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer      :: i

    ! Call Mol2Atom for each component
    do i = 1, this%NComponents
      call Mol2Atom( this%Component(i), this%Component(i)%NPart )
    end do

  end subroutine TEnsemble_Mol2Atom



!==============================================================!
!  Subroutine TEnsemble_Atom2Mol                               !
!==============================================================!

  subroutine TEnsemble_Atom2Mol( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Call Atom2Mol for each component
    do i = 1, this%NComponents
      call Atom2Mol( this%Component(i), this%Component(i)%NPart )
    end do

  end subroutine TEnsemble_Atom2Mol



!==============================================================!
!  Subroutine TEnsemble_Predict                                !
!==============================================================!

  subroutine TEnsemble_Predict( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Call predictor
    select case( IntegratorType )
    case( IntegratorTypeGear )
      call PredictGear( this )
    case( IntegratorTypeLeapFrog )
      call PredictLeapFrog( this )
    case( IntegratorTypeVerlet )
      call PredictVerlet( this )
    case( IntegratorTypeVV )
      call PredictVV( this )
    end select

  end subroutine TEnsemble_Predict



!==============================================================!
!  Subroutine TEnsemble_Correct                                !
!==============================================================!

  subroutine TEnsemble_Correct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Call corrector
    select case( IntegratorType )
    case( IntegratorTypeGear )
      call CorrectGear( this )
    case( IntegratorTypeLeapFrog )
      call CorrectLeapFrog( this )
    case( IntegratorTypeVerlet )
      call CorrectVerlet( this )
    case( IntegratorTypeVV )
      call CorrectVV( this )
    end select

  end subroutine TEnsemble_Correct



!==============================================================!
!  Subroutine TEnsemble_PredictGear                            !
!==============================================================!

  subroutine TEnsemble_PredictGear( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Call predictor for each component
    if( RootProc ) then
      do i = 1, this%NComponents
        call PredictGear( this%Component(i) )
      end do
    end if

    ! Predict volume of simulation box
    if( ConstantPressure .and. .not. NVTEquilibration ) then
      if( RootProc ) then
        this%Volume0 = this%Volume0 &
&                    + this%Volume1 &
&                    + this%Volume2 &
&                    + this%Volume3 &
&                    + this%Volume4 &
&                    + this%Volume5
        this%Volume1 = this%Volume1 &
&            + 2._RK * this%Volume2 &
&            + 3._RK * this%Volume3 &
&            + 4._RK * this%Volume4 &
&            + 5._RK * this%Volume5
        this%Volume2 = this%Volume2 &
&            + 3._RK * this%Volume3 &
&            + 6._RK * this%Volume4 &
&            +10._RK * this%Volume5
        this%Volume3 = this%Volume3 &
&            + 4._RK * this%Volume4 &
&            +10._RK * this%Volume5
        this%Volume4 = this%Volume4 &
&            + 5._RK * this%Volume5
      end if
#if MPI_VER > 0
      call MPI_Bcast( this%Volume0, 1, MPI_DOUBLE_PRECISION, NRootProc, &
&       MPI_COMM_WORLD, ierror )
#endif
      call UpdateBoxLength( this )
    end if

  end subroutine TEnsemble_PredictGear



!==============================================================!
!  Subroutine TEnsemble_CorrectGear                            !
!==============================================================!

  subroutine TEnsemble_CorrectGear( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer  :: i
    real(RK) :: dLogVolumeThird, Volume2, Corr

#ifdef ABL
    real(RK) :: vol
    real(RK) :: fac
    real(RK) :: denom,denom2
    real(RK) :: nen
    integer  :: j
#endif

    ! Call corrector for each component
    if( RootProc ) then
      dLogVolumeThird = this%Volume1 / (3._RK * this%Volume0)
      do i = 1, this%NComponents
        call CorrectGear( this%Component(i), dLogVolumeThird )
      end do
    end if

    ! Correct volume of simulation box
    if( ConstantPressure .and. .not. NVTEquilibration ) then
      if( RootProc ) then
        Volume2 = (this%Pressure - this%RefPressure) &
&         * TimeStepSquared2 / this%PistonMass
        Corr = Volume2 - this%Volume2
        this%Volume0 = this%Volume0 + Corr * Gear20
        this%Volume1 = this%Volume1 + Corr * Gear21
        this%Volume2 =      Volume2
        this%Volume3 = this%Volume3 + Corr * Gear23
        this%Volume4 = this%Volume4 + Corr * Gear24
        this%Volume5 = this%Volume5 + Corr * Gear25

#if ABL
        vol = this%Volume0 + this%Volume1 + this%Volume2 + this%Volume3 +&
&             this%Volume4 + this%Volume5
        fac = TimeStepSquared2*Gear20
        denom = fac*(this%Pressure - this%RefPressure) - &
&                       this%PistonMass*this%Volume2*Gear20
        denom2 = denom**2
        nen = this%PistonMass*fac / (vol * denom2)
        do i=1,this%NComponents
          do j=1,this%Component(i)%Molecule%NLJ126
            this%AblPS(i,j)   =  this%AblPS(i,j) + &
&              this%Interaction(1, 1)%PotLJ126LJ126(i, j)%AblSigCorr(i,j)
            this%AblPE(i,j)   =  this%AblPE(i,j) + &
&                this%Interaction(1, 1)%PotLJ126LJ126(i, j)%AblEpsCorr(i,j)
            this%AblRhoS(i,j) = nen * this%AblPS(i,j)
            this%AblRhoE(i,j) = nen * this%AblPE(i,j)
          end do
        end do
#endif

      end if
#if MPI_VER > 0
      call MPI_Bcast( this%Volume0, 1, MPI_DOUBLE_PRECISION, NRootProc, &
&       MPI_COMM_WORLD, ierror )
#endif
      call UpdateBoxLength( this )
    end if

  end subroutine TEnsemble_CorrectGear



!==============================================================!
!  Subroutine TEnsemble_PredictLeapFrog                        !
!==============================================================!

  subroutine TEnsemble_PredictLeapFrog( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Call predictor for each component
    if( RootProc ) then
      do i = 1, this%NComponents
        call PredictLeapFrog( this%Component(i), this%scale )
      end do
    end if

    ! Predict volume of simulation box
    if( ConstantPressure .and. .not. NVTEquilibration ) then
      if( RootProc ) then
        this%Volume1 = this%Volume1 + this%Volume2
        this%Volume0 = this%Volume0 + this%Volume1
      end if
#if MPI_VER > 0
      call MPI_Bcast( this%Volume0, 1, MPI_DOUBLE_PRECISION, NRootProc, &
&       MPI_COMM_WORLD, ierror )
#endif
      call UpdateBoxLength( this )
    end if

  end subroutine TEnsemble_PredictLeapFrog



!==============================================================!
!  Subroutine TEnsemble_CorrectLeapFrog                        !
!==============================================================!

  subroutine TEnsemble_CorrectLeapFrog( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer  :: i
    real(RK) :: dLogVolumeThird

    ! Call corrector for each component
    if( RootProc ) then
      dLogVolumeThird = this%Volume1 / (3._RK * this%Volume0)
      do i = 1, this%NComponents
        call CorrectLeapFrog( this%Component(i), dLogVolumeThird )
      end do
    end if

    ! Correct volume of simulation box
    if( ConstantPressure .and. .not. NVTEquilibration ) then
      if( RootProc ) then
        this%Volume2 = (this%Pressure - this%RefPressure) &
&         * TimeStepSquared2 / this%PistonMass
        this%Volume1 = this%Volume1 + this%Volume2
      end if
    end if

  end subroutine TEnsemble_CorrectLeapFrog



!==============================================================!
!  Subroutine TEnsemble_PredictVerlet                          !
!==============================================================!

  subroutine TEnsemble_PredictVerlet( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Call predictor for each component
    do i = 1, this%NComponents
      call PredictVerlet( this%Component(i) )
    end do

  end subroutine TEnsemble_PredictVerlet



!==============================================================!
!  Subroutine TEnsemble_CorrectVerlet                          !
!==============================================================!

  subroutine TEnsemble_CorrectVerlet( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer  :: i
    real(RK) :: dLogVolumeThird

    ! Check for root process
    if( .not. RootProc ) return

    ! Call corrector for each component
    dLogVolumeThird = this%Volume1 / (3._RK * this%Volume0)
    do i = 1, this%NComponents
      call CorrectVerlet( this%Component(i), dLogVolumeThird )
    end do

  end subroutine TEnsemble_CorrectVerlet



!==============================================================!
!  Subroutine TEnsemble_PredictVV                              !
!==============================================================!

  subroutine TEnsemble_PredictVV( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Call predictor for each component
    do i = 1, this%NComponents
      call PredictVV( this%Component(i) )
    end do

  end subroutine TEnsemble_PredictVV



!==============================================================!
!  Subroutine TEnsemble_CorrectVV                              !
!==============================================================!

  subroutine TEnsemble_CorrectVV( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer  :: i
    real(RK) :: dLogVolumeThird

    ! Check for root process
    if( .not. RootProc ) return

    ! Call corrector for each component
    dLogVolumeThird = this%Volume1 / (3._RK * this%Volume0)
    do i = 1, this%NComponents
      call CorrectVV( this%Component(i), dLogVolumeThird )
    end do

  end subroutine TEnsemble_CorrectVV



!==============================================================!
!  Subroutine TEnsemble_Force                                  !
!==============================================================!

  subroutine TEnsemble_Force( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK)                  :: EPot, Virial
    integer                   :: i, j
    type(TComponent), pointer :: pc
#ifdef ABL
    integer                   :: k,l
    integer                   :: numbi, numbj, numb
#endif

    ! Zero forces
    do i = 1, this%NComponents
      pc => this%Component(i)
      do j = 1, this%Component(i)%Molecule%NLJ126
        pc%Molecule%SiteLJ126(j)%FX(1:pc%NPart) = 0._RK
        pc%Molecule%SiteLJ126(j)%FY(1:pc%NPart) = 0._RK
        pc%Molecule%SiteLJ126(j)%FZ(1:pc%NPart) = 0._RK
      end do
      do j = 1, this%Component(i)%Molecule%NCharge
        pc%Molecule%SiteCharge(j)%FX(1:pc%NPart) = 0._RK
        pc%Molecule%SiteCharge(j)%FY(1:pc%NPart) = 0._RK
        pc%Molecule%SiteCharge(j)%FZ(1:pc%NPart) = 0._RK
      end do
      do j = 1, this%Component(i)%Molecule%NDipole
        pc%Molecule%SiteDipole(j)%FX(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%FY(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%FZ(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%TX(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%TY(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%TZ(1:pc%NPart) = 0._RK
      end do
      do j = 1, this%Component(i)%Molecule%NQuadrupole
        pc%Molecule%SiteQuadrupole(j)%FX(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%FY(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%FZ(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%TX(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%TY(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%TZ(1:pc%NPart) = 0._RK
      end do
      if( pc%Molecule%isElongated ) then
        pc%tRFX(:) = 0._RK
        pc%tRFY(:) = 0._RK
        pc%tRFZ(:) = 0._RK
      end if
    end do

    ! Zero potential
    EPot = this%Density * this%EPotCorrLJ + this%EPotCorrRF

    ! Zero virial
    Virial = this%Density * this%VirialCorrLJ

    ! Loop over components
    do i = 1, this%NComponents
      do j = i, this%NComponents
#ifndef ABL
        call Force( this%Interaction( i, j ), &
&                   EPot, Virial, this%BoxLength )
#else
        this%Interaction(i,j)%AblPS => this%AblPS
        this%Interaction(i,j)%AblPE => this%AblPE
        call Force( this%Interaction( i, j ), &
&                   EPot, Virial, this%BoxLength, i, j)
#endif
      end do
    end do

! #ifdef ABL
!     do i=1,this%NComponents
!       do j=i,this%NComponents
!         pli => this%Component(i)%Molecule
!         pli => this%Component(j)%Molecule
!         do k=1,pli%NSiteLJ126
!           do l=1,plj%NSiteLJ126
!             numbi = (k-1)*pli%NSiteLJ126
!             numbj = (l-1)*plj%NSiteLJ126
!             numb  = numbi + numbj + 1
!             if ( (i .eq. j) .AND. (k .eq.l) ) then
!               VirialAblSig = VirialAblSig + 2._RK*this%Interaction(i,j)%AblS(numb)
! !              VirialAblEps = VirialAblEps + this%Interaction(i,j)%AblE(i)
!             else
! !              epsi = this%Component(i)%Molecule%Site
!               VirialAblSig = VirialAblSig +       this%Interaction(i,j)%AblS(numb)
! !              VirialAblEps = VirialAblEps + 2._RK*this%Interaction(i,j)%AblE(i)
!             end if
!           end do
!         end do
!       end do
!     end do
! 
! 
!     posi = 0
!     do i=1,this%NComponents
!       pli => this%Component(i)%Molecule
!       do k=1,pli%NSiteLJ126
!       posi = posi + 1
!         do j=1,this%NComponents
!           VirialAblSig(posi) = VirialAblSig(posi) + this%Interaction(i,j)%AblS(numb)
!         pli => this%Component(j)%Molecule
!           do l=1,plj%NSiteLJ126
!             numbi = (k-1)*pli%NSiteLJ126
!             numbj = (l-1)*plj%NSiteLJ126
!             numb  = numbi + numbj + 1
!             if ( (i .eq. j) .AND. (k .eq.l) ) then
!               VirialAblSig = VirialAblSig + 2._RK*this%Interaction(i,j)%AblS(numb)
! !              VirialAblEps = VirialAblEps + this%Interaction(i,j)%AblE(i)
!             else
! !              epsi = this%Component(i)%Molecule%Site
!               VirialAblSig = VirialAblSig +       this%Interaction(i,j)%AblS(numb)
! !              VirialAblEps = VirialAblEps + 2._RK*this%Interaction(i,j)%AblE(i)
!             end if
!           end do
!         end do
!       end do
!     end do
! 
! 
! #endif


    if (LongRange .eq. Ewald) then
      call EwaldFourierTerm (this)
    end if
    if (LongRange .eq. PME) then
      call PMEFourierTerm (this)
    end if


    ! Collect sums from all processes
#if MPI_VER > 0
    call MPI_Reduce( EPot, this%EPot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
&     NRootProc, MPI_COMM_WORLD, ierror )
    call MPI_Reduce( Virial, this%Virial, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
&     NRootProc, MPI_COMM_WORLD, ierror )
#else
    this%EPot = EPot
    this%Virial = Virial
#endif

    if ((LongRange .eq. Ewald) .or. (LongRange .eq. PME))then
      if( RootProc ) then

       this%EPot   = this%EPot   + this%UFourier + this%USelbstTerm + this%UIntra
       this%Virial = this%Virial + this%EVirial
      end if
    end if

!     write(*,*) this%UFourier
    ! Calculate pressure
    this%Pressure = this%Density * this%Temperature + this%Virial / this%Volume0

  end subroutine TEnsemble_Force



!==============================================================!
!  Subroutine TEnsemble_ChemicalPotential                      !
!==============================================================!

  subroutine TEnsemble_ChemicalPotential( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK)                  :: ChemPot, qsum
    integer                   :: i, j
    integer                   :: ndf, ndfmove, ndfbiased, ndffluct, ndfchange, &
&                                ndfcp
    integer                   :: r, s, nc, np, ncf, npf
    type(TComponent), pointer :: pc
!DEBUG
!     type(TComponent), pointer :: pcf
    integer                   :: nstate( 0:this%NFluctMax ) !, counter
!DEBUG

    ! No calculation of chemical potential in equilibration
    if( Equilibration ) then
      do i = 1, this%NRealComponents
        this%Component(i)%CalcChemPot = .false.
        this%Component(i)%ChemPot = 0._RK
!DEBUG
        this%Component(i)%ChemPot1 = 0._RK
        this%Component(i)%ChemPot2 = 0._RK
!DEBUG
      end do

      if( NVTEquilibration ) return

    else
      ! Throw test particles
      do i = 1, this%NTestMax
        do j = 1, 3
          this%P0Test( i, j ) = tprnd( -.5_RK, .5_RK )
        end do
        do
          qsum = 0._RK
          do j = 1, 4
            this%Q0Test(i, j) = tprnd( -1._RK, 1._RK )
          end do
          qsum = sum( this%Q0Test(i, :)**2 )
          if( qsum <= 1._RK ) exit
        end do
#if ARCH == 3
        this%Q0Test(i, :) = this%Q0Test(i, :) * rsqrt( qsum )
#else
        this%Q0Test(i, :) = this%Q0Test(i, :) / sqrt( qsum )
#endif
      end do
    end if

    ! Outer loop over components
    do i = 1, this%NRealComponents

      pc => this%Component(i)
      if( Equilibration .and. pc%WFMethod .ne. WFMethodGuess ) cycle
      select case( pc%ChemPotMethod )

      ! Chemical potential by gradual insertion
      case( ChemPotMethodGradIns )

        ! Reset variables
        if( Step == 1 ) then
          pc%ProbW0 = 0._RK
          pc%ProbW1 = 0._RK
          pc%ProbW0V = 0._RK
          pc%ProbW1Rho = 0._RK
          pc%NStateWF(:) = 0
!           pc%NStateBF(:) = 0
        end if

        if( mod( Step, GradInsFrequency ) == 0 ) then
          pc%CalcChemPot = .true.

          ! Save current state
!           call SaveState( this )

          ndfmove = this%NDF
          ndfbiased = this%NDF * 50
          ndffluct = this%NDF * 10
          ndfchange = this%NDF * 10
          ndfcp = ndfmove + ndfbiased + ndffluct + ndfchange
          pc%NState(:) = 1
!DEBUG
  nstate = 0
!   counter = 0
!DEBUG

          ! Set fluctuating particle
          ncf = pc%NFluctComp( pc%NFluctState )
          if( ncf == i ) then
            npf = rnd( pc%NPart )
            call Move2End( this, ncf, npf )
          else
            npf = 1
          end if

giloop:   do j = 1, NFullFluct * ndfcp
! giloop:   do j = 1, 1

            ! Choose particle randomly
            s = 0
            r = rnd( ndfcp )

            if( r <= ndfmove ) then
loop1:        do nc = 1, this%NComponents
                s = s + this%Component(nc)%NDF
                if( r <= s ) exit loop1
              end do loop1
              ndf = this%Component(nc)%Molecule%NDF
              np = 1 + (s - r) / ndf

              ! Move or rotate
              if( mod( s - r, ndf ) < 3 ) then
                call Move( this, nc, np )
              else
                call Rotate( this, nc, np )
              end if

            else if( r <= (ndfmove + ndfbiased) ) then
              r = (r - ndfmove - 1) / 50 + 1
loop2:        do nc = 1, this%NComponents
                s = s + this%Component(nc)%NDF
                if( r <= s ) exit loop2
              end do loop2
              ndf = this%Component(nc)%Molecule%NDF
              np = 1 + (s - r) / ndf

              ! Move or rotate biased
              if( mod( s - r, ndf ) < 3 ) then
                call MoveBiased( this, nc, np, ncf, npf )
              else
                call RotateBiased( this, nc, np, ncf, npf )
              end if

            else if( r <= (ndfmove + ndfbiased + ndffluct) ) then
              r = r - ndfmove - ndfbiased
              ndf = this%Component(ncf)%Molecule%NDF

              ! Move or rotate fluctuating particle
              if( mod( r, ndf ) < 3 ) then
                call Move( this, ncf, npf )
              else
                call Rotate( this, ncf, npf )
              end if

            else
              ! Change fluctuating particle
              call ChangeFluct( this, i, ncf, npf )
!DEBUG
!               pc%NState(pc%NFluctState) = pc%NState(pc%NFluctState) + 1
!               pc%NStateWF(pc%NFluctState) = pc%NStateWF(pc%NFluctState) + 1
!   if( maxcounter > 0 .and. counter > maxcounter ) exit giloop
!   if( maxcounter > 0 ) exit giloop
!   nstate(pc%NFluctState) = nstate(pc%NFluctState) + 1
!DEBUG

            end if

          end do giloop

!DEBUG
  pc%NStateWF = pc%NStateWF + nstate(0:pc%NFluctMax)
!   if( maxcounter > 0 .and. counter > maxcounter ) then
  if( maxcounter > 0 ) then
    write( IOBuffer, &
&     '("GradIns abgebrochen. Zeitschritt: ", I0, "  FluctState: ", I0)') &
&     Step, pc%NFluctState
    call LogWrite
    write( IOBuffer, '("NStates:")' )
    call LogWrite
    do j = 0, pc%NFluctMax
      write( IOBuffer, '(I10)' ) nstate(j)
      call LogWrite
    end do
  else
    pc%NState = pc%NState + nstate(0:pc%NFluctMax)
  end if
!DEBUG

          ! Reset fluctuating particle
!           if( .not. ncf == i ) then
!             pcf => this%Component( ncf )
!             if( pcf%Molecule%IsElongated ) then
!               call AddParticle( pc, pcf%P0( npf, : ), pcf%Q0( npf, : ) )
!             else
!               call AddParticle( pc, pcf%P0( npf, : ) )
!             end if
!             call RemoveParticle( pcf, npf )
!             pc%NFluctState = 0
!           end if

          ! Restore saved state
!           call RestoreState( this )

          ! Calculate weighted propabilities
          pc%ProbW0 = pc%ProbW0 + real(pc%NState(0), RK)
          pc%ProbW1 = pc%ProbW1 &
&           + real( pc%NState(pc%NFluctMax), RK ) / pc%WF( pc%NFluctMax )
          pc%ProbW0V = pc%ProbW0V &
&           + real(pc%NState(0), RK) / this%Density
          pc%ProbW1Rho = pc%ProbW1Rho &
&           + real( pc%NState(pc%NFluctMax), RK ) / pc%WF(pc%NFluctMax) &
&             * this%Density

          ! Calculate chemical potential
          ! (long range correction already done in ChangeFluct)
          pc%ChemPot = pc%ProbW0 / pc%ProbW1Rho
!DEBUG
          pc%ChemPot1 = pc%ProbW0 / pc%ProbW1
          pc%ChemPot2 = pc%ProbW0V / pc%ProbW1
!DEBUG

        else
          pc%CalcChemPot = .false.
          pc%ChemPot = 0._RK
        end if

        if( mod( Step, ErrorsUpdateFrequency ) == 0 .or. &
&           ( Equilibration .and. Step == NStepsP ) ) then
          do j = 1, pc%NFluctMax
            pc%WF(j) = pc%WF(j) * real(pc%NStateWF(0) + 1, RK) &
&                               / real(pc%NStateWF(j) + 1, RK)
          end do
          write( IOBuffer, '("New weighting factors for ",A," calculated:")' ) &
&           trim( pc%PotModFileName )
          call LogWrite
          write( IOBuffer, &
&           '("   State      NState      new WF     up        down (%)")' )
          call LogWrite
          write( IOBuffer, &
&           '("   --------------------------------  --------  --------")' )
          call LogWrite
          write( IOBuffer, &
&           '(I8, I12, F15.2, 2F10.4)' ) 0, pc%NStateWF(0), pc%WF(0), &
&           real(pc%NFluctUpSuccesses(1), RK) / &
&             real(pc%NFluctUpAttempts(1), RK) * 100._RK, 0._RK
          call LogWrite
          do j = 1, pc%NFluctMax - 1
            write( IOBuffer, '(I8, I12, F15.2, 2F10.4)' ) j, pc%NStateWF(j), &
&             pc%WF(j), real(pc%NFluctUpSuccesses(j+1), RK) / &
&               real(pc%NFluctUpAttempts(j+1), RK) * 100._RK, &
&             real(pc%NFluctDownSuccesses(j), RK) / &
&               real(pc%NFluctDownAttempts(j), RK) * 100._RK
            call LogWrite
          end do
          j = pc%NFluctMax
          write( IOBuffer, '(I8, I12, F15.2, 2F10.4)' ) j, pc%NStateWF(j), &
&           pc%WF(j), 0._RK, real(pc%NFluctDownSuccesses(j), RK) / &
&             real(pc%NFluctDownAttempts(j), RK) * 100._RK
            call LogWrite
          call LogWriteBlank
          pc%NStateWF(:) = 0
        end if

      ! Chemical potential by Widom's test particle method
      case( ChemPotMethodWidom )
        pc%CalcChemPot = .true.
        call Mol2AtomTest( this%Component(i), this%Component(i)%NTest )

        if ((LongRange .eq. RField) .or. (LongRange .eq. ExtRField)) then
          this%EPotTest(:) = this%Density * pc%EPotTestCorrLJ &
&                                       + pc%EPotTestCorrRF
          do j = 1, this%NRealComponents
            call ChemicalPotential( this%Interaction( i, j ), &
&                                 this%EPotTest, this%BoxLength )
          end do


! Ewald Summation
        else           ! Ewald
! Chemical Potential  - BLODSINN, nur fuer restart wichtig
           write (*,*) 'Widom does not yet work with Ewald Summation'
           STOP
!          DO j=1,this%NComponents
           call Ewald_ChemPotSelf(this,i)
!          END DO

          this%EPotTest(:) = this%Density * pc%EPotTestCorrLJ + &
&                            pc%EPotTestSelf
          do j = 1, this%NRealComponents
            call ChemicalPotential( this%Interaction( i, j ), &
&                                 this%EPotTest, this%BoxLength )
          end do

! ! !           call Ewald_ChemPotFour(this,i)
        end if

        ChemPot = sum( exp( -( this%EPotTest(:) ) / this%Temperature ) ) &
&                   / pc%NTestAll
#if MPI_VER > 0
        call MPI_Reduce( ChemPot, pc%ChemPot, 1, &
&         MPI_DOUBLE_PRECISION, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
#else
        pc%ChemPot = ChemPot
#endif

      case default
        pc%CalcChemPot = .false.
        pc%ChemPot = 0._RK
      end select

    end do



    contains



    function tprnd( l_range, h_range ) result( rharvest )

      implicit none

      ! Include MPI header
#if MPI_VER > 1
      include 'mpif.h'
#endif

      ! Declare arguments
      real(RK), intent(in) :: l_range, h_range

      ! Declare result
      real(RK) :: rharvest

      ! Declare local variables
      integer(K4B), parameter :: IA=16807, IM=2147483647, &
&       IQ=127773, IR=2836

      ! Generate fast random number
      tpix = ieor(tpix, ishft(tpix, 13))
      tpix = ieor(tpix, ishft(tpix, -17))
      tpix = ieor(tpix, ishft(tpix, 5))
      rharvest = l_range + am * ior(iand(IM, tpix), 1) &
&       * (h_range - l_range)

    end function tprnd



  end subroutine TEnsemble_ChemicalPotential



!==============================================================!
!  Subroutine TEnsemble_UpdateEnergy                           !
!==============================================================!

  subroutine TEnsemble_UpdateEnergy( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this

    ! Declare local variables
    type(TInteraction), pointer :: pi
    integer                     :: n1, n2
    integer                     :: i, j

    ! Update potential energy and virial matrices
    do i = 1, this%NComponents
      do j = 1, this%NComponents
        pi => this%Interaction(j, i)
        n1 = pi%NPart1
        n2 = pi%NPart2
        pi%EPot(1:n1, 1:n2) = pi%EPotNew(1:n1, 1:n2)
        pi%Virial(1:n1, 1:n2) = pi%VirialNew(1:n1, 1:n2)
      end do
    end do

  end subroutine TEnsemble_UpdateEnergy



!==============================================================!
!  Subroutine TEnsemble_UpdateEnergy1                          !
!==============================================================!

  subroutine TEnsemble_UpdateEnergy1( this, nc, np )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np

    ! Declare local variables
    type(TInteraction), pointer :: pi
    integer                     :: n
    integer                     :: i

    ! Update potential energy and virial matrices for a particle
    do i = 1, this%NComponents
      pi => this%Interaction(nc, i)
      n = pi%NPart2
      pi%EPot(np, 1:n) = pi%EPot1(1:n)
      pi%Virial(np, 1:n) = pi%Virial1(1:n)
      this%Interaction(i, nc)%EPot(1:n, np) = pi%EPot1(1:n)
      this%Interaction(i, nc)%Virial(1:n, np) = pi%Virial1(1:n)
    end do

  end subroutine TEnsemble_UpdateEnergy1



!==============================================================!
!  Subroutine TEnsemble_Energy                                 !
!==============================================================!

  subroutine TEnsemble_Energy( this, E )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)       :: this
    real(RK), intent(out) :: E

    ! Declare local variables
    type(TInteraction), pointer :: pi
    integer                     :: nc, np
    integer                     :: i, n

    ! Initialize new energy
    E = 0._RK

    if (LongRange .eq. Ewald) then
      call EwaldSelfTerm_Energy ( this )
    else if (LongRange .eq. PME) then
      call PMESelfTermMC ( this )
    end if
    ! Loop over components
    do nc = 1, this%NComponents
      do i = 1, this%NComponents
        pi => this%Interaction(nc, i)
        n = pi%NPart2

        ! Loop over particles
        do np = 1, this%Component(nc)%NPart
          call Energy( pi, np, this%BoxLength )

!           if ((LongRange .eq. Ewald) .AND. (i .eq. nc)) then
!              call EwaldSelfTerm_Energy(this,nc,np)	! immer gleich - Beschleunigung möglich!!!!!!!!!!
!           end if
          ! Save new energy matrix
          pi%EPotNew(np, 1:n) = pi%EPot1(1:n)
          pi%VirialNew(np, 1:n) = pi%Virial1(1:n)

          ! Sum energy
          E = E + sum( pi%EPot1(1:n) )
        end do
      end do
    end do

    ! Calculate new energy
    E = .5_RK * E + this%Density * this%EPotCorrLJ + this%EPotCorrRF

! Ewald 
    if (LongRange .eq. Ewald) then
      call EwaldFourierEnergy(this)
      E = E + this%UFourier + this%UIntra + this%USelbstTerm
    else if (LongRange .eq. PME) then
      call charge_grid_MCall ( this )
      call PMEFourierTermMC ( this )
      E = E + this%UFourier + this%UIntra + this%USelbstTerm
    end if

  end subroutine TEnsemble_Energy



!==============================================================!
!  Subroutine TEnsemble_Energy1                                !
!==============================================================!

  subroutine TEnsemble_Energy1( this, nc, np, EPotNew )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)       :: this
    integer, intent(in)   :: nc, np
    real(RK), intent(out) :: EPotNew

    ! Declare local variables
    type(TInteraction), pointer :: pi
    integer                     :: n
    integer                     :: i

    ! Initialize new energy
    EPotNew = 0._RK

!     if ((LongRange .eq. Ewald) .AND. (i .eq. nc)) then
!        call EwaldSelfTerm_Energy(this)
!     end if

    ! Loop over components
    do i = 1, this%NComponents
      pi => this%Interaction(nc, i)
      n = pi%NPart2

      call Energy( pi, np, this%BoxLength )


      ! Calculate new energy
      EPotNew = EPotNew + sum( pi%EPot1(1:n) )
    end do

    if (LongRange .eq. Ewald) then
       call EwaldFourierEnergy(this,nc,np)
       EPotNew = EPotnew + this%UFourier
    else if (LongRange .eq. PME) then
       call PMEFourierTermMC ( this )
       EPotNew = EPotnew + this%UFourier
    end if
  end subroutine TEnsemble_Energy1



!==============================================================!
!  Subroutine TEnsemble_Energy_CF                              !
!==============================================================!

  subroutine TEnsemble_Energy1_CF( this, nc, np, ncold, npold, EPotNew )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)       :: this
    integer, intent(in)   :: nc, np
    integer, intent(in)   :: ncold, npold
    real(RK), intent(out) :: EPotNew

    ! Declare local variables
    type(TInteraction), pointer :: pi
    integer                     :: n
    integer                     :: i

    ! Initialize new energy
    EPotNew = 0._RK

    ! Loop over components
    do i = 1, this%NComponents
      pi => this%Interaction(nc, i)
      n = pi%NPart2

      call Energy( pi, np, this%BoxLength )

      ! Calculate new energy
      EPotNew = EPotNew + sum( pi%EPot1(1:n) )
    end do

    if (LongRange .eq. Ewald) then
       call EwaldFourierEnergy(this,nc,np,ncold,npold)
       EPotNew = EPotnew + this%UFourier + this%USelbstTerm + this%UIntra
    else if (LongRange .eq. PME) then
       call PMEFourierTermMC ( this )
       EPotNew = EPotnew + this%UFourier
    end if
  end subroutine TEnsemble_Energy1_CF



!==============================================================!
!  Function TEnsemble_GetEnergy                                !
!==============================================================!

  function TEnsemble_GetEnergy( this ) result(E)

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare result
    real(RK) :: E

    ! Declare local variables
    integer :: i, j
    integer :: n

    ! Calculate potential energy of a particle
    E = 0._RK
    do i = 1, this%NComponents
      n = this%Component(i)%NPart
      do j = 1, this%NComponents
        E = E + sum( this%Interaction(j, i)% &
&         EPot(1:this%Component(j)%NPart, 1:n) )
      end do
    end do
    E = .5_RK * E + this%Density * this%EPotCorrLJ + this%EPotCorrRF

! Ewald 
    if (LongRange .eq. Ewald) then
      call EwaldFourierEnergy(this)
      E = E + this%UFourier + this%UIntra + this%USelbstTerm
    else if (LongRange .eq. PME) then
      call charge_grid_MCall (this)
      call PMEFourierTermMC(this)
      E = E + this%UFourier + this%UIntra + this%USelbstTerm
    end if

  end function TEnsemble_GetEnergy



!==============================================================!
!  Function TEnsemble_GetEnergy1                               !
!==============================================================!

  function TEnsemble_GetEnergy1( this, nc, np ) result(E)

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np

    ! Declare result
    real(RK) :: E

    ! Declare local variables
    integer :: i

    ! Calculate potential energy of a particle
    E = 0._RK
    do i = 1, this%NComponents
      E = E &
&         + sum( this%Interaction(i, nc)%EPot(1:this%Component(i)%NPart, np) )
    end do

! Ewald 
    if ((LongRange .eq. Ewald) .or. (LongRange .eq. PME)) then
      E = E + this%UFourier
    end if

  end function TEnsemble_GetEnergy1



!==============================================================!
!  Function TEnsemble_GetVirial                                !
!==============================================================!

  function TEnsemble_GetVirial( this ) result(V)

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare result
    real(RK) :: V

    ! Declare local variables
    integer :: i, j
    integer :: n

    ! Calculate potential energy of a particle
    V = 0._RK
    do i = 1, this%NComponents
      n = this%Component(i)%NPart
      do j = 1, this%NComponents
        V = V + sum( this%Interaction(j, i)% &
&         Virial(1:this%Component(j)%NPart, 1:n) )
      end do
    end do
    V = .5_RK * V + this%Density * this%VirialCorrLJ

    if (LongRange .eq. Ewald) then
!       call EwaldFourierEnergy(this)
       V = V + this%EVirial
    else if (LongRange .eq. PME) then
       V = V + this%EVirial
    end if

  end function TEnsemble_GetVirial



!==============================================================!
!  Subroutine TEnsemble_Move                                   !
!==============================================================!

  subroutine TEnsemble_Move( this, nc, np )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np

    ! Declare local variables
    real(RK)                  :: r(3)
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier, EVirial
    type(TComponent), pointer :: pc
    integer                   :: i
#if MPI_VER > 0
    real(RK)                  :: EPotDeltaAll
#endif

    ! Assign local variables
    pc => this%Component(nc)

    ! Update number of move attempts
    pc%NMoveAttempts = pc%NMoveAttempts + 1

    ! Save current particle position and energy
    r(:) = pc%P0(np, :)
    EPotOld = GetEnergy( this, nc, np )

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO
!       EVirial  = this%EVirial
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
!       this%qgrida_old = this%qgrida
      call chargegrid_min  (this, nc, np)
    end if

    ! Generate a trial displacement
    do i = 1, 3
      pc%P0(np, i) = pc%P0(np, i) + rnd( -pc%DispTran, pc%DispTran )
    end do

    ! Apply periodic boundary conditions
    pc%P0(np, :) = pc%P0(np, :) - anint( pc%P0(np, :) )

    ! Convert molecular coordinates to atom positions
    call Mol2Atom1( pc, np )

    ! Calculate changes in the SPME grid
    if (LongRange .eq. PME) then
      call chargegrid_plus (this, nc, np)
!       call charge_grid_MCall (this )
    end if

    ! Calculate particle energy at trial position
    call Energy( this, nc, np, EPotNew )
    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    call MPI_Allreduce( EPotOld - EPotNew, EPotDeltaAll, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
    if( exp( EPotDeltaAll / this%Temperature ) .gt. rnd( 0._RK, 1._RK ) ) then
#else
    if( exp( (EPotOld - EPotNew) / this%Temperature ) &
&       .gt. rnd( 0._RK, 1._RK ) ) then
#endif

      ! Accept move
      pc%NMoveSuccesses = pc%NMoveSuccesses + 1
      call UpdateEnergy( this, nc, np )

    else

      ! Reject move
      if (LongRange .eq. Ewald) then
          this%UFourier = EFourier
!         this%EVirial  = EVirial
          DO i=1,pc%Molecule%NCharge
            this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
            this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
            this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
          END DO
          pc%P0(np, :) = r(:)
          call Mol2Atom1( pc, np )
          call EwaldFourierEnergy(this,nc,np)
!           this%sinfac_s(nc,1:this%Component(nc)%Molecule%NCharge,np) = this%sinfac_s_old
!           this%cosfac_s(nc,1:this%Component(nc)%Molecule%NCharge,np) = this%cosfac_s_old
      else if (LongRange .eq. PME) then
          this%UFourier = EFourier
          this%EVirial  = EVirial
          call chargegrid_min  (this, nc, np)
          pc%P0(np, :) = r(:)
          call Mol2Atom1( pc, np )
          call chargegrid_plus (this, nc, np)
!         this%qgrida   = this%qgrida_old
      else
          pc%P0(np, :) = r(:)
          call Mol2Atom1( pc, np )
      end if
    end if

  end subroutine TEnsemble_Move



!==============================================================!
!  Subroutine TEnsemble_Rotate                                 !
!==============================================================!

  subroutine TEnsemble_Rotate( this, nc, np )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments03
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np

    ! Declare local variables
    real(RK)                  :: q(4), dq(3)
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier, EVirial
    type(TComponent), pointer :: pc
    integer                   :: i
#if MPI_VER > 0
    real(RK)                  :: EPotDeltaAll
#endif

    ! Assign local variables
    pc => this%Component(nc)

    ! Update number of rotation attempts
    pc%NRotateAttempts = pc%NRotateAttempts + 1

    ! Save current particle orientation and energy
    q(:) = pc%Q0(np, :)
    EPotOld = GetEnergy( this, nc, np )

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO
!       EVirial  = this%EVirial
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
!       this%qgrida_old = this%qgrida
      call chargegrid_min  (this, nc, np)
    end if

    ! Generate a trial rotation
    do i = 1, 3
      dq(i) = rnd( -pc%DispRot, pc%DispRot )
    end do
    pc%Q0(np, 1) = q(1) - dq(1) * q(2) - dq(2) * q(3) - dq(3) * q(4)
    pc%Q0(np, 2) = q(2) + dq(1) * q(1) - dq(2) * q(4) + dq(3) * q(3)
    pc%Q0(np, 3) = q(3) + dq(1) * q(4) + dq(2) * q(1) - dq(3) * q(2)
    pc%Q0(np, 4) = q(4) - dq(1) * q(3) + dq(2) * q(2) + dq(3) * q(1)


    ! Convert molecular coordinates to atom positions
    call Mol2Atom1( pc, np )

    if (LongRange .eq. PME) then
      call chargegrid_plus (this, nc, np)
!       call charge_grid_MCall (this)
    end if

    ! Calculate particle energy with trial orientation
    call Energy( this, nc, np, EPotNew )

    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    call MPI_Allreduce( EPotOld - EPotNew, EPotDeltaAll, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
    if( exp( EPotDeltaAll / this%Temperature ) .gt. rnd( 0._RK, 1._RK ) ) then
#else
    if( exp( (EPotOld - EPotNew) / this%Temperature ) &
&       .gt. rnd( 0._RK, 1._RK ) ) then
#endif

      ! Accept rotation
      pc%NRotateSuccesses = pc%NRotateSuccesses + 1
      call UpdateEnergy( this, nc, np )

    else

      ! Reject move
      if (LongRange .eq. Ewald) then
        this%UFourier = EFourier
        DO i=1,pc%Molecule%NCharge
          this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
          this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
          this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
        END DO
        pc%Q0(np, :) = q(:)
        call Mol2Atom1( pc, np )
        call EwaldFourierEnergy(this,nc,np)
!         this%sinfac_s(nc,1:this%Component(nc)%Molecule%NCharge,np) = this%sinfac_s_old
!         this%cosfac_s(nc,1:this%Component(nc)%Molecule%NCharge,np) = this%cosfac_s_old
      else if (LongRange .eq. PME) then
        this%UFourier = EFourier
        this%EVirial  = EVirial
        call chargegrid_min  (this, nc, np)
        pc%Q0(np, :) = q(:)
        call Mol2Atom1( pc, np )
        call chargegrid_plus (this, nc, np)
!         this%qgrida   = this%qgrida_old
      else
        pc%Q0(np, :) = q(:)
        call Mol2Atom1( pc, np )
      end if

    end if

  end subroutine TEnsemble_Rotate



!==============================================================!
!  Subroutine TEnsemble_MoveBiased                             !
!==============================================================!

  subroutine TEnsemble_MoveBiased( this, nc, np, ncf, npf )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np, ncf, npf

    ! Declare local variables
    real(RK)                  :: r(3), dr(3), f1, f2
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier, EVirial
    type(TComponent), pointer :: pc, pcf
    integer                   :: i
#if MPI_VER > 0
    real(RK)                  :: EPotDeltaAll
#endif

    ! Test for fluctuating particle
    if( nc .eq. ncf .and. np .eq. npf ) return

    ! Assign local variables
    pc => this%Component(nc)
    pcf => this%Component(ncf)

    ! Update number of move attempts
    pc%NMoveBiasedAttempts = pc%NMoveBiasedAttempts + 1

    ! Save current particle position and energy
    r(:) = pc%P0(np, :)
    EPotOld = GetEnergy( this, nc, np )

    ! Apply distance criterion
    dr(:) = r(:) - pcf%P0(npf, :)
    dr(:) = ( dr(:) - anint( dr(:) ) ) * this%BoxLength
    f1 = 1._RK / ( dr(1)**2 + dr(2)**2 + dr(3)**2 )**2
    if( rnd(0._RK, 1._RK) < (1._RK - f1) ) return

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO
!       EVirial  = this%EVirial
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      this%qgrida_old = this%qgrida
      call chargegrid_min  (this, nc, np)
    end if

    ! Generate a trial displacement
    do i = 1, 3
      pc%P0(np, i) = pc%P0(np, i) + rnd( -pc%DispTran, pc%DispTran )
    end do

    ! Apply periodic boundary conditions
    pc%P0(np, :) = pc%P0(np, :) - anint( pc%P0(np, :) )

    ! Apply direction criterion
    dr(:) = pc%P0(np, :) - pcf%P0(npf, :)
    dr(:) = ( dr(:) - anint( dr(:) ) ) * this%BoxLength
    f2 = 1._RK / ( dr(1)**2 + dr(2)**2 + dr(3)**2 )**2
    if( rnd(0._RK, 1._RK) < (1._RK - f2/f1) ) then
      pc%P0(np, :) = r(:)
      return
    end if

    ! Convert molecular coordinates to atom positions
    call Mol2Atom1( pc, np )

    ! Save Energies, Virials for faster Rejection
    if (LongRange .eq. PME) then
      call chargegrid_plus (this, nc, np)
    end if

    ! Calculate particle energy at trial position
    call Energy( this, nc, np, EPotNew )

    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    call MPI_Allreduce( EPotOld - EPotNew, EPotDeltaAll, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
    if( exp( EPotDeltaAll / this%Temperature ) .gt. rnd( 0._RK, 1._RK ) ) then
#else
    if( exp( (EPotOld - EPotNew) / this%Temperature ) &
&       .gt. rnd( 0._RK, 1._RK ) ) then
#endif

      ! Accept move
      pc%NMoveBiasedSuccesses = pc%NMoveBiasedSuccesses + 1
      call UpdateEnergy( this, nc, np )

    else

      ! Reject move
      if (LongRange .eq. Ewald) then
          this%UFourier = EFourier
!         this%EVirial  = EVirial
          DO i=1,pc%Molecule%NCharge
            this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
            this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
            this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
          END DO
          pc%P0(np, :) = r(:)
          call Mol2Atom1( pc, np )
          call EwaldFourierEnergy(this,nc,np)
      else if (LongRange .eq. PME) then
        this%UFourier = EFourier
        this%EVirial  = EVirial
        this%qgrida   = this%qgrida_old
      else
        pc%P0(np, :) = r(:)
        call Mol2Atom1( pc, np )
      end if

    end if

  end subroutine TEnsemble_MoveBiased



!==============================================================!
!  Subroutine TEnsemble_RotateBiased                           !
!==============================================================!

  subroutine TEnsemble_RotateBiased( this, nc, np, ncf, npf )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments03
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np, ncf, npf

    ! Declare local variables
    real(RK)                  :: dr(3), q(4), dq(3), f1
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier, EVirial
    type(TComponent), pointer :: pc, pcf
    integer                   :: i
#if MPI_VER > 0
    real(RK)                  :: EPotDeltaAll
#endif

    ! Test for fluctuating particle
    if( nc .eq. ncf .and. np .eq. npf ) return

    ! Assign local variables
    pc => this%Component(nc)
    pcf => this%Component(ncf)

    ! Update number of rotation attempts
    pc%NRotateBiasedAttempts = pc%NRotateBiasedAttempts + 1

    ! Save current particle orientation and energy
    q(:) = pc%Q0(np, :)
    EPotOld = GetEnergy( this, nc, np )

    ! Apply distance criterion
    dr(:) = pc%P0(np, :) - pcf%P0(npf, :)
    dr(:) = ( dr(:) - anint( dr(:) ) ) * this%BoxLength
    f1 = 1._RK / ( dr(1)**2 + dr(2)**2 + dr(3)**2 )**2
    if( rnd(0._RK, 1._RK) < (1._RK - f1) ) return

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO
!       EVirial  = this%EVirial
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      this%qgrida_old = this%qgrida
      call chargegrid_min  (this, nc, np)
    end if

    ! Generate a trial rotation
    do i = 1, 3
      dq(i) = rnd( -pc%DispRot, pc%DispRot )
    end do
    pc%Q0(np, 1) = q(1) - dq(1) * q(2) - dq(2) * q(3) - dq(3) * q(4)
    pc%Q0(np, 2) = q(2) + dq(1) * q(1) - dq(2) * q(4) + dq(3) * q(3)
    pc%Q0(np, 3) = q(3) + dq(1) * q(4) + dq(2) * q(1) - dq(3) * q(2)
    pc%Q0(np, 4) = q(4) - dq(1) * q(3) + dq(2) * q(2) + dq(3) * q(1)

    ! Convert molecular coordinates to atom positions
    call Mol2Atom1( pc, np )

    ! Save Energies, Virials for faster Rejection
    if (LongRange .eq. PME) then
      call chargegrid_plus (this, nc, np)
    end if

    ! Calculate particle energy with trial orientation
    call Energy( this, nc, np, EPotNew )

    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    call MPI_Allreduce( EPotOld - EPotNew, EPotDeltaAll, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
    if( exp( EPotDeltaAll / this%Temperature ) .gt. rnd( 0._RK, 1._RK ) ) then
#else
    if( exp( (EPotOld - EPotNew) / this%Temperature ) &
&       .gt. rnd( 0._RK, 1._RK ) ) then
#endif

      ! Accept rotation
      pc%NRotateBiasedSuccesses = pc%NRotateBiasedSuccesses + 1
      call UpdateEnergy( this, nc, np )

    else

      ! Reject move
      if (LongRange .eq. Ewald) then
        this%UFourier = EFourier
        DO i=1,pc%Molecule%NCharge
          this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
          this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
          this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
        END DO
        pc%Q0(np, :) = q(:)
        call Mol2Atom1( pc, np )
        call EwaldFourierEnergy(this,nc,np)
      else if (LongRange .eq. PME) then
        this%UFourier = EFourier
        this%EVirial  = EVirial
        this%qgrida   = this%qgrida_old
      else
        pc%Q0(np, :) = q(:)
        call Mol2Atom1( pc, np )
      end if

    end if

  end subroutine TEnsemble_RotateBiased



!==============================================================!
!  Subroutine TEnsemble_ChangeFluct                            !
!==============================================================!

  subroutine TEnsemble_ChangeFluct( this, nc, ncf, npf )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)        :: this
    integer, intent(in)    :: nc
    integer, intent(inout) :: ncf, npf
    integer                :: i
!DEBUG
!     integer, intent(inout) :: counter
!DEBUG

    ! Declare local variables
    type(TComponent), pointer :: pc, pcf, pcfnew
    integer                   :: oldstate, newstate
    integer                   :: ncfnew, npfnew
    real(RK)                  :: EPotOld, EPotNew
#if MPI_VER > 0
    real(RK)                  :: EPotDeltaAll
#endif
    real(RK)                  :: EFourier, EVirial
!DEBUG
!   logical :: accepted, unequal
!DEBUG

    ! Assign local variables
    pc => this%Component(nc)
    pcf => this%Component(ncf)
    oldstate = pc%NFluctState

    if( oldstate .eq. 0 ) then
      if( rnd( 0._RK, 1._RK ) < .5_RK ) then
        npf = rnd( pcf%NPart )
        call Move2End( this, ncf, npf )
!DEBUG
!   counter = 0
!DEBUG
      end if
      newstate = 1
      pc%NFluctUpAttempts( newstate ) = pc%NFluctUpAttempts( newstate ) + 1
    elseif( oldstate .eq. pc%NFluctMax ) then
      newstate = oldstate - 1
      pc%NFluctDownAttempts( oldstate ) = pc%NFluctDownAttempts( oldstate ) + 1
    else
      if( rnd( 0._RK, 1._RK ) < .5_RK ) then
        newstate = oldstate + 1
        pc%NFluctUpAttempts( newstate ) = pc%NFluctUpAttempts( newstate ) + 1
      else
        newstate = oldstate - 1
        pc%NFluctDownAttempts( oldstate ) = pc%NFluctDownAttempts( oldstate ) + 1
      end if
    end if

!DEBUG
!   counter = counter + 1
!DEBUG

    ! Get old energy of fluctuating particle
    EPotOld = GetEnergy( this, ncf, npf )

    ! Change state of fluctuating particle
    ncfnew = pc%NFluctComp( newstate )
    pcfnew => this%Component( ncfnew )
    if( pcf%Molecule%IsElongated ) then
      call AddParticle( pcfnew, pcf%P0( npf, : ), pcf%Q0( npf, : ) )
    else
      call AddParticle( pcfnew, pcf%P0( npf, : ) )
    end if
    call RemoveParticle( pcf, npf )
    npfnew = pcfnew%NPart

! Save states for the Ewald Summation and/or derivates
    if (LongRange .eq. Ewald) then     ! Ewald Summation
       ! Save the initial state
       EFourier = this%UFourier
       EPotOld = EPotOld  + this%USelbstTerm + this%UIntra
       EVirial  = this%EVirial
       DO i=1,pcf%Molecule%NCharge
         this%rold(i,1) = pcf%Molecule%SiteCharge(i)%RX(npf)
         this%rold(i,2) = pcf%Molecule%SiteCharge(i)%RY(npf)
         this%rold(i,3) = pcf%Molecule%SiteCharge(i)%RZ(npf)
       END DO

       ! Calculate new energies
       call EwaldSelfTerm_Energy(this)

       ! Convert molecular coordinates to atom positions
       call Mol2Atom1( pcfnew, npfnew )

       ! Calculate particle energy at new fluctuating state
       call Energy( this, ncfnew, npfnew, ncf, npf, EPotNew )

       ! Acceptance Criteria
#if MPI_VER > 0
       call MPI_Allreduce( EPotOld - EPotNew, EPotDeltaAll, 1, &
&        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
       EPotDeltaAll = EPotDeltaAll &
&        + this%Density * ( pcf%EPotTestCorrLJ - pcfnew%EPotTestCorrLJ )
       if( rnd( 0._RK, 1._RK ) < pc%WF(newstate) / pc%WF(oldstate) * &
&        exp( ( EPotDeltaAll ) / this%Temperature ) ) then
#else
       EPotOld = EPotOld &
&        + this%Density * ( pcf%EPotTestCorrLJ - pcfnew%EPotTestCorrLJ )
       if( rnd( 0._RK, 1._RK ) < pc%WF(newstate) / pc%WF(oldstate) * &
&        exp( ( EPotOld - EPotNew ) / this%Temperature ) ) then
#endif

         ! Accept
         pc%NFluctState = newstate
         ncf = ncfnew
         npf = npfnew
         call UpdateEnergy( this, ncf, npf )
         if( newstate > oldstate ) then
           pc%NFluctUpSuccesses( newstate ) = pc%NFluctUpSuccesses( newstate ) + 1
         else
           pc%NFluctDownSuccesses( oldstate ) = pc%NFluctDownSuccesses( oldstate ) + 1
         end if

       else
         ! Reject
         if( pcf%Molecule%IsElongated ) then
           call AddParticle( pcf, pcfnew%P0( npfnew, : ), pcfnew%Q0( npfnew, : ) )
         else
           call AddParticle( pcf, pcfnew%P0( npfnew, : ) )
         end if
         call RemoveParticle( pcfnew, npfnew )
         call EwaldSelfTerm_Energy(this)
         DO i=1,pcfnew%Molecule%NCharge
           this%rold(i,1) = pcfnew%Molecule%SiteCharge(i)%RX(npfnew)
           this%rold(i,2) = pcfnew%Molecule%SiteCharge(i)%RY(npfnew)
           this%rold(i,3) = pcfnew%Molecule%SiteCharge(i)%RZ(npfnew)
         END DO
         call Energy( this, ncf, npf, ncfnew, npfnew, EPotNew )

       end if       ! Acceptance Criteria

! ----------------------------------------------------------------
    else if (LongRange .eq. PME) then ! PME 
      EFourier = this%UFourier
      EVirial  = this%EVirial
      call PMESetup(this)
      write (*,*) 'Gradual Insertion does not yet work with PME'
      STOP

! ----------------------------------------------------------------
    else   ! REACTION FIELD
       ! Convert molecular coordinates to atom positions
       call Mol2Atom1( pcfnew, npfnew )

       ! Calculate particle energy at new fluctuating state
       call Energy( this, ncfnew, npfnew, EPotNew )

    ! Apply acceptance criterion
#if MPI_VER > 0
       call MPI_Allreduce( EPotOld - EPotNew, EPotDeltaAll, 1, &
&        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
       EPotDeltaAll = EPotDeltaAll &
&        + this%Density * ( pcf%EPotTestCorrLJ - pcfnew%EPotTestCorrLJ ) &
&        + pcf%EPotTestCorrRF - pcfnew%EPotTestCorrRF
       if( rnd( 0._RK, 1._RK ) < pc%WF(newstate) / pc%WF(oldstate) * &
&        exp( ( EPotDeltaAll ) / this%Temperature ) ) then
#else
       EPotOld = EPotOld &
&        + this%Density * ( pcf%EPotTestCorrLJ - pcfnew%EPotTestCorrLJ ) &
&        + pcf%EPotTestCorrRF - pcfnew%EPotTestCorrRF
       if( rnd( 0._RK, 1._RK ) < pc%WF(newstate) / pc%WF(oldstate) * &
&        exp( ( EPotOld - EPotNew ) / this%Temperature ) ) then
#endif

         ! Accept
         pc%NFluctState = newstate
         ncf = ncfnew
         npf = npfnew
         call UpdateEnergy( this, ncf, npf )
         if( newstate > oldstate ) then
           pc%NFluctUpSuccesses( newstate ) = pc%NFluctUpSuccesses( newstate ) + 1
         else
           pc%NFluctDownSuccesses( oldstate ) = pc%NFluctDownSuccesses( oldstate ) + 1
         end if
!DEBUG
!      counter = 0
!   accepted = .TRUE.
!DEBUG

       else

         ! Reject
         if( pcf%Molecule%IsElongated ) then
           call AddParticle( pcf, pcfnew%P0( npfnew, : ), pcfnew%Q0( npfnew, : ) )
         else
           call AddParticle( pcf, pcfnew%P0( npfnew, : ) )
         end if
         call RemoveParticle( pcfnew, npfnew )

       end if

    end if      ! LongRange - Correction
!
!
!DEBUG
!  counter = i
!   call MPI_AllReduce( accepted, unequal, 1, MPI_LOGICAL, MPI_LXOR, MPI_COMM_WORLD, ierror )
!   if( unequal ) then
!     write(0, '(I2, ": ", A, " EPotDeltaAll=", F20.16, ", ix=", I0)') &
! &     NProc, merge( "acc.", "rej.", accepted ), EPotDeltaAll, ix
!     call MPI_Barrier( MPI_COMM_WORLD, ierror )
!     write(0, '(I2, ": WFnew=", F12.4, " WFold=", F12.4, " Temp=", F12.8)') &
! &     NProc, pc%WF(newstate), pc%WF(oldstate), this%Temperature
!     call MPI_Barrier( MPI_COMM_WORLD, ierror )
!     stop
!   end if
!DEBUG

  end subroutine TEnsemble_ChangeFluct



!==============================================================!
!  Subroutine TEnsemble_Insert                                 !
!==============================================================!

  subroutine TEnsemble_Insert( this, nc )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc

    ! Declare local variables
    real(RK)                  :: r(3)
    real(RK)                  :: q(4)
    real(RK)                  :: EPotIns
    type(TComponent), pointer :: pc
    integer                   :: i, np
    real(RK)                  :: s
    real(RK)                  :: EFourier, EVirial
#if MPI_VER > 0
    real(RK)                  :: EPotInsAll
!DEBUG
!  logical                   :: accepted, different
!DEBUG
#endif

    ! Assign local variables
    pc => this%Component(nc)

    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      this%qgrida_old = this%qgrida
    end if
    ! Update number of insert attempts
    this%NInsertAttempts = this%NInsertAttempts + 1

    ! Generate a random position and orientation
    do i = 1, 3
      r(i) = rnd( -.5_RK, .5_RK )
    end do
    do
      s = 0._RK
      do i = 1, 4
        q(i) = rnd( -1._RK, 1._RK )
      end do
      s = sum( q**2 )
      if( s <= 1._RK ) exit
    end do
#if ARCH == 3
    q = q * rsqrt( s )
#else
    q = q / sqrt( s )
#endif

    call AddParticle( pc, r, q )
    np = pc%NPart
    this%NPart = this%NPart + 1

    ! Convert molecular coordinates to atom positions
    call Mol2Atom1( pc, np )

    ! Save Energies, Virials for faster Rejection
    if (LongRange .eq. PME) then
      call chargegrid_plus (this, nc, np)
      call PMESelfTermMC( this )
    end if

    ! Calculate particle energy at trial position
    call Energy( this, nc, np, EPotIns )

    ! Apply acceptance criterion
#if MPI_VER > 0
    call MPI_Allreduce( EPotIns, EPotInsAll, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
    EPotInsAll = EPotInsAll + this%Density * pc%EPotTestCorrLJ &
&                           + pc%EPotTestCorrRF
    if ( LongRange .eq. PME ) then
      EPotInsAll = EPotInsAll + this%USelbstTerm + this%UIntra
    end if
!DEBUG
!  write(0, '(I2, ": EPotIns = ", F12.6)') NProc, EPotInsAll
!DEBUG
    if( rnd( 0._RK, 1._RK ) .lt. &
&       ( exp( pc%ChemPot - EPotInsAll / this%Temperature ) &
&         * this%Volume0 / np )) then
#else
    EPotIns = EPotIns + this%Density * pc%EPotTestCorrLJ &
&                     + pc%EPotTestCorrRF
    if (LongRange .eq. PME) then
      EPotIns = EPotIns + this%USelbstTerm + this%UIntra
    end if
!DEBUG
!  write(0, '(I2, ": EPotIns = ", F12.6)') NProc, EPotIns
!DEBUG
    if( rnd( 0._RK, 1._RK ) .lt. &
&       ( exp( pc%ChemPot - EPotIns / this%Temperature ) &
&         * this%Volume0 / np )) then
#endif

      ! Accept Insertion
      this%NInsertSuccesses = this%NInsertSuccesses + 1

      ! Update energy matrix
      call UpdateEnergy( this, nc, np )

      ! Update density
      this%Density = this%NPart / this%Volume0

      ! Update fractions and NDF
      call UpdateFractions( this )

      ! Update long range correction
      call CalculateCorr( this )

!DEBUG
!#if MPI_VER > 0
!  accepted = .TRUE.
!#endif
!DEBUG

    else

      ! Reject Insertion
      call RemoveParticle( pc, np )
      this%NPart = this%NPart - 1

      if (LongRange .eq. Ewald) then
        this%UFourier = EFourier
        this%EVirial  = EVirial
      else if (LongRange .eq. PME) then
        this%UFourier = EFourier
        this%EVirial  = EVirial
        this%qgrida   = this%qgrida_old
        call PMESelfTermMC (this)
      end if

!DEBUG
!#if MPI_VER > 0
!  accepted = .FALSE.
!#endif
!DEBUG

    end if

!DEBUG
!#if MPI_VER > 0
!  call MPI_Allreduce( accepted, different, 1, MPI_LOGICAL, MPI_LXOR, &
!&   MPI_COMM_WORLD, ierror )
!  if( different ) then
!    write(0, '(I2, ": Insert of comp. ", I0, A, " at step ", I0)') &
!&     NProc, nc, merge("    accepted", "not accepted", accepted), step
!    write(0, '(I2, ": Next random number = ", F12.10)') NProc, rnd(0._RK, 1._RK)
!    stop
!  end if
!#endif
!DEBUG

  end subroutine TEnsemble_Insert



!==============================================================!
!  Subroutine TEnsemble_Delete                                 !
!==============================================================!

  subroutine TEnsemble_Delete( this, nc, np )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np

    ! Declare local variables
    real(RK)                    :: EPotDel
    type(TComponent), pointer   :: pc
    type(TInteraction), pointer :: pi
    integer                     :: i, n1, n2
!     real(RK)                    :: s

! Ewald Parameter
    real(RK)                    :: EFourier, EPotNew
    real(RK)                    :: EVirial, EVirialIntra
    real(RK)                    :: USelf, UIntra
    real(RK)                    :: r(3)
    real(RK)                    :: q(4)

!DEBUG
!#if MPI_VER > 0
!  logical                     :: accepted, different
!#endif
!DEBUG

    ! Assign local variables
    pc => this%Component(nc)

    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      USelf    = this%USelbstTerm
      UIntra   = this%UIntra
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      EVirialIntra = this%EVirialIntra
      USelf    = this%USelbstTerm
      UIntra   = this%UIntra
      this%qgrida_old = this%qgrida
      call chargegrid_min ( this, nc,np )
      this%NPart = this%NPart - 1
      this%Component(nc)%NPart = this%Component(nc)%NPart - 1
      call PMESelfTermMC ( this )
! For further use of the following code
      this%NPart = this%NPart + 1
      this%Component(nc)%NPart = this%Component(nc)%NPart + 1
    end if

    ! Update number of delete attempts
    this%NDeleteAttempts = this%NDeleteAttempts + 1

    ! Calculate particle energy
#if MPI_VER > 0
    call MPI_Allreduce( GetEnergy( this, nc, np ), EPotDel, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
#else
    EPotDel = GetEnergy( this, nc, np )
#endif
    EPotDel = EPotDel + this%Density * pc%EPotTestCorrLJ &
&                     + pc%EPotTestCorrRF

    if (LongRange .eq. Ewald) then
    ! Save Coordinates
      do i=1,3,1
         r(i) = pc%P0(np,i)
      end do
      do i=1,4,1
         q(i) = pc%Q0(np,i)
      end do

      call RemoveParticle( pc, np )
      this%NPart = this%NPart - 1

!        call EwaldSelfTerm_Energy( this,nc,np )

       EPotNew = GetEnergy( this, nc, np )

       EPotDel = EPotDel + this%USelbstTerm + this%UIntra
    else if (LongRange .eq. PME) then
       EPotDel = EpotDel + (this%UIntra + this%USelbstTerm) - (USelf+UIntra)
    end if

!DEBUG
!  write(0, '(I2, ": EPotDel = ", F12.6)') NProc, EPotDel
!#if MPI_VER > 0
!  accepted = .FALSE.
!#endif
!DEBUG

    ! Apply acceptance criterion
    if( rnd( 0._RK, 1._RK ) .lt. &
&       ( exp( EPotDel / this%Temperature - pc%ChemPot ) &
&         * this%Density * pc%Fraction )) then

      ! Accept Deletion
      this%NDeleteSuccesses = this%NDeleteSuccesses + 1
      call RemoveParticle( pc, np )

      ! Copy energies and virial
      n1 = pc%NPart + 1
      do i = 1, this%NComponents
        pi => this%Interaction(nc, i)
        n2 = pi%NPart2
        pi%EPot(np, 1:n2) = pi%EPot(n1, 1:n2)
        pi%Virial(np, 1:n2) = pi%Virial(n1, 1:n2)
        this%Interaction(i, nc)%EPot(1:n2, np) = pi%EPot(n1, 1:n2)
        this%Interaction(i, nc)%Virial(1:n2, np) = pi%Virial(n1, 1:n2)
      end do

      ! Zero diagonal elements
      this%Interaction(nc, nc)%EPot(np, np) = 0._RK
      this%Interaction(nc, nc)%Virial(np, np) = 0._RK

      this%NPart = this%NPart - 1

      ! Update density
      this%Density = this%NPart / this%Volume0

      ! Update fractions and NDF
      call UpdateFractions( this )

      ! Update long range correction
      call CalculateCorr( this )

!DEBUG
!#if MPI_VER > 0
!  accepted = .TRUE.
!#endif
!DEBUG

    else 
      if (LongRange .eq. Ewald) then
        this%UFourier = EFourier
        call AddParticle(pc,r,q)
        this%NPart = this%NPart + 1
      else if (LongRange .eq. PME) then
        this%UFourier = EFourier
        this%EVirial  = EVirial
        this%qgrida   = this%qgrida_old
        call PMESelfTermMC (this)
      end if
    end if

!DEBUG
!#if MPI_VER > 0
!  call MPI_Allreduce( accepted, different, 1, MPI_LOGICAL, MPI_LXOR, &
!&   MPI_COMM_WORLD, ierror )
!  if( different ) then
!    write(0, '(I2, ": Delete of comp. ", I0, A, " at step ", I0)') &
!&     NProc, nc, merge("    accepted", "not accepted", accepted), step
!    write(0, '(I2, ": Next random number = ", F12.10)') NProc, rnd(0._RK, 1._RK)
!    stop
!  end if
!#endif
!DEBUG

  end subroutine TEnsemble_Delete



!==============================================================!
!  Subroutine TEnsemble_Move2End                               !
!==============================================================!

  subroutine TEnsemble_Move2End( this, nc, np )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)        :: this
    integer, intent(in)    :: nc
    integer, intent(inout) :: np

    ! Declare local variables
    type(TComponent), pointer   :: pc
    type(TInteraction), pointer :: pi
    integer                     :: i, n1, n2
    real(RK)                    :: PSave(3), QSave(4)
    real(RK)                    :: ESave(this%NPartMax), VSave(this%NPartMax)

    ! Assign local variables
    pc => this%Component(nc)
    n1 = pc%NPart

    ! Copy position and quaternions
    PSave(:) = pc%P0(np, :)
    pc%P0(np, :) = pc%P0(n1, :)
    pc%P0(n1, :) = PSave(:)
    if( pc%Molecule%IsElongated ) then
      QSave(:) = pc%Q0(np, :)
      pc%Q0(np, :) = pc%Q0(n1, :)
      pc%Q0(n1, :) = QSave(:)
    end if

    ! Convert molecular coordinates to atom positions
    call Mol2Atom1( pc, np )
    call Mol2Atom1( pc, n1 )

    ! Copy energies and virial
    do i = 1, this%NRealComponents
      pi => this%Interaction(nc, i)
      n2 = pi%NPart2
      ESave(1:n2) = pi%EPot(np, :)
      VSave(1:n2) = pi%Virial(np, :)
      if( i .eq. nc ) then
        ESave(np) = pi%EPot(np, n2)
        VSave(np) = pi%Virial(np, n2)
      end if
      pi%EPot(np, :) = pi%EPot(n1, :)
      pi%Virial(np, :) = pi%Virial(n1, :)
      this%Interaction(i, nc)%EPot(:, np) = pi%EPot(n1, :)
      this%Interaction(i, nc)%Virial(:, np) = pi%Virial(n1, :)
      pi%EPot(n1, :) = ESave(1:n2)
      pi%Virial(n1, :) = VSave(1:n2)
      this%Interaction(i, nc)%EPot(:, n1) = ESave(1:n2)
      this%Interaction(i, nc)%Virial(:, n1) = VSave(1:n2)
    end do

    ! Zero diagonal elements
    this%Interaction(nc, nc)%EPot(np, np) = 0._RK
    this%Interaction(nc, nc)%Virial(np, np) = 0._RK

    ! Set new particle number
    np = n1

  end subroutine TEnsemble_Move2End



!==============================================================!
!  Subroutine TEnsemble_Resize                                 !
!==============================================================!

  subroutine TEnsemble_Resize( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK) :: VolumeOld, EPotOld
    real(RK) :: EPotDelta
    real(RK) :: EVirial
    real(RK) :: UFourier
    real(RK) :: UIntra, EVirialintra
#if MPI_VER > 0
    real(RK) :: EPotNew
#endif

    ! Update number of resizing attempts
    this%NResizeAttempts = this%NResizeAttempts + 1

    ! Save current simulation box size, volume, energy, virial
    VolumeOld = this%Volume0
    EPotOld = this%EPot
    if (LongRange .eq. Ewald) then
       UIntra  = this%UIntra
       EVirialIntra = this%EVirialIntra
       UFourier= this%UFourier
       EVirial = this%EVirial
    else if (LongRange .eq. PME) then
       UIntra  = this%UIntra
       EVirialIntra = this%EVirialIntra
       UFourier= this%UFourier
       EVirial = this%EVirial
    end if

    ! Generate a trial volume change
    this%Volume0 = this%Volume0 * (1._RK + rnd( -this%DispVol, this%DispVol ))
    call UpdateBoxLength( this )

    ! Convert molecular coordinates to atom positions
    call Mol2Atom( this )

    ! Calculate potential energy and virial at trial position
#if MPI_VER > 0
    call Energy( this, EPotNew )
    call MPI_Allreduce( EPotNew, this%EPot, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
#else
    call Energy( this, this%EPot )
#endif

    ! Find potential change
    EPotDelta = this%RefPressure * (this%Volume0 - VolumeOld) &
&     + this%EPot - EPotOld &
&     + this%NPart * this%Temperature * log( VolumeOld / this%Volume0 )

    ! Apply Metropolis acceptance criterion
    if( exp( -EPotDelta / this%Temperature ) .gt. rnd( 0._RK, 1._RK ) ) then

      ! Accept volume change
      this%NResizeSuccesses = this%NResizeSuccesses + 1

      ! Update energy and virial matrices
      call UpdateEnergy( this )
#if MPI_VER > 0
      call MPI_Allreduce( GetVirial( this ), this%Virial, 1, &
&       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
#else
      this%Virial = GetVirial( this )
#endif

    else

      ! Reject volume change
      this%Volume0 = VolumeOld
      call UpdateBoxLength( this )
      call Mol2Atom( this )
      this%EPot = EPotOld
      if (LongRange .eq. Ewald) then
         this%UIntra = UIntra
!          this%EVirialIntra = EVirialIntra
         this%UFourier = UFourier
         call Energy(this,this%Epot)
!          this%EVirial = EVirial
#if MPI_VER > 0
         call MPI_Allreduce( GetEnergy( this ), this%EPot, 1 , &
&            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
         call MPI_Allreduce( GetVirial( this ), this%Virial, 1 , &
&            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
#else
         this%EPot = GetEnergy(this)
         this%Virial = GetVirial( this )
#endif
      else if (LongRange .eq. PME) then
         this%UIntra = UIntra
         this%EVirialIntra = EVirialIntra
         this%UFourier = UFourier
         this%EVirial = EVirial
         call charge_grid_MCall ( this )
      end if

    end if

  end subroutine TEnsemble_Resize



!==============================================================!
!  Subroutine TEnsemble_ZeroNAttempts                          !
!==============================================================!

  subroutine TEnsemble_ZeroNAttempts( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Zero number of MC attempts and successes
    do i = 1, this%NComponents
      call ZeroNAttempts( this%Component(i) )
    end do
    this%NResizeAttempts  = 0
    this%NResizeSuccesses = 0
    this%NInsertAttempts  = 0
    this%NInsertSuccesses = 0
    this%NDeleteAttempts  = 0
    this%NDeleteSuccesses = 0

  end subroutine TEnsemble_ZeroNAttempts



!==============================================================!
!  Subroutine TEnsemble_UpdateDisplacements                    !
!==============================================================!

  subroutine TEnsemble_UpdateDisplacements( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    type(TComponent), pointer :: pc
    integer                   :: i
    real(RK) :: AccRateTran, AccRateRot, AccRateVol

    ! Calculate acceptance rates
    do i = 1, this%NRealComponents
      pc => this%Component(i)
      AccRateTran = real(pc%NMoveSuccesses) / real(pc%NMoveAttempts)
      AccRateRot = real(pc%NRotateSuccesses) / real(pc%NRotateAttempts)

      ! Update translational displacement
      if(( AccRateTran .gt. AccUpperLimit) .and. &
&        ( pc%DispTran .lt. DispTranLimit )) then
        pc%DispTran = pc%DispTran * 1.05_RK
      else if( AccRateTran .lt. AccLowerLimit ) then
        pc%DispTran = pc%DispTran * .95_RK
      end if

      ! Update rotational displacement
      if(( AccRateRot .gt. AccUpperLimit ) .and. &
&       ( pc%DispRot .lt. DispRotLimit )) then
        pc%DispRot = pc%DispRot * 1.05_RK
      else if( AccRateRot .lt. AccLowerLimit ) then
        pc%DispRot = pc%DispRot * 0.95_RK
      end if
    end do

    if( ConstantPressure .and. .not. NVTEquilibration ) then
      AccRateVol = real(this%NResizeSuccesses) / real(this%NResizeAttempts)

      ! Update volume displacement
      if(( AccRateVol .gt. AccUpperLimit ) .and. &
&        ( this%DispVol .lt. DispVolLimit )) then
        this%DispVol = this%DispVol * 1.05_RK
      else if( AccRateVol .lt. AccLowerLimit ) then
        this%DispVol = this%DispVol * 0.95_RK
      end if
    end if

    ! Zero attempts
    call ZeroNAttempts( this )

  end subroutine TEnsemble_UpdateDisplacements



!==============================================================!
!  Subroutine TEnsemble_SaveState                              !
!==============================================================!

  subroutine TEnsemble_SaveState( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Save current state
    do i = 1, this%NRealComponents
      call SaveState( this%Component(i) )
    end do

    ! Save current random number
!     ixsave = ix
!     iysave = iy

  end subroutine TEnsemble_SaveState



!==============================================================!
!  Subroutine TEnsemble_RestoreState                           !
!==============================================================!

  subroutine TEnsemble_RestoreState( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Restore current state
    do i = 1, this%NRealComponents
      call RestoreState( this%Component(i) )
    end do

    call Energy( this, this%EPot )
    call UpdateEnergy( this )

    ! Restore current random number
!     ix = ixsave
!     iy = iysave

  end subroutine TEnsemble_RestoreState



!==============================================================!
!  Subroutine TEnsemble_ResultOpen                             !
!==============================================================!

  subroutine TEnsemble_ResultOpen( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    if( Restart ) then
      ! Open result file
      write( IOBuffer, '(I16)' ) this%EnsembleNumber
      call FileAppend( this%iounit_result, &
&       trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ResultFileExtension )

      if( .not. SimulationType .eq. SecondVirialCoeff ) then
        ! Open running average result file
        write( IOBuffer, '(I16)' ) this%EnsembleNumber
        call FileAppend( this%iounit_runave, &
&         trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//RunAveFileExtension )
      end if

    else
      ! Open result file
      write( IOBuffer, '(I16)' ) this%EnsembleNumber
      call FileRewrite( this%iounit_result, &
&       trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ResultFileExtension )

      if( .not. SimulationType .eq. SecondVirialCoeff ) then
        ! Open running average result file
        write( IOBuffer, '(I16)' ) this%EnsembleNumber
        call FileRewrite( this%iounit_runave, &
&         trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//RunAveFileExtension )
      end if

    end if

  end subroutine TEnsemble_ResultOpen



!==============================================================!
!  Subroutine TEnsemble_ResultUpdate                           !
!==============================================================!

  subroutine TEnsemble_ResultUpdate( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    type(TComponent), pointer :: pc
    integer                   :: i,err
    real(RK)                  :: value

    if( Step == 1 ) then
      ! Reset accumulators
      ! 1.) Basic sums
      call Reset( this%SumPressure )
      call Reset( this%SumDensity )
      call Reset( this%SumTemperature )
      call Reset( this%SumEPot )
      call Reset( this%SumEnthalpy )
      call Reset( this%SumVolume )
      call Reset( this%SumVirial )
      if( EnsembleType .eq. EnsembleTypeGE .or. &
&         EnsembleType .eq. EnsembleTypeHA ) then
        call Reset( this%SumNPart )
        do i = 1, this%NComponents
          call Reset( this%Component(i)%SumFraction )
        end do
      end if

      ! 2.) Combined sums
      call Reset( this%SumEPotSquared )
      call Reset( this%SumEPotV )
      call Reset( this%SumEPotVirial )
      call Reset( this%SumEnthalpySquared )
      call Reset( this%SumEnthalpyV )
      call Reset( this%SumVolumeSquared )

      ! 3.) Derived sums
      if( ConstantPressure ) then
        call Reset( this%SumBetaT )
        call Reset( this%SumdHdP )
        call Reset( this%SumCP )
        call Reset( this%SumAlphaP )
      else
        call Reset( this%SumdUdV )
        call Reset( this%SumCV )
      endif

      ! 4.) Chemical potential and partial molar volumes
      do i = 1, this%NRealComponents
        select case( this%Component(i)%ChemPotMethod )
        case( ChemPotMethodGradIns )
          call Reset( this%Component(i)%SumInvChemPotRho )
          call Reset( this%Component(i)%SumInvChemPot )
!DEBUG
          call Reset( this%Component(i)%SumInvChemPotRho1 )
          call Reset( this%Component(i)%SumInvChemPot1 )
          call Reset( this%Component(i)%SumInvChemPotRho2 )
          call Reset( this%Component(i)%SumInvChemPot2 )
!DEBUG
        case( ChemPotMethodWidom )
          call Reset( this%Component(i)%SumChemPotV )
          call Reset( this%Component(i)%SumChemPotVV )
        end select
      end do
!       if( ConstantPressure .and. this%NRealComponents > 1 ) then
        do i = 1, this%NRealComponents
          if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodNone ) then
            call Reset( this%Component(i)%SumVW )
!DEBUG
            if( this%Component(i)%ChemPotMethod .eq. ChemPotMethodGradIns ) then
              call Reset( this%Component(i)%SumVW1 )
              call Reset( this%Component(i)%SumVW2 )
            end if
!DEBUG
          end if
        end do
!       end if

      ! Update result header
      call FileWriteBlank( this%iounit_result )
      call FileWriteBlank( this%iounit_runave )

      ! Number of steps
      write( IOBuffer, '("     NR")' )
      call FileWriteNoAdvance( this%iounit_result )
      call FileWriteNoAdvance( this%iounit_runave )

      ! Displacement
      if( SimulationType .eq. MolecularDynamics ) then
        write( IOBuffer, '("    DISP")' )
        call FileWriteNoAdvance( this%iounit_runave )
      end if

      ! Pressure
      write( IOBuffer, '("     DRUCK")' )
      call FileWriteNoAdvance( this%iounit_result )
      call FileWriteNoAdvance( this%iounit_runave )

      ! Density
      write( IOBuffer, '("    DICHTE")' )
      call FileWriteNoAdvance( this%iounit_result )
      call FileWriteNoAdvance( this%iounit_runave )

      ! Temperature
      write( IOBuffer, '("      TEMP")' )
      call FileWriteNoAdvance( this%iounit_result )
      call FileWriteNoAdvance( this%iounit_runave )

      ! Potential energy
      write( IOBuffer, '("      EPOT")' )
      call FileWriteNoAdvance( this%iounit_result )
      call FileWriteNoAdvance( this%iounit_runave )

      ! Enthalpy
      write( IOBuffer, '("     ENTLP")' )
      call FileWriteNoAdvance( this%iounit_result )
      call FileWriteNoAdvance( this%iounit_runave )

      ! Chemical potential
      do i = 1, this%NRealComponents
        if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodNone ) then
          write( IOBuffer, '("     MUE", I2)' ) i
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )
        end if
      end do

      ! Partial molar volume
!       if( ConstantPressure .and. this%NRealComponents > 1 ) then
        do i = 1, this%NRealComponents
          if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodNone ) then
            write( IOBuffer, '("      VW", I2)' ) i
            call FileWriteNoAdvance( this%iounit_result )
            call FileWriteNoAdvance( this%iounit_runave )
          end if
        end do
!       end if

      ! Number of particles in ensemble
      if( EnsembleType .eq. EnsembleTypeGE .or. &
&         EnsembleType .eq. EnsembleTypeHA ) then
        write( IOBuffer, '("     NPART")' )
        call FileWriteNoAdvance( this%iounit_result )
        call FileWriteNoAdvance( this%iounit_runave )

        ! Molar fraction of each component
        do i = 1, this%NComponents
          write( IOBuffer, '("   FRACT", I2)' ) i
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )
        end do
      end if

#if CONSTR > 0
      do i=1, this%NCons
        write( IOBuffer, '("     PMF", I2)' ) i
        call FileWriteNoAdvance( this%iounit_runave )
        write( IOBuffer, '("     MF",  I2)' ) i
        call FileWriteNoAdvance( this%iounit_runave )
      end do
#endif

      call FileWriteBlank( this%iounit_result )
      call FileWriteBlank( this%iounit_runave )
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END IF of step ==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Update accumulators
    ! 1.) Basic sums
    call Update( this%SumPressure, this%Pressure )
    call Update( this%SumDensity, this%Density )
    call Update( this%SumTemperature, this%Temperature )
    call Update( this%SumEPot, this%EPot / real( this%NPart, RK ) )
    if( ConstantPressure ) then
      call Update( this%SumEnthalpy, this%EPot / real( this%NPart, RK ) + &
&       this%RefPressure / this%Density - this%RefTemperature )
    else
      call Update( this%SumEnthalpy, this%EPot / real( this%NPart, RK ) + &
&       this%Pressure / this%Density - this%RefTemperature )
    end if
    call Update( this%SumVolume, 1._RK / this%Density )
    call Update( this%SumVirial, -3._RK * this%Virial )
    if( EnsembleType .eq. EnsembleTypeGE .or. &
&       EnsembleType .eq. EnsembleTypeHA ) then
      call Update( this%SumNPart, real( this%NPart, RK ) )
      do i = 1, this%NComponents
        pc => this%Component(i)
        call Update( pc%SumFraction, pc%Fraction )
      end do
    end if

    ! 2.) Combined sums
    call Update( this%SumEPotSquared, &
&     ( this%EPot / real( this%NPart, RK ) )**2 )
    call Update( this%SumEPotV, &
&     this%EPot / ( real( this%NPart, RK ) * this%Density ) )
    call Update( this%SumEPotVirial, &
&     -3. * this%Virial * this%EPot / real( this%NPart, RK ) )
    call Update( this%SumEnthalpySquared, &
&     ( this%EPot / real( this%NPart, RK ) + &
&     this%RefPressure / this%Density - this%RefTemperature )**2 )
    call Update( this%SumEnthalpyV, &
&     ( this%EPot / real( this%NPart, RK ) + &
&     this%RefPressure / this%Density - this%RefTemperature ) / this%Density )
    call Update( this%SumVolumeSquared, 1._RK / this%Density**2 )

    ! 3.) Derived sums
    if( ConstantPressure ) then
      call Update( this%SumBetaT, &
&       real( this%NPart, RK ) / this%RefTemperature &
&       * ( this%SumVolumeSquared%Average / this%SumVolume%Average &
&         - this%SumVolume%Average ) )
      call Update( this%SumdHdP, this%SumVolume%Average &
&       - real( this%NPart, RK ) / this%RefTemperature &
&       * ( this%SumEPotV%Average - this%SumEPot%Average * this%SumVolume%Average &
&         + this%RefPressure &
&         * ( this%SumVolumeSquared%Average - this%SumVolume%Average**2 ) ) )
      call Update( this%SumCP, real( this%NPart, RK ) / this%RefTemperature**2 &
&       * ( this%SumEnthalpySquared%Average - this%SumEnthalpy%Average**2 ) )
      call Update( this%SumAlphaP, real( this%NPart, RK ) / this%RefTemperature**2 &
&       * this%SumDensity%Average * ( this%SumEnthalpyV%Average &
&         - this%SumEnthalpy%Average * this%SumVolume%Average ) )
    else
      call Update( this%SumdUdV, this%Density / ( 3. * real( this%NPart, RK ) ) &
&       * ( this%NPart / this%RefTemperature &
&       * ( this%SumVirial%Average * this%SumEPot%Average &
&         - this%SumEPotVirial%Average ) + this%SumVirial%Average ) )
      call Update( this%SumCV, real( this%NPart, RK ) / this%RefTemperature**2 &
&       * ( this%SumEPotSquared%Average - this%SumEPot%Average**2 ) )
    endif

    ! 4.) Chemical potential and partial molar volumes
    do i = 1, this%NRealComponents
      pc => this%Component(i)
      if( pc%CalcChemPot ) then
        select case( pc%ChemPotMethod )
        case( ChemPotMethodGradins )
          call Update( pc%SumInvChemPotRho, 1._RK / pc%ChemPot )
          call Update( pc%SumInvChemPot, 1._RK / pc%ChemPot1 )
!DEBUG
          call Update( pc%SumInvChemPotRho1, this%Density / pc%ChemPot1 )
          call Update( pc%SumInvChemPot1, 1._RK / pc%ChemPot1 )
          call Update( pc%SumInvChemPotRho2, 1._RK / pc%ChemPot2 )
          call Update( pc%SumInvChemPot2, 1._RK / pc%ChemPot1 )
!DEBUG
        case( ChemPotMethodWidom )
          call Update( pc%SumChemPotV, pc%ChemPot / this%Density )
          call Update( pc%SumChemPotVV, pc%ChemPot / this%Density**2 )
        end select
      end if
    end do
!     if( ConstantPressure .and. this%NRealComponents > 1 ) then
      do i = 1, this%NRealComponents
        pc => this%Component(i)
        if( pc%CalcChemPot ) then
          select case( pc%ChemPotMethod )
          case( ChemPotMethodGradIns )
            call Update( pc%SumVW, this%NPart &
&             * ( this%SumVolume%Average &
&               - pc%SumInvChemPot%Average / pc%SumInvChemPotRho%Average ) )
!DEBUG
          call Update( pc%SumVW1, this%NPart &
&             * ( this%SumVolume%Average &
&               - pc%SumInvChemPot1%Average / pc%SumInvChemPotRho1%Average ) )
          call Update( pc%SumVW2, this%NPart &
&             * ( this%SumVolume%Average &
&               - pc%SumInvChemPot2%Average / pc%SumInvChemPotRho2%Average ) )
!DEBUG
          case( ChemPotMethodWidom )
            call Update( pc%SumVW, this%NPart &
&             * ( pc%SumChemPotVV%Average / pc%SumChemPotV%Average &
&               - this%SumVolume%Average ) )
          end select
        end if
      end do
!     end if

    ! Update result files
    if( mod( Step, BlockSize ) == 0 ) then

      ! Number of steps
      write( IOBuffer, '(I7)' ) Step
      call FileWriteNoAdvance( this%iounit_result )
      call FileWriteNoAdvance( this%iounit_runave )

      ! Displacement
      if( SimulationType .eq. MolecularDynamics ) then
        value = 0._RK
        do i = 1, this%NComponents
          value = value + sum( this%Component(i)%Disp(:, :)**2 )
        end do
        value = value * this%BoxLength**2 &
&               / ( 6._RK * this%NPart * TimeStep * Step )
        write( IOBuffer, '(F8.3)' ) value
        call FileWriteNoAdvance( this%iounit_runave )
      end if

      ! Pressure
      write( IOBuffer, '(F10.5)' ) this%SumPressure%BlockAverage
      call FileWriteNoAdvance( this%iounit_result )
      write( IOBuffer, '(F10.5)' ) this%SumPressure%Average
      call FileWriteNoAdvance( this%iounit_runave )

      ! Density
      write( IOBuffer, '(F10.5)' ) this%SumDensity%BlockAverage
      call FileWriteNoAdvance( this%iounit_result )
      write( IOBuffer, '(F10.5)' ) this%SumDensity%Average
      call FileWriteNoAdvance( this%iounit_runave )

      ! Temperature
      write( IOBuffer, '(F10.5)' ) this%SumTemperature%BlockAverage
      call FileWriteNoAdvance( this%iounit_result )
      write( IOBuffer, '(F10.5)' ) this%SumTemperature%Average
      call FileWriteNoAdvance( this%iounit_runave )

      ! Potential energy
      write( IOBuffer, '(F10.5)' ) this%SumEPot%BlockAverage
      call FileWriteNoAdvance( this%iounit_result )
      write( IOBuffer, '(F10.5)' ) this%SumEPot%Average
      call FileWriteNoAdvance( this%iounit_runave )

      ! Enthalpy
      write( IOBuffer, '(F10.5)' ) this%SumEnthalpy%BlockAverage
      call FileWriteNoAdvance( this%iounit_result )
      write( IOBuffer, '(F10.5)' ) this%SumEnthalpy%Average
      call FileWriteNoAdvance( this%iounit_runave )

      ! Chemical potential
      do i = 1, this%NRealComponents
        pc => this%Component(i)
        if( pc%ChemPotMethod .ne. ChemPotMethodNone ) then
          if( Equilibration ) then
            write( IOBuffer, '(F10.5)' ) 0._RK
            call FileWriteNoAdvance( this%iounit_result )
            call FileWriteNoAdvance( this%iounit_runave )
          else
            if( pc%Fraction > 0._RK ) then
              select case( pc%ChemPotMethod )
              case( ChemPotMethodGradIns )
                write( IOBuffer, '(F10.5)' ) &
&                 log( pc%Fraction * pc%SumInvChemPotRho%BlockAverage )
                call FileWriteNoAdvance( this%iounit_result )
                write( IOBuffer, '(F10.5)' ) &
&                 log( pc%Fraction * pc%SumInvChemPotRho%Average )
                call FileWriteNoAdvance( this%iounit_runave )
              case( ChemPotMethodWidom )
                write( IOBuffer, '(F10.5)' ) &
&                 log( pc%Fraction / pc%SumChemPotV%BlockAverage )
                call FileWriteNoAdvance( this%iounit_result )
                write( IOBuffer, '(F10.5)' ) &
&                 log( pc%Fraction / pc%SumChemPotV%Average )
                call FileWriteNoAdvance( this%iounit_runave )
              end select
            else
              write( IOBuffer, '(F10.5)' ) &
&               log( 1._RK / pc%SumChemPotV%BlockAverage )
              call FileWriteNoAdvance( this%iounit_result )
              write( IOBuffer, '(F10.5)' ) &
&               log( 1._RK / pc%SumChemPotV%Average )
              call FileWriteNoAdvance( this%iounit_runave )
            end if
          end if
        end if
      end do

      ! Partial molar volume
!       if( ConstantPressure .and. this%NRealComponents > 1 ) then
        do i = 1, this%NRealComponents
          pc => this%Component(i)
          if( pc%ChemPotMethod .ne. ChemPotMethodNone ) then
            if( Equilibration ) then
              write( IOBuffer, '(F10.5)' ) 0._RK
              call FileWriteNoAdvance( this%iounit_result )
              call FileWriteNoAdvance( this%iounit_runave )
            else
              write( IOBuffer, '(F10.5)' ) pc%SumVW%BlockAverage
              call FileWriteNoAdvance( this%iounit_result )
              write( IOBuffer, '(F10.5)' ) pc%SumVW%Average
              call FileWriteNoAdvance( this%iounit_runave )
            end if
          end if
        end do
!       end if

      ! Number of particles in ensemble
      if( EnsembleType .eq. EnsembleTypeGE .or. &
&         EnsembleType .eq. EnsembleTypeHA ) then
        write( IOBuffer, '(F10.2)' ) this%SumNPart%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(F10.2)' ) this%SumNPart%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! Molar fraction of each component
        do i = 1, this%NComponents
          pc => this%Component(i)
          write( IOBuffer, '(F10.5)' ) pc%SumFraction%BlockAverage
          call FileWriteNoAdvance( this%iounit_result )
          write( IOBuffer, '(F10.5)' ) pc%SumFraction%Average
          call FileWriteNoAdvance( this%iounit_runave )
        end do
      end if

      call FileWriteBlank( this%iounit_result )
#if CONSTR == 0
      call FileWriteBlank( this%iounit_runave )
#else
    this%consup = .true.
#endif
#if ARCH == 2
      call flush( this%iounit_result )
      call flush( this%iounit_runave )
#endif

    end if

! Exit, if specific file is in the folder!
    open( 99, file = 'stop.txt', action = 'READ', status = 'OLD', &
&     iostat = err )
#ifdef __INTEL_COMPILER
    if( err .eq. 0 ) err = SetTerminateProgram( 1 ) 
#else
    if( err .eq. 0 ) call SetTerminateProgram
#endif


  end subroutine TEnsemble_ResultUpdate



!==============================================================!
!  Subroutine TEnsemble_ResultClose                            !
!==============================================================!

  subroutine TEnsemble_ResultClose( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Close running average result file
    if( .not. SimulationType .eq. SecondVirialCoeff ) &
&     call FileClose( this%iounit_runave )

    ! Close result file
    call FileClose( this%iounit_result )

  end subroutine TEnsemble_ResultClose



!==============================================================!
!  Subroutine TEnsemble_ErrorsUpdate                           !
!==============================================================!

  subroutine TEnsemble_ErrorsUpdate( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK)                  :: Average, Variance
    type(TComponent), pointer :: pc
    integer                   :: i, j

    ! Declare local variables for velocity of sound
    real(RK) :: molmass, cpid

    ! Declare local variables for phase equilibrium results
    real(RK) :: NN, yvi
    real(RK) :: dpdmu( this%NComponents ), dpdv( this%NComponents )
    real(RK) :: dydmu( this%NComponents, this%NComponents ), &
&               dydv( this%NComponents, this%NComponents )
    real(RK) :: varmu( this%NComponents ), varv( this%NComponents )
    real(RK) :: vary( this%NComponents - 1 )
    real(RK) :: VarPressure, DeltaHv, VarDeltaHv
#ifdef ABL
    integer  :: counter
#endif

    ! Calculate averages and errors
    call Error( this%SumPressure )
    call Error( this%SumDensity )
    call Error( this%SumTemperature )
    call Error( this%SumEPot )
    call Error( this%SumEnthalpy )
    call Error( this%SumVolume )
    if( EnsembleType .eq. EnsembleTypeGE .or. &
&       EnsembleType .eq. EnsembleTypeHA ) then
      do i = 1, this%NComponents
        pc => this%Component(i)
        call Error( pc%SumFraction )
      end do
    else
      if( ConstantPressure ) then
        call Error( this%SumBetaT )
        call Error( this%SumdHdP )
        call Error( this%SumCP )
        call Error( this%SumAlphaP )
      else
        call Error( this%SumdUdV )
        call Error( this%SumCV )
      end if
      do i = 1, this%NRealComponents
        pc => this%Component(i)
        if( pc%CalcChemPot ) then
          select case( pc%ChemPotMethod )
          case( ChemPotMethodGradIns )
            call Error( pc%SumInvChemPotRho )
!DEBUG
            call Error( pc%SumInvChemPotRho1 )
            call Error( pc%SumInvChemPotRho2 )
!DEBUG
          case( ChemPotMethodWidom )
            call Error( pc%SumChemPotV )
          end select
!           if( ConstantPressure .and. this%NRealComponents > 1 ) &
! &           call Error( pc%SumVW )
            call Error( pc%SumVW )
!DEBUG
            if( pc%ChemPotMethod .eq. ChemPotMethodGradIns ) then
              call Error( pc%SumVW1 )
              call Error( pc%SumVW2 )
            end if
!DEBUG
        end if
      end do
    end if

    ! Open final result file
    write( IOBuffer, '(I16)' ) this%EnsembleNumber
    call FileRewrite( this%iounit_errors, &
&     trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ErrorsFileExtension )

    ! Separator
    write( IOBuffer, '(76("="))' )
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )
    write( IOBuffer, '(T24, "SIMULATION RESULT FILE")' )
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T24, "----------------------")' )
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Simulation type
    write( IOBuffer, '("Simulation type", T36, ":", 9X, A)' ) &
&     trim( SimulationTypeString )
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Ensemble type", T36, ":", 9X, A)' ) &
&     trim( EnsembleTypeString )
    call FileWrite( this%iounit_errors )
    if( SimulationType .eq. MolecularDynamics ) then
      write( IOBuffer, '("Integrator type", T36, ":", 9X, A)' ) &
&       trim( IntegratorTypeString )
      call FileWrite( this%iounit_errors )
    end if
    call FileWriteBlank( this%iounit_errors )

    ! Number of steps
    write( IOBuffer, '("Number of NVT equilibration steps", T36, ":", I10)' ) &
&     NStepsV
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Number of NPT equilibration steps", T36, ":", I10)' ) &
&     NStepsP
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Number of production steps", T36, ":", I10)' ) &
&     Step
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Time step
    if( SimulationType .eq. MolecularDynamics ) then
      write( IOBuffer, '("Time step", T29, "reduced:", F20.9)' ) &
&       TimeStep
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T31, "in fs:", F20.9)' ) &
&       TimeStep * UnitTime * 1E15_RK
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )
    end if

    ! Acceptance rate
    if( SimulationType .eq. MonteCarlo ) then
      write( IOBuffer, '("Acceptance rate", T36, ":", F20.9)' ) &
&       Acceptance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )
    end if

    ! Mass of piston
    if( SimulationType .eq. MolecularDynamics .and. ConstantPressure ) then
      write( IOBuffer, '("Mass of piston", T36, ":", F20.9)' ) &
&       this%PistonMass
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )
    end if

    ! Number of particles
    write( IOBuffer, '("Number of particles", T36, ":", I10)' ) &
&     this%NPart
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Potential models
    if( EnsembleType .ne. EnsembleTypeGE .or. &
&       EnsembleType .ne. EnsembleTypeHA ) then
      do i = 1, this%NRealComponents
        write( IOBuffer, '("Molar fraction of ", A, T36, ":", F20.9)' ) &
&         trim( this%Component(i)%Molecule%PotModFileName ), &
&         this%Component(i)%Fraction
        call FileWrite( this%iounit_errors )
        select case( this%Component(i)%ChemPotMethod )
        case( ChemPotMethodGradIns )
          write( IOBuffer, '("Chemical potential calculated by gradual insertion")' )
          call FileWrite( this%iounit_errors )
        case( ChemPotMethodWidom )
          write( IOBuffer, '("Number of test particles", T36, ":", I10)' ) &
&           this%Component(i)%NTestAll
          call FileWrite( this%iounit_errors )
        end select
      end do
      call FileWriteBlank( this%iounit_errors )
    end if

    ! Initial pressure
    if( ConstantPressure ) then
      write( IOBuffer, '("Initial pressure", T29, "reduced:", F20.9)' ) &
&       this%RefPressure
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T30, "in MPa:", F20.9)' ) &
&       this%RefPressure * UnitPressure * 1E-6_RK
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )
    end if

    ! Initial density
    write( IOBuffer, '("Initial density", T29, "reduced:", F20.9)' ) &
&     this%RefDensity
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T28, "in mol/l:", F20.9)' ) &
&     this%RefDensity * UnitDensity
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Initial temperature
    write( IOBuffer, '("Initial temperature", T29, "reduced:", F20.9)' ) &
&     this%RefTemperature
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T32, "in K:", F20.9)' ) &
&     this%RefTemperature * UnitTemperature
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! System of units
    write( IOBuffer, '("Unit of length", T36, ":", F20.9, " A")' ) &
&     UnitLength / Angstroem
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Unit of energy", T36, ":", F20.9, " K")' ) &
&     UnitEnergy / kBoltzmann
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Unit of mass", T36, ":", F20.9, " a.u.")' ) &
&     UnitMass * NAvogadro * 1000._RK
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Cutoff radii
    if( this%NLJ126Max > 0 ) then
      write( IOBuffer, &
&       '("Lennard-Jones cutoff radius", T36, ":", F20.9, " sigma")' ) &
&       this%RCutoffLJ126LJ126
      call FileWrite( this%iounit_errors )
    end if
    if( this%NDipoleMax > 0 ) then
      write( IOBuffer, &
&       '("Dipole-dipole cutoff radius", T36, ":", F20.9, " A")' ) &
&       this%RCutoffDipoleDipole * UnitLength / Angstroem
      call FileWrite( this%iounit_errors )
      if( this%NQuadrupoleMax > 0 ) then
        write( IOBuffer, &
&         '("Dipole-quadrupole cutoff radius", T36, ":", F20.9, " A")' ) &
&         this%RCutoffDipoleQuadrupole * UnitLength / Angstroem
        call FileWrite( this%iounit_errors )
      end if
    end if
    if( this%NQuadrupoleMax > 0 ) then
      write( IOBuffer, &
&       '("Quadrupole-quadrupole cutoff radius", T36, ":", F20.9, " A")' ) &
&       this%RCutoffQuadrupoleQuadrupole * UnitLength / Angstroem
      call FileWrite( this%iounit_errors )
    end if
    call FileWriteBlank( this%iounit_errors )

    ! Dielectric constant
    if( this%NDipoleMax > 0 ) then
      write( IOBuffer, &
&       '("Dielectric constant:", F36.9)' ) &
&       this%RFEpsilon
      call FileWrite( this%iounit_errors )
    end if
    call FileWriteBlank( this%iounit_errors )

    ! Separator
    write( IOBuffer, '(76("="))' )
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )
    write( IOBuffer, '("VALUE", T31, "UNITS", T46, "AVERAGE", T66, "ERROR")' )
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("-----", T31, "-----", T46, "-------", T66, "-----")' )
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Pressure
    Average = this%SumPressure%Average
    Variance = this%SumPressure%Variance
    write( IOBuffer, '("Pressure", T29, "reduced:", 2F20.9)' ) &
&     Average, Variance
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T30, "in MPa:", 2F20.9)' ) &
&     Average * UnitPressure * 1E-6_RK, Variance * UnitPressure * 1E-6_RK
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Density
    Average = this%SumDensity%Average
    Variance = this%SumDensity%Variance
    write( IOBuffer, '("Density", T29, "reduced:", 2F20.9)' ) &
&     Average, Variance
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T28, "in mol/l:", 2F20.9)' ) &
&     Average * UnitDensity, Variance * UnitDensity
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Temperature
    Average = this%SumTemperature%Average
    Variance = this%SumTemperature%Variance
    write( IOBuffer, '("Temperature", T29, "reduced:", 2F20.9)' ) &
&     Average, Variance
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T32, "in K:", 2F20.9)' ) &
&     Average * UnitTemperature, Variance * UnitTemperature
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Potential energy
    Average = this%SumEPot%Average
    Variance = this%SumEPot%Variance
    write( IOBuffer, '("Potential energy", T29, "reduced:", 2F20.9)' ) &
&     Average, Variance
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) &
&     Average * UnitEnergy * NAvogadro, &
&     Variance * UnitEnergy * NAvogadro
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Enthalpy
    Average = this%SumEnthalpy%Average
    Variance = this%SumEnthalpy%Variance
    write( IOBuffer, '("Enthalpy", T29, "reduced:", 2F20.9)' ) &
&     Average, Variance
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) &
&     Average * UnitEnergy * NAvogadro, &
&     Variance * UnitEnergy * NAvogadro
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    if( EnsembleType .eq. EnsembleTypeGE .or. &
&       EnsembleType .eq. EnsembleTypeHA ) then
      ! Molar fraction
      do i = 1, this%NComponents
        pc => this%Component(i)
        Average = pc%SumFraction%Average
        Variance = pc%SumFraction%Variance
        write( IOBuffer, &
&         '("Molar fraction of ", A, T36, ":", 2F20.9)' ) &
&         trim( this%Component(i)%Molecule%PotModFileName ), Average, Variance
        call FileWrite( this%iounit_errors )
      end do
      call FileWriteBlank( this%iounit_errors )

    else
      ! Chemical potential
      do i = 1, this%NRealComponents
        pc => this%Component(i)
        select case( pc%ChemPotMethod )
        case( ChemPotMethodGradIns )
          Variance = pc%SumInvChemPotRho%Variance / pc%SumInvChemPotRho%Average
          Average = log( pc%Fraction * pc%SumInvChemPotRho%Average )
          write( IOBuffer, &
&           '("Chemical potential of ", A, T33, "r`d:", 2F20.9)' ) &
&           trim( this%Component(i)%Molecule%PotModFileName ), &
&           Average, Variance
          call FileWrite( this%iounit_errors )
!DEBUG
          Variance = pc%SumInvChemPotRho1%Variance / pc%SumInvChemPotRho1%Average
          Average = log( pc%Fraction * pc%SumInvChemPotRho1%Average )
          write( IOBuffer, &
&           '("Chemical potential 0 of ", A, T33, "r`d:", 2F20.9)' ) &
&           trim( this%Component(i)%Molecule%PotModFileName ), &
&           Average, Variance
          call FileWrite( this%iounit_errors )
          Variance = pc%SumInvChemPotRho2%Variance / pc%SumInvChemPotRho2%Average
          Average = log( pc%Fraction * pc%SumInvChemPotRho2%Average )
          write( IOBuffer, &
&           '("Chemical potential 1 of ", A, T33, "r`d:", 2F20.9)' ) &
&           trim( this%Component(i)%Molecule%PotModFileName ), &
&           Average, Variance
          call FileWrite( this%iounit_errors )
!DEBUG
        case( ChemPotMethodWidom )
          Variance = pc%SumChemPotV%Variance / pc%SumChemPotV%Average
          if( pc%Fraction > 0.0_RK ) then
            Average = log( pc%Fraction / pc%SumChemPotV%Average )
            write( IOBuffer, &
&             '("Chemical potential of ", A, T33, "r`d:", 2F20.9)' ) &
&             trim( this%Component(i)%Molecule%PotModFileName ), &
&             Average, Variance
          else
            Average = -log( pc%SumChemPotV%Average )
            write( IOBuffer, &
&             '("Chem. pot. at inf. dilution of ", A, T33, "r`d:", 2F20.9)' ) &
&             trim( this%Component(i)%Molecule%PotModFileName ), &
&             Average, Variance
            Average = this%Temperature / pc%SumChemPotV%Average
            write( IOBuffer, &
&             '("Henrys law constant of ", A, T33, "r`d:", 2F20.9)' ) &
&             trim( pc%Molecule%PotModFileName ), Average, Variance
            call FileWrite( this%iounit_errors )
            write( IOBuffer, '(T30, "in MPa:", 2F20.9)' ) &
&             Average * UnitPressure * 1E-6_RK, Variance * UnitPressure * 1E-6_RK
          end if
          call FileWrite( this%iounit_errors )
        end select
      end do
      if( any(this%Component(:)%ChemPotMethod .ne. ChemPotMethodNone) ) &
&       call FileWriteBlank( this%iounit_errors )

      ! Partial molar volume
!     if( ConstantPressure .and. this%NRealComponents > 1 ) then
      do i = 1, this%NRealComponents
        pc => this%Component(i)
        if( pc%ChemPotMethod .ne. ChemPotMethodNone ) then
          Average = pc%SumVW%Average
          Variance = pc%SumVW%Variance
          write( IOBuffer, &
&           '("Partial molar volume of ", A, T33, "r`d:", 2F20.9)' ) &
&           trim( this%Component(i)%Molecule%PotModFileName ), Average, Variance
          call FileWrite( this%iounit_errors )
          write( IOBuffer, &
&           '(T28, "in l/mol:", 2F20.9)' ) &
&           Average / UnitDensity, Variance / UnitDensity
          call FileWrite( this%iounit_errors )
!DEBUG
          if( pc%ChemPotMethod .eq. ChemPotMethodGradIns ) then
            Average = pc%SumVW1%Average
            Variance = pc%SumVW1%Variance
            write( IOBuffer, &
&             '("Partial molar volume0 of ", A, T33, "r`d:", 2F20.9)' ) &
&             trim( this%Component(i)%Molecule%PotModFileName ), Average, Variance
            call FileWrite( this%iounit_errors )
            write( IOBuffer, &
&             '(T28, "in l/mol:", 2F20.9)' ) &
&             Average / UnitDensity, Variance / UnitDensity
            call FileWrite( this%iounit_errors )
            Average = pc%SumVW2%Average
            Variance = pc%SumVW2%Variance
            write( IOBuffer, &
&             '("Partial molar volume1 of ", A, T33, "r`d:", 2F20.9)' ) &
&             trim( this%Component(i)%Molecule%PotModFileName ), Average, Variance
            call FileWrite( this%iounit_errors )
            write( IOBuffer, &
&             '(T28, "in l/mol:", 2F20.9)' ) &
&             Average / UnitDensity, Variance / UnitDensity
            call FileWrite( this%iounit_errors )
          end if
!DEBUG
        end if
      end do
      if( any(this%Component(:)%ChemPotMethod .ne. ChemPotMethodNone) ) &
&       call FileWriteBlank( this%iounit_errors )
!     end if

      if( ConstantPressure ) then
        ! Isothermal compressibility
        Average = this%SumBetaT%Average
        Variance = this%SumBetaT%Variance
        write( IOBuffer, &
&         '("Isothermal compressibility", T29, "reduced:", 2F20.9)' ) &
&         Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T28, "in 1/MPa:", 2F20.9)' ) &
&         Average / ( UnitPressure * 1E-6_RK ), &
&         Variance / ( UnitPressure * 1E-6_RK )
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )
        ! dH/dP
        Average = this%SumdHdP%Average
        Variance = this%SumdHdP%Variance
        write( IOBuffer, '("dH/dP", T29, "reduced:", 2F20.9)' ) &
&         Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T28, "in l/mol:", 2F20.9)' ) &
&         Average / UnitDensity, Variance / UnitDensity
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )
        ! CP
        Average = this%SumCP%Average
        Variance = this%SumCP%Variance
        write( IOBuffer, '("Isobaric heat capacity", T29, "reduced:", 2F20.9)' ) &
&         Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T24, "in J/(mol K):", 2F20.9)' ) &
&         Average * kBoltzmann * NAvogadro, &
&         Variance * kBoltzmann * NAvogadro
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )
        ! AlphaP
        Average = this%SumAlphaP%Average
        Variance = this%SumAlphaP%Variance
        write( IOBuffer, '("Volume expansivity", T29, "reduced:", 2F20.9)' ) &
&         Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T30, "in 1/K:", 2F20.9)' ) &
&         Average / UnitTemperature, Variance / UnitTemperature
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )

        ! Speed of sound
        molmass = 0._RK
        cpid = 0._RK
        do i = 1, this%NRealComponents
          pc => this%Component(i)
          molmass = molmass + pc%Fraction * pc%Molecule%Mass
          cpid = cpid + .5_RK * pc%Fraction * pc%Molecule%NDF
        end do
        Average = 1._RK / sqrt( molmass * ( this%SumBetaT%Average * &
&         this%SumDensity%Average - this%RefTemperature * &
&         this%SumAlphaP%Average**2 / ( this%SumCP%Average + cpid ) ) )
        Variance = .25_RK / molmass / ( this%SumBetaT%Average * &
&         this%SumDensity%Average - this%RefTemperature * &
&         this%SumAlphaP%Average**2 / ( this%SumCP%Average + cpid ) )**3 * &
&         ( this%SumDensity%Average**2 * this%SumBetaT%Variance**2 + &
&           this%SumBetaT%Average**2 * this%SumDensity%Variance**2 + &
&           this%RefTemperature**2 * this%SumAlphaP%Average**2 / &
&           ( this%SumCP%Average + cpid )**2 * &
&           ( 4._RK * this%SumAlphaP%Variance**2 + &
&             this%SumAlphaP%Average**2 / ( this%SumCP%Average + cpid )**2 * &
&             this%SumCP%Variance**2 ) )
        write( IOBuffer, '("Speed of sound", T29, "reduced:", 2F20.9)' ) &
&         Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T30, "in m/s:", 2F20.9)' ) &
&         Average * UnitLength / UnitTime, Variance * UnitLength / UnitTime
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )
      else
        ! dU/dV
        Average = this%SumdUdV%Average
        Variance = this%SumdUdV%Variance
        write( IOBuffer, '("dU/dV", T29, "reduced:", 2F20.9)' ) &
&         Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T30, "in MPa:", 2F20.9)' ) &
&         Average * UnitPressure * 1E-6_RK, &
&         Variance * UnitPressure * 1E-6_RK
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )
        ! Cv
        Average = this%SumCV%Average
        Variance = this%SumCV%Variance
        write( IOBuffer, '("Isochore heat capacity", T29, "reduced:", 2F20.9)' ) &
&         Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T24, "in J/(mol K):", 2F20.9)' ) &
&         Average * kBoltzmann * NAvogadro, &
&         Variance * kBoltzmann * NAvogadro
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )
      endif

    end if

    ! Separator
    write( IOBuffer, '(76("="))' )
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Too large cutoff radius
    write( IOBuffer, &
&     '("Cutoff radius is", I10, " times (", F6.2, "%) too large")' ) &
&     this%NRCutoffMax, ( 100._RK * this%NRCutoffMax ) / Step
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Separator
    write( IOBuffer, '(76("="))' )
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Phase equilibria data for GE-ensemble
    if( EnsembleType == EnsembleTypeGE ) then
      write( IOBuffer, '("PHASE EQUILIBRIA DATA")' )
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '("---------------------")' )
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Simulation temperature
      write( IOBuffer, '("Simulation temperature", T29, "reduced:", F20.9)' ) &
&       this%RefTemperature
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T32, "in K:", F20.9)' ) &
&       this%Temperature * UnitTemperature
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Molar fractions of liquid phase
      do i = 1, this%NComponents
        pc => this%Component(i)
        write( IOBuffer, &
&         '("Liquid molar fraction of ", A, T36, ":", F20.9)' ) &
&         trim( pc%Molecule%PotModFileName ), pc%LiqFraction
        call FileWrite( this%iounit_errors )
      end do
      call FileWriteBlank( this%iounit_errors )

      ! Simulation pressure of liquid phase
      write( IOBuffer, &
&       '("Liquid simulation pressure", T29, "reduced:", F20.9)' ) &
&       this%RefPressure
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T30, "in MPa:", F20.9)' ) &
&       this%RefPressure * UnitPressure * 1e-6_RK
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Vapor pressure
      Average = this%SumPressure%Average
      Variance = this%SumPressure%Variance
      NN = 0._RK
      do i = 1, this%NComponents
        NN = NN + this%Component(i)%SumFraction%Average * &
&         this%Component(i)%PartialMolarVolume
      end do
      NN = NN - this%Temperature / Average
      do i = 1, this%NComponents
        pc => this%Component(i)
        dpdmu(i) = -this%Temperature * pc%SumFraction%Average / NN
        dpdv(i) = -pc%SumFraction%Average * (Average-this%RefPressure) / NN
        varmu(i) = pc%VarChemPot
        varv(i) = pc%VarPartialMolarVolume
      end do
      VarPressure = sqrt( Variance**2 + sum( (dpdmu * varmu)**2 ) + &
&       sum( (dpdv * varv)**2 ) )
      write( IOBuffer, '("Vapor pressure", T29, "reduced:", 2F20.9)' ) &
&       Average, VarPressure
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T30, "in MPa:", 2F20.9)' ) &
&       Average * UnitPressure * 1E-6_RK, VarPressure * UnitPressure * 1E-6_RK
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Molar fractions of vapor phase
      do i = 1, this%NComponents
        pc => this%Component(i)
        yvi = pc%SumFraction%Average * &
&         ( pc%PartialMolarVolume / this%Temperature - &
&         1 / this%SumPressure%Average )
        do j = 1, this%NComponents
          dydmu(i, j) = yvi * dpdmu(j)
          dydv(i, j) = yvi * dpdv(j)
        end do
        dydmu(i, i) = dydmu(i, i) + pc%SumFraction%Average
        dydv(i, i) = dydv(i, i) + pc%SumFraction%Average * &
&         1 / this%Temperature * ( this%SumPressure%Average - this%RefPressure )
      end do

      do i = 1, (this%NComponents - 1)
        pc => this%Component(i)
        Average = pc%SumFraction%Average
        vary(i) = sqrt( pc%SumFraction%Variance**2 + &
&         sum( (dydmu(i, :) * varmu)**2 ) + sum( (dydv(i, :) * varv)**2 ) )
        write( IOBuffer, &
&         '("Vapor molar fraction of ", A, T36, ":", 2F20.9)' ) &
&         trim( pc%Molecule%PotModFileName ), Average, vary(i)
        call FileWrite( this%iounit_errors )
      end do
      pc => this%Component( this%NComponents )
      Average = pc%SumFraction%Average
      Variance = sqrt( sum( vary(1:(this%NComponents - 1))**2 ) )
      write( IOBuffer, &
&       '("Vapor molar fraction of ", A, T36, ":", 2F20.9)' ) &
&       trim( pc%Molecule%PotModFileName ), Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Saturated liquid density
      Average = this%LiqDensity + this%LiqDensity * this%LiqBetaT * &
&       ( this%SumPressure%Average - this%RefPressure )
      Variance = sqrt( this%VarLiqDensity**2 + ( this%VarLiqBetaT * &
&       ( this%SumPressure%Average - this%RefPressure ) + &
&       VarPressure * this%LiqBetaT )**2 )
      write( IOBuffer, '("Liquid density", T29, "reduced:", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in mol/l:", 2F20.9)' ) &
&       Average * UnitDensity, Variance * UnitDensity
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Saturated vapor density
      Average = this%SumDensity%Average
      Variance = this%SumDensity%Variance
      write( IOBuffer, '("Vapor density", T29, "reduced:", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in mol/l:", 2F20.9)' ) &
&       Average * UnitDensity, Variance * UnitDensity
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Saturated liquid enthalpy
      Average = this%LiqEnthalpy + this%LiqdHdP * &
&       ( this%SumPressure%Average - this%RefPressure )
      Variance = sqrt( this%VarLiqEnthalpy**2 + ( this%VarLiqdHdP * &
&       ( this%SumPressure%Average - this%RefPressure ) + &
&       VarPressure * this%LiqdHdP )**2 )
      write( IOBuffer, '("Liquid enthalpy", T29, "reduced:", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) &
&       Average * UnitEnergy * NAvogadro, &
&       Variance * UnitEnergy * NAvogadro
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      DeltaHv = Average
      VarDeltaHv = Variance

      ! Saturated vapor enthalpy
      Average = this%SumEnthalpy%Average
      Variance = this%SumEnthalpy%Variance
      write( IOBuffer, '("Vapor enthalpy", T29, "reduced:", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) &
&       Average * UnitEnergy * NAvogadro, &
&       Variance * UnitEnergy * NAvogadro
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      DeltaHv = Average - DeltaHv
      VarDeltaHv = Variance + VarDeltaHv

      ! Evaporation enthalpy
      write( IOBuffer, &
&       '("Enthalpy of vaporization", T29, "reduced:", 2F20.9)' ) &
&       DeltaHv, VarDeltaHv
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) &
&       DeltaHv * UnitEnergy * NAvogadro, &
&       VarDeltaHv * UnitEnergy * NAvogadro
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Separator
      write( IOBuffer, '(76("="))' )
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

    end if

    if( SimulationType .eq. MonteCarlo ) then
      ! Statistics section
      write( IOBuffer, '("Statistics")' )
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '("----------")' )
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Volume change acceptance rate and maximum displacement
      if( ConstantPressure ) then
        write( IOBuffer, '("Acceptance rate volume changes", T32, "in %:", F20.9)' ) &
&         100._RK * real( this%NResizeSuccesses, RK ) / &
&         real ( this%NResizeAttempts, RK )
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '("Maximum displacement volume", T33, "r`d:", F20.9)' ) &
          this%DispVol
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )
      end if

      do i = 1, this%NRealComponents
        pc => this%Component(i)

        ! Move and rotate acceptance rates
        write( IOBuffer, '("Component ", A)' ) pc%PotModFileName
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '("Acceptance rate moves", T32, "in %:", F20.9)' ) &
&         100._RK * real( pc%NMoveSuccesses, RK ) / real ( pc%NMoveAttempts, RK )
        call FileWrite( this%iounit_errors )
        if( pc%Molecule%IsElongated ) then
          write( IOBuffer, '(T17, "rotates", T32, "in %:", F20.9)' ) 100._RK &
&           * real( pc%NRotateSuccesses, RK ) / real ( pc%NRotateAttempts, RK )
          call FileWrite( this%iounit_errors )
        end if

        if( pc%NMoveBiasedAttempts > 0 ) then
          ! Biased move and rotate acceptance rates
          write( IOBuffer, '(T17, "biased moves", T32, "in %:", F20.9)' ) &
&          100._RK * real( pc%NMoveBiasedSuccesses, RK ) / &
&          real ( pc%NMoveBiasedAttempts, RK )
          call FileWrite( this%iounit_errors )
          if( pc%Molecule%IsElongated ) then
            write( IOBuffer, '(T17, "biased rotates", T32, "in %:", F20.9)' ) &
&             100._RK * real( pc%NRotateBiasedSuccesses, RK ) / &
&             real ( pc%NRotateBiasedAttempts, RK )
            call FileWrite( this%iounit_errors )
          end if
        end if

        ! Maximum translational and rotational displacements
        write( IOBuffer, '("Maximum displacement trans.", T33, "r`d:", F20.9)' ) &
          pc%DispTran
        call FileWrite( this%iounit_errors )
        if( pc%Molecule%IsElongated ) then
          write( IOBuffer, '(T22, "rotational", T33, "r`d:", F20.9)' ) &
            pc%DispRot
          call FileWrite( this%iounit_errors )
        end if
        call FileWriteBlank( this%iounit_errors )

        ! Gradual insertion change fluctuating particle acceptance rates
        if( pc%ChemPotMethod .eq. ChemPotMethodGradIns ) then
          write(IOBuffer, &
&           '("Acceptance rate gradual insertion change fluctuating particle moves:")')
          call FileWrite( this%iounit_errors )
          write(IOBuffer, '("  up        down (%)")')
          call FileWrite( this%iounit_errors )
          write(IOBuffer, '("  --------  --------")')
          call FileWrite( this%iounit_errors )
          do j = 1, pc%NFluctMax
            write(IOBuffer, '(2F10.4)') &
&             real(pc%NFluctUpSuccesses(j), RK) / &
&               real(pc%NFluctUpAttempts(j), RK) * 100._RK, &
&             real(pc%NFluctDownSuccesses(j), RK) / &
&               real(pc%NFluctDownAttempts(j), RK) * 100._RK
            call FileWrite( this%iounit_errors )
          end do
        end if
      end do

      ! Inserts and deletes acceptance rates
      if( EnsembleType .eq. EnsembleTypeGE .or. &
&         EnsembleType .eq. EnsembleTypeHA ) then
        write( IOBuffer, '("Acceptance rate inserts", T32, "in %:", F20.9)' ) &
&         100._RK * real( this%NInsertSuccesses, RK ) / real ( this%NInsertAttempts, RK )
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '("Acceptance rate deletes", T32, "in %:", F20.9)' ) &
&         100._RK * real( this%NDeleteSuccesses, RK ) / real ( this%NDeleteAttempts, RK )
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )
      end if
    end if

#ifdef ABL
    counter = 0
    do i=1,this%NComponents,1
      do j=1,this%Component(i)%Molecule%NLJ126
        counter = counter + 1
        write(IOBuffer, '("Molecule  ", T5, "Site  ", T5)' ), i, j
        call FileWrite( this%iounit_errors )
        write(IOBuffer, '("dp / dsigma", F20,9)' ), this%AblPS(i,j)
        call FileWrite( this%iounit_errors )
        write(IOBuffer, '("dp / deps", F20,9)' ), this%AblPE(i,j)
        call FileWrite( this%iounit_errors )
        write(IOBuffer, '("drho / dsigma", F20,9)' ), this%AblRhoS(i,j)
        call FileWrite( this%iounit_errors )
        write(IOBuffer, '("drho / deps", F20,9)' ), this%AblRhoE(i,j)
        call FileWrite( this%iounit_errors )
      end do
    end do
#endif
    ! Close final result file
    call FileClose( this%iounit_errors )





  end subroutine TEnsemble_ErrorsUpdate



!==============================================================!
!  Subroutine TEnsemble_SVCOutput                              !
!==============================================================!

  subroutine TEnsemble_SVCOutput( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK) :: value
    integer  :: i, j

    ! Open final result file
    write( IOBuffer, '(I16)' ) this%EnsembleNumber
    call FileRewrite( this%iounit_errors, &
&     trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ErrorsFileExtension )

    ! Separator
    write( IOBuffer, '(76("="))' )
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )
    write( IOBuffer, '(T24, "SIMULATION RESULT FILE")' )
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T24, "----------------------")' )
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Simulation type
    write( IOBuffer, '("Simulation type", T36, ":", 9X, A)' ) &
&     trim( SimulationTypeString )
    call FileWrite( this%iounit_errors )

    ! Number of orientations
    write( IOBuffer, '("Number of orientations", T36, ":", I10)' ) NOrient
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Number of radial steps", T36, ":", I10)' ) NSteps
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Minimum radius", T29, "reduced:", F20.9)' ) MinRadius
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T32, "in A:", F20.9)' ) &
&     MinRadius * UnitLength / Angstroem
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Maximum radius", T29, "reduced:", F20.9)' ) MaxRadius
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T32, "in A:", F20.9)' ) &
&     MaxRadius * UnitLength / Angstroem
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Temperature
    write( IOBuffer, '("Temperature", T29, "reduced:", F20.9)' ) &
&     this%Temperature
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T32, "in K:", F20.9)' ) &
&     this%Temperature * UnitTemperature
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! System of units
    write( IOBuffer, '("Unit of length", T36, ":", F20.9, " A")' ) &
&     UnitLength / Angstroem
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Unit of energy", T36, ":", F20.9, " K")' ) &
&     UnitEnergy / kBoltzmann
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Unit of mass", T36, ":", F20.9, " a.u.")' ) &
&     UnitMass * NAvogadro * 1000._RK
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Separator
    write( IOBuffer, '(76("="))' )
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )
    write( IOBuffer, '("VALUE", T31, "UNITS", T46, "AVERAGE")' )
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("-----", T31, "-----", T46, "-------")' )
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Second virial coefficient
    do i = 1, this%NComponents, 2
      do j = i + 1, this%NComponents, 2
        value = this%Interaction(i, j)%IntFFunction(NSteps) + &
&         .5_RK * this%Interaction(i, j)%EPotCorrLJ / this%Temperature
        write( IOBuffer, '("2. VC of ", A, "-", A, T29, "reduced:", F20.9)' ) &
&         trim( this%Component(i)%Molecule%PotModFileName ), &
&         trim( this%Component(j)%Molecule%PotModFileName ), value
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T28, "in l/mol:", F20.9)' ) value / UnitDensity
        call FileWrite( this%iounit_errors )
      end do
    end do
    call FileWriteBlank( this%iounit_errors )

    ! Temperature deviation of second virial coefficient
    do i = 1, this%NComponents, 2
      do j = i + 1, this%NComponents, 2
        value = ( this%Interaction(i, j)%IntFFunction2(NSteps) &
&         - this%Interaction(i,j)%IntFFunction1(NSteps) ) &
&         / ( .0002_RK * this%Temperature )
        write( IOBuffer, '("dB/dT of ", A, "-", A, T29, "reduced:", F20.9)' ) &
&         trim( this%Component(i)%Molecule%PotModFileName ), &
&         trim( this%Component(j)%Molecule%PotModFileName ), value
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T24, "in l/(mol K):", F20.9)' ) &
&         value / ( UnitDensity * UnitTemperature )
        call FileWrite( this%iounit_errors )
      end do
    end do
    call FileWriteBlank( this%iounit_errors )

    ! Separator
    write( IOBuffer, '(76("="))' )
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Close final result file
    call FileClose( this%iounit_errors )

  end subroutine TEnsemble_SVCOutput



!==============================================================!
!  Subroutine TEnsemble_VisualOpen                             !
!==============================================================!

  subroutine TEnsemble_VisualOpen( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer                   :: i, j
    type(TSiteLJ126), pointer :: psLJ126

    ! Open visualization file
    write( IOBuffer, '(I16)' ) this%EnsembleNumber
    call FileRewrite( this%iounit_visual, &
&     trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//VisualFileExtension )

    ! Create header
    do i = 1, this%NComponents
      do j = 1, this%Component(i)%Molecule%NLJ126
        psLJ126 => this%Component(i)%Molecule%SiteLJ126(j)
        write( IOBuffer, '("~", I3, " LJ", 4F8.4, "  1")' ) i, &
&         psLJ126%r(:) * UnitLength / Angstroem, &
&         psLJ126%sig  * UnitLength / Angstroem
        call FileWrite( this%iounit_visual )
      end do
    end do
    call FileWriteBlank( this%iounit_visual )

  end subroutine TEnsemble_VisualOpen



!==============================================================!
!  Subroutine TEnsemble_VisualUpdate                           !
!==============================================================!

  subroutine TEnsemble_VisualUpdate( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer  :: i, j
    logical  :: l
    real(RK) :: r(3), q(4)

    ! Update visualization file
    write( IOBuffer, '("#", F10.4, "  new Frame")' ) &
&     this%BoxLength * UnitLength / Angstroem
    call FileWrite( this%iounit_visual )
    do i = 1, this%NComponents
      l = this%Component(i)%Molecule%isElongated
      do j = 1, this%Component(i)%NPart
        r(:) = this%Component(i)%P0(j, :)
        if( l ) then
          q(:) = this%Component(i)%Q0(j, :)
        else
          q(1) = 1._RK
          q(2:4) = 0._RK
        end if
        where( r(:) < 0._RK ) r(:) = r(:) + 1._RK
        write( IOBuffer, '("!", I3, 3I4, 4I5)' ) &
&         i, nint( r(:) * 999 ), nint( q(:) * 999 )
        call FileWrite( this%iounit_visual )
      end do
    end do
    call FileWriteBlank( this%iounit_visual )

  end subroutine TEnsemble_VisualUpdate



!==============================================================!
!  Subroutine TEnsemble_VisualClose                            !
!==============================================================!

  subroutine TEnsemble_VisualClose( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Close visualization file
    write( IOBuffer, '("##")' )
    call FileWrite( this%iounit_visual )
    call FileClose( this%iounit_visual )

  end subroutine TEnsemble_VisualClose



!==============================================================!
!  Subroutine TEnsemble_RestartSave                            !
!==============================================================!

  subroutine TEnsemble_RestartSave( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    type(TComponent), pointer :: pc
    integer                   :: i

    ! Save contents to restart file
    write( iounit_restart, '(I10)' ) this%NPart
    write( iounit_restart, '(ES20.12E3)' ) this%Volume0

    if( SimulationType .eq. MolecularDynamics ) then
      write( iounit_restart, '(ES20.12E3)' ) this%Volume1
      write( iounit_restart, '(ES20.12E3)' ) this%Volume2

      if( IntegratorType .eq. IntegratorTypeGear ) then
        write( iounit_restart, '(ES20.12E3)' ) this%Volume3
        write( iounit_restart, '(ES20.12E3)' ) this%Volume4
        write( iounit_restart, '(ES20.12E3)' ) this%Volume5
      end if
    else
      write( iounit_restart, '(ES20.12E3)' ) this%DispVol
      write(iounit_restart, '(2I10)' ) &
&       this%NResizeAttempts, this%NResizeSuccesses
      if( EnsembleType .eq. EnsembleTypeGE .or. &
&         EnsembleType .eq. EnsembleTypeHA ) then
        write(iounit_restart, '(2I10)' ) &
&         this%NInsertAttempts, this%NInsertSuccesses
        write(iounit_restart, '(2I10)' ) &
&         this%NDeleteAttempts, this%NDeleteSuccesses
      end if
    end if
    write( iounit_restart, '(I10)' ) this%NRCutoffMax

    ! Save components
    do i = 1, this%NComponents
      call RestartSave( this%Component(i) )
    end do

    ! Save accumulators
    ! 1.) Basic sums
    call RestartSave( this%SumPressure )
    call RestartSave( this%SumDensity )
    call RestartSave( this%SumTemperature )
    call RestartSave( this%SumEPot )
    call RestartSave( this%SumEnthalpy )
    call RestartSave( this%SumVolume )
    call RestartSave( this%SumVirial )
    if( EnsembleType .eq. EnsembleTypeGE .or. &
&       EnsembleType .eq. EnsembleTypeHA ) then
      call RestartSave( this%SumNPart )
      do i = 1, this%NComponents
        pc => this%Component(i)
        call RestartSave( pc%SumFraction )
      end do
    end if

    ! 2.) Combined sums
    call RestartSave( this%SumEPotSquared )
    call RestartSave( this%SumEPotV )
    call RestartSave( this%SumEPotVirial )
    call RestartSave( this%SumEnthalpySquared )
    call RestartSave( this%SumEnthalpyV )
    call RestartSave( this%SumVolumeSquared )

    ! 3.) Derived sums
    if( ConstantPressure ) then
      call RestartSave( this%SumBetaT )
      call RestartSave( this%SumdHdP )
      call RestartSave( this%SumCP )
      call RestartSave( this%SumAlphaP )
    else
      call RestartSave( this%SumdUdV )
      call RestartSave( this%SumCV )
    endif

    ! 4.) Chemical potential and partial molar volumes
    do i = 1, this%NRealComponents
      pc => this%Component(i)
      select case( pc%ChemPotMethod )
      case( ChemPotMethodGradIns )
        call RestartSave( pc%SumInvChemPotRho )
        call RestartSave( pc%SumInvChemPot )
!DEBUG
        call RestartSave( pc%SumInvChemPotRho1 )
        call RestartSave( pc%SumInvChemPot1 )
        call RestartSave( pc%SumInvChemPotRho2 )
        call RestartSave( pc%SumInvChemPot2 )
!DEBUG
      case( ChemPotMethodWidom )
        call RestartSave( pc%SumChemPotV )
        call RestartSave( pc%SumChemPotVV )
      end select
      if( pc%ChemPotMethod .ne. ChemPotMethodNone .and. ConstantPressure &
&      .and. this%NRealComponents > 1 ) then
        call RestartSave( pc%SumVW )
!DEBUG
        if( pc%ChemPotMethod .eq. ChemPotMethodGradIns ) then
          call RestartSave( pc%SumVW1 )
          call RestartSave( pc%SumVW2 )
        end if
!DEBUG
      end if
    end do

  end subroutine TEnsemble_RestartSave



!==============================================================!
!  Subroutine TEnsemble_RestartRead                            !
!==============================================================!

  subroutine TEnsemble_RestartRead( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    type(TComponent), pointer :: pc
    integer                   :: i,j,stat

    if( RootProc ) then

      ! Read contents from restart file
      read( iounit_restart, '(I10)' ) this%NPart
      read( iounit_restart, '(ES20.12E3)' ) this%Volume0

      if( SimulationType .eq. MolecularDynamics ) then
        read( iounit_restart, '(ES20.12E3)' ) this%Volume1
        read( iounit_restart, '(ES20.12E3)' ) this%Volume2

        if( IntegratorType .eq. IntegratorTypeGear ) then
          read( iounit_restart, '(ES20.12E3)' ) this%Volume3
          read( iounit_restart, '(ES20.12E3)' ) this%Volume4
          read( iounit_restart, '(ES20.12E3)' ) this%Volume5
        end if
      else
        read( iounit_restart, '(ES20.12E3)' ) this%DispVol
        read(iounit_restart, '(2I10)' ) &
&         this%NResizeAttempts, this%NResizeSuccesses
        if( EnsembleType .eq. EnsembleTypeGE .or. &
&           EnsembleType .eq. EnsembleTypeHA ) then
          read(iounit_restart, '(2I10)' ) &
&           this%NInsertAttempts, this%NInsertSuccesses
          read(iounit_restart, '(2I10)' ) &
&           this%NDeleteAttempts, this%NDeleteSuccesses
        end if
      end if
      read( iounit_restart, '(I10)' ) this%NRCutoffMax

    end if

#if MPI_VER > 0
    call MPI_Bcast( this%NPart, 1, MPI_INTEGER, NRootProc, &
&     MPI_COMM_WORLD, ierror )
    call MPI_Bcast( this%Volume0, 1, MPI_DOUBLE_PRECISION, NRootProc, &
&     MPI_COMM_WORLD, ierror )
    if( SimulationType .eq. MonteCarlo ) then
      call MPI_Bcast( this%DispVol, 1, MPI_DOUBLE_PRECISION, NRootProc, &
&       MPI_COMM_WORLD, ierror )
      call MPI_Bcast( this%NResizeAttempts, 1, MPI_INTEGER, NRootProc, &
&       MPI_COMM_WORLD, ierror )
      call MPI_Bcast( this%NResizeSuccesses, 1, MPI_INTEGER, NRootProc, &
&       MPI_COMM_WORLD, ierror )
      if( EnsembleType .eq. EnsembleTypeGE .or. &
&         EnsembleType .eq. EnsembleTypeHA ) then
        call MPI_Bcast( this%NInsertAttempts, 1, MPI_INTEGER, NRootProc, &
&         MPI_COMM_WORLD, ierror )
        call MPI_Bcast( this%NInsertSuccesses, 1, MPI_INTEGER, NRootProc, &
&         MPI_COMM_WORLD, ierror )
        call MPI_Bcast( this%NDeleteAttempts, 1, MPI_INTEGER, NRootProc, &
&         MPI_COMM_WORLD, ierror )
        call MPI_Bcast( this%NDeleteSuccesses, 1, MPI_INTEGER, NRootProc, &
&         MPI_COMM_WORLD, ierror )
      end if
    end if
    call MPI_Bcast( this%NRCutoffMax, 1, MPI_INTEGER, NRootProc, &
&     MPI_COMM_WORLD, ierror )
#endif

    ! Read components
    do i = 1, this%NComponents
      call RestartRead( this%Component(i) )
    end do

    ! Read accumulators
    ! 1.) Basic sums
    call RestartRead( this%SumPressure )
    call RestartRead( this%SumDensity )
    call RestartRead( this%SumTemperature )
    call RestartRead( this%SumEPot )
    call RestartRead( this%SumEnthalpy )
    call RestartRead( this%SumVolume )
    call RestartRead( this%SumVirial )
    if( EnsembleType .eq. EnsembleTypeGE .or. &
&       EnsembleType .eq. EnsembleTypeHA ) then
      call RestartRead( this%SumNPart )
      do i = 1, this%NComponents
        pc => this%Component(i)
        call RestartRead( pc%SumFraction )
      end do
    end if

    ! 2.) Combined sums
    call RestartRead( this%SumEPotSquared )
    call RestartRead( this%SumEPotV )
    call RestartRead( this%SumEPotVirial )
    call RestartRead( this%SumEnthalpySquared )
    call RestartRead( this%SumEnthalpyV )
    call RestartRead( this%SumVolumeSquared )

    ! 3.) Derived sums
    if( ConstantPressure ) then
      call RestartRead( this%SumBetaT )
      call RestartRead( this%SumdHdP )
      call RestartRead( this%SumCP )
      call RestartRead( this%SumAlphaP )
    else
      call RestartRead( this%SumdUdV )
      call RestartRead( this%SumCV )
    endif

    ! 4.) Chemical potential and partial molar volumes
    do i = 1, this%NRealComponents
      pc => this%Component(i)
      select case( pc%ChemPotMethod )
      case( ChemPotMethodGradIns )
        call RestartRead( pc%SumInvChemPotRho )
        call RestartRead( pc%SumInvChemPot )
!DEBUG
        call RestartRead( pc%SumInvChemPotRho1 )
        call RestartRead( pc%SumInvChemPot1 )
        call RestartRead( pc%SumInvChemPotRho2 )
        call RestartRead( pc%SumInvChemPot2 )
!DEBUG
        pc%NFluctState = 0;
        do j = 1, this%NComponents
          if (this%Component(j)%NPart .eq. 1) then
            pc%NFluctState = j-1
          end if
        end do

      case( ChemPotMethodWidom )
        call RestartRead( pc%SumChemPotV )
        call RestartRead( pc%SumChemPotVV )
      end select
      if( pc%ChemPotMethod .ne. ChemPotMethodNone .and. ConstantPressure &
&       .and. this%NRealComponents > 1 ) then
        call RestartRead( pc%SumVW )
!DEBUG
        if( pc%ChemPotMethod .eq. ChemPotMethodGradIns ) then
          call RestartRead( pc%SumVW1 )
          call RestartRead( pc%SumVW2 )
        end if
!DEBUG
      end if
    end do

    ! Update density and box length
    call UpdateBoxLength( this )

    ! Update molar fractions
    call UpdateFractions( this )

    if (LongRange .eq. Ewald) then
    ! Calculate initial energies for the Ewald Summation
       this%Kappa = this%KappaL/this%BoxLength   !Boxlength bereits normiert
       do i=1,this%NComponents
         do j=1,this%NComponents
           this%Interaction(i,j)%Kappa = this%Kappa
         end do
       end do
      ! Calculate initial energies for the Ewald Summation



!       if (LongRange .eq. Ewald) then
!          write (*,*) 'this%BoxLength wurde festgelegt in Zeile ms2_ensemble 958'
        ! Memory Allocation for Ewald Summation
         allocate(this%Faktor(this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error Faktor'
         allocate(this%U_fourierLocal(this%BoxenAnzahl),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error U_fourier'
         allocate(this%SSin(this%BoxenAnzahl),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error SSin'
         allocate(this%SCos(this%BoxenAnzahl),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error SCos'
         allocate(this%sinfac(this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error sinfac'
         allocate(this%cosfac(this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error cosfac'
         allocate(this%SSin_Vec(this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error SSin_Vec'
         allocate(this%SCos_Vec(this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error SCos_Vec'
         allocate(this%sinfac_s(this%NComponents,5,this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error sin_facs'
         allocate(this%cosfac_s(this%NComponents,5,this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error cosfac_s'
         allocate(this%Vec2(this%BoxenAnzahl),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error Vec2'
         allocate(this%HFac(this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error HFac'
         allocate(this%distx(this%NComponents,5,this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error distx'
         allocate(this%disty(this%NComponents,5,this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error disty'
         allocate(this%distz(this%NComponents,5,this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error distz'
         allocate(this%VirIntra(this%NPart),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error VirIntra'
         allocate(this%rold(5,3),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error rold'
         DO i=1,5
           DO j=1,3
             this%rold(i,j) = 0._RK
             this%rold(i,j) = 0._RK
             this%rold(i,j) = 0._RK
           END DO
         END DO


    else if (LongRange .eq. PME) then
    ! Calculate initial energies for the Ewald Summation
       this%Kappa = this%KappaL/this%BoxLength   !Boxlength bereits normiert
       do i=1,this%NComponents
         do j=1,this%NComponents
           this%Interaction(i,j)%Kappa = this%Kappa
         end do
       end do

! Allocation needed!
         allocate(this%bsp_arr(this%gridx+5),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error bsp_arr'
         allocate(this%bsp_modx(this%gridx+5),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error bsp_modx'
         allocate(this%bsp_mody(this%gridx+5),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error bsp_mody'
         allocate(this%bsp_modz(this%gridx+5),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error bsp_modz'

       call PMESetup(this)
    end if

    if( SimulationType .eq. MolecularDynamics ) then

      ! Calculate temperature and set up ndf
      call CalculateEKin( this, .false. )

    else

      ! Set temperature
      this%Temperature = this%RefTemperature

      ! Initialize energy matrix
      call Mol2Atom( this )
      call Energy( this, this%EPot )
      call UpdateEnergy( this )

    end if

  end subroutine TEnsemble_RestartRead



!==============================================================!
!  Subroutine TEnsemble_DbgWrite                               !
!==============================================================!

  subroutine TEnsemble_DbgWrite( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare external procedures
    integer, external :: system

    ! Declare local variables
    type(TComponent), pointer :: pc
    type(TSiteLJ126), pointer :: plj
    integer                   :: nc, i, i1, i2, j, n, n2, n3
!     real(RK)                  :: a1, a2, a3, anorm, EPot
    real(RK)                  :: C(this%NPart * 3), Q(this%NPart * 4)

    if( .not. RootProc ) return

    n = this%NPart
    n2 = 2 * n
    n3 = 3 * n
    i1 = 1
    i2 = 0

    do nc = 1, this%NComponents
      pc => this%Component(nc)
      plj => pc%Molecule%SiteLJ126(1)
      i2 = i2 + pc%NPart
      do i = i1, i2
        j = i - i1 + 1
        C(i) = pc%P0(j, 1)
        C(i+n) = pc%P0(j, 2)
        C(i+n2) = pc%P0(j, 3)
        Q(i) = pc%Q0(j, 1)
        Q(i+n) = pc%Q0(j, 2)
        Q(i+n2) = pc%Q0(j, 3)
        Q(i+n3) = pc%Q0(j, 4)
      end do
      i1 = i2 + 1
    end do

    open(996, file='EIN_1', status='REPLACE', action='WRITE')
    write(996, '(3ES22.12E3)') C
    write(996, '(4ES22.12E3)') Q
    write(996,*) this%Component(1)%NFLuctMax - this%Component(1)%NFluctState + 1
    write(996,*) ix, iy
    close(996)

!     EPot = GetEnergy( this )
!     write(*,*) (EPot - this%Density * this%EPotCorrLJ - &
! &     this%EPotCorrRF) / real( this%NPart, 8 )
! 
!     i = system('./igi')
! 
!     stop

  end subroutine TEnsemble_DbgWrite




!==============================================================!
!  Subroutine TSimulation_Ewald_SelfTerm                       !
!==============================================================!

  subroutine TEnsemble_EwaldSelfTerm( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif
    include 'fftw3.f'
    ! Declare arguments
    type(TEnsemble)         :: this
    ! Declare local variables
    integer :: i,Si,Sj
    integer :: NX,NY,NZ
    real(RK):: USelbstTermKomp
!     real(RK):: ZweiPi_inv
    real(RK):: Faktor, NSQ
    real(RK):: RXi,RYi,RZi,RXj,RYj,RZj
    real(RK):: drxij,dryij,drzij,dr
    real(RK):: approx
    real(RK):: UIntraTermKomp
    real(RK):: twopi
    real(RK):: test

!     ZweiPi_inv = 2.0 / sqrt(PI)
    twopi = 2*PI


#if DEBUG_4Part > 0
! DEBUG REASONS
    this%Component(1)%Molecule%SiteCharge(1)%RX(1) = 0.5041221762148172
    this%Component(1)%Molecule%SiteCharge(1)%RY(1) = 0.8107396528595158
    this%Component(1)%Molecule%SiteCharge(1)%RZ(1) = 0.1894949962472363
    this%Component(1)%Molecule%SiteCharge(2)%RX(1) = 0.502268554151700
    this%Component(1)%Molecule%SiteCharge(2)%RY(1) = 0.8117808162647254
    this%Component(1)%Molecule%SiteCharge(2)%RZ(1) = 0.186927667845592
    this%Component(1)%Molecule%SiteCharge(1)%RX(2) = 2.099380236705095
    this%Component(1)%Molecule%SiteCharge(1)%RY(2) = 1.860649777419017
    this%Component(1)%Molecule%SiteCharge(1)%RZ(2) = 1.313235759816588
    this%Component(1)%Molecule%SiteCharge(2)%RX(2) = 2.097708401492369
    this%Component(1)%Molecule%SiteCharge(2)%RY(2) = 1.858592437418999
    this%Component(1)%Molecule%SiteCharge(2)%RZ(2) = 1.311215009892196
    this%Component(1)%Molecule%SiteCharge(1)%RX(3) = 1.317611713633559
    this%Component(1)%Molecule%SiteCharge(1)%RY(3) = 2.087579656401032
    this%Component(1)%Molecule%SiteCharge(1)%RZ(3) = 0.4484239832120850
    this%Component(1)%Molecule%SiteCharge(2)%RX(3) = 1.315500144816407
    this%Component(1)%Molecule%SiteCharge(2)%RY(3) = 2.085384351686846
    this%Component(1)%Molecule%SiteCharge(2)%RZ(3) = 0.4470700905082808
    this%Component(1)%Molecule%SiteCharge(1)%RX(4) = 1.482005925799646
    this%Component(1)%Molecule%SiteCharge(1)%RY(4) = 0.7714424297398743
    this%Component(1)%Molecule%SiteCharge(1)%RZ(4) = 1.556938840875683
    this%Component(1)%Molecule%SiteCharge(2)%RX(4) = 1.484582190214423
    this%Component(1)%Molecule%SiteCharge(2)%RY(4) = 0.7716945714654903
    this%Component(1)%Molecule%SiteCharge(2)%RZ(4) = 1.554838746276274
! 
! 
! 
! 
! DEBUG STEPHAN
    this%Component(1)%Molecule%SiteCharge(1)%RX(1) = 0.50412 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(1) = 0.81073 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(1) = 0.18949 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(1) = 0.50226 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(1) = 0.81178 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(1) = 0.18692 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(2) = 2.09938 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(2) = 1.86064 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(2) = 1.31323 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(2) = 2.09770 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(2) = 1.85859 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(2) = 1.31121 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(3) = 1.31761 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(3) = 2.08757 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(3) = 0.44842 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(3) = 1.31550 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(3) = 2.08538 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(3) = 0.44707 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(4) = 1.48200 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(4) = 0.77144 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(4) = 1.55693 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(4) = 1.48458 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(4) = 0.77169 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(4) = 1.55483 / this%BoxLength
#endif
    i = 0
    DO NX = 0, this%NMAX, 1
       IF (NX .EQ.0) THEN
         Faktor = 1.0
       ELSE
         Faktor = 2.0
       END IF
       DO NY = -this%NMAX,this%NMAX,1
         DO NZ = -this%NMAX,this%NMAX,1
           NSQ = NX*NX + NY*NY + NZ*NZ
           IF ( (NSQ .LE. this%NSQMAX) .AND. (NSQ .NE. 0) ) THEN
             i = i+1
             IF (i .GT. this%NVecMax) STOP   'BoxAnzahl zu gross'
!              test = (PI/this%KappaL)**2*NSQ
!              this%Ewald_Prefac(i) = Faktor*exp(-(PI/this%KappaL)**2 * NSQ) / (Pi2 * NSQ)
             this%Ewald_Prefac(i) = Faktor*exp(-(PI/this%KappaL)**2 * NSQ) / (Pi2 * NSQ * this%BoxLength)
             this%Ewald_Vec(1,i)  = twopi*NX
             this%Ewald_Vec(2,i)  = twopi*NY
             this%Ewald_Vec(3,i)  = twopi*NZ
           END IF
         END DO
       END DO
    END DO
    this%BoxenAnzahl = i

#if MPI_VER > 0
    call MPI_Bcast( this%Ewald_Vec(1,:), this%BoxenAnzahl, MPI_DOUBLE_PRECISION, &
&     NRootProc, MPI_COMM_WORLD, ierror )
    call MPI_Bcast( this%Ewald_Vec(2,:), this%BoxenAnzahl, MPI_DOUBLE_PRECISION, &
&     NRootProc, MPI_COMM_WORLD, ierror )
    call MPI_Bcast( this%Ewald_Vec(3,:), this%BoxenAnzahl, MPI_DOUBLE_PRECISION, &
&     NRootProc, MPI_COMM_WORLD, ierror )
    call MPI_Bcast( this%Ewald_Prefac, this%BoxenAnzahl, MPI_DOUBLE_PRECISION, &
&     NRootProc, MPI_COMM_WORLD, ierror )
#endif

! Selbstterm
    this%USelbstTerm = 0.0
    DO i=1,this%NComponents,1
       USelbstTermKomp = 0.0
       DO Si=1,this%Component(i)%Molecule%NCharge,1
         USelbstTermKomp = USelbstTermKomp + this%Component(i)%Molecule%SiteCharge(Si)%e**2
       END DO
       this%USelbstTerm = this%USelbstTerm + this%Component(i)%NPart * USelbstTermKomp
    END DO

!     this%USelbstTerm = -this%USelbstTerm * this%Kappa / sqrt(Pi) *UnitLength*UnitEnergy
    this%USelbstTerm = -this%USelbstTerm * this%Kappa / sqrt(Pi)


! intramolecular term
    this%UIntra = 0._RK
    DO i=1,this%NComponents,1
      UIntraTermKomp = 0.0
      DO Si = 1,this%component(i)%Molecule%NCharge-1
        DO Sj = Si+1,this%component(i)%Molecule%NCharge
!           mol = this%Molecule
!         Following lines are applicable for non-rigid molecules
!  Watch out: Before using these lines, they have to be shifted to some other spot
!             RX not assigned yet (done in Mol2Atom!)
!           RXi  = this%Component(i)%Molecule%SiteCharge(Si)%RX(1)
!           RYi  = this%Component(i)%Molecule%SiteCharge(Si)%RY(1)
!           RZi  = this%Component(i)%Molecule%SiteCharge(Si)%RZ(1)
!           RXj  = this%Component(i)%Molecule%SiteCharge(Sj)%RX(1)
!           RYj  = this%Component(i)%Molecule%SiteCharge(Sj)%RY(1)
!           RZj  = this%Component(i)%Molecule%SiteCharge(Sj)%RZ(1)
! 
!           drxij = (RXi-RXj)*this%BoxLength
!           dryij = (RYi-RYj)*this%BoxLength
!           drzij = (RZi-RZj)*this%BoxLength


          RXi  = this%Component(i)%Molecule%SiteCharge(Si)%r(1)
          RYi  = this%Component(i)%Molecule%SiteCharge(Si)%r(2)
          RZi  = this%Component(i)%Molecule%SiteCharge(Si)%r(3)
          RXj  = this%Component(i)%Molecule%SiteCharge(Sj)%r(1)
          RYj  = this%Component(i)%Molecule%SiteCharge(Sj)%r(2)
          RZj  = this%Component(i)%Molecule%SiteCharge(Sj)%r(3)

          drxij = (RXi-RXj)
          dryij = (RYi-RYj)
          drzij = (RZi-RZj)

! Debugging
!           drxij = (RXi-RXj)
!           dryij = (RYi-RYj)
!           drzij = (RZi-RZj)

          dr = sqrt(drxij*drxij + dryij*dryij + drzij*drzij)

          if (dr .ne. 0.0) then
            call erfc_approx (this%Kappa*dr, approx)

            UIntraTermKomp = UIntraTermKomp - this%Component(i)%Molecule%SiteCharge(Si)%e* &
&                   this%Component(i)%Molecule%SiteCharge(Sj)%e / dr * (1-approx)
          end if
        END DO
      END DO
      this%UIntra = this%UIntra + this%component(i)%NPart * UIntraTermKomp
    END DO

#if MPI_VER > 0
    this%NBox1 = 1 +(this%BoxenAnzahl - 1) / NProcs
    DO i=1,NProcs
      this%NBox0(i) = 1 + this%NBox1(i) * (i-1)
    END DO
    this%NBox2 = min(this%NBox0 + this%NBox1 - 1, this%BoxenAnzahl)
    this%NBox1 = this%NBox2 - this%NBox0
#endif
    end subroutine TEnsemble_EwaldSelfTerm



!==============================================================!
!  Subroutine TSimulation_Ewald_FourierTerm                    !
!==============================================================!

   subroutine TEnsemble_EwaldFourierTerm(this)

   implicit none

#if MPI_VER > 0
    include 'mpif.h'
#endif
   ! Declare arguments
   type(TEnsemble)   :: this

   integer :: i,j,l,m
   integer :: molec

!    real(RK),pointer:: U_fourier(:)
   real(RK),pointer:: RX(:),RY(:),RZ(:)
   real(RK),pointer:: PX(:),PY(:),PZ(:)
   real(RK),pointer:: FX(:),FY(:),FZ(:)
   real(RK),pointer:: q(:)
   real(RK) :: RXloc(this%NPart),RYloc(this%NPart),RZloc(this%NPart)
   real(RK) :: PXloc(this%NPart),PYloc(this%NPart),PZloc(this%NPart)

   real(RK):: KVec(3)
   real(RK):: EPotLocal
   real(RK):: Viriallocal, VirIntra
   real(RK):: SSinSum,SCosSum
   real(RK):: KappaL2, vorfac
   real(RK):: facx,facy,facz
!    real(RK),pointer:: SSin_Fac, SCos_fac
!   real(RK),pointer :: test(:)
#if MPI_VER > 0
   integer:: i0,i1
#endif

   type(TMolecule), pointer               :: mol
   integer:: stat



#if DEBUG_4Part > 0
! DEBUG STEPHAN
    this%Component(1)%Molecule%SiteCharge(1)%RX(1) = 0.50412 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(1) = 0.81073 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(1) = 0.18949 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(1) = 0.50226 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(1) = 0.81178 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(1) = 0.18692 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(2) = 2.09938 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(2) = 1.86064 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(2) = 1.31323 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(2) = 2.09770 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(2) = 1.85859 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(2) = 1.31121 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(3) = 1.31761 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(3) = 2.08757 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(3) = 0.44842 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(3) = 1.31550 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(3) = 2.08538 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(3) = 0.44707 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(4) = 1.48200 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(4) = 0.77144 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(4) = 1.55693 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(4) = 1.48458 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(4) = 0.77169 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(4) = 1.55483 / this%BoxLength
#endif

   EPotLocal = 0.0_RK
   VirialLocal = 0.0_RK
   KappaL2 = 1.0_RK/(2._RK*this%KappaL**2)
   vorfac = 2._RK/this%BoxLength

   this%VirIntra = 0._RK

   this%SSin=0._RK
   this%SCos=0._RK
#if MPI_VER > 0
   j=NProc + 1
   i0 = this%NBox0(j)
   i1 = this%NBox2(j)
!    i0 = 1
   DO i=i0,i1,1
#else
   DO i=1,this%BoxenAnzahl,1
# endif
     KVec = this%Ewald_Vec(:,i)
     this%SSin_Vec = 0._RK
     this%SCos_Vec = 0._RK
     DO j=1,this%NComponents,1
       mol => this%Component(j)%Molecule
       q => mol%SiteCharge(1:mol%NCharge)%e
       molec = this%Component(j)%NPart
       DO l=1,mol%NCharge
         RX => this%Component(j)%Molecule%SiteCharge(l)%RX(1:molec)
         RY => this%Component(j)%Molecule%SiteCharge(l)%RY(1:molec)
         RZ => this%Component(j)%Molecule%SiteCharge(l)%RZ(1:molec)
         PX => this%Component(j)%Molecule%SiteCharge(l)%PX(1:molec)
         PY => this%Component(j)%Molecule%SiteCharge(l)%PY(1:molec)
         PZ => this%Component(j)%Molecule%SiteCharge(l)%PZ(1:molec)
!           q => this%Component(j)%Molecule%SiteCharge(l)%e

         RXloc(1:molec) = RX(1:molec)
         RYloc(1:molec) = RY(1:molec)
         RZloc(1:molec) = RZ(1:molec)
         PXloc(1:molec) = PX(1:molec)
         PYloc(1:molec) = PY(1:molec)
         PZloc(1:molec) = PZ(1:molec)

         DO m=1,molec
         if (RX(m) < 0) RXloc(m) = RXloc(m) + 1._RK
         if (RY(m) < 0) RYloc(m) = RYloc(m) + 1._RK
         if (RZ(m) < 0) RZloc(m) = RZloc(m) + 1._RK
         if (PX(m) < 0) PXloc(m) = PXloc(m) + 1._RK
         if (PY(m) < 0) PYloc(m) = PYloc(m) + 1._RK
         if (PZ(m) < 0) PZloc(m) = PZloc(m) + 1._RK
         end DO
         this%distx(j,l,1:molec) = (RXloc - PXloc)*this%BoxLength
         this%disty(j,l,1:molec) = (RYloc - PYloc)*this%BoxLength
         this%distz(j,l,1:molec) = (RZloc - PZloc)*this%BoxLength

         this%Faktor(1:molec) = KVec(1) * RX + KVec(2)*RY + KVec(3)*RZ

         this%sinfac_s(j,l,1:molec) = sin(this%Faktor)
         this%cosfac_s(j,l,1:molec) = cos(this%Faktor)
         this%sinfac = q(l)*this%sinfac_s(j,l,1:molec)
         this%cosfac = q(l)*this%cosfac_s(j,l,1:molec)

         this%SSin_Vec = this%SSin_Vec + this%sinfac
         this%SCos_Vec = this%SCos_Vec + this%cosfac
       END DO
     END DO

     SSinSum = sum(this%SSin_Vec)
     SCosSum = sum(this%SCos_Vec)
     this%SSin(i) = SSinSum
     this%SCos(i) = SCosSum

! Forces
     Facx = KVec(1)*this%Ewald_Prefac(i)*vorfac
     Facy = KVec(2)*this%Ewald_Prefac(i)*vorfac
     Facz = KVec(3)*this%Ewald_Prefac(i)*vorfac

     DO j=1,this%NComponents,1
       mol => this%Component(j)%Molecule
       q => mol%SiteCharge(1:mol%NCharge)%e
       molec = this%Component(j)%NPart
       DO l=1,mol%NCharge
          FX => this%Component(j)%Molecule%SiteCharge(l)%FX
          FY => this%Component(j)%Molecule%SiteCharge(l)%FY
          FZ => this%Component(j)%Molecule%SiteCharge(l)%FZ


          this%HFac = q(l)*(this%sinfac_s(j,l,1:molec)*this%SCos(i) - this%cosfac_s(j,l,1:molec)*this%SSin(i))
!           test = Facy*HFac
          FX = FX + Facx*this%HFac
          FY = FY + Facy*this%HFac
          FZ = FZ + Facz*this%HFac

          this%VirIntra = this%VirIntra + Facx*this%HFac*this%distx(j,l,1:molec)+&
   &         Facy*this%HFac*this%disty(j,l,1:molec)+&
&            Facz*this%HFac*this%distz(j,l,1:molec)
       END DO
     END DO
   END DO ! Boxenschleife

! Finish Calculation

! Energy
   this%U_fourierLocal = this%Ewald_Prefac * (this%SSin*this%SSin + this%SCos*this%SCos)
! Virial
   this%Vec2 = this%Ewald_Vec(1,:)**2 + this%Ewald_Vec(2,:)**2 + this%Ewald_Vec(3,:)**2 


#if MPI_VER > 0
   call MPI_Reduce( sum(this%U_fourierLocal), EPotLocal, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
   call MPI_Reduce( EPotLocal - sum(this%U_fourierLocal*KappaL2*this%Vec2), VirialLocal, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
   call MPI_Reduce( sum(this%VirIntra), VirIntra, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
#else
   EPotLocal = sum(this%U_fourierLocal)
   VirialLocal = EPotLocal - sum(this%U_fourierLocal *KappaL2*this%Vec2)
   VirIntra = sum(this%VirIntra)
#endif
  if( RootProc ) then
   this%UFourier= EPotLocal
   this%EVirial = -(Viriallocal - VirIntra)*Third
!    this%EVirial = -(Viriallocal + VirIntra)*Third
  end if

   end subroutine TEnsemble_EwaldFourierTerm



!==============================================================!
!  Subroutine TSimulation_Ewald_SelfTerm_Energy                !
!==============================================================!

  subroutine TEnsemble_EwaldSelf_Energy( this )

    implicit none

    ! Declare arguments
    type(TEnsemble)            :: this
    type(TInteraction),pointer :: inter
    integer        :: np

    ! Declare local variables
    integer :: Si,Sj,i
    real(RK):: RXi,RYi,RZi,RXj,RYj,RZj
    real(RK):: drxij,dryij,drzij,dr
    real(RK):: approx
    real(RK):: UIntraTermKomp, USelbstTermKomp


#if DEBUG_4Part > 0
! DEBUG STEPHAN
    this%Component(1)%Molecule%SiteCharge(1)%RX(1) = 0.50412 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(1) = 0.81073 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(1) = 0.18949 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(1) = 0.50226 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(1) = 0.81178 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(1) = 0.18692 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(2) = 2.09938 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(2) = 1.86064 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(2) = 1.31323 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(2) = 2.09770 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(2) = 1.85859 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(2) = 1.31121 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(3) = 1.31761 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(3) = 2.08757 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(3) = 0.44842 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(3) = 1.31550 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(3) = 2.08538 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(3) = 0.44707 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(4) = 1.48200 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(4) = 0.77144 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(4) = 1.55693 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(4) = 1.48458 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(4) = 0.77169 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(4) = 1.55483 / this%BoxLength
#endif

! Selbstterm
    this%USelbstTerm = 0.0
    DO Sj=1,this%NComponents,1
       USelbstTermKomp = 0.0
       DO Si=1,this%Component(Sj)%Molecule%NCharge,1
         USelbstTermKomp = USelbstTermKomp + this%Component(Sj)%Molecule%SiteCharge(Si)%e**2
       END DO
       this%USelbstTerm = this%USelbstTerm + this%Component(Sj)%NPart * USelbstTermKomp
    END DO

    this%USelbstTerm = -this%USelbstTerm * this%Kappa / sqrt(Pi) / NProcs


! intramolecular term
    this%UIntra = 0._RK
    DO i=1,this%NComponents,1
      UIntraTermKomp = 0.0
      DO Si = 1,this%component(i)%Molecule%NCharge-1
        DO Sj = Si+1,this%component(i)%Molecule%NCharge

          RXi  = this%Component(i)%Molecule%SiteCharge(Si)%r(1)
          RYi  = this%Component(i)%Molecule%SiteCharge(Si)%r(2)
          RZi  = this%Component(i)%Molecule%SiteCharge(Si)%r(3)
          RXj  = this%Component(i)%Molecule%SiteCharge(Sj)%r(1)
          RYj  = this%Component(i)%Molecule%SiteCharge(Sj)%r(2)
          RZj  = this%Component(i)%Molecule%SiteCharge(Sj)%r(3)

          drxij = (RXi-RXj)
          dryij = (RYi-RYj)
          drzij = (RZi-RZj)

          dr = sqrt(drxij*drxij + dryij*dryij + drzij*drzij)

          if (dr .ne. 0.0) then
            call erfc_approx (this%Kappa*dr, approx)

            UIntraTermKomp = UIntraTermKomp - this%Component(i)%Molecule%SiteCharge(Si)%e* &
&                   this%Component(i)%Molecule%SiteCharge(Sj)%e / dr * (1-approx)
          end if
        END DO
      END DO
      this%UIntra = this%UIntra + this%component(i)%NPart * UIntraTermKomp / NProcs
    END DO


    end subroutine TEnsemble_EwaldSelf_Energy



!==============================================================!
!  Subroutine TSimulation_Ewald_FourierTermEnergy              !
!==============================================================!

   subroutine TEnsemble_EwaldFourierEnergy(this)

   implicit none

#if MPI_VER > 0
    include 'mpif.h'
#endif

   type(TMolecule), pointer               :: mol

   ! Declare arguments
   type(TEnsemble)   :: this

   integer :: i,j,l
   integer :: molec
# if MPI_VER > 0
   integer :: i0
   integer :: i1
!    integer :: counter
!    real(RK):: summe(NProcs,this%NPart)
# endif

   real(RK),pointer:: RX(:),RY(:),RZ(:)
   real(RK),pointer:: PX(:),PY(:),PZ(:)
   real(RK),pointer:: q(:)

   real(RK):: KVec(3)
   real(RK):: EPotLocal
   real(RK):: Viriallocal,VirIntra
   real(RK):: SSinSum,SCosSum
   real(RK):: KappaL2
   real(RK):: vorfac
   real(RK):: facx,facy,facz


! DEBUG STEPHAN
!     this%Component(1)%Molecule%SiteCharge(1)%RX(1) = 0.50412 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(1)%RY(1) = 0.81073 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(1)%RZ(1) = 0.18949 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(2)%RX(1) = 0.50226 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(2)%RY(1) = 0.81178 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(2)%RZ(1) = 0.18692 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(1)%RX(2) = 2.09938 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(1)%RY(2) = 1.86064 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(1)%RZ(2) = 1.31323 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(2)%RX(2) = 2.09770 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(2)%RY(2) = 1.85859 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(2)%RZ(2) = 1.31121 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(1)%RX(3) = 1.31761 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(1)%RY(3) = 2.08757 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(1)%RZ(3) = 0.44842 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(2)%RX(3) = 1.31550 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(2)%RY(3) = 2.08538 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(2)%RZ(3) = 0.44707 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(1)%RX(4) = 1.48200 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(1)%RY(4) = 0.77144 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(1)%RZ(4) = 1.55693 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(2)%RX(4) = 1.48458 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(2)%RY(4) = 0.77169 / this%BoxLength
!     this%Component(1)%Molecule%SiteCharge(2)%RZ(4) = 1.55483 / this%BoxLength

   EPotLocal = 0._RK
   VirialLocal = 0._RK
   KappaL2 = 1.0_RK/(2._RK*this%KappaL**2)
   vorfac = 2._RK / this%BoxLength

   this%VirIntra = 0._RK

   this%SSin = 0._RK
   this%SCos = 0._RK
# if MPI_VER > 0
!   summe = 0._RK
   j=NProc+1
   i0 = this%NBox0(j)
   i1 = this%NBox2(j)
   DO i=i0,i1,1
# else
   DO i=1,this%BoxenAnzahl,1
# endif
     KVec = this%Ewald_Vec(:,i)
     this%SSin_Vec = 0._RK
     this%SCos_Vec = 0._RK
     DO j=1,this%NComponents,1
       mol => this%Component(j)%Molecule
       q => mol%SiteCharge(1:mol%NCharge)%e
       molec = this%Component(j)%NPart
       DO l=1,mol%NCharge
         RX => this%Component(j)%Molecule%SiteCharge(l)%RX(1:molec)
         RY => this%Component(j)%Molecule%SiteCharge(l)%RY(1:molec)
         RZ => this%Component(j)%Molecule%SiteCharge(l)%RZ(1:molec)
         PX => this%Component(j)%Molecule%SiteCharge(l)%PX(1:molec)
         PY => this%Component(j)%Molecule%SiteCharge(l)%PY(1:molec)
         PZ => this%Component(j)%Molecule%SiteCharge(l)%PZ(1:molec)

         this%distx(j,l,1:molec) = (RX - PX)*this%BoxLength
         this%disty(j,l,1:molec) = (RY - PY)*this%BoxLength
         this%distz(j,l,1:molec) = (RZ - PZ)*this%BoxLength

         this%Faktor(1:molec) = KVec(1) * RX + KVec(2)*RY + KVec(3)*RZ

         this%sinfac_s(j,l,1:molec) = sin(this%Faktor)
         this%cosfac_s(j,l,1:molec) = cos(this%Faktor)
         this%sinfac = q(l)*this%sinfac_s(j,l,1:molec)
         this%cosfac = q(l)*this%cosfac_s(j,l,1:molec)

         this%SSin_Vec = this%SSin_Vec + this%sinfac
         this%SCos_Vec = this%SCos_Vec + this%cosfac
       END DO
     END DO

     SSinSum = sum(this%SSin_Vec)
     SCosSum = sum(this%SCos_Vec)
     this%SSin(i) = SSinSum
     this%SCos(i) = SCosSum

! Forces
     Facx = KVec(1)*this%Ewald_Prefac(i)*vorfac
     Facy = KVec(2)*this%Ewald_Prefac(i)*vorfac
     Facz = KVec(3)*this%Ewald_Prefac(i)*vorfac

     DO j=1,this%NComponents,1
       mol => this%Component(j)%Molecule
       q => mol%SiteCharge(1:mol%NCharge)%e
       molec = this%Component(j)%NPart
       DO l=1,mol%NCharge
          this%HFac = q(l)*(this%sinfac_s(j,l,1:molec)*this%SCos(i) - this%cosfac_s(j,l,1:molec)*this%SSin(i))
          this%VirIntra = this%VirIntra + Facx*this%HFac*this%distx(j,l,1:molec)+&
   &         Facy*this%HFac*this%disty(j,l,1:molec) + &
   &         Facz*this%HFac*this%distz(j,l,1:molec)
! # endif
       END DO
     END DO
   END DO ! Boxenschleife


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Finish Calculation

! Energy
   this%U_fourierLocal = this%Ewald_Prefac * (this%SSin*this%SSin + this%SCos*this%SCos)
! Virial
   this%Vec2 = this%Ewald_Vec(1,:)**2 + this%Ewald_Vec(2,:)**2 + this%Ewald_Vec(3,:)**2 

#if MPI_VER > 0
   call MPI_Allreduce( sum(this%U_fourierLocal), EPotLocal, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
   call MPI_Allreduce( sum(this%U_fourierLocal) - sum(this%U_fourierLocal*KappaL2*this%Vec2), VirialLocal, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
   call MPI_Allreduce( sum(this%VirIntra), VirIntra, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
#else
   EPotLocal = sum(this%U_fourierLocal)
   VirialLocal = EPotLocal - sum(this%U_fourierLocal *KappaL2*this%Vec2)
   VirIntra = sum(this%VirIntra)
#endif
!     call MPI_Allreduce( GetEnergy( this ), this%EPot, 1 , &
! &     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )

!   if( RootProc ) then
   this%UFourier= EPotLocal / NProcs
   this%EVirial = -(Viriallocal - VirIntra)*Third / NProcs
!   end if


  END subroutine TEnsemble_EwaldFourierEnergy


!==============================================================!
!  Subroutine TSimulation_Ewald_FourierTermEnergy1             !
!==============================================================!
   subroutine TEnsemble_EwaldFourierEnergy1(this,nc,np)

   implicit none

#if MPI_VER > 0
    include 'mpif.h'
#endif

   ! Declare arguments
   type(TEnsemble)          :: this
   type(TMolecule), pointer :: mol

   real(RK):: RX,RY,RZ
   real(RK),pointer:: q(:)
   real(RK):: KVec(3)
   real(RK):: EPotLocal
   real(RK):: SSinSum,SCosSum
   real(RK)::SSin_Vec,SCos_Vec
   real(RK):: KappaL2
   real(RK):: sinfac,cosfac
   real(RK)::Faktor,Faktor2

   integer :: i,j,l
   integer,intent(in)::nc,np

# if MPI_VER > 0
   integer :: i0
   integer :: i1
# endif


! Declarations
   KappaL2 = 1.0_RK/(2._RK*this%KappaL**2)

! Calculation
# if MPI_VER > 0
   j=NProc+1
   i0 = this%NBox0(j)
   i1 = this%NBox2(j)
   DO i=i0,i1,1
# else
   DO i=1,this%BoxenAnzahl,1
# endif

     KVec = this%Ewald_Vec(:,i)
     mol => this%Component(nc)%Molecule
     q => mol%SiteCharge(1:mol%NCharge)%e
!     this%sinfac_s_old = this%sinfac_s(nc,1:mol%NCharge,np)
!     this%cosfac_s_old = this%cosfac_s(nc,1:mol%NCharge,np)
     SSin_Vec =0._RK
     SCos_Vec =0._RK
       DO l=1,mol%NCharge
         RX = this%Component(nc)%Molecule%SiteCharge(l)%RX(np)
         RY = this%Component(nc)%Molecule%SiteCharge(l)%RY(np)
         RZ = this%Component(nc)%Molecule%SiteCharge(l)%RZ(np)

         Faktor = KVec(1) * RX + KVec(2)*RY + KVec(3)*RZ
         Faktor2 = KVec(1) * this%rold(l,1) + KVec(2)*this%rold(l,2) + KVec(3)*this%rold(l,3)

!         this%sinfac_s(nc,l,np) = sin(Faktor)
!         this%cosfac_s(nc,l,np) = cos(Faktor)
!         sinfac = q(l)*(this%sinfac_s(nc,l,np)-this%sinfac_s_old(l))
!         cosfac = q(l)*(this%cosfac_s(nc,l,np)-this%cosfac_s_old(l))
         sinfac = q(l)*(sin(Faktor)-sin(Faktor2))
         cosfac = q(l)*(cos(Faktor)-cos(Faktor2))

         SSin_Vec = SSin_Vec + sinfac
         SCos_Vec = SCos_Vec + cosfac
       END DO

!     SSinSum = (SSin_Vec)
!     SCosSum = (SCos_Vec)
     this%SSin(i) = this%SSin(i) + SSin_Vec
     this%SCos(i) = this%SCos(i) + SCos_Vec

   END DO ! Boxenschleife


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Finish Calculation

! Energy
   this%U_fourierLocal = this%Ewald_Prefac * (this%SSin*this%SSin + this%SCos*this%SCos)
! Virial
   this%Vec2 = this%Ewald_Vec(1,:)**2 + this%Ewald_Vec(2,:)**2 + this%Ewald_Vec(3,:)**2 

#if MPI_VER > 0
   call MPI_Allreduce( sum(this%U_fourierLocal), EPotLocal, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
#else
   EPotLocal = sum(this%U_fourierLocal)
#endif

!   if( RootProc ) then
   this%UFourier= EPotLocal  / NProcs
!   end if


  END subroutine TEnsemble_EwaldFourierEnergy1




!==============================================================!
!  Subroutine TSimulation_Ewald_FourierTermEnergy1             !
!==============================================================!
   subroutine TEnsemble_EwaldFourierEnergy_CF(this,nc,np,ncold,npold)

   implicit none

#if MPI_VER > 0
    include 'mpif.h'
#endif

   ! Declare arguments
   type(TEnsemble)          :: this
   type(TMolecule), pointer :: mol
   type(TMolecule), pointer :: mol2

   real(RK):: RX,RY,RZ
   real(RK),pointer:: q(:)
   real(RK):: KVec(3)
   real(RK):: EPotLocal
   real(RK):: SSinSum,SCosSum
   real(RK)::SSin_Vec,SCos_Vec
   real(RK):: KappaL2
   real(RK):: sinfac,cosfac
   real(RK)::Faktor,Faktor2

   integer :: i,j,l
   integer,intent(in)::nc,np
   integer,intent(in)::ncold,npold

# if MPI_VER > 0
   integer :: i0
   integer :: i1
# endif

! Declarations
   KappaL2 = 1.0_RK/(2._RK*this%KappaL**2)

! Calculation
# if MPI_VER > 0
   j=NProc+1
   i0 = this%NBox0(j)
   i1 = this%NBox2(j)
   DO i=i0,i1,1
# else
   DO i=1,this%BoxenAnzahl,1
# endif
     KVec = this%Ewald_Vec(:,i)
     mol => this%Component(nc)%Molecule
     mol2 => this%Component(ncold)%Molecule
     q => mol%SiteCharge(1:mol%NCharge)%e
     SSin_Vec =0._RK
     SCos_Vec =0._RK
       DO l=1,mol%NCharge
         RX = this%Component(nc)%Molecule%SiteCharge(l)%RX(np)
         RY = this%Component(nc)%Molecule%SiteCharge(l)%RY(np)
         RZ = this%Component(nc)%Molecule%SiteCharge(l)%RZ(np)

         Faktor = KVec(1) * RX + KVec(2)*RY + KVec(3)*RZ

         sinfac = q(l)*sin(Faktor)
         cosfac = q(l)*cos(Faktor)

         SSin_Vec = SSin_Vec + sinfac
         SCos_Vec = SCos_Vec + cosfac
       END DO
       DO l=1,mol2%NCharge
         RX = this%rold(l,1)
         RY = this%rold(l,2)
         RZ = this%rold(l,3)

         Faktor2 = KVec(1) * this%rold(l,1) + KVec(2)*this%rold(l,2) + KVec(3)*this%rold(l,3)

         sinfac = q(l)*sin(Faktor2)
         cosfac = q(l)*cos(Faktor2)

         SSin_Vec = SSin_Vec - sinfac
         SCos_Vec = SCos_Vec - cosfac
       END DO

     this%SSin(i) = this%SSin(i) + SSin_Vec
     this%SCos(i) = this%SCos(i) + SCos_Vec

   END DO ! Boxenschleife

! Energy
   this%U_fourierLocal = this%Ewald_Prefac * (this%SSin*this%SSin + this%SCos*this%SCos)

#if MPI_VER > 0
   call MPI_Allreduce( sum(this%U_fourierLocal), EPotLocal, 1, &
&     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror )
#else
   EPotLocal = sum(this%U_fourierLocal)
#endif
   this%UFourier= EPotLocal  / NProcs

  END subroutine TEnsemble_EwaldFourierEnergy_CF






!==============================================================!
!  Subroutine TEnsemble_Ewald_ChemicalPotential                !
!==============================================================!

  subroutine TEnsemble_Ewald_ChemPotSelf( this,nc1)

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)            :: this
    type(TMolecule), pointer   :: mol
    integer, intent(in)        :: nc1

    ! Declare local variables
    integer :: Si,Sj, np
    real(RK):: fac
    real(RK):: UTermKomp1

    fac = 1.0_RK / this%RCutoffLJ126LJ126
! Selfterm
    UTermKomp1 = 0._RK
    DO Si=1,this%Component(nc1)%Molecule%NCharge,1
         UTermKomp1 = UTermKomp1 - &
&                  this%Component(nc1)%Molecule%SiteCharge(Si)%e**2
    END DO


! Final Summation
    this%Component(nc1)%EPotTestSelf = UTermKomp1*fac
    this%Component(nc1)%EPotTestSelf = 0._RK

    end subroutine TEnsemble_Ewald_ChemPotSelf



! ! ! ! 
! ! ! ! 
! ! ! ! 
! ! ! ! 
! ! ! ! 
! ! ! !   subroutine TEnsemble_Ewald_ChemPotFour( this,nc1)
! ! ! ! 
! ! ! !     implicit none
! ! ! ! 
! ! ! !     ! Include MPI header
! ! ! ! #if MPI_VER > 0
! ! ! !     include 'mpif.h'
! ! ! ! #endif
! ! ! ! 
! ! ! !     ! Declare arguments
! ! ! !     type(TEnsemble)            :: this
! ! ! !     type(TMolecule), pointer   :: mol
! ! ! !     integer, intent(in)        :: nc1
! ! ! ! 
! ! ! !     ! Declare local variables
! ! ! !     integer :: j,i,l
! ! ! !     real(RK):: RXi,RYi,RZi
! ! ! !     real(RK):: KVec(3)
! ! ! !     real(RK):: SSinSum, SCosSum
! ! ! !     real(RK):: SSinTest, SCosTest
! ! ! !     real(RK),pointer :: q(:)
! ! ! !     real(RK):: Faktor
! ! ! ! 
! ! ! ! 
! ! ! ! 
! ! ! !    DO j=1,this%Component(nc1)%NTest
! ! ! ! 
! ! ! ! 
! ! ! ! ! #if MPI_VER > 0
! ! ! ! !      l=NProc + 1
! ! ! ! !      i0 = this%NBox0(l)
! ! ! ! !      i1 = this%NBox2(l)
! ! ! ! ! !    i0 = 1
! ! ! ! !      DO i=i0,i1,1
! ! ! ! ! #else
! ! ! !      DO i=1,this%BoxenAnzahl,1
! ! ! ! ! # endif
! ! ! !        SSinSum  = 0.0_RK
! ! ! !        SCosSum  = 0.0_RK
! ! ! ! 
! ! ! !        SSinTest = this%SSin(i)
! ! ! !        SCosTest = this%SCos(i)
! ! ! ! 
! ! ! !        KVec = this%Ewald_Vec(:,i)
! ! ! !        mol => this%Component(nc1)%Molecule
! ! ! !        q => mol%SiteCharge(1:mol%NCharge)%e
! ! ! ! !        molec = this%Component(j)%NPart
! ! ! !        DO l=1,mol%NCharge
! ! ! !          RXi = this%Component(nc1)%Molecule%SiteCharge(l)%RXTest(j)
! ! ! !          RYi = this%Component(nc1)%Molecule%SiteCharge(l)%RYTest(j)
! ! ! !          RZi = this%Component(nc1)%Molecule%SiteCharge(l)%RZTest(j)
! ! ! ! 
! ! ! !          Faktor = KVec(1) * RXi + KVec(2)*RYi + KVec(3)*RZi
! ! ! ! 
! ! ! !          SSinSum  = SSinSum + q(l) * sin(Faktor)
! ! ! !          SCosSum  = SCosSum + q(l) * cos(Faktor)
! ! ! ! 
! ! ! !          SSinTest = SSinTest + SSinSum 
! ! ! !          SCosTest = SCosTest + SCosSum
! ! ! !        END DO
! ! ! ! 
! ! ! ! ! Subtract the image of just the test particle immediately *Sum*-terms!
! ! ! !        this%EPotTest(j) = this%EPotTest(j) + this%Ewald_Prefac(i) * &
! ! ! ! &              (SSinTest*SSinTest + SCosTest*SCosTest - SSinSum*SSinSum - &
! ! ! ! &               SCosSum*SCosSum)
! ! ! ! 
! ! ! !      END DO ! Boxenschleife
! ! ! ! 
! ! ! !    END DO
! ! ! ! 
! ! ! ! 
! ! ! ! ! Finish Calculation
! ! ! ! ! Energy
! ! ! ! 
! ! ! ! ! #if MPI_VER > 0
! ! ! ! !    call MPI_Reduce( sum(this%U_fourierLocal), EPotLocal, 1, &
! ! ! ! ! &     MPI_DOUBLE_PRECISION, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
! ! ! ! ! #else
! ! ! ! !    EPotLocal = sum(this%U_fourierLocal)
! ! ! ! ! #endif
! ! ! ! !   if( RootProc ) then
! ! ! ! !    this%UFourier= EPotLocal
! ! ! ! !   end if
! ! ! ! 
! ! ! !     end subroutine TEnsemble_Ewald_ChemPotFour
! ! ! ! 






!==============================================================!
!  Subroutine TEnsemble_Constraints                            !
!  Calculate part of SHAKE                                     !
!==============================================================!
#if CONSTR > 0

  subroutine TEnsemble_Constraints( this)

    implicit none

   type(TEnsemble) :: this

   integer         :: maxit
   integer         :: i, j, aa, bb
   integer         :: aacomp, bbcomp
!    real(RK),pointer:: PX(:),PY(:),PZ(:)
   real(RK)        :: PX1, PY1, PZ1
   real(RK)        :: PX2, PY2, PZ2
   real(RK)        :: dx,dy,dz
   real(RK)        :: ddx,ddy,ddz
   real(RK)        :: dist2, dr2, dist
   real(RK)        :: fac
   real(RK)        :: Forc
   real(RK)        :: dLOgVolumeThird
   real(RK)        :: tol

   logical         :: cont

! Initialization of important variables
   cont  = .false.
! Initialization of the max number of iterations for SHAKE
   maxit = 300
   tol   = 1e-7
   dLogVolumeThird = this%Volume1 / (3._RK * this%Volume0)


   if (this%consup .eq. .true.) then
     DO j=1,this%NCons,1
          write( IOBuffer, '(F10.5)' ) this%UCons(j) / BlockSize
          call FileWriteNoAdvance( this%iounit_runave )
          write( IOBuffer, '(F10.5)' ) this%FCons(j) / BlockSize
          call FileWriteNoAdvance( this%iounit_runave )

          this%FCons(j) = 0._RK
          this%UCons(j) = 0._RK
     END DO

     call FileWriteBlank( this%iounit_runave )
#if ARCH == 2
     call flush( this%iounit_runave )
#endif
     this%consup = .false.
   end if
!     WHILE ((i<maxit) .AND. (cont .eq. .false.)) DO
    i=0
    DO WHILE ((i .le. maxit) .AND. (cont .eq. .false.))
       cont  = .true.
       DO j=1,this%NCons,1
         aacomp  = this%Cons1Comp(j)
         bbcomp  = this%Cons2Comp(j)
         aa      = this%Cons1(j)
         bb      = this%Cons2(j)
         dr2     = this%ConsR(j)

         PX1 = this%Component(aacomp)%P0(aa,1)
         PY1 = this%Component(aacomp)%P0(aa,2)
         PZ1 = this%Component(aacomp)%P0(aa,3)
         PX2 = this%Component(bbcomp)%P0(bb,1)
         PY2 = this%Component(bbcomp)%P0(bb,2)
         PZ2 = this%Component(bbcomp)%P0(bb,3)

         dx  = (PX2 - PX1)
         dy  = (PY2 - PY1)
         dz  = (PZ2 - PZ1)

         dx  = (dx - anint(dx))*this%BoxLength
         dy  = (dy - anint(dy))*this%BoxLength
         dz  = (dz - anint(dz))*this%BoxLength

         dist2 = dx*dx + dy*dy + dz*dz

         dist  = dist2 - dr2

         if (abs(dist) .gt. tol) then
            Forc = 0._RK
            cont = .false.
            fac = 1 - sqrt(dr2 / dist2)

            ddx = 0.5*fac * dx/this%BoxLength
            ddy = 0.5*fac * dy/this%BoxLength
            ddz = 0.5*fac * dz/this%BoxLength

            call CorrectGear_Constraint(this%Component(aacomp),aa,dLogVolumeThird,&
&                                        Forc, ddx,ddy,ddz)
            call CorrectGear_Constraint(this%Component(bbcomp),bb,dLogVolumeThird,&
&                                        Forc,-ddx,-ddy,-ddz)
         end if
         this%FCons(j) = this%FCons(j) + Forc
         this%UCons(j) = this%UCons(j) + Forc*dist

!         PX1 = this%Component(aacomp)%P0(aa,1)
!         PY1 = this%Component(aacomp)%P0(aa,2)
!         PZ1 = this%Component(aacomp)%P0(aa,3)
!         PX2 = this%Component(bbcomp)%P0(bb,1)
!         PY2 = this%Component(bbcomp)%P0(bb,2)
!         PZ2 = this%Component(bbcomp)%P0(bb,3)

!         dx  = (PX2 - PX1)*this%BoxLength
!         dy  = (PY2 - PY1)*this%BoxLength
!         dz  = (PZ2 - PZ1)*this%BoxLength

!         dist2 = dx*dx + dy*dy + dz*dz

!         dist  = dist2 - dr2

       END DO

!        if (cont = .true.) then
!           goto 1000
!        end if
    END DO

! 1000 

  end subroutine TEnsemble_Constraints

#endif
!==============================================================!
!  Subroutine TSimulation_Ewald_FourierTermEnergy              !
!==============================================================!
! 
!    subroutine TEnsemble_EwaldFourierEnergy(this, nc1,np)
! 
!    implicit none
! 
!    integer :: i,j,k,l
!    real(RK):: Faktor
!    real(RK):: fac1,fac2,fac3
!    real(RK):: KVec(3)
!    real(RK):: RX,RY,RZ
!    real(RK):: kappasqrt4, summe, KVec_LenInv
!    real(RK):: UTermKomp, UTerm, VirialLocal
!    real(RK):: complr(this%NMax),compll(this%NMax)
! 
!    ! Declare arguments
!    type(TEnsemble)                 :: this
!    type(TInteraction),pointer      :: inter
!    type(TMolecule), pointer        :: mol
!    integer, intent(in) :: nc1
!    integer, intent(in) :: np
! 
! ! initialization
!    summe = 0._RK
!    UTerm = 0._RK
!    VirialLocal = 0._RK
!    kappasqrt4 = 1.0_RK/(4*this%kappa**2)
! 
!    inter => this%Interaction(nc1,nc1)
! 
! 
! ! Preparation
! 
!   complr(1) = 1._RK
!   compll(1) = 0._RK
!   complr(2) = twopi * RX
!    DO i=2,this%NMax,1
!       compl(i) = compl(i-1) * compl(1)
!    DO i=1,this%BoxenAnzahl,1
! ! Preparation
!      Fac3  = 0._RK
!      summe=0._RK
!      KVec(1:3) = Pi2 * this%Ewald_Vec(1:3,i) / this%BoxLength
! 
! ! Virials
!      KVec_LenInv = 1._RK / (KVec(1)**2 + KVec(2)**2 + KVec(3)**2)
!      DO j=1,3,1
!          summe = summe + KVec(j)*KVec(j)
!      END DO
! 
! 
!      mol => this%Component(nc1)%Molecule
! 
!      DO l=1,mol%NCharge
!        RX = mol%SiteCharge(l)%RX(np)*this%BoxLength
!        RY = mol%SiteCharge(l)%RY(np)*this%BoxLength
!        RZ = mol%SiteCharge(l)%RZ(np)*this%BoxLength
! 
!        Faktor = KVec(1)*RX + KVec(2)*RY + KVec(3)*RZ
! 
!        Fac1 = (mol%SiteCharge(l)%e*sin(Faktor))**2
!        Fac2 = (mol%SiteCharge(l)%e*cos(Faktor))**2
!        Fac3 = Fac3 + (Fac1+Fac2)
! !        DO k=1,this%Component(nc2)%NPart,1
! !          inter%Epot1(k) = inter%Epot1(k) + Fac3
! !        END DO
!      END DO
!      UTermKomp = Fac3 * this%Ewald_Prefac(i)
!      UTerm = UTerm + UTermKomp
!      VirialLocal = VirialLocal + UTermKomp * (3-2*(KVec_LenInv+kappasqrt4)*summe)
! 
!    END DO      ! end loop over all Ewald boxes k
!    inter%Epot1(np) = inter%Epot1(np) +  UTerm / this%BoxLength
!    inter%Virial1(np) = inter%Virial1(np) &
! &                               + Third * VirialLocal / (this%BoxLength)
! 
!    end subroutine TEnsemble_EwaldFourierEnergy
! 

! !==============================================================!
! !  Subroutine TSimulation_Ewald_FourierTermEnergy              !
! !==============================================================!
! 
!    subroutine TEnsemble_EwaldFourierEnergy(this, nc1,np)
! 
!    implicit none
! 
!    integer :: i,j,k,l
!    real(RK):: Faktor
!    real(RK):: fac1,fac2,fac3
!    real(RK):: KVec(3)
!    real(RK):: RX,RY,RZ
!    real(RK):: kappasqrt4, summe, KVec_LenInv
!    real(RK):: UTermKomp, UTerm, VirialLocal
! 
!    ! Declare arguments
!    type(TEnsemble)                 :: this
!    type(TInteraction),pointer      :: inter
!    type(TMolecule), pointer        :: mol
!    integer, intent(in) :: nc1
!    integer, intent(in) :: np
! 
! ! initialization
!    summe = 0._RK
!    UTerm = 0._RK
!    VirialLocal = 0._RK
!    kappasqrt4 = 1.0_RK/(4*this%kappa**2)
! 
!    inter => this%Interaction(nc1,nc1)
!    DO i=1,this%BoxenAnzahl,1
! ! Preparation
!      Fac3  = 0._RK
!      summe=0._RK
!      KVec(1:3) = Pi2 * this%Ewald_Vec(1:3,i) / this%BoxLength
! 
! ! Virials
!      KVec_LenInv = 1._RK / (KVec(1)**2 + KVec(2)**2 + KVec(3)**2)
!      DO j=1,3,1
!          summe = summe + KVec(j)*KVec(j)
!      END DO
! 
! 
!      mol => this%Component(nc1)%Molecule
! 
!      DO l=1,mol%NCharge
!        RX = mol%SiteCharge(l)%RX(np)*this%BoxLength
!        RY = mol%SiteCharge(l)%RY(np)*this%BoxLength
!        RZ = mol%SiteCharge(l)%RZ(np)*this%BoxLength
! 
!        Faktor = KVec(1)*RX + KVec(2)*RY + KVec(3)*RZ
! 
!        Fac1 = (mol%SiteCharge(l)%e*sin(Faktor))**2
!        Fac2 = (mol%SiteCharge(l)%e*cos(Faktor))**2
!        Fac3 = Fac3 + (Fac1+Fac2)
! !        DO k=1,this%Component(nc2)%NPart,1
! !          inter%Epot1(k) = inter%Epot1(k) + Fac3
! !        END DO
!      END DO
!      UTermKomp = Fac3 * this%Ewald_Prefac(i)
!      UTerm = UTerm + UTermKomp
!      VirialLocal = VirialLocal + UTermKomp * (3-2*(KVec_LenInv+kappasqrt4)*summe)
! 
!    END DO      ! end loop over all Ewald boxes k
!    inter%Epot1(np) = inter%Epot1(np) +  UTerm / this%BoxLength
!    inter%Virial1(np) = inter%Virial1(np) &
! &                               + Third * VirialLocal / (this%BoxLength)
! 
! 
!    end subroutine TEnsemble_EwaldFourierEnergy








!==============================================================!
!  Subroutine TSimulation_Ewald_ChemPotFourier                 !
!==============================================================!

!    subroutine TSimulation_EwaldFourier_ChemicalPotential(this, nc1)
! 
!    implicit none
! 
!    integer :: i,l,np
!    real(RK):: Faktor
!    real(RK):: fac1,fac2,fac3
!    real(RK):: KVec(3)
!    real(RK):: RX,RY,RZ
!    real(RK):: UTermKomp
! 
!    ! Declare arguments
!    type(TEnsemble)             :: this
!    type(TInteraction),pointer  :: inter
!    type(TMolecule), pointer    :: mol
!    integer, intent(in) :: nc1
! 
! 
! ! initialization
!    UTermKomp = 0._RK
! 
!    DO i=1,this%BoxenAnzahl,1
!      KVec(1:3) = Pi2 * this%Ewald_Vec(1:3,i) / this%BoxLength
! 
!      mol => this%Component(nc1)%Molecule
! 
! 
!      DO np=1,this%Component(nc1)%NPart
!        Fac3 = 0._RK
!        DO l=1,mol%NCharge
!          RX = mol%SiteCharge(l)%RX(np)*this%BoxLength
!          RY = mol%SiteCharge(l)%RY(np)*this%BoxLength
!          RZ = mol%SiteCharge(l)%RZ(np)*this%BoxLength
! 
!          Faktor = KVec(1)*RX + KVec(2)*RY + KVec(3)*RZ
! 
!          Fac1 = (mol%SiteCharge(l)%e*sin(Faktor))**2
!          Fac2 = (mol%SiteCharge(l)%e*cos(Faktor))**2
!          Fac3 = Fac3 + (Fac1+Fac2)
!        END DO
!        this%EPotTest(np) = this%EPotTest(np) + Fac3*this%Ewald_Prefac(i) / this%BoxLength
!      END DO
! 
! 
!    END DO      ! end loop over all Ewald boxes k
! 
!    end subroutine TSimulation_EwaldFourier_ChemicalPotential





! ! ! 
! ! ! 
! ! ! 
! ! ! 
! ! ! !==============================================================!
! ! ! !  Subroutine TSimulation_Ewald_FourierTerm                    !
! ! ! !==============================================================!
! ! ! 
! ! !    subroutine TEnsemble_EwaldFourierTerm(this)
! ! ! 
! ! !    implicit none
! ! ! 
! ! !    integer :: i,j,k,l, counter
! ! !    real(RK):: Faktor, Faktor2, Faktor3(3)
! ! !    real(RK):: U_fourier
! ! !    real(RK):: KVec(3)
! ! !    real(RK):: SSin,SCos
! ! !    real(RK):: RX,RY,RZ
! ! !    real(RK):: FX,FY,FZ
! ! !    real(RK):: EPotLocal
! ! !    real(RK):: Viriallocal
! ! !    real(RK):: test
! ! !    real(RK):: summe, KappaLsqrt4, KVec_LenInv
! ! ! !    real(RK),pointer:: SSin_Fac, SCos_fac
! ! ! 
! ! !    type(TMolecule), pointer               :: mol
! ! ! 
! ! !    ! Declare arguments
! ! !    type(TEnsemble)   :: this
! ! ! 
! ! ! !    real(RK):: sin_test(this%BoxenAnzahl,this%Component(1)%NPart,this%Component(1)%molecule%NCharge)
! ! ! !    real(RK):: cos_test(this%BoxenAnzahl,this%Component(1)%NPart,this%Component(1)%molecule%NCharge)
! ! ! 
! ! ! ! DEBUG STEPHAN
! ! ! !     this%Component(1)%Molecule%SiteCharge(1)%RX(1) = 0.50412 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(1)%RY(1) = 0.81073 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(1)%RZ(1) = 0.18949 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(2)%RX(1) = 0.50226 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(2)%RY(1) = 0.81178 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(2)%RZ(1) = 0.18692 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(1)%RX(2) = 2.09938 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(1)%RY(2) = 1.86064 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(1)%RZ(2) = 1.31323 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(2)%RX(2) = 2.09770 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(2)%RY(2) = 1.85859 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(2)%RZ(2) = 1.31121 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(1)%RX(3) = 1.31761 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(1)%RY(3) = 2.08757 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(1)%RZ(3) = 0.44842 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(2)%RX(3) = 1.31550 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(2)%RY(3) = 2.08538 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(2)%RZ(3) = 0.44707 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(1)%RX(4) = 1.48200 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(1)%RY(4) = 0.77144 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(1)%RZ(4) = 1.55693 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(2)%RX(4) = 1.48458 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(2)%RY(4) = 0.77169 / this%BoxLength
! ! ! !     this%Component(1)%Molecule%SiteCharge(2)%RZ(4) = 1.55483 / this%BoxLength
! ! ! 
! ! !    EPotLocal = 0.0_RK
! ! !    VirialLocal = 0.0_RK
! ! !    KappaLsqrt4 = 1.0_RK/(4*this%KappaL**2)
! ! ! !    this%BoxLength = 3.6445654373
! ! ! 
! ! !    DO i=1,this%BoxenAnzahl,1
! ! !      KVec(1:3) = this%Ewald_Vec(1:3,i)
! ! ! !      KVec_LenInv = 1._RK / (KVec(1)**2 + KVec(2)**2 + KVec(3)**2)
! ! !      summe = KVec(1)*KVec(1) + KVec(2)*KVec(2) + KVec(3)*KVec(3)
! ! !      SSin = 0.0_RK
! ! !      SCos = 0.0_RK
! ! !      counter = 1
! ! ! 
! ! !      DO j=1,this%NComponents,1
! ! !        mol => this%Component(j)%Molecule
! ! !        DO k=1,this%Component(j)%NPart,1
! ! ! 
! ! !          DO l=1,mol%NCharge
! ! !            RX = mol%SiteCharge(l)%RX(k)
! ! !            RY = mol%SiteCharge(l)%RY(k)
! ! !            RZ = mol%SiteCharge(l)%RZ(k)
! ! ! 
! ! !            Faktor = KVec(1)*RX + KVec(2)*RY + KVec(3)*RZ
! ! ! 
! ! !            this%SSin_Fac(counter) = sin(Faktor)
! ! !            this%SCos_Fac(counter) = cos(Faktor)
! ! !            SSin = SSin + mol%SiteCharge(l)%e*this%SSin_Fac(counter)
! ! !            SCos = SCos + mol%SiteCharge(l)%e*this%SCos_Fac(counter)
! ! !            counter = counter + 1
! ! ! !            sin_test(i,k,l) = sin(Faktor)
! ! ! !            cos_test(i,k,l) = cos(Faktor)
! ! !          END DO
! ! !        END DO
! ! !      END DO
! ! ! 
! ! !      U_fourier = this%Ewald_Prefac(i)*(SSin*SSin + SCos*SCos)
! ! ! !      test = SSin*SSin + SCos*SCos
! ! !      EPotLocal = EPotLocal + U_fourier
! ! ! 
! ! ! ! Virial just Trace(Virial) of interest
! ! ! !      summe =0.0_RK
! ! ! !      DO j=1,3,1
! ! ! !          summe = summe + KVec(j)*KVec(j)
! ! ! !      END DO
! ! ! !      VirialLocal = VirialLocal + U_fourier * (3-2*(KVec_LenInv+KappaLsqrt4)*summe)
! ! !      VirialLocal = VirialLocal + U_fourier * (1-2*KappaLsqrt4*summe)
! ! ! !
! ! ! !
! ! ! ! Calculation of Forces
! ! !      counter = 1
! ! !      Faktor3(1:3) = 2._RK * this%Ewald_Prefac(i) / this%BoxLength * KVec(1:3)
! ! !      DO j=1,this%NComponents,1
! ! !        mol => this%Component(j)%Molecule
! ! !        DO k=1,this%Component(j)%NPart,1
! ! ! 
! ! !          DO l=1,mol%NCharge
! ! ! !            RX = mol%SiteCharge(l)%RX(k)
! ! ! !            RY = mol%SiteCharge(l)%RY(k)
! ! ! !            RZ = mol%SiteCharge(l)%RZ(k)
! ! ! ! 
! ! ! !            Faktor = KVec(1)*RX + KVec(2)*RY + KVec(3)*RZ
! ! ! ! Calculation of Forces
! ! ! !            Faktor2 = 2*mol%SiteCharge(l)%e*this%Ewald_Prefac(i)*&
! ! ! ! &                    (this%SSin_Fac(counter)*SCos - this%SCos_Fac(counter)*SSin)
! ! !            Faktor2 = mol%SiteCharge(l)%e*&
! ! ! &                    (this%SSin_Fac(counter)*SCos - this%SCos_Fac(counter)*SSin)
! ! ! 
! ! ! !            test = sin(Faktor)*SCos - cos(Faktor)*SSin
! ! ! 
! ! !            FX = Faktor2 * Faktor3(1)
! ! !            FY = Faktor2 * Faktor3(2)
! ! !            FZ = Faktor2 * Faktor3(3)
! ! !            mol%SiteCharge(l)%FX(k) = mol%SiteCharge(l)%FX(k) + FX
! ! !            mol%SiteCharge(l)%FY(k) = mol%SiteCharge(l)%FY(k) + FY
! ! !            mol%SiteCharge(l)%FZ(k) = mol%SiteCharge(l)%FZ(k) + FZ
! ! ! 
! ! !            counter = counter + 1
! ! !          END DO
! ! !        END DO
! ! !      END DO
! ! !    END DO      ! end loop over all Ewald boxes k
! ! ! 
! ! !    this%Epot   = this%Epot   + EPotLocal
! ! !    this%Virial = this%Virial + Viriallocal*Third
! ! ! 
! ! !    end subroutine TEnsemble_EwaldFourierTerm
! ! ! 
! ! ! 
! ! ! 














!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!















!==============================================================!
!  Subroutine TEnsemble_PME_SelfTerm                       !
!==============================================================!

!   subroutine TEnsemble_PMESelfTerm( this )
! 
!     implicit none
! 
!     ! Include MPI header
! #if MPI_VER > 0
!     include 'mpif.h'
! #endif
!     include 'fftw3.f'
!     ! Declare arguments
!     type(TEnsemble)         :: this
!     ! Declare local variables
!     integer :: i,Si,Sj
!     integer :: NX,NY,NZ
!     real(RK):: USelbstTermKomp, VirialLocal
! !     real(RK):: ZweiPi_inv
!     real(RK):: Faktor, NSQ
!     real(RK):: approx
!     real(RK):: twopi
!     real(RK):: err
!     real(RK)::bspx(this%splineorder)
!     real(RK)::bspy(this%splineorder)
!     real(RK)::bspz(this%splineorder)
! 
!     real(RK):: UIntraTermKomp
!     real(RK):: dr, drxij,dryij,drzij
!     real(RK):: RXi,RYi,RZi,RXj,RYj,RZj
!     real(RK):: q1,q2
! 
! ! DEBUG REASONS
! !     this%Component(1)%Molecule%SiteCharge(1)%RX(1) = 0.5041221762148172
! !     this%Component(1)%Molecule%SiteCharge(1)%RY(1) = 0.8107396528595158
! !     this%Component(1)%Molecule%SiteCharge(1)%RZ(1) = 0.1894949962472363
! !     this%Component(1)%Molecule%SiteCharge(2)%RX(1) = 0.502268554151700
! !     this%Component(1)%Molecule%SiteCharge(2)%RY(1) = 0.8117808162647254
! !     this%Component(1)%Molecule%SiteCharge(2)%RZ(1) = 0.186927667845592
! !     this%Component(1)%Molecule%SiteCharge(1)%RX(2) = 2.099380236705095
! !     this%Component(1)%Molecule%SiteCharge(1)%RY(2) = 1.860649777419017
! !     this%Component(1)%Molecule%SiteCharge(1)%RZ(2) = 1.313235759816588
! !     this%Component(1)%Molecule%SiteCharge(2)%RX(2) = 2.097708401492369
! !     this%Component(1)%Molecule%SiteCharge(2)%RY(2) = 1.858592437418999
! !     this%Component(1)%Molecule%SiteCharge(2)%RZ(2) = 1.311215009892196
! !     this%Component(1)%Molecule%SiteCharge(1)%RX(3) = 1.317611713633559
! !     this%Component(1)%Molecule%SiteCharge(1)%RY(3) = 2.087579656401032
! !     this%Component(1)%Molecule%SiteCharge(1)%RZ(3) = 0.4484239832120850
! !     this%Component(1)%Molecule%SiteCharge(2)%RX(3) = 1.315500144816407
! !     this%Component(1)%Molecule%SiteCharge(2)%RY(3) = 2.085384351686846
! !     this%Component(1)%Molecule%SiteCharge(2)%RZ(3) = 0.4470700905082808
! !     this%Component(1)%Molecule%SiteCharge(1)%RX(4) = 1.482005925799646
! !     this%Component(1)%Molecule%SiteCharge(1)%RY(4) = 0.7714424297398743
! !     this%Component(1)%Molecule%SiteCharge(1)%RZ(4) = 1.556938840875683
! !     this%Component(1)%Molecule%SiteCharge(2)%RX(4) = 1.484582190214423
! !     this%Component(1)%Molecule%SiteCharge(2)%RY(4) = 0.7716945714654903
! !     this%Component(1)%Molecule%SiteCharge(2)%RZ(4) = 1.554838746276274
! ! 
! ! 
! ! 
! ! 
! ! DEBUG STEPHAN
! !     this%Component(1)%Molecule%SiteCharge(1)%RX(1) = 0.50412 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(1)%RY(1) = 0.81073 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(1)%RZ(1) = 0.18949 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(2)%RX(1) = 0.50226 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(2)%RY(1) = 0.81178 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(2)%RZ(1) = 0.18692 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(1)%RX(2) = 2.09938 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(1)%RY(2) = 1.86064 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(1)%RZ(2) = 1.31323 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(2)%RX(2) = 2.09770 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(2)%RY(2) = 1.85859 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(2)%RZ(2) = 1.31121 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(1)%RX(3) = 1.31761 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(1)%RY(3) = 2.08757 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(1)%RZ(3) = 0.44842 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(2)%RX(3) = 1.31550 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(2)%RY(3) = 2.08538 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(2)%RZ(3) = 0.44707 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(1)%RX(4) = 1.48200 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(1)%RY(4) = 0.77144 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(1)%RZ(4) = 1.55693 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(2)%RX(4) = 1.48458 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(2)%RY(4) = 0.77169 / this%BoxLength
! !     this%Component(1)%Molecule%SiteCharge(2)%RZ(4) = 1.55483 / this%BoxLength
! 
! ! Selbstterm
!     this%USelbstTerm = 0.0
!     DO i=1,this%NComponents,1
!        USelbstTermKomp = 0.0
!        DO Si=1,this%Component(i)%Molecule%NCharge,1
!          USelbstTermKomp = USelbstTermKomp + &
! &            this%Component(i)%Molecule%SiteCharge(Si)%e**2
!        END DO
!        this%USelbstTerm = this%USelbstTerm + &
! &            this%Component(i)%NPart * USelbstTermKomp
!     END DO
!     this%USelbstTerm = -this%USelbstTerm * this%Kappa / sqrt(Pi)
! 
! 
! ! Intramolecular
!     this%UIntra = 0._RK
!     this%EVirialIntra = 0._RK
!     DO i=1,this%NComponents,1
!       UIntraTermKomp = 0.0
!       VirialLocal    = 0._RK
!       DO Si = 1,this%component(i)%Molecule%NCharge-1
!         q1 = this%Component(i)%Molecule%SiteCharge(Si)%e
!         DO Sj = Si+1,this%component(i)%Molecule%NCharge
!           q2 = this%Component(i)%Molecule%SiteCharge(Sj)%e
!           RXi  = this%Component(i)%Molecule%SiteCharge(Si)%r(1)
!           RYi  = this%Component(i)%Molecule%SiteCharge(Si)%r(2)
!           RZi  = this%Component(i)%Molecule%SiteCharge(Si)%r(3)
!           RXj  = this%Component(i)%Molecule%SiteCharge(Sj)%r(1)
!           RYj  = this%Component(i)%Molecule%SiteCharge(Sj)%r(2)
!           RZj  = this%Component(i)%Molecule%SiteCharge(Sj)%r(3)
! 
!           drxij = (RXi-RXj)
!           dryij = (RYi-RYj)
!           drzij = (RZi-RZj)
! 
!           dr = sqrt(drxij*drxij + dryij*dryij + drzij*drzij)
! 
!           call erfc_approx (this%Kappa*dr, approx)
! 
!           UIntraTermKomp = UIntraTermKomp - this%Component(i)%Molecule%SiteCharge(Si)%e* &
! &                   this%Component(i)%Molecule%SiteCharge(Sj)%e / dr * (1-approx)
! 
!           VirialLocal = VirialLocal + q1 * q2 / dr * (1-approx) -2._RK*this%Kappa&
! &                   / sqrt(pi) * exp(-(this%Kappa*dr)**2) * q1*q2
!         END DO
!       END DO
!       this%UIntra = this%UIntra + this%component(i)%NPart * UIntraTermKomp
!       this%EVirialIntra = this%EVirialIntra + VirialLocal *this%component(i)%NPart
!     END DO
! 
!     this%EVirialIntra = this%EVirialIntra * Third
! 
!     NX = this%gridx
!     NY = this%gridy
!     NZ = this%gridz
! 
!     call dfftw_plan_dft_3d(this%qgrid_forward,NX,NY,NZ,this%qgrida,this%qgrida,FFTW_FORWARD,FFTW_PATIENT)
!     call dfftw_plan_dft_3d(this%qgrid_backward,NX,NY,NZ,this%qgrida,this%qgrida,FFTW_BACKWARD,FFTW_PATIENT)
! 
!     err = pme_bspline(NX,NY,NZ)
! 
! 
! !     this%NBox1 = 1 +(this%BoxenAnzahl - 1) / NProcs
! !     DO i=1,NProcs
! !       this%NBox0(i) = 1 + this%NBox1(i) * (i-1)
! !     END DO
! !     this%NBox2 = min(this%NBox0 + this%NBox1 - 1, this%BoxenAnzahl)
! !     this%NBox1 = this%NBox2 - this%NBox0
! 
! 
! 
! contains
! 
! !!!!!!!!!!!!!!!!!!!!
!     real(RK) function pme_bspline(NX,NY,NZ)
! 
!     implicit none
!     integer,intent(in) :: NX,NY,NZ
! 
!     real(RK)    :: w
!     real(RK)    :: err
!     real(RK),dimension(this%splineorder+5):: spline,dspline
!     integer     :: ngrid
! 
!     w = 0.0_RK
!     ngrid = max(NX,NY,NZ) + 1
! 
! ! Fill the spline function
!     err =  fillspline(w,spline,dspline)
! 
! ! Initialize the bsp_arrays
!     DO i=1,ngrid
!       this%bsp_arr(i) = 0.0_RK
!     END DO
!     DO i=2,this%splineorder
!       this%bsp_arr(i) = spline(i-1)
!     END DO
! 
!     this%bsp_modx = pme_dftmod(NX)
!     this%bsp_mody = pme_dftmod(NY)
!     this%bsp_modZ = pme_dftmod(NZ)
! 
!     end function pme_bspline
! 
! 
! 
! 
!     function pme_dftmod (grid) result(bmod)
! 
!     implicit none
! 
!     real(RK),dimension(this%gridx) :: bmod
!     real(RK),pointer :: barr(:)
!     integer  :: grid
! 
!     real(RK) :: twopi
!     real(RK) :: tinys
!     real(RK) :: sum1,sum2
!     real(RK) :: arg
!     integer  :: i,j
! 
!     twopi = 2._RK*PI
!     tinys = 1.0E-7;
! 
!     barr => this%bsp_arr
!     DO i=1,grid,1
!       sum1 = 0._RK
!       sum2 = 0._RK
!       DO j=1,grid,1
!         arg = twopi * (i-1)*(j-1) / grid
!         sum1 = sum1 + barr(j) * cos(arg)
!         sum2 = sum2 + barr(j) * sin(arg)
!       END DO
!       bmod(i) = sum1*sum1 + sum2*sum2
!     END DO
! 
!     end function pme_dftmod
! 
! 
! 
! !!!!!!!!!!!!!!!!!!
! ! Calculate Spline function and Derivative of the spline function
!     real(RK) function fillspline(w,spline,dspline)
! 
!     implicit none
!       real(RK),intent(in) :: w
!       real(RK),dimension(this%splineorder+5),intent(in out) ::  spline
!       real(RK),dimension(this%splineorder+5),intent(in out) :: dspline
! 
!       real(RK)            :: div
!       real(RK)            :: dj,di,dorder
!       integer             :: i,j,order
! 
! ! Nullify all the arrays
! !       nullify ( this%spline  )
! !       nullify ( this%dspline )
!       order = this%splineorder
! 
!       spline    = 0._RK
!       dspline   = 0._RK
! 
! ! Calculate the single spline function contributions
!       spline(1) = 1._RK - w
!       spline(2) = w
! 
!       DO i=3,order-1,1
!         div = 1._RK / (i-1)
!         spline(i) = div * w * spline(i-1)
!         di = real(i)
! ! Second Contribution
!          DO j=1,i-2,1
!            dj = real(j)
!            spline(i-j) = div * ((w+dj)*spline(i-j-1) + &
! &                                   (di-dj-w)*spline(i-j))
!          END DO
!        spline(1) = spline(1) * div * (1._RK-w)
!        END DO
! 
! ! Differential of the spline function
!     dspline(1) = -spline(1)
!       DO i=2,order,1
!         dspline(i) = spline(i-1) - spline(i)
!       END DO
! 
! ! Generate order spline derivative
!     div = 1._RK / (order - 1._RK)
!     spline(order) = div * w * spline(order-1)
!     dorder = real(order)
!     DO i=1,order-2,1
!       di = real(i)
!       spline(order- i) = div*((w+di)*spline(order-i-1) + &
! &                                  (dorder-di-w)*spline(order-i))
!     END DO
!     spline(1) = spline(1) * div* (1._RK-w)
! 
!     end function fillspline
! 
! 
!     end subroutine TEnsemble_PMESelfTerm





!==============================================================!
!  Subroutine TSimulation_PME_FourierTerm                      !
!==============================================================!

   subroutine TEnsemble_PMEFourierTerm(this)

   implicit none

    include 'fftw3.f'
#if MPI_VER > 0
    include 'mpif.h'
#endif

   ! Declare arguments

   real(RK):: EPotLocal
   real(RK):: Viriallocal
   real(RK):: kap,fac
   real(RK):: boxl,vol
   real(RK):: manhx,manhy,manhz
   real(RK):: den
   real(RK):: mm,bb
   real(RK):: eterm,wterm
   real(RK):: qr,qi,struc
   real(RK):: energ
!    real(RK),pointer:: FX(:), FY(:), FZ(:)
   real(RK),pointer:: FX, FY, FZ
   real(RK):: FXi, FYi, FZi
   real(RK):: FXcum, FYcum, FZcum
   real(RK):: q
   real(RK):: fac2,fac2s
   real(RK):: factor
   real(RK):: facx,facy,facz

   integer :: order
   integer :: NX,NY,NZ
   integer :: ngrid,ngridyz
   integer :: k1,k2,k3,k1bck
   integer :: m1,m2,m3
   integer :: index_loc

   integer :: counter
   integer :: i,j,k
   integer :: inter

!    real(RK):: test

   type(TEnsemble)   :: this
   type(TMolecule),pointer   :: pm

!!!!!!!!!!!! Charge Grid
   real(RK),dimension(this%splineorder+5) ::  transf
   real(RK),dimension(this%splineorder+5) :: dtransf
   real(RK),dimension(this%splineorder+5) ::  splinex
   real(RK),dimension(this%splineorder+5) :: dsplinex
   real(RK),dimension(this%splineorder+5) ::  spliney
   real(RK),dimension(this%splineorder+5) :: dspliney
   real(RK),dimension(this%splineorder+5) ::  splineZ
   real(RK),dimension(this%splineorder+5) :: dsplinez

   real(RK)           :: RX,RY,RZ
   real(RK)           :: RXgit ,RYgit, RZgit
   real(RK)           :: RXgit1,RYgit1,RZgit1
   real(RK)           ::strucx,strucy,strucz
   real(RK)           :: VirLoc

   integer            :: xxo,yyo,zzo
   integer            :: yo,zo

   integer            :: xi,yi,zi
   integer            :: x,y,z

   real(RK)           :: err,err2

   real(RK)           :: qgrid_safe


#if MPI_VER > 0
   integer:: i0,i1
   real(RK),pointer :: qgrid(:,:)
   real(RK),dimension(this%gridx*this%gridx*this%gridx+1):: mult
   real(RK),dimension(this%gridx*this%gridx*this%gridx+1):: mult2
#endif

#if DEBUG_4Part > 0
! DEBUG STEPHAN
    this%Component(1)%Molecule%SiteCharge(1)%RX(1) = 0.50412 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(1) = 0.81073 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(1) = 0.18949 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(1) = 0.50226 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(1) = 0.81178 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(1) = 0.18692 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(2) = 2.09938 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(2) = 1.86064 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(2) = 1.31323 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(2) = 2.09770 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(2) = 1.85859 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(2) = 1.31121 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(3) = 1.31761 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(3) = 2.08757 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(3) = 0.44842 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(3) = 1.31550 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(3) = 2.08538 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(3) = 0.44707 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(4) = 1.48200 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(4) = 0.77144 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(4) = 1.55693 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(4) = 1.48458 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(4) = 0.77169 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(4) = 1.55483 / this%BoxLength
#endif

! Dimensions of charge
   fac2s= (ElementaryCharge / UnitCharge)
   fac2 = (fac2s)**2
   factor = 1._RK / (UnitLength/Angstroem)

   EPotLocal = 0.0_RK
   VirialLocal = 0.0_RK

   kap  = this%Kappa
   fac  = (PI/kap)**2
   boxl = this%boxlength
   vol  = boxl*boxl*boxl

   order = this%splineorder
   NX    = this%gridx
   NY    = this%gridY
   NZ    = this%gridz

   ngrid   = NX*NY*NZ
   ngridyz = NY*NZ

! Calculate the charges on the grid
   call charge_grid(this)

! Backward Transformation of the chargegrid
   call dfftw_execute(this%qgrid_backward)

! Parallelization
# if MPI_VER > 0
   qgrid => this%qgrida
   mult = 0._RK
   mult2 = 0._RK
   j=NProc+1
   i0 = this%NBox0(j)
   i1 = this%NBox2(j)
   err = 0
   DO i=i0,i1,1
# else
   DO i=2,ngrid
#endif
     k1    = int((i-1) / ngridyz)
     k1bck = int((i-1) - k1*ngridyz)
     k2    = int(k1bck / NZ)
     k3    = int(k1bck - k2*NZ)

     m1 = k1
     m2 = k2
     m3 = k3

     if (m1 > NX/2._RK) m1 = m1 - NX
     if (m2 > NY/2._RK) m2 = m2 - NY
     if (m3 > NZ/2._RK) m3 = m3 - NZ

     manhx = m1 / boxl
     manhy = m2 / boxl
     manhz = m3 / boxl

     mm = manhx*manhx + manhy*manhy + manhz*manhz

     bb = this%bsp_modx(k1+1)*this%bsp_mody(k2+1)*this%bsp_modz(k3+1)
     den = PI*vol*mm*bb

     eterm = exp(-fac*mm) / den
     wterm = -2._RK*(fac*mm + 1._RK)

     index_loc = k3 + k2*NZ + k1*NY*NZ + 1

     qr = this%qgrida(1,index_loc)
     qi = this%qgrida(2,index_loc)
     struc = qr*qr + qi*qi

!      test = eterm * struc
!      write (*,*) test
     EPotLocal = EPotLocal + eterm*struc
!      write(*,*) '\n', eterm*struc

     err = err + 1
     VirialLocal = VirialLocal + eterm*struc*(3._RK  + wterm)

# if MPI_VER > 0
     mult(index_loc)      =  eterm * factor
# else
     this%qgrida(1,index_loc) = this%qgrida(1,index_loc) * eterm * factor
     this%qgrida(2,index_loc) = this%qgrida(2,index_loc) * eterm * factor
#endif
   END DO

!     STOP
#if MPI_VER > 0
!    if (RootProc) then
      call MPI_Reduce( EPotLocal, this%UFourier, 1, &
&        MPI_DOUBLE_PRECISION, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
      call MPI_Reduce( VirialLocal, this%EVirial, 1, &
&        MPI_DOUBLE_PRECISION, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
      call MPI_Reduce( mult, mult2, ngrid+1, &
&        MPI_DOUBLE_PRECISION, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
      call MPI_Reduce( err, err2, 1, &
&        MPI_DOUBLE_PRECISION, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
      this%UFourier = 0.5*this%UFourier * fac2
!       this%EVirial  = -0.5*(this%EVirial)*fac2*Third + this%EVirialIntra
      this%EVirial  = -0.5*(this%EVirial)*fac2*Third

! Update of the charge vector on the grid, normally done in the loop above.
      this%qgrida(1,1:ngrid+1) = this%qgrida(1,1:ngrid+1)*mult2(1:ngrid+1)
      this%qgrida(2,1:ngrid+1) = this%qgrida(2,1:ngrid+1)*mult2(1:ngrid+1)
!     end if
#else
   this%UFourier=  0.5*EPotLocal * fac2
!    this%EVirial = -0.5*(Viriallocal)*fac2*Third + this%EVirialIntra
   this%EVirial = -0.5*(Viriallocal)*fac2*Third
#endif

   call dfftw_execute(this%qgrid_forward)

! Last factor cause of multiplication with unitlength/angstroem in multiplication of this%qgrida! (look just some lines above this comment!)
  facx = NX / boxl * fac2s / factor
  facy = NY / boxl * fac2s / factor
  facz = NZ / boxl * fac2s / factor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Calculation of the forces!!!!
!    counter = 1 
   VirLoc = 0._RK
   DO i=1,this%NComponents
     DO j=1,this%Component(i)%Molecule%NCharge
       q =  this%Component(i)%Molecule%SiteCharge(j)%e
       DO k=1,this%Component(i)%NPart

       strucx = 0._RK
       strucy = 0._RK
       strucz = 0._RK

!!!!!!!!!!!!!!!!!!!!!
       q  = this%Component(i)%Molecule%SiteCharge(j)%e
       RX = this%Component(i)%Molecule%SiteCharge(j)%RX(k)
       RY = this%Component(i)%Molecule%SiteCharge(j)%RY(k)
       RZ = this%Component(i)%Molecule%SiteCharge(j)%RZ(k)

       RX     = RX - anint(RX)
       RY     = RY - anint(RY)
       RZ     = RZ - anint(RZ)

       RXgit1 = this%gridx*(RX+0.5)
       RYgit1 = this%gridy*(RY+0.5)
       RZgit1 = this%gridz*(RZ+0.5)

       RXgit = (RXgit1 - floor( RXgit1 ))
       RYgit = (RYgit1 - floor( RYgit1 ))
       RZgit = (RZgit1 - floor( RZgit1 ))

! x-coordinate
       err=fillspline(RXgit,transf,dtransf)
         splinex  = transf
         dsplinex = dtransf

! y-coordinate
       err=fillspline(RYgit,transf,dtransf)
         spliney  = transf
         dspliney = dtransf

! z-coordinate
       err=fillspline(RZgit,transf,dtransf)
         splinez  = transf
         dsplinez = dtransf

       xxo = floor (RXgit1 - this%splineorder)
       yyo = floor (RYgit1 - this%splineorder)
       zzo = floor (RZgit1 - this%splineorder)

!        if (xxo > NX-1) xxo = NX-1
!        if (yyo > NY-1) yyo = NY-1
!        if (zzo > NZ-1) zzo = Nz-1

       DO xi=1,order,1
         xxo = xxo + 1
         x = xxo
         if (x < 0) x = x + NX
         yo   = yyo
         DO yi=1,order,1
           yo = yo + 1
           y   = yo
           if (y < 0) y = y + NY
           zo  = zzo
           DO zi=1,order,1
             zo = zo + 1
             z  = zo
             if (z < 0) z = z + NZ
             index_loc = z + y*NZ + x*NY*NZ + 1
             qgrid_safe = this%qgrida(1,index_loc)
               strucx = strucx - &
&                    dsplinex(xi)*spliney(yi)*splinez(zi)*qgrid_safe
               strucy = strucy - &
&                    splinex(xi)*dspliney(yi)*splinez(zi)*qgrid_safe
               strucz = strucz - &
&                    splinex(xi)*spliney(yi)*dsplinez(zi)*qgrid_safe
           END DO
         END DO
       END DO
!!!!!!!!!1



         FX =>  this%Component(i)%Molecule%SiteCharge(j)%FX(k)
         FY =>  this%Component(i)%Molecule%SiteCharge(j)%FY(k)
         FZ =>  this%Component(i)%Molecule%SiteCharge(j)%FZ(k)

         FXi = strucx*facx*q
         FYi = strucy*facy*q
         FZi = strucz*facz*q

         FX    = FX + FXi
         FY    = FY + FYi
         FZ    = FZ + FZi

         FXcum = FXcum + FXi
         FYcum = FYcum + FYi
         FZcum = FZcum + FZi

         VirLoc = VirLoc + FXi*(this%Component(i)%Molecule%SiteCharge(j)%RX(k) - &
&                    this%Component(i)%Molecule%SiteCharge(j)%PX(k))*this%BoxLength
         VirLoc = VirLoc + FYi*(this%Component(i)%Molecule%SiteCharge(j)%RY(k) - &
&                    this%Component(i)%Molecule%SiteCharge(j)%PY(k))*this%BoxLength
         VirLoc = VirLoc + FZi*(this%Component(i)%Molecule%SiteCharge(j)%RZ(k) - &
&                    this%Component(i)%Molecule%SiteCharge(j)%PZ(k))*this%BoxLength
!          counter = counter + 3

       END DO
     END DO
   END DO

   DO i=1,this%NComponents
     DO j=1,this%Component(i)%Molecule%NCharge
       pm=>this%Component(i)%Molecule
       DO k=1,this%Component(i)%NPart
         pm%SiteCharge(j)%FX(k) = pm%SiteCharge(j)%FX(k) - FXcum/this%NPart
         pm%SiteCharge(j)%FY(k) = pm%SiteCharge(j)%FY(k) - FYcum/this%NPart
         pm%SiteCharge(j)%FZ(k) = pm%SiteCharge(j)%FZ(k) - FZcum/this%NPart
       END DO
     END DO
   END DO

   call TEnsemble_VirialIntra(this)
   this%EVirial = this%EVirial + VirLoc*Third

contains

!!!!!!!!!!!!!!!!!!
! Calculate Spline function and Derivative of the spline function
    real(RK) function fillspline(w,spline,dspline)

    implicit none
      real(RK),intent(in) ::  w
      real(RK),dimension(this%splineorder+5),intent(in out) ::  spline
      real(RK),dimension(this%splineorder+5),intent(in out) :: dspline


      real(RK)            :: div
      real(RK)            :: dj,di,dorder
      integer             :: i,j,order

      order = this%splineorder

! Calculate the single spline function contributions
      spline    = 0._RK
      dspline   = 0._RK


      spline(1) = 1._RK - w
      spline(2) = w

      DO i=3,order-1,1
        div = 1._RK / (i-1)
        spline(i) = div * w * spline(i-1)
        di = real(i)
! Second Contribution
         DO j=1,i-2,1
           dj = real(j)
           spline(i-j) = div * ((w+dj)*spline(i-j-1) + &
&                                   (di-dj-w)*spline(i-j))
         END DO
       spline(1) = spline(1) * div * (1._RK-w)
       END DO

! Differential of the spline function
    dspline(1) = -spline(1)
      DO i=2,order,1
        dspline(i) = spline(i-1) - spline(i)
      END DO

! Generate order spline derivative
    div = 1._RK / (order - 1._RK)
    spline(order) = div * w * spline(order-1)
    dorder = real(order)
    DO i=1,order-2,1
      di = real(i)
      spline(order- i) = div*((w+di)*spline(order-i-1) + &
&                                  (dorder-di-w)*spline(order-i))
    END DO
    spline(1) = spline(1) * div* (1._RK-w)

    end function fillspline



   end subroutine TEnsemble_PMEFourierTerm



   subroutine charge_grid(this)

   implicit none
   type(TEnsemble)    :: this

   real(RK)           :: boxl
   real(RK),dimension(this%splineorder+5) ::  transf
   real(RK),dimension(this%splineorder+5) :: dtransf
   real(RK),dimension(this%splineorder+5) ::  splinex
   real(RK),dimension(this%splineorder+5) :: dsplinex
   real(RK),dimension(this%splineorder+5) ::  spliney
   real(RK),dimension(this%splineorder+5) :: dspliney
   real(RK),dimension(this%splineorder+5) ::  splineZ
   real(RK),dimension(this%splineorder+5) :: dsplinez

   real(RK)           :: RX,RY,RZ
   real(RK)           :: RXgit ,RYgit, RZgit
   real(RK)           :: RXgit1,RYgit1,RZgit1
   real(RK)           :: q

   integer            :: xxo,yyo,zzo
   integer            :: yo,zo

   integer            :: i,j,k
   integer            :: xi,yi,zi
   integer            :: x,y,z
   integer            :: order
   integer            :: NX,NY,NZ
   integer            :: index_loc

   real(RK)           :: err
   real(RK)           :: fac


!    nullify(this%qgrida)
!    write(*,*) 'nicht ausgenullt!'
   this%qgrida = 0._RK

   NX = this%gridx
   NY = this%gridy
   NZ = this%gridz

   boxl  = this%BoxLength
   order = this%splineorder

! Debug Stephan
!    this%Component(1)%molecule%SiteCharge(1)%RX(1) = -0.3604397088940393
!    this%Component(1)%molecule%SiteCharge(1)%RY(1) = -0.2647705271429488
!    this%Component(1)%molecule%SiteCharge(1)%RZ(1) = 0.2643293164661566
!    this%Component(1)%molecule%SiteCharge(2)%RX(1) = -0.1543060806967627
!    this%Component(1)%molecule%SiteCharge(2)%RY(1) = -0.2437836059502038
!    this%Component(1)%molecule%SiteCharge(2)%RZ(1) = 0.2319150387444959
!    this%Component(1)%molecule%SiteCharge(3)%RX(1) = -0.1114351099351288
!    this%Component(1)%molecule%SiteCharge(3)%RY(1) = -0.1269657257582112
!    this%Component(1)%molecule%SiteCharge(3)%RZ(1) = 0.2941404511947440
! 
!    this%Component(1)%molecule%SiteCharge(1)%RX(2) = -0.2553580327047170 
!    this%Component(1)%molecule%SiteCharge(1)%RY(2) = 0.1997507804063506
!    this%Component(1)%molecule%SiteCharge(1)%RZ(2) = -0.1536295576590148
!    this%Component(1)%molecule%SiteCharge(2)%RX(2) = 0.2790804600935731
!    this%Component(1)%molecule%SiteCharge(2)%RY(2) = 0.3321433974891049
!    this%Component(1)%molecule%SiteCharge(2)%RZ(2) = -0.3145390914850995
!    this%Component(1)%molecule%SiteCharge(3)%RX(2) = -0.3085680403258462
!    this%Component(1)%molecule%SiteCharge(3)%RY(2) = 0.2582111950813060
!    this%Component(1)%molecule%SiteCharge(3)%RZ(2) = -0.4286508795261004
! 
!    this%Component(1)%molecule%SiteCharge(1)%RX(3) =  0.1660540503174808
!    this%Component(1)%molecule%SiteCharge(1)%RY(3) = -0.3100534496043796
!    this%Component(1)%molecule%SiteCharge(1)%RZ(3) = -0.2347976151662044
!    this%Component(1)%molecule%SiteCharge(2)%RX(3) =0.3524229108967265
!    this%Component(1)%molecule%SiteCharge(2)%RY(3) = -0.2156196577483070
!    this%Component(1)%molecule%SiteCharge(2)%RZ(3) = -0.2529313739823603
!    this%Component(1)%molecule%SiteCharge(3)%RX(3) = 0.3638105584647310 
!    this%Component(1)%molecule%SiteCharge(3)%RY(3) =  -0.1205093428544866
!    this%Component(1)%molecule%SiteCharge(3)%RZ(3) = -0.1520264922000034
! 
!    this%Component(1)%molecule%SiteCharge(1)%RX(4) = 0.2965553442990528
!    this%Component(1)%molecule%SiteCharge(1)%RY(4) = 0.3442031559880894
!    this%Component(1)%molecule%SiteCharge(1)%RZ(4) = 0.2098397978823303
!    this%Component(1)%molecule%SiteCharge(2)%RX(4) = 0.1884454921646807 
!    this%Component(1)%molecule%SiteCharge(2)%RY(4) = 0.1653657616616851  
!    this%Component(1)%molecule%SiteCharge(2)%RZ(4) = 0.2274146898412163
!    this%Component(1)%molecule%SiteCharge(3)%RX(4) = 0.2750154390580618
!    this%Component(1)%molecule%SiteCharge(3)%RY(4) = 0.0654462503908305
!    this%Component(1)%molecule%SiteCharge(3)%RZ(4) = 0.2707268712796423
! 
! 
! 
! ! 5 Digits like in Spidey
!    this%Component(1)%molecule%SiteCharge(1)%RX(1) = -0.36043
!    this%Component(1)%molecule%SiteCharge(1)%RY(1) = -0.26477
!    this%Component(1)%molecule%SiteCharge(1)%RZ(1) =  0.26432
!    this%Component(1)%molecule%SiteCharge(2)%RX(1) = -0.15430
!    this%Component(1)%molecule%SiteCharge(2)%RY(1) = -0.24378
!    this%Component(1)%molecule%SiteCharge(2)%RZ(1) =  0.23191
!    this%Component(1)%molecule%SiteCharge(3)%RX(1) = -0.11143
!    this%Component(1)%molecule%SiteCharge(3)%RY(1) = -0.12696
!    this%Component(1)%molecule%SiteCharge(3)%RZ(1) =  0.29414
! 
!    this%Component(1)%molecule%SiteCharge(1)%RX(2) = -0.25535
!    this%Component(1)%molecule%SiteCharge(1)%RY(2) =  0.19975
!    this%Component(1)%molecule%SiteCharge(1)%RZ(2) = -0.15362
!    this%Component(1)%molecule%SiteCharge(2)%RX(2) =  0.27908
!    this%Component(1)%molecule%SiteCharge(2)%RY(2) =  0.33214
!    this%Component(1)%molecule%SiteCharge(2)%RZ(2) = -0.31453
!    this%Component(1)%molecule%SiteCharge(3)%RX(2) = -0.30856
!    this%Component(1)%molecule%SiteCharge(3)%RY(2) =  0.25821
!    this%Component(1)%molecule%SiteCharge(3)%RZ(2) = -0.42865
! 
!    this%Component(1)%molecule%SiteCharge(1)%RX(3) =  0.16605
!    this%Component(1)%molecule%SiteCharge(1)%RY(3) = -0.31005
!    this%Component(1)%molecule%SiteCharge(1)%RZ(3) = -0.23479
!    this%Component(1)%molecule%SiteCharge(2)%RX(3) =  0.35242
!    this%Component(1)%molecule%SiteCharge(2)%RY(3) = -0.21561
!    this%Component(1)%molecule%SiteCharge(2)%RZ(3) = -0.25293
!    this%Component(1)%molecule%SiteCharge(3)%RX(3) =  0.36381
!    this%Component(1)%molecule%SiteCharge(3)%RY(3) = -0.12050
!    this%Component(1)%molecule%SiteCharge(3)%RZ(3) = -0.15202
! 
!    this%Component(1)%molecule%SiteCharge(1)%RX(4) = 0.296555
!    this%Component(1)%molecule%SiteCharge(1)%RY(4) = 0.344203
!    this%Component(1)%molecule%SiteCharge(1)%RZ(4) = 0.209839
!    this%Component(1)%molecule%SiteCharge(2)%RX(4) = 0.188445
!    this%Component(1)%molecule%SiteCharge(2)%RY(4) = 0.165365
!    this%Component(1)%molecule%SiteCharge(2)%RZ(4) = 0.227414
!    this%Component(1)%molecule%SiteCharge(3)%RX(4) = 0.275015
!    this%Component(1)%molecule%SiteCharge(3)%RY(4) = 0.065446
!    this%Component(1)%molecule%SiteCharge(3)%RZ(4) = 0.270726


! Dimensions of charge
   fac = UnitCharge / ElementaryCharge

   DO i=1,this%NComponents
     DO j=1,this%Component(i)%Molecule%NCharge
      DO k=1,this%Component(i)%NPart

       q  = this%Component(i)%Molecule%SiteCharge(j)%e
       RX = this%Component(i)%Molecule%SiteCharge(j)%RX(k)
       RY = this%Component(i)%Molecule%SiteCharge(j)%RY(k)
       RZ = this%Component(i)%Molecule%SiteCharge(j)%RZ(k)

       RX     = RX - anint(RX)
       RY     = RY - anint(RY)
       RZ     = RZ - anint(RZ)

       RXgit1 = this%gridx*(RX+0.5)
       RYgit1 = this%gridy*(RY+0.5)
       RZgit1 = this%gridz*(RZ+0.5)

       RXgit = (RXgit1 - floor( RXgit1 ))
       RYgit = (RYgit1 - floor( RYgit1 ))
       RZgit = (RZgit1 - floor( RZgit1 ))

! x-coordinate
       err=fillspline(RXgit,transf,dtransf)
         splinex  = transf
         dsplinex = dtransf

! y-coordinate
       err=fillspline(RYgit,transf,dtransf)
         spliney  = transf
         dspliney = dtransf

! z-coordinate
       err=fillspline(RZgit,transf,dtransf)
         splinez  = transf
         dsplinez = dtransf

       xxo = floor (RXgit1 - this%splineorder)
       yyo = floor (RYgit1 - this%splineorder)
       zzo = floor (RZgit1 - this%splineorder)

       if (xxo > NX-1) xxo = NX-1
       if (yyo > NY-1) yyo = NY-1
       if (zzo > NZ-1) zzo = Nz-1

       DO xi=1,order,1
         xxo = xxo + 1
         x = xxo
         if (x < 0) x = x + NX
         yo   = yyo
         DO yi=1,order,1
           yo = yo + 1
           y   = yo
           if (y < 0) y = y + NY
           zo  = zzo
           DO zi=1,order,1
             zo = zo + 1
             z  = zo
             if (z < 0) z = z + NZ
             index_loc = z + y*NZ + x*NY*NZ + 1
!              if (index_loc .eq. 26) then
!                write(*,*) 'halt'
!              end if
!              write(*,*) i,j,k,RX,RY,RZ, RXgit1, RYgit1, RZgit1
             this%qgrida(1,index_loc) = this%qgrida(1,index_loc) + splinex(xi)*spliney(yi)*splinez(zi)*q*fac
           END DO
         END DO
       END DO

      END DO ! Particles
     END DO  ! Charges
   END DO ! Components


contains

!!!!!!!!!!!!!!!!!!
! Calculate Spline function and Derivative of the spline function
    real(RK) function fillspline(w,spline,dspline)

    implicit none
      real(RK),intent(in) ::  w
      real(RK),dimension(this%splineorder+5),intent(in out) ::  spline
      real(RK),dimension(this%splineorder+5),intent(in out) :: dspline


      real(RK)            :: div
      real(RK)            :: dj,di,dorder
      integer             :: i,j,order

      order = this%splineorder

! Calculate the single spline function contributions
      spline    = 0._RK
      dspline   = 0._RK


      spline(1) = 1._RK - w
      spline(2) = w

      DO i=3,order-1,1
        div = 1._RK / (i-1)
        spline(i) = div * w * spline(i-1)
        di = real(i)
! Second Contribution
         DO j=1,i-2,1
           dj = real(j)
           spline(i-j) = div * ((w+dj)*spline(i-j-1) + &
&                                   (di-dj-w)*spline(i-j))
         END DO
       spline(1) = spline(1) * div * (1._RK-w)
       END DO

! Differential of the spline function
    dspline(1) = -spline(1)
      DO i=2,order,1
        dspline(i) = spline(i-1) - spline(i)
      END DO

! Generate order spline derivative
    div = 1._RK / (order - 1._RK)
    spline(order) = div * w * spline(order-1)
    dorder = real(order)
    DO i=1,order-2,1
      di = real(i)
      spline(order- i) = div*((w+di)*spline(order-i-1) + &
&                                  (dorder-di-w)*spline(order-i))
    END DO
    spline(1) = spline(1) * div* (1._RK-w)

    end function fillspline

   end subroutine charge_grid


















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SPME for MonteCarlo Simulations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine TEnsemble_PME_Setup ( this )

    implicit none

    include 'fftw3.f'
    ! Declare arguments
    type(TEnsemble)         :: this

    real(RK)::bspx(this%splineorder)
    real(RK)::bspy(this%splineorder)
    real(RK)::bspz(this%splineorder)
    real(RK)::err

    integer :: NX,NY,NZ
    integer :: i
    integer :: k1,k1bck,k2,k3,ngrid,ngridyz
    integer :: stat
    integer :: manhx,manhy,manhz

    NX = this%gridx
    NY = this%gridy
    NZ = this%gridz

    if (SimulationType .eq. MolecularDynamics) then
    call dfftw_plan_dft_3d(this%qgrid_forward,NX,NY,NZ,this%qgrida,this%qgrida,FFTW_FORWARD,FFTW_PATIENT)
    call dfftw_plan_dft_3d(this%qgrid_backward,NX,NY,NZ,this%qgrida,this%qgrida,FFTW_BACKWARD,FFTW_PATIENT)
    else if (SimulationType .eq. MonteCarlo) then
    call dfftw_plan_dft_3d(this%qgrid_backward,NX,NY,NZ,this%qgrida,this%qgridb,FFTW_BACKWARD,FFTW_PATIENT)
    end if

    err = pme_bspline(NX,NY,NZ)

! Self Term
    call PMESelfTermMC (this )


! Allocation of bbtot
    allocate(this%bbtot(NX*NY*NZ+1),STAT=stat)
    if(stat >0) write(*,*) 'Allocation Error bbtot'
    allocate(this%mm2(NX*NY*NZ+1),STAT=stat)
    if(stat >0) write(*,*) 'Allocation Error mm2'
    allocate(this%EPotPME(NX*NY*NZ+1),STAT=stat)
    if(stat >0) write(*,*) 'Allocation Error EPotPME'
    allocate(this%VirialPME(NX*NY*NZ+1),STAT=stat)
    if(stat >0) write(*,*) 'Allocation Error VirialPME'

! Preparation of bb
   ngrid   = NX*NY*NZ
   ngridyz = NY*NZ

   DO i=2,ngrid
     k1    = int((i-1) / ngridyz)
     k1bck = int((i-1) - k1*ngridyz)
     k2    = int(k1bck / NZ)
     k3    = int(k1bck - k2*NZ)

! Further specifications
     manhx = k1
     manhy = k2
     manhz = k3

     if (manhx > NX/2._RK) manhx = manhx - NX ! manhx = manhx - ANINT(manhx/NX)*NX
     if (manhy > NY/2._RK) manhy = manhy - NY! manhy = manhy - ANINT(manhx/NY)*NY
     if (manhz > NZ/2._RK) manhz = manhz - NZ! manhz = manhz - ANINT(manhx/NZ)*NZ

     this%mm2(i) = manhx*manhx + manhy*manhy + manhz*manhz

     this%bbtot(i)= this%bsp_modx(k1+1)*this%bsp_mody(k2+1)*this%bsp_modz(k3+1)*this%mm2(i)

   END DO

#if MPI_VER > 0
    this%NBox1 = 1 +(ngrid - 1) / NProcs
    DO i=1,NProcs
      this%NBox0(i) = 1 + this%NBox1(i) * (i-1)
    END DO
    this%NBox2 = min(this%NBox0 + this%NBox1 - 1, ngrid)
    this%NBox0(1) = 2
    this%NBox1 = this%NBox2 - this%NBox0
#endif


contains

!!!!!!!!!!!!!!!!!!!!
    real(RK) function pme_bspline(NX,NY,NZ)

    implicit none
    integer,intent(in) :: NX,NY,NZ

    real(RK)    :: w
    real(RK)    :: err
    real(RK),dimension(this%splineorder+5):: spline,dspline
    integer     :: ngrid

    w = 0.0_RK
    ngrid = max(NX,NY,NZ) + 1

! Fill the spline function
    err =  fillspline(w,spline,dspline)

! Initialize the bsp_arrays
    DO i=1,ngrid
      this%bsp_arr(i) = 0.0_RK
    END DO
    DO i=2,this%splineorder
      this%bsp_arr(i) = spline(i-1)
    END DO

    this%bsp_modx = pme_dftmod(NX)
    this%bsp_mody = pme_dftmod(NY)
    this%bsp_modZ = pme_dftmod(NZ)

    end function pme_bspline




    function pme_dftmod (grid) result(bmod)

    implicit none

    real(RK),dimension(this%gridx) :: bmod
    real(RK),pointer :: barr(:)
    integer  :: grid

    real(RK) :: twopi
    real(RK) :: tinys
    real(RK) :: sum1,sum2
    real(RK) :: arg
    integer  :: i,j

    twopi = 2._RK*PI
    tinys = 1.0E-7;

    barr => this%bsp_arr
    DO i=1,grid,1
      sum1 = 0._RK
      sum2 = 0._RK
      DO j=1,grid,1
        arg = twopi * (i-1)*(j-1) / grid
        sum1 = sum1 + barr(j) * cos(arg)
        sum2 = sum2 + barr(j) * sin(arg)
      END DO
      bmod(i) = sum1*sum1 + sum2*sum2
    END DO

    end function pme_dftmod



!!!!!!!!!!!!!!!!!!
! Calculate Spline function and Derivative of the spline function
    real(RK) function fillspline(w,spline,dspline)

    implicit none
      real(RK),intent(in) :: w
      real(RK),dimension(this%splineorder+5),intent(in out) ::  spline
      real(RK),dimension(this%splineorder+5),intent(in out) :: dspline

      real(RK)            :: div
      real(RK)            :: dj,di,dorder
      integer             :: i,j,order

! Nullify all the arrays
!       nullify ( this%spline  )
!       nullify ( this%dspline )
      order = this%splineorder

      spline    = 0._RK
      dspline   = 0._RK

! Calculate the single spline function contributions
      spline(1) = 1._RK - w
      spline(2) = w

      DO i=3,order-1,1
        div = 1._RK / (i-1)
        spline(i) = div * w * spline(i-1)
        di = real(i)
! Second Contribution
         DO j=1,i-2,1
           dj = real(j)
           spline(i-j) = div * ((w+dj)*spline(i-j-1) + &
&                                   (di-dj-w)*spline(i-j))
         END DO
       spline(1) = spline(1) * div * (1._RK-w)
       END DO

! Differential of the spline function
    dspline(1) = -spline(1)
      DO i=2,order,1
        dspline(i) = spline(i-1) - spline(i)
      END DO

! Generate order spline derivative
    div = 1._RK / (order - 1._RK)
    spline(order) = div * w * spline(order-1)
    dorder = real(order)
    DO i=1,order-2,1
      di = real(i)
      spline(order- i) = div*((w+di)*spline(order-i-1) + &
&                                  (dorder-di-w)*spline(order-i))
    END DO
    spline(1) = spline(1) * div* (1._RK-w)

    end function fillspline

  end subroutine TEnsemble_PME_Setup





  subroutine TEnsemble_PMESelfTerm_MC( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif
    ! Declare arguments
    type(TEnsemble)         :: this

    ! Declare local variables
    integer :: i,Si,Sj
    real(RK):: USelbstTermKomp, VirialLocal
    real(RK):: Faktor, NSQ
    real(RK):: approx
    real(RK):: twopi
    real(RK):: UIntraTermKomp
    real(RK):: dr, drxij,dryij,drzij
    real(RK):: RXi,RYi,RZi,RXj,RYj,RZj
    real(RK):: q1,q2

! Selbstterm
    this%USelbstTerm = 0.0
    DO i=1,this%NComponents,1
       USelbstTermKomp = 0.0
       DO Si=1,this%Component(i)%Molecule%NCharge,1
         USelbstTermKomp = USelbstTermKomp + &
&            this%Component(i)%Molecule%SiteCharge(Si)%e**2
       END DO
       this%USelbstTerm = this%USelbstTerm + &
&            this%Component(i)%NPart * USelbstTermKomp
    END DO
    this%USelbstTerm = -this%USelbstTerm * this%Kappa / sqrt(Pi)


! Intramolecular
    this%UIntra = 0._RK
!     this%EVirialIntra = 0._RK
    DO i=1,this%NComponents,1
      UIntraTermKomp = 0.0
!       VirialLocal    = 0._RK
      DO Si = 1,this%component(i)%Molecule%NCharge-1
        q1 = this%Component(i)%Molecule%SiteCharge(Si)%e
        DO Sj = Si+1,this%component(i)%Molecule%NCharge
          q2 = this%Component(i)%Molecule%SiteCharge(Sj)%e
          RXi  = this%Component(i)%Molecule%SiteCharge(Si)%r(1)
          RYi  = this%Component(i)%Molecule%SiteCharge(Si)%r(2)
          RZi  = this%Component(i)%Molecule%SiteCharge(Si)%r(3)
          RXj  = this%Component(i)%Molecule%SiteCharge(Sj)%r(1)
          RYj  = this%Component(i)%Molecule%SiteCharge(Sj)%r(2)
          RZj  = this%Component(i)%Molecule%SiteCharge(Sj)%r(3)

          drxij = (RXi-RXj)
          dryij = (RYi-RYj)
          drzij = (RZi-RZj)

          dr = sqrt(drxij*drxij + dryij*dryij + drzij*drzij)

          call erfc_approx (this%Kappa*dr, approx)

          UIntraTermKomp = UIntraTermKomp - this%Component(i)%Molecule%SiteCharge(Si)%e* &
&                   this%Component(i)%Molecule%SiteCharge(Sj)%e / dr * (1-approx)

!           VirialLocal = VirialLocal + q1 * q2 / dr * (1-approx) -2._RK*this%Kappa&
! &                   / sqrt(pi) * exp(-(this%Kappa*dr)**2) * q1*q2
        END DO
      END DO
      this%UIntra = this%UIntra + this%component(i)%NPart * UIntraTermKomp
!       this%EVirialIntra = this%EVirialIntra + VirialLocal *this%component(i)%NPart
    END DO

!     this%EVirialIntra = this%EVirialIntra * Third


    end subroutine TEnsemble_PMESelfTerm_MC





!==============================================================!
!  Subroutine TSimulation_PME_FourierTerm MonteCarlo           !
!==============================================================!

   subroutine TEnsemble_PMEFourierTerm_MC(this)

   implicit none

    include 'fftw3.f'
#if MPI_VER > 0
    include 'mpif.h'
#endif

   ! Declare arguments

! Factors
   real(RK):: fac2, fac2s
   real(RK):: factor
   real(RK):: fac
   real(RK):: mult
! Energy
   real(RK):: EPotLocal
   real(RK):: Viriallocal
   real(RK):: boxl,boxl2,vol
!    real(RK):: manhx,manhy,manhz
!    real(RK)::mm2,den2
!    real(RK):: bb
   real(RK):: eterm,wterm
   real(RK):: struc
   real(RK):: energ
!Pointers
   real(RK),pointer::lad(:,:)
!    real(RK),pointer:: den(:)
!    real(RK),pointer:: mm(:)
!    real(RK),pointer:: qr,qi
   real(RK):: qr,qi
!Positioning
   integer :: order
   integer :: NX,NY,NZ
   integer :: ngrid,ngridyz
   integer :: k1,k2,k3,k1bck
   integer :: m1,m2,m3
   integer :: index_loc

   integer :: i,j,k

   type(TEnsemble)   :: this



#if MPI_VER > 0
   integer:: i0,i1
#endif

#if DEBUG_4Part > 0
! DEBUG STEPHAN
    this%Component(1)%Molecule%SiteCharge(1)%RX(1) = 0.50412 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(1) = 0.81073 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(1) = 0.18949 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(1) = 0.50226 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(1) = 0.81178 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(1) = 0.18692 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(2) = 2.09938 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(2) = 1.86064 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(2) = 1.31323 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(2) = 2.09770 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(2) = 1.85859 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(2) = 1.31121 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(3) = 1.31761 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(3) = 2.08757 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(3) = 0.44842 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(3) = 1.31550 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(3) = 2.08538 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(3) = 0.44707 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RX(4) = 1.48200 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RY(4) = 0.77144 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(1)%RZ(4) = 1.55693 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RX(4) = 1.48458 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RY(4) = 0.77169 / this%BoxLength
    this%Component(1)%Molecule%SiteCharge(2)%RZ(4) = 1.55483 / this%BoxLength
#endif

! Factors
   boxl = this%boxlength
   boxl2= boxl**2
   vol  = boxl**3

! Spline factors
   order = this%splineorder
   NX    = this%gridx
   NY    = this%gridY
   NZ    = this%gridz

   ngrid   = NX*NY*NZ
   ngridyz = NY*NZ

! Dimensions of charge
   fac  = (PI/this%Kappa)**2 / boxl2
   fac2s= (ElementaryCharge / UnitCharge)
   fac2 = (fac2s)**2 / (PI*boxl)
   factor = 1._RK / (UnitLength/Angstroem)

! Backward Transformation of the chargegrid
   call dfftw_execute(this%qgrid_backward)

   lad => this%qgridb
!    den => this%bbtot
!    mm  => this%mm2

   EPotLocal = 0._RK
   VirialLocal = 0._RK

! Summation over all the Energies
# if MPI_VER > 0
   j=NProc+1
   i0 = this%NBox0(j)
   i1 = this%NBox2(j)
   DO i=i0,i1,1
# else
   DO i=2,ngrid
#endif
! Positioning
     k1    = int((i-1) / ngridyz)
     k1bck = int((i-1) - k1*ngridyz)
     k2    = int(k1bck / NZ)
     k3    = int(k1bck - k2*NZ)

     index_loc = k3 + k2*NZ + k1*ngridyz + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      manhx = k1
!      manhy = k2
!      manhz = k3
! 
!      if (manhx > NX/2._RK) manhx = manhx-NX! manhx = manhx - ANINT(manhx/NX)*NX
!      if (manhy > NY/2._RK) manhy = manhy-NY! manhy = manhy - ANINT(manhx/NY)*NY
!      if (manhz > NZ/2._RK) manhz = manhz-NZ! manhz = manhz - ANINT(manhx/NZ)*NZ
! 
!      mm2 = manhx*manhx + manhy*manhy + manhz*manhz
! 
!      bb = this%bsp_modx(k1+1)*this%bsp_mody(k2+1)*this%bsp_modz(k3+1)
! !      den = this%bbtot(i)
!      den2 = mm2*bb
! 
!      mult = mm2*fac
!      eterm = exp(-mult) / den2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Charge contribution
     qr = lad(1,index_loc)
     qi = lad(2,index_loc)
     struc = qr*qr + qi*qi

! Distance contribution
     mult  = fac*this%mm2(i)
     eterm = exp(-mult) / this%bbtot(i)
!      wterm = -2._RK*(fac*mm(i) + 1._RK)
!      wterm = -2._RK*(fac*mm(i) + 1._RK)

! Energy
     energ = eterm*struc

!      this%EPotPME(i) = energ
     EPotLocal = EPotLocal + energ
!      this%VirialPME(i) = energ*(3._RK  + wterm)
!      this%VirialPME(i) = energ*(1._RK  - mult)
     VirialLocal = VirialLocal + energ*(1._RK - mult)
!      qr = qr * eterm * factor 
!      qi = qi * eterm * factor 
!      write(*,*) index_loc
   END DO

!    this%UFourier=  0.5*sum(this%EPotPME) * fac2
!    this%EVirial = -0.5*(sum(this%VirialPME))*fac2*Third + this%EVirialIntra
   this%UFourier =  0.5*EPotLocal*fac2
   this%EVirial  = -0.5*VirialLocal*fac2*Third + this%EVirialIntra
! STOP
! ForwardTransformation of the chargegrid

   end subroutine TEnsemble_PMEFourierTerm_MC





   subroutine charge_grid_MCall(this)

   implicit none
   type(TEnsemble)    :: this

   real(RK)           :: boxl
   real(RK),dimension(this%splineorder+5) ::  transf
   real(RK),dimension(this%splineorder+5) ::  splinex
   real(RK),dimension(this%splineorder+5) ::  spliney
   real(RK),dimension(this%splineorder+5) ::  splineZ

   real(RK)           :: RX,RY,RZ
   real(RK)           :: RXgit ,RYgit, RZgit
   real(RK)           :: RXgit1,RYgit1,RZgit1
   real(RK)           :: q

   integer            :: xxo,yyo,zzo
   integer            :: yo,zo

   integer            :: i,j,k
   integer            :: xi,yi,zi
   integer            :: x,y,z
   integer            :: order
   integer            :: NX,NY,NZ
   integer            :: index_loc

   real(RK)           :: err
   real(RK)           :: fac


!    nullify(this%qgrida)
!    write(*,*) 'nicht ausgenullt!'
   this%qgrida = 0._RK

   NX = this%gridx
   NY = this%gridy
   NZ = this%gridz

   boxl  = this%BoxLength
   order = this%splineorder

! 
! ! 5 Digits like in Spidey
!    this%Component(1)%molecule%SiteCharge(1)%RX(1) = -0.36043
!    this%Component(1)%molecule%SiteCharge(1)%RY(1) = -0.26477
!    this%Component(1)%molecule%SiteCharge(1)%RZ(1) =  0.26432
!    this%Component(1)%molecule%SiteCharge(2)%RX(1) = -0.15430
!    this%Component(1)%molecule%SiteCharge(2)%RY(1) = -0.24378
!    this%Component(1)%molecule%SiteCharge(2)%RZ(1) =  0.23191
!    this%Component(1)%molecule%SiteCharge(3)%RX(1) = -0.11143
!    this%Component(1)%molecule%SiteCharge(3)%RY(1) = -0.12696
!    this%Component(1)%molecule%SiteCharge(3)%RZ(1) =  0.29414
! 
!    this%Component(1)%molecule%SiteCharge(1)%RX(2) = -0.25535
!    this%Component(1)%molecule%SiteCharge(1)%RY(2) =  0.19975
!    this%Component(1)%molecule%SiteCharge(1)%RZ(2) = -0.15362
!    this%Component(1)%molecule%SiteCharge(2)%RX(2) =  0.27908
!    this%Component(1)%molecule%SiteCharge(2)%RY(2) =  0.33214
!    this%Component(1)%molecule%SiteCharge(2)%RZ(2) = -0.31453
!    this%Component(1)%molecule%SiteCharge(3)%RX(2) = -0.30856
!    this%Component(1)%molecule%SiteCharge(3)%RY(2) =  0.25821
!    this%Component(1)%molecule%SiteCharge(3)%RZ(2) = -0.42865
! 
!    this%Component(1)%molecule%SiteCharge(1)%RX(3) =  0.16605
!    this%Component(1)%molecule%SiteCharge(1)%RY(3) = -0.31005
!    this%Component(1)%molecule%SiteCharge(1)%RZ(3) = -0.23479
!    this%Component(1)%molecule%SiteCharge(2)%RX(3) =  0.35242
!    this%Component(1)%molecule%SiteCharge(2)%RY(3) = -0.21561
!    this%Component(1)%molecule%SiteCharge(2)%RZ(3) = -0.25293
!    this%Component(1)%molecule%SiteCharge(3)%RX(3) =  0.36381
!    this%Component(1)%molecule%SiteCharge(3)%RY(3) = -0.12050
!    this%Component(1)%molecule%SiteCharge(3)%RZ(3) = -0.15202
! 
!    this%Component(1)%molecule%SiteCharge(1)%RX(4) = 0.296555
!    this%Component(1)%molecule%SiteCharge(1)%RY(4) = 0.344203
!    this%Component(1)%molecule%SiteCharge(1)%RZ(4) = 0.209839
!    this%Component(1)%molecule%SiteCharge(2)%RX(4) = 0.188445
!    this%Component(1)%molecule%SiteCharge(2)%RY(4) = 0.165365
!    this%Component(1)%molecule%SiteCharge(2)%RZ(4) = 0.227414
!    this%Component(1)%molecule%SiteCharge(3)%RX(4) = 0.275015
!    this%Component(1)%molecule%SiteCharge(3)%RY(4) = 0.065446
!    this%Component(1)%molecule%SiteCharge(3)%RZ(4) = 0.270726


! Dimensions of charge
   fac = UnitCharge / ElementaryCharge

   DO i=1,this%NComponents
     DO j=1,this%Component(i)%Molecule%NCharge
      DO k=1,this%Component(i)%NPart

       q  = this%Component(i)%Molecule%SiteCharge(j)%e
       RX = this%Component(i)%Molecule%SiteCharge(j)%RX(k)
       RY = this%Component(i)%Molecule%SiteCharge(j)%RY(k)
       RZ = this%Component(i)%Molecule%SiteCharge(j)%RZ(k)

       RX     = RX - anint(RX)
       RY     = RY - anint(RY)
       RZ     = RZ - anint(RZ)

       RXgit1 = this%gridx*(RX+0.5)
       RYgit1 = this%gridy*(RY+0.5)
       RZgit1 = this%gridz*(RZ+0.5)

       RXgit = (RXgit1 - floor( RXgit1 ))
       RYgit = (RYgit1 - floor( RYgit1 ))
       RZgit = (RZgit1 - floor( RZgit1 ))

! x-coordinate
       err=fillspline(RXgit,transf)
         splinex  = transf

! y-coordinate
       err=fillspline(RYgit,transf)
         spliney  = transf

! z-coordinate
       err=fillspline(RZgit,transf)
         splinez  = transf

       xxo = floor (RXgit1 - this%splineorder)
       yyo = floor (RYgit1 - this%splineorder)
       zzo = floor (RZgit1 - this%splineorder)

       if (xxo > NX-1) xxo = NX-1
       if (yyo > NY-1) yyo = NY-1
       if (zzo > NZ-1) zzo = Nz-1

       DO xi=1,order,1
         xxo = xxo + 1
         x = xxo
         if (x < 0) x = x + NX
         yo   = yyo
         DO yi=1,order,1
           yo = yo + 1
           y   = yo
           if (y < 0) y = y + NY
           zo  = zzo
           DO zi=1,order,1
             zo = zo + 1
             z  = zo
             if (z < 0) z = z + NZ
             index_loc = z + y*NZ + x*NY*NZ + 1
             this%qgrida(1,index_loc) = this%qgrida(1,index_loc) + splinex(xi)*spliney(yi)*splinez(zi)*q*fac
           END DO
         END DO
       END DO

      END DO ! Particles
     END DO  ! Charges
   END DO ! Components


contains

!!!!!!!!!!!!!!!!!!
! Calculate Spline function and Derivative of the spline function
    real(RK) function fillspline(w,spline)

    implicit none
      real(RK),intent(in) ::  w
      real(RK),dimension(this%splineorder+5),intent(in out) ::  spline


      real(RK)            :: div
      real(RK)            :: dj,di,dorder
      integer             :: i,j,order

      order = this%splineorder

! Calculate the single spline function contributions
      spline    = 0._RK


      spline(1) = 1._RK - w
      spline(2) = w

      DO i=3,order-1,1
        div = 1._RK / (i-1)
        spline(i) = div * w * spline(i-1)
        di = real(i)
! Second Contribution
         DO j=1,i-2,1
           dj = real(j)
           spline(i-j) = div * ((w+dj)*spline(i-j-1) + &
&                                   (di-dj-w)*spline(i-j))
         END DO
       spline(1) = spline(1) * div * (1._RK-w)
       END DO

! Generate order spline
    div = 1._RK / (order - 1._RK)
    spline(order) = div * w * spline(order-1)
    dorder = real(order)
    DO i=1,order-2,1
      di = real(i)
      spline(order- i) = div*((w+di)*spline(order-i-1) + &
&                                  (dorder-di-w)*spline(order-i))
    END DO
    spline(1) = spline(1) * div* (1._RK-w)

    end function fillspline

   end subroutine charge_grid_MCall



   subroutine TEnsemble_PMChargeGrid_plus(this,nc,np)

   implicit none
   type(TEnsemble)    :: this
   integer,intent(in) :: nc
   integer,intent(in) :: np


   real(RK)           :: boxl
   real(RK),dimension(this%splineorder+5) ::  transf
   real(RK),dimension(this%splineorder+5) ::  splinex
   real(RK),dimension(this%splineorder+5) ::  spliney
   real(RK),dimension(this%splineorder+5) ::  splinez

   real(RK)           :: RX,RY,RZ
   real(RK)           :: RXgit ,RYgit, RZgit
   real(RK)           :: RXgit1,RYgit1,RZgit1
   real(RK)           :: q

   integer            :: xxo,yyo,zzo
   integer            :: yo,zo

   integer            :: i,j,k
   integer            :: xi,yi,zi
   integer            :: x,y,z
   integer            :: order
   integer            :: NX,NY,NZ
   integer            :: index_loc

   real(RK)           :: err
   real(RK)           :: fac


!    nullify(this%qgrida)
!    write(*,*) 'nicht ausgenullt!'
   NX = this%gridx
   NY = this%gridy
   NZ = this%gridz

   boxl  = this%BoxLength
   order = this%splineorder

!    this%qgrida(2,1:NX*NY*NZ+1) = 0._RK
! Dimensions of charge
   fac = UnitCharge / ElementaryCharge

   DO j=1,this%Component(nc)%Molecule%NCharge

       q  = this%Component(nc)%Molecule%SiteCharge(j)%e
       RX = this%Component(nc)%Molecule%SiteCharge(j)%RX(np)
       RY = this%Component(nc)%Molecule%SiteCharge(j)%RY(np)
       RZ = this%Component(nc)%Molecule%SiteCharge(j)%RZ(np)

       RX     = RX - anint(RX)
       RY     = RY - anint(RY)
       RZ     = RZ - anint(RZ)

       RXgit1 = this%gridx*(RX+0.5)
       RYgit1 = this%gridy*(RY+0.5)
       RZgit1 = this%gridz*(RZ+0.5)

       RXgit = (RXgit1 - floor( RXgit1 ))
       RYgit = (RYgit1 - floor( RYgit1 ))
       RZgit = (RZgit1 - floor( RZgit1 ))

! x-coordinate
       err=fillspline(RXgit,transf)
         splinex  = transf

! y-coordinate
       err=fillspline(RYgit,transf)
         spliney  = transf

! z-coordinate
       err=fillspline(RZgit,transf)
         splinez  = transf

       xxo = floor (RXgit1 - this%splineorder)
       yyo = floor (RYgit1 - this%splineorder)
       zzo = floor (RZgit1 - this%splineorder)

       if (xxo > NX-1) xxo = NX-1
       if (yyo > NY-1) yyo = NY-1
       if (zzo > NZ-1) zzo = Nz-1

       DO xi=1,order,1
         xxo = xxo + 1
         x = xxo
         if (x < 0) x = x + NX
         yo   = yyo
         DO yi=1,order,1
           yo = yo + 1
           y   = yo
           if (y < 0) y = y + NY
           zo  = zzo
           DO zi=1,order,1
             zo = zo + 1
             z  = zo
             if (z < 0) z = z + NZ
             index_loc = z + y*NZ + x*NY*NZ + 1
             this%qgrida(1,index_loc) = this%qgrida(1,index_loc) + splinex(xi)*spliney(yi)*splinez(zi)*q*fac
           END DO
         END DO
       END DO

     END DO  ! Charges


contains

!!!!!!!!!!!!!!!!!!
! Calculate Spline function and Derivative of the spline function
    real(RK) function fillspline(w,spline)

    implicit none
      real(RK),intent(in) ::  w
      real(RK),dimension(this%splineorder+5),intent(in out) ::  spline


      real(RK)            :: div
      real(RK)            :: dj,di,dorder
      integer             :: i,j,order

      order = this%splineorder

! Calculate the single spline function contributions
      spline    = 0._RK


      spline(1) = 1._RK - w
      spline(2) = w

      DO i=3,order-1,1
        div = 1._RK / (i-1)
        spline(i) = div * w * spline(i-1)
        di = real(i)
! Second Contribution
         DO j=1,i-2,1
           dj = real(j)
           spline(i-j) = div * ((w+dj)*spline(i-j-1) + &
&                                   (di-dj-w)*spline(i-j))
         END DO
       spline(1) = spline(1) * div * (1._RK-w)
       END DO

! Generate order spline
    div = 1._RK / (order - 1._RK)
    spline(order) = div * w * spline(order-1)
    dorder = real(order)
    DO i=1,order-2,1
      di = real(i)
      spline(order- i) = div*((w+di)*spline(order-i-1) + &
&                                  (dorder-di-w)*spline(order-i))
    END DO
    spline(1) = spline(1) * div* (1._RK-w)

    end function fillspline

   end subroutine TEnsemble_PMChargeGrid_plus



   subroutine TEnsemble_PMChargeGrid_min(this,nc,np)

   implicit none
   type(TEnsemble)    :: this
   integer,intent(in) :: nc
   integer,intent(in) :: np

   real(RK)           :: boxl
   real(RK),dimension(this%splineorder+5) ::  transf
   real(RK),dimension(this%splineorder+5) ::  splinex
   real(RK),dimension(this%splineorder+5) ::  spliney
   real(RK),dimension(this%splineorder+5) ::  splineZ

   real(RK)           :: RX,RY,RZ
   real(RK)           :: RXgit ,RYgit, RZgit
   real(RK)           :: RXgit1,RYgit1,RZgit1
   real(RK)           :: q

   integer            :: xxo,yyo,zzo
   integer            :: yo,zo

   integer            :: i,j,k
   integer            :: xi,yi,zi
   integer            :: x,y,z
   integer            :: order
   integer            :: NX,NY,NZ
   integer            :: index_loc

   real(RK)           :: err
   real(RK)           :: fac


!    nullify(this%qgrida)
!    write(*,*) 'nicht ausgenullt!'
   NX = this%gridx
   NY = this%gridy
   NZ = this%gridz

   boxl  = this%BoxLength
   order = this%splineorder

! Dimensions of charge
   fac = UnitCharge / ElementaryCharge

   DO j=1,this%Component(nc)%Molecule%NCharge

       q  = this%Component(nc)%Molecule%SiteCharge(j)%e
       RX = this%Component(nc)%Molecule%SiteCharge(j)%RX(np)
       RY = this%Component(nc)%Molecule%SiteCharge(j)%RY(np)
       RZ = this%Component(nc)%Molecule%SiteCharge(j)%RZ(np)

       RX     = RX - anint(RX)
       RY     = RY - anint(RY)
       RZ     = RZ - anint(RZ)

       RXgit1 = this%gridx*(RX+0.5)
       RYgit1 = this%gridy*(RY+0.5)
       RZgit1 = this%gridz*(RZ+0.5)

       RXgit = (RXgit1 - floor( RXgit1 ))
       RYgit = (RYgit1 - floor( RYgit1 ))
       RZgit = (RZgit1 - floor( RZgit1 ))

! x-coordinate
       err=fillspline(RXgit,transf)
         splinex  = transf

! y-coordinate
       err=fillspline(RYgit,transf)
         spliney  = transf

! z-coordinate
       err=fillspline(RZgit,transf)
         splinez  = transf

       xxo = floor (RXgit1 - this%splineorder)
       yyo = floor (RYgit1 - this%splineorder)
       zzo = floor (RZgit1 - this%splineorder)

       if (xxo > NX-1) xxo = NX-1
       if (yyo > NY-1) yyo = NY-1
       if (zzo > NZ-1) zzo = Nz-1

       DO xi=1,order,1
         xxo = xxo + 1
         x = xxo
         if (x < 0) x = x + NX
         yo   = yyo
         DO yi=1,order,1
           yo = yo + 1
           y   = yo
           if (y < 0) y = y + NY
           zo  = zzo
           DO zi=1,order,1
             zo = zo + 1
             z  = zo
             if (z < 0) z = z + NZ
             index_loc = z + y*NZ + x*NY*NZ + 1
             this%qgrida(1,index_loc) = this%qgrida(1,index_loc) - splinex(xi)*spliney(yi)*splinez(zi)*q*fac
           END DO
         END DO
       END DO

     END DO  ! Charges


contains

!!!!!!!!!!!!!!!!!!
! Calculate Spline function and Derivative of the spline function
    real(RK) function fillspline(w,spline)

    implicit none
      real(RK),intent(in) ::  w
      real(RK),dimension(this%splineorder+5),intent(in out) ::  spline


      real(RK)            :: div
      real(RK)            :: dj,di,dorder
      integer             :: i,j,order

      order = this%splineorder

! Calculate the single spline function contributions
      spline    = 0._RK


      spline(1) = 1._RK - w
      spline(2) = w

      DO i=3,order-1,1
        div = 1._RK / (i-1)
        spline(i) = div * w * spline(i-1)
        di = real(i)
! Second Contribution
         DO j=1,i-2,1
           dj = real(j)
           spline(i-j) = div * ((w+dj)*spline(i-j-1) + &
&                                   (di-dj-w)*spline(i-j))
         END DO
       spline(1) = spline(1) * div * (1._RK-w)
       END DO

! Generate order spline
    div = 1._RK / (order - 1._RK)
    spline(order) = div * w * spline(order-1)
    dorder = real(order)
    DO i=1,order-2,1
      di = real(i)
      spline(order- i) = div*((w+di)*spline(order-i-1) + &
&                                  (dorder-di-w)*spline(order-i))
    END DO
    spline(1) = spline(1) * div* (1._RK-w)

    end function fillspline

   end subroutine TEnsemble_PMChargeGrid_min


   subroutine TEnsemble_VirialIntra(this)

   implicit none
   type(TEnsemble)            :: this
   type(TMolecule),pointer    :: pm

   integer  :: i,j,jj,k
   real(RK) :: RX,RY,RZ
   real(RK) :: PX,PY,PZ
   real(RK) :: dx,dy,dz,dr
   real(RK) :: ex,ey,ez
   real(RK) :: qj,qjj
   real(RK) :: virloc
   real(RK) :: Kappa
   real(RK) :: kapparij,approx
   real(RK) :: fij
   real(RK) :: Faktor

   virloc = 0._RK
   Kappa  = this%Kappa
   Faktor = 2./sqrt(Pi) * Kappa


   DO i=1,this%NComponents
     pm => this%Component(i)%Molecule
     DO j=1,pm%NCharge
      qj = pm%SiteCharge(j)%e

      DO jj=1,j-1
       qjj = pm%SiteCharge(jj)%e
       DO k=1,this%Component(i)%NPart
         RX = pm%SiteCharge(j)%RX(k)
         RY = pm%SiteCharge(j)%RY(k)
         RZ = pm%SiteCharge(j)%RZ(k)

         PX = pm%SiteCharge(jj)%RX(k)
         PY = pm%SiteCharge(jj)%RY(k)
         PZ = pm%SiteCharge(jj)%RZ(k)

         dx = (RX-PX) * this%BoxLength
         dy = (RY-PY) * this%BoxLength
         dz = (RZ-PZ) * this%BoxLength

         dr = sqrt(dx**2 + dy**2 + dz**2)

         KappaRij = Kappa*dr
         call erfc_approx(KappaRij,approx)

         eX = dx / dr
         eY = dy / dr
         eZ = dz / dr

         Fij  = (qj*qjj/dr*(1._RK-approx) - Faktor*exp(-KappaRij**2)*qj*qjj) /dr

 	 Virloc = Virloc + Fij* (eX * dx + eY * dy + eZ * dz)
        END DO
      ENd DO

      DO jj=j+1,pm%NCharge
       qjj = pm%SiteCharge(jj)%e
       DO k=1,this%Component(i)%NPart
         RX = pm%SiteCharge(j)%RX(k)
         RY = pm%SiteCharge(j)%RY(k)
         RZ = pm%SiteCharge(j)%RZ(k)

         PX = pm%SiteCharge(jj)%RX(k)
         PY = pm%SiteCharge(jj)%RY(k)
         PZ = pm%SiteCharge(jj)%RZ(k)

         dx = (RX-PX) * this%BoxLength
         dy = (RY-PY) * this%BoxLength
         dz = (RZ-PZ) * this%BoxLength

         dr = sqrt(dx**2 + dy**2 + dz**2)

         KappaRij = Kappa*dr
         call erfc_approx(KappaRij,approx)

         eX = dx / dr
         eY = dy / dr
         eZ = dz / dr

         Fij  = (qj*qjj/dr*(1._RK-approx) - Faktor*exp(-KappaRij**2)*qj*qjj) /dr

 	 Virloc = Virloc + Fij* (eX * dx + eY * dy + eZ * dz)
        END DO
       END DO
     END DO
   END DO

   end subroutine TEnsemble_VirialIntra

end module ms2_ensemble
