!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 1.0                !
!  (c) 2011 by TU Kaiserslautern                               !
!      P.O. Box 67653                                          !
!      67653 Kaiserslautern                                    !
!==============================================================!
!  Module ms2_ensemble                                         !
!  Contains TEnsemble object                                   !
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
!DEC$ MESSAGE:'Compiling ms2_ensemble.F90...'
#endif

module ms2_ensemble

!#ifdev MPI_VER > 0
!  use mpi
!#endif

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
#if  TRANS == 1
    ! I/O unit for result ACF
    integer :: iounit_rescf   !TRANSPORT_thisline
#endif
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

    ! Optional calculation of pressure
    logical :: OptPressure

    ! Positions and orientations of test particles
    real(RK), pointer :: P0Test(:, :), Q0Test(:, :)

    ! Number of components in ensemble
    integer :: NComponents, NRealComponents, NGradInsComp

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
#if  TRANS == 1
!TRANSPORT_start
    ! Correlation functions

    integer :: Ncorr,Mmess,MmessMax
    integer :: NSpancf,Nviewcf
    real(RK), pointer :: cf_db(:), cf_vs(:), cf_vb(:), cf_c(:)
    real(RK), pointer :: lamda(:, :)
!    real(RK), pointer :: cf_vs_kk(:), cf_vs_kp(:), cf_vs_pp(:) ! FIX_ME these seem to be artifacts from an old
!    real(RK), pointer :: cf_ct_kk(:), cf_ct_kp(:), cf_ct_pp(:) ! and obsolete version
!    real(RK), pointer :: cf_cr_kk(:), cf_cr_kp(:), cf_cr_pp(:)
    real(RK), pointer :: sinte_i(:, :), sinte_lamda(:,:)
    real(RK), pointer :: sinte_db(:), sinte_vs(:), sinte_vb(:)
    real(RK), pointer :: sinte_c(:)
!    real(RK), pointer :: sinte_vs_kk(:), sinte_vs_kp(:), sinte_vs_pp(:) ! FIX_ME these seem to be artifacts from an old
!    real(RK), pointer :: sinte_ct_kk(:), sinte_ct_kp(:), sinte_ct_pp(:) ! and obsolete version
!    real(RK), pointer :: sinte_cr_kk(:), sinte_cr_kp(:), sinte_cr_pp(:)
    real(RK), pointer :: a(:, :), cf_d (:, :), vsk(:, :), vsp(:, :), vbk(:, :), vbp(:, :)
    real(RK), pointer :: vckt(:, :), vckr(:, :), vcpt(:, :), vcpr(:, :)
    real(RK)          :: sc(3),sp(3)

    real(RK),pointer :: selfd_i(:)
    real(RK)         :: ternary_a, ternary_b, ternary_c
    real(RK)         :: binary_d
    real(RK)         :: visco_s
!    real(RK)         :: visco_s_kk, visco_s_kp, visco_s_pp
    real(RK)         :: visco_b
    real(RK)         :: conduct
!    real(RK)         :: conduct_t_kk, conduct_t_kp, conduct_t_pp ! FIX_ME these seem to be artifacts from an old
!    real(RK)         :: conduct_r_kk, conduct_r_kp, conduct_r_pp ! and obsolete version

    ! 4.) Transport properties

    type(TAccumulatorCF),pointer :: Sumself_i(:)
    type(TAccumulatorCF)         :: SumTer_a, SumTer_b, SumTer_c
    type(TAccumulatorCF)         :: SumBin_d
    type(TAccumulatorCF)         :: SumVisco_s
!    type(TAccumulatorCF)         :: SumVisco_s_kk  ! FIX_ME these seem to be artifacts from an old
!    type(TAccumulatorCF)         :: SumVisco_s_kp  ! and obsolete version
!    type(TAccumulatorCF)         :: SumVisco_s_pp
    type(TAccumulatorCF)         :: SumVisco_b
    type(TAccumulatorCF)         :: SumConduct
!    type(TAccumulatorCF)         :: SumConduct_t_kk
!    type(TAccumulatorCF)         :: SumConduct_t_kp
!    type(TAccumulatorCF)         :: SumConduct_t_pp
!    type(TAccumulatorCF)         :: SumConduct_r_kk
!    type(TAccumulatorCF)         :: SumConduct_r_kp
!    type(TAccumulatorCF)         :: SumConduct_r_pp
!TRANSPORT_END
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
#if  TRANS == 1
!TRANSPORT_start
  interface CalCorrFun
    module procedure TEnsemble_CalCorrFun
  end interface

  interface IntCorrFun
    module procedure TEnsemble_IntCorrFun
  end interface
!TRANSPORT_END
#endif
contains



!==============================================================!
!  Subroutine TEnsemble_Construct                              !
!==============================================================!

  subroutine TEnsemble_Construct( this, ne )

    implicit none

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: ne

    ! Declare local variables
    integer :: i, j
    integer :: stat
    character( IOBufferLength ) :: str

    ! Allocate simulation box length
    allocate( this%BoxLength, STAT = stat )
    call AllocationError( stat, 'simulation box length' )

    ! Allocate maximum number of particles
    allocate( this%NPartMax, STAT = stat )
    call AllocationError( stat, 'maximum number of particles' )

    ! Set number of ensemble
    this%EnsembleNumber = ne
    call LogWriteBlank
    write( IOBuffer, '(72(1H-))')
    call LogWrite
    write( IOBuffer, '(T14, "Reading parameters of ensemble", I3)' ) &
&     this%EnsembleNumber
    call LogWrite

    ! Read temperature
    call FileReadParameter( this%RefTemperature, iounit_params , IdRefTemperature, .false. )
    if( .not. UseReducedUnits ) then
      this%RefTemperature = this%RefTemperature / UnitTemperature
    end if

    ! Read pressure
    if( EnsembleType .eq. EnsembleTypeGE ) then
      call FileReadParameter( this%RefPressure, iounit_params , IdPressure0, .false. )
      if( .not. UseReducedUnits ) then
        this%RefPressure = this%RefPressure * 1E6_RK / UnitPressure
      end if
    end if
    if( ConstantPressure ) then
      call FileReadParameter( this%RefPressure, iounit_params , IdRefPressure, .false. )
      if( .not. UseReducedUnits ) then
        this%RefPressure = this%RefPressure * 1E6_RK / UnitPressure
      end if
    end if

    ! Read liquid simulation data
    if( EnsembleType .eq. EnsembleTypeGE ) then
      call FileReadParameter( this%LiqDensity, iounit_params , IdLiqDensity, .false. )
      if( .not. UseReducedUnits ) then
        this%LiqDensity = this%LiqDensity / UnitDensity
      end if
      call FileReadParameter( this%VarLiqDensity, iounit_params , IdVarLiqDensity, .false. )
      if( .not. UseReducedUnits ) then
        this%VarLiqDensity = this%VarLiqDensity / UnitDensity
      end if
      call FileReadParameter( this%LiqEnthalpy, iounit_params , IdLiqEnthalpy, .false. )
      if( .not. UseReducedUnits ) then
        this%LiqEnthalpy = this%LiqEnthalpy / ( UnitEnergy * NAvogadro )
      end if
      call FileReadParameter( this%VarLiqEnthalpy, iounit_params , IdVarLiqEnthalpy, .false. )
      if( .not. UseReducedUnits ) then
        this%VarLiqEnthalpy = this%VarLiqEnthalpy / ( UnitEnergy * NAvogadro )
      end if
      call FileReadParameter( this%LiqBetaT, iounit_params , IdLiqBetaT, .false. )
      if( .not. UseReducedUnits ) then
        this%LiqBetaT = this%LiqBetaT * UnitPressure * 1E-6_RK
      end if
      call FileReadParameter( this%VarLiqBetaT, iounit_params , IdVarLiqBetaT, .false. )
      if( .not. UseReducedUnits ) then
        this%VarLiqBetaT = this%VarLiqBetaT * UnitPressure * 1E-6_RK
      end if
      call FileReadParameter( this%LiqdHdP, iounit_params , IdLiqdHdP, .false. )
      if( .not. UseReducedUnits ) then
        this%LiqdHdP = this%LiqdHdP * UnitDensity
      end if
      call FileReadParameter( this%VarLiqdHdP, iounit_params , IdVarLiqdHdP, .false. )
      if( .not. UseReducedUnits ) then
        this%VarLiqdHdP = this%VarLiqdHdP * UnitDensity
      end if
    end if

    ! Read density
    call FileReadParameter( this%RefDensity, iounit_params , IdRefDensity, .false. )
    if( .not. UseReducedUnits ) then
      this%RefDensity = this%RefDensity / UnitDensity
    end if

    ! Update log file
    write( IOBuffer, '("Temperature: ",T26, F9.3, " K")' ) &
&     this%RefTemperature * UnitTemperature
    call LogWrite
    if( ConstantPressure ) then
      write( IOBuffer, '("Pressure: ",T26, F9.3, " MPa")' ) &
&       this%RefPressure * UnitPressure * 1E-6_RK
      call LogWrite
    end if
    if( EnsembleType .eq. EnsembleTypeGE ) then
      write( IOBuffer, '("Pressure0: ",T26, F12.6, " MPa")' ) &
&       this%RefPressure * UnitPressure * 1E-6_RK
      call LogWrite
      write( IOBuffer, '("Liquid density: ",T26, F12.6, " (", F13.6, ") mol/l")' ) &
&       this%LiqDensity * UnitDensity, this%VarLiqDensity * UnitDensity
      call LogWrite
      write( IOBuffer, '("Liquid enthalpy: ",T22, F16.6, " (", F13.6, ") J/mol")' ) &
&       this%LiqEnthalpy * UnitEnergy * NAvogadro, &
&       this%VarLiqEnthalpy * UnitEnergy * NAvogadro
      call LogWrite
      write( IOBuffer, '("Liquid betaT: ",T26, F12.6, " (", F13.6, ") 1/MPa")' ) &
&       this%LiqBetaT / UnitPressure * 1E6_RK, &
&       this%VarLiqBetaT / UnitPressure * 1E6_RK
      call LogWrite
      write( IOBuffer, '("Liquid dHdP: ",T26, F12.6, " (", F13.6, ") l/mol")' ) &
&       this%LiqdHdP / UnitDensity, this%VarLiqdHdP / UnitDensity
      call LogWrite
    end if
    write( IOBuffer, '("Density: ", T26, F9.3, " mol/l")' ) &
&     this%RefDensity * UnitDensity
    call LogWrite
    write( IOBuffer, '("Reduced temperature: ", T26, F12.6)' ) this%RefTemperature
    call LogWrite
    if( ConstantPressure ) then
      write( IOBuffer, '("Reduced pressure: ",T26, F12.6)' ) this%RefPressure
      call LogWrite
    end if
    if( EnsembleType .eq. EnsembleTypeGE ) then
      write( IOBuffer, '("Reduced pressure0: ",T26, F12.6)' ) this%RefPressure
      call LogWrite
      write( IOBuffer, &
&       '("Reduced liquid density: ",T29, F9.6, " (", F13.6, ")")' ) &
&       this%LiqDensity, this%VarLiqDensity
      call LogWrite
      write( IOBuffer, &
&       '("Reduced liquid enthalpy: ", T28, F10.6, " (", F13.6, ")")' ) &
&       this%LiqEnthalpy, this%VarLiqEnthalpy
      call LogWrite
      write( IOBuffer, '("Reduced liquid betaT: ",T26, F12.6, " (", F13.6, ")")' ) &
&       this%LiqBetaT, this%VarLiqBetaT
      call LogWrite
      write( IOBuffer, '("Reduced liquid dHdP: ",T26, F12.6, " (", F13.4, ")")' ) &
&       this%LiqdHdP, this%VarLiqdHdP
      call LogWrite
    end if
    write( IOBuffer, '("Reduced density: ",T26, F12.6)' ) this%RefDensity
    call LogWrite

    ! Read mass of piston
    if( SimulationType .eq. MolecularDynamics .and. ConstantPressure ) then
      call FileReadParameter( this%PistonMass, iounit_params , IdPistonMass, .false. )
      if( .not. UseReducedUnits ) then
!        this%PistonMass = this%PistonMass / UnitMass * UnitLength**4
      end if
      write( IOBuffer, '("Mass of piston: ",T26, F15.9)' ) this%PistonMass
      call LogWrite
    end if

    ! Read optional pressure calculation
    this%OptPressure = .true.
    if( SimulationType .eq. MonteCarlo ) then
      call FileReadParameter( str, iounit_params , IdOptPressure, .false., "yes" )
      select case( str )
        case( 'YES', 'Yes', 'yes' )
          this%OptPressure = .true.
          write( IOBuffer, '("Pressure calculation: ",T30, A)' ) trim( str )
          call LogWrite
        case( 'NO', 'No', 'no')
          this%OptPressure = .false.
          write( IOBuffer, '("Pressure calculation: ",T30, A)' ) trim( str )
          call LogWrite
        case default
          call Error( 'Select yes/no for calculation of pressure '// &
&           ProgramFileName//ConfigFileExtension )
      end select
      if ( .not. ConstantPressure .and. .not. this%OptPressure) &
&          call Error( 'Pressure Calculation in NVT necessary' )
    end if
    if ( EnsembleType .eq. EnsembleTypeGE .and. .not. this%OptPressure ) then
      write(IOBuffer, '("For GE simulations, please set Logical OptPressure to yes")' )
      call LogWrite
      call Error( ' ms2 has to quit' )
    end if

   ! if( SimulationType .eq. MonteCarlo ) then
      ! Read whether to perform the MC equilibration in parallel
      call FileReadParameter( str, iounit_params , IdCommonEqui, .false., "yes" )
      select case( str )
        case( 'YES', 'Yes', 'yes' )
          CommonEqui = .true.
          if( SimulationType .eq. MonteCarlo ) then
            write( IOBuffer, '("Common equilibration: ",T30, A)' ) trim( str )
            call LogWrite
          endif
        case( 'NO', 'No', 'no')
          CommonEqui = .false.
          write( IOBuffer, '("Common equilibration: ",T30, A)' ) trim( str )
          call LogWrite
        case default
          call Error( 'Select yes/no for common equilibration '// &
&           ProgramFileName//ConfigFileExtension )
      end select
    !endif
  
    ! Read initial number of particles in ensemble
    call FileReadParameter( this%NPart, iounit_params , IdNPart, .false. )
    if( EnsembleType .eq. EnsembleTypeGE .or. &
&       EnsembleType .eq. EnsembleTypeHA ) then
      this%NPartInitial = this%NPart
      this%NPartLBound = int( real( this%NPart, RK ) / 1.2_RK )
      this%NPartUBound = int( real( this%NPart, RK ) * 1.2_RK )
    end if

    ! Read number of components in ensemble
    call FileReadParameter( this%NComponents, iounit_params , IdNComponents, .false. )
    write( IOBuffer, '("Number of components:",T28, I3)' ) this%NComponents
    call LogWrite
    if( this%NComponents <= 0 ) then
      write( ErrorBuffer, &
&       '("There must be at least 1 component in ensemble", I2)' ) &
&       this%EnsembleNumber
      call Error
    end if
    if( this%NComponents > 999 ) &
&     call Error( 'Cannot work with more than 999 components on '//Hardware )
#if  TRANS == 1
!TRANSPORT_start
    ! Read correlation function
    if( CorrfunMode .eq. active ) then

      ! Read legth of the correlation function
      call FileReadParameter( this%NCorr , iounit_params , IdCorrlength )
      ! Read time span between correlations
      call FileReadParameter( this%NSpanCF , iounit_params , IdSpanCF )
 
      if(mod(this%NCorr, this%NSpanCF) .eq. 0) then
        write( IOBuffer, '("Length of CorrFunction:",T26, I5)' ) this%NCorr
        call LogWrite
      else
        this%NCorr = (AINT(real( this%NCorr, RK )/real( this%NSpanCF, RK ))+1)*this%NSpanCF
        write( IOBuffer, '("Length of CorrFunction is extended to:",T40, I7)') this%NCorr
        call LogWrite
      endif

      write( IOBuffer, '("Time Span between cf:",T26, I5)' ) this%NSpanCF
      call LogWrite

      call FileReadParameter( this%Nviewcf , iounit_params , IdNviewcf )
      write( IOBuffer, '("Print cf each:",T26, I5)' ) this%Nviewcf
      call LogWrite

      ! Read frequency of updating result file CF
      call FileReadParameter( BlockSizeCF , iounit_params , IdBlockSizeCF, .false., 0 )
      if( BlockSize > 0 ) then
        write( IOBuffer, &
&         '("Result files will be updated each", I7, " Correlation Functions")' ) &
&         BlockSizeCF
      else
        write( IOBuffer, '("Result files will not be created")' )
      end if
      call LogWrite

      ! Calculate MmessMax
      this%MmessMax = int((NSteps-this%Ncorr)/this%NSpancf)

      if( BlockSizeCF > 0 ) then
        NBlocksMaxCF = this%MmessMax / BlockSizeCF
        NBlockSizesMaxCF = int( sqrt( real( NBlocksMaxCF, RK ) ) )
      else
        NBlocksMaxCF = 0
        NBlockSizesMaxCF = 0
      end if

    end if
!TRANSPORT_END
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
        call FileReadParameter( this%ScaleSigma(i, j), iounit_params , IdScaleSigma, .false. )
        if( i /= j ) this%ScaleSigma(j, i) = this%ScaleSigma(i, j)
        call FileReadParameter( this%ScaleEpsilon(i, j), iounit_params , IdScaleEpsilon, .false. )
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
      call FileReadParameter( this%RCutoffLJ126LJ126, iounit_params , IdRCutoffCOM, .false. )
      if (this%RCutoffLJ126LJ126 < 0._RK) then
        this%RCutoffLJ126LJ126 = 0.9*0.5*(this%NPart / &
&          (NAvogadro*this%RefDensity*UnitDensity*1000))**(1._RK/3._RK)/UnitLength
      end if
      call LogWriteBlank
      write( IOBuffer, '("Reduced center of mass cutoff radius: ",T45, F6.3)' ) &
&       this%RCutoffLJ126LJ126
      call LogWrite
      this%RCutoffDipoleDipole = this%RCutoffLJ126LJ126
      this%RCutoffDipoleQuadrupole = this%RCutoffLJ126LJ126
      this%RCutoffQuadrupoleQuadrupole = this%RCutoffLJ126LJ126
    else
      if( this%NLJ126Max > 0 ) then
        call FileReadParameter( this%RCutoffLJ126LJ126, iounit_params , IdRCutoffLJ126LJ126, .false. )
        write( IOBuffer, &
&         '("Lennard-Jones cutoff radius: ",T45, F6.3, " sigma")' ) &
&         this%RCutoffLJ126LJ126
        call LogWrite
      end if
      if( this%NDipoleMax > 0 ) then
        call FileReadParameter( this%RCutoffDipoleDipole, iounit_params , IdRCutoffDipoleDipole, .false. )
        write( IOBuffer, '("Reduced dipole-dipole cutoff radius: ",T42, F8.3)' ) &
&         this%RCutoffDipoleDipole
        call LogWrite
        if( this%NQuadrupoleMax > 0 ) then
          call FileReadParameter( this%RCutoffDipoleQuadrupole, iounit_params , IdRCutoffDipoleQuadrupole, .false. )
          write( IOBuffer, &
&           '("Reduced dipole-quadrupole cutoff radius: ",T42, F8.3)' ) &
&           this%RCutoffDipoleQuadrupole
          call LogWrite
        end if
      end if
      if( this%NQuadrupoleMax > 0 ) then
        call FileReadParameter( this%RCutoffQuadrupoleQuadrupole, iounit_params , IdRCutoffQuadrupoleQuadrupole, .false. )
        write( IOBuffer, &
&         '("Reduced quadrupole-quadrupole cutoff radius: ",T42, F8.3)' ) &
&         this%RCutoffQuadrupoleQuadrupole
        call LogWrite
      end if
    end if

    ! Read characteristic dielectric constant
    this%RFEpsilon = 0._RK
    if(( this%NDipoleMax > 0 ) .or. ( this%NChargeMax > 0 )) then
      call LogWriteBlank
      call FileReadParameter( this%RFEpsilon, iounit_params , IdRFEpsilon, .false. )
      write( IOBuffer, '("Characteristic dielectric constant: ",T41, E16.5)' ) &
&       this%RFEpsilon
      call LogWrite
    end if

    ! Create potentials
    call CreatePotentials( this )

    ! Calculate long-range corrections
    call CalculateCorr( this )
    call LogWriteBlank
    write( IOBuffer, '("Cutoff correction to")' )
    call LogWrite
    
    if ( SimulationType .eq. MonteCarlo .and. (.not.(Equilibration .and. CommonEqui)) ) then 
      write( IOBuffer, &
&       '("- potential energy from LJ",T44, F12.8)' ) &
&       this%EPotCorrLJ / this%NPart
    else
      write( IOBuffer, &
&       '("- potential energy from LJ",T44, F12.8)' ) &
&       this%EPotCorrLJ * Nprocs/ this%NPart
    endif
        
    call LogWrite
    if ( SimulationType .eq. MonteCarlo .and. (.not.(Equilibration .and. CommonEqui)) ) then  
      write( IOBuffer, &
&       '("- pressure from LJ ",T44, F12.8)' ) &
&       this%VirialCorrLJ / this%NPart
    else
      write( IOBuffer, &
&       '("- pressure from LJ ",T44, F12.8)' ) &
&       this%VirialCorrLJ * NProcs / this%NPart
    endif
    
    call LogWrite

    do i = 1, this%NRealComponents
      write( IOBuffer, &
&       '("- chem. pot. of ", A, " from LJ",T44, F12.8)' ) &
&       trim( this%Component(i)%PotModFileName ), &
&       this%Component(i)%EPotTestCorrLJ
      call LogWrite
    end do
    
    if ( SimulationType .eq. MonteCarlo .and. (.not.(Equilibration .and. CommonEqui)) ) then 
      write( IOBuffer, &
&       '("- potential energy from reaction field (RF)",T44, F12.8)' ) &
&       this%EPotCorrRF / this%NPart
    else
      write( IOBuffer, &
&       '("- potential energy from reaction field (RF)",T44, F12.8)' ) &
&       this%EPotCorrRF * NProcs / this%NPart
    endif
    call LogWrite

    do i = 1, this%NRealComponents
      write( IOBuffer, &
&       '("- chem. pot. of ", A, " from RF",T44, F12.8)' ) &
&       trim( this%Component(i)%PotModFileName ), &
&       this%Component(i)%EPotTestCorrRF
      call LogWrite
    end do

    ! Calculate maximum cutoff radius
!    this%RCutoffMax2 = 0._RK
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

      ! Abort, if maximum cutoff larger than boxlength 
      if (this%RCutoffMax2 > this%BoxLength) &
&       call Error('Cutoff is larger than the boxsize')

      ! Set initial positions of particles in simulation box
      call InitPositions( this )

      ! Set initial orientations of particles in simulation box
      call InitOrientations( this )

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

    end if

    ! Set I/O unit numbers
    i = FilesPerEnsemble * this%EnsembleNumber
    this%iounit_result = iounit_result + i
    this%iounit_runave = iounit_runave + i
    this%iounit_errors = iounit_errors + i
    this%iounit_visual = iounit_visual + i
#if  TRANS == 1
    this%iounit_rescf  = iounit_rescf  + i   !TRANSPORT_thisline
#endif

    write( IOBuffer, '(T15, "Reading ensemble ", I3, " successful")') &
&          this%EnsembleNumber
    call LogWrite
    write( IOBuffer, '(72(1H-))')
    call LogWrite

  end subroutine TEnsemble_Construct



!==============================================================!
!  Subroutine TEnsemble_ConstructSVC                           !
!==============================================================!

  subroutine TEnsemble_ConstructSVC( this, ne )

    implicit none

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: ne

    ! Declare local variables
    integer                   :: i, j
    integer                   :: stat
    type(TComponent), pointer :: pc
    character(FileNameLength) :: PotModFileName
    real(RK)                  :: scaleSigma, scaleEpsilon


#if MPI_VER > 0

    call Error('Up to now, SVC can only be used with the serial version.' )

#endif

    ! Allocate simulation box length
    allocate( this%BoxLength, STAT = stat )
    call AllocationError( stat, 'simulation box length' )

    ! Allocate maximum number of particles
    allocate( this%NPartMax, STAT = stat )
    call AllocationError( stat, 'maximum number of particles' )

    ! Set number of ensemble
    this%EnsembleNumber = ne
    call LogWriteBlank
    write( IOBuffer, '(72(1H-))')
    call LogWrite
    write( IOBuffer, '(T14, "Reading parameters of ensemble", I3)' ) &
&     this%EnsembleNumber
    call LogWrite

    ! Read temperature
    call FileReadParameter( this%Temperature, iounit_params , IdRefTemperature, .false. )
    if( .not. UseReducedUnits ) then
      this%Temperature = this%Temperature / UnitTemperature
    end if

    ! Update log file
    write( IOBuffer, '("Temperature: ",T26, F9.3, " K")' ) &
&     this%RefTemperature * UnitTemperature
    call LogWrite
    write( IOBuffer, '("Reduced temperature: ",T26, F12.6)' ) this%RefTemperature
    call LogWrite

    ! Read number of components in ensemble
    call FileReadParameter( this%NComponents, iounit_params , IdNComponents, .false. )
    write( IOBuffer, '("Number of components:",T28, I3)' ) this%NComponents
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
      call FileReadParameter( PotModFileName, iounit_params , IdPotModFileName, .false. )

      call Construct( this%Component(i), PotModFileName )
      call Construct( this%Component(i+1), PotModFileName )

    end do

    ! Calculate number of particles in each process
    do i = 1, this%NComponents
      pc => this%Component(i)
      pc%NPart = NOrient
      pc%NPart1 = ProcRange( pc%NPart, pc%NPart0, pc%NPart2 )
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
        call FileReadParameter( scaleSigma, iounit_params , IdScaleSigma, .false. )
        this%ScaleSigma(i:i+1, j:j+1) = scaleSigma
        if( i /= j ) this%ScaleSigma(j:j+1, i:i+1) = scaleSigma
        call FileReadParameter( scaleEpsilon, iounit_params , IdScaleEpsilon, .false. )
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

  end subroutine TEnsemble_Destruct



!==============================================================!
!  Subroutine TEnsemble_CreateComponents                       !
!==============================================================!

  subroutine TEnsemble_CreateComponents( this )

    implicit none

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

    this%NGradInsComp = 0
    ! Create components for fluctuating particle states
    this%NRealComponents = ncomp
    do i = 1, this%NRealComponents
      if( this%Component(i)%ChemPotMethod .eq. ChemPotMethodGradIns ) then
        this%Component(i)%NGradThis = this%NGradInsComp
        this%NGradInsComp = this%NGradInsComp + 1
        nfluct = this%Component(i)%Molecule%NFluct
        ncomp = ncomp + nfluct

        ! Reallocate component array
        allocate( reallocate(ncomp), STAT = stat )
        ! reallocate will be stored in this%Component, which thus has to be deallocated
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

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Destroy components
    do i = 1, this%NComponents
      if ( i .le. this%NRealComponents ) then
        call Destruct( this%Component(i) )
      else
        call DestructFluct( this%Component(i) )
      end if
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
          this%Interaction(i,j)%OptPressure = this%OptPressure
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
          this%Interaction(i,j)%OptPressure = this%OptPressure
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
#if  TRANS == 1
!TRANSPORT_start
    ! 4.) Transport properties
    if( CorrfunMode .eq. active ) then
      do i = 1, this%NComponents
        call ConstructCF( this%Sumself_i(i),  .true. )
        this%Sumself_i(i)%BLOCKSUM(:) = 0._RK
        this%Sumself_i(i)%TOTALSUM    = 0._RK
        this%Sumself_i(i)%AVERAGE     = 0._RK
        this%Sumself_i(i)%VARIANCE    = 0._RK
      end do

      call ConstructCF( this%SumBin_d,   .true. )
      this%SumBin_d%BLOCKSUM(:) = 0._RK
      this%SumBin_d%TOTALSUM    = 0._RK
      this%SumBin_d%AVERAGE     = 0._RK
      this%SumBin_d%VARIANCE    = 0._RK

      call ConstructCF( this%SumTer_a,  .true. )
      this%SumTer_a%BLOCKSUM(:) = 0._RK
      this%SumTer_a%TOTALSUM    = 0._RK
      this%SumTer_a%AVERAGE     = 0._RK
      this%SumTer_a%VARIANCE    = 0._RK
      call ConstructCF( this%SumTer_b,  .true. )
      this%SumTer_b%BLOCKSUM(:) = 0._RK
      this%SumTer_b%TOTALSUM    = 0._RK
      this%SumTer_b%AVERAGE     = 0._RK
      this%SumTer_b%VARIANCE    = 0._RK
      call ConstructCF( this%SumTer_c,  .true. )
      this%SumTer_c%BLOCKSUM(:) = 0._RK
      this%SumTer_c%TOTALSUM    = 0._RK
      this%SumTer_c%AVERAGE     = 0._RK
      this%SumTer_c%VARIANCE    = 0._RK


      call ConstructCF( this%SumVisco_s, .true. )
      this%SumVisco_s%BLOCKSUM(:) = 0._RK
      this%SumVisco_s%TOTALSUM    = 0._RK
      this%SumVisco_s%AVERAGE     = 0._RK
      this%SumVisco_s%VARIANCE    = 0._RK

      call ConstructCF( this%SumVisco_b, .true. )
      this%SumVisco_b%BLOCKSUM(:) = 0._RK
      this%SumVisco_b%TOTALSUM    = 0._RK
      this%SumVisco_b%AVERAGE     = 0._RK
      this%SumVisco_b%VARIANCE    = 0._RK

      call ConstructCF( this%SumConduct, .true. )
      this%SumConduct%BLOCKSUM(:) = 0._RK
      this%SumConduct%TOTALSUM    = 0._RK
      this%SumConduct%AVERAGE     = 0._RK
      this%SumConduct%VARIANCE    = 0._RK
    end if
!TRANSPORT_END
#endif

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
#if  TRANS == 1
!TRANSPORT_start
    if( CorrfunMode.eq.active ) then
      do i = 1, this%NComponents
         call DestructCF( this%Sumself_i(i) )
      end do

      call DestructCF( this%SumBin_d   )
      call DestructCF( this%SumTer_a   )
      call DestructCF( this%Sumter_b   )
      call DestructCF( this%Sumter_c   )
      call DestructCF( this%SumVisco_s )
      call DestructCF( this%SumVisco_b )
      call DestructCF( this%SumConduct )
    end if
!TRANSPORT_END
#endif
    do i = 1, this%NRealComponents
      call DestroyAccumulators( this%Component(i) )
    end do


  end subroutine TEnsemble_DestroyAccumulators



!==============================================================!
!  Subroutine TEnsemble_CalculateNPart                         !
!==============================================================!

  subroutine TEnsemble_CalculateNPart( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK)                  :: s
    type(TComponent), pointer :: pc
    integer                   :: i

    ! Adjust number of cells to cube of integer
    if( this%NPart < NPartInCell ) this%NPart = NPartInCell

    ! Set maximum number of particles
    if( EnsembleType .eq. EnsembleTypeGE .or. &
&       EnsembleType .eq. EnsembleTypeHA ) then
      this%NPartMax = 2 * this%NPart
    else
      this%NPartMax = this%NPart
    end if

    ! Normalize mole fractions
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
    ! according to mole fraction
    this%Component(this%NRealComponents)%NPart = this%NPart
    do i = 1, this%NRealComponents - 1
      pc => this%Component(i)
      pc%NPart = nint( this%NPart * pc%Fraction )
      this%Component(this%NRealComponents)%NPart = &
&       this%Component(this%NRealComponents)%NPart - pc%NPart
    end do

    ! Set mole fractions according to real number of particles
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
      pc%NPart1 = ProcRange( pc%NPart, pc%NPart0, pc%NPart2 )
      if( pc%NTest > 0 ) pc%NTest = 1 + (pc%NTest - 1) / NProcs
      pc%NTestAll = NProcs * pc%NTest
      this%NTestMax = max( pc%NTest, this%NTestMax )
      this%NFluctMax = max( pc%NFluctMax, this%NFluctMax )
    end do

    ! Calculate volume of simulation box
    this%Volume0 = this%NPart / this%RefDensity

    ! Update log file
    write( IOBuffer, '("Number of particles:",T36, I11)' ) this%NPart
    call LogWrite
    do i = 1, this%NRealComponents
      pc => this%Component(i)
      write( IOBuffer, '("Mole fraction of ", A, ":",T45, F6.3)' ) &
&       trim( pc%PotModFileName ), pc%Fraction
      call LogWrite
      write( IOBuffer, '(T10,"->  number of particles: ",T36, I11)' ) &
&       pc%NPart
      call LogWrite
      if( pc%NTestAll > 0 ) then
        write( IOBuffer, '("Number of test particles:",T36, I11)' ) pc%NTestAll
        call LogWrite
      end if
    end do

  end subroutine TEnsemble_CalculateNPart



!==============================================================!
!  Subroutine TEnsemble_UpdateBoxLength                        !
!==============================================================!

  subroutine TEnsemble_UpdateBoxLength( this )

    implicit none

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
    end if

    do i = 1, this%NComponents
      do j = 1, this%NComponents
        call UpdateBoxLength( this%Interaction( j, i ), this%BoxLength )
      end do
    end do

  end subroutine TEnsemble_UpdateBoxLength



!==============================================================!
!  Subroutine TEnsemble_UpdateFractions                        !
!==============================================================!

  subroutine TEnsemble_UpdateFractions( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer                   :: i
    type(TComponent), pointer :: pc

    ! Set mole fractions according to real number of particles
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

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i
    integer :: stat
#if TRANS ==1
    integer :: Npart3
    integer :: Ncomp2
#endif

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
      call LogWriteBlank
      write( IOBuffer, '("Memory for test particles allocated successfully")' )
      call LogWrite
    end if

    ! Allocate components
    do i = 1, this%NComponents
      this%Component(i)%NPartMax => this%NPartMax
      if( this%Component(i)%NTest > 0 ) then
        this%Component(i)%P0Test => this%P0Test
        this%Component(i)%Q0Test => this%Q0Test
      end if
      this%Component(i)%BoxLength => this%BoxLength
      call Allocate( this%Component(i) )
    end do

#if  TRANS == 1
!TRANSPORT_start
      Npart3 = 3*this%Npart
      Ncomp2 = this%NComponents*this%NComponents

    ! Allocate correlation fucntions
     if( CorrfunMode .eq. active ) then

      allocate( this%cf_vs(this%Ncorr), STAT = stat )
      call AllocationError( stat, 'viscosity_shear_cf_vs', this%Ncorr )

      allocate( this%cf_vb(this%Ncorr), STAT = stat )
      call AllocationError( stat, 'viscosity_bulk_cf_vb', this%Ncorr )

      allocate( this%cf_c(this%Ncorr), STAT = stat )
      call AllocationError( stat, 'conductivity_cf_c', this%Ncorr )

      allocate( this%cf_d( this%NComponents, this%Ncorr), STAT = stat )
      call AllocationError( stat, 'self_diffusion', this%Ncorr )

      allocate( this%cf_db( this%Ncorr), STAT = stat )
      call AllocationError( stat, 'binary_diffusion', this%Ncorr )

      allocate( this%lamda( NComp2, this%Ncorr ), STAT = stat )
      call AllocationError( stat, 'binary_diffusion', this%Ncorr )

      allocate( this%sinte_i( this%NComponents, this%Ncorr), STAT = stat )
      call AllocationError( stat, 'self_diffusion_integrated', this%Ncorr )

      allocate( this%sinte_db( this%Ncorr), STAT = stat )
      call AllocationError( stat, 'mutual diffusion integrated', this%Ncorr )

      allocate( this%sinte_lamda( Ncomp2, this%Ncorr), STAT = stat )
      call AllocationError( stat, 'mutual diffusion integrated', this%Ncorr )

      allocate( this%sinte_vs( this%Ncorr), STAT = stat )
      call AllocationError( stat, 'shear_viscosity_integrated', this%Ncorr )

      allocate( this%sinte_vb( this%Ncorr), STAT = stat )
      call AllocationError( stat, 'bulk_viscosity_integrated', this%Ncorr )

      allocate( this%sinte_c( this%Ncorr), STAT = stat )
      call AllocationError( stat, 'Thermal_conductivity_integrated', this%Ncorr )

      allocate( this%a( NPart3, this%Ncorr), STAT = stat  )
      call AllocationError( stat, 'diffusion_matrix', NPart3 )

      allocate( this%vsk(this%Ncorr, 3), STAT = stat  )
      call AllocationError( stat, 'vsk', this%Npart )

      allocate( this%vsp(this%Ncorr, 3), STAT = stat  )
      call AllocationError( stat, 'vsp', this%Npart )

      allocate( this%vbk(this%Ncorr, 3), STAT = stat  )
      call AllocationError( stat, 'vbk', this%Npart )

      allocate( this%vbp(this%Ncorr, 3), STAT = stat  )
      call AllocationError( stat, 'vbp', this%Npart )

      allocate( this%vckt(this%Ncorr, 3), STAT = stat  )
      call AllocationError( stat, 'vckt', this%Npart )

      allocate( this%vckr(this%Ncorr, 3), STAT = stat  )
      call AllocationError( stat, 'vckr', this%Npart )

      allocate( this%vcpt(this%Ncorr, 3), STAT = stat  )
      call AllocationError( stat, 'vcpt', this%Npart )

      allocate( this%vcpr(this%Ncorr, 3), STAT = stat  )
      call AllocationError( stat, 'vcpr', this%Npart )

      allocate( this%selfd_i(this%NComponents), STAT = stat  )
      call AllocationError( stat, 'selfd_i', this%NComponents )

      allocate( this%Sumself_i(this%NComponents), STAT = stat  )
      call AllocationError( stat, 'Sumselfd_i', this%NComponents )


      ! Set correlation-fucntion vectors
      this%cf_d(: , :)    = 0._RK
      this%cf_db( :  )    = 0._RK
      this%cf_vs( :  )    = 0._RK
      this%cf_vb( :  )    = 0._RK
      this%cf_c( :   )    = 0._RK

      this%a(:  ,   :   ) = 0._RK
      this%lamda( : , : ) = 0._RK
      this%sinte_i(: , :) = 0._RK
      this%sinte_db(:   ) = 0._RK
      this%sinte_lamda(:,:)=0._RK
      this%sinte_vs(:   ) = 0._RK
      this%sinte_vb(:   ) = 0._RK
      this%sinte_c(:    ) = 0._RK

      this%selfd_i(:)     = 0._RK
      this%vsk(: ,  :)    = 0._RK
      this%vsp(: ,  :)    = 0._RK
      this%vbk(: ,  :)    = 0._RK
      this%vbp(: ,  :)    = 0._RK
      this%vckt(: ,  :)   = 0._RK
      this%vckr(: ,  :)   = 0._RK
      this%vcpt(: ,  :)   = 0._RK
      this%vcpr(: ,  :)   = 0._RK

      this%sc(:) = 0._RK
      this%sp(:) = 0._RK
    end if
    !TRANSPORT_END
#endif

  end subroutine TEnsemble_Allocate



!==============================================================!
!  Subroutine TEnsemble_DeallocateEPot                         !
!==============================================================!

  subroutine TEnsemble_DeallocateEPot( this )

    implicit none

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
#if  TRANS == 1
!TRANSPORT_start
    ! Deallocate arrays for correlation fucntions
    if( associated( this%cf_d  ) )   then
      deallocate( this%cf_d  )
    end if
    if( associated( this%cf_db ) )   then
      deallocate( this%cf_db )
    end if
    if( associated( this%cf_vs ) )   then
      deallocate( this%cf_vs )
    end if
    if( associated( this%cf_vb ) )   then
      deallocate( this%cf_vb )
    end if
    if( associated( this%cf_c  ) )   then
      deallocate( this%cf_c  )
    end if
    if( associated( this%lamda ) )   then
      deallocate( this%lamda  )
    end if
    if( associated( this%sinte_i) )  then
      deallocate( this%sinte_i  )
    end if
    if( associated( this%sinte_db) ) then
      deallocate( this%sinte_db )
    end if
    if( associated( this%sinte_vs) ) then
      deallocate( this%sinte_vs )
    end if
    if( associated( this%sinte_vb) ) then
      deallocate( this%sinte_vb )
    end if
    if( associated( this%sinte_c)  ) then
      deallocate( this%sinte_c  )
    end if

    if( associated( this%a )  )   then
      deallocate( this%a  )
    end if
    if( associated( this%vsk ) )  then
      deallocate( this%vsk )
    end if
    if( associated( this%vsp ) )  then
      deallocate( this%vsp )
    end if
    if( associated( this%vbk ) )  then
      deallocate( this%vbk )
    end if
    if( associated( this%vbp ) )  then
      deallocate( this%vbp )
    end if
    if( associated( this%vckt ) ) then
      deallocate( this%vckt )
    end if
    if( associated( this%vckr ) ) then
      deallocate( this%vckr )
    end if
    if( associated( this%vcpt ) ) then
      deallocate( this%vcpt )
    end if
    if( associated( this%vcpr ) ) then
      deallocate( this%vcpr )
    end if
    if( associated( this%selfd_i ) ) then
      deallocate( this%selfd_i )
    end if
    if( associated( this%Sumself_i ) ) then
      deallocate( this%Sumself_i )
    end if
!TRANSPORT_END
#endif
  end subroutine TEnsemble_Deallocate



!==============================================================!
!  Subroutine TEnsemble_CalculateCorr                          !
!==============================================================!

  subroutine TEnsemble_CalculateCorr( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK)                        :: NPartInv, Scale, RFConst
    type(TComponent), pointer       :: pc
    type(TPotLJ126LJ126), pointer   :: plj
    integer                         :: i1, i2, j1, j2

    ! Assign local variables
    NPartInv = 1._RK / this%NPart
    RFConst = -1._RK / this%RCutoffDipoleDipole**3 &
&     * (this%RFEpsilon - 1._RK) / (2._RK * this%RFEpsilon + 1._RK)

    ! Set maximum cutoff radius
    this%NRCutoffMax = 0
    this%RCutoffMax2 = 0._RK

    ! Zero long-range corrections
    this%EPotCorrLJ   = 0._RK
    this%EPotCorrRF   = 0._RK
    this%VirialCorrLJ = 0._RK
    do i1 = 1, this%NComponents
      this%Component(i1)%EPotTestCorrLJ = 0._RK
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
    if( (this%NChargeMax > 0).or.(this%NDipoleMax > 0) ) then
      do i1 = 1, this%NComponents
        pc => this%Component(i1)
        pc%EPotTestCorrRF = pc%Molecule%MueSquared * 2._RK * RFConst
        this%EPotCorrRF = this%EPotCorrRF + pc%Molecule%MueSquared * pc%NPart
      end do
      this%EPotCorrRF = this%EPotCorrRF * RFConst / NProcs
    end if

#if MPI_VER >0
  if ( SimulationType .eq. MonteCarlo .and. (.not.(Equilibration .and. CommonEqui)) ) then 
    if( this%NLJ126Max > 0 ) then    
      this%EPotCorrLJ = this%EPotCorrLJ * NProcs
      this%VirialCorrLJ = this%VirialCorrLJ * NProcs
    endif
    if( (this%NChargeMax > 0).or.(this%NDipoleMax > 0) ) then
      this%EPotCorrRF = this%EPotCorrRF * NProcs
    endif
  endif
#endif

  end subroutine TEnsemble_CalculateCorr



!==============================================================!
!  Subroutine TEnsemble_InitPositions                          !
!==============================================================!

  subroutine TEnsemble_InitPositions( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer                   :: i, j, k, l, n, nc, nm
    real(RK), dimension(3)    :: xl
    integer                   :: comp(this%NComponents)
    type(TComponent), pointer :: pc
    real(RK)                  :: NCells                    ! Number of unit cells
    integer, dimension(3)     :: NCells1dim                ! Number of unit cells in one dimension of lattice


    NCells = ceiling( real( this%NPart, RK ) / real( NPartInCell, RK ) )
    NCells1dim = ceiling( NCells**Third )

    if( (NCells1dim(1)-1)*NCells1dim(2)*NCells1dim(3)>=NCells ) NCells1dim(1)=NCells1dim(1)-1
    if( NCells1dim(1)*(NCells1dim(2)-1)*NCells1dim(3)>=NCells ) NCells1dim(2)=NCells1dim(2)-1

    ! Initialize comp array
    do i = 1, this%NComponents
      comp(i) = this%Component(i)%NPart
    end do

    ! Create FCC lattice
    n = 0
    xl(1) = 1._RK / real( NCells1dim(1), RK )
    xl(2) = 1._RK / real( NCells1dim(2), RK )
    xl(3) = 1._RK / real( NCells1dim(3), RK )
    call LogWriteBlank
    write( IOBuffer, '("Initialize positions:")' )
    call LogWrite
    write( IOBuffer, '(T10,"FCC lattice: ",I3," *",I4,"x",I4,"x",I4," cells")' ) &
&          NPartInCell, ( NCells1dim(i), i=1,3 )
    call LogWrite
    write( IOBuffer, '(T10, "with",I3," molecules/cell")' ) NPartInCell
    call LogWrite
loop:do l = 1, NPartInCell
      do i = 1, NCells1dim(1)
        do j = 1, NCells1dim(2)
          do k = 1, NCells1dim(3)
            nc = select_component( comp )
            nm = comp(nc) + 1
            pc => this%Component(nc)
            pc%P0(nm, 1) = xl(1) * (CellX(l) + i - 1)
            pc%P0(nm, 2) = xl(2) * (CellY(l) + j - 1)
            pc%P0(nm, 3) = xl(3) * (CellZ(l) + k - 1)
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
    call MPI_Bcast( this%Temperature, 1, MPI_RK, &
&     NRootProc, Communicator, ierror )
#endif

  end subroutine TEnsemble_CalculateEKin



!==============================================================!
!  Subroutine TEnsemble_CheckNPart                             !
!==============================================================!

  subroutine TEnsemble_CheckNPart( this, NPartOk )

    implicit none

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
    call Force( this )
    call ChemicalPotential( this )
    call Atom2Mol( this )
    call Correct( this )
#if  TRANS == 1
!TRANSPORT_start
    if( .not.Equilibration.and.(CorrfunMode .eq. active) )then
      call CalCorrFun( this )
         if((this%Mmess.gt.0).and.(mod(Step, this%NSpancf).eq.0))then
           call IntCorrFun( this )
         end if
    end if
!TRANSPORT_END
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

    ! Zero number of MC attempts and successes
    if( Step == 1 ) call ZeroNAttempts( this )

    ! Outer loop
    do i = 1, this%NDF / 3

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
    ! in MC simulations we only communicate during common equilibration 
    if (Equilibration .and. CommonEqui) then
      call MPI_Allreduce( GetEnergy( this ), this%EPot, 1 , &
&       MPI_RK, MPI_SUM, Communicator, ierror )
        if ( this%OptPressure ) then
          call MPI_Allreduce( GetVirial( this ), this%Virial, 1 , &
&          MPI_RK, MPI_SUM, Communicator, ierror )
        endif
    else
      this%EPot = GetEnergy( this )
      if ( this%OptPressure ) then
        this%Virial = GetVirial( this )
      endif
    endif  

#else
    this%EPot = GetEnergy( this )
    if ( this%OptPressure ) then
      this%Virial = GetVirial( this )
    endif
#endif
    ! Resize simulation box
    if( ConstantPressure .and. .not. NVTEquilibration ) then
      call Resize( this )

      ! Check whether cutoff radius is too large
      if( this%RCutoffMax2 > this%BoxLength ) &
&       this%NRCutoffMax = this%NRCutoffMax + 1
    end if

    if ( this%OptPressure ) then
      ! Calculate pressure
      this%Pressure = this%Density * this%Temperature &
&                     + this%Virial / this%Volume0
    end if

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

      ! Initialize result
      integral = 0._RK

      ! Return if no values to integrate
      if( n < 1 ) return
!       if( n < 3 ) then
!         print *,"ERROR: TEnsemble_RunSVCStep simpson: n=",n,"<3" ! DEBUG
!         return
        ! could automatically use
        ! - a trapazoidal rule for n=2
        ! - a quadrilateral rule for n=1
!       end if

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
      call MPI_Bcast( this%Volume0, 1, MPI_RK, NRootProc, &
&       Communicator, ierror )
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
      end if
#if MPI_VER > 0
      call MPI_Bcast( this%Volume0, 1, MPI_RK, NRootProc, &
&       Communicator, ierror )
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
      call MPI_Bcast( this%Volume0, 1, MPI_RK, NRootProc, &
&       Communicator, ierror )
#endif
      call UpdateBoxLength( this )
    end if

  end subroutine TEnsemble_PredictLeapFrog



!==============================================================!
!  Subroutine TEnsemble_CorrectLeapFrog                        !
!==============================================================!

  subroutine TEnsemble_CorrectLeapFrog( this )

    implicit none

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

    ! Zero forces
    do i = 1, this%NComponents
      pc => this%Component(i)
      do j = 1, this%Component(i)%Molecule%NLJ126
        pc%Molecule%SiteLJ126(j)%FX(1:pc%NPart) = 0._RK
        pc%Molecule%SiteLJ126(j)%FY(1:pc%NPart) = 0._RK
        pc%Molecule%SiteLJ126(j)%FZ(1:pc%NPart) = 0._RK
#if  TRANS == 1
        !TRANSPORT_start
        pc%Molecule%SiteLJ126(j)%vsLJx(1:pc%Npart) = 0._RK
        pc%Molecule%SiteLJ126(j)%vsLJy(1:pc%Npart) = 0._RK
        pc%Molecule%SiteLJ126(j)%vsLJz(1:pc%Npart) = 0._RK
        pc%Molecule%SiteLJ126(j)%vsuLJx(1:pc%Npart)= 0._RK
        pc%Molecule%SiteLJ126(j)%vsuLJy(1:pc%Npart)= 0._RK
        pc%Molecule%SiteLJ126(j)%vsuLJz(1:pc%Npart)= 0._RK
        pc%Molecule%SiteLJ126(j)%vbLJx(1:pc%Npart) = 0._RK
        pc%Molecule%SiteLJ126(j)%vbLJy(1:pc%Npart) = 0._RK
        pc%Molecule%SiteLJ126(j)%vbLJz(1:pc%Npart) = 0._RK
        pc%Molecule%SiteLJ126(j)%cLJx(1:pc%Npart)  = 0._RK
        pc%Molecule%SiteLJ126(j)%cLJy(1:pc%Npart)  = 0._RK
        pc%Molecule%SiteLJ126(j)%cLJz(1:pc%Npart)  = 0._RK
        pc%Molecule%SiteLJ126(j)%tuLJx(1:pc%Npart) = 0._RK
        pc%Molecule%SiteLJ126(j)%tuLJy(1:pc%Npart) = 0._RK
        pc%Molecule%SiteLJ126(j)%tuLJz(1:pc%Npart) = 0._RK
        pc%Molecule%SiteLJ126(j)%tlLJx(1:pc%Npart) = 0._RK
        pc%Molecule%SiteLJ126(j)%tlLJy(1:pc%Npart) = 0._RK
        pc%Molecule%SiteLJ126(j)%tlLJz(1:pc%Npart) = 0._RK
        pc%Molecule%SiteLJ126(j)%tdLJx(1:pc%Npart) = 0._RK
        pc%Molecule%SiteLJ126(j)%tdLJy(1:pc%Npart) = 0._RK
        pc%Molecule%SiteLJ126(j)%tdLJz(1:pc%Npart) = 0._RK
        !TRANSPORT_END
#endif
      end do
      do j = 1, this%Component(i)%Molecule%NCharge
        pc%Molecule%SiteCharge(j)%FX(1:pc%NPart) = 0._RK
        pc%Molecule%SiteCharge(j)%FY(1:pc%NPart) = 0._RK
        pc%Molecule%SiteCharge(j)%FZ(1:pc%NPart) = 0._RK
#if  TRANS == 1
        !TRANSPORT_start
        pc%Molecule%SiteCharge(j)%vsCx(1:pc%Npart) = 0._RK
        pc%Molecule%SiteCharge(j)%vsCy(1:pc%Npart) = 0._RK
        pc%Molecule%SiteCharge(j)%vsCz(1:pc%Npart) = 0._RK
        pc%Molecule%SiteCharge(j)%vsuCx(1:pc%Npart)= 0._RK
        pc%Molecule%SiteCharge(j)%vsuCy(1:pc%Npart)= 0._RK
        pc%Molecule%SiteCharge(j)%vsuCz(1:pc%Npart)= 0._RK
        pc%Molecule%SiteCharge(j)%vbCx(1:pc%Npart) = 0._RK
        pc%Molecule%SiteCharge(j)%vbCy(1:pc%Npart) = 0._RK
        pc%Molecule%SiteCharge(j)%vbCz(1:pc%Npart) = 0._RK
        pc%Molecule%SiteCharge(j)%cCx(1:pc%Npart)  = 0._RK
        pc%Molecule%SiteCharge(j)%cCy(1:pc%Npart)  = 0._RK
        pc%Molecule%SiteCharge(j)%cCz(1:pc%Npart)  = 0._RK
        pc%Molecule%SiteCharge(j)%tuCx(1:pc%Npart) = 0._RK
        pc%Molecule%SiteCharge(j)%tuCy(1:pc%Npart) = 0._RK
        pc%Molecule%SiteCharge(j)%tuCz(1:pc%Npart) = 0._RK
        pc%Molecule%SiteCharge(j)%tlCx(1:pc%Npart) = 0._RK
        pc%Molecule%SiteCharge(j)%tlCy(1:pc%Npart) = 0._RK
        pc%Molecule%SiteCharge(j)%tlCz(1:pc%Npart) = 0._RK
        pc%Molecule%SiteCharge(j)%tdCx(1:pc%Npart) = 0._RK
        pc%Molecule%SiteCharge(j)%tdCy(1:pc%Npart) = 0._RK
        pc%Molecule%SiteCharge(j)%tdCz(1:pc%Npart) = 0._RK
        !TRANSPORT_END
#endif
      end do
      do j = 1, this%Component(i)%Molecule%NDipole
        pc%Molecule%SiteDipole(j)%FX(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%FY(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%FZ(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%TX(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%TY(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%TZ(1:pc%NPart) = 0._RK
#if  TRANS == 1
        !TRANSPORT_start
        pc%Molecule%SiteDipole(j)%vsDx(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%vsDy(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%vsDz(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%vsuDx(1:pc%NPart)= 0._RK
        pc%Molecule%SiteDipole(j)%vsuDy(1:pc%NPart)= 0._RK
        pc%Molecule%SiteDipole(j)%vsuDz(1:pc%NPart)= 0._RK
        pc%Molecule%SiteDipole(j)%vbDx(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%vbDy(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%vbDz(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%cDx(1:pc%NPart)  = 0._RK
        pc%Molecule%SiteDipole(j)%cDy(1:pc%NPart)  = 0._RK
        pc%Molecule%SiteDipole(j)%cDz(1:pc%NPart)  = 0._RK
        pc%Molecule%SiteDipole(j)%tuDx(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%tuDy(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%tuDz(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%tlDx(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%tlDy(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%tlDz(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%tdDx(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%tdDy(1:pc%NPart) = 0._RK
        pc%Molecule%SiteDipole(j)%tdDz(1:pc%NPart) = 0._RK
        !TRANSPORT_END
#endif
      end do
      do j = 1, this%Component(i)%Molecule%NQuadrupole
        pc%Molecule%SiteQuadrupole(j)%FX(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%FY(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%FZ(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%TX(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%TY(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%TZ(1:pc%NPart) = 0._RK
#if  TRANS == 1
        !TRANSPORT_start
        pc%Molecule%SiteQuadrupole(j)%vsQx(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%vsQy(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%vsQz(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%vsuQx(1:pc%NPart)= 0._RK
        pc%Molecule%SiteQuadrupole(j)%vsuQy(1:pc%NPart)= 0._RK
        pc%Molecule%SiteQuadrupole(j)%vsuQz(1:pc%NPart)= 0._RK
        pc%Molecule%SiteQuadrupole(j)%vbQx(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%vbQy(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%vbQz(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%cQx(1:pc%NPart)  = 0._RK
        pc%Molecule%SiteQuadrupole(j)%cQy(1:pc%NPart)  = 0._RK
        pc%Molecule%SiteQuadrupole(j)%cQz(1:pc%NPart)  = 0._RK
        pc%Molecule%SiteQuadrupole(j)%tuQx(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%tuQy(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%tuQz(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%tlQx(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%tlQy(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%tlQz(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%tdQx(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%tdQy(1:pc%NPart) = 0._RK
        pc%Molecule%SiteQuadrupole(j)%tdQz(1:pc%NPart) = 0._RK
        !TRANSPORT_END
#endif
      end do
      if( pc%Molecule%isElongated ) then
        pc%tRFX(:) = 0._RK
        pc%tRFY(:) = 0._RK
        pc%tRFZ(:) = 0._RK
      end if
#if  TRANS == 1
      !TRANSPORT_start
      do j = 1, this%Component(i)%Npart
        this%Component(i)%fs(j, 1)    = 0._RK
        this%Component(i)%fs(j, 2)    = 0._RK
        this%Component(i)%fs(j, 3)    = 0._RK
        this%Component(i)%fb(j, 1)    = 0._RK
        this%Component(i)%fb(j, 2)    = 0._RK
        this%Component(i)%fb(j, 3)    = 0._RK
        this%Component(i)%ftc(j, 1)   = 0._RK
        this%Component(i)%ftc(j, 2)   = 0._RK
        this%Component(i)%ftc(j, 3)   = 0._RK
        this%Component(i)%frc(j, 1)   = 0._RK
        this%Component(i)%frc(j, 2)   = 0._RK
        this%Component(i)%frc(j, 3)   = 0._RK
      end do
      !TRANSPORT_END
#endif
    end do

    ! Zero potential
    EPot = this%Density * this%EPotCorrLJ + this%EPotCorrRF

    ! Zero virial
    Virial = this%Density * this%VirialCorrLJ

    ! Loop over components
    do i = 1, this%NComponents
      do j = i, this%NComponents
        call Force( this%Interaction( i, j ), &
&                   EPot, Virial, this%BoxLength )
      end do
    end do

    ! Collect sums from all processes
#if MPI_VER > 0
    call MPI_Reduce( EPot, this%EPot, 1, MPI_RK, MPI_SUM, &
&     NRootProc, Communicator, ierror )
    call MPI_Reduce( Virial, this%Virial, 1, MPI_RK, MPI_SUM, &
&     NRootProc, Communicator, ierror )
#else
    this%EPot = EPot
    this%Virial = Virial
#endif

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
    integer                   :: nstate( 0:this%NFluctMax )
#if MPI_VER > 0
    integer                   :: tempComm
    integer                   :: tempVec(0:this%NFluctMax)
    integer                   :: tempVal, tempVal2
    integer                   :: tempVec1(this%NFluctMax), tempVec2(this%NFluctMax)
    integer                   :: tempVec3(this%NFluctMax), tempVec4(this%NFluctMax)
    
    !integer                   :: color
#endif
    ! No calculation of chemical potential in equilibration
    if( Equilibration ) then
      do i = 1, this%NRealComponents
        this%Component(i)%CalcChemPot = .false.
        this%Component(i)%ChemPot = 0._RK
        this%Component(i)%ChemPot1 = 0._RK
        this%Component(i)%ChemPot2 = 0._RK
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

#if MPI_VER > 0
    tempComm = Communicator
#endif

    ! Outer loop over components
componentLoop:    do i = 1, this%NRealComponents

      pc => this%Component(i)
      if( Equilibration .and. pc%WFMethod .ne. WFMethodGuess ) cycle
      select case( pc%ChemPotMethod )

      ! Chemical potential by gradual insertion
      case( ChemPotMethodGradIns )
        if( Equilibration) cycle componentLoop
        if( GradInsInitialization .and. (pc%WFMethod .ne. WFMethodGuess)) cycle componentLoop

        ! Reset variables
        if( Step == 1 ) then
          pc%ProbW0 = 0._RK
          pc%ProbW1 = 0._RK
          pc%ProbW0V = 0._RK
          pc%ProbW1Rho = 0._RK
          pc%NStateWF(:) = 0
!           pc%NStateBF(:) = 0
        end if

        ! determine, if chemical potential has to be calculated
        pc%CalcChemPot = .false.
        if( GradInsFrequency > 0 ) then
          if( mod( Step, GradInsFrequency ) == 0 ) pc%CalcChemPot = .true.
        end if
        if (GradInsInitialization) then
          pc%CalcChemPot = .true.
        end if

#if MPI_VER > 0
        ! Per Process we calculate GI only for one component 
       ! if (mod(NProc,this%NRealComponents)/=i-1) then
      !   color = 100000*i
     !   else
     !    color = 2000*i

      !  endif
        if (mod(NProc,this%NGradInsComp)/= pc%NGradThis) pc%CalcChemPot = .false.

      !  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,NProc,Communicator,ierror) 
           ! Careful, Nproc and NProcs are now specific for Communicator
      !  call SetCommunicator( Communicator )
#endif

        if( pc%CalcChemPot ) then
          !pc%CalcChemPot = .true.

          ! Save current state
!           call SaveState( this )

          ndfmove = this%NDF
          ndfbiased = this%NDF * 50
          ndffluct = this%NDF * 10
          ndfchange = this%NDF * 10
          ndfcp = ndfmove + ndfbiased + ndffluct + ndfchange
          pc%NState(:) = 1
          nstate = 0

          ! Set fluctuating particle
          ncf = pc%NFluctComp( pc%NFluctState )
          if( ncf == i ) then
            npf = rnd( pc%NPart )
            call Move2End( this, ncf, npf )
          else
            npf = 1
          end if

giloop:   do j = 1, NFullFluct * ndfcp

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

              ! Acceleration of MC Moves
              if (np .gt. this%Component(nc)%NPart) cycle

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

              ! Acceleration of MC Moves
              if (np .gt. this%Component(nc)%NPart) cycle

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
              nstate(pc%NFluctState) = nstate(pc%NFluctState) + 1

            end if

          end do giloop

          pc%NStateWF = pc%NStateWF + nstate(0:pc%NFluctMax)
          pc%NState = pc%NState + nstate(0:pc%NFluctMax)


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
          pc%ChemPot1 = pc%ProbW0 / pc%ProbW1
          pc%ChemPot2 = pc%ProbW0V / pc%ProbW1

        else
          !pc%CalcChemPot = .false.
          pc%ChemPot = 0._RK
        end if

        if( mod( Step, ErrorsUpdateFrequency ) == 0 .or. &
&           ( GradInsInitialization .and. mod(Step, max(NStepsMC,1)) ==0 ) ) then
    
          ! Here we sum up the NStateWF over all processes dealing with a specific component to improve statistics
#if MPI_VER > 0
          call MPI_Allreduce(pc%NStateWF, tempVec(0:pc%NFluctMax), size(pc%NStateWF), MPI_INTEGER, &
&           MPI_SUM, Communicator, ierror)
          pc%NStateWF = tempVec(0:pc%NFluctMax)
          
          call MPI_Reduce( pc%NFluctUpSuccesses(:),tempVec1(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, &
          & MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
          call MPI_Reduce( pc%NFluctUpAttempts(:),tempVec2(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, &
          & MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
          call MPI_Reduce( pc%NFluctDownSuccesses(:),tempVec3(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, &
          & MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
          call MPI_Reduce( pc%NFluctDownAttempts(:),tempVec4(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, &
          & MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
          
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
          j = pc%NFluctMax
          write( IOBuffer, '(I8, I12, F15.2, 2F10.4)' ) j, pc%NStateWF(j), &
&           pc%WF(j), 0._RK, real(tempVec3(j), RK) / &
&             real(tempVec4(j), RK) * 100._RK
            call LogWrite

          do j = pc%NFluctMax - 1, 1, -1
            write( IOBuffer, '(I8, I12, F15.2, 2F10.4)' ) j, pc%NStateWF(j), &
&             pc%WF(j), real(tempVec1(j+1), RK) / &
&               real(tempVec2(j+1), RK) * 100._RK, &
&             real(tempVec3(j), RK) / &
&               real(tempVec4(j), RK) * 100._RK
            call LogWrite
          end do
          write( IOBuffer, &
&           '(I8, I12, F15.2, 2F10.4)' ) 0, pc%NStateWF(0), pc%WF(0), &
&           real(tempVec1(1), RK) / &
&             real(tempVec2(1), RK) * 100._RK, 0._RK
          call LogWrite
          call LogWriteBlank
          pc%NStateWF(:) = 0
          
          
#else        
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
          j = pc%NFluctMax
          write( IOBuffer, '(I8, I12, F15.2, 2F10.4)' ) j, pc%NStateWF(j), &
&           pc%WF(j), 0._RK, real(pc%NFluctDownSuccesses(j), RK) / &
&             real(pc%NFluctDownAttempts(j), RK) * 100._RK
            call LogWrite

          do j = pc%NFluctMax - 1, 1, -1
            write( IOBuffer, '(I8, I12, F15.2, 2F10.4)' ) j, pc%NStateWF(j), &
&             pc%WF(j), real(pc%NFluctUpSuccesses(j+1), RK) / &
&               real(pc%NFluctUpAttempts(j+1), RK) * 100._RK, &
&             real(pc%NFluctDownSuccesses(j), RK) / &
&               real(pc%NFluctDownAttempts(j), RK) * 100._RK
            call LogWrite
          end do
          write( IOBuffer, &
&           '(I8, I12, F15.2, 2F10.4)' ) 0, pc%NStateWF(0), pc%WF(0), &
&           real(pc%NFluctUpSuccesses(1), RK) / &
&             real(pc%NFluctUpAttempts(1), RK) * 100._RK, 0._RK
          call LogWrite
          call LogWriteBlank
          pc%NStateWF(:) = 0
          
          
#endif          
        end if

      ! Chemical potential by Widom's test particle method
      case( ChemPotMethodWidom )
        pc%CalcChemPot = .true.
        call Mol2AtomTest( this%Component(i), this%Component(i)%NTest )
        this%EPotTest(:) = this%Density * pc%EPotTestCorrLJ &
&                                       + pc%EPotTestCorrRF
        do j = 1, this%NRealComponents
          call ChemicalPotential( this%Interaction( i, j ), &
&                                 this%EPotTest, this%BoxLength )
        end do
        ChemPot = sum( exp( -( this%EPotTest(:) ) / this%Temperature ) ) &
&                   / pc%NTestAll
#if MPI_VER > 0
        if ( SimulationType .ne. MonteCarlo .or. (Equilibration .and. CommonEqui) ) then
          call MPI_Reduce( ChemPot, pc%ChemPot, 1, &
&           MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
        else
            pc%ChemPot = ChemPot
        endif
#else
        pc%ChemPot = ChemPot
#endif

      case default
        pc%CalcChemPot = .false.
        pc%ChemPot = 0._RK
      end select

    end do  componentLoop
#if MPI_VER > 0
!    call SetCommunicator( tempComm )
#endif
    contains



    function tprnd( l_range, h_range ) result( rharvest )

      implicit none

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
        if ( this%OptPressure ) then
          pi%Virial(1:n1, 1:n2) = pi%VirialNew(1:n1, 1:n2)
        end if
      end do
    end do

  end subroutine TEnsemble_UpdateEnergy



!==============================================================!
!  Subroutine TEnsemble_UpdateEnergy1                          !
!==============================================================!

  subroutine TEnsemble_UpdateEnergy1( this, nc, np )

    implicit none

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
      if ( this%OptPressure ) then
        pi%Virial(np, 1:n) = pi%Virial1(1:n)
      end if
      this%Interaction(i, nc)%EPot(1:n, np) = pi%EPot1(1:n)
      if ( this%OptPressure ) then
        this%Interaction(i, nc)%Virial(1:n, np) = pi%Virial1(1:n)
      end if
    end do

  end subroutine TEnsemble_UpdateEnergy1



!==============================================================!
!  Subroutine TEnsemble_Energy                                 !
!==============================================================!

  subroutine TEnsemble_Energy( this, E )

    implicit none

    ! Declare arguments
    type(TEnsemble)       :: this
    real(RK), intent(out) :: E

    ! Declare local variables
    type(TInteraction), pointer :: pi
    integer                     :: nc, np
    integer                     :: i, n

    ! Initialize new energy
    E = 0._RK

    ! Loop over components
    do nc = 1, this%NComponents
      do i = 1, this%NComponents
        pi => this%Interaction(nc, i)
        n = pi%NPart2

        ! Loop over particles
        do np = 1, this%Component(nc)%NPart
          call Energy( pi, np, this%BoxLength )

          ! Save new energy matrix
          pi%EPotNew(np, 1:n) = pi%EPot1(1:n)
          if ( this%OptPressure ) then
            pi%VirialNew(np, 1:n) = pi%Virial1(1:n)
          end if

          ! Sum energy
          E = E + sum( pi%EPot1(1:n) )
        end do
      end do
    end do

    ! Calculate new energy
    E = .5_RK * E + this%Density * this%EPotCorrLJ + this%EPotCorrRF

  end subroutine TEnsemble_Energy



!==============================================================!
!  Subroutine TEnsemble_Energy1                                !
!==============================================================!

  subroutine TEnsemble_Energy1( this, nc, np, EPotNew )

    implicit none

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

    ! Loop over components
    do i = 1, this%NComponents
      pi => this%Interaction(nc, i)
      n = pi%NPart2

      call Energy( pi, np, this%BoxLength )

      ! Calculate new energy
      EPotNew = EPotNew + sum( pi%EPot1(1:n) )
    end do

  end subroutine TEnsemble_Energy1



!==============================================================!
!  Function TEnsemble_GetEnergy                                !
!==============================================================!

  function TEnsemble_GetEnergy( this ) result(E)

    implicit none

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

  end function TEnsemble_GetEnergy



!==============================================================!
!  Function TEnsemble_GetEnergy1                               !
!==============================================================!

  function TEnsemble_GetEnergy1( this, nc, np ) result(E)

    implicit none

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

  end function TEnsemble_GetEnergy1



!==============================================================!
!  Function TEnsemble_GetVirial                                !
!==============================================================!

  function TEnsemble_GetVirial( this ) result(V)

    implicit none

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
    type(TComponent), pointer :: pc
    integer                   :: i
    real(RK)                  :: EPotDelta
    logical                   :: accepted

    ! Assign local variables
    pc => this%Component(nc)

    ! Update number of move attempts
    pc%NMoveAttempts = pc%NMoveAttempts + 1

    ! Save current particle position and energy
    r(:) = pc%P0(np, :)
    EPotOld = GetEnergy( this, nc, np )

    ! Generate a trial displacement
    do i = 1, 3
      pc%P0(np, i) = pc%P0(np, i) + rnd( -pc%DispTran, pc%DispTran )
    end do

    ! Apply periodic boundary conditions
    pc%P0(np, :) = pc%P0(np, :) - anint( pc%P0(np, :) )

    ! Convert molecular coordinates to atom positions
    call Mol2Atom1( pc, np )

    ! Calculate particle energy at trial position
    call Energy( this, nc, np, EPotNew )

    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    if (Equilibration .and. CommonEqui) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, &
&       MPI_RK, MPI_SUM, Communicator, ierror )
    else
      EPotDelta = EPotOld - EPotNew
    endif 
#else
    EPotDelta = EPotOld - EPotNew

#endif

    accepted = EPotDelta > 0._RK
    if( .not. accepted ) accepted = exp( EPotDelta / this%Temperature ) > rnd( 0._RK, 1._RK )
    if( accepted ) then

      ! Accept move
      pc%NMoveSuccesses = pc%NMoveSuccesses + 1
      call UpdateEnergy( this, nc, np )

    else

      ! Reject move
      pc%P0(np, :) = r(:)
      call Mol2Atom1( pc, np )

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
    type(TComponent), pointer :: pc
    integer                   :: i
    real(RK)                  :: EPotDelta
    logical                   :: accepted

    ! Assign local variables
    pc => this%Component(nc)

    ! Update number of rotation attempts
    pc%NRotateAttempts = pc%NRotateAttempts + 1

    ! Save current particle orientation and energy
    q(:) = pc%Q0(np, :)
    EPotOld = GetEnergy( this, nc, np )

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

    ! Calculate particle energy with trial orientation
    call Energy( this, nc, np, EPotNew )

    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    if (Equilibration .and. CommonEqui) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, &
&       MPI_RK, MPI_SUM, Communicator, ierror )
    else
      EPotDelta = EPotOld - EPotNew
    endif 
#else
    EPotDelta = EPotOld - EPotNew

#endif

    accepted = EPotDelta > 0._RK
    if( .not. accepted ) accepted = exp( EPotDelta / this%Temperature ) > rnd( 0._RK, 1._RK )
    if( accepted ) then

      ! Accept rotation
      pc%NRotateSuccesses = pc%NRotateSuccesses + 1
      call UpdateEnergy( this, nc, np )

    else

      ! Reject rotation
      pc%Q0(np, :) = q(:)
      call Mol2Atom1( pc, np )

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
    type(TComponent), pointer :: pc, pcf
    integer                   :: i
    real(RK)                  :: EPotDelta
    logical                   :: accepted

    ! Test for fluctuating particle
    if( nc .eq. ncf .and. np .eq. npf ) return

    ! Assign local variables
    pc => this%Component(nc)
    pcf => this%Component(ncf)

    ! Save current particle position and energy
    r(:) = pc%P0(np, :)
    EPotOld = GetEnergy( this, nc, np )

    ! Apply distance criterion
    dr(:) = r(:) - pcf%P0(npf, :)
    dr(:) = ( dr(:) - anint( dr(:) ) ) * this%BoxLength
    f1 = 1._RK / ( dr(1)**2 + dr(2)**2 + dr(3)**2 )**2
    if( rnd(0._RK, 1._RK) < (1._RK - f1) ) return

    ! Update number of move attempts
    pc%NMoveBiasedAttempts = pc%NMoveBiasedAttempts + 1

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

    ! Calculate particle energy at trial position
    call Energy( this, nc, np, EPotNew )

    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    if (Equilibration .and. CommonEqui) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, &
&       MPI_RK, MPI_SUM, Communicator, ierror )
    else
        EPotDelta = EPotOld - EPotNew
    endif 
#else
    EPotDelta = EPotOld - EPotNew

#endif

    accepted = EPotDelta > 0._RK
    if( .not. accepted ) accepted = exp( EPotDelta / this%Temperature ) > rnd( 0._RK, 1._RK )
    if( accepted ) then

      ! Accept move
      pc%NMoveBiasedSuccesses = pc%NMoveBiasedSuccesses + 1
      call UpdateEnergy( this, nc, np )

    else

      ! Reject move
      pc%P0(np, :) = r(:)
      call Mol2Atom1( pc, np )

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
    type(TComponent), pointer :: pc, pcf
    integer                   :: i
    real(RK)                  :: EPotDelta
    logical                   :: accepted

    ! Test for fluctuating particle
    if( nc .eq. ncf .and. np .eq. npf ) return

    ! Assign local variables
    pc => this%Component(nc)
    pcf => this%Component(ncf)

    ! Save current particle orientation and energy
    q(:) = pc%Q0(np, :)
    EPotOld = GetEnergy( this, nc, np )

    ! Apply distance criterion
    dr(:) = pc%P0(np, :) - pcf%P0(npf, :)
    dr(:) = ( dr(:) - anint( dr(:) ) ) * this%BoxLength
    f1 = 1._RK / ( dr(1)**2 + dr(2)**2 + dr(3)**2 )**2
    if( rnd(0._RK, 1._RK) < (1._RK - f1) ) return
!     if( rnd(0._RK, 1._RK) > f1 ) return

    ! Update number of rotation attempts
    pc%NRotateBiasedAttempts = pc%NRotateBiasedAttempts + 1

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

    ! Calculate particle energy with trial orientation
    call Energy( this, nc, np, EPotNew )

    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    if (Equilibration .and. CommonEqui) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, &
&       MPI_RK, MPI_SUM, Communicator, ierror )
    else
      EPotDelta = EPotOld - EPotNew
    endif 
#else
    EPotDelta = EPotOld - EPotNew

#endif

    accepted = EPotDelta > 0._RK
    if( .not. accepted ) accepted = exp( EPotDelta / this%Temperature ) > rnd( 0._RK, 1._RK )
    if( accepted ) then

      ! Accept rotation
      pc%NRotateBiasedSuccesses = pc%NRotateBiasedSuccesses + 1
      call UpdateEnergy( this, nc, np )

    else

      ! Reject rotation
      pc%Q0(np, :) = q(:)
      call Mol2Atom1( pc, np )

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

    ! Declare local variables
    type(TComponent), pointer :: pc, pcf, pcfnew
    integer                   :: oldstate, newstate
    integer                   :: ncfnew, npfnew
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EPotDelta

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

    ! Convert molecular coordinates to atom positions
    call Mol2Atom1( pcfnew, npfnew )

    ! Calculate particle energy at new fluctuating state
    call Energy( this, ncfnew, npfnew, EPotNew )

    ! Apply acceptance criterion
#if MPI_VER > 0
    if (Equilibration .and. CommonEqui) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, &
&       MPI_RK, MPI_SUM, Communicator, ierror )
    else
      EPotDelta = EPotOld - EPotNew
    endif 
#else

     EPotDelta = EPotOld - EPotNew

#endif


    ! Apply long range corrections
    EPotDelta = EPotDelta &
&     + this%Density * ( pcf%EPotTestCorrLJ - pcfnew%EPotTestCorrLJ ) &
&     + pcf%EPotTestCorrRF - pcfnew%EPotTestCorrRF
    if( exp( ( EPotDelta ) / this%Temperature ) > &
&       pc%WF(oldstate) / pc%WF(newstate) * rnd( 0._RK, 1._RK ) ) then

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

    end if

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
#if MPI_VER > 0
    real(RK)                  :: EPotInsAll
#endif

    ! Assign local variables
    pc => this%Component(nc)

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
    if ( tooManyParticles ) return
    np = pc%NPart
    this%NPart = this%NPart + 1

    ! Convert molecular coordinates to atom positions
    call Mol2Atom1( pc, np )

    ! Calculate particle energy at trial position
    call Energy( this, nc, np, EPotIns )

    ! Apply acceptance criterion
#if MPI_VER > 0
    if ( SimulationType .ne. MonteCarlo .or. (Equilibration .and. CommonEqui) ) then
      call MPI_Allreduce( EPotIns, EPotInsAll, 1, &
&       MPI_RK, MPI_SUM, Communicator, ierror )
      EPotInsAll = EPotInsAll + this%Density * pc%EPotTestCorrLJ &
&                             + pc%EPotTestCorrRF
!DEBUG
!  write(0, '(I2, ": EPotIns = ", F12.6)') NProc, EPotInsAll
!DEBUG
    else
      EPotInsAll = EPotIns + this%Density * pc%EPotTestCorrLJ &
&                 + pc%EPotTestCorrRF
    endif 

    ! check if (pc%ChemPot-EPotInsAll)/this%Temperature>exp_arg_max ?
    if( exp( pc%ChemPot - EPotInsAll / this%Temperature ) > &
&       np / this%Volume0 * rnd( 0._RK, 1._RK ) ) then
#else
    EPotIns = EPotIns + this%Density * pc%EPotTestCorrLJ &
&                     + pc%EPotTestCorrRF

    ! check if (pc%ChemPot-EPotIns)/this%Temperature>exp_arg_max ?
    if( exp( pc%ChemPot - EPotIns  / this%Temperature ) > &
&       np / this%Volume0 * rnd( 0._RK, 1._RK ) ) then
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

    else

      ! Reject Insertion
      call RemoveParticle( pc, np )
      this%NPart = this%NPart - 1

    end if

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

    ! Assign local variables
    pc => this%Component(nc)

    ! Update number of delete attempts
    this%NDeleteAttempts = this%NDeleteAttempts + 1

    ! Calculate particle energy
#if MPI_VER > 0
    if ( SimulationType .ne. MonteCarlo .or. (Equilibration .and. CommonEqui) ) then
      call MPI_Allreduce( GetEnergy( this, nc, np ), EPotDel, 1, &
&       MPI_RK, MPI_SUM, Communicator, ierror )
    else
      EPotDel = GetEnergy( this, nc, np )
    endif 
#else
    EPotDel = GetEnergy( this, nc, np )
#endif
    EPotDel = EPotDel + this%Density * pc%EPotTestCorrLJ &
&                     + pc%EPotTestCorrRF

    ! Apply acceptance criterion
    ! check if EPotDel/this%Temperature-pc%ChemPot>exp_arg_max ?
    if( exp( EPotDel / this%Temperature - pc%ChemPot ) > &
&       rnd( 0._RK, 1._RK ) / this%Density / pc%Fraction ) then

      ! Accept Deletion
      this%NDeleteSuccesses = this%NDeleteSuccesses + 1
      call RemoveParticle( pc, np )

      ! Copy energies and virial
      n1 = pc%NPart + 1
      do i = 1, this%NComponents
        pi => this%Interaction(nc, i)
        n2 = pi%NPart2
        pi%EPot(np, 1:n2) = pi%EPot(n1, 1:n2)
        if ( this%OptPressure ) then
          pi%Virial(np, 1:n2) = pi%Virial(n1, 1:n2)
        end if
        this%Interaction(i, nc)%EPot(1:n2, np) = pi%EPot(n1, 1:n2)
        if ( this%OptPressure ) then
          this%Interaction(i, nc)%Virial(1:n2, np) = pi%Virial(n1, 1:n2)
        endif
      end do

      ! Zero diagonal elements
      this%Interaction(nc, nc)%EPot(np, np) = 0._RK
      if ( this%OptPressure ) then
        this%Interaction(nc, nc)%Virial(np, np) = 0._RK
      end if

      this%NPart = this%NPart - 1

      ! Update density
      this%Density = this%NPart / this%Volume0

      ! Update fractions and NDF
      call UpdateFractions( this )

      ! Update long range correction
      call CalculateCorr( this )

    end if

  end subroutine TEnsemble_Delete



!==============================================================!
!  Subroutine TEnsemble_Move2End                               !
!==============================================================!

  subroutine TEnsemble_Move2End( this, nc, np )

    implicit none

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
      if ( this%OptPressure ) then
        VSave(1:n2) = pi%Virial(np, :)
      end if
      if( i .eq. nc ) then
        ESave(np) = pi%EPot(np, n2)
        if ( this%OptPressure ) then
          VSave(np) = pi%Virial(np, n2)
        end if
      end if
      pi%EPot(np, :) = pi%EPot(n1, :)
      this%Interaction(i, nc)%EPot(:, np) = pi%EPot(n1, :)
      pi%EPot(n1, :) = ESave(1:n2)
      this%Interaction(i, nc)%EPot(:, n1) = ESave(1:n2)
      if ( this%OptPressure ) then
        pi%Virial(np, :) = pi%Virial(n1, :)
        this%Interaction(i, nc)%Virial(:, np) = pi%Virial(n1, :)
        pi%Virial(n1, :) = VSave(1:n2)
        this%Interaction(i, nc)%Virial(:, n1) = VSave(1:n2)
      end if
    end do

    ! Zero diagonal elements
    this%Interaction(nc, nc)%EPot(np, np) = 0._RK
    if ( this%OptPressure ) then
      this%Interaction(nc, nc)%Virial(np, np) = 0._RK
    end if

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
#if MPI_VER > 0
    real(RK) :: EPotNew
#endif
    logical  :: accepted

    ! Update number of resizing attempts
    this%NResizeAttempts = this%NResizeAttempts + 1

    ! Save current simulation box size, volume, energy, virial
    VolumeOld = this%Volume0
    EPotOld = this%EPot

    ! Generate a trial volume change
    this%Volume0 = this%Volume0 * (1._RK + rnd( -this%DispVol, this%DispVol ))
    call UpdateBoxLength( this )

    ! Convert molecular coordinates to atom positions
    call Mol2Atom( this )

    ! Calculate potential energy and virial at trial position
#if MPI_VER > 0
    if ( SimulationType .ne. MonteCarlo .or. (Equilibration .and. CommonEqui) ) then
      call Energy( this, EPotNew )
      call MPI_Allreduce( EPotNew, this%EPot, 1, &
&       MPI_RK, MPI_SUM, Communicator, ierror )
    else
     call Energy( this, this%EPot )
    endif
#else
    call Energy( this, this%EPot )
#endif

    ! Find potential change
    EPotDelta = this%RefPressure * (this%Volume0 - VolumeOld) &
&     + this%EPot - EPotOld &
&     + this%NPart * this%Temperature * log( VolumeOld / this%Volume0 )

    ! Apply Metropolis acceptance criterion
    accepted = EPotDelta < 0._RK
    if ( .not. accepted ) accepted = exp( -EPotDelta / this%Temperature ) > rnd( 0._RK, 1._RK )
    if( accepted ) then

      ! Accept volume change
      this%NResizeSuccesses = this%NResizeSuccesses + 1

      ! Update energy and virial matrices
      call UpdateEnergy( this )
      if ( this%OptPressure ) then
#if MPI_VER > 0
        if ( SimulationType .ne. MonteCarlo .or. (Equilibration .and. CommonEqui) ) then
          call MPI_Allreduce( GetVirial( this ), this%Virial, 1, &
&           MPI_RK, MPI_SUM, Communicator, ierror )
        else
          this%Virial = GetVirial( this )
        endif
#else
        this%Virial = GetVirial( this )
#endif
      end if

    else

      ! Reject volume change
      this%Volume0 = VolumeOld
      call UpdateBoxLength( this )
      call Mol2Atom( this )
      this%EPot = EPotOld

    end if

  end subroutine TEnsemble_Resize



!==============================================================!
!  Subroutine TEnsemble_ZeroNAttempts                          !
!==============================================================!

  subroutine TEnsemble_ZeroNAttempts( this )

    implicit none

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

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i

    ! Save current state
    do i = 1, this%NRealComponents
      call SaveState( this%Component(i) )
    end do

  end subroutine TEnsemble_SaveState



!==============================================================!
!  Subroutine TEnsemble_RestoreState                           !
!==============================================================!

  subroutine TEnsemble_RestoreState( this )

    implicit none

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

  end subroutine TEnsemble_RestoreState



!==============================================================!
!  Subroutine TEnsemble_ResultOpen                             !
!==============================================================!

  subroutine TEnsemble_ResultOpen( this )

    implicit none

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
#if TRANS ==1
      write( IOBuffer, '(I16)' ) this%EnsembleNumber
      call FileAppend( this%iounit_rescf, &
&       trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ResultTransportExtension )

#endif

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
#if  TRANS == 1
      ! Open result file for correlation function
      write( IOBuffer, '(I16)' ) this%EnsembleNumber
      call FileRewrite( this%iounit_rescf, &
&       trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ResultTransportExtension )
!TRANSPORT_END
#endif
    end if

  end subroutine TEnsemble_ResultOpen



!==============================================================!
!  Subroutine TEnsemble_ResultUpdate                           !
!==============================================================!

  subroutine TEnsemble_ResultUpdate( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    type(TComponent), pointer :: pc
    integer                   :: i
    real(RK)                  :: value
#if TRANS ==1
    integer                   :: j
#endif

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

      ! 4.) Chemical potential and partial mole volumes
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
      write( IOBuffer, '("     PRESS")' )
      call FileWriteNoAdvance( this%iounit_result )
      call FileWriteNoAdvance( this%iounit_runave )

      ! Density
      write( IOBuffer, '("   DENSITY")' )
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

      ! Partial mole volume
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

        ! Mole fraction of each component
        do i = 1, this%NComponents
          write( IOBuffer, '("   FRACT", I2)' ) i
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )
        end do
      end if

      call FileWriteBlank( this%iounit_result )
      call FileWriteBlank( this%iounit_runave )
    end if

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
#if  TRANS == 1
    ! 4.) Tranport properties !TRANSPORT_start
    if( mod( Step - this%Ncorr, BlockSizeCF * this%NSpancf ) == 0 .and. &
&     (this%Mmess > 0) ) then
      do i = 1, this%NComponents
        call UpdateCF( this%Sumself_i(i), this%selfd_i(i), this%Mmess  )
      end do

      if(this%Ncomponents == 2) then
        call UpdateCF( this%SumBin_d, this%binary_d, this%Mmess )
      end if

     if(this%Ncomponents == 3 ) then
        call UpdateCF( this%Sumter_a, this%ternary_a, this%Mmess )
        call UpdateCF( this%Sumter_b, this%ternary_b, this%Mmess )
        call UpdateCF( this%Sumter_c, this%ternary_c, this%Mmess )
      end if

      call UpdateCF( this%SumVisco_s, this%visco_s, this%Mmess )
      call UpdateCF( this%SumVisco_b, this%visco_b, this%Mmess )
      call UpdateCF( this%SumConduct, this%conduct, this%Mmess )
    end if
!TRANSPORT_END
#endif
    ! 4.) Chemical potential and partial mole volumes
    do i = 1, this%NRealComponents
      pc => this%Component(i)
      if( pc%CalcChemPot ) then
        select case( pc%ChemPotMethod )
        case( ChemPotMethodGradins )
#if MPI_VER > 0
          ! Per Process we calculate GI only for one component 
          if (mod(NProc,this%NGradInsComp)/=pc%NGradThis) cycle
#endif
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
#if MPI_VER > 0
          ! Per Process we calculate GI only for one component 
          if (mod(NProc,this%NGradInsComp)/=pc%NGradThis) cycle 
#endif
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
      if( SimulationType .eq. MolecularDynamics ) then
        write( IOBuffer, '(F10.5)' ) this%SumPressure%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(F10.5)' ) this%SumPressure%Average
        call FileWriteNoAdvance( this%iounit_runave )
      else
        if ( this%OptPressure ) then
          write( IOBuffer, '(F10.5)' ) this%SumPressure%BlockAverage
          call FileWriteNoAdvance( this%iounit_result )
          write( IOBuffer, '(F10.5)' ) this%SumPressure%Average
          call FileWriteNoAdvance( this%iounit_runave )
        else
          write( IOBuffer, '(F10.5)' ) this%RefPressure
          call FileWriteNoAdvance( this%iounit_result )
          write( IOBuffer, '(F10.5)' ) this%RefPressure
          call FileWriteNoAdvance( this%iounit_runave )
        end if
      end if
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
            if( pc%NPart > 1 ) then
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
              select case( pc%ChemPotMethod )
              case( ChemPotMethodGradIns )
                write( IOBuffer, '(F10.5)' ) &
&                 log( pc%SumInvChemPotRho%BlockAverage )
                call FileWriteNoAdvance( this%iounit_result )
                write( IOBuffer, '(F10.5)' ) &
&                 log( pc%SumInvChemPotRho%Average )
                call FileWriteNoAdvance( this%iounit_runave )
              case( ChemPotMethodWidom )
              write( IOBuffer, '(F10.5)' ) &
&               log( 1._RK / pc%SumChemPotV%BlockAverage )
              call FileWriteNoAdvance( this%iounit_result )
              write( IOBuffer, '(F10.5)' ) &
&               log( 1._RK / pc%SumChemPotV%Average )
              call FileWriteNoAdvance( this%iounit_runave )
              end select
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
              write( IOBuffer, '(F10.4)' ) 0._RK
              call FileWriteNoAdvance( this%iounit_result )
              call FileWriteNoAdvance( this%iounit_runave )
            else
              write( IOBuffer, '(F10.4)' ) pc%SumVW%BlockAverage
              call FileWriteNoAdvance( this%iounit_result )
              write( IOBuffer, '(F10.4)' ) pc%SumVW%Average
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

        ! Mole fraction of each component
        do i = 1, this%NComponents
          pc => this%Component(i)
          write( IOBuffer, '(F10.5)' ) pc%SumFraction%BlockAverage
          call FileWriteNoAdvance( this%iounit_result )
          write( IOBuffer, '(F10.5)' ) pc%SumFraction%Average
          call FileWriteNoAdvance( this%iounit_runave )
        end do
      end if

      call FileWriteBlank( this%iounit_result )
      call FileWriteBlank( this%iounit_runave )
#if ARCH == 2
      call flush( this%iounit_result )
      call flush( this%iounit_runave )
#endif
    end if
#if  TRANS == 1
    ! Transport properties !TRANSPORT_start
    if( ( this%Mmess > 0 ) .and. ( mod(this%Mmess, this%Nviewcf) == 0 ) ) then
      rewind( this%iounit_rescf )
      write( IOBuffer, '("  TIME[ps]")' )
      call FileWriteNoAdvance( this%iounit_rescf )
!      if(this%Ncomponents==2)then
!        write( IOBuffer, '(T11,"D_12")' )
!        call FileWriteNoAdvance( this%iounit_rescf )
!      end if
!      if(this%Ncomponents==3)then
!          write( IOBuffer, '(T10, "D_ijk", I2)') i
!          call FileWriteNoAdvance( this%iounit_rescf )
!      end if
      if(this%Ncomponents>1)then
        do i=1,this%NComponents*this%NComponents
            write( IOBuffer, '(T10, "L_ij", I1)') i
            call FileWriteNoAdvance( this%iounit_rescf )
        end do
      end if
      do i = 1, this%NComponents
        write( IOBuffer, '(T10,"D_i",I2)' ) i
        call FileWriteNoAdvance( this%iounit_rescf )
      end do

      write( IOBuffer, '(T13,"VS")' )
      call FileWriteNoAdvance( this%iounit_rescf )

      write( IOBuffer, '(T13,"VB")' )
      call FileWriteNoAdvance( this%iounit_rescf )

      write( IOBuffer, '(T13,"CO")' )
      call FileWriteNoAdvance( this%iounit_rescf )

      !if( this%Ncomponents == 2 ) then
       ! write( IOBuffer, '(T9,"IntD12")' )
       ! call FileWriteNoAdvance( this%iounit_rescf )
      !end if

      if( this%Ncomponents > 1 ) then
        do i=1,this%NComponents*this%NComponents
           write( IOBuffer, '(T7,"Int_Lij",I1)')i
           call FileWriteNoAdvance( this%iounit_rescf )
       end do
      end if

      do i = 1, this%NComponents
         write( IOBuffer, '(T7,"IntD_i",I2)' ) i
         call FileWriteNoAdvance( this%iounit_rescf )
      end do

      write( IOBuffer, '(T9,"Int VS")' )
      call FileWriteNoAdvance( this%iounit_rescf )

      write( IOBuffer, '(T9,"Int VB")' )
      call FileWriteNoAdvance( this%iounit_rescf )

      write( IOBuffer, '(T10,"Int C ")' )
      call FileWriteNoAdvance( this%iounit_rescf )

      call FileWriteBlank( this%iounit_rescf )

      ! integration time
      do i  = 1, this%Ncorr
        value = TimeStep*UnitTime/1E-12_RK
        write( IOBuffer, '(F10.5)' ) (i-1)*value
        call FileWriteNoAdvance( this%iounit_rescf )

        ! Binary diffusion coefficient
!        if(this%Ncomponents==2)then
!          write( IOBuffer, '(T5,F10.5)' ) this%cf_db(i)/this%cf_db(1)
!          call FileWriteNoAdvance( this%iounit_rescf )
!        end if

!        ! Ternary diffusion coefficient
!        if(this%Ncomponents==3)then
!              write( IOBuffer, '(T5, F10.5)' ) this%lamda(1, i)/this%lamda(1,1)
!              call FileWriteNoAdvance( this%iounit_rescf )
!        end if
        ! Ternary diffusion coefficient
        if(this%Ncomponents>1)then
          do j=1,this%NComponents*this%NComponents
              write( IOBuffer, '(T5, F10.5)' ) this%lamda(j, i)/this%lamda(j,1)
              call FileWriteNoAdvance( this%iounit_rescf )
          end do
        end if

        ! Self-diffusion coefficients
        do j = 1, this%NComponents
          write( IOBuffer, '(T5, F10.5)' ) this%cf_d(j,i)/this%cf_d(j,1)
          call FileWriteNoAdvance( this%iounit_rescf )
        end do

        ! Shear viscosity
        write( IOBuffer, '(T5, F10.5)' ) this%cf_vs(i)/this%cf_vs(1)
        call FileWriteNoAdvance( this%iounit_rescf )

        ! Bulk viscosity
        write( IOBuffer, '(T5, F10.5)' ) this%cf_vb(i)/this%cf_vb(1)
        call FileWriteNoAdvance( this%iounit_rescf )

        ! Thermal conductivity
        write( IOBuffer, '(T5, F10.5)' ) this%cf_c(i)/this%cf_c(1)
        call FileWriteNoAdvance( this%iounit_rescf )

        ! integral ======================================================!
        value = dsqrt(UnitEnergy/UnitMass)*UnitLength/1E-10_RK

        ! Binary diffusion coefficient
 !       if( this%Ncomponents == 2) then
  !        write( IOBuffer, '(T5, F10.4)' ) &
!&           this%sinte_db(i) / this%sinte_db(this%Ncorr) * this%binary_d * value
!          call FileWriteNoAdvance( this%iounit_rescf )
!        end if

 !       if( this%Ncomponents == 3) then
 !            write( IOBuffer, '(T5, F10.4)' ) &
!&            this%sinte_lamda(1,i) / this%sinte_lamda(1,this%Ncorr) 
!             call FileWriteNoAdvance( this%iounit_rescf )
!        end if

        if( this%Ncomponents > 1) then
          do j = 1, this%NComponents*this%NComponents
             write( IOBuffer, '(T5, F10.4)' ) &
&            this%sinte_lamda(j,i) / this%sinte_lamda(j,this%Ncorr)* value
             call FileWriteNoAdvance( this%iounit_rescf )
          end do
        end if

        ! Self-diffusion coefficient
        do j = 1, this%NComponents
          write( IOBuffer, '(T5, F10.4)' ) &
&           this%sinte_i(j,i) / this%sinte_i(j,this%Ncorr) * this%selfd_i(j) * value
          call FileWriteNoAdvance( this%iounit_rescf )
        end do

       !viscosity
        value = dsqrt(UnitEnergy*UnitMass)/UnitLength**2/1E-4_RK

       !shear
        write( IOBuffer, '(T5, F10.5)' ) &
&         this%sinte_vs(i) / this%sinte_vs(this%Ncorr) * this%visco_s * value
        call FileWriteNoAdvance( this%iounit_rescf )

       ! bulk
        write( IOBuffer, '(T5, F10.5)' ) &
&         this%sinte_vb(i) / this%sinte_vb(this%Ncorr) * this%visco_b * value
        call FileWriteNoAdvance( this%iounit_rescf )

       ! thermal conductivity
        value = dsqrt(UnitEnergy/UnitMass)*kBoltzmann/UnitLength**2

        write( IOBuffer, '(T5, F10.5)' ) &
&         this%sinte_c(i) / this%sinte_c(this%Ncorr) * this%conduct * value
        call FileWriteNoAdvance( this%iounit_rescf )
        call FileWriteBlank( this%iounit_rescf )

      end do

#if ARCH == 2
      call flush( this%iounit_rescf )
#endif
    end if
!TRANSPORT_END
#endif
  end subroutine TEnsemble_ResultUpdate



!==============================================================!
!  Subroutine TEnsemble_ResultClose                            !
!==============================================================!

  subroutine TEnsemble_ResultClose( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Close running average result file
    if( .not. SimulationType .eq. SecondVirialCoeff ) &
&     call FileClose( this%iounit_runave )

    ! Close result file
    call FileClose( this%iounit_result )

!TRANSPORT Hier feht Befehl schliessen *.rtr

  end subroutine TEnsemble_ResultClose



!==============================================================!
!  Subroutine TEnsemble_ErrorsUpdate                           !
!==============================================================!

  subroutine TEnsemble_ErrorsUpdate( this )

    implicit none
#if MPI_VER > 0
    include 'mpif.h'
#endif
    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK)                  :: Average, Variance
    type(TComponent), pointer :: pc
    integer                   :: i, j
#if  TRANS == 1
    real(RK)                  :: value
#endif

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

#if MPI_VER > 0
    integer :: tempVal, tempVal2, color
    integer :: tempVec1(this%NFluctMax), tempVec2(this%NFluctMax)
    integer :: tempVec3(this%NFluctMax), tempVec4(this%NFluctMax)
    real(RK) :: tempReal
    if ( SimulationType .eq. MonteCarlo) then
      NBlockSizes = int( sqrt( real( Step*NProcs / BlockSize, RK ) ) )
      tempVal = NBlocks
      NBlocks = tempVal*NProcs
    endif
    
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
  !      if( pc%CalcChemPot ) then
          select case( pc%ChemPotMethod )
          case( ChemPotMethodGradIns )
          
#if MPI_VER > 0          
           if (mod(NProc,this%NGradInsComp)==pc%NGradThis) then
             color = i
           else
             color = 100000
           endif 
            NProc_W = NProc
            RootProc_W = Rootproc
           call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,NProc,Communicator,ierror) 
           ! Careful, Nproc and NProcs are now specific for Communicator
           call SetCommunicator( Communicator )
           NBlockSizes = int( sqrt( real( Step*NProcs / BlockSize, RK ) ) )
           NBlocks = tempVal*NProcs   
           RootProc = NProc_W==(pc%NGradThis)     
           NRootProc_W = (pc%NGradThis)     

#endif          
          
            call ErrorGI( pc%SumInvChemPotRho )
!DEBUG
            call ErrorGI( pc%SumInvChemPotRho1 )
            call ErrorGI( pc%SumInvChemPotRho2 )
!DEBUG
            call ErrorGI( pc%SumVW1 )
            call ErrorGI( pc%SumVW2 )
            call ErrorGI( pc%SumVW )

#if MPI_VER > 0       
!            call MPI_Allreduce (pc%NStateWF(0:pc%NFluctMax), tempVec1, pc%NFluctMax, MPI_INTEGER, MPI_SUM, Communicator, ierror )
!            pc%NStateWF = tempVec1/NProcs
!            call MPI_Allreduce (pc%NState(0:pc%NFluctMax), tempVec1, pc%NFluctMax, MPI_INTEGER, MPI_SUM, Communicator, ierror)
!            pc%NState = tempVec1/NProcs
!            call MPI_Allreduce (pc%ProbW0, tempReal, 1, MPI_RK, MPI_SUM, Communicator, ierror)
!            pc%ProbW0 = tempReal/NProcs
!            call MPI_Allreduce (pc%ProbW1, tempReal, 1, MPI_RK, MPI_SUM, Communicator, ierror)
!            pc%ProbW1 = tempReal/NProcs
!            call MPI_Allreduce (pc%ProbW0V, tempReal, 1, MPI_RK, MPI_SUM, Communicator, ierror)
!            pc%ProbW0V = tempReal/NProcs
!            call MPI_Allreduce (pc%ProbW1Rho, tempReal, 1, MPI_RK, MPI_SUM, Communicator, ierror)
!            pc%ProbW1Rho = tempReal/NProcs
            
                     
            call SetCommunicator(MPI_COMM_WORLD)
            RootProc = Rootproc_W                          
#endif 

          case( ChemPotMethodWidom )
            call Error( pc%SumChemPotV )
            call Error( pc%SumVW )
          case default
           ! DO NOTHING
          end select
!           if( ConstantPressure .and. this%NRealComponents > 1 ) &
! &           call Error( pc%SumVW )

     !   end if
      end do
    end if


#if MPI_VER >0
    if ( SimulationType .eq. MonteCarlo) then
      NBlockSizes = int( sqrt( real( Step / BlockSize, RK ) ) )
      NBlocks = tempVal
    end if
#endif

    ! Open final result file
    write( IOBuffer, '(I16)' ) this%EnsembleNumber
    call FileRewrite( this%iounit_errors, &
&     trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ErrorsFileExtension )

    write( IOBuffer, '(76("="))')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("*                           Publishing with ms2                            *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* Every user agrees to cite ms2 upon usage as follows                      *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* ------------------------------------------------------------------------ *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* S. Deublein, B. Eckl, J. Stoll, S. Lishchuk, G. Guevara-Carrion,         *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* C.W. Glass, T. Merker, M. Bernreuther, H. Hasse, J. Vrabec               *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* Computer Physics Communications (2011)                                   *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* DOI:10.1016/j.cpc.2011.04.026                                        *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(76("="))')
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

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

    if ( SimulationType .eq. MonteCarlo .and. (Nproc == NRootProc)) then
      ! The RootProc receives data from all processes and therefore the # of 
      ! Step is increased accordingly
      write( IOBuffer, '("Number of production steps", T36, ":", I10)' ) &
&       Step*NProcs
    else 
      write( IOBuffer, '("Number of production steps", T36, ":", I10)' ) &
&       Step
    end if

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
      write( IOBuffer, '("Acceptance rate", T34, ":", F20.9)' ) &
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
        write( IOBuffer, '("Mole fraction of ", A, T36, ":", F20.9)' ) &
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
    if( SimulationType .eq. MolecularDynamics ) then
       Average = this%SumPressure%Average
       Variance = this%SumPressure%Variance
       write( IOBuffer, '("Pressure", T29, "reduced:", 2F20.9)' ) &
&        Average, Variance
       call FileWrite( this%iounit_errors )
       write( IOBuffer, '(T30, "in MPa:", 2F20.9)' ) &
&        Average * UnitPressure * 1E-6_RK, Variance * UnitPressure * 1E-6_RK
       call FileWrite( this%iounit_errors )
       call FileWriteBlank( this%iounit_errors )
    else
       if ( this%OptPressure ) then
         Average = this%SumPressure%Average
         Variance = this%SumPressure%Variance
       else
         Average = this%RefPressure
         Variance = 0._RK
       end if
       write( IOBuffer, '("Pressure", T29, "reduced:", 2F20.9)' ) &
&        Average, Variance
       call FileWrite( this%iounit_errors )
       write( IOBuffer, '(T30, "in MPa:", 2F20.9)' ) &
&        Average * UnitPressure * 1E-6_RK, Variance * UnitPressure * 1E-6_RK
       call FileWrite( this%iounit_errors )
       call FileWriteBlank( this%iounit_errors )
    end if

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
      ! Mole fraction
      do i = 1, this%NComponents
        pc => this%Component(i)
        Average = pc%SumFraction%Average
        Variance = pc%SumFraction%Variance
        write( IOBuffer, &
&         '("Mole fraction of ", A, T36, ":", 2F20.9)' ) &
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
          if( pc%NPart > 1 ) then

            Variance = pc%SumInvChemPotRho%Variance / pc%SumInvChemPotRho%Average
            Average = log( pc%Fraction * pc%SumInvChemPotRho%Average )
            write( IOBuffer, &
&             '("Chemical potential of ", A, T33, "r`d:", 2F20.9)' ) &
&             trim( this%Component(i)%Molecule%PotModFileName ), &
&             Average, Variance
            call FileWrite( this%iounit_errors )
!DEBUG
            Variance = pc%SumInvChemPotRho1%Variance / pc%SumInvChemPotRho1%Average
            Average = log( pc%Fraction * pc%SumInvChemPotRho1%Average )
            write( IOBuffer, &
&             '("Chemical potential 0 of ", A, T33, "r`d:", 2F20.9)' ) &
&             trim( this%Component(i)%Molecule%PotModFileName ), &
&             Average, Variance
            call FileWrite( this%iounit_errors )
            Variance = pc%SumInvChemPotRho2%Variance / pc%SumInvChemPotRho2%Average
            Average = log( pc%Fraction * pc%SumInvChemPotRho2%Average )
            write( IOBuffer, &
&             '("Chemical potential 1 of ", A, T33, "r`d:", 2F20.9)' ) &
&             trim( this%Component(i)%Molecule%PotModFileName ), &
&             Average, Variance
            call FileWrite( this%iounit_errors )
!DEBUG
!MERKER
          else
            Variance = pc%SumInvChemPotRho%Variance / pc%SumInvChemPotRho%Average
            Average = -log( 1/pc%SumInvChemPotRho%Average )
            write( IOBuffer, &
&             '("Chem. pot. at inf. dilution of ", A, T33, "r`d:", 2F20.9)' ) &
&             trim( this%Component(i)%Molecule%PotModFileName ), &
&             Average, Variance
            Average = this%Temperature*pc%SumInvChemPotRho%Average
            write( IOBuffer, &
&             '("Henrys law constant of ", A, T33, "r`d:", 2F20.9)' ) &
&             trim( pc%Molecule%PotModFileName ), Average, Variance
            call FileWrite( this%iounit_errors )
            write( IOBuffer, '(T30, "in MPa:", 2F20.9)' ) &
&           Average * UnitPressure * 1E-6_RK, Variance * UnitPressure * 1E-6_RK
          end if
          call FileWrite( this%iounit_errors )
!MERKER
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
        ! CP - subtract ideal gas contribution of the pressure
        Average = this%SumCP%Average - 1._RK
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
        write( IOBuffer, '("Isochoric heat capacity", T29, "reduced:", 2F20.9)' ) &
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
#if  TRANS == 1
    ! Transport properties !TRANSPORT_start
    if (CorrfunMode .eq. active) Then

      write( IOBuffer, '(T24, "TRANSPORT PROPERTIES")' )
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      write( IOBuffer, '("VALUE", T31, "UNITS", T46, "AVERAGE", T66, "ERROR")' )
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '("-----", T31, "-----", T46, "-------", T66, "-----")' )
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      write( IOBuffer, '("Number of ACF", T36, ":",T46, I5 )' ) this%Mmess
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      value = this%Ncorr*TimeStep
      write( IOBuffer, '("Length ACF  ", T29, "reduced:", F20.9)' ) value
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T31, "in ps:", F20.9)' )  value*UnitTime/1E-12_RK
      call FileWrite( this%iounit_errors )

      call FileWriteBlank( this%iounit_errors )

      value = this%NSpancf*TimeStep
      write( IOBuffer, '("Time span between ACF ", T29, "reduced:", F20.9)' ) value
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T31, "in ps:", F20.9)' )  value*UnitTime/1E-12_RK
      call FileWrite( this%iounit_errors )

      call FileWriteBlank( this%iounit_errors )

      if( this%Mmess > 0 ) then

        if( this%NComponents == 2  ) then

          if((NBlockSizesCF >= 2 ).and.(NBlocksCF.le.NBlocksMaxCF))then
          call ErrorCF(this%SumBin_d, this%Mmess)
          Average  = this%SumBin_d%Average
          Variance = this%SumBin_d%Variance
          else
          Average  = this%SumBin_d%Average
          Variance = this%SumBin_d%Variance
          end if
          value = dsqrt(UnitEnergy/UnitMass)*UnitLength/1E-10_RK
          write( IOBuffer, '("Binary diff. coeff.", T29, "reduced:", 2F20.9)' ) this%binary_d, Variance
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' )  this%binary_d*value, Variance*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
        end if


       if( this%NComponents == 3 ) then
          if((NBlockSizesCF >= 2 ).and.(NBlocksCF.le.NBlocksMaxCF))then
           call ErrorCF(this%SumTer_a, this%Mmess)
            Average  = this%SumTer_a%Average
            Variance = this%SumTer_a%Variance
           else
            Average  = this%SumTer_a%Average
            Variance = this%SumTer_a%Variance
          end if
          value = dsqrt(UnitEnergy/UnitMass)*UnitLength/1E-10_RK
          write( IOBuffer, '("Ternary diff. coeff. 1_3", T29, "reduced:", 2F20.9)' ) this%ternary_a, Variance
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' )  this%ternary_a*value, Variance*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
       end if

       if( this%NComponents == 3 ) then
          if((NBlockSizesCF >= 2 ).and.(NBlocksCF.le.NBlocksMaxCF))then
           call ErrorCF(this%SumTer_b, this%Mmess)
            Average  = this%SumTer_b%Average
            Variance = this%SumTer_b%Variance
           else
            Average  = this%SumTer_b%Average
            Variance = this%SumTer_b%Variance
          end if
          value = dsqrt(UnitEnergy/UnitMass)*UnitLength/1E-10_RK
          write( IOBuffer, '("Ternary diff. coeff. 1_2", T29, "reduced:", 2F20.9)' ) this%ternary_b, Variance
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' )  this%ternary_b*value, Variance*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
       end if

       if( this%NComponents == 3 ) then
          if((NBlockSizesCF >= 2 ).and.(NBlocksCF.le.NBlocksMaxCF))then
           call ErrorCF(this%SumTer_c, this%Mmess)
            Average  = this%SumTer_c%Average
            Variance = this%SumTer_c%Variance
           else
            Average  = this%SumTer_c%Average
            Variance = this%SumTer_c%Variance
          end if
          value = dsqrt(UnitEnergy/UnitMass)*UnitLength/1E-10_RK
          write( IOBuffer, '("Ternary diff. coeff. 2_3", T29, "reduced:", 2F20.9)' ) this%ternary_c, Variance
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' )  this%ternary_c*value, Variance*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
       end if



         do i = 1, this%NComponents
          if((NBlockSizesCF >= 2 ).and.(NBlocksCF.le.NBlocksMaxCF))then
          call ErrorCF(this%Sumself_i(i), this%Mmess)
          Average  = this%Sumself_i(i)%Average
          Variance = this%Sumself_i(i)%Variance
          else
          Average  = this%Sumself_i(i)%Average
          Variance = this%Sumself_i(i)%Variance
          end if
          value = dsqrt(UnitEnergy/UnitMass)*UnitLength/1E-10_RK
          write( IOBuffer, '("Self-diff. coeff.",A ,T29, "reduced:", 2F20.9)' )  &
&          trim( this%Component(i)%Molecule%PotModFileName ), this%selfd_i(i), Variance
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' )  this%selfd_i(i)*value, Variance*value
          call FileWrite( this%iounit_errors )
         end do
          call FileWriteBlank( this%iounit_errors )

          if((NBlockSizesCF >= 2 ).and.(NBlocksCF.le.NBlocksMaxCF))then
          call ErrorCF(this%SumVisco_s, this%Mmess)
          Average  = this%SumVisco_s%Average
          Variance = this%SumVisco_s%Variance
          else
          Average  = this%SumVisco_s%Average
          Variance = this%SumVisco_s%Variance
          end if
          value = dsqrt(UnitEnergy*UnitMass)/UnitLength**2/1E-4_RK
          write( IOBuffer, '("Shear-Viscosity    ", T29, "reduced:", 2F20.9)' ) this%visco_s, Variance
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T23, "in 10E-4 Pa s:", 2F20.9)' ) this%visco_s*value, Variance*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )

          if((NBlockSizesCF >= 2 ).and.(NBlocksCF.le.NBlocksMaxCF))then
          call ErrorCF(this%SumVisco_b, this%Mmess)
          Average  = this%SumVisco_b%Average
          Variance = this%SumVisco_b%Variance
          else
          Average  = this%SumVisco_b%Average
          Variance = this%SumVisco_b%Variance
          end if
          write( IOBuffer, '("Bulk-Viscosity    ", T29, "reduced:", 2F20.9)' ) this%visco_b, Variance
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T23, "in 10E-4 Pa s:", 2F20.9)' ) this%visco_b*value, Variance*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )

          if((NBlockSizesCF >= 2 ).and.(NBlocksCF.le.NBlocksMaxCF))then
          call ErrorCF(this%SumConduct, this%Mmess)
          Average  = this%SumConduct%Average
          Variance = this%SumConduct%Variance
          else
          Average  = this%SumConduct%Average
          Variance = this%SumConduct%Variance
           end if
          value = dsqrt(UnitEnergy/UnitMass)*kBoltzmann/UnitLength**2
          write( IOBuffer, '("Thermal conductivity ", T29, "reduced:", 2F20.9)' ) this%conduct, Variance
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T23, "in W / (m K) :", 2F20.9)' ) this%conduct*value, Variance*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )

      else

        if(this%NComponents==2)then
          write( IOBuffer, '("Binary diff. coeff.", T29, "reduced:", F20.9)' ) 0._8
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", F20.9)' )  0._8
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
        end if

         do i = 1, this%NComponents
          write( IOBuffer, '("Self-diff. coeff.",A ,T29, "reduced:", F20.9)' )  &
&          trim( this%Component(i)%Molecule%PotModFileName ), 0._8
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", F20.9)' )  0._8
          call FileWrite( this%iounit_errors )
         end do
          call FileWriteBlank( this%iounit_errors )

          write( IOBuffer, '("Shear-Viscosity    ", T29, "reduced:", F20.9)' )  0._8
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T23, "in 10E-4 Pa s:", F20.9)' ) 0._8
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )

          write( IOBuffer, '("Bulk-Viscosity     ", T29, "reduced:", F20.9)' )  0._8
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T23, "in 10E-4 Pa s:", F20.9)' ) 0._8
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )

          write( IOBuffer, '("Thermal conductivity", T29, "reduced:", F20.9)' )  0._8
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T23, "in W / (m K) :", F20.9)' ) 0._8
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )

      end if
!TRANSPORT_END
      ! Separator
      write( IOBuffer, '(76("="))' )
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

    end if
#endif
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
      write( IOBuffer, '("PHASE EQUILIBRIUM DATA")' )
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

      ! Mole fractions of liquid phase
      do i = 1, this%NComponents
        pc => this%Component(i)
        write( IOBuffer, &
&         '("Liquid mole fraction of ", A, T36, ":", F20.9)' ) &
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
    if ( this%OptPressure ) then
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

      ! Mole fractions of vapor phase
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
&         '("Vapor mole fraction of ", A, T36, ":", 2F20.9)' ) &
&         trim( pc%Molecule%PotModFileName ), Average, vary(i)
        call FileWrite( this%iounit_errors )
      end do
      pc => this%Component( this%NComponents )
      Average = pc%SumFraction%Average
      Variance = sqrt( sum( vary(1:(this%NComponents - 1))**2 ) )
      write( IOBuffer, &
&       '("Vapor mole fraction of ", A, T36, ":", 2F20.9)' ) &
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
&       Average, Average * VarPressure / this%SumPressure%Average
!&       Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in mol/l:", 2F20.9)' ) &
&       Average * UnitDensity, Average * VarPressure / this%SumPressure%Average * UnitDensity
!&       Average * UnitDensity, Variance * UnitDensity
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

#if MPI_VER > 0

        call MPI_Reduce( this%NResizeSuccesses,tempVal, 1, MPI_INTEGER, MPI_SUM, &
&         NRootProc, Communicator, ierror )
        call MPI_Reduce( this%NResizeAttempts,tempVal2, 1, MPI_INTEGER, MPI_SUM, &
&         NRootProc, Communicator, ierror )
        if ( Nproc == NRootProc) then
          write( IOBuffer, '("Acceptance rate volume changes", T32, "in %:", F20.9)' ) &
        &         100._RK * real(tempVal, RK ) / &
        &         real (tempVal2, RK )
        endif

#else
        write( IOBuffer, '("Acceptance rate volume changes", T32, "in %:", F20.9)' ) &
&         100._RK * real( this%NResizeSuccesses, RK ) / &
&         real ( this%NResizeAttempts, RK )

#endif          

        call FileWrite( this%iounit_errors )
        
#if MPI_VER > 0
        call MPI_Reduce( this%DispVol,tempReal, 1, MPI_RK, MPI_MAX, &
        &     NRootProc, Communicator, ierror )
        if (Nproc == NRootProc) then
          write( IOBuffer, '("Maximum displacement volume", T33, "r`d:", F20.9)' ) &
        &    tempReal
        endif

#else
          write( IOBuffer, '("Maximum displacement volume", T33, "r`d:", F20.9)' ) &
        &    this%DispVol  

#endif  
     
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )
      end if

      do i = 1, this%NRealComponents
        pc => this%Component(i)

        ! Move and rotate acceptance rates
        write( IOBuffer, '("Component ", A)' ) pc%PotModFileName
        call FileWrite( this%iounit_errors )
#if MPI_VER > 0
        call MPI_Reduce( pc%NMoveSuccesses,tempVal, 1, MPI_INTEGER, MPI_SUM, &
        &     NRootProc, Communicator, ierror )
        call MPI_Reduce( pc%NMoveAttempts,tempVal2, 1, MPI_INTEGER, MPI_SUM, &
        &     NRootProc, Communicator, ierror )
        if (Nproc == NRootProc) then
          write( IOBuffer, '("Acceptance rate trans.", T32, "in %:", F20.9)' ) &
&         100._RK * real( tempVal, RK ) / real ( tempVal2, RK ) 
        endif

#else
          write( IOBuffer, '("Acceptance rate trans.", T32, "in %:", F20.9)' ) &
&         100._RK * real( pc%NMoveSuccesses, RK ) / real ( pc%NMoveAttempts, RK )   

#endif            

        call FileWrite( this%iounit_errors )
        if( pc%Molecule%IsElongated ) then
        
#if MPI_VER > 0
          call MPI_Reduce( pc%NRotateSuccesses,tempVal, 1, MPI_INTEGER, MPI_SUM, &
          &     NRootProc, Communicator, ierror )
          call MPI_Reduce( pc%NRotateAttempts,tempVal2, 1, MPI_INTEGER, MPI_SUM, &
          &     NRootProc, Communicator, ierror )
          if (Nproc == NRootProc) then
            write( IOBuffer, '(T17, "rotates", T32, "in %:", F20.9)' ) 100._RK &
&           * real( tempVal, RK ) / real (tempVal2, RK )
        endif

#else
            write( IOBuffer, '(T17, "rotates", T32, "in %:", F20.9)' ) 100._RK &
&           * real( pc%NRotateSuccesses, RK ) / real ( pc%NRotateAttempts, RK )

#endif         

          call FileWrite( this%iounit_errors )
        end if

        if( pc%ChemPotMethod .eq. ChemPotMethodGradIns ) then
          ! Biased move and rotate acceptance rates
#if MPI_VER > 0
          call MPI_Reduce( pc%NMoveBiasedSuccesses,tempVal, 1, MPI_INTEGER, MPI_SUM, &
          &     NRootProc, Communicator, ierror )
          call MPI_Reduce( pc%NMoveBiasedAttempts,tempVal2, 1, MPI_INTEGER, MPI_SUM, &
          &     NRootProc, Communicator, ierror )
          if (Nproc == NRootProc) then
            write( IOBuffer, '(T17, "biased trans.", T32, "in %:", F20.9)' ) &
&           100._RK * real(tempVal, RK ) / &
&           real ( tempVal2, RK )
          endif

#else
            write( IOBuffer, '(T17, "biased trans.", T32, "in %:", F20.9)' ) &
&           100._RK * real( pc%NMoveBiasedSuccesses, RK ) / &
&           real ( pc%NMoveBiasedAttempts, RK )

#endif              
     
          call FileWrite( this%iounit_errors )
          if( pc%Molecule%IsElongated ) then
          
#if MPI_VER > 0
          call MPI_Reduce( pc%NRotateBiasedSuccesses,tempVal, 1, MPI_INTEGER, MPI_SUM, &
          &     NRootProc, Communicator, ierror )
          call MPI_Reduce( pc%NRotateBiasedAttempts,tempVal2, 1, MPI_INTEGER, MPI_SUM, &
          &     NRootProc, Communicator, ierror )
          if (Nproc == NRootProc) then
            write( IOBuffer, '(T17, "biased rotates", T32, "in %:", F20.9)' ) &
&             100._RK * real(tempVal, RK ) / &
&             real (tempVal2, RK )
        endif

#else
            write( IOBuffer, '(T17, "biased rotates", T32, "in %:", F20.9)' ) &
&             100._RK * real( pc%NRotateBiasedSuccesses, RK ) / &
&             real ( pc%NRotateBiasedAttempts, RK )

#endif          
          
            call FileWrite( this%iounit_errors )
          end if
        end if

        ! Maximum translational and rotational displacements
#if MPI_VER > 0
        call MPI_Reduce( pc%DispTran,tempReal, 1, MPI_RK, MPI_MAX, &
        &     NRootProc, Communicator, ierror )
        if (Nproc == NRootProc) then
          write( IOBuffer, '("Maximum displacement trans.", T33, "r`d:", F20.9)' ) &
&         tempReal
        endif

#else
          write( IOBuffer, '("Maximum displacement trans.", T33, "r`d:", F20.9)' ) &
&          pc%DispTran

#endif  

        call FileWrite( this%iounit_errors )
        if( pc%Molecule%IsElongated ) then
#if MPI_VER > 0
        call MPI_Reduce( pc%DispRot,tempReal, 1, MPI_RK, MPI_MAX, &
        &     NRootProc, Communicator, ierror )
        if (Nproc == NRootProc) then
          write( IOBuffer, '(T22, "rotational", T33, "r`d:", F20.9)' ) &
&            tempReal
        endif

#else
          write( IOBuffer, '(T22, "rotational", T33, "r`d:", F20.9)' ) &
&            pc%DispRot

#endif  
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
#if MPI_VER > 0
          call MPI_Reduce( pc%NFluctUpSuccesses(:),tempVec1(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, &
          & MPI_SUM, NRootProc, Communicator, ierror )
          call MPI_Reduce( pc%NFluctUpAttempts(:),tempVec2(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, &
          & MPI_SUM, NRootProc, Communicator, ierror )
          call MPI_Reduce( pc%NFluctDownSuccesses(:),tempVec3(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, &
          & MPI_SUM, NRootProc, Communicator, ierror )
          call MPI_Reduce( pc%NFluctDownAttempts(:),tempVec4(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, &
          & MPI_SUM, NRootProc, Communicator, ierror )
          if (Nproc == NRootProc) then
            write(IOBuffer, '(2F10.4)') &
&             0._RK, &
&             real(tempVec3(pc%NFluctMax), RK) / &
&             real(tempVec4(pc%NFluctMax), RK) * 100._RK
            call FileWrite( this%iounit_errors )           
           
           do j = pc%NFluctMax -1, 1, -1
             write(IOBuffer, '(2F10.4)') &
&            real(tempVec1(j+1), RK) / &
&            real(tempVec2(j+1), RK) * 100._RK, &
&            real(tempVec3(j), RK) / &
&            real(tempVec4(j), RK) * 100._RK
             call FileWrite( this%iounit_errors )
           end do
           write(IOBuffer, '(2F10.4)') &
&             real(tempVec1(1), RK) / &
&             real(tempVec2(1), RK) * 100._RK, &
&             0._RK
           call FileWrite( this%iounit_errors )
           call FileWriteBlank( this%iounit_errors )
           call FileWriteBlank( this%iounit_errors )
          endif
#else
          write(IOBuffer, '(2F10.4)') &
&             0._RK, &
&             real(pc%NFluctDownSuccesses(pc%NFluctMax), RK) / &
&               real(pc%NFluctDownAttempts(pc%NFluctMax), RK) * 100._RK
            call FileWrite( this%iounit_errors )
          do j = pc%NFluctMax -1, 1, -1
            write(IOBuffer, '(2F10.4)') &
&             real(pc%NFluctUpSuccesses(j+1), RK) / &
&               real(pc%NFluctUpAttempts(j+1), RK) * 100._RK, &
&             real(pc%NFluctDownSuccesses(j), RK) / &
&               real(pc%NFluctDownAttempts(j), RK) * 100._RK
            call FileWrite( this%iounit_errors )
          end do
          write(IOBuffer, '(2F10.4)') &
&             real(pc%NFluctUpSuccesses(1), RK) / &
&               real(pc%NFluctUpAttempts(1), RK) * 100._RK, &
&             0._RK
            call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )

#endif             
        end if
      end do

      ! Inserts and deletes acceptance rates
      if( EnsembleType .eq. EnsembleTypeGE .or. &
&         EnsembleType .eq. EnsembleTypeHA ) then
#if MPI_VER > 0
        call MPI_Reduce( this%NInsertSuccesses,tempVal, 1, MPI_INTEGER, MPI_SUM, &
        &     NRootProc, Communicator, ierror )
        call MPI_Reduce( this%NInsertAttempts,tempVal2, 1, MPI_INTEGER, MPI_SUM, &
        &     NRootProc, Communicator, ierror )
        if (Nproc == NRootProc) then
          write( IOBuffer, '("Acceptance rate inserts", T32, "in %:", F20.9)' ) &
&         100._RK * real( tempVal, RK ) / real ( tempVal2, RK )
        endif

#else
        write( IOBuffer, '("Acceptance rate inserts", T32, "in %:", F20.9)' ) &
&         100._RK * real( this%NInsertSuccesses, RK ) / real ( this%NInsertAttempts, RK )

#endif 

        call FileWrite( this%iounit_errors )
#if MPI_VER > 0
        call MPI_Reduce( this%NDeleteSuccesses,tempVal, 1, MPI_INTEGER, MPI_SUM, &
        &     NRootProc, Communicator, ierror )
        call MPI_Reduce( this%NDeleteAttempts,tempVal2, 1, MPI_INTEGER, MPI_SUM, &
        &     NRootProc, Communicator, ierror )
        if (Nproc == NRootProc) then
          write( IOBuffer, '("Acceptance rate deletes", T32, "in %:", F20.9)' ) &
&         100._RK * real(tempVal, RK ) / real ( tempVal2, RK )
        endif

#else
          write( IOBuffer, '("Acceptance rate deletes", T32, "in %:", F20.9)' ) &
&         100._RK * real( this%NDeleteSuccesses, RK ) / real ( this%NDeleteAttempts, RK )

#endif         
        
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )
      end if
    end if
    call FileWriteBlank( this%iounit_errors )

    ! Close final result file
    call FileClose( this%iounit_errors )

  end subroutine TEnsemble_ErrorsUpdate



!==============================================================!
!  Subroutine TEnsemble_SVCOutput                              !
!==============================================================!

  subroutine TEnsemble_SVCOutput( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK) :: value
    integer  :: i, j

    ! Open final result file
    write( IOBuffer, '(I16)' ) this%EnsembleNumber
    call FileRewrite( this%iounit_errors, &
&     trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ErrorsFileExtension )

    write( IOBuffer, '(76("="))')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("*                           Publishing with ms2                            *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* Every user agrees to cite ms2 upon usage as follows                      *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* ------------------------------------------------------------------------ *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* S. Deublein, B. Eckl, J. Stoll, S. Lishchuk, G. Guevara-Carrion,         *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* C.W. Glass, T. Merker, M. Bernreuther, H. Hasse, J. Vrabec               *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* Computer Physics Communications (2011)                                   *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* DOI:10.1016/j.cpc.2011.04.026                                        *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(76("="))')
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Separator
    call FileWriteBlank( this%iounit_errors )
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
!#if MPI_VER > 0
!    include 'mpif.h'
!#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    type(TComponent), pointer :: pc
    integer                   :: i
#if TRANS ==1
    integer                   :: j
#endif

    if( SimulationType .eq. MonteCarlo ) then
      if( NProc /= NRootProc ) return
    endif 
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

#if TRANS ==1
    if( RootProc .and. ( CorrfunMode .eq. active ) ) then
    write( iounit_restart, '(I10)' ) this%Ncorr
    write( iounit_restart, '(I10)' ) this%Mmess

    do i = 1, 3*this%Npart
        do j = 1, this%Ncorr
            write( iounit_restart, '(ES20.12E3)' )  this%a( i, j)
        end do
    end do

    do i = 1, this%Ncorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vsk(i,:)
    end do
    do i = 1, this%Ncorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vsp(i,:)
    end do
    do i = 1, this%Ncorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vbk(i,:)
    end do
    do i = 1, this%Ncorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vbp(i,:)
    end do
    do i = 1, this%Ncorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vckt(i,:)
    end do
    do i = 1, this%Ncorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vckr(i,:)
    end do
    do i = 1, this%Ncorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vcpt(i,:)
    end do
    do i = 1, this%Ncorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vcpr(i,:)
    end do

    do i = 1, this%Ncorr
        write( iounit_restart, '(ES20.12E3)' ) this%cf_c(i)
    end do
    do i = 1, this%Ncorr
        write( iounit_restart, '(ES20.12E3)' ) this%cf_vb(i)
    end do
    do i = 1, this%Ncorr
        write( iounit_restart, '(ES20.12E3)' ) this%cf_vs(i)
    end do

!    if (this%Ncomponents==2) then
!        do i = 1, this%Ncorr
!            write( iounit_restart, '(ES20.12E3)' ) this%cf_db(i)
!        end do
!    end if

    do i = 1, this%Ncomponents
      do j = 1, this%Ncorr
        write( iounit_restart, '(ES20.12E3)' ) this%cf_d(i , j)
      end do
    end do

    if (this%Ncomponents>1) then
      do i = 1, this%Ncomponents*this%Ncomponents
        do j = 1, this%Ncorr
          write( iounit_restart, '(ES20.12E3)' ) this%lamda(i , j)
        end do
      end do
    end if

    write( iounit_restart, '(I10)' ) NBlocksMaxCF

    do i = 1, this%NComponents
      call RestartSaveCF( this%Sumself_i(i) )
    end do

    if(this%Ncomponents == 2) then
      call RestartSaveCF( this%SumBin_d )
    end if

    if(this%Ncomponents == 3) then
    call RestartSaveCF( this%SumTer_a )
    call RestartSaveCF( this%SumTer_b )
    call RestartSaveCF( this%SumTer_c )
    end if

    call RestartSaveCF( this%SumVisco_s )
    call RestartSaveCF( this%SumVisco_b )
    call RestartSaveCF( this%SumConduct )

    do i = 1,3
      write( iounit_restart, '(ES20.12E3)' )  this%sp(i)
    end do

    do i = 1,3
      write( iounit_restart, '(ES20.12E3)' )  this%sc(i)
    end do

endif
#endif
!#if MPI_VER > 0
! call MPI_Barrier( Communicator, ierror )
!#endif


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
    integer                   :: i,j,counter

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
&     Communicator, ierror )
    call MPI_Bcast( this%Volume0, 1, MPI_RK, NRootProc, &
&     Communicator, ierror )
    if( SimulationType .eq. MonteCarlo ) then
      call MPI_Bcast( this%DispVol, 1, MPI_RK, NRootProc, &
&       Communicator, ierror )
      call MPI_Bcast( this%NResizeAttempts, 1, MPI_INTEGER, NRootProc, &
&       Communicator, ierror )
      call MPI_Bcast( this%NResizeSuccesses, 1, MPI_INTEGER, NRootProc, &
&       Communicator, ierror )
      if( EnsembleType .eq. EnsembleTypeGE .or. &
&         EnsembleType .eq. EnsembleTypeHA ) then
        call MPI_Bcast( this%NInsertAttempts, 1, MPI_INTEGER, NRootProc, &
&         Communicator, ierror )
        call MPI_Bcast( this%NInsertSuccesses, 1, MPI_INTEGER, NRootProc, &
&         Communicator, ierror )
        call MPI_Bcast( this%NDeleteAttempts, 1, MPI_INTEGER, NRootProc, &
&         Communicator, ierror )
        call MPI_Bcast( this%NDeleteSuccesses, 1, MPI_INTEGER, NRootProc, &
&         Communicator, ierror )
      end if
    end if
    call MPI_Bcast( this%NRCutoffMax, 1, MPI_INTEGER, NRootProc, &
&     Communicator, ierror )
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
    counter = this%NRealComponents+1
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
        pc%NFluctState = 0
        do j = counter,counter + this%Component(i)%Molecule%NFluct-1
          if (this%Component(j)%NPart .eq. 1) then
            pc%NFluctState = j-counter+1
          end if
        end do
        counter = counter + this%Component(i)%Molecule%NFluct

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

#if TRANS ==1
if ( CorrfunMode .eq. active ) then
  if( RootProc ) then
      read( iounit_restart, '(I10)' ) this%Ncorr
      read( iounit_restart, '(I10)' ) this%Mmess
  end if

#if MPI_VER > 0
    call MPI_Bcast( this%Ncorr, 1, MPI_INTEGER, NRootProc, &
&       Communicator, ierror )
    call MPI_Bcast( this%Mmess, 1, MPI_INTEGER, NRootProc, &
&     Communicator, ierror )
#endif

  if( RootProc ) then
    do i = 1, 3*this%Npart
        do j = 1, this%Ncorr
            read( iounit_restart, '(ES20.12E3)' )  this%a( i, j)
        end do
    end do

    do i = 1, this%Ncorr
        read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vsk(i,:)
    end do
    do i = 1, this%Ncorr
        read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vsp(i,:)
    end do
    do i = 1, this%Ncorr
        read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vbk(i,:)
    end do
    do i = 1, this%Ncorr
        read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vbp(i,:)
    end do
    do i = 1, this%Ncorr
        read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vckt(i,:)
    end do
    do i = 1, this%Ncorr
        read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vckr(i,:)
    end do
    do i = 1, this%Ncorr
        read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vcpt(i,:)
    end do
    do i = 1, this%Ncorr
        read( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vcpr(i,:)
    end do

    do i = 1, this%Ncorr
        read( iounit_restart, '(ES20.12E3)' ) this%cf_c(i)
    end do
    do i = 1, this%Ncorr
        read( iounit_restart, '(ES20.12E3)' ) this%cf_vb(i)
    end do
    do i = 1, this%Ncorr
        read( iounit_restart, '(ES20.12E3)' ) this%cf_vs(i)
    end do

!     if (this%Ncomponents==2) then
!         do i = 1, this%Ncorr
!             read( iounit_restart, '(ES20.12E3)' ) this%cf_db(i)
!         end do
!     end if

    do i = 1, this%Ncomponents
        do j = 1, this%Ncorr
            read( iounit_restart, '(ES20.12E3)' ) this%cf_d(i , j)
        end do
    end do

    if (this%Ncomponents>1) then
      do i = 1, this%Ncomponents*this%Ncomponents
        do j = 1, this%Ncorr
            read( iounit_restart, '(ES20.12E3)' ) this%lamda(i , j)
        end do
      end do
    end if

    read( iounit_restart, '(I10)' ) NBlocksMaxCF

    do i = 1, this%NComponents
      call RestartReadCF( this%Sumself_i(i) )
    end do

    if(this%Ncomponents == 2) then
      call RestartReadCF( this%SumBin_d )
    end if

    if(this%Ncomponents == 3) then
      call RestartReadCF( this%SumTer_a )
      call RestartReadCF( this%SumTer_b )
      call RestartReadCF( this%SumTer_c )
    end if

    call RestartReadCF( this%SumVisco_s )
    call RestartReadCF( this%SumVisco_b )
    call RestartReadCF( this%SumConduct )

    do i = 1,3
        read( iounit_restart, '(ES20.12E3)' )  this%sp(i)
    end do

    do i = 1,3
        read( iounit_restart, '(ES20.12E3)' )  this%sc(i)
    end do

  end if
end if
#endif


    ! Update density and box length
    call UpdateBoxLength( this )

    ! Update mole fractions
    call UpdateFractions( this )

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


! TRANSPORT ab hier
!==============================================================!
!  Subroutine TEnsemble_CalCorrFun                             !
!==============================================================!

    subroutine TEnsemble_CalCorrFun( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this
#if TRANS==1
    ! Declare local variables
    integer  :: nmess, i, j, j0, k, s
    integer  :: CFindex, Mindex
    integer  :: Npart3,Npart2, Ncomp2
    integer  :: np, nc
    real(RK) :: sx(this%Ncomponents, this%Ncorr ), sy(this%Ncomponents, this%Ncorr )
    real(RK) :: sz(this%Ncomponents, this%Ncorr )
    real(RK) :: EKinTran(this%Ncomponents,this%Npart)
    real(RK) :: EKinRot(this%Ncomponents,this%Npart)
    real(RK) :: BoxLength_dt,BoxLength_dt2
    real(RK) :: tempf(3), virf(3)
    real(RK), pointer :: pfb(:,:), pfs(:,:), pftc(:,:), pfrc(:,:)

    Npart3 = 3*this%Npart
    Npart2 = 2*this%Npart
    Ncomp2 = this%NComponents*this%NComponents
    BoxLength_dt       =  this%BoxLength/TimeStep
    BoxLength_dt2      =  BoxLength_dt**2

!Evaluate FTc and FRC components (parallel version)

      do i = 1, this%NComponents
        call ForceTransport( this%Component(i) )
#if MPI_VER > 0
        call MPI_Bcast( this%Component(i)%P1(:, :), size( this%Component(i)%P1 ), MPI_RK, &
&       NRootProc, Communicator, ierror )
#endif
      end do

! Calculate sum of terms of the pressure tensor (kinetic and potential)

    do i = 1, this%Ncomponents
       np = this%Component(i)%Npart
#if MPI_VER > 0
       pfb => this%Component(i)%fbAll(:,:)
#else
       pfb => this%Component(i)%fb(:,:)
#endif
       do k =1, 3
          do j = 1, np

             this%sc(k) = this%sc(k) + this%Component(i)%KinETran(j,k)
             this%sp(k) = this%sp(k) + pfb(j, k)
          end do
       end do
     end do

     tempf(:)  = this%sc(:)/Step
     virf(:)   = this%sp(:)/Step

! Calculate the kinetic Energy/particle and total

     EKinTran(:,:) = 0._RK
     EKinRot (:,:) = 0._RK

   do i = 1, this%Ncomponents
     np = this%Component(i)%Npart

     do k = 1, 3
          do j = 1, np
             EkinTran(i,j) = EkinTran(i,j)+ this%Component(i)%KinETran(j,k)*0.5d0

             if( this%Component(i)%Molecule%IsElongated ) then

               EKinRot(i,j)= EKinRot(i,j) + (this%Component(i)%W0(j,k)*this%Component(i)%W0(j,k)*    &
  &                                          this%Component(i)%Molecule%MOI(k))*0.5d0
             else
              EKinRot(i,j)= 0._RK
            end if

          end do
       end do
     end do


!Caculate matrix indexes

    if (Step.le.this%Ncorr) then
        Mindex = Step
        else
           if(mod(Step, this%Ncorr).eq.0)then
              Mindex = this%Ncorr
              else
                 Mindex = mod(Step, this%NCorr)
           end if
    end if


!Write transport properties Matrixes (root Processor)

    do k = 1, 3
       this%vsk(Mindex,  k) = 0._RK
       this%vsp(Mindex,  k) = 0._RK
       this%vbk(Mindex,  k) = 0._RK
       this%vbp(Mindex,  k) = 0._RK
       this%vckt(Mindex, k) = 0._RK
       this%vckr(Mindex, k) = 0._RK
       this%vcpt(Mindex, k) = 0._RK
       this%vcpr(Mindex, k) = 0._RK
    end do


! part calculated together with force

    do i = 1, this%Ncomponents
       np = this%Component(i)%Npart
#if MPI_VER > 0
       pfb => this%Component(i)%fbAll(:,:)
       pfs => this%Component(i)%fsAll(:,:)
#else
       pfb => this%Component(i)%fb(:,:)
       pfs => this%Component(i)%fs(:,:)
#endif
       pftc => this%Component(i)%ftc(:,:)
       pfrc => this%Component(i)%frc(:,:)

       do k = 1,3
          do j = 1, np
             this%vsp(Mindex, k) = this%vsp(Mindex, k) + pfs(j, k)
             this%vbp(Mindex, k) = this%vbp(Mindex, k) + pfb(j, k)
             this%vcpr(Mindex, k) = this%vcpr(Mindex, k) + pfrc(j, k)
             this%vcpt(Mindex, k) = this%vcpt(Mindex, k) + pftc(j, k)
          end do
       end do
     end do


! kinetic part

   !Diffusion matrix a

     do i = 1, this%Ncomponents
       if (i.eq.1) then
         j0 = 0
         k  = this%Component(i)%Npart
       else
         j0 = k
         k  = k + this%Component(i)%Npart
       end if
         do j = 1, this%Component(i)%Npart

         this%a(j+j0           , Mindex) = this%Component(i)%P1(j,1)*BoxLength_dt
         this%a(j+j0+this%Npart, Mindex) = this%Component(i)%P1(j,2)*BoxLength_dt
         this%a(j+j0+Npart2    , Mindex) = this%Component(i)%P1(j,3)*BoxLength_dt

         end do
     end do


   !parts of the stress and energy tensors

    ! shear off diagonal terms
     do i = 1, this%Ncomponents
        np = this%Component(i)%Npart

          do j = 1, np

         this%vsk(Mindex, 1) = this%vsk(Mindex, 1) + this%Component(i)%P1(j,1)*        &
      &                        this%Component(i)%P1(j,2)*this%Component(i)%Molecule%Mass*  &
      &                        BoxLength_dt2

         this%vsk(Mindex, 2) = this%vsk(Mindex, 2) + this%Component(i)%P1(j,1)*        &
      &                        this%Component(i)%P1(j,3)*this%Component(i)%Molecule%Mass*  &
      &                        BoxLength_dt2

         this%vsk(Mindex, 3) = this%vsk(Mindex, 3) + this%Component(i)%P1(j,2)*        &
      &                        this%Component(i)%P1(j,3)*this%Component(i)%Molecule%Mass*  &
      &                        BoxLength_dt2


       end do
     end do

  !bulk diagonal terms and energy tensor kinetic part

    do i = 1, this%Ncomponents
       np = this%Component(i)%Npart
       do k = 1, 3
          do j = 1, np

              this%vbk(Mindex, k) = this%vbk(Mindex, k) + this%Component(i)%KinETran(j,k)

              this%vckt(Mindex, k)= this%vckt(Mindex, k) + this%Component(i)%P1(j, k)*EKinTran(i,j)*  &
     &                              this%Component(i)%Molecule%Mass*BoxLength_dt

              this%vckr(Mindex, k)= this%vckr(Mindex, k) + this%Component(i)%P1(j, k)*EKinRot(i,j)*  &
     &                              this%Component(i)%Molecule%Mass*BoxLength_dt
         end do
      end do
    end do




! Calculate Auto Correlation Functions

 if((Step.gt.this%Ncorr).and.(mod(Step, this%NSpancf).eq.0))then


      if (Mindex .eq. this%Ncorr) then
         CFindex = 1                       !index of t = t0
         else
            CFindex = Mindex +1
      end if

      if(this%Ncomponents .gt. 1)then
        do nmess =1, this%Ncorr

           sx(:,nmess) = 0._RK
           sy(:,nmess) = 0._RK
           sz(:,nmess) = 0._RK

           do i = 1, this%Ncomponents
              if (i.eq.1) then
                 j0 = 0
                 k  = this%Component(i)%Npart
                 else
                  j0 = k
                  k  = k + this%Component(i)%Npart
              end if
              do j = 1, this%Component(i)%Npart
                 sx(i,nmess)  = sx(i, nmess) + this%a(j+j0            , nmess)
                 sy(i,nmess)  = sy(i, nmess) + this%a(j+j0+this%Npart , nmess)
                 sz(i,nmess)  = sz(i, nmess) + this%a(j+j0+Npart2     , nmess)
              end do
           end do
        end do
      end if


! Write the diffusion matrix

      k      = 0
      do nmess = 1, this%Ncorr
         if (nmess .le. (this%Ncorr-Mindex)) then
            s = Mindex + nmess                              ! index t= t
            else
               if (Mindex .eq. this%Ncorr) then
                  s = nmess
                  else
                  s = mod ((nmess + Mindex), this%Ncorr)
               end if
         end if

!Calculate auto-correlation functions
         !Nur selbstdiffusion
         do i = 1, this%Ncomponents
            if (i.eq.1) then
               j0 = 0
               k  = this%Component(i)%Npart
            else
                  j0 = k
                  k  = k + this%Component(i)%Npart
            end if
            do j = 1, this%Component(i)%Npart

               this%cf_d(i , nmess) = this%cf_d(i, nmess)                                       &
                                    + this%a(j+j0            ,CFindex)*this%a(j+j0            ,s) &
                                    + this%a(j+j0+this%Npart ,CFindex)*this%a(j+j0+this%Npart ,s) &
                                    + this%a(j+j0+Npart2     ,Cfindex)*this%a(j+j0+Npart2     ,s)
            end do
         end do
      end do

      do nmess= 1, this%Ncorr

          if (nmess .le. (this%Ncorr-Mindex)) then
             s = Mindex + nmess
            else
              if (Mindex .eq. this%Ncorr) then
                 s = nmess
                else
                  s = mod ((nmess + Mindex), this%Ncorr)
               end if
            end if

           if (this%Ncomponents .gt. 1)then
               nc = this%NComponents
               k = 1
               do i = 1, nc
                  do j = 1,nc
                  this%lamda(k, nmess) = this%lamda(k, nmess) + sx(i, CFindex)*sx(j, s) &
                                                              + sy(i, CFindex)*sy(j, s) &
                                                              + sz(i, CFindex)*sz(j, s)
                  k = k + 1
                  end do
               end do
            end if


  !         if(this%Ncomponents==2)then

  !         this%cf_db(nmess) = this%cf_db(nmess) + sx(1, CFindex)*sx(1, s) &
   !                                              + sy(1, CFindex)*sy(1, s) &
   !                                              + sz(1, CFindex)*sz(1, s)

    !       end if

           do k = 1, 3

             this%cf_vs(nmess) = this%cf_vs(nmess) + this%vsk(CFindex, k)*this%vsk(s, k) + &
&                                                    this%vsp(CFindex, k)*this%vsp(s, k) + &
&                                                    this%vsk(CFindex, k)*this%vsp(s, k) + &
&                                                    this%vsp(CFindex, k)*this%vsk(s, k)
             do j = 1, 3 ! FIXME: PROBLEM mit gemischten Termen j index zuviel (vermutlich)

                this%cf_vb(nmess) =   this%cf_vb(nmess) + &
&                                   ( this%vbk(CFindex, j)-tempf(j))*(this%vbk(s, k)-tempf(k)) + &
&                                   ( this%vbp(CFindex, j)-virf(j)) *(this%vbp(s, k)-virf(k) ) + &
&                                   ( this%vbk(CFindex, j)-tempf(j))*(this%vbp(s, k)-virf(k) ) + &
&                                   ( this%vbp(CFindex, j)-virf(j))*(this%vbk(s, k)-tempf(k) )
             end do

             this%cf_c(nmess) =  this%cf_c(nmess) + this%vckt(CFindex, k)*this%vckt(s, k) + &
&                                                   this%vckr(CFindex, k)*this%vckr(s, k) + &
&                                                   this%vcpt(CFindex, k)*this%vcpt(s, k) + &
&                                                   this%vcpr(CFindex, k)*this%vcpr(s, k) + &
&                                                   this%vckt(CFindex, k)*this%vcpt(s, k) + &
&                                                   this%vcpt(CFindex, k)*this%vckt(s, k) + &
&                                                   this%vckr(CFindex, k)*this%vcpr(s, k) + &
&                                                   this%vcpr(CFindex, k)*this%vckr(s, k) + &
&                                                   this%vckr(CFindex, k)*this%vckt(s, k) + &
&                                                   this%vckt(CFindex, k)*this%vckr(s, k) + &
&                                                   this%vcpt(CFindex, k)*this%vcpr(s, k) + &
&                                                   this%vcpr(CFindex, k)*this%vcpt(s, k) + &
&                                                   this%vckt(CFindex, k)*this%vcpr(s, k) + &
&                                                   this%vcpr(CFindex, k)*this%vckt(s, k) + &
&                                                   this%vckr(CFindex, k)*this%vcpt(s, k) + &
&                                                   this%vcpt(CFindex, k)*this%vckr(s, k)
       end do

      end do
    this%Mmess  = this%Mmess +1
    end if
#endif
   end subroutine TEnsemble_CalCorrFun



!==============================================================!
!  Subroutine TEnsemble_IntCorrFun                             !
!==============================================================!

  subroutine TEnsemble_IntCorrFun( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this
#if TRANS==1
    ! Declare local varibles
    integer  :: i, k
    integer  :: ncomp2
    real(RK) :: helpvar, det, deter1, deter2, deter3, deter4
    real(RK) :: x1, x2, x3
    real(RK) :: Inv_x1, Inv_x2, Inv_x3
    real(RK) :: B11, B12, B21, B22

    ncomp2 = this%NComponents*this%NComponents

    do i  = 1, this%NComponents

      helpvar =  1._RK /(3._RK *this%Component(i)%Npart * this%Mmess)

      this%sinte_i(i,:) = simpson( this%cf_d(i,:)/this%cf_d(i, 1), TimeStep, this%NCorr )
      this%selfd_i(i) = this%sinte_i(i, this%NCorr) * this%cf_d(i, 1) * helpvar

    end do


    if ( this%NComponents .gt. 1) then

      helpvar =  1._RK /(3._RK *this%Npart * this%Mmess)
      do k = 1, ncomp2
        this%sinte_lamda(k, :) = simpson(this%lamda(k,:)/this%lamda(k,1), TimeStep, this%NCorr)
      end do

      if( this%NComponents == 2 ) then

 !       this%sinte_db = simpson( this%cf_db(:)/this%cf_db(1), TimeStep, this%NCorr )

        this%binary_d = (((this%sinte_lamda(1,this%NCorr)*this%lamda(1,1)) * &
&                         (this%Component(2)%Fraction/this%Component(1)%Fraction)) + &
&                        ((this%sinte_lamda(4,this%NCorr)*this%lamda(4,1)) * &
&                         (this%Component(1)%Fraction/this%Component(2)%Fraction)) - &
&                         (this%sinte_lamda(2,this%NCorr)*this%lamda(2,1)) - &
&                         (this%sinte_lamda(3,this%NCorr)*this%lamda(3,1)))* helpvar

!      this%binary_d = this%sinte_db( this%NCorr ) * this%cf_db(1)* (1._RK /(3._RK *this%Component(1)%Npart * this%Mmess) * &
!&       ( (this%Component(1)%Fraction * this%Component(1)%Molecule%Mass + &
!&          this%Component(2)%Fraction * this%Component(2)%Molecule%Mass) / &
!&         (this%Component(2)%Fraction * this%Component(2)%Molecule%Mass) ) **2 * &
!&          this%Component(2)%Fraction
      end if

     if( this%NComponents == 3 ) then

        x1 = this%Component(1)%Fraction
        x2 = this%Component(2)%Fraction
        x3 = this%Component(3)%Fraction
        Inv_x1 = 1._RK / x1
        Inv_x2 = 1._RK / x2
        Inv_x3 = 1._RK / x3

      deter1 = (((1._RK - x1)*(((this%sinte_lamda(1,this%NCorr)*this%lamda(1,1))*Inv_x1) - &
&                              ((this%sinte_lamda(3,this%NCorr)*this%lamda(3,1))*Inv_x3)))- &
&                       ((x1)*(((this%sinte_lamda(4,this%NCorr)*this%lamda(4,1))*Inv_x1) - &
&                              ((this%sinte_lamda(6,this%NCorr)*this%lamda(6,1))*Inv_x3) + &
&                              ((this%sinte_lamda(7,this%NCorr)*this%lamda(7,1))*Inv_x1) - &
&                              ((this%sinte_lamda(9,this%NCorr)*this%lamda(9,1))*Inv_x3))))*helpvar

       deter2 = (((1._RK - x1)*(((this%sinte_lamda(2,this%NCorr)*this%lamda(2,1))*Inv_x2) - &
&                               ((this%sinte_lamda(3,this%NCorr)*this%lamda(3,1))*Inv_x3)))- &
&                        ((x1)*(((this%sinte_lamda(5,this%NCorr)*this%lamda(5,1))*Inv_x2) - &
&                               ((this%sinte_lamda(6,this%NCorr)*this%lamda(6,1))*Inv_x3) + &
&                               ((this%sinte_lamda(8,this%NCorr)*this%lamda(8,1))*Inv_x2) - &
&                               ((this%sinte_lamda(9,this%NCorr)*this%lamda(9,1))*Inv_x3))))*helpvar

       deter3 = (((1._RK - x2)*(((this%sinte_lamda(4,this%NCorr)*this%lamda(4,1))*Inv_x1) - &
&                               ((this%sinte_lamda(6,this%NCorr)*this%lamda(6,1))*Inv_x3)))- &
&                        ((x2)*(((this%sinte_lamda(1,this%NCorr)*this%lamda(1,1))*Inv_x1) - &
&                               ((this%sinte_lamda(3,this%NCorr)*this%lamda(3,1))*Inv_x3) + &
&                               ((this%sinte_lamda(7,this%NCorr)*this%lamda(7,1))*Inv_x1) - &
&                               ((this%sinte_lamda(9,this%NCorr)*this%lamda(9,1))*Inv_x3))))*helpvar

       deter4 = (((1._RK - x2)*(((this%sinte_lamda(5,this%NCorr)*this%lamda(5,1))*Inv_x2) - &
&                               ((this%sinte_lamda(6,this%NCorr)*this%lamda(6,1))*Inv_x3)))- &
&                        ((x2)*(((this%sinte_lamda(2,this%NCorr)*this%lamda(2,1))*Inv_x2) - &
&                               ((this%sinte_lamda(3,this%NCorr)*this%lamda(3,1))*Inv_x3) + &
&                               ((this%sinte_lamda(8,this%NCorr)*this%lamda(8,1))*Inv_x2) - &
&                               ((this%sinte_lamda(9,this%NCorr)*this%lamda(9,1))*Inv_x3))))*helpvar

      !obtain matrix [B] so that [B]=[D]-1

     ! determinat of matrix [B]

      det = (deter1*deter4)-(deter2*deter3)

      B11 =  deter4 * (1._RK/det)
      B12 = -deter2 * (1._RK/det)
      B21 = -deter3 * (1._RK/det)
      B22 =  deter1 * (1._RK/det)

      this%ternary_a =  1._RK  / ( (B11) + ( x2* B12 * Inv_x1) )
      this%ternary_b =  1._RK  / ( (B11) - ( (x1 + x3) * B12 *Inv_x1))
      this%ternary_c =  1._RK  / ( (B22) + ( x1* B21 * Inv_x2))

     end if

    end if


    helpvar =  this%Density /(3._RK *this%Npart * this%Mmess * this%Temperature)

    this%sinte_vs = simpson( this%cf_vs(:)/this%cf_vs(1), TimeStep, this%NCorr )
    this%visco_s = this%sinte_vs( this%NCorr ) * this%cf_vs(1) * helpvar

    this%sinte_vb = simpson( this%cf_vb(:)/this%cf_vb(1), TimeStep, this%NCorr )
    this%visco_b = this%sinte_vb( this%NCorr ) * this%cf_vb(1) * (helpvar / 3._RK)

    this%sinte_c = simpson( this%cf_c(:)/this%cf_c(1), TimeStep, this%NCorr )
    this%conduct = this%sinte_c( this%NCorr ) * this%cf_c(1) * (helpvar / this%Temperature)


  contains



    function trapezoid(values, step, n) result(integral)

      ! Declare arguments
      real(RK), intent(in) :: values(:), step
      integer, intent(in)  :: n

      ! Declare result
      real(RK) :: integral(n)

      ! Declare local variables
      integer :: i

      ! Return if no values to integrate
      if( n < 1 ) return

      ! Initialize result
      integral = 0._RK

      ! Calculate integral via trapezoid method
      do i = 2, n
        integral(i) = integral(i-1) + values(i) + values(i-1)
      end do
      integral = integral * .5_RK * step

    end function



    function simpson(values, step, n) result(integral)

      ! Declare arguments
      real(RK), intent(in) :: values(:), step
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
        integral(i) = integral(i-2) + values(i) + 4._RK * values(i-1) + values(i-2)
        integral(i-1) = .5 * (integral(i) + integral(i-2))
      end do
      integral = integral * step / 3._RK

      if( mod(n, 2) == 0 ) integral(n) = integral(n-1)

    end function


#endif
  end subroutine TEnsemble_IntCorrFun

!==============================================================!
!  Subroutine TEnsemble_DbgWrite                               !
!==============================================================!

  subroutine TEnsemble_DbgWrite( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare external procedures
    integer, external :: system

    ! Declare local variables
    type(TComponent), pointer :: pc
    type(TSiteLJ126), pointer :: plj
    integer                   :: nc, i, i1, i2, j, n, n2, n3
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



end module ms2_ensemble

