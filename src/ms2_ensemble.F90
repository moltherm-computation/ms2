!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 2.0                !
!  (c) 2014 by TU Kaiserslautern                               !
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


#ifndef OSMOP
#define OSMOP 0
#endif

#ifndef HBOND
#define HBOND 0
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

    ! I/O unit for RDF file
    integer :: iounit_rdf

    ! I/O unit for ThermoInt File
    integer :: iounit_thermoint

    ! I/O unit for result ACF
    integer :: iounit_rescf

    ! I/O unit for visualization H-bonding file
    integer :: iounit_visualHB

    ! I/O unit for Profile file
    integer :: iounit_dcp

#if  TRANS == 1
    logical :: Conductivity   !TRANSPORT_thisline
    logical :: EConductivity
    logical :: MolarEnthConduct
    logical :: Bulkviscosity
#endif

    ! Maximum number of particles
    integer, pointer :: NPartMax
    integer, pointer :: NPartMaxFluct

    ! Maximum number of units
    integer, pointer :: NUnitMax

    ! Number of particles in ensemble
    integer :: NPart, NPartInitial
    integer :: NPartLBound, NPartUBound

    ! Maximum number of test particles
    integer :: NTestMax

    ! Maximum number of fluctuating states
    integer :: NFluctMax

    ! Number of degrees of freedom
    integer :: NDFTran, NDFRot, NDF, constrNDF

    ! Mass of piston
    real(RK) :: PistonMass

    ! Optional calculation of pressure
    logical :: OptPressure

    ! Positions and orientations of test particles
    real(RK), pointer, contiguous :: P0Test(:, :, :), Q0Test(:, :, :)

    ! Number of unit cells in one dimension of lattice
    integer :: NCells

    ! Number of components in ensemble
    integer :: NComponents, NRealComponents, NGradInsComp

    ! Maximum numbers of sites in components
    integer :: NLJ126Max, NChargeMax, NDipoleMax, NQuadrupoleMax

    ! Total number of Units
    integer :: NUnitTotal

    ! Components
    type(TComponent), pointer, contiguous :: Component(:)

    ! Interactions
    type(TInteraction), pointer :: Interaction(:, :)

    ! Initial values of temperature, pressure, density, hamiltonian and enthalpy
    real(RK) :: RefTemperature, RefPressure, RefDensity, RefHamiltonian, RefEnthalpy

    ! Values of density, enthalpy, betaT, dHdp and their uncertainties
    ! in corresponding liquid simulation (for GE ensemble only)
    real(RK) :: LiqDensity, VarLiqDensity, LiqEnthalpy, VarLiqEnthalpy
    real(RK) :: LiqBetaT, VarLiqBetaT, LiqdHdP, VarLiqdHdP

    ! Current values of temperature, pressure, density
    real(RK) :: Temperature, Pressure, Density

    ! Velocity scaling factor for temperature control
    real(RK) :: scale

#if OSMOP > 0
    real(RK) :: OsmoticPressure
#if OSMOP == 2
    real(RK), pointer, contiguous :: VirialProfile(:)
    real(RK), pointer, contiguous :: PressureProfile(:)
#endif
#endif

    ! Virial
    real(RK) :: Virial

    ! Inter und Intra(Bond) Virial
    real(RK) :: VirialInter
    real(RK) :: VirialIntra
    real(RK) :: VirialShake

    ! Scale coefficients for LJ126 epsilon and sigma
    real(RK), pointer, contiguous :: ScaleEpsilon(:, :), ScaleSigma(:, :)

    ! Cutoff radii
    real(RK) :: RCutoffLJ126LJ126
    real(RK) :: RCutoffDipoleDipole
    real(RK) :: RCutoffDipoleQuadrupole
    real(RK) :: RCutoffQuadrupoleQuadrupole

    !RDF Hilfsvariable
    real(RK) :: RDFdr
    real(RK), pointer, contiguous :: RDFVSchale(:)
    real(RK), pointer, contiguous :: RDFValue(:)

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
    integer :: NTransferAttempts, NTransferSuccesses

    ! Kinetic energy
    real(RK) :: EKin, EKinTran, EKinRot

    ! Potential energy
    real(RK) :: EPot

	! d2EpotdV2
    real(RK) :: d2EpotdV2

    ! Intra and inter potential energy
    real(RK) :: EPotIntra
    real(RK) :: EPotInter

    ! Different Intra energies: bond, angle, dihedral, nonbonded
    real(RK) :: EPotIntra_Bond
    real(RK) :: EPotIntra_Angle
    real(RK) :: EPotIntra_Dihedral
    real(RK) :: EPotIntra_Nonbonded

    ! Potential energy of test particles
    real(RK), pointer, contiguous :: EPotTest(:)

    ! Long-range corrections
    real(RK) :: EPotCorrLJ, VirialCorrLJ, d2EpotdV2CorrLJ
    real(RK) :: EPotCorrRF
    real(RK) :: EPotCorrRFPart, EPotCorrRFVol
    real(RK) :: VirialCorrRF

    ! Accumulated sums, averages and errors
    ! 1.) Basic sums
    type(TAccumulator) :: SumPressure
    type(TAccumulator) :: SumDensity
    type(TAccumulator) :: SumTemperature
    type(TAccumulator) :: SumEPot
    type(TAccumulator) :: SumEnthalpy
    type(TAccumulator) :: SumConfEnthalpy
    type(TAccumulator) :: SumEPotIntra
    type(TAccumulator) :: SumEPotInter
    type(TAccumulator) :: SumEPotIntra_Bond
    type(TAccumulator) :: SumEPotIntra_Angle
    type(TAccumulator) :: SumEPotIntra_Dihedral
    type(TAccumulator) :: SumEPotIntra_Nonbonded
    type(TAccumulator) :: SumVolume
    type(TAccumulator) :: SumVirial
    type(TAccumulator) :: SumVirialIntra
    type(TAccumulator) :: SumVirialInter
#if OSMOP > 0
    type(TAccumulator) :: SumOsmoticPressure
#if OSMOP == 2
    type(TAccumulator),pointer, contiguous :: SumPressureProfile(:)
#endif
#endif
    type(TAccumulator) :: SumNPart
    type(TAccumulator) :: SumdEpotdV
    type(TAccumulator) :: Sumd2EpotdV2

    !if( this%EnsembleType .eq. this%EnsembleTypeNVE .and. this%LongRange .eq. this%Rfield) then
      type(TAccumulator) :: SumHmU
      type(TAccumulator) :: SumHmUm1
      type(TAccumulator) :: SumHmUm2
      type(TAccumulator) :: SumHmUm3
      type(TAccumulator) :: SumHmUm1dUdV
      type(TAccumulator) :: SumHmUm1dUdV2
      type(TAccumulator) :: SumHmUm1d2UdV2
      type(TAccumulator) :: SumHmUm2dUdV
      type(TAccumulator) :: SumHmUm2dUdV2
      type(TAccumulator) :: SumHmUm2d2UdV2
      type(TAccumulator) :: SumHmUm3dUdV
      type(TAccumulator) :: SumHmUm3dUdV2
    !end if

    ! 2.) Combined sums
    type(TAccumulator) :: SumEPotSquared
    type(TAccumulator) :: SumEPotV
    type(TAccumulator) :: SumEPotVirial
    type(TAccumulator) :: SumEnthalpySquared
    type(TAccumulator) :: SumEnthalpyV
    type(TAccumulator) :: SumVolumeSquared
    type(TAccumulator) :: SumEPotCubic
    type(TAccumulator) :: SumdEpotdVSquared
    type(TAccumulator) :: SumEPotdEpotdV
    type(TAccumulator) :: SumEPotSquareddEpotdV
    type(TAccumulator) :: SumEPotdEpotdVSquared
    type(TAccumulator) :: SumEPotd2EpotdV2
    !if( (EnsembleType .eq. EnsembleTypeNVT .or. EnsembleType .eq. EnsembleTypeNVE) .and. LongRange .eq. Rfield) then
      type(TAccumulator) :: SumA10resI
      type(TAccumulator) :: SumA01resI
      type(TAccumulator) :: SumA20resI
      type(TAccumulator) :: SumA11resI
      type(TAccumulator) :: SumA02resI
      type(TAccumulator) :: SumA30resI
      type(TAccumulator) :: SumA21resI
      type(TAccumulator) :: SumA12resI
      type(TAccumulator) :: SumA10resII
      type(TAccumulator) :: SumA01resII
      type(TAccumulator) :: SumA20resII
      type(TAccumulator) :: SumA11resII
      type(TAccumulator) :: SumA02resII
      type(TAccumulator) :: SumA30resII
      type(TAccumulator) :: SumA21resII
      type(TAccumulator) :: SumA12resII
    !end if

    ! 3.) Derived sums
    type(TAccumulator) :: SumBetaT
    type(TAccumulator) :: SumdHdP
    type(TAccumulator) :: SumdUdV
    type(TAccumulator) :: SumCV
    type(TAccumulator) :: SumCP
    type(TAccumulator) :: SumAlphaP

    ! BiasedPartners
    integer            :: NGradIns
    integer,pointer, contiguous    :: BiasedPartners(:)

    ! Ewald Parameters
    real(RK),pointer, contiguous :: Ewald_Prefac(:)
    real(RK),pointer, contiguous :: Ewald_Vec(:,:)
    integer          :: NVecMax, NSQMAX, NMAX
    integer          :: BoxenAnzahlMax
    real(RK)         :: Kappa, KappaL
    real(RK)         :: USelbstterm
    real(RK)         :: UIntra
    real(RK)         :: UFourier, EVirial
    integer,pointer ::  NBox0,NBox1,NBox2
    real(RK),pointer, contiguous :: U_fourierLocal(:)
    real(RK),pointer, contiguous :: SSin(:),SCos(:)
    real(RK),pointer, contiguous :: rold(:,:)
    real(RK),pointer, contiguous :: Vec2(:)
    real(RK),pointer, contiguous :: VirIntra(:)

#if SPME > 0
    ! SPME parameters
    integer         :: splineorder
    integer         :: gridx
    integer         :: gridy
    integer         :: gridz
    real(RK),pointer, contiguous::   qgrida(:,:), qgrida_old(:,:), qgridb(:,:)
    integer*8       :: qgrid_forward, qgrid_backward
    real(RK),pointer, contiguous:: bbtot(:)
    real(RK),pointer, contiguous:: bsp_arr(:), bsp_modx(:), bsp_mody(:), bsp_modz(:)
    real(RK)        :: EVirialIntra
    real(RK),pointer, contiguous:: EPotPME(:), VirialPME(:), mm2(:)
#endif

    ! Extended ReactionField Method
    real(RK)        :: DebyeLen
    
    ! Residence Time
    integer         :: ResidPairs
    integer         :: ResidComp1, ResidSite1
    integer         :: ResidComp2, ResidSite2
    integer         :: ResidCem
    integer,pointer, contiguous :: CompPair(:,:), CompPair_Old(:,:)
    integer,pointer, contiguous :: ResidTimesStart(:), ResidTimesStart_Old(:), ResidPairsCem(:,:)
    integer         :: ResidPeriod
    integer         :: ResidBreak
    real(RK)        :: ResidLength
    real(RK)        :: ResidenceDuration
    logical         :: ResidenceTime
    type(TAccumulator) :: SumResidenceDuration
    type(TAccumulator) :: SumResidencePairs

#if  TRANS == 1
!TRANSPORT_start
    ! Correlation functions
    logical  :: CorrFunMode
    
    integer  :: NStepCorr
    real(RK) :: TimeStepCorr
    integer  :: NCorr,Mmess
    integer  :: NSpanCF,Nviewcf
    real(RK), pointer, contiguous :: cf_db(:), cf_soret(:)
    real(RK), pointer, contiguous :: average_cf_db(:), average_cf_soret(:)
    real(RK), pointer, contiguous :: cf_vs(:), cf_vb(:), cf_c(:), cf_ec(:)
    real(RK), pointer, contiguous :: average_cf_vs(:), average_cf_vb(:), average_cf_c(:), average_cf_ec(:)
    real(RK), pointer, contiguous :: lamda(:, :)
    real(RK), pointer, contiguous :: average_lamda(:, :)
    real(RK), pointer, contiguous :: sinte_i(:, :), sinte_lamda(:,:)
    real(RK), pointer, contiguous :: average_sinte_i(:, :), average_sinte_lamda(:,:)
    real(RK), pointer, contiguous :: sinte_db(:), sinte_soret(:)
    real(RK), pointer, contiguous :: average_sinte_soret(:), average_sinte_db(:)
    real(RK), pointer, contiguous :: sinte_vs(:), sinte_vb(:)
    real(RK), pointer, contiguous :: average_sinte_vs(:), average_sinte_vb(:)
    real(RK), pointer, contiguous :: sinte_c(:), sinte_ec(:)
    real(RK), pointer, contiguous :: average_sinte_c(:), average_sinte_ec(:)
    real(RK), pointer, contiguous :: a(:, :), A_SpanCF(:,:)
    real(RK), pointer, contiguous :: cf_d (:, :),  average_cf_d (:, :), vsk(:, :)
    real(RK),pointer, contiguous  :: vsp(:, :), vbk(:, :), vbp(:, :)
    real(RK), pointer, contiguous :: vckt(:, :), vckr(:, :), vcpt(:, :), vcpr(:, :), vcmt(:,:)
    real(RK)          :: sc(3),sp(3)

    real(RK),pointer, contiguous :: selfd_i(:)
    real(RK),pointer, contiguous :: Onsager(:,:)
    real(RK)         :: visco_s
    real(RK)         :: visco_b
    real(RK)         :: conduct
    real(RK)         :: soret
    real(RK)         :: econduct

    ! 4.) Transport properties

    type(TAccumulator),pointer, contiguous :: Sumself_i(:)
    type(TAccumulator),pointer, contiguous :: SumOnsager(:,:)
    type(TAccumulator)         :: SumVisco_s
    type(TAccumulator)         :: SumVisco_b
    type(TAccumulator)         :: SumSoret
    type(TAccumulator)         :: SumConduct
    type(TAccumulator)         :: SumEConduct
!TRANSPORT_END
#endif

#if CONSTR > 0
   integer         :: NCons
   integer,pointer, contiguous :: Cons1Comp(:)
   integer,pointer, contiguous :: Cons2Comp(:)
   integer,pointer, contiguous :: Cons1(:)
   integer,pointer, contiguous :: Cons2(:)
   real(RK),pointer, contiguous:: ConsR(:)
   real(RK),pointer, contiguous:: FCons(:)
   real(RK),pointer, contiguous:: UCons(:)
   logical         :: consup
#endif

#ifdef ABL
   real(RK),pointer, contiguous:: AblPS(:,:)
   real(RK),pointer, contiguous:: AblPE(:,:)
   real(RK),pointer, contiguous:: AblRhoS(:,:)
   real(RK),pointer, contiguous:: AblRhoE(:,:)
#endif

#if HBOND > 0
   integer          :: NHBondCrit
   integer,pointer, contiguous  :: AccComp(:), AccAccSite(:), AccDonSite(:)
   integer,pointer, contiguous  :: DonComp(:), DonAccSite(:), DonDonSite(:)
   real(RK),pointer, contiguous :: DistCrit1(:), DistCrit2(:), AngleCrit(:)
   integer,pointer, contiguous  :: NHBond0(:), NHBond1(:,:), NHBond2(:,:,:), NHBond3(:,:,:,:), NHBondN(:)
   type(TAccumulator),pointer, contiguous :: SumHBond0(:)
   type(TAccumulator),pointer, contiguous :: SumHBond1(:,:)
   type(TAccumulator),pointer, contiguous :: SumHBond2(:,:,:)
   type(TAccumulator),pointer, contiguous :: SumHBond3(:,:,:,:)
   type(TAccumulator),pointer, contiguous :: SumHBondN(:)
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

  interface LongRangeCheck
    module procedure TEnsemble_LongRangeCheck
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

  interface Unit2Atom
    module procedure TEnsemble_Unit2Atom
  end interface

  interface Atom2Unit
    module procedure TEnsemble_Atom2Unit
  end interface

  interface Mol2Unit
    module procedure TEnsemble_Mol2Unit
  end interface

  interface Unit2Mol
    module procedure TEnsemble_Unit2Mol
  end interface

  interface PredictVol
    module procedure TEnsemble_PredictVol
  end interface

  interface CorrectVol
    module procedure TEnsemble_CorrectVol
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

  interface QShake
    module procedure TEnsemble_QShake
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
    module procedure TEnsemble_UpdateEnergy1Mol
  end interface

  interface Energy
    module procedure TEnsemble_Energy
    module procedure TEnsemble_Energy1
    module procedure TEnsemble_Energy1Mol
    module procedure TEnsemble_Energy1_CF
    module procedure TEnsemble_EwaldEnergy1
  end interface

  interface GetEnergy
    module procedure TEnsemble_GetEnergy
    module procedure TEnsemble_GetEnergy1
    module procedure TEnsemble_GetEnergy1Mol
  end interface

  interface GetEnergyIntra
    module procedure TEnsemble_GetEnergyIntra
    module procedure TEnsemble_GetEnergyIntra1Mol
  end interface
  
  interface GetEnergyIntra_Bond
    module procedure TEnsemble_GetEnergyIntra_Bond
  end interface
  
  interface GetEnergyIntra_Angle
    module procedure TEnsemble_GetEnergyIntra_Angle
  end interface
  
  interface GetEnergyIntra_Dihedral
    module procedure TEnsemble_GetEnergyIntra_Dihedral
  end interface

  interface GetVirial
    module procedure TEnsemble_GetVirial
  end interface

  interface GetVirialIntra
    module procedure TEnsemble_GetVirialIntra
  end interface

  interface Getd2EpotdV2
    module procedure TEnsemble_Getd2EpotdV2
  end interface

  interface Move
    module procedure TEnsemble_Move
    module procedure TEnsemble_MoveMol
  end interface

  interface Rotate
    module procedure TEnsemble_Rotate
    module procedure TEnsemble_RotateMol
  end interface

  interface Move_NVE
    module procedure TEnsemble_Move_NVE
    module procedure TEnsemble_MoveMol_NVE
  end interface

  interface Rotate_NVE
    module procedure TEnsemble_Rotate_NVE
    module procedure TEnsemble_RotateMol_NVE
  end interface
  
  interface Move_NPH
    module procedure TEnsemble_Move_NPH
    module procedure TEnsemble_MoveMol_NPH
  end interface

  interface Rotate_NPH
    module procedure TEnsemble_Rotate_NPH
    module procedure TEnsemble_RotateMol_NPH
  end interface

  interface MoveBiased
    module procedure TEnsemble_MoveBiased
  end interface

  interface RotateBiased
    module procedure TEnsemble_RotateBiased
  end interface

  interface BiasedPartners
    module procedure TEnsemble_PartnersBiased
  end interface

  interface ChangeFluct
    module procedure TEnsemble_ChangeFluct
    module procedure TEnsemble_ChangeFluctTI
  end interface

  interface ScaleInteractionThermoInt
    module procedure TEnsemble_ScaleInteractionThermoInt
  end interface

  interface ChangeLambda
    module procedure TEnsemble_ChangeLambda
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

  interface ResizeMol
    module procedure TEnsemble_ResizeMol
  end interface

  interface Resize_Gibbs
    module procedure TEnsemble_Resize_liq
    module procedure TEnsemble_Resize_vap
  end interface

  interface Update_Gibbs
    module procedure TEnsemble_ResizeLiquid_Update
    module procedure TEnsemble_PartChangeUpdate
  end interface

  interface Insert_Gibbs
    module procedure TEnsemble_GibbsInsert
  end interface

  interface Remove_Gibbs
    module procedure TEnsemble_GibbsRemove
  end interface

  interface Gibbs_Delete
    module procedure TEnsemble_GibbsDelete
  end interface

  interface ZeroNAttempts
    module procedure TEnsemble_ZeroNAttempts
  end interface

  interface UpdateDisplacements
    module procedure TEnsemble_UpdateDisplacements
  end interface

  interface Residence
    module procedure TEnsemble_Residence
  end interface

  interface ResidencePartners
    module procedure TEnsemble_ResidencePartners
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

  interface RDFOpen
    module procedure TEnsemble_RDFOpen
  end interface

  interface RDFUpdate
    module procedure TEnsemble_RDFUpdate
  end interface

  interface RDFClose
    module procedure TEnsemble_RDFClose
  end interface

  interface ErrorsUpdate
    module procedure TEnsemble_ErrorsUpdate
  end interface

  interface ErrorsUpdateThermoInt
    module procedure TEnsemble_ErrorsUpdateThermoInt
  end interface

  interface SVCOutput
    module procedure TEnsemble_SVCOutput
  end interface

  interface VisualOpen
    module procedure TEnsemble_VisualOpen
  end interface

  interface VisualUpdate
    module procedure TEnsemble_VisualUpdate
#if HBOND > 0
    module procedure TEnsemble_VisualUpdateHB
#endif
  end interface

  interface VisualClose
    module procedure TEnsemble_VisualClose
  end interface

#if OSMOP > 0
  interface ProfileOpen
    module procedure TEnsemble_ProfileOpen
  end interface

  interface ProfileUpdate
    module procedure TEnsemble_ProfileUpdate
  end interface

  interface ProfileClose
    module procedure TEnsemble_ProfileClose
  end interface
#endif

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
    module procedure TEnsemble_EwaldFourierAddDel
  end interface

  interface EwaldSelfTerm_Energy
    module procedure TEnsemble_EwaldSelf_Energy
  end interface

#if SPME > 0
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
#endif

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

#if CONSTR > 0
  interface Constraints
    module procedure TEnsemble_Constraints
  end interface
#endif

#if HBOND > 0
  interface HBonding
    module procedure TEnsemble_HBonding
  end interface
#endif
  
contains

!==============================================================!
!  Subroutine TEnsemble_Construct                              !
!==============================================================!

  subroutine TEnsemble_Construct( this, ne )

    implicit none

    ! Include MPI header
#if HBOND > 0
#if MPI_VER > 0
    include 'mpif.h'
#endif
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: ne

    ! Declare local variables
    integer :: i, j
    integer :: stat
    character( IOBufferLength ) :: str
    
    !Declare variable for walltime solution in ms2_global
    integer :: time_limit

    ! Allocate simulation box length
    allocate( this%BoxLength, STAT = stat )
    call AllocationError( stat, 'simulation box length' )

    ! Allocate maximum number of particles
    allocate( this%NPartMax, STAT = stat )
    call AllocationError( stat, 'maximum number of particles' )
    allocate( this%NPartMaxFluct, STAT = stat )
    call AllocationError( stat, 'maximum number of particles' )

    ! Allocate maximum number of units
    allocate( this%NUnitMax, STAT = stat )
    call AllocationError( stat, 'maximum number of units' )

    ! Set number of ensemble
    this%EnsembleNumber = ne
    call LogWriteBlank
    write( IOBuffer, '(72(1H-))')
    call LogWrite
    write( IOBuffer, '(T14, "Reading parameters of ensemble", I3)' ) this%EnsembleNumber
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

    ! Read hamiltonian
    if( EnsembleType .eq. EnsembleTypeNVE ) then !.and. SimulationType .eq. MonteCarlo
      call FileReadParameter( this%RefHamiltonian, iounit_params , IdRefHamiltonian, .false. )
      if( .not. UseReducedUnits ) then
        this%RefHamiltonian = this%RefHamiltonian / UnitEnergy / NAvogadro
      end if
    end if

    ! Read Enthalpy (NPH)
    if( EnsembleType .eq. EnsembleTypeNPH ) then
      call FileReadParameter( this%RefEnthalpy, iounit_params , IdRefEnthalpy, .false. )
      if( .not. UseReducedUnits ) then
        this%RefEnthalpy = this%RefEnthalpy / UnitEnergy / NAvogadro 
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
    write( IOBuffer, '("Temperature: ",T26, F9.3, " K")' ) this%RefTemperature * UnitTemperature
    call LogWrite

    if( ConstantPressure ) then
      write( IOBuffer, '("Pressure: ",T26, F9.3, " MPa")' ) this%RefPressure * UnitPressure * 1E-6_RK
      call LogWrite
    end if

    if( EnsembleType .eq. EnsembleTypeGE ) then

      write( IOBuffer, '("Pressure0: ",T26, F12.6, " MPa")' ) this%RefPressure * UnitPressure * 1E-6_RK
      call LogWrite
      write( IOBuffer, '("Liquid density: ",T26, F12.6, " (", F13.6, ") mol/l")' ) & 
&       this%LiqDensity * UnitDensity, this%VarLiqDensity * UnitDensity
      call LogWrite
      write( IOBuffer, '("Liquid enthalpy: ",T22, F16.6, " (", F13.6, ") J/mol")' ) this%LiqEnthalpy * UnitEnergy * NAvogadro, &
&       this%VarLiqEnthalpy * UnitEnergy * NAvogadro
      call LogWrite
      write( IOBuffer, '("Liquid betaT: ",T26, F12.6, " (", F13.6, ") 1/MPa")' ) this%LiqBetaT / UnitPressure * 1E6_RK, &
&       this%VarLiqBetaT / UnitPressure * 1E6_RK
      call LogWrite
      write( IOBuffer, '("Liquid dHdP: ",T26, F12.6, " (", F13.6, ") l/mol")' ) &
&       this%LiqdHdP / UnitDensity, this%VarLiqdHdP / UnitDensity
      call LogWrite
    end if

    write( IOBuffer, '("Density: ", T26, F9.3, " mol/l")' ) this%RefDensity * UnitDensity
    call LogWrite
    if( EnsembleType .eq. EnsembleTypeNPH ) then
      write( IOBuffer, '("Enthalpy: ",T26, F9.3, " J/mol")' ) this%RefEnthalpy * UnitEnergy * NAvogadro
      call LogWrite
    end if
    write( IOBuffer, '("Reduced temperature: ", T26, F14.8)' ) this%RefTemperature
    call LogWrite

    if( ConstantPressure ) then
      write( IOBuffer, '("Reduced pressure: ",T26, F14.8)' ) this%RefPressure
      call LogWrite
    end if

    if( EnsembleType .eq. EnsembleTypeGE ) then
      write( IOBuffer, '("Reduced pressure0: ",T26, F14.8)' ) this%RefPressure
      call LogWrite
      write( IOBuffer, '("Reduced liquid density: ",T29, F11.8, " (", F14.8, ")")' ) this%LiqDensity, this%VarLiqDensity
      call LogWrite
      write( IOBuffer, '("Reduced liquid enthalpy: ", T28, F11.8, " (", F14.8, ")")' ) &
&       this%LiqEnthalpy, this%VarLiqEnthalpy
      call LogWrite
      write( IOBuffer, '("Reduced liquid betaT: ",T26, F14.8, " (", F14.8, ")")' ) this%LiqBetaT, this%VarLiqBetaT
      call LogWrite
      write( IOBuffer, '("Reduced liquid dHdP: ",T26, F14.8, " (", F14.8, ")")' ) this%LiqdHdP, this%VarLiqdHdP
      call LogWrite
    end if

    write( IOBuffer, '("Reduced density: ",T26, F14.8)' ) this%RefDensity
    call LogWrite
    if( EnsembleType .eq. EnsembleTypeNPH ) then
      write( IOBuffer, '("Reduced Enthalpy: ",T26, F12.6)' ) this%RefEnthalpy
      call LogWrite 
    end if

    ! Read mass of piston
    if( SimulationType .eq. MolecularDynamics .and. ConstantPressure ) then
      call FileReadParameter( this%PistonMass, iounit_params , IdPistonMass, .false. )
      if ( (.not. UseReducedUnits) .and. (parVersionNr .ge. 2.0_RK) ) then
        this%PistonMass = this%PistonMass / UnitMass * UnitLength**4
      end if
      write( IOBuffer, '("Reduced mass of piston: ",T26, F14.10)' ) this%PistonMass
      call LogWrite
    end if

    ! Read optional pressure calculation
    if( SimulationType .eq. MonteCarlo ) then
      if ( .not. ConstantPressure ) then
        this%OptPressure = .true.
      else
        call FileReadParameter( str, iounit_params , IdOptPressure, .false., "no" )
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
            call Error( 'Select yes/no for calculation of pressure '// ProgramFileName//ConfigFileExtension )
        end select
      end if
    else
      this%OptPressure = .true.
    end if

    ! Read calculation of residence time
    this%ResidenceTime = .false.
    if( SimulationType .eq. MolecularDynamics ) then
      call FileReadParameter( str, iounit_params , IdResidTime, .false., "no" )
      select case( str )
        case( 'YES', 'Yes', 'yes' )
          this%ResidenceTime = .true.
          write( IOBuffer, '("Calculation of residence time: ",A)' ) trim( str )
          call LogWrite
          write( IOBuffer, '("Maximum molecules in hydration shell: 10 (default in code)")' ) 
          call LogWrite
          call FileReadParameter( this%ResidComp1, iounit_params , IdResidComp1, .false. )
          call FileReadParameter( this%ResidSite1, iounit_params , IdResidSite1, .false. )
          call FileReadParameter( this%ResidComp2, iounit_params , IdResidComp2, .false. )
          call FileReadParameter( this%ResidSite2, iounit_params , IdResidSite2, .false. )
          write( IOBuffer, '("For component",I2,", site:",I2)' ) this%ResidComp1, this%ResidSite1
          call LogWrite
          write( IOBuffer, '("and component",I2,", site:",I2)' ) this%ResidComp2, this%ResidSite2
          call LogWrite
          call FileReadParameter( this%ResidPeriod, iounit_params , IdResidPeriod, .false. )
          write( IOBuffer, '("Update period / steps: ",T25, I6)' ) this%ResidPeriod
          call LogWrite
          call FileReadParameter( this%ResidLength, iounit_params , IdResidLength, .false. )
          write( IOBuffer, '("Pairing at distances lower: ",T28, F9.5)' ) this%ResidLength
          call LogWrite
          this%ResidLength = this%ResidLength / UnitLength * Angstroem
          call FileReadParameter( this%ResidBreak, iounit_params , IdResidBreak, .false., 0 )
          write( IOBuffer, '("Pairing can be invalid for steps: ",T28, I7)' ) this%ResidBreak
          call LogWrite

        case( 'NO', 'No', 'no')
          this%ResidenceTime = .false.
          write( IOBuffer, '("Calculation of residence time: ",A)' ) trim( str )
          call LogWrite
        case default
          call Error( 'Select yes/no for calculation of residence time'// ProgramFileName//ConfigFileExtension )
      end select
    end if

    ! Read whether to perform the MC equilibration in parallel
    call FileReadParameter( str, iounit_params , IdCommonEqui, .false., "no" )
    select case( str )
      case( 'YES', 'Yes', 'yes' )
        CommonEqui = .true.
        if( SimulationType .eq. MonteCarlo ) then
          write( IOBuffer, '("Common equilibration: ",T30, A)' ) trim( str )
          call LogWrite
        endif
      case( 'NO', 'No', 'no')
        if( SimulationType .eq. MolecularDynamics ) then
          CommonEqui = .true.
          write(IOBuffer, '("MD simulation: Logical CommonEqui no is invalid, set to yes")' )
          call LogWrite
        else
          CommonEqui = .false.
          write( IOBuffer, '("Common equilibration: ",T30, A)' ) trim( str )
          call LogWrite
        endif
      case default
        call Error( 'Select yes/no for common equilibration '// ProgramFileName//ConfigFileExtension )
    end select
  
    ! Read initial number of particles in ensemble
    call FileReadParameter( this%NPart, iounit_params , IdNPart, .false. )
    if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then

      this%NPartInitial = this%NPart
      this%NPartLBound = int( real( this%NPart, RK ) / 1.2_RK )
      this%NPartUBound = int( real( this%NPart, RK ) * 1.2_RK )

    end if

    ! Read number of components in ensemble
    call FileReadParameter( this%NComponents, iounit_params , IdNComponents, .false. )
    write( IOBuffer, '("Number of components:",T28, I3)' ) this%NComponents
    call LogWrite
    if( this%NComponents <= 0 ) then
      write( ErrorBuffer, '("There must be at least 1 component in ensemble", I2)' ) this%EnsembleNumber
      call Error
    end if

    if( this%NComponents > 999 ) call Error( 'Cannot work with more than 999 components on '//Hardware )

#if  TRANS == 1
!TRANSPORT_start
    call LogWriteBlank
    if ( parVersionNr .ge. 2.0_RK ) then
      call FileReadParameter( str , iounit_params , IdCorrFun, .false. , 'no' )
      select case( str )

      case( 'yes' , 'ok', 'ja' )
        this%CorrFunMode = .true.
        str = 'Include transport properties'

      case( 'no', 'nein' )
        this%CorrFunMode = .false.
        str = 'No transport properties'

      case default
        call Error( 'Unknown transport properties ('//trim(IdCorrFun)//'='//trim(str)//')' )
      end select
      write( IOBuffer, '("Transport properties:",T26, A)' ) trim(str)
      call LogWrite
    endif

    ! Read correlation function
    if( this%CorrfunMode ) then

      ! Calculate correlation function every n-th time step
      call FileReadParameter( this%NStepCorr , iounit_params , IdNStepcf )

      ! Read legth of the correlation function
      call FileReadParameter( this%NCorr , iounit_params , IdCorrlength )

      ! Read time span between correlations
      call FileReadParameter( this%NSpanCF , iounit_params , IdSpanCF )

      ! Calculation of the correlation function every n-th time step

      if(mod(this%NSpanCF, this%NStepCorr) .eq. 0) then
        this%NSpanCF = this%NSpanCF/this%NStepCorr
        this%NCorr = this%NCorr/this%NStepCorr
        write( IOBuffer, '("Correlation Function (CF) is calculated every",I4,"-th time step")') this%NStepCorr
        call LogWrite
      else
        this%NStepCorr = 1
        write( IOBuffer, '("Correlation Function (CF) is calculated every time step")')
        call LogWrite
        write( IOBuffer, '("StepsCorrfun is set to 1. SpanCorrfun is not divisible by StepsCorrfun")') 
        call LogWrite
      endif

      this%TimeStepCorr = TimeStep * this%NStepCorr

      if(mod(this%NCorr, this%NSpanCF) .eq. 0) then
        write( IOBuffer, '("Length of CF:",T26, I5)' ) this%NCorr*this%NStepCorr
        call LogWrite
      else
        this%NCorr = (AINT(real( this%NCorr, RK )/real( this%NSpanCF, RK ))+1)*this%NSpanCF
        write( IOBuffer, '("Length of CF is extended to:",T40, I7)') this%NCorr*this%NStepCorr
        call LogWrite
      endif
      
      ! Correlation length output
      write( IOBuffer, '("Time Span between CF:",T26, I5)' ) this%NSpanCF*this%NStepCorr
      call LogWrite

      call FileReadParameter( this%Nviewcf , iounit_params , IdNviewcf )
      write( IOBuffer, '("Print CF each:",T26, I5)' ) this%Nviewcf
      call LogWrite

      if ( ((this%Nviewcf*this%NSpanCF*this%NStepCorr+this%NCorr*this%NStepCorr) > NSteps) .or. (this%Nviewcf .eq. 0) ) then
        write(IOBuffer, '("Warning: Updates of CF not sufficient - Output once at the end of simulation")')
        call LogWrite
        this%Nviewcf = int((NSteps-this%NCorr*this%NStepCorr)/(this%NSpanCF*this%NStepCorr))
        write( IOBuffer, '("Print after", I6," CF")' ) this%Nviewcf
        call LogWrite
      end if

      ! Read frequency of updating result file CF
      call FileReadParameter( BlockSizeCF , iounit_params , IdBlockSizeCF, .false., 1 )
      if( BlockSizeCF > 0 ) then
        write( IOBuffer, '("Result files will be updated each", I3, " CF")' ) BlockSizeCF
      else
        write( IOBuffer, '("Result files will not be created")' )
      end if
      call LogWrite

      ! Calculate max of Mmess for BlockSize determination
      this%Mmess = int((((NSteps+this%NStepCorr-1)/this%NStepCorr)-this%NCorr)/this%NSpanCF)

      if( BlockSizeCF > 0 ) then
        NBlocksMaxCF = int(this%Mmess / BlockSizeCF)
        NBlockSizesMaxCF = int( sqrt( real(NBlocksMaxCF,RK) ) )

      else
        NBlocksMaxCF = 0
        NBlockSizesMaxCF = 0
      end if

      ! Initialization of Mmess
      this%Mmess = 0
    end if
!TRANSPORT_END
#endif

#if MPI_VER > 0
    if ((LongRange .eq. Ewald) .or. (LongRange .eq. PME))then
      allocate(this%NBox0,STAT=stat)
      call AllocationError( stat, 'NProcs' )
      allocate(this%NBox1,STAT=stat)
      call AllocationError( stat, 'NProcs' )
      allocate(this%NBox2,STAT=stat)
      call AllocationError( stat, 'NProcs' )
    end if
#endif
    ! Create components
    call CreateComponents( this )

    ! Calculate initial number of particles of each component
    call CalculateNPart( this )

    ! LongRange-Corrections need electro-neutral system
    call LongRangeCheck ( this )

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
    ! Setting scale coefficients for ThermoInt-Components
    j = this%NRealComponents+1
    do i = 1, this%NRealComponents
      if (this%Component(i)%ChemPotMethod == ChemPotMethodThermoInt ) then
        this%ScaleSigma(j, :) = this%ScaleSigma(i, :)
        this%ScaleSigma(:, j) = this%ScaleSigma(:, i)
        this%ScaleEpsilon(j, :) = this%ScaleEpsilon(i, :)
        this%ScaleEpsilon(:, j) = this%ScaleEpsilon(:, i)
        j = j+1
      end if
    end do

#if HBOND > 0
    call FileReadParameter( this%NHBondCrit, iounit_params , IdNHBonds, .false. )
    call LogWriteBlank
    write( IOBuffer, '("Reading ", I2, " H-Bonding criteria:")' ) this%NHBondCrit
    call LogWrite
    write( IOBuffer, '("CritNo.  AccComp  AccAccSite AccDonSite  DonComp  DonAccSite  DonDonSite  DistCrit1  DistCrit2  AngleCrit")' )
    call LogWrite

    !Allocate H-Bond Arrays
    allocate( this%AccComp( this%NHBondCrit ), STAT = stat )
    call AllocationError( stat, 'components', this%NHBondCrit )
    allocate( this%AccAccSite( this%NHBondCrit ), STAT = stat )
    call AllocationError( stat, 'components', this%NHBondCrit )
    allocate( this%AccDonSite( this%NHBondCrit ), STAT = stat )
    call AllocationError( stat, 'components', this%NHBondCrit )
    allocate( this%DonComp( this%NHBondCrit ), STAT = stat )
    call AllocationError( stat, 'components', this%NHBondCrit )
    allocate( this%DonAccSite( this%NHBondCrit ), STAT = stat )
    call AllocationError( stat, 'components', this%NHBondCrit )
    allocate( this%DonDonSite( this%NHBondCrit ), STAT = stat )
    call AllocationError( stat, 'components', this%NHBondCrit )
    allocate( this%DistCrit1( this%NHBondCrit ), STAT = stat )
    call AllocationError( stat, 'components', this%NHBondCrit )
    allocate( this%DistCrit2( this%NHBondCrit ), STAT = stat )
    call AllocationError( stat, 'components', this%NHBondCrit )
    allocate( this%AngleCrit( this%NHBondCrit ), STAT = stat )
    call AllocationError( stat, 'components', this%NHBondCrit )

    if (RootProc) then
      do i = 1, this%NHBondCrit
        read( iounit_params, * ) this%AccComp(i), this%AccAccSite(i), this%AccDonSite(i), this%DonComp(i), &
  &           this%DonAccSite(i), this%DonDonSite(i), this%DistCrit1(i), this%DistCrit2(i), this%AngleCrit(i)
        if ( (this%AngleCrit(i) .le. 0._RK) .or. (this%AngleCrit(i) .gt. 180._RK) ) then
          call Error('Angle of the H-Bonding criteria(s) should be between 0 and 180.')
        else
          this%AngleCrit(i) = cos(this%AngleCrit(i)*PI/180._RK)
        end if
        write( IOBuffer, '("  ", I3, "      ", I2, "         ", I2, "         ", I2, "       ", I2, "         ", I2, "         ", I2, "     ", F9.4, "   ", F9.4, "   ", F9.4)' ) i, &
  &              this%AccComp(i), this%AccAccSite(i), this%AccDonSite(i), this%DonComp(i), this%DonAccSite(i), this%DonDonSite(i), &
  &              this%DistCrit1(i), this%DistCrit2(i), this%AngleCrit(i)
        call LogWrite
        if (.not. UseReducedUnits ) then
          this%DistCrit1(i) = this%DistCrit1(i)*Angstroem/UnitLength
          this%DistCrit2(i) = this%DistCrit2(i)*Angstroem/UnitLength
        end if
      end do
    end if
#if MPI_VER > 0
    call MPI_Bcast( this%AccComp,    this%NHBondCrit, MPI_INTEGER, NRootProc, Communicator, ierror )
    call MPI_Bcast( this%AccAccSite, this%NHBondCrit, MPI_INTEGER, NRootProc, Communicator, ierror )
    call MPI_Bcast( this%AccDonSite, this%NHBondCrit, MPI_INTEGER, NRootProc, Communicator, ierror )
    call MPI_Bcast( this%DonComp,    this%NHBondCrit, MPI_INTEGER, NRootProc, Communicator, ierror )
    call MPI_Bcast( this%DonAccSite, this%NHBondCrit, MPI_INTEGER, NRootProc, Communicator, ierror )
    call MPI_Bcast( this%DonDonSite, this%NHBondCrit, MPI_INTEGER, NRootProc, Communicator, ierror )
    call MPI_Bcast( this%DistCrit1,  this%NHBondCrit, MPI_RK, NRootProc, Communicator, ierror )
    call MPI_Bcast( this%DistCrit2,  this%NHBondCrit, MPI_RK, NRootProc, Communicator, ierror )
    call MPI_Bcast( this%AngleCrit,  this%NHBondCrit, MPI_RK, NRootProc, Communicator, ierror )
#endif
#endif

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
      write( IOBuffer, '("Reduced center of mass cutoff radius: ",T45, F6.3)' ) this%RCutoffLJ126LJ126
      call LogWrite
      this%RCutoffDipoleDipole = this%RCutoffLJ126LJ126
      this%RCutoffDipoleQuadrupole = this%RCutoffLJ126LJ126
      this%RCutoffQuadrupoleQuadrupole = this%RCutoffLJ126LJ126

    else

      if( this%NLJ126Max > 0 ) then
        call FileReadParameter( this%RCutoffLJ126LJ126, iounit_params , IdRCutoffLJ126LJ126, .false. )
        write( IOBuffer, '("Lennard-Jones cutoff radius: ",T45, F6.3, " sigma")' ) this%RCutoffLJ126LJ126
        call LogWrite
      end if

      if( this%NDipoleMax > 0 ) then
        call FileReadParameter( this%RCutoffDipoleDipole, iounit_params , IdRCutoffDipoleDipole, .false. )
        write( IOBuffer, '("Reduced dipole-dipole cutoff radius: ",T42, F8.3)' ) this%RCutoffDipoleDipole
        call LogWrite

        if( this%NQuadrupoleMax > 0 ) then
          call FileReadParameter( this%RCutoffDipoleQuadrupole, iounit_params , IdRCutoffDipoleQuadrupole, .false. )
          write( IOBuffer, '("Reduced dipole-quadrupole cutoff radius: ",T42, F8.3)' ) this%RCutoffDipoleQuadrupole
          call LogWrite
        end if

      end if

      if( this%NQuadrupoleMax > 0 ) then
        call FileReadParameter( this%RCutoffQuadrupoleQuadrupole, iounit_params , IdRCutoffQuadrupoleQuadrupole, .false. )
        write( IOBuffer, '("Reduced quadrupole-quadrupole cutoff radius: ",T42, F8.3)' ) this%RCutoffQuadrupoleQuadrupole
        call LogWrite
      end if
    end if

    ! Read characteristic dielectric constant
    this%RFEpsilon = 0._RK

    if(( this%NDipoleMax > 0 ) .or. ( this%NChargeMax > 0 )) then
      call LogWriteBlank
      call FileReadParameter( this%RFEpsilon, iounit_params , IdRFEpsilon, .false. )
      write( IOBuffer, '("Characteristic dielectric constant: ",T41, E16.5)' ) this%RFEpsilon
      call LogWrite
    end if

! Initialization of the transport property "Conductivity" and "EConductivity"
#if TRANS == 1
      
     this%Bulkviscosity = .false.
     this%MolarEnthConduct = .true.
     this%Conductivity = .false.
     this%EConductivity = .false.

     if (LongRange .eq. Rfield) this%Conductivity = .true.
     if (EnsembleType .eq. EnsembleTypeNVE) this%Bulkviscosity = .true.
 
      if ( this%NComponents > 1 .and. LongRange .eq. RField) then
       do i = 1, this%NComponents
            if (this%Component(i)%PartialMolarEnthalpy .eq. 0._RK) then
               this%Conductivity = .false.
               this%MolarEnthConduct = .false.
            end if
       end do
      end if

      
      if (LongRange .eq. Ewald) then
	do i = 1, this%NComponents
	  if ( abs(this%Component(i)%Molecule%Charge) .gt. 1e-7) then
             this%EConductivity = .true.
	  end if
	end do
      end if


#endif

#if CONSTR > 0
    write( IOBuffer, '("CONSTRAINED DYNAMICS")' )
    call LogWrite

    call FileReadParameter( iounit_params , IdNCons )
    read( IOBuffer, * ) this%NCons
    write( IOBuffer, '("Number of Constrained Molecules:", I3)' ) this%NCons
    call LogWrite

    allocate( this%Cons1Comp( this%NCons), STAT = stat )
       if(stat >0) write(*,*) 'Allocation Error Cons1Comp'
    allocate( this%Cons2Comp( this%NCons), STAT = stat )
       if(stat >0) write(*,*) 'Allocation Error Cons2Comp'
    allocate( this%Cons1( this%NCons), STAT = stat )
           if(stat >0) write(*,*) 'Allocation Error Cons1'
    allocate( this%Cons2( this%NCons), STAT = stat )
         if(stat >0) write(*,*) 'Allocation Error Cons2'
    allocate( this%ConsR( this%NCons), STAT = stat )
         if(stat >0) write(*,*) 'Allocation Error ConsR'
    allocate( this%FCons( this%NCons), STAT = stat )
         if(stat >0) write(*,*) 'Allocation Error FCons'
    allocate( this%UCons( this%NCons), STAT = stat )
         if(stat >0) write(*,*) 'Allocation Error UCons'

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
    call LogWriteBlank
    write( IOBuffer, '("Cutoff correction to")' )
    call LogWrite

    if ( SimulationType .eq. MonteCarlo .and. (.not.  CommonEqui))  then
      write( IOBuffer, '("- potential energy from LJ",T44, F12.8)' ) this%EPotCorrLJ  / this%NPart

    else
      write( IOBuffer, '("- potential energy from LJ",T44, F12.8)' ) this%EPotCorrLJ * Nprocs/ this%NPart
    endif

    call LogWrite

    if ( SimulationType .eq. MonteCarlo .and. (.not. CommonEqui))  then  
      write( IOBuffer, '("- pressure from LJ ",T44, F12.8)' ) this%VirialCorrLJ  / this%NPart

    else
      write( IOBuffer, '("- pressure from LJ ",T44, F12.8)' ) this%VirialCorrLJ * NProcs / this%NPart
    endif

    call LogWrite

    do i = 1, this%NRealComponents
      write( IOBuffer, '("- chem. pot. of ", A, " from LJ",T44, F12.8)' ) trim( this%Component(i)%PotModFileName ), &
&        this%Component(i)%EPotTestCorrLJ
      call LogWrite
    end do

    if ( SimulationType .eq. MonteCarlo .and. (.not. CommonEqui))  then 
      write( IOBuffer, '("- potential energy from reaction field (RF)",T44, F12.8)' ) &
&       this%EPotCorrRF  / this%NPart

    else
      write( IOBuffer, '("- potential energy from reaction field (RF)",T44, F12.8)' ) &
&       this%EPotCorrRF * NProcs / this%NPart

    endif

    call LogWrite

    write( IOBuffer, '("- pressure from reaction field:",T44, F12.8)' ) this%VirialCorrRF / this%NPart
    call LogWrite

    do i = 1, this%NRealComponents
      write( IOBuffer, '("- chem. pot. of ", A, " from RF",T44, F12.8)' ) trim( this%Component(i)%PotModFileName ), &
&       this%Component(i)%EPotTestCorrRF
      call LogWrite
    end do
    
    ! Calculate maximum cutoff radius
    if( this%NDipoleMax > 0 ) then
      this%RCutoffMax2 = max(this%RCutoffMax2, 2._RK * this%RCutoffDipoleDipole )

      if( this%NQuadrupoleMax > 0 ) then
        this%RCutoffMax2 = max(this%RCutoffMax2, 2._RK * this%RCutoffDipoleQuadrupole )
      end if

    end if

    if( this%NQuadrupoleMax > 0 ) then
      this%RCutoffMax2 = max(this%RCutoffMax2, 2._RK * this%RCutoffQuadrupoleQuadrupole )
    end if

    if( .not. Restart ) then

      ! Update all BoxLength-dependent constants
      call UpdateBoxLength( this )

      ! Abort, if maximum cutoff larger than boxlength 
      if (this%RCutoffMax2 > this%BoxLength) call Error('Cutoff is larger than the boxsize')

      ! Set initial positions of particles in simulation box
      call InitPositions( this )

      ! Set initial orientations of particles in simulation box
      call InitOrientations( this )

      ! Calculate initial energies for the Ewald Summation
      if (LongRange .eq. Ewald) then

         if (this%KappaL .eq. 0.) then
            this%Kappa = sqrt(PI) * (4.0_RK*this%NPart / this%Volume0**2)**(1._RK/6._RK)
            this%KappaL = this%Kappa*this%BoxLength
         else
            this%Kappa = this%KappaL/this%BoxLength   !Boxlength bereits normiert
         end if

         do i=1,this%NComponents
           do j=1,this%NComponents
             this%Interaction(i,j)%Kappa = this%Kappa
           end do
         end do

         call EwaldSelfTerm( this )

        ! Memory Allocation for Ewald Summation
         nullify ( this%U_fourierLocal )
         nullify ( this%SSin )
         nullify ( this%SCos )
         nullify ( this%rold )
         nullify ( this%Vec2 )
         nullify ( this%VirIntra )

         allocate(this%U_fourierLocal(this%BoxenAnzahlMax),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error U_fourier'
         allocate(this%SSin(this%BoxenAnzahlMax),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error SSin'
         allocate(this%SCos(this%BoxenAnzahlMax),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error SCos'
         allocate(this%rold(5,3),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error rold'
         allocate(this%Vec2(this%BoxenAnzahlMax),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error Vec2'
        allocate(this%VirIntra(this%NPartMax),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error VirIntra'

#if SPME > 0
      else if (LongRange .eq. PME) then
         if (this%KappaL .eq. 0.) then
            this%Kappa = sqrt(PI) * (4.0_RK*this%NPart / this%Volume0**2)**(1._RK/6._RK)
            this%KappaL = this%Kappa*this%BoxLength

         else
            this%Kappa = this%KappaL/this%BoxLength   !Boxlength bereits normiert
         end if

         do i=1,this%NComponents
           do j=1,this%NComponents
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
#endif
      else if (LongRange .eq. ExtRField) then
         do i=1,this%NComponents
           do j=1,this%NComponents
             this%Interaction(i,j)%DebyeLen = this%DebyeLen
           end do
         end do

      else if (LongRange .eq. Rodgers) then
         do i=1,this%NComponents
           do j=1,this%NComponents
             this%Interaction(i,j)%Kappa = this%Kappa
           end do
         end do
      end if
     
      if( SimulationType .eq. MolecularDynamics .and. .not. MCOverlapReduction ) then

        ! Calculate positions of units
        call Mol2Unit( this)   ! Calculate initial orientations and positions of units

        ! Initialize molecular dynamics simulation
        call InitMolecularDynamics( this, .false. )

      else

        ! Set temperature
        this%Temperature = this%RefTemperature

        ! Convert unit coordinates to atom positions
        call Unit2Atom( this )

        ! Set all potential energy matrices
        call Energy( this, this%EPot )
        call UpdateEnergy( this )

        ! Set initial values of maximum allowed MC displacements
        this%DispVol = DispVolStart
        do i = 1, this%NRealComponents
          if (this%Component(i)%Molecule%NUnit .eq. 1) then
            this%Component(i)%DispTran = DispMolTranStart
            this%Component(i)%DispRot = DispMolRotStart
          else
            this%Component(i)%DispTran = DispTranStart
            this%Component(i)%DispRot = DispRotStart
            this%Component(i)%DispMolTran = DispMolTranStart
            this%Component(i)%DispMolRot = DispMolRotStart
          end if
        end do

      end if
    else        ! Restart
      ! Calculate initial energies for the Ewald Summation
    end if         ! Restart

    ! Set I/O unit numbers
    i = FilesPerEnsemble * this%EnsembleNumber
    this%iounit_result    = iounit_result    + i
    this%iounit_runave    = iounit_runave    + i
    this%iounit_errors    = iounit_errors    + i
    this%iounit_visual    = iounit_visual    + i
    this%iounit_rdf       = iounit_rdf       + i
    this%iounit_thermoint = iounit_thermoint + i
    this%iounit_rescf     = iounit_rescf     + i
    this%iounit_visualHB  = iounit_visualHB  + i
    this%iounit_dcp       = iounit_dcp       + i

    ! Calculate RDF VSchale 
    this%RDFdr = this%RCutoffLJ126LJ126 / RDFNumberShells
    do i = 1, RDFNumberShells
      this%RDFVSchale(i) = 4./3.*pi* this%RDFdr**3 *(i**3 - (i-1)**3)
    end do

    write( IOBuffer, '(T15, "Reading ensemble ", I3, " successful")') this%EnsembleNumber
    call LogWrite
    write( IOBuffer, '(72(1H-))')
    call LogWrite

#if MPI_VER > 0
! Abortion of simulation run due to wall-time Constraints
    call time_left(time_limit)
#endif



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
    ! Allocate maximum number of units
    allocate( this%NUnitMax, STAT = stat )
    call AllocationError( stat, 'maximum number of units' )

    ! Set number of ensemble
    this%EnsembleNumber = ne
    call LogWriteBlank
    write( IOBuffer, '(72(1H-))')
    call LogWrite
    write( IOBuffer, '(T14, "Reading parameters of ensemble", I3)' ) this%EnsembleNumber
    call LogWrite

    ! Read temperature
    call FileReadParameter( this%Temperature, iounit_params , IdRefTemperature, .false. )

    if( .not. UseReducedUnits ) then
      this%Temperature = this%Temperature / UnitTemperature
    end if

    ! Update log file
    write( IOBuffer, '("Temperature: ",T26, F9.3, " K")' ) this%RefTemperature * UnitTemperature
    call LogWrite
    write( IOBuffer, '("Reduced temperature: ",T26, F12.6)' ) this%RefTemperature
    call LogWrite

    ! Read number of components in ensemble
    call FileReadParameter( this%NComponents, iounit_params , IdNComponents, .false. )
    write( IOBuffer, '("Number of components:",T28, I3)' ) this%NComponents
    call LogWrite

    if( this%NComponents <= 0 ) then
      write( ErrorBuffer, '("There must be at least 1 component in ensemble", I2)' ) &
&       this%EnsembleNumber
      call Error
    end if

    if( this%NComponents > 999 ) call Error( 'Cannot work with more than 999 components on '//Hardware )

    ! Create components
    this%NComponents = 2 * this%NComponents
    allocate( this%Component(this%NComponents), STAT = stat )
    call AllocationError( stat, 'components', this%NComponents )

    do i = 1, this%NComponents, 2

      ! Read file name for potential model
      call FileReadParameter( PotModFileName, iounit_params , IdPotModFileName, .false. )
      call Construct( this%Component(i), PotModFileName )
      call Construct( this%Component(i+1), PotModFileName )

      if( (this%Component(i)%ChemPotMethod .eq. ChemPotMethodWidom) .and. (LongRange .ne. RField) ) &
&     call Error( 'Widom cannot be used with Ewald Summation or its Deriavtes ' )

    end do

    ! Calculate number of particles in each process
    do i = 1, this%NComponents
      pc => this%Component(i)
      pc%NPart = NOrient
      pc%NPart1 = ProcRange( pc%NPart, pc%NPart0, pc%NPart2 )
    end do

    this%NPartMax = NOrient
    this%NTestMax = 0
    ! Calculate Max Number of Units in Component
    this%NUnitMax = 0
    do i = 1, this%NComponents
      pc => this%Component(i)
      pc%NUnitMax => pc%Molecule%NUnit
      if (pc%Molecule%NUnit > this%NUnitMax) this%NUnitMax = pc%Molecule%NUnit
    end do

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
        write( IOBuffer, '(A, "-", A, " Lennard-Jones interaction:  eta =", F6.3, ", xi =", F6.3)' ) &
&         trim( this%Component(i)%PotModFileName ), trim( this%Component(j)%PotModFileName ), &
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

        this%Interaction(i, j)%EPotCorrLJ = sum( this%Interaction(i, j)%PotLJ126LJ126(:, :)%EPotCorr )
        write( IOBuffer, '("Cutoff correction to SVC of ", A, "-", A, " from LJ:", F12.8)' ) &
&         trim( this%Component(i)%Molecule%PotModFileName ), trim( this%Component(j)%Molecule%PotModFileName ), &
&         .5_RK * this%Interaction(i, j)%EPotCorrLJ / this%Temperature
        call LogWrite

      end do
    end do

    ! Update all BoxLength-dependent constants
    call UpdateBoxLength( this )

    ! Set initial positions of particles
    do i = 1, this%NComponents
      this%Component(i)%Pm0 = 0._RK
    end do

    ! Set initial orientations of particles
    call InitOrientations( this )

    ! Convert molecular coordinates to atom positions
    call Unit2Atom( this )

    ! Set I/O unit numbers
    i = FilesPerEnsemble * this%EnsembleNumber
    this%iounit_result = iounit_result + i
    this%iounit_runave = iounit_runave + i
    this%iounit_errors = iounit_errors + i
    this%iounit_visual = iounit_visual + i
    this%iounit_visualHB = iounit_visualHB + i

  end subroutine TEnsemble_ConstructSVC


!==============================================================!
!  Subroutine TEnsemble_Destruct                               !
!==============================================================!

  subroutine TEnsemble_Destruct( this )

    implicit none
#if SPME > 0
    include 'fftw3.f'
#endif

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

    if( associated( this%NPartMaxFluct ) ) then
      deallocate( this%NPartMaxFluct )
    end if

    ! Deallocate maximum number of units
    if( associated( this%NUnitMax ) ) deallocate( this%NUnitMax )

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

      if ( associated ( this%U_fourierLocal ) ) then 
         deallocate( this%U_fourierLocal )
      end if

      if ( associated ( this%SSin ) ) then 
         deallocate( this%SSin )
      end if

      if ( associated ( this%SCos ) ) then 
         deallocate( this%SCos )
      end if

      if ( associated ( this%rold ) ) then 
         deallocate( this%rold )
      end if

      if ( associated ( this%Vec2 ) ) then 
         deallocate( this%Vec2 )
      end if

      if ( associated ( this%VirIntra ) ) then
         deallocate( this%VirIntra )
      end if

    end if

#if SPME > 0
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
#endif

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
    real(RK):: q
    type(TComponent), pointer, contiguous :: reallocate(:)

    ! Create components
    ncomp = this%NComponents
    allocate( this%Component(ncomp), STAT = stat )
    call AllocationError( stat, 'components', ncomp )
    do i = 1, ncomp
      call Construct( this%Component(i), i )
    end do

    ! Create components for fluctuating particle states
    this%NGradInsComp = 0
    this%NRealComponents = ncomp
    do i = 1, this%NRealComponents
      q = 0._RK

      if( this%Component(i)%ChemPotMethod .eq. ChemPotMethodGradIns ) then
        ! LongRange Check
        do j=1,this%Component(i)%Molecule%NCharge
          q = q + this%Component(i)%Molecule%SiteCharge(j)%e
        end do

        if (abs(q) .gt. 1e-7) call Error ('Gradual Insertion not possible for charged molecule! No Electroneutrality')

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
          call Construct( this%Component(ncomp - nfluct + j), this%Component(i), j )
          this%Component(i)%NFluctComp(j) = ncomp - nfluct + j
        end do
        this%Component(i)%NFluctComp(0) = i
      end if

      if( this%Component(i)%ChemPotMethod .eq. ChemPotMethodThermoInt ) then

        ! LongRange Check
        do j=1,this%Component(i)%Molecule%NCharge
          q = q + this%Component(i)%Molecule%SiteCharge(j)%e
        end do

        if (abs(q) .gt. 1e-7) call Error ('Thermodynamic Integration not possible for charged molecule! No Electroneutrality')

        ncomp = ncomp + 1

        ! Reallocate component array
        allocate( reallocate(ncomp), STAT = stat )
        ! reallocate will be stored in this%Component, which thus has to be deallocated
        call AllocationError( stat, 'components', ncomp )
        reallocate( 1:size(this%Component) ) = this%Component(:)
        deallocate( this%Component )
        this%Component => reallocate
        call Construct( this%Component(ncomp), this%Component(i))
      end if

    end do

#if OSMOP > 0
    if (SimulationType .eq. MolecularDynamics ) then
      j = 0._RK
      do i = 1, this%NRealComponents
        if (this%Component(i)%permeable) then
          j = 1._RK
        end if
      end do
      if (j == 0._RK ) call Error ('At least one component has to be permeable.')
    end if
#endif
#if OSMOP == 2
    do i = 1, this%NRealComponents
      if ( this%Component(i)%ChemPotMethod .ne. ChemPotMethodWidom ) then
        call LogWriteBlank
        write( IOBuffer, '("Chem. Pot.-Profile for", A, " cannot be calculated if Widom is not used!")' ) trim(this%Component(i)%PotModFileName)
        call LogWrite
        call LogWriteBlank
        this%Component(i)%ChemPotProfile(:) = 0._RK
      end if
    end do
#endif

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
    !if( associated( this%Component ) ) then  !Michael Sch.: fix me
    !  deallocate( this%Component )
    !end if

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
#if OSMOP == 2
    integer :: j1, j2
#endif

    ! Create interactions
    allocate( this%Interaction(this%NComponents, this%NComponents ), STAT = stat )
    call AllocationError( stat, 'components', this%NComponents )

    if( SimulationType .eq. SecondVirialCoeff ) then
      do i = 1, this%NComponents, 2
        do j = i + 1, this%NComponents, 2
          this%Interaction(i,j)%OptPressure = this%OptPressure
          call Construct(this%Interaction(i, j), i, j, &
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

          if (LongRange .eq. ExtRField) this%Interaction(i,j)%DebyeLen = this%DebyeLen
          this%Interaction(i,j)%OptPressure = this%OptPressure
          call Construct(this%Interaction(i, j), i, j, &
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

#if OSMOP == 2
    ! set pointers for virial profile
    do i = 1, this%NComponents
      do j = 1, this%NComponents
      
        if( this%Interaction(i,j)%N1LJ126 > 0 .and. this%Interaction(i,j)%N2LJ126 > 0 ) then
          do j1 = 1, this%Interaction(i,j)%N1LJ126
            do j2 = 1, this%Interaction(i,j)%N2LJ126
              this%Interaction(i,j)%PotLJ126LJ126(j1, j2)%VirialProfile => this%VirialProfile
            end do
          end do
        end if
        if( this%Interaction(i,j)%N1Charge > 0 .and. this%Interaction(i,j)%N2Charge > 0 ) then
          do j1 = 1, this%Interaction(i,j)%N1Charge
            do j2 = 1, this%Interaction(i,j)%N2Charge
              this%Interaction(i,j)%PotChargeCharge(j1, j2)%VirialProfile => this%VirialProfile
            end do
          end do
        end if
        if( this%Interaction(i,j)%N1Charge > 0 .and. this%Interaction(i,j)%N2Dipole > 0 ) then
          do j1 = 1, this%Interaction(i,j)%N1Charge
            do j2 = 1, this%Interaction(i,j)%N2Dipole
              this%Interaction(i,j)%PotChargeDipole(j1, j2)%VirialProfile => this%VirialProfile
            end do
          end do
        end if
        if( this%Interaction(i,j)%N1Charge > 0 .and. this%Interaction(i,j)%N2Quadrupole > 0 ) then
          do j1 = 1, this%Interaction(i,j)%N1Charge
            do j2 = 1, this%Interaction(i,j)%N2Quadrupole
              this%Interaction(i,j)%PotChargeQuadrupole(j1, j2)%VirialProfile => this%VirialProfile
            end do
          end do
        end if
        if( this%Interaction(i,j)%N1Dipole > 0 .and. this%Interaction(i,j)%N2Charge > 0 ) then
          do j1 = 1, this%Interaction(i,j)%N1Dipole
            do j2 = 1, this%Interaction(i,j)%N2Charge
              this%Interaction(i,j)%PotDipoleCharge(j1, j2)%VirialProfile => this%VirialProfile
            end do
          end do
        end if
        if( this%Interaction(i,j)%N1Dipole > 0 .and. this%Interaction(i,j)%N2Dipole > 0 ) then
          do j1 = 1, this%Interaction(i,j)%N1Dipole
            do j2 = 1, this%Interaction(i,j)%N2Dipole
              this%Interaction(i,j)%PotDipoleDipole(j1, j2)%VirialProfile => this%VirialProfile
            end do
          end do
        end if
        if( this%Interaction(i,j)%N1Dipole > 0 .and. this%Interaction(i,j)%N2Quadrupole > 0 ) then
          do j1 = 1, this%Interaction(i,j)%N1Dipole
            do j2 = 1, this%Interaction(i,j)%N2Quadrupole
              this%Interaction(i,j)%PotDipoleQuadrupole(j1, j2)%VirialProfile => this%VirialProfile
            end do
          end do
        end if
        if( this%Interaction(i,j)%N1Quadrupole > 0 .and. this%Interaction(i,j)%N2Charge > 0 ) then
          do j1 = 1, this%Interaction(i,j)%N1Quadrupole
            do j2 = 1, this%Interaction(i,j)%N2Charge
              this%Interaction(i,j)%PotQuadrupoleCharge(j1, j2)%VirialProfile => this%VirialProfile
            end do
          end do
        end if
        if( this%Interaction(i,j)%N1Quadrupole > 0 .and. this%Interaction(i,j)%N2Dipole > 0 ) then
          do j1 = 1, this%Interaction(i,j)%N1Quadrupole
            do j2 = 1, this%Interaction(i,j)%N2Dipole
              this%Interaction(i,j)%PotQuadrupoleDipole(j1, j2)%VirialProfile => this%VirialProfile
            end do
          end do
        end if
        if( this%Interaction(i,j)%N1Quadrupole > 0 .and. this%Interaction(i,j)%N2Quadrupole > 0 ) then
          do j1 = 1, this%Interaction(i,j)%N1Quadrupole
            do j2 = 1, this%Interaction(i,j)%N2Quadrupole
              this%Interaction(i,j)%PotQuadrupoleQuadrupole(j1, j2)%VirialProfile => this%VirialProfile
            end do
          end do
        end if

      end do
    end do
#endif

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
    integer   :: i, j
#if HBOND > 0
    integer   :: k, l,stat

    allocate( this%SumHBond0(this%NComponents ), STAT = stat )
    call AllocationError( stat, 'components', this%NComponents )
    allocate( this%SumHBond1(this%NComponents, this%NComponents ), STAT = stat )
    call AllocationError( stat, 'components', this%NComponents**2 )
    allocate( this%SumHBond2(this%NComponents, this%NComponents, this%NComponents ), STAT = stat )
    call AllocationError( stat, 'components', (this%NComponents**2)*this%NComponents )
    allocate( this%SumHBond3(this%NComponents, this%NComponents, this%NComponents, this%NComponents ), STAT = stat )
    call AllocationError( stat, 'components', (this%NComponents**2)*(this%NComponents**2) )
    allocate( this%SumHBondN(this%NComponents ), STAT = stat )
    call AllocationError( stat, 'components', this%NComponents )
#endif
#if OSMOP == 2
#if HBOND == 0
    integer   :: stat
#endif

    allocate( this%SumPressureProfile(NBinsDen ), STAT = stat )
    call AllocationError( stat, 'NBinsDen', NBinsDen )
#endif

    ! Construct accumulators
    if( .not. SimulationType .eq. SecondVirialCoeff ) then

      ! 1.) Basic sums
      call Construct( this%SumPressure, .false. )
      call Construct( this%SumDensity, .false. )
      call Construct( this%SumTemperature, .false. )
      call Construct( this%SumEPot, .false. )
      call Construct( this%SumEnthalpy, .false. )
      call Construct( this%SumConfEnthalpy, .false. )
      call Construct( this%SumVolume, .false. )
      call Construct( this%SumVirial, .false. )
      call Construct( this%SumEPotInter, .false. )
      call Construct( this%SumEPotIntra, .false. )
      if (printIDF) then
        call Construct( this%SumEPotIntra_Bond, .false. )
        call Construct( this%SumEPotIntra_Angle, .false. )
        call Construct( this%SumEPotIntra_Dihedral, .false. )
        call Construct( this%SumEPotIntra_Nonbonded, .false. )
        call Construct( this%SumVirialIntra, .false. )
        call Construct( this%SumVirialInter, .false. )
      end if
      call Construct( this%SumdEpotdV, .false. )
      call Construct( this%Sumd2EpotdV2, .false. )
#if OSMOP > 0
      call Construct( this%SumOsmoticPressure, .false. )
#if OSMOP == 2
      do i = 1, NBinsDen
         call Construct( this%SumPressureProfile(i), .false. )
      end do
#endif
#endif

#if HBOND > 0
      do i = 1, this%NComponents
        call Construct( this%SumHBond0(i), .false. )
        do j = 1, this%NComponents
          call Construct( this%SumHBond1(i,j), .false. )
          do k = j, this%NComponents
            call Construct( this%SumHBond2(i,j,k), .false. )
            do l = k, this%NComponents
              call Construct( this%SumHBond3(i,j,k,l), .false. )
            end do
          end do
        end do
        call Construct( this%SumHBondN(i), .false. )
      end do
#endif

      if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
        call Construct( this%SumNPart, .false. )
      end if

      ! 2.) Combined sums
      call Construct( this%SumEPotSquared, .false. )
      call Construct( this%SumEPotV, .false. )
      call Construct( this%SumEPotVirial, .false. )
      call Construct( this%SumEnthalpySquared, .false. )
      call Construct( this%SumEnthalpyV, .false. )
      call Construct( this%SumVolumeSquared, .false. )
      call Construct( this%SumEPotCubic, .false. )
      call Construct( this%SumdEpotdVSquared, .false. )
      call Construct( this%SumEPotdEpotdV, .false. )
      call Construct( this%SumEPotSquareddEpotdV, .false. )
      call Construct( this%SumEPotdEpotdVSquared, .false. )
      call Construct( this%SumEPotd2EpotdV2, .false. )
      if( EnsembleType .eq. EnsembleTypeNVE .and. LongRange .eq. Rfield) then
        call Construct( this%SumHmU, .false. )
        call Construct( this%SumHmUm1, .false. )
        call Construct( this%SumHmUm2, .false. )
        call Construct( this%SumHmUm3, .false. )
        call Construct( this%SumHmUm1dUdV, .false. )
        call Construct( this%SumHmUm1dUdV2, .false. )
        call Construct( this%SumHmUm1d2UdV2, .false. )
        call Construct( this%SumHmUm2dUdV, .false. )
        call Construct( this%SumHmUm2dUdV2, .false. )
        call Construct( this%SumHmUm2d2UdV2, .false. )
        call Construct( this%SumHmUm3dUdV, .false. )
        call Construct( this%SumHmUm3dUdV2, .false. )
      end if

      ! 3.) Derived sums
      call Construct( this%SumBetaT, .true. )
      call Construct( this%SumdHdP, .true. )
      call Construct( this%SumdUdV, .true. )
      call Construct( this%SumCV, .true. )
      call Construct( this%SumCP, .true. )
      call Construct( this%SumAlphaP, .true. )
      if( LongRange .eq. Rfield) then
        if ( EnsembleType .eq. EnsembleTypeNVT ) then
          call Construct( this%SumA10resI, .true. )
          call Construct( this%SumA01resI, .true. )
          call Construct( this%SumA20resI, .true. )
          call Construct( this%SumA11resI, .true. )
          call Construct( this%SumA02resI, .true. )
          call Construct( this%SumA30resI, .true. )
          call Construct( this%SumA21resI, .true. )
          call Construct( this%SumA12resI, .true. )
        elseif ( EnsembleType .eq. EnsembleTypeNVE ) then
          call Construct( this%SumA10resI, .true. )
          call Construct( this%SumA01resI, .true. )
          call Construct( this%SumA20resI, .true. )
          call Construct( this%SumA11resI, .true. )
          call Construct( this%SumA02resI, .true. )
          call Construct( this%SumA30resI, .true. )
          call Construct( this%SumA21resI, .true. )
          call Construct( this%SumA12resI, .true. )
          call Construct( this%SumA10resII, .true. )
          call Construct( this%SumA01resII, .true. )
          call Construct( this%SumA20resII, .true. )
          call Construct( this%SumA11resII, .true. )
          call Construct( this%SumA02resII, .true. )
          call Construct( this%SumA30resII, .true. )
          call Construct( this%SumA21resII, .true. )
          call Construct( this%SumA12resII, .true. )
        end if
      end if

#if  TRANS == 1
!TRANSPORT_start
    ! 4.) Transport properties
    if( this%CorrfunMode ) then
      do i = 1, this%NComponents
        call Construct( this%Sumself_i(i),  .false., .true. )
        do j = 1, this%NComponents
           call Construct( this%SumOnsager(i,j), .false., .true. )
        end do
     end do

      call Construct( this%SumVisco_s, .false., .true. )
      call Construct( this%SumVisco_b, .false., .true. )
      call Construct( this%SumConduct, .false., .true. )
      call Construct( this%SumSoret,   .false., .true. )
      call Construct( this%SumEConduct,.false., .true. )

    end if
!TRANSPORT_END
#endif

! Calculation of residence times
     if (this%ResidenceTime) then
      call Construct( this%SumResidenceDuration, .false. )
      call Construct( this%SumResidencePairs, .false. )
     end if

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
    integer :: i, j
#if HBOND > 0
    integer :: k, l
#endif

    ! Destruct accumulators
    ! 1.) Basic sums
    call Destruct( this%SumPressure )
    call Destruct( this%SumDensity )
    call Destruct( this%SumTemperature )
    call Destruct( this%SumEPot )
    call Destruct( this%SumEnthalpy )
    call Destruct( this%SumConfEnthalpy )
    call Destruct( this%SumVolume )
    call Destruct( this%SumVirial )
    call Destruct( this%SumEPotInter )
    call Destruct( this%SumEPotIntra )
    if (printIDF) then
      call Destruct( this%SumEPotIntra_Bond )
      call Destruct( this%SumEPotIntra_Angle )
      call Destruct( this%SumEPotIntra_Dihedral )
      call Destruct( this%SumEPotIntra_Nonbonded )
      call Destruct( this%SumVirialIntra )
      call Destruct( this%SumVirialInter )
    end if
    call Destruct( this%SumdEpotdV )
    call Destruct( this%Sumd2EpotdV2 )
#if OSMOP > 0
      call Destruct( this%SumOsmoticPressure )
#if OSMOP == 2
      do i = 1, NBinsDen
         call Destruct( this%SumPressureProfile(i) )
      end do
#endif
#endif

#if HBOND > 0
      do i = 1, this%NComponents
        call Destruct( this%SumHBond0(i) )
        do j = 1, this%NComponents
          call Destruct( this%SumHBond1(i,j) )
          do k = j, this%NComponents
            call Destruct( this%SumHBond2(i,j,k) )
            do l = k, this%NComponents
              call Destruct( this%SumHBond3(i,j,k,l) )
            end do
          end do
        end do
        call Destruct( this%SumHBondN(i) )
      end do
#endif

    if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
      call Destruct( this%SumNPart )
    end if

    ! 2.) Combined sums
    call Destruct( this%SumEPotSquared )
    call Destruct( this%SumEPotV )
    call Destruct( this%SumEPotVirial )
    call Destruct( this%SumEnthalpySquared )
    call Destruct( this%SumEnthalpyV )
    call Destruct( this%SumVolumeSquared )
    call Destruct( this%SumEPotCubic )
    call Destruct( this%SumdEpotdVSquared )
    call Destruct( this%SumEPotdEpotdV )
    call Destruct( this%SumEPotSquareddEpotdV )
    call Destruct( this%SumEPotdEpotdVSquared )
    call Destruct( this%SumEPotd2EpotdV2 )
    if( EnsembleType .eq. EnsembleTypeNVE .and. LongRange .eq. Rfield) then
      call Destruct( this%SumHmU )
      call Destruct( this%SumHmUm1)
      call Destruct( this%SumHmUm2 )
      call Destruct( this%SumHmUm3 )
      call Destruct( this%SumHmUm1dUdV )
      call Destruct( this%SumHmUm1dUdV2 )
      call Destruct( this%SumHmUm1d2UdV2 )
      call Destruct( this%SumHmUm2dUdV )
      call Destruct( this%SumHmUm2dUdV2 )
      call Destruct( this%SumHmUm2d2UdV2 )
      call Destruct( this%SumHmUm3dUdV )
      call Destruct( this%SumHmUm3dUdV2 )
    end if

    ! 3.) Derived sums
    call Destruct( this%SumBetaT )
    call Destruct( this%SumdHdP )
    call Destruct( this%SumdUdV )
    call Destruct( this%SumCV )
    call Destruct( this%SumCP )
    call Destruct( this%SumAlphaP )
    if( LongRange .eq. Rfield) then
      if ( EnsembleType .eq. EnsembleTypeNVT ) then
        call Destruct( this%SumA10resI )
        call Destruct( this%SumA01resI )
        call Destruct( this%SumA20resI )
        call Destruct( this%SumA11resI )
        call Destruct( this%SumA02resI )
        call Destruct( this%SumA30resI )
        call Destruct( this%SumA21resI )
        call Destruct( this%SumA12resI )
      elseif ( EnsembleType .eq. EnsembleTypeNVE ) then
        call Destruct( this%SumA10resI )
        call Destruct( this%SumA01resI )
        call Destruct( this%SumA20resI )
        call Destruct( this%SumA11resI )
        call Destruct( this%SumA02resI )
        call Destruct( this%SumA30resI )
        call Destruct( this%SumA21resI )
        call Destruct( this%SumA12resI )
        call Destruct( this%SumA10resII )
        call Destruct( this%SumA01resII )
        call Destruct( this%SumA20resII )
        call Destruct( this%SumA11resII )
        call Destruct( this%SumA02resII )
        call Destruct( this%SumA30resII )
        call Destruct( this%SumA21resII )
        call Destruct( this%SumA12resII )
      end if
    end if

#if  TRANS == 1
!TRANSPORT_start
    if( this%CorrfunMode ) then

      do i = 1, this%NComponents
        call Destruct( this%Sumself_i(i) )
        do j = 1, this%NComponents
          call Destruct( this%SumOnsager(i,j) )
        end do
      end do

      call Destruct( this%SumVisco_s )
      call Destruct( this%SumVisco_b )
      call Destruct( this%SumConduct )
      call Destruct( this%SumSoret )
      call Destruct( this%SumEConduct )

    end if

!TRANSPORT_END
#endif

! Calculation of residence times
    if ( this%ResidenceTime ) then
      call Destruct( this%SumResidenceDuration )
      call Destruct( this%SumResidencePairs )
    end if

    do i = 1, this%NRealComponents
      call DestroyAccumulators( this%Component(i) )
    end do

#if OSMOP == 2
    if( associated( this%SumPressureProfile ) ) then
      deallocate( this%SumPressureProfile )
    end if
#endif

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
    integer                   :: i,t
    integer                   :: k, j,  NUnitMax

    ! Adjust number of cells to cube of integer
    if( this%NPart < NPartInCell ) this%NPart = NPartInCell

    ! Set maximum number of particles
    if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
      this%NPartMax = 2 * this%NPart
! Max. number of particles of component i in a fluctuating state 
      this%NPartMaxFluct = 1
    else
      this%NPartMax = this%NPart
! Max. number of particles of component i in a fluctuating state 
      this%NPartMaxFluct = 1
    end if

    ! Calculate Max Number of Units in Component
    NUnitMax = 0
    do i = 1, this%NComponents
        k = this%Component(i)%Molecule%NUnit
        if (k>NUnitMax) NUnitMax = k
    end do
    this%NUnitMax=NUnitMax

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
      this%Component(this%NRealComponents)%NPart = this%Component(this%NRealComponents)%NPart - pc%NPart
    end do

    ! Setting Particles for ThermoInt
    t = this%NRealComponents+1
    do i = 1, this%NRealComponents
      pc => this%Component(i)
      if ( pc%ChemPotMethod == ChemPotMethodThermoInt ) then
        if ( pc%NPart >= 1) then
          pc%NPart = pc%NPart - 1
          this%Component(t)%NPart = 1
          t = t+1
        else
          write (ErrorBuffer,'("At least one particle of ", A, " is required for Thermodynamic Integration" ,T45)' ) trim( pc%PotModFileName )
          call Error
        end if
      end if
    end do

    ! Set mole fractions according to real number of particles
    ! and calculate number of degrees of freedom
    this%NDFTran = 0
    this%NDFRot = 0
    this%constrNDF = 0
    do i = 1, this%NComponents
      pc => this%Component(i)
      pc%Fraction = real( pc%NPart, RK ) / real( this%NPart, RK )
      pc%NDFTran = pc%NPart * pc%Molecule%NUnit * 3  ! all unit*3
      pc%NDFRot=0
      do j = 1, pc%Molecule%NUnit
        pc%NDFRot = pc%NDFRot + pc%Molecule%Unit(j)%NDFRot ! for one molecule
      end do
     ! Inner Degrees of Freedom of one particle
      if (UseIntDegFreed) then
        pc%Molecule%NDF = pc%Molecule%NUnit * 3
        pc%Molecule%NDF = pc%Molecule%NDF + pc%NDFRot
        if ( Shake > 0 ) then
          this%constrNDF = this%constrNDF + pc%NPart*pc%Molecule%NBond
        end if
      end if
      pc%NDFRot = pc%NPart * pc%NDFRot
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

!      if( pc%NTest > 0 ) pc%NTest = 1 + (pc%NTest - 1) / NProcs
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
      write( IOBuffer, '("Mole fraction of ", A, ":",T45, F6.3)' ) trim( pc%PotModFileName ), pc%Fraction
      call LogWrite
      write( IOBuffer, '(T10,"->  number of particles: ",T36, I11)' ) pc%NPart
      call LogWrite

      if( pc%NTest > 0 ) then
        write( IOBuffer, '("Number of test particles:",T36, I11)' ) pc%NTest
        call LogWrite
      end if
    end do

    this%NUnitTotal=0
    do i = 1, this%NComponents
        this%NUnitTotal = this%NUnitTotal + &
&                this%Component(i)%NPart*this%Component(i)%Molecule%NUnit
    end do

  end subroutine TEnsemble_CalculateNPart



!==============================================================!
!  Subroutine TEnsemble_LongRangeCheck                         !
!==============================================================!

   subroutine TEnsemble_LongRangeCheck ( this )

   implicit none
   type(TEnsemble)            :: this
   integer                    :: i
   real(RK)                   :: q

! Calculation
   q = 0._RK
   do i=1,this%NComponents
     call LongRangeCheck ( this%Component(i), q)
   end do

! Error Analysis
   if (abs(q) .ge. 1e-1) then
     write (ErrorBuffer,'("You have a non-neutral system.\n NetCharge normed = ", &
&                 F20.10, "\n Conflicts arise applying long range corrections")') q
     call Error
   end if

   end subroutine TEnsemble_LongRangeCheck



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

      if (LongRange .eq. Ewald) then
        this%Kappa = this%KappaL / this%BoxLength
        do i=1,this%NComponents
          do j=1,this%NComponents
            this%Interaction(i,j)%Kappa = this%Kappa
          end do
        end do
        call EwaldSelfTerm(this)

#if SPME > 0
      else if (LongRange .eq. PME ) then
        this%Kappa = this%KappaL / this%BoxLength
        do i=1,this%NComponents
          do j=1,this%NComponents
            this%Interaction(i,j)%Kappa = this%Kappa
          end do
        end do
        call PMESelfTermMC ( this )
#endif

      else if (LongRange .eq. Rodgers ) then
        this%EPotCorrRF = this%EPotCorrRFVol*this%Volume0 + this%EPotCorrRFPart
      end if

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
    integer                   :: i, j
    type(TComponent), pointer :: pc

    ! Set mole fractions according to real number of particles
    ! and calculate number of degrees of freedom
    this%NDFTran = 0
    this%NDFRot = 0
    this%constrNDF = 0
    do i = 1, this%NComponents
      pc => this%Component(i)
      this%Component(i)%Fraction = real( pc%NPart, RK ) / real( this%NPart, RK )
      pc%NDFTran = pc%NPart * pc%Molecule%NUnit * 3
      pc%NDFRot=0
      do j = 1, pc%Molecule%NUnit
        pc%NDFRot = pc%NDFRot + pc%Molecule%Unit(j)%NDFRot ! for one molecule
      end do
      if ( Shake > 0 .and. UseIntDegFreed) then
          this%constrNDF = this%constrNDF + pc%NPart*pc%Molecule%NBond
      end if
      pc%NDFRot = pc%NPart * pc%NDFRot
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
      if( pc%Molecule%NLJ126 > this%NLJ126Max ) this%NLJ126Max = pc%Molecule%NLJ126
      if( pc%Molecule%NCharge > this%NChargeMax ) this%NChargeMax = pc%Molecule%NCharge
      if( pc%Molecule%NDipole > this%NDipoleMax ) this%NDipoleMax = pc%Molecule%NDipole
      if( pc%Molecule%NQuadrupole > this%NQuadrupoleMax ) this%NQuadrupoleMax = pc%Molecule%NQuadrupole
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
    integer :: number
#if TRANS ==1
    integer :: NPart3
    integer :: NComp2
#endif


    ! Nullify pointers
    nullify( this%P0Test )
    nullify( this%Q0Test )
    nullify( this%EPotTest )
    nullify( this%BiasedPartners )
    nullify( this%RDFValue )
    nullify( this%RDFVSchale )

    ! Allocate scale coefficients for sigma and epsilon
    allocate( this%ScaleSigma(this%NComponents, this%NComponents), STAT = stat )
    call AllocationError( stat, 'components', this%NComponents )
    allocate( this%ScaleEpsilon(this%NComponents, this%NComponents), STAT = stat )
    call AllocationError( stat, 'components', this%NComponents )
    allocate( this%BiasedPartners(this%NPartMax), STAT = stat )
    call AllocationError( stat, 'NPartMax', this%NPartMax )

    ! Allocate RDF arrays
    if( RDFUpdateFrequency > 0 ) then
      allocate( this%RDFVSchale(RDFNumberShells), STAT = stat )
      call AllocationError( stat, 'components', RDFNumberShells )
      allocate( this%RDFValue(RDFNumberShells), STAT = stat )
      call AllocationError( stat, 'components', RDFNumberShells )    
    endif

    ! Allocate test particles
    if( this%NTestMax > 0 ) then
      allocate( this%P0Test( this%NTestMax, 3, this%NUnitMax ), STAT = stat )
      call AllocationError( stat, 'test particles', this%NTestMax )
      allocate( this%Q0Test( this%NTestMax, 4, this%NUnitMax ), STAT = stat )
      call AllocationError( stat, 'test particles', this%NTestMax )
      allocate( this%EPotTest( this%NTestMax ), STAT = stat )
      call AllocationError( stat, 'test particles', this%NTestMax )
      call LogWriteBlank
      write( IOBuffer, '("Memory for test particles allocated successfully")' )
      call LogWrite
    end if

    ! Allocate components
    if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. &
&       SimulationType .eq. Gibbs .or. SimulationType .eq. SecondVirialCoeff ) then       

       do i = 1, this%NComponents
         this%Component(i)%NPartMax => this%NPartMax
         this%Component(i)%NUnitMax => this%NUnitMax

         if( this%Component(i)%NTest > 0 ) then
           this%Component(i)%P0Test => this%P0Test
           this%Component(i)%Q0Test => this%Q0Test
         end if

         this%Component(i)%BoxLength => this%BoxLength
         call Allocate( this%Component(i) )
       end do

    else
       do i = 1, this%NComponents
         if (i .le. this%NRealComponents) then
           this%Component(i)%NPartMax => this%Component(i)%NPart
           this%Component(i)%NUnitMax => this%Component(i)%Molecule%NUnit

         else 
           this%Component(i)%NPartMax => this%NPartMaxFluct
           this%Component(i)%NUnitMax => this%NUnitMax
         end if

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

#if HBOND > 0
    !Allocate H-Bond counters
    allocate( this%NHBond0( this%NComponents ), STAT = stat )
    call AllocationError( stat, 'components', this%NComponents )
    allocate( this%NHBond1( this%NComponents, this%NComponents ), STAT = stat )
    call AllocationError( stat, 'components', this%NComponents**2 )
    allocate( this%NHBond2( this%NComponents, this%NComponents, this%NComponents ), STAT = stat )
    call AllocationError( stat, 'components', this%NComponents**2*this%NComponents )
    allocate( this%NHBond3( this%NComponents, this%NComponents, this%NComponents, this%NComponents ), STAT = stat )
    call AllocationError( stat, 'components', this%NComponents**2*this%NComponents**2 )
    allocate( this%NHBondN( this%NComponents ), STAT = stat )
    call AllocationError( stat, 'components', this%NComponents )
#endif

#if OSMOP == 2
    !Allocate pressure in the bins of the pressureprofile
    allocate( this%VirialProfile( NBinsDen ), STAT = stat )
    call AllocationError( stat, 'viral profile', NBinsDen )
    this%VirialProfile(:) = 0._RK
    allocate( this%PressureProfile( NBinsDen ), STAT = stat )
    call AllocationError( stat, 'pressure profile', NBinsDen )
    this%PressureProfile(:) = 0._RK
#endif

#if  TRANS == 1
!TRANSPORT_start
      NPart3 = 3*this%NPart
      NComp2 = this%NComponents*this%NComponents

    ! Allocate correlation fucntions
     if( this%CorrfunMode ) then

      allocate( this%cf_vs(this%NCorr), STAT = stat )
      call AllocationError( stat, 'viscosity_shear_cf_vs', this%NCorr )

      allocate( this%average_cf_vs(this%NCorr), STAT = stat )
      call AllocationError( stat, 'viscosity_shear_cf_vs', this%NCorr )

      allocate( this%cf_vb(this%NCorr), STAT = stat )
      call AllocationError( stat, 'viscosity_bulk_cf_vb', this%NCorr )

      allocate( this%average_cf_vb(this%NCorr), STAT = stat )
      call AllocationError( stat, 'viscosity_bulk_cf_vb', this%NCorr )

      allocate( this%cf_c(this%NCorr), STAT = stat )
      call AllocationError( stat, 'conductivity_cf_c', this%NCorr )

      allocate( this%average_cf_c(this%NCorr), STAT = stat )
      call AllocationError( stat, 'conductivity_cf_c', this%NCorr )

      allocate( this%cf_ec(this%NCorr), STAT = stat )
      call AllocationError( stat, 'conductivity_cf_ec', this%NCorr )

      allocate( this%average_cf_ec(this%NCorr), STAT = stat )
      call AllocationError( stat, 'conductivity_cf_ec', this%NCorr )

      allocate( this%cf_d( this%NComponents, this%NCorr), STAT = stat )
      call AllocationError( stat, 'self_diffusion', this%NCorr )

      allocate( this%average_cf_d( this%NComponents, this%NCorr), STAT = stat )
      call AllocationError( stat, 'self_diffusion', this%NCorr )

      allocate( this%cf_db( this%NCorr), STAT = stat )
      call AllocationError( stat, 'binary_diffusion', this%NCorr )

      allocate( this%average_cf_db( this%NCorr), STAT = stat )
      call AllocationError( stat, 'binary_diffusion', this%NCorr )

      allocate( this%cf_soret( this%NCorr), STAT = stat )
      call AllocationError( stat, 'thermal_diffusion', this%NCorr )

      allocate( this%average_cf_soret( this%NCorr), STAT = stat )
      call AllocationError( stat, 'thermal_diffusion', this%NCorr )

      allocate( this%lamda( NComp2, this%NCorr ), STAT = stat )
      call AllocationError( stat, 'onsager_coefficient', this%NCorr )

      allocate( this%average_lamda( NComp2, this%NCorr ), STAT = stat )
      call AllocationError( stat, 'onsager_coefficient', this%NCorr )

      allocate( this%sinte_i( this%NComponents, this%NCorr), STAT = stat )
      call AllocationError( stat, 'self_diffusion_integrated', this%NCorr )

      allocate( this%average_sinte_i( this%NComponents, this%NCorr), STAT = stat )
      call AllocationError( stat, 'self_diffusion_integrated', this%NCorr )

      allocate( this%sinte_db( this%NCorr), STAT = stat )
      call AllocationError( stat, 'mutual_diffusion integrated', this%NCorr )

      allocate( this%average_sinte_db( this%NCorr), STAT = stat )
      call AllocationError( stat, 'mutual_diffusion integrated', this%NCorr )

      allocate( this%sinte_soret( this%NCorr), STAT = stat )
      call AllocationError( stat, 'thermal_diffusion integrated', this%NCorr )

      allocate( this%average_sinte_soret( this%NCorr), STAT = stat )
      call AllocationError( stat, 'thermal_diffusion integrated', this%NCorr )

      allocate( this%sinte_lamda( NComp2, this%NCorr), STAT = stat )
      call AllocationError( stat, 'mutual diffusion integrated', this%NCorr )

      allocate( this%average_sinte_lamda( NComp2, this%NCorr), STAT = stat )
      call AllocationError( stat, 'mutual diffusion integrated', this%NCorr )      

      allocate( this%sinte_vs( this%NCorr), STAT = stat )
      call AllocationError( stat, 'shear_viscosity_integrated', this%NCorr )

      allocate( this%average_sinte_vs( this%NCorr), STAT = stat )
      call AllocationError( stat, 'shear_viscosity_integrated', this%NCorr )

      allocate( this%sinte_vb( this%NCorr), STAT = stat )
      call AllocationError( stat, 'bulk_viscosity_integrated', this%NCorr )

      allocate( this%average_sinte_vb( this%NCorr), STAT = stat )
      call AllocationError( stat, 'bulk_viscosity_integrated', this%NCorr )

      allocate( this%sinte_c( this%NCorr), STAT = stat )
      call AllocationError( stat, 'Thermal_conductivity_integrated', this%NCorr )

      allocate( this%average_sinte_c( this%NCorr), STAT = stat )
      call AllocationError( stat, 'Thermal_conductivity_integrated', this%NCorr )

      allocate( this%sinte_ec( this%NCorr), STAT = stat )
      call AllocationError( stat, 'Electric_conductivity_integrated', this%NCorr )

      allocate( this%average_sinte_ec( this%NCorr), STAT = stat )
      call AllocationError( stat, 'Electric_conductivity_integrated', this%NCorr )

      allocate( this%a( NPart3, this%NCorr), STAT = stat  )
      call AllocationError( stat, 'diffusion_matrix', NPart3 )

      allocate( this%A_SpanCF( NPart3, this%NSpanCF), STAT = stat  )
      call AllocationError( stat, 'diffusion_matrix', NPart3 )

      allocate( this%vsk(this%NCorr, 3), STAT = stat  )
      call AllocationError( stat, 'vsk', this%NPart )

      allocate( this%vsp(this%NCorr, 3), STAT = stat  )
      call AllocationError( stat, 'vsp', this%NPart )

      allocate( this%vbk(this%NCorr, 3), STAT = stat  )
      call AllocationError( stat, 'vbk', this%NPart )

      allocate( this%vbp(this%NCorr, 3), STAT = stat  )
      call AllocationError( stat, 'vbp', this%NPart )

      allocate( this%vckt(this%NCorr, 3), STAT = stat  )
      call AllocationError( stat, 'vckt', this%NPart )

      allocate( this%vckr(this%NCorr, 3), STAT = stat  )
      call AllocationError( stat, 'vckr', this%NPart )

      allocate( this%vcpt(this%NCorr, 3), STAT = stat  )
      call AllocationError( stat, 'vcpt', this%NPart )

      allocate( this%vcpr(this%NCorr, 3), STAT = stat  )
      call AllocationError( stat, 'vcpr', this%NPart )

      allocate( this%vcmt(this%NCorr, 3), STAT = stat  )
      call AllocationError( stat, 'vcmt', this%NPart )

      allocate( this%selfd_i(this%NComponents), STAT = stat  )
      call AllocationError( stat, 'selfd_i', this%NComponents )

      allocate( this%Sumself_i(this%NComponents), STAT = stat  )
      call AllocationError( stat, 'Sumselfd_i', this%NComponents )

       allocate( this%Onsager(this%NComponents,this%NComponents), STAT = stat  )
      call AllocationError( stat, 'Onsager', this%NComponents )

      allocate( this%SumOnsager(this%NComponents,this%NComponents), STAT = stat  )
      call AllocationError( stat, 'SumOnsager', this%NComponents )




      ! Set correlation-fucntion vectors
      this%cf_d(:,:)      = 0._RK
      this%cf_db(:)       = 0._RK
      this%cf_soret(:)    = 0._RK
      this%lamda(:,:)     = 0._RK
      this%cf_vs(:)       = 0._RK
      this%cf_vb(:)       = 0._RK
      this%cf_c(:)        = 0._RK
      this%cf_ec(:)       = 0._RK

      this%average_cf_d(:,:)   = 0._RK
      this%average_cf_db(:)    = 0._RK
      this%average_cf_soret(:) = 0._RK
      this%average_lamda(:,:)  = 0._RK
      this%average_cf_vs(:)    = 0._RK
      this%average_cf_vb(:)    = 0._RK
      this%average_cf_c(:)     = 0._RK
      this%average_cf_ec(:)    = 0._RK
  
      this%a(:,:)           = 0._RK
      this%A_SpanCF(:,:)    = 0._RK
      
      this%sinte_i(:,:)     = 0._RK
      this%sinte_lamda(:,:) = 0._RK
      this%sinte_db (:)     = 0._RK
      this%sinte_soret(:)   = 0._RK
      this%sinte_vs(:)      = 0._RK
      this%sinte_vb(:)      = 0._RK
      this%sinte_c(:)       = 0._RK
      this%sinte_ec(:)      = 0._RK
      
      this%average_sinte_i(:,:)     = 0._RK
      this%average_sinte_db (:)     = 0._RK
      this%average_sinte_soret(:)   = 0._RK
      this%average_sinte_lamda(:,:) = 0._RK
      this%average_sinte_vs(:)      = 0._RK
      this%average_sinte_vb(:)      = 0._RK
      this%average_sinte_c(:)       = 0._RK
      this%average_sinte_ec(:)      = 0._RK

      this%selfd_i(:)  = 0._RK
      this%vsk(:,:)    = 0._RK
      this%vsp(:,:)    = 0._RK
      this%vbk(:,:)    = 0._RK
      this%vbp(:,:)    = 0._RK
      this%vckt(:,:)   = 0._RK
      this%vckr(:,:)   = 0._RK
      this%vcpt(:,:)   = 0._RK
      this%vcpr(:,:)   = 0._RK
      this%vcmt(:,:)   = 0._RK

      this%sc(:) = 0._RK
      this%sp(:) = 0._RK
    end if
    !TRANSPORT_END
#endif

! Calculation of residence times
    if ( this%ResidenceTime ) then

      nullify( this%CompPair )
      nullify( this%CompPair_Old )
      nullify( this%ResidTimesStart )
      nullify( this%ResidTimesStart_Old )
      nullify( this%ResidPairsCem )

      allocate(this%CompPair(this%Component(this%ResidComp1)%NPart*10,2),STAT=stat)
      call AllocationError( stat, 'CompPair' )
      allocate(this%CompPair_Old(this%Component(this%ResidComp1)%NPart*10,2),STAT=stat)
      call AllocationError( stat, 'CompPair_Old' )
      allocate(this%ResidTimesStart(this%Component(this%ResidComp1)%NPart*10),STAT=stat)
      call AllocationError( stat, 'ResidTimesStart' )
      allocate(this%ResidTimesStart_Old(this%Component(this%ResidComp1)%NPart*10),STAT=stat)
      call AllocationError( stat, 'ResidTimesStart_Old' )
      allocate(this%ResidPairsCem(this%Component(this%ResidComp1)%NPart*10,4),STAT=stat)
      call AllocationError( stat, 'ResidPairsCem' )
    end if

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

    if( associated( this%BiasedPartners ) ) then
      deallocate( this%BiasedPartners )
    end if

    if( associated( this%RDFVSchale ) ) then
      deallocate( this%RDFVSchale )
    end if

    if( associated( this%RDFValue ) ) then
      deallocate( this%RDFValue )
    end if    

#if  TRANS == 1
!TRANSPORT_start
    ! Deallocate arrays for correlation fucntions

    if( associated( this%cf_d  ) )   then
      deallocate( this%cf_d  )
    end if
   
     if( associated( this%average_cf_d  ) )   then
      deallocate( this%average_cf_d  )
    end if

     if( associated( this%cf_db ) )   then
       deallocate( this%cf_db )
     end if

     if( associated( this%average_cf_db ) )   then
       deallocate( this%average_cf_db )
     end if

    if( associated( this%cf_soret ) )   then
       deallocate( this%cf_soret )
     end if

    if( associated( this%average_cf_soret ) )   then
       deallocate( this%average_cf_soret )
     end if

    if( associated( this%cf_vs ) )   then
      deallocate( this%cf_vs )
    end if

    if( associated( this%average_cf_vs ) )   then
      deallocate( this%average_cf_vs )
    end if

    if( associated( this%cf_vb ) )   then
      deallocate( this%cf_vb )
    end if

    if( associated( this%average_cf_vb ) )   then
      deallocate( this%average_cf_vb )
    end if

    if( associated( this%cf_c  ) )   then
      deallocate( this%cf_c  )
    end if

    if( associated( this%average_cf_c  ) )   then
      deallocate( this%average_cf_c  )
    end if

    if( associated( this%cf_ec  ) )   then
      deallocate( this%cf_ec  )
    end if

    if( associated( this%average_cf_ec  ) )   then
      deallocate( this%average_cf_ec  )
    end if

    if( associated( this%lamda ) )   then
      deallocate( this%lamda  )
    end if

    if( associated( this%average_lamda ) )   then
      deallocate( this%average_lamda  )
    end if

    if( associated( this%sinte_i) )  then
      deallocate( this%sinte_i  )
    end if

    if( associated( this%average_sinte_i) )  then
      deallocate( this%average_sinte_i  )
    end if

     if( associated( this%sinte_db) ) then
       deallocate( this%sinte_db )
    end if

    if( associated( this%average_sinte_db) ) then
       deallocate( this%average_sinte_db )
    end if

    if( associated( this%average_sinte_soret) ) then
       deallocate( this%average_sinte_soret )
    end if

    if( associated( this%sinte_vs) ) then
      deallocate( this%sinte_vs )
    end if

    if( associated( this%average_sinte_vs) ) then
      deallocate( this%average_sinte_vs )
    end if
    
    if( associated( this%sinte_vb) ) then
      deallocate( this%sinte_vb )
    end if

    if( associated( this%average_sinte_vb) ) then
      deallocate( this%average_sinte_vb )
    end if

    if( associated( this%sinte_c) ) then
      deallocate( this%sinte_c )
    end if

    if( associated( this%average_sinte_c) ) then
      deallocate( this%average_sinte_c )
    end if

    if( associated( this%sinte_ec) ) then
      deallocate( this%sinte_ec )
    end if

    if( associated( this%average_sinte_c)  ) then
      deallocate( this%average_sinte_c )
    end if

    if( associated( this%a )  )   then
      deallocate( this%a  )
    end if

    if( associated( this%A_SpanCF )  )   then
      deallocate( this%A_SpanCF  )
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

    if( associated( this%vcmt ) ) then
      deallocate( this%vcmt )
    end if

    if( associated( this%selfd_i ) ) then
      deallocate( this%selfd_i )
    end if

    if( associated( this%Sumself_i ) ) then
      deallocate( this%Sumself_i )
    end if

    if( associated( this%Onsager ) ) then
      deallocate( this%Onsager )
    end if

    if( associated( this%SumOnsager ) ) then
      deallocate( this%SumOnsager )
    end if

!TRANSPORT_END
#endif

! Calculation of residence times
    if ( this%ResidenceTime ) then
      if( associated( this%CompPair ) ) then
        deallocate( this%CompPair )
      end if

      if( associated( this%CompPair_Old ) ) then
        deallocate( this%CompPair_Old )
      end if

      if( associated( this%ResidTimesStart ) ) then
        deallocate( this%ResidTimesStart )
      end if

      if( associated( this%ResidTimesStart_Old ) ) then
        deallocate( this%ResidTimesStart_Old )
      end if

      if( associated( this%ResidPairsCem ) ) then
        deallocate( this%ResidPairsCem )
      end if
    end if

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
    real(RK)                        :: fac
    real(RK)                        :: fac_neutral, fac_charge1, fac_charge2
    real(RK)                        :: totcharge, sumcharge
    real(RK)                        :: drx, dry, drz, dr2

    ! Assign local variables
    NPartInv = 1._RK / this%NPart

    if (LongRange .eq. ExtRField) then
      fac = this%DebyeLen*this%RCutoffDipoleDipole
      RFConst = -1._RK / this%RCutoffDipoleDipole**3 * ((this%RFEpsilon - 1._RK)*(1._RK+fac) + 0.5*this%RFEpsilon*(fac)**2) &
&               / ( (2._RK * this%RFEpsilon+1._RK)*(1._RK+fac) + this%RFEpsilon*(fac)**2 )
    else
      RFConst = -1._RK / this%RCutoffDipoleDipole**3 * (this%RFEpsilon - 1._RK) / (2._RK * this%RFEpsilon + 1._RK)
    endif

    ! Set maximum cutoff radius
    this%NRCutoffMax = 0
    this%RCutoffMax2 = 0._RK

    ! Zero long-range corrections
    this%EPotCorrLJ   = 0._RK
    this%EPotCorrRF   = 0._RK
    this%VirialCorrLJ = 0._RK
    this%VirialCorrRF = 0._RK
    this%d2EpotdV2CorrLJ   = 0._RK

    do i1 = 1, this%NComponents
      this%Component(i1)%EPotTestCorrLJ = 0._RK
    end do

    ! Calculate Lennard-Jones long-range corrections
    if( this%NLJ126Max > 0 ) then
      do i1 = 1, this%NComponents
        do i2 = 1, this%NComponents
          Scale = this%Component(i1)%NPart * this%Component(i2)%NPart * NPartInv
          do j1 = 1, this%Component(i1)%Molecule%NLJ126
            do j2 = 1, this%Component(i2)%Molecule%NLJ126
              plj => this%Interaction(i1, i2)%PotLJ126LJ126(j1, j2)
              this%EPotCorrLJ = this%EPotCorrLJ + Scale * plj%EPotCorr
              this%VirialCorrLJ = this%VirialCorrLJ + Scale * plj%VirialCorr
              this%d2EpotdV2CorrLJ = this%d2EpotdV2CorrLJ + Scale * plj%d2EpotdV2Corr
              this%Component(i1)%EPotTestCorrLJ = this%Component(i1)%EPotTestCorrLJ &
&                 + this%Component(i2)%Fraction * plj%EPotTestCorr
              this%RCutoffMax2 = max( this%RCutoffMax2, 2._RK * sqrt( plj%RCutoffSquared ) )
            end do
          end do
        end do
      end do

      this%EPotCorrLJ = this%EPotCorrLJ / NProcs
      this%VirialCorrLJ = this%VirialCorrLJ / NProcs
      this%d2EpotdV2CorrLJ = this%d2EpotdV2CorrLJ / NProcs
    end if

    ! Calculate electrostatic long-range corrections
    ! This is the self term of the reaction field
    if( (this%NChargeMax > 0).or.(this%NDipoleMax > 0) ) then

      if (LongRange .ne. Rodgers ) then
        do i1 = 1, this%NComponents
          pc => this%Component(i1)
          pc%EPotTestCorrRF = 0._RK
          do j1 = 1, pc%Molecule%NUnit
            this%EPotCorrRF = this%EPotCorrRF + pc%Molecule%Unit(j1)%MueSquared * pc%NPart
            pc%EPotTestCorrRF = pc%EPotTestCorrRF + pc%Molecule%Unit(j1)%MueSquared * 2._RK * RFConst
          end do
        end do
        this%EPotCorrRF = this%EPotCorrRF * RFConst / NProcs

      else ! Rodgers
        this%Kappa = UnitLength * this%KappaL / Angstroem  ! = 1/sigma* aus Paper
        fac_charge1 = 1._RK /  sqrt(Pi) * this%Kappa
        fac_charge2 = fac_charge1 /  3._RK * this%Kappa*this%Kappa
        fac_neutral = 2._RK *fac_charge1 * this%Kappa*this%Kappa

        this%EPotCorrRFVol = fac_neutral/(4._RK*PI) * (1._RK-1._RK/this%RFEpsilon) * this%RefTemperature / NProcs
        this%EPotCorrRFPart = 0._RK

        do i1 = 1, this%NComponents
          pc => this%Component(i1)

          if ( pc%charged ) then

            ! Calculate total charge * distances
            totcharge = 0._RK
            sumcharge = 0._RK
            do j1 = 1, pc%Molecule%NCharge
              do j2 = 1, pc%Molecule%NCharge
                drx = pc%Molecule%SiteCharge(j1)%r(1) - pc%Molecule%SiteCharge(j2)%r(1)
                dry = pc%Molecule%SiteCharge(j1)%r(2) - pc%Molecule%SiteCharge(j2)%r(2)
                drz = pc%Molecule%SiteCharge(j1)%r(3) - pc%Molecule%SiteCharge(j2)%r(3)
                dr2  = drx*drx + dry*dry + drz*drz
                totcharge = totcharge + pc%Molecule%SiteCharge(j1)%e * pc%Molecule%SiteCharge(j2)%e * dr2
              end do
              sumcharge = sumcharge + pc%Molecule%SiteCharge(j1)%e
            end do

            this%EPotCorrRFPart = this%EPotCorrRFPart - fac_charge1 * sumcharge*sumcharge * pc%NPart / NProcs + &
&                   fac_charge2 * totcharge * pc%NPart / NProcs

          else ! charged
            do j1 = 1, pc%Molecule%NUnit
              this%EPotCorrRFPart = this%EPotCorrRFPart - fac_neutral / 3._RK * pc%Molecule%Unit(j1)%MueSquared * pc%NPart / NProcs
              this%VirialCorrRF   = this%VirialCorrRF - fac_neutral / (4._RK*PI) * this%RefTemperature * &
&                       (1._RK-1._RK/this%RFEpsilon) / NProcs 
            end do
          end if

        end do
        this%EPotCorrRF = this%EPotCorrRFVol* this%Volume0 + this%EPotCorrRFPart
      end if
    end if



#if MPI_VER >0
  if ( SimulationType .eq. MonteCarlo .and. (.not.(Equilibration .and. CommonEqui)) ) then
      if( this%NLJ126Max > 0 ) then
        this%EPotCorrLJ = this%EPotCorrLJ * NProcs
        this%VirialCorrLJ = this%VirialCorrLJ * NProcs
        this%d2EpotdV2CorrLJ = this%d2EpotdV2CorrLJ * NProcs
      endif

      if( (this%NChargeMax > 0).or.(this%NDipoleMax > 0) ) then
        this%EPotCorrRF = this%EPotCorrRF * NProcs
      endif

      if (LongRange .eq. Rodgers ) then
        this%EPotCorrRFPart = this%EPotCorrRFPart * NProcs
        this%EPotCorrRFVol  = this%EPotCorrRFVol  * NProcs
      end if
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
#if OSMOP > 0
    integer                   :: Npermeable, nc1, nm1
    real(RK), dimension(3)    :: swap
#endif

    NCells = ceiling( real( this%NPart, RK ) / real( NPartInCell, RK ) )
    NCells1dim = ceiling( NCells**Third )

#if OSMOP > 0
    ! Count permeable Particles and check NCells
    if (SimulationType .eq. MolecularDynamics) then
      Npermeable = 0
      do i = 1, this%NComponents
        if (this%Component(i)%permeable) then
          Npermeable = Npermeable + this%Component(i)%NPart
        end if
      end do

      do while ( (this%NPart-Npermeable) .gt. 0.5_RK*NCells*NPartInCell )
        NCells = 1.2_RK * ( ceiling( real( this%NPart-Npermeable, RK ) / real( NPartInCell, RK ) ) )
        NCells1dim = ceiling( NCells**Third )
      end do
    end if
#endif

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
    write( IOBuffer, '(T10,"FCC lattice: ",I3," *",I4,"x",I4,"x",I4," cells")' ) NPartInCell, ( NCells1dim(i), i=1,3 )
    call LogWrite
    write( IOBuffer, '(T10, "with",I3," molecules/cell")' ) NPartInCell
    call LogWrite

    ! Set all positions
loop:do l = 1, NPartInCell
      do i = 1, NCells1dim(1)
        do j = 1, NCells1dim(2)
          do k = 1, NCells1dim(3)
            nc = select_component( comp )
! #if OSMOP > 0
!             if (SimulationType .eq. MolecularDynamics ) then
!               if ( Npermeable == 0 .and. i < ceiling((NCells1dim(1)+1)/2._RK) ) then
!                 comp(nc)=comp(nc)+1
!                 cycle xloop
!               end if
!               if ( i < ceiling((NCells1dim(1)+1)/2._RK) ) then
!                 do while (.not. this%Component(nc)%permeable ) ! set permeable particles as long as some are still to be set
!                   comp(nc)=comp(nc)+1
!                   nc = select_component( comp )
!                 end do
!               end if
!               if ( this%Component(nc)%permeable ) Npermeable = Npermeable - 1
!             end if
! #endif
            nm = comp(nc) + 1
            pc => this%Component(nc)
            pc%Pm0(nm, 1) = xl(1) * (CellX(l) + i - 1)
            pc%Pm0(nm, 2) = xl(2) * (CellY(l) + j - 1)
            pc%Pm0(nm, 3) = xl(3) * (CellZ(l) + k - 1)
            n = n + 1

            if( n == this%NPart ) exit loop
          end do
        end do
      end do
    end do loop

#if OSMOP > 0
    do i = 1, this%NComponents
      if (.not. this%Component(i)%permeable) then
        do j = 1, this%Component(i)%NPart
          if ( this%Component(i)%Pm0(j,1) .le. 0.25_RK .or. this%Component(i)%Pm0(j,1) .ge. 0.75_RK ) then
            nc1 = rnd(this%NComponents)
            do while (.not. this%Component(nc1)%permeable )
              nc1 = rnd(this%NComponents)
            end do
            nm1 = rnd(this%Component(nc1)%NPart)
            do while ( this%Component(nc1)%Pm0(nm1,1) .le. 0.25_RK .or. this%Component(nc1)%Pm0(nm1,1) .ge. 0.75_RK )
              nm1 = rnd(this%Component(nc1)%NPart)
            end do
            swap(:) = this%Component(i)%Pm0(j,:)
            this%Component(i)%Pm0(j,:) = this%Component(nc1)%Pm0(nm1,:)
            this%Component(nc1)%Pm0(nm1,:) = swap(:)
          end if
        end do
      end if
    end do

!    do i = 1, this%NComponents
!      do j = 1, this%Component(i)%NPart
!        write( IOBuffer, '(T10,I3," ",F9.3," ",F9.3," ",F9.3)' ) i, this%Component(i)%Pm0(j, 1), this%Component(i)%Pm0(j, 2), this%Component(i)%Pm0(j, 3)
!        call LogWrite
!      end do
!    end do
#endif

    do i = 1, this%NComponents
      this%Component(i)%P0 = this%Component(i)%P0 - 0.5_RK
    end do

    ! Save old positions
    do i = 1, this%NComponents
      this%Component(i)%Pm0old = this%Component(i)%Pm0
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
              pc%Qm0(j, k) = rnd( -1._RK, 1._RK )
            end do
            r = sum( pc%Qm0(j, :)**2 )
            if( r <= 1._RK ) exit
          end do
          pc%Qm0(j, :) = pc%Qm0(j, :) / sqrt( r )
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
      this%Component(i)%Pm0old = this%Component(i)%Pm0
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
      call RemoveNetMomentum( this%Component(i), this%Component(i)%Molecule%NUnit )
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
    real(RK)                  :: scale, Reference
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
      this%Temperature = 2._RK * this%EKin / (this%NDF-this%constrNDF) ! constrNDF due to Shake

      if(ConstantTemperature .or. NVTEquilibration) then
        Reference=this%RefTemperature
      else if (EnsembleType .eq. EnsembleTypeNVE ) then
        Reference= 2._RK * (this%RefHamiltonian*this%NPart - this%Epot) / real (this%NDF-this%constrNDF, RK)
      else if (EnsembleType .eq. EnsembleTypeNPH .and. .not. NVTEquilibration ) then
        Reference= 2._RK * (this%RefEnthalpy*this%NPart - this%Epot - this%RefPressure * this%Volume0) / real (this%NDF, RK) 
      end if

      ! Rescale velocities
      if( rescale ) then
      scale = sqrt( Reference / this%Temperature )
        do i = 1, this%NComponents
          pc => this%Component(i)
          np = pc%NPart
          pc%P1(:, :, :) = pc%P1(:, :, :) * scale
          if( pc%Molecule%isElongated ) then
            pc%W0(:, :, :) = pc%W0(:, :, :) * scale
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
    call MPI_Bcast( this%Temperature, 1, MPI_RK, NRootProc, Communicator, ierror )
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
    NPartOk = NPartOk .and. this%NPart > this%NPartLBound .and. this%NPart < this%NPartUBound

  end subroutine TEnsemble_CheckNPart


!==============================================================!
!  Subroutine TEnsemble_ResetEnsemble                          !
!==============================================================!

  subroutine TEnsemble_ResetEnsemble( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Calculate new initial density
    this%RefDensity = this%RefDensity * real( this%NPart, RK ) / real( this%NPartInitial, RK )
    this%NPart = this%NPartInitial

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

    ! Convert unit coordinates to atom positions
    call Mol2Unit( this )
    call Unit2Atom( this )

    if( SimulationType .eq. MolecularDynamics .and. .not. MCOverlapReduction ) then
      call InitMolecularDynamics( this, .false. )
    else
    ! Set all potential energy matrices
    call Energy( this, this%EPot )
    call UpdateEnergy( this )
    end if

  end subroutine TEnsemble_ResetEnsemble


!==============================================================!
!  Subroutine TEnsemble_RunMDStep                              !
!==============================================================!

  subroutine TEnsemble_RunMDStep( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer :: i, k, nc, np
    integer :: j     !debug
    real(RK) :: diffpressure

    ! Zero displacement
    if( Step == 1 ) then

      do i = 1, this%NComponents
        this%Component(i)%Disp(:, :) = 0._RK
      end do
      
      ! Initiate calculation of residence times
      if (this%ResidenceTime .and. .not. Equilibration) then 
        call ResidencePartners ( this )
      end if
      
    end if

    ! Run MD simulation step
    call Predict( this )
    call Unit2Atom( this )

    if( EnsembleType .eq. EnsembleTypeGE .and. (.not. NVTEquilibration) ) then

      if( Step == 1 ) call ZeroNAttempts( this )

      ! Update chemical potentials
      diffpressure = ( this%Pressure - this%RefPressure ) / this%Temperature
      do nc = 1, this%NComponents
        call UpdateChemPot( this%Component(nc), diffpressure )
      end do

      if ( mod(Step,int(0.01/timestep)) .eq. 0 ) then
        j = 0._RK
        k = rnd( 0._RK, 1._RK )
loop4:  do nc = 1, this%NComponents
          j = j + this%Component(nc)%Fraction
          if( k <= j ) exit loop4
        end do loop4
        call Insert( this, nc )

        j = 0._RK
        k = rnd( this%NPart )
loop5:  do nc = 1, this%NComponents
          j = j + this%Component(nc)%NPart
          if( k <= j ) exit loop5
        end do loop5
        np = 1 + j - k
        call Delete( this, nc, np )
      end if
      
      if ( mod(Step,100) .eq. 0 .and. RootProc ) call RemoveNetMomentum( this )

    else

      call ChemicalPotential( this )

    end if

    call Force( this )
    call Atom2Unit( this )
    call Correct( this )

#if CONSTR > 0
    call Constraints(this)
#endif

#if  TRANS == 1

!TRANSPORT_start
    if(.not. Equilibration .and. (mod((Step+this%NStepCorr-1),this%NStepCorr) .eq. 0)) then
      call CalCorrFun( this )  
    end if
!TRANSPORT_END

#endif

    ! Calculation of residence time
    if ( .not. Equilibration .and. this%ResidenceTime ) then
      call Residence ( this )
    end if

    call CalculateEKin( this, .true. )
    if( .not. Equilibration .and. this%RCutoffMax2 > this%BoxLength ) this%NRCutoffMax = this%NRCutoffMax + 1

#if HBOND > 0
    call HBonding(this)
#endif

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
    integer  :: i, NPart2, t, nu
    real(RK) :: rx, sx
    real(RK) :: diffpressure

    ! Zero number of MC attempts and successes
    if( Step == 1 ) call ZeroNAttempts( this )

    ! Outer loop
    do i = 1, this%NDF / 15 !3 ! Michael Sch.: decreased by 5 to sample smaller molecules better

      ! Choose particle randomly
      s = 0
      r = rnd( this%NDF )

loop1:do nc = 1, this%NComponents
        s = s + this%Component(nc)%NDF
        if( r <= s ) exit loop1
      end do loop1

      ndf = this%Component(nc)%Molecule%NDF
      np = 1 + (s - r) / ndf
      t = 1 + ((s - r) - (np-1) * ndf)

      ! Assign Unit
loop2:do nu = 1, this%Component(nc)%Molecule%NUnit
        if (t <= sum(this%Component(nc)%Molecule%Unit(1:nu)%NDF)) exit loop2
      end do loop2

      ! Move or Rotate Unit
      if ( this%Component(nc)%Molecule%Unit(nu)%isElongated ) then
        if( EnsembleType .eq. EnsembleTypeNVE .and. .not. NVTEquilibration) then
          ! Move or rotate for NVE ensemble
          if( mod( s - r, 2 ) .eq. 0 ) then
            call Move_NVE( this, nc, np, nu )
          else
            call Rotate_NVE( this, nc, np, nu )
          end if
            
        else if( EnsembleType .eq. EnsembleTypeNPH .and. .not. NVTEquilibration) then
          ! Move or rotate for NPH ensemble
          if( mod( s - r, 2 ) .eq. 0 ) then
            call Move_NPH( this, nc, np, nu )
          else
            call Rotate_NPH( this, nc, np, nu )
          end if
      
        else
          ! Move or rotate for constant temperature ensembles
          if( mod( s - r, 2 ) .eq. 0 ) then
            call Move( this, nc, np, nu )
          else
            call Rotate( this, nc, np, nu )
          end if
        endif
      else
        if( EnsembleType .eq. EnsembleTypeNVE .and. .not. NVTEquilibration) then
          call Move_NVE( this, nc, np, nu )
        else if( EnsembleType .eq. EnsembleTypeNPH .and. .not. NVTEquilibration) then
          call Move_NPH( this, nc, np, nu )
        else
          call Move( this, nc, np, nu )
        endif
      end if
    end do

    ! Special Moves for the simulation of molecules with internal degrees of freedom
    if ( UseIntDegFreed ) then
      NPart2 = 2*this%NPart
      ! Move or Rotate entire Molecule
      do i = 1, NPart2
        s = 0
        r = rnd( NPart2 )
loop3:  do nc = 1, this%NComponents 
          s = s + 2*this%Component(nc)%NPart
          if( r <= s ) exit loop3
        end do loop3

        ! Specify molecule
        np = int((s-r)/2+1)

        ! Move or Rotate
        if ( this%Component(nc)%Molecule%NUnit > 1 ) then
          if( EnsembleType .eq. EnsembleTypeNVE .and. .not. NVTEquilibration) then
            if( mod( s - r, 2 ) .eq. 0 ) then
              call Move_NVE( this, nc, np )
            else
              call Rotate_NVE( this, nc, np )
            end if
          else if( EnsembleType .eq. EnsembleTypeNPH .and. .not. NVTEquilibration) then
            if( mod( s - r, 2 ) .eq. 0 ) then
              call Move_NPH( this, nc, np )
            else
              call Rotate_NPH( this, nc, np )
            end if
          else
            if( mod( s - r, 2 ) .eq. 0 ) then
              call Move( this, nc, np )
            else
              call Rotate( this, nc, np )
            end if
          endif
        end if

      end do
    end if

    ! Calculate potential energy and virial
#if MPI_VER > 0
    ! in MC simulations we only communicate during common equilibration 
    if (Equilibration .and. CommonEqui) then

      ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
      call MPI_Allreduce( GetEnergy( this ), this%EPot, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
      call MPI_Allreduce( GetEnergyIntra( this ), this%EPotIntra, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
      if (printIDF) then
        call MPI_Allreduce( GetEnergyIntra_Bond( this ), this%EPotIntra_Bond, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
        call MPI_Allreduce( GetEnergyIntra_Angle( this ), this%EPotIntra_Angle, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
        call MPI_Allreduce( GetEnergyIntra_Dihedral( this ), this%EPotIntra_Dihedral, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
        this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
      endif
      this%EPotInter = this%EPot - this%EPotIntra
      call MPI_Allreduce( Getd2EpotdV2( this ), this%d2EpotdV2, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
        if ( this%OptPressure ) then
          call MPI_Allreduce( GetVirial( this ), this%Virial, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
          call MPI_Allreduce( GetVirialIntra( this ), this%VirialIntra, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
          this%VirialInter = this%Virial - this%VirialIntra
        endif

    else

      this%EPot = GetEnergy( this )
      this%EPotIntra   = GetEnergyIntra( this )
      if (printIDF) then
        this%EpotIntra_Bond = GetEnergyIntra_Bond( this )
        this%EpotIntra_Angle = GetEnergyIntra_Angle( this )
        this%EpotIntra_Dihedral = GetEnergyIntra_Dihedral( this )
        this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
      endif
      this%EPotInter   = this%EPot - this%EPotIntra
      this%d2EpotdV2 = Getd2EpotdV2( this )
      if ( this%OptPressure ) then
        this%Virial = GetVirial( this )
        this%VirialIntra = GetVirialIntra( this )
        this%VirialInter = this%Virial - this%VirialIntra
      endif

    endif  

#else

    this%EPot = GetEnergy( this )
    this%EPotIntra   = GetEnergyIntra( this )
    if (printIDF) then
      this%EpotIntra_Bond = GetEnergyIntra_Bond( this )
      this%EpotIntra_Angle = GetEnergyIntra_Angle( this )
      this%EpotIntra_Dihedral = GetEnergyIntra_Dihedral( this )
      this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
    endif
    this%EPotInter   = this%EPot - this%EPotIntra
    this%d2EpotdV2 = Getd2EpotdV2( this )
    if ( this%OptPressure ) then
      this%Virial = GetVirial( this )
      this%VirialIntra = GetVirialIntra( this )
      this%VirialInter = this%Virial - this%VirialIntra
    endif

#endif
    ! Resize simulation box
    if( ConstantPressure .and. .not. NVTEquilibration ) then
      call Resize( this )
      ! Check whether cutoff radius is too large
      if( this%RCutoffMax2 > this%BoxLength ) this%NRCutoffMax = this%NRCutoffMax + 1
    end if

    if ( this%OptPressure ) then
      ! Calculate pressure
      this%Pressure = (this%NUnitTotal * this%Temperature + this%Virial) / this%Volume0
    end if

    if( EnsembleType .eq. EnsembleTypeGE ) then

      if( .not. NVTEquilibration ) then
        ! Update chemical potentials
        diffpressure = ( this%Pressure - this%RefPressure ) / this%Temperature
        do nc = 1, this%NComponents
          call UpdateChemPot( this%Component(nc), diffpressure )
        end do

        ! Attempt inserts and deletes
        ! Prevent insertion of ions - look at Gibbs Ensemble - Particle Transfer

        do i = 1, 2
          sx = 0._RK
          rx = rnd( 0._RK, 1._RK )

loop4:    do nc = 1, this%NComponents
            sx = sx + this%Component(nc)%Fraction
            if( rx <= sx ) exit loop4
          end do loop4

          call Insert( this, nc )
          sx = 0._RK
          rx = rnd( this%NPart )

loop5:    do nc = 1, this%NComponents
            sx = sx + this%Component(nc)%NPart
            if( rx <= sx ) exit loop5
          end do loop5

          np = 1 + sx - rx
          call Delete( this, nc, np )

        end do
      end if

    elseif( EnsembleType .eq. EnsembleTypeHA ) then

      if( .not. NVTEquilibration ) then
        ! Attempt inserts and deletes on phase changing component (first one)
        do i = 1, 2
          if( rnd( 0._RK, 1._RK ) < this%Component(1)%Fraction ) call Insert( this, 1 )
          np = rnd( this%NPart )
          if( np <= this%Component(1)%NPart ) call Delete( this, 1, np )
        end do

      end if

    else

      ! Calculate chemical potential
      call ChemicalPotential( this )

    end if

    ! Update MC displacements
    if( Equilibration .and. mod( Step, DispUpdateFrequency ) == 0 ) then
      call UpdateDisplacements( this )
    end if

#if HBOND > 0
    call HBonding(this)
#endif

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
    integer                     :: nu
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
      this%Component(i)%P0(:, 1, :) = this%Component(i)%P0(:, 1, :) + r / this%BoxLength ! (i+1,1:3) and (i,2:3) set to 0.0 in ConstructSVC
      call Unit2Atom(this%Component(i), this%Component(i)%NPart, this%Component(i)%Molecule%NUnit)
    end do

    ! Loop over components
    do i = 1, this%NComponents, 2
      do j = i + 1, this%NComponents, 2
        pi => this%Interaction(i, j)
        n = pi%NPart2*pi%NUnit2
        pi%MayerFFunction(Step) = 0._RK

        ! Loop over units
        do np = 1, this%Component(i)%NPart
          do nu = 1, this%Component(i)%Molecule%NUnit
            call Energy( pi, np, nu, this%BoxLength )
            if ( pi%SameComponent .and. UseIntDegFreed ) then
              call IntraEnergy( pi, np, nu, this%BoxLength )
            end if
          end do

          ! Sum Mayer f-function
          pi%MayerFFunction(Step) = pi%MayerFFunction(Step) + sum( exp( betaneg * pi%EPot1(1:n) ) - 1._RK )
          pi%MayerFFunction1(Step) = pi%MayerFFunction1(Step) + sum( exp( betaneg1 * pi%EPot1(1:n) ) - 1._RK )
          pi%MayerFFunction2(Step) = pi%MayerFFunction2(Step) + sum( exp( betaneg2 * pi%EPot1(1:n) ) - 1._RK )

        end do

        ! Average Mayer f-function
        pi%MayerFFunction(Step) = pi%MayerFFunction(Step) / real( n * this%Component(i)%NPart, RK )
        pi%MayerFFunction1(Step) = pi%MayerFFunction1(Step) / real( n * this%Component(i)%NPart, RK )
        pi%MayerFFunction2(Step) = pi%MayerFFunction2(Step) / real( n * this%Component(i)%NPart, RK )

        ! Integrate Mayer f-function
        pi%IntFFunction = Bij0 + Piminus2 * simpson( pi%MayerFFunction * rsquared, rdist, Step)
        pi%IntFFunction1 = Bij0 + Piminus2 * simpson( pi%MayerFFunction1 * rsquared, rdist, Step)
        pi%IntFFunction2 = Bij0 + Piminus2 * simpson( pi%MayerFFunction2 * rsquared, rdist, Step)

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
    write( IOBuffer, '(I9)' ) Step
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
      real(RK), intent(in) :: values(n), stepsize
      integer, intent(in)  :: n

      ! Declare result
      real(RK) :: integral(n)

      ! Declare local variables
      integer :: i

      ! Return if no values to integrate
      if( n < 1 ) return

      ! Initialize result
      integral(:) = 0._RK

      ! Calculate integral via Simpson's rule
      do i = 3, n, 2
        integral(i) = integral(i-2) + values(i) + 4._RK * values(i-1) + values(i-2)
        integral(i-1) = .5 * (integral(i) + integral(i-2))
      end do

      if( mod(n, 2) == 0 .and. n > 2 ) integral(n) = integral(n-1) + .5_RK * values(n) + 2._RK * values(n-1) + .5_RK * values(n-2)
      integral(:) = integral(:) * Third * stepsize

    end function

  end subroutine TEnsemble_RunSVCStep


!==============================================================!
!  Subroutine TEnsemble_Unit2Atom                               !
!==============================================================!

  subroutine TEnsemble_Unit2Atom( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer      :: i!, i0, i1

    ! Call Unit2Atom for each component
    do i = 1, this%NComponents
! #if MPI_VER > 0
!       i0 = this%Component(i)%NPart0
!       i1 = this%Component(i)%NPart2
! #else
!       i0 = 1
!       i1 = this%Component(i)%NPart1
! #endif
!       call Mol2Atom( this%Component(i), i0, i1-i0+1 )
      call Unit2Atom( this%Component(i), this%Component(i)%NPart, this%Component(i)%Molecule%NUnit )
    end do

  end subroutine TEnsemble_Unit2Atom


!==============================================================!
!  Subroutine TEnsemble_Atom2Unit                              !
!==============================================================!

  subroutine TEnsemble_Atom2Unit( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer   :: i!, i0, i1
#if OSMOP > 0
    integer   :: m
    real(RK)  :: TotalDenProfile(NBinsDen)

    ! Osmotic pressure
    this%OsmoticPressure = 0._RK
    TotalDenProfile(:) = 0._RK
#endif

    ! Call Atom2Unit for each component
    do i = 1, this%NComponents
! #if MPI_VER > 0
!       i0 = this%Component(i)%NPart0
!       i1 = this%Component(i)%NPart2
! #else
!       i0 = 1
!       i1 = this%Component(i)%NPart1
! #endif
! #if  TRANS == 1
!       if(.not. Equilibration .and. (mod((Step+this%NStepCorr-1),this%NStepCorr) .eq. 0)) then
!          call Atom2Mol_Trans( this%Component(i), i0, i1-i0+1 )
!       else
!          call Atom2Mol( this%Component(i), i0, i1-i0+1 )
!       end if
! #else
!       call Atom2Mol( this%Component(i), i0, i1-i0+1 )
! #endif
#if  TRANS == 1
      if(.not. Equilibration .and. (mod((Step+this%NStepCorr-1),this%NStepCorr) .eq. 0)) then
         if (this%Component(i)%Molecule%NUnit .ne. 1)  call Error( "!!!!!!Transportproperties only implemented for rigid molecules!!!!!!!" )
         call Atom2Unit_Trans( this%Component(i), this%Component(i)%NPart, this%Component(i)%Molecule%NUnit )
      else
         call Atom2Unit( this%Component(i), this%Component(i)%NPart, this%Component(i)%Molecule%NUnit )
      end if
#else
      call Atom2Unit( this%Component(i), this%Component(i)%NPart, this%Component(i)%Molecule%NUnit )
#endif

#if OSMOP == 0
    end do
#else
      ! Dichteprofil für die Bestimmung des chemischen Potentials
      call DensityProfile( this%Component(i) )

      !if(this%Component(i)%Molecule%Charge .ne. 0)then 
      if ( .not. this%Component(i)%permeable ) then   
        this%OsmoticPressure = this%OsmoticPressure + sum(this%Component(i)%FOsmoticPressure(:))
      end if
      do m = 1, NBinsDen
        TotalDenProfile(m) = TotalDenProfile(m) + this%Component(i)%SumDenProfile(m)%Average
      end do

    end do

#if OSMOP == 2
    !Correct virial profile
    do m = 1, NBinsDen
       this%VirialProfile(m) = this%VirialProfile(m) + (TotalDenProfile(m) * this%VirialCorrLJ * NProcs)/NBinsDen
    end do

    if (LongRange .eq. Ewald) then
      this%VirialProfile(:) = this%VirialProfile(:) + this%EVirial/NBinsDen
#if SPME > 0
    else if (LongRange .eq. PME) then
      this%VirialProfile(:) = this%VirialProfile(:) + this%EVirial/NBinsDen
#endif
    end if
 
    !Calculation of the pressure profile (real and ideal part)
    do m = 1, NBinsDen
      this%PressureProfile(m) = this%VirialProfile(m)*NBinsDen/this%Volume0 + TotalDenProfile(m) * this%Temperature
    end do
    this%VirialProfile(:) = 0._RK
#endif

    this%OsmoticPressure = 0.5_RK * this%OsmoticPressure / (this%BoxLength*this%BoxLength)
#endif

  end subroutine TEnsemble_Atom2Unit


!==============================================================!
!  Subroutine TEnsemble_Mol2Unit                               !
!==============================================================!

  subroutine TEnsemble_Mol2Unit( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer      :: i

    ! Call Mol2Unit for each component
    do i = 1, this%NComponents
      call Mol2Unit( this%Component(i), this%Component(i)%NPart, &
&                      this%Component(i)%Molecule%NUnit )
    end do

  end subroutine TEnsemble_Mol2Unit


!==============================================================!
!  Subroutine TEnsemble_Unit2Mol                               !
!==============================================================!

  subroutine TEnsemble_Unit2Mol( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer      :: i

    ! Call Unit2Mol for each component
    do i = 1, this%NComponents
      call Mol2Unit( this%Component(i), this%Component(i)%NPart, &
&                      this%Component(i)%Molecule%NUnit )
    end do

  end subroutine TEnsemble_Unit2Mol


!==============================================================!
!  Subroutine TEnsemble_ResizeMol                              !
!==============================================================!

  subroutine TEnsemble_ResizeMol( this, DelBoxFrac )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this
    real(RK), intent(in) :: DelBoxFrac

    ! Declare local variables
    integer      :: i

    ! Call ResizeMol for each component
    do i = 1, this%NComponents
      call ResizeMol( this%Component(i), DelBoxFrac )
    end do

  end subroutine TEnsemble_ResizeMol


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

    if( ConstantPressure .and. .not. NVTEquilibration ) then
      call PredictVol( this )
    end if

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

    if( ConstantPressure .and. .not. NVTEquilibration ) then
      call CorrectVol(this)
    end if

  end subroutine TEnsemble_Correct


!==============================================================!
!  Subroutine TEnsemble_PredictGear                            !
!==============================================================!

  subroutine TEnsemble_PredictGear( this )

    implicit none

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

  end subroutine TEnsemble_PredictGear



!==============================================================!
!  Subroutine TEnsemble_CorrectGear                            !
!==============================================================!

  subroutine TEnsemble_CorrectGear( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer  :: i
    real(RK) :: dLogVolumeThird

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

  end subroutine TEnsemble_CorrectGear



!==============================================================!
!  Subroutine TEnsemble_PredictLeapFrog                        !
!==============================================================!

  subroutine TEnsemble_PredictLeapFrog( this )

    implicit none

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
    type(TComponent), pointer :: pc
    real(RK)            :: EPot, Virial, d2EpotdV2
    integer             :: i, j
    real(RK)            :: EPotIntra, VirialIntra
    real(RK)            :: EPotIntra_Bond, EPotIntra_Angle, EPotIntra_Dihedral
    real(RK)            :: EPotIntra_Nonbonded
    real(RK)            :: EPotInter, VirialInter
#ifdef ABL
    integer             :: k,l
    integer             :: numbi, numbj, numb
#endif
! #if MPI_VER > 0
!     integer             :: i0, i1, np, np1, n
! #endif
#if OSMOP == 2
#if MPI_VER > 0
    real(RK),allocatable :: VirialProfile(:)

    allocate( VirialProfile(NBinsDen) )
    VirialProfile(:) = 0._RK
#endif
#endif

! Zero forces
    do i = 1, this%NComponents
      pc => this%Component(i)
      do j = 1, this%Component(i)%Molecule%NLJ126
        pc%Molecule%SiteLJ126(j)%FX(1:pc%NPart) = 0._RK
        pc%Molecule%SiteLJ126(j)%FY(1:pc%NPart) = 0._RK
        pc%Molecule%SiteLJ126(j)%FZ(1:pc%NPart) = 0._RK
#if  TRANS == 1
        !TRANSPORT_start
        if(mod((Step+this%NStepCorr-1),this%NStepCorr) .eq. 0) then
          pc%Molecule%SiteLJ126(j)%vsLJx(1:pc%NPart) = 0._RK
          pc%Molecule%SiteLJ126(j)%vsLJy(1:pc%NPart) = 0._RK
          pc%Molecule%SiteLJ126(j)%vsLJz(1:pc%NPart) = 0._RK
          pc%Molecule%SiteLJ126(j)%vbLJx(1:pc%NPart) = 0._RK
          pc%Molecule%SiteLJ126(j)%vbLJy(1:pc%NPart) = 0._RK
          pc%Molecule%SiteLJ126(j)%vbLJz(1:pc%NPart) = 0._RK
  !        if ( this%Conductivity ) then
            pc%Molecule%SiteLJ126(j)%vsuLJx(1:pc%NPart)= 0._RK
            pc%Molecule%SiteLJ126(j)%vsuLJy(1:pc%NPart)= 0._RK
            pc%Molecule%SiteLJ126(j)%vsuLJz(1:pc%NPart)= 0._RK
            pc%Molecule%SiteLJ126(j)%cLJx(1:pc%NPart)  = 0._RK
            pc%Molecule%SiteLJ126(j)%cLJy(1:pc%NPart)  = 0._RK
            pc%Molecule%SiteLJ126(j)%cLJz(1:pc%NPart)  = 0._RK
            pc%Molecule%SiteLJ126(j)%tuLJx(1:pc%NPart) = 0._RK
            pc%Molecule%SiteLJ126(j)%tuLJy(1:pc%NPart) = 0._RK
            pc%Molecule%SiteLJ126(j)%tuLJz(1:pc%NPart) = 0._RK
            pc%Molecule%SiteLJ126(j)%tlLJx(1:pc%NPart) = 0._RK
            pc%Molecule%SiteLJ126(j)%tlLJy(1:pc%NPart) = 0._RK
            pc%Molecule%SiteLJ126(j)%tlLJz(1:pc%NPart) = 0._RK
            pc%Molecule%SiteLJ126(j)%tdLJx(1:pc%NPart) = 0._RK
            pc%Molecule%SiteLJ126(j)%tdLJy(1:pc%NPart) = 0._RK
            pc%Molecule%SiteLJ126(j)%tdLJz(1:pc%NPart) = 0._RK
   !       end if
        end if
        !TRANSPORT_END
#endif
      end do
      do j = 1, this%Component(i)%Molecule%NCharge
        pc%Molecule%SiteCharge(j)%FX(1:pc%NPart) = 0._RK
        pc%Molecule%SiteCharge(j)%FY(1:pc%NPart) = 0._RK
        pc%Molecule%SiteCharge(j)%FZ(1:pc%NPart) = 0._RK
#if  TRANS == 1
        !TRANSPORT_start
        if(mod((Step+this%NStepCorr-1),this%NStepCorr) .eq. 0) then
          pc%Molecule%SiteCharge(j)%vsCx(1:pc%NPart) = 0._RK
          pc%Molecule%SiteCharge(j)%vsCy(1:pc%NPart) = 0._RK
          pc%Molecule%SiteCharge(j)%vsCz(1:pc%NPart) = 0._RK
          pc%Molecule%SiteCharge(j)%vbCx(1:pc%NPart) = 0._RK
          pc%Molecule%SiteCharge(j)%vbCy(1:pc%NPart) = 0._RK
          pc%Molecule%SiteCharge(j)%vbCz(1:pc%NPart) = 0._RK
    !      if ( this%Conductivity ) then
            pc%Molecule%SiteCharge(j)%vsuCx(1:pc%NPart)= 0._RK
            pc%Molecule%SiteCharge(j)%vsuCy(1:pc%NPart)= 0._RK
            pc%Molecule%SiteCharge(j)%vsuCz(1:pc%NPart)= 0._RK
            pc%Molecule%SiteCharge(j)%cCx(1:pc%NPart)  = 0._RK
            pc%Molecule%SiteCharge(j)%cCy(1:pc%NPart)  = 0._RK
            pc%Molecule%SiteCharge(j)%cCz(1:pc%NPart)  = 0._RK
            pc%Molecule%SiteCharge(j)%tuCx(1:pc%NPart) = 0._RK
            pc%Molecule%SiteCharge(j)%tuCy(1:pc%NPart) = 0._RK
            pc%Molecule%SiteCharge(j)%tuCz(1:pc%NPart) = 0._RK
            pc%Molecule%SiteCharge(j)%tlCx(1:pc%NPart) = 0._RK
            pc%Molecule%SiteCharge(j)%tlCy(1:pc%NPart) = 0._RK
            pc%Molecule%SiteCharge(j)%tlCz(1:pc%NPart) = 0._RK
            pc%Molecule%SiteCharge(j)%tdCx(1:pc%NPart) = 0._RK
            pc%Molecule%SiteCharge(j)%tdCy(1:pc%NPart) = 0._RK
            pc%Molecule%SiteCharge(j)%tdCz(1:pc%NPart) = 0._RK
     !     end if
        end if
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
        if(mod((Step+this%NStepCorr-1),this%NStepCorr) .eq. 0) then
          pc%Molecule%SiteDipole(j)%vsDx(1:pc%NPart) = 0._RK
          pc%Molecule%SiteDipole(j)%vsDy(1:pc%NPart) = 0._RK
          pc%Molecule%SiteDipole(j)%vsDz(1:pc%NPart) = 0._RK
          pc%Molecule%SiteDipole(j)%vbDx(1:pc%NPart) = 0._RK
          pc%Molecule%SiteDipole(j)%vbDy(1:pc%NPart) = 0._RK
          pc%Molecule%SiteDipole(j)%vbDz(1:pc%NPart) = 0._RK
      !    if ( this%Conductivity ) then
            pc%Molecule%SiteDipole(j)%vsuDx(1:pc%NPart)= 0._RK
            pc%Molecule%SiteDipole(j)%vsuDy(1:pc%NPart)= 0._RK
            pc%Molecule%SiteDipole(j)%vsuDz(1:pc%NPart)= 0._RK
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
       !   end if
        end if
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
        if(mod((Step+this%NStepCorr-1),this%NStepCorr) .eq. 0) then
          pc%Molecule%SiteQuadrupole(j)%vsQx(1:pc%NPart) = 0._RK
          pc%Molecule%SiteQuadrupole(j)%vsQy(1:pc%NPart) = 0._RK
          pc%Molecule%SiteQuadrupole(j)%vsQz(1:pc%NPart) = 0._RK
          pc%Molecule%SiteQuadrupole(j)%vbQx(1:pc%NPart) = 0._RK
          pc%Molecule%SiteQuadrupole(j)%vbQy(1:pc%NPart) = 0._RK
          pc%Molecule%SiteQuadrupole(j)%vbQz(1:pc%NPart) = 0._RK
        !  if ( this%Conductivity ) then
            pc%Molecule%SiteQuadrupole(j)%vsuQx(1:pc%NPart)= 0._RK
            pc%Molecule%SiteQuadrupole(j)%vsuQy(1:pc%NPart)= 0._RK
            pc%Molecule%SiteQuadrupole(j)%vsuQz(1:pc%NPart)= 0._RK
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
         ! end if
        end if
        !TRANSPORT_END
#endif
      end do

      do j = 1, this%Component(i)%Molecule%NBond
        pc%Molecule%IdfBond(j)%FX1(1:pc%NPart) = 0._RK
        pc%Molecule%IdfBond(j)%FY1(1:pc%NPart) = 0._RK
        pc%Molecule%IdfBond(j)%FZ1(1:pc%NPart) = 0._RK
        pc%Molecule%IdfBond(j)%FX2(1:pc%NPart) = 0._RK
        pc%Molecule%IdfBond(j)%FY2(1:pc%NPart) = 0._RK
        pc%Molecule%IdfBond(j)%FZ2(1:pc%NPart) = 0._RK
      end do
      do j = 1, this%Component(i)%Molecule%NAngle
        pc%Molecule%IdfAngle(j)%FX1(1:pc%NPart) = 0._RK
        pc%Molecule%IdfAngle(j)%FY1(1:pc%NPart) = 0._RK
        pc%Molecule%IdfAngle(j)%FZ1(1:pc%NPart) = 0._RK
        pc%Molecule%IdfAngle(j)%FX2(1:pc%NPart) = 0._RK
        pc%Molecule%IdfAngle(j)%FY2(1:pc%NPart) = 0._RK
        pc%Molecule%IdfAngle(j)%FZ2(1:pc%NPart) = 0._RK
        pc%Molecule%IdfAngle(j)%FX3(1:pc%NPart) = 0._RK
        pc%Molecule%IdfAngle(j)%FY3(1:pc%NPart) = 0._RK
        pc%Molecule%IdfAngle(j)%FZ3(1:pc%NPart) = 0._RK
      end do
      do j = 1, this%Component(i)%Molecule%NDihedral
        pc%Molecule%IdfDihedral(j)%FX1(1:pc%NPart) = 0._RK
        pc%Molecule%IdfDihedral(j)%FY1(1:pc%NPart) = 0._RK
        pc%Molecule%IdfDihedral(j)%FZ1(1:pc%NPart) = 0._RK
        pc%Molecule%IdfDihedral(j)%FX2(1:pc%NPart) = 0._RK
        pc%Molecule%IdfDihedral(j)%FY2(1:pc%NPart) = 0._RK
        pc%Molecule%IdfDihedral(j)%FZ2(1:pc%NPart) = 0._RK
        pc%Molecule%IdfDihedral(j)%FX3(1:pc%NPart) = 0._RK
        pc%Molecule%IdfDihedral(j)%FY3(1:pc%NPart) = 0._RK
        pc%Molecule%IdfDihedral(j)%FZ3(1:pc%NPart) = 0._RK
        pc%Molecule%IdfDihedral(j)%FX4(1:pc%NPart) = 0._RK
        pc%Molecule%IdfDihedral(j)%FY4(1:pc%NPart) = 0._RK
        pc%Molecule%IdfDihedral(j)%FZ4(1:pc%NPart) = 0._RK
      end do

      if( pc%Molecule%isElongated ) then
        pc%tRFX(:, :) = 0._RK
        pc%tRFY(:, :) = 0._RK
        pc%tRFZ(:, :) = 0._RK
      end if
#if  TRANS == 1
      !TRANSPORT_start
      if(mod((Step+this%NStepCorr-1),this%NStepCorr) .eq. 0) then
        do j = 1, this%Component(i)%NPart
          this%Component(i)%FS(j, 1)    = 0._RK
          this%Component(i)%FS(j, 2)    = 0._RK
          this%Component(i)%FS(j, 3)    = 0._RK
          this%Component(i)%FB(j, 1)    = 0._RK
          this%Component(i)%FB(j, 2)    = 0._RK
          this%Component(i)%FB(j, 3)    = 0._RK
    !      if ( this%Conductivity ) then
            this%Component(i)%FTC(j, 1)   = 0._RK
            this%Component(i)%FTC(j, 2)   = 0._RK
            this%Component(i)%FTC(j, 3)   = 0._RK
            this%Component(i)%FRC(j, 1)   = 0._RK
            this%Component(i)%FRC(j, 2)   = 0._RK
            this%Component(i)%FRC(j, 3)   = 0._RK
     !     end if
        end do
      end if
      !TRANSPORT_END
#endif
    end do

    ! Zero potential
    EPot = this%Density * this%EPotCorrLJ + this%EPotCorrRF
    EPotInter = this%Density * this%EPotCorrLJ + this%EPotCorrRF
    EPotIntra = 0._RK
    EPotIntra_Bond = 0._RK
    EPotIntra_Angle = 0._RK
    EPotIntra_Dihedral = 0._RK
    EPotIntra_Nonbonded = 0._RK

    ! Zero virial
    Virial = this%Density * this%VirialCorrLJ + this%VirialCorrRF*this%Volume0
    VirialInter = Virial
    VirialIntra = 0._RK

    ! Zero d2Epot/dV2
    d2EpotdV2 = this%Density * this%d2EpotdV2CorrLJ 

!     ! Calculate interactions partners within cutoff sphere
!     if( CutoffMode .eq. CenterofMass ) then
!       do i = 1, this%NComponents
!         do j = i, this%NComponents
!           call CalcCutoffPartners( this%Interaction( i, j ) )
!         end do
!       end do
!     end if
! 
! #if MPI_VER > 0
!     do i = 1, this%NComponents
!       this%Component(i)%NAdd(:) = .false.
!     end do
!     do i = 1, this%NComponents
!       i0 = this%Component(i)%NPart0
!       i1 = this%Component(i)%NPart2
!       do j = i, this%NComponents
!         ! Calculate SitePositions by demand
!         do np1 = i0, i1
!           do n = 1, this%Interaction( i, j )%NInCutoff(np1)
!             np = this%Interaction( i, j )%CutoffPartner(n, np1)
!             if ( np < this%Component(j)%NPart0 .or. np > this%Component(j)%NPart2 ) then
!               this%Component(j)%NAdd(np) = .true.
!               call Mol2Atom1( this%Component(j), np)
!             end if
!           end do
!         end do
!       end do
!     end do
! #endif

    ! Loop over components
    do i = 1, this%NComponents
      do j = i, this%NComponents
#if TRANS == 1
#ifndef ABL
        if(.not. Equilibration .and. (mod((Step+this%NStepCorr-1),this%NStepCorr) .eq. 0)) then
           call Force_Trans( this%Interaction( i, j ), EPot, Virial, EPotIntra, EPotIntra_Bond,  &
&                           EPotIntra_Angle, EPotIntra_Dihedral, EPotIntra_Nonbonded, EPotInter, &
&                           VirialIntra, VirialInter, d2EpotdV2, this%BoxLength )
        else
          call Force( this%Interaction( i, j ), EPot, Virial, EPotIntra, EPotIntra_Bond,  &
&                     EPotIntra_Angle, EPotIntra_Dihedral, EPotIntra_Nonbonded, EPotInter, &
&                     VirialIntra, VirialInter, d2EpotdV2, this%BoxLength )
        endif
#else
        this%Interaction(i,j)%AblPS => this%AblPS
        this%Interaction(i,j)%AblPE => this%AblPE
        if(.not. Equilibration .and. (mod((Step+this%NStepCorr-1),this%NStepCorr) .eq. 0)) then
           call Force_Trans( this%Interaction( i, j ), EPot, Virial, EPotIntra, EPotIntra_Bond,  &
&                           EPotIntra_Angle, EPotIntra_Dihedral, EPotIntra_Nonbonded, EPotInter, &
&                           VirialIntra, VirialInter, d2EpotdV2, this%BoxLength, i, j )
        else
          call Force( this%Interaction( i, j ), EPot, Virial, EPotIntra, EPotIntra_Bond,  &
&                     EPotIntra_Angle, EPotIntra_Dihedral, EPotIntra_Nonbonded, EPotInter, &
&                     VirialIntra, VirialInter, d2EpotdV2, this%BoxLength, i, j )
        endif
#endif
#else
#ifndef ABL
        call Force( this%Interaction( i, j ), EPot, Virial, EPotIntra, EPotIntra_Bond,  &
&                   EPotIntra_Angle, EPotIntra_Dihedral, EPotIntra_Nonbonded, EPotInter, &
&                   VirialIntra, VirialInter, d2EpotdV2, this%BoxLength )
#else
        this%Interaction(i,j)%AblPS => this%AblPS
        this%Interaction(i,j)%AblPE => this%AblPE
        call Force( this%Interaction( i, j ), EPot, Virial, EPotIntra, EPotIntra_Bond,  &
&                   EPotIntra_Angle, EPotIntra_Dihedral, EPotIntra_Nonbonded, EPotInter, &
&                   VirialIntra, VirialInter, d2EpotdV2, this%BoxLength, i, j )
#endif
#endif

      end do
    end do

    if (LongRange .eq. Ewald) then
      call EwaldFourierTerm (this)
#if SPME > 0
    else if (LongRange .eq. PME) then
      call PMEFourierTerm (this)
#endif
    end if


    ! Collect sums from all processes
#if MPI_VER > 0
    ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
    call MPI_Reduce( EPot, this%EPot, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    call MPI_Reduce( Virial, this%Virial, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    call MPI_Reduce( EPotInter, this%EPotInter, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    call MPI_Reduce( EPotIntra, this%EPotIntra, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    if (printIDF) then
      call MPI_Reduce( EPotIntra_Bond, this%EPotIntra_Bond, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( EPotIntra_Angle, this%EPotIntra_Angle, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( EPotIntra_Dihedral, this%EPotIntra_Dihedral, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( EPotIntra_Nonbonded, this%EPotIntra_Nonbonded, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    endif
    call MPI_Reduce( VirialInter, this%VirialInter, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    call MPI_Reduce( VirialIntra, this%VirialIntra, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    call MPI_Reduce( d2EpotdV2, this%d2EpotdV2, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
#if OSMOP == 2
    call MPI_Reduce( this%VirialProfile(:), VirialProfile(:), NBinsDen, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
    if (RootProc) this%VirialProfile(:) = VirialProfile(:)
#endif
#else
    this%EPot = EPot
    this%Virial = Virial
    this%EPotInter = EPotInter
    this%EPotIntra = EPotIntra
    if (printIDF) then
      this%EPotIntra_Bond = EPotIntra_Bond
      this%EPotIntra_Angle = EPotIntra_Angle
      this%EPotIntra_Dihedral = EPotIntra_Dihedral
      this%EPotIntra_Nonbonded = EPotIntra_Nonbonded
    endif
    this%VirialInter = VirialInter
    this%VirialIntra = VirialIntra
    this%d2EpotdV2 = d2EpotdV2
#endif

    if (RootProc) then
      if (LongRange .eq. Ewald) then
        this%EPot   = this%EPot   + this%UFourier + this%USelbstTerm + this%UIntra
        this%Virial = this%Virial + this%EVirial
#if SPME > 0
      else if (LongRange .eq. PME) then
        this%EPot   = this%EPot   + this%UFourier + this%USelbstTerm + this%UIntra
        this%Virial = this%Virial + this%EVirial
#endif
      end if
    end if


    this%Pressure = (this%NUnitTotal * this%Temperature + this%Virial) / this%Volume0

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
    real(RK)                  :: ChemPot, F(3,this%NUnitMax), rm(3), ExpMinusBetaEnLaMin, factor
    real(RK)                  :: HW_H_local, HW_V_local, HW_counter_local, HW_denom_local
    integer                   :: i, j, t
    real(RK)                  :: EPotTest(this%NTestMax)
    real(RK)                  :: EBin, E, EIntra, EBond, EAngle, EDihedral
    integer                   :: ndf, ndfmove, ndfbiased, ndffluct, ndfchange, ndfcp
    integer                   :: r, s, nc, np, ncf, npf
    integer                   :: ewald_h, ratio, sndf, nuh
    integer                   :: nu, nuh2, j0, j1, selected
    type(TComponent), pointer :: pc
    integer                   :: nstate( 0:this%NFluctMax )
#if MPI_VER > 0
    integer                   :: tempComm
    integer                   :: tempVec(0:this%NFluctMax)
    real(RK)                  :: EPot_h
    integer                   :: tempVal, tempVal2
    integer                   :: tempVec1(this%NFluctMax), tempVec2(this%NFluctMax)
    integer                   :: tempVec3(this%NFluctMax), tempVec4(this%NFluctMax)
#endif
#if OSMOP == 2
    real(RK)                  :: ChemPotProfile(NBinsDen)
    integer                   :: NTestBinDen, m
#if MPI_VER > 0
    integer                   :: NTestBinDenAll
#endif
#endif

    ! No calculation of chemical potential in equilibration
    if( Equilibration ) then
      if (Step == 1) then
        do i = 1, this%NRealComponents
          this%Component(i)%CalcChemPot = .false.
          this%Component(i)%ChemPot = 0._RK
        end do
      end if
      return
    end if

    if (this%NTestMax == 0) return

    if (Step == 1) then
      do i = 1, this%NComponents
        pc => this%Component(i)
        pc%NTest1 = ProcRange( pc%NTest, pc%NTest0, pc%NTest2 )
        select case( pc%ChemPotMethod )
        case (ChemPotMethodWidom)
          pc%CalcChemPot = .true.
        case (ChemPotMethodGradIns)
          pc%ProbW0      = 0._RK
          pc%ProbW1      = 0._RK
          pc%ProbW0V     = 0._RK
          pc%ProbW1Rho   = 0._RK
          pc%ChemPot1    = 0._RK
          pc%ChemPot2    = 0._RK
          pc%NStateWF(:) = 0
          pc%CalcChemPot= .true.
        case (ChemPotMethodThermoInt)
          pc%BinsVisit(:)=0
          pc%BinsEn(:)=0.0_RK
          pc%BinsdEndLa(:)=0.0_RK
          pc%BinsIntdEndLa(:)=0.0_RK
          pc%BinsdEndLaV(:)=0.0_RK
          pc%BinsdEndLaH(:)=0.0_RK
          pc%BinsIntVW(:)=0.0_RK
          pc%BinsIntHW(:)=0.0_RK
          pc%CalcChemPot= .true.
        case default
          pc%ChemPot = 0._RK
          pc%CalcChemPot = .false.
        end select
      end do
    end if

    ! Throw test particles
    if( mod( Step, BlockSize ) == 1 ) then
      do i = 1, this%NComponents
        if (this%Component(i)%NTest > 0) then
          pc => this%Component(i)
          pc%EPotTestIntra(:) = 0._RK
#if MPI_VER > 0
          if ( (SimulationType .ne. MonteCarlo) .or. (Equilibration .and. CommonEqui) ) then
            j0 = pc%NTest0
            j1 = pc%NTest2
          else
            j0 = 1
            j1 = pc%NTest
          end if
#else
          j0 = 1
          j1 = pc%NTest
#endif
          if (pc%NPart == 0) then ! for Henry NPart==0
            do j = j0, j1
              do t = 1, 3
                rm(t) = rnd( -.5_RK, .5_RK )
              end do
              do r = 1, pc%Molecule%NUnit
                pc%P0Test(j,:,r) = pc%Molecule%Unit(r)%P0(:) + rm(:)
              end do
              if (pc%Molecule%isElongated) then
                do r = 1, pc%Molecule%NUnit
                  pc%Q0Test(j,:,r) = pc%Molecule%Unit(r)%Q0(:)
                end do
                do t = 1, 3
                  rm(t) = rnd( -1._RK, 1._RK )
                end do
                call RotateTest( pc, j, rm)
              end if 
            end do
            call Unit2AtomTest( pc, pc%Ntest, pc%Molecule%NUnit )

          else ! not Henry or Fraction > 0
            call Unit2Mol(pc) ! needed? Michael Sch.
            do j = j0, j1
              do t = 1, 3
                rm(t) = rnd( -.5_RK, .5_RK )
              end do
              selected = rnd( pc%NPart )
              do r = 1, pc%Molecule%NUnit
                pc%P0Test(j,1:3,r) = pc%P0(selected,1:3,r) + rm(1:3)
              end do
              do r = 1, pc%Molecule%NUnit
                pc%P0Test(j,1:3,r) = pc%P0Test(j,1:3,r) - pc%Pm0(selected,1:3)
              end do

              if (pc%Molecule%isElongated) then
                pc%Q0Test(j,:,:) = pc%Q0(selected,:,:)
                do t = 1, 3
                  rm(t) = rnd( -1._RK, 1._RK )
                end do
                call RotateTest( pc, j, rm)
              end if
              if (SimulationType .eq. MolecularDynamics) then
                E = 0._RK; EIntra = 0._RK; EBond = 0._RK; EAngle = 0._RK; EDihedral = 0._RK; F(:,:) = 0._RK
                t = this%Component(i)%Molecule%NUnit
                call MDEnergy( this%Interaction(i,i), selected, t, F(:,1:t), E, EIntra, EBond, EAngle, EDihedral, this%BoxLength, .true. )
                pc%EPotTestIntra(j) = pc%EPotTestIntra(j) + EIntra
              else
                pc%EPotTestIntra(j) = pc%EPotTestIntra(j) + GetEnergyIntra(this, i, selected)
              end if
            end do
            call Unit2AtomTest( pc, pc%Ntest, pc%Molecule%NUnit )
          end if
        end if
      end do

    end if

#if MPI_VER > 0
    tempComm = Communicator
#endif

    t = this%NRealComponents+1  ! pseudo component identifier for ThermoInt (ThermoInt does not function together with GradIns)

    ! Outer loop over components
componentLoop:       do i = 1, this%NRealComponents
      pc => this%Component(i)
      if( Equilibration .and. pc%WFMethod .ne. WFMethodGuess ) cycle
      select case( pc%ChemPotMethod )

      ! Chemical potential by gradual insertion
      case( ChemPotMethodGradIns )
        if( (((pc%GradInsInit .eq. 0) .or. (Step .gt. pc%GradInsInit)) .and. GradInsInitialization) .or. (pc%WFMethod .ne. WFMethodGuess)) cycle componentLoop

        ! determine, if chemical potential has to be calculated
        pc%CalcChemPot = .false.

        if( GradInsFrequency > 0 .and. mod( Step, GradInsFrequency ) == 0 ) pc%CalcChemPot = .true.

        if (GradInsInitialization) pc%CalcChemPot = .true.

#if MPI_VER > 0
        ! Per Process we calculate GI only for one component
        if (mod(NProc,this%NGradInsComp)/= pc%NGradThis) pc%CalcChemPot = .false.
#endif

        if( pc%CalcChemPot ) then

          ! Set fluctuating particle
          ncf = pc%NFluctComp( pc%NFluctState )
          if( ncf == i ) then
            npf = rnd( pc%NPart )
            call Move2End( this, ncf, npf )
          else
            npf = 1
          end if

          ! Specify max amount of biased interaction partners
          call BiasedPartners(this, ncf, npf)

          ! Specify ratio of biased Moves
          ratio = (this%NGradIns * 50 / this%NDF )
          ndfmove = this%NDF
          ndfbiased = this%NDF * ratio
          ndffluct = this%NDF * 10
          ndfchange = this%NDF * 10
          ndfcp = ndfmove + ndfbiased + ndffluct + ndfchange
          pc%NState(:) = 1
          nstate = 0

giloop:   do j = 1, NFullFluct * ndfcp

            ! Update List for Biased Moves
            if (mod(j,ndfcp) == 0) call BiasedPartners(this, ncf, npf)

            ! Choose particle randomly
            s = 0
            sndf = 0
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

              nuh = 1 + ((s - r) - (np-1) * ndf)
loop2:        do nu = 1, this%Component(nc)%Molecule%NUnit
                if (nuh <= sum(this%Component(nc)%Molecule%Unit(1:nu)%NDF)) exit loop2
              end do loop2

              ! Move or rotate
              if( mod( s - r, ndf ) < 3 ) then
                call Move( this, nc, np, nu )
              else
                call Rotate( this, nc, np, nu )
              end if

            else if( r <= (ndfmove + ndfbiased) ) then

              nuh2 = 0
              r = (r - ndfmove - 1) / ratio + 1
              nuh = int( (r-1)* this%NGradIns / (ndfbiased/ratio) ) 
loop3:        do nc = 1, this%NComponents
                s = s + this%Component(nc)%BiasedPartners
                if( nuh < s ) exit loop3
                sndf = sndf + this%Component(nc)%BiasedPartners
              end do loop3
              ndf = this%Component(nc)%Molecule%NDF
              np = this%BiasedPartners(int((nuh-sndf)*this%Component(nc)%BiasedPartnersNum / this%Component(nc)%BiasedPartners)+1)
              nuh= int(( (r-1)/(ndfbiased/ratio)*this%NGradIns + 1 ) * ndf - nuh )
              do nu = 1, this%Component(nc)%Molecule%NUnit
                if (nuh <= sum(this%Component(nc)%Molecule%Unit(1:nu)%NDF)) exit
                nuh2 = nuh2 + this%Component(nc)%Molecule%Unit(nu)%NDF
              end do

              ! Acceleration of MC Moves
              if (np .gt. this%Component(nc)%NPart) cycle

              ! Move or rotate biased
              if( (mod( sndf - r, ndf)-nuh2) < 3 ) then
                call MoveBiased( this, nc, np, nu, ncf, npf )
              else
                call RotateBiased( this, nc, np, nu, ncf, npf )
              end if

            else if( r <= (ndfmove + ndfbiased + ndffluct) ) then
              r = r - ndfmove - ndfbiased
              ndf = this%Component(ncf)%Molecule%NDF

              np = 1 + (s - r) / ndf
              nuh = 1 + ((s - r) - (np-1) * ndf)
loop5:        do nu = 1, this%Component(ncf)%Molecule%NUnit
                if (nuh <= sum(this%Component(ncf)%Molecule%Unit(1:nu)%NDF)) exit loop5
              end do loop5


              ! Move or rotate fluctuating particle
              if( mod( r, ndf ) < 3 ) then
                call Move( this, ncf, npf, nu )
              else
                call Rotate( this, ncf, npf, nu )
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
          pc%ProbW1 = pc%ProbW1 + real( pc%NState(pc%NFluctMax), RK ) / pc%WF( pc%NFluctMax )
          pc%ProbW0V = pc%ProbW0V + real(pc%NState(0), RK) / this%Density
          pc%ProbW1Rho = pc%ProbW1Rho + real( pc%NState(pc%NFluctMax), RK ) / pc%WF(pc%NFluctMax) * this%Density

          ! Calculate chemical potential
          ! (long range correction already done in ChangeFluct)
          pc%ChemPot = pc%ProbW0 / pc%ProbW1Rho
          pc%ChemPot1 = pc%ProbW0 / pc%ProbW1
          pc%ChemPot2 = pc%ProbW0V / pc%ProbW1

        else
!           pc%CalcChemPot = .false.
          pc%ChemPot = 0._RK
        end if

         if( mod( Step, ErrorsUpdateFrequency ) == 0 .or. GradInsInitialization ) then
          ! Here we sum up the NStateWF over all processes 
          ! dealing with a specific component to improve statistics
#if MPI_VER > 0
          call MPI_Allreduce(pc%NStateWF,tempVec(0:pc%NFluctMax),size(pc%NStateWF),MPI_INTEGER, MPI_SUM, Communicator, ierror)
          pc%NStateWF = tempVec(0:pc%NFluctMax)
          call MPI_Reduce( pc%NFluctUpSuccesses(:),tempVec1(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
          call MPI_Reduce( pc%NFluctUpAttempts(:),tempVec2(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER,  MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
          call MPI_Reduce( pc%NFluctDownSuccesses(:),tempVec3(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
          call MPI_Reduce( pc%NFluctDownAttempts(:),tempVec4(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, MPI_SUM, NRootProc, MPI_COMM_WORLD, ierror )
          
           do j = 1, pc%NFluctMax
            pc%WF(j) = pc%WF(j) * real(pc%NStateWF(0) + 1, RK) / real(pc%NStateWF(j) + 1, RK)
           end do
          write( IOBuffer, '("New weighting factors for ",A," calculated:")' ) trim( pc%PotModFileName )
          call LogWrite
          write( IOBuffer, '("   State      NState      new WF     up        down (%)")' )
          call LogWrite
          write( IOBuffer, '("   --------------------------------  --------  --------")' )
          call LogWrite
          j = pc%NFluctMax
          write( IOBuffer, '(I8, I12, F15.2, 2F10.4)' ) j, pc%NStateWF(j), pc%WF(j), 0._RK, real(tempVec3(j), RK) / &
&             real(tempVec4(j), RK) * 100._RK
            call LogWrite

          do j = pc%NFluctMax - 1, 1, -1
            write( IOBuffer, '(I8, I12, F15.2, 2F10.4)' ) j, pc%NStateWF(j), pc%WF(j), real(tempVec1(j+1), RK) / &
&               real(tempVec2(j+1), RK) * 100._RK, real(tempVec3(j), RK) / real(tempVec4(j), RK) * 100._RK
            call LogWrite
          end do
          write( IOBuffer, '(I8, I12, F15.2, 2F10.4)' ) 0, pc%NStateWF(0), pc%WF(0), real(tempVec1(1), RK) / &
&             real(tempVec2(1), RK) * 100._RK, 0._RK
          call LogWrite
          call LogWriteBlank
          pc%NStateWF(:) = 0
          
#else        
          do j = 1, pc%NFluctMax
            pc%WF(j) = pc%WF(j) * real(pc%NStateWF(0) + 1, RK) / real(pc%NStateWF(j) + 1, RK)
          end do
          write( IOBuffer, '("New weighting factors for ",A," calculated:")' ) trim( pc%PotModFileName )
          call LogWrite
          write( IOBuffer, '("   State      NState      new WF     up        down (%)")' )
          call LogWrite
          write( IOBuffer, '("   --------------------------------  --------  --------")' )
          call LogWrite
          j = pc%NFluctMax
          write( IOBuffer, '(I8, I12, F15.2, 2F10.4)' ) j, pc%NStateWF(j), pc%WF(j), 0._RK, real(pc%NFluctDownSuccesses(j), RK) / &
&             real(pc%NFluctDownAttempts(j), RK) * 100._RK
            call LogWrite

          do j = pc%NFluctMax - 1, 1, -1
            write( IOBuffer, '(I8, I12, F15.2, 2F10.4)' ) j, pc%NStateWF(j), pc%WF(j), real(pc%NFluctUpSuccesses(j+1), RK) / &
&               real(pc%NFluctUpAttempts(j+1), RK) * 100._RK, real(pc%NFluctDownSuccesses(j), RK) / &
&               real(pc%NFluctDownAttempts(j), RK) * 100._RK
            call LogWrite
          end do

          write( IOBuffer, '(I8, I12, F15.2, 2F10.4)' ) 0, pc%NStateWF(0), pc%WF(0), &
&           real(pc%NFluctUpSuccesses(1), RK) / real(pc%NFluctUpAttempts(1), RK) * 100._RK, 0._RK
          call LogWrite
          call LogWriteBlank
          pc%NStateWF(:) = 0
#endif         
        end if


      ! Chemical potential by Widom's test particle method
      ! Just applicable for ReactionField Method. 
      ! - otherwise you cannot calculate the energy of only 1 particle
      ! Check earlier in ms2_ensemble (right after component construction

      case( ChemPotMethodWidom )

        call Unit2AtomTest( pc, pc%NTest, pc%Molecule%NUnit )

#if MPI_VER > 0
        if ( (SimulationType .ne. MonteCarlo) .or. (Equilibration .and. CommonEqui) ) then
          this%EPotTest(:) = 0._RK
          this%EPotTest(pc%NTest0:pc%NTest2) = this%Density * pc%EPotTestCorrLJ + pc%EPotTestCorrRF
        else
          this%EPotTest(:) = this%Density * pc%EPotTestCorrLJ + pc%EPotTestCorrRF
        end if
#else
        this%EPotTest(:) = this%Density * pc%EPotTestCorrLJ + pc%EPotTestCorrRF
#endif
        this%EPotTest(1:pc%NTest) = this%EPotTest(1:pc%NTest) + pc%EPotTestIntra(1:pc%NTest) ! EPotTest can be longer than Intra-Array
        do j = 1, this%NRealComponents
          call ChemicalPotential( this%Interaction( i, j ), this%EPotTest, this%BoxLength )
        end do
#if MPI_VER > 0
        if ( (SimulationType .ne. MonteCarlo) .or. (Equilibration .and. CommonEqui) ) then
          call MPI_Reduce( this%EPotTest, EPotTest, this%NTestMax, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
          this%EPotTest = EPotTest
        end if
#endif

        ChemPot = sum( exp( -( this%EPotTest(:) ) / this%Temperature ) ) / pc%NTest


! #if MPI_VER > 0
!         if ( SimulationType .ne. MonteCarlo .or. (Equilibration .and. CommonEqui) ) then
!           call MPI_Bcast( this%Density, 1, MPI_RK, NRootProc, Communicator, ierror )
!           call MPI_Bcast( this%EPot, 1, MPI_RK,  NRootProc, Communicator, ierror )
!         endif
! #endif

        ! partial molar enthalpy
       HW_H_local = this%EPot + this%RefPressure / this%Density * real( this%NPart, RK )
       HW_V_local = (1.0 / this%Density) * this%NPart
       HW_denom_local = 0
       HW_counter_local = 0

       do j=1, pc%NTest 
          HW_counter_local = HW_counter_local + ( HW_H_local + this%EPotTest(j) ) * &
&                            exp( - this%EPotTest(j) / this%RefTemperature )
          HW_denom_local = HW_denom_local + exp( - this%EPotTest(j) / this%RefTemperature )
       end do

       HW_counter_local = HW_V_local * HW_counter_local / pc%NTest
       HW_denom_local = HW_V_local * HW_denom_local / pc%NTest

! #if MPI_VER > 0
!         if ( SimulationType .ne. MonteCarlo .or. (Equilibration .and. CommonEqui) ) then
!           ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
!           call MPI_Reduce( ChemPot, pc%ChemPot, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
!           call MPI_Reduce( HW_counter_local, pc%HW_counter, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
!           call MPI_Reduce( HW_denom_local, pc%HW_denom, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
!             pc%ChemPot = pc%ChemPot/NProcs
!             pc%HW_counter = pc%HW_counter/NProcs
!             pc%HW_denom = pc%HW_denom/NProcs 
!         else
!             pc%ChemPot = ChemPot
!             pc%HW_counter = HW_counter_local
!             pc%HW_denom = HW_denom_local
!         endif
! #else
        pc%ChemPot = ChemPot
        pc%HW_counter = HW_counter_local
        pc%HW_denom = HW_denom_local
! #endif

#if OSMOP == 2
        ! Profile of  the chemical potential 
        ! Initialize local arrays
        if (SimulationType .eq. MolecularDynamics ) then
          ChemPotProfile(:) = 0._RK

          do j=1, NBinsDen

            NTestBinDen = 0
  
            do m=1, pc%NTest ! Michael Sch.: for IDF: P0Test(.,1,1) is not correct here...should be either Pm0 or unit specific!
              if (this%P0Test(m,1,1) .ge. real(j-1)/NBinsDen) then
                if (this%P0Test(m,1,1) < real(j)/NBinsDen) then
                  ChemPotProfile(j) = ChemPotProfile(j) + (exp( -( this%EPotTest(m)) /  this%Temperature))
                  NTestBinDen = NTestBinDen + 1
                end if
              end if
            end do

#if MPI_VER > 0
            call MPI_Allreduce(NTestBinDen,NTestBinDenAll,1,MPI_INTEGER, MPI_SUM, NRootProc, Communicator, ierror)
            if ( NTestBinDenAll > 1) then
              ChemPotProfile(j)= ChemPotProfile(j) / NTestBinDenAll
            end if
            call MPI_Reduce( ChemPotProfile(j), pc%ChemPotProfile(j), 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
#else
            if ( NTestBinDen > 1) then
              ChemPotProfile(j)= ChemPotProfile(j) / NTestBinDen
            end if 
            pc%ChemPotProfile(j) = ChemPotProfile(j)
#endif

          end do !NBinsDen
        end if
#endif

      case( ChemPotMethodThermoInt )

        if ( mod(Step,pc%changeLaFreq)==0 ) then
          call ChangeLambda( this, t, i )
        end if
        ! Calculating the energy of the fluctuating particle
        if (mod(Step,pc%changeLaFreq) .ge. pc%forfeitLaSampl ) then
          pc%currentBinsEn = (this%Density * pc%EPotTestCorrLJ + pc%EPotTestCorrRF)*this%Component(t)%Lambda**pc%LambdaExponent
          if (SimulationType .ne. MolecularDynamics ) then
            pc%currentBinsEn = pc%currentBinsEn + GetEnergy( this, t, 1 ) - GetEnergyIntra( this, t, 1 )
          else
            E = 0._RK; EIntra = 0._RK; EBond = 0._RK; EAngle = 0._RK; EDihedral = 0._RK; F(:,:) = 0._RK
            nu =this%Component(t)%Molecule%NUnit
            do j = 1, this%NComponents
              if (j > t) then
                call MDEnergy( this%Interaction(t,j), nu, F(:,1:nu), E, EIntra, EBond, EAngle, EDihedral, this%BoxLength, .true. )
              else
                call MDEnergy( this%Interaction(j,t), nu, F(:,1:nu), E, EIntra, EBond, EAngle, EDihedral, this%BoxLength, .false. )
              end if
            end do
            E = E - EIntra
#if MPI_VER > 0
            call MPI_Reduce( E, EBin, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
#else
            EBin = E
#endif
            pc%currentBinsEn = pc%currentBinsEn + EBin
          end if
        end if

        ! chemPot with LambdaMin by Widom
        call Unit2AtomTest( pc, pc%NTest, pc%Molecule%NUnit )
#if MPI_VER > 0
        if ( (SimulationType .ne. MonteCarlo) .or. (Equilibration .and. CommonEqui) ) then
          this%EPotTest(:) = 0._RK
          this%EPotTest(pc%NTest0:pc%NTest2) = this%Density * pc%EPotTestCorrLJ + pc%EPotTestCorrRF
        else
          this%EPotTest(:) = this%Density * pc%EPotTestCorrLJ + pc%EPotTestCorrRF
        end if
#else
        this%EPotTest(:) = this%Density * pc%EPotTestCorrLJ + pc%EPotTestCorrRF
#endif
        this%EPotTest(1:pc%NTest) = this%EPotTest(1:pc%NTest) + pc%EPotTestIntra(1:pc%NTest) ! EPotTest can be longer than Intra-Array
        do j = 1, this%NComponents
          call ChemicalPotential( this%Interaction( i, j ), this%EPotTest, this%BoxLength )
        end do
#if MPI_VER > 0
        if ( (SimulationType .ne. MonteCarlo) .or. (Equilibration .and. CommonEqui) ) then
          call MPI_Reduce( this%EPotTest, EPotTest, this%NTestMax, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
          this%EPotTest = EPotTest
        end if
#endif
        factor=pc%LaMin**pc%LambdaExponent
        ExpMinusBetaEnLaMin = sum( exp( -( factor*this%EPotTest(:) ) / this%Temperature ) ) / pc%NTest

        ! partial molar enthalpy
       HW_H_local = this%EPot + ( this%RefPressure / this%Density ) * real( this%NPart, RK )
       HW_V_local = (1.0 / this%Density) * this%NPart
       HW_denom_local = 0
       HW_counter_local = 0

       do j=1, pc%NTest 
          HW_counter_local= HW_counter_local + HW_V_local * ( HW_H_local + factor*this%EPotTest(j) ) * exp( - factor*this%EPotTest(j) / this%RefTemperature )
          HW_denom_local = HW_denom_local + HW_V_local * exp( - factor*this%EPotTest(j) / this%RefTemperature )
       end do

       HW_counter_local = HW_counter_local / pc%NTest
       HW_denom_local = HW_V_local * HW_denom_local / pc%NTest

! #if MPI_VER > 0
!         if ( SimulationType .eq. MolecularDynamics  ) then
!           ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
!           call MPI_Reduce( ExpMinusBetaEnLaMin, pc%ExpMinusBetaEnLaMin, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
!           call MPI_Reduce( HW_counter_local, pc%HW_counter, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
!           call MPI_Reduce( HW_denom_local, pc%HW_denom, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
!           pc%ExpMinusBetaEnLaMin = pc%ExpMinusBetaEnLaMin/NProcs
!           pc%HW_counter = pc%HW_counter/NProcs
!           pc%HW_denom = pc%HW_denom/NProcs
!         else
!           pc%ExpMinusBetaEnLaMin = ExpMinusBetaEnLaMin
!           pc%HW_counter = HW_counter_local
!           pc%HW_denom = HW_denom_local
!         endif
! #else
        pc%ExpMinusBetaEnLaMin = ExpMinusBetaEnLaMin
        pc%HW_counter = HW_counter_local
        pc%HW_denom = HW_denom_local
! #endif
        ! end of Widom for ThermoInt with LambdaMin

        t=t+1 ! end of ThermoInt

      case default

      end select

    end do  componentLoop

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
        n1 = pi%NPart1 * pi%NUnit1
        n2 = pi%NPart2 * pi%NUnit2
        pi%EPot(1:n1, 1:n2) = pi%EPotNew(1:n1, 1:n2)
        pi%d2EpotdV2(1:n1, 1:n2) = pi%d2EpotdV2New(1:n1, 1:n2)
        if ( this%OptPressure ) then
          pi%Virial(1:n1, 1:n2) = pi%VirialNew(1:n1, 1:n2)
        end if
        if ( j == i .and. UseIntDegFreed) then
          pi%EPotAngle(:) = pi%EPotAngleNew(:)
          pi%EPotTo(:) = pi%EPotToNew(:)
        end if
      end do
    end do

  end subroutine TEnsemble_UpdateEnergy



!==============================================================!
!  Subroutine TEnsemble_UpdateEnergy1                          !
!==============================================================!

  subroutine TEnsemble_UpdateEnergy1( this, nc, np, nu )

    implicit none

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np, nu

    ! Declare local variables
    type(TInteraction), pointer :: pi
    integer                     :: n
    integer                     :: i
    integer                     :: npu
    integer                     :: NBond, NAngle, NDihedral

    npu = (np-1) * this%Component(nc)%Molecule%NUnit + nu

    ! Update potential energy and virial matrices for a particle
    do i = 1, this%NComponents
      pi => this%Interaction(nc, i)
      n = pi%NPart2 * pi%NUnit2
      pi%EPot(npu, 1:n) = pi%EPot1(1:n)
      pi%d2EpotdV2(npu, 1:n) = pi%d2EpotdV21(1:n)

      if ( this%OptPressure ) then
        pi%Virial(npu, 1:n) = pi%Virial1(1:n)
      end if

      this%Interaction(i, nc)%EPot(1:n, npu) = pi%EPot1(1:n)
      this%Interaction(i, nc)%d2EpotdV2(1:n, npu) = pi%d2EpotdV21(1:n)

      if ( this%OptPressure ) then
        this%Interaction(i, nc)%Virial(1:n, npu) = pi%Virial1(1:n)
      end if
    end do

    if ( UseIntDegFreed ) then
      pi => this%Interaction(nc,nc)
      NAngle = pi%NAngle
      NDihedral = pi%NDihedral

      pi%EPotAngle((np-1)*NAngle+1:np*NAngle) = pi%EPot1Angle(:)
      pi%EPotTo((np-1)*NDihedral+1:np*Ndihedral) = pi%EPot1To(:)
    end if

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
    integer                     :: nu1, nu
    integer                     :: i, n
    real(RK)                    :: Intra

    ! Initialize new energy
    E = 0._RK
    Intra = 0._RK

    if (LongRange .eq. Ewald) then
      call EwaldSelfTerm_Energy ( this )

#if SPME > 0
    else if (LongRange .eq. PME) then
      call PMESelfTermMC ( this )
#endif
    end if

    ! Loop over components
    do nc = 1, this%NComponents
      do i = 1, this%NComponents
        pi => this%Interaction(nc, i)
        n = pi%NUnit2*pi%NPart2

        ! Loop over units
        do np = 1, this%Component(nc)%NPart
          do nu=1, this%Component(nc)%Molecule%NUnit

            nu1=(np-1)*pi%NUnit1+nu ! global number of unit

            call Energy( pi, np, nu, this%BoxLength )

            if ( pi%SameComponent .and. UseIntDegFreed ) then
              call IntraEnergy(pi, np, nu, this%BoxLength)
              pi%EPotAngleNew((np-1)*pi%NAngle+1:np*pi%NAngle) = pi%EPot1Angle(:)
              pi%EPotToNew((np-1)*pi%NDihedral+1:np*pi%NDihedral) = pi%EPot1To(:)
            end if

            ! Save new energy matrix
            pi%EPotNew(nu1, 1:n) = pi%EPot1(1:n)
            pi%d2EpotdV2New(nu1, 1:n) = pi%d2EpotdV21(1:n)

            if (this%OptPressure) then
              pi%VirialNew(nu1, 1:n) = pi%Virial1(1:n)
            end if

            ! Sum energy
            E = E + sum( pi%EPot1(1:n) )
          end do
        end do
      end do
      Intra = Intra + sum(this%Interaction(nc,nc)%EPotAngleNew(:)) + &
&                     sum(this%Interaction(nc,nc)%EPotToNew(:))
    end do

    ! Calculate new energy
    E = .5_RK * E + this%Density * this%EPotCorrLJ + this%EPotCorrRF + Intra

! Ewald 
    if (LongRange .eq. Ewald) then
      call EwaldFourierEnergy(this)
      E = E + this%UFourier + this%UIntra + this%USelbstTerm
#if SPME > 0

    else if (LongRange .eq. PME) then
      call charge_grid_MCall ( this )
      call PMEFourierTermMC ( this )
      E = E + this%UFourier + this%UIntra + this%USelbstTerm
#endif
    end if

  end subroutine TEnsemble_Energy



!==============================================================!
!  Subroutine TEnsemble_Energy1                                !
!==============================================================!

  subroutine TEnsemble_Energy1( this, nc, np, nu, EPotNew )

    implicit none

    ! Declare arguments
    type(TEnsemble)       :: this
    integer, intent(in)   :: nc, np, nu
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
      n = pi%NPart2 * pi%NUnit2

      call Energy( pi, np, nu, this%BoxLength )

      if ( pi%SameComponent .and. UseIntDegFreed ) then
        call IntraEnergy( pi, np,  nu, this%BoxLength )
        EPotNew = EPotNew + sum(pi%EPot1Angle) + sum(pi%EPot1To)
      end if

      ! Calculate new energy
      EPotNew = EPotNew + sum( pi%EPot1(1:n) )
    end do

    if (LongRange .eq. Ewald) then
       call EwaldFourierEnergy(this,nc,np)
       EPotNew = EPotnew + this%UFourier
#if SPME > 0

    else if (LongRange .eq. PME) then
       call PMEFourierTermMC ( this )
       EPotNew = EPotnew + this%UFourier
#endif
    end if
  end subroutine TEnsemble_Energy1


!==============================================================!
!  Subroutine TEnsemble_EwaldEnergy1                           !
!==============================================================!

  subroutine TEnsemble_EwaldEnergy1( this, nc, np, EPotNew, m )

    implicit none

    ! Declare arguments
    type(TEnsemble)       :: this
    integer, intent(in)   :: nc, np,m
    real(RK), intent(out) :: EPotNew

    ! Declare local variables
    type(TInteraction), pointer :: pi
    integer                     :: n
    integer                     :: i, nu

    ! Initialize new energy
    EPotNew = 0._RK

    ! Loop over components
    do nu=1,this%Component(nc)%Molecule%NUnit
      do i = 1, this%NComponents
        pi => this%Interaction(nc, i)
        n = pi%NPart2*pi%NUnit2
        call Energy( pi, np, nu, this%BoxLength )
        if ( pi%SameComponent .and. UseIntDegFreed ) then
          call IntraEnergy( pi, np, nu, this%BoxLength )
          EPotNew = EPotNew + sum(pi%EPot1Angle) + sum(pi%EPot1To)
        end if
        ! Calculate new energy
        EPotNew = EPotNew + sum( pi%EPot1(1:n) )
      end do
    end do

    if (LongRange .eq. Ewald) then
       call EwaldFourierEnergy(this,nc,np,m)
       EPotNew = EPotnew + this%UFourier
#if SPME > 0

    else if (LongRange .eq. PME) then
       call PMEFourierTermMC ( this )
       EPotNew = EPotnew + this%UFourier
#endif
    end if
  end subroutine TEnsemble_EwaldEnergy1




!==============================================================!
!  Subroutine TEnsemble_Energy_CF                              !
!==============================================================!

  subroutine TEnsemble_Energy1_CF( this, nc, np, ncold, npold, EPotNew )

    implicit none

    ! Declare arguments
    type(TEnsemble)       :: this
    integer, intent(in)   :: nc, np
    integer, intent(in)   :: ncold, npold
    real(RK), intent(out) :: EPotNew

    ! Declare local variables
    type(TInteraction), pointer :: pi
    integer                     :: n
    integer                     :: i, nu

    ! Initialize new energy
    EPotNew = 0._RK

    ! Loop over components
    do i = 1, this%NComponents
      pi => this%Interaction(nc, i)
      n = pi%NPart2*pi%NUnit2
      do nu = 1,pi%NUnit1

        call Energy( pi, np, nu, this%BoxLength )

        if ( pi%SameComponent .and. UseIntDegFreed ) then
          call IntraEnergy( pi, np, nu, this%BoxLength )
          EPotNew = EPotNew + sum(pi%EPot1Angle) + sum(pi%EPot1To)
        end if

        ! Calculate new energy
        EPotNew = EPotNew + sum( pi%EPot1(1:n) )
      end do
    end do

    if (LongRange .eq. Ewald) then
       call EwaldFourierEnergy(this,nc,np,ncold,npold)
       EPotNew = EPotnew + this%UFourier + this%USelbstTerm + this%UIntra
#if SPME > 0

    else if (LongRange .eq. PME) then
       call PMEFourierTermMC ( this )
       EPotNew = EPotnew + this%UFourier
#endif
    end if
  end subroutine TEnsemble_Energy1_CF



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
    integer :: n, n2
    real(RK):: Intra

    ! Calculate potential energy of a particle
    E = 0._RK
    Intra = 0._RK
    do i = 1, this%NComponents
      n = this%Component(i)%NPart*this%Component(i)%Molecule%NUnit
      do j = 1, this%NComponents
        n2 = this%Component(j)%NPart*this%Component(j)%Molecule%NUnit
        E = E + sum( this%Interaction(j, i)%EPot(1:n2, 1:n) )
      end do
      ! Kein Faktor 2, weil unten einfach aufaddiert wird
      Intra = Intra + sum(this%Interaction(i,i)%EPotAngle(:)) + &
&                   sum(this%Interaction(i,i)%EPotTo(:))
    end do
    E = .5_RK * E + this%Density * this%EPotCorrLJ + this%EPotCorrRF + Intra

! Ewald 
    if (LongRange .eq. Ewald) then
      call EwaldFourierEnergy(this)
      E = E + this%UFourier + this%UIntra + this%USelbstTerm
#if SPME > 0

    else if (LongRange .eq. PME) then
      call charge_grid_MCall (this)
      call PMEFourierTermMC(this)
      E = E + this%UFourier + this%UIntra + this%USelbstTerm
#endif
    end if

  end function TEnsemble_GetEnergy



!==============================================================!
!  Function TEnsemble_GetEnergy1                               !
!==============================================================!

  function TEnsemble_GetEnergy1( this, nc, np, nu ) result(E)

    implicit none

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np, nu

    ! Declare result
    real(RK) :: E

    ! Declare local variables
    type(TInteraction), pointer :: pi
    integer :: i
    integer :: NAngle, NDihedral
    integer :: NUnitPart
    integer :: nup1

    ! Calculate potential energy of a particle
    E = 0._RK
    nup1= this%Component(nc)%Molecule%NUnit * (np - 1) + nu
    do i = 1, this%NComponents
      NUnitPart = this%Component(i)%Molecule%NUnit * this%Component(i)%NPart
      E = E + sum( this%Interaction(i, nc)%EPot(1:NUnitPart, nup1) )
    end do

    if ( UseIntDegFreed ) then
      pi => this%Interaction(nc,nc)
      NAngle = pi%NAngle
      NDihedral = pi%NDihedral
      do i=1,NAngle
        if (this%Component(nc)%molecule%idfangle(i)%unitid1 == nu .or. &
&           this%Component(nc)%molecule%idfangle(i)%unitid2 == nu .or. &
&           this%Component(nc)%molecule%idfangle(i)%unitid3 == nu) then
          E = E + pi%EPotAngle((np-1)*NAngle+i)
        end if
      end do
      do i=1,NDihedral
        if (this%Component(nc)%molecule%idfdihedral(i)%unitid1 == nu .or. &
&           this%Component(nc)%molecule%idfdihedral(i)%unitid2 == nu .or. &
&           this%Component(nc)%molecule%idfdihedral(i)%unitid3 == nu .or. &
&           this%Component(nc)%molecule%idfdihedral(i)%unitid4 == nu) then
          E = E + pi%EPotTo((np-1)*NDihedral+i)
        end if
      end do
    end if

! Ewald
    if (LongRange .eq. Ewald) then
      E = E + this%UFourier
#if SPME > 0
    else if (LongRange .eq. PME) then
      E = E + this%UFourier
#endif
    end if

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
      n = this%Component(i)%NPart * this%Component(i)%Molecule%NUnit
      do j = 1, this%NComponents
        V = V + sum( this%Interaction(j, i)%Virial(1:this%Component(j)%NPart * this%Component(j)%Molecule%NUnit, 1:n) )
      end do
    end do
    V = .5_RK * V + this%Density * this%VirialCorrLJ + this%VirialCorrRF*this%Volume0

    if (LongRange .eq. Ewald) then
!       call EwaldFourierEnergy(this)
       V = V + this%EVirial
#if SPME > 0
    else if (LongRange .eq. PME) then
       V = V + this%EVirial
#endif
    end if

  end function TEnsemble_GetVirial

!==============================================================!
!  Function TEnsemble_Getd2EpotdV2                                !
!==============================================================!

  function TEnsemble_Getd2EpotdV2( this ) result(V)

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
      n = this%Component(i)%NPart * this%Component(i)%Molecule%NUnit
      do j = 1, this%NComponents
        V = V + sum( this%Interaction(j, i)% &
&         d2EpotdV2(1:this%Component(j)%NPart * this%Component(j)%Molecule%NUnit, 1:n) )
      end do
    end do
    V = .5_RK * V + this%Density * this%d2EpotdV2CorrLJ

  end function TEnsemble_Getd2EpotdV2


!==============================================================!
!  Subroutine TEnsemble_Move                                   !
!==============================================================!

  subroutine TEnsemble_Move( this, nc, np, nu )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np, nu

    ! Declare local variables
    real(RK)                  :: r(3)
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier
#if SPME > 0
    real(RK)                  :: EVirial
#endif
    real(RK)                  :: EPotDelta
    type(TComponent), pointer :: pc
    integer                   :: i
    logical                   :: accepted

    ! Assign local variables
    pc => this%Component(nc)

    ! Update number of move attempts
    pc%NMoveAttempts = pc%NMoveAttempts + 1

    ! Save current particle position and energy
    r(:) = pc%P0(np, :, nu)
    EPotOld = GetEnergy( this, nc, np, nu )

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO

#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      call chargegrid_min  (this, nc, np)
#endif
    end if

    ! Generate a trial displacement
    do i = 1, 3
      pc%P0(np, i, nu) = pc%P0(np, i, nu) + rnd( -pc%DispTran, pc%DispTran )
    end do

    !!! Implement fixed bond length here. 
    ! 1. check for bond partners and if forconst=0
    ! 2. stretch bond to original, calculate scaling factor
    ! 3. move all atoms which are connected to the displaced atom in any number of chains

    ! Apply periodic boundary conditions
    pc%P0(np, :, nu) = pc%P0(np, :, nu) - anint( pc%P0(np, :, nu) )

    ! Convert unit coordinates to atom positions
    call Unit2Atom1( pc, np )

#if SPME > 0
    ! Calculate changes in the SPME grid
    if (LongRange .eq. PME) then
      call chargegrid_plus (this, nc, np)
    end if
#endif

    ! Calculate particle energy at trial position
    call Energy( this, nc, np, nu, EPotNew )
    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, MPI_RK, MPI_SUM, Communicator, ierror )

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
      call UpdateEnergy( this, nc, np, nu )

      ! Calculate new COM
      call Unit2Mol( pc, np )
      pc%Pm0(np, :)    = pc%Pm0(np, :) - anint( pc%Pm0(np, :) )
    else

      ! Reject move
      if (LongRange .eq. Ewald) then
          this%UFourier = EFourier

          DO i=1,pc%Molecule%NCharge
            this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
            this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
            this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
          END DO

          pc%P0(np, :, nu) = r(:)
          call Unit2Atom1( pc, np, nu )
          call EwaldFourierEnergy(this,nc,np)

#if SPME > 0
      else if (LongRange .eq. PME) then
          this%UFourier = EFourier
          this%EVirial  = EVirial
          call chargegrid_min  (this, nc, np)
          pc%P0(np, :, nu) = r(:)
          call Unit2Atom1( pc, np, nu )
          call chargegrid_plus (this, nc, np)
#endif
      else
          pc%P0(np, :, nu) = r(:)
          call Unit2Atom1( pc, np, nu )
      end if
    end if

  end subroutine TEnsemble_Move



!==============================================================!
!  Subroutine TEnsemble_Rotate                                 !
!==============================================================!

  subroutine TEnsemble_Rotate( this, nc, np, nu )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments03
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np, nu

    ! Declare local variables
    real(RK)                  :: q(4), dq(3)
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier
#if SPME > 0
    real(RK)                  :: EVirial
#endif
    type(TComponent), pointer :: pc
    integer                   :: i
    real(RK)                  :: EPotDelta
    logical                   :: accepted

    ! Assign local variables
    pc => this%Component(nc)

    ! Update number of rotation attempts
    pc%NRotateAttempts = pc%NRotateAttempts + 1

    ! Save current particle orientation and energy
    q(:) = pc%Q0(np, :, nu)
    EPotOld = GetEnergy( this, nc, np, nu )

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO

#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      call chargegrid_min  (this, nc, np)
#endif
    end if

    ! Generate a trial rotation
    do i = 1, 3
      dq(i) = rnd( -pc%DispRot, pc%DispRot )
    end do
    ! rotate unit 
    pc%Q0(np, 1, nu) = q(1) - dq(1) * q(2) - dq(2) * q(3) - dq(3) * q(4)
    pc%Q0(np, 2, nu) = q(2) + dq(1) * q(1) - dq(3) * q(3) + dq(2) * q(4)
    pc%Q0(np, 3, nu) = q(3) + dq(2) * q(1) + dq(3) * q(2) - dq(1) * q(4)
    pc%Q0(np, 4, nu) = q(4) + dq(3) * q(1) - dq(2) * q(2) + dq(1) * q(3)

    ! Convert unit coordinates to atom positions
    call Unit2Atom1( pc, np, nu )

#if SPME > 0
    if (LongRange .eq. PME) then
      call chargegrid_plus (this, nc, np)
    end if
#endif

    ! Calculate particle energy with trial orientation
    call Energy( this, nc, np, nu, EPotNew )

    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, MPI_RK, MPI_SUM, Communicator, ierror )
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
      call UpdateEnergy( this, nc, np, nu )

    else

      ! Reject move
      if (LongRange .eq. Ewald) then
        this%UFourier = EFourier
        DO i=1,pc%Molecule%NCharge
          this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
          this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
          this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
        END DO
        pc%Q0(np, :, nu) = q(:)
        call Unit2Atom1( pc, np, nu )
        call EwaldFourierEnergy(this,nc,np)

#if SPME > 0
      else if (LongRange .eq. PME) then
        this%UFourier = EFourier
        this%EVirial  = EVirial
        call chargegrid_min  (this, nc, np)
        pc%Q0(np, :, nu) = q(:)
        call Unit2Atom1( pc, np, nu )
        call chargegrid_plus (this, nc, np)
#endif

      else
        pc%Q0(np, :, nu) = q(:)
        call Unit2Atom1( pc, np, nu )
      end if

    end if

  end subroutine TEnsemble_Rotate

!==============================================================!
!  Subroutine TEnsemble_Move_NVE                               !
!==============================================================!

  subroutine TEnsemble_Move_NVE( this, nc, np, nu )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np, nu

    ! Declare local variables
    real(RK)                  :: r(3), rm(3)
    real(RK)                  :: EPotOld, EPotNew, NewOmega
    real(RK)                  :: EFourier
#if SPME > 0
    real(RK)                  :: EVirial
#endif
    real(RK)                  :: EPotDelta
    type(TComponent), pointer :: pc
    integer                   :: i

    ! Assign local variables
    pc => this%Component(nc)

    ! Update number of move attempts
    pc%NMoveAttempts = pc%NMoveAttempts + 1

    ! Save current particle position and energy
    r(:) = pc%P0(np, :, nu)
    rm(:) = pc%Pm0(np, :)
    EPotOld = GetEnergy( this, nc, np, nu )

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO

#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      call chargegrid_min  (this, nc, np)
#endif
    end if

    ! Generate a trial displacement
    do i = 1, 3
      pc%P0(np, i, nu) = pc%P0(np, i, nu) + rnd( -pc%DispTran, pc%DispTran )
    end do

    !!! Implement fixed bond length here. 
    ! 1. check for bond partners and if forconst=0
    ! 2. stretch bond to original, calculate scaling factor
    ! 3. move all atoms which are connected to the displaced atom in any number of chains
    
    ! Calculate new COM
    call Unit2Mol( pc, np )

    ! Apply periodic boundary conditions
    pc%P0(np, :, nu) = pc%P0(np, :, nu) - anint( pc%P0(np, :, nu) )
    pc%Pm0(np, :)    = pc%Pm0(np, :) - anint( pc%Pm0(np, :) )

    ! Convert unit coordinates to atom positions
    call Unit2Atom1( pc, np, nu )

#if SPME > 0
    ! Calculate changes in the SPME grid
    if (LongRange .eq. PME) then
      call chargegrid_plus (this, nc, np)
    end if
#endif

    ! Calculate particle energy at trial position
    call Energy( this, nc, np, nu, EPotNew )
    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, MPI_RK, MPI_SUM, Communicator, ierror )

    else
          EPotDelta = EPotOld - EPotNew
    endif
#else
     EPotDelta = EPotOld - EPotNew
#endif

    if( (this%RefHamiltonian*this%NPart - this%Epot+EPotDelta) < 0._RK ) then 
      NewOmega = 0._RK
    else 
      NewOmega = 1._RK
    end if

    if( ((this%RefHamiltonian*this%NPart - this%Epot+EPotDelta)/(this%RefHamiltonian*this%NPart - this%Epot))**((real (this%NDF-this%constrNDF, RK)-2._RK)/2._RK) &
&         * NewOmega .ge. rnd( 0._RK, 1._RK ) ) then
      ! Accept move
      pc%NMoveSuccesses = pc%NMoveSuccesses + 1
      call UpdateEnergy( this, nc, np, nu )
    else

      ! Reject move
      if (LongRange .eq. Ewald) then
          this%UFourier = EFourier

          DO i=1,pc%Molecule%NCharge
            this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
            this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
            this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
          END DO

          pc%P0(np, :, nu) = r(:)
          pc%Pm0(np, :) = rm(:)
          call Unit2Atom1( pc, np, nu )
          call EwaldFourierEnergy(this,nc,np)

#if SPME > 0
      else if (LongRange .eq. PME) then
          this%UFourier = EFourier
          this%EVirial  = EVirial
          call chargegrid_min  (this, nc, np)
          pc%P0(np, :, nu) = r(:)
          pc%Pm0(np, :) = rm(:)
          call Unit2Atom1( pc, np, nu )
          call chargegrid_plus (this, nc, np)
#endif
      else
          pc%P0(np, :, nu) = r(:)
          pc%Pm0(np, :) = rm(:)
          call Unit2Atom1( pc, np, nu )
      end if
    end if

  end subroutine TEnsemble_Move_NVE



!==============================================================!
!  Subroutine TEnsemble_Rotate_NVE                             !
!==============================================================!

  subroutine TEnsemble_Rotate_NVE( this, nc, np, nu )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments03
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np, nu

    ! Declare local variables
    real(RK)                  :: q(4), dq(3)
    real(RK)                  :: EPotOld, EPotNew, NewOmega
    real(RK)                  :: EFourier
#if SPME > 0
    real(RK)                  :: EVirial
#endif
    type(TComponent), pointer :: pc
    integer                   :: i
    real(RK)                  :: EPotDelta
    logical                   :: accepted

    ! Assign local variables
    pc => this%Component(nc)

    ! Update number of rotation attempts
    pc%NRotateAttempts = pc%NRotateAttempts + 1

    ! Save current particle orientation and energy
    q(:) = pc%Q0(np, :, nu)
    EPotOld = GetEnergy( this, nc, np, nu )

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO

#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      call chargegrid_min  (this, nc, np)
#endif
    end if

    ! Generate a trial rotation
    do i = 1, 3
      dq(i) = rnd( -pc%DispRot, pc%DispRot )
    end do
    pc%Q0(np, 1, nu) = q(1) - dq(1) * q(2) - dq(2) * q(3) - dq(3) * q(4)
    pc%Q0(np, 2, nu) = q(2) + dq(1) * q(1) - dq(2) * q(4) + dq(3) * q(3)
    pc%Q0(np, 3, nu) = q(3) + dq(1) * q(4) + dq(2) * q(1) - dq(3) * q(2)
    pc%Q0(np, 4, nu) = q(4) - dq(1) * q(3) + dq(2) * q(2) + dq(3) * q(1)


    ! Convert molecular coordinates to atom positions
    call Unit2Atom1( pc, np, nu )

#if SPME > 0
    if (LongRange .eq. PME) then
      call chargegrid_plus (this, nc, np)
    end if
#endif

    ! Calculate particle energy with trial orientation
    call Energy( this, nc, np, nu, EPotNew )

    ! Apply acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, MPI_RK, MPI_SUM, Communicator, ierror )
    else
      EPotDelta = EPotOld - EPotNew
    endif

#else
    EPotDelta = EPotOld - EPotNew
#endif

    if( (this%RefHamiltonian*this%NPart - this%Epot+EPotDelta) < 0._RK ) then 
      NewOmega = 0._RK
    else 
      NewOmega = 1._RK
    end if

    if( ((this%RefHamiltonian*this%NPart - this%Epot+EPotDelta)/(this%RefHamiltonian*this%NPart - this%Epot))**((real (this%NDF-this%constrNDF, RK)-2._RK)/2._RK) * NewOmega .ge. rnd( 0._RK, 1._RK ) ) then
      ! Accept rotation
      pc%NRotateSuccesses = pc%NRotateSuccesses + 1
      call UpdateEnergy( this, nc, np, nu )

    else

      ! Reject move
      if (LongRange .eq. Ewald) then
        this%UFourier = EFourier
        DO i=1,pc%Molecule%NCharge
          this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
          this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
          this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
        END DO
        pc%Q0(np, :, nu) = q(:)
        call Unit2Atom1( pc, np, nu )
        call EwaldFourierEnergy(this,nc,np)

#if SPME > 0
      else if (LongRange .eq. PME) then
        this%UFourier = EFourier
        this%EVirial  = EVirial
        call chargegrid_min  (this, nc, np)
        pc%Q0(np, :, nu) = q(:)
        call Unit2Atom1( pc, np, nu )
        call chargegrid_plus (this, nc, np)
#endif

      else
        pc%Q0(np, :, nu) = q(:)
        call Unit2Atom1( pc, np, nu )
      end if

    end if

  end subroutine TEnsemble_Rotate_NVE
  
!==============================================================!
!  Subroutine TEnsemble_Move_NPH                               !
!==============================================================!

  subroutine TEnsemble_Move_NPH( this, nc, np, nu )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np, nu

    ! Declare local variables
    real(RK)                  :: r(3), rm(3)
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier, EVirial
    real(RK)                  :: EPotDelta
    type(TComponent), pointer :: pc
    integer                   :: i
    logical                   :: accepted

    ! Assign local variables
    pc => this%Component(nc)

    ! Update number of move attempts
    pc%NMoveAttempts = pc%NMoveAttempts + 1

    ! Save current particle position and energy
    r(:) = pc%P0(np, :, nu)
    rm(:) = pc%Pm0(np, :)
    EPotOld = GetEnergy( this, nc, np, nu )

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO

#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      call chargegrid_min  (this, nc, np)
#endif
    end if

    ! Generate a trial displacement
    do i = 1, 3
      pc%P0(np, i, nu) = pc%P0(np, i, nu) + rnd( -pc%DispTran, pc%DispTran )
    end do

    !!! Implement fixed bond length here. 
    ! 1. check for bond partners and if forconst=0
    ! 2. stretch bond to original, calculate scaling factor
    ! 3. move all atoms which are connected to the displaced atom in any number of chains
    
    ! Calculate new COM
    call Unit2Mol( pc, np )

    ! Apply periodic boundary conditions
    pc%P0(np, :, nu) = pc%P0(np, :, nu) - anint( pc%P0(np, :, nu) )
    pc%Pm0(np, :)    = pc%Pm0(np, :) - anint( pc%Pm0(np, :) )

    ! Convert unit coordinates to atom positions
    call Unit2Atom1( pc, np, nu )

#if SPME > 0
    ! Calculate changes in the SPME grid
    if (LongRange .eq. PME) then
      call chargegrid_plus (this, nc, np)
    end if
#endif

    ! Calculate particle energy at trial position
    call Energy( this, nc, np, nu, EPotNew )
    ! Apply acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, MPI_RK, MPI_SUM, Communicator, ierror )

    else
          EPotDelta = EPotOld - EPotNew
    endif
#else
     EPotDelta = EPotOld - EPotNew
#endif

     ! Acceptance criterion
    if( exp(( real (this%NDF, RK) / 2._RK  - 1._RK) * log((this%RefEnthalpy*this%NPart - this%Epot+EpotDelta - this%RefPressure * this%Volume0) &
&       / (this%RefEnthalpy*this%NPart - this%Epot - this%RefPressure * this%Volume0))) > rnd( 0._RK, 1._RK ) ) then
!print*, 'MOVE', real (this%NDF, RK), this%RefEnthalpy, this%Epot, EpotDelta, this%RefPressure, this%Volume0

     ! Accept move
      this%Temperature = 2._RK * (this%RefEnthalpy*this%NPart - this%Epot+EpotDelta - this%RefPressure * this%Volume0) / real (this%NDF, RK)
      pc%NMoveSuccesses = pc%NMoveSuccesses + 1
      call UpdateEnergy( this, nc, np, nu )
    else

      ! Reject move
      if (LongRange .eq. Ewald) then
          this%UFourier = EFourier

          DO i=1,pc%Molecule%NCharge
            this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
            this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
            this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
          END DO

          pc%P0(np, :, nu) = r(:)
          pc%Pm0(np, :) = rm(:)
          call Unit2Atom1( pc, np, nu )
          call EwaldFourierEnergy(this,nc,np)

#if SPME > 0
      else if (LongRange .eq. PME) then
          this%UFourier = EFourier
          this%EVirial  = EVirial
          call chargegrid_min  (this, nc, np)
          pc%P0(np, :, nu) = r(:)
          pc%Pm0(np, :) = rm(:)
          call Unit2Atom1( pc, np, nu )
          call chargegrid_plus (this, nc, np)
#endif
      else
          pc%P0(np, :, nu) = r(:)
          pc%Pm0(np, :) = rm(:)
          call Unit2Atom1( pc, np, nu )
      end if
    end if

  end subroutine TEnsemble_Move_NPH

!==============================================================!
!  Subroutine TEnsemble_Rotate_NPH                             !
!==============================================================!

  subroutine TEnsemble_Rotate_NPH( this, nc, np, nu )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments03
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np, nu

    ! Declare local variables
    real(RK)                  :: q(4), dq(3)
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier, EVirial
    type(TComponent), pointer :: pc
    integer                   :: i
    real(RK)                  :: EPotDelta
    logical                   :: accepted

    ! Assign local variables
    pc => this%Component(nc)

    ! Update number of rotation attempts
    pc%NRotateAttempts = pc%NRotateAttempts + 1

    ! Save current particle orientation and energy
    q(:) = pc%Q0(np, :, nu)
    EPotOld = GetEnergy( this, nc, np, nu )

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO

#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      call chargegrid_min  (this, nc, np)
#endif
    end if

    ! Generate a trial rotation
    do i = 1, 3
      dq(i) = rnd( -pc%DispRot, pc%DispRot )
    end do
    pc%Q0(np, 1, nu) = q(1) - dq(1) * q(2) - dq(2) * q(3) - dq(3) * q(4)
    pc%Q0(np, 2, nu) = q(2) + dq(1) * q(1) - dq(2) * q(4) + dq(3) * q(3)
    pc%Q0(np, 3, nu) = q(3) + dq(1) * q(4) + dq(2) * q(1) - dq(3) * q(2)
    pc%Q0(np, 4, nu) = q(4) - dq(1) * q(3) + dq(2) * q(2) + dq(3) * q(1)


    ! Convert molecular coordinates to atom positions
    call Unit2Atom1( pc, np, nu )

#if SPME > 0
    if (LongRange .eq. PME) then
      call chargegrid_plus (this, nc, np)
    end if
#endif

    ! Calculate particle energy with trial orientation
    call Energy( this, nc, np, nu, EPotNew )

    ! Apply acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, MPI_RK, MPI_SUM, Communicator, ierror )
    else
      EPotDelta = EPotOld - EPotNew
    endif

#else
    EPotDelta = EPotOld - EPotNew
#endif

     ! Acceptance criterion
    if( exp(( real (this%NDF, RK) / 2._RK - 1._RK) * log((this%RefEnthalpy*this%NPart - this%Epot+EpotDelta - this%RefPressure * this%Volume0) &
&       / (this%RefEnthalpy*this%NPart - this%Epot - this%RefPressure * this%Volume0))) > rnd( 0._RK, 1._RK ) ) then
!print*, 'ROTATE', real (this%NDF, RK), this%RefEnthalpy, this%Epot, EpotDelta, this%RefPressure, this%Volume0
     ! Accept rotation
      this%Temperature = 2._RK * (this%RefEnthalpy*this%NPart - this%Epot+EpotDelta - this%RefPressure * this%Volume0) / real (this%NDF, RK)
      pc%NRotateSuccesses = pc%NRotateSuccesses + 1
      call UpdateEnergy( this, nc, np, nu )
    else

      ! Reject move
      if (LongRange .eq. Ewald) then
        this%UFourier = EFourier
        DO i=1,pc%Molecule%NCharge
          this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
          this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
          this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
        END DO
        pc%Q0(np, :, nu) = q(:)
        call Unit2Atom1( pc, np, nu )
        call EwaldFourierEnergy(this,nc,np)

#if SPME > 0
      else if (LongRange .eq. PME) then
        this%UFourier = EFourier
        this%EVirial  = EVirial
        call chargegrid_min  (this, nc, np)
        pc%Q0(np, :, nu) = q(:)
        call Unit2Atom1( pc, np, nu )
        call chargegrid_plus (this, nc, np)
#endif

      else
        pc%Q0(np, :, nu) = q(:)
        call Unit2Atom1( pc, np, nu )
      end if

    end if

  end subroutine TEnsemble_Rotate_NPH

!==============================================================!
!  Subroutine TEnsemble_MoveBiased                             !
!==============================================================!

  subroutine TEnsemble_MoveBiased( this, nc, np, nu, ncf, npf )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np, nu, ncf, npf

    ! Declare local variables
    real(RK)                  :: r(3), rm(3), dr(3), f1, f2
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier
#if SPME > 0
    real(RK)                  :: EVirial
#endif
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
    r(:) = pc%P0(np, :, nu)
    rm(:) = pc%Pm0(np, :)
    EPotOld = GetEnergy( this, nc, np, nu )

    ! Apply distance criterion
    dr(:) = r(:) - pcf%P0(npf, :, nu)
    dr(:) = ( dr(:) - anint( dr(:) ) ) * this%BoxLength
    f1 = 1._RK / ( dr(1)**2 + dr(2)**2 + dr(3)**2 )**2
    if( rnd(0._RK, 1._RK) < (1._RK - f1) ) return

    ! Update number of move attempts
    pc%NMoveBiasedAttempts = pc%NMoveBiasedAttempts + 1

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO

#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      this%qgrida_old = this%qgrida
      call chargegrid_min  (this, nc, np)
#endif
    end if

    ! Generate a trial displacement
    do i = 1, 3
      pc%P0(np, i,nu) = pc%P0(np, i,nu) + rnd( -pc%DispTran, pc%DispTran )
    end do

    ! Calculate new COM
    call Unit2Mol( pc, np )

    ! Apply periodic boundary conditions
    pc%P0(np, :, nu) = pc%P0(np, :, nu) - anint( pc%Pm0(np, :) )
    pc%Pm0(np, :)    = pc%Pm0(np, :) - anint( pc%Pm0(np, :) )

    ! Apply direction criterion
    dr(:) = pc%P0(np, :, nu) - pcf%P0(npf, :, nu)
    dr(:) = ( dr(:) - anint( dr(:) ) ) * this%BoxLength
    f2 = 1._RK / ( dr(1)**2 + dr(2)**2 + dr(3)**2 )**2
    if( rnd(0._RK, 1._RK) < (1._RK - f2/f1) ) then
      pc%P0(np, :, nu) = r(:)
      pc%Pm0(np, :) = rm(:)
      return
    end if

    ! Convert unit coordinates to atom positions
    call Unit2Atom1( pc, np, nu )

#if SPME > 0
    ! Save Energies, Virials for faster Rejection
    if (LongRange .eq. PME) then
      call chargegrid_plus (this, nc, np)
    end if
#endif

    ! Calculate particle energy at trial position
    call Energy( this, nc, np, nu, EPotNew )

    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, MPI_RK, MPI_SUM, Communicator, ierror )
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
      call UpdateEnergy( this, nc, np, nu )

    else

      ! Reject move
      if (LongRange .eq. Ewald) then
          this%UFourier = EFourier
          DO i=1,pc%Molecule%NCharge
            this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
            this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
            this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
          END DO
          pc%P0(np, :, nu)  = r(:)
          pc%Pm0(np, :) = rm(:)
          call Unit2Atom1( pc, np, nu )
          call EwaldFourierEnergy(this,nc,np)

#if SPME > 0
      else if (LongRange .eq. PME) then
        pc%P0(np, :, nu)  = r(:)
        pc%Pm0(np, :) = rm(:)
        call Unit2Atom1( pc, np, nu )
        this%UFourier = EFourier
        this%EVirial  = EVirial
        this%qgrida   = this%qgrida_old
#endif
      else
        pc%P0(np, :, nu)  = r(:)
        pc%Pm0(np, :) = rm(:)
        call Unit2Atom1( pc, np, nu )
      end if

    end if

  end subroutine TEnsemble_MoveBiased



!==============================================================!
!  Subroutine TEnsemble_RotateBiased                           !
!==============================================================!

  subroutine TEnsemble_RotateBiased( this, nc, np, nu, ncf, npf )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments03
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np, nu, ncf, npf

    ! Declare local variables
    real(RK)                  :: dr(3), q(4), dq(3), f1
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier
#if SPME > 0
    real(RK)                  :: EVirial
#endif
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
    q(:) = pc%Q0(np, :, nu)
    EPotOld = GetEnergy( this, nc, np, nu )

    ! Apply distance criterion
    dr(:) = pc%P0(np, :, nu) - pcf%P0(npf, :, nu)
    dr(:) = ( dr(:) - anint( dr(:) ) ) * this%BoxLength
    f1 = 1._RK / ( dr(1)**2 + dr(2)**2 + dr(3)**2 )**2
    if( rnd(0._RK, 1._RK) < (1._RK - f1) ) return

    ! Update number of rotation attempts
    pc%NRotateBiasedAttempts = pc%NRotateBiasedAttempts + 1

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO

#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      this%qgrida_old = this%qgrida
      call chargegrid_min  (this, nc, np)
#endif
    end if

    ! Generate a trial rotation
    do i = 1, 3
      dq(i) = rnd( -pc%DispRot, pc%DispRot )
    end do
    ! rotate unit 
    pc%Q0(np, 1,nu) = q(1) - dq(1) * q(2) - dq(2) * q(3) - dq(3) * q(4)
    pc%Q0(np, 2,nu) = q(2) + dq(1) * q(1) - dq(3) * q(3) + dq(2) * q(4)
    pc%Q0(np, 3,nu) = q(3) + dq(2) * q(1) + dq(3) * q(2) - dq(1) * q(4)
    pc%Q0(np, 4,nu) = q(4) + dq(3) * q(1) - dq(2) * q(2) + dq(1) * q(3)

    ! Convert unit coordinates to atom positions
    call Unit2Atom1( pc, np, nu )

#if SPME > 0
    ! Save Energies, Virials for faster Rejection
    if (LongRange .eq. PME) then
      call chargegrid_plus (this, nc, np)
    end if
#endif

    ! Calculate particle energy with trial orientation
    call Energy( this, nc, np, nu, EPotNew )

    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, MPI_RK, MPI_SUM, Communicator, ierror )
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
      call UpdateEnergy( this, nc, np, nu )

    else

      ! Reject move
      if (LongRange .eq. Ewald) then
        this%UFourier = EFourier
        DO i=1,pc%Molecule%NCharge
          this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
          this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
          this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
        END DO
        pc%Q0(np, :,nu) = q(:)
        call Unit2Atom1( pc, np, nu )
        call EwaldFourierEnergy(this,nc,np)

#if SPME > 0
      else if (LongRange .eq. PME) then
        pc%Q0(np,:,nu) = q(:)
        call Unit2Atom1( pc, np, nu )
        this%UFourier = EFourier
        this%EVirial  = EVirial
        this%qgrida   = this%qgrida_old

#endif
      else
        pc%Q0(np,:,nu) = q(:)
        call Unit2Atom1( pc, np, nu )
      end if

    end if

  end subroutine TEnsemble_RotateBiased



!==============================================================!
!  Subroutine TEnsemble_PartnersBiased                         !
!==============================================================!

  subroutine TEnsemble_PartnersBiased ( this, ncf, npf )
  
    implicit none

    ! Declare arguments
    type(TEnsemble)         :: this
    integer, intent(in)     :: ncf, npf
    
    real(RK)            :: BoxLength
    real(RK)            :: dxf, dyf, dzf
    real(RK)            :: dx, dy, dz, dr2
    integer             :: NGradIns
    integer             :: counter, counter1
    integer             :: i, j
    type(TComponent),pointer :: pc, pcf
    
    BoxLength = this%BoxLength
    counter   = 0
    counter1  = 0
    NGradIns  = 0

    pcf => this%Component(ncf)
    dxf = pcf%Pm0(npf,1)
    dyf = pcf%Pm0(npf,2)
    dzf = pcf%Pm0(npf,3)
    
    do i=1, this%NComponents
      pc => this%Component(i)
      counter1=0
      do j=1,pc%NPart
        if (i .eq. ncf .and. npf .eq. j) cycle
        dx = pc%Pm0(j,1) - dxf
        dy = pc%Pm0(j,2) - dyf
        dz = pc%Pm0(j,3) - dzf

        dx  = (dx-anint(dx))*BoxLength
        dy  = (dy-anint(dy))*BoxLength
        dz  = (dz-anint(dz))*BoxLength
        dr2 = dx**2 + dy**2 + dz**2

        if ( dr2 .le. 25._RK ) then
          counter = counter + 1
          counter1= counter1+ 1
          this%BiasedPartners(counter) = j
        end if
      end do
      pc%BiasedPartners = counter1*pc%Molecule%NDF
      pc%BiasedPartnersNum = counter1
      NGradIns = NGradIns + counter1*pc%Molecule%NDF
    end do
    
    this%NGradIns = NGradIns
  
  end subroutine TEnsemble_PartnersBiased


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

    ! Declare local variables
    type(TComponent), pointer :: pc, pcf, pcfnew
    integer                   :: oldstate, newstate
    integer                   :: ncfnew, npfnew
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EPotDeltaAll
    real(RK)                  :: EFourier, EVirial

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
    call DuplicateParticle( pcfnew, pcf, npf )
    call RemoveParticle( pcf, npf )
    npfnew = pcfnew%NPart

! Save states for the Ewald Summation and/or derivates
    if (LongRange .eq. Ewald) then     ! Ewald Summation
       ! Save the initial state
       EFourier = this%UFourier
       EPotOld = EPotOld  + this%USelbstTerm + this%UIntra

       if ( this%OptPressure ) then
         EVirial  = this%EVirial
       end if
!  Sufficient, since no call to Mol2Atom1 yet

       do i=1,pcf%Molecule%NCharge
         this%rold(i,1) = pcf%Molecule%SiteCharge(i)%RX(npf)
         this%rold(i,2) = pcf%Molecule%SiteCharge(i)%RY(npf)
         this%rold(i,3) = pcf%Molecule%SiteCharge(i)%RZ(npf)
       end do

       ! Calculate new energies
       call EwaldSelfTerm_Energy(this)

       ! Convert unit coordinates to atom positions
       call Unit2Atom1( pcfnew, npfnew )

       ! Calculate particle energy at new fluctuating state
       call Energy( this, ncfnew, npfnew, ncf, npf, EPotNew )

       ! Acceptance Criteria
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
       call MPI_Allreduce( EPotOld - EPotNew, EPotDeltaAll, 1, MPI_RK, MPI_SUM, Communicator, ierror )
       EPotDeltaAll = EPotDeltaAll + this%Density * ( pcf%EPotTestCorrLJ - pcfnew%EPotTestCorrLJ )
    else
       EPotOld = EPotOld + this%Density * ( pcf%EPotTestCorrLJ - pcfnew%EPotTestCorrLJ )
       EPotDeltaAll = EPotOld - EPotNew
    end if

#else
       EPotOld = EPotOld + this%Density * ( pcf%EPotTestCorrLJ - pcfnew%EPotTestCorrLJ )
       EPotDeltaAll = EPotOld - EPotNew
#endif

       if( rnd( 0._RK, 1._RK ) < pc%WF(newstate) / pc%WF(oldstate) * exp( ( EPotDeltaAll ) / this%Temperature ) ) then
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
         call DuplicateParticle( pcf, pcfnew, npf )
         call RemoveParticle( pcfnew, npfnew )
         call Unit2Atom1( pcf, npf )
         call EwaldSelfTerm_Energy(this)
         do i=1,pcfnew%Molecule%NCharge
           this%rold(i,1) = pcfnew%Molecule%SiteCharge(i)%RX(npfnew)
           this%rold(i,2) = pcfnew%Molecule%SiteCharge(i)%RY(npfnew)
           this%rold(i,3) = pcfnew%Molecule%SiteCharge(i)%RZ(npfnew)
         end do
         call Energy( this, ncf, npf, ncfnew, npfnew, EPotNew )

       end if       ! Acceptance Criteria

#if SPME > 0
! ----------------------------------------------------------------
    else if (LongRange .eq. PME) then ! PME 
      EFourier = this%UFourier
      EVirial  = this%EVirial
      call PMESetup(this)
      write (*,*) 'Gradual Insertion does not yet work with PME'
      STOP
#endif
! ----------------------------------------------------------------
    else   ! REACTION FIELD
       ! Convert unit coordinates to atom positions
       call Unit2Atom1( pcfnew, npfnew )

       ! Calculate particle energy at new fluctuating state
       call Energy( this, ncfnew, npfnew, EPotNew )

    ! Apply acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
       call MPI_Allreduce( EPotOld - EPotNew, EPotDeltaAll, 1, MPI_RK, MPI_SUM, Communicator, ierror )
       EPotDeltaAll = EPotDeltaAll + this%Density * ( pcf%EPotTestCorrLJ - pcfnew%EPotTestCorrLJ ) &
&        + pcf%EPotTestCorrRF - pcfnew%EPotTestCorrRF

    else
       EPotOld = EPotOld + this%Density * ( pcf%EPotTestCorrLJ - pcfnew%EPotTestCorrLJ ) &
&        + pcf%EPotTestCorrRF - pcfnew%EPotTestCorrRF

       EPotDeltaAll = EPotOld - EPotNew
    end if
#else
       EPotOld = EPotOld + this%Density * ( pcf%EPotTestCorrLJ - pcfnew%EPotTestCorrLJ ) &
&        + pcf%EPotTestCorrRF - pcfnew%EPotTestCorrRF

       EPotDeltaAll = EPotOld - EPotNew
#endif
       if( rnd( 0._RK, 1._RK ) < pc%WF(newstate) / pc%WF(oldstate) * exp( ( EPotDeltaAll ) / this%Temperature ) ) then

         ! Accept
         pc%NFluctState = newstate
         ncf = ncfnew
         npf = npfnew
         call UpdateEnergy( this, ncf, npf )

         if( newstate > oldstate ) then
           pc%NFluctUpSuccesses(newstate) = pc%NFluctUpSuccesses(newstate)+1
         else
           pc%NFluctDownSuccesses(oldstate) = pc%NFluctDownSuccesses(oldstate)+1
         end if

        if ( newstate .eq. 0 ) call BiasedPartners ( this, ncf, npf )

       else

         ! Reject
         call DuplicateParticle( pcf, pcfnew, npf )
         call RemoveParticle( pcfnew, npfnew )
         call Unit2Atom1( pcf, npf )

       end if

    end if      ! LongRange - Correction

  end subroutine TEnsemble_ChangeFluct

!==============================================================!
!  Subroutine TEnsemble_ScaleInteractionThermoInt              !
!==============================================================!

subroutine TEnsemble_ScaleInteractionThermoInt( this, nt , factor)

   implicit none

    ! Declare arguments
    type(TEnsemble)        :: this
    integer, intent(in)    :: nt
    real(RK), intent(in)   :: factor

    ! Declare local variables
    integer                 :: i

    do i = 1, this%NComponents
      if (nt == i) cycle

      if( associated(this%Interaction(nt, i)%PotLJ126LJ126)) then
        this%Interaction(nt, i)%PotLJ126LJ126(:, :)%Epsilon      = this%Interaction(nt, i)%PotLJ126LJ126(:, :)%Epsilon * Factor
        this%Interaction(nt, i)%PotLJ126LJ126(:, :)%Epsilon4     = this%Interaction(nt, i)%PotLJ126LJ126(:, :)%Epsilon4 * Factor
        this%Interaction(nt, i)%PotLJ126LJ126(:, :)%Epsilon48    = this%Interaction(nt, i)%PotLJ126LJ126(:, :)%Epsilon48 * Factor
        this%Interaction(i, nt)%PotLJ126LJ126(:, :)%Epsilon      = this%Interaction(i, nt)%PotLJ126LJ126(:, :)%Epsilon * Factor
        this%Interaction(i, nt)%PotLJ126LJ126(:, :)%Epsilon4     = this%Interaction(i, nt)%PotLJ126LJ126(:, :)%Epsilon4 * Factor
        this%Interaction(i, nt)%PotLJ126LJ126(:, :)%Epsilon48    = this%Interaction(i, nt)%PotLJ126LJ126(:, :)%Epsilon48 * Factor
      endif
      if( associated(this%Interaction(nt, i)%PotChargeCharge)) then
        this%Interaction(nt, i)%PotChargeCharge(:, :)%Epsilon    = this%Interaction(nt, i)%PotChargeCharge(:, :)%Epsilon * Factor
        this%Interaction(i, nt)%PotChargeCharge(:, :)%Epsilon    = this%Interaction(i, nt)%PotChargeCharge(:, :)%Epsilon * Factor
        !do k=1,this%Interaction(nt,i)%N1Charge
        !  Shield1 = this%Component(nt)%Molecule%SiteCharge(k)%shield
        !  do l=1,this%Interaction(nt,i)%N2Charge
        !    Shield2 = this%Component(i)%Molecule%SiteCharge(l)%shield
        !    this%Interaction(nt, i)%PotChargeCharge(k, l)%RShieldSquared = .25_RK * ( Shield1 * Factor + Shield2 )**2
        !    this%Interaction(i, nt)%PotChargeCharge(l, k)%RShieldSquared = .25_RK * ( Shield2 + Shield1 * Factor )**2
        !  enddo
        !enddo
      endif
      if( associated(this%Interaction(nt, i)%PotChargeDipole)) then
        this%Interaction(nt, i)%PotChargeDipole(:, :)%Epsilon    = this%Interaction(nt, i)%PotChargeDipole(:, :)%Epsilon * Factor
        this%Interaction(i, nt)%PotDipoleCharge(:, :)%Epsilon    = this%Interaction(i, nt)%PotDipoleCharge(:, :)%Epsilon * Factor
        !do k=1,this%Interaction(nt,i)%N1Charge
        !  Shield1 = this%Component(nt)%Molecule%SiteCharge(k)%shield
        !  do l=1,this%Interaction(nt,i)%N2Dipole
        !    Shield2 = this%Component(i)%Molecule%SiteDipole(l)%shield
        !    this%Interaction(nt, i)%PotChargeDipole(k, l)%RShieldSquared = .25_RK * ( Shield1 * Factor + Shield2 )**2
        !  enddo
        !enddo
        !do k=1,this%Interaction(i,nt)%N1Dipole
        !  Shield2 = this%Component(i)%Molecule%SiteDipole(k)%shield
        !  do l=1,this%Interaction(i,nt)%N2Charge
        !    Shield1 = this%Component(nt)%Molecule%SiteCharge(l)%shield
        !    this%Interaction(i, nt)%PotDipoleCharge(k, l)%RShieldSquared = .25_RK * ( Shield2  + Shield1 * Factor)**2
        !  enddo
        !enddo
      endif
      if( associated(this%Interaction(nt, i)%PotChargeQuadrupole)) then       
        this%Interaction(nt, i)%PotChargeQuadrupole(:, :)%Epsilon    = this%Interaction(nt, i)%PotChargeQuadrupole(:, :)%Epsilon * Factor
        this%Interaction(i, nt)%PotQuadrupoleCharge(:, :)%Epsilon    = this%Interaction(i, nt)%PotQuadrupoleCharge(:, :)%Epsilon * Factor
        !do k=1,this%Interaction(nt,i)%N1Charge
        !  Shield1 = this%Component(nt)%Molecule%SiteCharge(k)%shield
        !  do l=1,this%Interaction(nt,i)%N2Quadrupole
        !    Shield2 = this%Component(i)%Molecule%SiteQuadrupole(l)%shield
        !    this%Interaction(nt, i)%PotChargeQuadrupole(k, l)%RShieldSquared = .25_RK * ( Shield1 * Factor + Shield2 )**2
        !  enddo
        !enddo
        !do k=1,this%Interaction(i,nt)%N1Quadrupole
        !  Shield2 = this%Component(i)%Molecule%SiteQuadrupole(k)%shield
        !  do l=1,this%Interaction(i,nt)%N2Charge
        !    Shield1 = this%Component(nt)%Molecule%SiteCharge(l)%shield
        !    this%Interaction(i, nt)%PotQuadrupoleCharge(k, l)%RShieldSquared = .25_RK * ( Shield2  + Shield1 * Factor)**2
        !  enddo
        !enddo
      endif
      if( associated(this%Interaction(nt, i)%PotDipoleCharge)) then
        this%Interaction(nt, i)%PotDipoleCharge(:, :)%Epsilon    = this%Interaction(nt, i)%PotDipoleCharge(:, :)%Epsilon * Factor
        this%Interaction(i, nt)%PotChargeDipole(:, :)%Epsilon    = this%Interaction(i, nt)%PotChargeDipole(:, :)%Epsilon * Factor
        !do k=1,this%Interaction(nt,i)%N1Dipole
        !  Shield1 = this%Component(nt)%Molecule%SiteDipole(k)%shield
        !  do l=1,this%Interaction(nt,i)%N2Charge
        !    Shield2 = this%Component(i)%Molecule%SiteCharge(l)%shield
        !    this%Interaction(nt, i)%PotDipoleCharge(k, l)%RShieldSquared = .25_RK * ( Shield1 * Factor + Shield2 )**2
        !  enddo
        !enddo
        !do k=1,this%Interaction(i,nt)%N1Charge
        !  Shield2 = this%Component(i)%Molecule%SiteCharge(k)%shield
        !  do l=1,this%Interaction(i,nt)%N2Dipole
        !    Shield1 = this%Component(nt)%Molecule%SiteDipole(l)%shield
        !    this%Interaction(i, nt)%PotChargeDipole(k, l)%RShieldSquared = .25_RK * ( Shield2  + Shield1 * Factor)**2
        !  enddo
        !enddo
      endif
      if( associated(this%Interaction(nt, i)%PotDipoleDipole)) then
        this%Interaction(nt, i)%PotDipoleDipole(:, :)%Epsilon    = this%Interaction(nt, i)%PotDipoleDipole(:, :)%Epsilon * Factor
        this%Interaction(i, nt)%PotDipoleDipole(:, :)%Epsilon    = this%Interaction(i, nt)%PotDipoleDipole(:, :)%Epsilon * Factor
        !do k=1,this%Interaction(nt,i)%N1Dipole
        !  Shield1 = this%Component(nt)%Molecule%SiteDipole(k)%shield
        !  do l=1,this%Interaction(nt,i)%N2Dipole
        !    Shield2 = this%Component(i)%Molecule%SiteDipole(l)%shield
        !    this%Interaction(nt, i)%PotDipoleDipole(k, l)%RShieldSquared = .25_RK * ( Shield1 * Factor + Shield2 )**2
        !    this%Interaction(i, nt)%PotDipoleDipole(l, k)%RShieldSquared = .25_RK * ( Shield2 + Shield1 * Factor )**2
        !  enddo
        !enddo
      endif
      if( associated(this%Interaction(nt, i)%PotDipoleQuadrupole)) then
        this%Interaction(nt, i)%PotDipoleQuadrupole(:, :)%Epsilon    = this%Interaction(nt, i)%PotDipoleQuadrupole(:, :)%Epsilon * Factor
        this%Interaction(i, nt)%PotQuadrupoleDipole(:, :)%Epsilon    = this%Interaction(i, nt)%PotQuadrupoleDipole(:, :)%Epsilon * Factor
        !do k=1,this%Interaction(nt,i)%N1Dipole
        !  Shield1 = this%Component(nt)%Molecule%SiteDipole(k)%shield
        !  do l=1,this%Interaction(nt,i)%N2Quadrupole
        !    Shield2 = this%Component(i)%Molecule%SiteQuadrupole(l)%shield
        !    this%Interaction(nt, i)%PotDipoleQuadrupole(k, l)%RShieldSquared = .25_RK * ( Shield1 * Factor + Shield2 )**2
        !  enddo
        !enddo
        !do k=1,this%Interaction(i,nt)%N1Quadrupole
        !  Shield2 = this%Component(i)%Molecule%SiteQuadrupole(k)%shield
        !  do l=1,this%Interaction(i,nt)%N2Dipole
        !    Shield1 = this%Component(nt)%Molecule%SiteDipole(l)%shield
        !    this%Interaction(i, nt)%PotQuadrupoleDipole(k, l)%RShieldSquared = .25_RK * ( Shield2  + Shield1 * Factor)**2
        !  enddo
        !enddo
      endif
      if( associated(this%Interaction(nt, i)%PotQuadrupoleCharge)) then
        this%Interaction(nt, i)%PotQuadrupoleCharge(:, :)%Epsilon    = this%Interaction(nt, i)%PotQuadrupoleCharge(:, :)%Epsilon * Factor
        this%Interaction(i, nt)%PotChargeQuadrupole(:, :)%Epsilon    = this%Interaction(i, nt)%PotChargeQuadrupole(:, :)%Epsilon * Factor
        !do k=1,this%Interaction(nt,i)%N1Quadrupole
        !  Shield1 = this%Component(nt)%Molecule%SiteQuadrupole(k)%shield
        !  do l=1,this%Interaction(nt,i)%N2Charge
        !    Shield2 = this%Component(i)%Molecule%SiteCharge(l)%shield
        !    this%Interaction(nt, i)%PotQuadrupoleCharge(k, l)%RShieldSquared = .25_RK * ( Shield1 * Factor + Shield2 )**2
        !  enddo
        !enddo
        !do k=1,this%Interaction(i,nt)%N1Charge
        !  Shield2 = this%Component(i)%Molecule%SiteCharge(k)%shield
        !  do l=1,this%Interaction(i,nt)%N2Quadrupole
        !    Shield1 = this%Component(nt)%Molecule%SiteQuadrupole(l)%shield
        !    this%Interaction(i, nt)%PotChargeQuadrupole(k, l)%RShieldSquared = .25_RK * ( Shield2  + Shield1 * Factor)**2
        !  enddo
        !enddo
      endif
      if( associated(this%Interaction(nt, i)%PotQuadrupoleDipole)) then
        this%Interaction(nt, i)%PotQuadrupoleDipole(:, :)%Epsilon    = this%Interaction(nt, i)%PotQuadrupoleDipole(:, :)%Epsilon * Factor
        this%Interaction(i, nt)%PotDipoleQuadrupole(:, :)%Epsilon    = this%Interaction(i, nt)%PotDipoleQuadrupole(:, :)%Epsilon * Factor
        !do k=1,this%Interaction(nt,i)%N1Quadrupole
        !  Shield1 = this%Component(nt)%Molecule%SiteQuadrupole(k)%shield
        !  do l=1,this%Interaction(nt,i)%N2Dipole
        !    Shield2 = this%Component(i)%Molecule%SiteDipole(l)%shield
        !    this%Interaction(nt, i)%PotQuadrupoleDipole(k, l)%RShieldSquared = .25_RK * ( Shield1 * Factor + Shield2 )**2
        !  enddo
        !enddo
        !do k=1,this%Interaction(i,nt)%N1Dipole
        !  Shield2 = this%Component(i)%Molecule%SiteDipole(k)%shield
        !  do l=1,this%Interaction(i,nt)%N2Quadrupole
        !    Shield1 = this%Component(nt)%Molecule%SiteQuadrupole(l)%shield
        !    this%Interaction(i, nt)%PotDipoleQuadrupole(k, l)%RShieldSquared = .25_RK * ( Shield2  + Shield1 * Factor)**2
        !  enddo
        !enddo
      endif
      if( associated(this%Interaction(nt, i)%PotQuadrupoleQuadrupole)) then
        this%Interaction(nt, i)%PotQuadrupoleQuadrupole(:, :)%Epsilon= this%Interaction(nt, i)%PotQuadrupoleQuadrupole(:, :)%Epsilon * Factor
        this%Interaction(i, nt)%PotQuadrupoleQuadrupole(:, :)%Epsilon= this%Interaction(i, nt)%PotQuadrupoleQuadrupole(:, :)%Epsilon * Factor
        !do k=1,this%Interaction(nt,i)%N1Quadrupole
        !  Shield1 = this%Component(nt)%Molecule%SiteQuadrupole(k)%shield
        !  do l=1,this%Interaction(nt,i)%N2Quadrupole
        !    Shield2 = this%Component(i)%Molecule%SiteQuadrupole(l)%shield
        !    this%Interaction(nt, i)%PotQuadrupoleQuadrupole(k, l)%RShieldSquared = .25_RK * ( Shield1 * Factor + Shield2 )**2
        !    this%Interaction(i, nt)%PotQuadrupoleQuadrupole(l, k)%RShieldSquared = .25_RK * ( Shield2 + Shield1 * Factor )**2
        !  enddo
        !enddo
      endif
    end do
    if( associated(this%Component(nt)%MueX)) then  ! if MueX then also MueY and Z
      do i=1,this%Component(nt)%Molecule%NUnit
        this%Component(nt)%Molecule%Unit(i)%Mue(:) = this%Component(nt)%Molecule%Unit(i)%Mue(:) * Factor
      end do
      !this%Component(nt)%Molecule%MueY(:) = this%Component(nt)%Molecule%MueY(:) * Factor
      !this%Component(nt)%Molecule%MueZ(:) = this%Component(nt)%Molecule%MueZ(:) * Factor
    endif

end subroutine TEnsemble_ScaleInteractionThermoInt


!==============================================================!
!  Subroutine TEnsemble_ChangeLambda                           !
!==============================================================!

  subroutine TEnsemble_ChangeLambda( this, nt , nc)

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)        :: this
    integer, intent(in)    :: nt
    integer, intent(in)    :: nc

    ! Declare local variables
    type(TComponent), pointer  :: pc, pt
    type(TInteraction), pointer:: plj
    integer                    :: i, j, k, l, currentbin
    real(RK)                   :: Shield1, Shield2
    real(RK)                   :: LambdaNew, Factor, FactorOld, ChempotDelta
    real(RK)                   :: EPotOld, EPotNew
    real(RK)                   :: EPotDeltaAll, Scale
    real(RK)                   :: EFourier, EVirial

    ! Assign local variables
    pt => this%Component(nt)
    pc => this%Component(nc)

    ! Get old energy of fluctuating particle
    if (SimulationType .ne. MolecularDynamics) then
      if ( Step .gt. pc%changeLaPart .and. nint(pt%Lambda/pc%deltaLa) .ge. pt%NBins) then
        pc%changeLaPart = pc%changeLaPart + pc%changeLaPart
        call ChangeFluct( this, nt, nc )
      end if
      EPotOld = (this%Density * pc%EPotTestCorrLJ + pc%EPotTestCorrRF)*pt%Lambda**pc%LambdaExponent
      EPotOld = EPotOld + GetEnergy( this, nt, 1 )

      ! Save states for the Ewald Summation and/or derivates
      if (LongRange .eq. Ewald) then     ! Ewald Summation
        ! Save the initial state
        EFourier = this%UFourier
        EPotOld = EPotOld  + this%USelbstTerm + this%UIntra

        if ( this%OptPressure ) then
          EVirial  = this%EVirial
        end if

        call EwaldSelfTerm_Energy(this)
      end if

      ! Change state of lambda
      LambdaNew=pt%Lambda+2.0_RK*pc%LaStepMax*(rnd(0.0_RK,1.0_RK)-0.5_RK)

      if (LambdaNew>=pc%LaMin .and. LambdaNew<=pc%LaMax) then 
        
        currentbin=int((LambdaNew-pc%LaMin)/pc%deltaLa)
        ChempotDelta=ChempotDelta+pc%BinsIntdEndLa(currentbin)
        ! Calculate energy of fluctuating particle
        Factor = (LambdaNew/pt%Lambda)**pc%LambdaExponent
        EPotNew=Factor*EPotOld

        ! Acceptance Criteria
        EPotDeltaAll = EPotOld - EPotNew
        if( rnd( 0._RK, 1._RK ) < exp( ( EPotDeltaAll + ChempotDelta) / this%Temperature ) ) then
          ! Accept
          ! Apply scaling factors
          call ScaleInteractionThermoInt(this, nt, Factor)
          !call Unit2Atom( this )
          call Unit2Atom1( pt, 1 )
          call Energy( this, nt, 1, EPotNew )
          call UpdateEnergy( this, nt, 1 )
          pt%Lambda=LambdaNew
        else
          ! Reject
          if (LongRange == Ewald) then
            call EwaldSelfTerm_Energy(this)
            call Energy( this, nt, 1, EPotNew )
          end if
        end if       ! Acceptance Criteria

      end if !Lamba change

    else ! MolecularDynamics

      if ( Step .gt. pc%changeLaPart .and. pt%Lambda+pc%LaStepMax > pc%LaMax) then
        pc%changeLaPart = pc%changeLaPart + pc%changeLaPart
        call ChangeFluct( this, nt, nc )
      end if
      if (RootProc) then
        LambdaNew=pt%Lambda+pc%LaStepMax
        ! should be 1/10 of MC-stepwidth for equl distribution (estimation by Gabor and Michael)
        if (LambdaNew<pc%LaMin .or. LambdaNew>pc%LaMax) then
          pc%LaStepMax = -pc%LaStepMax
          LambdaNew = pt%Lambda ! +pc%LaStepMax ! Minh test
        end if
        Factor = (LambdaNew/pt%Lambda)**pc%LambdaExponent
        pt%Lambda=LambdaNew
      end if
      ! Apply scaling factors
#if MPI_VER > 0
      call MPI_Bcast( Factor, 1, MPI_RK, NRootProc, Communicator, ierror )
      call MPI_Bcast( pt%Lambda, 1, MPI_RK, NRootProc, Communicator, ierror )
#endif
      call ScaleInteractionThermoInt(this, nt, Factor)

      call Unit2Atom1( pt, 1 )

    end if

  end subroutine TEnsemble_ChangeLambda


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
    real(RK)                  :: q(3)
    real(RK)                  :: EPotIns
    type(TComponent), pointer :: pc
    integer                   :: i, np, nu, dummy, j
    logical                   :: success, barrier
    real(RK)                  :: UIntra, USelbst, EFourier, EVirial
    real(RK)                  :: E, EIntra, EBond, EAngle, EDihedral, FIns(3,this%NUnitMax), InvDensityCorr
#if MPI_VER > 0
    real(RK)                  :: EPotInsAll
#endif

    ! Assign local variables
    pc => this%Component(nc)
    nu = pc%Molecule%NUnit
    success = .true.
    barrier = .true.

    do dummy = 1, 1
    if ( (SimulationType .eq. MonteCarlo) .or. ( (SimulationType .eq. MolecularDynamics) .and. RootProc ) ) then

      success = .false.
      ! Update number of insert attempts
      this%NInsertAttempts = this%NInsertAttempts + 1

      ! Generate a random position and orientation
      do i = 1, 3
        r(i) = rnd( -.5_RK, .5_RK )
      end do
      do i = 1, 3
        q(i) = rnd( -1._RK, 1._RK )
      end do
      !  Michael Sch.: Rotation problematic with IDF, especially for MD since velocities are not changed alongside rotation
      ! Instead of just duplicating and moving a particle, a velocity spin could replace the rotation....
      if (Shake > 0) q(:) = 0._RK 

      call AddParticle( pc, r, q )
      if ( tooManyParticles ) exit
      np = pc%NPart
      this%NPart = this%NPart + 1
      this%NUnitTotal = this%NUnitTotal + nu

      ! Force criteria for acceptance in MD Simulations
      ! derived from standard deviation of the velocity distribution ...3.57 is means 3.57 times the standard deviation 
      ! 3.57_RK * sqrt((3 * PI - 8._RK )/ PI) + sqrt(8/PI) ..= 4.0
      ! Currently Forces only implemented/calculated for LJ, bond, angle and dihedral potential!!!
      !Fbarrier(:) = 4._RK * sqrt( this%Temperature * this%Component(nc)%Molecule%Unit(:)%Mass )

      if (LongRange .eq. Ewald) then           ! EWALD-SUMMATION
        UIntra   = this%UIntra
        USelbst  = this%USelbstTerm
        EFourier = this%UFourier
        if ( this%OptPressure ) then
          EVirial  = this%EVirial
        end if
        ! Energy
        call EwaldSelfTerm_Energy(this)
        call Energy ( this, nc, np, EPotIns, 1 )

#if MPI_VER > 0
        if ( Equilibration .and. CommonEqui ) then
          ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
          call MPI_Allreduce( EPotIns, EPotInsAll, 1, MPI_RK, MPI_SUM, Communicator, ierror )
          EPotInsAll = EPotInsAll + this%Density * pc%EPotTestCorrLJ + this%UIntra-UIntra + this%USelbstTerm-USelbst-EFourier

        else
          EPotInsAll = EPotIns + this%Density * pc%EPotTestCorrLJ + this%UIntra-UIntra + this%USelbstTerm-USelbst-EFourier
        endif 
   
        if( rnd( 0._RK, 1._RK ) .lt. ( exp( pc%ChemPot - EPotInsAll / this%Temperature ) * this%Volume0 / np )) then
#else
          EPotIns = EPotIns + this%Density * pc%EPotTestCorrLJ + this%UIntra-UIntra + this%USelbstTerm-USelbst-EFourier

        ! Apply acceptance criterion - SINGLE
        if( rnd( 0._RK, 1._RK ) .lt. ( exp( pc%ChemPot - EPotIns / this%Temperature ) * this%Volume0 / np )) then
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
          this%NUnitTotal = this%NUnitTotal - nu
          call EwaldFourierEnergy ( this, nc, np, -1 )
          this%USelbstTerm = USelbst
          this%UIntra  = UIntra
        end if 

#if SPME > 0
      else if (LongRange .eq. PME) then           ! PME-SUMMATION
        EVirial  = this%EVirial
        this%qgrida_old = this%qgrida
        call chargegrid_plus (this, nc, np)
        call PMESelfTermMC( this )
        write (*,*) 'Insertion and Deletion is not supported for PME!'
        STOP
#endif

      else                                         ! REACTION FIELD
        ! Calculate particle energy at trial position
        if (SimulationType .eq. MonteCarlo) then
          call Energy( this, nc, np, EPotIns )
        else
          E = 0._RK; EIntra = 0._RK; EBond = 0._RK; EAngle = 0._RK; EDihedral = 0._RK; FIns(:,:) = 0._RK;
          do j = 1, this%NComponents
            if (j > nc) then
              call MDEnergy( this%Interaction(nc,j), np, nu, FIns(:,1:nu), E, EIntra, EBond, EAngle, EDihedral, this%BoxLength, .true. )
            else
              call MDEnergy( this%Interaction(j,nc), np, nu, FIns(:,1:nu), E, EIntra, EBond, EAngle, EDihedral, this%BoxLength, .false. )
            end if
          end do
          EPotIns = E - EIntra
        end if
        InvDensityCorr = this%Volume0 / np
        if (Shake > 0) InvDensityCorr =  this%Volume0 / (this%NUnitTotal-nu-this%constrNDF/3._RK)

        ! Apply acceptance criterion
#if MPI_VER > 0
        if ( (SimulationType .eq. MonteCarlo) .and. (Equilibration .and. CommonEqui) ) then
          ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
          call MPI_Allreduce( EPotIns, EPotInsAll, 1, MPI_RK, MPI_SUM, Communicator, ierror )
          EPotInsAll = EPotInsAll + this%Density * pc%EPotTestCorrLJ + pc%EPotTestCorrRF

        else
          EPotInsAll = EPotIns + this%Density * pc%EPotTestCorrLJ + pc%EPotTestCorrRF
        endif 
    
        if( rnd( 0._RK, 1._RK ) .lt. ( exp( pc%ChemPot - EPotInsAll / this%RefTemperature ) * InvDensityCorr )) then

#else
        EPotIns = EPotIns + this%Density * pc%EPotTestCorrLJ + pc%EPotTestCorrRF
        if( rnd( 0._RK, 1._RK ) .lt. ( exp( pc%ChemPot - EPotIns / this%RefTemperature ) * InvDensityCorr )) then
#endif

          ! Accept Insertion
          this%NInsertSuccesses = this%NInsertSuccesses + 1
          if (SimulationType .ne. MonteCarlo) then
            success = .true.
            this%Density = this%NPart / this%Volume0
          else
            ! Update energy matrix
            call UpdateEnergy( this, nc, np )
            ! Update density
            this%Density = this%NPart / this%Volume0
            ! Update fractions and NDF
            call UpdateFractions( this )
            ! Update long range correction
            call CalculateCorr( this )
          end if

        else
          ! Reject Insertion
          call RemoveParticle( pc, np )
          this%NPart = this%NPart - 1
          this%NUnitTotal = this%NUnitTotal - nu
        end if 

      end if

    end if
    end do

    if (SimulationType .ne. MonteCarlo) then
#if MPI_VER > 0
      call MPI_Bcast( success, 1, MPI_LOGICAL, NRootProc, Communicator, ierror )
      if (success) then
        call MPI_Bcast( this%NPart, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NUnitTotal, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( pc%NPart, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%Density, 1, MPI_RK, NRootProc, Communicator, ierror )
        np = pc%NPart
        call Unit2Atom1( pc, np)
        pc%NPart1 = ProcRange( pc%NPart, pc%NPart0, pc%NPart2 )
        ! Update fractions and NDF
        call UpdateFractions( this )
        ! Update long range correction
        call CalculateCorr( this )
      end if
#else
      if (success) then
        ! Update fractions and NDF
        call UpdateFractions( this )
        ! Update long range correction
        call CalculateCorr( this )
      end if
#endif
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
    logical                     :: success
    integer                     :: i, n1, n2, nu, nup, k
!     real(RK)                    :: s
    real(RK)                    :: E, EIntra, EBond, DensityCorr, EAngle, EDihedral, FDel(3,this%NUnitMax)

! Ewald Parameter
    real(RK)                    :: EFourier, EPotNew
    real(RK)                    :: USelf, UIntra
    real(RK)                    :: r(3)
    real(RK)                    :: q(4)

    ! Assign local variables
    pc => this%Component(nc)
    nu = pc%Molecule%NUnit
    success = .true.

    if ( (SimulationType .eq. MonteCarlo) .or. ( (SimulationType .eq. MolecularDynamics) .and. RootProc ) ) then

      success = .false.
      ! Update number of delete attempts
      this%NDeleteAttempts = this%NDeleteAttempts + 1


      if (LongRange .eq. Ewald) then
        EFourier = this%UFourier
        !if ( this%OptPressure ) then
        !  EVirial  = this%EVirial
        !end if

        USelf    = this%USelbstTerm
        UIntra   = this%UIntra
#if MPI_VER > 0
        if ( Equilibration .and. CommonEqui ) then
          ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
          call MPI_Allreduce( GetEnergy( this, nc, np ), EPotDel, 1,MPI_RK, MPI_SUM, Communicator, ierror )
        else
          EPotDel = GetEnergy( this, nc, np )
        endif 
#else
        EPotDel = GetEnergy( this, nc, np )
#endif
        EPotDel = EPotDel + this%Density * pc%EPotTestCorrLJ + this%UIntra-UIntra + this%USelbstTerm-USelf-EFourier

        ! Apply acceptance criterion
        if( rnd( 0._RK, 1._RK ) .lt. ( exp( EPotDel / this%Temperature - pc%ChemPot ) * this%Density * pc%Fraction )) then

          ! Accept Deletion
          this%NDeleteSuccesses = this%NDeleteSuccesses + 1
          call RemoveParticle( pc, np )

          ! Copy energies and virial
          nup = (np-1)*nu
          n1 = pc%NPart
          do k= 1, nu
            do i = 1, this%NComponents
              pi => this%Interaction(nc, i)
              n2 = pi%NPart2*pi%NUnit2
              pi%EPot(nup+k, 1:n2) = pi%EPot(n1+k, 1:n2)
              if ( this%OptPressure ) then
                pi%Virial(nup+k, 1:n2) = pi%Virial(n1+k, 1:n2)
              end if
              this%Interaction(i, nc)%EPot(1:n2, nup+k) = pi%EPot(n1+k, 1:n2)
              if ( this%OptPressure ) then
                this%Interaction(i, nc)%Virial(1:n2, nup+k) = pi%Virial(n1+k, 1:n2)
              end if
            end do
          end do
          ! Zero diagonal elements
          this%Interaction(nc, nc)%EPot(nup+1:nup+nu, nup+1:nup+nu) = 0._RK
          if ( this%OptPressure ) then
            this%Interaction(nc, nc)%Virial(nup+1:nup+nu, nup+1:nup+nu) = 0._RK
          end if

          this%NPart = this%NPart - 1
          this%NUnitTotal = this%NUnitTotal - nu
          ! Update density
          this%Density = this%NPart / this%Volume0

          ! Update fractions and NDF
          call UpdateFractions( this )
          ! Update long range correction
          call CalculateCorr( this )
        else        ! Rejection
          call EwaldSelfTerm_Energy (this)
          call EwaldFourierEnergy(this,nc,np,1)
        end if

#if SPME > 0
      else if (LongRange .eq. PME) then
        EFourier = this%UFourier
        !EVirial  = this%EVirial
        !EVirialIntra = this%EVirialIntra
        USelf    = this%USelbstTerm
        UIntra   = this%UIntra
        this%qgrida_old = this%qgrida
        call chargegrid_min( this, nc,np )
        this%NPart = this%NPart - 1
        this%NUnitTotal = this%NUnitTotal - nu
        this%Component(nc)%NPart = this%Component(nc)%NPart - 1
        call PMESelfTermMC( this )
  ! For further use of the following code
        this%NPart = this%NPart + 1
        this%NUnitTotal = this%NUnitTotal + nu
        this%Component(nc)%NPart = this%Component(nc)%NPart + 1
        write(*,*) 'Molecule Deletion is not supported yet with PME'
        STOP
#endif


! ReactionField
      else
        ! Calculate particle energy
        if (SimulationType .eq. MonteCarlo) then
#if MPI_VER > 0
          if ( Equilibration .and. CommonEqui ) then
            ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
            call MPI_Allreduce( GetEnergy( this, nc, np ), EPotDel, 1,MPI_RK, MPI_SUM, Communicator, ierror )
          else
            EPotDel = GetEnergy( this, nc, np )
          endif 
#else
          EPotDel = GetEnergy( this, nc, np )
#endif
        else
          E = 0._RK; EIntra = 0._RK; EBond = 0._RK; EAngle = 0._RK; EDihedral = 0._RK; FDel(:,:) = 0._RK;
          do k = 1, this%NComponents
            if (k > nc) then
              call MDEnergy( this%Interaction(nc,k), np, nu, FDel(:,1:nu), E, EIntra, EBond, EAngle, EDihedral, this%BoxLength, .true. )
            else
              call MDEnergy( this%Interaction(k,nc), np, nu, FDel(:,1:nu), E, EIntra, EBond, EAngle, EDihedral, this%BoxLength, .false. )
            end if
          end do
          EPotDel = E - EIntra
        end if

        EPotDel = EPotDel + this%Density * pc%EPotTestCorrLJ + pc%EPotTestCorrRF
        DensityCorr = this%Density
        if (Shake > 0) DensityCorr = (this%NUnitTotal-this%constrNDF/3._RK) / this%Volume0

        ! Apply acceptance criterion
        if( rnd( 0._RK, 1._RK ) .lt. ( exp( EPotDel / this%RefTemperature - pc%ChemPot ) * DensityCorr * pc%Fraction )) then

          ! Accept Deletion
          this%NDeleteSuccesses = this%NDeleteSuccesses + 1
          call RemoveParticle( pc, np )

          if (SimulationType .ne. MonteCarlo) then
            success = .true.
          else
            ! Copy energies and virial
            nup = (np-1)*nu
            n1 = pc%NPart
            do k= 1, nu
              do i = 1, this%NComponents
                pi => this%Interaction(nc, i)
                n2 = pi%NPart2*pi%NUnit2
                pi%EPot(nup+k, 1:n2) = pi%EPot(n1+k, 1:n2)
                if ( this%OptPressure ) then
                  pi%Virial(nup+k, 1:n2) = pi%Virial(n1+k, 1:n2)
                end if
                this%Interaction(i, nc)%EPot(1:n2, nup+k) = pi%EPot(n1+k, 1:n2)
                if ( this%OptPressure ) then
                  this%Interaction(i, nc)%Virial(1:n2, nup+k) = pi%Virial(n1+k, 1:n2)
                end if
              end do
            end do
            ! Zero diagonal elements
            this%Interaction(nc, nc)%EPot(nup+1:nup+nu, nup+1:nup+nu) = 0._RK
            if ( this%OptPressure ) then
              this%Interaction(nc, nc)%Virial(nup+1:nup+nu, nup+1:nup+nu) = 0._RK
            end if
          end if

          this%NPart = this%NPart - 1
          this%NUnitTotal = this%NUnitTotal - nu
          ! Update density
          this%Density = this%NPart / this%Volume0

          if (SimulationType .eq. MonteCarlo) then
            ! Update fractions and NDF
            call UpdateFractions( this )
            ! Update long range correction
            call CalculateCorr( this )
          end if
        else
        end if

      end if

    end if

    if (SimulationType .ne. MonteCarlo) then
#if MPI_VER > 0
      call MPI_Bcast( success, 1, MPI_LOGICAL, NRootProc, Communicator, ierror )
      if (success) then
        call MPI_Bcast( this%NPart, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NUnitTotal, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( pc%NPart, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%Density, 1, MPI_RK, NRootProc, Communicator, ierror )
        call Unit2Atom1( pc, np)
        pc%NPart1 = ProcRange( pc%NPart, pc%NPart0, pc%NPart2 )
        ! Update fractions and NDF
        call UpdateFractions( this )

        ! Update long range correction
        call CalculateCorr( this )
      end if
#else
      if (success) then
        ! Update fractions and NDF
        call UpdateFractions( this )
        ! Update long range correction
        call CalculateCorr( this )
      end if
#endif
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
    integer                     :: i, k, n1, n2, nu, nu1, nu2, nu1k, nu2k
    real(RK)                    :: PSave(3)
    real(RK)                    :: P0Save(3, 1:this%Component(nc)%Molecule%NUnit)
    real(RK)                    :: Q0Save(4, 1:this%Component(nc)%Molecule%NUnit)
    real(RK)                    :: ESave(this%NUnitMax,this%NPartMax*this%NUnitMax)
    real(RK)                    :: VSave(this%NUnitMax,this%NPartMax*this%NUnitMax)

    ! Assign local variables
    pc => this%Component(nc)
    n1 = pc%NPart
    nu = pc%Molecule%NUnit

    ! Copy position and quaternions
    PSave(:) = pc%Pm0(np, :)
    pc%Pm0(np, :) = pc%Pm0(n1, :)
    pc%Pm0(n1, :) = PSave(:)
    P0Save(:,1:nu) = pc%P0(np,:,1:nu )
    pc%P0(np,:,1:nu) = pc%P0(n1,:,1:nu)
    pc%P0(n1,:,1:nu) = P0Save(:,1:nu)

    if( pc%Molecule%IsElongated ) then
      Q0Save(:, 1:nu) = pc%Q0(np, :, 1:nu)
      pc%Q0(np, :, 1:nu) = pc%Q0(n1, :, 1:nu)
      pc%Q0(n1, :, 1:nu) = Q0Save(:, 1:nu)
    end if

    ! Convert molecular coordinates to atom positions
    call Unit2Atom1( pc, np )
    call Unit2Atom1( pc, n1 )

    ! Copy energies and virial
    do i = 1, this%NRealComponents
      pi => this%Interaction(nc, i)
      n2 = pi%NPart2*pi%NUnit2
      nu1 = (np-1)*pi%NUnit1
      nu2 = (n1-1)*pi%NUnit1
      do k=1, pi%NUnit1
        nu1k = nu1 + k
        nu2k = nu2 + k
        ESave(k,1:n2) = pi%EPot(nu1k, :)
        if ( this%OptPressure ) then
          VSave(k,1:n2) = pi%Virial(nu1k, :)
        end if
        if( i .eq. nc ) then
          ESave(k,nu1k) = pi%EPot(nu1k, n2)
          if ( this%OptPressure ) then
            VSave(k,nu1k) = pi%Virial(nu1k, n2)
          end if
        end if
        pi%EPot(nu1k, :) = pi%EPot(nu2k, :)
        this%Interaction(i, nc)%EPot(:, nu1k) = pi%EPot(nu2k, :)
        pi%EPot(nu2k, :) = ESave(k,1:n2)
        this%Interaction(i, nc)%EPot(:, nu2k) = ESave(k,1:n2)
        if ( this%OptPressure ) then
          pi%Virial(nu1k, :) = pi%Virial(nu2k, :)
          this%Interaction(i, nc)%Virial(:, nu1k) = pi%Virial(nu2k, :)
          pi%Virial(nu2k, :) = VSave(k,1:n2)
          this%Interaction(i, nc)%Virial(:, nu2k) = VSave(k,1:n2)
        end if
      end do
    end do

    ! Zero diagonal elements
    do i=1,nu
      this%Interaction(nc, nc)%EPot(nu1+i, nu1+i) = 0._RK
      if ( this%OptPressure ) then
        this%Interaction(nc, nc)%Virial(nu1+i, nu1+i) = 0._RK
      end if
    end do

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
#if SPME > 0
    real(RK) :: UIntra, EVirialintra
#endif
    real(RK) :: DelBoxL,BoxLengthOld
    logical  :: accepted
#if MPI_VER > 0
    real(RK) :: EPotNew
#endif

    ! Update number of resizing attempts
    this%NResizeAttempts = this%NResizeAttempts + 1

    ! Save current simulation box size, volume, energy, virial
    VolumeOld = this%Volume0
    EPotOld = this%EPot
    BoxLengthOld = this%BoxLength
    if (LongRange .eq. Ewald) then
       UFourier= this%UFourier
       if ( this%OptPressure ) then
         EVirial = this%EVirial
       end if

#if SPME > 0
    else if (LongRange .eq. PME) then
       EVirialIntra = this%EVirialIntra
       UFourier= this%UFourier
       EVirial = this%EVirial
       UIntra = this%UIntra
#endif
    end if

    ! Generate a trial volume change
    this%Volume0 = this%Volume0 * (1._RK + rnd( -this%DispVol, this%DispVol ))
    call UpdateBoxLength( this )

    ! Convert molecular coordinates to atom positions
    DelBoxL = this%BoxLength / BoxLengthOld
    call ResizeMol( this, DelBoxL )
    call Unit2Atom( this )

    ! Calculate potential energy and virial at trial position
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      call Energy( this, EPotNew )
      ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
      call MPI_Allreduce( EPotNew, this%EPot, 1, MPI_RK, MPI_SUM, Communicator, ierror )
    else
     call Energy( this, this%EPot )
    endif
#else
    call Energy( this, this%EPot )
#endif

    ! Find potential change
	
	! NPH
    if( EnsembleType .eq. EnsembleTypeNPH ) then
	  if( exp(( real (this%NDF, RK) / 2._RK - 1._RK) * log((this%RefEnthalpy*this%NPart - this%Epot - this%RefPressure * this%Volume0) &
&       / (this%RefEnthalpy*this%NPart - EPotOld - this%RefPressure * VolumeOld)) + this%NPart * log(this%Volume0 / VolumeOld)) > rnd( 0._RK, 1._RK )) then

	    ! Accept volume change
        this%Temperature = 2._RK * (this%RefEnthalpy*this%NPart - this%Epot - this%RefPressure * this%Volume0) / real (this%NDF, RK)
        this%NResizeSuccesses = this%NResizeSuccesses + 1
        call UpdateEnergy( this )

#if MPI_VER > 0
          if ( Equilibration .and. CommonEqui ) then         !Michael Sch.: move to after if clause!
            ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
            call MPI_Allreduce( GetEnergyIntra( this ), this%EPotIntra, 1, MPI_RK, MPI_SUM, Communicator, ierror )
            if (printIDF) then
              call MPI_Allreduce( GetEnergyIntra_Bond( this ), this%EPotIntra_Bond, 1, MPI_RK, MPI_SUM, Communicator, ierror )
              call MPI_Allreduce( GetEnergyIntra_Angle( this ), this%EPotIntra_Angle, 1, MPI_RK, MPI_SUM, Communicator, ierror )
              call MPI_Allreduce( GetEnergyIntra_Dihedral( this ), this%EPotIntra_Dihedral, 1, MPI_RK, MPI_SUM, Communicator, ierror )
              this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
            endif
            this%EPotInter   = this%EPot - this%EPotIntra
            call MPI_Allreduce( Getd2EpotdV2( this ), this%d2EpotdV2, 1, MPI_RK, MPI_SUM, Communicator, ierror )
            if ( this%OptPressure ) then
              call MPI_Allreduce( GetVirial( this ), this%Virial, 1, MPI_RK, MPI_SUM, Communicator, ierror )
              call MPI_Allreduce( GetVirialIntra(this), this%VirialIntra, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
              this%VirialInter = this%Virial - this%VirialIntra
            end if
          else
            this%EPotIntra   = GetEnergyIntra( this )
            if (printIDF) then
              this%EpotIntra_Bond = GetEnergyIntra_Bond( this )
              this%EpotIntra_Angle = GetEnergyIntra_Angle( this )
              this%EpotIntra_Dihedral = GetEnergyIntra_Dihedral( this )
              this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
            endif
            this%EPotInter   = this%EPot - this%EPotIntra
            this%d2EpotdV2 = Getd2EpotdV2( this )
            if ( this%OptPressure ) then
              this%Virial = GetVirial( this )
              this%VirialIntra = GetVirialIntra( this )
              this%VirialInter = this%Virial - this%VirialIntra
            end if
          endif
#else
          this%EPotIntra   = GetEnergyIntra( this )
          if (printIDF) then
            this%EpotIntra_Bond = GetEnergyIntra_Bond( this )
            this%EpotIntra_Angle = GetEnergyIntra_Angle( this )
            this%EpotIntra_Dihedral = GetEnergyIntra_Dihedral( this )
            this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
          endif
          this%EPotInter   = this%EPot - this%EPotIntra
          this%d2EpotdV2 = Getd2EpotdV2( this )
          if ( this%OptPressure ) then
            this%Virial = GetVirial( this )
            this%VirialIntra = GetVirialIntra( this )
            this%VirialInter = this%Virial - this%VirialIntra
          end if
#endif
	  else
        ! Reject volume change
        this%Volume0 = VolumeOld
        call UpdateBoxLength( this )
        call ResizeMol( this, 1._RK / DelBoxL )
        call Unit2Atom( this )
        this%EPot = EPotOld
        if (LongRange .eq. Ewald) then
          this%UFourier = UFourier
          call Energy(this,this%Epot)

#if MPI_VER > 0
          if ( Equilibration .and. CommonEqui ) then
            call MPI_Allreduce( GetEnergy( this ), this%EPot, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
            call MPI_Allreduce( GetEnergyIntra( this ), this%EPotIntra, 1, MPI_RK, MPI_SUM, Communicator, ierror )
            if (printIDF) then
              call MPI_Allreduce( GetEnergyIntra_Bond( this ), this%EPotIntra_Bond, 1, MPI_RK, MPI_SUM, Communicator, ierror )
              call MPI_Allreduce( GetEnergyIntra_Angle( this ), this%EPotIntra_Angle, 1, MPI_RK, MPI_SUM, Communicator, ierror )
              call MPI_Allreduce( GetEnergyIntra_Dihedral( this ), this%EPotIntra_Dihedral, 1, MPI_RK, MPI_SUM, Communicator, ierror )
              this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
            endif
            this%EPotInter   = this%EPot - this%EPotIntra
            call MPI_Allreduce( Getd2EpotdV2( this ), this%d2EpotdV2, 1, MPI_RK, MPI_SUM, Communicator, ierror )
            if ( this%OptPressure ) then
              call MPI_Allreduce( GetVirial( this ), this%Virial, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
              call MPI_Allreduce( GetVirialIntra( this ), this%VirialIntra, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
              this%VirialInter = this%Virial - this%VirialIntra
            end if
          else
            this%EPot = GetEnergy(this)
            this%EPotIntra   = GetEnergyIntra( this )
            if (printIDF) then
              this%EpotIntra_Bond = GetEnergyIntra_Bond( this )
              this%EpotIntra_Angle = GetEnergyIntra_Angle( this )
              this%EpotIntra_Dihedral = GetEnergyIntra_Dihedral( this )
              this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
            endif
            this%EPotInter   = this%EPot - this%EPotIntra
            this%d2EpotdV2 = Getd2EpotdV2( this )
            if ( this%OptPressure ) then
              this%Virial = GetVirial( this )
              this%VirialIntra = GetVirialIntra( this )
              this%VirialInter = this%Virial - this%VirialIntra
            end if
          end if

#else
          this%EPot = GetEnergy(this)
          this%EPotIntra   = GetEnergyIntra( this )
          if (printIDF) then
            this%EpotIntra_Bond = GetEnergyIntra_Bond( this )
            this%EpotIntra_Angle = GetEnergyIntra_Angle( this )
            this%EpotIntra_Dihedral = GetEnergyIntra_Dihedral( this )
            this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
          endif
          this%EPotInter   = this%EPot - this%EPotIntra
          this%d2EpotdV2 = Getd2EpotdV2( this )
          if ( this%OptPressure ) then
            this%Virial = GetVirial( this )
            this%VirialIntra = GetVirialIntra( this )
            this%VirialInter = this%Virial - this%VirialIntra
          end if
#endif

#if SPME > 0
        else if (LongRange .eq. PME) then
          this%UIntra = UIntra
          this%EVirialIntra = EVirialIntra
          this%UFourier = UFourier
          this%EVirial = EVirial
          call charge_grid_MCall ( this )
#endif
        end if
      end if
	
	else !NPT
      EPotDelta = this%RefPressure * (this%Volume0 - VolumeOld) + this%EPot - EPotOld &
&     + this%NPart * this%Temperature * log( VolumeOld / this%Volume0 )

      accepted = EPotDelta < 0._RK
      if ( .not. accepted ) accepted = exp( -EPotDelta / this%Temperature ) > rnd( 0._RK, 1._RK )

      if( accepted ) then

        ! Accept volume change
        this%NResizeSuccesses = this%NResizeSuccesses + 1

        ! Update energy and virial matrices
        call UpdateEnergy( this )

#if MPI_VER > 0
          if ( Equilibration .and. CommonEqui ) then         !Michael Sch.: move to after if clause!
            ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
            call MPI_Allreduce( GetEnergyIntra( this ), this%EPotIntra, 1, MPI_RK, MPI_SUM, Communicator, ierror )
            if (printIDF) then
              call MPI_Allreduce( GetEnergyIntra_Bond( this ), this%EPotIntra_Bond, 1, MPI_RK, MPI_SUM, Communicator, ierror )
              call MPI_Allreduce( GetEnergyIntra_Angle( this ), this%EPotIntra_Angle, 1, MPI_RK, MPI_SUM, Communicator, ierror )
              call MPI_Allreduce( GetEnergyIntra_Dihedral( this ), this%EPotIntra_Dihedral, 1, MPI_RK, MPI_SUM, Communicator, ierror )
              this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
            endif
            this%EPotInter   = this%EPot - this%EPotIntra
            call MPI_Allreduce( Getd2EpotdV2( this ), this%d2EpotdV2, 1, MPI_RK, MPI_SUM, Communicator, ierror )
            if ( this%OptPressure ) then
              call MPI_Allreduce( GetVirial( this ), this%Virial, 1, MPI_RK, MPI_SUM, Communicator, ierror )
              call MPI_Allreduce( GetVirialIntra(this), this%VirialIntra, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
              this%VirialInter = this%Virial - this%VirialIntra
            end if
          else
            this%EPotIntra   = GetEnergyIntra( this )
            if (printIDF) then
              this%EpotIntra_Bond = GetEnergyIntra_Bond( this )
              this%EpotIntra_Angle = GetEnergyIntra_Angle( this )
              this%EpotIntra_Dihedral = GetEnergyIntra_Dihedral( this )
              this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
            endif
            this%EPotInter   = this%EPot - this%EPotIntra
            this%d2EpotdV2 = Getd2EpotdV2( this )
            if ( this%OptPressure ) then
              this%Virial = GetVirial( this )
              this%VirialIntra = GetVirialIntra( this )
              this%VirialInter = this%Virial - this%VirialIntra
            end if
          endif
#else
          this%EPotIntra   = GetEnergyIntra( this )
          if (printIDF) then
            this%EpotIntra_Bond = GetEnergyIntra_Bond( this )
            this%EpotIntra_Angle = GetEnergyIntra_Angle( this )
            this%EpotIntra_Dihedral = GetEnergyIntra_Dihedral( this )
            this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
          endif
          this%EPotInter   = this%EPot - this%EPotIntra
          this%d2EpotdV2 = Getd2EpotdV2( this )
          if ( this%OptPressure ) then
            this%Virial = GetVirial( this )
            this%VirialIntra = GetVirialIntra( this )
            this%VirialInter = this%Virial - this%VirialIntra
          end if
#endif

      else

        ! Reject volume change
        this%Volume0 = VolumeOld
        call UpdateBoxLength( this )
        call ResizeMol( this, 1._RK / DelBoxL )
        call Unit2Atom( this )
        this%EPot = EPotOld
        if (LongRange .eq. Ewald) then
          this%UFourier = UFourier
          call Energy(this,this%Epot)

#if MPI_VER > 0
          if ( Equilibration .and. CommonEqui ) then
            call MPI_Allreduce( GetEnergy( this ), this%EPot, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
            call MPI_Allreduce( GetEnergyIntra( this ), this%EPotIntra, 1, MPI_RK, MPI_SUM, Communicator, ierror )
            if (printIDF) then
              call MPI_Allreduce( GetEnergyIntra_Bond( this ), this%EPotIntra_Bond, 1, MPI_RK, MPI_SUM, Communicator, ierror )
              call MPI_Allreduce( GetEnergyIntra_Angle( this ), this%EPotIntra_Angle, 1, MPI_RK, MPI_SUM, Communicator, ierror )
              call MPI_Allreduce( GetEnergyIntra_Dihedral( this ), this%EPotIntra_Dihedral, 1, MPI_RK, MPI_SUM, Communicator, ierror )
              this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
            endif
            this%EPotInter   = this%EPot - this%EPotIntra
            call MPI_Allreduce( Getd2EpotdV2( this ), this%d2EpotdV2, 1, MPI_RK, MPI_SUM, Communicator, ierror )
            if ( this%OptPressure ) then
              call MPI_Allreduce( GetVirial( this ), this%Virial, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
              call MPI_Allreduce( GetVirialIntra( this ), this%VirialIntra, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
              this%VirialInter = this%Virial - this%VirialIntra
            end if
          else
            this%EPot = GetEnergy(this)
            this%EPotIntra   = GetEnergyIntra( this )
            if (printIDF) then
              this%EpotIntra_Bond = GetEnergyIntra_Bond( this )
              this%EpotIntra_Angle = GetEnergyIntra_Angle( this )
              this%EpotIntra_Dihedral = GetEnergyIntra_Dihedral( this )
              this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
            endif
            this%EPotInter   = this%EPot - this%EPotIntra
            this%d2EpotdV2 = Getd2EpotdV2( this )
            if ( this%OptPressure ) then
              this%Virial = GetVirial( this )
              this%VirialIntra = GetVirialIntra( this )
              this%VirialInter = this%Virial - this%VirialIntra
            end if
          end if

#else
          this%EPot = GetEnergy(this)
          this%EPotIntra   = GetEnergyIntra( this )
          if (printIDF) then
            this%EpotIntra_Bond = GetEnergyIntra_Bond( this )
            this%EpotIntra_Angle = GetEnergyIntra_Angle( this )
            this%EpotIntra_Dihedral = GetEnergyIntra_Dihedral( this )
            this%EpotIntra_Nonbonded = this%EPotIntra - this%EPotIntra_Bond - this%EPotIntra_Angle - this%EPotIntra_Dihedral
          endif
          this%EPotInter   = this%EPot - this%EPotIntra
          this%d2EpotdV2 = Getd2EpotdV2( this )
          if ( this%OptPressure ) then
            this%Virial = GetVirial( this )
            this%VirialIntra = GetVirialIntra( this )
            this%VirialInter = this%Virial - this%VirialIntra
          end if
#endif

#if SPME > 0
        else if (LongRange .eq. PME) then
          this%UIntra = UIntra
          this%EVirialIntra = EVirialIntra
          this%UFourier = UFourier
          this%EVirial = EVirial
          call charge_grid_MCall ( this )
#endif
        end if
      end if
    end if

  end subroutine TEnsemble_Resize



!==============================================================!
!  Subroutine TEnsemble_Resize_LiquidPhase                     !
!==============================================================!

  subroutine TEnsemble_Resize_liq( this,dv,EPotDelta )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK), intent(in out) :: dv
    real(RK), intent(in out) :: EPotDelta
    real(RK) :: VolumeOld, EPotOld, BoxLengthOld, DelBoxL
#if MPI_VER > 0
    real(RK) :: EPotNew
#endif

    ! Update number of resizing attempts
    this%NResizeAttempts = this%NResizeAttempts + 1

    ! Save current simulation energy
    EPotOld = this%EPot
    VolumeOld = this%Volume0

    ! Generate a trial volume change
    dv = this%Volume0 * rnd( -this%DispVol, this%DispVol )
    this%Volume0 = this%Volume0 + dv
    BoxLengthOld = this%BoxLength
    call UpdateBoxLength( this )

    ! Convert molecular coordinates to atom positions
    DelBoxL = this%BoxLength / BoxLengthOld
    call ResizeMol( this, DelBoxL ) ! testing needed before old mol2unit was used (before rev388)
    call Unit2Atom( this )

    ! Calculate potential energy and virial at trial position
#if MPI_VER > 0
    call Energy( this, EPotNew )
    call MPI_Allreduce( EPotNew, this%EPot, 1, MPI_RK, MPI_SUM, Communicator, ierror )
#else
    call Energy( this, this%EPot )
#endif

    ! Find potential change
    EPotDelta = this%EPot - EPotOld + this%NPart * this%Temperature * log( VolumeOld / this%Volume0 )

  end subroutine TEnsemble_Resize_liq



!==============================================================!
!  Subroutine TEnsemble_Resize_LiquidPhaseUpdate               !
!==============================================================!
  subroutine TEnsemble_ResizeLiquid_Update(this,accept,EPotOldliq,VolumeOld)

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this
    real(RK),intent(in) :: EPotOldliq,VolumeOld
    logical  :: accept
    real(RK) :: BoxLengthOld, DelBoxL

    if ( accept ) then
      ! Accept volume change
      this%NResizeSuccesses = this%NResizeSuccesses + 1

      ! Update energy and virial matrices
      call UpdateEnergy( this )

#if MPI_VER > 0
      call MPI_Allreduce( GetVirial( this ), this%Virial, 1, MPI_RK, MPI_SUM, Communicator, ierror )
#else
      this%Virial = GetVirial( this )
#endif

    else
      ! Reject volume change
      this%Volume0 = VolumeOld
      BoxLengthOld = this%BoxLength ! correct
      call UpdateBoxLength( this )
      DelBoxL = this%BoxLength / BoxLengthOld
      call ResizeMol( this, DelBoxL ) ! testing needed before old mol2unit was used (before rev388)
      call Unit2Atom( this )
      this%EPot = EPotOldliq
      if (LongRange .eq. Ewald) then
         call Energy(this,this%Epot)

#if MPI_VER > 0
         call MPI_Allreduce( GetEnergy( this ), this%EPot, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
         call MPI_Allreduce( GetVirial( this ), this%Virial, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
#else
         this%EPot = GetEnergy(this)
         this%Virial = GetVirial( this )
#endif

#if SPME > 0
      else if (LongRange .eq. PME) then
        call charge_grid_MCall ( this )
#endif
      end if
    end if

  end subroutine TEnsemble_ResizeLiquid_Update


!==============================================================!
!  Subroutine TEnsemble_Resize_VaporPhase                      !
!==============================================================!

  subroutine TEnsemble_Resize_vap( this,dv,EPotDelta,accept )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK), intent(in) :: dv
    real(RK), intent(in out) :: EPotDelta
    real(RK) :: VolumeOld, EPotOld
    real(RK) :: EVirial
    real(RK) :: UFourier, BoxLengthOld, DelBoxL
#if SPME > 0
    real(RK) :: UIntra, EVirialintra
#endif
#if MPI_VER > 0
    real(RK) :: EPotNew
#endif
    logical  :: accept

    ! Save current simulation box size, volume, energy, virial
    VolumeOld = this%Volume0
    EPotOld = this%EPot
    if (LongRange .eq. Ewald) then
       UFourier= this%UFourier
       EVirial = this%EVirial

#if SPME > 0
    else if (LongRange .eq. PME) then
       EVirialIntra = this%EVirialIntra
       UFourier= this%UFourier
       EVirial = this%EVirial
       UIntra = this%UIntra
#endif
    end if

    ! Generate a trial volume change
    this%Volume0 = this%Volume0 - dv
    BoxLengthOld = this%BoxLength
    call UpdateBoxLength( this )

    ! Convert molecular coordinates to atom positions
    DelBoxL = this%BoxLength / BoxLengthOld
    call ResizeMol( this, DelBoxL ) ! testing needed before old mol2unit was used (before rev388)
    call Unit2Atom( this )

    ! Calculate potential energy and virial at trial position

#if MPI_VER > 0
    call Energy( this, EPotNew )
    call MPI_Allreduce( EPotNew, this%EPot, 1, MPI_RK, MPI_SUM, Communicator, ierror )
#else
    call Energy( this, this%EPot )
#endif

    ! Find potential change
    EPotDelta = EPotDelta + this%EPot - EPotOld + this%NPart * this%Temperature * log( VolumeOld / this%Volume0 )

    ! Acceptance criteria
    if( exp( -EPotDelta / this%Temperature ) .gt. rnd( 0._RK, 1._RK ) ) then

      accept = .true.

      ! Update energy and virial matrices
      call UpdateEnergy( this )

#if MPI_VER > 0
      call MPI_Allreduce( GetVirial( this ), this%Virial, 1, MPI_RK, MPI_SUM, Communicator, ierror )
#else
      this%Virial = GetVirial( this )
#endif

    else
      ! Reject volume change
      this%Volume0 = VolumeOld
      call UpdateBoxLength( this )
      call ResizeMol( this, 1._RK/DelBoxL ) ! testing needed before old mol2unit was used (before rev388)
      call Unit2Atom( this )
      this%EPot = EPotOld
      if (LongRange .eq. Ewald) then
         this%UFourier = UFourier
         call Energy(this,this%Epot)

#if MPI_VER > 0
         call MPI_Allreduce( GetEnergy( this ), this%EPot, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
         call MPI_Allreduce( GetVirial( this ), this%Virial, 1 , MPI_RK, MPI_SUM, Communicator, ierror )
#else
         this%EPot = GetEnergy(this)
         this%Virial = GetVirial( this )
#endif

#if SPME > 0
      else if (LongRange .eq. PME) then
         this%UIntra = UIntra
         this%EVirialIntra = EVirialIntra
         this%UFourier = UFourier
         this%EVirial = EVirial
         call charge_grid_MCall ( this )
#endif
      end if

    end if

  end subroutine TEnsemble_Resize_vap



!==============================================================!
!  Subroutine TEnsemble_GibbsRemoveParticle                    !
!==============================================================!

  subroutine TEnsemble_GibbsRemove( this, nc, np, EPotDelta )

    implicit none

    ! Declare arguments
    type(TEnsemble)         :: this

    real(RK)                :: rx,sx
    real(RK),intent(in out) :: EPotDelta
    integer,intent(in out)  :: nc,np

    real(RK) :: charge

     charge = 1._RK

     ! Chose component - exclude charged components
     do while  ( Charge > 0.1_RK )
       sx = 0._RK
       rx = rnd( 0._RK, 1._RK )
       do nc = 1, this%NComponents
          sx = sx + this%Component(nc)%Fraction
          if( rx <= sx ) exit 
       end do 
       Charge = this%Component(nc)%Molecule%Charge
     end do

     np = rnd( this%Component(nc)%NPart )

     call Gibbs_Delete( this, nc, np, EPotDelta )

   end subroutine TEnsemble_GibbsRemove


!==============================================================!
!  Subroutine TEnsemble_GibbsDelete                            !
!==============================================================!

  subroutine TEnsemble_GibbsDelete( this, nc, np, EPotDel )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np
    real(RK),intent(in out) :: EPotDel

    ! Declare local variables
    type(TComponent), pointer   :: pc
    type(TInteraction), pointer :: pi
    integer                     :: i, n1, n2

! Ewald Parameter
    real(RK)                    :: EFourier, EPotNew
    real(RK)                    :: EVirial
#if SPME > 0
    real(RK)                    :: EVirialIntra
#endif
    real(RK)                    :: USelf, UIntra
    real(RK)                    :: r(3)
    real(RK)                    :: q(4)

    ! Assign local variables
    pc => this%Component(nc)
    EPotDel = 0._RK

    ! Update number of delete attempts
    this%NDeleteAttempts = this%NDeleteAttempts + 1


    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      USelf    = this%USelbstTerm
      UIntra   = this%UIntra
      pc%NPart = pc%NPart - 1
      call EwaldSelfTerm_Energy (this)
      call EwaldFourierEnergy(this,nc,np,-1)
      pc%NPart = pc%NPart + 1
      ! Calculate particle energy

#if MPI_VER > 0
      call MPI_Allreduce( GetEnergy( this, nc, np ), EPotDel, 1, MPI_RK, MPI_SUM, Communicator, ierror )
#else
      EPotDel = GetEnergy( this, nc, np )
#endif

      EPotDel = EPotDel + this%Density * pc%EPotTestCorrLJ + NProcs*(this%UIntra-UIntra + this%USelbstTerm-USelf-EFourier) - &
&                  this%Temperature*log(this%Volume0/(this%NPart) )

#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      EVirialIntra = this%EVirialIntra
      USelf    = this%USelbstTerm
      UIntra   = this%UIntra
      this%qgrida_old = this%qgrida
      call chargegrid_min ( this, nc,np )
      this%NPart = this%NPart - 1
      this%NUnitTotal = this%NUnitTotal - pc%Molecule%NUnit
      this%Component(nc)%NPart = this%Component(nc)%NPart - 1
      call PMESelfTermMC ( this )
! For further use of the following code
      this%NPart = this%NPart + 1
      this%NUnitTotal = this%NUnitTotal + pc%Molecule%NUnit
      this%Component(nc)%NPart = this%Component(nc)%NPart + 1
      write(*,*) 'Molecule Deletion is not supported yet with PME'
      STOP
#endif


! ReactionField
    else
      ! Calculate particle energy
#if MPI_VER > 0
      call MPI_Allreduce( GetEnergy( this, nc, np ), EPotDel, 1, MPI_RK, MPI_SUM, Communicator, ierror )
#else
      EPotDel = GetEnergy( this, nc, np )
#endif
      EPotDel = EPotDel + this%Density * pc%EPotTestCorrLJ + pc%EPotTestCorrRF - this%Temperature*log(this%Volume0/(this%NPart) )

    end if

  end subroutine TEnsemble_GibbsDelete


!==============================================================!
!  Subroutine TEnsemble_GibbsInsertParticle                    !
!==============================================================!

  subroutine TEnsemble_GibbsInsert( this, nc, EPotDelta,accept )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc
    real(RK),intent(in out)   :: EPotDelta
    logical             :: accept

    ! Declare local variables
    real(RK)                  :: r(3)
    real(RK)                  :: q(3)
    type(TComponent), pointer :: pc
    integer                   :: i, np
    real(RK)                  :: EPotIns
    real(RK)                  :: UIntra, USelbst, EFourier, EVirial
#if MPI_VER > 0
    real(RK)                  :: EPotInsAll
#endif

    ! Assign local variables
    pc => this%Component(nc)

    ! Generate a random position and orientation
    do i = 1, 3
      r(i) = rnd( -.5_RK, .5_RK )
    end do
    do i = 1, 3
      q(i) = rnd( -1._RK, 1._RK )
    end do

    call AddParticle( pc, r, q )
    np = pc%NPart
    this%NPart = this%NPart + 1
    this%NUnitTotal = this%NUnitTotal + pc%Molecule%NUnit

    ! Convert unit coordinates to atom positions
    call Unit2Atom1( pc, np )

    if (LongRange .eq. Ewald) then           ! EWALD-SUMMATION
      UIntra   = this%UIntra
      USelbst  = this%USelbstTerm
      EFourier = this%UFourier
      EVirial  = this%EVirial
      ! Energy
      call EwaldSelfTerm_Energy(this)
      call Energy ( this, nc, np, EPotIns, 1 )
#if MPI_VER > 0
      call MPI_Allreduce( EPotIns, EPotInsAll, 1, MPI_RK, MPI_SUM, Communicator, ierror )
      EPotDelta = EpotDelta - EPotInsAll - this%Density * pc%EPotTestCorrLJ - this%Temperature * log((this%NPart)/this%Volume0 ) - &
&            NProcs * this%UIntra + NProcs * UIntra - NProcs * this%USelbstTerm + NProcs * USelbst + NProcs * EFourier

      if( rnd( 0._RK, 1._RK ) .lt. ( exp( EPotDelta / this%Temperature ) )) then
#else
      EPotDelta = EPotDelta - EPotIns - this%Density * pc%EPotTestCorrLJ  - this%Temperature * log((this%NPart)/this%Volume0 ) - &
&            this%UIntra + UIntra - this%USelbstTerm + USelbst + EFourier

      if( rnd( 0._RK, 1._RK ) .lt. ( exp( EPotDelta / this%Temperature ) )) then
#endif
        ! Accept Insertion
        accept = .true.
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
        this%NUnitTotal = this%NUnitTotal - pc%Molecule%NUnit
        call EwaldFourierEnergy ( this, nc, np, -1 )
        this%USelbstTerm = USelbst
        this%UIntra  = UIntra
      end if 

    else                                         ! REACTION FIELD
      ! Calculate particle energy at trial position
      call Energy( this, nc, np, EPotIns )
    ! Apply acceptance criterion
#if MPI_VER > 0
      call MPI_Allreduce( EPotIns, EPotInsAll, 1, MPI_RK, MPI_SUM, Communicator, ierror )

      EPotDelta = EpotDelta - EPotInsAll - this%Density * pc%EPotTestCorrLJ &
&        - pc%EPotTestCorrRF - this%Temperature*log((this%NPart)/this%Volume0 )

      if( rnd( 0._RK, 1._RK ) .lt. ( exp( EPotDelta / this%Temperature ) )) then
#else
      EPotDelta = EPotDelta - EPotIns - this%Density * pc%EPotTestCorrLJ &
&         - pc%EPotTestCorrRF - this%Temperature*log((this%NPart)/this%Volume0 )

      if( rnd( 0._RK, 1._RK ) .lt. ( exp( EPotDelta / this%Temperature ) )) then
#endif

        accept = .true.
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
        this%NUnitTotal = this%NUnitTotal - pc%Molecule%NUnit
      end if 

    end if

   end subroutine TEnsemble_GibbsInsert



!==============================================================!
!  Subroutine TEnsemble_PartChangeUpdate                       !
!==============================================================!
  subroutine TEnsemble_PartChangeUpdate(this,nc,np,TransferRate,accept)

    implicit none

    ! Declare arguments
    type(TEnsemble)       :: this
    integer, intent(in)   :: nc,np
    integer,intent(in out):: TransferRate
    logical               :: accept

    ! Declare variables
    type(TComponent)  ,pointer :: pc
    type(TInteraction),pointer :: pi

    integer  :: n1,n2,i
    real(RK) :: AccRateTransfer

    pc => this%Component(nc)
    this%NTransferAttempts = this%NTransferAttempts + 1
    if ( accept ) then
      this%NTransferSuccesses = this%NTransferSuccesses + 1
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
        this%NUnitTotal = this%NUnitTotal - pc%Molecule%NUnit

        ! Update density
        this%Density = this%NPart / this%Volume0

        ! Update fractions and NDF
        call UpdateFractions( this )

        ! Update long range correction
        call CalculateCorr( this )
    else
      if (LongRange .eq. Ewald) then
        call EwaldSelfTerm_Energy (this)
        call EwaldFourierEnergy(this,nc,np,1)
      end if
    endif

! Adjusting the amount of transfer moves
    AccRateTransfer = real(this%NTransferSuccesses) / real(this%NTransferAttempts)

      ! Update transfer rate
      if(( AccRateTransfer .gt. AccUpperLimit ) .and. ( TransferRate .lt. TransferRateLimit )) then
           TransferRate = int(TransferRate * 1.05_RK)

      else if( AccRateTransfer .lt. AccLowerLimit ) then
           TransferRate = int(TransferRate * 0.95_RK)
      end if



  end subroutine TEnsemble_PartChangeUpdate


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
    this%NTransferAttempts  = 0
    this%NTransferSuccesses = 0

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

      if (pc%Molecule%NUnit .eq. 1) then ! only one type of moves - molecular, no unit differentiation
        ! Update translational displacement
        if(( AccRateTran .gt. AccUpperLimit) .and. ( pc%DispTran .lt. DispMolTranUppLimit )) then
          pc%DispTran = pc%DispTran * 1.05_RK
        else if(( AccRateTran .lt. AccLowerLimit ) .and. ( pc%DispTran .gt. DispMolTranLowLimit )) then
          pc%DispTran = pc%DispTran * .95_RK
        end if

        ! Update rotational displacement
        if(( AccRateRot .gt. AccUpperLimit ) .and. ( pc%DispRot .lt. DispMolRotUppLimit )) then
          pc%DispRot = pc%DispRot * 1.05_RK

        else if(( AccRateRot .lt. AccLowerLimit ) .and. ( pc%DispRot .gt. DispMolRotLowLimit )) then
          pc%DispRot = pc%DispRot * 0.95_RK
        end if

      else
        ! Update translational displacement
        if(( AccRateTran .gt. AccUpperLimit) .and. ( pc%DispTran .lt. DispTranUppLimit )) then
          pc%DispTran = pc%DispTran * 1.05_RK
        else if(( AccRateTran .lt. AccLowerLimit ) .and. ( pc%DispTran .gt. DispTranLowLimit )) then
          pc%DispTran = pc%DispTran * .95_RK
        end if

        ! Update rotational displacement
        if(( AccRateRot .gt. AccUpperLimit ) .and. ( pc%DispRot .lt. DispRotUppLimit )) then
          pc%DispRot = pc%DispRot * 1.05_RK
        else if(( AccRateRot .lt. AccLowerLimit ) .and. ( pc%DispRot .gt. DispRotLowLimit )) then
          pc%DispRot = pc%DispRot * 0.95_RK
        end if

        ! Update Displacements for entire molecule
        AccRateTran = real(pc%NMoveMolSuccesses) / real(pc%NMoveMolAttempts)
        AccRateRot = real(pc%NRotateMolSuccesses) / real(pc%NRotateMolAttempts)

        ! Update molecular translational displacement
        if(( AccRateTran .gt. AccUpperLimit) .and. ( pc%DispMolTran .lt. DispMolTranUppLimit )) then
          pc%DispMolTran = pc%DispMolTran * 1.05_RK
        else if(( AccRateTran .lt. AccLowerLimit ) .and. ( pc%DispMolTran .gt. DispMolTranLowLimit )) then
          pc%DispMolTran = pc%DispMolTran * .95_RK
        end if

        ! Update molecular rotational displacement
        if(( AccRateRot .gt. AccUpperLimit ) .and. ( pc%DispMolRot .lt. DispMolRotUppLimit )) then
          pc%DispMolRot = pc%DispMolRot * 1.05_RK
        else if(( AccRateRot .lt. AccLowerLimit ) .and. ( pc%DispMolRot .gt. DispMolRotLowLimit )) then
          pc%DispMolRot = pc%DispMolRot * 0.95_RK
        end if
      end if
    end do

    if( ConstantPressure .and. .not. NVTEquilibration ) then
      AccRateVol = real(this%NResizeSuccesses) / real(this%NResizeAttempts)

      ! Update volume displacement
      if(( AccRateVol .gt. AccUpperLimit ) .and. ( this%DispVol .lt. DispVolUppLimit )) then
        this%DispVol = this%DispVol * 1.05_RK

      else if(( AccRateVol .lt. AccLowerLimit ) .and. ( this%DispVol .gt. DispVolLowLimit )) then
        this%DispVol = this%DispVol * 0.95_RK
      end if
    end if

    ! Zero attempts
    call ZeroNAttempts( this )

  end subroutine TEnsemble_UpdateDisplacements



!==============================================================!
!  Subroutine TEnsemble_ResidenceTime                          !
!==============================================================!

  subroutine TEnsemble_Residence( this )

    implicit none

    ! Declare arguments
    type(TEnsemble)    :: this

    ! Declare local variables
    real(RK) :: R1x, R1y, R1z
    real(RK) :: R2x, R2y, R2z
    real(RK) :: drx, dry, drz
    real(RK) :: CriticalLength
    integer  :: i
    integer  :: Numb1, Numb2
    integer  :: counter
    
!     CriticalLength = (this%ResidLength * this%BoxLength)**2
    CriticalLength = (this%ResidLength)**2
    counter = 0
    do i = 1, this%ResidPairs
      Numb1 = this%CompPair(i,1)
      Numb2 = this%CompPair(i,2)
      
      ! Do not evaluate the pair, since it seperated earlier
      if ( (Numb1 == 0) .or. (Numb2 == 0) ) then
        counter = counter + 1
        cycle
      end if
      
      ! Calculate distance
      R1x = this%Component(this%ResidComp1)%Molecule%SiteLJ126(this%ResidSite1)%RX(Numb1)
      R1y = this%Component(this%ResidComp1)%Molecule%SiteLJ126(this%ResidSite1)%RY(Numb1)
      R1z = this%Component(this%ResidComp1)%Molecule%SiteLJ126(this%ResidSite1)%RZ(Numb1)
      R2x = this%Component(this%ResidComp2)%Molecule%SiteLJ126(this%ResidSite2)%RX(Numb2)
      R2y = this%Component(this%ResidComp2)%Molecule%SiteLJ126(this%ResidSite2)%RY(Numb2)
      R2z = this%Component(this%ResidComp2)%Molecule%SiteLJ126(this%ResidSite2)%RZ(Numb2)
      
      drx = (R1x - R2x)
      drx = ( (drx -anint(drx))*this%BoxLength )**2
      dry = (R1y - R2y)
      dry = ( (dry -anint(dry))*this%BoxLength )**2
      drz = (R1z - R2z)
      drz = ( (drz -anint(drz))*this%BoxLength )**2

      if ( drx+dry+drz .gt. CriticalLength ) then
        this%ResidCem = this%ResidCem + 1
        this%ResidPairsCem(this%ResidCem,1) = this%CompPair(i,1)
        this%ResidPairsCem(this%ResidCem,2) = this%CompPair(i,2)
        this%ResidPairsCem(this%ResidCem,3) = this%ResidTimesStart(i)
        this%ResidPairsCem(this%ResidCem,4) = Step
        this%CompPair(i,1) = 0
        this%CompPair(i,2) = 0
        this%ResidTimesStart(i) = 0._RK
      end if
    end do
    
    ! Update list, if time is ready
    if (mod(Step,this%ResidPeriod) == 0) then
      call ResidencePartners ( this )
    end if
    ! Count number of pairs per time step
    call Update( this%SumResidencePairs, 1._RK*this%ResidPairs - 1._RK*counter )
    
   end subroutine TEnsemble_Residence

!==============================================================!
!  Subroutine TEnsemble_ResidencePartners                      !
!==============================================================!

  subroutine TEnsemble_ResidencePartners( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this
    
    ! Declare local variables
    integer :: i, j
    integer :: ResidPairs
    integer :: counter 
    type(TComponent),pointer :: pc1, pc2
    real(RK) :: R1x, R1y, R1z
    real(RK) :: R2x, R2y, R2z
    real(RK) :: drx, dry, drz
    real(RK) :: CriticalLength
    
    pc1 => this%Component(this%ResidComp1)
    pc2 => this%Component(this%ResidComp2)
    CriticalLength = (this%ResidLength)**2
    ResidPairs = 0
    this%ResidTimesStart_Old = this%ResidTimesStart
    this%CompPair_Old = this%CompPair
    

    do i=1, pc1%NPart
      do j=1, pc2%NPart

        ! Calculate distance
        R1x = pc1%Molecule%SiteLJ126(this%ResidSite1)%RX(i)
        R1y = pc1%Molecule%SiteLJ126(this%ResidSite1)%RY(i)
        R1z = pc1%Molecule%SiteLJ126(this%ResidSite1)%RZ(i)
        R2x = pc2%Molecule%SiteLJ126(this%ResidSite2)%RX(j)
        R2y = pc2%Molecule%SiteLJ126(this%ResidSite2)%RY(j)
        R2z = pc2%Molecule%SiteLJ126(this%ResidSite2)%RZ(j)
        
        drx = (R1x - R2x)
        drx = ( (drx -anint(drx))*this%BoxLength )**2
        dry = (R1y - R2y)
        dry = ( (dry -anint(dry))*this%BoxLength )**2
        drz = (R1z - R2z)
        drz = ( (drz -anint(drz))*this%BoxLength )**2
        if ( drx+dry+drz .le. CriticalLength ) then
          ResidPairs = ResidPairs + 1
          this%CompPair(ResidPairs,1) = i
          this%CompPair(ResidPairs,2) = j
          this%ResidTimesStart(ResidPairs) = Step
        end if 
      end do
    end do

! Update pairs, that are still grouped
    do i=1,ResidPairs
      do j=1, this%ResidPairs
        if ( (this%CompPair(i,1) .eq. this%CompPair_Old(j,1)) .and. (this%CompPair(i,2) .eq. this%CompPair_Old(j,2)) ) then
          this%ResidTimesStart(i) = this%ResidTimesStart_Old(j)
        end if
      end do
    end do

! Update pairs, that were grouped but were separated for less time than allowed (ResidBreak)
    counter  = 0
    do i=1,ResidPairs
      do j=1, this%ResidCem
        if ( (this%CompPair(i,1) .eq. this%ResidPairsCem(j,1)) .and. (this%CompPair(i,2) .eq. this%ResidPairsCem(j,2)) ) then
          this%ResidTimesStart(i) = this%ResidPairsCem(j,3)
          this%ResidPairsCem(this%ResidCem,1:4) = this%ResidPairsCem(j,1:4)
          counter = counter + 1
        end if
      end do
    end do
    this%ResidCem = this%ResidCem - counter

! Update Residence Time
    counter  = 0
    do i=1, this%ResidCem
      if ( (Step - this%ResidPairsCem(i,4)) .gt. this%ResidBreak ) then
         call Update( this%SumResidenceDuration, ( this%ResidPairsCem(i,4)-this%ResidPairsCem(i,3) )*TimeStep )
          this%ResidPairsCem(this%ResidCem,1:4) = this%ResidPairsCem(i,1:4)
          counter = counter + 1
      end if
    end do
    this%ResidCem = this%ResidCem - counter

! Update Number of pairs in the system (really grouped!)
    this%ResidPairs = ResidPairs

   end subroutine TEnsemble_ResidencePartners



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
#if MPI_VER > 0
      if (SimulationType .eq. MonteCarlo) then
        write( IOBuffer, '(I16)' ) this%EnsembleNumber
        call FileAppend_parallel( this%iounit_result,trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ResultFileExtension )

        if( .not. SimulationType .eq. SecondVirialCoeff ) then

          ! Open running average result file
          write( IOBuffer, '(I16)' ) this%EnsembleNumber
          call FileAppend_parallel( this%iounit_runave, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//RunAveFileExtension )
        end if
      else
        write( IOBuffer, '(I16)' ) this%EnsembleNumber
        call FileAppend( this%iounit_result,trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ResultFileExtension )

        if( .not. SimulationType .eq. SecondVirialCoeff ) then

          ! Open running average result file
          write( IOBuffer, '(I16)' ) this%EnsembleNumber
          call FileAppend( this%iounit_runave, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//RunAveFileExtension )
        end if
      endif
#else
      ! Open result file
      write( IOBuffer, '(I16)' ) this%EnsembleNumber
      call FileAppend( this%iounit_result,trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ResultFileExtension )

      if( .not. SimulationType .eq. SecondVirialCoeff ) then

        ! Open running average result file
        write( IOBuffer, '(I16)' ) this%EnsembleNumber
        call FileAppend( this%iounit_runave, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//RunAveFileExtension )
      end if
#endif
#if TRANS ==1
      write( IOBuffer, '(I16)' ) this%EnsembleNumber
      call FileAppend( this%iounit_rescf, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ResultTransportExtension )

#endif

    else
#if MPI_VER > 0
      if (SimulationType .eq. MonteCarlo) then
        ! Open result file
        write( IOBuffer, '(I16)' ) this%EnsembleNumber
        call FileRewrite_parallel( this%iounit_result, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ResultFileExtension )
        if( .not. SimulationType .eq. SecondVirialCoeff ) then

          ! Open running average result file
          write( IOBuffer, '(I16)' ) this%EnsembleNumber
          call FileRewrite_parallel( this%iounit_runave, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//RunAveFileExtension )
        end if
     else
        ! Open result file
        write( IOBuffer, '(I16)' ) this%EnsembleNumber
        call FileRewrite( this%iounit_result, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ResultFileExtension )

        if( .not. SimulationType .eq. SecondVirialCoeff ) then

          ! Open running average result file
          write( IOBuffer, '(I16)' ) this%EnsembleNumber
          call FileRewrite( this%iounit_runave, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//RunAveFileExtension )
        end if
     endif
#else

      ! Open result file
      write( IOBuffer, '(I16)' ) this%EnsembleNumber
      call FileRewrite( this%iounit_result, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ResultFileExtension )

      if( .not. SimulationType .eq. SecondVirialCoeff ) then

        ! Open running average result file
        write( IOBuffer, '(I16)' ) this%EnsembleNumber
        call FileRewrite( this%iounit_runave, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//RunAveFileExtension )
      end if
#endif
#if  TRANS == 1
      ! Open result file for correlation function
      write( IOBuffer, '(I16)' ) this%EnsembleNumber
      call FileRewrite( this%iounit_rescf, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ResultTransportExtension )
!TRANSPORT_END
#endif

    end if

  end subroutine TEnsemble_ResultOpen



!==============================================================!
!  Subroutine TEnsemble_ResultUpdate                           !
!==============================================================!

  subroutine TEnsemble_ResultUpdate( this )

    implicit none

#if MPI_VER > 0
  include 'mpif.h'
    integer         :: fields = 0
    integer         :: accumulate_step = 0
    integer         :: headers = 0
    integer(kind=MPI_OFFSET_KIND) :: offset = 0 
    integer                       :: ierr
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    type(TComponent), pointer :: pc
    integer                   :: i,j,t,err,currentbin
    real(RK)                  :: value
    real(RK)                  :: currentdEpotdV,currentd2EpotdV2
    real(RK)                  :: A10res, A01res, A20res, A11res, A02res, A20id, A30res, A21res, A12res
    real(RK)                  :: specv, specv2, Beta, Beta2, Beta3, Numb, U, U2, U3, dUdV, UdUdV, dUdV2, U2dUdV, UdUdV2, d2UdV2, Ud2UdV2
    real(RK)                  :: currentHmU, currentHmUm1, currentH
    real(RK)                  :: O10, O01, O20, O11, O02, O30, O21, O12, O40, O31, O22, O00
    real(RK)                  :: S10, S01, S20, S11, S02, S30, S21, S12
    real(RK)                  :: O00m1, O00m2, O00m3, O012, O20m1, S20m1, S20m2, S20m3 
    real(RK)                  :: F, invF, funcF, rho, rho2, HmU, HmUm1, HmUm2, HmUm3, HmUm1dUdV, HmUm1dUdV2, HmUm1d2UdV2, HmUm2dUdV, HmUm2dUdV2, HmUm2d2UdV2, HmUm3dUdV, HmUm3dUdV2
    real(RK)                  :: Momentum(3), Momentumd2Mass, Mass
#if HBOND > 0
    integer                   :: k, l
#endif

    if( Step == 1 ) then
      ! Reset accumulators
      ! 1.) Basic sums
      call Reset( this%SumPressure )
      call Reset( this%SumDensity )
      call Reset( this%SumTemperature )
      call Reset( this%SumEPot )
      call Reset( this%SumEPotInter )
      call Reset( this%SumEPotIntra )
      if (printIDF) then
        call Reset( this%SumEPotIntra_Bond )
        call Reset( this%SumEPotIntra_Angle )
        call Reset( this%SumEPotIntra_Dihedral )
        call Reset( this%SumEPotIntra_Nonbonded )
        call Reset( this%SumVirialIntra )
        call Reset( this%SumVirialInter )
      end if
      call Reset( this%SumEnthalpy )
      call Reset( this%SumConfEnthalpy )
      call Reset( this%SumVolume )
      call Reset( this%SumVirial )
      call Reset( this%SumdEpotdV )
      call Reset( this%Sumd2EpotdV2 )
#if OSMOP > 0
      call Reset( this%SumOsmoticPressure )
      do i = 1, this%NComponents
        do j = 1, NBinsDen
          call Reset( this%Component(i)%SumDenProfile(j) )
#if OSMOP == 2
          if( this%Component(i)%ChemPotMethod .eq. ChemPotMethodWidom ) then
            call Reset( this%Component(i)%SumChemPotProfile(j) )
          end if
#endif
        end do
      end do
#if OSMOP == 2
      do j = 1, NBinsDen
        call Reset( this%SumPressureProfile(j) )
      end do
#endif
#endif

#if HBOND > 0
      do i = 1, this%NComponents
        call Reset( this%SumHBond0(i) )
        do j = 1, this%NComponents
          call Reset( this%SumHBond1(i,j) )
          do k = j, this%NComponents
            call Reset( this%SumHBond2(i,j,k) )
            do l = k, this%NComponents
              call Reset( this%SumHBond3(i,j,k,l) )
            end do
          end do
        end do
        call Reset( this%SumHBondN(i) )
      end do
#endif

      if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
        call Reset( this%SumNPart )

        do i = 1, this%NComponents
          call Reset( this%Component(i)%SumFraction )
        end do
      end if

      ! Calculation of residence times
      if ( this%ResidenceTime ) then
        call Reset( this%SumResidenceDuration )
        call Reset( this%SumResidencePairs )
      end if
      
      ! 2.) Combined sums
      call Reset( this%SumEPotSquared )
      call Reset( this%SumEPotV )
      call Reset( this%SumEPotVirial )
      call Reset( this%SumEnthalpySquared )
      call Reset( this%SumEnthalpyV )
      call Reset( this%SumVolumeSquared )
      call Reset( this%SumEPotCubic )
      call Reset( this%SumdEpotdVSquared )
      call Reset( this%SumEPotdEpotdV )
      call Reset( this%SumEPotSquareddEpotdV )
      call Reset( this%SumEPotdEpotdVSquared )
      call Reset( this%SumEPotd2EpotdV2 )
      if( EnsembleType .eq. EnsembleTypeNVE .and. LongRange .eq. Rfield) then
        call Reset( this%SumHmU )
        call Reset( this%SumHmUm1)
        call Reset( this%SumHmUm2 )
        call Reset( this%SumHmUm3 )
        call Reset( this%SumHmUm1dUdV )
        call Reset( this%SumHmUm1dUdV2 )
        call Reset( this%SumHmUm1d2UdV2 )
        call Reset( this%SumHmUm2dUdV )
        call Reset( this%SumHmUm2dUdV2 )
        call Reset( this%SumHmUm2d2UdV2 )
        call Reset( this%SumHmUm3dUdV )
        call Reset( this%SumHmUm3dUdV2 )
      end if

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
      if( LongRange .eq. Rfield) then
        if ( EnsembleType .eq. EnsembleTypeNVT ) then
          call Reset( this%SumA10resI )
          call Reset( this%SumA01resI )
          call Reset( this%SumA20resI )
          call Reset( this%SumA11resI )
          call Reset( this%SumA02resI )
          call Reset( this%SumA30resI )
          call Reset( this%SumA21resI )
          call Reset( this%SumA12resI )
        elseif ( EnsembleType .eq. EnsembleTypeNVE ) then
          call Reset( this%SumA10resI )
          call Reset( this%SumA01resI )
          call Reset( this%SumA20resI )
          call Reset( this%SumA11resI )
          call Reset( this%SumA02resI )
          call Reset( this%SumA30resI )
          call Reset( this%SumA21resI )
          call Reset( this%SumA12resI )
          call Reset( this%SumA10resII )
          call Reset( this%SumA01resII )
          call Reset( this%SumA20resII )
          call Reset( this%SumA11resII )
          call Reset( this%SumA02resII )
          call Reset( this%SumA30resII )
          call Reset( this%SumA21resII )
          call Reset( this%SumA12resII )
        end if
      end if

      ! 4.) Chemical potential and partial molar volumes
      do i = 1, this%NRealComponents
        select case( this%Component(i)%ChemPotMethod )
        case( ChemPotMethodGradIns )
          call Reset( this%Component(i)%SumInvChemPotRho )
          call Reset( this%Component(i)%SumInvChemPot )
        case( ChemPotMethodWidom )
          call Reset( this%Component(i)%SumChemPotV )
          call Reset( this%Component(i)%SumChemPotVV )
          call Reset( this%Component(i)%SumHW_counter )
          call Reset( this%Component(i)%SumHW_denom )
        case( ChemPotMethodThermoInt )
          call Reset( this%Component(i)%SumChemPotV )
          call Reset( this%Component(i)%SumChemPotThermoIntWidom )
          call Reset( this%Component(i)%SumChemPotThermoIntWidomV )
          call Reset( this%Component(i)%SumHW_counter )
          call Reset( this%Component(i)%SumHW_denom )
        end select
      end do

        do i = 1, this%NRealComponents
          if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodNone ) then
            call Reset( this%Component(i)%SumVW )
            call Reset( this%Component(i)%SumHM )
          end if
        end do

      ! Update result header
      if (SimulationType .eq. MonteCarlo) then
#if MPI_VER > 0
         fields = 0
         headers = headers + 1
         fields = fields + 7
         do i = 1, this%NRealComponents
           if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodNone ) then 
             fields = fields + 1
             if( EnsembleType .eq. EnsembleTypeNPT) then 
               fields = fields + 1
               if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodWidom ) fields = fields + 1
             end if
           end if
         enddo
         if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) fields = fields + this%NComponents + 1
#if CONSTR > 0
         fields = fields + 2 *  this%NCons
#endif

         if (RootProc) then
           if (CommonEqui) then
             offset = (accumulate_step/BlockSize+headers-1) * (11 * fields + 1) + headers-1 
           else
             offset = (NProcs * (accumulate_step/BlockSize)+headers-1) * (11 * fields + 1) + headers-1
           endif
           call MPI_File_Seek((this%iounit_result), offset, MPI_SEEK_SET, ierr)
           call MPI_File_Seek((this%iounit_runave), offset, MPI_SEEK_SET, ierr)
           write( IOBuffer, '(A)' )new_line('a')
           call FileWriteNoAdvance_parallel( this%iounit_result )
           call FileWriteNoAdvance_parallel( this%iounit_runave )
           ! PROC
           write( IOBuffer, '("       PROC")' )
           call FileWriteNoAdvance_parallel( this%iounit_result )
           call FileWriteNoAdvance_parallel( this%iounit_runave )

           ! Number of steps
           write( IOBuffer, '("         NR")' )
           call FileWriteNoAdvance_parallel( this%iounit_result )
           call FileWriteNoAdvance_parallel( this%iounit_runave )

           ! Pressure
           write( IOBuffer, '("      PRESS")' )
           call FileWriteNoAdvance_parallel( this%iounit_result )
           call FileWriteNoAdvance_parallel( this%iounit_runave )

           ! Density
           write( IOBuffer, '("    DENSITY")' )
           call FileWriteNoAdvance_parallel( this%iounit_result )
           call FileWriteNoAdvance_parallel( this%iounit_runave )

           ! Temperature
           write( IOBuffer, '("       TEMP")' )
           call FileWriteNoAdvance_parallel( this%iounit_result )
           call FileWriteNoAdvance_parallel( this%iounit_runave )

           ! Potential energy
           write( IOBuffer, '("         EPOT")' )
           call FileWriteNoAdvance_parallel( this%iounit_result )
           call FileWriteNoAdvance_parallel( this%iounit_runave )

           ! Enthalpy
           write( IOBuffer, '("        ENTLP")' )
           call FileWriteNoAdvance_parallel( this%iounit_result )
           call FileWriteNoAdvance_parallel( this%iounit_runave )

          if (printIDF) then
            ! Inter Potential energy
            write( IOBuffer, '("     EP_Inter")' )
            call FileWriteNoAdvance_parallel( this%iounit_result )
            call FileWriteNoAdvance_parallel( this%iounit_runave )

            ! Intra Potential energy
            write( IOBuffer, '("     EP_Intra")' )
            call FileWriteNoAdvance_parallel( this%iounit_result )
            call FileWriteNoAdvance_parallel( this%iounit_runave )

            ! Intra Potential energy - Bonds
            write( IOBuffer, '("     EP_Bonds")' )
            call FileWriteNoAdvance_parallel( this%iounit_result )
            call FileWriteNoAdvance_parallel( this%iounit_runave )

            ! Intra Potential energy - Angles
            write( IOBuffer, '("    EP_Angles")' )
            call FileWriteNoAdvance_parallel( this%iounit_result )
            call FileWriteNoAdvance_parallel( this%iounit_runave )

            ! Intra Potential energy - Dihedral
            write( IOBuffer, '("     EP_Dihed")' )
            call FileWriteNoAdvance_parallel( this%iounit_result )
            call FileWriteNoAdvance_parallel( this%iounit_runave )

            ! Intra Potential energy - Nonbonded
            write( IOBuffer, '("     EP_14_15")' )
            call FileWriteNoAdvance_parallel( this%iounit_result )
            call FileWriteNoAdvance_parallel( this%iounit_runave )

            ! Intra Virial
            write( IOBuffer, '("    Vir_Intra")' )
            call FileWriteNoAdvance_parallel( this%iounit_result )
            call FileWriteNoAdvance_parallel( this%iounit_runave )

            ! Inter Virial
            write( IOBuffer, '("      Vir_Inter")' )
            call FileWriteNoAdvance_parallel( this%iounit_result )
            call FileWriteNoAdvance_parallel( this%iounit_runave )
          end if

           ! Chemical potential
           do i = 1, this%NRealComponents
             if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodNone ) then
               if( i < 10 ) then
                 write( IOBuffer, '("      MUE_", I1)' ) i
               else
                 write( IOBuffer, '("     MUE_", I2)' ) i
               end if
               call FileWriteNoAdvance_parallel( this%iounit_result )
               call FileWriteNoAdvance_parallel( this%iounit_runave )
             end if
           end do

           ! Partial molar volume
           do i = 1, this%NRealComponents
             if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT ) then
               if( i < 10 ) then
                 write( IOBuffer, '("       VW_", I1)' ) i
               else
                 write( IOBuffer, '("      VW_", I2)' ) i
               end if
               call FileWriteNoAdvance_parallel( this%iounit_result )
               call FileWriteNoAdvance_parallel( this%iounit_runave )
             end if
           end do

           ! Partial molar enthalpy
           do i = 1, this%NRealComponents
             if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT ) then
               if( i < 10 ) then
                 write( IOBuffer, '("       HM_", I1)' ) i
               else
                 write( IOBuffer, '("      HM_", I2)' ) i
               end if
               call FileWriteNoAdvance_parallel( this%iounit_result )
               call FileWriteNoAdvance_parallel( this%iounit_runave )
             end if
           end do

           ! Number of particles in ensemble
           if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
             write( IOBuffer, '("      NPART")' )
             call FileWriteNoAdvance_parallel( this%iounit_result )
             call FileWriteNoAdvance_parallel( this%iounit_runave )

             ! Mole fraction of each component
             do i = 1, this%NComponents
               if( i < 10 ) then
                 write( IOBuffer, '("    FRACT_", I1)' ) i
               else
                 write( IOBuffer, '("   FRACT_", I2)' ) i
               end if
               call FileWriteNoAdvance_parallel( this%iounit_result )
               call FileWriteNoAdvance_parallel( this%iounit_runave )
             end do
           end if

#if CONSTR > 0
           do i=1, this%NCons
             if ( i < 10 ) then
               write( IOBuffer, '("      PMF_", I1)' ) i
               call FileWriteNoAdvance_parallel( this%iounit_result )
               call FileWriteNoAdvance_parallel( this%iounit_runave )
               write( IOBuffer, '("      MF_",  I1)' ) i
               call FileWriteNoAdvance_parallel( this%iounit_result )
               call FileWriteNoAdvance_parallel( this%iounit_runave )
             else
               write( IOBuffer, '("     PMF_", I2)' ) i
               call FileWriteNoAdvance_parallel( this%iounit_result )
               call FileWriteNoAdvance_parallel( this%iounit_runave )
               write( IOBuffer, '("     MF_",  I2)' ) i
               call FileWriteNoAdvance_parallel( this%iounit_result )
               call FileWriteNoAdvance_parallel( this%iounit_runave )
             end if
           end do
#endif  
           write( IOBuffer, '(A)' )new_line('a')
           call FileWriteNoAdvance_parallel( this%iounit_result )
           call FileWriteNoAdvance_parallel( this%iounit_runave )
         endif
#else
         call FileWriteBlank( this%iounit_result )
         call FileWriteBlank( this%iounit_runave )
         ! Number of steps
         write( IOBuffer, '("     NR")' )
         call FileWriteNoAdvance( this%iounit_result )
         call FileWriteNoAdvance( this%iounit_runave )

         ! Pressure
         write( IOBuffer, '("      PRESS")' )
         call FileWriteNoAdvance( this%iounit_result )
         call FileWriteNoAdvance( this%iounit_runave )

         ! Density
         write( IOBuffer, '("    DENSITY")' )
         call FileWriteNoAdvance( this%iounit_result )
         call FileWriteNoAdvance( this%iounit_runave )

         ! Temperature
         write( IOBuffer, '("       TEMP")' )
         call FileWriteNoAdvance( this%iounit_result )
         call FileWriteNoAdvance( this%iounit_runave )

         ! Potential energy
         write( IOBuffer, '("         EPOT")' )
         call FileWriteNoAdvance( this%iounit_result )
         call FileWriteNoAdvance( this%iounit_runave )

         ! Enthalpy
         write( IOBuffer, '("        ENTLP")' )
         call FileWriteNoAdvance( this%iounit_result )
         call FileWriteNoAdvance( this%iounit_runave )

        if (printIDF) then
          ! Inter Potential energy
          write( IOBuffer, '("     EP_Inter")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Intra Potential energy
          write( IOBuffer, '("     EP_Intra")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Intra Potential energy - Bonds
          write( IOBuffer, '("     EP_Bonds")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Intra Potential energy - Angles
          write( IOBuffer, '("    EP_Angles")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Intra Potential energy - Dihedral
          write( IOBuffer, '("     EP_Dihed")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Intra Potential energy - Nonbonded
          write( IOBuffer, '("     EP_14_15")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Intra Virial
          write( IOBuffer, '("    Vir_Intra")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Inter Virial
          write( IOBuffer, '("      Vir_Inter")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )
        end if

         ! Chemical potential
         do i = 1, this%NRealComponents
           if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodNone ) then
             if( i < 10 ) then
               write( IOBuffer, '("      MUE_", I1)' ) i
             else
               write( IOBuffer, '("     MUE_", I2)' ) i
             end if
             call FileWriteNoAdvance( this%iounit_result )
             call FileWriteNoAdvance( this%iounit_runave )
           end if
         end do

         ! Partial molar volume
         do i = 1, this%NRealComponents
           if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT) then
             if( i < 10 ) then
               write( IOBuffer, '("       VW_", I1)' ) i
             else
               write( IOBuffer, '("      VW_", I2)' ) i
             end if
             call FileWriteNoAdvance( this%iounit_result )
             call FileWriteNoAdvance( this%iounit_runave )
           end if
         end do

         ! Partial molar enthalpy
         do i = 1, this%NRealComponents
           if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT) then
             if( i < 10 ) then
               write( IOBuffer, '("       HM_", I1)' ) i
             else
               write( IOBuffer, '("      HM_", I2)' ) i
             end if
             call FileWriteNoAdvance( this%iounit_result )
             call FileWriteNoAdvance( this%iounit_runave )
           end if
         end do

         ! Number of particles in ensemble
         if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
           write( IOBuffer, '("      NPART")' )
           call FileWriteNoAdvance( this%iounit_result )
           call FileWriteNoAdvance( this%iounit_runave )

           ! Mole fraction of each component
           do i = 1, this%NComponents
             if( i < 10 ) then
               write( IOBuffer, '("    FRACT_", I1)' ) i
             else
               write( IOBuffer, '("   FRACT_", I2)' ) i
             end if
             call FileWriteNoAdvance( this%iounit_result )
             call FileWriteNoAdvance( this%iounit_runave )
           end do
         end if

#if CONSTR > 0
         do i=1, this%NCons
           if ( i < 10 ) then
             write( IOBuffer, '("      PMF_", I1)' ) i
             call FileWriteNoAdvance( this%iounit_runave )
             write( IOBuffer, '("      MF_",  I1)' ) i
             call FileWriteNoAdvance( this%iounit_runave )
           else
             write( IOBuffer, '("     PMF_", I2)' ) i
             call FileWriteNoAdvance( this%iounit_runave )
             write( IOBuffer, '("     MF_",  I2)' ) i
             call FileWriteNoAdvance( this%iounit_runave )
           end if
         end do
#endif  
        call FileWriteBlank( this%iounit_result )
        call FileWriteBlank( this%iounit_runave )
#endif
      else !MD
        call FileWriteBlank( this%iounit_result )
        call FileWriteBlank( this%iounit_runave )
        ! Number of steps
        write( IOBuffer, '("     NR")' )
        call FileWriteNoAdvance( this%iounit_result )
        call FileWriteNoAdvance( this%iounit_runave )

        ! Displacement
          write( IOBuffer, '("     DISP")' )
          call FileWriteNoAdvance( this%iounit_runave )

        ! Pressure
        write( IOBuffer, '("      PRESS")' )
        call FileWriteNoAdvance( this%iounit_result )
        call FileWriteNoAdvance( this%iounit_runave )

        ! Density
        write( IOBuffer, '("    DENSITY")' )
        call FileWriteNoAdvance( this%iounit_result )
        call FileWriteNoAdvance( this%iounit_runave )

        ! Temperature
        write( IOBuffer, '("       TEMP")' )
        call FileWriteNoAdvance( this%iounit_result )
        call FileWriteNoAdvance( this%iounit_runave )

#if OSMOP > 0
           ! OsmoticPressure
           write( IOBuffer, '("      OSPR")' )
           call FileWriteNoAdvance( this%iounit_result )
           call FileWriteNoAdvance( this%iounit_runave )
#endif

        ! Potential energy
        write( IOBuffer, '("         EPOT")' )
        call FileWriteNoAdvance( this%iounit_result )
        call FileWriteNoAdvance( this%iounit_runave )

        ! Enthalpy
        write( IOBuffer, '("        ENTLP")' )
        call FileWriteNoAdvance( this%iounit_result )
        call FileWriteNoAdvance( this%iounit_runave )

        if (printIDF) then
          ! Inter Potential energy
          write( IOBuffer, '("     EP_Inter")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Intra Potential energy
          write( IOBuffer, '("     EP_Intra")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Intra Potential energy - Bonds
          write( IOBuffer, '("     EP_Bonds")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Intra Potential energy - Angles
          write( IOBuffer, '("    EP_Angles")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Intra Potential energy - Dihedral
          write( IOBuffer, '("     EP_Dihed")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Intra Potential energy - Nonbonded
          write( IOBuffer, '("     EP_14_15")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Intra Virial
          write( IOBuffer, '("    Vir_Intra")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Inter Virial
          write( IOBuffer, '("      Vir_Inter")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )
        end if

        ! Chemical potential
        do i = 1, this%NRealComponents
          if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodNone ) then
            if( i < 10 ) then
              write( IOBuffer, '("       MUE_", I1)' ) i
            else
              write( IOBuffer, '("      MUE_", I2)' ) i
            end if
            call FileWriteNoAdvance( this%iounit_result )
            call FileWriteNoAdvance( this%iounit_runave )
          end if
        end do

        ! Partial molar volume
        do i = 1, this%NRealComponents
          if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT) then
            if( i < 10 ) then
              write( IOBuffer, '("        VW_", I1)' ) i
            else
              write( IOBuffer, '("       VW_", I2)' ) i
            end if
            call FileWriteNoAdvance( this%iounit_result )
            call FileWriteNoAdvance( this%iounit_runave )
          end if
        end do

        ! Partial molar enthalpy
        do i = 1, this%NRealComponents
          if( this%Component(i)%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT) then
            if( i < 10 ) then
              write( IOBuffer, '("        HM_", I1)' ) i
            else
              write( IOBuffer, '("       HM_", I2)' ) i
            end if
            call FileWriteNoAdvance( this%iounit_result )
            call FileWriteNoAdvance( this%iounit_runave )
          end if
        end do

#if HBOND > 0
        do i = 1, this%NComponents
          write( IOBuffer, '("  HB0_(", I1, ")")' ) i
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )
        end do
        do i = 1, this%NComponents
          do  j = 1, this%NComponents
            write( IOBuffer, '("  HB1_(", I1, ",", I1, ")")' ) i, j
            call FileWriteNoAdvance( this%iounit_result )
            call FileWriteNoAdvance( this%iounit_runave )
          end do
        end do
        do i = 1, this%NComponents
          do  j = 1, this%NComponents
            do k = j, this%NComponents
              write( IOBuffer, '("  HB2_(", I1, ",", I1, ",", I1, ")")' ) i, j, k
              call FileWriteNoAdvance( this%iounit_result )
              call FileWriteNoAdvance( this%iounit_runave )
            end do
          end do
        end do
        do i = 1, this%NComponents
          do  j = 1, this%NComponents
            do k = j, this%NComponents
              do  l = k, this%NComponents
                write( IOBuffer, '("  HB3_(", I1, ",", I1, ",", I1, ",", I1, ")")' ) i, j, k, l
                call FileWriteNoAdvance( this%iounit_result )
                call FileWriteNoAdvance( this%iounit_runave )
              end do
            end do
          end do
        end do
        do i = 1, this%NComponents
          write( IOBuffer, '("  HB4+_(", I1, ")")' ) i
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )
        end do
#endif

#if OSMOP > 0
        !Density Profile
        do i = 1, this%NComponents
          do j = 1, NBinsDen
            if (j .le. 9) then 
              write( IOBuffer, '("   DP", I1, "B00", I1)' ) i, j
            elseif (j .le. 99) then 
              write( IOBuffer, '("    DP", I1, "B0", I2)' ) i, j
            else
              write( IOBuffer, '("     DP", I1, "B", I3)' ) i, j
            endif
            call FileWriteNoAdvance( this%iounit_result )
            call FileWriteNoAdvance( this%iounit_runave )
          end do
        end do

#if OSMOP == 2
        !Pressure Profile
        do j = 1, NBinsDen
          if (j .le. 9) then 
            write( IOBuffer, '(" PPB00", I1)' ) j
          elseif (j .le. 99) then 
            write( IOBuffer, '(" PPB0", I2)' ) j
          else
            write( IOBuffer, '(" PPB", I3)' ) j
          endif
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )
        end do

        !Chemical Potential Profile
        do i = 1, this%NRealComponents
          if( this%Component(i)%ChemPotMethod .eq. ChemPotMethodWidom ) then
            do j = 1, NBinsDen
              if (j .le. 9) then 
                write( IOBuffer, '("     CP", I1, "B00", I1)' ) i, j
              elseif (j .le. 99) then 
                write( IOBuffer, '("     CP", I1, "B0", I2)' ) i, j
              else
                write( IOBuffer, '("     CP", I1, "B", I3)' ) i, j
              endif
              call FileWriteNoAdvance( this%iounit_result )
              call FileWriteNoAdvance( this%iounit_runave )
            end do
          end if
        end do
#endif
#endif

        ! Number of particles in ensemble
        if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
          write( IOBuffer, '("      NPART")' )
          call FileWriteNoAdvance( this%iounit_result )
          call FileWriteNoAdvance( this%iounit_runave )

          ! Mole fraction of each component
          do i = 1, this%NComponents
            if( i < 10 ) then
              write( IOBuffer, '("     FRACT_", I1)' ) i
            else
              write( IOBuffer, '("    FRACT_", I2)' ) i
            end if
            call FileWriteNoAdvance( this%iounit_result )
            call FileWriteNoAdvance( this%iounit_runave )
          end do
        end if

#if CONSTR > 0
        do i=1, this%NCons
          if ( i < 10 ) then
            write( IOBuffer, '("      PMF_", I1)' ) i
            call FileWriteNoAdvance( this%iounit_runave )
            write( IOBuffer, '("      MF_",  I1)' ) i
            call FileWriteNoAdvance( this%iounit_runave )
          else
            write( IOBuffer, '("     PMF_", I2)' ) i
            call FileWriteNoAdvance( this%iounit_runave )
            write( IOBuffer, '("     MF_",  I2)' ) i
            call FileWriteNoAdvance( this%iounit_runave )
          end if
        end do
#endif   
        call FileWriteBlank( this%iounit_result )
        call FileWriteBlank( this%iounit_runave )
      end if
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END IF of step ==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Update accumulators
    ! 1.) Basic sums
    call Update( this%SumPressure, this%Pressure )
    call Update( this%SumDensity, this%Density )
    call Update( this%SumTemperature, this%Temperature )
    call Update( this%SumEPot, this%EPot / real( this%NPart, RK ) )
    call Update( this%SumVolume, 1._RK / this%Density )
    call Update( this%SumVirial, -3._RK * this%Virial )
    call Update( this%SumEPotInter, this%EPotInter / real( this%NPart, RK ) )
    call Update( this%SumEPotIntra, this%EPotIntra / real( this%NPart, RK ) )
    if (printIDF) then
      call Update( this%SumEPotIntra_Bond, this%EPotIntra_Bond / real( this%NPart, RK ) )
      call Update( this%SumEPotIntra_Angle, this%EPotIntra_Angle / real( this%NPart, RK ) )
      call Update( this%SumEPotIntra_Dihedral, this%EPotIntra_Dihedral / real( this%NPart, RK ) )
      call Update( this%SumEPotIntra_Nonbonded, this%EPotIntra_Nonbonded / real( this%NPart, RK ) )
      call Update( this%SumVirialIntra, -3._RK * this%VirialIntra )
      call Update( this%SumVirialInter, -3._RK * this%VirialInter )
    end if
#if OSMOP > 0
    call Update( this%SumOsmoticPressure, this%OsmoticPressure )
    do i = 1, this%NComponents
      do j = 1, NBinsDen
        call Update( this%Component(i)%SumDenProfile(j), real(this%Component(i)%DensityProfileN(j),RK)*NBinsDen/this%Volume0)
#if OSMOP == 2
        if( this%Component(i)%ChemPotMethod .eq. ChemPotMethodWidom ) then
          call Update( this%Component(i)%SumChemPotProfile(j), this%Component(i)%ChemPotProfile(j))
        end if
#endif
      end do
    end do
#if OSMOP == 2
    do j = 1, NBinsDen
      call Update( this%SumPressureProfile(j), this%PressureProfile(j))
    end do
#endif
#endif

#if HBOND > 0
    do i = 1, this%NComponents
      call Update( this%SumHBond0(i), real(this%NHBond0(i),RK) )
      do j = 1, this%NComponents
        call Update( this%SumHBond1(i,j), real(this%NHBond1(i,j),RK) )
        do k = j, this%NComponents
          call Update( this%SumHBond2(i,j,k), real(this%NHBond2(i,j,k),RK) )
          do l = k, this%NComponents
            call Update( this%SumHBond3(i,j,k,l), real(this%NHBond3(i,j,k,l),RK) )
          end do
        end do
      end do
      call Update( this%SumHBondN(i), real(this%NHBondN(i),RK) )
    end do
#endif
    if( ConstantPressure ) then
      call Update( this%SumEnthalpy, this%EPotInter / real( this%NPart, RK ) + this%RefPressure / this%Density - &
&      (1-this%NUnitTotal/this%Npart)*this%RefTemperature )
!       call Update( this%SumEnthalpy, this%EPot/real(this%NPart,RK) + this%Pressure/this%Density - this%RefTemperature) - refT to adjust H=U+pv with p_res, u already u_res
    else
      call Update( this%SumEnthalpy, this%EPot/real(this%NPart,RK) + this%Pressure/this%Density)
    end if

    call Update( this%SumVirialIntra, -3._RK * this%VirialIntra )
    call Update( this%SumVirialInter, -3._RK * this%VirialInter )

    currentdEpotdV   = -(this%Virial+(this%NUnitTotal-this%Npart)*this%RefTemperature)/this%Volume0
    currentd2EpotdV2 =  ((2._RK*this%Virial/3._RK + this%d2EpotdV2) + (this%NUnitTotal-this%Npart)*this%RefTemperature)/this%Volume0**2 ! diff to trunk...wrong! GABOR!!!
    call Update( this%SumdEpotdV,   currentdEpotdV)
    call Update( this%Sumd2EpotdV2, currentd2EpotdV2)

    if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
      call Update( this%SumNPart, real( this%NPart, RK ) )

      do i = 1, this%NComponents
        pc => this%Component(i)
        call Update( pc%SumFraction, pc%Fraction )
      end do
    end if

    ! 2.) Combined sums
    call Update( this%SumEPotSquared,      ( this%EPotInter / real( this%NPart, RK ) )**2 ) ! diff to trunk all 7 lines
    call Update( this%SumEPotCubic,          this%EPotInter**3 )
    call Update( this%SumdEpotdVSquared,                    currentdEpotdV**2 )
    call Update( this%SumEPotdEpotdV,        this%EPotInter    * currentdEpotdV    )             
    call Update( this%SumEPotSquareddEpotdV, this%EPotInter**2 * currentdEpotdV    )
    call Update( this%SumEPotdEpotdVSquared, this%EPotInter    * currentdEpotdV**2 )
    call Update( this%SumEPotd2EpotdV2,      this%EPotInter    * currentd2EpotdV2  )

    if( EnsembleType .eq. EnsembleTypeNVE .and. LongRange .eq. Rfield ) then
      !Following was part was commented, even if J.Chem.Phys.100(4)1994 prescribes it for NVEMom MD, because the results are identical with and without it.
      !if( SimulationType .eq. MolecularDynamics ) then  
      !  Momentum(:) = 0._RK
      !  Mass = 0._RK
      !  do j = 1, this%NComponents
      !    Mass=Mass+this%Component(j)%Molecule%Mass*real(this%Component(j)%NPart, RK)
      !    do i = 1, 3
      !      Momentum(i)=Momentum(i)+this%Component(j)%Molecule%Mass*sum(this%Component(j)%P1(1:this%Component(j)%NPart,i))
      !    end do
      !  end do
      !  Momentumd2Mass=(Momentum(1)*Momentum(1)+Momentum(2)*Momentum(2)+Momentum(3)*Momentum(3))/(2._RK*Mass)
      !  currentHmU = (this%RefHamiltonian*real( this%NPart, RK ) - this%EPot) - Momentumd2Mass
      !else
        currentHmU = (this%RefHamiltonian*real( this%NPart, RK ) - this%EPotInter)
      !endif
      currentHmUm1 = 1._RK/currentHmU
      call Update( this%SumHmU,            currentHmU )
      call Update( this%SumHmUm1,          currentHmUm1 )
      call Update( this%SumHmUm2,          currentHmUm1**2 )
      call Update( this%SumHmUm3,          currentHmUm1**3 )
      call Update( this%SumHmUm1dUdV,      currentHmUm1     * currentdEpotdV    )
      call Update( this%SumHmUm1dUdV2,     currentHmUm1     * currentdEpotdV**2 )
      call Update( this%SumHmUm1d2UdV2,    currentHmUm1     * currentd2EpotdV2  )
      call Update( this%SumHmUm2dUdV,      currentHmUm1**2  * currentdEpotdV    )
      call Update( this%SumHmUm2dUdV2,     currentHmUm1**2  * currentdEpotdV**2 )
      call Update( this%SumHmUm2d2UdV2,    currentHmUm1**2  * currentd2EpotdV2  )
      call Update( this%SumHmUm3dUdV,      currentHmUm1**3  * currentdEpotdV    )
      call Update( this%SumHmUm3dUdV2,     currentHmUm1**3  * currentdEpotdV**2 )
    endif

    call Update( this%SumEPotV, this%EPotInter / this%Volume0  )

    call Update( this%SumEPotVirial, -3. * this%Virial * this%EPot / real( this%NPart, RK ) )

    if( ConstantPressure ) then
       call Update( this%SumEnthalpySquared, ( this%EPotInter / real( this%NPart, RK ) + &
&                this%RefPressure / this%Density - (1-this%NUnitTotal/this%Npart)*this%RefTemperature )**2 )
   
       call Update( this%SumEnthalpyV, ( this%EPotInter / real( this%NPart, RK ) + &
&                this%RefPressure / this%Density - (1-this%NUnitTotal/this%Npart)*this%RefTemperature ) / this%Density )
    else
       call Update( this%SumEnthalpySquared, ( this%EPot / real( this%NPart, RK ) + &
&                this%Pressure / this%Density )**2 )
   
       call Update( this%SumEnthalpyV, ( this%EPotInter / real( this%NPart, RK ) + &
&                this%Pressure / this%Density  ) / this%Density )
    end if
    call Update( this%SumVolumeSquared, 1._RK / this%Density**2 )

    ! 3.) Derived sums
    if( ConstantPressure ) then

      call Update( this%SumBetaT, real( this%NPart, RK ) / this%RefTemperature &
&                * ( this%SumVolumeSquared%Average / this%SumVolume%Average - this%SumVolume%Average ) )

      call Update( this%SumdHdP, this%SumVolume%Average - real( this%NPart, RK ) / this%RefTemperature &
&                * ( this%SumEPotV%Average - this%SumEPotInter%Average * this%SumVolume%Average + this%RefPressure &
&                * ( this%SumVolumeSquared%Average - this%SumVolume%Average**2 ) ) )

      call Update( this%SumCP, real( this%NPart, RK ) / this%RefTemperature**2 &
&                * ( this%SumEnthalpySquared%Average - this%SumEnthalpy%Average**2 ) )

      call Update( this%SumAlphaP, real( this%NPart, RK ) / this%RefTemperature**2 &
&                * this%SumDensity%Average * ( this%SumEnthalpyV%Average &
&                - this%SumEnthalpy%Average * this%SumVolume%Average ) )

    else

      call Update( this%SumdUdV, (this%NPart / this%RefTemperature &
&                * ( this%SumVirial%Average * this%SumEPot%Average - this%SumEPotVirial%Average )&
&                + this%SumVirial%Average ) / 3._RK / this%Volume0 )

      call Update( this%SumCV, real( this%NPart, RK ) / this%RefTemperature**2 &
&                * ( this%SumEPotSquared%Average - this%SumEPot%Average**2 ) )
    endif

    if( EnsembleType .eq. EnsembleTypeNVT .and. LongRange .eq. Rfield ) then
      Beta    = 1._RK/this%RefTemperature
      Beta2   = Beta*Beta
      Beta3   = Beta*Beta2
      specv   = this%Volume0/this%NUnitTotal
      specv2  = specv*specv
      Numb    = real( this%NUnitTotal, RK )
      U       = this%SumEpot%Average*real( this%NPart, RK )
      U2      = this%SumEpotSquared%Average*real( this%NPart, RK )**2
      U3      = this%SumEpotCubic%Average
      dUdV    = this%SumdEpotdV%Average
      dUdV2   = this%SumdEpotdVSquared%Average
      UdUdV   = this%SumEpotdEpotdV%Average
      U2dUdV  = this%SumEPotSquareddEpotdV%Average
      UdUdV2  = this%SumEPotdEpotdVSquared%Average
      d2UdV2  = this%Sumd2EpotdV2%Average
      Ud2UdV2 = this%SumEPotd2EpotdV2%Average

      A10res =  Beta*U/Numb
      A01res =  (Numb-this%constrNDF/3._RK)/Numb - Beta*specv*dUdV
      A20res =  Beta2*(U*U-U2)/Numb
      A11res =  specv*(-Beta*dUdV + Beta2*UdUdV - Beta2*U*dUdV)
      A02res =  Numb*specv2*(Beta*d2UdV2 - Beta2*dUdV2 + Beta2*dUdV**2) + 2._RK*specv*Beta*dUdV - (Numb-this%constrNDF/3._RK)/Numb
      A30res =  Beta3*(U3 -3._RK*U*U2 + 2._RK*U**3)/Numb
      A21res =  specv*( Beta2*( 2._RK*UdUdV - 2._RK*U*dUdV) + Beta3*(U2*dUdV - U2dUdV + 2._RK*U*UdUdV - 2._RK*U**2*dUdV) )
      A12res =  Numb*specv2*Beta3*( UdUdV2 + 2._RK*U*dUdV**2 - U*dUdV2 - 2._RK*UdUdV*dUdV)+&
&               Numb*specv2*Beta2*( 2._RK*dUdV**2 + U*d2UdV2 - 2._RK*dUdV2 - Ud2UdV2)+&
&               Numb*specv2*Beta*d2UdV2+&
&               specv*Beta2*(2._RK*U*dUdV - 2._RK*UdUdV)+&
&               specv*Beta*2._RK*dUdV

      call Update( this%SumA10resI, A10res )
      call Update( this%SumA01resI, A01res )
      call Update( this%SumA20resI, A20res )
      call Update( this%SumA11resI, A11res )
      call Update( this%SumA02resI, A02res )
      call Update( this%SumA30resI, A30res )
      call Update( this%SumA21resI, A21res )
      call Update( this%SumA12resI, A12res )
    end if

    if( EnsembleType .eq. EnsembleTypeNVE .and. LongRange .eq. Rfield) then
      specv       = 1._RK/this%Density
      specv2      = specv*specv
      Numb        = real( this%NPart, RK )
      rho         = this%Density
      rho2        = rho*rho
      !Following was part was commented, even if J.Chem.Phys.100(4)1994 prescribes it for NVEMom MD, because it deteriorates the results considerably.
      !if( SimulationType .eq. MolecularDynamics ) then 
      ! F = (real(this%NDF, RK)- 3._RK)/2._RK
      !else
       F = real(this%NDF-this%constrNDF, RK)/2._RK
      !endif
      invF        = 2._RK/real(this%NDF-this%constrNDF, RK)
      dUdV        = this%SumdEpotdV%Average
      dUdV2       = this%SumdEpotdVSquared%Average
      d2UdV2      = this%Sumd2EpotdV2%Average
      HmU         =this%SumHmU%Average
      HmUm1       =this%SumHmUm1%Average
      HmUm2       =this%SumHmUm2%Average
      HmUm3       =this%SumHmUm3%Average
      HmUm1dUdV   =this%SumHmUm1dUdV%Average
      HmUm1dUdV2  =this%SumHmUm1dUdV2%Average
      HmUm1d2UdV2 =this%SumHmUm1d2UdV2%Average
      HmUm2dUdV   =this%SumHmUm2dUdV%Average
      HmUm2dUdV2  =this%SumHmUm2dUdV2%Average
      HmUm2d2UdV2 =this%SumHmUm2d2UdV2%Average
      HmUm3dUdV   =this%SumHmUm3dUdV%Average
      HmUm3dUdV2  =this%SumHmUm3dUdV2%Average

      O00 = HmU*invF
      O10 = 1._RK
      O01 = rho*HmU*invF - dUdV
      O20 = (F-1._RK)*HmUm1
      O11 = rho + (1._RK-F)*HmUm1dUdV
      O02 = rho*rho*(1._RK-1._RK/Numb)*HmU*invF - 2._RK*rho*dUdV - d2UdV2 + (F-1._RK)*HmUm1dUdV2
      O30 = (F*F-3._RK*F+2._RK)*HmUm2
      O21 = rho*(F-1._RK)*HmUm1+(3._RK*F-F*F-2._RK)*HmUm2dUdV
      O12 = rho*(rho*(1._RK-1._RK/Numb) + 2._RK*HmUm1dUdV*(1._RK-F)) + (1._RK-F)*HmUm1d2UdV2 + (2._RK-3._RK*F+F*F)*HmUm2dUdV2
      O40 = (F*F*F-6._RK*F*F+11._RK*F-6._RK)*HmUm3
      O31 = rho*(F*F-3._RK*F+2._RK)*HmUm2 + (6._RK-11._RK*F+6._RK*F*F-F*F*F)*HmUm3dUdV
      O22 = rho*rho*(F-1._RK-F/Numb+1._RK/Numb)*HmUm1 + rho*(6._RK*F-2._RK*F*F-4._RK)*HmUm2dUdV + (3._RK*F-F*F-2._RK)*HmUm2d2UdV2 + (F*F*F-6._RK*F*F+11._RK*F-6._RK)*HmUm3dUdV2

      O00m1 = 1._RK/O00
      O00m2 = O00m1*O00m1
      O00m3 = O00m2*O00m1
      O012  = O01*O01

      !Entropy definition I

      S01 =  O00m1*O01
      S10 =  O00m1
      S02 = -O00m2*O012+O00m1*O02
      S11 = -O00m2*O01 +O00m1*O11
      S20 = -O00m2     +O00m1*O20
      S12 = 2._RK*O00m3*O012-2._RK*O00m2*O01*O11      -O00m2*O02+O00m1*O12
      S21 = 2._RK*O00m3*O01       -O00m2*O01*O20-2._RK*O00m2*O11+O00m1*O21
      S30 = 2._RK*O00m3     -3._RK*O00m2*O20          +O00m1*O30

      S20m1 = 1._RK/S20
      S20m2 = S20m1*S20m1
      S20m3 = S20m2*S20m1

      A01res =  S01
      A10res = -this%RefHamiltonian*Numb
      A02res =  S02-S11*S11*S20m1
      A11res =  S11*S20m1
      A20res = -S20m1
      A12res =  S20m1*(S11*S20m1*(S11*S20m1*S30-2._RK*S21)+S12)
      A21res = -S20m2*(S11*S20m1*S30-S21)
      A30res =  S20m3*S30

      ! Substraction of the ideal part

      Beta    = O00m1
      Beta2   = Beta*Beta
      Beta3   = Beta*Beta2

      A10res =  A10res-(      -F*O00   ) !Beta=O00m1
      A20res =  A20res-(       F*O00**2)
      A30res =  A30res-(-2._RK*F*O00**3)
      A01res =  A01res-  rho
      A02res =  A02res-(-rho2/Numb)

      ! Final conversion

      call Update( this%SumA10resI, -A10res*Beta/Numb )
      call Update( this%SumA01resI,  A01res*specv ) !=-(-V*A01res)/Numb
      call Update( this%SumA20resI, -A20res*Beta2/Numb )
      call Update( this%SumA11resI,  A11res*Beta*specv ) !=-(-V*Beta*A11res)/Numb
      call Update( this%SumA02resI, -(specv2*A02res*Numb + 2._RK*specv*A01res) ) !=-(V^2*A02res+2V*A01res)/Numb
      call Update( this%SumA30resI, -A30res*Beta3/Numb )
      call Update( this%SumA21resI,  A21res*specv*Beta2 ) !=-(-V*Beta2*A21res)/Numb
      call Update( this%SumA12resI, -Beta*(specv2*A12res*Numb + 2._RK*specv*A11res) ) !=-Beta*(V^2*A12res+2V*A11res)/Numb

      !Entropy definition II

      S01 =  O11
      S10 =  O20
      S02 = -O11*O11+O12
      S11 = -O11*O20+O21
      S20 = -O20*O20+O30
      S12 = 2._RK*O11*O11*O20-2._RK*O11*O21      -O12*O20+O22
      S21 = 2._RK*O11*O20*O20      -O11*O30-2._RK*O20*O21+O31
      S30 = 2._RK*O20*O20*O20-3._RK*O20*O30      +O40

      S20m1 = 1._RK/S20
      S20m2 = S20m1*S20m1
      S20m3 = S20m2*S20m1

      A01res =  S01
      A10res = -this%RefHamiltonian*Numb
      A02res =  S02-S11*S11*S20m1
      A11res =  S11*S20m1
      A20res = -S20m1
      A12res =  S20m1*(S11*S20m1*(S11*S20m1*S30-2._RK*S21)+S12)
      A21res = -S20m2*(S11*S20m1*S30-S21)
      A30res =  S20m3*S30

      ! Substraction of the ideal part

      O20m1   = 1._RK/O20
      Beta    = O20
      Beta2   = Beta*Beta
      Beta3   = Beta*Beta2

      funcF=(1._RK-invF)

      A10res =  A10res-(      -F* funcF*O20m1    )       !Beta=O20
      A20res =  A20res-(       F*(funcF*O20m1)**2)/funcF
      A30res =  A30res-(-2._RK*F*(funcF*O20m1)**3)/funcF**2
      A01res =  A01res-rho
      A02res =  A02res-(-rho2/Numb)

      ! Final conversion

      call Update( this%SumA10resII, -A10res*Beta/Numb )
      call Update( this%SumA01resII,  A01res*specv ) !=-(-V*A01res)/Numb
      call Update( this%SumA20resII, -A20res*Beta2/Numb )
      call Update( this%SumA11resII,  A11res*Beta*specv ) !=-(-V*Beta*A11res)/Numb
      call Update( this%SumA02resII, -(specv2*A02res*Numb + 2._RK*specv*A01res) ) !=-(V^2*A02res+2V*A01res)/Numb
      call Update( this%SumA30resII, -A30res*Beta3/Numb )
      call Update( this%SumA21resII,  A21res*specv*Beta2 ) !=-(-V*Beta2*A21res)/Numb
      call Update( this%SumA12resII, -Beta*(specv2*A12res*Numb + 2._RK*specv*A11res) ) !=-Beta*(V^2*A12res+2V*A11res)/Numb
    end if

#if  TRANS == 1
    ! 4.) Tranport properties !TRANSPORT_start
    if( mod((Step-1),this%NStepCorr) .eq. 0 ) then ! Michael Sch.: this if needed?

      if( mod( (Step-1)/this%NStepCorr-this%NCorr+1, BlockSizeCF*this%NSpanCF ) == 0 .and. (this%Mmess > 0) ) then

        do i = 1, this%NComponents
          call Update( this%Sumself_i(i), this%selfd_i(i), this%Mmess )
        end do

        if(this%NComponents .gt. 1) then
          if(this%NComponents == 2) then
            call Update( this%SumSoret, this%soret, this%Mmess )
          end if
          do i = 1, this%NComponents
            do j = 1, this%NComponents
              call Update( this%SumOnsager(i,j),this%Onsager(i,j), this%Mmess )
            end do
          end do
        end if

        call Update( this%SumVisco_s, this%visco_s, this%Mmess )
        call Update( this%SumVisco_b, this%visco_b, this%Mmess )
        call Update( this%SumConduct, this%conduct, this%Mmess )
        call Update( this%SumEConduct, this%econduct, this%Mmess )
      end if
    end if
!TRANSPORT_END
#endif

    t = this%NRealComponents+1  ! pseudo component identifier for ThermoInt (ThermoInt does not function together with GradIns)

    ! 4.) Chemical potential and partial molar volumes
    do i = 1, this%NRealComponents
      pc => this%Component(i)
      if( pc%CalcChemPot ) then
        select case( pc%ChemPotMethod )
        case( ChemPotMethodGradins )
#if MPI_VER > 0
          ! Per Process we calculate GI only for one component 
          if (mod(NProc,this%NGradInsComp)/=pc%NGradThis) cycle
!           if (mod(NProc,this%NRealComponents)/=i-1) cycle
#endif
          call Update( pc%SumInvChemPotRho, 1._RK / pc%ChemPot )
          call Update( pc%SumInvChemPot, 1._RK / pc%ChemPot1 )
        case( ChemPotMethodWidom )
          call Update( pc%SumChemPotV, pc%ChemPot / this%Density )
          call Update( pc%SumChemPotVV, pc%ChemPot / this%Density**2 )
          call Update(pc%SumHW_counter, pc%HW_counter)
          call Update(pc%SumHW_denom, pc%HW_denom)
        case( ChemPotMethodThermoInt )
          if( .not. Equilibration .and. mod(Step,pc%changeLaFreq) .ge. pc%forfeitLaSampl ) then
            currentbin=int((this%Component(t)%Lambda-pc%LaMin)/pc%deltaLa)
            pc%BinsVisit(currentbin)=pc%BinsVisit(currentbin)+1

            currentH=this%EPot + this%RefPressure * real( this%NPart, RK ) / this%Density
            pc%BinsEn(currentbin)     =  (                  pc%currentBinsEn                                       + (pc%BinsVisit(currentbin)-1)*pc%BinsEn(currentbin)    )/pc%BinsVisit(currentbin)
            pc%BinsdEndLa(currentbin) =  (pc%LambdaExponent*pc%currentBinsEn/this%Component(t)%Lambda              + (pc%BinsVisit(currentbin)-1)*pc%BinsdEndLa(currentbin))/pc%BinsVisit(currentbin)
            pc%BinsdEndLaV(currentbin) = (pc%LambdaExponent*pc%currentBinsEn/this%Component(t)%Lambda/this%Density + (pc%BinsVisit(currentbin)-1)*pc%BinsdEndLaV(currentbin))/pc%BinsVisit(currentbin)
            pc%BinsdEndLaH(currentbin)=(pc%LambdaExponent*pc%currentBinsEn/this%Component(t)%Lambda*currentH     + (pc%BinsVisit(currentbin)-1)*pc%BinsdEndLaH(currentbin))/pc%BinsVisit(currentbin)

            pc%BinsIntdEndLa(0)=pc%BinsdEndLa(0)*pc%deltaLa
            pc%BinsIntVW(0)=(pc%BinsdEndLaV(0)-pc%BinsdEndLa(0)*this%SumVolume%Average)*pc%deltaLa
            pc%BinsIntHW(0)=(pc%BinsdEndLaH(0)-pc%BinsdEndLa(0)*this%SumConfEnthalpy%Average+pc%BinsdEndLa(0)/this%RefTemperature)*pc%deltaLa
            do j = 1, pc%NBins-1
              pc%BinsIntdEndLa(j)=pc%BinsIntdEndLa(j-1)+pc%BinsdEndLa(j)*pc%deltaLa
              pc%BinsIntVW(j)=pc%BinsIntVW(j-1)+(pc%BinsdEndLaV(j)-pc%BinsdEndLa(j)*this%SumVolume%Average)*pc%deltaLa
              pc%BinsIntHW(j)=pc%BinsIntHW(j-1)+(pc%BinsdEndLaH(j)-pc%BinsdEndLa(j)*this%SumConfEnthalpy%Average+pc%BinsdEndLa(j)/this%RefTemperature)*pc%deltaLa
            end do
          end if

            call Update( pc%SumChemPotThermoIntWidom,  pc%ExpMinusBetaEnLaMin/this%Density)
            call Update( pc%SumChemPotThermoIntWidomV, pc%ExpMinusBetaEnLaMin/this%Density/this%Density)
            call Update( pc%SumChemPotV, pc%BinsIntdEndLa(pc%NBins-1)/this%Temperature-log(pc%SumChemPotThermoIntWidom%Average/(pc%Fraction+1._RK/real( this%NPart, RK ))))
            call Update(pc%SumHW_counter, pc%HW_counter)
            call Update(pc%SumHW_denom, pc%HW_denom)
            t=t+1
        end select
      end if
    end do

      do i = 1, this%NRealComponents
        pc => this%Component(i)
        if( pc%CalcChemPot ) then
          select case( pc%ChemPotMethod )
          case( ChemPotMethodGradIns )

#if MPI_VER > 0
          ! Per Process we calculate GI only for one component 
          if (mod(NProc,this%NGradInsComp)/=pc%NGradThis) cycle 
#endif
          call Update( pc%SumVW, this%NPart * ( this%SumVolume%Average &
&                    - pc%SumInvChemPot%Average / pc%SumInvChemPotRho%Average ) )
            call Update( pc%SumHM, 0._RK )

          case( ChemPotMethodWidom )
            call Update( pc%SumVW, this%NPart * ( pc%SumChemPotVV%Average / pc%SumChemPotV%Average &
&                      - this%SumVolume%Average ) )
          ! partial molar enthalpy
            call Update( pc%SumHM, pc%SumHW_counter%Average / pc%SumHW_denom%Average &
&                      - this%SumConfEnthalpy%Average*this%NPart  - this%RefTemperature)

          case( ChemPotMethodThermoInt )
            call Update( pc%SumVW, this%NPart * ( pc%SumChemPotThermoIntWidomV%Average / pc%SumChemPotThermoIntWidom%Average &
&                      - this%SumVolume%Average - pc%BinsIntVW(pc%NBins-1)/this%RefTemperature) )
          ! partial molar enthalpy
            call Update( pc%SumHM, pc%BinsIntHW(pc%NBins-1)/this%RefTemperature + pc%SumHW_counter%Average / pc%SumHW_denom%Average &
&                      - this%SumConfEnthalpy%Average*this%NPart  - this%RefTemperature)

          end select
        end if
      end do

    ! Update result files
    if( mod( Step, BlockSize ) == 0 ) then
      if(SimulationType .eq. MonteCarlo) then
#if MPI_VER > 0
        if (Equilibration) then
          if(CommonEqui) then
            offset = (accumulate_step/BlockSize+headers) * (11 * fields + 1) + headers
            call MPI_File_Seek((this%iounit_result), offset, MPI_SEEK_SET, ierr)
            call MPI_File_Seek((this%iounit_runave), offset, MPI_SEEK_SET, ierr)
            if (RootProc) then
              ! PROC
              write( IOBuffer, '(I11)' ) NProc
              call FileWriteNoAdvance_parallel( this%iounit_result )
              call FileWriteNoAdvance_parallel( this%iounit_runave )
      
              ! Number of steps
              write( IOBuffer, '(I11)' ) Step
              call FileWriteNoAdvance_parallel( this%iounit_result )
              call FileWriteNoAdvance_parallel( this%iounit_runave )
      
              if ( this%OptPressure ) then
                write( IOBuffer, '(" ",F10.5)' ) this%SumPressure%BlockAverage
                call FileWriteNoAdvance_parallel( this%iounit_result )
                write( IOBuffer, '(" ",F10.5)' ) this%SumPressure%Average
                call FileWriteNoAdvance_parallel( this%iounit_runave )
              else
                write( IOBuffer, '(" ",F10.5)' ) this%RefPressure
                call FileWriteNoAdvance_parallel( this%iounit_result )
                write( IOBuffer, '(" ",F10.5)' ) this%RefPressure
                call FileWriteNoAdvance_parallel( this%iounit_runave )
              end if
      
              ! Density
              write( IOBuffer, '(" ",F10.5)' ) this%SumDensity%BlockAverage
              call FileWriteNoAdvance_parallel( this%iounit_result )
              write( IOBuffer, '(" ",F10.5)' ) this%SumDensity%Average
              call FileWriteNoAdvance_parallel( this%iounit_runave )
      
              ! Temperature
              write( IOBuffer, '(" ",F10.5)' ) this%SumTemperature%BlockAverage
              call FileWriteNoAdvance_parallel( this%iounit_result )
              write( IOBuffer, '(" ",F10.5)' ) this%SumTemperature%Average
              call FileWriteNoAdvance_parallel( this%iounit_runave )

              ! Potential energy
              write( IOBuffer, '(" ",F12.5)' ) this%SumEPot%BlockAverage
              call FileWriteNoAdvance_parallel( this%iounit_result )
              write( IOBuffer, '(" ",F12.5)' ) this%SumEPot%Average
              call FileWriteNoAdvance_parallel( this%iounit_runave )
      
              ! Enthalpy
              write( IOBuffer, '(" ",F12.5)' ) this%SumEnthalpy%BlockAverage
              call FileWriteNoAdvance_parallel( this%iounit_result )
              write( IOBuffer, '(" ",F12.5)' ) this%SumEnthalpy%Average
              call FileWriteNoAdvance_parallel( this%iounit_runave )
      
      if (printIDF) then
        ! EPotInter
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotInter%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotInter%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra_Bond
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Bond%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Bond%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra_Angle
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Angle%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Angle%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra_Dihedral
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Dihedral%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Dihedral%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra_Nonbonded
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Nonbonded%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Nonbonded%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! VirialIntra
        write( IOBuffer, '(" ",F12.5) ' ) this%SumVirialIntra%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumVirialIntra%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! VirialInter
        write( IOBuffer, '(" ",F14.5) ' ) this%SumVirialInter%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F14.5) ' ) this%SumVirialInter%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )
      end if

              ! Chemical potential
              do i = 1, this%NRealComponents
                pc => this%Component(i)
                if( pc%ChemPotMethod .ne. ChemPotMethodNone ) then
                    write( IOBuffer, '(" ",F10.5)' ) 0._RK
                    call FileWriteNoAdvance_parallel( this%iounit_result )
                    call FileWriteNoAdvance_parallel( this%iounit_runave )
                end if
              end do
      
            ! Partial molar volume
              do i = 1, this%NRealComponents
                pc => this%Component(i)
                if( pc%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT) then
                    write( IOBuffer, '(" ",F10.4)' ) 0._RK
                    call FileWriteNoAdvance_parallel( this%iounit_result )
                    call FileWriteNoAdvance_parallel( this%iounit_runave )
      
                end if
              end do
      
             ! Partial molar enthalphy 
              do i = 1, this%NRealComponents
                pc => this%Component(i)
                if( (pc%ChemPotMethod .eq. ChemPotMethodWidom .or. pc%ChemPotMethod .eq. ChemPotMethodThermoInt) .and. EnsembleType .eq. EnsembleTypeNPT) then
                    write( IOBuffer, '(" ",F10.4)' ) 0._RK
                    call FileWriteNoAdvance_parallel( this%iounit_result )
                    call FileWriteNoAdvance_parallel( this%iounit_runave )
                end if
              end do

            ! Number of particles in ensemble
              if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
                write( IOBuffer, '(" ",F10.2)' ) this%SumNPart%BlockAverage
                call FileWriteNoAdvance_parallel( this%iounit_result )
                write( IOBuffer, '(" ",F10.2)' ) this%SumNPart%Average
                call FileWriteNoAdvance_parallel( this%iounit_runave )
        
                ! Mole fraction of each component
                do i = 1, this%NComponents
                  pc => this%Component(i)
                  write( IOBuffer, '(" ",F10.5)' ) pc%SumFraction%BlockAverage
                  call FileWriteNoAdvance_parallel( this%iounit_result )
                  write( IOBuffer, '(" ",F10.5)' ) pc%SumFraction%Average
                  call FileWriteNoAdvance_parallel( this%iounit_runave )
                end do
              end if
        
#if CONSTR == 0
               write( IOBuffer, '()' )
               call FileWriteNoAdvance_parallel( this%iounit_result )
               call FileWriteNoAdvance_parallel( this%iounit_runave )
#else
               this%consup = .true.
#endif
               write( IOBuffer, '(A)' )new_line('a')
               call FileWriteNoAdvance_parallel( this%iounit_result )
               call FileWriteNoAdvance_parallel( this%iounit_runave )
            endif
          else ! No CommonEqui
            offset = (NProc  + ((accumulate_step/BlockSize)) * NProcs +headers) * (11 * fields + 1) + headers
            call MPI_File_Seek((this%iounit_result), offset, MPI_SEEK_SET, ierr)
            call MPI_File_Seek((this%iounit_runave), offset, MPI_SEEK_SET, ierr)
            ! PROC
            write( IOBuffer, '(I11)' ) NProc
            call FileWriteNoAdvance_parallel( this%iounit_result )
            call FileWriteNoAdvance_parallel( this%iounit_runave )
    
            ! Number of steps
            write( IOBuffer, '(I11)' ) ((Step/BlockSize) - 1) * (BlockSize * NProcs)  + (NProc + 1) * BlockSize
            call FileWriteNoAdvance_parallel( this%iounit_result )
            call FileWriteNoAdvance_parallel( this%iounit_runave )
    
            if ( this%OptPressure ) then
              write( IOBuffer, '(" ",F10.5)' ) this%SumPressure%BlockAverage
              call FileWriteNoAdvance_parallel( this%iounit_result )
              write( IOBuffer, '(" ",F10.5)' ) this%SumPressure%Average
              call FileWriteNoAdvance_parallel( this%iounit_runave )
            else
              write( IOBuffer, '(" ",F10.5)' ) this%RefPressure
              call FileWriteNoAdvance_parallel( this%iounit_result )
              write( IOBuffer, '(" ",F10.5)' ) this%RefPressure
              call FileWriteNoAdvance_parallel( this%iounit_runave )
            end if
    
            ! Density
            write( IOBuffer, '(" ",F10.5)' ) this%SumDensity%BlockAverage
            call FileWriteNoAdvance_parallel( this%iounit_result )
            write( IOBuffer, '(" ",F10.5)' ) this%SumDensity%Average
            call FileWriteNoAdvance_parallel( this%iounit_runave )
    
            ! Temperature
            write( IOBuffer, '(" ",F10.5)' ) this%SumTemperature%BlockAverage
            call FileWriteNoAdvance_parallel( this%iounit_result )
            write( IOBuffer, '(" ",F10.5)' ) this%SumTemperature%Average
            call FileWriteNoAdvance_parallel( this%iounit_runave )

            ! Potential energy
            write( IOBuffer, '(" ",F12.5)' ) this%SumEPot%BlockAverage
            call FileWriteNoAdvance_parallel( this%iounit_result )
            write( IOBuffer, '(" ",F12.5)' ) this%SumEPot%Average
            call FileWriteNoAdvance_parallel( this%iounit_runave )
    
            ! Enthalpy
            write( IOBuffer, '(" ",F12.5)' ) this%SumEnthalpy%BlockAverage
            call FileWriteNoAdvance_parallel( this%iounit_result )
            write( IOBuffer, '(" ",F12.5)' ) this%SumEnthalpy%Average
            call FileWriteNoAdvance_parallel( this%iounit_runave )

      if (printIDF) then
        ! EPotInter
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotInter%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotInter%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra_Bond
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Bond%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Bond%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra_Angle
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Angle%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Angle%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra_Dihedral
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Dihedral%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Dihedral%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra_Nonbonded
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Nonbonded%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Nonbonded%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! VirialIntra
        write( IOBuffer, '(" ",F12.5) ' ) this%SumVirialIntra%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumVirialIntra%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! VirialInter
        write( IOBuffer, '(" ",F14.5) ' ) this%SumVirialInter%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F14.5) ' ) this%SumVirialInter%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )
      end if
    
            ! Chemical potential
            do i = 1, this%NRealComponents
              pc => this%Component(i)
              if( pc%ChemPotMethod .ne. ChemPotMethodNone ) then          
                  write( IOBuffer, '(" ",F10.5)' ) 0._RK
                  call FileWriteNoAdvance_parallel( this%iounit_result )
                  call FileWriteNoAdvance_parallel( this%iounit_runave )
                  end if
            end do
    
          ! Partial molar volume
            do i = 1, this%NRealComponents
              pc => this%Component(i)
              if( pc%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT) then
                  write( IOBuffer, '(" ",F10.4)' ) 0._RK
                  call FileWriteNoAdvance_parallel( this%iounit_result )
                  call FileWriteNoAdvance_parallel( this%iounit_runave )
              end if
            end do
    
            ! Partial molar enthalphy
            do i = 1, this%NRealComponents
              pc => this%Component(i)
              if( (pc%ChemPotMethod .eq. ChemPotMethodWidom .or. pc%ChemPotMethod .eq. ChemPotMethodThermoInt) .and. EnsembleType .eq. EnsembleTypeNPT) then
                  write( IOBuffer, '(" ",F10.4)' ) 0._RK
                  call FileWriteNoAdvance_parallel( this%iounit_result )
                  call FileWriteNoAdvance_parallel( this%iounit_runave )
              end if
            end do

          ! Number of particles in ensemble
            if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
              write( IOBuffer, '(" ",F10.2)' ) this%SumNPart%BlockAverage
              call FileWriteNoAdvance_parallel( this%iounit_result )
              write( IOBuffer, '(" ",F10.2)' ) this%SumNPart%Average
              call FileWriteNoAdvance_parallel( this%iounit_runave )
      
              ! Mole fraction of each component
              do i = 1, this%NComponents
                pc => this%Component(i)
                write( IOBuffer, '(" ",F10.5)' ) pc%SumFraction%BlockAverage
                call FileWriteNoAdvance_parallel( this%iounit_result )
                write( IOBuffer, '(" ",F10.5)' ) pc%SumFraction%Average
                call FileWriteNoAdvance_parallel( this%iounit_runave )
              end do
            end if
      
#if CONSTR == 0
             write( IOBuffer, '()' )
             call FileWriteNoAdvance_parallel( this%iounit_result )
             call FileWriteNoAdvance_parallel( this%iounit_runave )
#else
             this%consup = .true.
#endif
             write( IOBuffer, '(A)' )new_line('a')
             call FileWriteNoAdvance_parallel( this%iounit_result )
             call FileWriteNoAdvance_parallel( this%iounit_runave )
          endif
          accumulate_step = accumulate_step + BlockSize
        else !production starts
          if (CommonEqui) then
            offset = (NProc  + ((Step/BlockSize)-1) * NProcs + (accumulate_step/BlockSize+headers)) * (11 * fields + 1) + headers 
          else
            offset = (NProc  + ((Step/BlockSize)-1) * NProcs + NProcs*(accumulate_step/BlockSize)+headers) * (11 * fields + 1) + headers 
          endif
          call MPI_File_Seek((this%iounit_result), offset, MPI_SEEK_SET, ierr)
          call MPI_File_Seek((this%iounit_runave), offset, MPI_SEEK_SET, ierr)
          ! PROC
          write( IOBuffer, '(I11)' ) NProc
          call FileWriteNoAdvance_parallel( this%iounit_result )
          call FileWriteNoAdvance_parallel( this%iounit_runave )
  
          ! Number of steps
          write( IOBuffer, '(I11)' ) ((Step/BlockSize) - 1) * (BlockSize * NProcs)  + (NProc + 1) * BlockSize
          call FileWriteNoAdvance_parallel( this%iounit_result )
          call FileWriteNoAdvance_parallel( this%iounit_runave )
  
          if ( this%OptPressure ) then
            write( IOBuffer, '(" ",F10.5)' ) this%SumPressure%BlockAverage
            call FileWriteNoAdvance_parallel( this%iounit_result )
            write( IOBuffer, '(" ",F10.5)' ) this%SumPressure%Average
            call FileWriteNoAdvance_parallel( this%iounit_runave )
          else
            write( IOBuffer, '(" ",F10.5)' ) this%RefPressure
            call FileWriteNoAdvance_parallel( this%iounit_result )
            write( IOBuffer, '(" ",F10.5)' ) this%RefPressure
            call FileWriteNoAdvance_parallel( this%iounit_runave )
          end if
  
          ! Density
          write( IOBuffer, '(" ",F10.5)' ) this%SumDensity%BlockAverage
          call FileWriteNoAdvance_parallel( this%iounit_result )
          write( IOBuffer, '(" ",F10.5)' ) this%SumDensity%Average
          call FileWriteNoAdvance_parallel( this%iounit_runave )
  
          ! Temperature
          write( IOBuffer, '(" ",F10.5)' ) this%SumTemperature%BlockAverage
          call FileWriteNoAdvance_parallel( this%iounit_result )
          write( IOBuffer, '(" ",F10.5)' ) this%SumTemperature%Average
          call FileWriteNoAdvance_parallel( this%iounit_runave )

          ! Potential energy
          write( IOBuffer, '(" ",F12.5)' ) this%SumEPot%BlockAverage
          call FileWriteNoAdvance_parallel( this%iounit_result )
          write( IOBuffer, '(" ",F12.5)' ) this%SumEPot%Average
          call FileWriteNoAdvance_parallel( this%iounit_runave )
  
          ! Enthalpy
          write( IOBuffer, '(" ",F12.5)' ) this%SumEnthalpy%BlockAverage
          call FileWriteNoAdvance_parallel( this%iounit_result )
          write( IOBuffer, '(" ",F12.5)' ) this%SumEnthalpy%Average
          call FileWriteNoAdvance_parallel( this%iounit_runave )

      if (printIDF) then
        ! EPotInter
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotInter%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotInter%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra_Bond
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Bond%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Bond%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra_Angle
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Angle%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Angle%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra_Dihedral
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Dihedral%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Dihedral%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! EPotIntra_Nonbonded
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Nonbonded%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Nonbonded%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! VirialIntra
        write( IOBuffer, '(" ",F12.5) ' ) this%SumVirialIntra%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumVirialIntra%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )

        ! VirialInter
        write( IOBuffer, '(" ",F14.5) ' ) this%SumVirialInter%BlockAverage
        call FileWriteNoAdvance_parallel( this%iounit_result )
        write( IOBuffer, '(" ",F14.5) ' ) this%SumVirialInter%Average
        call FileWriteNoAdvance_parallel( this%iounit_runave )
      end if
  
          ! Chemical potential 
          do i = 1, this%NRealComponents
            pc => this%Component(i)
            if( pc%ChemPotMethod .ne. ChemPotMethodNone ) then
                if( pc%NPart > 1 ) then
                  select case( pc%ChemPotMethod )
  
                  case( ChemPotMethodGradIns )
                      write( IOBuffer, '(" ",F10.5)' ) log( pc%Fraction * pc%SumInvChemPotRho%BlockAverage )
                      call FileWriteNoAdvance_parallel( this%iounit_result )
                      write( IOBuffer, '(" ",F10.5)' ) log( pc%Fraction * pc%SumInvChemPotRho%Average )
                      call FileWriteNoAdvance_parallel( this%iounit_runave )
                  case( ChemPotMethodWidom )
                    write( IOBuffer, '(" ",F10.5)' ) log( pc%Fraction / pc%SumChemPotV%BlockAverage )
                    call FileWriteNoAdvance_parallel( this%iounit_result )
                    write( IOBuffer, '(" ",F10.5)' ) log( pc%Fraction / pc%SumChemPotV%Average )
                    call FileWriteNoAdvance_parallel( this%iounit_runave )

                  case( ChemPotMethodThermoInt )
                    write( IOBuffer, '(" ",F10.5)' ) pc%SumChemPotV%BlockAverage
                    call FileWriteNoAdvance_parallel( this%iounit_result )
                    write( IOBuffer, '(" ",F10.5)' ) pc%SumChemPotV%Average
                    call FileWriteNoAdvance_parallel( this%iounit_runave )
                  end select
  
                else
                  select case( pc%ChemPotMethod )
  
                  case( ChemPotMethodGradIns )
                    write( IOBuffer, '(" ",F10.5)' ) log( pc%SumInvChemPotRho%BlockAverage )
                    call FileWriteNoAdvance_parallel( this%iounit_result )
                    write( IOBuffer, '(" ",F10.5)' ) log( pc%SumInvChemPotRho%Average )
                    call FileWriteNoAdvance_parallel( this%iounit_runave )
  
                  case( ChemPotMethodWidom )
                    write( IOBuffer, '(" ",F10.5)' ) -log( pc%SumChemPotV%BlockAverage )
                    call FileWriteNoAdvance_parallel( this%iounit_result )
                    write( IOBuffer, '(" ",F10.5)' ) -log( pc%SumChemPotV%Average )
                    call FileWriteNoAdvance_parallel( this%iounit_runave )

                  case( ChemPotMethodThermoInt )
                    write( IOBuffer, '(" ",F10.5)' ) pc%SumChemPotV%BlockAverage
                    call FileWriteNoAdvance_parallel( this%iounit_result )
                    write( IOBuffer, '(" ",F10.5)' ) pc%SumChemPotV%Average
                    call FileWriteNoAdvance_parallel( this%iounit_runave )
                  end select
                end if
            end if
          end do
  
        ! Partial molar volume
          do i = 1, this%NRealComponents
            pc => this%Component(i)
            if( pc%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT) then
                write( IOBuffer, '(" ",F10.4)' ) pc%SumVW%BlockAverage
                call FileWriteNoAdvance_parallel( this%iounit_result )
                write( IOBuffer, '(" ",F10.4)' ) pc%SumVW%Average
                call FileWriteNoAdvance_parallel( this%iounit_runave )
            end if
          end do
  
          ! Partial molar enthalphy
          do i = 1, this%NRealComponents
            pc => this%Component(i)
            if( pc%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT) then
                write( IOBuffer, '(" ",F10.4)' ) pc%SumHM%BlockAverage
                call FileWriteNoAdvance_parallel( this%iounit_result )     
                write( IOBuffer, '(" ",F10.4)' ) pc%SumHM%Average
                call FileWriteNoAdvance_parallel( this%iounit_runave )
            end if
          end do

        ! Number of particles in ensemble
          if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
            write( IOBuffer, '(" ",F10.2)' ) this%SumNPart%BlockAverage
            call FileWriteNoAdvance_parallel( this%iounit_result )
            write( IOBuffer, '(" ",F10.2)' ) this%SumNPart%Average
            call FileWriteNoAdvance_parallel( this%iounit_runave )
    
            ! Mole fraction of each component
            do i = 1, this%NComponents
              pc => this%Component(i)
              write( IOBuffer, '(" ",F10.5)' ) pc%SumFraction%BlockAverage
              call FileWriteNoAdvance_parallel( this%iounit_result )
              write( IOBuffer, '(" ",F10.5)' ) pc%SumFraction%Average
              call FileWriteNoAdvance_parallel( this%iounit_runave )
            end do
          end if
    
#if CONSTR == 0
           write( IOBuffer, '()' )
           call FileWriteNoAdvance_parallel( this%iounit_result )
           call FileWriteNoAdvance_parallel( this%iounit_runave )
#else
          this%consup = .true.
#endif
          write( IOBuffer, '(A)' )new_line('a')
          call FileWriteNoAdvance_parallel( this%iounit_result )
          call FileWriteNoAdvance_parallel( this%iounit_runave )
        endif
#else 
!MPI=0
        ! Number of steps
        write( IOBuffer, '(I9)' ) Step
        call FileWriteNoAdvance( this%iounit_result )
        call FileWriteNoAdvance( this%iounit_runave )

      ! Pressure
        if ( this%OptPressure ) then
          write( IOBuffer, '(" ",F10.5)' ) this%SumPressure%BlockAverage
          call FileWriteNoAdvance( this%iounit_result )
          write( IOBuffer, '(" ",F10.5)' ) this%SumPressure%Average
          call FileWriteNoAdvance( this%iounit_runave )
        else
          write( IOBuffer, '(" ",F10.5)' ) this%RefPressure
          call FileWriteNoAdvance( this%iounit_result )
          write( IOBuffer, '(" ",F10.5)' ) this%RefPressure
          call FileWriteNoAdvance( this%iounit_runave )
        end if

        ! Density
        write( IOBuffer, '(" ",F10.5)' ) this%SumDensity%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F10.5)' ) this%SumDensity%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! Temperature
        write( IOBuffer, '(" ",F10.5)' ) this%SumTemperature%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F10.5)' ) this%SumTemperature%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! Potential energy
        write( IOBuffer, '(" ",F12.5)' ) this%SumEPot%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5)' ) this%SumEPot%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! Enthalpy
        write( IOBuffer, '(" ",F12.5)' ) this%SumEnthalpy%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5)' ) this%SumEnthalpy%Average
        call FileWriteNoAdvance( this%iounit_runave )

      if (printIDF) then
        ! EPotInter
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotInter%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotInter%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! EPotIntra
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! EPotIntra_Bond
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Bond%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Bond%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! EPotIntra_Angle
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Angle%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Angle%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! EPotIntra_Dihedral
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Dihedral%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Dihedral%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! EPotIntra_Nonbonded
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Nonbonded%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Nonbonded%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! VirialIntra
        write( IOBuffer, '(" ",F12.5) ' ) this%SumVirialIntra%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumVirialIntra%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! VirialInter
        write( IOBuffer, '(" ",F14.5) ' ) this%SumVirialInter%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F14.5) ' ) this%SumVirialInter%Average
        call FileWriteNoAdvance( this%iounit_runave )
      end if

        ! Chemical potential
        do i = 1, this%NRealComponents
          pc => this%Component(i)
          if( pc%ChemPotMethod .ne. ChemPotMethodNone ) then
            if( Equilibration ) then
              write( IOBuffer, '(" ",F10.5)' ) 0._RK
              call FileWriteNoAdvance( this%iounit_result )
              call FileWriteNoAdvance( this%iounit_runave )
            else
              if( pc%NPart > 1 ) then
                select case( pc%ChemPotMethod )

                case( ChemPotMethodGradIns )
                  write( IOBuffer, '(" ",F10.5)' ) log( pc%Fraction * pc%SumInvChemPotRho%BlockAverage )
                  call FileWriteNoAdvance( this%iounit_result )
                  write( IOBuffer, '(" ",F10.5)' ) log( pc%Fraction * pc%SumInvChemPotRho%Average )
                  call FileWriteNoAdvance( this%iounit_runave )

                case( ChemPotMethodWidom )
                  write( IOBuffer, '(" ",F10.5)' ) log( pc%Fraction / pc%SumChemPotV%BlockAverage )
                  call FileWriteNoAdvance( this%iounit_result )
                  write( IOBuffer, '(" ",F10.5)' ) log( pc%Fraction / pc%SumChemPotV%Average )
                  call FileWriteNoAdvance( this%iounit_runave )

                case( ChemPotMethodThermoInt )
                  write( IOBuffer, '(" ",F10.5)' ) pc%SumChemPotV%BlockAverage
                  call FileWriteNoAdvance( this%iounit_result )
                  write( IOBuffer, '(" ",F10.5)' ) pc%SumChemPotV%Average
                  call FileWriteNoAdvance( this%iounit_runave )

                end select

              else
                select case( pc%ChemPotMethod )

                case( ChemPotMethodGradIns )
                  write( IOBuffer, '(" ",F10.5)' ) log( pc%SumInvChemPotRho%BlockAverage )
                  call FileWriteNoAdvance( this%iounit_result )
                  write( IOBuffer, '(" ",F10.5)' ) log( pc%SumInvChemPotRho%Average )
                  call FileWriteNoAdvance( this%iounit_runave )

                case( ChemPotMethodWidom )
                  write( IOBuffer, '(" ",F10.5)' ) -log( pc%SumChemPotV%BlockAverage )
                  call FileWriteNoAdvance( this%iounit_result )
                  write( IOBuffer, '(" ",F10.5)' ) -log( pc%SumChemPotV%Average )
                  call FileWriteNoAdvance( this%iounit_runave )

                case( ChemPotMethodThermoInt )
                  write( IOBuffer, '(" ",F10.5)' ) pc%SumChemPotV%BlockAverage
                  call FileWriteNoAdvance( this%iounit_result )
                  write( IOBuffer, '(" ",F10.5)' ) pc%SumChemPotV%Average
                  call FileWriteNoAdvance( this%iounit_runave )

                end select
              end if
            end if
          end if
        end do

      ! Partial molar volume
        do i = 1, this%NRealComponents
          pc => this%Component(i)
          if( pc%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT) then
            if( Equilibration ) then
              write( IOBuffer, '(" ",F10.4)' ) 0._RK
              call FileWriteNoAdvance( this%iounit_result )
              call FileWriteNoAdvance( this%iounit_runave )
            else
              write( IOBuffer, '(" ",F10.4)' ) pc%SumVW%BlockAverage
              call FileWriteNoAdvance( this%iounit_result )
              write( IOBuffer, '(" ",F10.4)' ) pc%SumVW%Average
              call FileWriteNoAdvance( this%iounit_runave )
           end if
          end if
        end do

        ! Partial molar enthalphy
        do i = 1, this%NRealComponents
          pc => this%Component(i)
          if( pc%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT) then
            if( Equilibration ) then
              write( IOBuffer, '(" ",F10.4)' ) 0._RK
              call FileWriteNoAdvance( this%iounit_result )
              call FileWriteNoAdvance( this%iounit_runave )
            else
              write( IOBuffer, '(" ",F10.4)' ) pc%SumHM%BlockAverage
              call FileWriteNoAdvance( this%iounit_result )     
              write( IOBuffer, '(" ",F10.4)' ) pc%SumHM%Average
              call FileWriteNoAdvance( this%iounit_runave )
            end if
          end if
        end do

      ! Number of particles in ensemble
        if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
          write( IOBuffer, '(" ",F10.2)' ) this%SumNPart%BlockAverage
          call FileWriteNoAdvance( this%iounit_result )
          write( IOBuffer, '(" ",F10.2)' ) this%SumNPart%Average
          call FileWriteNoAdvance( this%iounit_runave )
  
          ! Mole fraction of each component
          do i = 1, this%NComponents
            pc => this%Component(i)
            write( IOBuffer, '(" ",F10.5)' ) pc%SumFraction%BlockAverage
            call FileWriteNoAdvance( this%iounit_result )
            write( IOBuffer, '(" ",F10.5)' ) pc%SumFraction%Average
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
#endif
      else !MD
        ! Number of steps
        write( IOBuffer, '(I9)' ) Step
        call FileWriteNoAdvance( this%iounit_result )
        call FileWriteNoAdvance( this%iounit_runave )

        ! Displacement
        value = 0._RK
        do i = 1, this%NComponents
          value = value + sum( this%Component(i)%Disp(:, :)**2 )
        end do
        value = value * this%BoxLength**2 / ( 6._RK * this%NPart * TimeStep * Step )
        write( IOBuffer, '(" ",F8.3)' ) value
        call FileWriteNoAdvance( this%iounit_runave )

        ! Pressure
        write( IOBuffer, '(" ",F10.5)' ) this%SumPressure%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F10.5)' ) this%SumPressure%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! Density
        write( IOBuffer, '(" ",F10.5)' ) this%SumDensity%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F10.5)' ) this%SumDensity%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! Temperature
        write( IOBuffer, '(" ",F10.5)' ) this%SumTemperature%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F10.5)' ) this%SumTemperature%Average
        call FileWriteNoAdvance( this%iounit_runave )

#if OSMOP > 0
        ! OsmoticPressure
        write( IOBuffer, '(F10.5)' ) this%SumOsmoticPressure%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(F10.5)' ) this%SumOsmoticPressure%Average
        call FileWriteNoAdvance( this%iounit_runave )
#endif

        ! Potential energy
        write( IOBuffer, '(" ",F12.5)' ) this%SumEPot%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5)' ) this%SumEPot%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! Enthalpy
        write( IOBuffer, '(" ",F12.5)' ) this%SumEnthalpy%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5)' ) this%SumEnthalpy%Average
        call FileWriteNoAdvance( this%iounit_runave )

      if (printIDF) then
        ! EPotInter
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotInter%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotInter%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! EPotIntra
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! EPotIntra_Bond
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Bond%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Bond%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! EPotIntra_Angle
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Angle%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Angle%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! EPotIntra_Dihedral
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Dihedral%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Dihedral%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! EPotIntra_Nonbonded
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Nonbonded%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumEPotIntra_Nonbonded%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! VirialIntra
        write( IOBuffer, '(" ",F12.5) ' ) this%SumVirialIntra%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F12.5) ' ) this%SumVirialIntra%Average
        call FileWriteNoAdvance( this%iounit_runave )

        ! VirialInter
        write( IOBuffer, '(" ",F14.5) ' ) this%SumVirialInter%BlockAverage
        call FileWriteNoAdvance( this%iounit_result )
        write( IOBuffer, '(" ",F14.5) ' ) this%SumVirialInter%Average
        call FileWriteNoAdvance( this%iounit_runave )
      end if

        ! Chemical potential
        do i = 1, this%NRealComponents
          pc => this%Component(i)
          if( pc%ChemPotMethod .ne. ChemPotMethodNone ) then
            if( Equilibration ) then
              write( IOBuffer, '(" ",F10.5)' ) 0._RK
              call FileWriteNoAdvance( this%iounit_result )
              call FileWriteNoAdvance( this%iounit_runave )
            else
              if( pc%NPart > 1 ) then
                select case( pc%ChemPotMethod )

                case( ChemPotMethodGradIns )
                  write( IOBuffer, '(" ",F10.5)' ) log( pc%Fraction * pc%SumInvChemPotRho%BlockAverage )
                  call FileWriteNoAdvance( this%iounit_result )
                  write( IOBuffer, '(" ",F10.5)' ) log( pc%Fraction * pc%SumInvChemPotRho%Average )
                  call FileWriteNoAdvance( this%iounit_runave )

                case( ChemPotMethodWidom )
                  write( IOBuffer, '(" ",F10.5)' ) log( pc%Fraction / pc%SumChemPotV%BlockAverage )
                  call FileWriteNoAdvance( this%iounit_result )
                  write( IOBuffer, '(" ",F10.5)' ) log( pc%Fraction / pc%SumChemPotV%Average )
                  call FileWriteNoAdvance( this%iounit_runave )

                case( ChemPotMethodThermoInt )
                  write( IOBuffer, '(" ",F10.5)' ) pc%SumChemPotV%BlockAverage
                  call FileWriteNoAdvance( this%iounit_result )
                  write( IOBuffer, '(" ",F10.5)' ) pc%SumChemPotV%Average
                  call FileWriteNoAdvance( this%iounit_runave )
                end select

              else
                select case( pc%ChemPotMethod )

                case( ChemPotMethodGradIns )
                  write( IOBuffer, '(" ",F10.5)' ) log( pc%SumInvChemPotRho%BlockAverage )
                  call FileWriteNoAdvance( this%iounit_result )
                  write( IOBuffer, '(" ",F10.5)' ) log( pc%SumInvChemPotRho%Average )
                  call FileWriteNoAdvance( this%iounit_runave )

                case( ChemPotMethodWidom )
                  write( IOBuffer, '(" ",F10.5)' ) -log( pc%SumChemPotV%BlockAverage )
                  call FileWriteNoAdvance( this%iounit_result )
                  write( IOBuffer, '(" ",F10.5)' ) -log( pc%SumChemPotV%Average )
                  call FileWriteNoAdvance( this%iounit_runave )

                case( ChemPotMethodThermoInt )
                  write( IOBuffer, '(" ",F10.5)' ) pc%SumChemPotV%BlockAverage
                  call FileWriteNoAdvance( this%iounit_result )
                  write( IOBuffer, '(" ",F10.5)' ) pc%SumChemPotV%Average
                  call FileWriteNoAdvance( this%iounit_runave )
                end select
              end if
            end if
          end if
        end do

      ! Partial molar volume
        do i = 1, this%NRealComponents
          pc => this%Component(i)
          if( pc%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT) then
            if( Equilibration ) then
              write( IOBuffer, '(" ",F10.4)' ) 0._RK
              call FileWriteNoAdvance( this%iounit_result )
              call FileWriteNoAdvance( this%iounit_runave )

            else
              write( IOBuffer, '(" ",F10.4)' ) pc%SumVW%BlockAverage
              call FileWriteNoAdvance( this%iounit_result )
              write( IOBuffer, '(" ",F10.4)' ) pc%SumVW%Average
              call FileWriteNoAdvance( this%iounit_runave )
            end if
          end if
        end do

        ! Partial molar enthalphy
        do i = 1, this%NRealComponents
          pc => this%Component(i)
          if( pc%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT) then
            if( Equilibration ) then
              write( IOBuffer, '(" ",F10.4)' ) 0._RK
              call FileWriteNoAdvance( this%iounit_result )
              call FileWriteNoAdvance( this%iounit_runave )
            else
              write( IOBuffer, '(" ",F10.4)' ) pc%SumHM%BlockAverage
              call FileWriteNoAdvance( this%iounit_result )     
              write( IOBuffer, '(" ",F10.4)' ) pc%SumHM%Average
              call FileWriteNoAdvance( this%iounit_runave )
            end if
          end if
        end do

#if HBOND > 0
        do i = 1, this%NComponents
          write( IOBuffer, '(" ", F10.4)' ) this%SumHBond0(i)%BlockAverage
          call FileWriteNoAdvance( this%iounit_result )
          write( IOBuffer, '(" ", F10.4)' ) this%SumHBond0(i)%Average
          call FileWriteNoAdvance( this%iounit_runave )
        end do
        do i = 1, this%NComponents
          do  j = 1, this%NComponents
            write( IOBuffer, '("   ", F10.4)' ) this%SumHBond1(i,j)%BlockAverage
            call FileWriteNoAdvance( this%iounit_result )
            write( IOBuffer, '("   ", F10.4)' ) this%SumHBond1(i,j)%Average
            call FileWriteNoAdvance( this%iounit_runave )
          end do
        end do
        do i = 1, this%NComponents
          do  j = 1, this%NComponents
            do k = j, this%NComponents
              write( IOBuffer, '("      ", F10.4)' ) this%SumHBond2(i,j,k)%BlockAverage
              call FileWriteNoAdvance( this%iounit_result )
              write( IOBuffer, '("      ", F10.4)' ) this%SumHBond2(i,j,k)%Average
              call FileWriteNoAdvance( this%iounit_runave )
            end do
          end do
        end do
        do i = 1, this%NComponents
          do  j = 1, this%NComponents
            do k = j, this%NComponents
              do  l = k, this%NComponents
                write( IOBuffer, '("         ", F10.4)' ) this%SumHBond3(i,j,k,l)%BlockAverage
                call FileWriteNoAdvance( this%iounit_result )
                write( IOBuffer, '("         ", F10.4)' ) this%SumHBond3(i,j,k,l)%Average
                call FileWriteNoAdvance( this%iounit_runave )
              end do
            end do
          end do
        end do
        do i = 1, this%NComponents
          write( IOBuffer, '(" ", F10.4)' ) this%SumHBondN(i)%BlockAverage
          call FileWriteNoAdvance( this%iounit_result )
          write( IOBuffer, '(" ", F10.4)' ) this%SumHBondN(i)%Average
          call FileWriteNoAdvance( this%iounit_runave )
        end do
#endif

#if OSMOP > 0
        !Density Profile
        do i = 1, this%NComponents
          pc => this%Component(i)
          do j = 1, NBinsDen
            write( IOBuffer, '(F10.4)' ) pc%SumDenProfile(j)%BlockAverage
            call FileWriteNoAdvance( this%iounit_result )
            write( IOBuffer, '(F10.4)' ) pc%SumDenProfile(j)%Average
            call FileWriteNoAdvance( this%iounit_runave )
          end do
        end do

#if OSMOP == 2
        !Pressure Profile
        do j = 1, NBinsDen
            write( IOBuffer, '(F10.4)' ) this%SumPressureProfile(j)%BlockAverage
            call FileWriteNoAdvance( this%iounit_result )
            write( IOBuffer, '(F10.4)' ) this%SumPressureProfile(j)%Average 
            call FileWriteNoAdvance( this%iounit_runave )
        end do

        !Chemical Potential Profile
        do i = 1, this%NRealComponents
          pc => this%Component(i)
          if( pc%ChemPotMethod .eq. ChemPotMethodWidom ) then
            if( Equilibration ) then
              do j = 1, NBinsDen
                write( IOBuffer, '(F10.5)' ) 0._RK
                call FileWriteNoAdvance( this%iounit_result )
                call FileWriteNoAdvance( this%iounit_runave )
              end do
            else
              do j = 1, NBinsDen
                write( IOBuffer, '(F10.5)' ) &
&                      log( pc%SumDenProfile(j)%Average / pc%SumChemPotProfile(j)%BlockAverage )
                call FileWriteNoAdvance( this%iounit_result )
                write( IOBuffer, '(F10.5)' ) &
&                      log( pc%SumDenProfile(j)%Average / pc%SumChemPotProfile(j)%Average )
                call FileWriteNoAdvance( this%iounit_runave )
              end do
            end if
          end if
        end do
#endif
#endif

      ! Number of particles in ensemble
        if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
          write( IOBuffer, '(" ",F10.2)' ) this%SumNPart%BlockAverage
          call FileWriteNoAdvance( this%iounit_result )
          write( IOBuffer, '(" ",F10.2)' ) this%SumNPart%Average
          call FileWriteNoAdvance( this%iounit_runave )
  
          ! Mole fraction of each component
          do i = 1, this%NComponents
            pc => this%Component(i)
            write( IOBuffer, '(" ",F10.5)' ) pc%SumFraction%BlockAverage
            call FileWriteNoAdvance( this%iounit_result )
            write( IOBuffer, '(" ",F10.5)' ) pc%SumFraction%Average
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
    end if

#if  TRANS == 1
    ! Transport properties !TRANSPORT_start
    if( ( this%Mmess > 0 ) .and. ( mod(this%Mmess, this%Nviewcf) == 0 )&
&       .and. (mod(Step + this%NStepCorr -1, this%NSpanCF*this%NStepCorr) == 0) ) then

      rewind( this%iounit_rescf )
      write( IOBuffer, '("  TIME[ps]")' )
      call FileWriteNoAdvance( this%iounit_rescf )

      if(this%Ncomponents>1)then
        do i=1,this%NComponents*this%NComponents
          if( i < 10 ) then
            write( IOBuffer, '(T9, " L_ij_", I1)') i
          else
            write( IOBuffer, '(T9, "L_ij_", I2)') i
          end if
          call FileWriteNoAdvance( this%iounit_rescf )
        end do
      end if
      
      do i = 1, this%NComponents
        if( i < 10 ) then
          write( IOBuffer, '(T10," D_i_",I1)' ) i
        else
          write( IOBuffer, '(T10,"D_i_",I2)' ) i
        end if
        call FileWriteNoAdvance( this%iounit_rescf )
      end do

      write( IOBuffer, '(T13,"VS")' )
      call FileWriteNoAdvance( this%iounit_rescf )

      if (this%Bulkviscosity) then
        write( IOBuffer, '(T13,"VB")' )
        call FileWriteNoAdvance( this%iounit_rescf )
      end if

      if (this%Conductivity) then
        write( IOBuffer, '(T13,"CO")' )
        call FileWriteNoAdvance( this%iounit_rescf )
        if(this%NComponents==2)then
            write( IOBuffer, '(T10,"Th_Diff")' )
            call FileWriteNoAdvance( this%iounit_rescf )
        end if
      end if

      if (this%EConductivity) then
        write( IOBuffer, '(T13,"EC")' )
        call FileWriteNoAdvance( this%iounit_rescf )
      end if

      if( this%Ncomponents > 1 ) then
        do i=1,this%NComponents*this%NComponents
           if( i < 10 ) then
             write( IOBuffer, '(T7," Int_Lij_",I1)')i
           else
             write( IOBuffer, '(T7,"Int_Lij_",I1)')i
           end if
           call FileWriteNoAdvance( this%iounit_rescf )
       end do
      end if

      do i = 1, this%NComponents
         if( i < 10 ) then
           write( IOBuffer, '(T7," IntD_i_",I1)' ) i
         else
           write( IOBuffer, '(T7,"IntD_i_",I2)' ) i
         end if
         call FileWriteNoAdvance( this%iounit_rescf )
      end do

      write( IOBuffer, '(T9,"Int VS")' )
      call FileWriteNoAdvance( this%iounit_rescf )

      if (this%Bulkviscosity) then
        write( IOBuffer, '(T9,"Int VB")' )
        call FileWriteNoAdvance( this%iounit_rescf )
      end if

      if (this%Conductivity) then
        write( IOBuffer, '(T10,"Int C ")' )
        call FileWriteNoAdvance( this%iounit_rescf )
        if (this%NComponents == 2 ) then
          write( IOBuffer, '(T10,"Int Th_Diff")' )
          call FileWriteNoAdvance( this%iounit_rescf )
        end if       
      end if

      if (this%EConductivity) then
        write( IOBuffer, '(T9,"Int EC")' )
        call FileWriteNoAdvance( this%iounit_rescf )
      end if

      call FileWriteBlank( this%iounit_rescf )

      ! integration time
      do i  = 1, this%NCorr
        value = this%TimeStepCorr*UnitTime/1E-12_RK
        write( IOBuffer, '(" ",F10.5)' ) (i-1)*value
        call FileWriteNoAdvance( this%iounit_rescf )

!         ! Onsager Diffusion coefficients
        if(this%Ncomponents>1)then
          do j=1,this%NComponents*this%NComponents
              write( IOBuffer, '(T5, F10.5)' )  this%average_lamda(j, i)/this%average_lamda(j,1)
              call FileWriteNoAdvance( this%iounit_rescf )
          end do
        end if

        ! Self-diffusion coefficients
        do j = 1, this%NComponents
          write( IOBuffer, '(T5, F10.5)' ) this%average_cf_d(j,i)/this%average_cf_d(j,1)
          call FileWriteNoAdvance( this%iounit_rescf )
        end do

        ! Shear viscosity
        write( IOBuffer, '(T5, F10.5)' ) this%average_cf_vs(i)/this%average_cf_vs(1)
        call FileWriteNoAdvance( this%iounit_rescf )

        ! Bulk viscosity
        if (this%Bulkviscosity) then
          write( IOBuffer, '(T5, F10.5)' ) this%average_cf_vb(i)/this%average_cf_vb(1)
          call FileWriteNoAdvance( this%iounit_rescf )
        end if

        ! Thermal conductivity and thermal diffusion
        if (this%Conductivity) then
          write( IOBuffer, '(T5, F10.5)' )  this%average_cf_c(i)/this%average_cf_c(1)
          call FileWriteNoAdvance( this%iounit_rescf )
          if (this%NComponents==2)then
            write( IOBuffer, '(T5,F10.5)' ) this%average_cf_soret(i)/this%average_cf_soret(1)
            call FileWriteNoAdvance( this%iounit_rescf )
          end if
        end if

        ! Electric Conductivity
        if (this%EConductivity) then
          write( IOBuffer, '(T5, F10.5)' )  this%average_cf_ec(i)/this%average_cf_ec(1)
          call FileWriteNoAdvance( this%iounit_rescf )
        end if

        ! integral ======================================================!
        value = dsqrt(UnitEnergy/UnitMass)*UnitLength/1E-10_RK

        ! Onsager Diffusion coefficients
        if( this%Ncomponents > 1) then
          do j = 1, this%NComponents*this%NComponents
             write( IOBuffer, '(T5, F10.4)' )  this%average_sinte_lamda(j,i)*value !this%sinte_lamda(j,i) / this%sinte_lamda(j,this%Ncorr)* value
             call FileWriteNoAdvance( this%iounit_rescf )
          end do
        end if

        ! Self-diffusion coefficient
        do j = 1, this%NComponents
          write( IOBuffer, '(T5, F10.4)' ) this%average_sinte_i(j,i)* value !this%sinte_i(j,i) / this%sinte_i(j,this%NCorr) * this%selfd_i(j) * value
          call FileWriteNoAdvance( this%iounit_rescf )
        end do

       !viscosity
        value = dsqrt(UnitEnergy*UnitMass)/UnitLength**2/1E-4_RK

       !shear
        write( IOBuffer, '(T5, F10.5)' )  this%average_sinte_vs(i)* value !this%sinte_vs(i) / this%sinte_vs(this%NCorr) * this%visco_s * value
        call FileWriteNoAdvance( this%iounit_rescf )

       ! bulk
        write( IOBuffer, '(T5, F10.5)' ) this%average_sinte_vb(i)*value !this%sinte_vb(i) / this%sinte_vb(this%NCorr) * this%visco_b * value
        call FileWriteNoAdvance( this%iounit_rescf )

       ! thermal conductivity
        if (this%Conductivity) then
          value = dsqrt(UnitEnergy/UnitMass)*kBoltzmann/UnitLength**2
          write( IOBuffer, '(T5, F10.5)' ) this%average_sinte_c(i)*value !this%sinte_c(i) / this%sinte_c(this%NCorr) * this%conduct * value
          call FileWriteNoAdvance( this%iounit_rescf )
          if ( this%NComponents == 2) then
            value = dsqrt(UnitEnergy/UnitMass)*UnitLength*(kBoltzmann/UnitEnergy)/1E-12_RK
            write( IOBuffer, '(T5, F10.4)' ) this%average_sinte_soret(i)*value !/ this%average_sinte_soret(this%NCorr) * this%soret * value
            call FileWriteNoAdvance( this%iounit_rescf )
          end if
        end if

        ! electric conductivity
        if (this%EConductivity) then
           value = ElementaryCharge**2 /(dsqrt(UnitEnergy*UnitMass) * UnitLength**2)
           write( IOBuffer, '(T5, F10.5)' ) this%sinte_ec(i) / this%sinte_ec(this%NCorr) * this%econduct * value
           call FileWriteNoAdvance( this%iounit_rescf )
        end if

        call FileWriteBlank( this%iounit_rescf )
      end do

#if ARCH == 2
      call flush( this%iounit_rescf )
#endif
    end if
!TRANSPORT_END
#endif


! Exit, if specific file is in the folder!
    open( 99, file = 'stop.txt', action = 'READ', status = 'OLD', iostat = err )

#if MPI_VER > 0
! MPI Abortion
#ifdef __INTEL_COMPILER
    if( err .eq. 0 ) then
      call MPI_Bcast(err,1,MPI_INTEGER,NRootProc,Communicator,ierror)
      err = SetTerminateProgram( 1 ) 
    end if

#else
    if( err .eq. 0 ) then
      call MPI_Bcast(err,1,MPI_INTEGER,NRootProc,Communicator,ierror)
      call SetTerminateProgram
    end if
#endif

#else
! Single Abortion
#ifdef __INTEL_COMPILER
    if( err .eq. 0 ) err = SetTerminateProgram( 1 ) 
#else
    if( err .eq. 0 ) call SetTerminateProgram
#endif

#endif

#if MPI_VER > 0
! Abortion of simulation run in Karlsruhe
    call time_left(time_limit)
#endif

  end subroutine TEnsemble_ResultUpdate



!==============================================================!
!  Subroutine TEnsemble_ResultClose                            !
!==============================================================!

  subroutine TEnsemble_ResultClose( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

#if MPI_VER > 0    
    if (SimulationType .eq. MonteCarlo) then
      ! Close running average result file
      if( .not. SimulationType .eq. SecondVirialCoeff ) call FileClose_parallel( this%iounit_runave )
  
      ! Close result file
      call FileClose_parallel( this%iounit_result )
    else
      if( .not. SimulationType .eq. SecondVirialCoeff ) call FileClose( this%iounit_runave )
  
      ! Close result file
      call FileClose( this%iounit_result )
    endif
#else
    if( .not. SimulationType .eq. SecondVirialCoeff ) call FileClose( this%iounit_runave )
  
    ! Close result file
    call FileClose( this%iounit_result )
#endif
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
    integer                   :: i, j, t!, s, o
#if  TRANS == 1
    integer                   :: k, m
    real(RK)                  :: value
    real(RK)                  :: det, inv_det
    real(RK)                  :: x(this%NComponents)
    real(RK)                  :: Inv_x(this%NComponents)
    real(RK)                  :: L(this%NComponents, this%NComponents)
    real(RK)                  :: delta(this%NComponents-1, this%NComponents-1) 
    real(RK)                  :: err_delta(this%NComponents-1, this%NComponents-1) 
    real(RK)                  :: B(this%NComponents-1, this%NComponents-1)
    real(RK)                  :: err_B(this%NComponents-1, this%NComponents-1)
    real(RK)                  :: D_12, D_13, D_14, D_23, D_24, D_34 
    real(RK)                  :: err_D12, err_D13, err_D14, err_D23, err_D24, err_D34
    
#else
#if HBOND > 0
    integer                   :: k, l
#endif
#endif
    ! Declare local variables for velocity of sound
    real(RK) :: molmass, cpid

    ! Declare local variables for phase equilibrium results
    real(RK) :: NN, yvi
    real(RK) :: dpdmu( this%NComponents ), dpdv( this%NComponents )
    real(RK) :: dydmu( this%NComponents, this%NComponents ), dydv( this%NComponents, this%NComponents )
    real(RK) :: varmu( this%NComponents ), varv( this%NComponents )
    real(RK) :: vary( this%NComponents - 1 )
    real(RK) :: AvgPressure, VarPressure, DeltaHv, VarDeltaHv
#ifdef ABL
    integer  :: counter
#endif

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
    call Error( this%SumConfEnthalpy )
    call Error( this%SumVolume )
    call Error( this%SumEPotInter )
    call Error( this%SumEPotIntra )
    if (printIDF) then
      call Error( this%SumEPotIntra_Bond )
      call Error( this%SumEPotIntra_Angle )
      call Error( this%SumEPotIntra_Dihedral )
      call Error( this%SumEPotIntra_Nonbonded )
    end if
#if OSMOP > 0
    call Error( this%SumOsmoticPressure )
    do i = 1, this%NComponents
      do j = 1, NBinsDen
        call Error( this%Component(i)%SumDenProfile(j) )
#if OSMOP == 2
        if( this%Component(i)%ChemPotMethod .eq. ChemPotMethodWidom ) then
          call Error( this%Component(i)%SumChemPotProfile(j) )
        end if
#endif
      end do
    end do
#if OSMOP == 2
    do j = 1, NBinsDen
       call Error( this%SumPressureProfile(j) )
    end do
#endif
#endif

#if HBOND > 0
    do i = 1, this%NComponents
      call Error( this%SumHBond0(i) )
      do j = 1, this%NComponents
        call Error( this%SumHBond1(i,j) )
        do k = j, this%NComponents
          call Error( this%SumHBond2(i,j,k) )
          do l = k, this%NComponents
            call Error( this%SumHBond3(i,j,k,l) )
          end do
        end do
      end do
      call Error( this%SumHBondN(i) )
    end do
#endif

    call Error( this%SumdEpotdV )
    call Error( this%Sumd2EpotdV2 )

    call Error( this%SumEPotSquared)
    call Error( this%SumEPotCubic)
    call Error( this%SumdEpotdVSquared)
    call Error( this%SumEPotdEpotdV)
    call Error( this%SumEPotSquareddEpotdV)
    call Error( this%SumEPotdEpotdVSquared)
    call Error( this%SumEPotd2EpotdV2)

    if ( LongRange .eq. Rfield ) then
      if(EnsembleType .eq. EnsembleTypeNVT) then
        call Error( this%SumA10resI )
        call Error( this%SumA01resI )
        call Error( this%SumA20resI )
        call Error( this%SumA11resI )
        call Error( this%SumA02resI )
        call Error( this%SumA30resI )
        call Error( this%SumA21resI )
        call Error( this%SumA12resI )
      elseif (EnsembleType .eq. EnsembleTypeNVE ) then
        call Error( this%SumHmU )
        call Error( this%SumHmUm1)
        call Error( this%SumHmUm2 )
        call Error( this%SumHmUm3 )
        call Error( this%SumHmUm1dUdV )
        call Error( this%SumHmUm1dUdV2 )
        call Error( this%SumHmUm1d2UdV2 )
        call Error( this%SumHmUm2dUdV )
        call Error( this%SumHmUm2dUdV2 )
        call Error( this%SumHmUm2d2UdV2 )
        call Error( this%SumHmUm3dUdV )
        call Error( this%SumHmUm3dUdV2 )

        call Error( this%SumA10resI )
        call Error( this%SumA01resI )
        call Error( this%SumA20resI )
        call Error( this%SumA11resI )
        call Error( this%SumA02resI )
        call Error( this%SumA30resI )
        call Error( this%SumA21resI )
        call Error( this%SumA12resI )
        call Error( this%SumA10resII )
        call Error( this%SumA01resII )
        call Error( this%SumA20resII )
        call Error( this%SumA11resII )
        call Error( this%SumA02resII )
        call Error( this%SumA30resII )
        call Error( this%SumA21resII )
        call Error( this%SumA12resII )
      end if
    end if

    if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
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
          call MPI_COMM_SPLIT(Communicator,color,NProc,Communicator,ierror) 
           ! Careful, Nproc and NProcs are now specific for Communicator
          call SetCommunicator( Communicator )
          NBlockSizes = int( sqrt( real( Step*NProcs / BlockSize, RK ) ) )
          NBlocks = tempVal*NProcs   
          RootProc = NProc_W==(pc%NGradThis)     
          NRootProc_W = (pc%NGradThis)     

#endif          
          
          call ErrorGI( pc%SumInvChemPotRho )
          call ErrorGI( pc%SumVW )
          call ErrorGI( pc%SumHM )
            
#if MPI_VER > 0          
          call SetCommunicator(MPI_COMM_WORLD)
          RootProc = Rootproc_W                          
#endif 

        case( ChemPotMethodWidom )
          call Error( pc%SumChemPotV )
          call Error( pc%SumHW_counter )
          call Error( pc%SumHW_denom )
          call Error( pc%SumVW )
          call Error( pc%SumHM )
        case( ChemPotMethodThermoInt )
          call Error( pc%SumChemPotV )
          call Error( pc%SumChemPotThermoIntWidom )
          call Error( pc%SumChemPotThermoIntWidomV )
          call Error( pc%SumHW_counter )
          call Error( pc%SumHW_denom )
          call Error( pc%SumVW )
          call Error( pc%SumHM )
        case default
          ! DO NOTHING
        end select
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
    call FileRewrite( this%iounit_errors, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ErrorsFileExtension )

    write( IOBuffer, '(76("="))')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("*                           Publishing with ms2                                *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* Every user agrees to cite ms2 upon usage as follows                          *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* ---------------------------------------------------------------------------- *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* C.W. Glass, S. Reiser, G. Rutkai, S. Deublein, A. Koster, G. Guevara-Carrion *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* A. Wafai, M. Horsch, M. Bernreuther, T. Windmann, H. Hasse, J. Vrabec        *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* Computer Physics Communications (2014)                                       *")')
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
    write( IOBuffer, '("Simulation type", T36, ":", 9X, A)' ) trim( SimulationTypeString )
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Ensemble type", T36, ":", 9X, A)' ) trim( EnsembleTypeString )
    call FileWrite( this%iounit_errors )
    if( SimulationType .eq. MolecularDynamics ) then
      write( IOBuffer, '("Integrator type", T36, ":", 9X, A)' ) trim( IntegratorTypeString )
      call FileWrite( this%iounit_errors )
    end if
    call FileWriteBlank( this%iounit_errors )

    ! Number of steps
    write( IOBuffer, '("Number of NVT equilibration steps", T36, ":", I10)' ) NStepsV
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Number of NVE equilibration steps", T36, ":", I10)' ) NStepsE
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Number of NPT equilibration steps", T36, ":", I10)' ) NStepsP
    call FileWrite( this%iounit_errors )
	write( IOBuffer, '("Number of NPH equilibration steps", T36, ":", I10)' ) NStepsH
    call FileWrite( this%iounit_errors )


    if ( SimulationType .eq. MonteCarlo .and. (Nproc == NRootProc)) then
      ! The RootProc receives data from all processes and therefore the # of 
      ! Step is increased accordingly
      write( IOBuffer, '("Number of production steps", T34, ":", I12)' ) Step*NProcs

    else 
      write( IOBuffer, '("Number of production steps", T36, ":", I10)' ) Step
    end if
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Time step
    if( SimulationType .eq. MolecularDynamics ) then
      write( IOBuffer, '("Time step", T29, "reduced:", F20.9)' ) TimeStep
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T31, "in fs:", F20.9)' ) TimeStep * UnitTime * 1E15_RK
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )
    end if

    ! Acceptance rate
    if( (SimulationType .eq. MonteCarlo) .or. (SimulationType .eq. Gibbs)  ) then
      write( IOBuffer, '("Acceptance rate", T36, ":", F20.9)' ) Acceptance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )
    end if

    ! Mass of piston
    if( SimulationType .eq. MolecularDynamics .and. ConstantPressure ) then
      write( IOBuffer, '("Mass of piston", T29, "reduced:", F20.9)' ) this%PistonMass
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in kg/m⁴:", F20.9)' ) this%PistonMass * 0.001_RK * UnitMass / UnitLength**4
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )
    end if

    ! Number of particles
    write( IOBuffer, '("Number of particles", T36, ":", I10)' ) this%NPart
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Potential models
    if( EnsembleType .ne. EnsembleTypeGE .or. EnsembleType .ne. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
      do i = 1, this%NRealComponents
        write( IOBuffer, '("Mole fraction of ", A, T36, " :", F20.9)' )&
&              trim( this%Component(i)%Molecule%PotModFileName ), &
&              this%Component(i)%Fraction

        call FileWrite( this%iounit_errors )
        select case( this%Component(i)%ChemPotMethod )
        case( ChemPotMethodGradIns )
          write( IOBuffer, '("Chemical potential calculated by gradual insertion")' )
          call FileWrite( this%iounit_errors )
        case( ChemPotMethodWidom )
          write( IOBuffer, '("Number of test particles", T36, ":", I10)' ) this%Component(i)%NTest
          call FileWrite( this%iounit_errors )
        end select
      end do
      call FileWriteBlank( this%iounit_errors )
    end if

    ! Initial pressure
    if( ConstantPressure ) then
      write( IOBuffer, '("Initial pressure", T29, "reduced:", F20.9)' ) this%RefPressure
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T30, "in MPa:", F20.9)' ) this%RefPressure * UnitPressure * 1E-6_RK
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )
    end if

    ! Initial density
    write( IOBuffer, '("Initial density", T29, "reduced:", F20.9)' ) this%RefDensity
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T28, "in mol/l:", F20.9)' ) this%RefDensity * UnitDensity
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Initial temperature
    write( IOBuffer, '("Initial temperature", T29, "reduced:", F20.9)' ) this%RefTemperature
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T32, "in K:", F20.9)' ) this%RefTemperature * UnitTemperature
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! System of units
    write( IOBuffer, '("Unit of length", T36, ":", F20.9, " A")' ) UnitLength / Angstroem
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Unit of energy", T36, ":", F20.9, " K")' ) UnitEnergy / kBoltzmann
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Unit of mass", T36, ":", F20.9, " a.u.")' ) UnitMass * NAvogadro * 1000._RK
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Cutoff radii
    if( this%NLJ126Max > 0 ) then
      write( IOBuffer, '("Lennard-Jones cutoff radius", T36, ":", F20.9, " A")' ) &
&            this%RCutoffLJ126LJ126 * UnitLength / Angstroem
      call FileWrite( this%iounit_errors )
    end if

    if( this%NDipoleMax > 0 ) then
      write( IOBuffer, '("Dipole-dipole cutoff radius", T36, ":", F20.9, " A")' ) &
&            this%RCutoffDipoleDipole * UnitLength / Angstroem
      call FileWrite( this%iounit_errors )

      if( this%NQuadrupoleMax > 0 ) then
        write( IOBuffer, '("Dipole-quadrupole cutoff radius", T36, ":", F20.9, " A")' ) &
&              this%RCutoffDipoleQuadrupole * UnitLength / Angstroem
        call FileWrite( this%iounit_errors )
      end if

    end if

    if( this%NQuadrupoleMax > 0 ) then
      write( IOBuffer, '("Quadrupole-quadrupole cutoff radius", T36, ":", F20.9, " A")' ) &
&            this%RCutoffQuadrupoleQuadrupole * UnitLength / Angstroem
      call FileWrite( this%iounit_errors )
    end if

    call FileWriteBlank( this%iounit_errors )

    ! Dielectric constant
    if( this%NDipoleMax > 0 ) then
      write( IOBuffer, '("Dielectric constant:", F36.9)' ) this%RFEpsilon
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
    if ( this%OptPressure ) then
      Average = this%SumPressure%Average
      Variance = this%SumPressure%Variance
    else
      Average = this%RefPressure
      Variance = 0._RK
    end if
    write( IOBuffer, '("Pressure", T29, "reduced:", 2F20.9)' ) Average, Variance
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T30, "in MPa:", 2F20.9)' ) Average * UnitPressure * 1E-6_RK, Variance * UnitPressure * 1E-6_RK
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Density
    Average = this%SumDensity%Average
    Variance = this%SumDensity%Variance
    write( IOBuffer, '("Density", T29, "reduced:", 2F20.9)' ) Average, Variance
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T28, "in mol/l:", 2F20.9)' ) Average * UnitDensity, Variance * UnitDensity
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Temperature
    Average = this%SumTemperature%Average
    Variance = this%SumTemperature%Variance
    write( IOBuffer, '("Temperature", T29, "reduced:", 2F20.9)' ) Average, Variance
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T32, "in K:", 2F20.9)' ) Average * UnitTemperature, Variance * UnitTemperature
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

#if OSMOP > 0
    if (SimulationType .eq. MolecularDynamics) then
      ! OsmoticPressure
      Average = this%SumOsmoticPressure%Average
      Variance = this%SumOsmoticPressure%Variance
      write( IOBuffer, '("OsmoticPressure", T29, "reduced:", 2F20.9)' ) &
  &     Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T30, "in MPa:", 2F20.9)' ) &
  &     Average * UnitPressure * 1E-6_RK, Variance * UnitPressure * 1E-6_RK
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )
    end if
#endif

    if (printIDF) then
      !Intermolecular potential energy
      Average = this%SumEPotInter%Average
      Variance = this%SumEPotInter%Variance
      write( IOBuffer, '("Intermolecular energy", T29, "reduced:", 2F20.9)' ) Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) Average * UnitEnergy * NAvogadro, &
&         Variance * UnitEnergy * NAvogadro
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      !Intramolecular potential energy
      Average = this%SumEPotIntra%Average
      Variance = this%SumEPotIntra%Variance
      write( IOBuffer, '("Intramolecular energy", T29, "reduced:", 2F20.9)' ) Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) Average * UnitEnergy * NAvogadro, &
&          Variance * UnitEnergy * NAvogadro
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      !Intramolecular potential energy - Bonds
      Average = this%SumEPotIntra_Bond%Average
      Variance = this%SumEPotIntra_Bond%Variance
      write( IOBuffer, '("Bond energy", T29, "reduced:", 2F20.9)' ) Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) Average * UnitEnergy * NAvogadro, &
&          Variance * UnitEnergy * NAvogadro
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      !Intramolecular potential energy - Angles
      Average = this%SumEPotIntra_Angle%Average
      Variance = this%SumEPotIntra_Angle%Variance
      write( IOBuffer, '("Angle energy", T29, "reduced:", 2F20.9)' ) Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) Average * UnitEnergy * NAvogadro, &
&          Variance * UnitEnergy * NAvogadro
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      !Intramolecular potential energy - Dihedral
      Average = this%SumEPotIntra_Dihedral%Average
      Variance = this%SumEPotIntra_Dihedral%Variance
      write( IOBuffer, '("Dihedral energy", T29, "reduced:", 2F20.9)' ) Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) Average * UnitEnergy * NAvogadro, &
&         Variance * UnitEnergy * NAvogadro
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      !Intramolecular potential energy - Nonbonded
      Average = this%SumEPotIntra_Nonbonded%Average
      Variance = this%SumEPotIntra_Nonbonded%Variance
      write( IOBuffer, '("1-4, 1-5 energy", T29, "reduced:", 2F20.9)' ) Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) Average * UnitEnergy * NAvogadro, &
&          Variance * UnitEnergy * NAvogadro
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

    end if ! printIDF

    ! Potential energy
    Average = this%SumEPot%Average
    Variance = this%SumEPot%Variance
    write( IOBuffer, '("Potential energy", T29, "reduced:", 2F20.9)' ) Average, Variance
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) Average * UnitEnergy * NAvogadro, &
&          Variance * UnitEnergy * NAvogadro
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

! Some of the ensemble averages used in calculation of Amn's (they are not necessarily candidates for .res file output)

    ! EpotSuared
!    Average = this%SumEpotSquared%Average
!    Variance = this%SumEpotSquared%Variance
!    write( IOBuffer, '("EpotSquared", T29, "reduced:", 2F20.9)' ) &
!&     Average, Variance
!    call FileWrite( this%iounit_errors )
!    write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) &
!&     Average * UnitEnergy * NAvogadro, &
!&     Variance * UnitEnergy * NAvogadro
!    call FileWrite( this%iounit_errors )
!    call FileWriteBlank( this%iounit_errors )

    ! dEpot/dV
!    Average = this%SumdEpotdV%Average
!    Variance = this%SumdEpotdV%Variance
!    write( IOBuffer, '("dEpot/dV", T29, "reduced:", 2F20.9)' ) &
!&     Average, Variance
!    call FileWrite( this%iounit_errors )
!    write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) &
!&     Average * UnitEnergy * NAvogadro, &
!&     Variance * UnitEnergy * NAvogadro
!    call FileWrite( this%iounit_errors )
!    call FileWriteBlank( this%iounit_errors )

    ! dEpot/dVSquared
!    Average = this%SumdEpotdVSquared%Average
!    Variance = this%SumdEpotdVSquared%Variance
!    write( IOBuffer, '("dEpot/dVSquared", T29, "reduced:", 2F20.9)' ) &
!&     Average, Variance
!    call FileWrite( this%iounit_errors )
!    write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) &
!&     Average * UnitEnergy * NAvogadro, &
!&     Variance * UnitEnergy * NAvogadro
!    call FileWrite( this%iounit_errors )
!    call FileWriteBlank( this%iounit_errors )

    ! d2Epot/dV2
!    Average = this%Sumd2EpotdV2%Average
!    Variance = this%Sumd2EpotdV2%Variance
!    write( IOBuffer, '("d2Epot/dV2", T29, "reduced:", 2F20.9)' ) &
!&     Average, Variance
!    call FileWrite( this%iounit_errors )
!    write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) &
!&     Average * UnitEnergy * NAvogadro, &
!&     Variance * UnitEnergy * NAvogadro
!    call FileWrite( this%iounit_errors )
!    call FileWriteBlank( this%iounit_errors )

    ! EpotdEpot/dV
!    Average = this%SumEPotdEpotdV%Average
!    Variance = this%SumEPotdEpotdV%Variance
!    write( IOBuffer, '("EpotdEpot/dV", T29, "reduced:", 2F20.9)' ) &
!&     Average, Variance
!    call FileWrite( this%iounit_errors )
!    write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) &
!&     Average * UnitEnergy * NAvogadro, &
!&     Variance * UnitEnergy * NAvogadro
!    call FileWrite( this%iounit_errors )
!    call FileWriteBlank( this%iounit_errors )

    ! Enthalpy
    Average = this%SumEnthalpy%Average
    Variance = this%SumEnthalpy%Variance
    write( IOBuffer, '("Enthalpy", T29, "reduced:", 2F20.9)' ) Average, Variance
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) Average * UnitEnergy * NAvogadro, &
&          Variance * UnitEnergy * NAvogadro
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
      ! Mole fraction
      do i = 1, this%NComponents
        pc => this%Component(i)
        Average = pc%SumFraction%Average
        Variance = pc%SumFraction%Variance
        write( IOBuffer, '("Mole fraction of ", A, T36, ":", 2F20.9)' ) &
&              trim( this%Component(i)%Molecule%PotModFileName ), Average, Variance
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
          if( pc%NPart > 1 ) then
            Average = log( pc%Fraction * pc%SumInvChemPotRho%Average )
            write( IOBuffer, '("Chem. pot. of ", A, T33, "r`d:", 2F20.9)' ) &
&                  trim( this%Component(i)%Molecule%PotModFileName ), Average, Variance
!MERKER
          else
            Average = -log( 1/pc%SumInvChemPotRho%Average )
            write( IOBuffer, '("Chem. pot. at inf. dilution of ", A, T33, "r`d:", 2F20.9)' ) &
&                  trim( this%Component(i)%Molecule%PotModFileName ), Average, Variance
            call FileWrite( this%iounit_errors )
            Average = this%Temperature*pc%SumInvChemPotRho%Average
            Variance = this%Temperature * pc%SumInvChemPotRho%Average * pc%SumInvChemPotRho%Average / pc%SumInvChemPotRho%Variance
            write( IOBuffer, '("Henrys law constant of ", A, T33, "r`d:", 2F20.9)' ) &
&                  trim( pc%Molecule%PotModFileName ), Average, Variance
            call FileWrite( this%iounit_errors )
            write( IOBuffer, '(T30, "in MPa:", 2F20.9)' ) &
&                  Average * UnitPressure * 1E-6_RK, Variance * UnitPressure * 1E-6_RK
          end if
          call FileWrite( this%iounit_errors )
!MERKER

        case( ChemPotMethodWidom )
          Variance = pc%SumChemPotV%Variance / pc%SumChemPotV%Average

          if( pc%Fraction > 0.0_RK ) then
            Average = log( pc%Fraction / pc%SumChemPotV%Average )
            write( IOBuffer, '("Chem. pot. of ", A, T33, "r`d:", 2F20.9)' ) &
&                  trim( this%Component(i)%Molecule%PotModFileName ), Average, Variance

          else
            Average = this%Temperature / pc%SumChemPotV%Average
            Variance = this%Temperature * ( pc%SumChemPotV%Variance / (pc%SumChemPotV%Average * pc%SumChemPotV%Average))
            write( IOBuffer, '("Henrys law constant of ", A, T33, "r`d:", 2F20.9)' ) &
&                  trim( pc%Molecule%PotModFileName ), Average, Variance
            call FileWrite( this%iounit_errors )
            write( IOBuffer, '(T30, "in MPa:", 2F20.9)' ) &
&                  Average * UnitPressure * 1E-6_RK, Variance * UnitPressure * 1E-6_RK
          end if
          call FileWrite( this%iounit_errors )

        case( ChemPotMethodThermoInt )
          Average  = log( (pc%Fraction+1._RK/real( this%NPart, RK )) / pc%SumChemPotThermoIntWidom%Average )
          Variance =  pc%SumChemPotThermoIntWidom%Variance / pc%SumChemPotThermoIntWidom%Average
          write( IOBuffer, '("Chem. pot. at LambdaMin ", A, T33, "r`d:", 2F20.9)' ) &
&                trim( this%Component(i)%Molecule%PotModFileName ), Average , Variance
          call FileWrite( this%iounit_errors )
          Average  = pc%SumChemPotV%Average
          Variance = pc%SumChemPotV%Variance
          write( IOBuffer, '("Chem. pot. of ", A, T33, "r`d:", 2F20.9)' ) &
&                trim( this%Component(i)%Molecule%PotModFileName ), Average, Variance
          call FileWrite( this%iounit_errors )

          if( pc%Npart .eq. 0 ) then
            ! Actually: Average  = this%Temperature * exp(pc%SumChemPotV%Average)/(pc%Fraction+1._RK/real( this%NPart, RK )), but pc%Fraction=0.0
            Average  = this%Temperature * exp(pc%SumChemPotV%Average)*this%NPart
            Variance = Average * pc%SumChemPotV%Variance
            write( IOBuffer, '("Henrys law constant of ", A, T33, "r`d:", 2F20.9)' ) &
&                  trim( pc%Molecule%PotModFileName ), Average, Variance
            call FileWrite( this%iounit_errors )
            write( IOBuffer, '(T30, "in MPa:", 2F20.9)' ) &
&                  Average * UnitPressure * 1E-6_RK, Variance * UnitPressure * 1E-6_RK
            call FileWrite( this%iounit_errors )
          end if

        end select
      end do

      if( any(this%Component(:)%ChemPotMethod .ne. ChemPotMethodNone)) call FileWriteBlank( this%iounit_errors )

      ! Partial molar volume
      do i = 1, this%NRealComponents
        pc => this%Component(i)
        if( pc%ChemPotMethod .ne. ChemPotMethodNone .and. EnsembleType .eq. EnsembleTypeNPT) then
          Average = pc%SumVW%Average
          Variance = pc%SumVW%Variance
          write( IOBuffer, '("Partial molar volume of ", A, T33, "r`d:", 2F20.9)' ) &
&                trim( this%Component(i)%Molecule%PotModFileName ), Average, Variance
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T28, "in l/mol:", 2F20.9)' ) Average / UnitDensity, Variance / UnitDensity
          call FileWrite( this%iounit_errors ) 

          ! Partial molar enthalpy
          Average = pc%SumHM%Average
          Variance = pc%SumHM%Variance
          write( IOBuffer, '("Partial molar enthalpy of ", A, T33, "r`d:", 2F20.9)' ) &
&                trim( this%Component(i)%Molecule%PotModFileName ), Average, Variance
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) Average * UnitEnergy * NAvogadro, &
&                Variance * UnitEnergy * NAvogadro
          call FileWrite( this%iounit_errors )
        end if
      end do
      if( any(this%Component(:)%ChemPotMethod .ne. ChemPotMethodNone)) call FileWriteBlank( this%iounit_errors )

      if( ConstantPressure ) then
        ! Isothermal compressibility
        Average = this%SumBetaT%Average
        Variance = this%SumBetaT%Variance
        write( IOBuffer, '("Isothermal compressibility", T29, "reduced:", 2F20.9)' ) Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T28, "in 1/MPa:", 2F20.9)' ) Average / ( UnitPressure * 1E-6_RK ), &
&              Variance / ( UnitPressure * 1E-6_RK )
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )

        ! dH/dP
        Average = this%SumdHdP%Average
        Variance = this%SumdHdP%Variance
        write( IOBuffer, '("dH/dP", T29, "reduced:", 2F20.9)' ) Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T28, "in l/mol:", 2F20.9)' ) Average / UnitDensity, Variance / UnitDensity
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )

        ! CP - subtract ideal gas contribution of the pressure
        Average = this%SumCP%Average - 1._RK
        Variance = this%SumCP%Variance
        write( IOBuffer, '("Isobaric heat capacity", T29, "reduced:", 2F20.9)' ) Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T24, "in J/(mol K):", 2F20.9)' ) Average * kBoltzmann * NAvogadro, &
&              Variance * kBoltzmann * NAvogadro
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )

        ! AlphaP
        Average = this%SumAlphaP%Average
        Variance = this%SumAlphaP%Variance
        write( IOBuffer, '("Volume expansivity", T29, "reduced:", 2F20.9)' ) Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T30, "in 1/K:", 2F20.9)' ) Average / UnitTemperature, Variance / UnitTemperature
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

        Average = 1._RK / sqrt( molmass * ( this%SumBetaT%Average * this%SumDensity%Average - this%RefTemperature * &
&                 this%SumAlphaP%Average**2 / ( this%SumCP%Average + cpid ) ) )

        Variance = .25_RK / molmass / ( this%SumBetaT%Average * this%SumDensity%Average - this%RefTemperature * &
&                  this%SumAlphaP%Average**2 / ( this%SumCP%Average + cpid ) )**3 * ( this%SumDensity%Average**2 *&
&                  this%SumBetaT%Variance**2 + this%SumBetaT%Average**2 * this%SumDensity%Variance**2 +&
&                  this%RefTemperature**2 * this%SumAlphaP%Average**2 / ( this%SumCP%Average + cpid )**2 *&
&                  ( 4._RK * this%SumAlphaP%Variance**2 + this%SumAlphaP%Average**2 / ( this%SumCP%Average + cpid )**2 * &
&                  this%SumCP%Variance**2 ) )

        write( IOBuffer, '("Speed of sound", T29, "reduced:", 2F20.9)' ) Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T30, "in m/s:", 2F20.9)' ) Average * UnitLength / UnitTime, Variance * UnitLength / UnitTime
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )

      else
        ! dU/dV
        Average = this%SumdUdV%Average
        Variance = this%SumdUdV%Variance
        write( IOBuffer, '("dU/dV", T29, "reduced:", 2F20.9)' ) Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T30, "in MPa:", 2F20.9)' ) Average * UnitPressure * 1E-6_RK, Variance * UnitPressure * 1E-6_RK
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )

        ! Cv
        Average = this%SumCV%Average
        Variance = this%SumCV%Variance
        write( IOBuffer, '("Isochoric heat capacity", T29, "reduced:", 2F20.9)' ) Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T24, "in J/(mol K):", 2F20.9)' ) Average * kBoltzmann * NAvogadro, &
&              Variance * kBoltzmann * NAvogadro
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )
      endif

    end if

    if( EnsembleType .eq. EnsembleTypeNVT .and. LongRange .eq. Rfield ) then 
      ! A10
      Average = this%SumA10resI%Average
      Variance = this%SumA10resI%Variance
      write( IOBuffer, '("A10 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A01
      Average = this%SumA01resI%Average
      Variance = this%SumA01resI%Variance
      write( IOBuffer, '("A01 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A20
      Average = this%SumA20resI%Average
      Variance = this%SumA20resI%Variance
      write( IOBuffer, '("A20 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A11
      Average = this%SumA11resI%Average
      Variance = this%SumA11resI%Variance
      write( IOBuffer, '("A11 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A02
      Average = this%SumA02resI%Average
      Variance = this%SumA02resI%Variance
      write( IOBuffer, '("A02 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A30
      Average = this%SumA30resI%Average
      Variance = this%SumA30resI%Variance
      write( IOBuffer, '("A30 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A21
      Average = this%SumA21resI%Average
      Variance = this%SumA21resI%Variance
      write( IOBuffer, '("A21 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A12
      Average = this%SumA12resI%Average
      Variance = this%SumA12resI%Variance
      write( IOBuffer, '("A12 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )
    end if

    if( EnsembleType .eq. EnsembleTypeNVE .and. LongRange .eq. Rfield ) then 
      ! A10I
      Average = this%SumA10resI%Average
      Variance = this%SumA10resI%Variance
      write( IOBuffer, '("A10 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A01I
      Average = this%SumA01resI%Average
      Variance = this%SumA01resI%Variance
      write( IOBuffer, '("A01 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A20I
      Average = this%SumA20resI%Average
      Variance = this%SumA20resI%Variance
      write( IOBuffer, '("A20 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A11I
      Average = this%SumA11resI%Average
      Variance = this%SumA11resI%Variance
      write( IOBuffer, '("A11 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A02I
      Average = this%SumA02resI%Average
      Variance = this%SumA02resI%Variance
      write( IOBuffer, '("A02 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A30I
      Average = this%SumA30resI%Average
      Variance = this%SumA30resI%Variance
      write( IOBuffer, '("A30 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A21I
      Average = this%SumA21resI%Average
      Variance = this%SumA21resI%Variance
      write( IOBuffer, '("A21 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A12I
      Average = this%SumA12resI%Average
      Variance = this%SumA12resI%Variance
      write( IOBuffer, '("A12 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A10II
      Average = this%SumA10resII%Average
      Variance = this%SumA10resII%Variance
      write( IOBuffer, '("A10 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A01II
      Average = this%SumA01resII%Average
      Variance = this%SumA01resII%Variance
      write( IOBuffer, '("A01 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A20II
      Average = this%SumA20resII%Average
      Variance = this%SumA20resII%Variance
      write( IOBuffer, '("A20 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A11II
      Average = this%SumA11resII%Average
      Variance = this%SumA11resII%Variance
      write( IOBuffer, '("A11 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A02II
      Average = this%SumA02resII%Average
      Variance = this%SumA02resII%Variance
      write( IOBuffer, '("A02 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A30II
      Average = this%SumA30resII%Average
      Variance = this%SumA30resII%Variance
      write( IOBuffer, '("A30 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A21II
      Average = this%SumA21resII%Average
      Variance = this%SumA21resII%Variance
      write( IOBuffer, '("A21 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! A12II
      Average = this%SumA12resII%Average
      Variance = this%SumA12resII%Variance
      write( IOBuffer, '("A12 - Dimensionless, residual", T36,":", 2F20.9)' ) &
&       Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )
    end if

    ! Separator
    write( IOBuffer, '(76("="))' )
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

#if HBOND > 0
    do i = 1, this%NComponents
      Average = this%SumHBond0(i)%Average
      Variance = this%SumHBond0(i)%Variance
      write( IOBuffer, '("HBond0 of [", I2, "]", T36, ":", 2F20.9)' ) i, Average, Variance
      call FileWrite( this%iounit_errors )
    end do
    do i = 1, this%NComponents
      do  j = 1, this%NComponents
        Average = this%SumHBond1(i,j)%Average
        Variance = this%SumHBond1(i,j)%Variance
        write( IOBuffer, '("HBond1 of [", I2, "] with (", I2, ")", T36, ":", 2F20.9)' ) i, j, Average, Variance
        call FileWrite( this%iounit_errors )
      end do
    end do
    do i = 1, this%NComponents
      do  j = 1, this%NComponents
        do k = j, this%NComponents
          Average = this%SumHBond2(i,j,k)%Average
          Variance = this%SumHBond2(i,j,k)%Variance
          write( IOBuffer, '("HBond2 of [", I2, "] with (", I2, ",", I2, ")", T36, ":", 2F20.9)' ) i, j, k, Average, Variance
          call FileWrite( this%iounit_errors )
        end do
      end do
    end do
    do i = 1, this%NComponents
      do  j = 1, this%NComponents
        do k = j, this%NComponents
          do  l = k, this%NComponents
            Average = this%SumHBond3(i,j,k,l)%Average
            Variance = this%SumHBond3(i,j,k,l)%Variance
            write( IOBuffer, '("HBond3 of [", I2, "] with (", I2, ",", I2, ",", I2, ")", T36, ":", 2F20.9)' ) i, j, k, l, Average, Variance
            call FileWrite( this%iounit_errors )
          end do
        end do
      end do
    end do
    do i = 1, this%NComponents
      Average = this%SumHBondN(i)%Average
      Variance = this%SumHBondN(i)%Variance
      write( IOBuffer, '("HBond4+ of [", I2, "]", T36, ":", 2F20.9)' ) i, Average, Variance
      call FileWrite( this%iounit_errors )
    end do
    call FileWriteBlank( this%iounit_errors )

    ! Separator
    write( IOBuffer, '(76("="))' )
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )
#endif

#if  TRANS == 1
    ! Transport properties !TRANSPORT_start
    if ( this%CorrfunMode ) then

      write( IOBuffer, '(T24, "TRANSPORT PROPERTIES")' )
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      write( IOBuffer, '("VALUE", T31, "UNITS", T46, "AVERAGE", T66, "ERROR")' )
      call FileWrite( this%iounit_errors )

      write( IOBuffer, '("-----", T31, "-----", T46, "-------", T66, "-----")' )
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      write( IOBuffer, '("Number of ACF", T36, ":",T45, I6 )' ) this%Mmess
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      value = this%NCorr*this%TimeStepCorr
      write( IOBuffer, '("Length ACF  ", T29, "reduced:", F20.9)' ) value
      call FileWrite( this%iounit_errors )

      write( IOBuffer, '(T31, "in ps:", F20.9)' )  value*UnitTime/1E-12_RK
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      value = this%NSpanCF*this%TimeStepCorr
      write( IOBuffer, '("Time span between ACF ", T29, "reduced:", F20.9)' ) value
      call FileWrite( this%iounit_errors )

      write( IOBuffer, '(T31, "in ps:", F20.9)' )  value*UnitTime/1E-12_RK
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

        
      if ( this%Mmess > 0 ) then

        ! Error update for transport properties
        if( mod( ((Step-1)/this%NStepCorr) - this%NCorr + 1, BlockSizeCF * this%NSpanCF ) == 0 .and. NBlockSizesCF >= 2) then

          if ( this%NComponents > 1 ) then
            do i = 1, this%NComponents
              do j = 1, this%NComponents
                call Error(this%SumOnsager(i,j), .true.)
              end do
            end do 
          end if
          if( this%NComponents == 2  ) then
            if (this%MolarEnthConduct .eqv. .true.) then
              call Error(this%SumSoret, .true.)
            end if
          end if
          do i = 1, this%NComponents
            call Error(this%Sumself_i(i), .true.)
          end do
          call Error(this%SumVisco_s, .true.)
          call Error(this%SumVisco_b, .true.)
          call Error(this%SumConduct, .true.)
          call Error(this%SumEConduct,.true.)

        end if

        ! Onsager coefficients
        if ( this%NComponents > 1 ) then
          do i = 1, this%NComponents
            do j = 1, this%NComponents
              Average  = this%SumOnsager(i,j)%Average
              Variance = this%SumOnsager(i,j)%Variance
              value = dsqrt(UnitEnergy/UnitMass)*UnitLength/1E-10_RK
              write( IOBuffer, '("Onsager-diff. coeff.",2I2,T29, "reduced:", 2F20.9)' ) i,j,Average, Variance
              call FileWrite( this%iounit_errors )
              write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) Average*value, Variance*value
              call FileWrite( this%iounit_errors )     
            end do
          end do 
          call FileWriteBlank( this%iounit_errors )
        end if !this%NComponents

        !for multicomponent mixtures

        if( this%NComponents >= 2  ) then
 
           do i = 1, this%NComponents
              x(i) = this%Component(i)%Fraction
              Inv_x(i) = 1._RK/x(i)
           end do

           do i = 1, this%NComponents
              do  j = 1, this%NComponents
                 if (i ==j) then
                  L(i,j) = this%SumOnsager(i,j)%Average
                 else
                  L(i,j) = (this%SumOnsager(i,j)%Average + this%SumOnsager(j,i)%Average)/2._RK
                 end if
              end do
           end do
            
        end if          

        
        
        
        !binary diffusion and thermal diffusion
        if( this%NComponents == 2  ) then
          value = dsqrt(UnitEnergy/UnitMass)*UnitLength/1E-10_RK
          
          D_12 = L(1,1) * x(2)*Inv_x(1) + L(2,2) * x(1)*Inv_x(2) - L(1,2) - L(2,1)
          err_D12 = this%SumOnsager(1,1)%Variance * x(2)*Inv_x(1) + &
&                   this%SumOnsager(2,2)%Variance * x(1)* Inv_x(2) + &
&                   this%SumOnsager(1,2)%Variance + this%SumOnsager(2,1)%Variance
          
          write( IOBuffer, '("Binary diff. coeff.", T29, "reduced:", 2F20.9)' ) D_12, err_D12
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) D_12*value, err_D12*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )

          if (this%MolarEnthConduct .eqv. .true.) then
            Average  = this%SumSoret%Average
            Variance = this%SumSoret%Variance
            value = dsqrt(UnitEnergy/UnitMass)*UnitLength*(kBoltzmann/UnitEnergy)/1E-12_RK
            write( IOBuffer, '("Thermal diff. coeff",A, T29, "reduced:", 2F20.9)' ) trim(this%Component(2)%Molecule%PotModFileName), Average, Variance
            call FileWrite( this%iounit_errors )
            write( IOBuffer, '(T17, "in 10E-12 m^2/(K s):", 2F20.9)' ) Average*value, Variance*value
          else
            write( IOBuffer, '("Thermal diffusivity requires the partial molar enthalpies of all components")' )
          end if  !this%MolarEnthConduct
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )

        end if !this components = 2
     
     
        ! Ternary and Quaternary diffusion
        if(( this%NComponents == 3 ) .or. ( this%NComponents == 4 )) then

          !obtain matrix [delta] Equations 48 to 55 from Supplementary material 
          !Krishna and van Baten, Ind. Eng. Chem. Res., 2005, 44 (17), pp 6939

          delta(:,:) = 0._RK
          do i=1, (this%NComponents-1)
             do j =1, (this%NComponents-1)
                delta(i,j) = (1._RK-x(i))*(L(i,j)*Inv_x(j)-L(i,this%NComponents)*Inv_x(this%NComponents))
                do k = 1, this%NComponents
                   if (k /= i) then
                    delta(i,j) = delta(i,j) - x(i)* (L(k,j)*Inv_x(j)-L(k,this%NComponents)*Inv_x(this%NComponents))
                   end if
                end do
             end do
          end do

         !calculate variance by error propagation
          err_delta(:,:) = 0._RK
           do i=1, (this%NComponents-1)
             do j =1, (this%NComponents-1)
                err_delta(i,j) = (1._RK-x(i))*Inv_x(j)*this%SumOnsager(i,j)%Variance + (1._RK-x(i))*Inv_x(this%NComponents)*this%SumOnsager(i,this%NComponents)%Variance
                do k = 1, this%NComponents
                   if (k /= i) then
                    err_delta(i,j) = err_delta(i,j) + x(i)*Inv_x(j)*this%SumOnsager(k,j)%Variance + x(i)*Inv_x(this%NComponents)*this%SumOnsager(k,this%NComponents)%Variance
                   end if
                end do
             end do
          end do

        end if 


     
        !Ternary diffusion
        if( this%NComponents == 3 ) then
          value = dsqrt(UnitEnergy/UnitMass)*UnitLength/1E-10_RK
  
          ! determinat of matrix [delta]
          det = (delta(1,1)*delta(2,2))-(delta(1,2)*delta(2,1))
          inv_det = 1._RK/det

          !obtain matrix [B] so that [B]=[D]-1
          B(1,1) =  inv_det* delta(2,2) !B1
          B(1,2) =  inv_det*(-delta(1,2)) !B2
          B(2,1) =  inv_det*(-delta(2,1)) !B3
          B(2,2) =  inv_det* delta(1,1) !B4

          !Obtain Error matrix B (from Propagation of Errors for Matrix Inversion, 
          !Lefebvre et al.Nucl.Instrm.Meth. A451 (2000) 520-528)
          
           err_B(:,:) = 0._RK

           do k = 1, (this%NComponents-1)
              do m = 1, (this%NComponents-1)
                 do i = 1, (this%NComponents-1)
                    do  j = 1, (this%NComponents-1)
                       err_B(k,m) = err_B(k,m) + ABS(B(k,i)*B(j,m))*err_delta(i,j)
                    end do
                 end do
               end do
           end do
 
          !Calculate diffusion coefficients
          D_13 =  1._RK  / ( (B(1,1)) + ( x(2)* B(1,2) * Inv_x(1)) )            
          D_12 =  1._RK  / ( (B(1,1)) - ( (x(1) + x(3)) * B(1,2) *Inv_x(1)))  
          D_23 =  1._RK  / ( (B(2,2)) + ( x(1)* B(2,1) * Inv_x(2)))         
 
          !Obtain error of Diffusion coefficients
          err_D13 = ABS(1._RK/((x(2)*Inv_x(1)*B(1,2)+B(1,1))**2))*err_B(1,1) + &
&                   ABS(x(2)*Inv_x(1)/((B(1,1)+x(2)*Inv_x(1)*B(1,2))**2))*err_B(1,2)
          err_D12 = ABS(1._RK/((B(1,1)-((x(1)+x(3))*Inv_x(1)*B(1,2)))**2))*err_B(1,1) + &
&                   ABS(((x(1)+x(3))*Inv_x(1))/((B(1,1)-((x(1)+x(3))*Inv_x(1)*B(1,2)))**2))*err_B(1,2)
          err_D23 = ABS(1._RK/((x(1)*Inv_x(2)*B(2,1)+B(2,2))**2))*err_B(2,2) + &
&                   ABS(x(1)*Inv_x(2)/((B(2,2)+x(1)*Inv_x(2)*B(2,1))**2))*err_B(2,1)

          write( IOBuffer, '("Ternary diff. coeff. 1 2", T29, "reduced:", 2F20.9)' ) D_12, err_D12 
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) D_12*value, err_D12*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
          write( IOBuffer, '("Ternary diff. coeff. 1 3", T29, "reduced:", 2F20.9)' ) D_13, err_D13 
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) D_13*value, err_D13*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )      
          write( IOBuffer, '("Ternary diff. coeff. 2 3", T29, "reduced:", 2F20.9)' ) D_23, err_D23
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) D_23*value, err_D23*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
        end if !this%NComponents == 3 

        
        if ( this%NComponents == 4 ) then

          value = dsqrt(UnitEnergy/UnitMass)*UnitLength/1E-10_RK

          ! determinat of matrix [delta]
          det = (delta(1,1)*delta(2,2)*delta(3,3))+(delta(2,1)*delta(3,2)*delta(1,3))+(delta(3,1)*delta(1,2)*delta(2,3))-&
&               (delta(1,1)*delta(3,2)*delta(2,3))-(delta(3,1)*delta(2,2)*delta(1,3))-(delta(2,1)*delta(1,2)*delta(3,3))

          inv_det = 1._RK/det

          !obtain matrix [B] so that [B]=[D]-1
          B(1,1) =  inv_det* (delta(2,2)*delta(3,3)-delta(2,3)*delta(3,2)) !B1
          B(1,2) =  inv_det* (delta(1,3)*delta(3,2)-delta(1,2)*delta(3,3)) !B2
          B(1,3) =  inv_det* (delta(1,2)*delta(2,3)-delta(1,3)*delta(2,2))
          B(2,1) =  inv_det* (delta(2,3)*delta(3,1)-delta(2,1)*delta(3,3))
          B(2,2) =  inv_det* (delta(1,1)*delta(3,3)-delta(1,3)*delta(3,1))
          B(2,3) =  inv_det* (delta(1,3)*delta(2,1)-delta(1,1)*delta(2,3))
          B(3,1) =  inv_det* (delta(2,1)*delta(3,2)-delta(2,2)*delta(3,1))
          B(3,2) =  inv_det* (delta(1,2)*delta(3,1)-delta(1,1)*delta(3,2))
          B(3,3) =  inv_det* (delta(1,1)*delta(2,2)-delta(1,2)*delta(2,1)) 

          !Obtain Error matrix B (from Propagation of Errors for Matrix Inversion, 
          !Lefebvre et al.Nucl.Instrm.Meth. A451 (2000) 520-528)
          
          err_B(:,:) = 0._RK

           do k = 1, (this%NComponents-1)  
              do m = 1, (this%NComponents-1)  
                 do i = 1, (this%NComponents-1)
                    do  j = 1, (this%NComponents-1)
                       err_B(k,m) = err_B(k,m) + ABS(B(k,i)*B(j,m))*err_delta(i,j)
                    end do
                 end do
               end do
           end do


          !Calculate diffusion coefficients
          D_14 =  1._RK  / ( (B(1,1)) + ( x(2)* B(1,2) * Inv_x(1)) + (x(3) * B(1,3)* Inv_x(3)) )         
          D_24 =  1._RK  / ( (B(2,2)) + ( x(1)* B(2,1) * Inv_x(2)) + (x(3) * B(2,3)* Inv_x(2)) )        
          D_34 =  1._RK  / ( (B(3,3)) + ( x(1)* B(3,1) * Inv_x(3)) + (x(2) * B(3,2)* Inv_x(3)) )       
          D_12 =  1._RK  / ( (1._RK/D_24) - (B(2,1)*Inv_x(2)))
          D_13 =  1._RK  / ( (1._RK/D_14) - (B(1,3)*Inv_x(1)))
          D_23 =  1._RK  / ( (1._RK/D_24) - (B(2,3)*Inv_x(2)))


          !Obtain error of Diffusion coefficients
          err_D14 = ABS(1._RK/(((B(1,1)) + ( x(2)* B(1,2) * Inv_x(1)) + (x(3) * B(1,3)* Inv_x(3)) )**2))*err_B(1,1) + &
                    ABS(x(2)*Inv_x(1)/(((B(1,1)) + ( x(2)* B(1,2) * Inv_x(1)) + (x(3) * B(1,3)* Inv_x(3)) )**2))*err_B(1,2) + &
                    ABS(x(3)*Inv_x(1)/(((B(1,1)) + ( x(2)* B(1,2) * Inv_x(1)) + (x(3) * B(1,3)* Inv_x(3)) )**2))*err_B(1,3)

          err_D24 = ABS(1._RK/(((B(2,2)) + ( x(1)* B(2,1) * Inv_x(2)) + (x(3) * B(2,3)* Inv_x(2)) )**2))*err_B(2,2) + &
                    ABS(x(1)*Inv_x(2)/(((B(2,2)) + ( x(1)* B(2,1) * Inv_x(2)) + (x(3) * B(2,3)* Inv_x(2)) )**2))*err_B(2,1) + &
                    ABS(x(3)*Inv_x(2)/(((B(2,2)) + ( x(1)* B(2,1) * Inv_x(2)) + (x(3) * B(2,3)* Inv_x(2)) )**2))*err_B(2,3)

          err_D34 = ABS(1._RK/(((B(3,3)) + ( x(1)* B(3,1) * Inv_x(3)) + (x(2) * B(3,2)* Inv_x(3)) )**2))*err_B(3,3) + &
                    ABS(x(1)*Inv_x(3)/(((B(3,3)) + ( x(1)* B(3,1) * Inv_x(3)) + (x(2) * B(3,2)* Inv_x(3)) )**2))*err_B(3,1) + &
                    ABS(x(2)*Inv_x(3)/(((B(3,3)) + ( x(1)* B(3,1) * Inv_x(3)) + (x(2) * B(3,2)* Inv_x(3)) )**2))*err_B(3,2)

          err_D12 = ABS(1._RK/(((B(2,2)) + ( (x(1)-1._RK)* B(2,1) * Inv_x(2)) + (x(3) * B(2,3)* Inv_x(2)) )**2))*err_B(2,2) + &
                    ABS((x(1)-1._RK)*Inv_x(2)/(((B(2,2)) + ( (x(1)-1._RK)* B(2,1) * Inv_x(2)) + (x(3) * B(2,3)* Inv_x(2)) )**2))*err_B(2,1) + &
                    ABS(x(3)*Inv_x(2)/(((B(2,2)) + ( (x(1)-1._RK)* B(2,1) * Inv_x(2)) + (x(3) * B(2,3)* Inv_x(2)) )**2))*err_B(2,3)

          err_D13 = ABS(1._RK/(((B(1,1)) + ( x(2)* B(1,2) * Inv_x(1)) + ((x(3)-1._RK) * B(1,3)* Inv_x(3)) )**2))*err_B(1,1) + &
                    ABS(x(2)*Inv_x(1)/(((B(1,1)) + ( x(2)* B(1,2) * Inv_x(1)) + ((x(3)-1._RK) * B(1,3)* Inv_x(3)) )**2))*err_B(1,2) + &
                    ABS((x(3)-1._RK)*Inv_x(1)/(((B(1,1)) + ( x(2)* B(1,2) * Inv_x(1)) + ((x(3)-1._RK) * B(1,3)* Inv_x(3)) )**2))*err_B(1,3)

          err_D23 = ABS(1._RK/(((B(2,2)) + ( x(1)* B(2,1) * Inv_x(2)) + ((x(3)-1._RK) * B(2,3)* Inv_x(2)) )**2))*err_B(2,2) + &
                    ABS(x(1)*Inv_x(2)/(((B(2,2)) + ( x(1)* B(2,1) * Inv_x(2)) + ((x(3)-1._RK) * B(2,3)* Inv_x(2)) )**2))*err_B(2,1) + &
                    ABS((x(3)-1._RK)*Inv_x(2)/(((B(2,2)) + ( x(1)* B(2,1) * Inv_x(2)) + ((x(3)-1._RK) * B(2,3)* Inv_x(2)) )**2))*err_B(2,3)


          write( IOBuffer, '("Quat. diff. coeff. 1 2", T29, "reduced:", 2F20.9)' ) D_12, err_D12
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) D_12*value, err_D12*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
          write( IOBuffer, '("Quat. diff. coeff. 1 3", T29, "reduced:", 2F20.9)' ) D_13, err_D13
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) D_13*value, err_D13*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
          write( IOBuffer, '("Quat. diff. coeff. 1 4", T29, "reduced:", 2F20.9)' ) D_14, err_D14
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) D_14*value, err_D14*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
          write( IOBuffer, '("Quat. diff. coeff. 2 3", T29, "reduced:", 2F20.9)' ) D_23, err_D23
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) D_23*value, err_D23*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
          write( IOBuffer, '("Quat. diff. coeff. 2 4", T29, "reduced:", 2F20.9)' ) D_24, err_D24
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) D_24*value, err_D24*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
          write( IOBuffer, '("Quat. diff. coeff. 3 4", T29, "reduced:", 2F20.9)' ) D_34, err_D34
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) D_34*value, err_D34*value
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )

        end if !this%NComponents == 4


        !self-diffusion coefficient
        do i = 1, this%NComponents
          Average  = this%Sumself_i(i)%Average
          Variance = this%Sumself_i(i)%Variance
          value = dsqrt(UnitEnergy/UnitMass)*UnitLength/1E-10_RK
          write( IOBuffer, '("Self-diff. coeff.",A ,T29, "reduced:", 2F20.9)' )  &
&                trim( this%Component(i)%Molecule%PotModFileName ), Average, Variance
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) Average*value, Variance*value
          call FileWrite( this%iounit_errors )
        end do
        call FileWriteBlank( this%iounit_errors )

        !shear viscosity
        Average  = this%SumVisco_s%Average
        Variance = this%SumVisco_s%Variance
        value = dsqrt(UnitEnergy*UnitMass)/UnitLength**2/1E-4_RK
        write( IOBuffer, '("Shear-Viscosity    ", T29, "reduced:", 2F20.9)' ) Average, Variance
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T23, "in 10E-4 Pa s:", 2F20.9)' ) Average*value, Variance*value
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )

        !bulk viscosity
        if (this%Bulkviscosity ) then
          Average  = this%SumVisco_b%Average
           Variance = this%SumVisco_b%Variance
          write( IOBuffer, '("Bulk-Viscosity    ", T29, "reduced:", 2F20.9)' ) Average, Variance
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T23, "in 10E-4 Pa s:", 2F20.9)' ) Average*value, Variance*value
        else
          write( IOBuffer, '("Bulk viscosity only defined for the NVE ensemble")' )
        end if
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )

        !Thermal conductivity
        Average  = this%SumConduct%Average
        Variance = this%SumConduct%Variance
        value = dsqrt(UnitEnergy/UnitMass)*kBoltzmann/UnitLength**2
        if (LongRange .eq. Ewald) then
          write( IOBuffer, '("Thermal conductivity just implemented for reaction field")' )
        elseif (this%NComponents==1 .or. this%MolarEnthConduct .eqv. .true.) then
           write( IOBuffer, '("Thermal conductivity ", T29, "reduced:", 2F20.9)' ) Average, Variance
           call FileWrite( this%iounit_errors )
           write( IOBuffer, '(T23, "in W / (m K) :", 2F20.9)' ) Average*value, Variance*value
        else
          write( IOBuffer, '("Thermal conductivity requires the partial molar enthalpies of all components")' )
        end if
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )

        !Electric conductivity
        Average  = this%SumEConduct%Average
        Variance = this%SumEConduct%Variance
        value = ElementaryCharge**2 / (dsqrt(UnitEnergy*UnitMass) * UnitLength**2)
        if (this%EConductivity) then
          write( IOBuffer, '("Electric conductivity ", T29, "reduced:", 2F20.9)' ) Average, Variance
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T23, "in 1 / (Ohm m):", 2F20.9)' ) Average*value, Variance*value
        else
          write( IOBuffer, '("Electric conductivity only defined for charged particles")' )
        end if
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )

      else    ! ( this%Mmess > 0 )

         !Onsager coefficients
        if ( this%NComponents > 1 ) then
           do i = 1, this%NComponents
              do j = 1, this%NComponents
                 write( IOBuffer, '("Onsager-diff. coeff.",2I2,T29, "reduced:", 2F20.9)' ) i,j,0._RK
                 call FileWrite( this%iounit_errors )
                 write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) 0._RK
                 call FileWrite( this%iounit_errors )
              end do
           end do
           call FileWriteBlank( this%iounit_errors )
        end if

        if ( this%NComponents==2 ) then
          write( IOBuffer, '("Binary diff. coeff.", T29, "reduced:", F20.9)' ) 0._RK
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", F20.9)' )  0._RK
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
          if (this%MolarEnthConduct .eqv. .true.) then
            write( IOBuffer, '("Thermal diff. coeff.", A, T29, "reduced:", F20.9)' ) trim(this%Component(2)%Molecule%PotModFileName), 0._RK
            call FileWrite( this%iounit_errors )
            write( IOBuffer, '(T21, "in 10E-12 m^2/(K s):", F20.9)' ) 0._RK 
          else
            write( IOBuffer, '("Thermal diffusivity requires the partial molar enthalpies of all components")' )
          end if
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
          
        end if

         !ternary diffusion coefficient
        if( this%NComponents == 3 ) then
          write( IOBuffer, '("Ternary diff. coeff. 1 3", T29, "reduced:", 2F20.9)') 0._RK
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) 0._RK
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
          write( IOBuffer, '("Ternary diff. coeff. 1 2", T29, "reduced:", 2F20.9)' ) 0._RK
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) 0._RK
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
          write( IOBuffer, '("Ternary diff. coeff. 2 3", T29, "reduced:", 2F20.9)' ) 0._RK
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", 2F20.9)' ) 0._RK
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
        end if    

        do i = 1, this%NComponents
          write( IOBuffer, '("Self-diff. coeff.",A ,T29, "reduced:", F20.9)' ) trim( this%Component(i)%Molecule%PotModFileName ), 0._RK
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T21, "in 10E-10 m^2/s:", F20.9)' )  0._RK
          call FileWrite( this%iounit_errors )
        end do
        call FileWriteBlank( this%iounit_errors )

        write( IOBuffer, '("Shear-Viscosity    ", T29, "reduced:", F20.9)' )  0._RK
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T23, "in 10E-4 Pa s:", F20.9)' ) 0._RK
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )

        if (this%Bulkviscosity ) then
          write( IOBuffer, '("Bulk-Viscosity     ", T29, "reduced:", F20.9)' )  0._RK
          call FileWrite( this%iounit_errors )
          write( IOBuffer, '(T23, "in 10E-4 Pa s:", F20.9)' ) 0._RK
        else
          write( IOBuffer, '("Bulk viscosity only defined for the NVE ensemble")' )
        end if
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )

        if (LongRange .eq. Ewald) then
          write( IOBuffer, '("Thermal conductivity just implemented for reaction field")' )
        elseif (this%NComponents==1 .or. this%MolarEnthConduct .eqv. .true.) then
           write( IOBuffer, '("Thermal conductivity ", T29, "reduced:", 2F20.9)' ) 0._RK
           call FileWrite( this%iounit_errors )
           write( IOBuffer, '(T23, "in W / (m K) :", 2F20.9)' ) 0._RK
        else
          write( IOBuffer, '("Thermal conductivity requires the partial molar enthalpies of all components")' )
        end if
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )
      

        if (this%EConductivity) then
           write( IOBuffer, '("Electric Conductivity ", T29, "reduced:", F20.9)' )  0._RK
           call FileWrite( this%iounit_errors )
           write( IOBuffer, '(T23, "in 1 / (Ohm m):", F20.9)' ) 0._RK
        else
          write( IOBuffer, '("Electric conductivity only defined for pure charged particles")' )
        end if
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )

      end if
!TRANSPORT_END
      ! Separator
      write( IOBuffer, '(76("="))' )
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

    end if   ! (this%CorrfunMode)
#endif

    ! Too large cutoff radius
    write( IOBuffer, '("Cutoff radius is", I10, " times (", F6.2, "%) too large")' ) &
&          this%NRCutoffMax, ( 100._RK * this%NRCutoffMax ) / Step
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
      Average = this%SumTemperature%Average
      write( IOBuffer, '("Simulation temperature", T29, "reduced:", F20.9)' ) Average
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T32, "in K:", F20.9)' ) Average * UnitTemperature
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Mole fractions of liquid phase
      do i = 1, this%NComponents
        pc => this%Component(i)
        write( IOBuffer, '("Liquid mole fraction of ", A, T36, ":", F20.9)' ) &
&              trim( pc%Molecule%PotModFileName ), pc%LiqFraction
        call FileWrite( this%iounit_errors )
      end do
      call FileWriteBlank( this%iounit_errors )

      if ( this%OptPressure ) then
         Average = this%SumPressure%Average
         Variance = this%SumPressure%Variance
      else
         Average = this%RefPressure
         Variance = 0._RK
      end if
      ! Simulation pressure of liquid phase
      write( IOBuffer, '("Liquid simulation pressure", T29, "reduced:", F20.9)' ) Average
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T30, "in MPa:", F20.9)' ) Average * UnitPressure * 1e-6_RK
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Vapor pressure
    if ( this%OptPressure ) then
      AvgPressure = this%SumPressure%Average
      NN = 0._RK
      do i = 1, this%NComponents
        NN = NN + this%Component(i)%SumFraction%Average * this%Component(i)%PartialMolarVolume
      end do
      NN = NN - this%Temperature / Average
      do i = 1, this%NComponents
        pc => this%Component(i)
        dpdmu(i) = -this%Temperature * pc%SumFraction%Average / NN
        dpdv(i) = -pc%SumFraction%Average * (Average-this%RefPressure) / NN
        varmu(i) = pc%VarChemPot
        varv(i) = pc%VarPartialMolarVolume
      end do

      VarPressure = sqrt( this%SumPressure%Variance**2 + sum( (dpdmu * varmu)**2 ) + sum( (dpdv * varv)**2 ) )
    else
      AvgPressure = this%RefPressure
      NN = 0._RK
      do i = 1, this%NComponents
        NN = NN + this%Component(i)%SumFraction%Average * this%Component(i)%PartialMolarVolume
      end do
      NN = NN - this%Temperature / Average
      do i = 1, this%NComponents
        pc => this%Component(i)
        dpdmu(i) = -this%Temperature * pc%SumFraction%Average / NN
        dpdv(i) = 0._RK
        varmu(i) = pc%VarChemPot
        varv(i) = pc%VarPartialMolarVolume
      end do
      VarPressure = sqrt( sum( (dpdmu * varmu)**2 ) )
    end if

      write( IOBuffer, '("Vapor pressure", T29, "reduced:", 2F20.9)' ) AvgPressure, VarPressure
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T30, "in MPa:", 2F20.9)' ) AvgPressure * UnitPressure * 1E-6_RK, VarPressure * UnitPressure * 1E-6_RK
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Mole fractions of vapor phase
      do i = 1, this%NComponents
        pc => this%Component(i)
        yvi = pc%SumFraction%Average * ( pc%PartialMolarVolume / this%Temperature - 1 / Average )
        do j = 1, this%NComponents
          dydmu(i, j) = yvi * dpdmu(j)
          dydv(i, j) = yvi * dpdv(j)
        end do
        dydmu(i, i) = dydmu(i, i) + pc%SumFraction%Average
        dydv(i, i) = dydv(i, i) + pc%SumFraction%Average * 1 / this%Temperature * ( Average - this%RefPressure )
      end do

      do i = 1, (this%NComponents - 1)
        pc => this%Component(i)
        Average = pc%SumFraction%Average
        vary(i) = sqrt( pc%SumFraction%Variance**2 + sum( (dydmu(i, :) * varmu)**2 ) + sum( (dydv(i, :) * varv)**2 ) )
        write( IOBuffer, '("Vapor mole fraction of ", A, T36, ":", 2F20.9)' ) &
&              trim( pc%Molecule%PotModFileName ), Average, vary(i)
        call FileWrite( this%iounit_errors )
      end do

      pc => this%Component( this%NComponents )
      Average = pc%SumFraction%Average
      Variance = sqrt( sum( vary(1:(this%NComponents - 1))**2 ) )
      write( IOBuffer, '("Vapor mole fraction of ", A, T36, ":", 2F20.9)' ) &
&            trim( pc%Molecule%PotModFileName ), Average, Variance
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Saturated liquid density
      Average = this%LiqDensity + this%LiqDensity * this%LiqBetaT * ( AvgPressure - this%RefPressure)

      Variance = sqrt( this%VarLiqDensity**2 + ( this%VarLiqBetaT * ( AvgPressure - this%RefPressure )&
&                + VarPressure * this%LiqBetaT )**2 )

      write( IOBuffer, '("Liquid density", T29, "reduced:", 2F20.9)' ) Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in mol/l:", 2F20.9)' ) Average * UnitDensity, Variance * UnitDensity
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Saturated vapor density
      Average = this%SumDensity%Average
      Variance = this%SumDensity%Variance
      write( IOBuffer, '("Vapor density", T29, "reduced:", 2F20.9)' ) &
&            Average, Average * VarPressure / AvgPressure

      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in mol/l:", 2F20.9)' ) Average * UnitDensity, Average&
&            * VarPressure / AvgPressure * UnitDensity
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Saturated liquid enthalpy
      Average = this%LiqEnthalpy + this%LiqdHdP * ( AvgPressure - this%RefPressure )

      Variance = sqrt( this%VarLiqEnthalpy**2 + ( this%VarLiqdHdP * &
&                ( AvgPressure - this%RefPressure ) + VarPressure * this%LiqdHdP )**2 )

      write( IOBuffer, '("Liquid enthalpy", T29, "reduced:", 2F20.9)' ) Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) Average * UnitEnergy * NAvogadro, &
&            Variance * UnitEnergy * NAvogadro
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      DeltaHv = Average
      VarDeltaHv = Variance

      ! Saturated vapor enthalpy
      Average = this%SumEnthalpy%Average
      Variance = this%SumEnthalpy%Variance
      write( IOBuffer, '("Vapor enthalpy", T29, "reduced:", 2F20.9)' ) Average, Variance
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) Average * UnitEnergy * NAvogadro, &
&            Variance * UnitEnergy * NAvogadro
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      DeltaHv = Average - DeltaHv
      VarDeltaHv = Variance + VarDeltaHv

      ! Evaporation enthalpy
      write( IOBuffer, '("Enthalpy of vaporization", T29, "reduced:", 2F20.9)' ) DeltaHv, VarDeltaHv
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '(T28, "in J/mol:", 2F20.9)' ) DeltaHv * UnitEnergy * NAvogadro, &
&            VarDeltaHv * UnitEnergy * NAvogadro
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )

      ! Separator
      write( IOBuffer, '(76("="))' )
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )
    end if

    if( (SimulationType .eq. MolecularDynamics) .and. (EnsembleType .eq. EnsembleTypeGE) ) then
      ! Statistics section
      write( IOBuffer, '("Statistics")' )
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '("----------")' )
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )
    end if

    if( (SimulationType .eq. MonteCarlo) .or. (SimulationType .eq. Gibbs) ) then
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
&            NRootProc, Communicator, ierror )
        call MPI_Reduce( this%NResizeAttempts,tempVal2, 1, MPI_INTEGER, MPI_SUM, &
&            NRootProc, Communicator, ierror )
        if ( Nproc == NRootProc) then
          write( IOBuffer, '("Acceptance rate volume changes", T32, "in %:", F20.9)' ) &
&                100._RK * real(tempVal, RK ) / real (tempVal2, RK )
        endif

#else
        write( IOBuffer, '("Acceptance rate volume changes", T32, "in %:", F20.9)' ) &
&              100._RK * real( this%NResizeSuccesses, RK ) / real ( this%NResizeAttempts, RK )
#endif         
 
        call FileWrite( this%iounit_errors )
#if MPI_VER > 0
        call MPI_Reduce( this%DispVol,tempReal, 1, MPI_RK, MPI_MAX, &
        &    NRootProc, Communicator, ierror )
        if (Nproc == NRootProc) then
          write( IOBuffer, '("Maximum displacement volume", T33, "r`d:", F20.9)' ) tempReal
        endif

#else
          write( IOBuffer, '("Maximum displacement volume", T33, "r`d:", F20.9)' ) this%DispVol  
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
        &    NRootProc, Communicator, ierror )
        call MPI_Reduce( pc%NMoveAttempts,tempVal2, 1, MPI_INTEGER, MPI_SUM, &
        &    NRootProc, Communicator, ierror )

        if (Nproc == NRootProc) then
          write( IOBuffer, '("Acceptance rate trans.", T32, "in %:", F20.9)' ) &
&                100._RK * real( tempVal, RK ) / real ( tempVal2, RK ) 
        endif

#else
          write( IOBuffer, '("Acceptance rate trans.", T32, "in %:", F20.9)' ) &
&                100._RK * real( pc%NMoveSuccesses, RK ) / real ( pc%NMoveAttempts, RK )   
#endif

        call FileWrite( this%iounit_errors )

        if( this%NDFRot > 0 ) then
#if MPI_VER > 0
          call MPI_Reduce( pc%NRotateSuccesses,tempVal, 1, MPI_INTEGER, MPI_SUM, &
&              NRootProc, Communicator, ierror )
          call MPI_Reduce( pc%NRotateAttempts,tempVal2, 1, MPI_INTEGER, MPI_SUM, &
&              NRootProc, Communicator, ierror )

          if (Nproc == NRootProc) then
            write( IOBuffer, '(T17, "rotates", T32, "in %:", F20.9)' ) 100._RK &
&                  * real( tempVal, RK ) / real (tempVal2, RK )
          endif

#else
            write( IOBuffer, '(T17, "rotates", T32, "in %:", F20.9)' ) 100._RK &
&                  * real( pc%NRotateSuccesses, RK ) / real ( pc%NRotateAttempts, RK )
#endif
          call FileWrite( this%iounit_errors )
        end if

        if( pc%ChemPotMethod .eq. ChemPotMethodGradIns ) then
          ! Biased move and rotate acceptance rates
#if MPI_VER > 0
          call MPI_Reduce( pc%NMoveBiasedSuccesses,tempVal, 1, MPI_INTEGER, MPI_SUM, &
          &    NRootProc, Communicator, ierror )
          call MPI_Reduce( pc%NMoveBiasedAttempts,tempVal2, 1, MPI_INTEGER, MPI_SUM, &
          &    NRootProc, Communicator, ierror )

          if (Nproc == NRootProc) then
            write( IOBuffer, '(T17, "biased trans.", T32, "in %:", F20.9)' ) 100._RK * real(tempVal, RK ) / &
&                  real ( tempVal2, RK )
          endif

#else
          write( IOBuffer, '(T17, "biased trans.", T32, "in %:", F20.9)' ) &
&                100._RK * real( pc%NMoveBiasedSuccesses, RK ) / real ( pc%NMoveBiasedAttempts, RK )
#endif

          call FileWrite( this%iounit_errors )
          if( this%NDFRot > 0 ) then
#if MPI_VER > 0
            call MPI_Reduce( pc%NRotateBiasedSuccesses,tempVal, 1, MPI_INTEGER, MPI_SUM, &
&                NRootProc, Communicator, ierror )
            call MPI_Reduce( pc%NRotateBiasedAttempts,tempVal2, 1, MPI_INTEGER, MPI_SUM, &
&                NRootProc, Communicator, ierror )

            if (Nproc == NRootProc) then
              write( IOBuffer, '(T17, "biased rotates", T32, "in %:", F20.9)' ) &
&                    100._RK * real(tempVal, RK ) / real (tempVal2, RK )
            endif

#else
            write( IOBuffer, '(T17, "biased rotates", T32, "in %:", F20.9)' ) &
&                  100._RK * real( pc%NRotateBiasedSuccesses, RK ) / real ( pc%NRotateBiasedAttempts, RK )
#endif

            call FileWrite( this%iounit_errors )
          end if
        end if

        ! Maximum translational and rotational displacements
#if MPI_VER > 0
        call MPI_Reduce( pc%DispTran,tempReal, 1, MPI_RK, MPI_MAX, &
&            NRootProc, Communicator, ierror )
        if (Nproc == NRootProc) then
          write( IOBuffer, '("Maximum displacement trans.", T33, "r`d:", F20.9)' ) tempReal
        endif

#else
          write( IOBuffer, '("Maximum displacement trans.", T33, "r`d:", F20.9)' ) pc%DispTran
#endif  

        call FileWrite( this%iounit_errors )
        if( this%NDFRot > 0 ) then
#if MPI_VER > 0
        call MPI_Reduce( pc%DispRot,tempReal, 1, MPI_RK, MPI_MAX, &
&            NRootProc, Communicator, ierror )
        if (Nproc == NRootProc) then
          write( IOBuffer, '(T22, "rotational", T33, "r`d:", F20.9)' ) tempReal
        endif

#else
          write( IOBuffer, '(T22, "rotational", T33, "r`d:", F20.9)' ) pc%DispRot
#endif  
          call FileWrite( this%iounit_errors )
        end if

        ! Maximum molecular translational and rotational displacements
        if ( UseIntDegFreed ) then ! Michael Sch.: consider other processes in MC!!!
          write( IOBuffer, '("Acceptance rate mol. moves", T32, "in %:", F20.9)' ) &
&            100._RK * real( pc%NMoveMolSuccesses, RK ) / real ( pc%NMoveMolAttempts, RK )
          call FileWrite( this%iounit_errors )
          if( pc%Molecule%IsElongated ) then
            write( IOBuffer, '(T17, "mol. rotates", T32, "in %:", F20.9)' ) 100._RK &
&              * real( pc%NRotateMolSuccesses, RK ) / real ( pc%NRotateMolAttempts, RK )
            call FileWrite( this%iounit_errors )
          end if
          write( IOBuffer, '("Maximum displ. mol. trans.", T33, "r`d:", F20.9)' ) pc%DispMolTran
          call FileWrite( this%iounit_errors )
          if( pc%Molecule%IsElongated ) then
            write( IOBuffer, '(T17, "mol. rotational", T33, "r`d:", F20.9)' ) pc%DispMolRot
            call FileWrite( this%iounit_errors )
          end if
          call FileWriteBlank( this%iounit_errors )
        end if
        call FileWriteBlank( this%iounit_errors )

        ! Gradual insertion change fluctuating particle acceptance rates
        if( pc%ChemPotMethod .eq. ChemPotMethodGradIns ) then
          write(IOBuffer, '("Acceptance rate gradual insertion change fluctuating particle moves:")')
          call FileWrite( this%iounit_errors )
          write(IOBuffer, '("  up        down (%)")')
          call FileWrite( this%iounit_errors )
          write(IOBuffer, '("  --------  --------")')
          call FileWrite( this%iounit_errors )

#if MPI_VER > 0
          call MPI_Reduce( pc%NFluctUpSuccesses(:),tempVec1(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, &
&              MPI_SUM, NRootProc, Communicator, ierror )
          call MPI_Reduce( pc%NFluctUpAttempts(:),tempVec2(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, &
&              MPI_SUM, NRootProc, Communicator, ierror )
          call MPI_Reduce( pc%NFluctDownSuccesses(:),tempVec3(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, &
&              MPI_SUM, NRootProc, Communicator, ierror )
          call MPI_Reduce( pc%NFluctDownAttempts(:),tempVec4(1:pc%NFluctMax), pc%NFluctMax, MPI_INTEGER, &
&              MPI_SUM, NRootProc, Communicator, ierror )

          if (Nproc == NRootProc) then
            write(IOBuffer, '(2F10.4)') 0._RK, real(tempVec3(pc%NFluctMax), RK) / &
&                 real(tempVec4(pc%NFluctMax), RK) * 100._RK
            call FileWrite( this%iounit_errors )           
           
           do j = pc%NFluctMax -1, 1, -1
             write(IOBuffer, '(2F10.4)') real(tempVec1(j+1), RK) / &
&                  real(tempVec2(j+1), RK) * 100._RK, real(tempVec3(j), RK) / &
&                  real(tempVec4(j), RK) * 100._RK
             call FileWrite( this%iounit_errors )
           end do

           write(IOBuffer, '(2F10.4)') real(tempVec1(1), RK) / real(tempVec2(1), RK) * 100._RK, 0._RK
           call FileWrite( this%iounit_errors )
           call FileWriteBlank( this%iounit_errors )
           call FileWriteBlank( this%iounit_errors )
          endif

#else
          write(IOBuffer, '(2F10.4)') 0._RK, real(pc%NFluctDownSuccesses(pc%NFluctMax), RK) / &
&               real(pc%NFluctDownAttempts(pc%NFluctMax), RK) * 100._RK
          call FileWrite( this%iounit_errors )
          do j = pc%NFluctMax -1, 1, -1
            write(IOBuffer, '(2F10.4)') real(pc%NFluctUpSuccesses(j+1), RK) / &
&                 real(pc%NFluctUpAttempts(j+1), RK) * 100._RK, real(pc%NFluctDownSuccesses(j), RK) / &
&                 real(pc%NFluctDownAttempts(j), RK) * 100._RK
            call FileWrite( this%iounit_errors )
          end do

          write(IOBuffer, '(2F10.4)') real(pc%NFluctUpSuccesses(1), RK) / &
&               real(pc%NFluctUpAttempts(1), RK) * 100._RK, 0._RK
          call FileWrite( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )
          call FileWriteBlank( this%iounit_errors )

#endif
        end if
      end do

      ! Inserts and deletes acceptance rates
      if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA ) then
#if MPI_VER > 0
        call MPI_Reduce( this%NInsertSuccesses,tempVal, 1, MPI_INTEGER, MPI_SUM, &
&            NRootProc, Communicator, ierror )
        call MPI_Reduce( this%NInsertAttempts,tempVal2, 1, MPI_INTEGER, MPI_SUM, &
&            NRootProc, Communicator, ierror )

        if (Nproc == NRootProc) then
          write( IOBuffer, '("Acceptance rate inserts", T32, "in %:", F20.9)' ) &
&                100._RK * real( tempVal, RK ) / real ( tempVal2, RK )
        endif

#else
        write( IOBuffer, '("Acceptance rate inserts", T32, "in %:", F20.9)' ) &
&              100._RK * real( this%NInsertSuccesses, RK ) / real ( this%NInsertAttempts, RK )
#endif 

        call FileWrite( this%iounit_errors )
#if MPI_VER > 0
        call MPI_Reduce( this%NDeleteSuccesses,tempVal, 1, MPI_INTEGER, MPI_SUM, &
&            NRootProc, Communicator, ierror )
        call MPI_Reduce( this%NDeleteAttempts,tempVal2, 1, MPI_INTEGER, MPI_SUM, &
&            NRootProc, Communicator, ierror )

        if (Nproc == NRootProc) then
          write( IOBuffer, '("Acceptance rate deletes", T32, "in %:", F20.9)' ) &
&                100._RK * real(tempVal, RK ) / real ( tempVal2, RK )
        endif

#else
          write( IOBuffer, '("Acceptance rate deletes", T32, "in %:", F20.9)' ) &
&                100._RK * real( this%NDeleteSuccesses, RK ) / real ( this%NDeleteAttempts, RK )
#endif  
       
        call FileWrite( this%iounit_errors )
        call FileWriteBlank( this%iounit_errors )
      end if
    elseif ( SimulationType .eq. MolecularDynamics .and. EnsembleType .eq. EnsembleTypeGE ) then
      write( IOBuffer, '("Acceptance rate inserts", T32, "in %:", F20.9)' ) &
&            100._RK * real( this%NInsertSuccesses, RK ) / real ( this%NInsertAttempts, RK )
      call FileWrite( this%iounit_errors )
      write( IOBuffer, '("Acceptance rate deletes", T32, "in %:", F20.9)' ) &
&            100._RK * real( this%NDeleteSuccesses, RK ) / real ( this%NDeleteAttempts, RK )
      call FileWrite( this%iounit_errors )
      call FileWriteBlank( this%iounit_errors )
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


! Calculation of residence times
    if ( this%ResidenceTime ) then
      do i = 1, NBlockSizes
        do j = i, NBlocks, i
          if ( sum(this%SumResidenceDuration%NBlockSum(j - i + 1:j)) == 0 ) then
             this%SumResidenceDuration%NBlockSum(j - i + 1) = 1
             this%SumResidenceDuration%BlockSum(j - i + 1) = this%SumResidenceDuration%Average
          end if
        end do
      end do
      write(IOBuffer, '("Average pairs between")' )
      call FileWrite( this%iounit_errors )
      write(IOBuffer, '("Comp.",I2," Site",I2,"  and Comp.",I2," Site",I2," =", F14.5)' ) &
&           this%ResidComp1, this%ResidSite1, &
&           this%ResidComp2, this%ResidSite2, this%SumResidencePairs%Average/this%Component(this%ResidComp1)%NPart
      call FileWrite( this%iounit_errors )
      write(IOBuffer, '("Average residence time between")' )
      call FileWrite( this%iounit_errors )
      call Error (this%SumResidenceDuration)

      if ( (this%SumResidenceDuration%NTotalsum .eq. 0) .and. (this%ResidPairs .ne. 0) ) then
         write(IOBuffer, '("Comp.",I2," Site",I2,"  and Comp.",I2," Site",I2," =" F20.5" fs")' ) &
&           this%ResidComp1, this%ResidSite1, &
&           this%ResidComp2, this%ResidSite2, Step*TimeStep* UnitTime * 1E15_RK
        call FileWrite( this%iounit_errors )
        write(IOBuffer, '("No separation between the two components observed")' )

      else if ( (this%SumResidenceDuration%NTotalsum .eq. 0) .and. (this%ResidPairs .eq. 0) ) then
        write(IOBuffer, '("Comp.",I2," Site",I2,"  and Comp.",I2," Site",I2," =" F14.5" fs")' ) &
&           this%ResidComp1, this%ResidSite1, this%ResidComp2,this%ResidSite2,&
&           this%ResidenceDuration*UnitTime*1E15_RK
        call FileWrite( this%iounit_errors )
        write(IOBuffer, '("No pairing between the two components observed")' )

      else
        write(IOBuffer, '("Comp.",I2," Site",I2,"  and Comp.",I2," Site",I2," =" F14.5" fs +-",F10.5)' ) &
&         this%ResidComp1,this%ResidSite1, &
&         this%ResidComp2,this%ResidSite2, this%SumResidenceDuration%Average*UnitTime*1E15_RK ,&
&         this%SumResidenceDuration%Variance*UnitTime*1E15_RK
      end if

      call FileWrite( this%iounit_errors )
      write(IOBuffer, '("Critical distance: ",F10.5," A")' ) &
&           this%ResidLength*UnitLength/Angstroem
      call FileWrite( this%iounit_errors )
    end if
    call FileWriteBlank( this%iounit_errors )

    ! Close final result file
    call FileClose( this%iounit_errors )

    ! Open ThermoInt result file
    if ( any(this%Component(:)%ChemPotMethod .eq. ChemPotMethodThermoInt)) then
      if (RootProc) then
        write( IOBuffer, '(I16)' ) this%EnsembleNumber
        call FileRewrite( this%iounit_thermoint, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ThermoIntFileExtension)
      end if

      t = this%NRealComponents+1
      do i=1,this%NRealComponents
         pc => this%Component(i)
        if (pc%ChemPotMethod .eq. ChemPotMethodThermoInt) then

          !first two lines
          if (RootProc) then
            write( IOBuffer, '(" Component:", T15, I3)' ) i
            call FileWrite( this%iounit_thermoint )
            write( IOBuffer, '("currentlambda =", T20, F8.5)' ) this%Component(t)%lambda
            call FileWrite( this%iounit_thermoint )
            write( IOBuffer, '(" BINID")' )
            call FileWriteNoAdvance( this%iounit_thermoint )
            write( IOBuffer, '(" LAMBDA")' )
            call FileWriteNoAdvance( this%iounit_thermoint )
            write( IOBuffer, '("            EPOT")' )
            call FileWriteNoAdvance( this%iounit_thermoint )
            write( IOBuffer, '("    dEPOTdLAMBDA")' )
            call FileWriteNoAdvance( this%iounit_thermoint )
            write( IOBuffer, '("   dEPOTdLAMBDAV")' )
            call FileWriteNoAdvance( this%iounit_thermoint )
            write( IOBuffer, '("   dEPOTdLAMBDAH")' )
            call FileWriteNoAdvance( this%iounit_thermoint )
            write( IOBuffer, '(" INTdEPOTdLAMBDA")' )
            call FileWriteNoAdvance( this%iounit_thermoint )
            write( IOBuffer, '("       INTParVol")' )
            call FileWriteNoAdvance( this%iounit_thermoint )
            write( IOBuffer, '("       INTParEnt")' )
            call FileWriteNoAdvance( this%iounit_thermoint )
            write( IOBuffer, '("     VISITS")' )
            call FileWriteNoAdvance( this%iounit_thermoint )
            call FileWriteBlank( this%iounit_thermoint )
          end if

          call ErrorsUpdateThermoInt( this, i, pc%NBins )
          t = t+1
        end if

      end do

      ! Close final result file
      call FileClose( this%iounit_thermoint)
    end if

  end subroutine TEnsemble_ErrorsUpdate

!==============================================================!
!  Subroutine TEnsemble_ErrorsUpdateThermoInt                  !
!==============================================================!

  subroutine TEnsemble_ErrorsUpdateThermoInt( this, i, NBins )

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif
    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: i
    integer, intent(in) :: NBins

    ! Declare local variables
    type(TComponent), pointer  :: pc
    integer                    :: j, LocalVisit
    real(RK)                   :: BinsEn(0:NBins-1), BinsdEndLa(0:NBins-1), BinsIntdEndLa(0:NBins-1)
    real(RK)                   :: BinsdEndLaV(0:NBins-1), BinsdEndLaH(0:NBins-1), BinsIntVW(0:NBins-1), BinsIntHW(0:NBins-1)
    integer                    :: BinsVisit(0:NBins-1)

    pc => this%Component(i)
    ! Avearge of each MPI process's histogram is saved in the thi file 
#if MPI_VER > 0
    if (SimulationType .eq. MonteCarlo) then
      call MPI_Reduce( pc%BinsEn(0: NBins-1)       *pc%BinsVisit(0: NBins-1), BinsEn(0: NBins-1), NBins, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( pc%BinsdEndLa(0: NBins-1)   *pc%BinsVisit(0: NBins-1), BinsdEndLa(0: NBins-1), NBins, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( pc%BinsdEndLaV(0: NBins-1)  *pc%BinsVisit(0: NBins-1), BinsdEndLaV(0: NBins-1), NBins, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( pc%BinsdEndLaH(0: NBins-1)  *pc%BinsVisit(0: NBins-1), BinsdEndLaH(0: NBins-1), NBins, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( pc%BinsIntdEndLa(0: NBins-1)                         , BinsIntdEndLa(0: NBins-1), NBins, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( pc%BinsIntVW(0: NBins-1)                             , BinsIntVW(0: NBins-1), NBins, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( pc%BinsIntHW(0: NBins-1)                             , BinsIntHW(0: NBins-1), NBins, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( pc%BinsVisit(0: NBins-1)                             , BinsVisit(0: NBins-1), NBins, MPI_INTEGER, MPI_SUM, NRootProc, Communicator, ierror )
      do j=0,pc%NBins-1
        LocalVisit=BinsVisit(j)
        if (LocalVisit == 0) LocalVisit=1 
        BinsEn(j)        = BinsEn(j)/LocalVisit
        BinsdEndLa(j)    = BinsdEndLa(j)/LocalVisit
        BinsdEndLaV(j)   = BinsdEndLaV(j)/LocalVisit
        BinsdEndLaH(j)   = BinsdEndLaH(j)/LocalVisit
        BinsIntdEndLa(j) = BinsIntdEndLa(j)/NProcs
        BinsIntVW(j)     = BinsIntVW(j)/NProcs
        BinsIntHW(j)     = BinsIntHW(j)/NProcs
      end do
    else
      BinsEn(:)        = pc%BinsEn(:)
      BinsdEndLa(:)    = pc%BinsdEndLa(:)
      BinsdEndLaV(:)   = pc%BinsdEndLaV(:)
      BinsdEndLaH(:)   = pc%BinsdEndLaH(:)
      BinsIntdEndLa(:) = pc%BinsIntdEndLa(:)
      BinsIntVW(:)     = pc%BinsIntVW(:)
      BinsIntHW(:)     = pc%BinsIntHW(:)
      BinsVisit(:)     = pc%BinsVisit(:)
    endif
#else
      BinsEn(:)        = pc%BinsEn(:)
      BinsdEndLa(:)    = pc%BinsdEndLa(:)
      BinsdEndLaV(:)   = pc%BinsdEndLaV(:)
      BinsdEndLaH(:)   = pc%BinsdEndLaH(:)
      BinsIntdEndLa(:) = pc%BinsIntdEndLa(:)
      BinsIntVW(:)     = pc%BinsIntVW(:)
      BinsIntHW(:)     = pc%BinsIntHW(:)
      BinsVisit(:)     = pc%BinsVisit(:)
#endif

    ! Rest.
    if (RootProc) then
      do j=0,pc%NBins-1
        write( IOBuffer, '(I6)' ) j
        call FileWriteNoAdvance( this%iounit_thermoint )
        write( IOBuffer, '(" ",F6.3)' ) pc%LaMin+j*pc%deltaLa
        call FileWriteNoAdvance( this%iounit_thermoint )
        write( IOBuffer, '(" ",E15.6)' ) BinsEn(j)
        call FileWriteNoAdvance( this%iounit_thermoint )
        write( IOBuffer, '(" ",E15.6)' ) BinsdEndLa(j)
        call FileWriteNoAdvance( this%iounit_thermoint )
        write( IOBuffer, '(" ",E15.6)' ) BinsdEndLaV(j)
        call FileWriteNoAdvance( this%iounit_thermoint )
        write( IOBuffer, '(" ",E15.6)' ) BinsdEndLaH(j)
        call FileWriteNoAdvance( this%iounit_thermoint )
        write( IOBuffer, '(" ",E15.6)' ) BinsIntdEndLa(j)
        call FileWriteNoAdvance( this%iounit_thermoint )
        write( IOBuffer, '(" ",E15.6)' ) BinsIntVW(j)
        call FileWriteNoAdvance( this%iounit_thermoint )
        write( IOBuffer, '(" ",E15.6)' ) BinsIntHW(j)
        call FileWriteNoAdvance( this%iounit_thermoint )
        write( IOBuffer, '(" ",I10)' ) BinsVisit(j)
        call FileWriteNoAdvance( this%iounit_thermoint )
        call FileWriteBlank( this%iounit_thermoint )
      end do
    end if

  end subroutine TEnsemble_ErrorsUpdateThermoInt


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
    call FileRewrite( this%iounit_errors, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//ErrorsFileExtension )

    write( IOBuffer, '(76("="))')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("*                           Publishing with ms2                                *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* Every user agrees to cite ms2 upon usage as follows                          *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* ---------------------------------------------------------------------------- *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* C.W. Glass, S. Reiser, G. Rutkai, S. Deublein, A. Koster, G. Guevara-Carrion *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* A. Wafai, M. Horsch, M. Bernreuther, T. Windmann, H. Hasse, J. Vrabec        *")')
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("* Computer Physics Communications (2014)                                       *")')
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
    write( IOBuffer, '("Simulation type", T36, ":", 9X, A)' ) trim( SimulationTypeString )
    call FileWrite( this%iounit_errors )

    ! Number of orientations
    write( IOBuffer, '("Number of orientations", T36, ":", I10)' ) NOrient
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Number of radial steps", T36, ":", I10)' ) NSteps
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Minimum radius", T29, "reduced:", F20.9)' ) MinRadius
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T32, "in A:", F20.9)' ) MinRadius * UnitLength / Angstroem
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Maximum radius", T29, "reduced:", F20.9)' ) MaxRadius
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T32, "in A:", F20.9)' ) MaxRadius * UnitLength / Angstroem
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! Temperature
    write( IOBuffer, '("Temperature", T29, "reduced:", F20.9)' ) this%Temperature
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '(T32, "in K:", F20.9)' ) this%Temperature * UnitTemperature
    call FileWrite( this%iounit_errors )
    call FileWriteBlank( this%iounit_errors )

    ! System of units
    write( IOBuffer, '("Unit of length", T36, ":", F20.9, " A")' ) UnitLength / Angstroem
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Unit of energy", T36, ":", F20.9, " K")' ) UnitEnergy / kBoltzmann
    call FileWrite( this%iounit_errors )
    write( IOBuffer, '("Unit of mass", T36, ":", F20.9, " a.u.")' ) UnitMass * NAvogadro * 1000._RK
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
&               .5_RK * this%Interaction(i, j)%EPotCorrLJ / this%Temperature
        write( IOBuffer, '("2. VC of ", A, "-", A, T29, "reduced:", F20.9)' ) &
&              trim( this%Component(i)%Molecule%PotModFileName ), &
&              trim( this%Component(j)%Molecule%PotModFileName ), value
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T28, "in l/mol:", F20.9)' ) value / UnitDensity
        call FileWrite( this%iounit_errors )
      end do
    end do
    call FileWriteBlank( this%iounit_errors )

    ! Temperature deviation of second virial coefficient
    do i = 1, this%NComponents, 2
      do j = i + 1, this%NComponents, 2
        value = ( this%Interaction(i, j)%IntFFunction2(NSteps) - this%Interaction(i,j)%IntFFunction1(NSteps) ) &
&               / ( .0002_RK * this%Temperature )
        write( IOBuffer, '("dB/dT of ", A, "-", A, T29, "reduced:", F20.9)' ) &
&              trim( this%Component(i)%Molecule%PotModFileName ), &
&              trim( this%Component(j)%Molecule%PotModFileName ), value
        call FileWrite( this%iounit_errors )
        write( IOBuffer, '(T24, "in l/(mol K):", F20.9)' ) value / ( UnitDensity * UnitTemperature )
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
    integer                   :: i, j, k, num
    type(TSiteLJ126), pointer :: psLJ126
    type(TSiteCharge), pointer :: psCharge
    real(RK)                   :: ch_sig

    ! Open visualization file
    write( IOBuffer, '(I16)' ) this%EnsembleNumber
    call FileRewrite( this%iounit_visual, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//VisualFileExtension )

    ! Create header
    num = 0
    do i = 1, this%NComponents
      do k = 1, this%Component(i)%Molecule%NUnit
        if (this%Component(i)%Molecule%Unit(k)%NLJ126 > 0) then
          do j = 1, this%Component(i)%Molecule%Unit(k)%NLJ126
            psLJ126 => this%Component(i)%Molecule%Unit(k)%SiteLJ126(j)
            write( IOBuffer, '("~", I3, "     LJ", 4F8.4, "  1")' ) (num+k), psLJ126%r(:) * UnitLength / Angstroem, &
&                  psLJ126%sig  * UnitLength / Angstroem
            call FileWrite( this%iounit_visual )
          end do
        else  ! For visualisation of Units with no LJ sites
          ch_sig = UnitLength * 0.2
          do j = 1, this%Component(i)%Molecule%Unit(k)%NCharge
            psCharge => this%Component(i)%Molecule%Unit(k)%SiteCharge(j)
            write( IOBuffer, '("~", I3, " Charge", 4F8.4, "  1")' ) (num+k), &
&              psCharge%r(:) * UnitLength / Angstroem, ch_sig
            call FileWrite( this%iounit_visual )
          end do
        end if
      end do
      num = num+this%Component(i)%Molecule%NUnit
    end do
    call FileWriteBlank( this%iounit_visual )

#if HBOND > 0
    !Open visualization file for H-bondings
    write( IOBuffer, '(I16)' ) this%EnsembleNumber
    call FileRewrite( this%iounit_visualHB, &
&     trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//VisualHBFileExtension )
    write( IOBuffer, '("!"," Nr", " MH", " MO_1", " MO_2")' ) 
    call FileWrite( this%iounit_visualHB )
    call FileWriteBlank( this%iounit_visualHB )
#endif

  end subroutine TEnsemble_VisualOpen


!==============================================================!
!  Subroutine TEnsemble_VisualUpdate                           !
!==============================================================!

  subroutine TEnsemble_VisualUpdate( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer  :: i, j, k, num
    logical  :: l
    real(RK) :: r(3), q(4)

    ! Update visualization file
    num = 0
    write( IOBuffer, '("#", F10.4, "  new Frame")' ) this%BoxLength * UnitLength / Angstroem
    call FileWrite( this%iounit_visual )
    do i = 1, this%NComponents
      do j = 1, this%Component(i)%NPart
        do k = 1, this%Component(i)%Molecule%NUnit
          l = this%Component(i)%Molecule%Unit(k)%isElongated
          r(:) = this%Component(i)%P0(j, :, k) + .5_RK

          if( l ) then
            q(:) = this%Component(i)%Q0(j, :, k)
          else
            q(1) = 1._RK
            q(2:4) = .0_RK
          end if

          write( IOBuffer, '("!", I3,  3I5, 4I5)' ) (num+k),  nint( r(:) * 999.99_RK ), nint( q(:) * 999.99_RK )
          call FileWrite( this%iounit_visual )
        end do
      end do
      num = num+(i)*this%Component(i)%Molecule%NUnit
    end do
    call FileWriteBlank( this%iounit_visual )

  end subroutine TEnsemble_VisualUpdate


#if HBOND > 0
!==============================================================!
!  Subroutine TEnsemble_VisualUpdateHB                         !
!==============================================================!

  subroutine TEnsemble_VisualUpdateHB( this, np, MH, MO )

    implicit none

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: np
    integer, intent(in) :: MH(np)
    integer, intent(in) :: MO(2,np)

    ! Declare local variables
    integer  :: i

    ! Update visualization file
    write( IOBuffer, '("#", F10.4, "  new Frame")' ) this%BoxLength * UnitLength / Angstroem
    call FileWrite( this%iounit_visualHB )                                                       
    do i= 1, np
      write( IOBuffer, '("!", I5, I5, I5, I5)' ) i, MH(i), MO(1,i), MO(2,i) 
      call FileWrite( this%iounit_visualHB ) 
    end do 
    call FileWriteBlank( this%iounit_visualHB )

  end subroutine TEnsemble_VisualUpdateHB
#endif


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

#if HBOND > 0
    ! Close visualization H-bonding file
    write( IOBuffer, '("##")' )
    call FileWrite( this%iounit_visualHB )
    call FileClose( this%iounit_visualHB )
#endif

  end subroutine TEnsemble_VisualClose


#if OSMOP > 0
!==============================================================!
!  Subroutine TEnsemble_ProfileOpen                            !
!==============================================================!

  subroutine TEnsemble_ProfileOpen( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Open profile file
    write( IOBuffer, '(I16)' ) this%EnsembleNumber
    call FileRewrite( this%iounit_dcp, &
&     trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//DCPFileExtension )

    call FileWriteBlank( this%iounit_dcp )
    call FileClose( this%iounit_dcp )

  end subroutine TEnsemble_ProfileOpen


!==============================================================!
!  Subroutine TEnsemble_ProfileUpdate                          !
!==============================================================!

  subroutine TEnsemble_ProfileUpdate( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer  :: i, j
    real(RK) :: Variance, Average
    type(TComponent), pointer :: pc

    ! Open profile file
    write( IOBuffer, '(I16)' ) this%EnsembleNumber
    call FileRewrite( this%iounit_dcp, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//DCPFileExtension )

    ! Create header
    write( IOBuffer, '("#", F10.4, "  BoxLength / A")' ) this%BoxLength * UnitLength / Angstroem
    call FileWrite( this%iounit_dcp )
    call FileWriteBlank( this%iounit_dcp )
    write( IOBuffer, '(" Bin")')
    call FileWriteNoAdvance ( this%iounit_dcp )
    write( IOBuffer, '(" Position")')
    call FileWriteNoAdvance ( this%iounit_dcp )
    do j=1,this%NComponents
      pc => this%Component(j)
      write( IOBuffer, '(" rho_", A)') trim( pc%PotModFileName )
      call FileWriteNoAdvance( this%iounit_dcp )
#if OSMOP == 2
      if ( pc%ChemPotMethod .eq. ChemPotMethodWidom ) then
        write( IOBuffer, '(" mueAvg_", A, " mueVar_", A)') trim( pc%PotModFileName ), trim( pc%PotModFileName )
        call FileWriteNoAdvance ( this%iounit_dcp )
      end if    
    end do
    write( IOBuffer, '(" PressureAvg PressureVar")')
    call FileWriteNoAdvance ( this%iounit_dcp )
#else
    end do  
#endif
    call FileWriteBlank( this%iounit_dcp ) 

    ! Update profile file
    do i = 1, NBinsDen
      write( IOBuffer, '( I3, F8.4)') i, real(i)/NBinsDen
      call FileWriteNoAdvance ( this%iounit_dcp )
      do j = 1, this%NComponents
        pc => this%Component(j)
        write( IOBuffer, '(F10.4)') pc%SumDenProfile(i)%Average * UnitDensity
        call FileWriteNoAdvance( this%iounit_dcp )
#if OSMOP == 2
        if ( pc%ChemPotMethod .eq. ChemPotMethodWidom ) then
          Variance = pc%SumChemPotProfile(i)%Variance / pc%SumChemPotProfile(i)%Average
          Average = log( pc%SumDenProfile(i)%Average / pc%SumChemPotProfile(i)%Average )
          write( IOBuffer, '(F10.4, F10.4)') Average, Variance
          call FileWriteNoAdvance ( this%iounit_dcp )
        end if    
      end do
      Average = this%SumPressureProfile(i)%Average * UnitPressure * 1E-6_RK
      Variance = this%SumPressureProfile(i)%Variance * UnitPressure * 1E-6_RK
      write( IOBuffer, '(F11.4, F10.4)') Average, Variance 
      call FileWriteNoAdvance( this%iounit_dcp )
#else
      end do  
#endif
      call FileWriteBlank( this%iounit_dcp )
    end do
    call FileWriteBlank( this%iounit_dcp )

    call FileClose( this%iounit_dcp )

  end subroutine TEnsemble_ProfileUpdate


!==============================================================!
!  Subroutine TEnsemble_ProfileClose                           !
!==============================================================!

  subroutine TEnsemble_ProfileClose( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    write( IOBuffer, '("##")' )
    call FileWrite( this%iounit_dcp )
    call FileClose( this%iounit_dcp )

  end subroutine TEnsemble_ProfileClose
#endif

  
!==============================================================!
!  Subroutine TEnsemble_RDFOpen                                !
!==============================================================!

  subroutine TEnsemble_RDFOpen( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer                   :: i, j, s, t

    ! initialize RDFSum
    do i=1, this%NComponents
      do j=1, this%NComponents
        do s=1, this%component(i)%molecule%NLJ126
          do t=1, this%component(j)%molecule%NLJ126
            this%Interaction(i,j)%PotLJ126LJ126(s, t)%RDFSum(:) = 0
          end do
        end do
      end do
    end do

    ! Open visualization file
    write( IOBuffer, '(I16)' ) this%EnsembleNumber
    call FileRewrite( this%iounit_rdf, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//RDFFileExtension )
    call FileWriteBlank( this%iounit_rdf )
    call FileClose( this%iounit_rdf )

  end subroutine TEnsemble_RDFOpen


!==============================================================!
!  Subroutine TEnsemble_RDFUpdate                              !
!==============================================================!

  subroutine TEnsemble_RDFUpdate( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    integer  :: i, j, s, t, o
    real(RK) :: RDFRho, RDFRhoLocal

    ! Calculate RDF
    do i= 1, this%NComponents
      do j= i, this%NComponents
        call Get_RDF( this%Interaction(i,j), this%RDFdr/this%BoxLength )
      end do
    end do

    ! Open RDF file
    write( IOBuffer, '(I16)' ) this%EnsembleNumber
    call FileRewrite( this%iounit_rdf, trim( OutputNameTag )//'_'//trim( adjustl( IOBuffer ) )//RDFFileExtension )
    write(IOBuffer, '(T5," r [A]")')
    call FileWriteNoAdvance( this%iounit_rdf )
    do i= 1, this%NComponents
      do j= i, this%NComponents
        do s=1, this%Component(i)%molecule%NLJ126
          do t=1, this%Component(j)%molecule%NLJ126
            write(IOBuffer, '(I5,I5)') i, j
            call FileWriteNoAdvance( this%iounit_rdf )
          end do
        end do            
      end do
    end do
    call FileWriteBlank( this%iounit_rdf )
    write(IOBuffer, '(T5,"______")')
    call FileWriteNoAdvance( this%iounit_rdf )
 
    do i= 1, this%NComponents
      do j= i, this%NComponents
        do s=1, this%Component(i)%molecule%NLJ126
          do t=1, this%Component(j)%molecule%NLJ126 
            write(IOBuffer, '(I5,I5)') s, t
            call FileWriteNoAdvance( this%iounit_rdf )
          end do
        end do
      end do
    end do
    call FileWriteBlank( this%iounit_rdf )

    do o = 1, RDFNumberShells
      write(IOBuffer, '(F10.4)') (o*this%RDFdr*UnitLength/Angstroem)
      call FileWriteNoAdvance( this%iounit_rdf )
      do i= 1, this%NComponents
        do j= i, this%NComponents
          do s=1, this%Component(i)%molecule%NLJ126
            do t=1, this%Component(j)%molecule%NLJ126
              RDFRho = this%SumDensity%Average  * this%Component(j)%Fraction  
              if (i == j) then
                RDFRhoLocal = 2.0 * real(this%Interaction(i,j)%PotLJ126LJ126(s,t)%RDFSum(o),RK) & 
&                                       / (this%RDFVSchale(o) * ((Step-1)/RDFUpdateFrequency + 1) * this%Component(i)%NPart)
              else
               RDFRhoLocal = real(this%Interaction(i,j)%PotLJ126LJ126(s,t)%RDFSum(o),RK) & 
&                                 / (this%RDFVSchale(o) * ((Step-1)/RDFUpdateFrequency + 1) * this%Component(i)%NPart)
              end if
              this%RDFValue(o) = RDFRhoLocal / RDFRho  
              write(IOBuffer, '(F10.4)') this%RDFValue(o)
              call FileWriteNoAdvance( this%iounit_rdf )
            end do
          end do
        end do
      end do
      call FileWriteBlank( this%iounit_rdf )
    enddo

    ! Close RDF file
    call FileClose( this%iounit_rdf )

  end subroutine TEnsemble_RDFUpdate


!==============================================================!
!  Subroutine TEnsemble_RDFClose                               !
!==============================================================!

  subroutine TEnsemble_RDFClose( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Close visualization file
    write( IOBuffer, '("##")' )
    call FileWrite( this%iounit_rdf )
    call FileClose( this%iounit_rdf )

  end subroutine TEnsemble_RDFClose


!==============================================================!
!  Subroutine TEnsemble_RestartSave                            !
!==============================================================!

  subroutine TEnsemble_RestartSave( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    type(TComponent), pointer :: pc
    integer                   :: i
#if TRANS ==1
    integer                   :: j, k, Mindex, StepCorr
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
      write(iounit_restart, '(2I10)' ) this%NResizeAttempts, this%NResizeSuccesses

      if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA ) then
        write(iounit_restart, '(2I10)' ) this%NInsertAttempts, this%NInsertSuccesses
        write(iounit_restart, '(2I10)' ) this%NDeleteAttempts, this%NDeleteSuccesses
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
    call RestartSave( this%SumConfEnthalpy )
    call RestartSave( this%SumVolume )
    call RestartSave( this%SumVirial )
    call RestartSave( this%SumEPotInter )
    call RestartSave( this%SumEPotIntra )
    if (printIDF) then
      call RestartSave( this%SumEPotIntra_Bond )
      call RestartSave( this%SumEPotIntra_Angle )
      call RestartSave( this%SumEPotIntra_Dihedral )
      call RestartSave( this%SumEPotIntra_Nonbonded )
      call RestartSave( this%SumVirialIntra )
      call RestartSave( this%SumVirialInter )
    end if
    call RestartSave( this%SumdEpotdV )
    call RestartSave( this%Sumd2EpotdV2 )
    if( EnsembleType .eq. EnsembleTypeNVE .and. LongRange .eq. Rfield) then
      call RestartSave( this%SumHmU )
      call RestartSave( this%SumHmUm1)
      call RestartSave( this%SumHmUm2 )
      call RestartSave( this%SumHmUm3 )
      call RestartSave( this%SumHmUm1dUdV )
      call RestartSave( this%SumHmUm1dUdV2 )
      call RestartSave( this%SumHmUm1d2UdV2 )
      call RestartSave( this%SumHmUm2dUdV )
      call RestartSave( this%SumHmUm2dUdV2 )
      call RestartSave( this%SumHmUm2d2UdV2 )
      call RestartSave( this%SumHmUm3dUdV )
      call RestartSave( this%SumHmUm3dUdV2 )
    end if

    if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA ) then
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
    call RestartSave( this%SumEPotCubic )
    call RestartSave( this%SumdEpotdVSquared )
    call RestartSave( this%SumEPotdEpotdV )
    call RestartSave( this%SumEPotSquareddEpotdV )
    call RestartSave( this%SumEPotdEpotdVSquared )
    call RestartSave( this%SumEPotd2EpotdV2 )

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
    if( LongRange .eq. Rfield) then
      if ( EnsembleType .eq. EnsembleTypeNVT ) then
        call RestartSave( this%SumA10resI )
        call RestartSave( this%SumA01resI )
        call RestartSave( this%SumA20resI )
        call RestartSave( this%SumA11resI )
        call RestartSave( this%SumA02resI )
        call RestartSave( this%SumA30resI )
        call RestartSave( this%SumA21resI )
        call RestartSave( this%SumA12resI )
      elseif ( EnsembleType .eq. EnsembleTypeNVE ) then
        call RestartSave( this%SumA10resI )
        call RestartSave( this%SumA01resI )
        call RestartSave( this%SumA20resI )
        call RestartSave( this%SumA11resI )
        call RestartSave( this%SumA02resI )
        call RestartSave( this%SumA30resI )
        call RestartSave( this%SumA21resI )
        call RestartSave( this%SumA12resI )
        call RestartSave( this%SumA10resII )
        call RestartSave( this%SumA01resII )
        call RestartSave( this%SumA20resII )
        call RestartSave( this%SumA11resII )
        call RestartSave( this%SumA02resII )
        call RestartSave( this%SumA30resII )
        call RestartSave( this%SumA21resII )
        call RestartSave( this%SumA12resII )
      end if
    end if

    ! 4.) Chemical potential and partial molar volumes
    do i = 1, this%NRealComponents
      pc => this%Component(i)
      select case( pc%ChemPotMethod )
      case( ChemPotMethodGradIns )
        call RestartSave( pc%SumInvChemPotRho )
        call RestartSave( pc%SumInvChemPot )
      case( ChemPotMethodWidom )
        call RestartSave( pc%SumChemPotV )
        call RestartSave( pc%SumChemPotVV )
        call RestartSave( pc%SumHW_counter )
        call RestartSave( pc%SumHW_denom )
      case( ChemPotMethodThermoInt )
        call RestartSave( pc%SumChemPotV )
        call RestartSave( pc%SumChemPotThermoIntWidom )
        call RestartSave( pc%SumChemPotThermoIntWidomV )
        call RestartSave( pc%SumHW_counter )
        call RestartSave( pc%SumHW_denom )
      end select

      if( pc%ChemPotMethod .ne. ChemPotMethodNone .and. ConstantPressure .and. this%NRealComponents > 1 ) then
        call RestartSave( pc%SumVW )
       call RestartSave( pc%SumHM )
      end if
    end do

#if TRANS ==1
if( RootProc .and. this%CorrfunMode ) then

    !Aenderungen Koester, ASpan_CF Matrix wurde beim Restart nicht uebergeben
    !Reduced correlation steps
    StepCorr = (Step + this%NStepCorr -1) / this%NStepCorr

    !Calculate matrix indexes
    Mindex = mod(StepCorr, this%NCorr )
    if (Mindex .eq. 0) then
      Mindex = this%NCorr
    end if

    k=mod(Mindex,this%NSpanCF)
    this%a(:,Mindex + 1 - k:Mindex ) = this%A_SpanCF(:,1:k)
    !Aenderungen Koester

    write( iounit_restart, '(I10)' ) this%NCorr
    write( iounit_restart, '(I10)' ) this%Mmess

    do i = 1, 3*this%NPart
        do j = 1, this%NCorr
            write( iounit_restart, '(ES20.12E3)' )  this%a( i, j)
        end do
    end do

    do i = 1, this%NCorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vsk(i,:)
    end do
    do i = 1, this%NCorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vsp(i,:)
    end do
    do i = 1, this%NCorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vbk(i,:)
    end do
    do i = 1, this%NCorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vbp(i,:)
    end do
    do i = 1, this%NCorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vckt(i,:)
    end do
    do i = 1, this%NCorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vckr(i,:)
    end do
    do i = 1, this%NCorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vcpt(i,:)
    end do
    do i = 1, this%NCorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vcpr(i,:)
    end do
    do i = 1, this%NCorr
        write( iounit_restart, '(3(ES20.12E3, :, ";"))' ) this%vcmt(i,:)
    end do
    do i = 1, this%NCorr
        write( iounit_restart, '(ES20.12E3)' ) this%average_cf_c(i)
    end do
    do i = 1, this%NCorr
        write( iounit_restart, '(ES20.12E3)' )  this%average_cf_vb(i)
    end do
    do i = 1, this%NCorr
        write( iounit_restart, '(ES20.12E3)' )  this%average_cf_vs(i)
    end do
    do i = 1, this%NCorr
        write( iounit_restart, '(ES20.12E3)' ) this%average_cf_ec(i)
    end do

    if (this%NComponents==2) then
      do i = 1, this%NCorr
        write( iounit_restart, '(ES20.12E3)' ) this%average_cf_soret(i)
      end do
    end if

    do i = 1, this%NComponents
      do j = 1, this%NCorr
        write( iounit_restart, '(ES20.12E3)' )  this%average_cf_d(i , j)
      end do
    end do

    if (this%NComponents > 1) then
      do i = 1, this%NComponents*this%NComponents
        do j = 1, this%NCorr
          write( iounit_restart, '(ES20.12E3)' ) this%average_lamda(i , j)
        end do
      end do
    end if

    write( iounit_restart, '(I10)' ) NBlocksMaxCF

    do i = 1, this%NComponents
      call RestartSave( this%Sumself_i(i), .true. )
    end do

     if(this%NComponents > 1) then
      do i = 1, this%NComponents
         do j = 1, this%NComponents
           call RestartSave( this%SumOnsager(i,j), .true. )
         end do
      end do
    end if

    call RestartSave( this%SumVisco_s, .true. )
    call RestartSave( this%SumVisco_b, .true. )
    call RestartSave( this%SumConduct, .true. )
    call RestartSave( this%SumSoret,   .true. )
    call RestartSave( this%SumEConduct,.true. )

    do i = 1,3
      write( iounit_restart, '(ES20.12E3)' )  this%sp(i)
    end do

    do i = 1,3
      write( iounit_restart, '(ES20.12E3)' )  this%sc(i)
    end do

endif
#endif

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
    integer                   :: i,j,t,stat,counter,k
#if TRANS==1
    integer                   :: Mindex,StepCorr
#endif
    real(RK)                  :: dummy, Factor

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
        read(iounit_restart, '(2I10)' ) this%NResizeAttempts, this%NResizeSuccesses
        if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA ) then
          read(iounit_restart, '(2I10)' ) this%NInsertAttempts, this%NInsertSuccesses
          read(iounit_restart, '(2I10)' ) this%NDeleteAttempts, this%NDeleteSuccesses
        end if
      end if
      read( iounit_restart, '(I10)' ) this%NRCutoffMax

    end if

#if MPI_VER > 0
    call MPI_Bcast( this%NPart, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
    call MPI_Bcast( this%Volume0, 1, MPI_RK, NRootProc, Communicator, ierror )

    if( (SimulationType .eq. MonteCarlo) .or. (SimulationType .eq. Gibbs) ) then
      call MPI_Bcast( this%DispVol, 1, MPI_RK, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%NResizeAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%NResizeSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )

      if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA ) then
        call MPI_Bcast( this%NInsertAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NInsertSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NDeleteAttempts, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
        call MPI_Bcast( this%NDeleteSuccesses, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
      end if
    end if

    call MPI_Bcast( this%NRCutoffMax, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
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
    call RestartRead( this%SumConfEnthalpy )
    call RestartRead( this%SumVolume )
    call RestartRead( this%SumVirial )
    call RestartRead( this%SumEPotInter )
    call RestartRead( this%SumEPotIntra )
    if (printIDF) then
      call RestartRead( this%SumEPotIntra_Bond )
      call RestartRead( this%SumEPotIntra_Angle )
      call RestartRead( this%SumEPotIntra_Dihedral )
      call RestartRead( this%SumEPotIntra_Nonbonded )
      call RestartRead( this%SumVirialIntra )
      call RestartRead( this%SumVirialInter )
    end if
    call RestartRead( this%SumdEpotdV )
    call RestartRead( this%Sumd2EpotdV2 )
    if( EnsembleType .eq. EnsembleTypeNVE .and. LongRange .eq. Rfield) then
      call RestartRead( this%SumHmU )
      call RestartRead( this%SumHmUm1)
      call RestartRead( this%SumHmUm2 )
      call RestartRead( this%SumHmUm3 )
      call RestartRead( this%SumHmUm1dUdV )
      call RestartRead( this%SumHmUm1dUdV2 )
      call RestartRead( this%SumHmUm1d2UdV2 )
      call RestartRead( this%SumHmUm2dUdV )
      call RestartRead( this%SumHmUm2dUdV2 )
      call RestartRead( this%SumHmUm2d2UdV2 )
      call RestartRead( this%SumHmUm3dUdV )
      call RestartRead( this%SumHmUm3dUdV2 )
    end if

    if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA ) then
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
    call RestartRead( this%SumEPotCubic )
    call RestartRead( this%SumdEpotdVSquared )
    call RestartRead( this%SumEPotdEpotdV )
    call RestartRead( this%SumEPotSquareddEpotdV )
    call RestartRead( this%SumEPotdEpotdVSquared )
    call RestartRead( this%SumEPotd2EpotdV2 )

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
    if( LongRange .eq. Rfield) then
      if ( EnsembleType .eq. EnsembleTypeNVT ) then
        call RestartRead( this%SumA10resI )
        call RestartRead( this%SumA01resI )
        call RestartRead( this%SumA20resI )
        call RestartRead( this%SumA11resI )
        call RestartRead( this%SumA02resI )
        call RestartRead( this%SumA30resI )
        call RestartRead( this%SumA21resI )
        call RestartRead( this%SumA12resI )
      elseif ( EnsembleType .eq. EnsembleTypeNVE ) then
        call RestartRead( this%SumA10resI )
        call RestartRead( this%SumA01resI )
        call RestartRead( this%SumA20resI )
        call RestartRead( this%SumA11resI )
        call RestartRead( this%SumA02resI )
        call RestartRead( this%SumA30resI )
        call RestartRead( this%SumA21resI )
        call RestartRead( this%SumA12resI )
        call RestartRead( this%SumA10resII )
        call RestartRead( this%SumA01resII )
        call RestartRead( this%SumA20resII )
        call RestartRead( this%SumA11resII )
        call RestartRead( this%SumA02resII )
        call RestartRead( this%SumA30resII )
        call RestartRead( this%SumA21resII )
        call RestartRead( this%SumA12resII )
      end if
    end if

    ! 4.) Chemical potential and partial molar volumes
    counter = this%NRealComponents+1
    do i = 1, this%NRealComponents
      pc => this%Component(i)
      select case( pc%ChemPotMethod )
      case( ChemPotMethodGradIns )
        call RestartRead( pc%SumInvChemPotRho )
        call RestartRead( pc%SumInvChemPot )
        pc%NFluctState = 0
        do j = counter,counter + this%Component(i)%Molecule%NFluct-1
          if (this%Component(j)%NPart .eq. 1) pc%NFluctState = j-counter+1
        end do
        counter = counter + this%Component(i)%Molecule%NFluct

      case( ChemPotMethodWidom )
        call RestartRead( pc%SumChemPotV )
        call RestartRead( pc%SumChemPotVV )
        call RestartRead( pc%SumHW_counter )
        call RestartRead( pc%SumHW_denom )

      case( ChemPotMethodThermoInt )
        call RestartRead( pc%SumChemPotV )
        call RestartRead( pc%SumChemPotThermoIntWidom )
        call RestartRead( pc%SumChemPotThermoIntWidomV )
        call RestartRead( pc%SumHW_counter )
        call RestartRead( pc%SumHW_denom )
      end select
      if( pc%ChemPotMethod .ne. ChemPotMethodNone .and. ConstantPressure .and. this%NRealComponents > 1 ) then
        call RestartRead( pc%SumVW )
        call RestartRead( pc%SumHM )
      end if
    end do

#if TRANS ==1
    if( this%CorrfunMode ) then
      if( RootProc ) then
        read( iounit_restart, '(I10)' ) this%NCorr
        read( iounit_restart, '(I10)' ) this%Mmess
      end if
    end if

#if MPI_VER > 0
    call MPI_Bcast( this%NCorr, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
    call MPI_Bcast( this%Mmess, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
#endif

    if ( RootProc ) then

      !Aenderungen Koester, ASpan_CF Matrix wurde beim Restart nicht uebergeben
      !Reduced correlation steps
      StepCorr = (Step + this%NStepCorr -1) / this%NStepCorr

      !Calculate matrix indexes
      Mindex = mod(StepCorr, this%NCorr )
      if (Mindex .eq. 0) then
       Mindex = this%NCorr
      end if

      k=mod(Mindex,this%NSpanCF)
      !Aenderungen Koester

      do i = 1, 3*this%NPart
        do j = 1, this%NCorr
          read( iounit_restart, '(ES20.12E3)' )  this%a( i, j)
        end do
      end do

      !Aenderungen Koester
      this%A_SpanCF(:,1:k) =  this%a(:,Mindex + 1 - k:Mindex )
      !Aenderungen Koester

      do i = 1, this%NCorr
        read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%vsk(i,:)
      end do
      do i = 1, this%NCorr
        read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%vsp(i,:)
      end do
      do i = 1, this%NCorr
        read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%vbk(i,:)
      end do
      do i = 1, this%NCorr
        read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%vbp(i,:)
      end do
      do i = 1, this%NCorr
        read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%vckt(i,:)
      end do
      do i = 1, this%NCorr
        read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%vckr(i,:)
      end do
      do i = 1, this%NCorr
        read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%vcpt(i,:)
      end do
      do i = 1, this%NCorr
        read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%vcpr(i,:)
      end do
      do i = 1, this%NCorr
        read( iounit_restart, '(3(ES20.12E3, :, X))' ) this%vcmt(i,:)
      end do
      do i = 1, this%NCorr
        read( iounit_restart, '(ES20.12E3)' ) this%average_cf_c(i)
      end do
      do i = 1, this%NCorr
        read( iounit_restart, '(ES20.12E3)' )  this%average_cf_vb(i)
      end do
      do i = 1, this%NCorr
        read( iounit_restart, '(ES20.12E3)' )  this%average_cf_vs(i)
      end do
      do i = 1, this%NCorr
        read( iounit_restart, '(ES20.12E3)' )  this%average_cf_ec(i)
      end do

      if (this%NComponents==2) then
        do i = 1, this%NCorr
          read( iounit_restart, '(ES20.12E3)' ) this%average_cf_soret(i)
        end do
      end if

      do i = 1, this%NComponents
        do j = 1, this%NCorr
          read( iounit_restart, '(ES20.12E3)' )  this%average_cf_d(i , j)
        end do
      end do
      if (this%Ncomponents>1) then
        do i = 1, this%NComponents*this%NComponents
          do j = 1, this%NCorr
            read( iounit_restart, '(ES20.12E3)' )  this%average_lamda(i , j)
          end do
        end do
      end if

      read( iounit_restart, '(I10)' ) NBlocksMaxCF

      do i = 1, this%NComponents
      call RestartRead( this%Sumself_i(i) )
      end do

      if (this%NComponents > 1) then
        do i = 1, this%NComponents
           do j = 1, this%NComponents
             call RestartRead( this%SumOnsager(i,j) )
           end do
        end do
      end if

      call RestartRead( this%SumVisco_s )
      call RestartRead( this%SumVisco_b )
      call RestartRead( this%SumConduct )
      call RestartRead( this%SumSoret   )
      call RestartRead( this%SumEConduct)

      do i = 1,3
        read( iounit_restart, '(ES20.12E3)' )  this%sp(i)
      end do

      do i = 1,3
        read( iounit_restart, '(ES20.12E3)' )  this%sc(i)
      end do

    end if
#endif

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
         nullify ( this%U_fourierLocal )
         nullify ( this%SSin )
         nullify ( this%SCos )
         nullify ( this%rold )
         nullify ( this%Vec2 )
         allocate(this%U_fourierLocal(this%BoxenAnzahlMax),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error U_fourier'
         allocate(this%SSin(this%BoxenAnzahlMax),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error SSin'
         allocate(this%SCos(this%BoxenAnzahlMax),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error SCos'
         allocate(this%Vec2(this%BoxenAnzahlMax),STAT=stat)
         if(stat >0) write(*,*) 'Allocation Error Vec2'
         allocate(this%VirIntra(this%NPartMax),STAT=stat)
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

#if SPME > 0
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
#endif
    else if ( LongRange .eq. Rodgers ) then
       do i=1,this%NComponents
         do j=1,this%NComponents
           this%Interaction(i,j)%Kappa = this%Kappa
         end do
       end do
    end if

    ! Reading and broadcasting thi-file for ThermoInt
      t = this%NRealComponents+1
      write( IOBuffer, '(I16)' ) this%EnsembleNumber
      if ( any(this%Component(:)%ChemPotMethod .eq. ChemPotMethodThermoInt)) then
        call FileReset( this%iounit_thermoint, trim(OutputNameTag)//'_'//trim( adjustl(IOBuffer) )//ThermoIntFileExtension )
      end if
      do i=1,this%NRealComponents
        pc => this%Component(i)
        if (pc%ChemPotMethod .eq. ChemPotMethodThermoInt) then

          call FileReadParameter( this%Component(t)%lambda, this%iounit_thermoint , "currentlambda", .false. )
          pc%CalcChemPot = .true.
          if (RootProc) then
          !read empty line
          read( this%iounit_thermoint, * )

          ! read thermoint-profile
          do j = 0,pc%NBins-1
            read( this%iounit_thermoint, '(I6," ", F6.3,7(" ", E15.6)," ", I10)' )  k, dummy, pc%BinsEn(j), pc%BinsdEndLa(j), pc%BinsdEndLaV(j), pc%BinsdEndLaH(j), pc%BinsIntdEndLa(j), pc%BinsIntVW(j), pc%BinsIntHW(j), pc%BinsVisit(j)
          end do
          end if
          t = t+1
        end if 
      end do
      if ( any(this%Component(:)%ChemPotMethod .eq. ChemPotMethodThermoInt)) then
        call FileClose( this%iounit_thermoint )
      end if

#if MPI_VER > 0
    if (SimulationType .eq. MonteCarlo) then
      t = this%NRealComponents+1
      do i=1,this%NRealComponents
        pc => this%Component(i)
        if (pc%ChemPotMethod .eq. ChemPotMethodThermoInt) then
          !call MPI_Bcast( this%Component(t)%lambda, 1, MPI_RK, NRootProc, Communicator, ierror ) //done during the preceding call FileReadParameter 
          !The Broadcast of the following properties would not have been necessary for MD runs during the implementation of ThermoInt (they may be however important for future use)
          call MPI_Bcast( pc%BinsEn(0:pc%NBins-1), size( pc%BinsEn ), MPI_RK, NRootProc, Communicator, ierror )
          call MPI_Bcast( pc%BinsdEndLa(0:pc%NBins-1), size( pc%BinsdEndLa ), MPI_RK, NRootProc, Communicator, ierror )
          call MPI_Bcast( pc%BinsdEndLaV(0:pc%NBins-1), size( pc%BinsdEndLaV ), MPI_RK, NRootProc, Communicator, ierror )
          call MPI_Bcast( pc%BinsdEndLaH(0:pc%NBins-1), size( pc%BinsdEndLaH ), MPI_RK, NRootProc, Communicator, ierror )
          call MPI_Bcast( pc%BinsIntdEndLa(0:pc%NBins-1), size( pc%BinsIntdEndLa ), MPI_RK, NRootProc, Communicator, ierror )
          call MPI_Bcast( pc%BinsIntVW(0:pc%NBins-1), size( pc%BinsIntVW ), MPI_RK, NRootProc, Communicator, ierror )
          call MPI_Bcast( pc%BinsIntHW(0:pc%NBins-1), size( pc%BinsIntHW ), MPI_RK, NRootProc, Communicator, ierror )
          call MPI_Bcast( pc%BinsVisit(0:pc%NBins-1), size( pc%BinsVisit ), MPI_INTEGER, NRootProc, Communicator, ierror )
          t = t+1
        endif
      enddo
    end if
#endif
    t = this%NRealComponents+1
    do i=1,this%NRealComponents
      if (this%Component(i)%ChemPotMethod .eq. ChemPotMethodThermoInt) then
        Factor = this%Component(t)%lambda**this%Component(i)%LambdaExponent
        call ScaleInteractionThermoInt(this, t, Factor)
        t = t+1
      endif
    enddo
    ! End of reading and broadcasting thi-file for ThermoInt

    if( SimulationType .eq. MolecularDynamics ) then

      ! Calculate temperature and set up ndf
      call CalculateEKin( this, .false. )

    else

      ! Set temperature
      this%Temperature = this%RefTemperature

      ! Initialize energy matrix
      call Unit2Atom( this )
      call Energy( this, this%EPot )
      call UpdateEnergy( this )

    end if

  end subroutine TEnsemble_RestartRead



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
    integer                   :: nc, i, i1, i2, j, n, n2, n3, k, nu
    real(RK)                  :: C(this%NPart*this%NUnitMax* 3), Q(this%NPart*this%NUnitMax*4)

    if( .not. RootProc ) return

    nu = this%NUnitMax
    n = this%NPart*nu
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
        do k = 1, pc%Molecule%NUnit
          C(i*nu+k) = pc%P0(j, 1, k)
          C(i*nu+k+n) = pc%P0(j, 2, k)
          C(i*nu+k+n2) = pc%P0(j, 3, k)
          Q(i*nu+k) = pc%Q0(j, 1, k)
          Q(i*nu+k+n) = pc%Q0(j, 2, k)
          Q(i*nu+k+n2) = pc%Q0(j, 3, k)
          Q(i*nu+k+n3) = pc%Q0(j, 4, k)
        end do
      end do
      i1 = i2 + 1
    end do

    open(996, file='EIN_1', status='REPLACE', action='WRITE')
    write(996, '(3ES22.12E3)') C
    write(996, '(4ES22.12E3)') Q
    write(996,*) this%Component(1)%NFLuctMax - this%Component(1)%NFluctState + 1
    write(996,*) ix, iy
    close(996)


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

    ! Declare arguments
    type(TEnsemble)         :: this
    ! Declare local variables
    integer :: i,Si,Sj
    integer :: NX,NY,NZ
    real(RK):: USelbstTermKomp
    real(RK):: Faktor, NSQ
    real(RK):: RXi,RYi,RZi,RXj,RYj,RZj
    real(RK):: drxij,dryij,drzij,dr
    real(RK):: approx
    real(RK):: UIntraTermKomp
    real(RK):: twopi
    real(RK):: test

    twopi = 2._RK*PI

    i = 0
    DO NX = 0, this%NMAX, 1
       IF (NX .EQ.0) THEN
         Faktor = 1.0_RK
       ELSE
         Faktor = 2.0_RK
       END IF
       DO NY = -this%NMAX,this%NMAX,1
         DO NZ = -this%NMAX,this%NMAX,1
           NSQ = NX*NX + NY*NY + NZ*NZ
           IF ( (NSQ .LE. this%NSQMAX) .AND. (NSQ .NE. 0) ) THEN
             i = i+1
             IF (i .GT. this%NVecMax) STOP   'BoxAnzahl zu gross'
             this%Ewald_Prefac(i) = Faktor*exp(-(PI/this%KappaL)**2 * NSQ) / (Pi2 * NSQ * this%BoxLength)
             this%Ewald_Vec(1,i)  = twopi*NX
             this%Ewald_Vec(2,i)  = twopi*NY
             this%Ewald_Vec(3,i)  = twopi*NZ
           END IF
         END DO
       END DO
    END DO
    this%BoxenAnzahlMax = i

#if MPI_VER > 0
    this%NBox2 = i
    if ( (SimulationType .ne. MonteCarlo) .or. (Equilibration .and. CommonEqui) ) then
      call MPI_Bcast( this%Ewald_Vec(1,:), this%BoxenAnzahlMax, MPI_RK, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%Ewald_Vec(2,:), this%BoxenAnzahlMax, MPI_RK, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%Ewald_Vec(3,:), this%BoxenAnzahlMax, MPI_RK, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%Ewald_Prefac, this%BoxenAnzahlMax, MPI_RK, NRootProc, Communicator, ierror )
    end if
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

    this%USelbstTerm = -this%USelbstTerm * this%Kappa / sqrt(Pi)

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

          if (dr .ge. 0.0000001) then
            call ErrorApprox(this%interaction(i,i)%PotChargeCharge(Si,Sj), this%Kappa*dr, approx)

            UIntraTermKomp = UIntraTermKomp - this%Component(i)%Molecule%SiteCharge(Si)%e* &
&                            this%Component(i)%Molecule%SiteCharge(Sj)%e / dr * (1-approx)
          end if
        END DO
      END DO
      this%UIntra = this%UIntra + this%component(i)%NPart * UIntraTermKomp
    END DO

#if MPI_VER > 0
    this%NBox1  = ProcRange ( this%BoxenAnzahlMax, this%NBox0, this%NBox2 )
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

   real(RK),pointer, contiguous:: RX(:),RY(:),RZ(:)
   real(RK),pointer, contiguous:: PX(:),PY(:),PZ(:)
   real(RK),pointer, contiguous:: FX(:),FY(:),FZ(:)
   real(RK),pointer:: q(:)
   real(RK) :: RXloc(this%NPartMax),RYloc(this%NPartMax),RZloc(this%NPartMax)
   real(RK) :: PXloc(this%NPartMax),PYloc(this%NPartMax),PZloc(this%NPartMax)

   real(RK):: KVec(3), KVec05
   real(RK):: EPotLocal
   real(RK):: Viriallocal, VirIntra, VirIntraLocal
   real(RK):: SSinSum,SCosSum
   real(RK):: KappaL2, vorfac
   real(RK):: Facx,Facy,Facz
   real(RK):: BoxLength
   real(RK):: Faktor(1:this%NPart)
   real(RK):: HFac(1:this%NPart)
   real(RK):: HFacX(1:this%NPart), HFacY(1:this%NPart), HFacZ(1:this%NPart)
   real(RK):: distx(1:this%NPart),disty(1:this%NPart),distz(1:this%NPart)
   real(RK),allocatable:: qsinfac(:)
   real(RK),allocatable:: qcosfac(:)
   !real(RK):: qsinfac(1:this%NPart*this%NComponents*5)
   !real(RK):: qcosfac(1:this%NPart*this%NComponents*5)
   integer :: ChargeNumber
#if MPI_VER > 0
   integer, pointer:: i0,i1
#endif
#if  TRANS == 1
   real(RK),pointer, contiguous:: VSx(:),VSy(:),VSz(:)
   real(RK),pointer, contiguous:: VSux(:),VSuy(:),VSuz(:)
   real(RK),pointer, contiguous:: VBx(:),VBy(:),VBz(:)
   real(RK)        :: multiplicator
   real(RK)        :: Contrib
  
#endif

   type(TMolecule), pointer               :: mol
   integer:: stat

   EPotLocal = 0.0_RK
   VirialLocal = 0.0_RK
   KappaL2 = 1.0_RK/(2._RK*this%KappaL**2)
   BoxLength = this%BoxLength
   vorfac = 2._RK/BoxLength

   m = 0
   do i=1,this%NComponents
     m = m + this%Component(i)%NPart*this%Component(i)%Molecule%NCharge
   enddo
   allocate(qsinfac(1:m))
   allocate(qcosfac(1:m))
   qsinfac(1:m)=0
   qcosfac(1:m)=0

   ! Virial
   this%Vec2 = this%Ewald_Vec(1,:)**2 + this%Ewald_Vec(2,:)**2 + this%Ewald_Vec(3,:)**2 
   VirIntraLocal = 0._RK

   this%SSin=0._RK
   this%SCos=0._RK

#if MPI_VER > 0
   i0 => this%NBox0
   i1 => this%NBox2
   DO i=i0,i1,1
#else
   DO i=1,this%BoxenAnzahlMax,1
# endif

     ChargeNumber = 0
     KVec = this%Ewald_Vec(:,i)
     KVec05 = sum(KVec)/2._RK
     SSinSum = 0._RK
     SCosSum = 0._RK
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

! Alte Version
         RXloc(1:molec) = RX(1:molec)
         RYloc(1:molec) = RY(1:molec)
         RZloc(1:molec) = RZ(1:molec)
         PXloc(1:molec) = PX(1:molec)
         PYloc(1:molec) = PY(1:molec)
         PZloc(1:molec) = PZ(1:molec)

         DO m=1,molec
           if (PX(m) < 0) RXloc(m) = RXloc(m) + 1._RK
           if (PY(m) < 0) RYloc(m) = RYloc(m) + 1._RK
           if (PZ(m) < 0) RZloc(m) = RZloc(m) + 1._RK
           if (PX(m) < 0) PXloc(m) = PXloc(m) + 1._RK
           if (PY(m) < 0) PYloc(m) = PYloc(m) + 1._RK
           if (PZ(m) < 0) PZloc(m) = PZloc(m) + 1._RK
         end DO
         Faktor(1:molec) = KVec(1) * RXloc + KVec(2)*RYloc + KVec(3)*RZloc
         qsinfac(ChargeNumber+1:ChargeNumber+molec) = q(l)*sin(Faktor)
         qcosfac(ChargeNumber+1:ChargeNumber+molec) = q(l)*cos(Faktor)
         SSinSum = SSinSum + sum(qsinfac(ChargeNumber+1:ChargeNumber+molec))
         SCosSum = SCosSum + sum(qcosfac(ChargeNumber+1:ChargeNumber+molec))

         ChargeNumber = ChargeNumber + molec
       END DO
     END DO

     this%SSin(i) = SSinSum
     this%SCos(i) = SCosSum

! Forces
     Facx = KVec(1)*this%Ewald_Prefac(i)*vorfac
     Facy = KVec(2)*this%Ewald_Prefac(i)*vorfac
     Facz = KVec(3)*this%Ewald_Prefac(i)*vorfac

     ChargeNumber = 0
     DO j=1,this%NComponents,1
       mol => this%Component(j)%Molecule
       molec = this%Component(j)%NPart
       DO l=1,mol%NCharge
         RX => this%Component(j)%Molecule%SiteCharge(l)%RX(1:molec)
         RY => this%Component(j)%Molecule%SiteCharge(l)%RY(1:molec)
         RZ => this%Component(j)%Molecule%SiteCharge(l)%RZ(1:molec)
         PX => this%Component(j)%Molecule%SiteCharge(l)%PX(1:molec)
         PY => this%Component(j)%Molecule%SiteCharge(l)%PY(1:molec)
         PZ => this%Component(j)%Molecule%SiteCharge(l)%PZ(1:molec)
         FX => this%Component(j)%Molecule%SiteCharge(l)%FX
         FY => this%Component(j)%Molecule%SiteCharge(l)%FY
         FZ => this%Component(j)%Molecule%SiteCharge(l)%FZ

         distx(1:molec) = (RX-PX)*BoxLength
         disty(1:molec) = (RY-PY)*BoxLength
         distz(1:molec) = (RZ-PZ)*BoxLength

         HFac(1:molec) = qsinfac(ChargeNumber+1:ChargeNumber+molec)*SCosSum - qcosfac(ChargeNumber+1:ChargeNumber+molec)*SSinSum
         HFacX(1:molec) = Facx*HFac(1:molec)
         HFacy(1:molec) = Facy*HFac(1:molec)
         HFacZ(1:molec) = Facz*HFac(1:molec)
         
         FX = FX + HFacX(1:molec)
         FY = FY + HFacY(1:molec)
         FZ = FZ + HFacZ(1:molec)

         VirIntraLocal = VirIntraLocal + sum(HFacX(1:molec)*distx(1:molec) &
   &                                   +     HFacY(1:molec)*disty(1:molec) &
   &                                   +     HFacZ(1:molec)*distz(1:molec))

#if  TRANS == 1
         ! Preparation for Transport properties
         VSx  => mol%SiteCharge(l)%vsCx
         VSy  => mol%SiteCharge(l)%vsCy
         VSz  => mol%SiteCharge(l)%vsCz
         VBx  => mol%SiteCharge(l)%vbCx
         VBy  => mol%SiteCharge(l)%vbCy
         VBz  => mol%SiteCharge(l)%vbCz
       
          ! Intramolecular Forces
         VSx  = VSx  - HFacX(1:molec)*disty(1:molec)
         VSy  = VSy  - HFacX(1:molec)*distz(1:molec)
         VSz  = VSz  - HFacY(1:molec)*distz(1:molec)
         VBx  = VBx  - HFacX(1:molec)*distx(1:molec)
         VBy  = VBy  - HFacY(1:molec)*disty(1:molec)
         VBz  = VBz  - HFacZ(1:molec)*distz(1:molec)
       
#endif

         ChargeNumber = ChargeNumber + molec

       END DO
     END DO

#if  TRANS == 1
     ! Force contribution due to the intermolecular long-range forces
     ! Since these properties are specific for the entire solution and not valid for individual
     ! components, these quantities are calculated once and added onto the contribution of 
     ! the last component in the system
     multiplicator = (1._RK/this%Vec2(i) + 0.25_RK/this%KappaL**2)
     Contrib = -2._RK * multiplicator * this%Ewald_Prefac(i) * (SSinSum*SSinSum + SCosSum*SCosSum) / this%Volume0
     VSx(1) = VSx(1) + Contrib* (this%Ewald_Vec(1,i)*this%Ewald_Vec(2,i) )
     VSy(1) = VSy(1) + Contrib* (this%Ewald_Vec(1,i)*this%Ewald_Vec(3,i) )
     VSz(1) = VSz(1) + Contrib* (this%Ewald_Vec(2,i)*this%Ewald_Vec(3,i) )
     VBx(1) = VBx(1) + Contrib* (-0.5_RK/multiplicator + this%Ewald_Vec(1,i)*this%Ewald_Vec(1,i) )
     VBy(1) = VBy(1) + Contrib* (-0.5_RK/multiplicator + this%Ewald_Vec(2,i)*this%Ewald_Vec(2,i) )
     VBz(1) = VBz(1) + Contrib* (-0.5_RK/multiplicator + this%Ewald_Vec(3,i)*this%Ewald_Vec(3,i) )
#endif

   END DO ! Boxenschleife

! Finish Calculation
! STOP
! Energy
   this%U_fourierLocal = this%Ewald_Prefac * (this%SSin*this%SSin + this%SCos*this%SCos)


#if MPI_VER > 0
   call MPI_Reduce( sum(this%U_fourierLocal), EPotLocal, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
   call MPI_Reduce( EPotLocal - sum(this%U_fourierLocal*KappaL2*this%Vec2), VirialLocal, 1, MPI_RK, MPI_SUM, &
&                   NRootProc, Communicator, ierror )
   call MPI_Reduce( VirIntraLocal, VirIntra, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )

#else
   EPotLocal = sum(this%U_fourierLocal)
   VirialLocal = EPotLocal - sum(this%U_fourierLocal *KappaL2*this%Vec2)
   VirIntra = VirIntraLocal

#endif
   if( RootProc ) then
     this%UFourier= EPotLocal
     this%EVirial = -(Viriallocal - VirIntra)*Third
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
    integer        :: processes

    ! Declare local variables
    integer :: Si,Sj,i
    real(RK):: RXi,RYi,RZi,RXj,RYj,RZj
    real(RK):: drxij,dryij,drzij,dr
    real(RK):: approx
    real(RK):: UIntraTermKomp, USelbstTermKomp


    ! Setting of the right scaling for MC
    if (Equilibration .and. CommonEqui) then
      processes = NProcs
    else
      processes = 1
    end if


    ! Selbstterm
    this%USelbstTerm = 0.0_RK
    DO Sj=1,this%NComponents,1
       USelbstTermKomp = 0.0_RK
       DO Si=1,this%Component(Sj)%Molecule%NCharge,1
         USelbstTermKomp = USelbstTermKomp + this%Component(Sj)%Molecule%SiteCharge(Si)%e**2
       END DO
       this%USelbstTerm = this%USelbstTerm + this%Component(Sj)%NPart * USelbstTermKomp
    END DO

    this%USelbstTerm = -this%USelbstTerm * this%Kappa / sqrt(Pi) / processes


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

          if (dr .ge. 0.0000001) then
            call ErrorApprox(this%interaction(i,i)%PotChargeCharge(Si,Sj), this%Kappa*dr, approx)

            UIntraTermKomp = UIntraTermKomp - this%Component(i)%Molecule%SiteCharge(Si)%e* &
&                            this%Component(i)%Molecule%SiteCharge(Sj)%e / dr * (1-approx)
          end if
        END DO
      END DO
      this%UIntra = this%UIntra + this%component(i)%NPart * UIntraTermKomp / processes
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

   integer :: i,j,l,m
   integer :: molec
#if MPI_VER > 0
   integer,pointer :: i0
   integer,pointer :: i1
#endif

   real(RK),pointer, contiguous:: RX(:),RY(:),RZ(:)
   real(RK),pointer, contiguous:: PX(:),PY(:),PZ(:)
   real(RK),pointer:: q(:)
   real(RK) :: RXloc(this%NPartMax),RYloc(this%NPartMax),RZloc(this%NPartMax)
   real(RK) :: PXloc(this%NPartMax),PYloc(this%NPartMax),PZloc(this%NPartMax)

   real(RK):: KVec(3)
   real(RK):: EPotLocal
   real(RK):: Viriallocal,VirIntra, VirIntraLocal
   real(RK):: BoxLength
   real(RK):: SSinSum,SCosSum
   real(RK):: KappaL2
   real(RK):: vorfac
   real(RK):: facx,facy,facz
   real(RK):: Faktor(1:this%NPart)
   real(RK):: HFac(1:this%NPart)
   real(RK):: distx(1:this%NPart),disty(1:this%NPart),distz(1:this%NPart)
   real(RK),allocatable:: qsinfac(:)
   real(RK),allocatable:: qcosfac(:)
   !real(RK):: qsinfac(1:this%NPart*this%NComponents*5)
   !real(RK):: qcosfac(1:this%NPart*this%NComponents*5)
   integer :: ChargeNumber


   EPotLocal = 0._RK
   VirialLocal = 0._RK
   KappaL2 = 1.0_RK/(2._RK*this%KappaL**2)
   BoxLength = this%BoxLength
   vorfac = 2._RK / this%BoxLength

   m = 0
   do i=1,this%NComponents
     m = m + this%Component(i)%NPart*this%Component(i)%Molecule%NCharge
   enddo
   allocate(qsinfac(1:m))
   allocate(qcosfac(1:m))
   qsinfac(1:m)=0
   qcosfac(1:m)=0

   if (this%OptPressure) then
     VirIntraLocal = 0._RK
   end if

   this%SSin = 0._RK
   this%SCos = 0._RK
#if MPI_VER > 0
   i0 => this%NBox0
   i1 => this%NBox2
   DO i=i0,i1,1
#else
   DO i=1,this%BoxenAnzahlMax,1
#endif
     ChargeNumber = 0
     KVec = this%Ewald_Vec(:,i)
     SSinSum = 0._RK
     SCosSum = 0._RK
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

         RXloc(1:molec) = RX(1:molec)
         RYloc(1:molec) = RY(1:molec)
         RZloc(1:molec) = RZ(1:molec)
         PXloc(1:molec) = PX(1:molec)
         PYloc(1:molec) = PY(1:molec)
         PZloc(1:molec) = PZ(1:molec)

         DO m=1,molec
           if (PX(m) < 0) RXloc(m) = RXloc(m) + 1._RK
           if (PY(m) < 0) RYloc(m) = RYloc(m) + 1._RK
           if (PZ(m) < 0) RZloc(m) = RZloc(m) + 1._RK
           if (PX(m) < 0) PXloc(m) = PXloc(m) + 1._RK
           if (PY(m) < 0) PYloc(m) = PYloc(m) + 1._RK
           if (PZ(m) < 0) PZloc(m) = PZloc(m) + 1._RK
         end DO

         Faktor(1:molec) = KVec(1) * RXloc + KVec(2)*RYloc + KVec(3)*RZloc

         if ( this%OptPressure ) then
             qsinfac(ChargeNumber+1:ChargeNumber+molec) = q(l)*sin(Faktor)
             qcosfac(ChargeNumber+1:ChargeNumber+molec) = q(l)*cos(Faktor)
             SSinSum = SSinSum + sum(qsinfac(ChargeNumber+1:ChargeNumber+molec))
             SCosSum = SCosSum + sum(qcosfac(ChargeNumber+1:ChargeNumber+molec))
             ChargeNumber = ChargeNumber + molec

         else
             SSinSum = SSinSum + sum(q(l)*sin(Faktor))
             SCosSum = SCosSum + sum(q(l)*cos(Faktor))
         end if

       END DO
     END DO

     this%SSin(i) = SSinSum
     this%SCos(i) = SCosSum

! Forces
     if ( this%OptPressure ) then
       Facx = KVec(1)*this%Ewald_Prefac(i)*vorfac
       Facy = KVec(2)*this%Ewald_Prefac(i)*vorfac
       Facz = KVec(3)*this%Ewald_Prefac(i)*vorfac

       ChargeNumber = 0
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

           HFac(1:molec) = qsinfac(ChargeNumber+1:ChargeNumber+molec)*SCosSum - qcosfac(ChargeNumber+1:ChargeNumber+molec)*SSinSum

           distx(1:molec) = (RXloc - PXloc)*BoxLength
           disty(1:molec) = (RYloc - PYloc)*BoxLength
           distz(1:molec) = (RZloc - PZloc)*BoxLength
           VirIntraLocal = VirIntraLocal + sum(Facx*HFac(1:molec)*distx(1:molec)) &
&                                        + sum(Facy*HFac(1:molec)*disty(1:molec)) &
&                                        + sum(Facz*HFac(1:molec)*distz(1:molec))
           ChargeNumber = ChargeNumber + molec
         END DO
       END DO
     end if
   END DO ! Boxenschleife


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Finish Calculation

! Energy
   this%U_fourierLocal = this%Ewald_Prefac * (this%SSin*this%SSin + this%SCos*this%SCos)
! Virial
   this%Vec2 = this%Ewald_Vec(1,:)**2 + this%Ewald_Vec(2,:)**2 + this%Ewald_Vec(3,:)**2 

#if MPI_VER > 0
   if (Equilibration .and. CommonEqui) then
     ! Energy
      call MPI_Allreduce( sum(this%U_fourierLocal), EPotLocal, 1, MPI_RK, MPI_SUM, Communicator, ierror )

     ! Virial
      if (this%OptPressure) then
        call MPI_Allreduce( sum(this%U_fourierLocal) - sum(this%U_fourierLocal*KappaL2*this%Vec2), VirialLocal, 1, &
&                           MPI_RK, MPI_SUM, Communicator, ierror )
        call MPI_Allreduce( VirIntraLocal, VirIntra, 1, MPI_RK, MPI_SUM, Communicator, ierror )
        this%EVirial = -(Viriallocal - VirIntra)*Third / NProcs
      end if
     this%UFourier= EPotLocal / NProcs

   else
     ! Energy
     EPotLocal = sum(this%U_fourierLocal)
     this%UFourier= EPotLocal
     ! Virial
      if (this%OptPressure) then
        VirialLocal = EPotLocal - sum(this%U_fourierLocal *KappaL2*this%Vec2)
        VirIntra = VirIntraLocal
        this%EVirial = -(Viriallocal - VirIntra)*Third
      end if
   end if

#else
   ! Energy
   EPotLocal = sum(this%U_fourierLocal)
   this%UFourier= EPotLocal / NProcs
   ! Virial
   if (this%OptPressure) then
     VirialLocal = EPotLocal - sum(this%U_fourierLocal *KappaL2*this%Vec2)
     VirIntra = VirIntraLocal
     this%EVirial = -(Viriallocal - VirIntra)*Third
   end if
#endif

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
   real(RK):: RX2,RY2,RZ2
   real(RK),pointer:: q(:)
   real(RK):: KVec(3)
   real(RK):: EPotLocal
   real(RK)::SSin_Vec,SCos_Vec
   real(RK):: KappaL2
   real(RK):: sinfac,cosfac
   real(RK)::Faktor,Faktor2

   integer :: i,l
   integer,intent(in)::nc,np

#if MPI_VER > 0
   integer,pointer :: i0, i1
#endif


! Declarations
   KappaL2 = 1.0_RK/(2._RK*this%KappaL**2)

! Calculation
#if MPI_VER > 0
   i0 => this%NBox0
   i1 => this%NBox2
   DO i=i0,i1,1
#else
   DO i=1,this%BoxenAnzahlMax,1
#endif

     KVec = this%Ewald_Vec(:,i)
     mol => this%Component(nc)%Molecule
     q => mol%SiteCharge(1:mol%NCharge)%e
     SSin_Vec =0._RK
     SCos_Vec =0._RK
       DO l=1,mol%NCharge
         RX = this%Component(nc)%Molecule%SiteCharge(l)%RX(np)
         RY = this%Component(nc)%Molecule%SiteCharge(l)%RY(np)
         RZ = this%Component(nc)%Molecule%SiteCharge(l)%RZ(np)
         RX2 = this%rold(l,1)
         RY2 = this%rold(l,2)
         RZ2 = this%rold(l,3)

         if (RX < 0) RX = RX + 1._RK
         if (RY < 0) RY = RY + 1._RK
         if (RZ < 0) RZ = RZ + 1._RK
         if (RX2< 0) RX2= RX2+ 1._RK
         if (RY2< 0) RY2= RY2+ 1._RK
         if (RZ2< 0) RZ2= RZ2+ 1._RK

         Faktor = KVec(1) * RX + KVec(2)*RY + KVec(3)*RZ
         Faktor2 = KVec(1) * RX2 + KVec(2)*RY2 + KVec(3)*RZ2

         sinfac = q(l)*(sin(Faktor)-sin(Faktor2))
         cosfac = q(l)*(cos(Faktor)-cos(Faktor2))

         SSin_Vec = SSin_Vec + sinfac
         SCos_Vec = SCos_Vec + cosfac
       END DO

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
   if ( Equilibration .and. CommonEqui ) then
     call MPI_Allreduce( sum(this%U_fourierLocal), EPotLocal, 1, MPI_RK, MPI_SUM, Communicator, ierror )
     this%UFourier= EPotLocal / NProcs
   else
     EPotLocal = sum(this%U_fourierLocal)
     this%UFourier= EPotLocal
   end if

#else
   EPotLocal = sum(this%U_fourierLocal)
     this%UFourier= EPotLocal
#endif

  END subroutine TEnsemble_EwaldFourierEnergy1


!==============================================================!
!  Subroutine TSimulation_Ewald_FourierTermEnergyCF            !
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
   real(RK),pointer:: q2(:)
   real(RK):: KVec(3)
   real(RK):: EPotLocal
   real(RK)::SSin_Vec,SCos_Vec
   real(RK):: KappaL2
   real(RK):: sinfac,cosfac
   real(RK)::Faktor,Faktor2

   integer :: i,l
   integer,intent(in)::nc,np
   integer,intent(in)::ncold,npold

#if MPI_VER > 0
   integer,pointer :: i0
   integer,pointer :: i1
#endif

! Declarations
   KappaL2 = 1.0_RK/(2._RK*this%KappaL**2)

! Calculation
#if MPI_VER > 0
   i0 => this%NBox0
   i1 => this%NBox2
   DO i=i0,i1,1
#else
   DO i=1,this%BoxenAnzahlMax,1
#endif
     KVec = this%Ewald_Vec(:,i)
     mol => this%Component(nc)%Molecule
     mol2 => this%Component(ncold)%Molecule
     q => mol%SiteCharge(1:mol%NCharge)%e
     q2=> mol2%SiteCharge(1:mol2%NCharge)%e
     SSin_Vec =0._RK
     SCos_Vec =0._RK
       DO l=1,mol%NCharge
         RX = this%Component(nc)%Molecule%SiteCharge(l)%RX(np)
         RY = this%Component(nc)%Molecule%SiteCharge(l)%RY(np)
         RZ = this%Component(nc)%Molecule%SiteCharge(l)%RZ(np)

         if (RX < 0) RX = RX + 1._RK
         if (RY < 0) RY = RY + 1._RK
         if (RZ < 0) RZ = RZ + 1._RK

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

         if (RX < 0) RX = RX + 1._RK
         if (RY < 0) RY = RY + 1._RK
         if (RZ < 0) RZ = RZ + 1._RK

         Faktor2 = KVec(1) * RX + KVec(2)*RY + KVec(3)*RZ

         sinfac = q2(l)*sin(Faktor2)
         cosfac = q2(l)*cos(Faktor2)

         SSin_Vec = SSin_Vec - sinfac
         SCos_Vec = SCos_Vec - cosfac
       END DO

     this%SSin(i) = this%SSin(i) + SSin_Vec
     this%SCos(i) = this%SCos(i) + SCos_Vec

   END DO ! Boxenschleife

! Energy
   this%U_fourierLocal = this%Ewald_Prefac * (this%SSin*this%SSin + this%SCos*this%SCos)

#if MPI_VER > 0
   if ( Equilibration .and. CommonEqui ) then
     call MPI_Allreduce( sum(this%U_fourierLocal), EPotLocal, 1, MPI_RK, MPI_SUM, Communicator, ierror )
     this%UFourier= EPotLocal  / NProcs
   else
     EPotLocal = sum(this%U_fourierLocal)
     this%UFourier= EPotLocal
   end if

#else
   EPotLocal = sum(this%U_fourierLocal)
   this%UFourier= EPotLocal
#endif

  END subroutine TEnsemble_EwaldFourierEnergy_CF



!==============================================================!
!  Subroutine TSimulation_Ewald_FourierTermAddDel              !
!==============================================================!
   subroutine TEnsemble_EwaldFourierAddDel(this,nc,np,m)

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
   real(RK):: SSin_Vec,SCos_Vec
   real(RK):: KappaL2
   real(RK):: sinfac,cosfac
   real(RK):: Faktor

   integer :: i,l
   integer,intent(in)::nc,np,m

#if MPI_VER > 0
   integer,pointer :: i0
   integer,pointer :: i1
#endif


! Declarations
   KappaL2 = 1.0_RK/(2._RK*this%KappaL**2)

! Calculation
#if MPI_VER > 0
   i0 => this%NBox0
   i1 => this%NBox2
   DO i=i0,i1,1
#else
   DO i=1,this%BoxenAnzahlMax,1
#endif
     KVec = this%Ewald_Vec(:,i)
     mol => this%Component(nc)%Molecule
     q => mol%SiteCharge(1:mol%NCharge)%e
     SSin_Vec =0._RK
     SCos_Vec =0._RK
       DO l=1,mol%NCharge
         RX = this%Component(nc)%Molecule%SiteCharge(l)%RX(np)
         RY = this%Component(nc)%Molecule%SiteCharge(l)%RY(np)
         RZ = this%Component(nc)%Molecule%SiteCharge(l)%RZ(np)

         if (RX < 0) RX = RX + 1._RK
         if (RY < 0) RY = RY + 1._RK
         if (RZ < 0) RZ = RZ + 1._RK

         Faktor = KVec(1) * RX + KVec(2)*RY + KVec(3)*RZ

         sinfac = q(l)*sin(Faktor)
         cosfac = q(l)*cos(Faktor)

         SSin_Vec = SSin_Vec + sinfac*m
         SCos_Vec = SCos_Vec + cosfac*m
       END DO

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
   if ( Equilibration .and. CommonEqui ) then
      call MPI_Allreduce( sum(this%U_fourierLocal), EPotLocal, 1, MPI_RK, MPI_SUM, Communicator, ierror )
      this%UFourier= EPotLocal  / NProcs
   else
      EPotLocal = sum(this%U_fourierLocal)
      this%UFourier= EPotLocal
   end if

#else
   EPotLocal = sum(this%U_fourierLocal)
   this%UFourier= EPotLocal
#endif


  END subroutine TEnsemble_EwaldFourierAddDel


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

            ddx = 0.5_RK*fac * dx/this%BoxLength
            ddy = 0.5_RK*fac * dy/this%BoxLength
            ddz = 0.5_RK*fac * dz/this%BoxLength

            call CorrectGear_Constraint(this%Component(aacomp),aa,dLogVolumeThird,Forc, ddx,ddy,ddz)
            call CorrectGear_Constraint(this%Component(bbcomp),bb,dLogVolumeThird,Forc,-ddx,-ddy,-ddz)
         end if

         this%FCons(j) = this%FCons(j) + Forc
         this%UCons(j) = this%UCons(j) + Forc*dist

       END DO

    END DO


  end subroutine TEnsemble_Constraints

#endif



!==============================================================!
!  Subroutine TSimulation_PME_FourierTerm                      !
!==============================================================!
#if SPME > 0
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

   real(RK)           :: err

   real(RK)           :: qgrid_safe


#if MPI_VER > 0
   integer:: i0,i1
   real(RK) :: Virmpi, Fcum(3)
   real(RK),pointer, contiguous :: qgrid(:,:)
   real(RK) :: mult(this%gridx*this%gridx*this%gridx+3)
   real(RK) :: mult2(this%gridx*this%gridx*this%gridx+3)
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
#if MPI_VER > 0
   qgrid => this%qgrida
   mult = 0._RK
   mult2 = 0._RK
   i0 = this%NBox0
   i1 = this%NBox2
   DO i=i0,i1,1
     if (i==1) cycle
#else
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

     EPotLocal = EPotLocal + eterm*struc

     VirialLocal = VirialLocal + eterm*struc*(3._RK  + wterm)

#if MPI_VER > 0
     mult(index_loc)      =  eterm * factor
#else
     this%qgrida(1,index_loc) = this%qgrida(1,index_loc) * eterm * factor
     this%qgrida(2,index_loc) = this%qgrida(2,index_loc) * eterm * factor
#endif
   END DO

#if MPI_VER > 0
      call MPI_Reduce( EPotLocal, this%UFourier, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( VirialLocal, this%EVirial, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Reduce( mult, mult2, ngrid+1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )

      this%UFourier = 0.5*this%UFourier * fac2
      this%EVirial  = -0.5*(this%EVirial)*fac2*Third

! Update of the charge vector on the grid, normally done in the loop above.
      this%qgrida(1,1:ngrid+1) = this%qgrida(1,1:ngrid+1)*mult2(1:ngrid+1)
      this%qgrida(2,1:ngrid+1) = this%qgrida(2,1:ngrid+1)*mult2(1:ngrid+1)

#else
   this%UFourier=  0.5*EPotLocal * fac2
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
   FXcum = 0._RK
   FYcum = 0._RK
   FZcum = 0._RK
   DO i=1,this%NComponents
     pm=>this%Component(i)%Molecule
     DO j=1,pm%NCharge
       q = pm%SiteCharge(j)%e
#if MPI_VER > 0
       DO k=this%Component(i)%NPart0,this%Component(i)%NPart2
#else
       DO k=1,this%Component(i)%NPart ! can and should! be parallised using this%Component(i)%NPart0/2 (Michael Sch.)
#endif

       strucx = 0._RK
       strucy = 0._RK
       strucz = 0._RK

!!!!!!!!!!!!!!!!!!!!!
       q  = pm%SiteCharge(j)%e
       RX = pm%SiteCharge(j)%RX(k)
       RY = pm%SiteCharge(j)%RY(k)
       RZ = pm%SiteCharge(j)%RZ(k)

       RX     = RX - anint(RX)
       RY     = RY - anint(RY)
       RZ     = RZ - anint(RZ)

       RXgit1 = this%gridx*(RX+0.5_RK)
       RYgit1 = this%gridy*(RY+0.5_RK)
       RZgit1 = this%gridz*(RZ+0.5_RK)

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
               strucx = strucx - dsplinex(xi)*spliney(yi)*splinez(zi)*qgrid_safe
               strucy = strucy - splinex(xi)*dspliney(yi)*splinez(zi)*qgrid_safe
               strucz = strucz - splinex(xi)*spliney(yi)*dsplinez(zi)*qgrid_safe
           END DO
         END DO
       END DO
!!!!!!!!!1



         FX => pm%SiteCharge(j)%FX(k)
         FY => pm%SiteCharge(j)%FY(k)
         FZ => pm%SiteCharge(j)%FZ(k)

         FXi = strucx*facx*q
         FYi = strucy*facy*q
         FZi = strucz*facz*q

         FX    = FX + FXi
         FY    = FY + FYi
         FZ    = FZ + FZi

         FXcum = FXcum + FXi
         FYcum = FYcum + FYi
         FZcum = FZcum + FZi

         VirLoc = VirLoc + FXi*(pm%SiteCharge(j)%RX(k) - pm%SiteCharge(j)%PX(k))*this%BoxLength
         VirLoc = VirLoc + FYi*(pm%SiteCharge(j)%RY(k) - pm%SiteCharge(j)%PY(k))*this%BoxLength
         VirLoc = VirLoc + FZi*(pm%SiteCharge(j)%RZ(k) - pm%SiteCharge(j)%PZ(k))*this%BoxLength

       END DO
     END DO
   END DO

#if MPI_VER > 0
      Virmpi = VirLoc
      call MPI_Reduce( Virmpi, VirLoc, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
      call MPI_Allreduce( FXcum, Fcum(1), 1, MPI_RK, MPI_SUM, Communicator, ierror )
      call MPI_Allreduce( FYcum, Fcum(2), 1, MPI_RK, MPI_SUM, Communicator, ierror )
      call MPI_Allreduce( FZcum, Fcum(3), 1, MPI_RK, MPI_SUM, Communicator, ierror )
      FXcum = Fcum(1) /this%NPart
      FYcum = Fcum(2) /this%NPart
      FZcum = Fcum(3) /this%NPart
#else
      FXcum = FXcum /this%NPart
      FYcum = FYcum /this%NPart
      FZcum = FZcum /this%NPart
#endif

   DO i=1,this%NComponents
     DO j=1,this%Component(i)%Molecule%NCharge
       pm=>this%Component(i)%Molecule
#if MPI_VER > 0
       DO k=this%Component(i)%NPart0,this%Component(i)%NPart2 ! can also be rewritten to sequential...makes only sense for large systems
#else
       DO k=1,this%Component(i)%NPart
#endif
         pm%SiteCharge(j)%FX(k) = pm%SiteCharge(j)%FX(k) - FXcum
         pm%SiteCharge(j)%FY(k) = pm%SiteCharge(j)%FY(k) - FYcum
         pm%SiteCharge(j)%FZ(k) = pm%SiteCharge(j)%FZ(k) - FZcum
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
           spline(i-j) = div * ((w+dj)*spline(i-j-1) + (di-dj-w)*spline(i-j))
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
      spline(order- i) = div*((w+di)*spline(order-i-1) + (dorder-di-w)*spline(order-i))
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


   this%qgrida = 0._RK

   NX = this%gridx
   NY = this%gridy
   NZ = this%gridz

   boxl  = this%BoxLength
   order = this%splineorder


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

       RXgit1 = this%gridx*(RX+0.5_RK)
       RYgit1 = this%gridy*(RY+0.5_RK)
       RZgit1 = this%gridz*(RZ+0.5_RK)

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
           spline(i-j) = div * ((w+dj)*spline(i-j-1) + (di-dj-w)*spline(i-j))
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
      spline(order- i) = div*((w+di)*spline(order-i-1) + (dorder-di-w)*spline(order-i))
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

    else if (SimulationType .eq. MonteCarlo .or. SimulationType .eq. Gibbs ) then
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
    this%NBox1 = ProcRange( ngrid, this%NBox0, this%NBox2 )
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
    real(RK),pointer, contiguous :: barr(:)
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
           spline(i-j) = div * ((w+dj)*spline(i-j-1) + (di-dj-w)*spline(i-j))
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
      spline(order- i) = div*((w+di)*spline(order-i-1) + (dorder-di-w)*spline(order-i))
    END DO
    spline(1) = spline(1) * div* (1._RK-w)

    end function fillspline

  end subroutine TEnsemble_PME_Setup





  subroutine TEnsemble_PMESelfTerm_MC( this )

    implicit none

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
    this%USelbstTerm = 0.0_RK
    DO i=1,this%NComponents,1
       USelbstTermKomp = 0.0_RK
       DO Si=1,this%Component(i)%Molecule%NCharge,1
         USelbstTermKomp = USelbstTermKomp + this%Component(i)%Molecule%SiteCharge(Si)%e**2
       END DO
       this%USelbstTerm = this%USelbstTerm + this%Component(i)%NPart * USelbstTermKomp
    END DO
    this%USelbstTerm = -this%USelbstTerm * this%Kappa / sqrt(Pi)


! Intramolecular
    this%UIntra = 0._RK
    DO i=1,this%NComponents,1
      UIntraTermKomp = 0.0_RK
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

          if (dr .ge. 0.0000001) then
            call ErrorApprox(this%interaction(i,i)%PotChargeCharge(Si,Sj), this%Kappa*dr, approx)

            UIntraTermKomp = UIntraTermKomp - this%Component(i)%Molecule%SiteCharge(Si)%e* &
&                            this%Component(i)%Molecule%SiteCharge(Sj)%e / dr * (1-approx)
          end if
        END DO
      END DO
      this%UIntra = this%UIntra + this%component(i)%NPart * UIntraTermKomp
    END DO


    end subroutine TEnsemble_PMESelfTerm_MC



!==============================================================!
!  Subroutine TEnsemble_PME_FourierTerm MonteCarlo             !
!==============================================================!

   subroutine TEnsemble_PMEFourierTerm_MC(this)

   implicit none

    include 'fftw3.f'

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
   real(RK):: eterm,wterm
   real(RK):: struc
   real(RK):: energ
!Pointers
   real(RK),pointer, contiguous::lad(:,:)
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

   EPotLocal = 0._RK
   VirialLocal = 0._RK

! Summation over all the Energies
#if MPI_VER > 0
   i0 = this%NBox0
   i1 = this%NBox2
   DO i=i0,i1,1
#else
   DO i=2,ngrid
#endif
! Positioning
     k1    = int((i-1) / ngridyz)
     k1bck = int((i-1) - k1*ngridyz)
     k2    = int(k1bck / NZ)
     k3    = int(k1bck - k2*NZ)

     index_loc = k3 + k2*NZ + k1*ngridyz + 1

! Charge contribution
     qr = lad(1,index_loc)
     qi = lad(2,index_loc)
     struc = qr*qr + qi*qi

! Distance contribution
     mult  = fac*this%mm2(i)
     eterm = exp(-mult) / this%bbtot(i)

! Energy
     energ = eterm*struc
     EPotLocal = EPotLocal + energ
     VirialLocal = VirialLocal + energ*(1._RK - mult)
   END DO

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

   this%qgrida = 0._RK

   NX = this%gridx
   NY = this%gridy
   NZ = this%gridz

   boxl  = this%BoxLength
   order = this%splineorder


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
           spline(i-j) = div * ((w+dj)*spline(i-j-1) + (di-dj-w)*spline(i-j))
         END DO
       spline(1) = spline(1) * div * (1._RK-w)
       END DO

! Generate order spline
    div = 1._RK / (order - 1._RK)
    spline(order) = div * w * spline(order-1)
    dorder = real(order)
    DO i=1,order-2,1
      di = real(i)
      spline(order- i) = div*((w+di)*spline(order-i-1) + (dorder-di-w)*spline(order-i))
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

   NX = this%gridx
   NY = this%gridy
   NZ = this%gridz

   boxl  = this%BoxLength
   order = this%splineorder

   fac = UnitCharge / ElementaryCharge

   DO j=1,this%Component(nc)%Molecule%NCharge

       q  = this%Component(nc)%Molecule%SiteCharge(j)%e
       RX = this%Component(nc)%Molecule%SiteCharge(j)%RX(np)
       RY = this%Component(nc)%Molecule%SiteCharge(j)%RY(np)
       RZ = this%Component(nc)%Molecule%SiteCharge(j)%RZ(np)

       RX     = RX - anint(RX)
       RY     = RY - anint(RY)
       RZ     = RZ - anint(RZ)

       RXgit1 = this%gridx*(RX+0.5_RK)
       RYgit1 = this%gridy*(RY+0.5_RK)
       RZgit1 = this%gridz*(RZ+0.5_RK)

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
           spline(i-j) = div * ((w+dj)*spline(i-j-1) + (di-dj-w)*spline(i-j))
         END DO
       spline(1) = spline(1) * div * (1._RK-w)
       END DO

! Generate order spline
    div = 1._RK / (order - 1._RK)
    spline(order) = div * w * spline(order-1)
    dorder = real(order)
    DO i=1,order-2,1
      di = real(i)
      spline(order- i) = div*((w+di)*spline(order-i-1) + (dorder-di-w)*spline(order-i))
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

       RXgit1 = this%gridx*(RX+0.5_RK)
       RYgit1 = this%gridy*(RY+0.5_RK)
       RZgit1 = this%gridz*(RZ+0.5_RK)

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
           spline(i-j) = div * ((w+dj)*spline(i-j-1) + (di-dj-w)*spline(i-j))
         END DO
       spline(1) = spline(1) * div * (1._RK-w)
       END DO

! Generate order spline
    div = 1._RK / (order - 1._RK)
    spline(order) = div * w * spline(order-1)
    dorder = real(order)
    DO i=1,order-2,1
      di = real(i)
      spline(order- i) = div*((w+di)*spline(order-i-1) + (dorder-di-w)*spline(order-i))
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

   Virloc = 0._RK
   Kappa  = this%Kappa
   Faktor = 2._RK/sqrt(Pi) * Kappa


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
!          !TODO: Check if this should be changed to ErrorApprox
!          call erfc_approx(KappaRij,approx)
         !Michael: For SPME correct?
         if (dr .ge. 0.0000001) then
           call ErrorApprox(this%Interaction(i,i)%PotChargeCharge(j,jj), KappaRij, approx)
           Fij  = (qj*qjj/dr*(1._RK-approx) - Faktor*exp(-KappaRij**2)*qj*qjj) /dr
         else
           Fij  = (qj*qjj/dr - Faktor*exp(-KappaRij**2)*qj*qjj) /dr
         end if

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
!          !TODO: Check if this should be changed to ErrorApprox
!          call erfc_approx(KappaRij,approx)
         !Michael: For SPME correct?
         if (dr .ge. 0.0000001) then
           call ErrorApprox(this%Interaction(i,i)%PotChargeCharge(j,jj), KappaRij, approx)
           Fij  = (qj*qjj/dr*(1._RK-approx) - Faktor*exp(-KappaRij**2)*qj*qjj) /dr
         else
           Fij  = (qj*qjj/dr - Faktor*exp(-KappaRij**2)*qj*qjj) /dr
         end if

         eX = dx / dr
         eY = dy / dr
         eZ = dz / dr

         Fij  = (qj*qjj/dr*(1._RK-approx) - Faktor*exp(-KappaRij**2)*qj*qjj) /dr

         Virloc = Virloc + Fij* (eX * dx + eY * dy + eZ * dz)
        END DO
       END DO
     END DO
   END DO
!!! Michael Sch.: This routine does nothing...no result transferrred no type/class variable safed/changed!

   end subroutine TEnsemble_VirialIntra

#endif


! TRANSPORT ab hier
#if TRANS==1
!==============================================================!
!  Subroutine TEnsemble_CalCorrFun                             !
!==============================================================!

    subroutine TEnsemble_CalCorrFun( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this
    ! Declare local variables
    integer  :: nmess, i, j, j0, j1, j2, k, l, s, CFTMP
    integer  :: CFindex, Mindex
    integer  :: NPart, NPart2, StepCorr
    integer  :: np, nc, np1, np2, ncomp2
    real(RK) :: qi, qj
    real(RK) :: sx(this%NComponents), sy(this%NComponents)
    real(RK) :: sz(this%NComponents)
    real(RK) :: SXindex(this%NComponents),SYindex(this%NComponents),SZindex(this%NComponents)
    real(RK) :: Sindex(3), ss(3)
    real(RK) :: KinERot(this%NPart)
    real(RK) :: BoxLength_dt,BoxLength_dt2
    real(RK) :: tempf(3), virf(3)
    real(RK) :: Mass
    real(RK), pointer, contiguous :: pFB(:,:), pFS(:,:), pFTC(:,:), pFRC(:,:)
    type(TComponent),pointer :: pc
    logical  :: Conductivity, EConductivity, Bulkviscosity

    NPart  = this%NPart
    NPart2 = 2*this%NPart
    BoxLength_dt       =  this%BoxLength/TimeStep
    BoxLength_dt2      =  BoxLength_dt**2
    Conductivity = this%Conductivity
    EConductivity = this%EConductivity
    Bulkviscosity = this%Bulkviscosity
    ncomp2 = this%NComponents**2

    !Reduced correlation steps
    StepCorr = (Step + this%NStepCorr -1) / this%NStepCorr
 

    !Calculate matrix indexes
    Mindex = mod(StepCorr, this%NCorr )
    if (Mindex .eq. 0) then
      Mindex = this%NCorr
    end if

    !Write transport properties Matrixes (root Processor)
    this%vsk(Mindex,  :) = 0._RK
    this%vsp(Mindex,  :) = 0._RK
    this%vbk(Mindex,  :) = 0._RK
    this%vbp(Mindex,  :) = 0._RK
    this%vckt(Mindex, :) = 0._RK
    this%vcpt(Mindex, :) = 0._RK
    this%vckr(Mindex, :) = 0._RK
    this%vcpr(Mindex, :) = 0._RK
    this%vcmt(Mindex, :) = 0._RK

    !Evaluate FTC and FRC components (parallel version)
    do i = 1, this%NComponents
      call ForceTransport( this%Component(i) )
    end do


  if (RootProc) then
    ! Loop Variable
    j0 = 0
    do i = 1, this%NComponents
      pc => this%Component(i)
      np   = pc%NPart
      Mass =  pc%Molecule%Mass
#if MPI_VER > 0
      pFB => this%Component(i)%FBAll(:,:)
      pFS => this%Component(i)%FSAll(:,:)
#else
      pFB => this%Component(i)%FB(:,:)
      pFS => this%Component(i)%FS(:,:)
#endif

   
        pFTC => this%Component(i)%FTC(:,:)
        pFRC => this%Component(i)%FRC(:,:)

      !Michael Sch.: works only for 1unit per molecule, see also interaction.F90
        do j = 1, np
          if ( pc%Molecule%IsElongated ) then
            KinERot(j)= sum( pc%W0(j,1:3, 1) * pc%W0(j,1:3, 1) * pc%Molecule%Unit(1)%MOI(1:3))*0.5_RK
          end if
        end do

      do k =1, 3
        ! Calculate sum of terms of the pressure tensor (kinetic and potential)
        this%vsp(Mindex, k)  = this%vsp(Mindex, k) + sum(pFS (:, k)) ! potential part off-diagonal elements of pressure tensor
        this%vbp(Mindex, k)  = this%vbp(Mindex, k) + sum(pFB (:, k)) ! potential part diagonal elements of the pressure tensor
        this%vbk(Mindex, k)  = this%vbk(Mindex, k) + pc%KinETranTotal(k)!kinetic part diagonal elements of the pressure tensor
        
        if (Bulkviscosity) then
          this%sc(k) = this%sc(k) + pc%KinETranTotal(k)
          this%sp(k) = this%sp(k) + sum(pFB(:, k))
        end if
          
                  
        if (Conductivity) then
          this%vcpr(Mindex, k) = this%vcpr(Mindex, k) + sum(pFRC(:, k))  !Thermal conductivity for mixtures
          this%vcpt(Mindex, k) = this%vcpt(Mindex, k) + sum(pFTC(:, k))
          this%vckt(Mindex, k) = this%vckt(Mindex, k) + sum( pc%P1(:, k, 1)*sum( pc%KinETran(:,1:3),2 ) )*0.5_RK*BoxLength_dt
          this%vcmt(Mindex, k) = this%vcmt(Mindex, k) + pc%PartialMolarEnthalpy*sum(pc%P1(:, k, 1))*BoxLength_dt 
          if ( pc%Molecule%IsElongated ) then
            this%vckr(Mindex, k)= this%vckr(Mindex, k) + sum( pc%P1(:, k, 1) * KinERot(:) ) * BoxLength_dt
          end if
        end if
      end do !k =1, 3

      ! kinetic part
      !Diffusion matrix a
      k=mod(Mindex,this%NSpanCF)
      if (k==0) k=this%NSpanCF
      ! Multiplikation mit BoxLength_dt erst nach der Integration
      this%A_SpanCF(j0+1:j0+np              , k) = pc%P1(1:np,1,1)   ! P1 only works for 1unit per molecule, not defined for COM
      this%A_SpanCF(j0+NPart+1 :j0+NPart +np, k) = pc%P1(1:np,2,1)
      this%A_SpanCF(j0+NPart2+1:j0+NPart2+np, k) = pc%P1(1:np,3,1)

      !parts of the stress and energy tensors
      ! shear off diagonal terms
      this%vsk(Mindex, 1) = this%vsk(Mindex,1) + sum( pc%P1(:,1,1) * pc%P1(:,2,1) ) * Mass*BoxLength_dt2  ! P1 only works for 1unit per molecule, not defined for COM
      this%vsk(Mindex, 2) = this%vsk(Mindex,2) + sum( pc%P1(:,1,1) * pc%P1(:,3,1) ) * Mass*BoxLength_dt2
      this%vsk(Mindex, 3) = this%vsk(Mindex,3) + sum( pc%P1(:,2,1) * pc%P1(:,3,1) ) * Mass*BoxLength_dt2

      j0 = j0 + np
    end do    ! Component

    if (Bulkviscosity) then
      tempf(:)  = this%sc(:)/StepCorr
      virf(:)   = this%sp(:)/StepCorr
    end if

   
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate Auto Correlation Functions      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if ( mod(StepCorr, this%NSpanCF) .eq. 0 ) then
     if (StepCorr .gt. this%NCorr) then

      CFindex = Mindex +1
      this%a(:,CFindex - this%NSpanCF:CFindex-1) = this%A_SpanCF(:,1:this%NSpanCF)

      if (Mindex .eq. this%NCorr) CFindex = 1      !index of t = t0

     ! nullify the autocorrelation functions
	do i = 1, this%NComponents
	  this%cf_d(i, :) = 0._RK
        end do

        this%cf_vs(:) = 0._RK
        this%cf_c(:)  = 0._RK
        this%cf_vb(:) = 0._RK
        this%cf_ec(:) = 0._RK

        if (this%NComponents == 2) this%cf_soret = 0._RK
 
        if (this%NComponents .gt. 1) then 
          do k = 1, ncomp2
             this%lamda(k,:) = 0._RK
          end do
        end if


       ! Preparation of the Autocorrelation function - safe the Startpoints
        j0 = 0
        do i = 1, this%NComponents
          np = this%Component(i)%NPart
          SXindex(i)  = sum(this%a(j0       +1:j0+np        , CFindex))
          SYindex(i)  = sum(this%a(j0+NPart +1:j0+NPart +np , CFindex))
          SZindex(i)  = sum(this%a(j0+NPart2+1:j0+NPart2+np , CFindex))
          j0 = j0 + np
        end do

        Sindex(1)=SXindex(2)*BoxLength_dt  !mass flux for thermal diffusion
        Sindex(2)=SYindex(2)*BoxLength_dt 
        Sindex(3)=SZindex(2)*BoxLength_dt 

      ! Calculation of all transport properties 
      ! s .. matrix index of the corresponding values
!mn      s = CFindex    Has to be set inside of the loop for OMP parallelisation
      CFtmp=CFindex
! Directive inserted by Cray Reveal.  May be incomplete.
!$OMP  parallel do default(none)                                         &
!$OMP&   private (i,j,j0,j1,j2,k,l,nc,nmess,np,np1,np2,qi,qj,s)          &
!$OMP&   shared  (this,cfindex,econductivity,npart,npart2,sxindex,       &
!$OMP&            syindex,szindex,tempf,unitcharge,virf,CFtmp           &
!$OMP&            BoxLength_dt,Sindex,Bulkviscosity)           &
!$OMP&   private (sx,sy,sz,ss)

      
      do nmess= 1, this%NCorr
         ! Loop over particles      
        s=CFtmp+nmess-1
        if (s > this%NCorr) s = s-this%NCorr
        j0 = 0
        do i = 1, this%NComponents
          np = this%Component(i)%NPart
          this%cf_d(i, nmess) = this%cf_d(i, nmess) + DOT_PRODUCT( this%a(j0+1 : j0+np,CFindex) , &
&                                                     this%a(j0+1: j0+np,s) ) &
&                             + DOT_PRODUCT( this%a(j0+NPart +1 : j0+NPart +np,CFindex) , &
&                                                     this%a(j0+NPart +1 : j0+NPart +np,s) ) &
&                             + DOT_PRODUCT( this%a(j0+NPart2+1 : j0+NPart2+np,CFindex) , &
&                                                     this%a(j0+NPart2+1 : j0+NPart2+np,s) )

          if ( this%NComponents .gt. 1 ) then
            sx(i)  = sum(this%a(j0       +1:j0+np        , s))
            sy(i)  = sum(this%a(j0+NPart +1:j0+NPart +np , s))
            sz(i)  = sum(this%a(j0+NPart2+1:j0+NPart2+np , s))
          end if
          j0 = j0 + np
        end do

        ss(1) = sx(2)* BoxLength_dt 
        ss(2) = sy(2)* BoxLength_dt 
        ss(3) = sz(2)* BoxLength_dt  

        ! Just loops over components!
        if (this%NComponents .gt. 1) then
          nc = this%NComponents
          k = 1
          do i = 1, nc
            do j = 1,nc
              this%lamda(k, nmess) = this%lamda(k, nmess) + SXindex(i)*sx(j) &
&                                                         + SYindex(i)*sy(j) &
&                                                         + SZindex(i)*sz(j)
              k = k + 1
            end do
          end do
        end if

        ! Calculated in general
        do k = 1, 3
          ! shear viscosity (off-diagonal elements)
          this%cf_vs(nmess) = this%cf_vs(nmess) + this%vsk(CFindex, k)*this%vsk(s, k) + &
&                                                 this%vsp(CFindex, k)*this%vsp(s, k) + &
&                                                 this%vsk(CFindex, k)*this%vsp(s, k) + &
&                                                 this%vsp(CFindex, k)*this%vsk(s, k)

          ! bulk viscosity
          if (Bulkviscosity) then
            do j = 1, 3 
              this%cf_vb(nmess) =   this%cf_vb(nmess) + &
&                               ( this%vbk(CFindex, j)-tempf(j))*(this%vbk(s, k)-tempf(k)) + &
&                               ( this%vbp(CFindex, j)-virf(j)) *(this%vbp(s, k)-virf(k) ) + &
&                               ( this%vbk(CFindex, j)-tempf(j))*(this%vbp(s, k)-virf(k) ) + &
&                               ( this%vbp(CFindex, j)-virf(j))*(this%vbk(s, k)-tempf(k) )
            end do
          end if

          ! conductivity
            this%cf_c(nmess) =  this%cf_c(nmess) + this%vckt(CFindex, k)*this%vckt(s, k) + &
&                                                  this%vckt(CFindex, k)*this%vcpt(s, k) + &
&                                                  this%vckt(CFindex, k)*this%vckr(s, k) + &
&                                                  this%vckt(CFindex, k)*this%vcpr(s, k) + &
&                                                  this%vckt(CFindex, k)*this%vcmt(s, k) + &
&                                                  this%vckr(CFindex, k)*this%vckt(s, k) + &
&                                                  this%vckr(CFindex, k)*this%vcpt(s, k) + &
&                                                  this%vckr(CFindex, k)*this%vckr(s, k) + &
&                                                  this%vckr(CFindex, k)*this%vcpr(s, k) + &
&                                                  this%vckr(CFindex, k)*this%vcmt(s, k) + &
&                                                  this%vcpt(CFindex, k)*this%vckt(s, k) + &
&                                                  this%vcpt(CFindex, k)*this%vcpt(s, k) + &
&                                                  this%vcpt(CFindex, k)*this%vckr(s, k) + &
&                                                  this%vcpt(CFindex, k)*this%vcpr(s, k) + &
&                                                  this%vcpt(CFindex, k)*this%vcmt(s, k) + &
&                                                  this%vcpr(CFindex, k)*this%vckt(s, k) + &
&                                                  this%vcpr(CFindex, k)*this%vcpt(s, k) + &
&                                                  this%vcpr(CFindex, k)*this%vckr(s, k) + &
&                                                  this%vcpr(CFindex, k)*this%vcpr(s, k) + &
&                                                  this%vcpr(CFindex, k)*this%vcmt(s, k) + &
&                                                  this%vcmt(CFindex, k)*this%vckt(s, k) + &
&                                                  this%vcmt(CFindex, k)*this%vcpt(s, k) + &
&                                                  this%vcmt(CFindex, k)*this%vckr(s, k) + &
&                                                  this%vcmt(CFindex, k)*this%vcpr(s, k) + &
&                                                  this%vcmt(CFindex, k)*this%vcmt(s, k) 

        end do

        ! include the digonal elements to the shear viscosity (Pxx-Pyy)/2 and (Pyy-Pzz)/2
          this%cf_vs(nmess) = this%cf_vs(nmess) + (1._RK/4._RK)* ((this%vbk(CFindex, 1)-this%vbk(CFindex, 2))*(this%vbk(s, 1)-this%vbk(s, 2)) + &
&                                                              (this%vbp(CFindex, 1)-this%vbp(CFindex, 2))*(this%vbp(s, 1)-this%vbp(s, 2)) + &
&                                                              (this%vbk(CFindex, 1)-this%vbk(CFindex, 2))*(this%vbp(s, 1)-this%vbp(s, 2)) + &
&                                                              (this%vbp(CFindex, 1)-this%vbp(CFindex, 2))*(this%vbk(s, 1)-this%vbk(s, 2)) + &
&                                                              (this%vbk(CFindex, 2)-this%vbk(CFindex, 3))*(this%vbk(s, 2)-this%vbk(s, 3)) + &
&                                                              (this%vbp(CFindex, 2)-this%vbp(CFindex, 3))*(this%vbp(s, 2)-this%vbp(s, 3)) + &
&                                                              (this%vbk(CFindex, 2)-this%vbk(CFindex, 3))*(this%vbp(s, 2)-this%vbp(s, 3)) + &
&                                                              (this%vbp(CFindex, 2)-this%vbp(CFindex, 3))*(this%vbk(s, 2)-this%vbk(s, 3)))


         !Thermal diffusivity
	if (this%Ncomponents==2) then
           do k = 1, 3
	     this%cf_soret(nmess) =  this%cf_soret(nmess) + this%vckt(CFindex, k)*ss(k) + &
&	               					    this%vckr(CFindex, k)*ss(k) + &
&		  					    this%vcpt(CFindex, k)*ss(k) + &
&							    this%vcpr(CFindex, k)*ss(k) + &
&							    this%vcmt(CFindex, k)*ss(k) + &
&							    Sindex(k)*this%vckt(s, k) + &
&							    Sindex(k)*this%vckr(s, k) + &
&							    Sindex(k)*this%vcpt(s, k) + &
&							    Sindex(k)*this%vcpr(s, k) + &
&							    Sindex(k)*this%vcmt(s, k) 
            end do
         end if

        ! electric conductivity
        if (EConductivity) then
           j1 = 0
           do k = 1, this%NComponents
              np1 = this%Component(k)%NPart
              if((this%Component(k)%Molecule%Charge) .ne. 0._RK) then !Electric conductivity only defined for charged particles
!             if((this%Component(k)%Molecule%Charge) .ne. 0._RK) then !Electric conductivity only defined for charged particles
                 qi = this%Component(k)%Molecule%Charge*UnitCharge/ElementaryCharge
                 do i = 1, np1
                    j2 = 0
                    do l = 1, this%NComponents
                       np2 = this%Component(l)%NPart
                       if( (this%Component(l)%Molecule%Charge) .ne. 0._RK) then !Electric conductivity only defined for charged particles
                           qj = this%Component(l)%Molecule%Charge*UnitCharge/ElementaryCharge
                           do j = 1, np2
                              this%cf_ec(nmess) = this%cf_ec(nmess) + qi * qj * this%a(j1+i, CFindex) * this%a(j2+j, s)
                              this%cf_ec(nmess) = this%cf_ec(nmess) + qi * qj * this%a(j1+NPart+i, CFindex) * this%a(j2+NPart+j, s)
                              this%cf_ec(nmess) = this%cf_ec(nmess) + qi * qj * this%a(j1+NPart2+i, CFindex) * this%a(j2+NPart2+j, s)          
                           end do
                       endif
                       j2 = j2 + np2
                    end do
                 end do  
              endif 
              j1 = j1 + np1    
           end do
        end if !if(EConductivity)

!mn        if (s == this%NCorr) s = 0
!mn        s = s+1

      end do  ! NMess
      this%Mmess  = this%Mmess +1

      do i = 1, this%NComponents
	 this%average_cf_d(i, :) = (this%average_cf_d(i,:) + this%cf_d(i,:))
      end do

      this%average_cf_vs(:)= (this%average_cf_vs(:) + this%cf_vs(:))
      this%average_cf_vb(:)= (this%average_cf_vb(:) + this%cf_vb(:))
      this%average_cf_c(:) = (this%average_cf_c(:) + this%cf_c(:))
      this%average_cf_ec(:)= (this%average_cf_ec(:) + this%cf_ec(:))

      if (this%NComponents == 2) this%average_cf_soret(:)= (this%average_cf_soret(:) + this%cf_soret(:))

      if (this%NComponents .gt. 1) then 
        do k = 1, ncomp2
           this%average_lamda(k,:) = (this%average_lamda(k,:) + this%lamda(k,:))
        end do
      end if

      
      ! Call integration for ACF
      call IntCorrFun ( this )

! -----------------------------------
     else ! if (Step .gt. this%NCorr)
        CFindex = Mindex +1
        this%a(:,CFindex - this%NSpanCF:CFindex-1) = this%A_SpanCF(:,1:this%NSpanCF)
     end if ! if (Step .gt. this%NCorr)
    end if ! if (mod(Step, this%NSpanCF).eq.0)

 end if ! RootProc
   end subroutine TEnsemble_CalCorrFun
#endif



#if TRANS==1
!==============================================================!
!  Subroutine TEnsemble_IntCorrFun                             !
!==============================================================!

  subroutine TEnsemble_IntCorrFun( this )

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this
    ! Declare local varibles
    integer  :: i, j, k
    integer  :: ncomp2
    real(RK) :: helpvar!, det, deter1, deter2, deter3, deter4
    real(RK) :: x1, x2, x3, w1, w2, MM
!    real(RK) :: Inv_x1, Inv_x2, Inv_x3
!    real(RK) :: B11, B12, B21, B22
    real(RK) :: BoxLength_dt2

    BoxLength_dt2      =  (this%BoxLength/TimeStep)**2
    ncomp2 = this%NComponents*this%NComponents
    MM = this%Component(1)%Molecule%Mass*this%Component(1)%Fraction + this%Component(2)%Molecule%Mass*this%Component(2)%Fraction
    w1 = this%Component(1)%Molecule%Mass*this%Component(1)%Fraction/MM
    w2 = this%Component(2)%Molecule%Mass*this%Component(2)%Fraction/MM

    

    do i  = 1, this%NComponents
      helpvar =  1._RK /(3._RK *this%Component(i)%NPart) * BoxLength_dt2
      if (abs(this%cf_d(i, 1)) .gt. 1e-15) then 
         this%sinte_i(i,:) = simpson( this%cf_d(i,:)/this%cf_d(i, 1), this%TimeStepCorr, this%NCorr )
         this%average_sinte_i(i,:) = simpson( this%average_cf_d(i,:)/this%average_cf_d(i, 1),this%TimeStepCorr, this%NCorr )
         this%average_sinte_i(i,:) = this%average_sinte_i(i,:)*this%average_cf_d(i, 1)*helpvar/this%Mmess
         this%selfd_i(i) = this%sinte_i(i, this%NCorr) * this%cf_d(i, 1) * helpvar
      end if
    end do


    if ( this%NComponents .gt. 1) then
      helpvar =  1._RK /(3._RK *this%NPart) * BoxLength_dt2
      do k = 1, ncomp2
        if (abs(this%lamda(k, 1)) .gt. 1e-15) then 
           this%sinte_lamda(k, :) = simpson(this%lamda(k,:)/this%lamda(k,1),this%TimeStepCorr, this%NCorr)
           this%average_sinte_lamda(k,:) = simpson(this%average_lamda(k,:)/this%average_lamda(k,1),this%TimeStepCorr, this%NCorr)
           this%average_sinte_lamda(k,:) = this%average_sinte_lamda(k,:)* this%average_lamda(k,1)*helpvar/this%Mmess
        end if
      end do
      
      k = 1
       do i = 1, this%NComponents
         do j = 1, this%NComponents
             this%Onsager(i,j) = this%sinte_lamda(k,this%NCorr)*this%lamda(k,1)*helpvar
             k = k +1
         end do
      end do   

    end if


    helpvar =  this%Density /(5._RK *this%NPart * this%Temperature)
    this%sinte_vs = simpson( this%cf_vs(:)/this%cf_vs(1), this%TimeStepCorr, this%NCorr )
    this%average_sinte_vs = simpson( this%average_cf_vs(:)/this%average_cf_vs(1), this%TimeStepCorr, this%NCorr)
    this%average_sinte_vs = this%average_sinte_vs(:)*this%average_cf_vs(1)*helpvar/this%Mmess
    this%visco_s = this%sinte_vs( this%NCorr ) * this%cf_vs(1) * helpvar


    helpvar = 5._RK / 3._RK * helpvar
    if (this%Bulkviscosity) then
      this%sinte_vb = simpson( this%cf_vb(:)/this%cf_vb(1),this%TimeStepCorr, this%NCorr )
      this%average_sinte_vb = simpson( this%average_cf_vb(:)/this%average_cf_vb(1),this%TimeStepCorr, this%NCorr)
      this%average_sinte_vb = this%average_sinte_vb(:)*this%average_cf_vb(1)*helpvar/(3._RK*this%Mmess)
      this%visco_b = this%sinte_vb( this%NCorr ) * this%cf_vb(1) * (helpvar / 3._RK)
    end if

    if (this%Conductivity) then
      this%sinte_c = simpson( this%cf_c(:)/this%cf_c(1), this%TimeStepCorr, this%NCorr )
      this%average_sinte_c = simpson( this%average_cf_c(:)/this%average_cf_c(1),this%TimeStepCorr, this%NCorr)
      this%average_sinte_c = this%average_sinte_c(:)*this%average_cf_c(1)*(helpvar/(this%Temperature*this%Mmess))
      this%conduct = this%sinte_c( this%NCorr ) * this%cf_c(1) * (helpvar / this%Temperature)
    end if

    if ( this%NComponents == 2 .and. (abs(this%cf_soret(1)) .gt. 1e-15) ) then
      this%sinte_soret = simpson (this%cf_soret(:)/this%cf_soret(1), this%TimeStepCorr, this%NCorr )
      this%average_sinte_soret = simpson (this%average_cf_soret(:)/this%average_cf_soret(1), this%TimeStepCorr, this%NCorr)
      this%average_sinte_soret = this%average_sinte_soret(:)*this%average_cf_soret(1)*helpvar*this%Component(2)%Molecule%Mass/(2._RK*this%Mmess*this%Density*MM*this%Temperature*w1*w2)
      this%soret =  this%sinte_soret( this%NCorr)*this%cf_soret(1)*this%Component(2)%Molecule%Mass*helpvar/(this%Density*this%Temperature*w1*w2*MM*2._RK)
    end if


    if (this%EConductivity) then
      this%sinte_ec = simpson( this%cf_ec(:)/this%cf_ec(1), this%TimeStepCorr, this%NCorr )
      this%average_sinte_ec = simpson( this%average_cf_ec(:)/this%average_cf_ec(1),this%TimeStepCorr, this%NCorr)
      this%average_sinte_ec = this%average_sinte_ec(:)*this%average_cf_ec(1)*(helpvar * BoxLength_dt2/this%Mmess)
      this%econduct = this%sinte_ec( this%NCorr ) * this%cf_ec(1) * (helpvar * BoxLength_dt2)
    end if


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
        integral(i-1) = .5_RK * (integral(i) + integral(i-2))
      end do
      
      if( mod(n, 2) == 0 .and. n > 2 ) integral(n) = integral(n-1) + .5_RK * values(n) + 2._RK * values(n-1) + .5_RK * values(n-2)
      integral = integral * step / 3._RK
      

    end function

  end subroutine TEnsemble_IntCorrFun
#endif


#if HBOND > 0
!==============================================================!
!  Subroutine TEnsemble_HBonding                               !
!==============================================================!

  subroutine TEnsemble_HBonding( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif
    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    type(TComponent), pointer :: pacc, pdon
    type(TSiteCharge), pointer :: paccacc, pdonacc, pmixdon
    logical             :: MixTerm
    integer             :: h, i, i0, i1, j, k ,l, m
    real(RK)            :: BoxLengthInv
    real(RK)            :: LAA, LAD, LintraAD, CosAngle
    real(RK)            :: AngleCrit, DistCrit1, DistCrit2
    real(RK)            :: PXij, PYij, PZij
    real(RK)            :: drAA(3), drAD(3), drintraAD(3)
    integer,allocatable :: Counter(:,:), NHBAll(:)

    ! Initialize arrays
    this%NHBond0(:)       = 0
    this%NHBond1(:,:)     = 0
    this%NHBond2(:,:,:)   = 0
    this%NHBond3(:,:,:,:) = 0
    this%NHBondN(:)       = 0

    allocate( Counter(this%NComponents,this%NPart) )
    Counter(:,:) = 0
    allocate( NHBAll(this%NComponents) )
    NHBAll(:) = 0
    BoxLengthInv = 1._RK / this%BoxLength

    do h = 1, this%NHBondCrit

      !HBonding criteria
      DistCrit1 = this%DistCrit1(h)*BoxLengthInv
      DistCrit2 = this%DistCrit2(h)*BoxLengthInv
      AngleCrit = this%AngleCrit(h)

      !Definition of the H-Bonding components and sites
      pacc => this%Component(this%AccComp(h))
      pdon => this%Component(this%DonComp(h))

      paccacc => this%Component(this%AccComp(h))%Molecule%SiteCharge(this%AccAccSite(h))
      pdonacc => this%Component(this%DonComp(h))%Molecule%SiteCharge(this%DonAccSite(h))
      if ( this%AccDonSite(h) > 0 ) then
        pmixdon => this%Component(this%AccComp(h))%Molecule%SiteCharge(this%AccDonSite(h))
        MixTerm =.true.
      else
        pmixdon => this%Component(this%DonComp(h))%Molecule%SiteCharge(this%DonDonSite(h))
        MixTerm = .false.
      end if
 

      i0 = 1
      i1 = pacc%NPart
#if MPI_VER > 0
      if (SimulationType .eq. MolecularDynamics) then
        i0 = pacc%NPart0
        i1 = pacc%NPart2
      end if
#endif
      !Loop over all particles
      do i = i0, i1

        do j = 1,pdon%NPart

          !Minimum Image Convention ! ...old for rigid/COM based molecules
          !PXij = pacc%P0(i, 1) - pdon%P0(j, 1)
          !PYij = pacc%P0(i, 2) - pdon%P0(j, 2)
          !PZij = pacc%P0(i, 3) - pdon%P0(j, 3)

          !Calculation of the distance vector Acc of AccComp and Acc of DonComp 
          drAA(1)=(paccacc%RX(i)-pdonacc%RX(j)) !- anint( PXij )
          drAA(2)=(paccacc%RY(i)-pdonacc%RY(j)) !- anint( PYij )
          drAA(3)=(paccacc%RZ(i)-pdonacc%RZ(j)) !- anint( PZij )
          drAA(1)=drAA(1) - anint(drAA(1))
          drAA(2)=drAA(2) - anint(drAA(2))
          drAA(3)=drAA(3) - anint(drAA(3))
          LAA = SQRT( DOT_PRODUCT(drAA,drAA) )
  
          if (LAA .le. DistCrit1)then
            !Calculation of the distance vector Don of one Comp and  Acc of the other Comp
            if (MixTerm) then
              drAD(1)=(pmixdon%RX(i)-pdonacc%RX(j)) !- anint( PXij )
              drAD(2)=(pmixdon%RY(i)-pdonacc%RY(j)) !- anint( PYij )
              drAD(3)=(pmixdon%RZ(i)-pdonacc%RZ(j)) !- anint( PZij )
            else
              drAD(1)=(paccacc%RX(i)-pmixdon%RX(j)) !- anint( PXij )
              drAD(2)=(paccacc%RY(i)-pmixdon%RY(j)) !- anint( PYij )
              drAD(3)=(paccacc%RZ(i)-pmixdon%RZ(j)) !- anint( PZij )
            end if
            drAD(1)=drAD(1) - anint(drAD(1))
            drAD(2)=drAD(2) - anint(drAD(2))
            drAD(3)=drAD(3) - anint(drAD(3))
            !Hier die Minimum Image-Convention
            LAD = SQRT( DOT_PRODUCT(drAD,drAD) )

            if (LAD .le. DistCrit2) then
              !Calculation of the angle between dono and acceptors
              if (MixTerm) then
                drintraAD(1)=pmixdon%RX(i)-paccacc%RX(i)
                drintraAD(2)=pmixdon%RY(i)-paccacc%RY(i)
                drintraAD(3)=pmixdon%RZ(i)-paccacc%RZ(i)
              else
                drintraAD(1)=pmixdon%RX(j)-pdonacc%RX(j)
                drintraAD(2)=pmixdon%RY(j)-pdonacc%RY(j)
                drintraAD(3)=pmixdon%RZ(j)-pdonacc%RZ(j)
              end if
              drintraAD(1)=drintraAD(1) - anint(drintraAD(1))
              drintraAD(2)=drintraAD(2) - anint(drintraAD(2))
              drintraAD(3)=drintraAD(3) - anint(drintraAD(3))
              LintraAD = SQRT( DOT_PRODUCT(drintraAD,drintraAD) )
              CosAngle=abs( DOT_PRODUCT(drAA,drintraAD) / LAA / LintraAD )

              if (CosAngle .ge. AngleCrit) then
                ! not working for more than 99 species
                if ( Counter(this%AccComp(h),i) == 0 ) then
                  Counter(this%AccComp(h),i) = this%DonComp(h)
                elseif ( Counter(this%AccComp(h),i) < 100 ) then
                  Counter(this%AccComp(h),i) = Counter(this%AccComp(h),i) + this%DonComp(h)*100
                elseif ( Counter(this%AccComp(h),i) < 10000 ) then
                  Counter(this%AccComp(h),i) = Counter(this%AccComp(h),i) + this%DonComp(h)*10000
                else
                  Counter(this%AccComp(h),i) = 1000000
                end if
              end if

            end if

          end if

        end do !do j=1,npDon
      end do !do i=1,npAcc
    end do ! NHBondCrit

    !HBonding statistics
    do h = 1, this%NComponents
      i0 = 1
      i1 = this%Component(h)%NPart
#if MPI_VER > 0
      if (SimulationType .eq. MolecularDynamics) then
        i0 = this%Component(h)%NPart0
        i1 = this%Component(h)%NPart2
      end if
#endif
      do i = i0, i1
        m = 0
        if ( Counter(h,i) == 0 ) then
          this%NHBond0(h)=this%NHBond0(h) + 1
        elseif ( Counter(h,i) == 1000000 ) then
          this%NHBondN(h)=this%NHBondN(h) + 1
        elseif ( Counter(h,i) < 100 ) then
          do while (Counter(h,i) > 0 )
            m = m + 1
            Counter(h,i) = Counter(h,i) - 1
          end do
          j = m
          this%NHBond1(h,j)=this%NHBond1(h,j) + 1
        elseif ( Counter(h,i) < 10000 ) then
          do while (Counter(h,i) > 100 )
            m = m + 1
            Counter(h,i) = Counter(h,i) - 100
          end do
          j = m
          m = 0
          do while (Counter(h,i) > 0 )
            m = m + 1
            Counter(h,i) = Counter(h,i) - 1
          end do
          if (j .le. m) then
            k = m
          else
            k = j
            j = m
          end if
          this%NHBond2(h,j,k)=this%NHBond2(h,j,k) + 1
        else
          do while (Counter(h,i) > 10000 )
            m = m + 1
            Counter(h,i) = Counter(h,i) - 10000
          end do
          j = m
          m = 0
          do while (Counter(h,i) > 100 )
            m = m + 1
            Counter(h,i) = Counter(h,i) - 100
          end do
          if (j .le. m) then
            k = m
          else
            k = j
            j = m
          end if
          m = 0
          do while (Counter(h,i) > 0 )
            m = m + 1
            Counter(h,i) = Counter(h,i) - 1
          end do
          if (k .le. m) then
            l = m
          else
            l = k
            if (j .le. m) then
              k = m
            else
              k = j
              j = m
            end if
          end if
          this%NHBond3(h,j,k,l)=this%NHBond3(h,j,k,l) + 1
        end if

      end do
    end do

#if MPI_VER > 0
    if (SimulationType .eq. MolecularDynamics) then
      call MPI_Reduce( this%NHBond0(:), NHBAll(:), this%NComponents, MPI_INTEGER, MPI_SUM, NRootProc, Communicator, ierror )
      if (RootProc) this%NHBond0(:) = NHBAll(:)
      do j = 1, this%NComponents
        call MPI_Reduce( this%NHBond1(:,j), NHBAll(:), this%NComponents, MPI_INTEGER, MPI_SUM, NRootProc, Communicator, ierror )
        if (RootProc) this%NHBond1(:,j) = NHBAll(:)
        do k = j, this%NComponents
          call MPI_Reduce( this%NHBond2(:,j,k), NHBAll(:), this%NComponents, MPI_INTEGER, MPI_SUM, NRootProc, Communicator, ierror )
          if (RootProc) this%NHBond2(:,j,k) = NHBAll(:)
          do l = k, this%NComponents
            call MPI_Reduce( this%NHBond3(:,j,k,l), NHBAll(:), this%NComponents, MPI_INTEGER, MPI_SUM, NRootProc, Communicator, ierror )
            if (RootProc) this%NHBond3(:,j,k,l) = NHBAll(:)
          end do
        end do
      end do
      call MPI_Reduce( this%NHBondN(:), NHBAll(:), this%NComponents, MPI_INTEGER, MPI_SUM, NRootProc, Communicator, ierror )
      if (RootProc) this%NHBondN(:) = NHBAll(:)
    end if
#endif
    !this%NHBondN(1)=this%NPart-this%NHBond0(1)-this%NHBond1(1,1)-this%NHBond2(1,1,1)-this%NHBond3(1,1,1,1)

!      !Output of the H-bonded Molecules
!      if( (StepTotal > 1) .and. (mod( StepTotal - 1, VisualUpdateFrequency ) == 0) ) then  
!        call VisualUpdate( this, np, MH(:), MO(:,:) )
!      endif

  end subroutine TEnsemble_HBonding
#endif


!==============================================================!
!  Subroutine TEnsemble_PredictVol                             !
!==============================================================!

  subroutine TEnsemble_PredictVol( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK) :: BoxLengthOld, DelBoxL

    ! Predict volume of simulation box
    if ( RootProc ) then
      ! Call predictor
      select case( IntegratorType )
      case( IntegratorTypeGear )

        this%Volume0 = this%Volume0 + this%Volume1 + this%Volume2 + this%Volume3 + this%Volume4 + this%Volume5
        this%Volume1 = this%Volume1 + 2._RK * this%Volume2 + 3._RK * this%Volume3 &
&                    + 4._RK * this%Volume4 + 5._RK * this%Volume5
        this%Volume2 = this%Volume2 + 3._RK * this%Volume3 + 6._RK * this%Volume4 &
&                    + 10._RK * this%Volume5
        this%Volume3 = this%Volume3 + 4._RK * this%Volume4 + 10._RK * this%Volume5
        this%Volume4 = this%Volume4 + 5._RK * this%Volume5

      case( IntegratorTypeLeapFrog )
        this%Volume1 = this%Volume1 + this%Volume2
        this%Volume0 = this%Volume0 + this%Volume1

      case( IntegratorTypeVerlet )

      case( IntegratorTypeVV )

      end select
    end if

#if MPI_VER > 0
    ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
    call MPI_Bcast( this%Volume0, 1, MPI_RK, NRootProc, Communicator, ierror )
#endif
    BoxLengthOld = this%BoxLength
    call UpdateBoxLength( this )

    DelBoxL = this%BoxLength / BoxLengthOld

  end subroutine TEnsemble_PredictVol


!==============================================================!
!  Subroutine TEnsemble_ChangeFluctTI                          !
!==============================================================!

  subroutine TEnsemble_ChangeFluctTI( this, nt, nc )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)        :: this
    integer, intent(in)    :: nt
    integer, intent(in)    :: nc

    ! Declare local variables
    type(TComponent), pointer   :: pt, pc
    type(TInteraction), pointer :: pti, pci
    integer                     :: i, k, n, nu, n1, nu2, nu2k
    real(RK)                    :: PSave(3)
    real(RK)                    :: P0Save(3, 1:this%Component(nc)%Molecule%NUnit)
    real(RK)                    :: Q0Save(4, 1:this%Component(nc)%Molecule%NUnit)
    real(RK)                    :: ESave(this%NUnitMax,this%NPartMax*this%NUnitMax)
    real(RK)                    :: VSave(this%NUnitMax,this%NPartMax*this%NUnitMax)

    ! Assign local variables
    pt => this%Component(nt)
    pc => this%Component(nc)

    if ( SimulationType .eq. MonteCarlo .or. (SimulationType .eq. MolecularDynamics .and. RootProc) ) then
      n1 = rnd(pc%NPart)
      nu = pc%Molecule%NUnit
      write( IOBuffer, '("Exchanging fluctuating particle with particle ", I3, " of the corresponding TI component", I3)' ) n1, nc
      call LogWrite

      ! Copy position and quaternions
      PSave(:) = pt%Pm0(1, :)
      pt%Pm0(1, :) = pc%Pm0(n1, :)
      pc%Pm0(n1, :) = PSave(:)
      P0Save(:,1:nu) = pt%P0(1,:,1:nu )
      pt%P0(1,:,1:nu) = pc%P0(n1,:,1:nu)
      pc%P0(n1,:,1:nu) = P0Save(:,1:nu)

      if( pc%Molecule%IsElongated ) then
        Q0Save(:, 1:nu) = pt%Q0(1, :, 1:nu)
        pt%Q0(1, :, 1:nu) = pc%Q0(n1, :, 1:nu)
        pc%Q0(n1, :, 1:nu) = Q0Save(:, 1:nu)
      end if
    end if
#if MPI_VER > 0
    if (SimulationType .eq. MolecularDynamics) then
      call MPI_Bcast( n1, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
    end if
#endif

    ! Convert molecular coordinates to atom positions
    call Unit2Atom1( pt, 1 )
    call Unit2Atom1( pc, n1 )

    ! Copy energies and virial
    if (SimulationType .eq. MonteCarlo) then
      do i = 1, this%NRealComponents
        pti => this%Interaction(nt, i)
        pci => this%Interaction(nc, i)
        n = pci%NPart2*pci%NUnit2
        nu2 = (n1-1)*pci%NUnit1
        do k=1, pci%NUnit1
          nu2k = nu2 + k
          ESave(k,1:n) = pti%EPot(k, :)
          if ( this%OptPressure ) then
            VSave(k,1:n) = pti%Virial(k, :)
          end if
          pti%EPot(k, :) = pci%EPot(nu2k, :)
          this%Interaction(i, nt)%EPot(:, k) = pci%EPot(nu2k, :)
          pci%EPot(nu2k, :) = ESave(k,1:n)
          this%Interaction(i, nc)%EPot(:, nu2k) = ESave(k,1:n)
          if ( this%OptPressure ) then
            pti%Virial(k, :) = pci%Virial(nu2k, :)
            this%Interaction(i, nt)%Virial(:, k) = pci%Virial(nu2k, :)
            pci%Virial(nu2k, :) = VSave(k,1:n)
            this%Interaction(i, nc)%Virial(:, nu2k) = VSave(k,1:n)
          end if
        end do
      end do
    end if

  end subroutine TEnsemble_ChangeFluctTI


!==============================================================!
!  Subroutine TEnsemble_RotateMol_NPH                          !
!==============================================================!

  subroutine TEnsemble_RotateMol_NPH( this, nc, np )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments03
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np

    ! Declare local variables
    real(RK)                  :: p(3, this%Component(nc)%Molecule%NUnit)
    real(RK)                  :: q(4, this%Component(nc)%Molecule%NUnit)
    real(RK)                  :: dq(3), EPotOld, EPotNew
    real(RK)                  :: EFourier, EVirial
    type(TComponent), pointer :: pc
    integer                   :: i, NUnit
    real(RK)                  :: EPotDelta
    logical                   :: accepted

    ! Assign local variables
    pc => this%Component(nc)
    NUnit = pc%Molecule%NUnit

    ! Update number of rotation attempts
    pc%NRotateAttempts = pc%NRotateAttempts + 1

    ! Save old positions
    do i=1,NUnit
      p(:,i) = pc%P0(np, :, i)
      q(:,i) = pc%Q0(np, :, i)
    end do

    ! Calculate old Energies
    EPotOld = GetEnergy( this, nc, np )   ! IDF

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO

#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      call chargegrid_min(this, nc, np)
#endif
    end if

    ! Generate a trial rotation
    do i = 1, 3
      dq(i) = rnd( -pc%DispMolRot, pc%DispMolRot )
    end do

    ! Calculate new unit and atom positions
    call RotateMol(pc,np,dq)
    call Unit2Atom1(pc,np)

#if SPME > 0
    if (LongRange .eq. PME) then
      call chargegrid_plus(this, nc, np)
    end if
#endif

    ! Calculate particle energy with trial orientation
    call Energy( this, nc, np, EPotNew )

    ! Apply acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, MPI_RK, MPI_SUM, Communicator, ierror )
    else
      EPotDelta = EPotOld - EPotNew
    endif
#else
    EPotDelta = EPotOld - EPotNew
#endif

     ! Acceptance criterion
    if( exp(( real (this%NDF, RK) / 2._RK - 1._RK) * log((this%RefEnthalpy*this%NPart - this%Epot+EpotDelta - this%RefPressure * this%Volume0) &
&       / (this%RefEnthalpy*this%NPart - this%Epot - this%RefPressure * this%Volume0))) > rnd( 0._RK, 1._RK ) ) then

     ! Accept rotation
      this%Temperature = 2._RK * (this%RefEnthalpy*this%NPart - this%Epot+EpotDelta - this%RefPressure * this%Volume0) / real (this%NDF, RK)
      pc%NRotateMolSuccesses = pc%NRotateMolSuccesses + 1
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
        do i=1,NUnit
          pc%P0(np, :, i) = p(:,i)
          pc%Q0(np, :, i) = q(:,i)
          call Unit2Atom1( pc, np, i )
        end do
        call EwaldFourierEnergy(this,nc,np)

#if SPME > 0
      else if (LongRange .eq. PME) then
        this%UFourier = EFourier
        this%EVirial  = EVirial
        call chargegrid_min(this, nc, np)
        do i=1,NUnit
          pc%P0(np, :, i) = p(:,i)
          pc%Q0(np, :, i) = q(:,i)
          call Unit2Atom1( pc, np, i )
        end do
        call chargegrid_plus(this, nc, np)
#endif

      else
        do i=1,NUnit
          pc%P0(np, :, i) = p(:,i)
          pc%Q0(np, :, i) = q(:,i)
          call Unit2Atom1( pc, np, i )
        end do
      end if

    end if

  end subroutine TEnsemble_RotateMol_NPH


!==============================================================!
!  Subroutine TEnsemble_MoveMol_NPH                            !
!==============================================================!

  subroutine TEnsemble_MoveMol_NPH( this, nc, np )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np

    ! Declare local variables
    real(RK)                  :: rm(3), trans(3)
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier, EVirial
    real(RK)                  :: EPotDelta
    type(TComponent), pointer :: pc
    integer                   :: i, j, NUnit
    logical                   :: accepted

    ! Assign local variables
    pc => this%Component(nc)
    NUnit = pc%Molecule%NUnit

    ! Update number of move attempts
    pc%NMoveAttempts = pc%NMoveAttempts + 1

    ! Save current particle position and energy
    rm(:) = pc%Pm0(np, :)
    EPotOld = GetEnergy( this, nc, np )

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO

#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      call chargegrid_min(this, nc, np)
#endif
    end if

    ! Generate a trial displacement & Apply periodic boundary conditions
    do i = 1, 3
      trans(i) = rnd( -pc%DispMolTran, pc%DispMolTran )
      pc%Pm0(np, i) = pc%Pm0(np, i) + trans(i)
      pc%Pm0(np, i) = pc%Pm0(np, i) - anint( pc%Pm0(np, i) )
      do j=1, NUnit
        pc%P0(np, i, j ) = pc%P0(np, i, j ) + trans(i)
        pc%P0(np, i, j ) = pc%P0(np, i, j ) - anint( pc%P0(np, i, j) )
      end do
    end do
    
    ! Convert unit coordinates to atom positions and calculate Energies
    call Unit2Atom1( pc, np )

#if SPME > 0
    ! Calculate changes in the SPME grid
    if (LongRange .eq. PME) then
      call chargegrid_plus(this, nc, np)
    end if
#endif

    ! Calculate particle energy at trial position
    call Energy( this, nc, np, EPotNew )
    ! Apply acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, MPI_RK, MPI_SUM, Communicator, ierror )
    else
          EPotDelta = EPotOld - EPotNew
    endif
#else
     EPotDelta = EPotOld - EPotNew
#endif

     ! Acceptance criterion
    if( exp(( real (this%NDF, RK) / 2._RK  - 1._RK) * log((this%RefEnthalpy*this%NPart - this%Epot+EpotDelta - this%RefPressure * this%Volume0) &
&       / (this%RefEnthalpy*this%NPart - this%Epot - this%RefPressure * this%Volume0))) > rnd( 0._RK, 1._RK ) ) then
!print*, 'MOVE', real (this%NDF, RK), this%RefEnthalpy, this%Epot, EpotDelta, this%RefPressure, this%Volume0

     ! Accept move
      this%Temperature = 2._RK * (this%RefEnthalpy*this%NPart - this%Epot+EpotDelta - this%RefPressure * this%Volume0) / real (this%NDF, RK)
      pc%NMoveMolSuccesses = pc%NMoveMolSuccesses + 1
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

        pc%Pm0(np, :) = rm(:)
        do j=1, NUnit
          pc%P0(np, :, j) = pc%P0(np, :, j) - trans(:)
          pc%P0(np, :, j) = pc%P0(np, :, j) - anint( pc%P0(np, :, j) )
          call Unit2Atom1( pc, np, j )
        end do
        call EwaldFourierEnergy(this,nc,np)

#if SPME > 0
      else if (LongRange .eq. PME) then
        this%UFourier = EFourier
        this%EVirial  = EVirial
        call chargegrid_min(this, nc, np)
        pc%Pm0(np, :) = rm(:)
        do j=1, NUnit
          pc%P0(np, :, j) = pc%P0(np, :, j) - trans(:)
          pc%P0(np, :, j) = pc%P0(np, :, j) - anint( pc%P0(np, :, j) )
          call Unit2Atom1( pc, np, j )
        end do
        call chargegrid_plus(this, nc, np)
#endif
      else
        pc%Pm0(np, :) = rm(:)
        do j=1, NUnit
          pc%P0(np, :, j) = pc%P0(np, :, j) - trans(:)
          pc%P0(np, :, j) = pc%P0(np, :, j) - anint( pc%P0(np, :, j) )
          call Unit2Atom1( pc, np, j )
        end do
      end if

    end if

  end subroutine TEnsemble_MoveMol_NPH


!==============================================================!
!  Subroutine TEnsemble_RotateMol_NVE                          !
!==============================================================!

  subroutine TEnsemble_RotateMol_NVE( this, nc, np )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments03
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np

    ! Declare local variables
    real(RK)                  :: p(3, this%Component(nc)%Molecule%NUnit)
    real(RK)                  :: q(4, this%Component(nc)%Molecule%NUnit)
    real(RK)                  :: dq(3), NewOmega
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier
#if SPME > 0
    real(RK)                  :: EVirial
#endif
    type(TComponent), pointer :: pc
    integer                   :: i
    integer                   :: NUnit
    real(RK)                  :: EPotDelta

    ! Assign local variables
    pc => this%Component(nc)
    NUnit = pc%Molecule%NUnit

    ! Update number of rotation attempts
    pc%NRotateMolAttempts = pc%NRotateMolAttempts + 1

    ! Save old positions
    do i=1,NUnit
      p(:,i) = pc%P0(np, :, i)
      q(:,i) = pc%Q0(np, :, i)
    end do

    ! Calculate old Energies
    EPotOld = GetEnergy( this, nc, np )   ! IDF

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO
#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      call chargegrid_min(this, nc, np)
#endif
    end if

    ! Generate a trial rotation
    do i = 1, 3
      dq(i) = rnd( -pc%DispMolRot, pc%DispMolRot )
    end do

    ! Calculate new unit and atom positions
    call RotateMol(pc,np,dq)
    call Unit2Atom1(pc,np)

#if SPME > 0
    if (LongRange .eq. PME) then
      call chargegrid_plus(this, nc, np)
    end if
#endif

    ! Calculate particle energy with trial orientation
    call Energy( this, nc, np, EPotNew )

    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, MPI_RK, MPI_SUM, Communicator, ierror )
    else
      EPotDelta = EPotOld - EPotNew
    endif
#else
    EPotDelta = EPotOld - EPotNew
#endif

    if( (this%RefHamiltonian*this%NPart - this%Epot+EPotDelta) < 0._RK ) then
      NewOmega = 0._RK
    else
      NewOmega = 1._RK
    end if

    if( ((this%RefHamiltonian*this%NPart - this%Epot+EPotDelta)/(this%RefHamiltonian*this%NPart - this%Epot))**((real (this%NDF-this%constrNDF, RK)-2._RK)/2._RK) &
&         * NewOmega .ge. rnd( 0._RK, 1._RK ) ) then
      ! Accept rotation
      pc%NRotateMolSuccesses = pc%NRotateMolSuccesses + 1
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
        do i=1,NUnit
          pc%P0(np, :, i) = p(:,i)
          pc%Q0(np, :, i) = q(:,i)
          call Unit2Atom1( pc, np, i )
        end do
        call EwaldFourierEnergy(this,nc,np)

#if SPME > 0
      else if (LongRange .eq. PME) then
        this%UFourier = EFourier
        this%EVirial  = EVirial
        call chargegrid_min(this, nc, np)
        do i=1,NUnit
          pc%P0(np, :, i) = p(:,i)
          pc%Q0(np, :, i) = q(:,i)
          call Unit2Atom1( pc, np, i )
        end do
        call chargegrid_plus(this, nc, np)
#endif
      else
        do i=1,NUnit
          pc%P0(np, :, i) = p(:,i)
          pc%Q0(np, :, i) = q(:,i)
          call Unit2Atom1( pc, np, i )
        end do
      end if

    end if

  end subroutine TEnsemble_RotateMol_NVE


!==============================================================!
!  Subroutine TEnsemble_MoveMol_NVE                            !
!==============================================================!

  subroutine TEnsemble_MoveMol_NVE( this, nc, np )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np

    ! Declare local variables
    real(RK)                  :: rm(3), trans(3)
    real(RK)                  :: TransMove, NewOmega
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier
#if SPME > 0
    real(RK)                  :: EVirial
#endif
    real(RK)                  :: EPotDelta
    type(TComponent), pointer :: pc
    integer                   :: i, j
    integer                   :: NUnit

    ! Assign local variables
    pc => this%Component(nc)
    NUnit = pc%Molecule%NUnit

    ! Update number of move attempts
    pc%NMoveMolAttempts = pc%NMoveMolAttempts + 1

    ! Save current particle position and energy
    rm(:) = pc%Pm0(np, :)
    EPotOld = GetEnergy( this, nc, np )   ! IDF

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO
#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      call chargegrid_min(this, nc, np)
#endif
    end if

    ! Generate a trial displacement & Apply periodic boundary conditions
    do i = 1, 3
      trans(i) = rnd( -pc%DispMolTran, pc%DispMolTran )
      pc%Pm0(np, i) = pc%Pm0(np, i) + trans(i)
      pc%Pm0(np, i) = pc%Pm0(np, i) - anint( pc%Pm0(np, i) )
      do j=1, NUnit
        pc%P0(np, i, j ) = pc%P0(np, i, j ) + trans(i)
        pc%P0(np, i, j ) = pc%P0(np, i, j ) - anint( pc%P0(np, i, j) )
      end do

    end do

    ! Convert unit coordinates to atom positions and calculate Energies
    call Unit2Atom1( pc, np )


#if SPME > 0
    ! Calculate changes in the SPME grid
    if (LongRange .eq. PME) then
      call chargegrid_plus(this, nc, np)
    end if
#endif

    ! Calculate particle energy at trial position
    call Energy( this, nc, np, EPotNew )

    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, MPI_RK, MPI_SUM, Communicator, ierror )
    else
      EPotDelta = EPotOld - EPotNew
    endif
#else
    EPotDelta = EPotOld - EPotNew
#endif

    if( (this%RefHamiltonian*this%NPart - this%Epot+EPotDelta) < 0._RK ) then
      NewOmega = 0._RK
    else
      NewOmega = 1._RK
    end if

    if( ((this%RefHamiltonian*this%NPart - this%Epot+EPotDelta)/(this%RefHamiltonian*this%NPart - this%Epot))**((real (this%NDF-this%constrNDF, RK)-2._RK)/2._RK) &
&         * NewOmega .ge. rnd( 0._RK, 1._RK ) ) then
      ! Accept move
      pc%NMoveMolSuccesses = pc%NMoveMolSuccesses + 1
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
          pc%Pm0(np, :) = rm(:)
          do j=1, NUnit
            pc%P0(np, :, j) = pc%P0(np, :, j) - trans(:)
            pc%P0(np, :, j) = pc%P0(np, :, j) - anint( pc%P0(np, :, j) )
            call Unit2Atom1( pc, np, j )
          end do
          call EwaldFourierEnergy(this,nc,np)

#if SPME > 0
      else if (LongRange .eq. PME) then
          this%UFourier = EFourier
          this%EVirial  = EVirial
          call chargegrid_min(this, nc, np)
          pc%Pm0(np, :) = rm(:)
          do j=1, NUnit
            pc%P0(np, :, j) = pc%P0(np, :, j) - trans(:)
            pc%P0(np, :, j) = pc%P0(np, :, j) - anint( pc%P0(np, :, j) )
            call Unit2Atom1( pc, np, j )
          end do
          call chargegrid_plus(this, nc, np)
#endif
      else
        pc%Pm0(np, :) = rm(:)
        do j=1, NUnit
          pc%P0(np, :, j) = pc%P0(np, :, j) - trans(:)
          pc%P0(np, :, j) = pc%P0(np, :, j) - anint( pc%P0(np, :, j) )
          call Unit2Atom1( pc, np, j )
        end do
      end if

    end if

  end subroutine TEnsemble_MoveMol_NVE


!==============================================================!
!  Subroutine TEnsemble_RotateMol                              !
!==============================================================!

  subroutine TEnsemble_RotateMol( this, nc, np )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments03
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np

    ! Declare local variables
    real(RK)                  :: p(3, this%Component(nc)%Molecule%NUnit)
    real(RK)                  :: q(4, this%Component(nc)%Molecule%NUnit)
    real(RK)                  :: dq(3)
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier
#if SPME > 0
    real(RK)                  :: EVirial
#endif
    type(TComponent), pointer :: pc
    integer                   :: i
    integer                   :: NUnit
    real(RK)                  :: EPotDelta
    logical                   :: accepted

    ! Assign local variables
    pc => this%Component(nc)
    NUnit = pc%Molecule%NUnit

    ! Update number of rotation attempts
    pc%NRotateMolAttempts = pc%NRotateMolAttempts + 1

    ! Save old positions
    do i=1,NUnit
      p(:,i) = pc%P0(np, :, i)
      q(:,i) = pc%Q0(np, :, i)
    end do
    ! Calculate old Energies
    EPotOld = GetEnergy( this, nc, np )   ! IDF

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO
#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      call chargegrid_min(this, nc, np)
#endif
    end if

    ! Generate a trial rotation
    do i = 1, 3
      dq(i) = rnd( -pc%DispMolRot, pc%DispMolRot )
    end do

    ! Calculate new unit and atom positions
    call RotateMol(pc,np,dq)
    call Unit2Atom1( pc, np )

#if SPME > 0
    if (LongRange .eq. PME) then
      call chargegrid_plus(this, nc, np)
    end if
#endif

    ! Calculate particle energy with trial orientation
    call Energy( this, nc, np, EPotNew )

    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, MPI_RK, MPI_SUM, Communicator, ierror )
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
      pc%NRotateMolSuccesses = pc%NRotateMolSuccesses + 1
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
        do i=1,NUnit
          pc%P0(np, :, i) = p(:,i)
          pc%Q0(np, :, i) = q(:,i)
          call Unit2Atom1( pc, np, i )
        end do
        call EwaldFourierEnergy(this,nc,np)

#if SPME > 0
      else if (LongRange .eq. PME) then
        this%UFourier = EFourier
        this%EVirial  = EVirial
        call chargegrid_min(this, nc, np)
        do i=1,NUnit
          pc%P0(np, :, i) = p(:,i)
          pc%Q0(np, :, i) = q(:,i)
          call Unit2Atom1( pc, np, i )
        end do
        call chargegrid_plus(this, nc, np)
#endif
      else
        do i=1,NUnit
          pc%P0(np, :, i) = p(:,i)
          pc%Q0(np, :, i) = q(:,i)
          call Unit2Atom1( pc, np, i )
        end do
      end if

    end if

  end subroutine TEnsemble_RotateMol


!==============================================================!
!  Subroutine TEnsemble_MoveMol                                !
!==============================================================!

  subroutine TEnsemble_MoveMol( this, nc, np )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np

    ! Declare local variables
    real(RK)                  :: rm(3), trans(3)
    real(RK)                  :: r(3,this%Component(nc)%Molecule%NUnit)
    real(RK)                  :: TransMove
    real(RK)                  :: EPotOld, EPotNew
    real(RK)                  :: EFourier
#if SPME > 0
    real(RK)                  :: EVirial
#endif
    real(RK)                  :: EPotDelta
    type(TComponent), pointer :: pc
    integer                   :: i, j, NUnit
    logical                   :: accepted

    ! Assign local variables
    pc => this%Component(nc)
    NUnit = pc%Molecule%NUnit

    ! Update number of move attempts
    pc%NMoveMolAttempts = pc%NMoveMolAttempts + 1

    ! Save current particle position and energy
    rm(:) = pc%Pm0(np, :)
    r(:,:) = pc%P0(np,:,:)
    EPotOld = GetEnergy( this, nc, np )   ! IDF

    ! Save the Energies and Virials for a faster MoveRejction
    if (LongRange .eq. Ewald) then
      EFourier = this%UFourier
      DO i=1,pc%Molecule%NCharge
        this%rold(i,1) = pc%Molecule%SiteCharge(i)%RX(np)
        this%rold(i,2) = pc%Molecule%SiteCharge(i)%RY(np)
        this%rold(i,3) = pc%Molecule%SiteCharge(i)%RZ(np)
      END DO
#if SPME > 0
    else if (LongRange .eq. PME) then
      EFourier = this%UFourier
      EVirial  = this%EVirial
      call chargegrid_min(this, nc, np)
#endif
    end if

    ! Generate a trial displacement & Apply periodic boundary conditions
    do i = 1, 3
      trans(i) = rnd( -pc%DispMolTran, pc%DispMolTran )
      pc%Pm0(np, i) = pc%Pm0(np, i) + trans(i)
      do j=1, NUnit
        pc%P0(np, i, j ) = pc%P0(np, i, j ) + trans(i)
        pc%P0(np, i, j ) = pc%P0(np, i, j ) - anint( pc%P0(np, i, j) )
      end do
      pc%Pm0(np, i) = pc%Pm0(np, i) - anint( pc%Pm0(np, i) )
    end do

    ! Convert molecular coordinates to atom positions and calculate Energies
    call Unit2Atom1( pc, np )

#if SPME > 0
    ! Calculate changes in the SPME grid
    if (LongRange .eq. PME) then
      call chargegrid_plus(this, nc, np)
    end if
#endif

    ! Calculate particle energy at trial position
    call Energy( this, nc, np, EPotNew )

    ! Apply Metropolis acceptance criterion
#if MPI_VER > 0
    if ( Equilibration .and. CommonEqui ) then
      ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
      call MPI_Allreduce( EPotOld - EPotNew, EPotDelta, 1, MPI_RK, MPI_SUM, Communicator, ierror )
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
      pc%NMoveMolSuccesses = pc%NMoveMolSuccesses + 1
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
          pc%Pm0(np, :) = rm(:)
          pc%P0(np, :, :) = r(:,:)
          call Unit2Atom1( pc, np )
          call EwaldFourierEnergy(this,nc,np)
#if SPME > 0
      else if (LongRange .eq. PME) then
          this%UFourier = EFourier
          this%EVirial  = EVirial
          call chargegrid_min(this, nc, np)
          pc%Pm0(np, :) = rm(:)
          pc%P0(np, :, :) = r(:,:)
          call Unit2Atom1( pc, np )
          call chargegrid_plus(this, nc, np)
#endif
      else
        pc%Pm0(np, :) = rm(:)
        pc%P0(np, :, :) = r(:,:)
        call Unit2Atom1( pc, np )
      end if

    end if

  end subroutine TEnsemble_MoveMol


!==============================================================!
!  Function TEnsemble_GetVirialIntra                           !
!==============================================================!

  function TEnsemble_GetVirialIntra( this ) result(V)

    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare result
    real(RK) :: V

    ! Declare local variables
    integer :: i, j
    integer :: n
    integer :: NUnit, np

    ! Calculate potential energy of a particle
    V = 0._RK
    do i = 1, this%NComponents
      NUnit = this%Component(i)%Molecule%NUnit
      np = this%Component(i)%NPart
      n = np*NUnit
      do j = 1, np
        V = V + sum( this%Interaction(i, i)% &
&         Virial((j-1)*NUnit+1:j*NUnit,(j-1)*NUnit+1:j*NUnit) )
      end do
    end do
    V = .5_RK * V

  end function TEnsemble_GetVirialIntra


!==============================================================!
!  Function TEnsemble_GetEnergyIntra                           !
!==============================================================!

  function TEnsemble_GetEnergyIntra( this ) result(E)
    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare result
    real(RK) :: E

    ! Declare local variables
    integer :: i, j, np, nu
    real(RK):: Intra

    ! Calculate potential energy of a particle
    E = 0._RK
    if ( UseIntDegFreed ) then
      Intra = 0._RK
      do i = 1, this%NComponents
        nu = this%Component(i)%Molecule%NUnit
        np = this%Component(i)%NPart
        do j=1,np
          E = E + sum( this%Interaction(i, i)%EPot((j-1)*nu+1:j*nu,(j-1)*nu+1:j*nu) )
        end do

        ! Kein Faktor 2, weil unten einfach aufaddiert wird
        Intra = Intra + sum(this%Interaction(i,i)%EPotAngle(:)) + sum(this%Interaction(i,i)%EPotTo(:))
      end do
      E = .5_RK * E + Intra
    endif

  end function TEnsemble_GetEnergyIntra


!==============================================================!
!  Function TEnsemble_GetEnergyIntra1Mol (per molecule)        !
!==============================================================!

  function TEnsemble_GetEnergyIntra1Mol( this, nc, np ) result(E)
    implicit none

    ! Declare arguments
    type(TEnsemble) :: this
    integer, intent(in) :: nc, np

    ! Declare result
    real(RK) :: E

    ! Declare local variables
    integer :: nu

    ! Calculate potential energy of a particle
    E = 0._RK
    if ( UseIntDegFreed ) then
      nu = this%Component(nc)%Molecule%NUnit
      E = 0.5_RK * sum( this%Interaction(nc, nc)%EPot((np-1)*nu+1:np*nu,(np-1)*nu+1:np*nu) )
      if (associated(this%Component(nc)%Molecule%idfangle)) then   !Michael Sch.: assoicated terms only needed here,
        E = E + this%Interaction(nc,nc)%EPotAngle(np)              !          since : omits '0' entries/empty arrays
      endif
      if (associated(this%Component(nc)%Molecule%idfdihedral)) then
        E = E + this%Interaction(nc,nc)%EPotTo(np)
      endif
    endif

  end function TEnsemble_GetEnergyIntra1Mol


 !==============================================================!
 !  Function TEnsemble_GetEnergyIntra_Bond                      !
 !==============================================================!

  function TEnsemble_GetEnergyIntra_Bond( this ) result(E)
    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare result
    real(RK) :: E

    ! Declare local variables
    integer :: i, j, nu, np

    ! Calculate potential energy of a particle
    E = 0._RK
    do i = 1, this%NComponents
      nu = this%Component(i)%Molecule%NUnit
      np = this%Component(i)%NPart
      do j=1,np
        E = E + sum( this%Interaction(i, i)%EPot((j-1)*nu+1:j*nu,(j-1)*nu+1:j*nu) )
      end do
    end do
    E = .5_RK * E

  end function TEnsemble_GetEnergyIntra_Bond


 !==============================================================!
 !  Function TEnsemble_GetEnergyIntra_Angle                     !
 !==============================================================!

  function TEnsemble_GetEnergyIntra_Angle( this ) result(E)
    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare result
    real(RK) :: E

    ! Declare local variables
    integer :: i

    ! Calculate potential energy of a particle
    E = 0._RK
    do i = 1, this%NComponents
      E = E + sum(this%Interaction(i,i)%EPotAngle(:))
    end do

  end function TEnsemble_GetEnergyIntra_Angle


 !==============================================================!
 !  Function TEnsemble_GetEnergyIntra_Dihedral                  !
 !==============================================================!

  function TEnsemble_GetEnergyIntra_Dihedral( this ) result(E)
    implicit none

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare result
    real(RK) :: E

    ! Declare local variables
    integer :: i

    ! Calculate potential energy of a particle
    E = 0._RK
    do i = 1, this%NComponents
      E = E + sum(this%Interaction(i,i)%EPotTo(:))
    end do

  end function TEnsemble_GetEnergyIntra_Dihedral


!==============================================================!
!  Function TEnsemble_GetEnergy1Mol (per molecule)             !
!==============================================================!

  function TEnsemble_GetEnergy1Mol( this, nc, np ) result(E)

    implicit none

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np

    ! Declare result
    real(RK) :: E

    ! Declare local variables
    integer :: i, k
    integer :: NAngle, NDihedral
    integer :: NAngleNum, NDihedralNum
    integer :: NUnitPart
    integer :: nu, nup1

    ! Calculate potential energy of a particle
    E = 0._RK
    nu = this%Component(nc)%Molecule%NUnit
    nup1 = nu * (np - 1)
    do i = 1, this%NComponents
      NUnitPart = this%Component(i)%Molecule%NUnit*this%Component(i)%NPart
      do k=1, nu
        E = E + sum( this%Interaction(i, nc)%EPot(1:NUnitPart, nup1+k) )
      end do
    end do

    if ( UseIntDegFreed ) then
      E = E - 0.5_RK*sum( this%Interaction(nc,nc)%EPot(nup1+1:nup1+nu,nup1+1:nup1+nu) )
      NAngle = this%Interaction(nc,nc)%NAngle
      NDihedral = this%Interaction(nc,nc)%NDihedral
      NAngleNum = (np-1)*NAngle
      NDihedralNum = (np-1)*NDihedral
      E = E + sum(this%Interaction(nc,nc)%EPotAngle(NAngleNum+1:NAngleNum+NAngle)) + &
&       sum(this%Interaction(nc,nc)%EPotTo(NDihedralNum +1:NDihedralNum +NDihedral))
    end if


    ! Ewald
    if (LongRange .eq. Ewald) then
      E = E + this%UFourier
#if SPME > 0
    else if (LongRange .eq. PME) then
      E = E + this%UFourier
#endif
    end if

  end function TEnsemble_GetEnergy1Mol


!==============================================================!
!  Subroutine TEnsemble_Energy1 (per molecule)                 !
!==============================================================!

  subroutine TEnsemble_Energy1Mol( this, nc, np, EPotNew )

    implicit none

    ! Declare arguments
    type(TEnsemble)       :: this
    integer, intent(in)   :: nc, np
    real(RK), intent(out) :: EPotNew

    ! Declare local variables
    type(TInteraction), pointer :: pi
    integer                     :: n, nu, nup
    integer                     :: i

    ! Initialize new energy
    EPotNew = 0._RK
    nup = (np-1)*this%Component(nc)%Molecule%NUnit

    ! Loop over components
    do i = 1, this%NComponents
      pi => this%Interaction(nc, i)
      n = pi%NPart2*pi%NUnit2
      do nu=1, this%Component(nc)%Molecule%NUnit
          call Energy( pi, np, nu, this%BoxLength )
          if ( pi%SameComponent .and. UseIntDegFreed ) then
            call IntraEnergy( pi, np, nu, this%BoxLength )
            EPotNew = EPotNew - 0.5_RK*sum( pi%EPot1(nup+1:nup+this%Component(nc)%Molecule%NUnit) )
          end if
          ! Calculate new energy
          EPotNew = EPotNew + sum( pi%EPot1(1:n) )  !includes Bond energies
          pi%EPotMol(nu,:) = pi%Epot1
          pi%d2EpotdV2Mol(nu, :) = pi%d2EpotdV21
          if (this%OptPressure) then
            pi%VirialMol(nu,:) = pi%Virial1
          end if
      end do
    end do

    !Michael Sch.: new form
    if ( UseIntDegFreed ) then
      pi => this%Interaction(nc, nc)
      EPotNew = EPotNew  + ( sum(pi%EPot1Angle) + sum(pi%EPot1To) )
    end if

    if (LongRange .eq. Ewald) then
       call EwaldFourierEnergy(this,nc,np)
       EPotNew = EPotnew + this%UFourier
#if SPME > 0
    else if (LongRange .eq. PME) then
       call PMEFourierTermMC( this )
       EPotNew = EPotnew + this%UFourier
#endif
    end if

  end subroutine TEnsemble_Energy1Mol



!==============================================================!
!  Subroutine TEnsemble_UpdateEnergy1Mol                       !
!==============================================================!

  subroutine TEnsemble_UpdateEnergy1Mol( this, nc, np )

    implicit none

    ! Declare arguments
    type(TEnsemble)     :: this
    integer, intent(in) :: nc, np

    ! Declare local variables
    type(TInteraction), pointer :: pi
    integer                     :: n
    integer                     :: i, j
    integer                     :: NBond, NAngle, NDihedral
    integer                     :: npu, npu1

    ! Update potential energy and virial matrices for a particle
    npu = (np-1) * this%Component(nc)%Molecule%NUnit

    do i = 1, this%NComponents
      pi => this%Interaction(nc, i)
      n = pi%NPart2 * pi%NUnit2
      do j=1,pi%NUnit1
        npu1 = npu + j
        pi%EPot(npu1, 1:n) = pi%EPotMol(j, 1:n)
        pi%d2EpotdV2(npu1, 1:n) = pi%d2EpotdV2Mol(j, 1:n)

        if ( this%OptPressure ) then
          pi%Virial(npu1, 1:n) = pi%VirialMol(j, 1:n)
        end if

        this%Interaction(i, nc)%EPot(1:n, npu1) = pi%EPotMol(j, 1:n)
        this%Interaction(i, nc)%d2EpotdV2(1:n, npu1) = pi%d2EpotdV2Mol(j, 1:n)

        if ( this%OptPressure ) then
          this%Interaction(i, nc)%Virial(1:n, npu1) = pi%VirialMol(j, 1:n)
        end if
      end do
    end do

    if ( UseIntDegFreed ) then
      pi => this%Interaction(nc,nc)
      NAngle = pi%NAngle
      NDihedral = pi%NDihedral

      pi%EPotAngle((np-1)*NAngle+1:np*NAngle) = pi%EPot1Angle(:)
      pi%EPotTo((np-1)*NDihedral+1:np*Ndihedral) = pi%EPot1To(:)
    end if

  end subroutine TEnsemble_UpdateEnergy1Mol


!==============================================================!
!  Subroutine TEnsemble_QShake                                 !
!==============================================================!

  subroutine TEnsemble_QShake( this )

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
    real(RK)                  :: VirialShake, tempVirial
    real(RK)                  :: dLogVolumeThird
    real(RK)                  :: oldF(this%NPartmax,3,this%NUnitmax)

    select case( IntegratorType )
    case( IntegratorTypeGear )
      call Error( "QShake only valid for Verlet-Algorithms" )
      ! QShake could be used with Gear, but precision advantage of Gear (o5) is lost when using QShake (o2)
      ! also Virial-contribution of Shake can't be accounted for in either correction or prediction step
    case( IntegratorTypeVerlet )
      call Error( "QShake only implemented for LeapFrog" )
    case( IntegratorTypeVV )
      call Error( "QShake only implemented for LeapFrog" )
    end select

    VirialShake = 0._RK
    tempVirial = 0._RK
    dLogVolumeThird = this%Volume1 / (3._RK * this%Volume0)

    ! calculate unconstrained and unscaled(T) positions
    this%scale = 1._RK ! shutoff Thermostat for unconstrained timestep
    call CorrectLeapFrog( this )
    call PredictLeapFrog( this )
    do i =1, this%NComponents
      pc => this%Component(i)
      if (RootProc) then
        oldF(:,:,:) = 0._RK
#if MPI_VER > 0
        oldF(1:pc%NPart,1:3,1:pc%Molecule%NUnit) = pc%FAll(1:pc%NPart,1:3,1:pc%Molecule%NUnit)
#else
        oldF(1:pc%NPart,1:3,1:pc%Molecule%NUnit) = pc%F(1:pc%NPart,1:3,1:pc%Molecule%NUnit)
#endif
      end if
      ! calculate new forces and positions due to constraints (bonds)
      call Constraints( pc, tempVirial )
      ! reverse unconstrained timestep
      if (RootProc) then
        call ReverseLeapFrog( pc, oldF(1:pc%NPart,1:3,1:pc%Molecule%NUnit), dLogVolumeThird )
      end if
    end do

#if MPI_VER > 0
    call MPI_Reduce( tempVirial, VirialShake, 1, MPI_RK, MPI_SUM, NRootProc, Communicator, ierror )
#else
    VirialShake = tempVirial
#endif

    VirialShake = Third * VirialShake
    this%Virial = this%Virial + VirialShake
    this%VirialIntra = this%VirialIntra + VirialShake
    this%Pressure = this%Pressure + VirialShake/this%Volume0

  end subroutine TEnsemble_QShake


!==============================================================!
!  Subroutine TEnsemble_CorrectVol                             !
!==============================================================!

  subroutine TEnsemble_CorrectVol( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TEnsemble) :: this

    ! Declare local variables
    real(RK) :: Volume2, Corr
    real(RK) :: BoxLengthOld, DelBoxL

#ifdef ABL
    real(RK) :: vol
    real(RK) :: fac
    real(RK) :: denom,denom2
    real(RK) :: nen
    integer  :: j
#endif

    ! Correct volume of simulation box
    if( RootProc ) then
      ! Call corrector
      select case( IntegratorType )
      case( IntegratorTypeGear )

        Volume2 = (this%Pressure - this%RefPressure) * TimeStepSquared2 / this%PistonMass
        Corr = Volume2 - this%Volume2
        this%Volume0 = this%Volume0 + Corr * Gear20
        this%Volume1 = this%Volume1 + Corr * Gear21
        this%Volume2 =      Volume2
        this%Volume3 = this%Volume3 + Corr * Gear23
        this%Volume4 = this%Volume4 + Corr * Gear24
        this%Volume5 = this%Volume5 + Corr * Gear25

#if ABL
        vol = this%Volume0 + this%Volume1 + this%Volume2 + this%Volume3 + this%Volume4 + this%Volume5
        fac = TimeStepSquared2*Gear20
        denom = fac*(this%Pressure - this%RefPressure) - this%PistonMass*this%Volume2*Gear20 ! Michael Sch.: per def = 0, also obsolet...
        denom2 = denom**2
        nen = this%PistonMass*fac / (vol * denom2)
        do i=1,this%NComponents
          do j=1,this%Component(i)%Molecule%NLJ126
            this%AblPS(i,j)   =  this%AblPS(i,j) + this%Interaction(1, 1)%PotLJ126LJ126(i, j)%AblSigCorr(i,j)
            this%AblPE(i,j)   =  this%AblPE(i,j) + this%Interaction(1, 1)%PotLJ126LJ126(i, j)%AblEpsCorr(i,j)
            this%AblRhoS(i,j) = nen * this%AblPS(i,j)
            this%AblRhoE(i,j) = nen * this%AblPE(i,j)
          end do
        end do
#endif

      case( IntegratorTypeLeapFrog )
        this%Volume2 = (this%Pressure - this%RefPressure) * TimeStepSquared2 / this%PistonMass
        this%Volume1 = this%Volume1 + this%Volume2

      case( IntegratorTypeVerlet )

      case( IntegratorTypeVV )
      end select

    end if

    if ( IntegratorType .eq. IntegratorTypeGear ) then
#if MPI_VER > 0
      ! use MPI_RK (cmp. ms2_global.F90) instead of MPI_RK
      call MPI_Bcast( this%Volume0, 1, MPI_RK, NRootProc, Communicator, ierror )
#endif
      BoxLengthOld = this%BoxLength
      call UpdateBoxLength( this )

      DelBoxL = this%BoxLength / BoxLengthOld

    end if


  end subroutine TEnsemble_CorrectVol


end module ms2_ensemble
