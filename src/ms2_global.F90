!==============================================================!
!  MOLECULAR SIMULATION PROGRAM MS2 Version 1.1 v12            !
!  (c) 2001 by Sergey Lishchuk, ITT                            !
!  (c) 2007 by Bernhard Eckl, ITT                              !
!  (c) 2008 by Ekaterina Elts, TUM                                  !
!==============================================================!
!  Module ms2_global                                           !
!  Contains declarations of global constants and functions     !
!==============================================================!

#ifndef ARCH
#define ARCH    0
#define FORTRAN 90
#define MPI_VER 0
#endif

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2_global.F90...'
#endif

module ms2_global

#ifdef _WIN32
  use dfport
#endif

#ifdef __INTEL_COMPILER
  use IFPORT
#endif



!==============================================================!
!  Global constants and variables                              !
!==============================================================!

  ! Define kind of real type
  ! 4: single precision
  ! 8: double precision
#ifdef SINGLEPRECISION
  integer, parameter :: RK = 4
#else
  integer, parameter :: RK = 8
#endif

  ! Define maximum length of file names
  integer, parameter :: FileNameLength = 128

  ! Name of program
#if ARCH == 1 || ARCH == 2 || ARCH == 3
  character(28)           :: ProgramFileName
#else
  character(*), parameter :: ProgramFileName = 'ms2'
#endif

  ! Version of program
#ifndef _WIN32
#ifdef __DATE__
  character(*), parameter :: VersionString = 'v12 ' // __DATE__
#else
  character(*), parameter :: VersionString = 'v12'
#endif
#else
  character(*), parameter :: VersionString = 'v12 Windows'
#endif

  ! Name of platform
#if ARCH == 1
  character(*), parameter :: Hardware = 'alpha'
#elif ARCH == 2
#ifdef _WIN32
  character(*), parameter :: Hardware = 'Win32'
#elif defined __INTEL_COMPILER
  character(*), parameter :: Hardware = 'Linux/Ifort'
#elif defined _PGF
  character(*), parameter :: Hardware = 'Linux/PGF'
#elif defined __GNUC__
  character(*), parameter :: Hardware = 'Linux/GFortran'
#else
  character(*), parameter :: Hardware = 'Linux/any'
#endif
#elif ARCH == 3
  character(*), parameter :: Hardware = 'NEC SX-8'
#elif ARCH == 4
  character(*), parameter :: Hardware = 'IBM p690'
#else
  character(*), parameter :: Hardware = 'generic platform'
#endif

  ! Extension of configuration file.
  character(*), parameter :: ConfigFileExtension = '.cfg'

  ! Extension of parameter file.
  character(*), parameter :: ParameterFileExtension = '.par'

  ! Extension of log file.
  character(*), parameter :: LogFileExtension = '.log'

  ! Extension of result file.
  character(*), parameter :: ResultFileExtension = '.run'

  ! Extension of running average result file.
  character(*), parameter :: RunAveFileExtension = '.rav'

  ! Extension of final result file.
  character(*), parameter :: ErrorsFileExtension = '.res'

  ! Extension of visualisation file.
  character(*), parameter :: VisualFileExtension = '.vim'

  ! Extension of normalized potential model file
  character(*), parameter :: NormalizedPotModExtension = '.nrm'

  ! Extension of restart file
  character(*), parameter :: RestartFileExtension = '.rst'

  ! Name tag for output files
  character(FileNameLength) :: OutputNameTag

  ! Parameter file name
  character(FileNameLength) :: ParameterFileName

  ! Restart file name
  character(FileNameLength) :: RestartFileName

  ! Define minimum allowed i/o unit number
#if ARCH == 1
  integer, parameter :: iounit_start  = 7
#else
  integer, parameter :: iounit_start  = 1000
#endif

  ! Define i/o unit numbers
  integer, parameter :: iounit_log     = iounit_start + 0
  integer, parameter :: iounit_config  = iounit_start + 1
  integer, parameter :: iounit_params  = iounit_start + 2
  integer, parameter :: iounit_potmod  = iounit_start + 3
  integer, parameter :: iounit_normal  = iounit_start + 4
  integer, parameter :: iounit_restart = iounit_start + 5
  integer, parameter :: iounit_result  = iounit_start + 6
  integer, parameter :: iounit_runave  = iounit_start + 7
  integer, parameter :: iounit_errors  = iounit_start + 8
  integer, parameter :: iounit_visual  = iounit_start + 9

  ! Define number of output files for each ensemble
  integer, parameter :: FilesPerEnsemble = iounit_visual - iounit_result + 1

  ! Define maximum length of input/output buffer string
  integer, parameter :: IOBufferLength = 1024

  ! Declare input/output buffer strings
  character(IOBufferLength) :: IOBuffer
  character(IOBufferLength) :: ErrorBuffer

  ! Define comment character
  character, parameter :: CommentSign = '#'

  ! Define identifiers used in configuration file
  character(*), parameter :: IdRestart                     = 'Restart'
  character(*), parameter :: IdRestartFileName             = 'RestartFile'
  character(*), parameter :: IdParamsFileName              = 'ParamFile'

  ! Define identifiers used in parameters file
  character(*), parameter :: IdOutputNameTag               = 'OutputNameTag'
  character(*), parameter :: IdUseReducedUnits             = 'Units'
  character(*), parameter :: IdUnitLength                  = 'LengthUnit'
  character(*), parameter :: IdUnitEnergy                  = 'EnergyUnit'
  character(*), parameter :: IdUnitMass                    = 'MassUnit'
  character(*), parameter :: IdEnsembleType                = 'Ensemble'
  character(*), parameter :: IdPistonMass                  = 'PistonMass'
  character(*), parameter :: IdSimulationType              = 'Simulation'
  character(*), parameter :: IdIntegratorType              = 'Integrator'
  character(*), parameter :: IdTimeStep                    = 'TimeStep'
  character(*), parameter :: IdAcceptance                  = 'Acceptance'
  character(*), parameter :: IdNStepsMC                    = 'MCORSteps'
  character(*), parameter :: IdNStepsV                     = 'NVTSteps'
  character(*), parameter :: IdNStepsP                     = 'NPTSteps'
  character(*), parameter :: IdNStepsMue                   = 'mueVTSteps'
  character(*), parameter :: IdNStepsMueP                  = 'muePTSteps'
  character(*), parameter :: IdNSteps                      = 'RunSteps'
  character(*), parameter :: IdBlockSize                   = 'ResultFreq'
  character(*), parameter :: IdErrorsUpdateFrequency       = 'ErrorsFreq'
  character(*), parameter :: IdVisualUpdateFrequency       = 'VisualFreq'
  character(*), parameter :: IdCutoffMode                  = 'CutoffMode'
  character(*), parameter :: IdLongRange                   = 'LongRange'
  character(*), parameter :: IdKappa                       = 'Kappa'
  character(*), parameter :: Idnsqmax                      = 'NsqMax'
  character(*), parameter :: IdNVecMax                     = 'NVecMax'
  character(*), parameter :: IdNMax                        = 'NMax'
  character(*), parameter :: IdGrid                        = 'Grid'
  character(*), parameter :: IdSpline                      = 'Spline'
  character(*), parameter :: IdDebyeLen                    = 'DebyeLen'
  character(*), parameter :: IdNOrient                     = 'NOrient'
  character(*), parameter :: IdRSteps                      = 'RSteps'
  character(*), parameter :: IdMinRadius                   = 'RMinRadius'
  character(*), parameter :: IdMaxRadius                   = 'RMaxRadius'
  character(*), parameter :: IdNEnsembles                  = 'NEnsembles'
  character(*), parameter :: IdRefTemperature              = 'Temperature'
  character(*), parameter :: IdRefPressure                 = 'Pressure'
  character(*), parameter :: IdPressure0                   = 'Pressure0'
  character(*), parameter :: IdLiqDensity                  = 'LiqDensity'
  character(*), parameter :: IdVarLiqDensity               = 'VarDensity'
  character(*), parameter :: IdLiqEnthalpy                 = 'LiqEnthalpy'
  character(*), parameter :: IdVarLiqEnthalpy              = 'VarEnthalpy'
  character(*), parameter :: IdLiqBetaT                    = 'LiqBetaT'
  character(*), parameter :: IdVarLiqBetaT                 = 'VarBetaT'
  character(*), parameter :: IdLiqdHdP                     = 'LiqdHdP'
  character(*), parameter :: IdVarLiqdHdP                  = 'VardHdP'
  character(*), parameter :: IdRefDensity                  = 'Density'
  character(*), parameter :: IdNPart                       = 'NParticles'
  character(*), parameter :: IdNComponents                 = 'NComponents'
  character(*), parameter :: IdPotModFileName              = 'PotModel'
  character(*), parameter :: IdFraction                    = 'MolarFract'
  character(*), parameter :: IdChemPotMethod               = 'ChemPotMethod'
  character(*), parameter :: IdWeightFactors               = 'WeightFactors'
  character(*), parameter :: IdNTest                       = 'NTest'
  character(*), parameter :: IdLiqFraction                 = 'LiqMolarFract'
  character(*), parameter :: IdChemPot                     = 'ChemPot'
  character(*), parameter :: IdVarChemPot                  = 'VarChemPot'
  character(*), parameter :: IdPartialMolarVolume          = 'PartMolVol'
  character(*), parameter :: IdVarPartialMolarVolume       = 'VarPartMolVol'
  character(*), parameter :: IdScaleSigma                  = 'eta'
  character(*), parameter :: IdScaleEpsilon                = 'xi'
  character(*), parameter :: IdRCutoffCOM                  = 'Cutoff'
  character(*), parameter :: IdRCutoffLJ126LJ126           = 'CutoffLJ'
  character(*), parameter :: IdRCutoffDipoleDipole         = 'CutoffDD'
  character(*), parameter :: IdRCutoffDipoleQuadrupole     = 'CutoffDQ'
  character(*), parameter :: IdRCutoffQuadrupoleQuadrupole = 'CutoffQQ'
  character(*), parameter :: IdRFEpsilon                   = 'Epsilon'
  character(*), parameter :: IdFluctFreq                   = 'FluctFreq'
  character(*), parameter :: IdNFullFluct                  = 'NFullFluct'
  character(*), parameter :: IdMaxCounter                  = 'MaxCounter'
  character(*), parameter :: IdUseIntDegFreed              = 'IntDegFreed'   ! switch on internal degrees of freedom
  character(*), parameter :: IdShake                       = 'Shake'         ! QShake algorithm in the case of Leap-frog
  character(*), parameter :: IdIntraLJEl                   = 'IntraLJ_El'    ! all intramolecular 1-5 electrostatic & LJ interaction incl/excl
  character(*), parameter :: IdLJEl14                      = 'LJ_El_14'      ! all 1-4 pairs must be included/excluded from non-bonded inter.
!  character(*), parameter :: IdScaleEl14                   = 'ScaleEl_14'    ! Scale 1-4 electrostatic interactions by this factor
!  character(*), parameter :: IdScaleLJ14                   = 'ScaleLJ_14'    ! Scale 1-4 LJ interactions by this factor

  ! Define identifiers used in potential model file
  character(*), parameter :: IdSite_ntypes                 = 'NSiteTypes'
  character(*), parameter :: IdIdf_ntypes                  = 'NIdfTypes' !number of types for internel degrees of freedom
  character(*), parameter :: IdSite_stype                  = 'SiteType'
  character(*), parameter :: IdIdf_stype                   = 'IdfType'   ! type of internal degree of freedom: bond, angle, dihedral angle
  character(*), parameter :: IdSite_NLJ126                 = 'NSites'
  character(*), parameter :: IdSite_NCharge                = 'NSites'
  character(*), parameter :: IdSite_NDipole                = 'NSites'
  character(*), parameter :: IdSite_NQuadrupole            = 'NSites'
  character(*), parameter :: IdIdf_NBond                   = 'NIdfs'
  character(*), parameter :: IdIdf_NAngle                  = 'NIdfs'
  character(*), parameter :: IdIdf_NDihedral               = 'NIdfs'
  character(*), parameter :: IdUnit_NConstraint            = 'NConstraintUnits'
  character(*), parameter :: IdSite_NDFRot                 = 'NRotAxes'
  character(*), parameter :: IdSite_Mass                   = 'TotalMass'
  character(*), parameter :: IdSite_MOI1                   = 'InertMomX'
  character(*), parameter :: IdSite_MOI2                   = 'InertMomY'
  character(*), parameter :: IdSite_MOI3                   = 'InertMomZ'
  character(*), parameter :: IdLJ126_SiteId                = 'SiteID'
  character(*), parameter :: IdLJ126_r1                    = 'x'
  character(*), parameter :: IdLJ126_r2                    = 'y'
  character(*), parameter :: IdLJ126_r3                    = 'z'
  character(*), parameter :: IdLJ126_sig                   = 'sigma'
  character(*), parameter :: IdLJ126_eps                   = 'epsilon'
  character(*), parameter :: IdLJ126_mass                  = 'mass'
  character(*), parameter :: IdCharge_SiteId               = 'SiteID'
  character(*), parameter :: IdCharge_r1                   = 'x'
  character(*), parameter :: IdCharge_r2                   = 'y'
  character(*), parameter :: IdCharge_r3                   = 'z'
  character(*), parameter :: IdCharge_e                    = 'charge'
  character(*), parameter :: IdCharge_mass                 = 'mass'
  character(*), parameter :: IdCharge_shield               = 'shielding'
  character(*), parameter :: IdDipole_SiteId               = 'SiteID'
  character(*), parameter :: IdDipole_r1                   = 'x'
  character(*), parameter :: IdDipole_r2                   = 'y'
  character(*), parameter :: IdDipole_r3                   = 'z'
  character(*), parameter :: IdDipole_theta                = 'theta'
  character(*), parameter :: IdDipole_phi                  = 'phi'
  character(*), parameter :: IdDipole_D                    = 'dipole'
  character(*), parameter :: IdDipole_mass                 = 'mass'
  character(*), parameter :: IdDipole_shield               = 'shielding'
  character(*), parameter :: IdQuadrupole_SiteId           = 'SiteID'
  character(*), parameter :: IdQuadrupole_r1               = 'x'
  character(*), parameter :: IdQuadrupole_r2               = 'y'
  character(*), parameter :: IdQuadrupole_r3               = 'z'
  character(*), parameter :: IdQuadrupole_theta            = 'theta'
  character(*), parameter :: IdQuadrupole_phi              = 'phi'
  character(*), parameter :: IdQuadrupole_Q                = 'quadrupole'
  character(*), parameter :: IdQuadrupole_mass             = 'mass'
  character(*), parameter :: IdQuadrupole_shield           = 'shielding'
! Inner Degrees of Freedom
  character(*), parameter :: IdConstraint_NSites           = 'Constraint'
  character(*), parameter :: IdConstraint_SiteIds          = 'SiteIds'
  character(*), parameter :: IdConstraint_NDFRot           = 'NRotAxes'
  character(*), parameter :: IdConstraint_Mass             = 'TotalMass'
  character(*), parameter :: IdConstraint_MOI1             = 'InertMomX'
  character(*), parameter :: IdConstraint_MOI2             = 'InertMomY'
  character(*), parameter :: IdConstraint_MOI3             = 'InertMomZ'
  character(*), parameter :: IdBond_Sites                  = 'Bond'
  character(*), parameter :: IdAngle_Sites                 = 'Angle'
  character(*), parameter :: IdDihedral_Sites              = 'Dihedral'
  character(*), parameter :: IdBond_ForConst               = 'ForConst'
  character(*), parameter :: IdBond_R0                     = 'R0'
  character(*), parameter :: IdAngle_ForConst              = 'ForConst'
  character(*), parameter :: IdAngle_Angle0                = 'Angle0'
  character(*), parameter :: IdDihedral_PotBarrier         = 'PotBarrier'
  character(*), parameter :: IdDihedral_n                  = 'n'
  character(*), parameter :: IdDihedral_gamma              = 'gamma'
  character(*), parameter :: IdDihedral_ScaleLJ14          = 'ScaleLJ14'  ! Scale 1-4 LJ interactions by this factor
  character(*), parameter :: IdDihedral_ScaleEl14          = 'ScaleEl14' ! Scale 1-4 electrostatic interactions by this factor
! Chemical Potential
  character(*), parameter :: IdNFluct                      = 'NFluct'

  ! (Almost) zero for mass of inertia
  real(RK), parameter :: Zero = 1E-10_RK

  ! Contribution limit in chemical potential by Widom
  real(RK), parameter :: ContributionLimit = 708.3964185

  ! General mathematical constants
  real(RK), parameter :: Pi = 3.141592689793238462643_RK
  real(RK), parameter :: Pi2 = Pi * 2._RK
  real(RK), parameter :: Pi23 = Pi * 2._RK / 3._RK
  real(RK), parameter :: Pi8 = Pi * 8._RK
  real(RK), parameter :: Pi89 = Pi * 8._RK / 9._RK
  real(RK), parameter :: Pi329 = Pi * 32._RK / 9._RK
  real(RK), parameter :: Piminus2 = Pi * (-2._RK)
  real(RK), parameter :: Piminus83 = Pi * (-8._RK) / 3._RK
  real(RK), parameter :: Third = 1._RK / 3._RK
  real(RK), parameter :: FourThird = 4._RK / 3._RK
  real(RK), parameter :: FiveThird = 5._RK / 3._RK

  ! General physical constants
  real(RK), parameter :: NAvogadro = 6.022137E23_RK
  real(RK), parameter :: kBoltzmann = 1.380658E-23_RK
  real(RK), parameter :: ElementaryCharge = 1.602177E-19_RK
  real(RK), parameter :: VacuumPermittivity = 8.854188E-12_RK

  ! Some physical units
  real(RK), parameter :: Angstroem = 1E-10_RK
  real(RK), parameter :: DegreesInRadian = 180._RK / Pi
  real(RK)            :: DebyesInSI
  real(RK)            :: BuckinghamsInSI

  ! Use reduced units for temperature, pressure, density
  logical :: UseReducedUnits

  ! Use integral degrees of freedom
  logical :: UseIntDegFreed

  ! Shake tolerance for constraint bonds in the case of IDF
  real(RK) :: Shake
  integer, parameter :: Maxit = 100 ! Maximum allowed iterations for QShake

  ! Include intramolecular 1-5 nonbonded interactions
  logical :: IntraLJEl

  ! Include intramolecular 1-4 nonbonded interactions
  logical :: LJEl14

  ! Scaling factor for electrostatic terms in intramolecular 1,4 interactions
!  real(RK) :: ScaleEl14

  ! Scaling factor for Lennard-Jones terms in intramolecular 1,4 interactions
!  real(RK) :: ScaleLJ14

  ! Basic reduced units
  real(RK) :: UnitLength
  real(RK) :: UnitEnergy
  real(RK) :: UnitMass

  ! Derived reduced units
  real(RK) :: UnitVolume
  real(RK) :: UnitTemperature
  real(RK) :: UnitDensity
  real(RK) :: UnitTime
  real(RK) :: UnitForce
  real(RK) :: UnitTorque
  real(RK) :: UnitPressure
  real(RK) :: UnitInertia
  real(RK) :: UnitCharge
  real(RK) :: UnitDipole
  real(RK) :: UnitQuadrupole

  ! Gear corrector coefficients
  real(RK), parameter :: Gear20 = 3._RK / 20._RK
  real(RK), parameter :: Gear21 = 251._RK / 360._RK
  real(RK), parameter :: Gear23 = 11._RK / 18._RK
  real(RK), parameter :: Gear24 = 1._RK / 6._RK
  real(RK), parameter :: Gear25 = 1._RK / 60._RK
  real(RK), parameter :: Gear10 = 251._RK / 720._RK
  real(RK), parameter :: Gear12 = 11._RK / 12._RK
  real(RK), parameter :: Gear13 = 1._RK / 3._RK
  real(RK), parameter :: Gear14 = 1._RK / 24._RK

  ! Describe initial unit cell
  integer, parameter :: NPartInCell = 4
  real(RK), parameter :: &
&   CellX(NPartInCell) = (/.25_RK, .25_RK, .75_RK, .75_RK/), &
&   CellY(NPartInCell) = (/.25_RK, .75_RK, .25_RK, .75_RK/), &
&   CellZ(NPartInCell) = (/.25_RK, .75_RK, .75_RK, .25_RK/)

  ! Type of simulation
  character(80)      :: SimulationTypeString
  integer, parameter :: MolecularDynamics = 1
  integer, parameter :: MonteCarlo        = 2
  integer, parameter :: SecondVirialCoeff = 3
  integer, parameter :: Gibbs             = 4
  integer            :: SimulationType

  ! Type of integrator
  character(80)      :: IntegratorTypeString
  integer, parameter :: IntegratorTypeGear     = 1
  integer, parameter :: IntegratorTypeLeapFrog = 2
  integer, parameter :: IntegratorTypeVerlet   = 3
  integer, parameter :: IntegratorTypeVV       = 4
  integer            :: IntegratorType

  ! Type of ensembles
  character(80)      :: EnsembleTypeString
  integer, parameter :: EnsembleTypeNVE = 1
  integer, parameter :: EnsembleTypeNVT = 2
  integer, parameter :: EnsembleTypeNPH = 3
  integer, parameter :: EnsembleTypeNPT = 4
  integer, parameter :: EnsembleTypeGE  = 5                ! Grand Equilibrium muVT
  integer, parameter :: EnsembleTypeHA  = 6                ! Humid Air mupT
  integer            :: EnsembleType
  logical            :: ConstantTemperature, ConstantPressure

  ! Type of cutoff
  character(80)      :: CutoffModeString
  integer, parameter :: SiteSite     = 1
  integer, parameter :: CenterofMass = 2
  integer            :: CutoffMode
!   integer            :: 

  ! Type of LongRange Mode
  character(80)      :: LongRangeString
  integer, parameter :: Ewald         = 1
  integer, parameter :: RField        = 2
  integer, parameter :: PME           = 3
  integer, parameter :: extRField     = 4
  integer            :: LongRange

  ! Type of method for chemical potential
  integer, parameter :: ChemPotMethodNone    = 0
  integer, parameter :: ChemPotMethodWidom   = 1
  integer, parameter :: ChemPotMethodGradIns = 2

  ! Type of method for weighting factors
  integer, parameter :: WFMethodNone   = 0
  integer, parameter :: WFMethodAuto   = 1
  integer, parameter :: WFMethodGuess  = 2
  integer, parameter :: WFMethodOptSet = 3

  ! MD time step
  real(RK) :: TimeStep, TimeStep2
  real(RK) :: TimeStepSquared, TimeStepSquared2, TimeStepSquaredInv2

  ! MC acceptance rate
  real(RK) :: Acceptance
  real(RK) :: AccUpperLimit,AccLowerLimit

  ! MC acceptance rate auto-adjustment parameters
  real(RK), parameter :: DispTranStart = 0.020_RK
  real(RK), parameter :: DispTranLimit = 0.150_RK
  real(RK), parameter :: DispRotStart  = 0.050_RK
  real(RK), parameter :: DispRotLimit  = 0.150_RK
  real(RK), parameter :: DispVolStart  = 0.010_RK
  real(RK), parameter :: DispVolLimit  = 0.100_RK
  real(RK), parameter :: DispMolTranStart = 0.020_RK
  real(RK), parameter :: DispMolTranLimit = 0.150_RK
  real(RK), parameter :: DispMolRotStart  = 0.050_RK
  real(RK), parameter :: DispMolRotLimit  = 0.150_RK

  ! Frequency of updating MC displacements
  integer, parameter :: DispUpdateFrequency = 100

  ! Number of simulation time steps
  integer :: NSteps

  ! Number of MC overlap reduction steps
  integer :: NStepsMC

  ! Number of NVT equilibration time steps
  integer :: NStepsV

  ! Number of NPT equilibration time steps
  integer :: NStepsP

  ! Number of orientations for second virial coefficient
  integer :: NOrient

  ! Radii for second virial coefficient
  real(RK) :: MinRadius, MaxRadius

  ! Number of current time step
  integer :: Step, StepTotal

  ! Equilibration flags
  logical :: Equilibration, NVTEquilibration, MCOverlapReduction

  ! Restart flag
  logical :: Restart

  ! Too many particles flag (in GE runs)
  logical :: tooManyParticles

  ! Parameters of gradual insertion
  integer :: GradInsFrequency, NFullFluct, MaxCounter

  ! Maximum number of blocks
  integer :: NBlocksMax

  ! Frequency of updating result file
  integer :: BlockSize

  ! Maximum number of block sizes for error calculation
  integer :: NBlockSizesMax

  ! Number of block sizes for error calculation
  integer :: NBlockSizes

  ! Current number of blocks
  integer :: NBlocks

  ! Frequency of updating final result file
  integer :: ErrorsUpdateFrequency

  ! Frequency of updating visualisation file
  integer :: VisualUpdateFrequency

  ! Frequency of updating log file
  integer, parameter :: LogUpdateFrequency = 1000

  ! Internal variables of random number generator
  integer, parameter :: K4B = selected_int_kind(9)
  integer(K4B)       :: ix, iy, tpix
  real(RK)           :: am

  ! Internal variable of FileReadParameter
  integer :: FileReadParameter_LineNumber = 0

  ! MPI variables
#if MPI_VER > 0
  integer :: ierror
  integer :: NProcs
  integer :: NProc
  integer :: NRootProc
  logical :: RootProc
#else
  integer, parameter :: NProcs    = 1
  integer, parameter :: NProc     = 0
  integer, parameter :: NRootProc = NProc
  logical, parameter :: RootProc  = .true.
#endif

#if ARCH == 1 || ARCH == 2 || ARCH == 3
  ! Flag for catched terminate signal
  logical :: TerminateProgram

! PGF compiler version < 6.0 seems to need this
! #ifdef _PGF
!   ! External funtion for signal handling
!   external SetTerminateProgram
! #endif

#else
  logical, parameter :: TerminateProgram = .false.
#endif



!==============================================================!
!  Global procedure interfaces                                 !
!==============================================================!

  interface InitializeProgram
    module procedure Global_InitializeProgram
  end interface

  interface FinalizeProgram
    module procedure Global_FinalizeProgram
  end interface

  interface Warning
    module procedure Global_Warning
  end interface

  interface Error
    module procedure Global_Error
  end interface

  interface AllocationError
    module procedure Global_AllocationError
  end interface

  interface LogOpen
    module procedure Global_LogOpen
  end interface

  interface LogClose
    module procedure Global_LogClose
  end interface

  interface LogWrite
    module procedure Global_LogWrite
  end interface

  interface LogWriteNoAdvance
    module procedure Global_LogWriteNoAdvance
  end interface

  interface LogWriteBlank
    module procedure Global_LogWriteBlank
  end interface

  interface LogWriteTime
    module procedure Global_LogWriteTime
  end interface

  interface LogWriteStep
    module procedure Global_LogWriteStep
  end interface

  interface FileReset
    module procedure Global_FileReset
  end interface

  interface FileRewind
   module procedure Global_FileRewind
  end interface

  interface FileRewrite
    module procedure Global_FileRewrite
  end interface

  interface FileAppend
    module procedure Global_FileAppend
  end interface

  interface FileClose
    module procedure Global_FileClose
  end interface

  interface FileWrite
    module procedure Global_FileWrite
  end interface

  interface FileWriteNoAdvance
    module procedure Global_FileWriteNoAdvance
  end interface

  interface FileWriteBlank
    module procedure Global_FileWriteBlank
  end interface

  ! backward compatible version of FileReadParameter
  interface FileReadParameter_IOBuffer
    module procedure Global_FileReadParameter_buffer
  end interface

  interface FileReadParameter
    module procedure Global_FileReadParameter_String, &
&                    Global_FileReadParameter_Int, &
&                    Global_FileReadParameter_RK, &
&                    Global_FileReadParameter_RKdim1!,
!&                    Global_FileReadParameter_buffer
!                    ambiguous for SX compiler, collision of _buffer with _Int
  end interface

  interface FileWriteParameter
    module procedure Global_FileWriteParameter
  end interface

  interface Randomize
    module procedure Global_Randomize
  end interface

  interface rnd
    module procedure Global_Irnd
  end interface

  interface rnd
    module procedure Global_Rrnd
  end interface

  interface ProcRange
    module procedure Global_GetProcRange
  end interface



!==============================================================!
!  External (intrinsic) procedures                             !
!==============================================================!

#if !( defined _WIN32 || defined __GNUC__ )

  ! Command line arguments
#if ARCH == 1 || ARCH == 2 || ARCH == 3
  integer, external :: iargc
  external getarg
#endif

#ifndef __INTEL_COMPILER

  ! Flush of I/O units
#if ARCH == 1 || ARCH == 2 || ARCH == 3
  external flush
#endif

  ! User name from console
#if ARCH == 1 || defined _PGF
  character(256), external :: getlog
#elif ARCH == 2 || ARCH==3
  external getlog
#endif

  ! Signal handler
#if ARCH == 1 || ARCH == 2 || ARCH == 3
  integer, external :: signal
#endif

  ! Host name
#if ARCH == 1
  external getenv
#elif ARCH == 2 || ARCH == 3
  integer, external :: hostnm
#endif

#endif

#endif



contains



!==============================================================!
!  Subroutine Global_InitializeProgram                         !
!==============================================================!

  subroutine Global_InitializeProgram

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare local variables
#if ARCH == 1 || ARCH == 2 || ARCH == 3
    integer                   :: narg, dot, stat, i
    character(IOBufferLength) :: buffer
#endif

    ! Initialize MPI
#if MPI_VER > 0
    call MPI_Init( ierror )
    call MPI_Comm_size( MPI_COMM_WORLD, NProcs, ierror )
    call MPI_Comm_rank( MPI_COMM_WORLD, NProc, ierror )
    NRootProc = 0
    RootProc = NProc == NRootProc
#endif

!DEBUG
!   if( NProcs > 1 ) then
!     write(IOBuffer, '("debug.out.", I0)') NProc
!     open(999, file=trim(IOBuffer), status='REPLACE', action='WRITE')
!   else
!     open(999, file='debug.out', status='REPLACE', action='WRITE')
!   endif
!DEBUG

    ! Initialize restart flag
    Restart = .false.

    ! Read command line parameter
#if ARCH == 1 || ARCH == 2 || ARCH == 3
    if( RootProc ) then
      call getarg( 0, buffer )
#ifdef _WIN32
      i = index( buffer, '\', BACK=.true. )
#else
      i = index( buffer, '/', BACK=.true. )
#endif
      if( i > 0 ) then
        ProgramFileName = trim( buffer( i+1:len( buffer ) ) )
      else
        ProgramFileName = trim( buffer )         ! possible truncation?
      end if
      narg = iargc()
      !narg = command_argument_count()
      if( narg .lt. 1 ) then
        print *, trim( ProgramFileName ), ' Version: ', VersionString
        print *, 'usage: ', trim( ProgramFileName ), &
&         ' {<par-file[', ParameterFileExtension, ']|<rst-file>', &
&         RestartFileExtension, '}'

        ! Abort program
#if MPI_VER > 0
        call MPI_Abort( MPI_COMM_WORLD, 1, ierror )
#endif
        stop
      end if
      call getarg( 1, buffer )
      !call get_command_argument( 1, buffer )
      buffer = trim( buffer )
      dot = index( buffer, '.', BACK=.true. )
      if( dot > 0 ) then
        if( buffer( dot:len( buffer ) ) .eq. RestartFileExtension ) then
          Restart = .true.
          RestartFileName = trim( buffer )       ! possible truncation

          ! Open restart file for reading
          open( iounit_restart , file = RestartFileName, action = 'READ', &
&           status = 'OLD', iostat = stat )
          if( stat /= 0 ) then
            print *, 'Cannot open restart file ', trim( RestartFileName ), &
&             ' for reading'

            ! Abort program
#if MPI_VER > 0
            call MPI_Abort( MPI_COMM_WORLD, 1, ierror )
#endif
            stop
          end if

          ! Read parameter file name from restart file
          read( iounit_restart, '(A128)' ) ParameterFileName
          buffer = ParameterFileName( 1:len( trim( ParameterFileName ) ) - 4 )

        else if( buffer( dot:len( buffer ) ) .eq. ParameterFileExtension ) then
          buffer = buffer( 1:dot - 1 )
        end if
      end if
      OutputNameTag = trim( buffer )             ! possible truncation
    end if

#if MPI_VER > 0
    call MPI_Bcast( Restart, 1, MPI_LOGICAL, NRootProc, MPI_COMM_WORLD, ierror )
#endif
#endif

    ! Open log file
    call LogOpen

    ! Update log file
#if MPI_VER > 0
    write( IOBuffer, '(I4, " MPI processes initialized")' ) NProcs
#else
    write( IOBuffer, '("sequential Version")' )
#endif
    call LogWrite
    call LogWriteBlank

    ! Set signal handler
#if ARCH == 1 || ARCH == 2
#ifdef __GNUC__
    call signal( 1, IgnoreSignal )
    call signal( 2, SetTerminateProgram )
    call signal( 15, SetTerminateProgram )
#else
    i = signal( 1, SetTerminateProgram, 1 ) ! Ignore SIGHUP
    i = signal( 2, SetTerminateProgram, -1 ) ! Catch SIGINT
    i = signal( 15, SetTerminateProgram, -1 ) ! Catch SIGTERM
#endif
#elif ARCH == 3
    i = signal( 15, SetTerminateProgram )
#endif
#if ARCH == 1 || ARCH == 2 || ARCH == 3
    if( i < 0 ) then
      call Warning('Cannot set signal handler')
    else
      write( IOBuffer, '("Signal handler set successfully")' )
      call LogWrite
      call LogWriteBlank
    end if
#endif

    ! Initialize random number generator
    call Randomize( seed = 5333 )

    ! Define some constants
#ifdef SINGLEPRECISION
    DebyesInSI = real( sqrt( 1E49_8 / (4._8 * real(Pi, 8) &
&     * real(VacuumPermittivity, 8) ) ), RK )
    BuckinghamsInSI = real( sqrt( 1E69_8 / (4._8 * real(Pi, 8) &
&     * real(VacuumPermittivity, 8) ) ), RK )
#else
    DebyesInSI = sqrt( 1E49_8 / (4._8 * Pi * VacuumPermittivity) )
    BuckinghamsInSI = sqrt( 1E69_8 / (4._8 * Pi * VacuumPermittivity) )
#endif

  end subroutine Global_InitializeProgram



!==============================================================!
!  Subroutine Global_FinalizeProgram                           !
!==============================================================!

  subroutine Global_FinalizeProgram

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Close log file
    call LogClose

!DEBUG
!   close(999)
!DEBUG

    ! Finalize MPI
#if MPI_VER > 0
    call MPI_Barrier( MPI_COMM_WORLD, ierror )
    call MPI_Finalize( ierror )
#endif

  end subroutine Global_FinalizeProgram



!==============================================================!
!  Subroutine Global_Message                                   !
!==============================================================!
!
!  subroutine Global_Message( MessageString )
!
!    implicit none
!
!    ! Declare arguments
!    character(*), intent(in), optional :: MessageString
!
!    if( present( MessageString ) ) then
!      IOBuffer = 'MESSAGE: '// trim(MessageString)
!    else
!      IOBuffer = 'MESSAGE: '// trim(MessageBuffer)        ! possible truncation
!    end if
!    if( RootProc ) print *, trim( IOBuffer )
!    call LogWrite
!
!  end subroutine Global_MESSAGE



!==============================================================!
!  Subroutine Global_Warning                                   !
!==============================================================!

  subroutine Global_Warning( ErrorString )

    implicit none

    ! Declare arguments
    character(*), intent(in), optional :: ErrorString

    ! Issue warning
    if( present( ErrorString ) ) then
      IOBuffer = 'WARNING: '// trim(ErrorString)
    else
      IOBuffer = 'WARNING: '// trim(ErrorBuffer)           ! possible truncation
    end if
    if( RootProc ) print *, trim( IOBuffer )
    call LogWrite

  end subroutine Global_Warning



!==============================================================!
!  Subroutine Global_Error                                     !
!==============================================================!

  subroutine Global_Error( ErrorString )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    character(*), intent(in), optional :: ErrorString

    ! Output error message
    call LogWriteBlank
    if( present( ErrorString ) ) then
      IOBuffer = 'ERROR: '// trim( ErrorString )
    else
      IOBuffer = 'ERROR: '// trim( ErrorBuffer ) ! possible truncation
    end if
    if( RootProc ) print *, trim( IOBuffer )
    call LogWrite

    ! Close log file
    call LogClose

    ! Abort program
#if MPI_VER > 0
    call MPI_Abort( MPI_COMM_WORLD, 1, ierror )
#endif
    stop

  end subroutine Global_Error



!==============================================================!
!  Subroutine Global_AllocationError                           !
!==============================================================!

  subroutine Global_AllocationError( stat, str, NPart )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    integer, intent(in)           :: stat
    character(*), intent(in)      :: str
    integer, intent(in), optional :: NPart

    ! Declare local variables
#if MPI_VER > 0
    logical :: ok, okAll
#endif

    ! Check for allocation error
#if MPI_VER > 0
    ok = stat == 0
    call MPI_Allreduce( ok, okAll, 1, MPI_LOGICAL, MPI_LAND, &
&     MPI_COMM_WORLD, ierror )
    if( okAll ) return
#else
    if( stat == 0 ) return
#endif

    ! Terminate program
    if( present( NPart ) ) then
      write( ErrorBuffer, &
&       '("Cannot allocate memory for", I11, " ", A)' ) &
&       NPart, str
    else
      write( ErrorBuffer, '("Cannot allocate memory for ", A)' ) str
    end if
    call Error

  end subroutine Global_AllocationError


!==============================================================!
!  Subroutine Global_LogOpen                                   !
!==============================================================!

  subroutine Global_LogOpen()

    implicit none

    ! Declare local variables
#if ARCH == 1 || ARCH == 2 || ARCH == 3
    character(IOBufferLength) :: hostname = 'unknown host'
    character(IOBufferLength) :: username = 'unknown user'
#if ARCH == 2 || ARCH == 3
    integer                   :: i
#endif
#else
    character(*), parameter   :: hostname = 'unknown host'
    character(*), parameter   :: username = 'unknown user'
#endif

    ! Check for root process
    if( .not. RootProc ) return
    ! Get name of host
#if ARCH == 1
    call getenv( 'HOSTNAME', hostname )
#elif ARCH == 2 || ARCH == 3
#if defined _PGF || (defined __GNUC__ && !defined __PATHSCALE__) || defined __SUNPRO_F90 || ARCH == 3
    i = hostnm( hostname )
#else
    i = hostnam( hostname )
#endif
    if( i .ne. 0 ) hostname = 'unknown host'
#endif

#if ARCH == 1 || defined _PGF
    username = getlog()
#elif ARCH == 2 || ARCH == 3
    call getlog( username )
#endif

    ! Open log file
#if ARCH == 1 || ARCH == 2 || ARCH == 3
    call FileRewrite( iounit_log, trim( OutputNameTag )//LogFileExtension )
#else
    call FileRewrite( iounit_log, ProgramFileName//LogFileExtension )
#endif
    write( IOBuffer, &
&     '("Program ", A, " compiled for ", A)' ) &
&     trim( ProgramFileName ), Hardware
    call LogWrite
    write( IOBuffer, '("version ", A)' ) trim( VersionString )
    call LogWrite
    write( IOBuffer, '("started by ", A," on ", A)' ) &
&     trim( username ), &
&     trim( hostname )
    call LogWriteTime
    call LogWriteBlank

  end subroutine Global_LogOpen



!==============================================================!
!  Subroutine Global_LogClose                                  !
!==============================================================!

  subroutine Global_LogClose()

    implicit none

    ! Check for root process
    if( .not. RootProc ) return

    ! Close log file
    call LogWriteBlank
    write( IOBuffer, '("Program terminated")' )
    call LogWriteTime
    call FileClose( iounit_log )

  end subroutine Global_LogClose



!==============================================================!
!  Subroutine Global_LogWrite                                  !
!==============================================================!

  subroutine Global_LogWrite()

    implicit none

    ! Check for root process
    if( .not. RootProc ) return

    ! Write contents of buffer to log file
    call FileWrite( iounit_log )

    ! Update log file
#if ARCH == 1 || ARCH == 2 || ARCH == 3
    call flush( iounit_log )
#endif

  end subroutine Global_LogWrite



!==============================================================!
!  Subroutine Global_LogWriteNoAdvance                         !
!==============================================================!

  subroutine Global_LogWriteNoAdvance()

    implicit none

    ! Check for root process
    if( .not. RootProc ) return

    ! Write contents of buffer to log file
    call FileWriteNoAdvance( iounit_log )

  end subroutine Global_LogWriteNoAdvance



!==============================================================!
!  Subroutine Global_LogWriteBlank                             !
!==============================================================!

  subroutine Global_LogWriteBlank()

    implicit none

    ! Check for root process
    if( .not. RootProc ) return

    ! Write blank line to log file
    call FileWriteBlank( iounit_log )

  end subroutine Global_LogWriteBlank



!==============================================================!
!  Subroutine Global_LogWriteTime                              !
!==============================================================!

  subroutine Global_LogWriteTime()

    implicit none

    ! Declare local variables
    character(8)  :: date_string
    character(10) :: time_string

    ! Check for root process
    if( .not. RootProc ) return

    ! Update log file
    call LogWriteNoAdvance
    call date_and_time( date_string, time_string )
    write( IOBuffer, &
&     '(" on ", A, ".", A, ".", A, " at ", A, ":", A, ":", A)' ) &
&     date_string(7:8), date_string(5:6), date_string(1:4), &
&     time_string(1:2), time_string(3:4), time_string(5:6)
    call LogWrite

  end subroutine Global_LogWriteTime



!==============================================================!
!  Subroutine Global_LogWriteStep                              !
!==============================================================!

  subroutine Global_LogWriteStep()

    implicit none

    ! Check for root process
    if( .not. RootProc ) return

    ! Update log file
    write( IOBuffer, '(I7, " steps completed")' ) Step
    call LogWriteTime

  end subroutine Global_LogWriteStep



!==============================================================!
!  Subroutine Global_FileReset                                 !
!==============================================================!

  subroutine Global_FileReset( iounit, filename )

    implicit none

    ! Declare arguments
    integer, intent(in)      :: iounit
    character(*), intent(in) :: filename

    ! Declare local variables
    integer :: stat

    ! Check for root process
    if( .not. RootProc ) return

    ! Open file for reading
    write( IOBuffer, '("Opening file <", A, "> for reading")' ) &
&     trim( filename )
    call LogWrite
    open( iounit, file = filename, action = 'READ', status = 'OLD', &
&     iostat = stat )
    if( stat /= 0 ) &
&     call Error( 'Cannot open file '//trim( filename )//' for reading' )

  end subroutine Global_FileReset


!==============================================================!
!  Subroutine Global_FileRewind                                !
!==============================================================!

  subroutine Global_FileRewind( iounit, filename )

    implicit none

    ! Declare arguments
    integer, intent(in)      :: iounit
    character(*), intent(in) :: filename

    ! Declare local variables
    integer :: stat

    ! Check for root process
    if( .not. RootProc ) return

    ! Open file for reading
    write( IOBuffer, '("Opening file <", A, "> for reading")' ) &
&     trim( filename )
    call LogWrite
    rewind( iounit, iostat = stat )
    if( stat /= 0 ) &
&     call Error( 'Cannot rewind file '//trim( filename )//' for reading' )

  end subroutine Global_FileRewind




!==============================================================!
!  Subroutine Global_FileRewrite                               !
!==============================================================!

  subroutine Global_FileRewrite( iounit, filename )

    implicit none

    ! Declare arguments
    integer, intent(in)           :: iounit
    character(*), intent(in)      :: filename

    ! Check for root process
    if( .not. RootProc ) return

    ! Open file for writing
    if( iounit /= iounit_log ) then
      write( IOBuffer, '("Opening file <", A, "> for writing")' ) &
&       trim( filename )
      call LogWrite
    end if
    open( iounit, file = filename, action = 'WRITE', status = 'REPLACE' )

  end subroutine Global_FileRewrite



!==============================================================!
!  Subroutine Global_FileAppend                                !
!==============================================================!

  subroutine Global_FileAppend( iounit, filename )

    implicit none

    ! Declare arguments
    integer, intent(in)           :: iounit
    character(*), intent(in)      :: filename

    ! Declare local variables
    logical :: ex

    ! Check for root process
    if( .not. RootProc ) return

    ! Open file for writing
    if( iounit /= iounit_log ) then
      write( IOBuffer, '("Opening file <", A, "> for appending")' ) &
&       trim( filename )
      call LogWrite
    end if
    inquire( file = filename, exist = ex )
    if( ex ) then
      open( iounit, file = filename, action = 'WRITE', status = 'OLD', &
&       position = 'APPEND' )
    else
      write( IOBuffer, '("File does not exist. Creating new")' )
      call LogWrite
      open( iounit, file = filename, action = 'WRITE', status = 'REPLACE' )
    end if

  end subroutine Global_FileAppend



!==============================================================!
!  Subroutine Global_FileClose                                 !
!==============================================================!

  subroutine Global_FileClose( iounit )

    implicit none

    ! Declare arguments
    integer, intent(in) :: iounit

    ! Declare local variables
    character(FileNameLength) :: fn
#ifdef _WIN32
    integer :: i
#endif

    ! Check for root process
    if( .not. RootProc ) return

    ! Close file
    inquire( iounit, NAME = fn )
#ifdef _WIN32
    i = index( fn, '\', BACK=.true. )
    if( i > 0 ) fn = fn( i+1:len( fn ) )
#endif
    close( iounit )
    if( iounit /= iounit_log ) then
      write( IOBuffer, '("File <", A, "> closed")' ) trim( fn )
      call LogWrite
    end if

  end subroutine Global_FileClose



!==============================================================!
!  Subroutine Global_FileWrite                                 !
!==============================================================!

  subroutine Global_FileWrite( iounit )

    implicit none

    ! Declare arguments
    integer, intent(in) :: iounit

    ! Check for root process
    if( .not. RootProc ) return

    ! Write contents of buffer to file
    call FileWriteNoAdvance( iounit )
    call FileWriteBlank( iounit )

  end subroutine Global_FileWrite



!==============================================================!
!  Subroutine Global_FileWriteNoAdvance                        !
!==============================================================!

  subroutine Global_FileWriteNoAdvance( iounit )

    implicit none

    ! Declare arguments
    integer, intent(in) :: iounit

    ! Check for root process
    if( .not. RootProc ) return

    ! Write contents of buffer to file
    write( iounit, '(A)', advance = 'NO' ) trim( IOBuffer )

  end subroutine Global_FileWriteNoAdvance



!==============================================================!
!  Subroutine Global_FileWriteBlank                            !
!==============================================================!

  subroutine Global_FileWriteBlank( iounit )

    implicit none

    ! Declare arguments
    integer, intent(in) :: iounit

    ! Check for root process
    if( .not. RootProc ) return

    ! Write blank line to file
    write( iounit, '()' )

  end subroutine Global_FileWriteBlank


!==============================================================!
!  Function Global_FileReadParameter                           !
!==============================================================!

  function Global_FileReadParameter( iounit, parameterqualifier, &
&                                    rewind_before, status ) &
&          result (parametervalue)


    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    integer, intent(in)                :: iounit
    character(*), intent(in)           :: parameterqualifier
    logical, intent(in), optional      :: rewind_before
    integer, intent(out), optional     :: status

    character(IOBufferLength) :: parametervalue


    ! Declare local variables
    integer                   :: stat, comment_pos, linesread, i
    character(FileNameLength) :: fn

    ! determine filename
    inquire( iounit, NAME = fn )
    ! Only RootProc reads parameter from file
    if( RootProc ) then
      ! rewind file, if requested
      if( present(rewind_before) ) then
        if( rewind_before ) then
!          write( IOBuffer, '("(",A,":",I4,") rewind")' ) trim(fn),FileReadParameter_LineNumber; call LogWrite
          rewind( iounit )
          FileReadParameter_LineNumber = 0
        end if
      end if
      linesread = 0
      ! loop to read lines until parameter is found
      do
        read( iounit, '(A)', IOSTAT = stat ) parametervalue
        ! error reading from file?
        if( stat > 0 ) then
          call Error( "ERROR reading file "//trim(fn)// &
&                     " while searching for parameter <"//parameterqualifier//">" )
          !if( present(status) ) status = stat
          !return
        ! end of file reached?
        elseif( stat < 0 ) then
          !call Warning( trim(fn)//": Could not find parameter <"//parameterqualifier//">" )
          parametervalue=""
          if( present(status) ) status = stat
          ! (try to) restore position
          if( present(rewind_before) ) then
            if( rewind_before ) then
!              write( IOBuffer, '("(",A,":",I4,") rewind")' ) trim(fn),FileReadParameter_LineNumber; call LogWrite
              rewind( iounit )
              FileReadParameter_LineNumber = 0
              linesread=0
              exit    !not nice!
            end if
          end if
          ! rewind to the position, where the reading process was started
          backspace( iounit )   ! "undo" last read, where eof was encountered
          do i = 1,linesread
!            write( IOBuffer, '("(",A,":",I4,") backspace")' ) trim(fn),FileReadParameter_LineNumber; call LogWrite
            backspace( iounit )
            FileReadParameter_LineNumber = FileReadParameter_LineNumber - 1
          end do
          exit
        end if
        FileReadParameter_LineNumber = FileReadParameter_LineNumber + 1
        linesread = linesread + 1
!        write( IOBuffer, '("(",A,":",I4,") read:",A)' ) trim(fn),FileReadParameter_LineNumber,trim(parametervalue); call LogWrite
!         check for comment token
        comment_pos = index( parametervalue, CommentSign )
        if( comment_pos > 0 ) then
!          write( IOBuffer, '("(",A,":",I4,") comment:",A)' ) trim(fn),FileReadParameter_LineNumber, &
!&               trim(parametervalue(comment_pos:len(parametervalue))); call LogWrite
          !                eliminate comment part of line
          parametervalue = parametervalue(1:comment_pos - 1)
        end if
        if( index( adjustl( parametervalue ), trim( parameterqualifier ) ) == 1 ) then
          ! extract value part (after =)
          parametervalue = parametervalue( index( parametervalue, '=' )+1:len( parametervalue ) )
          parametervalue = trim( adjustl( parametervalue ) )
          if( present(status) ) status = 0
          exit
        end if
      end do

    end if

    ! Broadcast parameter to other processes
    ! (2 Broadcast are not very efficient, but it doesn't need to be efficient here.
    !  Better broadcast the integer, float parametervalues, instead of the string?)
#if MPI_VER > 0
    call MPI_Bcast( parametervalue, len(parametervalue), &
&     MPI_CHARACTER, NRootProc, MPI_COMM_WORLD, ierror )
    if( present(status) ) then
      call MPI_Bcast( status, 1, &
&       MPI_INTEGER, NRootProc, MPI_COMM_WORLD, ierror )
    end if
#endif

!    write( IOBuffer, '(I5," (",A,":",I4,") String ",A," =",A)' ) NProc,trim(fn),FileReadParameter_LineNumber, &
!&                      trim(parameterqualifier),trim(parametervalue); call LogWrite

  end function Global_FileReadParameter


!==============================================================!
!  Subroutine Global_FileReadParameter_buffer                  !
!==============================================================!

  subroutine Global_FileReadParameter_buffer( iounit, parameterqualifier, &
&                                            rewind_before, defaultvalue, status )
  ! this subroutine is for backward compatibily purposes

    implicit none

    ! Declare arguments
    integer, intent(in)                :: iounit
    character(*), intent(in)           :: parameterqualifier
    logical, intent(in), optional      :: rewind_before
    character(*), intent(in), optional :: defaultvalue
    integer, intent(out), optional     :: status

    ! Declare local variables
    integer                   :: stat

    IOBuffer = Global_FileReadParameter(iounit, parameterqualifier, rewind_before, stat)
    if ( stat < 0 ) then
      if ( present(defaultvalue) ) then
        write( IOBuffer, '("setting ",A," (IOBuffer) to default value ",A)' ) &
&             trim(parameterqualifier), trim(defaultvalue)
        call LogWrite
        IOBuffer = defaultvalue
      else
        call Error( "Could not find parameter <"//parameterqualifier//">" )
        !return
      end if
    end if
    if ( present(status) ) status=stat

    ! Broadcast parameter
!#if MPI_VER > 0
!    call MPI_Bcast( IOBuffer, IOBufferLength, &
!&     MPI_CHARACTER, NRootProc, MPI_COMM_WORLD, ierror )
!#endif

!    inquire( iounit, NAME = fn )
!    write( IOBuffer, '(I5," (",A,":",I4,";",I2,") IOBuffer ",A," =",A)' ) NProc,trim(fn),FileReadParameter_LineNumber, &
!&          stat,trim(parameterqualifier),trim(IOBuffer); call LogWrite

  end subroutine Global_FileReadParameter_buffer


!==============================================================!
!  Subroutine Global_FileReadParameter_String                  !
!==============================================================!

  subroutine Global_FileReadParameter_String( parametervariable, iounit, parameterqualifier, &
&                                            rewind_before, defaultvalue, status )
  ! setting up functions with result (parametervalue) for different data types is ambigious
  ! for a FileReadParameter polymorphism 

    implicit none

    ! Declare arguments
    character(*), intent(out)          :: parametervariable
    integer, intent(in)                :: iounit
    character(*), intent(in)           :: parameterqualifier
    logical, intent(in), optional      :: rewind_before
    character(*), intent(in), optional :: defaultvalue
    integer, intent(out), optional     :: status


    ! Declare local variables
    integer                   :: stat
    !character(FileNameLength) :: fn

    parametervariable = Global_FileReadParameter(iounit, parameterqualifier, rewind_before, stat)
    if ( stat < 0 ) then
      if ( present(defaultvalue) ) then
        write( IOBuffer, '("setting ",A," to default value ",A)' ) &
&             trim(parameterqualifier), trim(defaultvalue)
        call LogWrite
        parametervariable = defaultvalue
      else
        call Error( "Could not find parameter <"//parameterqualifier//">" )
        !return
      end if
    end if
    if ( present(status) ) status=stat

    ! Broadcast parameter to other processes
!#if MPI_VER > 0
!    call MPI_Bcast( parametervariable, len(parametervariable), &
!&     MPI_CHARACTER, NRootProc, MPI_COMM_WORLD, ierror )
!#endif

!    inquire( iounit, NAME = fn )
!    write( IOBuffer, '(I5," (",A,":",I4,";",I2,") String ",A," =",A)' ) NProc,trim(fn),FileReadParameter_LineNumber, &
!&          stat,trim(parameterqualifier),trim(parametervariable); call LogWrite

  end subroutine Global_FileReadParameter_String


!==============================================================!
!  Subroutine Global_FileReadParameter_Int                     !
!==============================================================!

  subroutine Global_FileReadParameter_Int( parametervariable, iounit, parameterqualifier, &
&                                         rewind_before, defaultvalue, status )
  ! Global_FileReadParameter_Integer has 32>31 characters!

    implicit none

    ! Declare arguments
    integer, intent(out)           :: parametervariable
    integer, intent(in)            :: iounit
    character(*), intent(in)       :: parameterqualifier
    logical, intent(in), optional  :: rewind_before
    integer, intent(in), optional  :: defaultvalue
    integer, intent(out), optional :: status

    ! Declare local variables
    character(IOBufferLength) :: buffer
    integer                   :: stat
    !character(FileNameLength) :: fn

    buffer = Global_FileReadParameter(iounit, parameterqualifier, rewind_before, stat)
    if ( stat == 0 ) then
      read( buffer, * ) parametervariable
    else if ( stat < 0 ) then
      if ( present(defaultvalue) ) then
        write( IOBuffer, '("setting ",A," to default value ",I7)' ) &
&             trim(parameterqualifier), defaultvalue
        call LogWrite
        parametervariable = defaultvalue
      else
        call Error( "Could not find parameter <"//parameterqualifier//">" )
        !return
      end if
    end if
    if ( present(status) ) status=stat

    ! Broadcast parameter to other processes
!#if MPI_VER > 0
!    call MPI_Bcast( parametervariable, 1, &
!&     MPI_INTEGER, NRootProc, MPI_COMM_WORLD, ierror )
!#endif

!    inquire( iounit, NAME = fn )
!    write( IOBuffer, '(I5," (",A,":",I4,";",I2,") Integer ",A," =",I7)' ) NProc,trim(fn),FileReadParameter_LineNumber, &
!&          stat,trim(parameterqualifier),parametervariable; call LogWrite

  end subroutine Global_FileReadParameter_Int


!==============================================================!
!  Subroutine Global_FileReadParameter_RK                      !
!==============================================================!

  subroutine Global_FileReadParameter_RK( parametervariable, iounit, parameterqualifier, &
&                                        rewind_before, defaultvalue, status )

    implicit none

    ! Declare arguments
    real(RK), intent(out)          :: parametervariable
    integer, intent(in)            :: iounit
    character(*), intent(in)       :: parameterqualifier
    logical, intent(in), optional  :: rewind_before
    real(RK), intent(in), optional :: defaultvalue
    integer, intent(out), optional :: status


    ! Declare local variables
    character(IOBufferLength) :: buffer
    integer                   :: stat
    !character(FileNameLength) :: fn

    buffer = Global_FileReadParameter(iounit, parameterqualifier, rewind_before, stat)
    if ( stat == 0 ) then
      read( buffer, * ) parametervariable
    else if ( stat < 0 ) then
      if ( present(defaultvalue) ) then
        write( IOBuffer, '("setting ",A," to default value ",G15.9)' ) &
&             trim(parameterqualifier), defaultvalue
        call LogWrite
        parametervariable = defaultvalue
      else
        call Error( "Could not find parameter <"//parameterqualifier//">" )
        !return
      end if
    end if
    if ( present(status) ) status=stat

    ! Broadcast parameter to other processes
!#if MPI_VER > 0
!    call MPI_Bcast( parametervariable, 1, &
!&     MPI_RK, NRootProc, MPI_COMM_WORLD, ierror )
!#endif

!    inquire( iounit, NAME = fn )
!    write( IOBuffer, '(I5," (",A,":",I4,";",I2,") Integer ",A," =",G15.9)' ) NProc,trim(fn),FileReadParameter_LineNumber, &
!&          stat,trim(parameterqualifier),parametervariable; call LogWrite

  end subroutine Global_FileReadParameter_RK


!==============================================================!
!  Subroutine Global_FileReadParameter_RKdim1                  !
!==============================================================!

  subroutine Global_FileReadParameter_RKdim1( parametervariable, iounit, parameterqualifier, &
&                                            rewind_before, defaultvalue, status )

    implicit none

    ! Declare arguments
    real(RK), dimension(:), intent(out)          :: parametervariable
    integer, intent(in)                          :: iounit
    character(*), intent(in)                     :: parameterqualifier
    logical, intent(in), optional                :: rewind_before
    real(RK), dimension(:), intent(in), optional :: defaultvalue
    integer, intent(out), optional               :: status

    ! Declare local variables
    character(IOBufferLength) :: buffer
    integer                   :: stat
    !character(FileNameLength) :: fn

    buffer = Global_FileReadParameter(iounit, parameterqualifier, rewind_before, stat)
    if ( stat == 0 ) then
      read( buffer, * ) parametervariable
    else if ( stat < 0 ) then
      if ( present(defaultvalue) ) then
        write( IOBuffer, '("setting ",A," to default value ")' ) trim(parameterqualifier)
        call LogWrite
        write( IOBuffer, * ) defaultvalue
        call LogWrite
        parametervariable = defaultvalue
      else
        call Error( "Could not find parameter <"//parameterqualifier//">" )
        !return
      end if
    end if
    if ( present(status) ) status=stat

    ! Broadcast parameter to other processes
!#if MPI_VER > 0
!    call MPI_Bcast( parametervariable, size(parametervariable), &
!&     MPI_INTEGER, NRootProc, MPI_COMM_WORLD, ierror )
!#endif

!    inquire( iounit, NAME = fn )
!    write( IOBuffer, '(I5," (",A,":",I4,";",I2,") Real Array ",A," =")' ) NProc,trim(fn),FileReadParameter_LineNumber, &
!&          stat,trim(parameterqualifier); call LogWrite
!    write( IOBuffer, * ) parametervariable; call LogWrite

  end subroutine Global_FileReadParameter_RKdim1





!==============================================================!
!  Subroutine Global_FileWriteParameter                        !
!==============================================================!

  subroutine Global_FileWriteParameter( iounit, parameter )

    implicit none

    ! Declare arguments
    integer, intent(in)      :: iounit
    character(*), intent(in) :: parameter

    ! Check for root process
    if( .not. RootProc ) return

    ! Write parameter to file
    write( iounit, '(A, T17, "=", A)' ) trim( parameter ), trim( IOBuffer )

  end subroutine Global_FileWriteParameter



!==============================================================!
!  Subroutine Global_Randomize                                 !
!==============================================================!

  subroutine Global_Randomize( seed )

    implicit none

    ! Declare arguments
    integer, intent(in) :: seed

    ! Initialize with given seed
    ix = ior(ieor(888889999, seed), 1)
    iy = ieor(777755555, seed)

    ! Initialize test particle random number generator
    tpix = NProc
    tpix = ieor( tpix, ishft(tpix,5) ) + 1422217823
    tpix = ieor( tpix, ishft(tpix,-16) ) + 1842055030
    tpix = ieor( tpix, ishft(tpix,9) ) + 80567781
    tpix = ior( tpix, 1 )

    ! Calculate normalization factor
    am = nearest(1._RK, -1._RK) / huge(ix)

    write( IOBuffer, '("Random number generator initialized")' )
    call LogWrite
    call LogWriteBlank

  end subroutine Global_Randomize



!==============================================================!
!  Function Global_Irnd                                        !
!==============================================================!

  function Global_Irnd( range ) result( iharvest )

    implicit none

    ! Declare arguments
    integer, intent(in) :: range

    ! Declare result
    integer :: iharvest

    ! Declare local variables
    integer(K4B), parameter :: IA=16807, IM=2147483647, &
&     IQ=127773, IR=2836
    integer(K4B), save      :: k

    ! Generate random number
    ix = ieor(ix, ishft(ix, 13))
    ix = ieor(ix, ishft(ix, -17))
    ix = ieor(ix, ishft(ix, 5))
    k = iy / IQ
    iy = IA * (iy - k * IQ) - IR * k
    if( iy < 0 ) iy = iy + IM
    iharvest = 1 + ishft(int(range, 8) * &
&     ior(iand(IM, ieor(ix, iy)), 1), -31)

  end function Global_Irnd



!==============================================================!
!  Function Global_Rrnd                                        !
!==============================================================!

  function Global_Rrnd( l_range, h_range ) result( rharvest )

    implicit none

    ! Declare arguments
    real(RK), intent(in) :: l_range, h_range

    ! Declare result
    real(RK) :: rharvest

    ! Declare local variables
    integer(K4B), parameter :: IA=16807, IM=2147483647, &
&     IQ=127773, IR=2836
    integer(K4B), save      :: k

    ! Generate random number
    ix = ieor(ix, ishft(ix, 13))
    ix = ieor(ix, ishft(ix, -17))
    ix = ieor(ix, ishft(ix, 5))
    k = iy / IQ
    iy = IA * (iy - k * IQ) - IR * k
    if( iy < 0 ) iy = iy + IM
    rharvest = l_range + am * ior(iand(IM,ieor(ix,iy)),1) &
&     * (h_range - l_range)

  end function Global_Rrnd


#if ARCH == 1 || ARCH == 2 || ARCH == 3
!==============================================================!
!  Subroutine SetTerminateProgram                              !
!==============================================================!

#ifdef __INTEL_COMPILER
  function SetTerminateProgram( signum ) result( sigout )
#else
  subroutine SetTerminateProgram
#endif

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

#ifdef __INTEL_COMPILER
    ! Declare arguments
    integer, intent(in) :: signum

    ! Declare result
    integer :: sigout
#endif

    ! Set flag to terminate program
    TerminateProgram = .true.

#ifdef __INTEL_COMPILER
    ! Set return value
    sigout = 0

  end function SetTerminateProgram
#else
  end subroutine SetTerminateProgram
#endif

#endif

#if __GNUC__
!==============================================================!
!  Subroutine IgnoreSignal                                     !
!==============================================================!

  subroutine IgnoreSignal

    implicit none

    continue

  end subroutine IgnoreSignal
#endif


!==============================================================!
!  Subroutine Write Restart File on xc2 in Karlsruhe, Germany  !
!==============================================================!

#if MPI_VER > 0
  subroutine time_left(time_limit)

    implicit none

    ! Include MPI header
    include 'mpif.h'

    real(RK) :: time_remaining
    real(RK) :: cputime,max_cpu_time
    character*10 string_max_time
    integer  :: ierror,max_time, err
    integer  :: time_limit

! Get CPU time consumed by each task and compute the maximum value
    call cpu_time(cputime)
!    call MPI_Allreduce(cputime,max_cpu_time,1,MPI_DOUBLE_PREICSION,MPI_MAX,MPI_COMM_WORLD,ierror)

    max_time = 1e7
#ifdef KARLS
! getenv delivers the value of the environment variable JMS_t
    call getenv('JMS_t',string_max_time)

! Convert to integer
    read(string_max_time,*) max_time
#endif
#ifdef ITWM
! getenv WALLTIME
    call getenv('WALLTIME',string_max_time)
! Convert to integer
    read(string_max_time,*) max_time
#endif

! Compute the remaining CPU time
    time_remaining = max_time - real(cputime)/60.

    if (time_remaining .le. time_limit) then
       write( IOBuffer, '("Simulation Abort due to Time Constraints on simulation cluster")' )
       call LogWrite
       call LogWriteBlank

#ifdef __INTEL_COMPILER
         err = SetTerminateProgram( 1 )
#else
         call SetTerminateProgram
#endif
    end if

  end subroutine time_left
#endif


!==============================================================!
!  Function Global_GetProcRange                                !
!==============================================================!

  function Global_GetProcRange( overall_size, first_index, last_index ) result( range_size )

    implicit none

    ! Declare arguments
    integer, intent(in) :: overall_size
    integer, intent(out) :: first_index, last_index

    ! Declare result
    integer :: range_size
    ! the function could return an array containing the indices, but NPart0..NPart2 are already scalar values.

    ! Declare local variables
    !integer :: range_size0             ! version 1 only: range size for the first process

#if MPI_VER > 0
    ! original version: last process might get smaller range_size
    range_size = 1 + (overall_size - 1) / NProcs
    first_index = 1 + NProc * range_size
    last_index = min( first_index + range_size - 1, overall_size )
    range_size = last_index - first_index + 1

    ! alternative version 1: first process ("master", NProc==0) might get smaller range_size
    !range_size = ceiling( real(overall_size)/NProcs )
    !range_size0 = mod( overall_size, range_size )
    !last_index = range_size0 + NProc*range_size
    !if ( NProc == 0 ) then
    !  range_size = range_size0
    !end if
    !first_index = last_index-range_size+1

    ! alternative version 2: distribute, use round instead of int?
    !first_index = int(real(NProc)/NProcs*overall_size)+1
    !last_index = int(real(NProc+1)/NProcs*overall_size)
    !range_size = last_index - first_index + 1

#else
    first_index=1
    last_index = overall_size
    range_size=overall_size
#endif

  end function Global_GetProcRange



end module ms2_global


