!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 3.0                !
!  (c) 2017 by TU Kaiserslautern / U Paderborn                 !
!      P.O. Box 67653                                          !
!      67653 Kaiserslautern                                    !
!==============================================================!
!  Module ms2_global                                           !
!  Contains declarations of global constants and functions     !
!==============================================================!

!****************************************************************
!* Updates and auxiliary routines are available from            *
!* http://www.ms-2.de                                           *
!****************************************************************

!#define USE_PRINTPROCSTATUS

#ifndef ARCH
#define ARCH    0
#define FORTRAN 90
#define MPI_VER 0
#endif

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2_global.F90...'
#endif

!           __GFORTRAN__
#if defined __GNUC__
! the gfortran preprocessor seems not to support the # operator
#define MACRODEF_STRINGIFY(x) "x"
#else
#define MACRODEF_STRINGIFY(x) #x
#endif
#define MACRODEF_TO_STRING(x)  MACRODEF_STRINGIFY(x)

#if defined(__GNUC__)
# if defined(__GNUC_PATCHLEVEL__)
#  define __GNUC_VERSION__ (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
# else
#  define __GNUC_VERSION__ (__GNUC__ * 10000 + __GNUC_MINOR__ * 100)
# endif
#endif

#ifndef TRANS
#define TRANS 0
#endif

#ifndef OSMOP
#define OSMOP 0
#endif

#ifndef HBOND
#define HBOND 0
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
  integer, parameter :: RK = KIND(1.0)
#else
  integer, parameter :: RK = KIND(1D0)
#endif

#if MPI_VER > 0
  integer :: MPI_RK
#endif

  ! Identifier for MC overlaps
   logical :: MCOverlapDetected

  ! limits
  real(RK)            :: limits_RK_MAX
  real(RK)            :: exp_arg_max  != log(limits_RK_MAX)

  ! Define maximum length of file names
  integer, parameter :: FileNameLength = 128

  ! Name of program
#if ARCH == 1 || ARCH == 2 || ARCH == 3
  character(28)           :: ProgramFileName
#else
  character(*), parameter :: ProgramFileName = 'ms2'
#endif

  ! Version of program
  character(*), parameter :: VersionString = 'v2.0'
  real(RK)                :: ms2VersionNr = 2.0_RK
#ifdef __DATE__
#ifdef __TIME__
  character(*), parameter :: CompileTime = __DATE__ // ',' // __TIME__
#else
  character(*), parameter :: CompileTime = __DATE__
#endif
#else
#ifdef __TIME__
  character(*), parameter :: CompileTime = __TIME__
#else
  character(*), parameter :: CompileTime = 'unknown compile time'
#endif
#endif

  ! Name of platform
#if ARCH == 1
  character(*), parameter :: Hardware = 'alpha'
#elif ARCH == 2
#ifdef _WIN32
  character(*), parameter :: Hardware = 'pc/win32'
#elif defined __INTEL_COMPILER
  character(*), parameter :: Hardware = 'pc/ifort'
#elif defined __SUNPRO_F90
  character(*), parameter :: Hardware = 'pc/sunF90'
#elif defined __PATHSCALE__
  character(*), parameter :: Hardware = 'pc/pathf9X'
#elif defined _PGF || defined __PGI
  character(*), parameter :: Hardware = 'pc/PGI'
#elif defined _CRAYFTN
  character(*), parameter :: Hardware = 'XE6/CRAY'
#elif defined __GNUC__
  character(*), parameter :: Hardware = 'pc/gfortran'
#else
  character(*), parameter :: Hardware = 'pc/any'
#endif
#elif ARCH == 3
  character(*), parameter :: Hardware = 'NEC SX'
#elif ARCH == 4
  character(*), parameter :: Hardware = 'IBM p690'
#else
  character(*), parameter :: Hardware = 'generic platform'
#endif


! define platform-specific path separator
#ifdef _WIN32
  character(*), parameter :: FileSep = '\'
#else
  character(*), parameter :: FileSep = '/'
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
  
  !DC NOTE- Extension of cluster related visualisation file.
  character(*), parameter :: VisualCCFileExtension = '.cvim'
  
  !DC NOTE- Extension of cluster criteria info file.
  character(*), parameter :: CCFileExtension = '.clust'

  !DC NOTE- Extension of cluster criteria grid position file.
  character(*), parameter :: GridFileExtension = '.grid'
  
  ! Extension of visualisation of h-bonding file.
  character(*), parameter :: VisualHBFileExtension = '.hbvim'

  ! Extension of normalized potential model file
  character(*), parameter :: NormalizedPotModExtension = '.nrm'

  ! Extension of restart file
  character(*), parameter :: RestartFileExtension = '.rst'

  ! Extension of RDF file
  character(*), parameter :: RDFFileExtension = '.rdf'
  
  ! Extension of ODF file 
  character(*), parameter :: ODFFileExtension = '.odf'

  ! Extension of KBI file (Kirkwood-Buff Integration)
  character(*), parameter :: KBIrdfFileExtension = '.kbirdf'
  character(*), parameter :: KBIravFileExtension = '.kbirav'

  ! Extension of alpha2 file (displacement correlation function)
  character(*), parameter :: ALPHA2ravFileExtension = '.a2rav'

  !EinsteinCoef data extension
  character(*), parameter :: EinsteinCoefFileExtension = '.ecoef'

  ! Extension of ThermoInt filename
  character(*), parameter :: ThermoIntFileExtension = '.thi'

  ! Extension fo result correlation fucntion
  character(*), parameter :: ResultTransportExtension = '.rtr'

  ! Extension of DCP file
  character(*), parameter :: DCPFileExtension = '.dcp'

  ! Marker within a result file for each ensemble data
  character(*), parameter :: RstEnsembleMarker = 'ENSEMBLE'

  ! Name tag for output files
  character(FileNameLength) :: OutputNameTag
  ! true, if OutputNameTag is set through the command line argument
  ! MPI_VER>0: only set on RootProc
  logical :: OutputNameTagfromCommandline

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
  integer, parameter :: iounit_log       = iounit_start +  0
  integer, parameter :: iounit_config    = iounit_start +  1
  integer, parameter :: iounit_params    = iounit_start +  2
  integer, parameter :: iounit_potmod    = iounit_start +  3
  integer, parameter :: iounit_normal    = iounit_start +  4
  integer, parameter :: iounit_restart   = iounit_start +  5
  integer, parameter :: iounit_result    = iounit_start +  6
  integer, parameter :: iounit_runave    = iounit_start +  7
  integer, parameter :: iounit_errors    = iounit_start +  8
  integer, parameter :: iounit_visual    = iounit_start +  9
  integer, parameter :: iounit_rdf       = iounit_start + 10
  integer, parameter :: iounit_thermoint = iounit_start + 11
  integer, parameter :: iounit_rescf     = iounit_start + 12
  integer, parameter :: iounit_visualHB  = iounit_start + 13
  integer, parameter :: iounit_dcp       = iounit_start + 14
  integer, parameter :: iounit_kbirdf    = iounit_start + 15
  integer, parameter :: iounit_kbirav    = iounit_start + 16
  integer, parameter :: iounit_a2rav     = iounit_start + 17
  integer, parameter :: iounit_proc      = iounit_start + 18
  integer, parameter :: iounit_ecoef     = iounit_start + 19   !EinsteinCoef
  integer, parameter :: iounit_ccpos     = iounit_start + 20 !DC TODO - this should be changed appropriate to the other output files
  integer, parameter :: iounit_cc        = iounit_start + 21 !DC TODO - this should be changed appropriate to the other output files
  integer, parameter :: iounit_ccgrid    = iounit_start + 22 !DC TODO - this should be changed appropriate to the other output files

#if MPI_VER > 0
  integer            :: iounit_result_parallel = iounit_start + 6
  integer            :: iounit_runave_parallel = iounit_start + 7
#endif

  ! Define number of output files for each ensemble
  integer, parameter :: FilesPerEnsemble = iounit_dcp - iounit_result + 1

  ! Define maximum length of input/output buffer string
  integer, parameter :: IOBufferLength = 1024

  ! Declare input/output buffer strings
  character(IOBufferLength) :: IOBuffer
  character(IOBufferLength) :: ErrorBuffer

  ! Define comment character
  character, parameter :: CommentSign = '#'
  ! Define whitespaces                            TAB       CR
  character(*), parameter :: Whitespaces=' '//char(9)//char(13)

  ! Define identifiers used in configuration file
  character(*), parameter :: IdRestart                     = 'Restart'
  character(*), parameter :: IdRestartFileName             = 'RestartFile'
  character(*), parameter :: IdParamsFileName              = 'ParamFile'

  ! Define identifiers used in parameters file
  character(*), parameter :: IdOutputNameTag               = 'OutputNameTag'
  character(*), parameter :: IdparVersionNr                = 'ms2Version'
  character(*), parameter :: IdWallTime                    = 'WallTime'
  character(*), parameter :: IdTimeLimit                   = 'TimeLimit'
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
  character(*), parameter :: IdNStepsE                     = 'NVESteps'
  character(*), parameter :: IdNStepsP                     = 'NPTSteps'
  character(*), parameter :: IdNStepsH                     = 'NPHSteps'
  character(*), parameter :: IdNStepsMue                   = 'mueVTSteps'
  character(*), parameter :: IdNStepsMueP                  = 'muePTSteps'
  character(*), parameter :: IdNSteps                      = 'RunSteps'
  character(*), parameter :: IdBlockSize                   = 'ResultFreq'
  character(*), parameter :: IdErrorsUpdateFrequency       = 'ErrorsFreq'
  character(*), parameter :: IdVisualUpdateFrequency       = 'VisualFreq'
  character(*), parameter :: IdRDFUpdateFrequency          = 'RDFFreq'
  character(*), parameter :: IdRDFNumberShells             = 'NumShells'
  character(*), parameter :: IdnR                          = 'NumShellsODF'
  character(*), parameter :: IdnPhi                        = 'nPhiODF'
  character(*), parameter :: IdnGamma                      = 'nGammaODF'
  character(*), parameter :: IdODFUpdateFrequency          = 'ODFRecordingFreq'
  character(*), parameter :: IdODFOutputFrequency          = 'ODFOutputFreq'
  character(*), parameter :: IdKBIUpdateFrequency          = 'KBIFreq' !Kirkwood-Buff Integration
  character(*), parameter :: IdKBINumberShells             = 'KBINumShells'
  character(*), parameter :: IdKBIResetFrequency           = 'KBIResetFreq'
  character(*), parameter :: IdALPHA2UpdateFrequency       = 'ALPHA2Freq' !Alpha2 correlation function
  character(*), parameter :: IdALPHA2Length                = 'ALPHA2Length'
  character(*), parameter :: IdALPHA2Shift                 = 'ALPHA2Span'
  character(*), parameter :: IdEinsteinCoefCalc            = 'EinsteinCoefCalc' !EinsteinCoef
  character(*), parameter :: IdNBinsDen                    = 'NumDenBins'
  character(*), parameter :: IdWallForce                   = 'Wallforce'
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
  character(*), parameter :: IdmpiEnsembleGroups           = 'mpiEnsembleGroups'
  character(*), parameter :: IdRefTemperature              = 'Temperature'
  character(*), parameter :: IdRefHamiltonian              = 'Hamiltonian'
  character(*), parameter :: IdRefEnthalpy                 = 'Enthalpy'
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
  character(*), parameter :: IdFraction                    = 'MolarFract:MoleFract'
  character(*), parameter :: IdChemPotMethod               = 'ChemPotMethod'
  character(*), parameter :: IdPermeability                = 'Permeability'
  character(*), parameter :: IdNHBonds                     = 'NHBondCriteria'

  !DC NOTE- cluster criteria relevant global Id
  character(*), parameter :: IdIsClusterCriteria           = 'ClusterIsCriteria'  
  character(*), parameter :: IdCCUpdateFrequency           = 'ClusterCriteriaFreq'  
  character(*), parameter :: IdCcrittype                   = 'ClusterCriteriaType'  
  character(*), parameter :: IdCcritdist                   = 'ClusterCriteriaDistance'  
  character(*), parameter :: IdCcount                      = 'ClusterMoleculeCount'  
  character(*), parameter :: IdCmax                        = 'ClusterMaximumAllowed'  
  character(*), parameter :: IdIsCvim                      = 'ClusterIsCvim'  

  !Koester
  character(*), parameter :: IdGradInsInit                 = 'GISteps'
  character(*), parameter :: IdWeightFactors               = 'WeightFactors'
  character(*), parameter :: IdNTest                       = 'NTest'
  character(*), parameter :: IdLiqFraction                 = 'LiqMolarFract:LiqMoleFract'
  character(*), parameter :: IdChemPot                     = 'ChemPot'
  character(*), parameter :: IdVarChemPot                  = 'VarChemPot'
  character(*), parameter :: IdPartialMolarVolume          = 'PartMolVol'
  character(*), parameter :: IdVarPartialMolarVolume       = 'VarPartMolVol'
  character(*), parameter :: IdPartialMolarEnthalpy        = 'PartMolEnt'
  character(*), parameter :: IdScaleSigma                  = 'eta'
  character(*), parameter :: IdScaleEpsilon                = 'xi'
  character(*), parameter :: IdRCutoffCOM                  = 'Cutoff'
  character(*), parameter :: IdRCutoffMIEnmMIEnm           = 'CutoffMIE'
  character(*), parameter :: IdRCutoffDipoleDipole         = 'CutoffDD'
  character(*), parameter :: IdRCutoffDipoleQuadrupole     = 'CutoffDQ'
  character(*), parameter :: IdRCutoffQuadrupoleQuadrupole = 'CutoffQQ'
  character(*), parameter :: IdRFEpsilon                   = 'Epsilon'
  character(*), parameter :: IdFluctFreq                   = 'FluctFreq'
  character(*), parameter :: IdNFullFluct                  = 'NFullFluct'
  character(*), parameter :: IdMaxCounter                  = 'MaxCounter'

  character(*), parameter :: IdLambdaMin                   = 'LambdaMin'
  character(*), parameter :: IdLambdaMax                   = 'LambdaMax'
  character(*), parameter :: IdNBins                       = 'NBins'
  character(*), parameter :: IdLambdaStepMax               = 'LambdaStepMax'
  character(*), parameter :: IdLambdaExponent              = 'LambdaExponent'

  character(*), parameter :: IdBlockSizeCF                 = 'ResultFreqCF'
  character(*), parameter :: IdCorrFun                     = 'CorrfunMode'
  character(*), parameter :: IdCorrlength                  = 'Corrlength'
  character(*), parameter :: IdNStepcf                     = 'StepsCorrfun'
  character(*), parameter :: IdSpancf                      = 'SpanCorrfun'
  character(*), parameter :: IdNviewcf                     = 'ViewCorrfun'

  ! Define identifiers used in potential model file
  character(*), parameter :: IdSite_ntypes                 = 'NSiteTypes'
  character(*), parameter :: IdSite_stype                  = 'SiteType' !Mie-Potential or LJ126
  character(*), parameter :: IdSite_NMIEnm                 = 'NSites' !Mie-Potential or LJ126
  character(*), parameter :: IdSite_NCharge                = 'NSites'
  character(*), parameter :: IdSite_NDipole                = 'NSites'
  character(*), parameter :: IdSite_NQuadrupole            = 'NSites'
  character(*), parameter :: IdSite_NDFRot                 = 'NRotAxes'
  character(*), parameter :: IdSite_Mass                   = 'TotalMass'
  character(*), parameter :: IdSite_MOI1                   = 'InertMomX'
  character(*), parameter :: IdSite_MOI2                   = 'InertMomY'
  character(*), parameter :: IdSite_MOI3                   = 'InertMomZ'
  character(*), parameter :: IdMIE_n                       = 'MIE_n'
  character(*), parameter :: IdMIE_m                       = 'MIE_m'
  character(*), parameter :: IdMIEnm_r1                    = 'x'
  character(*), parameter :: IdMIEnm_r2                    = 'y'
  character(*), parameter :: IdMIEnm_r3                    = 'z'
  character(*), parameter :: IdMIEnm_sig                   = 'sigma'
  character(*), parameter :: IdMIEnm_eps                   = 'epsilon'
  character(*), parameter :: IdMIEnm_mass                  = 'mass'
  character(*), parameter :: IdCharge_r1                   = 'x'
  character(*), parameter :: IdCharge_r2                   = 'y'
  character(*), parameter :: IdCharge_r3                   = 'z'
  character(*), parameter :: IdCharge_e                    = 'charge'
  character(*), parameter :: IdCharge_mass                 = 'mass'
  character(*), parameter :: IdCharge_shield               = 'shielding'
  character(*), parameter :: IdDipole_r1                   = 'x'
  character(*), parameter :: IdDipole_r2                   = 'y'
  character(*), parameter :: IdDipole_r3                   = 'z'
  character(*), parameter :: IdDipole_theta                = 'theta'
  character(*), parameter :: IdDipole_phi                  = 'phi'
  character(*), parameter :: IdDipole_D                    = 'dipole'
  character(*), parameter :: IdDipole_mass                 = 'mass'
  character(*), parameter :: IdDipole_shield               = 'shielding'
  character(*), parameter :: IdQuadrupole_r1               = 'x'
  character(*), parameter :: IdQuadrupole_r2               = 'y'
  character(*), parameter :: IdQuadrupole_r3               = 'z'
  character(*), parameter :: IdQuadrupole_theta            = 'theta'
  character(*), parameter :: IdQuadrupole_phi              = 'phi'
  character(*), parameter :: IdQuadrupole_Q                = 'quadrupole'
  character(*), parameter :: IdQuadrupole_mass             = 'mass'
  character(*), parameter :: IdQuadrupole_shield           = 'shielding'
  character(*), parameter :: IdNFluct                      = 'NFluct'

#if CONSTR > 0
  character(*), parameter :: IdNCons                       = 'NConstr'
  character(*), parameter :: IdCons1Comp                   = 'Constr1Typ'
  character(*), parameter :: IdCons2Comp                   = 'Constr2Typ'
  character(*), parameter :: IdCons1                       = 'Constr1'
  character(*), parameter :: IdCons2                       = 'Constr2'
  character(*), parameter :: IdConsR                       = 'ConstrDist'
#endif
  character(*), parameter :: IdOptPressure                 = 'CalcPressure'
  character(*), parameter :: IdCommonEqui                  = 'CommonEqui'

! Calculation of residence times
  character(*), parameter :: IdResidTime                   = 'ResidTime'
  character(*), parameter :: IdResidComp1                  = 'ResidComp1'
  character(*), parameter :: IdResidSite1                  = 'ResidSite1'
  character(*), parameter :: IdResidComp2                  = 'ResidComp2'
  character(*), parameter :: IdResidSite2                  = 'ResidSite2'
  character(*), parameter :: IdResidPeriod                 = 'ResidPeriod'
  character(*), parameter :: IdResidLength                 = 'ResidLength'
  character(*), parameter :: IdResidBreak                  = 'ResidBreak'

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
  real(RK), parameter :: Ninth = 1._RK / 9._RK

  ! General physical constants
  real(RK), parameter :: NAvogadro = 6.022137E23_RK
  real(RK), parameter :: kBoltzmann = 1.380658E-23_RK
  real(RK), parameter :: ElementaryCharge = 1.602177E-19_RK
  real(RK), parameter :: VacuumPermittivity = 8.854188E-12_RK

  ! Some physical units
  real(RK), parameter :: Angstroem = 1E-10_RK
  real(RK), parameter :: DegreesInRadian = 180._RK / Pi
  real(8)             :: DebyesInSI
  real(8)             :: BuckinghamsInSI
  real(RK)            :: kForceOsmoticPressure


  ! Version of the parameter file
  real(RK) :: parVersionNr
  
  ! Walltime settings
  integer :: max_time
  integer :: time_limit
  
  ! Use reduced units for temperature, pressure, density
  logical :: UseReducedUnits
  
  ! LJ126 or Mie-Potential
  character(16) :: LJorMIE

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

  ! Type of LongRange Mode
  character(80)      :: LongRangeString
  integer, parameter :: Ewald         = 1
  integer, parameter :: RField        = 2
  integer, parameter :: PME           = 3
  integer, parameter :: extRField     = 4
  integer, parameter :: Rodgers       = 5
  integer            :: LongRange

  ! Type of method for chemical potential
  integer, parameter :: ChemPotMethodNone    = 0
  integer, parameter :: ChemPotMethodWidom   = 1
  integer, parameter :: ChemPotMethodGradIns = 2
  integer, parameter :: ChemPotMethodThermoInt = 3

  ! Type of method for weighting factors
  integer, parameter :: WFMethodNone   = 0
  integer, parameter :: WFMethodAuto   = 1
  integer, parameter :: WFMethodGuess  = 2
  integer, parameter :: WFMethodOptSet = 3
  
  integer, parameter :: CCritTypeVapor  = 0  
  integer, parameter :: CCritTypeGridvap = 2
  integer, parameter :: CCritTypeGridliq = 3

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
  integer,  parameter :: TransferRateLimit = 50

  ! Frequency of updating MC displacements
  integer, parameter :: DispUpdateFrequency = 100

  ! Number of simulation time steps
  integer :: NSteps

  ! Number of MC overlap reduction steps
  integer :: NStepsMC

  ! Number of NVT equilibration time steps
  integer :: NStepsV

  ! Number of NVE equilibration time steps
  integer :: NStepsE

  ! Number of NPT equilibration time steps
  integer :: NStepsP

  ! Number of NPH equilibration time steps
  integer :: NStepsH

  ! Number of gradual insertion initialization steps
  integer :: GradInsInit
  
  ! Number of orientations for second virial coefficient
  integer :: NOrient

  ! Radii for second virial coefficient
  real(RK) :: MinRadius, MaxRadius

  ! Number of current time step
  integer :: Step, StepTotal

  ! Equilibration flags
  logical :: Equilibration, NVTEquilibration, MCOverlapReduction, GradInsInitialization

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

#if TRANS == 1
  ! Maximum number of blocks CF
  integer :: NBlocksMaxCF

  ! Frequency of updating result file CF
  integer :: BlockSizeCF

  ! Maximum number of block sizes for error calculation CF
  integer :: NBlockSizesMaxCF

  ! Number of block sizes for error calculation CF
  integer :: NBlockSizesCF

  ! Current number of blocks CF
  integer :: NBlocksCF
#endif
  
  ! Kirkwood-Buff integration parameters
  ! Maximum number of blocks KBI
  integer :: NBlocksMaxKBI

  ! Frequency of updating result file KBI
  integer :: BlockSizeKBI

  ! Maximum number of block sizes for error calculation KBI
  integer :: NBlockSizesMaxKBI

  ! Number of block sizes for error calculation KBI
  integer :: NBlockSizesKBI

  ! Current number of blocks KBI
  integer :: NBlocksKBI

  ! Frequency of updating final result file
  integer :: ErrorsUpdateFrequency

  ! Frequency of updating visualisation file
  integer :: VisualUpdateFrequency

  !DC NOTE- Frequency of updating visualisation file
  integer :: VisualCCUpdateFrequency

  ! Frequency of updating RDF file
  integer :: RDFUpdateFrequency
  
  ! Number of RDF shells
  integer :: RDFNumberShells

  ! Number of ODF shells
  integer :: nR
  
  ! Discretisation of angle phi for ODF
  integer :: nPhi
  
  ! Discretisation of angle gamma12 for ODF
  integer :: nGamma
  
  ! Frequency of updating ODF file
  integer :: ODFUpdateFrequency
  
  ! Frequency of creating ODF output files
  integer :: ODFOutputFrequency  
  
  ! Frequency of updating KBI file
  integer :: KBIUpdateFrequency
  
  ! Frequency of updating Alpha2 displacement
  integer :: ALPHA2UpdateFrequency
  integer :: ALPHA2Length
  integer :: ALPHA2Shift

  !EinsteinCoef variables
  logical :: EinsteinCoefCalc

  ! Number of KBI shells
  integer :: KBINumberShells
  integer :: KBINumberShellsMax
  integer :: KBINShellsCubeEdge

  ! Number of density profile bins
  integer :: NBinsDen

  ! Common equilibration flag for MC. Determines whether one shared
  ! equilibration is performed
  logical :: CommonEqui

 ! Frequency of updating log file
  integer, parameter :: LogUpdateFrequency = 1000

  ! Internal variables of random number generator
  integer, parameter :: K4B = selected_int_kind(9)
  integer(K4B)       :: ix, iy, tpix, randk
  real(RK)           :: am

  ! Internal variable of FileReadParameter
  integer :: FileReadParameter_LineNumber = 0

  ! MPI variables
#if MPI_VER > 0
  integer :: ierror
  integer :: Communicator   ! actual MPI communicator
  !integer :: Communicator_W    ! =MPI_COMM_WORLD
  integer :: Communicator_R ! MPI communicator containing all roots
  integer :: NProcs ! number of PEs within actual MPI communicator
  integer :: NProc  ! MPI rank of actual MPI communicator
  integer :: NRootProc  ! MPI rank of root of actual MPI communicator
  logical :: RootProc   ! is PE root within actual MPI communicator
  integer :: NProcs_W   ! number of PEs within MPI_COMM_WORLD
  integer :: NProc_W    ! MPI rank within MPI_COMM_WORLD
  integer :: NRootProc_W    ! MPI rank of root PE within MPI_COMM_WORLD
  logical :: RootProc_W     ! is PE root of MPI_COMM_WORLD?
  integer :: NProcs_R   ! number of PEs within actual Communicator_R
  integer :: NProc_R    ! MPI rank within actual Communicator_R
  integer :: NRootProc_R    ! MPI rank of root PE within actual Communicator_R
  logical :: RootProc_R     ! is PE root of actual Communicator_R?
  integer :: NCommunicators ! number of Communicators (useful after MPI_Comm_Split)
  integer :: NCommunicator  ! ID of the Communicator
  !
  !integer, parameter :: mpimsgtag_log    = 0
  integer, parameter :: mpimsgtag_simTerm = 1
#else
  integer, parameter :: NProcs       = 1
  integer, parameter :: NProc        = 0
  integer, parameter :: NRootProc    = 0
  logical, parameter :: RootProc     = .true.
  integer, parameter :: NProcs_W     = 1
  integer, parameter :: NProc_W      = 0
  integer, parameter :: NRootProc_W  = 0
  logical, parameter :: RootProc_W   = .true.
  integer, parameter :: NCommunicators = 0
  integer, parameter :: NCommunicator = 0
#endif

#if ARCH == 1 || ARCH == 2 || ARCH == 3
  ! Flag for catched terminate signal
  logical :: TerminateProgram = .false.

! PGF compiler version < 6.0 seems to need this
! #if defined _PGF || defined __PGI
!   ! External funtion for signal handling
!   external SetTerminateProgram
! #endif

#else
  logical, parameter :: TerminateProgram = .false.
#endif
  integer :: TerminateStatus = 0

  integer, parameter :: IdErrorCodeBase = b'1000000000000000'   !=32768
  ! e.g. 10000 would be better to read for pure addition, but
  ! bits might code error type, origin (module&function),...


!==============================================================!
!  Global procedure interfaces                                 !
!==============================================================!

#if MPI_VER > 0
  interface SetCommunicator
    module procedure Global_SetCommunicator
  end interface

  interface SplitCommunicator
    module procedure Global_SplitCommunicator
  end interface
#endif

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
!#if MPI_VER > 0
!    module procedure Global_LogWrite_MPI
!#else
    module procedure Global_LogWrite
!#endif
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

# if MPI_VER > 0
  interface FileRewrite_parallel
    module procedure Global_FileRewrite_parallel
  end interface

  interface FileWriteNoAdvance_parallel
    module procedure Global_FileWriteNoAdvance_parallel
  end interface

  interface FileAppend_parallel
    module procedure Global_FileAppend_parallel
  end interface

  interface FileClose_parallel
    module procedure Global_FileClose_parallel
  end interface
#endif

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
&                    Global_FileReadParameter_RKdim1
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

  interface strlen_trim
    module procedure Global_String_Len_Trim
  end interface

  interface strtrim
    module procedure Global_String_TrimR
  end interface
  
  interface strtrimr
    module procedure Global_String_TrimR
  end interface
  
  interface strtriml
    module procedure Global_String_TrimL
  end interface
  
  interface strtrimlr
    module procedure Global_String_TrimLR
  end interface

  interface ProcRange
    module procedure Global_GetProcRange
  end interface

#ifdef USE_PRINTPROCSTATUS
  interface printProcStatus
    module procedure Global_printProcStatus
  end interface
#endif

!==============================================================!
!  External (intrinsic) procedures                             !
!==============================================================!

#if !( defined _WIN32 || defined __GNUC__ )

#if ARCH == 1 || ARCH == 2 || ARCH == 3
  ! Command line arguments
  integer, external :: iargc
  external getarg
#endif

#ifndef __INTEL_COMPILER

#if ARCH == 1 || ARCH == 2 || ARCH == 3
  ! Flush of I/O units
  external flush
  
  ! get/set file position
  integer, external :: ftell
#ifdef __GNUC__
  external fseek
#else
  integer, external :: fseek
#endif
  
  ! change current directory
#if defined _PGF || defined __PGI
  integer, external :: chdir
!#elif defined
  !external chdir
#endif

  ! Signal handler
  integer, external :: signal
#endif

  ! User name from console
#if ARCH == 1 || defined _PGF || defined __PGI
  character(256), external :: getlog
#elif ARCH == 2 || ARCH==3
  external getlog
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

  subroutine Global_printVersion()
    implicit none
    if (RootProc) then
      print *, trim( ProgramFileName ), ' Version: ', VersionString &
&            , ' (compiled at ', CompileTime, ')'
    end if
  end subroutine Global_printVersion

  subroutine Global_printUsage()
    implicit none
    if (RootProc) then
      print *, 'usage: ' &
&            , trim(ProgramFileName), '[-V|--version] [-h|--help] [-r|--restart] <par-file>[' &
&            , ParameterFileExtension, '] [<OutputPrefix>]', '}'
    end if
  end subroutine Global_printUsage
  

!==============================================================!

#if MPI_VER > 0

!==============================================================!
!  Subroutine Global_SetCommunicator                           !
!==============================================================!
! setting Communicator, NProc, NProcs, NRootProc, RootProc

  subroutine Global_SetCommunicator(comm)

    implicit none

    ! Include MPI header
    include 'mpif.h'

    ! Declare arguments
    integer, intent(in) :: comm

    Communicator = comm
    if( Communicator /= MPI_COMM_NULL ) then
      call MPI_Comm_size( Communicator, NProcs, ierror )
      call MPI_Comm_rank( Communicator, NProc, ierror )
    else
      NProcs = 0
      NProc = -1
    end if
    NRootProc = 0
    RootProc = NProc == NRootProc
    
  end subroutine Global_SetCommunicator

!==============================================================!
!  Subroutine Global_SplitCommunicator                         !
!==============================================================!
! setting NCommunicator, NCommunicators (& calling SetCommunicator)

  subroutine Global_SplitCommunicator(ngroups)

    implicit none

    ! Include MPI header
    include 'mpif.h'

    ! Declare arguments
    integer, intent(in)         :: ngroups
    
    integer :: groupId
    integer :: oldCommunicator,newCommunicator
    
    oldCommunicator=Communicator
    
    if( ngroups > NProcs ) then
      NCommunicators=NProcs
    else
      NCommunicators=ngroups
    endif
    
    write( IOBuffer, '("splitting communicator with",I4," PEs to ",I3," subcommunicators")') NProcs, NCommunicators
    call LogWrite
    write( IOBuffer, '("closing (and reopening) logfile - opening ",I3," additional new logfile(s) ",A,"_*",A," ...")') &
&          NCommunicators-1,trim(OutputNameTag), LogFileExtension
    call LogWrite
    call LogWriteBlank

    write( IOBuffer, '(72("#"))')
    call LogWrite
    call LogWriteBlank
    ! close log file to reopen/open new ones
    call LogClose
    
    !NCommunicator=mod(NProc,NCommunicators)
    NCommunicator=NProc*NCommunicators/NProcs
    ! NCommunicator -> color, NProc -> key (NProc_W also could be used)
    call MPI_Comm_Split(oldCommunicator,NCommunicator,NProc,newCommunicator,ierror)
    ! MPI_Comm_Group + MPI_Group_Range_incl + MPI_Comm_Create might be more efficient
    ! (avoiding some internal communication within the MPI library)    
    call SetCommunicator(newCommunicator)   !   RootProc is now true for the root of the new communicator(s)
    ! (re)open log files
    call LogOpen
    
    ! creating a communicator for all the RootProc (resp. non-RootProc) within the old communicator
    if (RootProc) then
      groupId=0
    else
      groupId=1
    endif
    call MPI_Comm_Split(oldCommunicator,groupId,NProc_W,Communicator_R,ierror)
    call MPI_Comm_size( Communicator_R, NProcs_R, ierror )
    call MPI_Comm_rank( Communicator_R, NProc_R, ierror )
    NRootProc_R = 0
    RootProc_R = NProc_R == NRootProc_R ! =RootProc_W

    !write(IOBuffer, '("after MPI_Comm_Split: NProc_W RootProc_W=",I6,L2," NProc RootProc=",I6,L2," NProc_R RootProc_R=",I6,L2)') NProc_W, RootProc_W, NProc,RootProc, NProc_R,RootProc_R
    !call LogWrite

  end subroutine Global_SplitCommunicator

#endif


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
    integer :: stat
#if ARCH == 1 || ARCH == 2 || ARCH == 3
    integer                   :: narg, dot, i
    integer                   :: argpos
    character(IOBufferLength) :: buffer
#endif
#if MPI_VER > 0
    integer                                    :: mpiversion, mpisubversion
    character*(MPI_MAX_PROCESSOR_NAME)         :: procname
    integer                                    :: procnamelen
    character*(MPI_MAX_PROCESSOR_NAME),pointer, contiguous :: procnames(:)
    integer                                    :: hostrank = MPI_PROC_NULL
    integer                                    :: iorank = MPI_PROC_NULL
    integer,pointer, contiguous                            :: ioranks(:)
    logical                                    :: flag
#endif
#ifdef ENABLE_OMP
    integer                                    :: OMP_GET_MAX_THREADS, OMP_GET_NUM_PROCS
    integer                                    :: ompmaxnumthreads, ompnumprocs
#endif
#if ARCH == 1 || ARCH == 2 || ARCH == 3
    character(IOBufferLength) :: hostname = 'unknown host'
    character(IOBufferLength) :: username = 'unknown user'
#else
    character(*), parameter   :: hostname = 'unknown host'
    character(*), parameter   :: username = 'unknown user'
#endif
    integer :: length
    character :: Version*6
    integer :: color

    ! Initialize MPI
#if MPI_VER > 0
    call MPI_Init( ierror )
    call SetCommunicator( MPI_COMM_WORLD )
    NProcs_W = NProcs
    NProc_W = NProc
    NRootProc_W = NRootProc
    RootProc_W = RootProc
    NCommunicators=1
    NCommunicator=0
    
    ! just do something to initialize the communicator between the (root of the) communicators
    !Communicator_R=MPI_COMM_NULL
    if (RootProc) then
      color=0
    else
      color=1
    endif
    call MPI_Comm_Split(Communicator,color,NProc,Communicator_R,ierror)
    
    ! better define and initialize as parameter...
    if ( RK == 8 ) then
      !MPI_RK = MPI_DOUBLE_PRECISION
      MPI_RK = MPI_REAL8
    else if ( RK == 4 ) then
      !MPI_RK = MPI_REAL
      MPI_RK = MPI_REAL4
    else
      if( RootProc ) then
        print *,"ERROR: RK==",RK," not supported for MPI version"
      end if
      call MPI_Abort( MPI_COMM_WORLD, 1, ierror )
    end if
#endif

    ! Initialize flags
    Restart = .false.
    OutputNameTagfromCommandline = .false.

#if ARCH == 1 || ARCH == 2 || ARCH == 3
    ! Read command line parameter
    if( RootProc ) then
      ! first argument is binary file path, which was executed
      call getarg( 0, buffer )
      i = scan( buffer, FileSep, BACK=.true. )
      if( i > 0 ) then
        ProgramFileName = trim( buffer( i+1:len( buffer ) ) )
      else
        ProgramFileName = trim( buffer )         ! possible truncation?
      end if
#if defined __PATHSCALE__ || defined _CRAYFTN
! this should work for all Fortran2003 compilers
      narg = command_argument_count()
#else
      narg = iargc()
#endif
      if( narg .lt. 1 ) then
    call Global_printVersion()
    call Global_printUsage()
        ! Abort program
#if MPI_VER > 0
        call MPI_Abort( MPI_COMM_WORLD, 2, ierror )
#endif
        stop
      end if
      ! processing command line arguments
      do i = 1,narg
        argpos=i
        call getarg( argpos, buffer )
        !print *,"processing command line argument ",trim(buffer)
        if (trim(buffer).eq."-V" .or. trim(buffer).eq."--version") then
          call Global_printVersion()
#if MPI_VER > 0
          call MPI_Finalize( ierror )
#endif
          stop
        else if (trim(buffer).eq."-h" .or. trim(buffer).eq."--help") then
          call Global_printUsage()
#if MPI_VER > 0
          call MPI_Finalize( ierror )
#endif
          stop
        else if (trim(buffer).eq."-r" .or. trim(buffer).eq."--restart") then
          Restart = .true.
    else
    !  print *,"WARNING: command line argument not known and disregarded: ",trim(buffer)
      exit
        end if
      end do
      if (argpos>narg) then
#if MPI_VER > 0
          call MPI_Finalize( ierror )
#endif
          stop
      end if
      ! next argument should be the input file name
      call getarg( argpos, buffer )
      argpos=argpos+1
      ! 
      buffer = trim( buffer )
      ParameterFileName =  trim(buffer)

      ! separate directory and filename
      i = scan(buffer, FileSep, .true.)
      if( i>0 ) then
        ! path includes directory
#if defined __INTEL_COMPILER || defined _PGF || defined __PGI || defined __PATHSCALE__ 
        stat = chdir( buffer(:max(i-1,1)) )
#elif defined _CRAYFTN
        call PXFCHDIR( buffer(:max(i-1,1)), 0, stat)
#elif ARCH==3 || defined __GNUC__ 
        call chdir( buffer(:max(i-1,1)), stat )
#else
        print *, 'chdir not supported!'
        stat=-1
        i=0
#endif
        if( stat==0 ) then
          print *, 'chdir to ', trim(buffer(:max(i-1,1)))
        else
          print *, 'cannot change to ', trim(buffer(:max(i-1,1))), ' stat=', stat
        end if
        buffer=trim(buffer(i+1:))
      end if

      dot = index( buffer, '.', BACK=.true. )
      if( dot > 0 ) then
        if( buffer( dot:len( buffer ) ) .eq. ParameterFileExtension ) then
!           buffer = buffer( 1:dot - 1 )
          ParameterFileName =  trim( buffer )    ! possible truncation
    !else
        !  ParameterFileName =  trim(buffer)//ParameterFileExtension
        end if
        !RestartFileName=trim(buffer(1:dot-1))//RestartFileExtension
      end if

      if( narg .ge. argpos ) then
        ! if present, the third argument should be the output file name
        call getarg( argpos, buffer )
        OutputNameTagfromCommandline = .true.
      else
        ! otherwise use the input file name without extension
        buffer = buffer( 1:dot - 1 )
      end if
      OutputNameTag = trim( buffer )             ! possible truncation
    end if
#else
    ! Plattform does not support command line parameters
    ProgramFileName='ms2'
    ParameterFileName = 'ms2.par'
    OutputNameTag='ms2out'
#endif
    RestartFileName=trim(OutputNameTag)//RestartFileExtension

#if MPI_VER > 0
    call MPI_Bcast( Restart, 1, MPI_LOGICAL, NRootProc, Communicator, ierror )
    call MPI_Bcast( OutputNameTag, len(OutputNameTag), MPI_CHARACTER, NRootProc, Communicator, ierror )
#endif

    ! Open log file
    call LogOpen

    ! Open log file
    call LogWriteBlank
    write( IOBuffer, '(74("*"))')
    call LogWrite
    write( IOBuffer, '("*                         Molecular Simulation 2                         *")')
    call LogWrite
    write( IOBuffer, '(74("*"))')
    call LogWrite
    call LogWriteBlank
    write( IOBuffer, '(74("*"))')
    call LogWrite
    write( IOBuffer, '("*                         Publishing with ms2                            *")')
    call LogWrite
    write( IOBuffer, '("* Every user agrees to cite ms2 upon usage as follows                    *")')
    call LogWrite
    write( IOBuffer, '("* ---------------------------------------------------------------------- *")')
    call LogWrite
    write( IOBuffer, '("* G. Rutkai, A. Koester, G. Guevara-Carrion, T. Janzen, M. Schappal,     *")')
    call LogWrite
    write( IOBuffer, '("* C.W. Glass, M. Bernreuther, A. Wafai, S. Stephan, M. Kohns, S. Reiser, *")')
    call LogWrite
    write( IOBuffer, '("* S. Deublein, M. Horsch, H. Hasse, J. Vrabec                            *")')
    call LogWrite
    write( IOBuffer, '("* Computer Physics Communications (2017)                                 *")')
    call LogWrite
    write( IOBuffer, '(74("*"))')
    call LogWrite
    call LogWriteBlank
    write( IOBuffer, '(74("*"))')
    call LogWrite
    write( IOBuffer, '("* (c) by TU Kaiserslautern / U Paderborn                                 *")')
    call LogWrite
    write( IOBuffer, '("*     P.O. Box 67653                                                     *")')
    call LogWrite
    write( IOBuffer, '("*     67653 Kaiserslautern                                               *")')
    call LogWrite
    write( IOBuffer, '(74("*"))')
    call LogWrite
    call LogWriteBlank
    write( IOBuffer, '(74("*"))')
    call LogWrite
    write( IOBuffer, '("* Updates and auxiliary routines are available from                      *")')
    call LogWrite
    write( IOBuffer, '("* http://www.ms-2.de                                                     *")')
    call LogWrite
    write( IOBuffer, '(74("*"))')
    call LogWrite
    call LogWriteBlank
    write( IOBuffer, '("Program ", A, " version ", A)' ) trim( ProgramFileName ), trim( VersionString )
    call LogWrite
    write( IOBuffer, '("Hardware architecture: ", A)' ) Hardware
    call LogWrite

! cmp. http://predef.sourceforge.net/precomp.html
!           __GFORTRAN__
#if defined _CRAYFTN
    call GET_ENVIRONMENT_VARIABLE('CRAY_CC_VERSION',Version,length,stat,.FALSE.)
    write( IOBuffer, '("Compiler version     : CRAYFTN ftn ", A6)' )   Version 
#elif defined __GNUC__
    write( IOBuffer, '("Compiler version     : GNU gfortran", I6)' ) __GNUC_VERSION__
#elif defined __INTEL_COMPILER
    write( IOBuffer, '("Compiler version     : INTEL ", I4, ", build ", I8)' ) __INTEL_COMPILER, __INTEL_COMPILER_BUILD_DATE
#elif defined __PGI || defined _PGF
    write( IOBuffer, '("Compiler version     : PGI pgf")' )
#elif defined __SUNPRO_F95
    write( IOBuffer, '("Compiler version     : SUN studio sunf95 ", A)' ) MACRODEF_TO_STRING(__SUNPRO_F95)
#elif defined __SUNPRO_F90
    write( IOBuffer, '("Compiler version     : SUN studio sunf90 ", A)' ) MACRODEF_TO_STRING(__SUNPRO_F90)
#else
!                                                         __VERSION__
    write( IOBuffer, '("Compiler version     : unknown")' )
#endif
    call LogWrite

    write( IOBuffer, '("Compiler flags       :")' )
    call LogWriteNoAdvance
#if MPI_VER > 0
    write( IOBuffer, '(" MPI=1")' )
    call LogWriteNoAdvance
#endif
#if TRANS == 1
    write( IOBuffer, '(" TRANS=1")' )
    call LogWriteNoAdvance
#endif
#if HBOND == 1
    write( IOBuffer, '(" HBOND=1")' )
    call LogWriteNoAdvance
#endif
#if OSMOP == 1
    write( IOBuffer, '(" OSMOP=1")' )
    call LogWriteNoAdvance
#endif
#if OSMOP == 2
    write( IOBuffer, '(" OSMOP=2")' )
    call LogWriteNoAdvance
#endif
    ! new compiler flags should be added
    ! include target, omp and precision???
    write( IOBuffer, '(" ")' )
    call LogWrite

    write( IOBuffer, '("Compile time         : ", A)' ) CompileTime
    call LogWrite
    write( IOBuffer, '("Real Kind            :", I2)' ) RK
    call LogWrite
    call LogWriteBlank

    ! Get name of host and user
#if ARCH == 1  || defined _CRAYFTN
    call getenv( 'HOSTNAME', hostname )
#elif ARCH == 2 || ARCH == 3
#if defined _PGF || defined __PGI || defined __GNUC__ || defined __PATHSCALE__ || defined __SUNPRO_F90 || ARCH == 3
    i = hostnm( hostname )
#else
    i = hostnam( hostname )
#endif
    if( i .ne. 0 ) hostname = 'unknown host'
#endif
#ifdef _CRAYFTN
   username = 'Getlog is not supported'
#elif ARCH == 1 || defined _PGF || defined __PGI
    username = getlog()
#elif ARCH == 2 || ARCH == 3
    call getlog( username )
#endif
    write( IOBuffer, '("Hostname             : ", A)' ) trim( hostname )
    call LogWrite
    write( IOBuffer, '("started by user ", A)' ) trim( username )
    call LogWriteTime
    call LogWriteBlank
    write( IOBuffer, '(72("-"))')
    call LogWrite

    write( IOBuffer, '("Parallelization:")' )
    call LogWrite

    ! Update log file
#if MPI_VER > 0
    nullify( procnames )
    nullify( ioranks )
    if( RootProc ) then
      call MPI_Get_version(mpiversion, mpisubversion, ierror)
      write( IOBuffer, '("MPI version        :",I2,".",I1)' ) mpiversion, mpisubversion
      call LogWrite
      write( IOBuffer, '("Number of processes:",I4)' ) NProcs
      call LogWrite
      write( IOBuffer, '("Root process rank  :",I4)' ) NRootProc
      call LogWrite
      call MPI_Attr_get(Communicator, MPI_HOST, hostrank, flag, ierror)
      if(ierror==0 .and. flag .and. hostrank/=MPI_PROC_NULL ) then
        write( IOBuffer, '("MPI Host rank      :",I4)' ) hostrank
        call LogWrite
      end if
    end if
    allocate( procnames(NProcs), STAT = stat )
    allocate( ioranks(NProcs), STAT = stat )
    call MPI_Get_processor_name(procname, procnamelen, ierror)
    !                         procnamelen might be variable
    call MPI_Gather(procname, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER &
&                  ,procnames, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER &
&                  ,NRootProc, Communicator, ierror)
    call MPI_Attr_get(Communicator, MPI_IO, iorank, flag, ierror)
    call MPI_Gather(iorank, 1, MPI_INTEGER, ioranks, 1, MPI_INTEGER &
&                  ,NRootProc, Communicator, ierror)

    if( RootProc ) then
      write( IOBuffer, '("rank:  I/O: processor name:")' )
      call LogWrite
      do i = 1,NProcs
        if( ioranks(i) == MPI_ANY_SOURCE )  then
          write( IOBuffer, '(I5,"   +   ", A)' ) i-1, procnames(i)
        else
          write( IOBuffer, '(I5, 1X, I5, 1X, A)' ) i-1, ioranks(i), procnames(i)
        end if
        call LogWrite
      end do
    end if

    call LogWriteBlank
    if( associated( ioranks ) ) deallocate( ioranks )
    if( associated( procnames ) ) deallocate( procnames )
#else

#ifndef ENABLE_OMP
    write( IOBuffer, '("sequential Version")' )
    call LogWrite
#endif
#endif
#ifdef ENABLE_OMP
    write( IOBuffer, '("OpenMP enabled")' )
    call LogWrite
    ompmaxnumthreads = OMP_GET_MAX_THREADS()
    write( IOBuffer, '("Number of max. threads:",I4)' ) ompmaxnumthreads
    call LogWrite
    ompnumprocs = OMP_GET_NUM_PROCS()
    write( IOBuffer, '("Number of processors  :",I4)' ) ompnumprocs
    call LogWrite
#endif
    call LogWriteBlank

    ! Set signal handler
#if ARCH == 1 || ARCH == 2
#ifdef _CRAYFTN
#elif defined  __GNUC__
    call signal( 1, IgnoreSignal )  ! Ignore SIGHUP
    call signal( 2, SetTerminateProgram )   ! Catch SIGINT
    call signal( 3, SetTerminateProgram )   ! Catch SIGQUIT
    call signal( 15, SetTerminateProgram )  ! Catch SIGTERM
#else
    i = signal( 1, SetTerminateProgram, 1 ) ! Ignore SIGHUP (HangUP)
    i = signal( 2, SetTerminateProgram, -1 )    ! Catch SIGINT (INTerrupt)
    i = signal( 3, SetTerminateProgram, -1 )    ! Catch SIGQUIT (QUIT)
    i = signal( 15, SetTerminateProgram, -1 )   ! Catch SIGTERM (TERMinate)
#endif
#elif ARCH == 3
    i = signal( 15, SetTerminateProgram )   ! Catch SIGTERM
#endif
    write( IOBuffer, '(72("-"))')
    call LogWrite
#if ARCH == 1 || ARCH == 2 || ARCH == 3
    if( i < 0 ) then
      call Warning('Cannot set signal handler')
    else
      write( IOBuffer, '("Signal handler set successfully")' )
      call LogWrite
       !call LogWriteBlank
    end if
#endif

    ! Initialize random number generator
    call Randomize( seed = 5333 )

    ! Define some constants
    limits_RK_MAX = huge(limits_RK_MAX)
    exp_arg_max = log(limits_RK_MAX)

#ifdef SINGLEPRECISION
    DebyesInSI = real( sqrt( 1E49_8 / (4._RK * real(Pi, RK) * real(VacuumPermittivity, RK) ) ), 8 )
    BuckinghamsInSI = real( sqrt( 1E69_8 / (4._RK * real(Pi, RK) * real(VacuumPermittivity, RK) ) ), 8 )
#else
    DebyesInSI = sqrt( 1E49_RK / (4._RK * Pi * VacuumPermittivity) )
    BuckinghamsInSI = sqrt( 1E69_RK / (4._RK * Pi * VacuumPermittivity) )
#endif

#ifdef USE_PRINTPROCSTATUS
    call printProcStatus("end of InitializeProgram")
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

#ifdef USE_PRINTPROCSTATUS
    call printProcStatus("beginning of FinalizeProgram")
#endif

    call LogWriteBlank
    write( IOBuffer, '(72("*"))')
    call LogWrite
    write( IOBuffer, '("Program terminated")' )
    call LogWriteTime
    write( IOBuffer, '(72("*"))')
    call LogWrite
    
    ! Close log file
    call LogClose


    ! Finalize MPI
#if MPI_VER > 0
    call MPI_Finalize( ierror )
#endif

  end subroutine Global_FinalizeProgram



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

  subroutine Global_Error( ErrorString, ErrorCode )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    character(*), intent(in), optional :: ErrorString
    integer, intent(in), optional :: ErrorCode
    integer :: GlobalErrorCode = IdErrorCodeBase
    
    ! Output error message (might not show up in the MPI version if not initiated by NRootProc!)
    call LogWriteBlank
    if( present( ErrorString ) ) then
      IOBuffer = 'ERROR: '// trim( ErrorString )
    else
      IOBuffer = 'ERROR: '// trim( ErrorBuffer ) ! possible truncation
    end if
    if( RootProc ) print *, trim( IOBuffer )
    call LogWrite

    if( present( ErrorCode ) ) then
      GlobalErrorCode=IdErrorCodeBase+ErrorCode
      !GlobalErrorCode=ior(IdErrorCodeBase,ErrorCode)
    end if

    call LogWriteBlank
    write( IOBuffer, '(72("*"))')
    call LogWrite
    write( IOBuffer, '("Program terminated with Error (",I5,")")' ) GlobalErrorCode
    call LogWriteTime
    write( IOBuffer, '(72("*"))')
    call LogWrite
    
    ! Close log file
    call LogClose

    ! Abort program
#if MPI_VER > 0
    ! ErrorCode will be used (at least) by MPI...
    call MPI_Abort( MPI_COMM_WORLD, GlobalErrorCode, ierror )
#endif
    !    GlobalErrorCode is not a constant and therefore not accepted by older Fortran versions :-( ...
    stop IdErrorCodeBase
    !error stop IdErrorCodeBase ! this is an error, so error stop might be favorable
    !stop 4 ! very old Fortran versions only support char (0-255)
    ! should check for Fortran2008+ solution...

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
    call MPI_Allreduce( ok, okAll, 1, MPI_LOGICAL, MPI_LAND, Communicator, ierror )
    if( okAll ) return
#else
    if( stat == 0 ) return
#endif

    ! Terminate program
    if( present( NPart ) ) then
      write( ErrorBuffer, '("Cannot allocate memory for", I11, " ", A)' ) NPart, str
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
    character(FileNameLength) :: filename
    
    ! Check for root process
    if( .not. RootProc ) return

  
    ! using <OutputNameTag>.log, if only one communicator exists date_and_time
    ! and   <OutputNameTag>_<CommId>.log for several
    ! could be extended to <OutputNameTag>_<Phase>.<CommId>.log, for multiple communicator splits/phases
    
    ! generate filename
    if ( NCommunicators .gt. 1 .and. NCommunicator .gt. 0 ) then
      write( filename, '(A,"_",I0,A)' ) trim( OutputNameTag ),NCommunicator+1,LogFileExtension
    else
      write( filename, '(A,A)' ) trim( OutputNameTag ),LogFileExtension
    endif

    if ( NCommunicators .gt. 1 .and. NCommunicator .eq. 0 ) then
      call FileAppend( iounit_log, trim(filename) )
      write( IOBuffer, '("ms2 logfile ",A," reopened")' ) trim(filename)
    else
      call FileRewrite( iounit_log, trim(filename) )
      write( IOBuffer, '("ms2 logfile ",A," created")' ) trim(filename)
    endif
#if MPI_VER > 0
    write( IOBuffer(len_trim(IOBuffer)+1:), '(" by PE",I0)' ) NProc_W
#endif

    call LogWriteTime
    !call LogWriteBlank
      
  end subroutine Global_LogOpen



!==============================================================!
!  Subroutine Global_LogClose                                  !
!==============================================================!

  subroutine Global_LogClose()

    implicit none

    ! Check for root process
    if( .not. RootProc ) return

    ! Close log file
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

! #if MPI_VER > 0
! !==============================================================!
! !  Subroutine Global_LogWrite_MPI                              !
! !==============================================================!
! 
! subroutine Global_LogWrite_MPI(rank)
! 
!     implicit none
!     include 'mpif.h'
!     
!     ! Declare local variables
!     integer, intent(in), optional      :: rank
!
!     integer             :: mpistatus(MPI_STATUS_SIZE)
!
!
!     if( present( rank ) .and. (rank .ne. NRootProc) ) then
!       ! transfer IOBuffer to NRootProc
!       call MPI_Sendrecv( IOBuffer, IOBufferLength, MPI_CHARACTER, NRootProc, mpimsgtag_log, &
! &                        IOBuffer, IOBufferLength, MPI_CHARACTER, rank,      mpimsgtag_log, &
! &                        Communicator, mpistatus, ierror)
!       !call MPI_Barrier( Communicator, ierror )
!     endif
!     ! execute LogWrite on NRootProc
!     if( RootProc ) call Global_LogWrite()
! 
!   end subroutine Global_LogWrite_MPI
! #endif


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
    write( IOBuffer, '(I9, " steps completed")' ) Step
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
    write( IOBuffer, '("Opening file <", A, "> for reading (unit",I5,")")' ) trim( filename ), iounit
    call LogWrite
    open( iounit, file = filename, action = 'READ', status = 'OLD', iostat = stat )
    if( stat /= 0 ) call Error( 'Cannot open file '//trim( filename )//' for reading' )

  end subroutine Global_FileReset

#if MPI_VER > 0
!==============================================================!
!  Subroutine Global_FileClose_parallel                        !
!==============================================================!

  subroutine Global_FileClose_parallel( iounit )

    implicit none

    ! Declare arguments
    integer, intent(in) :: iounit

    call MPI_File_Close(iounit, ierror)

    if( RootProc )then
        write( IOBuffer, '("File <", A, "> closed")' )"*.run or *.rav"
        call LogWrite
    endif

  end subroutine Global_FileClose_parallel


!==============================================================!
!  Subroutine Global_FileRewrite_parallel                      !
!==============================================================!

  subroutine Global_FileRewrite_parallel( iounit, filename )

    implicit none
    include 'mpif.h'
    ! Declare arguments
    integer                       :: iounit 
    character(*), intent(in)      :: filename

    if(RootProc) then
      ! open file for writing
      if( iounit /= iounit_log ) then
        write( iobuffer, '("opening file <", a, "> for writing")' ) trim( filename )
        call logwrite
        open( iounit, file = filename, action = 'WRITE', status = 'REPLACE' )
        close(iounit)
      end if
    end if
    call MPI_File_Open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, iounit, ierror)
    if(RootProc) then
      if( ierror .ne. 0 ) then
        write( IOBuffer,'(a,a)') 'Can not create ',trim( filename )
        call logwrite
      end if
    end if

  end subroutine Global_FileRewrite_parallel

!==============================================================!
!  Subroutine Global_FileAppend_parallel                       !
!==============================================================!

  subroutine Global_FileAppend_parallel( iounit, filename )

    implicit none
    include 'mpif.h'
    ! Declare arguments
    integer, intent(in)           :: iounit
    character(*), intent(in)      :: filename

    ! Declare local variables

    logical :: ex

    ! Check for root process
    if( RootProc ) then

      ! Open file for writing
      if( iounit /= iounit_log ) then
        write( IOBuffer, '("Opening file <", A, "> for appending")' ) trim( filename )
        call LogWrite
      end if

      inquire( file = filename, exist = ex )
      if( ex ) then
        open( iounit, file = filename, action = 'WRITE', status = 'OLD', position = 'APPEND' )
      else
        write( IOBuffer, '("File does not exist. Creating new")' )
        call LogWrite
        open( iounit, file = filename, action = 'WRITE', status = 'REPLACE' )
      end if
    endif
    !!! ERRONEOUS
    ! don't mix Fortran POSIX IO with mpi IO; Fortran units != MPI units; mpi iounit is intend(out) here...
    call MPI_File_Open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, iounit, ierror)
    if(RootProc) then
      if( ierror /= 0 ) then
        write( IOBuffer,'(a,a)') 'Can not create ',trim( filename )
        call logwrite
      end if
    end if

  end subroutine Global_FileAppend_parallel

!==============================================================!
!  Subroutine Global_FileWriteNoAdvance_parallel               !
!==============================================================!

  subroutine Global_FileWriteNoAdvance_parallel( iounit )

    implicit none
    include 'mpif.h'
    ! Declare arguments
    integer             :: mpistatus(MPI_STATUS_SIZE)
    integer, intent(in) :: iounit

    ! Write contents of buffer to file
    call MPI_File_write(iounit,IOBuffer, len(trim(IOBuffer)), MPI_CHARACTER, mpistatus, ierror)


  end subroutine Global_FileWriteNoAdvance_parallel

#endif

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
      write( IOBuffer, '("Opening file <", A, "> for writing (unit",I5,")")' ) trim( filename ), iounit
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
      write( IOBuffer, '("Opening file <", A, "> for appending (unit",I5,")")' ) trim( filename ), iounit
      call LogWrite
    end if
    inquire( file = filename, exist = ex )
    if( ex ) then
      open( iounit, file = filename, action = 'WRITE', status = 'OLD', position = 'APPEND' )
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
      write( IOBuffer, '("File <", A, "> closed (unit",I5,")")' ) trim( fn ), iounit
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

  function Global_FileReadParameter( iounit, parameterqualifiers, &
&                                    rewind_before, status ) &
&          result (parametervalue)


    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    integer, intent(in)                :: iounit
    character(*), intent(in)           :: parameterqualifiers
    logical, intent(in), optional      :: rewind_before
    integer, intent(out), optional     :: status

    character(IOBufferLength) :: parametervalue


    ! Declare local variables
    integer                   :: stat, comment_pos, linesread, i
    character(FileNameLength) :: fn
    logical                   :: foundqualifier = .false.
    character(IOBufferLength) :: parameterqualifier
    integer                   :: delimiterpos1, delimiterpos2

    ! determine filename
    inquire( iounit, NAME = fn )
    ! Only RootProc reads parameter from file
    if( RootProc ) then
      ! rewind file, if requested
      if( present(rewind_before) ) then
        if( rewind_before ) then
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
          call Error( "ERROR reading file "//trim(fn)// " while searching for parameter <"//parameterqualifiers//">" )
        ! end of file reached?
        elseif( stat < 0 ) then
          !call Warning( trim(fn)//": Could not find parameter <"//parameterqualifiers//">" )
          parametervalue = ''
          if( present(status) ) status = stat
          ! (try to) restore position
          if( present(rewind_before) ) then
            if( rewind_before ) then
              rewind( iounit )
              FileReadParameter_LineNumber = 0
              linesread = 0
              exit    !not nice!
            end if
          end if
          ! rewind to the position, where the reading process was started
          backspace( iounit )   ! "undo" last read, where eof was encountered
          do i = 1,linesread
            backspace( iounit )
            FileReadParameter_LineNumber = FileReadParameter_LineNumber - 1
          end do
          exit
        end if
        FileReadParameter_LineNumber = FileReadParameter_LineNumber + 1
        linesread = linesread + 1
        comment_pos = index( parametervalue, CommentSign )
        if( comment_pos > 0 ) then
          parametervalue = parametervalue(1:comment_pos - 1)
        end if
        delimiterpos2 = 0
        do ! test all qualifier alternatives (if parameterqualifier is a list delimited with :)
          delimiterpos1 = delimiterpos2
          if ( delimiterpos1>=len(trim(parameterqualifiers)) ) exit
          delimiterpos2 = delimiterpos1 + scan(trim(parameterqualifiers(delimiterpos1+1:)),":")
          if( delimiterpos2>delimiterpos1 ) then
            parameterqualifier = parameterqualifiers(delimiterpos1+1:delimiterpos2-1)
          else
            parameterqualifier = parameterqualifiers(delimiterpos1+1:)
          end if
          foundqualifier = index( strtriml( parametervalue ), trim( parameterqualifier ) ) == 1
          if( foundqualifier ) then
            ! extract value part (after =)
            parametervalue = parametervalue( index( parametervalue, '=' )+1:len( parametervalue ) )
            parametervalue = strtrimlr( parametervalue )
            if( present(status) ) status = 0
            exit
          end if
          if ( delimiterpos2<=delimiterpos1 ) exit
        end do
        if ( foundqualifier ) exit
      end do
    end if

    ! Broadcast parameter to other processes
    ! (2 Broadcast are not very efficient, but it doesn't need to be efficient here.
    !  Better broadcast the integer, float parametervalues, instead of the string?)
#if MPI_VER > 0
    ! RootProc knows length (len_trim) of parametervalue, but it's easier to bcast the whole buffer
    call MPI_Bcast( parametervalue, len(parametervalue), MPI_CHARACTER, NRootProc, Communicator, ierror )
    if( present(status) ) then
      call MPI_Bcast( status, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
    end if
#endif

  end function Global_FileReadParameter


!==============================================================!
!  Subroutine Global_FileReadParameter_buffer                  !
!==============================================================!

  subroutine Global_FileReadParameter_buffer( iounit, parameterqualifier, rewind_before, defaultvalue, status )
  ! this subroutine is for backward compatibily purposes

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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
        write( IOBuffer, '("setting ",A," (IOBuffer) to default value ",A)' ) trim(parameterqualifier), trim(defaultvalue)
        call LogWrite
        IOBuffer = defaultvalue
      else
        call Error( "Could not find parameter <"//parameterqualifier//">" )
        !return
      end if
    end if
    if ( present(status) ) status=stat

  end subroutine Global_FileReadParameter_buffer


!==============================================================!
!  Subroutine Global_FileReadParameter_String                  !
!==============================================================!

  subroutine Global_FileReadParameter_String( parametervariable, iounit, parameterqualifiers, &
&                                            rewind_before, defaultvalue, status )
  ! setting up functions with result (parametervalue) for different data types is ambigious
  ! for a FileReadParameter polymorphism

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    character(*), intent(out)          :: parametervariable
    integer, intent(in)                :: iounit
    character(*), intent(in)           :: parameterqualifiers
    logical, intent(in), optional      :: rewind_before
    character(*), intent(in), optional :: defaultvalue
    integer, intent(out), optional     :: status


    ! Declare local variables
    character(IOBufferLength) :: buffer
    integer                   :: stat
    !character(FileNameLength) :: fn

    buffer = Global_FileReadParameter(iounit, parameterqualifiers, rewind_before, stat)
    if ( stat == 0 ) then
      parametervariable = buffer
    else
      ! parameter could not be read
      if ( present(defaultvalue) ) then
        ! set default value
        write( IOBuffer, '("setting ",A," to default value ",A)' ) trim(parameterqualifiers), trim(defaultvalue)
        call LogWrite
        parametervariable = defaultvalue
      else if ( .not. present(status) ) then
        ! Terminate with error, if error can not be returned through status
        call Error( "Could not find parameter <"//parameterqualifiers//">" )
        !return
      end if
    end if
    if ( present(status) ) status=stat

  end subroutine Global_FileReadParameter_String


!==============================================================!
!  Subroutine Global_FileReadParameter_Int                     !
!==============================================================!

  subroutine Global_FileReadParameter_Int( parametervariable, iounit, parameterqualifiers, &
&                                         rewind_before, defaultvalue, status )
  ! Global_FileReadParameter_Integer has 32>31 characters!

    implicit none

    ! Declare arguments
    integer, intent(out)           :: parametervariable
    integer, intent(in)            :: iounit
    character(*), intent(in)       :: parameterqualifiers
    logical, intent(in), optional  :: rewind_before
    integer, intent(in), optional  :: defaultvalue
    integer, intent(out), optional :: status

    ! Declare local variables
    character(IOBufferLength) :: buffer
    integer                   :: stat
    !character(FileNameLength) :: fn

    buffer = Global_FileReadParameter(iounit, parameterqualifiers, rewind_before, stat)
    if ( stat == 0 ) then
      read( buffer, * ) parametervariable
    else !if ( stat < 0 ) then

      ! parameter could not be read
      if ( present(defaultvalue) ) then
        ! set default value
        write( IOBuffer, '("setting ",A," to default value ",I7)' ) trim(parameterqualifiers), defaultvalue
        call LogWrite
        parametervariable = defaultvalue
      else if ( .not. present(status) ) then
        ! Terminate with error, if error can not be returned through status
        call Error( "Could not find parameter <"//parameterqualifiers//">" )
        !return
      end if
    end if
    if ( present(status) ) status=stat

  end subroutine Global_FileReadParameter_Int


!==============================================================!
!  Subroutine Global_FileReadParameter_RK                      !
!==============================================================!

  subroutine Global_FileReadParameter_RK( parametervariable, iounit, parameterqualifiers, &
&                                        rewind_before, defaultvalue, status )

    implicit none

    ! Declare arguments
    real(RK), intent(out)          :: parametervariable
    integer, intent(in)            :: iounit
    character(*), intent(in)       :: parameterqualifiers
    logical, intent(in), optional  :: rewind_before
    real(RK), intent(in), optional :: defaultvalue
    integer, intent(out), optional :: status


    ! Declare local variables
    character(IOBufferLength) :: buffer
    integer                   :: stat
    !character(FileNameLength) :: fn

    buffer = Global_FileReadParameter(iounit, parameterqualifiers, rewind_before, stat)
    if ( stat == 0 ) then
      read( buffer, * ) parametervariable

    else !if ( stat < 0 ) then
      ! parameter could not be read
      if ( present(defaultvalue) ) then
        ! set default value
        write( IOBuffer, '("setting ",A," to default value ",G16.9)' ) trim(parameterqualifiers), defaultvalue
        call LogWrite
        parametervariable = defaultvalue
      else if ( .not. present(status) ) then
        ! Terminate with error, if error can not be returned through status
        call Error( "Could not find parameter <"//parameterqualifiers//">" )
        !return
      end if
    end if
    if ( present(status) ) status=stat

  end subroutine Global_FileReadParameter_RK


!==============================================================!
!  Subroutine Global_FileReadParameter_RKdim1                  !
!==============================================================!

  subroutine Global_FileReadParameter_RKdim1( parametervariable, iounit, parameterqualifiers, &
&                                            rewind_before, defaultvalue, status )

    implicit none

    ! Declare arguments
    real(RK), dimension(:), intent(out)          :: parametervariable
    integer, intent(in)                          :: iounit
    character(*), intent(in)                     :: parameterqualifiers
    logical, intent(in), optional                :: rewind_before
    real(RK), dimension(:), intent(in), optional :: defaultvalue
    integer, intent(out), optional               :: status

    ! Declare local variables
    character(IOBufferLength) :: buffer
    integer                   :: stat
    !character(FileNameLength) :: fn

    buffer = Global_FileReadParameter(iounit, parameterqualifiers, rewind_before, stat)
    if ( stat == 0 ) then
      read( buffer, * ) parametervariable

    else !if ( stat < 0 ) then
      ! parameter could not be read
      if ( present(defaultvalue) ) then
        ! set default value
        write( IOBuffer, '("setting ",A," to default value ")' ) trim(parameterqualifiers)
        call LogWrite
        write( IOBuffer, * ) defaultvalue
        call LogWrite
        parametervariable = defaultvalue
      else if ( .not. present(status) ) then
        ! Terminate with error, if error can not be returned through status
        call Error( "Could not find parameter <"//parameterqualifiers//">" )
        !return
      end if
    end if
    if ( present(status) ) status=stat

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
    write( iounit, '(A, T12, "=", A)' ) trim( parameter ), trim( IOBuffer )

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
    write( IOBuffer, '(72("-"))')
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
    integer(K4B), parameter :: IA=16807, IM=2147483647, IQ=127773, IR=2836

    ! Generate random number
    ix = ieor(ix, ishft(ix, 13))
    ix = ieor(ix, ishft(ix, -17))
    ix = ieor(ix, ishft(ix, 5))
    randk = iy / IQ
    iy = IA * (iy - randk * IQ) - IR * randk
    if( iy < 0 ) iy = iy + IM
    iharvest = 1 + ishft(int(range, RK) * ior(iand(IM, ieor(ix, iy)), 1), -31)

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
    integer(K4B), parameter :: IA=16807, IM=2147483647, IQ=127773, IR=2836

    ! Generate random number
    ix = ieor(ix, ishft(ix, 13))
    ix = ieor(ix, ishft(ix, -17))
    ix = ieor(ix, ishft(ix, 5))
    randk = iy / IQ
    iy = IA * (iy - randk * IQ) - IR * randk
    if( iy < 0 ) iy = iy + IM
    rharvest = l_range + am * ior(iand(IM,ieor(ix,iy)),1) * (h_range - l_range)

  end function Global_Rrnd


!==============================================================!
!  Function Global_String_Len_Trim                                 !
!==============================================================!

  pure function Global_String_Len_Trim( string, trim_left, trim_right ) result( length )

    !> Get options
    !> \param string     ... string to trim  character(*)
    !> \param trim_left  ... trim left? (default: .false.)  logical
    !> \param trim_right ... trim right? (default: .true.)  logical
    !> \return length    ... length of trimmed string

    implicit none

    ! Declare arguments
    character(*), intent(in) :: string
    logical,optional,intent(in) :: trim_left
    logical,optional,intent(in) :: trim_right

    ! Declare local variables
    logical :: do_trim_left, do_trim_right
    integer :: pos1, pos2

    ! Declare result
    integer :: length

    do_trim_left = .false.
    do_trim_right = .true.

    if( present(trim_left) ) do_trim_left = trim_left
    if( present(trim_right) ) do_trim_right = trim_right

    pos1 = 1
    pos2 = len(string)
    if( do_trim_right ) pos2 = verify(string,Whitespaces,.true.)
    if( do_trim_left ) pos1 = verify(string,Whitespaces)

    if( pos1/=0 .and. pos2/=0 ) then
      length = pos2-pos1+1
    else
      length = 0
    end if

  end function Global_String_Len_Trim


!==============================================================!
!  Function Global_String_TrimR                                !
!==============================================================!

  pure function Global_String_TrimR( string ) result( trimmed_string )

    !> Get options
    !> \param string          ... string to trim  character(*)
    !> \return trimmed_string ... trimmed string  character()

    implicit none

    ! Declare arguments
    character(*), intent(in) :: string

    ! Declare local variables
    integer :: pos2

    ! Declare result
    !character(len_trim(string)) :: trimmed_string
    character(strlen_trim(string,.false.,.true.)) :: trimmed_string

    pos2 = verify(string,Whitespaces,.true.)

    if( pos2/=0 ) then
      trimmed_string = trim(string(:pos2))
    else
      trimmed_string = ""
    end if

  end function Global_String_TrimR

!==============================================================!
!  Function Global_String_TrimL                                !
!==============================================================!

  pure function Global_String_TrimL( string ) result( trimmed_string )

    !> Get options
    !> \param string          ... string to trim  character(*)
    !> \return trimmed_string ... trimmed string  character()

    implicit none

    ! Declare arguments
    character(*), intent(in) :: string

    ! Declare local variables
    integer :: pos1

    ! Declare result
    !character(len(string)) :: trimmed_string
    character(strlen_trim(string,.true.,.false.)) :: trimmed_string

    pos1 = verify(string,Whitespaces)

    if( pos1/=0 ) then
      trimmed_string = trim(string(pos1:))
    else
      trimmed_string = ""
    end if

  end function Global_String_TrimL

!==============================================================!
!  Function Global_String_TrimLR                               !
!==============================================================!

  pure function Global_String_TrimLR( string ) result( trimmed_string )

    !> Get options
    !> \param string          ... string to trim   character(*)
    !> \return trimmed_string ... trimmed string   character()

    implicit none

    ! Declare arguments
    character(*), intent(in) :: string

    ! Declare local variables
    integer :: pos1, pos2

    ! Declare result
    !character(len(string)) :: trimmed_string
    character(strlen_trim(string,.true.,.true.)) :: trimmed_string

    pos1 = verify(string,Whitespaces)
    pos2 = verify(string,Whitespaces,.true.)

    if( pos1/=0 .and. pos2/=0 .and. pos2>=pos1 ) then
      trimmed_string = trim(string(pos1:pos2))
    else
      trimmed_string = ""
    end if

  end function Global_String_TrimLR


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

!     ! Include MPI header
! #if MPI_VER > 0
!     include 'mpif.h'
! #endif

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

#if MPI_VER > 0
    if( NProcs > 0 ) then
      ! original version 0: last process might get smaller range_size
      ! The if-statement reads: 
      ! only do it if we are in the equilibration phase of a MC  simulation
      ! and common equilibration is active. It is a little complicated, but that cannot be helped
      if( (SimulationType .ne. MonteCarlo) .or. (CommonEqui .and. (Equilibration .or. Step==0))) then 
        range_size = 1 + (overall_size - 1) / NProcs
        first_index = 1 + NProc * range_size
        last_index = min( first_index + range_size - 1, overall_size )
        range_size = last_index - first_index + 1
      else
        first_index=1
        last_index = overall_size
        range_size=overall_size
      endif
    
    else
      first_index=0
      last_index = -1
      range_size=0
    end if

#else
    first_index=1
    last_index = overall_size
    range_size=overall_size
#endif

  end function Global_GetProcRange

!==============================================================!
!  Subroutine Write Restart File when more writing time needed !
!==============================================================!

subroutine time_left(time_limit)

    implicit none

    ! could also use (an extended version of) TStopwatch
    
#if MPI_VER > 0
    ! Include MPI header
    include 'mpif.h'
#endif

    real(RK) :: time_remaining
    integer  :: time_limit
    
!     integer  :: ierror
#ifdef __INTEL_COMPILER
    integer  :: err
#endif
#ifdef KARLS
    character*10 string_max_time
#endif

    real(RK) :: time_elapsed    ! [sec]
    real(RK), save :: first_time
    logical, save :: FirstCAll =.TRUE.
    !integer :: time
    integer(4) :: sysclkcount, sysclkcountrate, sysclkcountmax

    if (FirstCAll)then
#if MPI_VER > 0
       first_time = MPI_WTIME()
!#elif defined ENABLE_OMP ! comment put by simon -> otherwise omp error
!       first_time = omp_get_wtime()   !-"-
#else
       !first_time = real(time())
       !!first_time = rtc()
       ! call system_clock(count_rate=sysclkcountrate,count_max=sysclkcountmax)
       ! call system_clock(sysclkcount)
       call system_clock(sysclkcount, sysclkcountrate, sysclkcountmax)
       first_time = real(real(sysclkcount)/sysclkcountrate)
#endif       
       FirstCall = .FALSE.
    end if
#if MPI_VER > 0
    time_elapsed = MPI_WTIME() - first_time
!#elif defined ENABLE_OMP   ! comment put by simon -> otherwise omp error
!      first_time = omp_get_wtime() - first_time        ! -"-
#else
    !time_elapsed = real(time()) - first_time
    call system_clock(sysclkcount, sysclkcountrate, sysclkcountmax)
    time_elapsed = real(sysclkcount)/sysclkcountrate - first_time
#endif       

! Get CPU time consumed by each task and compute the maximum value
!    call cpu_time(cputime)
! CPU time (!= elapsed wallclock time) does not make much sense here! There are also problems with multithreaded programs and "wrap around".

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

! Compute the remaining walltime
    time_remaining = max_time - real(time_elapsed)/60.

    if (time_remaining .le. time_limit) then
       write( IOBuffer, '("Simulation Abort due to Time Constraints on simulation cluster (time remaining=",G8.1,"<",G8.1," min)")' ) time_remaining, real(time_elapsed)/60.
       call LogWrite
       call LogWriteBlank

#ifdef __INTEL_COMPILER
         err = SetTerminateProgram( 1 )
#else
         call SetTerminateProgram
#endif
    !else
    !   write( IOBuffer, '("time remaining [min]: ",I5)' ) time_remaining
    !   call LogWrite
    end if

  end subroutine time_left


#ifdef USE_PRINTPROCSTATUS
!==============================================================!
!  Subroutine Global_printprocStatus
!==============================================================!

subroutine Global_printprocStatus(tag_string)
      
      implicit none
      
#if MPI_VER > 0
      ! Include MPI header
      include 'mpif.h'
#endif
      
      ! Declare arguments
      character(*), intent(in), optional :: tag_string
      
      ! Declare local variables
      character(*), parameter   :: procfilename = '/proc/self/status'
      character(IOBufferLength) :: buffer,token,valbuffer
      integer                   :: stat, linenr, seppos, i, val, nvalread
      integer, parameter        :: numvalues = 3
      character(*),parameter,dimension(numvalues) :: tokens = (/ &
&                                                       "VmPeak" &
&                                                      ,"VmSize" &
&                                                      ,"VmRSS " &
&                                                             /)
      integer(8) values(numvalues)
#if MPI_VER > 0
      integer(8) values_minmaxsum(numvalues,3)
#endif
      
      open(unit=iounit_proc, file=procfilename, action='read', iostat=stat)
      !if ( stat /= 0 ) return  ! dangerous for MPI version if not all ranks do exit...
      
      !linenr = 0
      nvalread = 0
      values=0
      
      do ! endless loop
        read(iounit_proc, '(A)', iostat=stat) buffer
        if (stat /= 0) exit  ! exit if nothing to read (EOF)
        !linenr = linenr + 1
        seppos=scan(buffer,': ')
        token=trim(adjustl(buffer(1:seppos-1)))
        valbuffer=trim(adjustl(buffer(seppos+1:)))
        do i=1,numvalues
          if ( trim(token) .eq. trim(tokens(i)) ) then
            read(valbuffer,*,iostat=stat)  values(i)
            !if (stat /= 0) continue
            nvalread = nvalread + 1
            if (index(valbuffer,'kB',.true.).gt.0) then
              values(i)=values(i)*1024
            end if
            if (nvalread .ge. numvalues) exit  ! exit if all data already read
          end if
        end do
      end do
      
      close(iounit_proc)
      
      call LogWriteBlank
      if( present( tag_string ) ) then
        write( IOBuffer, '("( ",A," ",A)' ) trim(procfilename), trim(tag_string)
      else
        write( IOBuffer, '( "(",A)' ) trim(procfilename)
      end if
      call LogWriteTime
#if MPI_VER > 0
      values_minmaxsum(:,1)=-values
      values_minmaxsum(:,2)=values
      values_minmaxsum(:,3)=values
      if ( RootProc_W ) then
        call MPI_Reduce( MPI_IN_PLACE, values_minmaxsum, numvalues*2, MPI_INTEGER8, MPI_MAX, NRootProc_W, MPI_COMM_WORLD, ierror )
      else
        call MPI_Reduce( values_minmaxsum, values_minmaxsum, numvalues*2, MPI_INTEGER8, MPI_MAX, NRootProc_W, MPI_COMM_WORLD, ierror )
      end if
      call MPI_Reduce( values, values_minmaxsum(:,3), numvalues, MPI_INTEGER8, MPI_SUM, NRootProc_W, MPI_COMM_WORLD, ierror )
      values_minmaxsum(:,1)=-values_minmaxsum(:,1)
      do i=1,numvalues
        write( IOBuffer, '(" ",A,"(min max sum):",3I20)' ) trim(tokens(i)),values_minmaxsum(i,:)
        call LogWrite
      end do
#else
      do i=1,numvalues
        write( IOBuffer, '(" ",A,":",I16)' ) trim(tokens(i)),values(i)
        call LogWrite
      end do
#endif
      write( IOBuffer, '(") ",A)' ) trim(procfilename)
      call LogWriteTime
      call LogWriteBlank

  end subroutine Global_printprocStatus
#endif


end module ms2_global
