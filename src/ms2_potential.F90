!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 3.0                !
!  (c) 2017 by TU Kaiserslautern / U Paderborn                 !
!      P.O. Box 67653                                          !
!      67653 Kaiserslautern                                    !
!==============================================================!
!  Module ms2_potential                                        !
!  Contains TPot* objects                                      !
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
!DEC$ MESSAGE:'Compiling ms2_potential.F90...'
#endif

module ms2_potential

  use ms2_molecule
  use ms2_site


!==============================================================!
!  Type TPotMIEnmMIEnm                                         !
!==============================================================!

  type TPotMIEnmMIEnm

    type(TSiteMIEnm), pointer :: Site1, Site2
    real(RK)                  :: Sigma, Epsilon
    real(RK)                  :: Mie_n, Mie_m, Mie_a, Mie_nHalf, Mie_mHalf
    real(RK)                  :: RCutoffSquared, RCutoffSquaredScaled
    real(RK)                  :: EPotCorr, VirialCorr, d2EpotdV2Corr,EPotTestCorr
    logical                   :: SameComponent
    real(RK)                  :: SigmaSquared
    real(RK)                  :: EpsilonMie_a, EpsilonMie_aF
    real(RK)                  :: BoxlengthInv, BoxLengthThird
    integer, pointer, contiguous          :: NInCutoff(:), CutoffPartner(:, :)
    integer, pointer, contiguous          :: RDFSum(:)
#if OSMOP == 2
    real(RK), pointer, contiguous         :: VirialProfile(:)
#endif


  end type TPotMIEnmMIEnm

  interface Construct
    module procedure TPotMIEMIE_Construct
  end interface

  interface Destruct
    module procedure TPotMIEMIE_Destruct
  end interface

  interface Force
    module procedure TPotMIEMIE_Force
  end interface

  interface Get_RDF
    module procedure TPotMIEMIE_RDF
  end interface
  
  interface Force_Trans
    module procedure TPotMIEMIE_Force_Trans
  end interface

  interface ChemicalPotential
    module procedure TPotMIEMIE_ChemicalPotential
  end interface

  interface UpdateBoxLength
    module procedure TPotMIEMIE_UpdateBoxLength
  end interface



!==============================================================!
!  Type TPotChargeCharge                                       !
!==============================================================!

  type TPotChargeCharge

    type(TSiteCharge), pointer :: Site1, Site2
    real(RK)                   :: Epsilon
    real(RK)                   :: RShieldSquared
    real(RK)                   :: RCutoffSquared
    real(RK)                   :: RFConstant
    logical                    :: SameComponent
    integer, pointer, contiguous           :: NInCutoff(:), CutoffPartner(:, :)
#if OSMOP == 2
    real(RK), pointer, contiguous         :: VirialProfile(:)
#endif

  end type TPotChargeCharge

  interface Construct
    module procedure TPotCC_Construct
  end interface

  interface Destruct
    module procedure TPotCC_Destruct
  end interface

  interface Force
    module procedure TPotCC_Force
    module procedure TPotCC_Force_Ewald
  end interface

  interface Force_Trans
    module procedure TPotCC_Force_Trans
    module procedure TPotCC_Force_Ewald_Trans
  end interface

  interface ChemicalPotential
    module procedure TPotCC_ChemicalPotential
  end interface

  interface ErrorApprox
    module procedure TPoterfc_approx
  end interface



!==============================================================!
!  Type TPotChargeDipole                                       !
!==============================================================!

  type TPotChargeDipole

    type(TSiteCharge), pointer :: Site1
    type(TSiteDipole), pointer :: Site2
    real(RK)                   :: Epsilon
    real(RK)                   :: RShieldSquared
    real(RK)                   :: RCutoffSquared
    logical                    :: SameComponent
    integer, pointer, contiguous           :: NInCutoff(:), CutoffPartner(:, :)
#if OSMOP == 2
    real(RK), pointer, contiguous          :: VirialProfile(:)
#endif

  end type TPotChargeDipole

  interface Construct
    module procedure TPotCD_Construct
  end interface

  interface Destruct
    module procedure TPotCD_Destruct
  end interface

  interface Force
    module procedure TPotCD_Force
  end interface

  interface Force_Trans
    module procedure TPotCD_Force_Trans
  end interface

  interface ChemicalPotential
    module procedure TPotCD_ChemicalPotential
  end interface



!==============================================================!
!  Type TPotChargeQuadrupole                                   !
!==============================================================!

  type TPotChargeQuadrupole

    type(TSiteCharge), pointer     :: Site1
    type(TSiteQuadrupole), pointer :: Site2
    real(RK)                       :: Epsilon
    real(RK)                       :: RShieldSquared
    real(RK)                       :: RCutoffSquared
    logical                        :: SameComponent
    integer, pointer, contiguous               :: NInCutoff(:), CutoffPartner(:, :)
#if OSMOP == 2
    real(RK), pointer, contiguous              :: VirialProfile(:)
#endif

  end type TPotChargeQuadrupole

  interface Construct
    module procedure TPotCQ_Construct
  end interface

  interface Destruct
    module procedure TPotCQ_Destruct
  end interface

  interface Force
    module procedure TPotCQ_Force
  end interface

  interface Force_Trans
    module procedure TPotCQ_Force_Trans
  end interface

  interface ChemicalPotential
    module procedure TPotCQ_ChemicalPotential
  end interface


!==============================================================!
!  Type TPotDipoleCharge                                       !
!==============================================================!

  type TPotDipoleCharge

    type(TSiteDipole), pointer :: Site1
    type(TSiteCharge), pointer :: Site2
    real(RK)                   :: Epsilon
    real(RK)                   :: RShieldSquared
    real(RK)                   :: RCutoffSquared
    logical                    :: SameComponent
    integer, pointer, contiguous           :: NInCutoff(:), CutoffPartner(:, :)
#if OSMOP == 2
    real(RK), pointer, contiguous          :: VirialProfile(:)
#endif

  end type TPotDipoleCharge

  interface Construct
    module procedure TPotDC_Construct
  end interface

  interface Destruct
    module procedure TPotDC_Destruct
  end interface

  interface Force
    module procedure TPotDC_Force
  end interface

  interface Force_Trans
    module procedure TPotDC_Force_Trans
  end interface

  interface ChemicalPotential
    module procedure TPotDC_ChemicalPotential
  end interface


!==============================================================!
!  Type TPotDipoleDipole                                       !
!==============================================================!

  type TPotDipoleDipole

    type(TSiteDipole), pointer :: Site1, Site2
    real(RK)                   :: Epsilon
    real(RK)                   :: RCutoffSquared
    real(RK)                   :: RShieldSquared
    real(RK)                   :: RFConstant
    logical                    :: SameComponent
    integer, pointer, contiguous           :: NInCutoff(:), CutoffPartner(:, :)
#if OSMOP == 2
    real(RK), pointer, contiguous          :: VirialProfile(:)
#endif

  end type TPotDipoleDipole

  interface Construct
    module procedure TPotDD_Construct
  end interface

  interface Destruct
    module procedure TPotDD_Destruct
  end interface

  interface Force
    module procedure TPotDD_Force
  end interface

  interface Force_Trans
    module procedure TPotDD_Force_Trans
  end interface

  interface ChemicalPotential
    module procedure TPotDD_ChemicalPotential
  end interface


!==============================================================!
!  Type TPotDipoleQuadrupole                                   !
!==============================================================!

  type TPotDipoleQuadrupole

    type(TSiteDipole), pointer     :: Site1
    type(TSiteQuadrupole), pointer :: Site2
    real(RK)                       :: Epsilon
    real(RK)                       :: RCutoffSquared
    real(RK)                       :: RShieldSquared
    logical                        :: SameComponent
    integer, pointer, contiguous               :: NInCutoff(:), CutoffPartner(:, :)
#if OSMOP == 2
    real(RK), pointer, contiguous              :: VirialProfile(:)
#endif

  end type TPotDipoleQuadrupole

  interface Construct
    module procedure TPotDQ_Construct
  end interface

  interface Destruct
    module procedure TPotDQ_Destruct
  end interface

  interface Force
    module procedure TPotDQ_Force
  end interface

  interface Force_Trans
    module procedure TPotDQ_Force_Trans
  end interface

  interface ChemicalPotential
    module procedure TPotDQ_ChemicalPotential
  end interface

!==============================================================!
!  Type TPotQuadrupoleCharge                                   !
!==============================================================!

  type TPotQuadrupoleCharge

    type(TSiteQuadrupole), pointer :: Site1
    type(TSiteCharge), pointer     :: Site2
    real(RK)                       :: Epsilon
    real(RK)                       :: RShieldSquared
    real(RK)                       :: RCutoffSquared
    logical                        :: SameComponent
    integer, pointer, contiguous               :: NInCutoff(:), CutoffPartner(:, :)
#if OSMOP == 2
    real(RK), pointer, contiguous              :: VirialProfile(:)
#endif

  end type TPotQuadrupoleCharge

  interface Construct
    module procedure TPotQC_Construct
  end interface

  interface Destruct
    module procedure TPotQC_Destruct
  end interface

  interface Force
    module procedure TPotQC_Force
  end interface

  interface Force_Trans
    module procedure TPotQC_Force_Trans
  end interface

  interface ChemicalPotential
    module procedure TPotQC_ChemicalPotential
  end interface


!==============================================================!
!  Type TPotQuadrupoleDipole                                   !
!==============================================================!

  type TPotQuadrupoleDipole

    type(TSiteQuadrupole), pointer :: Site1
    type(TSiteDipole), pointer     :: Site2
    real(RK)                       :: Epsilon
    real(RK)                       :: RCutoffSquared
    real(RK)                       :: RShieldSquared
    logical                        :: SameComponent
    integer, pointer, contiguous               :: NInCutoff(:), CutoffPartner(:, :)
#if OSMOP == 2
    real(RK), pointer, contiguous              :: VirialProfile(:)
#endif

  end type TPotQuadrupoleDipole

  interface Construct
    module procedure TPotQD_Construct
  end interface

  interface Destruct
    module procedure TPotQD_Destruct
  end interface

  interface Force
    module procedure TPotQD_Force
  end interface

  interface Force_Trans
    module procedure TPotQD_Force_Trans
  end interface

  interface ChemicalPotential
    module procedure TPotQD_ChemicalPotential
  end interface


!==============================================================!
!  Type TPotQuadrupoleQuadrupole                               !
!==============================================================!

  type TPotQuadrupoleQuadrupole

    type(TSiteQuadrupole), pointer :: Site1, Site2
    real(RK)                       :: Epsilon
    real(RK)                       :: RCutoffSquared
    real(RK)                       :: RShieldSquared
    logical                        :: SameComponent
    integer, pointer, contiguous               :: NInCutoff(:), CutoffPartner(:, :)
#if OSMOP == 2
    real(RK), pointer, contiguous              :: VirialProfile(:)
#endif

  end type TPotQuadrupoleQuadrupole

  interface Construct
    module procedure TPotQQ_Construct
  end interface

  interface Destruct
    module procedure TPotQQ_Destruct
  end interface

  interface Force
    module procedure TPotQQ_Force
  end interface

  interface Force_Trans
    module procedure TPotQQ_Force_Trans
  end interface

  interface ChemicalPotential
    module procedure TPotQQ_ChemicalPotential
  end interface



contains



!==============================================================!
!  Subroutine TPotMIEMIE_Construct                               !
!==============================================================!

  subroutine TPotMIEMIE_Construct( this, i1, i2, j1, j2, Molecule1, Molecule2, RCutoff, ScaleSigma, ScaleEpsilon )

    implicit none

    ! Declare arguments
    type(TPotMIEnmMIEnm)        :: this
    integer, intent(in)         :: i1, i2, j1, j2
    type(TMolecule), intent(in) :: Molecule1, Molecule2
    real(RK), intent(in)        :: RCutoff
    real(RK), intent(in)        :: ScaleSigma, ScaleEpsilon

    ! Declare local variables
    real(RK) :: RCutoff3Inv, RCutoff9Inv
    real(RK) :: tau, tau1, tau2
    real(RK) :: Pi2mie_a, Piminus23mie_a, Pi29mie_a

        
    ! Construct potential
    this%Site1 => Molecule1%SiteMIEnm(j1)
    this%Site2 => Molecule2%SiteMIEnm(j2)
    this%SameComponent = i1 == i2
    this%Sigma = .5_RK * (this%Site1%sig + this%Site2%sig)
    this%Epsilon = sqrt(this%Site1%eps * this%Site2%eps)
    
    
    ! Calculate parameter for MIE-Potential -> R. Fingerhut -> Note: define mix rule for n and m
    !this%Mie_n = .5_RK * (this%Site1%mie_n + this%Site2%mie_n) !Später so einführen und Anpassungsparameter wie für sig und eps bei Mischung
    !this%Mie_m = .5_RK * (this%Site1%mie_m + this%Site2%mie_m)
    this%Mie_n = 3._RK+((this%Site1%mie_n-3)*(this%Site2%mie_n-3))**.5_RK !Von Erich-> Imperial College
    this%Mie_m = 3._RK+((this%Site1%mie_m-3)*(this%Site2%mie_m-3))**.5_RK
    
    ! prefactor in the mie-function: Mie_a=f(n,m)
    this%Mie_a = 1._RK / (this%Mie_n-this%Mie_m) * (this%Mie_n**this%Mie_n/(this%Mie_m**this%Mie_m))**(1._RK/(this%Mie_n-this%Mie_m))
    
    this%Mie_nHalf = .5_RK * this%Mie_n
    this%Mie_mHalf = .5_RK * this%Mie_m
    
    Pi2mie_a = Pi * 2._RK * this%Mie_a
    Piminus23mie_a = Pi * this%Mie_a * (-2._RK)/3._RK
    Pi29mie_a = Pi * this%Mie_a * (2._RK/9._RK)
    

    if( .not. this%SameComponent ) then
      this%Sigma = this%Sigma * ScaleSigma
      this%Epsilon = this%Epsilon * ScaleEpsilon
    end if
    this%EpsilonMie_a = this%Mie_a * this%Epsilon

    ! Calculate long-range corrections
    this%RCutoffSquared = RCutoff**2
    tau1 = sqrt( sum( this%Site1%r(:)**2 ))
    tau2 = sqrt( sum( this%Site2%r(:)**2 ))
    tau = max( tau1, tau2 )

    if( (CutoffMode .eq. CenterofMass) .and. (tau > 1E-10_RK) ) then
      if( (tau1 > 1E-10_RK) .and. (tau2 > 1E-10_RK) ) then

        this%EPotCorr = Pi2mie_a * this%Epsilon * &
&         ( TISSu(-this%Mie_nHalf, RCutoff, this%Sigma**2._RK, tau1, tau2) &
&         - TISSu(-this%Mie_mHalf, RCutoff, this%Sigma**2._RK, tau1, tau2) )

        this%VirialCorr = Piminus23mie_a * this%Epsilon * &
&         ( TISSp(-this%Mie_nHalf, RCutoff, this%Sigma**2._RK, tau1, tau2) &
&         - TISSp(-this%Mie_mHalf, RCutoff, this%Sigma**2._RK, tau1, tau2) )

        this%d2EpotdV2Corr = Pi29mie_a * this%Epsilon * &
&         ( TISSd2EpotdV2(-this%Mie_nHalf, RCutoff, this%Sigma**2._RK, tau1, tau2) &
&         - TISSd2EpotdV2(-this%Mie_mHalf, RCutoff, this%Sigma**2._RK, tau1, tau2) )

      else
        this%EPotCorr = Pi2mie_a * this%Epsilon * &
&         ( TICSu(-this%Mie_nHalf, RCutoff, this%Sigma**2._RK, tau) &
&         - TICSu(-this%Mie_mHalf, RCutoff, this%Sigma**2._RK, tau) )

        this%VirialCorr = Piminus23mie_a * this%Epsilon * &
&         ( TICSp(-this%Mie_nHalf, RCutoff, this%Sigma**2._RK, tau) &
&         - TICSp(-this%Mie_mHalf, RCutoff, this%Sigma**2._RK, tau) )

        this%d2EpotdV2Corr = Pi29mie_a * this%Epsilon * &
&         ( TICSd2EpotdV2(-this%Mie_nHalf, RCutoff, this%Sigma**2._RK, tau) &
&         - TICSd2EpotdV2(-this%Mie_mHalf, RCutoff, this%Sigma**2._RK, tau) )

      endif
    else ! Site-site cutoff or both sites in center of mass
     this%EPotCorr = Pi2mie_a * this%Epsilon * ( TICCu(-this%Mie_nHalf, RCutoff, this%Sigma**2._RK) - TICCu(-this%Mie_mHalf, RCutoff, this%Sigma**2._RK) )

     this%VirialCorr = Piminus23mie_a * this%Epsilon * ( TICCp(-this%Mie_nHalf, RCutoff, this%Sigma**2._RK) - TICCp(-this%Mie_mHalf, RCutoff, this%Sigma**2._RK) )

     this%d2EpotdV2Corr = Pi29mie_a * this%Epsilon *  ( TICCd2EpotdV2(-this%Mie_nHalf, RCutoff, this%Sigma**2._RK) - TICCd2EpotdV2(-this%Mie_mHalf, RCutoff, this%Sigma**2._RK) )

    end if
    this%EPotTestCorr = 2._RK * this%EPotCorr

 
  contains

    real(RK) function TISSu( n, rc, sigma2, tau1, tau2 )

      ! Declare arguments
      real(RK), intent(in) :: n
      real(RK), intent(in) :: rc, sigma2, tau1, tau2

      ! Declare local variables
      real(RK) :: tauPlus, tauMinus

      tauPlus = tau1 + tau2
      tauMinus = abs( tau1 - tau2 )

      ! Calculate angle averaged partial integral
      TISSu = - ( (rc+tauPlus)**(2._RK*n+4._RK) - (rc+tauMinus)**(2._RK*n+4._RK) &
&             - (rc-tauMinus)**(2._RK*n+4._RK) + (rc-tauPlus)**(2._RK*n+4._RK) ) * rc &
&             / ( 8._RK * sigma2**n * tau1 * tau2 &
&             * (n+1._RK) * (2._RK*n+3._RK) * (2._RK*n+4._RK) ) &
&             + ( (rc+tauPlus)**(2._RK*n+5._RK) - (rc+tauMinus)**(2._RK*n+5._RK) &
&             - (rc-tauMinus)**(2._RK*n+5._RK) + (rc-tauPlus)**(2._RK*n+5._RK) ) &
&             / ( 8._RK * sigma2**n * tau1 * tau2 &
&             * (n+1._RK) * (2._RK*n+3._RK) * (2._RK*n+4._RK) * (2._RK*n+5._RK) )

    end function TISSu



    real(RK) function TICSu( n, rc, sigma2, tau )

      ! Declare arguments
      real(RK), intent(in) :: n
      real(RK), intent(in) :: rc, sigma2, tau

      ! Calculate angle averaged partial integral
      TICSu = - ( (rc+tau)**(2._RK*n+3._RK) - (rc-tau)**(2._RK*n+3._RK) ) * rc &
&             / ( 4._RK * sigma2**n * tau * (n+1._RK) * (2._RK*n+3._RK) ) &
&             + ( (rc+tau)**(2._RK*n+4._RK) - (rc-tau)**(2._RK*n+4._RK) ) &
&             / ( 4._RK * sigma2**n * tau * (n+1._RK) * (2._RK*n+3._RK) * (2._RK*n+4._RK) )

    end function TICSu


        
    real(RK) function TICCu( n, rc, sigma2 )

      ! Declare arguments
      real(RK), intent(in) :: n
      real(RK), intent(in) :: rc, sigma2

      ! Calculate angle averaged partial integral
      TICCu = - ( rc**(2._RK*n+3._RK) ) / ( sigma2**n * (2._RK*n+3._RK) )

    end function TICCu



    real(RK) function TISSp( n, rc, sigma2, tau1, tau2 )

      ! Declare arguments
      real(RK), intent(in) :: n
      real(RK), intent(in) :: rc, sigma2, tau1, tau2

      ! Declare local variables
      real(RK) :: tauPlus, tauMinus

      tauPlus = tau1 + tau2
      tauMinus = abs( tau1 - tau2 )

      ! Calculate angle averaged partial integral
      TISSp = - ( (rc+tauPlus)**(2._RK*n+3._RK) - (rc+tauMinus)**(2._RK*n+3._RK) &
&             - (rc-tauMinus)**(2._RK*n+3._RK) + (rc-tauPlus)**(2._RK*n+3._RK) ) * rc**2._RK &
&             / ( 8._RK * sigma2**n * tau1 * tau2 * (n+1._RK) * (2._RK*n+3._RK) ) &
&             - 3._RK * TISSu(n,rc,sigma2,tau1,tau2)

    end function TISSp



    real(RK) function TICSp( n, rc, sigma2, tau )

      ! Declare arguments
      real(RK), intent(in) :: n
      real(RK), intent(in) :: rc, sigma2, tau

      ! Calculate angle averaged partial integral
      TICSp = - ( (rc+tau)**(2._RK*n+2._RK) - (rc-tau)**(2._RK*n+2._RK) ) * rc**2._RK &
&             / ( 4._RK * sigma2**n * tau * (n+1._RK) ) &
&             - 3._RK * TICSu(n,rc,sigma2,tau)

    end function TICSp


    real(RK) function TICCp( n, rc, sigma2 )

      ! Declare arguments
      real(RK), intent(in) :: n
      real(RK), intent(in) :: rc, sigma2

      ! Calculate angle averaged partial integral
      TICCp = 2._RK*n*TICCu(n,rc,sigma2)

    end function TICCp


    real(RK) function TISSd2EpotdV2( n, rc, sigma2, tau1, tau2 )

      ! Declare arguments
      real(RK), intent(in) :: n
      real(RK), intent(in) :: rc, sigma2, tau1, tau2

      ! Declare local variables
      real(RK) :: tauPlus, tauMinus, A, B, C, D, AN, BN, CN, DN

      tauPlus = tau1 + tau2
      tauMinus = abs( tau1 - tau2 )

      ! Calculate angle averaged partial integral
      A = 18._RK*n**2._RK*rc**3._RK          -6._RK *n**2._RK*rc**2._RK*tauPlus&
&        -3._RK *tauPlus**3._RK              +6._RK *n*tauPlus**2._RK*rc&
&        +4._RK *rc**3._RK*n**3._RK          +12._RK*rc**3._RK&
&        +26._RK*n*rc**3._RK                 -9._RK *rc**2._RK*tauPlus&
&        +6._RK *tauPlus**2._RK*rc           -15._RK*rc**2._RK*n*tauPlus

      B = 18._RK*n**2._RK*rc**3._RK          +6._RK *n**2._RK*rc**2._RK*tauPlus&
&        +3._RK *tauPlus**3._RK              +6._RK *n*tauPlus**2._RK*rc&
&        +4._RK *rc**3._RK*n**3._RK          +12._RK*rc**3._RK&
&        +26._RK*n*rc**3._RK                 +9._RK *rc**2._RK*tauPlus&
&        +6._RK *tauPlus**2._RK*rc           +15._RK*rc**2._RK*n*tauPlus

      C = 18._RK*n**2._RK*rc**3._RK          -15._RK*rc**2._RK*tauMinus*n&
&        +4._RK *rc**3._RK*n**3._RK          +12._RK*rc**3._RK&
&        -6._RK *n**2._RK*rc**2._RK*tauMinus +6._RK *n*tauMinus**2._RK*rc&
&        -3._RK *tauMinus**3._RK             +26._RK*n*rc**3._RK&
&        +6._RK *tauMinus**2._RK*rc          -9._RK *rc**2._RK*tauMinus

      D = 18._RK*n**2._RK*rc**3._RK          +15._RK*rc**2._RK*tauMinus*n&
&        +4._RK *rc**3._RK*n**3._RK          +12._RK*rc**3._RK&
&        +6._RK *n**2._RK*rc**2._RK*tauMinus +6._RK *n*tauMinus**2._RK*rc&
&        +3._RK *tauMinus**3._RK             +26._RK*n*rc**3._RK&
&        +6._RK *tauMinus**2._RK*rc          +9._RK *rc**2._RK*tauMinus

      AN = (1._RK/8._RK)*(rc+tauPlus)**(2._RK*n+3._RK)/(tau1*tau2*(n+1._RK)*(2._RK*n+3._RK)*(rc+tauPlus)*(5._RK+2._RK*n)*(2._RK+n))

      BN = (1._RK/8._RK)*(rc-tauPlus)**(2._RK*n+3._RK)/(tau1*tau2*(n+1._RK)*(2._RK*n+3._RK)*(rc-tauPlus)*(5._RK+2._RK*n)*(2._RK+n))

      CN =-(1._RK/8._RK)*(rc+tauMinus)**(2._RK*n+3._RK)/(tau1*tau2*(n+1._RK)*(2._RK*n+3._RK)*(rc+tauMinus)*(5._RK+2._RK*n)*(2._RK+n))

      DN =-(1._RK/8._RK)*(rc-tauMinus)**(2._RK*n+3._RK)/(tau1*tau2*(n+1._RK)*(2._RK*n+3._RK)*(rc-tauMinus)*(5._RK+2._RK*n)*(2._RK+n))

      TISSd2EpotdV2 =-(A*AN+B*BN+C*CN+D*DN)/sigma2**n -2._RK * TISSp(n,rc,sigma2,tau1,tau2)

    end function TISSd2EpotdV2


    real(RK) function TICSd2EpotdV2( n, rc, sigma2, tau )

      ! Declare arguments
      real(RK), intent(in) :: n
      real(RK), intent(in) :: rc, sigma2, tau

      ! Declare local variables
      real(RK) :: A, B, AN, BN

      ! Calculate angle averaged partial integral

      A = 3._RK *tau**3._RK             +3._RK *rc**2._RK*tau&
&        +12._RK*n**2._RK*rc**3._RK     +11._RK*n*rc**3._RK&
&        +6._RK *n**2._RK*rc**2._RK*tau +9._RK *n*rc**2._RK*tau&
&        +3._RK *rc**3._RK              +4._RK *n**3._RK*rc**3._RK&
&        +3._RK *rc*tau**2._RK          +6._RK *rc*n*tau**2._RK

      B =-3._RK *tau**3._RK             -3._RK *rc**2._RK*tau&
&        +12._RK*n**2._RK*rc**3._RK     +11._RK*n*rc**3._RK&
&        -6._RK *n**2._RK*rc**2._RK*tau -9._RK *n*rc**2._RK*tau&
&        +3._RK *rc**3._RK              +4._RK *n**3._RK*rc**3._RK&
&        +3._RK *rc*tau**2._RK          +6._RK *rc*n*tau**2._RK

      AN = -(1._RK/4._RK)*(rc-tau)**(2._RK*n+2._RK)/(tau*(n+1._RK)*(rc-tau)*(2._RK+n)*(3._RK+2._RK*n))

      BN =  (1._RK/4._RK)*(rc+tau)**(2._RK*n+2._RK)/(tau*(n+1._RK)*(rc+tau)*(2._RK+n)*(3._RK+2._RK*n))

      TICSd2EpotdV2 =-(A*AN+B*BN)/sigma2**n - 2._RK * TICSp(n,rc,sigma2,tau)

    end function TICSd2EpotdV2


    real(RK) function TICCd2EpotdV2( n, rc, sigma2 )

      ! Declare arguments
      real(RK), intent(in) :: n
      real(RK), intent(in) :: rc, sigma2

      ! Declare local variables
      real(RK) :: A, AN

      A = rc**(3._RK+2._RK*n)

      AN = 2._RK*n*(2._RK*n-1._RK)/(3._RK+2._RK*n)

      ! Calculate angle averaged partial integral
      TICCd2EpotdV2 =-(A*AN)/sigma2**n

    end function TICCd2EpotdV2 


    real(RK) function TISSuAbl( n, rc, sigma2, tau1, tau2 )

      ! Declare arguments
      real(RK), intent(in) :: n
      real(RK), intent(in) :: rc, sigma2, tau1, tau2

      ! Declare local variables
      real(RK) :: tauPlus, tauMinus

      tauPlus = tau1 + tau2
      tauMinus = abs( tau1 - tau2 )

      ! Calculate angle averaged partial integral
      TISSuAbl = - ( (rc+tauPlus)**(2._RK*n+4._RK) - (rc+tauMinus)**(2._RK*n+4._RK) &
&                - (rc-tauMinus)**(2._RK*n+4._RK) + (rc-tauPlus)**(2._RK*n+4._RK) ) * rc * (-2._RK*n) &
&                / ( 8._RK * sigma2**(n+1._RK) * tau1 * tau2 &
&                * (n+1._RK) * (2._RK*n+3._RK) * (2._RK*n+4._RK) ) &
&                + ( (rc+tauPlus)**(2._RK*n+5._RK) - (rc+tauMinus)**(2._RK*n+5._RK) &
&                - (rc-tauMinus)**(2._RK*n+5._RK) + (rc-tauPlus)**(2._RK*n+5._RK) ) * (-2._RK*n)&
&                / ( 8._RK * sigma2**(n+1._RK) * tau1 * tau2 &
&                * (n+1._RK) * (2._RK*n+3._RK) * (2._RK*n+4._RK) * (2._RK*n+5._RK) )

    end function TISSuAbl

    real(RK) function TISSpAbl( n, rc, sigma2, tau1, tau2 )

      ! Declare arguments
      real(RK), intent(in) :: n
      real(RK), intent(in) :: rc, sigma2, tau1, tau2

      ! Declare local variables
      real(RK) :: tauPlus, tauMinus

      tauPlus = tau1 + tau2
      tauMinus = abs( tau1 - tau2 )

      ! Calculate angle averaged partial integral
      TISSpAbl = - ( (rc+tauPlus)**(2._RK*n+3._RK) - (rc+tauMinus)**(2._RK*n+3._RK) &
&                - (rc-tauMinus)**(2._RK*n+3._RK) + (rc-tauPlus)**(2._RK*n+3._RK) ) *rc**2._RK * (-2._RK*n)&
&                / ( 8._RK * sigma2**(n+1._RK) * tau1 * tau2 * (n+1._RK) * (2._RK*n+3._RK) ) &
&                - 3._RK * TISSuAbl(n,rc,sigma2,tau1,tau2)

    end function TISSpAbl


  end subroutine TPotMIEMIE_Construct



!==============================================================!
!  Subroutine TPotMIEMIE_Destruct                                !
!==============================================================!

  subroutine TPotMIEMIE_Destruct( this )

    implicit none

    ! Declare arguments
    type(TPotMIEnmMIEnm) :: this

    ! Destroy potential
    continue

  end subroutine TPotMIEMIE_Destruct

!==============================================================!
!  Subroutine TPotMIEMIE_Force                                   !
!==============================================================!

  subroutine TPotMIEMIE_Force( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotMIEnmMIEnm)     :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK)          :: SigmaSquared
    real(RK)          :: EpsilonMie_a, EpsilonMie_aF
    real(RK)          :: Mie_n, Mie_m, Mie_n1, Mie_m1, Mie_nHalf, Mie_mHalf, Mie_nRijMie_n, Mie_mRijMie_m
    real(RK)          :: RCutoffSquared
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: FXij, FYij, FZij, Fij
    real(RK)          :: RijSquared, RijSquaredInv, RijMie_nInv, RijMie_mInv
    real(RK)          :: EPotLocal, EPotLocal1, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr
    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
#if MPI_VER > 0
    integer           :: i0, N1, N2, ji
    logical           :: EvenN
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif
        

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ


    ! Assign local variables
    SameComponent = this%SameComponent
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK
    SigmaSquared = this%SigmaSquared
    EpsilonMie_a = this%EpsilonMie_a
    EpsilonMie_aF = this%EpsilonMie_aF
    Mie_n = this%Mie_n
    Mie_m = this%Mie_m
    Mie_n1 = Mie_n+1._RK
    Mie_m1 = Mie_m+1._RK
    Mie_nHalf = this%Mie_nHalf
    Mie_mHalf = this%Mie_mHalf
    RCutoffSquared = this%RCutoffSquaredScaled
    
#if MPI_VER > 0
    N1 = this%Site2%NPart
    N2 = N1 / 2
    EvenN = mod( N1, 2 ) == 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
    j1 = 0
    ji = 0
#else
    i1 = this%Site1%NPart
    j1 = this%Site2%NPart
#endif

!$OMP PARALLEL DEFAULT(SHARED) & 
#if MPI_VER > 0
!$OMP FIRSTPRIVATE(i0, N1, N2, ji, EvenN) &
#endif
!$OMP FIRSTPRIVATE(i1, j1) &
!$OMP PRIVATE( i, j, k, j0) &
!$OMP PRIVATE(sitecorr, EPotLocal1) &
!$OMP PRIVATE(RXi, RYi, RZi,  PXi, PYi, PZi,  FXi, FYi, FZi) &
!$OMP PRIVATE(RXij, RYij, RZij, PXij, PYij, PZij) &
!$OMP PRIVATE(FXij, FYij, FZij, Fij, RijSquared, RijSquaredInv, RijMie_nInv, RijMie_mInv ) 

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
!CDIR NODEP
#if OSMOP == 2
loop0:  do m=1,NBinsDen
          if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PXi < real(m)/NBinsDen-0.5_RK) then
              Bin1=m
              exit loop0
            end if
          end if
        end do loop0
#endif
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = RXij - anint( PXij )
          RYij = RYij - anint( PYij )
          RZij = RZij - anint( PZij )
          PXij = PXij - anint( PXij )
          PYij = PYij - anint( PYij )
          PZij = PZij - anint( PZij )
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          RijSquaredInv = SigmaSquared / RijSquared
          RijMie_nInv = RijSquaredInv**Mie_nHalf
          RijMie_mInv = RijSquaredInv**Mie_mHalf
          Mie_nRijMie_n = Mie_n * RijMie_nInv
          Mie_mRijMie_m = Mie_m * RijMie_mInv
          EPotLocal1 = RijMie_nInv - RijMie_mInv
          EPotLocal = EPotLocal + EPotLocal1
          Fij = EpsilonMie_aF * (Mie_nRijMie_n - Mie_mRijMie_m) * RijSquaredInv
          FXij = Fij * RXij
          FYij = Fij * RYij
          FZij = Fij * RZij
          VirialLocal = VirialLocal + (PXij * FXij + PYij * FYij + PZij * FZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
                Bin2=m 
                exit loop2
              end if
            end if
          end do loop2
          tempMin = min(Bin1, Bin2)
          tempMax = max(Bin1, Bin2)
          if(abs(PXij) .le. 0.5_RK) then
             VirialPart = (PXij * FXij + PYij * FYij + PZij * FZij)/(tempMax-tempMin+1._RK) 
             do m = tempMin, tempMax
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
             end do
          else
             VirialPart = (PXij * FXij + PYij * FYij + PZij * FZij)/(NBinsDen-tempMax+tempMin+1._RK) 
             do m = 1, tempMin
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
             end do
             do m = tempMax, NBinsDen
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
             end do
          end if
#endif
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)/RijSquared
          d2EpotdV2Local = d2EpotdV2Local + EpsilonMie_a * Ninth * ((Mie_nRijMie_n - Mie_mRijMie_m)*(sitecorr*sitecorr-(PXij*PXij+PYij*PYij+PZij*PZij)/RijSquared) &
                           + (Mie_n1*Mie_nRijMie_n - Mie_m1*Mie_mRijMie_m)*sitecorr*sitecorr)     
          FXi = FXi + FXij
          FYi = FYi + FYij
          FZi = FZi + FZij
          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij

        end do loop1
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local)      
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, merge( i1 - 1, i1, SameComponent )
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)

#if MPI_VER > 0
        if( SameComponent ) then
          j0 = i + 1
          j1 = i + N2
          if( EvenN .and. i > N2 ) j1 = j1 - 1
        else
          j0 = 1
          j1 = N1
        end if
loop3:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )
#else
        j0 = merge( i + 1, 1, SameComponent )
loop3:  do j = j0, j1
#endif
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = PXij - anint( RXij )
          PYij = PYij - anint( RYij )
          PZij = PZij - anint( RZij )
          RXij = RXij - anint( RXij )
          RYij = RYij - anint( RYij )
          RZij = RZij - anint( RZij )
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          if( RijSquared >= RCutoffSquared ) cycle loop3
          RijSquaredInv = SigmaSquared / RijSquared
          RijMie_nInv = RijSquaredInv**Mie_nHalf
          RijMie_mInv = RijSquaredInv**Mie_mHalf
          Mie_nRijMie_n = Mie_n * RijMie_nInv
          Mie_mRijMie_m = Mie_m * RijMie_mInv
          EPotLocal = EPotLocal + (RijMie_nInv - RijMie_mInv)
          Fij = EpsilonMie_aF * (Mie_nRijMie_n - Mie_mRijMie_m) * RijSquaredInv
          FXij = Fij * RXij
          FYij = Fij * RYij
          FZij = Fij * RZij
          VirialLocal = VirialLocal + (PXij * FXij + PYij * FYij + PZij * FZij)
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)/RijSquared
          d2EpotdV2Local = d2EpotdV2Local + EpsilonMie_a * Ninth * ((Mie_nRijMie_n - Mie_mRijMie_m)*(sitecorr*sitecorr-(PXij*PXij+PYij*PYij+PZij*PZij)/RijSquared) &
                           + (Mie_n1*Mie_nRijMie_n - Mie_m1*Mie_mRijMie_m)*sitecorr*sitecorr)     
          FXi = FXi + FXij
          FYi = FYi + FYij
          FZi = FZi + FZij
          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij
        end do loop3
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
      end do
!$OMP END DO 
    end if

!$OMP END PARALLEL

    ! Update potential energy and virial
   FX2 = FX2 + forceTempX
   FY2 = FY2 + forceTempY
   FZ2 = FZ2 + forceTempZ
   EPot = EPot + this%EpsilonMie_a * EPotLocal
   Virial = Virial + Third * VirialLocal * BoxLength
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:) * BoxLength
#endif
   d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotMIEMIE_Force


!==============================================================!
!  Subroutine TPotMIEMIE_Force_Trans                             !
!==============================================================!

  subroutine TPotMIEMIE_Force_Trans( this, EPot, Virial, d2EpotdV2, BoxLength )


    implicit none

    ! Declare arguments
    type(TPotMIEnmMIEnm)     :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength


    ! Declare local variables
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK)          :: SigmaSquared
    real(RK)          :: EpsilonMie_a, EpsilonMie_aF
    real(RK)          :: Mie_n, Mie_m, Mie_n1, Mie_m1, Mie_nHalf, Mie_mHalf, Mie_nRijMie_n, Mie_mRijMie_m
    real(RK)          :: RCutoffSquared
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: FXij, FYij, FZij, Fij
    real(RK)          :: RijSquared, RijSquaredInv, RijMie_nInv, RijMie_mInv
    real(RK)          :: EPotLocal, EPotLocal1, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
#if MPI_VER > 0
    integer           :: i0, N1, N2, ji
    logical           :: EvenN
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)

#if TRANS == 1
    real(RK), pointer, contiguous :: VSx1(:), VSy1(:), VSz1(:)
    real(RK), pointer, contiguous :: VSx2(:), VSy2(:), VSz2(:) 
    real(RK), pointer, contiguous :: VSux1(:),VSuy1(:),VSuz1(:)
    real(RK), pointer, contiguous :: VSux2(:),VSuy2(:),VSuz2(:)
    real(RK), pointer, contiguous :: VBx1(:), VBy1(:), VBz1(:)
    real(RK), pointer, contiguous :: VBx2(:), VBy2(:), VBz2(:)
    real(RK), pointer, contiguous :: Cx1(:) , Cy1(:) , Cz1(:)
    real(RK), pointer, contiguous :: Cx2(:) , Cy2(:) , Cz2(:)
    real(RK), pointer, contiguous :: tux1(:) , tuy1(:) , tuz1(:)
    real(RK), pointer, contiguous :: tux2(:) , tuy2(:) , tuz2(:)
    real(RK), pointer, contiguous :: tlx1(:) , tly1(:) , tlz1(:)
    real(RK), pointer, contiguous :: tlx2(:) , tly2(:) , tlz2(:)
    real(RK), pointer, contiguous :: tdx1(:) , tdy1(:) , tdz1(:)
    real(RK), pointer, contiguous :: tdx2(:) , tdy2(:) , tdz2(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
    real(RK)          :: VSTempX(1:this%Site2%NPart), VSTempY(1:this%Site2%NPart), VSTempZ(1:this%Site2%NPart)
    real(RK)          :: VSuTempX(1:this%Site2%NPart), VSuTempY(1:this%Site2%NPart), VSuTempZ(1:this%Site2%NPart)
    real(RK)          :: VBTempX(1:this%Site2%NPart), VBTempY(1:this%Site2%NPart), VBTempZ(1:this%Site2%NPart)
    real(RK)          :: CTempX(1:this%Site2%NPart), CTempY(1:this%Site2%NPart), CTempZ(1:this%Site2%NPart)
    real(RK)          :: tuTempX(1:this%Site2%NPart), tuTempY(1:this%Site2%NPart), tuTempZ(1:this%Site2%NPart)
    real(RK)          :: tlTempX(1:this%Site2%NPart), tlTempY(1:this%Site2%NPart), tlTempZ(1:this%Site2%NPart)
    real(RK)          :: tdTempX(1:this%Site2%NPart), tdTempY(1:this%Site2%NPart), tdTempZ(1:this%Site2%NPart)
    real(RK)          :: SigmaInvEpsMie_a
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSxij, VSyij, VSzij
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VSuxij,VSuyij,VSuzij
    real(RK)          :: VBxi, VByi, VBzi
    real(RK)          :: VBxij, VByij, VBzij
    real(RK)          :: Cxi,  Cyi,  Czi
    real(RK)          :: tuxi,  tuyi,  tuzi
    real(RK)          :: tuxij,  tuyij,  tuzij
    real(RK)          :: tlxi,  tlyi,  tlzi
    real(RK)          :: tlxij,  tlyij,  tlzij
    real(RK)          :: tdxi,  tdyi,  tdzi
    real(RK)          :: tdxij,  tdyij,  tdzij
    real(RK)          :: txii,  tyii , tzii
    real(RK)          :: txi ,  tyi  , tzi
    real(RK)          :: UU, BoxLength2
    real(RK)          :: r1x, r1y, r1z
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33

    VSTempX(:)  = 0._RK
    VSTempY(:)  = 0._RK
    VSTempZ(:)  = 0._RK
    VSuTempX(:) = 0._RK
    VSuTempY(:) = 0._RK
    VSuTempZ(:) = 0._RK
    VBTempX(:)  = 0._RK
    VBTempY(:)  = 0._RK
    VBTempZ(:)  = 0._RK
    CTempX(:)   = 0._RK
    CTempY(:)   = 0._RK
    CTempZ(:)   = 0._RK
    tuTempX(:)  = 0._RK
    tuTempY(:)  = 0._RK
    tuTempz(:)  = 0._RK
    tlTempX(:)  = 0._RK
    tlTempY(:)  = 0._RK
    tlTempz(:)  = 0._RK
    tdTempX(:)  = 0._RK
    tdTempY(:)  = 0._RK
    tdTempz(:)  = 0._RK
#endif



    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK


!$OMP PARALLEL PRIVATE(i, j, k, i1, j0, j1) &
!$OMP PRIVATE( RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE( PX1, PY1, PZ1, PX2, PY2, PZ2, FX1, FY1, FZ1, FX2, FY2) &
!$OMP PRIVATE(FZ2, SigmaSquared, EpsilonMie_a, EpsilonMie_aF, RCutoffSquared,EPotLocal1) &
!$OMP PRIVATE(RXi, RYi, RZi,  PXi, PYi, PZi,  FXi, FYi, FZi,  RXij, RYij, RZij, PXij, PYij, PZij) &
!$OMP PRIVATE(FXij, FYij, FZij, Fij, RijSquared, RijSquaredInv, RijMie_nInv, RijMie_mInv ) &
#if MPI_VER > 0
!$OMP PRIVATE(i0, N1, N2, ji, EvenN) &
#endif
#if  TRANS == 1
!$OMP PRIVATE(VSx1, VSy1, VSz1, VSux1,VSuy1,VSuz1, VBx1, VBy1, VBz1, Cx1, Cy1, Cz1) &
!$OMP PRIVATE( tux1 , tuy1 , tuz1, tlx1 , tly1 , tlz1, tdx1 , tdy1, tdz1) &
!$OMP PRIVATE( q1, q2, q3, q4, SigmaInvEpsMie_a, VSxi, VSyi, VSzi, VSuxi,VSuyi,VSuzi) &
!$OMP PRIVATE( VBxi, VByi, VBzi, Cxi,  Cyi,  Czi, tuxi,  tuyi,  tuzi, tlxi,  tlyi,  tlzi) &
!$OMP PRIVATE(  tdxi,  tdyi,  tdzi, txii,  tyii , tzii, txi ,  tyi  , tzi ) &
!$OMP PRIVATE(  UU , BoxLength2, r1x, r1y, r1z) &
!$OMP PRIVATE( A11, A12, A13, A21, A22, A23, A31, A32, A33) &
#endif
!$OMP PRIVATE( SameComponent )

      ! Assign local variables
    SameComponent = this%SameComponent
#if MPI_VER > 0
    N1 = this%Site2%NPart
    N2 = N1 / 2
    EvenN = mod( N1, 2 ) == 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
    j1 = this%Site2%NPart
#endif

    SigmaSquared = this%SigmaSquared
    EpsilonMie_a = this%EpsilonMie_a
    EpsilonMie_aF = this%EpsilonMie_aF
    Mie_n = this%Mie_n
    Mie_m = this%Mie_m
    Mie_n1 = Mie_n+1._RK
    Mie_m1 = Mie_m+1._RK
    Mie_nHalf = this%Mie_nHalf
    Mie_mHalf = this%Mie_mHalf
    RCutoffSquared = this%RCutoffSquaredScaled

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ


#if TRANS == 1
    SigmaInvEpsMie_a = EpsilonMie_a/Sqrt(this%SigmaSquared)
    BoxLength2   = BoxLength**2
    VSx1 => this%Site1%vsMIEx
    VSy1 => this%Site1%vsMIEy
    VSz1 => this%Site1%vsMIEz
    VSx2 => this%Site2%vsMIEx
    VSy2 => this%Site2%vsMIEy
    VSz2 => this%Site2%vsMIEz
    VBx1 => this%Site1%vbMIEx
    VBy1 => this%Site1%vbMIEy
    VBz1 => this%Site1%vbMIEz
    VBx2 => this%Site2%vbMIEx
    VBy2 => this%Site2%vbMIEy
    VBz2 => this%Site2%vbMIEz
    VSux1=> this%Site1%vsuMIEx
    VSuy1=> this%Site1%vsuMIEy
    VSuz1=> this%Site1%vsuMIEz
    VSux2=> this%Site2%vsuMIEx
    VSuy2=> this%Site2%vsuMIEy
    VSuz2=> this%Site2%vsuMIEz
    Cx1  => this%Site1%cMIEx
    Cy1  => this%Site1%cMIEy
    Cz1  => this%Site1%cMIEz
    Cx2  => this%Site2%cMIEx
    Cy2  => this%Site2%cMIEy
    Cz2  => this%Site2%cMIEz
    tux1 => this%Site1%tuMIEx
    tuy1 => this%Site1%tuMIEy
    tuz1 => this%Site1%tuMIEz
    tux2 => this%Site2%tuMIEx
    tuy2 => this%Site2%tuMIEy
    tuz2 => this%Site2%tuMIEz
    tlx1 => this%Site1%tlMIEx
    tly1 => this%Site1%tlMIEy
    tlz1 => this%Site1%tlMIEz
    tlx2 => this%Site2%tlMIEx
    tly2 => this%Site2%tlMIEy
    tlz2 => this%Site2%tlMIEz
    tdx1 => this%Site1%tdMIEx
    tdy1 => this%Site1%tdMIEy
    tdz1 => this%Site1%tdMIEz
    tdx2 => this%Site2%tdMIEx
    tdy2 => this%Site2%tdMIEy
    tdz2 => this%Site2%tdMIEz
    q1  => this%Site1%Q0r(:, 1)
    q2  => this%Site1%Q0r(:, 2)
    q3  => this%Site1%Q0r(:, 3)
    q4  => this%Site1%Q0r(:, 4)

#endif


    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i) 
#if  TRANS == 1
        !TRANSPORT_start
        VSxi = 0._RK
        VSyi = 0._RK
        VSzi = 0._RK
        VBxi = 0._RK
        VByi = 0._RK
        VBzi = 0._RK
        VSuxi= 0._RK
        VSuyi= 0._RK
        VSuzi= 0._RK
        Cxi  = Cx1(i)
        Cyi  = Cy1(i)
        Czi  = Cz1(i)
        tuxi = tux1(i)
        tuyi = tuy1(i)
        tuzi = tuz1(i)
        tlxi = tlx1(i)
        tlyi = tly1(i)
        tlzi = tlz1(i)
        tdxi = tdx1(i)
        tdyi = tdy1(i)
        tdzi = tdz1(i)
        r1x  = ( RXi-PXi ) * BoxLength2
        r1y  = ( RYi-PYi ) * BoxLength2
        r1z  = ( RZi-PZi ) * BoxLength2

        A11 = q1(i)**2 + q2(i)**2 - q3(i)**2 - q4(i)**2
        A12 = 2._RK * (q2(i) * q3(i) + q1(i) * q4(i))
        A13 = 2._RK * (q2(i) * q4(i) - q1(i) * q3(i))
        A21 = 2._RK * (q2(i) * q3(i) - q1(i) * q4(i))
        A22 = q1(i)**2 - q2(i)**2 + q3(i)**2 - q4(i)**2
        A23 = 2._RK * (q3(i) * q4(i) + q1(i) * q2(i))
        A31 = 2._RK * (q2(i) * q4(i) + q1(i) * q3(i))
        A32 = 2._RK * (q3(i) * q4(i) - q1(i) * q2(i))
        A33 = q1(i)**2 - q2(i)**2 - q3(i)**2 + q4(i)**2
#endif

!CDIR NODEP
#if OSMOP == 2
loop0:  do m=1,NBinsDen
          if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PXi < real(m)/NBinsDen-0.5_RK) then
              Bin1=m
              exit loop0
            end if
          end if
        end do loop0
#endif
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = RXij - anint( PXij )
          RYij = RYij - anint( PYij )
          RZij = RZij - anint( PZij )
          PXij = PXij - anint( PXij )
          PYij = PYij - anint( PYij )
          PZij = PZij - anint( PZij )
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          RijSquaredInv = SigmaSquared / RijSquared
          RijMie_nInv = RijSquaredInv**Mie_nHalf
          RijMie_mInv = RijSquaredInv**Mie_mHalf
          Mie_nRijMie_n = Mie_n * RijMie_nInv
          Mie_mRijMie_m = Mie_m * RijMie_mInv
          EPotLocal1 = RijMie_nInv - RijMie_mInv
          EPotLocal = EPotLocal + EPotLocal1
          Fij = EpsilonMie_aF * (Mie_nRijMie_n - Mie_mRijMie_m) * RijSquaredInv
          FXij = Fij * RXij
          FYij = Fij * RYij
          FZij = Fij * RZij
          VirialLocal = VirialLocal + (PXij * FXij + PYij * FYij + PZij * FZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
                Bin2=m 
                exit loop2
              end if
            end if
          end do loop2
          tempMin = min(Bin1, Bin2)
          tempMax = max(Bin1, Bin2)
          if(abs(PXij) .le. 0.5_RK) then
             VirialPart = (PXij * FXij + PYij * FYij + PZij * FZij)/(tempMax-tempMin+1._RK) 
             do m = tempMin, tempMax
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
             end do
          else
             VirialPart = (PXij * FXij + PYij * FYij + PZij * FZij)/(NBinsDen-tempMax+tempMin+1._RK) 
             do m = 1, tempMin
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
             end do
             do m = tempMax, NBinsDen
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
             end do
          end if
#endif
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)/RijSquared
          d2EpotdV2Local = d2EpotdV2Local + EpsilonMie_a * Ninth * ((Mie_nRijMie_n - Mie_mRijMie_m)*(sitecorr*sitecorr-(PXij*PXij+PYij*PYij+PZij*PZij)/RijSquared) &
                           + (Mie_n1*Mie_nRijMie_n - Mie_m1*Mie_mRijMie_m)*sitecorr*sitecorr)!xxxx MIE T      
          FXi = FXi + FXij
          FYi = FYi + FYij
          FZi = FZi + FZij
          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij

#if  TRANS == 1
          !TRANSPORT_start
          VSxij   = 0.5 * FXij * PYij
          VSyij   = 0.5 * FXij * PZij
          VSzij   = 0.5 * FYij * PZij
          VBxij   = 0.5 * FXij * PXij
          VByij   = 0.5 * FYij * PYij
          VBzij   = 0.5 * FZij * PZij
          VSuxij  = 0.5 * FYij * PXij
          VSuyij  = 0.5 * FZij * PXij
          VSuzij  = 0.5 * FZij * PYij

          VSxi   = VSxi + VSxij
          VSyi   = VSyi + VSyij
          VSzi   = VSzi + VSzij
          VBxi   = VBxi + VBxij
          VByi   = VByi + VByij
          VBzi   = VBzi + VBzij
          VSuxi  = VSuxi+ VSuxij
          VSuyi  = VSuyi+ VSuyij
          VSuzi  = VSuzi+ VSuzij
          VSTempX(j) = VSTempX(j) + VSxij
          VSTempY(j) = VSTempY(j) + VSyij
          VSTempZ(j) = VSTempZ(j) + VSzij
          VBTempX(j) = VBTempX(j) + VBxij
          VBTempY(j) = VBTempY(j) + VByij
          VBTempZ(j) = VBTempZ(j) + VBzij
          VSuTempX(j)= VSuTempX(j)+ VSuxij
          VSuTempY(j)= VSuTempY(j)+ VSuyij
          VSuTempZ(j)= VSuTempZ(j)+ VSuzij

          UU     = 0.5 * EPotLocal1*this%EpsilonMie_a
          Cxi    = Cxi  + UU
          Cyi    = Cyi  + UU
          Czi    = Czi  + UU
          CTempX(j) = CTempX(j) + UU
          CTempY(j) = CTempY(j) + UU
          CTempZ(j) = CTempZ(j) + UU
          txii   = r1y * FZij - r1z * FYij
          tyii   = r1z * FXij - r1x * FZij
          tzii   = r1x * FYij - r1y * FXij
          txi    = A11 * txii + A12 * tyii + A13 * tzii
          tyi    = A21 * txii + A22 * tyii + A23 * tzii
          tzi    = A31 * txii + A32 * tyii + A33 * tzii
          tuxij  = 0.5 * PXij * tyi
          tuyij  = 0.5 * PXij * tzi
          tuzij  = 0.5 * PYij * tzi
          tlxij  = 0.5 * PYij * txi
          tlyij  = 0.5 * PZij * txi
          tlzij  = 0.5 * PZij * tyi
          tdxij  = 0.5 * PXij * txi
          tdyij  = 0.5 * PYij * tyi
          tdzij  = 0.5 * PZij * tzi

          tuxi   = tuxi + tuxij
          tuyi   = tuyi + tuyij
          tuzi   = tuzi + tuzij
          tlxi   = tlxi + tlxij
          tlyi   = tlyi + tlyij
          tlzi   = tlzi + tlzij
          tdxi   = tdxi + tdxij
          tdyi   = tdyi + tdyij
          tdzi   = tdzi + tdzij
          tuTempX(j)= tuTempX(j) + tuxij
          tuTempY(j)= tuTempY(j) + tuyij
          tuTempZ(j)= tuTempZ(j) + tuzij
          tlTempX(j)= tlTempX(j) + tlxij
          tlTempY(j)= tlTempY(j) + tlyij
          tlTempZ(j)= tlTempZ(j) + tlzij
          tdTempX(j)= tdTempX(j) + tdxij
          tdTempY(j)= tdTempY(j) + tdyij
          tdTempZ(j)= tdTempZ(j) + tdzij
#endif


        end do loop1
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
#if TRANS == 1
        !TRANSPORT_start
        VSx1(i) = VSx1(i) + VSxi * BoxLength
        VSy1(i) = VSy1(i) + VSyi * BoxLength
        VSz1(i) = VSz1(i) + VSzi * BoxLength
        VBx1(i) = VBx1(i) + VBxi * BoxLength
        VBy1(i) = VBy1(i) + VByi * BoxLength
        VBz1(i) = VBz1(i) + VBzi * BoxLength
        VSux1(i)= VSux1(i)+ VSuxi* BoxLength
        VSuy1(i)= VSuy1(i)+ VSuyi* BoxLength
        VSuz1(i)= VSuz1(i)+ VSuzi* BoxLength
        Cx1(i)  = Cxi
        Cy1(i)  = Cyi
        Cz1(i)  = Czi
        tux1(i) = tuxi
        tuy1(i) = tuyi
        tuz1(i) = tuzi
        tlx1(i) = tlxi
        tly1(i) = tlyi
        tlz1(i) = tlzi
        tdx1(i) = tdxi
        tdy1(i) = tdyi
        tdz1(i) = tdzi
        !TRANSPORT_END
#endif

      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local)      
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, merge( i1 - 1, i1, SameComponent )
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)

#if MPI_VER > 0
        if( SameComponent ) then
          j0 = i + 1
          j1 = i + N2
          if( EvenN .and. i > N2 ) j1 = j1 - 1
        else
          j0 = 1
          j1 = N1
        end if
loop3:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )

#else
          j0 = merge( i + 1, 1, SameComponent )
loop3:  do j = j0, j1
#endif
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = PXij - anint( RXij )
          PYij = PYij - anint( RYij )
          PZij = PZij - anint( RZij )
          RXij = RXij - anint( RXij )
          RYij = RYij - anint( RYij )
          RZij = RZij - anint( RZij )
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          if( RijSquared >= RCutoffSquared ) cycle loop3
          RijSquaredInv = SigmaSquared / RijSquared
          RijMie_nInv = RijSquaredInv**Mie_nHalf
          RijMie_mInv = RijSquaredInv**Mie_mHalf
          Mie_nRijMie_n = Mie_n * RijMie_nInv
          Mie_mRijMie_m = Mie_m * RijMie_mInv
          EPotLocal = EPotLocal + (RijMie_nInv - RijMie_mInv)
          Fij = EpsilonMie_aF * (Mie_nRijMie_n - Mie_mRijMie_m) * RijSquaredInv
          FXij = Fij * RXij
          FYij = Fij * RYij
          FZij = Fij * RZij
          VirialLocal = VirialLocal + (PXij * FXij + PYij * FYij + PZij * FZij)
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)/RijSquared
          d2EpotdV2Local = d2EpotdV2Local + EpsilonMie_a * Ninth * ((Mie_nRijMie_n - Mie_mRijMie_m)*(sitecorr*sitecorr-(PXij*PXij+PYij*PYij+PZij*PZij)/RijSquared) &
                           + (Mie_n1*Mie_nRijMie_n - Mie_m1*Mie_mRijMie_m)*sitecorr*sitecorr)     !xxxx MIE SS T
          FXi = FXi + FXij
          FYi = FYi + FYij
          FZi = FZi + FZij
          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij

        end do loop3
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
      end do
!$OMP END DO 
    end if

!$OMP END PARALLEL
 
    ! Update potential energy and virial
   FX2 = FX2 + forceTempX
   FY2 = FY2 + forceTempY
   FZ2 = FZ2 + forceTempZ

#if  TRANS == 1
   VSx2 = VSx2 + VSTempX*BoxLength
   VSy2 = VSy2 + VSTempY*BoxLength
   VSz2 = VSz2 + VSTempZ*BoxLength

   VSux2 = VSux2 + VSuTempX*BoxLength
   VSuy2 = VSuy2 + VSuTempY*BoxLength
   VSuz2 = VSuz2 + VSuTempZ*BoxLength

   VBx2 = VBx2 + VBTempX*BoxLength
   VBy2 = VBy2 + VBTempY*BoxLength
   VBz2 = VBz2 + VBTempZ*BoxLength

   Cx2 = Cx2 + CTempX
   Cy2 = Cy2 + CTempY
   Cz2 = Cz2 + CTempZ

   tux2 = tux2 + tuTempX
   tuy2 = tuy2 + tuTempY
   tuz2 = tuz2 + tuTempZ

   tdx2 = tdx2 + tdTempX
   tdy2 = tdy2 + tdTempY
   tdz2 = tdz2 + tdTempZ

   tlx2 = tlx2 + tlTempX
   tly2 = tly2 + tlTempY
   tlz2 = tlz2 + tlTempZ
#endif

   EPot = EPot + this%EpsilonMie_a * EPotLocal
   Virial = Virial + Third * VirialLocal * BoxLength

#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:) * BoxLength
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotMIEMIE_Force_Trans


!==============================================================!
!  Subroutine TPotMIEMIE_RDF                                   !
!==============================================================!

  subroutine TPotMIEMIE_RDF( this, RDFdr )

    implicit none

    ! Declare arguments
    type(TPotMIEnmMIEnm)     :: this
    real(RK), intent(in)     :: RDFdr

    !RDF RDFdr und RDFSchalenIndex
    real(RK)          :: distance
    integer           :: RDFSchalenIndex

    ! Declare local variables
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: RXi, RYi, RZi
    integer           :: i, j, k

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ

    ! Loop over molecules
    do i = 1, this%Site1%NPart
      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)

!CDIR NODEP
loop1:do k = 1, this%NInCutoff(i)
        j = this%CutoffPartner(k, i)
        RXij = RXi - RX2(j)
        RYij = RYi - RY2(j)
        RZij = RZi - RZ2(j)
        RXij = RXij - anint( RXij )
        RYij = RYij - anint( RYij )
        RZij = RZij - anint( RZij )

!RDF in Schalen sortieren
        distance = sqrt(RXij*RXij + RYij*RYij + RZij*RZij)
        RDFSchalenIndex = INT(distance/RDFdr) + 1
        if (RDFSchalenIndex .le. RDFNumberShells) then
          this%RDFSum(RDFSchalenIndex) = this%RDFSum(RDFSchalenIndex) + 1
        endif
      end do loop1
    end do

  end subroutine TPotMIEMIE_RDF
  

!==============================================================!
!  Subroutine TPotMIEMIE_ChemicalPotential                     !
!==============================================================!

  subroutine TPotMIEMIE_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotMIEnmMIEnm) :: this
    real(RK), pointer, contiguous    :: EPotTest(:)
    real(RK), intent(in) :: BoxLength

    ! Declare local variables
    real(RK)          :: SigmaSquared
    real(RK)          :: EpsilonMie_a
    real(RK)          :: Mie_nHalf, Mie_mHalf
    real(RK)          :: RCutoffSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: RijSquared, RijSquaredInv, RijMie_nInv, RijMie_mInv
    real(RK)          :: EPotLocal
    integer           :: N2
    integer           :: i, j, k

    ! Assign local variables
    N2 = this%Site2%NPart
    SigmaSquared = this%SigmaSquared
    EpsilonMie_a = this%EpsilonMie_a
    Mie_nHalf = this%Mie_nHalf
    Mie_mHalf = this%Mie_mHalf
    RCutoffSquared = this%RCutoffSquaredScaled

    ! Assign pointers
    RX1 => this%Site1%RXTest
    RY1 => this%Site1%RYTest
    RZ1 => this%Site1%RZTest
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    PX1 => this%Site1%PXTest
    PY1 => this%Site1%PYTest
    PZ1 => this%Site1%PZTest
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (RXi,RYi,RZi,PXi,PYi,PZi) &
!$OMP PRIVATE (RXij,RYij,RZij,PXij,PYij,PZij) &
!$OMP PRIVATE (RijSquared,RijSquaredInv,RijMie_nInv, RijMie_mInv) &
!$OMP PRIVATE (EpotLocal,i,j,k) 

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
!$OMP DO
      do i = 1, this%Site1%NTest
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        EPotLocal = 0._RK
!CDIR NODEP
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = RXij - anint( PXij )
          RYij = RYij - anint( PYij )
          RZij = RZij - anint( PZij )
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          RijSquaredInv = SigmaSquared / RijSquared
          RijMie_nInv = RijSquaredInv**Mie_nHalf
          RijMie_mInv = RijSquaredInv**Mie_mHalf
          EPotLocal = EPotLocal + (RijMie_nInv - RijMie_mInv)
        end do loop1
        EPotTest(i) = EPotTest(i) + EpsilonMie_a * EPotLocal
      end do
!$OMP END DO
    else

      ! Loop over test particles
!$OMP DO
      do i = 1, this%Site1%NTest
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        EPotLocal = 0._RK
!CDIR NODEP
loop2:  do j = 1, N2
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          RXij = RXij - anint( RXij )
          RYij = RYij - anint( RYij )
          RZij = RZij - anint( RZij )
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          if( RijSquared >= RCutoffSquared ) cycle loop2
          RijSquaredInv = SigmaSquared / RijSquared
          RijMie_nInv = RijSquaredInv**Mie_nHalf
          RijMie_mInv = RijSquaredInv**Mie_mHalf
          EPotLocal = EPotLocal + (RijMie_nInv - RijMie_mInv)
        end do loop2
        EPotTest(i) = EPotTest(i) + EpsilonMie_a * EPotLocal
      end do
!$OMP END DO
    end if
!$OMP END PARALLEL
  end subroutine TPotMIEMIE_ChemicalPotential



!==============================================================!
!  Subroutine TPotMIEMIE_UpdateBoxLength                       !
!==============================================================!

  subroutine TPotMIEMIE_UpdateBoxLength( this, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotMIEnmMIEnm) :: this
    real(RK), intent(in) :: BoxLength

    ! Declare local variables
    real(RK) :: BoxLengthInv

    ! Update constants
    BoxLengthInv = 1._RK / BoxLength
    this%BoxLengthInv = BoxLengthInv
    this%BoxLengthThird = Third * BoxLength
    this%SigmaSquared = (this%Sigma * BoxLengthInv)**2
    this%EpsilonMie_aF = this%EpsilonMie_a * BoxLengthInv / this%SigmaSquared
    this%RCutoffSquaredScaled = this%RCutoffSquared * BoxLengthInv**2

  end subroutine TPotMIEMIE_UpdateBoxLength



!==============================================================!
!  Subroutine TPotCC_Construct                                 !
!==============================================================!

  subroutine TPotCC_Construct( this, i1, i2, j1, j2, &
&                              Molecule1, Molecule2, RCutoff, RFEpsilon )

    implicit none

    ! Declare arguments
    type(TPotChargeCharge)      :: this
    integer, intent(in)         :: i1, i2, j1, j2
    type(TMolecule), intent(in) :: Molecule1, Molecule2
    real(RK), intent(in)        :: RCutoff
    real(RK), intent(in)        :: RFEpsilon

    ! Construct potential
    this%Site1 => Molecule1%SiteCharge(j1)
    this%Site2 => Molecule2%SiteCharge(j2)
    this%SameComponent = i1 == i2
    this%Epsilon = this%Site1%e * this%Site2%e
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2
    this%RFConstant = this%Epsilon / RCutoff**3 * (RFEpsilon - 1._RK) / (2._RK * RFEpsilon + 1._RK)

  end subroutine TPotCC_Construct



!==============================================================!
!  Subroutine TPotCC_Destruct                                  !
!==============================================================!

  subroutine TPotCC_Destruct( this )

    implicit none

    ! Declare arguments
    type(TPotChargeCharge) :: this

    ! Destroy potential
    continue

  end subroutine TPotCC_Destruct



!==============================================================!
!  Subroutine TPotCC_Force                                     !
!==============================================================!

  subroutine TPotCC_Force( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)

    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ      ! Site-Site-Einheitvektor
    real(RK)          :: RijInv
    real(RK)          :: EPotLocal, EPotLocal1, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: Rij2
    integer           :: i, j, k, i1
#if MPI_VER > 0
    integer           :: i0
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ

    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK
    Epsilon = this%Epsilon

!$OMP PARALLEL &
#if MPI_VER > 0
!$OMP FIRSTPRIVATE (i0) &
#endif
!$OMP FIRSTPRIVATE(i1) &
!$OMP PRIVATE (Plen2,sitecorr) &
!$OMP PRIVATE (RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE (RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE (eX, eY, eZ  , RijInv, EPotLocal1,  i, j, k)

    ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local)    
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)
      FXi = FX1(i)
      FYi = FY1(i)
      FZi = FZ1(i)
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)

!CDIR NODEP
#if OSMOP == 2
loop0:do m=1,NBinsDen
        if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
          if (PXi < real(m)/NBinsDen-0.5_RK) then
            Bin1=m
            exit loop0
          end if
        end if
      end do loop0
#endif
loop1:do k = 1, this%NInCutoff(i)
        j = this%CutoffPartner(k, i)
        RXij = RXi - RX2(j)
        RYij = RYi - RY2(j)
        RZij = RZi - RZ2(j)
        PXij = PXi - PX2(j)
        PYij = PYi - PY2(j)
        PZij = PZi - PZ2(j)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        Rij2   = RXij*RXij + RYij*RYij + RZij*RZij

        RijInv = 1._RK / sqrt( Rij2 )

        eX = RXij * RijInv
        eY = RYij * RijInv
        eZ = RZij * RijInv
        EPotLocal1 = Epsilon * RijInv
        EPotLocal  = EPotLocal + EPotLocal1
        VirialLocal = VirialLocal + (EPotLocal1 * RijInv * (eX * PXij + eY * PYij + eZ * PZij))
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
              Bin2=m 
              exit loop2
            end if
          end if
        end do loop2
        tempMin = min(Bin1, Bin2)
        tempMax = max(Bin1, Bin2)
        if(abs(PXij) .le. 0.5_RK) then
            VirialPart = (EPotLocal1 * RijInv * (eX * PXij + eY * PYij + eZ * PZij))/(tempMax-tempMin+1._RK) 
            do m = tempMin, tempMax
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        else
            VirialPart = (EPotLocal1 * RijInv * (eX * PXij + eY * PYij + eZ * PZij))/(NBinsDen-tempMax+tempMin+1._RK) 
            do m = 1, tempMin
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
            do m = tempMax, NBinsDen
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        end if
#endif
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (RXij * PXij + RYij * PYij + RZij * PZij)*RijInv*RijInv
        d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 * (3._RK*sitecorr*sitecorr - Plen2*RijInv*RijInv)*Ninth !XXXX CC
        FXij = EPotLocal1 * RijInv * eX
        FYij = EPotLocal1 * RijInv * eY
        FZij = EPotLocal1 * RijInv * eZ
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        forceTempX(j) = forceTempX(j) - FXij
        forceTempY(j) = forceTempY(j) - FYij
        forceTempZ(j) = forceTempZ(j) - FZij

      end do loop1
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi

    end do
!$OMP END DO
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotCC_Force


!==============================================================!
!  Subroutine TPotCC_Force_Ewald                               !
!==============================================================!

  subroutine TPotCC_Force_Ewald( this, EPot, Virial, d2EpotdV2, BoxLength, Kappa )

    implicit none

    ! Declare arguments
    type(TPotChargeCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength
    real(RK), intent(in)     :: Kappa

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ         ! Site-Site-Einheitvektor
    real(RK)          :: RijInv, Rij
    real(RK)          :: EPotLocal, EPotLocal1, VirialLocal
    real(RK)          :: approx, Faktor
    real(RK)          :: Fij,KappaRij
    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: Rij2
    integer           :: i, j, k, i1, i2
#if MPI_VER > 0
    integer           :: i0
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ


    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
    i2 = 0
#else
    i1 = this%Site1%NPart
    i2 = this%Site2%NPart
#endif
    Epsilon = this%Epsilon
    Faktor = 2._RK/sqrt(Pi) * Kappa
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    EPotLocal=0._RK
    VirialLocal=0._RK

!$OMP PARALLEL DEFAULT(SHARED) &
#if MPI_VER > 0
!$OMP FIRSTPRIVATE ( i0) &
#endif
!$OMP FIRSTPRIVATE (i1, i2) &
!$OMP PRIVATE ( approx, Fij,KappaRij,Rij2) &
!$OMP PRIVATE ( RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi) &
!$OMP PRIVATE ( RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE ( eX, eY, eZ  , RijInv,Rij, EPotLocal1,  i, j, k)


    ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal)    
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)
      FXi = FX1(i)
      FYi = FY1(i)
      FZi = FZ1(i)
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)

!CDIR NODEP
#if OSMOP == 2
loop0:do m=1,NBinsDen
        if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
          if (PXi < real(m)/NBinsDen-0.5_RK) then
            Bin1=m
            exit loop0
          end if
        end if
      end do loop0
#endif
 loop1:do k = 1, this%NInCutoff(i),1
         j = this%CutoffPartner(k, i)
        RXij = RXi - RX2(j)
        RYij = RYi - RY2(j)
        RZij = RZi - RZ2(j)
        PXij = PXi - PX2(j)
        PYij = PYi - PY2(j)
        PZij = PZi - PZ2(j)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength

        Rij2   = RXij*RXij + RYij*RYij + RZij*RZij
        Rij =  sqrt(Rij2)


        RijInv = 1._RK /  Rij 

        KappaRij = Kappa*Rij
        call ErrorApprox(this, KappaRij,approx)

        eX = RXij * RijInv
        eY = RYij * RijInv
        eZ = RZij * RijInv
        EPotLocal1 = Epsilon * RijInv * approx
        EPotLocal  = EPotLocal + EPotLocal1
        Fij  = (EPotLocal1 + Faktor*exp(-KappaRij**2)*Epsilon) * RijInv
        VirialLocal = VirialLocal + ( Fij * (eX * PXij + eY * PYij + eZ * PZij))
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
              Bin2=m 
              exit loop2
            end if
          end if
        end do loop2
        tempMin = min(Bin1, Bin2)
        tempMax = max(Bin1, Bin2)
        if(abs(PXij) .le. 0.5_RK) then
            VirialPart = ( Fij * (eX * PXij + eY * PYij + eZ * PZij))/(tempMax-tempMin+1._RK) 
            do m = tempMin, tempMax
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        else
            VirialPart = ( Fij * (eX * PXij + eY * PYij + eZ * PZij))/(NBinsDen-tempMax+tempMin+1._RK) 
            do m = 1, tempMin
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
            do m = tempMax, NBinsDen
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        end if
#endif
        FXij = Fij * eX
        FYij = Fij * eY
        FZij = Fij * eZ
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        forceTempX(j) = forceTempX(j) - FXij
        forceTempY(j) = forceTempY(j) - FYij
        forceTempZ(j) = forceTempZ(j) - FZij


      end do loop1
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
    end do
!$OMP END DO
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif

  end subroutine TPotCC_Force_Ewald


!==============================================================!
!  Subroutine TPotCC_Force_Trans                               !
!==============================================================!

  subroutine TPotCC_Force_Trans( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ      ! Site-Site-Einheitvektor
    real(RK)          :: RijInv
    real(RK)          :: EPotLocal, EPotLocal1, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    real(RK)          :: Rij2
    integer           :: i, j, k, i1
#if MPI_VER > 0
    integer           :: i0
#endif    
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)

#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer, contiguous :: VSx1(:), VSy1(:), VSz1(:)
    real(RK), pointer, contiguous :: VSx2(:), VSy2(:), VSz2(:)
    real(RK), pointer, contiguous :: VSux1(:), VSuy1(:), VSuz1(:)
    real(RK), pointer, contiguous :: VSux2(:), VSuy2(:), VSuz2(:)
    real(RK), pointer, contiguous :: VBx1(:), VBy1(:), VBz1(:)
    real(RK), pointer, contiguous :: VBx2(:), VBy2(:), VBz2(:)
    real(RK), pointer, contiguous :: Cx1(:) , Cy1(:) , Cz1(:)
    real(RK), pointer, contiguous :: Cx2(:) , Cy2(:) , Cz2(:)
    real(RK), pointer, contiguous :: tux1(:) , tuy1(:) , tuz1(:)
    real(RK), pointer, contiguous :: tux2(:) , tuy2(:) , tuz2(:)
    real(RK), pointer, contiguous :: tlx1(:) , tly1(:) , tlz1(:)
    real(RK), pointer, contiguous :: tlx2(:) , tly2(:) , tlz2(:)
    real(RK), pointer, contiguous :: tdx1(:) , tdy1(:) , tdz1(:)
    real(RK), pointer, contiguous :: tdx2(:) , tdy2(:) , tdz2(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
    real(RK)          :: VSTempX(1:this%Site2%NPart), VSTempY(1:this%Site2%NPart), VSTempZ(1:this%Site2%NPart)
    real(RK)          :: VSuTempX(1:this%Site2%NPart), VSuTempY(1:this%Site2%NPart), VSuTempZ(1:this%Site2%NPart)
    real(RK)          :: VBTempX(1:this%Site2%NPart), VBTempY(1:this%Site2%NPart), VBTempZ(1:this%Site2%NPart)
    real(RK)          :: CTempX(1:this%Site2%NPart), CTempY(1:this%Site2%NPart), CTempZ(1:this%Site2%NPart)
    real(RK)          :: tuTempX(1:this%Site2%NPart), tuTempY(1:this%Site2%NPart), tuTempZ(1:this%Site2%NPart)
    real(RK)          :: tlTempX(1:this%Site2%NPart), tlTempY(1:this%Site2%NPart), tlTempZ(1:this%Site2%NPart)
    real(RK)          :: tdTempX(1:this%Site2%NPart), tdTempY(1:this%Site2%NPart), tdTempZ(1:this%Site2%NPart)
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSxij, VSyij, VSzij
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VSuxij,VSuyij,VSuzij
    real(RK)          :: VBxi, VByi, VBzi
    real(RK)          :: VBxij, VByij, VBzij
    real(RK)          :: Cxi,  Cyi,  Czi
    real(RK)          :: tuxi,  tuyi,  tuzi
    real(RK)          :: tuxij,  tuyij,  tuzij
    real(RK)          :: tlxi,  tlyi,  tlzi
    real(RK)          :: tlxij,  tlyij,  tlzij
    real(RK)          :: tdxi,  tdyi,  tdzi
    real(RK)          :: tdxij,  tdyij,  tdzij
    real(RK)          :: txii,  tyii , tzii
    real(RK)          :: txi ,  tyi  , tzi
    real(RK)          :: r1x, r1y, r1z, UU
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33

    VSTempX(:) = 0._RK
    VSTempY(:) = 0._RK
    VSTempZ(:) = 0._RK
    VSuTempX(:)= 0._RK
    VSuTempY(:)= 0._RK
    VSuTempZ(:)= 0._RK
    VBTempX(:) = 0._RK
    VBTempY(:) = 0._RK
    VBTempZ(:) = 0._RK
    CTempX(:)  = 0._RK
    CTempY(:)  = 0._RK
    CTempZ(:)  = 0._RK
    tuTempX(:) = 0._RK
    tuTempY(:) = 0._RK
    tuTempz(:) = 0._RK
    tlTempX(:) = 0._RK
    tlTempY(:) = 0._RK
    tlTempz(:) = 0._RK
    tdTempX(:) = 0._RK
    tdTempY(:) = 0._RK
    tdTempz(:) = 0._RK

    !TRANSPORT_END
#endif

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

!$OMP PARALLEL &
!$OMP PRIVATE( Epsilon, RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE(  FX1, FY1, FZ1, FX2, FY2, FZ2) &
!$OMP PRIVATE(Plen2,sitecorr, PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE(   RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE(   RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
#if  TRANS == 1
!$OMP PRIVATE(VSx, VSy, VSz ,VSux,VSuy,VSuz, VBx, VBy, VBz, Cx , Cy , Cz) &
!$OMP PRIVATE( tux , tuy , tuz, tlx , tly , tlz, tdx , tdy , tdz) &
!$OMP PRIVATE( q1, q2, q3, q4, VSxi, VSyi, VSzi, VSuxi,VSuyi,VSuzi) &
!$OMP PRIVATE( VBxi, VByi, VBzi, Cxi,  Cyi,  Czi, tuxi,  tuyi,  tuzi, tlxi,  tlyi,  tlzi) &
!$OMP PRIVATE(  tdxi,  tdyi,  tdzi, txii,  tyii , tzii, txi ,  tyi  , tzi ) &
!$OMP PRIVATE(  UU ,  Uxi,  Uyi, Uzi, r1x, r1y, r1z) &
!$OMP PRIVATE( A11, A12, A13, A21, A22, A23, A31, A32, A33) &
#endif
#if MPI_VER > 0
!$OMP PRIVATE ( i0) &
#endif
!$OMP PRIVATE ( eX, eY, eZ  , RijInv, EPotLocal1,  i, j, k, i1)


    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    Epsilon = this%Epsilon

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
#if  TRANS == 1
    !TRANSPORT_start
    VSx1 => this%Site1%vsCx
    VSy1 => this%Site1%vsCy
    VSz1 => this%Site1%vsCz
    VSx2 => this%Site2%vsCx
    VSy2 => this%Site2%vsCy
    VSz2 => this%Site2%vsCz
    VBx1 => this%Site1%vbCx
    VBy1 => this%Site1%vbCy
    VBz1 => this%Site1%vbCz
    VBx2 => this%Site2%vbCx
    VBy2 => this%Site2%vbCy
    VBz2 => this%Site2%vbCz
    VSux1=> this%Site1%vsuCx
    VSuy1=> this%Site1%vsuCy
    VSuz1=> this%Site1%vsuCz
    VSux2=> this%Site2%vsuCx
    VSuy2=> this%Site2%vsuCy
    VSuz2=> this%Site2%vsuCz
    Cx1  => this%Site1%cCx
    Cy1  => this%Site1%cCy
    Cz1  => this%Site1%cCz
    Cx2  => this%Site2%cCx
    Cy2  => this%Site2%cCy
    Cz2  => this%Site2%cCz
    tux1 => this%Site1%tuCx
    tuy1 => this%Site1%tuCy
    tuz1 => this%Site1%tuCz
    tux2 => this%Site2%tuCx
    tuy2 => this%Site2%tuCy
    tuz2 => this%Site2%tuCz
    tlx1 => this%Site1%tlCx
    tly1 => this%Site1%tlCy
    tlz1 => this%Site1%tlCz
    tlx2 => this%Site2%tlCx
    tly2 => this%Site2%tlCy
    tlz2 => this%Site2%tlCz
    tdx1 => this%Site1%tdCx
    tdy1 => this%Site1%tdCy
    tdz1 => this%Site1%tdCz
    tdx2 => this%Site1%tdCx
    tdy2 => this%Site1%tdCy
    tdz2 => this%Site1%tdCz
    q1  => this%Site1%Q0r(:, 1)
    q2  => this%Site1%Q0r(:, 2)
    q3  => this%Site1%Q0r(:, 3)
    q4  => this%Site1%Q0r(:, 4)
!TRANSPORT_END
#endif

    ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local)    
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)
      FXi = FX1(i)
      FYi = FY1(i)
      FZi = FZ1(i)
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)
#if  TRANS == 1
      !TRANSPORT_start
      VSxi = VSx1(i)
      VSyi = VSy1(i)
      VSzi = VSz1(i)
      VBxi = VBx1(i)
      VByi = VBy1(i)
      VBzi = VBz1(i)
      VSuxi= VSux1(i)
      VSuyi= VSuy1(i)
      VSuzi= VSuz1(i)
      Cxi  = Cx1(i)
      Cyi  = Cy1(i)
      Czi  = Cz1(i)
      tuxi = tux1(i)
      tuyi = tuy1(i)
      tuzi = tuz1(i)
      tlxi = tlx1(i)
      tlyi = tly1(i)
      tlzi = tlz1(i)
      tdxi = tdx1(i)
      tdyi = tdy1(i)
      tdzi = tdz1(i)
      r1x  = ( RXi-PXi ) * BoxLength
      r1y  = ( RYi-PYi ) * BoxLength
      r1z  = ( RZi-PZi ) * BoxLength
      A11  = q1(i)**2 + q2(i)**2 - q3(i)**2 - q4(i)**2
      A12  = 2._RK * (q2(i) * q3(i) + q1(i) * q4(i))
      A13  = 2._RK * (q2(i) * q4(i) - q1(i) * q3(i))
      A21  = 2._RK * (q2(i) * q3(i) - q1(i) * q4(i))
      A22  = q1(i)**2 - q2(i)**2 + q3(i)**2 - q4(i)**2
      A23  = 2._RK * (q3(i) * q4(i) + q1(i) * q2(i))
      A31  = 2._RK * (q2(i) * q4(i) + q1(i) * q3(i))
      A32  = 2._RK * (q3(i) * q4(i) - q1(i) * q2(i))
      A33  = q1(i)**2 - q2(i)**2 - q3(i)**2 + q4(i)**2
      !TRANSPORT_END
#endif

!CDIR NODEP
#if OSMOP == 2
loop0:do m=1,NBinsDen
        if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
          if (PXi < real(m)/NBinsDen-0.5_RK) then
            Bin1=m
            exit loop0
          end if
        end if
      end do loop0
#endif
loop1:do k = 1, this%NInCutoff(i)
        j = this%CutoffPartner(k, i)
        RXij = RXi - RX2(j)
        RYij = RYi - RY2(j)
        RZij = RZi - RZ2(j)
        PXij = PXi - PX2(j)
        PYij = PYi - PY2(j)
        PZij = PZi - PZ2(j)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        Rij2   = RXij*RXij + RYij*RYij + RZij*RZij

        RijInv = 1._RK / sqrt( Rij2 )

        eX = RXij * RijInv
        eY = RYij * RijInv
        eZ = RZij * RijInv
        EPotLocal1 = Epsilon * RijInv
        EPotLocal  = EPotLocal + EPotLocal1
        VirialLocal = VirialLocal + (EPotLocal1 * RijInv * (eX * PXij + eY * PYij + eZ * PZij))
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
              Bin2=m 
              exit loop2
            end if
          end if
        end do loop2
        tempMin = min(Bin1, Bin2)
        tempMax = max(Bin1, Bin2)
        if(abs(PXij) .le. 0.5_RK) then
            VirialPart = (EPotLocal1 * RijInv * (eX * PXij + eY * PYij + eZ * PZij))/(tempMax-tempMin+1._RK) 
            do m = tempMin, tempMax
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        else
            VirialPart = (EPotLocal1 * RijInv * (eX * PXij + eY * PYij + eZ * PZij))/(NBinsDen-tempMax+tempMin+1._RK) 
            do m = 1, tempMin
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
            do m = tempMax, NBinsDen
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        end if
#endif
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (RXij * PXij + RYij * PYij + RZij * PZij)*RijInv*RijInv
        d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 * (3._RK*sitecorr*sitecorr - Plen2*RijInv*RijInv)*Ninth !xxxx CC T
        FXij = EPotLocal1 * RijInv * eX
        FYij = EPotLocal1 * RijInv * eY
        FZij = EPotLocal1 * RijInv * eZ
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        forceTempX(j) = forceTempX(j) - FXij
        forceTempY(j) = forceTempY(j) - FYij
        forceTempZ(j) = forceTempZ(j) - FZij

#if TRANS==1
        !TRANSPORT_start

        VSxij   = 0.5 * FXij * PYij
        VSyij   = 0.5 * FXij * PZij
        VSzij   = 0.5 * FYij * PZij
        VBxij   = 0.5 * FXij * PXij
        VByij   = 0.5 * FYij * PYij
        VBzij   = 0.5 * FZij * PZij
        VSuxij  = 0.5 * FYij * PXij
        VSuyij  = 0.5 * FZij * PXij
        VSuzij  = 0.5 * FZij * PYij

        VSxi   = VSxi + VSxij
        VSyi   = VSyi + VSyij
        VSzi   = VSzi + VSzij
        VBxi   = VBxi + VBxij
        VByi   = VByi + VByij
        VBzi   = VBzi + VBzij
        VSuxi  = VSuxi+ VSuxij
        VSuyi  = VSuyi+ VSuyij
        VSuzi  = VSuzi+ VSuzij        

        VSTempX(j) = VSTempX(j) + VSxij
        VSTempY(j) = VSTempY(j) + VSyij
        VSTempZ(j) = VSTempZ(j) + VSzij
        VBTempX(j) = VBTempX(j) + VBxij
        VBTempY(j) = VBTempY(j) + VByij
        VBTempZ(j) = VBTempZ(j) + VBzij
        VSuTempX(j)= VSuTempX(j)+ VSuxij
        VSuTempY(j)= VSuTempY(j)+ VSuyij
        VSuTempZ(j)= VSuTempZ(j)+ VSuzij
      
        UU     = 0.5 * EpotLocal1
        Cxi    = Cxi  + UU 
        Cyi    = Cyi  + UU
        Czi    = Czi  + UU
        CTempX(j) = CTempX(j) + UU
        CTempY(j) = CTempY(j) + UU
        CTempZ(j) = CTempZ(j) + UU

        txii   = r1y * FZij - r1z * FYij
        tyii   = r1z * FXij - r1x * FZij
        tzii   = r1x * FYij - r1y * FXij
        txi    = A11 * txii + A12 * tyii + A13 * tzii
        tyi    = A21 * txii + A22 * tyii + A23 * tzii
        tzi    = A31 * txii + A32 * tyii + A33 * tzii
        tuxij   = 0.5 * PXij* tyi
        tuyij   = 0.5 * PXij* tzi
        tuzij   = 0.5 * PYij* tzi
        tlxij   = 0.5 * PYij* txi
        tlyij   = 0.5 * PZij* txi
        tlzij   = 0.5 * PZij* tyi
        tdxij   = 0.5 * PXij* txi
        tdyij   = 0.5 * PYij* tyi
        tdzij   = 0.5 * PZij* tzi
        tuxi   = tuxi + tuxij
        tuyi   = tuyi + tuyij
        tuzi   = tuzi + tuzij
        tlxi   = tlxi + tlxij
        tlyi   = tlyi + tlyij
        tlzi   = tlzi + tlzij
        tdxi   = tdxi + tdxij
        tdyi   = tdyi + tdyij
        tdzi   = tdzi + tdzij
        tuTempX(j)= tuTempX(j) + tuxij
        tuTempY(j)= tuTempY(j) + tuyij
        tuTempZ(j)= tuTempZ(j) + tuzij
        tlTempX(j)= tlTempX(j) + tlxij
        tlTempY(j)= tlTempY(j) + tlyij
        tlTempZ(j)= tlTempZ(j) + tlzij
        tdTempX(j)= tdTempX(j) + tdxij
        tdTempY(j)= tdTempY(j) + tdyij
        tdTempZ(j)= tdTempZ(j) + tdzij
#endif

      end do loop1
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
#if  TRANS == 1
      !TRANSPORT_start
      VSx1(i) = VSxi
      VSy1(i) = VSyi
      VSz1(i) = VSzi
      VBx1(i) = VBxi
      VBy1(i) = VByi
      VBz1(i) = VBzi
      VSux1(i)= VSuxi
      VSuy1(i)= VSuyi
      VSuz1(i)= VSuzi
      Cx1(i)  = Cxi
      Cy1(i)  = Cyi
      Cz1(i)  = Czi
      tux1(i) = tuxi
      tuy1(i) = tuyi
      tuz1(i) = tuzi
      tlx1(i) = tlxi
      tly1(i) = tlyi
      tlz1(i) = tlzi
      tdx1(i) = tdxi
      tdy1(i) = tdyi
      tdz1(i) = tdzi
      !TRANSPORT_END

#endif
    end do
!$OMP END DO
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ

#if  TRANS == 1
   VSx2 = VSx2 + VSTempX
   VSy2 = VSy2 + VSTempY
   VSz2 = VSz2 + VSTempZ

   VSux2 = VSux2 + VSuTempX
   VSuy2 = VSuy2 + VSuTempY
   VSuz2 = VSuz2 + VSuTempZ

   VBx2 = VBx2 + VBTempX
   VBy2 = VBy2 + VBTempY
   VBz2 = VBz2 + VBTempZ

   Cx2 = Cx2 + CTempX
   Cy2 = Cy2 + CTempY
   Cz2 = Cz2 + CTempZ

   tux2 = tux2 + tuTempX
   tuy2 = tuy2 + tuTempY
   tuz2 = tuz2 + tuTempZ

   tdx2 = tdx2 + tdTempX
   tdy2 = tdy2 + tdTempY
   tdz2 = tdz2 + tdTempZ

   tlx2 = tlx2 + tlTempX
   tly2 = tly2 + tlTempY
   tlz2 = tlz2 + tlTempZ

#endif
    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotCC_Force_Trans



!==============================================================!
!  Subroutine TPotCC_Force_Ewald_Trans                         !
!==============================================================!

  subroutine TPotCC_Force_Ewald_Trans( this, EPot, Virial, d2EpotdV2, BoxLength, Kappa )

    implicit none

    ! Declare arguments
    type(TPotChargeCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength
    real(RK), intent(in)     :: Kappa

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ         ! Site-Site-Einheitvektor
    real(RK)          :: RijInv, Rij
    real(RK)          :: EPotLocal, EPotLocal1, VirialLocal
    real(RK)          :: approx, Faktor
    real(RK)          :: Fij,KappaRij
    real(RK)          :: Rij2
    integer           :: i, j, k, i1, i2
#if MPI_VER > 0
    integer           :: i0
#endif    
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)

#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer, contiguous :: VSx1(:), VSy1(:), VSz1(:)
    real(RK), pointer, contiguous :: VSx2(:), VSy2(:), VSz2(:)
    real(RK), pointer, contiguous :: VBx1(:), VBy1(:), VBz1(:)
    real(RK), pointer, contiguous :: VBx2(:), VBy2(:), VBz2(:)
    real(RK)          :: VSTempX(1:this%Site2%NPart), VSTempY(1:this%Site2%NPart), VSTempZ(1:this%Site2%NPart)
    real(RK)          :: VBTempX(1:this%Site2%NPart), VBTempY(1:this%Site2%NPart), VBTempZ(1:this%Site2%NPart)
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VBxi, VByi, VBzi


    VSTempX(:)=0._RK
    VSTempY(:)=0._RK
    VSTempZ(:)=0._RK
    VBTempX(:)=0._RK
    VBTempY(:)=0._RK
    VBTempZ(:)=0._RK

    !TRANSPORT_END
#endif

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    forceTempX(:)=0._RK   
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    EPotLocal=0._RK
    VirialLocal=0._RK
  
!$OMP PARALLEL &
!$OMP PRIVATE( Epsilon, RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE( approx, Faktor, Fij,KappaRij ) &
!$OMP PRIVATE(  FX1, FY1, FZ1, FX2, FY2, FZ2) &
!$OMP PRIVATE( PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE(   RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE(   RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
#if  TRANS == 1
!$OMP PRIVATE(VSx1, VSy1, VSz1 ,VBx1, VBy1, VBz1) &
!$OMP PRIVATE( VSxi, VSyi, VSzi) &
!$OMP PRIVATE( VBxi, VByi, VBzi) &
#endif
#if MPI_VER > 0
!$OMP PRIVATE ( i0) &
#endif
!$OMP PRIVATE ( eX, eY, eZ  , RijInv, Rij, EPotLocal1,  i, j, k, i1, i2)


    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
    i2 = this%Site2%NPart
#endif
    Epsilon = this%Epsilon
    Faktor = 2._RK/sqrt(Pi) * Kappa

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
#if  TRANS == 1
    !TRANSPORT_start

    VSx1 => this%Site1%vsCx
    VSy1 => this%Site1%vsCy
    VSz1 => this%Site1%vsCz
    VSx2 => this%Site2%vsCx
    VSy2 => this%Site2%vsCy
    VSz2 => this%Site2%vsCz
    VBx1 => this%Site1%vbCx
    VBy1 => this%Site1%vbCy
    VBz1 => this%Site1%vbCz
    VBx2 => this%Site2%vbCx
    VBy2 => this%Site2%vbCy
    VBz2 => this%Site2%vbCz

!TRANSPORT_END
#endif

    ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal)    
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)
      FXi = FX1(i)
      FYi = FY1(i)
      FZi = FZ1(i)
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)
#if  TRANS == 1
      !TRANSPORT_start
      VSxi= VSx1(i)
      VSyi= VSy1(i)
      VSzi= VSz1(i)
      VBxi= VBx1(i)
      VByi= VBy1(i)
      VBzi= VBz1(i)
      !TRANSPORT_END
#endif

!CDIR NODEP
#if OSMOP == 2
loop0:do m=1,NBinsDen
        if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
          if (PXi < real(m)/NBinsDen-0.5_RK) then
            Bin1=m
            exit loop0
          end if
        end if
      end do loop0
#endif
 loop1:do k = 1, this%NInCutoff(i),1
         j = this%CutoffPartner(k, i)
        RXij = RXi - RX2(j)
        RYij = RYi - RY2(j)
        RZij = RZi - RZ2(j)
        PXij = PXi - PX2(j)
        PYij = PYi - PY2(j)
        PZij = PZi - PZ2(j)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength

        Rij2   = RXij*RXij + RYij*RYij + RZij*RZij
        Rij =  sqrt(Rij2)


        RijInv = 1._RK /  Rij 

        KappaRij = Kappa*Rij
        call ErrorApprox(this, KappaRij,approx)

        eX = RXij * RijInv
        eY = RYij * RijInv
        eZ = RZij * RijInv
        EPotLocal1 = Epsilon * RijInv * approx
        EPotLocal  = EPotLocal + EPotLocal1
        Fij  = (EPotLocal1 + Faktor*exp(-KappaRij**2)*Epsilon) * RijInv
        VirialLocal = VirialLocal + (Fij * (eX * PXij + eY * PYij + eZ * PZij))
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
              Bin2=m 
              exit loop2
            end if
          end if
        end do loop2
        tempMin = min(Bin1, Bin2)
        tempMax = max(Bin1, Bin2)
        if(abs(PXij) .le. 0.5_RK) then
            VirialPart = (Fij * (eX * PXij + eY * PYij + eZ * PZij))/(tempMax-tempMin+1._RK) 
            do m = tempMin, tempMax
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        else
            VirialPart = (Fij * (eX * PXij + eY * PYij + eZ * PZij))/(NBinsDen-tempMax+tempMin+1._RK) 
            do m = 1, tempMin
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
            do m = tempMax, NBinsDen
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        end if
#endif
        FXij = Fij * eX
        FYij = Fij * eY
        FZij = Fij * eZ
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij

        forceTempX(j) = forceTempX(j) - FXij
        forceTempY(j) = forceTempY(j) - FYij
        forceTempZ(j) = forceTempZ(j) - FZij

#if TRANS==1
        !TRANSPORT_start vielleicht
        VSxi   = VSxi + 0.5 * FXij * PYij
        VSyi   = VSyi + 0.5 * FXij * PZij
        VSzi   = VSzi + 0.5 * FYij * PZij
        VBxi   = VBxi + 0.5 * FXij * PXij
        VByi   = VByi + 0.5 * FYij * PYij
        VBzi   = VBzi + 0.5 * FZij * PZij
        VSTempX(j) = VSTempX(j) + 0.5 * FXij * PYij
        VSTempY(j) = VSTempY(j) + 0.5 * FXij * PZij
        VSTempZ(j) = VSTempZ(j) + 0.5 * FYij * PZij
        VBTempX(j) = VBTempX(j) + 0.5 * FXij * PXij
        VBTempY(j) = VBTempY(j) + 0.5 * FYij * PYij
        VBTempZ(j) = VBTempZ(j) + 0.5 * FZij * PZij
#endif

      end do loop1
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
#if  TRANS == 1
      !TRANSPORT_start
      VSx1(i) = VSxi
      VSy1(i) = VSyi
      VSz1(i) = VSzi
      VBx1(i) = VBxi
      VBy1(i) = VByi
      VBz1(i) = VBzi
#endif

    end do
!$OMP END DO
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ

#if  TRANS == 1
    VSx2 = VSx2 + VSTempX*BoxLength
    VSy2 = VSy2 + VSTempY*BoxLength
    VSz2 = VSz2 + VSTempZ*BoxLength
    VBx2 = VBx2 + VBTempX*BoxLength
    VBy2 = VBy2 + VBTempY*BoxLength
    VBz2 = VBz2 + VBTempZ*BoxLength
#endif


    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif

  end subroutine TPotCC_Force_Ewald_Trans


!==============================================================!
!  Subroutine TPotCC_ChemicalPotential                         !
!==============================================================!

  subroutine TPotCC_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeCharge) :: this
    real(RK), pointer, contiguous      :: EPotTest(:)
    real(RK), intent(in)   :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: RijInv, RijSquared
    real(RK)          :: EPotLocal
    integer           :: i, j, k, i1


    ! Assign local variables
    i1 = this%Site1%NTest
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

    ! Assign pointers
    RX1 => this%Site1%RXTest
    RY1 => this%Site1%RYTest
    RZ1 => this%Site1%RZTest
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    PX1 => this%Site1%PXTest
    PY1 => this%Site1%PYTest
    PZ1 => this%Site1%PZTest
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

   ! Loop over test particles
   do i = 1, i1
     RXi = RX1(i)
     RYi = RY1(i)
     RZi = RZ1(i)
     PXi = PX1(i)
     PYi = PY1(i)
     PZi = PZ1(i)
     EPotLocal = 0._RK



!CDIR NODEP
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij

          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop1
          end if
           RijInv = 1._RK / sqrt( RijSquared )

           EPotLocal = EPotLocal + Epsilon * RijInv
        end do loop1

        EPotTest(i) = EPotTest(i) + EPotLocal

   end do

  end subroutine TPotCC_ChemicalPotential



!==============================================================!
!  Subroutine TPotCD_Construct                                 !
!==============================================================!

  subroutine TPotCD_Construct( this, i1, i2, j1, j2, Molecule1, Molecule2, RCutoff )

    implicit none

    ! Declare arguments
    type(TPotChargeDipole)      :: this
    integer, intent(in)         :: i1, i2, j1, j2
    type(TMolecule), intent(in) :: Molecule1, Molecule2
    real(RK), intent(in)        :: RCutoff

    ! Construct potential
    this%Site1 => Molecule1%SiteCharge(j1)
    this%Site2 => Molecule2%SiteDipole(j2)
    this%SameComponent = i1 == i2
    this%Epsilon = this%Site1%e * this%Site2%D
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

  end subroutine TPotCD_Construct


!==============================================================!
!  Subroutine TPotCD_Destruct                                  !
!==============================================================!

  subroutine TPotCD_Destruct( this )

    implicit none

    ! Declare arguments
    type(TPotChargeDipole) :: this

    ! Destroy potential
    continue

  end subroutine TPotCD_Destruct



!==============================================================!
!  Subroutine TPotCD_Force                                     !
!==============================================================!

  subroutine TPotCD_Force( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeDipole)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1, Epsilon2
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: TX2(:), TY2(:), TZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv
    real(RK)          :: CosTheta, CosTheta3
    real(RK)          :: EPotLocal, Viriallocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1
    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: momTempX(1:this%Site2%NPart)
    real(RK)          :: momTempY(1:this%Site2%NPart)
    real(RK)          :: momTempZ(1:this%Site2%NPart)
#if MPI_VER > 0
    integer           :: i0
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ    

    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    momTempX(:)=0._RK
    momTempY(:)=0._RK
    momTempZ(:)=0._RK    
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

!$OMP PARALLEL &
!$OMP PRIVATE ( Epsilon, Epsilon1, Epsilon2, RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (  FX1, FY1, FZ1, OX2, OY2, OZ2) &
!$OMP PRIVATE ( Plen2, PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (   RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE (   RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE ( OXj, OYj, OZj, eX, eY, eZ, RijSquaredInv, RijInv) &
#if MPI_VER > 0
!$OMP PRIVATE ( CosTheta, CosTheta3,  i, j, k, i1) &
!$OMP PRIVATE ( i0)
#else
!$OMP PRIVATE ( CosTheta, CosTheta3, i, j, k, i1)
#endif

    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    Epsilon = this%Epsilon

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ

    ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)       
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)
      FXi = FX1(i)
      FYi = FY1(i)
      FZi = FZ1(i)
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)

!CDIR NODEP
#if OSMOP == 2
loop0:do m=1,NBinsDen
        if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
          if (PXi < real(m)/NBinsDen-0.5_RK) then
            Bin1=m
            exit loop0
          end if
        end if
      end do loop0
#endif
loop1:do k = 1, this%NInCutoff(i)
        j = this%CutoffPartner(k, i)
        RXij = RXi - RX2(j)
        RYij = RYi - RY2(j)
        RZij = RZi - RZ2(j)
        PXij = PXi - PX2(j)
        PYij = PYi - PY2(j)
        PZij = PZi - PZ2(j)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        OXj = OX2(j)
        OYj = OY2(j)
        OZj = OZ2(j)
        RijSquaredInv = 1._RK / ( RXij*RXij + RYij*RYij + RZij*RZij )
        RijInv = sqrt( RijSquaredInv )
        eX = RXij * RijInv                                                      ! Einheitsabstandvektor nach Price
        eY = RYij * RijInv
        eZ = RZij * RijInv
        CosTheta  = OXj * ex + OYj * eY + OZj * eZ                              ! cos(alpha) nach Price
        CosTheta3 = 3._RK * CosTheta
        Epsilon1 = Epsilon * RijSquaredInv
        Epsilon2 = Epsilon1 * RijInv
        EPotLocal  = EPotLocal + Epsilon1 * CosTheta                            ! Uebereinstimmumg mit Price
        FXij = Epsilon2 * ( CosTheta3 * eX - OXj )                              ! F2 bei Price
        FYij = Epsilon2 * ( CosTheta3 * eY - OYj )
        FZij = Epsilon2 * ( CosTheta3 * eZ - OZj )
        VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F2*R_COM_Price; stimmt so
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
              Bin2=m 
              exit loop2
            end if
          end if
        end do loop2
        tempMin = min(Bin1, Bin2)
        tempMax = max(Bin1, Bin2)
        if(abs(PXij) .le. 0.5_RK) then
            VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
            do m = tempMin, tempMax
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        else
            VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
            do m = 1, tempMin
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
            do m = tempMax, NBinsDen
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        end if
#endif
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + Epsilon1*CosTheta*(8._RK*sitecorr*sitecorr-2._RK*Plen2*RijSquaredInv)*Ninth   !xxxx2 CD
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij

        forceTempX(j) = forceTempX(j) - FXij
        forceTempY(j) = forceTempY(j) - FYij
        forceTempZ(j) = forceTempZ(j) - FZij
        momTempX(j) = momTempX(j) - Epsilon1 * eX   
        momTempY(j) = momTempY(j) - Epsilon1 * eY   
        momTempZ(j) = momTempZ(j) - Epsilon1 * eZ   

      end do loop1
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
    end do
!$OMP END DO
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ
    TX2 = TX2 + momTempX                                 
    TY2 = TY2 + momTempY
    TZ2 = TZ2 + momTempZ

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotCD_Force

!==============================================================!
!  Subroutine TPotCD_Force_Trans                               !
!==============================================================!

  subroutine TPotCD_Force_Trans( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeDipole)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1, Epsilon2
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: TX2(:), TY2(:), TZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv
    real(RK)          :: CosTheta, CosTheta3
    real(RK)          :: EPotLocal, EPotLocal1, Viriallocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1
    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: momTempX(1:this%Site2%NPart)
    real(RK)          :: momTempY(1:this%Site2%NPart)
    real(RK)          :: momTempZ(1:this%Site2%NPart)

#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer, contiguous :: VSx1(:), VSy1(:), VSz1(:)
    real(RK), pointer, contiguous :: VSx2(:), VSy2(:), VSz2(:)
    real(RK), pointer, contiguous :: VSux1(:), VSuy1(:), VSuz1(:)
    real(RK), pointer, contiguous :: VSux2(:), VSuy2(:), VSuz2(:)
    real(RK), pointer, contiguous :: VBx1(:), VBy1(:), VBz1(:)
    real(RK), pointer, contiguous :: VBx2(:), VBy2(:), VBz2(:)
    real(RK), pointer, contiguous :: Cx1(:) , Cy1(:) , Cz1(:)
    real(RK), pointer, contiguous :: Cx2(:) , Cy2(:) , Cz2(:)
    real(RK), pointer, contiguous :: tux1(:) , tuy1(:) , tuz1(:)
    real(RK), pointer, contiguous :: tux2(:) , tuy2(:) , tuz2(:)
    real(RK), pointer, contiguous :: tlx1(:) , tly1(:) , tlz1(:)
    real(RK), pointer, contiguous :: tlx2(:) , tly2(:) , tlz2(:)
    real(RK), pointer, contiguous :: tdx1(:) , tdy1(:) , tdz1(:)
    real(RK), pointer, contiguous :: tdx2(:) , tdy2(:) , tdz2(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
    real(RK)          :: VSTempX(1:this%Site2%NPart), VSTempY(1:this%Site2%NPart), VSTempZ(1:this%Site2%NPart)
    real(RK)          :: VSuTempX(1:this%Site2%NPart), VSuTempY(1:this%Site2%NPart), VSuTempZ(1:this%Site2%NPart)
    real(RK)          :: VBTempX(1:this%Site2%NPart), VBTempY(1:this%Site2%NPart), VBTempZ(1:this%Site2%NPart)
    real(RK)          :: CTempX(1:this%Site2%NPart), CTempY(1:this%Site2%NPart), CTempZ(1:this%Site2%NPart)
    real(RK)          :: tuTempX(1:this%Site2%NPart), tuTempY(1:this%Site2%NPart), tuTempZ(1:this%Site2%NPart)
    real(RK)          :: tlTempX(1:this%Site2%NPart), tlTempY(1:this%Site2%NPart), tlTempZ(1:this%Site2%NPart)
    real(RK)          :: tdTempX(1:this%Site2%NPart), tdTempY(1:this%Site2%NPart), tdTempZ(1:this%Site2%NPart)
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VBxi, VByi, VBzi
    real(RK)          :: Cxi,  Cyi,  Czi
    real(RK)          :: tuxi,  tuyi,  tuzi
    real(RK)          :: tlxi,  tlyi,  tlzi
    real(RK)          :: tdxi,  tdyi,  tdzi
    real(RK)          :: txii,  tyii , tzii
    real(RK)          :: txi ,  tyi  , tzi
    real(RK)          :: UU
    real(RK)          :: r1x, r1y, r1z
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33
    !TRANSPORT_END
#endif
   
#if MPI_VER > 0
    integer           :: i0
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ    
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    momTempX(:)=0._RK
    momTempY(:)=0._RK
    momTempZ(:)=0._RK    
    EPotLocal =0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

#if  TRANS == 1
    VSTempX(:) = 0._RK
    VSTempY(:) = 0._RK
    VSTempZ(:) = 0._RK
    VSuTempX(:)= 0._RK
    VSuTempY(:)= 0._RK
    VSuTempZ(:)= 0._RK
    VBTempX(:) = 0._RK
    VBTempY(:) = 0._RK
    VBTempZ(:) = 0._RK
    CTempX(:)  = 0._RK
    CTempY(:)  = 0._RK
    CTempZ(:)  = 0._RK
    tuTempX(:) = 0._RK
    tuTempY(:) = 0._RK
    tuTempz(:) = 0._RK
    tlTempX(:) = 0._RK
    tlTempY(:) = 0._RK
    tlTempz(:) = 0._RK
    tdTempX(:) = 0._RK
    tdTempY(:) = 0._RK
    tdTempz(:) = 0._RK
#endif


!$OMP PARALLEL &
!$OMP PRIVATE ( Epsilon, Epsilon1, Epsilon2, RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (  FX1, FY1, FZ1, OX2, OY2, OZ2) &
!$OMP PRIVATE ( Plen2,PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (   RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE (   RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE ( OXj, OYj, OZj, eX, eY, eZ, RijSquaredInv, RijInv) &
#if  TRANS == 1
!$OMP PRIVATE(VSx1, VSy1, VSz1,VSux1,VSuy1,VSuz1, VBx1, VBy1, VBz1, Cx1, Cy1, Cz1) &
!$OMP PRIVATE(VSx2, VSy2, VSz2,VSux2,VSuy2,VSuz2, VBx2, VBy2, VBz2, Cx2, Cy2, Cz2) &
!$OMP PRIVATE( tux1, tuy1, tuz1, tlx1, tly1, tlz1, tdx1, tdy1, tdz1) &
!$OMP PRIVATE( tux2, tuy2, tuz2, tlx2, tly2, tlz2, tdx2, tdy2, tdz2) &
!$OMP PRIVATE( q1, q2, q3, q4, VSxi, VSyi, VSzi, VSuxi,VSuyi,VSuzi) &
!$OMP PRIVATE( VBxi, VByi, VBzi, Cxi,  Cyi,  Czi, tuxi,  tuyi,  tuzi, tlxi,  tlyi,  tlzi) &
!$OMP PRIVATE(  tdxi,  tdyi,  tdzi, txii,  tyii , tzii, txi ,  tyi  , tzi ) &
!$OMP PRIVATE(  UU , r1x, r1y, r1z) &
!$OMP PRIVATE( A11, A12, A13, A21, A22, A23, A31, A32, A33) &
#endif

#if MPI_VER > 0
!$OMP PRIVATE ( CosTheta, CosTheta3,  i, j, k, i1) &
!$OMP PRIVATE ( i0)
#else
!$OMP PRIVATE ( CosTheta, CosTheta3, i, j, k, i1)
#endif

    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    Epsilon = this%Epsilon

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
#if  TRANS == 1
    !TRANSPORT_start
    VSx1 => this%Site1%vsCx
    VSy1 => this%Site1%vsCy
    VSz1 => this%Site1%vsCz
    VSx2 => this%Site2%vsDx
    VSy2 => this%Site2%vsDy
    VSz2 => this%Site2%vsDz
    VBx1 => this%Site1%vbCx
    VBy1 => this%Site1%vbCy
    VBz1 => this%Site1%vbCz
    VBx2 => this%Site2%vbDx
    VBy2 => this%Site2%vbDy
    VBz2 => this%Site2%vbdz
    VSux1=> this%Site1%vsuCx
    VSuy1=> this%Site1%vsuCy
    VSuz1=> this%Site1%vsuCz
    VSux2=> this%Site2%vsuDx
    VSuy2=> this%Site2%vsuDy
    VSuz2=> this%Site2%vsuDz
    Cx1  => this%Site1%cCx
    Cy1  => this%Site1%cCy
    Cz1  => this%Site1%cCz
    Cx2  => this%Site2%cDx
    Cy2  => this%Site2%cDy
    Cz2  => this%Site2%cDz
    tux1 => this%Site1%tuCx
    tuy1 => this%Site1%tuCy
    tuz1 => this%Site1%tuCz
    tux2 => this%Site2%tuDx
    tuy2 => this%Site2%tuDy
    tuz2 => this%Site2%tuDz
    tlx1 => this%Site1%tlCx
    tly1 => this%Site1%tlCy
    tlz1 => this%Site1%tlCz
    tlx2 => this%Site2%tlDx
    tly2 => this%Site2%tlDy
    tlz2 => this%Site2%tlDz
    tdx1 => this%Site1%tdCx
    tdy1 => this%Site1%tdCy
    tdz1 => this%Site1%tdCz
    tdx2 => this%Site2%tdDx
    tdy2 => this%Site2%tdDy
    tdz2 => this%Site2%tdDz
    q1  => this%Site1%Q0r(:, 1)
    q2  => this%Site1%Q0r(:, 2)
    q3  => this%Site1%Q0r(:, 3)
    q4  => this%Site1%Q0r(:, 4)
!TRANSPORT_END
#endif

    ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)       
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)
      FXi = FX1(i)
      FYi = FY1(i)
      FZi = FZ1(i)
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)

#if  TRANS == 1
      !TRANSPORT_start
      VSxi = VSx1(i)
      VSyi = VSy1(i)
      VSzi = VSz1(i)
      VBxi = VBx1(i)
      VByi = VBy1(i)
      VBzi = VBz1(i)
      VSuxi= VSux1(i)
      VSuyi= VSuy1(i)
      VSuzi= VSuz1(i)
      Cxi  = Cx1(i)
      Cyi  = Cy1(i)
      Czi  = Cz1(i)
      tuxi = tux1(i)
      tuyi = tuy1(i)
      tuzi = tuz1(i)
      tlxi = tlx1(i)
      tlyi = tly1(i)
      tlzi = tlz1(i)
      tdxi = tdx1(i)
      tdyi = tdy1(i)
      tdzi = tdz1(i)
      r1x  = ( RXi-PXi ) * BoxLength
      r1y  = ( RYi-PYi ) * BoxLength
      r1z  = ( RZi-PZi ) * BoxLength
      A11 = q1(i)**2 + q2(i)**2 - q3(i)**2 - q4(i)**2
      A12 = 2._RK * (q2(i) * q3(i) + q1(i) * q4(i))
      A13 = 2._RK * (q2(i) * q4(i) - q1(i) * q3(i))
      A21 = 2._RK * (q2(i) * q3(i) - q1(i) * q4(i))
      A22 = q1(i)**2 - q2(i)**2 + q3(i)**2 - q4(i)**2
      A23 = 2._RK * (q3(i) * q4(i) + q1(i) * q2(i))
      A31 = 2._RK * (q2(i) * q4(i) + q1(i) * q3(i))
      A32 = 2._RK * (q3(i) * q4(i) - q1(i) * q2(i))
      A33 = q1(i)**2 - q2(i)**2 - q3(i)**2 + q4(i)**2
      !TRANSPORT_END
#endif

!CDIR NODEP
#if OSMOP == 2
loop0:do m=1,NBinsDen
        if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
          if (PXi < real(m)/NBinsDen-0.5_RK) then
            Bin1=m
            exit loop0
          end if
        end if
      end do loop0
#endif
loop1:do k = 1, this%NInCutoff(i)
        j = this%CutoffPartner(k, i)
        RXij = RXi - RX2(j)
        RYij = RYi - RY2(j)
        RZij = RZi - RZ2(j)
        PXij = PXi - PX2(j)
        PYij = PYi - PY2(j)
        PZij = PZi - PZ2(j)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        OXj = OX2(j)
        OYj = OY2(j)
        OZj = OZ2(j)
        RijSquaredInv = 1._RK / ( RXij*RXij + RYij*RYij + RZij*RZij )
        RijInv = sqrt( RijSquaredInv )
        eX = RXij * RijInv                                                      ! Einheitsabstandvektor nach Price
        eY = RYij * RijInv
        eZ = RZij * RijInv
        CosTheta  = OXj * ex + OYj * eY + OZj * eZ                              ! cos(alpha) nach Price
        CosTheta3 = 3._RK * CosTheta
        Epsilon1 = Epsilon * RijSquaredInv
        Epsilon2 = Epsilon1 * RijInv
        EPotlocal1 = Epsilon1 * CosTheta             !Define EPotLocal 1 
        EPotLocal  = EPotLocal + EPotlocal1                         ! Uebereinstimmumg mit Price
        FXij = Epsilon2 * ( CosTheta3 * eX - OXj )                              ! F2 bei Price
        FYij = Epsilon2 * ( CosTheta3 * eY - OYj )
        FZij = Epsilon2 * ( CosTheta3 * eZ - OZj )
        VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F2*R_COM_Price; stimmt so
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
              Bin2=m 
              exit loop2
            end if
          end if
        end do loop2
        tempMin = min(Bin1, Bin2)
        tempMax = max(Bin1, Bin2)
        if(abs(PXij) .le. 0.5_RK) then
            VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
            do m = tempMin, tempMax
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        else
            VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
            do m = 1, tempMin
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
            do m = tempMax, NBinsDen
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        end if
#endif
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + EPotlocal1*(8._RK*sitecorr*sitecorr-2._RK*Plen2*RijSquaredInv)*Ninth   !xxxx2 CD T
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        forceTempX(j) = forceTempX(j) - FXij
        forceTempY(j) = forceTempY(j) - FYij
        forceTempZ(j) = forceTempZ(j) - FZij
        momTempX(j) = momTempX(j) - Epsilon1 * eX   
        momTempY(j) = momTempY(j) - Epsilon1 * eY   
        momTempZ(j) = momTempZ(j) - Epsilon1 * eZ   

#if TRANS==1
        !TRANSPORT_start 
        VSxi   = VSxi + 0.5 * FXij * PYij
        VSyi   = VSyi + 0.5 * FXij * PZij
        VSzi   = VSzi + 0.5 * FYij * PZij
        VBxi   = VBxi + 0.5 * FXij * PXij
        VByi   = VByi + 0.5 * FYij * PYij
        VBzi   = VBzi + 0.5 * FZij * PZij
        VSuxi  = VSuxi+ 0.5 * FYij * PXij
        VSuyi  = VSuyi+ 0.5 * FZij * PXij
        VSuzi  = VSuzi+ 0.5 * FZij * PYij
        VSTempX(j) = VSTempX(j) + 0.5*FXij * PYij
        VSTempY(j) = VSTempY(j) + 0.5*FXij * PZij
        VSTempZ(j) = VSTempZ(j) + 0.5*FYij * PZij
        VBTempX(j) = VBTempX(j) + 0.5*FXij * PXij
        VBTempY(j) = VBTempY(j) + 0.5*FYij * PYij
        VBTempZ(j) = VBTempZ(j) + 0.5*FZij * PZij
        VSuTempX(j)= VSuTempX(j)+ 0.5*FYij * PXij
        VSuTempY(j)= VSuTempY(j)+ 0.5*FZij * PXij
        VSuTempZ(j)= VSuTempZ(j)+ 0.5*FZij * PYij
       
        UU    =  0.5 * EpotLocal1        
        Cxi    = Cxi  + UU
        Cyi    = Cyi  + UU
        Czi    = Czi  + UU
        CTempX(j) = CTempX(j) + UU
        CTempY(j) = CTempY(j) + UU
        CTempZ(j) = CTempZ(j) + UU

        txii   = r1y * FZij - r1z * FYij
        tyii   = r1z * FXij - r1x * FZij
        tzii   = r1x * FYij - r1y * FXij
        txi    = A11 * txii + A12 * tyii + A13 * tzii
        tyi    = A21 * txii + A22 * tyii + A23 * tzii
        tzi    = A31 * txii + A32 * tyii + A33 * tzii
        tuxi   = tuxi + 0.5 * PXij* tyi
        tuyi   = tuyi + 0.5 * PXij* tzi
        tuzi   = tuzi + 0.5 * PYij* tzi
        tlxi   = tlxi + 0.5 * PYij* txi
        tlyi   = tlyi + 0.5 * PZij* txi
        tlzi   = tlzi + 0.5 * PZij* tyi
        tdxi   = tdxi + 0.5 * PXij* txi
        tdyi   = tdyi + 0.5 * PYij* tyi
        tdzi   = tdzi + 0.5 * PZij* tzi
        tuTempX(j)= tuTempX(j) + 0.5*PXij*tyi
        tuTempY(j)= tuTempY(j) + 0.5*PXij*tzi
        tuTempZ(j)= tuTempZ(j) + 0.5*PYij*tzi
        tlTempX(j)= tlTempX(j) + 0.5*PYij*txi
        tlTempY(j)= tlTempY(j) + 0.5*PZij*txi
        tlTempZ(j)= tlTempZ(j) + 0.5*PZij*tyi
        tdTempX(j)= tdTempX(j) + 0.5*PXij*txi
        tdTempY(j)= tdTempY(j) + 0.5*PYij*tyi
        tdTempZ(j)= tdTempZ(j) + 0.5*PZij*tzi  
#endif

      end do loop1

      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi

#if TRANS == 1
      !TRANSPORT_start
      VSx1(i) = VSxi
      VSy1(i) = VSyi
      VSz1(i) = VSzi
      VBx1(i) = VBxi
      VBy1(i) = VByi
      VBz1(i) = VBzi
      VSux1(i)= VSuxi
      VSuy1(i)= VSuyi
      VSuz1(i)= VSuzi
      Cx1(i)  = Cxi
      Cy1(i)  = Cyi
      Cz1(i)  = Czi
      tux1(i) = tuxi
      tuy1(i) = tuyi
      tuz1(i) = tuzi
      tlx1(i) = tlxi
      tly1(i) = tlyi
      tlz1(i) = tlzi
      tdx1(i) = tdxi
      tdy1(i) = tdyi
      tdz1(i) = tdzi
#endif

    end do
!$OMP END DO
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ
    TX2 = TX2 + momTempX                                 
    TY2 = TY2 + momTempY
    TZ2 = TZ2 + momTempZ

#if  TRANS == 1
   VSx2 = VSx2 + VSTempX
   VSy2 = VSy2 + VSTempY
   VSz2 = VSz2 + VSTempZ

   VSux2 = VSux2 + VSuTempX
   VSuy2 = VSuy2 + VSuTempY
   VSuz2 = VSuz2 + VSuTempZ

   VBx2 = VBx2 + VBTempX
   VBy2 = VBy2 + VBTempY
   VBz2 = VBz2 + VBTempZ

   Cx2 = Cx2 + CTempX
   Cy2 = Cy2 + CTempY
   Cz2 = Cz2 + CTempZ

   tux2 = tux2 + tuTempX
   tuy2 = tuy2 + tuTempY
   tuz2 = tuz2 + tuTempZ

   tdx2 = tdx2 + tdTempX
   tdy2 = tdy2 + tdTempY
   tdz2 = tdz2 + tdTempZ

   tlx2 = tlx2 + tlTempX
   tly2 = tly2 + tlTempY
   tlz2 = tlz2 + tlTempZ

#endif


    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

 end subroutine TPotCD_Force_Trans


!==============================================================!
!  Subroutine TPotCD_ChemicalPotential                         !
!==============================================================!

  subroutine TPotCD_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeDipole) :: this
    real(RK), pointer, contiguous      :: EPotTest(:)
    real(RK), intent(in)   :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: OX2(:), OY2(:), OZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv, RijSquared
    real(RK)          :: CosTheta
    real(RK)          :: EPotLocal
    integer           :: i, j, k, i1


    ! Assign local variables
    i1 = this%Site1%NTest
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

    ! Assign pointers
    RX1 => this%Site1%RXTest
    RY1 => this%Site1%RYTest
    RZ1 => this%Site1%RZTest
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    PX1 => this%Site1%PXTest
    PY1 => this%Site1%PYTest
    PZ1 => this%Site1%PZTest
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ

   ! Loop over test particles
   do i = 1, i1
     RXi = RX1(i)
     RYi = RY1(i)
     RZi = RZ1(i)
     PXi = PX1(i)
     PYi = PY1(i)
     PZi = PZ1(i)
     EPotLocal = 0._RK


!CDIR NODEP
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij

          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop1
          end if
          RijSquaredInv = 1._RK / RijSquared
          RijInv = sqrt( RijSquaredInv )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosTheta  = OXj * ex + OYj * eY + OZj * eZ
          EPotLocal  = EPotLocal + Epsilon * RijSquaredInv * CosTheta
        end do loop1

        EPotTest(i) = EPotTest(i) + EPotLocal

   end do

  end subroutine TPotCD_ChemicalPotential


!==============================================================!
!  Subroutine TPotCQ_Construct                                 !
!==============================================================!

  subroutine TPotCQ_Construct( this, i1, i2, j1, j2, Molecule1, Molecule2, RCutoff )
    implicit none

    ! Declare arguments
    type(TPotChargeQuadrupole)  :: this
    integer, intent(in)         :: i1, i2, j1, j2
    type(TMolecule), intent(in) :: Molecule1, Molecule2
    real(RK), intent(in)        :: RCutoff

    ! Construct potential
    this%Site1 => Molecule1%SiteCharge(j1)
    this%Site2 => Molecule2%SiteQuadrupole(j2)
    this%SameComponent = i1 == i2
    this%Epsilon = 1.5_RK * this%Site1%e * this%Site2%Q
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

  end subroutine TPotCQ_Construct


!==============================================================!
!  Subroutine TPotCQ_Destruct                                  !
!==============================================================!

  subroutine TPotCQ_Destruct( this )

    implicit none

    ! Declare arguments
    type(TPotChargeQuadrupole) :: this

    ! Destroy potential
    continue

  end subroutine TPotCQ_Destruct


!==============================================================!
!  Subroutine TPotCQ_Force                                     !
!==============================================================!

  subroutine TPotCQ_Force( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeQuadrupole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: d2EpotdV2
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1, Epsilon2
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: TX2(:), TY2(:), TZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv
    real(RK)          :: CosTheta, CosTheta2, CosAux
    real(RK)          :: EPotLocal, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: momTempX(1:this%Site2%NPart)
    real(RK)          :: momTempY(1:this%Site2%NPart)
    real(RK)          :: momTempZ(1:this%Site2%NPart)
    
#if MPI_VER > 0
    integer           :: i0
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ    
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    momTempX(:)=0._RK
    momTempY(:)=0._RK
    momTempZ(:)=0._RK    
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

!$OMP PARALLEL &
!$OMP PRIVATE ( Epsilon, Epsilon1, Epsilon2, RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (  FX1, FY1, FZ1, OX2, OY2, OZ2, TX2, TY2, TZ2) &
!$OMP PRIVATE ( Plen2, PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (   RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE (   RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE ( OXj, OYj, OZj, eX, eY, eZ, RijSquaredInv, RijInv) &
#if MPI_VER > 0
!$OMP PRIVATE ( CosTheta, CosTheta2, CosAux,  i, j, k, i1) &
!$OMP PRIVATE ( i0)
#else
!$OMP PRIVATE ( CosTheta, CosTheta2, CosAux,  i, j, k, i1)
#endif

    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    Epsilon = this%Epsilon

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ

    ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)     
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)
      FXi = FX1(i)
      FYi = FY1(i)
      FZi = FZ1(i)
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)

!CDIR NODEP
#if OSMOP == 2
loop0:do m=1,NBinsDen
        if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
          if (PXi < real(m)/NBinsDen-0.5_RK) then
            Bin1=m
            exit loop0
          end if
        end if
      end do loop0
#endif
loop1:do k = 1, this%NInCutoff(i)
        j = this%CutoffPartner(k, i)
        RXij = RXi - RX2(j)
        RYij = RYi - RY2(j)
        RZij = RZi - RZ2(j)
        PXij = PXi - PX2(j)
        PYij = PYi - PY2(j)
        PZij = PZi - PZ2(j)
        RXij = (RXij - anint( PXij )) * BoxLength                               ! Abstandsvektor von Q nach C wie bei Price
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        OXj = OX2(j)                                                            ! Orientierungsvektor Quadrupol
        OYj = OY2(j)
        OZj = OZ2(j)
        RijSquaredInv = 1._RK / ( RXij*RXij + RYij*RYij + RZij*RZij )
        RijInv = sqrt( RijSquaredInv )
        eX = RXij * RijInv                                                      ! Normierter Abstandsvektor
        eY = RYij * RijInv
        eZ = RZij * RijInv
        CosTheta  = OXj * ex + OYj * eY + OZj * eZ
        Epsilon1 = Epsilon * RijSquaredInv * RijInv
        EPotLocal  = EPotLocal + Epsilon1 * ( CosTheta * CosTheta - Third )
        CosTheta2 = 2._RK * CosTheta
        CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
        Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
        FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXj )                     ! F2 nach Price bzw. Kraft auf Punktladung
        FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYj )
        FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZj )
        VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! Vorzeichen richtig so
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
              Bin2=m 
              exit loop2
            end if
          end if
        end do loop2
        tempMin = min(Bin1, Bin2)
        tempMax = max(Bin1, Bin2)
        if(abs(PXij) .le. 0.5_RK) then
            VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
            do m = tempMin, tempMax
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        else
            VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
            do m = 1, tempMin
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
            do m = tempMax, NBinsDen
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        end if
#endif
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + Epsilon1*(CosTheta*CosTheta-Third)*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijSquaredInv)*Ninth   !xxxx3 CQ

        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        
        forceTempX(j) = forceTempX(j) - FXij
        forceTempY(j) = forceTempY(j) - FYij
        forceTempZ(j) = forceTempZ(j) - FZij
        momTempX(j) = momTempX(j) - Epsilon1 * CosTheta2 * eX   
        momTempY(j) = momTempY(j) - Epsilon1 * CosTheta2 * eY   
        momTempZ(j) = momTempZ(j) - Epsilon1 * CosTheta2 * eZ   

      end do loop1
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
    end do
!$OMP END DO
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ
    TX2 = TX2 + momTempX                                 
    TY2 = TY2 + momTempY
    TZ2 = TZ2 + momTempZ

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotCQ_Force

!==============================================================!
!  Subroutine TPotCQ_Force_Trans                               !
!==============================================================!

  subroutine TPotCQ_Force_Trans( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeQuadrupole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: d2EpotdV2
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1, Epsilon2
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: TX2(:), TY2(:), TZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv
    real(RK)          :: CosTheta, CosTheta2, CosAux
    real(RK)          :: EPotLocal, EPotLocal1, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: momTempX(1:this%Site2%NPart)
    real(RK)          :: momTempY(1:this%Site2%NPart)
    real(RK)          :: momTempZ(1:this%Site2%NPart)
 
#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer, contiguous :: VSx1(:), VSy1(:), VSz1(:)
    real(RK), pointer, contiguous :: VSx2(:), VSy2(:), VSz2(:)
    real(RK), pointer, contiguous :: VSux1(:),VSuy1(:),VSuz1(:)
    real(RK), pointer, contiguous :: VSux2(:),VSuy2(:),VSuz2(:)
    real(RK), pointer, contiguous :: VBx1(:), VBy1(:), VBz1(:)
    real(RK), pointer, contiguous :: VBx2(:), VBy2(:), VBz2(:)
    real(RK), pointer, contiguous :: Cx1(:) , Cy1(:) , Cz1(:)
    real(RK), pointer, contiguous :: Cx2(:) , Cy2(:) , Cz2(:)
    real(RK), pointer, contiguous :: tux1(:) , tuy1(:) , tuz1(:)
    real(RK), pointer, contiguous :: tux2(:) , tuy2(:) , tuz2(:)
    real(RK), pointer, contiguous :: tlx1(:) , tly1(:) , tlz1(:)
    real(RK), pointer, contiguous :: tlx2(:) , tly2(:) , tlz2(:)
    real(RK), pointer, contiguous :: tdx1(:) , tdy1(:) , tdz1(:)
    real(RK), pointer, contiguous :: tdx2(:) , tdy2(:) , tdz2(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
    real(RK)          :: VSTempX(1:this%Site2%NPart), VSTempY(1:this%Site2%NPart), VSTempZ(1:this%Site2%NPart)
    real(RK)          :: VSuTempX(1:this%Site2%NPart), VSuTempY(1:this%Site2%NPart), VSuTempZ(1:this%Site2%NPart)
    real(RK)          :: VBTempX(1:this%Site2%NPart), VBTempY(1:this%Site2%NPart), VBTempZ(1:this%Site2%NPart)
    real(RK)          :: CTempX(1:this%Site2%NPart), CTempY(1:this%Site2%NPart), CTempZ(1:this%Site2%NPart)
    real(RK)          :: tuTempX(1:this%Site2%NPart), tuTempY(1:this%Site2%NPart), tuTempZ(1:this%Site2%NPart)
    real(RK)          :: tlTempX(1:this%Site2%NPart), tlTempY(1:this%Site2%NPart), tlTempZ(1:this%Site2%NPart)
    real(RK)          :: tdTempX(1:this%Site2%NPart), tdTempY(1:this%Site2%NPart), tdTempZ(1:this%Site2%NPart)
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VBxi, VByi, VBzi
    real(RK)          :: Cxi,  Cyi,  Czi
    real(RK)          :: tuxi,  tuyi,  tuzi
    real(RK)          :: tlxi,  tlyi,  tlzi
    real(RK)          :: tdxi,  tdyi,  tdzi
    real(RK)          :: txii,  tyii , tzii
    real(RK)          :: txi ,  tyi  , tzi
    real(RK)          :: UU
    real(RK)          :: r1x, r1y, r1z
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33
    !TRANSPORT_END
#endif  
 
#if MPI_VER > 0
    integer           :: i0
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ    

    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    momTempX(:)=0._RK
    momTempY(:)=0._RK
    momTempZ(:)=0._RK    
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

#if  TRANS == 1
    VSTempX(:) = 0._RK
    VSTempY(:) = 0._RK
    VSTempZ(:) = 0._RK
    VSuTempX(:)= 0._RK
    VSuTempY(:)= 0._RK
    VSuTempZ(:)= 0._RK
    VBTempX(:) = 0._RK
    VBTempY(:) = 0._RK
    VBTempZ(:) = 0._RK
    CTempX(:)  = 0._RK
    CTempY(:)  = 0._RK
    CTempZ(:)  = 0._RK
    tuTempX(:) = 0._RK
    tuTempY(:) = 0._RK
    tuTempz(:) = 0._RK
    tlTempX(:) = 0._RK
    tlTempY(:) = 0._RK
    tlTempz(:) = 0._RK
    tdTempX(:) = 0._RK
    tdTempY(:) = 0._RK
    tdTempz(:) = 0._RK
#endif

!$OMP PARALLEL &
!$OMP PRIVATE ( Epsilon, Epsilon1, Epsilon2, RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (  FX1, FY1, FZ1, OX2, OY2, OZ2, TX2, TY2, TZ2) &
!$OMP PRIVATE ( Plen2, PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (   RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE (   RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE ( OXj, OYj, OZj, eX, eY, eZ, RijSquaredInv, RijInv) &
#if  TRANS == 1
!$OMP PRIVATE(VSx1, VSy1, VSz1,VSux1,VSuy1,VSuz1, VBx1, VBy1, VBz1, Cx1, Cy1, Cz1) &
!$OMP PRIVATE(VSx2, VSy2, VSz2,VSux2,VSuy2,VSuz2, VBx2, VBy2, VBz2, Cx2, Cy2, Cz2) &
!$OMP PRIVATE( tux1 , tuy1, tuz1, tlx1, tly1, tlz1, tdx1, tdy1, tdz1) &
!$OMP PRIVATE( tux2 , tuy2, tuz2, tlx2, tly2, tlz2, tdx2, tdy2, tdz2) &
!$OMP PRIVATE( q1, q2, q3, q4, VSxi, VSyi, VSzi, VSuxi,VSuyi,VSuzi) &
!$OMP PRIVATE( VBxi, VByi, VBzi, Cxi,  Cyi,  Czi, tuxi,  tuyi,  tuzi, tlxi,  tlyi,  tlzi) &
!$OMP PRIVATE(  tdxi,  tdyi,  tdzi, txii,  tyii , tzii, txi ,  tyi  , tzi ) &
!$OMP PRIVATE(  UU , r1x, r1y, r1z) &
!$OMP PRIVATE( A11, A12, A13, A21, A22, A23, A31, A32, A33) &
#endif
#if MPI_VER > 0
!$OMP PRIVATE ( CosTheta, CosTheta2, CosAux,  i, j, k, i1) &
!$OMP PRIVATE ( i0)
#else
!$OMP PRIVATE ( CosTheta, CosTheta2, CosAux,  i, j, k, i1)
#endif

    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    Epsilon = this%Epsilon

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ

#if  TRANS == 1
    !TRANSPORT_start
    VSx1 => this%Site1%vsCx
    VSy1 => this%Site1%vsCy
    VSz1 => this%Site1%vsCz
    VSx2 => this%Site2%vsQx
    VSy2 => this%Site2%vsQy
    VSz2 => this%Site2%vsQz
    VBx1 => this%Site1%vbCx
    VBy1 => this%Site1%vbCy
    VBz1 => this%Site1%vbCz
    VBx2 => this%Site2%vbQx
    VBy2 => this%Site2%vbQy
    VBz2 => this%Site2%vbQz
    VSux1=> this%Site1%vsuCx
    VSuy1=> this%Site1%vsuCy
    VSuz1=> this%Site1%vsuCz
    VSux2=> this%Site2%vsuQx
    VSuy2=> this%Site2%vsuQy
    VSuz2=> this%Site2%vsuQz
    Cx1  => this%Site1%cCx
    Cy1  => this%Site1%cCy
    Cz1  => this%Site1%cCz
    Cx2  => this%Site2%cQx
    Cy2  => this%Site2%cQy
    Cz2  => this%Site2%cQz
    tux1 => this%Site1%tuCx
    tuy1 => this%Site1%tuCy
    tuz1 => this%Site1%tuCz
    tux2 => this%Site2%tuQx
    tuy2 => this%Site2%tuQy
    tuz2 => this%Site2%tuQz
    tlx1 => this%Site1%tlCx
    tly1 => this%Site1%tlCy
    tlz1 => this%Site1%tlCz
    tlx2 => this%Site2%tlQx
    tly2 => this%Site2%tlQy
    tlz2 => this%Site2%tlQz
    tdx1 => this%Site1%tdCx
    tdy1 => this%Site1%tdCy
    tdz1 => this%Site1%tdCz
    tdx2 => this%Site2%tdQx
    tdy2 => this%Site2%tdQy
    tdz2 => this%Site2%tdQz
    q1  => this%Site1%Q0r(:, 1)
    q2  => this%Site1%Q0r(:, 2)
    q3  => this%Site1%Q0r(:, 3)
    q4  => this%Site1%Q0r(:, 4)
!TRANSPORT_END
#endif
 
    ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)     
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)
      FXi = FX1(i)
      FYi = FY1(i)
      FZi = FZ1(i)
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)
#if  TRANS == 1
      !TRANSPORT_start
      VSxi = VSx1(i)
      VSyi = VSy1(i)
      VSzi = VSz1(i)
      VBxi = VBx1(i)
      VByi = VBy1(i)
      VBzi = VBz1(i)
      VSuxi= VSux1(i)
      VSuyi= VSuy1(i)
      VSuzi= VSuz1(i)
      Cxi  = Cx1(i)
      Cyi  = Cy1(i)
      Czi  = Cz1(i)
      tuxi = tux1(i)
      tuyi = tuy1(i)
      tuzi = tuz1(i)
      tlxi = tlx1(i)
      tlyi = tly1(i)
      tlzi = tlz1(i)
      tdxi = tdx1(i)
      tdyi = tdy1(i)
      tdzi = tdz1(i)
      r1x  = ( RXi-PXi ) * BoxLength
      r1y  = ( RYi-PYi ) * BoxLength
      r1z  = ( RZi-PZi ) * BoxLength
      A11 = q1(i)**2 + q2(i)**2 - q3(i)**2 - q4(i)**2
      A12 = 2._RK * (q2(i) * q3(i) + q1(i) * q4(i))
      A13 = 2._RK * (q2(i) * q4(i) - q1(i) * q3(i))
      A21 = 2._RK * (q2(i) * q3(i) - q1(i) * q4(i))
      A22 = q1(i)**2 - q2(i)**2 + q3(i)**2 - q4(i)**2
      A23 = 2._RK * (q3(i) * q4(i) + q1(i) * q2(i))
      A31 = 2._RK * (q2(i) * q4(i) + q1(i) * q3(i))
      A32 = 2._RK * (q3(i) * q4(i) - q1(i) * q2(i))
      A33 = q1(i)**2 - q2(i)**2 - q3(i)**2 + q4(i)**2
      !TRANSPORT_END
#endif

!CDIR NODEP
#if OSMOP == 2
loop0:do m=1,NBinsDen
        if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
          if (PXi < real(m)/NBinsDen-0.5_RK) then
            Bin1=m
            exit loop0
          end if
        end if
      end do loop0
#endif
loop1:do k = 1, this%NInCutoff(i)
        j = this%CutoffPartner(k, i)
        RXij = RXi - RX2(j)
        RYij = RYi - RY2(j)
        RZij = RZi - RZ2(j)
        PXij = PXi - PX2(j)
        PYij = PYi - PY2(j)
        PZij = PZi - PZ2(j)
        RXij = (RXij - anint( PXij )) * BoxLength                               ! Abstandsvektor von Q nach C wie bei Price
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        OXj = OX2(j)                                                            ! Orientierungsvektor Quadrupol
        OYj = OY2(j)
        OZj = OZ2(j)
        RijSquaredInv = 1._RK / ( RXij*RXij + RYij*RYij + RZij*RZij )
        RijInv = sqrt( RijSquaredInv )
        eX = RXij * RijInv                                                      ! Normierter Abstandsvektor
        eY = RYij * RijInv
        eZ = RZij * RijInv
        CosTheta  = OXj * ex + OYj * eY + OZj * eZ
        Epsilon1 = Epsilon * RijSquaredInv * RijInv
        EPotLocal1 = Epsilon1 * ( CosTheta * CosTheta - Third ) !Gabriela: Definition of EpotLocal1 
        EPotLocal = EPotLocal + EPotlocal1
        CosTheta2 = 2._RK * CosTheta
        CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
        Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
        FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXj )                     ! F2 nach Price bzw. Kraft auf Punktladung
        FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYj )
        FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZj )
        VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! Vorzeichen richtig so
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
              Bin2=m 
              exit loop2
            end if
          end if
        end do loop2
        tempMin = min(Bin1, Bin2)
        tempMax = max(Bin1, Bin2)
        if(abs(PXij) .le. 0.5_RK) then
            VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
            do m = tempMin, tempMax
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        else
            VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
            do m = 1, tempMin
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
            do m = tempMax, NBinsDen
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        end if
#endif
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijSquaredInv)*Ninth   !xxxx3 CQ T

        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        
        forceTempX(j) = forceTempX(j) - FXij
        forceTempY(j) = forceTempY(j) - FYij
        forceTempZ(j) = forceTempZ(j) - FZij
        momTempX(j) = momTempX(j) - Epsilon1 * CosTheta2 * eX   
        momTempY(j) = momTempY(j) - Epsilon1 * CosTheta2 * eY   
        momTempZ(j) = momTempZ(j) - Epsilon1 * CosTheta2 * eZ   

#if TRANS==1
        !TRANSPORT_start vielleicht
        VSxi   = VSxi + 0.5 * FXij * PYij
        VSyi   = VSyi + 0.5 * FXij * PZij
        VSzi   = VSzi + 0.5 * FYij * PZij
        VBxi   = VBxi + 0.5 * FXij * PXij
        VByi   = VByi + 0.5 * FYij * PYij
        VBzi   = VBzi + 0.5 * FZij * PZij
        VSuxi  = VSuxi+ 0.5 * FYij * PXij
        VSuyi  = VSuyi+ 0.5 * FZij * PXij
        VSuzi  = VSuzi+ 0.5 * FZij * PYij
        VSTempX(j) = VSTempX(j) + 0.5 * FXij * PYij
        VSTempY(j) = VSTempY(j) + 0.5 * FXij * PZij
        VSTempZ(j) = VSTempZ(j) + 0.5 * FYij * PZij
        VBTempX(j) = VBTempX(j) + 0.5 * FXij * PXij
        VBTempY(j) = VBTempY(j) + 0.5 * FYij * PYij
        VBTempZ(j) = VBTempZ(j) + 0.5 * FZij * PZij
        VSuTempX(j)= VSuTempX(j)+ 0.5 * FYij * PXij
        VSuTempY(j)= VSuTempY(j)+ 0.5 * FZij * PXij
        VSuTempZ(j)= VSuTempZ(j)+ 0.5 * FZij * PYij

        UU    =  0.5 * EpotLocal1  
        Cxi    = Cxi  + UU
        Cyi    = Cyi  + UU
        Czi    = Czi  + UU
        CTempX(j) = CTempX(j) + UU
        CTempY(j) = CTempY(j) + UU
        CTempZ(j) = CTempZ(j) + UU
   

        txii   = r1y * FZij - r1z * FYij
        tyii   = r1z * FXij - r1x * FZij
        tzii   = r1x * FYij - r1y * FXij
        txi    = A11 * txii + A12 * tyii + A13 * tzii
        tyi    = A21 * txii + A22 * tyii + A23 * tzii
        tzi    = A31 * txii + A32 * tyii + A33 * tzii
        tuxi   = tuxi + 0.5 * PXij* tyi
        tuyi   = tuyi + 0.5 * PXij* tzi
        tuzi   = tuzi + 0.5 * PYij* tzi
        tlxi   = tlxi + 0.5 * PYij* txi
        tlyi   = tlyi + 0.5 * PZij* txi
        tlzi   = tlzi + 0.5 * PZij* tyi
        tdxi   = tdxi + 0.5 * PXij* txi
        tdyi   = tdyi + 0.5 * PYij* tyi
        tdzi   = tdzi + 0.5 * PZij* tzi
#endif

      end do loop1
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
#if  TRANS == 1
      !TRANSPORT_start
      VSx1(i) = VSxi
      VSy1(i) = VSyi
      VSz1(i) = VSzi
      VBx1(i) = VBxi
      VBy1(i) = VByi
      VBz1(i) = VBzi
      VSux1(i)= VSuxi
      VSuy1(i)= VSuyi
      VSuz1(i)= VSuzi
      Cx1(i)  = Cxi
      Cy1(i)  = Cyi
      Cz1(i)  = Czi
      tux1(i) = tuxi
      tuy1(i) = tuyi
      tuz1(i) = tuzi
      tlx1(i) = tlxi
      tly1(i) = tlyi
      tlz1(i) = tlzi
      tdx1(i) = tdxi
      tdy1(i) = tdyi
      tdz1(i) = tdzi
      !TRANSPORT_END
#endif  
  end do

!$OMP END DO
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ
    TX2 = TX2 + momTempX                                 
    TY2 = TY2 + momTempY
    TZ2 = TZ2 + momTempZ
   
#if  TRANS == 1
   VSx2 = VSx2 + VSTempX*BoxLength
   VSy2 = VSy2 + VSTempY*BoxLength
   VSz2 = VSz2 + VSTempZ*BoxLength

   VSux2 = VSux2 + VSuTempX*BoxLength
   VSuy2 = VSuy2 + VSuTempY*BoxLength
   VSuz2 = VSuz2 + VSuTempZ*BoxLength

   VBx2 = VBx2 + VBTempX*BoxLength
   VBy2 = VBy2 + VBTempY*BoxLength
   VBz2 = VBz2 + VBTempZ*BoxLength

   Cx2 = Cx2 + CTempX
   Cy2 = Cy2 + CTempY
   Cz2 = Cz2 + CTempZ

   tux2 = tux2 + tuTempX
   tuy2 = tuy2 + tuTempY
   tuz2 = tuz2 + tuTempZ

   tdx2 = tdx2 + tdTempX
   tdy2 = tdy2 + tdTempY
   tdz2 = tdz2 + tdTempZ

   tlx2 = tlx2 + tlTempX
   tly2 = tly2 + tlTempY
   tlz2 = tlz2 + tlTempZ

#endif
  

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotCQ_Force_Trans


!==============================================================!
!  Subroutine TPotCQ_ChemicalPotential                         !
!==============================================================!

  subroutine TPotCQ_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeQuadrupole) :: this
    real(RK), pointer, contiguous          :: EPotTest(:)
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: OX2(:), OY2(:), OZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv, RijSquared
    real(RK)          :: CosTheta
    real(RK)          :: EPotLocal
    integer           :: i, j, k, i1


    ! Assign local variables
    i1 = this%Site1%NTest
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

    ! Assign pointers
    RX1 => this%Site1%RXTest
    RY1 => this%Site1%RYTest
    RZ1 => this%Site1%RZTest
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    PX1 => this%Site1%PXTest
    PY1 => this%Site1%PYTest
    PZ1 => this%Site1%PZTest
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ

   ! Loop over test particles
   do i = 1, i1
     RXi = RX1(i)
     RYi = RY1(i)
     RZi = RZ1(i)
     PXi = PX1(i)
     PYi = PY1(i)
     PZi = PZ1(i)
     EPotLocal = 0._RK

!CDIR NODEP
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij

          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop1
          end if

          RijSquaredInv = 1._RK / RijSquared
          RijInv = sqrt( RijSquaredInv )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosTheta  = OXj * ex + OYj * eY + OZj * eZ
          EPotLocal  = EPotLocal + Epsilon * RijSquaredInv * RijInv * ( CosTheta * CosTheta - Third )
        end do loop1

        EPotTest(i) = EPotTest(i) + EPotLocal
   end do

  end subroutine TPotCQ_ChemicalPotential


!==============================================================!
!  Subroutine TPotDC_Construct                                 !
!==============================================================!

  subroutine TPotDC_Construct( this, i1, i2, j1, j2, Molecule1, Molecule2, RCutoff )

    implicit none

    ! Declare arguments
    type(TPotDipoleCharge)      :: this
    integer, intent(in)         :: i1, i2, j1, j2
    type(TMolecule), intent(in) :: Molecule1, Molecule2
    real(RK), intent(in)        :: RCutoff

    ! Construct potential
    this%Site1 => Molecule1%SiteDipole(j1)
    this%Site2 => Molecule2%SiteCharge(j2)
    this%SameComponent = i1 == i2
    this%Epsilon = this%Site1%D * this%Site2%e
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

  end subroutine TPotDC_Construct


!==============================================================!
!  Subroutine TPotDC_Destruct                                  !
!==============================================================!

  subroutine TPotDC_Destruct( this )

    implicit none

    ! Declare arguments
    type(TPotDipoleCharge) :: this

    ! Destroy potential
    continue

  end subroutine TPotDC_Destruct



!==============================================================!
!  Subroutine TPotDC_Force                                     !
!==============================================================!

  subroutine TPotDC_Force( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1, Epsilon2
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:)
    real(RK), pointer, contiguous :: TX1(:), TY1(:), TZ1(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: TXi, TYi, TZi
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv
    real(RK)          :: CosTheta, CosTheta3
    real(RK)          :: EPotLocal, Viriallocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1
    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
#if MPI_VER > 0
    integer           :: i0
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

!$OMP PARALLEL &
!$OMP PRIVATE ( Epsilon,Epsilon1, Epsilon2, RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (  FX1, FY1, FZ1, OX1, OY1, OZ1, TX1, TY1, TZ1) &
!$OMP PRIVATE ( Plen2, PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (   RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE (   RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE ( OXi, OYi, OZi,  TXi, TYi, TZi,  eX, eY, eZ) &
!$OMP PRIVATE (  RijSquaredInv, RijInv,  CosTheta, CosTheta3) &
#if MPI_VER > 0
!$OMP PRIVATE (  i, j, k, i1) &
!$OMP PRIVATE ( i0)
#else
!$OMP PRIVATE (  i, j, k, i1) 
#endif

    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    Epsilon = this%Epsilon

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ

    ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local)        
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)
      FXi = FX1(i)
      FYi = FY1(i)
      FZi = FZ1(i)
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)
      OXi = OX1(i)
      OYi = OY1(i)
      OZi = OZ1(i)
      TXi = TX1(i)
      TYi = TY1(i)
      TZi = TZ1(i)

!CDIR NODEP
#if OSMOP == 2
loop0:do m=1,NBinsDen
        if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
          if (PXi < real(m)/NBinsDen-0.5_RK) then
            Bin1=m
            exit loop0
          end if
        end if
      end do loop0
#endif
loop1:do k = 1, this%NInCutoff(i)
        j = this%CutoffPartner(k, i)
        RXij = RXi - RX2(j)
        RYij = RYi - RY2(j)
        RZij = RZi - RZ2(j)
        PXij = PXi - PX2(j)
        PYij = PYi - PY2(j)
        PZij = PZi - PZ2(j)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        RijSquaredInv = 1._RK / ( RXij*RXij + RYij*RYij + RZij*RZij )
        RijInv = sqrt( RijSquaredInv )
        eX = RXij * RijInv
        eY = RYij * RijInv
        eZ = RZij * RijInv
        CosTheta  = OXi * ex + OYi * eY + OZi * eZ                              ! -cos(alpha) bei Price
        CosTheta3 = 3._RK * CosTheta
        Epsilon1 = Epsilon * RijSquaredInv
        Epsilon2 = Epsilon1 * RijInv
        EPotLocal  = EPotLocal - Epsilon1 * CosTheta                            ! Uebereinstimmumg mit Price
        FXij = Epsilon2 * ( OXi - CosTheta3 * eX )                              ! F1 bei Price
        FYij = Epsilon2 * ( OYi - CosTheta3 * eY )
        FZij = Epsilon2 * ( OZi - CosTheta3 * eZ )
        VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F1*(-R_COM_Price); stimmt so
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
              Bin2=m 
              exit loop2
            end if
          end if
        end do loop2
        tempMin = min(Bin1, Bin2)
        tempMax = max(Bin1, Bin2)
        if(abs(PXij) .le. 0.5_RK) then
            VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
            do m = tempMin, tempMax
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        else
            VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
            do m = 1, tempMin
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
            do m = tempMax, NBinsDen
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        end if
#endif
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local - Epsilon1*CosTheta*(8._RK*sitecorr*sitecorr-2._RK*Plen2*RijSquaredInv)*Ninth          !xxxx4  DC
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        TXi = TXi + Epsilon1 * eX   ! Uebereinstimmumg mit Price; Rest bei Atom2Mol in Component
        TYi = TYi + Epsilon1 * eY   ! Reaktionsfeldbeitrag in Interaction
        TZi = TZi + Epsilon1 * eZ
        forceTempX(j) = forceTempX(j) - FXij
        forceTempY(j) = forceTempY(j) - FYij
        forceTempZ(j) = forceTempZ(j) - FZij
      end do loop1
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
      TX1(i) = TXi
      TY1(i) = TYi
      TZ1(i) = TZi
    end do
!$OMP END DO
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ
    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotDC_Force

!==============================================================!
!  Subroutine TPotDC_Force_Trans                               !
!==============================================================!

  
  subroutine TPotDC_Force_Trans( this, EPot, Virial, d2EpotdV2, BoxLength )

!Hier noch die Transportgrößen rein!

    implicit none

    ! Declare arguments
    type(TPotDipoleCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1, Epsilon2
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:)
    real(RK), pointer, contiguous :: TX1(:), TY1(:), TZ1(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: TXi, TYi, TZi
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv
    real(RK)          :: CosTheta, CosTheta3
    real(RK)          :: EPotLocal, EPotLocal1, Viriallocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1
    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
#if MPI_VER > 0
    integer           :: i0
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer, contiguous :: VSx1(:), VSy1(:), VSz1(:)
    real(RK), pointer, contiguous :: VSux1(:),VSuy1(:),VSuz1(:)
    real(RK), pointer, contiguous :: VBx1(:), VBy1(:), VBz1(:)
    real(RK), pointer, contiguous :: Cx1(:) , Cy1(:) , Cz1(:)
    real(RK), pointer, contiguous :: tux1(:) , tuy1(:) , tuz1(:)
    real(RK), pointer, contiguous :: tlx1(:) , tly1(:) , tlz1(:)
    real(RK), pointer, contiguous :: tdx1(:) , tdy1(:) , tdz1(:)
    real(RK), pointer, contiguous :: VSx2(:) , VSy2(:), VSz2(:)
    real(RK), pointer, contiguous :: VSux2(:),VSuy2(:),VSuz2(:)
    real(RK), pointer, contiguous :: VBx2(:), VBy2(:), VBz2(:)
    real(RK), pointer, contiguous :: Cx2(:) , Cy2(:) , Cz2(:)
    real(RK), pointer, contiguous :: tux2(:) , tuy2(:) , tuz2(:)
    real(RK), pointer, contiguous :: tlx2(:) , tly2(:) , tlz2(:)
    real(RK), pointer, contiguous :: tdx2(:) , tdy2(:) , tdz2(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
    real(RK)          :: VSTempX(1:this%Site2%NPart), VSTempY(1:this%Site2%NPart), VSTempZ(1:this%Site2%NPart)
    real(RK)          :: VSuTempX(1:this%Site2%NPart), VSuTempY(1:this%Site2%NPart), VSuTempZ(1:this%Site2%NPart)
    real(RK)          :: VBTempX(1:this%Site2%NPart), VBTempY(1:this%Site2%NPart), VBTempZ(1:this%Site2%NPart)
    real(RK)          :: CTempX(1:this%Site2%NPart), CTempY(1:this%Site2%NPart), CTempZ(1:this%Site2%NPart)
    real(RK)          :: tuTempX(1:this%Site2%NPart), tuTempY(1:this%Site2%NPart), tuTempZ(1:this%Site2%NPart)
    real(RK)          :: tlTempX(1:this%Site2%NPart), tlTempY(1:this%Site2%NPart), tlTempZ(1:this%Site2%NPart)
    real(RK)          :: tdTempX(1:this%Site2%NPart), tdTempY(1:this%Site2%NPart), tdTempZ(1:this%Site2%NPart)
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VBxi, VByi, VBzi
    real(RK)          :: Cxi,  Cyi,  Czi
    real(RK)          :: tuxi,  tuyi,  tuzi
    real(RK)          :: tlxi,  tlyi,  tlzi
    real(RK)          :: tdxi,  tdyi,  tdzi
    real(RK)          :: txii,  tyii , tzii
    real(RK)          :: txir , tyir , tzir
    real(RK)          :: FTXi , FTYi , FTZi
    real(RK)          :: UU
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33
 
    VSTempX(:) = 0._RK
    VSTempY(:) = 0._RK
    VSTempZ(:) = 0._RK
    VSuTempX(:)= 0._RK
    VSuTempY(:)= 0._RK
    VSuTempZ(:)= 0._RK
    VBTempX(:) = 0._RK
    VBTempY(:) = 0._RK
    VBTempZ(:) = 0._RK
    CTempX(:)  = 0._RK
    CTempY(:)  = 0._RK
    CTempZ(:)  = 0._RK
    tuTempX(:) = 0._RK
    tuTempY(:) = 0._RK
    tuTempz(:) = 0._RK
    tlTempX(:) = 0._RK
    tlTempY(:) = 0._RK
    tlTempz(:) = 0._RK
    tdTempX(:) = 0._RK
    tdTempY(:) = 0._RK
    tdTempz(:) = 0._RK
#endif



    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

!$OMP PARALLEL &
!$OMP PRIVATE ( Epsilon,Epsilon1, Epsilon2, RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (  FX1, FY1, FZ1, OX1, OY1, OZ1, TX1, TY1, TZ1) &
!$OMP PRIVATE ( Plen2,PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (   RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE (   RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE ( OXi, OYi, OZi,  TXi, TYi, TZi,  eX, eY, eZ) &
!$OMP PRIVATE (  RijSquaredInv, RijInv,  CosTheta, CosTheta3) &
#if MPI_VER > 0
!$OMP PRIVATE (  i, j, k, i1) &
!$OMP PRIVATE ( i0)
#else
!$OMP PRIVATE (  i, j, k, i1) 
#endif

    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    Epsilon = this%Epsilon

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
#if  TRANS == 1
    !TRANSPORT_start
    VSx1 => this%Site1%vsDx
    VSy1 => this%Site1%vsDy
    VSz1 => this%Site1%vsDz
    VSx2 => this%Site2%vsCx
    VSy2 => this%Site2%vsCy
    VSz2 => this%Site2%vsCz
    VBx1 => this%Site1%vbDx
    VBy1 => this%Site1%vbDy
    VBz1 => this%Site1%vbDz
    VBx2 => this%Site2%vbCx
    VBy2 => this%Site2%vbCy
    VBz2 => this%Site2%vbCz
    VSux1=> this%Site1%vsuDx
    VSuy1=> this%Site1%vsuDy
    VSuz1=> this%Site1%vsuDz
    VSux2=> this%Site2%vsuCx
    VSuy2=> this%Site2%vsuCy
    VSuz2=> this%Site2%vsuCz
    Cx1  => this%Site1%cDx
    Cy1  => this%Site1%cDy
    Cz1  => this%Site1%cDz
    Cx2  => this%Site2%cCx
    Cy2  => this%Site2%cCy
    Cz2  => this%Site2%cCz
    tux1 => this%Site1%tuDx
    tuy1 => this%Site1%tuDy
    tuz1 => this%Site1%tuDz
    tux2 => this%Site2%tuCx
    tuy2 => this%Site2%tuCy
    tuz2 => this%Site2%tuCz
    tlx1 => this%Site1%tlDx
    tly1 => this%Site1%tlDy
    tlz1 => this%Site1%tlDz
    tlx2 => this%Site2%tlCx
    tly2 => this%Site2%tlCy
    tlz2 => this%Site2%tlCz
    tdx1 => this%Site1%tdDx
    tdy1 => this%Site1%tdDy
    tdz1 => this%Site1%tdDz
    tdx2 => this%Site2%tdCx
    tdy2 => this%Site2%tdCy
    tdz2 => this%Site2%tdCz
    q1  => this%Site1%Q0r(:, 1)
    q2  => this%Site1%Q0r(:, 2)
    q3  => this%Site1%Q0r(:, 3)
    q4  => this%Site1%Q0r(:, 4)
!TRANSPORT_END
#endif

    ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local)        
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)
      FXi = FX1(i)
      FYi = FY1(i)
      FZi = FZ1(i)
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)
      OXi = OX1(i)
      OYi = OY1(i)
      OZi = OZ1(i)
      TXi = TX1(i)
      TYi = TY1(i)
      TZi = TZ1(i)
#if  TRANS == 1
        !TRANSPORT_start
        VSxi= VSx1(i)
        VSyi= VSy1(i)
        VSzi= VSz1(i)
        VBxi= VBx1(i)
        VByi= VBy1(i)
        VBzi= VBz1(i)
        VSuxi= VSux1(i)
        VSuyi= VSuy1(i)
        VSuzi= VSuz1(i)
        Cxi = Cx1(i)
        Cyi = Cy1(i)
        Czi = Cz1(i)
        tuxi = tux1(i)
        tuyi = tuy1(i)
        tuzi = tuz1(i)
        tlxi = tlx1(i)
        tlyi = tly1(i)
        tlzi = tlz1(i)
        tdxi = tdx1(i)
        tdyi = tdy1(i)
        tdzi = tdz1(i)
        A11 = q1(i)**2 + q2(i)**2 - q3(i)**2 - q4(i)**2
        A12 = 2._RK * (q2(i) * q3(i) + q1(i) * q4(i))
        A13 = 2._RK * (q2(i) * q4(i) - q1(i) * q3(i))
        A21 = 2._RK * (q2(i) * q3(i) - q1(i) * q4(i))
        A22 = q1(i)**2 - q2(i)**2 + q3(i)**2 - q4(i)**2
        A23 = 2._RK * (q3(i) * q4(i) + q1(i) * q2(i))
        A31 = 2._RK * (q2(i) * q4(i) + q1(i) * q3(i))
        A32 = 2._RK * (q3(i) * q4(i) - q1(i) * q2(i))
        A33 = q1(i)**2 - q2(i)**2 - q3(i)**2 + q4(i)**2
        !TRANSPORT_END
#endif

!CDIR NODEP
#if OSMOP == 2
loop0:do m=1,NBinsDen
        if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
          if (PXi < real(m)/NBinsDen-0.5_RK) then
            Bin1=m
            exit loop0
          end if
        end if
      end do loop0
#endif
loop1:do k = 1, this%NInCutoff(i)
        j = this%CutoffPartner(k, i)
        RXij = RXi - RX2(j)
        RYij = RYi - RY2(j)
        RZij = RZi - RZ2(j)
        PXij = PXi - PX2(j)
        PYij = PYi - PY2(j)
        PZij = PZi - PZ2(j)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        RijSquaredInv = 1._RK / ( RXij*RXij + RYij*RYij + RZij*RZij )
        RijInv = sqrt( RijSquaredInv )
        eX = RXij * RijInv
        eY = RYij * RijInv
        eZ = RZij * RijInv
        CosTheta  = OXi * ex + OYi * eY + OZi * eZ                              ! -cos(alpha) bei Price
        CosTheta3 = 3._RK * CosTheta
        Epsilon1 = Epsilon * RijSquaredInv
        Epsilon2 = Epsilon1 * RijInv
        EPotlocal1 = -Epsilon1 * CosTheta
        EPotLocal  = EPotLocal + EPotLocal1                          ! Uebereinstimmumg mit Price
        FXij = Epsilon2 * ( OXi - CosTheta3 * eX )                              ! F1 bei Price
        FYij = Epsilon2 * ( OYi - CosTheta3 * eY )
        FZij = Epsilon2 * ( OZi - CosTheta3 * eZ )
        VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F1*(-R_COM_Price); stimmt so
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
              Bin2=m 
              exit loop2
            end if
          end if
        end do loop2
        tempMin = min(Bin1, Bin2)
        tempMax = max(Bin1, Bin2)
        if(abs(PXij) .le. 0.5_RK) then
            VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
            do m = tempMin, tempMax
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        else
            VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
            do m = 1, tempMin
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
            do m = tempMax, NBinsDen
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        end if
#endif
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(8._RK*sitecorr*sitecorr-2._RK*Plen2*RijSquaredInv)*Ninth         !xxxx4  DC T
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        TXi = TXi + Epsilon1 * eX   ! Uebereinstimmumg mit Price; Rest bei Atom2Mol in Component
        TYi = TYi + Epsilon1 * eY   ! Reaktionsfeldbeitrag in Interaction
        TZi = TZi + Epsilon1 * eZ

        forceTempX(j) = forceTempX(j) - FXij
        forceTempY(j) = forceTempY(j) - FYij
        forceTempZ(j) = forceTempZ(j) - FZij

#if TRANS==1
          !TRANSPORT_start
          VSxi   = VSxi + 0.5 * FXij * PYij
          VSyi   = VSyi + 0.5 * FXij * PZij
          VSzi   = VSzi + 0.5 * FYij * PZij
          VBxi   = VBxi + 0.5 * FXij * PXij
          VByi   = VByi + 0.5 * FYij * PYij
          VBzi   = VBzi + 0.5 * FZij * PZij
          VSuxi  = VSuxi+ 0.5 * FYij * PXij
          VSuyi  = VSuyi+ 0.5 * FZij * PXij
          VSuzi  = VSuzi+ 0.5 * FZij * PYij
          VSTempX(j) = VSTempX(j) + 0.5*FXij * PYij
          VSTempY(j) = VSTempY(j) + 0.5*FXij * PZij
          VSTempZ(j) = VSTempZ(j) + 0.5*FYij * PZij
          VBTempX(j) = VBTempX(j) + 0.5*FXij * PXij
          VBTempY(j) = VBTempY(j) + 0.5*FYij * PYij
          VBTempZ(j) = VBTempZ(j) + 0.5*FZij * PZij
          VSuTempX(j)= VSuTempX(j)+ 0.5*FYij * PXij
          VSuTempY(j)= VSuTempY(j)+ 0.5*FZij * PXij
          VSuTempZ(j)= VSuTempZ(j)+ 0.5*FZij * PYij

          UU     = 0.5 * EpotLocal1
          Cxi    = Cxi  + UU
          Cyi    = Cyi  + UU
          Czi    = Czi  + UU
          CTempX(j) = CTempX(j) + UU
          CTempY(j) = CTempY(j) + UU
          CTempZ(j) = CTempZ(j) + UU

          FTXi = Epsilon1 * eX       ! Chequear                    
          FTYi = Epsilon1 * eY       !Chequear
          FTZi = Epsilon1 * eZ       !Chequear
          txii   = OYi * FTZi - OZi * FTYi   ! Chequear 
          tyii   = OZi * FTXi - OXi * FTZi   ! Chequear 
          tzii   = OXi * FTYi - OYi * FTXi   ! Chequear 
          txir   = A11 * txii + A12 * tyii + A13 * tzii
          tyir   = A21 * txii + A22 * tyii + A23 * tzii
          tzir   = A31 * txii + A32 * tyii + A33 * tzii
          tuxi   = tuxi + 0.5 * PXij* tyir
          tuyi   = tuyi + 0.5 * PXij* tzir
          tuzi   = tuzi + 0.5 * PYij* tzir
          tlxi   = tlxi + 0.5 * PYij* txir
          tlyi   = tlyi + 0.5 * PZij* txir
          tlzi   = tlzi + 0.5 * PZij* tyir
          tdxi   = tdxi + 0.5 * PXij* txir
          tdyi   = tdyi + 0.5 * PYij* tyir
          tdzi   = tdzi + 0.5 * PZij* tzir
          tuTempX(j)= tuTempX(j) + 0.5*PXij*tyi
          tuTempY(j)= tuTempY(j) + 0.5*PXij*tzi
          tuTempZ(j)= tuTempZ(j) + 0.5*PYij*tzi
          tlTempX(j)= tlTempX(j) + 0.5*PYij*txi
          tlTempY(j)= tlTempY(j) + 0.5*PZij*txi
          tlTempZ(j)= tlTempZ(j) + 0.5*PZij*tyi
          tdTempX(j)= tdTempX(j) + 0.5*PXij*txi
          tdTempY(j)= tdTempY(j) + 0.5*PYij*tyi
          tdTempZ(j)= tdTempZ(j) + 0.5*PZij*tzi
#endif

      end do loop1
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
      TX1(i) = TXi
      TY1(i) = TYi
      TZ1(i) = TZi
#if  TRANS == 1
        !TRANSPORT_start
        VSx1(i) = VSxi
        VSy1(i) = VSyi
        VSz1(i) = VSzi
        VBx1(i) = VBxi
        VBy1(i) = VByi
        VBz1(i) = VBzi
        VSux1(i)= VSuxi
        VSuy1(i)= VSuyi
        VSuz1(i)= VSuzi
        Cx1(i)  = Cxi
        Cy1(i)  = Cyi
        Cz1(i)  = Czi
        tux1(i) = tuxi
        tuy1(i) = tuyi
        tuz1(i) = tuzi
        tlx1(i) = tlxi
        tly1(i) = tlyi
        tlz1(i) = tlzi
        tdx1(i) = tdxi
        tdy1(i) = tdyi
        tdz1(i) = tdzi
#endif

    end do
!$OMP END DO
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ

#if  TRANS == 1
   VSx2 = VSx2 + VSTempX
   VSy2 = VSy2 + VSTempY
   VSz2 = VSz2 + VSTempZ

   VSux2 = VSux2 + VSuTempX
   VSuy2 = VSuy2 + VSuTempY
   VSuz2 = VSuz2 + VSuTempZ

   VBx2 = VBx2 + VBTempX
   VBy2 = VBy2 + VBTempY
   VBz2 = VBz2 + VBTempZ

   Cx2 = Cx2 + CTempX
   Cy2 = Cy2 + CTempY
   Cz2 = Cz2 + CTempZ

   tux2 = tux2 + tuTempX
   tuy2 = tuy2 + tuTempY
   tuz2 = tuz2 + tuTempZ

   tdx2 = tdx2 + tdTempX
   tdy2 = tdy2 + tdTempY
   tdz2 = tdz2 + tdTempZ

   tlx2 = tlx2 + tlTempX
   tly2 = tly2 + tlTempY
   tlz2 = tlz2 + tlTempZ

#endif


    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotDC_Force_Trans


!==============================================================!
!  Subroutine TPotDC_ChemicalPotential                         !
!==============================================================!

  subroutine TPotDC_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleCharge) :: this
    real(RK), pointer, contiguous      :: EPotTest(:)
    real(RK), intent(in)   :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv, RijSquared
    real(RK)          :: CosTheta
    real(RK)          :: EPotLocal
    integer           :: i, j, k, i1

    ! Assign local variables
    i1 = this%Site1%NTest
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

    ! Assign pointers
    RX1 => this%Site1%RXTest
    RY1 => this%Site1%RYTest
    RZ1 => this%Site1%RZTest
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    PX1 => this%Site1%PXTest
    PY1 => this%Site1%PYTest
    PZ1 => this%Site1%PZTest
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
    OX1 => this%Site1%OXTest
    OY1 => this%Site1%OYTest
    OZ1 => this%Site1%OZTest

    ! Loop over test particles
    do i = 1, i1
      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)
      OXi = OX1(i)
      OYi = OY1(i)
      OZi = OZ1(i)
      EPotLocal = 0._RK

!CDIR NODEP
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij

          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop1
          end if
          RijSquaredInv = 1._RK / RijSquared
          RijInv = sqrt( RijSquaredInv )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosTheta  = OXi * ex + OYi * eY + OZi * eZ
          EPotLocal = EPotLocal - Epsilon * RijSquaredInv * CosTheta
        end do loop1

        EPotTest(i) = EPotTest(i) + EPotLocal

    end do

  end subroutine TPotDC_ChemicalPotential


!==============================================================!
!  Subroutine TPotDD_Construct                                 !
!==============================================================!

  subroutine TPotDD_Construct( this, i1, i2, j1, j2, Molecule1, Molecule2, RCutoff, RFEpsilon )

    implicit none

    ! Declare arguments
    type(TPotDipoleDipole)      :: this
    integer, intent(in)         :: i1, i2, j1, j2
    type(TMolecule), intent(in) :: Molecule1, Molecule2
    real(RK), intent(in)        :: RCutoff
    real(RK), intent(in)        :: RFEpsilon

    ! Construct potential
    this%Site1 => Molecule1%SiteDipole(j1)
    this%Site2 => Molecule2%SiteDipole(j2)
    this%SameComponent = i1 == i2
    this%Epsilon = this%Site1%D * this%Site2%D
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2
    this%RFConstant = this%Epsilon / RCutoff**3 * (RFEpsilon - 1._RK) / (2._RK * RFEpsilon + 1._RK)

  end subroutine TPotDD_Construct


!==============================================================!
!  Subroutine TPotDD_Destruct                                  !
!==============================================================!

  subroutine TPotDD_Destruct( this )

    implicit none

    ! Declare arguments
    type(TPotDipoleDipole) :: this

    ! Destroy potential
    continue

  end subroutine TPotDD_Destruct



!==============================================================!
!  Subroutine TPotDD_Force                                     !
!==============================================================!

  subroutine TPotDD_Force( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleDipole)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RFConstant2
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: TX1(:), TY1(:), TZ1(:), TX2(:), TY2(:), TZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)

    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: TXi, TYi, TZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij3Inv, Rij4Inv3
    real(RK)          :: CosThetai, CosThetaj, CosGammaij
    real(RK)          :: CosThetai3, CosThetaj3
    real(RK)          :: Tmp
    real(RK)          :: EPotLocal, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: momTempX(1:this%Site2%NPart)
    real(RK)          :: momTempY(1:this%Site2%NPart)
    real(RK)          :: momTempZ(1:this%Site2%NPart)
    
#if MPI_VER > 0
    integer           :: N1, N2, i0, ji
    logical           :: EvenN
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ    
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    momTempX(:)=0._RK
    momTempY(:)=0._RK
    momTempZ(:)=0._RK    
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

!$OMP PARALLEL &
!$OMP PRIVATE (Epsilon,  RCutoffSquared, RFConstant2) &
!$OMP PRIVATE (RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (OX1, OY1, OZ1, OX2, OY2, OZ2) &
!$OMP PRIVATE (FX1, FY1, FZ1, TX1, TY1, TZ1) &
!$OMP PRIVATE (Plen2,sitecorr,PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (RXi, RYi, RZi, OXi, OYi, OZi,  FXi, FYi, FZi) &
!$OMP PRIVATE (TXi, TYi, TZi, PXi, PYi, PZi,  RXij, RYij, RZij) &
!$OMP PRIVATE (OXj, OYj, OZj, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE (eX, eY, eZ, RijSquared, RijInv, Rij3Inv, Rij4Inv3) &
!$OMP PRIVATE (CosThetai, CosThetaj, CosGammaij) &
!$OMP PRIVATE (CosThetai3, CosThetaj3,  Tmp) &
!$OMP PRIVATE (  SameComponent) &
#if MPI_VER > 0
!$OMP PRIVATE (i, j, k, i1, j0, j1) &
!$OMP PRIVATE ( N1, N2, i0, ji, EvenN)
#else
!$OMP PRIVATE (i, j, k, i1, j0, j1)
#endif

    ! Assign local variables
    SameComponent = this%SameComponent
#if MPI_VER > 0
    N1 = this%Site2%NPart
    N2 = N1 / 2
    EvenN = mod( N1, 2 ) == 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
    j1 = this%Site2%NPart
#endif
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RFConstant2 = 2._RK * this%RFConstant

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)      
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)

!CDIR NODEP
#if OSMOP == 2
loop0:  do m=1,NBinsDen
          if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PXi < real(m)/NBinsDen-0.5_RK) then
              Bin1=m
              exit loop0
            end if
          end if
        end do loop0
#endif
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RXij*RXij + RYij*RYij + RZij*RZij )

          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          CosThetai3 = 3._RK * CosThetai
          CosThetaj3 = 3._RK * CosThetaj
          Tmp = CosGammaij - CosThetai * CosThetaj3
          Rij3Inv = Epsilon * RijInv**3
          Rij4Inv3 = 3._RK * Rij3Inv * RijInv
          EPotLocal = EPotLocal +  Rij3Inv * Tmp

          FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                     - (eX * CosThetaj - OXj) * CosThetai)
          FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                     - (eY * CosThetaj - OYj) * CosThetai)
          FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                     - (eZ * CosThetaj - OZj) * CosThetai)

          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
                Bin2=m 
                exit loop2
              end if
            end if
          end do loop2
          tempMin = min(Bin1, Bin2)
          tempMax = max(Bin1, Bin2)
          if(abs(PXij) .le. 0.5_RK) then
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
              do m = tempMin, tempMax
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          else
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
              do m = 1, tempMin
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
              do m = tempMax, NBinsDen
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          end if
#endif
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + Rij3Inv*Tmp*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijInv*RijInv)*Ninth         !xxxx5   DD

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij
          TXi    = TXi    + Rij3Inv * (eX * CosThetaj3 - OXj)
          TYi    = TYi    + Rij3Inv * (eY * CosThetaj3 - OYj)
          TZi    = TZi    + Rij3Inv * (eZ * CosThetaj3 - OZj)
          
          momTempX(j) = momTempX(j) + Rij3Inv * (eX * CosThetai3 - OXi)  
          momTempY(j) = momTempY(j) + Rij3Inv * (eY * CosThetai3 - OYi)  
          momTempZ(j) = momTempZ(j) + Rij3Inv * (eZ * CosThetai3 - OZi)  

        end do loop1
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi

      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, merge( i1 - 1, i1, SameComponent )
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)

#if MPI_VER > 0
        if( SameComponent ) then
          j0 = i + 1
          j1 = i + N2
          if( EvenN .and. i > N2 ) j1 = j1 - 1
        else
          j0 = 1
          j1 = N1
        end if
loop3:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )

#else
        j0 = merge( i + 1, 1, SameComponent )
loop3:  do j = j0, j1
#endif

          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = (PXij - anint( RXij )) * BoxLength
          PYij = (PYij - anint( RYij )) * BoxLength
          PZij = (PZij - anint( RZij )) * BoxLength
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )

          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          CosThetai3 = 3._RK * CosThetai
          CosThetaj3 = 3._RK * CosThetaj
          Tmp = CosGammaij - CosThetai * CosThetaj3
          Rij3Inv = Epsilon * RijInv**3
          Rij4Inv3 = 3._RK * Rij3Inv * RijInv
          EPotLocal = EPotLocal + (Rij3Inv * Tmp - RFConstant2 * CosGammaij)

          FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                     - (eX * CosThetaj - OXj) * CosThetai)
          FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                     - (eY * CosThetaj - OYj) * CosThetai)
          FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                     - (eZ * CosThetaj - OZj) * CosThetai)
          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)

          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + Rij3Inv*Tmp*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijInv*RijInv)*Ninth         !xxxx5   DD ss

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij

          TXi    = TXi    + Rij3Inv * (eX * CosThetaj3 - OXj) &
&                         + RFConstant2 * OXj
          TYi    = TYi    + Rij3Inv * (eY * CosThetaj3 - OYj) &
&                         + RFConstant2 * OYj
          TZi    = TZi    + Rij3Inv * (eZ * CosThetaj3 - OZj) &
&                         + RFConstant2 * OZj

          momTempX(j) = momTempX(j) + Rij3Inv * (eX * CosThetai3 - OXi) &
&                         + RFConstant2 * OXi 
          momTempY(j) = momTempY(j) + Rij3Inv * (eY * CosThetai3 - OYi) &
&                         + RFConstant2 * OYi   
          momTempZ(j) = momTempZ(j) + Rij3Inv * (eZ * CosThetai3 - OZi) &
&                         + RFConstant2 * OZi   

        end do loop3
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do
!$OMP END DO

    end if
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ
    TX2 = TX2 + momTempX                                 
    TY2 = TY2 + momTempY
    TZ2 = TZ2 + momTempZ

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotDD_Force



!==============================================================!
!  Subroutine TPotDD_Force_Trans                               !
!==============================================================!

  subroutine TPotDD_Force_Trans( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleDipole)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RFConstant2
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: TX1(:), TY1(:), TZ1(:), TX2(:), TY2(:), TZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: TXi, TYi, TZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij3Inv, Rij4Inv3
    real(RK)          :: CosThetai, CosThetaj, CosGammaij
    real(RK)          :: CosThetai3, CosThetaj3
    real(RK)          :: Tmp, EPotLocal1
    real(RK)          :: EPotLocal, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: momTempX(1:this%Site2%NPart)
    real(RK)          :: momTempY(1:this%Site2%NPart)
    real(RK)          :: momTempZ(1:this%Site2%NPart)
    
#if MPI_VER > 0
    integer           :: N1, N2, i0, ji
    logical           :: EvenN
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer, contiguous :: VSx1(:), VSy1(:), VSz1(:)
    real(RK), pointer, contiguous :: VSx2(:), VSy2(:), VSz2(:)
    real(RK), pointer, contiguous :: VSux1(:),VSuy1(:),VSuz1(:)
    real(RK), pointer, contiguous :: VSux2(:),VSuy2(:),VSuz2(:)
    real(RK), pointer, contiguous :: VBx1(:), VBy1(:), VBz1(:)
    real(RK), pointer, contiguous :: VBx2(:), VBy2(:), VBz2(:)
    real(RK), pointer, contiguous :: Cx1(:) , Cy1(:) , Cz1(:)
    real(RK), pointer, contiguous :: Cx2(:) , Cy2(:) , Cz2(:)
    real(RK), pointer, contiguous :: tux1(:) , tuy1(:) , tuz1(:)
    real(RK), pointer, contiguous :: tux2(:) , tuy2(:) , tuz2(:)
    real(RK), pointer, contiguous :: tlx1(:) , tly1(:) , tlz1(:)
    real(RK), pointer, contiguous :: tlx2(:) , tly2(:) , tlz2(:)
    real(RK), pointer, contiguous :: tdx1(:) , tdy1(:) , tdz1(:)
    real(RK), pointer, contiguous :: tdx2(:) , tdy2(:) , tdz2(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
    real(RK)          :: VSTempX(1:this%Site2%NPart), VSTempY(1:this%Site2%NPart), VSTempZ(1:this%Site2%NPart)
    real(RK)          :: VSuTempX(1:this%Site2%NPart), VSuTempY(1:this%Site2%NPart), VSuTempZ(1:this%Site2%NPart)
    real(RK)          :: VBTempX(1:this%Site2%NPart), VBTempY(1:this%Site2%NPart), VBTempZ(1:this%Site2%NPart)
    real(RK)          :: CTempX(1:this%Site2%NPart), CTempY(1:this%Site2%NPart), CTempZ(1:this%Site2%NPart)
    real(RK)          :: tuTempX(1:this%Site2%NPart), tuTempY(1:this%Site2%NPart), tuTempZ(1:this%Site2%NPart)
    real(RK)          :: tlTempX(1:this%Site2%NPart), tlTempY(1:this%Site2%NPart), tlTempZ(1:this%Site2%NPart)
    real(RK)          :: tdTempX(1:this%Site2%NPart), tdTempY(1:this%Site2%NPart), tdTempZ(1:this%Site2%NPart)
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VBxi, VByi, VBzi
    real(RK)          :: Cxi,  Cyi,  Czi
    real(RK)          :: tuxi,  tuyi,  tuzi
    real(RK)          :: tlxi,  tlyi,  tlzi
    real(RK)          :: tdxi,  tdyi,  tdzi
    real(RK)          :: txii,  tyii , tzii
    real(RK)          :: txir , tyir , tzir
    real(RK)          :: FTXi , FTYi , FTZi
    real(RK)          :: UU
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33

    VSTempX(:) = 0._RK
    VSTempY(:) = 0._RK
    VSTempZ(:) = 0._RK
    VSuTempX(:)= 0._RK
    VSuTempY(:)= 0._RK
    VSuTempZ(:)= 0._RK
    VBTempX(:) = 0._RK
    VBTempY(:) = 0._RK
    VBTempZ(:) = 0._RK
    CTempX(:)  = 0._RK
    CTempY(:)  = 0._RK
    CTempZ(:)  = 0._RK
    tuTempX(:) = 0._RK
    tuTempY(:) = 0._RK
    tuTempz(:) = 0._RK
    tlTempX(:) = 0._RK
    tlTempY(:) = 0._RK
    tlTempz(:) = 0._RK
    tdTempX(:) = 0._RK
    tdTempY(:) = 0._RK
    tdTempz(:) = 0._RK
#endif
 

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ    

    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    momTempX(:)=0._RK
    momTempY(:)=0._RK
    momTempZ(:)=0._RK    
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

!$OMP PARALLEL &
!$OMP PRIVATE (Epsilon,  RCutoffSquared, RFConstant2) &
!$OMP PRIVATE (RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (OX1, OY1, OZ1, OX2, OY2, OZ2) &
!$OMP PRIVATE (FX1, FY1, FZ1, TX1, TY1, TZ1) &
!$OMP PRIVATE (Plen2,sitecorr,PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (RXi, RYi, RZi, OXi, OYi, OZi,  FXi, FYi, FZi) &
!$OMP PRIVATE (TXi, TYi, TZi, PXi, PYi, PZi,  RXij, RYij, RZij) &
!$OMP PRIVATE (OXj, OYj, OZj, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE (eX, eY, eZ, RijSquared, RijInv, Rij3Inv, Rij4Inv3) &
!$OMP PRIVATE (CosThetai, CosThetaj, CosGammaij) &
!$OMP PRIVATE (CosThetai3, CosThetaj3,  Tmp) &
!$OMP PRIVATE (  SameComponent) &
#if  TRANS == 1
!$OMP PRIVATE(VSx1, VSy1, VSz1, VSux1,VSuy1, VSuz1, VBx1, VBy1, VBz1, Cx1, Cy1, Cz1) &
!$OMP PRIVATE(VSx2, VSy2, VSz2, VSux2,VSuy2, VSuz2, VBx2, VBy2, VBz2, Cx2, Cy2, Cz2) &
!$OMP PRIVATE( tux1, tuy1, tuz1, tlx1, tly1, tlz1, tdx1, tdy1, tdz1) &
!$OMP PRIVATE( tux2, tuy2, tuz2, tlx2, tly2, tlz2, tdx2, tdy2, tdz2) &
!$OMP PRIVATE( q1, q2, q3, q4, VSxi, VSyi, VSzi, VSuxi,VSuyi,VSuzi) &
!$OMP PRIVATE( VBxi, VByi, VBzi, Cxi,  Cyi,  Czi, tuxi,  tuyi,  tuzi, tlxi,  tlyi,  tlzi) &
!$OMP PRIVATE(  tdxi,  tdyi,  tdzi, txii,  tyii , tzii, txir ,  tyir  , tzir ) &
!$OMP PRIVATE(   Uxi,  FTXi , FTYi , FTZi) &
!$OMP PRIVATE( A11, A12, A13, A21, A22, A23, A31, A32, A33) &
#endif
#if MPI_VER > 0
!$OMP PRIVATE ( N1, N2, i0, ji, EvenN) &
#endif
!$OMP PRIVATE (i, j, k, i1, j0, j1)

    ! Assign local variables
    SameComponent = this%SameComponent
#if MPI_VER > 0
    N1 = this%Site2%NPart
    N2 = N1 / 2
    EvenN = mod( N1, 2 ) == 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
    j1 = this%Site2%NPart
#endif
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RFConstant2 = 2._RK * this%RFConstant

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

#if  TRANS == 1
    !TRANSPORT_start
    VSx1 => this%Site1%vsDx
    VSy1 => this%Site1%vsDy
    VSz1 => this%Site1%vsDz
    VSx2 => this%Site2%vsDx
    VSy2 => this%Site2%vsDy
    VSz2 => this%Site2%vsDz
    VBx1 => this%Site1%vbDx
    VBy1 => this%Site1%vbDy
    VBz1 => this%Site1%vbDz
    VBx2 => this%Site2%vbDx
    VBy2 => this%Site2%vbDy
    VBz2 => this%Site2%vbDz
    VSux1=> this%Site1%vsuDx
    VSuy1=> this%Site1%vsuDy
    VSuz1=> this%Site1%vsuDz
    VSux2=> this%Site2%vsuDx
    VSuy2=> this%Site2%vsuDy
    VSuz2=> this%Site2%vsuDz
    Cx1  => this%Site1%cDx
    Cy1  => this%Site1%cDy
    Cz1  => this%Site1%cDz
    Cx2  => this%Site2%cDx
    Cy2  => this%Site2%cDy
    Cz2  => this%Site2%cDz
    tux1 => this%Site1%tuDx
    tuy1 => this%Site1%tuDy
    tuz1 => this%Site1%tuDz
    tux2 => this%Site2%tuDx
    tuy2 => this%Site2%tuDy
    tuz2 => this%Site2%tuDz
    tlx1 => this%Site1%tlDx
    tly1 => this%Site1%tlDy
    tlz1 => this%Site1%tlDz
    tlx2 => this%Site2%tlDx
    tly2 => this%Site2%tlDy
    tlz2 => this%Site2%tlDz
    tdx1 => this%Site1%tdDx
    tdy1 => this%Site1%tdDy
    tdz1 => this%Site1%tdDz
    tdx2 => this%Site2%tdDx
    tdy2 => this%Site2%tdDy
    tdz2 => this%Site2%tdDz
    q1  => this%Site1%Q0r(:, 1)
    q2  => this%Site1%Q0r(:, 2)
    q3  => this%Site1%Q0r(:, 3)
    q4  => this%Site1%Q0r(:, 4)
!TRANSPORT_END
#endif

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)      
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
#if  TRANS == 1
        !TRANSPORT_start
        VSxi = VSx1(i)
        VSyi = VSy1(i)
        VSzi = VSz1(i)
        VBxi = VBx1(i)
        VByi = VBy1(i)
        VBzi = VBz1(i)
        VSuxi= VSux1(i)
        VSuyi= VSuy1(i)
        VSuzi= VSuz1(i)
        Cxi  = Cx1(i)
        Cyi  = Cy1(i)
        Czi  = Cz1(i)
        tuxi = tux1(i)
        tuyi = tuy1(i)
        tuzi = tuz1(i)
        tlxi = tlx1(i)
        tlyi = tly1(i)
        tlzi = tlz1(i)
        tdxi = tdx1(i)
        tdyi = tdy1(i)
        tdzi = tdz1(i)
        A11 = q1(i)**2 + q2(i)**2 - q3(i)**2 - q4(i)**2
        A12 = 2._RK * (q2(i) * q3(i) + q1(i) * q4(i))
        A13 = 2._RK * (q2(i) * q4(i) - q1(i) * q3(i))
        A21 = 2._RK * (q2(i) * q3(i) - q1(i) * q4(i))
        A22 = q1(i)**2 - q2(i)**2 + q3(i)**2 - q4(i)**2
        A23 = 2._RK * (q3(i) * q4(i) + q1(i) * q2(i))
        A31 = 2._RK * (q2(i) * q4(i) + q1(i) * q3(i))
        A32 = 2._RK * (q3(i) * q4(i) - q1(i) * q2(i))
        A33 = q1(i)**2 - q2(i)**2 - q3(i)**2 + q4(i)**2
        !TRANSPORT_END
#endif

!CDIR NODEP
#if OSMOP == 2
loop0:  do m=1,NBinsDen
          if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PXi < real(m)/NBinsDen-0.5_RK) then
              Bin1=m
              exit loop0
            end if
          end if
        end do loop0
#endif
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RXij*RXij + RYij*RYij + RZij*RZij )

          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          CosThetai3 = 3._RK * CosThetai
          CosThetaj3 = 3._RK * CosThetaj
          Tmp = CosGammaij - CosThetai * CosThetaj3
          Rij3Inv = Epsilon * RijInv**3
          Rij4Inv3 = 3._RK * Rij3Inv * RijInv
          EPotLocal1 = Rij3Inv * Tmp
          EPotLocal = EPotLocal +  EPotLocal1
          FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                     - (eX * CosThetaj - OXj) * CosThetai)
          FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                     - (eY * CosThetaj - OYj) * CosThetai)
          FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                     - (eZ * CosThetaj - OZj) * CosThetai)

          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
                Bin2=m 
                exit loop2
              end if
            end if
          end do loop2
          tempMin = min(Bin1, Bin2)
          tempMax = max(Bin1, Bin2)
          if(abs(PXij) .le. 0.5_RK) then
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
              do m = tempMin, tempMax
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          else
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
              do m = 1, tempMin
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
              do m = tempMax, NBinsDen
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          end if
#endif
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + Rij3Inv*Tmp*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijInv*RijInv)*Ninth         !xxxx5   DD T

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij

          TXi    = TXi    + Rij3Inv * (eX * CosThetaj3 - OXj)
          TYi    = TYi    + Rij3Inv * (eY * CosThetaj3 - OYj)
          TZi    = TZi    + Rij3Inv * (eZ * CosThetaj3 - OZj)
          
          momTempX(j) = momTempX(j) + Rij3Inv * (eX * CosThetai3 - OXi)  
          momTempY(j) = momTempY(j) + Rij3Inv * (eY * CosThetai3 - OYi)  
          momTempZ(j) = momTempZ(j) + Rij3Inv * (eZ * CosThetai3 - OZi)  
#if TRANS==1
          !TRANSPORT_start
          VSxi   = VSxi + 0.5 * FXij * PYij
          VSyi   = VSyi + 0.5 * FXij * PZij
          VSzi   = VSzi + 0.5 * FYij * PZij
          VBxi   = VBxi + 0.5 * FXij * PXij
          VByi   = VByi + 0.5 * FYij * PYij
          VBzi   = VBzi + 0.5 * FZij * PZij
          VSuxi  = VSuxi+ 0.5 * FYij * PXij
          VSuyi  = VSuyi+ 0.5 * FZij * PXij
          VSuzi  = VSuzi+ 0.5 * FZij * PYij
          VSTempX(j) = VSTempX(j) + 0.5*FXij * PYij
          VSTempY(j) = VSTempY(j) + 0.5*FXij * PZij
          VSTempZ(j) = VSTempZ(j) + 0.5*FYij * PZij
          VBTempX(j) = VBTempX(j) + 0.5*FXij * PXij
          VBTempY(j) = VBTempY(j) + 0.5*FYij * PYij
          VBTempZ(j) = VBTempZ(j) + 0.5*FZij * PZij
          VSuTempX(j)= VSuTempX(j)+ 0.5*FYij * PXij
          VSuTempY(j)= VSuTempY(j)+ 0.5*FZij * PXij
          VSuTempZ(j)= VSuTempZ(j)+ 0.5*FZij * PYij

          UU        = 0.5 * (EPotLocal1 - RFConstant2 * CosGammaij)
          Cxi    = Cxi  + UU
          Cyi    = Cyi  + UU
          Czi    = Czi  + UU
          CTempX(j) = CTempX(j) + UU
          CTempY(j) = CTempY(j) + UU
          CTempZ(j) = CTempZ(j) + UU
         

          FTXi   = Rij3Inv * (eX * CosThetaj3 - OXj) + OXj * RFConstant2
          FTYi   = Rij3Inv * (eY * CosThetaj3 - OYj) + OYj * RFConstant2
          FTZi   = Rij3Inv * (eZ * CosThetaj3 - OZj) + OZj * RFConstant2
          txii   = OYi * FTZi - OZi * FTYi
          tyii   = OZi * FTXi - OXi * FTZi
          tzii   = OXi * FTYi - OYi * FTXi
          txir   = A11 * txii + A12 * tyii + A13 * tzii
          tyir   = A21 * txii + A22 * tyii + A23 * tzii
          tzir   = A31 * txii + A32 * tyii + A33 * tzii
          tuxi   = tuxi + 0.5 * PXij* tyir
          tuyi   = tuyi + 0.5 * PXij* tzir
          tuzi   = tuzi + 0.5 * PYij* tzir
          tlxi   = tlxi + 0.5 * PYij* txir
          tlyi   = tlyi + 0.5 * PZij* txir
          tlzi   = tlzi + 0.5 * PZij* tyir
          tdxi   = tdxi + 0.5 * PXij* txir
          tdyi   = tdyi + 0.5 * PYij* tyir
          tdzi   = tdzi + 0.5 * PZij* tzir
          tuTempX(j)= tuTempX(j) + 0.5 * PXij*tyi
          tuTempY(j)= tuTempY(j) + 0.5 * PXij*tzi
          tuTempZ(j)= tuTempZ(j) + 0.5 * PYij*tzi
          tlTempX(j)= tlTempX(j) + 0.5 * PYij*txi
          tlTempY(j)= tlTempY(j) + 0.5 * PZij*txi
          tlTempZ(j)= tlTempZ(j) + 0.5 * PZij*tyi
          tdTempX(j)= tdTempX(j) + 0.5 * PXij*txi
          tdTempY(j)= tdTempY(j) + 0.5 * PYij*tyi
          tdTempZ(j)= tdTempZ(j) + 0.5 * PZij*tzi
          !TRANSPORT_END
#endif
        end do loop1

        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
#if  TRANS == 1
        !TRANSPORT_start
        VSx1(i) = VSxi
        VSy1(i) = VSyi
        VSz1(i) = VSzi
        VBx1(i) = VBxi
        VBy1(i) = VByi
        VBz1(i) = VBzi
        VSux1(i)= VSuxi
        VSuy1(i)= VSuyi
        VSuz1(i)= VSuzi
        Cx1(i)  = Cxi
        Cy1(i)  = Cyi
        Cz1(i)  = Czi
        tux1(i) = tuxi
        tuy1(i) = tuyi
        tuz1(i) = tuzi
        tlx1(i) = tlxi
        tly1(i) = tlyi
        tlz1(i) = tlzi
        tdx1(i) = tdxi
        tdy1(i) = tdyi
        tdz1(i) = tdzi
        !TRANSPORT_END
#endif
      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, merge( i1 - 1, i1, SameComponent )
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
#if MPI_VER > 0
        if( SameComponent ) then
          j0 = i + 1
          j1 = i + N2
          if( EvenN .and. i > N2 ) j1 = j1 - 1
        else
          j0 = 1
          j1 = N1
        end if
loop3:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )
#else
        j0 = merge( i + 1, 1, SameComponent )
loop3:  do j = j0, j1
#endif
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = (PXij - anint( RXij )) * BoxLength
          PYij = (PYij - anint( RYij )) * BoxLength
          PZij = (PZij - anint( RZij )) * BoxLength
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )

          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          CosThetai3 = 3._RK * CosThetai
          CosThetaj3 = 3._RK * CosThetaj
          Tmp = CosGammaij - CosThetai * CosThetaj3
          Rij3Inv = Epsilon * RijInv**3
          Rij4Inv3 = 3._RK * Rij3Inv * RijInv
          EPotLocal = EPotLocal + (Rij3Inv * Tmp - RFConstant2 * CosGammaij)

          FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                     - (eX * CosThetaj - OXj) * CosThetai)
          FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                     - (eY * CosThetaj - OYj) * CosThetai)
          FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                     - (eZ * CosThetaj - OZj) * CosThetai)
          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)

          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + Rij3Inv*Tmp*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijInv*RijInv)*Ninth         !xxxx5   DD ss T

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij
          TXi    = TXi    + Rij3Inv * (eX * CosThetaj3 - OXj) &
&                         + RFConstant2 * OXj
          TYi    = TYi    + Rij3Inv * (eY * CosThetaj3 - OYj) &
&                         + RFConstant2 * OYj
          TZi    = TZi    + Rij3Inv * (eZ * CosThetaj3 - OZj) &
&                         + RFConstant2 * OZj

          momTempX(j) = momTempX(j) + Rij3Inv * (eX * CosThetai3 - OXi) &
&                         + RFConstant2 * OXi 
          momTempY(j) = momTempY(j) + Rij3Inv * (eY * CosThetai3 - OYi) &
&                         + RFConstant2 * OYi   
          momTempZ(j) = momTempZ(j) + Rij3Inv * (eZ * CosThetai3 - OZi) &
&                         + RFConstant2 * OZi   

        end do loop3
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do
!$OMP END DO

    end if
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ
    TX2 = TX2 + momTempX                                 
    TY2 = TY2 + momTempY
    TZ2 = TZ2 + momTempZ

#if  TRANS == 1
   VSx2 = VSx2 + VSTempX*BoxLength
   VSy2 = VSy2 + VSTempY*BoxLength
   VSz2 = VSz2 + VSTempZ*BoxLength

   VSux2 = VSux2 + VSuTempX*BoxLength
   VSuy2 = VSuy2 + VSuTempY*BoxLength
   VSuz2 = VSuz2 + VSuTempZ*BoxLength

   VBx2 = VBx2 + VBTempX*BoxLength
   VBy2 = VBy2 + VBTempY*BoxLength
   VBz2 = VBz2 + VBTempZ*BoxLength

   Cx2 = Cx2 + CTempX
   Cy2 = Cy2 + CTempY
   Cz2 = Cz2 + CTempZ

   tux2 = tux2 + tuTempX
   tuy2 = tuy2 + tuTempY
   tuz2 = tuz2 + tuTempZ

   tdx2 = tdx2 + tdTempX
   tdy2 = tdy2 + tdTempY
   tdz2 = tdz2 + tdTempZ

   tlx2 = tlx2 + tlTempX
   tly2 = tly2 + tlTempY
   tlz2 = tlz2 + tlTempZ

#endif

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotDD_Force_Trans

!==============================================================!
!  Subroutine TPotDD_ChemicalPotential                         !
!==============================================================!

  subroutine TPotDD_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleDipole) :: this
    real(RK), pointer, contiguous      :: EPotTest(:)
    real(RK), intent(in)   :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK)          :: RFConstant2
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij3Inv
    real(RK)          :: CosThetai, CosThetaj, CosGammaij
    real(RK)          :: Tmp
    real(RK)          :: EPotLocal
    integer           :: i, j, k, i1, j1

    ! Assign local variables
    i1 = this%Site1%NTest
    j1 = this%Site2%NPart
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared
    RFConstant2 = 2._RK * this%RFConstant

    ! Assign pointers
    RX1 => this%Site1%RXTest
    RY1 => this%Site1%RYTest
    RZ1 => this%Site1%RZTest
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    OX1 => this%Site1%OXTest
    OY1 => this%Site1%OYTest
    OZ1 => this%Site1%OZTest
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    PX1 => this%Site1%PXTest
    PY1 => this%Site1%PYTest
    PZ1 => this%Site1%PZTest
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (RXi,RYi,RZi,PXi,PYi,PZi) &
!$OMP PRIVATE (OXi, OYi, OZi, OXj, OYj, OZj) &
!$OMP PRIVATE (RXij,RYij,RZij,PXij,PYij,PZij) &
!$OMP PRIVATE (CosThetai, CosThetaj, CosGammaij) &
!$OMP PRIVATE (Tmp,RijSquared,RijInv, Rij3Inv) &
!$OMP PRIVATE (eX,eY,eZ) &
!$OMP PRIVATE (EPotLocal,i,j,k) 
    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
!$OMP DO
      do i = 1, i1
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        EPotLocal = 0._RK

!CDIR NODEP
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij

          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop1
          end if

          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )

          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          Tmp = CosGammaij - 3._RK * CosThetai * CosThetaj
          Rij3Inv = Epsilon * RijInv**3
          EPotLocal = EPotLocal + Rij3Inv * Tmp
        end do loop1

        EPotTest(i) = EPotTest(i) + EPotLocal
      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over test particles
!$OMP DO
      do i = 1, i1
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        EPotLocal = 0._RK

!CDIR NODEP
loop2:  do j = 1, j1
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij**2 + RYij**2 + RZij**2
          if( RijSquared >= RCutoffSquared ) cycle loop2

          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop2
          end if
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)
          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          Tmp = CosGammaij - 3._RK * CosThetai * CosThetaj
          Rij3Inv = Epsilon * RijInv**3
          EPotLocal = EPotLocal + Rij3Inv * Tmp - RFConstant2 * CosGammaij
        end do loop2
        EPotTest(i) = EPotTest(i) + EPotLocal
      end do
!$OMP END DO
    end if
!$OMP END PARALLEL
  end subroutine TPotDD_ChemicalPotential


!==============================================================!
!  Subroutine TPotDQ_Construct                                 !
!==============================================================!

  subroutine TPotDQ_Construct( this, i1, i2, j1, j2, Molecule1, Molecule2, RCutoff )

    implicit none

    ! Declare arguments
    type(TPotDipoleQuadrupole)  :: this
    integer, intent(in)         :: i1, i2, j1, j2
    type(TMolecule), intent(in) :: Molecule1, Molecule2
    real(RK), intent(in)        :: RCutoff

    ! Construct potential
    this%Site1 => Molecule1%SiteDipole(j1)
    this%Site2 => Molecule2%SiteQuadrupole(j2)
    this%SameComponent = i1 == i2
    this%Epsilon = 1.5_RK * this%Site1%D * this%Site2%Q
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

  end subroutine TPotDQ_Construct



!==============================================================!
!  Subroutine TPotDQ_Destruct                                  !
!==============================================================!

  subroutine TPotDQ_Destruct( this )

    implicit none

    ! Declare arguments
    type(TPotDipoleQuadrupole) :: this

    ! Destroy potential
    continue

  end subroutine TPotDQ_Destruct

!==============================================================!
!  Subroutine TPotDQ_Force                                     !
!==============================================================!

  subroutine TPotDQ_Force( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleQuadrupole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: d2EpotdV2
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: TX1(:), TY1(:), TZ1(:), TX2(:), TY2(:), TZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: TXi, TYi, TZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij4Inv
    real(RK)          :: CosThetai, CosThetaj, CosThetaj2, CosGammaij
    real(RK)          :: dCosThetai, dCosThetaj, dCosGammaij
    real(RK)          :: Tmp
    real(RK)          :: EPotLocal1, EPotLocal, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: momTempX(1:this%Site2%NPart)
    real(RK)          :: momTempY(1:this%Site2%NPart)
    real(RK)          :: momTempZ(1:this%Site2%NPart)
    
#if MPI_VER > 0
    integer           :: N1, N2, i0, ji
    logical           :: EvenN
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ    
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    momTempX(:)=0._RK
    momTempY(:)=0._RK
    momTempZ(:)=0._RK    
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

!$OMP PARALLEL &
!$OMP PRIVATE (Epsilon,  RCutoffSquared) &
!$OMP PRIVATE (RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (OX1, OY1, OZ1, OX2, OY2, OZ2) &
!$OMP PRIVATE (FX1, FY1, FZ1, TX1, TY1, TZ1) &
!$OMP PRIVATE (Plen2,sitecorr,PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (RXi, RYi, RZi, OXi, OYi, OZi,  FXi, FYi, FZi) &
!$OMP PRIVATE (TXi, TYi, TZi, PXi, PYi, PZi,  RXij, RYij, RZij) &
!$OMP PRIVATE (OXj, OYj, OZj, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE (eX, eY, eZ, RijSquared, RijInv, Rij4Inv) &
!$OMP PRIVATE (CosThetai, CosThetaj, CosThetaj2, CosGammaij) &
!$OMP PRIVATE (dCosThetai, dCosThetaj, dCosGammaij) &
!$OMP PRIVATE (Tmp, EPotLocal1) &
!$OMP PRIVATE ( SameComponent) &
#if MPI_VER > 0
!$OMP PRIVATE (i, j, k, i1, j0, j1) &
!$OMP PRIVATE ( N1, N2, i0, ji, EvenN)
#else
!$OMP PRIVATE (i, j, k, i1, j0, j1)
#endif


    ! Assign local variables
    SameComponent = this%SameComponent
#if MPI_VER > 0
    N1 = this%Site2%NPart
    N2 = N1 / 2
    EvenN = mod( N1, 2 ) == 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
    j1 = this%Site2%NPart
#endif
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)

!CDIR NODEP
#if OSMOP == 2
loop0:  do m=1,NBinsDen
          if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PXi < real(m)/NBinsDen-0.5_RK) then
              Bin1=m
              exit loop0
            end if
          end if
        end do loop0
#endif
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosThetaj2 = CosThetaj**2
          CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal1 = Rij4Inv * (CosGammaij * CosThetaj - CosThetai * (5._RK * CosThetaj2 - 1))
          EPotLocal = EPotLocal + EPotLocal1
          dCosThetai = Rij4Inv * (1 - 5._RK * CosThetaj2)
          dCosThetaj = Rij4Inv * (CosGammaij - 10._RK * CosThetai * CosThetaj)
          dCosGammaij = 2._RK * Rij4Inv * CosThetaj
          Tmp = -4._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)

          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
                Bin2=m 
                exit loop2
              end if
            end if
          end do loop2
          tempMin = min(Bin1, Bin2)
          tempMax = max(Bin1, Bin2)
          if(abs(PXij) .le. 0.5_RK) then
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
              do m = tempMin, tempMax
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          else
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
              do m = 1, tempMin
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
              do m = tempMax, NBinsDen
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          end if
#endif
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 *(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Ninth    !xxxx6   DQ

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij
          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij

          momTempX(j) = momTempX(j) - eX * dCosThetaj - OXi * dCosGammaij
          momTempY(j) = momTempY(j) - eY * dCosThetaj - OYi * dCosGammaij  
          momTempZ(j) = momTempZ(j) - eZ * dCosThetaj - OZi * dCosGammaij   

        end do loop1
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, merge( i1 - 1, i1, SameComponent )
#endif

        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)

#if MPI_VER > 0
        if( SameComponent ) then
          j0 = i + 1
          j1 = i + N2
          if( EvenN .and. i > N2 ) j1 = j1 - 1
        else
          j0 = 1
          j1 = N1
        end if
loop3:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )
#else
        j0 = merge( i + 1, 1, SameComponent )
loop3:  do j = j0, j1
#endif

          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = (PXij - anint( RXij )) * BoxLength
          PYij = (PYij - anint( RYij )) * BoxLength
          PZij = (PZij - anint( RZij )) * BoxLength
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosThetaj2 = CosThetaj**2
          CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal1 = Rij4Inv * (CosGammaij * CosThetaj - CosThetai * (5._RK * CosThetaj2 - 1))
          EPotLocal = EPotLocal + EPotLocal1
          dCosThetai = Rij4Inv * (1 - 5._RK * CosThetaj2)
          dCosThetaj = Rij4Inv * (CosGammaij - 10._RK * CosThetai * CosThetaj)
          dCosGammaij = 2._RK * Rij4Inv * CosThetaj
          Tmp = -4._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)

          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)

          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 *(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Ninth    !xxxx6   DQ ss

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij
          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
          momTempX(j) = momTempX(j) - eX * dCosThetaj - OXi * dCosGammaij
          momTempY(j) = momTempY(j) - eY * dCosThetaj - OYi * dCosGammaij  
          momTempZ(j) = momTempZ(j) - eZ * dCosThetaj - OZi * dCosGammaij   

        end do loop3
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do
!$OMP END DO

    end if
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ
    TX2 = TX2 + momTempX                                 
    TY2 = TY2 + momTempY
    TZ2 = TZ2 + momTempZ

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotDQ_Force



!==============================================================!
!  Subroutine TPotDQ_Force_Trans                               !
!==============================================================!

  subroutine TPotDQ_Force_Trans( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleQuadrupole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: d2EpotdV2
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: TX1(:), TY1(:), TZ1(:), TX2(:), TY2(:), TZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: TXi, TYi, TZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij4Inv
    real(RK)          :: CosThetai, CosThetaj, CosThetaj2, CosGammaij
    real(RK)          :: dCosThetai, dCosThetaj, dCosGammaij
    real(RK)          :: Tmp
    real(RK)          :: EPotLocal1, EPotLocal, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: momTempX(1:this%Site2%NPart)
    real(RK)          :: momTempY(1:this%Site2%NPart)
    real(RK)          :: momTempZ(1:this%Site2%NPart)
    
#if MPI_VER > 0
    integer           :: N1, N2, i0, ji
    logical           :: EvenN
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer, contiguous :: VSx1(:), VSy1(:), VSz1(:)
    real(RK), pointer, contiguous :: VSux1(:),VSuy1(:),VSuz1(:)
    real(RK), pointer, contiguous :: VBx1(:), VBy1(:), VBz1(:)
    real(RK), pointer, contiguous :: Cx1(:) , Cy1(:) , Cz1(:)
    real(RK), pointer, contiguous :: tux1(:), tuy1(:), tuz1(:)
    real(RK), pointer, contiguous :: tlx1(:), tly1(:), tlz1(:)
    real(RK), pointer, contiguous :: tdx1(:), tdy1(:), tdz1(:)
    real(RK), pointer, contiguous :: VSx2(:), VSy2(:), VSz2(:)
    real(RK), pointer, contiguous :: VSux2(:),VSuy2(:),VSuz2(:)
    real(RK), pointer, contiguous :: VBx2(:), VBy2(:), VBz2(:)
    real(RK), pointer, contiguous :: Cx2(:) , Cy2(:) , Cz2(:)
    real(RK), pointer, contiguous :: tux2(:), tuy2(:), tuz2(:)
    real(RK), pointer, contiguous :: tlx2(:), tly2(:), tlz2(:)
    real(RK), pointer, contiguous :: tdx2(:), tdy2(:), tdz2(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
    real(RK)          :: VSTempX(1:this%Site2%NPart), VSTempY(1:this%Site2%NPart), VSTempZ(1:this%Site2%NPart)
    real(RK)          :: VSuTempX(1:this%Site2%NPart), VSuTempY(1:this%Site2%NPart), VSuTempZ(1:this%Site2%NPart)
    real(RK)          :: VBTempX(1:this%Site2%NPart), VBTempY(1:this%Site2%NPart), VBTempZ(1:this%Site2%NPart)
    real(RK)          :: CTempX(1:this%Site2%NPart), CTempY(1:this%Site2%NPart), CTempZ(1:this%Site2%NPart)
    real(RK)          :: tuTempX(1:this%Site2%NPart), tuTempY(1:this%Site2%NPart), tuTempZ(1:this%Site2%NPart)
    real(RK)          :: tlTempX(1:this%Site2%NPart), tlTempY(1:this%Site2%NPart), tlTempZ(1:this%Site2%NPart)
    real(RK)          :: tdTempX(1:this%Site2%NPart), tdTempY(1:this%Site2%NPart), tdTempZ(1:this%Site2%NPart)
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VBxi, VByi, VBzi
    real(RK)          :: Cxi,  Cyi,  Czi
    real(RK)          :: tuxi,  tuyi,  tuzi
    real(RK)          :: tlxi,  tlyi,  tlzi
    real(RK)          :: tdxi,  tdyi,  tdzi
    real(RK)          :: txii,  tyii , tzii
    real(RK)          :: txir , tyir , tzir
    real(RK)          :: FTXi , FTYi , FTZi
    real(RK)          :: UU
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33

    VSTempX(:) = 0._RK
    VSTempY(:) = 0._RK
    VSTempZ(:) = 0._RK
    VSuTempX(:)= 0._RK
    VSuTempY(:)= 0._RK
    VSuTempZ(:)= 0._RK
    VBTempX(:) = 0._RK
    VBTempY(:) = 0._RK
    VBTempZ(:) = 0._RK
    CTempX(:)  = 0._RK
    CTempY(:)  = 0._RK
    CTempZ(:)  = 0._RK
    tuTempX(:) = 0._RK
    tuTempY(:) = 0._RK
    tuTempz(:) = 0._RK
    tlTempX(:) = 0._RK
    tlTempY(:) = 0._RK
    tlTempz(:) = 0._RK
    tdTempX(:) = 0._RK
    tdTempY(:) = 0._RK
    tdTempz(:) = 0._RK
#endif

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ    
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    momTempX(:)=0._RK
    momTempY(:)=0._RK
    momTempZ(:)=0._RK    
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

!$OMP PARALLEL &
!$OMP PRIVATE (Epsilon,  RCutoffSquared) &
!$OMP PRIVATE (RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (OX1, OY1, OZ1, OX2, OY2, OZ2) &
!$OMP PRIVATE (FX1, FY1, FZ1, TX1, TY1, TZ1) &
!$OMP PRIVATE (Plen2,sitecorr,PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (RXi, RYi, RZi, OXi, OYi, OZi,  FXi, FYi, FZi) &
!$OMP PRIVATE (TXi, TYi, TZi, PXi, PYi, PZi,  RXij, RYij, RZij) &
!$OMP PRIVATE (OXj, OYj, OZj, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE (eX, eY, eZ, RijSquared, RijInv, Rij4Inv) &
!$OMP PRIVATE (CosThetai, CosThetaj, CosThetaj2, CosGammaij) &
!$OMP PRIVATE (dCosThetai, dCosThetaj, dCosGammaij) &
!$OMP PRIVATE (Tmp, EPotLocal1) &
!$OMP PRIVATE ( SameComponent) &
#if  TRANS == 1
!$OMP PRIVATE(VSx1, VSy1, VSz1, VSux1,VSuy1, VSuz1, VBx1, VBy1, VBz1, Cx1, Cy1, Cz1) &
!$OMP PRIVATE(VSx2, VSy2, VSz2, VSux2,VSuy2, VSuz2, VBx2, VBy2, VBz2, Cx2, Cy2, Cz2) &
!$OMP PRIVATE( tux1, tuy1, tuz1, tlx1, tly1, tlz1, tdx1, tdy1, tdz1) &
!$OMP PRIVATE( tux2, tuy2, tuz2, tlx2, tly2, tlz2, tdx2, tdy2, tdz2) &
!$OMP PRIVATE( q1, q2, q3, q4, VSxi, VSyi, VSzi, VSuxi,VSuyi,VSuzi) &
!$OMP PRIVATE( VBxi, VByi, VBzi, Cxi,  Cyi,  Czi, tuxi,  tuyi,  tuzi, tlxi,  tlyi,  tlzi) &
!$OMP PRIVATE(  tdxi,  tdyi,  tdzi, txii,  tyii , tzii, txir ,  tyir  , tzir ) &
!$OMP PRIVATE(   UU, FTXi , FTYi , FTZi) &
!$OMP PRIVATE( A11, A12, A13, A21, A22, A23, A31, A32, A33) &
#endif
#if MPI_VER > 0
!$OMP PRIVATE (i, j, k, i1, j0, j1) &
!$OMP PRIVATE ( N1, N2, i0, ji, EvenN)
#else
!$OMP PRIVATE (i, j, k, i1, j0, j1)
#endif


    ! Assign local variables
    SameComponent = this%SameComponent
#if MPI_VER > 0
    N1 = this%Site2%NPart
    N2 = N1 / 2
    EvenN = mod( N1, 2 ) == 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
    j1 = this%Site2%NPart
#endif
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

#if  TRANS == 1
    !TRANSPORT_start
   ! Conductivity = this%Conductivity
    VSx1 => this%Site1%vsDx
    VSy1 => this%Site1%vsDy
    VSz1 => this%Site1%vsDz
    VSx2 => this%Site2%vsQx
    VSy2 => this%Site2%vsQy
    VSz2 => this%Site2%vsQz
    VBx1 => this%Site1%vbDx
    VBy1 => this%Site1%vbDy
    VBz1 => this%Site1%vbDz
    VBx2 => this%Site2%vbQx
    VBy2 => this%Site2%vbQy
    VBz2 => this%Site2%vbQz
    VSux1=> this%Site1%vsuDx
    VSuy1=> this%Site1%vsuDy
    VSuz1=> this%Site1%vsuDz
    VSux2=> this%Site2%vsuQx
    VSuy2=> this%Site2%vsuQy
    VSuz2=> this%Site2%vsuQz
    Cx1  => this%Site1%cDx
    Cy1  => this%Site1%cDy
    Cz1  => this%Site1%cDz
    Cx2  => this%Site2%cQx
    Cy2  => this%Site2%cQy
    Cz2  => this%Site2%cQz
    tux1 => this%Site1%tuDx
    tuy1 => this%Site1%tuDy
    tuz1 => this%Site1%tuDz
    tux2 => this%Site2%tuQx
    tuy2 => this%Site2%tuQy
    tuz2 => this%Site2%tuQz
    tlx1 => this%Site1%tlDx
    tly1 => this%Site1%tlDy
    tlz1 => this%Site1%tlDz
    tlx2 => this%Site2%tlQx
    tly2 => this%Site2%tlQy
    tlz2 => this%Site2%tlQz
    tdx1 => this%Site1%tdDx
    tdy1 => this%Site1%tdDy
    tdz1 => this%Site1%tdDz
    tdx2 => this%Site2%tdQx
    tdy2 => this%Site2%tdQy
    tdz2 => this%Site2%tdQz

    q1  => this%Site1%Q0r(:, 1)
    q2  => this%Site1%Q0r(:, 2)
    q3  => this%Site1%Q0r(:, 3)
    q4  => this%Site1%Q0r(:, 4)
!TRANSPORT_END
#endif

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
#if  TRANS == 1
        !TRANSPORT_start
        VSxi = VSx1(i)
        VSyi = VSy1(i)
        VSzi = VSz1(i)
        VBxi = VBx1(i)
        VByi = VBy1(i)
        VBzi = VBz1(i)
        VSuxi= VSux1(i)
        VSuyi= VSuy1(i)
        VSuzi= VSuz1(i)
        Cxi  = Cx1(i)
        Cyi  = Cy1(i)
        Czi  = Cz1(i)
        tuxi = tux1(i)
        tuyi = tuy1(i)
        tuzi = tuz1(i)
        tlxi = tlx1(i)
        tlyi = tly1(i)
        tlzi = tlz1(i)
        tdxi = tdx1(i)
        tdyi = tdy1(i)
        tdzi = tdz1(i)
        A11 = q1(i)**2 + q2(i)**2 - q3(i)**2 - q4(i)**2
        A12 = 2._RK * (q2(i) * q3(i) + q1(i) * q4(i))
        A13 = 2._RK * (q2(i) * q4(i) - q1(i) * q3(i))
        A21 = 2._RK * (q2(i) * q3(i) - q1(i) * q4(i))
        A22 = q1(i)**2 - q2(i)**2 + q3(i)**2 - q4(i)**2
        A23 = 2._RK * (q3(i) * q4(i) + q1(i) * q2(i))
        A31 = 2._RK * (q2(i) * q4(i) + q1(i) * q3(i))
        A32 = 2._RK * (q3(i) * q4(i) - q1(i) * q2(i))
        A33 = q1(i)**2 - q2(i)**2 - q3(i)**2 + q4(i)**2
        !TRANSPORT_END
#endif

!CDIR NODEP
#if OSMOP == 2
loop0:  do m=1,NBinsDen
          if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PXi < real(m)/NBinsDen-0.5_RK) then
              Bin1=m
              exit loop0
            end if
          end if
        end do loop0
#endif
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij

          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosThetaj2 = CosThetaj**2
          CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal1 = Rij4Inv * (CosGammaij * CosThetaj - CosThetai * (5._RK * CosThetaj2 - 1))
          EPotLocal = EPotLocal + EPotLocal1
          dCosThetai = Rij4Inv * (1 - 5._RK * CosThetaj2)
          dCosThetaj = Rij4Inv * (CosGammaij - 10._RK * CosThetai * CosThetaj)
          dCosGammaij = 2._RK * Rij4Inv * CosThetaj
          Tmp = -4._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)

          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
                Bin2=m 
                exit loop2
              end if
            end if
          end do loop2
          tempMin = min(Bin1, Bin2)
          tempMax = max(Bin1, Bin2)
          if(abs(PXij) .le. 0.5_RK) then
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
              do m = tempMin, tempMax
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          else
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
              do m = 1, tempMin
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
              do m = tempMax, NBinsDen
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          end if
#endif
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 *(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Ninth    !xxxx6   DQ T

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij

          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij

          momTempX(j) = momTempX(j) - eX * dCosThetaj - OXi * dCosGammaij
          momTempY(j) = momTempY(j) - eY * dCosThetaj - OYi * dCosGammaij  
          momTempZ(j) = momTempZ(j) - eZ * dCosThetaj - OZi * dCosGammaij   
#if TRANS==1
          !TRANSPORT_start
          VSxi   = VSxi + 0.5 * FXij * PYij
          VSyi   = VSyi + 0.5 * FXij * PZij
          VSzi   = VSzi + 0.5 * FYij * PZij
          VBxi   = VBxi + 0.5 * FXij * PXij
          VByi   = VByi + 0.5 * FYij * PYij
          VBzi   = VBzi + 0.5 * FZij * PZij
          VSuxi  = VSuxi+ 0.5 * FYij * PXij
          VSuyi  = VSuyi+ 0.5 * FZij * PXij
          VSuzi  = VSuzi+ 0.5 * FZij * PYij
          VSTempX(j) = VSTempX(j) + 0.5*FXij * PYij
          VSTempY(j) = VSTempY(j) + 0.5*FXij * PZij
          VSTempZ(j) = VSTempZ(j) + 0.5*FYij * PZij
          VBTempX(j) = VBTempX(j) + 0.5*FXij * PXij
          VBTempY(j) = VBTempY(j) + 0.5*FYij * PYij
          VBTempZ(j) = VBTempZ(j) + 0.5*FZij * PZij
          VSuTempX(j)= VSuTempX(j)+ 0.5*FYij * PXij
          VSuTempY(j)= VSuTempY(j)+ 0.5*FZij * PXij
          VSuTempZ(j)= VSuTempZ(j)+ 0.5*FZij * PYij

          UU     =  0.5 * EpotLocal1 
          Cxi    = Cxi  + UU
          Cyi    = Cyi  + UU
          Czi    = Czi  + UU
          CTempX(j) = CTempX(j) + UU
          CTempY(j) = CTempY(j) + UU
          CTempZ(j) = CTempZ(j) + UU
                  
          FTXi   = - eX * dCosThetai - OXj * dCosGammaij 
          FTYi   = - eY * dCosThetai - OYj * dCosGammaij  
          FTZi   = - eZ * dCosThetaj - OZi * dCosGammaij 
          txii   = OYi * FTZi - OZi * FTYi
          tyii   = OZi * FTXi - OXi * FTZi
          tzii   = OXi * FTYi - OYi * FTXi
          txir   = A11 * txii + A12 * tyii + A13 * tzii
          tyir   = A21 * txii + A22 * tyii + A23 * tzii
          tzir   = A31 * txii + A32 * tyii + A33 * tzii
          tuxi   = tuxi + 0.5 * PXij*tyir
          tuyi   = tuyi + 0.5 * PXij*tzir
          tuzi   = tuzi + 0.5 * PYij*tzir
          tlxi   = tlxi + 0.5 * PYij*txir
          tlyi   = tlyi + 0.5 * PZij*txir
          tlzi   = tlzi + 0.5 * PZij*tyir
          tdxi   = tdxi + 0.5 * PXij*txir
          tdyi   = tdyi + 0.5 * PYij*tyir
          tdzi   = tdzi + 0.5 * PZij*tzir
          tuTempX(j)= tuTempX(j) + 0.5*PXij*tyi
          tuTempY(j)= tuTempY(j) + 0.5*PXij*tzi
          tuTempZ(j)= tuTempZ(j) + 0.5*PYij*tzi
          tlTempX(j)= tlTempX(j) + 0.5*PYij*txi
          tlTempY(j)= tlTempY(j) + 0.5*PZij*txi
          tlTempZ(j)= tlTempZ(j) + 0.5*PZij*tyi
          tdTempX(j)= tdTempX(j) + 0.5*PXij*txi
          tdTempY(j)= tdTempY(j) + 0.5*PYij*tyi
          tdTempZ(j)= tdTempZ(j) + 0.5*PZij*tzi
          !TRANSPORT_END
#endif
        end do loop1

        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
#if  TRANS == 1
        !TRANSPORT_start
        VSx1(i) = VSxi
        VSy1(i) = VSyi
        VSz1(i) = VSzi
        VBx1(i) = VBxi
        VBy1(i) = VByi
        VBz1(i) = VBzi
        VSux1(i)= VSuxi
        VSuy1(i)= VSuyi
        VSuz1(i)= VSuzi
        Cx1(i)  = Cxi
        Cy1(i)  = Cyi
        Cz1(i)  = Czi
        tux1(i) = tuxi
        tuy1(i) = tuyi
        tuz1(i) = tuzi
        tlx1(i) = tlxi
        tly1(i) = tlyi
        tlz1(i) = tlzi
        tdx1(i) = tdxi
        tdy1(i) = tdyi
        tdz1(i) = tdzi
        !TRANSPORT_END
#endif      
       end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, merge( i1 - 1, i1, SameComponent )
#endif
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)

#if MPI_VER > 0
        if( SameComponent ) then
          j0 = i + 1
          j1 = i + N2
          if( EvenN .and. i > N2 ) j1 = j1 - 1
        else
          j0 = 1
          j1 = N1
        end if
loop3:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )
#else
        j0 = merge( i + 1, 1, SameComponent )
loop3:  do j = j0, j1
#endif

          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = (PXij - anint( RXij )) * BoxLength
          PYij = (PYij - anint( RYij )) * BoxLength
          PZij = (PZij - anint( RZij )) * BoxLength
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosThetaj2 = CosThetaj**2
          CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal1 = Rij4Inv * (CosGammaij * CosThetaj - CosThetai * (5._RK * CosThetaj2 - 1))
          EPotLocal = EPotLocal + EPotLocal1
          dCosThetai = Rij4Inv * (1 - 5._RK * CosThetaj2)
          dCosThetaj = Rij4Inv * (CosGammaij - 10._RK * CosThetai * CosThetaj)
          dCosGammaij = 2._RK * Rij4Inv * CosThetaj
          Tmp = -4._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)

          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)

          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 *(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Ninth    !xxxx6   DQ ss T

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij
          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
          momTempX(j) = momTempX(j) - eX * dCosThetaj - OXi * dCosGammaij
          momTempY(j) = momTempY(j) - eY * dCosThetaj - OYi * dCosGammaij  
          momTempZ(j) = momTempZ(j) - eZ * dCosThetaj - OZi * dCosGammaij   

        end do loop3
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do
!$OMP END DO

    end if
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ
    TX2 = TX2 + momTempX                                 
    TY2 = TY2 + momTempY
    TZ2 = TZ2 + momTempZ

#if  TRANS == 1
   VSx2 = VSx2 + VSTempX
   VSy2 = VSy2 + VSTempY
   VSz2 = VSz2 + VSTempZ

   VSux2 = VSux2 + VSuTempX
   VSuy2 = VSuy2 + VSuTempY
   VSuz2 = VSuz2 + VSuTempZ

   VBx2 = VBx2 + VBTempX
   VBy2 = VBy2 + VBTempY
   VBz2 = VBz2 + VBTempZ

   Cx2 = Cx2 + CTempX
   Cy2 = Cy2 + CTempY
   Cz2 = Cz2 + CTempZ

   tux2 = tux2 + tuTempX
   tuy2 = tuy2 + tuTempY
   tuz2 = tuz2 + tuTempZ

   tdx2 = tdx2 + tdTempX
   tdy2 = tdy2 + tdTempY
   tdz2 = tdz2 + tdTempZ

   tlx2 = tlx2 + tlTempX
   tly2 = tly2 + tlTempY
   tlz2 = tlz2 + tlTempZ

#endif
 

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotDQ_Force_Trans



!==============================================================!
!  Subroutine TPotDQ_ChemicalPotential                         !
!==============================================================!

  subroutine TPotDQ_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleQuadrupole) :: this
    real(RK), pointer, contiguous          :: EPotTest(:)
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij4Inv
    real(RK)          :: CosThetai, CosThetaj, CosGammaij
    real(RK)          :: EPotLocal
    integer           :: i, j, k, i1, j1

    ! Assign local variables
    i1 = this%Site1%NTest
    j1 = this%Site2%NPart
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

    ! Assign pointers
    RX1 => this%Site1%RXTest
    RY1 => this%Site1%RYTest
    RZ1 => this%Site1%RZTest
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    OX1 => this%Site1%OXTest
    OY1 => this%Site1%OYTest
    OZ1 => this%Site1%OZTest
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    PX1 => this%Site1%PXTest
    PY1 => this%Site1%PYTest
    PZ1 => this%Site1%PZTest
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (RXi,RYi,RZi,PXi,PYi,PZi) &
!$OMP PRIVATE (OXi, OYi, OZi, OXj, OYj, OZj) &
!$OMP PRIVATE (RXij,RYij,RZij,PXij,PYij,PZij) &
!$OMP PRIVATE (CosThetai, CosThetaj, CosGammaij) &
!$OMP PRIVATE (RijSquared,RijInv, Rij4Inv) &
!$OMP PRIVATE (eX,eY,eZ) &
!$OMP PRIVATE (EPotLocal,i,j,k) 
    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
!$OMP DO
      do i = 1, i1
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        EPotLocal = 0._RK

!CDIR NODEP

loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij

          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop1
          end if
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal = EPotLocal + Rij4Inv * (2._RK * CosGammaij * CosThetaj - CosThetai * (5._RK * CosThetaj**2 - 1._RK))
        end do loop1

        EPotTest(i) = EPotTest(i) + EPotLocal
      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over test particles
!$OMP DO
      do i = 1, i1
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        EPotLocal = 0._RK

!CDIR NODEP
loop2:  do j = 1, j1
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          if( RijSquared >= RCutoffSquared ) cycle loop2
          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop2
          end if
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal = EPotLocal + Rij4Inv * (2._RK * CosGammaij * CosThetaj - CosThetai * (5._RK * CosThetaj**2 - 1._RK))
        end do loop2

        EPotTest(i) = EPotTest(i) + EPotLocal
      end do
!$OMP END DO
    end if
!$OMP END PARALLEL
  end subroutine TPotDQ_ChemicalPotential


!==============================================================!
!  Subroutine TPotQC_Construct                                 !
!==============================================================!

  subroutine TPotQC_Construct( this, i1, i2, j1, j2, Molecule1, Molecule2, RCutoff )
    implicit none

    ! Declare arguments
    type(TPotQuadrupoleCharge)  :: this
    integer, intent(in)         :: i1, i2, j1, j2
    type(TMolecule), intent(in) :: Molecule1, Molecule2
    real(RK), intent(in)        :: RCutoff

    ! Construct potential
    this%Site1 => Molecule1%SiteQuadrupole(j1)
    this%Site2 => Molecule2%SiteCharge(j2)
    this%SameComponent = i1 == i2
    this%Epsilon = 1.5_RK * this%Site1%Q * this%Site2%e
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

  end subroutine TPotQC_Construct



!==============================================================!
!  Subroutine TPotQC_Destruct                                  !
!==============================================================!

  subroutine TPotQC_Destruct( this )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleCharge) :: this

    ! Destroy potential
    continue

  end subroutine TPotQC_Destruct



!==============================================================!
!  Subroutine TPotQC_Force                                     !
!==============================================================!

  subroutine TPotQC_Force( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleCharge) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: d2EpotdV2
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1, Epsilon2
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:)
    real(RK), pointer, contiguous :: TX1(:), TY1(:), TZ1(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: TXi, TYi, TZi
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv
    real(RK)          :: CosTheta, CosTheta2, CosAux
    real(RK)          :: EPotLocal, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    
#if MPI_VER > 0
    integer           :: i0
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
  
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK   
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

!$OMP PARALLEL &
!$OMP PRIVATE ( Epsilon, Epsilon1, Epsilon2, RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (  FX1, FY1, FZ1, OX1, OY1, OZ1, TX1, TY1, TZ1) &
!$OMP PRIVATE (Plen2,sitecorr, PX1, PY1, PZ1, PX2, PY2, PZ2, TXi, TYi, TZi) &
!$OMP PRIVATE (   RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE (   RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE (  eX, eY, eZ, RijSquaredInv, RijInv) &
#if MPI_VER > 0
!$OMP PRIVATE ( CosTheta, CosTheta2, CosAux,  i, j, k, i1) &
!$OMP PRIVATE ( i0)
#else
!$OMP PRIVATE ( CosTheta, CosTheta2, CosAux,  i, j, k, i1)
#endif


    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    Epsilon = this%Epsilon

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

    ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local)    
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif

      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)
      FXi = FX1(i)
      FYi = FY1(i)
      FZi = FZ1(i)
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)
      OXi = OX1(i)
      OYi = OY1(i)
      OZi = OZ1(i)
      TXi = TX1(i)
      TYi = TY1(i)
      TZi = TZ1(i)
!CDIR NODEP
#if OSMOP == 2
loop0:do m=1,NBinsDen
        if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
          if (PXi < real(m)/NBinsDen-0.5_RK) then
            Bin1=m
            exit loop0
          end if
        end if
      end do loop0
#endif
loop1:do k = 1, this%NInCutoff(i)
        j = this%CutoffPartner(k, i)
        RXij = RXi - RX2(j)
        RYij = RYi - RY2(j)
        RZij = RZi - RZ2(j)
        PXij = PXi - PX2(j)
        PYij = PYi - PY2(j)
        PZij = PZi - PZ2(j)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        RijSquaredInv = 1._RK / ( RXij*RXij + RYij*RYij + RZij*RZij )
        RijInv = sqrt( RijSquaredInv )
        eX = - RXij * RijInv                          ! Normierter Abstandsvektor nach Price
        eY = - RYij * RijInv
        eZ = - RZij * RijInv
        CosTheta  = OXi * ex + OYi * eY + OZi * eZ          ! Scalarprodukt normierter 
!                                              Abstandsvektor mit Orientierungsvektor Quadrupol
        Epsilon1 = Epsilon * RijSquaredInv * RijInv
        EPotLocal  = EPotLocal + Epsilon1 * ( CosTheta * CosTheta - Third )
        CosTheta2 = 2._RK * CosTheta
        CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
        Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
        FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXi ) ! Kraft auf die Punktladung, sprich F2
        FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYi )
        FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZi )
        VirialLocal = VirialLocal - (FXij * PXij + FYij * PYij + FZij * PZij)     ! Vorzeichen richtig
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
              Bin2=m 
              exit loop2
            end if
          end if
        end do loop2
        tempMin = min(Bin1, Bin2)
        tempMax = max(Bin1, Bin2)
        if(abs(PXij) .le. 0.5_RK) then
            VirialPart = -(FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
            do m = tempMin, tempMax
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        else
            VirialPart = -(FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
            do m = 1, tempMin
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
            do m = tempMax, NBinsDen
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        end if
#endif
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + Epsilon1*(CosTheta*CosTheta-Third)*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijSquaredInv)*Ninth    !xxxx7  QC

        FXi    = FXi    - FXij
        FYi    = FYi    - FYij
        FZi    = FZi    - FZij

        forceTempX(j) = forceTempX(j) + FXij
        forceTempY(j) = forceTempY(j) + FYij
        forceTempZ(j) = forceTempZ(j) + FZij

        TXi    = TXi - Epsilon1*CosTheta2*eX  
        ! Drehmomentanteil auf Quadrupol wegen Punktladung. Kreuzprodukt
        TYi    = TYi - Epsilon1*CosTheta2*eY  ! in Atom2Mol von Component
        TZi    = TZi - Epsilon1*CosTheta2*eZ
      end do loop1

      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
      TX1(i) = TXi
      TY1(i) = TYi
      TZ1(i) = TZi
    end do
!$OMP END DO
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotQC_Force

!==============================================================!
!  Subroutine TPotQC_Force_Trans                               !
!==============================================================!

  subroutine TPotQC_Force_Trans( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleCharge) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: d2EpotdV2
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1, Epsilon2
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:)
    real(RK), pointer, contiguous :: TX1(:), TY1(:), TZ1(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: TXi, TYi, TZi
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv
    real(RK)          :: CosTheta, CosTheta2, CosAux
    real(RK)          :: EPotLocal, EPotLocal1, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    
#if MPI_VER > 0
    integer           :: i0
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer, contiguous :: VSx1(:), VSy1(:), VSz1(:)
    real(RK), pointer, contiguous :: VSux1(:),VSuy1(:),VSuz1(:)
    real(RK), pointer, contiguous :: VBx1(:), VBy1(:), VBz1(:)
    real(RK), pointer, contiguous :: Cx1(:) , Cy1(:) , Cz1(:)
    real(RK), pointer, contiguous :: tux1(:) , tuy1(:) , tuz1(:)
    real(RK), pointer, contiguous :: tlx1(:) , tly1(:) , tlz1(:)
    real(RK), pointer, contiguous :: tdx1(:) , tdy1(:) , tdz1(:)
    real(RK), pointer, contiguous :: VSx2(:), VSy2(:), VSz2(:)
    real(RK), pointer, contiguous :: VSux2(:),VSuy2(:),VSuz2(:)
    real(RK), pointer, contiguous :: VBx2(:), VBy2(:), VBz2(:)
    real(RK), pointer, contiguous :: Cx2(:) , Cy2(:) , Cz2(:)
    real(RK), pointer, contiguous :: tux2(:) , tuy2(:) , tuz2(:)
    real(RK), pointer, contiguous :: tlx2(:) , tly2(:) , tlz2(:)
    real(RK), pointer, contiguous :: tdx2(:) , tdy2(:) , tdz2(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
    real(RK)          :: VSTempX(1:this%Site2%NPart), VSTempY(1:this%Site2%NPart), VSTempZ(1:this%Site2%NPart)
    real(RK)          :: VSuTempX(1:this%Site2%NPart), VSuTempY(1:this%Site2%NPart), VSuTempZ(1:this%Site2%NPart)
    real(RK)          :: VBTempX(1:this%Site2%NPart), VBTempY(1:this%Site2%NPart), VBTempZ(1:this%Site2%NPart)
    real(RK)          :: CTempX(1:this%Site2%NPart), CTempY(1:this%Site2%NPart), CTempZ(1:this%Site2%NPart)
    real(RK)          :: tuTempX(1:this%Site2%NPart), tuTempY(1:this%Site2%NPart), tuTempZ(1:this%Site2%NPart)
    real(RK)          :: tlTempX(1:this%Site2%NPart), tlTempY(1:this%Site2%NPart), tlTempZ(1:this%Site2%NPart)
    real(RK)          :: tdTempX(1:this%Site2%NPart), tdTempY(1:this%Site2%NPart), tdTempZ(1:this%Site2%NPart)
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VBxi, VByi, VBzi
    real(RK)          :: Cxi,  Cyi,  Czi
    real(RK)          :: tuxi,  tuyi,  tuzi
    real(RK)          :: tlxi,  tlyi,  tlzi
    real(RK)          :: tdxi,  tdyi,  tdzi
    real(RK)          :: txii,  tyii , tzii
    real(RK)          :: txir , tyir , tzir
    real(RK)          :: FTXi , FTYi , FTZi
    real(RK)          :: UU
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33

    VSTempX(:) = 0._RK
    VSTempY(:) = 0._RK
    VSTempZ(:) = 0._RK
    VSuTempX(:)= 0._RK
    VSuTempY(:)= 0._RK
    VSuTempZ(:)= 0._RK
    VBTempX(:) = 0._RK
    VBTempY(:) = 0._RK
    VBTempZ(:) = 0._RK
    CTempX(:)  = 0._RK
    CTempY(:)  = 0._RK
    CTempZ(:)  = 0._RK
    tuTempX(:) = 0._RK
    tuTempY(:) = 0._RK
    tuTempz(:) = 0._RK
    tlTempX(:) = 0._RK
    tlTempY(:) = 0._RK
    tlTempz(:) = 0._RK
    tdTempX(:) = 0._RK
    tdTempY(:) = 0._RK
    tdTempz(:) = 0._RK
#endif


    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
  
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK   
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

!$OMP PARALLEL &
!$OMP PRIVATE ( Epsilon, Epsilon1, Epsilon2, RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (  FX1, FY1, FZ1, OX1, OY1, OZ1, TX1, TY1, TZ1) &
!$OMP PRIVATE (Plen2,sitecorr, PX1, PY1, PZ1, PX2, PY2, PZ2, TXi, TYi, TZi) &
!$OMP PRIVATE (   RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE (   RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE (  eX, eY, eZ, RijSquaredInv, RijInv) &
#if  TRANS == 1
!$OMP PRIVATE(VSx1, VSy1, VSz1,VSux1,VSuy1,VSuz1, VBx1, VBy1, VBz1, Cx1, Cy1, Cz1) &
!$OMP PRIVATE(VSx2, VSy2, VSz2,VSux2,VSuy2,VSuz2, VBx2, VBy2, VBz2, Cx2, Cy2, Cz2) &
!$OMP PRIVATE( tux1, tuy1, tuz1, tlx1, tly1, tlz1, tdx1, tdy1, tdz1) &
!$OMP PRIVATE( tux2, tuy2, tuz2, tlx2, tly2, tlz2, tdx2, tdy2, tdz2) &
!$OMP PRIVATE( q1, q2, q3, q4, VSxi, VSyi, VSzi, VSuxi,VSuyi,VSuzi) &
!$OMP PRIVATE( VBxi, VByi, VBzi, Cxi,  Cyi,  Czi, tuxi,  tuyi,  tuzi, tlxi,  tlyi,  tlzi) &
!$OMP PRIVATE(  tdxi,  tdyi,  tdzi, txii,  tyii , tzii, txir ,  tyir  , tzir ) &
!$OMP PRIVATE(   UU, FTXi , FTYi , FTZi) &
!$OMP PRIVATE( A11, A12, A13, A21, A22, A23, A31, A32, A33) &
#endif
#if MPI_VER > 0
!$OMP PRIVATE ( CosTheta, CosTheta2, CosAux,  i, j, k, i1) &
!$OMP PRIVATE ( i0)
#else
!$OMP PRIVATE ( CosTheta, CosTheta2, CosAux,  i, j, k, i1)
#endif


    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    Epsilon = this%Epsilon

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

#if  TRANS == 1
    !TRANSPORT_start
    VSx1 => this%Site1%vsQx
    VSy1 => this%Site1%vsQy
    VSz1 => this%Site1%vsQz
    VSx2 => this%Site2%vsCx
    VSy2 => this%Site2%vsCy
    VSz2 => this%Site2%vsCz
    VBx1 => this%Site1%vbQx
    VBy1 => this%Site1%vbQy
    VBz1 => this%Site1%vbQz
    VBx2 => this%Site2%vbCx
    VBy2 => this%Site2%vbCy
    VBz2 => this%Site2%vbCz
    VSux1=> this%Site1%vsuQx
    VSuy1=> this%Site1%vsuQy
    VSuz1=> this%Site1%vsuQz
    VSux2=> this%Site2%vsuCx
    VSuy2=> this%Site2%vsuCy
    VSuz2=> this%Site2%vsuCz
    Cx1  => this%Site1%cQx
    Cy1  => this%Site1%cQy
    Cz1  => this%Site1%cQz
    Cx2  => this%Site2%cCx
    Cy2  => this%Site2%cCy
    Cz2  => this%Site2%cCz
    tux1 => this%Site1%tuQx
    tuy1 => this%Site1%tuQy
    tuz1 => this%Site1%tuQz
    tux2 => this%Site2%tuCx
    tuy2 => this%Site2%tuCy
    tuz2 => this%Site2%tuCz
    tlx1 => this%Site1%tlQx
    tly1 => this%Site1%tlQy
    tlz1 => this%Site1%tlQz
    tlx2 => this%Site2%tlCx
    tly2 => this%Site2%tlCy
    tlz2 => this%Site2%tlCz
    tdx1 => this%Site1%tdQx
    tdy1 => this%Site1%tdQy
    tdz1 => this%Site1%tdQz
    tdx2 => this%Site2%tdCx
    tdy2 => this%Site2%tdCy
    tdz2 => this%Site2%tdCz
    q1  => this%Site1%Q0r(:, 1)
    q2  => this%Site1%Q0r(:, 2)
    q3  => this%Site1%Q0r(:, 3)
    q4  => this%Site1%Q0r(:, 4)
!TRANSPORT_END
#endif
    ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local)    

#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)
      FXi = FX1(i)
      FYi = FY1(i)
      FZi = FZ1(i)
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)
      OXi = OX1(i)
      OYi = OY1(i)
      OZi = OZ1(i)
      TXi = TX1(i)
      TYi = TY1(i)
      TZi = TZ1(i)

!CDIR NODEP
#if  TRANS == 1
        !TRANSPORT_start
        VSxi = VSx1(i)
        VSyi = VSy1(i)
        VSzi = VSz1(i)
        VBxi = VBx1(i)
        VByi = VBy1(i)
        VBzi = VBz1(i)
        VSuxi= VSux1(i)
        VSuyi= VSuy1(i)
        VSuzi= VSuz1(i)
        Cxi  = Cx1(i)
        Cyi  = Cy1(i)
        Czi  = Cz1(i)
        tuxi = tux1(i)
        tuyi = tuy1(i)
        tuzi = tuz1(i)
        tlxi = tlx1(i)
        tlyi = tly1(i)
        tlzi = tlz1(i)
        tdxi = tdx1(i)
        tdyi = tdy1(i)
        tdzi = tdz1(i)
        A11 = q1(i)**2 + q2(i)**2 - q3(i)**2 - q4(i)**2
        A12 = 2._RK * (q2(i) * q3(i) + q1(i) * q4(i))
        A13 = 2._RK * (q2(i) * q4(i) - q1(i) * q3(i))
        A21 = 2._RK * (q2(i) * q3(i) - q1(i) * q4(i))
        A22 = q1(i)**2 - q2(i)**2 + q3(i)**2 - q4(i)**2
        A23 = 2._RK * (q3(i) * q4(i) + q1(i) * q2(i))
        A31 = 2._RK * (q2(i) * q4(i) + q1(i) * q3(i))
        A32 = 2._RK * (q3(i) * q4(i) - q1(i) * q2(i))
        A33 = q1(i)**2 - q2(i)**2 - q3(i)**2 + q4(i)**2
#endif

#if OSMOP == 2
loop0:do m=1,NBinsDen
        if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
          if (PXi < real(m)/NBinsDen-0.5_RK) then
            Bin1=m
            exit loop0
          end if
        end if
      end do loop0
#endif
loop1:do k = 1, this%NInCutoff(i)
        j = this%CutoffPartner(k, i)
        RXij = RXi - RX2(j)
        RYij = RYi - RY2(j)
        RZij = RZi - RZ2(j)
        PXij = PXi - PX2(j)
        PYij = PYi - PY2(j)
        PZij = PZi - PZ2(j)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        RijSquaredInv = 1._RK / ( RXij*RXij + RYij*RYij + RZij*RZij )
        RijInv = sqrt( RijSquaredInv )
        eX = - RXij * RijInv                          ! Normierter Abstandsvektor nach Price
        eY = - RYij * RijInv
        eZ = - RZij * RijInv
        CosTheta  = OXi * ex + OYi * eY + OZi * eZ          ! Scalarprodukt normierter 
!                                              Abstandsvektor mit Orientierungsvektor Quadrupol
        Epsilon1 = Epsilon * RijSquaredInv * RijInv
        EPotLocal1 = Epsilon1 * ( CosTheta * CosTheta - Third )
        EPotLocal  = EPotLocal + EPotLocal1
        CosTheta2 = 2._RK * CosTheta
        CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
        Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
        FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXi ) ! Kraft auf die Punktladung, sprich F2
        FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYi )
        FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZi )
        VirialLocal = VirialLocal - (FXij * PXij + FYij * PYij + FZij * PZij)     ! Vorzeichen richtig
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
              Bin2=m 
              exit loop2
            end if
          end if
        end do loop2
        tempMin = min(Bin1, Bin2)
        tempMax = max(Bin1, Bin2)
        if(abs(PXij) .le. 0.5_RK) then
            VirialPart = -(FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
            do m = tempMin, tempMax
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        else
            VirialPart = -(FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
            do m = 1, tempMin
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
            do m = tempMax, NBinsDen
              this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
            end do
        end if
#endif
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + Epsilon1*(CosTheta*CosTheta-Third)*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijSquaredInv)*Ninth    !xxxx7  QC T
        FXi    = FXi    - FXij
        FYi    = FYi    - FYij
        FZi    = FZi    - FZij

        forceTempX(j) = forceTempX(j) + FXij
        forceTempY(j) = forceTempY(j) + FYij
        forceTempZ(j) = forceTempZ(j) + FZij

        TXi    = TXi - Epsilon1*CosTheta2*eX  
        ! Drehmomentanteil auf Quadrupol wegen Punktladung. Kreuzprodukt
        TYi    = TYi - Epsilon1*CosTheta2*eY  ! in Atom2Mol von Component
        TZi    = TZi - Epsilon1*CosTheta2*eZ
#if  TRANS == 1
!TRANSPORT_start
          VSxi   = VSxi + 0.5 * FXij * PYij
          VSyi   = VSyi + 0.5 * FXij * PZij
          VSzi   = VSzi + 0.5 * FYij * PZij
          VBxi   = VBxi + 0.5 * FXij * PXij
          VByi   = VByi + 0.5 * FYij * PYij
          VBzi   = VBzi + 0.5 * FZij * PZij
          VSuxi  = VSuxi+ 0.5 * FYij * PXij
          VSuyi  = VSuyi+ 0.5 * FZij * PXij
          VSuzi  = VSuzi+ 0.5 * FZij * PYij
          VSTempX(j) = VSTempX(j) + 0.5*FXij * PYij
          VSTempY(j) = VSTempY(j) + 0.5*FXij * PZij
          VSTempZ(j) = VSTempZ(j) + 0.5*FYij * PZij
          VBTempX(j) = VBTempX(j) + 0.5*FXij * PXij
          VBTempY(j) = VBTempY(j) + 0.5*FYij * PYij
          VBTempZ(j) = VBTempZ(j) + 0.5*FZij * PZij
          VSuTempX(j)= VSuTempX(j)+ 0.5*FYij * PXij
          VSuTempY(j)= VSuTempY(j)+ 0.5*FZij * PXij
          VSuTempZ(j)= VSuTempZ(j)+ 0.5*FZij * PYij

          UU     = 0.5 * EpotLocal1  
          Cxi    = Cxi  + UU
          Cyi    = Cyi  + UU
          Czi    = Czi  + UU
          CTempX(j) = CTempX(j) + UU
          CTempY(j) = CTempY(j) + UU
          CTempZ(j) = CTempZ(j) + UU

          FTXi   = - Epsilon1*CosTheta2*eX 
          FTYi   = - Epsilon1*CosTheta2*eY 
          FTZi   = - Epsilon1*CosTheta2*eZ 
          txii   = OYi * FTZi - OZi * FTYi
          tyii   = OZi * FTXi - OXi * FTZi
          tzii   = OXi * FTYi - OYi * FTXi
          txir   = A11 * txii + A12 * tyii + A13 * tzii
          tyir   = A21 * txii + A22 * tyii + A23 * tzii
          tzir   = A31 * txii + A32 * tyii + A33 * tzii
          tuxi   = tuxi + 0.5 * PXij *tyir
          tuyi   = tuyi + 0.5 * PXij *tzir
          tuzi   = tuzi + 0.5 * PYij *tzir
          tlxi   = tlxi + 0.5 * PYij *txir
          tlyi   = tlyi + 0.5 * PZij *txir
          tlzi   = tlzi + 0.5 * PZij *tyir
          tdxi   = tdxi + 0.5 * PXij *txir
          tdyi   = tdyi + 0.5 * PYij *tyir
          tdzi   = tdzi + 0.5 * PZij *tzir
          tuTempX(j)= tuTempX(j) + 0.5*PXij*tyi
          tuTempY(j)= tuTempY(j) + 0.5*PXij*tzi
          tuTempZ(j)= tuTempZ(j) + 0.5*PYij*tzi
          tlTempX(j)= tlTempX(j) + 0.5*PYij*txi
          tlTempY(j)= tlTempY(j) + 0.5*PZij*txi
          tlTempZ(j)= tlTempZ(j) + 0.5*PZij*tyi
          tdTempX(j)= tdTempX(j) + 0.5*PXij*txi
          tdTempY(j)= tdTempY(j) + 0.5*PYij*tyi
          tdTempZ(j)= tdTempZ(j) + 0.5*PZij*tzi          
          !TRANSPORT_END
#endif
      end do loop1

      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
      TX1(i) = TXi
      TY1(i) = TYi
      TZ1(i) = TZi
#if  TRANS == 1
        !TRANSPORT_start
        VSx1(i) = VSxi
        VSy1(i) = VSyi
        VSz1(i) = VSzi
        VBx1(i) = VBxi
        VBy1(i) = VByi
        VBz1(i) = VBzi
        VSux1(i)= VSuxi
        VSuy1(i)= VSuyi
        VSuz1(i)= VSuzi
        Cx1(i)  = Cxi
        Cy1(i)  = Cyi
        Cz1(i)  = Czi
        tux1(i) = tuxi
        tuy1(i) = tuyi
        tuz1(i) = tuzi
        tlx1(i) = tlxi
        tly1(i) = tlyi
        tlz1(i) = tlzi
        tdx1(i) = tdxi
        tdy1(i) = tdyi
        tdz1(i) = tdzi
        !TRANSPORT_END
#endif
    end do
!$OMP END DO
!$OMP END PARALLEL




    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ

#if  TRANS == 1
   VSx2 = VSx2 + VSTempX
   VSy2 = VSy2 + VSTempY
   VSz2 = VSz2 + VSTempZ

   VSux2 = VSux2 + VSuTempX
   VSuy2 = VSuy2 + VSuTempY
   VSuz2 = VSuz2 + VSuTempZ

   VBx2 = VBx2 + VBTempX
   VBy2 = VBy2 + VBTempY
   VBz2 = VBz2 + VBTempZ

   Cx2 = Cx2 + CTempX
   Cy2 = Cy2 + CTempY
   Cz2 = Cz2 + CTempZ

   tux2 = tux2 + tuTempX
   tuy2 = tuy2 + tuTempY
   tuz2 = tuz2 + tuTempZ

   tdx2 = tdx2 + tdTempX
   tdy2 = tdy2 + tdTempY
   tdz2 = tdz2 + tdTempZ

   tlx2 = tlx2 + tlTempX
   tly2 = tly2 + tlTempY
   tlz2 = tlz2 + tlTempZ

#endif

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotQC_Force_Trans



!==============================================================!
!  Subroutine TPotQC_ChemicalPotential                         !
!==============================================================!

  subroutine TPotQC_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleCharge) :: this
    real(RK), pointer, contiguous          :: EPotTest(:)
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv, RijSquared
    real(RK)          :: CosTheta
    real(RK)          :: EPotLocal
    integer           :: i, j, k, i1

    ! Assign local variables
    i1 = this%Site1%NTest
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

    ! Assign pointers
    RX1 => this%Site1%RXTest
    RY1 => this%Site1%RYTest
    RZ1 => this%Site1%RZTest
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    PX1 => this%Site1%PXTest
    PY1 => this%Site1%PYTest
    PZ1 => this%Site1%PZTest
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
    OX1 => this%Site1%OXTest
    OY1 => this%Site1%OYTest
    OZ1 => this%Site1%OZTest

    ! Loop over test particles
    do i = 1, i1
      RXi = RX1(i)
      RYi = RY1(i)
      RZi = RZ1(i)
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)
      OXi = OX1(i)
      OYi = OY1(i)
      OZi = OZ1(i)
      EPotLocal = 0._RK

!CDIR NODEP
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij

          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop1
          end if
          RijSquaredInv = 1._RK / RijSquared
          RijInv = sqrt( RijSquaredInv )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosTheta  = OXi * ex + OYi * eY + OZi * eZ
          EPotLocal = EPotLocal + Epsilon * RijSquaredInv * RijInv * ( CosTheta * CosTheta - Third )
        end do loop1

        EPotTest(i) = EPotTest(i) + EPotLocal
    end do

  end subroutine TPotQC_ChemicalPotential



!==============================================================!
!  Subroutine TPotQD_Construct                                 !
!==============================================================!

  subroutine TPotQD_Construct( this, i1, i2, j1, j2, Molecule1, Molecule2, RCutoff )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleDipole)  :: this
    integer, intent(in)         :: i1, i2, j1, j2
    type(TMolecule), intent(in) :: Molecule1, Molecule2
    real(RK), intent(in)        :: RCutoff

    ! Construct potential
    this%Site1 => Molecule1%SiteQuadrupole(j1)
    this%Site2 => Molecule2%SiteDipole(j2)
    this%SameComponent = i1 == i2
    this%Epsilon = 1.5_RK * this%Site1%Q * this%Site2%D
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

  end subroutine TPotQD_Construct



!==============================================================!
!  Subroutine TPotQD_Destruct                                  !
!==============================================================!

  subroutine TPotQD_Destruct( this )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleDipole) :: this

    ! Destroy potential
    continue

  end subroutine TPotQD_Destruct



!==============================================================!
!  Subroutine TPotQD_Force                                     !
!==============================================================!

  subroutine TPotQD_Force( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleDipole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: d2EpotdV2
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: TX1(:), TY1(:), TZ1(:), TX2(:), TY2(:), TZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: TXi, TYi, TZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij4Inv
    real(RK)          :: CosThetai, CosThetaj, CosThetai2, CosGammaij
    real(RK)          :: dCosThetai, dCosThetaj, dCosGammaij
    real(RK)          :: Tmp
    real(RK)          :: EPotLocal1, EPotLocal, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: momTempX(1:this%Site2%NPart)
    real(RK)          :: momTempY(1:this%Site2%NPart)
    real(RK)          :: momTempZ(1:this%Site2%NPart)
    
#if MPI_VER > 0
    integer           :: N1, N2, i0, ji
    logical           :: EvenN
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ    
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    momTempX(:)=0._RK
    momTempY(:)=0._RK
    momTempZ(:)=0._RK    
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

!$OMP PARALLEL &
!$OMP PRIVATE (Epsilon,  RCutoffSquared) &
!$OMP PRIVATE (RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (OX1, OY1, OZ1, OX2, OY2, OZ2) &
!$OMP PRIVATE (FX1, FY1, FZ1, TX1, TY1, TZ1) &
!$OMP PRIVATE (Plen2,sitecorr,PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (RXi, RYi, RZi, OXi, OYi, OZi,  FXi, FYi, FZi) &
!$OMP PRIVATE (TXi, TYi, TZi, PXi, PYi, PZi,  RXij, RYij, RZij) &
!$OMP PRIVATE (OXj, OYj, OZj, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE (eX, eY, eZ, RijSquared, RijInv, Rij4Inv) &
!$OMP PRIVATE (CosThetai, CosThetaj, CosThetai2, CosGammaij) &
!$OMP PRIVATE (dCosThetai, dCosThetaj, dCosGammaij) &
!$OMP PRIVATE (Tmp, EPotLocal1) &
!$OMP PRIVATE (SameComponent) &
#if MPI_VER > 0
!$OMP PRIVATE (i, j, k, i1, j0, j1) &
!$OMP PRIVATE ( N1, N2, i0, ji, EvenN)
#else
!$OMP PRIVATE (i, j, k, i1, j0, j1)
#endif

    ! Assign local variables
    SameComponent = this%SameComponent
#if MPI_VER > 0
    N1 = this%Site2%NPart
    N2 = N1 / 2
    EvenN = mod( N1, 2 ) == 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
    j1 = this%Site2%NPart
#endif
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)

!CDIR NODEP
#if OSMOP == 2
loop0:  do m=1,NBinsDen
          if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PXi < real(m)/NBinsDen-0.5_RK) then
              Bin1=m
              exit loop0
            end if
          end if
        end do loop0
#endif
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)
          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosThetai2 = CosThetai**2
          CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal1 = Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) - CosGammaij * CosThetai)
          EPotLocal = EPotLocal + EPotLocal1
          dCosThetai = Rij4Inv * (10._RK * CosThetai * CosThetaj - CosGammaij)
          dCosThetaj = Rij4Inv * (5._RK * CosThetai2 - 1._RK)
          dCosGammaij = -2._RK * Rij4Inv * CosThetai
          Tmp = -4._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)

          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
                Bin2=m 
                exit loop2
              end if
            end if
          end do loop2
          tempMin = min(Bin1, Bin2)
          tempMax = max(Bin1, Bin2)
          if(abs(PXij) .le. 0.5_RK) then
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
              do m = tempMin, tempMax
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          else
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
              do m = 1, tempMin
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
              do m = tempMax, NBinsDen
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          end if
#endif
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Ninth    !xxxx8   QD

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij
          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij

          momTempX(j) = momTempX(j) - eX * dCosThetaj - OXi * dCosGammaij
          momTempY(j) = momTempY(j) - eY * dCosThetaj - OYi * dCosGammaij  
          momTempZ(j) = momTempZ(j) - eZ * dCosThetaj - OZi * dCosGammaij   

        end do loop1
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, merge( i1 - 1, i1, SameComponent )
#endif
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)

#if MPI_VER > 0
        if( SameComponent ) then
          j0 = i + 1
          j1 = i + N2
          if( EvenN .and. i > N2 ) j1 = j1 - 1
        else
          j0 = 1
          j1 = N1
        end if

!CDIR NODEP
loop3:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )
#else
        j0 = merge( i + 1, 1, SameComponent )

!CDIR NODEP
loop3:  do j = j0, j1
#endif
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = (PXij - anint( RXij )) * BoxLength
          PYij = (PYij - anint( RYij )) * BoxLength
          PZij = (PZij - anint( RZij )) * BoxLength
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosThetai2 = CosThetai**2
          CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal1 = Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) - CosGammaij * CosThetai)
          EPotLocal = EPotLocal + EPotLocal1
          dCosThetai = Rij4Inv * (10._RK * CosThetai * CosThetaj - CosGammaij)
          dCosThetaj = Rij4Inv * (5._RK * CosThetai2 - 1._RK)
          dCosGammaij = -2._RK * Rij4Inv * CosThetai
          Tmp = -4._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)

          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)

          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Ninth    !xxxx8   QD ss

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij
          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
          momTempX(j) = momTempX(j) - eX * dCosThetaj - OXi * dCosGammaij
          momTempY(j) = momTempY(j) - eY * dCosThetaj - OYi * dCosGammaij  
          momTempZ(j) = momTempZ(j) - eZ * dCosThetaj - OZi * dCosGammaij   

        end do loop3
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do
!$OMP END DO

    end if
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ
    TX2 = TX2 + momTempX                                 
    TY2 = TY2 + momTempY
    TZ2 = TZ2 + momTempZ

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotQD_Force

!==============================================================!
!  Subroutine TPotQD_Force_Trans                               !
!==============================================================!

  subroutine TPotQD_Force_Trans( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleDipole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: d2EpotdV2
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: TX1(:), TY1(:), TZ1(:), TX2(:), TY2(:), TZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: TXi, TYi, TZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij4Inv
    real(RK)          :: CosThetai, CosThetaj, CosThetai2, CosGammaij
    real(RK)          :: dCosThetai, dCosThetaj, dCosGammaij
    real(RK)          :: Tmp
    real(RK)          :: EPotLocal1, EPotLocal, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: momTempX(1:this%Site2%NPart)
    real(RK)          :: momTempY(1:this%Site2%NPart)
    real(RK)          :: momTempZ(1:this%Site2%NPart)
    
#if MPI_VER > 0
    integer           :: N1, N2, i0, ji
    logical           :: EvenN
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

#if  TRANS == 1
    real(RK), pointer, contiguous :: VSx1(:), VSy1(:), VSz1(:)
    real(RK), pointer, contiguous :: VSux1(:),VSuy1(:),VSuz1(:)
    real(RK), pointer, contiguous :: VBx1(:), VBy1(:), VBz1(:)
    real(RK), pointer, contiguous :: Cx1(:) , Cy1(:) , Cz1(:)
    real(RK), pointer, contiguous :: tux1(:) , tuy1(:) , tuz1(:)
    real(RK), pointer, contiguous :: tlx1(:) , tly1(:) , tlz1(:)
    real(RK), pointer, contiguous :: tdx1(:) , tdy1(:) , tdz1(:)
    real(RK), pointer, contiguous :: VSx2(:), VSy2(:), VSz2(:)
    real(RK), pointer, contiguous :: VSux2(:),VSuy2(:),VSuz2(:)
    real(RK), pointer, contiguous :: VBx2(:), VBy2(:), VBz2(:)
    real(RK), pointer, contiguous :: Cx2(:) , Cy2(:) , Cz2(:)
    real(RK), pointer, contiguous :: tux2(:) , tuy2(:) , tuz2(:)
    real(RK), pointer, contiguous :: tlx2(:) , tly2(:) , tlz2(:)
    real(RK), pointer, contiguous :: tdx2(:) , tdy2(:) , tdz2(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
    real(RK)          :: VSTempX(1:this%Site2%NPart), VSTempY(1:this%Site2%NPart), VSTempZ(1:this%Site2%NPart)
    real(RK)          :: VSuTempX(1:this%Site2%NPart), VSuTempY(1:this%Site2%NPart), VSuTempZ(1:this%Site2%NPart)
    real(RK)          :: VBTempX(1:this%Site2%NPart), VBTempY(1:this%Site2%NPart), VBTempZ(1:this%Site2%NPart)
    real(RK)          :: CTempX(1:this%Site2%NPart), CTempY(1:this%Site2%NPart), CTempZ(1:this%Site2%NPart)
    real(RK)          :: tuTempX(1:this%Site2%NPart), tuTempY(1:this%Site2%NPart), tuTempZ(1:this%Site2%NPart)
    real(RK)          :: tlTempX(1:this%Site2%NPart), tlTempY(1:this%Site2%NPart), tlTempZ(1:this%Site2%NPart)
    real(RK)          :: tdTempX(1:this%Site2%NPart), tdTempY(1:this%Site2%NPart), tdTempZ(1:this%Site2%NPart)
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VBxi, VByi, VBzi
    real(RK)          :: Cxi,  Cyi,  Czi
    real(RK)          :: tuxi,  tuyi,  tuzi
    real(RK)          :: tlxi,  tlyi,  tlzi
    real(RK)          :: tdxi,  tdyi,  tdzi
    real(RK)          :: txii,  tyii , tzii
    real(RK)          :: txir , tyir , tzir
    real(RK)          :: FTXi , FTYi , FTZi
    real(RK)          :: UU
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33

    VSTempX(:) = 0._RK
    VSTempY(:) = 0._RK
    VSTempZ(:) = 0._RK
    VSuTempX(:)= 0._RK
    VSuTempY(:)= 0._RK
    VSuTempZ(:)= 0._RK
    VBTempX(:) = 0._RK
    VBTempY(:) = 0._RK
    VBTempZ(:) = 0._RK
    CTempX(:)  = 0._RK
    CTempY(:)  = 0._RK
    CTempZ(:)  = 0._RK
    tuTempX(:) = 0._RK
    tuTempY(:) = 0._RK
    tuTempz(:) = 0._RK
    tlTempX(:) = 0._RK
    tlTempY(:) = 0._RK
    tlTempz(:) = 0._RK
    tdTempX(:) = 0._RK
    tdTempY(:) = 0._RK
    tdTempz(:) = 0._RK
#endif

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ    

    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    momTempX(:)=0._RK
    momTempY(:)=0._RK
    momTempZ(:)=0._RK    
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

!$OMP PARALLEL &
!$OMP PRIVATE (Epsilon,  RCutoffSquared) &
!$OMP PRIVATE (RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (OX1, OY1, OZ1, OX2, OY2, OZ2) &
!$OMP PRIVATE (FX1, FY1, FZ1, TX1, TY1, TZ1) &
!$OMP PRIVATE (Plen2,sitecorr,PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (RXi, RYi, RZi, OXi, OYi, OZi,  FXi, FYi, FZi) &
!$OMP PRIVATE (TXi, TYi, TZi, PXi, PYi, PZi,  RXij, RYij, RZij) &
!$OMP PRIVATE (OXj, OYj, OZj, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE (eX, eY, eZ, RijSquared, RijInv, Rij4Inv) &
!$OMP PRIVATE (CosThetai, CosThetaj, CosThetai2, CosGammaij) &
!$OMP PRIVATE (dCosThetai, dCosThetaj, dCosGammaij) &
!$OMP PRIVATE (Tmp, EPotLocal1) &
!$OMP PRIVATE (SameComponent) &
#if  TRANS == 1
!$OMP PRIVATE(VSx1, VSy1, VSz1, VSux1,VSuy1,VSuz1, VBx1, VBy1, VBz1, Cx1, Cy1, Cz1) &
!$OMP PRIVATE(VSx2, VSy2, VSz2, VSux2,VSuy2,VSuz2, VBx2, VBy2, VBz2, Cx2, Cy2, Cz2) &
!$OMP PRIVATE( tux1, tuy1, tuz1, tlx1, tly1, tlz1, tdx1, tdy1, tdz1) &
!$OMP PRIVATE( tux2, tuy2, tuz2, tlx2, tly2, tlz2, tdx2, tdy2, tdz2) &
!$OMP PRIVATE( q1, q2, q3, q4, VSxi, VSyi, VSzi, VSuxi,VSuyi,VSuzi) &
!$OMP PRIVATE( VBxi, VByi, VBzi, Cxi,  Cyi,  Czi, tuxi,  tuyi,  tuzi, tlxi,  tlyi,  tlzi) &
!$OMP PRIVATE(  tdxi,  tdyi,  tdzi, txii,  tyii , tzii, txir ,  tyir  , tzir ) &
!$OMP PRIVATE( UU, FTXi , FTYi , FTZi) &
!$OMP PRIVATE( A11, A12, A13, A21, A22, A23, A31, A32, A33) &
#endif
#if MPI_VER > 0
!$OMP PRIVATE (i, j, k, i1, j0, j1) &
!$OMP PRIVATE ( N1, N2, i0, ji, EvenN)
#else
!$OMP PRIVATE (i, j, k, i1, j0, j1)
#endif


    ! Assign local variables
    SameComponent = this%SameComponent
#if MPI_VER > 0
    N1 = this%Site2%NPart
    N2 = N1 / 2
    EvenN = mod( N1, 2 ) == 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
    j1 = this%Site2%NPart
#endif
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
#if  TRANS == 1
 
    VSx1 => this%Site1%vsQx
    VSy1 => this%Site1%vsQy
    VSz1 => this%Site1%vsQz
    VSx2 => this%Site2%vsDx
    VSy2 => this%Site2%vsDy
    VSz2 => this%Site2%vsDz
    VBx1 => this%Site1%vbQx
    VBy1 => this%Site1%vbQy
    VBz1 => this%Site1%vbQz
    VBx2 => this%Site2%vbDx
    VBy2 => this%Site2%vbDy
    VBz2 => this%Site2%vbDz
    VSux1=> this%Site1%vsuQx
    VSuy1=> this%Site1%vsuQy
    VSuz1=> this%Site1%vsuQz
    VSux2=> this%Site2%vsuDx
    VSuy2=> this%Site2%vsuDy
    VSuz2=> this%Site2%vsuDz
    Cx1  => this%Site1%cQx
    Cy1  => this%Site1%cQy
    Cz1  => this%Site1%cQz
    Cx2  => this%Site2%cDx
    Cy2  => this%Site2%cDy
    Cz2  => this%Site2%cDz
    tux1 => this%Site1%tuQx
    tuy1 => this%Site1%tuQy
    tuz1 => this%Site1%tuQz
    tux2 => this%Site2%tuDx
    tuy2 => this%Site2%tuDy
    tuz2 => this%Site2%tuDz
    tlx1 => this%Site1%tlQx
    tly1 => this%Site1%tlQy
    tlz1 => this%Site1%tlQz
    tlx2 => this%Site2%tlDx
    tly2 => this%Site2%tlDy
    tlz2 => this%Site2%tlDz
    tdx1 => this%Site1%tdQx
    tdy1 => this%Site1%tdQy
    tdz1 => this%Site1%tdQz
    tdx2 => this%Site2%tdDx
    tdy2 => this%Site2%tdDy
    tdz2 => this%Site2%tdDz
    q1  => this%Site1%Q0r(:, 1)
    q2  => this%Site1%Q0r(:, 2)
    q3  => this%Site1%Q0r(:, 3)
    q4  => this%Site1%Q0r(:, 4)
!TRANSPORT_END
#endif
    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)

!CDIR NODEP
#if  TRANS == 1
        !TRANSPORT_start
        VSxi= VSx1(i)
        VSyi= VSy1(i)
        VSzi= VSz1(i)
        VBxi= VBx1(i)
        VByi= VBy1(i)
        VBzi= VBz1(i)
        VSuxi= VSux1(i)
        VSuyi= VSuy1(i)
        VSuzi= VSuz1(i)
        Cxi = Cx1(i)
        Cyi = Cy1(i)
        Czi = Cz1(i)
        tuxi = tux1(i)
        tuyi = tuy1(i)
        tuzi = tuz1(i)
        tlxi = tlx1(i)
        tlyi = tly1(i)
        tlzi = tlz1(i)
        tdxi = tdx1(i)
        tdyi = tdy1(i)
        tdzi = tdz1(i)
        A11 = q1(i)**2 + q2(i)**2 - q3(i)**2 - q4(i)**2
        A12 = 2._RK * (q2(i) * q3(i) + q1(i) * q4(i))
        A13 = 2._RK * (q2(i) * q4(i) - q1(i) * q3(i))
        A21 = 2._RK * (q2(i) * q3(i) - q1(i) * q4(i))
        A22 = q1(i)**2 - q2(i)**2 + q3(i)**2 - q4(i)**2
        A23 = 2._RK * (q3(i) * q4(i) + q1(i) * q2(i))
        A31 = 2._RK * (q2(i) * q4(i) + q1(i) * q3(i))
        A32 = 2._RK * (q3(i) * q4(i) - q1(i) * q2(i))
        A33 = q1(i)**2 - q2(i)**2 - q3(i)**2 + q4(i)**2
        !TRANSPORT_END
#endif

#if OSMOP == 2
loop0:  do m=1,NBinsDen
          if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PXi < real(m)/NBinsDen-0.5_RK) then
              Bin1=m
              exit loop0
            end if
          end if
        end do loop0
#endif
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosThetai2 = CosThetai**2
          CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal1 = Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) - CosGammaij * CosThetai)
          EPotLocal = EPotLocal + EPotLocal1
          dCosThetai = Rij4Inv * (10._RK * CosThetai * CosThetaj - CosGammaij)
          dCosThetaj = Rij4Inv * (5._RK * CosThetai2 - 1._RK)
          dCosGammaij = -2._RK * Rij4Inv * CosThetai
          Tmp = -4._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)

          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
                Bin2=m 
                exit loop2
              end if
            end if
          end do loop2
          tempMin = min(Bin1, Bin2)
          tempMax = max(Bin1, Bin2)
          if(abs(PXij) .le. 0.5_RK) then
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
              do m = tempMin, tempMax
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          else
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
              do m = 1, tempMin
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
              do m = tempMax, NBinsDen
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          end if
#endif
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Ninth    !xxxx8   QD T

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij
          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij

          momTempX(j) = momTempX(j) - eX * dCosThetaj - OXi * dCosGammaij
          momTempY(j) = momTempY(j) - eY * dCosThetaj - OYi * dCosGammaij  
          momTempZ(j) = momTempZ(j) - eZ * dCosThetaj - OZi * dCosGammaij   

#if  TRANS == 1
!TRANSPORT_start
          VSxi   = VSxi + 0.5 * FXij * PYij
          VSyi   = VSyi + 0.5 * FXij * PZij
          VSzi   = VSzi + 0.5 * FYij * PZij
          VBxi   = VBxi + 0.5 * FXij * PXij
          VByi   = VByi + 0.5 * FYij * PYij
          VBzi   = VBzi + 0.5 * FZij * PZij
          VSuxi  = VSuxi+ 0.5 * FYij * PXij
          VSuyi  = VSuyi+ 0.5 * FZij * PXij
          VSuzi  = VSuzi+ 0.5 * FZij * PYij
          VSTempX(j) = VSTempX(j) + 0.5*FXij * PYij
          VSTempY(j) = VSTempY(j) + 0.5*FXij * PZij
          VSTempZ(j) = VSTempZ(j) + 0.5*FYij * PZij
          VBTempX(j) = VBTempX(j) + 0.5*FXij * PXij
          VBTempY(j) = VBTempY(j) + 0.5*FYij * PYij
          VBTempZ(j) = VBTempZ(j) + 0.5*FZij * PZij
          VSuTempX(j)= VSuTempX(j)+ 0.5*FYij * PXij
          VSuTempY(j)= VSuTempY(j)+ 0.5*FZij * PXij
          VSuTempZ(j)= VSuTempZ(j)+ 0.5*FZij * PYij
          
          UU     = 0.5 * EpotLocal1  
          Cxi    = Cxi  + UU
          Cyi    = Cyi  + UU
          Czi    = Czi  + UU
          CTempX(j) = CTempX(j) + UU
          CTempY(j) = CTempY(j) + UU
          CTempZ(j) = CTempZ(j) + UU
          

          FTXi   = - eX * dCosThetai - OXj * dCosGammaij 
          FTYi   = - eY * dCosThetai - OYj * dCosGammaij
          FTZi   = - eZ * dCosThetai - OZj * dCosGammaij
          txii   = OYi * FTZi - OZi * FTYi
          tyii   = OZi * FTXi - OXi * FTZi
          tzii   = OXi * FTYi - OYi * FTXi
          txir   = A11 * txii + A12 * tyii + A13 * tzii
          tyir   = A21 * txii + A22 * tyii + A23 * tzii
          tzir   = A31 * txii + A32 * tyii + A33 * tzii
          tuxi   = tuxi + PXij*tyir
          tuyi   = tuyi + PXij*tzir
          tuzi   = tuzi + PYij*tzir
          tlxi   = tlxi + PYij*txir
          tlyi   = tlyi + PZij*txir
          tlzi   = tlzi + PZij*tyir
          tdxi   = tdxi + PXij*txir
          tdyi   = tdyi + PYij*tyir
          tdzi   = tdzi + PZij*tzir
          tuTempX(j)= tuTempX(j) + 0.5 * PXij*tyi
          tuTempY(j)= tuTempY(j) + 0.5 * PXij*tzi
          tuTempZ(j)= tuTempZ(j) + 0.5 * PYij*tzi
          tlTempX(j)= tlTempX(j) + 0.5 * PYij*txi
          tlTempY(j)= tlTempY(j) + 0.5 * PZij*txi
          tlTempZ(j)= tlTempZ(j) + 0.5 * PZij*tyi
          tdTempX(j)= tdTempX(j) + 0.5 * PXij*txi
          tdTempY(j)= tdTempY(j) + 0.5 * PYij*tyi
          tdTempZ(j)= tdTempZ(j) + 0.5 * PZij*tzi
  
          !TRANSPORT_END
#endif
        end do loop1

        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi

#if  TRANS == 1
        !TRANSPORT_start
        VSx1(i) = VSxi
        VSy1(i) = VSyi
        VSz1(i) = VSzi
        VBx1(i) = VBxi
        VBy1(i) = VByi
        VBz1(i) = VBzi
        VSux1(i)= VSuxi
        VSuy1(i)= VSuyi
        VSuz1(i)= VSuzi
        Cx1(i)  = Cxi
        Cy1(i)  = Cyi
        Cz1(i)  = Czi
        tux1(i) = tuxi
        tuy1(i) = tuyi
        tuz1(i) = tuzi
        tlx1(i) = tlxi
        tly1(i) = tlyi
        tlz1(i) = tlzi
        tdx1(i) = tdxi
        tdy1(i) = tdyi
        tdz1(i) = tdzi
        !TRANSPORT_END
#endif
      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, merge( i1 - 1, i1, SameComponent )
#endif
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)

#if MPI_VER > 0
        if( SameComponent ) then
          j0 = i + 1
          j1 = i + N2
          if( EvenN .and. i > N2 ) j1 = j1 - 1
        else
          j0 = 1
          j1 = N1
        end if

!CDIR NODEP
loop3:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )
#else
        j0 = merge( i + 1, 1, SameComponent )
!CDIR NODEP
loop3:  do j = j0, j1
#endif
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = (PXij - anint( RXij )) * BoxLength
          PYij = (PYij - anint( RYij )) * BoxLength
          PZij = (PZij - anint( RZij )) * BoxLength
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosThetai2 = CosThetai**2
          CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal1 = Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) - CosGammaij * CosThetai)
          EPotLocal = EPotLocal + EPotLocal1
          dCosThetai = Rij4Inv * (10._RK * CosThetai * CosThetaj - CosGammaij)
          dCosThetaj = Rij4Inv * (5._RK * CosThetai2 - 1._RK)
          dCosGammaij = -2._RK * Rij4Inv * CosThetai
          Tmp = -4._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)

          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)

          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Ninth    !xxxx8   QD ss T

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij
          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij

          momTempX(j) = momTempX(j) - eX * dCosThetaj - OXi * dCosGammaij
          momTempY(j) = momTempY(j) - eY * dCosThetaj - OYi * dCosGammaij  
          momTempZ(j) = momTempZ(j) - eZ * dCosThetaj - OZi * dCosGammaij   

        end do loop3
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do
!$OMP END DO

    end if
!$OMP END PARALLEL


    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ
    TX2 = TX2 + momTempX                                 
    TY2 = TY2 + momTempY
    TZ2 = TZ2 + momTempZ
#if  TRANS == 1
   VSx2 = VSx2 + VSTempX
   VSy2 = VSy2 + VSTempY
   VSz2 = VSz2 + VSTempZ

   VSux2 = VSux2 + VSuTempX
   VSuy2 = VSuy2 + VSuTempY
   VSuz2 = VSuz2 + VSuTempZ

   VBx2 = VBx2 + VBTempX
   VBy2 = VBy2 + VBTempY
   VBz2 = VBz2 + VBTempZ

   Cx2 = Cx2 + CTempX
   Cy2 = Cy2 + CTempY
   Cz2 = Cz2 + CTempZ

   tux2 = tux2 + tuTempX
   tuy2 = tuy2 + tuTempY
   tuz2 = tuz2 + tuTempZ

   tdx2 = tdx2 + tdTempX
   tdy2 = tdy2 + tdTempY
   tdz2 = tdz2 + tdTempZ

   tlx2 = tlx2 + tlTempX
   tly2 = tly2 + tlTempY
   tlz2 = tlz2 + tlTempZ

#endif



    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotQD_Force_Trans



!==============================================================!
!  Subroutine TPotQD_ChemicalPotential                         !
!==============================================================!

  subroutine TPotQD_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleDipole) :: this
    real(RK), pointer, contiguous          :: EPotTest(:)
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij4Inv
    real(RK)          :: CosThetai, CosThetaj, CosGammaij
    real(RK)          :: EPotLocal
    integer           :: i, j, k, i1, j1

    ! Assign local variables
    i1 = this%Site1%NTest
    j1 = this%Site2%NPart
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

    ! Assign pointers
    RX1 => this%Site1%RXTest
    RY1 => this%Site1%RYTest
    RZ1 => this%Site1%RZTest
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    OX1 => this%Site1%OXTest
    OY1 => this%Site1%OYTest
    OZ1 => this%Site1%OZTest
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    PX1 => this%Site1%PXTest
    PY1 => this%Site1%PYTest
    PZ1 => this%Site1%PZTest
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (RXi,RYi,RZi,PXi,PYi,PZi) &
!$OMP PRIVATE (OXi, OYi, OZi, OXj, OYj, OZj) &
!$OMP PRIVATE (RXij,RYij,RZij,PXij,PYij,PZij) &
!$OMP PRIVATE (CosThetai, CosThetaj, CosGammaij) &
!$OMP PRIVATE (RijSquared,RijInv, Rij4Inv) &
!$OMP PRIVATE (eX,eY,eZ) &
!$OMP PRIVATE (EPotLocal,i,j,k) 

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
!$OMP DO
      do i = 1, i1
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        EPotLocal = 0._RK


!CDIR NODEP
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij

          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop1
          end if
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal = EPotLocal + Rij4Inv * (CosThetaj * (5._RK * CosThetai**2 - 1._RK) - 2._RK * CosGammaij * CosThetai)
        end do loop1

        EPotTest(i) = EPotTest(i) + EPotLocal
      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over test particles
!$OMP DO
      do i = 1, i1
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        EPotLocal = 0._RK

!CDIR NODEP
loop2:  do j = 1, j1
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          if( RijSquared >= RCutoffSquared ) cycle loop2
          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop2
          end if
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)
          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal = EPotLocal + Rij4Inv * ( CosThetaj * (5._RK * CosThetai**2 - 1._RK) - 2._RK * CosGammaij * CosThetai )
        end do loop2

        EPotTest(i) = EPotTest(i) + EPotLocal
      end do
!$OMP END DO
    end if
!$OMP END PARALLEL
  end subroutine TPotQD_ChemicalPotential



!==============================================================!
!  Subroutine TPotQQ_Construct                                 !
!==============================================================!

  subroutine TPotQQ_Construct( this, i1, i2, j1, j2, Molecule1, Molecule2, RCutoff )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleQuadrupole) :: this
    integer, intent(in)            :: i1, i2, j1, j2
    type(TMolecule), intent(in)    :: Molecule1, Molecule2
    real(RK), intent(in)           :: RCutoff

    ! Construct potential
    this%Site1 => Molecule1%SiteQuadrupole(j1)
    this%Site2 => Molecule2%SiteQuadrupole(j2)
    this%SameComponent = i1 == i2
    this%Epsilon = .75_RK * this%Site1%Q * this%Site2%Q
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

  end subroutine TPotQQ_Construct



!==============================================================!
!  Subroutine TPotQQ_Destruct                                  !
!==============================================================!

  subroutine TPotQQ_Destruct( this )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleQuadrupole) :: this

    ! Destroy potential
    continue

  end subroutine TPotQQ_Destruct



!==============================================================!
!  Subroutine TPotQQ_Force                                     !
!==============================================================!

  subroutine TPotQQ_Force( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleQuadrupole) :: this
    real(RK), intent(in out)       :: EPot
    real(RK), intent(in out)       :: Virial
    real(RK), intent(in out)       :: d2EpotdV2
    real(RK), intent(in)           :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: TX1(:), TY1(:), TZ1(:), TX2(:), TY2(:), TZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: TXi, TYi, TZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij5Inv
    real(RK)          :: CosThetai, CosThetaj, CosGammaij
    real(RK)          :: CosThetaiSquared, CosThetajSquared
    real(RK)          :: dCosThetai, dCosThetaj, dCosGammaij
    real(RK)          :: Tmp
    real(RK)          :: EPotLocal1, EPotLocal, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: momTempX(1:this%Site2%NPart)
    real(RK)          :: momTempY(1:this%Site2%NPart)
    real(RK)          :: momTempZ(1:this%Site2%NPart)
    
#if MPI_VER > 0
    integer           :: N1, N2, i0, ji
    logical           :: EvenN
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ
    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ    

    ! Assign local variables
    SameComponent = this%SameComponent
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    momTempX(:)=0._RK
    momTempY(:)=0._RK
    momTempZ(:)=0._RK    
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK
    i = 0
    j = 0
    k = 0
    j0= 0
#if MPI_VER > 0
    N1 = this%Site2%NPart
    N2 = N1 / 2
    EvenN = mod( N1, 2 ) == 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
    ji = 0
    j1 = 0
#else
    i1 = this%Site1%NPart
    j1 = this%Site2%NPart
#endif

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE (i, j, k, i1, j0, j1) &
#if MPI_VER > 0
!$OMP FIRSTPRIVATE ( N1, N2, i0, ji, EvenN) &
#endif
!$OMP PRIVATE (RXi, RYi, RZi, OXi, OYi, OZi,  FXi, FYi, FZi) &
!$OMP PRIVATE (TXi, TYi, TZi, PXi, PYi, PZi,  RXij, RYij, RZij) &
!$OMP PRIVATE (OXj, OYj, OZj, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE (eX, eY, eZ, RijSquared, RijInv, Rij5Inv) &
!$OMP PRIVATE (CosThetai, CosThetaj, CosGammaij) &
!$OMP PRIVATE (CosThetaiSquared, CosThetajSquared) &
!$OMP PRIVATE (dCosThetai, dCosThetaj, dCosGammaij, Tmp) &
!$OMP PRIVATE (EPotLocal1, Plen2, sitecorr) 

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)

!CDIR NODEP
#if OSMOP == 2
loop0:  do m=1,NBinsDen
          if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PXi < real(m)/NBinsDen-0.5_RK) then
              Bin1=m
              exit loop0
            end if
          end if
        end do loop0
#endif
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          CosThetaiSquared = CosThetai**2
          CosThetajSquared = CosThetaj**2
          Tmp = CosGammaij - 5._RK * CosThetai * CosThetaj

#if ARCH == 1
          Rij5Inv = Epsilon * RijInv * (RijInv**2)**2
#else
          Rij5Inv = Epsilon * RijInv**5
#endif
          EPotLocal1 = Rij5Inv * (1._RK - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&                      - 15._RK * CosThetaiSquared * CosThetajSquared + 2._RK * Tmp**2)

          EPotLocal = EPotLocal + EPotLocal1

          dCosThetai = Rij5Inv * (-10._RK * CosThetai - 30._RK * CosThetai * CosThetajSquared &
&                                - 20._RK * CosThetaj * Tmp)
          dCosThetaj = Rij5Inv * (-10._RK * CosThetaj - 30._RK * CosThetaj * CosThetaiSquared &
&                                - 20._RK * CosThetai * Tmp)

          dCosGammaij = 4._RK * Rij5Inv * Tmp
          Tmp = -5._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)

          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
                Bin2=m 
                exit loop2
              end if
            end if
          end do loop2
          tempMin = min(Bin1, Bin2)
          tempMax = max(Bin1, Bin2)
          if(abs(PXij) .le. 0.5_RK) then
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
              do m = tempMin, tempMax
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          else
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
              do m = 1, tempMin
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
              do m = tempMax, NBinsDen
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          end if
#endif
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(35._RK*sitecorr*sitecorr-5._RK*Plen2*RijInv*RijInv)*Ninth      !xxxx9  QQ

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij

          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
          momTempX(j) = momTempX(j) - eX * dCosThetaj - OXi * dCosGammaij
          momTempY(j) = momTempY(j) - eY * dCosThetaj - OYi * dCosGammaij  
          momTempZ(j) = momTempZ(j) - eZ * dCosThetaj - OZi * dCosGammaij   

        end do loop1
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do
!$OMP END DO

    else ! Site-site cutoff

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, merge( i1 - 1, i1, SameComponent )
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)

#if MPI_VER > 0
        if( SameComponent ) then
          j0 = i + 1
          j1 = i + N2
          if( EvenN .and. i > N2 ) j1 = j1 - 1
        else
          j0 = 1
          j1 = N1
        end if

!CDIR NODEP
loop3:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )
#else
        j0 = merge( i + 1, 1, SameComponent )
!CDIR NODEP
loop3:  do j = j0, j1
#endif

          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = (PXij - anint( RXij )) * BoxLength
          PYij = (PYij - anint( RYij )) * BoxLength
          PZij = (PZij - anint( RZij )) * BoxLength
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          CosThetaiSquared = CosThetai**2
          CosThetajSquared = CosThetaj**2
          Tmp = CosGammaij - 5._RK * CosThetai * CosThetaj

#if ARCH == 1
          Rij5Inv = Epsilon * RijInv * (RijInv**2)**2
#else
          Rij5Inv = Epsilon * RijInv**5
#endif

          EPotLocal1 = Rij5Inv * (1._RK - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&                      - 15._RK * CosThetaiSquared * CosThetajSquared + 2._RK * Tmp**2)

          EPotLocal = EPotLocal + EPotLocal1

          dCosThetai = Rij5Inv * (-10._RK * CosThetai - 30._RK * CosThetai * CosThetajSquared &
&                                - 20._RK * CosThetaj * Tmp)
          dCosThetaj = Rij5Inv * (-10._RK * CosThetaj - 30._RK * CosThetaj * CosThetaiSquared &
&                                - 20._RK * CosThetai * Tmp)

          dCosGammaij = 4._RK * Rij5Inv * Tmp
          Tmp = -5._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)

          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)

          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(35._RK*sitecorr*sitecorr-5._RK*Plen2*RijInv*RijInv)*Ninth     !xxxx9  QQ ss

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij

          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij

          momTempX(j) = momTempX(j) - eX * dCosThetaj - OXi * dCosGammaij
          momTempY(j) = momTempY(j) - eY * dCosThetaj - OYi * dCosGammaij  
          momTempZ(j) = momTempZ(j) - eZ * dCosThetaj - OZi * dCosGammaij   

        end do loop3

        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do
!$OMP END DO

    end if
!$OMP END PARALLEL

    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ
    TX2 = TX2 + momTempX                                 
    TY2 = TY2 + momTempY
    TZ2 = TZ2 + momTempZ


    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotQQ_Force


!==============================================================!
!  Subroutine TPotQQ_Force_Trans                               !
!==============================================================!

  subroutine TPotQQ_Force_Trans( this, EPot, Virial, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleQuadrupole) :: this
    real(RK), intent(in out)       :: EPot
    real(RK), intent(in out)       :: Virial
    real(RK), intent(in out)       :: d2EpotdV2
    real(RK), intent(in)           :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer, contiguous :: TX1(:), TY1(:), TZ1(:), TX2(:), TY2(:), TZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: TXi, TYi, TZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij5Inv
    real(RK)          :: CosThetai, CosThetaj, CosGammaij
    real(RK)          :: CosThetaiSquared, CosThetajSquared
    real(RK)          :: dCosThetai, dCosThetaj, dCosGammaij
    real(RK)          :: Tmp
    real(RK)          :: EPotLocal1, EPotLocal, VirialLocal
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1

    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: momTempX(1:this%Site2%NPart)
    real(RK)          :: momTempY(1:this%Site2%NPart)
    real(RK)          :: momTempZ(1:this%Site2%NPart)
    
#if MPI_VER > 0
    integer           :: N1, N2, i0, ji
    logical           :: EvenN
#endif
#if OSMOP == 2
    integer           :: m
    real(RK)          :: VirialPart
    integer           :: Bin1, Bin2
    integer           :: tempMin, tempMax
#endif

#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer, contiguous :: VSx1(:), VSy1(:), VSz1(:)
    real(RK), pointer, contiguous :: VSux1(:),VSuy1(:),VSuz1(:)
    real(RK), pointer, contiguous :: VBx1(:), VBy1(:), VBz1(:)
    real(RK), pointer, contiguous :: Cx1(:) , Cy1(:) , Cz1(:)
    real(RK), pointer, contiguous :: tux1(:) , tuy1(:) , tuz1(:)
    real(RK), pointer, contiguous :: tlx1(:) , tly1(:) , tlz1(:)
    real(RK), pointer, contiguous :: tdx1(:) , tdy1(:) , tdz1(:)
    real(RK), pointer, contiguous :: VSx2(:), VSy2(:), VSz2(:)
    real(RK), pointer, contiguous :: VSux2(:),VSuy2(:),VSuz2(:)
    real(RK), pointer, contiguous :: VBx2(:), VBy2(:), VBz2(:)
    real(RK), pointer, contiguous :: Cx2(:) , Cy2(:) , Cz2(:)
    real(RK), pointer, contiguous :: tux2(:) , tuy2(:) , tuz2(:)
    real(RK), pointer, contiguous :: tlx2(:) , tly2(:) , tlz2(:)
    real(RK), pointer, contiguous :: tdx2(:) , tdy2(:) , tdz2(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
    real(RK)          :: VSTempX(1:this%Site2%NPart), VSTempY(1:this%Site2%NPart), VSTempZ(1:this%Site2%NPart)
    real(RK)          :: VSuTempX(1:this%Site2%NPart), VSuTempY(1:this%Site2%NPart), VSuTempZ(1:this%Site2%NPart)
    real(RK)          :: VBTempX(1:this%Site2%NPart), VBTempY(1:this%Site2%NPart), VBTempZ(1:this%Site2%NPart)
    real(RK)          :: CTempX(1:this%Site2%NPart), CTempY(1:this%Site2%NPart), CTempZ(1:this%Site2%NPart)
    real(RK)          :: tuTempX(1:this%Site2%NPart), tuTempY(1:this%Site2%NPart), tuTempZ(1:this%Site2%NPart)
    real(RK)          :: tlTempX(1:this%Site2%NPart), tlTempY(1:this%Site2%NPart), tlTempZ(1:this%Site2%NPart)
    real(RK)          :: tdTempX(1:this%Site2%NPart), tdTempY(1:this%Site2%NPart), tdTempZ(1:this%Site2%NPart)
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VBxi, VByi, VBzi
    real(RK)          :: Cxi,  Cyi,  Czi
    real(RK)          :: tuxi,  tuyi,  tuzi
    real(RK)          :: tlxi,  tlyi,  tlzi
    real(RK)          :: tdxi,  tdyi,  tdzi
    real(RK)          :: txii,  tyii , tzii
    real(RK)          :: txir , tyir , tzir
    real(RK)          :: FTXi , FTYi , FTZi
    real(RK)          :: UU
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33

    VSTempX(:)=0._RK
    VSTempY(:)=0._RK
    VSTempZ(:)=0._RK
    VSuTempX(:)=0._RK
    VSuTempY(:)=0._RK
    VSuTempZ(:)=0._RK
    VBTempX(:)=0._RK
    VBTempY(:)=0._RK
    VBTempZ(:)=0._RK
    CTempX(:)=0._RK
    CTempY(:)=0._RK
    CTempZ(:)=0._RK
    tuTempX(:)=0._RK
    tuTempY(:)=0._RK
    tuTempz(:)=0._RK
    tlTempX(:)=0._RK
    tlTempY(:)=0._RK
    tlTempz(:)=0._RK
    tdTempX(:)=0._RK
    tdTempY(:)=0._RK
    tdTempz(:)=0._RK
#endif

    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ    

    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    momTempX(:)=0._RK
    momTempY(:)=0._RK
    momTempZ(:)=0._RK    
    EPotLocal=0._RK
    VirialLocal=0._RK
    d2EpotdV2Local= 0._RK

!$OMP PARALLEL &
!$OMP PRIVATE (Epsilon,  RCutoffSquared) &
!$OMP PRIVATE (RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (OX1, OY1, OZ1, OX2, OY2, OZ2) &
!$OMP PRIVATE (FX1, FY1, FZ1, TX1, TY1, TZ1) &
!$OMP PRIVATE (Plen2,sitecorr, PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (RXi, RYi, RZi, OXi, OYi, OZi,  FXi, FYi, FZi) &
!$OMP PRIVATE (TXi, TYi, TZi, PXi, PYi, PZi,  RXij, RYij, RZij) &
!$OMP PRIVATE (OXj, OYj, OZj, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE (eX, eY, eZ, RijSquared, RijInv, Rij5Inv) &
!$OMP PRIVATE (CosThetai, CosThetaj, CosGammaij) &
!$OMP PRIVATE (CosThetaiSquared, CosThetajSquared) &
!$OMP PRIVATE (dCosThetai, dCosThetaj, dCosGammaij, Tmp) &
!$OMP PRIVATE (EPotLocal1, SameComponent) &
#if  TRANS == 1
!$OMP PRIVATE(VSx, VSy, VSz ,VSux,VSuy,VSuz, VBx, VBy, VBz, Cx , Cy , Cz) &
!$OMP PRIVATE( tux , tuy , tuz, tlx , tly , tlz, tdx , tdy , tdz) &
!$OMP PRIVATE( q1, q2, q3, q4, VSxi, VSyi, VSzi, VSuxi,VSuyi,VSuzi) &
!$OMP PRIVATE( VBxi, VByi, VBzi, Cxi,  Cyi,  Czi, tuxi,  tuyi,  tuzi, tlxi,  tlyi,  tlzi) &
!$OMP PRIVATE(  tdxi,  tdyi,  tdzi, txii,  tyii , tzii, txir ,  tyir  , tzir ) &
!$OMP PRIVATE(   Uxi,  Uyi, Uzi, FTXi , FTYi , FTZi) &
!$OMP PRIVATE( A11, A12, A13, A21, A22, A23, A31, A32, A33) &
#endif
#if MPI_VER > 0
!$OMP PRIVATE ( N1, N2, i0, ji, EvenN) &
#endif
!$OMP PRIVATE (i, j, k, i1, j0, j1)

    ! Assign local variables
    SameComponent = this%SameComponent
#if MPI_VER > 0
    N1 = this%Site2%NPart
    N2 = N1 / 2
    EvenN = mod( N1, 2 ) == 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
    j1 = this%Site2%NPart
#endif
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared

    ! Assign pointers
    RX1 => this%Site1%RX
    RY1 => this%Site1%RY
    RZ1 => this%Site1%RZ
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    FX1 => this%Site1%FX
    FY1 => this%Site1%FY
    FZ1 => this%Site1%FZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

#if  TRANS == 1
    !TRANSPORT_start
    VSx1 => this%Site1%vsQx
    VSy1 => this%Site1%vsQy
    VSz1 => this%Site1%vsQz
    VSx2 => this%Site2%vsQx
    VSy2 => this%Site2%vsQy
    VSz2 => this%Site2%vsQz
    VBx1 => this%Site1%vbQx
    VBy1 => this%Site1%vbQy
    VBz1 => this%Site1%vbQz
    VBx2 => this%Site2%vbQx
    VBy2 => this%Site2%vbQy
    VBz2 => this%Site2%vbQz
    VSux1=> this%Site1%vsuQx
    VSuy1=> this%Site1%vsuQy
    VSuz1=> this%Site1%vsuQz
    VSux2=> this%Site2%vsuQx
    VSuy2=> this%Site2%vsuQy
    VSuz2=> this%Site2%vsuQz
    Cx1  => this%Site1%cQx
    Cy1  => this%Site1%cQy
    Cz1  => this%Site1%cQz
    Cx2  => this%Site2%cQx
    Cy2  => this%Site2%cQy
    Cz2  => this%Site2%cQz
    tux1 => this%Site1%tuQx
    tuy1 => this%Site1%tuQy
    tuz1 => this%Site1%tuQz
    tux2 => this%Site2%tuQx
    tuy2 => this%Site2%tuQy
    tuz2 => this%Site2%tuQz
    tlx1 => this%Site1%tlQx
    tly1 => this%Site1%tlQy
    tlz1 => this%Site1%tlQz
    tlx2 => this%Site2%tlQx
    tly2 => this%Site2%tlQy
    tlz2 => this%Site2%tlQz
    tdx1 => this%Site1%tdQx
    tdy1 => this%Site1%tdQy
    tdz1 => this%Site1%tdQz
    tdx2 => this%Site2%tdQx
    tdy2 => this%Site2%tdQy
    tdz2 => this%Site2%tdQz
    q1  => this%Site1%Q0r(:, 1)
    q2  => this%Site1%Q0r(:, 2)
    q3  => this%Site1%Q0r(:, 3)
    q4  => this%Site1%Q0r(:, 4)
!TRANSPORT_END
#endif

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)

!CDIR NODEP
#if  TRANS == 1
        !TRANSPORT_start
        VSxi= VSx1(i)
        VSyi= VSy1(i)
        VSzi= VSz1(i)
        VBxi= VBx1(i)
        VByi= VBy1(i)
        VBzi= VBz1(i)
        VSuxi= VSux1(i)
        VSuyi= VSuy1(i)
        VSuzi= VSuz1(i)
        Cxi = Cx1(i)
        Cyi = Cy1(i)
        Czi = Cz1(i)
        tuxi = tux1(i)
        tuyi = tuy1(i)
        tuzi = tuz1(i)
        tlxi = tlx1(i)
        tlyi = tly1(i)
        tlzi = tlz1(i)
        tdxi = tdx1(i)
        tdyi = tdy1(i)
        tdzi = tdz1(i)
        A11 = q1(i)**2 + q2(i)**2 - q3(i)**2 - q4(i)**2
        A12 = 2._RK * (q2(i) * q3(i) + q1(i) * q4(i))
        A13 = 2._RK * (q2(i) * q4(i) - q1(i) * q3(i))
        A21 = 2._RK * (q2(i) * q3(i) - q1(i) * q4(i))
        A22 = q1(i)**2 - q2(i)**2 + q3(i)**2 - q4(i)**2
        A23 = 2._RK * (q3(i) * q4(i) + q1(i) * q2(i))
        A31 = 2._RK * (q2(i) * q4(i) + q1(i) * q3(i))
        A32 = 2._RK * (q3(i) * q4(i) - q1(i) * q2(i))
        A33 = q1(i)**2 - q2(i)**2 - q3(i)**2 + q4(i)**2
        !TRANSPORT_END
#endif

#if OSMOP == 2
loop0:  do m=1,NBinsDen
          if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PXi < real(m)/NBinsDen-0.5_RK) then
              Bin1=m
              exit loop0
            end if
          end if
        end do loop0
#endif
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          CosThetaiSquared = CosThetai**2
          CosThetajSquared = CosThetaj**2
          Tmp = CosGammaij - 5._RK * CosThetai * CosThetaj

#if ARCH == 1
          Rij5Inv = Epsilon * RijInv * (RijInv**2)**2
#else
          Rij5Inv = Epsilon * RijInv**5
#endif
          EPotLocal1 = Rij5Inv * (1._RK - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&                      - 15._RK * CosThetaiSquared * CosThetajSquared + 2._RK * Tmp**2)

          EPotLocal = EPotLocal + EPotLocal1
          dCosThetai = Rij5Inv * (-10._RK * CosThetai - 30._RK * CosThetai * CosThetajSquared &
&                                - 20._RK * CosThetaj * Tmp)

          dCosThetaj = Rij5Inv * (-10._RK * CosThetaj - 30._RK * CosThetaj * CosThetaiSquared &
&                                - 20._RK * CosThetai * Tmp)

          dCosGammaij = 4._RK * Rij5Inv * Tmp
          Tmp = -5._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)

          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(j) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(j) < real(m)/NBinsDen-0.5_RK) then
                Bin2=m 
                exit loop2
              end if
            end if
          end do loop2
          tempMin = min(Bin1, Bin2)
          tempMax = max(Bin1, Bin2)
          if(abs(PXij) .le. 0.5_RK) then
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(tempMax-tempMin+1._RK) 
              do m = tempMin, tempMax
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          else
              VirialPart = (FXij * PXij + FYij * PYij + FZij * PZij)/(NBinsDen-tempMax+tempMin+1._RK) 
              do m = 1, tempMin
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
              do m = tempMax, NBinsDen
                this%VirialProfile(m) = this%VirialProfile(m) + VirialPart
              end do
          end if
#endif
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(35._RK*sitecorr*sitecorr-5._RK*Plen2*RijInv*RijInv)*Ninth      !xxxx9  QQ T

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij

          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
          momTempX(j) = momTempX(j) - eX * dCosThetaj - OXi * dCosGammaij
          momTempY(j) = momTempY(j) - eY * dCosThetaj - OYi * dCosGammaij  
          momTempZ(j) = momTempZ(j) - eZ * dCosThetaj - OZi * dCosGammaij   

#if  TRANS == 1
!TRANSPORT_start
          VSxi   = VSxi + 0.5 * FXij * PYij
          VSyi   = VSyi + 0.5 * FXij * PZij
          VSzi   = VSzi + 0.5 * FYij * PZij
          VBxi   = VBxi + 0.5 * FXij * PXij
          VByi   = VByi + 0.5 * FYij * PYij
          VBzi   = VBzi + 0.5 * FZij * PZij
          VSuxi  = VSuxi+ 0.5 * FYij * PXij
          VSuyi  = VSuyi+ 0.5 * FZij * PXij
          VSuzi  = VSuzi+ 0.5 * FZij * PYij
          VSTempX(j) = VSTempX(j) + 0.5*FXij * PYij
          VSTempY(j) = VSTempY(j) + 0.5*FXij * PZij
          VSTempZ(j) = VSTempZ(j) + 0.5*FYij * PZij
          VBTempX(j) = VBTempX(j) + 0.5*FXij * PXij
          VBTempY(j) = VBTempY(j) + 0.5*FYij * PYij
          VBTempZ(j) = VBTempZ(j) + 0.5*FZij * PZij
          VSuTempX(j)= VSuTempX(j)+ 0.5*FYij * PXij
          VSuTempY(j)= VSuTempY(j)+ 0.5*FZij * PXij
          VSuTempZ(j)= VSuTempZ(j)+ 0.5*FZij * PYij

          UU     = 0.5 * EPotLocal1
          Cxi    = Cxi  + UU
          Cyi    = Cyi  + UU
          Czi    = Czi  + UU
          CTempX(j) = CTempX(j) + UU
          CTempY(j) = CTempY(j) + UU
          CTempZ(j) = CTempZ(j) + UU

          FTXi   = - eX * dCosThetai - OXj * dCosGammaij
          FTYi   = - eY * dCosThetai - OYj * dCosGammaij
          FTZi   = - eZ * dCosThetai - OZj * dCosGammaij

          txii   = OYi * FTZi - OZi * FTYi
          tyii   = OZi * FTXi - OXi * FTZi
          tzii   = OXi * FTYi - OYi * FTXi
          txir   = A11 * txii + A12 * tyii + A13 * tzii
          tyir   = A21 * txii + A22 * tyii + A23 * tzii
          tzir   = A31 * txii + A32 * tyii + A33 * tzii
          tuxi   = tuxi + 0.5 * PXij* tyir
          tuyi   = tuyi + 0.5 * PXij* tzir
          tuzi   = tuzi + 0.5 * PYij* tzir
          tlxi   = tlxi + 0.5 * PYij* txir
          tlyi   = tlyi + 0.5 * PZij* txir
          tlzi   = tlzi + 0.5 * PZij* tyir
          tdxi   = tdxi + 0.5 * PXij* txir
          tdyi   = tdyi + 0.5 * PYij* tyir
          tdzi   = tdzi + 0.5 * PZij* tzir
          tuTempX(j)= tuTempX(j) + 0.5*PXij*tyi
          tuTempY(j)= tuTempY(j) + 0.5*PXij*tzi
          tuTempZ(j)= tuTempZ(j) + 0.5*PYij*tzi
          tlTempX(j)= tlTempX(j) + 0.5*PYij*txi
          tlTempY(j)= tlTempY(j) + 0.5*PZij*txi
          tlTempZ(j)= tlTempZ(j) + 0.5*PZij*tyi
          tdTempX(j)= tdTempX(j) + 0.5*PXij*txi
          tdTempY(j)= tdTempY(j) + 0.5*PYij*tyi
          tdTempZ(j)= tdTempZ(j) + 0.5*PZij*tzi

      
          !TRANSPORT_END
#endif
        end do loop1

        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi

#if  TRANS == 1
        !TRANSPORT_start
        VSx1(i) = VSxi
        VSy1(i) = VSyi
        VSz1(i) = VSzi
        VBx1(i) = VBxi
        VBy1(i) = VByi
        VBz1(i) = VBzi
        VSux1(i)= VSuxi
        VSuy1(i)= VSuyi
        VSuz1(i)= VSuzi
        Cx1(i)  = Cxi
        Cy1(i)  = Cyi
        Cz1(i)  = Czi
        tux1(i) = tuxi
        tuy1(i) = tuyi
        tuz1(i) = tuzi
        tlx1(i) = tlxi
        tly1(i) = tlyi
        tlz1(i) = tlzi
        tdx1(i) = tdxi
        tdy1(i) = tdyi
        tdz1(i) = tdzi
        !TRANSPORT_END
#endif
      end do
!$OMP END DO

    else ! Site-site cutoff

      ! Loop over molecules
!$OMP DO REDUCTION(+:forceTempX,forceTempY,forceTempZ,EPotLocal,VirialLocal,d2EpotdV2Local) &
!$OMP REDUCTION(+:momTempX, momTempY, momTempZ)
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, merge( i1 - 1, i1, SameComponent )
#endif

        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)
        TXi = TX1(i)
        TYi = TY1(i)
        TZi = TZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)

#if MPI_VER > 0
        if( SameComponent ) then
          j0 = i + 1
          j1 = i + N2
          if( EvenN .and. i > N2 ) j1 = j1 - 1
        else
          j0 = 1
          j1 = N1
        end if

!CDIR NODEP
loop3:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )
#else
        j0 = merge( i + 1, 1, SameComponent )
!CDIR NODEP
loop3:  do j = j0, j1
#endif
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = (PXij - anint( RXij )) * BoxLength
          PYij = (PYij - anint( RYij )) * BoxLength
          PZij = (PZij - anint( RZij )) * BoxLength
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij

          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          CosThetaiSquared = CosThetai**2
          CosThetajSquared = CosThetaj**2
          Tmp = CosGammaij - 5._RK * CosThetai * CosThetaj
#if ARCH == 1
          Rij5Inv = Epsilon * RijInv * (RijInv**2)**2
#else
          Rij5Inv = Epsilon * RijInv**5
#endif
          EPotLocal1 = Rij5Inv * (1._RK - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&           - 15._RK * CosThetaiSquared * CosThetajSquared + 2._RK * Tmp**2)

          EPotLocal = EPotLocal + EPotLocal1

          dCosThetai = Rij5Inv * (-10._RK * CosThetai - 30._RK * CosThetai * CosThetajSquared &
&                                - 20._RK * CosThetaj * Tmp)

          dCosThetaj = Rij5Inv * (-10._RK * CosThetaj - 30._RK * CosThetaj * CosThetaiSquared &
&                                - 20._RK * CosThetai * Tmp)

          dCosGammaij = 4._RK * Rij5Inv * Tmp
          Tmp = -5._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)

          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)

          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(35._RK*sitecorr*sitecorr-5._RK*Plen2*RijInv*RijInv)*Ninth      !xxxx9  QQ ss T

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(j) = forceTempX(j) - FXij
          forceTempY(j) = forceTempY(j) - FYij
          forceTempZ(j) = forceTempZ(j) - FZij

          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
          momTempX(j) = momTempX(j) - eX * dCosThetaj - OXi * dCosGammaij
          momTempY(j) = momTempY(j) - eY * dCosThetaj - OYi * dCosGammaij  
          momTempZ(j) = momTempZ(j) - eZ * dCosThetaj - OZi * dCosGammaij   

        end do loop3
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do
!$OMP END DO

    end if
!$OMP END PARALLEL
    FX2 = FX2 + forceTempX
    FY2 = FY2 + forceTempY
    FZ2 = FZ2 + forceTempZ
    TX2 = TX2 + momTempX                                 
    TY2 = TY2 + momTempY
    TZ2 = TZ2 + momTempZ

#if  TRANS == 1
   VSx2 = VSx2 + VSTempX
   VSy2 = VSy2 + VSTempY
   VSz2 = VSz2 + VSTempZ

   VSux2 = VSux2 + VSuTempX
   VSuy2 = VSuy2 + VSuTempY
   VSuz2 = VSuz2 + VSuTempZ

   VBx2 = VBx2 + VBTempX
   VBy2 = VBy2 + VBTempY
   VBz2 = VBz2 + VBTempZ

   Cx2 = Cx2 + CTempX
   Cy2 = Cy2 + CTempY
   Cz2 = Cz2 + CTempZ

   tux2 = tux2 + tuTempX
   tuy2 = tuy2 + tuTempY
   tuz2 = tuz2 + tuTempZ

   tdx2 = tdx2 + tdTempX
   tdy2 = tdy2 + tdTempY
   tdz2 = tdz2 + tdTempZ

   tlx2 = tlx2 + tlTempX
   tly2 = tly2 + tlTempY
   tlz2 = tlz2 + tlTempZ

#endif



    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotQQ_Force_Trans



!==============================================================!
!  Subroutine TPotQQ_ChemicalPotential                         !
!==============================================================!

  subroutine TPotQQ_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleQuadrupole) :: this
    real(RK), pointer, contiguous              :: EPotTest(:)
    real(RK), intent(in)           :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij5Inv
    real(RK)          :: CosThetai, CosThetaj, CosGammaij
    real(RK)          :: CosThetaiSquared, CosThetajSquared
    real(RK)          :: Tmp
    real(RK)          :: EPotLocal
    integer           :: i, j, k, i1, j1

    ! Assign local variables
    i1 = this%Site1%NTest
    j1 = this%Site2%NPart
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

    ! Assign pointers
    RX1 => this%Site1%RXTest
    RY1 => this%Site1%RYTest
    RZ1 => this%Site1%RZTest
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    OX1 => this%Site1%OXTest
    OY1 => this%Site1%OYTest
    OZ1 => this%Site1%OZTest
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    PX1 => this%Site1%PXTest
    PY1 => this%Site1%PYTest
    PZ1 => this%Site1%PZTest
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (RXi,RYi,RZi,PXi,PYi,PZi) &
!$OMP PRIVATE (OXi, OYi, OZi, OXj, OYj, OZj) &
!$OMP PRIVATE (RXij,RYij,RZij,PXij,PYij,PZij) &
!$OMP PRIVATE (CosThetai, CosThetaj, CosGammaij) &
!$OMP PRIVATE (CosThetaiSquared, CosThetajSquared,RijSquared,RijInv, Rij5Inv) &
!$OMP PRIVATE (Tmp,eX,eY,eZ) &
!$OMP PRIVATE (EPotLocal,i,j,k) 

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
!$OMP DO 
      do i = 1, i1
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        EPotLocal = 0._RK


!CDIR NODEP
loop1:  do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij

          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop1
          end if

          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )

          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          CosThetaiSquared = CosThetai**2
          CosThetajSquared = CosThetaj**2
          Tmp = CosGammaij - 5._RK * CosThetai * CosThetaj

#if ARCH == 1
          Rij5Inv = Epsilon * RijInv * (RijInv**2)**2
#else
          Rij5Inv = Epsilon * RijInv**5
#endif
          EPotLocal = EPotLocal + Rij5Inv * (1._RK - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&                     - 15._RK * CosThetaiSquared * CosThetajSquared + 2._RK * Tmp**2)
        end do loop1

        EPotTest(i) = EPotTest(i) + EPotLocal
      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over test particles
!$OMP DO
      do i = 1, i1
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        EPotLocal = 0._RK
!CDIR NODEP
loop2:  do j = 1, j1
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij*RXij + RYij*RYij + RZij*RZij
          if( RijSquared >= RCutoffSquared ) cycle loop2
          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop2
          end if

          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

          RijInv = 1._RK / sqrt( RijSquared )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          CosThetaiSquared = CosThetai**2
          CosThetajSquared = CosThetaj**2
          Tmp = CosGammaij - 5._RK * CosThetai * CosThetaj

#if ARCH == 1
          Rij5Inv = Epsilon * RijInv * (RijInv**2)**2
#else
          Rij5Inv = Epsilon * RijInv**5
#endif
          EPotLocal = EPotLocal + (Rij5Inv * (1._RK - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&           - 15._RK * CosThetaiSquared * CosThetajSquared + 2._RK * Tmp**2))
        end do loop2

        EPotTest(i) = EPotTest(i) + EPotLocal
      end do
!$OMP END DO
    end if
!$OMP END PARALLEL

  end subroutine TPotQQ_ChemicalPotential



  subroutine TPoterfc_approx(this,in,approx_out)

  type(TPotChargeCharge)      :: this
  real(RK),intent(in)     :: in
  real(RK),intent(out)    :: approx_out

! Local variables
  real(RK)                :: C1,C2,C3,C4,C5,P

  C1 =  0.254829592
  C2 = -0.284496736
  C3 =  1.421413741
  C4 = -1.453152027
  C5 =  1.061405429
  P  =  0.3275911

  argu = 1._RK / (1._RK + P*in)
  approx_out = argu*(C1+argu*(C2+argu*(C3+argu*(C4+argu*C5))))*exp(-in**2)

  end subroutine TPoterfc_approx


end module ms2_potential
