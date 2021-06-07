!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 2.0                !
!  (c) 2014 by TU Kaiserslautern                               !
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

#include "mathMacros.F90"

module ms2_potential

  use ms2_molecule
  use ms2_site



!==============================================================!
!  Type TPotLJ126LJ126                                         !
!==============================================================!

  type TPotLJ126LJ126

    type(TSiteLJ126), pointer :: Site1, Site2
    integer, pointer          :: NUnit1, NUnit2
    real(RK)                  :: Sigma, Epsilon
    real(RK)                  :: RCutoffSquared, RCutoffSquaredScaled
    real(RK)                  :: EPotCorr, VirialCorr, d2EpotdV2Corr,EPotTestCorr
    logical                   :: SameComponent
    logical                   :: potintra15, potintra14
    real(RK)                  :: SigmaSquared
    real(RK)                  :: Epsilon4, Epsilon48
    real(RK)                  :: BoxlengthInv, BoxLengthThird
    real(RK)                  :: ScaleLJ14
    integer, pointer, contiguous          :: NInCutoff(:), CutoffPartner(:, :)
    integer, pointer, contiguous          :: RDFSum(:)
#if OSMOP == 2
    real(RK), pointer, contiguous         :: VirialProfile(:)
#endif
#ifdef ABL
    real(RK),pointer, contiguous          :: AblEpsCorr(:,:)
    real(RK),pointer, contiguous          :: AblSigCorr(:,:)
#endif

  end type TPotLJ126LJ126

  interface Construct
    module procedure TPotLJLJ_Construct
  end interface

  interface Destruct
    module procedure TPotLJLJ_Destruct
  end interface

  interface Force
    module procedure TPotLJLJ_Force
  end interface

  interface Get_RDF
    module procedure TPotLJLJ_RDF
  end interface

  interface Force_Trans
    module procedure TPotLJLJ_Force_Trans
  end interface

  interface ChemicalPotential
    module procedure TPotLJLJ_ChemicalPotential
  end interface

  interface Energy
    module procedure TPotLJLJ_Energy
  end interface

  interface UpdateBoxLength
    module procedure TPotLJLJ_UpdateBoxLength
  end interface



!==============================================================!
!  Type TPotChargeCharge                                       !
!==============================================================!

  type TPotChargeCharge

    type(TSiteCharge), pointer :: Site1, Site2
    integer, pointer           :: NUnit1, NUnit2
    real(RK)                   :: Epsilon
    real(RK)                   :: RShieldSquared
    real(RK)                   :: RCutoffSquared
    real(RK)                   :: RFConstant
    real(RK)                   :: ScaleEl14
    logical                    :: SameComponent
    logical                    :: potintra15, potintra14
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

  interface Energy
    module procedure TPotCC_Energy
    module procedure TPotCC_Energy_Ewald
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
    integer, pointer           :: NUnit1, NUnit2
    real(RK)                   :: Epsilon
    real(RK)                   :: RShieldSquared
    real(RK)                   :: RCutoffSquared
    real(RK)                   :: ScaleEl14
    logical                    :: SameComponent
    logical                    :: potintra15, potintra14
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

  interface Energy
    module procedure TPotCD_Energy
  end interface



!==============================================================!
!  Type TPotChargeQuadrupole                                   !
!==============================================================!

  type TPotChargeQuadrupole

    type(TSiteCharge), pointer     :: Site1
    type(TSiteQuadrupole), pointer :: Site2
    integer, pointer               :: NUnit1, NUnit2
    real(RK)                       :: Epsilon
    real(RK)                       :: RShieldSquared
    real(RK)                       :: RCutoffSquared
    real(RK)                       :: ScaleEl14
    logical                        :: SameComponent
    logical                        :: potintra15, potintra14
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

  interface Energy
    module procedure TPotCQ_Energy
  end interface


!==============================================================!
!  Type TPotDipoleCharge                                       !
!==============================================================!

  type TPotDipoleCharge

    type(TSiteDipole), pointer :: Site1
    type(TSiteCharge), pointer :: Site2
    integer, pointer           :: NUnit1, NUnit2
    real(RK)                   :: Epsilon
    real(RK)                   :: RShieldSquared
    real(RK)                   :: RCutoffSquared
    real(RK)                   :: ScaleEl14
    logical                    :: SameComponent
    logical                    :: potintra15, potintra14
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

  interface Energy
    module procedure TPotDC_Energy
  end interface


!==============================================================!
!  Type TPotDipoleDipole                                       !
!==============================================================!

  type TPotDipoleDipole

    type(TSiteDipole), pointer :: Site1, Site2
    integer, pointer           :: NUnit1, NUnit2
    real(RK)                   :: Epsilon
    real(RK)                   :: RCutoffSquared
    real(RK)                   :: RShieldSquared
    real(RK)                   :: RFConstant
    real(RK)                   :: ScaleEl14
    logical                    :: SameComponent
    logical                    :: potintra15, potintra14
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

  interface Energy
    module procedure TPotDD_Energy
  end interface


!==============================================================!
!  Type TPotDipoleQuadrupole                                   !
!==============================================================!

  type TPotDipoleQuadrupole

    type(TSiteDipole), pointer     :: Site1
    type(TSiteQuadrupole), pointer :: Site2
    integer, pointer               :: NUnit1, NUnit2
    real(RK)                       :: Epsilon
    real(RK)                       :: RCutoffSquared
    real(RK)                       :: RShieldSquared
    real(RK)                       :: ScaleEl14
    logical                        :: SameComponent
    logical                        :: potintra15, potintra14
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

  interface Energy
    module procedure TPotDQ_Energy
  end interface

!==============================================================!
!  Type TPotQuadrupoleCharge                                   !
!==============================================================!

  type TPotQuadrupoleCharge

    type(TSiteQuadrupole), pointer :: Site1
    type(TSiteCharge), pointer     :: Site2
    integer, pointer               :: NUnit1, NUnit2
    real(RK)                       :: Epsilon
    real(RK)                       :: RShieldSquared
    real(RK)                       :: RCutoffSquared
    real(RK)                       :: ScaleEl14
    logical                        :: SameComponent
    logical                        :: potintra15, potintra14
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

  interface Energy
    module procedure TPotQC_Energy
  end interface


!==============================================================!
!  Type TPotQuadrupoleDipole                                   !
!==============================================================!

  type TPotQuadrupoleDipole

    type(TSiteQuadrupole), pointer :: Site1
    type(TSiteDipole), pointer     :: Site2
    integer, pointer               :: NUnit1, NUnit2
    real(RK)                       :: Epsilon
    real(RK)                       :: RCutoffSquared
    real(RK)                       :: RShieldSquared
    real(RK)                       :: ScaleEl14
    logical                        :: SameComponent
    logical                        :: potintra15, potintra14
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

  interface Energy
    module procedure TPotQD_Energy
  end interface


!==============================================================!
!  Type TPotQuadrupoleQuadrupole                               !
!==============================================================!

  type TPotQuadrupoleQuadrupole

    type(TSiteQuadrupole), pointer :: Site1, Site2
    integer, pointer               :: NUnit1, NUnit2
    real(RK)                       :: Epsilon
    real(RK)                       :: RCutoffSquared
    real(RK)                       :: RShieldSquared
    real(RK)                       :: ScaleEl14
    logical                        :: SameComponent
    logical                        :: potintra15, potintra14
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

  interface Energy
    module procedure TPotQQ_Energy
  end interface


!==============================================================!
!  Type TPotBond                                               !
!==============================================================!

  type TPotBond

    type(TIdfBond), pointer   :: Bond
    integer                   :: Site1, Site2
    integer                   :: Unit1, Unit2
    real(RK),pointer          :: ForConst
    real(RK)                  :: R0

  end type TPotBond

   interface Construct
     module procedure TPotBond_Construct
   end interface

  interface Destruct
    module procedure TPotBond_Destruct
  end interface

  interface Force
    module procedure TPotBond_Force
  end interface

  interface Energy
    module procedure TPotBond_Energy
  end interface


!==============================================================!
!  Type TPotAngle                                              !
!==============================================================!

  type TPotAngle


    type(TIdfAngle), pointer  :: Angle
    integer                   :: Site1, Site2, Site3
    integer                   :: Unit1, Unit2, Unit3
    real(RK),pointer          :: ForConst
    real(RK)                  :: Angle0
    logical,pointer           :: orientation1, orientation2

  end type TPotAngle

   interface Construct
     module procedure TPotAngle_Construct
   end interface

  interface Destruct
    module procedure TPotAngle_Destruct
  end interface

  interface Force
    module procedure TPotAngle_Force
  end interface

  interface Energy
    module procedure TPotAngle_Energy
  end interface


!==============================================================!
!  Type TPotDihedral                                           !
!==============================================================!

  type TPotDihedral


    type(TIdfDihedral), pointer  :: Dihedral
    integer                      :: Site1, Site2, Site3, Site4
    integer                      :: Unit1, Unit2, Unit3, Unit4
    integer                      :: nmax
    real(RK),pointer, contiguous :: ForConst(:)
    real(RK),pointer, contiguous :: gamma0(:)

  end type TPotDihedral

   interface Construct
     module procedure TPotDihedral_Construct
   end interface

  interface Destruct
    module procedure TPotDihedral_Destruct
  end interface

  interface Force
    module procedure TPotDihedral_Force
  end interface

  interface Energy
    module procedure TPotDihedral_Energy
  end interface

  type idfPotentialEnergies

      real(RK) :: EPotIntra
      real(RK) :: EPotIntra_Bond
      real(RK) :: EPotIntra_Angle
      real(RK) :: EPotIntra_Dihedral
      real(RK) :: EPotIntra_Nonbonded
      real(RK) :: EPotInter

  end type idfPotentialEnergies

contains



!==============================================================!
!  Subroutine TPotLJLJ_Construct                               !
!==============================================================!

  subroutine TPotLJLJ_Construct( this, i1, i2, j1, j2, Molecule1, Molecule2, RCutoff, ScaleSigma, ScaleEpsilon )

    implicit none

    ! Declare arguments
    type(TPotLJ126LJ126)        :: this
    integer, intent(in)         :: i1, i2, j1, j2
    type(TMolecule), intent(in) :: Molecule1, Molecule2
    real(RK), intent(in)        :: RCutoff
    real(RK), intent(in)        :: ScaleSigma, ScaleEpsilon

    ! Declare local variables
    real(RK) :: RCutoff3Inv, RCutoff9Inv
    real(RK) :: tau, tau1, tau2
    integer :: k, ende


    ! Construct potential
    this%Site1 => Molecule1%SiteLJ126(j1)
    this%NUnit1 => Molecule1%NUnit
    this%Site2 => Molecule2%SiteLJ126(j2)
    this%NUnit2 => Molecule2%NUnit
    this%SameComponent = i1 == i2
    this%Sigma = .5_RK * (this%Site1%sig + this%Site2%sig)
    this%Epsilon = sqrt(this%Site1%eps * this%Site2%eps)

    if( .not. this%SameComponent ) then
      this%Sigma = this%Sigma * ScaleSigma
      this%Epsilon = this%Epsilon * ScaleEpsilon
    end if
    this%Epsilon4 = 4._RK * this%Epsilon

    ! if this potential is intra
    if (this%SameComponent .and. Molecule1%hasIntraLJEl ) then
      ende = size (Molecule1%IntLJ15(:,1))
      do k=1, ende
        if (Molecule1%IntLJ15(k,1)==this%Site1%SiteId .and. Molecule1%IntLJ15(k,2)==this%Site2%SiteId) then
           this%potintra15 = .true.
        end if
      end do
      if (LJEl14 .and. .not. this%potintra15 ) then
        ende = size (Molecule1%IntLJ14(:,1))
        do k=1, ende
          if (Molecule1%IntLJ14(k,1)==this%Site1%SiteId .and. Molecule1%IntLJ14(k,2)==this%Site2%SiteId) then
            this%potintra14 = .true.
            this%ScaleLJ14 = Molecule1%ScaleLJ14(k)
          end if
        end do
      end if
    else
      this%potintra15 = .false.
      this%potintra14 = .false.
    end if

    ! Calculate long-range corrections
    this%RCutoffSquared = RCutoff**2
    tau1 = sqrt( sum( this%Site1%r(:)**2 ))
    tau2 = sqrt( sum( this%Site2%r(:)**2 ))
    tau = max( tau1, tau2 )

    if( (CutoffMode .eq. CenterofMass) .and. (tau > 1E-10_RK) ) then
      if( (tau1 > 1E-10_RK) .and. (tau2 > 1E-10_RK) ) then

        this%EPotCorr = Pi8 * this%Epsilon * &
&         ( TISSu(-6, RCutoff, this%Sigma**2, tau1, tau2) &
&         - TISSu(-3, RCutoff, this%Sigma**2, tau1, tau2) )

        this%VirialCorr = Piminus83 * this%Epsilon * &
&         ( TISSp(-6, RCutoff, this%Sigma**2, tau1, tau2) &
&         - TISSp(-3, RCutoff, this%Sigma**2, tau1, tau2) )

        this%d2EpotdV2Corr = Pi89 * this%Epsilon * &
&         ( TISSd2EpotdV2(-6, RCutoff, this%Sigma**2, tau1, tau2) &
&         - TISSd2EpotdV2(-3, RCutoff, this%Sigma**2, tau1, tau2) )

#ifdef ABL
      this%AblEpsCorr(i1,j1) = this%VirialCorr / this%Epsilon
      this%AblEpsCorr(i2,j2) = this%VirialCorr / this%Epsilon

      if ( .not. this%SameComponent ) then
        this%AblEpsCorr(i1,j1) = this%AblEpsCorr(i1,j1) * this%Site2%eps / (2._RK*this%Epsilon)
        this%AblEpsCorr(i2,j2) = this%AblEpsCorr(i2,j2) * this%Site1%eps / (2._RK*this%Epsilon)
      end if

      this%AblSigCorr(i1,j1) = Piminus83 * this%Epsilon / 2._RK * &
&         ( TISSpAbl(-6, RCutoff, this%Sigma**2, tau1, tau2) &
&         - TISSpAbl(-3, RCutoff, this%Sigma**2, tau1, tau2) )

      this%AblSigCorr(i2,j2) = Piminus83 * this%Epsilon / 2._RK * &
&         ( TISSpAbl(-6, RCutoff, this%Sigma**2, tau1, tau2) &
&         - TISSpAbl(-3, RCutoff, this%Sigma**2, tau1, tau2) )
#endif
      else
        this%EPotCorr = Pi8 * this%Epsilon * &
&         ( TICSu(-6, RCutoff, this%Sigma**2, tau) &
&         - TICSu(-3, RCutoff, this%Sigma**2, tau) )

        this%VirialCorr = Piminus83 * this%Epsilon * &
&         ( TICSp(-6, RCutoff, this%Sigma**2, tau) &
&         - TICSp(-3, RCutoff, this%Sigma**2, tau) )

        this%d2EpotdV2Corr = Pi89 * this%Epsilon * &
&         ( TICSd2EpotdV2(-6, RCutoff, this%Sigma**2, tau) &
&         - TICSd2EpotdV2(-3, RCutoff, this%Sigma**2, tau) )

#ifdef ABL
      this%AblEpsCorr(i1,j1) = this%VirialCorr / this%Epsilon
      this%AblEpsCorr(i2,j2) = this%VirialCorr / this%Epsilon

      if ( .not. this%SameComponent ) then
        this%AblEpsCorr(i1,j1) = this%AblEpsCorr(i1,j1) * this%Site2%eps / (2._RK*this%Epsilon)
        this%AblEpsCorr(i2,j2) = this%AblEpsCorr(i2,j2) * this%Site1%eps / (2._RK*this%Epsilon)
      end if

      this%AblSigCorr(i1,j1) = Piminus83 * this%Epsilon / 2._RK * &
&         ( TISSpAbl(-6, RCutoff, this%Sigma**2, tau1, tau2) &
&         - TISSpAbl(-3, RCutoff, this%Sigma**2, tau1, tau2) )

      this%AblSigCorr(i2,j2) = Piminus83 * this%Epsilon / 2._RK * &
&         ( TISSpAbl(-6, RCutoff, this%Sigma**2, tau1, tau2) &
&         - TISSpAbl(-3, RCutoff, this%Sigma**2, tau1, tau2) )

#endif
      endif
    else ! Site-site cutoff or both sites in center of mass
     this%EPotCorr = Pi8 * this%Epsilon * ( TICCu(-6, RCutoff, this%Sigma**2) - TICCu(-3, RCutoff, this%Sigma**2) )

     this%VirialCorr = Piminus83 * this%Epsilon * ( TICCp(-6, RCutoff, this%Sigma**2) - TICCp(-3, RCutoff, this%Sigma**2) )

      this%d2EpotdV2Corr = Pi89 * this%Epsilon *  ( TICCd2EpotdV2(-6, RCutoff, this%Sigma**2) - TICCd2EpotdV2(-3, RCutoff, this%Sigma**2) )

    end if
    this%EPotTestCorr = 2._RK * this%EPotCorr

 
  contains

    real(RK) function TISSu( n, rc, sigma2, tau1, tau2 )

      ! Declare arguments
      integer, intent(in)  :: n
      real(RK), intent(in) :: rc, sigma2, tau1, tau2

      ! Declare local variables
      real(RK) :: tauPlus, tauMinus

      tauPlus = tau1 + tau2
      tauMinus = abs( tau1 - tau2 )

      ! Calculate angle averaged partial integral
      TISSu = - ( (rc+tauPlus)**(2*n+4) - (rc+tauMinus)**(2*n+4) &
&             - (rc-tauMinus)**(2*n+4) + (rc-tauPlus)**(2*n+4) ) * rc &
&             / ( 8._RK * sigma2**n * tau1 * tau2 &
&             * (n+1) * (2*n+3) * (2*n+4) ) &
&             + ( (rc+tauPlus)**(2*n+5) - (rc+tauMinus)**(2*n+5) &
&             - (rc-tauMinus)**(2*n+5) + (rc-tauPlus)**(2*n+5) ) &
&             / ( 8._RK * sigma2**n * tau1 * tau2 &
&             * (n+1) * (2*n+3) * (2*n+4) * (2*n+5) )

    end function TISSu



    real(RK) function TICSu( n, rc, sigma2, tau )

      ! Declare arguments
      integer, intent(in)  :: n
      real(RK), intent(in) :: rc, sigma2, tau

      ! Calculate angle averaged partial integral
      TICSu = - ( (rc+tau)**(2*n+3) - (rc-tau)**(2*n+3) ) * rc &
&             / ( 4._RK * sigma2**n * tau * (n+1) * (2*n+3) ) &
&             + ( (rc+tau)**(2*n+4) - (rc-tau)**(2*n+4) ) &
&             / ( 4._RK * sigma2**n * tau * (n+1) * (2*n+3) * (2*n+4) )

    end function TICSu


        
    real(RK) function TICCu( n, rc, sigma2 )

      ! Declare arguments
      integer, intent(in)  :: n
      real(RK), intent(in) :: rc, sigma2

      ! Calculate angle averaged partial integral
      TICCu = - ( rc**(2*n+3) ) / ( sigma2**n * (2*n+3) )

    end function TICCu



    real(RK) function TISSp( n, rc, sigma2, tau1, tau2 )

      ! Declare arguments
      integer, intent(in)  :: n
      real(RK), intent(in) :: rc, sigma2, tau1, tau2

      ! Declare local variables
      real(RK) :: tauPlus, tauMinus

      tauPlus = tau1 + tau2
      tauMinus = abs( tau1 - tau2 )

      ! Calculate angle averaged partial integral
      TISSp = - ( (rc+tauPlus)**(2*n+3) - (rc+tauMinus)**(2*n+3) &
&             - (rc-tauMinus)**(2*n+3) + (rc-tauPlus)**(2*n+3) ) * rc**2 &
&             / ( 8._RK * sigma2**n * tau1 * tau2 * (n+1) * (2*n+3) ) &
&             - 3._RK * TISSu(n,rc,sigma2,tau1,tau2)

    end function TISSp



    real(RK) function TICSp( n, rc, sigma2, tau )

      ! Declare arguments
      integer, intent(in)  :: n
      real(RK), intent(in) :: rc, sigma2, tau

      ! Calculate angle averaged partial integral
      TICSp = - ( (rc+tau)**(2*n+2) - (rc-tau)**(2*n+2) ) * rc**2 &
&             / ( 4._RK * sigma2**n * tau * (n+1) ) &
&             - 3._RK * TICSu(n,rc,sigma2,tau)

    end function TICSp


    real(RK) function TICCp( n, rc, sigma2 )

      ! Declare arguments
      integer, intent(in)  :: n
      real(RK), intent(in) :: rc, sigma2

      ! Calculate angle averaged partial integral
      TICCp = 2._RK*n*TICCu(n,rc,sigma2)

    end function TICCp


    real(RK) function TISSd2EpotdV2( n, rc, sigma2, tau1, tau2 )

      ! Declare arguments
      integer, intent(in)  :: n
      real(RK), intent(in) :: rc, sigma2, tau1, tau2

      ! Declare local variables
      real(RK) :: tauPlus, tauMinus, A, B, C, D, AN, BN, CN, DN

      tauPlus = tau1 + tau2
      tauMinus = abs( tau1 - tau2 )

      ! Calculate angle averaged partial integral
      A = 18._RK*n**2*rc**3          -6._RK *n**2*rc**2*tauPlus&
&        -3._RK *tauPlus**3          +6._RK *n*tauPlus**2*rc&
&        +4._RK *rc**3*n**3          +12._RK*rc**3&
&        +26._RK*n*rc**3             -9._RK *rc**2*tauPlus&
&        +6._RK *tauPlus**2*rc       -15._RK*rc**2*n*tauPlus

      B = 18._RK*n**2*rc**3          +6._RK *n**2*rc**2*tauPlus&
&        +3._RK *tauPlus**3          +6._RK *n*tauPlus**2*rc&
&        +4._RK *rc**3*n**3          +12._RK*rc**3&
&        +26._RK*n*rc**3             +9._RK *rc**2*tauPlus&
&        +6._RK *tauPlus**2*rc       +15._RK*rc**2*n*tauPlus

      C = 18._RK*n**2*rc**3          -15._RK*rc**2*tauMinus*n&
&        +4._RK *rc**3*n**3          +12._RK*rc**3&
&        -6._RK *n**2*rc**2*tauMinus +6._RK *n*tauMinus**2*rc&
&        -3._RK *tauMinus**3         +26._RK*n*rc**3&
&        +6._RK *tauMinus**2*rc      -9._RK *rc**2*tauMinus

      D = 18._RK*n**2*rc**3          +15._RK*rc**2*tauMinus*n&
&        +4._RK *rc**3*n**3          +12._RK*rc**3&
&        +6._RK *n**2*rc**2*tauMinus +6._RK *n*tauMinus**2*rc&
&        +3._RK *tauMinus**3         +26._RK*n*rc**3&
&        +6._RK *tauMinus**2*rc      +9._RK *rc**2*tauMinus

      AN = (1._RK/8._RK)*(rc+tauPlus)**(2*n+3)/(tau1*tau2*(n+1._RK)*(2._RK*n+3._RK)*(rc+tauPlus)*(5._RK+2._RK*n)*(2._RK+n))

      BN = (1._RK/8._RK)*(rc-tauPlus)**(2*n+3)/(tau1*tau2*(n+1._RK)*(2._RK*n+3._RK)*(rc-tauPlus)*(5._RK+2._RK*n)*(2._RK+n))

      CN =-(1._RK/8._RK)*(rc+tauMinus)**(2*n+3)/(tau1*tau2*(n+1._RK)*(2._RK*n+3._RK)*(rc+tauMinus)*(5._RK+2._RK*n)*(2._RK+n))

      DN =-(1._RK/8._RK)*(rc-tauMinus)**(2*n+3)/(tau1*tau2*(n+1._RK)*(2._RK*n+3._RK)*(rc-tauMinus)*(5._RK+2._RK*n)*(2._RK+n))

      TISSd2EpotdV2 =-(A*AN+B*BN+C*CN+D*DN)/sigma2**n -2._RK * TISSp(n,rc,sigma2,tau1,tau2)

    end function TISSd2EpotdV2


    real(RK) function TICSd2EpotdV2( n, rc, sigma2, tau )

      ! Declare arguments
      integer, intent(in)  :: n
      real(RK), intent(in) :: rc, sigma2, tau

      ! Declare local variables
      real(RK) :: A, B, AN, BN

      ! Calculate angle averaged partial integral

      A = 3._RK *tau**3         +3._RK *rc**2*tau&
&        +12._RK*n**2*rc**3     +11._RK*n*rc**3&
&        +6._RK *n**2*rc**2*tau +9._RK *n*rc**2*tau&
&        +3._RK *rc**3          +4._RK *n**3*rc**3&
&        +3._RK *rc*tau**2      +6._RK *rc*n*tau**2

      B =-3._RK *tau**3         -3._RK *rc**2*tau&
&        +12._RK*n**2*rc**3     +11._RK*n*rc**3&
&        -6._RK *n**2*rc**2*tau -9._RK *n*rc**2*tau&
&        +3._RK *rc**3          +4._RK *n**3*rc**3&
&        +3._RK *rc*tau**2      +6._RK *rc*n*tau**2

      AN = -(1._RK/4._RK)*(rc-tau)**(2*n+2)/(tau*(n+1._RK)*(rc-tau)*(2._RK+n)*(3._RK+2._RK*n))

      BN =  (1._RK/4._RK)*(rc+tau)**(2*n+2)/(tau*(n+1._RK)*(rc+tau)*(2._RK+n)*(3._RK+2._RK*n))

      TICSd2EpotdV2 =-(A*AN+B*BN)/sigma2**n - 2._RK * TICSp(n,rc,sigma2,tau)

    end function TICSd2EpotdV2


    real(RK) function TICCd2EpotdV2( n, rc, sigma2 )

      ! Declare arguments
      integer, intent(in)  :: n
      real(RK), intent(in) :: rc, sigma2

      ! Declare local variables
      real(RK) :: A, AN

      A = rc**(3+2*n)

      AN = 2._RK*n*(2._RK*n-1._RK)/(3._RK+2._RK*n)

      ! Calculate angle averaged partial integral
      TICCd2EpotdV2 =-(A*AN)/sigma2**n

    end function TICCd2EpotdV2 


    real(RK) function TISSuAbl( n, rc, sigma2, tau1, tau2 )

      ! Declare arguments
      integer, intent(in)  :: n
      real(RK), intent(in) :: rc, sigma2, tau1, tau2

      ! Declare local variables
      real(RK) :: tauPlus, tauMinus

      tauPlus = tau1 + tau2
      tauMinus = abs( tau1 - tau2 )

      ! Calculate angle averaged partial integral
      TISSuAbl = - ( (rc+tauPlus)**(2*n+4) - (rc+tauMinus)**(2*n+4) &
&                - (rc-tauMinus)**(2*n+4) + (rc-tauPlus)**(2*n+4) ) * rc * (-2*n) &
&                / ( 8._RK * sigma2**(n+1) * tau1 * tau2 &
&                * (n+1) * (2*n+3) * (2*n+4) ) &
&                + ( (rc+tauPlus)**(2*n+5) - (rc+tauMinus)**(2*n+5) &
&                - (rc-tauMinus)**(2*n+5) + (rc-tauPlus)**(2*n+5) ) * (-2*n)&
&                / ( 8._RK * sigma2**(n+1) * tau1 * tau2 &
&                * (n+1) * (2*n+3) * (2*n+4) * (2*n+5) )

    end function TISSuAbl

    real(RK) function TISSpAbl( n, rc, sigma2, tau1, tau2 )

      ! Declare arguments
      integer, intent(in)  :: n
      real(RK), intent(in) :: rc, sigma2, tau1, tau2

      ! Declare local variables
      real(RK) :: tauPlus, tauMinus

      tauPlus = tau1 + tau2
      tauMinus = abs( tau1 - tau2 )

      ! Calculate angle averaged partial integral
      TISSpAbl = - ( (rc+tauPlus)**(2*n+3) - (rc+tauMinus)**(2*n+3) &
&                - (rc-tauMinus)**(2*n+3) + (rc-tauPlus)**(2*n+3) ) *rc**2 * (-2*n)&
&                / ( 8._RK * sigma2**(n+1) * tau1 * tau2 * (n+1) * (2*n+3) ) &
&                - 3._RK * TISSuAbl(n,rc,sigma2,tau1,tau2)

    end function TISSpAbl


  end subroutine TPotLJLJ_Construct



!==============================================================!
!  Subroutine TPotLJLJ_Destruct                                !
!==============================================================!

  subroutine TPotLJLJ_Destruct( this )

    implicit none

    ! Declare arguments
    type(TPotLJ126LJ126) :: this

    ! Destroy potential
    continue

  end subroutine TPotLJLJ_Destruct

!==============================================================!
!  Subroutine TPotLJLJ_Force                                   !
!==============================================================!

#ifdef ABL
  subroutine TPotLJLJ_Force( this, EPot, Virial, EPotInter, &
&            VirialInter, EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength, &
&            VirAblSig, VirAblEps, eps1,eps2)
#else
  subroutine TPotLJLJ_Force( this, EPot, Virial, EPotInter, &
&            VirialInter, EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )
#endif

    implicit none

    ! Declare arguments
    type(TPotLJ126LJ126)     :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra_Nonbonded
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength
#ifdef ABL
    real(RK), intent(in out) :: VirAblSig
    real(RK), intent(in out) :: VirAblEps
    real(RK), intent(in out) :: eps1,eps2
#endif

    ! Declare local variables
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK)          :: SigmaSquared
    real(RK)          :: Epsilon4, Epsilon48
    real(RK)          :: RCutoffSquared
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: FXij, FYij, FZij, Fij
    real(RK)          :: RijSquared, RijSquaredInv, Rij6Inv
    real(RK)          :: EPotLocal, EPotLocal1, VirialLocal
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
    integer           :: jk, unit
    real(RK)          :: coeff

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

#ifdef ABL
    real(RK)          :: dr2Abl
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
#ifdef ABL
    VirAblSig = 0.0_RK
    VirAblEps = 0.0_RK
    eps1      = this%Site1%eps
    eps2      = this%Site2%eps
#endif
    SameComponent = this%SameComponent
    forceTempX(:)=0._RK
    forceTempY(:)=0._RK
    forceTempZ(:)=0._RK
    EPotLocal=0._RK
    VirialLocal=0._RK
    EPotLocalIntra   = 0._RK
    VirialLocalIntra = 0._RK
    EPotLocalInter   = 0._RK
    VirialLocalInter = 0._RK
    d2EpotdV2Local= 0._RK
    SigmaSquared = this%SigmaSquared
    Epsilon4 = this%Epsilon4
    Epsilon48 = this%Epsilon48
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
    
    if (this%potintra14) then
      coeff = this%ScaleLJ14
    else
      coeff = 1._RK
    end if

!$OMP PARALLEL DEFAULT(SHARED) & 
#if MPI_VER > 0
!$OMP FIRSTPRIVATE(i0, N1, N2, ji, EvenN) &
#endif
!$OMP FIRSTPRIVATE(i1, j1) &
!$OMP PRIVATE( i, j, k, j0) &
!$OMP PRIVATE(Plen2, sitecorr, EPotLocal1) &
!$OMP PRIVATE(RXi, RYi, RZi,  PXi, PYi, PZi,  FXi, FYi, FZi) &
!$OMP PRIVATE(RXij, RYij, RZij, PXij, PYij, PZij) &
!$OMP PRIVATE(FXij, FYij, FZij, Fij, RijSquared, RijSquaredInv, Rij6Inv ) 

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

        unit=this%NUnit1*(i-1)+this%Site1%UnitNumber ! Number of unit, to which this site corresponds

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit) ! Unit-partner of this unit
          if ( mod(j-this%Site2%UnitNumber, this%NUnit2)==0) then  ! choose only units, to which our Site2 correspond
            if (mod(j,this%NUnit2)==0) then
              jk = INT(j/this%NUnit2)   ! number of molecule, to which this unit correspond
            else
              jk = INT(j/this%NUnit2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = RXij - anint( PXij )
            RYij = RYij - anint( PYij )
            RZij = RZij - anint( PZij )
            PXij = PXij - anint( PXij )
            PYij = PYij - anint( PYij )
            PZij = PZij - anint( PZij )
            RijSquared = RXij**2 + RYij**2 + RZij**2
            RijSquaredInv = SigmaSquared / RijSquared
            Rij6Inv = RijSquaredInv**3
            EPotLocal1 = Rij6Inv * (Rij6Inv - 1._RK)
            EPotLocal = EPotLocal + EPotLocal1
            EPotLocalInter = EPotLocalInter + EPotLocal1
            Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
            FXij = Fij * RXij
            FYij = Fij * RYij
            FZij = Fij * RZij
            VirialLocal = VirialLocal + (PXij * FXij + PYij * FYij + PZij * FZij)
            VirialLocalInter = VirialLocalInter + (PXij * FXij + PYij * FYij + PZij * FZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)/RijSquared
          d2EpotdV2Local = d2EpotdV2Local + Epsilon4 * Rij6Inv * (12._RK*Rij6Inv  -  6._RK) * (sitecorr * sitecorr - Plen2/RijSquared)*Third*Third !xxxx LJ
          d2EpotdV2Local = d2EpotdV2Local + Epsilon4 * Rij6Inv * (156._RK*Rij6Inv - 42._RK) *  sitecorr * sitecorr*Third*Third
          FXi = FXi + FXij
          FYi = FYi + FYij
          FZi = FZi + FZij
          forceTempX(jk) = forceTempX(jk) - FXij
          forceTempY(jk) = forceTempY(jk) - FYij
          forceTempZ(jk) = forceTempZ(jk) - FZij

#ifdef ABL
          dr2Abl  = RXij**2 + RYij**2 + RZij**2
          VirAblSig = VirAblSig + Rij6Inv*(1._RK-4._RK*Rij6Inv)*(PXij*RXij+ &
&                     PYij*RYij + PZij*RZij) / dr2Abl
          VirAblEps = VirAblEps + Rij6Inv*(1._RK-2._RK*Rij6Inv)*(PXij*RXij+ &
&                     PYij*RYij + PZij*RZij) / dr2Abl
#endif
          end if
        end do loop1
        ! Include intramolecular interaction if need
        if (this%potintra15 .or. this%potintra14) then ! Michael Sch.: intra15/14 enough, .and. redundant (changed for all pot-classes)
        ! previous: if (SameComponent .and. (intra15 .or. intra14)) then
          RXij = RXi - RX2(i)
          RYij = RYi - RY2(i)
          RZij = RZi - RZ2(i)
          PXij = PXi - PX2(i)
          PYij = PYi - PY2(i)
          PZij = PZi - PZ2(i)
          RXij = RXij - anint( PXij )
          RYij = RYij - anint( PYij )
          RZij = RZij - anint( PZij )
          PXij = PXij - anint( PXij )
          PYij = PYij - anint( PYij )
          PZij = PZij - anint( PZij )
          RijSquaredInv = SigmaSquared / ( RXij**2 + RYij**2 + RZij**2 )
          Rij6Inv = RijSquaredInv**3
          EPotLocal1 = Rij6Inv * (Rij6Inv - 1._RK) * coeff
          EPotLocal = EPotLocal + EPotLocal1
          EPotLocalIntra = EPotLocalIntra + EPotLocal1
          Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv*coeff
          FXij = Fij * RXij
          FYij = Fij * RYij
          FZij = Fij * RZij
          VirialLocal = VirialLocal + (PXij * FXij + PYij * FYij + PZij * FZij)
          VirialLocalIntra = VirialLocalIntra + (PXij * FXij + PYij * FYij + PZij * FZij)
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)/RijSquared
          d2EpotdV2Local = d2EpotdV2Local + Epsilon4 * Rij6Inv * (12._RK*Rij6Inv  -  6._RK) * (sitecorr * sitecorr - Plen2/RijSquared)*Third*Third !xxxx LJ
          d2EpotdV2Local = d2EpotdV2Local + Epsilon4 * Rij6Inv * (156._RK*Rij6Inv - 42._RK) *  sitecorr * sitecorr*Third*Third
          FXi = FXi + FXij
          FYi = FYi + FYij
          FZi = FZi + FZij
          forceTempX(i) = forceTempX(i) - FXij
          forceTempY(i) = forceTempY(i) - FYij
          forceTempZ(i) = forceTempZ(i) - FZij
        end if
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
          RijSquared = RXij**2 + RYij**2 + RZij**2
          if( RijSquared >= RCutoffSquared ) cycle loop3
          RijSquaredInv = SigmaSquared / RijSquared
          Rij6Inv = RijSquaredInv**3
          EPotLocal = EPotLocal + (Rij6Inv * (Rij6Inv - 1._RK))
          EPotLocalInter = EPotLocalInter + (Rij6Inv * (Rij6Inv - 1._RK))
          Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
          FXij = Fij * RXij
          FYij = Fij * RYij
          FZij = Fij * RZij
          VirialLocal = VirialLocal + (PXij * FXij + PYij * FYij + PZij * FZij)
          VirialLocalInter = VirialLocalInter + (PXij * FXij + PYij * FYij + PZij * FZij)
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)/RijSquared
          d2EpotdV2Local = d2EpotdV2Local + Epsilon4 * Rij6Inv * (12._RK*Rij6Inv  -  6._RK) * (sitecorr * sitecorr - Plen2/RijSquared)*Third*Third !xxxx LJ SS
          d2EpotdV2Local = d2EpotdV2Local + Epsilon4 * Rij6Inv * (156._RK*Rij6Inv - 42._RK) *  sitecorr * sitecorr*Third*Third
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
   EPot = EPot + this%Epsilon4 * EPotLocal
   Virial = Virial + Third * VirialLocal * BoxLength
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:) * BoxLength
#endif
   EPotInter = EPotInter + Epsilon4 * EPotLocalInter
   VirialInter = VirialInter + Third * VirialLocalInter * BoxLength
   if (IntraLJEl) then
     EPotIntra_Nonbonded = EPotIntra_Nonbonded + Epsilon4 * EPotLocalIntra
     VirialIntra = VirialIntra + Third * VirialLocalIntra * BoxLength
   end if
   d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

#ifdef ABL
    VirAblSig = VirAblSig * Third * BoxLength * 18._RK * Epsilon4 / this%Sigma
    VirAblEps = VirAblEps * Third * BoxLength * 24._RK
#endif
  end subroutine TPotLJLJ_Force


!==============================================================!
!  Subroutine TPotLJLJ_Force_Trans                             !
!==============================================================!

#ifdef ABL
  subroutine TPotLJLJ_Force_Trans( this, EPot, Virial, EPotInter, &
&            VirialInter, EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength, &
&            VirAblSig, VirAblEps, eps1,eps2)
#else
  subroutine TPotLJLJ_Force_Trans( this, EPot, Virial, EPotInter, &
&            VirialInter, EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )
#endif

    implicit none

    ! Declare arguments
    type(TPotLJ126LJ126)     :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra_Nonbonded
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength

#ifdef ABL
    real(RK), intent(in out) :: VirAblSig
    real(RK), intent(in out) :: VirAblEps
    real(RK), intent(in out) :: eps1,eps2
#endif

    ! Declare local variables
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK)          :: SigmaSquared
    real(RK)          :: Epsilon4, Epsilon48
    real(RK)          :: RCutoffSquared
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: FXij, FYij, FZij, Fij
    real(RK)          :: RijSquared, RijSquaredInv, Rij6Inv
    real(RK)          :: EPotLocal, EPotLocal1, VirialLocal
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
    integer           :: jk, unit
    real(RK)          :: coeff
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

#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer, contiguous :: VSx(:), VSy(:), VSz(:) 
    real(RK), pointer, contiguous :: VSux(:),VSuy(:),VSuz(:)
    real(RK), pointer, contiguous :: VBx(:), VBy(:), VBz(:)
    real(RK), pointer, contiguous :: Cx(:) , Cy(:) , Cz(:)
    real(RK), pointer, contiguous :: tux(:) , tuy(:) , tuz(:)
    real(RK), pointer, contiguous :: tlx(:) , tly(:) , tlz(:)
    real(RK), pointer, contiguous :: tdx(:) , tdy(:) , tdz(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
    real(RK)          :: SigmaInvEps4
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VBxi, VByi, VBzi
    real(RK)          :: Cxi,  Cyi,  Czi
    real(RK)          :: tuxi,  tuyi,  tuzi
    real(RK)          :: tlxi,  tlyi,  tlzi
    real(RK)          :: tdxi,  tdyi,  tdzi
    real(RK)          :: txii,  tyii , tzii
    real(RK)          :: txi ,  tyi  , tzi
    real(RK)          :: UU ,  Uxi,  Uyi, Uzi, RijSInvNorm
    real(RK)          :: BoxLength2
    real(RK)          :: r1x, r1y, r1z
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33
   !TRANSPORT_END
#endif

#ifdef ABL
    real(RK)          :: dr2Abl
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
!$OMP PRIVATE( Plen2,PX1, PY1, PZ1, PX2, PY2, PZ2, FX1, FY1, FZ1 ) &
!$OMP PRIVATE(SigmaSquared, Epsilon4, Epsilon48, RCutoffSquared,EPotLocal1) &
!$OMP PRIVATE(RXi, RYi, RZi,  PXi, PYi, PZi,  FXi, FYi, FZi,  RXij, RYij, RZij, PXij, PYij, PZij) &
!$OMP PRIVATE(FXij, FYij, FZij, Fij, RijSquared, RijSquaredInv, Rij6Inv ) &
#if MPI_VER > 0
!$OMP PRIVATE(i0, N1, N2, ji, EvenN) &
#endif
#if  TRANS == 1
!$OMP PRIVATE(VSx, VSy, VSz ,VSux,VSuy,VSuz, VBx, VBy, VBz, Cx , Cy , Cz) &
!$OMP PRIVATE( tux , tuy , tuz, tlx , tly , tlz, tdx , tdy , tdz) &
!$OMP PRIVATE( q1, q2, q3, q4, SigmaInvEps4, VSxi, VSyi, VSzi, VSuxi,VSuyi,VSuzi) &
!$OMP PRIVATE( VBxi, VByi, VBzi, Cxi,  Cyi,  Czi, tuxi,  tuyi,  tuzi, tlxi,  tlyi,  tlzi) &
!$OMP PRIVATE(  tdxi,  tdyi,  tdzi, txii,  tyii , tzii, txi ,  tyi  , tzi ) &
!$OMP PRIVATE(  UU ,  Uxi,  Uyi, Uzi, RijSInvNorm, BoxLength2, r1x, r1y, r1z) &
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
#ifdef ABL
    VirAblSig = 0.0_RK
    VirAblEps = 0.0_RK
    eps1      = this%Site1%eps
    eps2      = this%Site2%eps
#endif
    SigmaSquared = this%SigmaSquared
    Epsilon4 = this%Epsilon4
    Epsilon48 = this%Epsilon48
    RCutoffSquared = this%RCutoffSquaredScaled
    EPotLocalIntra   = 0._RK
    VirialLocalIntra = 0._RK
    EPotLocalInter   = 0._RK
    VirialLocalInter = 0._RK

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

    if (this%potintra14) then
       coeff = this%ScaleLJ14
    else
       coeff = 1._RK
    end if

#if  TRANS == 1
    !TRANSPORT_start

    SigmaInvEps4 = Epsilon4/Sqrt(this%SigmaSquared)
    BoxLength2   = BoxLength**2
    VSx => this%Site1%vsLJx
    VSy => this%Site1%vsLJy
    VSz => this%Site1%vsLJz
    VBx => this%Site1%vbLJx
    VBy => this%Site1%vbLJy
    VBz => this%Site1%vbLJz
    VSux=> this%Site1%vsuLJx
    VSuy=> this%Site1%vsuLJy
    VSuz=> this%Site1%vsuLJz
    Cx  => this%Site1%cLJx
    Cy  => this%Site1%cLJy
    Cz  => this%Site1%cLJz
    tux => this%Site1%tuLJx
    tuy => this%Site1%tuLJy
    tuz => this%Site1%tuLJz
    tlx => this%Site1%tlLJx
    tly => this%Site1%tlLJy
    tlz => this%Site1%tlLJz
    tdx => this%Site1%tdLJx
    tdy => this%Site1%tdLJy
    tdz => this%Site1%tdLJz
    q1  => this%Site1%Qm0r(:, 1, 1)
    q2  => this%Site1%Qm0r(:, 2, 1)
    q3  => this%Site1%Qm0r(:, 3, 1)
    q4  => this%Site1%Qm0r(:, 4, 1)
!TRANSPORT_END
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
        VSxi= 0._RK
        VSyi= 0._RK
        VSzi= 0._RK
        VBxi= 0._RK
        VByi= 0._RK
        VBzi= 0._RK
        VSuxi= 0._RK
        VSuyi= 0._RK
        VSuzi= 0._RK
        Cxi = Cx(i)
        Cyi = Cy(i)
        Czi = Cz(i)
        tuxi = tux(i)
        tuyi = tuy(i)
        tuzi = tuz(i)
        tlxi = tlx(i)
        tlyi = tly(i)
        tlzi = tlz(i)
        tdxi = tdx(i)
        tdyi = tdy(i)
        tdzi = tdz(i)
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

        unit=this%NUnit1*(i-1)+this%Site1%UnitNumber ! Number of unit, to which this site corresponds

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit) ! Unit-partner of this unit
          if ( mod(j-this%Site2%UnitNumber, this%NUnit2)==0) then  ! choose only units, to which our Site2 correspond
            if (mod(j,this%NUnit2)==0) then
              jk = INT(j/this%NUnit2)   ! number of molecule, to which this unit correspond
            else
              jk = INT(j/this%NUnit2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = RXij - anint( PXij )
            RYij = RYij - anint( PYij )
            RZij = RZij - anint( PZij )
            PXij = PXij - anint( PXij )
            PYij = PYij - anint( PYij )
            PZij = PZij - anint( PZij )
            RijSquared = RXij**2 + RYij**2 + RZij**2
            RijSquaredInv = SigmaSquared / RijSquared
            Rij6Inv = RijSquaredInv**3
            EPotLocal1 = Rij6Inv * (Rij6Inv - 1._RK)
            EPotLocal = EPotLocal + EPotLocal1
            EPotLocalInter = EPotLocalInter + EPotLocal1
            Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
            FXij = Fij * RXij
            FYij = Fij * RYij
            FZij = Fij * RZij
            VirialLocal = VirialLocal + (PXij * FXij + PYij * FYij + PZij * FZij)
#if OSMOP == 2
loop2:      do m=1,NBinsDen
              if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
                if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
            VirialLocalInter = VirialLocalInter + (PXij * FXij + PYij * FYij + PZij * FZij)
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)/RijSquared
            d2EpotdV2Local = d2EpotdV2Local + Epsilon4 * Rij6Inv * (12._RK*Rij6Inv  -  6._RK) * (sitecorr * sitecorr - Plen2/RijSquared)*Third*Third  !xxxx LJ T
            d2EpotdV2Local = d2EpotdV2Local + Epsilon4 * Rij6Inv * (156._RK*Rij6Inv - 42._RK) *  sitecorr * sitecorr*Third*Third
            FXi = FXi + FXij
            FYi = FYi + FYij
            FZi = FZi + FZij
            forceTempX(jk) = forceTempX(jk) - FXij
            forceTempY(jk) = forceTempY(jk) - FYij
            forceTempZ(jk) = forceTempZ(jk) - FZij
#if  TRANS == 1
            !TRANSPORT_start
            VSxi   = VSxi + FXij * PYij
            VSyi   = VSyi + FXij * PZij
            VSzi   = VSzi + FYij * PZij
            VBxi   = VBxi + FXij * PXij
            VByi   = VByi + FYij * PYij
            VBzi   = VBzi + FZij * PZij
            VSuxi  = VSuxi+ FYij * PXij
            VSuyi  = VSuyi+ FZij * PXij
            VSuzi  = VSuzi+ FZij * PYij
            RijSInvNorm   = Sqrt(RijSquaredInv)
            UU   = RijSInvNorm*EPotLocal1*SigmaInvEps4
            Cxi    = Cxi  + UU*RXij
            Cyi    = Cyi  + UU*RYij
            Czi    = Czi  + UU*RZij
            txii   = r1y * FZij - r1z * FYij
            tyii   = r1z * FXij - r1x * FZij
            tzii   = r1x * FYij - r1y * FXij
            txi    = A11 * txii + A12 * tyii + A13 * tzii
            tyi    = A21 * txii + A22 * tyii + A23 * tzii
            tzi    = A31 * txii + A32 * tyii + A33 * tzii
            tuxi   = tuxi + PXij*tyi
            tuyi   = tuyi + PXij*tzi
            tuzi   = tuzi + PYij*tzi
            tlxi   = tlxi + PYij*txi
            tlyi   = tlyi + PZij*txi
            tlzi   = tlzi + PZij*tyi
            tdxi   = tdxi + PXij*txi
            tdyi   = tdyi + PYij*tyi
            tdzi   = tdzi + PZij*tzi
            !TRANSPORT_END
#endif

#ifdef ABL
          dr2Abl  = RXij**2 + RYij**2 + RZij**2
          VirAblSig = VirAblSig + Rij6Inv*(1._RK-4._RK*Rij6Inv)*(PXij*RXij+ &
&                     PYij*RYij + PZij*RZij) / dr2Abl
          VirAblEps = VirAblEps + Rij6Inv*(1._RK-2._RK*Rij6Inv)*(PXij*RXij+ &
&                     PYij*RYij + PZij*RZij) / dr2Abl
#endif
          end if
        end do loop1
        ! Include intramolecular interaction if need
        if (this%potintra15 .or. this%potintra14) then
          RXij = RXi - RX2(i)
          RYij = RYi - RY2(i)
          RZij = RZi - RZ2(i)
          PXij = PXi - PX2(i)
          PYij = PYi - PY2(i)
          PZij = PZi - PZ2(i)
          RXij = RXij - anint( PXij )
          RYij = RYij - anint( PYij )
          RZij = RZij - anint( PZij )
          PXij = PXij - anint( PXij )
          PYij = PYij - anint( PYij )
          PZij = PZij - anint( PZij )
          RijSquaredInv = SigmaSquared / ( RXij**2 + RYij**2 + RZij**2 )
          Rij6Inv = RijSquaredInv**3
          EPotLocal1 = Rij6Inv * (Rij6Inv - 1._RK) * coeff
          EPotLocal = EPotLocal + EPotLocal1
          EPotLocalIntra = EPotLocalIntra + EPotLocal1
          Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv*coeff
          FXij = Fij * RXij
          FYij = Fij * RYij
          FZij = Fij * RZij
          VirialLocal = VirialLocal + (PXij * FXij + PYij * FYij + PZij * FZij)
          VirialLocalIntra = VirialLocalIntra + (PXij * FXij + PYij * FYij + PZij * FZij)
          FXi = FXi + FXij
          FYi = FYi + FYij
          FZi = FZi + FZij
          forceTempX(i) = forceTempX(i) - FXij
          forceTempY(i) = forceTempY(i) - FYij
          forceTempZ(i) = forceTempZ(i) - FZij
        end if
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
#if  TRANS == 1
        !TRANSPORT_start
        VSx(i) = VSx(i) + VSxi *BoxLength
        VSy(i) = VSy(i) + VSyi *BoxLength
        VSz(i) = VSz(i) + VSzi *BoxLength
        VBx(i) = VBx(i) + VBxi*BoxLength
        VBy(i) = VBy(i) + VByi*BoxLength
        VBz(i) = VBz(i) + VBzi*BoxLength
        VSux(i)= VSux(i)+ VSuxi*BoxLength
        VSuy(i)= VSuy(i)+ VSuyi*BoxLength
        VSuz(i)= VSuz(i)+ VSuzi*BoxLength
        Cx(i)  = Cxi
        Cy(i)  = Cyi
        Cz(i)  = Czi
        ! Multiplication with Boxlength for the following terms already done in rx1, ...
        tux(i) = tuxi
        tuy(i) = tuyi
        tuz(i) = tuzi
        tlx(i) = tlxi
        tly(i) = tlyi
        tlz(i) = tlzi
        tdx(i) = tdxi
        tdy(i) = tdyi
        tdz(i) = tdzi
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
          RijSquared = RXij**2 + RYij**2 + RZij**2
          if( RijSquared >= RCutoffSquared ) cycle loop3
          RijSquaredInv = SigmaSquared / RijSquared
          Rij6Inv = RijSquaredInv**3
          EPotLocal = EPotLocal + (Rij6Inv * (Rij6Inv - 1._RK))
          EPotLocalInter = EPotLocalInter + (Rij6Inv * (Rij6Inv - 1._RK))
          Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
          FXij = Fij * RXij
          FYij = Fij * RYij
          FZij = Fij * RZij
          VirialLocal = VirialLocal + (PXij * FXij + PYij * FYij + PZij * FZij)
          VirialLocalInter = VirialLocalInter + (PXij * FXij + PYij * FYij + PZij * FZij)
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)/RijSquared
          d2EpotdV2Local = d2EpotdV2Local + Epsilon4 * Rij6Inv * (12._RK*Rij6Inv  -  6._RK) * (sitecorr * sitecorr - Plen2/RijSquared)*Third*Third  !xxxx LJ SS T
          d2EpotdV2Local = d2EpotdV2Local + Epsilon4 * Rij6Inv * (156._RK*Rij6Inv - 42._RK) *  sitecorr * sitecorr*Third*Third
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
   EPot = EPot + this%Epsilon4 * EPotLocal
   Virial = Virial + Third * VirialLocal * BoxLength
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:) * BoxLength
#endif
   EPotInter = EPotInter + this%Epsilon4 * EPotLocalInter
   VirialInter = VirialInter + Third * VirialLocalInter * BoxLength
   if (IntraLJEl) then
     EPotIntra_Nonbonded = EPotIntra_Nonbonded + this%Epsilon4 * EPotLocalIntra
     VirialIntra = VirialIntra + Third * VirialLocalIntra * BoxLength
   end if
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

#ifdef ABL
    VirAblSig = VirAblSig * Third * BoxLength * 18._RK * Epsilon4 / this%Sigma
    VirAblEps = VirAblEps * Third * BoxLength * 24._RK
#endif
  end subroutine TPotLJLJ_Force_Trans


!==============================================================!
!  Subroutine TPotLJLJ_RDF                                     !
!==============================================================!

  subroutine TPotLJLJ_RDF( this, RDFdr )

    implicit none

    ! Declare arguments
    type(TPotLJ126LJ126)     :: this
    real(RK), intent(in)     :: RDFdr

    !RDF RDFdr und RDFSchalenIndex
    real(RK)          :: distance
    integer           :: RDFSchalenIndex

    ! Declare local variables
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: RXi, RYi, RZi
    integer           :: i, j, k
    integer           :: jk, unit

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

      unit=this%NUnit1*(i-1)+this%Site1%UnitNumber

!CDIR NODEP
loop1:do k = 1, this%NInCutoff(unit)
        j = this%CutoffPartner(k, unit) ! Unit-partner of this unit
        if ( mod(j-this%Site2%UnitNumber, this%NUnit2)==0) then  ! choose only units, to which our Site2 correspond
          if (mod(j,this%NUnit2)==0) then
            jk = INT(j/this%NUnit2)   ! number of molecule, to which this unit correspond
          else
            jk = INT(j/this%NUnit2)+1
          end if
          RXij = RXi - RX2(jk)
          RYij = RYi - RY2(jk)
          RZij = RZi - RZ2(jk)
          RXij = RXij - anint( RXij )
          RYij = RYij - anint( RYij )
          RZij = RZij - anint( RZij )

!RDF in Schalen sortieren
          distance = sqrt(RXij**2 + RYij**2 + RZij**2)
          RDFSchalenIndex = INT(distance/RDFdr) + 1
          if (RDFSchalenIndex .le. RDFNumberShells) then
            this%RDFSum(RDFSchalenIndex) = this%RDFSum(RDFSchalenIndex) + 1
          end if
        end if
      end do loop1
    end do

  end subroutine TPotLJLJ_RDF


!==============================================================!
!  Subroutine TPotLJLJ_ChemicalPotential                       !
!==============================================================!

  subroutine TPotLJLJ_ChemicalPotential( this, EPotTest )

    implicit none

    ! Declare arguments
    type(TPotLJ126LJ126) :: this
    real(RK), pointer, contiguous    :: EPotTest(:)

    ! Declare local variables
    real(RK)          :: SigmaSquared
    real(RK)          :: Epsilon4
    real(RK)          :: RCutoffSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer:: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: RijSquared, RijSquaredInv, Rij6Inv
    real(RK)          :: EPotLocal
    integer           :: N2
    integer           :: i, j, k
    integer           :: i0, i1, jk, unit

    ! Assign local variables
    N2 = this%Site2%NPart
    SigmaSquared = this%SigmaSquared
    Epsilon4 = this%Epsilon4
    RCutoffSquared = this%RCutoffSquaredScaled
#if MPI_VER > 0
    i0 = this%Site1%NTest0
    i1 = this%Site1%NTest2
#else
    i0 = 1
    i1 = this%Site1%NTest
#endif

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
!$OMP PRIVATE (RijSquared,RijSquaredInv,Rij6Inv) &
!$OMP PRIVATE (EpotLocal,i,i0,i1,j,k) 

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
!$OMP DO
      do i = i0, i1
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        EPotLocal = 0._RK

        unit = this%NUnit1*(i-1)+this%Site1%UnitNumber

!CDIR NODEP
loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, this%NUnit2)==0) then  ! choose only units, to which our Site2 correspond
            if (mod(j,this%NUnit2)==0) then
              jk = INT(j/this%NUnit2)   ! number of molecule, to which this unit correspond
            else
              jk = INT(j/this%NUnit2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = RXij - anint( PXij )
            RYij = RYij - anint( PYij )
            RZij = RZij - anint( PZij )
            RijSquared = RXij**2 + RYij**2 + RZij**2
            RijSquaredInv = SigmaSquared / RijSquared
            Rij6Inv = RijSquaredInv**3
            EPotLocal = EPotLocal + Rij6Inv * (Rij6Inv - 1._RK)
          end if
        end do loop1
        EPotTest(i) = EPotTest(i) + Epsilon4 * EPotLocal
      end do
!$OMP END DO
    else

      ! Loop over test particles
!$OMP DO
      do i = i0, i1
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
          RijSquared = RXij**2 + RYij**2 + RZij**2

          if( RijSquared >= RCutoffSquared ) cycle loop2
          RijSquaredInv = SigmaSquared / RijSquared
          Rij6Inv = RijSquaredInv**3
          EPotLocal = EPotLocal + Rij6Inv * (Rij6Inv - 1._RK)
        end do loop2
        EPotTest(i) = EPotTest(i) + Epsilon4 * EPotLocal
      end do
!$OMP END DO
    end if
!$OMP END PARALLEL
  end subroutine TPotLJLJ_ChemicalPotential



!==============================================================!
!  Subroutine TPotLJLJ_Energy                                  !
!==============================================================!

  subroutine TPotLJLJ_Energy( this, np, nu, F, E, EIntra, BoxLength, CompIdent )

    implicit none

    ! Declare arguments
    type(TPotLJ126LJ126) :: this
    integer, intent(in)  :: np
    integer, intent(in)      :: nu
    real(RK), intent(in out) :: F(3,nu)
    real(RK), intent(in out) :: E
    real(RK), intent(in out) :: EIntra
    real(RK), intent(in) :: BoxLength
    logical, intent(in)      :: CompIdent

    ! Declare local variables
    real(RK)          :: SigmaSquared
    real(RK)          :: Epsilon48
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij, Fij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: RijSquared, RijSquaredInv, Rij6Inv
    real(RK)          :: tempF(3,nu)
    real(RK)          :: EPot, EIntra1, ELocal
    integer           :: j, k
    integer           :: nu2, unit, jk
    real(RK)          :: coeff

    ! Assign local variables
    SigmaSquared = this%SigmaSquared
    Epsilon48 = this%Epsilon48
    nu2 = this%NUnit2
    coeff = 1._RK
    if (this%potintra14) coeff = this%ScaleLJ14
    EPot = 0._RK
    EIntra1 = 0._RK
    tempF(:,:) = 0._RK

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

    ! Loop over molecules
    RXi = RX1(np)
    RYi = RY1(np)
    RZi = RZ1(np)
    PXi = PX1(np)
    PYi = PY1(np)
    PZi = PZ1(np)

    unit=this%NUnit1*(np-1)+this%Site1%UnitNumber ! Number of unit, to which this site corresponds

    do k = 1, this%NInCutoff(unit)
      j = this%CutoffPartner(k, unit) ! Unit-partner of this unit
      if ( mod(j-this%Site2%UnitNumber, nu2)==0) then  ! choose only units, to which our Site2 correspond
        if (mod(j,nu2)==0) then
          jk = INT(j/nu2)   ! number of molecule, to which this unit correspond
        else
          jk = INT(j/nu2)+1
        end if
        RXij = RXi - RX2(jk)
        RYij = RYi - RY2(jk)
        RZij = RZi - RZ2(jk)
        PXij = PXi - PX2(jk)
        PYij = PYi - PY2(jk)
        PZij = PZi - PZ2(jk)
        RXij = RXij - anint( PXij )
        RYij = RYij - anint( PYij )
        RZij = RZij - anint( PZij )
        RijSquared = RXij**2 + RYij**2 + RZij**2
        RijSquaredInv = SigmaSquared / RijSquared
        Rij6Inv = RijSquaredInv**3
        EPot = EPot + Rij6Inv * (Rij6Inv - 1._RK)
        Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
        FXij = Fij * RXij
        FYij = Fij * RYij
        FZij = Fij * RZij
        if (CompIdent) then
          tempF(1,this%Site1%UnitNumber) = tempF(1,this%Site1%UnitNumber) + FXij
          tempF(2,this%Site1%UnitNumber) = tempF(2,this%Site1%UnitNumber) + FYij
          tempF(3,this%Site1%UnitNumber) = tempF(3,this%Site1%UnitNumber) + FZij
        else
          tempF(1,this%Site2%UnitNumber) = tempF(1,this%Site2%UnitNumber) - FXij
          tempF(2,this%Site2%UnitNumber) = tempF(2,this%Site2%UnitNumber) - FYij
          tempF(3,this%Site2%UnitNumber) = tempF(3,this%Site2%UnitNumber) - FZij
        end if
      end if
    end do
    ! Include intramolecular interaction if need
    if (this%potintra15 .or. this%potintra14) then
      RXij = RXi - RX2(np)
      RYij = RYi - RY2(np)
      RZij = RZi - RZ2(np)
      PXij = PXi - PX2(np)
      PYij = PYi - PY2(np)
      PZij = PZi - PZ2(np)
      RXij = RXij - anint( PXij )
      RYij = RYij - anint( PYij )
      RZij = RZij - anint( PZij )
      RijSquared = RXij**2 + RYij**2 + RZij**2
      RijSquaredInv = SigmaSquared / RijSquared
      Rij6Inv = RijSquaredInv**3
      EIntra1 = EIntra1 + (Rij6Inv * (Rij6Inv - 1._RK) * coeff)
      Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
      tempF(1,this%Site1%UnitNumber) = tempF(1,this%Site1%UnitNumber) + Fij * RXij
      tempF(2,this%Site1%UnitNumber) = tempF(2,this%Site1%UnitNumber) + Fij * RYij
      tempF(3,this%Site1%UnitNumber) = tempF(3,this%Site1%UnitNumber) + Fij * RZij
      tempF(1,this%Site2%UnitNumber) = tempF(1,this%Site2%UnitNumber) - Fij * RXij
      tempF(2,this%Site2%UnitNumber) = tempF(2,this%Site2%UnitNumber) - Fij * RYij
      tempF(3,this%Site2%UnitNumber) = tempF(3,this%Site2%UnitNumber) - Fij * RZij
    end if

    F(:,:) = F(:,:) + tempF(:,:)
    E = E + this%Epsilon4 * (EPot + EIntra1)
    EIntra = EIntra + this%Epsilon4 * EIntra1

  end subroutine TPotLJLJ_Energy



!==============================================================!
!  Subroutine TPotLJLJ_UpdateBoxLength                         !
!==============================================================!

  subroutine TPotLJLJ_UpdateBoxLength( this, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotLJ126LJ126) :: this
    real(RK), intent(in) :: BoxLength

    ! Declare local variables
    real(RK) :: BoxLengthInv

    ! Update constants
    BoxLengthInv = 1._RK / BoxLength
    this%BoxLengthInv = BoxLengthInv
    this%BoxLengthThird = Third * BoxLength
    this%SigmaSquared = (this%Sigma * BoxLengthInv)**2
    this%Epsilon48 = 12._RK * this%Epsilon4 * BoxLengthInv / this%SigmaSquared
    this%RCutoffSquaredScaled = this%RCutoffSquared * BoxLengthInv**2

  end subroutine TPotLJLJ_UpdateBoxLength



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

    ! Declare local variables
    integer :: k, ende

    ! Construct potential
    this%Site1 => Molecule1%SiteCharge(j1)
    this%Site2 => Molecule2%SiteCharge(j2)
    this%NUnit1 => Molecule1%NUnit
    this%NUnit2 => Molecule2%NUnit
    this%SameComponent = i1 == i2
    this%Epsilon = this%Site1%e * this%Site2%e
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2
    this%RFConstant = this%Epsilon / RCutoff**3 * (RFEpsilon - 1._RK) / (2._RK * RFEpsilon + 1._RK)

    ! if this potential is intra
    if (this%SameComponent .and. Molecule1%hasIntraLJEl) then
      ende = size(Molecule1%IntCC15(:,1))
      do k=1, ende
        if (Molecule1%IntCC15(k,1)==this%Site1%SiteId .and. Molecule1%IntCC15(k,2)==this%Site2%SiteId) then
          this%potintra15 = .true.
        end if
      end do
      if (LJEL14 .and. .not. this%potintra15) then
        ende = size(Molecule1%IntCC14(:,1))
        do k=1, ende
          if (Molecule1%IntCC14(k,1)==this%Site1%SiteId .and. Molecule1%IntCC14(k,2)==this%Site2%SiteId) then
            this%potintra14=.true.
            this%ScaleEl14 = Molecule1%ScaleCC14(k)
          end if
        end do
      end if
    else
      this%potintra15 = .false.
      this%potintra14 = .false.
    end if    

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

  subroutine TPotCC_Force( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra_Nonbonded
    real(RK), intent(in out) :: VirialIntra
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
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: Rij2
    integer           :: i, j, k, i1
    integer           :: jk, unit
    real(RK)          :: coeff
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
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

!$OMP PARALLEL &
#if MPI_VER > 0
!$OMP FIRSTPRIVATE (i0) &
#endif
!$OMP PRIVATE(i1) &
!$OMP PRIVATE (Plen2,sitecorr) &
!$OMP PRIVATE (RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE (RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE (eX, eY, eZ  , RijInv, EPotLocal1,  i, j, k)
    
    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

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

      unit=this%NUnit1*(i-1)+this%Site1%UnitNumber

loop1:do k = 1, this%NInCutoff(unit)
        j = this%CutoffPartner(k, unit)
        if ( mod(j-this%Site2%UnitNumber, this%NUnit2)==0) then
          if (mod(j,this%NUnit2)==0) then
            jk = INT(j/this%NUnit2)
          else
            jk = INT(j/this%NUnit2)+1
          end if
          RXij = RXi - RX2(jk)
          RYij = RYi - RY2(jk)
          RZij = RZi - RZ2(jk)
          PXij = PXi - PX2(jk)
          PYij = PYi - PY2(jk)
          PZij = PZi - PZ2(jk)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          Rij2   = RXij**2 + RYij**2 + RZij**2
#if ARCH == 3
          RijInv = rsqrt( Rij2 )
#else
          RijInv = 1._RK / sqrt( Rij2 )

#endif
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          EPotLocal1 = Epsilon * RijInv
          EPotLocal  = EPotLocal + EPotLocal1
          EPotLocalInter  = EPotLocalInter + EPotLocal1
          VirialLocal = VirialLocal + (EPotLocal1 * RijInv * (eX * PXij + eY * PYij + eZ * PZij))
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
          VirialLocalInter = VirialLocalInter + (EPotLocal1 * RijInv * (eX * PXij + eY * PYij + eZ * PZij))
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (RXij * PXij + RYij * PYij + RZij * PZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 * (3._RK*sitecorr*sitecorr - Plen2*RijInv*RijInv)*Third*Third !XXXX CC
          FXij = EPotLocal1 * RijInv * eX
          FYij = EPotLocal1 * RijInv * eY
          FZij = EPotLocal1 * RijInv * eZ
          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          forceTempX(jk) = forceTempX(jk) - FXij
          forceTempY(jk) = forceTempY(jk) - FYij
          forceTempZ(jk) = forceTempZ(jk) - FZij
        end if
      end do loop1
      ! Include intramolecular interaction if need
      if (this%potintra15 .or. this%potintra14) then
        RXij = RXi - RX2(i)
        RYij = RYi - RY2(i)
        RZij = RZi - RZ2(i)
        PXij = PXi - PX2(i)
        PYij = PYi - PY2(i)
        PZij = PZi - PZ2(i)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
#if ARCH == 3
        RijInv = rsqrt( RXij**2 + RYij**2 + RZij**2 )
#else
        RijInv = 1._RK / sqrt( RXij**2 + RYij**2 + RZij**2 )
#endif
        eX = RXij * RijInv
        eY = RYij * RijInv
        eZ = RZij * RijInv
        EPotLocal1 = Epsilon * RijInv*coeff
        EPotLocal  = EPotLocal + EPotLocal1
        EPotLocalIntra  = EPotLocalIntra + EPotLocal1
        VirialLocal = VirialLocal + (EPotLocal1 * RijInv * (eX * PXij + eY * PYij + eZ * PZij))
        VirialLocalIntra = VirialLocalIntra + (EPotLocal1 * RijInv * (eX * PXij + eY * PYij + eZ * PZij))
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (RXij * PXij + RYij * PYij + RZij * PZij)*RijInv*RijInv
        d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 * (3._RK*sitecorr*sitecorr - Plen2*RijInv*RijInv)*Third*Third !XXXX CC
        FXij = EPotLocal1 * RijInv * eX
        FYij = EPotLocal1 * RijInv * eY
        FZij = EPotLocal1 * RijInv * eZ
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        forceTempX(i) = forceTempX(i) - FXij
        forceTempY(i) = forceTempY(i) - FYij
        forceTempZ(i) = forceTempZ(i) - FZij

      end if
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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotCC_Force


!==============================================================!
!  Subroutine TPotCC_Force_Ewald                               !
!==============================================================!

  subroutine TPotCC_Force_Ewald( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength, Kappa )

    implicit none

    ! Declare arguments
    type(TPotChargeCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra_Nonbonded
    real(RK), intent(in out) :: VirialIntra
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
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: approx, Faktor
    real(RK)          :: Fij,KappaRij
    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: Rij2
    integer           :: i, j, k, i1, i2
    integer           :: nu1, nu2, jk, unit
    logical           :: SameComponent
    real(RK)          :: coeff
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
    SameComponent = this%SameComponent
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

!$OMP PARALLEL DEFAULT(SHARED) &
#if MPI_VER > 0
!$OMP FIRSTPRIVATE ( i0) &
#endif
!$OMP PRIVATE (i1, i2) &
!$OMP PRIVATE ( approx, Fij,KappaRij,Rij2) &
!$OMP PRIVATE ( RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi) &
!$OMP PRIVATE ( RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE ( eX, eY, eZ  , RijInv,Rij, EPotLocal1,  i, j, k)

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

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

      unit=nu1*(i-1)+this%Site1%UnitNumber

 loop1:do k = 1, this%NInCutoff(unit)
         j = this%CutoffPartner(k, unit)
        if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
          if (mod(j,nu2)==0) then
            jk = INT(j/nu2)
          else
            jk = INT(j/nu2)+1
          end if
          RXij = RXi - RX2(jk)
          RYij = RYi - RY2(jk)
          RZij = RZi - RZ2(jk)
          PXij = PXi - PX2(jk)
          PYij = PYi - PY2(jk)
          PZij = PZi - PZ2(jk)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          Rij2   = RXij**2 + RYij**2 + RZij**2
          Rij =  sqrt(Rij2)
#if ARCH == 3
          RijInv = 1._RK /  Rij 
#else
          RijInv = 1._RK /  Rij 
#endif
          KappaRij = Kappa*Rij
          call ErrorApprox(this, KappaRij,approx)

          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          EPotLocal1 = Epsilon * RijInv * approx
          EPotLocal  = EPotLocal + EPotLocal1
          EPotLocalInter  = EPotLocalInter + EPotLocal1
          Fij  = (EPotLocal1 + Faktor*exp(-KappaRij**2)*Epsilon) * RijInv
          VirialLocal = VirialLocal + (Fij * (eX * PXij + eY * PYij + eZ * PZij))
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
          VirialLocalInter = VirialLocalInter + (Fij * (eX * PXij + eY * PYij + eZ * PZij))
          FXij = Fij * eX
          FYij = Fij * eY
          FZij = Fij * eZ
          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          forceTempX(jk) = forceTempX(jk) - FXij
          forceTempY(jk) = forceTempY(jk) - FYij
          forceTempZ(jk) = forceTempZ(jk) - FZij
        end if

      end do loop1
      ! Include intramolecular interaction if need
      if (this%potintra15 .or. this%potintra14) then
        RXij = RXi - RX2(i)
        RYij = RYi - RY2(i)
        RZij = RZi - RZ2(i)
        PXij = PXi - PX2(i)
        PYij = PYi - PY2(i)
        PZij = PZi - PZ2(i)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        Rij =  sqrt(RXij**2 + RYij**2 + RZij**2)
        if (PXij**2+PYij**2+PZij**2 .eq. 0.0) then
          Rij = 1.e33
        end if

#if ARCH == 3
        RijInv = 1._RK /  Rij
#else
        RijInv = 1._RK /  Rij
#endif
        KappaRij = Kappa*Rij
        call ErrorApprox(this, KappaRij,approx)

        eX = RXij * RijInv
        eY = RYij * RijInv
        eZ = RZij * RijInv
        EPotLocal1 = Epsilon * RijInv * approx*coeff
        EPotLocal  = EPotLocal + EPotLocal1
        EPotLocalIntra  = EPotLocalIntra + EPotLocal1
        Fij  = (EPotLocal1 + coeff*Faktor*exp(-KappaRij**2)*Epsilon) * RijInv
        VirialLocal = VirialLocal + (Fij * (eX * PXij + eY * PYij + eZ * PZij))
        VirialLocalIntra = VirialLocalIntra + (Fij * (eX * PXij + eY * PYij + eZ * PZij))
        FXij = Fij * eX
        FYij = Fij * eY
        FZij = Fij * eZ
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        forceTempX(i) = forceTempX(i) - FXij
        forceTempY(i) = forceTempY(i) - FYij
        forceTempZ(i) = forceTempZ(i) - FZij
      end if

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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if

  end subroutine TPotCC_Force_Ewald


!==============================================================!
!  Subroutine TPotCC_Force_Trans                               !
!==============================================================!

  subroutine TPotCC_Force_Trans( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra_Nonbonded
    real(RK), intent(in out) :: VirialIntra
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
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    real(RK)          :: Rij2
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, jk, unit
    real(RK)          :: coeff
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
    real(RK), pointer, contiguous :: VSx(:), VSy(:), VSz(:)
    real(RK), pointer, contiguous :: VSux(:), VSuy(:), VSuz(:)
    real(RK), pointer, contiguous :: VBx(:), VBy(:), VBz(:)
    real(RK), pointer, contiguous :: Cx(:) , Cy(:) , Cz(:)
    real(RK), pointer, contiguous :: tux(:) , tuy(:) , tuz(:)
    real(RK), pointer, contiguous :: tlx(:) , tly(:) , tlz(:)
    real(RK), pointer, contiguous :: tdx(:) , tdy(:) , tdz(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VBxi, VByi, VBzi
    real(RK)          :: Cxi,  Cyi,  Czi
    real(RK)          :: tuxi,  tuyi,  tuzi
    real(RK)          :: tlxi,  tlyi,  tlzi
    real(RK)          :: tdxi,  tdyi,  tdzi
    real(RK)          :: txii,  tyii , tzii
    real(RK)          :: txi ,  tyi  , tzi
    real(RK)          :: UU, Uxi,  Uyi, Uzi
    real(RK)          :: r1x, r1y, r1z
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33
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
!$OMP PRIVATE(  FX1, FY1, FZ1 ) &
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
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

#if  TRANS == 1
    !TRANSPORT_start
    VSx => this%Site1%vsCx
    VSy => this%Site1%vsCy
    VSz => this%Site1%vsCz
    VBx => this%Site1%vbCx
    VBy => this%Site1%vbCy
    VBz => this%Site1%vbCz
    VSux=> this%Site1%vsuCx
    VSuy=> this%Site1%vsuCy
    VSuz=> this%Site1%vsuCz
    Cx  => this%Site1%cCx
    Cy  => this%Site1%cCy
    Cz  => this%Site1%cCz
    tux => this%Site1%tuCx
    tuy => this%Site1%tuCy
    tuz => this%Site1%tuCz
    tlx => this%Site1%tlCx
    tly => this%Site1%tlCy
    tlz => this%Site1%tlCz
    tdx => this%Site1%tdCx
    tdy => this%Site1%tdCy
    tdz => this%Site1%tdCz
    q1  => this%Site1%Qm0r(:, 1, 1)
    q2  => this%Site1%Qm0r(:, 2, 1)
    q3  => this%Site1%Qm0r(:, 3, 1)
    q4  => this%Site1%Qm0r(:, 4, 1)
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
      VSxi= VSx(i)
      VSyi= VSy(i)
      VSzi= VSz(i)
      VBxi= VBx(i)
      VByi= VBy(i)
      VBzi= VBz(i)
      VSuxi= VSux(i)
      VSuyi= VSuy(i)
      VSuzi= VSuz(i)
      Cxi = Cx(i)
      Cyi = Cy(i)
      Czi = Cz(i)
      tuxi = tux(i)
      tuyi = tuy(i)
      tuzi = tuz(i)
      tlxi = tlx(i)
      tlyi = tly(i)
      tlzi = tlz(i)
      tdxi = tdx(i)
      tdyi = tdy(i)
      tdzi = tdz(i)
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

      unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:do k = 1, this%NInCutoff(unit)
        j = this%CutoffPartner(k, unit)
        if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
          if (mod(j,nu2)==0) then
            jk = INT(j/nu2)
          else
            jk = INT(j/nu2)+1
          end if
          RXij = RXi - RX2(jk)
          RYij = RYi - RY2(jk)
          RZij = RZi - RZ2(jk)
          PXij = PXi - PX2(jk)
          PYij = PYi - PY2(jk)
          PZij = PZi - PZ2(jk)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          Rij2   = RXij**2 + RYij**2 + RZij**2
#if ARCH == 3
          RijInv = rsqrt( Rij2 )
#else
          RijInv = 1._RK / sqrt( Rij2 )
#endif
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          EPotLocal1 = Epsilon * RijInv
          EPotLocal  = EPotLocal + EPotLocal1
          EPotLocalInter  = EPotLocalInter + EPotLocal1
          VirialLocal = VirialLocal + (EPotLocal1 * RijInv * (eX * PXij + eY * PYij + eZ * PZij))
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
        VirialLocalInter = VirialLocalInter + (EPotLocal1 * RijInv * (eX * PXij + eY * PYij + eZ * PZij))
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (RXij * PXij + RYij * PYij + RZij * PZij)*RijInv*RijInv
        d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 * (3._RK*sitecorr*sitecorr - Plen2*RijInv*RijInv)*Third*Third !xxxx CC T
        FXij = EPotLocal1 * RijInv * eX
        FYij = EPotLocal1 * RijInv * eY
        FZij = EPotLocal1 * RijInv * eZ
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        forceTempX(jk) = forceTempX(jk) - FXij
        forceTempY(jk) = forceTempY(jk) - FYij
        forceTempZ(jk) = forceTempZ(jk) - FZij

#if TRANS==1
        !TRANSPORT_start vielleicht
        VSxi   = VSxi + FXij * PYij
        VSyi   = VSyi + FXij * PZij
        VSzi   = VSzi + FYij * PZij
        VBxi   = VBxi + FXij * PXij
        VByi   = VByi + FYij * PYij
        VBzi   = VBzi + FZij * PZij
        VSuxi  = VSuxi+ FYij * PXij
        VSuyi  = VSuyi+ FZij * PXij
        VSuzi  = VSuzi+ FZij * PYij
        UU        = EpotLocal1 + this%RFConstant * Rij2
        Cxi    = Cxi  + UU * eX
        Cyi    = Cyi  + UU * eY
        Czi    = Czi  + UU * eZ
        txii   = r1y * FZij - r1z * FYij
        tyii   = r1z * FXij - r1x * FZij
        tzii   = r1x * FYij - r1y * FXij
        txi    = A11 * txii + A12 * tyii + A13 * tzii
        tyi    = A21 * txii + A22 * tyii + A23 * tzii
        tzi    = A31 * txii + A32 * tyii + A33 * tzii
        tuxi   = tuxi + PXij*tyi
        tuyi   = tuyi + PXij*tzi
        tuzi   = tuzi + PYij*tzi
        tlxi   = tlxi + PYij*txi
        tlyi   = tlyi + PZij*txi
        tlzi   = tlzi + PZij*tyi
        tdxi   = tdxi + PXij*txi
        tdyi   = tdyi + PYij*tyi
        tdzi   = tdzi + PZij*tzi
#endif

        end if
      end do loop1
      ! Include intramolecular interaction if need
      if (this%potintra15 .or. this%potintra14) then
        RXij = RXi - RX2(i)
        RYij = RYi - RY2(i)
        RZij = RZi - RZ2(i)
        PXij = PXi - PX2(i)
        PYij = PYi - PY2(i)
        PZij = PZi - PZ2(i)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
#if ARCH == 3
        RijInv = rsqrt( RXij**2 + RYij**2 + RZij**2 )
#else
        RijInv = 1._RK / sqrt( RXij**2 + RYij**2 + RZij**2 )
#endif
        eX = RXij * RijInv
        eY = RYij * RijInv
        eZ = RZij * RijInv
        EPotLocal1 = Epsilon * RijInv*coeff
        EPotLocal  = EPotLocal + EPotLocal1
        EPotLocalIntra  = EPotLocalIntra + EPotLocal1
        VirialLocal = VirialLocal + (EPotLocal1 * RijInv * (eX * PXij + eY * PYij + eZ * PZij))
        VirialLocalIntra = VirialLocalIntra + (EPotLocal1 * RijInv * (eX * PXij + eY * PYij + eZ * PZij))
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (RXij * PXij + RYij * PYij + RZij * PZij)*RijInv*RijInv
        d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 * (3._RK*sitecorr*sitecorr - Plen2*RijInv*RijInv)*Third*Third !xxxx CC T
        FXij = EPotLocal1 * RijInv * eX
        FYij = EPotLocal1 * RijInv * eY
        FZij = EPotLocal1 * RijInv * eZ
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        forceTempX(i) = forceTempX(i) - FXij
        forceTempY(i) = forceTempY(i) - FYij
        forceTempZ(i) = forceTempZ(i) - FZij


      end if

      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
#if  TRANS == 1
      !TRANSPORT_start
      VSx(i) = VSxi
      VSy(i) = VSyi
      VSz(i) = VSzi
      VBx(i) = VBxi
      VBy(i) = VByi
      VBz(i) = VBzi
      VSux(i)= VSuxi
      VSuy(i)= VSuyi
      VSuz(i)= VSuzi
      Cx(i)  = Cxi
      Cy(i)  = Cyi
      Cz(i)  = Czi
      tux(i) = tuxi
      tuy(i) = tuyi
      tuz(i) = tuzi
      tlx(i) = tlxi
      tly(i) = tlyi
      tlz(i) = tlzi
      tdx(i) = tdxi
      tdy(i) = tdyi
      tdz(i) = tdzi
      !TRANSPORT_END

#endif
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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotCC_Force_Trans



!==============================================================!
!  Subroutine TPotCC_Force_Ewald_Trans                         !
!==============================================================!

  subroutine TPotCC_Force_Ewald_Trans( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength, Kappa )

    implicit none

    ! Declare arguments
    type(TPotChargeCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra_Nonbonded
    real(RK), intent(in out) :: VirialIntra
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
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: EPotLocalInter, VirialLocalInter
    !real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    real(RK)          :: approx, Faktor
    real(RK)          :: Fij,KappaRij
    real(RK)          :: Rij2
    integer           :: i, j, k, i1, i2
    integer           :: nu1, nu2, jk, unit
    logical           :: SameComponent
    real(RK)          :: coeff
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
    real(RK), pointer, contiguous :: VSx(:), VSy(:), VSz(:)
    real(RK), pointer, contiguous :: VSux(:), VSuy(:), VSuz(:)
    real(RK), pointer, contiguous :: VBx(:), VBy(:), VBz(:)
    real(RK), pointer, contiguous :: Cx(:) , Cy(:) , Cz(:)
    real(RK), pointer, contiguous :: tux(:) , tuy(:) , tuz(:)
    real(RK), pointer, contiguous :: tlx(:) , tly(:) , tlz(:)
    real(RK), pointer, contiguous :: tdx(:) , tdy(:) , tdz(:)

    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VBxi, VByi, VBzi

    real(RK)          :: r1x, r1y, r1z
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33

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
!$OMP PRIVATE(  FX1, FY1, FZ1) &
!$OMP PRIVATE( PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE(   RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE(   RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
#if  TRANS == 1
!$OMP PRIVATE(VSx, VSy, VSz ,VSux,VSuy,VSuz, VBx, VBy, VBz, Cx , Cy , Cz) &
!$OMP PRIVATE( tux , tuy , tuz, tlx , tly , tlz, tdx , tdy , tdz) &
!$OMP PRIVATE( VSxi, VSyi, VSzi, VSuxi,VSuyi,VSuzi) &
!$OMP PRIVATE( VBxi, VByi, VBzi,  r1x, r1y, r1z) &
!$OMP PRIVATE( A11, A12, A13, A21, A22, A23, A31, A32, A33) &
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
    SameComponent = this%SameComponent
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

#if  TRANS == 1
    !TRANSPORT_start

    VSx => this%Site1%vsCx
    VSy => this%Site1%vsCy
    VSz => this%Site1%vsCz
    VBx => this%Site1%vbCx
    VBy => this%Site1%vbCy
    VBz => this%Site1%vbCz

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
      VSxi= VSx(i)
      VSyi= VSy(i)
      VSzi= VSz(i)
      VBxi= VBx(i)
      VByi= VBy(i)
      VBzi= VBz(i)
!       !TRANSPORT_END
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

      unit=nu1*(i-1)+this%Site1%UnitNumber

 loop1:do k = 1, this%NInCutoff(unit)
         j = this%CutoffPartner(k, unit)
        if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
          if (mod(j,nu2)==0) then
            jk = INT(j/nu2)
          else
            jk = INT(j/nu2)+1
          end if
          RXij = RXi - RX2(jk)
          RYij = RYi - RY2(jk)
          RZij = RZi - RZ2(jk)
          PXij = PXi - PX2(jk)
          PYij = PYi - PY2(jk)
          PZij = PZi - PZ2(jk)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          Rij2   = RXij**2 + RYij**2 + RZij**2
          Rij =  sqrt(Rij2)

#if ARCH == 3
          RijInv = 1._RK /  Rij 
#else
          RijInv = 1._RK /  Rij 
#endif
          KappaRij = Kappa*Rij
          call ErrorApprox(this, KappaRij,approx)
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          EPotLocal1 = Epsilon * RijInv * approx
          EPotLocal  = EPotLocal + EPotLocal1
          EPotLocalInter  = EPotLocalInter + EPotLocal1
          Fij  = (EPotLocal1 + Faktor*exp(-KappaRij**2)*Epsilon) * RijInv
          VirialLocal = VirialLocal + (Fij * (eX * PXij + eY * PYij + eZ * PZij))
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
        VirialLocalInter = VirialLocalInter + (Fij * (eX * PXij + eY * PYij + eZ * PZij))
        FXij = Fij * eX
        FYij = Fij * eY
        FZij = Fij * eZ
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij

        forceTempX(jk) = forceTempX(jk) - FXij
        forceTempY(jk) = forceTempY(jk) - FYij
        forceTempZ(jk) = forceTempZ(jk) - FZij

#if TRANS==1
        !TRANSPORT_start vielleicht
        VSxi   = VSxi + FXij * PYij
        VSyi   = VSyi + FXij * PZij
        VSzi   = VSzi + FYij * PZij
        VBxi   = VBxi + FXij * PXij
        VByi   = VByi + FYij * PYij
        VBzi   = VBzi + FZij * PZij

#endif
        end if
      end do loop1
      ! Include intramolecular interaction if need
      if (this%potintra15 .or. this%potintra14) then
        RXij = RXi - RX2(i)
        RYij = RYi - RY2(i)
        RZij = RZi - RZ2(i)
        PXij = PXi - PX2(i)
        PYij = PYi - PY2(i)
        PZij = PZi - PZ2(i)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        Rij =  sqrt(RXij**2 + RYij**2 + RZij**2)
        if (PXij**2+PYij**2+PZij**2 .eq. 0.0) then
          Rij = 1.e33
        end if

#if ARCH == 3
        RijInv = 1._RK /  Rij
#else
        RijInv = 1._RK /  Rij
#endif
        KappaRij = Kappa*Rij
        call ErrorApprox(this, KappaRij,approx)

        eX = RXij * RijInv
        eY = RYij * RijInv
        eZ = RZij * RijInv
        EPotLocal1 = Epsilon * RijInv * approx*coeff
        EPotLocal  = EPotLocal + EPotLocal1
        EPotLocalIntra  = EPotLocalIntra + EPotLocal1
        Fij  = (EPotLocal1 + coeff*Faktor*exp(-KappaRij**2)*Epsilon) * RijInv
        VirialLocal = VirialLocal + (Fij * (eX * PXij + eY * PYij + eZ * PZij))
        VirialLocalIntra = VirialLocalIntra + (Fij * (eX * PXij + eY * PYij + eZ * PZij))
        FXij = Fij * eX
        FYij = Fij * eY
        FZij = Fij * eZ
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        forceTempX(i) = forceTempX(i) - FXij
        forceTempY(i) = forceTempY(i) - FYij
        forceTempZ(i) = forceTempZ(i) - FZij

      end if
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
#if  TRANS == 1
      !TRANSPORT_start
      VSx(i) = VSxi
      VSy(i) = VSyi
      VSz(i) = VSzi
      VBx(i) = VBxi
      VBy(i) = VByi
      VBz(i) = VBzi

      !TRANSPORT_END
#endif
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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if

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
    real(RK), pointer:: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: RijInv, RijSquared
    real(RK)          :: EPotLocal
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, unit, i0, jk
#if ARCH == 3
    logical           :: hit
#endif

    ! Assign local variables
    i1 = this%Site1%NTest
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared
    nu1 = this%NUnit1
    nu2 = this%NUnit2
#if MPI_VER > 0
    i0 = this%Site1%NTest0
    i1 = this%Site1%NTest2
#else
    i0 = 1
    i1 = this%Site1%NTest
#endif

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
   do i = i0, i1
     RXi = RX1(i)
     RYi = RY1(i)
     RZi = RZ1(i)
     PXi = PX1(i)
     PYi = PY1(i)
     PZi = PZ1(i)
     EPotLocal = 0._RK

#if ARCH == 3
        hit = .false.
#endif

     unit = nu1*(i-1)+this%Site1%UnitNumber

!CDIR NODEP
loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then  ! choose only units, to which our Site2 correspond
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)   ! number of molecule, to which this unit correspond
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            PXij = (PXij - anint( PXij )) * BoxLength
            PYij = (PYij - anint( PYij )) * BoxLength
            PZij = (PZij - anint( PZij )) * BoxLength
            RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
            if( RijSquared <= RShieldSquared ) hit = .true.
#else
            if( RijSquared <= RShieldSquared ) then
              EPotLocal = 1E33_RK
              exit loop1
            end if
#endif

#if ARCH == 3
            RijInv = rsqrt( RijSquared )
#else
            RijInv = 1._RK / sqrt( RijSquared )
#endif
            EPotLocal = EPotLocal + Epsilon * RijInv
          end if
        end do loop1

#if ARCH == 3
        if( .not. hit ) then
          EPotTest(i) = EPotTest(i) + EPotLocal
        else
          EPotTest(i) = EPotTest(i) + 1E33_RK
        endif
#else
        EPotTest(i) = EPotTest(i) + EPotLocal
#endif
   end do

  end subroutine TPotCC_ChemicalPotential


!==============================================================!
!  Subroutine TPotCC_Energy                                    !
!==============================================================!

  subroutine TPotCC_Energy( this, np, nu, F,E, EIntra, BoxLength, CompIdent )

    implicit none

    ! Declare arguments
    type(TPotChargeCharge) :: this
    integer, intent(in)    :: np
    integer, intent(in)      :: nu
    real(RK), intent(in out) :: F(3,nu)
    real(RK), intent(in out) :: E
    real(RK), intent(in out) :: EIntra
    real(RK), intent(in)   :: BoxLength
    logical, intent(in)      :: CompIdent

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
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijInv, RijSquared
    real(RK)          :: EPot, EIntra1, EPotLocal, tempF(3,nu)
    integer           :: j, k
    integer           :: nu2, unit, su, jk
    real(RK)          :: coeff

    ! Assign local variables
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

    nu2 = this%NUnit2
    coeff = 1._RK
    if (this%potintra14) coeff = this%ScaleEl14
    EPot = 0._RK
    EIntra1 = 0._RK
    tempF(:,:) = 0._RK
    su = this%Site2%UnitNumber
    if (CompIdent) su = this%Site1%UnitNumber

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

    ! Loop over molecules
    RXi = RX1(np)
    RYi = RY1(np)
    RZi = RZ1(np)
    PXi = PX1(np)
    PYi = PY1(np)
    PZi = PZ1(np)

    unit=this%NUnit1*(np-1)+this%Site1%UnitNumber

    do k = 1, this%NInCutoff(unit)
      j = this%CutoffPartner(k, unit)
      if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
        if (mod(j,nu2)==0) then
          jk = INT(j/nu2)
        else
          jk = INT(j/nu2)+1
        end if
        RXij = RXi - RX2(jk)
        RYij = RYi - RY2(jk)
        RZij = RZi - RZ2(jk)
        PXij = PXi - PX2(jk)
        PYij = PYi - PY2(jk)
        PZij = PZi - PZ2(jk)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        RijSquared = RXij**2 + RYij**2 + RZij**2

        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
          tempF(:,su) = 1E33_RK
        else
#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          EPotLocal = Epsilon * RijInv
        end if
        EPot = EPot + EPotLocal
      end if
    end do
    ! Include intramolecular interaction if need
    if (this%potintra14 .or. this%potintra15) then
      RXij = RXi - RX2(np)
      RYij = RYi - RY2(np)
      RZij = RZi - RZ2(np)
      PXij = PXi - PX2(np)
      PYij = PYi - PY2(np)
      PZij = PZi - PZ2(np)
      RXij = (RXij - anint( PXij )) * BoxLength
      RYij = (RYij - anint( PYij )) * BoxLength
      RZij = (RZij - anint( PZij )) * BoxLength
#if ARCH == 3
      RijInv = rsqrt( RXij**2 + RYij**2 + RZij**2 )
#else
      RijInv = 1._RK / sqrt( RXij**2 + RYij**2 + RZij**2 )
#endif
      EPotLocal = Epsilon * RijInv * coeff
      eX = RXij * RijInv
      eY = RYij * RijInv
      eZ = RZij * RijInv

      EIntra1  = EIntra1 + EPotLocal
    end if

    ! Update potential energy and virial
    F(:,:) = F(:,:) + tempF(:,:)
    E = E + EPot + EIntra1
    EIntra = EIntra + EIntra1

  end subroutine TPotCC_Energy



!==============================================================!
!  Subroutine TPotCC_Energy_Ewald                              !
!==============================================================!

  subroutine TPotCC_Energy_Ewald( this, np, nu, F, E, EIntra, BoxLength, Kappa, CompIdent )

    implicit none

    ! Declare arguments
    type(TPotChargeCharge) :: this
    integer, intent(in)    :: np
    integer, intent(in)      :: nu
    real(RK), intent(in out) :: F(3,nu)
    real(RK), intent(in out) :: E
    real(RK), intent(in out) :: EIntra
    real(RK), intent(in)   :: BoxLength
    real(RK), intent(in)   :: Kappa
    logical, intent(in)      :: CompIdent

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
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijInv, RijSquared
    real(RK)          :: Epot, EIntra1, EPotLocal, coeff
    real(RK)          :: Fij, Faktor, tempF(3,nu)
    integer           :: j, k
    real(RK)          :: approx
    integer           :: nu2, unit, su, jk
    real(RK)          :: Rij, KappaRij

    ! Assign local variables
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared
    nu2 = this%NUnit2
    coeff = 1._RK
    if (this%potintra14) coeff = this%ScaleEl14
    Epot = 0._RK
    EIntra1 = 0._RK
    Faktor = 2._RK/sqrt(Pi) * Kappa
    su = this%Site2%UnitNumber
    if (CompIdent) su = this%Site1%UnitNumber

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

    ! Loop over molecules
    RXi = RX1(np)
    RYi = RY1(np)
    RZi = RZ1(np)
    PXi = PX1(np)
    PYi = PY1(np)
    PZi = PZ1(np)

    unit=this%NUnit1*(np-1)+this%Site1%UnitNumber

    do k = 1, this%NInCutoff(unit)
      j = this%CutoffPartner(k, unit)
      if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
        if (mod(j,nu2)==0) then
          jk = INT(j/nu2)
        else
          jk = INT(j/nu2)+1
        end if
        RXij = RXi - RX2(jk)
        RYij = RYi - RY2(jk)
        RZij = RZi - RZ2(jk)
        PXij = PXi - PX2(jk)
        PYij = PYi - PY2(jk)
        PZij = PZi - PZ2(jk)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        RijSquared = RXij**2 + RYij**2 + RZij**2

        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
          tempF(:,su) = 1E33_RK
        else
        Rij =  sqrt(RijSquared)
#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
          KappaRij = Kappa*Rij
          call ErrorApprox(this, KappaRij,approx)

          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          EPotLocal = Epsilon * RijInv * approx
          Fij  = (EPotLocal + Faktor*exp(-KappaRij**2)*Epsilon) * RijInv
        end if
        Epot = Epot + EPotLocal
      end if
    end do
    ! Include intramolecular interaction if need
    if (this%potintra15 .or. this%potintra14) then
      RXij = RXi - RX2(np)
      RYij = RYi - RY2(np)
      RZij = RZi - RZ2(np)
      PXij = PXi - PX2(np)
      PYij = PYi - PY2(np)
      PZij = PZi - PZ2(np)
      RXij = (RXij - anint( PXij )) * BoxLength
      RYij = (RYij - anint( PYij )) * BoxLength
      RZij = (RZij - anint( PZij )) * BoxLength
#if ARCH == 3
      RijInv = rsqrt( RXij**2 + RYij**2 + RZij**2 )
#else
      RijInv = 1._RK / sqrt( RXij**2 + RYij**2 + RZij**2 )
#endif
      Rij =  sqrt(RijSquared)
      KappaRij = Kappa*Rij
      call ErrorApprox(this, KappaRij,approx)
      EPotLocal = Epsilon * RijInv * approx*coeff
      eX = RXij * RijInv
      eY = RYij * RijInv
      eZ = RZij * RijInv
      Fij  = (EPotLocal + Faktor*exp(-KappaRij**2)*Epsilon) * RijInv

      EIntra1  = EIntra1 + EPotLocal
    end if

    F(:,:) = F(:,:) + tempF(:,:)
    E = E + Epot + EIntra1
    EIntra = EIntra + EIntra1

  end subroutine TPotCC_Energy_Ewald



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

    ! Declare local variables
    integer :: k, ende

    ! Construct potential
    this%Site1 => Molecule1%SiteCharge(j1)
    this%Site2 => Molecule2%SiteDipole(j2)
    this%NUnit1 => Molecule1%NUnit
    this%NUnit2 => Molecule2%NUnit
    this%SameComponent = i1 == i2
    this%Epsilon = this%Site1%e * this%Site2%D
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

    ! if this potential is intra
    if (this%SameComponent .and. Molecule1%hasIntraLJEl) then
      ende = size(Molecule1%IntCD15(:,1))
      do k=1, ende
        if (Molecule1%IntCD15(k,1)==this%Site1%SiteId .and. Molecule1%IntCD15(k,2)==this%Site2%SiteId) then
           this%potintra15 = .true.
        end if
      end do
      if (LJEL14 .and. .not. this%potintra15) then
        ende = size(Molecule1%IntCD14(:,1))
        do k=1, ende
          if (Molecule1%IntCD14(k,1)==this%Site1%SiteId .and. Molecule1%IntCD14(k,2)==this%Site2%SiteId) then
            this%potintra14=.true.
            this%ScaleEl14 = Molecule1%ScaleCD14(k)
          end if
        end do
      end if
    else
      this%potintra15 = .false.
      this%potintra14 = .false.
    end if

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

  subroutine TPotCD_Force( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeDipole)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra_Nonbonded
    real(RK), intent(in out) :: VirialIntra
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
    real(RK)          :: EPotLocalInter, ViriallocalInter
    real(RK)          :: EPotLocalIntra, ViriallocalIntra
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, jk, unit
    logical           :: SameComponent
    real(RK)          :: coeff
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
    SameComponent = this%SameComponent
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._Rk
    end if

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

      unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:do k = 1, this%NInCutoff(unit)
        j = this%CutoffPartner(k, unit)
        if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
          if (mod(j,nu2)==0) then
            jk = INT(j/nu2)
          else
            jk = INT(j/nu2)+1
          end if
          RXij = RXi - RX2(jk)
          RYij = RYi - RY2(jk)
          RZij = RZi - RZ2(jk)
          PXij = PXi - PX2(jk)
          PYij = PYi - PY2(jk)
          PZij = PZi - PZ2(jk)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          OXj = OX2(jk)
          OYj = OY2(jk)
          OZj = OZ2(jk)
          RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
          RijInv = sqrt( RijSquaredInv )
          eX = RXij * RijInv                                                      ! Einheitsabstandvektor nach Price
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosTheta  = OXj * ex + OYj * eY + OZj * eZ                              ! cos(alpha) nach Price
          CosTheta3 = 3._RK * CosTheta
          Epsilon1 = Epsilon * RijSquaredInv
          Epsilon2 = Epsilon1 * RijInv
          EPotLocal  = EPotLocal + Epsilon1 * CosTheta                            ! Uebereinstimmumg mit Price
          EPotLocalInter  = EPotLocalInter + Epsilon1 * CosTheta
          FXij = Epsilon2 * ( CosTheta3 * eX - OXj )                              ! F2 bei Price
          FYij = Epsilon2 * ( CosTheta3 * eY - OYj )
          FZij = Epsilon2 * ( CosTheta3 * eZ - OZj )
          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F2*R_COM_Price; stimmt so
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
          VirialLocalInter = VirialLocalInter + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F2*R_COM_Price; stimmt so
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
          d2EpotdV2Local = d2EpotdV2Local + Epsilon1 * CosTheta*(8._RK*sitecorr*sitecorr-2._RK*Plen2*RijSquaredInv)*Third*Third   !xxxx2 CD
          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(jk) = forceTempX(jk) - FXij
          forceTempY(jk) = forceTempY(jk) - FYij
          forceTempZ(jk) = forceTempZ(jk) - FZij
          momTempX(jk) = momTempX(jk) - Epsilon1 * eX   
          momTempY(jk) = momTempY(jk) - Epsilon1 * eY   
          momTempZ(jk) = momTempZ(jk) - Epsilon1 * eZ   

        end if

      end do loop1
      ! Include intramolecular interactions if need
      if (this%potintra15 .or. this%potintra14) then
        RXij = RXi - RX2(i)
        RYij = RYi - RY2(i)
        RZij = RZi - RZ2(i)
        PXij = PXi - PX2(i)
        PYij = PYi - PY2(i)
        PZij = PZi - PZ2(i)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        OXj = OX2(i)
        OYj = OY2(i)
        OZj = OZ2(i)
        RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
        RijInv = sqrt( RijSquaredInv )
        eX = RXij * RijInv                                                      ! Einheitsabstandvektor nach Price
        eY = RYij * RijInv
        eZ = RZij * RijInv
        CosTheta  = OXj * ex + OYj * eY + OZj * eZ                              ! cos(alpha) nach Price
        CosTheta3 = 3._RK * CosTheta
        Epsilon1 = Epsilon * RijSquaredInv*coeff                                ! 1-4 non-bonded interaction coeff
        Epsilon2 = Epsilon1 * RijInv
        EPotLocal  = EPotLocal + Epsilon1 * CosTheta                           ! Uebereinstimmumg mit Price
        EPotLocalIntra  = EPotLocalIntra + Epsilon1 * CosTheta
        FXij = Epsilon2 * ( CosTheta3 * eX - OXj )                              ! F2 bei Price
        FYij = Epsilon2 * ( CosTheta3 * eY - OYj )
        FZij = Epsilon2 * ( CosTheta3 * eZ - OZj )
        VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F2*R_COM_Price; stimmt so
        VirialLocalIntra = VirialLocalIntra + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F2*R_COM_Price; stimmt so
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + Epsilon1 * CosTheta*(8._RK*sitecorr*sitecorr-2._RK*Plen2*RijSquaredInv)*Third*Third   !xxxx2 CD
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij

        forceTempX(i) = forceTempX(i) - FXij
        forceTempY(i) = forceTempY(i) - FYij
        forceTempZ(i) = forceTempZ(i) - FZij
        momTempX(i) = momTempX(i) - Epsilon1 * eX   
        momTempY(i) = momTempY(i) - Epsilon1 * eY   
        momTempZ(i) = momTempZ(i) - Epsilon1 * eZ   

      end if
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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
       EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
       VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotCD_Force

!==============================================================!
!  Subroutine TPotCD_Force_Trans                               !
!==============================================================!

  subroutine TPotCD_Force_Trans( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeDipole)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra_Nonbonded
    real(RK), intent(in out) :: VirialIntra
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
    real(RK)          :: EPotLocalInter, ViriallocalInter
    real(RK)          :: EPotLocalIntra, ViriallocalIntra
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, jk, unit
    logical           :: SameComponent
    real(RK)          :: coeff
    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: momTempX(1:this%Site2%NPart)
    real(RK)          :: momTempY(1:this%Site2%NPart)
    real(RK)          :: momTempZ(1:this%Site2%NPart)

#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer, contiguous :: VSx(:), VSy(:), VSz(:)
    real(RK), pointer, contiguous :: VSux(:), VSuy(:), VSuz(:)
    real(RK), pointer, contiguous :: VBx(:), VBy(:), VBz(:)
    real(RK), pointer, contiguous :: Cx(:) , Cy(:) , Cz(:)
    real(RK), pointer, contiguous :: tux(:) , tuy(:) , tuz(:)
    real(RK), pointer, contiguous :: tlx(:) , tly(:) , tlz(:)
    real(RK), pointer, contiguous :: tdx(:) , tdy(:) , tdz(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VBxi, VByi, VBzi
    real(RK)          :: Cxi,  Cyi,  Czi
    real(RK)          :: tuxi,  tuyi,  tuzi
    real(RK)          :: tlxi,  tlyi,  tlzi
    real(RK)          :: tdxi,  tdyi,  tdzi
    real(RK)          :: txii,  tyii , tzii
    real(RK)          :: txi ,  tyi  , tzi
    real(RK)          :: UU, Uxi,  Uyi, Uzi
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

!$OMP PARALLEL &
!$OMP PRIVATE ( Epsilon, Epsilon1, Epsilon2, RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (  FX1, FY1, FZ1, OX2, OY2, OZ2) &
!$OMP PRIVATE ( Plen2,PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (   RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE (   RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE ( OXj, OYj, OZj, eX, eY, eZ, RijSquaredInv, RijInv) &
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
    SameComponent = this%SameComponent
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

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
    VSx => this%Site1%vsCx
    VSy => this%Site1%vsCy
    VSz => this%Site1%vsCz
    VBx => this%Site1%vbCx
    VBy => this%Site1%vbCy
    VBz => this%Site1%vbCz
    VSux=> this%Site1%vsuCx
    VSuy=> this%Site1%vsuCy
    VSuz=> this%Site1%vsuCz
    Cx  => this%Site1%cCx
    Cy  => this%Site1%cCy
    Cz  => this%Site1%cCz
    tux => this%Site1%tuCx
    tuy => this%Site1%tuCy
    tuz => this%Site1%tuCz
    tlx => this%Site1%tlCx
    tly => this%Site1%tlCy
    tlz => this%Site1%tlCz
    tdx => this%Site1%tdCx
    tdy => this%Site1%tdCy
    tdz => this%Site1%tdCz
    q1  => this%Site1%Qm0r(:, 1, 1)
    q2  => this%Site1%Qm0r(:, 2, 1)
    q3  => this%Site1%Qm0r(:, 3, 1)
    q4  => this%Site1%Qm0r(:, 4, 1)
!TRANSPORT_END
#endif

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._Rk
    end if

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
      VSxi= VSx(i)
      VSyi= VSy(i)
      VSzi= VSz(i)
      VBxi= VBx(i)
      VByi= VBy(i)
      VBzi= VBz(i)
      VSuxi= VSux(i)
      VSuyi= VSuy(i)
      VSuzi= VSuz(i)
      Cxi = Cx(i)
      Cyi = Cy(i)
      Czi = Cz(i)
      tuxi = tux(i)
      tuyi = tuy(i)
      tuzi = tuz(i)
      tlxi = tlx(i)
      tlyi = tly(i)
      tlzi = tlz(i)
      tdxi = tdx(i)
      tdyi = tdy(i)
      tdzi = tdz(i)
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

      unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:do k = 1, this%NInCutoff(unit)
        j = this%CutoffPartner(k, unit)
        if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
          if (mod(j,nu2)==0) then
            jk = INT(j/nu2)
          else
            jk = INT(j/nu2)+1
          end if
          RXij = RXi - RX2(jk)
          RYij = RYi - RY2(jk)
          RZij = RZi - RZ2(jk)
          PXij = PXi - PX2(jk)
          PYij = PYi - PY2(jk)
          PZij = PZi - PZ2(jk)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          OXj = OX2(jk)
          OYj = OY2(jk)
          OZj = OZ2(jk)
          RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
          RijInv = sqrt( RijSquaredInv )
          eX = RXij * RijInv                                                      ! Einheitsabstandvektor nach Price
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosTheta  = OXj * ex + OYj * eY + OZj * eZ                              ! cos(alpha) nach Price
          CosTheta3 = 3._RK * CosTheta
          Epsilon1 = Epsilon * RijSquaredInv
          Epsilon2 = Epsilon1 * RijInv
          EPotlocal1 = Epsilon1 * CosTheta              !Define EPotLocal 1 
          EPotLocal  = EPotLocal + EPotlocal1                         ! Uebereinstimmumg mit Price
          EPotLocalInter  = EPotLocalInter + EPotLocal1
          FXij = Epsilon2 * ( CosTheta3 * eX - OXj )                              ! F2 bei Price
          FYij = Epsilon2 * ( CosTheta3 * eY - OYj )
          FZij = Epsilon2 * ( CosTheta3 * eZ - OZj )
          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F2*R_COM_Price; stimmt so
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
        VirialLocalInter = VirialLocalInter + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F2*R_COM_Price; stimmt so
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + EPotlocal1*(8._RK*sitecorr*sitecorr-2._RK*Plen2*RijSquaredInv)*Third*Third   !xxxx2 CD T
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        forceTempX(jk) = forceTempX(jk) - FXij
        forceTempY(jk) = forceTempY(jk) - FYij
        forceTempZ(jk) = forceTempZ(jk) - FZij
        momTempX(jk) = momTempX(jk) - Epsilon1 * eX   
        momTempY(jk) = momTempY(jk) - Epsilon1 * eY   
        momTempZ(jk) = momTempZ(jk) - Epsilon1 * eZ   

#if TRANS==1
        !TRANSPORT_start vielleicht
        VSxi   = VSxi + FXij * PYij
        VSyi   = VSyi + FXij * PZij
        VSzi   = VSzi + FYij * PZij
        VBxi   = VBxi + FXij * PXij
        VByi   = VByi + FYij * PYij
        VBzi   = VBzi + FZij * PZij
        VSuxi  = VSuxi+ FYij * PXij
        VSuyi  = VSuyi+ FZij * PXij
        VSuzi  = VSuzi+ FZij * PYij
        UU    =  EpotLocal1        !EpotLocal1 not defined
        Uxi   =  UU * eX
        Uyi   =  UU * eY
        Uzi   =  UU * eZ
        Cxi    = Cxi  + Uxi !Why was this term left out?
        Cyi    = Cyi  + Uyi !Why was this term left out?
        Czi    = Czi  + Uzi !Why was this term left out?
        txii   = r1y * FZij - r1z * FYij
        tyii   = r1z * FXij - r1x * FZij
        tzii   = r1x * FYij - r1y * FXij
        txi    = A11 * txii + A12 * tyii + A13 * tzii
        tyi    = A21 * txii + A22 * tyii + A23 * tzii
        tzi    = A31 * txii + A32 * tyii + A33 * tzii
        tuxi   = tuxi + PXij*tyi
        tuyi   = tuyi + PXij*tzi
        tuzi   = tuzi + PYij*tzi
        tlxi   = tlxi + PYij*txi
        tlyi   = tlyi + PZij*txi
        tlzi   = tlzi + PZij*tyi
        tdxi   = tdxi + PXij*txi
        tdyi   = tdyi + PYij*tyi
        tdzi   = tdzi + PZij*tzi
#endif
        end if
      end do loop1
      ! Include intramolecular interactions if need
      if (this%potintra15 .or. this%potintra14) then
        RXij = RXi - RX2(i)
        RYij = RYi - RY2(i)
        RZij = RZi - RZ2(i)
        PXij = PXi - PX2(i)
        PYij = PYi - PY2(i)
        PZij = PZi - PZ2(i)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        OXj = OX2(i)
        OYj = OY2(i)
        OZj = OZ2(i)
        RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
        RijInv = sqrt( RijSquaredInv )
        eX = RXij * RijInv                                                      ! Einheitsabstandvektor nach Price
        eY = RYij * RijInv
        eZ = RZij * RijInv
        CosTheta  = OXj * ex + OYj * eY + OZj * eZ                              ! cos(alpha) nach Price
        CosTheta3 = 3._RK * CosTheta
        Epsilon1 = Epsilon * RijSquaredInv*coeff                                ! 1-4 non-bonded interaction coeff
        Epsilon2 = Epsilon1 * RijInv
        EPotLocal1 = Epsilon1 * CosTheta
        EPotLocal  = EPotLocal + EPotLocal1                          ! Uebereinstimmumg mit Price
        EPotLocalIntra  = EPotLocalIntra + EPotLocal1
        FXij = Epsilon2 * ( CosTheta3 * eX - OXj )                              ! F2 bei Price
        FYij = Epsilon2 * ( CosTheta3 * eY - OYj )
        FZij = Epsilon2 * ( CosTheta3 * eZ - OZj )
        VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F2*R_COM_Price; stimmt so
        VirialLocalIntra = VirialLocalIntra + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F2*R_COM_Price; stimmt so
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + EPotlocal1*(8._RK*sitecorr*sitecorr-2._RK*Plen2*RijSquaredInv)*Third*Third   !xxxx2 CD T
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        forceTempX(i) = forceTempX(i) - FXij
        forceTempY(i) = forceTempY(i) - FYij
        forceTempZ(i) = forceTempZ(i) - FZij
        momTempX(i) = momTempX(i) - Epsilon1 * eX   
        momTempY(i) = momTempY(i) - Epsilon1 * eY   
        momTempZ(i) = momTempZ(i) - Epsilon1 * eZ

      end if

      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi

#if  TRANS == 1
      !TRANSPORT_start
      VSx(i) = VSxi
      VSy(i) = VSyi
      VSz(i) = VSzi
      VBx(i) = VBxi
      VBy(i) = VByi
      VBz(i) = VBzi
      VSux(i)= VSuxi
      VSuy(i)= VSuyi
      VSuz(i)= VSuzi
      Cx(i)  = Cxi
      Cy(i)  = Cyi
      Cz(i)  = Czi
      tux(i) = tuxi
      tuy(i) = tuyi
      tuz(i) = tuzi
      tlx(i) = tlxi
      tly(i) = tlyi
      tlz(i) = tlzi
      tdx(i) = tdxi
      tdy(i) = tdyi
      tdz(i) = tdzi
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

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
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
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
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
    integer           :: nu1, nu2, unit, i0, jk
#if ARCH == 3
    logical           :: hit
#endif

    ! Assign local variables
    i1 = this%Site1%NTest
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared
    nu1 = this%NUnit1
    nu2 = this%NUnit2
#if MPI_VER > 0
    i0 = this%Site1%NTest0
    i1 = this%Site1%NTest2
#else
    i0 = 1
    i1 = this%Site1%NTest
#endif

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
   do i = i0, i1
     RXi = RX1(i)
     RYi = RY1(i)
     RZi = RZ1(i)
     PXi = PX1(i)
     PYi = PY1(i)
     PZi = PZ1(i)
     EPotLocal = 0._RK
#if ARCH == 3
     hit = .false.
#endif

     unit = nu1*(i-1)+this%Site1%UnitNumber

!CDIR NODEP
loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then  ! choose only units, to which our Site2 correspond
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)   ! number of molecule, to which this unit correspond
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            PXij = (PXij - anint( PXij )) * BoxLength
            PYij = (PYij - anint( PYij )) * BoxLength
            PZij = (PZij - anint( PZij )) * BoxLength
            OXj = OX2(jk)
            OYj = OY2(jk)
            OZj = OZ2(jk)
            RijSquared = RXij**2 + RYij**2 + RZij**2
#if ARCH == 3
            if( RijSquared <= RShieldSquared ) hit = .true.
#else
            if( RijSquared <= RShieldSquared ) then
              EPotLocal = 1E33_RK
              exit loop1
            end if
#endif
            RijSquaredInv = 1._RK / RijSquared
            RijInv = sqrt( RijSquaredInv )
            eX = RXij * RijInv
            eY = RYij * RijInv
            eZ = RZij * RijInv
            CosTheta  = OXj * ex + OYj * eY + OZj * eZ
            EPotLocal  = EPotLocal + Epsilon * RijSquaredInv * CosTheta
          end if
        end do loop1

#if ARCH == 3
        if( .not. hit ) then
          EPotTest(i) = EPotTest(i) + EPotLocal
        else
          EPotTest(i) = EPotTest(i) + 1E33_RK
        endif
#else
        EPotTest(i) = EPotTest(i) + EPotLocal
#endif
   end do

  end subroutine TPotCD_ChemicalPotential



!==============================================================!
!  Subroutine TPotCD_Energy                                    !
!==============================================================!

  subroutine TPotCD_Energy( this, np, nu, F, E, EIntra, BoxLength, CompIdent )

    implicit none

    ! Declare arguments
    type(TPotChargeDipole) :: this
    integer, intent(in)    :: np
    integer, intent(in)      :: nu
    real(RK), intent(in out) :: F(3,nu)
    real(RK), intent(in out) :: E
    real(RK), intent(in out) :: EIntra
    real(RK), intent(in)   :: BoxLength
    logical, intent(in)      :: CompIdent

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1
    real(RK)          :: RShieldSquared
    real(RK), pointer, contiguous :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer, contiguous :: OX2(:), OY2(:), OZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv, RijSquared
    real(RK)          :: EPotLocal, EPot, EIntra1, tempF(3,nu)
    real(RK)          :: CosTheta
    integer           :: j, k, nu2, jk, unit
    real(RK)          :: coeff

    ! Assign local variables
    Epsilon = this%Epsilon
    RShieldSquared = this%RShieldSquared
    nu2 = this%NUnit2
    coeff = 1._Rk
    if (this%potintra14) coeff = this%ScaleEl14
    EPot = 0._RK
    EIntra1 = 0._RK

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
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ

    ! Loop over molecules
    RXi = RX1(np)
    RYi = RY1(np)
    RZi = RZ1(np)
    PXi = PX1(np)
    PYi = PY1(np)
    PZi = PZ1(np)

    unit=this%NUnit1*(np-1)+this%Site1%UnitNumber

    do k = 1, this%NInCutoff(unit)
      j = this%CutoffPartner(k, unit)
      if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
        if (mod(j,nu2)==0) then
          jk = INT(j/nu2)
        else
          jk = INT(j/nu2)+1
        end if
        RXij = RXi - RX2(jk)
        RYij = RYi - RY2(jk)
        RZij = RZi - RZ2(jk)
        PXij = PXi - PX2(jk)
        PYij = PYi - PY2(jk)
        PZij = PZi - PZ2(jk)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        OXj = OX2(jk)
        OYj = OY2(jk)
        OZj = OZ2(jk)
        RijSquared = RXij**2 + RYij**2 + RZij**2

        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
        else
          RijSquaredInv = 1._RK / RijSquared
          RijInv = sqrt( RijSquaredInv )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosTheta  = OXj * ex + OYj * eY + OZj * eZ
          Epsilon1 = Epsilon * RijSquaredInv
          EPotLocal  = Epsilon1 * 3._RK * CosTheta
        end if
        EPot = EPot + EPotLocal
      end if
    end do
    ! Include intramolecular interactions if need
    if (this%potintra15 .or. this%potintra14) then
      RXij = RXi - RX2(np)
      RYij = RYi - RY2(np)
      RZij = RZi - RZ2(np)
      PXij = PXi - PX2(np)
      PYij = PYi - PY2(np)
      PZij = PZi - PZ2(np)
      RXij = (RXij - anint( PXij )) * BoxLength
      RYij = (RYij - anint( PYij )) * BoxLength
      RZij = (RZij - anint( PZij )) * BoxLength
      OXj = OX2(np)
      OYj = OY2(np)
      OZj = OZ2(np)
      RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
      RijInv = sqrt( RijSquaredInv )
      eX = RXij * RijInv
      eY = RYij * RijInv
      eZ = RZij * RijInv
      CosTheta  = OXj * ex + OYj * eY + OZj * eZ
      EPotLocal = Epsilon * RijSquaredInv*coeff * 3._RK * CosTheta
      EIntra1  = EIntra1 + EPotLocal
    end if

    F(:,:) = F(:,:) + tempF(:,:)
    E = E + EPot + EIntra1
    EIntra = EIntra + EIntra1

  end subroutine TPotCD_Energy


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

    ! Declare local variables
    integer :: k, ende

    ! Construct potential
    this%Site1 => Molecule1%SiteCharge(j1)
    this%Site2 => Molecule2%SiteQuadrupole(j2)
    this%NUnit1 => Molecule1%NUnit
    this%NUnit2 => Molecule2%NUnit
    this%SameComponent = i1 == i2
    this%Epsilon = 1.5_RK * this%Site1%e * this%Site2%Q
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

    ! if this potential is intra
    if (this%SameComponent .and. Molecule1%hasIntraLJEl) then
      ende = size(Molecule1%IntCQ15(:,1))
      do k=1, ende
        if (Molecule1%IntCQ15(k,1)==this%Site1%SiteId .and. Molecule1%IntCQ15(k,2)==this%Site2%SiteId) then
          this%potintra15 = .true.
        end if
      end do
      if (LJEL14 .and. .not. this%potintra15) then
        ende = size(Molecule1%IntCQ14(:,1))
        do k=1, ende
          if (Molecule1%IntCQ14(k,1)==this%Site1%SiteId .and. Molecule1%IntCQ14(k,2)==this%Site2%SiteId) then
            this%potintra14=.true.
            this%ScaleEl14 = Molecule1%ScaleCQ14(k)
          end if
        end do
      end if
   else
      this%potintra15 = .false.
      this%potintra14 = .false.
    end if

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

  subroutine TPotCQ_Force( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeQuadrupole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: EPotInter
    real(RK), intent(in out)   :: VirialInter
    real(RK), intent(in out)   :: EPotIntra_Nonbonded
    real(RK), intent(in out)   :: VirialIntra
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
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, jk, unit
    logical           :: SameComponent
    real(RK)          :: coeff
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
    SameComponent = this%SameComponent
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14 !Scale 1,4 El interactions
    else
      coeff = 1._RK
    end if

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

      unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:do k = 1, this%NInCutoff(unit)
        j = this%CutoffPartner(k, unit)
        if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
          if (mod(j,nu2)==0) then
            jk = INT(j/nu2)
          else
            jk = INT(j/nu2)+1
          end if
          RXij = RXi - RX2(jk)
          RYij = RYi - RY2(jk)
          RZij = RZi - RZ2(jk)
          PXij = PXi - PX2(jk)
          PYij = PYi - PY2(jk)
          PZij = PZi - PZ2(jk)
          RXij = (RXij - anint( PXij )) * BoxLength                               ! Abstandsvektor von Q nach C wie bei Price
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          OXj = OX2(jk)                                                            ! Orientierungsvektor Quadrupol
          OYj = OY2(jk)
          OZj = OZ2(jk)
          RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
          RijInv = sqrt( RijSquaredInv )
          eX = RXij * RijInv                                                      ! Normierter Abstandsvektor
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosTheta  = OXj * ex + OYj * eY + OZj * eZ
          Epsilon1 = Epsilon * RijSquaredInv * RijInv
          EPotLocal  = EPotLocal + Epsilon1 * ( CosTheta * CosTheta - Third )
          EPotLocalInter  = EPotLocalInter + Epsilon1 * ( CosTheta * CosTheta - Third )
          CosTheta2 = 2._RK * CosTheta
          CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
          Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
          FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXj )                     ! F2 nach Price bzw. Kraft auf Punktladung
          FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYj )
          FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZj )
          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! Vorzeichen richtig so
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
        VirialLocalInter = VirialLocalInter + (FXij * PXij + FYij * PYij + FZij * PZij)
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + Epsilon1 * ( CosTheta * CosTheta - Third )*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijSquaredInv)*Third*Third   !xxxx3 CQ

        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        
        forceTempX(jk) = forceTempX(jk) - FXij
        forceTempY(jk) = forceTempY(jk) - FYij
        forceTempZ(jk) = forceTempZ(jk) - FZij
        momTempX(jk) = momTempX(jk) - Epsilon1 * CosTheta2 * eX   
        momTempY(jk) = momTempY(jk) - Epsilon1 * CosTheta2 * eY   
        momTempZ(jk) = momTempZ(jk) - Epsilon1 * CosTheta2 * eZ   

        end if

      end do loop1
      ! Include intramolecular interactions if need
      if (SameComponent .and. (this%potintra14 .or. this%potintra15)) then
        RXij = RXi - RX2(i)
        RYij = RYi - RY2(i)
        RZij = RZi - RZ2(i)
        PXij = PXi - PX2(i)
        PYij = PYi - PY2(i)
        PZij = PZi - PZ2(i)
        RXij = (RXij - anint( PXij )) * BoxLength                               ! Abstandsvektor von Q nach C wie bei Price
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        OXj = OX2(i)                                                            ! Orientierungsvektor Quadrupol
        OYj = OY2(i)
        OZj = OZ2(i)
        RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
        RijInv = sqrt( RijSquaredInv )
        eX = RXij * RijInv                                                      ! Normierter Abstandsvektor
        eY = RYij * RijInv
        eZ = RZij * RijInv
        CosTheta  = OXj * ex + OYj * eY + OZj * eZ
        Epsilon1 = Epsilon * RijSquaredInv * RijInv * coeff
        EPotLocal  = EPotLocal + Epsilon1 * ( CosTheta * CosTheta - Third )
        EPotLocalIntra  = EPotLocalIntra + Epsilon1 * ( CosTheta * CosTheta - Third )
        CosTheta2 = 2._RK * CosTheta
        CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
        Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv*coeff
        FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXj )                      ! F2 nach Price bzw. Kraft auf Punktladung
        FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYj )
        FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZj )
        VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! Vorzeichen richtig so
        VirialLocalIntra = VirialLocalIntra + (FXij * PXij + FYij * PYij + FZij * PZij)
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + Epsilon1 * ( CosTheta * CosTheta - Third ) *(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijSquaredInv)*Third*Third   !xxxx3 CQ
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        forceTempX(i) = forceTempX(i) - FXij
        forceTempY(i) = forceTempY(i) - FYij
        forceTempZ(i) = forceTempZ(i) - FZij
        momTempX(i) = momTempX(i) - Epsilon1 * CosTheta2 * eX   
        momTempY(i) = momTempY(i) - Epsilon1 * CosTheta2 * eY   
        momTempZ(i) = momTempZ(i) - Epsilon1 * CosTheta2 * eZ
      end if

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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotCQ_Force

!==============================================================!
!  Subroutine TPotCQ_Force_Trans                               !
!==============================================================!

  subroutine TPotCQ_Force_Trans( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeQuadrupole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: EPotInter
    real(RK), intent(in out)   :: VirialInter
    real(RK), intent(in out)   :: EPotIntra_Nonbonded
    real(RK), intent(in out)   :: VirialIntra
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
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, jk, unit
    logical           :: SameComponent
    real(RK)          :: coeff
    real(RK)          :: forceTempX(1:this%Site2%NPart)
    real(RK)          :: forceTempY(1:this%Site2%NPart)
    real(RK)          :: forceTempZ(1:this%Site2%NPart)
    real(RK)          :: momTempX(1:this%Site2%NPart)
    real(RK)          :: momTempY(1:this%Site2%NPart)
    real(RK)          :: momTempZ(1:this%Site2%NPart)
 
#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer, contiguous :: VSx(:), VSy(:), VSz(:)
    real(RK), pointer, contiguous :: VSux(:), VSuy(:), VSuz(:)
    real(RK), pointer, contiguous :: VBx(:), VBy(:), VBz(:)
    real(RK), pointer, contiguous :: Cx(:) , Cy(:) , Cz(:)
    real(RK), pointer, contiguous :: tux(:) , tuy(:) , tuz(:)
    real(RK), pointer, contiguous :: tlx(:) , tly(:) , tlz(:)
    real(RK), pointer, contiguous :: tdx(:) , tdy(:) , tdz(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
    real(RK)          :: VSxi, VSyi, VSzi
    real(RK)          :: VSuxi,VSuyi,VSuzi
    real(RK)          :: VBxi, VByi, VBzi
    real(RK)          :: Cxi,  Cyi,  Czi
    real(RK)          :: tuxi,  tuyi,  tuzi
    real(RK)          :: tlxi,  tlyi,  tlzi
    real(RK)          :: tdxi,  tdyi,  tdzi
    real(RK)          :: txii,  tyii , tzii
    real(RK)          :: txi ,  tyi  , tzi
    real(RK)          :: UU, Uxi,  Uyi, Uzi
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

!$OMP PARALLEL &
!$OMP PRIVATE ( Epsilon, Epsilon1, Epsilon2, RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (  FX1, FY1, FZ1, OX2, OY2, OZ2) &
!$OMP PRIVATE ( Plen2, PX1, PY1, PZ1, PX2, PY2, PZ2) &
!$OMP PRIVATE (   RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE (   RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE ( OXj, OYj, OZj, eX, eY, eZ, RijSquaredInv, RijInv) &
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
    SameComponent = this%SameComponent
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14 !Scale 1,4 El interactions
    else
      coeff = 1._RK
    end if

#if  TRANS == 1
    !TRANSPORT_start
    VSx => this%Site1%vsCx
    VSy => this%Site1%vsCy
    VSz => this%Site1%vsCz
    VBx => this%Site1%vbCx
    VBy => this%Site1%vbCy
    VBz => this%Site1%vbCz
    VSux=> this%Site1%vsuCx
    VSuy=> this%Site1%vsuCy
    VSuz=> this%Site1%vsuCz
    Cx  => this%Site1%cCx
    Cy  => this%Site1%cCy
    Cz  => this%Site1%cCz
    tux => this%Site1%tuCx
    tuy => this%Site1%tuCy
    tuz => this%Site1%tuCz
    tlx => this%Site1%tlCx
    tly => this%Site1%tlCy
    tlz => this%Site1%tlCz
    tdx => this%Site1%tdCx
    tdy => this%Site1%tdCy
    tdz => this%Site1%tdCz
    q1  => this%Site1%Qm0r(:, 1, 1)
    q2  => this%Site1%Qm0r(:, 2, 1)
    q3  => this%Site1%Qm0r(:, 3, 1)
    q4  => this%Site1%Qm0r(:, 4, 1)
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
      VSxi= VSx(i)
      VSyi= VSy(i)
      VSzi= VSz(i)
      VBxi= VBx(i)
      VByi= VBy(i)
      VBzi= VBz(i)
      VSuxi= VSux(i)
      VSuyi= VSuy(i)
      VSuzi= VSuz(i)
      Cxi = Cx(i)
      Cyi = Cy(i)
      Czi = Cz(i)
      tuxi = tux(i)
      tuyi = tuy(i)
      tuzi = tuz(i)
      tlxi = tlx(i)
      tlyi = tly(i)
      tlzi = tlz(i)
      tdxi = tdx(i)
      tdyi = tdy(i)
      tdzi = tdz(i)
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

      unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:do k = 1, this%NInCutoff(unit)
        j = this%CutoffPartner(k, unit)
        if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
          if (mod(j,nu2)==0) then
            jk = INT(j/nu2)
          else
            jk = INT(j/nu2)+1
          end if
          RXij = RXi - RX2(jk)
          RYij = RYi - RY2(jk)
          RZij = RZi - RZ2(jk)
          PXij = PXi - PX2(jk)
          PYij = PYi - PY2(jk)
          PZij = PZi - PZ2(jk)
          RXij = (RXij - anint( PXij )) * BoxLength                               ! Abstandsvektor von Q nach C wie bei Price
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          OXj = OX2(jk)                                                            ! Orientierungsvektor Quadrupol
          OYj = OY2(jk)
          OZj = OZ2(jk)
          RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
          RijInv = sqrt( RijSquaredInv )
          eX = RXij * RijInv                                                      ! Normierter Abstandsvektor
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosTheta  = OXj * ex + OYj * eY + OZj * eZ
          Epsilon1 = Epsilon * RijSquaredInv * RijInv
          EPotLocal1 = Epsilon1 * ( CosTheta * CosTheta - Third ) !Gabriela: Definition of EpotLocal1 
          EPotLocal = EPotLocal + EPotlocal1
          EPotLocalInter  = EPotLocalInter + EPotLocal1
          CosTheta2 = 2._RK * CosTheta
          CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
          Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
          FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXj )                     ! F2 nach Price bzw. Kraft auf Punktladung
          FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYj )
          FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZj )
          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! Vorzeichen richtig so
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
        VirialLocalInter = VirialLocalInter + (FXij * PXij + FYij * PYij + FZij * PZij)
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijSquaredInv)*Third*Third   !xxxx3 CQ T

        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        
        forceTempX(jk) = forceTempX(jk) - FXij
        forceTempY(jk) = forceTempY(jk) - FYij
        forceTempZ(jk) = forceTempZ(jk) - FZij
        momTempX(jk) = momTempX(jk) - Epsilon1 * CosTheta2 * eX   
        momTempY(jk) = momTempY(jk) - Epsilon1 * CosTheta2 * eY   
        momTempZ(jk) = momTempZ(jk) - Epsilon1 * CosTheta2 * eZ   

#if TRANS==1
        !TRANSPORT_start vielleicht
        VSxi   = VSxi + FXij * PYij
        VSyi   = VSyi + FXij * PZij
        VSzi   = VSzi + FYij * PZij
        VBxi   = VBxi + FXij * PXij
        VByi   = VByi + FYij * PYij
        VBzi   = VBzi + FZij * PZij
        VSuxi  = VSuxi+ FYij * PXij
        VSuyi  = VSuyi+ FZij * PXij
        VSuzi  = VSuzi+ FZij * PYij
        UU    =  EpotLocal1  ! Change: EPotLocal1 was not defined before .... 
        Uxi   =  UU * eX
        Uyi   =  UU * eY
        Uzi   =  UU * eZ   
        Cxi    = Cxi  + Uxi !Why was this term left out?
        Cyi    = Cyi  + Uyi !Why was this term left out?
        Czi    = Czi  + Uzi !Why was this term left out?
        txii   = r1y * FZij - r1z * FYij
        tyii   = r1z * FXij - r1x * FZij
        tzii   = r1x * FYij - r1y * FXij
        txi    = A11 * txii + A12 * tyii + A13 * tzii
        tyi    = A21 * txii + A22 * tyii + A23 * tzii
        tzi    = A31 * txii + A32 * tyii + A33 * tzii
        tuxi   = tuxi + PXij*tyi
        tuyi   = tuyi + PXij*tzi
        tuzi   = tuzi + PYij*tzi
        tlxi   = tlxi + PYij*txi
        tlyi   = tlyi + PZij*txi
        tlzi   = tlzi + PZij*tyi
        tdxi   = tdxi + PXij*txi
        tdyi   = tdyi + PYij*tyi
        tdzi   = tdzi + PZij*tzi
#endif

        end if
      end do loop1
      ! Include intramolecular interactions if need
      if (SameComponent .and. (this%potintra14 .or. this%potintra15)) then
        RXij = RXi - RX2(i)
        RYij = RYi - RY2(i)
        RZij = RZi - RZ2(i)
        PXij = PXi - PX2(i)
        PYij = PYi - PY2(i)
        PZij = PZi - PZ2(i)
        RXij = (RXij - anint( PXij )) * BoxLength                               ! Abstandsvektor von Q nach C wie bei Price
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        OXj = OX2(i)                                                            ! Orientierungsvektor Quadrupol
        OYj = OY2(i)
        OZj = OZ2(i)
        RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
        RijInv = sqrt( RijSquaredInv )
        eX = RXij * RijInv                                                      ! Normierter Abstandsvektor
        eY = RYij * RijInv
        eZ = RZij * RijInv
        CosTheta  = OXj * ex + OYj * eY + OZj * eZ
        Epsilon1 = Epsilon * RijSquaredInv * RijInv * coeff
        EPotLocal1 = Epsilon1 * ( CosTheta * CosTheta - Third )
        EPotLocal  = EPotLocal + EPotLocal1
        EPotLocalIntra  = EPotLocalIntra + EPotLocal1
        CosTheta2 = 2._RK * CosTheta
        CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
        Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv*coeff
        FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXj )                      ! F2 nach Price bzw. Kraft auf Punktladung
        FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYj )
        FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZj )
        VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
        VirialLocalIntra = VirialLocalIntra + (FXij * PXij + FYij * PYij + FZij * PZij)
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijSquaredInv)*Third*Third   !xxxx3 CQ
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        forceTempX(i) = forceTempX(i) - FXij
        forceTempY(i) = forceTempY(i) - FYij
        forceTempZ(i) = forceTempZ(i) - FZij
        momTempX(i) = momTempX(i) - Epsilon1 * CosTheta2 * eX   
        momTempY(i) = momTempY(i) - Epsilon1 * CosTheta2 * eY   
        momTempZ(i) = momTempZ(i) - Epsilon1 * CosTheta2 * eZ

      end if

      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
#if  TRANS == 1
      !TRANSPORT_start
      VSx(i) = VSxi
      VSy(i) = VSyi
      VSz(i) = VSzi
      VBx(i) = VBxi
      VBy(i) = VByi
      VBz(i) = VBzi
      VSux(i)= VSuxi
      VSuy(i)= VSuyi
      VSuz(i)= VSuzi
      Cx(i)  = Cxi
      Cy(i)  = Cyi
      Cz(i)  = Czi
      tux(i) = tuxi
      tuy(i) = tuyi
      tuz(i) = tuzi
      tlx(i) = tlxi
      tly(i) = tlyi
      tlz(i) = tlzi
      tdx(i) = tdxi
      tdy(i) = tdyi
      tdz(i) = tdzi
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

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal
#if OSMOP == 2
    this%VirialProfile(:) = Third * this%VirialProfile(:)
#endif
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
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
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
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
    integer           :: nu1, nu2, unit, i0, jk
#if ARCH == 3
    logical           :: hit
#endif

    ! Assign local variables
    i1 = this%Site1%NTest
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared
    nu1 = this%NUnit1
    nu2 = this%NUnit2
#if MPI_VER > 0
    i0 = this%Site1%NTest0
    i1 = this%Site1%NTest2
#else
    i0 = 1
    i1 = this%Site1%NTest
#endif

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
   do i = i0, i1
     RXi = RX1(i)
     RYi = RY1(i)
     RZi = RZ1(i)
     PXi = PX1(i)
     PYi = PY1(i)
     PZi = PZ1(i)
     EPotLocal = 0._RK

#if ARCH == 3
        hit = .false.
#endif
!CDIR NODEP

     unit = nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then  ! choose only units, to which our Site2 correspond
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)   ! number of molecule, to which this unit correspond
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            PXij = (PXij - anint( PXij )) * BoxLength
            PYij = (PYij - anint( PYij )) * BoxLength
            PZij = (PZij - anint( PZij )) * BoxLength
            OXj = OX2(jk)
            OYj = OY2(jk)
            OZj = OZ2(jk)
            RijSquared = RXij**2 + RYij**2 + RZij**2
#if ARCH == 3
            if( RijSquared <= RShieldSquared ) hit = .true.
#else
            if( RijSquared <= RShieldSquared ) then
              EPotLocal = 1E33_RK
              exit loop1
            end if
#endif
            RijSquaredInv = 1._RK / RijSquared
            RijInv = sqrt( RijSquaredInv )
            eX = RXij * RijInv
            eY = RYij * RijInv
            eZ = RZij * RijInv
            CosTheta  = OXj * ex + OYj * eY + OZj * eZ
            EPotLocal  = EPotLocal + Epsilon * RijSquaredInv * RijInv * ( CosTheta * CosTheta - Third )
          end if
        end do loop1

#if ARCH == 3
        if( .not. hit ) then
          EPotTest(i) = EPotTest(i) + EPotLocal
        else
          EPotTest(i) = EPotTest(i) + 1E33_RK
        endif
#else
        EPotTest(i) = EPotTest(i) + EPotLocal
#endif
   end do

  end subroutine TPotCQ_ChemicalPotential


!==============================================================!
!  Subroutine TPotCQ_Energy                                    !
!==============================================================!

  subroutine TPotCQ_Energy( this, np, nu, F, E, EIntra, BoxLength, CompIdent )

    implicit none

    ! Declare arguments
    type(TPotChargeQuadrupole) :: this
    integer, intent(in)        :: np
    integer, intent(in)        :: nu
    real(RK), intent(in out)   :: F(3,nu)
    real(RK), intent(in out)   :: E
    real(RK), intent(in out)   :: EIntra
    real(RK), intent(in)       :: BoxLength
    logical, intent(in)        :: CompIdent

    ! Declare local variables
    real(RK)          :: Epsilon
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
    real(RK)          :: EPot, EIntra1, EPotLocal, tempF(3,nu)
    integer           :: j, k, nu2, jk, unit
    real(RK)          :: coeff

    ! Assign local variables
    Epsilon = this%Epsilon
    RShieldSquared = this%RShieldSquared
    nu2 = this%NUnit2
    coeff = 1._RK
    if (this%potintra14) coeff = this%ScaleEl14
    EPot = 0._RK
    EIntra1 = 0._RK

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
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ

    ! Loop over molecules
    RXi = RX1(np)
    RYi = RY1(np)
    RZi = RZ1(np)
    PXi = PX1(np)
    PYi = PY1(np)
    PZi = PZ1(np)

    unit=this%NUnit1*(np-1)+this%Site1%UnitNumber

    do k = 1, this%NInCutoff(unit)
      j = this%CutoffPartner(k, unit)
      if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
        if (mod(j,nu2)==0) then
          jk = INT(j/nu2)
        else
          jk = INT(j/nu2)+1
        end if
        RXij = RXi - RX2(jk)
        RYij = RYi - RY2(jk)
        RZij = RZi - RZ2(jk)
        PXij = PXi - PX2(jk)
        PYij = PYi - PY2(jk)
        PZij = PZi - PZ2(jk)
        RXij = (RXij - anint( PXij )) * BoxLength                                 ! Abstandsvektor von Q nach C wie bei Price
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        OXj = OX2(jk)
        OYj = OY2(jk)
        OZj = OZ2(jk)
        RijSquared = RXij**2 + RYij**2 + RZij**2

        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
        else
          RijSquaredInv = 1._RK / RijSquared
          RijInv = sqrt( RijSquaredInv )
          eX = RXij * RijInv                                                      ! Normierter Abstandsvektor
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosTheta  = OXj * ex + OYj * eY + OZj * eZ
          EPotLocal = Epsilon * RijSquaredInv * RijInv * ( CosTheta * CosTheta - Third )
        end if
        EPot = EPot + EPotLocal
      end if
    end do
    ! Include intramolecular interactions if need
    if (this%potintra14 .or. this%potintra15) then
      RXij = RXi - RX2(np)
      RYij = RYi - RY2(np)
      RZij = RZi - RZ2(np)
      PXij = PXi - PX2(np)
      PYij = PYi - PY2(np)
      PZij = PZi - PZ2(np)
      RXij = (RXij - anint( PXij )) * BoxLength
      RYij = (RYij - anint( PYij )) * BoxLength
      RZij = (RZij - anint( PZij )) * BoxLength
      OXj = OX2(np)
      OYj = OY2(np)
      OZj = OZ2(np)
      RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
      RijInv = sqrt( RijSquaredInv )
      eX = RXij * RijInv
      eY = RYij * RijInv
      eZ = RZij * RijInv
      CosTheta  = OXj * ex + OYj * eY + OZj * eZ
      EPotLocal = Epsilon * RijSquaredInv * RijInv * coeff * ( CosTheta * CosTheta - Third )
      EIntra1 = EIntra1 + EPotLocal
    end if

    F(:,:) = F(:,:) + tempF(:,:)
    E = E + EPot + EIntra1
    EIntra = EIntra + EIntra1

  end subroutine TPotCQ_Energy


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

    ! Declare local variables
    integer :: k, ende

    ! Construct potential
    this%Site1 => Molecule1%SiteDipole(j1)
    this%Site2 => Molecule2%SiteCharge(j2)
    this%NUnit1 => Molecule1%NUnit
    this%NUnit2 => Molecule2%NUnit
    this%SameComponent = i1 == i2
    this%Epsilon = this%Site1%D * this%Site2%e
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

    ! if this potential is intra
    if (this%SameComponent .and. Molecule1%hasIntraLJEl) then
      ende = size(Molecule1%IntDC15(:,1))
      do k=1, ende
        if (Molecule1%IntDC15(k,1)==this%Site1%SiteId .and. Molecule1%IntDC15(k,2)==this%Site2%SiteId) then
          this%potintra15 = .true.
        end if
      end do
      if (LJEL14 .and. .not. this%potintra15) then
        ende = size(Molecule1%IntDC14(:,1))
        do k=1, ende
          if (Molecule1%IntDC14(k,1)==this%Site1%SiteId .and. Molecule1%IntDC14(k,2)==this%Site2%SiteId) then
             this%potintra14=.true.
             this%ScaleEl14 = Molecule1%ScaleDC14(k)
          end if
        end do
      end if
    else
      this%potintra15 = .false.
      this%potintra14 = .false.
    end if

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

  subroutine TPotDC_Force( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra_Nonbonded
    real(RK), intent(in out) :: VirialIntra
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
    real(RK)          :: EPotLocalIntra, ViriallocalIntra
    real(RK)          :: EPotLocalInter, ViriallocalInter
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, jk, unit
    logical           :: SameComponent
    real(RK)          :: coeff
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
    SameComponent = this%SameComponent
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

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

      unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:do k = 1, this%NInCutoff(unit)
        j = this%CutoffPartner(k, unit)
        if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
          if (mod(j,nu2)==0) then
            jk = INT(j/nu2)
          else
            jk = INT(j/nu2)+1
          end if
          RXij = RXi - RX2(jk)
          RYij = RYi - RY2(jk)
          RZij = RZi - RZ2(jk)
          PXij = PXi - PX2(jk)
          PYij = PYi - PY2(jk)
          PZij = PZi - PZ2(jk)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
          RijInv = sqrt( RijSquaredInv )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosTheta  = OXi * ex + OYi * eY + OZi * eZ                              ! -cos(alpha) bei Price
          CosTheta3 = 3._RK * CosTheta
          Epsilon1 = Epsilon * RijSquaredInv
          Epsilon2 = Epsilon1 * RijInv
          EPotLocal  = EPotLocal - Epsilon1 * CosTheta                          ! Uebereinstimmumg mit Price
          EPotLocalInter  = EPotLocalInter - Epsilon1 * CosTheta
          FXij = Epsilon2 * ( OXi - CosTheta3 * eX )                              ! F1 bei Price
          FYij = Epsilon2 * ( OYi - CosTheta3 * eY )
          FZij = Epsilon2 * ( OZi - CosTheta3 * eZ )
          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F1*(-R_COM_Price); stimmt so
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
          VirialLocalInter = VirialLocalInter + (FXij * PXij + FYij * PYij + FZij * PZij)
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
          d2EpotdV2Local = d2EpotdV2Local - Epsilon1*CosTheta*(8._RK*sitecorr*sitecorr-2._RK*Plen2*RijSquaredInv)*Third*Third          !xxxx4  DC
          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          TXi = TXi + Epsilon1 * eX   ! Uebereinstimmumg mit Price; Rest bei Atom2Mol in Component
          TYi = TYi + Epsilon1 * eY   ! Reaktionsfeldbeitrag in Interaction
          TZi = TZi + Epsilon1 * eZ
          forceTempX(jk) = forceTempX(jk) - FXij
          forceTempY(jk) = forceTempY(jk) - FYij
          forceTempZ(jk) = forceTempZ(jk) - FZij
        end if
      end do loop1
      ! Include intramolecular interaction if need
      if (this%potintra15 .or. this%potintra14) then
        RXij = RXi - RX2(i)
        RYij = RYi - RY2(i)
        RZij = RZi - RZ2(i)
        PXij = PXi - PX2(i)
        PYij = PYi - PY2(i)
        PZij = PZi - PZ2(i)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
        RijInv = sqrt( RijSquaredInv )
        eX = RXij * RijInv
        eY = RYij * RijInv
        eZ = RZij * RijInv
        CosTheta  = OXi * ex + OYi * eY + OZi * eZ                              ! -cos(alpha) bei Price
        CosTheta3 = 3._RK * CosTheta
        Epsilon1 = Epsilon * RijSquaredInv * coeff
        Epsilon2 = Epsilon1 * RijInv
        EPotLocal  = EPotLocal - Epsilon1 * CosTheta                           ! Uebereinstimmumg mit Price
        EPotLocalIntra  = EPotLocalIntra - Epsilon1 * CosTheta
        FXij = Epsilon2 * ( OXi - CosTheta3 * eX )                      ! F1 bei Price, 1-4 Coeff is included in Epsilon1 -> Epsilon2
        FYij = Epsilon2 * ( OYi - CosTheta3 * eY )
        FZij = Epsilon2 * ( OZi - CosTheta3 * eZ )
        VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F1*(-R_COM_Price); stimmt so
        VirialLocalIntra = VirialLocalIntra + (FXij * PXij + FYij * PYij + FZij * PZij)
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local - Epsilon1 * CosTheta*(8._RK*sitecorr*sitecorr-2._RK*Plen2*RijSquaredInv)*Third*Third          !xxxx4  DC T
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        TXi = TXi + Epsilon1 * eX   ! Uebereinstimmumg mit Price; Rest bei Atom2Mol in Component
        TYi = TYi + Epsilon1 * eY   ! Reaktionsfeldbeitrag in Interaction
        TZi = TZi + Epsilon1 * eZ
        forceTempX(i) = forceTempX(i) - FXij
        forceTempY(i) = forceTempY(i) - FYij
        forceTempZ(i) = forceTempZ(i) - FZij
      end if
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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotDC_Force

!==============================================================!
!  Subroutine TPotDC_Force_Trans                               !
!==============================================================!

  
  subroutine TPotDC_Force_Trans( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

!Hier noch die Transportgrößen rein!

    implicit none

    ! Declare arguments
    type(TPotDipoleCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra_Nonbonded
    real(RK), intent(in out) :: VirialIntra
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
    real(RK)          :: EPotLocalIntra, ViriallocalIntra
    real(RK)          :: EPotLocalInter, ViriallocalInter
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, jk, unit
    logical           :: SameComponent
    real(RK)          :: coeff
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
    real(RK), pointer, contiguous :: VSx(:), VSy(:), VSz(:)
    real(RK), pointer, contiguous :: VSux(:),VSuy(:),VSuz(:)
    real(RK), pointer, contiguous :: VBx(:), VBy(:), VBz(:)
    real(RK), pointer, contiguous :: Cx(:) , Cy(:) , Cz(:)
    real(RK), pointer, contiguous :: tux(:) , tuy(:) , tuz(:)
    real(RK), pointer, contiguous :: tlx(:) , tly(:) , tlz(:)
    real(RK), pointer, contiguous :: tdx(:) , tdy(:) , tdz(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
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
    real(RK)          :: UU, Uxi,  Uyi, Uzi
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33
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
    SameComponent = this%SameComponent
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

#if  TRANS == 1
    !TRANSPORT_start
    VSx => this%Site1%vsDx
    VSy => this%Site1%vsDy
    VSz => this%Site1%vsDz
    VBx => this%Site1%vbDx
    VBy => this%Site1%vbDy
    VBz => this%Site1%vbDz
    VSux=> this%Site1%vsuDx
    VSuy=> this%Site1%vsuDy
    VSuz=> this%Site1%vsuDz
    Cx  => this%Site1%cDx
    Cy  => this%Site1%cDy
    Cz  => this%Site1%cDz
    tux => this%Site1%tuDx
    tuy => this%Site1%tuDy
    tuz => this%Site1%tuDz
    tlx => this%Site1%tlDx
    tly => this%Site1%tlDy
    tlz => this%Site1%tlDz
    tdx => this%Site1%tdDx
    tdy => this%Site1%tdDy
    tdz => this%Site1%tdDz
    q1  => this%Site1%Qm0r(:, 1, 1)
    q2  => this%Site1%Qm0r(:, 2, 1)
    q3  => this%Site1%Qm0r(:, 3, 1)
    q4  => this%Site1%Qm0r(:, 4, 1)
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
        VSxi= VSx(i)
        VSyi= VSy(i)
        VSzi= VSz(i)
        VBxi= VBx(i)
        VByi= VBy(i)
        VBzi= VBz(i)
        VSuxi= VSux(i)
        VSuyi= VSuy(i)
        VSuzi= VSuz(i)
        Cxi = Cx(i)
        Cyi = Cy(i)
        Czi = Cz(i)
        tuxi = tux(i)
        tuyi = tuy(i)
        tuzi = tuz(i)
        tlxi = tlx(i)
        tlyi = tly(i)
        tlzi = tlz(i)
        tdxi = tdx(i)
        tdyi = tdy(i)
        tdzi = tdz(i)
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

      unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:do k = 1, this%NInCutoff(unit)
        j = this%CutoffPartner(k, unit)
        if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
          if (mod(j,nu2)==0) then
            jk = INT(j/nu2)
          else
            jk = INT(j/nu2)+1
          end if
          RXij = RXi - RX2(jk)
          RYij = RYi - RY2(jk)
          RZij = RZi - RZ2(jk)
          PXij = PXi - PX2(jk)
          PYij = PYi - PY2(jk)
          PZij = PZi - PZ2(jk)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
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
        EPotLocalInter  = EPotLocalInter + EPotLocal1
        FXij = Epsilon2 * ( OXi - CosTheta3 * eX )                              ! F1 bei Price
        FYij = Epsilon2 * ( OYi - CosTheta3 * eY )
        FZij = Epsilon2 * ( OZi - CosTheta3 * eZ )
        VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F1*(-R_COM_Price); stimmt so
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
        VirialLocalInter = VirialLocalInter + (FXij * PXij + FYij * PYij + FZij * PZij)
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(8._RK*sitecorr*sitecorr-2._RK*Plen2*RijSquaredInv)*Third*Third          !xxxx4  DC T
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        TXi = TXi + Epsilon1 * eX   ! Uebereinstimmumg mit Price; Rest bei Atom2Mol in Component
        TYi = TYi + Epsilon1 * eY   ! Reaktionsfeldbeitrag in Interaction
        TZi = TZi + Epsilon1 * eZ

        forceTempX(jk) = forceTempX(jk) - FXij
        forceTempY(jk) = forceTempY(jk) - FYij
        forceTempZ(jk) = forceTempZ(jk) - FZij
#if TRANS==1
          !TRANSPORT_start
          VSxi   = VSxi + FXij * PYij
          VSyi   = VSyi + FXij * PZij
          VSzi   = VSzi + FYij * PZij
          VBxi   = VBxi + FXij * PXij
          VByi   = VByi + FYij * PYij
          VBzi   = VBzi + FZij * PZij
          VSuxi  = VSuxi+ FYij * PXij
          VSuyi  = VSuyi+ FZij * PXij
          VSuzi  = VSuzi+ FZij * PYij
          UU     = EpotLocal1  !What is EPotLocal1?
          Uxi     = UU * eX
          Uyi     = UU * eY
          Uzi     = UU * eZ   
          Cxi    = Cxi  + Uxi
          Cyi    = Cyi  + Uyi
          Czi    = Czi  + Uzi
          FTXi = Epsilon1 * eX       ! Chequear                    
          FTYi = Epsilon1 * eY       !Chequear
          FTZi = Epsilon1 * eZ       !Chequear
          txii   = OYi * FTZi - OZi * FTYi   ! Chequear 
          tyii   = OZi * FTXi - OXi * FTZi   ! Chequear 
          tzii   = OXi * FTYi - OYi * FTXi   ! Chequear 
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
          !TRANSPORT_END
#endif
        end if
      end do loop1
      ! Include intramolecular interaction if need
      if (this%potintra15 .or. this%potintra14) then
        RXij = RXi - RX2(i)
        RYij = RYi - RY2(i)
        RZij = RZi - RZ2(i)
        PXij = PXi - PX2(i)
        PYij = PYi - PY2(i)
        PZij = PZi - PZ2(i)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
        RijInv = sqrt( RijSquaredInv )
        eX = RXij * RijInv
        eY = RYij * RijInv
        eZ = RZij * RijInv
        CosTheta  = OXi * ex + OYi * eY + OZi * eZ                              ! -cos(alpha) bei Price
        CosTheta3 = 3._RK * CosTheta
        Epsilon1 = Epsilon * RijSquaredInv * coeff
        Epsilon2 = Epsilon1 * RijInv
        EPotLocal1 = - Epsilon1 * CosTheta
        EPotLocal  = EPotLocal + EPotLocal1                          ! Uebereinstimmumg mit Price
        EPotLocalIntra  = EPotLocalIntra + EPotLocal1
        FXij = Epsilon2 * ( OXi - CosTheta3 * eX )                      ! F1 bei Price, 1-4 Coeff is included in Epsilon1 -> Epsilon2
        FYij = Epsilon2 * ( OYi - CosTheta3 * eY )
        FZij = Epsilon2 * ( OZi - CosTheta3 * eZ )
        VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)     ! F1*(-R_COM_Price); stimmt so
        VirialLocalIntra = VirialLocalIntra + (FXij * PXij + FYij * PYij + FZij * PZij)
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(8._RK*sitecorr*sitecorr-2._RK*Plen2*RijSquaredInv)*Third*Third          !xxxx4  DC T
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        TXi = TXi + Epsilon1 * eX   ! Uebereinstimmumg mit Price; Rest bei Atom2Mol in Component
        TYi = TYi + Epsilon1 * eY   ! Reaktionsfeldbeitrag in Interaction
        TZi = TZi + Epsilon1 * eZ
        forceTempX(i) = forceTempX(i) - FXij
        forceTempY(i) = forceTempY(i) - FYij
        forceTempZ(i) = forceTempZ(i) - FZij

      end if
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
      TX1(i) = TXi
      TY1(i) = TYi
      TZ1(i) = TZi
#if  TRANS == 1
        !TRANSPORT_start
        VSx(i) = VSxi
        VSy(i) = VSyi
        VSz(i) = VSzi
        VBx(i) = VBxi
        VBy(i) = VByi
        VBz(i) = VBzi
        VSux(i)= VSuxi
        VSuy(i)= VSuyi
        VSuz(i)= VSuzi
        Cx(i)  = Cxi
        Cy(i)  = Cyi
        Cz(i)  = Czi
        tux(i) = tuxi
        tuy(i) = tuyi
        tuz(i) = tuzi
        tlx(i) = tlxi
        tly(i) = tlyi
        tlz(i) = tlzi
        tdx(i) = tdxi
        tdy(i) = tdyi
        tdz(i) = tdzi
        !TRANSPORT_END
#endif
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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
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
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
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
    integer           :: nu1, nu2, unit, i0, jk
#if ARCH == 3
    logical           :: hit
#endif

    ! Assign local variables
    i1 = this%Site1%NTest
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared
    nu1 = this%NUnit1
    nu2 = this%NUnit2
#if MPI_VER > 0
    i0 = this%Site1%NTest0
    i1 = this%Site1%NTest2
#else
    i0 = 1
    i1 = this%Site1%NTest
#endif

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
    do i = i0, i1
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

#if ARCH == 3
        hit = .false.
#endif
!CDIR NODEP

      unit = nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then  ! choose only units, to which our Site2 correspond
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)   ! number of molecule, to which this unit correspond
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            PXij = (PXij - anint( PXij )) * BoxLength
            PYij = (PYij - anint( PYij )) * BoxLength
            PZij = (PZij - anint( PZij )) * BoxLength
            RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
            if( RijSquared <= RShieldSquared ) hit = .true.
#else
            if( RijSquared <= RShieldSquared ) then
              EPotLocal = 1E33_RK
              exit loop1
            end if
#endif
            RijSquaredInv = 1._RK / RijSquared
            RijInv = sqrt( RijSquaredInv )
            eX = RXij * RijInv
            eY = RYij * RijInv
            eZ = RZij * RijInv
            CosTheta  = OXi * ex + OYi * eY + OZi * eZ
            EPotLocal = EPotLocal - Epsilon * RijSquaredInv * CosTheta
          end if
        end do loop1

#if ARCH == 3
        if( .not. hit ) then
          EPotTest(i) = EPotTest(i) + EPotLocal
        else
          EPotTest(i) = EPotTest(i) + 1E33_RK
        endif
#else
        EPotTest(i) = EPotTest(i) + EPotLocal
#endif
    end do

  end subroutine TPotDC_ChemicalPotential



!==============================================================!
!  Subroutine TPotDC_Energy                                    !
!==============================================================!

  subroutine TPotDC_Energy( this, np, nu, F, E, EIntra, BoxLength, CompIdent )

    implicit none

    ! Declare arguments
    type(TPotDipoleCharge) :: this
    integer, intent(in)    :: np
    integer, intent(in)      :: nu
    real(RK), intent(in out) :: F(3,nu)
    real(RK), intent(in out) :: E
    real(RK), intent(in out) :: EIntra
    real(RK), intent(in)   :: BoxLength
    logical, intent(in)      :: CompIdent

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1
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
    real(RK)          :: EPot, EIntra1, EPotLocal, tempF(3,nu)
    integer           :: j, k, nu2, jk, unit
    real(RK)          :: coeff

    ! Assign local variables
    Epsilon = this%Epsilon
    RShieldSquared = this%RShieldSquared
    nu2 = this%NUnit2
    coeff = 1._RK
    if (this%potintra14) coeff = this%ScaleEl14
    EPot = 0._RK
    EIntra1 = 0._RK

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
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ

    ! Loop over molecules
    RXi = RX1(np)
    RYi = RY1(np)
    RZi = RZ1(np)
    PXi = PX1(np)
    PYi = PY1(np)
    PZi = PZ1(np)
    OXi = OX1(np)
    OYi = OY1(np)
    OZi = OZ1(np)

    unit=this%NUnit1*(np-1)+this%Site1%UnitNumber

    do k = 1, this%NInCutoff(unit)
      j = this%CutoffPartner(k, unit)
      if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
        if (mod(j,nu2)==0) then
          jk = INT(j/nu2)
        else
          jk = INT(j/nu2)+1
        end if
        RXij = RXi - RX2(jk)
        RYij = RYi - RY2(jk)
        RZij = RZi - RZ2(jk)
        PXij = PXi - PX2(jk)
        PYij = PYi - PY2(jk)
        PZij = PZi - PZ2(jk)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        RijSquared = RXij**2 + RYij**2 + RZij**2

        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
        else
          RijSquaredInv = 1._RK / RijSquared
          RijInv = sqrt( RijSquaredInv )
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosTheta = OXi * ex + OYi * eY + OZi * eZ
          Epsilon1 = Epsilon * RijSquaredInv
          EPotLocal = - Epsilon1 * CosTheta
        end if
        EPot = EPot + EPotLocal
      end if
    end do
    ! Include intramolecular interaction if need
    if (this%potintra14 .or. this%potintra15) then
      RXij = RXi - RX2(np)
      RYij = RYi - RY2(np)
      RZij = RZi - RZ2(np)
      PXij = PXi - PX2(np)
      PYij = PYi - PY2(np)
      PZij = PZi - PZ2(np)
      RXij = (RXij - anint( PXij )) * BoxLength
      RYij = (RYij - anint( PYij )) * BoxLength
      RZij = (RZij - anint( PZij )) * BoxLength
      RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
      RijInv = sqrt( RijSquaredInv )
      eX = RXij * RijInv
      eY = RYij * RijInv
      eZ = RZij * RijInv
      CosTheta  = OXi * ex + OYi * eY + OZi * eZ
      EPotLocal = - Epsilon * RijSquaredInv * coeff * CosTheta
      EIntra1  = EIntra1 + EPotLocal
    end if

    F(:,:) = F(:,:) + tempF(:,:)
    E = E + EPot + EIntra1
    EIntra = EIntra + EIntra1

  end subroutine TPotDC_Energy


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

    ! Declare local variables
    integer :: k, ende

    ! Construct potential
    this%Site1 => Molecule1%SiteDipole(j1)
    this%Site2 => Molecule2%SiteDipole(j2)
    this%NUnit1 => Molecule1%NUnit
    this%NUnit2 => Molecule2%NUnit
    this%SameComponent = i1 == i2
    this%Epsilon = this%Site1%D * this%Site2%D
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2
    this%RFConstant = this%Epsilon / RCutoff**3 * (RFEpsilon - 1._RK) / (2._RK * RFEpsilon + 1._RK)

    ! if this potential is intra
    if (this%SameComponent .and. Molecule1%hasIntraLJEl) then
      ende = size(Molecule1%IntDD15(:,1))
      do k=1, ende
        if (Molecule1%IntDD15(k,1)==this%Site1%SiteId .and. Molecule1%IntDD15(k,2)==this%Site2%SiteId) then
          this%potintra15 = .true.
        end if
      end do
      if (LJEL14 .and. .not. this%potintra15) then
        ende = size(Molecule1%IntDD14(:,1))
        do k=1, ende
          if (Molecule1%IntDD14(k,1)==this%Site1%SiteId .and. Molecule1%IntDD14(k,2)==this%Site2%SiteId) then
            this%potintra14=.true.
            this%ScaleEl14 = Molecule1%ScaleDD14(k)
          end if
        end do
      end if
    else
      this%potintra15 = .false.
      this%potintra14 = .false.
    end if

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

  subroutine TPotDD_Force( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleDipole)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra_Nonbonded
    real(RK), intent(in out) :: VirialIntra
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
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
    integer           :: nu1, nu2, jk, unit
    real(RK)          :: coeff
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
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

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

        unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            PXij = (PXij - anint( PXij )) * BoxLength
            PYij = (PYij - anint( PYij )) * BoxLength
            PZij = (PZij - anint( PZij )) * BoxLength
            OXj = OX2(jk)
            OYj = OY2(jk)
            OZj = OZ2(jk)

#if ARCH == 3
            RijInv = rsqrt( RXij**2 + RYij**2 + RZij**2 )
#else
            RijInv = 1._RK / sqrt( RXij**2 + RYij**2 + RZij**2 )
#endif
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
            EPotLocalInter = EPotLocalInter + Rij3Inv * Tmp
            FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                       - (eX * CosThetaj - OXj) * CosThetai)
            FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                       - (eY * CosThetaj - OYj) * CosThetai)
            FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                       - (eZ * CosThetaj - OZj) * CosThetai)
            VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
            VirialLocalInter = VirialLocalInter + (FXij * PXij + FYij * PYij + FZij * PZij)
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
            d2EpotdV2Local = d2EpotdV2Local + Rij3Inv*Tmp*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijInv*RijInv)*Third*Third         !xxxx5   DD

            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            forceTempX(jk) = forceTempX(jk) - FXij
            forceTempY(jk) = forceTempY(jk) - FYij
            forceTempZ(jk) = forceTempZ(jk) - FZij
            TXi    = TXi    + Rij3Inv * (eX * CosThetaj3 - OXj)
            TYi    = TYi    + Rij3Inv * (eY * CosThetaj3 - OYj)
            TZi    = TZi    + Rij3Inv * (eZ * CosThetaj3 - OZj)
          
            momTempX(jk) = momTempX(jk) + Rij3Inv * (eX * CosThetai3 - OXi)  
            momTempY(jk) = momTempY(jk) + Rij3Inv * (eY * CosThetai3 - OYi)  
            momTempZ(jk) = momTempZ(jk) + Rij3Inv * (eZ * CosThetai3 - OZi)  

          end if
        end do loop1
        ! Include intramolecular interactions if need
        if (this%potintra15 .or. this%potintra14) then
          RXij = RXi - RX2(i)
          RYij = RYi - RY2(i)
          RZij = RZi - RZ2(i)
          PXij = PXi - PX2(i)
          PYij = PYi - PY2(i)
          PZij = PZi - PZ2(i)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          OXj = OX2(i)
          OYj = OY2(i)
          OZj = OZ2(i)
#if ARCH == 3
          RijInv = rsqrt( RXij**2 + RYij**2 + RZij**2 )
#else
          RijInv = 1._RK / sqrt( RXij**2 + RYij**2 + RZij**2 )
#endif
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          CosThetai3 = 3._RK * CosThetai
          CosThetaj3 = 3._RK * CosThetaj
          Tmp = CosGammaij - CosThetai * CosThetaj3
          Rij3Inv = Epsilon * RijInv**3 * coeff
          Rij4Inv3 = 3._RK * Rij3Inv * RijInv
          EPotLocal = EPotLocal +  Rij3Inv * Tmp
          EPotLocalIntra = EPotLocalIntra + Rij3Inv * Tmp
          FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                     - (eX * CosThetaj - OXj) * CosThetai)
          FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                     - (eY * CosThetaj - OYj) * CosThetai)
          FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                     - (eZ * CosThetaj - OZj) * CosThetai)
          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
          VirialLocalIntra = VirialLocalIntra + (FXij * PXij + FYij * PYij + FZij * PZij)
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + Rij3Inv*Tmp*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijInv*RijInv)*Third*Third         !xxxx5   DD

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          forceTempX(i) = forceTempX(i) - FXij
          forceTempY(i) = forceTempY(i) - FYij
          forceTempZ(i) = forceTempZ(i) - FZij
          TXi    = TXi    + Rij3Inv * (eX * CosThetaj3 - OXj)
          TYi    = TYi    + Rij3Inv * (eY * CosThetaj3 - OYj)
          TZi    = TZi    + Rij3Inv * (eZ * CosThetaj3 - OZj)
          
          momTempX(i) = momTempX(i) + Rij3Inv * (eX * CosThetai3 - OXi) 
          momTempY(i) = momTempY(i) + Rij3Inv * (eY * CosThetai3 - OYi)
          momTempZ(i) = momTempZ(i) + Rij3Inv * (eZ * CosThetai3 - OZi)

        end if
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
        unit=nu1*(i-1)+this%Site1%UnitNumber
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
          RijSquared = RXij**2 + RYij**2 + RZij**2
          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
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
          d2EpotdV2Local = d2EpotdV2Local + Rij3Inv*Tmp*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijInv*RijInv)*Third*Third         !xxxx5   DD ss

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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if ( IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotDD_Force



!==============================================================!
!  Subroutine TPotDD_Force_Trans                               !
!==============================================================!

  subroutine TPotDD_Force_Trans( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleDipole)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra_Nonbonded
    real(RK), intent(in out) :: VirialIntra
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
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
    integer           :: nu1, nu2, jk, unit
    real(RK)          :: coeff

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
    real(RK), pointer, contiguous :: VSx(:), VSy(:), VSz(:)
    real(RK), pointer, contiguous :: VSux(:),VSuy(:),VSuz(:)
    real(RK), pointer, contiguous :: VBx(:), VBy(:), VBz(:)
    real(RK), pointer, contiguous :: Cx(:) , Cy(:) , Cz(:)
    real(RK), pointer, contiguous :: tux(:) , tuy(:) , tuz(:)
    real(RK), pointer, contiguous :: tlx(:) , tly(:) , tlz(:)
    real(RK), pointer, contiguous :: tdx(:) , tdy(:) , tdz(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
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
    real(RK)          :: UU, Uxi,  Uyi, Uzi
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33
    !TRANSPORT_END
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
    RFConstant2 = 2._RK * this%RFConstant
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

#if  TRANS == 1
    !TRANSPORT_start
    VSx => this%Site1%vsDx
    VSy => this%Site1%vsDy
    VSz => this%Site1%vsDz
    VBx => this%Site1%vbDx
    VBy => this%Site1%vbDy
    VBz => this%Site1%vbDz
    VSux=> this%Site1%vsuDx
    VSuy=> this%Site1%vsuDy
    VSuz=> this%Site1%vsuDz
    Cx  => this%Site1%cDx
    Cy  => this%Site1%cDy
    Cz  => this%Site1%cDz
    tux => this%Site1%tuDx
    tuy => this%Site1%tuDy
    tuz => this%Site1%tuDz
    tlx => this%Site1%tlDx
    tly => this%Site1%tlDy
    tlz => this%Site1%tlDz
    tdx => this%Site1%tdDx
    tdy => this%Site1%tdDy
    tdz => this%Site1%tdDz
    q1  => this%Site1%Qm0r(:, 1, 1)
    q2  => this%Site1%Qm0r(:, 2, 1)
    q3  => this%Site1%Qm0r(:, 3, 1)
    q4  => this%Site1%Qm0r(:, 4, 1)
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
        VSxi= VSx(i)
        VSyi= VSy(i)
        VSzi= VSz(i)
        VBxi= VBx(i)
        VByi= VBy(i)
        VBzi= VBz(i)
        VSuxi= VSux(i)
        VSuyi= VSuy(i)
        VSuzi= VSuz(i)
        Cxi = Cx(i)
        Cyi = Cy(i)
        Czi = Cz(i)
        tuxi = tux(i)
        tuyi = tuy(i)
        tuzi = tuz(i)
        tlxi = tlx(i)
        tlyi = tly(i)
        tlzi = tlz(i)
        tdxi = tdx(i)
        tdyi = tdy(i)
        tdzi = tdz(i)
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

        unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            PXij = (PXij - anint( PXij )) * BoxLength
            PYij = (PYij - anint( PYij )) * BoxLength
            PZij = (PZij - anint( PZij )) * BoxLength
            OXj = OX2(jk)
            OYj = OY2(jk)
            OZj = OZ2(jk)
#if ARCH == 3
            RijInv = rsqrt( RXij**2 + RYij**2 + RZij**2 )
#else
            RijInv = 1._RK / sqrt( RXij**2 + RYij**2 + RZij**2 )
#endif
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
            EPotLocalInter = EPotLocalInter + Rij3Inv * Tmp
            FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                       - (eX * CosThetaj - OXj) * CosThetai)
            FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                       - (eY * CosThetaj - OYj) * CosThetai)
            FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                       - (eZ * CosThetaj - OZj) * CosThetai)
            VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
          VirialLocalInter = VirialLocalInter + (FXij * PXij + FYij * PYij + FZij * PZij)
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + Rij3Inv * Tmp*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijInv*RijInv)*Third*Third         !xxxx5   DD T

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          forceTempX(jk) = forceTempX(jk) - FXij
          forceTempY(jk) = forceTempY(jk) - FYij
          forceTempZ(jk) = forceTempZ(jk) - FZij
          TXi    = TXi    + Rij3Inv * (eX * CosThetaj3 - OXj)
          TYi    = TYi    + Rij3Inv * (eY * CosThetaj3 - OYj)
          TZi    = TZi    + Rij3Inv * (eZ * CosThetaj3 - OZj)
          
          momTempX(jk) = momTempX(jk) + Rij3Inv * (eX * CosThetai3 - OXi)  
          momTempY(jk) = momTempY(jk) + Rij3Inv * (eY * CosThetai3 - OYi)  
          momTempZ(jk) = momTempZ(jk) + Rij3Inv * (eZ * CosThetai3 - OZi)  
#if TRANS==1
          !TRANSPORT_start
          VSxi   = VSxi + FXij * PYij
          VSyi   = VSyi + FXij * PZij
          VSzi   = VSzi + FYij * PZij
          VBxi   = VBxi + FXij * PXij
          VByi   = VByi + FYij * PYij
          VBzi   = VBzi + FZij * PZij
          VSuxi  = VSuxi+ FYij * PXij
          VSuyi  = VSuyi+ FZij * PXij
          VSuzi  = VSuzi+ FZij * PYij
          UU        = Rij3Inv * Tmp - RFConstant2 * CosGammaij
          Uxi       = UU * eX
          Uyi       = UU * eY
          Uzi       = UU * eZ
          FTXi   = Rij3Inv * (eX * CosThetaj3 - OXj) + OXj * RFConstant2
          FTYi   = Rij3Inv * (eY * CosThetaj3 - OYj) + OYj * RFConstant2
          FTZi   = Rij3Inv * (eZ * CosThetaj3 - OZj) + OZj * RFConstant2
          Cxi    = Cxi  + Uxi
          Cyi    = Cyi  + Uyi
          Czi    = Czi  + Uzi
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
          !TRANSPORT_END
#endif
          end if
        end do loop1
        ! Include intramolecular interactions if need
        if (this%potintra15 .or. this%potintra14) then
          RXij = RXi - RX2(i)
          RYij = RYi - RY2(i)
          RZij = RZi - RZ2(i)
          PXij = PXi - PX2(i)
          PYij = PYi - PY2(i)
          PZij = PZi - PZ2(i)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          OXj = OX2(i)
          OYj = OY2(i)
          OZj = OZ2(i)
#if ARCH == 3
          RijInv = rsqrt( RXij**2 + RYij**2 + RZij**2 )
#else
          RijInv = 1._RK / sqrt( RXij**2 + RYij**2 + RZij**2 )
#endif
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          CosThetai3 = 3._RK * CosThetai
          CosThetaj3 = 3._RK * CosThetaj
          Tmp = CosGammaij - CosThetai * CosThetaj3
          Rij3Inv = Epsilon * RijInv**3 * coeff
          Rij4Inv3 = 3._RK * Rij3Inv * RijInv
          EPotLocal = EPotLocal + Rij3Inv * Tmp
          EPotLocalIntra = EPotLocalIntra + Rij3Inv * Tmp
          FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                     - (eX * CosThetaj - OXj) * CosThetai)
          FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                     - (eY * CosThetaj - OYj) * CosThetai)
          FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                     - (eZ * CosThetaj - OZj) * CosThetai)
          VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
          VirialLocalIntra = VirialLocalIntra + FXij * PXij + FYij * PYij + FZij * PZij
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + Rij3Inv * Tmp*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijInv*RijInv)*Third*Third         !xxxx5   DD T
          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          forceTempX(i) = forceTempX(i) - FXij
          forceTempY(i) = forceTempY(i) - FYij
          forceTempZ(i) = forceTempZ(i) - FZij
          TXi    = TXi    + Rij3Inv * (eX * CosThetaj3 - OXj)
          TYi    = TYi    + Rij3Inv * (eY * CosThetaj3 - OYj)
          TZi    = TZi    + Rij3Inv * (eZ * CosThetaj3 - OZj)
          momTempX(i) = momTempX(i) + Rij3Inv * (eX * CosThetai3 - OXi)
          momTempY(i) = momTempY(i) + Rij3Inv * (eY * CosThetai3 - OYi)
          momTempZ(i) = momTempZ(i) + Rij3Inv * (eZ * CosThetai3 - OZi)

        end if

        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
#if  TRANS == 1
        !TRANSPORT_start
        VSx(i) = VSxi
        VSy(i) = VSyi
        VSz(i) = VSzi
        VBx(i) = VBxi
        VBy(i) = VByi
        VBz(i) = VBzi
        VSux(i)= VSuxi
        VSuy(i)= VSuyi
        VSuz(i)= VSuzi
        Cx(i)  = Cxi
        Cy(i)  = Cyi
        Cz(i)  = Czi
        tux(i) = tuxi
        tuy(i) = tuyi
        tuz(i) = tuzi
        tlx(i) = tlxi
        tly(i) = tlyi
        tlz(i) = tlzi
        tdx(i) = tdxi
        tdy(i) = tdyi
        tdz(i) = tdzi
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
        unit=nu1*(i-1)+this%Site1%UnitNumber
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
          RijSquared = RXij**2 + RYij**2 + RZij**2
          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
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
          d2EpotdV2Local = d2EpotdV2Local + Rij3Inv*Tmp*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijInv*RijInv)*Third*Third         !xxxx5   DD ss T

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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if ( IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
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
    real(RK), pointer:: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
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
    integer           :: nu1, nu2, unit, jk, i0
#if ARCH == 3
    logical           :: hit
#endif

    ! Assign local variables
    i1 = this%Site1%NTest
    j1 = this%Site2%NPart
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared
    RFConstant2 = 2._RK * this%RFConstant
    nu1 = this%NUnit1
    nu2 = this%NUnit2
#if MPI_VER > 0
    i0 = this%Site1%NTest0
    i1 = this%Site1%NTest2
#else
    i0 = 1
    i1 = this%Site1%NTest
#endif

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
!$OMP PRIVATE (EPotLocal,i,i0,i1,j,k) 
    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
!$OMP DO
      do i = i0, i1
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
#if ARCH == 3
        hit = .false.
#endif

!CDIR NODEP

        unit = nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then  ! choose only units, to which our Site2 correspond
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)   ! number of molecule, to which this unit correspond
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
            if( RijSquared <= RShieldSquared ) hit = .true.
#else
            if( RijSquared <= RShieldSquared ) then
              EPotLocal = 1E33_RK
              exit loop1
            end if
#endif
            OXj = OX2(jk)
            OYj = OY2(jk)
            OZj = OZ2(jk)
#if ARCH == 3
            RijInv = rsqrt( RijSquared )
#else
            RijInv = 1._RK / sqrt( RijSquared )
#endif
            eX = RXij * RijInv
            eY = RYij * RijInv
            eZ = RZij * RijInv
            CosThetai = OXi * eX + OYi * eY + OZi * eZ
            CosThetaj = OXj * eX + OYj * eY + OZj * eZ
            CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
            Tmp = CosGammaij - 3._RK * CosThetai * CosThetaj
            Rij3Inv = Epsilon * RijInv**3
            EPotLocal = EPotLocal + Rij3Inv * Tmp
          end if
        end do loop1

#if ARCH == 3
        if( .not. hit ) then
          EPotTest(i) = EPotTest(i) + EPotLocal
        else
          EPotTest(i) = EPotTest(i) + 1E33_RK
        endif
#else
        EPotTest(i) = EPotTest(i) + EPotLocal
#endif
      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over test particles
!$OMP DO
      do i = i0, i1
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        EPotLocal = 0._RK
#if ARCH == 3
        hit = .false.
#endif

!CDIR NODEP
loop2:  do j = 1, j1
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij**2 + RYij**2 + RZij**2
#if ARCH == 3
          if( RijSquared <= RShieldSquared ) hit = .true.
#else
          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop2
          end if
#endif
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)
#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
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
#if ARCH == 3
        if( .not. hit ) then
          EPotTest(i) = EPotTest(i) + EPotLocal
        else
          EPotTest(i) = EPotTest(i) + 1E33_RK
        endif
#else
        EPotTest(i) = EPotTest(i) + EPotLocal
#endif
      end do
!$OMP END DO
    end if
!$OMP END PARALLEL
  end subroutine TPotDD_ChemicalPotential



!==============================================================!
!  Subroutine TPotDD_Energy                                    !
!==============================================================!

  subroutine TPotDD_Energy( this, np, nu, F, E, EIntra, BoxLength, CompIdent )

    implicit none

    ! Declare arguments
    type(TPotDipoleDipole) :: this
    integer, intent(in)    :: np
    integer, intent(in)      :: nu
    real(RK), intent(in out) :: F(3,nu)
    real(RK), intent(in out) :: E
    real(RK), intent(in out) :: EIntra
    real(RK), intent(in)   :: BoxLength
    logical, intent(in)      :: CompIdent

    ! Declare local variables
    real(RK), pointer, contiguous :: RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: PX2(:), PY2(:), PZ2(:)
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
    real(RK)          :: E1, EIntra1, ELocal, tempF(3,nu)
    integer           :: j, k, nu2, jk, unit
    real(RK)          :: coeff

    ! Assign local variables
    nu2 = this%NUnit2
    coeff = 1._RK
    if (this%potintra14) coeff = this%ScaleEl14
    E1 = 0._RK
    EIntra1 = 0._RK

    ! Assign pointers
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

    ! Loop over molecules
    RXi = this%Site1%RX(np)
    RYi = this%Site1%RY(np)
    RZi = this%Site1%RZ(np)
    OXi = this%Site1%OX(np)
    OYi = this%Site1%OY(np)
    OZi = this%Site1%OZ(np)
    PXi = this%Site1%PX(np)
    PYi = this%Site1%PY(np)
    PZi = this%Site1%PZ(np)

    unit=this%NUnit1*(np-1)+this%Site1%UnitNumber

    do k = 1, this%NInCutoff(unit)
      j = this%CutoffPartner(k, unit)
      if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
        if (mod(j,nu2)==0) then
          jk = INT(j/nu2)
        else
          jk = INT(j/nu2)+1
        end if
        RXij = RXi - RX2(jk)
        RYij = RYi - RY2(jk)
        RZij = RZi - RZ2(jk)
        PXij = PXi - PX2(jk)
        PYij = PYi - PY2(jk)
        PZij = PZi - PZ2(jk)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        RijSquared = RXij**2 + RYij**2 + RZij**2

        if( RijSquared <= this%RShieldSquared ) then
          ELocal = 1E33_RK
        else
          OXj = OX2(jk)
          OYj = OY2(jk)
          OZj = OZ2(jk)
#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          Tmp = CosGammaij - CosThetai * 3._RK * CosThetaj
          Rij3Inv = this%Epsilon * RijInv**3
          ELocal = Rij3Inv * Tmp
        end if
        E1 = E1 + ELocal
      end if
    end do
    ! Include intramolecular interactions if need
    if (this%potintra15 .or. this%potintra14) then
      RXij = RXi - RX2(np)
      RYij = RYi - RY2(np)
      RZij = RZi - RZ2(np)
      PXij = PXi - PX2(np)
      PYij = PYi - PY2(np)
      PZij = PZi - PZ2(np)
      RXij = (RXij - anint( PXij )) * BoxLength
      RYij = (RYij - anint( PYij )) * BoxLength
      RZij = (RZij - anint( PZij )) * BoxLength
      OXj = OX2(np)
      OYj = OY2(np)
      OZj = OZ2(np)
#if ARCH == 3
      RijInv = rsqrt( RXij**2 + RYij**2 + RZij**2 )
#else
      RijInv = 1._RK / sqrt( RXij**2 + RYij**2 + RZij**2 )
#endif
      eX = RXij * RijInv
      eY = RYij * RijInv
      eZ = RZij * RijInv
      CosThetai = OXi * eX + OYi * eY + OZi * eZ
      CosThetaj = OXj * eX + OYj * eY + OZj * eZ
      CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
      Tmp = CosGammaij - CosThetai * 3._RK * CosThetaj
      Rij3Inv = this%Epsilon * RijInv**3
      ELocal = Rij3Inv * Tmp
      EIntra1 = EIntra1 + ELocal
    end if

    F(:,:) = F(:,:) + tempF(:,:)
    E = E + E1 + EIntra1
    EIntra = EIntra + EIntra1

  end subroutine TPotDD_Energy


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

    ! Declare local variables
    integer :: k, ende

    ! Construct potential
    this%Site1 => Molecule1%SiteDipole(j1)
    this%Site2 => Molecule2%SiteQuadrupole(j2)
    this%NUnit1 => Molecule1%NUnit
    this%NUnit2 => Molecule2%NUnit
    this%SameComponent = i1 == i2
    this%Epsilon = 1.5_RK * this%Site1%D * this%Site2%Q
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

    ! if this potential is intra
    if (this%SameComponent .and. Molecule1%hasIntraLJEl) then
      ende = size(Molecule1%IntDQ15(:,1))
      do k=1, ende
        if (Molecule1%IntDQ15(k,1)==this%Site1%SiteId .and. Molecule1%IntDQ15(k,2)==this%Site2%SiteId) then
          this%potintra15 = .true.
        end if
      end do
      if (LJEL14 .and. .not. this%potintra15) then
        ende = size(Molecule1%IntDQ14(:,1))
        do k=1, ende
          if (Molecule1%IntDQ14(k,1)==this%Site1%SiteId .and. Molecule1%IntDQ14(k,2)==this%Site2%SiteId) then
            this%potintra14=.true.
            this%ScaleEl14 = Molecule1%ScaleDQ14(k)
          end if
        end do
      end if
    else
      this%potintra15 = .false.
      this%potintra14 = .false.
    end if

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

  subroutine TPotDQ_Force( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleQuadrupole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: EPotInter
    real(RK), intent(in out)   :: VirialInter
    real(RK), intent(in out)   :: EPotIntra_Nonbonded
    real(RK), intent(in out)   :: VirialIntra
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
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
    integer           :: nu1, nu2, jk, unit
    real(RK)          :: coeff

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
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter   = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra   = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

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

        unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            PXij = (PXij - anint( PXij )) * BoxLength
            PYij = (PYij - anint( PYij )) * BoxLength
            PZij = (PZij - anint( PZij )) * BoxLength
            OXj = OX2(jk)
            OYj = OY2(jk)
            OZj = OZ2(jk)
            RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
            RijInv = rsqrt( RijSquared )
#else
            RijInv = 1._RK / sqrt( RijSquared )
#endif
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
            EPotLocalInter = EPotLocalInter + EPotLocal1
            dCosThetai = Rij4Inv * (1 - 5._RK * CosThetaj2)
            dCosThetaj = Rij4Inv * (CosGammaij - 10._RK * CosThetai * CosThetaj)
            dCosGammaij = 2._RK * Rij4Inv * CosThetaj
            Tmp = -4._RK * RijInv * EPotLocal1

            FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                      + (eX * CosThetaj - OXj) * dCosThetaj)
            FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                      + (eY * CosThetaj - OYj) * dCosThetaj)
            FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                      + (eZ * CosThetaj - OZj) * dCosThetaj)

            VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
            VirialLocalInter = VirialLocalInter + (FXij * PXij + FYij * PYij + FZij * PZij)
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
            d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 *(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Third*Third    !xxxx6   DQ

            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            forceTempX(jk) = forceTempX(jk) - FXij
            forceTempY(j) = forceTempY(jk) - FYij
            forceTempZ(j) = forceTempZ(jk) - FZij
            TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
            TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
            TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
            momTempX(jk) = momTempX(jk) - eX * dCosThetaj - OXi * dCosGammaij
            momTempY(jk) = momTempY(jk) - eY * dCosThetaj - OYi * dCosGammaij  
            momTempZ(jk) = momTempZ(jk) - eZ * dCosThetaj - OZi * dCosGammaij
          end if

        end do loop1
        ! Include intramolecular interaction if need
        if (this%potintra15 .or. this%potintra14) then
          RXij = RXi - RX2(i)
          RYij = RYi - RY2(i)
          RZij = RZi - RZ2(i)
          PXij = PXi - PX2(i)
          PYij = PYi - PY2(i)
          PZij = PZi - PZ2(i)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          OXj = OX2(i)
          OYj = OY2(i)
          OZj = OZ2(i)
          RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosThetaj2 = CosThetaj**2
          CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal1 = coeff * Rij4Inv * (CosGammaij * CosThetaj - CosThetai * (5._RK * CosThetaj2 - 1))
          EPotLocal = EPotLocal + EPotLocal1
          EPotLocalIntra = EPotLocalIntra + EPotLocal1
          dCosThetai = Rij4Inv * (1 - 5._RK * CosThetaj2)
          dCosThetaj = Rij4Inv * (CosGammaij - 10._RK * CosThetai * CosThetaj)
          dCosGammaij = 2._RK * Rij4Inv * CosThetaj
          Tmp = -4._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FXij = FXij * coeff
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FYij = FYij * coeff
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)
          FZij = FZij * coeff
          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
          VirialLocalIntra = VirialLocalIntra + (FXij * PXij + FYij * PYij + FZij * PZij)
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 *(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Third*Third    !xxxx6   DQ T
          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(i) = forceTempX(i) - FXij
          forceTempY(i) = forceTempY(i) - FYij
          forceTempZ(i) = forceTempZ(i) - FZij
          TXi    = TXi    - (eX * dCosThetai + OXj * dCosGammaij) * coeff
          TYi    = TYi    - (eY * dCosThetai + OYj * dCosGammaij) * coeff
          TZi    = TZi    - (eZ * dCosThetai + OZj * dCosGammaij) * coeff
          momTempX(i) = momTempX(i) - eX * dCosThetaj - OXi * dCosGammaij * coeff
          momTempY(i) = momTempY(i) - eY * dCosThetaj - OYi * dCosGammaij * coeff
          momTempZ(i) = momTempZ(i) - eZ * dCosThetaj - OZi * dCosGammaij * coeff
        end if

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
          RijSquared = RXij**2 + RYij**2 + RZij**2
          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
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
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 *(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Third*Third    !xxxx6   DQ ss

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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotDQ_Force



!==============================================================!
!  Subroutine TPotDQ_Force_Trans                               !
!==============================================================!

  subroutine TPotDQ_Force_Trans( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleQuadrupole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: EPotInter
    real(RK), intent(in out)   :: VirialInter
    real(RK), intent(in out)   :: EPotIntra_Nonbonded
    real(RK), intent(in out)   :: VirialIntra
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
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
    integer           :: nu1, nu2, jk, unit
    real(RK)          :: coeff

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
    real(RK), pointer, contiguous :: VSx(:), VSy(:), VSz(:)
    real(RK), pointer, contiguous :: VSux(:),VSuy(:),VSuz(:)
    real(RK), pointer, contiguous :: VBx(:), VBy(:), VBz(:)
    real(RK), pointer, contiguous :: Cx(:) , Cy(:) , Cz(:)
    real(RK), pointer, contiguous :: tux(:) , tuy(:) , tuz(:)
    real(RK), pointer, contiguous :: tlx(:) , tly(:) , tlz(:)
    real(RK), pointer, contiguous :: tdx(:) , tdy(:) , tdz(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
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
    real(RK)          :: UU, Uxi,  Uyi, Uzi
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33
    !TRANSPORT_END
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
!$OMP PRIVATE(VSx, VSy, VSz ,VSux,VSuy,VSuz, VBx, VBy, VBz, Cx , Cy , Cz) &
!$OMP PRIVATE( tux , tuy , tuz, tlx , tly , tlz, tdx , tdy , tdz) &
!$OMP PRIVATE( q1, q2, q3, q4, VSxi, VSyi, VSzi, VSuxi,VSuyi,VSuzi) &
!$OMP PRIVATE( VBxi, VByi, VBzi, Cxi,  Cyi,  Czi, tuxi,  tuyi,  tuzi, tlxi,  tlyi,  tlzi) &
!$OMP PRIVATE(  tdxi,  tdyi,  tdzi, txii,  tyii , tzii, txir ,  tyir  , tzir ) &
!$OMP PRIVATE(   Uxi,  Uyi, Uzi, FTXi , FTYi , FTZi) &
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
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter   = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra   = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

#if  TRANS == 1
    !TRANSPORT_start
   ! Conductivity = this%Conductivity
    VSx => this%Site1%vsDx
    VSy => this%Site1%vsDy
    VSz => this%Site1%vsDz
    VBx => this%Site1%vbDx
    VBy => this%Site1%vbDy
    VBz => this%Site1%vbDz
    VSux=> this%Site1%vsuDx
    VSuy=> this%Site1%vsuDy
    VSuz=> this%Site1%vsuDz
    Cx  => this%Site1%cDx
    Cy  => this%Site1%cDy
    Cz  => this%Site1%cDz
    tux => this%Site1%tuDx
    tuy => this%Site1%tuDy
    tuz => this%Site1%tuDz
    tlx => this%Site1%tlDx
    tly => this%Site1%tlDy
    tlz => this%Site1%tlDz
    tdx => this%Site1%tdDx
    tdy => this%Site1%tdDy
    tdz => this%Site1%tdDz
    q1  => this%Site1%Qm0r(:, 1, 1)
    q2  => this%Site1%Qm0r(:, 2, 1)
    q3  => this%Site1%Qm0r(:, 3, 1)
    q4  => this%Site1%Qm0r(:, 4, 1)
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
        VSxi= VSx(i)
        VSyi= VSy(i)
        VSzi= VSz(i)
        VBxi= VBx(i)
        VByi= VBy(i)
        VBzi= VBz(i)
        VSuxi= VSux(i)
        VSuyi= VSuy(i)
        VSuzi= VSuz(i)
        Cxi = Cx(i)
        Cyi = Cy(i)
        Czi = Cz(i)
        tuxi = tux(i)
        tuyi = tuy(i)
        tuzi = tuz(i)
        tlxi = tlx(i)
        tlyi = tly(i)
        tlzi = tlz(i)
        tdxi = tdx(i)
        tdyi = tdy(i)
        tdzi = tdz(i)
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

        unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            PXij = (PXij - anint( PXij )) * BoxLength
            PYij = (PYij - anint( PYij )) * BoxLength
            PZij = (PZij - anint( PZij )) * BoxLength
            OXj = OX2(jk)
            OYj = OY2(jk)
            OZj = OZ2(jk)
            RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
            RijInv = rsqrt( RijSquared )
#else
            RijInv = 1._RK / sqrt( RijSquared )
#endif
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
            EPotLocalInter = EPotLocalInter + EPotLocal1
            dCosThetai = Rij4Inv * (1 - 5._RK * CosThetaj2)
            dCosThetaj = Rij4Inv * (CosGammaij - 10._RK * CosThetai * CosThetaj)
            dCosGammaij = 2._RK * Rij4Inv * CosThetaj
            Tmp = -4._RK * RijInv * EPotLocal1

            FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                      + (eX * CosThetaj - OXj) * dCosThetaj)
            FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                      + (eY * CosThetaj - OYj) * dCosThetaj)
            FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                      + (eZ * CosThetaj - OZj) * dCosThetaj)

            VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
            VirialLocalInter = VirialLocalInter + (FXij * PXij + FYij * PYij + FZij * PZij)
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
            d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 *(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Third*Third    !xxxx6   DQ T

            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij

            forceTempX(jk) = forceTempX(jk) - FXij
            forceTempY(jk) = forceTempY(jk) - FYij
            forceTempZ(jk) = forceTempZ(jk) - FZij

            TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
            TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
            TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij

            momTempX(jk) = momTempX(jk) - eX * dCosThetaj - OXi * dCosGammaij
            momTempY(jk) = momTempY(jk) - eY * dCosThetaj - OYi * dCosGammaij  
            momTempZ(jk) = momTempZ(jk) - eZ * dCosThetaj - OZi * dCosGammaij   
#if TRANS==1
            !TRANSPORT_start
            VSxi   = VSxi + FXij * PYij
            VSyi   = VSyi + FXij * PZij
            VSzi   = VSzi + FYij * PZij
            VBxi   = VBxi + FXij * PXij
            VByi   = VByi + FYij * PYij
            VBzi   = VBzi + FZij * PZij
            VSuxi  = VSuxi+ FYij * PXij
            VSuyi  = VSuyi+ FZij * PXij
            VSuzi  = VSuzi+ FZij * PYij
            UU     =  EpotLocal1
            Uxi    =  UU * eX
            Uyi    =  UU * eY
            Uzi    =  UU * eZ 
            Cxi    = Cxi  + Uxi
            Cyi    = Cyi  + Uyi
            Czi    = Czi  + Uzi          
            FTXi   = - eX * dCosThetai - OXj * dCosGammaij 
            FTYi   = - eY * dCosThetai - OYj * dCosGammaij  
            FTZi   = - eZ * dCosThetaj - OZi * dCosGammaij 
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
            !TRANSPORT_END
#endif
          end if
        end do loop1

        ! Include intramolecular interaction if need
        if (this%potintra15 .or. this%potintra14) then
          RXij = RXi - RX2(i)
          RYij = RYi - RY2(i)
          RZij = RZi - RZ2(i)
          PXij = PXi - PX2(i)
          PYij = PYi - PY2(i)
          PZij = PZi - PZ2(i)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          OXj = OX2(i)
          OYj = OY2(i)
          OZj = OZ2(i)
          RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosThetaj2 = CosThetaj**2
          CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal1 = coeff * Rij4Inv * (CosGammaij * CosThetaj &
&                      - CosThetai * (5._RK * CosThetaj2 - 1))
          EPotLocal = EPotLocal + EPotLocal1
          EPotLocalIntra = EPotLocalIntra + EPotLocal1
          dCosThetai = Rij4Inv * (1 - 5._RK * CosThetaj2)
          dCosThetaj = Rij4Inv * (CosGammaij - 10._RK * CosThetai * CosThetaj)
          dCosGammaij = 2._RK * Rij4Inv * CosThetaj
          Tmp = -4._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FXij = FXij * coeff
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FYij = FYij * coeff
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)
          FZij = FZij * coeff
          VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
          VirialLocalIntra = VirialLocalIntra + FXij * PXij + FYij * PYij + FZij * PZij
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 *(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Third*Third    !xxxx6   DQ ss T

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          forceTempX(i) = forceTempX(i) - FXij
          forceTempY(i) = forceTempY(i) - FYij
          forceTempZ(i) = forceTempZ(i) - FZij
          TXi    = TXi    - (eX * dCosThetai + OXj * dCosGammaij) * coeff
          TYi    = TYi    - (eY * dCosThetai + OYj * dCosGammaij) * coeff
          TZi    = TZi    - (eZ * dCosThetai + OZj * dCosGammaij) * coeff
          momTempX(i) = momTempX(i) - (eX * dCosThetaj + OXi * dCosGammaij) * coeff
          momTempY(i) = momTempY(i) - (eY * dCosThetaj + OYi * dCosGammaij) * coeff
          momTempZ(i) = momTempZ(i) - (eZ * dCosThetaj + OZi * dCosGammaij) * coeff

            end if

        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
#if  TRANS == 1
        !TRANSPORT_start
        VSx(i) = VSxi
        VSy(i) = VSyi
        VSz(i) = VSzi
        VBx(i) = VBxi
        VBy(i) = VByi
        VBz(i) = VBzi
        VSux(i)= VSuxi
        VSuy(i)= VSuyi
        VSuz(i)= VSuzi
        Cx(i)  = Cxi
        Cy(i)  = Cyi
        Cz(i)  = Czi
        tux(i) = tuxi
        tuy(i) = tuyi
        tuz(i) = tuzi
        tlx(i) = tlxi
        tly(i) = tlyi
        tlz(i) = tlzi
        tdx(i) = tdxi
        tdy(i) = tdyi
        tdz(i) = tdzi
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
          RijSquared = RXij**2 + RYij**2 + RZij**2
          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
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
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1 *(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Third*Third    !xxxx6   DQ ss T

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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
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
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
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
    integer           :: nu1, nu2, unit, jk, i0
#if ARCH == 3
    logical           :: hit
#endif

    ! Assign local variables
    i1 = this%Site1%NTest
    j1 = this%Site2%NPart
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared
    nu1 = this%NUnit1
    nu2 = this%NUnit2
#if MPI_VER > 0
    i0 = this%Site1%NTest0
    i1 = this%Site1%NTest2
#else
    i0 = 1
    i1 = this%Site1%NTest
#endif

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
!$OMP PRIVATE (EPotLocal,i,i0,i1,j,k) 
    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
!$OMP DO
      do i = i0, i1
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
#if ARCH == 3
        hit = .false.
#endif
!CDIR NODEP

        unit = nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then  ! choose only units, to which our Site2 correspond
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)   ! number of molecule, to which this unit correspond
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
            if( RijSquared <= RShieldSquared ) hit = .true.
#else
            if( RijSquared <= RShieldSquared ) then
              EPotLocal = 1E33_RK
              exit loop1
            end if
#endif
            OXj = OX2(jk)
            OYj = OY2(jk)
            OZj = OZ2(jk)

#if ARCH == 3
            RijInv = rsqrt( RijSquared )
#else
            RijInv = 1._RK / sqrt( RijSquared )
#endif
            eX = RXij * RijInv
            eY = RYij * RijInv
            eZ = RZij * RijInv
            CosThetai = OXi * eX + OYi * eY + OZi * eZ
            CosThetaj = OXj * eX + OYj * eY + OZj * eZ
            CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
            Rij4Inv = Epsilon / RijSquared**2
            EPotLocal = EPotLocal + Rij4Inv * (2._RK * CosGammaij * CosThetaj - CosThetai * (5._RK * CosThetaj**2 - 1._RK))
          end if
        end do loop1

#if ARCH == 3
        if( .not. hit ) then
          EPotTest(i) = EPotTest(i) + EPotLocal
        else
          EPotTest(i) = EPotTest(i) + 1E33_RK
        endif
#else
        EPotTest(i) = EPotTest(i) + EPotLocal
#endif
      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over test particles
!$OMP DO
      do i = i0, i1
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        EPotLocal = 0._RK

#if ARCH == 3
        hit = .false.
#endif
!CDIR NODEP
loop2:  do j = 1, j1
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
          if( RijSquared <= RShieldSquared ) hit = .true.
#else
          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop2
          end if
#endif
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal = EPotLocal + Rij4Inv * (2._RK * CosGammaij * CosThetaj - CosThetai * (5._RK * CosThetaj**2 - 1._RK))
        end do loop2

#if ARCH == 3
        if( .not. hit ) then
          EPotTest(i) = EPotTest(i) + EPotLocal
        else
          EPotTest(i) = EPotTest(i) + 1E33_RK
        endif
#else
        EPotTest(i) = EPotTest(i) + EPotLocal
#endif
      end do
!$OMP END DO
    end if
!$OMP END PARALLEL
  end subroutine TPotDQ_ChemicalPotential



!==============================================================!
!  Subroutine TPotDQ_Energy                                    !
!==============================================================!

  subroutine TPotDQ_Energy( this, np, nu, F, E, EIntra, BoxLength, CompIdent )

    implicit none

    ! Declare arguments
    type(TPotDipoleQuadrupole) :: this
    integer, intent(in)        :: np
    integer, intent(in)        :: nu
    real(RK), intent(in out)   :: F(3,nu)
    real(RK), intent(in out)   :: E
    real(RK), intent(in out)   :: EIntra
    real(RK), intent(in)       :: BoxLength
    logical, intent(in)        :: CompIdent

    ! Declare local variables
    real(RK), pointer, contiguous :: RX2(:), RY2(:), RZ2(:)
    real(RK), pointer, contiguous :: OX2(:), OY2(:), OZ2(:)
    real(RK), pointer, contiguous :: PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij4Inv
    real(RK)          :: CosThetai, CosThetaj, CosThetaj2, CosGammaij
    real(RK)          :: E1, EIntra1, ELocal, tempF(3,nu)
    integer           :: j, k, nu2, jk, unit
    real(RK)          :: coeff

    ! Assign local variables
    nu2 = this%NUnit2
    coeff = 1._RK
    if (this%potintra14) coeff = this%ScaleEl14
    E1   = 0._RK
    EIntra1   = 0._RK

    ! Assign pointers
    RX2 => this%Site2%RX
    RY2 => this%Site2%RY
    RZ2 => this%Site2%RZ
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

    ! Loop over molecules
    RXi = this%Site1%RX(np)
    RYi = this%Site1%RY(np)
    RZi = this%Site1%RZ(np)
    OXi = this%Site1%OX(np)
    OYi = this%Site1%OY(np)
    OZi = this%Site1%OZ(np)
    PXi = this%Site1%PX(np)
    PYi = this%Site1%PY(np)
    PZi = this%Site1%PZ(np)

    unit=this%NUnit1*(np-1)+this%Site1%UnitNumber

    do k = 1, this%NInCutoff(unit)
      j = this%CutoffPartner(k, unit)
      if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
        if (mod(j,nu2)==0) then
          jk = INT(j/nu2)
        else
          jk = INT(j/nu2)+1
        end if
        RXij = RXi - RX2(jk)
        RYij = RYi - RY2(jk)
        RZij = RZi - RZ2(jk)
        PXij = PXi - PX2(jk)
        PYij = PYi - PY2(jk)
        PZij = PZi - PZ2(jk)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        RijSquared = RXij**2 + RYij**2 + RZij**2

        if( RijSquared <= this%RShieldSquared ) then
          ELocal = 1E33_RK
        else
          OXj = OX2(jk)
          OYj = OY2(jk)
          OZj = OZ2(jk)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosThetaj2 = CosThetaj**2
          CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
          Rij4Inv = this%Epsilon / RijSquared**2
          ELocal = Rij4Inv * (CosGammaij * CosThetaj - CosThetai * (5._RK * CosThetaj2 - 1))
        end if
        E1 = E1 + ELocal
      end if
    end do
    ! Include intramolecular interaction if need
    if (this%potintra15 .or. this%potintra14) then
      RXij = RXi - RX2(np)
      RYij = RYi - RY2(np)
      RZij = RZi - RZ2(np)
      PXij = PXi - PX2(np)
      PYij = PYi - PY2(np)
      PZij = PZi - PZ2(np)
      RXij = (RXij - anint( PXij )) * BoxLength
      RYij = (RYij - anint( PYij )) * BoxLength
      RZij = (RZij - anint( PZij )) * BoxLength
      OXj = OX2(np)
      OYj = OY2(np)
      OZj = OZ2(np)
      RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
      RijInv = rsqrt( RijSquared )
#else
      RijInv = 1._RK / sqrt( RijSquared )
#endif
      eX = RXij * RijInv
      eY = RYij * RijInv
      eZ = RZij * RijInv
      CosThetai = OXi * eX + OYi * eY + OZi * eZ
      CosThetaj = OXj * eX + OYj * eY + OZj * eZ
      CosThetaj2 = CosThetaj**2
      CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
      Rij4Inv = this%Epsilon / RijSquared**2
      ELocal = coeff * Rij4Inv * (CosGammaij * CosThetaj - CosThetai * (5._RK * CosThetaj2 - 1))
      EIntra1 = EIntra1 + ELocal
    end if

    F(:,:) = F(:,:) + tempF(:,:)
    E = E + E1 + EIntra1
    EIntra = EIntra + EIntra1

  end subroutine TPotDQ_Energy


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

    ! Declare local variables
    integer :: k, ende

    ! Construct potential
    this%Site1 => Molecule1%SiteQuadrupole(j1)
    this%Site2 => Molecule2%SiteCharge(j2)
    this%NUnit1 => Molecule1%NUnit
    this%NUnit2 => Molecule2%NUnit
    this%SameComponent = i1 == i2
    this%Epsilon = 1.5_RK * this%Site1%Q * this%Site2%e
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

    ! if this potential is intra
    if (this%SameComponent .and. Molecule1%hasIntraLJEl) then
      ende = size(Molecule1%IntQC15(:,1))
      do k=1, ende
        if (Molecule1%IntQC15(k,1)==this%Site1%SiteId .and. Molecule1%IntQC15(k,2)==this%Site2%SiteId) then
          this%potintra15 = .true.
        end if
      end do
      if (LJEL14 .and. .not. this%potintra15) then
        ende = size(Molecule1%IntQC14(:,1))
        do k=1, ende
          if (Molecule1%IntQC14(k,1)==this%Site1%SiteId .and. Molecule1%IntQC14(k,2)==this%Site2%SiteId) then
            this%potintra14=.true.
            this%ScaleEl14 = Molecule1%ScaleQC14(k)
          end if
        end do
      end if
    else
      this%potintra15 = .false.
      this%potintra14 = .false.
    end if

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

  subroutine TPotQC_Force( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleCharge) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: EPotInter
    real(RK), intent(in out)   :: VirialInter
    real(RK), intent(in out)   :: EPotIntra_Nonbonded
    real(RK), intent(in out)   :: VirialIntra
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
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, jk, unit
    logical           :: SameComponent
    real(RK)          :: coeff

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
    SameComponent = this%SameComponent
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

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

      unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:do k = 1, this%NInCutoff(unit)
        j = this%CutoffPartner(k, unit)
        if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
          if (mod(j,nu2)==0) then
            jk = INT(j/nu2)
          else
            jk = INT(j/nu2)+1
          end if
          RXij = RXi - RX2(jk)
          RYij = RYi - RY2(jk)
          RZij = RZi - RZ2(jk)
          PXij = PXi - PX2(jk)
          PYij = PYi - PY2(jk)
          PZij = PZi - PZ2(jk)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
          RijInv = sqrt( RijSquaredInv )
          eX = - RXij * RijInv                          ! Normierter Abstandsvektor nach Price
          eY = - RYij * RijInv
          eZ = - RZij * RijInv
          CosTheta  = OXi * ex + OYi * eY + OZi * eZ          ! Scalarprodukt normierter 
!                                              Abstandsvektor mit Orientierungsvektor Quadrupol
          Epsilon1 = Epsilon * RijSquaredInv * RijInv
          EPotLocal  = EPotLocal + Epsilon1 * ( CosTheta * CosTheta - Third )
          EPotLocalInter  = EPotLocalInter + Epsilon1 * ( CosTheta * CosTheta - Third )
          CosTheta2 = 2._RK * CosTheta
          CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
          Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
          FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXi ) ! Kraft auf die Punktladung, sprich F2
          FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYi )
          FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZi )
          VirialLocal = VirialLocal - (FXij * PXij + FYij * PYij + FZij * PZij)     ! Vorzeichen richtig
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
          VirialLocalInter = VirialLocalInter - (FXij * PXij - FYij * PYij - FZij * PZij)
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
          d2EpotdV2Local = d2EpotdV2Local + Epsilon1 * ( CosTheta * CosTheta - Third )*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijSquaredInv)*Third*Third    !xxxx7  QC

          FXi    = FXi    - FXij
          FYi    = FYi    - FYij
          FZi    = FZi    - FZij

          forceTempX(jk) = forceTempX(jk) + FXij
          forceTempY(jk) = forceTempY(jk) + FYij
          forceTempZ(jk) = forceTempZ(jk) + FZij

          TXi    = TXi - Epsilon1*CosTheta2*eX  
          ! Drehmomentanteil auf Quadrupol wegen Punktladung. Kreuzprodukt
          TYi    = TYi - Epsilon1*CosTheta2*eY  ! in Atom2Mol von Component
          TZi    = TZi - Epsilon1*CosTheta2*eZ
        end if
      end do loop1

      ! Include intramolecular interaction if need
      if (this%potintra15 .or. this%potintra14) then
        RXij = RXi - RX2(i)
        RYij = RYi - RY2(i)
        RZij = RZi - RZ2(i)
        PXij = PXi - PX2(i)
        PYij = PYi - PY2(i)
        PZij = PZi - PZ2(i)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
        RijInv = sqrt( RijSquaredInv )
        eX = - RXij * RijInv                  ! Normierter Abstandsvektor nach Price
        eY = - RYij * RijInv
        eZ = - RZij * RijInv
        CosTheta  = OXi * ex + OYi * eY + OZi * eZ  ! Scalarprodukt normierter Abstandsvektor mit Orientierungsvektor Quadrupol
        Epsilon1 = Epsilon * RijSquaredInv * RijInv * coeff
        EPotLocal  = EPotLocal + Epsilon1 * ( CosTheta * CosTheta - Third )
        EPotLocalIntra  = EPotLocalIntra + Epsilon1 * ( CosTheta * CosTheta - Third )
        CosTheta2 = 2._RK * CosTheta
        CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
        Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv * coeff
        FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXi )   ! Kraft auf die Punktladung, sprich F2
        FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYi )
        FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZi )
        VirialLocal = VirialLocal - (FXij * PXij + FYij * PYij + FZij * PZij)     ! Vorzeichen richtig
        VirialLocalIntra = VirialLocalIntra - (FXij * PXij - FYij * PYij - FZij * PZij)
        Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
        sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
        d2EpotdV2Local = d2EpotdV2Local + Epsilon1 * ( CosTheta * CosTheta - Third )*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijSquaredInv)*Third*Third    !xxxx7  QC
        FXi    = FXi    - FXij
        FYi    = FYi    - FYij
        FZi    = FZi    - FZij

        forceTempX(i) = forceTempX(i) + FXij
        forceTempY(i) = forceTempY(i) + FYij
        forceTempZ(i) = forceTempZ(i) + FZij

        TXi    = TXi - Epsilon1*CosTheta2*eX  
        ! Drehmomentanteil auf Quadrupol wegen Punktladung. Kreuzprodukt
        TYi    = TYi - Epsilon1*CosTheta2*eY  ! in Atom2Mol von Component
        TZi    = TZi - Epsilon1*CosTheta2*eZ
      end if

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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotQC_Force

!==============================================================!
!  Subroutine TPotQC_Force_Trans                               !
!==============================================================!

  subroutine TPotQC_Force_Trans( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleCharge) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: EPotInter
    real(RK), intent(in out)   :: VirialInter
    real(RK), intent(in out)   :: EPotIntra_Nonbonded
    real(RK), intent(in out)   :: VirialIntra
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
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, jk, unit
    logical           :: SameComponent
    real(RK)          :: coeff

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
    real(RK), pointer, contiguous :: VSx(:), VSy(:), VSz(:)
    real(RK), pointer, contiguous :: VSux(:),VSuy(:),VSuz(:)
    real(RK), pointer, contiguous :: VBx(:), VBy(:), VBz(:)
    real(RK), pointer, contiguous :: Cx(:) , Cy(:) , Cz(:)
    real(RK), pointer, contiguous :: tux(:) , tuy(:) , tuz(:)
    real(RK), pointer, contiguous :: tlx(:) , tly(:) , tlz(:)
    real(RK), pointer, contiguous :: tdx(:) , tdy(:) , tdz(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
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
    real(RK)          :: UU, Uxi,  Uyi, Uzi
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33
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
!$OMP PRIVATE ( Epsilon, Epsilon1, Epsilon2, RX1, RY1, RZ1, RX2, RY2, RZ2) &
!$OMP PRIVATE (  FX1, FY1, FZ1, OX1, OY1, OZ1, TX1, TY1, TZ1) &
!$OMP PRIVATE (Plen2,sitecorr, PX1, PY1, PZ1, PX2, PY2, PZ2, TXi, TYi, TZi) &
!$OMP PRIVATE (   RXi, RYi, RZi, FXi, FYi, FZi, PXi, PYi, PZi)&
!$OMP PRIVATE (   RXij, RYij, RZij, FXij, FYij, FZij, PXij, PYij, PZij) &
!$OMP PRIVATE (  eX, eY, eZ, RijSquaredInv, RijInv) &
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
    SameComponent = this%SameComponent
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

#if  TRANS == 1
    !TRANSPORT_start
    VSx => this%Site1%vsQx
    VSy => this%Site1%vsQy
    VSz => this%Site1%vsQz
    VBx => this%Site1%vbQx
    VBy => this%Site1%vbQy
    VBz => this%Site1%vbQz
    VSux=> this%Site1%vsuQx
    VSuy=> this%Site1%vsuQy
    VSuz=> this%Site1%vsuQz
    Cx  => this%Site1%cQx
    Cy  => this%Site1%cQy
    Cz  => this%Site1%cQz
    tux => this%Site1%tuQx
    tuy => this%Site1%tuQy
    tuz => this%Site1%tuQz
    tlx => this%Site1%tlQx
    tly => this%Site1%tlQy
    tlz => this%Site1%tlQz
    tdx => this%Site1%tdQx
    tdy => this%Site1%tdQy
    tdz => this%Site1%tdQz
    q1  => this%Site1%Qm0r(:, 1, 1)
    q2  => this%Site1%Qm0r(:, 2, 1)
    q3  => this%Site1%Qm0r(:, 3, 1)
    q4  => this%Site1%Qm0r(:, 4, 1)
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
        VSxi= VSx(i)
        VSyi= VSy(i)
        VSzi= VSz(i)
        VBxi= VBx(i)
        VByi= VBy(i)
        VBzi= VBz(i)
        VSuxi= VSux(i)
        VSuyi= VSuy(i)
        VSuzi= VSuz(i)
        Cxi = Cx(i)
        Cyi = Cy(i)
        Czi = Cz(i)
        tuxi = tux(i)
        tuyi = tuy(i)
        tuzi = tuz(i)
        tlxi = tlx(i)
        tlyi = tly(i)
        tlzi = tlz(i)
        tdxi = tdx(i)
        tdyi = tdy(i)
        tdzi = tdz(i)
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
loop0:do m=1,NBinsDen
        if (PXi .ge. real(m-1)/NBinsDen-0.5_RK) then
          if (PXi < real(m)/NBinsDen-0.5_RK) then
            Bin1=m
            exit loop0
          end if
        end if
      end do loop0
#endif

      unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:do k = 1, this%NInCutoff(unit)
        j = this%CutoffPartner(k, unit)
        if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
          if (mod(j,nu2)==0) then
            jk = INT(j/nu2)
          else
            jk = INT(j/nu2)+1
          end if
          RXij = RXi - RX2(jk)
          RYij = RYi - RY2(jk)
          RZij = RZi - RZ2(jk)
          PXij = PXi - PX2(jk)
          PYij = PYi - PY2(jk)
          PZij = PZi - PZ2(jk)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
          RijInv = sqrt( RijSquaredInv )
          eX = - RXij * RijInv               ! Normierter Abstandsvektor nach Price
          eY = - RYij * RijInv
          eZ = - RZij * RijInv
          CosTheta  = OXi * ex + OYi * eY + OZi * eZ    ! Scalarprodukt normierter 
!                                              Abstandsvektor mit Orientierungsvektor Quadrupol
          Epsilon1 = Epsilon * RijSquaredInv * RijInv
          EPotLocal1 = Epsilon1 * ( CosTheta * CosTheta - Third )
          EPotLocal  = EPotLocal + EPotLocal1
          EPotLocalInter  = EPotLocalInter + EPotLocal1
          CosTheta2 = 2._RK * CosTheta
          CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
          Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
          FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXi ) ! Kraft auf die Punktladung, sprich F2
          FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYi )
          FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZi )
          VirialLocal = VirialLocal - (FXij * PXij + FYij * PYij + FZij * PZij)  ! Vorzeichen richtig
#if OSMOP == 2
loop2:  do m=1,NBinsDen
          if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
            if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
          VirialLocalInter = VirialLocalInter - (FXij * PXij - FYij * PYij - FZij * PZij)
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
          d2EpotdV2Local = d2EpotdV2Local + Epsilon1 * ( CosTheta * CosTheta - Third )*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijSquaredInv)*Third*Third    !xxxx7  QC T
          FXi    = FXi    - FXij
          FYi    = FYi    - FYij
          FZi    = FZi    - FZij

          forceTempX(jk) = forceTempX(jk) + FXij
          forceTempY(jk) = forceTempY(jk) + FYij
          forceTempZ(jk) = forceTempZ(jk) + FZij

          TXi    = TXi - Epsilon1*CosTheta2*eX  
          ! Drehmomentanteil auf Quadrupol wegen Punktladung. Kreuzprodukt
          TYi    = TYi - Epsilon1*CosTheta2*eY  ! in Atom2Mol von Component
          TZi    = TZi - Epsilon1*CosTheta2*eZ
#if  TRANS == 1
!TRANSPORT_start
          VSxi   = VSxi + FXij * PYij
          VSyi   = VSyi + FXij * PZij
          VSzi   = VSzi + FYij * PZij
          VBxi   = VBxi + FXij * PXij
          VByi   = VByi + FYij * PYij
          VBzi   = VBzi + FZij * PZij
          VSuxi  = VSuxi+ FYij * PXij
          VSuyi  = VSuyi+ FZij * PXij
          VSuzi  = VSuzi+ FZij * PYij
          UU     = EpotLocal1  
          Uxi    = UU * eX
          Uyi    = UU * eY
          Uzi    = UU * eZ
          Cxi    = Cxi  + Uxi
          Cyi    = Cyi  + Uyi
          Czi    = Czi  + Uzi
          FTXi   = - Epsilon1*CosTheta2*eX 
          FTYi   = - Epsilon1*CosTheta2*eY 
          FTZi   = - Epsilon1*CosTheta2*eZ 
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
          !TRANSPORT_END
#endif
        end if
      end do loop1

      ! Include intramolecular interaction if need
      if (this%potintra15 .or. this%potintra14) then
        RXij = RXi - RX2(i)
        RYij = RYi - RY2(i)
        RZij = RZi - RZ2(i)
        PXij = PXi - PX2(i)
        PYij = PYi - PY2(i)
        PZij = PZi - PZ2(i)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength
        RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
        RijInv = sqrt( RijSquaredInv )
        eX = - RXij * RijInv            ! Normierter Abstandsvektor nach Price
        eY = - RYij * RijInv
        eZ = - RZij * RijInv
        CosTheta  = OXi * ex + OYi * eY + OZi * eZ
        ! Scalarprodukt normierter Abstandsvektor mit Orientierungsvektor Quadrupol
        Epsilon1 = Epsilon * RijSquaredInv * RijInv * coeff
        EPotLocal1 = Epsilon1 * ( CosTheta * CosTheta - Third )
        EPotLocal  = EPotLocal + EPotLocal1
        EPotLocalIntra  = EPotLocalIntra + EPotLocal1
        CosTheta2 = 2._RK * CosTheta
        CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
        Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv * coeff
        FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXi )          ! Kraft auf die Punktladung, sprich F2
        FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYi )
        FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZi )
        VirialLocal = VirialLocal - (FXij * PXij + FYij * PYij + FZij * PZij)     ! Vorzeichen richtig
        VirialLocalIntra = VirialLocalIntra - (FXij * PXij - FYij * PYij - FZij * PZij)
        FXi    = FXi    - FXij
        FYi    = FYi    - FYij
        FZi    = FZi    - FZij
        forceTempX(i) = forceTempX(i) + FXij
        forceTempY(i) = forceTempY(i) + FYij
        forceTempZ(i) = forceTempZ(i) + FZij
        TXi    = TXi - Epsilon1*CosTheta2*eX    ! Drehmomentanteil auf Quadrupol wegen Punktladung. Kreuzprodukt
        TYi    = TYi - Epsilon1*CosTheta2*eY    ! in Atom2Mol von Component
        TZi    = TZi - Epsilon1*CosTheta2*eZ

      end if

      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
      TX1(i) = TXi
      TY1(i) = TYi
      TZ1(i) = TZi
#if  TRANS == 1
        !TRANSPORT_start
        VSx(i) = VSxi
        VSy(i) = VSyi
        VSz(i) = VSzi
        VBx(i) = VBxi
        VBy(i) = VByi
        VBz(i) = VBzi
        VSux(i)= VSuxi
        VSuy(i)= VSuyi
        VSuz(i)= VSuzi
        Cx(i)  = Cxi
        Cy(i)  = Cyi
        Cz(i)  = Czi
        tux(i) = tuxi
        tuy(i) = tuyi
        tuz(i) = tuzi
        tlx(i) = tlxi
        tly(i) = tlyi
        tlz(i) = tlzi
        tdx(i) = tdxi
        tdy(i) = tdyi
        tdz(i) = tdzi
        !TRANSPORT_END
#endif
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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
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
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
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
    integer           :: nu1, nu2, unit, jk, i0
#if ARCH == 3
    logical           :: hit
#endif

    ! Assign local variables
    i1 = this%Site1%NTest
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared
    nu1 = this%NUnit1
    nu2 = this%NUnit2
#if MPI_VER > 0
    i0 = this%Site1%NTest0
    i1 = this%Site1%NTest2
#else
    i0 = 1
    i1 = this%Site1%NTest
#endif

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
    do i = i0, i1
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
#if ARCH == 3
        hit = .false.
#endif

!CDIR NODEP

      unit = nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then  ! choose only units, to which our Site2 correspond
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)   ! number of molecule, to which this unit correspond
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            PXij = (PXij - anint( PXij )) * BoxLength
            PYij = (PYij - anint( PYij )) * BoxLength
            PZij = (PZij - anint( PZij )) * BoxLength
            RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
            if( RijSquared <= RShieldSquared ) hit = .true.
#else
            if( RijSquared <= RShieldSquared ) then
              EPotLocal = 1E33_RK
              exit loop1
            end if
#endif
            RijSquaredInv = 1._RK / RijSquared
            RijInv = sqrt( RijSquaredInv )
            eX = RXij * RijInv
            eY = RYij * RijInv
            eZ = RZij * RijInv
            CosTheta  = OXi * ex + OYi * eY + OZi * eZ
            EPotLocal = EPotLocal + Epsilon * RijSquaredInv * RijInv * ( CosTheta * CosTheta - Third )
          end if
        end do loop1

#if ARCH == 3
        if( .not. hit ) then
          EPotTest(i) = EPotTest(i) + EPotLocal
        else
          EPotTest(i) = EPotTest(i) + 1E33_RK
        endif
#else
        EPotTest(i) = EPotTest(i) + EPotLocal
#endif
    end do

  end subroutine TPotQC_ChemicalPotential



!==============================================================!
!  Subroutine TPotQC_Energy                                    !
!==============================================================!

  subroutine TPotQC_Energy( this, np, nu, F, E, EIntra, BoxLength, CompIdent )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleCharge) :: this
    integer, intent(in)        :: np
    integer, intent(in)        :: nu
    real(RK), intent(in out)   :: F(3,nu)
    real(RK), intent(in out)   :: E
    real(RK), intent(in out)   :: EIntra
    real(RK), intent(in)       :: BoxLength
    logical, intent(in)        :: CompIdent

    ! Declare local variables
    real(RK)          :: Epsilon
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
    real(RK)          :: CosTheta, CosTheta2, CosAux
    real(RK)          :: EPot, EIntra1, EPotLocal, tempF(3,nu)
    integer           :: j, k, nu2, jk, unit
    real(RK)          :: coeff

    ! Assign local variables
    Epsilon = this%Epsilon
    RShieldSquared = this%RShieldSquared
    nu2 = this%NUnit2
    coeff = 1._RK
    if (this%potintra14) coeff = this%ScaleEl14
    EPot   = 0._RK
    EIntra1   = 0._RK

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
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ

    ! Loop over molecules
    RXi = RX1(np)
    RYi = RY1(np)
    RZi = RZ1(np)
    PXi = PX1(np)
    PYi = PY1(np)
    PZi = PZ1(np)
    OXi = OX1(np)
    OYi = OY1(np)
    OZi = OZ1(np)

    unit=this%NUnit1*(np-1)+this%Site1%UnitNumber

    do k = 1, this%NInCutoff(unit)
      j = this%CutoffPartner(k, unit)
      if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
        if (mod(j,nu2)==0) then
          jk = INT(j/nu2)
        else
          jk = INT(j/nu2)+1
        end if
        RXij = RXi - RX2(jk)
        RYij = RYi - RY2(jk)
        RZij = RZi - RZ2(jk)
        PXij = PXi - PX2(jk)
        PYij = PYi - PY2(jk)
        PZij = PZi - PZ2(jk)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        RijSquared = RXij**2 + RYij**2 + RZij**2

        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
        else
          RijSquaredInv = 1._RK / RijSquared
          RijInv = sqrt( RijSquaredInv )
          eX = - RXij * RijInv        ! Normierter Abstandsvektor nach Price
          eY = - RYij * RijInv
          eZ = - RZij * RijInv
          CosTheta  = OXi * ex + OYi * eY + OZi * eZ  
          ! Scalarprodukt normierter Abstandsvektor mit Orientierungsvektor Quadrupol
          EPotLocal = Epsilon * RijSquaredInv * RijInv * ( CosTheta * CosTheta - Third )
        end if

        EPot = EPot + EPotLocal
      end if
    end do
    ! Include intramolecular interaction if need
    if (this%potintra15 .or. this%potintra14) then
      RXij = RXi - RX2(np)
      RYij = RYi - RY2(np)
      RZij = RZi - RZ2(np)
      PXij = PXi - PX2(np)
      PYij = PYi - PY2(np)
      PZij = PZi - PZ2(np)
      RXij = (RXij - anint( PXij )) * BoxLength
      RYij = (RYij - anint( PYij )) * BoxLength
      RZij = (RZij - anint( PZij )) * BoxLength
      RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
      RijInv = sqrt( RijSquaredInv )
      eX = - RXij * RijInv
      eY = - RYij * RijInv
      eZ = - RZij * RijInv
      CosTheta  = OXi * ex + OYi * eY + OZi * eZ
      EPotLocal = Epsilon * RijSquaredInv * RijInv * coeff * ( CosTheta * CosTheta - Third )
      EIntra1  = EIntra1 + EPotLocal
    end if

    F(:,:) = F(:,:) + tempF(:,:)
    E = E + EPot + EIntra1
    EIntra = EIntra + EIntra1

  end subroutine TPotQC_Energy



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

    ! Declare local variables
    integer :: k, ende

    ! Construct potential
    this%Site1 => Molecule1%SiteQuadrupole(j1)
    this%Site2 => Molecule2%SiteDipole(j2)
    this%NUnit1 => Molecule1%NUnit
    this%NUnit2 => Molecule2%NUnit
    this%SameComponent = i1 == i2
    this%Epsilon = 1.5_RK * this%Site1%Q * this%Site2%D
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

    ! if this potential is intra
   if (this%SameComponent .and. Molecule1%hasIntraLJEl ) then
      ende = size(Molecule1%IntQD15(:,1))
      do k=1, ende
        if (Molecule1%IntQD15(k,1)==this%Site1%SiteId .and. Molecule1%IntQD15(k,2)==this%Site2%SiteId) then
          this%potintra15 = .true.
        end if
      end do
      if (LJEL14 .and. .not. this%potintra15) then
        ende = size(Molecule1%IntQD14(:,1))
        do k=1, ende
          if (Molecule1%IntQD14(k,1)==this%Site1%SiteId .and. Molecule1%IntQD14(k,2)==this%Site2%SiteId) then
            this%potintra14=.true.
            this%ScaleEl14 = Molecule1%ScaleQD14(k)
          end if
        end do
      end if
    else
      this%potintra15 = .false.
      this%potintra14 = .false.
    end if

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

  subroutine TPotQD_Force( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleDipole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: EPotInter
    real(RK), intent(in out)   :: VirialInter
    real(RK), intent(in out)   :: EPotIntra_Nonbonded
    real(RK), intent(in out)   :: VirialIntra
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
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
    integer           :: nu1, nu2, jk, unit
    real(RK)          :: coeff
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
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter   = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra   = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

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

        unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            PXij = (PXij - anint( PXij )) * BoxLength
            PYij = (PYij - anint( PYij )) * BoxLength
            PZij = (PZij - anint( PZij )) * BoxLength
            RijSquared = RXij**2 + RYij**2 + RZij**2
            OXj = OX2(jk)
            OYj = OY2(jk)
            OZj = OZ2(jk)
#if ARCH == 3
            RijInv = rsqrt( RijSquared )
#else
            RijInv = 1._RK / sqrt( RijSquared )
#endif
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
            EPotLocalInter = EPotLocalInter + EPotLocal1
            dCosThetai = Rij4Inv * (10._RK * CosThetai * CosThetaj - CosGammaij)
            dCosThetaj = Rij4Inv * (5._RK * CosThetai2 - 1._RK)
            dCosGammaij = -2._RK * Rij4Inv * CosThetai
            Tmp = -4._RK * RijInv * EPotLocal1

            FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                      + (eX * CosThetaj - OXj) * dCosThetaj)
            FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                      + (eY * CosThetaj - OYj) * dCosThetaj)
            FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                      + (eZ * CosThetaj - OZj) * dCosThetaj)

            VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
            VirialLocalInter = VirialLocalInter + (FXij * PXij + FYij * PYij + FZij * PZij)
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
            d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Third*Third    !xxxx8   QD

            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij

            forceTempX(jk) = forceTempX(jk) - FXij
            forceTempY(jk) = forceTempY(jk) - FYij
            forceTempZ(jk) = forceTempZ(jk) - FZij
            TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
            TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
            TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij

            momTempX(jk) = momTempX(jk) - eX * dCosThetaj - OXi * dCosGammaij
            momTempY(jk) = momTempY(jk) - eY * dCosThetaj - OYi * dCosGammaij  
            momTempZ(jk) = momTempZ(jk) - eZ * dCosThetaj - OZi * dCosGammaij   
          end if

        end do loop1
        ! Include intramolecular interaction if need
        if (SameComponent .and. (this%potintra14 .or. this%potintra15)) then
          RXij = RXi - RX2(i)
          RYij = RYi - RY2(i)
          RZij = RZi - RZ2(i)
          PXij = PXi - PX2(i)
          PYij = PYi - PY2(i)
          PZij = PZi - PZ2(i)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          RijSquared = RXij**2 + RYij**2 + RZij**2
          OXj = OX2(i)
          OYj = OY2(i)
          OZj = OZ2(i)
#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosThetai2 = CosThetai**2
          CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal1 = Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) &
&                                  - CosGammaij * CosThetai)
          EPotLocal = EPotLocal + EPotLocal1*coeff
          EPotLocalIntra = EPotLocalIntra + EPotLocal1*coeff
          dCosThetai = Rij4Inv * (10._RK * CosThetai * CosThetaj - CosGammaij)
          dCosThetaj = Rij4Inv * (5._RK * CosThetai2 - 1._RK)
          dCosGammaij = -2._RK * Rij4Inv * CosThetai
          Tmp = -4._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FXij = FXij * coeff
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FYij = FYij * coeff
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)
          FZij = FZij * coeff

          VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
          VirialLocalIntra = VirialLocalIntra + FXij * PXij + FYij * PYij + FZij * PZij
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + coeff * EPotLocal1*(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Third*Third    !xxxx8   QD
          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(i) = forceTempX(i) - FXij
          forceTempY(i) = forceTempY(i) - FYij
          forceTempZ(i) = forceTempZ(i) - FZij
          TXi    = TXi    - (eX * dCosThetai + OXj * dCosGammaij) * coeff
          TYi    = TYi    - (eY * dCosThetai + OYj * dCosGammaij) * coeff
          TZi    = TZi    - (eZ * dCosThetai + OZj * dCosGammaij) * coeff

          momTempX(i) = momTempX(i) - (eX * dCosThetaj + OXi * dCosGammaij) * coeff
          momTempY(i) = momTempY(i) - (eY * dCosThetaj + OYi * dCosGammaij) * coeff
          momTempZ(i) = momTempZ(i) - (eZ * dCosThetaj + OZi * dCosGammaij) * coeff
        end if

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
          RijSquared = RXij**2 + RYij**2 + RZij**2
          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
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
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Third*Third    !xxxx8   QD ss

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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotQD_Force

!==============================================================!
!  Subroutine TPotQD_Force_Trans                               !
!==============================================================!

  subroutine TPotQD_Force_Trans( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleDipole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out)   :: EPotInter
    real(RK), intent(in out)   :: VirialInter
    real(RK), intent(in out)   :: EPotIntra_Nonbonded
    real(RK), intent(in out)   :: VirialIntra
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
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
    integer           :: nu1, nu2, jk, unit
    real(RK)          :: coeff

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
    real(RK), pointer, contiguous :: VSx(:), VSy(:), VSz(:)
    real(RK), pointer, contiguous :: VSux(:),VSuy(:),VSuz(:)
    real(RK), pointer, contiguous :: VBx(:), VBy(:), VBz(:)
    real(RK), pointer, contiguous :: Cx(:) , Cy(:) , Cz(:)
    real(RK), pointer, contiguous :: tux(:) , tuy(:) , tuz(:)
    real(RK), pointer, contiguous :: tlx(:) , tly(:) , tlz(:)
    real(RK), pointer, contiguous :: tdx(:) , tdy(:) , tdz(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
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
    real(RK)          :: UU, Uxi,  Uyi, Uzi
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33
    !TRANSPORT_END
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
!$OMP PRIVATE(VSx, VSy, VSz ,VSux,VSuy,VSuz, VBx, VBy, VBz, Cx , Cy , Cz) &
!$OMP PRIVATE( tux , tuy , tuz, tlx , tly , tlz, tdx , tdy , tdz) &
!$OMP PRIVATE( q1, q2, q3, q4, VSxi, VSyi, VSzi, VSuxi,VSuyi,VSuzi) &
!$OMP PRIVATE( VBxi, VByi, VBzi, Cxi,  Cyi,  Czi, tuxi,  tuyi,  tuzi, tlxi,  tlyi,  tlzi) &
!$OMP PRIVATE(  tdxi,  tdyi,  tdzi, txii,  tyii , tzii, txir ,  tyir  , tzir ) &
!$OMP PRIVATE(   Uxi,  Uyi, Uzi, FTXi , FTYi , FTZi) &
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
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter   = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra   = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

#if  TRANS == 1
    !TRANSPORT_start
    VSx => this%Site1%vsQx
    VSy => this%Site1%vsQy
    VSz => this%Site1%vsQz
    VBx => this%Site1%vbQx
    VBy => this%Site1%vbQy
    VBz => this%Site1%vbQz
    VSux=> this%Site1%vsuQx
    VSuy=> this%Site1%vsuQy
    VSuz=> this%Site1%vsuQz
    Cx  => this%Site1%cQx
    Cy  => this%Site1%cQy
    Cz  => this%Site1%cQz
    tux => this%Site1%tuQx
    tuy => this%Site1%tuQy
    tuz => this%Site1%tuQz
    tlx => this%Site1%tlQx
    tly => this%Site1%tlQy
    tlz => this%Site1%tlQz
    tdx => this%Site1%tdQx
    tdy => this%Site1%tdQy
    tdz => this%Site1%tdQz
    q1  => this%Site1%Qm0r(:, 1, 1)
    q2  => this%Site1%Qm0r(:, 2, 1)
    q3  => this%Site1%Qm0r(:, 3, 1)
    q4  => this%Site1%Qm0r(:, 4, 1)
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
        VSxi= VSx(i)
        VSyi= VSy(i)
        VSzi= VSz(i)
        VBxi= VBx(i)
        VByi= VBy(i)
        VBzi= VBz(i)
        VSuxi= VSux(i)
        VSuyi= VSuy(i)
        VSuzi= VSuz(i)
        Cxi = Cx(i)
        Cyi = Cy(i)
        Czi = Cz(i)
        tuxi = tux(i)
        tuyi = tuy(i)
        tuzi = tuz(i)
        tlxi = tlx(i)
        tlyi = tly(i)
        tlzi = tlz(i)
        tdxi = tdx(i)
        tdyi = tdy(i)
        tdzi = tdz(i)
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

        unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            PXij = (PXij - anint( PXij )) * BoxLength
            PYij = (PYij - anint( PYij )) * BoxLength
            PZij = (PZij - anint( PZij )) * BoxLength
            RijSquared = RXij**2 + RYij**2 + RZij**2
            OXj = OX2(jk)
            OYj = OY2(jk)
            OZj = OZ2(jk)

#if ARCH == 3
            RijInv = rsqrt( RijSquared )
#else
            RijInv = 1._RK / sqrt( RijSquared )
#endif
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
            EPotLocalInter = EPotLocalInter + EPotLocal1
            dCosThetai = Rij4Inv * (10._RK * CosThetai * CosThetaj - CosGammaij)
            dCosThetaj = Rij4Inv * (5._RK * CosThetai2 - 1._RK)
            dCosGammaij = -2._RK * Rij4Inv * CosThetai
            Tmp = -4._RK * RijInv * EPotLocal1

            FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                      + (eX * CosThetaj - OXj) * dCosThetaj)
            FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                      + (eY * CosThetaj - OYj) * dCosThetaj)
            FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                      + (eZ * CosThetaj - OZj) * dCosThetaj)

            VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
            VirialLocalInter = VirialLocalInter + (FXij * PXij + FYij * PYij + FZij * PZij)
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
            d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Third*Third    !xxxx8   QD T

            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij

            forceTempX(jk) = forceTempX(jk) - FXij
            forceTempY(jk) = forceTempY(jk) - FYij
            forceTempZ(jk) = forceTempZ(jk) - FZij
            TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
            TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
            TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij

            momTempX(jk) = momTempX(jk) - eX * dCosThetaj - OXi * dCosGammaij
            momTempY(jk) = momTempY(jk) - eY * dCosThetaj - OYi * dCosGammaij  
            momTempZ(jk) = momTempZ(jk) - eZ * dCosThetaj - OZi * dCosGammaij   

#if  TRANS == 1
            !TRANSPORT_start
            VSxi   = VSxi + FXij * PYij
            VSyi   = VSyi + FXij * PZij
            VSzi   = VSzi + FYij * PZij
            VBxi   = VBxi + FXij * PXij
            VByi   = VByi + FYij * PYij
            VBzi   = VBzi + FZij * PZij
            VSuxi  = VSuxi+ FYij * PXij
            VSuyi  = VSuyi+ FZij * PXij
            VSuzi  = VSuzi+ FZij * PYij          
            UU     = EpotLocal1  
            Uxi     = UU * eX
            Uyi     = UU * eY
            Uzi     = UU * eZ
            Cxi    = Cxi  + Uxi
            Cyi    = Cyi  + Uyi
            Czi    = Czi  + Uzi
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
            !TRANSPORT_END
#endif
          end if
        end do loop1

        ! Include intramolecular interaction if need
        if (SameComponent .and. (this%potintra14 .or. this%potintra15)) then
          RXij = RXi - RX2(i)
          RYij = RYi - RY2(i)
          RZij = RZi - RZ2(i)
          PXij = PXi - PX2(i)
          PYij = PYi - PY2(i)
          PZij = PZi - PZ2(i)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          RijSquared = RXij**2 + RYij**2 + RZij**2
          OXj = OX2(i)
          OYj = OY2(i)
          OZj = OZ2(i)
#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosThetai2 = CosThetai**2
          CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal1 = Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) &
&                                  - CosGammaij * CosThetai)
          EPotLocal = EPotLocal + EPotLocal1*coeff
          EPotLocalIntra = EPotLocalIntra + EPotLocal1*coeff
          dCosThetai = Rij4Inv * (10._RK * CosThetai * CosThetaj - CosGammaij)
          dCosThetaj = Rij4Inv * (5._RK * CosThetai2 - 1._RK)
          dCosGammaij = -2._RK * Rij4Inv * CosThetai
          Tmp = -4._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FXij = FXij * coeff
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                   + (eY * CosThetaj - OYj) * dCosThetaj)
          FYij = FYij * coeff
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)
          FZij = FZij * coeff
          VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
          VirialLocalIntra = VirialLocalIntra + FXij * PXij + FYij * PYij + FZij * PZij
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + coeff * EPotLocal1*(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Third*Third    !xxxx8   QD T
          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(i) = forceTempX(i) - FXij
          forceTempY(i) = forceTempY(i) - FYij
          forceTempZ(i) = forceTempZ(i) - FZij
          TXi    = TXi    - (eX * dCosThetai + OXj * dCosGammaij) * coeff
          TYi    = TYi    - (eY * dCosThetai + OYj * dCosGammaij) * coeff
          TZi    = TZi    - (eZ * dCosThetai + OZj * dCosGammaij) * coeff

          momTempX(i) = momTempX(i) - (eX * dCosThetaj + OXi * dCosGammaij) * coeff
          momTempY(i) = momTempY(i) - (eY * dCosThetaj + OYi * dCosGammaij) * coeff
          momTempZ(i) = momTempZ(i) - (eZ * dCosThetaj + OZi * dCosGammaij) * coeff

        end if

        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi

#if  TRANS == 1
        !TRANSPORT_start
        VSx(i) = VSxi
        VSy(i) = VSyi
        VSz(i) = VSzi
        VBx(i) = VBxi
        VBy(i) = VByi
        VBz(i) = VBzi
        VSux(i)= VSuxi
        VSuy(i)= VSuyi
        VSuz(i)= VSuzi
        Cx(i)  = Cxi
        Cy(i)  = Cyi
        Cz(i)  = Czi
        tux(i) = tuxi
        tuy(i) = tuyi
        tuz(i) = tuzi
        tlx(i) = tlxi
        tly(i) = tlyi
        tlz(i) = tlzi
        tdx(i) = tdxi
        tdy(i) = tdyi
        tdz(i) = tdzi
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
          RijSquared = RXij**2 + RYij**2 + RZij**2
          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
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
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv*RijInv)*Third*Third    !xxxx8   QD ss T

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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if

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
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
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
    integer           :: nu1, nu2, unit, i0, jk
#if ARCH == 3
    logical           :: hit
#endif

    ! Assign local variables
    i1 = this%Site1%NTest
    j1 = this%Site2%NPart
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared
    nu1 = this%NUnit1
    nu2 = this%NUnit2
#if MPI_VER > 0
    i0 = this%Site1%NTest0
    i1 = this%Site1%NTest2
#else
    i0 = 1
    i1 = this%Site1%NTest
#endif

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
!$OMP PRIVATE (EPotLocal,i,i0,i1,j,k) 

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
!$OMP DO
      do i = i0, i1
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

#if ARCH == 3
        hit = .false.
#endif

!CDIR NODEP
        unit = nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then  ! choose only units, to which our Site2 correspond
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)   ! number of molecule, to which this unit correspond
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
            if( RijSquared <= RShieldSquared ) hit = .true.
#else
            if( RijSquared <= RShieldSquared ) then
              EPotLocal = 1E33_RK
              exit loop1
            end if
#endif
            OXj = OX2(jk)
            OYj = OY2(jk)
            OZj = OZ2(jk)

#if ARCH == 3
            RijInv = rsqrt( RijSquared )
#else
            RijInv = 1._RK / sqrt( RijSquared )
#endif
            eX = RXij * RijInv
            eY = RYij * RijInv
            eZ = RZij * RijInv
            CosThetai = OXi * eX + OYi * eY + OZi * eZ
            CosThetaj = OXj * eX + OYj * eY + OZj * eZ
            CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
            Rij4Inv = Epsilon / RijSquared**2
            EPotLocal = EPotLocal + Rij4Inv * (CosThetaj * (5._RK * CosThetai**2 - 1._RK) - 2._RK * CosGammaij * CosThetai)
          end if
        end do loop1

#if ARCH == 3
        if( .not. hit ) then
          EPotTest(i) = EPotTest(i) + EPotLocal
        else
          EPotTest(i) = EPotTest(i) + 1E33_RK
        endif
#else
        EPotTest(i) = EPotTest(i) + EPotLocal
#endif
      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over test particles
!$OMP DO
      do i = i0, i1
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        EPotLocal = 0._RK

#if ARCH == 3
        hit = .false.
#endif
!CDIR NODEP
loop2:  do j = 1, j1
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
          if( RijSquared <= RShieldSquared ) hit = .true.
#else
          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop2
          end if
#endif
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)
#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosGammaij = OXi * OXj + OYi * OYj + OZi * OZj
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal = EPotLocal + Rij4Inv * ( CosThetaj * (5._RK * CosThetai**2 - 1._RK) - 2._RK * CosGammaij * CosThetai )
        end do loop2

#if ARCH == 3
        if( .not. hit ) then
          EPotTest(i) = EPotTest(i) + EPotLocal
        else
          EPotTest(i) = EPotTest(i) + 1E33_RK
        endif
#else
        EPotTest(i) = EPotTest(i) + EPotLocal
#endif
      end do
!$OMP END DO
    end if
!$OMP END PARALLEL
  end subroutine TPotQD_ChemicalPotential



!==============================================================!
!  Subroutine TPotQD_Energy                                    !
!==============================================================!

  subroutine TPotQD_Energy( this, np, nu, F, E, EIntra, BoxLength, CompIdent )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleDipole) :: this
    integer, intent(in)        :: np
    integer, intent(in)        :: nu
    real(RK), intent(in out)   :: F(3,nu)
    real(RK), intent(in out)   :: E
    real(RK), intent(in out)   :: EIntra
    real(RK), intent(in)       :: BoxLength
    logical, intent(in)        :: CompIdent

    ! Declare local variables
    real(RK)          :: Epsilon
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
    real(RK)          :: CosThetai, CosThetaj, CosThetai2, CosGammaij
    real(RK)          :: EPot, EIntra1, EPotLocal, tempF(3,nu)
    integer           :: j, k, nu2, jk, unit
    real(RK)          :: coeff

    ! Assign local variables
    Epsilon = this%Epsilon
    RShieldSquared = this%RShieldSquared
    nu2 = this%NUnit2
    coeff = 1._RK
    if (this%potintra14) coeff = this%ScaleEl14
    EPot   = 0._RK
    EIntra1   = 0._RK

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
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

    ! Loop over molecules
    RXi = RX1(np)
    RYi = RY1(np)
    RZi = RZ1(np)
    OXi = OX1(np)
    OYi = OY1(np)
    OZi = OZ1(np)
    PXi = PX1(np)
    PYi = PY1(np)
    PZi = PZ1(np)

    unit=this%NUnit1*(np-1)+this%Site1%UnitNumber

    do k = 1, this%NInCutoff(unit)
      j = this%CutoffPartner(k, unit)
      if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
        if (mod(j,nu2)==0) then
          jk = INT(j/nu2)
        else
          jk = INT(j/nu2)+1
        end if
        RXij = RXi - RX2(jk)
        RYij = RYi - RY2(jk)
        RZij = RZi - RZ2(jk)
        PXij = PXi - PX2(jk)
        PYij = PYi - PY2(jk)
        PZij = PZi - PZ2(jk)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        RijSquared = RXij**2 + RYij**2 + RZij**2

        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
        else
          OXj = OX2(jk)
          OYj = OY2(jk)
          OZj = OZ2(jk)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
          eX = RXij * RijInv
          eY = RYij * RijInv
          eZ = RZij * RijInv
          CosThetai = OXi * eX + OYi * eY + OZi * eZ
          CosThetaj = OXj * eX + OYj * eY + OZj * eZ
          CosThetai2 = CosThetai**2
          CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
          Rij4Inv = Epsilon / RijSquared**2
          EPotLocal = Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) - CosGammaij * CosThetai)
        end if
        EPot = EPot + EPotLocal
      end if
    end do
    ! Include intramolecular interaction if need
    if (this%potintra14 .or. this%potintra15) then
      RXij = RXi - RX2(np)
      RYij = RYi - RY2(np)
      RZij = RZi - RZ2(np)
      PXij = PXi - PX2(np)
      PYij = PYi - PY2(np)
      PZij = PZi - PZ2(np)
      RXij = (RXij - anint( PXij )) * BoxLength
      RYij = (RYij - anint( PYij )) * BoxLength
      RZij = (RZij - anint( PZij )) * BoxLength
      RijSquared = RXij**2 + RYij**2 + RZij**2
      OXj = OX2(np)
      OYj = OY2(np)
      OZj = OZ2(np)

#if ARCH == 3
      RijInv = rsqrt( RijSquared )
#else
      RijInv = 1._RK / sqrt( RijSquared )
#endif
      eX = RXij * RijInv
      eY = RYij * RijInv
      eZ = RZij * RijInv
      CosThetai = OXi * eX + OYi * eY + OZi * eZ
      CosThetaj = OXj * eX + OYj * eY + OZj * eZ
      CosThetai2 = CosThetai**2
      CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
      Rij4Inv = Epsilon / RijSquared**2
      EPotLocal = coeff * Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) - CosGammaij * CosThetai)
      EIntra1 = EIntra1 + EPotLocal
    end if

    F(:,:) = F(:,:) + tempF(:,:)
    E = E + EPot + EIntra1
    EIntra = EIntra + EIntra1

  end subroutine TPotQD_Energy



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

    ! Declare local variables
    integer :: k, ende

    ! Construct potential
    this%Site1 => Molecule1%SiteQuadrupole(j1)
    this%Site2 => Molecule2%SiteQuadrupole(j2)
    this%NUnit1 => Molecule1%NUnit
    this%NUnit2 => Molecule2%NUnit
    this%SameComponent = i1 == i2
    this%Epsilon = .75_RK * this%Site1%Q * this%Site2%Q
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

    ! if this potential is intra
    if (this%SameComponent .and. Molecule1%hasIntraLJEl) then
      ende = size(Molecule1%IntQQ15(:,1))
      do k=1, ende
        if (Molecule1%IntQQ15(k,1)==this%Site1%SiteId .and. Molecule1%IntQQ15(k,2)==this%Site2%SiteId) then
          this%potintra15 = .true.
        end if
      end do
      if (LJEL14 .and. .not. this%potintra15) then
        ende = size(Molecule1%IntQQ14(:,1))
        do k=1, ende
          if (Molecule1%IntQQ14(k,1)==this%Site1%SiteId .and. Molecule1%IntQQ14(k,2)==this%Site2%SiteId) then
            this%potintra14=.true.
            this%ScaleEl14 = Molecule1%ScaleQQ14(k)
          end if
        end do
      end if
    else
      this%potintra15 = .false.
      this%potintra14 = .false.
    end if

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

  subroutine TPotQQ_Force( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleQuadrupole) :: this
    real(RK), intent(in out)       :: EPot
    real(RK), intent(in out)       :: Virial
    real(RK), intent(in out)       :: EPotInter
    real(RK), intent(in out)       :: VirialInter
    real(RK), intent(in out)       :: EPotIntra_Nonbonded
    real(RK), intent(in out)       :: VirialIntra
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
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
    integer           :: nu1, nu2, jk, unit
    real(RK)          :: coeff
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
!$OMP PRIVATE (i, j, k, i1, j0, j1) &
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

    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter   = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra   = 0._RK
    VirialLocalIntra = 0._RK

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

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

        unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            PXij = (PXij - anint( PXij )) * BoxLength
            PYij = (PYij - anint( PYij )) * BoxLength
            PZij = (PZij - anint( PZij )) * BoxLength
            RijSquared = RXij**2 + RYij**2 + RZij**2
            OXj = OX2(jk)
            OYj = OY2(jk)
            OZj = OZ2(jk)

#if ARCH == 3
            RijInv = rsqrt( RijSquared )
#else
            RijInv = 1._RK / sqrt( RijSquared )
#endif
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
&                        - 15._RK * CosThetaiSquared * CosThetajSquared + 2._RK * Tmp**2)
            EPotLocal = EPotLocal + EPotLocal1
            EPotLocalInter = EPotLocalInter + EPotLocaL1
            dCosThetai = Rij5Inv * (-10._RK * CosThetai - 30._RK * CosThetai * CosThetajSquared &
&                                 - 20._RK * CosThetaj * Tmp)
            dCosThetaj = Rij5Inv * (-10._RK * CosThetaj - 30._RK * CosThetaj * CosThetaiSquared &
&                                  - 20._RK * CosThetai * Tmp)

            dCosGammaij = 4._RK * Rij5Inv * Tmp
            Tmp = -5._RK * RijInv * EPotLocal1
            FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                      + (eX * CosThetaj - OXj) * dCosThetaj)
            FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                      + (eY * CosThetaj - OYj) * dCosThetaj)
            FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                      + (eZ * CosThetaj - OZj) * dCosThetaj)

            VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
            VirialLocalInter = VirialLocalInter + (FXij * PXij + FYij * PYij + FZij * PZij)
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
            d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(35._RK*sitecorr*sitecorr-5._RK*Plen2*RijInv*RijInv)*Third*Third      !xxxx9  QQ

            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij

            forceTempX(jk) = forceTempX(jk) - FXij
            forceTempY(jk) = forceTempY(jk) - FYij
            forceTempZ(jk) = forceTempZ(jk) - FZij

            TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
            TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
            TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
            momTempX(jk) = momTempX(jk) - eX * dCosThetaj - OXi * dCosGammaij
            momTempY(jk) = momTempY(jk) - eY * dCosThetaj - OYi * dCosGammaij  
            momTempZ(jk) = momTempZ(jk) - eZ * dCosThetaj - OZi * dCosGammaij

          end if
        end do loop1

        ! Include intramolecular interaction if need
        if (SameComponent .and. (this%potintra14 .or. this%potintra15)) then
          RXij = RXi - RX2(i)
          RYij = RYi - RY2(i)
          RZij = RZi - RZ2(i)
          PXij = PXi - PX2(i)
          PYij = PYi - PY2(i)
          PZij = PZi - PZ2(i)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          RijSquared = RXij**2 + RYij**2 + RZij**2
          OXj = OX2(i)
          OYj = OY2(i)
          OZj = OZ2(i)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
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

          EPotLocal = EPotLocal + EPotLocal1*coeff

          EPotLocalIntra = EPotLocalIntra + EPotLocal1*coeff

          dCosThetai = Rij5Inv * (-10._RK * CosThetai - 30._RK * CosThetai * CosThetajSquared &
&                                - 20._RK * CosThetaj * Tmp)
          dCosThetaj = Rij5Inv * (-10._RK * CosThetaj - 30._RK * CosThetaj * CosThetaiSquared &
&                                - 20._RK * CosThetai * Tmp)

          dCosGammaij = 4._RK * Rij5Inv * Tmp
          Tmp = -5._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FXij = FXij * coeff
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FYij = FYij * coeff
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)
          FZij = FZij * coeff
          VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
          VirialLocalIntra = VirialLocalIntra + (FXij * PXij + FYij * PYij + FZij * PZij)
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + coeff * EPotLocal1*(35._RK*sitecorr*sitecorr-5._RK*Plen2*RijInv*RijInv)*Third*Third      !xxxx9  QQ
          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(i) = forceTempX(i) - FXij
          forceTempY(i) = forceTempY(i) - FYij
          forceTempZ(i) = forceTempZ(i) - FZij

          TXi    = TXi    - (eX * dCosThetai + OXj * dCosGammaij) * coeff
          TYi    = TYi    - (eY * dCosThetai + OYj * dCosGammaij) * coeff
          TZi    = TZi    - (eZ * dCosThetai + OZj * dCosGammaij) * coeff

          momTempX(i) = momTempX(i) - (eX * dCosThetaj + OXi * dCosGammaij) * coeff
          momTempY(i) = momTempY(i) - (eY * dCosThetaj + OYi * dCosGammaij) * coeff
          momTempZ(i) = momTempZ(i) - (eZ * dCosThetaj + OZi * dCosGammaij) * coeff
        end if

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
          RijSquared = RXij**2 + RYij**2 + RZij**2
          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
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
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(35._RK*sitecorr*sitecorr-5._RK*Plen2*RijInv*RijInv)*Third*Third      !xxxx9  QQ ss

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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
    d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local

  end subroutine TPotQQ_Force


!==============================================================!
!  Subroutine TPotQQ_Force_Trans                               !
!==============================================================!

  subroutine TPotQQ_Force_Trans( this, EPot, Virial, EPotInter, VirialInter, &
&            EPotIntra_Nonbonded, VirialIntra, d2EpotdV2, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleQuadrupole) :: this
    real(RK), intent(in out)       :: EPot
    real(RK), intent(in out)       :: Virial
    real(RK), intent(in out)       :: EPotInter
    real(RK), intent(in out)       :: VirialInter
    real(RK), intent(in out)       :: EPotIntra_Nonbonded
    real(RK), intent(in out)       :: VirialIntra
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
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: EPotLocalInter, VirialLocalInter
    real(RK)          :: d2EpotdV2Local, sitecorr, Plen2
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
    integer           :: nu1, nu2, jk, unit
    real(RK)          :: coeff

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
    real(RK), pointer, contiguous :: VSx(:), VSy(:), VSz(:)
    real(RK), pointer, contiguous :: VSux(:),VSuy(:),VSuz(:)
    real(RK), pointer, contiguous :: VBx(:), VBy(:), VBz(:)
    real(RK), pointer, contiguous :: Cx(:) , Cy(:) , Cz(:)
    real(RK), pointer, contiguous :: tux(:) , tuy(:) , tuz(:)
    real(RK), pointer, contiguous :: tlx(:) , tly(:) , tlz(:)
    real(RK), pointer, contiguous :: tdx(:) , tdy(:) , tdz(:)
    real(RK), pointer, contiguous :: q1(:), q2(:), q3(:), q4(:)
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
    real(RK)          :: Uxi,  Uyi, Uzi
    real(RK)          :: A11, A12, A13, A21, A22, A23, A31, A32, A33
    !TRANSPORT_END
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
    i = 0   !why needed here,in QQ, but not in all other force-routines
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
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    EPotLocalInter   = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra   = 0._RK
    VirialLocalIntra = 0._RK

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

    if (this%potintra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._RK
    end if

#if  TRANS == 1
    !TRANSPORT_start
    VSx => this%Site1%vsQx
    VSy => this%Site1%vsQy
    VSz => this%Site1%vsQz
    VBx => this%Site1%vbQx
    VBy => this%Site1%vbQy
    VBz => this%Site1%vbQz
    VSux=> this%Site1%vsuQx
    VSuy=> this%Site1%vsuQy
    VSuz=> this%Site1%vsuQz
    Cx  => this%Site1%cQx
    Cy  => this%Site1%cQy
    Cz  => this%Site1%cQz
    tux => this%Site1%tuQx
    tuy => this%Site1%tuQy
    tuz => this%Site1%tuQz
    tlx => this%Site1%tlQx
    tly => this%Site1%tlQy
    tlz => this%Site1%tlQz
    tdx => this%Site1%tdQx
    tdy => this%Site1%tdQy
    tdz => this%Site1%tdQz
    q1  => this%Site1%Qm0r(:, 1, 1)
    q2  => this%Site1%Qm0r(:, 2, 1)
    q3  => this%Site1%Qm0r(:, 3, 1)
    q4  => this%Site1%Qm0r(:, 4, 1)
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
        VSxi= VSx(i)
        VSyi= VSy(i)
        VSzi= VSz(i)
        VBxi= VBx(i)
        VByi= VBy(i)
        VBzi= VBz(i)
        VSuxi= VSux(i)
        VSuyi= VSuy(i)
        VSuzi= VSuz(i)
        Cxi = Cx(i)
        Cyi = Cy(i)
        Czi = Cz(i)
        tuxi = tux(i)
        tuyi = tuy(i)
        tuzi = tuz(i)
        tlxi = tlx(i)
        tlyi = tly(i)
        tlzi = tlz(i)
        tdxi = tdx(i)
        tdyi = tdy(i)
        tdzi = tdz(i)
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

        unit=nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            PXij = (PXij - anint( PXij )) * BoxLength
            PYij = (PYij - anint( PYij )) * BoxLength
            PZij = (PZij - anint( PZij )) * BoxLength
            RijSquared = RXij**2 + RYij**2 + RZij**2
            OXj = OX2(jk)
            OYj = OY2(jk)
            OZj = OZ2(jk)

#if ARCH == 3
            RijInv = rsqrt( RijSquared )
#else
            RijInv = 1._RK / sqrt( RijSquared )
#endif
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
&                        - 15._RK * CosThetaiSquared * CosThetajSquared + 2._RK * Tmp**2)
            EPotLocal = EPotLocal + EPotLocal1
            EPotLocalInter = EPotLocalInter + EPotLocaL1
            dCosThetai = Rij5Inv * (-10._RK * CosThetai - 30._RK * CosThetai * CosThetajSquared &
&                                  - 20._RK * CosThetaj * Tmp)
            dCosThetaj = Rij5Inv * (-10._RK * CosThetaj - 30._RK * CosThetaj * CosThetaiSquared &
&                                - 20._RK * CosThetai * Tmp)
            dCosGammaij = 4._RK * Rij5Inv * Tmp
            Tmp = -5._RK * RijInv * EPotLocal1

            FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                      + (eX * CosThetaj - OXj) * dCosThetaj)
            FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                      + (eY * CosThetaj - OYj) * dCosThetaj)
            FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                      + (eZ * CosThetaj - OZj) * dCosThetaj)

            VirialLocal = VirialLocal + (FXij * PXij + FYij * PYij + FZij * PZij)
#if OSMOP == 2
loop2:    do m=1,NBinsDen
            if (PX2(jk) .ge. real(m-1)/NBinsDen-0.5_RK) then
              if (PX2(jk) < real(m)/NBinsDen-0.5_RK) then
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
            VirialLocalInter = VirialLocalInter + (FXij * PXij + FYij * PYij + FZij * PZij)
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
            d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(35._RK*sitecorr*sitecorr-5._RK*Plen2*RijInv*RijInv)*Third*Third      !xxxx9  QQ T

            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij

            forceTempX(jk) = forceTempX(jk) - FXij
            forceTempY(jk) = forceTempY(jk) - FYij
            forceTempZ(jk) = forceTempZ(jk) - FZij

            TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
            TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
            TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
            momTempX(jk) = momTempX(jk) - eX * dCosThetaj - OXi * dCosGammaij
            momTempY(jk) = momTempY(jk) - eY * dCosThetaj - OYi * dCosGammaij  
            momTempZ(jk) = momTempZ(jk) - eZ * dCosThetaj - OZi * dCosGammaij   

#if  TRANS == 1
            !TRANSPORT_start
            VSxi   = VSxi + FXij * PYij
            VSyi   = VSyi + FXij * PZij
            VSzi   = VSzi + FYij * PZij
            VBxi   = VBxi + FXij * PXij
            VByi   = VByi + FYij * PYij
            VBzi   = VBzi + FZij * PZij
            VSuxi  = VSuxi+ FYij * PXij
            VSuyi  = VSuyi+ FZij * PXij
            VSuzi  = VSuzi+ FZij * PYij
            Uxi     = EPotLocal1 * eX
            Uyi     = EPotLocal1 * eY
            Uzi     = EPotLocal1 * eZ
            FTXi   = - eX * dCosThetai - OXj * dCosGammaij
            FTYi   = - eY * dCosThetai - OYj * dCosGammaij
            FTZi   = - eZ * dCosThetai - OZj * dCosGammaij
            Cxi    = Cxi  + Uxi
            Cyi    = Cyi  + Uyi
            Czi    = Czi  + Uzi
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
            !TRANSPORT_END
#endif
          end if
        end do loop1

        ! Include intramolecular interaction if need
        if (SameComponent .and. (this%potintra14 .or. this%potintra15)) then
          RXij = RXi - RX2(i)
          RYij = RYi - RY2(i)
          RZij = RZi - RZ2(i)
          PXij = PXi - PX2(i)
          PYij = PYi - PY2(i)
          PZij = PZi - PZ2(i)
          RXij = (RXij - anint( PXij )) * BoxLength
          RYij = (RYij - anint( PYij )) * BoxLength
          RZij = (RZij - anint( PZij )) * BoxLength
          PXij = (PXij - anint( PXij )) * BoxLength
          PYij = (PYij - anint( PYij )) * BoxLength
          PZij = (PZij - anint( PZij )) * BoxLength
          RijSquared = RXij**2 + RYij**2 + RZij**2
          OXj = OX2(i)
          OYj = OY2(i)
          OZj = OZ2(i)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
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

          EPotLocal = EPotLocal + EPotLocal1*coeff
          EPotLocalIntra = EPotLocalIntra + EPotLocal1*coeff
          dCosThetai = Rij5Inv * (-10._RK * CosThetai - 30._RK * CosThetai * CosThetajSquared &
&                                - 20._RK * CosThetaj * Tmp)

          dCosThetaj = Rij5Inv * (-10._RK * CosThetaj - 30._RK * CosThetaj * CosThetaiSquared &
&                                - 20._RK * CosThetai * Tmp)

          dCosGammaij = 4._RK * Rij5Inv * Tmp
          Tmp = -5._RK * RijInv * EPotLocal1

          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FXij = FXij * coeff
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FYij = FYij * coeff
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)
          FZij = FZij * coeff
          VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
          VirialLocalIntra = VirialLocalIntra + FXij * PXij + FYij * PYij + FZij * PZij
          Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
          sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv*RijInv
          d2EpotdV2Local = d2EpotdV2Local + coeff * EPotLocal1*(35._RK*sitecorr*sitecorr-5._RK*Plen2*RijInv*RijInv)*Third*Third      !xxxx9  QQ T

          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij

          forceTempX(i) = forceTempX(i) - FXij
          forceTempY(i) = forceTempY(i) - FYij
          forceTempZ(i) = forceTempZ(i) - FZij

          TXi    = TXi    - (eX * dCosThetai + OXj * dCosGammaij) * coeff
          TYi    = TYi    - (eY * dCosThetai + OYj * dCosGammaij) * coeff
          TZi    = TZi    - (eZ * dCosThetai + OZj * dCosGammaij) * coeff
          momTempX(i) = momTempX(i) - (eX * dCosThetaj + OXi * dCosGammaij) * coeff
          momTempY(i) = momTempY(i) - (eY * dCosThetaj + OYi * dCosGammaij) * coeff
          momTempZ(i) = momTempZ(i) - (eZ * dCosThetaj + OZi * dCosGammaij) * coeff

        end if

        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi

#if  TRANS == 1
        !TRANSPORT_start
        VSx(i) = VSxi
        VSy(i) = VSyi
        VSz(i) = VSzi
        VBx(i) = VBxi
        VBy(i) = VByi
        VBz(i) = VBzi
        VSux(i)= VSuxi
        VSuy(i)= VSuyi
        VSuz(i)= VSuzi
        Cx(i)  = Cxi
        Cy(i)  = Cyi
        Cz(i)  = Czi
        tux(i) = tuxi
        tuy(i) = tuyi
        tuz(i) = tuzi
        tlx(i) = tlxi
        tly(i) = tlyi
        tlz(i) = tlzi
        tdx(i) = tdxi
        tdy(i) = tdyi
        tdz(i) = tdzi
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
          RijSquared = RXij**2 + RYij**2 + RZij**2

          if( RijSquared >= RCutoffSquared ) cycle loop3
          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
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
          d2EpotdV2Local = d2EpotdV2Local + EPotLocal1*(35._RK*sitecorr*sitecorr-5._RK*Plen2*RijInv*RijInv)*Third*Third      !xxxx9  QQ ss T

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
    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter
    if (IntraLJEl) then
      EPotIntra_Nonbonded = EPotIntra_Nonbonded + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if
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
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
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
    integer           :: nu1, nu2, unit, i0, jk
#if ARCH == 3
    logical           :: hit
#endif

    ! Assign local variables
    i1 = this%Site1%NTest
    j1 = this%Site2%NPart
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared
    nu1 = this%NUnit1
    nu2 = this%NUnit2
#if MPI_VER > 0
    i0 = this%Site1%NTest0
    i1 = this%Site1%NTest2
#else
    i0 = 1
    i1 = this%Site1%NTest
#endif

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
!$OMP PRIVATE (EPotLocal,i,i0,i1,j,k) 

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
!$OMP DO 
      do i = i0, i1
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

#if ARCH == 3
        hit = .false.
#endif

!CDIR NODEP

        unit = nu1*(i-1)+this%Site1%UnitNumber

loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then  ! choose only units, to which our Site2 correspond
            if (mod(j,nu2)==0) then
              jk = INT(j/nu2)   ! number of molecule, to which this unit correspond
            else
              jk = INT(j/nu2)+1
            end if
            RXij = RXi - RX2(jk)
            RYij = RYi - RY2(jk)
            RZij = RZi - RZ2(jk)
            PXij = PXi - PX2(jk)
            PYij = PYi - PY2(jk)
            PZij = PZi - PZ2(jk)
            RXij = (RXij - anint( PXij )) * BoxLength
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
            if( RijSquared <= RShieldSquared ) hit = .true.
#else
            if( RijSquared <= RShieldSquared ) then
              EPotLocal = 1E33_RK
              exit loop1
            end if
#endif

            OXj = OX2(jk)
            OYj = OY2(jk)
            OZj = OZ2(jk)

#if ARCH == 3
            RijInv = rsqrt( RijSquared )
#else
            RijInv = 1._RK / sqrt( RijSquared )
#endif

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
&                       - 15._RK * CosThetaiSquared * CosThetajSquared + 2._RK * Tmp**2)
          end if
        end do loop1

#if ARCH == 3
        if( .not. hit ) then
          EPotTest(i) = EPotTest(i) + EPotLocal
        else
          EPotTest(i) = EPotTest(i) + 1E33_RK
        endif
#else
        EPotTest(i) = EPotTest(i) + EPotLocal
#endif
      end do
!$OMP END DO
    else ! Site-site cutoff

      ! Loop over test particles
!$OMP DO
      do i = i0, i1
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        OXi = OX1(i)
        OYi = OY1(i)
        OZi = OZ1(i)
        EPotLocal = 0._RK
#if ARCH == 3
        hit = .false.
#endif
!CDIR NODEP
loop2:  do j = 1, j1
          RXij = RXi - RX2(j)
          RYij = RYi - RY2(j)
          RZij = RZi - RZ2(j)
          RXij = (RXij - anint( RXij )) * BoxLength
          RYij = (RYij - anint( RYij )) * BoxLength
          RZij = (RZij - anint( RZij )) * BoxLength
          RijSquared = RXij**2 + RYij**2 + RZij**2

#if ARCH == 3
          if( RijSquared <= RShieldSquared ) hit = .true.
#else
          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop2
          end if
#endif

          OXj = OX2(j)
          OYj = OY2(j)
          OZj = OZ2(j)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
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

#if ARCH == 3
        if( .not. hit ) then
          EPotTest(i) = EPotTest(i) + EPotLocal
        else
          EPotTest(i) = EPotTest(i) + 1E33_RK
        endif
#else
        EPotTest(i) = EPotTest(i) + EPotLocal
#endif
      end do
!$OMP END DO
    end if
!$OMP END PARALLEL

  end subroutine TPotQQ_ChemicalPotential



!==============================================================!
!  Subroutine TPotQQ_Energy                                    !
!==============================================================!

  subroutine TPotQQ_Energy( this, np, nu, F, E, EIntra, BoxLength, CompIdent )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleQuadrupole) :: this
    integer, intent(in)            :: np
    integer, intent(in)            :: nu
    real(RK), intent(in out)       :: F(3,nu)
    real(RK), intent(in out)       :: E
    real(RK), intent(in out)       :: EIntra
    real(RK), intent(in)           :: BoxLength
    logical, intent(in)            :: CompIdent

    ! Declare local variables
    real(RK)          :: Epsilon
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
    real(RK)          :: EPot, EIntra1, EPotLocal, tempF(3,nu)
    integer           :: j, k, nu2, jk, unit
    real(RK)          :: coeff

    ! Assign local variables
    nu2 = this%NUnit2
    Epsilon = this%Epsilon
    RShieldSquared = this%RShieldSquared
    coeff = 1._RK
    if (this%potintra14) coeff = this%ScaleEl14
    EPot   = 0._RK
    EIntra1   = 0._RK

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
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

    ! Loop over molecules
    RXi = RX1(np)
    RYi = RY1(np)
    RZi = RZ1(np)
    OXi = OX1(np)
    OYi = OY1(np)
    OZi = OZ1(np)
    PXi = PX1(np)
    PYi = PY1(np)
    PZi = PZ1(np)

    unit=this%NUnit1*(np-1)+this%Site1%UnitNumber

    do k = 1, this%NInCutoff(unit)
      j = this%CutoffPartner(k, unit)
      if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
        if (mod(j,nu2)==0) then
          jk = INT(j/nu2)
        else
          jk = INT(j/nu2)+1
        end if
        RXij = RXi - RX2(jk)
        RYij = RYi - RY2(jk)
        RZij = RZi - RZ2(jk)
        PXij = PXi - PX2(jk)
        PYij = PYi - PY2(jk)
        PZij = PZi - PZ2(jk)
        RXij = (RXij - anint( PXij )) * BoxLength
        RYij = (RYij - anint( PYij )) * BoxLength
        RZij = (RZij - anint( PZij )) * BoxLength
        RijSquared = RXij**2 + RYij**2 + RZij**2

        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
        else
          OXj = OX2(jk)
          OYj = OY2(jk)
          OZj = OZ2(jk)

#if ARCH == 3
          RijInv = rsqrt( RijSquared )
#else
          RijInv = 1._RK / sqrt( RijSquared )
#endif
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
          EPotLocal = Rij5Inv * (1._RK - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&           - 15._RK * CosThetaiSquared * CosThetajSquared + 2._RK * Tmp**2)
        end if
        EPot = EPot + EPotLocal
      end if
    end do
    ! Include intramolecular interaction if need
    if (this%potintra15 .or. this%potintra14) then
      RXij = RXi - RX2(np)
      RYij = RYi - RY2(np)
      RZij = RZi - RZ2(np)
      PXij = PXi - PX2(np)
      PYij = PYi - PY2(np)
      PZij = PZi - PZ2(np)
      RXij = (RXij - anint( PXij )) * BoxLength
      RYij = (RYij - anint( PYij )) * BoxLength
      RZij = (RZij - anint( PZij )) * BoxLength
      RijSquared = RXij**2 + RYij**2 + RZij**2
      OXj = OX2(np)
      OYj = OY2(np)
      OZj = OZ2(np)
#if ARCH == 3
      RijInv = rsqrt( RijSquared )
#else
      RijInv = 1._RK / sqrt( RijSquared )
#endif
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
      EPotLocal = coeff * Rij5Inv * (1._RK - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&                    - 15._RK * CosThetaiSquared * CosThetajSquared + 2._RK * Tmp**2)
      EIntra1 = EIntra1 + EPotLocal
    end if

    F(:,:) = F(:,:) + tempF(:,:)
    E = E + EPot + EIntra1
    EIntra = EIntra + EIntra1

  end subroutine TPotQQ_Energy



  subroutine TPoterfc_approx(this,in,approx_out)

  type(TPotChargeCharge)      :: this
  real(RK),intent(in)     :: in
  real(RK),intent(out)    :: approx_out

! Local variables
  real(RK)                :: argu,C1,C2,C3,C4,C5,P

  C1 =  0.254829592
  C2 = -0.284496736
  C3 =  1.421413741
  C4 = -1.453152027
  C5 =  1.061405429
  P  =  0.3275911

  argu = 1._RK / (1._RK + P*in)
  approx_out = argu*(C1+argu*(C2+argu*(C3+argu*(C4+argu*C5))))*exp(-in**2)

  end subroutine TPoterfc_approx



!==============================================================!
!  Subroutine TPotBond_Construct                               !
!==============================================================!

  subroutine TPotBond_Construct( this, Molecule, j )


    implicit none

    ! Declare arguments
    type(TPotBond)              :: this
    type(TMolecule), intent(in) :: Molecule
    integer, intent(in)         :: j

    ! Construct potential

    this%Bond => Molecule%IdfBond(j)
    this%Site1 = this%Bond%SiteId1
    this%Site2 = this%Bond%SiteId2
    this%Unit1 = this%Bond%UnitId1
    this%Unit2 = this%Bond%UnitId2
    this%ForConst => this%Bond%ForConst
    this%R0 = this%Bond%R0

  end subroutine TPotBond_Construct


!==============================================================!
!  Subroutine TPotBond_Destruct                                !
!==============================================================!

  subroutine TPotBond_Destruct( this )

    implicit none

    ! Declare arguments
    type(TPotBond) :: this

    ! Destroy potential
    continue

  end subroutine TPotBond_Destruct


!==============================================================!
!  Subroutine TPotBond_Force                                   !
!==============================================================!

  subroutine TPotBond_Force( this, EPot, Virial, EPotIntra_Bond, VirialIntra, d2EpotdV2, BoxLength )

    use math_types

    implicit none

    ! Declare arguments
    type(TPotBond)     :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotIntra_Bond
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK)          :: R, RSquared
    real(RK)          :: Fijabs
    real(RK)          :: EPotLocal, VirialLocal
    real(RK)          :: d2EpotdV2Local , Plen2, sitecorr
    real(RK)          :: dR, F0, R0, ForConst
    type(vector)      :: P1, P2, R1, R2, Rij, Pij, Fij

    integer           :: i, i1
#if MPI_VER > 0
    integer           :: i0
#endif

#if MPI_VER > 0
     i0 = this%Bond%NPart0
     i1 = this%Bond%NPart2
#else
     i1 = this%Bond%NPart
#endif

     ForConst = this%ForConst
!     R0 = this%R0/BoxLength

     R0 = this%R0
     EPotLocal   = 0._RK
     VirialLocal = 0._RK
     d2EpotdV2Local = 0._RK

    ! Assign pointers
     FX1 => this%Bond%FX1
     FY1 => this%Bond%FY1
     FZ1 => this%Bond%FZ1
     FX2 => this%Bond%FX2
     FY2 => this%Bond%FY2
     FZ2 => this%Bond%FZ2

    ! Loop over molecules
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif

!CDIR NODEP

        R1 = vector(this%Bond%RX1(i), this%Bond%RY1(i), this%Bond%RZ1(i))
        R2 = vector(this%Bond%RX2(i), this%Bond%RY2(i), this%Bond%RZ2(i))

        ! Standard harmonic bond
        ! Energy and forces:
        ! formulae  E = ForConst*(R - R0)**2
        !           F = - 2*ForConst*(R-R0)/R - abs. value

        ! Calculate bond length
        Rij = SUB_VECTOR(R1, R2)
        Rij = SUB_ANINT_VECTOR(Rij)
        Rij = SCALE_VECTOR(Rij, BoxLength)

        RSquared=SQR_VECTOR_NORM(Rij)
        R=sqrt(RSquared) ! Bond length

        ! Deviation from equilibrium
        dR=R-R0

        ! Potential parameter
        F0 = dR*ForConst

        ! Energy of the bond
        EPotLocal = EPotLocal + dR*F0

        ! Force (abs. value)
        Fijabs=-2.0d0*F0/R

        ! Force components
        Fij = SCALE_VECTOR(Rij, Fijabs)

        ! For calculation of virial
        P1 = vector(this%Bond%PX1(i), this%Bond%PY1(i), this%Bond%PZ1(i))
        P2 = vector(this%Bond%PX2(i), this%Bond%PY2(i), this%Bond%PZ2(i))

        Pij = SUB_VECTOR(P1, P2)
        Pij = SUB_ANINT_VECTOR(Pij)
        Pij = SCALE_VECTOR(Pij, BoxLength)

        ! Contribution to virial
        VirialLocal = VirialLocal + SCALAR_PRODUCT(Pij, Fij)
        Plen2    = SQR_VECTOR_NORM(Pij)
        sitecorr = SCALAR_PRODUCT(Pij, Rij) / RSquared
        d2EpotdV2Local = d2EpotdV2Local - R * 2._RK * ForConst * dR * (sitecorr * sitecorr - Plen2/RSquared)*Third*Third !xxxx Bond CC
        d2EpotdV2Local = d2EpotdV2Local + RSquared * 2._RK * ForConst * sitecorr * sitecorr *Third*Third

         ! New Forces
         FX1(i) = FX1(i) + Fij%x
         FY1(i) = FY1(i) + Fij%y
         FZ1(i) = FZ1(i) + Fij%z
         FX2(i) = FX2(i) - Fij%x
         FY2(i) = FY2(i) - Fij%y
         FZ2(i) = FZ2(i) - Fij%z

       end do

     ! Update potential energy and virial
       EPot = EPot + EPotLocal
       Virial = Virial + Third * VirialLocal

     ! Update Intra potential energy and virial
       d2EpotdV2 = d2EpotdV2 + d2EpotdV2Local
       EPotIntra_Bond = EPotIntra_Bond + EPotLocal
       VirialIntra = VirialIntra + Third * VirialLocal

  end subroutine TPotBond_Force


!==============================================================!
!  Subroutine TPotBond_Energy                                  !
!==============================================================!

  subroutine TPotBond_Energy( this, np, nu, F, EBond, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotBond)           :: this
    integer, intent(in)      :: np
    integer, intent(in)      :: nu
    real(RK), intent(in out) :: F(3,nu)
    real(RK), intent(in out) :: EBond
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK)          :: dR, R0, R, RSquared
    real(RK)          :: ForConst, F0, Fij, FXij, FYij, FZij
    real(RK)          :: RXij, RYij, RZij, PXij, PYij, PZij
    real(RK)          :: EPotLocal

    ForConst = this%ForConst
    R0 = this%R0

    ! Calculate bond length
    RXij = this%Bond%RX1(np) - this%Bond%RX2(np)
    RYij = this%Bond%RY1(np) - this%Bond%RY2(np)
    RZij = this%Bond%RZ1(np) - this%Bond%RZ2(np)
    RXij = (RXij - anint( RXij )) * BoxLength
    RYij = (RYij - anint( RYij )) * BoxLength
    RZij = (RZij - anint( RZij )) * BoxLength

    RSquared=RXij**2+RYij**2+RZij**2
    R=sqrt(RSquared) ! Bond length
    ! Deviation from equilibrium
    dR=R-R0
    ! Potential parameter
    F0 = dR*ForConst
    ! Energy of the bond
    EBond = EBond + dR*F0
    ! Force (abs. value)
    Fij  = -2.0d0*F0/R
    FXij = Fij * RXij
    FYij = Fij * RYij
    FZij = Fij * RZij
    F(1,this%Bond%SiteId1) = F(1,this%Bond%SiteId1) + FXij
    F(2,this%Bond%SiteId1) = F(2,this%Bond%SiteId1) + FYij
    F(3,this%Bond%SiteId1) = F(3,this%Bond%SiteId1) + FZij
    F(1,this%Bond%SiteId2) = F(1,this%Bond%SiteId2) - FXij
    F(2,this%Bond%SiteId2) = F(2,this%Bond%SiteId2) - FYij
    F(3,this%Bond%SiteId2) = F(3,this%Bond%SiteId2) - FZij

  end subroutine TPotBond_Energy


!==============================================================!
!  Subroutine TPotAngle_Construct                              !
!==============================================================!

  subroutine TPotAngle_Construct( this, Molecule, j )


    implicit none

    ! Declare arguments
    type(TPotAngle)              :: this
    type(TMolecule), intent(in) :: Molecule
    integer, intent(in)         :: j

    ! Construct potential

    this%Angle => Molecule%IdfAngle(j)
    this%Site1 = this%Angle%SiteId1
    this%Site2 = this%Angle%SiteId2
    this%Site3 = this%Angle%SiteId3
    this%Unit1 = this%Angle%UnitId1
    this%Unit2 = this%Angle%UnitId2
    this%Unit3 = this%Angle%UnitId3
    this%ForConst => this%Angle%ForConst
    this%Angle0 = this%Angle%Angle0
    this%orientation1 => this%Angle%orientation1
    this%orientation2 => this%Angle%orientation2

  end subroutine TPotAngle_Construct


!==============================================================!
!  Subroutine TPotAngle_Destruct                                !
!==============================================================!

  subroutine TPotAngle_Destruct( this )

    implicit none

    ! Declare arguments
    type(TPotAngle) :: this

    ! Destroy potential
    continue

  end subroutine TPotAngle_Destruct


!==============================================================!
!  Subroutine TPotAngle_Force                                  !
!==============================================================!

  subroutine TPotAngle_Force( this, EPot, EPotIntra_Angle, BoxLength )

    use math_types

    implicit none

    ! Declare arguments
    type(TPotAngle)     :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: EPotIntra_Angle
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:), FX3(:), FY3(:), FZ3(:)
    real(RK)          :: RijRkj, RijSquared, RkjSquared
    type(vector)      :: R1, R2, R3, Rij, Rkj
    real(RK)          :: EPotLocal
    real(RK)          :: ForConst, Angle, Angle0, dAngle, cosa, sina
    real(RK)          :: abc, sab, cab, fab, fbb, faa, fax, fay, faz,  fbx, fby, fbz


    integer           :: i, i1
#if MPI_VER > 0
    integer           :: i0
#endif

#if MPI_VER > 0
     i0 = this%Angle%NPart0
     i1 = this%Angle%NPart2
#else
     i1 = this%Angle%NPart
#endif

     ForConst = this%ForConst
     Angle0 = this%Angle0
     EPotLocal = 0._RK

    ! Assign pointers
     FX1 => this%Angle%FX1 !           (1)    (3)
     FY1 => this%Angle%FY1 !             \    /
     FZ1 => this%Angle%FZ1 !            a \  / b
     FX2 => this%Angle%FX2 !               \/
     FY2 => this%Angle%FY2 !               (2)
     FZ2 => this%Angle%FZ2
     FX3 => this%Angle%FX3
     FY3 => this%Angle%FY3
     FZ3 => this%Angle%FZ3

    ! Standard harmonic potential of angle
    ! Energy:
    ! formulae  E = ForConst*(a - a0)**2

      ! Loop over molecules
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif

!CDIR NODEP
         R1 = vector(this%Angle%RX1(i), this%Angle%RY1(i), this%Angle%RZ1(i))
         R2 = vector(this%Angle%RX2(i), this%Angle%RY2(i), this%Angle%RZ2(i))
         R3 = vector(this%Angle%RX3(i), this%Angle%RY3(i), this%Angle%RZ3(i))

         Rij = SUB_VECTOR(R1, R2)
         Rkj = SUB_VECTOR(R3, R2)

         Rij = SUB_ANINT_VECTOR(Rij)
         Rij = SCALE_VECTOR(Rij, BoxLength)

         Rkj = SUB_ANINT_VECTOR(Rkj)
         Rkj = SCALE_VECTOR(Rkj, BoxLength)

         RijSquared = SQR_VECTOR_NORM(Rij)
         RkjSquared = SQR_VECTOR_NORM(Rkj)

         ! Calculate angle
         RijRkj = sqrt(RijSquared*RkjSquared)
         cosa = SCALAR_PRODUCT(Rij, Rkj) / RijRkj
         if( cosa .gt. 1._RK ) cosa = 1._RK
         if( cosa .lt.  -1._RK ) cosa = -1._RK
         Angle = acos(cosa)

         ! Calculate energy
         ! Deviation from equilibrium
         dAngle = Angle - Angle0

         ! Derivative of the energy
         abc = dAngle*ForConst
         EPotLocal = EPotLocal + abc*dAngle

         ! Force calculation
         sina = sqrt(1._RK-cosa**2)
         if( sina .lt. 1E-12_RK ) sina = 1E-12_RK
         sab = -2._RK*abc/sina
         cab = sab*cosa

         fab = sab/RijRkj
         faa = cab/RijSquared
         fbb = cab/RkjSquared

         fax = fab*Rkj%x-faa*Rij%x
         fay = fab*Rkj%y-faa*Rij%y
         faz = fab*Rkj%z-faa*Rij%z

         fbx = fab*Rij%x-fbb*Rkj%x
         fby = fab*Rij%y-fbb*Rkj%y
         fbz = fab*Rij%z-fbb*Rkj%z

         FX1(i) = FX1(i) - fax
         FY1(i) = FY1(i) - fay
         FZ1(i) = FZ1(i) - faz
         FX2(i) = FX2(i) + fax + fbx
         FY2(i) = FY2(i) + fay + fby
         FZ2(i) = FZ2(i) + faz + fbz
         FX3(i) = FX3(i) - fbx
         FY3(i) = FY3(i) - fby
         FZ3(i) = FZ3(i) - fbz
       end do

    ! Update potential energy, no contribution to virial!
     EPot = EPot + EPotLocal

    ! Update potential energy, no contribution to virial!
     EPotIntra_Angle = EPotIntra_Angle + EPotLocal

  end subroutine TPotAngle_Force


!==============================================================!
!  Subroutine TPotAngle_Energy                                 !
!==============================================================!

  subroutine TPotAngle_Energy( this, np, nu, F, EAngle, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotAngle)          :: this
    integer, intent(in)      :: np
    integer, intent(in)      :: nu
    real(RK), intent(in out) :: F(3,nu)
    real(RK), intent(in out) :: EAngle
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK)          :: RijRkj, RijSquared, RkjSquared
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: RXkj, RYkj, RZkj
    real(RK)          :: EPotLocal
    real(RK)          :: ForConst, Angle, Angle0, dAngle
    real(RK)          :: abc, cosa, sina, sab, cab
    real(RK)          :: faa, fab, fbb, fax, fay, faz, fbx, fby, fbz

    ForConst = this%ForConst
    Angle0 = this%Angle0

    RXij = this%Angle%RX1(np) - this%Angle%RX2(np)
    RYij = this%Angle%RY1(np) - this%Angle%RY2(np)
    RZij = this%Angle%RZ1(np) - this%Angle%RZ2(np)
    RXkj = this%Angle%RX3(np) - this%Angle%RX2(np)
    RYkj = this%Angle%RY3(np) - this%Angle%RY2(np)
    RZkj = this%Angle%RZ3(np) - this%Angle%RZ2(np)
    RXij = (RXij - anint( RXij )) * BoxLength
    RYij = (RYij - anint( RYij )) * BoxLength
    RZij = (RZij - anint( RZij )) * BoxLength
    RXkj = (RXkj - anint( RXkj )) * BoxLength
    RYkj = (RYkj - anint( RYkj )) * BoxLength
    RZkj = (RZkj - anint( RZkj )) * BoxLength

    RijSquared=RXij**2+RYij**2+RZij**2
    RkjSquared=RXkj**2+RYkj**2+RZkj**2

    ! Calculate angle
    RijRkj=sqrt(RijSquared*RkjSquared)
    cosa = (RXij*RXkj+RYij*RYkj+RZij*RZkj)/RijRkj
    if( cosa .gt. 1._RK ) cosa = 1._RK
    if( cosa .lt.  -1._RK ) cosa = -1._RK
    Angle = acos(cosa)

    ! Deviation from equilibrium
    dAngle = Angle - Angle0
    ! Derivative of the energy
    abc = dAngle*ForConst
    EAngle = EAngle + abc*dAngle
    ! Force calculation
    sina = sqrt(1._RK-cosa**2)
    if( sina .lt. 1E-12_RK ) sina = 1E-12_RK
    sab = -2._RK*abc/sina
    cab = sab*cosa

    fab = sab/RijRkj
    faa = cab/RijSquared
    fbb = cab/RkjSquared

    fax = fab*RXkj-faa*RXij
    fay = fab*RYkj-faa*RYij
    faz = fab*RZkj-faa*RZij

    fbx = fab*RXij-fbb*RXkj
    fby = fab*RYij-fbb*RYkj
    fbz = fab*RZij-fbb*RZkj

    F(1,this%Angle%SiteId1) = F(1,this%Angle%SiteId1) - fax
    F(2,this%Angle%SiteId1) = F(2,this%Angle%SiteId1) - fay
    F(3,this%Angle%SiteId1) = F(3,this%Angle%SiteId1) - faz
    F(1,this%Angle%SiteId2) = F(1,this%Angle%SiteId2) + fax + fbx
    F(2,this%Angle%SiteId2) = F(2,this%Angle%SiteId2) + fay + fby
    F(2,this%Angle%SiteId2) = F(3,this%Angle%SiteId2) + faz + fbz
    F(1,this%Angle%SiteId3) = F(1,this%Angle%SiteId3) - fbx
    F(2,this%Angle%SiteId3) = F(2,this%Angle%SiteId3) - fby
    F(3,this%Angle%SiteId3) = F(3,this%Angle%SiteId3) - fbz

  end subroutine TPotAngle_Energy


!==============================================================!
!  Subroutine TPotDihedral_Construct                           !
!==============================================================!

  subroutine TPotDihedral_Construct( this, Molecule, j )


    implicit none

    ! Declare arguments
    type(TPotDihedral)          :: this
    type(TMolecule), intent(in) :: Molecule
    integer, intent(in)         :: j

    ! Construct potential

    this%Dihedral => Molecule%IdfDihedral(j)
    this%Site1 = this%Dihedral%SiteId1
    this%Site2 = this%Dihedral%SiteId2
    this%Site3 = this%Dihedral%SiteId3
    this%Site4 = this%Dihedral%SiteId4
    this%Unit1 = this%Dihedral%UnitId1
    this%Unit2 = this%Dihedral%UnitId2
    this%Unit3 = this%Dihedral%UnitId3
    this%Unit4 = this%Dihedral%UnitId4
    this%nmax = this%Dihedral%nmax
    this%ForConst => this%Dihedral%ForConst
    this%gamma0 => this%Dihedral%gamma0

  end subroutine TPotDihedral_Construct


!==============================================================!
!  Subroutine TPotDihedral_Destruct                            !
!==============================================================!

  subroutine TPotDihedral_Destruct( this )

    implicit none

    ! Declare arguments
    type(TPotDihedral) :: this

    ! Destroy potential
    continue

  end subroutine TPotDihedral_Destruct


!==============================================================!
!  Subroutine TPotDihedral_Force                               !
!==============================================================!

  subroutine TPotDihedral_Force( this, EPot, EPotIntra_Dihedral,  BoxLength )

    use math_types

    implicit none

    ! Declare arguments
    type(TPotDihedral)     :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: EPotIntra_Dihedral
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK), pointer, contiguous :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:), FX3(:), FY3(:), FZ3(:), FX4(:), FY4(:), FZ4(:)
    real(RK)          :: EPotLocal, VirialLocal
    real(RK)          :: num, den, de1
    real(RK)          :: ab, bc, ac, aa, bb, cc, axb, bxc, co, si, signum, arg, earg
    real(RK)          :: deri,dnum,dden,ffi,ffj,ffk,ffl
    type(vector)      :: a, b, c, Rj, Ri, Rm, Rl

    integer           :: i, i1, j
#if MPI_VER > 0
     integer           :: i0
#endif

#if MPI_VER > 0
    i0 = this%Dihedral%NPart0
    i1 = this%Dihedral%NPart2
#else
     i1 = this%Dihedral%NPart
#endif

    EPotLocal   = 0._RK
    VirialLocal = 0._RK

    ! Assign pointers
     FX1 => this%Dihedral%FX1
     FY1 => this%Dihedral%FY1
     FZ1 => this%Dihedral%FZ1
     FX2 => this%Dihedral%FX2 !                  (i)            (l)
     FY2 => this%Dihedral%FY2 !                    \            /
     FZ2 => this%Dihedral%FZ2 !                  a  \          / c
     FX3 => this%Dihedral%FX3 !                      (j)-----(m)
     FY3 => this%Dihedral%FY3 !                            b
     FZ3 => this%Dihedral%FZ3
     FX4 => this%Dihedral%FX4
     FY4 => this%Dihedral%FY4
     FZ4 => this%Dihedral%FZ4

      ! Loop over molecules
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        Ri = vector(this%Dihedral%RX1(i), this%Dihedral%RY1(i), this%Dihedral%RZ1(i))
        Rj = vector(this%Dihedral%RX2(i), this%Dihedral%RY2(i), this%Dihedral%RZ2(i))
        Rm = vector(this%Dihedral%RX3(i), this%Dihedral%RY3(i), this%Dihedral%RZ3(i))
        Rl = vector(this%Dihedral%RX4(i), this%Dihedral%RY4(i), this%Dihedral%RZ4(i))

!CDIR NODEP

        deri = 0._RK
        if (this%nmax .eq. 0) then
           earg = 1._RK + cos(-this%gamma0(1))
           EPotLocal = EPotLocal + earg * this%ForConst(1)
        else
          ! Calculate vectors IJ, JK, KL
          a = SUB_VECTOR(Rj, Ri)
          b = SUB_VECTOR(Rm, Rj)
          c = SUB_VECTOR(Rl, Rm)

          a = SUB_ANINT_VECTOR(a)
          a = SCALE_VECTOR(a, BoxLength)

          b = SUB_ANINT_VECTOR(b)
          b = SCALE_VECTOR(b, BoxLength)

          c = SUB_ANINT_VECTOR(c)
          c = SCALE_VECTOR(c, BoxLength)

          ab = SCALAR_PRODUCT(a, b)
          bc = SCALAR_PRODUCT(b, c)
          ac = SCALAR_PRODUCT(a, c)
          aa = SCALAR_PRODUCT(a, a)
          bb = SCALAR_PRODUCT(b, b)
          cc = SCALAR_PRODUCT(c, c)

          ! Vector products
          axb = (aa*bb) - (ab*ab)
          bxc = (bb*cc) - (bc*bc)

          num = (ab*bc) - (ac*bb)
          den = axb*bxc

          ! Check, that any 3 atoms don't lie on one line and they define good dihedral angle:
          ! (Otherwise contribution is zero)

          if ( den .gt. 1E-10_RK ) then
            den = sqrt( den )

            ! cos of angle:
            co = num/den
            if ( co .gt. 1._RK ) co =   1._RK
            if ( co .lt. -1._RK ) co = -1._RK

            ! sign of angle:
            signum = a%x*(b%y*c%z-c%y*b%z)+a%y*(b%z*c%x-c%z*b%x)+a%z*(b%x*c%y-c%x*b%y)

            ! Value of angle:
            arg = sign( acos(co), signum)
            si = sin(arg)
            if( abs(si) .lt. 1E-10_RK ) si = sign( 1E-10_RK, si )


            if (this%nmax > 0) then
              ! Normal Amber-type torsion angle
              earg = 1._RK + cos(-this%gamma0(1))
              EPotLocal = EPotLocal + earg * this%ForConst(1)
              do j = 1,this%nmax
                earg= j*arg-this%gamma0(j+1)
                ! Energy and forces:
                ! formulae  E = ForConst*( 1 + cos(earg) )
                !           F = ForConst*n*sin(earg)
                EPotLocal = EPotLocal + this%ForConst(j+1)*(1._RK+cos(earg))
                deri = deri - this%ForConst(j+1)*j*sin(earg)
              end do

             else ! Improper dihedral angle
               earg= arg-this%gamma0(1)

               ! Energy and forces:
               ! formulae  E = ForConst*earg**2
               !           F = -2*ForConst*earg
                EPotLocal = EPotLocal + this%ForConst(1)*earg**2
                deri = 2._RK*this%ForConst(1)*earg
             end if

             ! Calculate Forces
             axb = axb/den*co
             bxc = bxc/den*co
             de1 = deri/den/si

            ! X components
            dnum = c%x*bb - b%x*bc
            dden = ( ab*b%x - a%x*bb )*bxc
            FFI = (dnum - dden) * de1
            dnum = ((b%x-a%x)*bc - ab*c%x ) + (2.0*ac*b%x - c%x*bb)
            dden = axb*(bc*c%x-b%x*cc) + (a%x*bb-aa*b%x-ab*(b%x-a%x))*bxc
            FFJ = (dnum - dden) * de1
            dnum = ab*b%x - a%x*bb
            dden = axb*( bb*c%x - bc*b%x )
            FFL = (dnum - dden) * de1
            FFK = -(ffi+ffj+ffl)

            ! Forces
            FX1(i) = FX1(i)+ffi
            FX2(i) = FX2(i)+ffj
            FX3(i) = FX3(i)+ffk
            FX4(i) = FX4(i)+ffl

            ! Y components
            dnum = c%y*bb - b%y*bc
            dden = ( ab*b%y - a%y*bb )*bxc
            FFI = (dnum - dden) * de1
            dnum = ((b%y-a%y)*bc - ab*c%y ) + (2.0*ac*b%y - c%y*bb)
            dden = axb*(bc*c%y-b%y*cc) + (a%y*bb-aa*b%y-ab*(b%y-a%y))*bxc
            FFJ = (dnum - dden) * de1
            dnum = ab*b%y - a%y*bb
            dden = axb*( bb*c%y - bc*b%y )
            FFL = (dnum - dden) * de1
            FFK = -(ffi+ffj+ffl)

            ! Forces
            FY1(i) = FY1(i)+ffi
            FY2(i) = FY2(i)+ffj
            FY3(i) = FY3(i)+ffk
            FY4(i) = FY4(i)+ffl

            ! Z components
            dnum = c%z*bb - b%z*bc
            dden = ( ab*b%z - a%z*bb )*bxc
            FFI = (dnum - dden) * de1
            dnum = ((b%z-a%z)*bc - ab*c%z ) + (2.0*ac*b%z - c%z*bb)
            dden = axb*(bc*c%z-b%z*cc) + (a%z*bb-aa*b%z-ab*(b%z-a%z))*bxc
            FFJ = (dnum - dden) * de1
            dnum = ab*b%z - a%z*bb
            dden = axb*( bb*c%z - bc*b%z )
            FFL = (dnum - dden) * de1
            FFK = -(ffi+ffj+ffl)

            ! Forces
            FZ1(i) = FZ1(i)+ffi
            FZ2(i) = FZ2(i)+ffj
            FZ3(i) = FZ3(i)+ffk
            FZ4(i) = FZ4(i)+ffl

          endif ! den>0
        endif ! nmax/=0
      enddo

    ! Update potential energy, no contribution to virial!
     EPot = EPot +  EPotLocal

    ! Update Intra potential energy
     EPotIntra_Dihedral = EPotIntra_Dihedral +  EPotLocal

  end subroutine TPotDihedral_Force


!==============================================================!
!  Subroutine TPotDihedral_Energy                              !
!==============================================================!

  subroutine TPotDihedral_Energy( this, np, nu, F, EDihedral,  BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDihedral)       :: this
    integer, intent(in)      :: np
    integer, intent(in)      :: nu
    real(RK), intent(in out) :: F(3,nu)
    real(RK), intent(in out) :: EDihedral
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    integer      :: j
    real(RK)     :: num, deri, den, ax, ay, az, bx, by, bz, cx, cy, cz
    real(RK)     :: ab, bc, ac, aa, bb, cc, axb, bxc, co, signum, arg, earg

    deri = 0._RK
    if (this%nmax .eq. 0) then
        earg = 1._RK + cos(-this%gamma0(1))
        EDihedral = EDihedral + earg * this%ForConst(1)
    else
      ! Calculate vectors IJ, JK, KL
      ax = this%Dihedral%RX2(np) - this%Dihedral%RX1(np)
      ay = this%Dihedral%RY2(np) - this%Dihedral%RY1(np)
      az = this%Dihedral%RZ2(np) - this%Dihedral%RZ1(np)
      bx = this%Dihedral%RX3(np) - this%Dihedral%RX2(np)
      by = this%Dihedral%RY3(np) - this%Dihedral%RY2(np)
      bz = this%Dihedral%RZ3(np) - this%Dihedral%RZ2(np)
      cx = this%Dihedral%RX4(np) - this%Dihedral%RX3(np)
      cy = this%Dihedral%RY4(np) - this%Dihedral%RY3(np)
      cz = this%Dihedral%RZ4(np) - this%Dihedral%RZ3(np)
      !
      ax = (ax - anint( ax )) * BoxLength
      ay = (ay - anint( ay )) * BoxLength
      az = (az - anint( az )) * BoxLength
      bx = (bx - anint( bx )) * BoxLength
      by = (by - anint( by )) * BoxLength
      bz = (bz - anint( bz )) * BoxLength
      cx = (cx - anint( cx )) * BoxLength
      cy = (cy - anint( cy )) * BoxLength
      cz = (cz - anint( cz )) * BoxLength

      ! Scalar products
      ab = ax*bx + ay*by + az*bz
      bc = bx*cx + by*cy + bz*cz
      ac = ax*cx + ay*cy + az*cz
      aa = ax*ax + ay*ay + az*az
      bb = bx*bx + by*by + bz*bz
      cc = cx*cx + cy*cy + cz*cz

      ! Vector products
      axb = (aa*bb) - (ab*ab)
      bxc = (bb*cc) - (bc*bc)

      num = (ab*bc) - (ac*bb)
      den = axb*bxc

      ! Check, that any 3 atoms don't lie on one line and they define good dihedral angle:
      ! (Otherwise contribution is zero)
      if ( den .gt. 1E-10_RK ) then
        den = sqrt( den )
        ! cos of angle:
        co = num/den
        if ( co .gt. 1._RK ) co =   1._RK
        if ( co .lt. -1._RK ) co = -1._RK
        ! sign of angle:
        signum = ax*(by*cz-cy*bz)+ay*(bz*cx-cz*bx)+az*(bx*cy-cx*by)
        ! Value of angle:
        arg = sign( acos(co), signum)

        if (this%nmax > 0) then
          ! Normal Amber-type torsion angle
          earg = 1._RK + cos(-this%gamma0(1))
          EDihedral = EDihedral + earg * this%ForConst(1)
          do j = 1,this%nmax
            earg= j*arg-this%gamma0(j+1)

            EDihedral = EDihedral + this%ForConst(j+1)*(1._RK+cos(earg))
          end do

        else ! Improper dihedral angle
          earg= arg-this%gamma0(1)
          EDihedral = EDihedral + this%ForConst(1)*earg**2
        end if

      endif ! den>0
    endif ! nmax/=0

  end subroutine TPotDihedral_Energy


end module ms2_potential
