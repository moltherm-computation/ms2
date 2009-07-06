!==============================================================!
!  MOLECULAR SIMULATION PROGRAM MS2 Version 1.1 v12            !
!  (c) 2001 by Sergey Lishchuk, ITT                            !
!  (c) 2007 by Bernhard Eckl, ITT                              !
!  (c) 2008 by Ekaterina Elts, TUM                             !
!==============================================================!
!  Module ms2_potential                                        !
!  Contains TPot* objects                                      !
!==============================================================!

!==============================================================!
! ChangeLog                                                    !
!==============================================================!
! 04/07/09  PotQuadrupoleDipole ergaenzt                       !
!           Subroutinen Force, Chempot und Energy angelegt     !
!           mit 2CLJDQ von Juergen verglichen                  !
!                                                              !
! 05/03/30  Charge Potentials ergaenzt (MC Energy noch falsch) !
!                                                              !
! 05/09/15  Thorsten hat alles bis auf Chempo gecheckt         !
!                                                              !
! 05/12/15  Schleifen mit neuer Cutoffberechnung konsistent    !
!                                                              !
!==============================================================!

#ifndef ARCH
#define ARCH    0
#define FORTRAN 90
#define MPI_VER 0
#endif

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2_potential.F90...'
#endif

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
    real(RK)                  :: EPotCorr, VirialCorr, EPotTestCorr
    logical                   :: SameComponent
    logical                   :: potintra15, potintra14
    real(RK)                  :: SigmaSquared
    real(RK)                  :: Epsilon4, Epsilon48
    real(RK)                  :: BoxlengthInv, BoxLengthThird
    real(RK)                  :: ScaleLJ14
    integer, pointer          :: NInCutoff(:), CutoffPartner(:, :)

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
    real(RK)                   :: ScaleEl14
    logical                    :: SameComponent
    logical                    :: potintra15, potintra14
    integer, pointer           :: NInCutoff(:), CutoffPartner(:, :)

  end type TPotChargeCharge

  interface Construct
    module procedure TPotCC_Construct
  end interface

  interface Destruct
    module procedure TPotCC_Destruct
  end interface

  interface Force
    module procedure TPotCC_Force
  end interface

  interface ChemicalPotential
    module procedure TPotCC_ChemicalPotential
  end interface

  interface Energy
    module procedure TPotCC_Energy
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
    integer, pointer           :: NInCutoff(:), CutoffPartner(:, :)

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
    integer, pointer               :: NInCutoff(:), CutoffPartner(:, :)

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
    integer, pointer           :: NInCutoff(:), CutoffPartner(:, :)

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
    integer, pointer           :: NInCutoff(:), CutoffPartner(:, :)

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
    integer, pointer               :: NInCutoff(:), CutoffPartner(:, :)

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
    integer, pointer               :: NInCutoff(:), CutoffPartner(:, :)

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
    integer, pointer               :: NInCutoff(:), CutoffPartner(:, :)

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
    integer, pointer               :: NInCutoff(:), CutoffPartner(:, :)

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
    real(RK)                  :: ForConst, R0
    real(RK)                  :: EPotCorr, VirialCorr, EPotTestCorr
    real(RK)                  :: BoxlengthInv, BoxLengthThird

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

  interface ChemicalPotential
    module procedure TPotBond_ChemicalPotential
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
    real(RK)                  :: ForConst, Angle0
    real(RK)                  :: EPotCorr, VirialCorr, EPotTestCorr
    real(RK)                  :: BoxlengthInv, BoxLengthThird

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

  interface ChemicalPotential
    module procedure TPotAngle_ChemicalPotential
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
    integer                      :: multi
    real(RK)                     :: ForConst, gamma
    real(RK)                     :: ScaleLJ14, ScaleEl14
    real(RK)                     :: EPotCorr, VirialCorr, EPotTestCorr
    real(RK)                     :: BoxlengthInv, BoxLengthThird
    real(RK)                     :: Sigma1, Sigma4, Epsilon1, Epsilon4

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

  interface ChemicalPotential
    module procedure TPotDihedral_ChemicalPotential
  end interface

  interface Energy
    module procedure TPotDihedral_Energy
  end interface

contains

!==============================================================!
!  Subroutine TPotLJLJ_Construct                               !
!==============================================================!

  subroutine TPotLJLJ_Construct( this, i1, i2, j1, j2, &
&                                Molecule1, Molecule2, &
&                                RCutoff, ScaleSigma, ScaleEpsilon )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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
!    if (this%SameComponent .and. IntraLJEL .and. this%Site1%SiteId<this%Site2%SiteId) then
    if (this%SameComponent .and. IntraLJEL ) then
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
      else
        this%EPotCorr = Pi8 * this%Epsilon * &
&         ( TICSu(-6, RCutoff, this%Sigma**2, tau) &
&         - TICSu(-3, RCutoff, this%Sigma**2, tau) )
        this%VirialCorr = Piminus83 * this%Epsilon * &
&         ( TICSp(-6, RCutoff, this%Sigma**2, tau) &
&         - TICSp(-3, RCutoff, this%Sigma**2, tau) )
      endif
    else ! Site-site cutoff or both sites in center of mass
      RCutoff3Inv = (this%Sigma / RCutoff)**3
      RCutoff9Inv = RCutoff3Inv**3
      this%EPotCorr = Pi89 * this%Epsilon &
&       * (RCutoff9Inv - 3._RK * RCutoff3Inv)
      this%VirialCorr = Pi329 * this%Epsilon &
&       * (RCutoff9Inv - 1.5_RK * RCutoff3Inv)
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
&               - (rc-tauMinus)**(2*n+4) + (rc-tauPlus)**(2*n+4) ) * rc &
&             / ( 8._RK * sigma2**n * tau1 * tau2 &
&               * (n+1) * (2*n+3) * (2*n+4) ) &
&             + ( (rc+tauPlus)**(2*n+5) - (rc+tauMinus)**(2*n+5) &
&               - (rc-tauMinus)**(2*n+5) + (rc-tauPlus)**(2*n+5) ) &
&             / ( 8._RK * sigma2**n * tau1 * tau2 &
&               * (n+1) * (2*n+3) * (2*n+4) * (2*n+5) )

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
&               - (rc-tauMinus)**(2*n+3) + (rc-tauPlus)**(2*n+3) ) * rc**2 &
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



  end subroutine TPotLJLJ_Construct



!==============================================================!
!  Subroutine TPotLJLJ_Destruct                                !
!==============================================================!

  subroutine TPotLJLJ_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotLJ126LJ126) :: this

    ! Destroy potential
    continue

  end subroutine TPotLJLJ_Destruct



!==============================================================!
!  Subroutine TPotLJLJ_Force                                   !
!==============================================================!

  subroutine TPotLJLJ_Force( this, EPot, Virial, EPotInter, VirialInter, EPotIntra, VirialIntra, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotLJ126LJ126)     :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
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
    real(RK)          :: EPotLocal, VirialLocal
    real(RK)          :: EPotLocalIntra, VirialLocalIntra
    real(RK)          :: EPotLocalInter, VirialLocalInter
    logical           :: SameComponent, noIntra, choice
    integer           :: i, j, k, i1, j0, j1, m, ende
    integer           :: nu1, nu2, jk, unit
    integer           :: unit1, unit2
    logical           :: ok, intra15, intra14
    real(RK)          :: coeff
#if MPI_VER > 0
    integer           :: i0, N1, N2, ji
    logical           :: EvenN
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
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    SigmaSquared = this%SigmaSquared
    Epsilon4 = this%Epsilon4
    Epsilon48 = this%Epsilon48
    RCutoffSquared = this%RCutoffSquaredScaled
    EPotLocal   = 0._RK
    VirialLocal = 0._RK
    EPotLocalIntra   = 0._RK
    VirialLocalIntra = 0._RK
    EPotLocalInter   = 0._RK
    VirialLocalInter = 0._RK
    intra15 = this%potintra15
    intra14 = this%potintra14

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

    if (intra14) then
       coeff = this%ScaleLJ14  !Scale 1,4 LJ interaction
    else
       coeff = 1._RK
    end if

   if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        unit=nu1*(i-1)+this%Site1%UnitNumber ! Number of unit, to which this site corresponds
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
loop1:  do k = 1, this%NInCutoff(unit)
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
            PXij = PXij - anint( PXij )
            PYij = PYij - anint( PYij )
            PZij = PZij - anint( PZij )
            RijSquaredInv = SigmaSquared / ( RXij**2 + RYij**2 + RZij**2 )
            Rij6Inv = RijSquaredInv**3
            EPotLocal = EPotLocal + Rij6Inv * (Rij6Inv - 1._RK)
            EPotLocalInter = EPotLocalInter + Rij6Inv * (Rij6Inv - 1._RK)
            Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
            FXij = Fij * RXij
            FYij = Fij * RYij
            FZij = Fij * RZij
            VirialLocal = VirialLocal + PXij * FXij + PYij * FYij + PZij * FZij
            VirialLocalInter = VirialLocalInter + PXij * FXij + PYij * FYij + PZij * FZij
            FXi = FXi + FXij
            FYi = FYi + FYij
            FZi = FZi + FZij
            FX2(jk) = FX2(jk) - FXij
            FY2(jk) = FY2(jk) - FYij
            FZ2(jk) = FZ2(jk) - FZij
          end if
        end do loop1
		! Include intramolecular interaction if need
        if (SameComponent .and. (intra15 .or. intra14)) then
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
            EPotLocal = EPotLocal + Rij6Inv * (Rij6Inv - 1._RK)
            EPotLocalIntra = EPotLocalIntra + Rij6Inv * (Rij6Inv - 1._RK)
            Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
            FXij = Fij * RXij * coeff
            FYij = Fij * RYij * coeff
            FZij = Fij * RZij * coeff
            VirialLocal = VirialLocal + PXij * FXij + PYij * FYij + PZij * FZij
            VirialLocalIntra = VirialLocalIntra + PXij * FXij + PYij * FYij + PZij * FZij
!            print*, 'VirialLocalIntra=', VirialLocalIntra
            FXi = FXi + FXij
            FYi = FYi + FYij
            FZi = FZi + FZij
            FX2(i) = FX2(i) - FXij
            FY2(i) = FY2(i) - FYij
            FZ2(i) = FZ2(i) - FZij
         end if
         FX1(i) = FXi
         FY1(i) = FYi
         FZ1(i) = FZi
      end do
	! Include all intramoleculare interactions if need



        else ! Site-site cutoff ! Should be corrected

!    noIntra = .not. IntraLJEl       ! no Intramolecular interaction
!    choice = SameComponent .and. noIntra ! SameComponent, but no Intramolecular interaction

    ! Loop over molecules
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, merge( i1 - 1, i1, choice )
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
loop2:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )
#else
          j0 = merge( i + 1, 1, SameComponent )

  loop2:  do j = j0, j1
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
            if( RijSquared >= RCutoffSquared ) cycle loop2
            RijSquaredInv = SigmaSquared / RijSquared
            Rij6Inv = RijSquaredInv**3
            EPotLocal = EPotLocal + Rij6Inv * (Rij6Inv - 1._RK)
            EPotLocalInter = EPotLocalInter + Rij6Inv * (Rij6Inv - 1._RK)
            Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
            FXij = Fij * RXij
            FYij = Fij * RYij
            FZij = Fij * RZij
            VirialLocal = VirialLocal + PXij * FXij + PYij * FYij + PZij * FZij
            VirialLocalInter = VirialLocalInter + PXij * FXij + PYij * FYij + PZij * FZij
            FXi = FXi + FXij
            FYi = FYi + FYij
            FZi = FZi + FZij
            FX2(j) = FX2(j) - FXij
            FY2(j) = FY2(j) - FYij
            FZ2(j) = FZ2(j) - FZij
          end do loop2
          FX1(i) = FXi
          FY1(i) = FYi
          FZ1(i) = FZi
        end do

    end if ! Cutoff=COM or site-site

    ! Update potential energy and virial
    EPot = EPot + Epsilon4 * EPotLocal
    Virial = Virial + Third * VirialLocal * BoxLength

    ! Update Inter potential energy and virial
    EPotInter = EPotInter + Epsilon4 * EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter * BoxLength


    if (IntraLJEl) then
      ! Update Intra potential energy and virial
      EPotIntra = EPotIntra + Epsilon4 * EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra * BoxLength
    end if

  end subroutine TPotLJLJ_Force

!==============================================================!
!  Subroutine TPotLJLJ_ChemicalPotential                       !
!==============================================================!

  subroutine TPotLJLJ_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotLJ126LJ126) :: this
    real(RK), pointer    :: EPotTest(:)
    real(RK), intent(in) :: BoxLength

    ! Declare local variables
    real(RK)          :: SigmaSquared
    real(RK)          :: Epsilon4
    real(RK)          :: RCutoffSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: RijSquared, RijSquaredInv, Rij6Inv
    real(RK)          :: EPotLocal
    integer           :: N2
    integer           :: i, j, k

    ! Assign local variables
    N2 = this%Site2%NPart
    SigmaSquared = this%SigmaSquared
    Epsilon4 = this%Epsilon4
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

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
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
          RijSquared = RXij**2 + RYij**2 + RZij**2
          RijSquaredInv = SigmaSquared / RijSquared
          Rij6Inv = RijSquaredInv**3
          EPotLocal = EPotLocal + Rij6Inv * (Rij6Inv - 1._RK)
        end do loop1
        EPotTest(i) = EPotTest(i) + Epsilon4 * EPotLocal
      end do

    else

      ! Loop over test particles
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
          RijSquared = RXij**2 + RYij**2 + RZij**2
          if( RijSquared >= RCutoffSquared ) cycle loop2
          RijSquaredInv = SigmaSquared / RijSquared
          Rij6Inv = RijSquaredInv**3
          EPotLocal = EPotLocal + Rij6Inv * (Rij6Inv - 1._RK)
        end do loop2
        EPotTest(i) = EPotTest(i) + Epsilon4 * EPotLocal
      end do
    end if

  end subroutine TPotLJLJ_ChemicalPotential



!==============================================================!
!  Subroutine TPotLJLJ_Energy                                  !
!==============================================================!

  subroutine TPotLJLJ_Energy( this, np, EPot, Virial, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotLJ126LJ126) :: this
    integer, intent(in)  :: np
    real(RK), pointer    :: EPot(:)
    real(RK), pointer    :: Virial(:)
    real(RK), intent(in) :: BoxLength

    ! Declare local variables
    real(RK)          :: SigmaSquared
    real(RK)          :: Epsilon4, Epsilon48
    real(RK)          :: RCutoffSquared
    real(RK)          :: BoxLengthThird
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij, Fij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: RijSquared, RijSquaredInv, Rij6Inv
    integer           :: N
    integer           :: j, k

    ! Assign local variables
    N = this%Site2%NPart
    SigmaSquared = this%SigmaSquared
    Epsilon4 = this%Epsilon4
    Epsilon48 = this%Epsilon48
    RCutoffSquared = this%RCutoffSquaredScaled
    BoxLengthThird = this%BoxLengthThird

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

    if( CutoffMode .eq. CenterofMass ) then

!CDIR NODEP
loop1:do k = 1, this%NInCutoff(np)
        j = this%CutoffPartner(k, np)
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
        RijSquared = RXij**2 + RYij**2 + RZij**2
        RijSquaredInv = SigmaSquared / RijSquared
        Rij6Inv = RijSquaredInv**3
        EPot(j) = EPot(j) + Epsilon4 * Rij6Inv * (Rij6Inv - 1._RK)
        Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
        FXij = Fij * RXij
        FYij = Fij * RYij
        FZij = Fij * RZij
        Virial(j) = Virial(j) &
&         + BoxLengthThird * (PXij * FXij + PYij * FYij + PZij * FZij)
      end do loop1

    else ! Site-site cutoff

#if MPI_VER > 0
!CDIR NODEP
loop2:do j = this%Site2%NPart0, this%Site2%NPart2
#else
!CDIR NODEP
loop2:do j = 1, N
#endif
        if( this%SameComponent .and. j == np ) cycle loop2
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
        if( RijSquared >= RCutoffSquared ) cycle loop2
        RijSquaredInv = SigmaSquared / RijSquared
        Rij6Inv = RijSquaredInv**3
        EPot(j) = EPot(j) + Epsilon4 * Rij6Inv * (Rij6Inv - 1._RK)
        Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
        FXij = Fij * RXij
        FYij = Fij * RYij
        FZij = Fij * RZij
        Virial(j) = Virial(j) &
&         + BoxLengthThird * (PXij * FXij + PYij * FYij + PZij * FZij)
      end do loop2
    end if

  end subroutine TPotLJLJ_Energy



!==============================================================!
!  Subroutine TPotLJLJ_UpdateBoxLength                         !
!==============================================================!

  subroutine TPotLJLJ_UpdateBoxLength( this, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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
    this%Epsilon48 = 12._RK * this%Epsilon4 * BoxLengthInv &
&      / this%SigmaSquared
    this%RCutoffSquaredScaled = this%RCutoffSquared * BoxLengthInv**2

  end subroutine TPotLJLJ_UpdateBoxLength



!==============================================================!
!  Subroutine TPotCC_Construct                                 !
!==============================================================!

  subroutine TPotCC_Construct( this, i1, i2, j1, j2, &
&                              Molecule1, Molecule2, RCutoff )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotChargeCharge)      :: this
    integer, intent(in)         :: i1, i2, j1, j2
    type(TMolecule), intent(in) :: Molecule1, Molecule2
    real(RK), intent(in)        :: RCutoff

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

    ! if this potential is intra

!    if (this%SameComponent .and. IntraLJEL .and. this%Site1%SiteId<this%Site2%SiteId) then
    if (this%SameComponent .and. IntraLJEL) then
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
    end if

  end subroutine TPotCC_Construct



!==============================================================!
!  Subroutine TPotCC_Destruct                                  !
!==============================================================!

  subroutine TPotCC_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotChargeCharge) :: this

    ! Destroy potential
    continue

  end subroutine TPotCC_Destruct



!==============================================================!
!  Subroutine TPotCC_Force                                     !
!==============================================================!

  subroutine TPotCC_Force( this, EPot, Virial, EPotInter, VirialInter,EPotIntra, VirialIntra, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotChargeCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ                                             ! Site-Site-Einheitvektor
    real(RK)          :: RijInv
    real(RK)          :: EPotLocal, EPotLocal1, VirialLocal
    real(RK)          :: EPotLocalIntra, EPotLocal1Intra, VirialLocalIntra
    real(RK)          :: EPotLocalInter, EPotLocal1Inter, VirialLocalInter
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, jk, unit
    logical           :: intra14, intra15, SameComponent
    real(RK)          :: coeff
#if MPI_VER > 0
    integer           :: i0
#endif

    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    Epsilon = this%Epsilon
    EPotLocal = 0._RK
    VirialLocal = 0._RK
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK
    intra15 = this%potintra15
    intra14 = this%potintra14
    SameComponent = this%SameComponent


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

    if (intra14) then
    		coeff = this%ScaleEl14 ! Scale 1,4 El interaction
	else
		coeff = 1._RK
	end if

    ! Loop over molecules
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
        unit=nu1*(i-1)+this%Site1%UnitNumber
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
#if ARCH == 3
            RijInv = rsqrt( RXij**2 + RYij**2 + RZij**2 )
#else
            RijInv = 1._RK / sqrt( RXij**2 + RYij**2 + RZij**2 )
#endif
            eX = RXij * RijInv
            eY = RYij * RijInv
            eZ = RZij * RijInv
            EPotLocal1 = Epsilon * RijInv
            EPotLocal  = EPotLocal + EPotLocal1
            VirialLocal = VirialLocal + EPotLocal1 * RijInv &
&                           * (eX * PXij + eY * PYij + eZ * PZij)
!            if (intra) then
!	       EPotLocal1Intra = Epsilon * RijInv
!               EPotLocalIntra  = EPotLocalIntra + EPotLocal1Intra
!               VirialLocalIntra = VirialLocalIntra + EPotLocal1Intra * RijInv &
!&                           * (eX * PXij + eY * PYij + eZ * PZij)
!	    else
	       EPotLocal1Inter = Epsilon * RijInv
               EPotLocalInter  = EPotLocalInter + EPotLocal1Inter
               VirialLocalInter = VirialLocalInter + EPotLocal1Inter * RijInv &
&                           * (eX * PXij + eY * PYij + eZ * PZij)

!            end if
            FXij = EPotLocal1 * RijInv * eX * coeff
            FYij = EPotLocal1 * RijInv * eY * coeff
            FZij = EPotLocal1 * RijInv * eZ * coeff
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            FX2(jk) = FX2(jk) - FXij
            FY2(jk) = FY2(jk) - FYij
            FZ2(jk) = FZ2(jk) - FZij
         end if
      end do loop1
      ! Include intramolecular interaction if need
      if (SameComponent .and. (intra15 .or. intra14)) then
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
            EPotLocal1 = Epsilon * RijInv
            EPotLocal  = EPotLocal + EPotLocal1
            VirialLocal = VirialLocal + EPotLocal1 * RijInv &
&                           * (eX * PXij + eY * PYij + eZ * PZij)
            EPotLocal1Intra = Epsilon * RijInv
            EPotLocalIntra  = EPotLocalIntra + EPotLocal1Intra
            VirialLocalIntra = VirialLocalIntra + EPotLocal1Intra * RijInv &
&                           * (eX * PXij + eY * PYij + eZ * PZij)
            FXij = EPotLocal1 * RijInv * eX * coeff
            FYij = EPotLocal1 * RijInv * eY * coeff
            FZij = EPotLocal1 * RijInv * eZ * coeff
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            FX2(i) = FX2(i) - FXij
            FY2(i) = FY2(i) - FYij
            FZ2(i) = FZ2(i) - FZij
         end if
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
    end do

    ! Update potential energy and virial



    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal

    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter

    if (IntraLJEl) then
       EPotIntra = EPotIntra + EPotLocalIntra
       VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if

  end subroutine TPotCC_Force



!==============================================================!
!  Subroutine TPotCC_ChemicalPotential                         !
!==============================================================!

  subroutine TPotCC_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotChargeCharge) :: this
    real(RK), pointer      :: EPotTest(:)
    real(RK), intent(in)   :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: RijInv, RijSquared
    real(RK)          :: EPotLocal
    integer           :: i, j, k, i1
#if ARCH == 3
    logical           :: hit
#endif

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
#if ARCH == 3
        hit = .false.
#endif
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

  subroutine TPotCC_Energy( this, np, EPot, Virial, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotChargeCharge) :: this
    integer, intent(in)    :: np
    real(RK), pointer      :: EPot(:)
    real(RK), pointer      :: Virial(:)
    real(RK), intent(in)   :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijInv, RijSquared
    real(RK)          :: EPotLocal, VirialLocal
    integer           :: j, k

    ! Assign local variables
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

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

!CDIR NODEP
    do k = 1, this%NInCutoff(np)
      j = this%CutoffPartner(k, np)
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
      RijSquared = RXij**2 + RYij**2 + RZij**2
      if( RijSquared <= RShieldSquared ) then
        EPotLocal = 1E33_RK
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
        VirialLocal = EPotLocal * RijInv * (eX * PXij + eY * PYij + eZ * PZij)
      end if
      EPot(j) = EPot(j) + EPotLocal
      Virial(j) = Virial(j) + Third * VirialLocal
    end do

  end subroutine TPotCC_Energy



!==============================================================!
!  Subroutine TPotCD_Construct                                 !
!==============================================================!

  subroutine TPotCD_Construct( this, i1, i2, j1, j2, &
&                              Molecule1, Molecule2, RCutoff )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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

!    if (this%SameComponent .and. IntraLJEL .and. this%Site1%SiteId<this%Site2%SiteId) then
    if (this%SameComponent .and. IntraLJEL) then
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
    end if


  end subroutine TPotCD_Construct


!==============================================================!
!  Subroutine TPotCD_Destruct                                  !
!==============================================================!

  subroutine TPotCD_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotChargeDipole) :: this

    ! Destroy potential
    continue

  end subroutine TPotCD_Destruct



!==============================================================!
!  Subroutine TPotCD_Force                                     !
!==============================================================!

  subroutine TPotCD_Force( this, EPot, Virial, EPotInter, VirialInter,EPotIntra, VirialIntra, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotChargeDipole)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1, Epsilon2
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer :: OX2(:), OY2(:), OZ2(:)
    real(RK), pointer :: TX2(:), TY2(:), TZ2(:)
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
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, jk, unit
    logical           :: intra14, intra15, SameComponent
    real(RK)          :: coeff
#if MPI_VER > 0
    integer           :: i0
#endif

    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    Epsilon = this%Epsilon
    EPotLocal = 0._RK
    VirialLocal = 0._RK
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

    intra14 = this%potintra14
    intra15 = this%potintra15
    SameComponent = this%SameComponent

    if (intra14) then
      coeff = this%ScaleEl14
    else
      coeff = 1._Rk
    end if

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
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ

    ! Loop over molecules
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
        unit=nu1*(i-1)+this%Site1%UnitNumber
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
loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
!          scale = this%ScaleCoeff(k, unit)
!          intra = this%Intra(k,unit)
!          if (scale) then
!             coeff = ScaleEl14
!          else
             coeff = 1._RK
!          end if
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
!            if (intra) then
!              EPotLocalIntra  = EPotLocalIntra + Epsilon1 * CosTheta
!            else
              EPotLocalInter  = EPotLocalInter + Epsilon1 * CosTheta
!            endif
            FXij = Epsilon2 * ( CosTheta3 * eX - OXj ) * coeff                      ! F2 bei Price
            FYij = Epsilon2 * ( CosTheta3 * eY - OYj ) * coeff
            FZij = Epsilon2 * ( CosTheta3 * eZ - OZj ) * coeff
            VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij     ! F2*R_COM_Price; stimmt so
            VirialLocalInter = VirialLocalInter + FXij * PXij + FYij * PYij + FZij * PZij     ! F2*R_COM_Price; stimmt so
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            FX2(jk) = FX2(jk) - FXij
            FY2(jk) = FY2(jk) - FYij
            FZ2(jk) = FZ2(jk) - FZij
            TX2(jk) = TX2(jk) - Epsilon1 * eX                                         ! Uebereinstimmung mit Price
            TY2(jk) = TY2(jk) - Epsilon1 * eY
            TZ2(jk) = TZ2(jk) - Epsilon1 * eZ
          end if
        end do loop1
        ! Include intramolecular interactions if need
        if (SameComponent .and. (intra15 .or. intra14)) then
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
            Epsilon1 = Epsilon * RijSquaredInv
            Epsilon2 = Epsilon1 * RijInv
            EPotLocal  = EPotLocal + Epsilon1 * CosTheta                            ! Uebereinstimmumg mit Price
            EPotLocalIntra  = EPotLocalIntra + Epsilon1 * CosTheta
            FXij = Epsilon2 * ( CosTheta3 * eX - OXj ) * coeff                      ! F2 bei Price
            FYij = Epsilon2 * ( CosTheta3 * eY - OYj ) * coeff
            FZij = Epsilon2 * ( CosTheta3 * eZ - OZj ) * coeff
            VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij     ! F2*R_COM_Price; stimmt so
            VirialLocalIntra = VirialLocalIntra + FXij * PXij + FYij * PYij + FZij * PZij     ! F2*R_COM_Price; stimmt so
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            FX2(i) = FX2(i) - FXij
            FY2(i) = FY2(i) - FYij
            FZ2(i) = FZ2(i) - FZij
            TX2(i) = TX2(i) - Epsilon1 * eX                                         ! Uebereinstimmung mit Price
            TY2(i) = TY2(i) - Epsilon1 * eY
            TZ2(i) = TZ2(i) - Epsilon1 * eZ
        end if
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
      end do

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal

    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter

    if (IntraLJEl) then
        EPotIntra = EPotIntra + EPotLocalIntra
        VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if


  end subroutine TPotCD_Force


!==============================================================!
!  Subroutine TPotCD_ChemicalPotential                         !
!==============================================================!

  subroutine TPotCD_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotChargeDipole) :: this
    real(RK), pointer      :: EPotTest(:)
    real(RK), intent(in)   :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer :: OX2(:), OY2(:), OZ2(:)
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
#if ARCH == 3
    logical           :: hit
#endif

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
#if ARCH == 3
     hit = .false.
#endif
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

  subroutine TPotCD_Energy( this, np, EPot, Virial, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotChargeDipole) :: this
    integer, intent(in)    :: np
    real(RK), pointer      :: EPot(:)
    real(RK), pointer      :: Virial(:)
    real(RK), intent(in)   :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1, Epsilon2
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer :: OX2(:), OY2(:), OZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv, RijSquared
    real(RK)          :: EPotLocal, VirialLocal
    real(RK)          :: CosTheta, CosTheta3
    integer           :: j, k

    ! Assign local variables
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

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

!CDIR NODEP
    do k = 1, this%NInCutoff(np)
      j = this%CutoffPartner(k, np)
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
        CosTheta3 = 3._RK * CosTheta
        Epsilon1 = Epsilon * RijSquaredInv
        Epsilon2 = Epsilon1 * RijInv
        EPotLocal  = Epsilon1 * CosTheta
        VirialLocal =  Epsilon2 * ( ( CosTheta3 * eX - OXj ) * PXij &
&                                 + ( CosTheta3 * eY - OYj ) * PYij &
&                                 + ( CosTheta3 * eZ - OZj ) * PZij )
      end if
      EPot(j) = EPot(j) + EPotLocal
      Virial(j) = Virial(j) + Third * VirialLocal
    end do

  end subroutine TPotCD_Energy


!==============================================================!
!  Subroutine TPotCQ_Construct                                 !
!==============================================================!

  subroutine TPotCQ_Construct( this, i1, i2, j1, j2, &
&                              Molecule1, Molecule2, RCutoff )
    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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
!    if (this%SameComponent .and. IntraLJEL .and. this%Site1%SiteId<this%Site2%SiteId) then
    if (this%SameComponent .and. IntraLJEL) then
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
    end if


  end subroutine TPotCQ_Construct


!==============================================================!
!  Subroutine TPotCQ_Destruct                                  !
!==============================================================!

  subroutine TPotCQ_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotChargeQuadrupole) :: this

    ! Destroy potential
    continue

  end subroutine TPotCQ_Destruct


!==============================================================!
!  Subroutine TPotCQ_Force                                     !
!==============================================================!

  subroutine TPotCQ_Force( this, EPot, Virial, EPotInter, VirialInter,EPotIntra, VirialIntra, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotChargeQuadrupole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1, Epsilon2
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer :: OX2(:), OY2(:), OZ2(:)
    real(RK), pointer :: TX2(:), TY2(:), TZ2(:)
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
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, jk, unit
    logical           :: intra14, intra15, SameComponent
    real(RK)          :: coeff
#if MPI_VER > 0
    integer           :: i0
#endif

    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    Epsilon = this%Epsilon
    EPotLocal = 0._RK
    VirialLocal = 0._RK
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

    intra14 = this%potintra14
    intra15 = this%potintra15
    SameComponent = this%SameComponent


    if (intra14) then
      coeff = this%ScaleEl14 !Scale 1,4 El interactions
    else
      coeff = 1._RK
    end if

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
    OX2 => this%Site2%OX
    OY2 => this%Site2%OY
    OZ2 => this%Site2%OZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ

    ! Loop over molecules
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
        unit=nu1*(i-1)+this%Site1%UnitNumber
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
loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
!          intra = this%Intra(k, unit)
!          scale = this%ScaleCoeff(k, unit)
!          if (scale) then
!             coeff = ScaleEl14
!          else
             coeff = 1._RK
!          end if
          if ( mod(j-this%Site2%UnitNumber, nu2)==0) then
!            print *, 'j =', j
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
!            if (intra) then
!              EPotLocalIntra  = EPotLocalIntra + Epsilon1 * ( CosTheta * CosTheta - Third )
!	    else
              EPotLocalInter  = EPotLocalInter + Epsilon1 * ( CosTheta * CosTheta - Third )
!	    end if
            CosTheta2 = 2._RK * CosTheta
            CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
            Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
            FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXj ) * coeff                     ! F2 nach Price bzw. Kraft auf Punktladung
            FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYj ) * coeff
            FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZj ) * coeff
            VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij     ! Vorzeichen richtig so
            VirialLocalInter = VirialLocalInter + FXij * PXij + FYij * PYij + FZij * PZij     ! Vorzeichen richtig so
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            FX2(jk) = FX2(jk) - FXij
            FY2(jk) = FY2(jk) - FYij
            FZ2(jk) = FZ2(jk) - FZij
            TX2(jk) = TX2(jk) - Epsilon1 * CosTheta2 * eX
            TY2(jk) = TY2(jk) - Epsilon1 * CosTheta2 * eY
            TZ2(jk) = TZ2(jk) - Epsilon1 * CosTheta2 * eZ
          end if
        end do loop1

        ! Include intramolecular interactions if need
        if (SameComponent .and. (intra14 .or. intra15)) then
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
            Epsilon1 = Epsilon * RijSquaredInv * RijInv
            EPotLocal  = EPotLocal + Epsilon1 * ( CosTheta * CosTheta - Third )
            EPotLocalIntra  = EPotLocalIntra + Epsilon1 * ( CosTheta * CosTheta - Third )
            CosTheta2 = 2._RK * CosTheta
            CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
            Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
            FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXj ) * coeff                     ! F2 nach Price bzw. Kraft auf Punktladung
            FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYj ) * coeff
            FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZj ) * coeff
            VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij     ! Vorzeichen richtig so
            VirialLocalIntra = VirialLocalIntra + FXij * PXij + FYij * PYij + FZij * PZij     ! Vorzeichen richtig so
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            FX2(i) = FX2(i) - FXij
            FY2(i) = FY2(i) - FYij
            FZ2(i) = FZ2(i) - FZij
            TX2(i) = TX2(i) - Epsilon1 * CosTheta2 * eX
            TY2(i) = TY2(i) - Epsilon1 * CosTheta2 * eY
            TZ2(i) = TZ2(i) - Epsilon1 * CosTheta2 * eZ
          end if
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
      end do

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal

    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter

    if (IntraLJEl) then
      EPotIntra = EPotIntra + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if

  end subroutine TPotCQ_Force



!==============================================================!
!  Subroutine TPotCQ_ChemicalPotential                         !
!==============================================================!

  subroutine TPotCQ_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotChargeQuadrupole) :: this
    real(RK), pointer          :: EPotTest(:)
    real(RK), intent(in)       :: BoxLength


    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer :: OX2(:), OY2(:), OZ2(:)
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
#if ARCH == 3
    logical           :: hit
#endif

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
#if ARCH == 3
        hit = .false.
#endif
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
          EPotLocal  = EPotLocal + Epsilon * RijSquaredInv * RijInv &
&                        * ( CosTheta * CosTheta - Third )
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

  subroutine TPotCQ_Energy( this, np, EPot, Virial, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotChargeQuadrupole) :: this
    integer, intent(in)        :: np
    real(RK), pointer          :: EPot(:)
    real(RK), pointer          :: Virial(:)
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon2
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer :: OX2(:), OY2(:), OZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv, RijSquared
    real(RK)          :: CosTheta, CosTheta2, CosAux
    real(RK)          :: EPotLocal, VirialLocal
    integer           :: j, k

    ! Assign local variables
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

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

!CDIR NODEP
    do k = 1, this%NInCutoff(np)
      j = this%CutoffPartner(k, np)
      RXij = RXi - RX2(j)
      RYij = RYi - RY2(j)
      RZij = RZi - RZ2(j)
      PXij = PXi - PX2(j)
      PYij = PYi - PY2(j)
      PZij = PZi - PZ2(j)
      RXij = (RXij - anint( PXij )) * BoxLength                                 ! Abstandsvektor von Q nach C wie bei Price
      RYij = (RYij - anint( PYij )) * BoxLength
      RZij = (RZij - anint( PZij )) * BoxLength
      PXij = (PXij - anint( PXij )) * BoxLength
      PYij = (PYij - anint( PYij )) * BoxLength
      PZij = (PZij - anint( PZij )) * BoxLength                                 ! Orientierungsvektor Quadrupol
      OXj = OX2(j)
      OYj = OY2(j)
      OZj = OZ2(j)
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
        EPotLocal  = Epsilon * RijSquaredInv * RijInv &
&                      * ( CosTheta * CosTheta - Third )
        CosTheta2 = 2._RK * CosTheta
        CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
        Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
        VirialLocal =  Epsilon2 * ( ( CosAux * eX - CosTheta2 * OXj ) * PXij &
&                                 + ( CosAux * eY - CosTheta2 * OYj ) * PYij &
&                                 + ( CosAux * eZ - CosTheta2 * OZj ) * PZij )
      end if
      EPot(j) = EPot(j) + EPotLocal
      Virial(j) = Virial(j) + Third * VirialLocal
    end do

  end subroutine TPotCQ_Energy


!==============================================================!
!  Subroutine TPotDC_Construct                                 !
!==============================================================!

  subroutine TPotDC_Construct( this, i1, i2, j1, j2, &
&                              Molecule1, Molecule2, RCutoff )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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
!    if (this%SameComponent .and. IntraLJEL .and. this%Site1%SiteId<this%Site2%SiteId) then
    if (this%SameComponent .and. IntraLJEL) then
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
    end if
  end subroutine TPotDC_Construct


!==============================================================!
!  Subroutine TPotDC_Destruct                                  !
!==============================================================!

  subroutine TPotDC_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDipoleCharge) :: this

    ! Destroy potential
    continue

  end subroutine TPotDC_Destruct



!==============================================================!
!  Subroutine TPotDC_Force                                     !
!==============================================================!

  subroutine TPotDC_Force( this, EPot, Virial, EPotInter, VirialInter,EPotIntra, VirialIntra, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDipoleCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1, Epsilon2
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:)
    real(RK), pointer :: TX1(:), TY1(:), TZ1(:)
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
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, jk, unit
    logical           :: intra14, intra15, SameComponent
    real(RK)          :: coeff
#if MPI_VER > 0
    integer           :: i0
#endif

    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    Epsilon = this%Epsilon
    EPotLocal = 0._RK
    VirialLocal = 0._RK
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

    intra14 = this%potintra14
    intra15 = this%potintra15
    SameComponent = this%SameComponent
    if (intra14) then
      coeff = this%ScaleEl14 !Scale 1,4 El interactions
    else
      coeff = 1._RK
    end if


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
    OX1 => this%Site1%OX
    OY1 => this%Site1%OY
    OZ1 => this%Site1%OZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ

    ! Loop over molecules
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
        unit=nu1*(i-1)+this%Site1%UnitNumber
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
loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
!          intra = this%Intra(k, unit)
!          scale = this%ScaleCoeff(k, unit)
!          if (scale) then
!             coeff = ScaleEl14
!          else
             coeff = 1._RK
!          end if
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
            EPotLocal  = EPotLocal - Epsilon1 * CosTheta                            ! Uebereinstimmumg mit Price
!            if (intra) then
!              EPotLocalIntra  = EPotLocalIntra - Epsilon1 * CosTheta                ! Uebereinstimmumg mit Price
!            else
              EPotLocalInter  = EPotLocalInter - Epsilon1 * CosTheta                ! Uebereinstimmumg mit Price
!            end if
            FXij = Epsilon2 * ( OXi - CosTheta3 * eX ) * coeff                      ! F1 bei Price
            FYij = Epsilon2 * ( OYi - CosTheta3 * eY ) * coeff
            FZij = Epsilon2 * ( OZi - CosTheta3 * eZ ) * coeff
            VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij     ! F1*(-R_COM_Price); stimmt so
            VirialLocalInter = VirialLocalInter + FXij * PXij + FYij * PYij + FZij * PZij
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            TXi = TXi + Epsilon1 * eX                                               ! Uebereinstimmumg mit Price; Rest bei Atom2Mol in Component
            TYi = TYi + Epsilon1 * eY                                               ! Reaktionsfeldbeitrag in Interaction
            TZi = TZi + Epsilon1 * eZ
            FX2(jk) = FX2(jk) - FXij
            FY2(jk) = FY2(jk) - FYij
            FZ2(jk) = FZ2(jk) - FZij
          end if
        end do loop1
        ! Include intramolecular interaction if need
        if (SameComponent .and. (intra15 .or. intra14)) then
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
            Epsilon1 = Epsilon * RijSquaredInv
            Epsilon2 = Epsilon1 * RijInv
            EPotLocal  = EPotLocal - Epsilon1 * CosTheta                            ! Uebereinstimmumg mit Price
            EPotLocalIntra  = EPotLocalIntra - Epsilon1 * CosTheta                ! Uebereinstimmumg mit Price
            FXij = Epsilon2 * ( OXi - CosTheta3 * eX ) * coeff                      ! F1 bei Price
            FYij = Epsilon2 * ( OYi - CosTheta3 * eY ) * coeff
            FZij = Epsilon2 * ( OZi - CosTheta3 * eZ ) * coeff
            VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij     ! F1*(-R_COM_Price); stimmt so
            VirialLocalIntra = VirialLocalIntra + FXij * PXij + FYij * PYij + FZij * PZij
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            TXi = TXi + Epsilon1 * eX                                               ! Uebereinstimmumg mit Price; Rest bei Atom2Mol in Component
            TYi = TYi + Epsilon1 * eY                                               ! Reaktionsfeldbeitrag in Interaction
            TZi = TZi + Epsilon1 * eZ
            FX2(i) = FX2(i) - FXij
            FY2(i) = FY2(i) - FYij
            FZ2(i) = FZ2(i) - FZij
          end if
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
    end do

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal

    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter

    if (IntraLJEl) then
      EPotIntra = EPotIntra + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if

  end subroutine TPotDC_Force



!==============================================================!
!  Subroutine TPotDC_ChemicalPotential                         !
!==============================================================!

  subroutine TPotDC_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDipoleCharge) :: this
    real(RK), pointer      :: EPotTest(:)
    real(RK), intent(in)   :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:)
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
#if ARCH == 3
    logical           :: hit
#endif

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
#if ARCH == 3
        hit = .false.
#endif
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

  subroutine TPotDC_Energy( this, np, EPot, Virial, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDipoleCharge) :: this
    integer, intent(in)    :: np
    real(RK), pointer      :: EPot(:)
    real(RK), pointer      :: Virial(:)
    real(RK), intent(in)   :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1, Epsilon2
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv, RijSquared
    real(RK)          :: CosTheta, CosTheta3
    real(RK)          :: EPotLocal, VirialLocal
    integer           :: j, k

    ! Assign local variables
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

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

!CDIR NODEP
    do k = 1, this%NInCutoff(np)
      j = this%CutoffPartner(k, np)
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
        CosTheta3 = 3._RK * CosTheta
        Epsilon1 = Epsilon * RijSquaredInv
        Epsilon2 = Epsilon1 * RijInv
        EPotLocal = - Epsilon1 * CosTheta
        VirialLocal = Epsilon2 * ( ( OXi - CosTheta3 * eX ) * PXij &
&                                + ( OYi - CosTheta3 * eY ) * PYij &
&                                + ( OZi - CosTheta3 * eZ ) * PZij )
      end if
      EPot(j) = EPot(j) + EPotLocal
      Virial(j) = Virial(j) + Third * VirialLocal
    end do

  end subroutine TPotDC_Energy


!==============================================================!
!  Subroutine TPotDD_Construct                                 !
!==============================================================!

  subroutine TPotDD_Construct( this, i1, i2, j1, j2, &
&                              Molecule1, Molecule2, RCutoff, RFEpsilon )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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
    this%RFConstant = this%Epsilon / RCutoff**3 &
&     * (RFEpsilon - 1._RK) / (2._RK * RFEpsilon + 1._RK)

    ! if this potential is intra
!    if (this%SameComponent .and. IntraLJEL .and. this%Site1%SiteId<this%Site2%SiteId) then
    if (this%SameComponent .and. IntraLJEL) then
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
    end if


  end subroutine TPotDD_Construct



!==============================================================!
!  Subroutine TPotDD_Destruct                                  !
!==============================================================!

  subroutine TPotDD_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDipoleDipole) :: this

    ! Destroy potential
    continue

  end subroutine TPotDD_Destruct



!==============================================================!
!  Subroutine TPotDD_Force                                     !
!==============================================================!

  subroutine TPotDD_Force( this, EPot, Virial, EPotInter, VirialInter,EPotIntra, VirialIntra, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDipoleDipole)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RFConstant2
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer :: TX1(:), TY1(:), TZ1(:), TX2(:), TY2(:), TZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
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
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
    integer           :: nu1, nu2, jk, unit
    logical           :: intra14, intra15
    real(RK)          :: coeff
#if MPI_VER > 0
    integer           :: N1, N2, i0, ji
    logical           :: EvenN
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
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RFConstant2 = 2._RK * this%RFConstant
    EPotLocal = 0._RK
    VirialLocal = 0._RK
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

    intra14 = this%potintra14
    intra15 = this%potintra15
    if (intra14) then
      coeff = this%ScaleEl14 !Scale 1,4 El interactions
    else
      coeff = 1._RK
    end if


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
    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        unit=nu1*(i-1)+this%Site1%UnitNumber
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
loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
!          intra = this%Intra(k, unit)
!          scale = this%ScaleCoeff(k, unit)
!          if (scale) then
!             coeff = ScaleEl14
!          else
             coeff = 1._RK
!          end if
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
!            if (intra) then
!              EPotLocalIntra = EPotLocalIntra +  Rij3Inv * Tmp
!            else
              EPotLocalInter = EPotLocalInter +  Rij3Inv * Tmp
!            end if
!           EPotLocal = EPotLocal +  Rij3Inv * Tmp - RFConstant2 * CosGammaij
            FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                       - (eX * CosThetaj - OXj) * CosThetai)
            FXij = FXij * coeff
            FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                       - (eY * CosThetaj - OYj) * CosThetai)
            FYij = FYij * coeff
            FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                       - (eZ * CosThetaj - OZj) * CosThetai)
            FZij = FZij * coeff
            VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
            VirialLocalInter = VirialLocalInter + FXij * PXij + FYij * PYij + FZij * PZij
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            FX2(jk) = FX2(jk) - FXij
            FY2(jk) = FY2(jk) - FYij
            FZ2(jk) = FZ2(jk) - FZij
            TXi    = TXi    + Rij3Inv * (eX * CosThetaj3 - OXj)
            TYi    = TYi    + Rij3Inv * (eY * CosThetaj3 - OYj)
            TZi    = TZi    + Rij3Inv * (eZ * CosThetaj3 - OZj)
            TX2(jk) = TX2(jk) + Rij3Inv * (eX * CosThetai3 - OXi)
            TY2(jk) = TY2(jk) + Rij3Inv * (eY * CosThetai3 - OYi)
            TZ2(jk) = TZ2(jk) + Rij3Inv * (eZ * CosThetai3 - OZi)
!           TXi    = TXi    + Rij3Inv * (eX * CosThetaj3 - OXj) &
! &                         + RFConstant2 * OXj
!           TYi    = TYi    + Rij3Inv * (eY * CosThetaj3 - OYj) &
! &                         + RFConstant2 * OYj
!           TZi    = TZi    + Rij3Inv * (eZ * CosThetaj3 - OZj) &
! &                         + RFConstant2 * OZj
!           TX2(j) = TX2(j) + Rij3Inv * (eX * CosThetai3 - OXi) &
! &                         + RFConstant2 * OXi
!           TY2(j) = TY2(j) + Rij3Inv * (eY * CosThetai3 - OYi) &
! &                         + RFConstant2 * OYi
!           TZ2(j) = TZ2(j) + Rij3Inv * (eZ * CosThetai3 - OZi) &
! &                         + RFConstant2 * OZi
          end if
        end do loop1
        ! Include intramolecular interactions if need
        if (SameComponent .and. (intra15 .or. intra14)) then
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
            Rij3Inv = Epsilon * RijInv**3
            Rij4Inv3 = 3._RK * Rij3Inv * RijInv
            EPotLocal = EPotLocal +  Rij3Inv * Tmp
            EPotLocalIntra = EPotLocalIntra +  Rij3Inv * Tmp
            FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                       - (eX * CosThetaj - OXj) * CosThetai)
            FXij = FXij * coeff
            FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                       - (eY * CosThetaj - OYj) * CosThetai)
            FYij = FYij * coeff
            FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                       - (eZ * CosThetaj - OZj) * CosThetai)
            FZij = FZij * coeff
            VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
            VirialLocalIntra = VirialLocalIntra + FXij * PXij + FYij * PYij + FZij * PZij
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            FX2(i) = FX2(i) - FXij
            FY2(i) = FY2(i) - FYij
            FZ2(i) = FZ2(i) - FZij
            TXi    = TXi    + Rij3Inv * (eX * CosThetaj3 - OXj)
            TYi    = TYi    + Rij3Inv * (eY * CosThetaj3 - OYj)
            TZi    = TZi    + Rij3Inv * (eZ * CosThetaj3 - OZj)
            TX2(i) = TX2(i) + Rij3Inv * (eX * CosThetai3 - OXi)
            TY2(i) = TY2(i) + Rij3Inv * (eY * CosThetai3 - OYi)
            TZ2(i) = TZ2(i) + Rij3Inv * (eZ * CosThetai3 - OZi)
!           TXi    = TXi    + Rij3Inv * (eX * CosThetaj3 - OXj) &
! &                         + RFConstant2 * OXj
!           TYi    = TYi    + Rij3Inv * (eY * CosThetaj3 - OYj) &
! &                         + RFConstant2 * OYj
!           TZi    = TZi    + Rij3Inv * (eZ * CosThetaj3 - OZj) &
! &                         + RFConstant2 * OZj
!           TX2(j) = TX2(j) + Rij3Inv * (eX * CosThetai3 - OXi) &
! &                         + RFConstant2 * OXi
!           TY2(j) = TY2(j) + Rij3Inv * (eY * CosThetai3 - OYi) &
! &                         + RFConstant2 * OYi
!           TZ2(j) = TZ2(j) + Rij3Inv * (eZ * CosThetai3 - OZi) &
! &                         + RFConstant2 * OZi
          end if
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do

    else ! Site-site cutoff

      ! Loop over molecules
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
loop2:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )
#else
        j0 = merge( i + 1, 1, SameComponent )
loop2:  do j = j0, j1
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
          if( RijSquared >= RCutoffSquared ) cycle loop2
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
          EPotLocal = EPotLocal + Rij3Inv * Tmp - RFConstant2 * CosGammaij
          FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                     - (eX * CosThetaj - OXj) * CosThetai)
          FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                     - (eY * CosThetaj - OYj) * CosThetai)
          FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                     - (eZ * CosThetaj - OZj) * CosThetai)
          VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          FX2(j) = FX2(j) - FXij
          FY2(j) = FY2(j) - FYij
          FZ2(j) = FZ2(j) - FZij
          TXi    = TXi    + Rij3Inv * (eX * CosThetaj3 - OXj) &
&                         + RFConstant2 * OXj
          TYi    = TYi    + Rij3Inv * (eY * CosThetaj3 - OYj) &
&                         + RFConstant2 * OYj
          TZi    = TZi    + Rij3Inv * (eZ * CosThetaj3 - OZj) &
&                         + RFConstant2 * OZj
          TX2(j) = TX2(j) + Rij3Inv * (eX * CosThetai3 - OXi) &
&                         + RFConstant2 * OXi
          TY2(j) = TY2(j) + Rij3Inv * (eY * CosThetai3 - OYi) &
&                         + RFConstant2 * OYi
          TZ2(j) = TZ2(j) + Rij3Inv * (eZ * CosThetai3 - OZi) &
&                         + RFConstant2 * OZi
        end do loop2
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do

    end if

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal

    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter

    if ( IntraLJEl) then
      EPotIntra = EPotIntra + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if


  end subroutine TPotDD_Force



!==============================================================!
!  Subroutine TPotDD_ChemicalPotential                         !
!==============================================================!

  subroutine TPotDD_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDipoleDipole) :: this
    real(RK), pointer      :: EPotTest(:)
    real(RK), intent(in)   :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK)          :: RFConstant2
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
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

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
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
#if ARCH == 3
        hit = .false.
#endif
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
          RijSquared = RXij**2 + RYij**2 + RZij**2
#if ARCH == 3
          if( RijSquared <= RShieldSquared ) hit = .true.
#else
          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop1
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
!           EPotLocal = EPotLocal + Rij3Inv * Tmp - RFConstant2 * CosGammaij
          EPotLocal = EPotLocal + Rij3Inv * Tmp
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

    else ! Site-site cutoff

      ! Loop over test particles
      do i = 1, i1
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

    end if

  end subroutine TPotDD_ChemicalPotential



!==============================================================!
!  Subroutine TPotDD_Energy                                    !
!==============================================================!

  subroutine TPotDD_Energy( this, np, EPot, Virial, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDipoleDipole) :: this
    integer, intent(in)    :: np
    real(RK), pointer      :: EPot(:)
    real(RK), pointer      :: Virial(:)
    real(RK), intent(in)   :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK)          :: RFConstant2
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij3Inv, Rij4Inv3
    real(RK)          :: CosThetai, CosThetaj, CosGammaij
    real(RK)          :: Tmp
    real(RK)          :: EPotLocal
    integer           :: j, k, j1

    ! Assign local variables
    j1 = this%Site2%NPart
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared
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

    if( CutoffMode .eq. CenterofMass ) then

!CDIR NODEP
loop1:do k = 1, this%NInCutoff(np)
        j = this%CutoffPartner(k, np)
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
        RijSquared = RXij**2 + RYij**2 + RZij**2
        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
        else
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
          Tmp = CosGammaij -  3._RK * CosThetai * CosThetaj
          Rij3Inv = Epsilon * RijInv**3
          Rij4Inv3 = 3._RK * Rij3Inv * RijInv
          EPotLocal = Rij3Inv * Tmp
!           EPotLocal = Rij3Inv * Tmp - RFConstant2 * CosGammaij
          FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                     - (eX * CosThetaj - OXj) * CosThetai)
          FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                     - (eY * CosThetaj - OYj) * CosThetai)
          FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                     - (eZ * CosThetaj - OZj) * CosThetai)
        end if
        EPot(j) = EPot(j) + EPotLocal
        Virial(j) = Virial(j) + Third &
&                     * ( FXij * PXij + FYij * PYij + FZij * PZij )
      end do loop1

    else ! Site-site cutoff

!CDIR NODEP
loop2:do j = 1, j1
        if( this%SameComponent .and. j == np ) cycle loop2
        PXij = PXi - PX2(j)
        PYij = PYi - PY2(j)
        PZij = PZi - PZ2(j)
        RXij = RXi - RX2(j)
        RYij = RYi - RY2(j)
        RZij = RZi - RZ2(j)
        PXij = (PXij - anint( RXij )) * BoxLength
        PYij = (PYij - anint( RYij )) * BoxLength
        PZij = (PZij - anint( RZij )) * BoxLength
        RXij = (RXij - anint( RXij )) * BoxLength
        RYij = (RYij - anint( RYij )) * BoxLength
        RZij = (RZij - anint( RZij )) * BoxLength
        RijSquared = RXij**2 + RYij**2 + RZij**2
        if( RijSquared >= RCutoffSquared ) cycle loop2
        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
        else
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
          Tmp = CosGammaij -  3._RK * CosThetai * CosThetaj
          Rij3Inv = Epsilon * RijInv**3
          Rij4Inv3 = 3._RK * Rij3Inv * RijInv
          EPotLocal = Rij3Inv * Tmp - RFConstant2 * CosGammaij
          FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                     - (eX * CosThetaj - OXj) * CosThetai)
          FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                     - (eY * CosThetaj - OYj) * CosThetai)
          FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                     - (eZ * CosThetaj - OZj) * CosThetai)
        end if
        EPot(j) = EPot(j) + EPotLocal
        Virial(j) = Virial(j) + Third &
&                     * ( FXij * PXij + FYij * PYij + FZij * PZij )
      end do loop2

    end if


  end subroutine TPotDD_Energy


!==============================================================!
!  Subroutine TPotDQ_Construct                                 !
!==============================================================!

  subroutine TPotDQ_Construct( this, i1, i2, j1, j2, &
&                              Molecule1, Molecule2, RCutoff )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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
!    if (this%SameComponent .and. IntraLJEL .and. this%Site1%SiteId<this%Site2%SiteId) then
    if (this%SameComponent .and. IntraLJEL) then
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
    end if


  end subroutine TPotDQ_Construct



!==============================================================!
!  Subroutine TPotDQ_Destruct                                  !
!==============================================================!

  subroutine TPotDQ_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDipoleQuadrupole) :: this

    ! Destroy potential
    continue

  end subroutine TPotDQ_Destruct



!==============================================================!
!  Subroutine TPotDQ_Force                                     !
!==============================================================!

  subroutine TPotDQ_Force( this, EPot, Virial, EPotInter, VirialInter,EPotIntra, VirialIntra,BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDipoleQuadrupole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer :: TX1(:), TY1(:), TZ1(:), TX2(:), TY2(:), TZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
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
    real(RK)          :: EPotLocal1Inter, EPotLocalInter, VirialLocalInter
    real(RK)          :: EPotLocal1Intra, EPotLocalIntra, VirialLocalIntra
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
    integer           :: nu1, nu2, jk, unit
    logical           :: intra14, intra15
    real(RK)          :: coeff
#if MPI_VER > 0
    integer           :: N1, N2, i0, ji
    logical           :: EvenN
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
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    EPotLocal   = 0._RK
    VirialLocal = 0._RK
    EPotLocalInter   = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra   = 0._RK
    VirialLocalIntra = 0._RK

    intra14 = this%potintra14
    intra15 = this%potintra15
    if (intra14) then
      coeff = this%ScaleEl14 !Scale 1,4 El interactions
    else
      coeff = 1._RK
    end if

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
    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        unit=nu1*(i-1)+this%Site1%UnitNumber
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
loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
!          intra = this%Intra(k, unit)
!          scale = this%ScaleCoeff(k, unit)
 !         if (scale) then
 !            coeff = ScaleEl14
 !         else
             coeff = 1._RK
!          end if
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
            EPotLocal1 = Rij4Inv * (CosGammaij * CosThetaj &
&                        - CosThetai * (5._RK * CosThetaj2 - 1))
            EPotLocal = EPotLocal + EPotLocal1
!            if (intra) then
!               EPotLocal1Intra = Rij4Inv * (CosGammaij * CosThetaj &
!&                        - CosThetai * (5._RK * CosThetaj2 - 1))
!               EPotLocalIntra = EPotLocalIntra + EPotLocal1Intra
!            else
               EPotLocal1Inter = Rij4Inv * (CosGammaij * CosThetaj &
&                        - CosThetai * (5._RK * CosThetaj2 - 1))
               EPotLocalInter = EPotLocalInter + EPotLocal1Inter
!            end if
            dCosThetai = Rij4Inv * (1 - 5._RK * CosThetaj2)
            dCosThetaj = Rij4Inv * (CosGammaij - 10._RK * CosThetai * CosThetaj)
            dCosGammaij = 2._RK * Rij4Inv * CosThetaj
            Tmp = -4._RK * RijInv * EPotLocal1
            FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                      + (eX * CosThetaj - OXj) * dCosThetaj)
            FXij = FXij * coeff
            FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                      + (eY * CosThetaj - OYj) * dCosThetaj)
            FYij = FYij * coeff
            FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                      + (eZ * CosThetaj - OZj) * dCosThetaj)
            FZij = FZij * coeff
            VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
            VirialLocalInter = VirialLocalInter + FXij * PXij + FYij * PYij + FZij * PZij
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            FX2(jk) = FX2(jk) - FXij
            FY2(jk) = FY2(jk) - FYij
            FZ2(jk) = FZ2(jk) - FZij
            TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
            TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
            TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
            TX2(jk) = TX2(jk) - eX * dCosThetaj - OXi * dCosGammaij
            TY2(jk) = TY2(jk) - eY * dCosThetaj - OYi * dCosGammaij
            TZ2(jk) = TZ2(jk) - eZ * dCosThetaj - OZi * dCosGammaij
          end if
        end do loop1
        ! Include intramolecular interaction if need
        if (SameComponent .and. (intra15 .or. intra14)) then
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
            EPotLocal1 = Rij4Inv * (CosGammaij * CosThetaj &
&                        - CosThetai * (5._RK * CosThetaj2 - 1))
            EPotLocal = EPotLocal + EPotLocal1
            EPotLocal1Intra = Rij4Inv * (CosGammaij * CosThetaj &
&                        - CosThetai * (5._RK * CosThetaj2 - 1))
            EPotLocalIntra = EPotLocalIntra + EPotLocal1Intra
            dCosThetai = Rij4Inv * (1 - 5._RK * CosThetaj2)
            dCosThetaj = Rij4Inv * (CosGammaij - 10._RK * CosThetai * CosThetaj)
            dCosGammaij = 2._RK * Rij4Inv * CosThetaj
            Tmp = -4._RK * RijInv * EPotLocal1
            FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                      + (eX * CosThetaj - OXj) * dCosThetaj)
            FXij = FXij * coeff
            FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                      + (eY * CosThetaj - OYj) * dCosThetaj)
            FYij = FYij * coeff
            FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                      + (eZ * CosThetaj - OZj) * dCosThetaj)
            FZij = FZij * coeff
            VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
            VirialLocalIntra = VirialLocalIntra + FXij * PXij + FYij * PYij + FZij * PZij
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            FX2(i) = FX2(i) - FXij
            FY2(i) = FY2(i) - FYij
            FZ2(i) = FZ2(i) - FZij
            TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
            TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
            TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
            TX2(i) = TX2(i) - eX * dCosThetaj - OXi * dCosGammaij
            TY2(i) = TY2(i) - eY * dCosThetaj - OYi * dCosGammaij
            TZ2(i) = TZ2(i) - eZ * dCosThetaj - OZi * dCosGammaij
          end if
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do

    else ! Site-site cutoff

      ! Loop over molecules
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
loop2:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )
#else
        j0 = merge( i + 1, 1, SameComponent )
loop2:  do j = j0, j1
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
          if( RijSquared >= RCutoffSquared ) cycle loop2
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
          EPotLocal1 = Rij4Inv * (CosGammaij * CosThetaj &
&                      - CosThetai * (5._RK * CosThetaj2 - 1))
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
          VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          FX2(j) = FX2(j) - FXij
          FY2(j) = FY2(j) - FYij
          FZ2(j) = FZ2(j) - FZij
          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
          TX2(j) = TX2(j) - eX * dCosThetaj - OXi * dCosGammaij
          TY2(j) = TY2(j) - eY * dCosThetaj - OYi * dCosGammaij
          TZ2(j) = TZ2(j) - eZ * dCosThetaj - OZi * dCosGammaij
        end do loop2
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do

    end if

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal

    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter

    if (IntraLJEl) then
      EPotIntra = EPotIntra + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if

  end subroutine TPotDQ_Force



!==============================================================!
!  Subroutine TPotDQ_ChemicalPotential                         !
!==============================================================!

  subroutine TPotDQ_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDipoleQuadrupole) :: this
    real(RK), pointer          :: EPotTest(:)
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
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
#if ARCH == 3
    logical           :: hit
#endif

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

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
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
#if ARCH == 3
        hit = .false.
#endif
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
          RijSquared = RXij**2 + RYij**2 + RZij**2
#if ARCH == 3
          if( RijSquared <= RShieldSquared ) hit = .true.
#else
          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop1
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
          EPotLocal = EPotLocal + Rij4Inv * (2._RK * CosGammaij * CosThetaj &
&                               - CosThetai * (5._RK * CosThetaj**2 - 1._RK))
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

    else ! Site-site cutoff

      ! Loop over test particles
      do i = 1, i1
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
          EPotLocal = EPotLocal + Rij4Inv * (2._RK * CosGammaij * CosThetaj &
&                               - CosThetai * (5._RK * CosThetaj**2 - 1._RK))
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

    end if

  end subroutine TPotDQ_ChemicalPotential



!==============================================================!
!  Subroutine TPotDQ_Energy                                    !
!==============================================================!

  subroutine TPotDQ_Energy( this, np, EPot, Virial, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDipoleQuadrupole) :: this
    integer, intent(in)        :: np
    real(RK), pointer          :: EPot(:)
    real(RK), pointer          :: Virial(:)
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij4Inv
    real(RK)          :: CosThetai, CosThetaj, CosThetaj2, CosGammaij
    real(RK)          :: dCosThetai, dCosThetaj, dCosGammaij
    real(RK)          :: Tmp
    real(RK)          :: EPotLocal
    integer           :: j, k, j1

    ! Assign local variables
    j1 = this%Site2%NPart
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

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

    if( CutoffMode .eq. CenterofMass ) then

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
!CDIR NODEP
loop1:do k = 1, this%NInCutoff(np)
        j = this%CutoffPartner(k, np)
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
        RijSquared = RXij**2 + RYij**2 + RZij**2
        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
        else
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
          EPotLocal = Rij4Inv * (CosGammaij * CosThetaj &
&                      - CosThetai * (5._RK * CosThetaj2 - 1._RK))
          dCosThetai = Rij4Inv * (1._RK - 5._RK * CosThetaj2)
          dCosThetaj = Rij4Inv * (CosGammaij - 10._RK * CosThetai * CosThetaj)
          dCosGammaij = 2._RK * Rij4Inv * CosThetaj
          Tmp = -4._RK * RijInv * EPotLocal
          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)
        end if
        EPot(j) = EPot(j) + EPotLocal
        Virial(j) = Virial(j) + Third &
&                     * ( FXij * PXij + FYij * PYij + FZij * PZij )
      end do loop1

    else ! Site-site cutoff

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
!CDIR NODEP
loop2:do j = 1, j1
        if( this%SameComponent .and. j == np ) cycle loop2
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
        if( RijSquared >= RCutoffSquared ) cycle loop2
        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
        else
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
          EPotLocal = Rij4Inv * (CosGammaij * CosThetaj &
&                      - CosThetai * (5._RK * CosThetaj2 - 1._RK))
          dCosThetai = Rij4Inv * (1._RK - 5._RK * CosThetaj2)
          dCosThetaj = Rij4Inv * (CosGammaij - 10._RK * CosThetai * CosThetaj)
          dCosGammaij = 2._RK * Rij4Inv * CosThetaj
          Tmp = -4._RK * RijInv * EPotLocal
          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)
        end if
        EPot(j) = EPot(j) + EPotLocal
        Virial(j) = Virial(j) + Third &
&                     * ( FXij * PXij + FYij * PYij + FZij * PZij )
      end do loop2

    end if

  end subroutine TPotDQ_Energy


!==============================================================!
!  Subroutine TPotQC_Construct                                 !
!==============================================================!

  subroutine TPotQC_Construct( this, i1, i2, j1, j2, &
&                              Molecule1, Molecule2, RCutoff )
    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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

!    ! if this potential is intra
!    if (this%SameComponent .and. IntraLJEL .and. this%Site1%SiteId<this%Site2%SiteId) then
    if (this%SameComponent .and. IntraLJEL) then
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
    end if


  end subroutine TPotQC_Construct



!==============================================================!
!  Subroutine TPotQC_Destruct                                  !
!==============================================================!

  subroutine TPotQC_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotQuadrupoleCharge) :: this

    ! Destroy potential
    continue

  end subroutine TPotQC_Destruct



!==============================================================!
!  Subroutine TPotQC_Force                                     !
!==============================================================!

  subroutine TPotQC_Force( this, EPot, Virial, EPotInter, VirialInter,EPotIntra, VirialIntra, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotQuadrupoleCharge) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in)       :: BoxLength


    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon1, Epsilon2
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:)
    real(RK), pointer :: TX1(:), TY1(:), TZ1(:)
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
    integer           :: i, j, k, i1
    integer           :: nu1, nu2, jk, unit
    logical           :: intra14, intra15, SameComponent
    real(RK)          :: coeff
#if MPI_VER > 0
    integer           :: i0
#endif

    ! Assign local variables
#if MPI_VER > 0
    i0 = this%Site1%NPart0
    i1 = this%Site1%NPart2
#else
    i1 = this%Site1%NPart
#endif
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    Epsilon = this%Epsilon
    EPotLocal = 0._RK
    VirialLocal = 0._RK
    EPotLocalInter = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra = 0._RK
    VirialLocalIntra = 0._RK

    intra14 = this%potintra14
    intra15 = this%potintra15
    SameComponent = this%SameComponent
    if (intra14) then
      coeff = this%ScaleEl14 !Scale 1,4 El interactions
    else
      coeff = 1._RK
    end if



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
    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

    ! Loop over molecules
#if MPI_VER > 0
    do i = i0, i1
#else
    do i = 1, i1
#endif
        unit=nu1*(i-1)+this%Site1%UnitNumber
!        print *, 'Unit =', unit
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
loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
!          intra = this%Intra(k, unit)
!          scale = this%ScaleCoeff(k, unit)
!          if (scale) then
!             coeff = ScaleEl14
!          else
             coeff = 1._RK
!          end if
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
            eX = - RXij * RijInv                                                    ! Normierter Abstandsvektor nach Price
            eY = - RYij * RijInv
            eZ = - RZij * RijInv
            CosTheta  = OXi * ex + OYi * eY + OZi * eZ                              ! Scalarprodukt normierter Abstandsvektor mit Orientierungsvektor Quadrupol
            Epsilon1 = Epsilon * RijSquaredInv * RijInv
            EPotLocal  = EPotLocal + Epsilon1 * ( CosTheta * CosTheta - Third )
!            if (intra) then
!              EPotLocalIntra  = EPotLocalIntra + Epsilon1 * ( CosTheta * CosTheta - Third )!
!	    else
              EPotLocalInter  = EPotLocalInter + Epsilon1 * ( CosTheta * CosTheta - Third )
!            end if
            CosTheta2 = 2._RK * CosTheta
            CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
            Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
            FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXi ) * coeff                     ! Kraft auf die Punktladung, sprich F2
            FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYi ) * coeff
            FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZi ) * coeff
            VirialLocal = VirialLocal - FXij * PXij - FYij * PYij - FZij * PZij     ! Vorzeichen richtig
            VirialLocalInter = VirialLocalInter - FXij * PXij - FYij * PYij - FZij * PZij     ! Vorzeichen richtig
            FXi    = FXi    - FXij
            FYi    = FYi    - FYij
            FZi    = FZi    - FZij
            FX2(jk) = FX2(jk) + FXij
            FY2(jk) = FY2(jk) + FYij
            FZ2(jk) = FZ2(jk) + FZij
            TXi    = TXi - Epsilon1*CosTheta2*eX                                    ! Drehmomentanteil auf Quadrupol wegen Punktladung. Kreuzprodukt
            TYi    = TYi - Epsilon1*CosTheta2*eY                                    ! in Atom2Mol von Component
            TZi    = TZi - Epsilon1*CosTheta2*eZ
          end if
        end do loop1
        ! Include intramolecular interaction if need
        if (SameComponent .and. (intra15 .or. intra14)) then
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
            eX = - RXij * RijInv                                                    ! Normierter Abstandsvektor nach Price
            eY = - RYij * RijInv
            eZ = - RZij * RijInv
            CosTheta  = OXi * ex + OYi * eY + OZi * eZ                              ! Scalarprodukt normierter Abstandsvektor mit Orientierungsvektor Quadrupol
            Epsilon1 = Epsilon * RijSquaredInv * RijInv
            EPotLocal  = EPotLocal + Epsilon1 * ( CosTheta * CosTheta - Third )
            EPotLocalIntra  = EPotLocalIntra + Epsilon1 * ( CosTheta * CosTheta - Third )!
            CosTheta2 = 2._RK * CosTheta
            CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
            Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
            FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXi ) * coeff                     ! Kraft auf die Punktladung, sprich F2
            FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYi ) * coeff
            FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZi ) * coeff
            VirialLocal = VirialLocal - FXij * PXij - FYij * PYij - FZij * PZij     ! Vorzeichen richtig
            VirialLocalIntra = VirialLocalIntra - FXij * PXij - FYij * PYij - FZij * PZij     ! Vorzeichen richtig
            FXi    = FXi    - FXij
            FYi    = FYi    - FYij
            FZi    = FZi    - FZij
            FX2(i) = FX2(i) + FXij
            FY2(i) = FY2(i) + FYij
            FZ2(i) = FZ2(i) + FZij
            TXi    = TXi - Epsilon1*CosTheta2*eX                                    ! Drehmomentanteil auf Quadrupol wegen Punktladung. Kreuzprodukt
            TYi    = TYi - Epsilon1*CosTheta2*eY                                    ! in Atom2Mol von Component
            TZi    = TZi - Epsilon1*CosTheta2*eZ
          end if
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
    end do

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal

    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter

    if (IntraLJEl) then
      EPotIntra = EPotIntra + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if

  end subroutine TPotQC_Force



!==============================================================!
!  Subroutine TPotQC_ChemicalPotential                         !
!==============================================================!

  subroutine TPotQC_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotQuadrupoleCharge) :: this
    real(RK), pointer          :: EPotTest(:)
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:)
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
#if ARCH == 3
    logical           :: hit
#endif

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
#if ARCH == 3
        hit = .false.
#endif
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
          EPotLocal = EPotLocal + Epsilon * RijSquaredInv * RijInv &
&                       * ( CosTheta * CosTheta - Third )
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

  subroutine TPotQC_Energy( this, np, EPot, Virial, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotQuadrupoleCharge) :: this
    integer, intent(in)        :: np
    real(RK), pointer          :: EPot(:)
    real(RK), pointer          :: Virial(:)
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon, Epsilon2
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquaredInv, RijInv, RijSquared
    real(RK)          :: CosTheta, CosTheta2, CosAux
    real(RK)          :: EPotLocal
    integer           :: j, k

    ! Assign local variables
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

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

!CDIR NODEP
    do k = 1, this%NInCutoff(np)
      j = this%CutoffPartner(k, np)
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
      RijSquared = RXij**2 + RYij**2 + RZij**2
      if( RijSquared <= RShieldSquared ) then
        EPotLocal = 1E33_RK
      else
        RijSquaredInv = 1._RK / RijSquared
        RijInv = sqrt( RijSquaredInv )
        eX = - RXij * RijInv                                                    ! Normierter Abstandsvektor nach Price
        eY = - RYij * RijInv
        eZ = - RZij * RijInv
        CosTheta  = OXi * ex + OYi * eY + OZi * eZ                              ! Scalarprodukt normierter Abstandsvektor mit Orientierungsvektor Quadrupol
        EPotLocal  = Epsilon * RijSquaredInv * RijInv &
&                      * ( CosTheta * CosTheta - Third )
        CosTheta2 = 2._RK * CosTheta
        CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
        Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
        FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXi )                     ! Kraft auf die Punktladung, sprich F2
        FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYi )
        FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZi )
      end if
      EPot(j) = EPot(j) + EPotLocal
      Virial(j) = Virial(j) - Third &
&                   * ( FXij * PXij + FYij * PYij + FZij * PZij )
    end do

  end subroutine TPotQC_Energy



!==============================================================!
!  Subroutine TPotQD_Construct                                 !
!==============================================================!

  subroutine TPotQD_Construct( this, i1, i2, j1, j2, &
&                              Molecule1, Molecule2, RCutoff )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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

!   ! if this potential is intra
!   if (this%SameComponent .and. IntraLJEL .and. this%Site1%SiteId<this%Site2%SiteId) then
   if (this%SameComponent .and. IntraLJEL ) then
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
    end if

  end subroutine TPotQD_Construct



!==============================================================!
!  Subroutine TPotQD_Destruct                                  !
!==============================================================!

  subroutine TPotQD_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotQuadrupoleDipole) :: this

    ! Destroy potential
    continue

  end subroutine TPotQD_Destruct



!==============================================================!
!  Subroutine TPotQD_Force                                     !
!==============================================================!

  subroutine TPotQD_Force( this, EPot, Virial, EPotInter, VirialInter,EPotIntra, VirialIntra, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotQuadrupoleDipole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer :: TX1(:), TY1(:), TZ1(:), TX2(:), TY2(:), TZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
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
    real(RK)          :: EPotLocal1Inter, EPotLocalInter, VirialLocalInter
    real(RK)          :: EPotLocal1Intra, EPotLocalIntra, VirialLocalIntra
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
    integer           :: nu1, nu2, jk, unit
    logical           :: intra14, intra15
    real(RK)          :: coeff
#if MPI_VER > 0
    integer           :: N1, N2, i0, ji
    logical           :: EvenN
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
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    EPotLocal   = 0._RK
    VirialLocal = 0._RK
    EPotLocalInter   = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra   = 0._RK
    VirialLocalIntra = 0._RK

    intra14 = this%potintra14
    intra15 = this%potintra15
    if (intra14) then
      coeff = this%ScaleEl14 !Scale 1,4 El interactions
    else
      coeff = 1._RK
    end if



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
    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        unit=nu1*(i-1)+this%Site1%UnitNumber
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
loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
!          intra = this%Intra(k, unit)
!          scale = this%ScaleCoeff(k, unit)
!          if (scale) then
!             coeff = ScaleEl14
!          else
             coeff = 1._RK
!          end if
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
            EPotLocal1 = Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) &
&                                    - CosGammaij * CosThetai)
            EPotLocal = EPotLocal + EPotLocal1
!            if (intra) then
!              EPotLocal1Intra = Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) &
!&                                    - CosGammaij * CosThetai)
!              EPotLocalIntra = EPotLocalIntra + EPotLocal1Intra
!	    else
              EPotLocal1Inter = Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) &
&                                    - CosGammaij * CosThetai)
              EPotLocalInter = EPotLocalInter + EPotLocal1Inter
!	    end if
            dCosThetai = Rij4Inv * (10._RK * CosThetai * CosThetaj - CosGammaij)
            dCosThetaj = Rij4Inv * (5._RK * CosThetai2 - 1._RK)
            dCosGammaij = -2._RK * Rij4Inv * CosThetai
            Tmp = -4._RK * RijInv * EPotLocal1
            FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                      + (eX * CosThetaj - OXj) * dCosThetaj)
            FXij = FXij * coeff
            FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                      + (eY * CosThetaj - OYj) * dCosThetaj)
            FYij = FYij * coeff
            FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                      + (eZ * CosThetaj - OZj) * dCosThetaj)
            FZij = FZij * coeff
            VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
            VirialLocalInter = VirialLocalInter + FXij * PXij + FYij * PYij + FZij * PZij
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            FX2(jk) = FX2(jk) - FXij
            FY2(jk) = FY2(jk) - FYij
            FZ2(jk) = FZ2(jk) - FZij
            TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
            TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
            TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
            TX2(jk) = TX2(jk) - eX * dCosThetaj - OXi * dCosGammaij
            TY2(jk) = TY2(jk) - eY * dCosThetaj - OYi * dCosGammaij
            TZ2(jk) = TZ2(jk) - eZ * dCosThetaj - OZi * dCosGammaij
          end if
        end do loop1
        ! Include intramolecular interaction if need
        if (SameComponent .and. (intra14 .or. intra15)) then
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
&                                    - CosGammaij * CosThetai)
            EPotLocal = EPotLocal + EPotLocal1
            EPotLocal1Intra = Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) &
&                                    - CosGammaij * CosThetai)
            EPotLocalIntra = EPotLocalIntra + EPotLocal1Intra
            dCosThetai = Rij4Inv * (10._RK * CosThetai * CosThetaj - CosGammaij)
            dCosThetaj = Rij4Inv * (5._RK * CosThetai2 - 1._RK)
            dCosGammaij = -2._RK * Rij4Inv * CosThetai
            Tmp = -4._RK * RijInv * EPotLocal1
            FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                      + (eX * CosThetaj - OXj) * dCosThetaj)
            FXij = FXij * coeff
            FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                      + (eY * CosThetaj - OYj) * dCosThetaj)
            FYij = FYij * coeff
            FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                      + (eZ * CosThetaj - OZj) * dCosThetaj)
            FZij = FZij * coeff
            VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
            VirialLocalIntra = VirialLocalIntra + FXij * PXij + FYij * PYij + FZij * PZij
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            FX2(i) = FX2(i) - FXij
            FY2(i) = FY2(i) - FYij
            FZ2(i) = FZ2(i) - FZij
            TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
            TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
            TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
            TX2(i) = TX2(i) - eX * dCosThetaj - OXi * dCosGammaij
            TY2(i) = TY2(i) - eY * dCosThetaj - OYi * dCosGammaij
            TZ2(i) = TZ2(i) - eZ * dCosThetaj - OZi * dCosGammaij
          end if
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do

    else ! Site-site cutoff

      ! Loop over molecules
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
loop2:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )
#else
        j0 = merge( i + 1, 1, SameComponent )
!CDIR NODEP
loop2:  do j = j0, j1
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
          if( RijSquared >= RCutoffSquared ) cycle loop2
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
          EPotLocal1 = Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) &
&                                  - CosGammaij * CosThetai)
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
          VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          FX2(j) = FX2(j) - FXij
          FY2(j) = FY2(j) - FYij
          FZ2(j) = FZ2(j) - FZij
          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
          TX2(j) = TX2(j) - eX * dCosThetaj - OXi * dCosGammaij
          TY2(j) = TY2(j) - eY * dCosThetaj - OYi * dCosGammaij
          TZ2(j) = TZ2(j) - eZ * dCosThetaj - OZi * dCosGammaij
        end do loop2
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do

    end if

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal

    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter

    if (IntraLJEl) then
      EPotIntra = EPotIntra + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
    end if

  end subroutine TPotQD_Force



!==============================================================!
!  Subroutine TPotQD_ChemicalPotential                         !
!==============================================================!

  subroutine TPotQD_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotQuadrupoleDipole) :: this
    real(RK), pointer          :: EPotTest(:)
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
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
#if ARCH == 3
    logical           :: hit
#endif

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

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
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
#if ARCH == 3
        hit = .false.
#endif
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
          RijSquared = RXij**2 + RYij**2 + RZij**2
#if ARCH == 3
          if( RijSquared <= RShieldSquared ) hit = .true.
#else
          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop1
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
          EPotLocal = EPotLocal + Rij4Inv * &
&           (CosThetaj * (5._RK * CosThetai**2 - 1._RK) - &
&           2._RK * CosGammaij * CosThetai)
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

    else ! Site-site cutoff

      ! Loop over test particles
      do i = 1, i1
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
          EPotLocal = EPotLocal + Rij4Inv &
&                       * ( CosThetaj * (5._RK * CosThetai**2 - 1._RK) &
&                            - 2._RK * CosGammaij * CosThetai )
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

    end if

  end subroutine TPotQD_ChemicalPotential



!==============================================================!
!  Subroutine TPotQD_Energy                                    !
!==============================================================!

  subroutine TPotQD_Energy( this, np, EPot, Virial, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotQuadrupoleDipole) :: this
    integer, intent(in)        :: np
    real(RK), pointer          :: EPot(:)
    real(RK), pointer          :: Virial(:)
    real(RK), intent(in)       :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij4Inv
    real(RK)          :: CosThetai, CosThetaj, CosThetai2, CosGammaij
    real(RK)          :: dCosThetai, dCosThetaj, dCosGammaij
    real(RK)          :: tmp
    real(RK)          :: EPotLocal
    integer           :: j, k, j1

    ! Assign local variables
    j1 = this%Site2%NPart
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

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

    if( CutoffMode .eq. CenterofMass ) then

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
!CDIR NODEP
loop1:do k = 1, this%NInCutoff(np)
        j = this%CutoffPartner(k, np)
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
        RijSquared = RXij**2 + RYij**2 + RZij**2
        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
        else
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
          EPotLocal = Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) &
&                              - CosGammaij * CosThetai)
          dCosThetai = Rij4Inv * (10._RK * CosThetai * CosThetaj - CosGammaij)
          dCosThetaj = Rij4Inv * (5._RK * CosThetai2 - 1._RK)
          dCosGammaij = -2._RK * Rij4Inv * CosThetai
          Tmp = -4._RK * RijInv * EPotLocal
          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)
        end if
        EPot(j) = EPot(j) + EPotLocal
        Virial(j) = Virial(j) + Third &
&                     * (FXij * PXij + FYij * PYij + FZij * PZij)
      end do loop1

    else ! Site-site cutoff

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
!CDIR NODEP
loop2:do j = 1, j1
        if( this%SameComponent .and. j == np ) cycle loop2
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
        if( RijSquared >= RCutoffSquared ) cycle loop2
        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
        else
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
          EPotLocal = Rij4Inv * (CosThetaj * (5._RK * CosThetai2 - 1._RK) &
&                              - CosGammaij * CosThetai)
          dCosThetai = Rij4Inv * (10._RK * CosThetai * CosThetaj - CosGammaij)
          dCosThetaj = Rij4Inv * (5._RK * CosThetai2 - 1._RK)
          dCosGammaij = -2._RK * Rij4Inv * CosThetai
          Tmp = -4._RK * RijInv * EPotLocal
          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)
        end if
        EPot(j) = EPot(j) + EPotLocal
        Virial(j) = Virial(j) + Third &
&                     * (FXij * PXij + FYij * PYij + FZij * PZij)
      end do loop2

    end if

  end subroutine TPotQD_Energy



!==============================================================!
!  Subroutine TPotQQ_Construct                                 !
!==============================================================!

  subroutine TPotQQ_Construct( this, i1, i2, j1, j2, &
&                              Molecule1, Molecule2, &
&                              RCutoff )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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
!    if (this%SameComponent .and. IntraLJEL .and. this%Site1%SiteId<this%Site2%SiteId) then
    if (this%SameComponent .and. IntraLJEL) then
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
    end if


  end subroutine TPotQQ_Construct



!==============================================================!
!  Subroutine TPotQQ_Destruct                                  !
!==============================================================!

  subroutine TPotQQ_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotQuadrupoleQuadrupole) :: this

    ! Destroy potential
    continue

  end subroutine TPotQQ_Destruct



!==============================================================!
!  Subroutine TPotQQ_Force                                     !
!==============================================================!

  subroutine TPotQQ_Force( this, EPot, Virial, EPotInter, VirialInter, EPotIntra, VirialIntra, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotQuadrupoleQuadrupole) :: this
    real(RK), intent(in out)       :: EPot
    real(RK), intent(in out)       :: Virial
    real(RK), intent(in out) :: EPotInter
    real(RK), intent(in out) :: VirialInter
    real(RK), intent(in out) :: EPotIntra
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in)           :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer :: TX1(:), TY1(:), TZ1(:), TX2(:), TY2(:), TZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
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
    real(RK)          :: EPotLocal1Intra, EPotLocalIntra, VirialLocalIntra
    real(RK)          :: EPotLocal1Inter, EPotLocalInter, VirialLocalInter
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
    integer           :: nu1, nu2, jk, unit
    logical           :: intra14, intra15
    real(RK)          :: coeff
#if MPI_VER > 0
    integer           :: N1, N2, i0, ji
    logical           :: EvenN
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
    nu1 = this%NUnit1
    nu2 = this%NUnit2
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    EPotLocal   = 0._RK
    VirialLocal = 0._RK
    EPotLocalInter   = 0._RK
    VirialLocalInter = 0._RK
    EPotLocalIntra   = 0._RK
    VirialLocalIntra = 0._RK

    intra14 = this%potintra14
    intra15 = this%potintra15
    if (intra14) then
      coeff = this%ScaleEl14 !Scale 1,4 El interactions
    else
      coeff = 1._RK
    end if

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
    FX2 => this%Site2%FX
    FY2 => this%Site2%FY
    FZ2 => this%Site2%FZ
    TX1 => this%Site1%TX
    TY1 => this%Site1%TY
    TZ1 => this%Site1%TZ
    TX2 => this%Site2%TX
    TY2 => this%Site2%TY
    TZ2 => this%Site2%TZ
    PX1 => this%Site1%PX
    PY1 => this%Site1%PY
    PZ1 => this%Site1%PZ
    PX2 => this%Site2%PX
    PY2 => this%Site2%PY
    PZ2 => this%Site2%PZ

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over molecules
#if MPI_VER > 0
      do i = i0, i1
#else
      do i = 1, i1
#endif
        unit=nu1*(i-1)+this%Site1%UnitNumber
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
loop1:  do k = 1, this%NInCutoff(unit)
          j = this%CutoffPartner(k, unit)
!          intra = this%Intra(k, unit)
!          scale = this%ScaleCoeff(k, unit)
!          if (scale) then
!             coeff = ScaleEl14
!          else
             coeff = 1._RK
!          end if
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
            EPotLocal1 = Rij5Inv * (1._RK &
&             - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&             - 15._RK * CosThetaiSquared * CosThetajSquared &
&             + 2._RK * Tmp**2)
            EPotLocal = EPotLocal + EPotLocal1
!            if (intra) then
!              EPotLocal1Intra = Rij5Inv * (1._RK &
!&             - 5._RK * (CosThetaiSquared + CosThetajSquared) &
!&             - 15._RK * CosThetaiSquared * CosThetajSquared &
!&             + 2._RK * Tmp**2)
!              EPotLocalIntra = EPotLocalIntra + EPotLocal1Intra
!	    else
              EPotLocal1Inter = Rij5Inv * (1._RK &
&             - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&             - 15._RK * CosThetaiSquared * CosThetajSquared &
&             + 2._RK * Tmp**2)
              EPotLocalInter = EPotLocalInter + EPotLocal1Inter
!	    end if
            dCosThetai = Rij5Inv * (-10._RK * CosThetai &
&                                  - 30._RK * CosThetai * CosThetajSquared &
&                                  - 20._RK * CosThetaj * Tmp)
            dCosThetaj = Rij5Inv * (-10._RK * CosThetaj &
&                                  - 30._RK * CosThetaj * CosThetaiSquared &
&                                  - 20._RK * CosThetai * Tmp)
            dCosGammaij = 4._RK * Rij5Inv * Tmp
            Tmp = -5._RK * RijInv * EPotLocal1
            FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                      + (eX * CosThetaj - OXj) * dCosThetaj)
            FXij = FXij * coeff
            FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                      + (eY * CosThetaj - OYj) * dCosThetaj)
            FYij = FYij * coeff
            FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                      + (eZ * CosThetaj - OZj) * dCosThetaj)
            FZij = FZij * coeff
            VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
            VirialLocalInter = VirialLocalInter + FXij * PXij + FYij * PYij + FZij * PZij
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            FX2(jk) = FX2(jk) - FXij
            FY2(jk) = FY2(jk) - FYij
            FZ2(jk) = FZ2(jk) - FZij
            TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
            TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
            TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
            TX2(jk) = TX2(jk) - eX * dCosThetaj - OXi * dCosGammaij
            TY2(jk) = TY2(jk) - eY * dCosThetaj - OYi * dCosGammaij
            TZ2(jk) = TZ2(jk) - eZ * dCosThetaj - OZi * dCosGammaij
          end if
        end do loop1
        ! Include intramolecular interaction if need
        if (SameComponent .and. (intra14 .or. intra15)) then
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
            EPotLocal1 = Rij5Inv * (1._RK &
&             - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&             - 15._RK * CosThetaiSquared * CosThetajSquared &
&             + 2._RK * Tmp**2)
            EPotLocal = EPotLocal + EPotLocal1
            EPotLocal1Intra = Rij5Inv * (1._RK &
&             - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&             - 15._RK * CosThetaiSquared * CosThetajSquared &
&             + 2._RK * Tmp**2)
            EPotLocalIntra = EPotLocalIntra + EPotLocal1Intra
            dCosThetai = Rij5Inv * (-10._RK * CosThetai &
&                                  - 30._RK * CosThetai * CosThetajSquared &
&                                  - 20._RK * CosThetaj * Tmp)
            dCosThetaj = Rij5Inv * (-10._RK * CosThetaj &
&                                  - 30._RK * CosThetaj * CosThetaiSquared &
&                                  - 20._RK * CosThetai * Tmp)
            dCosGammaij = 4._RK * Rij5Inv * Tmp
            Tmp = -5._RK * RijInv * EPotLocal1
            FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                      + (eX * CosThetaj - OXj) * dCosThetaj)
            FXij = FXij * coeff
            FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                      + (eY * CosThetaj - OYj) * dCosThetaj)
            FYij = FYij * coeff
            FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                      + (eZ * CosThetaj - OZj) * dCosThetaj)
            FZij = FZij * coeff
            VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
            VirialLocalIntra = VirialLocalIntra + FXij * PXij + FYij * PYij + FZij * PZij
            FXi    = FXi    + FXij
            FYi    = FYi    + FYij
            FZi    = FZi    + FZij
            FX2(i) = FX2(i) - FXij
            FY2(i) = FY2(i) - FYij
            FZ2(i) = FZ2(i) - FZij
            TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
            TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
            TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
            TX2(i) = TX2(i) - eX * dCosThetaj - OXi * dCosGammaij
            TY2(i) = TY2(i) - eY * dCosThetaj - OYi * dCosGammaij
            TZ2(i) = TZ2(i) - eZ * dCosThetaj - OZi * dCosGammaij
          end if
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do
   else ! Site-site cutoff

      ! Loop over molecules
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
loop2:  do ji = j0, j1
          j = 1 + mod( ji - 1, N1 )
#else
        j0 = merge( i + 1, 1, SameComponent )
!CDIR NODEP
loop2:  do j = j0, j1
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
          if( RijSquared >= RCutoffSquared ) cycle loop2
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
          EPotLocal1 = Rij5Inv * (1._RK &
&           - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&           - 15._RK * CosThetaiSquared * CosThetajSquared &
&           + 2._RK * Tmp**2)
          EPotLocal = EPotLocal + EPotLocal1
          dCosThetai = Rij5Inv * (-10._RK * CosThetai &
&                                - 30._RK * CosThetai * CosThetajSquared &
&                                - 20._RK * CosThetaj * Tmp)
          dCosThetaj = Rij5Inv * (-10._RK * CosThetaj &
&                                - 30._RK * CosThetaj * CosThetaiSquared &
&                                - 20._RK * CosThetai * Tmp)
          dCosGammaij = 4._RK * Rij5Inv * Tmp
          Tmp = -5._RK * RijInv * EPotLocal1
          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)
          VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij
          FXi    = FXi    + FXij
          FYi    = FYi    + FYij
          FZi    = FZi    + FZij
          FX2(j) = FX2(j) - FXij
          FY2(j) = FY2(j) - FYij
          FZ2(j) = FZ2(j) - FZij
          TXi    = TXi    - eX * dCosThetai - OXj * dCosGammaij
          TYi    = TYi    - eY * dCosThetai - OYj * dCosGammaij
          TZi    = TZi    - eZ * dCosThetai - OZj * dCosGammaij
          TX2(j) = TX2(j) - eX * dCosThetaj - OXi * dCosGammaij
          TY2(j) = TY2(j) - eY * dCosThetaj - OYi * dCosGammaij
          TZ2(j) = TZ2(j) - eZ * dCosThetaj - OZi * dCosGammaij
        end do loop2
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
        TX1(i) = TXi
        TY1(i) = TYi
        TZ1(i) = TZi
      end do

    end if

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal

    EPotInter = EPotInter + EPotLocalInter
    VirialInter = VirialInter + Third * VirialLocalInter

    if (IntraLJEl) then
      EPotIntra = EPotIntra + EPotLocalIntra
      VirialIntra = VirialIntra + Third * VirialLocalIntra
   end if

  end subroutine TPotQQ_Force



!==============================================================!
!  Subroutine TPotQQ_ChemicalPotential                         !
!==============================================================!

  subroutine TPotQQ_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotQuadrupoleQuadrupole) :: this
    real(RK), pointer              :: EPotTest(:)
    real(RK), intent(in)           :: BoxLength

    ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
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
#if ARCH == 3
    logical           :: hit
#endif

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

    if( CutoffMode .eq. CenterofMass ) then

      ! Loop over test particles
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
#if ARCH == 3
        hit = .false.
#endif
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
          RijSquared = RXij**2 + RYij**2 + RZij**2
#if ARCH == 3
          if( RijSquared <= RShieldSquared ) hit = .true.
#else
          if( RijSquared <= RShieldSquared ) then
            EPotLocal = 1E33_RK
            exit loop1
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
          EPotLocal = EPotLocal + Rij5Inv * (1._RK &
&           - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&           - 15._RK * CosThetaiSquared * CosThetajSquared &
&           + 2._RK * Tmp**2)
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

    else ! Site-site cutoff

      ! Loop over test particles
      do i = 1, i1
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
          EPotLocal = EPotLocal + Rij5Inv * (1._RK &
&           - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&           - 15._RK * CosThetaiSquared * CosThetajSquared &
&           + 2._RK * Tmp**2)
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

    end if

  end subroutine TPotQQ_ChemicalPotential



!==============================================================!
!  Subroutine TPotQQ_Energy                                    !
!==============================================================!

  subroutine TPotQQ_Energy( this, np, EPot, Virial, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotQuadrupoleQuadrupole) :: this
    integer, intent(in)            :: np
    real(RK), pointer              :: EPot(:)
    real(RK), pointer              :: Virial(:)
    real(RK), intent(in)           :: BoxLength

        ! Declare local variables
    real(RK)          :: Epsilon
    real(RK)          :: RCutoffSquared
    real(RK)          :: RShieldSquared
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: FXij, FYij, FZij
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, Rij5Inv
    real(RK)          :: CosThetai, CosThetaj, CosGammaij
    real(RK)          :: CosThetaiSquared, CosThetajSquared
    real(RK)          :: dCosThetai, dCosThetaj, dCosGammaij
    real(RK)          :: Tmp
    real(RK)          :: EPotLocal
    integer           :: j, k, j1

    ! Assign local variables
    j1 = this%Site2%NPart
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RShieldSquared = this%RShieldSquared

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

    if( CutoffMode .eq. CenterofMass ) then

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
!CDIR NODEP
loop1:do k = 1, this%NInCutoff(np)
        j = this%CutoffPartner(k, np)
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
        RijSquared = RXij**2 + RYij**2 + RZij**2
        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
        else
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
          EPotLocal = Rij5Inv * (1._RK &
&           - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&           - 15._RK * CosThetaiSquared * CosThetajSquared &
&           + 2._RK * Tmp**2)

          dCosThetai = Rij5Inv * (-10._RK * CosThetai &
&                                - 30._RK * CosThetai * CosThetajSquared &
&                                - 20._RK * CosThetaj * Tmp)
          dCosThetaj = Rij5Inv * (-10._RK * CosThetaj &
&                                - 30._RK * CosThetaj * CosThetaiSquared &
&                                - 20._RK * CosThetai * Tmp)
          dCosGammaij = 4._RK * Rij5Inv * Tmp
          Tmp = -5._RK * RijInv * EPotLocal
          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)
        end if
        EPot(j) = EPot(j) + EPotLocal
        Virial(j) = Virial(j) + Third &
&                     * ( FXij * PXij + FYij * PYij + FZij * PZij )
      end do loop1

    else ! Site-site cutoff

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
!CDIR NODEP
loop2:do j = 1, j1
        if( this%SameComponent .and. j == np ) cycle loop2
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
        if( RijSquared >= RCutoffSquared ) cycle loop2
        if( RijSquared <= RShieldSquared ) then
          EPotLocal = 1E33_RK
        else
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
          EPotLocal = Rij5Inv * (1._RK &
&           - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&           - 15._RK * CosThetaiSquared * CosThetajSquared &
&           + 2._RK * Tmp**2)
          dCosThetai = Rij5Inv * (-10._RK * CosThetai &
&                                - 30._RK * CosThetai * CosThetajSquared &
&                                - 20._RK * CosThetaj * Tmp)
          dCosThetaj = Rij5Inv * (-10._RK * CosThetaj &
&                                - 30._RK * CosThetaj * CosThetaiSquared &
&                                - 20._RK * CosThetai * Tmp)
          dCosGammaij = 4._RK * Rij5Inv * Tmp
          Tmp = -5._RK * RijInv * EPotLocal
          FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                    + (eX * CosThetaj - OXj) * dCosThetaj)
          FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                    + (eY * CosThetaj - OYj) * dCosThetaj)
          FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                    + (eZ * CosThetaj - OZj) * dCosThetaj)
        end if
        EPot(j) = EPot(j) + EPotLocal
        Virial(j) = Virial(j) + Third &
&                     * ( FXij * PXij + FYij * PYij + FZij * PZij )
      end do loop2

    end if


    ! Issue error
    call Error( 'Subroutine TPotQQ_Energy is not implemented' )

  end subroutine TPotQQ_Energy




!==============================================================!
!  Subroutine TPotBond_Construct                               !
!==============================================================!

  subroutine TPotBond_Construct( this, Molecule, j )


    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotBond)              :: this
    type(TMolecule), intent(in) :: Molecule
    integer, intent(in)         :: j

    ! Construct potential

    this%Bond => Molecule%IdfBond(j)
    this%Site1 = this%Bond%SiteId1
    this%Site2 = this%Bond%SiteId2
    this%ForConst = this%Bond%ForConst
    this%R0 = this%Bond%R0

  end subroutine TPotBond_Construct



!==============================================================!
!  Subroutine TPotBond_Destruct                                !
!==============================================================!

  subroutine TPotBond_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotBond) :: this

    ! Destroy potential
    continue

  end subroutine TPotBond_Destruct



!==============================================================!
!  Subroutine TPotBond_Force                                   !
!==============================================================!

  subroutine TPotBond_Force( this, EPot, Virial, EPotIntra, VirialIntra, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotBond)     :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotIntra
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: R, RSquared
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: FXij, FYij, FZij, Fij
    real(RK)          :: EPotLocal, VirialLocal
    real(RK)          :: dR, F0, R0, ForConst

    integer           :: i, j, k, i1, j0, j1
#if MPI_VER > 0
    integer           :: i0
!    integer           :: N1, N2, ji
!    logical           :: EvenN
#endif

#if MPI_VER > 0
!    N1 = this%Site2%NPart
!    N2 = N1 / 2
!    EvenN = mod( N1, 2 ) == 0
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

    ! Assign pointers

     RX1 => this%Bond%RX1
     RY1 => this%Bond%RY1
     RZ1 => this%Bond%RZ1
     RX2 => this%Bond%RX2
     RY2 => this%Bond%RY2
     RZ2 => this%Bond%RZ2

     PX1 => this%Bond%PX1
     PY1 => this%Bond%PY1
     PZ1 => this%Bond%PZ1
     PX2 => this%Bond%PX2
     PY2 => this%Bond%PY2
     PZ2 => this%Bond%PZ2

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
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        FXi = FX1(i)
        FYi = FY1(i)
        FZi = FZ1(i)

!CDIR NODEP

        ! Standard harmonic bond
        ! Energy and forces:
        ! formulae  E = ForConst*(R - R0)**2
        !           F = - 2*ForConst*(R-R0)/R - abs. value

        ! Calculate bond length
        RXij = RXi - RX2(i)
        RYij = RYi - RY2(i)
        RZij = RZi - RZ2(i)
        !
        RXij = (RXij - anint( RXij )) * BoxLength
        RYij = (RYij - anint( RYij )) * BoxLength
        RZij = (RZij - anint( RZij )) * BoxLength
        !
        RSquared=RXij**2+RYij**2+RZij**2
        R=dsqrt(RSquared) ! Bond length

        ! Deviation from equilibrium
        dR=R-R0

        ! Potential parameter
        F0 = dR*ForConst

        ! Energy of the bond
        EPotLocal = EPotLocal + dR*F0

        ! Force (abs. value)
        Fij=-2.0d0*F0/R

        ! Force components
        FXij = Fij * RXij
        FYij = Fij * RYij
        FZij = Fij * RZij

        ! For calculation of virial
        PXij = PXi - PX2(i)
        PYij = PYi - PY2(i)
        PZij = PZi - PZ2(i)
        !
        PXij = (PXij - anint( PXij )) * BoxLength
        PYij = (PYij - anint( PYij )) * BoxLength
        PZij = (PZij - anint( PZij )) * BoxLength

        ! Contribution to virial
        VirialLocal = VirialLocal + PXij * FXij + PYij * FYij + PZij * FZij

         ! New Forces
         FX1(i) = FX1(i) + FXij
         FY1(i) = FY1(i) + FYij
         FZ1(i) = FZ1(i) + FZij
         FX2(i) = FX2(i) - FXij
         FY2(i) = FY2(i) - FYij
         FZ2(i) = FZ2(i) - FZij

       end do

     ! Update potential energy and virial
       EPot = EPot + EPotLocal
       Virial = Virial + Third * VirialLocal

     ! Update Intra potential energy and virial
       EPotIntra = EPotIntra + EPotLocal
       VirialIntra = VirialIntra + Third * VirialLocal

  end subroutine TPotBond_Force



!==============================================================!
!  Subroutine TPotBond_ChemicalPotential                       !
!==============================================================!

  subroutine TPotBond_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotBond) :: this
    real(RK), pointer    :: EPotTest(:)
    real(RK), intent(in) :: BoxLength

    ! Issue error
    call Error( 'Subroutine TPotBond_ChemicalPotential is not implemented' )

  end subroutine TPotBond_ChemicalPotential



!==============================================================!
!  Subroutine TPotBond_Energy                                  !
!==============================================================!

  subroutine TPotBond_Energy( this, np, EPot, Virial, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotBond) :: this
    integer, intent(in)  :: np
    real(RK), pointer    :: EPot(:)
    real(RK), pointer    :: Virial(:)
    real(RK), intent(in) :: BoxLength

    ! Issue error
    call Error( 'Subroutine TPotBond_Energy is not implemented' )

  end subroutine TPotBond_Energy



!==============================================================!
!  Subroutine TPotAngle_Construct                               !
!==============================================================!

  subroutine TPotAngle_Construct( this, Molecule, j )


    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotAngle)              :: this
    type(TMolecule), intent(in) :: Molecule
    integer, intent(in)         :: j

    ! Construct potential

    this%Angle => Molecule%IdfAngle(j)
    this%Site1 = this%Angle%SiteId1
    this%Site2 = this%Angle%SiteId2
    this%Site3 = this%Angle%SiteId3
    this%ForConst = this%Angle%ForConst
    this%Angle0 = this%Angle%Angle0

  end subroutine TPotAngle_Construct



!==============================================================!
!  Subroutine TPotAngle_Destruct                                !
!==============================================================!

  subroutine TPotAngle_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotAngle) :: this

    ! Destroy potential
    continue

  end subroutine TPotAngle_Destruct



!==============================================================!
!  Subroutine TPotAngle_Force                                   !
!==============================================================!

  subroutine TPotAngle_Force( this, EPot, Virial, EPotIntra, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotAngle)     :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotIntra
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:), RX3(:), RY3(:), RZ3(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:), FX3(:), FY3(:), FZ3(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: RXk, RYk, RZk
    real(RK)          :: RijRkj, RijSquared, RkjSquared
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: FXk, FYk, FZk
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: RXkj, RYkj, RZkj
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: FXij, FYij, FZij, Fij
    real(RK)          :: EPotLocal, VirialLocal
    real(RK)          :: ForConst, Angle, Angle0, dAngle, cosa, sina
    real(RK)          :: abc, sab, cab, fab, fbb, faa, fax, fay, faz,  fbx, fby, fbz


    integer           :: i, j, k, i1, j0, j1
#if MPI_VER > 0
    integer           :: i0
!    integer           :: N1, N2, ji
!    logical           :: EvenN
#endif

#if MPI_VER > 0
!    N1 = this%Site2%NPart
!    N2 = N1 / 2
!    EvenN = mod( N1, 2 ) == 0
     i0 = this%Angle%NPart0
     i1 = this%Angle%NPart2
#else
     i1 = this%Angle%NPart
#endif

     ForConst = this%ForConst
     Angle0 = this%Angle0
     EPotLocal   = 0._RK
     VirialLocal = 0._RK

    ! Assign pointers
     RX1 => this%Angle%RX1
     RY1 => this%Angle%RY1
     RZ1 => this%Angle%RZ1
     RX2 => this%Angle%RX2
     RY2 => this%Angle%RY2
     RZ2 => this%Angle%RZ2
     RX3 => this%Angle%RX3
     RY3 => this%Angle%RY3
     RZ3 => this%Angle%RZ3

     FX1 => this%Angle%FX1
     FY1 => this%Angle%FY1
     FZ1 => this%Angle%FZ1
     FX2 => this%Angle%FX2
     FY2 => this%Angle%FY2
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
        RXi = RX1(i)
        RYi = RY1(i)
        RZi = RZ1(i)
        FXi = FX1(i)
        FYi = FY1(i) !           (i)    (k)
        FZi = FZ1(i) !             \    /
        RXk = RX3(i) !            a \  / b
        RYk = RY3(i) !               \/
        RZk = RZ3(i) !               (j)
        FXk = FX3(i)
        FYk = FY3(i)
        FZk = FZ3(i)

!CDIR NODEP
         RXij = RXi - RX2(i)
         RYij = RYi - RY2(i)
         RZij = RZi - RZ2(i)
         RXkj = RXk - RX2(i)
         RYkj = RYk - RY2(i)
         RZkj = RZk - RZ2(i)
         !
         RXij = (RXij - anint( RXij )) * BoxLength
         RYij = (RYij - anint( RYij )) * BoxLength
         RZij = (RZij - anint( RZij )) * BoxLength
         RXkj = (RXkj - anint( RXkj )) * BoxLength
         RYkj = (RYkj - anint( RYkj )) * BoxLength
         RZkj = (RZkj - anint( RZkj )) * BoxLength
         !

         RijSquared=RXij**2+RYij**2+RZij**2
         RkjSquared=RXkj**2+RYkj**2+RZkj**2

         ! Calculate angle
         RijRkj=dsqrt(RijSquared*RkjSquared)
         cosa = (RXij*RXkj+RYij*RYkj+RZij*RZkj)/RijRkj
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

         fax = fab*RXkj-faa*RXij
         fay = fab*RYkj-faa*RYij
         faz = fab*RZkj-faa*RZij

         fbx = fab*RXij-fbb*RXkj
         fby = fab*RYij-fbb*RYkj
         fbz = fab*RZij-fbb*RZkj

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
     EPotIntra = EPotIntra + EPotLocal

  end subroutine TPotAngle_Force



!==============================================================!
!  Subroutine TPotAngle_ChemicalPotential                       !
!==============================================================!

  subroutine TPotAngle_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotAngle) :: this
    real(RK), pointer    :: EPotTest(:)
    real(RK), intent(in) :: BoxLength

    ! Issue error
    call Error( 'Subroutine TPotAngle_ChemicalPotential is not implemented' )

  end subroutine TPotAngle_ChemicalPotential



!==============================================================!
!  Subroutine TPotAngle_Energy                                  !
!==============================================================!

  subroutine TPotAngle_Energy( this, np, EPot, Virial, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotAngle) :: this
    integer, intent(in)  :: np
    real(RK), pointer    :: EPot(:)
    real(RK), pointer    :: Virial(:)
    real(RK), intent(in) :: BoxLength

    ! Issue error
    call Error( 'Subroutine TPotAngle_Energy is not implemented' )

  end subroutine TPotAngle_Energy



!==============================================================!
!  Subroutine TPotDihedral_Construct                               !
!==============================================================!

  subroutine TPotDihedral_Construct( this, Molecule, j )


    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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
    this%ForConst = this%Dihedral%ForConst
    this%gamma = this%Dihedral%gamma
    this%multi = this%Dihedral%multi
    this%ScaleLJ14=this%Dihedral%ScaleLJ14
    this%ScaleEl14=this%Dihedral%ScaleEl14
    this%Sigma1=this%Dihedral%Sigma1
    this%Sigma4=this%Dihedral%Sigma4
    this%Epsilon1=this%Dihedral%Epsilon1
    this%Epsilon4=this%Dihedral%Epsilon4



  end subroutine TPotDihedral_Construct



!==============================================================!
!  Subroutine TPotDihedral_Destruct                                !
!==============================================================!

  subroutine TPotDihedral_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDihedral) :: this

    ! Destroy potential
    continue

  end subroutine TPotDihedral_Destruct



!==============================================================!
!  Subroutine TPotDihedral_Force                                   !
!==============================================================!

  subroutine TPotDihedral_Force( this, EPot, Virial, EPotIntra, VirialIntra,  BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDihedral)     :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: EPotIntra
    real(RK), intent(in out) :: VirialIntra
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:), RX3(:), RY3(:), RZ3(:), RX4(:), RY4(:), RZ4(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:), FX3(:), FY3(:), FZ3(:), FX4(:), FY4(:), FZ4(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX4(:), PY4(:), PZ4(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: RXj, RYj, RZj
    real(RK)          :: RXk, RYk, RZk
    real(RK)          :: RXl, RYl, RZl
    real(RK)          :: FXi, FYi, FZi
    real(RK)          :: FXk, FYk, FZk
    real(RK)          :: FXl, FYl, FZl
    real(RK)          :: FXj, FYj, FZj
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: PXl, PYl, PZl
    real(RK)          :: RilSquaredInv, SigmaSquared, Ril6Inv, Sigma,Epsilon, Epsilon48, BoxLengthInv
    real(RK)          :: RXil, RYil, RZil, PXil, PYil, PZil, FXil, FYil, FZil, Fil
    real(RK)          :: EPotLocal, VirialLocal
    real(RK)          :: num, den, de1, ax, ay, az, bx, by, bz, cx, cy, cz
    real(RK)          :: ab, bc, ac, aa, bb, cc, axb, bxc, co, si, signum, arg, cos1, cosn, earg
    real(RK)          :: deri,dnum,dden,ffi,ffj,ffk,ffl, ForConst, gamma

    integer           :: i, j, k, i1, j0, j1, multi

#if MPI_VER > 0
     integer           :: i0
!    integer           :: N1, N2, ji
!    logical           :: EvenN
#endif

#if MPI_VER > 0
!    N1 = this%Site2%NPart
!    N2 = N1 / 2
!    EvenN = mod( N1, 2 ) == 0
    i0 = this%Dihedral%NPart0
    i1 = this%Dihedral%NPart2
#else
     i1 = this%Dihedral%NPart
#endif

    gamma = this%gamma
    ForConst = this%ForConst
    multi =this%multi
    Sigma = .5_RK * (this%Sigma1 + this%Sigma4)
    Epsilon = sqrt(this%Epsilon1 * this%Epsilon4)
    BoxLengthInv = 1._RK / BoxLength
    SigmaSquared = (Sigma * BoxLengthInv)**2
    Epsilon48 = 48._RK * Epsilon * BoxLengthInv &
&      / SigmaSquared

    EPotLocal   = 0._RK
    VirialLocal = 0._RK

    ! Assign pointers

     RX1 => this%Dihedral%RX1
     RY1 => this%Dihedral%RY1
     RZ1 => this%Dihedral%RZ1
     RX2 => this%Dihedral%RX2 !                  (i)            (l)
     RY2 => this%Dihedral%RY2 !                    \            /
     RZ2 => this%Dihedral%RZ2 !                  a  \          / c
     RX3 => this%Dihedral%RX3 !                      (j)-----(k)
     RY3 => this%Dihedral%RY3 !                            b
     RZ3 => this%Dihedral%RZ3
     RX4 => this%Dihedral%RX4
     RY4 => this%Dihedral%RY4
     RZ4 => this%Dihedral%RZ4

     FX1 => this%Dihedral%FX1
     FY1 => this%Dihedral%FY1
     FZ1 => this%Dihedral%FZ1
     FX2 => this%Dihedral%FX2
     FY2 => this%Dihedral%FY2
     FZ2 => this%Dihedral%FZ2
     FX3 => this%Dihedral%FX3
     FY3 => this%Dihedral%FY3
     FZ3 => this%Dihedral%FZ3
     FX4 => this%Dihedral%FX4
     FY4 => this%Dihedral%FY4
     FZ4 => this%Dihedral%FZ4

!     if (LJEl14) then
!       PX1 => this%Dihedral%PX1  ! For 1-4 intramolecular interaction - calculation of virial
!       PY1 => this%Dihedral%PY1
!       PZ1 => this%Dihedral%PZ1
!       PX4 => this%Dihedral%PX4
!       PY4 => this%Dihedral%PY4
!       PZ4 => this%Dihedral%PZ4
!     end if

      ! Loop over molecules
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
        RXj = RX2(i)
        RYj = RY2(i)
        RZj = RZ2(i)
        FXj = FX2(i)
        FYj = FY2(i)
        FZj = FZ2(i)
        RXk = RX3(i)
        RYk = RY3(i)
        RZk = RZ3(i)
        FXk = FX3(i)
        FYk = FY3(i)
        FZk = FZ3(i)
        RXl = RX4(i)
        RYl = RY4(i)
        RZl = RZ4(i)
        FXl = FX4(i)
        FYl = FY4(i)
        FZl = FZ4(i)

!     if (LJEl14) then
!        PXi = PX1(i)
!        PYi = PY1(i)
!        PZi = PZ1(i)
!        PXl = PX4(i)
!        PYl = PY4(i)
!        PZl = PZ4(i)
!     end if

!CDIR NODEP

        if (multi .eq. 0) then
           EPotLocal = EPotLocal+ForConst*2._RK
        else
          ! Calculate vectors IJ, JK, KL
          ax = (RXj - RXi)
          ay = (RYj - RYi)
          az = (RZj - RZi)
          bx = (RXk - RXj)
          by = (RYk - RYj)
          bz = (RZk - RZj)
          cx = (RXl - RXk)
          cy = (RYl - RYk)
          cz = (RZl - RZk)
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
            si = sin(arg)
            if( abs(si) .lt. 1E-10_RK ) si = sign( 1E-10_RK, si )


            if (multi > 0) then
               ! Normal Amber-type torsion angle
               earg= multi*arg-gamma

               ! Energy and forces:
               ! formulae  E = ForConst*( 1 + cos(earg) )
               !           F = ForConst*n*sin(earg)

                EPotLocal  = EPotLocal + ForConst*(1.d0+dcos(earg))
                deri= -ForConst*dble(multi)*dsin(earg)

             else ! Improper dihedral angle
               earg= arg-gamma

               ! Energy and forces:
               ! formulae  E = ForConst*earg**2
               !           F = -2*ForConst*earg
                EPotLocal  = EPotLocal + ForConst*earg**2
                deri= 2.d0*ForConst*earg
             end if

             ! Calculate Forces
             axb = axb/den*co
             bxc = bxc/den*co
             de1 = deri/den/si

! X components
            dnum = cx*bb - bx*bc
            dden = ( ab*bx - ax*bb )*bxc
            FFI = (dnum - dden) * de1
            dnum = ((bx-ax)*bc - ab*cx ) + (2.0*ac*bx - cx*bb)
            dden = axb*(bc*cx-bx*cc) + (ax*bb-aa*bx-ab*(bx-ax))*bxc
            FFJ = (dnum - dden) * de1
            dnum = ab*bx - ax*bb
            dden = axb*( bb*cx - bc*bx )
            FFL = (dnum - dden) * de1
            FFK = -(ffi+ffj+ffl)

! Forces
            FX1(i) = FXi+ffi
            FX2(i) = FXj+ffj
            FX3(i) = FXk+ffk
            FX4(i) = FXl+ffl

! Y components
            dnum = cy*bb - by*bc
            dden = ( ab*by - ay*bb )*bxc
            FFI = (dnum - dden) * de1
            dnum = ((by-ay)*bc - ab*cy ) + (2.0*ac*by - cy*bb)
            dden = axb*(bc*cy-by*cc) + (ay*bb-aa*by-ab*(by-ay))*bxc
            FFJ = (dnum - dden) * de1
            dnum = ab*by - ay*bb
            dden = axb*( bb*cy - bc*by )
            FFL = (dnum - dden) * de1
            FFK = -(ffi+ffj+ffl)

! Forces
            FY1(i) = FYi+ffi
            FY2(i) = FYj+ffj
            FY3(i) = FYk+ffk
            FY4(i) = FYl+ffl

! Z components
            dnum = cz*bb - bz*bc
            dden = ( ab*bz - az*bb )*bxc
            FFI = (dnum - dden) * de1
            dnum = ((bz-az)*bc - ab*cz ) + (2.0*ac*bz - cz*bb)
            dden = axb*(bc*cz-bz*cc) + (az*bb-aa*bz-ab*(bz-az))*bxc
            FFJ = (dnum - dden) * de1
            dnum = ab*bz - az*bb
            dden = axb*( bb*cz - bc*bz )
            FFL = (dnum - dden) * de1
            FFK = -(ffi+ffj+ffl)

! Forces
            FZ1(i) = FZi+ffi
            FZ2(i) = FZj+ffj
            FZ3(i) = FZk+ffk
            FZ4(i) = FZl+ffl

          endif ! den>0
        endif ! multi/=0

!        ! Calculate Intramolecular 1-4 interactions if need
!         if ( LJEl14 .and. (multi > 0)) then
!           !Calculate LJ and El interactions, don't forget corrections!!
!            RXil = RXi - RXl
!            RYil = RYi - RYl
!            RZil = RZi - RZl
!            PXil = PXi - PXl ! the same molecule, but different units
!            PYil = PYi - PYl
!            PZil = PZi - PZl
!            RXil = RXil - anint( RXil )
!            RYil = RYil - anint( RYil )
!            RZil = RZil - anint( RZil )
!            PXil = PXil - anint( PXil )
!            PYil = PYil - anint( PYil )
!            PZil = PZil - anint( PZil )
!            RilSquaredInv = SigmaSquared / ( RXil**2 + RYil**2 + RZil**2 )
!            Ril6Inv = RilSquaredInv**3
!            EPotLocal = EPotLocal + this%ScaleLJ14*Ril6Inv * (Ril6Inv - 1._RK)/multi
!            Fil = Epsilon48 * Ril6Inv * (Ril6Inv - .5_RK) * RilSquaredInv/multi
!            FXil = Fil * RXil * this%ScaleLJ14
!            FYil = Fil * RYil * this%ScaleLJ14
!            FZil = Fil * RZil * this%ScaleLJ14
!            VirialLocal = VirialLocal + PXil * FXil + PYil * FYil + PZil * FZil
!            FX1(i) = FXi + FXil
!            FY1(i) = FYi + FYil
!            FZ1(i) = FZi + FZil
!            FX4(i) = FXl - FXil
!            FY4(i) = FYl - FYil
!            FZ4(i) = FZl - FZil
!          end if
    enddo


    ! Update potential energy, no contribution to virial!
     EPot = EPot +  EPotLocal

    ! Update Intra potential energy
     EPotIntra = EPotIntra +  EPotLocal

!    ! Contribution to virial from 1-4 interactions
!    Virial = Virial+VirialLocal
!    VirialIntra = VirialIntra+VirialLocal
!    print *, 'VirialIntra from Dihedral=', VirialLocal


  end subroutine TPotDihedral_Force

!==============================================================!
!  Subroutine TPotDihedral_ChemicalPotential                       !
!==============================================================!

  subroutine TPotDihedral_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDihedral) :: this
    real(RK), pointer    :: EPotTest(:)
    real(RK), intent(in) :: BoxLength

    ! Issue error
    call Error( 'Subroutine TPotDihedral_ChemicalPotential is not implemented' )

  end subroutine TPotDihedral_ChemicalPotential



!==============================================================!
!  Subroutine TPotDihedral_Energy                                  !
!==============================================================!

  subroutine TPotDihedral_Energy( this, np, EPot, Virial, BoxLength )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TPotDihedral) :: this
    integer, intent(in)  :: np
    real(RK), pointer    :: EPot(:)
    real(RK), pointer    :: Virial(:)
    real(RK), intent(in) :: BoxLength

    ! Issue error
    call Error( 'Subroutine TPotDihedral_Energy is not implemented' )

  end subroutine TPotDihedral_Energy




end module ms2_potential
