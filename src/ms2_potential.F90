!==============================================================!
!  MOLECULAR SIMULATION PROGRAM MS2 Version 1.1 v12            !
!  (c) 2001 by Sergey Lishchuk, ITT                            !
!  (c) 2007 by Bernhard Eckl, ITT                              !
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

!#ifdev MPI_VER > 0
!  use mpi
!#endif

  use ms2_molecule
  use ms2_site



!==============================================================!
!  Type TPotLJ126LJ126                                         !
!==============================================================!

  type TPotLJ126LJ126

    type(TSiteLJ126), pointer :: Site1, Site2
    real(RK)                  :: Sigma, Epsilon
    real(RK)                  :: RCutoffSquared, RCutoffSquaredScaled
    real(RK)                  :: EPotCorr, VirialCorr, EPotTestCorr
    logical                   :: SameComponent
    real(RK)                  :: SigmaSquared
    real(RK)                  :: Epsilon4, Epsilon48
    real(RK)                  :: BoxlengthInv, BoxLengthThird
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
    real(RK)                   :: Epsilon
    real(RK)                   :: RShieldSquared
    real(RK)                   :: RCutoffSquared
    logical                    :: SameComponent
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
    real(RK)                   :: Epsilon
    real(RK)                   :: RShieldSquared
    real(RK)                   :: RCutoffSquared
    logical                    :: SameComponent
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
    real(RK)                       :: Epsilon
    real(RK)                       :: RShieldSquared
    real(RK)                       :: RCutoffSquared
    logical                        :: SameComponent
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
    real(RK)                   :: Epsilon
    real(RK)                   :: RShieldSquared
    real(RK)                   :: RCutoffSquared
    logical                    :: SameComponent
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
    real(RK)                   :: Epsilon
    real(RK)                   :: RCutoffSquared
    real(RK)                   :: RShieldSquared
    real(RK)                   :: RFConstant
    logical                    :: SameComponent
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
    real(RK)                       :: Epsilon
    real(RK)                       :: RCutoffSquared
    real(RK)                       :: RShieldSquared
    logical                        :: SameComponent
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
    real(RK)                       :: Epsilon
    real(RK)                       :: RShieldSquared
    real(RK)                       :: RCutoffSquared
    logical                        :: SameComponent
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
    real(RK)                       :: Epsilon
    real(RK)                       :: RCutoffSquared
    real(RK)                       :: RShieldSquared
    logical                        :: SameComponent
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
    real(RK)                       :: Epsilon
    real(RK)                       :: RCutoffSquared
    real(RK)                       :: RShieldSquared
    logical                        :: SameComponent
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



contains



!==============================================================!
!  Subroutine TPotLJLJ_Construct                               !
!==============================================================!

  subroutine TPotLJLJ_Construct( this, i1, i2, j1, j2, &
&                                Molecule1, Molecule2, &
&                                RCutoff, ScaleSigma, ScaleEpsilon )

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

    ! Construct potential
    this%Site1 => Molecule1%SiteLJ126(j1)
    this%Site2 => Molecule2%SiteLJ126(j2)
    this%SameComponent = i1 == i2
    this%Sigma = .5_RK * (this%Site1%sig + this%Site2%sig)
    this%Epsilon = sqrt(this%Site1%eps * this%Site2%eps)
    if( .not. this%SameComponent ) then
      this%Sigma = this%Sigma * ScaleSigma
      this%Epsilon = this%Epsilon * ScaleEpsilon
    end if
    this%Epsilon4 = 4._RK * this%Epsilon

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

    ! Declare arguments
    type(TPotLJ126LJ126) :: this

    ! Destroy potential
    continue

  end subroutine TPotLJLJ_Destruct



!==============================================================!
!  Subroutine TPotLJLJ_Force                                   !
!==============================================================!

  subroutine TPotLJLJ_Force( this, EPot, Virial, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotLJ126LJ126)     :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
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
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
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
    SigmaSquared = this%SigmaSquared
    Epsilon4 = this%Epsilon4
    Epsilon48 = this%Epsilon48
    RCutoffSquared = this%RCutoffSquaredScaled
    EPotLocal   = 0._RK
    VirialLocal = 0._RK

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

    if( CutoffMode .eq. CenterofMass ) then

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
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
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
          PXij = PXij - anint( PXij )
          PYij = PYij - anint( PYij )
          PZij = PZij - anint( PZij )
          RijSquaredInv = SigmaSquared / ( RXij**2 + RYij**2 + RZij**2 )
          Rij6Inv = RijSquaredInv**3
          EPotLocal = EPotLocal + Rij6Inv * (Rij6Inv - 1._RK)
          Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
          FXij = Fij * RXij
          FYij = Fij * RYij
          FZij = Fij * RZij
          VirialLocal = VirialLocal + PXij * FXij + PYij * FYij + PZij * FZij
          FXi = FXi + FXij
          FYi = FYi + FYij
          FZi = FZi + FZij
          FX2(j) = FX2(j) - FXij
          FY2(j) = FY2(j) - FYij
          FZ2(j) = FZ2(j) - FZij
        end do loop1
        FX1(i) = FXi
        FY1(i) = FYi
        FZ1(i) = FZi
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
          Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
          FXij = Fij * RXij
          FYij = Fij * RYij
          FZij = Fij * RZij
          VirialLocal = VirialLocal + PXij * FXij + PYij * FYij + PZij * FZij
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

    end if

    ! Update potential energy and virial
    EPot = EPot + Epsilon4 * EPotLocal
    Virial = Virial + Third * VirialLocal * BoxLength

  end subroutine TPotLJLJ_Force



!==============================================================!
!  Subroutine TPotLJLJ_ChemicalPotential                       !
!==============================================================!

  subroutine TPotLJLJ_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

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

    ! Declare arguments
    type(TPotChargeCharge)      :: this
    integer, intent(in)         :: i1, i2, j1, j2
    type(TMolecule), intent(in) :: Molecule1, Molecule2
    real(RK), intent(in)        :: RCutoff

    ! Construct potential
    this%Site1 => Molecule1%SiteCharge(j1)
    this%Site2 => Molecule2%SiteCharge(j2)
    this%SameComponent = i1 == i2
    this%Epsilon = this%Site1%e * this%Site2%e
    this%RCutoffSquared = RCutoff**2
    this%RShieldSquared = .25_RK * ( this%Site1%shield + this%Site2%shield )**2

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

  subroutine TPotCC_Force( this, EPot, Virial, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
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
    integer           :: i, j, k, i1
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
    Epsilon = this%Epsilon
    EPotLocal = 0._RK
    VirialLocal = 0._RK

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
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)
!CDIR NODEP
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
&                       * (eX * PXij + eY * PYij + eZ * PZij)
        FXij = EPotLocal1 * RijInv * eX
        FYij = EPotLocal1 * RijInv * eY
        FZij = EPotLocal1 * RijInv * eZ
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        FX2(j) = FX2(j) - FXij
        FY2(j) = FY2(j) - FYij
        FZ2(j) = FZ2(j) - FZij
      end do loop1
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
    end do

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal

  end subroutine TPotCC_Force



!==============================================================!
!  Subroutine TPotCC_ChemicalPotential                         !
!==============================================================!

  subroutine TPotCC_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

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

  subroutine TPotCD_Force( this, EPot, Virial, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeDipole)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
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
    integer           :: i, j, k, i1
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
    Epsilon = this%Epsilon
    EPotLocal = 0._RK
    VirialLocal = 0._RK

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
        FXij = Epsilon2 * ( CosTheta3 * eX - OXj )                              ! F2 bei Price
        FYij = Epsilon2 * ( CosTheta3 * eY - OYj )  
        FZij = Epsilon2 * ( CosTheta3 * eZ - OZj )
        VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij     ! F2*R_COM_Price; stimmt so
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        FX2(j) = FX2(j) - FXij
        FY2(j) = FY2(j) - FYij
        FZ2(j) = FZ2(j) - FZij
        TX2(j) = TX2(j) - Epsilon1 * eX                                         ! Uebereinstimmung mit Price
        TY2(j) = TY2(j) - Epsilon1 * eY
        TZ2(j) = TZ2(j) - Epsilon1 * eZ
      end do loop1
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
    end do

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal

  end subroutine TPotCD_Force


!==============================================================!
!  Subroutine TPotCD_ChemicalPotential                         !
!==============================================================!

  subroutine TPotCD_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

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

  subroutine TPotCQ_Force( this, EPot, Virial, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotChargeQuadrupole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
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
    integer           :: i, j, k, i1
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
    Epsilon = this%Epsilon
    EPotLocal = 0._RK
    VirialLocal = 0._RK

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
        RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
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
        VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij     ! Vorzeichen richtig so
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        FX2(j) = FX2(j) - FXij
        FY2(j) = FY2(j) - FYij
        FZ2(j) = FZ2(j) - FZij
        TX2(j) = TX2(j) - Epsilon1 * CosTheta2 * eX 
        TY2(j) = TY2(j) - Epsilon1 * CosTheta2 * eY
        TZ2(j) = TZ2(j) - Epsilon1 * CosTheta2 * eZ
      end do loop1
      FX1(i) = FXi
      FY1(i) = FYi
      FZ1(i) = FZi
    end do

    ! Update potential energy and virial
    EPot = EPot + EPotLocal
    Virial = Virial + Third * VirialLocal

  end subroutine TPotCQ_Force



!==============================================================!
!  Subroutine TPotCQ_ChemicalPotential                         !
!==============================================================!

  subroutine TPotCQ_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

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

  subroutine TPotDC_Force( this, EPot, Virial, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleCharge)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
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
    integer           :: i, j, k, i1
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
    Epsilon = this%Epsilon
    EPotLocal = 0._RK
    VirialLocal = 0._RK

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
        FXij = Epsilon2 * ( OXi - CosTheta3 * eX )                              ! F1 bei Price
        FYij = Epsilon2 * ( OYi - CosTheta3 * eY )
        FZij = Epsilon2 * ( OZi - CosTheta3 * eZ )
        VirialLocal = VirialLocal + FXij * PXij + FYij * PYij + FZij * PZij     ! F1*(-R_COM_Price); stimmt so
        FXi    = FXi    + FXij
        FYi    = FYi    + FYij
        FZi    = FZi    + FZij
        TXi = TXi + Epsilon1 * eX                                       ! Uebereinstimmumg mit Price; Rest bei Atom2Mol in Component
        TYi = TYi + Epsilon1 * eY                                       ! Reaktionsfeldbeitrag in Interaction
        TZi = TZi + Epsilon1 * eZ
        FX2(j) = FX2(j) - FXij
        FY2(j) = FY2(j) - FYij
        FZ2(j) = FZ2(j) - FZij
      end do loop1
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

  end subroutine TPotDC_Force



!==============================================================!
!  Subroutine TPotDC_ChemicalPotential                         !
!==============================================================!

  subroutine TPotDC_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

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
    this%RFConstant = this%Epsilon / RCutoff**3 &
&     * (RFEpsilon - 1._RK) / (2._RK * RFEpsilon + 1._RK)

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

  subroutine TPotDD_Force( this, EPot, Virial, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleDipole)   :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
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
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
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
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    RFConstant2 = 2._RK * this%RFConstant
    EPotLocal = 0._RK
    VirialLocal = 0._RK

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
!           EPotLocal = EPotLocal +  Rij3Inv * Tmp - RFConstant2 * CosGammaij
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
          TXi    = TXi    + Rij3Inv * (eX * CosThetaj3 - OXj)
          TYi    = TYi    + Rij3Inv * (eY * CosThetaj3 - OYj)
          TZi    = TZi    + Rij3Inv * (eZ * CosThetaj3 - OZj)
          TX2(j) = TX2(j) + Rij3Inv * (eX * CosThetai3 - OXi)
          TY2(j) = TY2(j) + Rij3Inv * (eY * CosThetai3 - OYi)
          TZ2(j) = TZ2(j) + Rij3Inv * (eZ * CosThetai3 - OZi)
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
        end do loop1
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

  end subroutine TPotDD_Force



!==============================================================!
!  Subroutine TPotDD_ChemicalPotential                         !
!==============================================================!

  subroutine TPotDD_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

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

  subroutine TPotDQ_Force( this, EPot, Virial, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotDipoleQuadrupole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
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
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
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
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    EPotLocal   = 0._RK
    VirialLocal = 0._RK

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
        end do loop1
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

  end subroutine TPotDQ_Force



!==============================================================!
!  Subroutine TPotDQ_ChemicalPotential                         !
!==============================================================!

  subroutine TPotDQ_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

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

  subroutine TPotQC_Force( this, EPot, Virial, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleCharge) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
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
    integer           :: i, j, k, i1
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
    Epsilon = this%Epsilon
    EPotLocal = 0._RK
    VirialLocal = 0._RK

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
        RijSquaredInv = 1._RK / ( RXij**2 + RYij**2 + RZij**2 )
        RijInv = sqrt( RijSquaredInv )
        eX = - RXij * RijInv                                                    ! Normierter Abstandsvektor nach Price
        eY = - RYij * RijInv
        eZ = - RZij * RijInv
        CosTheta  = OXi * ex + OYi * eY + OZi * eZ       ! Scalarprodukt normierter Abstandsvektor mit Orientierungsvektor Quadrupol
        Epsilon1 = Epsilon * RijSquaredInv * RijInv
        EPotLocal  = EPotLocal + Epsilon1 * ( CosTheta * CosTheta - Third )
        CosTheta2 = 2._RK * CosTheta
        CosAux = 5._RK *  CosTheta * CosTheta - 1._RK
        Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
        FXij = Epsilon2 * ( CosAux * eX - CosTheta2 * OXi )                     ! Kraft auf die Punktladung, sprich F2
        FYij = Epsilon2 * ( CosAux * eY - CosTheta2 * OYi )
        FZij = Epsilon2 * ( CosAux * eZ - CosTheta2 * OZi )
        VirialLocal = VirialLocal - FXij * PXij - FYij * PYij - FZij * PZij     ! Vorzeichen richtig
        FXi    = FXi    - FXij
        FYi    = FYi    - FYij
        FZi    = FZi    - FZij
        FX2(j) = FX2(j) + FXij
        FY2(j) = FY2(j) + FYij
        FZ2(j) = FZ2(j) + FZij
        TXi    = TXi - Epsilon1*CosTheta2*eX                        ! Drehmomentanteil auf Quadrupol wegen Punktladung. Kreuzprodukt
        TYi    = TYi - Epsilon1*CosTheta2*eY                        ! in Atom2Mol von Component
        TZi    = TZi - Epsilon1*CosTheta2*eZ
      end do loop1
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

  end subroutine TPotQC_Force



!==============================================================!
!  Subroutine TPotQC_ChemicalPotential                         !
!==============================================================!

  subroutine TPotQC_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

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
        CosTheta  = OXi * ex + OYi * eY + OZi * eZ       ! Scalarprodukt normierter Abstandsvektor mit Orientierungsvektor Quadrupol
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

  subroutine TPotQD_Force( this, EPot, Virial, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleDipole) :: this
    real(RK), intent(in out)   :: EPot
    real(RK), intent(in out)   :: Virial
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
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
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
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    EPotLocal   = 0._RK
    VirialLocal = 0._RK

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
        end do loop1
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

  end subroutine TPotQD_Force



!==============================================================!
!  Subroutine TPotQD_ChemicalPotential                         !
!==============================================================!

  subroutine TPotQD_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

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

  subroutine TPotQQ_Force( this, EPot, Virial, BoxLength )

    implicit none

    ! Declare arguments
    type(TPotQuadrupoleQuadrupole) :: this
    real(RK), intent(in out)       :: EPot
    real(RK), intent(in out)       :: Virial
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
    logical           :: SameComponent
    integer           :: i, j, k, i1, j0, j1
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
    Epsilon = this%Epsilon
    RCutoffSquared = this%RCutoffSquared
    EPotLocal   = 0._RK
    VirialLocal = 0._RK

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
        end do loop1
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

  end subroutine TPotQQ_Force



!==============================================================!
!  Subroutine TPotQQ_ChemicalPotential                         !
!==============================================================!

  subroutine TPotQQ_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

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

  end subroutine TPotQQ_Energy



end module ms2_potential
