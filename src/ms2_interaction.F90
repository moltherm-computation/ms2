!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 1.0                !
!  (c) 2011 by TU Kaiserslautern                               !
!      P.O. Box 67653                                          !
!      67653 Kaiserslautern                                    !
!==============================================================!
!  Module ms2_interaction                                      !
!  Contains TInteraction object                                !
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
!DEC$ MESSAGE:'Compiling ms2_interaction.F90...'
#endif

module ms2_interaction

  use ms2_potential
  use ms2_component
  use ms2_site



!==============================================================!
!  Type TInteraction                                           !
!==============================================================!

  type TInteraction

    ! Site-site potentials
    type(TPotLJ126LJ126), pointer           :: PotLJ126LJ126(:, :)
    type(TPotChargeCharge), pointer         :: PotChargeCharge(:, :)
    type(TPotChargeDipole), pointer         :: PotChargeDipole(:, :)
    type(TPotChargeQuadrupole), pointer     :: PotChargeQuadrupole(:, :)
    type(TPotDipoleCharge), pointer         :: PotDipoleCharge(:, :)
    type(TPotDipoleDipole), pointer         :: PotDipoleDipole(:, :)
    type(TPotDipoleQuadrupole), pointer     :: PotDipoleQuadrupole(:, :)
    type(TPotQuadrupoleCharge), pointer     :: PotQuadrupoleCharge(:, :)
    type(TPotQuadrupoleDipole), pointer     :: PotQuadrupoleDipole(:, :)
    type(TPotQuadrupoleQuadrupole), pointer :: PotQuadrupoleQuadrupole(:, :)

    ! Potential energy
    real(RK), pointer :: EPot(:, :), EPot1(:), EPotNew(:, :)

    ! Mayer f-function for second virial coefficient
    real(RK), pointer :: MayerFFunction(:), IntFFunction(:)
    real(RK), pointer :: MayerFFunction1(:), IntFFunction1(:)
    real(RK), pointer :: MayerFFunction2(:), IntFFunction2(:)

    ! Virial
    real(RK), pointer :: Virial(:, :), Virial1(:), VirialNew(:, :)
    logical           :: OptPressure

    real(RK), pointer :: d2EpotdV2(:, :), d2EpotdV21(:), d2EpotdV2New(:, :)

    ! Arrays for center of mass cutoff
    integer, pointer :: NInCutoff(:), CutoffPartner(:, :)!, RDFSum (:)

    ! Center of mass positions
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)

    ! Total dipole moments of molecules for reaction field
    real(RK), pointer :: MueX1(:), MueY1(:), MueZ1(:)
    real(RK), pointer :: MueX2(:), MueY2(:), MueZ2(:)

    ! Torques from reaction field
    real(RK), pointer :: tRFX1(:), tRFY1(:), tRFZ1(:)
    real(RK), pointer :: tRFX2(:), tRFY2(:), tRFZ2(:)

    ! Center of mass positions of test particles
    real(RK), pointer :: PX1Test(:), PY1Test(:), PZ1Test(:)

    ! Total dipole moments of test particles for reaction field
    real(RK), pointer :: MueX1Test(:), MueY1Test(:), MueZ1Test(:)

    ! Maximum number of particles per component
    integer :: NPartMax

    ! Numbers of particles
    integer, pointer :: NPart1, NPart2
#if MPI_VER > 0
    integer, pointer :: NPart10, NPart12
    integer, pointer :: NPart20, NPart22
#endif

    ! Numbers of test particles in component 1
    integer :: NTest1

    ! Numbers of sites
    integer :: N1LJ126, N2LJ126
    integer :: N1Charge, N2Charge
    integer :: N1Dipole, N2Dipole
    integer :: N1Quadrupole, N2Quadrupole

    ! Squared cutoff radius
    real(RK) :: RCutoffSquared, RCutoffSquaredScaled

    ! Cutoff correction to LJ-interaction
    real(RK) :: EPotCorrLJ

    ! Flag for reaction field
    logical :: ReactionField

    ! (2*eps-1)/(2*eps+1)
    real(RK) :: RFConst2, RFConst3

    ! Same component
    logical :: SameComponent

    ! Ewald Summation
    real(RK) :: Kappa
    real(RK) :: DebyeLen, RFConstant
    real(RK) :: lad1,lad2

#ifdef ABL
    real(RK),pointer :: AblS(:)
    real(RK),pointer :: AblE(:)
    real(RK),pointer :: AblPS(:,:)
    real(RK),pointer :: AblPE(:,:)
#endif

  end type TInteraction

  interface Construct
    module procedure TInteraction_Construct
  end interface

  interface Destruct
    module procedure TInteraction_Destruct
  end interface

  interface Allocate
    module procedure TInteraction_Allocate
  end interface

  interface DeallocateEPot
    module procedure TInteraction_DeallocateEPot
  end interface

  interface Deallocate
    module procedure TInteraction_Deallocate
  end interface

  interface Force
    module procedure TInteraction_Force
  end interface
  
  interface GET_RDF
   module procedure TInteraction_RDF
  end interface

  interface ChemicalPotential
    module procedure TInteraction_ChemicalPotential
  end interface

  interface Energy
    module procedure TInteraction_Energy
  end interface

  interface UpdateBoxLength
    module procedure TInteraction_UpdateBoxLength
  end interface

  interface CalcCutoffPartners
    module procedure TInteraction_CalcPartners
    module procedure TInteraction_CalcPartners1
  end interface
  
  interface CalcCutoffPartnersRDF
    module procedure TInteraction_CalcPartnersRDF
  end interface  

  interface CalcCutoffPartnersTest
    module procedure TInteraction_CalcPartnersTest
  end interface


contains


!==============================================================!
!  Subroutine TInteraction_Construct                           !
!==============================================================!

  subroutine TInteraction_Construct( this, i1, i2, &
&                                    Component1, Component2, &
&                                    RCutoffLJ126LJ126, &
&                                    RCutoffDipoleDipole, &
&                                    RCutoffDipoleQuadrupole, &
&                                    RCutoffQuadrupoleQuadrupole, &
&                                    ScaleSigma, ScaleEpsilon, &
&                                    RFEpsilon )

    implicit none

    ! Declare arguments
    type(TInteraction)           :: this
    integer, intent(in)          :: i1, i2
    type(TComponent), intent(in) :: Component1, Component2
    real(RK), intent(in)         :: RCutoffLJ126LJ126
    real(RK), intent(in)         :: RCutoffDipoleDipole
    real(RK), intent(in)         :: RCutoffDipoleQuadrupole
    real(RK), intent(in)         :: RCutoffQuadrupoleQuadrupole
    real(RK), intent(in)         :: ScaleSigma, ScaleEpsilon
    real(RK), intent(in)         :: RFEpsilon

    ! Declare local variables
    integer :: j1, j2
    integer :: stat
    real    :: fac

    ! RFConstant2
    if (LongRange .eq. RField) then
      this%RFConst2 = -2._RK / RCutoffDipoleDipole**3 * (RFEpsilon - 1._RK) / (2._RK * RFEpsilon + 1._RK)

    else
      fac = this%DebyeLen*RCutoffDipoleDipole
      this%RFConst2 = - 2._RK / RCutoffDipoleDipole**3 &
&                     * ( (RFEpsilon - 1._RK)*(1._RK+fac)+ 0.5*RFEpsilon*(fac)**2 )   &
&                     / ( (2._RK * RFEpsilon+1._RK)*(1._RK+fac) + RFEpsilon*(fac)**2 )

      this%RFConst3 = -3._RK / RCutoffDipoleDipole * RFEpsilon*(1._RK+fac+0.5*(fac)**2) &
&                     / ( (2._RK * RFEpsilon + 1._RK)*(1+fac) + RFEpsilon*(fac)**2 )
      this%RFConstant=RCutoffDipoleDipole
    end if

    ! Set SameComponent flag
    this%SameComponent = i1 == i2

    ! Set number of particles
    this%NPart1 => Component1%NPart
    this%NPart2 => Component2%NPart
    this%NPartMax = max( Component1%NPartMax, Component2%NPartMax )
#if MPI_VER > 0
    this%NPart10 => Component1%NPart0
    this%NPart12 => Component1%NPart2
    this%NPart20 => Component2%NPart0
    this%NPart22 => Component2%NPart2
#endif
    this%NTest1 = Component1%NTest

    ! Set number of sites
    this%N1LJ126 = Component1%Molecule%NLJ126
    this%N2LJ126 = Component2%Molecule%NLJ126
    this%N1Charge = Component1%Molecule%NCharge
    this%N2Charge = Component2%Molecule%NCharge
    this%N1Dipole = Component1%Molecule%NDipole
    this%N2Dipole = Component2%Molecule%NDipole
    this%N1Quadrupole = Component1%Molecule%NQuadrupole
    this%N2Quadrupole = Component2%Molecule%NQuadrupole

    ! Set center of mass positions
    this%PX1 => Component1%P0(:, 1)
    this%PY1 => Component1%P0(:, 2)
    this%PZ1 => Component1%P0(:, 3)
    this%PX2 => Component2%P0(:, 1)
    this%PY2 => Component2%P0(:, 2)
    this%PZ2 => Component2%P0(:, 3)
!    write(*,*) "after PX1 ",LOC(this%PX1)
!    write(*,*) "after PY1 ",LOC(this%PY1)
!    write(*,*) "after PZ1 ",LOC(this%PZ1)
!    write(*,*) "after PX2 ",LOC(this%PX2)
!    write(*,*) "after PY2 ",LOC(this%PY2)
!    write(*,*) "after PZ2 ",LOC(this%PZ2)

    ! Total dipole moments of molecules for reaction field
    this%MueX1 => Component1%MueX(:)
    this%MueY1 => Component1%MueY(:)
    this%MueZ1 => Component1%MueZ(:)
    this%MueX2 => Component2%MueX(:)
    this%MueY2 => Component2%MueY(:)
    this%MueZ2 => Component2%MueZ(:)

    ! Torques from reaction field
    if( SimulationType .eq. MolecularDynamics ) then
      this%tRFX1 => Component1%tRFX(:)
      this%tRFY1 => Component1%tRFY(:)
      this%tRFZ1 => Component1%tRFZ(:)
      this%tRFX2 => Component2%tRFX(:)
      this%tRFY2 => Component2%tRFY(:)
      this%tRFZ2 => Component2%tRFZ(:)
    end if

    ! Charge for extended reactionField
    if ((this%N1Charge .gt. 0) .and. (.not. Component1%Molecule%isElongated) )then
       this%lad1 = Component1%Molecule%SiteCharge(1)%e
    end if
    if ((this%N1Charge .gt. 0) .and. (.not. Component2%Molecule%isElongated)) then
       this%lad2 = Component2%Molecule%SiteCharge(1)%e
    end if

    ! Set center of mass positions of test particles
    if( this%NTest1 > 0 ) then
      this%PX1Test => Component1%P0Test(:, 1)
      this%PY1Test => Component1%P0Test(:, 2)
      this%PZ1Test => Component1%P0Test(:, 3)

      ! Total dipole moments of test particles for reaction field
      this%MueX1Test => Component1%MueXTest(:)
      this%MueY1Test => Component1%MueYTest(:)
      this%MueZ1Test => Component1%MueZTest(:)

    else
      nullify( this%PX1Test )
      nullify( this%PY1Test )
      nullify( this%PZ1Test )
      nullify( this%MueX1Test )
      nullify( this%MueY1Test )
      nullify( this%MueZ1Test )
    end if

    ! Set squared cutoff radius
    this%RCutoffSquared = RCutoffLJ126LJ126**2

    ! Create arrays
    call Allocate( this )
    if( (SimulationType .eq. MonteCarlo) .or. (SimulationType .eq. Gibbs) .or. MCOverlapReduction ) then
      this%EPot = 0._RK
      this%EPotNew = 0._RK
      if (this%OptPressure ) then
        this%Virial = 0._RK
        this%VirialNew = 0._RK
      end if
    end if

    ! Nullify pointers
    nullify( this%PotLJ126LJ126 )
    nullify( this%PotChargeCharge )
    nullify( this%PotChargeDipole )
    nullify( this%PotChargeQuadrupole )
    nullify( this%PotDipoleCharge )
    nullify( this%PotDipoleDipole )
    nullify( this%PotDipoleQuadrupole )
    nullify( this%PotQuadrupoleCharge )
    nullify( this%PotQuadrupoleDipole )
    nullify( this%PotQuadrupoleQuadrupole )

    ! Construct Lennard-Jones potentials
    if( this%N1LJ126 > 0 .and. this%N2LJ126 > 0 ) then
      allocate( this%PotLJ126LJ126(this%N1LJ126, this%N2LJ126), STAT = stat )
      call AllocationError( stat, 'sites', this%N1LJ126 + this%N2LJ126 )
      do j1 = 1, this%N1LJ126
        do j2 = 1, this%N2LJ126
          call Construct( this%PotLJ126LJ126(j1, j2), &
&              i1, i2, j1, j2, Component1%Molecule, Component2%Molecule, &
&              RCutoffLJ126LJ126, ScaleSigma, ScaleEpsilon )
        
          this%PotLJ126LJ126(j1, j2)%NInCutoff => this%NInCutoff
          this%PotLJ126LJ126(j1, j2)%CutoffPartner => this%CutoffPartner

          if( RDFUpdateFrequency>0 ) then
            allocate( this%PotLJ126LJ126(j1, j2)%RDFSum(RDFNumberShells+10), STAT = stat )
            call AllocationError( stat, 'RDFSum', RDFNumberShells+10)
          end if
        end do
      end do
    end if

    ! Construct charge-charge potentials
    if( this%N1Charge > 0 .and. this%N2Charge > 0 ) then
      allocate(this%PotChargeCharge(this%N1Charge, this%N2Charge), STAT = stat )
      call AllocationError( stat, 'sites', this%N1Charge + this%N2Charge )
      do j1 = 1, this%N1Charge
        do j2 = 1, this%N2Charge
          call Construct( this%PotChargeCharge(j1, j2), &
&              i1, i2, j1, j2, Component1%Molecule, &
&              Component2%Molecule, RCutoffDipoleDipole, RFEpsilon )

          this%PotChargeCharge(j1, j2)%NInCutoff => this%NInCutoff
          this%PotChargeCharge(j1, j2)%CutoffPartner => this%CutoffPartner
        end do
      end do
    end if

    ! Construct charge-dipole potentials
    if( this%N1Charge > 0 .and. this%N2Dipole > 0 ) then
      allocate(this%PotChargeDipole(this%N1Charge, this%N2Dipole), STAT = stat )
      call AllocationError( stat, 'sites', this%N1Charge + this%N2Dipole )
      do j1 = 1, this%N1Charge
        do j2 = 1, this%N2Dipole
          call Construct( this%PotChargeDipole(j1, j2), &
&              i1, i2, j1, j2, Component1%Molecule, &
&              Component2%Molecule, RCutoffDipoleDipole )

          this%PotChargeDipole(j1, j2)%NInCutoff => this%NInCutoff
          this%PotChargeDipole(j1, j2)%CutoffPartner => this%CutoffPartner
        end do
      end do
    end if

    ! Construct charge-quadrupole potentials
    if( this%N1Charge > 0 .and. this%N2Quadrupole > 0 ) then
      allocate(this%PotChargeQuadrupole(this%N1Charge, this%N2Quadrupole), STAT = stat )
      call AllocationError( stat, 'sites', this%N1Charge + this%N2Quadrupole )
      do j1 = 1, this%N1Charge
        do j2 = 1, this%N2Quadrupole
          call Construct( this%PotChargeQuadrupole(j1, j2), &
&              i1, i2, j1, j2, Component1%Molecule, &
&              Component2%Molecule, RCutoffDipoleDipole )

          this%PotChargeQuadrupole(j1, j2)%NInCutoff => this%NInCutoff
          this%PotChargeQuadrupole(j1, j2)%CutoffPartner => this%CutoffPartner
        end do
      end do
    end if

    ! Construct dipole-charge potentials
    if( this%N1Dipole > 0 .and. this%N2Charge > 0 ) then
      allocate(this%PotDipoleCharge(this%N1Dipole, this%N2Charge), STAT = stat )
      call AllocationError( stat, 'sites', this%N1Dipole + this%N2Charge )
      do j1 = 1, this%N1Dipole
        do j2 = 1, this%N2Charge
          call Construct( this%PotDipoleCharge(j1, j2), &
&              i1, i2, j1, j2, Component1%Molecule, &
&              Component2%Molecule, RCutoffDipoleDipole )

          this%PotDipoleCharge(j1, j2)%NInCutoff => this%NInCutoff
          this%PotDipoleCharge(j1, j2)%CutoffPartner => this%CutoffPartner
        end do
      end do
    end if

    ! Construct dipole-dipole potentials
    if( this%N1Dipole > 0 .and. this%N2Dipole > 0 ) then
      allocate(this%PotDipoleDipole(this%N1Dipole, this%N2Dipole), STAT = stat )
      call AllocationError( stat, 'sites', this%N1Dipole + this%N2Dipole )
      do j1 = 1, this%N1Dipole
        do j2 = 1, this%N2Dipole
          call Construct( this%PotDipoleDipole(j1, j2), &
&              i1, i2, j1, j2, Component1%Molecule, Component2%Molecule, &
&              RCutoffDipoleDipole, RFEpsilon )

          this%PotDipoleDipole(j1, j2)%NInCutoff => this%NInCutoff
          this%PotDipoleDipole(j1, j2)%CutoffPartner => this%CutoffPartner
        end do
      end do
    end if

    ! Construct dipole-quadrupole potentials
    if( this%N1Dipole > 0 .and. this%N2Quadrupole > 0 ) then
      allocate(this%PotDipoleQuadrupole(this%N1Dipole, this%N2Quadrupole), STAT = stat )
      call AllocationError( stat, 'sites', this%N1Dipole + this%N2Quadrupole )
      do j1 = 1, this%N1Dipole
        do j2 = 1, this%N2Quadrupole
          call Construct( this%PotDipoleQuadrupole(j1, j2), &
&              i1, i2, j1, j2, Component1%Molecule, Component2%Molecule, &
&              RCutoffDipoleQuadrupole )

          this%PotDipoleQuadrupole(j1, j2)%NInCutoff => this%NInCutoff
          this%PotDipoleQuadrupole(j1, j2)%CutoffPartner => this%CutoffPartner
        end do
      end do
    end if

    ! Construct quadrupole-charge potentials
    if( this%N1Quadrupole > 0 .and. this%N2Charge > 0 ) then
      allocate(this%PotQuadrupoleCharge(this%N1Quadrupole, this%N2Charge), STAT = stat )
      call AllocationError( stat, 'sites', this%N1Quadrupole + this%N2Charge )
      do j1 = 1, this%N1Quadrupole
        do j2 = 1, this%N2Charge
          call Construct( this%PotQuadrupoleCharge(j1, j2), &
&              i1, i2, j1, j2, Component1%Molecule, &
&              Component2%Molecule, RCutoffDipoleDipole )

          this%PotQuadrupoleCharge(j1, j2)%NInCutoff => this%NInCutoff
          this%PotQuadrupoleCharge(j1, j2)%CutoffPartner => this%CutoffPartner
        end do
      end do
    end if

    ! Construct quadrupole-dipole potentials
    if( this%N1Quadrupole > 0 .and. this%N2Dipole > 0 ) then
      allocate(this%PotQuadrupoleDipole(this%N1Quadrupole, this%N2Dipole), STAT = stat )
      call AllocationError( stat, 'sites', this%N1Quadrupole + this%N2Dipole )
      do j1 = 1, this%N1Quadrupole
        do j2 = 1, this%N2Dipole
          call Construct( this%PotQuadrupoleDipole(j1, j2), &
&              i1, i2, j1, j2, Component1%Molecule, Component2%Molecule, &
&              RCutoffDipoleQuadrupole )

          this%PotQuadrupoleDipole(j1, j2)%NInCutoff => this%NInCutoff
          this%PotQuadrupoleDipole(j1, j2)%CutoffPartner => this%CutoffPartner
        end do
      end do
    end if

    ! Construct quadrupole-quadrupole potentials
    if( this%N1Quadrupole > 0 .and. this%N2Quadrupole > 0 ) then
      allocate(this%PotQuadrupoleQuadrupole(this%N1Quadrupole, this%N2Quadrupole), STAT = stat )
      call AllocationError(stat, 'sites', this%N1Quadrupole + this%N2Quadrupole )
      do j1 = 1, this%N1Quadrupole
        do j2 = 1, this%N2Quadrupole
          call Construct( this%PotQuadrupoleQuadrupole(j1, j2), &
&              i1, i2, j1, j2, Component1%Molecule, Component2%Molecule, &
&              RCutoffQuadrupoleQuadrupole )

          this%PotQuadrupoleQuadrupole(j1, j2)%NInCutoff => this%NInCutoff
          this%PotQuadrupoleQuadrupole(j1, j2)%CutoffPartner => this%CutoffPartner
        end do
      end do
    end if

    ! Set reaction field flag
    if ( (LongRange .eq. RField) .or. (LongRange .eq. ExtRField) ) then
    this%ReactionField = ( CutoffMode .eq. CenterofMass ) .and. &
&                        ((this%N1Charge > 0) .or. (this%N1Dipole > 0)) .and. &
&                        ((this%N2Charge > 0) .or. (this%N2Dipole > 0)) .and. &
&                        .not. ( SimulationType .eq. SecondVirialCoeff )

     else 
        this%ReactionField = .false.
     end if


  end subroutine TInteraction_Construct



!==============================================================!
!  Subroutine TInteraction_Destruct                            !
!==============================================================!

  subroutine TInteraction_Destruct( this )

    implicit none

    ! Declare arguments
    type(TInteraction) :: this

    ! Declare local variables
    integer :: i, j

    ! Destroy Lennard-Jones potentials
    do i = 1, this%N1LJ126
      do j = 1, this%N2LJ126
        call Destruct( this%PotLJ126LJ126(i, j) )
      end do
    end do
    if( associated( this%PotLJ126LJ126 ) ) then
      deallocate( this%PotLJ126LJ126 )
    end if

    ! Destroy charge-charge potentials
    do i = 1, this%N1Charge
      do j = 1, this%N2Charge
        call Destruct( this%PotChargeCharge(i, j) )
      end do
    end do
    if( associated( this%PotChargeCharge ) ) then
      deallocate( this%PotChargeCharge )
    end if

    ! Destroy charge-dipole potentials
    do i = 1, this%N1Charge
      do j = 1, this%N2Dipole
        call Destruct( this%PotChargeDipole(i, j) )
      end do
    end do
    if( associated( this%PotChargeDipole ) ) then
      deallocate( this%PotChargeDipole )
    end if

    ! Destroy charge-quadrupole potentials
    do i = 1, this%N1Charge
      do j = 1, this%N2Quadrupole
        call Destruct( this%PotChargeQuadrupole(i, j) )
      end do
    end do
    if( associated( this%PotChargeQuadrupole ) ) then
      deallocate( this%PotChargeQuadrupole )
    end if

    ! Destroy dipole-charge potentials
    do i = 1, this%N1Dipole
      do j = 1, this%N2Charge
        call Destruct( this%PotDipoleCharge(i, j) )
      end do
    end do
    if( associated( this%PotDipoleCharge ) ) then
      deallocate( this%PotDipoleCharge )
    end if

    ! Destroy dipole-dipole potentials
    do i = 1, this%N1Dipole
      do j = 1, this%N2Dipole
        call Destruct( this%PotDipoleDipole(i, j) )
      end do
    end do
    if( associated( this%PotDipoleDipole ) ) then
      deallocate( this%PotDipoleDipole )
    end if

    ! Destroy dipole-quadrupole potentials
    do i = 1, this%N1Dipole
      do j = 1, this%N2Quadrupole
        call Destruct( this%PotDipoleQuadrupole(i, j) )
      end do
    end do
    if( associated( this%PotDipoleQuadrupole ) ) then
      deallocate( this%PotDipoleQuadrupole )
    end if

    ! Destroy quadrupole-charge potentials
    do i = 1, this%N1Quadrupole
      do j = 1, this%N2Charge
        call Destruct( this%PotQuadrupoleCharge(i, j) )
      end do
    end do
    if( associated( this%PotQuadrupoleCharge ) ) then
      deallocate( this%PotQuadrupoleCharge )
    end if

    ! Destroy quadrupole-dipole potentials
    do i = 1, this%N1Quadrupole
      do j = 1, this%N2Dipole
        call Destruct( this%PotQuadrupoleDipole(i, j) )
      end do
    end do
    if( associated( this%PotQuadrupoleDipole ) ) then
      deallocate( this%PotQuadrupoleDipole )
    end if

    ! Destroy quadrupole-quadrupole potentials
    do i = 1, this%N1Quadrupole
      do j = 1, this%N2Quadrupole
        call Destruct( this%PotQuadrupoleQuadrupole(i, j) )
      end do
    end do
    if( associated( this%PotQuadrupoleQuadrupole ) ) then
      deallocate( this%PotQuadrupoleQuadrupole )
    end if

    ! Destroy arrays
    call Deallocate( this )

  end subroutine TInteraction_Destruct



!==============================================================!
!  Subroutine TInteraction_Allocate                            !
!==============================================================!

  subroutine TInteraction_Allocate( this )

    implicit none

    ! Declare arguments
    type(TInteraction) :: this

    ! Declare local variables
    integer :: N1, N2, stat

    ! Nullify pointers
    nullify( this%EPot )
    nullify( this%EPot1 )
    nullify( this%EPotNew )
    nullify( this%d2EpotdV2 )
    nullify( this%d2EpotdV21 )
    nullify( this%d2EpotdV2New )

    if (this%OptPressure) then
      nullify( this%Virial )
      nullify( this%Virial1 )
      nullify( this%VirialNew )
    end if
    nullify( this%NInCutoff )
    nullify( this%CutoffPartner )
    
    ! allocated only for SimulationType .eq. SecondVirialCoeff
    nullify( this%MayerFFunction )
    nullify( this%MayerFFunction1 )
    nullify( this%MayerFFunction2 )
    nullify( this%IntFFunction )
    nullify( this%IntFFunction1 )
    nullify( this%IntFFunction2 )
    

    ! Calculate dimension of arrays
    if( EnsembleType .eq. EnsembleTypeGE .or. EnsembleType .eq. EnsembleTypeHA .or. SimulationType .eq. Gibbs) then
      N1 = this%NPartMax
      N2 = this%NPartMax

    else
      N1 = max(this%NPart1, 1)
      N2 = max(this%NPart2, 1)
    end if

    ! Allocate arrays
    if( (SimulationType .eq. MonteCarlo) .or. (SimulationType .eq. Gibbs) .or. MCOverlapReduction ) then
      allocate( this%EPot(N1, N2), STAT = stat )
      call AllocationError( stat, 'particles', N1 * N2 )
      allocate( this%EPot1(N2), STAT = stat )
      call AllocationError( stat, 'particles', N2 )
      allocate( this%EPotNew(N1, N2), STAT = stat )
      call AllocationError( stat, 'particles', N1 * N2 )

     allocate( this%d2EpotdV2(N1, N2), STAT = stat )
      call AllocationError( stat, 'particles', N1 * N2 )
      allocate( this%d2EpotdV21(N2), STAT = stat )
      call AllocationError( stat, 'particles', N2 )
      allocate( this%d2EpotdV2New(N1, N2), STAT = stat )
      call AllocationError( stat, 'particles', N1 * N2 )

      if ( this%OptPressure ) then
        allocate( this%Virial(N1, N2), STAT = stat )
        call AllocationError( stat, 'particles', N1 * N2 )
        allocate( this%Virial1(N2), STAT = stat )
        call AllocationError( stat, 'particles', N2 )
        allocate( this%VirialNew(N1, N2), STAT = stat )
        call AllocationError( stat, 'particles', N1 * N2 )
      end if
    end if

    if( SimulationType .eq. SecondVirialCoeff ) then
      allocate( this%EPot1(this%NPartMax), STAT = stat )
      call AllocationError( stat, 'particles', this%NPartMax )
      if ( this%OptPressure ) then
        allocate( this%Virial1(this%NPartMax), STAT = stat )
        call AllocationError( stat, 'particles', this%NPartMax )
      end if

      allocate( this%MayerFFunction(NSteps), STAT = stat )
      call AllocationError( stat, 'Mayer f-function' )
      allocate( this%MayerFFunction1(NSteps), STAT = stat )
      call AllocationError( stat, 'Mayer f-function' )
      allocate( this%MayerFFunction2(NSteps), STAT = stat )
      call AllocationError( stat, 'Mayer f-function' )
      allocate( this%IntFFunction(NSteps), STAT = stat )
      call AllocationError( stat, 'Mayer f-function' )
      allocate( this%IntFFunction1(NSteps), STAT = stat )
      call AllocationError( stat, 'Mayer f-function' )
      allocate( this%IntFFunction2(NSteps), STAT = stat )
      call AllocationError( stat, 'Mayer f-function' )
    end if

    if(( CutoffMode .eq. CenterofMass ) .or. ( CutoffMode .eq. SiteSite ))  then
      N1 = max( N1, this%NTest1 )
      allocate( this%NInCutoff(N1), STAT = stat )
      call AllocationError( stat, 'particles', N1 )
      allocate( this%CutoffPartner(N2, N1), STAT = stat )
      call AllocationError( stat, 'particles', N1 * N2 )
    end if

  end subroutine TInteraction_Allocate



!==============================================================!
!  Subroutine TInteraction_DeallocateEPot                      !
!==============================================================!

  subroutine TInteraction_DeallocateEPot( this )

    implicit none

    ! Declare arguments
    type(TInteraction) :: this

    ! Deallocate arrays
    if( associated( this%EPot ) ) then
      deallocate( this%EPot )
    end if
    if( associated( this%EPot1 ) ) then
      deallocate( this%EPot1 )
    end if
    if( associated( this%EPotNew ) ) then
      deallocate( this%EPotNew )
    end if

    if( associated( this%d2EpotdV2 ) ) then
      deallocate( this%d2EpotdV2 )
    end if
    if( associated( this%d2EpotdV21 ) ) then
      deallocate( this%d2EpotdV21 )
    end if
    if( associated( this%d2EpotdV2New ) ) then
      deallocate( this%d2EpotdV2New )	  
    end if

    if ( this%OptPressure ) then
      if( associated( this%Virial ) ) then
        deallocate( this%Virial )
      end if
      if( associated( this%Virial1 ) ) then
        deallocate( this%Virial1 )
      end if
      if( associated( this%VirialNew ) ) then
        deallocate( this%VirialNew )
      end if
    end if

  end subroutine TInteraction_DeallocateEPot



!==============================================================!
!  Subroutine TInteraction_Deallocate                          !
!==============================================================!

  subroutine TInteraction_Deallocate( this )

    implicit none

    ! Declare arguments
    type(TInteraction) :: this

    ! Deallocate arrays
    call DeallocateEPot( this )

    ! allocated only for SimulationType .eq. SecondVirialCoeff
    if( associated( this%MayerFFunction ) ) then
      deallocate( this%MayerFFunction )
    end if
    if( associated( this%MayerFFunction1 ) ) then
      deallocate( this%MayerFFunction1 )
    end if
    if( associated( this%MayerFFunction2 ) ) then
      deallocate( this%MayerFFunction2 )
    end if
    if( associated( this%IntFFunction ) ) then
      deallocate( this%IntFFunction )
    end if
    if( associated( this%IntFFunction1 ) ) then
      deallocate( this%IntFFunction1 )
    end if
    if( associated( this%IntFFunction2 ) ) then
      deallocate( this%IntFFunction2 )
    end if
    if( associated( this%NInCutoff ) ) then
      deallocate( this%NInCutoff )
    end if
    if( associated( this%CutoffPartner ) ) then
      deallocate( this%CutoffPartner )
    end if

  end subroutine TInteraction_Deallocate

!==============================================================!
!  Subroutine TInteraction_RDF                               !
!==============================================================!

  subroutine TInteraction_RDF( this, BoxLength,RDFdr )

    implicit none

    ! Declare arguments
    type(TInteraction)       :: this
    real(RK), intent(in)     :: RDFdr
    real(RK), intent(in)     :: BoxLength

    ! Declare local variables

    integer           :: i, j


    ! Calculate interactions partners within cutoff sphere
      call CalcCutoffPartnersRDF( this )

    ! Calculate Lennard-Jones forces
    do i = 1, this%N1LJ126
      do j = 1, this%N2LJ126
        call GET_RDF( this%PotLJ126LJ126( i, j ), BoxLength,RDFdr )
      end do
    end do
    
 end subroutine TInteraction_RDF


!==============================================================!
!  Subroutine TInteraction_Force                               !
!==============================================================!
#ifndef ABL
  subroutine TInteraction_Force( this, EPot, Virial, d2EpotdV2, BoxLength )
#else
  subroutine TInteraction_Force( this, EPot, Virial, d2EpotdV2, BoxLength,C1,C2)
#endif

    implicit none

    ! Declare arguments
    type(TInteraction)       :: this
    real(RK), intent(in out) :: EPot
    real(RK), intent(in out) :: Virial
    real(RK), intent(in out) :: d2EpotdV2
    real(RK), intent(in)     :: BoxLength
#ifdef ABL
    integer, intent(in)      :: C1,C2
#endif

    ! Declare local variables
    real(RK), pointer :: MueX1(:), MueY1(:), MueZ1(:)
    real(RK), pointer :: MueX2(:), MueY2(:), MueZ2(:)
    real(RK), pointer :: TX1(:), TY1(:), TZ1(:), TX2(:), TY2(:), TZ2(:)
    real(RK)          :: mueXi, mueYi, mueZi, mueXj, mueYj, mueZj
    real(RK)          :: RFTX, RFTY, RFTZ
    real(RK)          :: EPotLocal, TXi, TYi, TZi
    integer           :: i, j, k, i1
#if MPI_VER > 0
    integer           :: i0
#endif
#ifdef ABL
    real(RK)          :: AblSig, AblEps
    real(RK)          :: eps1,eps2,fac
#endif

    ! Calculate interactions partners within cutoff sphere
    if( CutoffMode .eq. CenterofMass ) then
      call CalcCutoffPartners( this )
    end if

    ! Calculate Lennard-Jones forces
    do i = 1, this%N1LJ126
      do j = 1, this%N2LJ126
#ifndef ABL
#if  TRANS == 1
       if(.not. Equilibration .and. (mod((Step+NStepCorr-1),NStepCorr) .eq. 0)) then
         call Force_Trans( this%PotLJ126LJ126( i, j ), EPot, Virial, d2EpotdV2, BoxLength )

       else
         call Force( this%PotLJ126LJ126( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       end if

#else
       call Force( this%PotLJ126LJ126( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
#endif

#else
#if  TRANS == 1
       if(.not. Equilibration .and. (mod((Step+NStepCorr-1),NStepCorr) .eq. 0)) then
         call Force_Trans( this%PotLJ126LJ126( i, j ), &
&             EPot, Virial, d2EpotdV2, BoxLength,  AblSig, AblEps,eps1,eps2)

       else
         call Force( this%PotLJ126LJ126( i, j ), &
&             EPot, Virial, d2EpotdV2, BoxLength,  AblSig, AblEps,eps1,eps2)
       end if
#else
       call Force( this%PotLJ126LJ126( i, j ), &
&           EPot, Virial, d2EpotdV2, BoxLength,  AblSig, AblEps,eps1,eps2)
#endif

        this%AblPS(C1,i)    = this%AblPS(C1,i) + AblSig
        this%AblPS(C2,j)    = this%AblPS(C2,j) + AblSig

        if ( (C1.eq. C2) .AND. (i.eq.j) ) then 
          this%AblPE(C1,i)    = this%AblPE(C1,i) + AblEps

        else
          fac = 2._RK*sqrt(eps1*eps2)
          this%AblPE(C1,i)    = this%AblPE(C1,i) + AblEps*eps2 / fac
          this%AblPE(C2,j)    = this%AblPE(C2,j) + AblEps*eps1 / fac
        end if
#endif
      end do
    end do

    ! Calculate point charge forces
    do i = 1, this%N1Charge
      if ( .not. this%ReactionField ) then
        do j = 1, this%N2Charge

#if  TRANS == 1
          if(.not. Equilibration .and. (mod((Step+NStepCorr-1),NStepCorr) .eq. 0)) then
               call Force_Trans( this%PotChargeCharge( i, j ), EPot, Virial, d2EpotdV2, BoxLength, this%Kappa )

          else
               call Force( this%PotChargeCharge( i, j ), EPot, Virial, d2EpotdV2, BoxLength, this%Kappa )
          end if

#else
          call Force( this%PotChargeCharge( i, j ), EPot, Virial, d2EpotdV2, BoxLength, this%Kappa )
#endif
        end do

      else
        do j = 1, this%N2Charge

#if  TRANS == 1
         if(.not. Equilibration .and. (mod((Step+NStepCorr-1),NStepCorr) .eq. 0)) then
            call Force_Trans( this%PotChargeCharge( i, j ), EPot, Virial, d2EpotdV2, BoxLength )

         else
            call Force( this%PotChargeCharge( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
         end if

#else
            call Force( this%PotChargeCharge( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
#endif
        end do
      end if

      do j = 1, this%N2Dipole
#if  TRANS == 1
       if(.not. Equilibration .and. (mod((Step+NStepCorr-1),NStepCorr) .eq. 0)) then
          call Force_Trans( this%PotChargeDipole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )

       else
          call Force( this%PotChargeDipole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       end if

#else
          call Force( this%PotChargeDipole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
#endif
      end do

      do j = 1, this%N2Quadrupole
#if  TRANS == 1
       if(.not. Equilibration .and. (mod((Step+NStepCorr-1),NStepCorr) .eq. 0)) then
          call Force_Trans( this%PotChargeQuadrupole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       else
          call Force( this%PotChargeQuadrupole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       end if

#else
          call Force( this%PotChargeQuadrupole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
#endif
      end do
    end do


    ! Calculate dipolar forces
    do i = 1, this%N1Dipole
      do j = 1, this%N2Charge
#if  TRANS == 1
       if(.not. Equilibration .and. (mod((Step+NStepCorr-1),NStepCorr) .eq. 0)) then
          call Force_Trans( this%PotDipoleCharge( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       else
          call Force( this%PotDipoleCharge( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       end if
#else
          call Force( this%PotDipoleCharge( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
#endif
      end do

      do j = 1, this%N2Dipole
#if  TRANS == 1
       if(.not. Equilibration .and. (mod((Step+NStepCorr-1),NStepCorr) .eq. 0)) then
          call Force_Trans( this%PotDipoleDipole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       else
          call Force( this%PotDipoleDipole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       end if
#else
          call Force( this%PotDipoleDipole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
#endif
      end do

      do j = 1, this%N2Quadrupole
#if  TRANS == 1
       if(.not. Equilibration .and. (mod((Step+NStepCorr-1),NStepCorr) .eq. 0)) then
          call Force_Trans( this%PotDipoleQuadrupole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       else
          call Force( this%PotDipoleQuadrupole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       end if
#else
          call Force( this%PotDipoleQuadrupole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
#endif
      end do
    end do


    ! Calculate quadrupolar forces
    do i = 1, this%N1Quadrupole
      do j = 1, this%N2Charge
#if  TRANS == 1
       if(.not. Equilibration .and. (mod((Step+NStepCorr-1),NStepCorr) .eq. 0)) then
          call Force_Trans( this%PotQuadrupoleCharge( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       else
          call Force( this%PotQuadrupoleCharge( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       end if
#else
          call Force( this%PotQuadrupoleCharge( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
#endif
      end do

      do j = 1, this%N2Dipole
#if  TRANS == 1
       if(.not. Equilibration .and. (mod((Step+NStepCorr-1),NStepCorr) .eq. 0)) then
          call Force_Trans( this%PotQuadrupoleDipole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       else
          call Force( this%PotQuadrupoleDipole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       end if
#else
          call Force( this%PotQuadrupoleDipole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
#endif
      end do

      do j = 1, this%N2Quadrupole
#if  TRANS == 1
       if(.not. Equilibration .and. (mod((Step+NStepCorr-1),NStepCorr) .eq. 0)) then
          call Force_Trans( this%PotQuadrupoleQuadrupole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       else
          call Force( this%PotQuadrupoleQuadrupole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
       end if
#else
          call Force( this%PotQuadrupoleQuadrupole( i, j ), EPot, Virial, d2EpotdV2, BoxLength )
#endif
      end do
    end do

    ! Explicit reaction field contribution
    if ( this%ReactionField ) then
      MueX1 => this%MueX1
      MueY1 => this%MueY1
      MueZ1 => this%MueZ1
      MueX2 => this%MueX2
      MueY2 => this%MueY2
      MueZ2 => this%MueZ2
      TX1 => this%tRFX1
      TY1 => this%tRFY1
      TZ1 => this%tRFZ1
      TX2 => this%tRFX2
      TY2 => this%tRFY2
      TZ2 => this%tRFZ2
      EPotLocal = 0._RK

#if MPI_VER > 0
      i0 = this%NPart10
      i1 = this%NPart12
      do i = i0, i1
#else
      i1 = this%NPart1
      do i = 1, i1
#endif
        TXi = 0._RK
        TYi = 0._RK
        TZi = 0._RK
        mueXi = MueX1(i)
        mueYi = MueY1(i)
        mueZi = MueZ1(i)
        do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          mueXj = MueX2(j)
          mueYj = MueY2(j)
          mueZj = MueZ2(j)
          RFTX = this%RFConst2 * (mueZi * mueYj - mueYi * mueZj)
          RFTY = this%RFConst2 * (mueXi * mueZj - mueZi * mueXj)
          RFTZ = this%RFConst2 * (mueYi * mueXj - mueXi * mueYj)
          TXi = TXi + RFTX
          TYi = TYi + RFTY
          TZi = TZi + RFTZ
          TX2(j) = TX2(j) - RFTX
          TY2(j) = TY2(j) - RFTY
          TZ2(j) = TZ2(j) - RFTZ

          EPotLocal = EPotLocal + (mueXi * mueXj + mueYi * mueYj + mueZi * mueZj)
        end do
        TX1(i) = TX1(i) + TXi
        TY1(i) = TY1(i) + TYi
        TZ1(i) = TZ1(i) + TZi
      end do

      EPot = EPot + this%RFConst2 * EPotlocal
    end if

  end subroutine TInteraction_Force



!==============================================================!
!  Subroutine TInteraction_ChemicalPotential                   !
!==============================================================!

  subroutine TInteraction_ChemicalPotential( this, EPotTest, BoxLength )

    implicit none

    ! Declare arguments
    type(TInteraction)   :: this
    real(RK), pointer    :: EPotTest(:)
    real(RK), intent(in) :: BoxLength

    ! Declare local variables
    real(RK), pointer :: MueX1(:), MueY1(:), MueZ1(:)
    real(RK), pointer :: MueX2(:), MueY2(:), MueZ2(:)
    real(RK)          :: mueXi, mueYi, mueZi
    real(RK)          :: EPotLocal
    integer           :: i, j, k

    ! Calculate interactions partners within cutoff sphere
    if( CutoffMode .eq. CenterofMass ) then
      call CalcCutoffPartnersTest( this )
    end if

    ! Calculate Lennard-Jones chemical potential
    do i = 1, this%N1LJ126
      do j = 1, this%N2LJ126
        call ChemicalPotential( this%PotLJ126LJ126( i, j ), EPotTest, BoxLength )
      end do
    end do

    ! Calculate point charge chemical potential
    do i = 1, this%N1Charge
        do j = 1, this%N2Charge
          call ChemicalPotential( this%PotChargeCharge( i, j ), EPotTest, BoxLength )
        end do
      do j = 1, this%N2Dipole
        call ChemicalPotential( this%PotChargeDipole( i, j ), EPotTest, BoxLength )
      end do
      do j = 1, this%N2Quadrupole
        call ChemicalPotential( this%PotChargeQuadrupole( i, j ), EPotTest, BoxLength )
      end do
    end do

    ! Calculate dipolar chemical potential
    do i = 1, this%N1Dipole
      do j = 1, this%N2Charge
        call ChemicalPotential( this%PotDipoleCharge( i, j ), EPotTest, BoxLength )
      end do
      do j = 1, this%N2Dipole
        call ChemicalPotential( this%PotDipoleDipole( i, j ), EPotTest, BoxLength )
      end do
      do j = 1, this%N2Quadrupole
        call ChemicalPotential( this%PotDipoleQuadrupole( i, j ), EPotTest, BoxLength )
      end do
    end do

    ! Calculate quadrupolar chemical potential
    do i = 1, this%N1Quadrupole
      do j = 1, this%N2Charge
        call ChemicalPotential( this%PotQuadrupoleCharge( i, j ), EPotTest, BoxLength )
      end do
      do j = 1, this%N2Dipole
        call ChemicalPotential( this%PotQuadrupoleDipole( i, j ), EPotTest, BoxLength )
      end do
      do j = 1, this%N2Quadrupole
        call ChemicalPotential( this%PotQuadrupoleQuadrupole( i, j ), EPotTest, BoxLength )
      end do
    end do

    ! Explicit reaction field contribution
    if ( this%ReactionField ) then
      MueX1 => this%MueX1Test
      MueY1 => this%MueY1Test
      MueZ1 => this%MueZ1Test
      MueX2 => this%MueX2
      MueY2 => this%MueY2
      MueZ2 => this%MueZ2

      do i = 1, this%NTest1
        mueXi = MueX1(i)
        mueYi = MueY1(i)
        mueZi = MueZ1(i)
        EPotLocal = 0._RK
        do k = 1, this%NInCutoff(i)
          j = this%CutoffPartner(k, i)
          EPotLocal = EPotLocal + ( mueXi * MueX2(j) + mueYi * MueY2(j) + mueZi * MueZ2(j) )
        end do
        EPotTest(i) = EPotTest(i) + this%RFConst2 * EPotLocal
      end do
    end if

  end subroutine TInteraction_ChemicalPotential



!==============================================================!
!  Subroutine TInteraction_Energy                              !
!==============================================================!

  subroutine TInteraction_Energy( this, np, BoxLength )

    implicit none

    ! Declare arguments
    type(TInteraction)   :: this
    integer, intent(in)  :: np
    real(RK), intent(in) :: BoxLength

    ! Declare local variables
    type(TPotLJ126LJ126), pointer           :: plj
    type(TPotChargeCharge), pointer         :: pcc
    type(TPotChargeDipole), pointer         :: pcd
    type(TPotChargeQuadrupole), pointer     :: pcq
    type(TPotDipoleCharge), pointer         :: pdc
    type(TPotDipoleDipole), pointer         :: pdd
    type(TPotDipoleQuadrupole), pointer     :: pdq
    type(TPotQuadrupoleCharge), pointer     :: pqc
    type(TPotQuadrupoleDipole), pointer     :: pqd
    type(TPotQuadrupoleQuadrupole), pointer :: pqq
    real(RK) :: EPot(this%NPart2)
    real(RK), pointer :: Virial(:)
    real(RK), pointer :: d2EpotdV2(:)
    real(RK)          :: EPotLocal
    real(RK)          :: VirialLocal

    real(RK)          :: SigmaSquared
    real(RK)          :: Epsilon, Epsilon2, Epsilon4, Epsilon48
    real(RK)          :: RCutoffSquared, RCutoffSquaredScaled, RShieldSquared
    real(RK)          :: BoxLengthThird
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK), pointer :: OX1(:), OY1(:), OZ1(:), OX2(:), OY2(:), OZ2(:)
    real(RK)          :: RXi, RYi, RZi
    real(RK)          :: PXi, PYi, PZi
    real(RK)          :: OXi, OYi, OZi
    real(RK)          :: RXij, RYij, RZij
    real(RK)          :: FXij, FYij, FZij, Fij
    real(RK)          :: PXij, PYij, PZij
    real(RK)          :: OXj, OYj, OZj
    real(RK)          :: eX, eY, eZ
    real(RK)          :: RijSquared, RijInv, RijSquaredInv, Rij3Inv
    real(RK)          :: RijInv4, RijInv2
    real(RK)          :: Rij4Inv, Rij4Inv3, Rij5Inv, Rij6Inv
    real(RK)          :: CosThetai, CosThetaj
    real(RK)          :: CosThetaiSquared, CosThetajSquared
    real(RK)          :: CosAux, CosGammaij
    real(RK)          :: dCosThetai, dCosThetaj, dCosGammaij
    real(RK)          :: Tmp, RFConst2
    real(RK), pointer :: MueX2(:), MueY2(:), MueZ2(:)
    real(RK)          :: mueXi, mueYi, mueZi
    real(RK)          :: sitecorr, Plen2
    real(RK)          :: KappaRij, Rij, approx, Faktor, q
    integer           :: N
    integer           :: s1, s2, j, k
    logical           :: OptPressure

    this%EPot1(:)=0._RK

    ! Zero energy
    EPot(:)= 0._RK

    ! Calculate interactions partners within cutoff sphere
    if( CutoffMode .eq. CenterofMass ) then
      call CalcCutoffPartners( this, np )
    end if
      
    d2EpotdV2 => this%d2EpotdV21


    ! Assign local variables
    OptPressure = this%OptPressure
    if ( OptPressure ) then
      Virial => this%Virial1
      VirialLocal = 1E33_RK
    end if

    N = this%NPart2
    RCutoffSquared = this%RCutoffSquared
    RCutoffSquaredScaled = this%RCutoffSquaredScaled
    BoxLengthThird = Third * BoxLength
    PXi = this%PX1(np)
    PYi = this%PY1(np)
    PZi = this%PZ1(np)

    ! Assign pointers to COM positions
    PX2 => this%PX2
    PY2 => this%PY2
    PZ2 => this%PZ2

    ! d2Epot/dV2
    d2EpotdV2(:) = 0._RK

    if ( OptPressure ) then
      ! Zero virial
      Virial(:) = 0._RK
    end if

    ! Initialization Ewald Summation
    if ( .not. this%ReactionField ) then
       Faktor = 2._RK/sqrt(Pi) * this%Kappa
    end if 

    if( CutoffMode .eq. CenterofMass ) then


      ! Calculate Lennard-Jones energy
      do s1 = 1, this%N1LJ126
        do s2 = 1, this%N2LJ126

          ! Set site specific variables
          plj => this%PotLJ126LJ126(s1, s2)
          SigmaSquared = plj%SigmaSquared
          Epsilon4 = plj%Epsilon4

          if ( OptPressure ) then
            Epsilon48 = plj%Epsilon48
          end if

          ! Assign pointers to site positions
          RX1 => plj%Site1%RX
          RY1 => plj%Site1%RY
          RZ1 => plj%Site1%RZ
          RX2 => plj%Site2%RX
          RY2 => plj%Site2%RY
          RZ2 => plj%Site2%RZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)

          ! Loop over molecules
!CDIR NODEP
          do k = 1, this%NInCutoff(np)
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
            if ( OptPressure ) then
              Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
              FXij = Fij * RXij
              FYij = Fij * RYij
              FZij = Fij * RZij
              Virial(j) = Virial(j) + BoxLengthThird * (PXij * FXij + PYij * FYij + PZij * FZij)
            end if
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)/RijSquared
            d2EpotdV2(j) = d2EpotdV2(j) + Epsilon4 * Rij6Inv *(12._RK *Rij6Inv -  6._RK) * (sitecorr * sitecorr - Plen2/RijSquared)*Third*Third !xxxx2 LJ
            d2EpotdV2(j) = d2EpotdV2(j) + Epsilon4 * Rij6Inv *(156._RK*Rij6Inv - 42._RK) *  sitecorr * sitecorr *Third*Third
          end do
        end do
      end do

      ! Calculate point charge energy
      do s1 = 1, this%N1Charge
! Ewald-Summation
       if ( .not. this%ReactionField ) then
        do s2 = 1, this%N2Charge
          pcc => this%PotChargeCharge(s1, s2)
          Epsilon = pcc%Epsilon
          RShieldSquared = pcc%RShieldSquared

          ! Assign pointers to site positions
          RX1 => pcc%Site1%RX
          RY1 => pcc%Site1%RY
          RZ1 => pcc%Site1%RZ
          RX2 => pcc%Site2%RX
          RY2 => pcc%Site2%RY
          RZ2 => pcc%Site2%RZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)

          ! Loop over molecules
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
                Rij =  sqrt(RijSquared)
                RijInv = 1._RK /  Rij 
                KappaRij = this%Kappa*Rij
                call ErrorApprox(this%PotChargeCharge(s1,s2), KappaRij,approx)
                EPotLocal = Epsilon * RijInv * approx
                if ( OptPressure ) then
                  eX = RXij * RijInv
                  eY = RYij * RijInv
                  eZ = RZij * RijInv
                  VirialLocal = (EPotLocal + Faktor*exp(-KappaRij**2) * Epsilon) &
&                               * RijInv * (eX * PXij + eY * PYij + eZ * PZij)
                end if
            end if
            EPot(j) = EPot(j) + EPotLocal
            if ( OptPressure ) then
              Virial(j) = Virial(j) + Third * VirialLocal
            end if
          end do
        end do
! Reaction Field
       else
        do s2 = 1, this%N2Charge
          pcc => this%PotChargeCharge(s1, s2)
          Epsilon = pcc%Epsilon
          RShieldSquared = pcc%RShieldSquared

          ! Assign pointers to site positions
          RX1 => pcc%Site1%RX
          RY1 => pcc%Site1%RY
          RZ1 => pcc%Site1%RZ
          RX2 => pcc%Site2%RX
          RY2 => pcc%Site2%RY
          RZ2 => pcc%Site2%RZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)

          ! Loop over molecules
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
              EPotLocal = Epsilon * RijInv
              if ( OptPressure ) then
                eX = RXij * RijInv
                eY = RYij * RijInv
                eZ = RZij * RijInv
                VirialLocal = EPotLocal * RijInv * (eX * PXij + eY * PYij + eZ * PZij)
              end if
            end if
            EPot(j) = EPot(j) + EPotLocal
            if ( OptPressure ) then
              Virial(j) = Virial(j) + Third * VirialLocal
            end if
            RijInv2  =  RijInv*RijInv
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (RXij*PXij+RYij*PYij+RZij*PZij)*RijInv2
            d2EpotdV2(j) = d2EpotdV2(j) + EPotLocal * (3._RK * sitecorr*sitecorr - Plen2*RijInv2)*Third*Third !xxxx2 CC
          end do
        end do
       end if ! ReactionField - Ewald-Summation
!
        do s2 = 1, this%N2Dipole
          pcd => this%PotChargeDipole(s1, s2)
          Epsilon = pcd%Epsilon
          RShieldSquared = pcd%RShieldSquared

          ! Assign pointers to site positions
          RX1 => pcd%Site1%RX
          RY1 => pcd%Site1%RY
          RZ1 => pcd%Site1%RZ
          RX2 => pcd%Site2%RX
          RY2 => pcd%Site2%RY
          RZ2 => pcd%Site2%RZ
          OX2 => pcd%Site2%OX
          OY2 => pcd%Site2%OY
          OZ2 => pcd%Site2%OZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)

          ! Loop over molecules
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
              CosThetaj = OXj * ex + OYj * eY + OZj * eZ
              EPotLocal = Epsilon * RijSquaredInv * CosThetaj
              if ( OptPressure ) then
                Tmp = 3._RK * CosThetaj
                VirialLocal = Epsilon * RijSquaredInv * RijInv &
&                              * ( ( Tmp * eX - OXj ) * PXij &
&                                + ( Tmp * eY - OYj ) * PYij &
&                                + ( Tmp * eZ - OZj ) * PZij )
              end if
            end if
            EPot(j) = EPot(j) + EPotLocal
            if ( OptPressure ) then
              Virial(j) = Virial(j) + Third * VirialLocal
            end if
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
            d2EpotdV2(j) = d2EpotdV2(j) + EPotLocal*(8._RK*sitecorr*sitecorr-2._RK*Plen2*RijSquaredInv)*Third*Third !xxxx2 CD

          end do
        end do
        do s2 = 1, this%N2Quadrupole
          pcq => this%PotChargeQuadrupole(s1, s2)
          Epsilon = pcq%Epsilon
          RShieldSquared = pcq%RShieldSquared

          ! Assign pointers to site positions
          RX1 => pcq%Site1%RX
          RY1 => pcq%Site1%RY
          RZ1 => pcq%Site1%RZ
          RX2 => pcq%Site2%RX
          RY2 => pcq%Site2%RY
          RZ2 => pcq%Site2%RZ
          OX2 => pcq%Site2%OX
          OY2 => pcq%Site2%OY
          OZ2 => pcq%Site2%OZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)

          ! Loop over molecules
!CDIR NODEP
          do k = 1, this%NInCutoff(np)
            j = this%CutoffPartner(k, np)
            RXij = RXi - RX2(j)
            RYij = RYi - RY2(j)
            RZij = RZi - RZ2(j)
            PXij = PXi - PX2(j)
            PYij = PYi - PY2(j)
            PZij = PZi - PZ2(j)
            RXij = (RXij - anint( PXij )) * BoxLength                           ! Abstandsvektor von Q nach C wie bei Price
            RYij = (RYij - anint( PYij )) * BoxLength
            RZij = (RZij - anint( PZij )) * BoxLength
            PXij = (PXij - anint( PXij )) * BoxLength
            PYij = (PYij - anint( PYij )) * BoxLength
            PZij = (PZij - anint( PZij )) * BoxLength                           ! Orientierungsvektor Quadrupol
            OXj = OX2(j)
            OYj = OY2(j)
            OZj = OZ2(j)
            RijSquared = RXij**2 + RYij**2 + RZij**2
            if( RijSquared <= RShieldSquared ) then
              EPotLocal = 1E33_RK
            else
              RijSquaredInv = 1._RK / RijSquared
              RijInv = sqrt( RijSquaredInv )
              eX = RXij * RijInv                                                ! Normierter Abstandsvektor
              eY = RYij * RijInv
              eZ = RZij * RijInv
              CosThetaj = OXj * ex + OYj * eY + OZj * eZ
              EPotLocal = Epsilon * RijSquaredInv * RijInv * ( CosThetaj * CosThetaj - Third )

              if ( OptPressure ) then
                Tmp = 2._RK * CosThetaj
                CosAux = 5._RK * CosThetaj * CosThetaj - 1._RK
                VirialLocal =  Epsilon * RijSquaredInv * RijSquaredInv &
&                              * ( ( CosAux * eX - Tmp * OXj ) * PXij &
&                              + ( CosAux * eY - Tmp * OYj ) * PYij &
&                              + ( CosAux * eZ - Tmp * OZj ) * PZij )
              end if
            end if
            EPot(j) = EPot(j) + EPotLocal
            if ( OptPressure ) then
              Virial(j) = Virial(j) + Third * VirialLocal
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
            d2EpotdV2(j) = d2EpotdV2(j) + EPotLocal*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijSquaredInv)*Third*Third !xxxx3 CQ
            end if
          end do
        end do
      end do

      ! Calculate dipolar energy
      do s1 = 1, this%N1Dipole
        do s2 = 1, this%N2Charge
          pdc => this%PotDipoleCharge(s1, s2)
          Epsilon = pdc%Epsilon
          RShieldSquared = pdc%RShieldSquared

          ! Assign pointers to site positions
          RX1 => pdc%Site1%RX
          RY1 => pdc%Site1%RY
          RZ1 => pdc%Site1%RZ
          OX1 => pdc%Site1%OX
          OY1 => pdc%Site1%OY
          OZ1 => pdc%Site1%OZ
          RX2 => pdc%Site2%RX
          RY2 => pdc%Site2%RY
          RZ2 => pdc%Site2%RZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)
          OXi = OX1(np)
          OYi = OY1(np)
          OZi = OZ1(np)

          ! Loop over molecules
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
              CosThetai = OXi * ex + OYi * eY + OZi * eZ
              Tmp = 3._RK * CosThetai
              EPotLocal = - Epsilon * RijSquaredInv * CosThetai
              if ( OptPressure ) then
                VirialLocal = Epsilon * RijSquaredInv * RijInv &
&                              * ( ( OXi - Tmp * eX ) * PXij &
&                                + ( OYi - Tmp * eY ) * PYij &
&                                + ( OZi - Tmp * eZ ) * PZij )
              end if
            end if
            EPot(j) = EPot(j) + EPotLocal
            if ( OptPressure ) then
              Virial(j) = Virial(j) + Third * VirialLocal
            end if
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
            d2EpotdV2(j) = d2EpotdV2(j) + EPotLocal*(8._RK*sitecorr*sitecorr-2._RK*Plen2*RijSquaredInv)*Third*Third !xxxx4 DC
          end do
        end do
        do s2 = 1, this%N2Dipole
          pdd => this%PotDipoleDipole(s1, s2)
          Epsilon = pdd%Epsilon
          RShieldSquared = pdd%RShieldSquared

          ! Assign pointers to site positions
          RX1 => pdd%Site1%RX
          RY1 => pdd%Site1%RY
          RZ1 => pdd%Site1%RZ
          OX1 => pdd%Site1%OX
          OY1 => pdd%Site1%OY
          OZ1 => pdd%Site1%OZ
          RX2 => pdd%Site2%RX
          RY2 => pdd%Site2%RY
          RZ2 => pdd%Site2%RZ
          OX2 => pdd%Site2%OX
          OY2 => pdd%Site2%OY
          OZ2 => pdd%Site2%OZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)
          OXi = OX1(np)
          OYi = OY1(np)
          OZi = OZ1(np)

          ! Loop over molecules
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
              EPotLocal = Rij3Inv * Tmp
              if ( OptPressure ) then
                Rij4Inv3 = 3._RK * Rij3Inv * RijInv
                FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                           - (eX * CosThetaj - OXj) * CosThetai)
                FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                           - (eY * CosThetaj - OYj) * CosThetai)
                FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                           - (eZ * CosThetaj - OZj) * CosThetai)
              end if
            end if
            EPot(j) = EPot(j) + EPotLocal
            if ( OptPressure ) then
              Virial(j) = Virial(j) + Third * ( FXij * PXij + FYij * PYij + FZij * PZij )
            end if
            RijInv2  =  RijInv*RijInv
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv2
            d2EpotdV2(j) = d2EpotdV2(j) + EPotLocal*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijInv2)*Third*Third !xxxx5 DD
          end do
        end do
        do s2 = 1, this%N2Quadrupole
          pdq => this%PotDipoleQuadrupole(s1, s2)
          Epsilon = pdq%Epsilon
          RShieldSquared = pdq%RShieldSquared

          ! Assign pointers to site positions
          RX1 => pdq%Site1%RX
          RY1 => pdq%Site1%RY
          RZ1 => pdq%Site1%RZ
          OX1 => pdq%Site1%OX
          OY1 => pdq%Site1%OY
          OZ1 => pdq%Site1%OZ
          RX2 => pdq%Site2%RX
          RY2 => pdq%Site2%RY
          RZ2 => pdq%Site2%RZ
          OX2 => pdq%Site2%OX
          OY2 => pdq%Site2%OY
          OZ2 => pdq%Site2%OZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)
          OXi = OX1(np)
          OYi = OY1(np)
          OZi = OZ1(np)

          ! Loop over molecules
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
              CosAux = 1._RK - 5._RK * CosThetaj**2
              CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
              Rij4Inv = Epsilon / RijSquared**2
              EPotLocal = Rij4Inv * ( CosGammaij * CosThetaj &
&                                   + CosThetai * CosAux )

              if ( OptPressure ) then
                dCosThetai = Rij4Inv * CosAux
                dCosThetaj = Rij4Inv * (CosGammaij - 10._RK * CosThetai * CosThetaj)
                dCosGammaij = 2._RK * Rij4Inv * CosThetaj
                Tmp = -4._RK * RijInv * EPotLocal
                FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                          + (eX * CosThetaj - OXj) * dCosThetaj)
                FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                          + (eY * CosThetaj - OYj) * dCosThetaj)
                FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                          + (eZ * CosThetaj - OZj) * dCosThetaj)
              end if
            end if
            EPot(j) = EPot(j) + EPotLocal
            if ( OptPressure ) then
              Virial(j) = Virial(j) + Third * ( FXij * PXij + FYij * PYij + FZij * PZij )
            end if
            RijInv2  =  RijInv*RijInv
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv2
            d2EpotdV2(j) = d2EpotdV2(j) + EPotLocal*(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv2)*Third*Third !xxxx6 DQ
          end do
        end do
      end do

      ! Calculate quadrupolar energy
      do s1 = 1, this%N1Quadrupole
        do s2 = 1, this%N2Charge
          pqc => this%PotQuadrupoleCharge(s1, s2)
          Epsilon = pqc%Epsilon
          RShieldSquared = pqc%RShieldSquared

          ! Assign pointers to site positions
          RX1 => pqc%Site1%RX
          RY1 => pqc%Site1%RY
          RZ1 => pqc%Site1%RZ
          OX1 => pqc%Site1%OX
          OY1 => pqc%Site1%OY
          OZ1 => pqc%Site1%OZ
          RX2 => pqc%Site2%RX
          RY2 => pqc%Site2%RY
          RZ2 => pqc%Site2%RZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)
          OXi = OX1(np)
          OYi = OY1(np)
          OZi = OZ1(np)

          ! Loop over molecules
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
              eX = - RXij * RijInv                                              ! Normierter Abstandsvektor nach Price
              eY = - RYij * RijInv
              eZ = - RZij * RijInv
              CosThetai = OXi * ex + OYi * eY + OZi * eZ        ! Scalarprodukt normierter Abstandsvektor mit 
!                                                              Orientierungsvektor Quadrupol
              EPotLocal = Epsilon * RijSquaredInv * RijInv * ( CosThetai * CosThetai - Third )

              if ( OptPressure ) then
                Tmp = 2._RK * CosThetai
                CosAux = 5._RK *  CosThetai * CosThetai - 1._RK
                Epsilon2 = Epsilon * RijSquaredInv * RijSquaredInv
                FXij = Epsilon2 * ( CosAux * eX - Tmp * OXi )               ! Kraft auf die Punktladung, sprich F2
                FYij = Epsilon2 * ( CosAux * eY - Tmp * OYi )
                FZij = Epsilon2 * ( CosAux * eZ - Tmp * OZi )
              end if
            end if
            EPot(j) = EPot(j) + EPotLocal
            if ( OptPressure ) then
              Virial(j) = Virial(j) - Third * ( FXij * PXij + FYij * PYij + FZij * PZij )
            end if
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijSquaredInv
            d2EpotdV2(j) = d2EpotdV2(j) + EPotLocal*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijSquaredInv)*Third*Third !xxxx7 QC
          end do
        end do
        do s2 = 1, this%N2Dipole
          pqd => this%PotQuadrupoleDipole(s1, s2)
          Epsilon = pqd%Epsilon
          RShieldSquared = pqd%RShieldSquared

          ! Assign pointers to site positions
          RX1 => pqd%Site1%RX
          RY1 => pqd%Site1%RY
          RZ1 => pqd%Site1%RZ
          OX1 => pqd%Site1%OX
          OY1 => pqd%Site1%OY
          OZ1 => pqd%Site1%OZ
          RX2 => pqd%Site2%RX
          RY2 => pqd%Site2%RY
          RZ2 => pqd%Site2%RZ
          OX2 => pqd%Site2%OX
          OY2 => pqd%Site2%OY
          OZ2 => pqd%Site2%OZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)
          OXi = OX1(np)
          OYi = OY1(np)
          OZi = OZ1(np)

          ! Loop over molecules
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
              CosAux = 5._RK * CosThetai**2 - 1._RK
              CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
              Rij4Inv = Epsilon / RijSquared**2
              EPotLocal = Rij4Inv * ( CosThetaj * CosAux - CosGammaij * CosThetai )

              if ( OptPressure ) then
                dCosThetai = Rij4Inv * (10._RK * CosThetai * CosThetaj - CosGammaij)
                dCosThetaj = Rij4Inv * CosAux
                dCosGammaij = -2._RK * Rij4Inv * CosThetai
                Tmp = -4._RK * RijInv * EPotLocal
                FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                          + (eX * CosThetaj - OXj) * dCosThetaj)
                FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                          + (eY * CosThetaj - OYj) * dCosThetaj)
                FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                          + (eZ * CosThetaj - OZj) * dCosThetaj)
              end if
            end if
            EPot(j) = EPot(j) + EPotLocal
            if ( OptPressure ) then
              Virial(j) = Virial(j) + Third * (FXij * PXij + FYij * PYij + FZij * PZij)
            end if
            RijInv2  =  RijInv*RijInv
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv2
            d2EpotdV2(j) = d2EpotdV2(j) + EPotLocal*(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv2)*Third*Third !xxxx8 QD
          end do
        end do

        do s2 = 1, this%N2Quadrupole
          pqq => this%PotQuadrupoleQuadrupole(s1, s2)
          Epsilon = pqq%Epsilon
          RShieldSquared = pqq%RShieldSquared

          ! Assign pointers to site positions
          RX1 => pqq%Site1%RX
          RY1 => pqq%Site1%RY
          RZ1 => pqq%Site1%RZ
          OX1 => pqq%Site1%OX
          OY1 => pqq%Site1%OY
          OZ1 => pqq%Site1%OZ
          RX2 => pqq%Site2%RX
          RY2 => pqq%Site2%RY
          RZ2 => pqq%Site2%RZ
          OX2 => pqq%Site2%OX
          OY2 => pqq%Site2%OY
          OZ2 => pqq%Site2%OZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)
          OXi = OX1(np)
          OYi = OY1(np)
          OZi = OZ1(np)

          ! Loop over molecules
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
              EPotLocal = Rij5Inv * (1._RK - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&                         - 15._RK * CosThetaiSquared * CosThetajSquared + 2._RK * Tmp**2)

              if ( OptPressure ) then
                dCosThetai = Rij5Inv * (-10._RK * CosThetai &
&                                      - 30._RK * CosThetai * CosThetajSquared &
&                                      - 20._RK * CosThetaj * Tmp)
                dCosThetaj = Rij5Inv * (-10._RK * CosThetaj &
&                                      - 30._RK * CosThetaj * CosThetaiSquared &
&                                      - 20._RK * CosThetai * Tmp)
                dCosGammaij = 4._RK * Rij5Inv * Tmp
                Tmp = -5._RK * RijInv * EPotLocal
                FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                          + (eX * CosThetaj - OXj) * dCosThetaj)
                FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                          + (eY * CosThetaj - OYj) * dCosThetaj)
                FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                          + (eZ * CosThetaj - OZj) * dCosThetaj)
              end if
            end if
            EPot(j) = EPot(j) + EPotLocal
            if ( OptPressure ) then
              Virial(j) = Virial(j) + Third * ( FXij * PXij + FYij * PYij + FZij * PZij )
            end if
            RijInv2  =  RijInv*RijInv
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv2
            d2EpotdV2(j) = d2EpotdV2(j) + EPotLocal*(35._RK*sitecorr*sitecorr-5._RK*Plen2*RijInv2)/9._RK !xxxx9 QQ
          end do
        end do
      end do

      ! Explicit reaction field contribution
      if ( (this%ReactionField) .or. (LongRange .eq. ExtRField) ) then
        if ( LongRange .eq. RField) then    ! Normal ReactionField
          MueX2 => this%MueX2
          MueY2 => this%MueY2
          MueZ2 => this%MueZ2

          mueXi = this%MueX1(np)
          mueYi = this%MueY1(np)
          mueZi = this%MueZ1(np)

          do k = 1, this%NInCutoff(np)
            j = this%CutoffPartner(k, np)
            EPot(j) = EPot(j) + this%RFConst2 * ( mueXi * MueX2(j) + mueYi * MueY2(j) + mueZi * MueZ2(j) )
          end do

        else         ! Extended ReactionField
          if ( ((this%N1Charge > 1) .and. (this%N2Charge > 1) ) .or. (this%N1Charge+this%N2Charge .eq. 0)) then 
            MueX2 => this%MueX2
            MueY2 => this%MueY2
            MueZ2 => this%MueZ2

            mueXi = this%MueX1(np)
            mueYi = this%MueY1(np)
            mueZi = this%MueZ1(np)

            do k = 1, this%NInCutoff(np)
              j = this%CutoffPartner(k, np)
              EPot(j) = EPot(j) + this%RFConst2 * ( mueXi * MueX2(j) + mueYi * MueY2(j) + mueZi * MueZ2(j) )
            end do

          else if ( (this%N1Charge .eq. 1) .and. (this%N2Charge .ne. 1) ) then 
          ! Assign pointers to site positions
            if (this%N2Charge > 0) then
              pcc => this%PotChargeCharge(1,1)
              PX2 => pcc%Site2%PX
              PY2 => pcc%Site2%PY
              PZ2 => pcc%Site2%PZ
              RX1 => pcc%Site1%RX
              RY1 => pcc%Site1%RY
              RZ1 => pcc%Site1%RZ
            else
              pcd => this%PotChargeDipole(1,1)
              PX2 => pcd%Site2%PX
              PY2 => pcd%Site2%PY
              PZ2 => pcd%Site2%PZ
              RX1 => pcd%Site1%RX
              RY1 => pcd%Site1%RY
              RZ1 => pcd%Site1%RZ
            end if
!!!!!!!!!!!!!!!!!!!!!!!!!!
            muexi = 0.0_RK
            mueyi = 0.0_RK
            muezi = 0.0_RK
            q = this%lad1
            MueX2 => this%MueX2
            MueY2 => this%MueY2
            MueZ2 => this%MueZ2
            do k = 1, this%NInCutoff(np)
              j = this%CutoffPartner(k, np)
!!!!!!!!!!!!!!!!!!!!!!!!!
              RXij = PX2(j)-RX1(np)
              RYij = PY2(j)-RY1(np)
              RZij = PZ2(j)-RZ1(np)
              RXij = (RXij - anint(RXij))*BoxLength
              RYij = (RYij - anint(RYij))*BoxLength
              RZij = (RZij - anint(RZij))*BoxLength
              EPot(j) = EPot(j) -this%RFConst2 * q * ( RXij*MueX2(j) + RYij*MueY2(j) + RZij*MueZ2(j) )

            end do

          else if ( (this%N1Charge > 1) .and. (this%N2Charge .eq. 1) ) then 
          ! Assign pointers to site positions
           if (this%N1Charge > 0) then
            pcc => this%PotChargeCharge(1,1)
            PX1 => pcc%Site1%PX
            PY1 => pcc%Site1%PY
            PZ1 => pcc%Site1%PZ
            RX2 => pcc%Site2%RX
            RY2 => pcc%Site2%RY
            RZ2 => pcc%Site2%RZ
           else
            pdc => this%PotDipoleCharge(1,1)
            PX1 => pdc%Site1%PX
            PY1 => pdc%Site1%PY
            PZ1 => pdc%Site1%PZ
            RX2 => pdc%Site2%RX
            RY2 => pdc%Site2%RY
            RZ2 => pdc%Site2%RZ
           end if
            q = this%lad2

            muexi = 0.0_RK
            mueyi = 0.0_RK
            muezi = 0.0_RK

            do k = 1, this%NInCutoff(np)
              j = this%CutoffPartner(k, np)
              RXij = RX2(j)-PX1(np)
              RYij = RY2(j)-PY1(np)
              RZij = RZ2(j)-PZ1(np)
              RXij = (RXij - anint(RXij))*BoxLength
              RYij = (RYij - anint(RYij))*BoxLength
              RZij = (RZij - anint(RZij))*BoxLength
              muexi = (RXij)*q
              mueyi = (RYij)*q
              muezi = (RZij)*q
              EPot(j) = EPot(j) +this%RFConst2 * ( muexi * this%MueX1(np) + mueyi * this%MueY1(np) + muezi * this%MueZ1(np) )
            end do

          else if ( (this%N1Charge .eq. 1) .and. (this%N2Charge .eq. 1) ) then 
            pcc => this%PotChargeCharge(1, 1)
            Epsilon = pcc%Epsilon
            RShieldSquared = pcc%RShieldSquared

          ! Assign pointers to site positions
            RX1 => pcc%Site1%RX
            RY1 => pcc%Site1%RY
            RZ1 => pcc%Site1%RZ
            RX2 => pcc%Site2%RX
            RY2 => pcc%Site2%RY
            RZ2 => pcc%Site2%RZ
            do k = 1, this%NInCutoff(np)
              j = this%CutoffPartner(k, np)
              RXij = RX2(j)-RX1(np)
              RYij = RY2(j)-RY1(np)
              RZij = RZ2(j)-RZ1(np)
              RXij = (RXij - anint(RXij))*BoxLength
              RYij = (RYij - anint(RYij))*BoxLength
              RZij = (RZij - anint(RZij))*BoxLength
              Rij = (RXij**2+RYij**2+RZij**2)
            end do
          end if
        end if 
      end if 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else ! Site-site cutoff

      ! Calculate Lennard-Jones energy
      do s1 = 1, this%N1LJ126
        do s2 = 1, this%N2LJ126

          ! Set site specific variables
          plj => this%PotLJ126LJ126(s1, s2)
          SigmaSquared = plj%SigmaSquared
          Epsilon4 = plj%Epsilon4
          if ( OptPressure ) then
            Epsilon48 = plj%Epsilon48
          end if

          ! Assign pointers to site positions
          RX1 => plj%Site1%RX
          RY1 => plj%Site1%RY
          RZ1 => plj%Site1%RZ
          RX2 => plj%Site2%RX
          RY2 => plj%Site2%RY
          RZ2 => plj%Site2%RZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)

          ! Loop over molecules
#if MPI_VER > 0
!CDIR NODEP
          do j = this%NPart20, this%NPart22
#else
!CDIR NODEP
          do j = 1, N
#endif
            if( this%SameComponent .and. j == np ) cycle
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
            if( RijSquared >= RCutoffSquared ) cycle
            RijSquaredInv = SigmaSquared / RijSquared
            Rij6Inv = RijSquaredInv**3
            EPot(j) = EPot(j) + Epsilon4 * Rij6Inv * (Rij6Inv - 1._RK)
            if ( OptPressure ) then
              Fij = Epsilon48 * Rij6Inv * (Rij6Inv - .5_RK) * RijSquaredInv
              FXij = Fij * RXij
              FYij = Fij * RYij
              FZij = Fij * RZij
              Virial(j) = Virial(j) + BoxLengthThird * (PXij * FXij + PYij * FYij + PZij * FZij)
            end if
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)/RijSquared
            d2EpotdV2(j) = d2EpotdV2(j) + Epsilon4 * Rij6Inv *(12._RK *Rij6Inv -  6._RK) * (sitecorr * sitecorr - Plen2/RijSquared)*Third*Third !xxxx2 LJ
            d2EpotdV2(j) = d2EpotdV2(j) + Epsilon4 * Rij6Inv *(156._RK*Rij6Inv - 42._RK) *  sitecorr * sitecorr *Third*Third
          end do
        end do
      end do

      ! No point charges allowed with site-site cutoff

      ! Calculate dipolar energy
      do s1 = 1, this%N1Dipole
        do s2 = 1, this%N2Dipole
          pdd => this%PotDipoleDipole(s1, s2)
          Epsilon = pdd%Epsilon
          RShieldSquared = pdd%RShieldSquared
          RFConst2 = Epsilon * this%RFConst2

          ! Assign pointers to site positions
          RX1 => pdd%Site1%RX
          RY1 => pdd%Site1%RY
          RZ1 => pdd%Site1%RZ
          OX1 => pdd%Site1%OX
          OY1 => pdd%Site1%OY
          OZ1 => pdd%Site1%OZ
          RX2 => pdd%Site2%RX
          RY2 => pdd%Site2%RY
          RZ2 => pdd%Site2%RZ
          OX2 => pdd%Site2%OX
          OY2 => pdd%Site2%OY
          OZ2 => pdd%Site2%OZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)
          OXi = OX1(np)
          OYi = OY1(np)
          OZi = OZ1(np)

          ! Loop over molecules
#if MPI_VER > 0
!CDIR NODEP
          do j = this%NPart20, this%NPart22
#else
!CDIR NODEP
          do j = 1, N
#endif
            if( this%SameComponent .and. j == np ) cycle
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

            if( RijSquared >= RCutoffSquared ) cycle
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
              EPotLocal = Rij3Inv * Tmp + RFConst2 * CosGammaij
              if ( OptPressure ) then
                Rij4Inv3 = 3._RK * Rij3Inv * RijInv
                FXij = Rij4Inv3 * (eX * Tmp - (eX * CosThetai - OXi) * CosThetaj &
&                                           - (eX * CosThetaj - OXj) * CosThetai)
                FYij = Rij4Inv3 * (eY * Tmp - (eY * CosThetai - OYi) * CosThetaj &
&                                           - (eY * CosThetaj - OYj) * CosThetai)
                FZij = Rij4Inv3 * (eZ * Tmp - (eZ * CosThetai - OZi) * CosThetaj &
&                                           - (eZ * CosThetaj - OZj) * CosThetai)
              end if
            end if
            EPot(j) = EPot(j) + EPotLocal
            if ( OptPressure ) then
              Virial(j) = Virial(j) + Third * ( FXij * PXij + FYij * PYij + FZij * PZij )
            end if
            RijInv2  =  RijInv*RijInv
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv2
            d2EpotdV2(j) = d2EpotdV2(j) + EPotLocal*(15._RK*sitecorr*sitecorr-3._RK*Plen2*RijInv2)/9._RK !xxxxss5 DD
          end do
        end do

        do s2 = 1, this%N2Quadrupole
          pdq => this%PotDipoleQuadrupole(s1, s2)
          Epsilon = pdq%Epsilon
          RShieldSquared = pdq%RShieldSquared

          ! Assign pointers to site positions
          RX1 => pdq%Site1%RX
          RY1 => pdq%Site1%RY
          RZ1 => pdq%Site1%RZ
          OX1 => pdq%Site1%OX
          OY1 => pdq%Site1%OY
          OZ1 => pdq%Site1%OZ
          RX2 => pdq%Site2%RX
          RY2 => pdq%Site2%RY
          RZ2 => pdq%Site2%RZ
          OX2 => pdq%Site2%OX
          OY2 => pdq%Site2%OY
          OZ2 => pdq%Site2%OZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)
          OXi = OX1(np)
          OYi = OY1(np)
          OZi = OZ1(np)

          ! Loop over molecules
#if MPI_VER > 0
!CDIR NODEP
          do j = this%NPart20, this%NPart22
#else
!CDIR NODEP
          do j = 1, N
#endif
            if( this%SameComponent .and. j == np ) cycle
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

            if( RijSquared >= RCutoffSquared ) cycle
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
              CosAux = 1._RK - 5._RK * CosThetaj**2
              CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
              Rij4Inv = Epsilon / RijSquared**2
              EPotLocal = Rij4Inv * ( CosGammaij * CosThetaj + CosThetai * CosAux )

              if ( OptPressure ) then
                dCosThetai = Rij4Inv * CosAux
                dCosThetaj = Rij4Inv * (CosGammaij - 10._RK * CosThetai * CosThetaj)
                dCosGammaij = 2._RK * Rij4Inv * CosThetaj
                Tmp = -4._RK * RijInv * EPotLocal
                FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                          + (eX * CosThetaj - OXj) * dCosThetaj)
                FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                          + (eY * CosThetaj - OYj) * dCosThetaj)
                FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                          + (eZ * CosThetaj - OZj) * dCosThetaj)
              end if
            end if
            EPot(j) = EPot(j) + EPotLocal
            if ( OptPressure ) then
              Virial(j) = Virial(j) + Third * ( FXij * PXij + FYij * PYij + FZij * PZij )
            end if
            RijInv2  =  RijInv*RijInv
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv2
            d2EpotdV2(j) = d2EpotdV2(j) + EPotLocal*(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv2)/9._RK !xxxxss6 DQ
          end do
        end do
      end do

      ! Calculate quadrupolar energy
      do s1 = 1, this%N1Quadrupole
        do s2 = 1, this%N2Dipole
          pqd => this%PotQuadrupoleDipole(s1, s2)
          Epsilon = pqd%Epsilon
          RShieldSquared = pqd%RShieldSquared

          ! Assign pointers to site positions
          RX1 => pqd%Site1%RX
          RY1 => pqd%Site1%RY
          RZ1 => pqd%Site1%RZ
          OX1 => pqd%Site1%OX
          OY1 => pqd%Site1%OY
          OZ1 => pqd%Site1%OZ
          RX2 => pqd%Site2%RX
          RY2 => pqd%Site2%RY
          RZ2 => pqd%Site2%RZ
          OX2 => pqd%Site2%OX
          OY2 => pqd%Site2%OY
          OZ2 => pqd%Site2%OZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)
          OXi = OX1(np)
          OYi = OY1(np)
          OZi = OZ1(np)

          ! Loop over molecules
#if MPI_VER > 0
!CDIR NODEP
          do j = this%NPart20, this%NPart22
#else
!CDIR NODEP
          do j = 1, N
#endif
            if( this%SameComponent .and. j == np ) cycle
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

            if( RijSquared >= RCutoffSquared ) cycle
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
              CosAux = 5._RK * CosThetai**2 - 1._RK
              CosGammaij = 2._RK * (OXi * OXj + OYi * OYj + OZi * OZj)
              Rij4Inv = Epsilon / RijSquared**2
              EPotLocal = Rij4Inv * ( CosThetaj * CosAux &
&                                   - CosGammaij * CosThetai )
              if ( OptPressure ) then
                dCosThetai = Rij4Inv * (10._RK * CosThetai * CosThetaj - CosGammaij)
                dCosThetaj = Rij4Inv * CosAux
                dCosGammaij = -2._RK * Rij4Inv * CosThetai
                Tmp = -4._RK * RijInv * EPotLocal
                FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                          + (eX * CosThetaj - OXj) * dCosThetaj)
                FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                          + (eY * CosThetaj - OYj) * dCosThetaj)
                FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                          + (eZ * CosThetaj - OZj) * dCosThetaj)
              end if
            end if
            EPot(j) = EPot(j) + EPotLocal
            if ( OptPressure ) then
              Virial(j) = Virial(j) + Third * (FXij * PXij + FYij * PYij + FZij * PZij)
            end if
            RijInv2  =  RijInv*RijInv
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv2
            d2EpotdV2(j) = d2EpotdV2(j) + EPotLocal*(24._RK*sitecorr*sitecorr-4._RK*Plen2*RijInv2)/9._RK !xxxxss8 QD
          end do
        end do
        do s2 = 1, this%N2Quadrupole
          pqq => this%PotQuadrupoleQuadrupole(s1, s2)
          Epsilon = pqq%Epsilon
          RShieldSquared = pqq%RShieldSquared

          ! Assign pointers to site positions
          RX1 => pqq%Site1%RX
          RY1 => pqq%Site1%RY
          RZ1 => pqq%Site1%RZ
          OX1 => pqq%Site1%OX
          OY1 => pqq%Site1%OY
          OZ1 => pqq%Site1%OZ
          RX2 => pqq%Site2%RX
          RY2 => pqq%Site2%RY
          RZ2 => pqq%Site2%RZ
          OX2 => pqq%Site2%OX
          OY2 => pqq%Site2%OY
          OZ2 => pqq%Site2%OZ

          RXi = RX1(np)
          RYi = RY1(np)
          RZi = RZ1(np)
          OXi = OX1(np)
          OYi = OY1(np)
          OZi = OZ1(np)

          ! Loop over molecules
#if MPI_VER > 0
!CDIR NODEP
          do j = this%NPart20, this%NPart22
#else
!CDIR NODEP
          do j = 1, N
#endif
            if( this%SameComponent .and. j == np ) cycle
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

            if( RijSquared >= RCutoffSquared ) cycle
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
              EPotLocal = Rij5Inv * (1._RK - 5._RK * (CosThetaiSquared + CosThetajSquared) &
&                         - 15._RK * CosThetaiSquared * CosThetajSquared + 2._RK * Tmp**2)

              if ( OptPressure ) then
                dCosThetai = Rij5Inv * (-10._RK * CosThetai &
&                                      - 30._RK * CosThetai * CosThetajSquared &
&                                      - 20._RK * CosThetaj * Tmp)
                dCosThetaj = Rij5Inv * (-10._RK * CosThetaj &
&                                      - 30._RK * CosThetaj * CosThetaiSquared &
&                                      - 20._RK * CosThetai * Tmp)

                dCosGammaij = 4._RK * Rij5Inv * Tmp
                Tmp = -5._RK * RijInv * EPotLocal

                FXij = -eX * Tmp + RijInv * ((eX * CosThetai - OXi) * dCosThetai &
&                                          + (eX * CosThetaj - OXj) * dCosThetaj)
                FYij = -eY * Tmp + RijInv * ((eY * CosThetai - OYi) * dCosThetai &
&                                          + (eY * CosThetaj - OYj) * dCosThetaj)
                FZij = -eZ * Tmp + RijInv * ((eZ * CosThetai - OZi) * dCosThetai &
&                                          + (eZ * CosThetaj - OZj) * dCosThetaj)
              end if
            end if
            EPot(j) = EPot(j) + EPotLocal
            if ( OptPressure ) then
              Virial(j) = Virial(j) + Third * ( FXij * PXij + FYij * PYij + FZij * PZij )
            end if
            RijInv2  =  RijInv*RijInv
            Plen2    =  PXij*PXij+PYij*PYij+PZij*PZij
            sitecorr = (PXij*RXij+PYij*RYij+PZij*RZij)*RijInv2
            d2EpotdV2(j) = d2EpotdV2(j) + EPotLocal*(35._RK*sitecorr*sitecorr-5._RK*Plen2*RijInv2)/9._RK !xxxxss9 QQ
          end do
        end do
      end do

    end if

    this%EPot1 = EPot
end subroutine TInteraction_Energy



!==============================================================!
!  Subroutine TInteraction_UpdateBoxLength                     !
!==============================================================!

  subroutine TInteraction_UpdateBoxLength( this, BoxLength )

    implicit none

    ! Declare arguments
    type(TInteraction)   :: this
    real(RK), intent(in) :: BoxLength

    ! Declare local variables
    integer :: i, j

    ! Update constants
    this%RCutoffSquaredScaled = this%RCutoffSquared / BoxLength**2

    ! Update BoxLength in Potentials
    do i = 1, this%N1LJ126
      do j = 1, this%N2LJ126
        call UpdateBoxLength( this%PotLJ126LJ126( i, j ), BoxLength )
      end do
    end do

  end subroutine TInteraction_UpdateBoxLength



!==============================================================!
!  Subroutine TInteraction_CalcPartners                        !
!==============================================================!

  subroutine TInteraction_CalcPartners( this )

    implicit none

    ! Declare arguments
    type(TInteraction) :: this

    ! Declare local variables
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: PXi, PYi, PZi, PXij, PYij, PZij
    real(RK)          :: RijSquared
    real(RK)          :: RCutoff
    integer           :: i, j, N, N2, NInCutoff
 
    ! Set cutoff radius
    RCutoff = this%RCutoffSquaredScaled
    N = this%NPart1
    this%NInCutoff(:) = 0

!$OMP PARALLEL PRIVATE(PX1, PY1, PZ1, PX2, PY2, PZ2, i, j ,NInCutoff, N2, RijSquared,PXi, PYi, PZi, PXij, PYij, PZij)
    ! Assign local pointers
    PX1 => this%PX1
    PY1 => this%PY1
    PZ1 => this%PZ1
    PX2 => this%PX2
    PY2 => this%PY2
    PZ2 => this%PZ2

!    write(*,*) "-- after PX1 ",LOC(this%PX1)
!    write(*,*) "-- after PY1 ",LOC(this%PY1)
!    write(*,*) "-- after PZ1 ",LOC(this%PZ1)
!    write(*,*) "-- after PX2 ",LOC(this%PX2)
!    write(*,*) "-- after PY2 ",LOC(this%PY2)
!    write(*,*) "-- after PZ2 ",LOC(this%PZ2)
    
    ! Calculate partners within cutoff sphere
    if( this%SameComponent ) then
#if MPI_VER > 0
      if( this%NPart10 <= (N+1)/2 ) then
        if( this%NPart12 > (N+1)/2 ) then
!$OMP DO        
          do i = this%NPart10, (N+1) / 2
#else
!$OMP DO        
      do i = 1, (N+1) / 2
#endif
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        NInCutoff = this%NInCutoff(i)
        do j = i + 1, (N/2) + i
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = PXij - anint( PXij )
          PYij = PYij - anint( PYij )
          PZij = PZij - anint( PZij )
          RijSquared = PXij**2 + PYij**2 + PZij**2

          if( RijSquared < RCutoff ) then
            NInCutoff = NInCutoff + 1
            this%CutoffPartner(NInCutoff, i) = j
          end if
        end do
        this%NInCutoff(i) = NInCutoff
      end do
!$OMP END DO      

!$OMP DO
#if MPI_VER > 0
          do i = (N+1) / 2 + 1, this%NPart12
#else
      do i = (N+1) / 2 + 1, N
#endif
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        NInCutoff = this%NInCutoff(i)
        do j = 1, i - N/2 - 1
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = PXij - anint( PXij )
          PYij = PYij - anint( PYij )
          PZij = PZij - anint( PZij )
          RijSquared = PXij**2 + PYij**2 + PZij**2

          if( RijSquared < RCutoff ) then
            NInCutoff = NInCutoff + 1
            this%CutoffPartner(NInCutoff, i) = j
          end if
        end do
        do j = i + 1, N
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = PXij - anint( PXij )
          PYij = PYij - anint( PYij )
          PZij = PZij - anint( PZij )
          RijSquared = PXij**2 + PYij**2 + PZij**2

          if( RijSquared < RCutoff ) then
            NInCutoff = NInCutoff + 1
            this%CutoffPartner(NInCutoff, i) = j
          end if
        end do
        this%NInCutoff(i) = NInCutoff
      end do
!$OMP END DO      

#if MPI_VER > 0
        else
!$OMP DO         
          do i = this%NPart10, this%NPart12
            PXi = PX1(i)
            PYi = PY1(i)
            PZi = PZ1(i)
            NInCutoff = this%NInCutoff(i)
            do j = i + 1, N/2 + i
              PXij = PXi - PX2(j)
              PYij = PYi - PY2(j)
              PZij = PZi - PZ2(j)
              PXij = PXij - anint( PXij )
              PYij = PYij - anint( PYij )
              PZij = PZij - anint( PZij )
              RijSquared = PXij**2 + PYij**2 + PZij**2

              if( RijSquared < RCutoff ) then
                NInCutoff = NInCutoff + 1
                this%CutoffPartner(NInCutoff, i) = j
              end if
            end do
            this%NInCutoff(i) = NInCutoff
          end do
!$OMP END DO          
        end if

      else
!$OMP DO       
        do i = this%NPart10, this%NPart12
          PXi = PX1(i)
          PYi = PY1(i)
          PZi = PZ1(i)
          NInCutoff = this%NInCutoff(i)
          do j = 1, i - N/2 - 1
            PXij = PXi - PX2(j)
            PYij = PYi - PY2(j)
            PZij = PZi - PZ2(j)
            PXij = PXij - anint( PXij )
            PYij = PYij - anint( PYij )
            PZij = PZij - anint( PZij )
            RijSquared = PXij**2 + PYij**2 + PZij**2

            if( RijSquared < RCutoff ) then
              NInCutoff = NInCutoff + 1
              this%CutoffPartner(NInCutoff, i) = j
            end if
          end do
          do j = i + 1, N
            PXij = PXi - PX2(j)
            PYij = PYi - PY2(j)
            PZij = PZi - PZ2(j)
            PXij = PXij - anint( PXij )
            PYij = PYij - anint( PYij )
            PZij = PZij - anint( PZij )
            RijSquared = PXij**2 + PYij**2 + PZij**2
            if( RijSquared < RCutoff ) then
              NInCutoff = NInCutoff + 1
              this%CutoffPartner(NInCutoff, i) = j
            end if
          end do
          this%NInCutoff(i) = NInCutoff
        end do
!$OMP END DO        
        
      end if
#endif
    else
      N2 = this%NPart2

!$OMP DO      
#if MPI_VER > 0
      do i = this%NPart10, this%NPart12
#else
      do i = 1, N
#endif
        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        NInCutoff = this%NInCutoff(i)
        do j = 1, N2
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = PXij - anint( PXij )
          PYij = PYij - anint( PYij )
          PZij = PZij - anint( PZij )
          RijSquared = PXij**2 + PYij**2 + PZij**2

          if( RijSquared < RCutoff ) then
            NInCutoff = NInCutoff + 1
            this%CutoffPartner(NInCutoff, i) = j
          end if
        end do
        this%NInCutoff(i) = NInCutoff
      end do
!$OMP END DO      
    end if
!$OMP END PARALLEL
  end subroutine TInteraction_CalcPartners


!==============================================================!
!  Subroutine TInteraction_CalcPartners                        !
!==============================================================!

  subroutine TInteraction_CalcPartnersRDF( this )

    implicit none

    ! Declare arguments
    type(TInteraction) :: this

    ! Declare local variables
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: PXi, PYi, PZi, PXij, PYij, PZij
    real(RK)          :: RijSquared
    real(RK)          :: RCutoff
    integer           :: i, j, N, N2, NInCutoff

    ! Set cutoff radius
    RCutoff = this%RCutoffSquaredScaled
    N = this%NPart1
    this%NInCutoff(:) = 0

    ! Assign local pointers
    PX1 => this%PX1
    PY1 => this%PY1
    PZ1 => this%PZ1
    PX2 => this%PX2
    PY2 => this%PY2
    PZ2 => this%PZ2

    ! Calculate partners within cutoff sphere
    if( this%SameComponent ) then

      do i = 1, (N+1) / 2

        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        NInCutoff = this%NInCutoff(i)

        do j = i + 1, (N/2) + i
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = PXij - anint( PXij )
          PYij = PYij - anint( PYij )
          PZij = PZij - anint( PZij )
          RijSquared = PXij**2 + PYij**2 + PZij**2

          if( RijSquared < RCutoff ) then
            NInCutoff = NInCutoff + 1
            this%CutoffPartner(NInCutoff, i) = j
          end if
        end do
        this%NInCutoff(i) = NInCutoff
      end do

      do i = (N+1) / 2 + 1, N

        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        NInCutoff = this%NInCutoff(i)
        do j = 1, i - N/2 - 1
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = PXij - anint( PXij )
          PYij = PYij - anint( PYij )
          PZij = PZij - anint( PZij )
          RijSquared = PXij**2 + PYij**2 + PZij**2

          if( RijSquared < RCutoff ) then
            NInCutoff = NInCutoff + 1
            this%CutoffPartner(NInCutoff, i) = j
          end if
        end do

        do j = i + 1, N
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = PXij - anint( PXij )
          PYij = PYij - anint( PYij )
          PZij = PZij - anint( PZij )
          RijSquared = PXij**2 + PYij**2 + PZij**2
          if( RijSquared < RCutoff ) then
            NInCutoff = NInCutoff + 1
            this%CutoffPartner(NInCutoff, i) = j
          end if
        end do
        this%NInCutoff(i) = NInCutoff
      end do

    else
      N2 = this%NPart2

      do i = 1, N

        PXi = PX1(i)
        PYi = PY1(i)
        PZi = PZ1(i)
        NInCutoff = this%NInCutoff(i)
        do j = 1, N2
          PXij = PXi - PX2(j)
          PYij = PYi - PY2(j)
          PZij = PZi - PZ2(j)
          PXij = PXij - anint( PXij )
          PYij = PYij - anint( PYij )
          PZij = PZij - anint( PZij )
          RijSquared = PXij**2 + PYij**2 + PZij**2
          if( RijSquared < RCutoff ) then
            NInCutoff = NInCutoff + 1
            this%CutoffPartner(NInCutoff, i) = j
          end if
        end do
        this%NInCutoff(i) = NInCutoff
      end do
    end if

  end subroutine TInteraction_CalcPartnersRDF



!==============================================================!
!  Subroutine TInteraction_CalcPartners1                       !
!==============================================================!

  subroutine TInteraction_CalcPartners1( this, np )

    implicit none

    ! Declare arguments
    type(TInteraction)  :: this
    integer, intent(in) :: np

    ! Declare local variables
    real(RK), pointer :: PX2(:), PY2(:), PZ2(:)
    real(RK)          :: PXi, PYi, PZi, PXij, PYij, PZij
    real(RK)          :: RijSquared
    real(RK)          :: RCutoffSquaredScaled
    integer           :: j, NInCutoff

    ! Set cutoff radius
    RCutoffSquaredScaled = this%RCutoffSquaredScaled

    ! Assign local pointers
    PX2 => this%PX2
    PY2 => this%PY2
    PZ2 => this%PZ2

    ! Calculate partners within cutoff sphere
    PXi = this%PX1(np)
    PYi = this%PY1(np)
    PZi = this%PZ1(np)
    NInCutoff = 0
#if MPI_VER > 0
    do j = this%NPart20, this%NPart22
#else
    do j = 1, this%NPart2
#endif
      if( this%SameComponent .and. j == np ) cycle
      PXij = PXi - PX2(j)
      PYij = PYi - PY2(j)
      PZij = PZi - PZ2(j)
      PXij = PXij - anint( PXij )
      PYij = PYij - anint( PYij )
      PZij = PZij - anint( PZij )
      RijSquared = PXij**2 + PYij**2 + PZij**2

      if( RijSquared < RCutoffSquaredScaled ) then
        NInCutoff = NInCutoff + 1
        this%CutoffPartner(NInCutoff, np) = j
      end if
    end do
    this%NInCutoff(np) = NInCutoff

  end subroutine TInteraction_CalcPartners1



!==============================================================!
!  Subroutine TInteraction_CalcPartnersTest                    !
!==============================================================!

  subroutine TInteraction_CalcPartnersTest( this )

    implicit none

    ! Declare arguments
    type(TInteraction) :: this

    ! Declare local variables
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)
    real(RK)          :: PXi, PYi, PZi, PXij, PYij, PZij
    real(RK)          :: RijSquared
    real(RK)          :: RCutoff
    integer           :: i, j, NInCutoff

    ! Set cutoff radius
    RCutoff = this%RCutoffSquaredScaled

    ! Assign local pointers
    PX1 => this%PX1Test
    PY1 => this%PY1Test
    PZ1 => this%PZ1Test
    PX2 => this%PX2
    PY2 => this%PY2
    PZ2 => this%PZ2

    ! Calculate partners within cutoff sphere
    do i = 1, this%NTest1
      PXi = PX1(i)
      PYi = PY1(i)
      PZi = PZ1(i)
      NInCutoff = 0

      do j = 1, this%NPart2
        PXij = PXi - PX2(j)
        PYij = PYi - PY2(j)
        PZij = PZi - PZ2(j)
        PXij = PXij - anint( PXij )
        PYij = PYij - anint( PYij )
        PZij = PZij - anint( PZij )
        RijSquared = PXij**2 + PYij**2 + PZij**2

        if( RijSquared < RCutoff ) then
          NInCutoff = NInCutoff + 1
          this%CutoffPartner(NInCutoff, i) = j
        end if
      end do
      this%NInCutoff(i) = NInCutoff
    end do

  end subroutine TInteraction_CalcPartnersTest



end module ms2_interaction
