!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 1.0                !
!  (c) 2011 by TU Kaiserslautern                               !
!      P.O. Box 67653                                          !
!      67653 Kaiserslautern                                    !
!==============================================================!
!  Module ms2_site                                             !
!  Contains TSite* objects                                     !
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
!DEC$ MESSAGE:'Compiling ms2_site.F90...'
#endif

module ms2_site

  use ms2_global



!==============================================================!
!  Type TSiteLJ126                                             !
!==============================================================!

  type TSiteLJ126

    integer           :: SiteId
    integer           :: UnitNumber
    real(RK),pointer  :: r(:)
    real(RK)          :: sig, eps
    real(RK)          :: mass
    integer, pointer  :: NPartMax, NPart, NTest
    integer, pointer  :: NPart0, NPart1, NPart2
    real(RK), pointer :: RX(:), RY(:), RZ(:)
    real(RK), pointer :: FX(:), FY(:), FZ(:)
    real(RK), pointer :: PX(:), PY(:), PZ(:)
    real(RK), pointer :: RXTest(:), RYTest(:), RZTest(:)
    real(RK), pointer :: PXTest(:), PYTest(:), PZTest(:)

  end type TSiteLJ126

  interface Construct
    module procedure TSiteLJ126_Construct
  end interface

  interface Destruct
    module procedure TSiteLJ126_Destruct
  end interface

  interface Allocate
    module procedure TSiteLJ126_Allocate
  end interface

  interface Deallocate
    module procedure TSiteLJ126_Deallocate
  end interface

  interface Save
    module procedure TSiteLJ126_Save
  end interface



!==============================================================!
!  Type TSiteCharge                                            !
!==============================================================!

  type TSiteCharge

    integer           :: SiteId
    integer           :: UnitNumber
    real(RK),pointer  :: r(:)
    real(RK)          :: e
    real(RK)          :: mass
    real(RK)          :: shield
    integer, pointer  :: NPartMax, NPart, NTest
    integer, pointer  :: NPart0, NPart1, NPart2
    real(RK), pointer :: RX(:), RY(:), RZ(:)
    real(RK), pointer :: FX(:), FY(:), FZ(:)
    real(RK), pointer :: PX(:), PY(:), PZ(:)
    real(RK), pointer :: RXTest(:), RYTest(:), RZTest(:)
    real(RK), pointer :: PXTest(:), PYTest(:), PZTest(:)

  end type TSiteCharge

  interface Construct
    module procedure TSiteCharge_Construct
  end interface

  interface Destruct
    module procedure TSiteCharge_Destruct
  end interface

  interface Allocate
    module procedure TSiteCharge_Allocate
  end interface

  interface Deallocate
    module procedure TSiteCharge_Deallocate
  end interface

  interface Save
    module procedure TSiteCharge_Save
  end interface



!==============================================================!
!  Type TSiteDipole                                            !
!==============================================================!

  type TSiteDipole

    integer           :: SiteId
    integer           :: UnitNumber
    real(RK),pointer  :: r(:), or(:)
    real(RK)          :: D
    real(RK)          :: mass
    real(RK)          :: shield
    integer, pointer  :: NPartMax, NPart, NTest
    integer, pointer  :: NPart0, NPart1, NPart2
    real(RK), pointer :: RX(:), RY(:), RZ(:)
    real(RK), pointer :: OX(:), OY(:), OZ(:)
    real(RK), pointer :: FX(:), FY(:), FZ(:)
    real(RK), pointer :: TX(:), TY(:), TZ(:)
    real(RK), pointer :: PX(:), PY(:), PZ(:)
    real(RK), pointer :: RXTest(:), RYTest(:), RZTest(:)
    real(RK), pointer :: OXTest(:), OYTest(:), OZTest(:)
    real(RK), pointer :: PXTest(:), PYTest(:), PZTest(:)

  end type TSiteDipole

  interface Construct
    module procedure TSiteDipole_Construct
  end interface

  interface Destruct
    module procedure TSiteDipole_Destruct
  end interface

  interface Allocate
    module procedure TSiteDipole_Allocate
  end interface

  interface Deallocate
    module procedure TSiteDipole_Deallocate
  end interface

  interface Save
    module procedure TSiteDipole_Save
  end interface



!==============================================================!
!  Type TSiteQuadrupole                                        !
!==============================================================!

  type TSiteQuadrupole

    integer           :: SiteId
    integer           :: UnitNumber
    real(RK),pointer  :: r(:), or(:)
    real(RK)          :: Q
    real(RK)          :: mass
    real(RK)          :: shield
    integer, pointer  :: NPartMax, NPart, NTest
    integer, pointer  :: NPart0, NPart1, NPart2
    real(RK), pointer :: RX(:), RY(:), RZ(:)
    real(RK), pointer :: OX(:), OY(:), OZ(:)
    real(RK), pointer :: FX(:), FY(:), FZ(:)
    real(RK), pointer :: TX(:), TY(:), TZ(:)
    real(RK), pointer :: PX(:), PY(:), PZ(:)
    real(RK), pointer :: RXTest(:), RYTest(:), RZTest(:)
    real(RK), pointer :: OXTest(:), OYTest(:), OZTest(:)
    real(RK), pointer :: PXTest(:), PYTest(:), PZTest(:)

  end type TSiteQuadrupole

  interface Construct
    module procedure TSiteQuadrupole_Construct
  end interface

  interface Destruct
    module procedure TSiteQuadrupole_Destruct
  end interface

  interface Allocate
    module procedure TSiteQuadrupole_Allocate
  end interface

  interface Deallocate
    module procedure TSiteQuadrupole_Deallocate
  end interface

  interface Save
    module procedure TSiteQuadrupole_Save
  end interface


contains


!==============================================================!
!  Subroutine TSiteLJ126_Construct                             !
!==============================================================!

  subroutine TSiteLJ126_Construct( this )

    implicit none

    ! Declare arguments
    type(TSiteLJ126) :: this
    
    ! Declare local variables
    integer          :: stat

    ! Read site parameters
    if( UseIntDegFreed ) then
        call FileReadParameter( this%SiteId, iounit_potmod, IdLJ126_SiteId, .false. )
    end if
    
    nullify ( this%r )
    allocate( this%r( 3 ), STAT = stat )
    call AllocationError( stat, 'coordinates', 3 )

    call FileReadParameter( this%r(1), iounit_potmod, IdLJ126_r1, .false. )
    call FileReadParameter( this%r(2), iounit_potmod, IdLJ126_r2, .false. )
    call FileReadParameter( this%r(3), iounit_potmod, IdLJ126_r3, .false. )
    call FileReadParameter( this%sig, iounit_potmod, IdLJ126_sig, .false. )
    call FileReadParameter( this%eps, iounit_potmod, IdLJ126_eps, .false. )
    call FileReadParameter( this%mass, iounit_potmod, IdLJ126_mass, .false. )

    ! Convert to SI units
    this%r(:) = this%r(:) * Angstroem
    this%sig = this%sig * Angstroem
    this%eps = this%eps * kBoltzmann
    this%mass = this%mass * .001_RK / NAvogadro

    ! Convert to derived units
    this%r(:) = this%r(:) / UnitLength
    this%sig = this%sig / UnitLength
    this%eps = this%eps / UnitEnergy
    this%mass = this%mass / UnitMass

  end subroutine TSiteLJ126_Construct


!==============================================================!
!  Subroutine TSiteLJ126_Destruct                              !
!==============================================================!

  subroutine TSiteLJ126_Destruct( this )

    implicit none

    ! Declare arguments
    type(TSiteLJ126) :: this

    ! Destroy site
    continue

  end subroutine TSiteLJ126_Destruct


!==============================================================!
!  Subroutine TSiteLJ126_Allocate                              !
!==============================================================!

  subroutine TSiteLJ126_Allocate( this )

    implicit none

    ! Declare arguments
    type(TSiteLJ126) :: this

    ! Declare local variables
    integer :: np, nt
    integer :: stat

    ! Assign local variables
    np = this%NPartMax
    nt = this%NTest

    ! Nullify pointers
    nullify( this%RX )
    nullify( this%RY )
    nullify( this%RZ )
    nullify( this%FX )
    nullify( this%FY )
    nullify( this%FZ )
    nullify( this%RXTest )
    nullify( this%RYTest )
    nullify( this%RZTest )

    ! Allocate arrays
    allocate( this%RX( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%RY( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%RZ( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    if( SimulationType .eq. MolecularDynamics ) then
      allocate( this%FX( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%FY( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%FZ( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
    end if

    if( nt > 0 ) then
      allocate( this%RXTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
      allocate( this%RYTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
      allocate( this%RZTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
    end if

  end subroutine TSiteLJ126_Allocate



!==============================================================!
!  Subroutine TSiteLJ126_Deallocate                            !
!==============================================================!

  subroutine TSiteLJ126_Deallocate( this )

    implicit none

    ! Declare arguments
    type(TSiteLJ126) :: this

    ! Deallocate arrays
    if( associated( this%r ) ) deallocate( this%r )
    if( associated( this%RX ) ) deallocate( this%RX )
    if( associated( this%RY ) ) deallocate( this%RY )
    if( associated( this%RZ ) ) deallocate( this%RZ )
    if( associated( this%FX ) ) deallocate( this%FX )
    if( associated( this%FY ) ) deallocate( this%FY )
    if( associated( this%FZ ) ) deallocate( this%FZ )
    if( associated( this%RXTest ) ) deallocate( this%RXTest )
    if( associated( this%RYTest ) ) deallocate( this%RYTest )
    if( associated( this%RZTest ) ) deallocate( this%RZTest )

  end subroutine TSiteLJ126_Deallocate



!==============================================================!
!  Subroutine TSiteLJ126_Save                                  !
!==============================================================!

  subroutine TSiteLJ126_Save( this )

    implicit none

    ! Declare arguments
    type(TSiteLJ126) :: this

    ! Save site parameters
    if( UseIntDegFreed ) then 
        write( IOBuffer, '(I3)' ), this%SiteId
        call FileWriteParameter( iounit_normal, IdLJ126_SiteId )
    end if
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%r(1) * UnitLength / Angstroem, this%r(1)
    call FileWriteParameter( iounit_normal, IdLJ126_r1 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%r(2) * UnitLength / Angstroem, this%r(2)
    call FileWriteParameter( iounit_normal, IdLJ126_r2 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%r(3) * UnitLength / Angstroem, this%r(3)
    call FileWriteParameter( iounit_normal, IdLJ126_r3 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%sig * UnitLength / Angstroem, this%sig
    call FileWriteParameter( iounit_normal, IdLJ126_sig )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%eps * UnitEnergy / kBoltzmann, this%eps
    call FileWriteParameter( iounit_normal, IdLJ126_eps )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%mass * UnitMass * 1000._RK * NAvogadro, this%mass
    call FileWriteParameter( iounit_normal, IdLJ126_mass )

  end subroutine TSiteLJ126_Save



!==============================================================!
!  Subroutine TSiteCharge_Construct                            !
!==============================================================!

  subroutine TSiteCharge_Construct( this )

    implicit none

    ! Declare arguments
    type(TSiteCharge) :: this
    
    ! Declare local variables
    integer          :: stat

    ! Read site parameters
    if( UseIntDegFreed ) then
      call FileReadParameter( this%SiteId, iounit_potmod, IdCharge_SiteId, .false. )
    end if
    
    nullify ( this%r )
    allocate( this%r( 3 ), STAT = stat )
    call AllocationError( stat, 'coordinates', 3 )

    call FileReadParameter( this%r(1), iounit_potmod, IdCharge_r1, .false. )
    call FileReadParameter( this%r(2), iounit_potmod, IdCharge_r2, .false. )
    call FileReadParameter( this%r(3), iounit_potmod, IdCharge_r3, .false. )
    call FileReadParameter( this%e, iounit_potmod, IdCharge_e, .false. )
    call FileReadParameter( this%mass, iounit_potmod, IdCharge_mass, .false. )
    call FileReadParameter( this%shield, iounit_potmod, IdCharge_shield, .false. )

    ! Convert to SI units
    this%r(:) = this%r(:) * Angstroem
    this%e = this%e * ElementaryCharge
    this%mass = this%mass * .001_RK / NAvogadro
    this%shield = this%shield * Angstroem

    ! Convert to derived units
    this%r(:) = this%r(:) / UnitLength
    this%e = this%e / UnitCharge
    this%mass = this%mass / UnitMass
    this%shield = this%shield / UnitLength

  end subroutine TSiteCharge_Construct



!==============================================================!
!  Subroutine TSiteCharge_Destruct                             !
!==============================================================!

  subroutine TSiteCharge_Destruct( this )

    implicit none

    ! Declare arguments
    type(TSiteCharge) :: this

    ! Destroy site
    continue

  end subroutine TSiteCharge_Destruct



!==============================================================!
!  Subroutine TSiteCharge_Allocate                             !
!==============================================================!

  subroutine TSiteCharge_Allocate( this )

    implicit none

    ! Declare arguments
    type(TSiteCharge) :: this

    ! Declare local variables
    integer :: np, nt
    integer :: stat

    ! Assign local variables
    np = this%NPartMax
    nt = this%NTest

    ! Nullify pointers
    nullify( this%RX )
    nullify( this%RY )
    nullify( this%RZ )
    nullify( this%FX )
    nullify( this%FY )
    nullify( this%FZ )
    nullify( this%RXTest )
    nullify( this%RYTest )
    nullify( this%RZTest )

    ! Allocate arrays
    allocate( this%RX( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%RY( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%RZ( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    if( SimulationType .eq. MolecularDynamics ) then
      allocate( this%FX( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%FY( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%FZ( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
    end if

    if( nt > 0 ) then
      allocate( this%RXTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
      allocate( this%RYTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
      allocate( this%RZTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
    end if

  end subroutine TSiteCharge_Allocate


!==============================================================!
!  Subroutine TSiteCharge_Deallocate                           !
!==============================================================!

  subroutine TSiteCharge_Deallocate( this )

    implicit none

    ! Declare arguments
    type(TSiteCharge) :: this

    ! Deallocate arrays
    if( associated( this%r ) ) deallocate( this%r )
    if( associated( this%RX ) ) deallocate( this%RX )
    if( associated( this%RY ) ) deallocate( this%RY )
    if( associated( this%RZ ) ) deallocate( this%RZ )
    if( associated( this%FX ) ) deallocate( this%FX )
    if( associated( this%FY ) ) deallocate( this%FY )
    if( associated( this%FZ ) ) deallocate( this%FZ )
    if( associated( this%RXTest ) ) deallocate( this%RXTest )
    if( associated( this%RYTest ) ) deallocate( this%RYTest )
    if( associated( this%RZTest ) ) deallocate( this%RZTest )

  end subroutine TSiteCharge_Deallocate



!==============================================================!
!  Subroutine TSiteCharge_Save                                 !
!==============================================================!

  subroutine TSiteCharge_Save( this )

    implicit none

    ! Declare arguments
    type(TSiteCharge) :: this

    ! Save site parameters
    if( UseIntDegFreed ) then
        write( IOBuffer, '(I3)' ), this%SiteId
        call FileWriteParameter( iounit_normal, IdCharge_SiteId ) 
    end if

    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%r(1) * UnitLength / Angstroem, this%r(1)
    call FileWriteParameter( iounit_normal, IdCharge_r1 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%r(2) * UnitLength / Angstroem, this%r(2)
    call FileWriteParameter( iounit_normal, IdCharge_r2 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%r(3) * UnitLength / Angstroem, this%r(3)
    call FileWriteParameter( iounit_normal, IdCharge_r3 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%e * UnitCharge / ElementaryCharge,  this%e
    call FileWriteParameter( iounit_normal, IdCharge_e )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%mass * UnitMass * 1000._RK * NAvogadro, this%mass
    call FileWriteParameter( iounit_normal, IdCharge_mass )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%shield * UnitLength / Angstroem, this%shield
    call FileWriteParameter( iounit_normal, IdCharge_shield )

  end subroutine TSiteCharge_Save



!==============================================================!
!  Subroutine TSiteDipole_Construct                            !
!==============================================================!

  subroutine TSiteDipole_Construct( this )

    implicit none

    ! Declare arguments
    type(TSiteDipole) :: this

    ! Declare local variables
    real(RK) :: theta, phi
    integer :: stat

    ! Read site parameters
    if( UseIntDegFreed ) then
        call FileReadParameter( this%SiteId, iounit_potmod, IdDipole_SiteId, .false. )
    end if

    nullify ( this%r )
    nullify ( this%or )
    allocate( this%r( 3 ), STAT = stat )
    call AllocationError( stat, 'coordinates particles', 3 )
    allocate( this%or( 3 ), STAT = stat )
    call AllocationError( stat, 'coordinates particles', 3 )
    call FileReadParameter( this%r(1), iounit_potmod, IdDipole_r1, .false. )
    call FileReadParameter( this%r(2), iounit_potmod, IdDipole_r2, .false. )
    call FileReadParameter( this%r(3), iounit_potmod, IdDipole_r3, .false. )
    call FileReadParameter( theta, iounit_potmod, IdDipole_theta, .false. )
    call FileReadParameter( phi, iounit_potmod, IdDipole_phi, .false. )
    call FileReadParameter( this%D, iounit_potmod, IdDipole_D, .false. )
    call FileReadParameter( this%mass, iounit_potmod, IdDipole_mass, .false. )
    call FileReadParameter( this%shield, iounit_potmod, IdDipole_shield, .false. )

    ! Convert to SI units
    this%r(:) = this%r(:) * Angstroem
    this%D = this%D / DebyesInSI
    this%mass = this%mass * .001_RK / NAvogadro
    theta = theta / DegreesInRadian
    phi = phi / DegreesInRadian
    this%shield = this%shield * Angstroem

    ! Convert to derived units
    this%r(:) = this%r(:) / UnitLength
    this%or(1) = sin( theta ) * cos( phi )
    this%or(2) = sin( theta ) * sin( phi )
    this%or(3) = cos( theta )
    this%D = this%D / UnitDipole
    this%mass = this%mass / UnitMass
    this%shield = this%shield / UnitLength

  end subroutine TSiteDipole_Construct



!==============================================================!
!  Subroutine TSiteDipole_Destruct                             !
!==============================================================!

  subroutine TSiteDipole_Destruct( this )

    implicit none

    ! Declare arguments
    type(TSiteDipole) :: this

    ! Destroy site
    continue

  end subroutine TSiteDipole_Destruct



!==============================================================!
!  Subroutine TSiteDipole_Allocate                             !
!==============================================================!

  subroutine TSiteDipole_Allocate( this )

    implicit none

    ! Declare arguments
    type(TSiteDipole) :: this

    ! Declare local variables
    integer :: np, nt
    integer :: stat

    ! Assign local variables
    np = this%NPartMax
    nt = this%NTest

    ! Nullify pointers
    nullify( this%RX )
    nullify( this%RY )
    nullify( this%RZ )
    nullify( this%OX )
    nullify( this%OY )
    nullify( this%OZ )
    nullify( this%FX )
    nullify( this%FY )
    nullify( this%FZ )
    nullify( this%TX )
    nullify( this%TY )
    nullify( this%TZ )
    nullify( this%RXTest )
    nullify( this%RYTest )
    nullify( this%RZTest )
    nullify( this%OXTest )
    nullify( this%OYTest )
    nullify( this%OZTest )

    ! Allocate arrays
    allocate( this%RX( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%RY( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%RZ( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%OX( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%OY( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%OZ( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    if( SimulationType .eq. MolecularDynamics ) then
      allocate( this%FX( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%FY( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%FZ( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%TX( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%TY( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%TZ( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
    end if

    if( nt > 0 ) then
      allocate( this%RXTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
      allocate( this%RYTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
      allocate( this%RZTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
      allocate( this%OXTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
      allocate( this%OYTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
      allocate( this%OZTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
    end if

  end subroutine TSiteDipole_Allocate



!==============================================================!
!  Subroutine TSiteDipole_Deallocate                           !
!==============================================================!

  subroutine TSiteDipole_Deallocate( this )

    implicit none

    ! Declare arguments
    type(TSiteDipole) :: this

    ! Deallocate arrays
    if( associated( this%r ) ) deallocate( this%r )
    if( associated( this%or ) ) deallocate( this%or )
    if( associated( this%RX ) ) deallocate( this%RX )
    if( associated( this%RY ) ) deallocate( this%RY )
    if( associated( this%RZ ) ) deallocate( this%RZ )
    if( associated( this%OX ) ) deallocate( this%OX )
    if( associated( this%OY ) ) deallocate( this%OY )
    if( associated( this%OZ ) ) deallocate( this%OZ )
    if( associated( this%FX ) ) deallocate( this%FX )
    if( associated( this%FY ) ) deallocate( this%FY )
    if( associated( this%FZ ) ) deallocate( this%FZ )
    if( associated( this%TX ) ) deallocate( this%TX )
    if( associated( this%TY ) ) deallocate( this%TY )
    if( associated( this%TZ ) ) deallocate( this%TZ )
    if( associated( this%RXTest ) ) deallocate( this%RXTest )
    if( associated( this%RYTest ) ) deallocate( this%RYTest )
    if( associated( this%RZTest ) ) deallocate( this%RZTest )
    if( associated( this%OXTest ) ) deallocate( this%OXTest )
    if( associated( this%OYTest ) ) deallocate( this%OYTest )
    if( associated( this%OZTest ) ) deallocate( this%OZTest )

  end subroutine TSiteDipole_Deallocate



!==============================================================!
!  Subroutine TSiteDipole_Save                                 !
!==============================================================!

  subroutine TSiteDipole_Save( this )

    implicit none

    ! Declare arguments
    type(TSiteDipole) :: this

    ! Save site parameters
    if( UseIntDegFreed ) then 
        write( IOBuffer, '(I3)' ), this%SiteId
        call FileWriteParameter( iounit_normal, IdDipole_SiteId )
    end if
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%r(1) * UnitLength / Angstroem, this%r(1)
    call FileWriteParameter( iounit_normal, IdDipole_r1 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%r(2) * UnitLength / Angstroem, this%r(2)
    call FileWriteParameter( iounit_normal, IdDipole_r2 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%r(3) * UnitLength / Angstroem, this%r(3)
    call FileWriteParameter( iounit_normal, IdDipole_r3 )
    write( IOBuffer, '(G20.10)' ) acos( this%or(3) ) * DegreesInRadian
    call FileWriteParameter( iounit_normal, IdDipole_theta )
    if( abs( this%or(1) ) > Zero .or. abs( this%or(2) ) > Zero ) then
      write( IOBuffer, '(G20.10)' ) &
&       atan2( this%or(2), this%or(1) ) * DegreesInRadian
    else
      write( IOBuffer, '(G20.10)' ) 0._RK
    end if
    call FileWriteParameter( iounit_normal, IdDipole_phi )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%D * UnitDipole * DebyesInSI, this%D
    call FileWriteParameter( iounit_normal, IdDipole_D )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%mass * UnitMass * 1000._RK * NAvogadro, this%mass
    call FileWriteParameter( iounit_normal, IdDipole_mass )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%shield * UnitLength / Angstroem, this%shield
    call FileWriteParameter( iounit_normal, IdDipole_shield )

  end subroutine TSiteDipole_Save



!==============================================================!
!  Subroutine TSiteQuadrupole_Construct                        !
!==============================================================!

  subroutine TSiteQuadrupole_Construct( this )

    implicit none

    ! Declare arguments
    type(TSiteQuadrupole) :: this

    ! Declare local variables
    real(RK) :: theta, phi
    integer  :: stat

    ! Read site parameters
    if( UseIntDegFreed ) then
        call FileReadParameter( this%SiteId, iounit_potmod, IdQuadrupole_SiteId, .false. )
    end if

    nullify ( this%r )
    nullify ( this%or )
    allocate( this%r( 3 ), STAT = stat )
    call AllocationError( stat, 'coordinates particles', 3 )
    allocate( this%or( 3 ), STAT = stat )
    call AllocationError( stat, 'coordinates particles', 3 )
    call FileReadParameter( this%r(1), iounit_potmod, IdQuadrupole_r1, .false. )
    call FileReadParameter( this%r(2), iounit_potmod, IdQuadrupole_r2, .false. )
    call FileReadParameter( this%r(3), iounit_potmod, IdQuadrupole_r3, .false. )
    call FileReadParameter( theta, iounit_potmod, IdQuadrupole_theta, .false. )
    call FileReadParameter( phi, iounit_potmod, IdQuadrupole_phi, .false. )
    call FileReadParameter( this%Q, iounit_potmod, IdQuadrupole_Q, .false. )
    call FileReadParameter( this%mass, iounit_potmod, IdQuadrupole_mass, .false. )
    call FileReadParameter( this%shield, iounit_potmod, IdQuadrupole_shield, .false. )

    ! Convert to SI units
    this%r(:) = this%r(:) * Angstroem
    this%Q = this%Q / BuckinghamsInSI
    this%mass = this%mass * .001_RK / NAvogadro
    theta = theta / DegreesInRadian
    phi = phi / DegreesInRadian
    this%shield = this%shield * Angstroem

    ! Convert to derived units
    this%r(:) = this%r(:) / UnitLength
    this%or(1) = sin( theta ) * cos( phi )
    this%or(2) = sin( theta ) * sin( phi )
    this%or(3) = cos( theta )
    this%Q = this%Q / UnitQuadrupole
    this%mass = this%mass / UnitMass
    this%shield = this%shield / UnitLength

  end subroutine TSiteQuadrupole_Construct



!==============================================================!
!  Subroutine TSiteQuadrupole_Destruct                         !
!==============================================================!

  subroutine TSiteQuadrupole_Destruct( this )

    implicit none

    ! Declare arguments
    type(TSiteQuadrupole) :: this

    ! Destroy site
    continue

  end subroutine TSiteQuadrupole_Destruct



!==============================================================!
!  Subroutine TSiteQuadrupole_Allocate                         !
!==============================================================!

  subroutine TSiteQuadrupole_Allocate( this )

    implicit none

    ! Declare arguments
    type(TSiteQuadrupole) :: this

    ! Declare local variables
    integer :: np, nt
    integer :: stat

    ! Assign local variables
    np = this%NPartMax
    nt = this%NTest

    ! Nullify pointers
    nullify( this%RX )
    nullify( this%RY )
    nullify( this%RZ )
    nullify( this%OX )
    nullify( this%OY )
    nullify( this%OZ )
    nullify( this%FX )
    nullify( this%FY )
    nullify( this%FZ )
    nullify( this%TX )
    nullify( this%TY )
    nullify( this%TZ )
    nullify( this%RXTest )
    nullify( this%RYTest )
    nullify( this%RZTest )
    nullify( this%OXTest )
    nullify( this%OYTest )
    nullify( this%OZTest )

    ! Allocate arrays
    allocate( this%RX( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%RY( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%RZ( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%OX( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%OY( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%OZ( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    if( SimulationType .eq. MolecularDynamics ) then
      allocate( this%FX( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%FY( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%FZ( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%TX( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%TY( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%TZ( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
    end if
    if( nt > 0 ) then
      allocate( this%RXTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
      allocate( this%RYTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
      allocate( this%RZTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
      allocate( this%OXTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
      allocate( this%OYTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
      allocate( this%OZTest( nt ), STAT = stat )
      call AllocationError( stat, 'test particles', nt )
    end if

  end subroutine TSiteQuadrupole_Allocate



!==============================================================!
!  Subroutine TSiteQuadrupole_Deallocate                       !
!==============================================================!

  subroutine TSiteQuadrupole_Deallocate( this )

    implicit none

    ! Declare arguments
    type(TSiteQuadrupole) :: this

    ! Deallocate arrays
    if( associated( this%r ) ) deallocate( this%r )
    if( associated( this%or ) ) deallocate( this%or )
    if( associated( this%RX ) ) deallocate( this%RX )
    if( associated( this%RY ) ) deallocate( this%RY )
    if( associated( this%RZ ) ) deallocate( this%RZ )
    if( associated( this%OX ) ) deallocate( this%OX )
    if( associated( this%OY ) ) deallocate( this%OY )
    if( associated( this%OZ ) ) deallocate( this%OZ )
    if( associated( this%FX ) ) deallocate( this%FX )
    if( associated( this%FY ) ) deallocate( this%FY )
    if( associated( this%FZ ) ) deallocate( this%FZ )
    if( associated( this%TX ) ) deallocate( this%TX )
    if( associated( this%TY ) ) deallocate( this%TY )
    if( associated( this%TZ ) ) deallocate( this%TZ )
    if( associated( this%RXTest ) ) deallocate( this%RXTest )
    if( associated( this%RYTest ) ) deallocate( this%RYTest )
    if( associated( this%RZTest ) ) deallocate( this%RZTest )
    if( associated( this%OXTest ) ) deallocate( this%OXTest )
    if( associated( this%OYTest ) ) deallocate( this%OYTest )
    if( associated( this%OZTest ) ) deallocate( this%OZTest )

  end subroutine TSiteQuadrupole_Deallocate



!==============================================================!
!  Subroutine TSiteQuadrupole_Save                             !
!==============================================================!

  subroutine TSiteQuadrupole_Save( this )

    implicit none

    ! Declare arguments
    type(TSiteQuadrupole) :: this

    ! Save site parameters
    if( UseIntDegFreed ) then 
        write( IOBuffer, '(I3)' ), this%SiteId
        call FileWriteParameter( iounit_normal, IdQuadrupole_SiteId )
    end if
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%r(1) * UnitLength / Angstroem, this%r(1)
    call FileWriteParameter( iounit_normal, IdQuadrupole_r1 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%r(2) * UnitLength / Angstroem, this%r(2)
    call FileWriteParameter( iounit_normal, IdQuadrupole_r2 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%r(3) * UnitLength / Angstroem, this%r(3)
    call FileWriteParameter( iounit_normal, IdQuadrupole_r3 )
    write( IOBuffer, '(G20.10)' ) acos( this%or(3) ) * DegreesInRadian
    call FileWriteParameter( iounit_normal, IdQuadrupole_theta )
    if( abs( this%or(1) ) > Zero .or. abs( this%or(2) ) > Zero ) then
      write( IOBuffer, '(G20.10)' ) &
&       atan2( this%or(2), this%or(1) ) * DegreesInRadian
    else
      write( IOBuffer, '(G20.10)' ) 0._RK
    end if
    call FileWriteParameter( iounit_normal, IdQuadrupole_phi )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%Q * UnitQuadrupole * BuckinghamsInSI, this%Q
    call FileWriteParameter( iounit_normal, IdQuadrupole_Q )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%mass * UnitMass * 1000._RK * NAvogadro, this%mass
    call FileWriteParameter( iounit_normal, IdQuadrupole_mass )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%shield * UnitLength / Angstroem, this%shield
    call FileWriteParameter( iounit_normal, IdQuadrupole_shield )

  end subroutine TSiteQuadrupole_Save



end module ms2_site
