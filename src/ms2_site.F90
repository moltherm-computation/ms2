!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 2.0 + IDF          !
!  (c) 2014 by TU Kaiserslautern                               !
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
    integer, pointer  :: RDFSum(:)

#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer :: vsLJx(:) , vsLJy(:), vsLJz(:)
    real(RK), pointer :: vsuLJx(:), vsuLJy(:), vsuLJz(:)
    real(RK), pointer :: vbLJx(:) , vbLJy(:), vbLJz(:)
    real(RK), pointer :: cLJx(:)  , cLJy(:),  cLJz(:)
    real(RK), pointer :: tuLJx(:),  tuLJy(:),  tuLJz(:)
    real(RK), pointer :: tlLJx(:),  tlLJy(:),  tlLJz(:)
    real(RK), pointer :: tdLJx(:),  tdLJy(:),  tdLJz(:)
    real(RK), pointer :: Qm0r(:,:)
!TRANSPORT_END
#endif

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

#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer :: vsCx(:), vsCy(:), vsCz(:)
    real(RK), pointer :: vsuCx(:), vsuCy(:), vsuCz(:)
    real(RK), pointer :: vbCx(:), vbCy(:), vbCz(:)
    real(RK), pointer :: cCx(:),  cCy(:),  cCz(:)
    real(RK), pointer :: tuCx(:),  tuCy(:),  tuCz(:)
    real(RK), pointer :: tlCx(:),  tlCy(:),  tlCz(:)
    real(RK), pointer :: tdCx(:),  tdCy(:),  tdCz(:)
    real(RK), pointer :: Qm0r(:,:)
!TRANSPORT_END
#endif

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

#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer :: vsDx(:), vsDy(:), vsDz(:)
    real(RK), pointer :: vsuDx(:), vsuDy(:), vsuDz(:)
    real(RK), pointer :: vbDx(:), vbDy(:), vbDz(:)
    real(RK), pointer :: cDx(:),  cDy(:),  cDz(:)
    real(RK), pointer :: tuDx(:),  tuDy(:),  tuDz(:)
    real(RK), pointer :: tlDx(:),  tlDy(:),  tlDz(:)
    real(RK), pointer :: tdDx(:),  tdDy(:),  tdDz(:)
    real(RK), pointer :: Qm0r(:,:)
!TRANSPORT_END
#endif

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

#if  TRANS == 1
    !TRANSPORT_start
    real(RK), pointer :: vsQx(:) , vsQy(:), vsQz(:)
    real(RK), pointer :: vsuQx(:), vsuQy(:), vsuQz(:)
    real(RK), pointer :: vbQx(:), vbQy(:), vbQz(:)
    real(RK), pointer :: cQx(:) ,  cQy(:),  cQz(:)
    real(RK), pointer :: tuQx(:),  tuQy(:),  tuQz(:)
    real(RK), pointer :: tlQx(:),  tlQy(:),  tlQz(:)
    real(RK), pointer :: tdQx(:),  tdQy(:),  tdQz(:)
    real(RK), pointer :: Qm0r(:,:)
!TRANSPORT_END
#endif

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
    nullify( this%RDFSum)

#if  TRANS == 1
    !TRANSPORT_start
    nullify( this%vsLJx )
    nullify( this%vsLJy )
    nullify( this%vsLJz )
    nullify( this%vsuLJx )
    nullify( this%vsuLJy )
    nullify( this%vsuLJz )
    nullify( this%vbLJx )
    nullify( this%vbLJy )
    nullify( this%vbLJz )
    nullify( this%cLJx )
    nullify( this%cLJy )
    nullify( this%cLJz )
    nullify( this%tuLJx )
    nullify( this%tuLJy )
    nullify( this%tuLJz )
    nullify( this%tlLJx )
    nullify( this%tlLJy )
    nullify( this%tlLJz )
    nullify( this%tdLJx )
    nullify( this%tdLJy )
    nullify( this%tdLJz )
!TRANSPORT_END
#endif

    ! Allocate arrays
    allocate( this%RX( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%RY( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    allocate( this%RZ( np ), STAT = stat )
    call AllocationError( stat, 'particles', np )
    if( RDFUpdateFrequency > 0 ) then
      allocate( this%RDFSum(RDFNumberShells+10), STAT = stat )
      call AllocationError( stat, 'RDFSum', RDFNumberShells+10 )
    endif     
    
    call AllocationError( stat, 'particles', np )    
    if( SimulationType .eq. MolecularDynamics ) then
      allocate( this%FX( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%FY( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%FZ( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )

#if  TRANS == 1
!TRANSPORT_start
      allocate( this%vsLJx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsLJy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsLJz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsuLJx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsuLJy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsuLJz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vbLJx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vbLJy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vbLJz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%cLJx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%cLJy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%cLJz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tuLJx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tuLJy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tuLJz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tlLJx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tlLJy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tlLJz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tdLJx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tdLJy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tdLJz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
!TRANSPORT_END
#endif
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
    if( associated( this%RX ) ) then
      deallocate( this%RX )
    end if
    if( associated( this%RY ) ) then
      deallocate( this%RY )
    end if
    if( associated( this%RZ ) ) then
      deallocate( this%RZ )
    end if
    if( associated( this%RDFSum ) ) then
      deallocate( this%RDFSum )
    end if    
    if( associated( this%FX ) ) then
      deallocate( this%FX )
    end if
    if( associated( this%FY ) ) then
      deallocate( this%FY )
    end if
    if( associated( this%FZ ) ) then
      deallocate( this%FZ )
    end if
    if( associated( this%RXTest ) ) then
      deallocate( this%RXTest )
    end if
    if( associated( this%RYTest ) ) then
      deallocate( this%RYTest )
    end if
    if( associated( this%RZTest ) ) then
      deallocate( this%RZTest )
    end if

#if  TRANS == 1
    !TRANSPORT_start
    if( associated( this%vsLJx ) ) then
      deallocate( this%vsLJx )
    end if
    if( associated( this%vsLJy ) ) then
      deallocate( this%vsLJy )
    end if
    if( associated( this%vsLJz ) ) then
      deallocate( this%vsLJz )
    end if
    if( associated( this%vsuLJx ) ) then
      deallocate( this%vsuLJx )
    end if
    if( associated( this%vsuLJy ) ) then
      deallocate( this%vsuLJy )
    end if
    if( associated( this%vsuLJz ) ) then
      deallocate( this%vsuLJz )
    end if
    if( associated( this%vbLJx ) ) then
     deallocate( this%vbLJx )
    end if
    if( associated( this%vbLJy ) ) then
      deallocate( this%vbLJy )
    end if
    if( associated( this%vbLJz ) ) then
      deallocate( this%vbLJz )
    end if
    if( associated( this%cLJx ) ) then
      deallocate( this%cLJx )
    end if
    if( associated( this%cLJy ) ) then
      deallocate( this%cLJy )
    end if
    if( associated( this%cLJz ) ) then
      deallocate( this%cLJz )
    end if
    if( associated( this%tuLJx ) ) then
      deallocate( this%tuLJx )
    end if
    if( associated( this%tuLJy ) ) then
      deallocate( this%tuLJy )
    end if
    if( associated( this%tuLJz ) ) then
      deallocate( this%tuLJz )
    end if
    if( associated( this%tlLJx ) ) then
      deallocate( this%tlLJx )
    end if
    if( associated( this%tlLJy ) ) then
      deallocate( this%tlLJy )
    end if
    if( associated( this%tlLJz ) ) then
      deallocate( this%tlLJz )
    end if
    if( associated( this%tdLJx ) ) then
      deallocate( this%tdLJx )
    end if
    if( associated( this%tdLJy ) ) then
      deallocate( this%tdLJy )
    end if
    if( associated( this%tdLJz ) ) then
      deallocate( this%tdLJz )
    end if
!TRANSPORT_END
#endif
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
        write( IOBuffer, '(I3)' ) this%SiteId
        call FileWriteParameter( iounit_normal, IdLJ126_SiteId )
    end if
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%r(1) * UnitLength / Angstroem, this%r(1)
    call FileWriteParameter( iounit_normal, IdLJ126_r1 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%r(2) * UnitLength / Angstroem, this%r(2)
    call FileWriteParameter( iounit_normal, IdLJ126_r2 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%r(3) * UnitLength / Angstroem, this%r(3)
    call FileWriteParameter( iounit_normal, IdLJ126_r3 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%sig * UnitLength / Angstroem, this%sig
    call FileWriteParameter( iounit_normal, IdLJ126_sig )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%eps * UnitEnergy / kBoltzmann, this%eps
    call FileWriteParameter( iounit_normal, IdLJ126_eps )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%mass * UnitMass * 1000._RK * NAvogadro, this%mass
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

#if  TRANS == 1
    !TRANSPORT_start
    nullify( this%vsCx )
    nullify( this%vsCy )
    nullify( this%vsCz )
    nullify( this%vsuCx )
    nullify( this%vsuCy )
    nullify( this%vsuCz )
    nullify( this%vbCx )
    nullify( this%vbCy )
    nullify( this%vbCz )
    nullify( this%cCx )
    nullify( this%cCy )
    nullify( this%cCz )
    nullify( this%tuCx )
    nullify( this%tuCy )
    nullify( this%tuCz )
    nullify( this%tlCx )
    nullify( this%tlCy )
    nullify( this%tlCz )
    nullify( this%tdCx )
    nullify( this%tdCy )
    nullify( this%tdCz )
!TRANSPORT_END
#endif

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

#if  TRANS == 1
      !TRANSPORT_start
      allocate( this%vsCx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsCy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsCz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsuCx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsuCy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsuCz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vbCx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vbCy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vbCz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%cCx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%cCy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%cCz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tuCx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tuCy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tuCz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tlCx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tlCy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tlCz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tdCx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tdCy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tdCz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      !TRANSPORT_END
#endif
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
    if( associated( this%RX ) ) then
      deallocate( this%RX )
    end if
    if( associated( this%RY ) ) then
      deallocate( this%RY )
    end if
    if( associated( this%RZ ) ) then
      deallocate( this%RZ )
    end if
    if( associated( this%FX ) ) then
      deallocate( this%FX )
    end if
    if( associated( this%FY ) ) then
      deallocate( this%FY )
    end if
    if( associated( this%FZ ) ) then
      deallocate( this%FZ )
    end if
    if( associated( this%RXTest ) ) then
      deallocate( this%RXTest )
    end if
    if( associated( this%RYTest ) ) then
      deallocate( this%RYTest )
    end if
    if( associated( this%RZTest ) ) then
      deallocate( this%RZTest )
    end if

#if  TRANS == 1
    !TRANSPORT_start
    if( associated( this%vsCx ) )  then
      deallocate( this%vsCx )
    end if
    if( associated( this%vsCy ) ) then
      deallocate( this%vsCy )
    end if
    if( associated( this%vsCz ) ) then
      deallocate( this%vsCz )
    end if
    if( associated( this%vbCx ) ) then
      deallocate( this%vbCx )
    end if
    if( associated( this%vbCy ) ) then
      deallocate( this%vbCy )
    end if
    if( associated( this%vbCz ) ) then
      deallocate( this%vbCz )
    end if
    if( associated( this%cCx ) ) then
      deallocate( this%cCx )
    end if
    if( associated( this%cCy ) ) then
      deallocate( this%cCy )
    end if
    if( associated( this%cCz ) ) then
      deallocate( this%cCz )
    end if
    if( associated( this%tuCx ) ) then
      deallocate( this%tuCx )
    end if
    if( associated( this%tuCy ) ) then
      deallocate( this%tuCy )
    end if
    if( associated( this%tuCz ) ) then
     deallocate( this%tuCz )
    end if
    if( associated( this%tlCx ) ) then
      deallocate( this%tlCx )
    end if
    if( associated( this%tlCy ) ) then
      deallocate( this%tlCy )
    end if
    if( associated( this%tlCz ) ) then
      deallocate( this%tlCz )
    end if
    if( associated( this%tdCx ) ) then
      deallocate( this%tdCx )
    end if
    if( associated( this%tdCy ) ) then
      deallocate( this%tdCy )
    end if
    if( associated( this%tdCz ) ) then
      deallocate( this%tdCz )
    end if
!TRANSPORT_END
#endif

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
        write( IOBuffer, '(I3)' ) this%SiteId
        call FileWriteParameter( iounit_normal, IdCharge_SiteId ) 
    end if
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%r(1) * UnitLength / Angstroem, this%r(1)
    call FileWriteParameter( iounit_normal, IdCharge_r1 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%r(2) * UnitLength / Angstroem, this%r(2)
    call FileWriteParameter( iounit_normal, IdCharge_r2 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%r(3) * UnitLength / Angstroem, this%r(3)
    call FileWriteParameter( iounit_normal, IdCharge_r3 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%e * UnitCharge / ElementaryCharge,  this%e
    call FileWriteParameter( iounit_normal, IdCharge_e )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%mass * UnitMass * 1000._RK * NAvogadro, this%mass
    call FileWriteParameter( iounit_normal, IdCharge_mass )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%shield * UnitLength / Angstroem, this%shield
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
    integer  :: stat

    ! Read site parameters
    if( UseIntDegFreed ) then
        call FileReadParameter( this%SiteId, iounit_potmod, IdDipole_SiteId, .false. )
    end if

    nullify ( this%r )
    allocate( this%r( 3 ), STAT = stat )
    call AllocationError( stat, 'coordinates', 3 )
    nullify ( this%or )
    allocate( this%or( 3 ), STAT = stat )
    call AllocationError( stat, 'coordinates', 3 )

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

#if  TRANS == 1
    !TRANSPORT_start
    nullify( this%vsDx )
    nullify( this%vsDy )
    nullify( this%vsDz )
    nullify( this%vsuDx )
    nullify( this%vsuDy )
    nullify( this%vsuDz )
    nullify( this%vbDx )
    nullify( this%vbDy )
    nullify( this%vbDz )
    nullify( this%cDx )
    nullify( this%cDy )
    nullify( this%cDz )
    nullify( this%tuDx )
    nullify( this%tuDy )
    nullify( this%tuDz )
    nullify( this%tlDx )
    nullify( this%tlDy )
    nullify( this%tlDz )
    nullify( this%tdDx )
    nullify( this%tdDy )
    nullify( this%tdDz )
!TRANSPORT_END
#endif

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

#if  TRANS == 1
      !TRANSPORT_start
      allocate( this%vsDx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsDy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsDz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsuDx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsuDy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsuDz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vbDx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vbDy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vbDz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%cDx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%cDy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%cDz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tuDx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tuDy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tuDz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tlDx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tlDy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tlDz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tdDx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tdDy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tdDz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
!TRANSPORT_END
#endif
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
    if( associated( this%RX ) ) then
      deallocate( this%RX )
    end if
    if( associated( this%RY ) ) then
      deallocate( this%RY )
    end if
    if( associated( this%RZ ) ) then
      deallocate( this%RZ )
    end if
    if( associated( this%OX ) ) then
      deallocate( this%OX )
    end if
    if( associated( this%OY ) ) then
      deallocate( this%OY )
    end if
    if( associated( this%OZ ) ) then
      deallocate( this%OZ )
    end if
    if( associated( this%FX ) ) then
      deallocate( this%FX )
    end if
    if( associated( this%FY ) ) then
      deallocate( this%FY )
    end if
    if( associated( this%FZ ) ) then
      deallocate( this%FZ )
    end if
    if( associated( this%TX ) ) then
      deallocate( this%TX )
    end if
    if( associated( this%TY ) ) then
      deallocate( this%TY )
    end if
    if( associated( this%TZ ) ) then
      deallocate( this%TZ )
    end if
    if( associated( this%RXTest ) ) then
      deallocate( this%RXTest )
    end if
    if( associated( this%RYTest ) ) then
      deallocate( this%RYTest )
    end if
    if( associated( this%RZTest ) ) then
      deallocate( this%RZTest )
    end if
    if( associated( this%OXTest ) ) then
      deallocate( this%OXTest )
    end if
    if( associated( this%OYTest ) ) then
      deallocate( this%OYTest )
    end if
    if( associated( this%OZTest ) ) then
      deallocate( this%OZTest )
    end if

#if  TRANS == 1
    !TRANSPORT_start
    if( associated( this%vsDx ) ) then
      deallocate( this%vsDx )
    end if
    if( associated( this%vsDy ) ) then
      deallocate( this%vsDy )
    end if
    if( associated( this%vsDz ) ) then
      deallocate( this%vsDz )
    end if
    if( associated( this%vsuDx ) ) then
      deallocate( this%vsuDx )
    end if
    if( associated( this%vsuDy ) ) then
      deallocate( this%vsuDy )
    end if
    if( associated( this%vsuDz ) ) then
      deallocate( this%vsuDz )
    end if
    if( associated( this%vbDx ) ) then
      deallocate( this%vbDx )
    end if
    if( associated( this%vbDy ) ) then
      deallocate( this%vbDy )
    end if
    if( associated( this%vbDz ) ) then
      deallocate( this%vbDz )
    end if
    if( associated( this%cDx ) ) then
      deallocate( this%cDx )
    end if
    if( associated( this%cDy ) ) then
      deallocate( this%cDy )
    end if
    if( associated( this%cDz ) ) then
      deallocate( this%cDz )
    end if
    if( associated( this%tuDx ) ) then
      deallocate( this%tuDx )
    end if
    if( associated( this%tuDy ) ) then
      deallocate( this%tuDy )
    end if
    if( associated( this%tuDz ) ) then
      deallocate( this%tuDz )
    end if
    if( associated( this%tlDx ) ) then
      deallocate( this%tlDx )
    end if
    if( associated( this%tlDy ) ) then
      deallocate( this%tlDy )
    end if
    if( associated( this%tlDz ) ) then
      deallocate( this%tlDz )
    end if
    if( associated( this%tdDx ) ) then
      deallocate( this%tdDx )
    end if
    if( associated( this%tdDy ) ) then
      deallocate( this%tdDy )
    end if
    if( associated( this%tdDz ) ) then
      deallocate( this%tdDz )
    end if
!TRANSPORT_END
#endif

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
        write( IOBuffer, '(I3)' ) this%SiteId
        call FileWriteParameter( iounit_normal, IdDipole_SiteId )
    end if
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%r(1) * UnitLength / Angstroem, this%r(1)
    call FileWriteParameter( iounit_normal, IdDipole_r1 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%r(2) * UnitLength / Angstroem, this%r(2)
    call FileWriteParameter( iounit_normal, IdDipole_r2 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%r(3) * UnitLength / Angstroem, this%r(3)
    call FileWriteParameter( iounit_normal, IdDipole_r3 )
    write( IOBuffer, '(G20.10)' ) acos( this%or(3) ) * DegreesInRadian
    call FileWriteParameter( iounit_normal, IdDipole_theta )
    if( abs( this%or(1) ) > Zero .or. abs( this%or(2) ) > Zero ) then
      write( IOBuffer, '(G20.10)' ) atan2( this%or(2), this%or(1) ) * DegreesInRadian
    else
      write( IOBuffer, '(G20.10)' ) 0._RK
    end if
    call FileWriteParameter( iounit_normal, IdDipole_phi )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%D * UnitDipole * DebyesInSI, this%D
    call FileWriteParameter( iounit_normal, IdDipole_D )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%mass * UnitMass * 1000._RK * NAvogadro, this%mass
    call FileWriteParameter( iounit_normal, IdDipole_mass )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%shield * UnitLength / Angstroem, this%shield
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
    allocate( this%r( 3 ), STAT = stat )
    call AllocationError( stat, 'coordinates', 3 )
    nullify ( this%or )
    allocate( this%or( 3 ), STAT = stat )
    call AllocationError( stat, 'coordinates', 3 )

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

#if  TRANS == 1
    !TRANSPORT_start
    nullify( this%vsQx )
    nullify( this%vsQy )
    nullify( this%vsQz )
    nullify( this%vsuQx )
    nullify( this%vsuQy )
    nullify( this%vsuQz )
    nullify( this%vsQx )
    nullify( this%vsQy )
    nullify( this%vsQz )
    nullify( this%vbQx )
    nullify( this%vbQy )
    nullify( this%vbQz )
    nullify( this%cQx )
    nullify( this%cQy )
    nullify( this%cQz )
    nullify( this%tuQx )
    nullify( this%tuQy )
    nullify( this%tuQz )
    nullify( this%tlQx )
    nullify( this%tlQy )
    nullify( this%tlQz )
    nullify( this%tdQx )
    nullify( this%tdQy )
    nullify( this%tdQz )
!TRANSPORT_END
#endif

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

#if  TRANS == 1
      !TRANSPORT_start
      allocate( this%vsQx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsQy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsQz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsuQx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsuQy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vsuQz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vbQx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vbQy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%vbQz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%cQx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%cQy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%cQz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tuQx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tuQy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tuQz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tlQx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tlQy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tlQz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tdQx( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tdQy( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      allocate( this%tdQz( np ), STAT = stat )
      call AllocationError( stat, 'particles', np )
      !TRANSPORT_END
#endif
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
    if( associated( this%RX ) ) then
      deallocate( this%RX )
    end if
    if( associated( this%RY ) ) then
      deallocate( this%RY )
    end if
    if( associated( this%RZ ) ) then
      deallocate( this%RZ )
    end if
    if( associated( this%OX ) ) then
      deallocate( this%OX )
    end if
    if( associated( this%OY ) ) then
      deallocate( this%OY )
    end if
    if( associated( this%OZ ) ) then
      deallocate( this%OZ )
    end if
    if( associated( this%FX ) ) then
      deallocate( this%FX )
    end if
    if( associated( this%FY ) ) then
      deallocate( this%FY )
    end if
    if( associated( this%FZ ) ) then
      deallocate( this%FZ )
    end if
    if( associated( this%TX ) ) then
      deallocate( this%TX )
    end if
    if( associated( this%TY ) ) then
      deallocate( this%TY )
    end if
    if( associated( this%TZ ) ) then
      deallocate( this%TZ )
    end if
    if( associated( this%RXTest ) ) then
      deallocate( this%RXTest )
    end if
    if( associated( this%RYTest ) ) then
      deallocate( this%RYTest )
    end if
    if( associated( this%RZTest ) ) then
      deallocate( this%RZTest )
    end if
    if( associated( this%OXTest ) ) then
      deallocate( this%OXTest )
    end if
    if( associated( this%OYTest ) ) then
      deallocate( this%OYTest )
    end if
    if( associated( this%OZTest ) ) then
      deallocate( this%OZTest )
    end if

#if  TRANS == 1
    !TRANSPORT_start
    if( associated( this%vsQx ) ) then
      deallocate( this%vsQx )
    end if
    if( associated( this%vsQy ) ) then
      deallocate( this%vsQy )
    end if
    if( associated( this%vsQz ) ) then
      deallocate( this%vsQz )
    end if
    if( associated( this%vsuQx ) ) then
      deallocate( this%vsuQx )
    end if
    if( associated( this%vsuQy ) ) then
      deallocate( this%vsuQy )
    end if
    if( associated( this%vsuQz ) ) then
      deallocate( this%vsuQz )
    end if
    if( associated( this%vbQx ) ) then
      deallocate( this%vbQx )
    end if
    if( associated( this%vbQy ) ) then
      deallocate( this%vbQy )
    end if
    if( associated( this%vbQz ) ) then
      deallocate( this%vbQz )
    end if
    if( associated( this%cQx ) ) then
      deallocate( this%cQx )
    end if
    if( associated( this%cQy ) ) then
      deallocate( this%cQy )
    end if
    if( associated( this%cQz ) ) then
      deallocate( this%cQz )
    end if
    if( associated( this%tuQx ) ) then
      deallocate( this%tuQx )
    end if
    if( associated( this%tuQy ) ) then
      deallocate( this%tuQy )
    end if
    if( associated( this%tuQz ) ) then
      deallocate( this%tuQz )
    end if
    if( associated( this%tlQx ) ) then
      deallocate( this%tlQx )
    end if
    if( associated( this%tlQy ) ) then
      deallocate( this%tlQy )
    end if
    if( associated( this%tlQz ) ) then
      deallocate( this%tlQz )
    end if
    if( associated( this%tdQx ) ) then
      deallocate( this%tdQx )
    end if
    if( associated( this%tdQy ) ) then
      deallocate( this%tdQy )
    end if
    if( associated( this%tdQz ) ) then
      deallocate( this%tdQz )
    end if
!TRANSPORT_END
#endif


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
        write( IOBuffer, '(I3)' ) this%SiteId
        call FileWriteParameter( iounit_normal, IdQuadrupole_SiteId )
    end if
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%r(1) * UnitLength / Angstroem, this%r(1)
    call FileWriteParameter( iounit_normal, IdQuadrupole_r1 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%r(2) * UnitLength / Angstroem, this%r(2)
    call FileWriteParameter( iounit_normal, IdQuadrupole_r2 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%r(3) * UnitLength / Angstroem, this%r(3)
    call FileWriteParameter( iounit_normal, IdQuadrupole_r3 )
    write( IOBuffer, '(G20.10)' ) acos( this%or(3) ) * DegreesInRadian
    call FileWriteParameter( iounit_normal, IdQuadrupole_theta )
    if( abs( this%or(1) ) > Zero .or. abs( this%or(2) ) > Zero ) then
      write( IOBuffer, '(G20.10)' ) atan2( this%or(2), this%or(1) ) * DegreesInRadian
    else
      write( IOBuffer, '(G20.10)' ) 0._RK
    end if
    call FileWriteParameter( iounit_normal, IdQuadrupole_phi )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%Q * UnitQuadrupole * BuckinghamsInSI, this%Q
    call FileWriteParameter( iounit_normal, IdQuadrupole_Q )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%mass * UnitMass * 1000._RK * NAvogadro, this%mass
    call FileWriteParameter( iounit_normal, IdQuadrupole_mass )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) this%shield * UnitLength / Angstroem, this%shield
    call FileWriteParameter( iounit_normal, IdQuadrupole_shield )

  end subroutine TSiteQuadrupole_Save



end module ms2_site
