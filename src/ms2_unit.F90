!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 2.0 + IDF          !
!  (c) 2014 by TU Kaiserslautern                               !
!      P.O. Box 67653                                          !
!      67653 Kaiserslautern                                    !
!==============================================================!
!  Program ms2                                                 !
!  This file contains the main routine                         !
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

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2_unit.F90...'
#endif

module ms2_unit

  use ms2_global
  use ms2_site

!==============================================================!
!  Type TUnit                                                  !
!==============================================================!

  type TUnit

    ! NSites in a unit
    integer  :: NSites

    ! SiteIds in a unit
    integer, allocatable :: SiteIds(:) ! Michael Sch.: make allocatables contiguous?

    ! Geometry of unit
    logical :: isElongated, is3D

    ! Type of unit
    logical :: isConstraint

    ! Number of degrees of freedom
    integer :: NDFRot, NDF

    ! Total mass of a unit
    real(RK) :: Mass

    ! Initial position of COM in molecular-fixed system
    real(RK) :: P0(3)

    ! position of COM in space-fixed system
    real(RK), pointer, contiguous :: PX(:), PY(:), PZ(:)

    ! Principal moments of inertia
    real(RK) :: MOI(3)

    ! matrix of rotation - initial orientation of units in molecular-fixed system
    real(RK) :: Q0(4)

    ! 12-6 Lennard-Jones sites
    integer :: NMIEnm
    type(TSiteMIEnm), pointer, contiguous :: SiteMIEnm(:)

    ! Coulomb sites
    integer :: NCharge
    type(TSiteCharge), pointer, contiguous :: SiteCharge(:)

    ! Dipole sites
    integer :: NDipole
    type(TSiteDipole), pointer, contiguous :: SiteDipole(:)

    ! Quadrupole sites
    integer :: NQuadrupole
    type(TSiteQuadrupole), pointer, contiguous :: SiteQuadrupole(:)

   ! Body fixed dipole vector for reaction field
    real(RK) :: Mue(3), MueSquared

    integer, pointer :: NPartMax, NPart ! Michael Sch. pointer needed here?
    integer, pointer :: NPart0, NPart1, NPart2

  end type TUnit

  interface Construct
    module procedure TUnit_Construct
  end interface

  interface Destruct
    module procedure TUnit_Destruct
  end interface

  interface Save
    module procedure TUnit_Save
  end interface

  interface FindCOM
    module procedure TUnit_FindCOM
  end interface


  type Constraint

    integer :: nSites
    integer, allocatable :: siteIDs(:)

  end type Constraint

contains



!==============================================================!
!  Subroutine TUnit_Construct                                  !
!==============================================================!

  subroutine TUnit_Construct( this, constraint, NSites )

    implicit none

    ! Declare arguments
    type(TUnit)            :: this
    logical, intent(in)    :: constraint
    integer, intent(in)    :: NSites


    ! Declare local variables
    character(16) :: stype
    integer       :: stat

    ! Set type of Unit and number of Sites
    this%isConstraint = constraint
    this%NSites = NSites

    ! Nullify pointers.

    nullify( this%SiteMIEnm )
    nullify( this%SiteCharge )
    nullify( this%SiteDipole )
    nullify( this%SiteQuadrupole )

    ! Zero number of sites

    this%NMIEnm = 0
    this%NCharge = 0
    this%NDipole = 0
    this%NQuadrupole = 0

    allocate (this%SiteIds(this%NSites), STAT = stat)
    call AllocationError( stat, 'constraints', this%NSites )
    allocate (this%SiteMIEnm(this%NSites), STAT = stat)
    call AllocationError( stat, 'unitLJ', this%NSites )
    allocate (this%SiteCharge(this%NSites), STAT = stat)
    call AllocationError( stat, 'unitCharge', this%NSites )
    allocate (this%SiteDipole(this%NSites), STAT = stat)
    call AllocationError( stat, 'unitDipole', this%NSites )
    allocate (this%SiteQuadrupole(this%NSites), STAT = stat)
    call AllocationError( stat, 'unitQuadrupole', this%NSites )


    if ( UseIntDegFreed ) then
      if ( this%isConstraint ) then
          call FileReadParameter( this%NSites, potmodFile%iounit, IdConstraint_NSites, .false. )
          call FileReadParameter_IOBuffer( potmodFile%iounit, IdConstraint_SiteIds )
          read( IOBuffer, * ) this%SiteIds

          ! Read number of rotation axes
          call FileReadParameter( stype, potmodFile%iounit, IdConstraint_NDFRot, .false. )
          select case( stype )
            case( '0' )
               this%NDFRot = 0
            case( '2' )
               this%NDFRot = 2
            case( '3' )
               this%NDFRot = 3
            case( 'AUTO', 'Auto', 'auto' )
               this%NDFRot = -1
            case default
            call Error( IdConstraint_NDFRot//' cannot be equal to '//trim( stype ) )
          end select
      else
         this%NDFRot = -1
      end if
    end if

  end subroutine TUnit_Construct


!==============================================================!
!  Subroutine TUnit_Destruct                               !
!==============================================================!

  subroutine TUnit_Destruct( this )

    implicit none

    ! Declare arguments
    type(TUnit) :: this

    ! Declare local variables
    integer :: i

    ! Deallocate arrays
    if( associated( this%SiteMIEnm ) ) then
      do i = 1, this%NMIEnm
        call Destruct( this%SiteMIEnm(i) )
      end do
      deallocate( this%SiteMIEnm )
    end if
    if( associated( this%SiteCharge ) ) then
      do i = 1, this%NCharge
        call Destruct( this%SiteCharge(i) )
      end do
      deallocate( this%SiteCharge )
    end if
    if( associated( this%SiteDipole ) ) then
      do i = 1, this%NDipole
        call Destruct( this%SiteDipole(i) )
      end do
      deallocate( this%SiteDipole )
    end if
    if( associated( this%SiteQuadrupole ) ) then
      do i = 1, this%NDipole
        call Destruct( this%SiteQuadrupole(i) )
      end do
      deallocate( this%SiteQuadrupole )
    end if

  end subroutine TUnit_Destruct

!==============================================================!
!      Subroutine TUnit_Save                                   !
!==============================================================!

  subroutine TUnit_Save( this )

    implicit none

    ! Declare arguments
    type(TUnit) :: this

    ! Declare local variables
    integer     :: i, n

    n = this%NSites

   ! Save Constraint Unit parameters
    write( IOBuffer, '(I3)' ) this%NSites
    call FileWriteParameter( normalFile%iounit, IdConstraint_NSites )
    write( IOBuffer, '(20I3)' ) (this%SiteIds(i),i=1,n)
    call FileWriteParameter( normalFile%iounit, IdConstraint_SiteIds )

    ! Save number of rotation axes
    write( IOBuffer, '(I2)' ) this%NDFRot
    call FileWriteParameter( normalFile%iounit, IdConstraint_NDFRot )

    ! Save total mass of the constraint unit
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%Mass * UnitMass * 1000._RK * NAvogadro, this%Mass
    call FileWriteParameter( normalFile%iounit, IdConstraint_Mass )

    ! Save moments of inertia
    if( this%NDFRot > 0 ) then
      write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&       this%MOI(1) * UnitInertia * 1000._RK * NAvogadro / Angstroem**2, &
&       this%MOI(1)
      call FileWriteParameter( normalFile%iounit, IdConstraint_MOI1 )
      write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&       this%MOI(2) * UnitInertia * 1000._RK * NAvogadro / Angstroem**2, &
&       this%MOI(2)
      call FileWriteParameter( normalFile%iounit, IdConstraint_MOI2 )
      write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&       this%MOI(3) * UnitInertia * 1000._RK * NAvogadro / Angstroem**2, &
&       this%MOI(3)
      call FileWriteParameter( normalFile%iounit, IdConstraint_MOI3 )
    end if

  end subroutine TUnit_Save

!==============================================================!
!      Subroutine TUnit_FindCOM                                !
!==============================================================!

  subroutine TUnit_FindCOM( this )

    implicit none

    ! Declare arguments
    type(TUnit) :: this

    ! Declare local variables
    integer  :: i, j
    real(RK) :: r(3)

    ! Calculate mass of Unit and COM position
    this%Mass = 0._RK
    r(:) = 0._RK
    do i = 1, this%NMIEnm
      this%Mass = this%Mass + this%SiteMIEnm(i)%mass
      r(:) = r(:) + this%SiteMIEnm(i)%mass * this%SiteMIEnm(i)%r(:)
    end do
    do i = 1, this%NCharge
      this%Mass = this%Mass + this%SiteCharge(i)%mass
      r(:) = r(:) + this%SiteCharge(i)%mass * this%SiteCharge(i)%r(:)
    end do
    do i = 1, this%NDipole
      this%Mass = this%Mass + this%SiteDipole(i)%mass
      r(:) = r(:) + this%SiteDipole(i)%mass * this%SiteDipole(i)%r(:)
    end do
    do i = 1, this%NQuadrupole
      this%Mass = this%Mass + this%SiteQuadrupole(i)%mass
      r(:) = r(:) + this%SiteQuadrupole(i)%mass * this%SiteQuadrupole(i)%r(:)
    end do
    r(:) = r(:) / this%Mass
    this%P0(:) = r(:) ! position of COM of unit in molecule-fixed system
    ! Move COM to zero
    do i = 1, this%NMIEnm
      do j = 1, 3
        this%SiteMIEnm(i)%r(j) = this%SiteMIEnm(i)%r(j) - r(j)
      end do
    end do
    do i = 1, this%NCharge
      do j = 1, 3
        this%SiteCharge(i)%r(j) = this%SiteCharge(i)%r(j) - r(j)
      end do
    end do
    do i = 1, this%NDipole
      do j = 1, 3
        this%SiteDipole(i)%r(j) = this%SiteDipole(i)%r(j) - r(j)
      end do
    end do
    do i = 1, this%NQuadrupole
      do j = 1, 3
        this%SiteQuadrupole(i)%r(j) = this%SiteQuadrupole(i)%r(j) - r(j)
      end do
    end do

  end subroutine TUnit_FindCOM


end module ms2_unit
