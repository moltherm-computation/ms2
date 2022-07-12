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
    integer :: NLJ126
    type(TSiteLJ126), pointer, contiguous :: SiteLJ126(:)

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

  interface FindMOI
    module procedure TUnit_FindMOI
  end interface

  interface ReadMOI
    module procedure TUnit_ReadMOI
  end interface

  interface FindNDF
    module procedure TUnit_FindNDF
  end interface

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

    nullify( this%SiteLJ126 )
    nullify( this%SiteCharge )
    nullify( this%SiteDipole )
    nullify( this%SiteQuadrupole )

    ! Zero number of sites

    this%NLJ126 = 0
    this%NCharge = 0
    this%NDipole = 0
    this%NQuadrupole = 0

    allocate (this%SiteIds(this%NSites), STAT = stat)
    call AllocationError( stat, 'constraints', this%NSites )
    allocate (this%SiteLJ126(this%NSites), STAT = stat)
    call AllocationError( stat, 'unitLJ', this%NSites )
    allocate (this%SiteCharge(this%NSites), STAT = stat)
    call AllocationError( stat, 'unitCharge', this%NSites )
    allocate (this%SiteDipole(this%NSites), STAT = stat)
    call AllocationError( stat, 'unitDipole', this%NSites )
    allocate (this%SiteQuadrupole(this%NSites), STAT = stat)
    call AllocationError( stat, 'unitQuadrupole', this%NSites )


    if ( UseIntDegFreed ) then
      if ( this%isConstraint ) then
          call FileReadParameter( this%NSites, iounit_potmod, IdConstraint_NSites, .false. )
          call FileReadParameter_IOBuffer( iounit_potmod, IdConstraint_SiteIds )
          read( IOBuffer, * ) this%SiteIds

          ! Read number of rotation axes
          call FileReadParameter( stype, iounit_potmod, IdConstraint_NDFRot, .false. )
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
    if( associated( this%SiteLJ126 ) ) then
      do i = 1, this%NLJ126
        call Destruct( this%SiteLJ126(i) )
      end do
      deallocate( this%SiteLJ126 )
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
    call FileWriteParameter( iounit_normal, IdConstraint_NSites )
    write( IOBuffer, '(20I3)' ) (this%SiteIds(i),i=1,n)
    call FileWriteParameter( iounit_normal, IdConstraint_SiteIds )

    ! Save number of rotation axes
    write( IOBuffer, '(I2)' ) this%NDFRot
    call FileWriteParameter( iounit_normal, IdConstraint_NDFRot )

    ! Save total mass of the constraint unit
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%Mass * UnitMass * 1000._RK * NAvogadro, this%Mass
    call FileWriteParameter( iounit_normal, IdConstraint_Mass )

    ! Save moments of inertia
    if( this%NDFRot > 0 ) then
      write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&       this%MOI(1) * UnitInertia * 1000._RK * NAvogadro / Angstroem**2, &
&       this%MOI(1)
      call FileWriteParameter( iounit_normal, IdConstraint_MOI1 )
      write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&       this%MOI(2) * UnitInertia * 1000._RK * NAvogadro / Angstroem**2, &
&       this%MOI(2)
      call FileWriteParameter( iounit_normal, IdConstraint_MOI2 )
      write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&       this%MOI(3) * UnitInertia * 1000._RK * NAvogadro / Angstroem**2, &
&       this%MOI(3)
      call FileWriteParameter( iounit_normal, IdConstraint_MOI3 )
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
    do i = 1, this%NLJ126
      this%Mass = this%Mass + this%SiteLJ126(i)%mass
      r(:) = r(:) + this%SiteLJ126(i)%mass * this%SiteLJ126(i)%r(:)
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
    do i = 1, this%NLJ126
      do j = 1, 3
        this%SiteLJ126(i)%r(j) = this%SiteLJ126(i)%r(j) - r(j)
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

!==============================================================!
!  Subroutine TUnit_FindMOI                                !
!==============================================================!

  subroutine TUnit_FindMOI( this )

    implicit none

    ! Declare arguments
    type(TUnit) :: this

    ! Declare local variables
    integer  :: i
    real(RK) :: moi(3, 3), rotation(3, 3), Rot2(3, 3)
    real(RK) :: qu1,qu2,qu3,qu4,quinv
    real(RK) :: T,S,SInv

    ! Calculate moment-of-inertia tensor
    moi(:, :) = 0._RK
    do i = 1, this%NLJ126
      moi(1, 1) = moi(1, 1) + this%SiteLJ126(i)%mass &
&       * ( this%SiteLJ126(i)%r(2)**2 + this%SiteLJ126(i)%r(3)**2 )
      moi(1, 2) = moi(1, 2) - this%SiteLJ126(i)%mass &
&       * this%SiteLJ126(i)%r(1) * this%SiteLJ126(i)%r(2)
      moi(1, 3) = moi(1, 3) - this%SiteLJ126(i)%mass &
&       * this%SiteLJ126(i)%r(1) * this%SiteLJ126(i)%r(3)
      moi(2, 2) = moi(2, 2) + this%SiteLJ126(i)%mass &
&       * ( this%SiteLJ126(i)%r(1)**2 + this%SiteLJ126(i)%r(3)**2 )
      moi(2, 3) = moi(2, 3) - this%SiteLJ126(i)%mass &
&       * this%SiteLJ126(i)%r(2) * this%SiteLJ126(i)%r(3)
      moi(3, 3) = moi(3, 3) + this%SiteLJ126(i)%mass &
&       * ( this%SiteLJ126(i)%r(1)**2 + this%SiteLJ126(i)%r(2)**2 )
    end do
    do i = 1, this%NCharge
      moi(1, 1) = moi(1, 1) + this%SiteCharge(i)%mass &
&       * ( this%SiteCharge(i)%r(2)**2 + this%SiteCharge(i)%r(3)**2 )
      moi(1, 2) = moi(1, 2) - this%SiteCharge(i)%mass &
&       * this%SiteCharge(i)%r(1) * this%SiteCharge(i)%r(2)
      moi(1, 3) = moi(1, 3) - this%SiteCharge(i)%mass &
&       * this%SiteCharge(i)%r(1) * this%SiteCharge(i)%r(3)
      moi(2, 2) = moi(2, 2) + this%SiteCharge(i)%mass &
&       * ( this%SiteCharge(i)%r(1)**2 + this%SiteCharge(i)%r(3)**2 )
      moi(2, 3) = moi(2, 3) - this%SiteCharge(i)%mass &
&       * this%SiteCharge(i)%r(2) * this%SiteCharge(i)%r(3)
      moi(3, 3) = moi(3, 3) + this%SiteCharge(i)%mass &
&       * ( this%SiteCharge(i)%r(1)**2 + this%SiteCharge(i)%r(2)**2 )
    end do
    do i = 1, this%NDipole
      moi(1, 1) = moi(1, 1) + this%SiteDipole(i)%mass &
&       * ( this%SiteDipole(i)%r(2)**2 + this%SiteDipole(i)%r(3)**2 )
      moi(1, 2) = moi(1, 2) - this%SiteDipole(i)%mass &
&       * this%SiteDipole(i)%r(1) * this%SiteDipole(i)%r(2)
      moi(1, 3) = moi(1, 3) - this%SiteDipole(i)%mass &
&       * this%SiteDipole(i)%r(1) * this%SiteDipole(i)%r(3)
      moi(2, 2) = moi(2, 2) + this%SiteDipole(i)%mass &
&       * ( this%SiteDipole(i)%r(1)**2 + this%SiteDipole(i)%r(3)**2 )
      moi(2, 3) = moi(2, 3) - this%SiteDipole(i)%mass &
&       * this%SiteDipole(i)%r(2) * this%SiteDipole(i)%r(3)
      moi(3, 3) = moi(3, 3) + this%SiteDipole(i)%mass &
&       * ( this%SiteDipole(i)%r(1)**2 + this%SiteDipole(i)%r(2)**2 )
    end do
    do i = 1, this%NQuadrupole
      moi(1, 1) = moi(1, 1) + this%SiteQuadrupole(i)%mass &
&       * ( this%SiteQuadrupole(i)%r(2)**2 + this%SiteQuadrupole(i)%r(3)**2 )
      moi(1, 2) = moi(1, 2) - this%SiteQuadrupole(i)%mass &
&       * this%SiteQuadrupole(i)%r(1) * this%SiteQuadrupole(i)%r(2)
      moi(1, 3) = moi(1, 3) - this%SiteQuadrupole(i)%mass &
&       * this%SiteQuadrupole(i)%r(1) * this%SiteQuadrupole(i)%r(3)
      moi(2, 2) = moi(2, 2) + this%SiteQuadrupole(i)%mass &
&       * ( this%SiteQuadrupole(i)%r(1)**2 + this%SiteQuadrupole(i)%r(3)**2 )
      moi(2, 3) = moi(2, 3) - this%SiteQuadrupole(i)%mass &
&       * this%SiteQuadrupole(i)%r(2) * this%SiteQuadrupole(i)%r(3)
      moi(3, 3) = moi(3, 3) + this%SiteQuadrupole(i)%mass &
&       * ( this%SiteQuadrupole(i)%r(1)**2 + this%SiteQuadrupole(i)%r(2)**2 )
    end do

    ! Transform to principal axes
    call eigen_find( moi(:,:), this%MOI(:), rotation(:,:) )
    call eigen_sort( this%MOI(:), rotation(:,:))
    do i = 1, this%NLJ126
      this%SiteLJ126(i)%r(:) = &
&       matmul( this%SiteLJ126(i)%r(:), rotation(:, :) )
    end do
    do i = 1, this%NCharge
      this%SiteCharge(i)%r(:) = &
&       matmul( this%SiteCharge(i)%r(:), rotation(:, :) )
    end do
    do i = 1, this%NDipole
      this%SiteDipole(i)%r(:) = &
&       matmul( this%SiteDipole(i)%r(:), rotation(:, :) )
      this%SiteDipole(i)%or(:) = &
&       matmul( this%SiteDipole(i)%or(:), rotation(:, :) )
    end do
    do i = 1, this%NQuadrupole
      this%SiteQuadrupole(i)%r(:) = &
&       matmul( this%SiteQuadrupole(i)%r(:), rotation(:, :) )
      this%SiteQuadrupole(i)%or(:) = &
&       matmul( this%SiteQuadrupole(i)%or(:), rotation(:, :) )
    end do
    if( (this%NCharge > 0).or.(this%NDipole > 0) ) &
&     this%Mue(:) = matmul( this%Mue(:), rotation(:, :) )
 
    ! Calculate inverse of rotation matrix - from body coordinate to space axes
    Rot2(:,:) = rotation


    ! Implemented according to Bronstein et al. 2008, Revision 7
    T = Rot2(1,1)+Rot2(2,2)+Rot2(3,3)+1._RK
    if (T>0) then
       S = 0.5_RK/sqrt(T)
       qu1 = 0.25_RK/S
       qu2 = (Rot2(3,2)-Rot2(2,3))*S
       qu3 = (Rot2(1,3)-Rot2(3,1))*S
       qu4 = (Rot2(2,1)-Rot2(1,2))*S
    else if ( (Rot2(1,1)>Rot2(2,2)) .and. (Rot2(1,1)>Rot2(3,3)) ) then
       S = 2._RK*sqrt(1._RK + Rot2(1,1) - Rot2(2,2) - Rot2(3,3)) ! S = 4*qu2
       SInv = 1._RK/S
       qu1 = (Rot2(3,2) - Rot2(2,3))*SInv
       qu2 = 0.25_RK*S
       qu3 = (Rot2(1,2) + Rot2(2,1))*SInv
       qu4 = (Rot2(1,3) + Rot2(3,1))*SInv
    else if (Rot2(2,2)>Rot2(3,3)) then
       S = 2._RK*sqrt(1._RK + Rot2(2,2) - Rot2(1,1) - Rot2(3,3)) ! S = 4*qu3
       SInv = 1._RK/S
       qu1 = (Rot2(1,3)-Rot2(3,1))*SInv
       qu2 = (Rot2(1,2)+Rot2(2,1))*SInv
       qu3 = 0.25_RK*S
       qu4 = (Rot2(2,3)+Rot2(3,2))*SInv
    else
       S = 2._RK*sqrt(1._RK + Rot2(3,3) - Rot2(1,1) - Rot2(2,2)) ! S = 4*qu4
       SInv = 1._RK/S
       qu1 = (Rot2(1,2)-Rot2(2,1))*SInv
       qu2 = (Rot2(1,3)+Rot2(3,1))*SInv
       qu3 = (Rot2(2,3)+Rot2(3,2))*SInv
       qu4 = 0.25_RK*S
    end if
    quinv = 1._RK / sqrt( qu1**2 + qu2**2 + qu3**2 + qu4**2 )
    qu1 = qu1 * quinv
    qu2 = qu2 * quinv
    qu3 = qu3 * quinv
    qu4 = qu4 * quinv
    this%Q0(1) = qu1
    this%Q0(2) = qu2
    this%Q0(3) = qu3
    this%Q0(4) = qu4


  contains

    subroutine jrotate( a1, a2, s, tau )

      ! Declare arguments
      real(RK), intent(in out) :: a1(:), a2(:)
      real(RK), intent(in)     :: s, tau

      ! Declare local variables
      real(RK) :: a3(size( a1 ))

      ! Rotate
      a3(:) = a1(:)
      a1(:) = a1(:) - s * (a2(:) + a1(:) * tau)
      a2(:) = a2(:) + s * (a3(:) - a2(:) * tau)

    end subroutine jrotate

    subroutine eigen_find( a, d, v )

      ! Declare arguments
      real(RK), intent(in out) :: a(3, 3)
      real(RK), intent(out)    :: d(3), v(3, 3)

      ! Declare local variables
      integer  :: i, ip, iq
      real(RK) :: c, g, h, s, sm, t, tau, theta, thresh, b(3), z(3)

      ! Compute eigenvalues and eigenvectors using Jacobi rotations
      v(:, :) = 0._RK
      do i = 1, 3
        v(i, i) = 1._RK
        b(i) = a(i, i)
      end do
      d(:) = b(:)
      do i = 1, 50
        z(:) = 0._RK
        sm = 0._RK
        do ip = 1, 2
          do iq = ip + 1, 3
            sm = sm + abs( a(ip, iq) )
          end do
        end do
        if( sm == 0._RK ) return
        thresh = merge( sm / 45._RK, 0._RK, i < 4 )
        do ip = 1, 2
          do iq = ip + 1, 3
            g = 100._RK * abs( a(ip, iq ) )
            if( &
&             (i > 4) &
&             .and. (abs( d(ip) ) + g == abs( d(ip) )) &
&             .and. (abs( d(iq) ) + g == abs( d(iq) )) &
&           ) then
              a(ip, iq) = 0._RK
            else if( abs( a(ip, iq) ) > thresh ) then
              h = d(iq) - d(ip)
              if( abs( h ) + g == abs( h ) ) then
                t = a(ip, iq) / h
              else
                theta = .5_RK * h / a(ip, iq)
                t = 1._RK / (abs( theta ) + sqrt( 1._RK + theta**2 ))
                if( theta < 0._RK ) t = -t
              end if
              c = 1._RK / sqrt( 1._RK + t**2 )
              s = t * c
              tau = s / (1._RK + c)
              h = t * a(ip, iq)
              z(ip) = z(ip) - h
              z(iq) = z(iq) + h
              d(ip) = d(ip) - h
              d(iq) = d(iq) + h
              a(ip, iq) = 0._RK
              call jrotate( a(1:ip - 1, ip), a(1:ip - 1, iq), s, tau )
              call jrotate( a(ip, ip + 1:iq - 1), a(ip + 1:iq - 1, iq), s, tau )
              call jrotate( a(ip, iq + 1:3), a(iq, iq + 1:3), s, tau )
              call jrotate( v(:, ip), v(:, iq), s, tau )
            end if
          end do
        end do
        b(:) = b(:) + z(:)
        d(:) = b(:)
      end do

    end subroutine eigen_find

    subroutine eigen_sort( d, v )

      ! Declare arguments
      real(RK), intent(in out) :: d(3), v(3, 3)
      ! Declare local variables
      integer     :: i
      real(RK)    :: p, q(3)
      integer     :: j, j1(1)
      equivalence (j, j1)

      ! Sort eigenvalues into descending order
      ! and rearrange eigenvectors correspondingly
      do i = 1, 2
        j1(:) = maxloc( d(i:3) )
        j = j + i - 1
        if( j /= i ) then
          p = d(j)
          d(j) = d(i)
          d(i) = p
          q(:) = v(:, i)
          v(:, i) = v(:, j)
          v(:, j) = q(:)
        end if
      end do

    end subroutine eigen_sort

  end subroutine TUnit_FindMOI


!==============================================================!
!  Subroutine TUnit_ReadMOI                                    !
!==============================================================!

  subroutine TUnit_ReadMOI( this )

    implicit none

    ! Declare arguments
    type(TUnit) :: this

    ! Declare local variables
    integer :: i

    ! Read moments of inertia
    this%MOI(:) = 0._RK
    if( this%NDFRot > 0 ) then
      call FileReadParameter( this%MOI(1), iounit_potmod, IdConstraint_MOI1, .false. )
      call FileReadParameter( this%MOI(2), iounit_potmod, IdConstraint_MOI2, .false. )
      if( this%NDFRot == 3 ) then
        call FileReadParameter( this%MOI(3), iounit_potmod, IdConstraint_MOI3, .false. )
      end if
    end if

    ! Convert to derived units
    do i = 1, 3
      this%MOI(i) = this%MOI(i) * .001_RK / NAvogadro * Angstroem**2
      this%MOI(i) = this%MOI(i) / UnitInertia
    end do

    if ( this%NDipole .gt. 0 .or. this%NQuadrupole .gt. 0 ) then
      this%Q0(1) = 1._RK
      this%Q0(2) = 0._RK
      this%Q0(3) = 0._RK
      this%Q0(4) = 0._RK
    end if

  end subroutine TUnit_ReadMOI

!==============================================================!
!  Subroutine TUnit_FindNDF                                    !
!==============================================================!

  subroutine TUnit_FindNDF( this )

    implicit none

    ! Declare arguments
    type(TUnit) :: this

    ! Declare local variables
    logical :: disoriented
    integer :: i

    ! Calculate number of rotation axes
    if( this%NDFRot < 0 ) then
      if( maxval( abs( this%MOI(:) ) ) > Zero ) then
        if( abs( this%MOI(3) ) > Zero ) then
          this%NDFRot = 3
        else
          this%NDFRot = 2
          this%MOI(3) = 0._RK
        end if
      else
        this%NDFRot = 0
        this%MOI(:) = 0._RK
      end if
    end if

    ! Check orientation of dipoles and quadrupoles
    if( this%NDFRot < 3 ) then
      disoriented = this%NDFRot < 2 &
&       .and. (this%NDipole > 0 .or. this%NQuadrupole > 0)
      do i = 1, this%NDipole
        disoriented = disoriented &
&         .or. ( maxval( abs( this%SiteDipole(i)%or(1:2) ) ) > Zero )
      end do
      do i = 1, this%NQuadrupole
        disoriented = disoriented &
&         .or. ( maxval( abs( this%SiteQuadrupole(i)%or(1:2) ) ) > Zero )
      end do
      if( disoriented ) &
&       call Error( 'Must specify moments of inertia manually for unit' )
    end if

    ! Calculate total number of degrees of freedom
    this%NDF = 3 + this%NDFRot

    ! Set logical flags according to the number of rotation axes
    this%isElongated = this%NDFRot > 0
    this%is3D = this%NDFRot == 3

  end subroutine TUnit_FindNDF

end module ms2_unit
