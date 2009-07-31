!==============================================================!
!  MOLECULAR SIMULATION PROGRAM MS2 Version 1.1 v12            !
!  (c) 2001 by Sergey Lishchuk, ITT                            !
!  (c) 2007 by Bernhard Eckl, ITT
!  (c) 2007 by Ekaterina Elts, TUM
!==============================================================!
!  Module ms2_molecule                                         !
!  Contains TMolecule object                                   !
!==============================================================!

#ifndef ARCH
#define ARCH    0
#define FORTRAN 90
#define MPI_VER 0
#endif

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2_molecule.F90...'
#endif

module ms2_molecule

  use ms2_global
  use ms2_site
  use ms2_idf
  use ms2_unit

!==============================================================!
!  Type TMolecule                                              !
!==============================================================!

  type TMolecule

    ! Geometry of molecule
    logical :: isElongated, is3D

    ! Number of degrees of freedom
    integer :: NDFRot, NDF

    ! Total mass of a molecule
    real(RK) :: Mass

    ! Principal moments of inertia
    real(RK) :: MOI(3)

    ! 12-6 Lennard-Jones sites
    integer :: NLJ126
    type(TSiteLJ126), pointer :: SiteLJ126(:)

    ! Coulomb sites
    integer :: NCharge
    type(TSiteCharge), pointer :: SiteCharge(:)

    ! Dipole sites
    integer :: NDipole
    type(TSiteDipole), pointer :: SiteDipole(:)

    ! Quadrupole sites
    integer :: NQuadrupole
    type(TSiteQuadrupole), pointer :: SiteQuadrupole(:)

    ! All sites
    integer :: NSite

    ! SiteIds
    integer, allocatable :: SiteIds(:)
    integer, allocatable :: LJSiteIds(:)
    integer, allocatable :: ChargeSiteIds(:)
    integer, allocatable :: QuadrupoleSiteIds(:)
    integer, allocatable :: DipoleSiteIds(:)
    integer, allocatable :: ConstraintSiteIds(:)
    integer, allocatable :: NotConstraintSiteIds(:)

    ! Bond for internal degree of freedom
    integer :: NBond
    type(TIdfBond), pointer ::IdfBond(:)

    ! Angle for internal degree of freedom
    integer :: NAngle
    type(TIdfAngle), pointer ::IdfAngle(:)

    ! Dihedral for internal degree of freedom
    integer :: NDihedral
    type(TIdfDihedral), pointer ::IdfDihedral(:)

    ! Constraint for internal degree of freedom
    integer :: NConstraint
    integer :: NNotConstraint

    ! Units of molecule
    integer, pointer :: NUnit
    integer :: NEUnit           ! number of elongated Units
    type(TUnit), pointer ::Unit(:)

    ! File name for potential model
    character(FileNameLength) :: PotModFileName

    ! Bonded Units
    integer, allocatable :: BondedUnits(:, :)

    !1-4 and 1-5 Sites in molecule for intramolecular 1-4, 1-5 nonbonded interactions
     integer, allocatable :: Int14(:, :)
     integer, allocatable :: IntLJ14(:, :), IntLJ15(:, :)
     integer, allocatable :: IntCC14(:, :), IntCC15(:, :)
     integer, allocatable :: IntCD14(:, :), IntCD15(:, :)
     integer, allocatable :: IntCQ14(:, :), IntCQ15(:, :)
     integer, allocatable :: IntDC14(:, :), IntDC15(:, :)
     integer, allocatable :: IntDD14(:, :), IntDD15(:, :)
     integer, allocatable :: IntDQ14(:, :), IntDQ15(:, :)
     integer, allocatable :: IntQC14(:, :), IntQC15(:, :)
     integer, allocatable :: IntQD14(:, :), IntQD15(:, :)
     integer, allocatable :: IntQQ14(:, :), IntQQ15(:, :)



   !Scale Factors for 1-4 nonbonded interactions
    real, allocatable :: ScaleLJ14(:), ScaleCC14(:)
    real, allocatable :: ScaleCD14(:), ScaleCQ14(:)
    real, allocatable :: ScaleDD14(:), ScaleQQ14(:)
    real, allocatable :: ScaleDC14(:), ScaleQC14(:)
    real, allocatable :: ScaleQD14(:), ScaleDQ14(:)



    ! Number of fluctuating states
    integer :: NFluct

  end type TMolecule

  interface Construct
    module procedure TMolecule_Construct
  end interface

  interface Destruct
    module procedure TMolecule_Destruct
  end interface

  interface Save
    module procedure TMolecule_Save
  end interface

  interface SaveIDF
    module procedure TMolecule_SaveIDF
  end interface

  interface FindCOM
    module procedure TMolecule_FindCOM
  end interface

  interface FindMOI
    module procedure TMolecule_FindMOI
  end interface

  interface ReadMOI
    module procedure TMolecule_ReadMOI
  end interface

  interface FindNDF
    module procedure TMolecule_FindNDF
  end interface


! Find bond  distance
  interface FindBondR
    module procedure TMolecule_FindBondR
  end interface

! Find  bond angle
  interface FindAngle
    module procedure TMolecule_FindAngle
  end interface

contains



!==============================================================!
!  Subroutine TMolecule_Construct                              !
!==============================================================!

  subroutine TMolecule_Construct( this, filename, fluctstate )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TMolecule)          :: this
    character(*), intent(in) :: filename
    integer, intent(in)      :: fluctstate

    ! Declare local variables
    integer                :: i, j, k, index, u, l, Bondend, m, n, index1, index2
    integer                :: ntypes
    integer                :: nidftypes  !number of internal degree of freedom types
    character(16)          :: stype
    character(16)          :: sidftype  !type of internal degree of freedom
    integer                :: stat
    integer                :: ncs        ! number of all constraint sites
    integer, allocatable   :: ncspu(:)   ! number of constraint sites pro unit
    logical, allocatable   :: mask(:)
    logical                :: ok, ok1, ok2, ok3, done, LJ1, LJ2, same
    logical                :: charge1, charge2, dipole1, dipole2, quadrupole1, quadrupole2
    integer                :: cc, cd, cq, dd, dq, qq, qc, qd, dc, lj
    integer, allocatable   :: S(:, :), E( :, :), Units(:, :)
    integer                :: Site1, Site2, Site3, Site4
    integer, allocatable   :: AllSites(:, :), Int14(:,:), IntLJ15(:,:)
    integer, allocatable   :: IntLJ14(:, :), SameCoord(:,:)
    integer, allocatable   :: IntCC14(:, :), IntCD14(:,:), IntCQ14(:,:)
    integer, allocatable   :: IntQQ14(:, :), IntDD14(:,:), IntDQ14(:,:)
    integer, allocatable   :: IntQC14(:, :), IntQD14(:,:), IntDC14(:,:)
    integer, allocatable   :: IntCC15(:, :), IntCD15(:,:), IntCQ15(:,:)
    integer, allocatable   :: IntQQ15(:, :), IntDD15(:,:), IntDQ15(:,:)
    integer, allocatable   :: IntQC15(:, :), IntQD15(:,:), IntDC15(:,:)
    real, allocatable      :: ScaleLJ14(:), ScaleCC14(:), ScaleCD14(:), ScaleCQ14(:)
    real, allocatable      :: ScaleQQ14(:), ScaleDD14(:), ScaleDQ14(:)
    real, allocatable      :: ScaleDC14(:), ScaleQD14(:), ScaleQC14(:)
    real, allocatable      :: CoeffLJ14(:), CoeffEl14(:)
    integer                :: Charge1Id, Charge2Id, Dipole1Id, Dipole2Id
    integer                :: Quadrupole1Id, Quadrupole2Id
    real(RK)               :: scalegeo, scalesig, scaleeps, scaleest

    ! Nullify pointers.
    nullify( this%SiteLJ126 )
    nullify( this%SiteCharge )
    nullify( this%SiteDipole )
    nullify( this%SiteQuadrupole )

    nullify( this%IdfBond )
    nullify( this%IdfAngle )
    nullify( this%IdfDihedral )

    nullify( this%Unit )
    nullify( this%NUnit )

    ! Open potential model file
    this%PotModFileName = filename
    call FileReset( iounit_potmod, this%PotModFileName )

    ! Read number of potential types
    call FileReadParameter( iounit_potmod, IdSite_ntypes )
    read( IOBuffer, * ) ntypes


    ! Zero number of sites
    this%NLJ126 = 0
    this%NCharge = 0
    this%NDipole = 0
    this%NQuadrupole = 0


    ! Zero number of idf
    this%NBond = 0
    this%NAngle = 0
    this%NDihedral = 0

    ! Zero number of Units
    this%NEUnit = 0
    allocate( this%NUnit, STAT = stat )
    call AllocationError( stat, 'number of units' )


    ! Zero number of constraint and notconstrained Units
    this%NConstraint = 0
    this%NNotConstraint = 0

    ! Zero number of Sites
    this%NSite = 0


    ! Zero number of  constraint sites and not oriented unites
    ncs = 0

    ! Loop over potential types
    do i = 1, ntypes
      call FileReadParameter( iounit_potmod, IdSite_stype )
      read( IOBuffer, * ) stype
      select case( stype )
      case( 'LJ126', 'lj126', 'LJ', 'lj' )
        call FileReadParameter( iounit_potmod, IdSite_NLJ126 )
        read( IOBuffer, * ) this%NLJ126
        if( this%NLJ126 > 0 ) then
          allocate( this%SiteLJ126(this%NLJ126), STAT = stat )
          call AllocationError( stat, 'Lennard-Jones sites', this%NLJ126 )
          do j = 1, this%NLJ126
            call Construct( this%SiteLJ126(j) )
          end do
        end if
      case( 'CHARGE', 'Charge', 'charge', 'E', 'e' )
        call FileReadParameter( iounit_potmod, IdSite_NCharge )
        read( IOBuffer, * ) this%NCharge
        if( this%NCharge > 0 ) then
          allocate( this%SiteCharge(this%NCharge), STAT = stat )
          call AllocationError( stat, 'point charge sites', this%NCharge )
          do j = 1, this%NCharge
            call Construct( this%SiteCharge(j) )
          end do
        end if
      case( 'DIPOLE', 'Dipole', 'dipole', 'D', 'd' )
        call FileReadParameter( iounit_potmod, IdSite_NDipole )
        read( IOBuffer, * ) this%NDipole
        if( this%NDipole > 0 ) then
          allocate( this%SiteDipole(this%NDipole), STAT = stat )
          call AllocationError( stat, 'dipolar sites', this%NDipole )
          do j = 1, this%NDipole
            call Construct( this%SiteDipole(j) )
          end do
        end if
      case( 'QUADRUPOLE', 'Quadrupole', 'quadrupole', 'Q', 'q' )
        call FileReadParameter( iounit_potmod, IdSite_NQuadrupole )
        read( IOBuffer, * ) this%NQuadrupole
        if( this%NQuadrupole > 0 ) then
          allocate( this%SiteQuadrupole(this%NQuadrupole), STAT = stat )
          call AllocationError( stat, 'quadrupolar sites', this%NQuadrupole )
          do j = 1, this%NQuadrupole
            call Construct( this%SiteQuadrupole(j) )
          end do
        end if
      case default
        call Error( trim( stype )//' potential is not implemented' )
      end select
    end do

    ! Read number of rotation axes
    call FileReadParameter( iounit_potmod, IdSite_NDFRot )
    read( IOBuffer, * ) stype
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
      call Error( IdSite_NDFRot//' cannot be equal to '//trim( stype ) )
    end select


    ! Find center of mass position
    call FindCOM( this )

    ! Find moments of inertia for molecule
    if( this%NDFRot < 0 ) then
      call FindMOI( this )
    else
      call ReadMOI( this )
    end if

    ! Find number of degrees of freedom
    call FindNDF( this )

    !**************************************************************************************************************
    !***************************For Chemical Potential Grad. insertion*********************************************

    ! For fluctuating particle scale parameters
    if( fluctstate > 0 ) then
      call FileReadParameter( iounit_potmod, IdNFluct )

      ! Scaling factors start in next line
      if( RootProc ) then
        do i = 1, fluctstate
          read( iounit_potmod, * ) scalegeo, scalesig, scaleeps, scaleest
        end do
      end if
#if MPI_VER > 0
      call MPI_Bcast( scalegeo, 1, MPI_DOUBLE_PRECISION, NRootProc, &
&       MPI_COMM_WORLD, ierror )
      call MPI_Bcast( scalesig, 1, MPI_DOUBLE_PRECISION, NRootProc, &
&       MPI_COMM_WORLD, ierror )
      call MPI_Bcast( scaleeps, 1, MPI_DOUBLE_PRECISION, NRootProc, &
&       MPI_COMM_WORLD, ierror )
      call MPI_Bcast( scaleest, 1, MPI_DOUBLE_PRECISION, NRootProc, &
&       MPI_COMM_WORLD, ierror )
#endif
      if( scalegeo > 1._RK .or. scalesig > 1._RK .or. &
&         scaleeps > 1._RK .or. scaleest > 1._RK ) &
&       call Error( 'Scaling factors for fluctuating particle must be lower or equal 1' )

      ! Apply scaling factors
      do i = 1, this%NLJ126
        this%SiteLJ126(i)%r = this%SiteLJ126(i)%r * scalegeo
        this%SiteLJ126(i)%sig = this%SiteLJ126(i)%sig * scalesig
        this%SiteLJ126(i)%eps = this%SiteLJ126(i)%eps * scaleeps
      end do
      do i = 1, this%NCharge
        this%SiteCharge(i)%r = this%SiteCharge(i)%r * scalegeo
        this%SiteCharge(i)%shield = this%SiteCharge(i)%shield * scalegeo
        this%SiteCharge(i)%e = this%SiteCharge(i)%e * scaleest
      end do
      do i = 1, this%NDipole
        this%SiteDipole(i)%r = this%SiteDipole(i)%r * scalegeo
        this%SiteDipole(i)%shield = this%SiteDipole(i)%shield * scalegeo
        this%SiteDipole(i)%D = this%SiteDipole(i)%D * scaleest
      end do
      do i = 1, this%NQuadrupole
        this%SiteQuadrupole(i)%r = this%SiteQuadrupole(i)%r * scalegeo
        this%SiteQuadrupole(i)%shield = this%SiteQuadrupole(i)%shield * scalegeo
        this%SiteQuadrupole(i)%Q = this%SiteQuadrupole(i)%Q * scaleest
      end do

    else if( fluctstate .eq. 0 ) then

      call FileReadParameter( iounit_potmod, IdNFluct )
      read( IOBuffer, * ) this%NFluct

    else

      this%NFluct = 0

    end if

    !**************************************************************************************************************


    ! Save used potential model
    call Save( this, fluctstate )

    ! Calculate the total number of sites
    this%NSite = this%NLJ126+this%NCharge+this%NDipole+this%NQuadrupole

    ! Create SiteIds array, if use IDF
    if (UseIntDegFreed) then
       allocate (this%SiteIds(this%NSite), STAT = stat)
       call AllocationError( stat, 'MoleculeSiteIds', this%NSite )
       allocate (this%LJSiteIds(this%NLJ126), STAT = stat)
       call AllocationError( stat, 'LJSiteIds', this%NLJ126 )
       allocate (this%ChargeSiteIds(this%NCharge), STAT = stat)
       call AllocationError( stat, 'ChargeSiteIds', this%NCharge )
       allocate (this%DipoleSiteIds(this%NDipole), STAT = stat)
       call AllocationError( stat, 'DipoleSiteIds', this%NDipole )
       allocate (this%QuadrupoleSiteIds(this%NQuadrupole), STAT = stat)
       call AllocationError( stat, 'QuadrupoleSiteIds', this%NQuadrupole )

       if( this%NLJ126 > 0 ) then
          do j = 1, this%NLJ126
            this%SiteIds(j)=this%SiteLJ126(j)%SiteId
            this%LJSiteIds(j)=this%SiteLJ126(j)%SiteId
          end do
       end if
       i=this%NLJ126
       if( this%NCharge > 0 ) then
          do j = 1+i, this%NCharge+i
            this%SiteIds(j)=this%SiteCharge(j-i)%SiteId
            this%ChargeSiteIds(j-i)=this%SiteCharge(j-i)%SiteId
          end do
       end if
       i=i+this%NCharge
       if( this%NDipole > 0 ) then
          do j = 1+i, this%NDipole+i
            this%SiteIds(j)=this%SiteDipole(j-i)%SiteId
            this%DipoleSiteIds(j-i)=this%SiteDipole(j-i)%SiteId
          end do
       end if
       i=i+this%NDipole
       if( this%NQuadrupole > 0 ) then
          do j = 1+i, this%NQuadrupole+i
            this%SiteIds(j)=this%SiteQuadrupole(j-i)%SiteId
            this%QuadrupoleSiteIds(j-i)=this%SiteQuadrupole(j-i)%SiteId
          end do
       end if
    end if


    ! Read number of internal degree of freedom types
    if (UseIntDegFreed) then
      call FileReadParameter( iounit_potmod, IdIdf_ntypes )
      read( IOBuffer, * ) nidftypes

    ! Loop over internal degree of freedom types
      do i =  1, nidftypes
        call FileReadParameter( iounit_potmod, IdIdf_stype )
        read( IOBuffer, * ) sidftype
        select case( sidftype )
        case( 'BOND', 'Bond', 'bond', 'Bonds', 'BONDS' )
          call FileReadParameter( iounit_potmod, IdIdf_NBond )
          read( IOBuffer, * ) this%NBond
          if( this%NBond > 0 ) then
            allocate( this%IdfBond(this%NBond), STAT = stat )
            call AllocationError( stat, 'Bonds for integral degrees of freedom', this%NBond )
            do j = 1, this%NBond
              call Construct( this%IdfBond(j) )
              call FindBondR(this,this%IdfBond(j), j)
            end do
          end if
        case( 'ANGLE', 'Angle', 'angle', 'Angles', 'ANGLES' )
          call FileReadParameter( iounit_potmod, IdIdf_NAngle )
          read( IOBuffer, * ) this%NAngle
          if( this%NAngle > 0 ) then
            allocate( this%IdfAngle(this%NAngle), STAT = stat )
            call AllocationError( stat, 'angles for internal degrees of freedom', this%NAngle )
            do j = 1, this%NAngle
              call Construct( this%IdfAngle(j) )
              call FindAngle(this,this%IdfAngle(j), j)
	    end do
          end if
        case( 'DIHEDRAL', 'Dihedral', 'dihedral', 'Dihedrals', 'DIHEDRALS' )
         call FileReadParameter( iounit_potmod, IdIdf_NDihedral )
          read( IOBuffer, * ) this%NDihedral
          if( this%NDihedral > 0 ) then
            allocate( this%IdfDihedral(this%NDihedral), STAT = stat )
            call AllocationError( stat, 'dihedrals for internal degrees of freedom', this%NDihedral )
            do j = 1, this%NDihedral
              call Construct( this%IdfDihedral(j) )
            end do
          end if
        case default
          call Error( trim( sidftype )//' this internal degree of freedom is not implemented' )
        end select
      end do

    ! Calculate total number of Units
      call FileReadParameter( iounit_potmod, IdUnit_NConstraint )
      read( IOBuffer, * ) this%NConstraint
      if (this%NConstraint > 0) then
        allocate (ncspu(this%NConstraint), STAT = stat)
        call AllocationError( stat, 'ncspu', this%NConstraint )
        do j = 1,this%NConstraint
            call FileReadParameter( iounit_potmod, IdConstraint_NSites )
            read( IOBuffer, * ) ncspu(j)
            ncs = ncs + ncspu(j)  ! number of sites in all constraint units
        end do
        allocate (this%ConstraintSiteIds(ncs), STAT = stat)
        call AllocationError( stat, 'ConstraintSiteIds', ncs )
      end if
      this%NNotConstraint = this%NSite-ncs  ! number of not constraint units
      this%NUnit   = this%NNotConstraint+this%NConstraint ! total number of units

    else !  No IDF,  rigid molecule
        this%NUnit = 1
        this%NConstraint = 1 ! Only one constrained unit for the whole molecule
        this%NNotConstraint = 0
    end if

    allocate(this%Unit(this%NUnit),STAT = stat)
    call AllocationError( stat, 'Units', this%NUnit)

    ! Rewind File for reading Constraints
    call FileRewind( iounit_potmod, this%PotModFileName )

    ! Construct Units

    if (UseIntDegFreed) then
       if (this%NConstraint > 0) then
          ! Construct constrained Units and create array this%ConstraintSiteIds
          k = 0
          do i = 1, this%NConstraint
              call Construct(this%Unit(i), .true., ncspu(i))
              do j=1, ncspu(i)
                this%ConstraintSiteIds(j+k)=this%Unit(i)%SiteIds(j)
                call binar_search(this%SiteLJ126%SiteId, this%Unit(i)%SiteIds(j), ok, index )
                if (ok) then
                    this%Unit(i)%NLJ126=this%Unit(i)%NLJ126+1
                    this%Unit(i)%SiteLJ126(this%Unit(i)%NLJ126)=this%SiteLJ126(index)
                    this%SiteLJ126(index)%UnitNumber = i
                end if
                if  ( .not. ok .and. this%NCharge > 0) then
                    call binar_search(this%SiteCharge%SiteId, this%Unit(i)%SiteIds(j), ok, index )
                    if (ok) then
                      this%Unit(i)%NCharge=this%Unit(i)%NCharge+1
                      this%Unit(i)%SiteCharge(this%Unit(i)%NCharge)=this%SiteCharge(index)
                      this%SiteCharge(index)%UnitNumber = i
                    end if
                end if
                if  ( .not. ok .and. this%NDipole > 0) then
                      call binar_search(this%SiteDipole%SiteId, this%Unit(i)%SiteIds(j), ok, index )
                      if (ok) then
                         this%Unit(i)%NDipole=this%Unit(i)%NDipole+1
                          this%Unit(i)%SiteDipole(this%Unit(i)%NDipole)=this%SiteDipole(index)
                          this%SiteDipole(index)%UnitNumber = i
                      end if
                end if
                if  ( .not. ok .and. this%NQuadrupole > 0) then
                      call binar_search(this%SiteQuadrupole%SiteId, this%Unit(i)%SiteIds(j), ok, index )
                      if (ok) then
                         this%Unit(i)%NQuadrupole=this%Unit(i)%NQuadrupole+1
                         this%Unit(i)%SiteQuadrupole(this%Unit(i)%NQuadrupole)=this%SiteQuadrupole(index)
                        this%SiteQuadrupole(index)%UnitNumber = i
                      end if
                end if
             end do
             k=k+ncspu(i)
          end do
 !         print *, 'Number of constraint sites=', k
       end if
       if (this%NNotConstraint > 0) then
          allocate (this%NotConstraintSiteIds(this%NNotConstraint), STAT = stat)
          call AllocationError( stat, 'NotConstraintSiteIds', this%NNotConstraint )
          k=1
          if (this%NConstraint > 0) then
            call sort_array(this%ConstraintSiteIds)
            do i = 1, this%NSite
               call binar_search(this%ConstraintSiteIds, this%SiteIds(i), ok, index)
               if (.not. ok) then
                  this%NotConstraintSiteIds(k) = this%SiteIds(i)
                  if (k<this%NNotConstraint) then
                    k=k+1
                  else
                     exit
                  end if
               end if
            end do
          else
            do i = 1, this%NSite
                  this%NotConstraintSiteIds(i) = this%SiteIds(i)
            end do
          end if
          ! Construct not constrained Units
          do i = (this%NConstraint+1), this%NUnit
             call Construct(this%Unit(i), .false., 1)
             this%Unit(i)%SiteIds=this%NotConstraintSiteIds(i-this%NConstraint)
             ! To know about this site parameters like in Constraint Unit
                call binar_search(this%SiteLJ126%SiteId, this%Unit(i)%SiteIds(1), ok, index )
                if (ok) then
                    this%Unit(i)%NLJ126=1
                    this%Unit(i)%SiteLJ126(1)=this%SiteLJ126(index)
                    this%SiteLJ126(index)%UnitNumber = i
                 end if
                if  ( .not. ok .and. this%NCharge > 0) then
                    call binar_search(this%SiteCharge%SiteId, this%Unit(i)%SiteIds(1), ok, index )
                    if (ok) then
                      this%Unit(i)%NCharge=1
                      this%Unit(i)%SiteCharge(1)=this%SiteCharge(index)
                      this%SiteCharge(index)%UnitNumber = i
                     end if
                end if
                if  ( .not. ok .and. this%NDipole > 0) then
                      call binar_search(this%SiteDipole%SiteId, this%Unit(i)%SiteIds(1), ok, index )
                      if (ok) then
                         this%Unit(i)%NDipole=1
                         this%Unit(i)%SiteDipole(1)=this%SiteDipole(index)
                         this%SiteDipole(index)%UnitNumber = i
                      end if
                end if
                if  ( .not. ok .and. this%NQuadrupole > 0) then
                      call binar_search(this%SiteQuadrupole%SiteId, this%Unit(i)%SiteIds(1), ok, index )
                      if (ok) then
                         this%Unit(i)%NQuadrupole=1
                         this%Unit(i)%SiteQuadrupole(1)=this%SiteQuadrupole(index)
                         this%SiteQuadrupole(index)%UnitNumber = i
                      end if
                end if
                ! Finish to know Unit Site's parameters
          end do
       end if
    else ! For rigid molecules
      ! construct one Constraint Unit for the whole molecule
         call Construct(this%Unit(1), .true., this%NSite)
         this%Unit(1)%NLJ126 = this%NLJ126
         do j = 1, this%NLJ126
           this%Unit(1)%SiteLJ126(j) = this%SiteLJ126(j)
           this%SiteLJ126(j)%UnitNumber = 1
         end do
         this%Unit(1)%NCharge= this%NCharge
         do j = 1, this%NCharge
           this%Unit(1)%SiteCharge(j) = this%SiteCharge(j)
           this%SiteCharge(j)%UnitNumber = 1
         end do
         this%Unit(1)%NDipole= this%NDipole
         do j = 1, this%NDipole
           this%Unit(1)%SiteDipole(j) = this%SiteDipole(j)
           this%SiteDipole(j)%UnitNumber = 1
         end do
         this%Unit(1)%NQuadrupole= this%NQuadrupole
         do j = 1, this%NQuadrupole
            this%Unit(1)%SiteQuadrupole(j) = this%SiteQuadrupole(j)
            this%SiteQuadrupole(j)%UnitNumber = 1
         end do
         this%Unit(1)%NDFRot = -1
    end if


   ! For all Units find mass, COM, moment of inertia, number of degree of freedom
   do i = 1, this%NUnit
     ! Reduction of point charges to body fixed dipole vector
     this%Unit(i)%Mue(:) = 0._RK
     if( (this%Unit(i)%NCharge > 0).or.(this%Unit(i)%NDipole > 0) ) then
      do j =1, this%Unit(i)%NCharge
        this%Unit(i)%Mue(:) = this%Unit(i)%Mue(:) + this%Unit(i)%SiteCharge(j)%r(:) * this%Unit(i)%SiteCharge(j)%e
      end do
      do j =1, this%Unit(i)%NDipole
        this%Unit(i)%Mue(:) = this%Unit(i)%Mue(:) + this%Unit(i)%SiteDipole(j)%or(:) * this%Unit(i)%SiteDipole(j)%D
      end do
!      print *, 'this%Unit(i)%Mue(:)=',this%Unit(i)%Mue(:)
    end if
    this%Unit(i)%MueSquared = sum( this%Unit(i)%Mue(:)**2 )
      call FindCOM ( this%Unit(i) )
       if( this%Unit(i)%NDFRot < 0 ) then
          call FindMOI( this%Unit(i) )
        else
          call ReadMOI( this%Unit(i) )
      end if
      call FindNDF( this%Unit(i) )
      if  (this%Unit(i)%isElongated) this%NEUnit = this%NEunit + 1
   end do

    ! Close potential model file
    call FileClose( iounit_potmod )

    if (UseIntDegFreed) then
      ! Save  used potential model with IDF
       call SaveIDF( this )
    end if


    !Consider Intramolecular interactions
    if (IntraLJEl .and. (this%NLJ126 < 4) .and. (this%NUnit < 2)) then
       call Error( 'Check *.par file, molecule too small, no intramolecular interactions can be used' )
    end if


    if (UseIntDegFreed .and. IntraLJEl) then
!      do k = 1, this%NUnit
!        call sort_array(this%Unit(k)%SiteIds)
!      end do
       do k = 1, this%NLJ126
         call sort_array(this%LJSiteIds)
       end do
       do k = 1, this%NCharge
         call sort_array(this%ChargeSiteIds)
       end do
       do k = 1, this%NDipole
         call sort_array(this%DipoleSiteIds)
       end do
       do k = 1, this%NQuadrupole
         call sort_array(this%QuadrupoleSiteIds)
       end do
   end if

   ! create list of 1-4, 1-5 interactions

   if (IntraLJEl) then
     allocate (AllSites(this%NSite, this%NSite))
     call AllocationError( stat, 'AllSites', this%NSite*this%NSite )
     allocate (SameCoord(this%NLJ126, 3))
     call AllocationError( stat, 'SameCoord', this%NLJ126*3 )
     allocate (IntLJ15((this%NSite-3)*2, 2), STAT = stat)
     call AllocationError( stat, 'Int15', (this%NSite-3)*4 )
     if (this%NCharge>0) then
       allocate (IntCC15((this%NSite-3)*2, 2), STAT = stat)
       call AllocationError( stat, 'IntCC15', (this%NSite-3)*4 )
       if (this%NDipole>0) then
         allocate (IntCD15((this%NSite-3)*2, 2), STAT = stat)
         call AllocationError( stat, 'IntCD15', (this%NSite-3)*4 )
         allocate (IntDC15((this%NSite-3)*2, 2), STAT = stat)
         call AllocationError( stat, 'IntDC15', (this%NSite-3)*4 )
       end if
       if ( this%NQuadrupole>0) then
         allocate (IntCQ15((this%NSite-3)*2, 2), STAT = stat)
         call AllocationError( stat, 'IntCQ15', (this%NSite-3)*4 )
         allocate (IntQC15((this%NSite-3)*2, 2), STAT = stat)
         call AllocationError( stat, 'IntQC15', (this%NSite-3)*4 )
       end if
     end if
     if (this%NDipole>0) then
       allocate (IntDD15((this%NSite-3)*2, 2), STAT = stat)
       call AllocationError( stat, 'IntDD15', (this%NSite-3)*4 )
       if ( this%NQuadrupole >0) then
         allocate (IntQD15((this%NSite-3)*2, 2), STAT = stat)
         call AllocationError( stat, 'IntQD15', (this%NSite-3)*4 )
         allocate (IntDQ15((this%NSite-3)*2, 2), STAT = stat)
         call AllocationError( stat, 'IntDQ15', (this%NSite-3)*4 )
       end if
     end if
     if ( this%NQuadrupole>0) then
       allocate (IntQQ15((this%NSite-3)*2, 2), STAT = stat)
       call AllocationError( stat, 'IntQQ15', (this%NSite-3)*4 )
     end if
     if (LJEl14) then
       allocate (Int14(this%NDihedral, 2), STAT = stat)
       call AllocationError( stat, 'Int14', this%NDihedral*2 )
       allocate (CoeffLJ14(this%NDihedral), STAT = stat)
       call AllocationError( stat, 'CoeffLJ14', this%NDihedral )
       allocate (CoeffEl14(this%NDihedral), STAT = stat)
       call AllocationError( stat, 'CoeffEl14', this%NDihedral )
       allocate (IntLJ14(this%NDihedral, 2), STAT = stat)
       call AllocationError( stat, 'IntLJ14', this%NDihedral*2 )
       allocate (ScaleLJ14(this%NDihedral), STAT = stat)
       call AllocationError( stat, 'ScaleLJ14', this%NDihedral )
       if (this%NCharge>0) then
         allocate (IntCC14(this%NDihedral, 2), STAT = stat)
         call AllocationError( stat, 'IntCC14', this%NDihedral*2 )
         allocate (ScaleCC14(this%NDihedral), STAT = stat)
         call AllocationError( stat, 'ScaleCC14', this%NDihedral )
         if (this%NDipole>0) then
           allocate (IntCD14(this%NDihedral, 2), STAT = stat)
           call AllocationError( stat, 'IntCD14', this%NDihedral*2 )
           allocate (ScaleCD14(this%NDihedral), STAT = stat)
           call AllocationError( stat, 'ScaleCD14', this%NDihedral )
           allocate (IntDC14(this%NDihedral, 2), STAT = stat)
           call AllocationError( stat, 'IntDC14', this%NDihedral*2 )
           allocate (ScaleDC14(this%NDihedral), STAT = stat)
           call AllocationError( stat, 'ScaleDC14', this%NDihedral )
         end if
         if ( this%NQuadrupole>0) then
           allocate (IntCQ14(this%NDihedral, 2), STAT = stat)
           call AllocationError( stat, 'IntCQ14', this%NDihedral*2 )
           allocate (ScaleCQ14(this%NDihedral), STAT = stat)
           call AllocationError( stat, 'ScaleCQ14', this%NDihedral )
           allocate (IntQC14(this%NDihedral, 2), STAT = stat)
           call AllocationError( stat, 'IntQC14', this%NDihedral*2 )
           allocate (ScaleQC14(this%NDihedral), STAT = stat)
           call AllocationError( stat, 'ScaleQC14', this%NDihedral )
         end if
       end if
       if (this%NDipole>0) then
         allocate (IntDD14(this%NDihedral, 2), STAT = stat)
         call AllocationError( stat, 'IntDD14', this%NDihedral*2 )
         allocate (ScaleDD14(this%NDihedral), STAT = stat)
         call AllocationError( stat, 'ScaleDD14', this%NDihedral )
         if ( this%NQuadrupole>0) then
           allocate (IntDQ14(this%NDihedral, 2), STAT = stat)
           call AllocationError( stat, 'IntDQ14', this%NDihedral*2 )
           allocate (ScaleDQ14(this%NDihedral), STAT = stat)
           call AllocationError( stat, 'ScaleDQ14', this%NDihedral )
           allocate (IntQD14(this%NDihedral, 2), STAT = stat)
           call AllocationError( stat, 'IntQD14', this%NDihedral*2 )
           allocate (ScaleQD14(this%NDihedral), STAT = stat)
           call AllocationError( stat, 'ScaleQD14', this%NDihedral )
         end if
       end if
       if ( this%NQuadrupole>0) then
         allocate (IntQQ14(this%NDihedral, 2), STAT = stat)
         call AllocationError( stat, 'IntQQ14', this%NDihedral*2 )
         allocate (ScaleQQ14(this%NDihedral), STAT = stat)
         call AllocationError( stat, 'ScaleQQ14', this%NDihedral )
       end if
     end if

! Initialisierung
     do i=1, this%NSite
       do j=1, this%NSite
       if (i==j) then
         AllSites(i,j)=0
       else
         AllSites(i,j) = 1
       end if
      end do
     end do

    do i=1, this%NLJ126
      do j=1,3
        SameCoord(i,j)=0
      end do
    end do


     do i=1, this%NDihedral
       Site1=this%IdfDihedral(i)%SiteId1
       Site2=this%IdfDihedral(i)%SiteId2
       Site3=this%IdfDihedral(i)%SiteId3
       Site4=this%IdfDihedral(i)%SiteId4
       AllSites(Site1,Site2)=0
       AllSites(Site1,Site3)=0
       AllSites(Site1,Site4)=0
       AllSites(Site2,Site3)=0
       AllSites(Site2,Site4)=0
       AllSites(Site3,Site4)=0
       AllSites(Site2,Site1)=0
       AllSites(Site3,Site1)=0
       AllSites(Site4,Site1)=0
       AllSites(Site3,Site2)=0
       AllSites(Site4,Site2)=0
       AllSites(Site4,Site3)=0
     end do

     do i = 1, this%NUnit
        if (this%Unit(i)%NSites>1) then
           do j=1, this%Unit(i)%NSites
              AllSites(this%Unit(i)%SiteIds(j),this%Unit(i)%SiteIds(:))=0
           end do
         end if
      end do


     ! Find LJ and Charge Sites with the same coordinates (for1,4-1,5 interactions)
     do i=1, this%NLJ126
       do j=1, this%NCharge
         call compare_coord(this%SiteLJ126(i)%r(:), this%SiteCharge(j)%r(:), same)
         if (same) then
           SameCoord(this%SiteLJ126(i)%SiteId,1)=this%SiteCharge(j)%SiteId
           AllSites(this%SiteLJ126(i)%SiteId, this%SiteCharge(j)%SiteId)=0
           do k=1, this%NSite
             if (AllSites(this%SiteLJ126(i)%SiteId, k)== 0) then
                AllSites(this%SiteCharge(j)%SiteId, k)=0
                AllSites(k, this%SiteCharge(j)%SiteId)=0
             end if
           end do
        end if
       end do
     end do

    ! Find LJ and Dipole Sites with the same coordinates (for1,4-1,5 interactions)
    do i=1, this%NLJ126
       do j=1, this%NDipole
         call compare_coord(this%SiteLJ126(i)%r(:), this%SiteDipole(j)%r(:), same)
         if (same) then
           SameCoord(this%SiteLJ126(i)%SiteId,2)=this%SiteDipole(j)%SiteId
           AllSites(this%SiteLJ126(i)%SiteId, this%SiteDipole(j)%SiteId)=0
           do k=1, this%NSite
             if (AllSites(this%SiteLJ126(i)%SiteId, k)== 0) then
                AllSites(this%SiteDipole(j)%SiteId, k)=0
                AllSites(k, this%SiteDipole(j)%SiteId)=0
             end if
           end do
         end if
       end do
     end do

    ! Find LJ and Quadrupole Sites with the same coordinates (for1,4-1,5 interactions)
    do i=1, this%NLJ126
       do j=1, this%NQuadrupole
         call compare_coord(this%SiteLJ126(i)%r(:), this%SiteQuadrupole(j)%r(:), same )
         if (same) then
           SameCoord(this%SiteLJ126(i)%SiteId,3)=this%SiteQuadrupole(j)%SiteId
           AllSites(this%SiteLJ126(i)%SiteId, this%SiteQuadrupole(j)%SiteId)=0
           do k=1, this%NSite
             if (AllSites(this%SiteLJ126(i)%SiteId, k)== 0) then
                AllSites(this%SiteQuadrupole(j)%SiteId, k)=0
                AllSites(k, this%SiteQuadrupole(j)%SiteId)=0
             end if
           end do
          end if
       end do
     end do


     lj=1
     cc=1
     dd=1
     qq=1
     cd=1
     cq=1
     dq=1
     dc=1
     qc=1
     qd=1
     do i=1, this%NSite
       do j=i+1, this%NSite
         if (AllSites(i,j)==1) then
           LJ1 =.false.
           LJ2 =.false.
           charge1 = .false.
           charge2 = .false.
           dipole1 = .false.
           dipole2 = .false.
           quadrupole1 = .false.
           quadrupole2 = .false.
           call binar_search(this%LJSiteIds(:), i, LJ1, index)
           call binar_search(this%LJSiteIds(:), j, LJ2, index)
           if (LJ1 .and. LJ2) then
             IntLJ15(lj,1)=i
!             print *, 'IntLJ15(lj,1)=', i
             IntLJ15(lj,2)=j
!             print *, 'IntLJ15(lj,2)=', j
             lj=lj+1
           else if (.not. LJ1 .and. .not. LJ2) then
             call binar_search(this%ChargeSiteIds(:), i, charge1, index)
              if (.not. charge1) then
                call binar_search(this%DipoleSiteIds(:), i, dipole1, index)
                if (.not. dipole1) then
                  call binar_search(this%QuadrupoleSiteIds(:), i, quadrupole1, index)
                 end if
               end if
             call binar_search(this%ChargeSiteIds(:), j, charge2, index)
              if (.not. charge2) then
                call binar_search(this%DipoleSiteIds(:), i, dipole2, index)
                if (.not. dipole2) then
                  call binar_search(this%QuadrupoleSiteIds(:), i, quadrupole2, index)
                 end if
               end if
             if (charge1) then
               if (charge2) then
                 IntCC15(cc,1)= i
                 IntCC15(cc,2)= j
!                print *, 'CC15, 1=', IntCC15(cc,1)
!                print *, 'CC15, 2=', IntCC15(cc,2)
                 cc=cc+1
               else if (dipole2) then
                 IntCD15(cd,1)= i
                 IntCD15(cd,2)= j
!                print *, 'CD15, 1=', IntCD15(cd,1)
!                print *, 'CD15, 2=', IntCD15(cd,2)
                 cd=cd+1
               else !quadrupole2
                 IntCQ15(cq,1)=i
                 IntCQ15(cq,2)=j
                 cq=cq+1
               end if
             else if (dipole1) then
               if (charge2) then
                 IntDC15(dc,1)= i
                 IntDC15(dc,2)= j
                 dc=dc+1
               else if (dipole2) then
                 IntDD15(dd,1)= i
                 IntDD15(dd,2)= j
                 dd=dd+1
               else !quadrupole2
                 IntDQ15(dq,1)=i
                 IntDQ15(dq,2)=j
                 dq=dq+1
               end if
             else ! quadrupole1
               if (charge2)then
                 IntQC15(qc,1)= i
                 IntQC15(qc,2)= j
                 qc=qc+1
               else if (dipole2)then
                 IntQD15(qd,1)= i
                 IntQD15(qd,2)= j
                 qd=qd+1
               else ! quadrupole2
                 IntQQ15(qq,1)=i
                 IntQQ15(qq,2)=j
                 qq=qq+1
               end if
             end if
           else
             AllSites(i,j)=0
           end if
         end if
       end do
     end do



!     print *, 'Number of 1,5 LJ=', lj-1
!     print *, 'Number of 1,5 CC=', cc-1

     allocate (this%IntLJ15(lj-1, 2), STAT = stat)
     call AllocationError( stat, 'this%IntLJ15', (lj-1)*2 )
     this%IntLJ15 = IntLJ15(1:lj-1,:)

     if (this%NCharge>0) then
       allocate (this%IntCC15(cc-1, 2), STAT = stat)
       call AllocationError( stat, 'this%IntCC15', (cc-1)*2 )
       this%IntCC15 = IntCC15(1:cc-1,:)
!       print *, "this%IntCC15=", this%IntCC15(:,:)
       if (this%NDipole>0) then
         allocate (this%IntCD15(cd-1, 2), STAT = stat)
         call AllocationError( stat, 'this%IntCD15', (cd-1)*2 )
         this%IntCD15 = IntCD15(1:cd-1,:)
         allocate (this%IntDC15(dc-1, 2), STAT = stat)
         call AllocationError( stat, 'this%IntDC15', (dc-1)*2 )
         this%IntDC15 = IntDC15(1:dc-1,:)
       end if
       if (this%NQuadrupole>0) then
         allocate (this%IntCQ15(cq-1, 2), STAT = stat)
         call AllocationError( stat, 'this%IntCQ15', (cq-1)*2 )
         this%IntCQ15 = IntCQ15(1:cq-1,:)
         allocate (this%IntQC15(qc-1, 2), STAT = stat)
         call AllocationError( stat, 'this%IntQC15', (qc-1)*2 )
         this%IntQC15 = IntQC15(1:qc-1,:)
        end if
     end if
     if (this%NDipole>0) then
       allocate (this%IntDD15(dd-1, 2), STAT = stat)
       call AllocationError( stat, 'this%IntDD15', (dd-1)*2 )
       this%IntDD15 = IntDD15(1:dd-1,:)
       if (this%NQuadrupole>0) then
         allocate (this%IntDQ15(dq-1, 2), STAT = stat)
         call AllocationError( stat, 'this%IntDQ15', (dq-1)*2 )
         this%IntDQ15 = IntDQ15(1:dq-1,:)
         allocate (this%IntQD15(qd-1, 2), STAT = stat)
         call AllocationError( stat, 'this%IntQD15', (qd-1)*2 )
         this%IntQD15 = IntQD15(1:qd-1,:)
       end if
     end if
     if (this%NQuadrupole>0) then
       allocate (this%IntQQ15(qq-1, 2), STAT = stat)
       call AllocationError( stat, 'this%IntQQ15', (qq-1)*2 )
       this%IntQQ15 = IntQQ15(1:qq-1,:)
      end if

!     print *, "this%IntLJ15=", this%IntLJ15(:,:)


     if (LJEl14) then
       k=1
       do i=1, this%NDihedral
         if (this%IdfDihedral(i)%multi>=0) then  !If proper dihedral
           Site1=this%IdfDihedral(i)%SiteId1
           Site4=this%IdfDihedral(i)%SiteId4
           if (Site1>Site4) then
             Site1=this%IdfDihedral(i)%SiteId4
             Site4=this%IdfDihedral(i)%SiteId1
           end if
           if (k>1) then
             call binar_search(Int14(1:k,1), Site1, ok1, index)
             if (.not. ok1) then
                Int14(k,1)=Site1
                Int14(k,2)=Site4
                CoeffLJ14(k)=this%IdfDihedral(i)%ScaleLJ14
                CoeffEl14(k)=this%IdfDihedral(i)%ScaleEl14
                k=k+1
             else if (Int14(index, 2) .ne. Site4) then
                Int14(k,1)=Site1
                Int14(k,2)=Site4
                CoeffLJ14(k)=this%IdfDihedral(i)%ScaleLJ14
                CoeffEl14(k)=this%IdfDihedral(i)%ScaleEl14
                k=k+1
             end if
           else
             Int14(k,1)=Site1
             Int14(k,2)=Site4
             CoeffLJ14(k)=this%IdfDihedral(i)%ScaleLJ14
             CoeffEl14(k)=this%IdfDihedral(i)%ScaleEl14
             k=k+1
           end if
         end if
       end do

       lj = 1
       cc = 1
       cq = 1
       cd = 1
       dc = 1
       dd = 1
       dq = 1
       qc = 1
       qd = 1
       qq = 1
       do i =1, k-1
         call binar_search(this%LJSiteIds(:), Int14(i,1), LJ1, index1)
         call binar_search(this%LJSiteIds(:), Int14(i,2), LJ2, index2)
         if (LJ1 .and. LJ2) then
           IntLJ14(lj,1)=Int14(i,1)
           IntLJ14(lj,2)=Int14(i,2)
           ScaleLJ14(lj)=CoeffLJ14(i)
!          ScaleLJ14(lj)=this%IdfDihedral(i)%ScaleLJ14
!          print *, 'IntLJ14, 1=', IntLJ14(lj, 1)
!          print *, 'IntLJ14, 2=', IntLJ14(lj, 2)
!          print *, 'ScaleLJ14=', ScaleLJ14(lj)
           lj = lj+1
           Charge1Id = SameCoord(this%LJSiteIds(index1),1)
           Charge2Id = SameCoord(this%LJSiteIds(index2),1)
           Dipole1Id = SameCoord(this%LJSiteIds(index1),2)
           Dipole2Id = SameCoord(this%LJSiteIds(index2),2)
           Quadrupole1Id = SameCoord(this%LJSiteIds(index1),3)
           Quadrupole2Id = SameCoord(this%LJSiteIds(index1),3)
           if (Charge1Id > 0) then
             if ( Charge2Id > 0) then
               IntCC14(cc,1) = Charge1Id
               IntCC14(cc,2) = Charge2Id
               ScaleCC14(cc)=CoeffEl14(i)
!              ScaleCC14(cc)=this%IdfDihedral(i)%ScaleEl14
!              print *, 'IntCC14, 1=', IntCC14(cc, 1)
!              print *, 'IntCC14, 2=', IntCC14(cc, 2)
!              print *, 'ScaleCC14=', ScaleCC14(cc)
               cc = cc+1
             else if ( Dipole2Id > 0) then
!               IntCD14(cd,1) = min(Charge1Id, Dipole2Id)
!               IntCD14(cd,2) = max(Charge1Id, Dipole2Id)
               IntCD14(cd,1) = Charge1Id
               IntCD14(cd,2) = Dipole2Id
               ScaleCD14(cd)=CoeffEl14(i)
!               ScaleCD14(cd)=this%IdfDihedral(i)%ScaleEl14
               cd = cd+1
             else
                if ( Quadrupole2Id > 0) then
                  IntCQ14(cq,1) = Charge1Id
                  IntCQ14(cq,2) = Quadrupole2Id
                  ScaleCQ14(cq)=CoeffEl14(i)
!                  ScaleCQ14(cq)=this%IdfDihedral(i)%ScaleEl14
                  cq = cq+1
                end if
             end if
           end if
           if (Dipole1Id > 0) then
             if ( Charge2Id > 0) then
               IntDC14(dc,1) = Dipole1Id
               IntDC14(dc,2) = Charge2Id
               ScaleDC14(dc)=CoeffEl14(i)
!               ScaleDC14(dc)=this%IdfDihedral(i)%ScaleEl14
               dc = dc+1
             else if ( Dipole2Id > 0) then
               IntDD14(dd,1) = Dipole1Id
               IntDD14(dd,2) = Dipole2Id
               ScaleDD14(dd)=CoeffEl14(i)
!               ScaleDD14(dd)=this%IdfDihedral(i)%ScaleEl14
               dd = dd+1
             else
                if (Quadrupole2Id > 0) then
                  IntDQ14(dq,1) = Dipole1Id
                  IntDQ14(dq,2) = Quadrupole2Id
                  ScaleDQ14(dq)=CoeffEl14(i)
!                  ScaleDQ14(dq)=this%IdfDihedral(i)%ScaleEl14
                  dq = dq+1
                end if
             end if
           end if
           if (Quadrupole1Id) then
             if ( Charge2Id > 0) then
               IntQC14(qc,1) = Quadrupole1Id
               IntQC14(qc,2) = Charge2Id
               ScaleQC14(qc) = CoeffEl14(i)
!               ScaleQC14(qc)=this%IdfDihedral(i)%ScaleEl14
               qc = qc+1
             else if ( Dipole2Id > 0) then
               IntQD14(qd,1) = Quadrupole1Id
               IntQD14(qd,2) = Dipole2Id
               ScaleQD14(qd) = CoeffEl14(i)
!               ScaleQD14(qd)=this%IdfDihedral(i)%ScaleEl14
               qd = qd+1
             else
                if ( Quadrupole2Id > 0) then
                  IntQQ14(qq,1) = Quadrupole1Id
                  IntQQ14(qq,2) = Quadrupole2Id
                  ScaleQQ14(qq)= CoeffEl14(i)
!                  ScaleQQ14(qq)=this%IdfDihedral(i)%ScaleEl14
                  qq = qq+1
                end if
             end if
           end if
         end if
       end do


       allocate (this%IntLJ14(lj-1, 2), STAT = stat)
       call AllocationError( stat, 'this%IntLJ14', (lj-1)*2 )
       this%IntLJ14 = IntLJ14(1:lj-1,:)
!       print *, 'IntLJ14=', this%IntLJ14(:,:)
       allocate (this%ScaleLJ14(lj-1), STAT = stat)
       call AllocationError( stat, 'ScaleLJ14', lj-1 )
       this%ScaleLJ14 = ScaleLJ14(1:lj-1)
!       print *, 'ScaleLJ14=', this%ScaleLJ14(:)


       if (this%NCharge>0) then
         allocate (this%IntCC14(cc-1, 2), STAT = stat)
         call AllocationError( stat, 'this%IntCC14', (cc-1)*2 )
         this%IntCC14 = IntCC14(1:cc-1,:)
!         print *, 'IntCC14=', this%IntCC14(:,:)
         allocate (this%ScaleCC14(cc-1), STAT = stat)
         call AllocationError( stat, 'ScaleCC14', cc-1 )
         this%ScaleCC14 = ScaleCC14(1:cc-1)
!         print *, 'ScaleCC14=', this%ScaleCC14(:)
         if (this%NDipole>0) then
           allocate (this%IntCD14(cd-1, 2), STAT = stat)
           call AllocationError( stat, 'this%IntCD14', (cd-1)*2 )
           this%IntCD14 = IntCD14(1:cd-1,:)
           allocate (this%ScaleCD14(cd-1), STAT = stat)
           call AllocationError( stat, 'ScaleCD14', cd-1 )
           this%ScaleCD14 = ScaleCD14(1:cd-1)
           allocate (this%IntDC14(dc-1, 2), STAT = stat)
           call AllocationError( stat, 'this%IntDC14', (dc-1)*2 )
           this%IntDC14 = IntDC14(1:dc-1,:)
           allocate (this%ScaleDC14(dc-1), STAT = stat)
           call AllocationError( stat, 'ScaleDC14', dc-1 )
           this%ScaleDC14 = ScaleDC14(1:dc-1)
         end if
         if (this%NQuadrupole>0) then
           allocate (this%IntQC14(qc-1, 2), STAT = stat)
           call AllocationError( stat, 'this%IntQC14', (qc-1)*2 )
           this%IntQC14 = IntQC14(1:qc-1,:)
           allocate (this%ScaleQC14(qc-1), STAT = stat)
           call AllocationError( stat, 'ScaleQC14', qc-1 )
           this%ScaleQC14 = ScaleQC14(1:qc-1)
           allocate (this%IntCQ14(cq-1, 2), STAT = stat)
           call AllocationError( stat, 'this%IntCQ14', (cq-1)*2 )
           this%IntCQ14 = IntCQ14(1:cq-1,:)
           allocate (this%ScaleCQ14(cq-1), STAT = stat)
           call AllocationError( stat, 'ScaleCQ14', cq-1 )
           this%ScaleCQ14 = ScaleCQ14(1:cq-1)
         end if
       end if
       if (this%NDipole>0) then
         allocate (this%IntDD14(dd-1, 2), STAT = stat)
         call AllocationError( stat, 'this%IntDD14', (dd-1)*2 )
         this%IntDD14 = IntDD14(1:dd-1,:)
         allocate (this%ScaleDD14(dd-1), STAT = stat)
         call AllocationError( stat, 'ScaleDD14', dd-1 )
         this%ScaleDD14 = ScaleDD14(1:dd-1)
         if (this%NQuadrupole>0) then
           allocate (this%IntDQ14(dq-1, 2), STAT = stat)
           call AllocationError( stat, 'this%IntDQ14', (dq-1)*2 )
           this%IntDQ14 = IntDQ14(1:dq-1,:)
           allocate (this%ScaleDQ14(dq-1), STAT = stat)
           call AllocationError( stat, 'ScaleDQ14', dq-1 )
           this%ScaleDQ14 = ScaleDQ14(1:dq-1)
           allocate (this%IntQD14(qd-1, 2), STAT = stat)
           call AllocationError( stat, 'this%IntQD14', (qd-1)*2 )
           this%IntQD14 = IntQD14(1:qd-1,:)
           allocate (this%ScaleQD14(qd-1), STAT = stat)
           call AllocationError( stat, 'ScaleQD14', qd-1 )
           this%ScaleQD14 = ScaleQD14(1:qd-1)
         end if
       end if
       if (this%NQuadrupole>0) then
         allocate (this%IntQQ14(qq-1, 2), STAT = stat)
         call AllocationError( stat, 'this%IntQQ14', (qq-1)*2 )
         this%IntQQ14 = IntQQ14(1:qq-1,:)
         allocate (this%ScaleQQ14(qq-1), STAT = stat)
         call AllocationError( stat, 'ScaleQQ14', qq-1 )
         this%ScaleQQ14 = ScaleQQ14(1:qq-1)
       end if

     end if
   end if

   contains

    subroutine sort_array(array)
      ! Declare arguments
      integer, dimension(:), intent( inout ) :: array
      ! Declare local variables
      integer :: i,j

      do i = 1, size(array)
         do j = 1, size(array)
            if (array(j)>array(i))  call change(array(i), array(j))
         end do
      end do
    end subroutine sort_array

   subroutine change(a, b)
      ! Declare arguments
      integer, intent( inout ) :: a, b
      ! Declare local variables
      integer :: temp

      temp = a
      a = b
      b = temp
   end subroutine change

    subroutine binar_search (array, Id, treffer, index)

      ! Declare arguments
      integer, dimension(:), intent( in ) :: array
      integer, intent( in )               :: Id
      logical, intent( out )              :: treffer
      integer, intent( out )              :: index

      ! Declare local variables
      integer                             :: anfang, ende, mitte

      anfang = 1
      ende = size (array)
      do
         if ( anfang == ende ) exit
         mitte = (anfang + ende)*0.5
         if ( id <= array(mitte) ) then
           ende = mitte
         else
           anfang = mitte + 1
         end if
      end do
      index = anfang
      treffer = (id == array(index))

    end subroutine binar_search

    subroutine compare_coord(array1, array2, same)

    ! Declare arguments
    real(RK), dimension(3),intent(in) ::array1
    real(RK), dimension(3), intent(in) ::array2
    logical, intent(out)              :: same

    ! Declare local variables

     if (array1(1)==array2(1) .and. array1(2)==array2(2) .and. array1(3)==array2(3) ) then
     	same = .true.
     else
     	same = .false.
    end if

    end subroutine compare_coord

    subroutine FindEdgeFrom( unit, bondedunit, E, n )

      ! Declare arguments
      integer, dimension(this%NBond, 3), intent( in) :: bondedunit
      integer, intent( in )                           :: unit
      integer, intent( inout )                        :: n
      integer, dimension(this%NBond, 2),  intent( out ):: E

      ! Declare local variables
      integer                                 :: index
      logical                                 :: ok

      do index = 1, this%NBond
         if ( unit == bondedunit(index, 1)) then
            n = n + 1
            E(n, 1) = index
            E(n, 2) = bondedunit(index,2)
         else if (unit == bondedunit(index, 2)) then
            n = n + 1
            E(n, 1) = index
            E(n, 2) = bondedunit(index, 1)
         end if
      end do

    end subroutine FindEdgeFrom

  end subroutine TMolecule_Construct


!==============================================================!
!  Subroutine TMolecule_Destruct                               !
!==============================================================!

  subroutine TMolecule_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TMolecule) :: this

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

  end subroutine TMolecule_Destruct



!==============================================================!
!  Subroutine TMolecule_Save                                   !
!==============================================================!

  subroutine TMolecule_Save( this, fluctstate )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TMolecule)     :: this
    integer, intent(in) :: fluctstate

    ! Declare local variables
    character(FileNameLength) :: filename
    integer                   :: ntypes
    integer                   :: i

    ! Open file
    if( fluctstate < 1 ) then
      filename = trim( this%PotModFileName )//NormalizedPotModExtension
    else
      write( filename, '(A, ".", I0)') &
&       trim( this%PotModFileName )//NormalizedPotModExtension, fluctstate
    end if
    call FileRewrite( iounit_normal, filename )

    ! Save number of potential types
    ntypes = 0
    if( this%NLJ126 > 0 ) ntypes = ntypes + 1
    if( this%NCharge > 0 ) ntypes = ntypes + 1
    if( this%NDipole > 0 ) ntypes = ntypes + 1
    if( this%NQuadrupole > 0 ) ntypes = ntypes + 1
    write( IOBuffer, '(I2)' ) ntypes
    call FileWriteParameter( iounit_normal, IdSite_ntypes )

    ! Save Lennard-Jones sites
    if( this%NLJ126 > 0 ) then
      call FileWriteBlank( iounit_normal )
      write( IOBuffer, '(X, A)' ) 'LJ126'
      call FileWriteParameter( iounit_normal, IdSite_stype )
      write( IOBuffer, '(I2)' ) this%NLJ126
      call FileWriteParameter( iounit_normal, IdSite_NLJ126 )
      do i = 1, this%NLJ126
        call FileWriteBlank( iounit_normal )
        call Save( this%SiteLJ126(i) )
      end do
    end if

    ! Save point charge sites
    if( this%NCharge > 0 ) then
      call FileWriteBlank( iounit_normal )
      write( IOBuffer, '(X, A)' ) 'Charge'
      call FileWriteParameter( iounit_normal, IdSite_stype )
      write( IOBuffer, '(I2)' ) this%NCharge
      call FileWriteParameter( iounit_normal, IdSite_NCharge )
      do i = 1, this%NCharge
        call FileWriteBlank( iounit_normal )
        call Save( this%SiteCharge(i) )
      end do
    end if

    ! Save point dipole sites
    if( this%NDipole > 0 ) then
      call FileWriteBlank( iounit_normal )
      write( IOBuffer, '(X, A)' ) 'Dipole'
      call FileWriteParameter( iounit_normal, IdSite_stype )
      write( IOBuffer, '(I2)' ) this%NDipole
      call FileWriteParameter( iounit_normal, IdSite_NDipole )
      do i = 1, this%NDipole
        call FileWriteBlank( iounit_normal )
        call Save( this%SiteDipole(i) )
      end do
    end if

    ! Save point quadrupole sites
    if( this%NQuadrupole > 0 ) then
      call FileWriteBlank( iounit_normal )
      write( IOBuffer, '(X, A)' ) 'Quadrupole'
      call FileWriteParameter( iounit_normal, IdSite_stype )
      write( IOBuffer, '(I2)' ) this%NQuadrupole
      call FileWriteParameter( iounit_normal, IdSite_NQuadrupole )
      do i = 1, this%NQuadrupole
        call FileWriteBlank( iounit_normal )
        call Save( this%SiteQuadrupole(i) )
      end do
    end if


    ! Save number of rotation axes of molecule
    call FileWriteBlank( iounit_normal )
    write( IOBuffer, '(I2)' ) this%NDFRot
    call FileWriteParameter( iounit_normal, IdSite_NDFRot )

    ! Save total mass of the molecule
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%Mass * UnitMass * 1000._RK * NAvogadro, this%Mass
    call FileWriteParameter( iounit_normal, IdSite_Mass )

    ! Save moments of inertia of molecule
    if( this%NDFRot > 0 ) then
      write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&       this%MOI(1) * UnitInertia * 1000._RK * NAvogadro / Angstroem**2, &
&       this%MOI(1)
      call FileWriteParameter( iounit_normal, IdSite_MOI1 )
      write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&       this%MOI(2) * UnitInertia * 1000._RK * NAvogadro / Angstroem**2, &
&       this%MOI(2)
      call FileWriteParameter( iounit_normal, IdSite_MOI2 )
      write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&       this%MOI(3) * UnitInertia * 1000._RK * NAvogadro / Angstroem**2, &
&       this%MOI(3)
      call FileWriteParameter( iounit_normal, IdSite_MOI3 )
    end if



    ! Close file
    call FileClose( iounit_normal )

    ! Update log file
    write( IOBuffer, '("Normalized potential model for ", A, &
&     " saved to file <", A, ">")' ) &
&     trim( this%PotModFileName ), trim( filename )
    call LogWrite

  end subroutine TMolecule_Save


!==============================================================!
!  Subroutine TMolecule_SaveIDF                                !
!==============================================================!

  subroutine TMolecule_SaveIDF( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TMolecule) :: this

    ! Declare local variables
    character(FileNameLength)      :: filename
    integer                        :: nidftypes
    integer                        :: i

    ! Open file
    filename = trim( this%PotModFileName )//NormalizedPotModExtension
    call FileAppend( iounit_normal, filename )

    ! Save information about Idf
    ! Save number of potential types
    call FileWriteBlank( iounit_normal )
    nidftypes = 0
    if( this%NBond > 0 ) nidftypes = nidftypes + 1
    if( this%NAngle > 0 ) nidftypes = nidftypes + 1
    if( this%NDihedral > 0 ) nidftypes = nidftypes + 1
    write( IOBuffer, '(I2)' ) nidftypes
    call FileWriteParameter( iounit_normal, IdIdf_ntypes )

    ! Save Bonds
    if( this%NBond > 0 ) then
      call FileWriteBlank( iounit_normal )
      write( IOBuffer, '(X, A)' ) 'Bond'
      call FileWriteParameter( iounit_normal, IdIdf_stype )
      write( IOBuffer, '(I2)' ) this%NBond
      call FileWriteParameter( iounit_normal, IdIdf_NBond )
      do i = 1, this%NBond
        call FileWriteBlank( iounit_normal )
        call Save( this%IdfBond(i) )
      end do
    end if

   ! Save Angles
   if( this%NAngle > 0 ) then
      call FileWriteBlank( iounit_normal )
      write( IOBuffer, '(X, A)' ) 'Angle'
      call FileWriteParameter( iounit_normal, IdIdf_stype )
      write( IOBuffer, '(I2)' ) this%NAngle
      call FileWriteParameter( iounit_normal, IdIdf_NAngle )
      do i = 1, this%NAngle
        call FileWriteBlank( iounit_normal )
        call Save( this%IdfAngle(i) )
      end do
   end if

   ! Save Dihedrals
   if( this%NDihedral > 0 ) then
      call FileWriteBlank( iounit_normal )
      write( IOBuffer, '(X, A)' ) 'Dihedral'
      call FileWriteParameter( iounit_normal, IdIdf_stype )
      write( IOBuffer, '(I2)' ) this%NDihedral
      call FileWriteParameter( iounit_normal, IdIdf_NDihedral )
      do i = 1, this%NDihedral
        call FileWriteBlank( iounit_normal )
        call Save( this%IdfDihedral(i) )
      end do
    end if

   ! Save information about Constraint Units
   ! Save number of constraint unites
     call FileWriteBlank( iounit_normal )
     write( IOBuffer, '(I2)' ) this%NConstraint
     call FileWriteParameter( iounit_normal, IdUnit_NConstraint )
     if( this%NConstraint > 0 ) then
       call FileWriteBlank( iounit_normal )
       do i = 1, this%NConstraint
         call FileWriteBlank( iounit_normal )
         call Save( this%Unit(i) )
       end do
     end if


    ! Close file
    call FileClose( iounit_normal )

    ! Update log file
    write( IOBuffer, '("Normalized potential model with IDF for ", A, &
&     " saved to file <", A, ">")' ) &
&     trim( this%PotModFileName ), trim( filename )
    call LogWrite

  end subroutine TMolecule_SaveIDF


!==============================================================!
!  Subroutine TMolecule_FindCOM                                !
!==============================================================!

  subroutine TMolecule_FindCOM( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TMolecule) :: this

    ! Declare local variables
    integer  :: i, j
    real(RK) :: r(3)

    ! Calculate mass of molecule and COM position
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

  end subroutine TMolecule_FindCOM

!==============================================================!
!  Subroutine TMolecule_FindMOI                                !
!==============================================================!

  subroutine TMolecule_FindMOI( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TMolecule) :: this

    ! Declare local variables
    integer  :: i
    real(RK) :: moi(3, 3), rotation(3, 3), determinant

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
    call eigen_sort( this%MOI(:), rotation(:,:) )
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
          v(:, i) = -v(:, j)
          v(:, j) = q(:)
        end if
      end do

    end subroutine eigen_sort



  end subroutine TMolecule_FindMOI



!==============================================================!
!  Subroutine TMolecule_ReadMOI                                !
!==============================================================!

  subroutine TMolecule_ReadMOI( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TMolecule) :: this

    ! Declare local variables
    integer :: i

    ! Read moments of inertia
    this%MOI(:) = 0._RK
    if( this%NDFRot > 0 ) then
      call FileReadParameter( iounit_potmod, IdSite_MOI1 )
      read( IOBuffer, * ) this%MOI(1)
      call FileReadParameter( iounit_potmod, IdSite_MOI2 )
      read( IOBuffer, * ) this%MOI(2)
      if( this%NDFRot == 3 ) then
        call FileReadParameter( iounit_potmod, IdSite_MOI3 )
        read( IOBuffer, * ) this%MOI(3)
      end if
    end if

    ! Convert to derived units
    do i = 1, 3
      this%MOI(i) = this%MOI(i) * .001_RK / NAvogadro * Angstroem**2
      this%MOI(i) = this%MOI(i) / UnitInertia
    end do



  end subroutine TMolecule_ReadMOI



!==============================================================!
!  Subroutine TMolecule_FindNDF                                !
!==============================================================!

  subroutine TMolecule_FindNDF( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TMolecule) :: this

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
&       call Error( 'Must specify moments of inertia manually' )
    end if

    ! Calculate total number of degrees of freedom
    this%NDF = 3 + this%NDFRot

    ! Set logical flags according to the number of rotation axes
    this%isElongated = this%NDFRot > 0
    this%is3D = this%NDFRot == 3

  end subroutine TMolecule_FindNDF


!==============================================================!
!  Subroutine TMolecule_FindBondR                                !
!==============================================================!

  subroutine TMolecule_FindBondR( this, Bond, j)

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TMolecule)     :: this
    type(TIdfBond)      :: Bond
    integer, intent(in) :: j

    ! Declare local variables
    integer:: SiteId1, SiteId2
    integer           :: i
    logical           :: Site1, Site2
    real(RK)          :: r1(3),r2(3)
    real(RK)          :: RX, RY, RZ
    character(10)      ::stb

    SiteId1 = Bond%SiteId1
    SiteId2 = Bond%SiteId2

    Site1 = .false.
    Site2 = .false.

     if( this%NLJ126 > 0 ) then
     do i = 1, this%NLJ126
        if (this%SiteLJ126(i)%SiteId==SiteId1) then
            r1(1)=this%SiteLJ126(i)%r(1)
            r1(2)=this%SiteLJ126(i)%r(2)
            r1(3)=this%SiteLJ126(i)%r(3)
            Site1 = .true.
        else if (this%SiteLJ126(i)%SiteId==SiteId2) then
            r2(1)=this%SiteLJ126(i)%r(1)
            r2(2)=this%SiteLJ126(i)%r(2)
            r2(3)=this%SiteLJ126(i)%r(3)
            Site2 = .true.
        end if
        if (Site1 .and. Site2) then
            exit
        end if
     end do
end if

if((.not. Site1 .or. .not. Site2) .and. (this%NCharge > 0) ) then
     do i = 1, this%NCharge
        if (this%SiteCharge(i)%SiteId==SiteId1) then
            r1(1)=this%SiteCharge(i)%r(1)
            r1(2)=this%SiteCharge(i)%r(2)
            r1(3)=this%SiteCharge(i)%r(3)
            Site1 = .true.
        else if (this%SiteCharge(i)%SiteId==SiteId2) then
            r2(1)=this%SiteCharge(i)%r(1)
            r2(2)=this%SiteCharge(i)%r(2)
            r2(3)=this%SiteCharge(i)%r(3)
            Site2 = .true.
        end if
        if (Site1 .and. Site2) then
           exit
        end if
  end do
end if


if (.not. Site1 .or. .not. Site2) then
write (stb, '(i10)') j
call Error('Uncorrect sites for bond' // stb)
end if


RX=r2(1)-r1(1)
RY=r2(2)-r1(2)
RZ=r2(3)-r1(3)
Bond%R0=dsqrt(RX**2+RY**2+RZ**2)

end subroutine TMolecule_FindBondR

!==============================================================!
!  Subroutine TMolecule_FindAngle                                !
!==============================================================!

  subroutine TMolecule_FindAngle( this, Angle, j )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TMolecule)     :: this
    type(TIdfAngle)     :: Angle
    integer, intent(in) :: j

    ! Declare local variables

    integer           :: i
    integer           :: SiteId1, SiteId2, SiteId3
    logical           :: Site1, Site2, Site3
    real(RK)          :: r1(3),r2(3),r3(3)
    real(RK)          :: R1X, R1Y, R1Z, R1S
    real(RK)          :: R2X, R2Y, R2Z, R2S
    real(RK)          ::cosa, R1R2
    character(10)     ::sta

    SiteId1 = Angle%SiteId1
    SiteId2 = Angle%SiteId2
    SiteId3 = Angle%SiteId3

    Site1 = .false.   !    (Site1) (Site3)
    Site2 = .false.   !         \  /
    Site3 = .false.   !        (Site2)


     if( this%NLJ126 > 0 ) then
     do i = 1, this%NLJ126
        if (this%SiteLJ126(i)%SiteId==SiteId1) then
            r1(1)=this%SiteLJ126(i)%r(1)
            r1(2)=this%SiteLJ126(i)%r(2)
            r1(3)=this%SiteLJ126(i)%r(3)
            Site1 = .true.
        else if (this%SiteLJ126(i)%SiteId==SiteId2) then
            r2(1)=this%SiteLJ126(i)%r(1)
            r2(2)=this%SiteLJ126(i)%r(2)
            r2(3)=this%SiteLJ126(i)%r(3)
            Site2 = .true.
        else if (this%SiteLJ126(i)%SiteId==SiteId3) then
            r3(1)=this%SiteLJ126(i)%r(1)
            r3(2)=this%SiteLJ126(i)%r(2)
            r3(3)=this%SiteLJ126(i)%r(3)
            Site3=.true.
        end if
      if (Site1 .and. Site2 .and. Site3) then
          exit
      end if
  end do
end if
if((.not. Site1 .or. .not. Site2 .or. .not. Site3) .and. (this%NCharge > 0) ) then
     do i = 1, this%NCharge
        if (this%SiteCharge(i)%SiteId==SiteId1) then
            r1(1)=this%SiteCharge(i)%r(1)
            r1(2)=this%SiteCharge(i)%r(2)
            r1(3)=this%SiteCharge(i)%r(3)
            Site1 = .true.
        else if (this%SiteCharge(i)%SiteId==SiteId2) then
            r2(1)=this%SiteCharge(i)%r(1)
            r2(2)=this%SiteCharge(i)%r(2)
            r2(3)=this%SiteCharge(i)%r(3)
            Site2 = .true.
        else if (this%SiteCharge(i)%SiteId==SiteId3) then
            r3(1)=this%SiteCharge(i)%r(1)
            r3(2)=this%SiteCharge(i)%r(2)
            r3(3)=this%SiteCharge(i)%r(3)
            Site3 = .true.
        end if
      if (Site1 .and. Site2 .and. Site3) then
          exit
      end if
  end do
end if

if (.not. Site1 .or. .not. Site2 .or. .not. Site3) then
write (sta, '(i10)') j
call Error('Uncorrect sites for angle' // sta)
end if

R1X=r1(1)-r2(1)
R1Y=r1(2)-r2(2)
R1Z=r1(3)-r2(3)
R2X=r3(1)-r2(1)
R2Y=r3(2)-r2(2)
R2Z=r3(3)-r2(3)

R1S=R1X**2+R1Y**2+R1Z**2
R2S=R2X**2+R2Y**2+R2Z**2
R1R2=dsqrt(R1S*R2S)

cosa=(R1X*R2X+R1Y*R2Y+R1Z*R2Z)/R1R2
if ( cosa .gt. 1.0d0 ) cosa = 1.0d0
if ( cosa .lt. -1.0d0) cosa = -1.0d0

Angle%Angle0=dacos(cosa)

end subroutine TMolecule_FindAngle


end module ms2_molecule

