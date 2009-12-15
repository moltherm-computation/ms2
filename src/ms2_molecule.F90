!==============================================================!
!  MOLECULAR SIMULATION PROGRAM MS2 Version 1.1 v12            !
!  (c) 2001 by Sergey Lishchuk, ITT                            !
!  (c) 2007 by Bernhard Eckl, ITT                              !
!==============================================================!
!  Module ms2_molecule                                         !
!  Contains TMolecule object                                   !
!==============================================================!

#ifndef ARCH
#define ARCH    0
#define FORTRAN 90
#define MPI_VER 0
#define FVM_VER 0
#endif

#ifndef TRANS
#define TRANS 0
#endif

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2_molecule.F90...'
#endif

module ms2_molecule

!#ifdev MPI_VER > 0
!  use mpi
!#endif

  use ms2_global
  use ms2_site

#if defined PAR_PROF
  use ms2_profiler
#endif

#if FVM_VER > 0
  use libfvmf2003
  use fvmf2003extensions
  use, intrinsic :: iso_c_binding
#endif


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

    ! File name for potential model
    character(FileNameLength) :: PotModFileName

    ! Body fixed dipole vector for reaction field
    real(RK) :: Mue(3), MueSquared

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
    integer       :: i, j
    integer       :: ntypes
    character(16) :: stype
    integer       :: stat
#if FVM_VER == 0
    real(RK)      :: scalegeo, scalesig, scaleeps, scaleest
#endif

    ! Nullify pointers.
    nullify( this%SiteLJ126 )
    nullify( this%SiteCharge )
    nullify( this%SiteDipole )
    nullify( this%SiteQuadrupole )

    ! Open potential model file
    this%PotModFileName = filename
    call FileReset( iounit_potmod, this%PotModFileName )

    ! Read number of potential types
    call FileReadParameter( ntypes, iounit_potmod, IdSite_ntypes, .false. )

    ! Zero number of sites
    this%NLJ126 = 0
    this%NCharge = 0
    this%NDipole = 0
    this%NQuadrupole = 0

    ! Loop over potential types
    do i = 1, ntypes
      call FileReadParameter( stype, iounit_potmod, IdSite_stype, .false. )
      select case( stype )
      case( 'LJ126', 'lj126', 'LJ', 'lj' )
        call FileReadParameter( this%NLJ126, iounit_potmod, IdSite_NLJ126, .false. )
        if( this%NLJ126 > 0 ) then
          allocate( this%SiteLJ126(this%NLJ126), STAT = stat )
          call AllocationError( stat, 'Lennard-Jones sites', this%NLJ126 )
          do j = 1, this%NLJ126
            call Construct( this%SiteLJ126(j) )
          end do
        end if
      case( 'CHARGE', 'Charge', 'charge', 'E', 'e' )
        call FileReadParameter( this%NCharge, iounit_potmod, IdSite_NCharge, .false. )
        if( this%NCharge > 0 ) then
          allocate( this%SiteCharge(this%NCharge), STAT = stat )
          call AllocationError( stat, 'point charge sites', this%NCharge )
          do j = 1, this%NCharge
            call Construct( this%SiteCharge(j) )
          end do
        end if
      case( 'DIPOLE', 'Dipole', 'dipole', 'D', 'd' )
        call FileReadParameter( this%NDipole, iounit_potmod, IdSite_NDipole, .false. )
        if( this%NDipole > 0 ) then
          allocate( this%SiteDipole(this%NDipole), STAT = stat )
          call AllocationError( stat, 'dipolar sites', this%NDipole )
          do j = 1, this%NDipole
            call Construct( this%SiteDipole(j) )
          end do
        end if
      case( 'QUADRUPOLE', 'Quadrupole', 'quadrupole', 'Q', 'q' )
        call FileReadParameter( this%NQuadrupole, iounit_potmod, IdSite_NQuadrupole, .false. )
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
    call FileReadParameter( stype, iounit_potmod, IdSite_NDFRot, .false. )
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

    ! Find moments of inertia
    if( this%NDFRot < 0 ) then
      call FindMOI( this )
    else
      call ReadMOI( this )
    end if

    ! Find number of degrees of freedom
    call FindNDF( this )

    ! For fluctuating particle scale parameters
    if( fluctstate > 0 ) then
      call FileReadParameter_IOBuffer( iounit_potmod, IdNFluct, .false. )

      ! Scaling factors start in next line
      if( RootProc ) then
        do i = 1, fluctstate
          read( iounit_potmod, * ) scalegeo, scalesig, scaleeps, scaleest
        end do
      end if

#if defined PAR_PROF

      ! Parallel Profiling added by Hendrik Adorf (ITWM)
      call profileTagBefore( Profiler, &
&       'TMolecule_Construct: Bcast(several scales)' )

#endif

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
#if FVM_VER > 0

      fvmret = pv4dBarrier()

      !FVM_Bcast
      fvmret = readdma(fvmByteOffScalegeo, fvmByteOffScalegeo, &
&       sizeof(scalegeo), NRootProc, 0)
      fvmret = readdma(fvmByteOffScalesig, fvmByteOffScalesig, &
&       sizeof(scalesig), NRootProc, 1)
      fvmret = readdma(fvmByteOffScaleeps, fvmByteOffScaleeps, &
&       sizeof(scaleeps), NRootProc, 2)
      fvmret = readdma(fvmByteOffScaleest, fvmByteOffScaleest, &
&       sizeof(scaleest), NRootProc, 3)
      fvmret = waitonqueue(0)
      fvmret = waitonqueue(1)
      fvmret = waitonqueue(2)
      fvmret = waitonqueue(3)

      fvmret = pv4dBarrier()

#endif

#if defined PAR_PROF

      ! Parallel Profiling added by Hendrik Adorf (ITWM)
      call profileTagAfter( Profiler, &
&       'TMolecule_Construct: Bcast(several scales)' )

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

      call FileReadParameter( this%NFluct, iounit_potmod, IdNFluct, .false. )

    else

      this%NFluct = 0

    end if

    ! Close potential model file
    call FileClose( iounit_potmod )

    ! Reduction of point charges and dipoles to body fixed dipole vector
    this%Mue(:) = 0._RK
    if( (this%NCharge > 0).or.(this%NDipole > 0) ) then
      do i =1, this%NCharge
        this%Mue(:) = this%Mue(:) + &
&         this%SiteCharge(i)%r(:) * this%SiteCharge(i)%e
      end do
      do i =1, this%NDipole
        this%Mue(:) = this%Mue(:) + &
&         this%SiteDipole(i)%or(:) * this%SiteDipole(i)%D
      end do
    end if
    this%MueSquared = sum( this%Mue(:)**2 )

    ! Save used potential model
    call Save( this, fluctstate )

  end subroutine TMolecule_Construct



!==============================================================!
!  Subroutine TMolecule_Destruct                               !
!==============================================================!

  subroutine TMolecule_Destruct( this )

    implicit none

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

    ! Save number of rotation axes
    call FileWriteBlank( iounit_normal )
    write( IOBuffer, '(I2)' ) this%NDFRot
    call FileWriteParameter( iounit_normal, IdSite_NDFRot )

    ! Save total mass of the molecule
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%Mass * UnitMass * 1000._RK * NAvogadro, this%Mass
    call FileWriteParameter( iounit_normal, IdSite_Mass )

    ! Save moments of inertia
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
!  Subroutine TMolecule_FindCOM                                !
!==============================================================!

  subroutine TMolecule_FindCOM( this )

    implicit none

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

    ! Declare arguments
    type(TMolecule) :: this

    ! Declare local variables
    integer  :: i
    real(RK) :: moi(3, 3), rotation(3, 3)

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
    if( (this%NCharge > 0).or.(this%NDipole > 0) ) &
&     this%Mue(:) = matmul( this%Mue(:), rotation(:, :) )



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



  end subroutine TMolecule_FindMOI



!==============================================================!
!  Subroutine TMolecule_ReadMOI                                !
!==============================================================!

  subroutine TMolecule_ReadMOI( this )

    implicit none

    ! Declare arguments
    type(TMolecule) :: this

    ! Declare local variables
    integer :: i

    ! Read moments of inertia
    this%MOI(:) = 0._RK
    if( this%NDFRot > 0 ) then
      call FileReadParameter( this%MOI(1), iounit_potmod, IdSite_MOI1, .false. )
      call FileReadParameter( this%MOI(2), iounit_potmod, IdSite_MOI2, .false. )
      if( this%NDFRot == 3 ) then
        call FileReadParameter( this%MOI(3), iounit_potmod, IdSite_MOI3, .false. )
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



end module ms2_molecule

