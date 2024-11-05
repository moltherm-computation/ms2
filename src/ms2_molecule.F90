!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 4.0                !
!  (c) 2020 by TU Kaiserslautern / TU Berlin                   !
!      P.O. Box 67653                                          !
!      67653 Kaiserslautern                                    !
!==============================================================!
!  Module ms2_molecule                                         !
!  Contains TMolecule object                                   !
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
!DEC$ MESSAGE:'Compiling ms2_molecule.F90...'
#endif


module ms2_molecule

#if MPI_VER > 0
  use mpi_f08
#endif

  use ms2_global
  use ms2_site



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

    ! MIE sites
    integer :: NMIEnm
    type(TSiteMIEnm), pointer, contiguous :: SiteMIEnm(:)

    ! TT sites
    integer :: NTT
    type(TSiteTT), pointer, contiguous :: SiteTT(:)

    ! EATM sites
    integer :: NEATM
    type(TSiteEATM), pointer, contiguous :: SiteEATM(:)

    ! Coulomb sites
    integer :: NCharge
    type(TSiteCharge), pointer, contiguous :: SiteCharge(:)

    ! Dipole sites
    integer :: NDipole
    type(TSiteDipole), pointer, contiguous :: SiteDipole(:)

    ! Quadrupole sites
    integer :: NQuadrupole
    type(TSiteQuadrupole), pointer, contiguous :: SiteQuadrupole(:)

    ! File name for potential model
    character(FileNameLength) :: PotModFileName

    ! Body fixed dipole vector for reaction field
    real(RK) :: Mue(3), MueSquared

    ! Number of fluctuating states
    integer :: NFluct

    ! Total charge of the molecule
    real(RK) :: Charge

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

    ! Declare arguments
    type(TMolecule)          :: this
    character(*), intent(in) :: filename
    integer, intent(in)      :: fluctstate

    ! Declare local variables
    integer       :: i, j
    integer       :: ntypes
    character(16) :: stype
    integer       :: stat
    real(RK)      :: scalegeo, scalesig, scaleeps, scaleest

    ! Nullify pointers.
    nullify( this%SiteMIEnm )
    nullify( this%SiteTT )
    nullify( this%SiteEATM )
    nullify( this%SiteCharge )
    nullify( this%SiteDipole )
    nullify( this%SiteQuadrupole )

    ! Open potential model file
    this%PotModFileName = filename
    call FileReset(potmodFile, this%PotModFileName)

    ! Read number of potential types
    call FileReadParameter( ntypes, potmodFile%iounit, IdSite_ntypes, .false. )

    ! Zero number of sites
    this%NMIEnm = 0
    this%NTT = 0
    this%NEATM = 0
    this%NCharge = 0
    this%Charge = 0._RK
    this%NDipole = 0
    this%NQuadrupole = 0



    ! Loop over potential types
    do i = 1, ntypes
      call FileReadParameter( stype, potmodFile%iounit, IdSite_stype, .false. )
      select case( stype )
      case( 'MIEnm', 'mienm', 'MIE', 'mie', 'Mie' ) !Case: Mie-Potential
      LJorMIE = 'MIE'
        call FileReadParameter( this%NMIEnm, potmodFile%iounit, IdSite_NSites, .false. )
        if( this%NMIEnm > 0 ) then
          allocate( this%SiteMIEnm(this%NMIEnm), STAT = stat )
          call AllocationError( stat, 'MIE sites', this%NMIEnm )
          do j = 1, this%NMIEnm
            call Construct( this%SiteMIEnm(j) )
          end do
        end if

      case( 'LJ126', 'lj126', 'LJ', 'lj', 'Lj' ) !Case: LJ126-Potential
      LJorMIE = 'LJ'
        call FileReadParameter( this%NMIEnm, potmodFile%iounit, IdSite_NSites, .false. )
        if( this%NMIEnm > 0 ) then
          allocate( this%SiteMIEnm(this%NMIEnm), STAT = stat )
          call AllocationError( stat, 'LJ sites', this%NMIEnm )
          do j = 1, this%NMIEnm
            call Construct( this%SiteMIEnm(j) )
          end do
        end if

      case( 'TT68', 'tt68', 'tt' ) !Case: Tang-Tönnies-Potential
        TT68orEXT = 'TT68'
          call FileReadParameter( this%NTT, potmodFile%iounit, IdSite_NSites, .false. )
          if( this%NTT > 0 ) then
            allocate( this%SiteTT(this%NTT), STAT = stat )
            call AllocationError( stat, 'TT sites', this%NTT )
            do j = 1, this%NTT
              call Construct( this%SiteTT(j) )
            end do
          end if
  
      case( 'TTExt', 'TText', 'ttext' ) !Case: Jäger-Potential
        TT68orEXT = 'TTExt'
          call FileReadParameter( this%NTT, potmodFile%iounit, IdSite_NSites, .false. )
          if( this%NTT > 0 ) then
            allocate( this%SiteTT(this%NTT), STAT = stat )
            call AllocationError( stat, 'TT sites', this%NTT )
            do j = 1, this%NTT
              call Construct( this%SiteTT(j) )
            end do
          end if

        case( 'EATM', 'eATM', 'eatm' ) !Case: extended Axilrod-Teller-Muto potential
            call FileReadParameter( this%NEATM, potmodFile%iounit, IdSite_NSites, .false. )
            if( this%NEATM > 0 ) then
              allocate( this%SiteEATM(this%NEATM), STAT = stat )
              call AllocationError( stat, 'EATM sites', this%NEATM )
              do j = 1, this%NEATM
                call Construct( this%SiteEATM(j) )
              end do
            end if

      case( 'CHARGE', 'Charge', 'charge', 'E', 'e' )
        call FileReadParameter( this%NCharge, potmodFile%iounit, IdSite_NSites, .false. )
        if( this%NCharge > 0 ) then
          allocate( this%SiteCharge(this%NCharge), STAT = stat )
          call AllocationError( stat, 'point charge sites', this%NCharge )
          do j = 1, this%NCharge
            call Construct( this%SiteCharge(j) )
            this%Charge = this%Charge + this%SiteCharge(j)%e
          end do
        end if

      case( 'DIPOLE', 'Dipole', 'dipole', 'D', 'd' )
        call FileReadParameter( this%NDipole, potmodFile%iounit, IdSite_NSites, .false. )
        if( this%NDipole > 0 ) then
          allocate( this%SiteDipole(this%NDipole), STAT = stat )
          call AllocationError( stat, 'dipolar sites', this%NDipole )
          do j = 1, this%NDipole
            call Construct( this%SiteDipole(j) )
          end do
        end if

      case( 'QUADRUPOLE', 'Quadrupole', 'quadrupole', 'Q', 'q' )
        call FileReadParameter( this%NQuadrupole, potmodFile%iounit, IdSite_NSites, .false. )
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
    call FileReadParameter( stype, potmodFile%iounit, IdSite_NDFRot, .false. )
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
      call FileReadParameter_IOBuffer( potmodFile%iounit, IdNFluct, .false. )

      ! Scaling factors start in next line
      if( RootProc ) then
        do i = 1, fluctstate
          read( potmodFile%iounit, * ) scalegeo, scalesig, scaleeps, scaleest
        end do
      end if
#if MPI_VER > 0
      call MPI_Bcast( scalegeo, 1, MPI_RK, NRootProc, Communicator, ierror )
      call MPI_Bcast( scalesig, 1, MPI_RK, NRootProc, Communicator, ierror )
      call MPI_Bcast( scaleeps, 1, MPI_RK, NRootProc, Communicator, ierror )
      call MPI_Bcast( scaleest, 1, MPI_RK, NRootProc, Communicator, ierror )
#endif
      if( scalegeo > 1._RK .or. scalesig > 1._RK .or. scaleeps > 1._RK .or. scaleest > 1._RK ) &
&         call Error( 'Scaling factors for fluctuating particle must be lower or equal 1' )

      ! Apply scaling factors
      do i = 1, this%NMIEnm
        this%SiteMIEnm(i)%r = this%SiteMIEnm(i)%r * scalegeo
        this%SiteMIEnm(i)%sig = this%SiteMIEnm(i)%sig * scalesig
        this%SiteMIEnm(i)%eps = this%SiteMIEnm(i)%eps * scaleeps
      end do

      do i = 1, this%NTT
        this%SiteTT(i)%r = this%SiteTT(i)%r * scalegeo
        this%SiteTT(i)%shield = this%SiteTT(i)%shield * scalegeo
      end do

      do i = 1, this%NEATM
        this%SiteEATM(i)%r = this%SiteEATM(i)%r * scalegeo
        this%SiteEATM(i)%shield = this%SiteEATM(i)%shield * scalegeo
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

      call FileReadParameter( this%NFluct, potmodFile%iounit, IdNFluct, .false. )

    else

      this%NFluct = 0

    end if

    ! Close potential model file
    call FileClose(potmodFile)

    ! Reduction of point charges and dipoles to body fixed dipole vector
    this%Mue(:) = 0._RK
    if( (this%NCharge > 0).or.(this%NDipole > 0) ) then
      if ((LongRange .ne. Ewald).or.(ODFUpdateFrequency > 0)) then
        if ((LongRange .ne. PME).or.(ODFUpdateFrequency > 0)) then
          do i =1, this%NCharge
            this%Mue(:) = this%Mue(:) + this%SiteCharge(i)%r(:) * this%SiteCharge(i)%e
          end do
        end if
      end if
      do i =1, this%NDipole
        this%Mue(:) = this%Mue(:) + this%SiteDipole(i)%or(:) * this%SiteDipole(i)%D
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
    if( associated( this%SiteMIEnm ) ) then
      do i = 1, this%NMIEnm
        call Destruct( this%SiteMIEnm(i) )
      end do
      deallocate( this%SiteMIEnm )
    end if
    if( associated( this%SiteTT ) ) then
      do i = 1, this%NTT
        call Destruct( this%SiteTT(i) )
      end do
      deallocate( this%SiteTT )
    end if
    if( associated( this%SiteEATM ) ) then
      do i = 1, this%NEATM
        call Destruct( this%SiteEATM(i) )
      end do
      deallocate( this%SiteEATM )
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

    i = index(this%PotModFileName,'.',.true.)
    if ( i<=1 ) then
      i = len(trim(this%PotModFileName))
    end if

    if( fluctstate < 1 ) then
      write( filename, '(A,".",A,A)') trim(OutputNameTag),trim( this%PotModFileName(1:i-1) ),trim(NormalizedPotModExtension)

    else
      write( filename, '(A,".",A,"_",I0,A)') trim(OutputNameTag),trim( this%PotModFileName(1:i-1) ),fluctstate &
&           ,trim(NormalizedPotModExtension)
    end if
    call FileRewrite(normalFile, filename)

    ! Save number of potential types
    ntypes = 0
    if( this%NMIEnm > 0 ) ntypes = ntypes + 1
    if( this%NTT > 0 ) ntypes = ntypes + 1
    if( this%NEATM > 0 ) ntypes = ntypes + 1
    if( this%NCharge > 0 ) ntypes = ntypes + 1
    if( this%NDipole > 0 ) ntypes = ntypes + 1
    if( this%NQuadrupole > 0 ) ntypes = ntypes + 1
    write( IOBuffer, '(I2)' ) ntypes
    call FileWriteParameter( normalFile%iounit, IdSite_ntypes )

    ! Save MIE sites
    if( this%NMIEnm > 0 ) then
      call FileWriteBlank(normalFile)
      write( IOBuffer, '(1X, A)' ) LJorMIE !'MIEnm'
      call FileWriteParameter( normalFile%iounit, IdSite_stype )
      write( IOBuffer, '(I2)' ) this%NMIEnm
      call FileWriteParameter( normalFile%iounit, IdSite_NSites )
      do i = 1, this%NMIEnm
        call FileWriteBlank(normalFile)
        call Save( this%SiteMIEnm(i) )
      end do
    end if

    ! Save TT sites
    if( this%NTT > 0 ) then
      call FileWriteBlank( normalFile )
      write( IOBuffer, '(1X, A)' ) TT68orEXT
      call FileWriteParameter( normalFile%iounit, IdSite_stype )
      write( IOBuffer, '(I2)' ) this%NTT
      call FileWriteParameter( normalFile%iounit, IdSite_NSites )
      do i = 1, this%NTT
        call FileWriteBlank(normalFile)
        call Save( this%SiteTT(i) )
      end do
    end if

    ! Save EATM sites
    if( this%NEATM > 0 ) then
      call FileWriteBlank( normalFile )
      write( IOBuffer, '(1X, A)' ) 'EATM'
      call FileWriteParameter( normalFile%iounit, IdSite_stype )
      write( IOBuffer, '(I2)' ) this%NEATM
      call FileWriteParameter( normalFile%iounit, IdSite_NSites )
      do i = 1, this%NEATM
        call FileWriteBlank(normalFile)
        call Save( this%SiteEATM(i) )
      end do
    end if

    ! Save point charge sites
    if( this%NCharge > 0 ) then
      call FileWriteBlank(normalFile)
      write( IOBuffer, '(1X, A)' ) 'Charge'
      call FileWriteParameter( normalFile%iounit, IdSite_stype )
      write( IOBuffer, '(I2)' ) this%NCharge
      call FileWriteParameter( normalFile%iounit, IdSite_NSites )
      do i = 1, this%NCharge
        call FileWriteBlank(normalFile)
        call Save( this%SiteCharge(i) )
      end do
    end if

    ! Save point dipole sites
    if( this%NDipole > 0 ) then
      call FileWriteBlank(normalFile)
      write( IOBuffer, '(1X, A)' ) 'Dipole'
      call FileWriteParameter( normalFile%iounit, IdSite_stype )
      write( IOBuffer, '(I2)' ) this%NDipole
      call FileWriteParameter( normalFile%iounit, IdSite_NSites )
      do i = 1, this%NDipole
        call FileWriteBlank(normalFile)
        call Save( this%SiteDipole(i) )
      end do
    end if

    ! Save point quadrupole sites
    if( this%NQuadrupole > 0 ) then
      call FileWriteBlank(normalFile)
      write( IOBuffer, '(1X, A)' ) 'Quadrupole'
      call FileWriteParameter( normalFile%iounit, IdSite_stype )
      write( IOBuffer, '(I2)' ) this%NQuadrupole
      call FileWriteParameter( normalFile%iounit, IdSite_NSites )
      do i = 1, this%NQuadrupole
        call FileWriteBlank(normalFile)
        call Save( this%SiteQuadrupole(i) )
      end do
    end if

    ! Save number of rotation axes
    call FileWriteBlank(normalFile)
    write( IOBuffer, '(I2)' ) this%NDFRot
    call FileWriteParameter( normalFile%iounit, IdSite_NDFRot )

    ! Save total mass of the molecule
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&          this%Mass * UnitMass * 1000._RK * NAvogadro, this%Mass
    call FileWriteParameter( normalFile%iounit, IdMolecule_mass )

    ! Save moments of inertia
    if( this%NDFRot > 0 ) then
      write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&            this%MOI(1) * UnitInertia * 1000._RK * NAvogadro / Angstroem**2, &
&            this%MOI(1)

      call FileWriteParameter( normalFile%iounit, IdSite_MOI1 )
      write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&            this%MOI(2) * UnitInertia * 1000._RK * NAvogadro / Angstroem**2, &
&            this%MOI(2)

      call FileWriteParameter( normalFile%iounit, IdSite_MOI2 )
      write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&            this%MOI(3) * UnitInertia * 1000._RK * NAvogadro / Angstroem**2, &
&            this%MOI(3)
      call FileWriteParameter( normalFile%iounit, IdSite_MOI3 )
    end if

    ! Close file
    call FileClose(normalFile)

    ! Update log file
    write( IOBuffer, '("Normalized potential model for ", A, &
&          " saved to file <", A, ">")' ) trim( this%PotModFileName ), trim( filename )
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
    do i = 1, this%NMIEnm
      this%Mass = this%Mass + this%SiteMIEnm(i)%mass
      r(:) = r(:) + this%SiteMIEnm(i)%mass * this%SiteMIEnm(i)%r(:)
    end do
    do i = 1, this%NTT
      this%Mass = this%Mass + this%SiteTT(i)%mass
      r(:) = r(:) + this%SiteTT(i)%mass * this%SiteTT(i)%r(:)
    end do
    do i = 1, this%NEATM
      this%Mass = this%Mass + this%SiteEATM(i)%mass
      r(:) = r(:) + this%SiteEATM(i)%mass * this%SiteEATM(i)%r(:)
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
    do i = 1, this%NMIEnm
      do j = 1, 3
        this%SiteMIEnm(i)%r(j) = this%SiteMIEnm(i)%r(j) - r(j)
      end do
    end do
    do i = 1, this%NTT
      do j = 1, 3
        this%SiteTT(i)%r(j) = this%SiteTT(i)%r(j) - r(j)
      end do
    end do
    do i = 1, this%NEATM
      do j = 1, 3
        this%SiteEATM(i)%r(j) = this%SiteEATM(i)%r(j) - r(j)
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
    do i = 1, this%NMIEnm
      moi(1, 1) = moi(1, 1) + this%SiteMIEnm(i)%mass * ( this%SiteMIEnm(i)%r(2)**2 + this%SiteMIEnm(i)%r(3)**2 )
      moi(1, 2) = moi(1, 2) - this%SiteMIEnm(i)%mass * this%SiteMIEnm(i)%r(1) * this%SiteMIEnm(i)%r(2)
      moi(1, 3) = moi(1, 3) - this%SiteMIEnm(i)%mass * this%SiteMIEnm(i)%r(1) * this%SiteMIEnm(i)%r(3)
      moi(2, 2) = moi(2, 2) + this%SiteMIEnm(i)%mass * ( this%SiteMIEnm(i)%r(1)**2 + this%SiteMIEnm(i)%r(3)**2 )
      moi(2, 3) = moi(2, 3) - this%SiteMIEnm(i)%mass * this%SiteMIEnm(i)%r(2) * this%SiteMIEnm(i)%r(3)
      moi(3, 3) = moi(3, 3) + this%SiteMIEnm(i)%mass * ( this%SiteMIEnm(i)%r(1)**2 + this%SiteMIEnm(i)%r(2)**2 )
    end do

    do i = 1, this%NTT
      moi(1, 1) = moi(1, 1) + this%SiteTT(i)%mass * ( this%SiteTT(i)%r(2)**2 + this%SiteTT(i)%r(3)**2 )
      moi(1, 2) = moi(1, 2) - this%SiteTT(i)%mass * this%SiteTT(i)%r(1) * this%SiteTT(i)%r(2)
      moi(1, 3) = moi(1, 3) - this%SiteTT(i)%mass * this%SiteTT(i)%r(1) * this%SiteTT(i)%r(3)
      moi(2, 2) = moi(2, 2) + this%SiteTT(i)%mass * ( this%SiteTT(i)%r(1)**2 + this%SiteTT(i)%r(3)**2 )
      moi(2, 3) = moi(2, 3) - this%SiteTT(i)%mass * this%SiteTT(i)%r(2) * this%SiteTT(i)%r(3)
      moi(3, 3) = moi(3, 3) + this%SiteTT(i)%mass * ( this%SiteTT(i)%r(1)**2 + this%SiteTT(i)%r(2)**2 )
    end do

    do i = 1, this%NEATM
      moi(1, 1) = moi(1, 1) + this%SiteEATM(i)%mass * ( this%SiteEATM(i)%r(2)**2 + this%SiteEATM(i)%r(3)**2 )
      moi(1, 2) = moi(1, 2) - this%SiteEATM(i)%mass * this%SiteEATM(i)%r(1) * this%SiteEATM(i)%r(2)
      moi(1, 3) = moi(1, 3) - this%SiteEATM(i)%mass * this%SiteEATM(i)%r(1) * this%SiteEATM(i)%r(3)
      moi(2, 2) = moi(2, 2) + this%SiteEATM(i)%mass * ( this%SiteEATM(i)%r(1)**2 + this%SiteEATM(i)%r(3)**2 )
      moi(2, 3) = moi(2, 3) - this%SiteEATM(i)%mass * this%SiteEATM(i)%r(2) * this%SiteEATM(i)%r(3)
      moi(3, 3) = moi(3, 3) + this%SiteEATM(i)%mass * ( this%SiteEATM(i)%r(1)**2 + this%SiteEATM(i)%r(2)**2 )
    end do

    do i = 1, this%NCharge
      moi(1, 1) = moi(1, 1) + this%SiteCharge(i)%mass * ( this%SiteCharge(i)%r(2)**2 + this%SiteCharge(i)%r(3)**2 )
      moi(1, 2) = moi(1, 2) - this%SiteCharge(i)%mass * this%SiteCharge(i)%r(1) * this%SiteCharge(i)%r(2)
      moi(1, 3) = moi(1, 3) - this%SiteCharge(i)%mass * this%SiteCharge(i)%r(1) * this%SiteCharge(i)%r(3)
      moi(2, 2) = moi(2, 2) + this%SiteCharge(i)%mass * ( this%SiteCharge(i)%r(1)**2 + this%SiteCharge(i)%r(3)**2 )
      moi(2, 3) = moi(2, 3) - this%SiteCharge(i)%mass * this%SiteCharge(i)%r(2) * this%SiteCharge(i)%r(3)
      moi(3, 3) = moi(3, 3) + this%SiteCharge(i)%mass * ( this%SiteCharge(i)%r(1)**2 + this%SiteCharge(i)%r(2)**2 )
    end do

    do i = 1, this%NDipole
      moi(1, 1) = moi(1, 1) + this%SiteDipole(i)%mass * ( this%SiteDipole(i)%r(2)**2 + this%SiteDipole(i)%r(3)**2 )
      moi(1, 2) = moi(1, 2) - this%SiteDipole(i)%mass * this%SiteDipole(i)%r(1) * this%SiteDipole(i)%r(2)
      moi(1, 3) = moi(1, 3) - this%SiteDipole(i)%mass * this%SiteDipole(i)%r(1) * this%SiteDipole(i)%r(3)
      moi(2, 2) = moi(2, 2) + this%SiteDipole(i)%mass * ( this%SiteDipole(i)%r(1)**2 + this%SiteDipole(i)%r(3)**2 )
      moi(2, 3) = moi(2, 3) - this%SiteDipole(i)%mass * this%SiteDipole(i)%r(2) * this%SiteDipole(i)%r(3)
      moi(3, 3) = moi(3, 3) + this%SiteDipole(i)%mass * ( this%SiteDipole(i)%r(1)**2 + this%SiteDipole(i)%r(2)**2 )
    end do

    do i = 1, this%NQuadrupole
      moi(1, 1) = moi(1, 1) + this%SiteQuadrupole(i)%mass * ( this%SiteQuadrupole(i)%r(2)**2 + this%SiteQuadrupole(i)%r(3)**2 )
      moi(1, 2) = moi(1, 2) - this%SiteQuadrupole(i)%mass * this%SiteQuadrupole(i)%r(1) * this%SiteQuadrupole(i)%r(2)
      moi(1, 3) = moi(1, 3) - this%SiteQuadrupole(i)%mass * this%SiteQuadrupole(i)%r(1) * this%SiteQuadrupole(i)%r(3)
      moi(2, 2) = moi(2, 2) + this%SiteQuadrupole(i)%mass * ( this%SiteQuadrupole(i)%r(1)**2 + this%SiteQuadrupole(i)%r(3)**2 )
      moi(2, 3) = moi(2, 3) - this%SiteQuadrupole(i)%mass * this%SiteQuadrupole(i)%r(2) * this%SiteQuadrupole(i)%r(3)
      moi(3, 3) = moi(3, 3) + this%SiteQuadrupole(i)%mass * ( this%SiteQuadrupole(i)%r(1)**2 + this%SiteQuadrupole(i)%r(2)**2 )
    end do

    ! Transform to principal axes
    call eigen_find( moi(:,:), this%MOI(:), rotation(:,:) )
    call eigen_sort( this%MOI(:), rotation(:,:) )
    do i = 1, this%NMIEnm
      this%SiteMIEnm(i)%r(:) = matmul( this%SiteMIEnm(i)%r(:), rotation(:, :) )
    end do

    do i = 1, this%NTT
      this%SiteTT(i)%r(:) = matmul( this%SiteTT(i)%r(:), rotation(:, :) )
    end do

    do i = 1, this%NEATM
      this%SiteEATM(i)%r(:) = matmul( this%SiteEATM(i)%r(:), rotation(:, :) )
    end do

    do i = 1, this%NCharge
      this%SiteCharge(i)%r(:) = matmul( this%SiteCharge(i)%r(:), rotation(:, :) )
    end do

    do i = 1, this%NDipole
      this%SiteDipole(i)%r(:) = matmul( this%SiteDipole(i)%r(:), rotation(:, :) )
      this%SiteDipole(i)%or(:) = matmul( this%SiteDipole(i)%or(:), rotation(:, :) )
    end do

    do i = 1, this%NQuadrupole
      this%SiteQuadrupole(i)%r(:) = matmul( this%SiteQuadrupole(i)%r(:), rotation(:, :) )
      this%SiteQuadrupole(i)%or(:) = matmul( this%SiteQuadrupole(i)%or(:), rotation(:, :) )
    end do

    if( (this%NCharge > 0).or.(this%NDipole > 0) ) this%Mue(:) = matmul( this%Mue(:), rotation(:, :) )

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
            if((i > 4) .and. (abs( d(ip) ) + g == abs( d(ip) )) .and. (abs( d(iq) ) + g == abs( d(iq) ))) then
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
      call FileReadParameter( this%MOI(1), potmodFile%iounit, IdSite_MOI1, .false. )
      call FileReadParameter( this%MOI(2), potmodFile%iounit, IdSite_MOI2, .false. )
      if( this%NDFRot == 3 ) then
        call FileReadParameter( this%MOI(3), potmodFile%iounit, IdSite_MOI3, .false. )
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
      disoriented = this%NDFRot < 2 .and. (this%NDipole > 0 .or. this%NQuadrupole > 0)

      do i = 1, this%NDipole
        disoriented = disoriented .or. ( maxval( abs( this%SiteDipole(i)%or(1:2) ) ) > Zero )
      end do

      do i = 1, this%NQuadrupole
        disoriented = disoriented .or. ( maxval( abs( this%SiteQuadrupole(i)%or(1:2) ) ) > Zero )
      end do

      if( disoriented ) call Error( 'Must specify moments of inertia manually' )
    end if

    ! Calculate total number of degrees of freedom
    this%NDF = 3 + this%NDFRot

    ! Set logical flags according to the number of rotation axes
    this%isElongated = this%NDFRot > 0
    this%is3D = this%NDFRot == 3

  end subroutine TMolecule_FindNDF



end module ms2_molecule
