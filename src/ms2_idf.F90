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
!DEC$ MESSAGE:'Compiling ms2_site.F90...'
#endif

module ms2_idf

  use ms2_global


  type coordinatesPointer

    real(RK), pointer, contiguous :: X(:), Y(:), Z(:)

  end type coordinatesPointer


!==============================================================!
!  Type TIdfBond                                             !
!==============================================================!

  type TIdfBond

    real(RK)          :: ForConst
    real(RK)          :: R0
    integer, dimension(2) :: SiteId, UnitId
    integer, pointer  :: NPartMax, NPart
    integer, pointer  :: NPart0, NPart1, NPart2
    type(coordinatesPointer), dimension(3) :: P, R, F

  end type TIdfBond

  interface Construct
    module procedure TIdfBond_Construct
  end interface

  interface Destruct
    module procedure TIdfBond_Destruct
  end interface

  interface Save
    module procedure TIdfBond_Save
  end interface


!==============================================================!
!  Type TIdfAngle                                              !
!==============================================================!

  type TIdfAngle

    real(RK)          :: ForConst
    real(RK)          :: Angle0, Angle
    integer, dimension(3) :: SiteId, UnitId
    integer, pointer  :: NPartMax, NPart
    integer, pointer  :: NPart0, NPart1, NPart2
    type(coordinatesPointer), dimension(3) :: R, F
    logical           :: orientation1, orientation2


  end type TIdfAngle

  interface Construct
    module procedure TIdfAngle_Construct
  end interface

  interface Destruct
    module procedure TIdfAngle_Destruct
  end interface

  interface Save
    module procedure TIdfAngle_Save
  end interface



!==============================================================!
!  Type TIdfDihedral                                              !
!==============================================================!

  type TIdfDihedral

    real(RK)          :: phi
    integer           :: nmax    ! multiplicity
    integer, dimension(4) :: SiteId, UnitId
    real(RK)          :: ScaleLJ14
    real(RK)          :: ScaleEl14
    integer, pointer  :: NPartMax, NPart
    integer, pointer  :: NPart0, NPart1, NPart2
    real(RK), pointer, contiguous :: ForConst(:), gamma0(:)
    type(coordinatesPointer), dimension(4) :: R, F
    logical           :: orientation1, orientation2

  end type TIdfDihedral

  interface Construct
    module procedure TIdfDihedral_Construct
  end interface

  interface Destruct
    module procedure TIdfDihedral_Destruct
  end interface

  interface Save
    module procedure TIdfDihedral_Save
  end interface


contains


!==============================================================!
!  Subroutine TIdfBond_Construct                               !
!==============================================================!

  subroutine TIdfBond_Construct( this )

    implicit none

    ! Declare arguments
    type(TIdfBond) :: this

    ! Read site parameters
    call FileReadParameter_IOBuffer( potmodFile%iounit, IdBond_Sites )
    read( IOBuffer, * ) this%SiteId(1), this%SiteId(2)
    call FileReadParameter( this%R0, potmodFile%iounit, IdBond_R0, .false.)
    if (Shake > 0) then
      this%ForConst = 1e07_RK
    else
      call FileReadParameter( this%ForConst, potmodFile%iounit, IdBond_ForConst, .false.)
      if (this%ForConst < 100) then
        call LogWriteBlank
        write( IOBuffer, '("WARNING: Check your bond definition. Currently a/some connections are faulty.")' )
        call LogWrite
        write( IOBuffer, '("Either use Shake to keep bonds rigd or choose a reasonably high force constant.")' )
        call LogWrite
      end if
    end if

    ! Convert to SI units
    this%R0 = this%R0 * Angstroem
    this%ForConst = this%ForConst * kBoltzmann

    ! Convert to derived units
    this%R0 = this%R0 / UnitLength
    this%ForConst = this%ForConst / UnitEnergy


end subroutine TIdfBond_Construct


!==============================================================!
!  Subroutine TIdfBond_Destruct                              !
!==============================================================!

  subroutine TIdfBond_Destruct( this )

    implicit none

    ! Declare arguments
    type(TIdfBond) :: this

    ! Destroy site
    continue

  end subroutine TIdfBond_Destruct

!==============================================================!
!  Subroutine TIdfBond_Save                                    !
!==============================================================!

  subroutine TIdfBond_Save( this )

    implicit none

    ! Declare arguments
    type(TIdfBond) :: this

    ! Save site parameters

    write( IOBuffer, '(2I3)' ) this%SiteId(1), this%SiteId(2)
    call FileWriteParameter( normalFile%iounit, IdBond_Sites )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%R0 * UnitLength / Angstroem, this%R0
    call FileWriteParameter( normalFile%iounit, IdBond_R0 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%ForConst * UnitEnergy * Angstroem * Angstroem / kBoltzmann / UnitLength / UnitLength, this%ForConst
    call FileWriteParameter( normalFile%iounit, IdBond_ForConst )

  end subroutine TIdfBond_Save


!==============================================================!
!  Subroutine TIdfAngle_Construct                             !
!==============================================================!

  subroutine TIdfAngle_Construct( this )

    implicit none

    ! Declare arguments
    type(TIdfAngle) :: this

    ! Read site parameters
    call FileReadParameter_IOBuffer( potmodFile%iounit, IdAngle_Sites )
    read( IOBuffer, * ) this%SiteId(1), this%SiteId(2), this%SiteId(3)
    call FileReadParameter( this%Angle0, potmodFile%iounit, IdAngle_Angle0, .false.)
    call FileReadParameter( this%ForConst, potmodFile%iounit, IdAngle_ForConst, .false.)

    ! Convert to SI units
    this%ForConst = this%ForConst * kBoltzmann

    ! Convert to derived units
    this%Angle0 = this%Angle0*Pi/180
    this%ForConst = this%ForConst / UnitEnergy

end subroutine TIdfAngle_Construct


!==============================================================!
!  Subroutine TIdfAngle_Destruct                              !
!==============================================================!

  subroutine TIdfAngle_Destruct( this )

    implicit none

    ! Declare arguments
    type(TIdfAngle) :: this

    ! Destroy site
    continue

  end subroutine TIdfAngle_Destruct


!==============================================================!
!  Subroutine TIdfAngle_Save                                    !
!==============================================================!

  subroutine TIdfAngle_Save( this )

    implicit none

    ! Declare arguments
    type(TIdfAngle) :: this

    ! Save site parameters

    write( IOBuffer, '(3I3)' ) this%SiteId(1), this%SiteId(2), this%SiteId(3)
    call FileWriteParameter( normalFile%iounit, IdAngle_Sites )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%Angle0*180/PI, this%Angle0
    call FileWriteParameter( normalFile%iounit, IdAngle_Angle0 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%ForConst * UnitEnergy / kBoltzmann, this%ForConst
    call FileWriteParameter( normalFile%iounit, IdAngle_ForConst )

  end subroutine TIdfAngle_Save

!==============================================================!
!  Subroutine TIdfDihedral_Construct                             !
!==============================================================!

  subroutine TIdfDihedral_Construct( this )

    implicit none

    ! Declare arguments
    type(TIdfDihedral) :: this

    ! Declare local variables
    integer           :: i
    integer           :: stat
    character(2)      :: integ

    ! Read site parameters
    call FileReadParameter_IOBuffer( potmodFile%iounit, IdDihedral_Sites )
    read( IOBuffer, * ) this%SiteId(1), this%SiteId(2), this%SiteId(3), this%SiteId(4)
    call FileReadParameter( this%nmax, potmodFile%iounit, IdDihedral_nmax, .false. )

    if (this%nmax > 0 ) then
      allocate( this%ForConst(1:this%nmax+1), STAT = stat )
      call AllocationError( stat, 'dihedral ForConst for internal degrees of freedom', this%nmax )
      allocate( this%gamma0(1:this%nmax+1), STAT = stat )
      call AllocationError( stat, 'dihedral gamm0 for internal degrees of freedom', this%nmax )
    else
      allocate( this%ForConst(1), STAT = stat )
      call AllocationError( stat, 'dihedral ForConst for internal degrees of freedom', this%nmax )
      allocate( this%gamma0(1), STAT = stat )
      call AllocationError( stat, 'dihedral gamm0 for internal degrees of freedom', this%nmax )
    end if

    call FileReadParameter( this%ForConst(1), potmodFile%iounit, IdDihedral_ForConst//'0', .false. )
    call FileReadParameter(this%gamma0(1), potmodFile%iounit, IdDihedral_gamma0//'0',.false. )
    if (this%nmax > 0) then
      do i = 1,(this%nmax)
        if (i < 10) then
          write(integ,'(I1)') i
        else
          write(integ,'(I2)') i
        end if
        call FileReadParameter( this%ForConst(i+1), potmodFile%iounit, IdDihedral_ForConst//trim(integ), .false. )
        call FileReadParameter(this%gamma0(i+1), potmodFile%iounit, IdDihedral_gamma0//trim(integ), .false. )
      end do
    end if

    if (LJEl14 .and. (this%nmax > 0)) then
      this%ScaleLJ14 = 0.0
      this%ScaleEl14 = 0.0
      call FileReadParameter( this%ScaleLJ14, potmodFile%iounit, IdDihedral_ScaleLJ14, .false. )
      call FileReadParameter( this%ScaleEl14, potmodFile%iounit, IdDihedral_ScaleEl14, .false. )
    end if

    ! Convert to SI units
    this%ForConst(:) = this%ForConst(:) * kBoltzmann

    ! Convert to derived units
    this%gamma0(:) = this%gamma0(:)*Pi/180
    this%ForConst(:) = this%ForConst(:) / UnitEnergy

end subroutine TIdfDihedral_Construct


!==============================================================!
!  Subroutine TIdfDihedral_Destruct                              !
!==============================================================!

  subroutine TIdfDihedral_Destruct( this )

    implicit none

    ! Declare arguments
    type(TIdfDihedral) :: this

    ! Destroy site
    continue

  end subroutine TIdfDihedral_Destruct


!==============================================================!
!  Subroutine TIdfDihedral_Save                                    !
!==============================================================!

  subroutine TIdfDihedral_Save( this )

    implicit none

    ! Declare arguments
    type(TIdfDihedral) :: this

    ! Declare local variables
    integer           :: i

    ! Save site parameters
    write( IOBuffer, '(4I3)' ) this%SiteId(1), this%SiteId(2), this%SiteId(3), this%SiteId(4)
    call FileWriteParameter( normalFile%iounit, IdDihedral_Sites )
    write( IOBuffer, '(G20.10)' ) this%nmax
    call FileWriteParameter( normalFile%iounit, IdDihedral_nmax )
    do i= 1,this%nmax+1
      write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&       this%ForConst(i) * UnitEnergy / kBoltzmann, this%ForConst(i)
      call FileWriteParameter( normalFile%iounit, IdDihedral_ForConst)
      write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&       this%gamma0(i)*180/Pi, this%gamma0(i)
      call FileWriteParameter( normalFile%iounit, IdDihedral_gamma0 )
    end do
    if (LJEl14 .and. (this%nmax .ge. 0)) then
      write( IOBuffer, '(G20.10)' ) this%ScaleLJ14
      call FileWriteParameter( normalFile%iounit, IdDihedral_ScaleLJ14 )
      write( IOBuffer, '(G20.10)' ) this%ScaleEl14
      call FileWriteParameter( normalFile%iounit, IdDihedral_ScaleEl14 )
    end if

  end subroutine TIdfDihedral_Save





end module ms2_idf
