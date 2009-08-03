!==============================================================!
!  MOLECULAR SIMULATION PROGRAM MS2 Version 1.1 v11            !
!  (c) 2001 by Sergey Lishchuk, ITT                            !
!  (c) 2005 by Bernhard Eckl, ITT
!  (c) 2007 by Ekaterina Elts, TUM
!==============================================================!
!  Module ms2_idf                                             !
!  Contains TIdf* objects                                     !
!==============================================================!

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



!==============================================================!
!  Type TIdfBond                                             !
!==============================================================!

  type TIdfBond

    real(RK)          :: ForConst
    real(RK)          :: R0
    integer           :: SiteId1, SiteId2
    integer, pointer  :: NPartMax, NPart
    integer, pointer  :: NPart0, NPart1, NPart2
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:)
    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX2(:), PY2(:), PZ2(:)

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
    integer           :: SiteId1, SiteId2, SiteId3
    integer, pointer  :: NPartMax, NPart
    integer, pointer  :: NPart0, NPart1, NPart2
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:), RX3(:), RY3(:), RZ3(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:), FX3(:), FY3(:), FZ3(:)


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

    real(RK)          :: ForConst
    real(RK)          :: phi
    real(RK)          :: gamma ! phase factor
    integer           :: multi     ! multiplicity
    integer           :: SiteId1, SiteId2, SiteId3, SiteId4
    real(RK)          :: ScaleLJ14
    real(RK)          :: ScaleEl14
    integer, pointer  :: NPartMax, NPart
    integer, pointer  :: NPart0, NPart1, NPart2
    real(RK), pointer :: RX1(:), RY1(:), RZ1(:), RX2(:), RY2(:), RZ2(:), RX3(:), RY3(:), RZ3(:), RX4(:), RY4(:), RZ4(:)
    real(RK), pointer :: FX1(:), FY1(:), FZ1(:), FX2(:), FY2(:), FZ2(:), FX3(:), FY3(:), FZ3(:), FX4(:), FY4(:), FZ4(:)
!    real(RK), pointer :: PX1(:), PY1(:), PZ1(:), PX4(:), PY4(:), PZ4(:)
!    real(RK)          :: Sigma1, Sigma4, Epsilon1, Epsilon4


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
!  Subroutine TIdfBond_Construct                             !
!==============================================================!

  subroutine TIdfBond_Construct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif


    ! Declare arguments
    type(TIdfBond) :: this

    ! Read site parameters
    call FileReadParameter( iounit_potmod, IdBond_Sites )
    read( IOBuffer, * ) this%SiteId1, this%SiteId2
    call FileReadParameter( iounit_potmod, IdBond_ForConst)
    read( IOBuffer, * ) this%ForConst


    ! Convert to SI units
    this%ForConst = this%ForConst * kBoltzmann

    ! Convert to derived units
    this%ForConst = this%ForConst / UnitEnergy


end subroutine TIdfBond_Construct


!==============================================================!
!  Subroutine TIdfBond_Destruct                              !
!==============================================================!

  subroutine TIdfBond_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TIdfBond) :: this

    ! Save site parameters

    write( IOBuffer, '(2I3)' ) this%SiteId1, this%SiteId2
    call FileWriteParameter( iounit_normal, IdBond_Sites )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%R0 * UnitLength / Angstroem, this%R0
    call FileWriteParameter( iounit_normal, IdBond_R0 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%ForConst * UnitEnergy / kBoltzmann, this%ForConst
    call FileWriteParameter( iounit_normal, IdBond_ForConst )

  end subroutine TIdfBond_Save


!==============================================================!
!  Subroutine TIdfAngle_Construct                             !
!==============================================================!

  subroutine TIdfAngle_Construct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif


    ! Declare arguments
    type(TIdfAngle) :: this

    ! Read site parameters
    call FileReadParameter( iounit_potmod, IdAngle_Sites )
    read( IOBuffer, * ) this%SiteId1, this%SiteId2, this%SiteId3
    call FileReadParameter( iounit_potmod, IdAngle_ForConst)
    read( IOBuffer, * ) this%ForConst

    ! Convert to SI units
    this%ForConst = this%ForConst * kBoltzmann

    ! Convert to derived units
    this%ForConst = this%ForConst / UnitEnergy

end subroutine TIdfAngle_Construct


!==============================================================!
!  Subroutine TIdfAngle_Destruct                              !
!==============================================================!

  subroutine TIdfAngle_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TIdfAngle) :: this

    ! Save site parameters

    write( IOBuffer, '(3I3)' ) this%SiteId1, this%SiteId2, this%SiteId3
    call FileWriteParameter( iounit_normal, IdAngle_Sites )
    write( IOBuffer, '(G20.10)' ) this%Angle0
    call FileWriteParameter( iounit_normal, IdAngle_Angle0 )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%ForConst * UnitEnergy / kBoltzmann, this%ForConst
    call FileWriteParameter( iounit_normal, IdAngle_ForConst )

  end subroutine TIdfAngle_Save

!==============================================================!
!  Subroutine TIdfDihedral_Construct                             !
!==============================================================!

  subroutine TIdfDihedral_Construct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif


    ! Declare arguments
    type(TIdfDihedral) :: this

    ! Read site parameters
    call FileReadParameter( iounit_potmod, IdDihedral_Sites )
    read( IOBuffer, * ) this%SiteId1, this%SiteId2, this%SiteId3, this%SiteId4
    call FileReadParameter( iounit_potmod, IdDihedral_PotBarrier)
    read( IOBuffer, * ) this%ForConst
    call FileReadParameter( iounit_potmod, IdDihedral_gamma )
    read( IOBuffer, * ) this%gamma
    call FileReadParameter( iounit_potmod, IdDihedral_n )
    read( IOBuffer, * ) this%multi
    if (LJEl14 .and. this%multi) then
      call FileReadParameter( iounit_potmod, IdDihedral_ScaleLJ14 )
      read( IOBuffer, * ) this%ScaleLJ14
      call FileReadParameter( iounit_potmod, IdDihedral_ScaleEl14 )
      read( IOBuffer, * ) this%ScaleEl14
    end if

    ! Convert to SI units
    this%ForConst = this%ForConst * kBoltzmann

    ! Convert to derived units
    this%ForConst = this%ForConst / UnitEnergy

end subroutine TIdfDihedral_Construct


!==============================================================!
!  Subroutine TIdfDihedral_Destruct                              !
!==============================================================!

  subroutine TIdfDihedral_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

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

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TIdfDihedral) :: this

    ! Save site parameters

    write( IOBuffer, '(4I3)' ) this%SiteId1, this%SiteId2, this%SiteId3, this%SiteId4
    call FileWriteParameter( iounit_normal, IdDihedral_Sites )
    write( IOBuffer, '(G20.10, T32, "# reduced value: ", G20.10)' ) &
&     this%ForConst * UnitEnergy / kBoltzmann, this%ForConst
    call FileWriteParameter( iounit_normal, IdDihedral_PotBarrier)
    write( IOBuffer, '(G20.10)' ) this%gamma
    call FileWriteParameter( iounit_normal, IdDihedral_gamma )
    write( IOBuffer, '(G20.10)' ) this%multi
    call FileWriteParameter( iounit_normal, IdDihedral_n )
    if (LJEl14 .and. this%multi) then
      write( IOBuffer, '(G20.10)' ) this%ScaleLJ14
      call FileWriteParameter( iounit_normal, IdDihedral_ScaleLJ14 )
      write( IOBuffer, '(G20.10)' ) this%ScaleEl14
      call FileWriteParameter( iounit_normal, IdDihedral_ScaleEl14 )
    end if

  end subroutine TIdfDihedral_Save





end module ms2_idf
