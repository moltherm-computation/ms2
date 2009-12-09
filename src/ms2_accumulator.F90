!==============================================================!
!  MOLECULAR SIMULATION PROGRAM MS2 Version 1.1 v12            !
!  (c) 2001 by Sergey Lishchuk, ITT                            !
!  (c) 2007 by Bernhard Eckl, ITT                              !
!==============================================================!
!  Module ms2_accumulator                                      !
!  Contains TAccumulator object                                !
!==============================================================!

#ifndef ARCH
#define ARCH    0
#define FORTRAN 90
#define MPI_VER 0
#endif

#ifndef TRANS
#define TRANS 0
#endif

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2_accumulator.F90...'
#endif

module ms2_accumulator

  use ms2_global



!==============================================================!
!  Type TAccumulator                                           !
!==============================================================!

  type TAccumulator

    ! Block sum
    real(RK), pointer :: BlockSum(:)

    ! Number of summed values in block
    integer , pointer :: NBlockSum(:)

    ! Total sum and average
    real(RK) :: TotalSum, Average, BlockAverage

    ! Total number of summed values
    integer :: NTotalSum

    ! Variance
    real(RK) :: Variance

    ! Method of updating
    logical :: UpdateByAverage

  end type TAccumulator

  interface Construct
    module procedure TAccumulator_Construct
  end interface

  interface Destruct
    module procedure TAccumulator_Destruct
  end interface

  interface Allocate
    module procedure TAccumulator_Allocate
  end interface

  interface Deallocate
    module procedure TAccumulator_Deallocate
  end interface

  interface Reset
    module procedure TAccumulator_Reset
  end interface

  interface Update
    module procedure TAccumulator_Update
  end interface

  interface Error
    module procedure TAccumulator_Error
  end interface

  interface RestartSave
    module procedure TAccumulator_RestartSave
  end interface

  interface RestartRead
    module procedure TAccumulator_RestartRead
  end interface


!===============================================================
!TRANSPORT_start
!==============================================================!
!  Type TAccumulatorCF                                         !
!==============================================================!
  type TAccumulatorCF
    ! Block sum
    real(RK), pointer :: BlockSum(:)

    ! Total sum
    real(RK) :: TotalSum

    ! Average and variance
    real(RK) :: Average, Variance

    ! Method of updating
    logical :: UpdateByAverage
  end type TAccumulatorCF

  interface ConstructCF
    module procedure TAccumulatorCF_Construct
  end interface

  interface DestructCF
    module procedure TAccumulatorCF_Destruct
  end interface

  interface AllocateCF
    module procedure TAccumulatorCF_Allocate
  end interface

  interface DeallocateCF
    module procedure TAccumulatorCF_Deallocate
  end interface

  interface UpdateCF
    module procedure TAccumulatorCF_Update
  end interface

  interface ErrorCF
    module procedure TAccumulatorCF_Error
  end interface
  
  interface RestartSaveCF
    module procedure TAccumulator_RestartSaveCF
  end interface

  interface RestartReadCF
    module procedure TAccumulator_RestartReadCF
  end interface

contains

!TRANSPORT_END

!==============================================================!
!  Subroutine TAccumulator_Construct                           !
!==============================================================!

  subroutine TAccumulator_Construct( this, UpdateByAverage )

    implicit none

    ! Declare arguments
    type(TAccumulator)  :: this
    logical, intent(in) :: UpdateByAverage

    ! Set method of updating
    this%UpdateByAverage = UpdateByAverage

    ! Allocate arrays
    call Allocate( this )

  end subroutine TAccumulator_Construct
!TRANSPORT_start
!==============================================================!
!  Subroutine TAccumulatorCF_Construct                         !
!==============================================================!

  subroutine TAccumulatorCF_Construct( this, UpdateByAverage )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TAccumulatorCF):: this
    logical, intent(in) :: UpdateByAverage 

    ! Set method of updating
    this%UpdateByAverage = UpdateByAverage

    ! Allocate arrays
    call AllocateCF( this )

  end subroutine TAccumulatorCF_Construct

!TRANSPORT_END
!==============================================================!
!  Subroutine TAccumulator_Destruct                            !
!==============================================================!

  subroutine TAccumulator_Destruct( this )

    implicit none

    ! Declare arguments
    type(TAccumulator) :: this

    ! Deallocate arrays
    call Deallocate( this )

  end subroutine TAccumulator_Destruct

!TRANSPORT_start
!==============================================================!
!  Subroutine TAccumulatorCF_Destruct                            !
!==============================================================!

  subroutine TAccumulatorCF_Destruct( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TAccumulatorCF) :: this

    ! Deallocate arrays
    call DeallocateCF( this )

  end subroutine TAccumulatorCF_Destruct


!TRANSPORT_END
!==============================================================!
!  Subroutine TAccumulator_Allocate                            !
!==============================================================!

  subroutine TAccumulator_Allocate( this )

    implicit none

    ! Declare arguments
    type(TAccumulator) :: this

    ! Declare local variables
    integer :: stat

    ! Allocate arrays
    allocate( this%BlockSum( NBlocksMax ), STAT = stat )
    call AllocationError( stat, 'output blocks', NBlocksMax )
    allocate( this%NBlockSum( NBlocksMax ), STAT = stat )
    call AllocationError( stat, 'output blocks', NBlocksMax )

  end subroutine TAccumulator_Allocate
!TRANSPORT_start
!==============================================================!
!  Subroutine TAccumulatorCF_Allocate                            !
!==============================================================!

  subroutine TAccumulatorCF_Allocate( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TAccumulatorCF) :: this

    ! Declare local variables
    integer :: stat
#if TRANS==1
    ! Allocate arrays
    allocate( this%BlockSum( NBlocksMaxCF ), STAT = stat )
    call AllocationError( stat, 'output blocks', NBlocksMaxCF )
#endif
  end subroutine TAccumulatorCF_Allocate
!TRANSPORT_END
!==============================================================!
!  Subroutine TAccumulator_Deallocate                          !
!==============================================================!

  subroutine TAccumulator_Deallocate( this )

    implicit none

    ! Declare arguments
    type(TAccumulator) :: this

    ! Deallocate arrays
    if( associated( this%BlockSum ) ) then
      deallocate( this%BlockSum )
    end if
    if( associated( this%NBlockSum ) ) then
      deallocate( this%NBlockSum )
    end if

  end subroutine TAccumulator_Deallocate
!TRANSPORT_start
!==============================================================!
!  Subroutine TAccumulatorCF_Deallocate                          !
!==============================================================!

  subroutine TAccumulatorCF_Deallocate( this )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TAccumulatorCF) :: this


    ! Deallocate arrays
    if( associated( this%BlockSum ) ) then
       deallocate( this%BlockSum )
    end if

  end subroutine TAccumulatorCF_Deallocate

!TRANSPORT_END
!==============================================================!
!  Subroutine TAccumulator_Reset                               !
!==============================================================!

  subroutine TAccumulator_Reset( this )

    implicit none

    ! Declare arguments
    type(TAccumulator) :: this

    ! Zero sums
    this%BlockSum(:) = 0._RK
    this%TotalSum = 0._RK
    this%NBlockSum(:) = 0
    this%NTotalSum = 0

    ! Zero average
    this%Average = 0._RK
    this%BlockAverage = 0._RK

  end subroutine TAccumulator_Reset



!==============================================================!
!  Subroutine TAccumulator_Update                              !
!==============================================================!

  subroutine TAccumulator_Update( this, Value )

    implicit none

    ! Declare arguments
    type(TAccumulator)   :: this
    real(RK), intent(in) :: Value

    ! Update sums and calculate average
    if( this%UpdateByAverage ) then
      if( mod( Step, BlockSize ) == 0) then
        this%BlockSum(NBlocks) = NBlocks * BlockSize * Value - this%TotalSum
        this%NBlockSum(NBlocks) = BlockSize
        this%TotalSum = NBlocks * BlockSize * Value
        this%NTotalSum = NBlocks * BlockSize
        this%Average = this%TotalSum / real( this%NTotalSum, RK )
        this%BlockAverage = this%BlockSum(NBlocks) &
&         / real( this%NBlockSum(NBlocks), RK )
      end if
    else
      this%BlockSum(NBlocks) = this%BlockSum(NBlocks) + Value
      this%NBlockSum(NBlocks) = this%NBlockSum(NBlocks) + 1
      this%TotalSum = this%TotalSum + Value
      this%NTotalSum = this%NTotalSum + 1
      this%Average = this%TotalSum / real( this%NTotalSum, RK )
      this%BlockAverage = this%BlockSum(NBlocks) &
&       / real( this%NBlockSum(NBlocks), RK )
    end if

  end subroutine TAccumulator_Update

!TRANSPORT_start
!==============================================================!
!  Subroutine TAccumulatorCF_Update                              !
!==============================================================!

  subroutine TAccumulatorCF_Update( this, Value , Mmess )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TAccumulatorCF)    :: this
    real(RK), intent(in)    :: Value
    integer, intent(in)     :: Mmess
#if TRANS==1
    ! Nullify total sum
    if( Mmess == 1 ) this%TotalSum = 0._RK
    if( Mmess == BlockSizeCF ) this%TotalSum = 0._RK

    ! Update sums
    if( this%UpdateByAverage ) then
      if( mod( Mmess, BlockSizeCF ) == 0) then
        this%BlockSum(NBlocksCF) = NBlocksCF * BlockSizeCF * Value - this%TotalSum
        this%TotalSum = NBlocksCF * BlockSizeCF * Value
      end if
    else
      if( mod( Mmess, BlockSizeCF ) == 0 ) this%BlockSum(NBlocksCF) = 0._RK
      this%BlockSum(NBlocksCF) = this%BlockSum(NBlocksCF) + Value
      this%TotalSum = this%TotalSum + Value
    end if
#endif
  end subroutine TAccumulatorCF_Update

!TRANSPORT_END
!==============================================================!
!  Subroutine TAccumulator_Error                               !
!==============================================================!

  subroutine TAccumulator_Error( this )

    implicit none

    ! Declare arguments
    type(TAccumulator) :: this

    ! Declare local variables
    real(RK) :: Tau(NBlockSizes)
    real(RK) :: BlockAverage
    real(RK) :: sx1, sx2, sxy
    real(RK) :: TauSum, TauInf
    integer :: i, j

    ! Calculate variance
    Tau = 0._RK
    do i = 1, NBlockSizes
      do j = i, NBlocks, i
        BlockAverage = sum( this%BlockSum(j - i + 1:j) ) &
&         / real( sum( this%NBlockSum(j - i + 1:j) ), RK )
        Tau(i) = Tau(i) + (BlockAverage - this%Average)**2
      end do
#ifdef _PGF
      ! Call write to prevent vectorization of loop (a bug in pgi compiler)
      write( IOBuffer, '("Prevent loop vectorization")' )
#endif
      Tau(i) = Tau(i) / real( (NBlocks / i), RK )
    end do
    this%Variance = Tau(1)
    if( this%Variance == 0._RK ) return
    Tau(1:NBlockSizes) = Tau(1:NBlockSizes) / this%Variance
    sx1 = 0._RK
    sx2 = 0._RK
    sxy = 0._RK
    do i = 1, NBlockSizes
      sx1 = sx1 + 1._RK / i
      sx2 = sx2 + 1._RK / i**2
      sxy = sxy + Tau(i)
      Tau(i) = i * Tau(i)
    end do
    TauSum = sum( Tau(1:NBlockSizes) )
    TauInf = (TauSum * sx2 - sx1* sxy) / (NBlockSizes * sx2 - sx1**2)
    TauInf = max( TauInf, TauSum / NBlockSizes )
    this%Variance = sqrt( this%Variance / NBlocks * TauInf )

  end subroutine TAccumulator_Error

!TRANSPORT_start
!==============================================================!
!  Subroutine TAccumulatorCF_Error                             !
!==============================================================!

  subroutine TAccumulatorCF_Error( this, Mmess )

    implicit none

    ! Include MPI header
#if MPI_VER > 0
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TAccumulatorCF) :: this
#if TRANS==1
    ! Declare local variables
    real(RK) :: Tau(NBlockSizesCF)
#endif
    real(RK) :: BlockAverage
    real(RK) :: sx1, sx2, sxy
    real(RK) :: TauSum, TauInf
    integer :: i, j, k, Mmess
    ! Calculate average
#if TRANS==1
    this%Average = this%TotalSum / Mmess

    ! Calculate variance
    Tau = 0._RK
    do i = 1, NBlockSizesCF
      do j = i, NBlocksCF, i
        BlockAverage = sum( this%BlockSum(j - i + 1:j) ) / real( (i * BlockSizeCF), RK )
        Tau(i) = Tau(i) + (BlockAverage - this%Average)**2
      end do
#ifdef _PGF
      ! Call write to prevent vectorization of loop (a bug in pgi compiler)
      write( IOBuffer, '("Prevent loop vectorization")' )
#endif
      Tau(i) = Tau(i) / real( (NBlocksCF / i), RK )
    end do
    this%Variance = Tau(1)
    if( this%Variance == 0._RK ) return
    Tau(1:NBlockSizesCF) = Tau(1:NBlockSizesCF) / this%Variance
    sx1 = 0._RK
    sx2 = 0._RK
    sxy = 0._RK
    do i = 1, NBlockSizesCF
      sx1 = sx1 + 1._RK / i
      sx2 = sx2 + 1._RK / i**2
      sxy = sxy + Tau(i)
      Tau(i) = i * Tau(i)
    end do
    TauSum = sum( Tau(1:NBlockSizesCF) )
    TauInf = (TauSum * sx2 - sx1* sxy) / (NBlockSizesCF * sx2 - sx1**2)
    TauInf = max( TauInf, TauSum / NBlockSizesCF )
    this%Variance = sqrt( this%Variance / NBlocksCF * TauInf )
#endif
  end subroutine TAccumulatorCF_Error
!TRANSPORT_END
!==============================================================!
!  Subroutine TAccumulator_RestartSave                         !
!==============================================================!

  subroutine TAccumulator_RestartSave( this )

    implicit none

    ! Declare arguments
    type(TAccumulator) :: this

    ! Declare local variables
    integer :: i

    ! Check for root process
    if( .not. RootProc ) return

    ! Save contents to restart file
    write( iounit_restart, '(I10)' ) NBlocks
    if( NBlocks > 0 ) write( iounit_restart, '(ES20.12E3, ";", I10)' ) &
&     ( this%BlockSum(i), this%NBlockSum(i), i = 1, NBlocks )

  end subroutine TAccumulator_RestartSave



!==============================================================!
!  Subroutine TAccumulator_RestartRead                         !
!==============================================================!

  subroutine TAccumulator_RestartRead( this )

    implicit none

    ! Declare arguments
    type(TAccumulator) :: this

    ! Declare local variables
    integer :: i, j

    ! Check for root process
    if( .not. RootProc ) return

    ! Read contents from restart file
    read( iounit_restart, '(I10)' ) i
    read( iounit_restart, '(ES20.12E3, ";", I10)' ) &
&     ( this%BlockSum(j), this%NBlockSum(j), j = 1, i )

    this%TotalSum = sum( this%BlockSum(1:i) )
    this%NTotalSum = sum( this%NBlockSum(1:i) )

    ! Calculate average
    this%Average = this%TotalSum / real( this%NTotalSum, RK )
    this%BlockAverage = this%BlockSum(NBlocks) / real( this%NBlockSum(NBlocks), RK )

  end subroutine TAccumulator_RestartRead

!==============================================================!
!  Subroutine TAccumulator_RestartSaveCF                         !
!==============================================================!

  subroutine TAccumulator_RestartSaveCF( this )

    implicit none

    ! Declare arguments
    type(TAccumulatorCF) :: this
#if  TRANS == 1
    ! Declare local variables
    integer :: j

    ! Check for root process
    if( .not. RootProc ) return

    ! Save contents to restart file
   ! write( iounit_restart, '(I10)' ) NBlocksMaxCF
     
      do j = 1, NBlocksMaxCF
       write( iounit_restart, '(ES20.12E3)' )  this%BLOCKSUM(j)
      end do 
      write( iounit_restart, '(ES20.12E3)' )  this%TOTALSUM   
      write( iounit_restart, '(ES20.12E3)' )  this%AVERAGE 
      write( iounit_restart, '(ES20.12E3)' )  this%VARIANCE 
#endif   

  end subroutine TAccumulator_RestartSaveCF
  
!==============================================================!
!  Subroutine TAccumulator_RestartReadCF                        !
!==============================================================!

  subroutine TAccumulator_RestartReadCF( this )

    implicit none

    ! Declare arguments
    type(TAccumulatorCF) :: this
#if  TRANS == 1
    ! Declare local variables
    integer :: j

    ! Check for root process
    if( .not. RootProc ) return

    ! Read contents from restart file
!    read( iounit_restart, '(I10)' ) NBlocksMaxCF
 

    do j = 1, NBlocksMaxCF
       read( iounit_restart, '(ES20.12E3)' )  this%BLOCKSUM(j)
      end do 
      read( iounit_restart, '(ES20.12E3)' )  this%TOTALSUM   
      read( iounit_restart, '(ES20.12E3)' )  this%AVERAGE 
      read( iounit_restart, '(ES20.12E3)' )  this%VARIANCE 
#endif 
  end subroutine TAccumulator_RestartReadCF

end module ms2_accumulator
