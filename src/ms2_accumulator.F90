!==============================================================!
!  MOLECULAR SIMULATION PROGRAM ms2 Version 5.0                !
!  (c) 2025 by RPTU Kaiserslautern / TU Berlin                 !
!      P.O. Box 67653                                          !
!      67653 Kaiserslautern                                    !
!==============================================================!
!  Module ms2_accumulator                                      !
!  Contains TAccumulator object                                !
!==============================================================!

!****************************************************************
!* Updates and auxiliary routines are available from            *   
!* http://www.ms-2.de                                           *   
!****************************************************************

#ifndef ARCH
#define ARCH    0
#define MPI_VER 0
#endif

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2_accumulator.F90...'
#endif


module ms2_accumulator

#if MPI_VER > 0 && defined(MPI_USE_MODULE)
  use mpi_f08
#endif

  use ms2_global



!==============================================================!
!  Type TAccumulator                                           !
!==============================================================!

  type TAccumulator

    ! Block sum
    real(RK), pointer, contiguous :: BlockSum(:) => NULL()

    ! Number of summed values in block
    integer , pointer, contiguous :: NBlockSum(:) => NULL()

#if MPI_VER > 0        
    ! MC communication COL_DEBUG
    real(RK), pointer, contiguous :: BlockSumGathered(:) => NULL()

    ! MC communication COL_DEBUG
    integer , pointer, contiguous :: NBlockSumGathered(:) => NULL()
#endif  

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

contains


!==============================================================!
!  Subroutine TAccumulator_Construct                           !
!==============================================================!

  subroutine TAccumulator_Construct( this, UpdateByAverage, trans, kbi )

    implicit none

    ! Declare arguments
    type(TAccumulator)            :: this
    logical, intent(in)           :: UpdateByAverage
    logical, intent(in), optional :: trans
    logical, intent(in), optional :: kbi

    ! Set method of updating
    this%UpdateByAverage = UpdateByAverage

    ! Initialize
    this%TotalSum = 0._RK
    this%NTotalSum = 0

    ! Allocate arrays
    !DC EDIT- rework to preserv functionality while agreeing with FORTRAN standard on unambiguous evaluation order of "if" statement
    if (present(trans)) then
      if (trans .eqv. .true.) then
        call Allocate( this, trans )
        return
      end if
    end if
    if (present(kbi)) then
      if (kbi .eqv. .true.) then
        call Allocate( this, trans, kbi )
        return
      end if 
    end if

    call Allocate( this )
    

  end subroutine TAccumulator_Construct


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


!==============================================================!
!  Subroutine TAccumulator_Allocate                            !
!==============================================================!

  subroutine TAccumulator_Allocate( this, trans, kbi )

    implicit none

    ! Include MPI header
#if MPI_VER > 0 && !defined(MPI_USE_MODULE)
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TAccumulator)            :: this
    logical, intent(in), optional :: trans
    logical, intent(in), optional :: kbi

    ! Declare local variables
    integer :: stat, i

    i = NBlocksMax
    !DC EDIT- rework to preserv functionality while agreeing with FORTRAN standard on unambiguous evaluation order of "if" statement
#if TRANS == 1
    if (present(trans)) then
      if (trans .eqv. .true.) then
        i = NBlocksMaxCF
      end if
    end if
#endif
    if (present(kbi)) then
      if (kbi .eqv. .true.) then
        i = NBlocksMaxKBI
      end if
    end if
        
    ! Allocate arrays
    allocate( this%BlockSum( i ), STAT = stat )
    call AllocationError( stat, 'output blocks', i )
    allocate( this%NBlockSum( i ), STAT = stat )
    call AllocationError( stat, 'output blocks', i )
    this%BlockSum  = 0._RK
    this%NBlockSum = 0

#if MPI_VER > 0
    if ( SimulationType .eq. MonteCarlo ) then
      ! Allocate arrays for MC communication COL_DEBUG
      if ( mpiMCCommonGroups > 0 ) then !Production cycles are distributed over mpiMCCommonGroups
        allocate( this%BlockSumGathered( i*mpiMCCommonGroups ), STAT = stat )
      else
        allocate( this%BlockSumGathered( i*NProcs ), STAT = stat )
      endif     
      call AllocationError( stat, 'output blocks', i )
      this%BlockSumGathered = 0._RK
   
      if ( mpiMCCommonGroups > 0 ) then
        allocate( this%NBlockSumGathered( i*mpiMCCommonGroups ), STAT = stat )
      else
        allocate( this%NBlockSumGathered( i*NProcs ), STAT = stat )
      endif
      call AllocationError( stat, 'output blocks', i )
      this%NBlockSumGathered = 0._RK
    endif
#endif   

  end subroutine TAccumulator_Allocate


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

  subroutine TAccumulator_Update( this, Value, Mmess, kbi )

    implicit none

    ! Declare arguments
    type(TAccumulator)            :: this
    real(RK), intent(in)          :: Value
    integer, intent(in), optional :: Mmess
    logical, intent(in), optional :: kbi

    ! Declare local variables
    integer :: i, j, k
    

    i = Step
    j = BlockSize
    k = NBlocks
#if TRANS == 1
    if (present(Mmess)) then
      i = Mmess
      j = BlockSizeCF
      k = NBlocksCF
    end if
#endif
    !DC EDIT- rework to preserv functionality while agreeing with FORTRAN standard on unambiguous evaluation order of "if" statement
    if (present(kbi)) then
      if (kbi .eqv. .true.) then
        i = Step
        j = BlockSizeKBI
        k = NBlocksKBI
      end if
    end if

    ! Update sums and calculate average
    if( this%UpdateByAverage ) then
      if( mod( i, j ) == 0) then
        this%BlockSum(k) = k * j * Value - this%TotalSum
        this%NBlockSum(k) = j
        this%TotalSum = k * j * Value
        this%NTotalSum = k * j
        this%Average = this%TotalSum / real( this%NTotalSum, RK )
        this%BlockAverage = this%BlockSum(k) / real( this%NBlockSum(k), RK )
      end if
    else
      this%BlockSum(k) = this%BlockSum(k) + Value
      this%NBlockSum(k) = this%NBlockSum(k) + 1
      this%TotalSum = this%TotalSum + Value
      this%NTotalSum = this%NTotalSum + 1
      this%Average = this%TotalSum / real( this%NTotalSum, RK )
      this%BlockAverage = this%BlockSum(k) / real( this%NBlockSum(k), RK )
    end if

  end subroutine TAccumulator_Update


!==============================================================!
!  Subroutine TAccumulator_Error                               !
!==============================================================!

  subroutine TAccumulator_Error( this, trans, kbi )

    implicit none
    
#if MPI_VER > 0 && !defined(MPI_USE_MODULE)
  include 'mpif.h'
#endif

    ! Declare arguments
    type(TAccumulator)            :: this
    logical, intent(in), optional :: trans
    logical, intent(in), optional :: kbi

    ! Declare local variables
    real(RK), dimension(:), allocatable :: Tau
    real(RK) :: BlockAverage
    real(RK) :: sx1, sx2, sxy
    real(RK) :: TauSum, TauInf
    integer :: i, j, m, n, stat
#if MPI_VER > 0
    real(RK) :: ReducedAverage
#endif

    m = NBlockSizes
    n = NBlocks
    !DC EDIT- rework to preserv functionality while agreeing with FORTRAN standard on unambiguous evaluation order of "if" statement
#if TRANS == 1
    if (present(trans)) then
      if (trans .eqv. .true.) then
        m = NBlockSizesCF
        n = NBlocksCF
      end if
    end if

#endif
    if (present(kbi)) then
      if (kbi .eqv. .true.) then
        m = NBlockSizesKBI
        n = NBlocksKBI
      end if
    end if

    allocate(Tau(max(m,1)),STAT=stat)


#if MPI_VER > 0
    if ( SimulationType .eq. MonteCarlo .and. .not. present(kbi)) then
      if ( mpiMCCommonGroups > 0 ) then !gather and reduce accumulated values of each group to the head (RootProc_MCCom) of each group heads (RootProc)
        if (RootProc) then
          call MPI_Gather(this%BlockSum(1:(n/mpiMCCommonGroups)),n/mpiMCCommonGroups, MPI_RK , &
&           this%BlockSumGathered(1:n), n/mpiMCCommonGroups,MPI_RK,NRootProc_MCCom,MCCommonGroups_R,ierror )
          call MPI_Gather(this%NBlockSum(1:(n/mpiMCCommonGroups)),n/mpiMCCommonGroups, MPI_INTEGER , &
&           this%NBlockSumGathered(1:n), n/mpiMCCommonGroups,MPI_INTEGER,NRootProc_MCCom,MCCommonGroups_R,ierror )
          call MPI_Reduce( this%Average,ReducedAverage, 1, MPI_RK, MPI_SUM, &
&           NRootProc_MCCom, MCCommonGroups_R, ierror )
        endif
      else
        call MPI_Gather(this%BlockSum(1:(n/NProcs)),n/NProcs, MPI_RK , &
&         this%BlockSumGathered(1:n), n/NProcs,MPI_RK,NRootProc,Communicator,ierror )
        call MPI_Gather(this%NBlockSum(1:(n/NProcs)),n/NProcs, MPI_INTEGER , this%NBlockSumGathered(1:n), & 
&         n/NProcs,MPI_INTEGER,NRootProc,Communicator,ierror )
        call MPI_Reduce( this%Average,ReducedAverage, 1, MPI_RK, MPI_SUM, &
&         NRootProc, Communicator, ierror )
      endif

      !Carefull: This if statement should remain as is, 
      ! because in the MC parallelization, every processor is treated as root
      ! this means that NRootProc=0 for each process
      if (Nproc /= NRootProc) return
      
      if ( mpiMCCommonGroups > 0 ) then
        if ( .not. RootProc_MCCom ) return !=RootProc_W
        this%Average=ReducedAverage/mpiMCCommonGroups
      else        
        this%Average=ReducedAverage/NProcs
      endif
      ! Calculate variance
      Tau = 0._RK
      do i = 1, m
        do j = i, n, i
          BlockAverage = sum( this%BlockSumGathered(j - i + 1:j) ) / real( sum(this%NBlockSumGathered (j - i + 1:j) ), RK )
          Tau(i) = Tau(i) + (BlockAverage - this%Average)**2
        end do
#ifdef _PGF
        ! Call write to prevent vectorization of loop (a bug in pgi compiler)
        write( IOBuffer, '("Prevent loop vectorization")' )
#endif
        Tau(i) = Tau(i) / real( (n / i), RK )
      end do

    else
      Tau = 0._RK
      do i = 1, m
        do j = i, n, i
          BlockAverage = sum( this%BlockSum(j - i + 1:j) ) / real( sum( this%NBlockSum(j - i + 1:j) ), RK )
          Tau(i) = Tau(i) + (BlockAverage - this%Average)**2
        end do
#ifdef _PGF
        ! Call write to prevent vectorization of loop (a bug in pgi compiler)
        write( IOBuffer, '("Prevent loop vectorization")' )
#endif
        Tau(i) = Tau(i) / real( (n / i), RK )
      end do
    endif
#else
    ! Calculate variance
    Tau = 0._RK
    do i = 1, m
      do j = i, n, i
        BlockAverage = sum( this%BlockSum(j - i + 1:j) ) / real( sum( this%NBlockSum(j - i + 1:j) ), RK ) ! Michael Sch.: for trans NBlockSum was with i*BlockSizeCF before...
        Tau(i) = Tau(i) + (BlockAverage - this%Average)**2
      end do
#ifdef _PGF
      ! Call write to prevent vectorization of loop (a bug in pgi compiler)
      write( IOBuffer, '("Prevent loop vectorization")' )
#endif
      Tau(i) = Tau(i) / real( (n / i), RK )
    end do
#endif

    this%Variance = Tau(1)
    if( this%Variance == 0._RK ) return
    Tau(1:m) = Tau(1:m) / this%Variance
    sx1 = 0._RK
    sx2 = 0._RK
    sxy = 0._RK
    do i = 1, m
      sx1 = sx1 + 1._RK / i
      sx2 = sx2 + 1._RK / i**2
      sxy = sxy + Tau(i)
      Tau(i) = i * Tau(i)
    end do
    TauSum = sum( Tau(1:m) )
    TauInf = (TauSum * sx2 - sx1* sxy) / (m * sx2 - sx1**2)
    TauInf = max( TauInf, TauSum / m )
    this%Variance = sqrt( this%Variance / n * TauInf )
    
    deallocate(Tau,STAT=stat)

  end subroutine TAccumulator_Error


  subroutine writeAverages(this, resultFile, runaveFile, optionalFormatString, parallelMC)

    implicit none

    type(TAccumulator) :: this
    type(TFile)        :: resultFile, runaveFile
    logical            :: parallelMC
    character(len=*), intent(in), optional :: optionalFormatString
    character(:), allocatable :: formatString

    if (PRESENT(optionalFormatString)) then
        formatString = optionalFormatString
    else
        formatString = '(" ",F10.5)'
    end if

#if MPI_VER > 0
    if (parallelMC) then

        write( IOBuffer, formatString) this%BlockAverage
        call FileWriteNoAdvance_parallel(resultFile)

        write( IOBuffer, formatString) this%Average
        call FileWriteNoAdvance_parallel(runaveFile)

    else
#endif

    write( IOBuffer, formatString) this%BlockAverage
    call FileWriteNoAdvance(resultFile)

    write( IOBuffer, formatString) this%Average
    call FileWriteNoAdvance(runaveFile)

#if MPI_VER > 0
    end if
#endif

  end subroutine writeAverages


  subroutine writeDimlessValueWithError(errorsFile, variableName, average, error, reducedTitle)

    implicit none

    type(TFile), intent(in)        :: errorsFile
    character(:), allocatable    :: formatString
    character(len=*), intent(in)    :: variableName
    logical, optional :: reducedTitle
    real(RK) :: average, error


    if (present(reducedTitle) .and. reducedTitle) then

        formatString = '("'//variableName//'", T13, "Dimensionless:", 2F20.9)'
    else

        formatString = '("'//variableName//'", T13, "Dimensionless, residual:", 2F20.9)'
    end if

    write( IOBuffer, formatString) average, error
    call FileWrite(errorsFile)
    call FileWriteBlank(errorsFile)

  end subroutine writeDimlessValueWithError


  subroutine writeAverageAndVariance(this, variableName, errorsFile, reducedTitle)

    implicit none

    type(TAccumulator) :: this
    type(TFile), intent(in)        :: errorsFile
    character(len=*), intent(in)    :: variableName
    logical, optional :: reducedTitle

    call writeDimlessValueWithError(errorsFile, variableName, this%Average, this%Variance, present(reducedTitle))

  end subroutine writeAverageAndVariance




!==============================================================!
!  Subroutine TAccumulator_RestartSave                         !
!==============================================================!

  subroutine TAccumulator_RestartSave( this, trans, kbi )

    implicit none

    ! Declare arguments
    type(TAccumulator)             :: this
    logical, intent(in), optional  :: trans
    logical, intent(in), optional  :: kbi

    ! Declare local variables
    integer :: i, j

    ! Check for root process
    if( .not. RootProc ) return

    j = NBlocks
#if TRANS == 1
    !DC EDIT- rework to preserv functionality while agreeing with FORTRAN standard on unambiguous evaluation order of "if" statement
    if (present(trans)) then
      if (trans .eqv. .true.) then
        j = NBlocksCF         
      end if
    end if
#endif
    if (present(kbi) ) then
      if (kbi .eqv. .true.) then
         j = NBlocksKBI
      end if
    end if
    ! Save contents to restart file
    write( restartFile%iounit, '(I10)' ) j
    if( j > 0 ) write( restartFile%iounit, '(ES20.12E3, ";", I10)' ) &
&     ( this%BlockSum(i), this%NBlockSum(i), i = 1, j )

  end subroutine TAccumulator_RestartSave


!==============================================================!
!  Subroutine TAccumulator_RestartRead                         !
!==============================================================!

  subroutine TAccumulator_RestartRead( this, kbi )

    implicit none

    ! Include MPI header
#if MPI_VER > 0 && !defined(MPI_USE_MODULE)
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TAccumulator) :: this
    logical, intent(in), optional  :: kbi

    ! Declare local variables
    integer :: i, j

    ! Check for root process
    if( RootProc ) then
      ! Read contents from restart file
      read( restartFile%iounit, '(I10)' ) i
      !read( restartFile%iounit, '(ES20.12E3, X, I10)' ) ( this%BlockSum(j), this%NBlockSum(j), j = 1, i )
      do j = 1, i ! should be equivalent to the previous line, which produced an "input conversion error"
        read( restartFile%iounit, '(ES20.12E3, 1X, I10)' ) this%BlockSum(j), this%NBlockSum(j)
      end do
    endif
    
#if MPI_VER >0
    if( SimulationType .eq. MonteCarlo .and. .not. present(kbi)) then
      call MPI_Bcast( this%BlockSum(:), size( this%BlockSum ), MPI_RK, NRootProc, Communicator, ierror )
      call MPI_Bcast( this%NBlockSum(:), size( this%NBlockSum ), MPI_INTEGER, NRootProc, Communicator, ierror )
      call MPI_Bcast( i, 1, MPI_INTEGER, NRootProc, Communicator, ierror )
    elseif ( .not. RootProc) then
      return
    endif
#endif     

    this%TotalSum = sum( this%BlockSum(1:i) )
    this%NTotalSum = sum( this%NBlockSum(1:i) )

    ! Calculate average
    this%Average = this%TotalSum / real( this%NTotalSum, RK )
    this%BlockAverage = this%BlockSum(i) / real( this%NBlockSum(i), RK )

  end subroutine TAccumulator_RestartRead



end module ms2_accumulator
