!==============================================================!
!  MOLECULAR SIMULATION PROGRAM MS2 Version 1.1 v12            !
!==============================================================!
!> Module: ms2_stopwatch                                       !
!          contains TStopwatch class                           !
!> \author Martin Bernreuther <bernreuther@hlrs.de>            !
!> \date   06.2009                                             !
!==============================================================!
! alternative: http://math.nist.gov/StopWatch/                 !
!==============================================================!

#ifndef ARCH
#define ARCH    0
#define FORTRAN 90
#define MPI_VER 0
#endif

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2_stopwatch.F90...'
#endif

#define STOPWATCH_USE_DATIME

#define STOPWATCH_USE_SYSCLK
!#if ( ARCH == 2 ) && ( defined __GNUC__ || defined __INTEL_COMPILER || defined __PGI )
!#define STOPWATCH_USE_ETIME
!! bad accurancy? -> better use cpu_time?
!! __PGI: only SYSVR4 environments
!! IBM XL: etime_ instead
!! for ARCH==3 (NEC SX) etime(D) returns elapsed time since 1.1.1970!
!#else
#define STOPWATCH_USE_CPUTIME
!#endif

#if MPI_VER
#define USE_MPI
! the MPI2 standard requires the support of a mpi module (see 16.2),
! but some MPI1 distributions don't offer a Fortran90 binding and won't compile
! uncomment the next line, if you have a MPI distribution for MPI>=version 2
!#define MPI_USE_MODULE
#define STOPWATCH_USE_MPIWTIME
#endif

!----------

#if USE_PAPI
! Note, that PAPI calls outside ms2_stopwatch might interfere...
#define STOPWATCH_USE_PAPI
#endif

#ifdef STOPWATCH_USE_MPIWTIME
#ifndef USE_MPI
! STOPWATCH_USE_MPIWTIME requires MPI!
#undef STOPWATCH_USE_MPIWTIME
#else
#ifdef STOPWATCH_USE_CPUTIME
! CPU time not relevant for parallel, distributed memory application
#undef STOPWATCH_USE_CPUTIME
#endif
#endif
#endif


module ms2_stopwatch

#ifdef MPI_USE_MODULE
  use mpi
#endif

  use ms2_global

#if defined STOPWATCH_USE_ETIME && defined __INTEL_COMPILER
  use IFPORT
! __INTEL_COMPILER also offers DTIME and DCLOCK via IFPORT
#endif

!#ifdef STOPWATCH_USE_PAPI
!  #include "fpapi.h"
!#endif

  integer, parameter :: tag_string_length = 64


  !integer, parameter :: CStopwatch_omitDATIME = B'1'
  integer, parameter :: CStopwatch_omitDATIME = 1
  !integer, parameter :: CStopwatch_omitCPUTIME = B'10'
  integer, parameter :: CStopwatch_omitCPUTIME = 2
  !integer, parameter :: CStopwatch_omitSYSCLK = B'100'
  integer, parameter :: CStopwatch_omitSYSCLK = 4
  !integer, parameter :: CStopwatch_omitETIME = B'1000'
  integer, parameter :: CStopwatch_omitETIME = 8
  !integer, parameter :: CStopwatch_omitPAPI = B'10000'
  integer, parameter :: CStopwatch_omitPAPI = 16
  !integer, parameter :: CStopwatch_omitMPIWTIME = B'100000'
  integer, parameter :: CStopwatch_omitMPIWTIME = 32

  !integer, parameter :: CStopwatch_doMPIStartBarrier = B'100000000'
  integer, parameter :: CStopwatch_doMPIStartBarrier = 256
  !integer, parameter :: CStopwatch_doMPIStopBarrier = B'1000000000'
  integer, parameter :: CStopwatch_doMPIStopBarrier = 512
  !integer, parameter :: CStopwatch_doMPIReduce = B'10000000000'
  integer, parameter :: CStopwatch_doMPIReduce = 1024


#ifdef USE_MPI
  integer :: mpi_defaultCommunicator
#endif

!==============================================================!
!  Type TStopwatch                                             !
!==============================================================!

  type TStopwatch

    ! tag
    character(tag_string_length) :: tag_string
    ! options
    integer :: options
    ! Start date and time
#ifdef STOPWATCH_USE_DATIME
    integer, dimension(8) :: datime_array_start, datime_array_stop
#endif
#ifdef STOPWATCH_USE_SYSCLK
    integer :: sysclk_cnt_start, sysclk_cnt_stop
#endif
#if defined STOPWATCH_USE_CPUTIME
    real :: cputime_start, cputime_stop
#endif
#ifdef STOPWATCH_USE_ETIME
    real, dimension(2) :: etime_array_start, etime_array_stop
    real :: etime_sum_start, etime_sum_stop
#endif

#ifdef STOPWATCH_USE_MPIWTIME
    double precision :: wtime_start, wtime_stop
    double precision, dimension(2) :: wtime_diff
#endif
#ifdef USE_MPI
    integer mpi_communicator
    ! actually only used to reduce the wtime diff, but could be used also to reduce other data in the future
    logical mpi_diff_reduced
#endif

#ifdef STOPWATCH_USE_PAPI
    real(kind=4) papi_real_time_start, papi_real_time_stop, papi_proc_time_start, papi_proc_time_stop
    integer(kind=8) papi_flpops_start, papi_flpops_stop
#endif

  end type TStopwatch


  interface Construct
    module procedure TStopwatch_Construct
  end interface

!  interface Destruct
!    module procedure TStopwatch_Destruct
!  end interface

!  interface Timer_Tag
!    module procedure TStopwatch_Tag
!    !module procedure TStopwatch_GetTag, TStopwatch_SetTag
!    ! is it allowed to mix subroutines and functions?
!    ! the second version won't compile on NEC SX-9, but on other platforms...
!  end interface

  interface Timer_setTag
    module procedure TStopwatch_SetTag
  end interface

  interface Timer_getTag
    module procedure TStopwatch_GetTag
  end interface

  interface Timer_setOptions
    module procedure TStopwatch_SetOptions
  end interface

!  interface Timer_getOptions
!    module procedure TStopwatch_GetOptions
!  end interface

  interface Timer_activateOptions
    module procedure TStopwatch_ActivateOptions
  end interface

  interface Timer_deactivateOptions
    module procedure TStopwatch_DeactivateOptions
  end interface

#ifdef STOPWATCH_USE_DATIME
  interface Timer_getDatimeStart
    module procedure TStopwatch_GetDatimeStart
  end interface

  interface Timer_getDatimeStop
    module procedure TStopwatch_GetDatimeStop
  end interface
#endif

#ifdef STOPWATCH_USE_SYSCLK
  interface Timer_getSysClkDiff
    module procedure TStopwatch_GetSysClkDiff
  end interface
#endif

#ifdef USE_MPI
  interface Timer_SetMPIcommunicator
    module procedure TStopwatch_SetMPIcommunicator
  end interface
#endif

  interface start_Timer
    module procedure TStopwatch_Start
  end interface

  interface stop_Timer
    module procedure TStopwatch_Stop
  end interface

  interface logwritestart_Timer
    module procedure TStopwatch_LogWriteStart
  end interface

  interface logwritestop_Timer
    module procedure TStopwatch_LogWriteStop
  end interface

contains


!==============================================================!
!  Subroutine TStopwatch_Construct                            !
!==============================================================!
  !> Constructor (not really needed except to set tag/options, but since every class offers a constructor...)
  !> \param this       ... object  TStopwatch
  !> \param tag_string ... tag     string
  subroutine TStopwatch_Construct( this, tag_string, options )

    implicit none

    ! Declare arguments
    type(TStopwatch) :: this
    character(*), intent(in), optional :: tag_string
    integer, intent(in), optional :: options

    ! Declare local variables
#ifdef STOPWATCH_USE_PAPI
    integer papi_check
    real(kind=4) papi_mflops
#endif

    if( present( tag_string ) ) then
      this%tag_string = tag_string
    else
      this%tag_string = ""
    end if

    if( present( options ) ) then
      call TStopwatch_SetOptions( this, options )
    else
      this%options = 0
#ifdef USE_MPI
      call TStopwatch_SetOptions( this, &
&                                 CStopwatch_omitCPUTIME+CStopwatch_omitSYSCLK &
&                                +CStopwatch_doMPIReduce )
#endif
    end if

    !call TStopwatch_Start(this)

    this%datime_array_start = 0
    this%datime_array_stop = 0

#ifdef STOPWATCH_USE_SYSCLK
    this%sysclk_cnt_start = 0
    this%sysclk_cnt_stop = 0
#endif
#ifdef STOPWATCH_USE_CPUTIME
    this%cputime_start = 0
    this%cputime_stop = 0
#endif
#ifdef STOPWATCH_USE_ETIME
    this%etime_array_start = 0
    this%etime_array_stop = 0
    this%etime_sum_start = 0
    this%etime_sum_stop = 0
#endif
#ifdef STOPWATCH_USE_MPIWTIME
    this%wtime_start = 0
    this%wtime_stop = 0
    this%wtime_diff = 0
#endif
#ifdef USE_MPI
  !this%mpi_communicator=MPI_COMM_WORLD
  this%mpi_communicator=Communicator
  this%mpi_diff_reduced=.FALSE.
#endif
#ifdef STOPWATCH_USE_PAPI
    this%papi_real_time_start = 0.
    this%papi_real_time_stop = 0.
    this%papi_proc_time_start = 0.
    this%papi_proc_time_stop = 0.
    this%papi_flpops_start = 0
    this%papi_flpops_stop = 0
    ! do a first call to initialize PAPI...
!    call PAPIF_flops(this%papi_real_time_start, this%papi_proc_time_start, this%papi_flpops_start &
!&                   , papi_mflops, papi_check)
    ! PAPIF_flops calls of multiple concurrent stopwatches interfere...
    ! better use papif_library_init, papif_start_counters...?
#endif

  end subroutine TStopwatch_Construct



!!==============================================================!
!!  Subroutine TStopwatch_Destruct                              !
!!==============================================================!
!
!  !> Destructor (empty! not needed)
!  !> \param this       ... object  TStopwatch
!  subroutine TStopwatch_Destruct( this )
!
!    implicit none
!
!    ! Declare arguments
!    type(TStopwatch) :: this
!
!  end subroutine TStopwatch_Destruct



!!==============================================================!
!!  Function TStopwatch_Tag                                     !
!!==============================================================!
!
!  !> Set/get tag
!  !> \param this     ... object           TStopwatch
!  !> \param new_tag  ... new tag to set   string
!  !> \return act_tag ... tag already set  string
!  function TStopwatch_Tag( this, new_tag ) result ( act_tag )
!
!    implicit none
!
!    ! Declare arguments
!    type(TStopwatch) :: this
!    character(*), intent(in), optional :: new_tag
!    character(tag_string_length) :: act_tag
!
!    act_tag = trim(this%tag_string)
!    if( present( new_tag ) ) then
!      this%tag_string = trim( new_tag )
!    end if
!
!  end function TStopwatch_Tag

!==============================================================!
!  Function TStopwatch_GetTag                                  !
!==============================================================!

  !> Get tag
  !> \param this        ... object      TStopwatch
  !> \return tag_string ... actual tag  string
  pure function TStopwatch_GetTag( this ) result ( tag_string )

    implicit none

    ! Declare arguments
    type(TStopwatch), intent(in) :: this

    character(tag_string_length) :: tag_string

    tag_string = trim( this%tag_string )

  end function TStopwatch_GetTag

!==============================================================!
!  Function TStopwatch_SetTag                                  !
!==============================================================!

  !> Set tag
  !> \param this     ... object          TStopwatch
  !> \param new_tag  ... new tag to set  string
  subroutine TStopwatch_SetTag( this, new_tag )

    implicit none

    ! Declare arguments
    type(TStopwatch) :: this
    character(*), intent(in) :: new_tag

    this%tag_string = trim( new_tag )

  end subroutine TStopwatch_SetTag



!==============================================================!
!  Function TStopwatch_GetOptions                              !
!==============================================================!

  !> Get options
  !> \param this        ... object       TStopwatch
  !> \return options ... actual options  integer
  pure function TStopwatch_GetOptions( this ) result ( options )

    implicit none

    ! Declare arguments
    type(TStopwatch), intent(in) :: this

    integer :: options

    options = this%options

  end function TStopwatch_GetOptions

!==============================================================!
!  Function TStopwatch_SetOptions                              !
!==============================================================!

  !> Set options
  !> \param this         ... object                     TStopwatch
  !> \param new_options  ... new options bitset to set  integer
  subroutine TStopwatch_SetOptions( this, new_options )

    implicit none

    ! Declare arguments
    type(TStopwatch) :: this
    integer, intent(in) :: new_options

    this%options = new_options

#ifndef STOPWATCH_USE_DATIME
    this%options = IOR( this%options, CStopwatch_omitDATIME )
#endif
#ifndef STOPWATCH_USE_CPUTIME
    this%options = IOR( this%options, CStopwatch_omitCPUTIME )
#endif
#ifndef STOPWATCH_USE_SYSCLK
    this%options = IOR( this%options, CStopwatch_omitSYSCLK )
#endif
#ifndef STOPWATCH_USE_ETIME
    this%options = IOR( this%options, CStopwatch_omitETIME )
#endif
#ifndef STOPWATCH_USE_PAPI
    this%options = IOR( this%options, CStopwatch_omitPAPI )
#endif
#ifndef STOPWATCH_USE_MPIWTIME
    this%options = IOR( this%options, CStopwatch_omitMPIWTIME )
#endif

  end subroutine TStopwatch_SetOptions

!==============================================================!
!  Function TStopwatch_ActivateOptions                         !
!==============================================================!

  !> Activate options
  !> \param this         ... object                                   TStopwatch
  !> \param act_options  ... options bitset to activate (bitwise OR)  integer
  subroutine TStopwatch_ActivateOptions( this, act_options )

    implicit none

    ! Declare arguments
    type(TStopwatch) :: this
    integer, intent(in) :: act_options

    this%options = IOR( this%options, act_options )

  end subroutine TStopwatch_ActivateOptions

!==============================================================!
!  Function TStopwatch_DeactivateOptions                       !
!==============================================================!

  !> Deactivate options
  !> \param this           ... object                                     TStopwatch
  !> \param deact_options  ... options bitset to deactivate (bitwise OR)  integer
  subroutine TStopwatch_DeactivateOptions( this, deact_options )

    implicit none

    ! Declare arguments
    type(TStopwatch) :: this
    integer, intent(in) :: deact_options

    this%options = IAND( this%options, NOT( deact_options ) )

  end subroutine TStopwatch_DeactivateOptions


#ifdef STOPWATCH_USE_DATIME
!==============================================================!
!  Function TStopwatch_GetDatimeStart                          !
!==============================================================!

  !> Get starting date and time
  !> \param this          ... object           TStopwatch
  !> \return datime_start ... start date&time  integer(8)
  pure function TStopwatch_GetDatimeStart( this ) result ( datime_start )

    implicit none

    ! Declare arguments
    type(TStopwatch), intent(in) :: this

    integer, dimension(8) :: datime_start

    datime_start = this%datime_array_start

  end function TStopwatch_GetDatimeStart

!==============================================================!
!  Function TStopwatch_GetDatimeStop                           !
!==============================================================!

  !> Get stopping date and time
  !> \param this          ... object         TStopwatch
  !> \return datime_stop ... stop date&time  integer(8)
  pure function TStopwatch_GetDatimeStop( this ) result ( datime_stop )

    implicit none

    ! Declare arguments
    type(TStopwatch), intent(in) :: this

    integer, dimension(8) :: datime_stop

    datime_stop = this%datime_array_stop

  end function TStopwatch_GetDatimeStop
#endif



#ifdef STOPWATCH_USE_SYSCLK
!==============================================================!
!  Function TStopwatch_GetSysClkDiff                           !
!==============================================================!

  !> Get system clock time difference
  !> \param this         ... object TStopwatch
  !> \return sysclk_cnt_diff ... system clock difference integer
  pure function TStopwatch_GetSysClkDiff( this ) result ( sysclk_cnt_diff )

    implicit none

    ! Declare arguments
    type(TStopwatch), intent(in) :: this

    integer :: sysclk_cnt_diff

    sysclk_cnt_diff = this%sysclk_cnt_stop-this%sysclk_cnt_start

  end function TStopwatch_GetSysClkDiff
#endif

#ifdef USE_MPI
!==============================================================!
!  Function TStopwatch_SetMPIcommunicator                      !
!==============================================================!

  !> Set MPI communicator
  !> \param this     ... object          TStopwatch
  !> \param new_mpicommunicator  ... new MPI communicator to set     integer
  subroutine TStopwatch_SetMPIcommunicator( this, new_mpicommunicator )

    implicit none

    ! Declare arguments
    type(TStopwatch) :: this
    integer, intent(in) :: new_mpicommunicator

    this%mpi_communicator = new_mpicommunicator

  end subroutine TStopwatch_SetMPIcommunicator
#endif



!==============================================================!
!  Subroutine TStopwatch_Start                                 !
!==============================================================!

  !> start the stopwatch timer
  !> \param this     ... object  TStopwatch
  subroutine TStopwatch_Start( this )

    implicit none

#if defined(USE_MPI) && !defined(MPI_USE_MODULE)
    include 'mpif.h'
#endif
#if defined STOPWATCH_USE_ETIME && defined __PGI
  include 'lib3f.h'
#endif
#ifdef STOPWATCH_USE_PAPI
  include "f90papi.h"
#endif

    ! Declare arguments
    type(TStopwatch) :: this

    ! Declare local variables
#ifdef STOPWATCH_USE_PAPI
    integer papi_check
    real(kind=4) papi_mflops
#endif

#ifdef STOPWATCH_USE_DATIME
    if (IAND(this%options,CStopwatch_omitDATIME) == 0) &
&     call date_and_time( values=this%datime_array_start )
#endif

#ifdef STOPWATCH_USE_SYSCLK
    if (IAND(this%options,CStopwatch_omitSYSCLK) == 0) &
&     call system_clock( this%sysclk_cnt_start )
#endif
#ifdef STOPWATCH_USE_MPIWTIME
    if (IAND(this%options,CStopwatch_omitMPIWTIME) == 0) then
      !if (BTEST(this%options,1)) then
      if (IAND(this%options,CStopwatch_doMPIStartBarrier) /= 0) then
          call MPI_Barrier( Communicator, ierror )
      end if
      this%mpi_diff_reduced=.FALSE.
      this%wtime_start = MPI_WTIME()
    endif
#endif
#ifdef STOPWATCH_USE_CPUTIME
    if (IAND(this%options,CStopwatch_omitCPUTIME) == 0) &
&     call cpu_time( this%cputime_start )
#endif
#ifdef STOPWATCH_USE_ETIME
    if (IAND(this%options,CStopwatch_omitDATIME) == 0) then
      !call etime(this%etime_array_start, this%etime_sum_start)  ! not supported by PGF
      this%etime_sum_start = etime( this%etime_array_start )
    endif
#endif
#ifdef STOPWATCH_USE_PAPI
    if (IAND(this%options,CStopwatch_omitDATIME) == 0) then
      call PAPIF_flops(this%papi_real_time_start, this%papi_proc_time_start, this%papi_flpops_start &
&                     , papi_mflops, papi_check)
    endif
#endif

  end subroutine TStopwatch_Start



!==============================================================!
!  Subroutine TStopwatch_Stop                                  !
!==============================================================!

  !> stop the stopwatch timer (multiple stops are possible)
  !> \param this     ... object  TStopwatch
  subroutine TStopwatch_Stop( this )

    implicit none

    ! Include headers
#if defined(USE_MPI) && !defined(MPI_USE_MODULE)
    include 'mpif.h'
#endif
#if defined STOPWATCH_USE_ETIME && defined __PGI
  include 'lib3f.h'
#endif
#ifdef STOPWATCH_USE_PAPI
  include "f90papi.h"
#endif

    ! Declare arguments
    type(TStopwatch) :: this

    ! Declare local variables
#ifdef STOPWATCH_USE_PAPI
    integer papi_check
    real(kind=4) papi_mflops
#endif

#ifdef STOPWATCH_USE_PAPI
    if (IAND(this%options,CStopwatch_omitPAPI) == 0) then
      ! do PAPIF_flops calls of multiple concurrent stopwatches interfere?
      call PAPIF_flops(this%papi_real_time_stop, this%papi_proc_time_stop, this%papi_flpops_stop &
&                     , papi_mflops, papi_check)
    endif
#endif
#ifdef STOPWATCH_USE_ETIME
    if (IAND(this%options,CStopwatch_omitETIME) == 0) &
&     this%etime_sum_stop = etime( this%etime_array_stop )
#endif
#ifdef STOPWATCH_USE_CPUTIME
    if (IAND(this%options,CStopwatch_omitCPUTIME) == 0) &
&     call cpu_time( this%cputime_stop )
#endif
#ifdef STOPWATCH_USE_SYSCLK
    if (IAND(this%options,CStopwatch_omitSYSCLK) == 0) &
&     call system_clock( this%sysclk_cnt_stop )
#endif
#ifdef STOPWATCH_USE_MPIWTIME
    if (IAND(this%options,CStopwatch_omitMPIWTIME) == 0) then
      !if (BTEST(this%options,2)) then
      if (IAND(this%options,CStopwatch_doMPIStopBarrier) /= 0) then
         call MPI_Barrier( Communicator, ierror )
      end if
      this%wtime_stop = MPI_WTIME()
      this%wtime_diff(1)=this%wtime_stop-this%wtime_start
      this%wtime_diff(2)=0.
      this%mpi_diff_reduced=.FALSE.
    endif
#endif
#ifdef STOPWATCH_USE_DATIME
    if (IAND(this%options,CStopwatch_omitDATIME) == 0) &
&     call date_and_time( values=this%datime_array_stop )
#endif

  end subroutine TStopwatch_Stop



#ifdef STOPWATCH_USE_MPIWTIME
!==============================================================!
!  Subroutine TStopwatch_MPIReduceDiff                         !
!==============================================================!

  !> reduce the wtime differences to get a max and min
  !> \param this     ... object  TStopwatch
  subroutine TStopwatch_MPIReduceDiff( this )

    implicit none

    ! Include headers
#if defined(USE_MPI) && !defined(MPI_USE_MODULE)
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TStopwatch) :: this

    ! Declare local variables
    double precision :: wtime_diff(2)

    !wtime_diff(1) = this%wtime_diff(1)
    wtime_diff(1) = this%wtime_stop-this%wtime_start
    wtime_diff(2) = -wtime_diff(1)

    !if (BTEST(this%options,3)) then
    if (IAND(this%options,CStopwatch_doMPIReduce) /= 0) then
       call MPI_Reduce( wtime_diff, this%wtime_diff, 2, MPI_DOUBLE_PRECISION, MPI_MAX &
&                     , NRootProc, Communicator, ierror )
       this%wtime_diff(2)=-this%wtime_diff(2)
    else
       this%wtime_diff(1)=wtime_diff(1)
       this%wtime_diff(2)=0.
    end if
    this%mpi_diff_reduced=.TRUE.

  end subroutine TStopwatch_MPIReduceDiff
#endif



!==============================================================!
!  Subroutine TStopwatch_LogWriteStart                         !
!==============================================================!
! TODO: argument to set options for writing short/detailed information...

  !> write start information to log; this means additional (I)O within the tested code part!
  !> \param this     ... object  TStopwatch
  subroutine TStopwatch_LogWriteStart( this )

    implicit none

    ! Declare arguments
    type(TStopwatch), intent(in) :: this

    call LogWriteBlank
    write( IOBuffer, '(72(1H=))')
    call LogWrite
#ifdef STOPWATCH_USE_DATIME
    if (IAND(this%options,CStopwatch_omitDATIME) == 0) then
      write( IOBuffer, &
&       '(A," timer start: ",I2.2,".",I2.2,".",I4.4,  &
&         ", ",I2.2,":",I2.2,":",I2.2,".",I3.3," (UTC",SP,I4,")")' ) &
&       trim(this%tag_string), &
&       this%datime_array_start(3), this%datime_array_start(2), this%datime_array_start(1), &
&       this%datime_array_start(5), this%datime_array_start(6), this%datime_array_start(7), &
&       this%datime_array_start(8), this%datime_array_start(4)
      call LogWrite
    end if
#endif
#ifdef STOPWATCH_USE_CPUTIME
    if (IAND(this%options,CStopwatch_omitCPUTIME) == 0) then
      write( IOBuffer, &
&       '(T2,A," cpu_time start:",G16.9," sec")' ) &
&       trim(this%tag_string), &
&       this%cputime_start
      call LogWrite
    end if
#endif
#ifdef STOPWATCH_USE_ETIME
    if (IAND(this%options,CStopwatch_omitETIME) == 0) then
      write( IOBuffer, &
&       '(T2,A," etime start:",F12.4,"(user) +",F12.4,"(system) =",G16.9," sec")' ) &
&       trim(this%tag_string), &
&       this%etime_array_start, this%etime_sum_start
      call LogWrite
    end if
#endif
    write( IOBuffer, '(72(1H-))')
    call LogWrite

  end subroutine TStopwatch_LogWriteStart



!==============================================================!
!  Subroutine TStopwatch_LogWriteStop                          !
!==============================================================!
! TODO: argument to set options for writing short/detailed...

  !> write stop information to log
  !> \param this     ... object  TStopwatch
  subroutine TStopwatch_LogWriteStop( this )

    implicit none

    ! Include MPI header
#if defined(USE_MPI) && !defined(MPI_USE_MODULE)
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TStopwatch) :: this

    ! Declare local variables
#ifdef STOPWATCH_USE_SYSCLK
    integer :: sysclk_cnt_rate, sysclk_cnt_max
    integer :: sysclk_cnt_diff
    real :: sysclk_diff_sec, sysclk_max_sec
#ifdef STOPWATCH_USE_MPIWTIME
    integer :: i
#endif
#endif
#ifdef STOPWATCH_USE_CPUTIME
    real :: cputime_diff
#endif
#ifdef STOPWATCH_USE_ETIME
    real, dimension(2) :: etime_array_diff
    real :: etime_sum_diff
#endif
#ifdef STOPWATCH_USE_PAPI
    real(kind=4) papi_mflops
    real(kind=4) papi_real_time_diff, papi_proc_time_diff, papi_flpops_diff
#endif

    write( IOBuffer, '(72(1H-))')
    call LogWrite
    write( IOBuffer, '("Duration: ")' )
    call LogWrite

#ifdef STOPWATCH_USE_SYSCLK
    if (IAND(this%options,CStopwatch_omitSYSCLK) == 0) then
      call system_clock(count_max=sysclk_cnt_max, count_rate=sysclk_cnt_rate)
      sysclk_cnt_diff = this%sysclk_cnt_stop-this%sysclk_cnt_start
      !sysclk_cnt_diff=mod(sysclk_cnt_diff+sysclk_cnt_max, sysclk_cnt_max) !sum needs long integer
      if( sysclk_cnt_diff<0 ) sysclk_cnt_diff = sysclk_cnt_diff+sysclk_cnt_max
      sysclk_diff_sec = float(sysclk_cnt_diff)/float(sysclk_cnt_rate)
      sysclk_max_sec = float(sysclk_cnt_max)/float(sysclk_cnt_rate)
#ifdef STOPWATCH_USE_MPIWTIME
      if (IAND(this%options,CStopwatch_omitMPIWTIME) == 0) then
        i = NINT((this%wtime_diff(1)-sysclk_diff_sec)/sysclk_max_sec)
        sysclk_diff_sec = sysclk_diff_sec+i*sysclk_max_sec
      end if
#endif
      write( IOBuffer, &
&       '(T2,A," system_clock diff:", G16.9,"sec")' ) &
&       trim(this%tag_string), sysclk_diff_sec
      call LogWrite
      write( IOBuffer, &
&       '(T18,"=",I5,"h",I3,"min",F9.5,"sec (+i*",I5,"h",I3,"min",F9.5,"sec)")' ) &
&       int(sysclk_diff_sec)/3600, mod(int(sysclk_diff_sec),3600)/60, amod(sysclk_diff_sec,60.), &
&       int(sysclk_max_sec)/3600, mod(int(sysclk_max_sec),3600)/60, amod(sysclk_max_sec,60.)
      call LogWrite
    end if
#endif

#ifdef STOPWATCH_USE_ETIME
    if (IAND(this%options,CStopwatch_omitETIME) == 0) then
      etime_array_diff = this%etime_array_stop-this%etime_array_start
      etime_sum_diff=this%etime_sum_stop-this%etime_sum_start
      write( IOBuffer, &
&       '(T2,A," etime        diff:",F12.4,"(user) +",F12.4,"(system) =",G16.9," sec")' ) &
&       trim(this%tag_string), &
&       etime_array_diff, etime_sum_diff
      call LogWrite
    end if
#endif
#ifdef STOPWATCH_USE_CPUTIME
    if (IAND(this%options,CStopwatch_omitCPUTIME) == 0) then
      cputime_diff = this%cputime_stop-this%cputime_start
      write( IOBuffer, &
&       '(T2,A," cpu_time     diff:",G16.9," sec")' ) &
&       trim(this%tag_string), &
&       cputime_diff
!    write( IOBuffer, &
!&     '(A," cpu_time     diff:",I5," h",I3," min",F9.5," sec =",G16.9," sec")' ) &
!&     trim(this%tag_string), &
!&     int(this%cputime_diff)/3600, mod(int(cputime_diff),3600)/60, dmod(cputime_diff,60.), &
!&     cputime_diff
      call LogWrite
    end if
#endif
#ifdef STOPWATCH_USE_MPIWTIME
    if (IAND(this%options,CStopwatch_omitMPIWTIME) == 0) then
      if (.NOT. this%mpi_diff_reduced) then
          call TStopwatch_MPIReduceDiff( this )
      end if
      if (IAND(this%options,CStopwatch_doMPIReduce) /= 0) then
        ! min/max reduction available
        write( IOBuffer, &
&         '(T2,A," wtime        diff:",G16.9,"-",G16.9)' ) &
&         trim(this%tag_string), this%wtime_diff(2), this%wtime_diff(1)
          call LogWrite
        write( IOBuffer,'(T31,"<=",I5,"h",I3,"min",F9.5,"sec (+-",E8.2,"sec)")' ) &
&         int(this%wtime_diff(1))/3600, mod(int(this%wtime_diff(1)),3600)/60, dmod(this%wtime_diff(1),60.D0), &
&         MPI_WTICK()
      else
        ! no min/max reduction available
        write( IOBuffer, &
&         '(T2,A," wtime   root diff:",G16.9)' ) &
&         trim(this%tag_string), this%wtime_diff(1)
          call LogWrite
        write( IOBuffer,'(T32,"=",I5,"h",I3,"min",F9.5,"sec (+-",E8.2,"sec)")' ) &
&         int(this%wtime_diff(1))/3600, mod(int(this%wtime_diff(1)),3600)/60, dmod(this%wtime_diff(1),60.D0), &
&         MPI_WTICK()
      end if
      call LogWrite
    endif
#endif
#ifdef STOPWATCH_USE_PAPI
    if (IAND(this%options,CStopwatch_omitPAPI) == 0) then
      papi_real_time_diff = this%papi_real_time_stop - this%papi_real_time_start
      papi_proc_time_diff = this%papi_proc_time_stop - this%papi_proc_time_start
      papi_flpops_diff = this%papi_flpops_stop - this%papi_flpops_start
      papi_mflops = papi_flpops_diff/1.0E6/papi_proc_time_diff
      write( IOBuffer, &
&       '(T2,A," PAPI             : real time",G16.9,"sec, proc time",G16.9,"sec ;", &
&         I16," FlOps,",G16.9," MFlOps/sec")' ) &
&       trim(this%tag_string), &
&       papi_real_time_diff, papi_proc_time_diff, &
&       papi_flpops_diff, papi_mflops
      call LogWrite
    end if
#endif
#ifdef STOPWATCH_USE_DATIME
    if (IAND(this%options,CStopwatch_omitDATIME) == 0) then
      write( IOBuffer, &
&       '(A," timer stop: ", I2.2,".",I2.2,".",I4.4,  &
&         ", ",I2.2,":",I2.2,":",I2.2,".",I3.3," (UTC",SP,I4,")")' ) &
&       trim(this%tag_string), &
&       this%datime_array_stop(3), this%datime_array_stop(2), this%datime_array_stop(1), &
&       this%datime_array_stop(5), this%datime_array_stop(6), this%datime_array_stop(7), &
&       this%datime_array_stop(8), this%datime_array_stop(4)
      call LogWrite
    end if
#endif
! #ifdef STOPWATCH_USE_CPUTIME
!     write( IOBuffer, &
! &     '(T22,"Cpu_time:",T35,G16.9," sec")' ) &
! &     this%cputime_stop
!     call LogWrite
! #endif
! #ifdef STOPWATCH_USE_ETIME
!     write( IOBuffer, &
! &     '("Etime start:",F12.4,"(user) +",F12.4,"(system) =",G16.9," sec")' ) &
! &     this%etime_array_stop, this%etime_sum_stop
!     call LogWrite
! #endif

    write( IOBuffer, '(72(1H=))')
    call LogWrite
    call LogWriteBlank

  end subroutine TStopwatch_LogWriteStop


end module ms2_stopwatch
