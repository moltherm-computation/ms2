!==============================================================!
!  MOLECULAR SIMULATION PROGRAM MS2 Version 1.1 v12            !
!==============================================================!
!> Module: ms2_stopwatch                                       !
!          contains TStopwatch class                           !
!> \author Martin Bernreuther <bernreuther@hlrs.de>            !
!> \date   06.2009                                             !
!==============================================================!
! better use http://math.nist.gov/StopWatch/?                  !
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

#if MPI_VER
#define STOPWATCH_USE_WTIME
#else
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
#endif

module ms2_stopwatch

  use ms2_global
  !use mpi

#if defined STOPWATCH_USE_ETIME && defined __INTEL_COMPILER
  use IFPORT
! __INTEL_COMPILER also offers DTIME and DCLOCK via IFPORT
#endif


  integer, parameter :: tag_string_length = 64


  !integer, parameter :: CStopwatch_useStartBarrier = B'1'
  integer, parameter :: CStopwatch_useStartBarrier = 1
  !integer, parameter :: CStopwatch_useStopBarrier = B'10'
  integer, parameter :: CStopwatch_useStopBarrier = 2
  !integer, parameter :: CStopwatch_doReduceMax = B'100'
  integer, parameter :: CStopwatch_doReduceMax = 4
  !integer, parameter :: CStopwatch_doReduceMin = B'1000'
  integer, parameter :: CStopwatch_doReduceMin = 8


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
#ifdef STOPWATCH_USE_WTIME
    double precision :: wtime_start, wtime_stop
    double precision, dimension(2) :: wtime_diff
    logical wtime_diff_reduced
#endif
#if defined STOPWATCH_USE_CPUTIME
    real :: cputime_start, cputime_stop
#endif
#ifdef STOPWATCH_USE_ETIME
    real, dimension(2) :: etime_array_start, etime_array_stop
    real :: etime_sum_start, etime_sum_stop
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
  !> \param this       ... object	TStopwatch
  !> \param tag_string ... tag	string
  subroutine TStopwatch_Construct( this, tag_string, options )

    implicit none

    ! Declare arguments
    type(TStopwatch) :: this
    character(*), intent(in), optional :: tag_string
    integer, intent(in), optional :: options

    if( present( tag_string ) ) then
      this%tag_string = tag_string
    else
      this%tag_string = ""
    end if

    if( present( options ) ) then
      this%options = options
    else
      this%options = CStopwatch_doReduceMax + CStopwatch_doReduceMin
    end if

    !call TStopwatch_Start(this)

    this%datime_array_start = 0
    this%datime_array_stop = 0

#ifdef STOPWATCH_USE_SYSCLK
    this%sysclk_cnt_start = 0
    this%sysclk_cnt_stop = 0
#endif
#ifdef STOPWATCH_USE_WTIME
    this%wtime_start = 0
    this%wtime_stop = 0
    this%wtime_diff = 0
    this%wtime_diff_reduced=.FALSE.
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

  end subroutine TStopwatch_Construct



!!==============================================================!
!!  Subroutine TStopwatch_Destruct                              !
!!==============================================================!
!
!  !> Destructor (empty! not needed)
!  !> \param this       ... object	TStopwatch
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
!  !> \param this     ... object	TStopwatch
!  !> \param new_tag  ... new tag to set	string
!  !> \return act_tag ... tag already set	string
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
  !> \param this        ... object	TStopwatch
  !> \return tag_string ... actual tag	string
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
  !> \param this     ... object	TStopwatch
  !> \param new_tag  ... new tag to set	string
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
  !> \param this        ... object	TStopwatch
  !> \return options ... actual options	integer
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
  !> \param this         ... object	TStopwatch
  !> \param new_options  ... new options bitset to set	integer
  subroutine TStopwatch_SetOptions( this, new_options )

    implicit none

    ! Declare arguments
    type(TStopwatch) :: this
    integer, intent(in) :: new_options

    this%options = new_options

  end subroutine TStopwatch_SetOptions

!==============================================================!
!  Function TStopwatch_ActivateOptions                         !
!==============================================================!

  !> Activate options
  !> \param this         ... object	TStopwatch
  !> \param act_options  ... options bitset to activate (bitwise OR)	integer
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
  !> \param this           ... object	TStopwatch
  !> \param deact_options  ... options bitset to deactivate (bitwise OR)	integer
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
  !> \param this          ... object	TStopwatch
  !> \return datime_start ... start date&time	integer(8)
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
  !> \param this          ... object	TStopwatch
  !> \return datime_stop ... stop date&time	integer(8)
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
  !> \param this         ... object	TStopwatch
  !> \return sysclk_diff ... system clock difference	integer
  pure function TStopwatch_GetSysClkDiff( this ) result ( sysclk_diff )

    implicit none

    ! Declare arguments
    type(TStopwatch), intent(in) :: this

    integer :: sysclk_diff

    sysclk_diff = this%sysclk_cnt_stop-this%sysclk_cnt_start

  end function TStopwatch_GetSysClkDiff
#endif



!==============================================================!
!  Subroutine TStopwatch_Start                                 !
!==============================================================!

  !> start the stopwatch timer
  !> \param this     ... object	TStopwatch
  subroutine TStopwatch_Start( this )

    implicit none

    ! Include MPI header
#if MPI_VER
    include 'mpif.h'
#endif
#if defined STOPWATCH_USE_ETIME && defined __PGI
  include 'lib3f.h'
#endif

    ! Declare arguments
    type(TStopwatch) :: this

#ifdef STOPWATCH_USE_DATIME
    call date_and_time( values=this%datime_array_start )
#endif

#ifdef STOPWATCH_USE_SYSCLK
    call system_clock( this%sysclk_cnt_start )
#endif
#ifdef STOPWATCH_USE_WTIME
    !if (BTEST(this%options,1)) then
    if (IAND(this%options,CStopwatch_useStartBarrier) /= 0) then
        call MPI_Barrier( MPI_COMM_WORLD, ierror )
    end if
    this%wtime_diff_reduced=.FALSE.
    this%wtime_start = MPI_WTIME()
#endif
#ifdef STOPWATCH_USE_CPUTIME
    call cpu_time( this%cputime_start )
#endif
#ifdef STOPWATCH_USE_ETIME
    !call etime(this%etime_array_start, this%etime_sum_start)	! not supported by PGF
    this%etime_sum_start = etime( this%etime_array_start )
#endif

  end subroutine TStopwatch_Start



!==============================================================!
!  Subroutine TStopwatch_Stop                                  !
!==============================================================!

  !> stop the stopwatch timer (multiple stops are possible)
  !> \param this     ... object	TStopwatch
  subroutine TStopwatch_Stop( this )

    implicit none

    ! Include headers
#if MPI_VER
    include 'mpif.h'
#endif
#if defined STOPWATCH_USE_ETIME && defined __PGI
  include 'lib3f.h'
#endif

    ! Declare arguments
    type(TStopwatch) :: this

#ifdef STOPWATCH_USE_ETIME
    this%etime_sum_stop = etime( this%etime_array_stop )
#endif
#ifdef STOPWATCH_USE_CPUTIME
    call cpu_time( this%cputime_stop )
#endif
#ifdef STOPWATCH_USE_SYSCLK
    call system_clock( this%sysclk_cnt_stop )
#endif
#ifdef STOPWATCH_USE_WTIME
    !if (BTEST(this%options,2)) then
    if (IAND(this%options,CStopwatch_useStopBarrier) /= 0) then
       call MPI_Barrier( MPI_COMM_WORLD, ierror )
    end if
    this%wtime_stop = MPI_WTIME()
    this%wtime_diff(1)=this%wtime_stop-this%wtime_start
    this%wtime_diff(2)=0.
    this%wtime_diff_reduced=.FALSE.
#endif
#ifdef STOPWATCH_USE_DATIME
    call date_and_time( values=this%datime_array_stop )
#endif

  end subroutine TStopwatch_Stop



#ifdef STOPWATCH_USE_WTIME
!==============================================================!
!  Subroutine TStopwatch_ReduceWtimeDiff                       !
!==============================================================!

  !> reduce the wtime differences to get a max and min
  !> \param this     ... object	TStopwatch
  subroutine TStopwatch_ReduceWtimeDiff( this )

    implicit none

    ! Include headers
#if MPI_VER
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TStopwatch) :: this

    ! Declare local variables
    double precision :: wtime_diff

    !wtime_diff = this%wtime_diff(1)
    wtime_diff = this%wtime_stop-this%wtime_start
    !if (BTEST(this%options,3)) then
    if (IAND(this%options,CStopwatch_doReduceMax) /= 0) then
       call MPI_Reduce( wtime_diff, this%wtime_diff(1), 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
&                       NRootProc, MPI_COMM_WORLD, ierror )
    else
       this%wtime_diff(1)=wtime_diff
    end if
    !if (BTEST(this%options,4)) then
    if (IAND(this%options,CStopwatch_doReduceMin) /= 0) then
       call MPI_Reduce( wtime_diff, this%wtime_diff(2), 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
&                       NRootProc, MPI_COMM_WORLD, ierror )
    else
       this%wtime_diff(2)=0.
    end if
    this%wtime_diff_reduced=.TRUE.

  end subroutine TStopwatch_ReduceWtimeDiff
#endif



!==============================================================!
!  Subroutine TStopwatch_LogWriteStart                         !
!==============================================================!
! TODO: argument to set options for writing short/detailed...

  !> write start information to log; this means additional (I)O within the tested code part!
  !> \param this     ... object	TStopwatch
  subroutine TStopwatch_LogWriteStart( this )

    implicit none

    ! Declare arguments
    type(TStopwatch), intent(in) :: this

    call LogWriteBlank
#ifdef STOPWATCH_USE_DATIME
    write( IOBuffer, &
&     '("Timer ",A," start: ", &
&       I2.2,".",I2.2,".",I4.4,"T",I2.2,":",I2.2,":",I2.2,".",I3.3," (",SP,I4,")")' ) &
&     trim(this%tag_string), &
&     this%datime_array_start(3), this%datime_array_start(2), this%datime_array_start(1), &
&     this%datime_array_start(5), this%datime_array_start(6), this%datime_array_start(7), &
&     this%datime_array_start(8), this%datime_array_start(4)
    call LogWrite
#endif
#ifdef STOPWATCH_USE_CPUTIME
    write( IOBuffer, &
&     '(A," cpu_time start:",G15.9,"sec ")' ) &
&     trim(this%tag_string), &
&     this%cputime_start
    call LogWrite
#endif
#ifdef STOPWATCH_USE_ETIME
    write( IOBuffer, &
&     '(A," etime start:",F12.4,"(user) +",F12.4,"(system) = ",G15.9,"sec ")' ) &
&     trim(this%tag_string), &
&     this%etime_array_start, this%etime_sum_start
    call LogWrite
#endif

  end subroutine TStopwatch_LogWriteStart



!==============================================================!
!  Subroutine TStopwatch_LogWriteStop                          !
!==============================================================!
! TODO: argument to set options for writing short/detailed...

  !> write stop information to log
  !> \param this     ... object	TStopwatch
  subroutine TStopwatch_LogWriteStop( this )

    implicit none

    ! Include MPI header
#if MPI_VER
    include 'mpif.h'
#endif

    ! Declare arguments
    type(TStopwatch) :: this

    ! Declare local variables
#ifdef STOPWATCH_USE_SYSCLK
    integer :: sysclk_cnt_rate, sysclk_cnt_max
    integer :: sysclk_diff
    real :: sysclk_diff_sec
#endif
#ifdef STOPWATCH_USE_CPUTIME
    real :: cputime_diff
#endif
#ifdef STOPWATCH_USE_ETIME
    real, dimension(2) :: etime_array_diff
    real :: etime_sum_diff
#endif

#ifdef STOPWATCH_USE_DATIME
    write( IOBuffer, &
&     '("Timer ",A," stop: ", &
&       I2.2,".",I2.2,".",I4.4,"T",I2.2,":",I2.2,":",I2.2,".",I3.3," (",SP,I4,")")' ) &
&     trim(this%tag_string), &
&     this%datime_array_stop(3), this%datime_array_stop(2), this%datime_array_stop(1), &
&     this%datime_array_stop(5), this%datime_array_stop(6), this%datime_array_stop(7), &
&     this%datime_array_stop(8), this%datime_array_stop(4)
    call LogWrite
#endif

#ifdef STOPWATCH_USE_SYSCLK
    call system_clock(count_max=sysclk_cnt_max, count_rate=sysclk_cnt_rate)
    sysclk_diff = this%sysclk_cnt_stop-this%sysclk_cnt_start
    !sysclk_diff=mod(sysclk_diff+sysclk_cnt_max, sysclk_cnt_max) !sum needs long integer
    if( sysclk_diff<0 ) sysclk_diff = sysclk_diff+sysclk_cnt_max
    sysclk_diff_sec = float(sysclk_diff)/float(sysclk_cnt_rate)
    write( IOBuffer, &
&     '(A," system_clock diff:",I5," h",I3," min",F9.5, " sec (",G15.9,"+i*",E15.9," sec)")' ) &
&     trim(this%tag_string), &
&     int(sysclk_diff_sec)/3600, mod(int(sysclk_diff_sec),3600)/60, amod(sysclk_diff_sec,60.), &
&     sysclk_diff_sec, float(sysclk_cnt_max)/float(sysclk_cnt_rate)
    call LogWrite
#endif
#ifdef STOPWATCH_USE_WTIME
    if (.NOT. this%wtime_diff_reduced) then
        call TStopwatch_ReduceWtimeDiff( this )
    end if
    write( IOBuffer, &
&     '(A," wtime    max diff:",I5," h",I3," min",F9.5," sec =",G15.9," (min.",G15.9,") +-",E8.2," sec")' ) &
&     trim(this%tag_string), &
&     int(this%wtime_diff(1))/3600, mod(int(this%wtime_diff(1)),3600)/60, dmod(this%wtime_diff(1),60.), &
&     this%wtime_diff(1), this%wtime_diff(2), MPI_WTICK()
    call LogWrite
#endif
#ifdef STOPWATCH_USE_CPUTIME
    cputime_diff = this%cputime_stop-this%cputime_start
    write( IOBuffer, &
&     '(A," cpu_time     diff:",G15.9," sec ")' ) &
&     trim(this%tag_string), &
&     cputime_diff
    call LogWrite
#endif
#ifdef STOPWATCH_USE_ETIME
    etime_array_diff = this%etime_array_stop-this%etime_array_start
    !etime_sum_diff=this%etime_sum_stop-this%etime_sum_start
    write( IOBuffer, &
&     '(A," etime        diff:",F12.4,"(user) +",F12.4,"(system) = ",G15.9," sec ")' ) &
&     trim(this%tag_string), &
&     etime_array_diff, etime_sum_diff
    call LogWrite
#endif
    call LogWriteBlank

  end subroutine TStopwatch_LogWriteStop


end module ms2_stopwatch
