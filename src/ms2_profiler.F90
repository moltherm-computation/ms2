!==============================================================!
!> Module: ms2_profiler                                        !
!          contains TProfiler class                            !
!> \author Hendrik Adorf (ITWM)                                !
!> \date   March 2009                                          !
!> Comment: A Microsecond Timer                                !
!==============================================================!

#ifndef ARCH
#define ARCH    0
#define FORTRAN 90
#define MPI_VER 0
#define FVM_VER 0
#endif

#if ARCH == 1 || defined __INTEL_COMPILER
!DEC$ MESSAGE:'Compiling ms2_component.F90...'
#endif


module ms2_profiler

#if FVM_VER > 0
  use fvmf2003extensions
#endif

!==============================================================!
!  Type TProfiler                                              !
!==============================================================!

  type TProfiler

    real(8) :: StartTime
    real(8) :: PreviousTimeTag
    real(8) :: OverallCommunicationTime
  
    integer :: iounitRuntime
    integer :: iounitTrace
  
  end type TProfiler

  interface constructProfiler
    module procedure TProfiler_Construct
  end interface

  interface destructProfiler
    module procedure TProfiler_Destruct
  end interface

  interface profileTagBefore
    module procedure TProfiler_Before
  end interface

  interface profileTagAfter
    module procedure TProfiler_After
  end interface

contains

!==============================================================!
!  Subroutine TProfiler_Construct                              !
!==============================================================!

  subroutine TProfiler_Construct( this, TraceFile, RuntimeFile, &
&   iounitNumberTrace, iounitNumberRuntime )

    implicit none

    !arguments
    type(TProfiler) :: this
    character(*)  :: TraceFile
    character(*)  :: RuntimeFile
    integer         :: iounitNumberTrace
    integer         :: iounitNumberRuntime

    !local variables
    real(8)         :: timeTag, time
!    character(8)    :: dummyDate
!    character(10)   :: dummyTime
!    character(5)    :: dummyZone
!    integer         :: time(8)
!    real(4)         :: millisec
!    real(4)         :: timeSec
 
    this%iounitTrace   = iounitNumberTrace
    this%iounitRuntime = iounitNumberRuntime
    
    open(unit = this%iounitRuntime, file = trim(RuntimeFile), &
&      action = 'readwrite', status = 'replace', position = 'append')
    open(unit = this%iounitTrace, file = trim(TraceFile), &
&      action = 'write', status = 'replace')
    
    this%OverallCommunicationTime = 0.0
    
!    call date_and_time(dummyDate, dummyTime, dummyZone, time)
!    millisec = time(8)
!    timeSec = time(3)*3600*24 + time(5)*3600 + &
!&     time(6)*60 + time(7) + millisec/1000 
!
!    this%StartTime = timeSec
!    timeTag = timeSec - this%StartTime

    call vmTimer(time)
    this%StartTime = time
    timeTag = time - this%StartTime

    write(this%iounitTrace, '(F18.6, A)') timeTag, &
&     ' COMPLETED TProfiler_Construct'

  end subroutine TProfiler_Construct

!==============================================================!
!  Subroutine TProfiler_Destruct                               !
!==============================================================!

  subroutine TProfiler_Destruct( this )

    implicit none

    !arguments
    type(TProfiler) :: this

    close(this%iounitRuntime)
    close(this%iounitTrace)

  end subroutine TProfiler_Destruct

!==============================================================!
!  Subroutine TProfiler_Before                                 !
!==============================================================!

  subroutine TProfiler_Before( this, tagStr )

    implicit none

    !arguments
    type(TProfiler) :: this
    character(*)  :: tagStr

    !local variables
    real(8)         :: timeTag, time
!    character(8)    :: dummyDate
!    character(10)   :: dummyTime
!    character(5)    :: dummyZone
!    integer         :: time(8)
!    real(4)         :: millisec
!    real(4)         :: timeSec

!    call date_and_time(dummyDate, dummyTime, dummyZone, time)
!    millisec = time(8)
!    timeSec = time(3)*3600*24 + time(5)*3600 + &
!&     time(6)*60 + time(7) + millisec/1000
!
!    timeTag = timeSec - this%StartTime
!    this%PreviousTimeTag = timeTag

    call vmTimer(time)
    timeTag = time - this%StartTime
    this%PreviousTimeTag = timeTag

    write(this%iounitTrace, '(F18.6, A)') timeTag, &
&     trim(' REACHED '//trim(tagStr) )

  end subroutine TProfiler_Before

!==============================================================!
!  Subroutine TProfiler_After                                  !
!==============================================================!

  subroutine TProfiler_After( this, tagStr )

    implicit none

    !arguments
    type(TProfiler) :: this
    character(*)  :: tagStr

    !local variables
    real(8)         :: timeTag, time
!    character(8)    :: dummyDate
!    character(10)   :: dummyTime
!    character(5)    :: dummyZone
!    integer         :: time(8)
!    real(4)         :: millisec
!    real(4)         :: timeSec

!    call date_and_time(dummyDate, dummyTime, dummyZone, time)
!    millisec = time(8)
!    timeSec = time(3)*3600*24 + time(5)*3600 + &
!&     time(6)*60 + time(7) + millisec/1000
!
!    TimeTag = timeSec - this%StartTime
!    this%OverallCommunicationTime = this%OverallCommunicationTime + &
!&     (TimeTag - this%PreviousTimeTag)

    call vmTimer(time)
    timeTag = time - this%StartTime
    this%OverallCommunicationTime = this%OverallCommunicationTime + &
&     (timeTag - this%PreviousTimeTag)

    write(this%iounitRuntime, '(F18.6, A, F18.6, A, F6.2)') &
&     timeTag, ' ', this%OverallCommunicationTime, ' ', &
&     (this%OverallCommunicationTime/timeTag)*100
!    backspace this%iounitRuntime ! why does it not work?

    write(this%iounitTrace, '(F18.6, A)') timeTag, &
&     trim(' COMPLETED '//trim(tagStr) )

  end subroutine TProfiler_After

end module ms2_profiler
