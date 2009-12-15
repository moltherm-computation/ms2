!==============================================================!
!> Module: ms2_profiler                                        !
!          contains TProfiler class                            !
!> \author Hendrik Adorf (ITWM)                                !
!> \date   March 2009                                          !
!> Comment: Timer is awkward; do not use across months (!)     !
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

!==============================================================!
!  Type TProfiler                                              !
!==============================================================!

  type TProfiler

    real(4) :: StartTime
    real(4) :: PreviousTimeTag
    real(4) :: OverallCommunicationTime
  
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
    real(4)         :: TimeTag
    character(8)    :: dummyDate
    character(10)   :: dummyTime
    character(5)    :: dummyZone
    integer         :: time(8)
    real(4)         :: millisec
    real(4)         :: timeSec
 
    this%iounitTrace   = iounitNumberTrace
    this%iounitRuntime = iounitNumberRuntime
    
    open(unit = this%iounitRuntime, file = trim(RuntimeFile), &
&      action = 'readwrite', status = 'replace', position = 'append')
    open(unit = this%iounitTrace, file = trim(TraceFile), &
&      action = 'write', status = 'replace')
    
    this%OverallCommunicationTime = 0.0
    
    call date_and_time(dummyDate, dummyTime, dummyZone, time)
    millisec = time(8)
    timeSec = time(3)*3600*24 + time(5)*3600 + &
&     time(6)*60 + time(7) + millisec/1000 

    this%StartTime = timeSec
    TimeTag = timeSec - this%StartTime
    
    write(this%iounitTrace, '(F10.3, A)') TimeTag, &
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
    real(4)         :: TimeTag
    character(8)    :: dummyDate
    character(10)   :: dummyTime
    character(5)    :: dummyZone
    integer         :: time(8)
    real(4)         :: millisec
    real(4)         :: timeSec

    call date_and_time(dummyDate, dummyTime, dummyZone, time)
    millisec = time(8)
    timeSec = time(3)*3600*24 + time(5)*3600 + &
&     time(6)*60 + time(7) + millisec/1000

    TimeTag = timeSec - this%StartTime
    this%PreviousTimeTag = TimeTag

    write(this%iounitTrace, '(F10.3, A)') TimeTag, &
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
    real(4)         :: TimeTag
    character(8)    :: dummyDate
    character(10)   :: dummyTime
    character(5)    :: dummyZone
    integer         :: time(8)
    real(4)         :: millisec
    real(4)         :: timeSec

    call date_and_time(dummyDate, dummyTime, dummyZone, time)
    millisec = time(8)
    timeSec = time(3)*3600*24 + time(5)*3600 + &
&     time(6)*60 + time(7) + millisec/1000

    TimeTag = timeSec - this%StartTime
    this%OverallCommunicationTime = this%OverallCommunicationTime + &
&     (TimeTag - this%PreviousTimeTag)

    write(this%iounitRuntime, '(F10.3, A, F10.3, A, F6.2)') &
&     TimeTag, ' ', this%OverallCommunicationTime, ' ', &
&     (this%OverallCommunicationTime/TimeTag)*100
!    backspace this%iounitRuntime ! why does it not work?

    write(this%iounitTrace, '(F10.3, A)') TimeTag, &
&     trim(' COMPLETED '//trim(tagStr) )

  end subroutine TProfiler_After

end module ms2_profiler
