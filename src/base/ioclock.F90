!Generic output clock for various models
module ioclock
    use kdefs
    use xml_input
    implicit none

    type IOClock_T
        !Output/restart (cadence/next output)
        real(rp) :: dtOut,tOut
        real(rp) :: dtRes,tRes
        integer :: tsOut,tsNext !Timestep console output
        integer :: nOut=0,nRes=0 !Output numbering for output/restarts

        logical :: doResOut,doConOut,doDataOut,doTimerOut !Logical flags to do various outputs
        logical :: doFat,doSlim !Not yet used generically
        logical :: isTop = .true. !Not used yet

        contains

            procedure :: init => IOClockInit
            procedure :: doConsole => doConsoleIOClock
            procedure :: doOutput  => doOutputIOClock
            procedure :: doRestart => doRestartIOClock
            procedure :: doTimer => doTimerIOClock

    end type IOClock_T

    contains

    !Initialize IOClock from an already created XML reader
    subroutine IOClockInit(this,iXML,time,ts)
        class(IOClock_T), intent(inout) :: this
        type(XML_Input_T), intent(in)   :: iXML
        real(rp), intent(in) :: time
        integer, intent(in) :: ts

        call iXML%Set_Val(this%tsOut,'output/tsOut' ,10)
        call iXML%Set_Val(this%dtOut,'output/dtOut' ,10.0)
        call iXML%Set_Val(this%dtRes,'restart/dtRes',100.0)

        !Get booleans
        call iXML%Set_Val(this%doTimerOut,   'output/doTimer'   ,.false.)
        call iXML%Set_Val(this%doFat,        'output/doFat'     ,.false.)
        call iXML%Set_Val(this%doSlim,       'output/doSlim'    ,.false.)

        if (this%tsOut<0) then
            this%doConOut = .false.
        else
            this%doConOut = .true.
            this%tsNext = ts
        endif

        !NOTE: Setting it so that first output is @ time
        if (this%dtOut<0) then
            this%doDataOut = .false.
            this%tOut = HUGE
        else
            this%doDataOut = .true.
            this%tOut = time
        endif

        if (this%dtRes<0) then
            this%doResOut = .false.
            this%tRes = HUGE
        else
            this%doResOut = .true.
            this%tRes = time
        endif

    end subroutine IOClockInit
    
    !Check for various output times
    function doConsoleIOClock(this,ts)
        class(IOClock_T), intent(in) :: this
        integer, intent(in) :: ts
        logical :: doConsoleIOClock

        if (ts>=this%tsNext .and. this%doConOut ) then
            doConsoleIOClock = .true.
        else
            doConsoleIOClock = .false.
        endif

    end function doConsoleIOClock

    function doRestartIOClock(this,t)
        class(IOClock_T), intent(in) :: this
        real(rp), intent(in) :: t
        logical :: doRestartIOClock

        if ( (t>=this%tRes) .and. this%doResOut ) then
            doRestartIOClock = .true.
        else
            doRestartIOClock = .false.
        endif
    end function doRestartIOClock

    function doOutputIOClock(this,t)
        class(IOClock_T), intent(in) :: this
        real(rp), intent(in) :: t
        logical :: doOutputIOClock

        if ( (t>=this%tOut) .and. this%doDataOut ) then
            doOutputIOClock = .true.
        else
            doOutputIOClock = .false.
        endif
    end function doOutputIOClock

    !Only check for tsOut because we need to clean clocks even if not printing
    function doTimerIOClock(this,ts)
        class(IOClock_T), intent(in) :: this
        integer, intent(in) :: ts
        logical :: doTimerIOClock

        if (ts>=this%tsNext) then
            doTimerIOClock = .true.
        else
            doTimerIOClock = .false.
        endif

    end function doTimerIOClock

    !Copy ioA=>ioB using tScl scaling
    subroutine IOSync(ioA,ioB,tSclO)
        class(IOClock_T), intent(in   ) :: ioA
        class(IOClock_T), intent(inout) :: ioB
        real(rp), intent(in), optional :: tSclO

        real(rp) :: tScl

        if (present(tSclO)) then
            tScl = tSclO
        else
            tScl = 1.0
        endif

        ioB%tOut = ioA%tOut*tScl
        ioB%tRes = ioA%tRes*tScl

        ioB%dtOut = ioA%dtOut*tScl
        ioB%dtRes = ioA%dtRes*tScl

        ioB%tsOut  = ioA%tsOut
        ioB%tsNext = ioA%tsNext
        
    end subroutine IOSync

end module ioclock
