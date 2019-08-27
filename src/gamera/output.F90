module output
    use gdefs
    use gamtypes
    use clocks
    use gioH5
    use gridutils

    implicit none

    logical :: primOut = .true.
    integer :: nFloors = 0
    real(rp) :: dt0 = 0

    character :: TAB = char(9)
    character(len=strLen) :: zcsClk = "Gamera" !Clock ID to use for ZCS calculation

    !ConOut_T
    !Console output function pointer
    abstract interface
        subroutine ConsoleOut_T(Model,Grid,State)
            Import :: Model_T, Grid_T, State_T
            type(Model_T), intent(in) :: Model
            type(Grid_T), intent(in) :: Grid
            type(State_T), intent(in) :: State
        end subroutine ConsoleOut_T
    end interface

    !Set console output function here
    procedure(ConsoleOut_T), pointer :: consoleOutput => consoleOutput_STD

    !tStr_T
    !Time string function
    abstract interface
        subroutine tStr_T(T,tStr)
            Import :: rp,strLen
            real(rp), intent(in) :: T
            character(len=strLen), intent(out) :: tStr
        end subroutine tStr_T
    end interface

    procedure(tStr_T), pointer :: timeString => tStr_STD

contains

    subroutine consoleOutput_STD(Model, Grid, State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: State

        real(rp) :: ZCs, wTime
        integer :: nTh
        character(len=strLen) :: tStr
        if (dt0 < TINY) dt0 = Model%dt
        wTime = readClock(zcsClk)

        !Calculate zone-cycles per second
        if (Model%ts > 0) then
            ZCs = Model%tsOut*Grid%Nip*Grid%Njp*Grid%Nkp/wTime
        else
            ZCs = 0.0
        endif

        if (verbose > 0) then
            write(*,*) ANSICYAN
            call timeString(Model%t,tStr)
            write(*,'(a,a)')        'Sim Time   = ', trim(tStr)
            write (*, '(a,I8)')     '        ts = ', Model%ts
            call timeString(Model%dt,tStr)
            write (*, '(a,a)')      '        dt = ', trim(tStr)
            write (*, '(a,f8.3,a)') '    dt/dt0 = ', 100*Model%dt/dt0, '%'
            
            nTh = 0
#ifdef _OPENMP
            nTh = Model%nTh
#endif
            write (*, '(a,f9.2,a,I0,a)') '      kZCs = ', ZCs/1000.0, ' (',nTh,' threads)'
            if (nFloors > 0) then
                write (*, '(a,I8)') '      nFloors = ', nFloors
            endif
            write (*, *) ANSIRESET, ''

        endif

    end subroutine consoleOutput_STD

    subroutine tStr_STD(T,tStr)
        real(rp), intent(in) :: T
        character(len=strLen), intent(out) :: tStr

        if (T>1.0e-2) then
            write(tStr,'(f9.3,a)' ) T, ' [code]'
        else
            write(tStr,'(es9.2,a)') T, ' [code]'
        endif
    end subroutine tStr_STD

    subroutine fOutput(Model, Grid, State)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: State

        character(len=strLen) :: gStr,tStr

        write (gStr, '(A,I0)') "Step#", Model%nOut

        if (verbose > 0) then
            call timeString(Model%t,tStr)
            write (*, '(a,a,a,a,a)') ANSIPURPLE, '<Writing HDF5 output @ t = ', trim(tStr), ' >', ANSIRESET
        endif

        !call writeH5Slc(Model,Grid,State,gStr)
        call writeSlc(Model, Grid, State, gStr)

        !Setup for next output
        Model%tOut = Model%tOut + Model%dtOut
        Model%nOut = Model%nOut + 1
    end subroutine fOutput

    subroutine resOutput(Model, Grid, State)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: State

        character(len=strLen) :: ResF, lnResF !Name of restart file
        logical :: fExist

        write (ResF, '(A,A,I0.5,A)') trim(Model%RunID), ".Res.", Model%nRes, ".h5"

        if (verbose > 0) then
            write (*, *) '-----------------------'
            inquire (file=ResF, exist=fExist)
            if (fExist) then
                write (*, *) 'Restart file already exists, deleting file.'
                call EXECUTE_COMMAND_LINE('rm '//trim(ResF), wait=.true.)
            endif

            write (*, '(a,f8.3,a)') '<Writing HDF5 RESTART @ t = ', Model%t, ' >'
            write (*, *) '-----------------------'
        endif

        call writeH5Res(Model, Grid, State, ResF)

        !Setup for next restart
        Model%tRes = Model%tRes + Model%dtRes
        Model%nRes = Model%nRes + 1

        write (lnResF, '(A,A,A,A)') trim(Model%RunID), ".Res.", "XXXXX", ".h5"

        ! make a link to the default "XXXXX" restart file
        call EXECUTE_COMMAND_LINE('ln -sf '//trim(ResF)//' '//trim(lnResF), wait=.false.)

    end subroutine resOutput

end module output
