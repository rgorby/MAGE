module raijuOut
    use raijuIO

    implicit none

    contains

    subroutine raijuOutput(Model, Grid, State)
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        character(len=strLen) :: gStr,tStr

        write (gStr, '(A,I0)') "Step#", State%IO%nOut

        if (Model%isLoud) then
            !if (State%t>1.0e-2) then
            !    write(tStr,'(f12.3,a)' ) State%t, ' [code]'
            !else
            !    write(tStr,'(es12.2,a)') State%t, ' [code]'
            !endif
            call timeString(State%t, tStr)
            write (*, '(a,a,a,a,a)') ANSIGREEN, '<Writing RAIJU HDF5 DATA @ t = ', trim(tStr), ' >', ANSIRESET
        endif

        call WriteRAIJU(Model, Grid, State, gStr, Model%writeGhosts)

        !Setup for next output
        State%IO%nOut = State%IO%nOut + 1
        State%IO%tOut = State%t + State%IO%dtOut
    end subroutine raijuOutput


    subroutine raijuResOutput(Model, Grid, State)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        character(len=strLen) :: ResF, tStr,lnResF !Name of restart file
        logical :: fExist

        write (ResF, '(A,A,I0.5,A)') trim(Model%RunID), ".raiju.Res.", State%IO%nRes, ".h5"

        call CheckAndKill(ResF)

        if (Model%isLoud) then
            call timeString(State%t, tStr)
            write (*, '(a,a,a,a,a)') ANSIGREEN, '<Writing HDF5 RESTART @ t = ', trim(tStr), ' >', ANSIRESET
        endif

        !!!! Write here !!!!

        ! Prep for next restart
        Sate%IO%tRes = State%IO%tRes + State%IO%dtRes
        Sate%IO%nRes = State%IO%nRes + 1

        write (lnResF, '(A,A,A,A)') trim(Model%RunID), ".raiju.Res.", "XXXXX", ".h5"

        call MapSymLink(ResF,lnResF)

    end subroutine raijuResOutput


    subroutine raijuResInput(Model, Grid, State)
        !! Mirror of raijuResOutput to make sure state is restored to same exact state as above
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        character(len=strLen) :: ResF

        if (Model%)

    end subroutine raijuResInput

end module raijuOut