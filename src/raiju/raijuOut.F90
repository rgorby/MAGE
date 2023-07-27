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
            if (State%t>1.0e-2) then
                write(tStr,'(f9.3,a)' ) State%t, ' [code]'
            else
                write(tStr,'(es9.2,a)') State%t, ' [code]'
            endif
            write (*, '(a,a,a,a,a)') ANSIGREEN, '<Writing RAIJU HDF5 DATA @ t = ', trim(tStr), ' >', ANSIRESET
        endif

        call WriteRAIJU(Model, Grid, State, gStr)

        !Setup for next output
        State%IO%nOut = State%IO%nOut + 1
        State%IO%tOut = State%t + State%IO%dtOut
    end subroutine raijuOutput

end module raijuOut