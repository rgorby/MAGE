module sifOut
    use sifIO

    implicit none

    contains

    subroutine sifOutput(Model, Grid, State)
        type(sifModel_T), intent(inout) :: Model
        type(sifGrid_T), intent(in) :: Grid
        type(sifState_T), intent(in) :: State

        character(len=strLen) :: gStr,tStr

        write (gStr, '(A,I0)') "Step#", Model%SIFIO%nOut

        if (Model%isLoud) then
            if (State%time>1.0e-2) then
                write(tStr,'(f9.3,a)' ) State%time, ' [code]'
            else
                write(tStr,'(es9.2,a)') State%time, ' [code]'
            endif
            write (*, '(a,a,a,a,a)') ANSIGREEN, '<Writing SIF HDF5 DATA @ t = ', trim(tStr), ' >', ANSIRESET
        endif

        call WriteSIF(Model, Grid, State, gStr)

        !Setup for next output
        Model%SIFIO%nOut = Model%SIFIO%nOut + 1
    end subroutine sifOutput

end module sifOut