module voltappHelper

    use kdefs
    use voltTypes
    use imagtubes

    implicit none

    contains

    subroutine initVoltState(vApp)
        type(voltApp_T), intent(inout) :: vApp

        associate(sh=>vApp%shGrid, State=>vApp%State)
            allocate(State%ijTubes(sh%Nt+1, sh%Np+1))
            !write(*,*) "tubeShell"
            !call init_IMAGTubeShell(sh, State%tubeShell)
            call initShellVar(sh, SHGR_CORNER, State%potential)

        end associate 
    end subroutine initVoltState

end module voltappHelper