module shellGridIO

    use kdefs
    use ioH5
    use shellGrid

    implicit none

    integer, parameter, private :: MAXIOVAR = 50
    
    contains

    subroutine writeShellGrid(sg, outH5, gStrO)
        !! Writes ShellGrid data to file outH5
        type(ShellGrid_T), intent(in) :: sg
            !! Shell Grid object to write
        character(len=strLen), intent(in) :: outH5
            !! Output file name
        character(len=strLen), optional, intent(in) :: gStrO
            !! Optional group to write to, default = /ShellGrid

        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars


        call ClearIO(IOVars)

        call AddOutVar(IOVars, "name", sg%name)

        if (present(gStrO)) then
            call WriteVars(IOVars, .false., outH5, gStrO)
        else
            call WriteVars(IOVars, .false., outH5)
        endif


    end subroutine writeShellGrid

end module shellGridIO