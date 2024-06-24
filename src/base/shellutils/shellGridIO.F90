module shellGridIO

    use kdefs
    use ioH5
    use shellGrid

    implicit none

    integer, parameter, private :: MAXIOVAR = 50
    
    contains

    subroutine writeShellGrid(sg, outH5, gStrO)
        !! Writes ShellGrid data to file outH5
        !! Note: This does not save every variable,
        !!  just ones needed for restarts and generally useful output
        type(ShellGrid_T), intent(in) :: sg
            !! Shell Grid object to write
        character(len=strLen), intent(in) :: outH5
            !! Output file name
        character(len=strLen), optional, intent(in) :: gStrO
            !! Optional group to write to, default = /ShellGrid

        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars


        call ClearIO(IOVars)

        ! Attrs
        call AddOutVar(IOVars, "name", sg%name)
        call AddOutVar(IOVars, "radius", sg%radius,uStr="Rp")
        call AddOutVar(IOVars, "Nt", sg%Nt)
        call AddOutVar(IOVars, "Np", sg%Np)
        call AddOutVar(IOVars, "nGhosts_n", sg%Ngn)
        call AddOutVar(IOVars, "nGhosts_s", sg%Ngs)
        call AddOutVar(IOVars, "nGhosts_e", sg%Nge)
        call AddOutVar(IOVars, "nGhosts_w", sg%Ngw)
        call AddOutVar(IOVars, "isChild", sg%isChild*1.0_rp)
        call AddOutVar(IOVars, "parentName", trim(sg%parentName))
        call AddOutVar(IOVars, "bndis", sg%bndis)
        call AddOutVar(IOVars, "bndie", sg%bndie)
        call AddOutVar(IOVars, "bndjs", sg%bndjs)
        call AddOutVar(IOVars, "bndje", sg%bndje)

        ! Arrays
        call AddOutVar(IOVars, "theta", sg%th, uStr="radians")
        call AddOutVar(IOVars, "phi"  , sg%ph, uStr="radians")


        if (present(gStrO)) then
            call WriteVars(IOVars, .false., outH5, gStrO)
        else
            call WriteVars(IOVars, .false., outH5)
        endif


    end subroutine writeShellGrid

end module shellGridIO