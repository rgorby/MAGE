module glio5
    use glutils
    use ioH5
    use dates
    use files

    implicit none

    integer, parameter, private :: MAXIOVAR = 50
    type(IOVAR_T), dimension(MAXIOVAR), private :: IOVars
    logical, private :: doRoot = .true. !Whether root variables need to be written
    logical, private :: doFat = .false. !Whether to output lots of extra data

    !Necessary for IO routines
    character(len=strLen), public:: GLH5File

    contains

        !> Initialize h5 file and write standalone model grid
        !> 
        subroutine writeH5GridInit(Model, State)
            type(glModel_T), intent(in) :: Model
            type(glState_T),  intent(in) :: State
            logical :: isExist

            character(len=strLen) :: vID
            !Test if root variables (grid/force info is already there)
            vID = "X" !Value to test for
            isExist = ioExist(GLH5File,vID)

            if (isExist) then
                return
            endif

            doRoot = .false.

            !Reset IO chain
            call ClearIO(IOVars)

            call AddOutVar(IOVars,"R", State%r)
            call AddOutVar(IOVars,"THETA", State%thpb)
            call AddOutVar(IOVars,"PHI", State%phpb)

            call AddOutVar(IOVars, "X", State%xyz(:,:,:,XDIR))
            call AddOutVar(IOVars, "Y", State%xyz(:,:,:,YDIR))
            call AddOutVar(IOVars, "Z", State%xyz(:,:,:,ZDIR))
            call WriteVars(IOVars,.true.,GLH5File)
        end subroutine writeH5GridInit

        !> Write solution variables to h5 file
        !>
        subroutine writeSolution(Model, State, Solution, gStr)
            type(glModel_T), intent(in) :: Model
            type(glState_T),  intent(in) :: State
            type(glSolution_T), intent(in) :: Solution
            character(len=*), intent(in) :: gStr
            
            type(glSolution_T) :: SolutionCartesian

            !Check if root variables need to be written on first write
            if (doRoot) then
                call writeH5GridInit(Model,State)
            endif

            !Reset IO chain
            call ClearIO(IOVars)

            SolutionCartesian = solutionSphereToCartesian(State%xyz, Model, State, Solution)
            call AddOutVar(IOVars,"dens",Solution%dens)
            call AddOutVar(IOVars,"pressure",Solution%pres)
            call AddOutVar(IOVars,"temp",Solution%temp)
            call AddOutVar(IOVars,"Br",Solution%b(:,:,:,XDIR))
            call AddOutVar(IOVars,"Bt",Solution%b(:,:,:,YDIR))
            call AddOutVar(IOVars,"Bp",Solution%b(:,:,:,ZDIR))           
            call AddOutVar(IOVars,"Jr",Solution%j(:,:,:,XDIR))
            call AddOutVar(IOVars,"Jt",Solution%j(:,:,:,YDIR))
            call AddOutVar(IOVars,"Jp",Solution%j(:,:,:,ZDIR))            
            call AddOutVar(IOVars,"Vr",Solution%v(:,:,:,XDIR))
            call AddOutVar(IOVars,"Bx",SolutionCartesian%b(:,:,:,XDIR))
            call AddOutVar(IOVars,"By",SolutionCartesian%b(:,:,:,YDIR))
            call AddOutVar(IOVars,"Bz",SolutionCartesian%b(:,:,:,ZDIR))
            call AddOutVar(IOVars,"Jx",SolutionCartesian%j(:,:,:,XDIR))
            call AddOutVar(IOVars,"Jy",SolutionCartesian%j(:,:,:,YDIR))
            call AddOutVar(IOVars,"Jz",SolutionCartesian%j(:,:,:,ZDIR))            
            call AddOutVar(IOVars,"Vx",SolutionCartesian%v(:,:,:,XDIR))
            call AddOutVar(IOVars,"Vy",SolutionCartesian%v(:,:,:,XDIR))
            call AddOutVar(IOVars,"Vz",SolutionCartesian%v(:,:,:,XDIR))
            call AddOutVar(IOVars,"insideMask",Solution%inside_mask)
            call AddOutVar(IOVars,"time",Model%time)
            call AddOutVar(IOVars,"phiss",Model%phiss)
            call AddOutVar(IOVars,"timestep",Model%ts)
            !------------------
            !Finalize
            call WriteVars(IOVars,.true.,GLH5File,trim(gStr))

        end subroutine writeSolution

        !>
        !>
        subroutine setH5File(filename)
            character(len=strLen), intent(in) :: filename
            GLH5File = filename
        end subroutine setH5File

end module glio5