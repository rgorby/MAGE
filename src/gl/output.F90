module gloutput
    use gldefs
    use gltypes
    use glio5
    implicit none

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

        subroutine fOutput(Model, State, Solution)
            type(glModel_T), intent(inout) :: Model
            type(glSolution_T), intent(in) :: Solution
            type(glState_T), intent(in) :: State

            character(len=strLen) :: gStr,tStr

            write (gStr, '(A,I0)') "Step#", Model%IO%nOut

            if (Model%isLoud) then
                call timeString(Model%time,tStr)
                write (*, '(a,a,a,a,a)') ANSIGREEN, '<Writing HDF5 DATA @ t = ', trim(tStr), ' >', ANSIRESET
            endif

            call writeSolution(Model, State, Solution, gStr)

            !Setup for next output
            Model%IO%tOut = Model%IO%tOut + Model%IO%dtOut
            Model%IO%nOut = Model%IO%nOut + 1
        end subroutine fOutput
        
        subroutine tStr_STD(T,tStr)
            real(rp), intent(in) :: T
            character(len=strLen), intent(out) :: tStr
    
            if (T>1.0e-2) then
                write(tStr,'(f9.3,a)' ) T, ' [code]'
            else
                write(tStr,'(es9.2,a)') T, ' [code]'
            endif
        end subroutine tStr_STD

end module gloutput