module cmetypes
    use gdefs
    use ioclock

    implicit none
    !procedure(cmeSolution_I), pointer :: generateCMESolution => NULL()

    type, abstract :: baseCMEModel_T
        character(len=strLen) :: RunID
        !> Timer for Clocks
        integer :: ts
        !> Whether you can write to console
        logical :: isLoud = .true. 
        !> Whether you can write to debug to console
        logical :: isDebug = .false. 
        !> 
        real(rp) :: Tstart_transient
        !> 
        real(rp) :: time = 0.0
        !> Output info
        type (IOClock_T) :: IO
        procedure(), pointer, nopass :: updateModelTime
    end type

    type, abstract :: baseCMEState_T
      
    end type

    type, abstract :: baseCMESolution_T 
        real(rp), dimension(:, :, :), allocatable :: pres, dens, temp, inside_mask
        real(rp), dimension(:, :, :, :), allocatable :: b, v, j
        !DIR$ attributes align : ALIGN :: pres, dens, temp, inside_mask
        !DIR$ attributes align: ALIGN :: b, v, j
        procedure(), pointer, nopass :: generateSolution
    end type

    abstract interface
        subroutine cmeSolution_I(Solution, Model, State)
            Import baseCMESolution_T, baseCMEModel_T, baseCMEState_T
            class(baseCMESolution_T), intent(inout) :: Solution
            class(baseCMEModel_T), intent(inout) :: Model
            class(baseCMEState_T), intent(inout) :: State
        end subroutine
    end interface

end module