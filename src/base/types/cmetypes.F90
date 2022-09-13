module cmetypes
    use gdefs
    use ioclock

    procedure(cmeSolution_I), pointer :: generateCMESolution => NULL()

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
        real(rp) :: time
        !> Output info
        type (IOClock_T) :: IO
    end type

    type, abstract :: baseCMEState_t
      
    end type

    type, abstract :: baseCMESolution_T 
        real(rp), dimension(:, :, :), allocatable :: pres, dens, temp, inside_mask
        real(rp), dimension(:, :, :, :), allocatable :: b, v, j
        !DIR$ attributes align : ALIGN :: pres, dens, temp
        !DIR$ attributes align: ALIGN :: b, v, j
    end type

    abstract interface
        subroutine cmeSolution_I(Solution, Model, State)
            Import baseCMESolution_T
            class(*), intent(inout) :: Solution
            class(*), intent(inout) :: Model
            class(*), intent(inout) :: State
        end subroutine
    end interface

end module