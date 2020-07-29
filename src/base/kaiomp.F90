!Routines to handle OMP logic, default settings, XML reading
module kaiomp

    use kdefs
    use xml_input
    use strings
    
#ifdef _OPENMP
    use omp_lib
#endif    
    implicit none

#ifdef _OPENMP
    logical, parameter :: isOMP = .true.
#else
    logical, parameter :: isOMP = .false.
#endif

    contains

    !Read parameters from xxx/omp block
    subroutine SetOMP(xmlInp,doLoudO)
        type(XML_Input_T), intent(inout) :: xmlInp
        logical, intent(in), optional :: doLoudO

        logical :: doLoud
        integer :: NumTh,MaxTh,NumP
        character(len=strLen) :: rStr

        if (present(doLoudO)) then
            doLoud = doLoudO
        else
            doLoud = .false.
        endif

        call xmlInp%GetRootStr(rStr)

        if (.not. isOMP) then
            !OMP not enabled
            if (doLoud) then
                write(*,*) trim(toUpper(rStr)) // " running without threading ..."
            endif
            return
        endif
        
        !If you're still here, let's get our omp on
        !Find max possible threads
#ifdef _OPENMP
        MaxTh = omp_get_max_threads()
        NumP  = omp_get_num_procs() 
#else
        MaxTh = 0 !This shouldn't happen
        NumP  = 0
#endif
        call xmlInp%Set_Val(NumTh,"threading/NumTh",MaxTh)
        !Check for secret values (-x => x/core)
        if (NumTh<0) then
            NumTh = (-NumTh)*NumP
        endif

        !Now set threads
#ifdef _OPENMP
        call omp_set_num_threads(NumTh)
#endif
        !Talk about it
        if (doLoud) then
            write(*,*) trim(toUpper(rStr)) // " running threaded ..."
            write(*,*) '   # Threads = ', NumTh
            write(*,*) '   # Cores   = ', NumP

        endif

    end subroutine SetOMP

    !Returns current number of threads
    function NumOMP()
        integer :: NumOMP
#ifdef _OPENMP
        NumOMP = omp_get_num_threads()
#else
        NumOMP = 0
#endif
    end function NumOMP

    !Returns thread ID within group threads
    !Returns 0 in non-threaded code
    function ThreadID()
        integer :: ThreadID

#ifdef _OPENMP
        ThreadID = omp_get_thread_num()+1 !Start w/ 1
#else
        ThreadID = 0
#endif
    end function ThreadID


end module kaiomp
