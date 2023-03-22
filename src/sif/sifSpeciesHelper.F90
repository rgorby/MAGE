
module sifSpeciesHelper
    use ioh5
    use xml_input
    use shellgrid

    use sifdefs
    use siftypes

    implicit none 

    integer, parameter, private :: MAXIOVAR = 50

    contains
    
    !------
    ! Helper helpers
    !------

    function spcIdx(Grid, flav) result(idx)
        !! Get the index corresponding to a certain species
        type(SIFGrid_T), intent(in) :: Grid
        integer, intent(in) :: flav
        
        integer :: i, idx

        do i=1,Grid%nSpc
            if (Grid%spc(i)%flav .eq. flav) then
                idx = i
                return
            endif
        enddo

        ! If we exited the loop that means flav not in species list
        idx = -1
        
    end function

    function spcExists(Grid, flav) result(fExist)
        !! Determine if a certain species is in the species list, using flavor value
        type(SIFGrid_T), intent(in) :: Grid
        integer, intent(in) :: flav
        
        integer :: i
        logical :: fExist

        fExist = .false.
        if(spcIdx(Grid, flav) .ne. -1) then
            fExist = .true.
        endif
    end function spcExists


    !------
    ! Do-stuff helpers
    !------

    subroutine populateSpeciesFromConfig(Model, Grid, configfname)
        !! Reads species specification from config file
        !! Populates Grid%spc object
        type(sifModel_T)  , intent(inout) :: Model
        type(sifGrid_T), intent(inout) :: Grid
        character(len=strLen), intent(in) :: configfname

        integer :: i
        character(len=strLen) :: gStr
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars ! Just grabbing lami and we're done


        ! Make sure Species group exists with at least 1 entry
        if(.not. ioExist(configfname, "0", "Species")) then
            write(*,*) "This config file not structured for SIF, use genSIF.py. Good day."
            stop
        endif
        
        ! If still here, start populating Grid%spc
        allocate(Grid%spc(Grid%nSpc))

        
        ! TODO: Think through edge cases that will cause errors
        
        do i=1,Grid%nSpc 
            associate(spc=>Grid%spc(i))

            write(gStr, '(A,I0)') "Species/",i-1  ! Species indexing in config starts at 0
            if(.not. ioExist(configfname, gStr)) then
                write(*,*) "ERROR in sifGrids.F90:populateSpeciesFromConfig"
                write(*,'(A,I0,A,I0)') "  Expected ",Grid%nSpc," species but only found ",i-1
                stop
            endif
            
            ! If still here, continue with reading and assigning
            
            write(*,*)"Found spc group ",trim(gStr)
                
            ! Read
            call ClearIO(IOVars)
            call AddInVar(IOVars, "flav" )  ! Attr
            call AddInVar(IOVars, "N"    )  ! Attr
            call AddInVar(IOVars, "amu")  ! Attr
            call AddInVar(IOVars, "fudge")  ! Attr
            call AddInVar(IOVars, "alami")  ! Dataset
            call ReadVars(IOVars, .false., configfname, gStr)

            ! Assign
            spc%flav  = IOVars(FindIO(IOVars, "flav" ))%data(1)
            spc%N     = IOVars(FindIO(IOVars, "N"    ))%data(1)
            spc%amu   = IOVars(FindIO(IOVars, "amu"  ))%data(1)
            spc%fudge = IOVars(FindIO(IOVars, "fudge"))%data(1)
            spc%isElectron = (spc%amu*100.0 < 1) ! Gonna assume this will always work, but is dependent on config having the right amu value

            allocate(spc%alami(spc%N+1))
            call IOArray1DFill(IOVars,"alami",spc%alami)

            write(*,*)" Flav: ", spc%flav
            write(*,*)" N:    ", spc%N
            write(*,*)" amu:  ", spc%amu
            write(*,*)" Fudge:", spc%fudge
            !write(*,*) spc%alami

                

            end associate
        enddo
        

        ! Check if there are more species in the file than what we read
        write(gStr, '(A,I0)') "Species/",Grid%nSpc  ! Species indexing in config starts at 0
        if(ioExist(configfname, gStr)) then
            if (.not. Grid%ignoreConfigMismatch) then
                write(*,*) "SIF ERROR: More species defined in ",trim(configfname)," than SIF expected."
                write(*,*) " If you want to ignore this, set config/ignoreMismatch=T. But for now, we die."
                stop
            else
                write(*,*)"SIF WARNING: More species in config than we expect, but ignoring because config/ignoreMismatch=T"
            endif
        endif

        ! Ensure that we have the core species we expect to have
        if (Model%doPlasmasphere) then
            call assertSpcExists(F_PSPH) ! Plasmasphere
        endif
        call assertSpcExists(F_HOTE)  ! Hot electrons
        call assertSpcExists(F_HOTP)  ! Hot protons


        contains

        subroutine assertSpcExists(flav)
            integer, intent(in) :: flav

            ! TODO: Make this say the name of the missing species, e.g. plasmasphere, hot protons/electrons
            if(.not. spcExists(Grid, flav)) then
                write(*,'(A,I0,A,A)')" SIF ERROR: Expected a species with flav ",flav,", but ",trim(configfname)," did not contain one"
                stop
            endif

        end subroutine assertSpcExists

    end subroutine populateSpeciesFromConfig

    subroutine initAlamc(Grid)
        !! Allocates and Populates Grid%alamc using species
        !!  as well as tells each species what its range is inside alamc
        type(sifGrid_T), intent(inout) :: Grid

        integer :: i, kPos

        ! Calc Nk
        Grid%Nk = 0
        do i=1,Grid%nSpc
            Grid%Nk = Grid%Nk + Grid%spc(i)%N
        enddo

        ! Allocate
        allocate(Grid%alamc(Grid%Nk))
        Grid%alamc = 0.0

        ! Populate, tell each species its bounds
        kPos = 1
        do i=1,Grid%nSpc
            associate(spc=>Grid%spc(i))
            spc%kStart = kPos
            kPos = kPos + spc%N
            spc%kEnd = kPos-1

            ! N = # channels, and size(alami)=N+1
            Grid%alamc(spc%kStart:spc%kEnd) = 0.5*(spc%alami(2:spc%N+1) + spc%alami(1:spc%N))
            end associate
        enddo

        ! TODO: Write unit test to ensure this is working as expected

    end subroutine initAlamc

end module sifSpeciesHelper