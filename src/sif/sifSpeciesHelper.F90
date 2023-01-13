
module sifSpeciesHelper
    use ioh5
    use xml_input
    use shellgrid

    use sifdefs
    use siftypes

    implicit none 
    contains


    !------
    ! Helper helpers
    !------

    function spcExist(Grid, flav) result(fExist)
        type(SIFGrid_T), intent(in) :: Grid
        integer, intent(in) :: flav
        
        integer :: i
        logical :: fExist

        fExist = .false.

        do i=1,Grid%nSpc
            if (Grid%spc(i)%flav .eq. flav) then
                fExist = .true.
                exit
            endif
        enddo
    end function spcExist


    !------
    ! Do-stuff helpers
    !------

    subroutine populateSpeciesFromConfig(Model, Grid, configfname)
        !! Reads species specification from config file
        !! Populates Grid%spc object
        type(sifModel_T)  , intent(inout) :: Model
        type(sifGrid_T), intent(inout) :: Grid
        character(len=strLen), intent(in) :: configfname

        integer(HID_T) :: h5fId
        integer :: herr, i
        character(len=strLen) :: gStr, tmpStr
        integer(HID_T) :: gId
        logical :: gExist, isEnd
        type(IOVAR_T), dimension(2) :: IOVars ! Just grabbing lami and we're done
        type(IOVAR_T) :: IOV


        call CheckFileOrDie(configfname,"Unable to open file")
        call h5open_f(herr) !Setup Fortran interface
        ! Open file
        call h5fopen_f(trim(configfname), H5F_ACC_RDONLY_F, h5fId, herr)
        
        ! Make sure Species group is there
        gStr = "Species"
        call h5lexists_f(h5fId, gStr, gExist, herr)
        if (.not. gExist) then
            write(*,*) "This config file not structured for SIF, use genSIF.py. Good day."
            stop
        endif
        
        ! If still here, start populating Grid%spc
        allocate(Grid%spc(Grid%nSpc))

        
        associate(spc=>Grid%spc)
        
        ! TODO: Think through edge cases that will cause errors

        isEnd = .false.
        do i=1,Grid%nSpc 
            ! We do an extra iteration to see if there are more species in the config than what SIF expected
            ! If so, we can warn the user
            write(gStr, '(A,I0)') "Species/",i-1  ! Species indexing in config starts at 0
            call h5lexists_f(h5fId,gStr,gExist,herr)
            if (.not. gExist) then
                write(*,*) "ERROR in sifGrids.F90:populateSpeciesFromConfig"
                write(*,'(A,I0,A,I0)') "  Expected ",Grid%nSpc," species but only found ",i-1
                stop
            else
                !write(*,*)"Found spc group ",trim(gStr)
                
                call h5gopen_f(h5fId,trim(gStr),gId,herr)
                ! TODO: Figure out how to get Name string from attrs
                spc(i)%N = readIntHDF(gId, "N")
                spc(i)%flav = readIntHDF(gId, "flav")
                spc(i)%fudge = readRealHDF(gId, "fudge")
                
                !write(*,*)spc(i)%N
                !write(*,*)spc(i)%flav
                !write(*,*)spc(i)%fudge

                ! Now get our alami values
                allocate(spc(i)%alami(spc(i)%N+1))
                
                ! Do this manually, can't use IOVARS/ReadVars cause we already have the group open
                IOV%toRead = .true.
                IOV%idStr = "alami"
                IOV%vType = IONULL
                IOV%scale = 1.0
                call ReadHDFVar(IOV, gId)
                spc(i)%alami = IOV%data
                !write(*,*)spc(i)%alami

                call h5gclose_f(gId,herr)
                
            endif
        enddo

        end associate


        ! Check if there are more species in the file than what we read
        write(gStr, '(A,I0)') "Species/",Grid%nSpc  ! Species indexing in config starts at 0
        call h5lexists_f(h5fId,gStr,gExist,herr)
        if(gExist) then
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
            call assertSpcExist(PSPH) ! Plasmasphere
        endif
        call assertSpcExist(HOTE)  ! Hot electrons
        call assertSpcExist(HOTP)  ! Hot protons

        
        

        call h5fclose_f(h5fId,herr)
        call h5close_f(herr) !Close intereface

        contains

        subroutine assertSpcExist(flav)
            integer, intent(in) :: flav

            ! TODO: Make this say the name of the missing species, e.g. plasmasphere, hot protons/electrons
            if(.not. spcExist(Grid, flav)) then
                write(*,'(A,I0,A,A)')" SIF ERROR: Expected a species with flav ",flav,", but ",trim(configfname)," did not contain one"
                stop
            endif

        end subroutine assertSpcExist

    end subroutine populateSpeciesFromConfig

end module sifSpeciesHelper