
module raijuSpeciesHelper
    use ioh5
    use xml_input
    use shellgrid

    use raijudefs
    use raijutypes

    implicit none 

    integer, parameter, private :: MAXIOVAR = 50

    contains
    
!------
! Helper helpers
!------

    function spcIdx(Grid, flav) result(idx)
        !! Get the index corresponding to a certain species
        type(raijuGrid_T), intent(in) :: Grid
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
        type(raijuGrid_T), intent(in) :: Grid
        integer, intent(in) :: flav
        
        integer :: i
        logical :: fExist

        fExist = .false.
        if(spcIdx(Grid, flav) .ne. -1) then
            fExist = .true.
        endif
    end function spcExists


    function SpcAmu(spc) result (amu)
        !! Calculates a species's mass in amu/daltons
        type(raijuSpecies_T), intent(in) :: spc

        integer :: num_e
        real(rp) :: amu

        !! This is technically wrong. 
        !  1 mass_proton = 1.00727647 amu, but He = 4.0026 amu = 3.9736 mass_proton
        !  Because binding energy is a thing
        !  Also mass_proton != mass_neutron, but that's an even smaller difference
        !  If this becomes our biggest error we're doing pretty well
        
        num_e = spc%numNuc_p - spc%q
        amu = ((spc%numNuc_p + spc%numNuc_n)*Mp_cgs + num_e*Me_cgs)*1.e-3/dalton

    end function SpcAmu


    function SpcType(spc) result(sType)
        !! Determine the species type (e-, H+, O+, etc.)
        type(raijuSpecies_T), intent(in) :: spc

        integer :: sType

        select case(spc%numNuc_p)
            case(0)
                sType = RAIJUELE
            case(1)
                sType = RAIJUHPLUS
            case(8)
                sType = RAIJUOPLUS
            case DEFAULT
                write(*,*)"WARNING: RAIJU can't determine type for species with nProton=",spc%numNuc_p
                sType = RAIJUNSPC
        end select


    end function SpcType

    !------
    ! Do-stuff helpers
    !------

    subroutine populateSpeciesFromConfig(Model, Grid, configfname)
        !! TODO: Rewrite species determination. xml config sshould determine which species we want, and we get the details from raijuconfig.h5
        !! Reads species specification from config file
        !! Populates Grid%spc object
        type(raijuModel_T)  , intent(inout) :: Model
        type(raijuGrid_T), intent(inout) :: Grid
        character(len=strLen), intent(in) :: configfname

        integer :: i, kPos, f, ferr, psphIdx
        character(len=strLen) :: gStr
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars ! Just grabbing lami and we're done


        ! Make sure Species group exists with at least 1 entry
        if(.not. ioExist(configfname, "0", "Species")) then
            write(*,*) "This config file not structured for RAIJU, use genRAIJU.py. Good day."
            stop
        endif
        
        ! If still here, start populating Grid%spc
        allocate(Grid%spc(Grid%nSpc))
        
        ! TODO: Think through edge cases that will cause errors
        kPos = 1
        do i=1,Grid%nSpc 
            associate(spc=>Grid%spc(i))

            write(gStr, '(A,I0)') "Species/",i-1  ! Species indexing in config starts at 0
            if(.not. ioExist(configfname, gStr)) then
                write(*,*) "ERROR in raijuSpeciesHelper.F90:populateSpeciesFromConfig"
                write(*,'(A,I0,A,I0)') "  Expected ",Grid%nSpc," species but only found ",i-1
                stop
            endif
            
            ! If still here, continue with reading and assigning
            write(*,*)"Found spc group ",trim(gStr)
                
            ! Read
            call ClearIO(IOVars)
            call AddInVar(IOVars, "flav"    )  ! Attr
            call AddInVar(IOVars, "N"       )  ! Attr
            call AddInVar(IOVars, "numNuc_p")  ! Attr
            call AddInVar(IOVars, "numNuc_n")  ! Attr
            call AddInVar(IOVars, "q"       )  ! Attr
            call AddInVar(IOVars, "fudge"   )  ! Attr
            call AddInVar(IOVars, "alami"   )  ! Dataset
            call ReadVars(IOVars, .false., configfname, gStr)

            ! Assign
            spc%flav     = IOVars(FindIO(IOVars, "flav"    ))%data(1)
            spc%N        = IOVars(FindIO(IOVars, "N"       ))%data(1)
            spc%numNuc_p = IOVars(FindIO(IOVars, "numNuc_p"))%data(1)
            spc%numNuc_n = IOVars(FindIO(IOVars, "numNuc_n"))%data(1)
            spc%q        = IOVars(FindIO(IOVars, "q"       ))%data(1)
            spc%fudge    = IOVars(FindIO(IOVars, "fudge"   ))%data(1)
            
            spc%amu = SpcAmu(spc)
            spc%spcType = SpcType(spc)

            !! If species is H+ but not plasmasphere, we can map etas below its lambda bounds to plasmasphere
            !if (spc%spcType .eq. RAIJUHPLUS .and. spc%flav .ne. F_PSPH) then
            !    spc%mapExtraToPsph = .true.
            !else
            !    spc%mapExtraToPsph = .false.
            !endif

            ! Calc start and end bounds, use it to set alami index range
            spc%kStart = kPos
            kPos = kPos + spc%N
            spc%kEnd = kPos-1
        
            allocate(spc%alami(spc%kStart:spc%kEnd+1))
            call IOArray1DFill(IOVars,"alami",spc%alami)

            write(*,*)" Flav:         ", spc%flav
            write(*,*)" N:            ", spc%N
            write(*,*)" # nuc_p:      ", spc%numNuc_p
            write(*,*)" # nuc_n:      ", spc%numNuc_n
            write(*,*)" q:            ", spc%q
            write(*,*)" spcType:      ", spc%spcType
            write(*,*)" Fudge:        ", spc%fudge
            write(*,*)" kStart:       ", spc%kStart
            write(*,*)" kEnd:         ", spc%kEnd
            !write(*,*)" ExcesstoPsph: ", spc%mapExtraToPsph
            !write(*,*)" alami:   ", spc%alami
            

            end associate
        enddo
        

        ! Check if there are more species in the file than what we read
        ! TODO: Maybe instead:
        !  Read desired flavors from xml file (default to psph, hote, hotp)
        !  Read config file and load in flavors that match
        !  If desired flavors are missing from config, complain and die
        write(gStr, '(A,I0)') "Species/",Grid%nSpc  ! Species indexing in config starts at 0
        if(ioExist(configfname, gStr)) then
            if (.not. Grid%ignoreConfigMismatch) then
                write(*,*) "RAIJU ERROR: More species defined in ",trim(configfname)," than RAIJU expected."
                write(*,*) " If you want to ignore this, set config/ignoreMismatch=T. But for now, we die."
                stop
            else
                write(*,*)"RAIJU WARNING: More species in config than we expect, but ignoring because config/ignoreMismatch=T"
            endif
        endif

        ! Ensure that we have the core species we expect to have
        if (Model%doPlasmasphere) then
            call assertSpcExists(F_PSPH) ! Plasmasphere
        endif
        call assertSpcExists(F_HOTE)  ! Hot electrons
        call assertSpcExists(F_HOTP)  ! Hot protons


        ! Set coupling info
        if (Model%nFluidIn > 0) then
            do f=1,Model%nFluidIn
                ! First, make sure mapping info is good on our side
                if (.not. spcExists(Grid, Model%fluidInMaps(f)%flav)) then
                    write(*,*)"ERROR initLambdaGrid: one of the fluidInMaps has a flavor that doesn't exist!"
                    write(*,*)"Existing species (see raijudefs for enum names):"
                    do ferr=1,Grid%nSpc
                        write(*,*)"flavor=",Grid%spc(ferr)%flav,", spc =",Grid%spc(ferr)%spcType
                    enddo
                    write(*,*)"FluidIn's flav = ",Model%fluidInMaps(f)%flav
                    stop
                endif

                ! If good, determine which species expect to get written to
                i = spcIdx(Grid, Model%fluidInMaps(f)%flav)
                Grid%spc(i)%isMappedTo = .true.
                if (Model%doPlasmasphere .and. Model%fluidInMaps(f)%doExcessToPsph) then
                    psphIdx = spcIdx(Grid, F_PSPH)
                    Grid%spc(psphIdx)%isMappedTo = .true.
                endif
            enddo
        endif

        do i=1,Grid%nSpc
            write(*,*)Grid%spc(i)%flav,Grid%spc(i)%isMappedTo
        enddo

        contains

        subroutine assertSpcExists(flav)
            integer, intent(in) :: flav

            ! TODO: Make this say the name of the missing species, e.g. plasmasphere, hot protons/electrons
            if(.not. spcExists(Grid, flav)) then
                write(*,'(A,I0,A,A)')" RAIJU ERROR: Expected a species with flav ",flav,", but ",trim(configfname)," did not contain one"
                stop
            endif

        end subroutine assertSpcExists

    end subroutine populateSpeciesFromConfig

    subroutine initAlamc(Grid)
        !! Allocates and Populates Grid%alamc using species
        !!  as well as tells each species what its range is inside alamc
        type(raijuGrid_T), intent(inout) :: Grid

        integer :: i

        ! Calc Nk
        Grid%Nk = 0
        do i=1,Grid%nSpc
            Grid%Nk = Grid%Nk + Grid%spc(i)%N
        enddo

        ! Allocate
        allocate(Grid%alamc(Grid%Nk))
        allocate(Grid%k2spc(Grid%Nk))
        Grid%alamc = 0.0
        Grid%k2spc = 0.0

        ! Populate alamc, fill in k2spc map
        do i=1,Grid%nSpc
            associate(spc=>Grid%spc(i))
            ! N = # channels, and size(alami)=N+1
            Grid%alamc(spc%kStart:spc%kEnd) = 0.5*(spc%alami(spc%kStart+1:spc%kEnd+1) + spc%alami(spc%kStart:spc%kEnd))
            Grid%k2spc(spc%kStart:spc%kEnd) = i
            end associate
        enddo

        ! TODO: Write unit test to ensure this is working as expected

    end subroutine initAlamc

end module raijuSpeciesHelper