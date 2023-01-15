module sifIO
    use ioh5
    
    use siftypes

    implicit none

    contains

    subroutine sifInitIO(Model, Grid)
        type(sifModel_T), intent(inout) :: Model
        type(sifGrid_T), intent(in) :: Grid

        integer :: i
        logical :: fExist
        real(rp), dimension(:,:), allocatable :: lat2D, lon2D
        type(IOVAR_T), dimension(5) :: IOVars
        character(len=strLen) :: gStr

        associate(SIO => Model%SIFIO, sh => Grid%shGrid, spc=>Grid%spc)
        SIO%SIFH5 = trim(Model%RunID) // ".sif.h5"

        fExist = CheckFile(SIO%SIFH5)
        write(*,*) "SIF outputting to ",trim(SIO%SIFH5)

        if(.not. Model%isRestart) then
            ! Remove all old files, start fresh
            call CheckAndKill(SIO%SIFH5)
        endif

        if (Model%isRestart .and. fExist) then
            ! No need to do anything because file already exist
            !!TODO: When more output files are active, maybe put in a check for each of them
            return
        endif

        ! If still here, proceed to init

        ! Add spatial grid as 2D array
        allocate(lat2D(sh%Nt+1,sh%Np+1))  ! +1 because we're doing corners
        allocate(lon2D(sh%Nt+1,sh%Np+1))

        do i=1,sh%Np+1
            lat2D(:,i) = sh%th(sh%is:sh%ie)
        enddo

        do i=1,sh%Nt+1
            lon2D(i,:) = sh%ph(sh%js:sh%je)
        enddo

        ! Ready for output
        !! TODO: Stamp with git hash and branch
        call ClearIO(IOVars)

        call AddOutVar(IOVars,"X",lat2D)
        call AddOutVar(IOVars,"Y",lon2D)
        call AddOutVar(IOVars,"alamc",Grid%alamc)
        call AddOutVar(IOVars,"UnitsID","Rad")  ! Attribute
        call WriteVars(IOVars,.true.,SIO%SIFH5)

        ! Output detailed lambda grid info
        do i=1,Grid%nSpc
            call ClearIO(IOVars)
            call AddOutVar(IOVars,"flav" ,spc(i)%flav )
            call AddOutVar(IOVars,"N"    ,spc(i)%N    )
            call AddOutVar(IOVars,"fudge",spc(i)%fudge)
            call AddOutVar(IOVars,"alami",spc(i)%alami)
            write(gStr,'(I0)') spc(i)%flav  ! Idk if this is the easiest way to format ints as strings
            call WriteVars(IOVars,.true.,SIO%SIFH5,"Species",gStr)

        enddo

        end associate

    end subroutine sifInitIO

end module sifIO