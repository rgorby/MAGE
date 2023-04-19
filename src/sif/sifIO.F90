module sifIO
    use ioh5
    use planethelper
    
    use siftypes
    use sifetautils

    implicit none

    integer, parameter, private :: MAXIOVAR = 50
    type(IOVAR_T), dimension(MAXIOVAR), private :: IOVars
    logical, private :: doRoot = .true. !Whether root variables need to be written
    logical, private :: doFat = .false. !Whether to output lots of extra datalogical, private :: doRoot = .true. !Whether root variables need to be written

    contains

    subroutine sifInitIO(Model, Grid)
        type(sifModel_T), intent(inout) :: Model
        type(sifGrid_T), intent(in) :: Grid

        integer :: i
        logical :: fExist
        real(rp), dimension(:,:), allocatable :: lat2D, lon2D
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        character(len=strLen) :: gStr

        doRoot = .false. ! Don't call again

        associate(sh => Grid%shGrid, spc=>Grid%spc)
        Model%SIFH5 = trim(Model%RunID) // ".sif.h5"

        fExist = CheckFile(Model%SIFH5)
        write(*,*) "SIF outputting to ",trim(Model%SIFH5)

        if(.not. Model%isRestart) then
            ! Remove all old files, start fresh
            call CheckAndKill(Model%SIFH5)
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

        ! Grid data
        call AddOutVar(IOVars,"X",lat2D,uStr="radians")
        call AddOutVar(IOVars,"Y",lon2D,uStr="radians")
        call AddOutVar(IOVars,"alamc",Grid%alamc,uStr="eV * (Rx/nT)^(2/3)")
        call WriteVars(IOVars,.true.,Model%SIFH5)

        ! Output detailed lambda grid info
        do i=1,Grid%nSpc
            call ClearIO(IOVars)
            ! Attrs
            call AddOutVar(IOVars,"flav"     ,spc(i)%flav    )
            call AddOutVar(IOVars,"N"        ,spc(i)%N       )
            call AddOutVar(IOVars,"fudge"    ,spc(i)%fudge   )
            call AddOutVar(IOVars,"numNuc_p" ,spc(i)%numNuc_p)
            call AddOutVar(IOVars,"numNuc_n" ,spc(i)%numNuc_n)
            call AddOutVar(IOVars,"q       " ,spc(i)%q       )
            call AddOutVar(IOVars,"amu"      ,spc(i)%amu     )
            call AddOutVar(IOVars,"kStart"   ,spc(i)%kStart-1)  ! Change to assume zero-based
            call AddOutVar(IOVars,"kEnd"     ,spc(i)%kEnd  -1)  ! Change to assume zero-based
            ! Datasets
            call AddOutVar(IOVars,"alami",spc(i)%alami,uStr="eV * (Rx/nT)^(2/3)")
            write(gStr,'(I0)') spc(i)%flav  ! Idk if this is the easiest way to format ints as strings
            call WriteVars(IOVars,.true.,Model%SIFH5,"Species",gStr)
        enddo

        ! Output planet info
        call writePlanetParams(Model%planet, .true., Model%SIFH5)

        end associate

    end subroutine sifInitIO


    subroutine WriteSIF(Model, Grid, State, gStr)
        type(sifModel_T), intent(inout) :: Model
        type(sifGrid_T ), intent(in) :: Grid
        type(sifState_T), intent(in) :: State
        character(len=strLen), intent(in) :: gStr

        integer :: i,j,s
        real(rp), dimension(:,:,:), allocatable :: outDen, outIntensity


        ! First, make sure root variables are there
        if (doRoot) then
            call sifInitIO(Model, Grid)
            doRoot = .false.
        endif
        !Reset IO chain
        call ClearIO(IOVars)

        ! Add attributes
        call AddOutVar(IOVars,"time",State%t)

        ! Add State variables
        call AddOutVar(IOVars,"bmin",State%Bmin,uStr="nT")
        call AddOutVar(IOVars,"xmin",State%xyzmin(:,:,XDIR),uStr="Rx")
        call AddOutVar(IOVars,"ymin",State%xyzmin(:,:,YDIR),uStr="Rx")
        call AddOutVar(IOVars,"zmin",State%xyzmin(:,:,ZDIR),uStr="Rx")

        call AddOutVar(IOVars,"eta",State%eta,uStr="#/cm^3 * Rx/T") !! TODO: Maybe swap with just intensity instead

        ! Calc intensity
        allocate(outIntensity(Grid%shGrid%Nt,Grid%shGrid%Np,Grid%Nk))
        outIntensity = 0.0
        do s=1,Grid%nSpc
            if (Grid%spc(s)%flav==F_PSPH) then
                cycle  ! Skip plasmasphere since it has zero energy
            endif
            do i=1,Grid%shGrid%Nt
                do j=1,Grid%shGrid%Np
                    outIntensity(i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd) = &
                        eta2intensity(Grid%spc(s)%amu,   &
                                      State%bVol(i,j),   &
                                      Grid%spc(s)%alami, &
                                      State%eta (i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd))
                enddo
            enddo
        enddo
        call AddOutVar(IOVars,"intensity",outIntensity,uStr="1/(s*sr*keV*cm^2)")
        deallocate(outIntensity)
        

        call AddOutVar(IOVars,"topo",State%topo*1.0_rp,uStr="0=Open, 1=Closed")
        call AddOutVar(IOVars,"active",State%active*1.0_rp,uStr="-1=Inactive, 0=Buffer, 1=Active")
        call AddOutVar(IOVars,"espot",State%espot,uStr="kV")
        call AddOutVar(IOVars,"colatc",State%thc,uStr="radians")
        call AddOutVar(IOVars,"lonc"  ,State%phc,uStr="radians")
        call AddOutVar(IOVars,"bVol",State%bvol,uStr="Rx/nT")

    ! Moments
        call AddOutVar(IOVars,"Pressure",State%Press,uStr="nPa")
        ! Add density moment as #/cc instead of amu/cc
        allocate(outDen(Grid%shGrid%Nt,Grid%shGrid%Np,Grid%nSpc+1))
        outDen = 0.0
        do s=1, Grid%nSpc
            outDen(:,:,s+1) = State%Den(:,:,s+1)/Grid%spc(s)%amu
            ! Don't include electrons to total number density
            if(Grid%spc(s)%isElectron .eq. .false.) then
                outDen(:,:,1) = outDen(:,:,1) + outDen(:,:,s+1)
            endif
        enddo
        call AddOutVar(IOVars,"Density",outDen,uStr="#/cc")
        deallocate(outDen)


        call WriteVars(IOVars,.true.,Model%SIFH5, gStr)

    end subroutine WriteSIF

end module sifIO