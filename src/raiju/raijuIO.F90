module raijuIO
    use ioh5
    use files
    use planethelper
    use kai2geo
    
    use raijutypes
    use raijuetautils
    use raijuPreAdvancer, only : calcEffectivePotential
    use raijuELossWM, only : eWMOutput

    implicit none

    integer, parameter, private :: MAXIOVAR = 50
    type(IOVAR_T), dimension(MAXIOVAR), private :: IOVars
    logical, private :: doRoot = .true. !Whether root variables need to be written
    logical, private :: doFat = .false. !Whether to output lots of extra datalogical, private :: doRoot = .true. !Whether root variables need to be written

    contains

    subroutine raijuInitIO(Model, Grid, doGhostsO)
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        logical, optional, intent(in) :: doGhostsO

        integer :: is, ie, js, je
        logical :: doGhosts
        integer :: i
        logical :: fExist
        real(rp), dimension(:,:), allocatable :: lat2D, lon2D
        !type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        character(len=strLen) :: gStr

        doRoot = .false. ! Don't call again

        if (present(doGhostsO)) then
            doGhosts = doGhostsO
        else
            doGhosts = .false.
        endif

        if (doGhosts) then
            is = Grid%shGrid%isg
            ie = Grid%shGrid%ieg
            js = Grid%shGrid%jsg
            je = Grid%shGrid%jeg
        else
            is = Grid%shGrid%is
            ie = Grid%shGrid%ie
            js = Grid%shGrid%js
            je = Grid%shGrid%je
        endif

        associate(sh => Grid%shGrid, spc=>Grid%spc)
        Model%raijuH5 = trim(Model%RunID) // ".raiju.h5"

        fExist = CheckFile(Model%raijuH5)
        write(*,*) "RAIJU outputting to ",trim(Model%raijuH5)

        if(.not. Model%isRestart) then
            ! Remove all old files, start fresh
            call CheckAndKill(Model%raijuH5)
        endif

        if (Model%isRestart .and. fExist) then
            ! No need to do anything because file already exist
            !!TODO: When more output files are active, maybe put in a check for each of them
            return
        endif

        ! If still here, proceed to init

        ! Add spatial grid as 2D array
        allocate(lat2D(is:ie+1,js:je+1))  ! +1 because we're doing corners
        allocate(lon2D(is:ie+1,js:je+1))

        do i=js,je+1
            lat2D(:,i) = sh%th(is:ie)
        enddo

        do i=is,ie+1
            lon2D(i,:) = sh%ph(js:je)
        enddo

        ! Ready for output
        !! TODO: Stamp with git hash and branch
        call ClearIO(IOVars)

        ! Grid data
        call AddOutVar(IOVars,"X",lat2D,uStr="radians")
        call AddOutVar(IOVars,"Y",lon2D,uStr="radians")
        call AddOutVar(IOVars,"Bmag",Grid%Bmag(is:ie,js:je),uStr="nT")
        call AddOutVar(IOVars,"areaCC",Grid%areaCC  (is:ie,js:je),uStr="Rp^2")
        call AddOutVar(IOVars,"areaFace_th",Grid%areaFace(is:ie,js:je,1),uStr="Rp^2")
        call AddOutVar(IOVars,"areaFace_ph",Grid%areaFace(is:ie,js:je,2),uStr="Rp^2")
        call AddOutVar(IOVars,"alamc",Grid%alamc,uStr="eV * (Rx/nT)^(2/3)")
        call WriteVars(IOVars,.true.,Model%raijuH5)

        ! Output detailed lambda grid info
        do i=1,Grid%nSpc
            call ClearIO(IOVars)
            ! Attrs
            call AddOutVar(IOVars,"flav"     ,spc(i)%flav    )
            call AddOutVar(IOVars,"N"        ,spc(i)%N       )
            call AddOutVar(IOVars,"fudge"    ,spc(i)%fudge   )
            call AddOutVar(IOVars,"numNuc_p" ,spc(i)%numNuc_p)
            call AddOutVar(IOVars,"numNuc_n" ,spc(i)%numNuc_n)
            call AddOutVar(IOVars,"spcType"  ,spc(i)%spcType )
            call AddOutVar(IOVars,"q"        ,spc(i)%q       )
            call AddOutVar(IOVars,"amu"      ,spc(i)%amu     )
            call AddOutVar(IOVars,"kStart"   ,spc(i)%kStart-1)  ! Change to assume zero-based
            call AddOutVar(IOVars,"kEnd"     ,spc(i)%kEnd  -1)  ! Change to assume zero-based
            ! Datasets
            call AddOutVar(IOVars,"alami",spc(i)%alami,uStr="eV * (Rx/nT)^(2/3)")
            write(gStr,'(I0)') spc(i)%flav  ! Idk if this is the easiest way to format ints as strings
            call WriteVars(IOVars,.true.,Model%raijuH5,"Species",gStr)
        enddo

        ! Output planet info
        call writePlanetParams(Model%planet, .true., Model%raijuH5)

        end associate

    end subroutine raijuInitIO


    subroutine WriteRaiju(Model, Grid, State, gStr, doGhostsO)
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        character(len=strLen), intent(in) :: gStr
        logical, optional, intent(in) :: doGhostsO

        integer :: i,j,s
        integer :: is, ie, js, je, ks, ke
        logical :: doGhosts
        real(rp) :: axT, axP
        real(rp), dimension(:,:), allocatable :: outActiveShell
        real(rp), dimension(:,:,:), allocatable :: outDen, outIntensity
        real(rp), dimension(:,:,:), allocatable :: outPrecipN, outPrecipE, outPrecipAvgE
        real(rp), dimension(:,:,:), allocatable :: outPEff  ! effective Potential

        if (present(doGhostsO)) then
            doGhosts = doGhostsO
        else
            doGhosts = .false.
        endif

        if (doGhosts) then
            is = Grid%shGrid%isg
            ie = Grid%shGrid%ieg
            js = Grid%shGrid%jsg
            je = Grid%shGrid%jeg
        else
            is = Grid%shGrid%is
            ie = Grid%shGrid%ie
            js = Grid%shGrid%js
            je = Grid%shGrid%je
        endif

        ! First, make sure root variables are there
        if (doRoot) then
            call raijuInitIO(Model, Grid, doGhosts)
            doRoot = .false.
        endif
        !Reset IO chain
        call ClearIO(IOVars)

        ! Add attributes
        call AddOutVar(IOVars,"time",State%t)
        call AddOutVar(IOVars,"MJD",State%mjd)

        ! Save axis of rotation as attributes
        if (Model%doGeoCorot) then
            call GeoAxisTP(axT, axP)
        else
            axT = 0.0
            axP = 0.0
        endif
        call AddOutVar(IOVars,"rotAxisT",axT)  ! Radians
        call AddOutVar(IOVars,"rotAxisP",axP)  ! Radians

        ! Add State variables
        call AddOutVar(IOVars,"bminX",State%Bmin(is:ie,js:je,XDIR),uStr="nT")
        call AddOutVar(IOVars,"bminY",State%Bmin(is:ie,js:je,YDIR),uStr="nT")
        call AddOutVar(IOVars,"bminZ",State%Bmin(is:ie,js:je,ZDIR),uStr="nT")

        call AddOutVar(IOVars,"eta",State%eta(is:ie,js:je, :),uStr="#/cm^3 * Rx/T") !! TODO: Maybe swap with just intensity instead

        ! Calc intensity
        allocate(outIntensity(is:ie,js:je,Grid%Nk))
        outIntensity = 0.0
        do s=1,Grid%nSpc
            if (Grid%spc(s)%flav==F_PSPH) then
                cycle  ! Skip plasmasphere since it has zero energy
            endif
            do i=is,ie
                do j=js,je
                    outIntensity(i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd) = &
                        eta2intensity(Grid%spc(s),   &
                                      State%bVol(i,j),   &
                                      State%eta (i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd))
                enddo
            enddo
        enddo
        call AddOutVar(IOVars,"intensity",outIntensity(is:ie,js:je, :),uStr="1/(s*sr*keV*cm^2)")
        !deallocate(outIntensity)
        

        ! Coupling things
        call AddOutVar(IOVars,"xmin"   ,State%xyzMin (is:ie+1,js:je+1,XDIR),uStr="Rx")
        call AddOutVar(IOVars,"ymin"   ,State%xyzMin (is:ie+1,js:je+1,YDIR),uStr="Rx")
        call AddOutVar(IOVars,"zmin"   ,State%xyzMin (is:ie+1,js:je+1,ZDIR),uStr="Rx")
        call AddOutVar(IOVars,"topo"   ,State%topo   (is:ie+1,js:je+1)*1.0_rp,uStr="0=Open, 1=Closed")
        call AddOutVar(IOVars,"colatc" ,State%thcon  (is:ie+1,js:je+1),uStr="radians")
        call AddOutVar(IOVars,"lonc"   ,State%phcon  (is:ie+1,js:je+1),uStr="radians")
        call AddOutVar(IOVars,"active" ,State%active (is:ie,js:je)*1.0_rp,uStr="-1=Inactive, 0=Buffer, 1=Active")
        call AddOutVar(IOVars,"OCBDist",State%OCBDist(is:ie,js:je)*1.0_rp,uStr="Cell distance from an open closed boundary")
        call AddOutVar(IOVars,"espot"  ,State%espot  (is:ie,js:je),uStr="kV")
        call AddOutVar(IOVars,"bVol"   ,State%bvol   (is:ie,js:je),uStr="Rx/nT")
        call AddOutVar(IOVars,"Pavg_in",State%Pavg   (is:ie,js:je, :),uStr="nPa")
        call AddOutVar(IOVars,"Davg_in",State%Davg   (is:ie,js:je, :),uStr="#/cc")

        ! Idk about you but I did not expect true to equal -1
        allocate(outActiveShell(is:ie, Grid%Nk))
        where (State%activeShells)
            outActiveShell = 1.0
        elsewhere
            outActiveShell = 0.0
        end where
        call AddOutVar(IOVars,"activeShells",outActiveShell,uStr="[Ni, Nk]")

    ! Moments
        call AddOutVar(IOVars,"Pressure",State%Press(is:ie,js:je, :),uStr="nPa")
        ! Calculate flux tube entropy using bulk pressure
        call AddOutVar(IOVars,"FTEntropy",State%Press(is:ie,js:je,1)*State%bVol(is:ie,js:je)**(5./3.),uStr="nPa*(Rp/nT)^(5/3)")
        ! Add density moment as #/cc instead of amu/cc
        allocate(outDen(is:ie,js:je,Grid%nSpc+1))
        ! Convert amu/cc to #/cc
        outDen = 0.0
        do s=1, Grid%nSpc
            outDen(:,:,s+1) = State%Den(is:ie,js:je,s+1)/Grid%spc(s)%amu
            write(*,*)"Davg_in ",s,maxval(State%Davg(is:ie,js:je,s))
            write(*,*)"Den out ",s,maxval(outDen(:,:,s+1))
            ! Don't include electrons to total number density
            if(Grid%spc(s)%spcType .ne. RAIJUELE) then
                outDen(:,:,1) = outDen(:,:,1) + outDen(:,:,s+1)
            endif
        enddo
        call AddOutVar(IOVars,"Density",outDen(is:ie,js:je, :),uStr="#/cc")
        !call AddOutVar(IOVars,"Density",State%Den,uStr="#/cc")

    ! Precipitation

        ! Calculate accumulated precipitation for each species
        allocate(outPrecipN   (is:ie,js:je,Grid%nSpc))
        allocate(outPrecipE   (is:ie,js:je,Grid%nSpc))
        allocate(outPrecipAvgE(is:ie,js:je,Grid%nSpc))
        do s=1,Grid%nSpc
            ks = Grid%spc(s)%kStart
            ke = Grid%spc(s)%kEnd
            outPrecipN(:,:,s) = sum(State%precipNFlux(is:ie,js:je,kS:kE), dim=3)
            outPrecipE(:,:,s) = sum(State%precipEFlux(is:ie,js:je,kS:kE), dim=3)

            where (outPrecipN(:,:,s) > TINY)
                outPrecipAvgE(:,:,s) = outPrecipE(:,:,s)/outPrecipN(:,:,s) * erg2kev  ! Avg E [keV]
            elsewhere
                outPrecipAvgE(:,:,s) = 0.0
            end where
        enddo
        call AddOutVar(IOVars, "precipNFlux", outPrecipN   , uStr="#/cm^2/s")
        call AddOutVar(IOVars, "precipEFlux", outPrecipE   , uStr="erg/cm^2/s")
        call AddOutVar(IOVars, "precipAvgE" , outPrecipAvgE, uStr="keV")

        
        if (Model%doFatOutput) then
            ! Calc pEffective based on current state
            ! Make full ghost size since that's what the subroutine expects
            allocate(outPEff   (Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg,Grid%Nk))
            call calcEffectivePotential(Model, Grid, State, outPEff)
            call AddOutVar(IOVars, "pEffective", outPEff(is:ie,js:je,:)*1e-3, uStr="kV")
            call AddOutVar(IOVars, "gradPotE"    , State%gradPotE    (is:ie,js:je,:), uStr="V/m")
            call AddOutVar(IOVars, "gradPotCorot", State%gradPotCorot(is:ie,js:je,:), uStr="V/m")
            call AddOutVar(IOVars, "gradVM"      , State%gradVM      (is:ie,js:je,:), uStr="V/m/lambda")
            call AddOutVar(IOVars, "cVel_th", State%cVel(is:ie,js:je,:,1), uStr="m/s")
            call AddOutVar(IOVars, "cVel_ph", State%cVel(is:ie,js:je,:,2), uStr="m/s")
            call AddOutVar(IOVars, "preciplossRates_Nk", State%lossRates  (is:ie,js:je,:), uStr="1/s")
            call AddOutVar(IOVars, "precipNFlux_Nk"    , State%precipNFlux(is:ie,js:je,:), uStr="#/cm^2/s")
            call AddOutVar(IOVars, "precipEFlux_Nk"    , State%precipEFlux(is:ie,js:je,:), uStr="erg/cm^2/s")
        endif

        call WriteVars(IOVars,.true.,Model%raijuH5, gStr)

        ! Let sub-models add stuff if they want
        if (Model%eLossWM%doOutput) then
            call eWMOutput(Model, Grid, State, gStr, doGhostsO)
        endif

    end subroutine WriteRaiju

end module raijuIO