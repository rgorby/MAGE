! Collection of data and objects for the voltron middle man

module voltapp
    use mixtypes
    use ebtypes
    use chmpdefs
    use starter
    use mhd2mix_interface
    use mix2mhd_interface
    use mhd2chmp_interface
    use chmp2mhd_interface
    use eqmap

    implicit none

    !Projection types
    enum, bind(C) 
        !L-phi (equatorial), Lat-Lon (northern ionospheric)
        enumerator :: LPPROJ,LLPROJ
    endenum

    type voltApp_T
        type(mixApp_T) :: remixApp
        type(mhd2Mix_T) :: mhd2mix
        type(mix2Mhd_T) :: mix2mhd

        type(ebTrcApp_T)  :: ebTrcApp
        type(mhd2Chmp_T)  :: mhd2chmp
        type(chmp2Mhd_T)  :: chmp2mhd

        real(rp) :: ShallowT
        real(rp) :: ShallowDT

        !Deep coupling information
        real(rp) :: DeepT
        real(rp) :: DeepDT
        logical  :: doDeep = .false. !Whether to do deep coupling
        real(rp) :: rDeep !Radius (in code units) to do deep coupling
        integer  :: iDeep = 0 !Index of max i shell containing deep coupling radius
        logical  :: doEQ = .false. !Do equatorial pressure mapping


        real(rp) :: tilt

    end type voltApp_T

    contains

    subroutine initVoltron(vApp, gApp,optFilename)
        type(gamApp_T) , intent(inout) :: gApp
        type(voltApp_T), intent(inout) :: vApp
        character(len=*), optional, intent(in) :: optFilename

        character(len=strLen) :: inpXML
        type(XML_Input_T) :: xmlInp

        if(present(optFilename)) then
            ! read from the prescribed file
            inpXML = optFilename
            call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")
        else
            !Find input deck
            call getIDeckStr(inpXML)

        endif

        write(*,*) 'Reading input deck from ', trim(inpXML)
        inquire(file=inpXML,exist=fExist)
        if (.not. fExist) then
            write(*,*) 'Error opening input deck, exiting ...'
            write(*,*) ''
            stop
        endif

        !Create XML reader
        xmlInp = New_XML_Input(trim(inpXML),'Voltron',.true.)

        vApp%tilt = 0.0_rp

    !Shallow coupling
        vApp%ShallowT = 0.0_rp
        call xmlInp%Set_Val(vApp%ShallowDT ,"coupling/dt" , 0.1_rp)

    !Deep coupling
        vApp%DeepT = 0.0_rp
        call xmlInp%Set_Val(vApp%DeepDT, "coupling/dtDeep", -1.0_rp)
        call xmlInp%Set_Val(vApp%rDeep,  "coupling/rDeep" , 10.0_rp)

        if (vApp%DeepDT>0) then
            vApp%doDeep = .true.
        else
            vApp%doDeep = .false.
        endif

        if (vApp%doDeep) then
            !Initialize deep coupling type
            !For now just assuming we're doing empirical equatorial pressure
            call initEQMap(xmlInp)
            vApp%doEQ = .true.
        endif

        if(present(optFilename)) then
            ! read from the prescribed file
            call initializeFromGamera(vApp, gApp, optFilename)
        else
            call initializeFromGamera(vApp, gApp)
        endif

    end subroutine initVoltron

    subroutine initializeFromGamera(vApp, gApp, optFilename)
        type(voltApp_T), intent(inout) :: vApp
        type(gamApp_T), intent(inout) :: gApp
        character(len=*), optional, intent(in) :: optFilename

    !Remix from Gamera
        if(present(optFilename)) then
            ! read from the prescribed file
            call init_mix(vApp%remixApp%ion,[NORTH, SOUTH],optFilename)
        else
            call init_mix(vApp%remixApp%ion,[NORTH, SOUTH])
        endif

        call init_mhd2Mix(vApp%mhd2mix, gApp, vApp%remixApp)
        call init_mix2Mhd(vApp%mix2mhd, vApp%remixApp, gApp)

    !CHIMP (TRC) from Gamera
        if (vApp%doDeep) then
            !Verify that there's some place to put deep coupling info
            if (.not. gApp%Model%doSource) then
                write(*,*) 'For deep coupling, GAMERA/source/doSouce="T" must be set!'
                write(*,*) 'Exiting ...'
                stop
            endif

            ! initialize chimp
            associate(ebTrcApp=>vApp%ebTrcApp)
            if (present(optFilename)) then
                call init_volt2Chmp(ebTrcApp,gApp,optFilename=optFilename)
            else
                call init_volt2Chmp(ebTrcApp,gApp)
            endif

            call init_mhd2Chmp(vApp%mhd2chmp, gApp, ebTrcApp)
            call init_chmp2Mhd(vApp%chmp2mhd, ebTrcApp, gApp)

            end associate
            vApp%iDeep = ShellBoundary(gApp%Model,gApp%Grid,vApp%rDeep)
        endif !doDeep

    end subroutine initializeFromGamera


!----------
!Shallow coupling stuff
    subroutine ShallowUpdate(vApp, gApp, time)
        type(gamApp_T), intent(inout) :: gApp
        type(voltApp_T), intent(inout) :: vApp
        real(rp) :: time

        ! convert gamera data to mixInput
        call Tic("G2R")
        call convertGameraToRemix(vApp%mhd2mix, gApp, vApp%remixApp)
        call Toc("G2R")

        ! run remix
        call Tic("ReMIX")
        call runRemix(vApp, time)
        call Toc("ReMIX")

        ! convert mixOutput to gamera data
        call Tic("R2G")
        call convertRemixToGamera(vApp%mix2mhd, vApp%remixApp, gApp)
        call Toc("R2G")

        vApp%ShallowT = time + vApp%ShallowDT

    end subroutine ShallowUpdate

    subroutine runRemix(vApp, time)
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: time

        ! convert gamera inputs to remix
        call mapGameraToRemix(vApp%mhd2mix, vApp%remixApp)

        ! solve for remix output
        call run_mix(vApp%remixApp%ion,vApp%tilt)

        ! get stuff from mix to gamera
        call mapRemixToGamera(vApp%mix2mhd, vApp%remixApp)

        ! output remix info
        call mix_mhd_output(vApp%remixApp%ion,vApp%mix2mhd%mixOutput,time)

    end subroutine runRemix

!----------
!Deep coupling stuff
    subroutine DeepUpdate(vApp, gApp, time)
        type(gamApp_T) , intent(inout) :: gApp
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: time

        real(rp) :: t

        t = time*gApp%Model%Units%gT0 !Time scaled back to seconds

        if (.not. vApp%doDeep) then
            !Why are you even here?
            return
        endif
        
        if (vApp%doEQ) then
            call updateEQMap(t)
        endif
        
        ! convert gamera data to chimp
        call Tic("G2C")
        call convertGameraToChimp(vApp%mhd2chmp,gApp,vApp%ebTrcApp)
        call Toc("G2C")
        ! run chimp

        call Tic("CHIMP")
        call runChimp(vApp, time)
        call Toc("CHIMP")

        ! convert chimp to gamera data
        call Tic("C2G")
        call convertChimpToGamera(vApp%chmp2mhd,vApp%ebTrcApp,gApp)
        call Toc("C2G")
        
        !Setup next coupling
        vApp%DeepT = time + vApp%DeepDT

    end subroutine DeepUpdate

    subroutine runChimp(vApp, time)
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: time

        integer :: projType
        projType = LPPROJ

        call Squish(vApp,projType)

        !Create 2D mapping to stretch from x1,x2=>xyz

    end subroutine runChimp

    !Find i-index of outer boundary of coupling domain
    function ShellBoundary(gModel,Gr,R) result(iMax)
        type(Model_T), intent(in) :: gModel
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(in) :: R
        integer :: iMax
        integer :: i
        !Just using dayside line
        !Avoiding findloc because it's not well supported by GNU
        do i=Gr%is,Gr%ie
            iMax = i
            if (Gr%xyz(i,Gr%js,Gr%ks,XDIR)>=R) exit
        enddo

    end function ShellBoundary

    subroutine Squish(vApp,projType)
        type(voltApp_T), intent(inout) :: vApp
        integer, intent(in) :: projType
        
        integer :: i,j,k
        real(rp) :: t,x1,x2
        real(rp), dimension(NDIM) :: xyz,xy0
        
        associate(ebModel=>vApp%ebTrcApp%ebModel,ebGr=>vApp%ebTrcApp%ebState%ebGr,ebState=>vApp%ebTrcApp%ebState)
        t = ebState%eb1%time
        
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,xyz,x1,x2)
        do k=ebGr%ks,ebGr%ke
            do j=ebGr%js,ebGr%je
                do i=ebGr%is,ebGr%is+vApp%iDeep
                    xyz = ebGr%xyzcc(i,j,k,:)
                    call Proj2LP(ebModel,ebState,xyz,t,x1,x2)

                    vApp%chmp2mhd%xyzSquish(i,j,k,1) = x1
                    vApp%chmp2mhd%xyzSquish(i,j,k,2) = x2

                enddo
            enddo
        enddo

        vApp%chmp2mhd%iMax = vApp%iDeep

        end associate

        contains
            subroutine Proj2LP(ebModel,ebState,xyz,t,x1,x2)
                type(chmpModel_T), intent(in) :: ebModel
                type(ebState_T)  , intent(in) :: ebState
                real(rp), dimension(NDIM), intent(in) :: xyz
                real(rp), intent(in) :: t
                real(rp), intent(out) :: x1,x2

                real(rp) :: L,phi,z
                real(rp), dimension(NDIM) :: xy0

                call getProjection(ebModel,ebState,xyz,t,xy0)
                !Map projection to L,phi

                z = abs(xy0(ZDIR))
                L = sqrt(xy0(XDIR)**2.0+xy0(YDIR)**2.0)
                phi = atan2(xy0(YDIR),xy0(XDIR))
                if (phi<0) phi = phi+2*PI

                if (z/L > 1.0e-3) then
                    !Probably failed to get to equator, set L=0
                    x1 = 0.0
                    x2 = phi
                else
                    x1 = L
                    x2 = phi
                endif
            end subroutine Proj2LP

    end subroutine Squish

    !Initialize CHIMP data structure
    subroutine init_volt2Chmp(ebTrcApp,gApp,optFilename)
        type(ebTrcApp_T), intent(inout) :: ebTrcApp
        type(gamApp_T), intent(in) :: gApp
        character(len=*), intent(in), optional     :: optFilename

        character(len=strLen) :: xmlStr
        type(XML_Input_T) :: inpXML
        
    !Create input XML object
        if (present(optFilename)) then
            xmlStr = trim(optFilename)
        else
            call getIDeckStr(xmlStr)
        endif
        inpXML = New_XML_Input(trim(xmlStr),"Chimp",.true.)

    !Initialize model
        associate(Model=>ebTrcApp%ebModel,ebState=>ebTrcApp%ebState,Gr=>gApp%Grid)
        call setUnits (Model,inpXML)
        Model%T0   = 0.0
        Model%tFin = 0.0
        Model%dt   = 0.0
        Model%t    = 0.0
        call setInterpolation(Model,inpXML)
        Model%doMHD = .true.
        call inpXML%Set_Val(Model%epsds,'tracer/epsds',1.0e-2)    
        call setBackground(Model,inpXML)

    !Initialize ebState
        !CHIMP grid is initialized from Gamera's active corners
        call ebInit_fromMHDGrid(Model,ebState,inpXML,Gr%xyz(Gr%is:Gr%ie+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1,1:NDIM))
        call InitLoc(Model,ebState%ebGr,inpXML)
        end associate

    end subroutine init_volt2Chmp

end module voltapp

