!Various routines to handle inner magnetosphere models coupled to Voltron

module innermagsphere
    use kdefs
    use gamtypes
    use ebtypes
    use volttypes
    use gamapp
!    use sstimag
    use sstLLimag
    !use rcmimag
    use msphutils, only : RadIonosphere
    use rcmXimag
    use cmiutils, only : SquishCorners
    
    implicit none

    contains

    !Figure out which inner magnetosphere model we're using and initialize it
    !subroutine InitInnerMag(vApp,isRestart,iXML)
    subroutine InitInnerMag(vApp,gApp,iXML)
        type(voltApp_T)  , intent(inout) :: vApp
        !logical, intent(in) :: isRestart
        class(gamApp_T), intent(in) :: gApp
        type(XML_Input_T), intent(inout) :: iXML
        character(len=strLen) :: imStr

        real(rp) :: mhd_Rin = 2.0

        if (.not. vApp%doDeep) return !Why are you even here?

        call iXML%Set_Val(imStr,"coupling/imType","RAIJU")

        !NOTE: Using the fact that x2 is longitude and 2P periodic for both inner mag models
        select case (trim(toUpper(imStr)))
        case("SST","TS07")
            write(*,*)"ERROR: SST not updated for new imag plumbing."
            stop
            !vApp%imType = IMAGSST
            !vApp%prType = LPPROJ !R-phi
            !allocate(empData_T :: vApp%imagApp)
        ! case("SSTLL")  ! on lon-lat grid in the ionosphere -- like RCM
        !     vApp%imType = IMAGSSTLL
        !     vApp%prType = LLPROJ !R-phi
        !     allocate(empData_T :: vApp%imagApp)
        case("RCM")
            write(*,*)"ERROR: RCM not updated for new imag plumbing."
            stop
            !vApp%imType = IMAGRCM
            !vApp%prType = LLPROJ !Lat-lon
            !allocate(rcmIMAG_T :: vApp%imagApp)
        case("RCMX")
            write(*,*)"ERROR: RCMX not updated for new imag plumbing."
            stop
            !vApp%imType = IMAGRCMX
            !vApp%prType = LLPROJ !Lat-lon
            !allocate(rcmXIMAG_T :: vApp%imagApp)
        case("RAIJU")
            vApp%imType = IMAGRAIJU
            vApp%prType = LLPROJ
            allocate(raijuCoupler_T :: vApp%imagApp)
            allocate(imagOptions_T :: vApp%imagApp%opt)
        case DEFAULT
            write(*,*) 'Unkown imType, bailing ...'
            stop
        end select

        ! Set options
        associate(opt=>vApp%imagApp%opt, Gr=>gApp%Grid)
        opt%swF  = vApp%symh%wID
        opt%mjd0 = gApp%Model%MJD0
        opt%mhdRin  = norm2(Gr%xyz(Gr%is,Gr%js,Gr%ks,:))
        opt%mhdRinG = norm2(Gr%xyz(Gr%isg,Gr%js,Gr%ks,:)) ! Calc lowlat BC from Gamera
        opt%voltGrid = vApp%shGrid
        call iXML%Set_Val(opt%doColdStart,"/Kaiju/voltron/imag/doInit",.false.) ! Whether or not IMAG should coldStart at volt%t = 0
        end associate

        !call vApp%imagApp%doInit(iXML,gApp%Model%isRestart,vApp)
        call vApp%imagApp%InitModel(iXML)
        call vApp%imagApp%InitIO(iXML)

    end subroutine InitInnerMag

end module innermagsphere
