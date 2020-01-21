!Various routines to handle inner magnetosphere models coupled to Voltron

module innermagsphere
    use kdefs
    use gamtypes
    use ebtypes
    use volttypes
    use gamapp
    use sstimag
    use rcmimag
    
    implicit none

    !IMag eval type
    abstract interface
        subroutine IMagEval_T(x1,x2,t,imW)
            Import :: rp,NVARIMAG
            real(rp), intent(in) :: x1,x2,t
            real(rp), dimension(NVARIMAG), intent(out) :: imW
        end subroutine IMagEval_T
    end interface

    ! !K: This part is draft class for polymorphic approach to inner magnetosphere
    ! type, abstract :: InnerMag_T
    ! contains
    !     procedure(InitInnerMag_T), deferred :: InitIM
    !     procedure(AdvInnerMag_T) , deferred :: AdvIM
    !     procedure(IMagEval_T)    , deferred :: EvalIM
    !     procedure(IOInnerMag_T)  , deferred :: WriteIM
    ! end type InnerMag_T

    ! abstract interface
    !     subroutine InitInnerMag_T(vApp,isRestart,iXML)
    !         Import :: voltApp_T,XML_Input_T
    !         type(voltApp_T)  , intent(inout) :: vApp
    !         logical, intent(in) :: isRestart
    !         type(XML_Input_T), intent(inout) :: iXML
    !     end subroutine InitInnerMag_T
    ! end interface

    ! abstract interface
    !     subroutine AdvInnerMag_T(vApp,tAdv)
    !         Import :: rp,voltApp_T
    !         type(voltApp_T), intent(inout) :: vApp
    !         real(rp), intent(in) :: tAdv
    !     end subroutine AdvInnerMag_T
    ! end interface

    ! abstract interface
    !     subroutine IOInnerMag_T(vApp,nOut)
    !         Import :: voltApp_T
    !     type(voltApp_T), intent(inout) :: vApp
    !     integer, intent(in) :: nOut
    !     end subroutine IOInnerMag_T
    ! end interface

    contains

    !Figure out which inner magnetosphere model we're using and initialize it
    subroutine InitInnerMag(vApp,isRestart,iXML)
        type(voltApp_T)  , intent(inout) :: vApp
        logical, intent(in) :: isRestart
        type(XML_Input_T), intent(inout) :: iXML
        
        character(len=strLen) :: imStr

        if (.not. vApp%doDeep) return !Why are you even here?

        call iXML%Set_Val(imStr,"coupling/imType","SST")

        select case (trim(toUpper(imStr)))
        case("SST","TS07")
            vApp%imType = IMAGSST
            vApp%prType = LPPROJ !R-phi
            call InitSST(iXML,isRestart)
        case("RCM")
            vApp%imType = IMAGRCM
            vApp%prType = LLPROJ !Lat-lon
            call InitRCM(iXML,isRestart,vApp%imag2mix,vApp%time,vApp%DeepDT,vApp%IO%nRes)
        case DEFAULT
            write(*,*) 'Unkown imType, bailing ...'
            stop
        end select

    end subroutine InitInnerMag


    !Advance inner magnetosphere model to tAdv
    subroutine AdvanceInnerMag(vApp,tAdv)
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: tAdv

        if (.not. vApp%doDeep) return !Why are you even here?

        select case (vApp%imType)
        case(IMAGSST)
            call AdvanceSST(tAdv)            
        case(IMAGRCM)
            call AdvanceRCM(vApp,tAdv)
        case DEFAULT
            write(*,*) 'Unkown imType, bailing ...'
            stop
        end select

    end subroutine AdvanceInnerMag

    !Use inner mag model to prepare Gamera source terms
    subroutine InnerMag2Gamera(vApp,gApp)
        type(voltApp_T), intent(inout) :: vApp
        type(gamApp_T) , intent(inout) :: gApp

        integer :: i,j,k
        real(rp) :: x1,x2,t
        real(rp) :: imW(NVARIMAG)
        procedure(IMagEval_T), pointer :: IMagEval

        !Set evaluation routine
        IMagEval => NULL()

        select case (vApp%imType)
        case(IMAGSST)
            IMagEval => EvalSST
        case(IMAGRCM)
            IMagEval => EvalRCM
        case DEFAULT
            write(*,*) 'Unkown imType, bailing ...'
            stop
        end select

        !TODO: Think about what time to evaluate at
        t = gApp%Model%t*gApp%Model%Units%gT0

        associate(Gr=>gApp%Grid,chmp2mhd=>vApp%chmp2mhd)
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,k,x1,x2,imW)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is,Gr%is+chmp2mhd%iMax
                    x1 = chmp2mhd%xyzSquish(i,j,k,1) 
                    x2 = chmp2mhd%xyzSquish(i,j,k,2)
                    if (norm2([x1,x2])>TINY) then
                        call IMagEval(x1,x2,t,imW)
                    else
                        !Both x1/x2 are 0, projection failure
                        imW = 0.0
                    endif
                    !Assuming density/pressure coming in #/cc and nPa
                    !Lengthscale is in Rx
                    Gr%Gas0(i,j,k,:,:) = 0.0

                    Gr%Gas0(i,j,k,IMDEN ,BLK) = imW(IMDEN)
                    Gr%Gas0(i,j,k,IMPR  ,BLK) = imW(IMPR)/gApp%Model%Units%gP0
                    Gr%Gas0(i,j,k,IMX1  ,BLK) = imW(IMX1)
                    Gr%Gas0(i,j,k,IMX2  ,BLK) = imW(IMX2)

                    !Interpreting IMTSCL as units of coupling timescale
                    Gr%Gas0(i,j,k,IMTSCL,BLK) = vApp%DeepDT*imW(IMTSCL)/gApp%Model%Units%gT0
                enddo
            enddo
        enddo


        end associate
        

    end subroutine InnerMag2Gamera

    subroutine InnerMagIO(vApp,nOut)
        type(voltApp_T), intent(inout) :: vApp
        integer, intent(in) :: nOut

        select case (vApp%imType)
        case(IMAGSST)
            
        case(IMAGRCM)
            call WriteRCM(nOut,vApp%MJD,vApp%time)
        case DEFAULT
            write(*,*) 'Unkown imType, bailing ...'
            stop
        end select

    end subroutine InnerMagIO

    subroutine InnerMagRestart(vApp,nRes)
        type(voltApp_T), intent(inout) :: vApp
        integer, intent(in) :: nRes

        select case (vApp%imType)
        case(IMAGSST)
            
        case(IMAGRCM)
            call WriteRCMRestart(nRes,vApp%MJD,vApp%time)
        case DEFAULT
            write(*,*) 'Unkown imType, bailing ...'
            stop
        end select

    end subroutine InnerMagRestart
end module innermagsphere