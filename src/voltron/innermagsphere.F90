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
        !Mapping of cell corners
        !x12C = chmp2mhd%xyzSquish(i:i+1,j:j+1,k:k+1,1:2) 
        subroutine IMagEval_T(imagApp,x1,x2,x12C,t,imW)
            Import :: rp,NVARIMAG
            class(*), intent(inout) :: imagApp
            real(rp), intent(in) :: x1,x2,t
            real(rp), intent(in) :: x12C(2,2,2,2)
            real(rp), dimension(NVARIMAG), intent(out) :: imW
        end subroutine IMagEval_T

    end interface


    contains

    !Figure out which inner magnetosphere model we're using and initialize it
    subroutine InitInnerMag(vApp,gApp,iXML)
        type(voltApp_T)  , intent(inout) :: vApp
        type(gamApp_T), intent(in) :: gApp
        type(XML_Input_T), intent(inout) :: iXML
        character(len=strLen) :: imStr
        real(rp) :: iono_rad_m

        if (.not. vApp%doDeep) return !Why are you even here?

        call iXML%Set_Val(imStr,"coupling/imType","SST")

        !NOTE: Using the fact that x2 is longitude and 2P periodic for both inner mag models
        select case (trim(toUpper(imStr)))
        case("SST","TS07")
            vApp%imType = IMAGSST
            vApp%prType = LPPROJ !R-phi
            allocate(eqData_T :: vApp%imagApp)
        case("RCM")
            vApp%imType = IMAGRCM
            vApp%prType = LLPROJ !Lat-lon
            allocate(rcmIMAG_T :: vApp%imagApp)
        case DEFAULT
            write(*,*) 'Unkown imType, bailing ...'
            stop
        end select

        !iono_rad_m = gApp%Model%Units%Rion*gApp%Model%Units%gx0 !Convert Rp to m
        !Assume for now that ionosphere is always 1.01*planet radius
        iono_rad_m = 1.01*gApp%Model%Units%gx0 ! m
        call vApp%imagApp%doInit(iXML,gApp%Model%isRestart,gApp%Model%Units%gx0,iono_rad_m,vApp)

    end subroutine InitInnerMag


    !Advance inner magnetosphere model to tAdv
    subroutine AdvanceInnerMag(vApp,tAdv)
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: tAdv

        if (.not. vApp%doDeep) return !Why are you even here?

        call vApp%imagApp%doAdvance(vApp,tAdv)

    end subroutine AdvanceInnerMag

    !Use inner mag model to prepare Gamera source terms
    subroutine InnerMag2Gamera(vApp,gApp)
        type(voltApp_T), intent(inout) :: vApp
        type(gamApp_T) , intent(inout) :: gApp

        integer :: i,j,k,Nk
        real(rp) :: x1,x2,t
        real(rp) :: imW(NVARIMAG)
        real(rp) :: x12C(2,2,2,2)
        real(rp) :: xAng(8)
        real(rp) :: xMag

        !TODO: Think about what time to evaluate at
        t = gApp%Model%t*gApp%Model%Units%gT0

        associate(Gr=>gApp%Grid,chmp2mhd=>vApp%chmp2mhd)
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,k,x1,x2,imW,x12C,xMag,xAng)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is,Gr%is+chmp2mhd%iMax
                    x12C = chmp2mhd%xyzSquish(i:i+1,j:j+1,k:k+1,1:2)
                    xMag = minval( norm2(x12C,dim=4) )
                    if (xMag > TINY) then
                        !All projected corners are good
                        x1 = sum(x12C(:,:,:,1))/8.0
                        xAng = reshape( x12C(1:2,1:2,1:2,2), [8] )
                        x2 = CircMean(xAng)
                        call vApp%imagApp%doEval(x1,x2,x12C,t,imW)
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

                    !Use IMTSCL if set, otherwise set to coupling timescale
                    if (imW(IMTSCL) > TINY) then
                        Gr%Gas0(i,j,k,IMTSCL,BLK) = max(imW(IMTSCL),vApp%DeepDT)/gApp%Model%Units%gT0
                    else
                        Gr%Gas0(i,j,k,IMTSCL,BLK) = vApp%DeepDT/gApp%Model%Units%gT0
                    endif
                    
                enddo
            enddo
        enddo

        !Do averaging for first cell next to singularity
        !Do for +/- X pole and density/pressure
        Nk = Gr%ke-Gr%ks+1
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,imW)
        do i=Gr%is,Gr%is+chmp2mhd%iMax
            !+X pole
            imW(IMDEN) = AvgOverGood(Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMDEN,BLK),Nk)
            imW(IMPR ) = AvgOverGood(Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMPR ,BLK),Nk)
            Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMDEN,BLK) = imW(IMDEN)
            Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMPR ,BLK) = imW(IMPR )

            !-X pole
            imW(IMDEN) = AvgOverGood(Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMDEN,BLK),Nk)
            imW(IMPR ) = AvgOverGood(Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMPR ,BLK),Nk)
            Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMDEN,BLK) = imW(IMDEN)
            Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMPR ,BLK) = imW(IMPR )
        enddo

        end associate

        contains
            function AvgOverGood(Q,Nk) result(Qavg)
                real(rp), intent(in), dimension(Nk) :: Q
                integer , intent(in) :: Nk

                real(rp) :: Qavg
                integer :: Nkg

                if ( any(Q>TINY) ) then
                    Nkg = count(Q>TINY)
                    Qavg = sum(Q,mask=(Q>TINY))/Nkg
                else
                    Qavg = 0.0
                endif

            end function AvgOverGood

    end subroutine InnerMag2Gamera

    subroutine InnerMagIO(vApp,nOut)
        type(voltApp_T), intent(inout) :: vApp
        integer, intent(in) :: nOut

        call vApp%imagApp%doIO(nOut,vApp%MJD,vApp%time)

    end subroutine InnerMagIO

    subroutine InnerMagRestart(vApp,nRes)
        type(voltApp_T), intent(inout) :: vApp
        integer, intent(in) :: nRes

        call vApp%imagApp%doRestart(nRes,vApp%MJD,vApp%time)

    end subroutine InnerMagRestart

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

end module innermagsphere
