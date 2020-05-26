!Various routines to handle inner magnetosphere models coupled to Voltron

module innermagsphere
    use kdefs
    use gamtypes
    use ebtypes
    use volttypes
    use gamapp
    use sstimag
    use rcmimag
    use msphutils, only : RadIonosphere
    use cmiutils, only : SquishCorners
    
    implicit none

    contains

    !Figure out which inner magnetosphere model we're using and initialize it
    !subroutine InitInnerMag(vApp,isRestart,iXML)
    subroutine InitInnerMag(vApp,gApp,iXML)
        type(voltApp_T)  , intent(inout) :: vApp
        !logical, intent(in) :: isRestart
        type(gamApp_T), intent(in) :: gApp
        type(XML_Input_T), intent(inout) :: iXML
        real(rp) :: rad_planet_m, rad_iono_m, Rion, M0g
        character(len=strLen) :: imStr

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

        rad_planet_m = gApp%Model%Units%gx0
        !rad_iono_m = 1.01*rad_planet_m !CHANGE to get iono radius from gamera
        Rion = RadIonosphere() !Units of rp
        rad_iono_m = Rion*rad_planet_m ! m
        M0g = -gApp%Model%MagM0*gApp%Model%Units%gB0*1.e-5 ! Convert whatever units MagM0 are back to Gauss
        call vApp%imagApp%doInit(iXML,gApp%Model%isRestart,rad_planet_m,rad_iono_m,M0g,vApp)

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
        real(rp) :: imW(NVARIMAG),Qs(8)
        logical :: isTasty
        real(rp), dimension(:,:,:,:), allocatable :: SrcNC !Node-centered source terms 


    !Proceed in two steps
    ! 1) Get ingestion values at each node (cell corner)
    ! 2) Loop over cells and average from corners to cell centers

        !TODO: Think about what time to evaluate at
        t = gApp%Model%t*gApp%Model%Units%gT0
        
        associate(Gr=>gApp%Grid,chmp2mhd=>vApp%chmp2mhd)

        !Create local storage for cell corner imW's
        allocate(SrcNC(Gr%is:Gr%is+chmp2mhd%iMax+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1,1:NVARIMAG))
        SrcNC = 0.0
        chmp2mhd%isEdible = .false.
      
    ! 1) Cell corner ingestion
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,k,x1,x2,imW,isTasty)
        do k=Gr%ks,Gr%ke+1
            do j=Gr%js,Gr%je+1
                do i=Gr%is,Gr%is+chmp2mhd%iMax+1

                    if (chmp2mhd%isGood(i,j,k)) then
                        !Good projection, let's get some values
                        x1 = chmp2mhd%xyzSquish(i,j,k,1)
                        x2 = chmp2mhd%xyzSquish(i,j,k,2)
                        call vApp%imagApp%doEval(x1,x2,t,imW,isTasty)
                    else
                        !Projection wasn't good, nothing good to eat
                        imW = 0.0
                        isTasty = .false.
                        
                    endif !isGood
                    SrcNC(i,j,k,:) = imW
                    chmp2mhd%isEdible(i,j,k) = isTasty
                enddo !i loop
            enddo
        enddo


    ! 2) Corners => Centers
        Gr%Gas0 = 0.0

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,k,imW,Qs)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is,Gr%is+chmp2mhd%iMax

                    if ( all(chmp2mhd%isEdible(i:i+1,j:j+1,k:k+1)) ) then
                    !Density and pressure
                        call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMDEN),Qs)
                        imW(IMDEN) = ArithMean(Qs)
                        call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMPR) ,Qs)
                        imW(IMPR)  = ArithMean(Qs)
                    !x1 and x2
                        call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMX1),Qs)
                        imW(IMX1) = ArithMean(Qs)
                        call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMX2),Qs)
                        imW(IMX2) = CircMean(Qs)
                    !Timescale
                        call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMTSCL),Qs)
                        if ( all(Qs>TINY) ) then
                            imW(IMTSCL) = ArithMean(Qs)
                        else
                            imW(IMTSCL) = vApp%DeepDT
                        endif
                    else
                        !Not good to eat
                        imW = 0.0
                    endif

                    imW(IMTSCL) = max(imW(IMTSCL),vApp%DeepDT)

                    !Do scaling and store
                    !density/pressure coming back in #/cc and nPa
                    !ingestion timescale coming back in seconds
                    Gr%Gas0(i,j,k,IMDEN ,BLK) = imW(IMDEN)
                    Gr%Gas0(i,j,k,IMPR  ,BLK) = imW(IMPR)/gApp%Model%Units%gP0
                    Gr%Gas0(i,j,k,IMX1  ,BLK) = imW(IMX1)
                    Gr%Gas0(i,j,k,IMX2  ,BLK) = imW(IMX2)
                    Gr%Gas0(i,j,k,IMTSCL,BLK) = imW(IMTSCL)/gApp%Model%Units%gT0
                
                enddo !i loop
            enddo
        enddo

    !Now do some touch up at the axis and get outta here

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

end module innermagsphere
