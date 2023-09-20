! implementation of voltron coupling between gamera and other models

module gamCouple
    use gamtypes
    use step
    use init
    use mhdgroup
    use output

    implicit none

    integer, parameter :: mix2mhd_varn = 1  ! for now just the potential is sent back

    type, extends(gamApp_T) :: gamCoupler_T

        ! data for coupling remix to gamera
        real(rp), dimension(:,:,:,:,:), allocatable :: mixOutput
        real(rp), dimension(:,:,:), allocatable :: gPsi
        type(Map_T), allocatable, dimension(:,:) :: PsiMaps
        real(rp) :: rm2g
        real(rp) :: Rion

        ! data for coupling imag to gamera
        real(rp), dimension(:,:,:,:), allocatable :: SrcNC !Node-centered source terms

        contains

        ! only over-riding specific functions
        ! actually none for local non-mpi gamera?
        !procedure :: InitModel => gamCplInitModel
        !procedure :: InitIO => gamInitIO
        procedure :: WriteRestart => gamCplWriteRestart
        procedure :: ReadRestart => gamCplReadRestart
        !procedure :: WriteConsoleOutput => gamWriteConsoleOutput
        !procedure :: WriteFileOutput => gamWriteFileOutput
        !procedure :: WriteSlimFileOutput => gamWriteSlimFileOutput
        !procedure :: AdvanceModel => gamCplAdvanceModel

        ! add new coupling function which can be over-ridden by children
        procedure :: InitCoupler => gamInitCoupler
        procedure :: CoupleModel => gamCoupleModel

    end type gamCoupler_T

    contains

    ! procedures for gamCoupler_T

    subroutine gamCplWriteRestart(App)
        class(gamCoupler_T), intent(inout) :: App

        ! write parent's restart data
        call gamWriteRestart(App)

        ! then my own
        call writeGamCouplerRestart(App)

    end subroutine

    subroutine gamCplReadRestart(App, resId, nRes)
        class(gamCoupler_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        ! read parent's restart data
        call gamReadRestart(App, resId, nRes)

        ! then my own
        call readGamCouplerRestart(App, resId, nRes)

    end subroutine

    subroutine gamInitCoupler(App, voltApp)
        class(gamCoupler_T), intent(inout) :: App
        class(voltApp_T), intent(in) :: voltApp

        call init_mix2MhdCoupler(App, voltApp%remixApp)

        allocate(SrcNC(App%Grid%is:voltApp%chmp2mhd%iMax+1,App%Grid%js:App%Grid%je+1,App%Grid%ks:App%Grid%ke+1,1:NVARIMAG))

    end subroutine

    subroutine gamCoupleModel(App, voltApp)
        class(gamCoupler_T), intent(inout) :: App
        class(voltApp_T), intent(in) :: voltApp

        call mapRemixToGamera(App, voltApp%remixApp)
        call convertRemixToGamera(App, voltApp%remixApp)

        call InnerMag2Gamera(App, voltApp)

    end subroutine gamCoupleModel

    ! actual procedures

    subroutine writeGamCouplerRestart(App)
        class(gamCoupler_T), intent(inout) :: App

        character(len=strLen) :: ResF
        integer, parameter :: MAXGCIOVAR = 20
        type(IOVAR_T), dimension(MAXGCIOVAR) :: IOVars

        !Restart Filename
        write (ResF, '(A,A,I0.5,A)') trim(App%Model%RunID), ".mix2gam.Res.", Model%IO%nRes, ".h5"

        !Reset IO chain
        call ClearIO(IOVars)

        !Remix Coupling Variables
        call AddOutVar(IOVars,"mixOutput",App%mixOutput)
        call AddOutVar(IOVars,"gPsi"     ,App%gPsi)

        !Imag Coupling Variables
        call AddOutVar(IOVars,"SrcNC"    ,App%SrcNC)

        !Write out, force real precision
        call WriteVars(IOVars,.false.,ResF)

    end subroutine

    subroutine readGamCouplerRestart(App, resId, nRes)
        class(gamCoupler_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        character(len=strLen) :: ResF
        logical :: fExist
        integer, parameter :: MAXGCIOVAR = 20
        type(IOVAR_T), dimension(MAXGCIOVAR) :: IOVars

        !Restart Filename
        write (ResF, '(A,A,I0.5,A)') trim(resId), ".mix2gam.Res.", nRes, ".h5"

        inquire(file=ResF,exist=fExist)
        if (.not. fExist) then
            !Error out and leave
            write(*,*) 'Unable to open input restart file for mix2gam, exiting'
            stop
        endif

        !Reset IO chain
        call ClearIO(IOVars)

        !Read Remix Coupling Variables
        call AddInVar(IOVars,"mixOutput")
        call AddInVar(IOVars,"gPsi")

        !Read Imag Coupling Variables
        call AddInVar(IOVars,"SrcNC")

        !Get data
        call ReadVars(IOVars,.false.,ResF)

        call IOArray5DFill(IOVars,"mixOutput",App%mixOutput)
        call IOArray3DFill(IOVars,"gPsi",App%gPsi)
        call IOArray4DFill(IOVars,"SrcNC",App%SrcNC)

    end subroutine

     subroutine init_mix2MhdCoupler(gameraApp, remixApp)
        class(gamCoupler_T), intent(inout) :: gameraApp
        class(mixApp_T), intent(in) :: remixApp

        integer :: i,j,k,iG,h,l
        real(rp) :: mhd_Rin
        real(rp), allocatable, dimension(:,:,:,:,:) :: mhdPsiGrid
        real(rp), allocatable, dimension(:,:) :: mhdt, mhdp, mhdtFpd, mhdpFpd
        type(mixGrid_T) :: mhdG
        type(Map_T) :: Map
        real(rp) :: gB0,gv0,gx0
        logical :: isRestart

        gB0 = gameraApp%Model%Units%gB0
        gv0 = gameraApp%Model%Units%gv0
        gx0 = gameraApp%Model%Units%gx0

        gameraApp%rm2g = gB0*gV0*gx0*1.0e-12 !Scaling factor for remix potential [kV]
        gameraApp%Rion = RadIonosphere()

        ! allocate remix arrays
        allocate(gamerApp%gPsi(1:PsiSh+1,gameraApp%Grid%js:gameraApp%Grid%je+1,gameraApp%Grid%ks:gameraApp%Grid%ke+1))
        allocate(mhdPsiGrid(1:PsiSh+1, gameraApp%Grid%js:gameraApp%Grid%je+1, gameraApp%Grid%ks:gameraApp%Grid%ke/2+1, 1:3, 1:2))
        allocate(gameraApp%mixOutput(1:PsiSh+1, gameraApp%Grid%js:gameraApp%Grid%je+1, gameraApp%Grid%ks:gameraApp%Grid%ke/2+1, 1:mix2mhd_varn, 1:2))
        allocate(gameraApp%PsiMaps(PsiSh,size(remixApp%ion)))
        gameraApp%gPsi = 0.0

        ! get those grid coordinates (corner centers for Psi)
        do k=gameraApp%Grid%ks,gameraApp%Grid%ke+1
            do j=gameraApp%Grid%js,gameraApp%Grid%je+1
                ! note, PsiShells give shell numbers based on cell centers per our
                ! convenion
                ! thus, no -1 below
                do i=1,PsiSh+1
                    iG = PsiSt+i-1

                    ! note conversion to Rion units which are expected on the remix
                    ! side
                    if (k<=gameraApp%Grid%ke/2+1) then
                        mhdPsiGrid(i,j,k,:,NORTH) = gameraApp%Grid%xyz(iG,j,k,:)/gameraApp%Rion
                    endif
                    if (k>=gameraApp%Grid%ke/2+1) then
                        mhdPsiGrid(i,j,k-gameraApp%Grid%ke/2,:,SOUTH) = gameraApp%Grid%xyz(iG,j,k,:)/gameraApp%Rion
                    endif
                enddo
            enddo
        enddo

        do h=1,size(remixApp%ion)
            ! set up interpolation map(s) for mix2mhd
            do l=1,PsiSh
                call mix_mhd_grid(mhdPsiGrid(l,:,:,:,h),mhdt,mhdp,mhdtFpd,mhdpFpd,mhd_Rin)
                call init_grid_fromTP(mhdG,mhdt,mhdp,isSolverGrid=.false.)
                call mix_set_map(remixApp%ion(h)%G,mhdG,Map)
                gameraApp%PsiMaps(l,h) = Map
            enddo
        enddo

        !Initialize the mixOutput mapping
        gameraApp%mixOutput = 0.0
        isRestart = gameraApp%Model%isRestart
        if (isRestart) then
            !We have data from mix restart file
            write(*,*) "mapRemixToGamera"
            call mapRemixToGamera(mix2mhd, remixApp)
        end if

        deallocate(mhdPsiGrid, mhdt, mhdp, mhdtFpd, mhdpFpd)

    end subroutine init_mix2Mhd

    subroutine mapRemixToGamera(gameraApp, remixApp)
        class(gamCoupler_T), intent(inout) :: gameraApp
        class(mixApp_T), intent(inout) :: remixApp
        integer :: l,h ! hemisphere
        integer :: v ! mhd var

        real(rp), dimension(:,:), allocatable, save :: gPsi_tmp  ! gamera potential

        ! map to potential shells
        do h=1,size(remixApp%ion)
            do v=1,size(gameraApp%mixOutput,4)  ! mirroring mix2mhd here, but for now only one variable
                do l=1,PsiSh
                    call mix_map_grids(gameraApp%PsiMaps(l,h),remixApp%ion(h)%St%Vars(:,:,POT),gPsi_tmp)
                    ! this is going back to gamera
                    select case (v)
                    case (MHDPSI)
                        gameraApp%mixOutput(l,:,:,v,h) = gPsi_tmp   ! north
                    end select
                enddo
            enddo
        enddo
    end subroutine mapRemixToGamera

    subroutine convertRemixToGamera(gameraApp, remixApp, doCorotO)
        class(gamCoupler_T), intent(inout) :: gameraApp
        class(mixApp_T), intent(inout) :: remixApp
        logical, intent(in), optional :: doCorotO

        ! convert the "remixOutputs" variable to inEijk and inExyz, which are in
        ! Gamera coordinates
        integer :: i,nbc
        logical :: doCorot

         ! populate potential on gamera grid
         gameraApp%gPsi = 0.0
         do i=1,PsiSh+1
            gameraApp%gPsi(i,:,gameraApp%Grid%ks:gameraApp%Grid%ke/2+1)   = gameraApp%mixOutput(i,:,:,MHDPSI,NORTH)
            gameraApp%gPsi(i,:,gameraApp%Grid%ke/2+1:gameraApp%Grid%ke+1) = gameraApp%mixOutput(i,:,:,MHDPSI,SOUTH)
         enddo

        if (present(doCorotO)) then
          doCorot = doCorotO
        else
          doCorot = .true.
        endif

        ! add corotation
        if (doCorot) call CorotationPot(gameraApp%Model, gameraApp%Grid, gameraApp%gPsi)

        ! find the remix BC to write data into
        nbc = FindBC(gameraApp%Model,gameraApp%Grid,INI)
        SELECT type(iiBC=>gameraApp%Grid%externalBCs(nbc)%p)
            TYPE IS (IonInnerBC_T)
                call Ion2MHD(gameraApp%Model,gameraApp%Grid,gameraApp%gPsi,iiBC%inEijk,iiBC%inExyz,gameraApp%rm2g)
            CLASS DEFAULT
                write(*,*) 'Could not find Ion Inner BC in remix IC'
                stop
        END SELECT

    end subroutine convertRemixToGamera

    subroutine convertImagToGamera(gApp, vApp)
        class(gamCoupler_T), intent(inout) :: gApp
        type(voltApp_T), intent(inout) :: vApp

        integer :: i,j,k,Nk
        real(rp) :: x1,x2,t
        real(rp) :: imW(NVARIMAG),Qs(8)
        logical :: isTasty

    !Proceed in two steps
    ! 1) Get ingestion values at each node (cell corner)
    ! 2) Loop over cells and average from corners to cell centers

        !TODO: Think about what time to evaluate at
        t = gApp%Model%t*gApp%Model%Units%gT0

        associate(Gr=>gApp%Grid,chmp2mhd=>vApp%chmp2mhd,SrcNC=>gApp%SrcNC)

        !Create local storage for cell corner imW's
        SrcNC = 0.0
        chmp2mhd%isEdible = .false.

    ! 1) Cell corner ingestion
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,k,x1,x2,imW,isTasty)
        do k=Gr%ks,Gr%ke+1
            do j=Gr%js,Gr%je+1
                do i=Gr%is,chmp2mhd%iMax+1

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
                do i=Gr%is,chmp2mhd%iMax

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
                        imW(IMX2) = CircMeanDeg(Qs)

                    !Timescale
                        call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMTSCL),Qs)
                        if ( all(Qs>TINY) ) then
                            imW(IMTSCL) = ArithMean(Qs)
                        else if (any(Qs>TINY)) then
                            imW(IMTSCL) = sum(Qs,mask=(Qs>TINY))/count(Qs>TINY)
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
        do i=Gr%is,chmp2mhd%iMax
            !+X pole
            imW(IMDEN) = AvgOverGood(Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMDEN,BLK),Nk)
            imW(IMPR ) = AvgOverGood(Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMPR ,BLK),Nk)
            Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMDEN,BLK) = imW(IMDEN)
            Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMPR ,BLK) = imW(IMPR )

            !-X pole
            imW(IMDEN) = AvgOverGood(Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMDEN,BLK),Nk)
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

    end subroutine

end module

