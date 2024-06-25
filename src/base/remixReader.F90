!> Reads REMIX data from output files, performs interpolation in time, etc.
module remixReader

    use kdefs
    use ioH5
    use XML_Input
    use iotable
    use shellGrid
    use shellUtils

    use hdf5

    implicit none

    integer :: MAXIOVARS = 20

!------
! Types
!------
    type rmHemi_T
        integer :: nStp
        real(rp) :: time
        type(ShellGrid_T) :: shGr
            !! Copy of ShellGrid as defined in rmReader_T
        type(ShellGridVar_T) :: Fac,SigP,SigH,Pot
            !! Vars defined relative to ShellGrid
    end type rmHemi_T


    type rmReader_T
        ! Model stuff
        character(len=strLen) :: rmF
            !! Remix filename we are reading form
        
        ! Grid stuff
        real(rp), dimension(:,:,:), allocatable :: XY
            !! X/Y coordinates in 2D. Convenient since we need it for lots of 2D arrays
        type(ShellGrid_T) :: shGr
        ! Its kinda like a grid of file steps if you really think about it
        type(ioTab_T) :: rmTab
            !! Table of steps in a mix.h5 file

        ! State stuff
        real(rp) :: time
            !! Current sim time according to whoever's in charge of rmReader
        integer :: i1=-1,i2=-1
            !! Input file step numbers bracketing time
        type(rmHemi_T) :: rmN1,rmN2,rmS1,rmS2
            !! File data at step numbers bracketing time
        logical :: doStatic
            !! Whether or not we are out of valid remix time bounds and should switch to static operation
        type(ShellGridVar_T), dimension(NORTH:SOUTH) :: nsFac, nsSigP, nsSigH, nsPot
            !! Time-interpolated North/South hemisphere shellVar objects
    end type rmReader_T


    contains

!------
! Init
!------

    subroutine initRM(ftag,inpXML,rmReader)
        character(len=*), intent(in) :: ftag
            !! Filename tag in this format: <ftag>.mix.h5
        type(XML_Input_T), intent(in) :: inpXML
        type(rmReader_T), intent(inout) :: rmReader

        type(IOVAR_T), dimension(MAXIOVARS) :: IOVars
        real(rp), dimension(:), allocatable :: th1D, ph1D
            !! 1D arrays of theta and phi coordinates, derived from XY locations
        real(rp), dimension(:,:,:), allocatable :: tmpXY
            !! Remix data is coming in at (Np,Nt), we use this to read and then convert to (Nt,Np)
        real(rp), dimension(:,:), allocatable :: tmpXcn, tmpYcn
            !! Temporary X/Y in real mix corner format
        integer :: Nt,Np
            !! # of non-ghost cell centers
        integer :: j, h
            !! Loop indices
        real(rp) :: Rp_m, Ri_rp
        
        write(rmReader%rmF,'(2a)') trim(adjustl(ftag)),'.mix.h5'
        write(*,*) 'Initializing rmReader w/ ', trim(rmReader%rmF)
        
        ! Get file data loaded in
        call InitIOTab(rmReader%rmTab, inpXML, rmReader%rmF)
        rmReader%rmTab%bStr = trim(adjustl(ftag))
        rmReader%doStatic = (rmReader%rmTab%N == 1)

        ! Start with grid
        call ClearIO(IOVars)
        call AddInVar(IOVars,"X")
        call AddInVar(IOVars,"Y")
        call AddInVar(IOVars,"Ri_m")  ! Attribute
        call AddInVar(IOVars,"Rp_m")  ! Attribute
        call ReadVars(IOVars,.true.,rmReader%rmF)


        !! Grid conversion
        !! X/Y in file are grid corners around the data
        !! To get original X/Y coords, we need to take the average and retrieve the cell-centered coordinates
        !! This means that the variable data is now stored at cell corners
        !! But, the outputted variable data only goes from phi [0,2Pi-dphi]
        !! So in order to comply with shellGrid expectations, we need to copy first phi column to end of array as well
        !! So the final size of the arrays we care about are (Nt-1,Np)
        !! Also, keep in mind that mix.h5 drops the pole theta points, so we have 1 less theta than run-time mix. Doesn't matter here though
        Np = IOVars(FindIO(IOVars, "X"))%dims(1)
        Nt = IOVars(FindIO(IOVars, "X"))%dims(2)

        ! Get ionosphere radius in units of Rp if possible
        if (ioExist(rmReader%rmF, "Rp_m") .and. ioExist(rmReader%rmF, "Ri_m")) then
            Rp_m  = IOVars(FindIO(IOVars, "Rp_m"))%data(1)
            Ri_rp = IOVars(FindIO(IOVars, "Ri_m"))%data(1) / Rp_m
        else
            ! If info not present in file, default to hard-coded Earth values
            Ri_rp = RIonE*1e6/REarth
        endif

        ! One day we will get grid in this format directly from the mix.h5 file
        ! But today is not that day
        allocate(rmReader%XY(Nt,Np,XDIR:YDIR))
        allocate(tmpXY(Np,Nt,XDIR:YDIR))
        call IOArray2DFill(IOVars,"X",tmpXY(:,:,XDIR))
        call IOArray2DFill(IOVars,"Y",tmpXY(:,:,YDIR))
        call genInGrid(tmpXY(:,:,XDIR), tmpXY(:,:,YDIR), tmpXcn, tmpYcn)
        rmReader%XY(2:Nt,:Np-1,XDIR) = transpose(tmpXcn)
        rmReader%XY(2:Nt,:Np-1,YDIR) = transpose(tmpYcn)

        ! Wrap in j
        rmReader%XY(:,Np,:) = rmReader%XY(:,1,:)

        ! Fill the pole
        rmReader%XY(1,:,:) = 0.0
        
        ! Now XY is in the desired format, time for shellGrid

        ! XY are in units of ionospheric radii (r = 1 Ri)
        ! So the theta and phi we calculate are also for the ionospheric grid
        allocate(th1D(Nt))
        allocate(ph1D(Np))
        th1D = asin(rmReader%XY(:,1,XDIR))
        ph1D = atan2(rmReader%XY(Nt-1,:,YDIR), rmReader%XY(Nt-1,:,XDIR))  ! note, first index doesn't matter here

        ! Clean up phi for shellGrid generation
        do j=Np/2-1,Np
            if (ph1D(j) < 0) then
                ph1D(j) = ph1D(j) + 2.0_rp*PI
            endif
        enddo
        if (abs(ph1D(1)) < TINY) then
            ph1D(1) = 0.0
        endif
        if (ph1D(Np) < PI) then
            ph1D(Np) = ph1D(Np) + 2.0_rp*PI
        endif


        associate(sh=>rmReader%shGr)
        
        call GenShellGrid(sh, th1D, ph1D, "remixReader", radO=Ri_rp)

        ! Hooray we have a shellGrid now
        ! Init our vars
        do h=NORTH,SOUTH
            call initShellVar(sh, SHGR_CORNER, rmReader%nsFac(h))
            rmReader%nsFac(h)%mask(sh%is:sh%ie+1, sh%js:sh%je+1) = .true.
            call initShellVar(sh, SHGR_CORNER, rmReader%nsPot(h) , rmReader%nsFac(h)%mask)
            call initShellVar(sh, SHGR_CORNER, rmReader%nsSigP(h), rmReader%nsFac(h)%mask)
            call initShellVar(sh, SHGR_CORNER, rmReader%nsSigH(h), rmReader%nsFac(h)%mask)
        enddo
        
        ! Now init hemispheres
        call initHemi(rmReader%rmN1, sh, "_N1")
        call initHemi(rmReader%rmN2, sh, "_N2")
        call initHemi(rmReader%rmS1, sh, "_S1")
        call initHemi(rmReader%rmS2, sh, "_S2")

        end associate


        contains


        subroutine initHemi(rmHemi,shGr, nameSuffix)
            type(rmHemi_T), intent(inout) :: rmHemi
            type(ShellGrid_T), intent(in) :: shGr
                !! Reference shellGrid, already defined
            character(len=*), intent(in) :: nameSuffix

            character(len=strLen) :: fullName

            rmHemi%nStp = -1 !Not yet set
            rmHemi%time = 0.0
            
            associate(hsg=>rmHemi%shGr)
            
            ! Generate our own copy of the parent shellGrid
            write(fullName,"(A,A)") trim(shGr%name), nameSuffix
            call GenChildShellGrid(shGr, hsg, fullName)

            ! Init our variables
            call initShellVar(hsg, SHGR_CORNER, rmHemi%Fac)
            rmHemi%Fac%mask(hsg%is:hsg%ie+1, hsg%js:hsg%je+1) = .true.
            call initShellVar(hsg, SHGR_CORNER, rmHemi%Pot , rmHemi%Fac%mask) 
            call initShellVar(hsg, SHGR_CORNER, rmHemi%SigP, rmHemi%Fac%mask)
            call initShellVar(hsg, SHGR_CORNER, rmHemi%SigH, rmHemi%Fac%mask)
            
            end associate

        end subroutine initHemi


    end subroutine initRM


!------
! Update
!------

    subroutine updateRM(rmReader, t)
        ! Update rmReader to time t
        ! TODO: This should do fancier stuff, like ebICstd:updateFields
        type(rmReader_T), intent(inout) :: rmReader
        real(rp), intent(in) :: t

        real(rp) :: w1, w2
        
        ! Update tabSlices
        call GetTabSlc(rmReader%rmTab,t,rmReader%i1,rmReader%i2)
        if ( t >= maxval(rmReader%rmTab%times) ) then
            rmReader%doStatic = .true.
        endif

        ! Read the 4 hemispheres
        call readHemi(rmReader, rmReader%rmN1, rmReader%i1, NORTH)
        call readHemi(rmReader, rmReader%rmS1, rmReader%i1, SOUTH)

        call readHemi(rmReader, rmReader%rmN2, rmReader%i2, NORTH)
        call readHemi(rmReader, rmReader%rmS2, rmReader%i2, SOUTH)

        !Now fill in remix main state for this t
        call GetTWgts(rmReader%rmN1%time, rmReader%rmN2%time, t, rmReader%doStatic, w1, w2)
        call hemi2rm(rmReader, w1, w2)

        !write(*,*)"-----"
        !write(*,*)t
        !write(*,*)rmReader%rmN1%time,",",rmReader%rmN2%time
        !write(*,*)w1,",",w2

        contains

        subroutine readHemi(rmReader, rmHemi, nStp, nsID)
            !! Reads hemisphere information from remix file step
            type(rmReader_T), intent(in) :: rmReader
            type(rmHemi_T), intent(inout) :: rmHemi
            integer, intent(in) :: nStp, nsID

            character(len=strLen) :: hID,gStr
            type(IOVAR_T), dimension(MAXIOVARS) :: IOVars

            ! Check to see if we need to read
            if (rmHemi%nStp == nStp) then
                ! We've already read this data
                return
            endif
            ! Otherwise, get the date
            ! Which hemisphere?
            if (nsID == NORTH) then
                hID = "NORTH"
            else
                hID = "SOUTH"
            endif

            gStr = trim(rmReader%rmTab%gStrs(nStp))
            write(*,'(5a)') '<Reading hemisphere from ', trim(rmReader%rmTab%bStr), '/', trim(gStr), '>'

            rmHemi%time = rmReader%rmTab%times(nStp)
            rmHemi%nStp = nStp

            call ClearIO(IOVars)
            call AddInVar(IOVars,"Field-aligned current " // hID)
            call AddInVar(IOVars, "Pedersen conductance " // hID)
            call AddInVar(IOVars,     "Hall conductance " // hID)
            call AddInVar(IOVars,            "Potential " // hID)
            call ReadVars(IOVars,.true.,rmReader%rmF,gStr)

            ! Abstract away the janky mapping, eventually replace with ShellGrid-friendly variable reading
            call readVarJank(IOVars, "Field-aligned current " // hID, rmHemi%shGr, rmHemi%Fac )
            call readVarJank(IOVars,  "Pedersen conductance " // hID, rmHemi%shGr, rmHemi%SigP)
            call readVarJank(IOVars,      "Hall conductance " // hID, rmHemi%shGr, rmHemi%SigH)
            call readVarJank(IOVars,             "Potential " // hID, rmHemi%shGr, rmHemi%Pot )

        end subroutine readHemi

    end subroutine updateRM


!------
! Helpers
!------

    subroutine GetTWgts(t1, t2, t, doStatic, w1, w2)
        real(rp), intent(in)  :: t1, t2, t
        logical, intent(in) :: doStatic
        real(rp), intent(out) :: w1,w2

        real(rp) :: dt

        if (doStatic) then
            w1 = 1.0
            w2 = 0.0
        else
            dt = t2-t1
            w1 = (t2 -  t)/dt
            w2 = (t  - t1)/dt
        endif
    end subroutine GetTWgts


    subroutine hemi2rm(rmReader, w1, w2)
        type(rmReader_T), intent(inout) :: rmReader
        real(rp), intent(in) :: w1, w2

        rmReader%nsFac (NORTH)%data = w1*rmReader%rmN1%Fac %data + w2*rmReader%rmN2%Fac %data
        rmReader%nsPot (NORTH)%data = w1*rmReader%rmN1%Pot %data + w2*rmReader%rmN2%Pot %data
        rmReader%nsSigP(NORTH)%data = w1*rmReader%rmN1%SigP%data + w2*rmReader%rmN2%SigP%data
        rmReader%nsSigH(NORTH)%data = w1*rmReader%rmN1%SigH%data + w2*rmReader%rmN2%SigH%data
        rmReader%nsFac (SOUTH)%data = w1*rmReader%rmS1%Fac %data + w2*rmReader%rmS2%Fac %data
        rmReader%nsPot (SOUTH)%data = w1*rmReader%rmS1%Pot %data + w2*rmReader%rmS2%Pot %data
        rmReader%nsSigP(SOUTH)%data = w1*rmReader%rmS1%SigP%data + w2*rmReader%rmS2%SigP%data
        rmReader%nsSigH(SOUTH)%data = w1*rmReader%rmS1%SigH%data + w2*rmReader%rmS2%SigH%data
    end subroutine hemi2rm

!------
! Temporary Helpers
!------

    subroutine genInGrid(xc,yc,x,y)
        !! Mix h5 grid to "real" mix grid
        !! NOTE: This has been modified from the version in mixio
        real(rp), dimension(:,:),intent(in) :: xc,yc ! with corners 1/2-cell shifted from original    
        real(rp), dimension(:,:),allocatable,intent(out) :: x,y
    
        integer, dimension(2) :: dims
        integer :: Np, Nt
    
        dims = shape(xc)
        Np = dims(1)-1; Nt = dims(2)
    
        if (.not.allocated(x)) allocate(x(Np,Nt-1))
        if (.not.allocated(y)) allocate(y(Np,Nt-1))
    
        !x(:,2:Nt) = 0.25*(xc(2:Np+1,2:Nt)+xc(1:Np,2:Nt)+xc(2:Np+1,1:Nt-1)+xc(1:Np,1:Nt-1))
        !y(:,2:Nt) = 0.25*(yc(2:Np+1,2:Nt)+yc(1:Np,2:Nt)+yc(2:Np+1,1:Nt-1)+yc(1:Np,1:Nt-1))
        x = 0.25*(xc(2:Np+1,2:Nt)+xc(1:Np,2:Nt)+xc(2:Np+1,1:Nt-1)+xc(1:Np,1:Nt-1))
        y = 0.25*(yc(2:Np+1,2:Nt)+yc(1:Np,2:Nt)+yc(2:Np+1,1:Nt-1)+yc(1:Np,1:Nt-1))
    
        ! fix pole
        !x(:,1) = 0
        !y(:,1) = 0
      end subroutine genInGrid


      subroutine readVarJank(IOV, vName, shGr, v)
        !! Placeholder readVar to make current remix output file work will ShellGrid format
        !! TODO: One day we will remove this once we implement a more standardized ShellGrid-friendly output file format
        type(IOVAR_T), dimension(MAXIOVARS), intent(in) :: IOV
        character(len=strLen), intent(in) :: vName
        type(ShellGrid_T), intent(in) :: shGr
        type(ShellGridVar_T), intent(inout) :: v

        !real(rp), dimension(shGr%Np, shGr%Nt+1) :: tmpVar
        real(rp), dimension(shGr%js:shGr%je, shGr%is:shGr%ie) :: tmpVar
        integer :: idx

        associate (isg=>shGr%isg, ieg=>shGr%ieg, jsg=>shGr%jsg, jeg=>shGr%jeg,\
            is =>shGr%is , ie =>shGr%ie , js =>shGr%js , je =>shGr%je )

        call IOArray2DFill(IOV,vName,tmpVar)
        v%data(is+1:ie,js:je) = transpose(tmpVar)

        ! fix pole value
        ! NOTE, THIS IS ASSUMING IS=1 (POLE)
        ! NOTE, STUPIDLY ASSUMING UNIFORM SPACING IN LONGITUDE
        ! BOTH OF THESE ASSUMPTIONS ARE HERE BECAUSE THIS ENTIRE CODE IS TEMPORARY
        ! AND IS ONLY USED FOR TESTING OF SHELLUTILS AND SHELLINTERP ON STANDARD REMIX FILES
        v%data(is,js:je) = sum(v%data(is+1,js:je))/size(v%data(is+1,js:je))

        ! fix periodic boundary (note, going through the ghosts in the first dimension)
        v%data(:,je+1)  = v%data(:,1)

        ! Copy last good cell through ghosts just so there's a value there
        ! note this doesn't get executed if isg=is (no ghosts)
        do idx=isg,is-1
            v%data(idx,:) = v%data(is,:)
        enddo

        ! note this doesn't get executed if ieg=ie (no ghosts)
        do idx=ie+1,ieg
            v%data(idx,:) = v%data(ie+1,:)
        enddo

        call wrapJ_SGV(shGr, v)
        
        ! But, we're gonna say that only our real domain is valid
        v%mask = .false.
        v%mask(is:ie+1, js:je+1) = .true.

        end associate

    end subroutine readVarJank


    subroutine outputRMSG(rmReader, fname, isFirst, gStrO)
        !! Write rmReader stuff to file
        !! Pretty much just for debugging
        type(rmReader_T), intent(in) :: rmReader
        character(len=*), intent(in) :: fname
        logical, intent(in) :: isFirst
        character(len=*), optional, intent(in) :: gStrO

        logical :: gExist

        type(IOVAR_T), dimension(10) :: IOVars

        if (isFirst) then

            call CheckAndKill(fname, .true.)
        endif
        
        if (.not. ioExist(fname, "sh_th")) then
            call ClearIO(IOVars)
            call AddOutVar(IOVars,"sh_th",rmReader%shGr%th)
            call AddOutVar(IOVars,"sh_ph",rmReader%shGr%ph)
            call WriteVars(IOVars,.true.,fname)
        endif


        ! If still here, we are good and need to write more stuff to file
        ! If still here, varO and vNameO present
        if (present(gStrO)) then
            call ClearIO(IOVars)
            call AddOutVar(IOVars, "time", rmReader%time)
            call AddOutVar(IOVars, "N1_time", rmReader%rmN1%time)
            call AddOutVar(IOVars, "nsPot NORTH data", rmReader%nsPot(NORTH)%data)
            call AddOutVar(IOVars, "nsPot NORTH mask", merge(1.0_rp, 0.0_rp, rmReader%nsPot(NORTH)%mask))
            call WriteVars(IOVars,.true.,fname,gStrO=gStrO)
        endif
    end subroutine outputRMSG

end module remixReader
