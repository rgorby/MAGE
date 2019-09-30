!Various routines to read/write HDF5 files

module gioH5
    use gamtypes
    use gamutils
    use gridutils
    use ioH5
    use multifluid
    
    implicit none

    integer, parameter :: MAXIOVAR = 50
    type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
    logical :: doRoot = .true. !Whether root variables need to be written

    !Necessary for IO routines
    character(len=strLen) ,public:: h5File
    logical :: fExist
    logical :: doWriteGhost  = .false.
    integer, parameter :: maxPlotVar = 25
    integer :: is,ie,js,je,ks,ke !Variable bounds for output
    integer :: GhostCells(3,1,1)

    contains

    subroutine readH5Grid(Model,Grid,inH5)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Grid
        character(len=*), intent(in) :: inH5

        integer :: Nd
        logical :: fExist
        integer, dimension(NDIM) :: dims        

        !Reset IO chain
        call ClearIO(IOVars)

        inquire(file=inH5,exist=fExist)
        if (.not. fExist) then
            !Error out and leave
            write(*,*) 'Unable to open input mesh, exiting'
            stop
        endif

        !Setup input chain
        call AddInVar(IOVars,"X")
        call AddInVar(IOVars,"Y")
        call AddInVar(IOVars,"Z")

        call ReadVars(IOVars,.false.,inH5) !Don't use io precision

        Nd = IOVars(1)%Nr !Dimension
        if (Nd <3) then
            write(*,*) "Number of dimensions not supported"
            stop
        endif
        dims = IOVars(1)%dims(1:Nd)

        !Start with main indices, convert # of ghost corners to active
        Grid%Nip = dims(1) - 2*Model%nG - 1
        Grid%Njp = dims(2) - 2*Model%nG - 1
        Grid%Nkp = dims(3) - 2*Model%nG - 1
    
        !Derived quantities
        Grid%Ni = Grid%Nip + 2*Model%nG
        Grid%Nj = Grid%Njp + 2*Model%nG
        Grid%Nk = Grid%Nkp + 2*Model%nG
        Grid%is = 1; Grid%ie = Grid%Nip
        Grid%js = 1; Grid%je = Grid%Njp
        Grid%ks = 1; Grid%ke = Grid%Nkp
    
        Grid%isg = Grid%is-Model%nG
        Grid%ieg = Grid%ie+Model%nG
        Grid%jsg = Grid%js-Model%nG
        Grid%jeg = Grid%je+Model%nG
        Grid%ksg = Grid%ks-Model%nG
        Grid%keg = Grid%ke+Model%nG
    
        allocate(Grid%x(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1))
        allocate(Grid%y(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1))
        allocate(Grid%z(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1))

        Grid%x = reshape(IOVars(XDIR)%data,[Grid%Ni+1,Grid%Nj+1,Grid%Nk+1])
        Grid%y = reshape(IOVars(YDIR)%data,[Grid%Ni+1,Grid%Nj+1,Grid%Nk+1])
        Grid%z = reshape(IOVars(ZDIR)%data,[Grid%Ni+1,Grid%Nj+1,Grid%Nk+1])

    
    end subroutine readH5Grid

    !Write initial grid info to root of H5 output file
    subroutine writeH5GridInit(Model,Gr)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr

        real (rp), dimension(:,:,:),   allocatable :: gQ !Grid quality
        logical :: isExist

        character(len=strLen) :: vID

        !Don't call this function again
        doRoot = .false.

        !Test if root variables (grid/force info is already there)
        vID = "X" !Value to test for
        isExist = ioExist(h5File,vID)

        if (isExist) then
            return
        endif
        
        !Calculate grid quality
        call allocGridVar(Model,Gr,gQ)
        call GridQuality(Model,Gr,gQ) 

        !Reset IO chain
        call ClearIO(IOVars)

        !Get bounds (use global doWriteGhost)
        call getBds(Gr)

        !Fill IO chain, start with coordinates
        if (Gr%Nkp > 1) then
            !3D problem
            call AddOutVar(IOVars,"X",Gr%x(is:ie+1,js:je+1,ks:ke+1))            
            call AddOutVar(IOVars,"Y",Gr%y(is:ie+1,js:je+1,ks:ke+1))
            call AddOutVar(IOVars,"Z",Gr%z(is:ie+1,js:je+1,ks:ke+1))
        else
            !2D problem
            !Squash corner arrays to 2D
            call AddOutVar(IOVars,"X",reshape(Gr%x(is:ie+1,js:je+1,ks:ks),[ie-is+2,je-js+2]))
            call AddOutVar(IOVars,"Y",reshape(Gr%y(is:ie+1,js:je+1,ks:ks),[ie-is+2,je-js+2]))
        endif

        call AddOutVar(IOVars,"dV",Gr%volume(is:ie,js:je,ks:ke))
        call AddOutVar(IOVars,"gQ",       gQ(is:ie,js:je,ks:ke))
        if (Model%doMHD .and. Model%doBackground) then
            !Write out background field and force density
            call AddOutVar(IOVars,"Bx0"  ,Gr%B0  (is:ie,js:je,ks:ke,XDIR))
            call AddOutVar(IOVars,"By0"  ,Gr%B0  (is:ie,js:je,ks:ke,YDIR))
            call AddOutVar(IOVars,"Bz0"  ,Gr%B0  (is:ie,js:je,ks:ke,ZDIR))
            call AddOutVar(IOVars,"dPxB0",Gr%dpB0(is:ie,js:je,ks:ke,XDIR))
            call AddOutVar(IOVars,"dPyB0",Gr%dpB0(is:ie,js:je,ks:ke,YDIR))
            call AddOutVar(IOVars,"dPzB0",Gr%dpB0(is:ie,js:je,ks:ke,ZDIR))

        endif
        if (Model%doGrav) then
            !Write out grav accelerations
            call AddOutVar(IOVars,"gx",Gr%gxyz(is:ie,js:je,ks:ke,XDIR))
            call AddOutVar(IOVars,"gy",Gr%gxyz(is:ie,js:je,ks:ke,YDIR))
            call AddOutVar(IOVars,"gz",Gr%gxyz(is:ie,js:je,ks:ke,ZDIR))

        endif

        !Add information about time scaling/units
        call AddOutVar(IOVars,"tScl"   ,Model%gamOut%tScl)
        call AddOutVar(IOVars,"timeID" ,Model%gamOut%tID)
        call AddOutVar(IOVars,"UnitsID",Model%gamOut%uID)
        !Write out the chain
        call WriteVars(IOVars,.true.,h5File)
        
    end subroutine writeH5GridInit

    !Write specific output slice
    subroutine writeSlc(Model,Gr,State,gStr)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        type(State_T), intent(in) :: State
        character(len=*), intent(in) :: gStr

        integer :: i,j,k,s
        character(len=strLen) :: dID,VxID,VyID,VzID,PID

        !Fill base data
        real(rp), dimension(:,:,:),   allocatable :: gVar,DivBcc
        real(rp), dimension(:,:,:,:), allocatable :: gVec
        real (rp), dimension(:,:,:,:), allocatable :: VecA,VecB !Full-sized arrays
        real(rp) :: totDivB

        !Check if root variables need to be written
        if (doRoot) then
            call writeH5GridInit(Model,Gr)
            doRoot = .false.
        endif
        !Reset IO chain
        call ClearIO(IOVars)

        !Get bounds (use global doWriteGhost)
        call getBds(Gr)

        !Allocate holders
        allocate(gVar(is:ie,js:je,ks:ke))
        allocate(gVec(is:ie,js:je,ks:ke,1:NDIM))

        associate(gamOut=>Model%gamOut)

        do s=0,Model%nSpc
            if (s == 0) then
                dID = "D"
                VxID = "Vx"
                VyID = "Vy"
                VzID = "Vz"
                PID = "P"
            else
                write(dID ,'(A,I0)') "D" , s
                write(VxID,'(A,I0)') "Vx", s
                write(VyID,'(A,I0)') "Vy", s
                write(VzID,'(A,I0)') "Vz", s
                write(PID ,'(A,I0)') "P" , s
            endif

            !Density
            
            call GameraOut(dID,gamOut%dID,gamOut%dScl,State%Gas(is:ie,js:je,ks:ke,DEN,s))

            !---------------------
            !Calculate Velocities/Pressure
            !$OMP PARALLEL DO default(shared) collapse(2)
            do k=ks,ke
                do j=js,je
                    do i=is,ie
                        if (State%Gas(i,j,k,DEN,s)>TINY) then
                            gVec(i,j,k,:) = State%Gas(i,j,k,MOMX:MOMZ,s)/State%Gas(i,j,k,DEN,s)
                            gVar(i,j,k) = (Model%gamma-1)*( State%Gas(i,j,k,ENERGY,s) - &
                                          0.5*(State%Gas(i,j,k,MOMX,s)**2 + &
                                               State%Gas(i,j,k,MOMY,s)**2 + &
                                               State%Gas(i,j,k,MOMZ,s)**2)/ &
                                               State%Gas(i,j,k,DEN ,s))
                        else
                            gVec(i,j,k,:) = 0.0
                            gVar(i,j,k)   = 0.0
                        endif
                    enddo
                enddo
            enddo


            !Add V/P to chain
            call GameraOut(VxID,gamOut%vID,gamOut%vScl,gVec(is:ie,js:je,ks:ke,XDIR))
            call GameraOut(VyID,gamOut%vID,gamOut%vScl,gVec(is:ie,js:je,ks:ke,YDIR))
            call GameraOut(VzID,gamOut%vID,gamOut%vScl,gVec(is:ie,js:je,ks:ke,ZDIR))
            call GameraOut(PID ,gamOut%pID,gamOut%pScl,gVar(is:ie,js:je,ks:ke))

        enddo !Species loop
        !---------------------
        !Write MHD variables
        if (Model%doMHD) then
            call allocGridVec(Model,Gr,VecA)
            call allocGridVec(Model,Gr,VecB)

            !For current, use VecA to hold total Bxyz, VecB for Jxyz
            if (Model%doBackground) then
                gVec(:,:,:,:) = Gr%B0(is:ie,js:je,ks:ke,XDIR:ZDIR) + State%Bxyz(is:ie,js:je,ks:ke,XDIR:ZDIR)
                VecA = State%Bxyz + Gr%B0
            else
                gVec(:,:,:,:) = State%Bxyz(is:ie,js:je,ks:ke,XDIR:ZDIR)
                VecA = State%Bxyz
            endif

            !Add magnetic fields
            call GameraOut("Bx",gamOut%bID,gamOut%bScl,gVec(is:ie,js:je,ks:ke,XDIR))
            call GameraOut("By",gamOut%bID,gamOut%bScl,gVec(is:ie,js:je,ks:ke,YDIR))
            call GameraOut("Bz",gamOut%bID,gamOut%bScl,gVec(is:ie,js:je,ks:ke,ZDIR))

            !Write current
            call bFld2Jxyz(Model,Gr,VecA,VecB)
            gVec(:,:,:,:) = VecB(is:ie,js:je,ks:ke,XDIR:ZDIR)
            call FixRAVec(gVec)

            call AddOutVar(IOVars,"Jx",gVec(:,:,:,XDIR))
            call AddOutVar(IOVars,"Jy",gVec(:,:,:,YDIR))
            call AddOutVar(IOVars,"Jz",gVec(:,:,:,ZDIR))

            !Calculate/Write xyz electric fields
            !Divide by edge-length to go from potential to field
            !$OMP PARALLEL DO default(shared) collapse(2)
            do k=Gr%ksg,Gr%keg
                do j=Gr%jsg,Gr%jeg
                    do i=Gr%isg,Gr%ieg
                        if(any(Gr%Edge(i,j,k,:) == 0)) then
                             VecA(i,j,k,:) = 0.0
                        else
                             VecA(i,j,k,:) = State%Efld(i,j,k,:)/Gr%Edge(i,j,k,:)
                        end if
                    enddo
                enddo
            enddo
            call Eijk2xyz(Model,Gr,VecA,VecB)
            gVec(:,:,:,:) = VecB(is:ie,js:je,ks:ke,XDIR:ZDIR)
            call FixRAVec(gVec)
            call AddOutVar(IOVars,"Ex",gVec(:,:,:,XDIR))
            call AddOutVar(IOVars,"Ey",gVec(:,:,:,YDIR))
            call AddOutVar(IOVars,"Ez",gVec(:,:,:,ZDIR))
            

            if(Model%doResistive) then
                gVec(:,:,:,:) = State%Deta(is:ie,js:je,ks:ke,XDIR:ZDIR)
                call AddOutVar(IOVars,"Etax",gVec(:,:,:,XDIR))
                call AddOutVar(IOVars,"Etay",gVec(:,:,:,YDIR))
                call AddOutVar(IOVars,"Etaz",gVec(:,:,:,ZDIR))
            end if

            !Write divergence if necessary
            if (Model%doDivB) then
                call allocGridVar(Model,Gr,DivBcc)
                call DivB(Model,Gr,State,totDivB,DivBcc)
                gVar = DivBcc(is:ie,js:je,ks:ke)
                call AddOutVar(IOVars,"DivB",gVar)
                deallocate(DivBcc)
            endif
            deallocate(VecA,VecB)
        endif
        deallocate(gVec,gVar)

        !---------------------
        !Do attributes
        call AddOutVar(IOVars,"time",gamOut%tScl*Model%t)
        call AddOutVar(IOVars,"timestep",Model%ts)
        call AddOutVar(IOVars,"dt",gamOut%tScl*Model%dt)

        !---------------------
        !Call user routine
        !FIXME: Add this

        !------------------
        !Finalize
        call WriteVars(IOVars,.true.,h5File,trim(gStr))

        end associate

        contains
            subroutine GameraOut(vID,uID,vScl,V)
                character(len=*), intent(in) :: vID,uID
                real(rp), intent(in) :: vScl
                real(rp), intent(in) :: V(is:ie,js:je,ks:ke)

                integer :: n0
                call AddOutVar(IOVars,vID,V)
                n0 = FindIO(IOVars,vID)
                IOVars(n0)%scale = vScl
                IOVars(n0)%unitStr = uID
            end subroutine GameraOut

            !Fix up cell-centered vector (like current or electric field) to deal with axis
            !Note, this is only for output purposes since we don't have proper ghost information

            !Estimate Qxyz at pole by averaging about second ring
            !Calculate around first ring to interpolate pole/ring-2 values
            subroutine FixRAVec(Qxyz)
                real(rp), intent(inout) :: Qxyz(is:ie,js:je,ks:ke,1:NDIM)

                real(rp) :: Q0(NDIM)
                integer :: i,k
                real(rp) :: w1,w2
                !Check if we're doing ring avg and have either pole
                if ( Model%doRing .and. (Model%Ring%doS .or. Model%Ring%doE) ) then

                    select case (Model%Ring%GridID)
                        !------------------
                        case ("lfm")
                            !Move along axis
                            do i=is,ie
                                if (Model%Ring%doS) then
                                    call FixRAVec_S(Model,Gr,Qxyz(i,js:js+1,ks:ke,1:NDIM))
                                endif
                                if (Model%Ring%doE) then
                                    call FixRAVec_E(Model,Gr,Qxyz(i,je-1:je,ks:ke,1:NDIM))
                                endif
                            enddo !I-loop
                    end select

                endif
            end subroutine FixRAVec
    end subroutine writeSlc

    subroutine GridQuality(Model,Gr,gQ)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        real(rp), intent(inout) :: gQ(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg)

        integer :: i,j,k
        real(rp), dimension(NDIM) :: dXcc, ddV,ddx,gdV

        gQ = 0.0
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,dXcc,ddV,ddx,gdV)
        do k=Gr%ks, Gr%ke
            do j=Gr%js, Gr%je
                do i=Gr%is, Gr%ie
                    dXcc = [Gr%di(i,j,k),Gr%dj(i,j,k),Gr%dk(i,j,k)]
                    
                    ddV = [Gr%volume(i+1,j,k)-Gr%volume(i-1,j,k), &
                           Gr%volume(i,j+1,k)-Gr%volume(i,j-1,k), &
                           Gr%volume(i,j,k+1)-Gr%volume(i,j,k-1) ]
                    ddx = [norm2( Gr%xyzcc(i+1,j,k,:)-Gr%xyzcc(i-1,j,k,:) ), &
                           norm2( Gr%xyzcc(i,j+1,k,:)-Gr%xyzcc(i,j-1,k,:) ), &
                           norm2( Gr%xyzcc(i,j,k+1,:)-Gr%xyzcc(i,j,k-1,:) ) ]
                    gdV(IDIR) = ddV(IDIR)/ddx(IDIR) !Gradient dV
                    gdV(JDIR) = ddV(JDIR)/ddx(JDIR)
                    gdV(KDIR) = ddV(KDIR)/ddx(KDIR)

                    gQ(i,j,k) = norm2(gdV)*norm2(dXcc)/Gr%volume(i,j,k)
                enddo
            enddo
        enddo
    end subroutine GridQuality

    !Write restart dump to "ResF" output file
    subroutine writeH5Res(Model,Gr,State,ResF)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        type(State_T), intent(in) :: State
        character(len=*), intent(in) :: ResF

        !Reset IO chain
        call ClearIO(IOVars)

        !Main attributes
        call AddOutVar(IOVars,"nOut",Model%nOut)
        call AddOutVar(IOVars,"nRes",Model%nRes)
        call AddOutVar(IOVars,"ts"  ,Model%ts)
        call AddOutVar(IOVars,"t"   ,Model%t)

        !Coordinates of corners
        call AddOutVar(IOVars,"X",Gr%x)
        call AddOutVar(IOVars,"Y",Gr%y)
        call AddOutVar(IOVars,"Z",Gr%z)

        !Set bounds for active cell centers
        call getBds(Gr,.false.)

        !State variable
        call AddOutVar(IOVars,"Gas",State%Gas(is:ie,js:je,ks:ke,:,:))
        if (Model%doMHD) then
            call AddOutVar(IOVars,"magFlux",State%magFlux(is:ie+1,js:je+1,ks:ke+1,:))
        endif

        !Write out, force real precision
        call WriteVars(IOVars,.false.,ResF)

    end subroutine writeH5Res
    
    subroutine readH5Restart(Model,Gr,State,inH5,doResetO,tResetO)
        type(Model_T), intent(inout) :: Model
        type(Grid_T),  intent(inout) :: Gr
        type(State_T), intent(inout) :: State
        character(len=*), intent(in) :: inH5
        logical, intent(in), optional :: doResetO
        real(rp), intent(in), optional :: tResetO

        logical :: doReset
        real(rp) :: tReset
        integer :: wDims(5),bDims(4)
        integer :: rSpc

        !Test for resetting
        if (present(doResetO)) then
            doReset = doResetO
            if (present(tResetO)) then
                tReset = tResetO
            else
                tReset = 0.0
            endif
        else
            doReset = .false.
            tReset = 0.0
        endif

        write(*,*) 'Reading restart from ', trim(inH5)
        inquire(file=inH5,exist=fExist)
        if (.not. fExist) then
            !Error out and leave
            write(*,*) 'Unable to open input restart file, exiting'
            stop
        endif

        !Reset IO chain
        call ClearIO(IOVars)

        call AddInVar(IOVars,"Gas")    
        call AddInVar(IOVars,"magFlux")
        call AddInVar(IOVars,"nOut",vTypeO=IOINT)
        call AddInVar(IOVars,"nRes",vTypeO=IOINT)
        call AddInVar(IOVars,"ts"  ,vTypeO=IOINT)
        call AddInVar(IOVars,"t"   ,vTypeO=IOREAL)

        !Get data
        call ReadVars(IOVars,.false.,inH5)

        !Set sizes/bounds
        call getBds(Gr,.false.)

        !Find number of species in restart
        rSpc = IOVars(1)%dims(5)-1

        if (Model%nSpc == rSpc) then
            !Restart and State variable agree
            wDims = [Gr%Nip,Gr%Njp,Gr%Nkp,NVAR,Model%nSpc+1]
            State%Gas(is:ie,js:je,ks:ke,:,:) = reshape(IOVars(1)%data,wDims)
        else if (Model%nSpc > rSpc) then
            !Not enough species in restart, fill as many as possible
            wDims = [Gr%Nip,Gr%Njp,Gr%Nkp,NVAR,rSpc+1]
            State%Gas(is:ie,js:je,ks:ke,:,0:rSpc) = reshape(IOVars(1)%data,wDims)
            !Now initialize to empty remaining species
            State%Gas(:,:,:,DEN,rSpc+1:Model%nSpc) = dFloor
            State%Gas(:,:,:,MOMX:MOMZ,rSpc+1:Model%nSpc) = 0.0
            State%Gas(:,:,:,ENERGY,rSpc+1:Model%nSpc) = pFloor/(Model%gamma-1)
            !Now reaccumulate
            call State2Bulk(Model,Gr,State)
        else
            !Too many species in restart, this isn't good
            write(*,*) 'Restart error, more species in restart than room in State!'
            stop
        endif

        !Now handle magnetic fields
        bDims = [Gr%Nip+1,Gr%Njp+1,Gr%Nkp+1,3]
        !NOTE: For now, lazily assuming order
        !Should use FindIO routine
        State%magFlux(is:ie+1,js:je+1,ks:ke+1,:) = reshape(IOVars(2)%data,bDims)

        !Get main attributes
        if (doReset) then
            Model%nOut = 0
            Model%nRes = int(IOVars(4)%data(1)) + 1
            Model%ts = 0
            Model%t = tReset
        else
            Model%nOut = int(IOVars(3)%data(1))
            Model%nRes = int(IOVars(4)%data(1)) + 1
            Model%ts   = int(IOVars(5)%data(1))
            Model%t = IOVars(6)%data(1)
        endif        
    !Do touchup to data structures
        State%time = Model%t
        Model%tOut = floor(Model%t/Model%dtOut)*Model%dtOut
        Model%tRes = Model%t + Model%dtRes

    end subroutine readH5Restart

    !Output black box from crash
    subroutine WriteBlackBox(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        type(State_T), intent(in) :: State

        !Reset output file
        h5File = "CRASH" // trim(Model%RunID) // ".h5"
        doRoot = .true. !Make sure root vars are rewritten
        call CheckAndKill(h5File)

        call writeSlc(Model,Gr,State,"Step#0")

    end subroutine WriteBlackBox

    !Gets variable bounds, uses optional input logical or global doWriteGhost
    subroutine getBds(Grid,doGhostOpt)
        type(Grid_T), intent(in) :: Grid
        logical, intent(in), optional :: doGhostOpt

        logical :: doIncludeG

        if (present(doGhostOpt)) then
            doIncludeG = doGhostOpt
        else
            doIncludeG = doWriteGhost
        endif

        if (doIncludeG) then
            is = Grid%isg
            ie = Grid%ieg
            js = Grid%jsg
            je = Grid%jeg
            ks = Grid%ksg
            ke = Grid%keg
        else
            is = Grid%is
            ie = Grid%ie
            js = Grid%js
            je = Grid%je
            ks = Grid%ks
            ke = Grid%ke
        endif
        if (Grid%Nkp == 1) then
            ke = ks
        endif

    end subroutine getBds


end module gioH5
