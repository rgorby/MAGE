! Main data objects and functions to perform a gamera simulation

module gamapp_mpi
    use gamtypes
    use step
    use init
    use mhdgroup
    use gamapp
    use bcs_mpi
    use mpidefs
    use mpi

    implicit none

    type, extends(GamApp_T) :: gamAppMpi_T
        integer :: gamMpiComm = MPI_COMM_NULL
        integer, dimension(:), allocatable :: sendRanks, recvRanks
        integer, dimension(:), allocatable :: sendCountsGas, sendTypesGas
        integer, dimension(:), allocatable :: recvCountsGas, recvTypesGas
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: sendDisplsGas, recvDisplsGas
    end type gamAppMpi_T

    contains

    subroutine initGamera_mpi(gamAppMpi, userInitFunc, gamComm, optFilename, doIO)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi
        procedure(StateIC_T), pointer, intent(in) :: userInitFunc
        integer, intent(in) :: gamComm
        character(len=*), optional, intent(in) :: optFilename
        logical, optional, intent(in) :: doIO

        character(len=strLen) :: inpXML
        type(XML_Input_T) :: xmlInp
        logical :: doIOX

        if(present(optFilename)) then
            ! read from the prescribed file
            inpXML = optFilename
        else
            !Find input deck
            call getIDeckStr(inpXML)
        endif
        call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")

        if (present(doIO)) then
            doIOX = doIO
        else
            doIOX = .true.
        endif

        !Create XML reader
        write(*,*) 'Reading input deck from ', trim(inpXML)
        xmlInp = New_XML_Input(trim(inpXML),'Gamera',.true.)

        !Initialize Grid/State/Model (Hatch Gamera)
        !Will enforce 1st BCs, caculate 1st timestep, set oldState
        call Hatch_mpi(gamAppMpi,xmlInp,userInitFunc,gamComm)
        call cleanClocks()

        if (doIOX) then
            if (.not. gamAppMpi%Model%isRestart) call fOutput(gamAppMpi%Model,gamAppMpi%Grid,gamAppMpi%State)
            call consoleOutput(gamAppMpi%Model,gamAppMpi%Grid,gamAppMpi%State)
        endif

    end subroutine initGamera_mpi

    subroutine Hatch_mpi(gamAppMpi,xmlInp,userInitFunc,gamComm,endTime)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi
        type(XML_Input_T), intent(inout) :: xmlInp
        procedure(StateIC_T), pointer, intent(in) :: userInitFunc
        integer, intent(in) :: gamComm
        real(rp), optional, intent(in) :: endTime

        integer :: numNeighbors, ierr, length, commSize, rank, ic, jc, kc, targetRank, sendNumber, numInNeighbors, numOutNeighbors, listIndex, localIndexIn, localIndexOut, sendDataOffset, recvDataOffset
        integer, dimension(27) :: sourceRanks, sourceData
        logical :: reorder,periodicI,periodicJ,periodicK,wasWeighted
        character(len=strLen) :: message, bcType
        real(rp), dimension(:,:,:), allocatable :: tempX,tempY,tempZ
        real(rp) :: tmpDT
        integer :: dataSize, iG, iP, iGjG, iGjP, iPjG, iPjP, transDataType
        integer :: cornerMpiType,iEdgeMpiType,jEdgeMpiType,kEdgeMpiType,iFaceMpiType,jFaceMpiType,kFaceMpiType
        integer :: corner4MpiType,iEdge4MpiType,jEdge4MpiType,kEdge4MpiType,iFace4MpiType,jFace4MpiType,kFace4MpiType
        integer :: corner5MpiType,iEdge5MpiType,jEdge5MpiType,kEdge5MpiType,iFace5MpiType,jFace5MpiType,kFace5MpiType
        integer(kind=MPI_ADDRESS_KIND) :: tempOffsets(2)

        associate(Grid=>gamAppMpi%Grid,Model=>gamAppMpi%Model)

        ! call appropriate subroutines to read corner info and mesh size data
        call ReadCorners(Model,Grid,xmlInp,endTime)

        ! setup the distributed graph topology communicator for future MPI work
        call xmlInp%Set_Val(Grid%NumRi,'iPdir/N',1)
        call xmlInp%Set_Val(Grid%NumRj,'jPdir/N',1)
        call xmlInp%Set_Val(Grid%NumRk,'kPdir/N',1)

        if(Grid%NumRi*Grid%NumRj*Grid%NumRk > 1) then
            Grid%isTiled = .true.
        endif

        if(Grid%isTiled) then
            call mpi_comm_size(gamComm, commSize, ierr)
            if(commSize /= Grid%NumRi*Grid%NumRj*Grid%NumRk) then
                print *,'Expected ',Grid%NumRi*Grid%NumRj*Grid%NumRk,' ranks, but this job has ',commSize,' gamera ranks'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif

            if(modulo(Grid%Nip,Grid%NumRi) /= 0) then
                print *,'Number of cells in i dimension is ',Grid%Nip,' but it must be a multiple of the number of MPI ranks in that dimension, which was ',Grid%NumRi
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif
            if(modulo(Grid%Njp,Grid%NumRj) /= 0) then
                print *,'Number of cells in j dimension is ',Grid%Njp,' but it must be a multiple of the number of MPI ranks in that dimension, which was ',Grid%NumRj
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif
            if(modulo(Grid%Nkp,Grid%NumRk) /= 0) then
                print *,'Number of cells in k dimension is ',Grid%Nkp,' but it must be a multiple of the number of MPI ranks in that dimension, which was ',Grid%NumRk
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif

            ! rank assignment. These are now 0 to N-1, matching the way that MPI counts it
            call mpi_comm_rank(gamComm, rank, ierr)
            Grid%Rk = modulo(rank, Grid%NumRk)
            rank = (rank-Grid%Rk)/Grid%NumRk
            Grid%Rj = modulo(rank, Grid%NumRj)
            rank = (rank-Grid%Rj)/Grid%NumRj
            Grid%Ri = rank

            ! adjust corner information to reflect this individual node's grid data
            Grid%Nip = Grid%Nip/Grid%NumRi
            Grid%Njp = Grid%Njp/Grid%NumRj
            Grid%Nkp = Grid%Nkp/Grid%NumRk

            Grid%Ni = Grid%Nip + 2*Model%nG
            Grid%Nj = Grid%Njp + 2*Model%nG
            Grid%Nk = Grid%Nkp + 2*Model%nG

            ! check which dimensions are using MPI for periodicity
            call xmlInp%Set_Val(bcType,'ibc/bc','')
            if (bcType == 'periodic') then
                periodicI = .true.
            else
                periodicI = .false.
            endif
            call xmlInp%Set_Val(bcType,'jbc/bc','')
            if (bcType == 'periodic') then
                periodicJ = .true.
            else
                periodicJ = .false.
            endif
            call xmlInp%Set_Val(bcType,'kbc/bc','')
            if (bcType == 'periodic') then
                periodicK = .true.
            else
                periodicK = .false.
            endif

            ! create MPI topology
            sourceRanks(:) = -1
            sourceData(:) = 0
            numNeighbors = 0
            call mpi_comm_rank(gamComm, rank, ierr)
            do ic=-1,1
                if((Grid%Ri+ic >= 0 .and. Grid%Ri+ic < Grid%NumRi) .or. periodicI) then
                    do jc=-1,1
                        if((Grid%Rj+jc >= 0 .and. Grid%Rj+jc < Grid%NumRj) .or. periodicJ) then
                            do kc=-1,1
                                if((Grid%Rk+kc >= 0 .and. Grid%Rk+kc < Grid%NumRk) .or. periodicK)then
                                    targetRank = modulo(Grid%Ri+ic,Grid%NumRi)*Grid%NumRk*Grid%NumRj + &
                                                 modulo(Grid%Rj+jc,Grid%NumRj)*Grid%NumRk + &
                                                 modulo(Grid%Rk+kc,Grid%NumRk)
                                    if(targetRank /= rank) then ! ensure I'm not talking to myself
                                        listIndex = findloc(sourceRanks, targetRank, 1)
                                        if(listIndex == 0) then ! this rank not in the list yet
                                            numNeighbors = numNeighbors+1
                                            listIndex = numNeighbors
                                            sourceRanks(listIndex) = targetRank
                                        endif
                                        ! calculate the size of this region that will be transmitted
                                        sendNumber = 1
                                        if (ic == 0) then
                                            sendNumber = sendNumber * Grid%Nip
                                        else
                                            sendNumber = sendNumber * Model%nG
                                        endif
                                        if (jc == 0) then
                                            sendNumber = sendNumber * Grid%Njp
                                        else
                                            sendNumber = sendNumber * Model%nG
                                        endif
                                        if (kc == 0) then
                                            sendNumber = sendNumber * Grid%Nkp
                                        else
                                            sendNumber = sendNumber * Model%nG
                                        endif
                                        ! add these cells to whatever we are already sending to that rank
                                        sourceData(listIndex) = sourceData(listIndex) + sendNumber
                                    endif
                                endif
                            enddo
                        endif
                    enddo
                endif
            enddo

            reorder = .true. ! allow MPI to reorder the ranks
            call mpi_dist_graph_create_adjacent(gamComm, &
                numNeighbors,sourceRanks,sourceData, &
                numNeighbors,sourceRanks,sourceData, & ! comms symmetrical
                MPI_INFO_NULL, reorder, gamAppMpi%gamMpiComm, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if

            ! get our final ranks from the new MPI topology
            call mpi_comm_rank(gamAppMpi%gamMpiComm, rank, ierr)
            Grid%Rk = modulo(rank, Grid%NumRk)
            rank = (rank-Grid%Rk)/Grid%NumRk
            Grid%Rj = modulo(rank, Grid%NumRj)
            rank = (rank-Grid%Rj)/Grid%NumRj
            Grid%Ri = rank

            ! whether this rank has external BCs
            Grid%hasLowerBC(IDIR) = Grid%Ri == 0
            Grid%hasLowerBC(JDIR) = Grid%Rj == 0
            Grid%hasLowerBC(KDIR) = Grid%Rk == 0
            Grid%hasUpperBC(IDIR) = Grid%Ri == (Grid%NumRi-1)
            Grid%hasUpperBC(JDIR) = Grid%Rj == (Grid%NumRj-1)
            Grid%hasUpperBC(KDIR) = Grid%Rk == (Grid%NumRk-1)

            ! adjust grid info for these ranks
            Grid%ijkShift(IDIR) = Grid%Nip*Grid%Ri
            Grid%ijkShift(JDIR) = Grid%Njp*Grid%Rj
            Grid%ijkShift(KDIR) = Grid%Nkp*Grid%Rk

            Grid%is = 1; Grid%ie = Grid%Nip
            Grid%js = 1; Grid%je = Grid%Njp
            Grid%ks = 1; Grid%ke = Grid%Nkp

            Grid%isg = Grid%is-Model%nG
            Grid%ieg = Grid%ie+Model%nG

            Grid%jsg = Grid%js-Model%nG
            Grid%jeg = Grid%je+Model%nG

            Grid%ksg = Grid%ks-Model%nG
            Grid%keg = Grid%ke+Model%nG

            ! create temporary arrays to hold this rank's subset of the full corner array
            allocate(tempX(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1))
            allocate(tempY(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1))
            allocate(tempZ(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1))

            ! pull out this rank's relevant corner info
            tempX = Grid%x(Grid%isg+Grid%ijkShift(IDIR):Grid%ieg+1+Grid%ijkShift(IDIR), &
                           Grid%jsg+Grid%ijkShift(JDIR):Grid%jeg+1+Grid%ijkShift(JDIR), &
                           Grid%ksg+Grid%ijkShift(KDIR):Grid%keg+1+Grid%ijkShift(KDIR))
            tempY = Grid%y(Grid%isg+Grid%ijkShift(IDIR):Grid%ieg+1+Grid%ijkShift(IDIR), &
                           Grid%jsg+Grid%ijkShift(JDIR):Grid%jeg+1+Grid%ijkShift(JDIR), &
                           Grid%ksg+Grid%ijkShift(KDIR):Grid%keg+1+Grid%ijkShift(KDIR))
            tempZ = Grid%z(Grid%isg+Grid%ijkShift(IDIR):Grid%ieg+1+Grid%ijkShift(IDIR), &
                           Grid%jsg+Grid%ijkShift(JDIR):Grid%jeg+1+Grid%ijkShift(JDIR), &
                           Grid%ksg+Grid%ijkShift(KDIR):Grid%keg+1+Grid%ijkShift(KDIR))

            ! delete the old corner arrays
            deallocate(Grid%x)
            deallocate(Grid%y)
            deallocate(Grid%z)

            ! move the new arrays to the grid corner arrays
            call move_alloc(tempX, Grid%x)
            call move_alloc(tempY, Grid%y)
            call move_alloc(tempZ, Grid%z)

            ! now create the arrays that MPI will use to send and receive the data
            call mpi_dist_graph_neighbors_count(gamAppMpi%gamMpiComm,numInNeighbors,numOutNeighbors,wasWeighted,ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            if (numInNeighbors /= numOutNeighbors) then
                print *,'Number of in edges and out edges did not match for rank ', rank
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif

            allocate(gamAppMpi%sendRanks(numOutNeighbors))
            allocate(gamAppMpi%recvRanks(numInNeighbors))
            ! don't care about the weights, dump them into an existing array
            call mpi_dist_graph_neighbors(gamAppMpi%gamMpiComm, numInNeighbors, gamAppMpi%recvRanks, sourceData, numOutNeighbors, gamAppMpi%sendRanks, sourceData, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if

            allocate(gamAppMpi%sendCountsGas(numOutNeighbors))
            allocate(gamAppMpi%sendDisplsGas(numOutNeighbors))
            allocate(gamAppMpi%sendTypesGas(numOutNeighbors))
            allocate(gamAppMpi%recvCountsGas(numInNeighbors))
            allocate(gamAppMpi%recvDisplsGas(numInNeighbors))
            allocate(gamAppMpi%recvTypesgas(numInNeighbors))

            ! assemble the different datatypes
            call mpi_type_extent(MPI_MYFLOAT, dataSize, ierr) ! number of bytes per array entry

            ! I dimension
            call mpi_type_contiguous(Model%nG, MPI_MYFLOAT, iG, ierr) ! ghosts i
            call mpi_type_contiguous(Grid%Nip, MPI_MYFLOAT, iP, ierr) ! physical i

            ! J dimension
            call mpi_type_hvector(Model%nG, 1, Grid%Ni*dataSize, iG, iGjG, ierr) ! ghosts i   - ghosts j
            call mpi_type_hvector(Grid%Njp, 1, Grid%Ni*dataSize, iG, iGjP, ierr) ! ghosts i   - physical j
            call mpi_type_hvector(Model%nG, 1, Grid%Ni*dataSize, iP, iPjG, ierr) ! physical i - ghosts j
            call mpi_type_hvector(Grid%Njp, 1, Grid%Ni*dataSize, iP, iPjP, ierr) ! physical i - physical j

            ! K dimension
            call mpi_type_hvector(Model%nG, 1, Grid%Ni*Grid%Nj*dataSize, iGjG, cornerMpiType, ierr) ! ghosts i   - ghosts j   - ghosts k
            call mpi_type_hvector(Model%nG, 1, Grid%Ni*Grid%Nj*dataSize, iGjP, jEdgeMpiType,  ierr) ! ghosts i   - physical j - ghosts k
            call mpi_type_hvector(Model%nG, 1, Grid%Ni*Grid%Nj*dataSize, iPjG, iEdgeMpiType,  ierr) ! physical i - ghosts j   - ghosts k
            call mpi_type_hvector(Grid%Nkp, 1, Grid%Ni*Grid%Nj*dataSize, iGjG, kEdgeMpiType,  ierr) ! ghosts i   - ghosts j   - physical k
            call mpi_type_hvector(Model%nG, 1, Grid%Ni*Grid%Nj*dataSize, iPjP, kFaceMpiType,  ierr) ! physical i - physical j - ghosts k
            call mpi_type_hvector(Grid%Nkp, 1, Grid%Ni*Grid%Nj*dataSize, iGjP, iFaceMpiType,  ierr) ! ghosts i   - physical j - physical k
            call mpi_type_hvector(Grid%Nkp, 1, Grid%Ni*Grid%Nj*dataSize, iPjG, jFaceMpiType,  ierr) ! physical i - ghosts j   - physical k

            ! 4th dimension
            call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, cornerMpiType, corner4MpiType, ierr)
            call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iEdgeMpiType,  iEdge4MpiType,  ierr)
            call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, jEdgeMpiType,  jEdge4MpiType,  ierr)
            call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, kEdgeMpiType,  kEdge4MpiType,  ierr)
            call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iFaceMpiType,  iFace4MpiType,  ierr)
            call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, jFaceMpiType,  jFace4MpiType,  ierr)
            call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, kFaceMpiType,  kFace4MpiType,  ierr)

            ! 5th dimension
            call mpi_type_hvector(Model%nSpc+1, 1, NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, corner4MpiType, corner5MpiType, ierr)
            call mpi_type_hvector(Model%nSpc+1, 1, NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iEdge4MpiType,  iEdge5MpiType,  ierr)
            call mpi_type_hvector(Model%nSpc+1, 1, NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, jEdge4MpiType,  jEdge5MpiType,  ierr)
            call mpi_type_hvector(Model%nSpc+1, 1, NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, kEdge4MpiType,  kEdge5MpiType,  ierr)
            call mpi_type_hvector(Model%nSpc+1, 1, NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iFace4MpiType,  iFace5MpiType,  ierr)
            call mpi_type_hvector(Model%nSpc+1, 1, NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, jFace4MpiType,  jFace5MpiType,  ierr)
            call mpi_type_hvector(Model%nSpc+1, 1, NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, kFace4MpiType,  kFace5MpiType,  ierr)

            ! counts are always 1 because we're sending a single (complicated) mpi datatype
            gamAppMpi%sendCountsGas(:) = 1
            gamAppMpi%recvCountsGas(:) = 1

            ! displacements are always 0 because the displacements are baked into each mpi datatype
            gamAppMpi%sendDisplsGas(:) = 0
            gamAppMpi%recvDisplsGas(:) = 0

            ! set all datatypes to null by default
            gamAppMpi%sendTypesGas(:) = MPI_DATATYPE_NULL
            gamAppMpi%recvTypesGas(:) = MPI_DATATYPE_NULL

            ! figure out exactly what data needs to be sent to (and received from) each neighbor
            ! create custom MPI datatypes to perform these transfers
            call mpi_comm_rank(gamAppMpi%gamMpiComm, rank, ierr)

            ! calculate receive types
            do ic=-1,1
                if((Grid%Ri+ic >= 0 .and. Grid%Ri+ic < Grid%NumRi) .or. periodicI) then
                    do jc=-1,1
                        if((Grid%Rj+jc >= 0 .and. Grid%Rj+jc < Grid%NumRj) .or. periodicJ) then
                            do kc=-1,1
                                if((Grid%Rk+kc >= 0 .and. Grid%Rk+kc < Grid%NumRk) .or. periodicK)then
                                    targetRank = modulo(Grid%Ri+ic,Grid%NumRi)*Grid%NumRk*Grid%NumRj + &
                                                 modulo(Grid%Rj+jc,Grid%NumRj)*Grid%NumRk + &
                                                 modulo(Grid%Rk+kc,Grid%NumRk)
                                    if(targetRank /= rank) then ! ensure I'm not talking to myself
                                        localIndexIn = findloc(gamAppMpi%recvRanks, targetRank, 1)

                                        recvDataOffset = 0
                                        SELECT CASE (ic)
                                            case (-1)
                                                ! min i side
                                                ! no change to recvDataOffset
                                            case (0)
                                                ! central in i dimension
                                                recvDataOffset = recvDataOffset + Model%nG
                                            case (1)
                                                ! max i side
                                                recvDataOffset = recvDataOffset + Model%nG + Grid%Nip
                                            CASE DEFAULT
                                                print *, 'Unexpected ic direction'
                                                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                                        END SELECT

                                        SELECT CASE (jc)
                                            case (-1)
                                                ! min j side
                                                ! no change to recvDataOffset
                                            case (0)
                                                ! central in j dimension
                                                recvDataOffset = recvDataOffset + Model%nG*Grid%Ni
                                            case (1)
                                                ! max j side
                                                recvDataOffset = recvDataOffset + (Model%nG + Grid%Njp)*Grid%Ni
                                            CASE DEFAULT
                                                print *, 'Unexpected jc direction'
                                                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                                        END SELECT

                                        SELECT CASE (kc)
                                            case (-1)
                                                ! min k side
                                                ! no change to recvDataOffset
                                            case (0)
                                                ! central in k dimension
                                                recvDataOffset = recvDataOffset + Model%nG*Grid%Ni*Grid%Nj
                                            case (1)
                                                ! max k side
                                                recvDataOffset = recvDataOffset + (Model%nG+Grid%Nkp)*Grid%Ni*Grid%Nj
                                            CASE DEFAULT
                                                print *, 'Unexpected kc direction'
                                                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                                        END SELECT

                                        ! determine which of the previously created datatypes is the correct one
                                        transDataType = MPI_DATATYPE_NULL
                                        SELECT CASE (abs(ic)+abs(jc)+abs(kc))
                                            case (1) ! face
                                                if(ic /= 0) then
                                                    transDataType = iFace5MpiType
                                                elseif(jc /= 0) then
                                                    transDataType = jFace5MpiType
                                                else
                                                    transDataType = kFace5MpiType
                                                endif
                                            case (2) ! edge
                                                if(ic == 0) then
                                                    transDataType = iEdge5MpiType
                                                elseif(jc == 0) then
                                                    transDataType = jEdge5MpiType
                                                else
                                                    transDataType = kEdge5MpiType
                                                endif
                                            case (3) ! corner
                                                transDataType = corner5MpiType
                                            CASE DEFAULT
                                                print *, 'Sum of ic+jc+kc is nonsense'
                                                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                                        END SELECT

                                        if(transDataType /= MPI_DATATYPE_NULL) then
                                            if(gamAppMpi%recvTypesGas(localIndexIn) == MPI_DATATYPE_NULL) then
                                                ! not receiving any data from this rank yet, just add this datatype
                                                call mpi_type_hindexed(1, (/ 1 /), recvDataOffset*dataSize, transDataType, gamAppMpi%recvTypesGas(localIndexIn), ierr)
                                            else
                                                ! we're already receivng other data from this rank
                                                !  merge the datatypes into a struct
                                                ! need to use a temporary array so that the ints are of type MPI_ADDRESS_KIND
                                                tempOffsets = (/ 0, recvDataOffset*dataSize /)
                                                call mpi_type_create_struct(2, (/ 1, 1 /), tempOffsets, (/ gamAppMpi%recvTypesGas(localIndexIn), transDataType /),  gamAppMpi%recvTypesGas(localIndexIn),  ierr)
                                            endif
                                        endif
                                    endif
                                endif
                            enddo
                        endif
                    enddo
                endif
            enddo

            ! now calculate the send types
            ! split into separate loops because the order of iterating must be inverted
            do ic=1,-1,-1
                if((Grid%Ri+ic >= 0 .and. Grid%Ri+ic < Grid%NumRi) .or. periodicI) then
                    do jc=1,-1,-1
                        if((Grid%Rj+jc >= 0 .and. Grid%Rj+jc < Grid%NumRj) .or. periodicJ) then
                            do kc=1,-1,-1
                                if((Grid%Rk+kc >= 0 .and. Grid%Rk+kc < Grid%NumRk) .or. periodicK)then
                                    targetRank = modulo(Grid%Ri+ic,Grid%NumRi)*Grid%NumRk*Grid%NumRj + &
                                                 modulo(Grid%Rj+jc,Grid%NumRj)*Grid%NumRk + &
                                                 modulo(Grid%Rk+kc,Grid%NumRk)
                                    if(targetRank /= rank) then ! ensure I'm not talking to myself
                                        localIndexOut = findloc(gamAppMpi%sendRanks, targetRank, 1)

                                        sendDataOffset = 0
                                        SELECT CASE (ic)
                                            case (-1)
                                                ! min i side
                                                sendDataOffset = sendDataOffset + Model%nG
                                            case (0)
                                                ! central in i dimension
                                                sendDataOffset = sendDataOffset + Model%nG
                                            case (1)
                                                ! max i side
                                                sendDataOffset = sendDataOffset + Grid%Nip ! ghosts+physical-ghosts
                                            CASE DEFAULT
                                                print *, 'Unexpected ic direction'
                                                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                                        END SELECT

                                        SELECT CASE (jc)
                                            case (-1)
                                                ! min j side
                                                sendDataOffset = sendDataOffset + Model%nG*Grid%Ni
                                            case (0)
                                                ! central in j dimension
                                                sendDataOffset = sendDataOffset + Model%nG*Grid%Ni
                                            case (1)
                                                ! max j side
                                                sendDataOffset = sendDataOffset + Grid%Njp*Grid%Ni
                                            CASE DEFAULT
                                                print *, 'Unexpected jc direction'
                                                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                                        END SELECT

                                        SELECT CASE (kc)
                                            case (-1)
                                                ! min k side
                                                sendDataOffset = sendDataOffset + Model%nG*Grid%Ni*Grid%Nj
                                            case (0)
                                                ! central in k dimension
                                                sendDataOffset = sendDataOffset + Model%nG*Grid%Ni*Grid%Nj
                                            case (1)
                                                ! max k side
                                                sendDataOffset = sendDataOffset + Grid%Nkp*Grid%Ni*Grid%Nj
                                            CASE DEFAULT
                                                print *, 'Unexpected kc direction'
                                                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                                        END SELECT

                                        ! determine which of the previously created datatypes is the correct one
                                        transDataType = MPI_DATATYPE_NULL
                                        SELECT CASE (abs(ic)+abs(jc)+abs(kc))
                                            case (1) ! face
                                                if(ic /= 0) then
                                                    transDataType = iFace5MpiType
                                                elseif(jc /= 0) then
                                                    transDataType = jFace5MpiType
                                                else
                                                    transDataType = kFace5MpiType
                                                endif
                                            case (2) ! edge
                                                if(ic == 0) then
                                                    transDataType = iEdge5MpiType
                                                elseif(jc == 0) then
                                                    transDataType = jEdge5MpiType
                                                else
                                                    transDataType = kEdge5MpiType
                                                endif
                                            case (3) ! corner
                                                transDataType = corner5MpiType
                                            CASE DEFAULT
                                                print *, 'Sum of ic+jc+kc is nonsense'
                                                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                                        END SELECT

                                        if(transDataType /= MPI_DATATYPE_NULL) then
                                            if(gamAppMpi%sendTypesGas(localIndexOut) == MPI_DATATYPE_NULL) then
                                                ! not sending any data to this rank yet, just add this datatype
                                                call mpi_type_hindexed(1, (/ 1 /), sendDataOffset*dataSize, transDataType, gamAppMpi%sendTypesGas(localIndexOut), ierr)
                                            else
                                                ! we're already sending other data to this rank
                                                !  merge the datatypes into a struct
                                                ! need to use a temporary array so that the ints are of type MPI_ADDRESS_KIND
                                                tempOffsets = (/ 0, sendDataOffset*dataSize /)
                                                call mpi_type_create_struct(2, (/ 1, 1 /), tempOffsets, (/ gamAppMpi%sendTypesGas(localIndexOut), transDataType /), gamAppMpi%sendTypesGas(localIndexOut), ierr)
                                            endif
                                        endif
                                    endif
                                endif
                            enddo
                        endif
                    enddo
                endif
            enddo

            ! commit the created MPI datatypes
            do localIndexOut=1,numOutNeighbors
                call mpi_type_commit(gamAppMpi%sendTypesGas(localIndexOut), ierr)
            enddo
            do localIndexIn=1,numInNeighbors
                call mpi_type_commit(gamAppMpi%recvTypesGas(localIndexIn), ierr)
            enddo

        endif

        ! call appropriate subroutines to calculate all appropriate grid data from the corner data
        call CalcGridInfo(Model,Grid,gamAppMpi%State,gamAppMpi%oState,gamAppMpi%Solver,xmlInp,userInitFunc)

        if(Grid%isTiled) then
            ! correct boundary conditions if necessary
            if(Grid%NumRi > 1) then
                ! MPI decomposed in I dimension

                ! check min I BC
                SELECT type(iiBC=>Grid%externalBCs(INI)%p)
                    TYPE IS (mpiNullBc_T)
                        ! this BC is already an MPI BC
                        if(Grid%hasLowerBC(IDIR) .and. .not. periodicI) then
                            print *, 'Min I BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                            stop
                        endif
                    CLASS DEFAULT
                        ! this BC is something other than an MPI BC
                        if(.not. Grid%hasLowerBC(IDIR) .or. periodicI) then
                            print *, 'Over-writing min I BC to be an MPI BC'
                            deallocate(Grid%externalBCs(INI)%p)
                            allocate(mpiNullBc_T :: Grid%externalBCs(INI)%p)
                        endif
                END SELECT

                ! check max I BC
                SELECT type(iiBC=>Grid%externalBCs(OUTI)%p)
                    TYPE IS (mpiNullBc_T)
                        ! this BC is already an MPI BC
                        if(Grid%hasUpperBC(IDIR) .and. .not. periodicI) then
                            print *, 'Max I BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                            stop
                        endif
                    CLASS DEFAULT
                        ! this BC is something other than an MPI BC
                        if(.not. Grid%hasUpperBC(IDIR) .or. periodicI) then
                            print *, 'Over-writing max I BC to be an MPI BC'
                            deallocate(Grid%externalBCs(OUTI)%p)
                            allocate(mpiNullBc_T :: Grid%externalBCs(OUTI)%p)
                        endif
                END SELECT

            endif

            if(Grid%NumRj > 1) then
                ! MPI decomposed in J dimension

                ! check min J BC
                SELECT type(iiBC=>Grid%externalBCs(INJ)%p)
                    TYPE IS (mpiNullBc_T)
                        ! this BC is already an MPI BC
                        if(Grid%hasLowerBC(JDIR) .and. .not. periodicJ) then
                            print *, 'Min J BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                            stop
                        endif
                    CLASS DEFAULT
                        ! this BC is something other than an MPI BC
                        if(.not. Grid%hasLowerBC(JDIR) .or. periodicJ) then
                            print *, 'Over-writing min J BC to be an MPI BC'
                            deallocate(Grid%externalBCs(INJ)%p)
                            allocate(mpiNullBc_T :: Grid%externalBCs(INJ)%p)
                        endif
                END SELECT

                ! check max J BC
                SELECT type(iiBC=>Grid%externalBCs(OUTJ)%p)
                    TYPE IS (mpiNullBc_T)
                        ! this BC is already an MPI BC
                        if(Grid%hasUpperBC(JDIR) .and. .not. periodicJ) then
                            print *, 'Max J BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                            stop
                        endif
                    CLASS DEFAULT
                        ! this BC is something other than an MPI BC
                        if(.not. Grid%hasUpperBC(JDIR) .or. periodicJ) then
                            print *, 'Over-writing max J BC to be an MPI BC'
                            deallocate(Grid%externalBCs(OUTJ)%p)
                            allocate(mpiNullBc_T :: Grid%externalBCs(OUTJ)%p)
                        endif
                END SELECT

            endif

            if(Grid%NumRk > 1) then
                ! MPI decomposed in K dimension

                ! check min K BC
                SELECT type(iiBC=>Grid%externalBCs(INK)%p)
                    TYPE IS (mpiNullBc_T)
                        ! this BC is already an MPI BC
                        if(Grid%hasLowerBC(KDIR) .and. .not. periodicK) then
                            print *, 'Min K BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                            stop
                        endif
                    CLASS DEFAULT
                        ! this BC is something other than an MPI BC
                        if(.not. Grid%hasLowerBC(KDIR) .or. periodicK) then
                            print *, 'Over-writing min K BC to be an MPI BC'
                            deallocate(Grid%externalBCs(INK)%p)
                            allocate(mpiNullBc_T :: Grid%externalBCs(INK)%p)
                        endif
                END SELECT

                ! check max K BC
                SELECT type(iiBC=>Grid%externalBCs(OUTK)%p)
                    TYPE IS (mpiNullBc_T)
                        ! this BC is already an MPI BC
                        if(Grid%hasUpperBC(KDIR) .and. .not. periodicK) then
                            print *, 'Max K BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                            stop
                        endif
                    CLASS DEFAULT
                        ! this BC is something other than an MPI BC
                        if(.not. Grid%hasUpperBC(KDIR) .or. periodicK) then
                            print *, 'Over-writing max K BC to be an MPI BC'
                            deallocate(Grid%externalBCs(OUTK)%p)
                            allocate(mpiNullBc_T :: Grid%externalBCs(OUTK)%p)
                        endif
                END SELECT

            endif

            ! perform a halo update before the sim starts to ensure that the ghost cells have correct values
            call haloUpdate(gamAppMpi)

            !Update the old state
            gamAppMpi%oState = gamAppMpi%State

            !Ensure all processes have the same starting timestep
            tmpDT = Model%dt
            call MPI_AllReduce(MPI_IN_PLACE, tmpDT, 1, MPI_MYFLOAT, MPI_MIN, gamAppMpi%gamMpiComm,ierr)
            Model%dt = tmpDT
            gamAppMpi%oState%time = gamAppMpi%State%time-Model%dt !Initial old state

        endif

        end associate
    end subroutine Hatch_mpi

    subroutine stepGamera_mpi(gamAppMpi)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi

        integer :: ierr
        real(rp) :: tmp

        !update the state variables to the next timestep
        call UpdateStateData(gamAppMpi)

        !Update ghost cells
        call Tic("Halos")
        call HaloUpdate(gamAppMpi)
        call Toc("Halos")

        !Calculate new timestep
        call Tic("DT")
        gamAppMpi%Model%dt = CalcDT(gamAppMpi%Model,gamAppMpi%Grid,gamAppMpi%State)

        !All MPI ranks take the lowest dt
        call Tic("mpiDT")
        tmp = gamAppMpi%Model%dt
        call MPI_AllReduce(MPI_IN_PLACE, tmp, 1, MPI_MYFLOAT, MPI_MIN, gamAppMpi%gamMpiComm,ierr)
        gamAppMpi%Model%dt = tmp
        call Toc("mpiDT")

        call Toc("DT")

        !Enforce BCs
        call Tic("Halos")
        call EnforceBCs(gamAppMpi%Model,gamAppMpi%Grid,gamAppMpi%State)
        call Toc("Halos")

    end subroutine stepGamera_mpi

    subroutine haloUpdate(gamAppMpi)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi

        integer :: ierr, length
        character(len=strLen) :: message

        if(gamAppMpi%gamMpiComm /= MPI_COMM_NULL) then
            ! just tell MPI to use the arrays we defined during initialization to send and receive data!
            ! GAS data
            call mpi_neighbor_alltoallw(gamAppMpi%State%Gas, gamAppMpi%sendCountsGas, gamAppMpi%sendDisplsGas, gamAppMpi%sendTypesGas, &
                                        gamAppMpi%State%Gas, gamAppMpi%recvCountsGas, gamAppMpi%recvDisplsGas, gamAppMpi%recvTypesGas, &
                                        gamAppMpi%gamMpiComm, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif
        endif

    end subroutine haloUpdate

end module gamapp_mpi

