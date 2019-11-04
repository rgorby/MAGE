! Main data objects and functions to perform a gamera simulation

module gamapp_mpi
    use gamtypes
    use step
    use init
    use mhdgroup
    use gamapp
    use bcs_mpi
    use mpi

    implicit none

    type, extends(GamApp_T) :: gamAppMpi_T
        integer :: NumRi=1,NumRj=1,NumRk=1
        integer :: Ri=1,Rj=1,Rk=1
        type(MPI_Comm) :: gamMpiComm
        integer :: sendRanks(:), recvRanks(:)
        integer :: sendCountsGas(:), sendDisplsGas(:), sendTypesGas(:)
        integer :: recvCountsGas(:), recvDisplsGas(:), recvTypesGas(:)
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
            call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")
        else
            !Find input deck
            call getIDeckStr(inpXML)
        endif

        if (present(doIO)) then
            doIOX = doIO
        else
            doIOX = .true.
        endif

        write(*,*) 'Reading input deck from ', trim(inpXML)
        inquire(file=inpXML,exist=fExist)
        if (.not. fExist) then
            write(*,*) 'Error opening input deck, exiting ...'
            write(*,*) ''
            stop
        endif
        
        !Create XML reader
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

        integer :: numNeighbors, ierr, length, commSize, rank, ic, jc, kc, targetRank, sendNumber, numInNeighbors, numOutNeighbors, localIndex
        integer, dimension(27), allocatable :: sourceRanks, sourceData, destRanks, destData
        logical :: reorder,periodicI,periodicJ,periodicK,wasWeighted
        character(len=strLen) :: message, bcType
        real(rp), dimension(:,:,:), allocatable :: tempX,tempY,tempZ

        ! call appropriate subroutines to read corner info and mesh size data
        call ReadCorners(gamAppMpi%Model,gamAppMpi%Grid,xmlInp,endTime)

        ! setup the distributed graph topology communicator for future MPI work
        call xmlInp%Set_Val(gamAppMpi%NumRi,'iPdir/N',1)
        call xmlInp%Set_Val(gamAppMpi%NumRj,'jPdir/N',1)
        call xmlInp%Set_Val(gamAppMpi%NumRk,'kPdir/N',1)

        call mpi_comm_size(gamComm, commSize, ierr)
        if(commSize /= gamAppMpi%NumRi*gamAppMpi%NumRj*gamAppMpi%NumRk) then
            print *,'Expected ',gamAppMpi%NumRi*gamAppMpi%NumRj*gamAppMpi%NumRk,' ranks, but this job only has ',commSize,' gamera ranks'
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif

        if(modulo(gamAppMpi%Grid%Nip,gamAppMpi%NumRi) /= 0) then
            print *,'Number of cells in i dimension is ',gamAppMpi%Grid%Nip,' but it must be a multiple of the number of MPI ranks in that dimension, which was ',gamAppMpi%NumRi
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif
        if(modulo(gamAppMpi%Grid%Njp,gamAppMpi%NumRj) /= 0) then
            print *,'Number of cells in j dimension is ',gamAppMpi%Grid%Njp,' but it must be a multiple of the nu
mber of MPI ranks in that dimension, which was ',gamAppMpi%NumRj
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif
        if(modulo(gamAppMpi%Grid%Nkp,gamAppMpi%NumRk) /= 0) then
            print *,'Number of cells in k dimension is ',gamAppMpi%Grid%Nkp,' but it must be a multiple of the nu
mber of MPI ranks in that dimension, which was ',gamAppMpi%NumRk
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif

        ! rank assignment. These are now 0 to N-1, matching the way that MPI counts it
        call mpi_comm_rank(gamComm, rank, ierr)
        gamAppMpi%Rk = modulo(rank, gamAppMpi%NumRk)
        rank = (rank-gamAppMpi%Rk)/gamAppMpi%NumRk
        gamAppMpi%Rj = modulo(rank, gamAppMpi%NumRj)
        rank = (rank-gamAppMpi%Rj)/gamAppMpi%NumRj
        gamAppMpi%Ri = rank

        ! adjust corner information to reflect this individual node's grid data
        gamAppMpi%Grid%Nip = gamAppMpi%Grid%Nip/gamAppMpi%NumRi
        gamAppMpi%Grid%Njp = gamAppMpi%Grid%Njp/gamAppMpi%NumRj
        gamAppMpi%Grid%Nkp = gamAppMpi%Grid%Nkp/gamAppMpi%NumRk

        gamAppMpi%Grid%Ni = gamAppMpi%Grid%Nip + 2*gamAppMpi%Model%nG
        gamAppMpi%Grid%Nj = gamAppMpi%Grid%Njp + 2*gamAppMpi%Model%nG
        gamAppMpi%Grid%Nk = gamAppMpi%Grid%Nkp + 2*gamAppMpi%Model%nG

        ! check which dimensions are using MPI for periodicity
        call xmlInp%Set_Val(bcType,'ibc/bc')
        if (bcType == 'periodic') then
            periodicI = .true.
        else
            periodicI = .false.
        endif
        call xmlInp%Set_Val(bcType,'jbc/bc')
        if (bcType == 'periodic') then
            periodicJ = .true.
        else
            periodicJ = .false.
        endif
        call xmlInp%Set_Val(bcType,'kbc/bc')
        if (bcType == 'periodic') then
            periodicK = .true.
        else
            periodicK = .false.
        endif

        ! create MPI topology
        numNeighbors = 0
        do ic=-1,1
            if((gamAppMpi%Ri+ic >= 0 .and. gamAppMpi%Ri+ic < GamAppMpi%NumRi) .or. periodicI) then
                do jc=-1,1
                    if((gamAppMpi%Rj+jc >= 0 .and. gamAppMpi%Rj+jc < GamAppMpi%NumRj) .or. periodicJ) then
                        do kc=-1,1
                            if((gamAppMpi%Rk+kc >= 0 .and. gamAppMpi%Rk+kc < GamAppMpi%NumRk) .or. periodicK)then
                                if(ic /= 0 .or. jc /= 0 .or. kc /= 0) then ! ensure I'm not talking to myself
                                    targetRank = modulo(gamAppMpi%Ri+ic,gamAppMpi%NumRi)*gamAppMpi%NumRk*gamAppMpi%NumRj + 
                                                 modulo(gamAppMpi%Rj+jc,gamAppMpi%NumRj)*gamAppMpi%NumRk +
                                                 modulo(gamAppMpi%Rk+kc,gamAppMpi%NumRk)
                                    listIndex = findloc(sourceRanks, targetRank)
                                    if(listIndex == 0) then ! this rank not in the list yet
                                        numNeighbors = numNeighbors+1
                                        listIndex = numNeighbors
                                        sourceRanks(listIndex) = targetRank
                                    endif
                                    ! calculate the size of this region that will be transmitted
                                    sendNumber = 1
                                    if (ic == 0) then
                                        sendNumber = sendNumber * gamAppMpi%Grid%Nip
                                    else
                                        sendNumber = sendNumber * gamAppMpi%Model%nG
                                    endif
                                    if (jc == 0) then
                                        sendNumber = sendNumber * gamAppMpi%Grid%Njp
                                    else
                                        sendNumber = sendNumber * gamAppMpi%Model%nG
                                    endif
                                    if (kc == 0) then
                                        sendNumber = sendNumber * gamAppMpi%Grid%Nkp
                                    else
                                        sendNumber = sendNumber * gamAppMpi%Model%nG
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
        call mpi_dist_graph_create_adjacent(gamComm,
            numNeighbors,sourceRanks,sourceData,
            numNeighbors,sourceRanks,sourceData, ! comms symmetrical
            mpiInfo, reorder, gamAppMpi%gamMpiComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! get our final ranks from the new MPI topology
        call mpi_comm_rank(gamAppMpi%gamMpiComm, rank, ierr)
        gamAppMpi%Rk = modulo(rank, gamAppMpi%NumRk)
        rank = (rank-gamAppMpi%Rk)/gamAppMpi%NumRk
        gamAppMpi%Rj = modulo(rank, gamAppMpi%NumRj)
        rank = (rank-gamAppMpi%Rj)/gamAppMpi%NumRj
        gamAppMpi%Ri = rank

        ! adjust grid info for these ranks
        gamAppMpi%Grid%ijkShift(1) = gamAppMpi%Grid%Nip*gamAppMpi%Ri
        gamAppMpi%Grid%ijkShift(2) = gamAppMpi%Grid%Njp*gamAppMpi%Rj
        gamAppMpi%Grid%ijkShift(3) = gamAppMpi%Grid%Nkp*gamAppMpi%Rk

        gamAppMpi%Grid%is = 1; gamAppMpi%Grid%ie = gamAppMpi%Grid%Nip
        gamAppMpi%Grid%js = 1; gamAppMpi%Grid%je = gamAppMpi%Grid%Njp
        gamAppMpi%Grid%ks = 1; gamAppMpi%Grid%ke = gamAppMpi%Grid%Nkp

        gamAppMpi%Grid%isg = gamAppMpi%Grid%is-gamAppMpi%Model%nG
        gamAppMpi%Grid%ieg = gamAppMpi%Grid%ie+gamAppMpi%Model%nG

        gamAppMpi%Grid%jsg = gamAppMpi%Grid%js-gamAppMpi%Model%nG
        gamAppMpi%Grid%jeg = gamAppMpi%Grid%je+gamAppMpi%Model%nG

        gamAppMpi%Grid%ksg = gamAppMpi%Grid%ks-gamAppMpi%Model%nG
        gamAppMpi%Grid%keg = gamAppMpi%Grid%ke+gamAppMpi%Model%nG

        ! create temporary arrays to hold this rank's subset of the full corner array
        allocate(tempX(gamAppMpi%Grid%isg:gamAppMpi%Grid%ieg+1,gamAppMpi%Grid%jsg:gamAppMpi%Grid%jeg+1,gamAppMpi%Grid%ksg:gamAppMpi%Grid%keg+1))
        allocate(tempY(gamAppMpi%Grid%isg:gamAppMpi%Grid%ieg+1,gamAppMpi%Grid%jsg:gamAppMpi%Grid%jeg+1,gamAppMpi%Grid%ksg:gamAppMpi%Grid%keg+1))
        allocate(tempZ(gamAppMpi%Grid%isg:gamAppMpi%Grid%ieg+1,gamAppMpi%Grid%jsg:gamAppMpi%Grid%jeg+1,gamAppMpi%Grid%ksg:gamAppMpi%Grid%keg+1))

        ! pull out this rank's relevant corner info
        tempX = gamAppMpi%Grid%x(gamAppMpi%Grid%isg+gamAppMpi%Grid%ijkShift(1):gamAppMpi%Grid%ieg+gamAppMpi%Grid%ijkShift(1), gamAppMpi%Grid%jsg+gamAppMpi%Grid%ijkShift(2):gamAppMpi%Grid%jeg+gamAppMpi%Grid%ijkShift(2), gamAppMpi%Grid%ksg+gamAppMpi%Grid%ijkShift(3):gamAppMpi%Grid%keg+gamAppMpi%Grid%ijkShift(3))
        tempY = gamAppMpi%Grid%y(gamAppMpi%Grid%isg+gamAppMpi%Grid%ijkShift(1):gamAppMpi%Grid%ieg+gamAppMpi%Grid%ijkShift(1), gamAppMpi%Grid%jsg+gamAppMpi%Grid%ijkShift(2):gamAppMpi%Grid%jeg+gamAppMpi%Grid%ijkShift(2), gamAppMpi%Grid%ksg+gamAppMpi%Grid%ijkShift(3):gamAppMpi%Grid%keg+gamAppMpi%Grid%ijkShift(3))
        tempZ = gamAppMpi%Grid%z(gamAppMpi%Grid%isg+gamAppMpi%Grid%ijkShift(1):gamAppMpi%Grid%ieg+gamAppMpi%Grid%ijkShift(1), gamAppMpi%Grid%jsg+gamAppMpi%Grid%ijkShift(2):gamAppMpi%Grid%jeg+gamAppMpi%Grid%ijkShift(2), gamAppMpi%Grid%ksg+gamAppMpi%Grid%ijkShift(3):gamAppMpi%Grid%keg+gamAppMpi%Grid%ijkShift(3))

        ! delete the old corner arrays
        deallocate(gamAppMpi%Grid%x)
        deallocate(gamAppMpi%Grid%y)
        deallocate(gamAppMpi%Grid%z)

        ! move the new arrays to the grid corner arrays
        call move_alloc(tempX, gamAppMpi%Grid%x)
        call move_alloc(tempY, gamAppMpi%Grid%y)
        call move_alloc(tempZ, gamAppMpi%Grid%z)

        ! now create the arrays that MPI will use to send and receive the data
        call mpi_dist_graph_neighbors_count(gamAppMpi%gamMpiComm, numInNeighbors, numOutNeighbors, wasWeighted, ierr)
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
        call mpi_dist_graph_neighbors(gamAppMpi%gamMpiComm, numInNeighbors, gamAppMpi%recvRanks, sendData, numOutNeighbors, gamAppMpi%sendRanks, sendData, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        allocate(gamAppMpi%sendCountsGas(numOutNeighbors))
        allocate(gamAppMpi%sendDisplsGas(numOutNeighbors))
        allocate(gamAppMpi%sendTypesGas(numOutNeighbors))
        allocate(gamAppMpi$recvCountsGas(numInNeighbors))
        allocate(gamAppMpi%recvDisplsGas(numInNeighbors))
        allocate(gamAppMpi%recvTypesgas(numInNeighbors))

        ! assemble the different datatypes for each rank, and then save the resulting final datatype
        call mpi_type_extent(MPI_MYFLOAT, dataSize, ierr) ! number of bytes per array entry

        ! I dimension
        call mpi_type_contiguous(gamAppMpi%Model%nG, MPI_MYFLOAT, iG, ierr) ! ghosts i
        call mpi_type_contiguous(gamAppMpi%Model%Nip, MPI_MYFLOAT, iP, ierr) ! physical i

        ! J dimension
        call mpi_type_hvector(gamAppMpi%Model%nG, 1, gamAppMpi%Model%Ni*dataSize, iG, iGjG, ierr) ! ghosts i - ghosts j
        call mpi_type_hvector(gamAppMpi%Model%Njp, 1, gamAppMpi%Model%Ni*dataSize, iG, iGjP, ierr) ! ghosts i - physical j
        call mpi_type_hvector(gamAppMpi%Model%nG, 1, gamAppMpi%Model%Ni*dataSize, iP, iPjG, ierr) ! physical i - ghosts j
        call mpi_type_hvector(gamAppMpi%Model%Njp, 1, gamAppMpi%Model%Ni*dataSize, iP, iPjP, ierr) ! physical i - physical j

        ! K dimension
        call mpi_type_hvector(gamAppMpi%Model%nG, 1, gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*dataSize, iGjG, cornerMpiType, ierr) ! ghosts i - ghosts j - ghosts k
        call mpi_type_hvector(gamAppMpi%Model%nG, 1,  gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*dataSize, iGjP, jEdgeMpiType, ierr) ! ghosts i - physical j - ghosts k
        call mpi_type_hvector(gamAppMpi%Model%nG, 1,  gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*dataSize, iPjG, iEdgeMpiType, ierr) ! physical i - ghosts j - ghosts k
        call mpi_type_hvector(gamAppMpi%Model%Nkp, 1,  gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*dataSize, iGjG, kEdgeMpiType, ierr) ! ghosts i - ghosts j - physical k
        call mpi_type_hvector(gamAppMpi%Model%nG, 1,  gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*dataSize, iPjP, kFaceMpiType, ierr) ! physical i - physical j - ghosts k
        call mpi_type_hvector(gamAppMpi%Model%Nkp, 1,  gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*dataSize, iGjP, iFaceMpiType, ierr) ! ghosts i - physical j - physical k
        call mpi_type_hvector(gamAppMpi%Model%Nkp, 1,  gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*dataSize, iPjG, jFaceMpiType, ierr) ! physical i - ghosts j - physical k

        ! 4th dimension
        call mpi_type_hvector(NVAR, 1, gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*gamAppMpi%Model%Nk*dataSize, cornerMpiType, corner4MpiType, ierr)
        call mpi_type_hvector(NVAR, 1, gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*gamAppMpi%Model%Nk*dataSize, iEdgeMpiType, iEdge4MpiType, ierr)
        call mpi_type_hvector(NVAR, 1, gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*gamAppMpi%Model%Nk*dataSize, jEdgeMpiType, jEdge4MpiType, ierr)
        call mpi_type_hvector(NVAR, 1, gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*gamAppMpi%Model%Nk*dataSize, kEdgeMpiType, kEdge4MpiType, ierr)
        call mpi_type_hvector(NVAR, 1, gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*gamAppMpi%Model%Nk*dataSize, iFaceMpiType, iFace4MpiType, ierr)
        call mpi_type_hvector(NVAR, 1, gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*gamAppMpi%Model%Nk*dataSize, jFaceMpiType, jFace4MpiType, ierr)
        call mpi_type_hvector(NVAR, 1, gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*gamAppMpi%Model%Nk*dataSize, kFaceMpiType, kFace4MpiType, ierr)

        ! 5th dimension
        call mpi_type_hvector(gamAppMpi%Model%nSpc+1, 1, NVAR*gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*gamAppMpi%Model%Nk*dataSize, corner4MpiType, corner5MpiType, ierr)
        call mpi_type_hvector(gamAppMpi%Model%nSpc+1, 1, NVAR*gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*gamAppMpi%Model%Nk*dataSize, iEdge4MpiType, iEdge5MpiType, ierr)
        call mpi_type_hvector(gamAppMpi%Model%nSpc+1, 1, NVAR*gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*gamAppMpi%Model%Nk*dataSize, jEdge4MpiType, jEdge5MpiType, ierr)
        call mpi_type_hvector(gamAppMpi%Model%nSpc+1, 1, NVAR*gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*gamAppMpi%Model%Nk*dataSize, kEdge4MpiType, kEdge5MpiType, ierr)
        call mpi_type_hvector(gamAppMpi%Model%nSpc+1, 1, NVAR*gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*gamAppMpi%Model%Nk*dataSize, iFace4MpiType, iFace5MpiType, ierr)
        call mpi_type_hvector(gamAppMpi%Model%nSpc+1, 1, NVAR*gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*gamAppMpi%Model%Nk*dataSize, jFace4MpiType, jFace5MpiType, ierr)
        call mpi_type_hvector(gamAppMpi%Model%nSpc+1, 1, NVAR*gamAppMpi%Model%Ni*gamAppMpi%Model%Nj*gamAppMpi%Model%Nk*dataSize, kFace4MpiType, kFace5MpiType, ierr)

        ! counts are always 1 because we're sending a single (complicated) mpi datatype
        sendCountsGas(:) = 1
        recvCountsGas(:) = 1

        sourceRanks(:) = 0
        sourceData(:) = MPI_DATATYPE_NULL
        destRanks(:) = 0
        destData(:) = MPI_DATATYPE_NULL
        ! figure out exactly what data needs to be sent to (and received from) each neighbor, and create custom MPI datatypes to perform these transfers
        do ic=-1,1
            if((gamAppMpi%Ri+ic >= 0 .and. gamAppMpi%Ri+ic < GamAppMpi%NumRi) .or. periodicI) then
                do jc=-1,1
                    if((gamAppMpi%Rj+jc >= 0 .and. gamAppMpi%Rj+jc < GamAppMpi%NumRj) .or. periodicJ) then
                        do kc=-1,1
                            if((gamAppMpi%Rk+kc >= 0 .and. gamAppMpi%Rk+kc < GamAppMpi%NumRk) .or. periodicK)then
                                if(ic /= 0 .or. jc /= 0 .or. kc /= 0) then ! ensure I'm not talking to myself
                                    targetRank = modulo(gamAppMpi%Ri+ic,gamAppMpi%NumRi)*gamAppMpi%NumRk*gamAppMpi%NumRj +
                                                 modulo(gamAppMpi%Rj+jc,gamAppMpi%NumRj)*gamAppMpi%NumRk +
                                                 modulo(gamAppMpi%Rk+kc,gamAppMpi%NumRk)

                                    localIndex = 1 + (ic+1)*9 + (jc+1)*3 + (kc+1)
                                    sourceRanks(localIndex) = targetRank
                                    destRanks(localIndex) = targetRank

                                    sendDataOffset = 0
                                    recvDataOffset = 0
                                    SELECT CASE (ic)
                                        case (-1)
                                            ! min i side
                                            ! no change to sendDataOffset
                                            recvDataOffset = recvDataOffset + gamAppMpi%Model%nG
                                        case (0)
                                            ! central in i dimension
                                            sendDataOffset = sendDataOffset + gamAppMpi%Model%nG
                                            recvDataOffset = recvDataOffset + gamAppMpi%Model%nG
                                        case (1)
                                            ! max i side
                                            sendDataOffset = sendDataOffset + gamAppMpi%Model%nG + gamAppMpi%Grid%Nip
                                            recvDataOffset = recvDataOffset + gamAppMpi%Grid%Nip ! ghosts+physical-ghosts
                                        CASE DEFAULT
                                            print *, 'Unexpected ic direction'
                                            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                                    END SELECT

                                    SELECT CASE (jc)
                                        case (-1)
                                            ! min j side
                                            ! no change to sendDataOffset
                                            recvDataOffset = recvDataOffset + gamAppMpi%Model%nG*gamAppMpi%Grid%Ni
                                        case (0)
                                            ! central in j dimension
                                            sendDataOffset = sendDataOffset + gamAppMpi%Model%nG*gamAppMpi%Grid%Ni
                                            recvDataOffset = recvDataOffset + gamAppMpi%Model%nG*gamAppMpi%Grid%Ni
                                        case (1)
                                            ! max j side
                                            sendDataOffset = sendDataOffset + (gamAppMpi%Model%nG + gamAppMpi%Grid%Njp)*gamAppMpi%Grid%Ni
                                            recvDataOffset = recvDataOffset + gamAppMpi%Grid%Njp*gamAppMpi%Grid%Ni
                                        CASE DEFAULT
                                            print *, 'Unexpected jc direction'
                                            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                                    END SELECT

                                    SELECT CASE (kc)
                                        case (-1)
                                            ! min k side
                                            ! no change to sendDataOffset
                                            recvDataOffset = recvDataOffset + gamAppMpi%Model%nG*gamAppMpi%Grid%Ni*gamAppMpi%Grid%Nj
                                        case (0)
                                            ! central in k dimension
                                            sendDataOffset = sendDataOffset + gamAppMpi%Model%nG*gamAppMpi%Grid%Ni*gamAppMpi%Grid%Nj
                                            recvDataOffset = recvDataOffset + gamAppMpi%Model%nG*gamAppMpi%Grid%N
i*gamAppMpi%Grid%Nj
                                        case (1)
                                            ! max k side
                                            sendDataOffset = sendDataOffset + (gamAppMpi%Model%nG + gamAppMpi%Grid%Nkp)*gamAppMpi%Grid%Ni*gamAppMpi%Grid%Nj
                                            recvDataOffset = recvDataOffset + gamAppMpi%Grid%Nkp*gamAppMpi%Grid%Ni*gamAppMpi%Grid%Nj
                                        CASE DEFAULT
                                            print *, 'Unexpected kc direction'
                                            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                                    END SELECT

                                    ! determine which of the previously created datatypes is the correct one
                                    SELECT CASE (abs(ic)+abs(jc)+abs(kc))
                                        case (1) ! face
                                            if(ic ~= 0) then
                                                transDataType = iFace5MpiType
                                            elseif(jc ~= 0) then
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
                                            transDataType = corner
                                        CASE DEFAULT
                                            print *, 'Sum of ic+jc+kc is nonsense'
                                            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                                    END SELECT

                                    if(sourceData(localIndex) == MPI_DATATYPE_NULL) then
                                        ! not sending any data to this rank yet, just add this datatype
                                        call mpi_type_hindexed(1, 1, sendDataOffset*dataSize, transDataType, sourceData(localIndex), ierr)
                                        call mpi_type_hindexed(1, 1, recvDataOffset*dataSize, transDataType, recvData(localIndex), ierr)
                                    else
                                        ! we're already sending other data to this rank
                                        !  merge the datatypes into a struct
                                        call mpi_type_create_struct(2, [1 1], [0 sendDataOffset*dataSize] [sourceData(localIndex) transDataType], sourceData(localIndex), ierr)
                                        call mpi_type_create_struct(2, [1 1], [0 recvDataOffset*dataSize] [destData(localIndex) transDataType], destData(localIndex), ierr)
                                    endif

                                endif
                            endif
                        enddo
                    endif
                enddo
            endif
        enddo

        ! call appropriate subroutines to calculate all appropriate grid data from the corner data
        call CalcGridInfo(gamAppMpi%Model,gamAppMpi%Grid,gamAppMpi%State,gamAppMpi%oState,gamAppMpi%Solver,xmlInp,userInitFunc)

        ! correct boundary conditions if necessary
        if(gamAppMpi%NumRi > 1) then
            ! MPI decomposed in I dimension

            ! check min I BC
            SELECT type(iiBC=>gamAppMpi%Grid%externalBCs(INI)%p)
                TYPE IS (mpiNullBc_T)
                    ! this BC is already an MPI BC
                    if((gamAppMpi%Ri == 0 .and. .not. periodicI) then
                        print *, 'Min I BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                        stop
                    endif
                CLASS DEFAULT
                    ! this BC is something other than an MPI BC
                    if(gamAppMpi%Ri > 0 .or. periodicI) then
                        print *, 'Over-writing min I BC to be an MPI BC'
                        deallocate(gamAppMpi%Grid%externalBCs(INI)%p)
                        allocate(mpiNullBc_T :: gamAppMpi%Grid%externalBCs(INI)%p)
                    endif
            END SELECT

            ! check max I BC
            SELECT type(iiBC=>gamAppMpi%Grid%externalBCs(OUTI)%p)
                TYPE IS (mpiNullBc_T)
                    ! this BC is already an MPI BC
                    if((gamAppMpi%Ri == (gamAppMpi%NumRi-1) .and. .not. periodicI) then
                        print *, 'Max I BC was set to be an MPI BC by the user, but this is not an MPI periodic c
ase'
                        stop
                    endif
                CLASS DEFAULT
                    ! this BC is something other than an MPI BC
                    if(gamAppMpi%Ri < (gamAppMpi%NumRi-1) .or. periodicI) then
                        print *, 'Over-writing max I BC to be an MPI BC'
                        deallocate(gamAppMpi%Grid%externalBCs(OUTI)%p)
                        allocate(mpiNullBc_T :: gamAppMpi%Grid%externalBCs(OUTI)%p)
                    endif
            END SELECT

        endif

        if(gamAppMpi%NumRj > 1) then
            ! MPI decomposed in J dimension

            ! check min J BC
            SELECT type(iiBC=>gamAppMpi%Grid%externalBCs(INJ)%p)
                TYPE IS (mpiNullBc_T)
                    ! this BC is already an MPI BC
                    if((gamAppMpi%Rj == 0 .and. .not. periodicJ) then
                        print *, 'Min J BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                        stop
                    endif
                CLASS DEFAULT
                    ! this BC is something other than an MPI BC
                    if(gamAppMpi%Rj > 0 .or. periodicJ) then
                        print *, 'Over-writing min J BC to be an MPI BC'
                        deallocate(gamAppMpi%Grid%externalBCs(INJ)%p)
                        allocate(mpiNullBc_T :: gamAppMpi%Grid%externalBCs(INJ)%p)
                    endif
            END SELECT

            ! check max J BC
            SELECT type(iiBC=>gamAppMpi%Grid%externalBCs(OUTJ)%p)
                TYPE IS (mpiNullBc_T)
                    ! this BC is already an MPI BC
                    if((gamAppMpi%Rj == (gamAppMpi%NumRj-1) .and. .not. periodicJ) then
                        print *, 'Max J BC was set to be an MPI BC by the user, but this is not an MPI periodic c
ase'
                        stop
                    endif
                CLASS DEFAULT
                    ! this BC is something other than an MPI BC
                    if(gamAppMpi%Rj < (gamAppMpi%NumRj-1) .or. periodicJ) then
                        print *, 'Over-writing max J BC to be an MPI BC'
                        deallocate(gamAppMpi%Grid%externalBCs(OUTJ)%p)
                        allocate(mpiNullBc_T :: gamAppMpi%Grid%externalBCs(OUTJ)%p)
                    endif
            END SELECT

        endif

        if(gamAppMpi%NumRk > 1) then
            ! MPI decomposed in K dimension

            ! check min K BC
            SELECT type(iiBC=>gamAppMpi%Grid%externalBCs(INK)%p)
                TYPE IS (mpiNullBc_T)
                    ! this BC is already an MPI BC
                    if((gamAppMpi%Rk == 0 .and. .not. periodicK) then
                        print *, 'Min K BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                        stop
                    endif
                CLASS DEFAULT
                    ! this BC is something other than an MPI BC
                    if(gamAppMpi%Rk > 0 .or. periodicK) then
                        print *, 'Over-writing min K BC to be an MPI BC'
                        deallocate(gamAppMpi%Grid%externalBCs(INK)%p)
                        allocate(mpiNullBc_T :: gamAppMpi%Grid%externalBCs(INK)%p)
                    endif
            END SELECT

            ! check max K BC
            SELECT type(iiBC=>gamAppMpi%Grid%externalBCs(OUTK)%p)
                TYPE IS (mpiNullBc_T)
                    ! this BC is already an MPI BC
                    if((gamAppMpi%Rk == (gamAppMpi%NumRk-1) .and. .not. periodicK) then
                        print *, 'Max K BC was set to be an MPI BC by the user, but this is not an MPI periodic c
ase'
                        stop
                    endif
                CLASS DEFAULT
                    ! this BC is something other than an MPI BC
                    if(gamAppMpi%Rk < (gamAppMpi%NumRk-1) .or. periodicK) then
                        print *, 'Over-writing max K BC to be an MPI BC'
                        deallocate(gamAppMpi%Grid%externalBCs(OUTK)%p)
                        allocate(mpiNullBc_T :: gamAppMpi%Grid%externalBCs(OUTK)%p)
                    endif
            END SELECT

        endif

    end subroutine Hatch_mpi

    subroutine stepGamera_mpi(gamAppMpi)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi

        !update the state variables to the next timestep
        call UpdateStateData(gamAppMpi)

        !Update ghost cells
        call Tic("Halos")
        call HaloUpdate(gamAppMpi)
        call Toc("Halos")

        ! calculate new DT, update BCs, write output data, etc...
        call FinishStep(gamAppMpi)

    end subroutine stepGamera_mpi

    subroutine haloUpdate(gamAppMpi)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi

        integer :: ierr

        ! just tell MPI to use the arrays we defined during initialization to send and receive data!

        ! GAS data
        call mpi_neighbor_alltoallw(gamAppMpi%State%Gas, gamAppMpi%sendCountsGas, gamAppMpi%sendDisplsGas, gamAppMpi%sendTypesGas,
                                    gamAppMpi%State%Gas, gamAppMpi%recvCountsGas, gamAppMpi%recvDisplsGas, gamAppMpi%recvTypesGas,
                                    gamAppMpi%gamMpiComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine haloUpdate

end module gamapp_mpi

