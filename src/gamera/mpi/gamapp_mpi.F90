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

        ! Gas Data Transfer Variables
        integer, dimension(:), allocatable :: sendCountsGas, sendTypesGas
        integer, dimension(:), allocatable :: recvCountsGas, recvTypesGas
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: sendDisplsGas, recvDisplsGas

        ! Magnetic Flux Data Transfer Variables
        integer, dimension(:), allocatable :: sendCountsMagFlux, sendTypesMagFlux
        integer, dimension(:), allocatable :: recvCountsMagFlux, recvTypesMagFlux
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: sendDisplsMagFlux, recvDisplsMagFlux

        ! Magnetic Field Data Transfer Variables
        integer, dimension(:), allocatable :: sendCountsBxyz, sendTypesBxyz
        integer, dimension(:), allocatable :: recvCountsBxyz, recvTypesBxyz
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: sendDisplsBxyz, recvDisplsBxyz
        
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

        integer :: numNeighbors, ierr, length, commSize, rank, ic, jc, kc
        integer :: targetRank, sendNumber, numInNeighbors, numOutNeighbors, listIndex
        integer, dimension(27) :: sourceRanks, sourceData
        logical :: reorder,periodicI,periodicJ,periodicK,wasWeighted
        character(len=strLen) :: message, bcType
        real(rp), dimension(:,:,:), allocatable :: tempX,tempY,tempZ
        real(rp) :: tmpDT

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
                do jc=-1,1
                    do kc=-1,1
                        targetRank = 0

                        ! I rank offset
                        if((Grid%Ri+ic >= 0 .and. Grid%Ri+ic < Grid%NumRi) .or. periodicI) then
                            ! talking to a neighbor, or wrapping around periodic MPI boundaries
                            targetRank = targetRank + modulo(Grid%Ri+ic,Grid%NumRi)*Grid%NumRk*Grid%NumRj
                        else
                            ! talking to a non-periodic boundary, so this is my own I rank
                            targetRank = targetRank + Grid%Ri*Grid%NumRk*Grid%NumRj
                        endif

                        ! J rank offset
                        if((Grid%Rj+jc >= 0 .and. Grid%Rj+jc < Grid%NumRj) .or. periodicJ) then
                            ! talking to a neighbor, or wrapping around periodic MPI boundaries
                            targetRank = targetRank + modulo(Grid%Rj+jc,Grid%NumRj)*Grid%NumRk
                        else
                            ! talking to a non-periodic boundary, so this is my own J rank
                            targetRank = targetRank + Grid%Rj*Grid%NumRk
                        endif

                        ! K rank offset
                        if((Grid%Rk+kc >= 0 .and. Grid%Rk+kc < Grid%NumRk) .or. periodicK) then
                            ! talking to a neighbor, or wrapping around periodic MPI boundaries
                            targetRank = targetRank + modulo(Grid%Rk+kc,Grid%NumRk)
                        else
                            ! talking to a non-periodic boundary, so this is my own K rank
                            targetRank = targetRank + Grid%Rk
                        endif
                       
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
                    enddo
                enddo
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

            ! Create MPI datatypes for cell centered arrays
            call createCellCenteredDatatypes(gamAppMpi, periodicI, periodicJ, periodicK)

            ! Create MPI datatypes for face centered arrays
            call createFaceCenteredDatatypes(gamAppMpi, periodicI, periodicJ, periodicK)

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
            Model%dt = CalcDT(Model,Grid,gamAppMpi%State)
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
        call Tic("BCs")
        call EnforceBCs(gamAppMpi%Model,gamAppMpi%Grid,gamAppMpi%State)
        call Toc("BCs")

    end subroutine stepGamera_mpi

    subroutine haloUpdate(gamAppMpi)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi

        integer :: ierr, length
        character(len=strLen) :: message

        if(gamAppMpi%gamMpiComm /= MPI_COMM_NULL) then
            ! just tell MPI to use the arrays we defined during initialization to send and receive data!

            ! Gas Cell Data
            call mpi_neighbor_alltoallw(gamAppMpi%State%Gas, gamAppMpi%sendCountsGas, &
                                        gamAppMpi%sendDisplsGas, gamAppMpi%sendTypesGas, &
                                        gamAppMpi%State%Gas, gamAppMpi%recvCountsGas, &
                                        gamAppMpi%recvDisplsGas, gamAppMpi%recvTypesGas, &
                                        gamAppMpi%gamMpiComm, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif

            if(gamAppMpi%Model%doMHD) then
                ! Magnetic Face Flux Data
                call mpi_neighbor_alltoallw(gamAppMpi%State%magFlux, gamAppMpi%sendCountsMagFlux, &
                                            gamAppMpi%sendDisplsMagFlux, gamAppMpi%sendTypesMagFlux, &
                                            gamAppMpi%State%magFlux, gamAppMpi%recvCountsMagFlux, &
                                            gamAppMpi%recvDisplsMagFlux, gamAppMpi%recvTypesMagFlux, &
                                            gamAppMpi%gamMpiComm, ierr)
                if(ierr /= MPI_Success) then
                    call MPI_Error_string( ierr, message, length, ierr)
                    print *,message(1:length)
                    call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                endif

                ! Magnetic Field Cell Data
                call mpi_neighbor_alltoallw(gamAppMpi%State%Bxyz, gamAppMpi%sendCountsBxyz, &
                                            gamAppMpi%sendDisplsBxyz, gamAppMpi%sendTypesBxyz, &
                                            gamAppMpi%State%Bxyz, gamAppMpi%recvCountsBxyz, &
                                            gamAppMpi%recvDisplsBxyz, gamAppMpi%recvTypesBxyz, &
                                            gamAppMpi%gamMpiComm, ierr)
                if(ierr /= MPI_Success) then
                    call MPI_Error_string( ierr, message, length, ierr)
                    print *,message(1:length)
                    call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                endif
            endif

        endif

    end subroutine haloUpdate

    subroutine createCellCenteredDatatypes(gamAppMpi, periodicI, periodicJ, periodicK)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi
        logical, intent(in) :: periodicI, periodicJ, periodicK

        integer :: ierr, rank, ic, jc, kc, localIndexIn, localIndexOut, sendDataOffset, recvDataOffset
        integer :: dataSize, iG, iP, iGjG, iGjP, iPjG, iPjP, targetSendRank, targetRecvRank, gasType, fluxType
        integer :: cornerMpiType,iEdgeMpiType,jEdgeMpiType,kEdgeMpiType,iFaceMpiType,jFaceMpiType,kFaceMpiType
        integer :: corner4Gas,iEdge4Gas,jEdge4Gas,kEdge4Gas,iFace4Gas,jFace4Gas,kFace4Gas
        integer :: corner5Gas,iEdge5Gas,jEdge5Gas,kEdge5Gas,iFace5Gas,jFace5Gas,kFace5Gas
        integer(kind=MPI_ADDRESS_KIND) :: tempOffsets(2)

        associate(Grid=>gamAppMpi%Grid,Model=>gamAppMpi%Model)

        allocate(gamAppMpi%sendCountsGas(1:SIZE(gamAppMpi%sendRanks)))
        allocate(gamAppMpi%sendDisplsGas(1:SIZE(gamAppMpi%sendRanks)))
        allocate(gamAppMpi%sendTypesGas(1:SIZE(gamAppMpi%sendRanks)))
        allocate(gamAppMpi%recvCountsGas(1:SIZE(gamAppMpi%recvRanks)))
        allocate(gamAppMpi%recvDisplsGas(1:SIZE(gamAppMpi%recvRanks)))
        allocate(gamAppMpi%recvTypesGas(1:SIZE(gamAppMpi%recvRanks)))

        ! counts are always 1 because we're sending a single (complicated) mpi datatype
        gamAppMpi%sendCountsGas(:) = 1
        gamAppMpi%recvCountsGas(:) = 1

        ! displacements are always 0 because the displacements are baked into each mpi datatype
        gamAppMpi%sendDisplsGas(:) = 0
        gamAppMpi%recvDisplsGas(:) = 0

        ! set all datatypes to null by default
        gamAppMpi%sendTypesGas(:) = MPI_DATATYPE_NULL
        gamAppMpi%recvTypesGas(:) = MPI_DATATYPE_NULL

        ! get my rank
        call mpi_comm_rank(gamAppMpi%gamMpiComm, rank, ierr)

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
        call mpi_type_hvector(Model%nG, 1, Grid%Ni*Grid%Nj*dataSize, &
                              iGjG, cornerMpiType, ierr) ! ghosts i   - ghosts j   - ghosts k
        call mpi_type_hvector(Model%nG, 1, Grid%Ni*Grid%Nj*dataSize, &
                              iGjP, jEdgeMpiType,  ierr) ! ghosts i   - physical j - ghosts k
        call mpi_type_hvector(Model%nG, 1, Grid%Ni*Grid%Nj*dataSize, &
                              iPjG, iEdgeMpiType,  ierr) ! physical i - ghosts j   - ghosts k
        call mpi_type_hvector(Grid%Nkp, 1, Grid%Ni*Grid%Nj*dataSize, &
                              iGjG, kEdgeMpiType,  ierr) ! ghosts i   - ghosts j   - physical k
        call mpi_type_hvector(Model%nG, 1, Grid%Ni*Grid%Nj*dataSize, &
                              iPjP, kFaceMpiType,  ierr) ! physical i - physical j - ghosts k
        call mpi_type_hvector(Grid%Nkp, 1, Grid%Ni*Grid%Nj*dataSize, &
                              iGjP, iFaceMpiType,  ierr) ! ghosts i   - physical j - physical k
        call mpi_type_hvector(Grid%Nkp, 1, Grid%Ni*Grid%Nj*dataSize, &
        iPjG, jFaceMpiType,  ierr) ! physical i - ghosts j   - physical k

        ! 4th dimension for Gas
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, cornerMpiType, corner4Gas, ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iEdgeMpiType,  iEdge4Gas,  ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, jEdgeMpiType,  jEdge4Gas,  ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, kEdgeMpiType,  kEdge4Gas,  ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iFaceMpiType,  iFace4Gas,  ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, jFaceMpiType,  jFace4Gas,  ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, kFaceMpiType,  kFace4Gas,  ierr)

        ! 5th dimension for Gas
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, corner4Gas, corner5Gas, ierr)
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iEdge4Gas,  iEdge5Gas,  ierr)
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, jEdge4Gas,  jEdge5Gas,  ierr)
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, kEdge4Gas,  kEdge5Gas,  ierr)
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iFace4Gas,  iFace5Gas,  ierr)
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, jFace4Gas,  jFace5Gas,  ierr)
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, kFace4Gas,  kFace5Gas,  ierr)

        ! figure out exactly what data needs to be sent to (and received from) each neighbor
        ! create custom MPI datatypes to perform these transfers
        do ic=-1,1
            do jc=-1,1
                do kc=-1,1
                    ! send and receive types must be calculated in reverse order so that the data goes to the
                    ! correct location. Therefore, use ic,jc,kc for recv and -ic,-jc,-kc for send

                    targetRecvRank = 0
                    targetSendRank = 0

                    ! I rank offset
                    if((Grid%Ri+ic >= 0 .and. Grid%Ri+ic < Grid%NumRi) .or. periodicI) then
                        ! talking to a neighbor, or wrapping around periodic MPI boundaries
                        targetRecvRank = targetRecvRank + modulo(Grid%Ri+ic,Grid%NumRi)*Grid%NumRk*Grid%NumRj
                    else
                        ! talking to a non-periodic boundary, so this is my own I rank
                        targetRecvRank = targetRecvRank + Grid%Ri*Grid%NumRk*Grid%NumRj
                    endif
                    if((Grid%Ri-ic >= 0 .and. Grid%Ri-ic < Grid%NumRi) .or. periodicI) then
                        ! talking to a neighbor, or wrapping around periodic MPI boundaries
                        targetSendRank = targetSendRank + modulo(Grid%Ri-ic,Grid%NumRi)*Grid%NumRk*Grid%NumRj
                    else
                        ! talking to a non-periodic boundary, so this is my own I rank
                        targetSendRank = targetSendRank + Grid%Ri*Grid%NumRk*Grid%NumRj
                    endif

                    ! J rank offset
                    if((Grid%Rj+jc >= 0 .and. Grid%Rj+jc < Grid%NumRj) .or. periodicJ) then
                        ! talking to a neighbor, or wrapping around periodic MPI boundaries
                        targetRecvRank = targetRecvRank + modulo(Grid%Rj+jc,Grid%NumRj)*Grid%NumRk
                    else
                        ! talking to a non-periodic boundary, so this is my own J rank
                        targetRecvRank = targetRecvRank + Grid%Rj*Grid%NumRk
                    endif
                    if((Grid%Rj-jc >= 0 .and. Grid%Rj-jc < Grid%NumRj) .or. periodicJ) then
                        ! talking to a neighbor, or wrapping around periodic MPI boundaries
                        targetSendRank = targetSendRank + modulo(Grid%Rj-jc,Grid%NumRj)*Grid%NumRk
                    else
                        ! talking to a non-periodic boundary, so this is my own J rank
                        targetSendRank = targetSendRank + Grid%Rj*Grid%NumRk
                    endif

                    ! K rank offset
                    if((Grid%Rk+kc >= 0 .and. Grid%Rk+kc < Grid%NumRk) .or. periodicK) then
                        ! talking to a neighbor, or wrapping around periodic MPI boundaries
                        targetRecvRank = targetRecvRank + modulo(Grid%Rk+kc,Grid%NumRk)
                    else
                        ! talking to a non-periodic boundary, so this is my own K rank
                        targetRecvRank = targetRecvRank + Grid%Rk
                    endif
                    if((Grid%Rk-kc >= 0 .and. Grid%Rk-kc < Grid%NumRk) .or. periodicK) then
                        ! talking to a neighbor, or wrapping around periodic MPI boundaries
                        targetSendRank = targetSendRank + modulo(Grid%Rk-kc,Grid%NumRk)
                    else
                        ! talking to a non-periodic boundary, so this is my own K rank
                        targetSendRank = targetSendRank + Grid%Rk
                    endif

                    ! calculate the byte offset to the start of the data
                    ! remember send uses opposite order
                    recvDataOffset = 0
                    sendDataOffset = 0
                    SELECT CASE (ic)
                        case (-1)
                            ! receiving min i
                            ! no change to recvDataOffset
                            if(Grid%Ri+1 < Grid%NumRi .or. periodicI) then
                                ! sending to adjacent rank, send upper physical
                                sendDataOffset = sendDataOffset + Grid%Nip ! ghosts+physical-ghosts
                            else
                                ! no adjacent rank, send lower ghosts
                                ! no change to sendDataOffset
                            endif
                        case (0)
                            ! receiving central in i dimension
                            recvDataOffset = recvDataOffset + Model%nG
                            sendDataOffset = sendDataOffset + Model%nG
                        case (1)
                            ! receiving max i
                            recvDataOffset = recvDataOffset + Model%nG + Grid%Nip
                            if(Grid%Ri-1 >= 0 .or. periodicI) then
                                ! sending to adjacent rank, send lower physical
                                sendDataOffset = sendDataOffset + Model%nG
                            else
                                ! no adjacent rank, send upper ghosts
                                sendDataOffset = sendDataOffset + Model%nG + Grid%Nip
                            endif
                        CASE DEFAULT
                            print *, 'Unexpected ic direction'
                            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                    END SELECT

                    SELECT CASE (jc)
                        case (-1)
                            ! receiving min j
                            ! no change to recvDataOffset
                            if(Grid%Rj+1 < Grid%NumRj .or. periodicJ) then
                                ! sending to adjacent rank, send upper physical
                                sendDataOffset = sendDataOffset + Grid%Njp*Grid%Ni ! ghosts+physical-ghosts
                            else
                                ! no adjacent rank, send lower ghosts
                                ! no change to sendDataOffset
                            endif
                        case (0)
                            ! receiving central in j dimension
                            recvDataOffset = recvDataOffset + Model%nG*Grid%Ni
                            sendDataOffset = sendDataOffset + Model%nG*Grid%Ni
                        case (1)
                            ! receiving max j
                            recvDataOffset = recvDataOffset + (Model%nG + Grid%Njp)*Grid%Ni
                            if(Grid%Rj-1 >= 0 .or. periodicJ) then
                                ! sending to adjacent rank, send lower physical
                                sendDataOffset = sendDataOffset + Model%nG*Grid%Ni
                            else
                                ! no adjacent rank, send upper ghosts
                                sendDataOffset = sendDataOffset + (Model%nG + Grid%Njp)*Grid%Ni
                            endif
                        CASE DEFAULT
                            print *, 'Unexpected jc direction'
                            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                    END SELECT

                    SELECT CASE (kc)
                        case (-1)
                            ! receiving min k
                            ! no change to recvDataOffset
                            if(Grid%Rk+1 < Grid%NumRk .or. periodicK) then
                                ! sending to adjacent rank, send upper physical
                                sendDataOffset = sendDataOffset+Grid%Nkp*Grid%Ni*Grid%Nj ! ghosts+physical-ghosts
                            else
                                ! no adjacent rank, send lower ghosts
                                ! no change to sendDataOffset
                            endif
                        case (0)
                            ! receiving central in k dimension
                            recvDataOffset = recvDataOffset + Model%nG*Grid%Ni*Grid%Nj
                            sendDataOffset = sendDataOffset + Model%nG*Grid%Ni*Grid%Nj 
                        case (1)
                            ! receiving max k
                            recvDataOffset = recvDataOffset + (Model%nG + Grid%Nkp)*Grid%Ni*Grid%Nj
                            if(Grid%Rk-1 >= 0 .or. periodicK) then
                                ! sending to adjacent rank, send lower physical
                                sendDataOffset = sendDataOffset + Model%nG*Grid%Ni*Grid%Nj
                            else
                                ! no adjacent rank, send upper ghosts
                                sendDataOffset = sendDataOffset + (Model%nG + Grid%Nkp)*Grid%Ni*Grid%Nj
                            endif
                        CASE DEFAULT
                            print *, 'Unexpected kc direction'
                            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                    END SELECT

                    ! determine which of the previously created datatypes is the correct one
                    gasType = MPI_DATATYPE_NULL
                    SELECT CASE (abs(ic)+abs(jc)+abs(kc))
                        case (0) ! all physical cells. this is never actually copied
                            ! do nothing
                        case (1) ! face
                            if(ic /= 0) then
                                gasType = iFace5Gas
                            elseif(jc /= 0) then
                                gasType = jFace5Gas
                            else
                                gasType = kFace5Gas
                            endif
                        case (2) ! edge
                            if(ic == 0) then
                                gasType = iEdge5Gas
                            elseif(jc == 0) then
                                gasType = jEdge5Gas
                            else
                                gasType = kEdge5Gas
                            endif
                        case (3) ! corner
                            gasType = corner5Gas
                        CASE DEFAULT
                            print *, 'Sum of ic+jc+kc is nonsense'
                            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                    END SELECT

                    if(targetRecvRank /= rank) then
                        ! make sure I'm not talking to myself
                        localIndexIn = findloc(gamAppMpi%recvRanks, targetRecvRank, 1)
                        if(gamAppMpi%recvTypesGas(localIndexIn) == MPI_DATATYPE_NULL) then
                            ! not receiving any data from this rank yet, just add this datatype
                            call mpi_type_hindexed(1, (/ 1 /), recvDataOffset*dataSize, gasType, &
                                                   gamAppMpi%recvTypesGas(localIndexIn), ierr)
                        else
                            ! we're already receivng other data from this rank
                            !  merge the datatypes into a struct
                            ! need to use a temporary array so that the ints are of type MPI_ADDRESS_KIND
                            tempOffsets = (/ 0, recvDataOffset*dataSize /)
                            call mpi_type_create_struct(2, (/ 1, 1 /), tempOffsets, &
                                                        (/ gamAppMpi%recvTypesGas(localIndexIn), gasType /), &
                                                        gamAppMpi%recvTypesGas(localIndexIn),  ierr)
                        endif
                    endif
                    if(targetSendRank /= rank) then
                        ! make sure I'm not talking to myself
                        localIndexOut = findloc(gamAppMpi%sendRanks, targetSendRank, 1)
                        if(gamAppMpi%sendTypesGas(localIndexOut) == MPI_DATATYPE_NULL) then
                            ! not sending any data to this rank yet, just add this datatype
                            call mpi_type_hindexed(1, (/ 1 /), sendDataOffset*dataSize, gasType, &
                                                   gamAppMpi%sendTypesGas(localIndexOut), ierr)
                        else
                            ! we're already sending other data to this rank
                            !  merge the datatypes into a struct
                            ! need to use a temporary array so that the ints are of type MPI_ADDRESS_KIND
                            tempOffsets = (/ 0, sendDataOffset*dataSize /)
                            call mpi_type_create_struct(2, (/ 1, 1 /), tempOffsets, &
                                                        (/ gamAppMpi%sendTypesGas(localIndexOut), gasType /), &
                                                        gamAppMpi%sendTypesGas(localIndexOut), ierr)
                        endif
                    endif
                enddo
            enddo
        enddo

        ! commit the created MPI datatypes
        do localIndexOut=1,SIZE(gamAppMpi%sendTypesGas)
            call mpi_type_commit(gamAppMpi%sendTypesGas(localIndexOut), ierr)
        enddo
        do localIndexIn=1,SIZE(gamAppMpi%recvTypesGas)
            call mpi_type_commit(gamAppMpi%recvTypesGas(localIndexIn), ierr)
        enddo

        end associate

    end subroutine createCellCenteredDatatypes

    subroutine createFaceCenteredDatatypes(gamAppMpi, periodicI, periodicJ, periodicK)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi
        logical, intent(in) :: periodicI, periodicJ, periodicK

    end subroutine createFaceCenteredDatatypes

end module gamapp_mpi

