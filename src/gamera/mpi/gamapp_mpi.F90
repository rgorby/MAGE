! Main data objects and functions to perform a gamera simulation

module gamapp_mpi
    use gamtypes
    use step
    use init
    use mhdgroup
    use gamapp
    use mpi

    implicit none

    type, extends(GamApp_T) :: gamAppMpi_T
        integer :: NumRi=1,NumRj=1,NumRk=1
        integer :: Ri=1,Rj=1,Rk=1
        type(MPI_Comm) :: gamMpiComm
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

        integer :: numNeighbors, ierr, length, commSize, rank, ic, jc, kc, rankOffset
        integer, dimension(26), allocatable :: sourceRanks, sourceData
        logical :: reorder,periodicI,periodicJ,periodicK
        character(len=strLen) :: message
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

        gamAppMpi%Grid%ijkShift(1) = gamAppMpi%Grid%Nip*gamAppMpi%Ri
        gamAppMpi%Grid%ijkShift(2) = gamAppMpi%Grid%Njp*gamAppMpi%Rj
        gamAppMpi%Grid%ijkShift(3) = gamAppMpi%Grid%Nkp*gamAppMpi%Rk

        gamAppMpi%Grid%Ni = gamAppMpi%Grid%Nip + 2*gamAppMpi%Model%nG
        gamAppMpi%Grid%Nj = gamAppMpi%Grid%Njp + 2*gamAppMpi%Model%nG
        gamAppMpi%Grid%Nk = gamAppMpi%Grid%Nkp + 2*gamAppMpi%Model%nG

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

        ! move the new arryas to the grid corner arrays
        call move_alloc(tempX, gamAppMpi%Grid%x)
        call move_alloc(tempY, gamAppMpi%Grid%y)
        call move_alloc(tempZ, gamAppMpi%Grid%z)

        ! check which dimensions are using MPI for periodicity
        SELECT type(iiBC=>gamAppMpi%Grid%externalBCs(INI)%p)
            TYPE IS (periodicInnerIBCMpi_T)
                periodicI = .true.
            CLASS DEFAULT
                periodicI = .false.
        END SELECT

        SELECT type(iiBC=>gamAppMpi%Grid%externalBCs(INJ)%p)
            TYPE IS (periodicInnerJBCMpi_T)
                periodicJ = .true.
            CLASS DEFAULT
                periodicJ = .false.
        END SELECT

        SELECT type(iiBC=>gamAppMpi%Grid%externalBCs(INK)%p)
            TYPE IS (periodicInnerKBCMpi_T)
                periodicK = .true.
            CLASS DEFAULT
                periodicK = .false.
        END SELECT

        ! create MPI topology
        numNeighbors = 0
        do ic=-1,1
            if((gamAppMpi%Ri+ic >= 0 .and. gamAppMpi%Ri+ic < GamAppMpi%NumRi) .or. periodicI) then
                do jc=-1,1
                    if((gamAppMpi%Rj+jc >= 0 .and. gamAppMpi%Rj+jc < GamAppMpi%NumRj) .or. periodicJ) then
                        do kc=-1,1
                            if((gamAppMpi%Rk+kc >= 0 .and. gamAppMpi%Rk+kc < GamAppMpi%NumRk) .or. periodicK)then
                                rankOffset = abs(ic)+abs(jc)+abs(kc)
                                targetRank = ...
                                listIndex = findloc(sourceRanks,targetRank)
                                
                                if(rankOffset == 0) then ! don't talk to yourself
                                elseif (rankOffset == 3) then ! corner
                                    numNeighbors=numNeighbors+1
                                    sourceRanks(numNeighbors) = ...
                                    sourceData(numNeibhors) = ...
                                elseif (rankOffset == 2) then ! edge
                                    numNeighbors=numNeighbors+1
                                else ! face
                                    numNeighbors=numNeighbors+1
                                endif
                            endif
                        enddo
                    endif
                enddo
            endif
        enddo

        reorder = .false. ! DO NOT allow MPI to reorder the ranks
        call mpi_dist_graph_create_adjacent(gamComm,
            numNeighbors,sourceRanks,sourceData,
            numNeighbors,sourceRanks,sourceData, ! comms symmetrical
            mpiInfo, reorder, gamAppMpi%gamMpiComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

         ! call appropriate subroutines to calculate all appropriate grid data from the corner data
        call CalcGridInfo(gamAppMpi%Model,gamAppMpi%Grid,gamAppMpi%State,gamAppMpi%oState,gamAppMpi%Solver,xmlInp,userInitFunc)

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
    end subroutine haloUpdate

end module gamapp_mpi

