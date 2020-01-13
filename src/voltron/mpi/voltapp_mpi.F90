! Collection of data and objects for the voltron middle man
! MPI version

module voltapp_mpi
    use voltapp
    use gamapp_mpi
    use gamapp
    use mpi
    
    implicit none

    type, extends(voltApp_T) :: voltAppMpi_T
        integer :: voltMpiComm = MPI_COMM_NULL
        type(gamApp_T) :: gAppLocal
    end type voltAppMpi_T

    contains

    !Initialize Voltron (after Gamera has already been initialized)
    subroutine initVoltron_mpi(vApp, voltComm, optFilename)
        type(voltAppMpi_T), intent(inout) :: vApp
        integer, intent(in) :: voltComm
        character(len=*), optional, intent(in) :: optFilename
        integer :: commSize, ierr

        ! create a new communicator using MPI Topology

        ! create a stripped down gamApp object which will be used locally
        ! voltron uses Grid%is,js,ks,ie,je,ke - Grid%xyz - Grid%B0
        

        ! use standard voltron with local gamApp object
        if(present(optFilename)) then
            call initVoltron(vApp, vApp%gAppLocal, optFilename)
        else
            call initVoltron(vApp, vApp%gAppLocal)
        endif

    end subroutine initVoltron_mpi

    ! MPI version of updating voltron variables
    subroutine stepVoltron_mpi(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

    end subroutine stepVoltron_mpi

!----------
!Shallow coupling stuff
    subroutine ShallowUpdate_mpi(vApp, time)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp) :: time

        ! fetch data from Gamera ranks
        ! remix uses State%Bxyz
        ! remix uses State%Gas(:,:,:,:,BLK)

        ! call base update function with local data
        call ShallowUpdate(vApp, vApp%gAppLocal, time)

        ! send updated data to Gamera ranks
        ! remix updates inEijk and inExyz in the IonInnerBC_T

    end subroutine ShallowUpdate_mpi

!----------
!Deep coupling stuff (time coming from vApp%time, so in seconds)
    subroutine DeepUpdate_mpi(vApp, time)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), intent(in) :: time

        real(rp) :: tAdv

        if (.not. vApp%doDeep) then
            !Why are you even here?
            return
        endif

        ! fetch data from Gamera ranks
        ! chimp uses State%Gas(:,:,:,:,BLK)
        ! chimp uses State%Bxyz(:,:,:,:)

        ! call base update function with local data
        call DeepUpdate(vApp, vApp%gAppLocal, time)

        ! send updated data to Gamera ranks

    end subroutine DeepUpdate_mpi

end module voltapp_mpi

