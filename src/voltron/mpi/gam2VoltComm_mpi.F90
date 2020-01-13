! Functions and data for MPI Gamera to talk to Voltron running on its own rank

module gam2VoltComm_mpi
    use gamapp_mpi

    implicit none

    type :: gam2VoltCommMpi_T
        integer :: voltMpiComm = MPI_COMM_NULL

        real(rp) :: time, tFin, DeepT, ShallowT
        integer :: ts
        logical :: doDeep

        type (IOClock_T) :: IO

    end type gam2VoltCommMpi_T

    contains

    ! setup the MPI communicator to talk to voltron, and send grid data
    subroutine initGam2Volt(g2vComm, gApp, voltComm)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(in) :: gApp
        integer, intent(in) :: voltComm

    end subroutine initGam2Volt

    ! transmit data to voltron so that it can keep up with the simulation
    subroutine performStepVoltron(g2vComm, gApp)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(in) :: gApp

    end subroutine performStepVoltron

    ! transmit state data to voltron over MPI and receive new data
    subroutine performShallowUpdate(g2vComm, gApp, time)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(inout) :: gApp
        real(rp), intent(in) :: time

    end subroutine performShallowUpdate

    ! transmit state data to voltron over MPI and receive new data
    subroutine performDeepUpdate(g2vComm, gApp, time)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(inout) :: gApp
        real(rp), intent(in) :: time

    end subroutine performDeepUpdate

end module gam2VoltComm_mpi

