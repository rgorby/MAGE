!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mhd_scalars
  !> Indices for status scalars to be sent from the MHD
  !! Note: These must match up with the scalars on the sending side.
!  USE Rcm_mod_subs, ONLY : rprec,iprec
  USE rcm_precision
  implicit none
  integer(iprec), parameter, public :: KILL_SIGNAL = 1 !> 111=kil; 222=last exchange; else=keep running
  integer(iprec), parameter, public :: YEAR = 2
  integer(iprec), parameter, public :: MONTH = 3
  integer(iprec), parameter, public :: DAY = 4
  integer(iprec), parameter, public :: HOUR = 5
  integer(iprec), parameter, public :: MINUTE = 6
  integer(iprec), parameter, public :: SECOND = 7
  integer(iprec), parameter, public :: DELTA_T = 8 !> # of seconds between exchange with coupled models
  integer(iprec), parameter, public :: LABEL = 9   !> MHD Simulation Time step
  integer(iprec), parameter, public :: NUMBER_OF_SCALARS = 9
  ! Scalars at index KILL_SIGNAL accept the following:
  integer(iprec), parameter, public :: KILL_SIGNAL_CONTINUE = 0
  integer(iprec), parameter, public :: KILL_SIGNAL_SHUTDOWN = 111
  integer(iprec), parameter, public :: KILL_SIGNAL_LAST_EXCHANGE = 222
   integer(iprec), public :: iaScalars(NUMBER_OF_SCALARS)      !> Status scalars sent from MHD

 end module mhd_scalars