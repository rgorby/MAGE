module rcmdefs
    use kdefs, ONLY: ip,rp,PI,strLen

    implicit none

    ! use gamera precision
    INTEGER, PARAMETER :: iprec = ip
    INTEGER, PARAMETER :: rprec = rp

    INTEGER(iprec), parameter :: RCMIOVARS = 50

    INTEGER(iprec), parameter :: ICONWRITERESTART = 31337
    INTEGER(iprec), parameter :: ICONWRITEOUTPUT = ICONWRITERESTART + 1
    INTEGER(iprec), parameter :: ICONRESTART     = ICONWRITERESTART - 1

    integer(ip), parameter :: RCMINIT=0,RCMADVANCE=1,RCMRESTART=2,RCMWRITERESTART=-2,RCMWRITEOUTPUT=-3,RCMWRITETIMING=-1
    logical :: doRCMVerbose = .false.
    integer(ip), parameter :: RCMTOPCLOSED=-1,RCMTOPOPEN=+1

end module rcmdefs

