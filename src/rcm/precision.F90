module rcm_precision
    use kdefs, ONLY: ip,rp,PI,strLen,Me_cgs,Mp_cgs
  ! use gamera precision
    INTEGER, PARAMETER :: iprec = ip
    INTEGER, PARAMETER :: rprec = rp

    INTEGER(iprec), parameter :: RCMIOVARS = 50

    INTEGER(iprec), parameter :: ICONWRITERESTART = 31337
    INTEGER(iprec), parameter :: ICONWRITEOUTPUT  = ICONWRITERESTART + 1
    INTEGER(iprec), parameter :: ICONRESTART      = ICONWRITERESTART - 1
    INTEGER(iprec), parameter :: RCMELECTRON = 1
    INTEGER(iprec), parameter :: RCMPROTON   = 2
    INTEGER(iprec), parameter :: RCMOXYGEN   = 3
    INTEGER(iprec), parameter :: RCMNUMFLAV = 2 !Number of RCM flavors
end module rcm_precision
