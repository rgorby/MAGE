module rcm_precision
    use kdefs, ONLY: ip,rp,PI,strLen
  ! use gamera precision
    INTEGER, PARAMETER :: iprec = ip
    INTEGER, PARAMETER :: rprec = rp

    INTEGER(iprec), parameter :: RCMIOVARS = 50

    INTEGER(iprec), parameter :: ICONWRITERESTART = 31337
    INTEGER(iprec), parameter :: ICONWRITEOUTPUT = ICONWRITERESTART + 1
    INTEGER(iprec), parameter :: ICONRESTART     = ICONWRITERESTART - 1
end module rcm_precision
