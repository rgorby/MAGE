module rcm_precision
    use kdefs, ONLY: ip,rp,PI,strLen
  ! use gamera precision
    INTEGER, PARAMETER :: iprec = ip
    INTEGER, PARAMETER :: rprec = rp

    INTEGER(iprec), parameter :: RCMIOVARS = 50

    INTEGER(iprec), parameter :: ICONWRITERESTART = 31337
    INTEGER(iprec), parameter :: ICONWRITEOUTPUT = ICONWRITERESTART + 1
    INTEGER(iprec), parameter :: ICONRESTART     = ICONWRITERESTART - 1
    INTEGER(iprec), parameter :: RCMELECTRON = 1
    INTEGER(iprec), parameter :: RCMPROTON   = 2
    INTEGER(iprec), parameter :: RCMCLAWLIM = 1
	!Limiters: 1=minmod, 2=superbee, 3=Van Leer, 4=MC, 5=Beam Warming
	!           6=Fromm, 7=Albada2, 8=Albada3, 9=Ultrabee, 10=UltrabeeMod, 11=Arora Roe

end module rcm_precision
