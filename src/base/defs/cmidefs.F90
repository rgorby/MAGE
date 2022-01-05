! Constants for CMI simulations

module cmidefs
    use kdefs

    implicit none

    integer, parameter :: JpSh = 1 !Number of cc current shells for CMI
    integer, parameter :: JpSt = 2 !First shell (i) to calculate current
    

    !Get 5 shells, all ghosts plus first active
    integer, parameter :: PsiSh = 5 !Number of *SHELLS* getting nodes at, ie PsiSh+1 = # i nodes
    integer, parameter :: PsiSt = -3 !Starting shell of potential

    real(rp), parameter :: defD0 = dalton !kg
    real(rp), parameter :: defV0 = 100.e3 !m/s
    real(rp), parameter :: defB0 = sqrt(Mu0*defD0)*defV0*1.0e+9 !T->nT
    real(rp), parameter :: defP0 = defD0*defV0*defV0*1.0e+9 !P->nPa
    

end module cmidefs