! Constants for CMI simulations

module cmidefs
    use kdefs

    implicit none

    integer, parameter :: JpSh = 1 !Number of cc current shells for CMI
    integer, parameter :: JpSt = 2 !First shell (i) to calculate current
    

    !Get 5 shells, all ghosts plus first active
    integer, parameter :: PsiSh = 5 !Number of *SHELLS* getting nodes at, ie PsiSh+1 = # i nodes
    integer, parameter :: PsiSt = -3 !Starting shell of potential

    !Default spinup time [s]
    real(rp), parameter :: tSpinDef = 14400.0
    
    real(rp), parameter :: defD0 = dalton*1.0D+6 !AMU/cc -> kg/m3
    real(rp), parameter :: defV0 = 1.0D+5 ! 100 km/s in m/s
    real(rp), parameter :: defB0 = sqrt(Mu0*defD0)*defV0*1.0D+9 !T->nT
    real(rp), parameter :: defP0 = defD0*defV0*defV0*1.0D+9 !P->nPa

end module cmidefs