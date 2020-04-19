! Constants for CMI simulations

module cmidefs
    use kdefs

    implicit none

    integer, parameter :: JpSh = 1 !Number of cc current shells for CMI
    integer, parameter :: JpSt = 2 !First shell (i) to calculate current
    

    !Get 5 shells, all ghosts plus first active
    integer, parameter :: PsiSh = 5 !Number of *SHELLS* getting nodes at, ie PsiSh+1 = # i nodes
    integer, parameter :: PsiSt = -3 !Starting shell of potential

end module cmidefs