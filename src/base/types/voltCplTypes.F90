module voltCplTypes
    !! Lower-level types useful for voltron coupling
    !! Putting here to remove the basic types from the voltron-dependent routines that use then
    !! So other modules can rely on the types without caring about how they're used by others

    use kdefs
    use shellGrid
    use ebtypes

    implicit none

    integer, parameter :: MAXTUBEFLUIDS = 5
    type IMAGTube_T
        real(rp) :: Vol,bmin,beta_average
            !! Flux tube volume, minimum b along FL
            !! Average plasma beta
        real(rp), dimension(0:MAXTUBEFLUIDS) :: Pave, Nave
            !! average pressure, average density
        real(rp), dimension(0:MAXTUBEFLUIDS) :: Pstd, Nstd
            !! standard deviation of pressure and densiy along the field line
        real(rp) :: X_bmin(NDIM)
        integer(ip) :: topo
        real(rp) :: latc,lonc !Conjugate lat/lon
        real(rp) :: Lb, Tb !Arc length/bounce time
        real(rp) :: losscone,rCurv,wIMAG,TioTe0=4.0_rp
        real(rp) :: Veb  ! [km/s]
    end type IMAGTube_T


    type IMAGTubeShell_T
        type(ShellGrid_T) :: sh
            !! Home grid

        type(magLine_T), dimension(:,:), allocatable :: bTrc2D

        type(ShellGridVar_T) :: vol
            !! Flux tube volume
        type(ShellGridVar_T) :: bmin
            !! Magnitude B at bmin location
        type(ShellGridVar_T) :: beta_ave
            !! Average plasma beta along field line
        type(ShellGridVar_T), dimension(0:MAXTUBEFLUIDS) :: Pave
            !! (N_MHDspc) average plasma pressure [nPa]
        type(ShellGridVar_T), dimension(0:MAXTUBEFLUIDS) :: Nave
            !! (N_MHDspc) average plasma density [#/cc]
        type(ShellGridVar_T), dimension(0:MAXTUBEFLUIDS) :: Pstd
            !! (N_MHDspc) Standard deviation of pressure along field line
        type(ShellGridVar_T), dimension(0:MAXTUBEFLUIDS) :: Nstd
            !! (N_MHDspc) Standard deviation of density along field line
        type(ShellGridVar_T), dimension(NDIM) :: X_bmin
            !! (NDIM) xyz-SM location of bmin surface
        type(ShellGridVar_T) :: latc, lonc
            !! Conjugate lat/lon
        type(ShellGridVar_T) :: Lb, Tb
            !! Arc length, bounce time
        type(ShellGridVar_T) :: losscone, rCurv, wIMAG, TioTe
        type(ShellGridVar_T), dimension(2) :: Veb  ! (Theta, phi dir ExB velocity)
    end type IMAGTubeShell_T

end module voltCplTypes