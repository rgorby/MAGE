module voltCplTypes
    !! Lower-level types useful for voltron coupling
    !! Putting here to remove the basic types from the voltron-dependent routines that use then
    !! So other modules can rely on the types without caring about how they're used by others

    use kdefs
    use shellGrid
    use ebtypes

    implicit none

    integer, parameter :: MAXTUBEFLUIDS = 5

    enum, bind(C)
        enumerator :: TUBE_CLOSED=0, TUBE_OPEN, TUBE_UNDEF, TUBE_DISASTER
    endenum

    !New tube type for global Voltron grid
    type Tube_T
        real(rp) :: xyz0(NDIM)
        !! XYZ of seed point [Rx]
        real(rp) :: lat0,lon0,invlat
        !! Lat/lon of tube footpoint and invariant latitude [rad]
        real(rp) :: latc,lonc
        !! Lat/lon of conjugate point [rad] if topo=TUBE_CLOSED
        integer :: topo
        !! topo = TUBE_CLOSED,OPEN,UNDEF
        real(rp) :: bmin, X_bmin(NDIM)
        !! Minimum B along field line [nT] and coordinates [Rx]
        real(rp) :: bVol
        !! Flux-tube volume [Rx/nT]
        real(rp) :: Lb
        !! Length of field line [Rx]
        real(rp) :: Tb
        !! Alfven bounce time [s] if CLOSED
        real(rp) :: pot,crpot,potc
        !! Electrostatic potential (TOTAL and corotation only), and conjugate potential [kV]
        real(rp) :: wMAG
        !! Energy partition, Mag E / (Kin + Mag + Thermal). 0 <= wMAG <= 1
        !! wMAG = 1 => magnetically dominated, wMAG << 1 => not strongly magnetized
        real(rp) :: rCurv
        !!Curvature radius [Rx] @ bmin if CLOSED
        real(rp) :: avgBeta
        !! Flux-tube volume-weighted average beta
        real(rp), dimension(0:MAXTUBEFLUIDS) :: avgP,avgN,stdP,stdN
        !! Average and standard deviation of pressure [nPa] and number density [#/cc] on field line
        real(rp) :: losscone,lossconec
        !! Size of losscone [RAD] at footpoint and conjugate (if CLOSED)
        real(rp) :: TioTe0
        !! Empirical Ti/Te estimate if it exists
        integer :: nTrc
        !! Number of points traced for this tube
    end type Tube_T

    type TubeShell_T
        type(ShellGridVar_T), dimension(NDIM) :: xyz0
            !! XYZ of seed point [Rx]
        type(ShellGridVar_T) :: lat0, lon0, invlat
            !! Lat/lon of tube footpoint and invariant latitude [rad]
        type(ShellGridVar_T) :: latc, lonc
            !! Conjugate lat/lon
        type(ShellGridVar_T) :: topo
            !! topo = TUBE_CLOSED,OPEN,UNDEF
            !! Note: integer data cast to real precision
        type(ShellGridVar_T) :: bmin
            !! Magnitude B at bmin location
        type(ShellGridVar_T), dimension(NDIM) :: X_bmin
            !! (NDIM) xyz-SM location of bmin surface
        type(ShellGridVar_T) :: bVol
            !! Flux tube volume [Rx/nT]
        type(ShellGridVar_T) :: Lb
            !! Length of field line [Rx]
        type(ShellGridVar_T) :: Tb
            !! Alfven bounce time [s] if CLOSED
        type(ShellGridVar_T) :: wMAG
            !! Energy partition, Mag E / (Kin + Mag + Thermal). 0 <= wMAG <= 1
            !! wMAG = 1 => magnetically dominated, wMAG << 1 => not strongly magnetized
        type(ShellGridVar_T) :: rCurv
            !! Curvature radius [Rx] @ bmin if CLOSED
        type(ShellGridVar_T) :: avgBeta
            !! Flux-tube volume-weighted average beta
        type(ShellGridVar_T), dimension(0:MAXTUBEFLUIDS) :: avgP
            !! (N_MHDspc) average plasma pressure [nPa]
        type(ShellGridVar_T), dimension(0:MAXTUBEFLUIDS) :: avgN
            !! (N_MHDspc) average plasma density [#/cc]
        type(ShellGridVar_T), dimension(0:MAXTUBEFLUIDS) :: stdP
            !! (N_MHDspc) Standard deviation of pressure along field line
        type(ShellGridVar_T), dimension(0:MAXTUBEFLUIDS) :: stdN
            !! (N_MHDspc) Standard deviation of density along field line
        type(ShellGridVar_T) :: losscone
            !! Size of losscone [RAD] at footpoint (if CLOSED)
        type(ShellGridVar_T) :: lossconec
            !! Size of losscone [RAD] at conjugate footpoint (if CLOSED)
        type(ShellGridVar_T) :: TioTe0
            !! Empirical Ti/Te estimate if it exists
    end type TubeShell_T

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

end module voltCplTypes