
!! Gibson-Low defined constants
module gldefs
    use kdefs
    implicit none
    ! Variable labels
    enum, bind(C)
        ! Density, Vr, Vtheta, Vphi, Pressure, Br, Btheta, Bphi, inside_mask
        ! note, theta and phi are polar and azimuthal angles respectively
        ! phi=atan2(y,x)
        enumerator :: GLDEN = 1, GLVR, GLVT, GLVP, GLPRESSURE, GLBR, GLBT, GLBP, GLMASK
    end enum

    enum, bind(C)
        enumerator :: SPHERICAL = 1, CARTESIAN
    end enum

    real(rp), parameter :: alnotrbub0 = 5.763854
    real(rp), dimension(6), parameter :: alnotrbubuse1 = &
                                            (/5.763459, 9.095011, 12.322941, 15.514603, 18.689036, 21.853874/)
    real(rp), parameter :: &
        mdtor = pi/180, & ! TODO: remove, already in kdefs
        Rsun = 6.96e10, & ! TODO: remove, already in kdefs as KM
        GMm = 221.3, &
        mprot = 1.674e-24, & !TODO: remove, already in kdefs Mp_cgs
        mp = 1.6696e-24, & !TODO: where is this from? different from kdefs
        gas_R = 1./6.07e-9, &
        kboltz = 1.3807e-16, & !TODO: remove, already in kdefs, Kbltz
        BIGTINY = 1e-4, &
        eta0 = 1.38d-7

end module