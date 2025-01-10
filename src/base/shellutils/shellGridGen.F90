!> Helper module for generating theta and phi arrays to use within a ShellGrid
module shellGridGen
    use math
    use kdefs

    implicit none

    contains

    subroutine genThetaPhi_uniform(Nt, Np, thetaL, thetaU, theta, phi)
        integer, intent(in) :: Nt
            !! Number of active cells in theta direction
        integer, intent(in) :: Np
            !! Number of active cells in phi direction
        real(rp), intent(in) :: thetaL
            !! [deg] Lower theta bound (0 = pole, 90 = equator)
        real(rp), intent(in) :: thetaU
            !! [deg] Upper theta bound (0 = pole, 90 = equator)
        real(rp), dimension(Nt+1), intent(inout) :: theta
            !! Populated active theta grid we return
        real(rp), dimension(Np+1), intent(inout) :: phi
            !! Populated active phi grid we return

        real(rp) :: thetaL_rad, thetaU_rad
        real(rp) :: dTheta, dPhi
        integer :: i

        ! Get theta bounds in radians
        thetaL_rad = thetaL*deg2rad
        thetaU_rad = thetaU*deg2rad

        dTheta = (thetaU_rad-thetaL_rad)/Nt
        dPhi = 2.0*PI/Np

        do i=1,Nt+1
            theta(i) = thetaL_rad + (i-1)*dTheta
        enddo

        do i=1,Np+1
            phi(i) = (i-1)*dPhi
        enddo
        ! Catch for slight overshoot
        if ((phi(Np+1) > 2*PI) .and. (phi(Np+1) - 2*PI < TINY) ) then
            phi(Np+1) = 2.0*PI
        endif

    end subroutine genThetaPhi_Uniform


end module shellGridGen