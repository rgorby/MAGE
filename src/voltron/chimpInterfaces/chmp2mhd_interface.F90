! Converting chimp data to gamera data

module chmp2mhd_interface
    use gamtypes
    use math
    use gamapp
    use ebtypes
    use eqmap
    
    implicit none

    ! data for chimp -> gamera conversion
    type chmp2Mhd_T

        !Ni,Nj,Nk,2 array
        !Holds mapping from cell-centered xyz => x1,x1 (projection coordinates)
        !Projection coordinates can be R,phi (cylindrical) or lat,lon (ionospheric)
        real(rp), dimension(:,:,:,:), allocatable :: xyzSquish
        integer :: iMax !Possibly changing i-boundary of squish mapping

    end type chmp2Mhd_T

    contains

    subroutine init_chmp2Mhd(chmp2mhd, ebTrcApp, gamApp)
        type(chmp2Mhd_T), intent(inout) :: chmp2mhd
        type(ebTrcApp_T), intent(inout) :: ebTrcApp
        type(gamApp_T)  , intent(in)    :: gamApp

        associate(Gr=>gamApp%Grid)
        allocate(chmp2mhd%xyzSquish(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,2))
        chmp2mhd%xyzSquish = 0.0
        end associate

    end subroutine init_chmp2Mhd

    subroutine convertChimpToGamera(chmp2mhd,ebTrcApp,gamApp)
        type(chmp2Mhd_T), intent(inout) :: chmp2mhd
        type(ebTrcApp_T), intent(in)    :: ebTrcApp
        type(gamApp_T)  , intent(inout) :: gamApp

        integer :: i,j,k
        real(rp) :: D,P,x1,x2,t

        t = gamApp%Model%t*gamApp%Model%Units%gT0 !Get time scaled to seconds

        associate(Gr=>gamApp%Grid)

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,x1,x2,D,P)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is,Gr%is+chmp2mhd%iMax

                x1 = chmp2mhd%xyzSquish(i,j,k,1) 
                x2 = chmp2mhd%xyzSquish(i,j,k,2)
                call evalEQMap(x1,x2,t,D,P)

                Gr%Gas0(i,j,k,DEN     ) = D
                Gr%Gas0(i,j,k,PRESSURE) = P/gamApp%Model%Units%gP0 !Convert from nPa to code units
                Gr%Gas0(i,j,k,VELX:VELZ) = 0.0
                enddo
            enddo
        enddo

        end associate
        
    end subroutine convertChimpToGamera


    function lazyP(L,phi) result(P)
        real(rp), intent(in) :: L,phi
        real(rp) :: P

        real(rp) :: P0,L0,phi0,dL,dPhi
        real(rp) :: Lm,Pm

        P0 = 20.0 !nPa
        L0 = 5.0
        dL = 1.0
        phi0 = PI
        dPhi = 0.5*PI

        Lm = Mollify( abs(L-L0)    , dL  )
        Pm = Mollify( abs(phi-phi0), dPhi)
        P = P0*Lm*Pm

    end function lazyP

    !Wang++ 2013, empirical 2D pressure model
    function lesslazyP(L,phi) result(P)
        real(rp), intent(in) :: L,phi
        real(rp) :: P

        real(rp), dimension(9) :: B = [-0.7888,73.3651,38.6285,99.9520,-2.0651,69.7594,0.2883,9.8675,0.0787]
        real(rp) :: Sp,A1,A2

        Sp = sin(phi)
        A1 = B(2) + B(3)*Sp + B(4)*Sp*Sp
        A2 = B(6) + B(7)*Sp + B(8)*Sp*Sp

        P = exp(B(1)*L)*A1 + ( L**B(5) )*A2 + B(9)

    end function lesslazyP

end module chmp2mhd_interface
