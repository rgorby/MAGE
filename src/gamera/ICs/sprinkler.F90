module usergamic
    use gamtypes
    use gamutils
    use math
    use gridutils
    use xml_input
    use bcs
    use background

    implicit none

    !Various global would go here
    real(rp) :: Rho0, P0, V0, wScl, tOrb, Cs0, B0
    contains

    subroutine initUser(Model,Grid,State,inpXML)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML
        procedure(GasIC_T), pointer :: Wxyz

        real(rp) :: vMach
        
        !Get defaults from input deck
        call inpXML%Set_Val(Cs0  ,"prob/Cs0"  ,1.0_rp)
        call inpXML%Set_Val(vMach,"prob/vMach",10.0_rp)
        call inpXML%Set_Val(wScl ,"prob/wScl" ,PI/8.0)
        call inpXML%Set_Val(tOrb ,"prob/tOrb" ,10.0_rp)
        call inpXML%Set_Val(B0   ,"prob/B0"   ,1.0_rp)


        Rho0 = 1.0_rp
        P0 = Rho0*Cs0**2.0/Model%gamma
        Cs0 = sqrt(Model%gamma*P0/Rho0)

        V0 = vMach*Cs0

        !Set BCs (spherical, RPT)
        !Rearrange order to do periodic first
        Grid%HaloUps(5)%ApplyBC => periodic_ibcK
        Grid%HaloUps(6)%ApplyBC => periodic_obcK

        Grid%HaloUps(3)%ApplyBC => zeroGrad_ibcJ 
        Grid%HaloUps(4)%ApplyBC => zeroGrad_obcJ

        Grid%HaloUps(1)%ApplyBC => SprinklerBC
        Grid%HaloUps(2)%ApplyBC => zeroGrad_obcI


        !Map IC to grid
        Wxyz => GasIC
        call GasIC2State(Model,Grid,State,Wxyz)

        !Local functions
        !NOTE: Don't put BCs here as they won't be visible after the initialization call

        contains
            subroutine GasIC(x,y,z,D,Vx,Vy,Vz,P)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: D,Vx,Vy,Vz,P

                D = Rho0
                P = P0
                Vx = 0.0
                Vy = 0.0
                Vz = 0.0
            end subroutine GasIC

    end subroutine initUser

    !Put BCs here for global access
    !Inner-I BC for rotating inflow region
    subroutine SprinklerBC(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        integer :: n,ig,j,k
    
        real(rp) :: rHat(NDIM), n0(NDIM)
        real(rp) :: D,P,Vx,Vy,Vz, KinE,TotE
        real(rp) :: R,Phi,Theta,xc,yc,zc
        real(rp) :: vMag,dSig, Phi0
    
        !Find normal vector of center of sprinkler
        Phi0 = PI/2.0 + 2*PI*State%time/tOrb
        n0 = [cos(Phi0),sin(Phi0),0.0_rp] 

        !i-boundaries (IN)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg

                !Save D/P
                D = Rho0
                P = P0
                do n=1,Model%Ng
                    !Find cell center in this ghost
                    ig = Grid%is-n
                    
                    call cellCenter(Grid,ig,j,k,xc,yc,zc)

                    !Calculate radial hat vector
                    R = sqrt(xc**2.0 + yc**2.0 + zc**2.0)
                    Theta = acos(zc/R)
                    Phi = atan2(yc,xc)
    
                    rHat = [xc,yc,zc]/R
                    !Find angular distance between this outward normal and n0
                    dSig = abs(acos(dot_product(n0,rHat)))

                    vMag = v0*Mollify(dSig,wScl)
                    !vMag = v0

                    !Set primitives (already have D/P)
                    ! Vx = n0(XDIR)*vMag
                    ! Vy = n0(YDIR)*vMag
                    ! Vz = n0(ZDIR)*vMag
                    Vx = rHat(XDIR)*vMag
                    Vy = rHat(YDIR)*vMag
                    Vz = rHat(ZDIR)*vMag

    
                    KinE = 0.5*D*(Vx**2.0+Vy**2.0+Vz**2.0)
                    TotE = KinE + P/(Model%gamma-1)
    
                    State%Gas(ig,j,k,DEN,BLK) = D
                    State%Gas(ig,j,k,MOMX,BLK) = D*Vx
                    State%Gas(ig,j,k,MOMY,BLK) = D*Vy
                    State%Gas(ig,j,k,MOMZ,BLK) = D*Vz
                    State%Gas(ig,j,k,ENERGY,BLK) = TotE
                    
                    State%magFlux(ig,j,k,:) = 0.0
                    State%magFlux(ig,j,k,KDIR) = B0*Mollify(abs(PI/2 - Theta),wScl)*Grid%face(ig,j,k,KDIR)
                    
                enddo
            enddo
        enddo
    
    end subroutine SprinklerBC

end module usergamic
