! Data for converting between Gamera and Remix on a single OpenMP node

module mixmhdlocalinterface
    use mixapp
    use gamapp
    
    implicit none

    type mixMhdLocalInterface_T
    end type mixMhdLocalInterface_T

    contains

    subroutine convertGameraToRemix(gameraApp, remixApp)
        type(remixApp_T), intent(inout) :: remixApp
        type(gameraApp_T), intent(in) :: gameraApp

        ! convert incoming gamera data to the "remixInputs" variable
        real(rp) ::  B0mag,Bi2m
        integer :: i,j,k,iG
        real(rp) :: xc,yc,zc

        real(rp) :: Con(NVAR)
        real(rp) :: Cs


        !write(*,*) 'Prepping remix data at T = ',Model%t

        !Only working on Bxyz from perturbation
        !B0 in inner region (where we care for remix)
        !Assumd to be current-free
        call GetShellJ(gamGridXyz, gamGridXfc,gamGridVolume, gamStateBxyz, gJ)

        !Now loop over only shells and populate mhdvars to be sent to mix
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,iG,j,k,B0mag,Bi2m,xc,yc,zc,Con,Cs)
        do k=ks,ke
            do j=js,je
                do i=1,JShells
                    iG = JStart+i-1
                    B0mag = sqrt(dot_product(gamGridB0(iG,j,k,:),gamGridB0(iG,j,k,:)))
                    ! get cell centers
                    call cellCenter(gamGridXyz,iG,j,k,xc,yc,zc)
                    !Get ion2mag scaling factor for field
                    Bi2m = BIon2Mag(xc,yc,zc)

                    if (k<=ke/2) then
                        !!! NOTE: assuming gB0 in nT and gx0 in m and gv0 in m/s

                        ! note conversion to microA/m^2
                        remixInputs(i,j,k,MHDJ,NORTH) = dot_product(gJ(i,j,k,:),gamGridB0(iG,j,k,:)/B0mag)*Bi2m*(gB0/gx0*1.e4/4/PI)
                        ! note conversion to g/cm^3
                        remixInputs(i,j,k,MHDD,NORTH) = gamStateGas(iG,j,k,DEN,BLK)*(gB0*1.e-7/gv0)**2/4/pi
                        ! get sound speed first
                        Con = gamStateGas(iG,j,k,:,BLK)
                        call CellPress2Cs(gamma,Con,Cs)
                        remixInputs(i,j,k,MHDC,NORTH) = Cs*gv0*1.e2
                    else
                        remixInputs(i,j,k-ke/2,MHDJ,SOUTH) = dot_product(gJ(i,j,k,:),gamGridB0(iG,j,k,:)/B0mag)*Bi2m*(gB0/gx0*1.e4/4/PI)
                        remixInputs(i,j,k-ke/2,MHDD,SOUTH) = gamStateGas(iG,j,k,DEN,BLK)*(gB0*1.e-7/gv0)**2/4/pi
                        ! get sound speed first
                        Con = gamStateGas(iG,j,k,:,BLK)
                        call CellPress2Cs(gamma,Con,Cs)
                        remixInputs(i,j,k-ke/2,MHDC,SOUTH) = Cs*gv0*1.e2
                    endif
                enddo
            enddo
        enddo

    end subroutine convertGameraToRemix

    subroutine convertRemixToGamera(gameraApp, remixApp)
        type(remixApp_T), intent(inout) :: remixApp
        type(gameraApp_T), intent(in) :: gameraApp

        ! convert remix results to the "remixOutputs" variable
        integer :: i

         ! populate potential on gamera grid
         gPsiRemix = 0.0
         do i=1,PsiShells+1
            gPsiRemix(i,:,ks:ke/2+1)   = remixOutputs(i,:,:,MHDPSI,NORTH)
            gPsiRemix(i,:,ke/2+1:ke+1) = remixOutputs(i,:,:,MHDPSI,SOUTH)
         enddo

         ! add corotation
        call CorotationPot(gamGridXyz,gPsiRemix)

        call Ion2MHD(gamGridXyz,gamGridEdge,gPsiRemix,inEijkRemix,inExyzRemix,rm2g)

    end subroutine convertRemixToGamera

end module mixmhdlocalinterface

