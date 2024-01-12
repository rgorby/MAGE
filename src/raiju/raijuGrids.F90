

module raijugrids
    use ioh5
    use xml_input
    use shellgrid

    use raijudefs
    use raijutypes
    use raijuRecon
    use raijuSpeciesHelper

    implicit none

    contains

!------
! Spatial grid setup stuff
!------

    subroutine raijuGenUniSphGrid(Grid, iXML)
        type(raijuGrid_T)  , intent(inout) :: Grid
        type(XML_Input_T), intent(in)    :: iXML

        real(rp), dimension(:), allocatable :: Theta
        real(rp), dimension(:), allocatable :: Phi
        real(rp) :: dTheta, dPhi, thetaL, thetaU
        integer :: Nt,Np,Ng,Ngn,Ngs,Nge,Ngw
        integer :: i


        call iXML%Set_Val(thetaL , "grid/ThetaL", 15.0)
            !! Lower colat boundary [deg], ~15 Re in dipole
        call iXML%Set_Val(thetaU , "grid/ThetaU", 45.0)
            !! Upper colat boundary [deg], 2 Re in dipole
        call iXML%Set_Val(Nt, "grid/Nt", 31 )  ! 1 deg resolution
        call iXML%Set_Val(Np, "grid/Np", 361)  ! 1 deg resolution
        call iXML%Set_Val(Ng, "grid/Ng", 4  )  ! Number of ghosts, in every direction for now

        ! Turn degrees into radians
        thetaL = thetaL*deg2rad
        thetaU = thetaU*deg2rad

        ! Probably need to change
        Ngn = Ng
        Ngs = Ng
        Nge = Ng
        Ngw = Ng

        ! Allocate arrays
        allocate(Theta(Nt))
        allocate(Phi(Np))

        ! Create uniform grids
        dTheta = (thetaU-thetaL)/(Nt-1)
        dPhi = 2*PI/(Np-1)

        do i=1,Nt
            Theta(i) = thetaL + (i-1)*dTheta
        enddo

        do i=1,Np
            Phi(i) = (i-1)*dPhi
        enddo
        !write(*,*)"Orig theta:",Theta
        !write(*,*)"Orig phi:",Phi

        call GenShellGrid(Grid%shGrid,Theta,Phi,Ngn,Ngs,Nge,Ngw)

    end subroutine raijuGenUniSphGrid


    subroutine raijuGenGridFromShGrid(Grid, shGrid)
        type(raijuGrid_T)  , intent(inout) :: Grid
        type(ShellGrid_T), intent(in) :: shGrid

        !!TODO
        write(*,*) "You never should have come here"
        stop
    end subroutine raijuGenGridFromShGrid


    subroutine finalizeLLGrid(Grid, planet)
        !! Use a fully-created shell grid to allocate and populate the rest of the grid parameters
        !! NOTE: All lengths are either in units of meters or plentary radii, as specified below
        !!  no units of ionospheric radius (even though many of these parameters "live" on the ionospheric grid)
        type(raijuGrid_T), intent(inout) :: Grid
        type(planet_T), intent(in) :: planet

        integer :: i,j
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg) :: cosTh, BMagTh

        associate(shGr=>Grid%shGrid)
            ! First allocate remaining arrays

            allocate( Grid%delTh(shGr%isg:shGr%ieg+1) )
            allocate( Grid%delPh(shGr%jsg:shGr%jeg+1) )
            allocate( Grid%areaCC  (shGr%isg:shGr%ieg  ,shGr%jsg:shGr%jeg) )
            allocate( Grid%areaFace(shGr%isg:shGr%ieg+1,shGr%jsg:shGr%jeg+1, 2) )
            allocate( Grid%lenFace (shGr%isg:shGr%ieg+1,shGr%jsg:shGr%jeg+1, 2) )
            allocate( Grid%Bmag  (shGr%isg:shGr%ieg,shGr%jsg:shGr%jeg) )
            allocate( Grid%cosdip(shGr%isg:shGr%ieg,shGr%jsg:shGr%jeg) )
            allocate( Grid%Br(shGr%isg:shGr%ieg+1,shGr%jsg:shGr%jeg+1, 2) )
            

            ! Calculate theta/phi delta between cells [radians]
            ! First and last elements should be zero since there's no cells below/above isg/ieg
            Grid%delTh(shGr%isg) = 0.0
            Grid%delTh(shGr%ieg+1) = 0.0  
            do i=shGr%isg+1,shGr%ieg
                ! Define delta between th(i-1) and th(i)
                Grid%delTh(i) = abs(shGr%thc(i) - shGr%thc(i-1))
            enddo

            ! [radians]
            Grid%delPh(shGr%jsg) = 0.0
            Grid%delPh(shGr%jeg+1) = 0.0  
            do j=shGr%jsg+1,shGr%jeg
                ! Define delta between ph(i-1) and ph(i)
                Grid%delPh(j) = abs(shGr%phc(j) - shGr%phc(j-1))
            enddo

            ! Calc areas at cell centers [Rp^2]
            do i=shGr%isg,shGr%ieg
                do j=shGr%jsg,shGr%jeg
                    ! r^2 * sin(th) * dTh * dPh
                    Grid%areaCC(i,j) = (planet%ri_m/planet%rp_m)**2 &
                                        * sin(shGr%thc(i)) &
                                        * (shGr%th(i+1) - shGr%th(i)) &
                                        * (shGr%ph(j+1) - shGr%ph(j))
                enddo
            enddo

            ! Area at faces [Rp^2]
            ! Kinda overkill, but just in case: use 8-centered reconstruction for each face
            ! Only do for faces of active cells, because we don't need them in ghosts
            Grid%areaFace = 0.0
            do i=shGr%is,shGr%ie+1
                do j=shGr%js,shGr%je+1
                    ! Theta dir
                    Grid%areaFace(i,j,RAI_TH) = Central8(Grid%areaCC(i-4:i+3, j))
                    ! Phi dir
                    Grid%areaFace(i,j,RAI_PH) = Central8(Grid%areaCC(i, j-4:j+3))
                enddo
            enddo

            ! Arc length of faces [Rp]
            do i=shGr%isg,shGr%ieg
                do j=shGr%jsg,shGr%jeg
                    ! Faces in theta dir are i +/- 1/2 faces, each of which spans the phi direction
                    
                    ! Theta dir
                    Grid%lenFace(i,j,RAI_TH) = (planet%ri_m/planet%rp_m) &
                                        * sin(shGr%th(i)) &
                                        * (shGr%ph(j+1) - shGr%ph(j))  ! Line segment in phi direction
                    ! Phi dir
                    Grid%lenFace(i,j,RAI_PH) = (planet%ri_m/planet%rp_m) &
                                        * (shGr%th(i+1) - shGr%th(i))  ! Line segment in theta direction
                enddo
            enddo
            ! Do last row's top edges too
            i = shGr%ieg+1
            do j=shGr%jsg,shGr%jeg
                Grid%lenFace(i, j, RAI_TH) = (planet%ri_m/planet%rp_m) &
                                            * sin(shGr%th(shGr%ie+1)) &
                                            * (shGr%ph(shGr%jeg+1) - shGr%ph(shGr%jeg))
            enddo
            ! And last column's side edges
            j = shGr%jeg+1
            do i=shGr%isg,shGr%ieg
                Grid%lenFace(i,j, RAI_PH) = (planet%ri_m/planet%rp_m) &
                                            * (shGr%th(shGr%ieg+1) - shGr%th(shGr%ieg))
            enddo


            ! Calc ionospheric Bmag [nT] and cos of dip angle at cell centers
            cosTh = cos(shGr%thc)
            BMagTh = planet%magMoment*G2nT &
                    /(planet%ri_m/planet%rp_m)**3.0 &
                    * sqrt(1.0+3.0*cosTh**2.0)  ! [nT]
            do j=shGr%jsg,shGr%jeg
                Grid%Bmag(:,j) = BMagTh
                Grid%cosdip(:,j) = 2.0*cosTh/sqrt(1.0 + 3.0*cosTh**2.0) 
            enddo

            
            ! Radial component of Bmag [nT] at cell edges, used for velocity calculation
            Grid%Br = 0.0
            ! Theta-dir faces have theta locations at cell interfaces
            cosTh = cos(shGr%th)  ! Ni+1
            do j=shGr%jsg,shGr%jeg
                Grid%Br(:,j,RAI_TH) = planet%magMoment*G2nT &
                                    /(planet%ri_m/planet%rp_m)**3.0 &
                                    * 2*cosTh  ! [nT]
            enddo
            ! Phi-dir faces have theta locations at cell centers
            cosTh = cos(shGr%thc)  ! Ni
            do j=shGr%jsg,shGr%jeg+1
                Grid%Br(:,j,RAI_PH) = planet%magMoment*G2nT &
                                    /(planet%ri_m/planet%rp_m)**3.0 &
                                    * 2*cosTh  ! [nT]
            enddo

        end associate

    end subroutine finalizeLLGrid


!------
! Spatial grid operations
!------
    subroutine wrapJcc(sh, Q)
        !! Populates overlapping j ghosts with active cell-centered data
        !! Assumes that all values within (:, js:je) are valid to wrap
        !! And assumes that j ghosts overlap with last active cells on the other side
        type(ShellGrid_T), intent(in) :: sh
        real(rp), dimension(sh%isg:sh%ieg, sh%jsg:sh%jeg), intent(inout) :: Q

        ! Starting ghost cells
        Q(:, sh%jsg:sh%js-1) = Q(:, sh%je-sh%Ngw+1:sh%je)
        ! Ending ghosts cells
        Q(:, sh%je+1:sh%jeg) = Q(:, sh%js:sh%js+sh%Nge-1)

    end subroutine wrapJcc


    subroutine wrapJcorners(sh, Q)
        !! Populates overlapping j ghosts with active corner data
        !! Assumes that all values within (:, js:je) are valid to wrap
        !! And assumes that j ghosts overlap with last active cells on the other side
        type(ShellGrid_T), intent(in) :: sh
        real(rp), dimension(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1), intent(inout) :: Q

        ! Starting ghost cells
        Q(:, sh%jsg:sh%js) = Q(:, sh%je-sh%Ngw+1:sh%je+1)
        ! Ending ghosts cells
        Q(:, sh%je+1:sh%jeg+1) = Q(:, sh%js:sh%js+sh%Nge)

    end subroutine wrapJcorners


!------
! Lambda grid stuff
!------

!! Maybe should leave just spatial grid stuff in raijuGrids and move lambda stuff to a lambdaUtils

    subroutine initLambdaGrid(Model, Grid, configfname)
        type(raijuModel_T)  , intent(inout) :: Model
        type(raijuGrid_T), intent(inout) :: Grid
        character(len=strLen), intent(in) :: configfname

        ! First read in species information from config file
        call populateSpeciesFromConfig(Model, Grid, configfname)

        ! Now prepare our alamc grid, tell each Species what its range is in alamc array
        call initAlamc(Grid)
        
    end subroutine initLambdaGrid

end module raijugrids