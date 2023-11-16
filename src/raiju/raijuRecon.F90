module raijuRecon
    !! Reconstruction routines
    use raijudefs
    use raijuTypes

    implicit none

    real(rp), dimension(recLen), parameter :: interpWgt_4c = [ 0, 0,  -1,  7,  7,  -1, 0, 0]/60.0_rp
    real(rp), dimension(recLen), parameter :: interpWgt_6c = [ 0, 1,  -8, 37, 37,  -8, 1, 0]/60.0_rp
    real(rp), dimension(recLen), parameter :: interpWgt_8c = [-3,29,-139,533,533,-139,29,-3]/840.0_rp

    contains


!------
! Reconstruction stencil routines
!------


    !> 4th order central interpolation for a single interface
    function Central4(Qcc) result(Qi)
        !! Note: still takes 8-element stencil
        real(rp), dimension(recLen), intent(in) :: Qcc
        
        real(rp) :: Qi

        Qi = dot_product(interpWgt_4c,Qcc(1:recLen))

    end function Central4


    !> 6th order central interpolation for a single interface
    function Central6(Qcc) result(Qi)
        !! Note: still takes 8-element stencil
        real(rp), dimension(recLen), intent(in) :: Qcc
        
        real(rp) :: Qi

        Qi = dot_product(interpWgt_6c,Qcc(1:recLen))

    end function Central6


    !> 8th order central interpolation for a single interface
    function Central8(Qcc) result(Qi)
        real(rp), dimension(recLen), intent(in) :: Qcc
        
        real(rp) :: Qi

        Qi = dot_product(interpWgt_8c,Qcc(1:recLen))

    end function Central8


    !> Reconstruct cell-centered variable at all faces bordering active cells
    subroutine ReconFaces(Grid, isG, Qcc, Qfaces, Qcc_phO)
        !! If just Qcc is provided, we interpolate Qcc to faces in both theta and phi directions
        !! If Qcc_phO is provided, we assume we interpolate Qcc to JUST theta-direction faces and Qcc_phO to JUST phi-direction faces
        !!   This is helpful when doing velocities because we don't care about the theta direction at phi faces and visa versa
        !! Note, we are still returning just one (Nig, Njg, 2) array
        type(raijuGrid_T), intent(in) :: Grid
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: isG
                           !! Whether a cell is safe to use in reconstruction
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: Qcc
                           !! Cell-centered variable to interpolate to faces
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2), intent(out) :: Qfaces
                           !! Face-interpolated variable
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in), optional :: Qcc_phO
                           !! Can optionally use a variable other than Qcc for phi direction

        integer :: i,j

        !! Note: We only populate from is:ie+1, js:je+1 because that's the only place we'll use face values
        !!  But full size includes ghost faces for output purposes

        !! Note: Loop is determining face i-1/2 and j-1/2 for given i,j
        Qfaces = 0.0

        do i=Grid%shGrid%is,Grid%shGrid%ie+1
            do j=Grid%shGrid%js,Grid%shGrid%je+1
                ! Theta dir
                Qfaces(i,j,1) = reconOrder( isG(i-4:i+3, j), Qcc(i-4:i+3, j) )
                ! Phi dir
                if( present(Qcc_phO) ) then
                    Qfaces(i,j,2) = reconOrder( isG(i, j-4:j+3), Qcc_phO(i, j-4:j+3) )
                else
                    Qfaces(i,j,2) = reconOrder( isG(i, j-4:j+3),     Qcc(i, j-4:j+3) )
                endif
            enddo
        enddo


        contains

        function reconOrder(isG, Qcc) result(Qface)
            !! Takes an 8-element stencil and determines reconstruction order based on isG
            logical , dimension(recLen), intent(in) :: isG
            real(rp), dimension(recLen), intent(in) :: Qcc

            real(rp) :: Qface
            
            !! TODO: We should add 4-th order asymmetric options before dropping to 2nd order
            !!   Will need to re-work ordering a bit once we do

            if ( all(isG) ) then
                !! Yay we can do full 8-th order centered reconstruction
                Qface = Central8(Qcc)
            else if ( all(isG(2:7)) ) then
                !! See if we can do 6-th order centered reconstruction
                Qface = Central6(Qcc)
            else if ( all(isG(3:6)) ) then
                !! 4-th order centered, hopefully do not need to drop below this
                Qface = Central4(Qcc)
            else if ( all(isG(4:5)) ) then
                !! Sad. Implement asymmetric 4-th order later to try and avoid this
                Qface = 0.5_rp*(Qcc(4) + Qcc(5))
            else
                !! If we are still here, one of the adjacent cells is bad
                !! We will set to zero, and then a BC function later on will set the flux of this face to whatever it needs to be
                Qface = 0.0_rp
            endif
            
        end function reconOrder

    end subroutine ReconFaces


    subroutine raijuPDM(Qcc, Qi, pdmb, Qpdm)
        !! Calculate PDM-limited value of Q at interface
        !! Using notation where interface is i+1/2
        !! Qcc is 4 elements long, i-1, i, i+1, i+2
        !!   Note, this is independent from indexing of arrays when they're handed to this subroutine
        !! Left is < i+1/2, Right is > i+1/2
        real(rp), dimension(4), intent(in) :: Qcc
            !! Cell-centered values of quantity. First two are left of interface, last 2 are right of interface
        real(rp), intent(in) :: Qi
            !! High-order reconstruction of quantity at interface
        real(rp), intent(in) :: pdmb
            !! PDM limiter value
        real(rp), intent(out) :: Qpdm
            !! PDM-limited quantity to return

        real(rp) :: minQ, maxQ, Qn
        real(rp) :: dQL, dQi, dQR, sL, si, sR, qL, qR
        real(rp) :: dQnL, dQnR, QLrecon, QRrecon

        ! Determine which is the intermediate value between Qi, Qi+1/2, Qi+1
        minQ = min(Qcc(2), Qcc(3))
        maxQ = max(Qcc(2), Qcc(3))
        Qn   = max(minQ, min(Qi, maxQ))

        ! Calculate differences to either side of interface, and at interface
        dQL = pdmb*( Qcc(2) - Qcc(1) )
        dQi = pdmb*( Qcc(3) - Qcc(2) )
        dQR = pdmb*( Qcc(4) - Qcc(3) )

        ! Signs of differences
        sL = sign(1.0_rp, dQL)
        si = sign(1.0_rp, dQi)
        sR = sign(1.0_rp, dQR)

        ! Will be 0 if local extrema, or 2 if monotonic
        qL = abs(sL + si)
        qR = abs(si + sR)

        ! Slope between median Qn and L/R cell centers
        ! Note, may be 0 in the case where Qi or Qi+1 was the median value
        dQnL = Qn     - Qcc(2)
        dQnR = Qcc(3) - Qn

        ! Finally, calculate L/R reconstructed values
        QLrecon = Qn - si*max(0.0, abs(dQnL) - qL*abs(dQL))
        QRrecon = Qn + si*max(0.0, abs(dQnR) - qR*abs(dQR))

        ! Final PDM-value is the one that's non-negative, so we can add them both to get our value
        Qpdm = QLrecon + QRrecon

    end subroutine raijuPDM


    subroutine PDMFaces(sh, isG, Qcc, Qfaces, Qpdms, pdmb)
        !! Calculates PDM-limited value for all faces for providd quantity
        type(ShellGrid_T), intent(in) :: sh
        logical , dimension(sh%isg:sh%ieg, &
                            sh%jsg:sh%jeg)     , intent(in) :: isG
        real(rp), dimension(sh%isg:sh%ieg, &
                            sh%jsg:sh%jeg)     , intent(in) :: Qcc
        real(rp), dimension(sh%isg:sh%ieg+1, &
                            sh%jsg:sh%jeg+1, 2), intent(in) :: Qfaces
        real(rp), dimension(sh%isg:sh%ieg+1, &
                            sh%jsg:sh%jeg+1, 2), intent(out) :: Qpdms
        real(rp), intent(in) :: pdmb

        integer :: i,j

        Qpdms = 0.0

        ! Only calculate for faces touching non-ghost cells
        do j=sh%js,sh%je+1
            do i=sh%is,sh%ie+1

                ! Theta dir first
                if (any(isG(i-2:i+1, j)) .eq. .false.) then
                    ! If full stencil isn't valid, we didn't use a high-order method to calculate the interface value anyways
                    !  so we can just use the value that's there
                    Qpdms(i,j,1) = Qfaces(i,j,1)
                else
                    call raijuPDM(Qcc(i-2:i+1, j), Qfaces(i, j, 1), pdmb, Qpdms(i,j,1))
                endif

                ! Phi dir
                if (any(isG(i, j-2:j+1)) .eq. .false.) then
                    ! If full stencil isn't valid, we didn't use a high-order method to calculate the interface value anyways
                    !  so we can just use the value that's there
                    Qpdms(i,j,2) = Qfaces(i,j,2)
                else
                    call raijuPDM(Qcc(i, j-2:j+1), Qfaces(i, j, 2), pdmb, Qpdms(i,j,2))
                endif

            enddo
        enddo

    end subroutine PDMFaces


    subroutine calcBoundaryFluxes(sh, active, Qflux)
        !! Sets fluxes at invalid/buffer boundaries
        type(ShellGrid_T), intent(in) :: sh
        integer , dimension(sh%isg:sh%ieg, &
                            sh%jsg:sh%jeg)     , intent(in) :: active
        real(rp), dimension(sh%isg:sh%ieg+1, &
                            sh%jsg:sh%jeg+1, 2), intent(inout) :: Qflux

        integer :: i,j

        do j=sh%js,sh%je+1
            do i=sh%is,sh%ie+1
                ! Theta dir
                if (  (active(i-1,j) .eq. RAIJUINACTIVE) &
                .and. (active(i  ,j) .eq. RAIJUBUFFER  ) &
                .and. (Qflux (i+1,j,1) < 0.0 ) ) then
                ! If current cell is good, lower i is bad, only copy i+1 flux if we have outflow
                    Qflux(i,j,1) = Qflux(i+1,j,1)
                else if (  (active(i-1,j) .eq. RAIJUBUFFER   ) &
                     .and. (active(i  ,j) .eq. RAIJUINACTIVE ) &
                     .and. (Qflux (i-1,j,1) > 0.0 ) ) then
                ! If current cell is bad, lower i is good, only copy i-1 flux if we have outflow
                    Qflux(i,j,1) = Qflux(i-1,j,1)
                endif


                ! Psi dir
                if (  (active(i,j-1) .eq. RAIJUINACTIVE) &
                .and. (active(i,j  ) .eq. RAIJUBUFFER  ) &
                .and. (Qflux (i,j+1,2) < 0.0 ) ) then
                ! If current cell is good, lower i is bad, only copy j+1 flux if we have outflow
                    Qflux(i,j,2) = Qflux(i,j+1,2)
                else if (  (active(i,j-1) .eq. RAIJUBUFFER   ) &
                     .and. (active(i,j  ) .eq. RAIJUINACTIVE ) &
                     .and. (Qflux (i,j-1,2) > 0.0 ) ) then
               ! If current cell is bad, lower i is good, only copy i-1 flux if we have outflow
                   Qflux(i,j,2) = Qflux(i,j-1,2)
               endif

            enddo
        enddo

    end subroutine calcBoundaryFluxes


    subroutine calcFluxes(Model, Grid, State, k, Qcc, Qflux)
        !! Takes cell-centered quantity (Qcc) and uses state information to calculate flux of Q through cell faces (Qflux)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: k
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg),      intent(in) :: Qcc
            !! Cell-centered quantity (e.x.: eta)
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2), intent(out) :: Qflux
            !! Flux of Q through faces

        ! Make some needed arrays here, figure out how to optimize later
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: isG
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: QA
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2) :: QAface, Qface, Qpdm

        Qflux = 0.0

        where (State%active .ne. RAIJUINACTIVE)
            isG = .true.
        elsewhere
            isG = .false.
        end where

        QA = Qcc * Grid%areaCC  ! Get total quantity within cell
        call ReconFaces(Grid, isG, QA, QAface)  ! Interpolate to face positions
        Qface = QAface / Grid%areaFace  ! Return to area density at face position

        ! PDM to get our final interface values
        !call PDMFaces(Grid%shGrid, isG, Qcc, Qface, Qpdm, Model%pdmb)

        ! Calculate face fluxes
        !Qflux = Qpdm*State%iVel(:,:,k,:)*Grid%lenFace &
        Qflux = Qface*State%iVel(:,:,k,:)*Grid%lenFace &
                * Model%planet%ri_m  ! [Q * m^2 / s]

        ! Thus far we have ignored fluxes of faces at invalid/buffer boundary
        !  (ReconFaces set them to zero)
        ! Now that valid faces are decided, we can apply our invalid/buffer boundary conditions
        !call calcBoundaryFluxes(Grid%shGrid, State%active, Qflux)
        
    end subroutine calcFluxes

end module raijuRecon