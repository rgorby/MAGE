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

        !! Note: We only populate from is:ie+1, js:je+1 because that's the only place we'll ues face values
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

end module raijuRecon