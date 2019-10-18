! Converting chimp data to gamera data

module chmp2mhd_interface
    use gamtypes
    use volttypes
    use gamapp
    use ebtypes
    
    implicit none

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

end module chmp2mhd_interface
